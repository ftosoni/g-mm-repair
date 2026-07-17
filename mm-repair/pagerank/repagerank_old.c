/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 * RePageRank (old version using 3 vectors)
 * A new version using only 2 vectors is in repagerank.c
 * 
 * Compute pagerank of a graph given its adjacency matrix
 * repair compressed, dandling nodes and teleporting included
 * 
 * Overview
 * In the original formulation PageRank requires the left multiplication 
 * of the current rank vector, with entries divided by the outdegree
 * of each node, times the binary adjacency matrix.
 * 
 * MatRepair does suport left and right matrix-vector multiplication 
 * but since other compressed matrix formats only support right
 * multiplication we are assuming that the graph matrix has been 
 * transposed (and all self-loops already removed).  
 * 
 * The main iteration goes as follows:
 *   1 each rank_i entry is divided by the outdegree of node i
 *     which is now the # of nonzero elements in column i
 *   2 if column i has no nonzero then i is a dandling node 
 *     and rank_i is added to the sum of dandling nodes ranks
 *   3 the normalized rank vector is right multiplied by the 
 *     (transpose adjacency matrix)
 *   4 the new rank vector is obtained for the above product
 *     + contribution of teleporting and dandling nodes.
 * In mathematical terms, let
 *   N = # nodes
 *   d = damping factor (from the command line)
 *   X = current rank vector
 *   dnr = 0  # dandling nodes rank sum
 *   for i in range(N):
 *     if col_count[i]==0: dnr[i] += X[i]
 *     else Y[i] = X[i]/col_count[i]
 *   Z = M*Y
 *   for i in range(N):
 *     Z[i] = d*Z[i] + (d/N)*dnr + (1-d)/N
 *   X = Z  # prepare for next iteration    
 * 
 * We start the iterations with X = (1/N ... 1/N)^T
 * we stop after maxiter (from the command line) iterations or
 * when the sum of the abs differences of the ranks between two
 * consecutive iterations is smaller than eps (from the command line)   
 *  
 * Note: in this program we execute the above iterations
 * as describe above, using three vectors X, Y and Z. However,
 * it is possible to use only two vectors X and Y at the expense
 * of some loss of accuracy in the error computation. The trick is
 * that in Y there is enough information to retrieve the previous X:
 *   1. when computing Y if col_count[i]==0 set Y[i] = X[i]
 *   2. instead of computing the new rank Z store it in X
 *   3. compute the error and Y for the next iteration as follows:
 *      error = 0
 *      for i in range(N):
 *        if col_count[i]==0: error += abs(X[i]-Y[i]); Y[i] = X[i]
 *        else error += abs(X[i]-Y[i]*col_count[i]);    Y[i] = X[i]/col_count[i]
 * 
 * Copyright 2024- giovanni.manzini@unipi.it
 * >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#ifdef CSR_MATRIX
#include "../csrmatrix.h"
#elif defined(USE_INTVEC) || defined(USE_ANSIV)
#include "../rematrix.hpp"
#else
#include "../rematrix.h"
#endif
#ifdef MALLOC_COUNT
#include "../tools/malloc_count.h"
#endif
#include <time.h>
#ifdef DETAILED_TIMING
#include <sys/times.h>
#endif
#include <pthread.h>
#include <semaphore.h>
#include <math.h>
#include "../tools/xerrors.h"

// input/output data for each thread 
typedef struct {
  rematrix *m;   // matrix block
  vector *rv;    // row vector    (same size as a matrix row)
  vector *cv;    // column vector (      "               column)
  int op;        // operation to be performed
  sem_t *in;     // semaphore for input shared with the main thread
  sem_t *out;    // semaphore for output shared with the main thread
} tdata;

static rematrix **remat_create_multipart(int, int,const char *, int blocks);
static void remat_destroy_multipart(rematrix **b,int n);
static void *block_main(void *v);
static void remat_mult_mth(rematrix *m,vector *x,vector *yv, tdata *td, int n);
static void minHeapify(double v[], int arr[], int n, int i);
static void kLargest(double v[], int arr[], int n, int k);



static void usage_and_exit(char *name)
{
    fprintf(stderr,"Usage:\n\t  %s [options] matrix msize col_count_file\n",name);
    fprintf(stderr,"\t\t-v             verbose\n");
    fprintf(stderr,"\t\t-b num         number of row blocks, def. 1\n");
    fprintf(stderr,"\t\t-m maxiter     maximum number of iteration, def. 100\n");
    fprintf(stderr,"\t\t-e eps         stop if error<eps (default ignore error)\n");
    fprintf(stderr,"\t\t-d df          damping factor (default 0.9)\n");
    fprintf(stderr,"\t\t-k K           show top K nodes (default 3)n\n");
    exit(1);
}



int main (int argc, char **argv) { 
  extern char *optarg;
  extern int optind, opterr, optopt;
  int verbose=0;
  int size,c;
  time_t start_wc = time(NULL);
  #ifdef DETAILED_TIMING
  struct tms ignored;
  clock_t t1,t2;
  long m1=0;
  #endif
  // default values for command line parameters 
  int maxiter=100,nblocks=1,topk=3;
  double dampf = 0.9, eps = -1;
  
  /* ------------- read options from command line ----------- */
  opterr = 0;
  while ((c=getopt(argc, argv, "b:m:e:d:k:v")) != -1) {
    switch (c) 
      {
      case 'v':
        verbose++; break;
      case 'm':
        maxiter=atoi(optarg); break;
      case 'b':
        nblocks=atoi(optarg); break;
      case 'e':
        eps=atof(optarg); break;
      case 'd':
        dampf=atof(optarg); break;
      case 'k':
        topk=atoi(optarg); break;
      case '?':
        fprintf(stderr,"Unknown option: %c\n", optopt);
        exit(1);
      }
  }
  if(verbose>0) {
    fputs("==== Command line:\n",stderr);
    for(int i=0;i<argc;i++)
     fprintf(stderr," %s",argv[i]);
    fputs("\n",stderr);  
  }
  // check command line
  if(maxiter<1 || nblocks < 1 || topk<1) {
    fprintf(stderr,"Error! Options -b -m and -k must be at least one\n");
    usage_and_exit(argv[0]);
  }
  if(dampf<0 || dampf>1) {
    fprintf(stderr,"Error! Options -d must be in the range [0,1]\n");
    usage_and_exit(argv[0]);
  }
  
  
  // virtually get rid of options from the command line 
  optind -=1;
  if (argc-optind != 4) usage_and_exit(argv[0]); 
  argv += optind; argc -= optind;
  
  // ----------- read and check matrix size 
  size  = atoi(argv[2]);
  if(size<1) die("Invalid matrix size");
  // ----------- init outdegree vector from file
  u_int32_t *outd = (u_int32_t *) malloc(size*sizeof(*outd));
  {
    FILE *ccol_file  = fopen(argv[3],"rb");
    if(ccol_file==NULL) die("Cannot open col_count_file");
    if(outd==NULL) die("Cannot allocate out_degree vector");
    size_t e = fread(outd,sizeof(*outd),size,ccol_file);
    if(e!=size) die("cannot read out_degree vector from col_count_file");
    if(fclose(ccol_file)!=0) die("Error closing col_count_file");
  }
  if(verbose>0) {
    fprintf(stderr,"Number of nodes: %d\n",size);
    long dn=0,arcs=0;
    for(int i=0;i<size;i++) {
      if(outd[i]==0) dn++;
      else arcs += outd[i];
    }
    fprintf(stderr,"Number of dandling nodes: %ld\n",dn);
    fprintf(stderr,"Number of arcs: %ld\n",arcs);
  }


  // ------------ read matrix or row blocks
  rematrix *m = NULL;
  rematrix **rblocks = NULL; 
  if(nblocks==1)
    m = remat_create(size,size,argv[1],true); 
  else 
    rblocks = remat_create_multipart(size,size,argv[1],nblocks);
    
  // ------------ init rank and aux vectors
  vector *x = vector_create_value(size,1.0/size);
  vector *y = vector_create_value(size,0);
  vector *z = vector_create_value(size,0);
  
  // data structures for multithread computation (nblocks>1)
  vector *zv = NULL;  // array of subvectors of z initialized below
  tdata td[nblocks];
  pthread_t t[nblocks];
  sem_t tsem_in[nblocks];
  sem_t tsem_out[nblocks];
  
  // initialize thread data
  if(nblocks>1) {
    // zv entries coincide with those of z  
    zv = vector_split(z,nblocks);
    for(int i=0;i<nblocks;i++) {
      td[i].m = rblocks[i];
      td[i].in = &tsem_in[i];
      td[i].out = &tsem_out[i];
      xsem_init(&tsem_in[i],0,0,__LINE__,__FILE__);
      xsem_init(&tsem_out[i],0,0,__LINE__,__FILE__);
      xpthread_create(&t[i],NULL,&block_main,&td[i],__LINE__,__FILE__);
    }
  }
    
  // compute products  
  int iter=0;
  double delta=11+eps;
  while(iter<maxiter && delta>=eps) {
    // normalize rank vector and compute dandling nodes rank
    double dnr = 0;
    for(int i=0;i<size;i++) {
      if(outd[i]==0) dnr += x->v[i];
      else y->v[i] = x->v[i]/outd[i];
    }
    #ifdef DETAILED_TIMING
    t1 = times(&ignored);
    #endif 
    if(nblocks==1) remat_mult(m,y,z);       // z = M*y
    else remat_mult_mth(m,y,zv,td,nblocks); // z = M*y with each thread computing a portion of z
    #ifdef DETAILED_TIMING
    t2 = times(&ignored);
    m1 += (t2-t1);       // measure time for matrix multiplication only
    #endif
    // compute contribution of teleporting and dandling nodes
    double teleport = (dnr*dampf+1-dampf)/size;
    // compute new rank vector and L1 norm of the difference
    delta = 0;
    for(int i=0;i<size;i++) {
      double nextri = dampf*z->v[i] + teleport;
      delta += fabs(nextri-x->v[i]);
      x->v[i] = nextri;
    }
    iter++; // iteration complete
    if(verbose>1) fprintf(stderr,"Iteration %d, delta=%g \n",iter,delta);
  }
  if(verbose>0) {
    if (delta>eps) fprintf(stderr,"Stopped after %d iterations, delta=%g\n",iter,delta);
    else           fprintf(stderr,"Converged after %d iterations, delta=%g\n",iter,delta);
    double sum=0;
    for(int i=0;i<size;i++) sum += x->v[i];
    fprintf(stderr,"Sum of ranks: %f (should be 1)\n",sum);
  }
  // deallocate y, z and outd: we may need space for the topk array
  vector_destroy(z);
  vector_destroy(y);
  free(outd);

  // retrieve topk nodes
  if(topk>size) topk = size;
  int *top = (int *) malloc(topk*sizeof(*top));
  int *aux = (int *) malloc(topk*sizeof(*top));
  if(top==NULL || aux==NULL) die("Cannot allocate topk/aux array");
  kLargest(x->v,aux,size,topk);
  // get sorted nodes in top
  for(int i=topk-1;i>=0;i--) {
    top[i] = aux[0];
    aux[0] = aux[i];
    minHeapify(x->v,aux,i,0);
  }  
  // report topk nodes sorted by decreasing rank
  if (verbose>0) {
    fprintf(stderr, "Top %d ranks:\n",topk);
    for(int i=0;i<topk;i++) fprintf(stderr,"  %d %lf\n",top[i],x->v[top[i]]);
  }
  // report topk nodes id's only on stdout
  fprintf(stdout,"Top:");
  for(int i=0;i<topk;i++) fprintf(stdout," %d",top[i]);
  fprintf(stdout,"\n");
  // destroy everything
  free(top); free(aux);
  vector_destroy(x);
  if(nblocks==1) 
    remat_destroy(m,true);
  else {
    free(zv);
    remat_destroy_multipart(rblocks,nblocks);
    for(int i=0;i<nblocks;i++) {
      td[i].op = -1; // stop thread
      xsem_post(td[i].in,__LINE__,__FILE__);
      pthread_join(t[i],NULL);
      xsem_destroy(td[i].in,__LINE__,__FILE__);
      xsem_destroy(td[i].out,__LINE__,__FILE__);
    }
  }
  #ifdef MALLOC_COUNT
    fprintf(stderr,"Peak memory allocation: %zu bytes, %.4lf bytes/node\n",
           malloc_count_peak(), (double)malloc_count_peak()/(size));
    fprintf(stderr,"Current memory allocation: %zu bytes\n", malloc_count_current());
  #endif
  #ifdef DETAILED_TIMING
  fprintf(stderr,"Total mult time (secs): %lf  Average: %lf\n", ((double)m1)/sysconf(_SC_CLK_TCK), ((double)m1/iter)/sysconf(_SC_CLK_TCK));
  #endif
  fprintf(stderr,"Elapsed time: %.0lf secs\n",(double) (time(NULL)-start_wc));  
  return 0;
}


// ------- code unchanged from remm.c ---------

// function executed by each thread 
// wait on a semaphore for a new operation to execute
// on its given matrix
static void *block_main(void *v)
{
  tdata *td = (tdata *) v;
  vector *auxrow = vector_create();
  vector_set_zero(auxrow,td->m->cols);  
  
  while(true) {
    // wait for input 
    xsem_wait(td->in,__LINE__,__FILE__);
    if(td->op<0) break;
    else if(td->op==0) { //left mult
      assert(td->cv!=NULL); // the input is a column vector
      if(td->rv==NULL) td->rv = auxrow;
      remat_left_mult(td->cv,td->m,td->rv);
    }
    else if(td->op==1) { //right mult
      assert(td->rv!=NULL); // the input is a row vector
      assert(td->cv!=NULL); // output must be given here 
      remat_mult(td->m,td->rv,td->cv);
    }
    else die("Unknown operation");
    // output ready 
    xsem_post(td->out,__LINE__,__FILE__);
  }
  vector_destroy(auxrow);
  return NULL;
}


// read matrix consisting of n blocks  
static rematrix **remat_create_multipart(int rows,int cols,const char *base, int n)
{
  assert(n>1); // there must be at least 2 blocks 
  
  rematrix **b = (rematrix **) malloc(n*sizeof *b);
  if(b==NULL) die("Not enough memory");
  int maxblock = (rows+n-1)/n;
  assert(maxblock>=1);
  
  // read everything except values
  #ifdef U_MATRIX
  FILE *fum = fopen(base,"r");
  if(fum==NULL) die("Unable to open matrix file")
  #else
  char fname[PATH_MAX];
  #endif  
  int remaining = rows;
  for(int i=0;i<n;i++) {
    assert(remaining>0);
    int r = (remaining>maxblock? maxblock : remaining);
    assert(r>0);
    #ifdef U_MATRIX
    b[i] = umat_create(r,col,fum);
    #else
    snprintf(fname,PATH_MAX,"%s.%d.%d",base,n,i);
    b[i] = remat_create(r,cols,fname,false);// false=> do not read .val file
    #endif
    remaining -= r;
  }
  assert(remaining==0);
  #ifdef U_MATRIX
  fclose(fum);
  #endif
  
  // read values ad assign them to all matrices in  b[] 
  snprintf(fname,PATH_MAX,"%s%s",base,VFILE_EXT);
  FILE *f = fopen(fname,"rb");
  if(f==NULL) die("Cannot open matrix values (" VFILE_EXT ") file");
  b[0]->Mval = read_vals(f,&b[0]->Mnum);
  // copy Mval/Mnum to the other blocks
  if(fclose(f)!=0) die("Error closing values (" VFILE_EXT ") file");
  for(int i=1;i<n;i++) {
    b[i]->Mval = b[0]->Mval; b[i]->Mnum = b[0]->Mnum;
  }  
  return b;
}


static void remat_destroy_multipart(rematrix **b,int n)
{

  free(b[0]->Mval); // free the unique copy of Mval we have   
  for(int i=0;i<n;i++)
    remat_destroy(b[i],false);
  free(b);
}

// note: x is a pointer to the left operand, yz is an array of n vectors 
// which are subvectors of the result
static void remat_mult_mth(rematrix *m,vector *x,vector *yv, tdata *td, int n)
{
  // start the block multiplications
  for(int i=0;i<n;i++) {
    td[i].rv = x;      // input  
    td[i].cv = &yv[i]; // output here 
    td[i].op = 1;      // right mult
    xsem_post(td[i].in, __LINE__,__FILE__);
  }
  // wait for all blocks
  for(int i=0;i<n;i++)
    xsem_wait(td[i].out, __LINE__,__FILE__); 
}



// heap based algorithm for finding the k largest ranks
// in heap order

// A utility function to swap two elements
static void swap(int* a, int* b) {
    int t = *a;
    *a = *b;
    *b = t;
}

// Heapify a subtree rooted with node i which is an index in arr[]. 
// n is size of heap. the key associated to entry arr[i] is v[arr[i]]  
static void minHeapify(double v[], int arr[], int n, int i) {
    int smallest = i;  // Initialize smallest as root
    int left = 2*i + 1;
    int right = 2*i + 2;

    // If left child is smaller than root
    if (left < n && v[arr[left]] < v[arr[smallest]])
        smallest = left;

    // If right child is smaller than smallest so far
    if (right < n && v[arr[right]] < v[arr[smallest]])
        smallest = right;

    // If smallest is not root
    if (smallest != i) {
        swap(&arr[i], &arr[smallest]);
        // Recursively heapify the affected sub-tree
        minHeapify(v, arr, n, smallest);
    }
}

// Function to find the k'th largest elements in an array
// v[0..n-1], arr[0..k-1] is the output array already allocated
static void kLargest(double v[], int arr[], int n, int k) {
  assert(k<=n);
  assert(k>0);
  // init arr[] with the first k elements
  for(int i=0;i<k;i++) arr[i] = i;
  // Build a min heap of first (k) elements in arr[]
  for (int i = k / 2 - 1; i >= 0; i--)
    minHeapify(v, arr, k, i);
  // Iterate through the rest of the array elements
  for (int i = k; i < n; i++) {
    // If current element is larger than root of the heap
    if (v[i] > v[arr[0]]) {
      // Replace root with current element
      arr[0] = i;
      // Heapify the root
      minHeapify(v, arr, k, 0);
    }
  }
}

