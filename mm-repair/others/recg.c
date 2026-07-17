/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 * ReCG
 * 
 * coniugate gradient on repair compressed matrices
 * Given rxc matrix A and 1xc vector y computes rx1 vector x
 * that minimizes ||Ax-y||_2 by computing n steps 
 * of CG descent algorithm from the CLA paper
 * 
 * input matrix A must be in repair compressed format
 * input vector y must be in binary format
 * >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#define _GNU_SOURCE
#ifdef CSR_MATRIX
#include "csrmatrix.h"
#else
#include "rematrix.h"
#endif
#ifdef MALLOC_COUNT
#include "mc/malloc_count.h"
#endif
#include <math.h>

static void usage_and_exit(char *name)
{
    fprintf(stderr,"Usage:\n\t  %s [-n mul] matrix rows cols yvector\n",name);
    fprintf(stderr,"\t\t-n mul         number of multiplications, def. 1\n\n");
    exit(1);
}


int main (int argc, char **argv) { 
  extern char *optarg;
  extern int optind, opterr, optopt;
  FILE *f;
  int rows,cols,c,iter=1;
  // ----------- check input
  fputs("==== Command line:\n",stderr);
  for(int i=0;i<argc;i++)
   fprintf(stderr," %s",argv[i]);
  fputs("\n",stderr);  
  
  /* ------------- read options from command line ----------- */
  opterr = 0;
  while ((c=getopt(argc, argv, "n:")) != -1) {
    switch (c) 
      {
      case 'n':
        iter=atoi(optarg); break;
      case '?':
        fprintf(stderr,"Unknown option: %c\n", optopt);
        exit(1);
      }
  }
  // check command line
  if(iter<1) {
    fprintf(stderr,"Error! Option -n must be at least one\n");
    usage_and_exit(argv[0]);
  }  
  // virtually get rid of options from the command line 
  optind -=1;
  if (argc-optind != 5) usage_and_exit(argv[0]); 
  argv += optind; argc -= optind;
  
  // ----------- read and check # rows and cols 
  rows  = atoi(argv[2]);
  if(rows<1) die("Invalid number of rows");
  cols  = atoi(argv[3]);
  if(cols<1) die("Invalid number of columns");
  // ------------ read input vector
  f = fopen(argv[4],"rb");
  if(f==NULL) die("Cannot open input vector file");
  vector *y = vector_create();
  y->v = read_vals(f,&y->size);
  if(y->size!=rows) die("Input vector size should be equal to # of columns");
  fclose(f);
  // ------------ read matrix
  rematrix *m = remat_create(rows,cols,argv[1]); 
  
  // init CG
  double lambda = 0; // in the CLA paper 0.001 is used 
  vector *r = vector_create();
  vector *p = vector_create();
  vector *w = vector_create();
  vector *q = vector_create();
  vector *z = vector_create();
  
  vector_set_zero(w,cols); // w=0
  vector_set_zero(p,cols); // p=0
    
  remat_left_mult(y,m,r);         // note: y is no longer used
  vector_scalar_update(r,-1);     // r = - (y^t m)
  xmatval n2 = vector_norm2sq(r); // |r|_2^2
  vector_update(p,-1,r);          // p = -r  
  
  for(int i=0;i<iter;i++) {
    // compute conjugate gradient 
    remat_mult(m,p,y);      // y = m p 
    remat_left_mult(y,m,q); // q = m^t m p
    vector_update(q,lambda,p); // q = q + lambda p
    // compute step size
    xmatval alpha = (n2) / (vector_scalar_prod(p,q));
    //printf("(%lf %lf) %lf %lf %lf\n",n2,vector_scalar_prod(p,q),w->v[0],w->v[4],w->v[20]);
    // update model and residual
    vector_update(w,alpha,p);
    vector_update(r,alpha,q);
    xmatval old_n2 = n2;
    n2 = vector_norm2sq(r);
    vector_scalar_update(p,n2/old_n2);
    vector_update(p,-1,r);
    printf("%lg -> %lg\n",old_n2,n2);
    printf("(norm2=%lg alpha=%lg)   w[0]=%lg\n",n2,alpha,w->v[0]);
  }
  
  // display residual and write solution to file
  printf("Residual after %d iterations: %lf\n",iter,(double) n2);
  char *fnam;
  c = asprintf(&fnam,"%s.%d.sol",argv[4],iter);
  if(c<0) die("Cannot create .sol file name");
  f= fopen(fnam,"wb");  
  if (f == NULL) die("Cannot open output vector file");
  size_t e = fwrite(w->v,sizeof(matval),w->size,f);
  if(e!=w->size) die("Cannot write to output file");
  if(fclose(f)!=0) die("Cannot close output file");
  free(fnam);
  
  // destroy 
  vector_destroy(z);
  vector_destroy(q);
  vector_destroy(w);
  vector_destroy(p);
  vector_destroy(r);
  vector_destroy(y);
  remat_destroy(m);
  #ifdef MALLOC_COUNT
    printf("Peak memory allocation: %zu bytes, %.4lf bytes/entries\n",
           malloc_count_peak(), (double)malloc_count_peak()/(rows*cols));
  #endif
  return 0;
}

