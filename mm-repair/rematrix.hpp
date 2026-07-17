/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 * ReMatrix
 * 
 * operations on repair compressed matrices
 * 
 * The code here uses an sdsl int_vector for the R (rule) array
 * while the C sequence is compressed with ANS (if USE_ANSIV) 
 * or again with an int_vector.
 * Used by the executables ansivremm and ivremm
 *  
 * Copyright (C) 2021-2099   giovanni.manzini@unipi.it
 * >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>
#include <errno.h>

#include <sdsl/int_vector.hpp>
#include <iostream>
#include <string>

#ifdef USE_ANSIV
  #include <algorithm>
  #include <cstdint>
  #define ANSf 1
  #define CFILE_EXT ".vc.C.ansf.1"
  #include "ans/decode.hpp"
  #define RFILE_EXT ".vc.R.iv"
  #define BUF_LOG2 20                  // log of (size decompression buffer)  
#else
  #define CFILE_EXT ".vc.C.iv"
  #define RFILE_EXT ".vc.R.iv"
  #define BUF_LOG2 18                  // log of (size decompression buffer)  
#endif

#define BUF_MASK ((1<<BUF_LOG2)-1)     // mask to recognize beginning of buffer

#define VFILE_EXT ".val"

// set to 1 to print a lot of debug information 
#define DEBUG 0


// the use of float's has not been fully tested 
//  matval: type representing a matrix/vector entry
// xmatval: type representing a matrix entry with larger precision (for intermediate results)
#ifndef INT_VALS 
  #ifndef FLOAT_VALS
  typedef double matval;    // type representing a matrix/vector entry
  #else
  typedef float  matval;    // type representing a matrix/vector entry
  #endif     
  typedef double xmatval;   // type representing a matrix entry with larger precision   
#else
  typedef int matval;       // type representing a matrix entry   
  typedef int xmatval;      // type representing a matrix entry with larger precision   
#endif


// report error message and terminate
static void quit(const char *msg, int line, const char *file);
#define die(s) quit((s),__LINE__,__FILE__)


// include definitions for dense uncompressed vectors 
#include "vector.h"


// matrix represented as a re-pair grammar
typedef struct {
  int rows,cols;  // # rows and columns of input matrix 
  int Alpha;      // alphabet size of input matrix: coincides with smallest non terminal symbol 
  matval *NTval;  // values associated to non-terminals during computation 
  size_t NTnum;   // number of non terminals 
  sdsl::int_vector<> NTrules; // rule array from repair: length = 2*NTnum+1
  size_t Clen;    // len of C array
  int *Cseq;      // uncompressed array C from repair; here is only a buffer of size 1<<BUF_LOG2
#ifdef USE_ANSIV
  uint8_t *Ccseq; // ans-compressed array C from repair
  size_t Cclen;   // length of ans-compressed array
#else
  sdsl::int_vector<> Ccseq; // int-vector C from repair
#endif  
  size_t Mnum;    // number of distinct non zero matrix values
  matval *Mval;   // set of distinct nonzero matrix values
} rematrix;  



// main prototypes
rematrix *remat_create(int r, int c, char *basename, bool read_values);
void remat_destroy(rematrix *v, bool free_vals);
void remat_mult(rematrix *m, vector *x, vector *y);
void remat_left_mult(vector *y, rematrix *m, vector *x);
matval *read_vals(FILE *f, size_t* size);
xmatval decode_mult_entry(int p, rematrix *m, vector *x);
xmatval decode_entry(int p, rematrix *m, size_t *c);

// local functions 
static void fill_NTval(rematrix *m, vector *x, bool share);
static void clear_NTval(rematrix *m);
static void propagate_NTval(rematrix *m, vector *x);

// load compressed matrix information from files
rematrix *remat_create(int r, int c, char *basename, bool read_values)
{
  char fname[PATH_MAX];
  FILE *f; struct stat s;
  rematrix *m= new rematrix;  
  m->rows=r; m->cols=c;

  // ------------ open and read rule file
  if(strlen(basename)+20>PATH_MAX) die("Illegal base name");
  strcpy(fname,basename);
  strcat(fname,RFILE_EXT);
  if (stat (fname,&s) != 0)die("Cannot stat rule (" RFILE_EXT ") file");
  load_from_file(m->NTrules, std::string(fname));
  assert(m->NTrules.size()%2==1); // the first entry is the alpha size
  m->Alpha = m->NTrules[0];
  m->NTnum = m->NTrules.size()/2;
  
  // values are used only for right multiplication, no need to allocate now   
  m->NTval = NULL;
  
  // --- open and read C file
  strcpy(fname,basename);
  strcat(fname,CFILE_EXT);
  if (stat (fname,&s) != 0) die("Cannot stat C (" CFILE_EXT ") file");
#ifdef USE_ANSIV  
  f = fopen (fname,"r");
  if (f == NULL) die("Cannot open C (" CFILE_EXT ") file");
  m->Cclen = s.st_size;
  m->Ccseq = new uint8_t[m->Cclen];
  if(fread(m->Ccseq,sizeof(uint8_t),m->Cclen,f)!=m->Cclen)
   die("Cannot read " CFILE_EXT " file");    
  if(fclose(f)) 
    die("Error closing C (" CFILE_EXT ") file");
  // extract decompressed length and remove it from compressed data 
  m->Clen = *( (size_t *) m->Ccseq); // size of the decompressed C sequence 
  m->Ccseq += sizeof(size_t);
  m->Cclen -= sizeof(size_t);  
#else
  // much simpler if Ccseq is just an int_vector
  load_from_file(m->Ccseq,std::string(fname));
  m->Clen = m->Ccseq.size();
#endif
  // allocate decompressed buffer
  m->Cseq = (int *) malloc((1<<BUF_LOG2)*sizeof(int));
  if(m->Cseq==NULL) die("Cannot allocate buffer for ANS decompression");

  // ------------ read matrix values 
  if(read_values) {
    strcpy(fname,basename);
    strcat(fname,VFILE_EXT);
    f = fopen(fname,"rb");
    if(f==NULL) die("Cannot open matrix values (" VFILE_EXT ") file");
    m->Mval = read_vals(f,&m->Mnum);
    if(fclose(f)!=0) die("Error closing values (" VFILE_EXT ") file");
  }
  else {
    m->Mval=NULL; m->Mnum=0;
  }

  return m;
}

// right multiplication 
void remat_mult(rematrix *m, vector *x, vector *y)
{
  if(m->cols!=x->size) die("Dimension mismatch (remat_mult)");   
  if(m->rows!=y->size) die("Dimension mismatch (remat_mult)");   
  // compute NT values according to x
  fill_NTval(m,x,false);

  // --- compute output 
  #ifdef USE_ANSIV
  // create and initialize decoder
  auto ans_dec = ANSf_decoder<ANSf>();
  ans_dec.init(m->Ccseq, m->Cclen, m->Clen);
  #endif
  int ycur = 0;
  xmatval sum=0;
  for(size_t j=0; j < m->Clen; j++) {
    if((j & BUF_MASK) ==0) {
      // fill the buffer 
      size_t to_read = std::min((size_t)(1<<BUF_LOG2), m->Clen -j);
      #ifdef USE_ANSIV
      size_t d = ans_dec.decode((uint32_t *)m->Cseq,to_read);
      if(d==0) die("Illegal decode call");
      #else
      for(size_t d = 0; d < to_read; d++)
        m->Cseq[d] = m->Ccseq[d+j];  // decode a bunch of packed ints
      #endif
    }
    int i = m->Cseq[j&BUF_MASK]; // read a single int from buffer    
    if(i>=m->Alpha) { // non terminal 
     if( (i = i-m->Alpha)>= (int)m->NTnum ) die("Illegal non terminal in C file");
     sum += m->NTval[i];
    }
    else if(i>0) {// terminal representing a matrix entry
     sum += decode_mult_entry(i-1,m,x);
    }
    else { // i==0 row completed
     y->v[ycur] = (matval) sum;
     sum = 0;
     if(++ycur==y->size) assert(j+1==m->Clen);
    }
  }
  assert(ycur==y->size);
  assert(sum==0);
}


// left multiply the (rows x cols) matrix m by the
// vector y^T of size (1 x rows), obtaining x^T of size (1 x cols)
void remat_left_mult(vector *y, rematrix *m, vector *x)
{
  // make sure the rules are available and dimensions agree
  if(m->rows!=y->size) die("Dimension mismatch (remat_left_mult)");   
  if(m->cols!=x->size) die("Dimension mismatch (remat_left_mult)");   
  // clean x
  for(size_t i=0;i<x->size;i++) x->v[i]=0;
  // clean NT values
  clear_NTval(m);

  // propagate y-values down the tree  
  // variables used by decode_entry 
  #ifdef USE_ANSIV
  // create and initialize decoder
  auto ans_dec = ANSf_decoder<ANSf>();
  ans_dec.init(m->Ccseq, m->Cclen, m->Clen);
  #endif
  xmatval a; size_t col;   
  // propagate y-values to symbols in C
  int ycur=0;
  for(size_t j=0; j<m->Clen;j++) {  
    if((j & BUF_MASK) ==0) {
      size_t to_read = std::min((size_t)(1<<BUF_LOG2), m->Clen -j);
      #ifdef USE_ANSIV
      size_t d = ans_dec.decode((uint32_t *)m->Cseq,to_read);
      if(d==0) die("Illegal decode call");
      #else
      for(size_t d = 0; d < to_read; d++)
        m->Cseq[d] = m->Ccseq[d+j]; // decode a bunch of packed ints
      #endif
    }
    int i = m->Cseq[j&BUF_MASK];    // read a single int from buffer
    if(i>=m->Alpha) { // non terminal 
      if( (i = i-m->Alpha)>= (int)m->NTnum ) die("Illegal non terminal in C file (left-mult)");
      assert(ycur < y->size);
      assert(i<m->NTnum);
      m->NTval[i] += y->v[ycur];
    }
    else if(i>0) {// terminal representing a matrix entry
      a = decode_entry(i-1,m,&col);
      assert(col<x->size);
      x->v[col] += a * y->v[ycur];
    }
    else { // i==0 row completed
      if(++ycur==y->size) assert(j+1==m->Clen);
    }
  }
  assert(ycur==y->size);
  // propagate y-values to x through the set of rules
  propagate_NTval(m,x);
}


void remat_destroy(rematrix *m, bool free_vals)
{  
  if(free_vals) free(m->Mval);
  if(m->NTval)   {free(m->NTval); m->NTval=NULL;}
  sdsl::util::clear(m->NTrules);
#ifdef USE_ANSIV
  if(m->Ccseq!=NULL) {
    delete[] (m->Ccseq - sizeof(size_t)); // see initialization of m->Ccseq
    m->Ccseq = NULL;
  } 
#else
  sdsl::util::clear(m->Ccseq);
#endif
  if(m->Cseq)    {free(m->Cseq); m->Cseq=NULL;}
  delete m;
}


// get value and column from terminal representing a matrix entry
xmatval decode_entry(int p, rematrix *m, size_t *c)
{
  *c = p % m->cols;
  size_t pval = p/m->cols;
  if(pval>=m->Mnum) die("Illegal value reference found in terminal symbol (decode_entry)");
  return m->Mval[pval];  
}


// decodes a terminal representing a matrix entry
// does not return the column index, instead the matrix
// value is multiplied by the corresponding X entry 
xmatval decode_mult_entry(int p, rematrix *m, vector *x)
{
  size_t pcol = p % m->cols;
  size_t pval = p/m->cols;
  if(pval>=m->Mnum) die("Illegal value reference found in terminal symbol (decode_mult_entry)");
  assert(pcol<x->size);
  return ((xmatval) x->v[pcol])*m->Mval[pval];
}  

// read a set of matval values from file f
// return number of items in *n and pointer to array with values
matval *read_vals(FILE *f, size_t *n)
{
  // get files size
  if(fseek(f,0,SEEK_END))
    die("Error in read_vals:fseek");
  long size = ftell(f);
  if(size<0)
    die("Error in read_vals:ftell");
  // get array size  
  *n = (int) (size/sizeof(matval));
  matval *a  = (matval *)malloc(*n * sizeof(matval));
  if(a==NULL)
    die("Error in read_vals:malloc");
  rewind(f);
  size_t e = fread(a,sizeof(matval),*n,f);
  if(e!= *n)
    die("Error in read_vals:fread");
  return a;
}

// init to zero all entries of NTval
// possibly allocate the array, but otherwise just clean it 
// loosing the original content 
static void clear_NTval(rematrix *m)
{
  if(m->NTval==NULL) {
    m->NTval = (matval *) malloc(m->NTnum*sizeof(matval));
    if(m->NTval==NULL) die("Cannot allocate NTval (clear)");
  }
  for(size_t i=0;i<m->NTnum;i++) m->NTval[i] = 0;
}

// propagate NTvalues up to the terminals and the x array
// the NTrules array is read in right to left order
// used in left multiplication 
static void propagate_NTval(rematrix *m, vector *x) 
{
  // variables used by decode_entry 
  xmatval a; size_t col;   
  // scan rules in reverse order
  int pos = 1 + 2*m->NTnum;  // pointer beyond the last rule 
  for(ssize_t i=m->NTnum-1;i>=0;i--) {
    xmatval s = m->NTval[i]; // value associated to rule i
    pos -= 2;                // go back one rule in m->NTrules
    for(int j=0;j<2;j++) {   // propagate s to rule i right-hand-size   
      int p = m->NTrules[pos+j];
      if(p>=m->Alpha) { // non terminal
        p -= m->Alpha;
        if(p>=i) die("Fatal Error: Forward rule (propagate_NTval)");
        m->NTval[p] += s;    // add s to value of p-th non-terminal
      }
      else { // terminal symbol
        if(p==0) die("Unique row separator found in rule (propagate_NTval)");
        a = decode_entry(p-1,m,&col);
        assert(col<x->size);
        x->v[col] += a * s;
      }
    }
  }
  assert(pos==1);
}


// compute the value associated to each non-terminal
// cannot overwrite the rule array so share must be false
// used by right multiplication 
static void fill_NTval(rematrix *m, vector *x, bool share) 
{
  assert(!share); // share not possible if R is compressed 

  if(m->NTval==NULL) {
      m->NTval = (matval *) malloc(m->NTnum*sizeof(matval));
      if(m->NTval==NULL) die("Cannot allocate NTval (fill)");
  }
  
  // compute values 
  int pos = 1; // start of first rule
  for(int i=0; i<(ssize_t) m->NTnum;i++) {  // i is number of NT computed so far
    xmatval sum = 0;
    for(int j=0;j<2;j++) {
      int p = m->NTrules[pos+j];
      if(p>=m->Alpha) { // non terminal
        p -= m->Alpha;
        if(p>=i) die("Fatal Error: Forward rule (fill_NTval)");
        sum += m->NTval[p];
        if(DEBUG) printf("%d-nt:%d  ",j,p);//!!!!!!!!!!!!
      }
      else { // terminal symbol
        if(p==0) {
          if(DEBUG) printf("%d-sep: %d  ",j,p);//!!!!!!!!
          die("Unique row separator found in rule");
        }
        sum += decode_mult_entry(p-1,m,x);
        #ifndef INT_VALS 
        if(DEBUG) printf("%d-t: col:%d val:%f ",j,(p-1)%m->cols,m->Mval[(p-1)/m->cols]);//!!!!!!!
        #else
        if(DEBUG) printf("%d-t: col:%d val:%d ",j,(p-1)%m->cols,m->Mval[(p-1)/m->cols]);//!!!!!!!
        #endif
      }
    }
    pos += 2;   // advance to next rule
    m->NTval[i]=sum;
    #ifndef INT_VALS 
    if(DEBUG) printf("  NT[%d]: %f\n",i,sum); //!!!!!!!!!
    #else 
    if(DEBUG) printf("  NT[%d]: %d\n",i,sum); //!!!!!!!!!
    #endif
  }
  return;  
}

// write error message + extra info and exit
static void quit(const char *msg, int line, const char *file) {
  if(errno==0)  fprintf(stderr,"== %d == %s\n",getpid(), msg);
  else fprintf(stderr,"== %d == %s: %s\n",getpid(), msg,
               strerror(errno));
  fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
  exit(1);
}

