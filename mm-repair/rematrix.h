/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 * ReMatrix
 * 
 * operations on repair compressed matrices
 * 
 * The code here uses an int array for the R (rule array)
 * while the C sequence is compressed with ANS (if USE_ANS) 
 * or again with an int array.
 * Used by the executables ansremm and remm
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

#ifdef USE_ANS
  #include <iostream>
  #include <algorithm>
  #include <cstdint>
  #define ANSf 1
  #define CFILE_EXT ".vc.C.ansf.1"
  #define BUF_LOG2 20                  // log size decompression buffer  
  #define BUF_MASK ((1<<BUF_LOG2)-1)   // mask to recognize begining of buffer
  #include "ans/decode.hpp"
#else 
  #define CFILE_EXT ".vc.C"
#endif
#define RFILE_EXT ".vc.R"
#define VFILE_EXT ".val"

// set to 1 to print a lot of debug information 
#define DEBUG 0



// the use of float's has not been fully tested 
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


// report error message and terminates
//static void die(const char *s);
static void quit(const char *msg, int line, const char *file);
#define die(s) quit((s),__LINE__,__FILE__)


// include definitions for dense uncompressed vectors 
#include "vector.h"



// matrix represented as a re-pair grammar
typedef struct {
  int rows,cols;  // # rows and columns of input matrix 
  int Alpha;      // alphabet size of input matrix coincides with smallest non terminal symbol 
  matval *NTval;  // values associated to non-terminals (possibly shared with NTrules)
  size_t NTnum;   // number of non terminals 
  int *NTrules;   // right hand side of non-terminal 
  int *Cseq;      // array C of repair grammar
  size_t Clen;    // len of C array
#ifdef USE_ANS
  uint8_t *Ccseq; // ans-compressed array C of repair grammar
  size_t Cclen;   // length ans-compressed array
#endif  
  size_t Mnum;    // number of distinct non zero matrix values
  matval *Mval;   // set of distinct nonzero matrix values
  FILE *Cf;       // C file 
  FILE *Rf;       // R file  
} rematrix;  



// main prototypes
rematrix *remat_create(int r, int c, char *basename, bool read_vals);
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
// if read_vals==false then the actual matrix entries are not 
// read from the .val file and Mval and Mnum are not initialized
rematrix *remat_create(int r, int c, char *basename, bool read_values)
{
  char fname[PATH_MAX];
  FILE *f; struct stat s;
  rematrix *m=(rematrix *) malloc(sizeof(rematrix));
  if(m==NULL) die("Cannot alloc matrix");
  
  m->rows=r; m->cols=c;

  // ------------ open and read rule file
  if(strlen(basename)+20>PATH_MAX) die("Illegal base name");
  strcpy(fname,basename);
  strcat(fname,RFILE_EXT);
  if (stat (fname,&s) != 0)die("Cannot stat rule (" RFILE_EXT ") file");
  off_t len = s.st_size;
  m->Rf = fopen (fname,"rb");
  if (m->Rf == NULL) die("Cannot open rules (" RFILE_EXT ") file"); 

  // read alphabet size for the original input string 
  if (fread(&m->Alpha,sizeof(int),1,m->Rf) != 1)  
   die("Cannot read rules (" RFILE_EXT ") file (1)"); 
  m->NTnum = (len-sizeof(int))/(2*sizeof(int)); // number of non terminal (rules) 
  m->NTrules = (int *) malloc(m->NTnum*2*sizeof(int));
  if(m->NTnum>0 && m->NTrules==NULL) die("Cannot allocate R array");
  if(fread(m->NTrules,2*sizeof(int),m->NTnum,m->Rf)!=m->NTnum)
    die("Cannot read rules (" RFILE_EXT ") file (2)");
  if(fclose(m->Rf)) die("Error closing rules (" RFILE_EXT ") file");
    m->Rf=NULL;
  // values are used only for right multiplication, no need to allocate now   
  m->NTval = NULL;
  
  // --- open and read C file
  strcpy(fname,basename);
  strcat(fname,CFILE_EXT);
  if (stat (fname,&s) != 0) die("Cannot stat C (" CFILE_EXT ") file");
  m->Cf = fopen (fname,"r");
  if (m->Cf == NULL) die("Cannot open C (" CFILE_EXT ") file");
  #ifdef USE_ANS
  m->Cclen = s.st_size;
  m->Ccseq = new uint8_t[m->Cclen];
  if(fread(m->Ccseq,sizeof(uint8_t),m->Cclen,m->Cf)!=m->Cclen)
   die("Cannot read " CFILE_EXT " file");    
  // extract decompressed length and remove it from compressed data 
  m->Clen = *( (size_t *) m->Ccseq); // size of the decompressed C sequence 
  m->Ccseq += sizeof(size_t);
  m->Cclen -= sizeof(size_t);  
  // allocate decompressed buffer
  m->Cseq = (int *) malloc((1<<BUF_LOG2)*sizeof(int));
  if(m->Cseq==NULL) die("Cannot allocate buffer for ANS decompression");
  #else
  m->Clen = (s.st_size)/sizeof(int);
  m->Cseq = (int *) malloc(m->Clen*sizeof(int));
  if(fread(m->Cseq,sizeof(int),m->Clen,m->Cf)!=m->Clen)
   die("Cannot read " CFILE_EXT " file");   
  #endif 
  // close the file
  if(fclose(m->Cf)!=0) die("Error closing " CFILE_EXT " file"); 
  m->Cf = NULL;
  
  
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




// compute right multiplication: y = m x 
// right multiply the (rows x cols) matrix m by the
// vector x of size (cols x 1), obtaining y of size (rows x 1)
void remat_mult(rematrix *m, vector *x, vector *y)
{
  if(m->NTrules==NULL) die("Rules array missing (remat_mult)");
  if(m->cols!=x->size) die("Dimension mismatch (remat_mult x)");   
  if(m->rows!=y->size) die("Dimension mismatch (remat_mult y)");   

  // compute NT values according to x
  fill_NTval(m,x,false);
  // --- compute output 
  #ifdef USE_ANS
  // create and initialize decoder
  auto ans_dec = ANSf_decoder<ANSf>();
  ans_dec.init(m->Ccseq, m->Cclen, m->Clen);
  #endif
  int ycur = 0;
  xmatval sum=0;
  for(size_t j=0; j < m->Clen; j++) {
    #ifdef USE_ANS
    if((j & BUF_MASK) ==0) {
      size_t d = ans_dec.decode((uint32_t *)m->Cseq,std::min((size_t)(1<<BUF_LOG2), m->Clen -j));
      if(d==0) die("Illegal decode call");
    }
    int i = m->Cseq[j&BUF_MASK];
    #else
    int i = m->Cseq[j];
    #endif
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
  if(m->NTrules==NULL) die("Rules array missing (remat_left_mult)");
  if(m->rows!=y->size) die("Dimension mismatch (remat_left_mult y)");   
  if(m->cols!=x->size) die("Dimension mismatch (remat_left_mult x)");   
  // clean x
  for(size_t i=0;i<x->size;i++) x->v[i]=0;
  // clean NT values
  clear_NTval(m);

  // propagate y-values down the tree  

  // variables used by decode_entry 
  #ifdef USE_ANS
  // create and initialize decoder
  auto ans_dec = ANSf_decoder<ANSf>();
  ans_dec.init(m->Ccseq, m->Cclen, m->Clen);
  #endif
  xmatval a; size_t col;   
  // propagate y-values to symbols in C
  int ycur=0;
  for(size_t j=0; j<m->Clen;j++) {  
    #ifdef USE_ANS
    if((j & BUF_MASK) ==0) {
      size_t d = ans_dec.decode((uint32_t *)m->Cseq, std::min((size_t) (1<<BUF_LOG2), m->Clen -j));
      if(d==0) die("Illegal decode call");
    }
    int i = m->Cseq[j&BUF_MASK];
    #else
    int i = m->Cseq[j];
    #endif
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
  // these are usually accessed sequentially, so one could keep them on file
  // this feature is curenntly not used!
  if(m->NTrules) {free(m->NTrules); m->NTrules=NULL;}
  if(m->NTval)   {free(m->NTval); m->NTval=NULL;}
  if(m->Cseq)    {free(m->Cseq); m->Cseq=NULL;}
  #ifdef USE_ANS
  if(m->Ccseq!=NULL) {
    delete[] (m->Ccseq - sizeof(size_t)); // see initialization of m->Ccseq
    m->Ccseq = NULL;
  } 
  #endif

  // if the files were left open in case their data have to be reloaded close them now
  // currently we are not suing this and close the files immediately  
  if(m->Rf!=NULL) { 
    if(fclose(m->Rf)) die("Error closing rules (" RFILE_EXT ") file");
    m->Rf=NULL;
  }
  if(m->Cf!=NULL) {
    if(fclose(m->Cf)) die("Error closing " CFILE_EXT " file");
    m->Cf=NULL;
  }
  free(m);
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
// used in left multiplication 
static void propagate_NTval(rematrix *m, vector *x) 
{
  // variables used by decode_entry 
  xmatval a; size_t col;   
  // scan rules in reverse order
  int *pair = m->NTrules + 2*m->NTnum;  // pointer beyond the last rule 
  for(ssize_t i=m->NTnum-1;i>=0;i--) {
    xmatval s = m->NTval[i]; // value associated to rule i
    pair -= 2;               // go back one rule in m->NTrules
    for(int j=0;j<2;j++) {   // propagate s to rule i right-hand-size   
      int p = pair[j];
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
  assert(pair==m->NTrules);
}


// compute the value associated to each non-terminal
// possibly overwriting the rule array if share==true
// used by right multiplication 
// Note: currently we are always using share=false since the rule
// array will be needed in any additional operation  
static void fill_NTval(rematrix *m, vector *x, bool share) 
{
  if(m->NTrules==NULL) die("Invalid call to fill_NTval");

  if(share) {
    // make sure the overwriting is possible 
    assert(sizeof(matval)<=2*sizeof(int));
    // share the same memory space 
    m->NTval = (matval *) m->NTrules;
  }
  else {
    if(m->NTval==NULL) {
      m->NTval = (matval *) malloc(m->NTnum*sizeof(matval));
      if(m->NTval==NULL) die("Cannot allocate NTval (fill)");
    }
  } 
  
  // compute values 
  int *pair = m->NTrules;
  for(int i=0; i<(ssize_t) m->NTnum;i++) {  // i is number of NT computed so far
    xmatval sum = 0;
    for(int j=0;j<2;j++) {
      int p = pair[j];
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
    pair += 2;   // advance to next rule
    m->NTval[i]=sum;
    #ifndef INT_VALS 
    if(DEBUG) printf("  NT[%d]: %f\n",i,sum); //!!!!!!!!!
    #else 
    if(DEBUG) printf("  NT[%d]: %d\n",i,sum); //!!!!!!!!!
    #endif
  }
  // make NTrules[] invalid
  if(share) m->NTrules = NULL;
  return;  
}


// write error message + extra info and and exit
static void quit(const char *msg, int line, const char *file) {
  if(errno==0)  fprintf(stderr,"== %d == %s\n",getpid(), msg);
  else fprintf(stderr,"== %d == %s: %s\n",getpid(), msg,
               strerror(errno));
  fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
  exit(1);
}

