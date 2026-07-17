/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 * CsrMatrix
 * 
 * operations on Compressed Sparse Row Representation 
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

#define VFILE_EXT ".val"

// set to 1 to print a lot of debug information 
#define DEBUG 0


// the use of float's has not been fully tested 
#ifndef INT_VALS 
  #ifndef FLOAT_VALS
  typedef double matval;  // type representing a matrix/vector entry
  #else
  typedef float  matval;  // type representing a matrix/vector entry
  #endif     
  typedef double xmatval;   // type representing a matrix entry with larger precision   
#else
  typedef int matval;       // type representing a matrix entry   
  typedef int xmatval;      // type representing a matrix entry with larger precision   
#endif

typedef double  matval;    // type representing a matrix/vector entry
typedef double xmatval;   // type representing a matrix entry with larger precision   

// report error message and terminates
static void quit(const char *msg, int line, char *file);
#define die(s) quit((s),__LINE__,__FILE__)


// include definitions for vectors dense uncompressed vectors 
#include "vector.h"



// matrix represented as a re-pair grammar
typedef struct {
  int rows,cols;  // # rows and columns of input matrix 
  int *CSRseq;      // array C of repair grammar
  size_t CSRlen;    // len of C array
  size_t Mnum;    // number of distinct non zero matrix values
  matval *Mval;   // set of distinct nonzero matrix values
  FILE *CSRf;       // C file 
} rematrix;  



// main prototypes
rematrix *remat_create(int r, int c, char *basename, bool read_vals);
void remat_destroy(rematrix *v, bool free_vals);
void remat_mult(rematrix *m, vector *x, vector *y);
matval *read_vals(FILE *f, size_t* size);
xmatval decode_mult_entry(int p, rematrix *m, vector *x);
xmatval decode_entry(int p, rematrix *m, size_t *c);


rematrix *remat_create(int r, int c, char *basename,bool read_values)
{
  char fname[PATH_MAX];
  FILE *f; struct stat s;
  rematrix *m=malloc(sizeof(rematrix));
  if(m==NULL) die("Cannot allocate matrix");
  
  m->rows=r; m->cols=c;

  // ------------ read csr values
  if(strlen(basename)+10>PATH_MAX) die("Illegal base name");
  strcpy(fname,basename);
  strcat(fname,".vc");
  if (stat (fname,&s) != 0) die("Cannot stat csr (.vc) file");
  m->CSRf = fopen (fname,"r");
  if (m->CSRf == NULL)      die("Cannot open csr (.vc) file");
  m->CSRlen = (s.st_size)/sizeof(int);
  m->CSRseq = (int *) malloc(m->CSRlen*sizeof(int));
  if(fread(m->CSRseq,sizeof(int),m->CSRlen,m->CSRf)!=m->CSRlen)
   die("Cannot read .vc file");
  
  // ------------ read matrix values 
  if(read_values) {
    strcpy(fname,basename);
    strcat(fname,".val");
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
  if(m->cols!=x->size) die("Dimension mismatch (remat_mult x)");   
  if(m->rows!=y->size) die("Dimension mismatch (remat_mult y)");   

  // --- compute output 
  int ycur = 0;
  xmatval sum=0;
  for(size_t j=0; j < m->CSRlen; j++) {
    int i = m->CSRseq[j];
    if(i>0) {// symbol representing a matrix entry
     sum += decode_mult_entry(i-1,m,x);
    }
    else { // i==0 row completed
     y->v[ycur] = (matval) sum;
     sum = 0;
     if(++ycur==y->size) assert(j+1==m->CSRlen);
    }
  }
  assert(ycur==y->size);
  assert(sum==0);
}


// left multiply the (rows x cols) matrix m by the
// vector y^T of size (1 x rows), obtaining x^T of size (1 x cols)
void remat_left_mult(vector *y, rematrix *m, vector *x)
{
  // make sure dimensions agree
  if(m->rows!=y->size) die("Dimension mismatch (remat_left_mult y)");   
  if(m->cols!=x->size) die("Dimension mismatch (remat_left_mult x)");   
  // clean x
  for(size_t i=0;i<x->size;i++) x->v[i]=0;

  // variables used by decode_entry 
  xmatval a; size_t col;   
  // propagate y-values to symbols in C
  int ycur=0; // ycur is the current row index 
  for(size_t j=0; j<m->CSRlen;j++) {  
    int i = m->CSRseq[j];
    if(i>0) {         // symbol representing a matrix entry
      a = decode_entry(i-1,m,&col);
      assert(col<x->size);
      x->v[col] += a * y->v[ycur];
    }
    else { // i==0 row completed
      if(++ycur==y->size) assert(j+1==m->CSRlen);
    }
  }
  assert(ycur==y->size);
}


void remat_destroy(rematrix *m, bool free_vals)
{  
  if(free_vals) free(m->Mval);
  // these are usually accessed sequentially, so one could keep them on file
  if(m->CSRseq)    {free(m->CSRseq); m->CSRseq=NULL;}

  // the file .vc were left open in case data have to be reloaded or read form file
  if(m->CSRf!=NULL) {
    if(fclose(m->CSRf)) die("Error closing .vc file");
    m->CSRf=NULL;
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

// write error message + extra info and and exit
static void quit(const char *msg, int line, char *file) {
  if(errno==0)  fprintf(stderr,"== %d == %s\n",getpid(), msg);
  else fprintf(stderr,"== %d == %s: %s\n",getpid(), msg,
               strerror(errno));
  fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
  exit(1);
}

