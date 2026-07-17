/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 * ReMult
 * 
 * matrix multiplication using a repair compressed matrix
 * first very primitive prototype
 * No longer tested or maintained DO NOT USE!
 * 
 * Copyright (C) 2021-2099   giovanni.manzini@uniupo.it
 * >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>

#ifndef INT_VALS 
typedef float   matval;  // type representing a matrix entry   
typedef double xmatval;  // type representing a matrix entry with larger precision   
#else
typedef int matval;     // type representing a matrix entry   
typedef int xmatval;    // type representing a matrix entry with larger precision   
#endif


// some orrible global variables

int rows,cols;  // # rows and columns of input matrix 
int Alpha;      // alphabet size of input matrix representation/smallest non terminal symbol 
int NTnum;      // number of non terminals 
matval *NTval;  // values associated to non-terminals 

matval *Xval;    // input vector
//matval *Yval;    // output vector 

int Mnum;        // number of distinct non zero matrix values
matval *Mval;    // set of distinct nonzero matrix values

int Debug=0;



// write error message and exit
static void die(const char *s)
{
  perror(s);
  exit(1);
}    

// read a set of matval values from file f
// return number of items in *n and pointer to array with values
static matval *read_vals(FILE *f, int *n)
{
  // get files size
  if(fseek(f,0,SEEK_END))
    die("Error in read_vals:fseek");
  long size = ftell(f);
  if(size<0)
    die("Error in read_vals:ftell");
  // get array size  
  *n = (int) (size/sizeof(matval));
  matval *a  = malloc(*n * sizeof(matval));
  if(a==NULL)
    die("Error in read_vals:malloc");
  rewind(f);
  int e = fread(a,sizeof(matval),*n,f);
  if(e!= *n)
    die("Error in read_vals:fread");
  return a;
}

// decode a terminal representing a matrix entry
// the matrix value is multiplied by the corresponding X entry 
xmatval decode_entry(int p)
{
  int pcol = p % cols;
  int pval = p/cols;
  if(pval>=Mnum) die("Illegal value reference found in terminal");
  return ((xmatval) Xval[pcol])*Mval[pval];
}  

static void fill_NTval(FILE *f) 
{
  int pair[2]; // rule right hand size 
  for(int i=0; i<NTnum;i++) {  // i is number of NT computed so far
    int e = fread(pair,sizeof(int),2,f);
    if(e<2) die("Error reading rule file");
    xmatval sum = 0;
    if(Debug) fprintf(stderr,"%d -- ", i); //!!!!!!!!!!
    for(int j=0;j<2;j++) {
      int p = pair[j];
      if(p>=Alpha) { // non terminal
        p -= Alpha;
        if(p>=i) die("Fatal Error: Forward rule");
        sum += NTval[p];
        if(Debug) fprintf(stderr,"nt:%d  ",p);//!!!!!!!!!!!!
      }
      else { // terminal symbol
        if(p<1) {
          if(Debug) fprintf(stderr,"sep: %d  ",p);//!!!!!!!!
          die("Unique row separator found in rule");
        }
        sum += decode_entry(p-1);
        #ifndef INT_VALS 
        if(Debug) fprintf(stderr,"t: col:%d val:%f ",(p-1)%cols,Mval[(p-1)/cols]);//!!!!!!!111
        #else
        if(Debug) fprintf(stderr,"t: col:%d val:%d ",(p-1)%cols,Mval[(p-1)/cols]);//!!!!!!!111
        #endif
      }
    }
    NTval[i]=sum;
    #ifndef INT_VALS 
    if(Debug) fprintf(stderr,"NT[%d]: %f\n",i,sum); //!!!!!!!!!
    #else 
    if(Debug) fprintf(stderr,"NT[%d]: %d\n",i,sum); //!!!!!!!!!
    #endif
  }
  if(fread(pair,sizeof(int),2,f)>0) die("Unexpect trailing rule");
}



int main (int argc, char **argv) { 
   char fname[PATH_MAX];
   FILE *f,*Rf,*Cf;
   int i,len,e;
   struct stat s;

   // ----------- check input
   fputs("==== Command line:\n",stderr);
   for(int i=0;i<argc;i++)
     fprintf(stderr," %s",argv[i]);
   fputs("\n",stderr);     
   if (argc != 6) { 
     fprintf (stderr,"Usage:\n\t %s matrix rows cols invector outvector \n",argv[0]);
     exit(1);
   }
   // ----------- read and check # rows and cols 
   rows  = atoi(argv[2]);
   if(rows<1) die("Invalid number of rows");
   cols  = atoi(argv[3]);
   if(cols<1) die("Invalid number of columns");
   // ------------ read input vector
   f = fopen(argv[4],"rb");
   if(f==NULL) die("Cannot open input vector file");
   Xval = read_vals(f,&i);
   if(i!=cols) die("Input vector size should be equal to # of columns");
   fclose(f);
   // ------------ read matrix values 
   strcpy(fname,argv[1]);
   strcat(fname,".val");
   f = fopen(fname,"rb");
   if(f==NULL) die("Cannot open matrix values file");
   Mval = read_vals(f,&Mnum);
   
   // ------------ open rules file
   strcpy(fname,argv[1]);
   strcat(fname,".il.R");
   if (stat (fname,&s) != 0) { 
     fprintf (stderr,"Error: cannot stat file %s\n",fname);
     exit(1);
   }
   len = s.st_size;
   Rf = fopen (fname,"rb");
   if (Rf == NULL) die("Cannot open rule file"); 
   // read alphabet size for the original input string 
   if (fread(&Alpha,sizeof(int),1,Rf) != 1) { 
     fprintf (stderr,"Error: cannot read file %s\n",fname);
     exit(1);
   }
   NTnum = (len-sizeof(int))/(2*sizeof(int)); // number of non terminal (rules) 
   NTval =  (matval *) malloc(NTnum*sizeof(matval));  // values of non terminal
   fill_NTval(Rf);
   fclose(Rf);

   // --- open C file
   strcpy(fname,argv[1]);
   strcat(fname,".il.C");
   if (stat (fname,&s) != 0) { 
     fprintf (stderr,"Error: cannot stat file %s\n",fname);
     exit(1);
   }
   Cf = fopen (fname,"r");
   if (Cf == NULL) die("Cannot open C file");

   // --- open output vector file 
   f = fopen (argv[5],"w");
   if (f == NULL) die("Cannot open output vector file");
   // --- compute output 
   int ywritten = 0;
   xmatval sum=0;
   while(ywritten<rows) {
     e = fread(&i,sizeof(int),1,Cf);
     if(e!=1) die("Error reading C file");
     if(i>=Alpha) { // non terminal 
       if( (i = i-Alpha)>= NTnum ) die("Illegal non terminal in C file");
       sum += NTval[i];
     }
     else if(i>0) {// terminal representing a matrix entry
       sum += decode_entry(i-1);
     }
     else { // row completed
       // if(i!=ywritten) die("Incorrect end of row separator");
       matval tmp=sum; sum=0;
       e = fwrite(&tmp,sizeof(matval),1,f);
       if(e!=1) die("Error writing to the output file");
       if(++ywritten>=rows) break;
     }
  }
  assert(ywritten==rows);
  assert(sum==0);
  if(fclose(f)!=0) die("Cannot close output file");
  if(fread(&i,sizeof(int),1,Cf)>0) die("Unexpected trailing symbols in C");
  if(fclose(Cf)!=0) die("Cannot close C file");
  exit(0);
}

