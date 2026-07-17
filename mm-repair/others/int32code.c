/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   int32code.c
   Estimate cost of gamma and 7x8 encoding of a sequence of integers
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>

void estimate(FILE *f,size_t s);

int main(int argc, char *argv[])
{
  void count(FILE *, int);
  extern char *optarg;
  extern int optind, opterr, optopt;
  char *fnam=NULL;
  FILE *f;
  int c,num_opt;
  size_t size = 4;

  /* ------------- read options from command line ----------- */
  num_opt = opterr = 0;
  while ((c=getopt(argc, argv, "s:f:")) != -1) {
    switch (c) 
      {
      case 's':
        size=atoi(optarg); break;
      case 'f':
        fnam = optarg;  break;
      case '?':
        fprintf(stderr,"Unknown option: %c\n", optopt);
        exit(1);
      }
    num_opt++;
  }
 
 if(optind!=argc) {
    fprintf(stderr, "Usage:  %s [-s n] [-f filename]\n",argv[0]);
    fprintf(stderr,"\t-f n        file to process, def: stdin\n\n");    
    fprintf(stderr,"\t-s n        size of each uint, def: 4\n\n");    
    return 0;
  }
  if(size<1 || size>8) {
    fprintf(stderr,"Invalid size parameter: must be between 1 and 8\n");
    exit(EXIT_FAILURE);
  }
    

  /* -------- open input file ------------- */
  if(fnam!=NULL) {
    if (! (f=fopen(fnam, "rb"))) {
      perror(fnam);
      return 1;
    }
  }
  else {
    f = stdout;
  }

  estimate(f,size);

  if(fnam!=NULL) fclose(f);
   
  return 0;
}

int bitsize(unsigned long n)
{
  assert(n>0);
  int tot=0;
  while(n>0) {
    tot += 1;
    n = n/2;
  }
  return tot;
}

int v7x8enc(unsigned long n)
{
  if(n==0) return 1;
  int size = bitsize(n);
  return (size + 6)/7;
}
  
int gammaenc(unsigned long n)
{
  assert(n>0);
  int size = bitsize(n);
  return size + size-1;
}
  

void estimate(FILE *f,size_t s)
{
  unsigned long b;
  
  size_t gamma=0;
  size_t v7x8 = 0;
  size_t z,u=0;
  while( (z=fread(&b,s,1,f))==1) {
    gamma += gammaenc(b+1);
    v7x8 += v7x8enc(b);
    u +=1;
  }
  if(z!=EOF) 
    fprintf(stderr,"Unexpected trailing bytes\n");

  // --- print stats
  printf("# of input symbols:  %zd\n",u);
  printf("Gamma code: %zu bytes (%.2f bits x int)\n",(gamma+7)/8,1.0*gamma/u);
  printf("v7x8 code: %zu bytes (%.2f bits x int)\n",v7x8,8.0*v7x8/u);
}




