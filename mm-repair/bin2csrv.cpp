/*  bin2csrv[if]
    (here and in the following [if] stands for either nothing, i or f)
    Convert a binary file contaning int32/float/double entries into the 
      CSRV format (.vc & .[if]val files) 
    or the 
      DRV  format (.dv & .[if]vald files)
    In the CSRV format we only store nonzero entries and each entry
    is represented by an id identifying the value in the .[if]val file
    and the column number (hence the .vc extension)
    In the DRV format we store zero and nonzero entries and each entry is
    represented by an id identifying it in the .[if]vald file 
    (hence the .dv extension) 
    In both formats id's are different from zero and the zero value is used 
    to mark the end of a matrix row 

    This source file is compiled in 3 versions:
      bin2csrv  (default),     matrix entries are 8 byte floats (eg doubles) 
      bin2csrvi (Typecode==1), matrix entries are 4 byte ints 
      bin2csrvf (Typecode==2), matrix entries are 4 byte floats 
    By default the output values stored in the [if]val file are in the same 
    type as in the input file, if option -d is used, output values are 
    stored as doubles. The extension [if]val is chosen according to the 
    type of the file values.  
    The option -c can be used to specify that the input values are complex
    with the real and immaginary parts either int/float or double.
    The type of the matrix entries does not affect the format of the 
    .vc/dv file  
    */
#include <cassert>
#include <iostream>
#include <unordered_map>
#include <map>
#include <cstring>
#include <limits.h>
#include <unistd.h>
#include <time.h>
#include <complex>
#include <array>
#include <cstdint>

#if Typecode==1  // int32_t    
#define Type int32_t
#define valext ".ival"
#define scnf "%d"
#elif Typecode==2 // float
#define Type float
#define valext ".fval"
#define scnf "%f"
#else            // default type is double 
#define Type double
#define valext ".val"
#define scnf "%lf"
#endif


// bit representing types of input/output values 
// and general options of the algorithm 
//#define INT32_OUTPUT 1
//#define FLOAT_OUTPUT 2
#define COMPLEX_INPUT 4
#define NO_COL_ID 8
#define DOUBLE_OUTPUT 16
#define BOOLEAN_INPUT 32
#define SPARSE_INPUT 64


static void usage_and_exit(char *name)
{
    fprintf(stderr,"Usage:\n\t  %s [options] matrix rows cols\n",name);
    fprintf(stderr,"\t\t-b num         number of row blocks, def. 1\n");    
    fprintf(stderr,"\t\t-c             input are complex numbers\n");
    fprintf(stderr,"\t\t-B             input is a boolean matrix in textual sparse format: row col\n");
    fprintf(stderr,"\t\t-S             input is in textual sparse format: row col value\n");
    fprintf(stderr,"\t\t-d             save matrix entries as doubles\n");
    fprintf(stderr,"\t\t-n             don't store col id (drv format, debug only)\n");
    fprintf(stderr,"\t\t-v             verbose\n");
    fprintf(stderr,"The option -c can be combined with -d\n\n");
    exit(1);
}

// write error message and exit
static void quit(const char *s)
{
  if(errno==0) fprintf(stderr,"%s\n",s);
  else  perror(s);
  exit(1);
}

// number of bits to represent longs up to n
int bits (unsigned long n) {
  int b=1;
  while (n>1) { n>>=1; b++; }
  return b;
}


// write a x as an int32/float/double in binary format  
// if there is loss of information a warning is sent to stderr
void write_bin(Type x, int out_type, FILE *f)
{
  int e;
  if(out_type&DOUBLE_OUTPUT) {
    double xd = x;
    e = fwrite(&xd,sizeof(xd),1,f);
    if(e!=1) 
      quit("Error writing to binary .val file (DOUBLE_OUTPUT)");
  }
  else {
    e = fwrite(&x,sizeof(x),1,f);
    if(e!=1) 
      quit("Error writing to binary " valext " file");
  }
}


// provide hash function for a complex<Type> 
struct ComplexHasher
{
  size_t operator()(const std::complex<Type>& k) const
  {
    using std::hash;

    return ((hash<Type>()(k.real())
             ^ (hash<Type>()(k.imag()) << 1)) >> 1);
  }
};



// copy row :rindex from a sparse boolean matrix file :f
// to file :fvc in csrv format working only with the nonzero 
// elements and therefore in time O(nz) (only used for the csrv format)
// In the input file row indexes must be non-decreasing,
// some rows can be empty and there must be no duplicates
// (we cannot check for duplicates in this function)
// It is not necessary that elements are in strict row major order
// however, sorting will usually help grammar compression and it is used to 
// check for duplicates, so this condition should be satisfied as well.
// We take advantage of the fact that the id of the unique nonzero
// (the value 1) is 0, so we do not use the values[] dictionary
// Update :maxcode with the largest code written to the .vc file
// and :nonz with the total number of processed values (they are all nonzeros)
// used for boolean matrices in csrv format 
void copy_sparse_bool_row(int rindex, int cols, FILE *fvc, 
                          unsigned long *maxcode, unsigned long *nonz, FILE *f) {
  // keep track of last pair read
  static int r=-1, c = -1;
  const unsigned long id = 0; // id of the unique nonzero
  while (true) {
    if(r<0) {
      assert(c<0);
      int e = fscanf(f,"%d %d",&r,&c);
      if(e!=2) {
        if(feof(f)) return;    // from now on nothing to read, all elements are 0
        quit("Error reading sparse input file: invalid data format (copy_sparse_bool_row)");
      }
      if(r<0 or c<0) quit("Error! Sparse input file has negative row or column index (copy_sparse_bool_row)");
    }
    // now we have a row and a column (either from previous iteration or from fscanf)
    if(r==rindex) {
      if(c<0 or c>=cols) 
        quit("Error! Sparse input file has column index out of range (copy_sparse_bool_row)");
      unsigned long code;
      *nonz += 1;               // found a new nonzero value
      code = id*cols + c;       // we are using the same instructions as in the general case
      code += 1;                // we could simpy write code = c+1 since id==0 
      if (code>= 1ul<<30) 
        quit("Code larger than 2**30. We are in trouble (copy_sparse_bool_row)");
      if (code>*maxcode) *maxcode=code;
      uint32_t wcode = (uint32_t) code; // convert to 32 bit value
      if(fwrite(&wcode,sizeof(wcode),1,fvc)!=1) 
        quit("Error writing to .vc file (copy_sparse_bool_row)");
      r=c=-1;      // reset r,c pair
    }
    else {
      if(r<rindex) quit("Error! Sparse input file is not in row-major order (copy_sparse_bool_row)");
      return;    // row completed, save r,c pair for next row
    }
  }
}

// as above, but now the text file contains also the values of the nonzero elements
void copy_sparse_row(int rindex, int cols, FILE *fvc, unsigned long *maxcode, unsigned long *nonz, FILE *f, 
                     unsigned long *dnonz, std::unordered_map<Type,unsigned long> &values, FILE *fval, int vtype) 
{ 
  assert(!(vtype&NO_COL_ID));      // this function is only for csrv format
  // keep track of last pair read
  static int r=-1, c = -1;
  static Type v=0;
  while (true) {
    if(r<0) {
      assert(c<0);
      int e = fscanf(f,"%d %d" scnf,&r,&c, &v); // read row/column/value
      if(e!=3) {
        if(feof(f)) return;    // from now on nothing to read, all elements are 0
        quit("Error reading sparse input file: invalid data format (copy_sparse_row)");
      }
      if(r<0 or c<0) quit("Error! Sparse input file has negative row or column index (copy_sparse_row)");
    }
    // now we have a row a column and a value(either from previous iteration or from fscanf)
    if(r==rindex) {
      if(c<0 or c>=cols) 
        quit("Error! Sparse input file has column index out of range (copy_sparse_row)");
      if(v!=0) { // process non zero value
        *nonz += 1;               // found a new nonzero value
        unsigned long id, code;
        if(values.count(v)==0) {
          id = values[v] = (*dnonz)++;
          write_bin(v,vtype,fval);
        } else id = values[v];
        code = id*cols + c;       // we are using the same instructions as in the general case
        code += 1;                // we could simpy write code = c+1 since id==0 
        if (code>= 1ul<<30) 
          quit("Code larger than 2**30. We are in trouble (copy_sparse_bool_row)");
        if (code>*maxcode) *maxcode=code;
        uint32_t wcode = (uint32_t) code; // convert to 32 bit value
        if(fwrite(&wcode,sizeof(wcode),1,fvc)!=1) 
          quit("Error writing to .vc file (copy_sparse_bool_row)");
      }
      r=c=-1;      // reset r,c pair
    }
    else {
      if(r<rindex) quit("Error! Sparse input file is not in row-major order (copy_sparse_row)");
      return;    // row completed, save r,c,v triplet for next row
    }
  }
}


// read row :rindex from a sparse boolean matrix file :f
// the row is initialized to zero and the elements are set to 1
// every time a pair (r,c) is read from the file with r=rindex
// row indexes must be non-decreasing, some rows can be empty;
// for convenience we require that entries are in row-major order
// and without duplicates
// used for boolean matrices in drv format 
void read_sparse_bool_row(int rindex, int cols, Type *row, FILE *f)
{
  // keep track of last pair read
  static int r=-1, c = -1;
  // set all row elements to zero
  for(int i=0;i<cols;i++) row[i]=0;
  while (true) {
    if(r<0) {
      assert(c<0);
      int e = fscanf(f,"%d %d",&r,&c);
      if(e!=2) {
        if(feof(f)) return;    // from now on nothing to read, all elements are 0
        quit("Error reading sparse input file: invalid data format (read_sparse_bool_row)");
      }
      if(r<0 or c<0) quit("Error! Sparse input file has negative row or column index (read_sparse_bool_row)");
    }
    // now we have a row and a column (either from previous iteration or from fscanf)
    if(r==rindex) {
      if(c<0 or c>=cols) quit("Error! Sparse input file has column index out of range (read_sparse_bool_row)");
      if(row[c]!=0) {
        assert(row[c]==1);
        fprintf(stderr,"Warning! Duplicate entry \"%d %d\" in sparse input file\n",r,c);
        quit("Error! Duplicate entry in sparse input file (read_sparse_bool_row)");
      }
      row[c]=1;  // set element to 1  
      r=c=-1;      // reset r,c pair
    }
    else {
      if(r<rindex) quit("Error! Sparse input file is not in row-major order (read_sparse_bool_row)");
      return;    // row completed, save r,c pair for next row
    }
  }
}


int main (int argc, char **argv) { 
  extern char *optarg;
  extern int optind, opterr, optopt;
  int verbose=0,vtype=0; 
  int c,rows,cols,nblocks=1;
  char fname[PATH_MAX];
  time_t start_wc = time(NULL);
  uint32_t wcode;  // int written to the .vc files

  /* ------------- read options from command line ----------- */
  opterr = 0;
  while ((c=getopt(argc, argv, "b:cdnvBS")) != -1) {
    switch (c) 
      {
      case 'v':
        verbose++; break;
      case 'B':
        vtype |= BOOLEAN_INPUT; break;
      case 'S':
        vtype |= SPARSE_INPUT; break;        
      case 'c':
        vtype |= COMPLEX_INPUT; break;
      case 'd':
        vtype |= DOUBLE_OUTPUT; break;
      case 'n':
        vtype |= NO_COL_ID; break;
      case 'b':
        nblocks=atoi(optarg); break;
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
  if(nblocks < 1) {
    fprintf(stderr,"Error! Option -b must be at least one\n");
    usage_and_exit(argv[0]);
  }
  if((vtype & BOOLEAN_INPUT) and (vtype & COMPLEX_INPUT)) {
    fprintf(stderr,"Error! Boolean input cannot be complex\n");
    usage_and_exit(argv[0]);
  }
  if((vtype & SPARSE_INPUT) and (vtype & COMPLEX_INPUT)) {
    fprintf(stderr,"Error! Sparse input cannot be complex\n");
    usage_and_exit(argv[0]);
  }
  if((vtype & SPARSE_INPUT) and (vtype & NO_COL_ID) ) {
    fprintf(stderr,"Error! Sparse input not supported for drv format\n");
    usage_and_exit(argv[0]);
  }


  // virtually get rid of options from the command line 
  optind -=1;
  if (argc-optind != 4) usage_and_exit(argv[0]); 
  argv += optind; argc -= optind;
  
  // ----------- read and check # rows and cols 
  rows  = atoi(argv[2]);
  if(rows<1) quit("Invalid number of rows");
  cols  = atoi(argv[3]);
  if(cols<1) quit("Invalid number of columns");
  size_t rowsize = cols*(vtype&COMPLEX_INPUT?2:1);

  if(nblocks>rows) quit("Too many row blocks!");
  // compute size of each block 
  int block_size = (rows+nblocks-1)/nblocks;
  fprintf(stderr,"Block size: %d\n", block_size);

  // open input file and check number of rows/cols
  FILE *f = fopen(argv[1],"rb");
  if(f==NULL) quit("Cannot open infile");
  if(! ((vtype&BOOLEAN_INPUT)or(vtype&SPARSE_INPUT)) ) {// binary input: check file size 
    if(fseek(f,0,SEEK_END)!=0) quit("Cannot seek input file");
    size_t fsize = ftell(f);
    if(fsize<0) quit("Cannot tell input file size");
    if(fsize!=sizeof(Type)*rowsize*rows) quit("Invalid input file size");
    rewind(f);
  }
  // open output .[if]val[d] file
  strncpy(fname,argv[1],PATH_MAX-10);
  strncat(fname,vtype & DOUBLE_OUTPUT ? ".val" : valext,6);
  if(vtype & NO_COL_ID) strncat(fname,"d",2);
  FILE *fval = fopen(fname,"wb");
  if(fval==NULL) quit("Cannot open values file for writing");

  // init counters
  int  wr = 0;  // number of written rows
  unsigned long    nonz = 0;      // total number of processed values (nonzeros for crsv format)  
  unsigned long   dnonz = 0;      // distinct values in .val file (nonzeros for crsv format)  
  unsigned long maxcode = 0;      // largest code in a .vc/.dv file
  // dictionary for storing input values
  std::unordered_map<Type,unsigned long> values; // dictionary of distinct nonzero
  Type v;   // values read from file
  // same for complex values
  std::unordered_map<std::complex<Type>,unsigned long, ComplexHasher> covalues; // dictionary of distinct nonzero
  std::complex<Type> cov;
  Type re,im;
  // extension of the matrix file .vc or .dv
  const char *mext = (vtype &NO_COL_ID) ? ".dv" : ".vc";
  // array cointaining a single row of the matrix
  Type *row = new Type[rowsize];
  // special case of boolean input and csrv format 
  if(vtype&BOOLEAN_INPUT and !(vtype&NO_COL_ID)) {
    // since we only store 1 in fval and values[] we do it there
    // dnonz values[] and fval will not change
    values[1] = dnonz++;
    write_bin(1,vtype,fval);
  } 


  // main loop reading binary file 
  size_t n=0;
  char *buffer=NULL;
  for(int bn=0;bn<nblocks;bn++) {
    if(nblocks==1) snprintf(fname,PATH_MAX,"%s%s",argv[1],mext);
    else snprintf(fname,PATH_MAX,"%s.%d.%d%s",argv[1],nblocks,bn,mext);
    FILE *fvc = fopen(fname,"wb");
    if(fvc==NULL) quit("Cannot open a .vc/.dv file");
    while(true) {
      // process input matrix one row at a time
      if(vtype&BOOLEAN_INPUT and !(vtype&NO_COL_ID)) { // do not use row array & values dictionary
        // for boolean csrv matrices we read the row from a sparse file
        // and fill the vc file with the column indexes in overall O(nz) time
        copy_sparse_bool_row(wr,cols,fvc,&maxcode,&nonz,f); // copy from sparse file to fvc
      }
      else if(vtype&SPARSE_INPUT) { // do not use row array & values dictionary
        // for sparse matrices we read the row from a sparse file
        // and fill the vc file with the column indexes in time O(nz)
        assert(!(vtype&NO_COL_ID));    // drv matrices not supported
        assert(!(vtype&COMPLEX_INPUT)); // complex matrices not supported
        copy_sparse_row(wr,cols,fvc,&maxcode,&nonz,f,&dnonz,values,fval,vtype); // copy from sparse file to fvc
      }
      else { // copy row from input file to row array
        if(vtype & BOOLEAN_INPUT) {
          assert(vtype&NO_COL_ID);    // csrv matrices handled above
          read_sparse_bool_row(wr,cols,row,f); // copy from sparse file to row array
        }
        else {
          // copy from binary file to row array
          size_t e = fread(row,sizeof(Type),rowsize,f);
          if(e!=rowsize)
            if(ferror(f)) quit("Error reading input file");
            else quit("Error reading input file: unexpected end of file");
        }
        // process row regardless of type bool/real/complex
        for(int c=0;c<cols;c++) {
          if(vtype&COMPLEX_INPUT) {
            cov.real(row[2*c]);
            cov.imag(row[2*c+1]);
            v= (cov.real()==0 and cov.imag()==0) ? 0 : 1; // v==0 iff cov==0
          }
          else v = row[c];
          // process element if nonzero or NO_COL_ID is set
          if(v!=0 or (vtype&NO_COL_ID)) { // process non zero value
            unsigned long id, code;
            nonz += 1;      // note: v could be 0 so nonz is technically the wrong name
            // get id of entry and possibly write new values to .val file
            if(!(vtype&COMPLEX_INPUT)) { // real entry
              if(values.count(v)==0) {
                id = values[v] = dnonz++;
                write_bin(v,vtype,fval);
              } else id = values[v];
            } else { // complex entry 
              if(covalues.count(cov)==0) {
                id = covalues[cov] = dnonz++;
                write_bin(cov.real(),vtype, fval);
                write_bin(cov.imag(),vtype, fval);
              } else id = covalues[cov];
            }
            // generate code and write it to .vc file 
            if (vtype&NO_COL_ID) code = id;
            else code = id*cols + c;
            code += 1;                // shift by 1 to allow code for endrow code 0 
            if (code>= 1ul<<30) 
              quit("Code larger than 2**30. We are in trouble");
            if (code>maxcode) maxcode=code;
            wcode = (uint32_t) code; // convert to 32 bit value
            if(fwrite(&wcode,sizeof(wcode),1,fvc)!=1) 
              quit("Error writing to .vc/.dv file");
          }
        }
      } // end if(vtype&BOOLEAN_INPUT and !(vtype&NO_COL_ID)
      // end row
      wcode=0;
      if(fwrite(&wcode,sizeof(wcode),1,fvc)!=1) 
        quit("Error writing to .vc/.dv file");
      wr += 1;
      // if(wr%1000==0) fprintf(stdout,"Block %d/%d: row %d/%d\n",bn+1,nblocks,wr,rows);
      // if a block is full, stop and start the next 
      if (wr>=rows or ((wr%block_size) == 0) )
        break; // end while true
    }
    fclose(fvc);
  }
  delete[] row; // deallocate row array
  fclose(fval);
  fclose(f);
  if(vtype&COMPLEX_INPUT) assert(dnonz==covalues.size());
  else assert(dnonz==values.size());

  fprintf(stderr,"Elapsed time: %.0lf secs\n",(double) (time(NULL)-start_wc));  
  fprintf(stderr,"Number of stored values: %ld   Stored ratio: %.4f\n", nonz, (nonz/((double) wr*cols)));  
  fprintf(stderr, "%zd distinct values in .[if]val file (nonzeros in crsv format) \n", dnonz);
  fprintf(stderr,"Largest codeword: %lu   Bits x codeword: %d\n", maxcode, bits(maxcode));
  fprintf(stderr,"==== Done\n");
  
  return EXIT_SUCCESS;
}
