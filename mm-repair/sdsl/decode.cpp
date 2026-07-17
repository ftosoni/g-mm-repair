/*
 * Function to test the decompression of an int_vector
 */
#include <sdsl/util.hpp>
#include <sdsl/int_vector.hpp>
#include <iostream>
#include <string>
#ifdef MALLOC_COUNT
#include "../tools/malloc_count.h"
#endif



using namespace sdsl;
using namespace std;

int main(int argc, char* argv[])
{
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " iv_file" << endl;
        return 1;
    }

    string ifile = string(argv[1]);
    size_t size;
    {
      int_vector<> v;
      load_from_file(v, ifile);
      size = v.size();
      cout<<"v.size()="<<v.size()<<endl;
      cout<<"v.width()="<<(int)v.width()<<endl;
      cout<<"v[0]="<<v[0]<<endl;
      long sum=0;
      for (size_t i=0; i<size; i++) 
          sum+=v[i];
          // TODO: save the decoded array to file
      cout<<"sum="<<sum<<endl;
    }
  #ifdef MALLOC_COUNT
    // we write %'zu to get the thousands separators  
    fprintf(stderr,"Peak memory allocation: %'zu kbs, %.4lf bytes/entry\n",
           malloc_count_peak()/1024, (double)malloc_count_peak()/(size));
    fprintf(stderr,"Current memory allocation: %'zu bytes\n", malloc_count_current());
  #endif
    
}
