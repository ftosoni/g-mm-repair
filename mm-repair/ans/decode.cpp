// Licensed to the Apache Software Foundation (ASF) under one
// or more contributor license agreements.  See the NOTICE file
// distributed with this work for additional information
// regarding copyright ownership.  The ASF licenses this file
// to you under the Apache License, Version 2.0 (the
// "License"); you may not use this file except in compliance
// with the License.  You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing,
// software distributed under the License is distributed on an
// "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, either express or implied.  See the License for the
// specific language governing permissions and limitations
// under the License.
#include <cassert>
#include <iostream>
#include <vector>
#ifdef MALLOC_COUNT
#include "../tools/malloc_count.h"
#endif

#include "util.hpp"
#include "decode.hpp"



// decode a file produced by encode. create a file with extension .dec
// fold parameter must match the one used for encoding (default is 1)

// number of int32 decoded at each call of decode() 
const size_t Buf_size = 100000;  // must be multiple of 4 


template <uint32_t fidelity>
void run_ansf(std::string input_name)
{

    // read compressed data file 
    std::vector<uint8_t> in_u8;    
    in_u8 = read_file_u8(input_name);
    // size of compressed data
    size_t cSrcSize = in_u8.size()-sizeof(size_t);
    // retrieve decompressed file size
    size_t original_size = *((size_t *) in_u8.data() );
    std::cerr << "# uint32_t in original file: " << original_size << std::endl;
    
    // create decoder
    auto ans_dec = ANSf_decoder<fidelity>();
    ans_dec.init(in_u8.data()+sizeof(size_t), cSrcSize, original_size);
    // decode and save to file
    FILE *outfile = fopen_or_fail(input_name + ".dec", "w");
    uint32_t *buf = new uint32_t[Buf_size];
    size_t decoded = 0;
    while(decoded < original_size) {
      // std::cerr << decoded << std::endl;
      size_t d =  ans_dec.decode(buf,std::min(Buf_size,original_size-decoded));
      if(d==0) quit("Illegal decode result");
      int e = fwrite(buf,sizeof(*buf),d,outfile);
      if(e!=d) quit("Error writing to output file");
      decoded += d;
      assert(decoded <= original_size);
    }
    fclose_or_fail(outfile);
    delete[] buf;
}


int main(int argc, char const* argv[])
{
    if(argc!=3 && argc!=2) {
      std::cerr << "Usage:\n\t"<< argv[0] << " infile [fidelity, def=1]\n";
      exit(1);
    }
    int fidelity = 1; 
    if(argc==3) fidelity = atoi(argv[2]);
    if(fidelity<1 || fidelity>5) {
      std::cerr << "fidelity must be between 1 and 5\n";
      exit(1);
    }
    
    switch(fidelity) {
      case 1:  run_ansf<1>(argv[1]); break;
      case 2:  run_ansf<2>(argv[1]); break;
      case 3:  run_ansf<3>(argv[1]); break;
      case 4:  run_ansf<4>(argv[1]); break;
      case 5:  run_ansf<5>(argv[1]); break;
      default: std::cerr<<"Invalid parameter: " << fidelity << std::endl; exit(3);
    }
    
    #ifdef MALLOC_COUNT
    fprintf(stderr,"Peak memory allocation: %zu bytes, current: %zu bytes\n", 
            malloc_count_peak(),malloc_count_current());
    #endif
        
    return EXIT_SUCCESS;
}
