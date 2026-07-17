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

// encode a file of uint32_t using ans-fold k (default k=1)
// create a file with extension ansf.k
// the first 8 bytes of the outfile is the number of input uint32_t


#include <cassert>
#include <iostream>
#include <vector>
#ifdef MALLOC_COUNT
#include "../tools/malloc_count.h"
#endif

#include "util.hpp"
#include "methods.hpp"

// number of times we compress (for computing average speed) 
const int NUM_RUNS = 1;


template <class t_compressor>
void run(std::vector<uint32_t>& input, std::string input_name)
{
  // (0) compute entropy
  auto [input_entropy, sigma] = compute_entropy(input);

  // (1) encode
  std::vector<uint8_t> encoded_data(input.size() * 8 * 2);
  std::vector<uint8_t> tmp_buf(input.size() * 8 * 2);

  size_t encoded_bytes;
  size_t encoding_time_ns_min = std::numeric_limits<size_t>::max();
  for (int i = 0; i < NUM_RUNS; i++) {
    auto start_encode = std::chrono::high_resolution_clock::now();
    encoded_bytes = t_compressor::encode(input.data(), input.size(),
        encoded_data.data(), encoded_data.size(), tmp_buf.data());
    auto stop_encode = std::chrono::high_resolution_clock::now();
    auto encoding_time_ns = stop_encode - start_encode;
    encoding_time_ns_min
      = std::min((size_t)encoding_time_ns.count(), encoding_time_ns_min);
  }

  double BPI = double(encoded_bytes * 8) / double(input.size());
  double encode_IPS = compute_ips(input.size(), encoding_time_ns_min);
  double enc_ns_per_int = double(encoding_time_ns_min) / double(input.size());

  // (4) output stats
  printf("(%s, obytes: %zd, BPI: %2.4f, iint32: %zd)\n", t_compressor::name().c_str(), encoded_bytes, BPI,input.size());
  fflush(stdout);

  // save file
  encoded_data.resize(encoded_bytes);
  std:: string outfile = input_name + ".ansf." + std::to_string(t_compressor::fold_fidelity);
  auto fd = fopen_or_fail(outfile, "w");
  // write input size in sizeof(size_t) bytes  
  size_t isize = input.size(); 
  size_t e = fwrite(&isize,sizeof(isize),1,fd);
  if(e!=1) quit("Error writing input file size");    
  // write compressed data
  e = fwrite(encoded_data.data(),1,encoded_bytes, fd); 
  if(e!=encoded_bytes) quit("Error writing compressed data");
  fclose_or_fail(fd);
}

int main(int argc, char const* argv[])
{
  if(argc!=3 && argc!=2) {
    std::cerr << "Usage:\n\t"<< argv[0] << " file_of_u32 [fidelity def=1] \n";
    exit(1);
  }

  int fidelity = 1; 
  if(argc==3) fidelity = atoi(argv[2]);
  if(fidelity<1 || fidelity>5) {
    std::cerr << "fidelity must be between 1 and 5\n";
    exit(1);
  }
  // read file 
  std::vector<uint32_t> input_u32s;
  input_u32s = read_file_u32(argv[1]);
  switch(fidelity) {
    case 1:  run<ANSfold<1>>(input_u32s, argv[1]); break;
    case 2:  run<ANSfold<2>>(input_u32s, argv[1]); break;
    case 3:  run<ANSfold<3>>(input_u32s, argv[1]); break;
    case 4:  run<ANSfold<4>>(input_u32s, argv[1]); break;
    case 5:  run<ANSfold<5>>(input_u32s, argv[1]); break;
    default: std::cerr<<"Invalid parameter: " << fidelity << std::endl; exit(3);
  }

#ifdef MALLOC_COUNT
  fprintf(stderr,"Peak memory allocation: %zu bytes, current: %zu bytes\n", 
      malloc_count_peak(),malloc_count_current());
#endif

  return EXIT_SUCCESS;
}
