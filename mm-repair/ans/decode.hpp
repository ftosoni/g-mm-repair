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

#include <algorithm> 
#include "methods.hpp"



// main decoding routines: functions init() and decode() 
// together do the work of  ans_fold_decompress() in ans_fold.hpp
template <uint32_t fidelity> struct ANSf_decoder {
    static const uint32_t fold_fidelity = fidelity;

    ans_fold_decode<fidelity> ans_frame;
    const uint8_t *in_u8;    
    size_t in_u8_size;
    std::array<uint64_t, 4> states;
    size_t cur_idx, to_decode; 

    // call this to start from scratch the decoding of a new sequence
    void init(const uint8_t* cSrc, size_t cSrcSize, size_t original_size) {
      in_u8 = cSrc;
      in_u8_size = cSrcSize;
      ans_frame = ans_fold_decode<fidelity>::load(in_u8);
      in_u8 += cSrcSize; // go to the end and proceed backward
      states[3] = ans_frame.init_state(in_u8);
      states[2] = ans_frame.init_state(in_u8);
      states[1] = ans_frame.init_state(in_u8);
      states[0] = ans_frame.init_state(in_u8);
      cur_idx = 0;
      to_decode = original_size;   
    }
    
    // decode up to size uint32_t writing them to dst
    // return number of decoded uint32_t 
    size_t decode(uint32_t* odst, size_t size) {
      assert(size>0);   // don't waste time with empty requests
      // size must be multiple of 4 unless at the end of data stream
      assert(size%4==0 || size+cur_idx==to_decode);
      if(size+cur_idx>to_decode) 
        return 0; // illegal request forces illegal answer
      
      uint32_t *dst = odst;  
      size_t decoded = 0;
      size_t fast_decode = to_decode - (to_decode % 4);
      while (cur_idx != fast_decode) { // decode 4 ints per iteration 
        *dst++ = ans_frame.decode_sym(states[3], in_u8);
        *dst++ = ans_frame.decode_sym(states[2], in_u8);
        *dst++ = ans_frame.decode_sym(states[1], in_u8);
        *dst++ = ans_frame.decode_sym(states[0], in_u8);
        cur_idx += 4;
        decoded += 4;
        if(decoded>=size) {
          assert(decoded == size); // because of previous assert
          return decoded;
        }
      }
      // because of previous assert the last <4 bytes are read in a single shot
      assert(size-decoded == to_decode-cur_idx);
      while (cur_idx != to_decode) {  // last < 4 ints 
        *dst++ = ans_frame.decode_sym(states[0], in_u8);
        cur_idx++;
        decoded++;
      }
      if(decoded!=size) quit("Something is very wrong in decoding");
      return decoded;
    }
};

