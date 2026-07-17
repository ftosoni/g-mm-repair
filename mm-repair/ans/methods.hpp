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

#pragma once


#include "ans_byte.hpp"
#include "ans_fold.hpp"
//~ #include "ans_int.hpp"
//~ #include "ans_msb.hpp"
#include "ans_reorder_fold.hpp"

//~ #include "ans_sint.hpp"
//~ #include "ans_smsb.hpp"


template <uint32_t fidelity> struct ANSfold {
    static const uint32_t fold_fidelity = fidelity;
  
    static std::string name()
    {
        return std::string("ANSfold-") + std::to_string(fidelity);
    }

    static size_t encode(const uint32_t* in_ptr, size_t in_size_u32,
        uint8_t* out_ptr, size_t out_size_u8, uint8_t* buf = NULL)
    {
        return ans_fold_compress<fidelity>(
            out_ptr, out_size_u8, in_ptr, in_size_u32);
    }
    static void decode(const uint8_t* in_ptr, size_t in_size_u8,
        uint32_t* out_ptr, size_t out_size_u32, uint8_t* buf = NULL)
    {
        ans_fold_decompress<fidelity>(
            out_ptr, out_size_u32, in_ptr, in_size_u8);
    }
};

template <uint32_t fidelity> struct ANSrfold {
    static std::string name()
    {
        return std::string("ANSrfold-") + std::to_string(fidelity);
    }

    static size_t encode(const uint32_t* in_ptr, size_t in_size_u32,
        uint8_t* out_ptr, size_t out_size_u8, uint8_t* buf = NULL)
    {
        return ans_reorder_fold_compress<fidelity>(
            out_ptr, out_size_u8, in_ptr, in_size_u32);
    }
    static void decode(const uint8_t* in_ptr, size_t in_size_u8,
        uint32_t* out_ptr, size_t out_size_u32, uint8_t* buf = NULL)
    {
        ans_reorder_fold_decompress<fidelity>(
            out_ptr, out_size_u32, in_ptr, in_size_u8);
    }
};


struct entropy_only {
    static std::string name() { return std::string("entropy"); }

    static size_t encode(const uint32_t* in_ptr, size_t in_size_u32,
        uint8_t* out_ptr, size_t out_size_u8, uint8_t* buf = NULL)
    {
        auto [input_entropy, sigma] = compute_entropy(in_ptr, in_size_u32);
        return ceil((input_entropy * in_size_u32) / 8.0);
    }
    static void decode(const uint8_t* in_ptr, size_t in_size_u8,
        uint32_t* out_ptr, size_t out_size_u32, uint8_t* buf = NULL)
    {
        std::cerr << "don't do this!!" << std::endl;
        exit(EXIT_FAILURE);
    }
};
