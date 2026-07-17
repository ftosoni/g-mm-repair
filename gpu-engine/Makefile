# This file is part of g-mm-repair
# <https://github.com/ftosoni/g-mm-repair>.
# Copyright (c) 2026 Francesco Tosoni.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

CXX = g++
CXXFLAGS = -O3 -std=c++17 -Wall -fopenmp
NVCC = /usr/local/cuda/bin/nvcc
NVCCFLAGS = -O3 -std=c++17 -arch=sm_121 -Xcompiler -Wall -Xcompiler -fopenmp -lgomp

# Shared host-side source used by every driver (CUDA-free): the schedule builder.
HOST_SRC = build_schedule.cpp
HOST_HDR = grammar.h build_schedule.h layered_types.h

all: host_scheduler gpu_test cpu_selftest cusparse_test

host_scheduler: host_scheduler.cpp
	$(CXX) $(CXXFLAGS) -o host_scheduler host_scheduler.cpp

# Unified GPU driver for the RePair construction the paper evaluates.
gpu_test: gpu_engine.cu gpu_engine_test.cu gpu_engine.h $(HOST_SRC) $(HOST_HDR)
	$(NVCC) $(NVCCFLAGS) -o gpu_test gpu_engine.cu gpu_engine_test.cu $(HOST_SRC) -lnvidia-ml

# Unified CUDA-free correctness test for the RePair pipeline; this is what CI runs.
cpu_selftest: cpu_selftest.cpp $(HOST_SRC) $(HOST_HDR)
	$(CXX) $(CXXFLAGS) -o cpu_selftest cpu_selftest.cpp $(HOST_SRC)

# cuSPARSE baseline benchmark
cusparse_test: cusparse_test.cu $(HOST_SRC)
	$(NVCC) $(NVCCFLAGS) -o cusparse_test cusparse_test.cu $(HOST_SRC) -lcusparse -lnvidia-ml

clean:
	rm -f host_scheduler gpu_test cpu_selftest cusparse_test