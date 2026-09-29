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
# nvcc: $(CUDA_HOME)/bin/nvcc if CUDA_HOME is set, else the one on PATH, else /usr/local/cuda.
# Override on the command line, e.g. `make NVCC=/opt/cuda/bin/nvcc` or `make CUDA_HOME=/opt/cuda`.
ifdef CUDA_HOME
NVCC ?= $(CUDA_HOME)/bin/nvcc
else
NVCC ?= $(or $(shell command -v nvcc 2>/dev/null),/usr/local/cuda/bin/nvcc)
endif
# GPU compute capability; the paper's GB10 is sm_121. Use e.g. CUDA_ARCH=sm_90 (Hopper),
# sm_89 (Ada), sm_80 (Ampere), or CUDA_ARCH=native to target the GPU of the build machine.
CUDA_ARCH ?= sm_121
NVCCFLAGS = -O3 -std=c++17 -arch=$(CUDA_ARCH) -Xcompiler -Wall -Xcompiler -fopenmp -lgomp

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