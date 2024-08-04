#
# cudesa - GPU CUDESA Makefile
#
# For use with GNU Make
#
# Copyright (C) 2007 David Dynerman, Julie Mitchell
#
# Visit the Mitchell Lab: http://www.mitchell-lab.org
#
# Contact information: dynerman@cs.wisc.edu, mitchell@math.wisc.edu
#
#	This file is part of CUDESA
#
# CUDESA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CUDESA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with CUDESA  If not, see <http://www.gnu.org/licenses/>.

# CUDA_INSTALL_PATH
#
# Point this to the base of your CUDA installation
CUDA_INSTALL_PATH = /usr/local/cuda

# BINARY
#
# What we're creating
BINARY = bin/cudesa

# BINARY_VERSION
#
# What version we're building
BINARY_VERSION = 1.0

# SRC
#
# What we're compiling
C_SRC = src/cudesa_transform.c src/cudesa_dome.c src/cudesa_util.c src/cudesa.c src/main.c src/vector.c
CU_SRC = src/cuda/cudesa_complex.cu #src/cuda/cudesa.cu
C_OBJECTS = $(addsuffix .o, $(basename $(C_SRC)))
CU_OBJECTS = $(addsuffix .o, $(basename $(CU_SRC)))

#
# Compiler specific sections
#

# GCC
GCC = g++

# GCC_DEBUG_CFLAGS
#
# Debug mode CFLAGS for g++
GCC_DEBUG_CFLAGS = -g -Wall -Wshadow -D_DEBUG=1 -D _cudesa_DEBUG=1  

# GCC_CFLAGS
#
GCC_LINUX_CFLAGS = -fstack-protector -O2 -ffast-math -fno-exceptions -fomit-frame-pointer -fno-math-errno 
GCC_OSX_CFLAGS = -fstack-protector -O2 -ffast-math -fno-exceptions -fomit-frame-pointer -fno-math-errno -malign-double

# GCC_LFLAGS
#
# Linker flags for g++
GCC_LFLAGS = -L$(CUDA_INSTALL_PATH)/lib -lcuda -lcudart -lm -lpthread

# NVCC
NVCC = /usr/local/cuda/bin/nvcc

# NVCC_DEBUG_CFLAGS
#
# Debug mode CFLAGS for nvcc
NVCC_DEBUG_CFLAGS = -g -deviceemu -D _DEBUG=1 -D _cudesa_DEBUG=1

# NVCC_CFLAGS
#
# Release mode CFLAGS for nvcc 
NVCC_LINUX_CFLAGS = -O2 -DMULTI_THREAD  
NVCC_OSX_CFLAGS = -O2 -use_fast_math  

C_COMPILER = $(GCC)

# Targets

# linux 
#
# compile all targets using gcc for linux (assumes 32 bit, pthreads)
linux:
	make build C_COMPILER="$(GCC)" C_LFLAGS="$(GCC_LFLAGS)" C_CFLAGS="$(GCC_LINUX_CFLAGS)" CU_CFLAGS="$(NVCC_LINUX_CFLAGS)" BINARY_VERSION="$(BINARY_VERSION) gcc linux"

# osx 
#
# compile all targets using gcc for osx (assumes 64 bit, no pthreads)
osx:
	make build C_COMPILER="$(GCC)" C_LFLAGS="$(GCC_LFLAGS)" C_CFLAGS="$(GCC_OSX_CFLAGS)" CU_CFLAGS="$(NVCC_OSX_CFLAGS)" BINARY_VERSION="$(BINARY_VERSION) gcc osx"

# debug
#
# compile all targets using gcc in debug mode (for linux)
debug:
	make build C_COMPILER="$(GCC)" C_LFLAGS="$(GCC_LFLAGS)" C_CFLAGS="$(GCC_DEBUG_CFLAGS)" CU_CFLAGS="$(NVCC_DEBUG_CFLAGS)" BINARY_VERSION="$(BINARY_VERSION) gcc debug"

build: $(BINARY)

$(BINARY):  $(C_OBJECTS) $(CU_OBJECTS)
	$(C_COMPILER) $(C_LFLAGS) $(C_OBJECTS) $(CU_OBJECTS) -o $@

clean:
	rm -rf $(C_OBJECTS) $(CU_OBJECTS) $(BINARY)

%.o : %.c
	$(C_COMPILER) $(C_CFLAGS) -DBINARY_VERSION="\"$(BINARY_VERSION)\"" -DBINARY_NAME="\"$(notdir $(BINARY))\"" -c $< -o $@

%.o : %.cu
	$(NVCC) $(CU_CFLAGS) -I $(@D)/ -c $< -o $@
