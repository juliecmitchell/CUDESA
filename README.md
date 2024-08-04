# CUDESA
This is release 1.0 of the CUDESA GPU desolvation code. This work
appeared in the April 2009 issue of the Journal of Computational
Biology. Please acknowledge your use of this work by citing:

David Dynerman, Erick Butzlaff, Julie C. Mitchell. Journal of
Computational Biology. April 2009, 16(4):
523-537. doi:10.1089/cmb.2008.0157.

0. License
==========
CUDESA is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

CUDESA is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

1. Information
==============
The source code is contained in src/, and the CUDA-specific bits are
in src/cuda/. 

This code was compiled and tested on RHEL 5, with gcc version 4.1.2
and CUDA version 2.1. The code also compiles on Mac OS X.  

2. Compilation
==============
To compile the program, just run 

	make osx
	make linux
	
according to your architecture.

If you type 
	
	make debug
	
both gcc and nvcc will compile with debugging symbols, and
nvcc will compile for device emulation, NOT hardware. In this
compilation stage, _CUDESA_DEBUG is defined to be 1, so you can use
that as an #ifdef gate for debug-only code.

3. Running
==========
By default, binaries are placed in bin/. To run cudesa, type

	bin/cudesa

and examine the options. At minimum, you want to specify a ligand and
receptor PQR file on which to perform the desolvation. Some sample PQR
files are in sample/

	bin/cudesa -l sample/1CHO_l.pqr -r sample/1CHO_r.pqr

By default CUDESA runs the CPU algorithm. Use '-g' to run on the GPU
and '-b' to run both versions and compare.

4. PQR Files
============
PQR Files were generated using pdb2pqr:

	http://pdb2pqr.sourceforge.net/

Dolinsky TJ, Nielsen JE, McCammon JA, Baker NA. PDB2PQR: an automated
pipeline for the setup, execution, and analysis of Poisson-Boltzmann
electrostatics calculations. Nucleic Acids Research, 32, W665-W667
(2004).

5. Timing
=========
The '-t' flag can be used to run a timing analysis. By default this
will time how long it takes to compute desolvation energy for 2000
random ligand transformations. Specifying both '-b' (both CPU/GPU) and
'-t' (timing) will run the timing analysis on both the CPU and GPU
algorithm and print some comparison information.

Additionally, the '-x' flag can be used as a multiplier to increase
the number of transformations. For example, '-x 3' will perform 6000
ligand transformations, and so on.

On Linux both the CPU and GPU algorithms compile for 2 threads by
default. To increase the number of CPU threads used, modify the
NUM_THREADS define in src/cudesa.c. To increase the number of GPU
threads used (each GPU thread uses one CUDA card) modify NUM_THREADS
in src/cuda/cudesa_complex.cu

On OS X both the CPU and GPU algorithms are single thread/single GPU
due to limitations in the OS X version of pthreads.


=====================================================
