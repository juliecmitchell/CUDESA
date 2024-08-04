/*********************************************************************** 
 * cudesa_common_kernel.cu - CUDA device code common for all kernels
 *
 * Copyright (C) 2007-2009 Regents of the University of Wisconsin
 * 
 * Authors: David Dynerman and Julie Mitchell
 *
 * Visit the Mitchell Lab: http://www.mitchell-lab.org
 * Contact information: admin@mitchell-lab.org
 *
 * This file is part of CUDESA.  If you use CUDESA as an application or
 * adapt this code per the terms of the GNU license below, this license 
 * requires you to cite the following publication:
 *
 *     David Dynerman, Erick Butzlaff, Julie C. Mitchell (2009) "CUSA and CUDE: 
 *     GPU-Accelerated Methods for Estimating Solvent Accessible Surface Area 
 *     and Desolvation," Journal of Computational Biology 16(4): 523-537.
 *
 * CUDESA is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CUDESA is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with cudesa.  If not, see <http://www.gnu.org/licenses/>.
 *
 **********************************************************************/

#ifndef _CUDESA_COMMON_KERNEL_CU_
#define _CUDESA_COMMON_KERNEL_CU_

typedef unsigned int uint;

//
// compute_b() - This helper function just computes the b_ij term for
//               the SASA computation. See Equation (7) in the paper.
// input       - r_i, r_j: atomic radii for the two atoms in question
//               d_ij: distance between atom i and atom j
//               r_s: solvent probe radius
// output      - b_ij SASA factor as per Equation (7)
//
__device__ float d_compute_b( float r_i, float r_j, float d_ij, float r_s )
{
  return M_PI*(r_i + r_s)*(r_i + r_j + 2.0*r_s - d_ij)*(1.0 + (r_j - r_i)/d_ij);
}

//
// shared memory arrays for block-based iteration, as described in
// Section (2.2.3)
//
__shared__ float s_atomic_radii[ BLOCK_SIZE ];
__shared__ float s_p_i[ BLOCK_SIZE ];
__shared__ float s_S_i[ BLOCK_SIZE ];
__shared__ float3 s_origin[ BLOCK_SIZE ];

// rotation and translation matrices for cudesa_transform_kernel.cu
__constant__ float rotation[ 3 ][ 3 ];
__constant__ float translation[ 3 ];

//24-bit multiplication is faster on G80,
//but we must be sure to multiply integers
//only within [-8M, 8M - 1] range
#define IMUL(a, b) __mul24(a, b)

#define sl_scale 3.0

// for cudesa_print, etc
#ifdef _DEBUG
#include "../cudesa_util.h"
#endif

//
// calculate_solvent_accessibility_molecule_kernel() - GPU kernel definition
//
__global__ void calculate_solvent_accessibility_molecule_kernel( float* d_solvent_accessibilities,
#ifdef _DEVICE_DEBUG
                                                                 float* dbgOut,
#endif
                                                                 uint start,
                                                                 uint num_molecule_blocks,
                                                                 uint num_molecule_atoms,
                                                                 float* d_atomic_radii,
                                                                 float* d_p_i,
                                                                 float* d_S_i,
                                                                 float3* d_origin
                                                                 );


#endif

