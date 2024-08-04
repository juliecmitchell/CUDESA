/*********************************************************************** 
 * cudesa_transform_kernel.cu - CUDA device code for transforming
 *                              atom origins
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

#include <cudesa_common_kernel.cu>

//
// transform_origins() - Transforms the given list of origins by the
//                       transformation given in the arrays rotation
//                       and translation(defined in
//                       cudesa_common_kernel.cu)
//
// input               - d_origins: input origins to transform
//                       num_origins: size of input array
//
// output              - d_output_origins: transformed origins
//
__global__ void transform_origins( float3* d_origins, float3* d_output_origins, int num_origins )
{
  // which atom we process in a given block
  int atom_offset = IMUL( threadIdx.y, blockDim.x ) + threadIdx.x;
  
  // block index within the grid
  int block_index = IMUL( blockIdx.y, gridDim.x ) + blockIdx.x;

  // overall ligand atom index that this thread will process
  int atom_id = IMUL( block_index, BLOCK_SIZE ) + atom_offset;

  if( atom_id >= num_origins )
  {
    return;
  }
  
  float3 temp_origin = d_origins[ atom_id ];
  float3 out_origin;
  
  // transform ligand origin based on matrix
  out_origin.x = translation[ 0 ];
  out_origin.x += (rotation[ 0 ][ 0 ] * temp_origin.x) + (rotation[ 0 ][ 1 ] * temp_origin.y) + (rotation[ 0 ][ 2 ] * temp_origin.z);
  
  out_origin.y = translation[ 1 ];
  out_origin.y += (rotation[ 1 ][ 0 ] * temp_origin.x) + (rotation[ 1 ][ 1 ] * temp_origin.y) + (rotation[ 1 ][ 2 ] * temp_origin.z);
  
  out_origin.z = translation[ 2 ];
  out_origin.z += (rotation[ 2 ][ 0 ] * temp_origin.x) + (rotation[ 2 ][ 1 ] * temp_origin.y) + (rotation[ 2 ][ 2 ] * temp_origin.z);
  
  d_output_origins[ atom_id ] = out_origin;
}
  
