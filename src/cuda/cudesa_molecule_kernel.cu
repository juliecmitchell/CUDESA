/*********************************************************************** 
 * cudesa_molecule.cu - CUDA device code for intra-molecular cudesa
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
// calculate_solvent_accessibility_molecule_kernel() - Calculate the
//                                                     solvent
//                                                     accessiblity
//                                                     for a given
//                                                     molecule
//
// input -
//   start                     - which receptor atom to start processing with,
//                               for parallel kernels running on multiple GPU's
//   num_molecule_blocks       - number of BLOCK_SIZE-size blocks of atoms
//                               to process
//   num_molecule_atoms        - number of molecule atoms to process
//  
//   The following are parallel arrays indexed on atom # in the range [0, n]
//   where n is num_ligand_atoms or num_receptor_atoms aligned to the next
//   highest multiple of BLOCK_SIZE
// 
//   d_atomic_radii            - molecule atomic radii
//   d_p_i                     - molecule p_i parameters
//   d_S_i                     - molecule S_i parameters
//   d_origin                  - molecule atom position
//
//
// output - 
//
//   d_solvent_accessibilities - solvent accessibilities for molecule
//   dbgOut                    - optional debugging output, device code can
//                               put intermediate values, etc in here
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
                                                                 )
{
  // atom offset that this thread will process within each block
  int atom_offset = IMUL( threadIdx.y, blockDim.x ) + threadIdx.x;
  
  // block index within the grid
  int block_index = IMUL( blockIdx.y, gridDim.x ) + blockIdx.x;

  // overall molecule atom index that this thread will process
  int molecule_atom_id = IMUL( block_index, BLOCK_SIZE ) + atom_offset + start; 

  // which atom data we'll load each loop iteration
  int molecule_block_atom_id;
  
#ifdef _DEVICE_DEBUG  
  int dbg_counter = 0;
#endif

  // this atom's accessibility
  float molecule_solvent_accessibility = 1.0;

  // loop counters for looping over the receptor atoms
  int molecule_grid_index;
  int molecule_block_index;

  // local register cache for the ligand atom we're processing
  float atomic_radius, p_i, S_i;
  float3 origin;

  // partial solvation access term
  float f_ij;

  // Load the local data for the atom we're computing
  atomic_radius = d_atomic_radii[ molecule_atom_id ];
  p_i           = d_p_i[ molecule_atom_id ];
  S_i           = d_S_i[ molecule_atom_id ];
  origin        = d_origin[ molecule_atom_id ];

  // iterate over receptor blocks
  for( molecule_grid_index = 0; molecule_grid_index < num_molecule_blocks; molecule_grid_index++ )
  {
    __syncthreads();
    molecule_block_atom_id = IMUL( molecule_grid_index, BLOCK_SIZE ) + atom_offset;
    
    // Load this receptor atom's data into shared memory
    s_atomic_radii[ atom_offset ] = d_atomic_radii[ molecule_block_atom_id ];
    s_p_i[ atom_offset ]          = d_p_i[ molecule_block_atom_id ];
    s_S_i[ atom_offset ]          = d_S_i[ molecule_block_atom_id ];
    s_origin[ atom_offset ]       = d_origin[ molecule_block_atom_id ];
    
    __syncthreads();

    // compare ligand atom against all receptor atoms in this block
    for( molecule_block_index = 0; molecule_block_index < BLOCK_SIZE; molecule_block_index++ )
    {
      // This check isn't usually neccesary - bogus atoms beyond the real number are placed at
      // BOGUS_ORIGIN_COORD_X, _Y, _Z which should almost always cause the distance check below to fail
      // In _DEBUG this fact is checked by code in cudesa.cu, and if it fails (e.g. real atoms are
      // close to BOGUS_ORIGIN_COORD_X, _Y, _Z) then an error is printed suggesting use of _CHECK_BOUNDS
      int molecule_global_index = IMUL( molecule_grid_index, BLOCK_SIZE ) + molecule_block_index;

      if( molecule_global_index == molecule_atom_id )
      {
        continue;
      }      
      
#ifdef _CHECK_BOUNDS      
      if( molecule_global_index >= num_molecule_atoms )
      {
        break;
      }
#endif
      float distance = sqrt( (s_origin[ molecule_block_index ].x - origin.x) * (s_origin[ molecule_block_index ].x - origin.x) +
                             (s_origin[ molecule_block_index ].y - origin.y) * (s_origin[ molecule_block_index ].y - origin.y) +
                             (s_origin[ molecule_block_index ].z - origin.z) * (s_origin[ molecule_block_index ].z - origin.z) );

      if( molecule_atom_id == 0 && molecule_block_index + molecule_grid_index * BLOCK_SIZE == 1934 )
      {
        //printf( "1934 radius: %f + %f\n", s_atomic_radii[ molecule_block_index ], atomic_radius );
      }
      
      if( distance < ( s_atomic_radii[ molecule_block_index ] + atomic_radius + (2.0 * SOLVATION_RADIUS) ) )
      {
        if( distance < 1.8 )
        {
          f_ij = 1.0 - (0.8875 *
                        p_i *
                        d_compute_b( atomic_radius, s_atomic_radii[ molecule_block_index ], distance, SOLVATION_RADIUS) )
                        / S_i;
        }
        else
        {
          f_ij = 1.0 - (0.3516 *
                        p_i *
                        d_compute_b( atomic_radius, s_atomic_radii[ molecule_block_index ], distance, SOLVATION_RADIUS) )
                        / S_i;
        }
#ifdef _DEBUG
        if( molecule_atom_id == DEBUG_ATOM )
        {
          DEBUG_PRINT( CUDESA_PRINT_DEBUG_EXTRA, "*= %f [p_i: %f, S_i: %f d: %f r_i: %f r_j: %f] [%d->%d] [%d: <%.3f, %.3f, %.3f> %d: <%.3f, %.3f, %.3f>]\n", f_ij, p_i, S_i, distance, atomic_radius,
                       s_atomic_radii[ molecule_block_index ], molecule_block_index + molecule_grid_index * BLOCK_SIZE, molecule_atom_id, molecule_block_index + molecule_grid_index * BLOCK_SIZE, s_origin[ molecule_block_index ].x, s_origin[ molecule_block_index ].y, s_origin[ molecule_block_index ].z, molecule_global_index, origin.x, origin.y, origin.z );
        }
#endif


#ifdef _DEVICE_DEBUG
        if( molecule_atom_id == DEBUG_ATOM )
        {
          dbgOut[ (dbgCounter++) + block_index ] = f_ij;
        }
#endif

        molecule_solvent_accessibility *= f_ij;
      }
#ifdef _DEBUG
      else
      {
        if( molecule_atom_id == DEBUG_ATOM )
        {
          DEBUG_PRINT( CUDESA_PRINT_DEBUG_SPEW, "Rejected [%d->%d] [d=%f cutoff=%f] [%d: <%.3f, %.3f, %.3f> %d: <%.3f, %.3f, %.3f>]\n",
                       molecule_block_index + molecule_grid_index * BLOCK_SIZE, molecule_atom_id, distance,
                       s_atomic_radii[ molecule_block_index ] + atomic_radius + (2.0 * SOLVATION_RADIUS),
                       molecule_block_index + molecule_grid_index * BLOCK_SIZE, s_origin[ molecule_block_index ].x,
                       s_origin[ molecule_block_index ].y, s_origin[ molecule_block_index ].z, molecule_global_index,
                       origin.x, origin.y, origin.z );
        }
      }
#endif
    }
  }

#if _DEBUG
  if( molecule_atom_id == DEBUG_ATOM )
  {
    DEBUG_PRINT( CUDESA_PRINT_DEBUG_EXTRA, "Storing accessibility %f for atom %d [%.3f * %.3f]\n", molecule_solvent_accessibility * S_i,
                 molecule_atom_id, molecule_solvent_accessibility, S_i );
  }
#endif
  
  d_solvent_accessibilities[ molecule_atom_id - start ] = molecule_solvent_accessibility * S_i;
}
