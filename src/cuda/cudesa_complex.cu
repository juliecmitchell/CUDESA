/*********************************************************************** 
 * cudesa_complex.cu - CUDA host file for cudesaating complexes without
 *                     restrictions on ligand/receptor properties
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

#include "../cudesa_dome.h"
#include "../cudesa_transform.h"
#include "../cudesa_util.h"
#include "cutil.h"

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>


// number of threads per block
#ifdef _DEBUG
#define BLOCK_GRID_HEIGHT 2
#define BLOCK_GRID_WIDTH 2
#else
#define BLOCK_GRID_HEIGHT 8
#define BLOCK_GRID_WIDTH 8
#endif

/* // TODO rename these to avoid confusion (e.g. THREAD_BLOCK_HEIGHT) */

#define BLOCK_SIZE (BLOCK_GRID_HEIGHT * BLOCK_GRID_WIDTH)

/* #define DBG_SIZE 128 */

/* #define INTERNAL_ONLY 1 */

#define BOGUS_ORIGIN_COORD_X 1000000000.0f
#define BOGUS_ORIGIN_COORD_Y 1000000000.0f
#define BOGUS_ORIGIN_COORD_Z 1000000000.0f

#include <cudesa_transform_kernel.cu>
#include <cudesa_molecule_kernel.cu>
#include <cudesa_complex_kernel.cu>

#define ALIGN_INT( n, align ) ( ((n) % (align) == 0) ? (n) : (((n) - ((n) % ((align))) + (align))) )

//#define _DEBUG 3

//#define MULTI_THREAD

//
// NUM_THREADS - Number of GPU threads to spawn, each thread requires
//               one GPU card
//
#ifdef MULTI_THREAD
#define NUM_THREADS 2
#else
#define NUM_THREADS 1
#endif

//
// cudesa_step_arg_t
//
// This structure is passed to each thread working on cudesa
//
struct cudesa_complex_step_arg_t
{
  // receptor - receptor molecule
  PQRData_t* receptor;

  // ligand - ligand molecule
  PQRData_t* ligand;

  // transformations - linked list of transformations to perform
  protein_transformation_t* transformations;

  // id - id of a particular thread
  int id;

  // stride - stride to use when evaluating transformations
  int stride;
};

//
// calculate_cudesa_step_gpu() - Performs cudesa calculations on one GPU
//
static void* calculate_complex_cudesa_step_gpu( void* arg )
{
  cudesa_complex_step_arg_t* thread_arg = (cudesa_complex_step_arg_t*)arg;
  PQRData_t* receptor = thread_arg->receptor;
  PQRData_t* ligand = thread_arg->ligand;
  double cudesa_energy;
  
  protein_transformation_t* transform = thread_arg->transformations;
  int i;
  
  // Bind this thread to a CUDA device
  int deviceCount = 0;
  
  cudaGetDeviceCount( &deviceCount );

  if( thread_arg->id >= deviceCount )
  {
    cudesa_error( "%s(): Not enough CUDA devices for threads! Devices: %d requested thread: %d\n", __FUNCTION__, deviceCount, thread_arg->id );
    return NULL;
  }

  if( NULL == transform )
  {
    cudesa_print( CUDESA_PRINT_ALWAYS, "%s(): [thread %d] dummy thread exiting\n", __FUNCTION__, thread_arg->id );
    return NULL;
  }
  cudesa_print( CUDESA_PRINT_ALWAYS, "%s(): [thread %d] Setting CUDA device to %d\n", __FUNCTION__, thread_arg->id, thread_arg->id );

  cudaSetDevice( thread_arg->id );

  
  dim3 thread_block( BLOCK_GRID_HEIGHT, BLOCK_GRID_WIDTH );

  int ligand_aligned_atoms = ALIGN_INT( ligand->natoms, BLOCK_SIZE );
  int ligand_num_blocks = ligand_aligned_atoms / BLOCK_SIZE;

  int receptor_aligned_atoms = ALIGN_INT( receptor->natoms, BLOCK_SIZE );
  int receptor_num_blocks = receptor_aligned_atoms / BLOCK_SIZE;
 
  int complex_aligned_atoms = receptor_aligned_atoms + ligand_aligned_atoms;
  int complex_num_blocks = complex_aligned_atoms / BLOCK_SIZE;
  
  dim3 complex_grid( complex_num_blocks, 1 );
  dim3 ligand_grid( ligand_num_blocks, 1 );
  dim3 receptor_grid( receptor_num_blocks, 1 );
  
  cudesa_print( CUDESA_PRINT_DEBUG, "%s(): [thread %d] Processing complex (ligand: %s, receptor: %s), %d ligand atoms, %d receptor atoms, %d total atoms %d aligned atoms\n",
                __FUNCTION__, thread_arg->id, ligand->filename, receptor->filename, ligand->natoms, receptor->natoms, ligand->natoms + receptor->natoms, complex_aligned_atoms );
  
  srand( time( NULL ) );

  //
  // set up GPU data
  //
  
  //
  // output data
  //

  float* d_solvent_accessibilities;
  float* h_solvent_accessibilities;

  float* d_base_ligand_accessibilities;
  float* h_base_ligand_accessibilities;
  float* d_base_receptor_accessibilities;
  float* h_base_receptor_accessibilities;
  
  //
  // input data
  //

  // we need to translate our origins to CUDA origins on the CPU
  float3* h_origins = NULL;
  
  // parameters used in SASA computation
  float* d_atomic_radii = NULL;
  float* d_p_i = NULL;
  float* d_S_i = NULL;
  float3* d_origins = NULL;
  float3* d_transform_origins = NULL;

  double surface_area;
  double buried_area;

  double ligand_area;
  double receptor_area;
  
  // device debugging data
#ifdef _DEVICE_DEBUG
  float* d_dbgOut = NULL;
  float* h_dbgOut = NULL;

  CUDA_SAFE_CALL( cudaMalloc( (void**)&d_dbgOut, sizeof(float) * ligand_aligned_atoms ) );
  h_dbgOut = (float*)malloc(sizeof(float) * totalAtoms );
#endif

  //
  // allocate memory
  //
  
  CUDA_SAFE_CALL( cudaMalloc( (void**)&d_solvent_accessibilities, sizeof(float) * complex_aligned_atoms ) );

  h_solvent_accessibilities = (float*)malloc( sizeof(float) * complex_aligned_atoms );

  CUDA_SAFE_CALL( cudaMalloc( (void**)&d_base_ligand_accessibilities, sizeof(float) * ligand_aligned_atoms ) );
  h_base_ligand_accessibilities = (float*)malloc( sizeof(float) * ligand_aligned_atoms );

  CUDA_SAFE_CALL( cudaMalloc( (void**)&d_base_receptor_accessibilities, sizeof(float) * receptor_aligned_atoms ) );
  h_base_receptor_accessibilities = (float*)malloc( sizeof(float) * receptor_aligned_atoms );
  
  h_origins = (float3*)malloc( sizeof(float3) * complex_aligned_atoms );
  
  CUDA_SAFE_CALL( cudaMalloc( (void**)&d_atomic_radii, sizeof(float) * complex_aligned_atoms ) );
  CUDA_SAFE_CALL( cudaMalloc( (void**)&d_p_i, sizeof(float) * complex_aligned_atoms ) );
  CUDA_SAFE_CALL( cudaMalloc( (void**)&d_S_i, sizeof(float) * complex_aligned_atoms ) );
  CUDA_SAFE_CALL( cudaMalloc( (void**)&d_origins, sizeof(float3) * complex_aligned_atoms ) );
  CUDA_SAFE_CALL( cudaMalloc( (void**)&d_transform_origins, sizeof(float3) * complex_aligned_atoms ) );

#ifdef _DEVICE_DEBUG
  for( i = 0; i < ligand_aligned_atoms; i++ )
  {
    h_dbgOut[ i ] = -1.0f;
  }
  cudaMemcpy( d_dbgOut, h_dbgOut, sizeof(float) * ligand_aligned_atoms, cudaMemcpyHostToDevice );
#endif

  // initialize origins
  for( int i = 0; i < complex_aligned_atoms; i++ )
  {
    if( i < ligand->natoms )
    {
      if( ligand->is_hydrogen[ i ] )
      {
        h_origins[ i ].x = BOGUS_ORIGIN_COORD_X;
        h_origins[ i ].y = BOGUS_ORIGIN_COORD_Y;
        h_origins[ i ].z = BOGUS_ORIGIN_COORD_Z;
      }
      else
      {
        h_origins[ i ].x = ligand->atoms[ i ].xyz[ 0 ];
        h_origins[ i ].y = ligand->atoms[ i ].xyz[ 1 ];
        h_origins[ i ].z = ligand->atoms[ i ].xyz[ 2 ];
      }
    }
    else if( i < ligand_aligned_atoms )
    {
      h_origins[ i ].x = BOGUS_ORIGIN_COORD_X;
      h_origins[ i ].y = BOGUS_ORIGIN_COORD_Y;
      h_origins[ i ].z = BOGUS_ORIGIN_COORD_Z;
    }
    else if( (i - ligand_aligned_atoms) < receptor->natoms )
    {
      if( receptor->is_hydrogen[ i - ligand_aligned_atoms ] )
      {
        h_origins[ i ].x = BOGUS_ORIGIN_COORD_X;
        h_origins[ i ].y = BOGUS_ORIGIN_COORD_Y;
        h_origins[ i ].z = BOGUS_ORIGIN_COORD_Z;
      }
      else
      {
        h_origins[ i ].x = receptor->atoms[ i - ligand_aligned_atoms ].xyz[ 0 ];
        h_origins[ i ].y = receptor->atoms[ i - ligand_aligned_atoms ].xyz[ 1 ];
        h_origins[ i ].z = receptor->atoms[ i - ligand_aligned_atoms ].xyz[ 2 ];
      }
    }
    else
    {
      h_origins[ i ].x = BOGUS_ORIGIN_COORD_X;
      h_origins[ i ].y = BOGUS_ORIGIN_COORD_Y;
      h_origins[ i ].z = BOGUS_ORIGIN_COORD_Z;
    }
  }

  cudaMemset( d_solvent_accessibilities, 0, sizeof(float) * complex_aligned_atoms );
  cudaMemset( d_base_ligand_accessibilities, 0, sizeof(float) * ligand_aligned_atoms );
  cudaMemset( d_base_receptor_accessibilities, 0, sizeof(float) * receptor_aligned_atoms );
  
  // copy all the data to the graphics card
  // copy origins
  CUDA_SAFE_CALL( cudaMemcpy( d_origins, h_origins, sizeof(float3) * complex_aligned_atoms, cudaMemcpyHostToDevice ) );
  CUDA_SAFE_CALL( cudaMemcpy( d_transform_origins, h_origins, sizeof(float3) * complex_aligned_atoms, cudaMemcpyHostToDevice ) );
  
  // copy atomic radii
  cudaMemset( d_atomic_radii, 0, sizeof(float) * complex_aligned_atoms );
  CUDA_SAFE_CALL( cudaMemcpy( d_atomic_radii, ligand->atomicRadii, sizeof(float) * ligand->natoms, cudaMemcpyHostToDevice ) );
  CUDA_SAFE_CALL( cudaMemcpy( d_atomic_radii + ligand_aligned_atoms, receptor->atomicRadii, sizeof(float) * receptor->natoms, cudaMemcpyHostToDevice ) );

  // copy p_i parameters
  CUDA_SAFE_CALL( cudaMemcpy( d_p_i, ligand->p_iParams, sizeof(float) * ligand->natoms, cudaMemcpyHostToDevice ) );
  CUDA_SAFE_CALL( cudaMemcpy( d_p_i + ligand_aligned_atoms, receptor->p_iParams, sizeof(float) * receptor->natoms, cudaMemcpyHostToDevice ) );  
  
  // copy S_i parameters
  CUDA_SAFE_CALL( cudaMemcpy( d_S_i, ligand->S_i, sizeof(float) * ligand->natoms, cudaMemcpyHostToDevice ) );
  CUDA_SAFE_CALL( cudaMemcpy( d_S_i + ligand_aligned_atoms, receptor->S_i, sizeof(float) * receptor->natoms, cudaMemcpyHostToDevice ) );

  // Calculate the base solvation area
  calculate_solvent_accessibility_molecule_kernel<<<receptor_grid, thread_block>>>( d_base_receptor_accessibilities,
#ifdef _DEVICE_DEBUG
                                                                                    dbgOut,
#endif
                                                                                    0,
                                                                                    receptor_num_blocks,
                                                                                    receptor->natoms,
                                                                                    d_atomic_radii + ligand_aligned_atoms,
                                                                                    d_p_i + ligand_aligned_atoms,
                                                                                    d_S_i + ligand_aligned_atoms,
                                                                                    d_origins + ligand_aligned_atoms );

  calculate_solvent_accessibility_molecule_kernel<<<ligand_grid, thread_block>>>( d_base_ligand_accessibilities,
#ifdef _DEVICE_DEBUG
                                                                                  dbgOut,
#endif
                                                                                  0,
                                                                                  ligand_num_blocks,
                                                                                  ligand->natoms,
                                                                                  d_atomic_radii,
                                                                                  d_p_i,
                                                                                  d_S_i,
                                                                                  d_origins );

  CUDA_SAFE_CALL( cudaMemcpy( h_base_receptor_accessibilities, d_base_receptor_accessibilities, sizeof(float) * receptor_aligned_atoms, cudaMemcpyDeviceToHost ) );
  CUDA_SAFE_CALL( cudaMemcpy( h_base_ligand_accessibilities, d_base_ligand_accessibilities, sizeof(float) * ligand_aligned_atoms, cudaMemcpyDeviceToHost ) );  

  //
  // Process each transformation
  //

  float3* h_transform_origins = (float3*)malloc( sizeof(float3) * complex_aligned_atoms );
  
  while( transform )
  {
    buried_area = 0.0;
    surface_area = 0.0;
    ligand_area = 0.0;
    receptor_area = 0.0;
    
    // Load this tarnsformations rotation and translation
    CUDA_SAFE_CALL( cudaMemcpyToSymbol( rotation, transform->single_rotation, sizeof(float) * 9 ) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol( translation, transform->single_translation, sizeof(float) * 3 ) );

    // Transform the ligand origins
    transform_origins<<<complex_grid, thread_block>>>( d_origins, d_transform_origins, ligand->natoms );

    CUDA_SAFE_CALL( cudaMemcpy( h_transform_origins, d_transform_origins, sizeof(float3) * complex_aligned_atoms, cudaMemcpyDeviceToHost ) );
    
    // Calculate complex accessibility contributions
    calculate_solvent_accessibility_complex_kernel<<<complex_grid, thread_block>>>( d_solvent_accessibilities,
#ifdef _DEVICE_DEBUG
                                                                                    dbgOut,
#endif
                                                                                    complex_num_blocks,
                                                                                    complex_aligned_atoms,
                                                                                    d_atomic_radii,
                                                                                    d_p_i,
                                                                                    d_S_i,
                                                                                    d_transform_origins );
    
    // Retrieve cudesa contributions
    CUDA_SAFE_CALL( cudaMemcpy( h_solvent_accessibilities, d_solvent_accessibilities, sizeof(float) * complex_aligned_atoms, cudaMemcpyDeviceToHost ) );  

#ifdef _DEVICE_DEBUG
    CUDA_SAFE_CALL( cudaMemcpy( h_dbgOut, d_dbgOut, sizeof(float) * totalAtoms, cudaMemcpyDeviceToHost ) );

    cudesa_print( CUDESA_PRINT_ALWAYS, "\n\nDevice debug output:\n\n" );
    for( i = 0; i < totalAtoms; i++ )
    {
      if( h_dbgOut[ i ] > -0.5f )
      {
        cudesa_print( CUDESA_PRINT_ALWAYS, "%f\n", h_dbgOut[ i ] );
      }
    }    
#endif

    // sum cudesa
    cudesa_energy = 0.0;
    //    cudesa_print( CUDESA_PRINT_ALWAYS, "GPU delta A's:\n" );
    for( i = 0; i < complex_aligned_atoms; i++ )
    {
      if( i < ligand->natoms )
      {
/*         if( ligand->is_hydrogen[ i ] ) */
/*         { */
/*           continue; */
/*         } */
        cudesa_energy += (h_solvent_accessibilities[ i ] - h_base_ligand_accessibilities[ i ]) * ligand->asp[ i ];
        //cudesa_print( CUDESA_PRINT_ALWAYS, "ligand delta a[%d]: %f (%d)\n", i, fabs(h_solvent_accessibilities[ i ] - h_base_ligand_accessibilities[ i ]) < 0.00001 ? 0.0 : (h_solvent_accessibilities[ i ] - h_base_ligand_accessibilities[ i ]), ligand->is_hydrogen[ i ] );
        buried_area += (h_solvent_accessibilities[ i ] - h_base_ligand_accessibilities[ i ]);
        surface_area += h_solvent_accessibilities[ i ];
	ligand_area += h_base_ligand_accessibilities[ i ];
      }
      else if( i < ligand_aligned_atoms )
      {
        i = ligand_aligned_atoms - 1;
        continue;
      }
      else if( (i - ligand_aligned_atoms) < receptor->natoms )
      {
/*         if( receptor->is_hydrogen[ i - ligand_aligned_atoms ] ) */
/*         { */
/*           continue; */
/*         } */
        cudesa_energy += (h_solvent_accessibilities[ i ] - h_base_receptor_accessibilities[ i - ligand_aligned_atoms ]) * receptor->asp[ i - ligand_aligned_atoms ];
        //cudesa_print( CUDESA_PRINT_ALWAYS, "receptor delta a[%d]: %f (%d)\n", i - ligand_aligned_atoms, fabs(h_solvent_accessibilities[ i ] - h_base_receptor_accessibilities[ i - ligand_aligned_atoms ]) < 0.00001 ? 0.0 : (h_solvent_accessibilities[ i ] - h_base_receptor_accessibilities[ i - ligand_aligned_atoms ]), receptor->is_hydrogen[ i - ligand_aligned_atoms ] );        
        buried_area += (h_solvent_accessibilities[ i ] - h_base_receptor_accessibilities[ i - ligand_aligned_atoms ]);
        surface_area += h_solvent_accessibilities[ i ];
	receptor_area += h_base_receptor_accessibilities[ i - ligand_aligned_atoms ];
      }
      else
      {
        break;
      }
    }

    transform->calculated_gpu_complex_energy = cudesa_energy;
    transform->calculated_gpu_buried_area = buried_area;
    transform->calculated_gpu_sasa = surface_area;
    transform->calculated_gpu_ligand_sasa = ligand_area;
    transform->calculated_gpu_receptor_sasa = receptor_area;
    //cudesa_print( CUDESA_PRINT_ALWAYS, "GPU cudesa: %f\n", cudesa_energy );
    for( int i = 0; i < thread_arg->stride; i++ )
    {
      if( transform )
      {
        transform = transform->next;
      }
    }
  }

  cudaFree( d_solvent_accessibilities );
  cudaFree( d_base_ligand_accessibilities );
  cudaFree( d_base_receptor_accessibilities );
  cudaFree( d_atomic_radii );
  cudaFree( d_p_i );
  cudaFree( d_S_i );
  cudaFree( d_origins );
  cudaFree( d_transform_origins );

  free( h_solvent_accessibilities );
  free( h_base_ligand_accessibilities );
  free( h_base_receptor_accessibilities );
  free( h_origins );
  free( h_transform_origins );
#ifdef _DEVICE_DEBUG
  cudaFree( d_dbgOut );
  free( h_dbgOut );
#endif

  
  return NULL;
} 

//
// calculate_complex_cudesa_gpu() - Main entrypoint to compute
//                                  desolvation energy using the
//                                  CUDESA GPU algorithm
//
// input                          - receptor: receptor molecule
//                                  ligand: ligand molecule
//                                  transformations: list of
//                                  transformations for which to
//                                  compute desolvation energies
//
// output                         - output energies are stored in the
//                                  transformations array
//
void calculate_complex_cudesa_gpu( PQRData_t* receptor, PQRData_t* ligand, protein_transformation_t* transformations )
{
  cudesa_complex_step_arg_t thread_args[ NUM_THREADS ];
  pthread_t threads[ NUM_THREADS ];
  for( int i = 0; i < NUM_THREADS; i++ )
  {
    thread_args[ i ].receptor = receptor;
    thread_args[ i ].ligand = ligand;
    thread_args[ i ].transformations = transformations;

    for( int j = 0; j < i; j++ )
    {
      if( thread_args[ i ].transformations )
      {
        thread_args[ i ].transformations = thread_args[ i ].transformations->next;
      }
      else
      {
        cudesa_error( "%s(): Multi-threading stride is greater than number of transformations, creating dummy thread\n", __FUNCTION__ );
        thread_args[ i ].transformations = NULL;
        break;
      }
    }

    thread_args[ i ].id = i;
    thread_args[ i ].stride = NUM_THREADS;
  
    pthread_create( &threads[ i ], NULL, calculate_complex_cudesa_step_gpu, (void*)&thread_args[ i ] );
  }

  for( int i = 0; i < NUM_THREADS; i++ )
  {
    pthread_join( threads[ i ], NULL );
  }
}
