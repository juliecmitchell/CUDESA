/*********************************************************************** 
 * cudesa.c - CPU cudesa kernel
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

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <string.h>

#include "vector.h"
#include "cudesa_dome.h"
#include "cudesa_transform.h"
#include "cudesa_util.h"
#include "cudesa.h"

//
// cpu_cudesa_step_arg_t - Threading information for the CPU CUDESA algorithm.
//                         Passed as argument to each CUDESA thread
//
struct cpu_cudesa_step_arg_t
{
  // receptor - receptor molecule
  PQRData_t* receptor;

  // ligand - ligand molecule
  PQRData_t* ligand;

  // transformations - list of transformations that we're running
  //                   CUDESA on, each CUDESA CPU thread will operate
  //                   on some subset of this list
  protein_transformation_t* transformations;

  // stride - thread stride, determines which subset of
  //          transformations a particular CUDESA CPU thread will
  //          operate on
  int stride;
};

//
// NUM_CPU_THREADS - Number of threads to use for the CPU CUDESA algorithm
//
#define NUM_CPU_THREADS 2

//
// compute_b() - This helper function just computes the b_ij term for
//               the SASA computation. See Equation (7) in the paper.
// input       - r_i, r_j: atomic radii for the two atoms in question
//               d_ij: distance between atom i and atom j
//               r_s: solvent probe radius
// output      - b_ij SASA factor as per Equation (7)
//
inline double compute_b( double r_i, double r_j, double d_ij, double r_s )
{
  return M_PI*(r_i + r_s)*(r_i + r_j + 2.0*r_s - d_ij)*(1.0 + (r_j - r_i)/d_ij);
}

//
// calculate_intramolecular_sasa() - Computes SASA contributions on
//                                   each atom in molecule from other
//                                   atoms in molecule
// input                           - molecule: input molecule
//
//                                   molecule_origins: atomic
//                                   positions, use this information
//                                   rather than the positions in
//                                   molecule since we may be applying
//                                   various transformations
//
// output                          - sasa: array containing SASA
//                                   values in sq. ang
//
// NOTE: All arrays are parallel arrays indexed accoring to the atom
//       list in molecule
//
void calculate_intramolecular_sasa( PQRData_t* molecule, vec3_t* molecule_origins, double* sasa )
{
  int i, i_res, ifirst, ilast;
	int j, j_res, jfirst, jlast;
	int k;
	double d, f_ij, p_i, dA_i[3], A_i, S_i, r_s = 1.4;
  
	for( i_res = 0; i_res < molecule->nres; i_res++)
  {
    ifirst = molecule->residues[i_res].firstatom;
    if( i_res == (molecule->nres) - 1 )
    {
      ilast = molecule->natoms;
    }
    else
    {
      ilast = molecule->residues[i_res + 1].firstatom;
    }

    for( i = ifirst; i < ilast; i++)
    {
      A_i = 1.0;

      for( k = 0; k < 3; k++ )
      {
        dA_i[k] = 0.0;
      }
		  
		  S_i = molecule->atoms[i].r + r_s;
		  S_i = 4.0 * M_PI * S_i * S_i;
		  p_i = molecule->p_iParams[ i ]; 
		  
		  for( j_res = 0; j_res < molecule->nres; j_res++)
      {
        jfirst = molecule->residues[j_res].firstatom;
        if( j_res == (molecule->nres) - 1 )
        {
          jlast = molecule->natoms;
        }
        else
        {
          jlast = molecule->residues[j_res + 1].firstatom;
        }

        for( j = jfirst; j < jlast; j++ )
        {
          if( j == i )
          {
            continue;
          }
          
			    d = 0.0;
          
			    for( k = 0; k < 3; k++ )
          {
            d += (molecule_origins[ i ][ k ] - molecule_origins[ j ][ k ]) * (molecule_origins[ i ][ k ] - molecule_origins[ j ][ k ]);
          }
          d = sqrt(d);
			    
			    if( d < ( molecule->atoms[j].r + molecule->atoms[i].r + 2.0 * r_s ) )
          {
            if( d < 1.8 )
            {
              f_ij = 1.0 - (0.8875 * p_i * compute_b( molecule->atoms[i].r, molecule->atoms[j].r, d, r_s) ) / S_i;
            }
            else
            {
              f_ij = 1.0 - (0.3516 * p_i * compute_b( molecule->atoms[i].r, molecule->atoms[j].r, d, r_s) ) / S_i;
            }
           
            A_i *= f_ij;
			    }
        } //end j loop
		  } //end j_res loop

		  sasa[ i ] = S_i * A_i;
    } //end i loop
	} //end i_res loop
}

//
// calculate_intermolecular_sasa() - Calculates surface area
//                                   contributions on molecule_obs due
//                                   to molecule_src
// input                           - molecule_obs: the molecule whose
//                                   atoms we're computing SASA on
//
//                                   molecule_src: the molecule
//                                   contributing SASA effects to
//                                   atoms in molecule_obs
//
//                                   obs_origins: atomic positions for
//                                   atoms in molecule_obs
//
//                                   src_origins: atomic positions for
//                                   atoms in molecule_src
//
// output                          - sasa: computed SASA contributions
//
void calculate_intermolecular_sasa( PQRData_t* molecule_obs, PQRData_t* molecule_src, vec3_t* obs_origins, vec3_t* src_origins, double* sasa )
{
  int i, i_res, ifirst, ilast;
	int j, j_res, jfirst, jlast;
	int k;
	double d, f_ij, p_i, A_i, S_i, r_s = 1.4;
  
	for( i_res = 0; i_res < molecule_obs->nres; i_res++)
  {
    ifirst = molecule_obs->residues[i_res].firstatom;
    if( i_res == (molecule_obs->nres) - 1 )
    {
      ilast = molecule_obs->natoms;
    }
    else
    {
      ilast = molecule_obs->residues[i_res + 1].firstatom;
    }

    for( i = ifirst; i < ilast; i++)
    {
		A_i = 1.0;
		  
		  S_i = molecule_obs->atoms[i].r + r_s;
		  S_i = 4.0 * M_PI * S_i * S_i;
		  p_i = molecule_obs->p_iParams[ i ];  
		  
		  for( j_res = 0; j_res < molecule_src->nres; j_res++)
      {
        jfirst = molecule_src->residues[j_res].firstatom;
        if( j_res == (molecule_src->nres) - 1 )
        {
          jlast = molecule_src->natoms;
        }
        else
        {
          jlast = molecule_src->residues[j_res + 1].firstatom;
        }

        for( j = jfirst; j < jlast; j++ )
        {

          d = 0.0;
          
			    for( k = 0; k < 3; k++ )
          {
            d += (obs_origins[ i ][ k ] - src_origins[ j ][ k ]) * (obs_origins[ i ][ k ] - src_origins[ j ][ k ]);
          }
          d = sqrt(d);
			    
			    if( d < ( molecule_src->atoms[j].r + molecule_obs->atoms[i].r + 2.0 * r_s ) )
          {
            if( d < 1.8 )
            {
              f_ij = 1.0 - (0.8875 * p_i * compute_b( molecule_obs->atoms[i].r, molecule_src->atoms[j].r, d, r_s) ) / S_i;              
            }
            else
            {
              f_ij = 1.0 - (0.3516 * p_i * compute_b( molecule_obs->atoms[i].r, molecule_src->atoms[j].r, d, r_s) ) / S_i;
            }
            A_i *= f_ij;
			    }
        } //end j loop
		  } //end j_res loop

		  sasa[ i ] = A_i;
    } //end i loop
	} //end i_res loop
}

//
// calculate_complex_cudesa_step_cpu() - Calculates desolvation energy
//                                       using the CUDESA CPU
//                                       algorithm. This is the
//                                       worker function that runs on
//                                       each CUDESA CPU thread
//
// input                               - arg: cpu_cudesa_step_arg_t
//                                       thread argument structure
//
// output                              - output desolvation energies
//                                       are stored in
//                                       arg->transformations
//
static void* calculate_complex_cudesa_step_cpu( void* arg )
{
  cpu_cudesa_step_arg_t* thread_arg = (cpu_cudesa_step_arg_t*)arg;
  PQRData_t* ligand = thread_arg->ligand;
  PQRData_t* receptor = thread_arg->receptor;
  
  
  protein_transformation_t* transform = thread_arg->transformations;

  double cudesa_energy;
  
  srand( time( NULL ) );                 

  double* ligand_inter_access = (double*)malloc( sizeof(double) * ligand->natoms );
  double* receptor_inter_access = (double*)malloc( sizeof(double) * receptor->natoms );

  double* ligand_intra_access = (double*)malloc( sizeof(double) * ligand->natoms );
  double* receptor_intra_access = (double*)malloc( sizeof(double) * receptor->natoms );

  double surface_area;
  double buried_area;
  
  vec3_t* ligand_transform_origins = (vec3_t*)malloc( sizeof(vec3_t) * ligand->natoms );
  vec3_t* receptor_transform_origins = (vec3_t*)malloc( sizeof(vec3_t) * receptor->natoms );

  for( int receptor_atom = 0; receptor_atom < receptor->natoms; receptor_atom++ )
  {
    for( int i = 0; i < 3; i++)
    {
      // identity transform for now for receptor
      receptor_transform_origins[ receptor_atom ][ i ] = receptor->atoms[ receptor_atom ].xyz[ i ];
    }
  }
  
  while( transform )
  {
    cudesa_energy = 0.0;
    
    surface_area = 0.0;
    buried_area = 0.0;

    // Translate the complex atoms
    for( int ligand_atom = 0; ligand_atom < ligand->natoms; ligand_atom++ )
    {    
      for( int i = 0; i < 3; i++)
      {
        ligand_transform_origins[ ligand_atom ][ i ] = transform->translation[ i ];
        for( int j = 0; j < 3; j++)
        {
          ligand_transform_origins[ ligand_atom ][ i ] += transform->rotation[ i ][ j ] * (ligand->atoms[ ligand_atom ].xyz[ j ]);
        }
      }
    }

    // Calculate the solvent accessible surface area
    calculate_intramolecular_sasa( receptor, receptor_transform_origins, receptor_intra_access );
    calculate_intramolecular_sasa( ligand, ligand_transform_origins, ligand_intra_access );
    
    calculate_intermolecular_sasa( ligand, receptor, ligand_transform_origins, receptor_transform_origins, ligand_inter_access );
    calculate_intermolecular_sasa( receptor, ligand, receptor_transform_origins, ligand_transform_origins, receptor_inter_access );    
                  
    for( int i = 0; i < receptor->natoms; i++ )
    {
      cudesa_energy += ((receptor_inter_access[ i ] * receptor_intra_access[ i ]) - receptor->atoms[ i ].access) * receptor->asp[ i ];
      surface_area += receptor_inter_access[ i ] * receptor_intra_access[ i ];
      buried_area += (receptor_inter_access[ i ] * receptor_intra_access[ i ]) - receptor->atoms[ i ].access;
    }

    for( int i = 0; i < ligand->natoms; i++ )
    {
      cudesa_energy += ((ligand_inter_access[ i ] * ligand_intra_access[ i ]) - ligand->atoms[ i ].access) * ligand->asp[ i ];
      surface_area += ligand_inter_access[ i ] * ligand_intra_access[ i ];
      buried_area += (ligand_inter_access[ i ] * ligand_intra_access[ i ]) - ligand->atoms[ i ].access;
    }
    
    transform->calculated_cpu_complex_energy = cudesa_energy;
    transform->calculated_sasa = surface_area;
    transform->calculated_buried_area = buried_area;
    transform->calculated_ligand_sasa = ligand->sasa;
    transform->calculated_receptor_sasa = receptor->sasa;
    for( int i = 0; i < thread_arg->stride; i++ )
    {
      if( transform )
      {
        transform = transform->next;
      }
    }
  }

  free( ligand_transform_origins );
  free( receptor_transform_origins );
  free( ligand_inter_access );
  free( ligand_intra_access );
  free( receptor_inter_access );
  free( receptor_intra_access );
  
  return NULL;
} 

//
// calculate_complex_cudesa_cpu() - Main entry point to run the CUDESA
//                                  CPU algorithm. Divides work
//                                  amongst various worker threads.
// input                          - receptor: receptor molecule
//                                  ligand: ligand molecule
//                                  transformations: list of
//                                  transformations for which to
//                                  compute desolvation energies
//
// output                         - helper functions will store output
//                                  energies in the transformations
//                                  array
//
void calculate_complex_cudesa_cpu( PQRData_t* receptor, PQRData_t* ligand, protein_transformation_t* transformations )
{
  pthread_t threads[ NUM_CPU_THREADS ];
  cpu_cudesa_step_arg_t thread_args[ NUM_CPU_THREADS ];

  assign_access( receptor );
  assign_access( ligand );

  for( int i = 0; i < NUM_CPU_THREADS; i++ )
  {
    thread_args[ i ].receptor = receptor;
    thread_args[ i ].ligand = ligand;
    thread_args[ i ].stride = NUM_CPU_THREADS;
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

    DEBUG_PRINT( CUDESA_PRINT_DEBUG, "%s(): Creating thread %d\n", __FUNCTION__, i );
    pthread_create( &threads[ i ], NULL, calculate_complex_cudesa_step_cpu, (void*)&thread_args[ i ] );                    
  }

  for( int i = 0; i < NUM_CPU_THREADS; i++ )
  {
    pthread_join( threads[ i ], NULL );
  }
}
