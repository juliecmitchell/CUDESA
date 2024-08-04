/*********************************************************************** 
 * cudesa_transform.h - Declarations for ligand/receptor transformations
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

//
// Transformation structures and methods
//

// protein_transformation_t
// a linked list of ligand transformations relative to the receptor
//
typedef struct protein_transformation_s
{
  double rotation[3][3];
  double rotation_angles[3];
  double translation[3];

  // single precision version for the GPU
  float single_rotation[3][3];
  float single_translation[3];
  
  double reference_energy; // reference energy of cudesa
  double calculated_cpu_energy;
  double calculated_gpu_energy;
  double calculated_cpu_complex_energy;
  double calculated_gpu_complex_energy;
  double calculated_sasa;
  double calculated_buried_area;
  double calculated_ligand_sasa;
  double calculated_receptor_sasa;
  double calculated_gpu_sasa;
  double calculated_gpu_buried_area;
  double calculated_gpu_ligand_sasa;
  double calculated_gpu_receptor_sasa;
  struct protein_transformation_s* next;
} protein_transformation_t;

//
// free_transformations
//
// de-allocate the linked lists's memory
//
void free_transformations( protein_transformation_t* transformations );

//
// zero_transformation_calculated_energies
//
// zero the calculated energies of all transformations in the linked
// list
//
void zero_transformation_calculated_energies( protein_transformation_t* transformations );

//
// load_transformations
//
// load the transformations from the specified file, optionally
// performing each one more than once
//
protein_transformation_t* load_transformations( char* transformations_filename, int transform_multiplier );

//
// load_random_transformations
//
// loads the specified number of random transformations
//
protein_transformation_t* load_random_transformations( int num_transforms, int transformation_multiplier );

//
// load_identity_transformation
//
// returns an identity transformation
//
protein_transformation_t* load_identity_transformation( void );
