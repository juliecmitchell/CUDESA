/*********************************************************************** 
 * cudesa_transform.c - Ligand/receptor transformation code
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "cudesa_util.h"
#include "cudesa_transform.h"

// TRANSLATION_RANDOM_RANGE - random transformations will translate within +/-
//                            this value / 2
static const double TRANSLATION_RANDOM_RANGE = 10.0;

// ROTATION_RANDOM_RANGE - random rotations will rotat within [0, this value]
static const double ROTATION_RANDOM_RANGE = 2.0 * M_PI;

//
// create_rotation_matrix() - this function evaluates the rotation
//                            matrix r
//
static void create_rotation_matrix( double phi, double theta, double psi, double r[3][3] )
{
  double sin_phi   = sin( phi );
  double cos_phi   = cos( phi );
  double sin_theta = sin( theta );
  double cos_theta = cos( theta );
  double sin_psi   = sin( psi );
  double cos_psi   = cos( psi );
  
  r[0][0] =  cos_psi * cos_phi - sin_psi * cos_theta * sin_phi;
  r[1][0] =  sin_psi * cos_phi + cos_psi * cos_theta * sin_phi;
  r[2][0] =  sin_theta * sin_phi;
  
  r[0][1] = -cos_psi * sin_phi - sin_psi * cos_theta * cos_phi;
  r[1][1] = -sin_psi * sin_phi + cos_psi * cos_theta * cos_phi;
  r[2][1] =  sin_theta * cos_phi;
  
  r[0][2] =  sin_psi * sin_theta;
  r[1][2] = -cos_psi * sin_theta;
  r[2][2] =  cos_theta;

  return;
}

//
// zero_transformations() - initializes the given transformation
//
static void zero_transformation( protein_transformation_t* transformation )
{
  int i, j;

  for( i = 0; i < 3; i++ )
  {
    for( j = 0; j < 3; j++ )
    {
      transformation->rotation[ i ][ j ] = 0.0;
    }
    transformation->translation[ i ] = 0.0;
    transformation->rotation_angles[ i ] = 0.0;
  }

  transformation->reference_energy = 0.0;
  
  transformation->next = NULL;
}

//
// zero_transformation_calculated_energies() - re-initializes energy
//                                             information for a given
//                                             transformation
//
void zero_transformation_calculated_energies( protein_transformation_t* transformations )
{
  protein_transformation_t* transform = transformations;

  while( transform )
  {
    transform->calculated_cpu_energy = 0.0;
    transform->calculated_gpu_energy = 0.0;
    transform = transform->next;
  }
}

//
// free_transformations() - frees memory associated with the specified
//                          list of transformations
//
void free_transformations( protein_transformation_t* transformations )
{
  protein_transformation_t* next = NULL;

  while( transformations )
  {
    next = transformations->next;
    free( transformations );
    transformations = next;
  }
}

//
// load_transformations() - This function reads in a series of ligand
//                          transformations from a standard .dome
//                          text file.
//
// Line format:
// %f %f %f %f %f x y z phi theta psi %d
//
protein_transformation_t* load_transformations( char* transformations_filename, int transformation_multiplier )
{
  FILE* transformation_file = fopen( transformations_filename, "r" );
  double in_energy, in_x, in_y, in_z, in_phi, in_theta, in_psi, value;
  protein_transformation_t* transformations = NULL;
  protein_transformation_t* newTransformation = NULL;

  protein_transformation_t* newCopyTransformation = NULL;
  
  char line_buffer[ PATH_MAX ] = { '\0' };
  char* token = NULL;
  unsigned int parse_index = 0;
  int i, j;

  
  if( NULL == transformation_file )
  {
    cudesa_error( "%s(): I/O error opening file %s\n", __FUNCTION__, transformations_filename );
    return NULL;
  }

  while( NULL != fgets( line_buffer, PATH_MAX, transformation_file ) )
  {    
    token = strtok( line_buffer, " " );
    while( token )
    {
      value = atof( token );
      switch( parse_index )
      {
        case 0:
        case 1:
        case 2:
        case 3:
        case 11:
        {
          // ignored values
          break;
        }
        case 4:
        {
          in_energy = value;
          break;
        }
        case 5:
        {
          in_x = value;
          break;
        }
        case 6:
        {
          in_y = value;
          break;
        }
        case 7:
        {
          in_z = value;
          break;
        }
        case 8:
        {
          in_phi = value;
          break;
        }
        case 9:
        {
          in_theta = value;
          break;
        }
        case 10:
        {
          in_psi = value;
          break;
        }
        
      }
      parse_index++;
      token = strtok( NULL, " " );
    }

    if( parse_index != 12 )
    {
      cudesa_error( "%s(): Parsing error while reading file %s\n", __FUNCTION__, transformations_filename );
      free_transformations( transformations );
      fclose( transformation_file );
      return NULL;
    }
    
    DEBUG_PRINT( CUDESA_PRINT_DEBUG_SPEW, "%s(): Got xyz: <%.3f, %.3f, %.3f> rot: <%.3f, %.3f, %.3f> energy: %f from input\n", __FUNCTION__,
                 in_x, in_y, in_z, in_phi, in_theta, in_psi, in_energy );

    newTransformation = (protein_transformation_t*)malloc( sizeof( protein_transformation_t ) );
    if( NULL == newTransformation )
    {
      cudesa_error( "%s(): Error allocating memory for transformation object\n", __FUNCTION__ );
      free_transformations( transformations );
      fclose( transformation_file );
      return NULL;
    }
    
    zero_transformation( newTransformation );

    create_rotation_matrix( in_phi, in_theta, in_psi, newTransformation->rotation );
    
    newTransformation->rotation_angles[ 0 ] = in_phi;
    newTransformation->rotation_angles[ 1 ] = in_theta;
    newTransformation->rotation_angles[ 2 ] = in_psi;
    newTransformation->translation[ 0 ] = in_x;
    newTransformation->translation[ 1 ] = in_y;
    newTransformation->translation[ 2 ] = in_z;

    for( i = 0; i < 3; i++ )
    {
      newTransformation->single_translation[ i ] = (float)newTransformation->translation[ i ];
      for( j = 0; j < 3; j++ )
      {
        newTransformation->single_rotation[ i ][ j ] = (float)newTransformation->rotation[ i ][ j ];
      }
    }
    
    newTransformation->reference_energy = in_energy;

    for( i = 0; i < transformation_multiplier - 1; i++ )
    {
      newCopyTransformation = (protein_transformation_t*)malloc( sizeof( protein_transformation_t ) );
      memcpy( newCopyTransformation, newTransformation, sizeof(protein_transformation_t) );

      newCopyTransformation->next = transformations;
      transformations = newCopyTransformation;
    }
    
    newTransformation->next = transformations;
    transformations = newTransformation;

    // reset parsing state
    memset( line_buffer, 0, sizeof(char) * PATH_MAX );
    parse_index = 0;
  }

  fclose( transformation_file );
  
  return transformations;
}

//
// load_random_transformations() - Generates random transformations
//                                 for timing purposes
//
protein_transformation_t* load_random_transformations( int num_transforms, int transformation_multiplier )
{
  protein_transformation_t* transformations = NULL;
  protein_transformation_t* newTransformation = NULL;
  double in_x, in_y, in_z, in_phi, in_theta, in_psi;
  srand( time( NULL ) );
  
  for( int i = 0; i < num_transforms * transformation_multiplier; i++ )
  {
    newTransformation = (protein_transformation_t*)malloc( sizeof( protein_transformation_t ) );
    if( NULL == newTransformation )
    {
      cudesa_error( "%s(): Error allocating memory for transformation object\n", __FUNCTION__ );
      free_transformations( transformations );
      return NULL;
    }

    zero_transformation( newTransformation );

    in_x = (TRANSLATION_RANDOM_RANGE * ((double)rand() / ((double)RAND_MAX + 1.0))) - (TRANSLATION_RANDOM_RANGE / 2.0);
    in_y = (TRANSLATION_RANDOM_RANGE * ((double)rand() / ((double)RAND_MAX + 1.0))) - (TRANSLATION_RANDOM_RANGE / 2.0);
    in_z = (TRANSLATION_RANDOM_RANGE * ((double)rand() / ((double)RAND_MAX + 1.0))) - (TRANSLATION_RANDOM_RANGE / 2.0);    

    in_phi = (ROTATION_RANDOM_RANGE * ((double)rand() / ((double)RAND_MAX + 1.0)));
    in_theta = (ROTATION_RANDOM_RANGE * ((double)rand() / ((double)RAND_MAX + 1.0)));
    in_psi = (ROTATION_RANDOM_RANGE * ((double)rand() / ((double)RAND_MAX + 1.0)));    
    
    create_rotation_matrix( in_phi, in_theta, in_psi, newTransformation->rotation );
    
    newTransformation->rotation_angles[ 0 ] = in_phi;
    newTransformation->rotation_angles[ 1 ] = in_theta;
    newTransformation->rotation_angles[ 2 ] = in_psi;
    newTransformation->translation[ 0 ] = in_x;
    newTransformation->translation[ 1 ] = in_y;
    newTransformation->translation[ 2 ] = in_z;

    for( int j = 0; j < 3; j++ )
    {
      newTransformation->single_translation[ j ] = (float)newTransformation->translation[ j ];
      for( int k = 0; k < 3; k++ )
      {
        newTransformation->single_rotation[ j ][ k ] = (float)newTransformation->rotation[ j ][ k ];
      }
    }
    
    newTransformation->reference_energy = 0;

    newTransformation->next = transformations;
    transformations = newTransformation;
  }

  return transformations;
}

//
// load_identity_transformation() - Returns the identity transformation
//
protein_transformation_t* load_identity_transformation( void )
{
  protein_transformation_t* newTransformation = NULL;

  newTransformation = (protein_transformation_t*)malloc( sizeof( protein_transformation_t ) );
  if( NULL == newTransformation )
  {
    cudesa_error( "%s(): Error allocating memory for transformation object\n", __FUNCTION__ );
    return NULL;
  }

  zero_transformation( newTransformation );

  for( int j = 0; j < 3; j++ )
  {
    newTransformation->translation[ j ] = 0.0;
    newTransformation->rotation_angles[ j ] = 0.0;
    newTransformation->single_translation[ j ] = (float)newTransformation->translation[ j ];
    for( int k = 0; k < 3; k++ )
    {
      if( j == k )
      {
        newTransformation->rotation[ j ][ k ] = 1.0;
      }
      else
      {
        newTransformation->rotation[ j ][ k ] = 0.0;
      }
      
      newTransformation->single_rotation[ j ][ k ] = (float)newTransformation->rotation[ j ][ k ];
    }
  }
    
  newTransformation->reference_energy = 0;

  newTransformation->next = NULL;

  return newTransformation; 
}
