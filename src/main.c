/*********************************************************************** 
 * main.c - Main program routines
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
#include <sys/time.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "cudesa_dome.h"
#include "cudesa_transform.h"
#include "cudesa_util.h"
#include "cudesa.h"

// DEFAULT_NUM_TRANSFORMS - default number of random transforms to generate
//                          if a transformation file isn't specified - the
//                          multiplier parameter is applied to this
#ifdef _DEBUG
static const int DEFAULT_NUM_TRANSFORMS = 3;
#else
static const int DEFAULT_NUM_TRANSFORMS = 2000;
#endif

//
// main_cudesa() - The main program routine for both CPU and GPU
//                 versions of the CUDESA algorithm
//
int main_cudesa( PQRData_t* ligand, PQRData_t* receptor, protein_transformation_t* transformations, char* transformations_filename )
{
  struct timeval time1, time2;
  unsigned long long cpu_time, gpu_time, cpu_complex_time, gpu_complex_time;
  int num_transforms = 0;
  protein_transformation_t* transform = NULL;

  memset( &time1, 0, sizeof( struct timeval ) );
  memset( &time2, 0, sizeof( struct timeval ) );

  cpu_time = gpu_time = cpu_complex_time = gpu_complex_time = 0;

  zero_transformation_calculated_energies( transformations );

  if( cudesa_settings.use_cpu )
  {
    cudesa_print( CUDESA_PRINT_ALWAYS, "Beginning CPU cudesa calculations\n" );
    if( cudesa_settings.track_timings )
    {
      gettimeofday( &time1, NULL );
    }
    
    if( cudesa_settings.track_timings )
    {
      gettimeofday( &time2, NULL );
      cpu_time = ((time2.tv_sec - time1.tv_sec) * 1000000) + (time2.tv_usec - time1.tv_usec);
    }
    
    if( cudesa_settings.track_timings )
    {
      gettimeofday( &time1, NULL );
    }

    //
    // Compute desolvation energies using CUDESA CPU
    //
    calculate_complex_cudesa_cpu( receptor, ligand, transformations );  
    
    if( cudesa_settings.track_timings )
    {
      gettimeofday( &time2, NULL );
      cpu_complex_time = ((time2.tv_sec - time1.tv_sec) * 1000000) + (time2.tv_usec - time1.tv_usec);
    } 
  } // if( cudesa_settings.use_gpu )
  
  if( cudesa_settings.use_gpu )
  {
    cudesa_print( CUDESA_PRINT_ALWAYS, "Beginning GPU cudesa calculations\n" );

    if( cudesa_settings.track_timings )
    {
      gettimeofday( &time1, NULL );
    }

    if( cudesa_settings.track_timings )
    {
      gettimeofday( &time2, NULL );
      gpu_time = ((time2.tv_sec - time1.tv_sec) * 1000000) + (time2.tv_usec - time1.tv_usec);    
    }
    
    if( cudesa_settings.track_timings )
    {
      gettimeofday( &time1, NULL );
    }

    //
    // Compute desolvation energies using CUDESA GPU
    //
    calculate_complex_cudesa_gpu( receptor, ligand, transformations );

    if( cudesa_settings.track_timings )
    {
      gettimeofday( &time2, NULL );
      gpu_complex_time = ((time2.tv_sec - time1.tv_sec) * 1000000) + (time2.tv_usec - time1.tv_usec);    
    }
  } // if( cudesa_settings.use_gpu )
  
  transform = transformations;
 
  cudesa_print( CUDESA_PRINT_ALWAYS, "Desolvation results:\n" );
  while( transform )
  {
    cudesa_print( CUDESA_PRINT_ALWAYS, "Transform %d:\n", num_transforms );
    
    if (cudesa_settings.use_cpu) 
    {
      cudesa_print( CUDESA_PRINT_ALWAYS, "\t\tCPU Complex Desolvation Energy (cal): %10.3f\n\t\tBound SASA (sq. ang): %10.3f [Unbound SASA (sq. ang) ligand: %10.3f, receptor: %10.3f]\n\t\tBA (sq. ang): %10.3f\n", 
                    transform->calculated_cpu_complex_energy,  
                    transform->calculated_sasa, transform->calculated_ligand_sasa, transform->calculated_receptor_sasa, transform->calculated_buried_area );
    }
    
    if (cudesa_settings.use_gpu) 
    {
      cudesa_print( CUDESA_PRINT_ALWAYS, "\t\tGPU Complex Desolvation Energy (cal): %10.3f\n\t\tSASA (sq. ang): %10.3f [ligand: %10.3f, receptor: %10.3f]\n\t\tBA (sq. ang): %10.3f\n", 
                    transform->calculated_gpu_complex_energy,
                    transform->calculated_gpu_sasa, transform->calculated_gpu_ligand_sasa, transform->calculated_gpu_receptor_sasa, transform->calculated_gpu_buried_area );
    }

    if( cudesa_settings.use_gpu && cudesa_settings.use_cpu )
    {
      cudesa_print( CUDESA_PRINT_ALWAYS, "\t\tCPU/GPU Comparison: CPU/GPU Energy Delta (cal): %10.3f\n", 
                    fabs( transform->calculated_gpu_complex_energy - transform->calculated_cpu_complex_energy ) );

      if( fabs( transform->calculated_gpu_complex_energy - transform->calculated_cpu_complex_energy ) > 0.1)
      {
        cudesa_error( "\t\tDifference between CPU and GPU Solvation energy exceeded 0.1: CPU: %f, GPU: %f Delta: %f\n",
                      num_transforms, transform->calculated_cpu_complex_energy, transform->calculated_gpu_complex_energy,
                      fabs( transform->calculated_gpu_complex_energy - transform->calculated_cpu_complex_energy ) );
      }
      
    }
    cudesa_print( CUDESA_PRINT_ALWAYS, "\n" );
    num_transforms++;
    transform = transform->next;
  }


  if( cudesa_settings.track_timings )
  {
    cudesa_print( CUDESA_PRINT_ALWAYS, "=======================\n" );
    cudesa_print( CUDESA_PRINT_ALWAYS, "Timing Analysis:\n" );
    cudesa_print( CUDESA_PRINT_ALWAYS, "Ligand size: %d atoms, Receptor size: %d, Total: %d\n", ligand->natoms, receptor->natoms, ligand->natoms + receptor->natoms );
    cudesa_print( CUDESA_PRINT_ALWAYS, "Number of transformations: %d (source: %s)\n", num_transforms, (transformations_filename == NULL) ? transformations_filename : "randomly generated" );
    
    if( cudesa_settings.use_cpu )
    {
      cudesa_print( CUDESA_PRINT_ALWAYS, "CPU Timings: %llu usec (%.2f sec, %.1f mins)\n", cpu_complex_time, (double)cpu_complex_time / 1000000.0, (double)cpu_complex_time / (1000000.0 * 60.0) );
    }

    if( cudesa_settings.use_gpu )
    {
      cudesa_print( CUDESA_PRINT_ALWAYS, "GPU Timings: %llu usec (%.2f sec, %.1f mins)\n", gpu_complex_time, (double)gpu_complex_time / 1000000.0, (double)gpu_complex_time / (1000000.0 * 60.0) );      
    }

    if( cudesa_settings.use_gpu && cudesa_settings.use_cpu )
    {
      cudesa_print( CUDESA_PRINT_ALWAYS, "CPU Speed/GPU Speed = %llu / %llu = %.2f\n", cpu_complex_time, gpu_complex_time, (double)cpu_complex_time / (double)gpu_complex_time );
    }
    else
    {
      cudesa_print( CUDESA_PRINT_ALWAYS, "To automatically compare CPU/GPU timings, use both the '-b' and '-t' flags\n" );
    }
    cudesa_print( CUDESA_PRINT_ALWAYS, "=======================\n" );    
  }
  
  return 0;
}

//
// main() - Program entry point
//
int main( int argc, char* argv[] )
{
  PQRData_t* receptor = NULL;
  PQRData_t* ligand = NULL;
  protein_transformation_t* transformations;
  
  int ret_code = 0;
    
  int option;
  int option_index;

  char ligand_filename[ PATH_MAX ] = { '\0' };
  char receptor_filename[ PATH_MAX ] = { '\0' };
  char transformations_filename[ PATH_MAX ] = { '\0' };

  cudesa_settings.run_mode = CUDESA_RUN_CUDESA;
  cudesa_settings.transform_multiplier = 1;
  cudesa_settings.use_cpu = 1;
  cudesa_settings.use_gpu = 0;
  
  while((option = getopt_long( argc, argv, short_options, long_options, &option_index)) != -1)
  {
    //cudesa_print( CUDESA_PRINT_ALWAYS, "option: %c\n", option );
    switch( option )
    {
      // --gpu
      case 'b':
      {
        cudesa_settings.use_cpu = 1;
        cudesa_settings.use_gpu = 1;
        break;
      }
      // --gpu
      case 'g':
      {
        cudesa_settings.use_cpu = 0;
        cudesa_settings.use_gpu = 1;
        break;
      }
      // --receptor
      case 'r':
      {
        strncpy( receptor_filename, optarg, PATH_MAX - 1 );
        receptor_filename[ PATH_MAX - 1 ] = '\0';
        break;
      }
      // --ligand
      case 'l':
      {
        //cudesa_print( CUDESA_PRINT_ALWAYS, "l optarg: %s\n", optarg );        
        strncpy( ligand_filename, optarg, PATH_MAX - 1 );
        ligand_filename[ PATH_MAX - 1 ] = '\0';
        break;
      }
      // --verbose
      case 'v':
      {
        cudesa_settings.verbose_level = 1;
        if( NULL != optarg )
        {
          cudesa_settings.verbose_level = atoi( optarg );
        }
        break;
      }
      // --help
      case 'h':
      {
        cudesa_print_version();
        cudesa_print_usage();
        exit( 0 );
        break;
      }
      // --timings
      case 't':
      {
        cudesa_settings.use_cpu = 1;
        cudesa_settings.use_gpu = 1;
        cudesa_settings.track_timings = 1;
        break;
      }
      // --transformations
      case 'm':
      {
        strncpy( transformations_filename, optarg, PATH_MAX - 1 );
        transformations_filename[ PATH_MAX - 1 ] = '\0';
        break;
      }
      // --multiplier
      case 'x':
      {
        cudesa_settings.transform_multiplier = atoi( optarg );
        break;
      }
      // --p_i
      case 'p':
      {
        if( optind > argc - 1 )
        {
          cudesa_error( "%s(): Specified p_i index without specifying what p_i value to set\n", __FUNCTION__ );
        }
        else
        {
          p_i_values[ atoi(optarg) - 1 ] = atof( argv[ optind ] );
          DEBUG_PRINT( CUDESA_PRINT_DEBUG, "%s(): Setting %d p_i to %f\n", __FUNCTION__, atoi(optarg), atof( argv[ optind ] ) );
          optind++;
        }
        break;
      }
    }
  }

  cudesa_print_version();

  /*
  cudesa_print( CUDESA_PRINT_ALWAYS, "Using the following p_i values:\nAtom type\t\tp_i\n" );
  for( int i = 0; i < NUM_P_I_VALUES; i++ )
  {
    cudesa_print( CUDESA_PRINT_ALWAYS, "%d\t\t\t%f\n", i + 1, p_i_values[ i ] );
  }
  */
  
  switch( cudesa_settings.run_mode )
  {
    case CUDESA_RUN_CUDESA:
    {
      if( '\0' == *ligand_filename || '\0' == *receptor_filename )
      {
        cudesa_print_usage();
        exit( 1 );
      }

      receptor = ReadPQR( receptor_filename, true );
      
      if( NULL == receptor )
      {
        cudesa_error( "%s(): Error opening receptor file %s\n", __FUNCTION__, receptor_filename );
        exit( 1 );
      }
      
      ligand = ReadPQR( ligand_filename, true );
      
      if( NULL == ligand )
      {
        FreePQR( receptor );
        cudesa_error( "%s(): Error opening ligand file %s\n", __FUNCTION__, ligand_filename );
        exit( 1 );
      }

      if( '\0' != *transformations_filename )
      {
        transformations = load_transformations( transformations_filename, cudesa_settings.transform_multiplier );
        
        if( NULL == transformations )
        {
          cudesa_error( "%s(): Error reading transformations from file %s\n", __FUNCTION__, transformations_filename );
          FreePQR( ligand );
          FreePQR( receptor );
          exit( 1 );
        }
      }
      else
      {
        // If a transformations filename isn't specified

        // If we're timing, generate random transformations
        if( cudesa_settings.track_timings )
        {
          transformations = load_random_transformations( DEFAULT_NUM_TRANSFORMS, cudesa_settings.transform_multiplier );
          
          if( NULL == transformations )
          {
            cudesa_error( "%s(): Error generating random transformations\n", __FUNCTION__ );
            FreePQR( ligand );
            FreePQR( receptor );
            exit( 1 );
          }
        }
        else
        {
          // If we aren't timing, just compute in-place
          transformations = load_identity_transformation();
          
          if( NULL == transformations )
          {
            cudesa_error( "%s(): Error loading identity transformation\n", __FUNCTION__ );
            FreePQR( ligand );
            FreePQR( receptor );
            exit( 1 );
          }
        }
      }

      ret_code = main_cudesa( ligand, receptor, transformations, transformations_filename );

      FreePQR( ligand );
      FreePQR( receptor );
      
      break;
    }
    default:
    {
      cudesa_error( "%s(): Unknown cudesa run_mode: %d\n", __FUNCTION__, cudesa_settings.run_mode );
      FreePQR( receptor );
      exit( 1 );
    }
  }

  return ret_code;
}
