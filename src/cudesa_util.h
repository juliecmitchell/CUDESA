/*********************************************************************** 
 * cudesa_util.h - Utility declarations
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

#ifndef _CUDESA_UTIL_H_
#define _CUDESA_UTIL_H_

#define DEBUG_ATOM 0

#include <getopt.h>
#include <limits.h>

// Print levels - used with CUDESA print routines to determine when a
//                message will be printed
//
enum cudesa_print_levels_t
{
  CUDESA_PRINT_ALWAYS = -1,
  CUDESA_PRINT_DEBUG,
  CUDESA_PRINT_DEBUG_EXTRA,
  CUDESA_PRINT_DEBUG_SPEW,
  CUDESA_PRINT_DEBUG_SUPER_SPEW
};

enum cudesa_run_mode_t
{
  CUDESA_RUN_CUDESA
};

// CUDESA settings
struct cudesa_settings_t
{
  int verbose_level;
  const char* binary_name;
  const char* binary_version;
  int transform_multiplier;
  int track_timings;
  cudesa_run_mode_t run_mode;
  int use_cpu;
  int use_gpu;
};

extern struct cudesa_settings_t cudesa_settings;
extern char* short_options;
extern struct option long_options[];
extern char* optarg;

//
// External utility functions & macros
//
void cudesa_print_usage( void );
void cudesa_print_version( void );
void cudesa_error( char* fmt, ... );
extern "C"{
void cudesa_print( int message_verbose_level, char* fmt, ... );
}
#ifdef _CUDESA_DEBUG
#define DEBUG_PRINT( level, format, ... ) cudesa_print( level, format, __VA_ARGS__ )
#define DEBUG_PRINT_NO_ARGS( level, format ) cudesa_print( level, format )
#else
#define DEBUG_PRINT( level, format, ... )
#define DEBUG_PRINT_NO_ARGS( level, format )
#endif

#endif
