/*********************************************************************** 
 * cudesa_util.c - Utility code
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

#include "cudesa_util.h"
#include <stdarg.h>
#include <stdio.h>

//
// Authors, copyright notice
//
#define CUDESA_AUTHORS "David Dynerman, Erick Butzlaff, Julie C. Mitchell"
#define CUDESA_RELEASE_YEARS "2007"

//
// Global settings
//
// This struct contains global settings for the execution of the program
//
struct cudesa_settings_t cudesa_settings =
{
  0,
  BINARY_NAME,
  BINARY_VERSION,
  false
};

//
// getopt settings
//
// Settings for getopt command line parsing
//
// Also see cudesa_print_usage() below
//
struct option long_options[] =
{
  {"gpu", no_argument, 0, 'g' },
  {"both",no_argument, 0, 'b' },
  {"receptor", required_argument, 0, 'r' },
  {"ligand", required_argument, 0, 'l' },
  {"complex", required_argument, 0, 'c' },
  {"help", no_argument, 0, 'h' },
  {"verbose", optional_argument, 0, 'v' },
  {"timings", no_argument, 0, 't' },
  {"transformations", required_argument, 0, 'm' },
  {"multiplier", required_argument, 0, 'x' },
  {"p_i", required_argument, 0, 'p' },
  { 0, 0, 0, 0 }
};

char* short_options = "bgvm:thr:l:c:b:x:p:";

//
// Utility functions
//

//
// cudesa_print()
//
// Prints the specified message if the program's verbose level matches
//
void cudesa_print( int message_verbose_level, char* fmt, ... )
{
  if( message_verbose_level >= cudesa_settings.verbose_level )
  {
    return;
  }

	va_list ap;

	va_start( ap, fmt );
	vfprintf( stdout, fmt, ap );
	va_end( ap ); 
}

//
// cudesa_error()
//
// Indicates an error has occured
//
void cudesa_error( char* fmt, ... )
{
	va_list ap;

	va_start( ap, fmt );
	fprintf( stderr, "[ERROR] " );
	vfprintf( stderr, fmt, ap );
	va_end( ap );
}

//
// cudesa_print_version()
//
// Displays version information
//
void cudesa_print_version( void )
{
  cudesa_print( CUDESA_PRINT_ALWAYS, "\n");
  cudesa_print( CUDESA_PRINT_ALWAYS, "%s v%s\n\n", cudesa_settings.binary_name, cudesa_settings.binary_version );
  cudesa_print( CUDESA_PRINT_ALWAYS, "Copyright (C) %s, Regents of the University of Wisconsin\n", CUDESA_RELEASE_YEARS );
  cudesa_print( CUDESA_PRINT_ALWAYS, "Authors: %s\n\n", CUDESA_AUTHORS );
  cudesa_print( CUDESA_PRINT_ALWAYS, "This software is free for non-commercial purposes only.  Please see the \n" );
  cudesa_print( CUDESA_PRINT_ALWAYS, "licensing information before copying, redistributing or modifying CUDESA. \n\n" );
}

//
// cudesa_print_usage()
//
// Displays usage information
//
void cudesa_print_usage( void )
{
  static char* usage =
    "Usage:\t%s\t[-r|--receptor] [receptor filename]\n" \
    "\t\t[-l|--ligand] [ligand filename]\n" \
    "\t\t[-m|--transformations] [transformations filename]\n" \
    "\t\t[-g|--gpu] use GPU calculation\n" \
    "\t\t[-b|--both] use both CPU and GPU calculation\n" \
    "\t\t[-t|--timings]\n" \
    "\t\t[-v|--verbose] {verbose level}}\n"    \
    "\t\t[-x|--multiplier] [multiplier]\n"     \
    "\t\t[-p|--p_i] [p_i index] [p_i value]\n" \
    "\n" \
    "\t%s\t[-h|--help]\n" \
/*    "\n" \
    "\t--receptor [receptor filename]\n" \
    "\tSpecifies the PQR file that will act as the receptor protein\n"
    "\n" \
    "\t--ligand [ligand filename]\n" \
    "\tSpecifies the PQR file that will act as the ligand protein\n"
    "\n" \
    "\t--transformations [transformations filename]\n" \
    "\tdome-style file specifying iterative transformations of ligand\n" \
    "\n" \
    "\t--timings\n"
    "\tPrint wall-clock timing information\n" \
    "\n" \
    "\t--verbose {verbose level}\n" \
    "\tPrint more output, optionally specifying how much more (default: 1)\n" \
    "\n" \
    "\t--help\n" \
    "\tPrint this screen and exit\n" \ */
    "\n" ;
    
  cudesa_print( CUDESA_PRINT_ALWAYS, usage, cudesa_settings.binary_name, cudesa_settings.binary_name, cudesa_settings.binary_name, cudesa_settings.binary_name );
}

