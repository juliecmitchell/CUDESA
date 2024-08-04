/*********************************************************************** 
 * vector.c - Vector utility functions
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

#include "vector.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

//
// Standard direction vectors
// 
vec3_t vec3_x = { 1.0f, 0.0f, 0.0f };
vec3_t vec3_y = { 1.0f, 0.0f, 0.0f };
vec3_t vec3_z = { 1.0f, 0.0f, 0.0f };

//
// vector_set() - sets all components of the vector to the specified
//                value
//
void vector_set( vec3_t v, double value )
{
  v[ 0 ] = value;
  v[ 1 ] = value;
  v[ 2 ] = value;
}

//
// vector_set_xyz() - sets the x, y, and z components of the vector to
//                    the specified values
//
void vector_set_xyz( vec3_t v, double x, double y, double z )
{
  v[ 0 ] = x;
  v[ 1 ] = y;
  v[ 2 ] = z;
}

//
// vector_dot() - computes the scalar product of the specified vectors
//
double vector_dot( vec3_t a, vec3_t b )
{
  return a[ 0 ] * b[ 0 ] + a[ 1 ] * b[ 1 ] + a[ 2 ] * b[ 2 ];
}

//
// vector_subtract() - performs a - b and stores it in res
//
void vector_subtract( vec3_t a, vec3_t b, vec3_t res )
{
  res[ 0 ] = a[ 0 ] - b[ 0 ];
  res[ 1 ] = a[ 1 ] - b[ 1 ];
  res[ 2 ] = a[ 2 ] - b[ 2 ];
}

//
// vector_add() - performs a + b and stores it in res
//
void vector_add( vec3_t a, vec3_t b, vec3_t res )
{
  res[ 0 ] = a[ 0 ] + b[ 0 ];
  res[ 1 ] = a[ 1 ] + b[ 1 ];
  res[ 2 ] = a[ 2 ] + b[ 2 ];
}

//
// vector_scale() - scales the vector by scale
//
void vector_scale( vec3_t v, double scale )
{
  v[ 0 ] *= scale;
  v[ 1 ] *= scale;
  v[ 2 ] *= scale;
}

//
// vector_copy() - copies vector a to b
//
void vector_copy( vec3_t a, vec3_t b )
{
  b[ 0 ] = a[ 0 ];
  b[ 1 ] = a[ 1 ];
  b[ 2 ] = a[ 2 ];
}

//
// vector_zero() - zeros the specified vector
//
void vector_zero( vec3_t v )
{
  v[ 0 ] = 0.0f;
  v[ 1 ] = 0.0f;
  v[ 2 ] = 0.0f;
}

//
// vector_normalize() - normalizes the specified vector
//
void vector_normalize( vec3_t v )
{
  double d = sqrtf( v[ 0 ] * v[ 0 ] + v[ 1 ] * v[ 1 ] + v[ 2 ] * v[ 2 ] );
  v[ 0 ] /= d;
  v[ 1 ] /= d;
  v[ 2 ] /= d;
}

//
// vector_length() - computes the length of the specified vector
//
double vector_length( vec3_t v )
{
  return sqrtf( v[ 0 ] * v[ 0 ] + v[ 1 ] * v[ 1 ] + v[ 2 ] * v[ 2 ] );
}

//
// Vector printing
//

const static int NUM_VEC_PRINT_BUFFERS = 16;
static int vec_print_buffer_index = 0;
static char vec_print_buffer[ NUM_VEC_PRINT_BUFFERS ][ 256 ];

//
// vector_print() - Prints the specified vector to a static buffer,
//                  which is returned
//
char* vector_print( vec3_t v )
{
  memset( vec_print_buffer[ vec_print_buffer_index % NUM_VEC_PRINT_BUFFERS ], 0, 256 );
  sprintf( vec_print_buffer[ vec_print_buffer_index % NUM_VEC_PRINT_BUFFERS ], "<%.5f, %.5f, %.5f>", v[ 0 ], v[ 1 ], v[ 2 ] );
  return vec_print_buffer[ vec_print_buffer_index++ % NUM_VEC_PRINT_BUFFERS ];
}

