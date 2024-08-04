/*********************************************************************** 
 * vector.h - Vector utility functions
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

#ifndef _VECTOR_H_
#define _VECTOR_H_

typedef double vec3_t[3];

extern vec3_t vec3_x;
extern vec3_t vec3_y;
extern vec3_t vec3_z;

//
// Vector utility functions
//

void vector_set( vec3_t v, double value );
void vector_set_xyz( vec3_t v, double x, double y, double z );
double vector_dot( vec3_t a, vec3_t b );
void vector_subtract( vec3_t a, vec3_t b, vec3_t res );
void vector_add( vec3_t a, vec3_t b, vec3_t res );
char* vector_print( vec3_t v );
void vector_normalize( vec3_t v );
void vector_scale( vec3_t v, double scale );
void vector_copy( vec3_t a, vec3_t b );
void vector_zero( vec3_t a );
double vector_length( vec3_t v );
#endif
