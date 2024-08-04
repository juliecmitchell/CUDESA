/*********************************************************************** 
 * cudesa.h - Program declarations
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

#include "cudesa_dome.h"

//
// Main program entrypoints
//

//
// calculate_complex_cudesa_cpu() - Run the CPU CUDESA algorithm
//
void calculate_complex_cudesa_cpu( PQRData_t* receptor, PQRData_t* ligand, protein_transformation_t* transformations );

//
// calculate_complex_cudesa_gpu() - Run the GPU CUDESA algorithm
//
void calculate_complex_cudesa_gpu( PQRData_t* receptor, PQRData_t* ligand, protein_transformation_t* transformations );

