/*********************************************************************** 
 * cudesa_dome.h
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
 
#ifndef _CUDESA_DOME_H_
#define _CUDESA_DOME_H_

#define VERBOSE 0
#define PQR_FILENAME_MAX 256

//
// atom_t - structure representing an individual atom
//
typedef struct
{
	char	aname[5];	/* atom name */
	double	xyz[3];		/* original xyz coordinates */
	double	XYZ[3];		/* transformed xyz coordinates */
	double	q;		/* atom charge */
	double	r;		/* atom radius */
	int     s;    /* ACE solvation type */
	double	access;		/* accessibility */
	double  accessgrad[3];      
	double  abc[3];
} atom_t;

//
// residue_t - structure representing a residue
//
typedef struct
{
	char		rname[5];	/* residue name */
	int			rnum;		/* residue number */
	int			firstatom;	/* index into atom array */
	int			ca;			/* position of CA atom */
	int			surface;	/* surface residue */
	double 		phi; 		/* phi angle */
	double		psi; 		/* psi angle */
} residue_t;

//
// PQRData_t - main structure representing a loaded molecular
//             structure
//
typedef struct
{
	char            filename[ PQR_FILENAME_MAX ];
	int		          natoms;		/* number of atoms in the molecule */
	int		          nres;		/* number of residues in the molecule */
	atom_t*		      atoms;		/* array of atom_t's */
	residue_t*	    residues;	/* array of residue_t's */
	double		      tq;		/* total charge on the molecule */
	double          molradius;      /* radius of the molecule */   
	double          smolradius;     /* radius of the maximal sphere inside molecule */   
	double          Emat[3][3];     /* ellipsoid matrix */
	double          sasa; /* total sasa */
	float*          atomicRadii; /* array of atomic radii */
	float*          p_iParams; /* array of p_i parameters per Section
                                (1.1.2) */
	int*           is_hydrogen; 
	int*           polarity;
	float*          S_i; /* s_i parameters per Section (1.1.2) */
	double*         asp; /* ASP parameters per Section (1.2) */
} PQRData_t;

//
// Function declarations for cudesa
//

PQRData_t* ReadPQR( char *pqrfile, int prep_for_gpu );
void FreePQR( PQRData_t* molecule );
double determine_asp( char *aname, char *rname );
void assign_access( PQRData_t *molecule );
inline double compute_b( double r_i, double r_j, double d_ij, double r_s );
// Other cudesa requirements
#define FMIN( a, b ) ( (a) < (b) ? (a) : (b) )

// This can be converted to a variable, but is unchanged right now, so we just define it
#define SOLVATION_RADIUS 1.4

//#define DEBUG_ATOM 0

#define ITERATIONS 10000

#define NUM_P_I_VALUES 23
extern double p_i_values[ NUM_P_I_VALUES ];

//#define _DEBUG 2
#endif
