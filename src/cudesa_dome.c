/*****************************************************************************
 * cudesa_dome.c
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
#include <string.h>
#include <math.h>

#include "cudesa_dome.h"
#include "cudesa_util.h"

//
// isDigit() - Determine if a character is a digit
//
inline int isDigit( char A )
{
  if ( (A == '0') || (A == '1') || (A == '2') || (A == '3') || (A == '4') ||
       (A == '5') || (A == '6') || (A == '7') || (A == '8') || (A == '9') )
  {
    return 1;
  }

  return 0;
}

//
// isH()  - Determine if a given atom is a hydrogen atom
//
// input  - sAtom: given atom name
//
// output - true/false
//
inline int isH( char *sAtom )
{
  if ( (sAtom[0] == 'H') ||
       ( (sAtom[0] == ' ') && (sAtom[1] == 'H') ) ||
       ( (isDigit(sAtom[0]) == 1) && (sAtom[1] == 'H') ) )
  {
    return 1;
  }

  return 0;
}

//
// determine_type() - Assigns p_i parameter types to atoms. See Table (1)
//                    for details.
//
// input            - aname: atom name
//                    rname: residue name
//
// output           - numerical p_i type, used for indexing into
//                    p_i_values[] array
//
inline int determine_type( char* aname, char* rname )
{
  if( !strcmp( aname, "N" ) && strncmp( rname, "PRO", 3) )
  {
    return 17;
  }
  else if( !strcmp( aname, "CA") && strncmp( rname, "GLY", 3) )
  {
    return 2;
  }
  else if( !strcmp( aname, "C" ) )
  {
    return 5;
  }
  else if( !strcmp( aname, "O" ) )
  {
    return 14;
  }
  else if( !strcmp( aname, "H" ) )
  {
    return 10;
  }
  else if( !strcmp( aname, "OXT" ) )
  {
    return 15;
  }
  //standard amino acid residues
  else if( !strncmp( rname, "ALA", 3) )
  {
    return 4; //CB
  }
  else if( !strncmp( rname, "ARG", 3) )
  {
	    if( !strcmp( aname, "CB") )
      {
        return 3;
      }
	    else if( !strcmp( aname, "CG" ) )
      {
        return 3;
      }
	    else if( !strcmp( aname, "CD" ) )
      {
        return 3;
      }
	    else if( !strcmp( aname, "NE" ) )
      {
        return 17;
      }
	    else if( !strcmp( aname, "CZ" ) )
      {
        return 2;
      }
	    else if( !strcmp( aname, "NH1" ) )
      {
        return 17;
      }
	    else if( !strcmp( aname, "NH2" ) )
      {
        return 17;
      }
	    else if( !strcmp( aname, "HE" ) )
      {
        return 10;
      }
	    else
      {
        return 11; //all HH's
      }
  }
  else if( !strncmp( rname, "ASN", 3) )
  {
	    if( !strcmp( aname, "CB" ) )
      {
        return 3;
      }
	    else if( !strcmp( aname, "CG" ) )
      {
        return 5;
      }
	    else if( !strcmp( aname, "OD1" ) )
      {
        return 14;
      }
	    else if( !strcmp( aname, "ND2" ) )
      {
        return 16;
      }
	    else
      {
        return 10;
      }
  }  
  else if( !strncmp( rname, "ASP", 3) )
  {
	    if( !strcmp( aname, "CB" ) )
      {
        return 3;
      }
	    else if( !strcmp( aname, "CG" ) )
      {
        return 5;
      }
	    else
      {
        return 14; //OD1, OD2
      }
  }    
  else if( !strncmp( rname, "CYS", 3) || !strncmp( rname, "CYX", 3) )
  {
	    if( !strcmp( aname, "CB" ) )
      {
        return 3;
      }
	    else
      {
        return 19; //SG
      }
  }
  else if( !strncmp( rname, "GLN", 3) )
  {
	    if( !strcmp( aname, "CB") )
      {
        return 3;
      }
	    else if( !strcmp( aname, "CG" ) )
      {
        return 3;
      }
	    else if( !strcmp( aname, "CD" ) )
      {
        return 5;
      }
	    else if( !strcmp( aname, "OE1" ) )
      {
        return 14;
      }
	    else if( !strcmp( aname, "NE2" ) )
      {
        return 16;
      }
	    else
      {
        return 10; //1HE1, 1HE2
      }
  }
  else if( !strncmp( rname, "GLU", 3) )
  {
	    if( !strcmp( aname, "CB") )
      {
        return 3;
      }
	    else if( !strcmp( aname, "CG" ) )
      {
        return 3;
      }
	    else if( !strcmp( aname, "CD" ) )
      {
        return 5;
      }
	    else
      {
        return 14;
      }
  }
  else if( !strncmp( rname, "GLY", 3) )
  {
    return 3; //CA
  }
  else if( !strncmp( rname, "HIS", 3) )
  {
	    if( !strcmp( aname, "CB") )
      {
        return 3;
      }
	    else if( !strcmp( aname, "CG" ) )
      {
        return 5;
      }
	    else if( !strcmp( aname, "ND1" ) )
      {
        return 18;
      }
	    else if( !strcmp( aname, "CD2" ) )
      {
        return 6;
      }
	    else if( !strcmp( aname, "CE1" ) )
      {
        return 6;
      }
	    else if( !strcmp( aname, "NE2" ) )
      {
        return 18;
      }
	    else
      {
        return 10; //HD1
      }
  }
  else if( !strncmp( rname, "ILE", 3) )
  {
	    if( !strcmp( aname, "CB" ) )
      {
        return 2;
      }
	    else if( !strcmp( aname, "CG1" ) )
      {
        return 3;
      }
	    else
      {
        return 4; //CG2, CD1
      }
  }     
  else if( !strncmp( rname, "LEU", 3) )
  {
	    if( !strcmp( aname, "CB" ) )
      {
        return 3;
      }
	    else if( !strcmp( aname, "CG" ) )
      {
        return 2;
      }
	    else
      {
        return 4; //CD1, CD2
      }
  }
  else if( !strncmp( rname, "LYS", 3) )
  {
	    if( isH( aname) )
      {
        return 11; //HZ1, HZ2, HZ3
      }
	    else if( !strcmp( aname, "NZ" ) )
      {
        return 16;
      }
	    else
      {
        return 3; //all sidechain carbons
      }
  }     
  else if( !strncmp( rname, "MET", 3) )
  {
	    if( !strcmp( aname, "SD" ) )
      {
        return 19;
      }
	    else if( !strcmp( aname, "CE" ) )
      {
        return 4;
      }
	    else
      {
        return 3; //CB, CG
      }
  }
  else if( !strncmp( rname, "PHE", 3) )
  {
	    if( !strcmp( aname, "CB" ) )
      {
        return 3;
      }
	    else if( !strcmp( aname, "CG" ) )
      {
        return 5;
      }
	    else
      {
        return 6; //all other sidechain carbons
      }
  }
  else if( !strncmp( rname, "PRO", 3) )
  {
	    if( !strcmp( aname, "N" ) )
      {
        return 17;
      }
	    else
      {
        return 3; //all other carbons;
      }
  }
  else if( !strncmp( rname, "SER", 3) )
  {
	    if( !strcmp( aname, "CB" ) )
      {
        return 3;
      }
	    if( !strcmp( aname, "HG" ) )
      {
        return 9;
      }
	    else
      {
        return 13; //OG
      }
  }
  else if( !strncmp( rname, "THR", 3) )
  {
	    if( !strcmp( aname, "CB") )
      {
        return 2;
      }
	    else if( !strcmp( aname, "OG1") )
      {
        return 13;
      }
	    else if( !strcmp( aname, "HG1") )
      {
        return 9;
      }
	    else
      {
        return 4; //CG2
      }
  }     
  else if( !strncmp( rname, "TRP", 3) )
  {
	    if( !strcmp( aname, "CB") )
      {
        return 3;
      }
	    else if( !strcmp( aname, "CG") )
      {
        return 5;
      }
	    else if( !strcmp( aname, "CD1") )
      {
        return 6;
      }
	    else if( !strcmp( aname, "CD2") )
      {
        return 5;
      }
	    else if( !strcmp( aname, "NE1") )
      {
        return 17;
      }
	    else if( !strcmp( aname, "CE2") )
      {
        return 5;
      }
	    else if( !strcmp( aname, "HE1") )
      {
        return 10;
      }
	    else
      {
        return 6; //all other members of ph
      }
  }
  else if( !strncmp( rname, "TYR", 3) )
  {
	    if( !strcmp( aname, "CB") )
      {
        return 3;
      }
	    else if( !strcmp( aname, "CG" ) )
      {
        return 5;
      }
	    else if( !strcmp( aname, "CZ" ) )
      {
        return 5;
      }
	    else if( !strcmp( aname, "OH" ) )
      {
        return 13;
      }
	    else if( !strcmp( aname, "HH" ) )
      {
        return 9;
      }
	    else
      {
        return 6; //all other carbons
      }
  }
  else if( !strncmp( rname, "VAL", 3) )
  {
	    if( !strcmp( aname, "CB") )
      {
        return 4;
      }
	    else
      {
        return 4;//CG1, CG2
      }
  }     

  cudesa_error( "%s(): Cannot determine p_i value for %s:%s\n", __FUNCTION__, rname, aname );
  return 0; //figure out what to do here
}

//
// NUM_P_I_VALUES - number of p_i parameter groups - these are
//                  parameters used to scale SASA based on
//                  multiple overlaps. For more information,
//                  see Section (2)
//
#define NUM_P_I_VALUES 23

//
// p_i_values - actual p_i values, 0.0 values are those of unused p_i groups
//
double p_i_values[ NUM_P_I_VALUES ] =
{
  0.0,
  0.416000,
  0.472000,
  0.544000,
  0.356000,
  0.458000,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.518000,
  0.507000,
  0.577000,
  0.620000,
  0.411000,
  0.411000,
  0.568000,
  0.0,
  0.0,
  0.0,
  0.0
};

//
// radii - atomic radii, indexing corresponds to p_i parameter group indexing
//
double radii[ NUM_P_I_VALUES ] =
{
  0.0, // 1
  2.0, // 2
  2.0, // 3
  2.0, // 4
  1.5, // 5
  1.85,// 6
  0.0, // 7
  0.0, // 8
  0.0, // 9
  0.0, // 10
  0.0, // 11
  0.0, // 12
  1.4, // 13
  1.4, // 14
  1.4, // 15
  1.5, // 16
  1.5, // 17
  1.5, // 18
  1.85,// 19
  0.0, // 20
  0.0, // 21
  0.0, // 22
  0.0  // 23
};

//
// determine_p_i() - Maps numerical p_i type to actual p_i value from
//                   the p_i_values array
//
inline double determine_p_i( int type )
{
#ifdef _DEBUG
  if( (type - 1) < 0 || (type - 1) >= NUM_P_I_VALUES )
  {
    return 1.0;
  }
#endif
  return p_i_values[ type - 1 ];
}

//
// determine_asp_em86() - Determines the atomic solvation parameter
//                        for a given atom according to Eisenberg et
//                        al. For more details, see Section (1.2.1)
//
// input                - aname: atom name
//                        rname: residue name for given atom
//
// output               - ASP value
//
inline double determine_asp_em86( char* aname, char* rname )
{

  static const float EM86_ASP_C = 30.5f;
  static const float EM86_ASP_N_O = -0.9f;
  static const float EM86_ASP_O_CHARGED = -15.0f;
  static const float EM86_ASP_N_CHARGED = -38.5f;
  static const float EM86_ASP_S = 21.0f;
  
  //
  // Charged Oxygen Atoms
  //
  if( !strcmp( rname, "GLU" ) )
  {
    // Resonance: OE2 atoms are assigned "uncharged" Oxygen ASP
    // EM86 indicates that charged ASP should go on most exposed atom
    // but this does not meaningfully change the results, so we simply
    // always assign OE1 the charged ASP
    if( !strcmp( aname, "OE1" ) )
    {
      return EM86_ASP_O_CHARGED;
    }
    if( !strcmp( aname, "OE2" ) )
    {
      return EM86_ASP_N_O;
    }
  }
  if( !strcmp( rname, "ASP" ) )
  {
    // Resonance: OD2 atoms are assigned "uncharged" Oxygen ASP
    // EM86 indicates that charged ASP should go on most exposed atom
    // but this does not meaningfully change the results, so we simply
    // always assign OD1 the charged ASP
    if( !strcmp( aname, "OD1" ) )
    {
      return EM86_ASP_O_CHARGED;
    }
    if( !strcmp( aname, "OD2" ) )
    {
      return EM86_ASP_N_O;
    }
  }

  //
  // Charged Nitrogen Atoms
  //
  if( !strcmp( rname, "LYS" ) )
  {
    if( !strcmp( aname, "NZ" ) )
    {
      return EM86_ASP_N_CHARGED;
    }
  }
  if( !strcmp( rname, "ARG" ) )
  {
    // Resonance: NE and NH1 atoms are assigned "uncharged" Nitrogen ASP
    // EM86 indicates that charged ASP should go on most exposed atom
    // but this does not meaningfully change the results, so we simply
    // always assign NH2 the charged ASP
    if( !strcmp( aname, "NH2" ) )
    {
      return EM86_ASP_N_CHARGED;
    }
    if( !strcmp( aname, "NH1" ) )
    {

      return EM86_ASP_N_O;
    }
    if( !strcmp( aname, "NE" ) )
    {
      return EM86_ASP_N_O;
    }
  }
  if( !strcmp( rname, "HIS" ) )
  {
    // Resonance: ND1 atoms are assigned "uncharged" Nitrogen ASP
    // EM86 indicates that charged ASP should go on most exposed atom
    // but this does not meaningfully change the results, so we simply
    // always assign NE2 the charged ASP
    if( !strcmp( aname, "NE2" ) )
    {
      return EM86_ASP_N_CHARGED;
    }
    if( !strcmp( aname, "ND1" ) )
    {
      return EM86_ASP_N_O;
    }
  }

  //
  // Uncharged atoms
  //

  if( aname[0] == 'N' )
  {
    return EM86_ASP_N_O;
  }
  if( aname[0] == 'O' )
  {
    return EM86_ASP_N_O;
  }

  if( aname[0] == 'C' )
  {
    return EM86_ASP_C;
  }

  if( aname[0] == 'S' )
  {
    return EM86_ASP_S;
  }

  cudesa_error( "%s(): Cannot determine ASP for %s:%s\n", __FUNCTION__, rname, aname );
  return 0.0;
}

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
// myfopen() - This function opens a file for reading or writing and returns
//             a file pointer.  The advantage of using this routine is that
//             it does not crash if the user accidentally types in a non-existent
//             file name.  Instead, it will prompt the user for the correct
//             name.  In addition, if the user tries to open an existing file
//             for writing, s/he will be asked whether to overwrite the existing
//             file with new data.
// 
//             NOTE: No user interaction in parallel mode.  Files open for writing
//             are overwritten if they exist.  Files open for reading that do not exist
//             will cause the program to abort.
//
// input     - filename: file to open
//             iotype: r for reading, w for writing
//
// output    - returns file handle
//
static FILE* myfopen( char *filename, char *iotype )
{
#ifndef PARALLEL
	char newfilename[132];
#endif

	FILE *fp;

	/* open file for writing */
	if( strncmp(iotype,"w",1)==0 )
  {
#ifdef PARALLEL
		/* if we're running parallel, just overwrite the existing filename - no user interaction */
		fp = fopen( filename,"w" );
#else
		fp = fopen( filename,"r" );
		if( fp )
    {
			(void)fclose(fp);
			cudesa_print( CUDESA_PRINT_ALWAYS, "%s(): WARNING:  the file %s exists. Enter a new filename or 'o' to overwrite:  ", __FUNCTION__, filename );
			fscanf( stdin,"%s",newfilename );
			if ( (strncmp(newfilename,"o",1) == 0) && (strlen(newfilename) == 1) )
      {
				fp = fopen( filename, "w" );
			}
      else
      {
				fp = myfopen( newfilename, "w" );
			}
		}
    else
    {
			fp = fopen( filename, "w" );
		}
#endif
	}
  else if( strncmp(iotype,"r",1) == 0 )
  {
		fp = fopen( filename,"r" );
#ifdef PARALLEL
		if (fp==NULL)
    {
			cudesa_error( "%s(): Cannot open file %s... EXITING!", __FUNCTION__, filename );
			ExitandKill();
		}
#else
		if (fp == NULL)
    {
			cudesa_print( CUDESA_PRINT_ALWAYS, "%s(): ERROR:  the file %s cannot be found. Enter a new filename or 'q' to quit:  ", __FUNCTION__, filename );
			fscanf( stdin,"%s",newfilename );
			if ( (strncmp(newfilename,"q",1) == 0) && (strlen(newfilename) == 1) )
      {
				exit( 1 );
			}
      else
      {
				fp = myfopen( newfilename,"r" );
			}
		} 
#endif
	}

	return fp;
}

//
// ExtractString() - extract a substring of a charater string
//
// NOTE: this function was taken from the Rasmol distribution
//
//       infile.c
//       RasMol2 Molecular Graphics
//       Roger Sayle, August 1995
//       Version 2.6
//
static void ExtractString( int len, char *src, char *dst )
{
  register char *ptr;
  register char ch;
  register int i;

  ptr = dst;
  for( i=0; i<len; i++ )
  {
    if( *src )
    {
      ch = *src++;
      *dst++ = ch;
      if( ch != ' ' )
      {
        ptr = dst;
      }
    }
    else
    {
      break;
    }
  }
  *ptr = 0;   
}

//
// assign_access() - Assigns baseline SASA for atoms in molecule. Uses
//                   the method of Hasel et al. described in Section
//                   (1.1.2)
//
// input           - molecule: input molecule
//
// output          - SASA is stored in molecule's atom array
//
void assign_access( PQRData_t *molecule ) 
{
  int i, i_res, ifirst, ilast;
	int j, j_res, jfirst, jlast;
	int k;
	double d, f_ij, p_i, A_i, S_i, r_s = 1.4;
  
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

		  S_i = molecule->atoms[i].r + r_s;
		  S_i = 4.0 * M_PI * S_i * S_i;
		  p_i = determine_p_i( determine_type( molecule->atoms[i].aname, molecule->residues[i_res].rname) );
		  
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

        for( j = jfirst; ((j < jlast) /*&& (j != i)*/); j++ )
        {
          if( j == i )
          {
            continue;
          }

			    d = 0.0;
          
			    for( k = 0; k < 3; k++ )
          {
            d += ( molecule->atoms[i].xyz[k] - molecule->atoms[j].xyz[k] ) * ( molecule->atoms[i].xyz[k] - molecule->atoms[j].xyz[k] );
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

#ifdef _DEBUG
            if( i == DEBUG_ATOM )
            {
              DEBUG_PRINT( CUDESA_PRINT_DEBUG, "*= %f [%d->%d d=%f, cutoff=%f] [ri=%f rj=%f, d=%f, si=%f]"
                           "[%d: <%.3f, %.3f, %.3f> %d: <%.3f, %.3f, %.3f>]\n", f_ij,
                           j, i, d, ( molecule->atoms[j].r + molecule->atoms[i].r + 2.0 * r_s ),
                           molecule->atoms[i].r, molecule->atoms[j].r, d, S_i, j,
                           molecule->atoms[j].xyz[0], molecule->atoms[j].xyz[1], molecule->atoms[j].xyz[2],
                           i,
                           molecule->atoms[i].xyz[0], molecule->atoms[i].xyz[1], molecule->atoms[i].xyz[2] );
            }
#endif
            A_i *= f_ij;
			    }
        } //end j loop
		  } //end j_res loop

		  molecule->atoms[i].access = S_i * A_i;
#ifdef _DEBUG
      if( i == DEBUG_ATOM )
      {
        DEBUG_PRINT( CUDESA_PRINT_DEBUG_EXTRA, "*= %f [S_i]\n", S_i );
        DEBUG_PRINT( CUDESA_PRINT_DEBUG, "= [CPU] %f\n", molecule->atoms[i].access );
      }
#endif      
    } //end i loop
	} //end i_res loop

  molecule->sasa = 0.0;
  for(i = 0; i < molecule->natoms; i++ )
  {
    molecule->sasa += molecule->atoms[ i ].access;
  }
  
} //end assign_access

//
// assign_asp() - Load the atomic solvation parameters for each atom
//                in molecule. For more details on ASPs, see
//                Section (1.2.1)
//
// input        - molecule: input molecule
//
// output       - ASPs are stored in molecule's atom array
//
static void assign_asp( PQRData_t* molecule )
{
  int i_res, ifirst, ilast, i;
  
  molecule->asp = (double*)malloc( sizeof(double) * molecule->natoms );

  if( NULL == molecule->asp )
  {
    cudesa_error( "%s(): Error allocating memory for GPU ASP store\n", __FUNCTION__ );
    return;
  }
  
	for( i_res = 0; i_res < molecule->nres; i_res++ )
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
      molecule->asp[i] = determine_asp_em86( molecule->atoms[i].aname, molecule->residues[i_res].rname );
      if( molecule->asp[ i ] == 0.0 )
      {
        cudesa_error( "%s(): Unknown ASP (0.0) for %s:%s\n", __FUNCTION__, molecule->residues[ i_res ].rname, molecule->atoms[ i ].aname );
      }
    }
  }
}

//
// ReadPQR() - Reads a .pqr file into the PQR data struct. All memory
//             is dynamically allocated as the data is read in.
//
// input     - pqrfile: PQR filename to read in
//             prep_for_gpu: store parameters in parallel arrays in
//             preparation for GPU computation
//
// output    - returns PQRData_t struct representing loaded molecule
//
PQRData_t * ReadPQR( char *pqrfile, int prep_for_gpu )
{
  double  radius;
	double  xmin, ymin, zmin, xmax, ymax, zmax;

	FILE *fp = NULL;
	int acount;                        /* atom counter */
	int rcount;                        /* residue counter */
	int lastresnum = 0;                /* last residue number stored */
  char s_lastresnum[6];
	char pqrstring[132];

	PQRData_t *pqr_p;

  int i, i_res, ifirst, ilast;
  
	/* allocate and initialize new pqr struct  */
	pqr_p = (PQRData_t *) malloc(sizeof(PQRData_t));

  // check for out of memory
  if( NULL == pqr_p )
  {
    return NULL;
  }
  
	memset(pqr_p, 0, sizeof(PQRData_t));

	/* open the .pqr file */
	fp = myfopen(pqrfile, "r");

  // I/O error
  if( NULL == fp )
  {
    return NULL;
  }
  
	/* Fill the structure with data read in from the .pqr file */

	acount = 0;        /* initialize atom & residue counters */
	rcount = 0;

	pqr_p->molradius = 0.0;   /* initialize radius of molecule */
	xmax = -999.9;
	ymax = -999.9;
	zmax = -999.9;
	xmin =  999.9;
	ymin =  999.9;
	zmin =  999.9;


	while (fgets(pqrstring, 132, fp))
  {
		int resnum;
		char s_atomnum[5], s_resnum[6];
		char type[5], atomname[5], resname[5];
		double coordx, coordy, coordz, q, r;
		char s_coordx[10], s_coordy[10], s_coordz[10], s_q[10], s_r[10];
    
		ExtractString(5, &pqrstring[0], type);
		ExtractString(5, &pqrstring[6], s_atomnum);
		ExtractString(4, &pqrstring[12], atomname);
		ExtractString(5, &pqrstring[16], resname);
		ExtractString(5, &pqrstring[22], s_resnum);
		ExtractString(8, &pqrstring[30], s_coordx);
		ExtractString(8, &pqrstring[38], s_coordy);
		ExtractString(8, &pqrstring[46], s_coordz);
		ExtractString(7, &pqrstring[54], s_q);
		ExtractString(6, &pqrstring[62], s_r);

		resnum = atoi(s_resnum);
		coordx = atof(s_coordx);
		coordy = atof(s_coordy);
		coordz = atof(s_coordz);
		q = atof(s_q);
		r = atof(s_r);

		if ((strncmp(type, "ATOM", 4) == 0) || (strncmp(type, "HETATM", 6) == 0))
    {
			atom_t *thisatom;
			char*        p;

			/* Allocate, initialize, then store an atom */
			pqr_p->atoms = (atom_t*) realloc (pqr_p->atoms, sizeof(atom_t) * (acount+1));

			thisatom = &pqr_p->atoms[acount];

			/* get rid of leading & trailing whitespace */
			/* note: may want to use strsep function instead */
			p = strtok (atomname, " \n\t");

			/* Get the atom name */
			strcpy(thisatom->aname, p);
			strcat (thisatom->aname, "\0");

			/* figure out if we need to allocate a new residue_t for this guy */
			//if ((rcount == 0) || (resnum != lastresnum))
      if((rcount == 0) || strcmp( s_resnum, s_lastresnum ) )
      {
				char* res;

				pqr_p->residues = (residue_t*) realloc (pqr_p->residues, sizeof(residue_t) * (rcount+1));
				pqr_p->residues[rcount].firstatom = acount;

				res = strtok (resname, " \n\t");
				strcpy(pqr_p->residues[rcount].rname, res);

				pqr_p->residues[rcount].rnum = resnum;
				/* pqr_p->residues[rcount].ca = -1; */
				pqr_p->residues[rcount].ca = 0;
				pqr_p->residues[rcount].surface = 1;

				pqr_p->residues[rcount].phi = 0.;
				pqr_p->residues[rcount].psi = 0.;

				lastresnum = resnum;
        strcpy( s_lastresnum, s_resnum );
				rcount++;
			}

			/* get the atom coordinates, charge & radius */
			thisatom->xyz[0] = coordx;
			thisatom->xyz[1] = coordy;
			thisatom->xyz[2] = coordz;
			thisatom->XYZ[0] = coordx;
			thisatom->XYZ[1] = coordy;
			thisatom->XYZ[2] = coordz;
			thisatom->q = q;
      // Use indexed Pauling radii instead of what's contained in the file

      thisatom->r = radii[ determine_type( thisatom->aname, pqr_p->residues[ rcount - 1 ].rname ) - 1 ];
			thisatom->access = 0.;

      
			/* store coordinates of the CA atom */

			if ((strncmp(atomname, "CA", 2) == 0) || (strncmp(atomname, " CA", 3) == 0) )
      {
				pqr_p->residues[rcount-1].ca = acount;
			}

			/* sum up the charges */
			pqr_p->tq += thisatom->q;
			acount++;
		}
	}

	pqr_p->natoms = acount;
	pqr_p->nres = rcount;

	fclose(fp);

	cudesa_print( CUDESA_PRINT_DEBUG, "file: %s\tatoms: %d\tresidues: %d\t\n", pqrfile, pqr_p->natoms, pqr_p->nres);

  if( prep_for_gpu )
  {
    assign_asp( pqr_p );
    pqr_p->atomicRadii = (float*)malloc( sizeof(float) * pqr_p->natoms );
    pqr_p->p_iParams = (float*)malloc( sizeof(float) * pqr_p->natoms );
    pqr_p->is_hydrogen = (int*)malloc( sizeof(int) * pqr_p->natoms );
    pqr_p->S_i = (float*)malloc( sizeof(float) * pqr_p->natoms );
    pqr_p->polarity = (int*)malloc( sizeof(int) * pqr_p->natoms );
    
    if( NULL == pqr_p->atomicRadii || NULL == pqr_p->p_iParams || NULL == pqr_p->is_hydrogen )
    {
      cudesa_error( "%s(): Error allocating memory for GPU storage\n", __FUNCTION__ );
      FreePQR( pqr_p );
      return NULL;
    }

    for( i_res = 0; i_res < pqr_p->nres; i_res++ )
    {
      ifirst = pqr_p->residues[i_res].firstatom;
      if( i_res == (pqr_p->nres) - 1 )
      {
        ilast = pqr_p->natoms;
      }
      else
      {
        ilast = pqr_p->residues[i_res + 1].firstatom;
      }
      
      for( i = ifirst; i < ilast; i++ )
      {
        pqr_p->atomicRadii[ i ] = (float)pqr_p->atoms[ i ].r;
        pqr_p->p_iParams[ i ] = (float)determine_p_i( determine_type( pqr_p->atoms[ i ].aname, pqr_p->residues[ i_res ].rname ) );
        pqr_p->is_hydrogen[ i ] = isH( pqr_p->atoms[ i ].aname );
        pqr_p->S_i[ i ] = (float)((pqr_p->atoms[ i ].r + SOLVATION_RADIUS) * (pqr_p->atoms[ i ].r + SOLVATION_RADIUS) * 4.0 * M_PI);

        if( pqr_p->p_iParams[ i ] == 0.0 )
        {
          cudesa_error( "%s(): Invalid p_i for atom %d:%s\n", __FUNCTION__, i, pqr_p->atoms[ i ].aname );
        }

        if( pqr_p->atomicRadii[ i ] == 0.0 )
        {
          cudesa_error( "%s(): Invalid radii for atom %d:%s determine_type( %s, %s ) = %d\n", __FUNCTION__, i, pqr_p->atoms[ i ].aname, pqr_p->atoms[ i ].aname, pqr_p->residues[ i_res ].rname, determine_type( pqr_p->atoms[ i ].aname, pqr_p->residues[ i_res ].rname ) );
        }
      }
    }

  }
  else
  {
    pqr_p->atomicRadii = NULL;
    pqr_p->p_iParams = NULL;
    pqr_p->is_hydrogen = NULL;
    pqr_p->S_i = NULL;
    pqr_p->polarity = NULL;
  }

  strncpy( pqr_p->filename, pqrfile, PQR_FILENAME_MAX - 1 );
  pqr_p->filename[ PQR_FILENAME_MAX - 1 ] = '\0';
  
	return pqr_p;
}

//
// FreePQR() - Frees the memory of the given molecule
//
// input     - molecule: input molecule to free
//
void FreePQR( PQRData_t* molecule )
{
  free( molecule->atoms );
  free( molecule->asp );
  free( molecule->residues );

  free( molecule->atomicRadii );
  free( molecule->p_iParams );
  free( molecule->is_hydrogen );
  free( molecule->S_i );
  free( molecule->polarity );
  
  free( molecule );
}
