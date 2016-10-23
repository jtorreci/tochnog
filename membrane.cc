/*
    Copyright (C) 1998  Dennis Roddeman
    email: dennis.roddeman@feat.nl

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software Foundation 
    59 Temple Place, Suite 330, Boston, MA, 02111-1307, USA
*/

#include "tochnog.h"

#define EPS_STRAIN_CORRECTION 1.e-10

long int membrane_apply( long int element, long int gr, 
  double memmat[MDIM][MDIM], double inc_ept[], double inc_epe[], 
  double new_ept[], double new_epe[], double new_sig[] )

  /* Assumption: no shear stresses due to normal strains */

{
  long int swit=0, membrane_found=1, membrane=-NO, axisymmetric=-NO, ldum=0;
  double rdum=0., ddum[1], strain_correction[MDIM*MDIM], rh[2], mat[2*2], inv_mat[2*2], 
    work[MDIM];

      // elastic correction for zero perpendicular stress(es)
  db( GROUP_MATERI_MEMBRANE, gr, &membrane, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  if ( membrane==-YES ) {
    swit = set_swit(element,-1, "membrane");
    if ( swit ) pri( "In routine MEMBRANE" );
    db( GROUP_AXISYMMETRIC, gr, &axisymmetric, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    if ( axisymmetric==-YES ) {
      pri( "Error: either use GROUP_MATERI_MEMBRANE or GROUP_AXISYMMETRIC." );
      exit(TN_EXIT_STATUS);
    }
    if ( db_active_index( GROUP_MATERI_ELASTI_COMPRESSIBILITY, gr, VERSION_NORMAL ) ) {
      pri( "Error: GROUP_MATERI_ELASTI_COMPRESSIBILITY cannot be used with MEMBRANE." );
      exit(TN_EXIT_STATUS);
    }
    if ( db_active_index( GROUP_MATERI_ELASTI_COMPRESSIBILITY, gr, VERSION_NORMAL ) ) {
      pri( "Error: GROUP_MATERI_ELASTI_COMPRESSIBILITY cannot be used with MEMBRANE." );
      exit(TN_EXIT_STATUS);
    }
    mat[0] = memmat[1][1];
    mat[1] = memmat[1][2];
    mat[2] = memmat[2][1];
    mat[3] = memmat[2][2];
    array_set( work, 0., MDIM );
    array_set( strain_correction, 0., MDIM*MDIM );
    if ( ndim==2 ) {
      if ( mat[3]!=0. ) work[2] = -new_sig[8]/mat[3];
    }
    else if ( ndim==1 ) {
      rh[0] = -new_sig[4];
      rh[1] = -new_sig[8];
      if ( matrix_inverse( mat, inv_mat, rdum, 2 ) )
        matrix_ab( inv_mat, rh, &work[1], 2, 2, 1 );
    }
    strain_correction[4] = work[1];
    strain_correction[8] = work[2];
    array_add( inc_ept, strain_correction, inc_ept, MDIM*MDIM );
    array_add( inc_epe, strain_correction, inc_epe, MDIM*MDIM );
    array_add( new_ept, strain_correction, new_ept, MDIM*MDIM );
    array_add( new_epe, strain_correction, new_epe, MDIM*MDIM );
    if ( array_size( strain_correction, MDIM*MDIM ) > EPS_STRAIN_CORRECTION ) 
      membrane_found = 0;
    if ( swit ) {
      pri( "strain correction", strain_correction, MDIM, MDIM );
      pri( "Out routine MEMBRANE" );
    }
  }

  return membrane_found;
}
