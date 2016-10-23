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

void visco_elasticity( long int element, long int gr, 
  long int formulation, double new_unknowns[],
  double inc_epe[], double rotated_old_msig[], double new_sig[], 
  double new_msig[], double memmat[MDIM][MDIM] )

{

  long int m=0, n=0, swit=0, indx=0, ldum=0, idum[1], task[2];
  double em=0., tm=0., pois=0., dtime=0., ddum[1],
    C[MDIM][MDIM][MDIM][MDIM], maxwell_chain[DATA_ITEM_SIZE], 
    inc_msig[MDIM*MDIM];

  if ( db_active_index(GROUP_MATERI_MAXWELL_CHAIN,gr, VERSION_NORMAL) ) {
    swit = set_swit(element,-1,"visco_elasticity");
    if ( swit ) pri( "In routine VISCO_ELASTICITY" );
    db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET );
    get_group_data( GROUP_MATERI_MAXWELL_CHAIN, gr, element,
      new_unknowns, maxwell_chain, ldum, GET );
    get_group_data( GROUP_MATERI_ELASTI_POISSON, gr, element,
      new_unknowns, &pois, ldum, GET_IF_EXISTS );
    n = db_len( GROUP_MATERI_MAXWELL_CHAIN, gr, VERSION_NORMAL ) / 2;
    for ( m=0; m<n; m++ ) {
      indx = m*MDIM*MDIM;
      em = maxwell_chain[2*m]; tm = maxwell_chain[2*m+1];
      if ( swit ) {
        pri( "maxwell chain number", m );
        pri( "em ", em );
        pri( "tm ", tm );
      }
      visco_elastiticity_nonlinear( m, gr, element, new_unknowns, 
        new_sig, new_msig, em, tm );
      task[0] = GROUP_MATERI_ISOTROPY;
      task[1] = -NO;
      array_set( &C[0][0][0][0], 0., MDIM*MDIM*MDIM*MDIM );
      C_matrix( em, pois, ddum, C, task );
      matrix_a4b( C, inc_epe, inc_msig );
      array_multiply( inc_msig, inc_msig, tm/dtime, MDIM*MDIM );
      array_subtract( inc_msig, &rotated_old_msig[indx], inc_msig, MDIM*MDIM );
      array_multiply( inc_msig, inc_msig, 1.-exp(-dtime/tm), MDIM*MDIM );
      if ( swit ) pri( "inc_msig", inc_msig, MDIM, MDIM );
      array_add( &rotated_old_msig[indx], inc_msig, 
        &new_msig[indx], MDIM*MDIM );
      if ( formulation==TOTAL )
        array_add( new_sig, &new_msig[indx], new_sig, MDIM*MDIM );
      else {
        assert ( formulation==INCREMENTAL );
        array_add( new_sig, inc_msig, new_sig, MDIM*MDIM );
      }
      memmat[0][0] += C[0][0][0][0];
      memmat[0][1] += C[0][0][1][1];
      memmat[0][2] += C[0][0][2][2];
      memmat[1][0] += C[1][1][0][0];
      memmat[1][1] += C[1][1][1][1];
      memmat[1][2] += C[1][1][2][2];
      memmat[2][0] += C[2][2][0][0];
      memmat[2][1] += C[2][2][1][1];
      memmat[2][2] += C[2][2][2][2];
    }
    if ( swit ) {
      pri( "new_sig", new_sig, MDIM, MDIM );
      pri( "Out routine VISCO_ELASTICITY" );
    }
  }

}
