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

void wave( long int element, long int gr, long int nnol,
  double h[], double d[], double volume, double new_unknowns[], 
  double grad_new_unknowns[], double element_lhside[], double element_matrix[],
  double element_rhside[], double element_residue[] )

{

  long int swit=0, jdim=0, inol=0, jnol=0, ipuknwn=0, indx=0, 
    indxi=0, indxj=0, ldum=0, idum[1];
  double tmp=0., c=0., dtime=0.;

  swit = set_swit(element,-1,"wave");
  if ( swit ) pri( "In routine WAVE." );

  db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET );
  get_group_data( GROUP_WAVE_SPEED_OF_SOUND, gr, element,
    new_unknowns, &c, ldum, GET_IF_EXISTS );

  for ( inol=0; inol<nnol; inol++ ) {
    ipuknwn = fscal_indx/nder;
    indx = inol*npuknwn + ipuknwn;
    for ( jdim=0; jdim<ndim; jdim++ ) {
        // diffusive terms (with green partial integration)
      tmp = d[jdim*nnol+inol] * c * c *
        grad_new_unknowns[jdim*nuknwn+scal_indx];
      element_rhside[indx] -= volume * tmp;
      if ( residue ) element_residue[indx] += h[inol] * c * c *
        grad_new_unknowns[jdim*nuknwn+scal_indx+jdim+1];
      indxi = inol*npuknwn + ipuknwn;
      for ( jnol=0; jnol<nnol; jnol++ ) {
        indxj = jnol*npuknwn+ipuknwn;
        element_matrix[indxi*nnol*npuknwn+indxj] += 
          volume * d[jdim*nnol+inol] * c * c * d[jdim*nnol+jnol] * dtime;
        if ( inol==jnol ) element_lhside[indxi] +=
          volume * d[jdim*nnol+inol] * c * c * d[jdim*nnol+jnol] * dtime;
      }
    }
  }

  if ( swit ) {
    pri( "element_rhside", element_rhside, nnol, npuknwn );
    pri( "element_residue", element_residue, nnol, npuknwn );
  }

  if ( swit ) pri( "Out routine WAVE." );

}
