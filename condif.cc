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

void condif( long int element, long int gr, long int nnol, double h[], 
  double volume, double new_unknowns[], double element_lhside[],
  double element_matrix[], double element_rhside[], double element_residue[] )

{

  long int swit=0, inol=0, ipuknwn=0, indx=0, ldum=0;
  double condif_absorption=0.;

  if ( get_group_data( GROUP_CONDIF_ABSORPTION, gr, element, new_unknowns, 
      &condif_absorption, ldum, GET_IF_EXISTS ) ) {
    swit = set_swit(element,-1,"condif");
    if ( swit ) pri( "In routine CONDIF" );
      // absorption
    for ( inol=0; inol<nnol; inol++ ) {
      ipuknwn = temp_indx/nder;
      indx = inol*npuknwn + ipuknwn;
      element_rhside[indx] -= volume * h[inol] * condif_absorption *
        new_unknowns[temp_indx];
      if ( residue ) element_residue[indx] -= h[inol] * condif_absorption *
        new_unknowns[temp_indx];
      element_lhside[indx] += volume * h[inol] * condif_absorption;
      element_matrix[indx*nnol*npuknwn+indx] += 
        volume * h[inol] * condif_absorption;
    }
    if ( swit ) {
      pri( "element_lhside", element_lhside, nnol, npuknwn );
      pri( "element_rhside", element_rhside, nnol, npuknwn );
      if ( residue ) pri( "element_residue", element_residue, nnol, npuknwn );
      pri( "Out routine CONDIF" );
    }
  }

}
