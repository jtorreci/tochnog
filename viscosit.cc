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

void viscous_stress( long int element, long int gr, double user_data[], 
  double unknowns[], double grad_unknowns[], double new_sig[],
  double &viscosity, double &viscosity_heatgeneration )

{
  long int idim=0, jdim=0, ldum=0, group_materi_viscosity_user=-NO,
    group_materi_viscosity_heatgeneration=-NO;
  double ddum[1], D[MDIM*MDIM];

  get_group_data( GROUP_MATERI_VISCOSITY, gr, element,
    unknowns, &viscosity, ldum, GET_IF_EXISTS );
  db( GROUP_MATERI_VISCOSITY_USER, gr, 
    &group_materi_viscosity_user, ddum, ldum, 
    VERSION_NORMAL, GET_IF_EXISTS );
  if ( group_materi_viscosity_user==-YES )
    user_viscosity( user_data, unknowns, viscosity );

  array_set( D, 0., MDIM*MDIM );
  for ( idim=0; idim<ndim; idim++ ) {
    for ( jdim=0; jdim<ndim; jdim++ ) {
      D[idim*MDIM+jdim] = 0.5 *
        ( grad_unknowns[idim*nuknwn+vel_indx+jdim*nder] +
          grad_unknowns[jdim*nuknwn+vel_indx+idim*nder] );
    }
  }                

  array_multiply( D, new_sig, 2.*viscosity, MDIM*MDIM );

  if ( db( GROUP_MATERI_VISCOSITY_HEATGENERATION, gr, 
      &group_materi_viscosity_heatgeneration, ddum, ldum, 
      VERSION_NORMAL, GET_IF_EXISTS ) ) {
    if ( group_materi_viscosity_heatgeneration==-YES )
      viscosity_heatgeneration = 2. * viscosity * array_inproduct( D, D, MDIM*MDIM );
  }

}
