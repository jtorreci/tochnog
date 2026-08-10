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

void groundflow_data( long int element, long int gr, double old_unknowns[], 
  double new_unknowns[], double coord_ip[], double pe[], double &C,
  double h[], long int nnol )

{
  long int ldum=0, idum=0, inol=0, idim=0, jdim=0, vert=0, nuknwn=0;
  double ddum[1], pvs[5], sigv=0., tmp=0., sig=0., pres=0., vertmax=0.;
  double force_gravity[MDIM];

  nuknwn = npuknwn * nder;

  get_group_data( GROUP_GROUNDFLOW_PERMEABILITY, gr, element, new_unknowns,
    pe, ldum, GET_IF_EXISTS );
  get_group_data( GROUP_GROUNDFLOW_CAPACITY, gr, element, new_unknowns, 
    &C, ldum, GET_IF_EXISTS );

  // group_groundflow_permeability_vertical_stress: kp = a / (sigv/sig0)^b,
  // clamped to [minimum, maximum]. sigv = vertical EFFECTIVE stress.
  // Only used in 2D/3D; falls back to group_groundflow_permeability when the
  // vertical stress is extremely small (avoid division by zero).
  if ( db_active_index( GROUP_GROUNDFLOW_PERMEABILITY_VERTICAL_STRESS, gr, VERSION_NORMAL ) ) {
    ldum = 5;
    get_group_data( GROUP_GROUNDFLOW_PERMEABILITY_VERTICAL_STRESS, gr, element,
      new_unknowns, pvs, ldum, GET_AND_CHECK );
    // pvs[0]=a, pvs[1]=b, pvs[2]=sig0, pvs[3]=minimum, pvs[4]=maximum

    force_gravity_calculate( force_gravity );

    // vertical direction = component of gravity with largest magnitude
    vert = 0; vertmax = fabs(force_gravity[0]);
    for ( idim=1; idim<ndim; idim++ ) {
      if ( fabs(force_gravity[idim]) > vertmax ) { vertmax = fabs(force_gravity[idim]); vert = idim; }
    }

    // interpolate effective vertical stress to the integration point
    sigv = 0.;
    for ( inol=0; inol<nnol; inol++ ) {
      sig = new_unknowns[ inol*npuknwn*nder + stres_indx + stress_indx(vert,vert)*nder ];
      pres = new_unknowns[ inol*npuknwn*nder + pres_indx ];
      sigv += h[inol] * ( sig - pres );   // effective vertical stress
    }
    if ( fabs(sigv) > 1.e-10 ) {
      tmp = pvs[0] / pow( fabs(sigv)/pvs[2], pvs[1] );
      if ( tmp < pvs[3] ) tmp = pvs[3];
      if ( tmp > pvs[4] ) tmp = pvs[4];
      for ( idim=0; idim<ndim; idim++ ) pe[idim] = tmp;
    }
    // else: keep pe from group_groundflow_permeability (tiny vertical stress)
  }

  (void)old_unknowns; (void)coord_ip; (void)ddum; (void)idum;
}
