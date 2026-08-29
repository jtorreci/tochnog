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

void volume_factor( long int element_group, double coord[], double &volfac )

{
  int i=0, level=0;
  long int j=0, length=0, ldum=0, idum[1];
  double x=0., y=0., group_volume_factor=1.,
    factor=1., *ptr=NULL;

  if ( db_active_index( VOLUME_FACTOR, 0, VERSION_NORMAL ) ) {
    length = db_len( VOLUME_FACTOR, 0, VERSION_NORMAL );
    ptr = db_dbl( VOLUME_FACTOR, 0, VERSION_NORMAL );
    if ( ndim==1 ) {
      x = coord[0];
      factor = ptr[0];
      for ( j=level=1; j<length; level++ ) {
        factor += ptr[j] * scalar_power(x,level); j++;
      }
    }
    else {
      assert( ndim==2 );
      x = coord[0];
      y = coord[1];
      j = 0;
      factor = ptr[j]; j++;
      for ( level=1; j<length; level++ ) {
        if ( j<length ) {
          factor += ptr[j] * scalar_power(x,level); j++;
        }
        for ( i=level-1; i>0 && j<length; i-- ) {
          factor += ptr[j] * scalar_power(x,i) * scalar_power(y,level-i); j++;
        }
        if ( j<length ) {
          factor += ptr[j] * scalar_power(y,level); j++;
        }
      }
    }
  }
  volfac = factor;

  // volume_factor_x (manual Professional 6.1094): an in x-direction
  // piecewise-constant volume factor: left from x0 the factor is 1,
  // from x0 to x1 it is fac01, etcetera, and right from the last x the
  // factor is 1 again. The record is x0 fac01 x1 fac12 ... xn (an odd
  // number of values). It multiplies the polynomial volume_factor and
  // the group_volume_factor.
  if ( db_active_index( VOLUME_FACTOR_X, 0, VERSION_NORMAL ) ) {
    length = db_len( VOLUME_FACTOR_X, 0, VERSION_NORMAL );
    ptr = db_dbl( VOLUME_FACTOR_X, 0, VERSION_NORMAL );
    if ( length>=3 && length%2==1 ) {
      x = coord[0];
      factor = 1.;
      // inside [x0, xn): the factor of the containing interval; left
      // from x0 and right from xn: 1 (manual 6.1094)
      if ( x>=ptr[0] && x<ptr[length-1] ) {
        for ( j=0; j<length-1; j+=2 ) {
          if ( x>=ptr[j] ) factor = ptr[j+1];
          else break;
        }
      }
      volfac *= factor;
    }
    else {
      pri( "Warning: volume_factor_x needs an odd number of values "
           "(x0 fac01 x1 fac12 ... xn) - the record is neglected" );
    }
  }

  db( GROUP_VOLUME_FACTOR, element_group, idum, &group_volume_factor, 
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  volfac *= group_volume_factor;

}
