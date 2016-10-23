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


void ordered_list_apply( long int inn, long int ordered_coords[], 
  long int tmp_max_coord, double tmp_coord[], 
  double eps_coord, long int &equal, long int task ) 

{
  long int inod=0, jnod=0, ready=0, lower=0, higher=0, smaller=0;
  double tmp=0., radius=0., tmp_radius=0., diff_coord[MDIM], *tmp_co=NULL;

    // initialise
  lower = -1;
  equal = -1;
  higher = -1;

    // order by radius
  if ( tmp_max_coord==-1 )
    ordered_coords[0] = 0; 
  else {
    ready = 0; lower = -1; higher = tmp_max_coord + 1;
    tmp_radius = array_size( tmp_coord, ndim );
    while ( !ready ) {
      inod = ( lower + higher ) / 2;
      jnod = ordered_coords[inod];
      tmp_co = db_dbl( NODE, jnod, VERSION_TMP );
      radius = array_size( tmp_co, ndim );
      if ( lower==higher-1 ) {
        ready = 1;
      }
      else if ( radius<0.9*tmp_radius-10.*eps_coord ) {
        lower = inod;
      }
      else if ( radius>=1.1*tmp_radius+10.*eps_coord ) {
        higher = inod;
      }
      else {
        ready = 1;
      }
    }

      // just to be sure
    lower -= 200;
    higher += 200;

      // check if coordinate already exists
    for ( inod=lower; inod<=higher; inod++ ) {
      if ( inod>=0 && inod<=tmp_max_coord ) {
        jnod = ordered_coords[inod];
        tmp_co = db_dbl( NODE, jnod, VERSION_TMP );
        radius = array_size( tmp_co, ndim );
        if ( radius<tmp_radius-eps_coord ) smaller = inod;
        array_subtract( tmp_co, tmp_coord, diff_coord, ndim );
        tmp = array_size( diff_coord, ndim );
        if ( tmp<eps_coord ) equal = jnod;
      }
    }

      // add to ordered list
    if      ( task==ADD_ALWAYS ) {
      for ( inod=tmp_max_coord; inod>smaller; inod-- )
        ordered_coords[inod+1] = ordered_coords[inod];
      ordered_coords[smaller+1] = tmp_max_coord + 1;
    }
    else if ( task==ADD ) {
      if ( equal<0 ) {
        for ( inod=tmp_max_coord; inod>smaller; inod-- ) {
          ordered_coords[inod+1] = ordered_coords[inod];
        }
        ordered_coords[smaller+1] = tmp_max_coord + 1;
      }
    }
    else {
      assert( task==CHECK );
    }

  }

}
