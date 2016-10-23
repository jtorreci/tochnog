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

void slide( void )

{
  long int inod=0, max_node=0, idim=0, swit=0, islide=0, max_slide=0,
    in_geometry=0, ldum=0, idum[1], slide_geometry[2];
  double normal_velocity=0., slide_force=0., normal_force=0., slide_friction=0., 
    dtime=0., slide_penalty=1.e15, tmp=0., rdum=0., ddum[MDIM], 
    velocity[MDIM], slide_velocity[MDIM], normal[MDIM], 
    *new_node_dof=NULL, *node_lhside=NULL, *node_rhside=NULL;

  if ( db_max_index( SLIDE_GEOMETRY, max_slide, VERSION_NORMAL, GET ) > 0 ) {
    swit = set_swit(-1,-1,"slide");
    if ( swit ) pri( "In routine SLIDE" );
    for ( islide=0; islide<max_slide; islide++ ) {
      if ( db_active_index( SLIDE_GEOMETRY, islide, VERSION_NORMAL ) ) {
        db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET );
        db( SLIDE_GEOMETRY, islide, slide_geometry, ddum, 
          ldum, VERSION_NORMAL, GET );
        db( SLIDE_FRICTION, inod, idum, &slide_friction, 
          ldum, VERSION_NORMAL, GET_IF_EXISTS );
        db( SLIDE_PENALTY, islide, idum, &slide_penalty, ldum, 
          VERSION_NORMAL, GET_IF_EXISTS );
        db_max_index( NODE, max_node, VERSION_NORMAL, GET );
        for ( inod=0; inod<=max_node; inod++ ) {
          if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
            geometry( inod, ddum, slide_geometry, in_geometry, rdum, normal, 
              rdum, ddum, NODE_START_REFINED, PROJECT_EXACT, VERSION_NORMAL );
            if ( in_geometry ) {
              new_node_dof = db_dbl( NODE_DOF, inod, VERSION_NEW );
              for ( idim=0; idim<ndim; idim++ ) velocity[idim] = 
                new_node_dof[vel_indx+idim*nder];
              normal_velocity = array_inproduct( normal, velocity, ndim );
              for ( idim=0; idim<ndim; idim++ ) {
                tmp = velocity[idim] - normal_velocity * normal[idim];
                slide_velocity[idim] = tmp;
              }
              node_lhside = db_dbl( NODE_LHSIDE, inod, VERSION_NORMAL );
              node_rhside = db_dbl( NODE_RHSIDE, inod, VERSION_NORMAL );
              normal_force = array_inproduct( node_rhside, normal, ndim );
              slide_force = slide_friction * normal_force;
              if ( array_normalize( slide_velocity, ndim ) ) {
                for ( idim=0; idim<ndim; idim++ ) {
                  node_lhside[vel_indx+idim*nder] += slide_penalty * dtime;
                  node_rhside[vel_indx+idim*nder] += 
                    - slide_penalty * normal_velocity * dtime * normal[idim]
                    - slide_force * slide_velocity[idim] ;
                }
              }
            }
          }
        }
      }
    }
    if ( swit ) pri( "Out routine SLIDE" );
  }

}
