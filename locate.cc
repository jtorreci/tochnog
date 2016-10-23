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

void locate( void )

{
  long int inod=0, max_node=0, idim=0, ldum=0, length=0, swit=0, 
    idum[1], options_mesh[MDIM];
  double dtime=0., ddum[1], flow[MDIM],
    *new_node_dof=NULL, *tmp_node=NULL;

  db_max_index( NODE, max_node, VERSION_NORMAL, GET );
  if ( materi_velocity && !materi_displacement && max_node>=0 ) {
    swit = set_swit(-1,-1,"locate");
    if ( swit ) pri( "In routine LOCATE" );
    db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET_IF_EXISTS );
    db( OPTIONS_MESH, 0, options_mesh, ddum, ldum, VERSION_NORMAL, GET);
    length = ( max_node + 1 ) * ndim;
    tmp_node = get_new_dbl( length );
    array_set( tmp_node, DBL_MAX, length );
    array_move( db_dbl(NODE,VERSION_NORMAL), tmp_node, length );
    for ( inod=0; inod<=max_node; inod++ ) {
      if ( db_active_index( NODE_DOF, inod, VERSION_NEW ) ) {
        new_node_dof = db_dbl( NODE_DOF, inod, VERSION_NEW );
        for ( idim=0; idim<ndim; idim++ ) 
          flow[idim] = new_node_dof[vel_indx+idim*nder];
        for ( idim=0; idim<ndim; idim++ ) {
          if ( options_mesh[idim]==-FOLLOW_MATERIAL )
            tmp_node[inod*ndim+idim] += flow[idim] * dtime;
        }
      }
    }
    array_move( tmp_node, db_dbl(NODE,VERSION_NEW), length );
    delete[] tmp_node;
    if ( swit ) pri( "Out routine LOCATE" );
  }

}
