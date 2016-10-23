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

void remesh( long int version )

{
  long int inod=0, ino=0, idim=0, max_node=0, length=0,
    swit=0, swit_node=0, in=0, icontrol=0, nnod=0, max_node_boundary=0,
    ldum=0, idum[1], node_remesh_allowed[MDIM],
    *node_node=NULL, *next_of_loop=NULL;
  double dx=0., dtime=0., tmp=0., smallest_distance=0., ddum[1], dx_geometry[MDIM], 
    dx_residue[MDIM], average_coord[MDIM], residue_gradient[MDIM], 
    coord[MDIM], neighbour_coord[MDIM], work[MDIM], 
    node_remesh_velocity[MDIM], control_remesh_factor[2], *node_dof=NULL;

  swit = set_swit(-1,-1,"parallel_remesh");
  if ( swit ) pri( "In routine PARALLEL_REMESH" );

  length = db_data_length(NODE_NODE);
  node_node = get_new_int(length);
  node_dof = get_new_dbl(MUKNWN);

  if ( db_max_index( NODE_BOUNDARY, max_node_boundary, VERSION_NORMAL, GET ) < 0 ) {
    pri( "Error: node_boundary should be specified if you use control_mesh_remesh." );
    exit(TN_EXIT_STATUS);
  }

  db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET_IF_EXISTS );
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  if ( db_active_index( CONTROL_MESH_REMESH_FACTOR, icontrol, VERSION_NORMAL ) ) {
    db( CONTROL_MESH_REMESH_FACTOR, icontrol, idum, 
      control_remesh_factor, ldum, VERSION_NORMAL, GET );
    if ( control_remesh_factor[0]<0. || control_remesh_factor[1]<0. ) 
      db_error( CONTROL_MESH_REMESH_FACTOR, icontrol );
    if ( control_remesh_factor[1]>0. && !residue ) 
      db_error( CONTROL_MESH_REMESH_FACTOR, icontrol );
  }
  else {
    control_remesh_factor[0] = 1.;
    control_remesh_factor[1] = 0.;
  }

  db_max_index( NODE, max_node, version, GET );
  db_version_copy( version, VERSION_TMP );

  array_set( dx_geometry, 0., ndim );
  array_set( dx_residue, 0., ndim );
  for ( inod=0; inod<=max_node; inod++ ) {
    if ( db_active_index( NODE, inod, VERSION_TMP ) ) {
      swit_node = swit; swit = swit && set_swit(-1,inod,"");
      db( NODE, inod, idum, coord, ldum, VERSION_TMP, GET );
      db( NODE_NODE, inod, node_node, ddum, nnod, VERSION_TMP, GET );
      if ( materi_displacement ) {
        db( NODE_DOF, inod, idum, node_dof, nuknwn, VERSION_TMP, GET );
        for ( idim=0; idim<ndim; idim++ )
          coord[idim] += node_dof[dis_indx+idim*nder];
      }
      if ( swit ) {
        pri( "inod", inod );
        pri( "old coord", coord, ndim );
      }
         // relocate
      if ( nnod>0 ) {
          // geometric remeshing displacement
        if      ( control_remesh_factor[0]>0. ) {
          array_set( average_coord, 0., ndim );
          for ( in=0; in<nnod; in++ ) {
            ino = labs(node_node[in]);
            db( NODE, ino, idum, neighbour_coord, ldum, VERSION_TMP, GET );
            if ( materi_displacement ) {
              db( NODE_DOF, ino, idum, node_dof, nuknwn, VERSION_TMP, GET );
              for ( idim=0; idim<ndim; idim++ )
                neighbour_coord[idim] += node_dof[dis_indx+idim*nder];
            }
            if ( swit ) pri( "neighbour_coord-b", neighbour_coord, ndim );
            array_add( average_coord, neighbour_coord, average_coord, 
              ndim );
          }
          array_multiply( average_coord, average_coord, 1./nnod, ndim );
          if ( swit ) pri( "average_coord", average_coord, ndim );
          array_subtract( average_coord, coord, dx_geometry, ndim );
          array_multiply( dx_geometry, dx_geometry, control_remesh_factor[0], ndim );
          if ( swit ) pri( "dx_geometry 1", dx_geometry, ndim );
        }
          // residue remeshing displacement
        if ( control_remesh_factor[1]>0. ) {
          db( NODE_DOF, inod, idum, node_dof, nuknwn, VERSION_TMP, GET );
          array_move( &node_dof[res_indx+1], residue_gradient, ndim );
          array_normalize( residue_gradient, ndim );
          smallest_distance = DBL_MAX;
          for ( in=0; in<nnod; in++ ) {
            ino = labs(node_node[in]);
            db( NODE, ino, idum, neighbour_coord, ldum, VERSION_TMP, GET);
            tmp = array_distance( coord, neighbour_coord, work, ndim );
            if ( tmp<smallest_distance ) smallest_distance = tmp;
          }
          array_multiply( residue_gradient, dx_residue, 0.01*smallest_distance, 
            ndim );
          array_multiply( dx_residue, dx_residue, control_remesh_factor[1], ndim );
        }
        if      ( db_active_index( NODE_REMESH_ALLOWED, inod, VERSION_TMP ) ) {
          db( NODE_REMESH_ALLOWED, inod, node_remesh_allowed, 
            ddum, ldum, VERSION_TMP, GET );
          for ( idim=0; idim<ndim; idim++ ) {
            if ( node_remesh_allowed[idim]==-NO ) {
              dx_geometry[idim] = 0.;
              dx_residue[idim] = 0.;
            }
          }
        }
        else if ( db_active_index( NODE_BOUNDARY, inod, VERSION_TMP ) ) {
          array_set( dx_geometry, 0., ndim );
          array_set( dx_residue, 0., ndim );
        }
        for ( idim=0; idim<ndim; idim++ ) {
          dx = dx_geometry[idim] + dx_residue[idim];
          if ( dtime!=0. ) 
            node_remesh_velocity[idim] = dx/dtime;
          else 
            node_remesh_velocity[idim] = 0.;
          coord[idim] += dx;
          if ( swit ) pri( "new coord becoming", coord, ndim );
        }
        db( NODE, inod, idum, coord, ndim, VERSION_TMP, PUT );
        db( NODE_REMESH_VELOCITY, inod, idum, node_remesh_velocity, ndim, 
          VERSION_TMP, PUT );
        if ( swit ) pri( "new coord", coord, ndim );
      }
      swit = swit_node;
    }
    delete[] next_of_loop;
  }

  db_version_copy( VERSION_TMP, version );
  db_version_delete( VERSION_TMP );

  delete[] node_node;
  delete[] node_dof;

  if ( swit ) pri( "Out routine PARALLEL_REMESH" );
}
