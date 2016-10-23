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

#define EPS_WEIGHT 1.e-3

void parallel_map_element( void )

{
  long int element=0, max_map_element=0, iloop=0, nloop=0, swit=0,
    ithread=0, *next_of_loop=NULL;

  swit = set_swit(-1,-1,"parallel_map_element");
  if ( swit ) pri( "In routine PARALLEL_MAP_ELEMENT" );

    // loop over elements
  db_max_index( ELEMENT, max_map_element, map_version_to, GET );
  if ( max_map_element>=0 ) {
    next_of_loop = get_new_int(1+max_map_element);
    parallel_sys_next_of_loop( next_of_loop, max_map_element, nloop, ithread );
    for ( iloop=0; iloop<nloop; iloop++ ) {
      element = next_of_loop[iloop];
      if ( element>max_map_element )
        break;
      else if ( db_active_index( ELEMENT, element, map_version_to ) )
        map_element( element );
    }
    delete[] next_of_loop;
  }

  if ( swit ) pri( "Out routine PARALLEL_MAP_ELEMENT" );
}

void parallel_map_node( void )

{
  long int inod=0, max_map_node=0, iloop=0, nloop=0, swit=0,
    ithread=0, *next_of_loop=NULL;

  swit = set_swit(-1,-1,"parallel_map_node");
  if ( swit ) pri( "In routine PARALLEL_MAP_NODE" );

    // loop over elements
  db_max_index( NODE, max_map_node, map_version_to, GET );
  if ( max_map_node>=0 ) {
    next_of_loop = get_new_int(1+max_map_node);
    parallel_sys_next_of_loop( next_of_loop, max_map_node, nloop, ithread );
    for ( iloop=0; iloop<nloop; iloop++ ) {
      inod = next_of_loop[iloop];
      if ( inod>max_map_node )
        break;
      else if ( db_active_index( NODE, inod, map_version_to ) )
        map_node( inod );
    }
    delete[] next_of_loop;
  }

  if ( swit ) pri( "Out routine PARALLEL_MAP_NODE" );
}

void map_element( long int element )

{
  long int max_old_element=0, found=0, length=0, name=0, in=0, 
    iel=0, inol=0, old_nnol=0, tmp_nnol=0, idim=0,
    element_smallest_distance=0, use_smallest_distance=0,
    swit=0, ldum=0, idum[1], old_el[MNOL+1], old_nodes[MNOL], 
    tmp_el[MNOL+1], tmp_nodes[MNOL];
  double distance=0., smallest_distance=DBL_MAX, 
    ddum[1], coord[MDIM], coords[MDIM*MNOL],
    average_coord[MDIM], weight[MNOL], tmp_coord[MDIM],
    work[MDIM];

  swit = set_swit(-1,-1,"map_element");
  if ( swit ) pri( "In routine MAP_ELEMENT" );

  start_of_map_element:

  db( ELEMENT, element, tmp_el, ddum, length, map_version_to, GET );
  tmp_nnol = length - 1; array_move( &tmp_el[1], tmp_nodes, tmp_nnol );
  array_set( average_coord, 0., ndim );
  for ( inol=0; inol<tmp_nnol; inol++ ) {
    in = tmp_nodes[inol];
    db( NODE, in, idum, coord, ldum, map_version_to, GET );
    for ( idim=0; idim<ndim; idim++ ) 
      average_coord[idim] += coord[idim]/tmp_nnol;
  }
  found = 0;
  db_max_index( ELEMENT, max_old_element, map_version_from, GET );
  for ( iel=0; iel<=max_old_element && !found; iel++ ) {
    if ( db_active_index( ELEMENT, iel, map_version_from ) ) {
      db( ELEMENT, iel, old_el, ddum, length, map_version_from, GET );
      name = old_el[0];
      old_nnol = length - 1; array_move( &old_el[1], old_nodes, old_nnol );
      for ( inol=0; inol<old_nnol; inol++ ) {
        in = old_nodes[inol];
        db( NODE, in, idum, &coords[inol*ndim], ldum, map_version_from, GET );
      }
      if ( (point_el(average_coord,coords,weight,name,old_nnol) && 
            !use_smallest_distance) ||
           (iel==element_smallest_distance && use_smallest_distance) 
         ) {
        found = 1;
        length = 1 + tmp_nnol;
        create_element( iel, element, tmp_el, length,
          map_version_from, map_version_to );
      }
      else {
        array_set( tmp_coord, 0., ndim );
        for ( inol=0; inol<old_nnol; inol++ ) {
          for ( idim=0; idim<ndim; idim++ ) tmp_coord[idim] += 
            coords[inol*ndim+idim]/old_nnol;
        }
        distance = array_distance( coord, tmp_coord, work, ndim );
        if ( distance<smallest_distance ) {
          smallest_distance = distance;
          element_smallest_distance = iel;
        }
      }
    }
  }
  if ( !found && map_always ) {
    use_smallest_distance = 1;
    goto start_of_map_element;
  }

  if ( swit ) pri( "Out routine MAP_ELEMENT" );
}

void map_node( long int inod )

{
  long int max_old_element=0, found=0, length=0, name=0, in=0, 
    iel=0, inol=0, old_nnol=0, idim=0, iuknwn=0, tmp_node_number=0,
    element_smallest_distance=0, use_smallest_distance=0,
    ldum=0, weight_nnol=0, idum[1], old_el[MNOL+1],
    old_nodes[MNOL], weight_nodes[MNOL];
  double total_weight=0., distance=0., smallest_distance=DBL_MAX,
    ddum[1], weight[MNOL], we[MNOL], coord[MDIM], coords[MDIM*MNOL],
    node[MDIM], node_dof[MUKNWN], 
    node_start_refined[MDIM], node_dof_start_refined[MUKNWN],
    tmp_node[MDIM], tmp_node_dof[MUKNWN], 
    tmp_node_start_refined[MDIM], tmp_node_dof_start_refined[MUKNWN],
    tmp_coord[MDIM], work[MDIM];

  start_of_map_node:

  found = 0;
  db( NODE, inod, idum, coord, ndim, map_version_to, GET );
  db_max_index( ELEMENT, max_old_element, map_version_from, GET );
  for ( iel=0; iel<=max_old_element && !found; iel++ ) {
    if ( db_active_index( ELEMENT, iel, map_version_from ) ) {
      db( ELEMENT, iel, old_el, ddum, length, map_version_from, GET );
      name = old_el[0];
      old_nnol = length - 1; array_move( &old_el[1], old_nodes, old_nnol );
      for ( inol=0; inol<old_nnol; inol++ ) {
        in = old_nodes[inol];
        db( NODE, in, idum, &coords[inol*ndim], ldum, map_version_from, GET );
      }
      if ( (point_el(coord,coords,we,name,old_nnol) && !use_smallest_distance) ||
           (iel==element_smallest_distance && use_smallest_distance) 
         ) {
        found = 1;
          // only use nodes which really contribute something
        total_weight = 0.; weight_nnol = 0;
        array_set( weight, 0., old_nnol );
        for ( inol=0; inol<old_nnol; inol++ ) {
          if (scalar_dabs(we[inol])>EPS_WEIGHT ) {
            weight[weight_nnol] = we[inol];
            weight_nodes[weight_nnol] = old_nodes[inol];
            weight_nnol++;
            total_weight += we[inol];
          }
        }
        if ( total_weight>0. ) {
          for ( inol=0; inol<weight_nnol; inol++ ) weight[inol] /= total_weight;
        }
        array_set( tmp_node, 0., ndim );
        array_set( tmp_node_dof, 0., nuknwn );
        array_set( tmp_node_start_refined, 0., ndim );
        array_set( tmp_node_dof_start_refined, 0., nuknwn );
        for ( inol=0; inol<weight_nnol; inol++ ) {
          in = weight_nodes[inol];
          db( NODE, in, idum, node, ldum, map_version_from, GET );
          db( NODE_START_REFINED, in, idum, node_start_refined, 
            ldum, map_version_from, GET );
          for ( idim=0; idim<ndim; idim++ ) {
            tmp_node[idim] += weight[inol] * node[idim];
            tmp_node_start_refined[idim] += weight[inol] * 
              node_start_refined[idim];
          }
          if ( nuknwn>0 ) {
            db( NODE_DOF, in, idum, node_dof, ldum, map_version_from, GET );
            db( NODE_DOF_START_REFINED, in, idum, node_dof_start_refined, 
              ldum, map_version_from, GET );
            for ( iuknwn=0; iuknwn<nuknwn; iuknwn++ ) {
              tmp_node_dof[iuknwn] += weight[inol] * node_dof[iuknwn];
              tmp_node_dof_start_refined[iuknwn] += weight[inol] * 
                node_dof_start_refined[iuknwn];
            }
          }
        }
        tmp_node_number = inod;
        create_node( weight_nodes, weight_nnol, tmp_node_number, tmp_node,
          tmp_node_dof, tmp_node_start_refined,
          tmp_node_dof_start_refined,
          map_version_from, map_version_to );
      }
      else {
        array_set( tmp_coord, 0., ndim );
        for ( inol=0; inol<old_nnol; inol++ ) {
          for ( idim=0; idim<ndim; idim++ ) tmp_coord[idim] += 
            coords[inol*ndim+idim]/old_nnol;
        }
        distance = array_distance( coord, tmp_coord, work, ndim );
        if ( distance<smallest_distance ) {
          smallest_distance = distance;
          element_smallest_distance = iel;
        }
      }
    }
  }
  if ( !found && map_always ) {
    use_smallest_distance = 1;
    goto start_of_map_node;
  }

}
