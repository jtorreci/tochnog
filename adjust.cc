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

void adjust_geom( long int geometry_entity[], 
  long int geometry_entity_edge[] )

{
  long int element=0, max_element=0, inol=0, nnol=0, inod=0, length=0,
    all_in_geometry=0, some_in_geometry=0, in_geometry=0,
    node_boundary=-YES, swit=0, ldum=0, 
    idum[1], *nodes=NULL, *el=NULL, *node_in_geometry=NULL;
  double rdum=0., ddum[MDIM], coord[MDIM], 
    diff_coord[MDIM], node_start_refined[MDIM], projection[MDIM], 
    normal_dum[MDIM];

  area_node_dataitem();

  swit = set_swit(-1,-1,"adjust_geom");
  if ( swit ) pri( "In routine ADJUST_GEOM" );

  nodes = get_new_int(MNOL);
  el = get_new_int(MNOL+1);
  node_in_geometry = get_new_int(MNOL);

  db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
  if ( geometry_entity[0]==-GEOMETRY_SET ) {
    cout << "Sorry, a GEOMETRY_SET cannot be used for adjusting geometry.";
    exit(TN_EXIT_STATUS);
  }

  for ( element=0; element<=max_element; element++ ) {
    if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) &&
         !db_active_index( ELEMENT_ADJUST, element, VERSION_NORMAL ) ) {
      db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
      nnol = length - 1; array_move( &el[1], nodes, nnol );
      some_in_geometry = 0; all_in_geometry = 1; 
      array_set( node_in_geometry, 0., nnol );
      for ( inol=0; inol<nnol; inol++ ) {
        inod = nodes[inol];
        geometry( inod, ddum, geometry_entity, in_geometry, rdum, 
          normal_dum, rdum, projection, NODE_START_REFINED, PROJECT_EXACT, 
          VERSION_NORMAL );
        if ( in_geometry ) {
          some_in_geometry = 1;
          node_in_geometry[inol] = 1;
        }
        else
          all_in_geometry = 0;
      }
      if ( some_in_geometry && !all_in_geometry ) {
        for ( inol=0; inol<nnol; inol++ ) {
          inod = nodes[inol];
          if ( node_in_geometry[inol] &&
               !db_active_index( NODE_ADJUST, inod, VERSION_NORMAL ) ) {
            geometry( inod, ddum, geometry_entity_edge, in_geometry, rdum, 
              normal_dum, rdum, projection, NODE_START_REFINED, PROJECT_EXACT,
              VERSION_NORMAL );
            db( NODE_START_REFINED, inod, idum, node_start_refined, 
              ldum, VERSION_NORMAL, GET );
            db( NODE, inod, idum, coord, ldum, VERSION_NORMAL, GET );
            array_subtract( projection, node_start_refined, diff_coord, ndim );
            array_add( coord, diff_coord, coord, ndim );
            length = 1; 
            db( NODE, inod, idum, coord, ndim, VERSION_NORMAL, PUT );
            db( NODE_START_REFINED, inod, idum, projection, 
              ndim, VERSION_NORMAL, PUT );
            db( NODE_BOUNDARY, inod, &node_boundary, ddum, length,
               VERSION_NORMAL, PUT );
          }
        }
      }
    }
  }
  mesh_has_changed( VERSION_NORMAL );

  delete[] nodes;
  delete[] el;
  delete[] node_in_geometry;

  if ( swit ) pri( "Out routine ADJUST_GEOM" );

}
