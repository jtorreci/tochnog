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

void merge( void )

{
  long int element=0, merge=0, max_element=0, max_node=0, tmp_max_node=-1,
    tmp_node_number=0, swit=0, inol=0, inod=0, nnol=0, length=0, 
    icontrol=0, in_geometry=0, length_macro_generate=0,
    node_macro_generate=0, equal, ldum=0, 
    idum[1], el[MNOL+1], nodes[MNOL], 
    geometry_entity[2], macro_generate[DATA_ITEM_SIZE], 
    *old_node_numbers=NULL, *ordered_nodes=NULL, 
    *new_node_numbers=NULL, *node_merge_not=NULL;
  double eps_coord=EPS_COORD, rdum=0., ddum[MDIM], node[MDIM];

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  if ( db_active_index( CONTROL_MESH_MERGE, icontrol, VERSION_NORMAL ) ) {
    db( CONTROL_MESH_MERGE, icontrol, &merge, ddum, ldum, VERSION_NORMAL, GET );
    db( CONTROL_MESH_MERGE_EPSCOORD, icontrol, idum, &eps_coord, 
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
    if ( merge==-YES ) {
      swit = set_swit(-1,-1,"merge");
      if ( swit ) pri( "In routine MERGE" );

      db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET ); 
      db_max_index( NODE, max_node, VERSION_NORMAL, GET ); 

      ordered_nodes = get_new_int(1+max_node);
      old_node_numbers = get_new_int(1+max_node);
      new_node_numbers = get_new_int(1+max_node);
      node_merge_not = get_new_int(1+max_node);
      array_set( node_merge_not, 0, 1+max_node );

      if ( db( CONTROL_MESH_MERGE_NOT, icontrol, geometry_entity, 
          ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS ) ) {
        for ( inod=0; inod<=max_node; inod++ ) {
          if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
            geometry( inod, ddum, geometry_entity, in_geometry, rdum, ddum, rdum,
              ddum, NODE_START_REFINED, PROJECT_EXACT, VERSION_NORMAL );
            if ( in_geometry ) node_merge_not[inod] = 1;
          }
        }
      }

      if ( db( CONTROL_MESH_MERGE_MACRO_GENERATE, icontrol, macro_generate, 
          ddum, length_macro_generate, VERSION_NORMAL, GET_IF_EXISTS ) ) {
        for ( inod=0; inod<=max_node; inod++ ) {
          if ( db_active_index( NODE, inod, VERSION_NORMAL ) && node_merge_not[inod]==0 ) {
            node_macro_generate = -1;
            db( NODE_MACRO_GENERATE, inod, &node_macro_generate, 
              ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
            if ( node_macro_generate==-1 ) 
              node_merge_not[inod] = 1;
            else if ( !array_member(macro_generate,node_macro_generate,length_macro_generate,ldum) ) 
              node_merge_not[inod] = 1;
            else
              node_merge_not[inod] = 0;
          }
        }
      }

        // generate nodes; initialise
      array_set( new_node_numbers, -1, 1+max_node );

        // first generate nodes which can be merged
      for ( inod=0; inod<=max_node; inod++ ) {
        if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
          if ( !node_merge_not[inod] ) {
            db( NODE, inod, idum, node, ldum, VERSION_NORMAL, GET );
            ordered_list_apply( inod, ordered_nodes, tmp_max_node, 
              node, eps_coord, tmp_node_number, ADD );
            if ( tmp_node_number<0 ) {
              tmp_max_node++;
              old_node_numbers[tmp_max_node] = inod;
              db( NODE, tmp_max_node, idum, node, ndim, VERSION_TMP, PUT );
            }
            else {
              new_node_numbers[inod] = old_node_numbers[tmp_node_number];
              delete_node( inod, VERSION_NORMAL );
            }
          }
        }
      }

        // second generate nodes which can not be merged
      for ( inod=0; inod<=max_node; inod++ ) {
        if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
          if ( node_merge_not[inod] ) {
            db( NODE, inod, idum, node, ldum, VERSION_NORMAL, GET );
            tmp_node_number = -1;
            ordered_list_apply( inod, ordered_nodes, tmp_max_node, node, 
              eps_coord, equal, ADD_ALWAYS );
            if ( tmp_node_number<0 ) {
              tmp_max_node++;
              old_node_numbers[tmp_max_node] = inod;
              db( NODE, tmp_max_node, idum, node, ndim, VERSION_TMP, PUT );
            }
            else {
              new_node_numbers[inod] = old_node_numbers[tmp_node_number];
              delete_node( inod, VERSION_NORMAL );
            }
          }
        }
      }

      for ( element=0; element<=max_element; element++ ) {
        if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) ) {
          db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
          nnol = length - 1; array_move( &el[1], nodes, nnol );
          for ( inol=0; inol<nnol; inol++ ) {
            inod = nodes[inol];
            if ( new_node_numbers[inod]>=0 )
             nodes[inol] = new_node_numbers[inod];
          }
          array_move( nodes, &el[1], nnol );
          db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, PUT );
        }
      }

      delete[] ordered_nodes;
      delete[] old_node_numbers;
      delete[] new_node_numbers;
      delete[] node_merge_not;

      db_version_delete( VERSION_TMP );

      mesh_has_changed( VERSION_NORMAL );

      if ( swit ) pri( "Out routine MERGE" );
    }
  }
}
