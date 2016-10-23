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

#define MSPLIT_ELEMENT 12
#define MSPLIT_NODE 10
#define MNEW_NODE 4
#define MNEW_NODE_LIST 8

void mesh_split( long int version )

{
  long int swit=0, max_node=0, max_element=0,
    name=0, inod=0, inol=0, jnol=0, nnol=0, length=0, length_element=0, 
    tmp_max_element=0, isplit_element=0, isplit_node=0,
    nsplit_element=0, nsplit_node=0, nnew_node=0, 
    nnew_number=0, element=0, icontrol=0, in_geometry=0,
    all_in_geometry=0, ldum=0, iwork[1], control_mesh_split_only[2], 
    new_node[MNEW_NODE][MNEW_NODE_LIST], *el=NULL, *nodes=NULL, 
    *split_nodes=NULL, *new_node_numbers=NULL, *in_geometry_list=NULL;
  double rdum=0., tmp_node[MDIM], tmp_node_start_refined[MDIM], 
    *ddum=NULL, *tmp_node_dof=NULL, *tmp_node_dof_start_refined=NULL,
    *node_dof=NULL, *node=NULL, *node_start_refined=NULL, 
    *node_dof_start_refined=NULL;

  if ( ndim==1 ) return;

  swit = set_swit(-1,-1,"mesh_split");
  if ( swit ) pri( "In routine MESH_SPLIT" );

  db_max_index( NODE, max_node, version, GET );
  db_max_index( ELEMENT, max_element, version, GET );

  el = get_new_int(MNOL+1);
  nodes = get_new_int(MNOL);
  split_nodes = get_new_int(MSPLIT_ELEMENT*MSPLIT_NODE);
  new_node_numbers = get_new_int(MNEW_NODE);
  ddum = get_new_dbl(MDIM+MUKNWN);
  array_set( ddum, 0., MUKNWN );
  tmp_node_dof = get_new_dbl(MUKNWN);
  tmp_node_dof_start_refined = get_new_dbl(MUKNWN);
  in_geometry_list = get_new_int(1+max_node);
  array_set( in_geometry_list, 1, 1+max_node );

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  if ( db_active_index( CONTROL_MESH_SPLIT_ONLY, icontrol, VERSION_NORMAL ) ) {
    db( CONTROL_MESH_SPLIT_ONLY, icontrol, control_mesh_split_only, ddum, 
      ldum, version, GET );
    for ( inod=0; inod<=max_node; inod++ ) {
      if ( db_active_index( NODE_START_REFINED, inod, VERSION_NORMAL ) ) {
        geometry( inod, ddum, control_mesh_split_only, in_geometry, rdum, ddum, rdum,
          ddum, NODE_START_REFINED, PROJECT_EXACT, VERSION_NORMAL );
        if ( !in_geometry ) in_geometry_list[inod] = 0;
      }
    }
  }

  for ( element=0; element<=max_element; element++ ) {
    if ( db_active_index( ELEMENT, element, version ) ) {
      db( ELEMENT, element, el, ddum, length, version, GET );
      name = el[0];
      nnol = length - 1; array_move( &el[1], nodes, nnol );
      all_in_geometry = 1;
      for ( inol=0; inol<nnol; inol++ ) {
        inod = nodes[inol];
        if ( !in_geometry_list[inod] ) all_in_geometry = 0;
      }
      if ( all_in_geometry && ( name==-QUAD4 || name==-QUAD9 || name==-HEX8 || name==-HEX27 ) ) {
        nnew_node = -1;
        nnew_number = -1;
        if      ( name==-QUAD4 ) {
          split_nodes[0] = 0; split_nodes[1] = 1;
          split_nodes[2] = 1; split_nodes[3] = 3;
          split_nodes[4] = 3; split_nodes[5] = 2;
          split_nodes[6] = 2; split_nodes[7] = 0;
          nsplit_element = 4;
          nsplit_node = 2;
          nnew_node = 1;
          nnew_number = nnol;
          for ( inol=0; inol<nnol; inol++ ) new_node[0][inol] = inol;
        }
        else if ( name==-QUAD9 ) {
          nsplit_element = 4;
            /* number of new nodes */
          nnew_node = 4; 
            /* number of old nodes at the middle of which a new node is made */
          nnew_number = 2; 
          new_node[0][0] = 0;
          new_node[0][1] = 4;
          new_node[1][0] = 2;
          new_node[1][1] = 4;
          new_node[2][0] = 8;
          new_node[2][1] = 4;
          new_node[3][0] = 6;
          new_node[3][1] = 4;
        }
        else if ( name==-HEX8 ) {
          split_nodes[0]  = 0; split_nodes[1]  = 1; split_nodes[2]  = 3;
          split_nodes[3]  = 0; split_nodes[4]  = 3; split_nodes[5]  = 2;
          split_nodes[6]  = 4; split_nodes[7]  = 5; split_nodes[8]  = 7;
          split_nodes[9]  = 4; split_nodes[10] = 7; split_nodes[11] = 6;
          split_nodes[12] = 1; split_nodes[13] = 3; split_nodes[14] = 7;
          split_nodes[15] = 1; split_nodes[16] = 7; split_nodes[17] = 5;
          split_nodes[18] = 0; split_nodes[19] = 2; split_nodes[20] = 6;
          split_nodes[21] = 0; split_nodes[22] = 6; split_nodes[23] = 4;
          split_nodes[24] = 0; split_nodes[25] = 1; split_nodes[26] = 5;
          split_nodes[27] = 0; split_nodes[28] = 5; split_nodes[29] = 4;
          split_nodes[30] = 2; split_nodes[31] = 3; split_nodes[32] = 7;
          split_nodes[33] = 2; split_nodes[34] = 7; split_nodes[35] = 6;
          nsplit_element = 12;
          nsplit_node = 3;
          nnew_node = 1;
          nnew_number = nnol;
          for ( inol=0; inol<nnol; inol++ ) new_node[0][inol] = inol;
        }
        else if ( name==-HEX27 ) {
          split_nodes[0]  = 0; split_nodes[1]  = 1; split_nodes[2]  = 2;
          split_nodes[3]  = 3; split_nodes[4]  = 4; split_nodes[5]  = 6;
          split_nodes[6]  = 9; split_nodes[7]  = 10; split_nodes[8]  = 12;
          split_nodes[9]  = 18;   
          split_nodes[10] = 18; split_nodes[11] = 21; split_nodes[12] = 24; 
          split_nodes[13] = 19; split_nodes[14] = 22; split_nodes[15] = 20; 
          split_nodes[16] = 12; split_nodes[17] = 15; split_nodes[18] = 13; 
          split_nodes[19] = 6;
          split_nodes[20] = 18; split_nodes[21] = 19; split_nodes[22] = 20;
          split_nodes[23] = 10; split_nodes[24] = 11; split_nodes[25] = 2;
          split_nodes[26] = 12; split_nodes[27] = 13; split_nodes[28] = 4;
          split_nodes[29] = 6;
          split_nodes[30] = 8; split_nodes[31] = 7; split_nodes[32] = 6; 
          split_nodes[33] = 5; split_nodes[34] = 4; split_nodes[35] = 2; 
          split_nodes[36] = 17; split_nodes[37] = 16; split_nodes[38] = 14; 
          split_nodes[39] = 26; 
          split_nodes[40] = 24; split_nodes[41] = 25; split_nodes[42] = 26;
          split_nodes[43] = 22; split_nodes[44] = 23; split_nodes[45] = 20;
          split_nodes[46] = 15; split_nodes[47] = 16; split_nodes[48] = 13;
          split_nodes[49] = 6;
          split_nodes[50] = 20; split_nodes[51] = 23; split_nodes[52] = 26;
          split_nodes[53] = 11; split_nodes[54] = 14; split_nodes[55] = 2;
          split_nodes[56] = 13; split_nodes[57] = 16; split_nodes[58] = 4;
          split_nodes[59] = 6;
          nsplit_element = 6;
          nsplit_node = 10;
        }
        assert( nsplit_element<=MSPLIT_ELEMENT );
        assert( nsplit_node<=MSPLIT_NODE );
        assert( nnew_node<=MNEW_NODE );
          // generate old nodes
        for ( inol=0; inol<nnol; inol++ ) {
          inod = nodes[inol];
          node = db_dbl( NODE, inod, version );
          if ( db_active_index( NODE_START_REFINED, inod, version ) )
            node_start_refined = db_dbl( NODE_START_REFINED, inod, version );
          else
            node_start_refined = ddum;
          if ( nuknwn>0 ) {
            if ( db_active_index( NODE_DOF_START_REFINED, inod, version ) )
              node_dof_start_refined = 
              db_dbl( NODE_DOF_START_REFINED, inod, version );
            else
              node_dof_start_refined = ddum;
            node_dof = db_dbl( NODE_DOF, inod, version );
          }
          else {
            node_dof_start_refined = ddum;
            node_dof = ddum;
          }
          iwork[0] = inod; create_node( iwork, 1, inod, node, node_dof,
            node_start_refined, node_dof_start_refined,
            version, VERSION_TMP ); 
        }
          // generate new nodes
        for ( inol=0; inol<nnew_node; inol++ ) {
          array_set( tmp_node, 0., ndim );
          array_set( tmp_node_start_refined, 0., ndim );
          if ( nuknwn>0 ) {
            array_set( tmp_node_dof_start_refined, 0., nuknwn );
            array_set( tmp_node_dof, 0., nuknwn );
          }
          for ( jnol=0; jnol<nnew_number; jnol++ ) {
            inod = nodes[new_node[inol][jnol]];
            node = db_dbl( NODE, inod, version );
            if ( db_active_index( NODE_START_REFINED, inod, version ) )
              node_start_refined = db_dbl( NODE_START_REFINED, inod, version );
            else
              node_start_refined = ddum;
            array_add( node, tmp_node, tmp_node, ndim );
            array_add( node_start_refined, tmp_node_start_refined, 
              tmp_node_start_refined, ndim );
            if ( nuknwn>0 ) {
              if ( db_active_index( NODE_DOF_START_REFINED, inod, version ) )
                node_dof_start_refined = 
                db_dbl( NODE_DOF_START_REFINED, inod, version );
              else
                node_dof_start_refined = ddum;
              array_add( node_dof_start_refined, tmp_node_dof_start_refined, 
                tmp_node_dof_start_refined, nuknwn );
              node_dof = db_dbl( NODE_DOF, inod, version );
              array_add( node_dof, tmp_node_dof, tmp_node_dof, nuknwn );
            }
            else {
              node_dof_start_refined = ddum;
              node_dof = ddum;
            }
          }
          array_multiply( tmp_node, tmp_node, 1./nnew_number, ndim );
          array_multiply( tmp_node_start_refined, tmp_node_start_refined, 
            1./nnew_number, ndim );
          if ( nuknwn>0 ) {
            array_multiply( tmp_node_dof_start_refined, tmp_node_dof_start_refined, 
              1./nnew_number, nuknwn );
            array_multiply( tmp_node_dof, tmp_node_dof, 1./nnew_number, nuknwn );
          }
          max_node++;
          new_node_numbers[inol] = max_node;
          create_node( nodes, nnol, max_node, tmp_node, tmp_node_dof,
            tmp_node_start_refined, tmp_node_dof_start_refined,
            version, VERSION_TMP );
        }
        if      ( name==-QUAD4 ) {
          el[0] = -TRIA3;
          length_element = 4;
        }
        else if ( name==-QUAD9 ) {
          el[0] = -TRIA6;
          length_element = 7;
        }
        else if ( name==-HEX8 ) {
          el[0] = -TET4;
          length_element = 5;
        }
        else {
          assert( name==-HEX27 );
          el[0] = -TET10;
          length_element = 11;
        }
          // generate new elements
        for ( isplit_element=0; isplit_element<nsplit_element; isplit_element++ ) {
          if ( name==-QUAD9 ) {
            if      ( isplit_element==0 ) {
              el[1+0] = nodes[0];
              el[1+1] = nodes[1];
              el[1+2] = nodes[2];
              el[1+3] = new_node_numbers[0];
              el[1+4] = new_node_numbers[1];
              el[1+5] = nodes[4];
            }
            else if ( isplit_element==1 ) {
              el[1+0] = nodes[2];
              el[1+1] = nodes[5];
              el[1+2] = nodes[8];
              el[1+3] = new_node_numbers[1];
              el[1+4] = new_node_numbers[2];
              el[1+5] = nodes[4];
            }
            else if ( isplit_element==2 ) {
              el[1+0] = nodes[8];
              el[1+1] = nodes[7];
              el[1+2] = nodes[6];
              el[1+3] = new_node_numbers[2];
              el[1+4] = new_node_numbers[3];
              el[1+5] = nodes[4];
            }
            else {
              assert( isplit_element==3 );
              el[1+0] = nodes[6];
              el[1+1] = nodes[3];
              el[1+2] = nodes[0];
              el[1+3] = new_node_numbers[3];
              el[1+4] = new_node_numbers[0];
              el[1+5] = nodes[4];
            }
          }
          else {
            for ( isplit_node=0; isplit_node<nsplit_node; isplit_node++ )
              el[1+isplit_node] = 
            nodes[split_nodes[isplit_element*nsplit_node+isplit_node]];
            if ( name==-QUAD4 || name==-HEX8 ) 
              el[1+nsplit_node] = new_node_numbers[0];
          }
          create_element( element, tmp_max_element, el, length_element,
            version, VERSION_TMP );
          tmp_max_element++;
        }
      }
      else {
          // just recreate old element
        for ( inol=0; inol<nnol; inol++ ) {
          inod = nodes[inol];
          node = db_dbl( NODE, inod, version );
          if ( db_active_index( NODE_START_REFINED, inod, version ) )
            node_start_refined = db_dbl( NODE_START_REFINED, inod, version );
          else
            node_start_refined = ddum;
          if ( nuknwn>0 ) {
            if ( db_active_index( NODE_DOF_START_REFINED, inod, version ) )
              node_dof_start_refined = 
              db_dbl( NODE_DOF_START_REFINED, inod, version );
            else
              node_dof_start_refined = ddum;
            node_dof = db_dbl( NODE_DOF, inod, version );
          }
          iwork[0] = inod; create_node( iwork, 1, inod, node, node_dof,
            node_start_refined, node_dof_start_refined,
            version, VERSION_TMP ); 
        }
        create_element( element, tmp_max_element, el, length,
          version, VERSION_TMP );
        tmp_max_element++;
      }
    }
  }

  db_version_copy( VERSION_TMP, version );
  db_version_delete( VERSION_TMP );

  mesh_has_changed( version );

  delete[] el;
  delete[] nodes;
  delete[] split_nodes;
  delete[] new_node_numbers;
  delete[] ddum;
  delete[] tmp_node_dof;
  delete[] tmp_node_dof_start_refined;
  delete[] in_geometry_list;

  if ( swit ) pri( "Out routine MESH_SPLIT" );

}
