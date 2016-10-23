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

#define INITIAL_SIZE 100

void nod_nod( long int version )

{
  long int i=0, inol=0, inod=0, jnol=0, jnod=0, nnol=0, indx=0,
    ixi=0, jxi=0, nxi=0, ieta=0, jeta=0, neta=0, izeta=0, jzeta=0, nzeta=0,
    name=0, element=0, max_element=0, max_node=0, length=0, old_length=0,
    max_length_node_node=0, max_length_node_element=0, swit=0, ldum=0,
    allready_present=0, direct_neighbours[MNOL][MNOL],
    el[1+MNOL], nodes[MNOL], *active_node=NULL, 
    **node_element=NULL, *length_node_element=NULL, *nelement=NULL, 
    **node_node=NULL, *length_node_node=NULL, *nnode=NULL, *int_ptr=NULL;
  double ddum[1];

  swit = set_swit(-1,-1,"nod_nod");
  if ( swit ) pri( "In routine NOD_NOD" );

  db_max_index( NODE, max_node, version, GET );
  db_max_index( ELEMENT, max_element, version, GET );

  length = 1+max_node;
  active_node = get_new_int( length );
  nelement = get_new_int( length );
  nnode = get_new_int( length );
  length_node_element = get_new_int( length );
  length_node_node = get_new_int( length );
  node_element = new long int * [length];
  node_node = new long int * [length];
  for ( inod=0; inod<=max_node; inod++ ) {
    active_node[inod] = 0;
    nelement[inod] = 0;
    nnode[inod] = 0;
    length_node_element[inod] = INITIAL_SIZE;
    length_node_node[inod] = INITIAL_SIZE;
    node_element[inod] = get_new_int(INITIAL_SIZE);
    node_node[inod] = get_new_int(INITIAL_SIZE);
    array_set( node_element[inod], 0, INITIAL_SIZE );
    array_set( node_node[inod], 0, INITIAL_SIZE );
  }                            

  for ( element=0; element<=max_element; element++ ) {
    if ( db_active_index( ELEMENT, element, version ) ) {
      db( ELEMENT, element, el, ddum, length, version, GET );
      name = el[0];
      nnol = length - 1; array_move( &el[1], nodes, nnol );
      array_set( &direct_neighbours[0][0], 0, MNOL*MNOL );
      if      ( name==-BAR2 || name==-BAR3 || name==-BAR4 ||
                name==-QUAD4 || name==-QUAD9 || name==-QUAD16 ||
                name==-HEX8 || name==-HEX27 || name==-HEX64 ) {
        nxi = neta = nzeta = 1;
        if      ( name==-BAR2 ) {
          nxi = 2;
        }
        else if ( name==-BAR3 ) {
          nxi = 3;
        }
        else if ( name==-BAR4 ) {
          nxi = 4;
        }
        else if ( name==-QUAD4 ) {
          nxi = 2;
          neta = 2;
        }
        else if ( name==-QUAD9 ) {
          nxi = 3;
          neta = 3;
        }
        else if ( name==-QUAD16 ) {
          nxi = 4;
          neta = 4;
        }
        else if ( name==-HEX8 ) {
          nxi = 2;
          neta = 2;
          nzeta = 2;
        }
        else if ( name==-HEX27 ) {
          nxi = 3;
          neta = 3;
          nzeta = 3;
        }
        else if ( name==-HEX64 ) {
          nxi = 4;
          neta = 4;
          nzeta = 4;
        }
        inol = 0;
        for ( izeta=0; izeta<nzeta; izeta++ ) {
          for ( ieta=0; ieta<neta; ieta++ ) {
            for ( ixi=0; ixi<nxi; ixi++ ) {
              jnol = 0;
              for ( jzeta=0; jzeta<nzeta; jzeta++ ) {
                for ( jeta=0; jeta<neta; jeta++ ) {
                  for ( jxi=0; jxi<nxi; jxi++ ) {
                    if ( labs(jxi-ixi)>1 ||  
                         labs(jeta-ieta)>1 ||
                         labs(jzeta-izeta)>1 ) {
                      direct_neighbours[inol][jnol] = 
                        direct_neighbours[jnol][inol] = 0;
                    }
                    else {
                      direct_neighbours[inol][jnol] = 
                        direct_neighbours[jnol][inol] = 1;
                    }
                    jnol++;
                  }
                }
              }
              inol++;
            }
          }
        }
      }
      else if ( name==-TRIA6 ) {
        array_set( &direct_neighbours[0][0], 0, MNOL*MNOL );
        direct_neighbours[0][1] = 1;
        direct_neighbours[1][0] = 1;
        direct_neighbours[0][3] = 1;
        direct_neighbours[3][0] = 1;
        direct_neighbours[0][4] = 1;
        direct_neighbours[4][0] = 1;
        direct_neighbours[2][1] = 1;
        direct_neighbours[1][2] = 1;
        direct_neighbours[2][4] = 1;
        direct_neighbours[4][2] = 1;
        direct_neighbours[2][5] = 1;
        direct_neighbours[5][2] = 1;
        direct_neighbours[6][3] = 1;
        direct_neighbours[3][6] = 1;
        direct_neighbours[6][4] = 1;
        direct_neighbours[4][6] = 1;
        direct_neighbours[6][5] = 1;
        direct_neighbours[5][6] = 1;
      }
      else {
        array_set( &direct_neighbours[0][0], 1, MNOL*MNOL );
      }
      for ( inol=0; inol<nnol; inol++ ) {
        inod = nodes[inol];
        if ( !db_active_index( NODE, inod, version ) ) {
          pri( "Error detected for element ", element );
          pri( "Non existing node ", inod );
          exit_tn_on_error();
        }
        active_node[inod] = 1;
        if ( !array_member(node_element[inod],element,nelement[inod],ldum) ) {
          if ( nelement[inod]>length_node_element[inod]-1 ) {
            int_ptr = get_new_int( length_node_element[inod] );
            for ( i=0; i<length_node_element[inod]; i++ ) {
              int_ptr[i] = node_element[inod][i];
            }
            delete[] node_element[inod];
            old_length = length_node_element[inod];
            length_node_element[inod] += INITIAL_SIZE;
            node_element[inod] = get_new_int( length_node_element[inod] );
            array_set( node_element[inod], 0, length_node_element[inod] );
            for ( i=0; i<old_length; i++ ) {
              node_element[inod][i] = int_ptr[i];
            }
            delete[] int_ptr;
          }
          indx = nelement[inod];
          node_element[inod][indx] = element;
          nelement[inod]++;
          if ( nelement[inod]>max_length_node_element ) 
            max_length_node_element = nelement[inod];
        }
        for ( jnol=0; jnol<nnol; jnol++ ) {
          jnod = nodes[jnol];
          allready_present = 0;
          for ( i=0; i<nnode[inod]; i++ ) {
            if ( jnod==labs(node_node[inod][i]) ) allready_present = 1;
          }
          if ( jnod!=inod && !allready_present ) {
            if ( nnode[inod]>length_node_node[inod]-1 ) {
              int_ptr = get_new_int( length_node_node[inod] );
              for ( i=0; i<length_node_node[inod]; i++ ) {
                int_ptr[i] = node_node[inod][i];
              }
              delete[] node_node[inod];
              old_length = length_node_node[inod];
              length_node_node[inod] += INITIAL_SIZE;
              node_node[inod] = get_new_int( length_node_node[inod] );
              array_set( node_node[inod], 0, length_node_node[inod] );
              for ( i=0; i<old_length; i++ ) {
                node_node[inod][i] = int_ptr[i];
              }
              delete[] int_ptr;
            }
            indx = nnode[inod];
            if ( direct_neighbours[inol][jnol] ) {
              node_node[inod][indx] = +jnod;
            }
            else {
              node_node[inod][indx] = -jnod;
            }
            nnode[inod]++;
            if ( nnode[inod]>max_length_node_node ) 
              max_length_node_node = nnode[inod];
          }
        }
      }
    }
  }
    // this deletes all versions!!
  db_data_length_put( NODE_ELEMENT, max_length_node_element );
  db_data_length_put( NODE_NODE, max_length_node_node );
  for ( inod=0; inod<=max_node; inod++ ) {
    if ( db_active_index( NODE, inod, version ) ) {
      if ( !active_node[inod] ) 
        delete_node( inod, version );
      else {
        length = nnode[inod];
        if ( length>0 ) db( NODE_NODE, inod, node_node[inod], ddum, 
          nnode[inod], version, PUT );
        length = nelement[inod];
        if ( length>0 ) db( NODE_ELEMENT, inod, node_element[inod], ddum, 
          nelement[inod], version, PUT );
      }
    }
  }
  delete[] active_node;
  delete[] nelement;
  delete[] nnode;
  delete[] length_node_element;
  delete[] length_node_node;
  for ( inod=0; inod<=max_node; inod++ ) {
    delete[] node_element[inod];
    delete[] node_node[inod];
  }                            
  delete[] node_element;
  delete[] node_node;

     // since all versions are deleted we restore VERSION_NORMAL again!!
  if ( version!=VERSION_NORMAL ) nod_nod( VERSION_NORMAL );

  if ( swit ) pri( "Out routine NOD_NOD" );
}
