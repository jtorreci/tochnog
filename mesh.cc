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

void mesh_has_changed( long int version )

{ 
  long int idat=0, swit=0;

  swit = set_swit(-1,-1,"mesh_has_changed");
  if ( swit ) pri( "In routine MESH_HAS_CHANGED" );

  db_delete( CONTROL_PRINT_TECPLOT_MESH, VERSION_NORMAL );
  db_delete( NODE_ELEMENT, version );
  db_delete( NODE_NODE, version );

  if ( version==VERSION_NORMAL ) {
    for ( idat=0; idat<MDAT; idat++ ) {
      if ( db_data_class(idat)==NODE && !db_version( idat, VERSION_TMP ) )
        db_delete( idat, version );
    }
  }
  area_element_group( version );
  tendon_distribute();
  area_node_dataitem();
  nod_nod(version);
  nonlocal_first_set=0;
 
  if ( swit ) pri( "Out routine MESH_HAS_CHANGED" );
}

void mesh_add( long int version_from, long int version_to )

{

  long int inol=0, inod=0, indx=0, idat=0, element=0, max_node=0, max_element=0, 
    length=0, nnol=0, data_class=0, swit=0, idum[1], *nodes=NULL, *el=NULL, 
    *ival=NULL, *new_nodes=NULL;
  double ddum[1], *dval=NULL;

  swit = set_swit(-1,-1,"mesh_add");
  if ( swit ) pri( "In routine MESH_ADD" );

  db_max_index( NODE, max_node, version_from, GET );
  db_max_index( ELEMENT, max_element, version_from, GET );
  nodes = get_new_int( MNOL );
  el = get_new_int( 1+MNOL );
  new_nodes = get_new_int( 1+max_node );
  for ( inod=0; inod<=max_node; inod++ ) {
    if ( db_active_index( NODE, inod, version_from ) ) {
      indx = inod;
      while ( db_active_index( NODE, indx, version_to ) ) indx++;
      new_nodes[inod] = indx;
      for ( idat=0; idat<MDAT; idat++ ) {
        data_class = db_data_class( idat );
        if ( data_class==NODE && 
          db_version( idat, version_from ) &&
          db_version( idat, version_to ) ) {
          if ( db_active_index( idat, inod, version_from ) ) {
            length = db_len( idat, inod, version_from );
            if ( db_type(idat)==DOUBLE_PRECISION ) {
              dval = db_dbl( idat, inod, version_from );
              db( idat, indx, idum, dval, length, version_to, PUT );
            }
            else {
              ival = db_int( idat, inod, version_from );
              db( idat, indx, ival, ddum, length, version_to, PUT );
            }
          }
        }
      }
    }
  }
  for ( element=0; element<=max_element; element++ ) {
    if ( db_active_index( ELEMENT, element, version_from ) ) {
      indx = element;
      while ( db_active_index( ELEMENT, indx, version_to ) ) indx++;
      for ( idat=0; idat<MDAT; idat++ ) {
        data_class = db_data_class( idat );
        if ( data_class==ELEMENT && 
          db_version( idat, version_from ) &&
          db_version( idat, version_to ) ) {
          if ( db_active_index( idat, element, version_from ) ) {
            length = db_len( idat, element, version_from );
            if ( db_type(idat)==DOUBLE_PRECISION ) {
              dval = db_dbl( idat, element, version_from );
              db( idat, indx, idum, dval, length, version_to, PUT );
            }
            else {
              ival = db_int( idat, element, version_from );
              db( idat, indx, ival, ddum, length, version_to, PUT );
            }
          }
        }
      }
      db( ELEMENT, element, el, ddum, length, version_from, GET );
      nnol = length - 1; array_move( &el[1], nodes, nnol );
      for ( inol=0; inol<nnol; inol++ ) {
        inod = nodes[inol];
        el[1+inol] = new_nodes[inod];
      }
      db( ELEMENT, indx, el, ddum, length, version_to, PUT );
    }
  }
  delete[] nodes;
  delete[] el;
  delete[] new_nodes;
  db_version_delete( version_from );
  mesh_has_changed( version_to );

  if ( swit ) pri( "Out routine MESH_ADD" );
}
