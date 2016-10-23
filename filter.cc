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

long int filter( long int print_filter_index[], long int length_print_filter_index,
  long int data_number, long int index, long int number, long int task )

{
  long int ifil=0, found=1, max=0, filter_length=0, ind=0,
    any=0, all=0, range_length=0, data_class=0, inol=0, nnol=0,
    inod=0, itmp=0, length=0, node_macro=0, element_macro=0,
    ldum=0, *el=NULL, *nodes=NULL, 
    *filter=NULL, *integer_range=NULL, *dof_label=NULL;
  double rdum=0., ddum[MDIM];

  db_max_index( PRINT_FILTER, max, VERSION_NORMAL, GET );
  if ( max>=0 ) {
    el = get_new_int(MNOL+1);
    nodes = get_new_int(MNOL);
    filter = get_new_int(DATA_ITEM_SIZE);
    integer_range = get_new_int(MRANGE);
    dof_label = get_new_int(MUKNWN);
    db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    for ( ifil=0; ifil<=max; ifil++ ) {
      if ( db_active_index( PRINT_FILTER, ifil, VERSION_NORMAL ) ) {
        if ( print_filter_index[0]!=-NONE ) {
          if ( print_filter_index[0]==-ALL || 
              array_member(print_filter_index,ifil,length_print_filter_index,ldum) ) {
            db( PRINT_FILTER, ifil, filter, ddum, filter_length, VERSION_NORMAL, GET );
            if ( data_number==filter[0] ) {
              if      ( filter[1]==-ALL ) 
                length = 1;
              else if ( filter[1]==-RA ) {
                range_expand( &filter[1], integer_range, length, range_length );
                if ( task==CHECK_INDEX &&
                  !array_member( integer_range, index, range_length, ldum ) ) found=0;
              }
              else if ( filter[1]==-MACRO ) {
                if ( task==CHECK_INDEX ) {
                  data_class = db_data_class( data_number );
                  if      ( data_class==NODE ) {
                    node_macro = -NONE; db( NODE_MACRO_GENERATE, index, &node_macro, ddum, ldum, 
                      VERSION_NORMAL, GET_IF_EXISTS );
                    if ( node_macro!=filter[2] ) found = 0;
                  }
                  else if ( data_class==ELEMENT ) {
                    element_macro = -NONE; db( ELEMENT_MACRO_GENERATE, index, &element_macro, ddum, ldum, 
                      VERSION_NORMAL, GET_IF_EXISTS );
                    if ( element_macro!=filter[2] ) found = 0;
                  }
                  else 
                    db_error( PRINT_FILTER, ifil );
                }
                length = 2;
              }
              else if ( filter[1]<0 && db_data_class(filter[1])==GEOMETRY ) {
                if ( task==CHECK_INDEX ) {
                  data_class = db_data_class( data_number );
                  if      ( data_class==NODE ) {
                    if ( db_active_index( NODE, index, VERSION_NORMAL ) ) {
                      geometry( index, ddum, &filter[1], found, rdum, ddum, rdum,
                        ddum, NODE_START_REFINED, PROJECT_EXACT, VERSION_NORMAL );
                    }
                    else found = 0;
                  }
                  else if ( data_class==ELEMENT ) {
                    if ( db_active_index( ELEMENT, index, VERSION_NORMAL ) ) {
                      db( ELEMENT, index, el, ddum, length, VERSION_NORMAL, GET );
                      nnol = length - 1; array_move( &el[1], nodes, nnol );
                      all = 1;
                      for ( inol=0; inol<nnol; inol++ ) {
                        inod = nodes[inol];
                        geometry( inod, ddum, &filter[1], itmp, rdum, ddum, rdum,
                          ddum, NODE_START_REFINED, PROJECT_EXACT, VERSION_NORMAL );
                        if ( !itmp ) all = 0;
                      }
                      if ( all ) found = 1;
                      else found = 0;
                    }
                    else found = 0;
                  }
                  else db_error( PRINT_FILTER, ifil );
                }
                length = 2;
              }
              else {
                if ( task==CHECK_INDEX && !(index==filter[1]) ) found = 0;
                length = 1;
              }
              if ( task==CHECK_NUMBER ) {
                for ( ind=1+length, any=0; ind<filter_length; ind++ ) {
                  if ( filter[ind]<0 ) {
                    if ( filter[ind]==-ALL ) 
                      any = 1;
                    else {
                      if      ( db_len(data_number,index,VERSION_NORMAL)==npuknwn ) {
                        if ( number<npuknwn && filter[ind]==dof_label[number*nder] ) any = 1;
                      }
                      else if ( db_len(data_number,index,VERSION_NORMAL)==nuknwn ) {
                        if ( number<nuknwn && filter[ind]==dof_label[number] ) any = 1;
                      }
                    }
                  }
                  else if ( number==filter[ind] ) any = 1;
                }
                if ( !any ) found = 0;
              }
            }
          }
        }
      }
    }
    delete[] el;
    delete[] nodes;
    delete[] filter;
    delete[] integer_range;
    delete[] dof_label;
  }

  return found;

}
