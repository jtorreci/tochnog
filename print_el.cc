/*
    Copyright (C) 2000  Dennis Roddeman
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

void print_element( long int data_item_name )

{
  long int element=0, max_element=0, ival=0, nval=0, inol=0, nnol=2,
    icontrol=0, test=0, idim=0, inod=0, ok=0, length_print_filter_index=0, 
    length=0, ldum=0, idum[1], nodes[MNOL], el[1+MNOL], 
    print_filter_index[DATA_ITEM_SIZE];
  double ddum[1], dval[DATA_ITEM_SIZE], *coord=NULL;
  char filename[MCHAR], str[MCHAR];

  set_swit(-1,-1,"print_element");

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  if ( !db( CONTROL_PRINT_FILTER, icontrol, print_filter_index, ddum, 
      length_print_filter_index, VERSION_NORMAL, GET_IF_EXISTS ) ) {
    print_filter_index[0] = -ALL; length_print_filter_index = 1;
  }

  nval = db_data_length( data_item_name );
  db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
  for ( ival=0; ival<nval; ival++ ) {
    test = filter( print_filter_index, length_print_filter_index,
      data_item_name, element, ival, CHECK_NUMBER );
    if ( test ) {
      strcpy( filename, db_name(data_item_name) );
      strcat( filename, "_" );
      long_to_a( ival, str );
      strcat( filename, str );
      strcat( filename, "." );
      long_to_a( icontrol, str );
      strcat( filename, str );
      ofstream out( filename );
      out.precision(TN_PRECISION);
      for ( element=0; element<=max_element; element++ ) {
        if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) ) {
          test = filter( print_filter_index, length_print_filter_index,
            data_item_name, element, ival, CHECK_INDEX );
          if ( test ) {
            if ( db_active_index( data_item_name, element, VERSION_NORMAL ) ) {
              db( data_item_name, element, idum, dval, ldum, VERSION_NORMAL, GET );
              db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
              nnol = length - 1;
              array_move( &el[1], nodes, nnol );
              for ( inol=0; inol<nnol; inol++ ) {
                ok = 0;
                if      ( data_item_name==-ELEMENT_TRUSS_FORCE )
                  ok = 1;
                else if ( data_item_name==-ELEMENT_BEAM_MOMENT ) {
                  if      ( inol==0 && ival<nval/2 )
                    ok = 1;
                  else if ( inol==1 && ival>=nval/2 )
                    ok = 1;
                }
                else
                  db_error( CONTROL_PRINT_ELEMENT, icontrol );
                if ( ok ) {
                  inod = nodes[inol];
                  coord = db_dbl( NODE, inod, VERSION_NORMAL );     
                  for ( idim=0; idim<ndim; idim++ ) out << coord[idim] << " ";
                  out << dval[ival] << "\n";
                }
              }
            }
          }
        }
      }
      out.close();
    }
  }

}
