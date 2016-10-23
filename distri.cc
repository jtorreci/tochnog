/*
    Copyright (C) 1998  Dennis Roddeman
    email: dennis.roddeman@feat.nl

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is element_edge in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.


    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software Foundation 
    59 Temple Place, Suite 330, Boston, MA, 02111-1307, USA
*/

#include "tochnog.h"

void distribute( void )

{
  int distribute_ran=-1;
  long int distribution_type=0, data_item_name=0, data_item_number=0, length=0,
    icontrol=0, index=0, max_index=0, group_data=0, name=0, 
    length_control_distribute=0, 
    length_control_distribute_values=0,
    length_element_distribute=0, 
    length_element_distribute_values=0,
    ldum=0, idistribute=0, ndistribute=0,
    idum[1], element_distribute[2],
    *control_distribute=NULL, *dof_label=NULL;
  double range=0., variance=0., ran=0., tmp=0., ddum[1], 
    *element_distribute_values=NULL, *control_distribute_values=NULL, 
    *dval=NULL;
  char str[MCHAR];

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  if ( db_active_index( CONTROL_DISTRIBUTE, icontrol, VERSION_NORMAL ) ) {
    control_distribute = get_new_int(DATA_ITEM_SIZE);
    dof_label = get_new_int(MUKNWN);
    element_distribute_values = get_new_dbl(DATA_ITEM_SIZE);
    control_distribute_values = get_new_dbl(DATA_ITEM_SIZE);
    db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( CONTROL_DISTRIBUTE_VALUES, icontrol, idum, control_distribute_values, 
      length_control_distribute_values, VERSION_NORMAL, GET );
    db( CONTROL_DISTRIBUTE, icontrol, control_distribute, ddum, 
      length_control_distribute, VERSION_NORMAL, GET );
    if ( length_control_distribute!=3*length_control_distribute_values ) {
      pri( "Error: lengths of CONTROL_DISTRIBUTE and CONTROL_DISTRIBUTE_VALUES do not match." );
      exit(TN_EXIT_STATUS);
    }
    db_delete( ELEMENT_DISTRIBUTE, VERSION_NORMAL );
    db_delete( ELEMENT_DISTRIBUTE_VALUES, VERSION_NORMAL );
    ndistribute = length_control_distribute_values;
    for ( idistribute=0; idistribute<ndistribute; idistribute++ ) {
      distribution_type = control_distribute[idistribute*3+0];
      data_item_name = control_distribute[idistribute*3+1];
      data_item_number = control_distribute[idistribute*3+2];
      strcpy( str, db_name( data_item_name ) );
      group_data = ( str[0]=='g' && str[1]=='r' && str[2]=='o' && str[3]=='u' && 
        str[4]=='p' );
      if ( group_data )
        name = -ELEMENT;
      else
        name = data_item_name;
      if      ( distribution_type==-UNIFORM )
        range = control_distribute_values[idistribute];
      else if ( distribution_type==-NORMAL )
        variance = control_distribute_values[idistribute];
      else
        db_error(  CONTROL_DISTRIBUTE, icontrol );
      distribute_ran = -1;
      db_max_index( name, max_index, VERSION_NORMAL, GET );
      for ( index=0; index<=max_index; index++ ) {
        if ( db_active_index( name, index, VERSION_NORMAL ) ) {
          if ( distribution_type==-UNIFORM ) {
            ran = scalar_ran_uniform( distribute_ran );
            tmp = range*ran - 0.5*range;
          }
          else if ( distribution_type==-NORMAL ) {
            ran = scalar_ran_normal( distribute_ran );
            tmp = variance*ran;
          }
          else
            db_error(  CONTROL_DISTRIBUTE, icontrol );
          if ( name==-ELEMENT ) {
            if ( db_active_index( ELEMENT_DISTRIBUTE, index, VERSION_NORMAL ) ) {
              db( ELEMENT_DISTRIBUTE, index, element_distribute, ddum, 
                length_element_distribute, VERSION_NORMAL, GET );
              db( ELEMENT_DISTRIBUTE_VALUES, index, idum, element_distribute_values,
                length_element_distribute_values, VERSION_NORMAL, GET );
            }
            else {
              length_element_distribute = 0;
              length_element_distribute_values = 0;
            }
            element_distribute[length_element_distribute*2+0] = data_item_name;
            element_distribute[length_element_distribute*2+1] = data_item_number;
            element_distribute_values[length_element_distribute_values] = tmp;
            length_element_distribute += 2;
            length_element_distribute_values += 1;
            db( ELEMENT_DISTRIBUTE, index, element_distribute, ddum, 
              length_element_distribute, VERSION_NORMAL, PUT );
            db( ELEMENT_DISTRIBUTE_VALUES, index, idum, element_distribute_values, 
              length_element_distribute_values, VERSION_NORMAL, PUT );
          }
          else {
            if ( data_item_name==-NODE ) {
              if ( db_active_index( NODE_BOUNDARY, index, VERSION_NORMAL ) ) tmp = 0.;
            }
            dval = db_dbl( data_item_name, index, VERSION_NORMAL );
            length = db_len( data_item_name, index, VERSION_NORMAL );
            data_item_number = control_distribute[idistribute*3+2];
            if ( data_item_number<0 ) {
              array_member(dof_label,data_item_number,nuknwn,data_item_number);
              if ( db_len(data_item_name,index,VERSION_NORMAL)==npuknwn )
                data_item_number /= nder;
            }
            if ( data_item_number<0 || data_item_number>length-1 ) 
              db_error(  CONTROL_DISTRIBUTE, icontrol );
            dval[data_item_number] += tmp;
          }
        }
      }
    }
    delete[] control_distribute;
    delete[] dof_label;
    delete[] element_distribute_values;
    delete[] control_distribute_values;
  }
}
