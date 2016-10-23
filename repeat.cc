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

#define EPS 1.e-10

long int repeat( long int &start_control )

{
  long int icontrol=0, length=1, repeat_return=0, repeat_criterium=0,
    ldum=0, data_item_name=0, data_item_index=0, data_item_number=0, number=0, 
    idum[1], control_repeat[2], control_repeat_until_item[5],
    *dof_label=NULL;
  double control_repeat_until_tolerance=0., control_repeat_until_value=0.,
    old_control_repeat_until_value=0., ddum[1], *val=NULL;

  dof_label = get_new_int(MUKNWN);
  val = get_new_dbl(DATA_ITEM_SIZE);
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  if      ( db_active_index( CONTROL_REPEAT, icontrol, VERSION_NORMAL ) ) {
    db( CONTROL_REPEAT, icontrol, control_repeat, ddum, ldum, VERSION_NORMAL, GET );
    start_control = control_repeat[1];
    if ( icontrol<start_control ) db_error( CONTROL_REPEAT, icontrol );
    if ( control_repeat[0]>0 ) {
      repeat_return = 1;
      control_repeat[0]--; 
      length=2; db( CONTROL_REPEAT, icontrol, control_repeat, ddum, length, VERSION_NORMAL, PUT );
    }
  }
  else if ( db_active_index( CONTROL_REPEAT_UNTIL_ITEM, icontrol, VERSION_NORMAL ) ) {
    db( CONTROL_REPEAT_UNTIL_ITEM, icontrol, control_repeat_until_item, ddum,
      ldum, VERSION_NORMAL, GET );
    db( CONTROL_REPEAT_UNTIL_TOLERANCE, icontrol, idum, &control_repeat_until_tolerance, 
      ldum, VERSION_NORMAL, GET );
    start_control = control_repeat_until_item[0];
    repeat_criterium = control_repeat_until_item[1];
    data_item_name = control_repeat_until_item[2];
    data_item_index = control_repeat_until_item[3];
    data_item_number = control_repeat_until_item[4];
    if ( db_active_index( data_item_name, data_item_index, VERSION_NORMAL ) ) {
      db( data_item_name, data_item_index, idum, val, length, VERSION_NORMAL, GET );
      if ( data_item_number<0 ) {
        db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
        array_member(dof_label,data_item_number,nuknwn,number);
        if ( db_len(data_item_name,data_item_index,VERSION_NORMAL)==npuknwn ) 
          number /= nder;
      }
      else
        number = data_item_number;
      if ( number>=0 && number<=length ) {
        control_repeat_until_value = val[number];
      }
      else
        db_error( CONTROL_REPEAT_UNTIL_ITEM, icontrol );
      if ( repeat_criterium==-CHANGE ) {
        if ( db_active_index( CONTROL_REPEAT_UNTIL_VALUE, icontrol, VERSION_NORMAL ) ) {
           db( CONTROL_REPEAT_UNTIL_VALUE, icontrol, idum, &old_control_repeat_until_value, 
             ldum, VERSION_NORMAL, GET );
          if ( old_control_repeat_until_value>EPS &&
               ( (control_repeat_until_value-old_control_repeat_until_value)/
                 old_control_repeat_until_value < control_repeat_until_tolerance/100. ) ) {
            repeat_return = 0;
            db_delete( CONTROL_REPEAT_UNTIL_VALUE, VERSION_NORMAL );
          }
          else {
            repeat_return = 1;
            db( CONTROL_REPEAT_UNTIL_VALUE, icontrol, idum, &control_repeat_until_value, 
              ldum, VERSION_NORMAL, PUT );
          }
        }
        else {
          repeat_return = 1;
          db( CONTROL_REPEAT_UNTIL_VALUE, icontrol, idum, &control_repeat_until_value, 
            ldum, VERSION_NORMAL, PUT );
        }
      }
      else if ( repeat_criterium==-VALUE ) {
        if ( scalar_dabs(control_repeat_until_value) < control_repeat_until_tolerance )
          repeat_return = 0;
        else
          repeat_return = 1;
      }
      else
        db_error( CONTROL_REPEAT_UNTIL_ITEM, icontrol );
    }
    else
      repeat_return = 1;
  }

  delete[] dof_label;
  delete[] val;
  return repeat_return;
}
