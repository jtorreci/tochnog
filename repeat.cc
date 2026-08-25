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

// control_repeat_save (manual Professional 6.348): capture the current
// value of each data item selected by the CONTROL_REPEAT_SAVE record at
// icontrol and store it in REPEAT_SAVE_RESULT. The layout of the record
// is a list of triplets (data_item_name, data_item_index, data_item_number),
// one per item to save. Each repeat writes one index of REPEAT_SAVE_RESULT
// (first repeat -> index 0, second -> index 1, ...) holding the values of
// all selected items in triplet order. Called from repeat() when a jump
// back to the repeat start is about to happen.
static void control_repeat_save_data( long int icontrol )

{
  long int ldum=0, length=0, ndata=0, isave=0, max_save=0, idata=0,
    data_item_name=0, data_item_index=0, data_item_number=0, number=0,
    idum[1], *dof_label=NULL, *control_repeat_save=NULL;
  double ddum[1], *val=NULL, *repeat_save_values=NULL;

  length = db_len( CONTROL_REPEAT_SAVE, icontrol, VERSION_NORMAL );
  ndata = length/3;
  if ( ndata<1 ) db_error( CONTROL_REPEAT_SAVE, icontrol );
  control_repeat_save = get_new_int(length);
  db( CONTROL_REPEAT_SAVE, icontrol, control_repeat_save, ddum, ldum,
    VERSION_NORMAL, GET );
  repeat_save_values = get_new_dbl(ndata);
  dof_label = get_new_int(MUKNWN);
  // next free index of REPEAT_SAVE_RESULT: the number of ACTIVE indices
  // (db_max_index returns the allocated maximum, which includes the
  // heuristic extra of db_allocate).
  isave = 0;
  db_max_index( REPEAT_SAVE_RESULT, max_save, VERSION_NORMAL, GET );
  for ( long int i=0; i<=max_save; i++ )
    if ( db_active_index( REPEAT_SAVE_RESULT, i, VERSION_NORMAL ) ) isave++;
  for ( idata=0; idata<ndata; idata++ ) {
    data_item_name   = control_repeat_save[3*idata];
    data_item_index  = control_repeat_save[3*idata+1];
    data_item_number = control_repeat_save[3*idata+2];
    if ( db_active_index( data_item_name, data_item_index, VERSION_NORMAL ) ) {
      length = db_len( data_item_name, data_item_index, VERSION_NORMAL );
      val = get_new_dbl(length);
      db( data_item_name, data_item_index, idum, val, length,
        VERSION_NORMAL, GET );
      if ( data_item_number<0 ) {
        db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
        array_member(dof_label,data_item_number,nuknwn,number);
        if ( db_len(data_item_name,data_item_index,VERSION_NORMAL)==npuknwn ) 
          number /= nder;
      }
      else
        number = data_item_number;
      if ( number>=0 && number<length ) {
        repeat_save_values[idata] = val[number];
      }
      else
        db_error( CONTROL_REPEAT_SAVE, icontrol );
      delete[] val;
    }
      else
        db_error( CONTROL_REPEAT_SAVE, icontrol );
  }
  db( REPEAT_SAVE_RESULT, isave, idum, repeat_save_values, ndata,
    VERSION_NORMAL, PUT );
  delete[] control_repeat_save;
  delete[] repeat_save_values;
  delete[] dof_label;
}

// control_repeat_save_calculate (manual Professional 6.349): statistical
// analysis of the data saved by control_repeat_save once the repeat has
// completed (the CONTROL_REPEAT counter reached 0). For every data item
// the average and the (population) variance over all repeats are computed
// and stored in REPEAT_CALCULATE_RESULT: index = data item number, the
// record holds [average, variance].
static void control_repeat_save_calculate( long int icontrol )

{
  long int length=0, ndata=0, nsave=0, isave=0, idata=0,
    idum[1];
  double *saved=NULL, *repeat_calculate_values=NULL;

  length = db_len( CONTROL_REPEAT_SAVE, icontrol, VERSION_NORMAL );
  ndata = length/3;
  if ( ndata<1 ) return;
  // number of ACTIVE indices of REPEAT_SAVE_RESULT (db_max_index
  // returns the allocated maximum, which includes the heuristic extra
  // of db_allocate).
  nsave = 0;
  {
    long int max_save=0;
    db_max_index( REPEAT_SAVE_RESULT, max_save, VERSION_NORMAL, GET );
    for ( isave=0; isave<=max_save; isave++ )
      if ( db_active_index( REPEAT_SAVE_RESULT, isave, VERSION_NORMAL ) ) nsave++;
  }
  if ( nsave<1 ) return;
  saved = get_new_dbl(nsave*ndata);
  for ( isave=0; isave<nsave; isave++ ) {
    length = 1;
    db( REPEAT_SAVE_RESULT, isave, idum, &saved[isave*ndata], length,
      VERSION_NORMAL, GET );
  }
  repeat_calculate_values = get_new_dbl(2);
  for ( idata=0; idata<ndata; idata++ ) {
    double average=0., variance=0., diff=0.;
    for ( isave=0; isave<nsave; isave++ )
      average += saved[isave*ndata+idata];
    average /= nsave;
    for ( isave=0; isave<nsave; isave++ ) {
      diff = saved[isave*ndata+idata] - average;
      variance += diff*diff;
    }
    variance /= nsave;
    repeat_calculate_values[0] = average;
    repeat_calculate_values[1] = variance;
    length = 2;
    db( REPEAT_CALCULATE_RESULT, idata, idum, repeat_calculate_values,
      length, VERSION_NORMAL, PUT );
  }
  delete[] saved;
  delete[] repeat_calculate_values;
}

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
      // control_repeat_save (manual Professional 6.348): save the current
      // values of the selected data items before jumping back to the
      // repeat start (each repeat writes the next index of
      // REPEAT_SAVE_RESULT).
      if ( db_active_index( CONTROL_REPEAT_SAVE, icontrol, VERSION_NORMAL ) )
        control_repeat_save_data( icontrol );
      length=2; db( CONTROL_REPEAT, icontrol, control_repeat, ddum, length, VERSION_NORMAL, PUT );
    }
    // the repeat completed (counter reached 0): control_repeat_save_calculate
    // (manual Professional 6.349) performs the statistical analysis of the
    // saved data and stores average + variance in REPEAT_CALCULATE_RESULT.
    else if ( db_active_index( CONTROL_REPEAT_SAVE_CALCULATE, icontrol,
        VERSION_NORMAL ) )
      control_repeat_save_calculate( icontrol );
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
