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

void print_history( long int ival[], long int nval )

{
  long int data_item_name=0, data_item_index=0, iset=0, nset=0,
    number=0, len=0, icontrol=0, swit=0, ldum=0, 
    idum[1], *idat=NULL, *dof_label=NULL;
  double time_current=0., ddum[1], *ddat=NULL, factor=1., *factor_d=NULL;
  char str[MCHAR], filename[MCHAR];

  swit = set_swit(-1,-1,"print_history");
  if ( swit ) pri( "In routine PRINT_HISTORY" );

  idat = get_new_int(DATA_ITEM_SIZE);
  dof_label = get_new_int(MUKNWN);
  ddat = get_new_dbl(DATA_ITEM_SIZE);

  db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( TIME_CURRENT, 0, idum, &time_current, ldum, VERSION_NORMAL, GET );
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  nset = nval / 3;

  // multiplication factors for the printed data values
  if ( db_active_index( CONTROL_PRINT_HISTORY_FACTOR, icontrol,
      VERSION_NORMAL ) ) {
    factor_d = get_new_dbl(DATA_ITEM_SIZE);
    array_set( factor_d, 1., DATA_ITEM_SIZE );
    db( CONTROL_PRINT_HISTORY_FACTOR, icontrol, idum, factor_d,
      ldum, VERSION_NORMAL, GET );
  }

  for ( iset=0; iset<nset; iset++ ) {

    data_item_name = ival[iset*3+0];
    data_item_index = labs(ival[iset*3+1]);
    if ( ival[iset*3+2]<0 ) {
      array_member(dof_label,ival[iset*3+2],nuknwn,number);
      if ( db_active_index(data_item_name,data_item_index,VERSION_NORMAL) &&
           db_len(data_item_name,data_item_index,VERSION_NORMAL)==npuknwn ) 
        number /= nder;
    }
    else
      number = ival[iset*3+2];

    if ( ival[iset*3+2]<0 ) 
      strcpy( filename, db_name(labs(ival[iset*3+2])) );
    else 
      strcpy( filename, long_to_a(ival[iset*3+2],str) );
    strcat( filename, long_to_a(data_item_index,str) );
    strcat( filename, ".his" );

    if ( db_active_index( data_item_name, data_item_index, VERSION_NORMAL ) ) {
      ofstream out( filename, ios::app );
      out.precision(TN_PRECISION);
      factor = ( factor_d ? factor_d[iset] : 1. );
      if ( db_type(data_item_name)==INTEGER ) {
        db( data_item_name, data_item_index, idat, ddum, len, VERSION_NORMAL, GET );
        if ( number<0 || number>len-1 )
          db_error( CONTROL_PRINT_HISTORY, icontrol );
        out << time_current << " " << (long int)(factor*idat[number]) << "\n";
      }
      else {
        db( data_item_name, data_item_index, idum, ddat, len, VERSION_NORMAL, GET );
        if ( number<0 || number>len-1 )
          db_error( CONTROL_PRINT_HISTORY, icontrol );
        out << time_current << " " << factor*ddat[number] << "\n";
      }
      out.close();
    }

  }

  delete[] idat;
  delete[] dof_label;
  delete[] ddat;
  if ( factor_d ) delete[] factor_d;

  if ( swit ) pri( "Out routine PRINT_HISTORY" );
}

// print_history_smooth - control_print_history_smooth.
// Smooths the data values printed by control_print_history (same index)
// by averaging the last N values (moving average). Each smooth_N value
// applies to the corresponding data value of control_print_history; a
// single value applies to all. Results are written to a separate history
// file whose name starts with "smooth".
void print_history_smooth( long int ival[], long int nval )

{
  long int data_item_name=0, data_item_index=0, iset=0, nset=0,
    number=0, len=0, icontrol=0, swit=0, ldum=0, smooth=1, nsmooth=0,
    smooth_size=0, idum[1], *idat=NULL, *dof_label=NULL, *smooth_d=NULL,
    *history_d=NULL;
  double time_current=0., ddum[1], *ddat=NULL;
  char str[MCHAR], filename[MCHAR];

  // persistent buffers of the last values per history set
  static double* smooth_buf[DATA_ITEM_SIZE];
  static long int smooth_count[DATA_ITEM_SIZE];
  static long int smooth_n[DATA_ITEM_SIZE];
  static int smooth_buf_init=0;

  swit = set_swit(-1,-1,"print_history_smooth");
  if ( swit ) pri( "In routine PRINT_HISTORY_SMOOTH" );

  idat = get_new_int(DATA_ITEM_SIZE);
  dof_label = get_new_int(MUKNWN);
  ddat = get_new_dbl(DATA_ITEM_SIZE);

  db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( TIME_CURRENT, 0, idum, &time_current, ldum, VERSION_NORMAL, GET );
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  nset = nval / 3;

  // smooth window sizes: one per data value, or a single value for all
  if ( db_active_index( CONTROL_PRINT_HISTORY_SMOOTH, icontrol,
      VERSION_NORMAL ) ) {
    smooth_d = get_new_int(DATA_ITEM_SIZE);
    db( CONTROL_PRINT_HISTORY_SMOOTH, icontrol, smooth_d, ddum,
      ldum, VERSION_NORMAL, GET );
    nsmooth = ( ldum>0 ) ? ldum : 0;
  }
  if ( nsmooth<=0 ) {
    if ( swit ) pri( "control_print_history_smooth: no data" );
    delete[] idat; delete[] dof_label; delete[] ddat;
    if ( smooth_d ) delete[] smooth_d;
    return;
  }
  smooth_size = ( nsmooth>=nset ) ? nset : nsmooth;

  if ( !smooth_buf_init ) {
    for ( long int k=0; k<DATA_ITEM_SIZE; k++ ) {
      smooth_buf[k] = NULL;
      smooth_count[k] = 0;
      smooth_n[k] = 1;
    }
    smooth_buf_init = 1;
  }

  for ( iset=0; iset<nset; iset++ ) {

    // window for this set: explicit value or the single value for all
    smooth = ( nsmooth>=nset ) ? smooth_d[iset] : smooth_d[0];
    if ( smooth<1 ) smooth = 1;
    if ( smooth>DATA_ITEM_SIZE ) smooth = DATA_ITEM_SIZE;

    data_item_name = ival[iset*3+0];
    data_item_index = labs(ival[iset*3+1]);
    if ( ival[iset*3+2]<0 ) {
      array_member(dof_label,ival[iset*3+2],nuknwn,number);
      if ( db_active_index(data_item_name,data_item_index,VERSION_NORMAL) &&
           db_len(data_item_name,data_item_index,VERSION_NORMAL)==npuknwn ) 
        number /= nder;
    }
    else
      number = ival[iset*3+2];

    if ( ival[iset*3+2]<0 ) 
      strcpy( filename, "smooth" );
    else 
      strcpy( filename, "smooth" );
    if ( ival[iset*3+2]<0 ) 
      strcat( filename, db_name(labs(ival[iset*3+2])) );
    else 
      strcat( filename, long_to_a(ival[iset*3+2],str) );
    strcat( filename, long_to_a(data_item_index,str) );
    strcat( filename, ".his" );

    if ( !db_active_index( data_item_name, data_item_index, VERSION_NORMAL ) )
      continue;

    // buffer for this set (allocate/reallocate on window change)
    if ( smooth_buf[iset]==NULL || smooth_n[iset]!=smooth ) {
      if ( smooth_buf[iset] ) delete[] smooth_buf[iset];
      smooth_buf[iset] = get_new_dbl(smooth);
      smooth_n[iset] = smooth;
      smooth_count[iset] = 0;
    }

    double value = 0.;
    if ( db_type(data_item_name)==INTEGER ) {
      db( data_item_name, data_item_index, idat, ddum, len, VERSION_NORMAL, GET );
      if ( number<0 || number>len-1 )
        db_error( CONTROL_PRINT_HISTORY_SMOOTH, icontrol );
      value = (double)idat[number];
    }
    else {
      db( data_item_name, data_item_index, idum, ddat, len, VERSION_NORMAL, GET );
      if ( number<0 || number>len-1 )
        db_error( CONTROL_PRINT_HISTORY_SMOOTH, icontrol );
      value = ddat[number];
    }

    // ring buffer
    long int idx = smooth_count[iset] % smooth;
    smooth_buf[iset][idx] = value;
    smooth_count[iset]++;
    long int nbuf = ( smooth_count[iset]<smooth ) ? smooth_count[iset] : smooth;
    double sum=0.;
    for ( long int k=0; k<nbuf; k++ ) sum += smooth_buf[iset][k];
    double avg = sum / ((double) nbuf);

    ofstream out( filename, ios::app );
    out.precision(TN_PRECISION);
    out << time_current << " " << avg << "\n";
    out.close();
  }

  delete[] idat;
  delete[] dof_label;
  delete[] ddat;
  if ( smooth_d ) delete[] smooth_d;

  if ( swit ) pri( "Out routine PRINT_HISTORY_SMOOTH" );
}
