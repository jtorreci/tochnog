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

// print_dof - control_print_dof: print the primary dofs with the
// coordinates at which they hold. Lines like "x y z dof" per node; in 1D
// only x, etc. The coordinates themselves are also printed in separate
// files. Filenames: dof.<index> (-separate_index) or dof.<n>
// (-separate_sequential).
void print_dof( long int icontrol, long int task )

{
  long int inod=0, idim=0, ipuknwn=0, iuknwn=0, nder_=0, nuknwn_=0,
    swit=0, ldum=0, nval=0, seq=0;
  long int idum[1], *dof_label=NULL, *dof_scal_vec_mat=NULL;
  double ddum[1], coord[MDIM], *node_dof=NULL;
  char filename[MCHAR], str[MCHAR];

  swit = set_swit(-1,-1,"print_dof");
  if ( swit ) pri( "In routine PRINT_DOF" );

  db_version_copy( VERSION_NORMAL, VERSION_PRINT );
  renumbering( VERSION_PRINT, NO, 0, 0, idum, idum );
  db_highest_index( NODE, inod, VERSION_PRINT );
  long int max_node = inod;
  if ( max_node<0 ) return;
  nder_ = nder;
  nuknwn_ = nuknwn;

  dof_label = get_new_int(MUKNWN);
  dof_scal_vec_mat = get_new_int(MUKNWN);
  db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( DOF_SCAL_VEC_MAT, 0, dof_scal_vec_mat, ddum, ldum, VERSION_NORMAL,
    GET_IF_EXISTS );

  // file name: dof.<index> or dof.<seq>
  strcpy( filename, "dof." );
  if      ( task==-SEPARATE_INDEX && icontrol>=0 ) {
    long_to_a( icontrol, str );
    strcat( filename, str );
  }
  else if ( task==-SEPARATE_SEQUENTIAL ) {
    static long int dof_seq=0;
    long_to_a( dof_seq++, str );
    strcat( filename, str );
  }
  else {
    long_to_a( icontrol, str );
    strcat( filename, str );
  }

  ofstream out( filename, ios::app );
  out.precision(TN_PRECISION);

  // coordinates file: coord.<index> (written only the first time)
  {
    char cfname[MCHAR];
    strcpy( cfname, "coord." );
    if      ( task==-SEPARATE_INDEX && icontrol>=0 ) {
      long_to_a( icontrol, str );
      strcat( cfname, str );
    }
    else if ( task==-SEPARATE_SEQUENTIAL ) {
      long_to_a( seq, str );
      strcat( cfname, str );
    }
    else {
      long_to_a( icontrol, str );
      strcat( cfname, str );
    }
    {
      std::ifstream fexists( cfname );
      if ( !fexists.is_open() ) {
        ofstream outcoord( cfname, ios::app );
        outcoord.precision(TN_PRECISION);
        for ( inod=0; inod<=max_node; inod++ ) {
          db( NODE, inod, idum, coord, ldum, VERSION_PRINT, GET );
          for ( idim=0; idim<ndim; idim++ )
            outcoord << coord[idim] << " ";
          outcoord << "\n";
        }
        outcoord.close();
      }
      fexists.close();
    }
  }

  for ( ipuknwn=0; ipuknwn<nuknwn_; ipuknwn++ ) {
    if ( dof_scal_vec_mat[ipuknwn]!=-SCALAR &&
         dof_scal_vec_mat[ipuknwn]!=-VECTOR &&
         dof_scal_vec_mat[ipuknwn]!=-MATRIX ) continue;
    iuknwn = ipuknwn*nder_;
    nval = 1;
    if      ( dof_scal_vec_mat[ipuknwn]==-VECTOR ) nval = ndim;
    else if ( dof_scal_vec_mat[ipuknwn]==-MATRIX ) nval = 6;
    for ( long int k=0; k<nval; k++ ) {
      for ( inod=0; inod<=max_node; inod++ ) {
        db( NODE, inod, idum, coord, ldum, VERSION_PRINT, GET );
        node_dof = db_dbl( NODE_DOF, inod, VERSION_PRINT );
        long int indx = iuknwn;
        if      ( dof_scal_vec_mat[ipuknwn]==-VECTOR )
          indx = iuknwn + k*nder_;
        else if ( dof_scal_vec_mat[ipuknwn]==-MATRIX ) {
          long int kk, ll;
          if      ( k==0 ) { kk=0; ll=0; }
          else if ( k==1 ) { kk=1; ll=1; }
          else if ( k==2 ) { kk=2; ll=2; }
          else if ( k==3 ) { kk=0; ll=1; }
          else if ( k==4 ) { kk=0; ll=2; }
          else              { kk=1; ll=2; }
          indx = iuknwn + stress_indx(kk,ll)*nder_;
        }
        for ( idim=0; idim<ndim; idim++ )
          out << coord[idim] << " ";
        out << node_dof[indx] << "\n";
      }
    }
  }

  out.close();

  db_version_delete( VERSION_PRINT );
  delete[] dof_label;
  delete[] dof_scal_vec_mat;

  if ( swit ) pri( "Out routine PRINT_DOF" );
}
