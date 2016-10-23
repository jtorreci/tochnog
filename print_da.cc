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

void print_data_versus_data( long int ival[], long int length )

{
  long int data_item_name=0, data_item_index=0, number=0,
    len=0, icontrol=0, swit=0, idata=0, ndata=0, ldum=0, idum[1], 
    *idat=NULL, *dof_label=NULL;
  double ddum[1], *ddat=NULL;
  char filename[MCHAR];

  swit = set_swit(-1,-1,"print_data_versus_data");
  if ( swit ) pri( "In routine PRINT_DATA_VERSUS_DATA" );

  idat = get_new_int(DATA_ITEM_SIZE);
  dof_label = get_new_int(MUKNWN);
  ddat = get_new_dbl(DATA_ITEM_SIZE);

  db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );

  if ( length%3!=0 ) db_error( CONTROL_PRINT_DATA_VERSUS_DATA, icontrol );
  ndata = length / 3;

  strcpy( filename, "tn.dvd" );
  ofstream out( filename, ios::app );
  out.precision(TN_PRECISION);

  for ( idata=0; idata<ndata; idata++ ) {

    data_item_name = ival[idata*3+0];
    data_item_index = ival[idata*3+1];
    if ( ival[idata*3+2]<0 ) {
      array_member(dof_label,ival[idata*3+2],nuknwn,number);
      if ( db_len(data_item_name,data_item_index,VERSION_NORMAL)==npuknwn ) 
        number /= nder;
    }
    else 
      number = ival[idata*3+2];

    if ( db_active_index( data_item_name, data_item_index, VERSION_NORMAL ) ) {
      if ( db_type(data_item_name)==INTEGER ) {
        db( data_item_name, data_item_index, idat, ddum, len, VERSION_NORMAL, GET );
        if ( number<0 || number>len-1 )
          db_error( CONTROL_PRINT_DATA_VERSUS_DATA, icontrol );
        out << idat[number] << "  ";
      }
      else {
        db( data_item_name, data_item_index, idum, ddat, len, VERSION_NORMAL, GET );
        if ( number<0 || number>len-1 )
          db_error( CONTROL_PRINT_DATA_VERSUS_DATA, icontrol );
        out << ddat[number] << "  ";
      }
    }

  }

  out << "\n";

  out.close();

  delete[] idat;
  delete[] dof_label;
  delete[] ddat;

  if ( swit ) pri( "Out routine PRINT_DATA_VERSUS_DATA" );
}
