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

void print_restart( long int icontrol )

{
  long int idat=0, swit=0, max=0, indx=0, icontrol_current=0, ldum=0;
  double ddum[1];

  swit = set_swit(-1,-1,"print_restart");
  if ( swit ) pri( "In routine PRINT_RESTART" );

  db( ICONTROL, 0, &icontrol_current, ddum, ldum, VERSION_NORMAL, GET );

  for ( idat=0; idat<MDAT; idat++ ) {
    if ( db_data_class(idat)==CONTROL ) {
      db_max_index( idat, max, VERSION_NORMAL, GET );
      for ( indx=0; indx<=max && indx<=icontrol_current; indx++ ) {
        db_delete_index( idat, indx, VERSION_NORMAL );
      }
    }
  }

  print_database( icontrol, VERSION_NORMAL, -YES );

  if ( swit ) pri( "Out routine PRINT_RESTART" );
}
