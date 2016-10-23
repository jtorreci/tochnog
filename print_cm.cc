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

void print_cmd( void )

  /* 
     Print in tngid.cmd all possible data commands.
     Print in tngid.ind 0/1 for data indices.
     The files will be used by the graphical layer tngid.
  */

{
  long int idat=0;

    // write data part commands
  ofstream out_cmd( "tngid.cmd" );
  for ( idat=0; idat<MDAT; idat++ ) {
    if ( db_external( idat ) && !db_print_only(idat) && check( idat, CHECK_USAGE ) &&
         ( db_type(idat)==INTEGER || db_type(idat)==DOUBLE_PRECISION ) ) {
      out_cmd << db_name(idat) << "\n";
    }
  }
  out_cmd.close();

    // write data part commands index info
  ofstream out_ind( "tngid.ind" );
  for ( idat=0; idat<MDAT; idat++ ) {
    if ( db_external( idat ) && !db_print_only(idat) && check( idat, CHECK_USAGE ) &&
         ( db_type(idat)==INTEGER || db_type(idat)==DOUBLE_PRECISION ) ) {
       if ( db_no_index(idat) )
          out_ind << "0\n";
       else
          out_ind << "1\n";
    }
  }
  out_ind.close();

}
