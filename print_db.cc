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

void print_database( long int icontrol, long int version, long int data_number )

{
  long int i=0, l=0, idat=0, total_size=0, printed_any_index=0, 
    length_print_filter_index=0, max=0, index=0, end_of_data=0, ready=0, 
    test=0, ninitia=0, size=0, print_filter_index[DATA_ITEM_SIZE],
    initialization_values[DATA_ITEM_SIZE], *ival=NULL;
  double ddum[1], *dval=NULL;
  char str[MCHAR], filename[MCHAR];

  set_swit(-1,-1,"print_database");

  db( INITIALIZATION_VALUES, 0, initialization_values, ddum, 
    ninitia, VERSION_NORMAL, GET_IF_EXISTS );

  if ( data_number==-YES || data_number==-EVERYTHING || data_number==-SIZETOT ) {
    strcpy( filename, data_file_base );
    if ( icontrol!=-1 ) {
      long_to_a( icontrol, str );
      strcat( filename, str );
    }
    strcat( filename, ".dbs" );
  }
  else if ( data_number==-PRINT_LASTDATABASE ) {
    strcpy( filename, "last.dbs" );
  }
  else {
    strcpy( filename, "tn.tmp" );
  }
  ofstream out( filename );
  out.precision(TN_PRECISION);

  print_filter_index[0] = -ALL; length_print_filter_index = 1;
  if ( icontrol!=-1 ) {
    if ( !db( CONTROL_PRINT_FILTER, icontrol, print_filter_index, ddum, 
        length_print_filter_index, VERSION_NORMAL, GET_IF_EXISTS ) ) {
    }
  }

  if ( data_number==-YES || data_number==-EVERYTHING || data_number==-PRINT_LASTDATABASE ) {
    for ( i=0; i<ninitia; i++ ) {
      out << initialization_names[i];
      if      ( initialization_values[i]==-YES )
        out << " -yes";
      else if ( initialization_values[i]==-NO )
        out << " -no";
      else if ( initialization_values[i]!=-EMPTY )
        out << " " << initialization_values[i];
      out << "\n";
    }
    out << "\n";
  }


  while ( !ready ) {
    test = ( idat==labs(data_number) );
    test = test || ( data_number==-YES && db_external(idat) );
    test = test || ( data_number==-PRINT_LASTDATABASE && db_external(idat) );
    test = test || ( data_number==-EVERYTHING );
    test = test || ( data_number==-SIZETOT );
    if ( test ) {
      if ( db_version(idat,version) ) {
        db_max_index( idat, max, version, GET );
        if ( data_number==-SIZETOT ) {
          if ( max>=0 ) {
            size = max * db_data_length( idat );
            if ( db_type(idat)==INTEGER )
              size *= sizeof(long);
            else {
              assert( db_type(idat)==DOUBLE_PRECISION );
              size *= sizeof(double);
            }
            out << "Size of " <<  db_name(idat) << " is " << size << "\n";
            total_size += size;
          }
        }
        else {
          index = printed_any_index = 0;
          do {
            test = ( filter( print_filter_index, length_print_filter_index,
              data_number, index, i, CHECK_INDEX ) );
            test = test || ( data_number==-YES );
            test = test || ( data_number==-PRINT_LASTDATABASE );
            test = test || ( data_number==-EVERYTHING );
            test = test && ( db_active_index(idat,index,version) );
            if ( test ) {
              if ( data_number==-YES || data_number==-EVERYTHING ||
                  data_number==-PRINT_LASTDATABASE ) {
                out << db_name(idat) << "  ";
                if ( !db_no_index(idat) ) out << index << " ";
              }
              else {
                cout << db_name(idat) << "  ";
                if ( !db_no_index(idat) ) cout << index << " ";
              }
              printed_any_index = 1;
              if ( db_type(idat)==INTEGER )
                ival = db_int( idat, index, version );
              else
                dval = db_dbl( idat, index, version );
              end_of_data = 0;
              for ( i=0; !end_of_data; i++ ) {
                test = ( filter( print_filter_index, length_print_filter_index,
                  data_number, index, i, CHECK_NUMBER ) );
                test = test || ( data_number==-YES );
                test = test || ( data_number==-PRINT_LASTDATABASE );
                test = test || ( data_number==-EVERYTHING );
                if ( test ) {
                  if ( db_type(idat)==INTEGER ) {
                    if ( ival[i]<0 ) {
                      if ( data_number==-YES || data_number==-EVERYTHING ||
                          data_number==-PRINT_LASTDATABASE )
                        out  << " " << "-" << db_name(ival[i]);
                      else
                        cout << " " << "-" << db_name(ival[i]);
                    }
                    else {
                      if ( data_number==-YES || data_number==-EVERYTHING ||
                          data_number==-PRINT_LASTDATABASE )
                        out  << " " << ival[i];
                      else
                        cout << " " << ival[i];
                    }
                  }
                  else if ( db_type(idat)==DOUBLE_PRECISION ) {
                    if ( data_number==-YES || data_number==-EVERYTHING ||
                        data_number==-PRINT_LASTDATABASE )
                      out  << " " << dval[i];
                    else
                      cout << " " << dval[i];
                  }
                }
                l = db_len( idat, index, version );
                if ( (i+1)==l ) end_of_data = 1;
              }
              if ( data_number==-YES || data_number==-EVERYTHING ||
                  data_number==-PRINT_LASTDATABASE )
                out  << " \n";
              else
                cout << " \n";
            }
            index++;
          }
          while ( index<=max );
          if ( printed_any_index ) {
            if ( data_number==-YES || data_number==-EVERYTHING ||
                data_number==-PRINT_LASTDATABASE )
              out  << " \n";
            else
              cout << " \n";
          }
        }
      }
    }
    idat++;
    ready = ( idat==MDAT-1 );
  }

  if ( data_number==-SIZETOT ) 
    out << "\n\nTotal size is " << total_size << "\n";

  if ( data_number==-YES || data_number==-EVERYTHING ||
      data_number==-PRINT_LASTDATABASE ) {
    out << "end_data\n";
    out.close();
  }
  else cout << flush;

}
