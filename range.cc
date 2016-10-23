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

void range_scan( long int echo, double d, long int d_is_set, long int ival[], long int &length )

{
  long int i=0, from=0, to=0, step=0;
  char str[MCHAR];

  ival[0] = -RA;
  do {

    input_read_string( echo, str, d, d_is_set );
    string_convert_to_lower_case( str );
    if ( echo ) cout << str << " ";

    if      ( !strcmp(str,"-from") ) {
      i++; range_test_scan( i );
      ival[i] = -FROM;
      i++; range_test_scan( i );
      if ( !(cin >> ival[i]) ) {
        pri( "Error: in data part." );
        exit(TN_EXIT_STATUS);
      }
      if ( echo ) cout << ival[i] << " ";
      from = ival[i];
      if ( !(cin >> str) ) {
        pri( "Error: in data part." );
        exit(TN_EXIT_STATUS);
      }
      string_convert_to_lower_case( str );
      if ( echo ) cout << str << " ";
      if ( strcmp(str,"-to") ) {
        pri( "Error: in data part." );
        exit(TN_EXIT_STATUS);
      }
      i++; range_test_scan( i );
      ival[i] = -TO;
      i++; range_test_scan( i );
      if ( !(cin >> ival[i]) ) {
        pri( "Error: in data part." );
        exit(TN_EXIT_STATUS);
      }
      if ( echo ) cout << ival[i] << " ";
      to = ival[i];
      if ( from>=to ) {
        pri( "Error in data file." );
        exit(TN_EXIT_STATUS);
      }
      if ( !(cin >> str) ) {
        pri( "Error: in data part." );
        exit(TN_EXIT_STATUS);
      }
      string_convert_to_lower_case( str );
      if ( echo ) cout << str << " ";
      if ( !strcmp(str,"-step") ) {
        i++; range_test_scan( i );
        ival[i] = STEP;
        i++; range_test_scan( i );
        if ( !(cin >> ival[i]) ) {
          pri( "Error: in data part." );
          exit(TN_EXIT_STATUS);
        }
        if ( echo ) cout << ival[i] << " ";
        step = ival[i];
        if ( step<=0 ) {
          pri( "Error in data file." );
          exit(TN_EXIT_STATUS);
        }
      }
      else if ( !strcmp(str,"-ra") ) {
        i++; range_test_scan( i );
        ival[i] = -RA;
        length = i+1;
      }
      else {
        if ( !string_isinteger( str ) ) {
          pri( "Error: in data part." );
          exit(TN_EXIT_STATUS);
        }
        i++; range_test_scan( i );
        ival[i] = atoi( str );
      }
    }
    else if ( !strcmp(str,"-ra") ) {
      i++; range_test_scan( i );
      ival[i] = -RA;
      length = i+1;
    }
    else {
      if ( !string_isinteger( str ) ) {
        pri( "Error: in data part." );
        exit(TN_EXIT_STATUS);
      }
      i++; range_test_scan( i );
      ival[i] = atoi(str);
    }
  }
  while ( ival[i]!=-RA );

}

void range_expand( long int ival[], long int integer_range[],
  long int &length, long int &range_length )

{

  long int j=1, i=-1, k=0, from=0, to=0, step=0;

  do {
    if      ( ival[j]==-FROM ) {
      j++;
      from = ival[j];
      j++; j++;
      to = ival[j];
      j++;
      if ( ival[j]==-STEP ) {
        j++; step = ival[j];
        for ( k=from; k<=to; k+=step ) {
          i++; range_test_expand( i );
          integer_range[i] = k;
        }
        j++;
      }
      else {
        for ( k=from; k<=to; k++ ) {
          i++; range_test_expand( i );
          integer_range[i] = k;
        }
      }
    }
    else {
      i++; range_test_expand( i );
      integer_range[i] = ival[j];
      j++;
    }
  }
  while ( ival[j]!=-RA );

  length = j+1;
  range_length = i+1;

}

void range_test_scan( long int i )
 
{
  if ( i>MRANGE-1 ) {
    cout << "\nError: maximum MRANGE in tochnog.h exceeded.\n";
    cout << "\nIncrease it.\n";
    exit(TN_EXIT_STATUS);
   }
}            

void range_test_expand( long int i )

{
  if ( i>MRANGE-1 ) {
    cout << "\nError: maximum range length MRANGE in tochnog.h exceeded.\n";
    cout << "\nIncrease it.\n";
    exit(TN_EXIT_STATUS);
   }
}
