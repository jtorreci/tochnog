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

// Ranges "-ra ... -ra" (Professional user manual "Ranges"): a list of
// integers, -from a -to b blocks (a<=b) and an optional -step s (s>0) after
// a -to, in any order and number, closed by -ra. "-ra a b -ra" is a LIST of
// the two numbers a and b, never the range a..b.

static void range_scan_error( const char *msg, const char *word )

{
  pri( "\nError in data part." );
  if ( word==NULL ) pri( msg );
  else pri( msg, word );
  exit(TN_EXIT_STATUS);
}

  // Next word of an index range through the normal reader (bounded word
  // length, include files, defines, comments; BACKLOG C56 / F4 T6f).
static void range_read_word( long int echo, char str[], double &d,
  long int &d_is_set )

{
  input_read_string( echo, str, d, d_is_set );
  input_skip_comment( str );
  string_convert_to_lower_case( str );
  if ( echo ) cout << str << " ";
}

static long int range_read_integer( long int echo, char str[], double &d,
  long int &d_is_set )

{
  range_read_word( echo, str, d, d_is_set );
  if ( !string_isinteger( str ) )
    range_scan_error( "Error in range reading", str );
  return atoi( str );
}

void range_scan( long int echo, double d, long int d_is_set, long int ival[], long int &length )

{
  long int i=0, from=0, to=0, step=0, ready=0;
  char str[MCHAR];

  ival[0] = -RA;
  range_read_word( echo, str, d, d_is_set );
  while ( !ready ) {
    if      ( !strcmp(str,"-from") ) {
      from = range_read_integer( echo, str, d, d_is_set );
      range_read_word( echo, str, d, d_is_set );
      if ( strcmp(str,"-to") )
        range_scan_error( "Error in range, a -from should be followed by a -to", str );
      to = range_read_integer( echo, str, d, d_is_set );
      if ( from>to )
        range_scan_error( "Error in range, a -from is larger than a -to", NULL );
      i++; range_test_scan( i ); ival[i] = -FROM;
      i++; range_test_scan( i ); ival[i] = from;
      i++; range_test_scan( i ); ival[i] = -TO;
      i++; range_test_scan( i ); ival[i] = to;
      range_read_word( echo, str, d, d_is_set );
      if ( !strcmp(str,"-step") ) {
        step = range_read_integer( echo, str, d, d_is_set );
        if ( step<=0 )
          range_scan_error( "Error in range, a -step should be larger than 0", str );
        i++; range_test_scan( i ); ival[i] = -STEP;
        i++; range_test_scan( i ); ival[i] = step;
        range_read_word( echo, str, d, d_is_set );
      }
    }
    else if ( !strcmp(str,"-ra") ) {
      if ( i==0 ) range_scan_error( "Error in range, the range is empty.", NULL );
      i++; range_test_scan( i );
      ival[i] = -RA;
      length = i+1;
      ready = 1;
    }
    else {
      if ( !string_isinteger( str ) )
        range_scan_error( "Error in range reading", str );
      i++; range_test_scan( i );
      ival[i] = atoi(str);
      range_read_word( echo, str, d, d_is_set );
    }
  }

}

void range_expand( long int ival[], long int integer_range[],
  long int &length, long int &range_length )

{

  long int j=1, i=-1, k=0, from=0, to=0, step=0;

  while ( ival[j]!=-RA ) {
    if      ( ival[j]==-FROM ) {
        // index ranges are checked by range_scan; value ranges come
        // from the generic record reader and are checked here
      from = ival[j+1];
      if ( ival[j+1]<0 || ival[j+2]!=-TO || ival[j+3]<0 ) {
        pri( "\nError in range, a -from should be followed by a number, -to and a number." );
        exit(TN_EXIT_STATUS);
      }
      to = ival[j+3];
      if ( from>to ) {
        pri( "\nError in range, a -from is larger than a -to" );
        exit(TN_EXIT_STATUS);
      }
      j += 4;
      step = 1;
      if ( ival[j]==-STEP ) {
        step = ival[j+1];
        if ( step<=0 ) {
          pri( "\nError in range, a -step should be a number larger than 0." );
          exit(TN_EXIT_STATUS);
        }
        j += 2;
      }
      for ( k=from; k<=to; k+=step ) {
        i++; range_test_expand( i );
        integer_range[i] = k;
      }
    }
    else {
      i++; range_test_expand( i );
      integer_range[i] = ival[j];
      j++;
    }
  }

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
