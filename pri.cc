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

void pri( const char *s )
{
  cout << s;
  cout << "\n";
  cout << flush;
}

void pri( const char *s, const char *st )
{
  cout << s << ": " << st;
  cout << "\n";
  cout << flush;
}

void pri( const char *s, int i ) // print short integer
{
  cout << s;
  cout << ": ";
  if ( i<0 )
    cout << db_name(i) << "\n";
  else
    cout << i << "\n";
  cout << flush;
}

void pri( const char *s, long int i ) // print integer
{
  cout << s;
  cout << ": ";
  if ( i<0 )
    cout << db_name(i) << "\n";
  else
    cout << i << "\n";
  cout << flush;
}

void pri( const char *s, double d )  // print double
{
  cout << s;
  cout << ": ";
  cout << d << "\n";
  cout << flush;
}

void pri( const char *s, int *i, int n )  // print short integer vector
{
  long int j=0;

  cout << s;
  cout << ":\n";
  for ( j=0; j<n; j++ )  {
    if ( *(i+j)<0 )
      cout << db_name(*(i+j)) << " ";
    else
      cout << *(i+j) << " ";
  }
  cout << "\n";
  cout << flush;
}


void pri( const char *s, long int *i, long int n )  // print integer vector
{
  long int j=0;

  cout << s;
  cout << ":\n";
  for ( j=0; j<n; j++ )  {
    if ( *(i+j)<0 ) {
      cout << db_name(*(i+j)) << " ";
    }
    else
      cout << *(i+j) << " ";
  }
  cout << "\n";
  cout << flush;
}

void pri( const char *s, double *d, long int n ) // print double vector
{
  long int j=0;

  cout << s;
  cout << ":\n";
  for ( j=0; j<n; j++ ) {
    cout << *(d+j) << " ";
  }
  cout << "\n";
  cout << flush;
}

void pri( const char *s, long int *iar, long int n, long int m )
{
  long int i=0, j=0;

  cout << s;
  cout << ":\n";
  for ( i=0; i<n; i++ ) {
    for ( j=0; j<m; j++ ) {
      if ( *(iar+i*m+j)<0 )
        cout << db_name(*(iar+i*m+j)) << " ";
      else
        cout << *(iar+i*m+j) << " ";
    }
    cout << "\n";
  }
  cout << flush;
}


void pri( const char *s, double *d, long int n, long int m )
{
  long int i=0, j=0;

  cout << s;
  cout << ":\n";
  for ( i=0; i<n; i++ ) {
    for ( j=0; j<m; j++ ) {
      cout << *(d+i*m+j) << " ";
    }
    cout << "\n";
  }
  cout << flush;
}
