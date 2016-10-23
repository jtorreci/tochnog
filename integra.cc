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

long int integration_gauss( long int niso, double iso[], double weight_iso[] )
{
  long int ok=0;

  if      ( niso==1 ) {
    ok = 1;
    iso[0] =  0.00000000000; weight_iso[0] = 2.00000000000000000000000/2.;
  }
  else if ( niso==2 ) {
    ok = 1;
    iso[0] = -0.57735026919; weight_iso[0] = 1.0000000000000000000000/2.;
    iso[1] = +0.57735026919; weight_iso[1] = 1.0000000000000000000000/2.;
  }
  else if ( niso==3 ) {
    ok = 1;
    iso[0] = -0.77459666924; weight_iso[0] = 0.55555555555555555555555/2.;
    iso[1] =  0.00000000000; weight_iso[1] = 0.88888888888888888888888/2.;
    iso[2] = +0.77459666924; weight_iso[2] = 0.55555555555555555555555/2.;
  }
  else if ( niso==4 ) {
    ok = 1;
    iso[0] = -0.86113631159; weight_iso[0] = 0.34785484514/2.;
    iso[1] = -0.33998104358; weight_iso[1] = 0.65214515486/2.;
    iso[2] = +0.33998104358; weight_iso[2] = 0.65214515486/2.;
    iso[3] = +0.86113631159; weight_iso[3] = 0.34785484514/2.;
  }
  else if ( niso==5 ) {
    ok = 1;
    iso[0] = -0.90617984594; weight_iso[0] = 0.23692688506/2.;
    iso[1] = -0.53846931011; weight_iso[1] = 0.47862867050/2.;
    iso[2] =  0.00000000000; weight_iso[2] = 0.56888888888/2.;
    iso[3] = +0.53846931011; weight_iso[3] = 0.47862867050/2.;
    iso[4] = +0.90617984594; weight_iso[4] = 0.23692688506/2.;
  }

  return ok;
}

long int integration_lobatto( long int niso, double iso[], 
  double weight_iso[])
{
  long int ok=0;

  if      ( niso==1 ) {
    ok = 1;
    iso[0] =  0.00000000000; weight_iso[0] = 2.000000000000/2.;
  }
  else if ( niso==2 ) {
    ok = 1;
    iso[0] = -1.00000000000; weight_iso[0] = 0.5000000000000000000000;
    iso[1] = +1.00000000000; weight_iso[1] = 0.5000000000000000000000;
  }
  else if ( niso==3 ) {
    ok = 1;
    iso[0] = -1.00000000000; weight_iso[0] = 0.1666666666666666666666;
    iso[1] =  0.00000000000; weight_iso[1] = 0.6666666666666666666666;
    iso[2] = +1.00000000000; weight_iso[2] = 0.1666666666666666666666;
  }
  else if ( niso==4 ) {
    ok = 1;
    iso[0] = -1.00000000000; weight_iso[0] = 0.0833333333333333333333;
    iso[1] = -0.44721359014; weight_iso[1] = 0.4166666666666666666666;
    iso[2] = +0.44721359014; weight_iso[2] = 0.4166666666666666666666;
    iso[3] = +1.00000000000; weight_iso[3] = 0.0833333333333333333333;
  }
  else {
    ok = 1;
    assert( niso==5 );
    iso[0] = -1.00000000000; weight_iso[0] = 0.050000000000000000000;
    iso[1] = -0.65465366840; weight_iso[1] = 0.272222222222222222222;
    iso[2] =  0.00000000000; weight_iso[2] = 0.355555555555555555555;
    iso[3] = +0.65465366840; weight_iso[3] = 0.272222222222222222222;
    iso[4] = +1.00000000000; weight_iso[4] = 0.050000000000000000000;
  }

  return ok;
}
