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

long int intersect_line_with_line( double line_a0[], double line_a1[], 
  double line_b0[], double line_b1[], double &iso_line_a,
  double &iso_line_b )

{
  double rdum=0., vec[MDIM], mat[MDIM*MDIM], 
    inv_mat[MDIM*MDIM], rh[MDIM], iso[MDIM];

  assert( ndim==2 );

  array_subtract( line_a1, line_a0, vec, ndim );
  mat[0*ndim+0] = vec[0];
  mat[1*ndim+0] = vec[1];
  array_subtract( line_b1, line_b0, vec, ndim );
  mat[0*ndim+1] = -vec[0];
  mat[1*ndim+1] = -vec[1];
  if ( matrix_inverse( mat, inv_mat, rdum, ndim ) ) {
    array_subtract( line_b0, line_a0, rh, ndim );
    matrix_ab( inv_mat, rh, iso, ndim, ndim, 1 ); 
    iso_line_a = iso[0]; 
    iso_line_b = iso[1];
    return 1;
  }
  else
    return 0;
}

long int intersect_line_with_point( double line0[], double line1[], 
  double point[], double &iso_line )

{
  assert( ndim==1 );

  iso_line = (point[0]-line0[0]) / (line1[0]-line0[0]);
  return 1;
}

long int intersect_line_with_triangle( double line0[], double line1[], 
  double triangle0[], double triangle1[], double triangle2[], 
  double &iso_line, double iso_triangle[] )

{
  double rdum=0., vec[MDIM], mat[MDIM*MDIM], 
    inv_mat[MDIM*MDIM], rh[MDIM], iso[MDIM];

  assert( ndim==3 );

  array_subtract( line1, line0, vec, ndim );
  mat[0*ndim+0] = vec[0];
  mat[1*ndim+0] = vec[1];
  mat[2*ndim+0] = vec[2];
  array_subtract( triangle0, triangle2, vec, ndim );
  mat[0*ndim+1] = -vec[0];
  mat[1*ndim+1] = -vec[1];
  mat[2*ndim+1] = -vec[2];
  array_subtract( triangle1, triangle2, vec, ndim );
  mat[0*ndim+2] = -vec[0];
  mat[1*ndim+2] = -vec[1];
  mat[2*ndim+2] = -vec[2];
  if ( matrix_inverse( mat, inv_mat, rdum, ndim ) ) {
    array_subtract( triangle2, line0, rh, ndim );
    matrix_ab( inv_mat, rh, iso, ndim, ndim, 1 ); 
    iso_line = iso[0]; 
    iso_triangle[0] = iso[1]; 
    iso_triangle[1] = iso[2]; 
    return 1;
  }
  else
    return 0;
}
