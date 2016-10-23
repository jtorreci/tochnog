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

#define EPS_MAT 1.e-10

long int project_point_exactly_on_line( double coord[], double coord0[], 
  double coord1[], double weight[] )

    /* Project point on the line dictated by coord0 and coord1.  */

{
  long int return_value=0;

  if ( project_point_on_line( coord, coord0, coord1, weight ) ) {
    return_value = 1;
    if ( weight[0]<0. ) weight[0] = 0.;
    if ( weight[0]>1. ) weight[0] = 1.;
    weight[1] = 1. - weight[0];
  }
  else
    return_value = 0;

  return return_value;
}

long int project_point_exactly_on_quad( double coord[], double coord0[], 
  double coord1[], double coord2[], double coord3[], double weight[] )

{
  long int idim=0, return_value=1;
  double co0[MDIM], co1[MDIM], co2[MDIM], co3[MDIM], wa[2], wb[2];

  for ( idim=0; idim<ndim; idim++ ) {
    co0[idim] = 0.5*(coord0[idim] + coord1[idim]);
    co1[idim] = 0.5*(coord1[idim] + coord2[idim]);
    co2[idim] = 0.5*(coord2[idim] + coord3[idim]);
    co3[idim] = 0.5*(coord3[idim] + coord0[idim]);
  }

  project_point_exactly_on_line( coord, co3, co1, wa );
  project_point_exactly_on_line( coord, co0, co2, wb );

  weight[0] = wa[0] * wb[0];
  weight[1] = wa[1] * wb[0];
  weight[2] = wa[0] * wb[1];
  weight[3] = wa[1] * wb[1];

  return return_value;
}

long int project_point_exactly_on_triangle( double coord[], double coord0[], 
  double coord1[], double coord2[], double weight[] )

{
  long int return_value=0;
  double work[MDIM];

  if ( project_point_on_triangle( coord, coord0, coord1, coord2, work ) ) {
      // projection on an edge of the triangle if required
    return_value = 1;
    if      ( work[0]<0. ) {
      project_point_on_line( coord, coord1, coord2, work );
      if ( work[0]<0. ) work[0] = 0.;
      if ( work[0]>1. ) work[0] = 1.;
      weight[0] = 0.;
      weight[1] = work[0];
      weight[2] = 1. - weight[1];
    }
    else if ( work[1]<0. ) {
      project_point_on_line( coord, coord0, coord2, work );
      if ( work[0]<0. ) work[0] = 0.;
      if ( work[0]>1. ) work[0] = 1.;
      weight[0] = work[0];
      weight[1] = 0.;
      weight[2] = 1. - weight[0];
    }
    else if ( work[2]<0. ) {
      project_point_on_line( coord, coord0, coord1, work );
      if ( work[0]<0. ) work[0] = 0.;
      if ( work[0]>1. ) work[0] = 1.;
      weight[0] = work[0];
      weight[1] = 1. - weight[0];
      weight[2] = 0.;
    }
    else
      array_move( work, weight, ndim );
  }
  else {
    return_value = 0;
    array_move( work, weight, ndim );
  }

  return return_value;
}

long int project_point_on_line( double coord[], double coord0[], double coord1[],
  double weight[] )

    /* Project point on the line dictated by coord0 and coord1.
       The projection is not necessarily inside the line segment. */

{
  double tmp1=0., tmp2=0., vec1[MDIM], vec2[MDIM];

  array_subtract( coord1, coord0, vec1, ndim );
  array_subtract( coord, coord0, vec2, ndim );
  tmp1 = array_inproduct( vec1, vec1, ndim );
  if ( tmp1!=0. ) {
    tmp2 = array_inproduct( vec1, vec2, ndim ) / tmp1;
    weight[0] = 1. - tmp2; 
    weight[1] = 1. - weight[0];
    return 1;
  }
  else
    return 0;
}

long int project_point_on_triangle( double coord[], double coord0[], 
  double coord1[], double coord2[], double weight[] )

    /* Project point on the plane dictated by coord0, coord1 and coord2.
       The projection is not necessarily inside the triangle. */

{
  double rdum=0., vec01[MDIM], vec02[MDIM], vec12[MDIM],
    vec1[MDIM], mat[4], inv_mat[4], rh[2], work[2];

  array_subtract(coord0, coord1, vec01, ndim);
  array_subtract(coord0, coord2, vec02, ndim);
  array_subtract(coord1, coord2, vec12, ndim);
  mat[0] = array_inproduct( vec02, vec01, ndim );
  mat[1] = array_inproduct( vec12, vec01, ndim );
  mat[2] = array_inproduct( vec02, vec02, ndim );
  mat[3] = array_inproduct( vec12, vec02, ndim );
  array_subtract( coord, coord2, vec1, ndim );
  rh[0] = array_inproduct( vec1, vec01, ndim );
  rh[1] = array_inproduct( vec1, vec02, ndim );
  if ( matrix_inverse(mat,inv_mat,rdum,2) ) {
    matrix_ab( inv_mat, rh, work, 2, 2, 1 );
    weight[0] = work[0];
    weight[1] = work[1];
    weight[2] = 1. - weight[0] - weight[1];
    return 1;
  }
  else
    return 0;

}
