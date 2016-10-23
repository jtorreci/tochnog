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
// This file was modified and corrected so that the contact with an ellipse geometry
// will work. Modification made by Fernando Lorenzo on June 22, 2012.

#include "tochnog.h"

void geometry( long int inod, double co[], long int geometry_entity[],
  long int &in_geometry, double &factor, double normal[], 
  double &penetration, double projection[],
  long int node_type, long int projection_type, long int version )

{
  int i=0, level=0;
  long j=0, itest=0, idim=0, index=0, ind=0, length_geometry_bounda_factor=0,
    entity=0, itriangle=0, ntriangle=0, length=0, iset=0, nset=0,
    ok=0, ldum=0, idum[1], geometry_set[DATA_ITEM_SIZE];
  double xi=0., tolerance=0., tmp=0., tmp1=0., a=0., b=0.,
    xyint[2], ellilength=0.,
    l0=0., l1=0., l2=0., x0=0., x1=0., y0=0., y1=0., l=0.,
    radius=0., cylinder_length=0., x=0., y=0., z=0., eps_iso=EPS_ISO,
    side[MDIM], geometry_bounda_sine_x[2], geometry_bounda_sine_y[2],
    geometry_bounda_sine_z[2], coord[MDIM], dydx[MDIM], 
    ddum[1], weight[MNOL], tmp_vec[MDIM], point_first[MDIM], point_second[MDIM],
    tmp_vec0[MDIM], tmp_vec1[MDIM], tmp_vec2[MDIM], tmp_vec11[MDIM],
    tmp_vec3[MDIM], centre[MDIM], vec01[MDIM], vec02[MDIM], coord0[MDIM],
    coord1[MDIM], coord2[MDIM], geometry_point[1*MDIM+1], 
    geometry_line[2*MDIM+1], geometry_triangle[2*(3*MDIM+1)], 
    geometry_circle[MDIM+2], geometry_circle_segment[MDIM+1+MDIM+1], 
    geometry_circle_smallsegment[MDIM+1+2*MDIM+1], geometry_quadrilateral[4*MDIM+1], 
    geometry_sphere[MDIM+2], geometry_sphere_segment[MDIM+1+MDIM+1], 
    geometry_cylinder[2*MDIM+2], geometry_cylinder_segment[MDIM+MDIM+1+MDIM+1], 
    geometry_ellipse[MDIM+3], geometry_bounda_factor[4], 
    geometry_brick[2*MDIM+1], work[MDIM], 
    *node_dof=NULL, *geometry_polynomial=NULL;

  factor = 1.;
  in_geometry = 0;
  array_set( normal, 0., MDIM );
  array_set( vec01, 0., MDIM );
  array_set( vec02, 0., MDIM );
  entity = geometry_entity[0];
  index  = geometry_entity[1];

  if ( inod>=0 ) {
    if      ( node_type==PLUS_DISPLACEMENT ) {
      db( NODE, inod, idum, coord, ldum, version, GET );
      if ( materi_displacement ) {
        node_dof = db_dbl( NODE_DOF, inod, version );
        for ( idim=0; idim<ndim; idim++ )
          coord[idim] += node_dof[dis_indx+idim*nder];
      }
    }
    else
      db( node_type, inod, idum, coord, ldum, version, GET );
  }
  else
    array_move( co, coord, ndim );

  if ( entity==-GEOMETRY_SET ) {
    db( GEOMETRY_SET, index, geometry_set, ddum, nset, VERSION_NORMAL, GET );
    nset = nset / 2;
  }
  else {
    geometry_set[0] = entity;
    geometry_set[1] = index;
    nset = 1;
  }

  for ( iset=0; iset<nset && !in_geometry; iset++ ) {
    entity = geometry_set[iset*2];
    index  = geometry_set[iset*2+1];
    if ( entity==-GEOMETRY_BRICK ) {
      db( GEOMETRY_BRICK, index, idum, geometry_brick, 
        ldum, VERSION_NORMAL, GET );
      tolerance = geometry_brick[2*MDIM];
      ok = 1;
      for ( idim=0; idim<MDIM; idim++ ) {
        if ( coord[idim] < geometry_brick[idim]-
             geometry_brick[MDIM+idim]/2.-tolerance ) 
          ok = 0;
        if ( coord[idim] > geometry_brick[idim]+
             geometry_brick[MDIM+idim]/2.+tolerance )
          ok = 0;
      }
      if ( ok ) in_geometry = 1;
      for ( idim=0; idim<MDIM; idim++ ) {
        if ( coord[idim] < geometry_brick[idim] )
          tmp = geometry_brick[idim]-geometry_brick[MDIM+idim];
        if ( coord[idim] > geometry_brick[idim] )
          tmp = geometry_brick[idim]+geometry_brick[MDIM+idim];
        projection[idim] = tmp;
      }
    }
    if ( entity==-GEOMETRY_CIRCLE ) {
      db( GEOMETRY_CIRCLE, index, idum, geometry_circle, 
        ldum, VERSION_NORMAL, GET );
      array_move( geometry_circle, centre, ndim );
      radius = geometry_circle[ndim];
      tolerance = geometry_circle[ndim+1];
      array_subtract( coord, centre, tmp_vec1, ndim );
      tmp = array_size( tmp_vec1, ndim );
      if ( projection_type==CONTROL_MESH_DELETE_GEOMETRY ) {
        if ( tmp<=(radius+tolerance+EPS_COORD) ) in_geometry = 1;
      }
      else {
        if ( tmp>=(radius-tolerance-EPS_COORD) && 
            tmp<=(radius+tolerance+EPS_COORD) ) in_geometry = 1;
      }
      array_move( tmp_vec1, normal, ndim );
      array_normalize( normal, ndim );
      penetration = tmp - radius;
      for ( idim=0; idim<ndim; idim++ ) {
        if      ( tmp<EPS_COORD )
           projection[idim] = centre[idim];
        else
           projection[idim] = centre[idim] + 
             (radius/tmp)*tmp_vec1[idim];
      }
    }
    else if ( entity==-GEOMETRY_CIRCLE_SEGMENT ) {
      db( GEOMETRY_CIRCLE_SEGMENT, index, idum, geometry_circle_segment, 
        ldum, VERSION_NORMAL, GET );
      array_move( geometry_circle_segment, centre, ndim );
      radius = geometry_circle_segment[ndim];
      side[0] = geometry_circle_segment[ndim+1];
      side[1] = geometry_circle_segment[ndim+2];
      tolerance = geometry_circle_segment[ndim+3];
      array_subtract( coord, centre, tmp_vec1, ndim );
      tmp = array_size( tmp_vec1, ndim );
      ok = 0;
      if ( projection_type==CONTROL_MESH_DELETE_GEOMETRY ) {
        if ( tmp<=(radius+tolerance+EPS_COORD) )
          ok = 1;
      }
      else {
        if ( tmp>=(radius-tolerance-EPS_COORD) && 
             tmp<=(radius+tolerance+EPS_COORD) )
          ok = 1;
      }
      for ( idim=0; idim<ndim; idim++ ) {
        if ( side[idim]>0. && tmp_vec1[idim]<0. ) ok = 0;
        if ( side[idim]<0. && tmp_vec1[idim]>0. ) ok = 0;
      }         
      if ( ok ) in_geometry = 1;
      array_move( tmp_vec1, normal, ndim );
      array_normalize( normal, ndim );
      penetration = tmp - radius;
      for ( idim=0; idim<ndim; idim++ ) {
        if      ( tmp<EPS_COORD )
           projection[idim] = centre[idim];
        else
           projection[idim] = centre[idim] + 
             (radius/tmp)*tmp_vec1[idim];
      }
    }
    else if ( entity==-GEOMETRY_CIRCLE_SMALLSEGMENT ) {
      db( GEOMETRY_CIRCLE_SMALLSEGMENT, index, idum, geometry_circle_smallsegment, 
        ldum, VERSION_NORMAL, GET );
      array_move( geometry_circle_smallsegment, centre, ndim );
      radius = geometry_circle_smallsegment[ndim];
      array_move( &geometry_circle_smallsegment[ndim+1], point_first, ndim );
      array_move( &geometry_circle_smallsegment[ndim+1+ndim], point_second, ndim );
      tolerance = geometry_circle_smallsegment[ndim+1+2*ndim];
      array_subtract( coord, centre, tmp_vec1, ndim );
      tmp = array_size( tmp_vec1, ndim );
      ok = 0;
      if ( projection_type==CONTROL_MESH_DELETE_GEOMETRY ) {
        if ( tmp<=(radius+tolerance+EPS_COORD) )
          ok = 1;
      }
      else {
        if ( tmp>=(radius-tolerance-EPS_COORD) && 
             tmp<=(radius+tolerance+EPS_COORD) )
          ok = 1;
      }
      for ( idim=0; idim<ndim; idim++ ) {
        if ( coord[idim]>point_first[idim] && coord[idim]>point_second[idim] ) ok = 0;
        if ( coord[idim]<point_first[idim] && coord[idim]<point_second[idim] ) ok = 0;
      }         
      if ( ok ) in_geometry = 1;
      array_move( tmp_vec1, normal, ndim );
      array_normalize( normal, ndim );
      penetration = tmp - radius;
      for ( idim=0; idim<ndim; idim++ ) {
        if      ( tmp<EPS_COORD )
           projection[idim] = centre[idim];
        else
           projection[idim] = centre[idim] + 
             (radius/tmp)*tmp_vec1[idim];
      }
    }
    else if ( entity==-GEOMETRY_CYLINDER ) {
      db( GEOMETRY_CYLINDER, index, idum, geometry_cylinder, 
        ldum, VERSION_NORMAL, GET );
      array_subtract( &geometry_cylinder[3], 
        &geometry_cylinder[0], tmp_vec0, ndim);
      if ( array_null( tmp_vec0, ndim ) ) 
        db_error( GEOMETRY_CYLINDER, index );
      cylinder_length = array_size( tmp_vec0, ndim );
      array_normalize( tmp_vec0, ndim );
      array_subtract( coord, &geometry_cylinder[0], tmp_vec1, ndim );
      l = array_inproduct( tmp_vec1, tmp_vec0, ndim );
      radius = geometry_cylinder[6];
      tolerance = geometry_cylinder[7];
      array_multiply( tmp_vec0, tmp_vec0, l, ndim );
      array_subtract( tmp_vec1, tmp_vec0, tmp_vec2, ndim );
      tmp = array_size( tmp_vec2, ndim );
      if ( l>=-EPS_COORD && l<=cylinder_length+EPS_COORD ) {
        if ( projection_type==CONTROL_MESH_DELETE_GEOMETRY ) {
          if ( tmp<=(radius+tolerance+EPS_COORD) )
            in_geometry = 1;
        }
        else {
          if ( tmp>=(radius-tolerance-EPS_COORD) && 
               tmp<=(radius+tolerance+EPS_COORD) )
            in_geometry = 1;
        }
      }
      array_move( tmp_vec2, normal, ndim );
      array_normalize( normal, ndim );
      penetration = tmp - radius;
      for ( idim=0; idim<ndim; idim++ ) {
        if      ( tmp<EPS_COORD )
          projection[idim] = tmp_vec0[idim]+geometry_cylinder[idim];
        else
          projection[idim] = tmp_vec0[idim] +geometry_cylinder[idim]+
         (radius/tmp)*tmp_vec2[idim];
      }
    }
    else if ( entity==-GEOMETRY_CYLINDER_SEGMENT ) {
      db( GEOMETRY_CYLINDER_SEGMENT, index, idum, geometry_cylinder_segment, 
        ldum, VERSION_NORMAL, GET );
      array_subtract( &geometry_cylinder_segment[3], 
        &geometry_cylinder_segment[0], tmp_vec0, ndim);
      if ( array_null( tmp_vec0, ndim ) ) 
        db_error( GEOMETRY_CYLINDER_SEGMENT, index );
      cylinder_length = array_size( tmp_vec0, ndim );
      array_normalize( tmp_vec0, ndim );
      array_subtract( coord, &geometry_cylinder_segment[0], tmp_vec1, ndim );
      l = array_inproduct( tmp_vec1, tmp_vec0, ndim );
      radius = geometry_cylinder_segment[6];
      side[0] = geometry_cylinder_segment[7];
      side[1] = geometry_cylinder_segment[8];
      side[2] = geometry_cylinder_segment[9];
      tolerance = geometry_cylinder_segment[10];
      array_multiply( tmp_vec0, tmp_vec0, l, ndim );
      array_subtract( tmp_vec1, tmp_vec0, tmp_vec2, ndim );
      tmp = array_size( tmp_vec2, ndim );
      if ( l>=-EPS_COORD && l<=cylinder_length+EPS_COORD ) {
        ok = 0;
        if ( projection_type==CONTROL_MESH_DELETE_GEOMETRY ) {
          if ( tmp<=(radius+tolerance+EPS_COORD) )
            ok = 1;
        }
        else {
          if ( tmp>=(radius-tolerance-EPS_COORD) &&
               tmp<=(radius+tolerance+EPS_COORD) )
            ok = 1;
        }                  
        for ( idim=0; idim<ndim; idim++ ) {
          if ( side[idim]>0. && tmp_vec1[idim]<0. ) ok = 0;
          if ( side[idim]<0. && tmp_vec1[idim]>0. ) ok = 0;
        }         
        if ( ok ) in_geometry = 1;
      }
      array_move( tmp_vec2, normal, ndim );
      array_normalize( normal, ndim );
      penetration = tmp - radius;
      for ( idim=0; idim<ndim; idim++ ) {
        if      ( tmp<EPS_COORD )
          projection[idim] = tmp_vec0[idim]+geometry_cylinder_segment[idim];
        else
          projection[idim] = tmp_vec0[idim] +geometry_cylinder_segment[idim]+
           (radius/tmp)*tmp_vec2[idim];
      }
    }
      
// Ellipse contact analysis was modified and corrected by Fernando Lorenzo June 22, 2012
      
    else if ( entity==-GEOMETRY_ELLIPSE ) {
      db( GEOMETRY_ELLIPSE, index, idum, geometry_ellipse, 
        ldum, VERSION_NORMAL, GET );
      array_move( geometry_ellipse, centre, ndim );
// Parameters of the ellipse
      a = geometry_ellipse[ndim];
      b = geometry_ellipse[ndim+1];
      if ( a<=0. || b<=0. ) db_error( GEOMETRY_ELLIPSE, index );
      tolerance = geometry_ellipse[ndim+2]; 
// Coordinates of the nodal point with respect to the center of the ellipse
// Equivalent to displacing the origin of the ellipse to ease the calculations
      x = coord[0] - centre[0];
      y = coord[1] - centre[1];
  
// xyint[0], xyint[1] contain the intersection of a line from the center of the ellipse
// to the ellipse
        xyint[0]=(a*b*x)/(sqrt((a*a)*(y*y)+(b*b)*(x*x)))+coord[0];
        xyint[1]=(a*b*y)/(sqrt((a*a)*(y*y)+(b*b)*(x*x)))+coord[1];
// ellilength is the distance equivalent to the radius of a circle
        ellilength=sqrt(xyint[0]*xyint[0]+xyint[1]*xyint[1]);
// tmp contains the distance from the center of the ellipse to the node
        array_subtract( coord, centre, tmp_vec1, ndim );
        tmp = array_size( tmp_vec1, ndim );
// tmp11 contains the distance from the center of the ellipse to point on the ellipse
        array_subtract( xyint, centre, tmp_vec11, ndim );
        tmp1 = array_size( tmp_vec1, ndim );
        
      if ( projection_type==CONTROL_MESH_DELETE_GEOMETRY ) {
        if ( tmp<=(tmp1+tolerance +EPS_COORD) ) in_geometry = 1;
      }
      else {
          if ( tmp>=(tmp1 -tolerance-EPS_COORD) && 
              tmp<=(tmp1+tolerance+EPS_COORD) ) in_geometry = 1;
      }
       if ( y==0. ) {
         if ( x>0. ) {
          normal[0] = 1.;
          normal[1] = 0.;
        }
        else {
          normal[0] = -1.;
          normal[1] = 0.;
        }
      }
      else if ( x==0. ) {
        if ( y>0. ) {
          normal[0] = 0.;
          normal[1] = 1.;
        }
        else {
              normal[0] = 0.;
              normal[1] = -1.;
        }
      }
      else {
        normal[0] = (b*xyint[0])/a;
        normal[1] = (a*xyint[1])/b;
        array_normalize( normal, ndim );
      }

        array_move( tmp_vec1, normal, ndim );
        array_normalize( normal, ndim );
        penetration = tmp - ellilength;
        for ( idim=0; idim<ndim; idim++ ) {
            if      ( tmp<EPS_COORD )
                projection[idim] = centre[idim];
            else
                projection[idim] = centre[idim] + 
                (ellilength/tmp)*tmp_vec1[idim]; 
        }  
    }
    else if ( entity==-GEOMETRY_LINE ) {
      db( GEOMETRY_LINE, index, idum, geometry_line, 
        ldum, VERSION_NORMAL, GET );
      array_subtract( coord, &geometry_line[0], tmp_vec1, ndim );
      array_subtract( &geometry_line[ndim], &geometry_line[0], vec01, ndim );
      if ( array_null(vec01,ndim) ) db_error( GEOMETRY_LINE, index );
      xi = array_inproduct( tmp_vec1, vec01, ndim ) /
        array_inproduct( vec01, vec01, ndim );
      array_multiply( vec01, tmp_vec1, xi, ndim );
      array_add( &geometry_line[0], tmp_vec1, tmp_vec2, ndim );
      array_subtract( coord, tmp_vec2, tmp_vec1, ndim );
      tmp = array_size( tmp_vec1, ndim );
      tolerance = geometry_line[2*ndim];
      if ( xi>=-EPS_ISO && xi<=(1.+EPS_ISO) &&
          tmp<=tolerance+EPS_COORD ) in_geometry = 1;
      if ( ndim==2 ) {
        array_outproduct_2D( vec01, normal );
        array_normalize( normal, ndim );
        penetration = array_inproduct( normal, tmp_vec1, ndim );
      }
      if ( db_active_index( GEOMETRY_BOUNDA_FACTOR, index, VERSION_NORMAL ) ) {
        db( GEOMETRY_BOUNDA_FACTOR, index, idum, 
          geometry_bounda_factor, length_geometry_bounda_factor, 
          VERSION_NORMAL, GET );
        if ( length_geometry_bounda_factor==2 ) 
          factor = geometry_bounda_factor[0] * (1.-xi) +
          geometry_bounda_factor[1] * xi;
        else {
          assert( length_geometry_bounda_factor==3 );
          xi = 2.*xi - 1.;
          factor = geometry_bounda_factor[0] * 0.5 *( xi*xi - xi ) +
            geometry_bounda_factor[1] * ( 1. - xi*xi ) +
            geometry_bounda_factor[2] * 0.5 * ( xi*xi + xi );
        }
      }
      for ( idim=0; idim<ndim; idim++ ) {
        if ( projection_type==PROJECT_EXACT || tmp<EPS_COORD)
          projection[idim] = tmp_vec2[idim];
        else {
          projection[idim] = tmp_vec2[idim] + (tolerance/tmp)*tmp_vec1[idim];
        }
      }
    }
    else if ( entity==-GEOMETRY_POINT ) {
      db( GEOMETRY_POINT, index, idum, geometry_point, 
        ldum, VERSION_NORMAL, GET );
      array_subtract( coord, geometry_point, tmp_vec1, ndim );
      tmp = array_size( tmp_vec1, ndim );
      tolerance = geometry_point[1*ndim];
      if ( tmp<=tolerance+EPS_COORD ) in_geometry = 1;
      array_move( tmp_vec1, normal, ndim );
      array_normalize( normal, ndim );
      penetration = tmp - tolerance;
      for ( idim=0; idim<ndim; idim++ ) {
        if ( projection_type==PROJECT_EXACT || tmp<EPS_COORD )
          projection[idim] = geometry_point[idim];
        else
          projection[idim] = geometry_point[idim] + 
            (tolerance/tmp)*tmp_vec1[idim];
      }
    }
    else if ( entity==-GEOMETRY_POLYNOMIAL ) {
      length = db_len( GEOMETRY_POLYNOMIAL, index, VERSION_NORMAL );
      if ( length<((ndim-1)*2+2) ) db_error( GEOMETRY_POLYNOMIAL, index );
      geometry_polynomial = 
        db_dbl( GEOMETRY_POLYNOMIAL, index, VERSION_NORMAL );
      if ( ndim==2 ) {
        x  = coord[0];
        x0 = geometry_polynomial[length-3];
        x1 = geometry_polynomial[length-2];
        tolerance = geometry_polynomial[length-1];
        if ( x>=x0 && x<=x1 ) {
          y = geometry_polynomial[0];
          for ( j=level=1; j<(length-3); level++ ) {
            y += geometry_polynomial[j] * scalar_power(x,level); j++;
          }
          dydx[0] = 1.; dydx[1] = geometry_polynomial[1]; dydx[2] = 0.;
          for ( j=level=2; j<(length-3); level++ ) {
            dydx[1] += 
              level * geometry_polynomial[j] * scalar_power(x,level-1); j++;
          }
          if ( coord[1]>=(y-tolerance-EPS_COORD) && 
               coord[1]<=(y+tolerance+EPS_COORD) ) in_geometry = 1;
          normal[1] = 1.;
          penetration = coord[1] - y;
          if ( projection_type==PROJECT_EXACT ) {
            projection[0] = x; projection[1] = y;
          }
          else if ( coord[1]>y ) {
            projection[0] = x; projection[1] = y+tolerance;
          }
          else {
            projection[0] = x; projection[1] = y-tolerance;
          }
        }
      }
      else {
        assert( ndim==3 );
        x  = coord[0];
        y  = coord[1];
        x0 = geometry_polynomial[length-5];
        x1 = geometry_polynomial[length-4];
        y0 = geometry_polynomial[length-3];
        y1 = geometry_polynomial[length-2];
        tolerance = geometry_polynomial[length-1];
        if ( x>=x0 && x<=x1 && y>=y0 && y<=y1 ) {
          j = 0;
          z = geometry_polynomial[j]; j++;
          for ( level=1; j<(length-5); level++ ) {
            if ( j<(length-5) ) {
              z += geometry_polynomial[j] * scalar_power(x,level); j++;
            }
            for ( i=level-1; i>0 && j<(length-5); i-- ) {
              z += geometry_polynomial[j] * scalar_power(x,i) * 
                scalar_power(y,level-i); j++;
            }
            if ( j<(length-5) ) {
              z+= geometry_polynomial[j] * scalar_power(y,level); j++;
            }
          }
          if ( coord[2]>=(z-tolerance-EPS_COORD) && 
               coord[2]<=(z+tolerance+EPS_COORD) ) in_geometry = 1;
          normal[2] = 1.;
          penetration = coord[2] - z;
          if ( projection_type==PROJECT_EXACT ) {
            projection[0] = x; projection[1] = y; projection[2] = z;
          }
          else if ( coord[2]>z ) {
            projection[0] = x; projection[1] = y; projection[2] = z+tolerance;
          }
          else {
            projection[0] = x; projection[1] = y; projection[2] = z-tolerance;
          }
        }
      }
    }
    else if ( entity==-GEOMETRY_SPHERE ) {
      db( GEOMETRY_SPHERE, index, idum, geometry_sphere, 
        ldum, VERSION_NORMAL, GET );
      radius = geometry_sphere[3];
      tolerance = geometry_sphere[4];
      array_move( geometry_sphere, centre, ndim );
      array_subtract( coord, centre, tmp_vec1, ndim );
      tmp = array_size( tmp_vec1, ndim );
      if ( projection_type==CONTROL_MESH_DELETE_GEOMETRY ) {
        if ( tmp<=(radius+tolerance+EPS_COORD) )
          in_geometry = 1;
      }
      else {
        if ( tmp>=(radius-tolerance-EPS_COORD) && 
             tmp<=(radius+tolerance+EPS_COORD) )
          in_geometry = 1;
      }
      array_move( tmp_vec1, normal, ndim );
      array_normalize( normal, ndim );
      penetration = tmp - radius;
      for ( idim=0; idim<ndim; idim++ ) {
        if      ( tmp<EPS_COORD )
          projection[idim] = centre[idim];
        else
          projection[idim] = centre[idim] +
            (radius/tmp)*tmp_vec1[idim];
      }
    }
    else if ( entity==-GEOMETRY_SPHERE_SEGMENT ) {
      db( GEOMETRY_SPHERE_SEGMENT, index, idum, geometry_sphere_segment, 
        ldum, VERSION_NORMAL, GET );
      radius = geometry_sphere_segment[3];
      side[0] = geometry_sphere_segment[4];
      side[1] = geometry_sphere_segment[5];
      side[2] = geometry_sphere_segment[6];
      tolerance = geometry_sphere_segment[7];
      array_move( geometry_sphere_segment, centre, ndim );
      array_subtract( coord, centre, tmp_vec1, ndim );
      tmp = array_size( tmp_vec1, ndim );
      ok = 0;
      if ( projection_type==CONTROL_MESH_DELETE_GEOMETRY ) {
        if ( tmp<=(radius+tolerance+EPS_COORD) )
          ok = 1;
      }
      else {
        if ( tmp>=(radius-tolerance-EPS_COORD) && 
             tmp<=(radius+tolerance+EPS_COORD) )
          ok = 1;
      }
      for ( idim=0; idim<ndim; idim++ ) {
        if ( side[idim]>0. && tmp_vec1[idim]<0. ) ok = 0;
        if ( side[idim]<0. && tmp_vec1[idim]>0. ) ok = 0;
      }
      if ( ok ) in_geometry = 1;
      array_move( tmp_vec1, normal, ndim );
      array_normalize( normal, ndim );
      penetration = tmp - radius;
      for ( idim=0; idim<ndim; idim++ ) {
        if      ( tmp<EPS_COORD )
          projection[idim] = centre[idim];
        else
          projection[idim] = centre[idim] +
            (radius/tmp)*tmp_vec1[idim];
      }
    }
    else if ( entity==-GEOMETRY_TRIANGLE || entity==-GEOMETRY_QUADRILATERAL ) {
      if ( entity==-GEOMETRY_TRIANGLE ) {
        ntriangle = 1;
        db( GEOMETRY_TRIANGLE, index, idum, geometry_triangle, 
          ldum, VERSION_NORMAL, GET );
        eps_iso = EPS_ISO;
        db( GEOMETRY_TRIANGLE_EPSISO, index, idum, &eps_iso, 
          ldum, VERSION_NORMAL, GET_IF_EXISTS );
      }
      else {
        ntriangle = 2;
        db( GEOMETRY_QUADRILATERAL, index, idum, geometry_quadrilateral,
          ldum, VERSION_NORMAL, GET );
        for ( idim=0; idim<ndim; idim++ ) {
            // first triangle
          ind = 0*(3*ndim+1);
          geometry_triangle[ind+0*ndim+idim] =
            geometry_quadrilateral[0*ndim+idim];
          geometry_triangle[ind+1*ndim+idim] =
            geometry_quadrilateral[1*ndim+idim];
          geometry_triangle[ind+2*ndim+idim] =
            geometry_quadrilateral[2*ndim+idim];
            // second triangle
          ind = 1*(3*ndim+1);
          geometry_triangle[ind+0*ndim+idim] =
            geometry_quadrilateral[1*ndim+idim];
          geometry_triangle[ind+1*ndim+idim] =
            geometry_quadrilateral[2*ndim+idim];
          geometry_triangle[ind+2*ndim+idim] =
            geometry_quadrilateral[3*ndim+idim];
        }
        ind = 0*(3*ndim+1);
        geometry_triangle[ind+3*ndim] = geometry_quadrilateral[4*ndim];
        ind = 1*(3*ndim+1);
        geometry_triangle[ind+3*ndim] = geometry_quadrilateral[4*ndim];
      }
      for ( itriangle=0; itriangle<ntriangle && !in_geometry; itriangle++ ) {
        ind = itriangle*(3*ndim+1);
        itest = project_point_on_triangle( coord, &geometry_triangle[ind+0*ndim],
          &geometry_triangle[ind+1*ndim], &geometry_triangle[ind+2*ndim],
          weight );
        l0 = weight[0];
        l1 = weight[1];
        l2 = weight[2];
        array_multiply( &geometry_triangle[ind+0*ndim], tmp_vec1,
          l0, ndim );
        array_multiply( &geometry_triangle[ind+1*ndim], tmp_vec2,
          l1, ndim );
        array_add( tmp_vec1, tmp_vec2, tmp_vec1, ndim );
        array_multiply( &geometry_triangle[ind+2*ndim], tmp_vec2,
          l2, ndim );
        array_add( tmp_vec1, tmp_vec2, tmp_vec3, ndim );
        array_subtract( coord, tmp_vec3, tmp_vec2, ndim );
        tmp = array_size( tmp_vec2, ndim );
        tolerance = geometry_triangle[ind+3*ndim];
        if ( itest &&
             l0>=-eps_iso && l0<(1.+eps_iso) &&
             l1>=-eps_iso && l1<(1.+eps_iso) &&
             l2>=-eps_iso && l2<(1.+eps_iso) ) {
          if ( tmp<=tolerance ) in_geometry = 1;
        }
        array_move( &geometry_triangle[ind+0*ndim], coord0, ndim );
        array_move( &geometry_triangle[ind+1*ndim], coord1, ndim );
        array_move( &geometry_triangle[ind+2*ndim], coord2, ndim );
        array_subtract( coord0, coord1, vec01, ndim );
        array_subtract( coord0, coord2, vec02, ndim );
        array_outproduct_3D( vec01, vec02, normal );
        array_normalize( normal, ndim );
        penetration = array_inproduct( tmp_vec2, normal, ndim );
        if ( db_active_index( GEOMETRY_BOUNDA_FACTOR, index, 
            VERSION_NORMAL ) ) {
          if ( entity==-GEOMETRY_TRIANGLE ) 
            length = 3;
          else {
            assert( entity==-GEOMETRY_QUADRILATERAL );
            length = 4;
          }
          db( GEOMETRY_BOUNDA_FACTOR, index, idum, 
            geometry_bounda_factor, length, VERSION_NORMAL, GET_AND_CHECK );
          if ( itriangle==0 ) {
            factor = geometry_bounda_factor[0] * l0 +
              geometry_bounda_factor[1] * l1 +
              geometry_bounda_factor[2] * l2;
          }
          else {
            assert( itriangle==1 );
            factor = geometry_bounda_factor[1] * l0 +
              geometry_bounda_factor[2] * l1 +
              geometry_bounda_factor[3] * l2;
          }
        }
        for ( idim=0; idim<ndim; idim++ ) {
          if ( projection_type==PROJECT_EXACT || tmp<EPS_COORD )
            projection[idim] = tmp_vec3[idim];
          else {
            projection[idim] = tmp_vec3[idim] + 
              (tolerance/tmp)*tmp_vec2[idim];
          }
        }
      }
    }
    else if ( entity==-NODE_BOUNDARY ) {
      area_node_dataitem();
      if ( db_active_index( NODE_BOUNDARY, inod, version ) )
        in_geometry = 1;
      else
        in_geometry = 0;
    }
    if ( in_geometry ) {
      if ( db_active_index( GEOMETRY_BOUNDA_SINE_X, index, VERSION_NORMAL ) ) {
        db( GEOMETRY_BOUNDA_SINE_X, index, idum, geometry_bounda_sine_x, 
          ldum, VERSION_NORMAL, GET );
        a = geometry_bounda_sine_x[0];
        b = geometry_bounda_sine_x[1];
        factor *= sin(a+b*coord[0]);
      }
      if ( db_active_index( GEOMETRY_BOUNDA_SINE_Y, index, VERSION_NORMAL ) ) {
        db( GEOMETRY_BOUNDA_SINE_Y, index, idum, geometry_bounda_sine_y, 
          ldum, VERSION_NORMAL, GET );
        a = geometry_bounda_sine_y[0];
        b = geometry_bounda_sine_y[1];
        factor *= sin(a+b*coord[1]);
      }
      if ( db_active_index( GEOMETRY_BOUNDA_SINE_Z, index, VERSION_NORMAL ) ) {
        db( GEOMETRY_BOUNDA_SINE_Z, index, idum, geometry_bounda_sine_z, 
          ldum, VERSION_NORMAL, GET );
        a = geometry_bounda_sine_z[0];
        b = geometry_bounda_sine_z[1];
        factor *= sin(a+b*coord[2]);
      }
    }
  }

}

void interpolate_geometry( long int geometry_entity[],
  long int node_numbers[], long int n, double test_coord[], 
  double new_coord[], double test_coord_start_refined[], 
  double new_coord_start_refined[], 
  long int project_type, long int version )

{

  long int i=0, inod=0, part_of_same_geometry=0, in_geometry=0, 
    i_geometry_entity=0, n_geometry_entity=0, geometry_entities[DATA_ITEM_SIZE];
  double rdum=0., ddum[MDIM], test_co[MDIM], new_co[MDIM], diff_co[MDIM];

  if ( geometry_entity[0]==-GEOMETRY_SET ) {
    db( GEOMETRY_SET, geometry_entity[1], geometry_entities, ddum, 
      n_geometry_entity, VERSION_NORMAL, GET );
    n_geometry_entity /= 2;
  }
  else {
    array_move( geometry_entity, geometry_entities, 2 );
    n_geometry_entity = 1;
  }

  for ( i_geometry_entity=0; i_geometry_entity<n_geometry_entity; i_geometry_entity++ ) {
    if ( n>0 ) {
      array_move( test_coord_start_refined, test_co, ndim );
      part_of_same_geometry = 1;
      for ( i=0; i<n && part_of_same_geometry; i++ ) {
        inod = node_numbers[i];
        geometry( inod, ddum, &geometry_entities[i_geometry_entity*2], 
          in_geometry, rdum, ddum, rdum, ddum, NODE_START_REFINED, project_type, version );
        part_of_same_geometry = part_of_same_geometry && in_geometry;
      }
    }
    if ( part_of_same_geometry ) {
      geometry( -1, test_co, &geometry_entities[i_geometry_entity*2], 
        in_geometry, rdum, ddum, rdum, new_co, NODE_START_REFINED, project_type, version );
      array_subtract( new_co, test_co, diff_co, ndim );
      array_move( new_co, new_coord_start_refined, ndim );
      array_add( test_coord, diff_co, new_coord, ndim );
    }
  }

}

void parallel_geometry( void )

{
  long int inod=0, max_node=0, found=0, iloop=0, nloop=0, swit=0,
    ithread=0, *next_of_loop=NULL;
  double factor=0., rdum=0., ddum[MDIM];

  swit = set_swit(-1,-1,"parallel_geometry");
  if ( swit ) pri( "In routine PARALLEL_GEOMETRY" );

  db_max_index( NODE, max_node, VERSION_NORMAL, GET );
  if ( max_node>=0 ) {
    next_of_loop = get_new_int(100+max_node);
    parallel_sys_next_of_loop( next_of_loop, max_node, nloop, ithread );
    for ( iloop=0; iloop<nloop; iloop++ ) {
      inod = next_of_loop[iloop];
      if ( inod>max_node ) {
        break;
      }
      else {
        if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
          geometry( inod, ddum, geometry_ent, found, factor, ddum, rdum,
            ddum, NODE_START_REFINED, PROJECT_EXACT, VERSION_NORMAL );
          if ( found ) nodes_in_geometry[inod] = 1;
        }
      }
    }
    delete[] next_of_loop;
  }

  if ( swit ) pri( "Out routine PARALLEL_GEOMETRY" );
}
