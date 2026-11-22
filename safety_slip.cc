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

// Slip-surface safety factor (control_safety_slip, manual Professional
// 6.369; F5 T6). For every slip surface, sampled in points k with weight
// A_k, unit normal n_k and slip direction s_k:
//   Fs = sum A_k tau_mc_k / sum A_k |tau_k|,
//   tau_mc = c - sigma_n tan(phi) (Mohr-Coulomb of the element group),
//   sigma_n = n . sigma . n, tau = s . sigma . n.
// The stress in a point is the average, over the elements containing the
// point, of the element nodal stresses (element_dof) interpolated in that
// element. Points outside the mesh or in tension (sigma_n > 0) do not
// count; a surface without counting points gets no result. Measured on
// the Professional 02-08-2026 (odd/tasks/f5-safety-slip.md L-T6-measure).

#include "tochnog.h"

#define SAFETY_SLIP_EPS_ISO 1.e-3 // point_el tolerance (as post_point)
#define SAFETY_SLIP_LINE_SEGMENT_N 10 // default points per line piece
#define SAFETY_SLIP_CIRCLE_SEGMENT_N 90 // default points per circle

  // one sampled slip surface
struct safety_slip_surface_type {
  long int npoint, data_name, data_index, valid;
  double factor;
  double *point, *normal, *area, *stress, *limit; // 3, 3, 1, 2, 1 per point
  long int *pair; // npair point pairs
  long int npair;
};

static void safety_slip_surface_allocate( safety_slip_surface_type &s )

{
  s.npoint = 0; s.npair = 0; s.valid = 0; s.factor = 0.;
  s.point = get_new_dbl( 3*SAFETY_SLIP_MPOINT );
  s.normal = get_new_dbl( 3*SAFETY_SLIP_MPOINT );
  s.area = get_new_dbl( SAFETY_SLIP_MPOINT );
  s.stress = get_new_dbl( 2*SAFETY_SLIP_MPOINT );
  s.limit = get_new_dbl( SAFETY_SLIP_MPOINT );
  s.pair = get_new_int( 2*SAFETY_SLIP_MPOINT );
}

static void safety_slip_surface_free( safety_slip_surface_type &s )

{
  delete[] s.point; delete[] s.normal; delete[] s.area;
  delete[] s.stress; delete[] s.limit; delete[] s.pair;
}

static void safety_slip_add_point( safety_slip_surface_type &s,
  double x, double y, double nx, double ny, double area )

{
  if ( s.npoint>=SAFETY_SLIP_MPOINT ) {
    pri( "Error: too many points on a safety slip surface; maximum ",
      SAFETY_SLIP_MPOINT );
    exit_tn_on_error();
  }
  s.point[3*s.npoint+0] = x; s.point[3*s.npoint+1] = y;
  s.point[3*s.npoint+2] = 0.;
  s.normal[3*s.npoint+0] = nx; s.normal[3*s.npoint+1] = ny;
  s.normal[3*s.npoint+2] = 0.;
  s.area[s.npoint] = area;
  s.npoint++;
}

static void safety_slip_add_piece( safety_slip_surface_type &s,
  double x0, double y0, double x1, double y1, long int segment_n )

// segment_n equidistant points on the piece, trapezoid weights; the
// points of a piece are connected in pairs (Professional layout)

{
  long int i=0, first=s.npoint;
  double dx=x1-x0, dy=y1-y0, L=sqrt(dx*dx+dy*dy), h=0.;

  if ( L<=0. ) return;
  h = L/(segment_n-1);
  for ( i=0; i<segment_n; i++ ) {
    double t = (double) i/(segment_n-1);
    safety_slip_add_point( s, x0+t*dx, y0+t*dy, -dy/L, dx/L,
      ( i==0 || i==segment_n-1 ) ? 0.5*h : h );
    if ( i>0 ) {
      s.pair[2*s.npair+0] = first+i-1;
      s.pair[2*s.npair+1] = first+i;
      s.npair++;
    }
  }
}

static void safety_slip_add_circle( safety_slip_surface_type &s,
  double xm, double ym, double r, long int segment_n )

// segment_n points at theta = 0, 360/n, ... (counterclockwise from +x),
// outward normal, weight 2 pi r / n, closed ring of point pairs

{
  long int i=0, first=s.npoint;
  double theta=0.;

  for ( i=0; i<segment_n; i++ ) {
    theta = 2.*PIRAD*i/segment_n;
    safety_slip_add_point( s, xm+r*cos(theta), ym+r*sin(theta), cos(theta),
      sin(theta), 2.*PIRAD*r/segment_n );
    s.pair[2*s.npair+0] = first+i;
    s.pair[2*s.npair+1] = first+(i+1)%segment_n;
    s.npair++;
  }
}

static void safety_slip_evaluate( safety_slip_surface_type &s )

// stresses, Mohr-Coulomb limit and safety factor of the surface points

{
  long int k=0, element=0, max_element=-1, length=0, name=0, nnol=0,
    inol=0, inod=0, iuknwn=0, nfound=0, group=0, ldum=0, idum[1],
    el[1+MNOL];
  double sum_limit=0., sum_tau=0., ddum[1], coords[MNOL*MDIM],
    weight[MNOL], unknowns[MUKNWN], data[DATA_ITEM_SIZE], sig[MDIM*MDIM],
    t[MDIM], sn=0., tau=0., limit=0., sn_sum=0., tau_sum=0., limit_sum=0.,
    *element_dof=NULL;

  db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
  for ( k=0; k<s.npoint; k++ ) {
    double *point = &s.point[3*k], *n = &s.normal[3*k], d[3];
      // slip direction: the normal turned clockwise (the piece direction)
    d[0] = n[1]; d[1] = -n[0]; d[2] = 0.;
    nfound = 0; sn_sum = 0.; tau_sum = 0.; limit_sum = 0.;
    for ( element=0; element<=max_element; element++ ) {
      if ( !db_active_index( ELEMENT, element, VERSION_NORMAL ) ) continue;
      if ( !db_active_index( ELEMENT_DOF, element, VERSION_NORMAL ) ) continue;
      db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
      name = el[0]; nnol = length - 1;
      for ( inol=0; inol<nnol; inol++ ) {
        inod = el[1+inol];
        db( NODE, inod, idum, &coords[inol*ndim], ldum, VERSION_NORMAL, GET );
      }
      if ( !point_el( point, coords, weight, name, nnol, SAFETY_SLIP_EPS_ISO ) )
        continue;
      group = 0;
      db( ELEMENT_GROUP, element, &group, ddum, ldum, VERSION_NORMAL,
        GET_IF_EXISTS );
      element_dof = db_dbl( ELEMENT_DOF, element, VERSION_NORMAL );
      array_set( unknowns, 0., nuknwn );
      for ( inol=0; inol<nnol; inol++ ) {
        for ( iuknwn=0; iuknwn<nuknwn; iuknwn++ )
          unknowns[iuknwn] += weight[inol]*element_dof[inol*nuknwn+iuknwn];
      }
      if ( !get_group_data( GROUP_MATERI_PLASTI_MOHR_COUL, group, element,
             unknowns, data, ldum, GET_IF_EXISTS ) &&
           !get_group_data( GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT, group,
             element, unknowns, data, ldum, GET_IF_EXISTS ) ) continue;
        // stress tensor (sxx sxy sxz syy syz szz)
      sig[0] = unknowns[stres_indx+0*nder];
      sig[1] = sig[3] = unknowns[stres_indx+1*nder];
      sig[2] = sig[6] = unknowns[stres_indx+2*nder];
      sig[4] = unknowns[stres_indx+3*nder];
      sig[5] = sig[7] = unknowns[stres_indx+4*nder];
      sig[8] = unknowns[stres_indx+5*nder];
      matrix_ab( sig, n, t, 3, 3, 1 );
      sn = array_inproduct( n, t, 3 );
      tau = array_inproduct( d, t, 3 );
      limit = data[1] - sn*tan(data[0]);
      sn_sum += sn; tau_sum += tau; limit_sum += limit;
      nfound++;
    }
    if ( nfound>0 ) {
      sn = sn_sum/nfound; tau = tau_sum/nfound; limit = limit_sum/nfound;
    }
    else {
      sn = 0.; tau = 0.; limit = 0.;
    }
    s.stress[2*k+0] = sn;
    s.stress[2*k+1] = scalar_dabs(tau);
    s.limit[k] = limit;
    if ( nfound>0 && sn<=0. ) {
      sum_limit += s.area[k]*limit;
      sum_tau += s.area[k]*scalar_dabs(tau);
    }
  }
  s.valid = ( sum_tau>0. );
  s.factor = s.valid ? sum_limit/sum_tau : 0.;
}

static long int safety_slip_n( long int idat, long int index, long int deflt )

{
  long int n=deflt, ldum=0;
  double ddum[1];

  db( idat, index, &n, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  return n;
}

static void safety_slip_interpolate( double first[], double last[],
  long int nvar, long int ivar, long int nvalue, double result[] )

// variation ivar of nvar between first and last (equidistant)

{
  long int i=0;
  double t = ( nvar>1 ) ? (double) ivar/(nvar-1) : 0.;

  for ( i=0; i<nvalue; i++ ) result[i] = first[i] + t*( last[i]-first[i] );
}

static void safety_slip_put( long int isurface, safety_slip_surface_type &s )

// the internal safety_slip_surface_* records of one evaluated surface

{
  long int k=0, length=0, idum[1], *elements=get_new_int(3*SAFETY_SLIP_MPOINT);
  double ddum[1], *nodes=s.point;

  length = s.npoint;
  db( SAFETY_SLIP_SURFACE_AREAS, isurface, idum, s.area, length,
    VERSION_NORMAL, PUT );
  for ( k=0; k<s.npair; k++ ) {
    elements[3*k+0] = s.pair[2*k+0];
    elements[3*k+1] = s.pair[2*k+1];
    elements[3*k+2] = -NONE;
  }
  length = 3*s.npair;
  db( SAFETY_SLIP_SURFACE_ELEMENTS, isurface, elements, ddum, length,
    VERSION_NORMAL, PUT );
  length = 3*s.npoint;
  db( SAFETY_SLIP_SURFACE_NODES, isurface, idum, nodes, length,
    VERSION_NORMAL, PUT );
  length = 1;
  db( SAFETY_SLIP_SURFACE_FACTOR, isurface, idum, &s.factor, length,
    VERSION_NORMAL, PUT );
  length = 2*s.npoint;
  db( SAFETY_SLIP_SURFACE_STRESSES, isurface, idum, s.stress, length,
    VERSION_NORMAL, PUT );
  length = s.npoint;
  db( SAFETY_SLIP_SURFACE_STRESS_LIMIT, isurface, idum, s.limit, length,
    VERSION_NORMAL, PUT );
  length = 1;
  db( SAFETY_SLIP_SURFACE_DATA_INDEX, isurface, &s.data_index, ddum,
    length, VERSION_NORMAL, PUT );
  db( SAFETY_SLIP_SURFACE_DATA_NAME, isurface, &s.data_name, ddum,
    length, VERSION_NORMAL, PUT );
  delete[] elements;
}

static void safety_slip_put_critical( long int isurface )

{
  long int yes=-YES, length=1, idum[1];
  double ddum[1];

  idum[0] = yes;
  db( SAFETY_SLIP_SURFACE_CRITICAL, isurface, idum, ddum, length,
    VERSION_NORMAL, PUT );
}

static void safety_slip_lines( long int idat, long int idat_n,
  long int idat_segment_n, long int idat_result, long int combined,
  safety_slip_surface_type &s, long int &isurface )

// safety_slip_multi_linear (polylines, 6.1034) and
// safety_slip_combined_linear (independent point pairs, 6.1021): every
// index gives its variations; the critical one of the index is flagged
// and written to the result record

{
  long int index=0, max_index=-1, length=0, nvar=0, ivar=0, segment_n=0,
    nvalue=0, ipt=0, best=-1, ibest=-1, has_n=0, idum[1];
  double fs_best=0., data[DATA_ITEM_SIZE], line[DATA_ITEM_SIZE],
    result[DATA_ITEM_SIZE+1];

  db_max_index( idat, max_index, VERSION_NORMAL, GET );
  for ( index=0; index<=max_index; index++ ) {
    if ( !db_active_index( idat, index, VERSION_NORMAL ) ) continue;
    db( idat, index, idum, data, length, VERSION_NORMAL, GET );
    has_n = db_active_index( idat_n, index, VERSION_NORMAL );
    nvar = safety_slip_n( idat_n, index, 1 );
    segment_n = safety_slip_n( idat_segment_n, index,
      SAFETY_SLIP_LINE_SEGMENT_N );
      // values of one line: first and last line when _n is given
    nvalue = has_n ? length/2 : length;
    if ( nvar<1 || segment_n<2 || nvalue<4 || nvalue%2!=0 ||
         ( combined && nvalue%4!=0 ) || ( has_n && 2*nvalue!=length ) )
      db_error( idat, index );
    best = -1; ibest = -1; fs_best = 0.;
    for ( ivar=0; ivar<nvar; ivar++ ) {
      safety_slip_interpolate( data, has_n ? &data[nvalue] : data,
        nvar, ivar, nvalue, line );
      s.npoint = 0; s.npair = 0;
      if ( combined ) {
        for ( ipt=0; ipt+3<nvalue; ipt+=4 )
          safety_slip_add_piece( s, line[ipt], line[ipt+1], line[ipt+2],
            line[ipt+3], segment_n );
      }
      else {
        for ( ipt=0; ipt+3<nvalue; ipt+=2 )
          safety_slip_add_piece( s, line[ipt], line[ipt+1], line[ipt+2],
            line[ipt+3], segment_n );
      }
      safety_slip_evaluate( s );
      if ( !s.valid ) continue;
      s.data_name = -idat_result;
      s.data_index = index;
      safety_slip_put( isurface, s );
      if ( best<0 || s.factor<fs_best ) {
        best = ivar; ibest = isurface; fs_best = s.factor;
        array_move( line, result, nvalue );
      }
      isurface++;
    }
    if ( best>=0 ) {
      result[nvalue] = fs_best;
      length = nvalue+1;
      db( idat_result, index, idum, result, length, VERSION_NORMAL, PUT );
      safety_slip_put_critical( ibest );
    }
  }
}

static void safety_slip_circles( long int idat_middle, long int idat_middle_n,
  long int idat_radius, long int idat_radius_n, long int idat_segment_n,
  long int idat_result, long int grid, safety_slip_surface_type &s,
  long int &isurface )

// safety_slip_circle_grid_* (6.1009-6.1014: n x n middles on the
// rectangle first..last, x fastest) and safety_slip_circle_line_*
// (6.1015-6.1020: n middles on the line first..last); radii r_first ..
// r_last (radius_n); surfaces ordered radius, then middle

{
  long int index=0, max_index=-1, length=0, nm=0, nr=0, im=0, ix=0, iy=0,
    ir=0, segment_n=0, ibest=-1, idum[1];
  double middle[4], radius[2], xm=0., ym=0., r=0., tx=0., ty=0.,
    result[4];

  db_max_index( idat_middle, max_index, VERSION_NORMAL, GET );
  for ( index=0; index<=max_index; index++ ) {
    if ( !db_active_index( idat_middle, index, VERSION_NORMAL ) ) continue;
    db( idat_middle, index, idum, middle, length, VERSION_NORMAL, GET );
    nm = safety_slip_n( idat_middle_n, index, 1 );
    if ( ( length!=2 && length!=4 ) || nm<1 || ( nm>1 && length!=4 ) )
      db_error( idat_middle, index );
    if ( length==2 ) { middle[2] = middle[0]; middle[3] = middle[1]; }
    if ( !db_active_index( idat_radius, index, VERSION_NORMAL ) )
      db_error( idat_radius, index );
    db( idat_radius, index, idum, radius, length, VERSION_NORMAL, GET );
    nr = safety_slip_n( idat_radius_n, index, 1 );
    if ( ( length!=1 && length!=2 ) || nr<1 || ( nr>1 && length!=2 ) )
      db_error( idat_radius, index );
    if ( length==1 ) radius[1] = radius[0];
    segment_n = safety_slip_n( idat_segment_n, index,
      SAFETY_SLIP_CIRCLE_SEGMENT_N );
    if ( segment_n<3 ) db_error( idat_segment_n, index );
    ibest = -1;
    for ( ir=0; ir<nr; ir++ ) {
      r = radius[0] + ( nr>1 ? (double) ir/(nr-1) : 0. )*( radius[1]-radius[0] );
      for ( im=0; im<( grid ? nm*nm : nm ); im++ ) {
        ix = grid ? im%nm : im;
        iy = grid ? im/nm : im;
        tx = ( nm>1 ) ? (double) ix/(nm-1) : 0.;
        ty = ( nm>1 ) ? (double) iy/(nm-1) : 0.;
        xm = middle[0] + tx*( middle[2]-middle[0] );
        ym = middle[1] + ty*( middle[3]-middle[1] );
        s.npoint = 0; s.npair = 0;
        safety_slip_add_circle( s, xm, ym, r, segment_n );
        safety_slip_evaluate( s );
        if ( !s.valid ) continue;
        s.data_name = -idat_result;
        s.data_index = index;
        safety_slip_put( isurface, s );
        if ( ibest<0 || s.factor<result[3] ) {
          ibest = isurface;
          result[0] = xm; result[1] = ym; result[2] = r;
          result[3] = s.factor;
        }
        isurface++;
      }
    }
    if ( ibest>=0 ) {
      length = 4;
      db( idat_result, index, idum, result, length, VERSION_NORMAL, PUT );
      safety_slip_put_critical( ibest );
    }
  }
}

void safety_slip( void )

// evaluate all slip surfaces of the current control index; the surfaces
// are numbered over the types in alphabetical order (as the
// Professional): circle_grid, circle_line, combined_linear, multi_linear

{
  long int icontrol=0, control=-NO, isurface=0, ldum=0;
  double ddum[1];
  safety_slip_surface_type s;

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  if ( !db( CONTROL_SAFETY_SLIP, icontrol, &control, ddum, ldum,
      VERSION_NORMAL, GET_IF_EXISTS ) ) return;
  if ( control!=-YES ) {
    if ( control==-NO ) return;
    db_error( CONTROL_SAFETY_SLIP, icontrol );
  }

  db_delete( SAFETY_SLIP_SURFACE_AREAS, VERSION_NORMAL );
  db_delete( SAFETY_SLIP_SURFACE_CRITICAL, VERSION_NORMAL );
  db_delete( SAFETY_SLIP_SURFACE_DATA_INDEX, VERSION_NORMAL );
  db_delete( SAFETY_SLIP_SURFACE_DATA_NAME, VERSION_NORMAL );
  db_delete( SAFETY_SLIP_SURFACE_ELEMENTS, VERSION_NORMAL );
  db_delete( SAFETY_SLIP_SURFACE_FACTOR, VERSION_NORMAL );
  db_delete( SAFETY_SLIP_SURFACE_NODES, VERSION_NORMAL );
  db_delete( SAFETY_SLIP_SURFACE_STRESSES, VERSION_NORMAL );
  db_delete( SAFETY_SLIP_SURFACE_STRESS_LIMIT, VERSION_NORMAL );

  safety_slip_surface_allocate( s );
  safety_slip_circles( SAFETY_SLIP_CIRCLE_GRID_MIDDLE,
    SAFETY_SLIP_CIRCLE_GRID_MIDDLE_N, SAFETY_SLIP_CIRCLE_GRID_RADIUS,
    SAFETY_SLIP_CIRCLE_GRID_RADIUS_N, SAFETY_SLIP_CIRCLE_GRID_SEGMENT_N,
    SAFETY_SLIP_CIRCLE_GRID_RESULT, 1, s, isurface );
  safety_slip_circles( SAFETY_SLIP_CIRCLE_LINE_MIDDLE,
    SAFETY_SLIP_CIRCLE_LINE_MIDDLE_N, SAFETY_SLIP_CIRCLE_LINE_RADIUS,
    SAFETY_SLIP_CIRCLE_LINE_RADIUS_N, SAFETY_SLIP_CIRCLE_LINE_SEGMENT_N,
    SAFETY_SLIP_CIRCLE_LINE_RESULT, 0, s, isurface );
  safety_slip_lines( SAFETY_SLIP_COMBINED_LINEAR, SAFETY_SLIP_COMBINED_LINEAR_N,
    SAFETY_SLIP_COMBINED_LINEAR_SEGMENT_N, SAFETY_SLIP_COMBINED_LINEAR_RESULT,
    1, s, isurface );
  safety_slip_lines( SAFETY_SLIP_MULTI_LINEAR, SAFETY_SLIP_MULTI_LINEAR_N,
    SAFETY_SLIP_MULTI_LINEAR_SEGMENT_N, SAFETY_SLIP_MULTI_LINEAR_RESULT,
    0, s, isurface );
  safety_slip_surface_free( s );
}
