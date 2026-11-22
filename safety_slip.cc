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

void safety_slip( void )

// evaluate all slip surfaces of the current control index

{
  long int icontrol=0, control=-NO, index=0, max_index=-1, length=0,
    nvar=0, ivar=0, segment_n=0, npt=0, ipt=0, isurface=0, icritical=-1,
    ldum=0, idum[1];
  double ddum[1], fs_critical=0., *data=NULL, *line=NULL, *result=NULL;
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
  data = get_new_dbl( DATA_ITEM_SIZE );
  line = get_new_dbl( DATA_ITEM_SIZE );
  result = get_new_dbl( DATA_ITEM_SIZE+1 );

    // safety_slip_multi_linear (6.1034): polylines
  db_max_index( SAFETY_SLIP_MULTI_LINEAR, max_index, VERSION_NORMAL, GET );
  for ( index=0; index<=max_index; index++ ) {
    if ( !db_active_index( SAFETY_SLIP_MULTI_LINEAR, index, VERSION_NORMAL ) )
      continue;
    db( SAFETY_SLIP_MULTI_LINEAR, index, idum, data, length, VERSION_NORMAL,
      GET );
    nvar = safety_slip_n( SAFETY_SLIP_MULTI_LINEAR_N, index, 1 );
    segment_n = safety_slip_n( SAFETY_SLIP_MULTI_LINEAR_SEGMENT_N, index,
      SAFETY_SLIP_LINE_SEGMENT_N );
    npt = db_active_index( SAFETY_SLIP_MULTI_LINEAR_N, index,
      VERSION_NORMAL ) ? length/4 : length/2;
    if ( nvar<1 || segment_n<2 || npt<2 ||
         length!=( db_active_index( SAFETY_SLIP_MULTI_LINEAR_N, index,
           VERSION_NORMAL ) ? 4*npt : 2*npt ) )
      db_error( SAFETY_SLIP_MULTI_LINEAR, index );
    long int best=-1;
    double fs_best=0.;
    for ( ivar=0; ivar<nvar; ivar++ ) {
      safety_slip_interpolate( data, ( length==4*npt ) ? &data[2*npt] : data,
        nvar, ivar, 2*npt, line );
      s.npoint = 0; s.npair = 0;
      for ( ipt=0; ipt+1<npt; ipt++ )
        safety_slip_add_piece( s, line[2*ipt], line[2*ipt+1],
          line[2*ipt+2], line[2*ipt+3], segment_n );
      safety_slip_evaluate( s );
      if ( !s.valid ) continue;
      s.data_name = -SAFETY_SLIP_MULTI_LINEAR_RESULT;
      s.data_index = index;
      safety_slip_put( isurface, s );
      if ( best<0 || s.factor<fs_best ) {
        best = ivar; fs_best = s.factor;
        array_move( line, result, 2*npt );
      }
      if ( icritical<0 || s.factor<fs_critical ) {
        icritical = isurface; fs_critical = s.factor;
      }
      isurface++;
    }
    if ( best>=0 ) {
      result[2*npt] = fs_best;
      length = 2*npt+1;
      db( SAFETY_SLIP_MULTI_LINEAR_RESULT, index, idum, result, length,
        VERSION_NORMAL, PUT );
    }
  }

  if ( icritical>=0 ) {
    long int yes=-YES;
    length = 1;
    db( SAFETY_SLIP_SURFACE_CRITICAL, icritical, &yes, ddum, length,
      VERSION_NORMAL, PUT );
  }

  safety_slip_surface_free( s );
  delete[] data; delete[] line; delete[] result;
}
