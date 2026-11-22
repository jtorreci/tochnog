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
// 3D (L-T6-U6-measure): safety_slip_ellipsoide and safety_slip_grd sample
// surfaces of triangles; a point with a zero normal (an ellipsoid pole,
// zero area) counts with sigma_n 0, tau 0 and limit c.

#include "tochnog.h"

#define SAFETY_SLIP_EPS_ISO 1.e-3 // point_el tolerance (as post_point)
#define SAFETY_SLIP_LINE_SEGMENT_N 10 // default points per line piece
#define SAFETY_SLIP_CIRCLE_SEGMENT_N 90 // default points per circle
#define SAFETY_SLIP_ELLIPSOIDE_DEFAULT_N 90 // default angles per ellipsoid
#define SAFETY_SLIP_GRD_DEFAULT_N 10 // default points per grd part side
#define SAFETY_SLIP_MVARIATION 100000 // maximum variations of one record

  // one sampled slip surface
struct safety_slip_surface_type {
  long int npoint, mpoint, data_name, data_index, valid, dof_indx;
  double factor;
  double *point, *normal, *direction, *area, *stress, *limit; // 3, 3, 3,
    // 1, 2, 1 per point; direction: unit slip direction (the dof
    // directions of dof_indx>=0 are set in safety_slip_evaluate)
  long int *pair; // npair point pairs (2D)
  long int npair;
  long int *triangle; // ntriangle point triples (3D)
  long int ntriangle;
};

static void safety_slip_surface_allocate( safety_slip_surface_type &s,
  long int mpoint )

{
  s.npoint = 0; s.npair = 0; s.ntriangle = 0; s.valid = 0; s.factor = 0.;
  s.dof_indx = -1; s.mpoint = mpoint;
  s.point = get_new_dbl( 3*mpoint );
  s.normal = get_new_dbl( 3*mpoint );
  s.direction = get_new_dbl( 3*mpoint );
  s.area = get_new_dbl( mpoint );
  s.stress = get_new_dbl( 2*mpoint );
  s.limit = get_new_dbl( mpoint );
  s.pair = get_new_int( 2*mpoint );
  s.triangle = get_new_int( 6*mpoint );
}

static void safety_slip_surface_free( safety_slip_surface_type &s )

{
  delete[] s.point; delete[] s.normal; delete[] s.direction;
  delete[] s.area; delete[] s.stress; delete[] s.limit; delete[] s.pair;
  delete[] s.triangle;
}

static void safety_slip_too_many_points( void )

{
  pri( "Error: too many points on a safety slip surface; maximum ",
    SAFETY_SLIP_MPOINT );
  exit_tn_on_error();
}

static void safety_slip_add_point( safety_slip_surface_type &s,
  double x, double y, double nx, double ny, double area )

{
  if ( s.npoint>=s.mpoint ) safety_slip_too_many_points();
  s.point[3*s.npoint+0] = x; s.point[3*s.npoint+1] = y;
  s.point[3*s.npoint+2] = 0.;
  s.normal[3*s.npoint+0] = nx; s.normal[3*s.npoint+1] = ny;
  s.normal[3*s.npoint+2] = 0.;
    // slip direction: the normal turned clockwise (the piece direction)
  s.direction[3*s.npoint+0] = ny; s.direction[3*s.npoint+1] = -nx;
  s.direction[3*s.npoint+2] = 0.;
  s.area[s.npoint] = area;
  s.npoint++;
}

static void safety_slip_set_point( safety_slip_surface_type &s, long int k,
  double x[], double n[], double d[], double area )

// point k of a 3D surface: coordinates, unit normal (zero for a
// zero-area point), unit slip direction and weight

{
  if ( k>=s.mpoint ) safety_slip_too_many_points();
  array_move( x, &s.point[3*k], 3 );
  array_move( n, &s.normal[3*k], 3 );
  array_move( d, &s.direction[3*k], 3 );
  s.area[k] = area;
  if ( k>=s.npoint ) s.npoint = k+1;
}

static void safety_slip_tangent( double v[], double n[], double d[] )

// d: v projected on the plane normal to the unit vector n, normalized;
// zero when v is (nearly) normal to the plane or zero

{
  long int i=0;
  double vn=array_inproduct( v, n, 3 ), size=array_size( v, 3 ), length=0.;

  for ( i=0; i<3; i++ ) d[i] = v[i] - vn*n[i];
  length = array_size( d, 3 );
  if ( size>0. && length>1.e-10*size ) {
    for ( i=0; i<3; i++ ) d[i] /= length;
  }
  else array_set( d, 0., 3 );
}

static long int safety_slip_ellipsoide_base( double p[], double e[] )

// e1 e2 e3 (e[0..8]) of the ellipsoid parameters p (middle, base1,
// base2, a b c): e1 = base1 normalized, e2 = base2 orthogonalized to e1
// (Gram-Schmidt) and normalized, e3 = e1 x e2; 0 if degenerate

{
  long int i=0;
  double length=0., b2e1=0., size2=array_size( &p[6], 3 );

  array_move( &p[3], &e[0], 3 );
  length = array_size( &e[0], 3 );
  if ( length<=0. ) return 0;
  for ( i=0; i<3; i++ ) e[i] /= length;
  b2e1 = array_inproduct( &p[6], &e[0], 3 );
  for ( i=0; i<3; i++ ) e[3+i] = p[6+i] - b2e1*e[i];
  length = array_size( &e[3], 3 );
  if ( size2<=0. || length<=1.e-10*size2 ) return 0;
  for ( i=0; i<3; i++ ) e[3+i] /= length;
  array_outproduct_3D( &e[0], &e[3], &e[6] );
  return ( p[9]>0. && p[10]>0. && p[11]>0. );
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

static void safety_slip_add_ellipsoide( safety_slip_surface_type &s,
  double p[], double e[], long int n, long int geometric )

// ellipsoid (manual 6.1025-6.1029, L-T6-U6-measure): n angles theta_i =
// 2 pi i/n and phi_j = 2 pi j/n, point i*n+j at local r = (a sin t cos p,
// b sin t sin p, c cos t) in the base e; weight |r_t x r_p| (2 pi/n)^2;
// normal: the normalized gradient (r1/a^2, r2/b^2, r3/c^2); slip
// direction (geometric): e1 projected on the tangent plane. The surface
// is covered twice (theta over 2 pi), which cancels in the safety factor.
// Triangles (i,j),(i,j+1),(i+1,j) and (i,j+1),(i+1,j+1),(i+1,j), the
// indices wrapping in i and j.

{
  long int i=0, j=0, l=0, i1=0, j1=0;
  double h=2.*PIRAD/n, a=p[9], b=p[10], c=p[11], scale=a, t=0., ph=0.,
    r[3], rt[3], rp[3], cross[3], g[3], x[3], normal[3], d[3], cr=0.;

  if ( b>scale ) scale = b;
  if ( c>scale ) scale = c;
  for ( i=0; i<n; i++ ) {
    t = h*i;
    for ( j=0; j<n; j++ ) {
      ph = h*j;
      r[0] = a*sin(t)*cos(ph); r[1] = b*sin(t)*sin(ph); r[2] = c*cos(t);
      rt[0] = a*cos(t)*cos(ph); rt[1] = b*cos(t)*sin(ph); rt[2] = -c*sin(t);
      rp[0] = -a*sin(t)*sin(ph); rp[1] = b*sin(t)*cos(ph); rp[2] = 0.;
      array_outproduct_3D( rt, rp, cross );
      cr = array_size( cross, 3 );
      for ( l=0; l<3; l++ ) {
        x[l] = p[l] + r[0]*e[l] + r[1]*e[3+l] + r[2]*e[6+l];
        g[l] = r[0]/(a*a)*e[l] + r[1]/(b*b)*e[3+l] + r[2]/(c*c)*e[6+l];
      }
      array_set( normal, 0., 3 ); array_set( d, 0., 3 );
      if ( cr>1.e-12*scale*scale ) {
        array_move( g, normal, 3 );
        array_normalize( normal, 3 );
        if ( geometric ) safety_slip_tangent( &e[0], normal, d );
      }
      safety_slip_set_point( s, i*n+j, x, normal, d, cr*h*h );
    }
  }
  for ( i=0; i<n; i++ ) {
    i1 = (i+1)%n;
    for ( j=0; j<n; j++ ) {
      j1 = (j+1)%n;
      long int *tr = &s.triangle[3*s.ntriangle];
      tr[0] = i*n+j; tr[1] = i*n+j1; tr[2] = i1*n+j;
      tr[3] = i*n+j1; tr[4] = i1*n+j1; tr[5] = i1*n+j;
      s.ntriangle += 2;
    }
  }
}

static void safety_slip_add_grd( safety_slip_surface_type &s, long int nx,
  long int ny, double range[], double z[], long int n, double direction[],
  long int geometric )

// Surfer grid z(x,y) (manual 6.1030-6.1033, L-T6-U6-measure): every part
// between four grid points gets n x n points (bilinear in the part, parts
// do not share points); the part normal is the normalized sum of the
// cross products of its triangles 00-10-01 and 10-11-01, the part area
// (the sum of the two triangle areas) is distributed with 2D trapezoid
// weights. Points numbered over the whole grid, x fastest; slip direction
// (geometric): direction projected on the tangent plane.

{
  long int ip=0, jp=0, ii=0, jj=0, l=0, width=(nx-1)*n, k00=0;
  double dx=(range[1]-range[0])/(nx-1), dy=(range[3]-range[2])/(ny-1),
    u=0., v=0., area=0., w[4], c1[3], c2[3], e1[3], e2[3], normal[3], d[3],
    x[3], p00[3], p10[3], p01[3], p11[3];

  for ( jp=0; jp<ny-1; jp++ ) {
    for ( ip=0; ip<nx-1; ip++ ) {
      p00[0] = range[0]+ip*dx; p00[1] = range[2]+jp*dy;
      p00[2] = z[jp*nx+ip];
      p10[0] = p00[0]+dx; p10[1] = p00[1]; p10[2] = z[jp*nx+ip+1];
      p01[0] = p00[0]; p01[1] = p00[1]+dy; p01[2] = z[(jp+1)*nx+ip];
      p11[0] = p10[0]; p11[1] = p01[1]; p11[2] = z[(jp+1)*nx+ip+1];
      array_subtract( p10, p00, e1, 3 ); array_subtract( p01, p00, e2, 3 );
      array_outproduct_3D( e1, e2, c1 );
      array_subtract( p11, p10, e1, 3 ); array_subtract( p01, p10, e2, 3 );
      array_outproduct_3D( e1, e2, c2 );
      area = 0.5*( array_size( c1, 3 ) + array_size( c2, 3 ) );
      for ( l=0; l<3; l++ ) normal[l] = c1[l] + c2[l];
      array_normalize( normal, 3 );
      array_set( d, 0., 3 );
      if ( geometric ) safety_slip_tangent( direction, normal, d );
      for ( jj=0; jj<n; jj++ ) {
        v = (double) jj/(n-1);
        for ( ii=0; ii<n; ii++ ) {
          u = (double) ii/(n-1);
          w[0] = (1.-u)*(1.-v); w[1] = u*(1.-v); w[2] = (1.-u)*v; w[3] = u*v;
          for ( l=0; l<3; l++ )
            x[l] = w[0]*p00[l] + w[1]*p10[l] + w[2]*p01[l] + w[3]*p11[l];
          safety_slip_set_point( s, (jp*n+jj)*width+ip*n+ii, x, normal, d,
            area*( ( ii==0 || ii==n-1 ) ? 0.5 : 1. )*
            ( ( jj==0 || jj==n-1 ) ? 0.5 : 1. )/((n-1)*(n-1)) );
        }
      }
      for ( jj=0; jj<n-1; jj++ ) {
        for ( ii=0; ii<n-1; ii++ ) {
          k00 = (jp*n+jj)*width + ip*n+ii;
          long int *tr = &s.triangle[3*s.ntriangle];
          tr[0] = k00; tr[1] = k00+1; tr[2] = k00+width;
          tr[3] = k00+1; tr[4] = k00+width+1; tr[5] = k00+width;
          s.ntriangle += 2;
        }
      }
    }
  }
}

static void safety_slip_evaluate( safety_slip_surface_type &s )

// stresses, Mohr-Coulomb limit and safety factor of the surface points

{
  long int k=0, element=0, max_element=-1, length=0, name=0, nnol=0,
    inol=0, inod=0, iuknwn=0, nfound=0, group=0, ldum=0, idim=0, idum[1],
    el[1+MNOL];
  double sum_limit=0., sum_tau=0., ddum[1], coords[MNOL*MDIM],
    weight[MNOL], unknowns[MUKNWN], data[DATA_ITEM_SIZE], sig[MDIM*MDIM],
    t[MDIM], sn=0., tau=0., limit=0., sn_sum=0., tau_sum=0., limit_sum=0.,
    t_sum[3], v_sum[3], *element_dof=NULL, *node_dof=NULL;

  db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
  for ( k=0; k<s.npoint; k++ ) {
    double *point = &s.point[3*k], *n = &s.normal[3*k],
      *d = &s.direction[3*k];
      // zero normal: a zero-area point (sigma_n 0, tau 0, limit c)
    long int zero = ( n[0]==0. && n[1]==0. && n[2]==0. );
    nfound = 0; sn_sum = 0.; tau_sum = 0.; limit_sum = 0.;
    array_set( t_sum, 0., 3 ); array_set( v_sum, 0., 3 );
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
      if ( zero ) {
        sn = 0.; tau = 0.; limit = data[1];
      }
      else {
        matrix_ab( sig, n, t, 3, 3, 1 );
        sn = array_inproduct( n, t, 3 );
        tau = array_inproduct( d, t, 3 );
        limit = data[1] - sn*tan(data[0]);
        if ( s.dof_indx>=0 ) {
            // dof slip direction: the displacement or velocity at the
            // point, interpolated from node_dof (these dofs are not in
            // element_dof)
          for ( inol=0; inol<nnol; inol++ ) {
            node_dof = db_dbl( NODE_DOF, el[1+inol], VERSION_NORMAL );
            for ( idim=0; idim<ndim; idim++ ) v_sum[idim] +=
              weight[inol]*node_dof[s.dof_indx+idim*nder];
          }
          for ( idim=0; idim<3; idim++ ) t_sum[idim] += t[idim];
        }
      }
      sn_sum += sn; tau_sum += tau; limit_sum += limit;
      nfound++;
    }
    if ( nfound>0 ) {
      sn = sn_sum/nfound; tau = tau_sum/nfound; limit = limit_sum/nfound;
      if ( s.dof_indx>=0 && !zero ) {
        safety_slip_tangent( v_sum, n, d );
        tau = array_inproduct( d, t_sum, 3 )/nfound;
      }
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

// the internal safety_slip_surface_* records of one evaluated surface;
// elements: point pairs with -none (2D) or point triples (3D)

{
  long int k=0, length=0, idum[1],
    *elements=get_new_int(3*(s.npair+s.ntriangle)+1);
  double ddum[1], *nodes=s.point;

  length = s.npoint;
  db( SAFETY_SLIP_SURFACE_AREAS, isurface, idum, s.area, length,
    VERSION_NORMAL, PUT );
  for ( k=0; k<s.npair; k++ ) {
    elements[3*k+0] = s.pair[2*k+0];
    elements[3*k+1] = s.pair[2*k+1];
    elements[3*k+2] = -NONE;
  }
  array_move( s.triangle, &elements[3*s.npair], 3*s.ntriangle );
  length = 3*(s.npair+s.ntriangle);
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

  // one slip-surface definition (type and index) with its current zone of
  // variations; the zone shrinks around the critical surface in the
  // refinement passes (control_safety_slip N)
struct safety_slip_zone_type {
  long int kind, index, idat_result, nvalue, nvar, nm, nr, segment_n,
    found;
  double first[DATA_ITEM_SIZE], last[DATA_ITEM_SIZE],
    best[DATA_ITEM_SIZE], fs_best;
  long int ibest;
  long int dof_indx; // dof of the slip direction (-1: geometric)
  double direction[3]; // grd: safety_slip_grd_method_direction
  long int grd_nx, grd_ny; // grd: grid points in x and y
  double grd_range[4], *grd; // grd: xmin xmax ymin ymax, z values
};

#define SAFETY_SLIP_CIRCLE_GRID 0
#define SAFETY_SLIP_CIRCLE_LINE 1
#define SAFETY_SLIP_COMBINED 2
#define SAFETY_SLIP_MULTI 3
#define SAFETY_SLIP_ELLIPSOIDE_KIND 4
#define SAFETY_SLIP_GRD_KIND 5
#define SAFETY_SLIP_MZONE 1000

static long int safety_slip_nzone = 0;
static safety_slip_zone_type *safety_slip_zones = NULL;

static void safety_slip_zone_lines( long int kind, long int idat,
  long int idat_n, long int idat_segment_n, long int idat_result )

// zones of safety_slip_multi_linear (polylines, 6.1034) and
// safety_slip_combined_linear (independent point pairs, 6.1021)

{
  long int index=0, max_index=-1, length=0, has_n=0, idum[1];
  double data[DATA_ITEM_SIZE];

  db_max_index( idat, max_index, VERSION_NORMAL, GET );
  for ( index=0; index<=max_index; index++ ) {
    if ( !db_active_index( idat, index, VERSION_NORMAL ) ) continue;
    if ( safety_slip_nzone>=SAFETY_SLIP_MZONE ) db_error( idat, index );
    safety_slip_zone_type &z = safety_slip_zones[safety_slip_nzone];
    db( idat, index, idum, data, length, VERSION_NORMAL, GET );
    has_n = db_active_index( idat_n, index, VERSION_NORMAL );
    z.kind = kind; z.index = index; z.idat_result = idat_result;
    z.dof_indx = -1; z.grd = NULL;
    z.nvar = safety_slip_n( idat_n, index, 1 );
    z.segment_n = safety_slip_n( idat_segment_n, index,
      SAFETY_SLIP_LINE_SEGMENT_N );
      // values of one line: first and last line when _n is given
    z.nvalue = has_n ? length/2 : length;
    if ( z.nvar<1 || z.segment_n<2 || z.nvalue<4 || z.nvalue%2!=0 ||
         ( kind==SAFETY_SLIP_COMBINED && z.nvalue%4!=0 ) ||
         ( has_n && 2*z.nvalue!=length ) )
      db_error( idat, index );
    array_move( data, z.first, z.nvalue );
    array_move( has_n ? &data[z.nvalue] : data, z.last, z.nvalue );
    z.nm = z.nvar; z.nr = 1;
    safety_slip_nzone++;
  }
}

static void safety_slip_zone_circles( long int kind, long int idat_middle,
  long int idat_middle_n, long int idat_radius, long int idat_radius_n,
  long int idat_segment_n, long int idat_result )

// zones of safety_slip_circle_grid_* (6.1009-6.1014) and
// safety_slip_circle_line_* (6.1015-6.1020); first and last hold
// x y r of the first and the last middle and radius

{
  long int index=0, max_index=-1, length=0, idum[1];
  double middle[4], radius[2];

  db_max_index( idat_middle, max_index, VERSION_NORMAL, GET );
  for ( index=0; index<=max_index; index++ ) {
    if ( !db_active_index( idat_middle, index, VERSION_NORMAL ) ) continue;
    if ( safety_slip_nzone>=SAFETY_SLIP_MZONE ) db_error( idat_middle, index );
    safety_slip_zone_type &z = safety_slip_zones[safety_slip_nzone];
    z.kind = kind; z.index = index; z.idat_result = idat_result;
    z.dof_indx = -1; z.grd = NULL;
    db( idat_middle, index, idum, middle, length, VERSION_NORMAL, GET );
    z.nm = safety_slip_n( idat_middle_n, index, 1 );
    if ( ( length!=2 && length!=4 ) || z.nm<1 || ( z.nm>1 && length!=4 ) )
      db_error( idat_middle, index );
    if ( length==2 ) { middle[2] = middle[0]; middle[3] = middle[1]; }
    if ( !db_active_index( idat_radius, index, VERSION_NORMAL ) )
      db_error( idat_radius, index );
    db( idat_radius, index, idum, radius, length, VERSION_NORMAL, GET );
    z.nr = safety_slip_n( idat_radius_n, index, 1 );
    if ( ( length!=1 && length!=2 ) || z.nr<1 || ( z.nr>1 && length!=2 ) )
      db_error( idat_radius, index );
    if ( length==1 ) radius[1] = radius[0];
    z.segment_n = safety_slip_n( idat_segment_n, index,
      SAFETY_SLIP_CIRCLE_SEGMENT_N );
    if ( z.segment_n<3 ) db_error( idat_segment_n, index );
    z.first[0] = middle[0]; z.first[1] = middle[1]; z.first[2] = radius[0];
    z.last[0] = middle[2]; z.last[1] = middle[3]; z.last[2] = radius[1];
    z.nvalue = 3;
    z.nvar = ( kind==SAFETY_SLIP_CIRCLE_GRID ? z.nm*z.nm : z.nm )*z.nr;
    safety_slip_nzone++;
  }
}

static long int safety_slip_dof_indx( long int idat_method, long int index,
  long int geometric )

// the dof of the slip direction of method record idat_method: -1 for the
// geometric method, else materi_displacement (last calculated
// displacement), materi_velocity or materi_velocity_integrated

{
  long int method=geometric, ldum=0;
  double ddum[1];

  db( idat_method, index, &method, ddum, ldum, VERSION_NORMAL,
    GET_IF_EXISTS );
  if ( method==geometric ) return -1;
  else if ( method==-MATERI_DISPLACEMENT && materi_displacement )
    return dis_indx;
  else if ( method==-MATERI_VELOCITY && materi_velocity )
    return vel_indx;
  else if ( method==-MATERI_VELOCITY_INTEGRATED && materi_velocity_integrated )
    return veli_indx;
  db_error( idat_method, index );
  return -1;
}

static void safety_slip_zone_ellipsoide( void )

// zones of safety_slip_ellipsoide (6.1025-6.1029, 3D): first and last
// hold the 12 parameters of the first and the last ellipsoid; with
// safety_slip_ellipsoide_n m every parameter with first != last varies
// over m values, all combinations (m^k surfaces)

{
  long int index=0, max_index=-1, length=0, i=0, idum[1];
  double data[DATA_ITEM_SIZE], e[9];

  db_max_index( SAFETY_SLIP_ELLIPSOIDE, max_index, VERSION_NORMAL, GET );
  for ( index=0; index<=max_index; index++ ) {
    if ( !db_active_index( SAFETY_SLIP_ELLIPSOIDE, index, VERSION_NORMAL ) )
      continue;
    if ( ndim!=3 || safety_slip_nzone>=SAFETY_SLIP_MZONE )
      db_error( SAFETY_SLIP_ELLIPSOIDE, index );
    safety_slip_zone_type &z = safety_slip_zones[safety_slip_nzone];
    db( SAFETY_SLIP_ELLIPSOIDE, index, idum, data, length, VERSION_NORMAL,
      GET );
    z.kind = SAFETY_SLIP_ELLIPSOIDE_KIND; z.index = index;
    z.idat_result = SAFETY_SLIP_ELLIPSOIDE_RESULT; z.grd = NULL;
    z.nm = safety_slip_n( SAFETY_SLIP_ELLIPSOIDE_N, index, 1 );
    if ( ( length!=12 && length!=24 ) || z.nm<1 || ( z.nm>1 && length!=24 ) )
      db_error( SAFETY_SLIP_ELLIPSOIDE, index );
    z.nvalue = 12; z.nr = 1;
    array_move( data, z.first, 12 );
    array_move( length==24 ? &data[12] : data, z.last, 12 );
    if ( !safety_slip_ellipsoide_base( z.first, e ) ||
         !safety_slip_ellipsoide_base( z.last, e ) )
      db_error( SAFETY_SLIP_ELLIPSOIDE, index );
    z.nvar = 1;
    for ( i=0; i<12; i++ ) {
      if ( z.first[i]==z.last[i] ) continue;
      z.nvar *= z.nm;
      if ( z.nvar>SAFETY_SLIP_MVARIATION )
        db_error( SAFETY_SLIP_ELLIPSOIDE_N, index );
    }
    z.segment_n = safety_slip_n( SAFETY_SLIP_ELLIPSOIDE_SEGMENT_N, index,
      SAFETY_SLIP_ELLIPSOIDE_DEFAULT_N );
    if ( z.segment_n<3 || z.segment_n*z.segment_n>SAFETY_SLIP_MPOINT )
      db_error( SAFETY_SLIP_ELLIPSOIDE_SEGMENT_N, index );
    z.dof_indx = safety_slip_dof_indx( SAFETY_SLIP_ELLIPSOIDE_METHOD, index,
      -SAFETY_SLIP_ELLIPSOIDE );
    safety_slip_nzone++;
  }
}

static void safety_slip_grd_read( safety_slip_zone_type &z )

// the Surfer ascii grid <index>.grd of the working directory: DSAA, nx
// ny, xmin xmax, ymin ymax, zmin zmax, then z for all x of the first y,
// of the second y, ...

{
  long int i=0, n=0;
  double zrange[2];
  char str[MCHAR], filename[MCHAR], header[MCHAR];

  strcpy( filename, long_to_a( z.index, str ) );
  strcat( filename, ".grd" );
  ifstream in( filename );
  if ( !in ) {
    pri( "Error: cannot open the Surfer grid file", filename );
    db_error( SAFETY_SLIP_GRD, z.index );
  }
  in.width( MCHAR );
  if ( !(in>>header) || strcmp( header, "DSAA" ) ||
       !(in>>z.grd_nx>>z.grd_ny) || z.grd_nx<2 || z.grd_ny<2 ||
       z.grd_nx>SAFETY_SLIP_MPOINT || z.grd_ny>SAFETY_SLIP_MPOINT ||
       !(in>>z.grd_range[0]>>z.grd_range[1]>>z.grd_range[2]>>
         z.grd_range[3]>>zrange[0]>>zrange[1]) ||
       z.grd_range[1]<=z.grd_range[0] || z.grd_range[3]<=z.grd_range[2] ) {
    pri( "Error: not a Surfer ascii grid (DSAA, nx ny, xmin xmax, ymin ymax, zmin zmax, z values)",
      filename );
    db_error( SAFETY_SLIP_GRD, z.index );
  }
  if ( (z.grd_nx-1)*(z.grd_ny-1)*z.segment_n*z.segment_n>SAFETY_SLIP_MPOINT ) {
    pri( "Error: too many points on the safety_slip_grd surface of",
      filename );
    pri( "The maximum is", SAFETY_SLIP_MPOINT );
    db_error( SAFETY_SLIP_GRD, z.index );
  }
  n = z.grd_nx*z.grd_ny;
  z.grd = get_new_dbl( n );
  for ( i=0; i<n; i++ ) {
    if ( !(in>>z.grd[i]) ) {
      pri( "Error: the Surfer grid file has less z values than nx*ny",
        filename );
      db_error( SAFETY_SLIP_GRD, z.index );
    }
  }
  in.close();
}

static void safety_slip_zone_grd( void )

// zones of safety_slip_grd (6.1030-6.1033, 3D): one surface per index

{
  long int index=0, max_index=-1, length=0, swit=-NO, idum[1];
  double ddum[1];

  db_max_index( SAFETY_SLIP_GRD, max_index, VERSION_NORMAL, GET );
  for ( index=0; index<=max_index; index++ ) {
    if ( !db_active_index( SAFETY_SLIP_GRD, index, VERSION_NORMAL ) )
      continue;
    db( SAFETY_SLIP_GRD, index, &swit, ddum, length, VERSION_NORMAL, GET );
    if ( swit==-NO ) continue;
    if ( swit!=-YES || ndim!=3 || safety_slip_nzone>=SAFETY_SLIP_MZONE )
      db_error( SAFETY_SLIP_GRD, index );
    safety_slip_zone_type &z = safety_slip_zones[safety_slip_nzone];
    z.kind = SAFETY_SLIP_GRD_KIND; z.index = index;
    z.idat_result = SAFETY_SLIP_GRD_RESULT; z.grd = NULL;
    z.nvalue = 0; z.nvar = 1; z.nm = 1; z.nr = 1;
    z.segment_n = safety_slip_n( SAFETY_SLIP_GRD_SEGMENT_N, index,
      SAFETY_SLIP_GRD_DEFAULT_N );
    if ( z.segment_n<2 ) db_error( SAFETY_SLIP_GRD_SEGMENT_N, index );
    z.dof_indx = safety_slip_dof_indx( SAFETY_SLIP_GRD_METHOD, index,
      -SAFETY_SLIP_GRD_DIRECTION );
    if ( z.dof_indx<0 ) {
      if ( !db( SAFETY_SLIP_GRD_METHOD_DIRECTION, index, idum, z.direction,
             length, VERSION_NORMAL, GET_IF_EXISTS ) ||
           array_size( z.direction, 3 )<=0. )
        db_error( SAFETY_SLIP_GRD_METHOD_DIRECTION, index );
    }
    safety_slip_nzone++;
    safety_slip_grd_read( z );
  }
}

static void safety_slip_zone_size( safety_slip_zone_type &z,
  long int &npoint, long int &nelement )

// upper bounds of the points and of the safety_slip_surface_elements
// values of one surface of zone z

{
  long int npiece=0, n=z.segment_n;

  if ( z.kind==SAFETY_SLIP_COMBINED || z.kind==SAFETY_SLIP_MULTI ) {
    npiece = ( z.kind==SAFETY_SLIP_COMBINED ) ? z.nvalue/4 : z.nvalue/2-1;
    npoint = npiece*n; nelement = 3*npiece*(n-1);
  }
  else if ( z.kind==SAFETY_SLIP_ELLIPSOIDE_KIND ) {
    npoint = n*n; nelement = 6*n*n;
  }
  else if ( z.kind==SAFETY_SLIP_GRD_KIND ) {
    npiece = (z.grd_nx-1)*(z.grd_ny-1);
    npoint = npiece*n*n; nelement = 6*npiece*(n-1)*(n-1);
  }
  else {
    npoint = n; nelement = 3*n;
  }
}

static void safety_slip_variation( safety_slip_zone_type &z, long int ivar,
  safety_slip_surface_type &s, double geometry[] )

// the points of variation ivar of zone z; geometry gets its definition
// (line values, or x y r of the circle). Circles: radius outer loop,
// then the middles (x fastest in the grid).

{
  long int ipt=0, step=0, im=0, ir=0, ix=0, iy=0;
  double tx=0., ty=0., tr=0.;

  s.npoint = 0; s.npair = 0; s.ntriangle = 0; s.dof_indx = z.dof_indx;
  if ( z.kind==SAFETY_SLIP_ELLIPSOIDE_KIND ) {
      // the last varied parameter varies fastest
    long int i=0, rest=ivar, digit=0;
    double e[9];
    array_move( z.first, geometry, z.nvalue );
    for ( i=z.nvalue-1; i>=0; i-- ) {
      if ( z.first[i]==z.last[i] ) continue;
      digit = rest%z.nm; rest /= z.nm;
      tx = ( z.nm>1 ) ? (double) digit/(z.nm-1) : 0.;
      geometry[i] = z.first[i] + tx*( z.last[i]-z.first[i] );
    }
    if ( !safety_slip_ellipsoide_base( geometry, e ) )
      db_error( SAFETY_SLIP_ELLIPSOIDE, z.index );
    safety_slip_add_ellipsoide( s, geometry, e, z.segment_n, z.dof_indx<0 );
  }
  else if ( z.kind==SAFETY_SLIP_GRD_KIND )
    safety_slip_add_grd( s, z.grd_nx, z.grd_ny, z.grd_range, z.grd,
      z.segment_n, z.direction, z.dof_indx<0 );
  else if ( z.kind==SAFETY_SLIP_COMBINED || z.kind==SAFETY_SLIP_MULTI ) {
    safety_slip_interpolate( z.first, z.last, z.nvar, ivar, z.nvalue,
      geometry );
    step = ( z.kind==SAFETY_SLIP_COMBINED ) ? 4 : 2;
    for ( ipt=0; ipt+3<z.nvalue; ipt+=step )
      safety_slip_add_piece( s, geometry[ipt], geometry[ipt+1],
        geometry[ipt+2], geometry[ipt+3], z.segment_n );
  }
  else {
    long int nmiddle = ( z.kind==SAFETY_SLIP_CIRCLE_GRID ) ? z.nm*z.nm : z.nm;
    ir = ivar/nmiddle; im = ivar%nmiddle;
    ix = ( z.kind==SAFETY_SLIP_CIRCLE_GRID ) ? im%z.nm : im;
    iy = ( z.kind==SAFETY_SLIP_CIRCLE_GRID ) ? im/z.nm : im;
    tx = ( z.nm>1 ) ? (double) ix/(z.nm-1) : 0.;
    ty = ( z.nm>1 ) ? (double) iy/(z.nm-1) : 0.;
    tr = ( z.nr>1 ) ? (double) ir/(z.nr-1) : 0.;
    geometry[0] = z.first[0] + tx*( z.last[0]-z.first[0] );
    geometry[1] = z.first[1] + ty*( z.last[1]-z.first[1] );
    geometry[2] = z.first[2] + tr*( z.last[2]-z.first[2] );
    safety_slip_add_circle( s, geometry[0], geometry[1], geometry[2],
      z.segment_n );
  }
}

static void safety_slip_shrink( safety_slip_zone_type &z, long int pass )

// the zone of the next refinement pass: the critical surface so far plus
// and minus a part of the current step of every parameter (Professional
// 02-08-2026: half a step after the first pass, a quarter of a step
// after the later passes)

{
  long int i=0, n=0;
  double h = ( pass==1 ) ? 0.5 : 0.25, step=0.;

  if ( !z.found ) return;
  for ( i=0; i<z.nvalue; i++ ) {
    if ( z.kind==SAFETY_SLIP_COMBINED || z.kind==SAFETY_SLIP_MULTI ) n = z.nvar;
    else if ( z.kind==SAFETY_SLIP_ELLIPSOIDE_KIND ) n = z.nm;
    else n = ( i<2 ) ? z.nm : z.nr;
    step = ( n>1 ) ? ( z.last[i]-z.first[i] )/(n-1) : 0.;
    z.first[i] = z.best[i] - h*step;
    z.last[i] = z.best[i] + h*step;
  }
}

void safety_slip( void )

// evaluate all slip surfaces of the current control index. Every pass
// evaluates the zones of all types in the Professional order: in 2D
// circle_grid, circle_line, combined_linear, multi_linear; in 3D grd,
// ellipsoide (the 2D types are not evaluated in 3D); the surfaces of all
// passes are kept and numbered in that order.

{
  long int icontrol=0, control=-NO, npass=1, ipass=0, izone=0, ivar=0,
    isurface=0, nsurface=0, nevaluated=0, length=0, iset=0, max_set=-1, nset=0, i=0,
    mpoint=1, melement=1, npoint=0, nelement=0,
    ldum=0, idum[1], set[DATA_ITEM_SIZE], *surface_index=NULL;
  double ddum[1], geometry[DATA_ITEM_SIZE+1], e[9], *surface_factor=NULL;
  safety_slip_surface_type s;

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  if ( !db( CONTROL_SAFETY_SLIP, icontrol, &control, ddum, ldum,
      VERSION_NORMAL, GET_IF_EXISTS ) ) return;
  if      ( control==-NO ) return;
  else if ( control==-YES ) npass = 1;
  else if ( control>=1 ) npass = control;
  else db_error( CONTROL_SAFETY_SLIP, icontrol );

  db_delete( SAFETY_SLIP_SURFACE_AREAS, VERSION_NORMAL );
  db_delete( SAFETY_SLIP_SURFACE_CRITICAL, VERSION_NORMAL );
  db_delete( SAFETY_SLIP_SURFACE_DATA_INDEX, VERSION_NORMAL );
  db_delete( SAFETY_SLIP_SURFACE_DATA_NAME, VERSION_NORMAL );
  db_delete( SAFETY_SLIP_SURFACE_ELEMENTS, VERSION_NORMAL );
  db_delete( SAFETY_SLIP_SURFACE_FACTOR, VERSION_NORMAL );
  db_delete( SAFETY_SLIP_SURFACE_NODES, VERSION_NORMAL );
  db_delete( SAFETY_SLIP_SURFACE_STRESSES, VERSION_NORMAL );
  db_delete( SAFETY_SLIP_SURFACE_STRESS_LIMIT, VERSION_NORMAL );

  safety_slip_zones = new safety_slip_zone_type[SAFETY_SLIP_MZONE];
  safety_slip_nzone = 0;
    // the line and circle types are 2D surfaces, the grd and ellipsoide
    // types 3D surfaces (the Professional 02-08-2026 skips the line types
    // in 3D and numbers grd before ellipsoide)
  if ( ndim==2 ) {
  safety_slip_zone_circles( SAFETY_SLIP_CIRCLE_GRID,
    SAFETY_SLIP_CIRCLE_GRID_MIDDLE, SAFETY_SLIP_CIRCLE_GRID_MIDDLE_N,
    SAFETY_SLIP_CIRCLE_GRID_RADIUS, SAFETY_SLIP_CIRCLE_GRID_RADIUS_N,
    SAFETY_SLIP_CIRCLE_GRID_SEGMENT_N, SAFETY_SLIP_CIRCLE_GRID_RESULT );
  safety_slip_zone_circles( SAFETY_SLIP_CIRCLE_LINE,
    SAFETY_SLIP_CIRCLE_LINE_MIDDLE, SAFETY_SLIP_CIRCLE_LINE_MIDDLE_N,
    SAFETY_SLIP_CIRCLE_LINE_RADIUS, SAFETY_SLIP_CIRCLE_LINE_RADIUS_N,
    SAFETY_SLIP_CIRCLE_LINE_SEGMENT_N, SAFETY_SLIP_CIRCLE_LINE_RESULT );
  safety_slip_zone_lines( SAFETY_SLIP_COMBINED, SAFETY_SLIP_COMBINED_LINEAR,
    SAFETY_SLIP_COMBINED_LINEAR_N, SAFETY_SLIP_COMBINED_LINEAR_SEGMENT_N,
    SAFETY_SLIP_COMBINED_LINEAR_RESULT );
  safety_slip_zone_lines( SAFETY_SLIP_MULTI, SAFETY_SLIP_MULTI_LINEAR,
    SAFETY_SLIP_MULTI_LINEAR_N, SAFETY_SLIP_MULTI_LINEAR_SEGMENT_N,
    SAFETY_SLIP_MULTI_LINEAR_RESULT );
  }
  else if ( ndim==3 ) {
    safety_slip_zone_grd();
    safety_slip_zone_ellipsoide();
  }

  for ( izone=0; izone<safety_slip_nzone; izone++ ) {
    safety_slip_zones[izone].found = 0;
    safety_slip_zones[izone].ibest = -1;
    nsurface += safety_slip_zones[izone].nvar;
    safety_slip_zone_size( safety_slip_zones[izone], npoint, nelement );
    if ( npoint>mpoint ) mpoint = npoint;
    if ( nelement>melement ) melement = nelement;
  }
  nsurface *= npass;
  surface_index = get_new_int( nsurface+1 );
  surface_factor = get_new_dbl( nsurface+1 );

    // the surface records get the length of the largest surface (they
    // are allocated with that length for every surface)
  if ( mpoint>SAFETY_SLIP_MPOINT ) safety_slip_too_many_points();
  db_data_length_put( SAFETY_SLIP_SURFACE_AREAS, mpoint );
  db_data_length_put( SAFETY_SLIP_SURFACE_ELEMENTS, melement );
  db_data_length_put( SAFETY_SLIP_SURFACE_NODES, 3*mpoint );
  db_data_length_put( SAFETY_SLIP_SURFACE_STRESSES, 2*mpoint );
  db_data_length_put( SAFETY_SLIP_SURFACE_STRESS_LIMIT, mpoint );

  safety_slip_surface_allocate( s, mpoint );
  for ( ipass=1; ipass<=npass; ipass++ ) {
    for ( izone=0; izone<safety_slip_nzone; izone++ ) {
      safety_slip_zone_type &z = safety_slip_zones[izone];
      if ( ipass>1 ) safety_slip_shrink( z, ipass-1 );
      for ( ivar=0; ivar<z.nvar; ivar++ ) {
        safety_slip_variation( z, ivar, s, geometry );
        safety_slip_evaluate( s );
        if ( !s.valid ) continue;
        s.data_name = -z.idat_result;
        s.data_index = z.index;
        safety_slip_put( isurface, s );
        surface_index[isurface] = z.index;
        surface_factor[isurface] = s.factor;
        if ( !z.found || s.factor<z.fs_best ) {
          z.found = 1; z.fs_best = s.factor; z.ibest = isurface;
          array_move( geometry, z.best, z.nvalue );
        }
        isurface++;
      }
    }
  }
  nevaluated = isurface;
  safety_slip_surface_free( s );

    // results per index, critical flags per index or, with
    // safety_slip_set, per set (6.1038-6.1039: the minimum over the
    // surfaces of the indices of the set, all types)
  db_max_index( SAFETY_SLIP_SET, max_set, VERSION_NORMAL, GET );
  for ( izone=0; izone<safety_slip_nzone; izone++ ) {
    safety_slip_zone_type &z = safety_slip_zones[izone];
    if ( !z.found ) continue;
    array_move( z.best, geometry, z.nvalue );
    if ( z.kind==SAFETY_SLIP_ELLIPSOIDE_KIND ) {
        // the normalized bases e1 e2
      safety_slip_ellipsoide_base( z.best, e );
      array_move( e, &geometry[3], 6 );
    }
    geometry[z.nvalue] = z.fs_best;
    length = z.nvalue+1;
    db( z.idat_result, z.index, idum, geometry, length, VERSION_NORMAL, PUT );
    if ( max_set<0 ) safety_slip_put_critical( z.ibest );
  }
  for ( iset=0; iset<=max_set; iset++ ) {
    if ( !db_active_index( SAFETY_SLIP_SET, iset, VERSION_NORMAL ) ) continue;
    db( SAFETY_SLIP_SET, iset, set, ddum, nset, VERSION_NORMAL, GET );
    long int icritical = -1;
    for ( isurface=0; isurface<nevaluated; isurface++ ) {
      for ( i=0; i<nset; i++ ) {
        if ( set[i]==surface_index[isurface] ) break;
      }
      if ( i==nset ) continue;
      if ( icritical<0 || surface_factor[isurface]<surface_factor[icritical] )
        icritical = isurface;
    }
    if ( icritical>=0 ) {
      length = 1;
      db( SAFETY_SLIP_SET_RESULT, iset, idum, &surface_factor[icritical],
        length, VERSION_NORMAL, PUT );
      safety_slip_put_critical( icritical );
    }
  }

  delete[] surface_index; delete[] surface_factor;
  for ( izone=0; izone<safety_slip_nzone; izone++ ) {
    if ( safety_slip_zones[izone].grd!=NULL ) delete[] safety_slip_zones[izone].grd;
  }
  delete[] safety_slip_zones; safety_slip_zones = NULL;
}

static long int safety_slip_print_surface( long int isurface,
  long int critical_only )

// 1 if surface isurface is plotted in GiD (all surfaces, or with
// control_print_gid_safety_slip_critical -yes only the critical ones)

{
  if ( !db_active_index( SAFETY_SLIP_SURFACE_NODES, isurface,
      VERSION_NORMAL ) ) return 0;
  if ( critical_only && !db_active_index( SAFETY_SLIP_SURFACE_CRITICAL,
      isurface, VERSION_NORMAL ) ) return 0;
  return 1;
}

static long int safety_slip_print_critical_only( void )

{
  long int icontrol=0, critical=-NO, ldum=0;
  double ddum[1];

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_PRINT_GID_SAFETY_SLIP_CRITICAL, icontrol, &critical, ddum,
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  return critical==-YES;
}

static long int safety_slip_print_triangles( long int isurface,
  long int *elements, long int &length )

// 1 if the elements of surface isurface are triangles (3D surfaces), 0 for
// point pairs (2D surfaces, -none in every third value)

{
  double ddum[1];

  db( SAFETY_SLIP_SURFACE_ELEMENTS, isurface, elements, ddum, length,
    VERSION_NORMAL, GET );
  return ( length>=3 && elements[2]!=-NONE );
}

void safety_slip_print_gid_mesh( ofstream &outmesh, long int max_node,
  long int max_element )

// GiD mesh of the slip surfaces (manual Professional 6.310): line elements
// between the points of the 2D surfaces, triangles on the 3D surfaces (a
// mesh per element type; the line mesh is SafetySlip, the triangle mesh
// SafetySlip too unless both exist); node and element numbers follow the
// highest node and element of the mesh, the nodes numbered over the
// plotted surfaces in surface order (as the results)

{
  long int isurface=0, max_surface=-1, length=0, nelement=0, k=0, inod=0,
    element=max_element, critical_only=0, itype=0, nline=0, ntriangle=0,
    triangles=0, npoint=0, idum[1],
    *elements=get_new_int(db_data_length(SAFETY_SLIP_SURFACE_ELEMENTS)+1);
  double *nodes=get_new_dbl(db_data_length(SAFETY_SLIP_SURFACE_NODES)+1);

  db_max_index( SAFETY_SLIP_SURFACE_NODES, max_surface, VERSION_NORMAL, GET );
  critical_only = safety_slip_print_critical_only();
  for ( isurface=0; isurface<=max_surface; isurface++ ) {
    if ( !safety_slip_print_surface( isurface, critical_only ) ) continue;
    if ( safety_slip_print_triangles( isurface, elements, nelement ) )
      ntriangle++;
    else
      nline++;
  }
  for ( itype=0; itype<2; itype++ ) {
    if ( ( itype==0 && nline==0 ) || ( itype==1 && ntriangle==0 ) ) continue;
    if ( itype==0 )
      outmesh << "MESH SafetySlip dimension " << ndim
        << " ElemType Linear Nnode 2" << endl;
    else
      outmesh << "MESH " << ( nline>0 ? "SafetySlipTriangles" : "SafetySlip" )
        << " dimension " << ndim << " ElemType Triangle Nnode 3" << endl;
    outmesh << "Coordinates" << endl;
    inod = max_node;
    for ( isurface=0; isurface<=max_surface; isurface++ ) {
      if ( !safety_slip_print_surface( isurface, critical_only ) ) continue;
      triangles = safety_slip_print_triangles( isurface, elements, nelement );
      db( SAFETY_SLIP_SURFACE_NODES, isurface, idum, nodes, length,
        VERSION_NORMAL, GET );
      npoint = length/3;
      if ( triangles==itype ) {
        for ( k=0; k<npoint; k++ ) {
          outmesh << inod+1+k << "\t" << nodes[3*k] << "\t" << nodes[3*k+1];
          if ( ndim==3 ) outmesh << "\t" << nodes[3*k+2];
          outmesh << endl;
        }
      }
      inod += npoint;
    }
    outmesh << "End Coordinates" << endl << "Elements" << endl;
    inod = max_node;
    for ( isurface=0; isurface<=max_surface; isurface++ ) {
      if ( !safety_slip_print_surface( isurface, critical_only ) ) continue;
      triangles = safety_slip_print_triangles( isurface, elements, nelement );
      db( SAFETY_SLIP_SURFACE_NODES, isurface, idum, nodes, length,
        VERSION_NORMAL, GET );
      npoint = length/3;
      if ( triangles==itype ) {
        for ( k=0; k+2<nelement; k+=3 ) {
          element++;
          outmesh << element << "\t" << inod+1+elements[k] << "\t"
            << inod+1+elements[k+1];
          if ( triangles ) outmesh << "\t" << inod+1+elements[k+2];
          outmesh << "\t" << isurface << endl;
        }
      }
      inod += npoint;
    }
    outmesh << "End Elements" << endl;
  }
  delete[] elements; delete[] nodes;
}

void safety_slip_print_gid_results( ofstream &outres, double step_val,
  long int max_node )

// GiD results on the slip surface points: safety_factor (of the surface),
// safety_factor_local (tau_mc / |tau|), safety_sign (sigma_n), safety_sigt
// (|tau|)

{
  long int isurface=0, max_surface=-1, length=0, k=0, inod=0, iresult=0,
    critical_only=0, idum[1];
  double factor=0., value=0.,
    *stresses=get_new_dbl(db_data_length(SAFETY_SLIP_SURFACE_STRESSES)+1),
    *limit=get_new_dbl(db_data_length(SAFETY_SLIP_SURFACE_STRESS_LIMIT)+1);
  const char *result_name[4] = { "safety_factor", "safety_factor_local",
    "safety_sign", "safety_sigt" };

  db_max_index( SAFETY_SLIP_SURFACE_NODES, max_surface, VERSION_NORMAL, GET );
  critical_only = safety_slip_print_critical_only();
  if ( max_surface>=0 ) {
    for ( iresult=0; iresult<4; iresult++ ) {
      outres << "Result  \"" << result_name[iresult] << "\"     \"Load Analysis\" \t"
        << step_val << " Scalar OnNodes" << endl << "Values" << endl;
      inod = max_node;
      for ( isurface=0; isurface<=max_surface; isurface++ ) {
        if ( !safety_slip_print_surface( isurface, critical_only ) ) continue;
        db( SAFETY_SLIP_SURFACE_FACTOR, isurface, idum, &factor, length,
          VERSION_NORMAL, GET );
        db( SAFETY_SLIP_SURFACE_STRESSES, isurface, idum, stresses, length,
          VERSION_NORMAL, GET );
        db( SAFETY_SLIP_SURFACE_STRESS_LIMIT, isurface, idum, limit, length,
          VERSION_NORMAL, GET );
        for ( k=0; k<length; k++ ) {
          if      ( iresult==0 ) value = factor;
          else if ( iresult==1 ) value = ( stresses[2*k+1]>0. ) ?
            limit[k]/stresses[2*k+1] : 0.;
          else if ( iresult==2 ) value = stresses[2*k];
          else value = stresses[2*k+1];
          inod++;
          outres << inod << "\t" << value << endl;
        }
      }
      outres << "End Values" << endl;
    }
  }
  delete[] stresses; delete[] limit;
}
