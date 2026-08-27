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

// print_beam_force_moment - control_print_beam_force_moment (manual
// Professional 6.262) with the companion records
// control_print_beam_force_moment_coordinates (6.263) and
// control_print_beam_force_moment_switch (6.264).
//
// The forces and moments of the beam / truss / truss-beam elements
// crossed by the cut segment (xstart..xend from
// control_print_beam_force_moment_coordinates; in 2D only x and y) are
// written to the file beam_force_moment.<index> (or
// beam_force_moment.<seq> with -separate_sequential, the static
// counter pattern of print_interface_stress). One line per element.
// The FIRST column is the distance from the cut start point (projection
// of the closest point between the element axis and the cut on the cut
// direction); the following columns are the 12 components in the LOCAL
// beam axes: force_x/y/z and moment_x/y/z of the FIRST node (element
// record node 1) and of the SECOND node, exactly as the manual lists.
//
// Local axes (the SAME frame beam.cc uses for ELEMENT_BEAM_MOMENT):
// local x = element direction node 1 -> node 2; local y = in-plane
// perpendicular of the 2D beam (the rotation matrix of beam_2d);
// local z = out of the beam plane. The out-of-plane components
// (force_z, moment_x, moment_y) are identically zero because the beam
// element is a 2D element (2 in-plane forces + 1 out-of-plane moment
// per node; ELEMENT_BEAM_MOMENT stores NNOL*NDOF = 6 values, in the
// axes of the beam plane).
//
// Selection criterion (decision, documented): an element is printed if
// the minimum distance between its axis segment and the cut segment
// (3D segment-segment closest-point) is smaller than
// BEAM_FORCE_MOMENT_CUT_TOL * max(1., cut_length). The lines are
// sorted by ascending distance (the cut is traversed from xstart to
// xend; equal distances keep element order, stable sort).
//
// Truss (manual: "if the element contains a truss (either a truss
// element or a truss-beam element), the truss force will be used for
// the axial force"): for -truss and -trussbeam elements the axial
// columns are taken from ELEMENT_TRUSS_FORCE, +N at the first node and
// -N at the second node (the element nodal force vector, same
// anti-symmetric pattern as the beam transverse components). A pure
// -truss element has no beam moment: its non-axial columns are zero.
//
// control_print_beam_force_moment_switch -yes multiplies all 12
// components by -1 (manual 6.264). If no element crosses the cut no
// file is written (decision). The file is opened in append mode: each
// step_close appends its lines (pattern of print_interface_stress).

#define BEAM_FORCE_MOMENT_CUT_TOL 1.e-6 // cut tolerance, relative to the cut length

// index_plane of a beam group (pattern of beam_3d): the two in-plane
// axes and the out-of-plane axis (0,1,2 = x,y,z). 2D beams always lie
// in the x-y plane; in 3D the plane comes from group_beam_plane
// (default -x -y).
static void beam_plane_index( long int element_group, long int index_plane[] )

{
  long int idim=0, group_beam_plane[2], ldum=0;
  double ddum[1];

  if ( ndim==2 ) {
    index_plane[0] = 0; index_plane[1] = 1; index_plane[2] = 2;
    return;
  }
  group_beam_plane[0] = -X;
  group_beam_plane[1] = -Y;
  db( GROUP_BEAM_PLANE, element_group, group_beam_plane, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  for ( idim=0; idim<2; idim++ ) {
    if      ( group_beam_plane[idim] == -X ) index_plane[idim] = 0;
    else if ( group_beam_plane[idim] == -Y ) index_plane[idim] = 1;
    else if ( group_beam_plane[idim] == -Z ) index_plane[idim] = 2;
    else db_error( GROUP_BEAM_PLANE, element_group );
  }
  if      ( index_plane[0]==0 && index_plane[1]==1 ) index_plane[2] = 2;
  else if ( index_plane[0]==0 && index_plane[1]==2 ) index_plane[2] = 1;
  else if ( index_plane[0]==1 && index_plane[1]==2 ) index_plane[2] = 0;
  else db_error( GROUP_BEAM_PLANE, element_group );
}

// Closest point between the segment (p1,p2) and the segment (q1,q2)
// (Eberly / Real-Time Collision Detection, closest point between two
// line segments): returns the distance and t_cut in [0,1] such that
// p1 + t_cut*(p2-p1) is the point of the FIRST segment (the cut)
// closest to the second one (the beam). GOTCHA: in the reference
// algorithm the first parameter (s) belongs to the first segment and
// the second one (u) to the second segment; mixing them up swaps the
// closest points.
static double segment_segment_distance( double p1[], double p2[],
  double q1[], double q2[], double &t_cut )

{
  long int idim=0;
  double d1[3], d2[3], r[3], a=0., b=0., c=0., e=0., f=0., s=0., u=0.,
    denom=0., dist=0., closest[3];

  for ( idim=0; idim<3; idim++ ) {
    d1[idim] = p2[idim]-p1[idim];
    d2[idim] = q2[idim]-q1[idim];
    r[idim]  = p1[idim]-q1[idim];
  }
  a = array_inproduct( d1, d1, 3 );
  e = array_inproduct( d2, d2, 3 );
  if ( a<=0. || e<=0. ) return 1.e30; // degenerate segment
  f = array_inproduct( d2, r, 3 );
  b = array_inproduct( d1, d2, 3 );
  c = array_inproduct( d1, r, 3 );
  denom = a*e - b*b;
  if ( denom!=0. )
    s = ( b*f - c*e ) / denom; // s: parameter on the FIRST segment (cut)
  if ( s<0. ) s = 0.;
  else if ( s>1. ) s = 1.;
  u = ( b*s + f ) / e;         // u: parameter on the SECOND segment (beam)
  if ( u<0. ) {
    u = 0.;
    s = -c/a;
    if ( s<0. ) s = 0.;
    else if ( s>1. ) s = 1.;
  }
  else if ( u>1. ) {
    u = 1.;
    s = ( b-c )/a;
    if ( s<0. ) s = 0.;
    else if ( s>1. ) s = 1.;
  }
  t_cut = s;
  for ( idim=0; idim<3; idim++ ) {
    closest[idim] = p1[idim] + s*d1[idim];
    dist += ( closest[idim] - q1[idim] - u*d2[idim] ) *
            ( closest[idim] - q1[idim] - u*d2[idim] );
  }
  return sqrt( dist );
}

void print_beam_force_moment( long int icontrol, long int task )

{
  long int swit=0, ldum=0, idum[1], ncoord=0, element=0, max_element=0,
    element_group=0, length=0, nnol=0, inol=0, idim=0, i=0, j=0,
    n=0, index_plane[3], ip0=0, ip1=0, seq=0, switch_record=-NO,
    *el=NULL, *nodes=NULL;
  double ddum[1], bm[6], truss_force=0., a=0., b=0., f0=0., f1=0.,
    m=0., dist=0., t=0., cut_len=0., tol=0., factor=1., p1[3], p2[3],
    q1[3], q2[3], len=0., *coordinates=NULL, *key=NULL, *val=NULL;
  long int *idx=NULL;
  char filename[MCHAR], str[MCHAR];
  static long int beam_force_moment_seq=0;

  swit = set_swit(-1,-1,"print_beam_force_moment");
  if ( swit ) pri( "In routine PRINT_BEAM_FORCE_MOMENT" );

  if ( task!=-SEPARATE_INDEX && task!=-SEPARATE_SEQUENTIAL )
    db_error( CONTROL_PRINT_BEAM_FORCE_MOMENT, icontrol );

  // cut segment: control_print_beam_force_moment_coordinates (6.263):
  // 2D: xstart ystart xend yend; 3D: xstart ystart zstart xend yend
  // zend.
  coordinates = get_new_dbl(DATA_ITEM_SIZE);
  db( CONTROL_PRINT_BEAM_FORCE_MOMENT_COORDINATES, icontrol, idum,
    coordinates, ncoord, VERSION_NORMAL, GET );
  if ( ncoord!=2*ndim )
    db_error( CONTROL_PRINT_BEAM_FORCE_MOMENT_COORDINATES, icontrol );
  array_set( p1, 0., 3 );
  array_set( p2, 0., 3 );
  for ( idim=0; idim<ndim; idim++ ) {
    p1[idim] = coordinates[idim];
    p2[idim] = coordinates[ndim+idim];
  }
  delete[] coordinates;
  cut_len = array_distance( p1, p2, ddum, 3 );
  if ( cut_len<=0. )
    db_error( CONTROL_PRINT_BEAM_FORCE_MOMENT_COORDINATES, icontrol );
  tol = BEAM_FORCE_MOMENT_CUT_TOL * ( cut_len>1. ? cut_len : 1. );

  // switch (6.264): -yes inverts the sign of all 12 components
  if ( db_active_index( CONTROL_PRINT_BEAM_FORCE_MOMENT_SWITCH, icontrol,
       VERSION_NORMAL ) ) {
    db( CONTROL_PRINT_BEAM_FORCE_MOMENT_SWITCH, icontrol, &switch_record,
      ddum, ldum, VERSION_NORMAL, GET );
    if ( switch_record==-YES ) factor = -1.;
    else if ( switch_record!=-NO )
      db_error( CONTROL_PRINT_BEAM_FORCE_MOMENT_SWITCH, icontrol );
  }

  el = get_new_int(MAXIMUM_NODE+1);
  nodes = get_new_int(MAXIMUM_NODE);
  db_highest_index( ELEMENT, max_element, VERSION_NORMAL );

  // pass 1: count the crossed elements (no file is written if none)
  for ( element=0; element<=max_element; element++ ) {
    if ( !db_active_index( ELEMENT, element, VERSION_NORMAL ) ) continue;
    db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
    if ( el[0]!=-BEAM && el[0]!=-TRUSS && el[0]!=-TRUSSBEAM ) continue;
    nnol = length - 1;
    if ( nnol!=2 ) continue;
    array_move( &el[1], nodes, nnol );
    array_set( q1, 0., 3 );
    array_set( q2, 0., 3 );
    for ( inol=0; inol<nnol; inol++ ) {
      double *c = db_dbl( NODE, nodes[inol], VERSION_NORMAL );
      for ( idim=0; idim<ndim; idim++ ) {
        if ( inol==0 ) q1[idim] = c[idim];
        else           q2[idim] = c[idim];
      }
    }
    if ( segment_segment_distance( p1, p2, q1, q2, t ) < tol ) n++;
  }
  if ( n==0 ) {
    delete[] el;
    delete[] nodes;
    if ( swit )
      pri( "Out function PRINT_BEAM_FORCE_MOMENT (no crossed elements)" );
    return;
  }

  // file name: beam_force_moment.<index> / beam_force_moment.<seq>
  strcpy( filename, "beam_force_moment." );
  if ( task==-SEPARATE_SEQUENTIAL ) {
    seq = beam_force_moment_seq++;
    long_to_a( seq, str );
    strcat( filename, str );
  }
  else {
    long_to_a( icontrol, str );
    strcat( filename, str );
  }

  key = get_new_dbl(n);
  val = get_new_dbl(13*n);
  idx = get_new_int(n);

  // pass 2: collect the lines, key = distance from the cut start point
  n = 0;
  for ( element=0; element<=max_element; element++ ) {
    if ( !db_active_index( ELEMENT, element, VERSION_NORMAL ) ) continue;
    db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
    if ( el[0]!=-BEAM && el[0]!=-TRUSS && el[0]!=-TRUSSBEAM ) continue;
    nnol = length - 1;
    if ( nnol!=2 ) continue;
    array_move( &el[1], nodes, nnol );
    array_set( q1, 0., 3 );
    array_set( q2, 0., 3 );
    for ( inol=0; inol<nnol; inol++ ) {
      double *c = db_dbl( NODE, nodes[inol], VERSION_NORMAL );
      for ( idim=0; idim<ndim; idim++ ) {
        if ( inol==0 ) q1[idim] = c[idim];
        else           q2[idim] = c[idim];
      }
    }
    dist = segment_segment_distance( p1, p2, q1, q2, t );
    if ( dist>=tol ) continue;

    element_group = 0;
    db( ELEMENT_GROUP, element, &element_group, ddum, ldum,
      VERSION_NORMAL, GET_IF_EXISTS );

    // local axes: local x = element direction node 1 -> node 2;
    // (a,b) = the direction components in the beam-plane axes
    // (index_plane[0], index_plane[1]).
    beam_plane_index( element_group, index_plane );
    ip0 = index_plane[0];
    ip1 = index_plane[1];
    len = array_distance( q1, q2, ddum, 3 );
    a = 0.; b = 0.;
    if ( len>0. ) {
      for ( idim=0; idim<ndim; idim++ ) {
        if ( idim==ip0 ) a += ( q2[idim]-q1[idim] ) / len;
        if ( idim==ip1 ) b += ( q2[idim]-q1[idim] ) / len;
      }
    }

    // axial force: truss force (+N at the first node, -N at the second
    // node). Transverse forces and moments from ELEMENT_BEAM_MOMENT,
    // rotated from the beam-plane axes to the local axes with the
    // beam_2d rotation matrix (local x = a*axis0 + b*axis1, local y =
    // -b*axis0 + a*axis1). The out-of-plane components (force_z,
    // moment_x, moment_y) are identically zero for the 2D beam.
    if ( el[0]==-TRUSS ) {
      if ( !db( ELEMENT_TRUSS_FORCE, element, idum, &truss_force,
           ldum, VERSION_NORMAL, GET_IF_EXISTS ) ) continue;
      val[n*13+0] = +truss_force;   // force_x first node
      val[n*13+1] = 0.;             // force_y first node
      val[n*13+2] = 0.;             // force_z first node
      val[n*13+3] = 0.;             // moment_x first node
      val[n*13+4] = 0.;             // moment_y first node
      val[n*13+5] = 0.;             // moment_z first node
      val[n*13+6] = -truss_force;   // force_x second node
      val[n*13+7] = 0.;             // force_y second node
      val[n*13+8] = 0.;             // force_z second node
      val[n*13+9] = 0.;             // moment_x second node
      val[n*13+10] = 0.;            // moment_y second node
      val[n*13+11] = 0.;            // moment_z second node
    }
    else {
      if ( !db( ELEMENT_BEAM_MOMENT, element, idum, bm,
           length, VERSION_NORMAL, GET_IF_EXISTS ) ) continue;
      if ( length<6 ) continue;
      if ( el[0]==-TRUSSBEAM ) {
        if ( !db( ELEMENT_TRUSS_FORCE, element, idum, &truss_force,
             ldum, VERSION_NORMAL, GET_IF_EXISTS ) ) continue;
      }
      else
        truss_force = 0.;
      for ( inol=0; inol<nnol; inol++ ) {
        f0 = bm[inol*3+0];
        f1 = bm[inol*3+1];
        m  = bm[inol*3+2];
        val[n*13+inol*6+0] = ( inol==0 ? +truss_force : -truss_force );
        val[n*13+inol*6+1] = -b*f0 + a*f1;
        val[n*13+inol*6+2] = 0.;
        val[n*13+inol*6+3] = 0.;
        val[n*13+inol*6+4] = 0.;
        val[n*13+inol*6+5] = m;
      }
    }
    // sign switch (6.264): multiply all 12 components by -1
    for ( j=0; j<12; j++ ) val[n*13+j] *= factor;

    key[n] = t * cut_len; // distance from the cut start point
    idx[n] = n;
    n++;
  }

  // sort by ascending distance (stable; equal keys keep element order)
  for ( i=1; i<n; i++ ) {
    j = i;
    while ( j>0 ) {
      long int a_idx = idx[j-1], b_idx = idx[j];
      if ( key[b_idx]>=key[a_idx] ) break;
      idx[j] = a_idx; idx[j-1] = b_idx; j--;
    }
  }

  // write the file (append: one line per element per call)
  ofstream out( filename, ios::app );
  out.precision(TN_PRECISION);
  for ( i=0; i<n; i++ ) {
    long int k = idx[i];
    // snap floating point / solver-residual noise (the equilibrium
    // solve converges to ~1e-6 relative) to zero: the out-of-plane
    // components are identically zero for the 2D beam and the tip
    // moment of the analytic tests is exactly 0
    double maxv = 0.;
    for ( j=0; j<12; j++ )
      if ( fabs(val[k*13+j])>maxv ) maxv = fabs(val[k*13+j]);
    double snap = 1.e-6 * ( 1. + maxv );
    out << key[k];
    for ( j=0; j<12; j++ ) {
      double v = val[k*13+j];
      if ( fabs(v)<snap ) v = 0.;
      out << " " << v;
    }
    out << "\n";
  }
  out.close();

  delete[] el;
  delete[] nodes;
  delete[] key;
  delete[] val;
  delete[] idx;

  if ( swit ) pri( "Out function PRINT_BEAM_FORCE_MOMENT" );
}
