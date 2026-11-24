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

// print_dof_line / print_dof_point - control_print_dof_line (manual
// Professional 6.273-6.280) and control_print_dof_point (6.281-6.283):
// print the values of the node_dof and node_dof_calcul records along a
// line in space (a polyline) or in a point, to plain ASCII files with
// one "x y z <dof>" line per point (in 1D only x). One file per dof
// label, e.g. temp.10 / velx.10 with -separate_index (the "dof" of the
// manual is the dof LABEL, same naming as print_unknowns).

#define DOF_LINE_DEFAULT_N 5
#define DOF_LINE_DEFAULT_EPS_ISO 1.e-3

// Interpolate node_dof (and node_dof_calcul, when present) at one point
// in space. The point is accepted as part of the first element that
// contains it within the eps_iso tolerance (manual 6.276); only elements
// of the optionally given element groups are searched (manual 6.275).
// The geometry used for the point-in-element test is the node_start_refined
// or the node coordinates (manual 6.277). Returns 1 when the point is
// accepted, 0 otherwise. point_dof_calcul may be NULL.
static long int dof_line_interpolate( double point[], long int method,
  long int group_list[], long int ngroups, double eps_iso,
  double point_dof[], double point_dof_calcul[] )

{
  long int element=0, max_element=0, length=0, name=0, inol=0, nnol=0,
    inod=0, idum[1], el[1+MNOL], nodes[MNOL], element_group=0, igroup=0,
    in_group=0, ncalcul=0, icalcul=0, ldum=0;
  double ddum[1], coords[MNOL*MDIM], tmp_node_dof[MUKNWN],
    tmp_node_dof_calcul[MCALCUL], weight[MNOL];

  db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
  for ( element=0; element<=max_element; element++ ) {
    if ( !db_active_index( ELEMENT, element, VERSION_NORMAL ) ) continue;
    if ( ngroups>0 ) {
      db( ELEMENT_GROUP, element, &element_group, ddum, ldum,
        VERSION_NORMAL, GET_IF_EXISTS );
      in_group = 0;
      for ( igroup=0; igroup<ngroups; igroup++ )
        if ( group_list[igroup]==element_group ) in_group = 1;
      if ( !in_group ) continue;
    }
    db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
    name = el[0]; nnol = length - 1; array_move( &el[1], nodes, nnol );
    for ( inol=0; inol<nnol; inol++ ) {
      inod = nodes[inol];
      // -node_start_refined falls back to the stored node coordinates
      // when the record does not exist (geometrically linear analysis)
      if ( method==-NODE_START_REFINED && db_active_index(
           NODE_START_REFINED, inod, VERSION_NORMAL ) )
        db( NODE_START_REFINED, inod, idum, &coords[inol*ndim], ldum,
          VERSION_NORMAL, GET );
      else
        db( NODE, inod, idum, &coords[inol*ndim], ldum, VERSION_NORMAL, GET );
    }
    if ( point_el( point, coords, weight, name, nnol, eps_iso ) ) {
      array_set( point_dof, 0., nuknwn );
      if ( point_dof_calcul ) {
        ncalcul = db_len( POST_CALCUL_SCAL_VEC_MAT, 0, VERSION_NORMAL );
        array_set( point_dof_calcul, 0., MCALCUL );
      }
      for ( inol=0; inol<nnol; inol++ ) {
        inod = nodes[inol];
        if ( db_active_index( NODE_DOF, inod, VERSION_NEW ) )
          db( NODE_DOF, inod, idum, tmp_node_dof, ldum, VERSION_NEW, GET );
        else
          db( NODE_DOF, inod, idum, tmp_node_dof, ldum, VERSION_NORMAL, GET );
        array_multiply( tmp_node_dof, tmp_node_dof, weight[inol], nuknwn );
        array_add( tmp_node_dof, point_dof, point_dof, nuknwn );
        if ( point_dof_calcul && db_active_index( NODE_DOF_CALCUL, inod,
             VERSION_NORMAL ) ) {
          db( NODE_DOF_CALCUL, inod, idum, tmp_node_dof_calcul, ldum,
            VERSION_NORMAL, GET );
          for ( icalcul=0; icalcul<ncalcul; icalcul++ )
            point_dof_calcul[icalcul] += weight[inol]*
              tmp_node_dof_calcul[icalcul];
        }
      }
      return 1;
    }
  }
  return 0;
}

// Position of point ipoint (0..n-1) of a line of n points, distributed
// with equal spacing over the TOTAL length of the polyline defined by
// coordinates (ncoord values = nvertex*ndim). n==1 -> the start point.
// A degenerate polyline (zero total length) collapses to its first
// vertex.
static void dof_line_point_position( double coordinates[], long int ncoord,
  long int n, long int ipoint, double point[] )

{
  long int nvertex=0, iv=0, idim=0;
  double total=0., s=0., target=0., seglen=0., dir[MDIM];

  nvertex = ncoord/ndim;
  if ( nvertex<1 ) return;
  if ( n<=1 ) {
    array_move( &coordinates[0], point, ndim );
    return;
  }
  for ( iv=0; iv<nvertex-1; iv++ ) {
    array_subtract( &coordinates[(iv+1)*ndim], &coordinates[iv*ndim],
      dir, ndim );
    total += array_size( dir, ndim );
  }
  if ( total<=0. ) {
    array_move( &coordinates[0], point, ndim );
    return;
  }
  target = ( (double) ipoint ) * total / ( (double) (n-1) );
  s = 0.;
  for ( iv=0; iv<nvertex-1; iv++ ) {
    array_subtract( &coordinates[(iv+1)*ndim], &coordinates[iv*ndim],
      dir, ndim );
    seglen = array_size( dir, ndim );
    if ( s+seglen >= target-1.e-12*total ) {
      if ( seglen<=0. ) {
        array_move( &coordinates[iv*ndim], point, ndim );
        return;
      }
      for ( idim=0; idim<ndim; idim++ )
        point[idim] = coordinates[iv*ndim+idim] +
          dir[idim]*(target-s)/seglen;
      return;
    }
    s += seglen;
  }
  array_move( &coordinates[(nvertex-1)*ndim], point, ndim );
}

// Append the file name extension: .<icontrol> for -separate_index and
// -yes, .<seq> for -separate_sequential (one number per print call).
static void dof_line_extension( char filename[], long int icontrol,
  long int task, long int seq )

{
  char str[MCHAR];

  if ( task==-SEPARATE_INDEX || task==-YES ) {
    long_to_a( icontrol, str );
    strcat( filename, "." );
    strcat( filename, str );
  }
  else {
    assert( task==-SEPARATE_SEQUENTIAL );
    long_to_a( seq, str );
    strcat( filename, "." );
    strcat( filename, str );
  }
}

static void print_dof_line_point( long int icontrol, long int task,
  long int is_line )

{
  long int swit=0, ldum=0, idum[1], ipuknwn=0, iuknwn=0, nder_=0,
    npuknwn_=0, npoints=0, n=DOF_LINE_DEFAULT_N, ipoint=0,
    idim=0, ncalcul=0, icalcul=0, method=-NODE_START_REFINED,
    move=-NO, time=-NO, ncoord=0, ngroups=0, seq=0,
    nfound=0, *dof_label=NULL, *group_list=NULL, *point_found=NULL;
  double ddum[1], eps_iso=DOF_LINE_DEFAULT_EPS_ISO, time_current=0.,
    dtime=0., *coordinates=NULL, *point=NULL, *point_dof=NULL,
    *point_dof_calcul=NULL, *point_coords=NULL;
  char filename[MCHAR];
  static long int dof_line_seq=0, dof_point_seq=0;

  if ( is_line ) {
    swit = set_swit(-1,-1,"print_dof_line");
    if ( swit ) pri( "In routine PRINT_DOF_LINE" );
  }
  else {
    swit = set_swit(-1,-1,"print_dof_point");
    if ( swit ) pri( "In routine PRINT_DOF_POINT" );
  }

  if ( task!=-YES && task!=-SEPARATE_INDEX && task!=-SEPARATE_SEQUENTIAL ) {
    if ( is_line ) db_error( CONTROL_PRINT_DOF_LINE, icontrol );
    else db_error( CONTROL_PRINT_DOF_POINT, icontrol );
  }

  if ( nuknwn<=0 ) return;

  nder_ = nder;
  npuknwn_ = npuknwn;

  // --- read the records of this icontrol ---------------------------
  if ( is_line ) {
    coordinates = get_new_dbl(DATA_ITEM_SIZE);
    db( CONTROL_PRINT_DOF_LINE_COORDINATES, icontrol, idum, coordinates,
      ncoord, VERSION_NORMAL, GET );
    if ( ncoord%ndim!=0 || ncoord<2*ndim )
      db_error( CONTROL_PRINT_DOF_LINE_COORDINATES, icontrol );
    if ( db_active_index( CONTROL_PRINT_DOF_LINE_N, icontrol,
         VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_DOF_LINE_N, icontrol, &n, ddum, ldum,
        VERSION_NORMAL, GET );
      if ( n<1 ) db_error( CONTROL_PRINT_DOF_LINE_N, icontrol );
    }
    else
      n = DOF_LINE_DEFAULT_N;
    if ( db_active_index( CONTROL_PRINT_DOF_LINE_EPS_ISO, icontrol,
         VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_DOF_LINE_EPS_ISO, icontrol, idum, &eps_iso, ldum,
        VERSION_NORMAL, GET );
      if ( eps_iso<0. ) db_error( CONTROL_PRINT_DOF_LINE_EPS_ISO, icontrol );
    }
    if ( db_active_index( CONTROL_PRINT_DOF_LINE_METHOD, icontrol,
         VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_DOF_LINE_METHOD, icontrol, &method, ddum, ldum,
        VERSION_NORMAL, GET );
      if ( method!=-NODE && method!=-NODE_START_REFINED )
        db_error( CONTROL_PRINT_DOF_LINE_METHOD, icontrol );
    }
    if ( db_active_index( CONTROL_PRINT_DOF_LINE_ELEMENT_GROUP, icontrol,
         VERSION_NORMAL ) ) {
      group_list = get_new_int(DATA_ITEM_SIZE);
      db( CONTROL_PRINT_DOF_LINE_ELEMENT_GROUP, icontrol, group_list,
        ddum, ngroups, VERSION_NORMAL, GET );
    }
    if ( db_active_index( CONTROL_PRINT_DOF_LINE_MOVE, icontrol,
         VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_DOF_LINE_MOVE, icontrol, &move, ddum, ldum,
        VERSION_NORMAL, GET );
      if ( move!=-YES && move!=-NO )
        db_error( CONTROL_PRINT_DOF_LINE_MOVE, icontrol );
    }
    if ( db_active_index( CONTROL_PRINT_DOF_LINE_TIME, icontrol,
         VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_DOF_LINE_TIME, icontrol, &time, ddum, ldum,
        VERSION_NORMAL, GET );
      if ( time!=-YES && time!=-NO )
        db_error( CONTROL_PRINT_DOF_LINE_TIME, icontrol );
    }
    npoints = n;
  }
  else {
    coordinates = get_new_dbl(DATA_ITEM_SIZE);
    db( CONTROL_PRINT_DOF_POINT_COORDINATES, icontrol, idum, coordinates,
      ncoord, VERSION_NORMAL, GET );
    if ( ncoord<ndim ) db_error( CONTROL_PRINT_DOF_POINT_COORDINATES, icontrol );
    if ( db_active_index( CONTROL_PRINT_DOF_POINT_TIME, icontrol,
         VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_DOF_POINT_TIME, icontrol, &time, ddum, ldum,
        VERSION_NORMAL, GET );
      if ( time!=-YES && time!=-NO )
        db_error( CONTROL_PRINT_DOF_POINT_TIME, icontrol );
    }
    npoints = 1;
  }

  db( TIME_CURRENT, 0, idum, &time_current, ldum, VERSION_NORMAL, GET );

  // --- interpolate the dofs at every point -------------------------
  point = get_new_dbl(MDIM);
  point_coords = get_new_dbl(npoints*MDIM);
  point_dof = get_new_dbl(npoints*MUKNWN);
  point_found = get_new_int(npoints);
  if ( db_active_index( POST_CALCUL, 0, VERSION_NORMAL ) )
    point_dof_calcul = get_new_dbl(npoints*MCALCUL);

  dof_label = get_new_int(MUKNWN);
  db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET );

  ncalcul = 0;
  if ( point_dof_calcul )
    ncalcul = db_len( POST_CALCUL_SCAL_VEC_MAT, 0, VERSION_NORMAL );

  for ( ipoint=0; ipoint<npoints; ipoint++ ) {
    if ( is_line )
      dof_line_point_position( coordinates, ncoord, npoints, ipoint,
        &point_coords[ipoint*MDIM] );
    else
      array_move( &coordinates[0], &point_coords[ipoint*MDIM], ndim );
    if ( dof_line_interpolate( &point_coords[ipoint*MDIM], method,
         group_list, ngroups, eps_iso, &point_dof[ipoint*MUKNWN],
         point_dof_calcul ? &point_dof_calcul[ipoint*MCALCUL] : NULL ) ) {
      point_found[ipoint] = 1;
      nfound++;
    }
    else {
      // point outside the mesh (or outside the element groups): the
      // point is omitted from the printed files (documented decision:
      // consistent with "the point is accepted to be part of an element")
      point_found[ipoint] = 0;
    }
  }

  // --- write one file per dof label --------------------------------
  if ( task==-SEPARATE_SEQUENTIAL ) {
    if ( is_line ) seq = dof_line_seq++;
    else seq = dof_point_seq++;
  }

  for ( ipuknwn=0; ipuknwn<npuknwn_; ipuknwn++ ) {
    iuknwn = ipuknwn*nder_;
    strcpy( filename, db_name(dof_label[iuknwn]) );
    dof_line_extension( filename, icontrol, task, seq );
    ofstream out( filename, ios::app );
    out.precision(TN_PRECISION);
    if ( time==-YES )
      out << "# time " << time_current << "\n";
    for ( ipoint=0; ipoint<npoints; ipoint++ ) {
      if ( !point_found[ipoint] ) continue;
      for ( idim=0; idim<ndim; idim++ )
        out << point_coords[ipoint*MDIM+idim] << " ";
      out << point_dof[ipoint*MUKNWN+iuknwn] << "\n";
    }
    out.close();
  }

  // --- one file per node_dof_calcul item ---------------------------
  for ( icalcul=0; icalcul<ncalcul; icalcul++ ) {
    strcpy( filename, post_calcul_names_without_extension[icalcul] );
    dof_line_extension( filename, icontrol, task, seq );
    ofstream out( filename, ios::app );
    out.precision(TN_PRECISION);
    if ( time==-YES )
      out << "# time " << time_current << "\n";
    for ( ipoint=0; ipoint<npoints; ipoint++ ) {
      if ( !point_found[ipoint] ) continue;
      for ( idim=0; idim<ndim; idim++ )
        out << point_coords[ipoint*MDIM+idim] << " ";
      out << point_dof_calcul[ipoint*MCALCUL+icalcul] << "\n";
    }
    out.close();
  }

  // --- move: follow the material particles with the velocity field --
  // (manual 6.278; only meaningful when materi_velocity is initialized)
  if ( is_line && move==-YES && materi_velocity && nfound>0 ) {
    if ( db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET_IF_EXISTS ) ||
         db( DTIME, 0, idum, &dtime, ldum, VERSION_NORMAL, GET_IF_EXISTS ) ) {
      for ( ipoint=0; ipoint<npoints; ipoint++ ) {
        if ( !point_found[ipoint] ) continue;
        for ( idim=0; idim<ndim; idim++ )
          coordinates[ipoint*ndim+idim] +=
            point_dof[ipoint*MUKNWN+vel_indx+idim*nder_]*dtime;
      }
      db( CONTROL_PRINT_DOF_LINE_COORDINATES, icontrol, idum, coordinates,
        ncoord, VERSION_NORMAL, PUT );
    }
  }

  delete[] coordinates;
  delete[] point;
  delete[] point_coords;
  delete[] point_dof;
  delete[] point_found;
  if ( point_dof_calcul ) delete[] point_dof_calcul;
  delete[] dof_label;
  if ( group_list ) delete[] group_list;

  if ( swit ) {
    if ( is_line ) pri( "Out routine PRINT_DOF_LINE" );
    else pri( "Out routine PRINT_DOF_POINT" );
  }
}

void print_dof_line( long int icontrol, long int task )

{
  print_dof_line_point( icontrol, task, 1 );
}

void print_dof_point( long int icontrol, long int task )

{
  print_dof_line_point( icontrol, task, 0 );
}
