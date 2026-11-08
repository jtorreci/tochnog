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

// ---------------------------------------------------------------------
// geometry_moving family (manual Professional 6.533-6.538) and
// control_geometry_moving (6.135), F5 T4 U8.
//
// A geometrical entity (-geometry_point with a radius, -geometry_triangle,
// -geometry_quadrilateral, -geometry_hexahedral, -geometry_tetrahedral)
// comes into existence at start_time, translates with a constant
// velocity until end_time and then stays where it is. The elements it
// passes are excavated: an element excavated for a part p ends with the
// delete factor 1 - p (the element_delete_factor machinery of
// control_mesh_delete_geometry_factor: frozen released-stress load
// scaled by the factor, no stiffness), an element excavated completely
// is deleted.
//
// control_geometry_moving index -initialise (geometry_moving_initialise,
// once at its control index): for every element, nspace sample points
// per parametric direction (manual 6.538: nspace, nspace^2, nspace^3
// points for 1D/2D/3D elements) are tested against the entity at the
// ntime+1 time moments t_i = start_time + i*(end_time-start_time)/ntime
// with geometry() (a scratch geometry record of the entity at its start
// position; the point is moved back over the translation instead of
// moving the entity). The part is the volume-weighted fraction of the
// sample points that are inside the entity at any of these moments
// (excavation is permanent). Measured on the Professional 02-08-2026 and
// followed where the manual is silent:
//   - per sample point the entry time is the time moment BEFORE the
//     first moment it is inside (t_{i-1}; t_0 when inside at t_0), the
//     leave time the last moment it is inside before it is outside;
//   - element_geometry_moving_in_time = first and last entry time,
//     element_geometry_moving_out_time = first and last leave time (only
//     when a point leaves);
//   - geometry_moving_apply (every step_start): the delete factor ramps
//     linearly from 1 at in_time[0] to 1 - part at in_time[1] (at once
//     when they are equal) and stays there; a factor 0 deletes the
//     element (and the nodes no element uses any more).
// The Professional tests the element NODES instead of the sample points
// (the part is the fraction of nodes ever inside, nspace has no effect);
// the GNU follows the manual (DIVERGENCES GEOMETRY-MOVING-SAMPLING).
//
// Sample points (hardcoded rule, see the developer manual):
//   - tensor-product elements (-bar2/3/4, -truss, -beam, -trussbeam with
//     2 nodes, -quad4/9/16, -hex8/27): the midpoints of nspace equal
//     cells per isoparametric direction, weighted by the Jacobian measure;
//   - simplices (-tria3/6, -tet4/10): collapsed (Duffy) square/cube cell
//     midpoints on the corner nodes (straight-sided element), weighted by
//     the collapse Jacobian;
//   - any other element (prism, interface, spring, ...): its nodes, equal
//     weights (the Professional rule).
// ---------------------------------------------------------------------

#include "tochnog.h"

#define GM_EPS_FACTOR 1.e-6     // as EPS_ELEMENT_DELETE_FACTOR (delete.cc)
#define GM_EPS_TIME 1.e-12
#define GM_NTIME_DEFAULT 100
#define GM_NSPACE_DEFAULT 5
#define GM_BOX_MARGIN 1.e-2     // relative margin of the bounding-box prefilter

struct gm_entity_data {
  long int index, entity, nvertex, ntime, nspace;
  double t_start, t_end, velocity[MDIM], box_min[MDIM], box_max[MDIM];
};

static long int gm_entity_vertices( long int entity )
{
  if ( entity==-GEOMETRY_POINT ) return 1;
  if ( entity==-GEOMETRY_TRIANGLE ) return 3;
  if ( entity==-GEOMETRY_QUADRILATERAL ) return 4;
  if ( entity==-GEOMETRY_TETRAHEDRAL ) return 4;
  if ( entity==-GEOMETRY_HEXAHEDRAL ) return 8;
  return 0;
}

// read and check the records of geometry_moving index
static void gm_read( long int index, gm_entity_data &gm, double parameter[] )
{
  long int length=0, ldum=0, idim=0, ivertex=0, operat=0, n[2],
    idum[1];
  double ddum[1], time[2], velocity[MDIM], radius=0., margin=0.;

  gm.index = index;
  db( GEOMETRY_MOVING, index, &gm.entity, ddum, ldum, VERSION_NORMAL, GET );
  gm.nvertex = gm_entity_vertices( gm.entity );
  if ( gm.nvertex==0 ||
       ( ( gm.entity==-GEOMETRY_TRIANGLE ||
           gm.entity==-GEOMETRY_QUADRILATERAL ) && ndim<2 ) ||
       ( ( gm.entity==-GEOMETRY_TETRAHEDRAL ||
           gm.entity==-GEOMETRY_HEXAHEDRAL ) && ndim<3 ) ) {
    pri( "Error: geometry_moving needs -geometry_point, -geometry_triangle, "
         "-geometry_quadrilateral (2D, 3D), -geometry_hexahedral or "
         "-geometry_tetrahedral (3D)." );
    db_error( GEOMETRY_MOVING, index );
  }

  if ( !db_active_index( GEOMETRY_MOVING_PARAMETER, index, VERSION_NORMAL ) ) {
    pri( "Error: geometry_moving_parameter is missing." );
    db_error( GEOMETRY_MOVING_PARAMETER, index );
  }
  db( GEOMETRY_MOVING_PARAMETER, index, idum, parameter, length,
    VERSION_NORMAL, GET );
  if ( length!=gm.nvertex*ndim+( gm.entity==-GEOMETRY_POINT ? 1 : 0 ) ) {
    pri( "Error: geometry_moving_parameter needs the coordinates of every "
         "vertex (plus the radius for -geometry_point)." );
    db_error( GEOMETRY_MOVING_PARAMETER, index );
  }

  if ( !db_active_index( GEOMETRY_MOVING_OPERAT, index, VERSION_NORMAL ) ) {
    pri( "Error: geometry_moving_operat is missing." );
    db_error( GEOMETRY_MOVING_OPERAT, index );
  }
  db( GEOMETRY_MOVING_OPERAT, index, &operat, ddum, ldum, VERSION_NORMAL,
    GET );
  if ( operat!=-TRANSLATE ) {
    pri( "Error: geometry_moving_operat only allows -translate." );
    db_error( GEOMETRY_MOVING_OPERAT, index );
  }

  if ( !db_active_index( GEOMETRY_MOVING_OPERAT_PARAMETER, index,
      VERSION_NORMAL ) ) {
    pri( "Error: geometry_moving_operat_parameter is missing." );
    db_error( GEOMETRY_MOVING_OPERAT_PARAMETER, index );
  }
  db( GEOMETRY_MOVING_OPERAT_PARAMETER, index, idum, velocity, length,
    VERSION_NORMAL, GET );
  if ( length!=ndim ) {
    pri( "Error: geometry_moving_operat_parameter of -translate needs one "
         "velocity per space direction." );
    db_error( GEOMETRY_MOVING_OPERAT_PARAMETER, index );
  }
  array_set( gm.velocity, 0., MDIM );
  array_move( velocity, gm.velocity, ndim );

  if ( !db_active_index( GEOMETRY_MOVING_OPERAT_TIME, index,
      VERSION_NORMAL ) ) {
    pri( "Error: geometry_moving_operat_time is missing." );
    db_error( GEOMETRY_MOVING_OPERAT_TIME, index );
  }
  db( GEOMETRY_MOVING_OPERAT_TIME, index, idum, time, ldum, VERSION_NORMAL,
    GET );
  if ( time[1]<time[0] ) {
    pri( "Error: geometry_moving_operat_time needs end_time >= start_time." );
    db_error( GEOMETRY_MOVING_OPERAT_TIME, index );
  }
  gm.t_start = time[0];
  gm.t_end = time[1];

  n[0] = GM_NTIME_DEFAULT;
  n[1] = GM_NSPACE_DEFAULT;
  db( GEOMETRY_MOVING_N, index, n, ddum, ldum, VERSION_NORMAL,
    GET_IF_EXISTS );
  if ( n[0]<1 || n[1]<1 ) {
    pri( "Error: geometry_moving_n needs ntime >= 1 and nspace >= 1." );
    db_error( GEOMETRY_MOVING_N, index );
  }
  gm.ntime = n[0];
  gm.nspace = n[1];

    // bounding box of the entity at its start position (prefilter only)
  for ( idim=0; idim<ndim; idim++ ) {
    gm.box_min[idim] = gm.box_max[idim] = parameter[idim];
    for ( ivertex=1; ivertex<gm.nvertex; ivertex++ ) {
      gm.box_min[idim] = scalar_dmin( gm.box_min[idim],
        parameter[ivertex*ndim+idim] );
      gm.box_max[idim] = scalar_dmax( gm.box_max[idim],
        parameter[ivertex*ndim+idim] );
    }
  }
  if ( gm.entity==-GEOMETRY_POINT ) radius = scalar_dabs( parameter[ndim] );
  margin = radius;
  for ( idim=0; idim<ndim; idim++ )
    margin = scalar_dmax( margin, gm.box_max[idim]-gm.box_min[idim] );
  margin = radius + GM_BOX_MARGIN*margin + EPS_COORD;
  for ( idim=0; idim<ndim; idim++ ) {
    gm.box_min[idim] -= margin;
    gm.box_max[idim] += margin;
  }
}

// the scratch geometry record (entity at its start position) that
// geometry() tests the sample points against
static void gm_put_scratch( gm_entity_data &gm, double parameter[],
  long int scratch )
{
  long int length=0, idum[1];
  double record[8*MDIM+1];

  length = gm.nvertex*ndim;
  array_move( parameter, record, length );
  if ( gm.entity==-GEOMETRY_POINT ) {
    record[ndim] = parameter[ndim];
    length = ndim+1;
    db( GEOMETRY_POINT, scratch, idum, record, length, VERSION_NORMAL, PUT );
  }
  else if ( gm.entity==-GEOMETRY_TRIANGLE ) {
      // a filled triangle: the in-plane distance must vanish
    record[length] = EPS_COORD;
    length++;
    db( GEOMETRY_TRIANGLE, scratch, idum, record, length, VERSION_NORMAL,
      PUT );
  }
  else if ( gm.entity==-GEOMETRY_QUADRILATERAL ) {
    record[length] = EPS_COORD;
    length++;
    db( GEOMETRY_QUADRILATERAL, scratch, idum, record, length,
      VERSION_NORMAL, PUT );
  }
  else if ( gm.entity==-GEOMETRY_TETRAHEDRAL )
    db( GEOMETRY_TETRAHEDRAL, scratch, idum, record, length, VERSION_NORMAL,
      PUT );
  else
    db( GEOMETRY_HEXAHEDRAL, scratch, idum, record, length, VERSION_NORMAL,
      PUT );
}

static void gm_delete_scratch( gm_entity_data &gm, long int scratch )
{
  if ( gm.entity==-GEOMETRY_POINT )
    db_delete_index( GEOMETRY_POINT, scratch, VERSION_NORMAL );
  else if ( gm.entity==-GEOMETRY_TRIANGLE )
    db_delete_index( GEOMETRY_TRIANGLE, scratch, VERSION_NORMAL );
  else if ( gm.entity==-GEOMETRY_QUADRILATERAL )
    db_delete_index( GEOMETRY_QUADRILATERAL, scratch, VERSION_NORMAL );
  else if ( gm.entity==-GEOMETRY_TETRAHEDRAL )
    db_delete_index( GEOMETRY_TETRAHEDRAL, scratch, VERSION_NORMAL );
  else
    db_delete_index( GEOMETRY_HEXAHEDRAL, scratch, VERSION_NORMAL );
}

// a geometry index no geometry_* record uses, so that no per-index
// option (eps_iso, node type, filters) of a user geometry applies
static long int gm_scratch_index( void )
{
  long int idat=0, max=0, scratch=0;

  for ( idat=0; idat<MDAT; idat++ ) {
    if ( !strncmp( db_name(idat), "geometry_", 9 ) ) {
      db_max_index( idat, max, VERSION_NORMAL, GET );
      if ( max+1>scratch ) scratch = max+1;
    }
  }
  return scratch;
}

static void gm_node_coord( long int inod, double coord[] )
{
  long int ldum=0, idum[1];

  array_set( coord, 0., MDIM );
  if ( db_active_index( NODE_START_REFINED, inod, VERSION_NORMAL ) )
    db( NODE_START_REFINED, inod, idum, coord, ldum, VERSION_NORMAL, GET );
  else
    db( NODE, inod, idum, coord, ldum, VERSION_NORMAL, GET );
}

// Jacobian measure sqrt(det(J^T J)) of the ndim x pdim matrix jac
static double gm_measure( double jac[MDIM][MDIM], long int pdim )
{
  long int a=0, b=0, idim=0;
  double g[MDIM][MDIM], det=0.;

  for ( a=0; a<pdim; a++ ) {
    for ( b=0; b<pdim; b++ ) {
      g[a][b] = 0.;
      for ( idim=0; idim<ndim; idim++ ) g[a][b] += jac[idim][a]*jac[idim][b];
    }
  }
  if ( pdim==1 ) det = g[0][0];
  else if ( pdim==2 ) det = g[0][0]*g[1][1]-g[0][1]*g[1][0];
  else
    det = g[0][0]*(g[1][1]*g[2][2]-g[1][2]*g[2][1]) -
          g[0][1]*(g[1][0]*g[2][2]-g[1][2]*g[2][0]) +
          g[0][2]*(g[1][0]*g[2][1]-g[1][1]*g[2][0]);
  return det>0. ? sqrt(det) : 0.;
}

// sample points (coord, ndim per point) and weights of one element;
// returns the number of points (the arrays hold nspace^3 or MNOL points)
static long int gm_element_samples( long int name, long int nnol,
  long int nodes[], long int nspace, double coord[], double weight[] )
{
  long int npoint=0, pdim=0, npol=0, i=0, j=0, k=0, inol=0, idim=0, d=0,
    n[MDIM], ip[MDIM], corner[4];
  double h[MDIM][MNOL], p[MDIM][MNOL], iso=0., u=0., v=0., w=0.,
    L[4], fac=0., node_coord[MNOL][MDIM], jac[MDIM][MDIM];

  for ( inol=0; inol<nnol; inol++ )
    gm_node_coord( nodes[inol], node_coord[inol] );

    // tensor-product isoparametric elements: pdim directions, npol nodes
    // per direction, nodes numbered with the first direction fastest
  if ( ( name==-BAR2 && nnol==2 ) || ( name==-TRUSS && nnol==2 ) ||
       ( name==-BEAM && nnol==2 ) || ( name==-TRUSSBEAM && nnol==2 ) ) {
    pdim = 1; npol = 2;
  }
  else if ( name==-BAR3 && nnol==3 ) { pdim = 1; npol = 3; }
  else if ( name==-BAR4 && nnol==4 ) { pdim = 1; npol = 4; }
  else if ( name==-QUAD4 && nnol==4 ) { pdim = 2; npol = 2; }
  else if ( name==-QUAD9 && nnol==9 ) { pdim = 2; npol = 3; }
  else if ( name==-QUAD16 && nnol==16 ) { pdim = 2; npol = 4; }
  else if ( name==-HEX8 && nnol==8 ) { pdim = 3; npol = 2; }
  else if ( name==-HEX27 && nnol==27 ) { pdim = 3; npol = 3; }
  if ( pdim>ndim ) pdim = 0;
  if ( pdim>0 ) {
    for ( d=0; d<MDIM; d++ ) n[d] = ( d<pdim ? nspace : 1 );
    for ( k=0; k<n[2]; k++ ) {
      for ( j=0; j<n[1]; j++ ) {
        for ( i=0; i<n[0]; i++ ) {
          ip[0] = i; ip[1] = j; ip[2] = k;
          for ( d=0; d<pdim; d++ ) {
            iso = -1. + ( 2.*ip[d] + 1. ) / nspace;
            interpolation_polynomial( iso, npol, h[d], p[d] );
          }
          array_set( &coord[npoint*MDIM], 0., MDIM );
          for ( idim=0; idim<MDIM; idim++ )
            for ( d=0; d<MDIM; d++ ) jac[idim][d] = 0.;
          for ( inol=0; inol<nnol; inol++ ) {
            long int id[MDIM], rest=inol;
            double shape=1., dshape[MDIM];
            for ( d=0; d<pdim; d++ ) { id[d] = rest%npol; rest /= npol; }
            for ( d=0; d<pdim; d++ ) shape *= h[d][id[d]];
            for ( d=0; d<pdim; d++ ) {
              dshape[d] = 1.;
              for ( long int e=0; e<pdim; e++ )
                dshape[d] *= ( e==d ? p[e][id[e]] : h[e][id[e]] );
            }
            for ( idim=0; idim<ndim; idim++ ) {
              coord[npoint*MDIM+idim] += shape*node_coord[inol][idim];
              for ( d=0; d<pdim; d++ )
                jac[idim][d] += dshape[d]*node_coord[inol][idim];
            }
          }
          weight[npoint] = gm_measure( jac, pdim );
          npoint++;
        }
      }
    }
    return npoint;
  }

    // simplices: collapsed square / cube cell midpoints on the corners
  pdim = 0;
  if ( name==-TRIA3 && nnol==3 ) {
    pdim = 2; corner[0] = 0; corner[1] = 1; corner[2] = 2;
  }
  else if ( name==-TRIA6 && nnol==6 ) {
    pdim = 2; corner[0] = 0; corner[1] = 2; corner[2] = 5;
  }
  else if ( name==-TET4 && nnol==4 ) {
    pdim = 3; corner[0] = 0; corner[1] = 1; corner[2] = 2; corner[3] = 3;
  }
  else if ( name==-TET10 && nnol==10 ) {
    pdim = 3; corner[0] = 0; corner[1] = 2; corner[2] = 5; corner[3] = 9;
  }
  if ( pdim>ndim ) pdim = 0;
  if ( pdim>0 ) {
    for ( k=0; k<( pdim==3 ? nspace : 1 ); k++ ) {
      for ( j=0; j<nspace; j++ ) {
        for ( i=0; i<nspace; i++ ) {
          u = ( i + 0.5 ) / nspace;
          v = ( j + 0.5 ) / nspace;
          w = ( k + 0.5 ) / nspace;
          if ( pdim==2 ) {
            L[1] = u*(1.-v); L[2] = v; L[0] = 1.-L[1]-L[2];
            fac = 1.-v;
          }
          else {
            L[3] = w; L[2] = v*(1.-w); L[1] = u*(1.-v)*(1.-w);
            L[0] = 1.-L[1]-L[2]-L[3];
            fac = (1.-v)*(1.-w)*(1.-w);
          }
          array_set( &coord[npoint*MDIM], 0., MDIM );
          for ( inol=0; inol<=pdim; inol++ )
            for ( idim=0; idim<ndim; idim++ )
              coord[npoint*MDIM+idim] += L[inol]*node_coord[corner[inol]][idim];
          weight[npoint] = fac;
          npoint++;
        }
      }
    }
    return npoint;
  }

    // any other element: its nodes (the Professional rule)
  for ( inol=0; inol<nnol; inol++ ) {
    array_set( &coord[inol*MDIM], 0., MDIM );
    array_move( node_coord[inol], &coord[inol*MDIM], ndim );
    weight[inol] = 1.;
  }
  return nnol;
}

// part, entry and leave times of one element for one entity; returns 1
// when a sample point is ever inside
static long int gm_element_entity( gm_entity_data &gm, long int scratch,
  long int npoint, double coord[], double weight[], double &part,
  double in_time[], double out_time[], long int &any_out )
{
  long int ipoint=0, itime=0, idim=0, overlap=0, in_geometry=0,
    any_in=0, first_in=0, was_in=0, out_index=0, geometry_entity[2];
  double t=0., dt=0., total=0., inside=0., rdum=0., entry=0., leave=0.,
    shift[MDIM], point[MDIM], ddum[MDIM], box_min[MDIM], box_max[MDIM];

  geometry_entity[0] = gm.entity;
  geometry_entity[1] = scratch;
  dt = ( gm.t_end - gm.t_start ) / gm.ntime;

  for ( idim=0; idim<ndim; idim++ ) {
    box_min[idim] = box_max[idim] = coord[idim];
    for ( ipoint=1; ipoint<npoint; ipoint++ ) {
      box_min[idim] = scalar_dmin( box_min[idim], coord[ipoint*MDIM+idim] );
      box_max[idim] = scalar_dmax( box_max[idim], coord[ipoint*MDIM+idim] );
    }
  }

  long int *first = get_new_int( npoint );
  long int *now = get_new_int( npoint );
  long int *leave_index = get_new_int( npoint );
  array_set( first, -1, npoint );
  array_set( now, 0, npoint );
  array_set( leave_index, -1, npoint );

  for ( itime=0; itime<=gm.ntime; itime++ ) {
    t = gm.t_start + itime*dt;
    for ( idim=0; idim<ndim; idim++ )
      shift[idim] = gm.velocity[idim]*( t - gm.t_start );
    overlap = 1;
    for ( idim=0; idim<ndim; idim++ ) {
      if ( box_max[idim]-shift[idim]<gm.box_min[idim] ||
           box_min[idim]-shift[idim]>gm.box_max[idim] ) overlap = 0;
    }
    for ( ipoint=0; ipoint<npoint; ipoint++ ) {
      was_in = now[ipoint];
      in_geometry = 0;
      if ( overlap ) {
        array_set( point, 0., MDIM );
        for ( idim=0; idim<ndim; idim++ )
          point[idim] = coord[ipoint*MDIM+idim] - shift[idim];
        geometry( -1, point, geometry_entity, in_geometry, rdum, ddum, rdum,
          ddum, NODE_START_REFINED, PROJECT_INSIDE, VERSION_NORMAL );
      }
      now[ipoint] = in_geometry ? 1 : 0;
      if ( now[ipoint] ) {
        if ( first[ipoint]<0 ) first[ipoint] = itime;
      }
      else if ( was_in )
        leave_index[ipoint] = itime-1;
    }
  }

  part = 0.;
  any_out = 0;
  for ( ipoint=0; ipoint<npoint; ipoint++ ) {
    total += weight[ipoint];
    if ( first[ipoint]<0 ) continue;
    inside += weight[ipoint];
    first_in = first[ipoint]>0 ? first[ipoint]-1 : 0;
    entry = gm.t_start + first_in*dt;
    if ( !any_in ) { in_time[0] = in_time[1] = entry; }
    in_time[0] = scalar_dmin( in_time[0], entry );
    in_time[1] = scalar_dmax( in_time[1], entry );
    any_in = 1;
    out_index = leave_index[ipoint];
    if ( out_index>=0 ) {
      leave = gm.t_start + out_index*dt;
      if ( !any_out ) { out_time[0] = out_time[1] = leave; }
      out_time[0] = scalar_dmin( out_time[0], leave );
      out_time[1] = scalar_dmax( out_time[1], leave );
      any_out = 1;
    }
  }
  if ( any_in && total>0. ) part = inside / total;

  delete[] first;
  delete[] now;
  delete[] leave_index;

  return any_in && part>0.;
}

void geometry_moving_initialise( long int first_of_control )
{
  long int icontrol=0, control_value=0, index=0, max_index=0, ngm=0,
    igm=0, element=0, max_element=0, nnol=0, npoint=0, length=0,
    scratch=0, any_out=0, best_out=0, best_index=0, swit=0,
    ldum=0, idum[1], el[1+MNOL];
  double part=0., best_part=0., ddum[1], in_time[2], out_time[2],
    best_in[2], best_out_time[2], parameter[8*MDIM+1];

  if ( !first_of_control ) return;
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  if ( !db_active_index( CONTROL_GEOMETRY_MOVING, icontrol, VERSION_NORMAL ) )
    return;
  db( CONTROL_GEOMETRY_MOVING, icontrol, &control_value, ddum, ldum,
    VERSION_NORMAL, GET );
  if ( control_value!=-INITIALISE ) {
    pri( "Error: control_geometry_moving only allows -initialise." );
    db_error( CONTROL_GEOMETRY_MOVING, icontrol );
  }

  swit = set_swit(-1,-1,"geometry_moving_initialise");
  if ( swit ) pri( "In routine GEOMETRY_MOVING_INITIALISE" );

  db_max_index( GEOMETRY_MOVING, max_index, VERSION_NORMAL, GET );
  for ( index=0; index<=max_index; index++ )
    if ( db_active_index( GEOMETRY_MOVING, index, VERSION_NORMAL ) ) ngm++;

  gm_entity_data *gm = new gm_entity_data[ngm>0 ? ngm : 1];
  double *gm_parameter = get_new_dbl( (ngm>0 ? ngm : 1)*(8*MDIM+1) );
  igm = 0;
  for ( index=0; index<=max_index; index++ ) {
    if ( db_active_index( GEOMETRY_MOVING, index, VERSION_NORMAL ) ) {
      gm_read( index, gm[igm], parameter );
      array_move( parameter, &gm_parameter[igm*(8*MDIM+1)], 8*MDIM+1 );
      igm++;
    }
  }

    // the previous initialisation is replaced
  db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
  for ( element=0; element<=max_element; element++ ) {
    db_delete_index( ELEMENT_GEOMETRY_MOVING, element, VERSION_NORMAL );
    db_delete_index( ELEMENT_GEOMETRY_MOVING_PART, element, VERSION_NORMAL );
    db_delete_index( ELEMENT_GEOMETRY_MOVING_IN_TIME, element,
      VERSION_NORMAL );
    db_delete_index( ELEMENT_GEOMETRY_MOVING_OUT_TIME, element,
      VERSION_NORMAL );
  }

  long int max_points = MNOL;
  for ( igm=0; igm<ngm; igm++ ) {
    long int np = gm[igm].nspace*gm[igm].nspace*gm[igm].nspace;
    if ( np>max_points ) max_points = np;
  }
  double *coord = get_new_dbl( max_points*MDIM );
  double *weight = get_new_dbl( max_points );

  scratch = gm_scratch_index();
  for ( igm=0; igm<ngm; igm++ ) {
    gm_put_scratch( gm[igm], &gm_parameter[igm*(8*MDIM+1)], scratch );
    for ( element=0; element<=max_element; element++ ) {
      if ( !db_active_index( ELEMENT, element, VERSION_NORMAL ) ) continue;
      db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
      nnol = length - 1;
      npoint = gm_element_samples( el[0], nnol, &el[1], gm[igm].nspace,
        coord, weight );
      if ( !gm_element_entity( gm[igm], scratch, npoint, coord, weight,
          part, in_time, out_time, any_out ) ) continue;
        // several entities: the one with the largest part (the first
        // index on a tie)
      best_part = 0.;
      db( ELEMENT_GEOMETRY_MOVING_PART, element, idum, &best_part, ldum,
        VERSION_NORMAL, GET_IF_EXISTS );
      if ( part<=best_part ) continue;
      best_index = gm[igm].index;
      best_in[0] = in_time[0]; best_in[1] = in_time[1];
      best_out = any_out;
      best_out_time[0] = out_time[0]; best_out_time[1] = out_time[1];
      length = 1;
      db( ELEMENT_GEOMETRY_MOVING, element, &best_index, ddum, length,
        VERSION_NORMAL, PUT );
      db( ELEMENT_GEOMETRY_MOVING_PART, element, idum, &part, length,
        VERSION_NORMAL, PUT );
      length = 2;
      db( ELEMENT_GEOMETRY_MOVING_IN_TIME, element, idum, best_in, length,
        VERSION_NORMAL, PUT );
      if ( best_out )
        db( ELEMENT_GEOMETRY_MOVING_OUT_TIME, element, idum, best_out_time,
          length, VERSION_NORMAL, PUT );
      else
        db_delete_index( ELEMENT_GEOMETRY_MOVING_OUT_TIME, element,
          VERSION_NORMAL );
    }
    gm_delete_scratch( gm[igm], scratch );
  }

  delete[] coord;
  delete[] weight;
  delete[] gm_parameter;
  delete[] gm;

  if ( swit ) pri( "Out routine GEOMETRY_MOVING_INITIALISE" );
}

void geometry_moving_apply( double time_current )
{
  long int element=0, max_element=0, max_element_all=0,
    any_element_deleted=0, length=0, idum[1];
  double part=0., factor=0., ramp=0., in_time[2];
  long int ldum=0;

  if ( db_max_index( ELEMENT_GEOMETRY_MOVING_PART, max_element,
      VERSION_NORMAL, GET )<0 ) return;
  for ( element=0; element<=max_element; element++ ) {
    if ( !db_active_index( ELEMENT_GEOMETRY_MOVING_PART, element,
        VERSION_NORMAL ) ) continue;
    if ( !db_active_index( ELEMENT, element, VERSION_NORMAL ) ) continue;
    db( ELEMENT_GEOMETRY_MOVING_PART, element, idum, &part, ldum,
      VERSION_NORMAL, GET );
    db( ELEMENT_GEOMETRY_MOVING_IN_TIME, element, idum, in_time, ldum,
      VERSION_NORMAL, GET );
    if ( time_current<=in_time[0]+GM_EPS_TIME ) continue;
    if ( in_time[1]-in_time[0]>GM_EPS_TIME )
      ramp = scalar_dmin( 1., ( time_current - in_time[0] ) /
        ( in_time[1] - in_time[0] ) );
    else
      ramp = 1.;
    factor = 1. - part*ramp;
    if ( factor>GM_EPS_FACTOR ) {
        // elem() stores the frozen released load (ELEMENT_RHSIDE_DELETE)
        // in the parallel element loop: allocate its range here (a PUT
        // beyond the allocated range aborts with "Data is allocated in a
        // parallel loop"; the Professional 02-08-2026 stops on it with
        // threads for its own geometry_moving1 example)
      db_max_index( ELEMENT, max_element_all, VERSION_NORMAL, GET );
      db_allocate( ELEMENT_RHSIDE_DELETE, max_element_all, VERSION_NORMAL,
        -MINIMAL );
      length = 1;
      db( ELEMENT_DELETE_FACTOR, element, idum, &factor, length,
        VERSION_NORMAL, PUT );
    }
    else {
      delete_element( element, VERSION_NORMAL );
      any_element_deleted = 1;
    }
  }
  if ( any_element_deleted ) mesh_has_changed( VERSION_NORMAL );
}
