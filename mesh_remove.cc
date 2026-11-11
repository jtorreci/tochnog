/*
    Copyright (C) 2026  Tochnog GNU contributors

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

// control_mesh_remove (manual Professional 6.229-6.234), evaluated on every
// step of a control index that has the record (step_start, top.cc).
//
//  -method1: an element of element_group_0 is removed when all its nodes
//    are nodes of one element of element_group_1, element_group_2, ...
//    (node numbers, not positions: BACKLOG C126). An element removed earlier
//    in the same pass is no container any more.
//  control_mesh_remove_geometry (6.230): only elements with every node in
//    the geometry are removed.
//  control_mesh_remove_really (6.232, default -yes): -yes deletes the
//    removed elements; -no makes them inactive: they become active again
//    on a later step whose control_mesh_remove no longer removes them.
//  control_mesh_remove_keep_geometry (6.231): an inactive element with
//    every node in this geometry of the current index stays removed.
//  control_mesh_remove_really_activate_factor (6.234, default 1) and
//    _activate_all (6.233, default -no): the dofs of the nodes of a
//    reactivated element, see mesh_remove_activate.
//
// GNU representation of an inactive element (DIVERGENCES REMOVE-INACTIVE-
// RECORDS): the element is deleted from the mesh like a really removed one
// and kept in the records element_remove_inactive (name and nodes),
// element_remove_inactive_group and element_remove_inactive_node (the
// coordinates of its nodes). Its nodes that no active element uses are
// deleted by the mesh change (nod_nod) and created again on the
// reactivation. The Professional 02-08-2026 keeps the element in the mesh
// (element_remove_active -yes, element_empty -empty) and its unused nodes
// fully bounded.
//
// Measured on the Professional 02-08-2026 (F5 T4 U9e probes): the default of
// _really is -yes; a step of an index without control_mesh_remove changes
// nothing; the factor and _activate_all of the reactivating index apply; the
// _keep_geometry of the current index applies (an element is kept when all
// its nodes are in the geometry; the _really of the index then applies, so
// with the default -yes the kept element is deleted); a reactivated element
// is "removed once" no longer.

static long int mesh_remove_dof_scaled( long int iuknwn )

{
  // the stresses and the pore pressure (manual 6.234)
  if ( stres_indx>=0 && iuknwn>=stres_indx && iuknwn<stres_indx+6*nder )
    return 1;
  if ( pres_indx>=0 && iuknwn>=pres_indx && iuknwn<pres_indx+nder )
    return 1;
  return 0;
}

// Every node of the element (node numbers el[1..nnol], node_start_refined
// coordinates start_refined[] for the nodes that no longer exist) in the
// geometry: the membership test of element_in_geometry -all with
// PROJECT_EXACT (geometry.cc), used for the inactive elements.
static long int mesh_remove_nodes_in_geometry( long int el[], long int nnol,
  double start_refined[], long int geometry_entity[] )

{
  long int inol=0, found=0;
  double rdum=0., ddum[MDIM], coord[MDIM];

  for ( inol=0; inol<nnol; inol++ ) {
    long int inod = el[1+inol];
    if ( db_active_index( NODE, inod, VERSION_NORMAL ) )
      geometry( inod, ddum, geometry_entity, found, rdum, ddum, rdum, ddum,
        NODE_START_REFINED, PROJECT_EXACT, VERSION_NORMAL );
    else {
      array_move( &start_refined[inol*ndim], coord, ndim );
      geometry( -1, coord, geometry_entity, found, rdum, ddum, rdum, ddum,
        NODE_START_REFINED, PROJECT_EXACT, VERSION_NORMAL );
    }
    if ( !found ) return 0;
  }
  return nnol>0;
}

// Make an active element inactive: keep it in the element_remove_inactive
// records and delete it from the mesh.
static void mesh_remove_deactivate( long int element )

{
  long int length=0, nnol=0, inol=0, ldum=0, one=1, group=0,
    idum[1], el[1+MNOL];
  double ddum[1], node_coord[2*MNOL*MDIM];

  db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
  nnol = length - 1;
  long int has_group = db( ELEMENT_GROUP, element, &group, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  for ( inol=0; inol<nnol; inol++ ) {
    db( NODE, el[1+inol], idum, &node_coord[inol*ndim], ldum,
      VERSION_NORMAL, GET );
    db( NODE_START_REFINED, el[1+inol], idum, &node_coord[(nnol+inol)*ndim],
      ldum, VERSION_NORMAL, GET );
  }
  delete_element( element, VERSION_NORMAL );
  db( ELEMENT_REMOVE_INACTIVE, element, el, ddum, length, VERSION_NORMAL,
    PUT );
  if ( has_group )
    db( ELEMENT_REMOVE_INACTIVE_GROUP, element, &group, ddum, one,
      VERSION_NORMAL, PUT );
  ldum = 2*nnol*ndim;
  db( ELEMENT_REMOVE_INACTIVE_NODE, element, idum, node_coord, ldum,
    VERSION_NORMAL, PUT );
}

static void mesh_remove_forget( long int element )

{
  db_delete_index( ELEMENT_REMOVE_INACTIVE, element, VERSION_NORMAL );
  db_delete_index( ELEMENT_REMOVE_INACTIVE_GROUP, element, VERSION_NORMAL );
  db_delete_index( ELEMENT_REMOVE_INACTIVE_NODE, element, VERSION_NORMAL );
}

// Reactivate an inactive element (manual 6.234, Professional 02-08-2026):
//  - the active nodes are those of the elements active at the start of the
//    pass (node_active); the average is taken over the active nodes of this
//    element, of every node dof (node_dof and node_dof_start_refined);
//    without active nodes the average is 0;
//  - a node that is not active gets the average: the stresses and the pore
//    pressure times the factor, the other dofs (velocity, displacement,
//    temperature, ...) unscaled; a node that does not exist any more is
//    created again with its kept coordinates;
//  - with _activate_all -yes the active nodes of the element get the
//    factor times the average for the stresses and the pore pressure too
//    (manual 6.233; their other dofs are left as they are: DIVERGENCES
//    REMOVE-ACTIVATE-ALL);
//  - the elements are reactivated in element order and a node gets its
//    dofs from the first of them only (assigned).
static void mesh_remove_activate( long int element, long int node_active[],
  long int assigned[], double factor, long int activate_all )

{
  long int length=0, nnol=0, inol=0, iuknwn=0, nactive=0, ldum=0, one=1,
    group=0, idum[1], el[1+MNOL];
  double ddum[1], node_coord[2*MNOL*MDIM], *avg=NULL,
    *avg_start_refined=NULL, *dof=NULL;

  if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) ) {
    cout << "Error: control_mesh_remove cannot activate element " << element
         << " again: another element has its number.\n";
    exit(TN_EXIT_STATUS);
  }
  db( ELEMENT_REMOVE_INACTIVE, element, el, ddum, length, VERSION_NORMAL,
    GET );
  nnol = length - 1;
  long int has_group = db( ELEMENT_REMOVE_INACTIVE_GROUP, element, &group,
    ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( ELEMENT_REMOVE_INACTIVE_NODE, element, idum, node_coord, ldum,
    VERSION_NORMAL, GET );

  if ( nuknwn>0 ) {
    avg = get_new_dbl( nuknwn );
    avg_start_refined = get_new_dbl( nuknwn );
    dof = get_new_dbl( nuknwn );
    array_set( avg, 0., nuknwn );
    array_set( avg_start_refined, 0., nuknwn );
    for ( inol=0; inol<nnol; inol++ ) {
      long int inod = el[1+inol];
      if ( !node_active[inod] ) continue;
      array_add( db_dbl( NODE_DOF, inod, VERSION_NORMAL ), avg, avg,
        nuknwn );
      array_add( db_dbl( NODE_DOF_START_REFINED, inod, VERSION_NORMAL ),
        avg_start_refined, avg_start_refined, nuknwn );
      nactive++;
    }
    if ( nactive>0 ) {
      array_multiply( avg, avg, 1./nactive, nuknwn );
      array_multiply( avg_start_refined, avg_start_refined, 1./nactive,
        nuknwn );
    }
  }

  for ( inol=0; inol<nnol; inol++ ) {
    long int inod = el[1+inol];
    if ( assigned[inod] ) continue;
    if ( node_active[inod] ) {
      if ( activate_all!=-YES ) continue;
      if ( nuknwn>0 ) {
        array_move( db_dbl( NODE_DOF, inod, VERSION_NORMAL ), dof, nuknwn );
        for ( iuknwn=0; iuknwn<nuknwn; iuknwn++ )
          if ( mesh_remove_dof_scaled( iuknwn ) )
            dof[iuknwn] = factor * avg[iuknwn];
        db( NODE_DOF, inod, idum, dof, nuknwn, VERSION_NORMAL, PUT );
      }
    }
    else {
      if ( nuknwn>0 ) {
        for ( iuknwn=0; iuknwn<nuknwn; iuknwn++ ) {
          dof[iuknwn] = avg[iuknwn];
          if ( mesh_remove_dof_scaled( iuknwn ) ) dof[iuknwn] *= factor;
        }
      }
      if ( !db_active_index( NODE, inod, VERSION_NORMAL ) ) {
        // the node again with the records a node keeps over a mesh change
        // (as mesh_map.cc; create_node copies through db_dbl pointers that
        // a PUT beyond the highest index reallocates)
        db( NODE, inod, idum, &node_coord[inol*ndim], ndim, VERSION_NORMAL,
          PUT );
        db( NODE_START_REFINED, inod, idum, &node_coord[(nnol+inol)*ndim],
          ndim, VERSION_NORMAL, PUT );
        if ( nuknwn>0 )
          db( NODE_DOF_START_REFINED, inod, idum, avg_start_refined, nuknwn,
            VERSION_NORMAL, PUT );
      }
      if ( nuknwn>0 )
        db( NODE_DOF, inod, idum, dof, nuknwn, VERSION_NORMAL, PUT );
    }
    assigned[inod] = 1;
  }

  // the element again, without its old element data (as a remeshed element:
  // integration point dofs initialised from the nodes)
  mesh_remove_forget( element );
  create_element( element, element, el, length, VERSION_NORMAL,
    VERSION_NORMAL );
  if ( has_group )
    db( ELEMENT_GROUP, element, &group, ddum, one, VERSION_NORMAL, PUT );
  if ( options_element_dof==-YES ) {
    long int initialised = 1;
    db( ELEMENT_DOF_INITIALISED, element, &initialised, ddum, one,
      VERSION_NORMAL, PUT );
  }

  if ( avg!=NULL ) delete[] avg;
  if ( avg_start_refined!=NULL ) delete[] avg_start_refined;
  if ( dof!=NULL ) delete[] dof;
}

void mesh_remove( long int icontrol )

{
  long int method=0, length=0, ldum=0, ielem=0, max_elem=0, max_inactive=0,
    max_node=0, jelem=0, inol=0, jnol=0, nnol=0, ngrp=0, i=0, inlist=0,
    elem_group=0, really=-YES, activate_all=-NO, use_geometry=0,
    use_keep=0, nchanged=0, list[DATA_ITEM_SIZE], remove_geometry[2],
    keep_geometry[2], idum[1], el[1+MNOL], elj[1+MNOL];
  long int *active=NULL, *inactive=NULL, *removed=NULL, *node_active=NULL,
    *assigned=NULL;
  double factor=1., ddum[1], node_coord[2*MNOL*MDIM];

  db( CONTROL_MESH_REMOVE, icontrol, list, ddum, length, VERSION_NORMAL, GET );
  if ( length<1 ) return;
  method = list[0];
  ngrp = length - 1;

  if ( db_active_index( CONTROL_MESH_REMOVE_GEOMETRY, icontrol,
       VERSION_NORMAL ) ) {
    db( CONTROL_MESH_REMOVE_GEOMETRY, icontrol, remove_geometry, ddum, ldum,
      VERSION_NORMAL, GET );
    use_geometry = 1;
  }
  if ( db_active_index( CONTROL_MESH_REMOVE_KEEP_GEOMETRY, icontrol,
       VERSION_NORMAL ) ) {
    db( CONTROL_MESH_REMOVE_KEEP_GEOMETRY, icontrol, keep_geometry, ddum,
      ldum, VERSION_NORMAL, GET );
    use_keep = 1;
  }
  db( CONTROL_MESH_REMOVE_REALLY, icontrol, &really, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  if ( really!=-YES && really!=-NO )
    db_error( CONTROL_MESH_REMOVE_REALLY, icontrol );
  db( CONTROL_MESH_REMOVE_REALLY_ACTIVATE_ALL, icontrol, &activate_all, ddum,
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  if ( activate_all!=-YES && activate_all!=-NO )
    db_error( CONTROL_MESH_REMOVE_REALLY_ACTIVATE_ALL, icontrol );
  db( CONTROL_MESH_REMOVE_REALLY_ACTIVATE_FACTOR, icontrol, idum, &factor,
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  if ( factor<0. || factor>1. ) {
    // the Professional 02-08-2026 message
    cout << "Error, CONTROL_MESH_REMOVE_REALLY_ACTIVATE_FACTOR factor should "
            "be between 0 and 1.\n";
    exit(TN_EXIT_STATUS);
  }

  db_max_index( ELEMENT, max_elem, VERSION_NORMAL, GET );
  db_max_index( ELEMENT_REMOVE_INACTIVE, max_inactive, VERSION_NORMAL, GET );
  if ( max_inactive>max_elem ) max_elem = max_inactive;
  if ( max_elem<0 ) return;

  if ( method!=-METHOD1 ) {
    cout << "Error: control_mesh_remove method must be -method1 or -method3.\n";
    exit(TN_EXIT_STATUS);
  }

  // the elements: active (in the mesh) or inactive (kept in the
  // element_remove_inactive records)
  active = get_new_int( 1+max_elem );
  inactive = get_new_int( 1+max_elem );
  removed = get_new_int( 1+max_elem );
  db_max_index( NODE, max_node, VERSION_NORMAL, GET );
  for ( ielem=0; ielem<=max_elem; ielem++ ) {
    active[ielem] = db_active_index( ELEMENT, ielem, VERSION_NORMAL );
    inactive[ielem] = !active[ielem] &&
      db_active_index( ELEMENT_REMOVE_INACTIVE, ielem, VERSION_NORMAL );
    removed[ielem] = 0;
    if ( inactive[ielem] ) {
      db( ELEMENT_REMOVE_INACTIVE, ielem, el, ddum, length, VERSION_NORMAL,
        GET );
      for ( inol=1; inol<length; inol++ )
        if ( el[inol]>max_node ) max_node = el[inol];
    }
  }

  // the active nodes before the pass (Professional node_active_last)
  node_active = get_new_int( 1+max_node );
  assigned = get_new_int( 1+max_node );
  array_set( node_active, 0, 1+max_node );
  array_set( assigned, 0, 1+max_node );
  for ( ielem=0; ielem<=max_elem; ielem++ ) {
    if ( !active[ielem] ) continue;
    db( ELEMENT, ielem, el, ddum, length, VERSION_NORMAL, GET );
    for ( inol=1; inol<length; inol++ ) node_active[el[inol]] = 1;
  }

  // which elements are removed in this pass
  for ( ielem=0; ielem<=max_elem && ngrp>=1; ielem++ ) {
    if ( !active[ielem] && !inactive[ielem] ) continue;
    elem_group = 0;
    if ( active[ielem] ) {
      db( ELEMENT, ielem, el, ddum, length, VERSION_NORMAL, GET );
      db( ELEMENT_GROUP, ielem, &elem_group, ddum, ldum, VERSION_NORMAL,
        GET_IF_EXISTS );
    }
    else {
      db( ELEMENT_REMOVE_INACTIVE, ielem, el, ddum, length, VERSION_NORMAL,
        GET );
      db( ELEMENT_REMOVE_INACTIVE_GROUP, ielem, &elem_group, ddum, ldum,
        VERSION_NORMAL, GET_IF_EXISTS );
      db( ELEMENT_REMOVE_INACTIVE_NODE, ielem, idum, node_coord, ldum,
        VERSION_NORMAL, GET );
    }
    nnol = length - 1;
    long int in_geometry = 1;
    // control_mesh_remove_geometry (6.230): only elements with every
    // node in the geometry (Professional 02-08-2026)
    if ( use_geometry ) {
      if ( active[ielem] )
        in_geometry = element_in_geometry( ielem, remove_geometry, -ALL,
          NODE_START_REFINED, PROJECT_EXACT, NULL, VERSION_NORMAL );
      else
        in_geometry = mesh_remove_nodes_in_geometry( el, nnol,
          &node_coord[nnol*ndim], remove_geometry );
    }
    if ( in_geometry && elem_group==list[1] ) {
      // all nodes of ielem are nodes of an active element of the groups
      // list[2..] that is not removed in this pass
      long int inside_other = 0;
      for ( jelem=0; jelem<=max_elem && !inside_other; jelem++ ) {
        if ( jelem==ielem || !active[jelem] || removed[jelem] ) continue;
        long int jgrp = 0;
        db( ELEMENT_GROUP, jelem, &jgrp, ddum, ldum, VERSION_NORMAL,
          GET_IF_EXISTS );
        inlist = 0;
        for ( i=1; i<ngrp; i++ ) {
          if ( jgrp==list[1+i] ) inlist = 1;
        }
        if ( !inlist ) continue;
        long int length_j = 0;
        db( ELEMENT, jelem, elj, ddum, length_j, VERSION_NORMAL, GET );
        long int all_in = 1;
        for ( inol=0; inol<nnol && all_in; inol++ ) {
          long int found = 0;
          for ( jnol=1; jnol<length_j; jnol++ )
            if ( el[1+inol]==elj[jnol] ) found = 1;
          if ( !found ) all_in = 0;
        }
        if ( all_in ) inside_other = 1;
      }
      if ( inside_other ) removed[ielem] = 1;
    }
    // control_mesh_remove_keep_geometry (6.231): an inactive element
    // stays removed when all its nodes are in the geometry
    if ( !removed[ielem] && inactive[ielem] && use_keep &&
         mesh_remove_nodes_in_geometry( el, nnol, &node_coord[nnol*ndim],
           keep_geometry ) )
      removed[ielem] = 1;
  }

  // remove: delete, or make inactive
  for ( ielem=0; ielem<=max_elem; ielem++ ) {
    if ( !removed[ielem] ) continue;
    if ( active[ielem] ) {
      if ( really==-YES ) delete_element( ielem, VERSION_NORMAL );
      else mesh_remove_deactivate( ielem );
      nchanged++;
    }
    else if ( really==-YES ) {
      mesh_remove_forget( ielem );
      nchanged++;
    }
  }

  // activate again the inactive elements that are no longer removed
  for ( ielem=0; ielem<=max_elem; ielem++ ) {
    if ( !inactive[ielem] || removed[ielem] ) continue;
    mesh_remove_activate( ielem, node_active, assigned, factor,
      activate_all );
    nchanged++;
  }

  delete[] active;
  delete[] inactive;
  delete[] removed;
  delete[] node_active;
  delete[] assigned;

  // as before U9e: the mesh change after every -method1 pass with groups
  if ( ngrp>=1 || nchanged>0 ) mesh_has_changed( VERSION_NORMAL );
}
