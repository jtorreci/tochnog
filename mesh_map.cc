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

// control_mesh_map index -yes (manual Professional 6.210, F5 T4 U9d).
//
// Each map toggles the quadratic volume elements between their quadratic
// and their linear form, so that a large linear calculation is first done
// with the linear elements (pardiso) and the quadratic calculation (bicg)
// starts from that solution:
//   map to linear: every -hex27 / -tet10 / -prism15 becomes the -hex8 /
//     -tet4 / -prism6 of its corner nodes; the quadratic element is kept
//     in the element_map record, the coordinates of its nodes in the GNU
//     record element_map_node. The mid nodes no element uses any more
//     are deleted by the mesh change (nod_nod, as after any GNU mesh
//     change); the Professional keeps them, unused, with their dofs
//     unchanged (DIVERGENCES MESH-MAP-MID-NODES).
//   map back: every linear element with an element_map record becomes its
//     quadratic element again; its mid nodes are created again
//     (coordinates from element_map_node) and their dofs
//     (node_dof, all of them) are interpolated from the corner nodes with
//     the linear shape functions at the natural coordinates of the mid
//     node (edge mid node: mean of the 2 edge corners, hex27 face centre:
//     mean of 4, hex27 body centre: mean of 8), = Professional.
// The integration point data of a mapped element are re-initialised from
// the nodal values (as create_element does after a remesh).
// Measured on the Professional 02-08-2026 (test/other/map1 and probes):
// mid nodes kept, element_map kept after the map back, each further map
// toggles again, a -hex20 maps to -hex8 (the GNU converts every -hex20 to
// -hex27 first, mesh_convert_hex20), the other elements (-hex8, -prism18,
// -prism12, 2D -quad9 / -tria6, ...) are left as they are. The
// Professional takes the wrong corners of the -tet10 (1 4 6 10) and the
// -prism15 (1 3 3 4 5 6); the GNU takes the corners (DIVERGENCES
// MESH-MAP-TET10, MESH-MAP-PRISM15).

// Quadratic element -> linear element: name, number of nodes and the
// positions (0-based) of the corner nodes in the quadratic node list, in
// the node order of the linear element. Returns 0 for other elements.
static long int mesh_map_linear( long int name, long int nnol,
  long int &lin_name, long int &lin_nnol, long int corner[] )

{
  if ( name==-HEX27 && nnol==27 ) {
    long int c[8] = { 0, 2, 6, 8, 18, 20, 24, 26 };
    lin_name = -HEX8; lin_nnol = 8;
    for ( long int i=0; i<8; i++ ) corner[i] = c[i];
    return 1;
  }
  else if ( name==-TET10 && nnol==10 ) {
    long int c[4] = { 0, 2, 5, 9 };
    lin_name = -TET4; lin_nnol = 4;
    for ( long int i=0; i<4; i++ ) corner[i] = c[i];
    return 1;
  }
  else if ( name==-PRISM15 && nnol==15 ) {
    lin_name = -PRISM6; lin_nnol = 6;
    for ( long int i=0; i<6; i++ ) corner[i] = i;
    return 1;
  }
  return 0;
}

// Weights of the linear element nodes (lin_nnol of them) at quadratic node
// inol: the linear shape functions at the natural coordinates of that
// node.
static void mesh_map_weights( long int name, long int inol, double w[] )

{
  for ( long int i=0; i<8; i++ ) w[i] = 0.;
  if ( name==-HEX27 ) {
    // lattice node (i,j,k), i fastest; corner c = (ci,cj,ck) in {0,2}
    long int ijk[3] = { inol%3, (inol/3)%3, inol/9 };
    for ( long int c=0; c<8; c++ ) {
      long int cijk[3] = { 2*(c%2), 2*((c/2)%2), 2*(c/4) };
      double f = 1.;
      for ( long int d=0; d<3; d++ ) {
        long int dist = ijk[d] - cijk[d];
        if ( dist<0 ) dist = -dist;
        f *= 1. - 0.5*dist;
      }
      w[c] = f;
    }
  }
  else if ( name==-TET10 ) {
    // tet10: 0 c1, 1 (c1,c2), 2 c2, 3 (c1,c3), 4 (c2,c3), 5 c3,
    // 6 (c1,c4), 7 (c2,c4), 8 (c3,c4), 9 c4 (c = tet4 corners 1..4)
    long int a[10] = { 0, 0, 1, 0, 1, 2, 0, 1, 2, 3 },
      b[10] = { 0, 1, 1, 2, 2, 2, 3, 3, 3, 3 };
    w[a[inol]] += 0.5; w[b[inol]] += 0.5;
  }
  else {
    assert( name==-PRISM15 );
    // prism15: 0-5 corners, 6-8 vertical edges (k,k+3), 9-11 base edges
    // (0,1) (1,2) (2,0), 12-14 top edges (3,4) (4,5) (5,3)
    long int a[15] = { 0, 1, 2, 3, 4, 5, 0, 1, 2, 0, 1, 2, 3, 4, 5 },
      b[15] = { 0, 1, 2, 3, 4, 5, 3, 4, 5, 1, 2, 0, 4, 5, 3 };
    w[a[inol]] += 0.5; w[b[inol]] += 0.5;
  }
}

// Give element its new node list el (length entries) in place. The data
// that belong to the old node list or integration points are dropped
// (element_dof, element_intpnt_*, element_matrix_*, the rhside / nonlocal
// records) and the integration point dofs are initialised again from the
// nodes, as create_element does for a remeshed element
// (element_intpnt_dof zero, element_dof_initialised 0). create_element
// itself is not used in place: it copies every element record through a
// DATA_ITEM_SIZE buffer, too short for the element_matrix_* records of a
// solved element. The other records (group, geometry, activity, ...) stay.
static void mesh_map_element_put( long int element, long int el[],
  long int length )

{
  long int idat=0, zero=0, one=1, mnolnuknwn=npointmax*nuknwn;
  double ddum[1];

  for ( idat=0; idat<MDAT; idat++ ) {
    if ( db_data_class(idat)!=ELEMENT || idat==ELEMENT ||
         idat==ELEMENT_MAP || idat==ELEMENT_MAP_NODE ) continue;
    if ( !db_active_index( idat, element, VERSION_NORMAL ) ) continue;
    const char *nm = db_name( idat );
    if ( idat==ELEMENT_DOF || idat==ELEMENT_RHSIDE_DELETE ||
         idat==ELEMENT_NONLOCAL || idat==ELEMENT_NONLOCAL_IPOINT ||
         idat==ELEMENT_NONLOCAL_WEIGHT || idat==NONLOCAL_ELEMENT_INFO ||
         !strncmp( nm, "element_intpnt_", 15 ) ||
         !strncmp( nm, "element_matrix_", 15 ) )
      db_delete_index( idat, element, VERSION_NORMAL );
  }
  // db_len finds a length by bisection on the LONG_MIN end mark: delete
  // first, so a shorter list does not leave the old nodes behind the mark
  db_delete_index( ELEMENT, element, VERSION_NORMAL );
  db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, PUT );
  if ( mnolnuknwn>0 ) {
    double *zeros = get_new_dbl( mnolnuknwn );
    array_set( zeros, 0., mnolnuknwn );
    db( ELEMENT_INTPNT_DOF, element, &zero, zeros, mnolnuknwn,
      VERSION_NORMAL, PUT );
    delete[] zeros;
  }
  db( ELEMENT_DOF_INITIALISED, element, &zero, ddum, one, VERSION_NORMAL,
    PUT );
}

void mesh_map( long int first_of_control )

{
  long int icontrol=0, ldum=0, swit_yes=-NO, element=0, max_element=0,
    max_node=0, max_used=0, length=0, map_length=0, nnol=0, lin_name=0,
    lin_nnol=0, inol=0, jnol=0, iuknwn=0, nchanged=0, quadratic=0,
    idum[1], el[1+MNOL], el_map[1+MNOL], lin[1+MNOL], corner[MNOL],
    lin_nodes[MNOL], *used=NULL;
  double ddum[1], w[8], coord[MDIM], start_refined[MDIM],
    map_node[2*MNOL*MDIM], *dof=NULL, *dof_start_refined=NULL,
    *corner_dof=NULL, *corner_dof_start_refined=NULL;

  if ( !first_of_control ) return;
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  if ( !db_active_index( CONTROL_MESH_MAP, icontrol, VERSION_NORMAL ) )
    return;
  db( CONTROL_MESH_MAP, icontrol, &swit_yes, ddum, ldum, VERSION_NORMAL,
    GET );
  if ( swit_yes!=-YES ) return;

  db_highest_index( ELEMENT, max_element, VERSION_NORMAL );
  db_highest_index( NODE, max_node, VERSION_NORMAL );
  if ( max_element<0 || max_node<0 ) return;

  // the direction follows the mesh: quadratic elements map to linear,
  // linear elements with an element_map record map back
  max_used = max_node;
  for ( element=0; element<=max_element; element++ ) {
    if ( !db_active_index( ELEMENT, element, VERSION_NORMAL ) ) continue;
    db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
    if ( mesh_map_linear( el[0], length-1, lin_name, lin_nnol, corner ) )
      quadratic = 1;
    if ( db_active_index( ELEMENT_MAP, element, VERSION_NORMAL ) ) {
      db( ELEMENT_MAP, element, el_map, ddum, map_length, VERSION_NORMAL,
        GET );
      for ( inol=1; inol<map_length; inol++ )
        if ( el_map[inol]>max_used ) max_used = el_map[inol];
    }
  }

  // the nodes the current elements use: they keep their dofs on a map back
  used = get_new_int( 1+max_used );
  array_set( used, 0, 1+max_used );
  for ( element=0; element<=max_element; element++ ) {
    if ( !db_active_index( ELEMENT, element, VERSION_NORMAL ) ) continue;
    db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
    for ( inol=1; inol<length; inol++ ) used[el[inol]] = 1;
  }

  if ( nuknwn>0 ) {
    dof = get_new_dbl( nuknwn );
    dof_start_refined = get_new_dbl( nuknwn );
    corner_dof = get_new_dbl( 8*nuknwn );
    corner_dof_start_refined = get_new_dbl( 8*nuknwn );
  }
  for ( element=0; element<=max_element; element++ ) {
    if ( !db_active_index( ELEMENT, element, VERSION_NORMAL ) ) continue;
    db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
    nnol = length - 1;
    if ( quadratic ) {
      if ( !mesh_map_linear( el[0], nnol, lin_name, lin_nnol, corner ) )
        continue;
      // keep the quadratic element, and the coordinates of its nodes: the
      // mid nodes are deleted with the mesh change (nod_nod deletes the
      // nodes no element uses) and created again on the map back
      db( ELEMENT_MAP, element, el, ddum, length, VERSION_NORMAL, PUT );
      for ( inol=0; inol<nnol; inol++ ) {
        db( NODE, el[1+inol], idum, &map_node[inol*ndim], ldum,
          VERSION_NORMAL, GET );
        db( NODE_START_REFINED, el[1+inol], idum,
          &map_node[(nnol+inol)*ndim], ldum, VERSION_NORMAL, GET );
      }
      ldum = 2*nnol*ndim;
      db( ELEMENT_MAP_NODE, element, idum, map_node, ldum, VERSION_NORMAL,
        PUT );
      lin[0] = lin_name;
      for ( inol=0; inol<lin_nnol; inol++ ) lin[1+inol] = el[1+corner[inol]];
      mesh_map_element_put( element, lin, 1+lin_nnol );
      nchanged++;
    }
    else {
      if ( !db_active_index( ELEMENT_MAP, element, VERSION_NORMAL ) ||
           !db_active_index( ELEMENT_MAP_NODE, element, VERSION_NORMAL ) )
        continue;
      db( ELEMENT_MAP, element, el_map, ddum, map_length, VERSION_NORMAL,
        GET );
      if ( !mesh_map_linear( el_map[0], map_length-1, lin_name, lin_nnol,
           corner ) ) continue;
      // only an element that is still the linear element of its map
      long int same = ( el[0]==lin_name && nnol==lin_nnol );
      for ( inol=0; same && inol<lin_nnol; inol++ )
        if ( el[1+inol]!=el_map[1+corner[inol]] ) same = 0;
      if ( !same ) continue;
      db( ELEMENT_MAP_NODE, element, idum, map_node, ldum, VERSION_NORMAL,
        GET );
      long int map_nnol = map_length - 1;
      for ( inol=0; inol<lin_nnol; inol++ ) {
        lin_nodes[inol] = el[1+inol];
        if ( nuknwn==0 ) continue;
        array_move( db_dbl( NODE_DOF, lin_nodes[inol], VERSION_NORMAL ),
          &corner_dof[inol*nuknwn], nuknwn );
        array_move( db_dbl( NODE_DOF_START_REFINED, lin_nodes[inol],
          VERSION_NORMAL ), &corner_dof_start_refined[inol*nuknwn], nuknwn );
      }
      for ( inol=0; inol<map_nnol; inol++ ) {
        long int inod = el_map[1+inol];
        if ( used[inod] ) continue;
        // the dofs of the linear solution at the mid node
        mesh_map_weights( el_map[0], inol, w );
        if ( nuknwn>0 ) {
          array_set( dof, 0., nuknwn );
          array_set( dof_start_refined, 0., nuknwn );
          for ( jnol=0; jnol<lin_nnol; jnol++ ) {
            for ( iuknwn=0; iuknwn<nuknwn; iuknwn++ ) {
              dof[iuknwn] += w[jnol]*corner_dof[jnol*nuknwn+iuknwn];
              dof_start_refined[iuknwn] +=
                w[jnol]*corner_dof_start_refined[jnol*nuknwn+iuknwn];
            }
          }
        }
        if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
          if ( nuknwn>0 ) db( NODE_DOF, inod, idum, dof, nuknwn,
            VERSION_NORMAL, PUT );
        }
        else {
          // the node again with the records a node keeps over a mesh
          // change (mesh_has_changed deletes the other node records);
          // create_node is not used in place: it copies the records of
          // the first corner through db_dbl pointers, which a PUT beyond
          // the highest index of a record reallocates
          array_move( &map_node[inol*ndim], coord, ndim );
          array_move( &map_node[(map_nnol+inol)*ndim], start_refined, ndim );
          db( NODE, inod, idum, coord, ndim, VERSION_NORMAL, PUT );
          db( NODE_START_REFINED, inod, idum, start_refined, ndim,
            VERSION_NORMAL, PUT );
          if ( nuknwn>0 ) {
            db( NODE_DOF, inod, idum, dof, nuknwn, VERSION_NORMAL, PUT );
            db( NODE_DOF_START_REFINED, inod, idum, dof_start_refined,
              nuknwn, VERSION_NORMAL, PUT );
          }
        }
        used[inod] = 1; // a shared mid node is set once
      }
      mesh_map_element_put( element, el_map, map_length );
      nchanged++;
    }
  }

  delete[] used;
  if ( dof!=NULL ) delete[] dof;
  if ( dof_start_refined!=NULL ) delete[] dof_start_refined;
  if ( corner_dof!=NULL ) delete[] corner_dof;
  if ( corner_dof_start_refined!=NULL ) delete[] corner_dof_start_refined;
  if ( nchanged>0 ) mesh_has_changed( VERSION_NORMAL );
}
