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

// ---------------------------------------------------------------------------
// Multi point constraints (MPC) - manual Professional 6.856-6.875.
//
// mpc_node_number / mpc_node_factor (6.874/6.875): the dof dof_0 of node
// node_0 (the slave) is constrained to a linear combination of dofs of the
// master nodes:
//
//     dof_0(node_0) = sum_i factor_i * dof_i(master_i)
//
// with the factors from mpc_node_factor (default 1). Semantics verified
// against the Professional binary (25-10-2023, .dbs): the slave dof is
// treated as a KNOWN (bounded) dof whose value is recomputed every
// equilibrium iteration from the current master values. The slave's
// equilibrium equation is dropped (the row is excluded from the global
// system) and NO force redistribution to the masters takes place; the
// master's equation is untouched. Boundary conditions must not be given
// on slave nodes (manual 6.875).
//
// mpc_linear_quadratic (6.873): when a quadratic element (quad9, hex27,
// ...) borders a linear element (quad4, hex8, ...) at a common interface,
// the extra nodes of the quadratic element (mid-edge, mid-face) dangle:
// they are not attached to the linear element, producing non-compatible
// solution fields. The option generates mpc_node_number/mpc_node_factor
// records that tie each dangling node to the linear element with the shape
// functions of the linear element evaluated at the node position. The
// generated records are identical to the ones the Professional writes to
// the .dbs (verified on mpc3/mpc4/mpc5: mid-edge slaves get 2 masters with
// factor 0.5, mid-face slaves 4 masters with factor 0.25).
//
// The generated records are re-created whenever the mesh changes
// (refinement, deletion, splitting, ...): the mesh fingerprint is stored
// in the internal record MPC_LINEAR_QUADRATIC_MESH_FINGERPRINT together
// with the index range of the generated records.
// ---------------------------------------------------------------------------

// tolerance on the isoparametric coordinates below which a dangling node
// is considered to be located in the linear element (same default as
// mpc_element_group_eps_iso, manual 6.865)
#define MPC_EPS_ISO 1.e-4

// weights below this value are not written to the generated records
// (the Professional lists only the non-zero masters)
#define MPC_WEIGHT_ZERO 1.e-6

static long int mpc_element_is_linear( long int name )
{
  return ( name==-BAR2 || name==-TRIA3 || name==-QUAD4 ||
           name==-TET4 || name==-HEX8 );
}

static long int mpc_element_is_quadratic( long int name )
{
  return ( name==-BAR3 || name==-TRIA6 || name==-QUAD9 ||
           name==-TET10 || name==-HEX27 );
}

// mesh fingerprint: a cheap hash over the element topology. Any mesh
// change (refine/delete/split) alters the element set and thus the hash.
static long int mpc_mesh_fingerprint( void )

{
  long int element=0, max_element=0, length=0, nnol=0, inol=0, inod=0,
    max_node=0, hash=0;
  double ddum[1];
  long int el[MNOL+1];

  db_max_index( NODE, max_node, VERSION_NORMAL, GET );
  db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
  hash = 100003 * ( 1 + max_node ) + 1009 * ( 1 + max_element );
  for ( element=0; element<=max_element; element++ ) {
    if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) ) {
      db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
      nnol = length - 1;
      hash ^= element * ( nnol + 7 );
      for ( inol=1; inol<=nnol; inol++ ) {
        inod = el[inol];
        hash ^= inod * ( inol + 3 );
      }
    }
  }
  return hash;

}

// mpc_linear_quadratic (manual Professional 6.873): generate the tie
// records for the dangling quadratic nodes. Called from mpc_node_apply()
// (every equilibrium iteration, but the generation runs only when the
// mesh fingerprint changed).
static void mpc_linear_quadratic_generate( void )

{
  long int lq=-NO, ldum=0, idum[1], length=0;
  double ddum[1];

  db( MPC_LINEAR_QUADRATIC, 0, &lq, ddum, ldum, VERSION_NORMAL,
    GET_IF_EXISTS );
  if ( lq != -YES ) return;

  long int fingerprint = mpc_mesh_fingerprint();
  long int stored[3] = { -1, -1, 0 };
  long int stored_len = 0;
  db( MPC_LINEAR_QUADRATIC_MESH_FINGERPRINT, 0, stored, ddum, stored_len,
    VERSION_NORMAL, GET_IF_EXISTS );
  if ( stored_len>=3 && stored[0]==fingerprint ) return;

  // delete the previously generated records
  if ( stored_len>=3 && stored[2]>0 ) {
    long int ig=0;
    for ( ig=0; ig<stored[2]; ig++ ) {
      db_delete_index( MPC_NODE_NUMBER, stored[1]+ig, VERSION_NORMAL );
      db_delete_index( MPC_NODE_FACTOR, stored[1]+ig, VERSION_NORMAL );
    }
  }

  long int max_node=0, max_element=0, element=0, ile=0, inol=0, jnol=0,
    inod=0, nnol=0, nnol_lin=0, name=0, name_lin=0, length_lin=0,
    nlinear=0, nmaster=0, i=0, ipuknwn=0, iuknwn=0, start_index=0,
    count=0, km=0, length_put=0;
  long int el[MNOL+1], el_lin[MNOL+1];
  double coord[MDIM], coords_el[MNOL*MDIM], weight[MNOL];

  db_max_index( NODE, max_node, VERSION_NORMAL, GET );
  db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
  if ( max_node<0 || max_element<0 ) {
    // nothing to tie; remember the fingerprint so the (empty) state is
    // not regenerated on every iteration
    stored[0] = fingerprint; stored[1] = -1; stored[2] = 0;
    length_put = 3;
    db( MPC_LINEAR_QUADRATIC_MESH_FINGERPRINT, 0, stored, ddum, length_put,
      VERSION_NORMAL, PUT );
    return;
  }

  // marker: nodes that belong to at least one LINEAR element
  long int *in_linear = get_new_int( 1 + max_node );
  array_set( in_linear, 0, 1 + max_node );
  // list of linear elements (master candidates)
  long int *linear_elems = get_new_int( 1 + max_element );

  for ( element=0; element<=max_element; element++ ) {
    if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) ) {
      db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
      name = el[0];
      if ( mpc_element_is_linear( name ) ) {
        linear_elems[nlinear++] = element;
        for ( inol=1; inol<length; inol++ )
          in_linear[el[inol]] = 1;
      }
    }
  }

  long int dof_label[MUKNWN], dof_principal[MUKNWN];
  db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( DOF_PRINCIPAL, 0, dof_principal, ddum, ldum, VERSION_NORMAL,
    GET_IF_EXISTS );

  // start index of the generated records: after the highest ACTIVE user
  // mpc_node_number record
  db_highest_index( MPC_NODE_NUMBER, start_index, VERSION_NORMAL );
  start_index++;

  long int vals[DATA_ITEM_SIZE];

  for ( element=0; element<=max_element; element++ ) {
    if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) ) {
      db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
      name = el[0];
      if ( !mpc_element_is_quadratic( name ) ) continue;
      nnol = length - 1;
      for ( inol=1; inol<=nnol; inol++ ) {
        inod = el[inol];
        // only dangling nodes: not attached to any linear element
        if ( in_linear[inod] ) continue;
        db( NODE_START_REFINED, inod, idum, coord, ldum, VERSION_NORMAL,
          GET );
        // find the linear element that contains the node
        for ( ile=0; ile<nlinear; ile++ ) {
          db( ELEMENT, linear_elems[ile], el_lin, ddum, length_lin,
            VERSION_NORMAL, GET );
          name_lin = el_lin[0];
          nnol_lin = length_lin - 1;
          for ( jnol=0; jnol<nnol_lin; jnol++ ) {
            db( NODE_START_REFINED, el_lin[1+jnol], idum,
              &coords_el[jnol*ndim], ldum, VERSION_NORMAL, GET );
          }
          if ( point_el( coord, coords_el, weight, name_lin, nnol_lin,
              MPC_EPS_ISO ) ) {
            // tie the node to the linear element: one mpc record per
            // principal dof, master weights = the linear shape functions
            // evaluated at the node position
            nmaster = 0;
            for ( jnol=0; jnol<nnol_lin; jnol++ ) {
              if ( scalar_dabs( weight[jnol] ) > MPC_WEIGHT_ZERO ) {
                el_lin[1+nmaster] = el_lin[1+jnol];
                weight[nmaster] = weight[jnol];
                nmaster++;
              }
            }
            if ( nmaster==0 ) continue;
            for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
              iuknwn = ipuknwn*nder;
              if ( dof_principal[iuknwn]>=0 ) {
                // record: [node_0 dof_0 node_1 dof_0 node_2 dof_0 ...]
                i = 0;
                vals[i++] = inod;
                vals[i++] = dof_label[iuknwn];
                for ( km=0; km<nmaster; km++ ) {
                  vals[i++] = el_lin[1+km];
                  vals[i++] = dof_label[iuknwn];
                }
                length_put = i;
                db( MPC_NODE_NUMBER, start_index+count, vals, ddum,
                  length_put, VERSION_NORMAL, PUT );
                length_put = nmaster;
                db( MPC_NODE_FACTOR, start_index+count, idum, weight,
                  length_put, VERSION_NORMAL, PUT );
                count++;
              }
            }
            break;
          }
        }
      }
    }
  }

  delete[] in_linear;
  delete[] linear_elems;

  stored[0] = fingerprint; stored[1] = start_index; stored[2] = count;
  length_put = 3;
  db( MPC_LINEAR_QUADRATIC_MESH_FINGERPRINT, 0, stored, ddum, length_put,
    VERSION_NORMAL, PUT );

}

// ---------------------------------------------------------------------------
// mpc_element_group (manual Professional 6.860): each node of an element of
// group element_group_0 that is also located inside an element of group
// element_group_1 is tied to that element with multi point constraints
// consistent with the shape functions at the node's isoparametric
// location. mpc_element_group_dof (6.864) selects the dofs that are set
// equal (default: all principal dofs); mpc_element_group_geometry (6.866)
// restricts the nodes of element_group_0 to those lying on a geometry.
// The generated mpc_node_number/mpc_node_factor records are the same ones
// the Professional writes to the .dbs, and are re-created when the mesh
// changes (fingerprint bookkeeping, mirroring mpc_linear_quadratic).
//
// Nodes that are already members of the node list of the containing
// element of element_group_1 (i.e. the two meshes are connected there by
// shared nodes) are NOT tied: the tie would reduce to the node itself and
// freeze its equation (mpc_element_group_always -no semantics, 6.861).
// ---------------------------------------------------------------------------

static void mpc_element_group_generate( void )

{
  long int ldum=0, idum[1], length=0, max_meg=0, ieg=0, i=0;
  double ddum[1];

  db_max_index( MPC_ELEMENT_GROUP, max_meg, VERSION_NORMAL, GET );
  if ( max_meg<0 ) return;

  long int fingerprint = mpc_mesh_fingerprint();
  long int stored[3] = { -1, -1, 0 };
  long int stored_len = 0;
  db( MPC_ELEMENT_GROUP_MESH_FINGERPRINT, 0, stored, ddum, stored_len,
    VERSION_NORMAL, GET_IF_EXISTS );
  if ( stored_len>=3 && stored[0]==fingerprint ) return;

  // delete the previously generated records
  if ( stored_len>=3 && stored[2]>0 ) {
    for ( i=0; i<stored[2]; i++ ) {
      db_delete_index( MPC_NODE_NUMBER, stored[1]+i, VERSION_NORMAL );
      db_delete_index( MPC_NODE_FACTOR, stored[1]+i, VERSION_NORMAL );
    }
  }

  long int dof_label[MUKNWN], dof_principal[MUKNWN];
  db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( DOF_PRINCIPAL, 0, dof_principal, ddum, ldum, VERSION_NORMAL,
    GET_IF_EXISTS );

  // start index of the generated records: after the highest ACTIVE
  // mpc_node_number record (user records and any mpc_linear_quadratic
  // generated records)
  long int start_index = 0;
  db_highest_index( MPC_NODE_NUMBER, start_index, VERSION_NORMAL );
  start_index++;
  long int count = 0;
  long int vals[DATA_ITEM_SIZE];

  long int max_element=0, max_node=0;
  db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
  db_max_index( NODE, max_node, VERSION_NORMAL, GET );

  for ( ieg=0; ieg<=max_meg; ieg++ ) {
    if ( !db_active_index( MPC_ELEMENT_GROUP, ieg, VERSION_NORMAL ) )
      continue;
    long int *eg_vals = db_int( MPC_ELEMENT_GROUP, ieg, VERSION_NORMAL );
    long int eg_length = db_len( MPC_ELEMENT_GROUP, ieg, VERSION_NORMAL );
    if ( eg_length<2 ) db_error( MPC_ELEMENT_GROUP, ieg );
    long int gr0 = eg_vals[0];

    // master element groups: eg_vals[1..] or -all (all other groups)
    long int master_groups[DATA_ITEM_SIZE];
    long int nmaster_groups = 0, all_groups = 0;
    for ( i=1; i<eg_length; i++ ) {
      if ( eg_vals[i]==-ALL ) all_groups = 1;
      else master_groups[nmaster_groups++] = eg_vals[i];
    }

    // selected dofs (negative dof labels); empty -> all principal dofs
    long int dof_sel[MUKNWN];
    long int ndof_sel = 0;
    if ( db_active_index( MPC_ELEMENT_GROUP_DOF, ieg, VERSION_NORMAL ) ) {
      long int *dof_vals = db_int( MPC_ELEMENT_GROUP_DOF, ieg,
        VERSION_NORMAL );
      long int dof_length = db_len( MPC_ELEMENT_GROUP_DOF, ieg,
        VERSION_NORMAL );
      for ( i=0; i<dof_length; i++ ) dof_sel[ndof_sel++] = dof_vals[i];
    }

    // restricting geometry pairs (geometry entity + index per pair)
    long int geom_vals[DATA_ITEM_SIZE];
    long int geom_length = 0;
    if ( db_active_index( MPC_ELEMENT_GROUP_GEOMETRY, ieg,
        VERSION_NORMAL ) ) {
      long int *gvals = db_int( MPC_ELEMENT_GROUP_GEOMETRY, ieg,
        VERSION_NORMAL );
      geom_length = db_len( MPC_ELEMENT_GROUP_GEOMETRY, ieg,
        VERSION_NORMAL );
      for ( i=0; i<geom_length; i++ ) geom_vals[i] = gvals[i];
    }

    // visited marker: each node is tied at most once per record
    long int *visited = NULL;
    if ( max_node>=0 ) {
      visited = get_new_int( 1+max_node );
      array_set( visited, 0, 1+max_node );
    }

    long int element=0, inod=0, inol=0, nnol=0, jelem=0, jnol=0,
      nnol_m=0, name_m=0, length_m=0, nmaster=0, k=0,
      km=0, ipuknwn=0, iuknwn=0, length_put=0, iuse=0, found=0,
      in_geom=0;
    long int el[MNOL+1], el_m[MNOL+1], el_grp=0;
    double coord[MDIM], coords_m[MNOL*MDIM], weight[MNOL], rdum=0.;
    long int geometry_entity[2];

    for ( element=0; element<=max_element; element++ ) {
      if ( !db_active_index( ELEMENT, element, VERSION_NORMAL ) ) continue;
      el_grp = 0;
      db( ELEMENT_GROUP, element, &el_grp, ddum, length, VERSION_NORMAL,
        GET_IF_EXISTS );
      if ( el_grp!=gr0 ) continue;
      db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
      nnol = length - 1;
      for ( inol=1; inol<=nnol; inol++ ) {
        inod = el[inol];
        if ( visited && visited[inod] ) continue;
        // geometry restriction
        if ( geom_length>0 ) {
          in_geom = 0;
          for ( k=0; k+1<geom_length; k+=2 ) {
            geometry_entity[0] = geom_vals[k];
            geometry_entity[1] = geom_vals[k+1];
            geometry( inod, ddum, geometry_entity, in_geom, rdum, ddum,
              rdum, ddum, NODE_START_REFINED, PROJECT_EXACT,
              VERSION_NORMAL );
            if ( in_geom ) break;
          }
          if ( !in_geom ) continue;
        }
        db( NODE_START_REFINED, inod, idum, coord, ldum, VERSION_NORMAL,
          GET );
        // find a master element containing the node (first match)
        found = 0;
        for ( jelem=0; jelem<=max_element && !found; jelem++ ) {
          if ( !db_active_index( ELEMENT, jelem, VERSION_NORMAL ) )
            continue;
          el_grp = 0;
          db( ELEMENT_GROUP, jelem, &el_grp, ddum, length, VERSION_NORMAL,
            GET_IF_EXISTS );
          iuse = 0;
          if ( all_groups ) { if ( el_grp!=gr0 ) iuse = 1; }
          else {
            for ( k=0; k<nmaster_groups; k++ )
              if ( master_groups[k]==el_grp ) iuse = 1;
          }
          if ( !iuse ) continue;
          db( ELEMENT, jelem, el_m, ddum, length_m, VERSION_NORMAL, GET );
          name_m = el_m[0]; nnol_m = length_m - 1;
          // skip master elements that already contain the node as a
          // member (the meshes are connected there; a tie would reduce
          // to the node itself)
          for ( jnol=1; jnol<=nnol_m; jnol++ )
            if ( el_m[jnol]==inod ) iuse = 0;
          if ( !iuse ) continue;
          for ( jnol=0; jnol<nnol_m; jnol++ ) {
            db( NODE_START_REFINED, el_m[1+jnol], idum,
              &coords_m[jnol*ndim], ldum, VERSION_NORMAL, GET );
          }
          if ( point_el( coord, coords_m, weight, name_m, nnol_m,
              MPC_EPS_ISO ) ) {
            // masters: the master element nodes with non-zero shape
            // function weight at the node position
            nmaster = 0;
            for ( jnol=0; jnol<nnol_m; jnol++ ) {
              if ( scalar_dabs( weight[jnol] ) > MPC_WEIGHT_ZERO ) {
                el_m[1+nmaster] = el_m[1+jnol];
                weight[nmaster] = weight[jnol];
                nmaster++;
              }
            }
            if ( nmaster==0 ) { iuse = 0; continue; }
            found = 1;
          }
          if ( !found ) continue;
          // generate one mpc record per selected dof
          for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
            iuknwn = ipuknwn*nder;
            if ( dof_principal[iuknwn]<0 ) continue;
            if ( ndof_sel>0 ) {
              long int iuse_dof = 0;
              for ( k=0; k<ndof_sel; k++ )
                if ( dof_sel[k]==dof_label[iuknwn] ) iuse_dof = 1;
              if ( !iuse_dof ) continue;
            }
            // record: [node_0 dof_0 node_1 dof_0 node_2 dof_0 ...]
            long int iv = 0;
            vals[iv++] = inod;
            vals[iv++] = dof_label[iuknwn];
            for ( km=0; km<nmaster; km++ ) {
              vals[iv++] = el_m[1+km];
              vals[iv++] = dof_label[iuknwn];
            }
            length_put = iv;
            db( MPC_NODE_NUMBER, start_index+count, vals, ddum,
              length_put, VERSION_NORMAL, PUT );
            length_put = nmaster;
            db( MPC_NODE_FACTOR, start_index+count, idum, weight,
              length_put, VERSION_NORMAL, PUT );
            count++;
          }
          if ( visited ) visited[inod] = 1;
        }
      }
    }
    delete[] visited;
  }

  stored[0] = fingerprint; stored[1] = start_index; stored[2] = count;
  length = 3;
  db( MPC_ELEMENT_GROUP_MESH_FINGERPRINT, 0, stored, ddum, length,
    VERSION_NORMAL, PUT );

}

// consume the mpc_node_number / mpc_node_factor records. Called from
// bounda() AFTER the bounda records, once per equilibrium iteration: the
// slave dofs are bounded and get the value sum(factor*master dof) from
// the master values of the current iteration (VERSION_NEW).
void mpc_node_apply( void )

{
  long int ldum=0, idum[1], max_mpc=0, impc=0;
  double ddum[1];

  // mpc_apply (manual 6.859, no index) / control_mpc_apply (6.255,
  // current control index): switch the whole mpc machinery off
  long int mpc_apply = -YES, icontrol = 0;
  db( MPC_APPLY, 0, &mpc_apply, ddum, ldum, VERSION_NORMAL,
    GET_IF_EXISTS );
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_MPC_APPLY, icontrol, &mpc_apply, ddum, ldum, VERSION_NORMAL,
    GET_IF_EXISTS );
  if ( mpc_apply==-NO ) return;

  // mpc_linear_quadratic and mpc_element_group generation (mesh-change
  // detection inside)
  mpc_linear_quadratic_generate();
  mpc_element_group_generate();

  db_max_index( MPC_NODE_NUMBER, max_mpc, VERSION_NORMAL, GET );
  if ( max_mpc<0 ) return;

  long int dof_label[MUKNWN];
  db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );

  for ( impc=0; impc<=max_mpc; impc++ ) {
    if ( db_active_index( MPC_NODE_NUMBER, impc, VERSION_NORMAL ) ) {
      long int length_mpc = db_len( MPC_NODE_NUMBER, impc, VERSION_NORMAL );
      long int *mpc = db_int( MPC_NODE_NUMBER, impc, VERSION_NORMAL );
      if ( length_mpc<2 ) db_error( MPC_NODE_NUMBER, impc );

      // slave: node_0 dof_0
      long int node_0 = mpc[0];
      long int iuknwn_0 = -1;
      array_member( dof_label, mpc[1], nuknwn, iuknwn_0 );
      if ( iuknwn_0<0 ) db_error( MPC_NODE_NUMBER, impc );
      long int ipuknwn_0 = iuknwn_0 / nder;

      // masters: the record layout is
      // [node_1 dof_1a dof_1b ... node_2 dof_2a ...] with the node
      // numbers positive and the dof keywords negative
      long int nterms = ( length_mpc - 2 ) / 2 + 1;
      long int *master_nodes = get_new_int( nterms );
      long int *master_dofs  = get_new_int( nterms );
      long int nmaster = 0, i = 2;
      while ( i<length_mpc ) {
        long int mnode = mpc[i++];
        if ( mnode<=0 ) db_error( MPC_NODE_NUMBER, impc );
        while ( i<length_mpc && mpc[i]<0 ) {
          long int iuknwn_m = -1;
          array_member( dof_label, mpc[i++], nuknwn, iuknwn_m );
          if ( iuknwn_m<0 ) db_error( MPC_NODE_NUMBER, impc );
          master_nodes[nmaster] = mnode;
          master_dofs[nmaster]  = iuknwn_m;
          nmaster++;
        }
      }
      if ( nmaster==0 ) db_error( MPC_NODE_NUMBER, impc );

      // factors (manual 6.874): same order as the master dofs; missing
      // factors default to 1
      double *master_factors = get_new_dbl( nmaster );
      for ( i=0; i<nmaster; i++ ) master_factors[i] = 1.;
      long int length_fac = 0;
      db( MPC_NODE_FACTOR, impc, idum, master_factors, length_fac,
        VERSION_NORMAL, GET_IF_EXISTS );
      if ( length_fac>nmaster ) length_fac = nmaster;

      // apply: bound the slave dof and set its value from the current
      // master values (VERSION_NEW). Nodes that no longer exist (mesh
      // deletion) are skipped.
      if ( db_active_index( NODE, node_0, VERSION_NORMAL ) ) {
        long int *node_bounded = db_int( NODE_BOUNDED, node_0,
          VERSION_NORMAL );
        double *new_node_dof = db_dbl( NODE_DOF, node_0, VERSION_NEW );
        node_bounded[ipuknwn_0] = 1;
        double value = 0.;
        for ( i=0; i<nmaster; i++ ) {
          if ( db_active_index( NODE, master_nodes[i], VERSION_NORMAL ) ) {
            double *mnd = db_dbl( NODE_DOF, master_nodes[i], VERSION_NEW );
            double factor = 1.;
            if ( i<length_fac ) factor = master_factors[i];
            value += factor * mnd[master_dofs[i]];
          }
        }
        new_node_dof[iuknwn_0] = value;
      }

      delete[] master_nodes;
      delete[] master_dofs;
      delete[] master_factors;
    }
  }

}
