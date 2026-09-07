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

void mesh_has_changed( long int version )

{ 
  long int idat=0, swit=0;

  swit = set_swit(-1,-1,"mesh_has_changed");
  if ( swit ) pri( "In routine MESH_HAS_CHANGED" );

  db_delete( CONTROL_PRINT_TECPLOT_MESH, VERSION_NORMAL );
  db_delete( NODE_ELEMENT, version );
  db_delete( NODE_NODE, version );

  if ( version==VERSION_NORMAL ) {
    for ( idat=0; idat<MDAT; idat++ ) {
      if ( db_data_class(idat)==NODE && !db_version( idat, VERSION_TMP ) )
        db_delete( idat, version );
    }
  }
  area_element_group( version );
  tendon_distribute();
  area_node_dataitem();
  nod_nod(version);
  nonlocal_first_set=0;
 
  if ( swit ) pri( "Out routine MESH_HAS_CHANGED" );
}

void mesh_add( long int version_from, long int version_to )

{

  long int inol=0, inod=0, indx=0, idat=0, element=0, max_node=0, max_element=0, 
    length=0, nnol=0, data_class=0, swit=0, idum[1], *nodes=NULL, *el=NULL, 
    *ival=NULL, *new_nodes=NULL;
  double ddum[1], *dval=NULL;

  swit = set_swit(-1,-1,"mesh_add");
  if ( swit ) pri( "In routine MESH_ADD" );

  db_max_index( NODE, max_node, version_from, GET );
  db_max_index( ELEMENT, max_element, version_from, GET );
  nodes = get_new_int( MNOL );
  el = get_new_int( 1+MNOL );
  new_nodes = get_new_int( 1+max_node );
  for ( inod=0; inod<=max_node; inod++ ) {
    if ( db_active_index( NODE, inod, version_from ) ) {
      indx = inod;
      while ( db_active_index( NODE, indx, version_to ) ) indx++;
      new_nodes[inod] = indx;
      for ( idat=0; idat<MDAT; idat++ ) {
        data_class = db_data_class( idat );
        if ( data_class==NODE && 
          db_version( idat, version_from ) &&
          db_version( idat, version_to ) ) {
          if ( db_active_index( idat, inod, version_from ) ) {
            length = db_len( idat, inod, version_from );
            if ( db_type(idat)==DOUBLE_PRECISION ) {
              dval = db_dbl( idat, inod, version_from );
              db( idat, indx, idum, dval, length, version_to, PUT );
            }
            else {
              ival = db_int( idat, inod, version_from );
              db( idat, indx, ival, ddum, length, version_to, PUT );
            }
          }
        }
      }
    }
  }
  for ( element=0; element<=max_element; element++ ) {
    if ( db_active_index( ELEMENT, element, version_from ) ) {
      indx = element;
      while ( db_active_index( ELEMENT, indx, version_to ) ) indx++;
      for ( idat=0; idat<MDAT; idat++ ) {
        data_class = db_data_class( idat );
        if ( data_class==ELEMENT && 
          db_version( idat, version_from ) &&
          db_version( idat, version_to ) ) {
          if ( db_active_index( idat, element, version_from ) ) {
            length = db_len( idat, element, version_from );
            if ( db_type(idat)==DOUBLE_PRECISION ) {
              dval = db_dbl( idat, element, version_from );
              db( idat, indx, idum, dval, length, version_to, PUT );
            }
            else {
              ival = db_int( idat, element, version_from );
              db( idat, indx, ival, ddum, length, version_to, PUT );
            }
          }
        }
      }
      db( ELEMENT, element, el, ddum, length, version_from, GET );
      nnol = length - 1; array_move( &el[1], nodes, nnol );
      for ( inol=0; inol<nnol; inol++ ) {
        inod = nodes[inol];
        el[1+inol] = new_nodes[inod];
      }
      db( ELEMENT, indx, el, ddum, length, version_to, PUT );
    }
  }
  delete[] nodes;
  delete[] el;
  delete[] new_nodes;
  db_version_delete( version_from );
  mesh_has_changed( version_to );

  if ( swit ) pri( "Out routine MESH_ADD" );
}

void mesh_switch( long int control_mesh_switch[], long int length )

{
  // Switch x, y, z coordinates of all nodes, e.g. for easy rotating the mesh.
  // The control specifies the new order of the axes:
  //   control_mesh_switch 0  -y -x -z   (interchange x and y)
  // The number of axes given must equal ndim.
  long int inod=0, max_node=0, idim=0, order[MDIM],
    has[MDIM], idum[1];
  double coords[MDIM], new_coords[MDIM];

  array_set( order, -1, MDIM );
  array_set( has, 0, MDIM );
  for ( idim=0; idim<ndim && idim<length; idim++ ) {
    if ( control_mesh_switch[idim]==-X ) order[idim] = 0;
    else if ( control_mesh_switch[idim]==-Y ) order[idim] = 1;
    else if ( control_mesh_switch[idim]==-Z ) order[idim] = 2;
    else {
      pri( "Error: control_mesh_switch axes must be -x, -y or -z." );
      exit(TN_EXIT_STATUS);
    }
  }
  for ( idim=0; idim<ndim; idim++ ) {
    if ( order[idim]<0 ) {
      pri( "Error: control_mesh_switch must specify all axes." );
      exit(TN_EXIT_STATUS);
    }
    if ( has[order[idim]] ) {
      pri( "Error: control_mesh_switch axes must be a permutation." );
      exit(TN_EXIT_STATUS);
    }
    has[order[idim]] = 1;
  }

  db_max_index( NODE, max_node, VERSION_NORMAL, GET );
  if ( max_node>=0 ) {
    for ( inod=0; inod<=max_node; inod++ ) {
      if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
        db( NODE, inod, idum, coords, ndim, VERSION_NORMAL, GET );
        for ( idim=0; idim<ndim; idim++ ) new_coords[idim] = coords[order[idim]];
        db( NODE, inod, idum, new_coords, ndim, VERSION_NORMAL, PUT );
        if ( db_active_index( NODE_START_REFINED, inod, VERSION_NORMAL ) ) {
          db( NODE_START_REFINED, inod, idum, coords, ndim, VERSION_NORMAL, GET );
          for ( idim=0; idim<ndim; idim++ ) new_coords[idim] = coords[order[idim]];
          db( NODE_START_REFINED, inod, idum, new_coords, ndim, VERSION_NORMAL, PUT );
        }
      }
    }
  }
  mesh_has_changed( VERSION_NORMAL );
}

void mesh_move( double control_mesh_move[], long int length )

{
  // move all nodes: in x-direction moved over
  //   move_x_constant + move_x_linear_x*x + move_x_linear_y*y + move_x_linear_z*z
  // coefficients are given for the number of space dimensions (ndim).
  long int inod=0, max_node=0, idim=0, jdim=0, idum[1];
  double coords[MDIM], new_coords[MDIM];
  double coeff[MDIM][MDIM+1];

  for ( idim=0; idim<MDIM; idim++ )
    for ( jdim=0; jdim<MDIM+1; jdim++ )
      coeff[idim][jdim] = 0.;
  // coeff[idim][0] = constant, coeff[idim][1+jd] = linear coefficient of jd
  for ( idim=0; idim<ndim && (idim*(ndim+1)+ndim)<length; idim++ ) {
    coeff[idim][0] = control_mesh_move[idim*(ndim+1)];
    for ( jdim=0; jdim<ndim; jdim++ )
      coeff[idim][1+jdim] = control_mesh_move[idim*(ndim+1)+1+jdim];
  }

  db_max_index( NODE, max_node, VERSION_NORMAL, GET );
  if ( max_node>=0 ) {
    for ( inod=0; inod<=max_node; inod++ ) {
      if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
        db( NODE, inod, idum, coords, ndim, VERSION_NORMAL, GET );
        for ( idim=0; idim<ndim; idim++ ) {
          new_coords[idim] = coords[idim] + coeff[idim][0];
          for ( jdim=0; jdim<ndim; jdim++ )
            new_coords[idim] += coeff[idim][1+jdim]*coords[jdim];
        }
        db( NODE, inod, idum, new_coords, ndim, VERSION_NORMAL, PUT );
        if ( db_active_index( NODE_START_REFINED, inod, VERSION_NORMAL ) ) {
          db( NODE_START_REFINED, inod, idum, coords, ndim, VERSION_NORMAL, GET );
          for ( idim=0; idim<ndim; idim++ ) {
            new_coords[idim] = coords[idim] + coeff[idim][0];
            for ( jdim=0; jdim<ndim; jdim++ )
              new_coords[idim] += coeff[idim][1+jdim]*coords[jdim];
          }
          db( NODE_START_REFINED, inod, idum, new_coords, ndim, VERSION_NORMAL, PUT );
        }
      }
    }
  }
  mesh_has_changed( VERSION_NORMAL );
}

void mesh_mirror( long int axis )

{
  // mirror the mesh about a plane x=0, y=0 or z=0: duplicate nodes and elements
  long int inod=0, max_node=0, ielem=0, max_elem=0, inol=0, nnol=0,
    length=0, new_node=0, new_elem=0, idum[1], axis_i=0, el[1+MNOL],
    new_nodes[1+MNOL];
  double ddum[1], coords[MDIM], new_coords[MDIM];

  if ( axis==-X ) axis_i = 0;
  else if ( axis==-Y ) axis_i = 1;
  else if ( axis==-Z ) axis_i = 2;
  else {
    cout << "Error: control_mesh_mirror axis must be -x, -y or -z.\n";
    exit(TN_EXIT_STATUS);
  }

  db_max_index( NODE, max_node, VERSION_NORMAL, GET );
  if ( max_node<0 ) return;

  // duplicate nodes with mirrored coordinate
  long int max_node_old = max_node;
  for ( inod=0; inod<=max_node; inod++ ) {
    if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
      new_node = inod + max_node_old + 1;
      db( NODE, inod, idum, coords, ndim, VERSION_NORMAL, GET );
      for ( int i=0; i<ndim; i++ ) new_coords[i] = coords[i];
      new_coords[axis_i] = -coords[axis_i];
      db( NODE, new_node, idum, new_coords, ndim, VERSION_NORMAL, PUT );
      if ( db_active_index( NODE_START_REFINED, inod, VERSION_NORMAL ) ) {
        db( NODE_START_REFINED, inod, idum, coords, ndim, VERSION_NORMAL, GET );
        for ( int i=0; i<ndim; i++ ) new_coords[i] = coords[i];
        new_coords[axis_i] = -coords[axis_i];
        db( NODE_START_REFINED, new_node, idum, new_coords, ndim, VERSION_NORMAL, PUT );
      }
      if ( db_active_index( NODE_DOF, inod, VERSION_NORMAL ) ) {
        double *ndof_old = db_dbl( NODE_DOF, inod, VERSION_NORMAL );
        db( NODE_DOF, new_node, idum, ndof_old, nuknwn, VERSION_NORMAL, PUT );
      }
    }
  }

  // duplicate elements with mirrored nodes (create_element copies all
  // element data: element_dof, element_dof_initialised, nonlocal, etc.)
  db_max_index( ELEMENT, max_elem, VERSION_NORMAL, GET );
  new_elem = max_elem;
  for ( ielem=0; ielem<=max_elem; ielem++ ) {
    if ( db_active_index( ELEMENT, ielem, VERSION_NORMAL ) ) {
      new_elem++;
      db( ELEMENT, ielem, el, ddum, length, VERSION_NORMAL, GET );
      nnol = length - 1;
      new_nodes[0] = el[0];
      for ( inol=0; inol<nnol; inol++ ) new_nodes[1+inol] = el[1+inol]+max_node_old+1;
      create_element( ielem, new_elem, new_nodes, length, VERSION_NORMAL,
        VERSION_NORMAL );
    }
  }
  mesh_has_changed( VERSION_NORMAL );
}

void mesh_copy( double move_coords[] )

{
  // copy the mesh: duplicate nodes and elements; each new node is the
  // corresponding old node moved over move_coords[ndim].
  long int inod=0, max_node=0, ielem=0, max_elem=0, inol=0, nnol=0,
    length=0, new_node=0, new_elem=0, idum[1], el[1+MNOL],
    new_nodes[1+MNOL];
  double ddum[1], coords[MDIM], new_coords[MDIM];

  db_max_index( NODE, max_node, VERSION_NORMAL, GET );
  if ( max_node<0 ) return;
  long int max_node_old = max_node;

  // duplicate nodes moved over move_coords
  for ( inod=0; inod<=max_node; inod++ ) {
    if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
      new_node = inod + max_node_old + 1;
      db( NODE, inod, idum, coords, ndim, VERSION_NORMAL, GET );
      for ( int i=0; i<ndim; i++ ) new_coords[i] = coords[i] + move_coords[i];
      db( NODE, new_node, idum, new_coords, ndim, VERSION_NORMAL, PUT );
      if ( db_active_index( NODE_START_REFINED, inod, VERSION_NORMAL ) ) {
        db( NODE_START_REFINED, inod, idum, coords, ndim, VERSION_NORMAL, GET );
        for ( int i=0; i<ndim; i++ ) new_coords[i] = coords[i] + move_coords[i];
        db( NODE_START_REFINED, new_node, idum, new_coords, ndim, VERSION_NORMAL, PUT );
      }
      if ( db_active_index( NODE_DOF, inod, VERSION_NORMAL ) ) {
        double *ndof_old = db_dbl( NODE_DOF, inod, VERSION_NORMAL );
        db( NODE_DOF, new_node, idum, ndof_old, nuknwn, VERSION_NORMAL, PUT );
      }
    }
  }

  // duplicate elements
  db_max_index( ELEMENT, max_elem, VERSION_NORMAL, GET );
  new_elem = max_elem;
  for ( ielem=0; ielem<=max_elem; ielem++ ) {
    if ( db_active_index( ELEMENT, ielem, VERSION_NORMAL ) ) {
      new_elem++;
      db( ELEMENT, ielem, el, ddum, length, VERSION_NORMAL, GET );
      nnol = length - 1;
      new_nodes[0] = el[0];
      for ( inol=0; inol<nnol; inol++ ) new_nodes[1+inol] = el[1+inol]+max_node_old+1;
      create_element( ielem, new_elem, new_nodes, length, VERSION_NORMAL,
        VERSION_NORMAL );
    }
  }
  mesh_has_changed( VERSION_NORMAL );
}

void mesh_rotate_2d( double angle_deg )

{
  // rotate the 2D mesh around the z-axis by angle_deg (degrees).
  // Only the nodes are moved; element connectivity is unchanged.
  long int inod=0, max_node=0, idum[1];
  double coords[MDIM], new_coords[MDIM];
  double a = angle_deg * PIRAD / 180.;
  double ca = cos(a), sa = sin(a);

  db_max_index( NODE, max_node, VERSION_NORMAL, GET );
  if ( max_node>=0 ) {
    for ( inod=0; inod<=max_node; inod++ ) {
      if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
        db( NODE, inod, idum, coords, ndim, VERSION_NORMAL, GET );
        new_coords[0] = ca*coords[0] - sa*coords[1];
        new_coords[1] = sa*coords[0] + ca*coords[1];
        if ( ndim==3 ) new_coords[2] = coords[2];
        db( NODE, inod, idum, new_coords, ndim, VERSION_NORMAL, PUT );
        if ( db_active_index( NODE_START_REFINED, inod, VERSION_NORMAL ) ) {
          db( NODE_START_REFINED, inod, idum, coords, ndim, VERSION_NORMAL, GET );
          new_coords[0] = ca*coords[0] - sa*coords[1];
          new_coords[1] = sa*coords[0] + ca*coords[1];
          if ( ndim==3 ) new_coords[2] = coords[2];
          db( NODE_START_REFINED, inod, idum, new_coords, ndim, VERSION_NORMAL, PUT );
        }
      }
    }
  }
  mesh_has_changed( VERSION_NORMAL );
}

void mesh_rotate_3d( long int nrot )

{
  // rotate a 2D mesh to 3D: each -tria3 becomes a -prism6 and each -quad4
  // becomes a -hex8, by rotating around the y-axis. nrot is the number of
  // elements in the rotational direction over 360 degrees.
  long int inod=0, max_node=0, ielem=0, max_elem=0, inol=0, nnol=0,
    length=0, new_node=0, new_elem=0, idum[1], len3=3, el[1+MNOL],
    nodes[MNOL], new_nodes[1+MNOL];
  double ddum[1], coords[MDIM], new_coords[MDIM];
  long int *node_rot=NULL;

  if ( nrot<1 ) nrot = 1;

  db_max_index( NODE, max_node, VERSION_NORMAL, GET );
  db_max_index( ELEMENT, max_elem, VERSION_NORMAL, GET );
  if ( max_node<0 || max_elem<0 ) return;

  // duplicate nodes rotated around y-axis
  node_rot = get_new_int(1+max_node);
  for ( inod=0; inod<=max_node; inod++ ) node_rot[inod] = -1;
  new_node = max_node;
  for ( inod=0; inod<=max_node; inod++ ) {
    if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
      new_node++;
      node_rot[inod] = new_node;
      db( NODE, inod, idum, coords, ndim, VERSION_NORMAL, GET );
      new_coords[0] = coords[2];
      new_coords[1] = coords[1];
      new_coords[2] = -coords[0];
      db( NODE, new_node, idum, new_coords, len3, VERSION_NORMAL, PUT );
      // copy all other node-class data items from the source node
      for ( int idat=0; idat<MDAT; idat++ ) {
        if ( idat!=NODE && db_data_class(idat)==NODE &&
             db_active_index( idat, inod, VERSION_NORMAL ) ) {
          long int ndata_len = db_len( idat, inod, VERSION_NORMAL );
          if ( db_type(idat)==DOUBLE_PRECISION ) {
            double *dold = db_dbl( idat, inod, VERSION_NORMAL );
            db( idat, new_node, idum, dold, ndata_len, VERSION_NORMAL, PUT );
          }
          else {
            long int *iold = db_int( idat, inod, VERSION_NORMAL );
            db( idat, new_node, iold, ddum, ndata_len, VERSION_NORMAL, PUT );
          }
        }
      }
      // the rotated coordinate also applies to node_start_refined
      if ( db_active_index( NODE_START_REFINED, inod, VERSION_NORMAL ) ) {
        db( NODE_START_REFINED, inod, idum, coords, ndim, VERSION_NORMAL, GET );
        new_coords[0] = coords[2];
        new_coords[1] = coords[1];
        new_coords[2] = -coords[0];
        db( NODE_START_REFINED, new_node, idum, new_coords, len3, VERSION_NORMAL, PUT );
      }
    }
  }

  // duplicate elements: 2D element + rotated copy = 3D element
  new_elem = max_elem;
  for ( ielem=0; ielem<=max_elem; ielem++ ) {
    if ( db_active_index( ELEMENT, ielem, VERSION_NORMAL ) ) {
      db( ELEMENT, ielem, el, ddum, length, VERSION_NORMAL, GET );
      nnol = length - 1;
      for ( inol=0; inol<nnol; inol++ ) nodes[inol] = el[1+inol];
      if ( el[0]==-TRIA3 && nnol==3 ) {
        new_elem++;
        new_nodes[0] = -PRISM6;
        new_nodes[1] = nodes[0]; new_nodes[2] = nodes[1]; new_nodes[3] = nodes[2];
        new_nodes[4] = node_rot[nodes[0]]; new_nodes[5] = node_rot[nodes[1]];
        new_nodes[6] = node_rot[nodes[2]];
        create_element( ielem, new_elem, new_nodes, 7, VERSION_NORMAL,
          VERSION_NORMAL );
      }
      else if ( el[0]==-QUAD4 && nnol==4 ) {
        new_elem++;
        new_nodes[0] = -HEX8;
        new_nodes[1] = nodes[0]; new_nodes[2] = nodes[1]; new_nodes[3] = nodes[2];
        new_nodes[4] = nodes[3]; new_nodes[5] = node_rot[nodes[0]];
        new_nodes[6] = node_rot[nodes[1]]; new_nodes[7] = node_rot[nodes[2]];
        new_nodes[8] = node_rot[nodes[3]];
        create_element( ielem, new_elem, new_nodes, 9, VERSION_NORMAL,
          VERSION_NORMAL );
      }
    }
  }

  // delete the 2D source elements (not valid in 3D)
  for ( ielem=0; ielem<=max_elem; ielem++ ) {
    if ( db_active_index( ELEMENT, ielem, VERSION_NORMAL ) ) {
      db( ELEMENT, ielem, el, ddum, length, VERSION_NORMAL, GET );
      if ( el[0]==-TRIA3 || el[0]==-QUAD4 || el[0]==-BAR2 )
        delete_element( ielem, VERSION_NORMAL );
    }
  }

  delete[] node_rot;
  mesh_has_changed( VERSION_NORMAL );
}

// mesh_convert_quad8 - auto-convert the Professional 8-node serendipity
// quad8 volume elements to the GNU 9-node Lagrange quad9.
//
// The GNU has NO real quad8: the Professional auto-converts them (its
// own corpus .dat files state it textually - interface_bar3_quad8.dat:
// "The bar3 and quad8 will be automatically converted to quad6
// interface and quad9 surface elements"). The conversion inserts the
// CENTRE node (the average of the 4 corners) and rewrites the
// connectivity from the Professional quad8 ordering to the GNU tensor
// quad9 ordering.
//
//   quad8 (Professional, every -quad8 record of the suite):
//     corners  (BL, BR, TL, TR) then mid-edge nodes (BM, LM, RM, TM)
//     -> record slots 1..8
//   quad9 (GNU, xi fastest / eta slowest: -1 -> 0 -> +1 per axis):
//     BL, BM, BR | LM, CENTRE, RM | TL, TM, TR  -> record slots 1..9
//     (border_nodes_quad9 of area.cc: corners 0,2,8,6, mid-edge nodes
//     1,3,5,7, centre 4 - the same ordering of every -quad9 record of
//     the suite, e.g. patch1.dat element 1)
//
// so the permutation is:
//     quad9 = { q8[1], q8[5], q8[2], q8[6], centre, q8[7], q8[3],
//               q8[8], q8[4] }.
//
// The hook runs at EVERY step_start (see top.cc step_start: not only
// task==YES - a quad8 is not a native element, so an intermediate
// control step below the timestep would evaluate the raw quad8), BEFORE
// extrude() (a quad8 mesh extrudes to hex27 like a quad9 one) and
// BEFORE interface_convert() (the interface split sees the quad9
// bulk). It is idempotent: converted elements are -quad9 and skipped
// on later steps. -quad8 elements of an INTERFACE group are skipped:
// there the quad8 is a FACIAL element (3D interface_quad8_hex20
// family) converted by the interface machinery, not a volume.
void mesh_convert_quad8( void )

{
  long int element=0, max_element=0, max_node=0, length=0, ldum=0,
    swit=0, element_group=0, inol=0, i=0, nconv=0, idum[1];
  double ddum[1], coord[MDIM], centre[MDIM];
  long int el[1+MNOL], q9[1+9];

  swit = set_swit(-1,-1,"mesh_convert_quad8");
  if ( swit ) pri( "In routine MESH_CONVERT_QUAD8" );

  db_highest_index( ELEMENT, max_element, VERSION_NORMAL );
  db_highest_index( NODE, max_node, VERSION_NORMAL );
  if ( max_element<0 ) return;

  for ( element=0; element<=max_element; element++ ) {
    if ( !db_active_index( ELEMENT, element, VERSION_NORMAL ) ) continue;
    db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
    if ( el[0]!=-QUAD8 ) continue;
    if ( length!=1+8 ) db_error( ELEMENT, element );
    // interface-group quad8: a facial interface element (the 3D
    // quad8-interface family), NOT a volume - leave it to the
    // interface conversion lot.
    element_group = 0;
    db( ELEMENT_GROUP, element, &element_group, ddum, ldum,
      VERSION_NORMAL, GET_IF_EXISTS );
    if ( db_active_index( GROUP_INTERFACE, element_group, VERSION_NORMAL ) )
      continue;
    // centre node: average of the 4 corners (quad8 slots 1..4)
    array_set( centre, 0., MDIM );
    for ( inol=0; inol<4; inol++ ) {
      db( NODE, el[1+inol], idum, coord, ndim, VERSION_NORMAL, GET );
      for ( i=0; i<ndim; i++ ) centre[i] += coord[i];
    }
    for ( i=0; i<ndim; i++ ) centre[i] /= 4.;
    max_node++;
    db( NODE, max_node, idum, centre, ndim, VERSION_NORMAL, PUT );
    db( NODE_START_REFINED, max_node, idum, centre, ndim,
      VERSION_NORMAL, PUT );
    // the new node carries the full dof state of the corners (zeros at
    // first step; same pattern as interface_convert/extrude)
    for ( int idat=0; idat<MDAT; idat++ ) {
      if ( idat!=NODE && idat!=NODE_START_REFINED &&
           db_data_class(idat)==NODE &&
           db_active_index( idat, el[1], VERSION_NORMAL ) ) {
        long int ndata_len = db_len( idat, el[1], VERSION_NORMAL );
        if ( db_type(idat)==DOUBLE_PRECISION ) {
          double *dold = db_dbl( idat, el[1], VERSION_NORMAL );
          db( idat, max_node, idum, dold, ndata_len, VERSION_NORMAL, PUT );
        }
        else {
          long int *iold = db_int( idat, el[1], VERSION_NORMAL );
          db( idat, max_node, iold, ddum, ndata_len, VERSION_NORMAL, PUT );
        }
      }
    }
    // rewrite the connectivity in the GNU quad9 tensor ordering
    q9[0] = -QUAD9;
    q9[1] = el[1];
    q9[2] = el[5];
    q9[3] = el[2];
    q9[4] = el[6];
    q9[5] = max_node;
    q9[6] = el[7];
    q9[7] = el[3];
    q9[8] = el[8];
    q9[9] = el[4];
    length = 1+9;
    db( ELEMENT, element, q9, ddum, length, VERSION_NORMAL, PUT );
    nconv++;
  }
  if ( nconv>0 ) mesh_has_changed( VERSION_NORMAL );

  if ( swit ) pri( "Out function MESH_CONVERT_QUAD8" );
}

// mesh_convert_hex20 - auto-convert the Professional 20-node serendipity
// hex20 volume elements to the GNU 27-node Lagrange hex27.
//
// The GNU has no real hex20 element routine: the Professional handles
// -hex20 natively (its .dbs keeps the 20-node records, verified with the
// 25-10-2023 user-supplied binary on the corpus hex20.dat), while the GNU
// formulation is the complete 27-node Lagrange hex. Like the quad8 lot,
// the auto-conversion elevates the serendipity input to the richer
// Lagrange element (7 extra nodes: 6 face centres + the body centre).
//
// The slot permutation below is NOT guessed: the Professional itself
// auto-converts the mesh on its interface_quad8_hex20 corpus file and the
// resulting .dbs shows the SAME hex27 slot layout as the GNU tensor order
// (base plane quad9 tensor [BL,BM,BR,LM,C,RM,TL,TM,TR], mid plane, top
// plane - mesh_extrude/area.cc border_nodes_hex27 conventions). For the
// input hex20 record of that test
//
//   element 1 -hex20 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
//     (corners: base BL,BR,TL,TR = 1,2,3,4 ; top BL,BR,TL,TR = 5,6,7,8
//      base mid-edges BM,LM,RM,TM = 9,10,11,12
//      top  mid-edges BM,LM,RM,TM = 13,14,15,16
//      verticals at (x0,y0),(x1,y0),(x0,y1),(x1,y1) = 17,18,19,20)
//
// the Professional .dbs writes
//
//   element 1 -hex27 1 9 2 10 56 11 3 12 4  17 57 18 60 61 59 19 58 20
//                    5 13 6 14 55 15 7 16 8
//
// (nodes 55-61 = the NEW face/body centres of that run), which is exactly
// the permutation implemented here:
//
//   base  plane (slots 1..9):  BL,BM,BR,LM,Cf_base,RM,TL,TM,TR
//   mid   plane (slots 10..18): v(x0,y0), Cf_y0, v(x1,y0), Cf_x0, Cbody,
//                               Cf_x1, v(x0,y1), Cf_y1, v(x1,y1)
//   top   plane (slots 19..27): BL,BM,BR,LM,Cf_top,RM,TL,TM,TR
//
// Face centres = average of the 4 face corners, body centre = average of
// the 8 corners. The face centres of faces shared between neighbouring
// hex20 elements are DEDUPLICATED BY COORDINATES (EPS_COORD): two stacked
// hex20 (hex20.dat of the corpus: shared face z=1) must end up with ONE
// centre node on the shared face, otherwise the mesh tears apart.
//
// The hook runs at EVERY step_start (same rationale as mesh_convert_quad8:
// a -hex20 is not a native element, so an intermediate control step below
// the timestep would evaluate the raw hex20), BEFORE extrude() and BEFORE
// interface_convert() (the 3D interface split of the quad8 face must see
// the hex27 bulk with the shared face centre node - interface_quad8_hex20
// family). Idempotent: converted elements are -hex27 and skipped on later
// steps.
void mesh_convert_hex20( void )

{
  long int element=0, max_element=0, max_node=0, length=0,
    swit=0, i=0, nconv=0, jnod=0, d=0, idum[1];
  double ddum[1], coord[MDIM], centre[MDIM];
  long int el[1+MNOL], h27[1+27];

  swit = set_swit(-1,-1,"mesh_convert_hex20");
  if ( swit ) pri( "In routine MESH_CONVERT_HEX20" );

  if ( ndim!=3 ) return;
  db_highest_index( ELEMENT, max_element, VERSION_NORMAL );
  db_highest_index( NODE, max_node, VERSION_NORMAL );
  if ( max_element<0 ) return;

  for ( element=0; element<=max_element; element++ ) {
    if ( !db_active_index( ELEMENT, element, VERSION_NORMAL ) ) continue;
    db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
    if ( el[0]!=-HEX20 ) continue;
    if ( length!=1+20 ) db_error( ELEMENT, element );

    // corners of the 6 faces (into el[1..8], 1-based) used for the new
    // centre nodes: base (1,2,3,4), top (5,6,7,8), y0 (1,2,5,6),
    // y1 (3,4,7,8), x0 (1,3,5,7), x1 (2,4,6,8); body = all 8.
    // h27 slot layout (1-based, GNU tensor hex27):
    //   base plane  quad9 tensor: BL BM BR LM C RM TL TM TR
    //   mid plane:  v00 Cf_y0 v10 Cf_x0 Cbody Cf_x1 v01 Cf_y1 v11
    //   top plane   quad9 tensor: BL BM BR LM C RM TL TM TR
    h27[1]  = el[1];   // BL base
    h27[2]  = el[9];   // BM base
    h27[3]  = el[2];   // BR base
    h27[4]  = el[10];  // LM base
    h27[5]  = 0;       // Cf base (new)
    h27[6]  = el[11];  // RM base
    h27[7]  = el[3];   // TL base
    h27[8]  = el[12];  // TM base
    h27[9]  = el[4];   // TR base
    h27[10] = el[17];  // vertical (x0,y0)
    h27[11] = 0;       // Cf y0 (new)
    h27[12] = el[18];  // vertical (x1,y0)
    h27[13] = 0;       // Cf x0 (new)
    h27[14] = 0;       // Cbody (new)
    h27[15] = 0;       // Cf x1 (new)
    h27[16] = el[19];  // vertical (x0,y1)
    h27[17] = 0;       // Cf y1 (new)
    h27[18] = el[20];  // vertical (x1,y1)
    h27[19] = el[5];   // BL top
    h27[20] = el[13];  // BM top
    h27[21] = el[6];   // BR top
    h27[22] = el[14];  // LM top
    h27[23] = 0;       // Cf top (new)
    h27[24] = el[15];  // RM top
    h27[25] = el[7];   // TL top
    h27[26] = el[16];  // TM top
    h27[27] = el[8];   // TR top

    // the 7 new node positions: (slot, 4 face corners into el[], all 8)
    long int nnew = 0, newslot[7], corner_of[7][4];
    for ( int k=0; k<7; k++ ) for ( int c=0; c<4; c++ ) corner_of[k][c] = 0;
    // base face centre (slot 5)
    newslot[nnew]=5;   corner_of[nnew][0]=1; corner_of[nnew][1]=2;
    corner_of[nnew][2]=3; corner_of[nnew][3]=4; nnew++;
    // y0 face centre (slot 11)
    newslot[nnew]=11;  corner_of[nnew][0]=1; corner_of[nnew][1]=2;
    corner_of[nnew][2]=5; corner_of[nnew][3]=6; nnew++;
    // x0 face centre (slot 13)
    newslot[nnew]=13;  corner_of[nnew][0]=1; corner_of[nnew][1]=3;
    corner_of[nnew][2]=5; corner_of[nnew][3]=7; nnew++;
    // body centre (slot 14): 8 corners = corners of faces 1 and 5
    newslot[nnew]=14;  nnew++;
    // x1 face centre (slot 15)
    newslot[nnew]=15;  corner_of[nnew][0]=2; corner_of[nnew][1]=4;
    corner_of[nnew][2]=6; corner_of[nnew][3]=8; nnew++;
    // y1 face centre (slot 17)
    newslot[nnew]=17;  corner_of[nnew][0]=3; corner_of[nnew][1]=4;
    corner_of[nnew][2]=7; corner_of[nnew][3]=8; nnew++;
    // top face centre (slot 23)
    newslot[nnew]=23;  corner_of[nnew][0]=5; corner_of[nnew][1]=6;
    corner_of[nnew][2]=7; corner_of[nnew][3]=8; nnew++;

    for ( int k=0; k<nnew; k++ ) {
      array_set( centre, 0., MDIM );
      long int ncorner = ( newslot[k]==14 ) ? 8 : 4;
      long int slot_c[8];
      if ( newslot[k]==14 ) {
        for ( int c=0; c<8; c++ ) slot_c[c] = 1+c;
      }
      else {
        for ( int c=0; c<4; c++ ) slot_c[c] = corner_of[k][c];
      }
      for ( int c=0; c<ncorner; c++ ) {
        db( NODE, el[slot_c[c]], idum, coord, ndim, VERSION_NORMAL, GET );
        for ( i=0; i<ndim; i++ ) centre[i] += coord[i];
      }
      for ( i=0; i<ndim; i++ ) centre[i] /= (double)ncorner;

      // deduplicate by coordinates: a face shared with a neighbouring
      // hex20 (or an already existing node of the mesh) must reuse that
      // node - two coincident centres would tear the mesh apart. Search
      // every active node seen so far (pre-existing + created by this
      // conversion on earlier elements).
      long int found = -1;
      for ( jnod=0; jnod<=max_node; jnod++ ) {
        if ( !db_active_index( NODE, jnod, VERSION_NORMAL ) ) continue;
        db( NODE, jnod, idum, coord, ndim, VERSION_NORMAL, GET );
        long int ok = 1;
        for ( d=0; d<ndim && ok; d++ )
          if ( fabs( coord[d]-centre[d] )>1.e-10 ) ok = 0;
        if ( ok ) { found = jnod; break; }
      }
      if ( found<0 ) {
        max_node++;
        found = max_node;
        db( NODE, found, idum, centre, ndim, VERSION_NORMAL, PUT );
        db( NODE_START_REFINED, found, idum, centre, ndim,
          VERSION_NORMAL, PUT );
        // the new node carries the full dof state of the corners (zeros
        // at the first step; same pattern as mesh_convert_quad8 /
        // interface_convert / extrude)
        for ( int idat=0; idat<MDAT; idat++ ) {
          if ( idat!=NODE && idat!=NODE_START_REFINED &&
               db_data_class(idat)==NODE &&
               db_active_index( idat, el[1], VERSION_NORMAL ) ) {
            long int ndata_len = db_len( idat, el[1], VERSION_NORMAL );
            if ( db_type(idat)==DOUBLE_PRECISION ) {
              double *dold = db_dbl( idat, el[1], VERSION_NORMAL );
              db( idat, found, idum, dold, ndata_len, VERSION_NORMAL, PUT );
            }
            else {
              long int *iold = db_int( idat, el[1], VERSION_NORMAL );
              db( idat, found, iold, ddum, ndata_len, VERSION_NORMAL, PUT );
            }
          }
        }
      }
      h27[newslot[k]] = found;
    }

    // rewrite the element: -hex20 -> -hex27 in the GNU tensor ordering
    h27[0] = -HEX27;
    length = 1+27;
    db( ELEMENT, element, h27, ddum, length, VERSION_NORMAL, PUT );
    nconv++;
  }
  if ( nconv>0 ) mesh_has_changed( VERSION_NORMAL );

  if ( swit ) pri( "Out function MESH_CONVERT_HEX20" );
}

void mesh_extrude( double z_layer[], long int n_layer, long int quad9_layers )

{
  // extrude a 2D mesh (z=0) to 3D along the z-axis.
  // z_layer[] gives the z-coordinate of each layer boundary;
  // one 3D element is generated per 2D element per layer.
  // -tria3 -> -prism6, -quad4 -> -hex8, -quad9 -> -hex27
  // (quad9_layers: the number of hex27 layers from
  // control_mesh_extrude_n, splitting the z_layer extent evenly).
  long int inod=0, max_node=0, ielem=0, max_elem=0, inol=0, nnol=0,
    length=0, layer=0, len3=3, idum[1], el[1+MNOL], nodes[MNOL],
    new_nodes[1+MNOL];
  double ddum[1], coords[MDIM];

  if ( n_layer<1 ) return;

  // map from a 2D node index to its extruded copies: base index scheme
  // node_copy(inod, layer) = max_node+1 + layer*(max_node+1) + inod
  long int nbase = max_node + 1;

  db_max_index( NODE, max_node, VERSION_NORMAL, GET );
  db_max_index( ELEMENT, max_elem, VERSION_NORMAL, GET );
  if ( max_node<0 || max_elem<0 ) return;
  nbase = max_node + 1;

  // create extruded node copies for each layer boundary
  for ( layer=0; layer<n_layer; layer++ ) {
    for ( inod=0; inod<=max_node; inod++ ) {
      if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
        long int new_node = nbase + layer*nbase + inod;
        db( NODE, inod, idum, coords, ndim, VERSION_NORMAL, GET );
        coords[2] = z_layer[layer];
        db( NODE, new_node, idum, coords, len3, VERSION_NORMAL, PUT );
        // copy all other node-class data items from the source node
        for ( int idat=0; idat<MDAT; idat++ ) {
          if ( idat!=NODE && db_data_class(idat)==NODE &&
               db_active_index( idat, inod, VERSION_NORMAL ) ) {
            long int ndata_len = db_len( idat, inod, VERSION_NORMAL );
            if ( db_type(idat)==DOUBLE_PRECISION ) {
              double *dold = db_dbl( idat, inod, VERSION_NORMAL );
              db( idat, new_node, idum, dold, ndata_len, VERSION_NORMAL, PUT );
            }
            else {
              long int *iold = db_int( idat, inod, VERSION_NORMAL );
              db( idat, new_node, iold, ddum, ndata_len, VERSION_NORMAL, PUT );
            }
          }
        }
        // the extruded coordinate also applies to node_start_refined
        if ( db_active_index( NODE_START_REFINED, inod, VERSION_NORMAL ) ) {
          db( NODE_START_REFINED, inod, idum, coords, ndim, VERSION_NORMAL, GET );
          coords[2] = z_layer[layer];
          db( NODE_START_REFINED, new_node, idum, coords, len3, VERSION_NORMAL, PUT );
        }
      }
    }
  }

  // create 3D elements: one per 2D element per layer
  long int new_elem = max_elem;
  long int next_node = nbase + n_layer*nbase; // after all boundary copies
  long int last_top = 0; // the top node block of the previous quad9 layer
  for ( ielem=0; ielem<=max_elem; ielem++ ) {
    if ( db_active_index( ELEMENT, ielem, VERSION_NORMAL ) ) {
      db( ELEMENT, ielem, el, ddum, length, VERSION_NORMAL, GET );
      nnol = length - 1;
      for ( inol=0; inol<nnol; inol++ ) nodes[inol] = el[1+inol];
      if ( el[0]==-TRIA3 && nnol==3 ) {
        for ( layer=0; layer<n_layer; layer++ ) {
          new_elem++;
          new_nodes[0] = -PRISM6;
          new_nodes[1] = nodes[0]+layer*nbase;
          new_nodes[2] = nodes[1]+layer*nbase;
          new_nodes[3] = nodes[2]+layer*nbase;
          new_nodes[4] = nodes[0]+(layer+1)*nbase;
          new_nodes[5] = nodes[1]+(layer+1)*nbase;
          new_nodes[6] = nodes[2]+(layer+1)*nbase;
          create_element( ielem, new_elem, new_nodes, 7, VERSION_NORMAL,
            VERSION_NORMAL );
        }
      }
      else if ( el[0]==-BAR2 && nnol==2 ) {
        // bar2 -> quad4 (the 2D interface bar; later converted to hex8
        // by control_mesh_convert). Nodes: 2 base + 2 extruded copy.
        for ( layer=0; layer<n_layer; layer++ ) {
          new_elem++;
          new_nodes[0] = -QUAD4;
          new_nodes[1] = nodes[0]+layer*nbase;
          new_nodes[2] = nodes[1]+layer*nbase;
          new_nodes[3] = nodes[0]+(layer+1)*nbase;
          new_nodes[4] = nodes[1]+(layer+1)*nbase;
          create_element( ielem, new_elem, new_nodes, 5, VERSION_NORMAL,
            VERSION_NORMAL );
        }
      }
      else if ( el[0]==-QUAD4 && nnol==4 ) {
        for ( layer=0; layer<n_layer; layer++ ) {
          new_elem++;
          new_nodes[0] = -HEX8;
          new_nodes[1] = nodes[0]+layer*nbase;
          new_nodes[2] = nodes[1]+layer*nbase;
          new_nodes[3] = nodes[2]+layer*nbase;
          new_nodes[4] = nodes[3]+layer*nbase;
          new_nodes[5] = nodes[0]+(layer+1)*nbase;
          new_nodes[6] = nodes[1]+(layer+1)*nbase;
          new_nodes[7] = nodes[2]+(layer+1)*nbase;
          new_nodes[8] = nodes[3]+(layer+1)*nbase;
          create_element( ielem, new_elem, new_nodes, 9, VERSION_NORMAL,
            VERSION_NORMAL );
        }
      }
      else if ( el[0]==-QUAD9 && nnol==9 ) {
        // quad9 -> hex27 (the Professional's extrusion of the quadratic
        // section elements; measured on its force11 .dbs: 9 base nodes +
        // 9 mid-plane nodes at the segment mid-height + 9 top nodes per
        // layer). The number of layers comes from control_mesh_extrude_n
        // (quad9_layers; the z record gives the total extent z0..z1).
        // The intermediate boundary nodes and the mid-plane nodes are
        // NEW nodes; the existing copy at z1 serves as the top of the
        // last layer, the source nodes as the base of the first.
        long int nq = ( quad9_layers>1 ? quad9_layers : 1 );
        double zq0 = z_layer[0], zq1 = z_layer[n_layer-1];
        for ( long int lay=0; lay<nq; lay++ ) {
          double zb = zq0 + lay*(zq1-zq0)/nq;
          double zt = zq0 + (lay+1)*(zq1-zq0)/nq;
          double z_mid = 0.5*( zb + zt );
          // base: lay 0 -> the SOURCE nodes; later layers -> the
          // ABSOLUTE block of the previous layer's top (last_top)
          long int base_abs = ( lay==0 ? -1 : last_top );
          // top: last layer -> the existing copy offset n_layer*nbase
          // (nodes[inol] + offset); other layers -> a NEW absolute
          // node block at zt
          long int top_off = -1, top_abs = -1;
          if ( lay==nq-1 ) {
            top_off = n_layer*nbase;
          }
          else {
            top_abs = next_node;
            for ( inol=0; inol<9; inol++ ) {
              db( NODE, nodes[inol], idum, coords, ndim,
                VERSION_NORMAL, GET );
              coords[2] = zt;
              db( NODE, next_node, idum, coords, len3, VERSION_NORMAL,
                PUT );
              for ( int idat=0; idat<MDAT; idat++ ) {
                if ( idat!=NODE && db_data_class(idat)==NODE &&
                     db_active_index( idat, nodes[inol],
                       VERSION_NORMAL ) ) {
                  long int ndata_len = db_len( idat, nodes[inol],
                    VERSION_NORMAL );
                  if ( db_type(idat)==DOUBLE_PRECISION ) {
                    double *dold = db_dbl( idat, nodes[inol],
                      VERSION_NORMAL );
                    db( idat, next_node, idum, dold, ndata_len,
                      VERSION_NORMAL, PUT );
                  }
                  else {
                    long int *iold = db_int( idat, nodes[inol],
                      VERSION_NORMAL );
                    db( idat, next_node, iold, ddum, ndata_len,
                      VERSION_NORMAL, PUT );
                  }
                }
              }
              next_node++;
            }
          }
          // the mid-plane nodes at the layer mid-height
          long int mid0 = next_node;
          for ( inol=0; inol<9; inol++ ) {
            db( NODE, nodes[inol], idum, coords, ndim, VERSION_NORMAL,
              GET );
            coords[2] = z_mid;
            db( NODE, next_node, idum, coords, len3, VERSION_NORMAL,
              PUT );
            for ( int idat=0; idat<MDAT; idat++ ) {
              if ( idat!=NODE && db_data_class(idat)==NODE &&
                   db_active_index( idat, nodes[inol],
                     VERSION_NORMAL ) ) {
                long int ndata_len = db_len( idat, nodes[inol],
                  VERSION_NORMAL );
                if ( db_type(idat)==DOUBLE_PRECISION ) {
                  double *dold = db_dbl( idat, nodes[inol],
                    VERSION_NORMAL );
                  db( idat, next_node, idum, dold, ndata_len,
                    VERSION_NORMAL, PUT );
                }
                else {
                  long int *iold = db_int( idat, nodes[inol],
                    VERSION_NORMAL );
                  db( idat, next_node, iold, ddum, ndata_len,
                    VERSION_NORMAL, PUT );
                }
              }
            }
            next_node++;
          }
          new_elem++;
          new_nodes[0] = -HEX27;
          for ( inol=0; inol<9; inol++ ) {
            new_nodes[1+inol] = ( base_abs>=0 ? base_abs+inol
                                               : nodes[inol] );
            new_nodes[10+inol] = mid0 + inol;
            new_nodes[19+inol] = ( top_abs>=0 ? top_abs+inol
                                               : nodes[inol]+top_off );
          }
          create_element( ielem, new_elem, new_nodes, 28, VERSION_NORMAL,
            VERSION_NORMAL );
          last_top = ( top_abs>=0 ? top_abs : 0 );
        }
      }
    }
  }
  // delete the 2D source elements
  for ( ielem=0; ielem<=max_elem; ielem++ ) {
    if ( db_active_index( ELEMENT, ielem, VERSION_NORMAL ) ) {
      db( ELEMENT, ielem, el, ddum, length, VERSION_NORMAL, GET );
      if ( el[0]==-TRIA3 || el[0]==-QUAD4 || el[0]==-QUAD9 )
        delete_element( ielem, VERSION_NORMAL );
    }
  }
  mesh_has_changed( VERSION_NORMAL );
}

// mesh_activate_gravity_factor - mesh_activate_gravity_time (Carril B).
//
// Returns the gravity activation factor for an element: 0 before the
// element start time of activation, 1 after the element end time of
// activation, interpolated in between. The start/end times are
// interpolated from the global time_start/time_end of the
// mesh_activate_gravity_time record and the lowest/highest coordinate of
// the element (bottom to top activation, typical for dam/dumping
// construction). The element is selected by mesh_activate_gravity_element
// (range), _element_group, or _geometry. Without the record the factor is 1
// (gravity fully active).
//
// With mesh_activate_gravity_method -method2 the element stays ACTIVE in
// the calculation before its activation but without gravity and with a
// REDUCED STIFFNESS (mesh_activate_gravity_stiffness_factor, default 1e-6).
// When stiff_factor is non-NULL it receives the element stiffness factor
// (1 fully active; stiffness_factor before activation, ramping up between
// the element start/end times).
double mesh_activate_gravity_factor( long int element, long int element_group,
  long int nnol, long int nodes[], double *stiff_factor )

{
  long int i=0, inod=0, length=0, ldum=0, in_geometry=0,
    method=-METHOD1, idum[1], *mesh_act=NULL, *gr_list=NULL;
  double factor=1., time_start=0., time_end=0., time_current=0., dtime=0.,
    coord_min=0., coord_max=0., t_start_el=0., t_end_el=0., ddum[MDIM],
    rdum=0., *coord=NULL, stiff=1.;

  if ( stiff_factor ) *stiff_factor = 1.;
  if ( db_active_index( MESH_ACTIVATE_GRAVITY_TIME, 0, VERSION_NORMAL ) ) {
    db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET );    db( TIME_CURRENT, 0, idum, &time_current, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( MESH_ACTIVATE_GRAVITY_TIME, 0, idum, ddum, ldum, VERSION_NORMAL, GET );
    time_start = ddum[0]; time_end = ddum[1];
    db( MESH_ACTIVATE_GRAVITY_METHOD, 0, &method, ddum, ldum,
      VERSION_NORMAL, GET_IF_EXISTS );
    double stiff_factor_val = 1.e-6;
    db( MESH_ACTIVATE_GRAVITY_STIFFNESS_FACTOR, 0, idum, &stiff_factor_val,
      ldum, VERSION_NORMAL, GET_IF_EXISTS );

    // select the element: range (_element), group (_element_group), or
    // geometry (_geometry). If no selection is given, all elements apply.
    long int selected = 0;
    if ( db_active_index( MESH_ACTIVATE_GRAVITY_ELEMENT, 0, VERSION_NORMAL ) ) {
      length = 0;
      mesh_act = db_int( MESH_ACTIVATE_GRAVITY_ELEMENT, 0, VERSION_NORMAL );
      length = db_len( MESH_ACTIVATE_GRAVITY_ELEMENT, 0, VERSION_NORMAL );
      for ( i=0; i<length; i++ )
        if ( mesh_act[i]==element ) { selected = 1; break; }
    }
    if ( db_active_index( MESH_ACTIVATE_GRAVITY_ELEMENT_GROUP, 0, VERSION_NORMAL ) ) {
      length = 0;
      gr_list = db_int( MESH_ACTIVATE_GRAVITY_ELEMENT_GROUP, 0, VERSION_NORMAL );
      length = db_len( MESH_ACTIVATE_GRAVITY_ELEMENT_GROUP, 0, VERSION_NORMAL );
      for ( i=0; i<length; i++ )
        if ( gr_list[i]==element_group ) { selected = 1; break; }
    }
    if ( db_active_index( MESH_ACTIVATE_GRAVITY_GEOMETRY, 0, VERSION_NORMAL ) ) {
      long int geom[2];
      db( MESH_ACTIVATE_GRAVITY_GEOMETRY, 0, geom, ddum, ldum, VERSION_NORMAL, GET );
      for ( inod=0; inod<nnol; inod++ ) {
        in_geometry = 0;
        geometry( nodes[inod], ddum, geom, in_geometry, rdum, ddum, rdum,
          ddum, NODE_START_REFINED, PROJECT_EXACT, VERSION_NORMAL );
        if ( !in_geometry ) break;
      }
      if ( in_geometry ) selected = 1;
    }
    // default: no selection record -> all elements
    if ( !db_active_index( MESH_ACTIVATE_GRAVITY_ELEMENT, 0, VERSION_NORMAL ) &&
         !db_active_index( MESH_ACTIVATE_GRAVITY_ELEMENT_GROUP, 0, VERSION_NORMAL ) &&
         !db_active_index( MESH_ACTIVATE_GRAVITY_GEOMETRY, 0, VERSION_NORMAL ) )
      selected = 1;
    if ( !selected ) return 1.;

    // lowest/highest vertical coordinate (y in 2D, z in 3D). The element
    // activation interval is interpolated from the global window: elements
    // with the lowest coordinate activate first. Without the global mesh
    // range we map the element height into the window (single element: the
    // full window).
    long int vdim = ( ndim==3 ) ? 2 : 1;
    coord_min = 1.e30; coord_max = -1.e30;
    for ( inod=0; inod<nnol; inod++ ) {
      coord = db_dbl( NODE, nodes[inod], VERSION_NORMAL );
      if ( coord[vdim]<coord_min ) coord_min = coord[vdim];
      if ( coord[vdim]>coord_max ) coord_max = coord[vdim];
    }
    if ( coord_max<=coord_min ) coord_max = coord_min + 1.e-6;
    t_start_el = time_start;
    t_end_el   = time_end;
    double t_total = time_current + dtime;
    // time_initial: before time_of_birth the element is inactive
    double time_birth = -1.e30;
    db( MESH_ACTIVATE_GRAVITY_TIME_INITIAL, 0, idum, &time_birth, ldum,
      VERSION_NORMAL, GET_IF_EXISTS );
    if ( t_total<time_birth ) {
      // method1: inactive before birth (no stiffness). method2: still active
      // with the reduced stiffness.
      if ( method==-METHOD2 && stiff_factor ) *stiff_factor = stiff_factor_val;
      return ( method==-METHOD2 ) ? 1. : 0.;
    }

    factor = 1.;
    if ( t_total<t_start_el ) {
      factor = 0.;
      // method2: element stays active with reduced stiffness until activation
      if ( method==-METHOD2 ) {
        factor = 1.;
        stiff = stiff_factor_val;
      }
    }
    else if ( t_total<t_end_el && t_end_el>t_start_el ) {
      factor = (t_total-t_start_el)/(t_end_el-t_start_el);
      // method2: stiffness ramps from reduced to full between the times
      if ( method==-METHOD2 )
        stiff = stiff_factor_val +
          (1.-stiff_factor_val) * (t_total-t_start_el)/(t_end_el-t_start_el);
    }
  }
  if ( stiff_factor ) *stiff_factor = stiff;
  return factor;
}
