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

void mesh_extrude( double z_layer[], long int n_layer )

{
  // extrude a 2D mesh (z=0) to 3D along the z-axis.
  // z_layer[] gives the z-coordinate of each layer boundary;
  // one 3D element is generated per 2D element per layer.
  // -tria3 -> -prism6, -quad4 -> -hex8.
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
    }
  }
  db_max_index( ELEMENT, max_elem, VERSION_NORMAL, GET );
  for ( ielem=0; ielem<=max_elem; ielem++ ) {
    if ( db_active_index( ELEMENT, ielem, VERSION_NORMAL ) ) {
      db( ELEMENT, ielem, el, ddum, length, VERSION_NORMAL, GET );
      for ( int kk=0; kk<length; kk++ ) cout << " " << el[kk];
      cout << endl;
    }
  }
  // delete the 2D source elements
  for ( ielem=0; ielem<=max_elem; ielem++ ) {
    if ( db_active_index( ELEMENT, ielem, VERSION_NORMAL ) ) {
      db( ELEMENT, ielem, el, ddum, length, VERSION_NORMAL, GET );
      if ( el[0]==-TRIA3 || el[0]==-QUAD4 )
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
  long int i=0, idim=0, inod=0, icontrol=0, length=0, ldum=0, in_geometry=0,
    found=0, method=-METHOD1, idum[1], *mesh_act=NULL, *gr_list=NULL;
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
