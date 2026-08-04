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
