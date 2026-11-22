/*
    Copyright (C) 2000  Dennis Roddeman
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

// Direction from the side of the first spring node to the side of the
// second one: centroid of an element of node 1 minus centroid of an element
// of node 0, taking the first element of each node that is not a contact
// spring (the spring's own centroid is degenerate - both nodes coincide);
// a node without another element uses the spring itself.
static void contactspring_side_direction( long int nodes[], double coord[],
  double direction[] )

{
  long int idim=0;
  double cen[2][MDIM], ddum[1];

  for ( long int inol=0; inol<2; inol++ ) {
    array_set( cen[inol], 0., MDIM );
    for ( idim=0; idim<ndim; idim++ ) cen[inol][idim] = coord[inol*ndim+idim];
    long int *nel = db_int( NODE_ELEMENT, nodes[inol], VERSION_NORMAL );
    long int l = db_len( NODE_ELEMENT, nodes[inol], VERSION_NORMAL );
    if ( l<=0 ) continue;
    long int el = -1;
    long int *buf = get_new_int(MAXIMUM_NODE+1);
    for ( long int kk=0; kk<l && el<0; kk++ ) {
      long int lcand = 0;
      db( ELEMENT, nel[kk], buf, ddum, lcand, VERSION_NORMAL, GET );
      if ( buf[0]!=-CONTACTSPRING ) el = nel[kk];
    }
    if ( el<0 ) el = nel[0];
    long int le = 0, n = 0;
    db( ELEMENT, el, buf, ddum, le, VERSION_NORMAL, GET );
    array_set( cen[inol], 0., MDIM );
    for ( long int k=1; k<le; k++ ) {
      double *cn = db_dbl( NODE, buf[k], VERSION_NORMAL );
      for ( idim=0; idim<ndim; idim++ ) cen[inol][idim] += cn[idim];
      n++;
    }
    if ( n>0 ) for ( idim=0; idim<ndim; idim++ ) cen[inol][idim] /= n;
    delete[] buf;
  }
  array_set( direction, 0., MDIM );
  for ( idim=0; idim<ndim; idim++ ) direction[idim] = cen[1][idim] - cen[0][idim];
}

void contactspring( long int element, long int /*name*/, long int element_group, 
  long int nnol, long int nodes[], double coord[], double old_dof[], double new_dof[], 
  double element_lhside[], double element_matrix[],
  double element_rhside[] )

{
  long int idim=0, jdim=0, kdim=0, ipuknwn=0, iuknwn=0, jpuknwn=0, juknwn=0,
    inol=0, jnol=0, indx=0, inod=0, swit=0, length=0, memory=-UPDATED,
    group_contactspring_friction_automatic=-NO, iel=0, elnum=0, nel=0, gr=0,
    found=0, friction_unlimited=0, icontrol=0, ldum=0, options_convection=-NO, 
    options_skip_plasticity=-NO, idum[1], options_mesh[MDIM], *node_element=NULL;
  double k_n=0., k_t=0., f=0., cohesion=1.e20, dtime=0., 
    fac=0., gamma=0., phi_approximate=0.,
    tmp=0., total_tangential_force=0., group_contactspring_friction=0., 
    stiffness[MDIM], old_length[MDIM], new_length[MDIM], incremental_length[MDIM],
    group_contactspring_stiffness[2], group_contactspring_direction[MDIM], 
    element_contactspring_force[MDIM], element_contactspring_direction[MDIM*MDIM],
    initial_coord[MNOL*MDIM], old_coord[MNOL*MDIM], new_coord[MNOL*MDIM], 
    old_coord_diff[MDIM], new_coord_diff[MDIM], work[MDIM],
    group_materi_plasti_model[DATA_ITEM_SIZE], ddum[1];

  swit = set_swit(element,-1,"contactspring");
  if ( swit ) pri( "In routine CONTACTSPRING." );

  db( OPTIONS_SKIP_PLASTICITY, 0, &options_skip_plasticity, ddum,
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_OPTIONS_SKIP_PLASTICITY, icontrol, &options_skip_plasticity, ddum,
    ldum, VERSION_NORMAL, GET_IF_EXISTS );

  array_set( options_mesh, -FOLLOW_MATERIAL, ndim );
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  options_convection = convection_apply_switch();
  db( OPTIONS_MESH, 0, options_mesh, ddum,
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  // The OPTIONS_CONVECTION -no restriction applies ONLY to the
  // single-node surface contact (nnol==1), where the spring node sits on
  // a solid element and the mesh must stay fixed. The 2-node
  // contact_spring2 (a spring BETWEEN two nodes, nnol==2) does not need
  // it - it transmits force between the two nodes like a spring2.
  if ( nnol==1 ) {
    for ( idim=0; idim<ndim; idim++ ) {
      if ( options_mesh[idim]==-FIXED_IN_SPACE ) {
        if ( options_convection==-YES ) {
          pri( "Error: set OPTIONS_CONVECTION to -NO in analysis with contact springs." );
          exit_tn_on_error();
        }
      }
    }
  }


  db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET );
  db( GROUP_CONTACTSPRING_COHESION, element_group, idum, &cohesion, 
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( GROUP_CONTACTSPRING_STIFFNESS, element_group, idum, group_contactspring_stiffness, 
    ldum, VERSION_NORMAL, GET );
  db( GROUP_CONTACTSPRING_FRICTION, element_group, idum, &group_contactspring_friction, 
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( GROUP_CONTACTSPRING_FRICTION_AUTOMATIC, element_group,
    &group_contactspring_friction_automatic, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  k_n = group_contactspring_stiffness[0];
  k_t = group_contactspring_stiffness[1];
  stiffness[0] = k_n;
  stiffness[1] = k_t;
  stiffness[2] = k_t;

  if ( db( GROUP_CONTACTSPRING_MEMORY, element_group, &memory, ddum,
      ldum, VERSION_NORMAL, GET_IF_EXISTS ) ) {
    if ( memory!=-UPDATED && memory!=-UPDATED_WITHOUT_ROTATION
         && memory!=-TOTAL_LINEAR )
      db_error( GROUP_CONTACTSPRING_MEMORY, element_group );
  }

  if ( group_contactspring_friction_automatic==-YES ) {
    for ( inol=0; inol<nnol; inol++ ) {
      inod = nodes[inol];
      node_element = db_int( NODE_ELEMENT, inod, VERSION_NORMAL );
      nel = db_len( NODE_ELEMENT, inod, VERSION_NORMAL );
      for ( iel=0; iel<nel; iel++ ) {
        elnum = node_element[iel];
        gr = 0;
        db( ELEMENT_GROUP, elnum, &gr, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
        if      ( db_active_index( GROUP_MATERI_PLASTI_MOHRCOUL, gr, VERSION_NORMAL ) ) {
          found = 1;
          db( GROUP_MATERI_PLASTI_MOHRCOUL, gr, idum, group_materi_plasti_model, 
            ldum, VERSION_NORMAL, GET );
          f = BOUNDARY_REDUCTION_FACTOR * tan( group_materi_plasti_model[0] );
        }
        else if ( db_active_index( GROUP_MATERI_PLASTI_MATSUOKANAKAI, gr, VERSION_NORMAL ) ) {
          found = 1;
          db( GROUP_MATERI_PLASTI_MATSUOKANAKAI, gr, idum, group_materi_plasti_model, 
            ldum, VERSION_NORMAL, GET );
          f = BOUNDARY_REDUCTION_FACTOR * tan( group_materi_plasti_model[0] );
        }
        else if ( db_active_index( GROUP_MATERI_PLASTI_DIPRISCO, gr, VERSION_NORMAL ) ) {
          found = 1;
          db( GROUP_MATERI_PLASTI_DIPRISCO, gr, idum, group_materi_plasti_model, 
            ldum, VERSION_NORMAL, GET );
          gamma = group_materi_plasti_model[0];
          if      ( gamma>0. )
            phi_approximate = 1. * PIRAD / 180.;
          else if ( gamma>1. )
            phi_approximate = 1. * PIRAD / 180.;
          else
            phi_approximate = 1. * PIRAD / 180.;
          f = BOUNDARY_REDUCTION_FACTOR * tan( phi_approximate );
        }
        else if ( db_active_index( GROUP_MATERI_PLASTI_DIPRISCO_DENSITY, gr, VERSION_NORMAL ) ) {
          found = 1;
          phi_approximate = 1. * PIRAD / 180.;
          f = BOUNDARY_REDUCTION_FACTOR * tan( phi_approximate );
        }
      }
    }
    if ( !found ) f = 0.2;
  }
  else if ( db_active_index( GROUP_CONTACTSPRING_FRICTION, element_group,
      VERSION_NORMAL ) ) {
    f = group_contactspring_friction;
  }
  else {
    // no friction record: "a very large value for f will be applied"
    // (manual 6.604), so no frictional slip at all. Professional
    // 02-08-2026, conspr7 with a sideways load: tangential force -0.5 in
    // each spring, also in the one with normal force 0; disx 8.5. GNU
    // used f=0, which removed the tangential force and stiffness (disx
    // 1.5e12).
    friction_unlimited = 1;
  }
  if ( swit ) pri( "friction f", f );

    // geometry
  for ( idim=0; idim<ndim; idim++ ) {
    if      ( materi_velocity_integrated ) {
      for ( inol=0;inol<nnol; inol++ ) {
        initial_coord[inol*ndim+idim] = 
          coord[inol*ndim+idim];
        if ( options_mesh[idim]==-FOLLOW_MATERIAL )
          initial_coord[inol*ndim+idim] -= 
          old_dof[inol*nuknwn+veli_indx+idim*nder];
        old_coord[inol*ndim+idim] = 
          initial_coord[inol*ndim+idim] + 
          old_dof[inol*nuknwn+veli_indx+idim*nder];
        new_coord[inol*ndim+idim] = 
          old_coord[inol*ndim+idim] + 
          new_dof[inol*nuknwn+vel_indx+idim*nder]*dtime;
      }
    }
    else {
      assert( materi_displacement );
      for ( inol=0;inol<nnol; inol++ ) {
        initial_coord[inol*ndim+idim] = 
          coord[inol*ndim+idim];
        old_coord[inol*ndim+idim] = 
          initial_coord[inol*ndim+idim] + 
          old_dof[inol*nuknwn+dis_indx+idim*nder];
        new_coord[inol*ndim+idim] = 
          old_coord[inol*ndim+idim] + 
          new_dof[inol*nuknwn+vel_indx+idim*nder]*dtime;
      }
    }
  }

    // normal direction
  long int direction_automatic = -NO;
  db( GROUP_CONTACTSPRING_DIRECTION_AUTOMATIC, element_group, &direction_automatic,
    ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  array_set( element_contactspring_direction, 0., MDIM*MDIM );
  long int user_direction = 1;
  if ( direction_automatic==-YES ||
       !db( GROUP_CONTACTSPRING_DIRECTION, element_group, idum,
         group_contactspring_direction, ldum, VERSION_NORMAL, GET_IF_EXISTS ) ) {
    // automatic: the normal points FROM the first node's side TOWARDS the
    // second node's side, computed from the centroids of the elements the
    // two spring nodes belong to. This is stable through the solver
    // iterations: the spring nodes themselves coincide, and their raw
    // coordinate difference degenerates as soon as the mesh coordinates
    // update (second iteration gives a spurious (0,1e-6) that normalizes
    // to the WRONG sign). Verified against conspr7: node 5 belongs to the
    // upper block (centroid y=1.5), node 3 to the lower block (centroid
    // y=0.5); centroid(node1) - centroid(node0) = (0,-1) -> force = +1.0.
    user_direction = 0;
    array_set( group_contactspring_direction, 0., MDIM );
    if ( ndim>=2 ) {
       // direction = centroid(node1 element) - centroid(node0 element)
       // (the Professional convention: the normal points FROM the first
       // node's side TOWARDS the second node's side; conspr7 force = +1.0)
       contactspring_side_direction( nodes, coord, group_contactspring_direction );
       // group_contact_spring_direction_automatic_planes (manual 6.606):
       // a -no switch removes that component, then the direction is
       // normalized again. Professional 02-08-2026: centroid direction
       // (-0.447,-0.894) gives (0,-1) with switch_x -no and (-1,0) with
       // switch_y -no; nothing left is an input error.
       long int planes[MDIM];
       if ( db( GROUP_CONTACTSPRING_DIRECTION_AUTOMATIC_PLANES, element_group,
           planes, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS ) ) {
         for ( idim=0; idim<ndim; idim++ )
           if ( planes[idim]==-NO ) group_contactspring_direction[idim] = 0.;
         if ( !array_normalize( group_contactspring_direction, MDIM ) ) {
           pri( "Error: please specify a valid contact spring direction for element ", element );
           pri( "(group_contact_spring_direction_automatic_planes removes every component)" );
           exit(TN_EXIT_STATUS);
         }
       }
       if ( !array_normalize( group_contactspring_direction, MDIM ) )
         group_contactspring_direction[0] = 1.;
     }
     else
       group_contactspring_direction[0] = 1.;
  }
  for ( idim=ndim; idim<MDIM; idim++ ) {
    if ( group_contactspring_direction[idim]!=0. ) {
       pri( "Error: please specify a valid direction in GROUP_CONTACTSPRING_DIRECTION." );
       exit(TN_EXIT_STATUS);
    }
  }
  if ( user_direction ) {
    // user direction: normalised, and turned to point from the side of the
    // first node to the side of the second one, so that FN is positive in
    // extension (manual 6.603). Professional 02-08-2026: conspr6 with
    // 1 1 1 gives -0.57735 (x3), 0 0 1 and 0 0 -1 give the same normal;
    // conspr7 with 0 2 0 gives 0 -1 and force 1.
    if ( ndim>=2 ) {
      if ( !array_normalize( group_contactspring_direction, MDIM ) ) {
        pri( "Error: please specify a valid direction in GROUP_CONTACTSPRING_DIRECTION." );
        exit(TN_EXIT_STATUS);
      }
      double side[MDIM];
      contactspring_side_direction( nodes, coord, side );
      if ( array_inproduct( group_contactspring_direction, side, MDIM )<0. )
        array_multiply( group_contactspring_direction, group_contactspring_direction, -1., MDIM );
    }
  }
  array_move( group_contactspring_direction, 
    &element_contactspring_direction[0*MDIM], MDIM );

    // first tangential direction
  if ( ndim==2 ) {
    // 2D: the first tangential direction lies in the plane, the normal
    // turned a quarter: t1 = (-n_y, n_x) (Professional 02-08-2026,
    // element_contact_spring_direction: n (0,-1) gives t1 (1,0), n (-1,0)
    // gives t1 (0,-1)). The 3D rule below gave t1 = (0,0,-1) in 2D, so
    // the in-plane tangential stiffness was lost (a sideways load on
    // conspr7 gave disx 1.5e12; the Professional gives 8.5).
    element_contactspring_direction[1*MDIM+0] =
      -element_contactspring_direction[0*MDIM+1];
    element_contactspring_direction[1*MDIM+1] =
      element_contactspring_direction[0*MDIM+0];
    element_contactspring_direction[1*MDIM+2] = 0.;
  }
  else {
    // 3D: t1 = e_z x n, the 2D rule above extended (Professional
    // 02-08-2026, conspr6: n (-1,0,0) gives t1 (0,-1,0), n -(1,1,1)/sqrt(3)
    // gives t1 (0.7071,-0.7071,0)); for a normal along z the fallback
    // below gives t1 = e_y x n (n (0,0,1) gives t1 (1,0,0)).
    work[0] = 0.; work[1] = 0.; work[2] = 1.;
    array_outproduct_3D( work, &element_contactspring_direction[0*MDIM],
      &element_contactspring_direction[1*MDIM] );
  }
  if ( !array_normalize( &element_contactspring_direction[1*MDIM], MDIM ) ) {
    work[0] = 0.; work[1] = 1.; work[2] = 0.;
    array_outproduct_3D( work, &element_contactspring_direction[0*MDIM],
      &element_contactspring_direction[1*MDIM] );
    array_normalize( &element_contactspring_direction[1*MDIM], MDIM );
  }

    // second tangential direction
  array_outproduct_3D( &element_contactspring_direction[0*MDIM], 
    &element_contactspring_direction[1*MDIM], &element_contactspring_direction[2*MDIM] );
  array_normalize( &element_contactspring_direction[2*MDIM], MDIM );

    // deformation
  array_set( old_coord_diff, 0., MDIM );
  array_set( new_coord_diff, 0., MDIM );
  array_subtract( &old_coord[1*ndim], &old_coord[0*ndim], 
    old_coord_diff, ndim );
  array_subtract( &new_coord[1*ndim], &new_coord[0*ndim], 
    new_coord_diff, ndim );
  for ( kdim=0; kdim<MDIM; kdim++ ) {
    old_length[kdim] = array_inproduct( old_coord_diff, 
      &element_contactspring_direction[kdim*MDIM], ndim );
    new_length[kdim] = array_inproduct( new_coord_diff, 
      &element_contactspring_direction[kdim*MDIM], ndim );
    incremental_length[kdim] = new_length[kdim] - old_length[kdim];
  }

    // elastic contactspring forces
  array_set( element_contactspring_force, 0., MDIM );
  db( ELEMENT_CONTACTSPRING_FORCE, element, idum, 
    element_contactspring_force, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  for ( kdim=0; kdim<MDIM; kdim++ ) {
    element_contactspring_force[kdim] += stiffness[kdim] * incremental_length[kdim];
  }

  if ( element_contactspring_force[0]>cohesion && options_skip_plasticity==-NO ) {

      // cohesion exceeded
    array_set( element_contactspring_force, 0., MDIM );
    array_set( stiffness, 0., MDIM );

  }
  else {

      // frictional slip test
    if ( options_skip_plasticity==-NO && !friction_unlimited ) {
      total_tangential_force = sqrt( scalar_square(element_contactspring_force[1]) +
        scalar_square(element_contactspring_force[2]) );
      if ( scalar_dabs(total_tangential_force)>f*scalar_dabs(element_contactspring_force[0]) ) {
        fac = f*scalar_dabs(element_contactspring_force[0]/total_tangential_force);
        element_contactspring_force[1] *= fac;
        element_contactspring_force[2] *= fac;
        stiffness[1] *= fac;
        stiffness[2] *= fac;
      }
    }

      // nodal force
    for ( kdim=0; kdim<MDIM; kdim++ ) {
      for ( idim=0; idim<ndim; idim++ ) {
        iuknwn = vel_indx + idim*nder;
        ipuknwn = iuknwn/nder;
        for ( inol=0; inol<nnol; inol++ ) {
          if ( inol==0 )
            fac = +1.;
          else
            fac = -1.;
          indx = inol*npuknwn+ipuknwn;
          tmp = fac*element_contactspring_force[kdim]*
            element_contactspring_direction[kdim*MDIM+idim];
          element_rhside[indx] += tmp;
          for ( jdim=0; jdim<ndim; jdim++ ) {
            juknwn = vel_indx + jdim*nder;
            jpuknwn = juknwn/nder;
            for ( jnol=0; jnol<nnol; jnol++ ) {
              if ( jnol==inol )
                fac = +1.;
              else
                fac = -1.;
              indx = inol*npuknwn*nnol*npuknwn+ipuknwn*nnol*npuknwn+
                jnol*npuknwn+jpuknwn;
              tmp = fac*stiffness[kdim]*dtime*
                element_contactspring_direction[kdim*MDIM+idim]*
                element_contactspring_direction[kdim*MDIM+jdim];
              element_matrix[indx] += tmp;
              if ( jnol==inol && jdim==idim ) 
                element_lhside[inol*npuknwn+ipuknwn] += tmp;
            }
          }
        }
      }
    }

  }

  length = MDIM*MDIM;
  db( ELEMENT_CONTACTSPRING_DIRECTION, element, idum, 
    element_contactspring_direction, length, VERSION_NEW, PUT );
  length = MDIM;
  db( ELEMENT_CONTACTSPRING_FORCE, element, idum, 
    element_contactspring_force, length, VERSION_NEW, PUT );

  if ( swit ) {
    pri( "initial_coord", initial_coord, nnol, ndim );
    pri( "new_coord", new_coord, nnol, ndim );
    pri( "group_contactspring_direction", group_contactspring_direction, ndim );
    pri( "element_contactspring_direction", element_contactspring_direction, MDIM, MDIM );
    pri( "stiffness", stiffness, MDIM );
    pri( "old_length", old_length, MDIM );
    pri( "new_length", new_length, MDIM );
    pri( "incremental_length", incremental_length, MDIM );
    pri( "element_contactspring_force", element_contactspring_force, MDIM );
    pri( "element_lhside", element_lhside, nnol, npuknwn );
    pri( "element_rhside", element_rhside, nnol, npuknwn );
  }


  if ( swit ) pri( "Out function CONTACTSPRING" );
}
