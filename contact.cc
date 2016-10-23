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

#define SLIP 1
#define STICK 2

void parallel_contact( void )

{
  long int icont=0, inum=0, itar=0, idim=0, inol=0, inod=0,
    nnol=0, iside=0, nside=0, nnol_side=0, name=0, 
    max_node=0, max_element=0, swit=0, ipuknwn=0, ldum=0, 
    length=0, use_element_target=0, use_contact_geometry=0, 
    max_contact_geometry=0, max_total=0, inside_target=0,
    contact_geometry_switch=0, contact_stick=-YES, status=0, iteration=0, 
    number_of_side_boundary_nodes=0, iloop=0,
    nloop=0, idat=0, max=0, any_contact_data=0, ithread=0, 
    max_node_boundary=0, idum[1], contact_geometry[10], 
    *el=NULL, *nodes=NULL, *side_nodes=NULL, *next_of_loop=NULL;
  double penetration=0.,
    tmp=0., temp=0., pres=0., pressure_penalty=0., temperature_penalty=0., 
    velocity_penalty=0., dtime=0., normal_force=0., friction_force=0., 
    pressure_force=0., temperature_force=0., friction_energy=0., 
    contact_heat_generation=0., slip_size=0., contact_friction=0.,
    fac=0., contact_relaxation=1., rdum=0., 
    ddum[MNOL], *tar_coord=NULL, *node_force=NULL,
    *normal_dir=NULL, *vec1=NULL, *vec2=NULL, *average_tar_coord=NULL, 
    *average_tar_side_coord=NULL, *cont_coord=NULL, *tar_side_coord=NULL, 
    *cont_vel=NULL, *average_tar_vel=NULL, *tar_side_vel=NULL, 
    *slip_dir=NULL, *penetration_coord=NULL, *weight=NULL, 
    *node_lhside=NULL, *node_rhside=NULL, *new_node_dof=NULL;

  for ( idat=0; idat<MDAT; idat++ ) {
    db_highest_index( idat, max, VERSION_NORMAL );
    if ( db_data_class(idat)==CONTACT && max>=0 ) any_contact_data = 1;
  }

  db_max_index( NODE, max_node, VERSION_NORMAL, GET );
  if ( max_node >= 0 && any_contact_data ) {

    db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
    db_max_index( NODE_BOUNDARY, max_node_boundary, VERSION_NORMAL, GET );
    db_max_index( CONTACT_GEOMETRY, max_contact_geometry, VERSION_NORMAL, GET );
    max_total = 1;
    max_total += max_contact_geometry;
    if ( max_node_boundary>=0 ) max_total += max_element;

    el = get_new_int(MNOL+1);
    nodes = get_new_int(MNOL);
    side_nodes = get_new_int(MNOL);
    tar_coord = get_new_dbl( MNOL*MDIM );
    node_force = get_new_dbl( MDIM );
    normal_dir = get_new_dbl( MDIM );
    vec1 = get_new_dbl( MDIM );
    vec2 = get_new_dbl( MDIM );
    average_tar_coord = get_new_dbl( MDIM );
    average_tar_side_coord = get_new_dbl( MDIM );
    cont_coord = get_new_dbl( MDIM );
    tar_side_coord = get_new_dbl( MDIM*MNOL );
    cont_vel = get_new_dbl( MDIM );
    average_tar_vel = get_new_dbl( MDIM );
    tar_side_vel = get_new_dbl( MDIM );
    slip_dir = get_new_dbl( MDIM );
    penetration_coord = get_new_dbl( MDIM );
    weight = get_new_dbl( MDIM );

    area_node_dataitem();

    swit = set_swit(-1,-1,"contact");
    if ( swit ) pri( "In routine CONTACT" );

    db( NUMBER_ITERATIONS, 0, &iteration, ddum, ldum, VERSION_NEW, GET );
    db( CONTACT_HEATGENERATION, 0, idum, 
      &contact_heat_generation, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( CONTACT_PENALTY_PRESSURE, 0, idum, 
      &pressure_penalty, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( CONTACT_PENALTY_TEMPERATURE, 0, idum, 
      &temperature_penalty, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( CONTACT_PENALTY_VELOCITY, 0, idum, 
      &velocity_penalty, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( CONTACT_FRICTION, 0, idum, &contact_friction, 
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( CONTACT_STICK, 0, &contact_stick, ddum,
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( CONTACT_RELAXATION, 0, idum, &contact_relaxation, 
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET );
    if ( swit ) {
      pri( "pressure_penalty", pressure_penalty );
      pri( "temperature_penalty", temperature_penalty );
      pri( "velocity_penalty", velocity_penalty );
      pri( "contact_friction", contact_friction );
    }

      // loop over contacters
    if ( max_node>=0 ) {
      next_of_loop = get_new_int(1+max_node);
      parallel_sys_next_of_loop( next_of_loop, max_node, nloop, ithread );
      for ( iloop=0; iloop<nloop; iloop++ ) {
        icont = next_of_loop[iloop];
        if ( icont>max_node )
          break;
        else if ( db_active_index( NODE, icont, VERSION_NORMAL ) ) {
            // contacter coordinates
          db( NODE, icont, idum, cont_coord, ldum, VERSION_NEW, GET );
          new_node_dof = db_dbl( NODE_DOF, icont, VERSION_NEW );
          node_rhside = db_dbl( NODE_RHSIDE, icont, VERSION_NORMAL );
          for ( idim=0; idim<ndim; idim++ ) {
            if ( materi_displacement ) cont_coord[idim] += 
              new_node_dof[dis_indx+idim*nder];
              // explicit prediction
            if ( iteration==1 ) cont_coord[idim] += 
              new_node_dof[vel_indx+idim*nder]*dtime;
            cont_vel[idim] = new_node_dof[vel_indx+idim*nder];
          }
          if ( swit ) {
            pri( "node_contacter", icont );
            pri( "cont_coord", cont_coord, ndim );
            pri( "cont_vel", cont_vel, ndim );
          }
            // loop over target geometries and elements
          for ( inum=0; inum<=max_total; inum++ ) {
            inside_target = 0;
            if ( inum<=max_contact_geometry ) {
              use_contact_geometry = 1;
              use_element_target = 0;
              itar = inum;
            }
            else {
              use_contact_geometry = 0;
              use_element_target = 1;
              itar = inum - max_contact_geometry - 1;
            }
            if ( (use_element_target&&db_active_index(ELEMENT,itar, 
                  VERSION_NORMAL))||
                 (use_contact_geometry&&db_active_index(CONTACT_GEOMETRY,itar, 
                  VERSION_NORMAL)) ) {
              if ( use_element_target ) {
                if ( swit ) pri( "itar", itar );
                db( ELEMENT, itar, el, ddum, length, VERSION_NORMAL, GET );
                name = el[0]; nnol = length - 1; 
                if  ( name!=-BAR2 && name!=-TRIA3 && name!=-QUAD4 &&
                      name!=-TET4 && name!=-HEX8 ) {
                  cout << "\nError: " << db_name( name );
                  cout << " is not available for contact analysis.\n"; 
                  exit(TN_EXIT_STATUS);
                }
                array_move( &el[1], nodes, nnol );
                  // target coordinates and average
                array_set( average_tar_coord, 0., ndim );
                array_set( average_tar_vel, 0., ndim );
                for ( inol=0; inol<nnol; inol++ ) {
                  inod = nodes[inol];
                  new_node_dof = db_dbl( NODE_DOF, inod, VERSION_NEW );
                  db( NODE, inod, idum, &tar_coord[inol*ndim], 
                    ldum, VERSION_NEW, GET );
                  for ( idim=0; idim<ndim; idim++ ) {
                    if ( materi_displacement ) tar_coord[inol*ndim+idim] += 
                      new_node_dof[dis_indx+idim*nder];
                    if ( iteration==1 ) tar_coord[inol*ndim+idim] += 
                      new_node_dof[vel_indx+idim*nder]*dtime;
                    average_tar_coord[idim] += tar_coord[inol*ndim+idim]/nnol;
                    average_tar_vel[idim] += new_node_dof[vel_indx+idim*nder]/nnol;
                  }
                }
                  // inside target element?
                if ( !array_member(nodes,icont,nnol,ldum) &&
                    point_el( cont_coord, tar_coord, ddum, name, nnol ) ) {
                  if ( swit ) {
                    pri( "tar_coord", tar_coord, nnol, ndim );
                    pri( "average_tar_coord", average_tar_coord, ndim );
                    pri( "average_tar_vel", average_tar_vel, ndim );
                  }
                    // loop over target sides
                  if      ( name==-BAR2 ) {
                    nside = 2;
                    nnol_side = 1;
                  }
                  else if ( name==-TRIA3 ) {
                    nside = 3;
                    nnol_side = 2;
                  }
                  else if ( name==-QUAD4 ) {
                    nside = 4;
                    nnol_side = 2;
                  }
                  else if ( name==-TET4 ) {
                    nside = 4;
                    nnol_side = 3;
                  }
                  else if ( name==-HEX8 ) {
                    nside = 6;
                    nnol_side = 4;
                  }
                  else {
                    nside = 0;
                    nnol_side = 0;
                  }
                  for ( iside=0; iside<nside && !inside_target; iside++ ) {
                    if ( swit ) pri( "iside", iside );
                    if      ( name==-BAR2 ) {
                      if ( iside==0 )
                        side_nodes[0] = 0;
                      else {
                        assert( iside==1 );
                        side_nodes[0] = 1;
                      }
                    }
                    else if ( name==-TRIA3 ) {
                      if      ( iside==0 ) {
                        side_nodes[0] = 0;
                        side_nodes[1] = 1;
                      }
                      else if ( iside==1 ) {
                        side_nodes[0] = 1;
                        side_nodes[1] = 2;
                      }
                      else {
                        assert ( iside==2 );
                        side_nodes[0] = 2;
                        side_nodes[1] = 0;
                      }
                    }
                    else if ( name==-QUAD4 ) {
                      if      ( iside==0 ) {
                        side_nodes[0] = 0;
                        side_nodes[1] = 1;
                      }
                      else if ( iside==1 ) {
                        side_nodes[0] = 1;
                        side_nodes[1] = 3;
                      }
                      else if ( iside==2 ) {
                        side_nodes[0] = 3;
                        side_nodes[1] = 2;
                      }
                      else {
                        assert( iside==3 );
                        side_nodes[0] = 2;
                        side_nodes[1] = 0;
                      }
                    }
                    else if ( name==-TET4 ) {
                      if      ( iside==0 ) {
                        side_nodes[0] = 0;
                        side_nodes[1] = 1;
                        side_nodes[2] = 2;
                      }
                      else if ( iside==1 ) {
                        side_nodes[0] = 0;
                        side_nodes[1] = 1;
                        side_nodes[2] = 3;
                      }
                      else if ( iside==2 ) {
                        side_nodes[0] = 2;
                        side_nodes[1] = 3;
                        side_nodes[2] = 0;
                      }
                      else {
                        assert ( iside==3 );
                        side_nodes[0] = 2;
                        side_nodes[1] = 3;
                        side_nodes[2] = 1;
                      }
                    }
                    else if ( name==-HEX8 ) {
                      if      ( iside==0 ) {
                        side_nodes[0] = 0;
                        side_nodes[1] = 1;
                        side_nodes[2] = 2;
                        side_nodes[3] = 3;
                      }
                      else if ( iside==1 ) {
                        side_nodes[0] = 0;
                        side_nodes[1] = 1;
                        side_nodes[2] = 4;
                        side_nodes[3] = 5;
                      }
                      else if ( iside==2 ) {
                        side_nodes[0] = 1;
                        side_nodes[1] = 2;
                        side_nodes[2] = 5;
                        side_nodes[3] = 6;
                      }
                      else if ( iside==3 ) {
                        side_nodes[0] = 2;
                        side_nodes[1] = 3;
                        side_nodes[2] = 6;
                        side_nodes[3] = 7;
                      }
                      else if ( iside==4 ) {
                        side_nodes[0] = 3;
                        side_nodes[1] = 0;
                        side_nodes[2] = 7;
                        side_nodes[3] = 4;
                      }
                      else if ( iside==5 ) {
                        side_nodes[0] = 4;
                        side_nodes[1] = 5;
                        side_nodes[2] = 6;
                        side_nodes[3] = 7;
                      }
                    }
                    number_of_side_boundary_nodes = 0;
                    for ( inol=0; inol<nnol_side; inol++ ) {
                      inod = nodes[side_nodes[inol]];
                      if ( db_active_index( NODE_BOUNDARY, inod, VERSION_NORMAL ) )
                        number_of_side_boundary_nodes++;
                    }
                    if ( swit ) pri( "number_of_side_boundary_nodes",
                      number_of_side_boundary_nodes );
                    if ( number_of_side_boundary_nodes==ndim ) {
                      inside_target = 1;
                        // side properties
                      array_set( average_tar_side_coord, 0., ndim );
                      array_set( tar_side_vel, 0., ndim );
                      for ( inol=0; inol<nnol_side; inol++ ) {
                        inod = nodes[side_nodes[inol]];
                        new_node_dof = db_dbl( NODE_DOF, inod, VERSION_NEW );
                        db( NODE, inod, idum, &tar_side_coord[inol*ndim], ldum, 
                          VERSION_NEW, GET );
                        for ( idim=0; idim<ndim; idim++ ) {
                          if ( materi_displacement ) 
                            tar_side_coord[inol*ndim+idim] +=
                            new_node_dof[dis_indx+idim*nder];
                          if ( iteration==1 ) 
                            tar_side_coord[inol*ndim+idim] +=
                            new_node_dof[vel_indx+idim*nder]*dtime;
                          average_tar_side_coord[idim] += 
                            tar_side_coord[inol*ndim+idim]/nnol_side;
                          tar_side_vel[idim] += 
                            new_node_dof[vel_indx+idim*nder] / nnol_side;
                        }
                      }
                      if ( swit ) {
                        pri( "average_tar_side_coord", 
                          average_tar_side_coord, ndim );
                        pri( "tar_side_vel", tar_side_vel, ndim );
                      }
                      if      ( nnol_side==1 ) {
                        weight[0] = 1.;
                      }
                      else if ( nnol_side==2 ) {
                        project_point_exactly_on_line( cont_coord, 
                          &tar_side_coord[0*ndim], 
                          &tar_side_coord[1*ndim], weight );
                      }
                      else if ( nnol_side==3 ) {
                        project_point_exactly_on_triangle( cont_coord, 
                          &tar_side_coord[0*ndim], 
                          &tar_side_coord[1*ndim], 
                          &tar_side_coord[2*ndim], 
                          weight );
                      }
                      else {
                        assert( nnol_side==4 );
                        project_point_exactly_on_quad( cont_coord, 
                          &tar_side_coord[0*ndim], 
                          &tar_side_coord[1*ndim], 
                          &tar_side_coord[2*ndim], 
                          &tar_side_coord[3*ndim], 
                          weight );
                      }
                        // normal
                      if      ( nnol_side==1 ) {
                        normal_dir[0] = -1.;
                      }
                      else if ( nnol_side==2 ) {
                        array_subtract( &tar_side_coord[1*ndim], 
                          &tar_side_coord[0*ndim], vec1, ndim );
                        normal_dir[0] = -vec1[1];
                        normal_dir[1] = +vec1[0];
                      }
                      else if ( nnol_side==3 ) {
                        array_subtract( &tar_side_coord[1*ndim], 
                          &tar_side_coord[0*ndim], vec1, ndim );
                        array_subtract( &tar_side_coord[2*ndim], 
                          &tar_side_coord[0*ndim], vec2, ndim );
                        array_outproduct_3D( vec1, vec2, normal_dir );
                      }
                      else {
                        assert( nnol_side==4 );
                        array_subtract( &tar_side_coord[1*ndim], 
                          &tar_side_coord[0*ndim], vec1, ndim );
                        array_subtract( &tar_side_coord[2*ndim], 
                          &tar_side_coord[0*ndim], vec2, ndim );
                        array_outproduct_3D( vec1, vec2, normal_dir );
                      }
                        // check if normal is outward
                      array_subtract( average_tar_side_coord, 
                        average_tar_coord, vec1, ndim );
                      tmp = array_inproduct( vec1, normal_dir, ndim );
                      if ( tmp<0. ) 
                        array_multiply( normal_dir, normal_dir, -1., ndim );
                      array_normalize( normal_dir, ndim );
                        // penetration
                      array_set( penetration_coord, 0., ndim );
                      for ( inol=0; inol<nnol_side; inol++ ) {
                        for ( idim=0; idim<ndim; idim++ ) {
                          penetration_coord[idim] += weight[inol] *
                            tar_side_coord[inol*ndim+idim];
                        }
                      }
                      array_subtract( cont_coord, penetration_coord, vec1, ndim );
                      penetration = array_inproduct( vec1, normal_dir, ndim );
                      if ( swit ) {
                        pri( "normal", normal_dir, ndim );
                        pri( "penetration_coord", penetration_coord, ndim );
                        pri( "penetration vector", vec1, ndim );
                        pri( "penetration normal", normal_dir, ndim );
                        pri( "penetration", penetration );
                      }
                    }
                  }
                }
              }
              else {
                assert( use_contact_geometry );
                db( CONTACT_GEOMETRY, itar, contact_geometry, ddum, 
                  ldum, VERSION_NORMAL, GET );
                geometry( icont, ddum, contact_geometry, inside_target, rdum, 
                  normal_dir, penetration, ddum, PLUS_DISPLACEMENT, 
                  PROJECT_ON_EDGE, VERSION_NEW );
                if ( inside_target ) {
                  if ( db_active_index( CONTACT_GEOMETRY_SWITCH, 
                    itar, VERSION_NORMAL ) ) {
                    db( CONTACT_GEOMETRY_SWITCH, itar, 
                      &contact_geometry_switch, ddum, ldum, VERSION_NORMAL, GET );
                    if ( contact_geometry_switch==-YES ) {
                      array_multiply( normal_dir, normal_dir, -1., ndim );
                      penetration = -penetration;
                    }
                  }
                  if ( penetration<0. )
                    inside_target = 1;
                  else
                    inside_target = 0;            
                  if ( swit ) {
                    pri( "normal", normal_dir, ndim );
                    pri( "penetration", penetration );
                    pri( "inside_target", inside_target );
                  }
                }
              }
              if ( inside_target ) {
                if ( use_element_target ) {
                  fac = 2.;
                  if ( groundflow_pressure ) {
                    new_node_dof = db_dbl( NODE_DOF, icont, VERSION_NEW );
                    pres = -new_node_dof[pres_indx];
                    for ( inol=0; inol<nnol_side; inol++ ) {
                      inod = nodes[side_nodes[inol]];
                      new_node_dof = db_dbl( NODE_DOF, inod, VERSION_NEW );
                      pres += weight[inol] * new_node_dof[pres_indx];
                    } 
                    pressure_force = pressure_penalty * pres;
                  }
                  if ( condif_temperature ) {
                    new_node_dof = db_dbl( NODE_DOF, icont, VERSION_NEW );
                    temp = -new_node_dof[temp_indx];
                    for ( inol=0; inol<nnol_side; inol++ ) {
                      inod = nodes[side_nodes[inol]];
                      new_node_dof = db_dbl( NODE_DOF, inod, VERSION_NEW );
                      temp += weight[inol] * new_node_dof[temp_indx];
                    } 
                    temperature_force = temperature_penalty * temp;
                  }
                }
                else {
                  assert( use_contact_geometry );
                  fac = 1.;
                  array_set( tar_side_vel, 0., ndim );
                }
                  // slip direction
                array_subtract( tar_side_vel, cont_vel, slip_dir, ndim );
                tmp = array_inproduct( normal_dir, slip_dir, ndim );
                array_multiply( normal_dir, vec1, tmp, ndim );
                array_subtract( slip_dir, vec1, slip_dir, ndim );
                slip_size = array_size( slip_dir, ndim );
                array_normalize( slip_dir, ndim );
                  // penetration force
                normal_force = velocity_penalty * scalar_dabs(penetration);
                  // slip force
                for ( idim=0; idim<ndim; idim++ ) {
                  ipuknwn = vel_indx/nder + idim;
                  node_force[idim] = node_rhside[ipuknwn] * slip_dir[idim];
                }
                if ( contact_stick==-NO || 
                    array_size(node_force,ndim)>=0.5*contact_friction*normal_force ) {
                  status = SLIP;
                  friction_force = contact_friction * normal_force;
                  friction_energy = contact_heat_generation * friction_force * slip_size;
                }
                else {
                  status = STICK;
                  friction_force = velocity_penalty * slip_size * dtime;
                  friction_energy = 0.;
                }
                if ( swit ) {
                  pri( "normal_force", normal_force );
                  pri( "slip_dir", slip_dir, ndim );
                  pri( "node_force", node_force, ndim );
                  if ( status==STICK ) pri( "status is STICK" );
                  else pri( "status is SLIP" );
                  pri( "friction_force", friction_force );
                  pri( "friction_energy", friction_energy );
                  pri( "pressure_force", pressure_force );
                  pri( "temperature_force", temperature_force );
                }
                  // contacter contributions
                node_lhside = db_dbl( NODE_LHSIDE, icont, VERSION_NORMAL );
                node_rhside = db_dbl( NODE_RHSIDE, icont, VERSION_NORMAL );
                if ( groundflow_pressure ) {
                  ipuknwn = pres_indx/nder;
                  node_lhside[ipuknwn] += pressure_penalty * fac;
                  node_rhside[ipuknwn] += pressure_force;
                }
                if ( condif_temperature ) {
                  ipuknwn = temp_indx/nder;
                  node_lhside[ipuknwn] += temperature_penalty * fac;
                  node_rhside[ipuknwn] += temperature_force;
                  if ( use_contact_geometry ) 
                    node_rhside[ipuknwn] += friction_energy;
                  else 
                    node_rhside[ipuknwn] += 0.5 * friction_energy;
                }
                for ( idim=0; idim<ndim; idim++ ) {
                  ipuknwn = vel_indx/nder + idim;
                    // normal force on contacter
                  node_lhside[ipuknwn] += velocity_penalty *  contact_relaxation *
                    scalar_dabs(normal_dir[idim]) * dtime * fac;
                  node_rhside[ipuknwn] += normal_force * normal_dir[idim];
                    // friction force on contacter
                  if ( status==STICK ) {
                    node_lhside[ipuknwn] += velocity_penalty * 
                      scalar_dabs(slip_dir[idim]) * dtime * fac;
                    node_rhside[ipuknwn] += friction_force * slip_dir[idim];
                  }
                  else {
                    assert( status==SLIP );
                    node_rhside[ipuknwn] += friction_force * slip_dir[idim];
                  }
                }
                if ( swit ) {
                  pri( "cont_inod", icont );
                  pri( "node_lhside", node_lhside, npuknwn );
                  pri( "node_rhside", node_rhside, npuknwn );
                }
                  // target contributions
                if ( use_element_target ) {
                  parallel_sys_lock();
                  for ( inol=0; inol<nnol_side; inol++ ) {
                    inod = nodes[side_nodes[inol]];
                    node_lhside = db_dbl( NODE_LHSIDE, inod, VERSION_NORMAL );
                    node_rhside = db_dbl( NODE_RHSIDE, inod, VERSION_NORMAL );
                    if ( groundflow_pressure ) {
                      ipuknwn = pres_indx/nder;
                      node_lhside[ipuknwn] += pressure_penalty * weight[inol] * fac;
                      node_rhside[ipuknwn] -= pressure_force * weight[inol];
                    }
                    if ( condif_temperature ) {
                      ipuknwn = temp_indx/nder;
                      node_lhside[ipuknwn] += temperature_penalty * weight[inol] * fac;
                      node_rhside[ipuknwn] -= ( temperature_force +
                        0.5 * friction_energy ) * weight[inol];
                    }
                    for ( idim=0; idim<ndim; idim++ ) {
                      ipuknwn = vel_indx/nder + idim;
                        // normal force on target
                      node_lhside[ipuknwn] +=  contact_relaxation *
                        velocity_penalty * weight[inol] * 
                        scalar_dabs(normal_dir[idim]) * dtime * fac;
                      node_rhside[ipuknwn] -= 
                        normal_force * weight[inol] * normal_dir[idim];
                        // friction force on target
                      if ( status==STICK ) {
                        node_lhside[ipuknwn] += 
                          velocity_penalty * weight[inol] * 
                          scalar_dabs(slip_dir[idim]) * dtime * fac;
                        node_rhside[ipuknwn] -= 
                          friction_force * weight[inol] * slip_dir[idim];
                      }
                      else {
                        assert( status==SLIP );
                        node_rhside[ipuknwn] -= 
                          friction_force * weight[inol] * slip_dir[idim];
                      }
                    }
                    if ( swit ) {
                      pri( "tar_inod", inod );
                      pri( "node_lhside", node_lhside, npuknwn );
                      pri( "node_rhside", node_rhside, npuknwn );
                    }
                  }
                  parallel_sys_unlock();
                }
              }
            }
          }
        }
      }
      delete[] next_of_loop;
    }

    delete[] el;
    delete[] nodes;
    delete[] side_nodes;
    delete[] tar_coord;
    delete[] node_force;
    delete[] normal_dir;
    delete[] vec1;
    delete[] vec2;
    delete[] average_tar_coord;
    delete[] average_tar_side_coord;
    delete[] cont_coord;
    delete[] tar_side_coord;
    delete[] cont_vel;
    delete[] average_tar_vel;
    delete[] tar_side_vel;
    delete[] slip_dir;
    delete[] penetration_coord;
    delete[] weight;

    if ( swit ) pri( "Out routine CONTACT" );
  }

}
