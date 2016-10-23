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

void truss( long int element, long int element_group, 
  double coord[], double old_dof[], double new_dof[], 
  double element_lhside[], double element_matrix[],
  double element_rhside[] )

{
  long int idim=0, jdim=0, ipuknwn=0, iuknwn=0, jpuknwn=0, juknwn=0,
    inol=0, jnol=0, indx=0, nnol=2, swit=0, length=0,
    options_inertia=-YES, group_truss_rope=-NO, 
    memory=-UPDATED, options_convection=-YES,
    icontrol=0, ldum=0, idum[1], options_mesh[MDIM];
  double mass=0., dtime=0., 
    group_truss_young=0., group_truss_area=0., 
    group_truss_density=0., group_truss_plasti=1.e20,
    truss_stiffness=0., incremental_length=0., 
    old_truss_force=0., new_truss_force=0., fac=0., tmp=0., 
    old_length=0., new_length=0., initial_length=0., ddum[1],
    work[MDIM], truss_direction[MDIM], initial_coord[MNOL*MDIM], 
    old_coord[MNOL*MDIM], new_coord[MNOL*MDIM], diff_coord[MDIM],
    force_gravity[MDIM];

  swit = set_swit(element,-1,"truss");
  if ( swit ) pri( "In routine TRUSS." );

  if ( db_active_index( GROUP_TRUSS_AREA, element_group, VERSION_NORMAL ) ) {

    db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    if ( materi_velocity_integrated ) {
      db( OPTIONS_CONVECTION, 0, &options_convection, ddum,
        ldum, VERSION_NORMAL, GET_IF_EXISTS );
      db( CONTROL_OPTIONS_CONVECTION, icontrol, &options_convection, ddum,
        ldum, VERSION_NORMAL, GET_IF_EXISTS );
      db( OPTIONS_MESH, 0, options_mesh, ddum,
        ldum, VERSION_NORMAL, GET );
      for ( idim=0; idim<ndim; idim++ ) {
        if ( options_mesh[idim]==-FIXED_IN_SPACE ) {
          if ( options_convection==-YES ) {
            pri( "Error: set OPTIONS_CONVECTION to -NO in analysis with trusses." );
            exit_tn_on_error();
          }
        }
      }
    }

    db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET );
    db( GROUP_TRUSS_AREA, element_group, idum, &group_truss_area, 
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( GROUP_TRUSS_DENSITY, element_group, idum, &group_truss_density, 
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( GROUP_TRUSS_YOUNG, element_group, idum, &group_truss_young, 
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( GROUP_TRUSS_PLASTI, element_group, idum, &group_truss_plasti, 
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( GROUP_TRUSS_ROPE, element_group, &group_truss_rope, ddum,
       ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( GROUP_TRUSS_MEMORY, element_group, &memory, ddum,
       ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( OPTIONS_INERTIA, 0, &options_inertia, ddum, 
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( CONTROL_OPTIONS_INERTIA, icontrol, &options_inertia, ddum, 
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
    force_gravity_calculate( force_gravity );
      // geometry
    for ( idim=0; idim<ndim; idim++ ) {
      if      ( materi_velocity_integrated ) {
        for ( inol=0; inol<nnol; inol++ ) {
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

      // direction + deformations + stiffness
    initial_length = array_distance( &initial_coord[ndim], &initial_coord[0], work, ndim );
    if      ( memory==-UPDATED ) {
      array_subtract( &new_coord[ndim], &new_coord[0], truss_direction, ndim );
      if ( !array_normalize( truss_direction, ndim ) ) {
         pri( "Error: truss length zero for element ", element );
         exit(TN_EXIT_STATUS);
      }
        // deformation
      new_length = array_distance( &new_coord[0], &new_coord[ndim], work, ndim );
      old_length = array_distance( &old_coord[0], &old_coord[ndim], work, ndim );
      incremental_length = new_length - old_length;
      truss_stiffness = group_truss_young * group_truss_area / new_length;
    }
    else {
      array_subtract( &initial_coord[ndim], &initial_coord[0], truss_direction, ndim );
      if ( !array_normalize( truss_direction, ndim ) ) {
         pri( "Error: truss length zero for element ", element );
         exit(TN_EXIT_STATUS);
      }
        // deformation
      array_subtract( &new_coord[ndim], &new_coord[0], diff_coord, ndim );
      new_length = scalar_dabs( array_inproduct( diff_coord, truss_direction, ndim ) );
      array_subtract( &old_coord[ndim], &old_coord[0], diff_coord, ndim );
      old_length = scalar_dabs( array_inproduct( diff_coord, truss_direction, ndim ) );
      incremental_length = new_length - old_length;
        // stiffness
      truss_stiffness = group_truss_young * group_truss_area / initial_length;
    }

      // mass
    mass = initial_length * group_truss_area * group_truss_density;

      // forces
    db( ELEMENT_TRUSS_FORCE, element, idum, &old_truss_force, 
      length, VERSION_NORMAL, GET_IF_EXISTS );
    new_truss_force = old_truss_force + truss_stiffness * incremental_length;
    if      ( new_truss_force>group_truss_plasti*group_truss_area )
      new_truss_force = group_truss_plasti*group_truss_area;
    else if ( new_truss_force<-group_truss_plasti*group_truss_area )
      new_truss_force = -group_truss_plasti*group_truss_area;
    if ( group_truss_rope==-YES && new_truss_force<0. ) new_truss_force = 0.;
    for ( idim=0; idim<ndim; idim++ ) {
      iuknwn = vel_indx + idim*nder;
      ipuknwn = iuknwn/nder;
      for ( inol=0; inol<nnol; inol++ ) {
        if ( inol==0 )
          fac = +1.;
        else
          fac = -1.;
        indx = inol*npuknwn+ipuknwn;
        tmp = fac*truss_direction[idim]*new_truss_force;
        element_rhside[indx] += tmp + (0.5*mass*force_gravity[idim]);
        if ( options_inertia==-YES ) {
          element_rhside[indx] += - (0.5*mass)*
            (new_dof[inol*nuknwn+iuknwn]-old_dof[inol*nuknwn+iuknwn])/dtime;
        }
        for ( jdim=0; jdim<ndim; jdim++ ) {
          juknwn = vel_indx + jdim*nder;
          jpuknwn = juknwn/nder;
          for ( jnol=0; jnol<nnol; jnol++ ) {
            if ( jnol==inol ) {
              fac = +1.;
            }
            else {
              fac = -1.;
            }
            indx = inol*npuknwn*nnol*npuknwn+ipuknwn*nnol*npuknwn+
              jnol*npuknwn+jpuknwn;
            tmp = fac*truss_stiffness*dtime*
              truss_direction[idim]*truss_direction[jdim];
            if ( options_inertia==-YES ) {
              if ( jdim==idim && jnol==inol ) tmp += (mass/2) / dtime;
            }
            element_matrix[indx] += tmp;
            if ( jnol==inol && jdim==idim ) 
              element_lhside[inol*npuknwn+ipuknwn] += tmp;
          }
        }
      }
    }

    length = ndim;
    db( ELEMENT_TRUSS_DIRECTION, element, idum, truss_direction, 
      length, VERSION_NEW, PUT );
    length = 1;
    db( ELEMENT_TRUSS_FORCE, element, idum, &new_truss_force, 
      length, VERSION_NEW, PUT );
    if ( swit ) {
      pri( "coord", coord, nnol, ndim );
      pri( "initial_coord", initial_coord, nnol, ndim );
      pri( "old_coord", old_coord, nnol, ndim );
      pri( "new_coord", new_coord, nnol, ndim );
      pri( "truss_direction", truss_direction, ndim );
      pri( "initial_length", initial_length );
      pri( "old_length", old_length );
      pri( "new_length", new_length );
      pri( "incremental_length", incremental_length );
      pri( "old_truss_force", old_truss_force );
      pri( "new_truss_force", new_truss_force );
      pri( "element_matrix", element_matrix, nnol*npuknwn, nnol*npuknwn );
      pri( "element_lhside", element_lhside, nnol, npuknwn );
      pri( "element_rhside", element_rhside, nnol, npuknwn );
    }

  }

  if ( swit ) pri( "Out function TRUSS" );
}
