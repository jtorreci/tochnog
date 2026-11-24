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

void spring( long int element, long int name, long int element_group, 
  double coord[], double old_dof[], double new_dof[], 
  double element_lhside[], double element_matrix[],
  double element_rhside[] )

{
  long int idim=0, jdim=0, ipuknwn=0, iuknwn=0, jpuknwn=0, juknwn=0, inol=0, jnol=0,
    indx=0, nnol=0, swit=0, length=0, ldum=0, 
    icontrol=0, options_convection=-NO, idum[1], options_mesh[MDIM];
  double dtime=0., group_spring_stiffness=0., 
    group_spring_plasti=1.e20, spring_force=0.,
    fac=0., tmp=0., old_length=0., new_length=0., incremental_length=0.,
    initial_length=0., spring_strain_old=0., spring_strain_total=0.,
    ddum[1], group_spring_direction[MDIM], initial_coord[MNOL*MDIM], 
    old_coord[MNOL*MDIM], new_coord[MNOL*MDIM], work[MDIM];

  swit = set_swit(element,-1,"spring");
  if ( swit ) pri( "In routine SPRING." );

  if ( db_active_index( GROUP_SPRING_STIFFNESS, element_group, VERSION_NORMAL ) ||
       db_active_index( GROUP_SPRING_STIFFNESS_NONLINEAR, element_group,
         VERSION_NORMAL ) ) {
    db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
    db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET );
    db( GROUP_SPRING_STIFFNESS, element_group, idum, &group_spring_stiffness, 
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( GROUP_SPRING_PLASTI, element_group, idum, &group_spring_plasti, 
      ldum, VERSION_NORMAL, GET_IF_EXISTS );

    // group_spring_stiffness_nonlinear (manual Professional 6.768):
    // diagram epsilon0 k0 epsilon1 k1 ... of the spring stiffness vs
    // the total spring strain (= total spring elongation). When the
    // record is present (and the linear stiffness is absent) the
    // stiffness of the increment is read from the diagram at the
    // MIDPOINT total strain of the increment; with a piecewise linear
    // diagram this integrates the spring force exactly along a
    // piecewise linear strain path (the corpus spring6 test checks
    // F = 0.5 EXACT for k(eps) = eps over an elongation of 1).
    double *nl_diagram = NULL;
    long int nl_length = 0;
    if ( db_active_index( GROUP_SPRING_STIFFNESS_NONLINEAR, element_group,
         VERSION_NORMAL ) ) {
      nl_diagram = db_dbl( GROUP_SPRING_STIFFNESS_NONLINEAR, element_group,
        VERSION_NORMAL );
      nl_length = db_len( GROUP_SPRING_STIFFNESS_NONLINEAR, element_group,
        VERSION_NORMAL );
      if ( nl_length<4 || (nl_length%2)!=0 )
        db_error( GROUP_SPRING_STIFFNESS_NONLINEAR, element_group );
    }

    if ( materi_velocity_integrated ) {
      options_convection = convection_apply_switch();
      db( OPTIONS_MESH, 0, options_mesh, ddum, ldum, VERSION_NORMAL, GET );
      for ( idim=0; idim<ndim; idim++ ) {
        if ( options_mesh[idim]==-FIXED_IN_SPACE ) {
          if ( options_convection==-YES ) {
            pri( "Error: set OPTIONS_CONVECTION to -NO in analysis with springs." );
            exit_tn_on_error();
          }
        }
      }
    }

    if       ( name==-SPRING1 )
      nnol = 1;
    else {
      assert ( name==-SPRING2 );
      nnol = 2;
    }

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

      // deformation (C88, F4 T6f 2026-10-19): the elongation is the
      // relative displacement PROJECTED on the spring direction,
      // group_spring_direction or, without the record, the INITIAL
      // axis of a -spring2 (first -> second node); a -spring1 measures
      // the displacement of its node from the initial position. The
      // force and the tangent act along the same fixed direction, so
      // the spring is linear and step-independent (Professional
      // 02-08-2026, measured for every group_spring_memory: an oblique
      // spring2 moved (0.3,0.1) along (1,1) gives k*0.4/sqrt(2) in 1 or
      // 10 steps, a spring2 rotated by 90 degrees gives 0, a spring1
      // moved (-0.2,0.5) along (1,0) gives -0.2*k). Before, the length
      // was the node distance (the full motion, unsigned for spring1).
    if ( !db( GROUP_SPRING_DIRECTION, element_group, idum,
         group_spring_direction, ldum, VERSION_NORMAL, GET_IF_EXISTS ) ) {
      if ( name==-SPRING1 )
        db_error( GROUP_SPRING_DIRECTION, element_group );
      else {
        assert ( name==-SPRING2 );
        array_subtract( &initial_coord[ndim], &initial_coord[0], group_spring_direction, ndim );
      }
    }
    if ( !array_normalize( group_spring_direction, ndim ) ) {
       pri( "Error: please specify a direction in GROUP_SPRING_DIRECTION." );
       exit(TN_EXIT_STATUS);
    }
    if ( name==-SPRING1 ) {
      array_subtract( &old_coord[0], &initial_coord[0], work, ndim );
      old_length = array_inproduct( work, group_spring_direction, ndim );
      array_subtract( &new_coord[0], &initial_coord[0], work, ndim );
      new_length = array_inproduct( work, group_spring_direction, ndim );
      initial_length = 0.;
    }
    else {
      array_subtract( &initial_coord[ndim], &initial_coord[0], work, ndim );
      initial_length = array_inproduct( work, group_spring_direction, ndim );
      array_subtract( &old_coord[ndim], &old_coord[0], work, ndim );
      old_length = array_inproduct( work, group_spring_direction, ndim );
      array_subtract( &new_coord[ndim], &new_coord[0], work, ndim );
      new_length = array_inproduct( work, group_spring_direction, ndim );
    }
    incremental_length = new_length - old_length;

      // spring force
    db( ELEMENT_SPRING_FORCE, element, idum, &spring_force, 
      length, VERSION_NORMAL, GET_IF_EXISTS );
    if ( nl_diagram ) {
      // total spring strain (= total elongation, manual Professional
      // 6.768) at the start and the end of the increment; the stiffness
      // of the increment is read at the midpoint strain so a piecewise
      // linear diagram is integrated exactly along a linear strain path.
      spring_strain_old = old_length - initial_length;
      spring_strain_total = new_length - initial_length;
      if ( nl_length>=4 ) {
        long int ipair=0;
        double eps=0., k_here=0.;
        eps = 0.5 * ( spring_strain_old + spring_strain_total );
        if ( eps<=nl_diagram[0] ) {
          group_spring_stiffness = nl_diagram[1];
        }
        else {
          for ( ipair=2; ipair+1<nl_length; ipair+=2 ) {
            if ( eps<nl_diagram[ipair] ) {
              k_here = nl_diagram[ipair-1] +
                ( eps - nl_diagram[ipair-2] ) *
                ( nl_diagram[ipair+1] - nl_diagram[ipair-1] ) /
                ( nl_diagram[ipair] - nl_diagram[ipair-2] );
              group_spring_stiffness = k_here;
              break;
            }
          }
          if ( ipair+1>=nl_length )
            group_spring_stiffness = nl_diagram[nl_length-1];
        }
      }
      length = 1;
      db( ELEMENT_SPRING_STRAIN, element, idum, &spring_strain_total,
        length, VERSION_NEW, PUT );
    }
    spring_force += group_spring_stiffness * incremental_length;
    if      (  spring_force>group_spring_plasti ) {
       spring_force = group_spring_plasti;
       group_spring_stiffness = 0.;
    }
    else if ( spring_force<-group_spring_plasti ) {
       spring_force = -group_spring_plasti;
       group_spring_stiffness = 0.;
    }

      // nodal force: -spring1 -F*dir at its node, -spring2 +F*dir at
      // the first and -F*dir at the second node; the tangent is
      // k*dtime*dir*dir^T (signed: an oblique direction couples the
      // components with the sign of dir_i*dir_j; C88: the former
      // |dir_i*dir_j| made a (1,-1) spring stiff along (1,1), and the
      // -spring1 branch wrote the right-hand side and the matrix at a
      // stale index with a negative stiffness)
    for ( idim=0; idim<ndim; idim++ ) {
      iuknwn = vel_indx + idim*nder;
      ipuknwn = iuknwn/nder;
      for ( inol=0; inol<nnol; inol++ ) {
        if ( name==-SPRING1 )
          fac = -1.;
        else if ( inol==0 )
          fac = +1.;
        else
          fac = -1.;
        indx = inol*npuknwn+ipuknwn;
        tmp = fac*spring_force*group_spring_direction[idim];
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
            tmp = fac*group_spring_stiffness*dtime*
              group_spring_direction[idim]*group_spring_direction[jdim];
            element_matrix[indx] += tmp;
            if ( jnol==inol && jdim==idim )
              element_lhside[inol*npuknwn+ipuknwn] += tmp;
          }
        }
      }
    }

    length = ndim;
    db( ELEMENT_SPRING_DIRECTION, element, idum, group_spring_direction, 
      length, VERSION_NEW, PUT );
    length = 1;
    db( ELEMENT_SPRING_FORCE, element, idum, &spring_force, 
      length, VERSION_NEW, PUT );

    if ( swit ) {
      pri( "initial_coord", initial_coord, nnol, ndim );
      pri( "new_coord", new_coord, nnol, ndim );
      pri( "group_spring_direction", group_spring_direction, ndim );
      pri( "old_length", old_length );
      pri( "new_length", new_length );
      pri( "incremental_length", incremental_length );
      pri( "spring_force", spring_force );
      pri( "element_lhside", element_lhside, nnol, npuknwn );
      pri( "element_rhside", element_rhside, nnol, npuknwn );
    }

  }

  if ( swit ) pri( "Out function SPRING" );
}
