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

#define NDOF 3  // number of dof's 2d-beam per node
#define NNOL 2  // number of nodes 2d-beam
#define NDIM 2  // number of dimensions 2d-beam

void beam_3d( long int element, long int element_group, double coord_3d[], 
  double old_dof_3d[], double new_dof_3d[], double element_lhside_3d[], 
  double element_matrix_3d[], double element_rhside_3d[] ) 

// Beam element: the classic GNU beam in the x-y plane in 2D, the beam with
// local beam axes (beam_local) when a local z record applies, and in 3D
// always the native 3D beam (beam_local) with default local axes (F5 T2
// V4d, owner decision C37: bending about both local axes and torsion,
// as the Professional; the former 3D beam mapped the element onto the
// 2D beam in the group_beam_plane plane).

{
  double record_vec[3];
  long int record = beam_direction_z_record( element, element_group,
    record_vec );

  if ( record || ndim==3 ) {
    beam_local( element, element_group, record, record_vec, coord_3d,
      old_dof_3d, new_dof_3d, element_lhside_3d, element_matrix_3d,
      element_rhside_3d );
    return;
  }

  beam_2d( element, element_group, coord_3d, old_dof_3d, new_dof_3d, 
    element_lhside_3d, element_matrix_3d, element_rhside_3d );
}


void beam_2d( long int element, long int element_group, double coord[], 
  double old_dof[], double new_dof[], double element_lhside[], 
  double element_matrix[], double element_rhside[] )

// 2d beam in x-y plane

{
  long int idim=0, jdim=0, inol=0, jnol=0, nnol=2, iuknwn=0, juknwn=0,
    ipuknwn=0, jpuknwn=0, idof=0, jdof=0, indx=0, swit=0, ldum=0, 
    displacement_indx=0, memory=-TOTAL_LINEAR, options_convection=-NO,
    length=0, icontrol=0, idum[1], options_mesh[MDIM];
  double E=0., I=0., initial_L=0., dtime=0., ddum[NDIM],
    beam_inertia[3]={0.,0.,0.},
    old_beam_direction[NDIM], new_beam_direction[NDIM],
    old_beam_dof[NNOL*NDOF], new_beam_dof[NNOL*NDOF], 
    incremental_beam_dof[NNOL*NDOF], incremental_beam_moment[NNOL*NDOF], 
    old_beam_moment[NNOL*NDOF], new_beam_moment[NNOL*NDOF],
    backrotated_beam_moment[NNOL*NDOF], forwardrotated_beam_moment[NNOL*NDOF],
    old_coord[NNOL*NDIM], new_coord[NNOL*NDIM], initial_coord[NNOL*NDIM],
    old_rotation_matrix[NNOL*NDOF*NNOL*NDOF], new_rotation_matrix[NNOL*NDOF*NNOL*NDOF], 
    local_beam_matrix[NNOL*NDOF*NNOL*NDOF], new_beam_matrix[NNOL*NDOF*NNOL*NDOF], 
    work[NNOL*NDOF*NNOL*NDOF];

  swit = set_swit(element,-1,"beam");
  if ( swit ) pri( "In routine BEAM." );

  db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET );
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  db( GROUP_BEAM_YOUNG, element_group, idum, &E, 
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  // group_beam_inertia (manual Professional 6.590): index Iyy Izz J.
  // The GNU beam is the 2D x-y beam bending about the LOCAL z axis, so
  // it consumes Izz = value 2. The legacy GNU single-value form (only
  // Izz given) is still read for backward compatibility.
  {
    long int inertia_length=0;
    if ( db( GROUP_BEAM_INERTIA, element_group, idum, beam_inertia,
        inertia_length, VERSION_NORMAL, GET_IF_EXISTS ) ) {
      if ( inertia_length==1 ) I = beam_inertia[0];
      else if ( inertia_length>=2 ) I = beam_inertia[1];
      else db_error( GROUP_BEAM_INERTIA, element_group );
    }
  }
  // manual Professional 6.591: group_beam_memory accepts -updated,
  // -updated_without_rotation and -total_linear; the default is
  // -total_linear (Professional 02-08-2026; F5 T2 V4d owner decision,
  // the GNU default was -updated). The GNU beam integrates a
  // geometrically linear law; -total_linear is mapped onto
  // -updated_without_rotation (initial coordinates fixed).
  db( GROUP_BEAM_MEMORY, element_group, &memory, ddum,
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  if ( memory==-TOTAL_LINEAR ) memory = -UPDATED_WITHOUT_ROTATION;
  if ( memory!=-UPDATED && memory!=-UPDATED_WITHOUT_ROTATION )
    db_error( GROUP_BEAM_MEMORY, element_group );

  if ( materi_velocity_integrated ) {
    options_convection = convection_apply_switch();
    db( OPTIONS_MESH, 0, options_mesh, ddum, ldum, VERSION_NORMAL, GET );
    for ( idim=0; idim<ndim; idim++ ) {
      if ( options_mesh[idim]==-FIXED_IN_SPACE ) {
        if ( options_convection==-YES ) {
          pri( "Error: set OPTIONS_CONVECTION to -NO in analysis with beams." );
          exit_tn_on_error();
        }
      }
    }
  }

  if ( materi_velocity_integrated )
    displacement_indx = veli_indx;
  else {
    assert( materi_displacement );
    displacement_indx = dis_indx;
  }

    // geometry
  for ( idim=0; idim<NDIM; idim++ ) {
    if      ( materi_velocity_integrated ) {
      for ( inol=0;inol<nnol; inol++ ) {
        initial_coord[inol*NDIM+idim] = 
          coord[inol*NDIM+idim];
        if ( options_mesh[idim]==-FOLLOW_MATERIAL )
          initial_coord[inol*NDIM+idim] -= 
          old_dof[inol*nuknwn+veli_indx+idim*nder];
        old_coord[inol*NDIM+idim] = 
          initial_coord[inol*NDIM+idim] + 
          coord[inol*NDIM+idim];
        new_coord[inol*NDIM+idim] = 
          old_coord[inol*NDIM+idim] +
          new_dof[inol*nuknwn+vel_indx+idim*nder]*dtime;
      }
    }
    else {
      assert( materi_displacement );
      for ( inol=0;inol<nnol; inol++ ) {
        initial_coord[inol*NDIM+idim] = 
          coord[inol*NDIM+idim];
        old_coord[inol*NDIM+idim] = 
          initial_coord[inol*NDIM+idim] + 
          old_dof[inol*nuknwn+dis_indx+idim*nder];
        new_coord[inol*NDIM+idim] = 
          old_coord[inol*NDIM+idim] + 
          new_dof[inol*nuknwn+vel_indx+idim*nder]*dtime;
      }
    }
  }

  if ( memory==-UPDATED_WITHOUT_ROTATION ) {
    array_move( initial_coord, old_coord, NNOL*NDIM );
    array_move( initial_coord, new_coord, NNOL*NDIM );
  }
  initial_L = array_distance( &initial_coord[0], 
    &initial_coord[NDIM], ddum, NDIM );
  if ( swit ) {
    pri( "coord", coord, NNOL, NDIM );
    pri( "initial_coord", initial_coord, NNOL, NDIM );
    pri( "old_coord", old_coord, NNOL, NDIM );
    pri( "new_coord", new_coord, NNOL, NDIM );
    pri( "initial_length", initial_L );
  }

    // beam dofs in global x-y axes
  for ( inol=0;inol<NNOL; inol++ ) {
    for ( idim=0; idim<NDIM; idim++ ) {
      idof = idim;
      indx = inol*NDOF+idof;
      iuknwn = displacement_indx + idim*nder;
      old_beam_dof[indx] = old_dof[inol*nuknwn+iuknwn];
      new_beam_dof[indx] = old_dof[inol*nuknwn+iuknwn] +
        new_dof[inol*nuknwn+vel_indx+idim*nder]*dtime;
      incremental_beam_dof[indx] = 
        new_beam_dof[indx] - old_beam_dof[indx];
    }
    idof = 2;
    indx = inol*NDOF+idof;
    if ( ndim==2 )
      iuknwn = rot_indx;
    else
      iuknwn = rot_indx + 2*nder;
    old_beam_dof[indx] = old_dof[inol*nuknwn+iuknwn];
    new_beam_dof[indx] = new_dof[inol*nuknwn+iuknwn];
    incremental_beam_dof[indx] = 
      new_beam_dof[indx] - old_beam_dof[indx];
  }
  if ( swit ) {
    pri( "old_beam_dof", old_beam_dof, NNOL, NDOF );
    pri( "new_beam_dof", new_beam_dof, NNOL, NDOF );
    pri( "incremental_beam_dof", incremental_beam_dof, NNOL, NDOF );
  }

    // beam direction
  array_subtract( &old_coord[NDIM], &old_coord[0], 
    old_beam_direction, NDIM );
  if ( !array_normalize( old_beam_direction, NDIM ) ) {
     pri( "Error: zero length detected in beam ", element );
     exit(TN_EXIT_STATUS);
  }
  array_subtract( &new_coord[NDIM], &new_coord[0], 
    new_beam_direction, NDIM );
  if ( !array_normalize( new_beam_direction, NDIM ) ) {
     pri( "Error: zero length detected in beam ", element );
     exit(TN_EXIT_STATUS);
  }
  if ( swit ) {
    pri( "old_beam_direction", old_beam_direction, NDIM );
    pri( "new_beam_direction", new_beam_direction, NDIM );
  }

    // rotation matrix
  array_set( old_rotation_matrix, 0., NNOL*NDOF*NNOL*NDOF );
  array_set( new_rotation_matrix, 0., NNOL*NDOF*NNOL*NDOF );
  old_rotation_matrix[0*NNOL*NDOF+0] = +old_beam_direction[0];
  old_rotation_matrix[0*NNOL*NDOF+1] = +old_beam_direction[1];
  old_rotation_matrix[1*NNOL*NDOF+0] = -old_beam_direction[1];
  old_rotation_matrix[1*NNOL*NDOF+1] = +old_beam_direction[0];
  old_rotation_matrix[2*NNOL*NDOF+2] = +1.;
  old_rotation_matrix[3*NNOL*NDOF+3] = +old_beam_direction[0];
  old_rotation_matrix[3*NNOL*NDOF+4] = +old_beam_direction[1];
  old_rotation_matrix[4*NNOL*NDOF+3] = -old_beam_direction[1];
  old_rotation_matrix[4*NNOL*NDOF+4] = +old_beam_direction[0];
  old_rotation_matrix[5*NNOL*NDOF+5] = +1.;
  new_rotation_matrix[0*NNOL*NDOF+0] = +new_beam_direction[0];
  new_rotation_matrix[0*NNOL*NDOF+1] = +new_beam_direction[1];
  new_rotation_matrix[1*NNOL*NDOF+0] = -new_beam_direction[1];
  new_rotation_matrix[1*NNOL*NDOF+1] = +new_beam_direction[0];
  new_rotation_matrix[2*NNOL*NDOF+2] = +1.;
  new_rotation_matrix[3*NNOL*NDOF+3] = +new_beam_direction[0];
  new_rotation_matrix[3*NNOL*NDOF+4] = +new_beam_direction[1];
  new_rotation_matrix[4*NNOL*NDOF+3] = -new_beam_direction[1];
  new_rotation_matrix[4*NNOL*NDOF+4] = +new_beam_direction[0];
  new_rotation_matrix[5*NNOL*NDOF+5] = +1.;
  if ( swit ) {
    pri( "old_rotation_matrix", old_rotation_matrix, NNOL*NDOF, NNOL*NDOF );
    pri( "new_rotation_matrix", new_rotation_matrix, NNOL*NDOF, NNOL*NDOF );
  }

    // beam matrix in local coordinates, dofs: UY1 ROTZ1 UY2 ROTZ2
  array_set( local_beam_matrix, 0., NNOL*NDOF*NNOL*NDOF );
    // UY1
  local_beam_matrix[1*NNOL*NDOF+1] = +(E*I/initial_L)*12./(initial_L*initial_L);
  local_beam_matrix[1*NNOL*NDOF+2] = +(E*I/initial_L)*6./initial_L;
  local_beam_matrix[1*NNOL*NDOF+4] = -(E*I/initial_L)*12./(initial_L*initial_L);
  local_beam_matrix[1*NNOL*NDOF+5] = +(E*I/initial_L)*6./initial_L;
    // ROTZ1
  local_beam_matrix[2*NNOL*NDOF+2] = +(E*I/initial_L)*4;
  local_beam_matrix[2*NNOL*NDOF+4] = -(E*I/initial_L)*6./initial_L;
  local_beam_matrix[2*NNOL*NDOF+5] = +(E*I/initial_L)*2.;
    // UY2
  local_beam_matrix[4*NNOL*NDOF+4] = +(E*I/initial_L)*12./(initial_L*initial_L);
  local_beam_matrix[4*NNOL*NDOF+5] = -(E*I/initial_L)*6./initial_L;
    // ROTZ2
  local_beam_matrix[5*NNOL*NDOF+5] = +(E*I/initial_L)*4.;
  matrix_symmetric( local_beam_matrix, NNOL*NDOF );
  if ( swit ) pri( "local_beam_matrix", local_beam_matrix, NNOL*NDOF, NNOL*NDOF );

    // global beam matrix
  matrix_atba( new_rotation_matrix, local_beam_matrix, 
    new_beam_matrix, work, NNOL*NDOF, NNOL*NDOF );
  if ( swit ) pri( "new_beam_matrix", new_beam_matrix, 
    NNOL*NDOF, NNOL*NDOF );

    // old beam forces in global axes
  array_set( old_beam_moment, 0., NNOL*NDOF );
  db( ELEMENT_BEAM_MOMENT, element, idum, old_beam_moment,
    length, VERSION_NORMAL, GET_IF_EXISTS );
  if ( swit ) pri( "old_beam_moment", old_beam_moment, NNOL, NDOF );

    // backrotated old beam forces in global axes
  matrix_atb( old_rotation_matrix, old_beam_moment, 
    backrotated_beam_moment, NNOL*NDOF, NNOL*NDOF, 1 );
  if ( swit ) pri( "backrotated_beam_moment", backrotated_beam_moment, NNOL, NDOF );

    // forwardrotated old beam forces in global axes
  matrix_ab( new_rotation_matrix, backrotated_beam_moment, 
    forwardrotated_beam_moment, NNOL*NDOF, NNOL*NDOF, 1 );
  if ( swit ) pri( "forwardrotated_beam_moment", forwardrotated_beam_moment, NNOL, NDOF );

    // incremental beam forces in global axes
  matrix_ab( new_beam_matrix, incremental_beam_dof, 
    incremental_beam_moment, NNOL*NDOF, NNOL*NDOF, 1 );
  if ( swit ) pri( "incremental_beam_moment", incremental_beam_moment, NNOL, NDOF );

    // new beam moment in global axes
  array_add( incremental_beam_moment, forwardrotated_beam_moment, new_beam_moment, NNOL*NDOF );
  if ( swit ) pri( "new_beam_moment", new_beam_moment, NNOL, NDOF );

    // fill in force vector and matrix
  for ( inol=0; inol<NNOL; inol++ ) {
    for ( idim=0; idim<NDIM; idim++ ) {
      idof = idim;
      iuknwn = vel_indx + idim*nder;
      ipuknwn = iuknwn/nder;
      element_rhside[inol*npuknwn+ipuknwn] -= new_beam_moment[inol*NDOF+idof];
      for ( jnol=0; jnol<NNOL; jnol++ ) {
        for ( jdim=0; jdim<NDIM; jdim++ ) {
          jdof = jdim;
          juknwn = vel_indx + jdim*nder;
          jpuknwn = juknwn/nder;
          indx = inol*npuknwn*NNOL*npuknwn+ipuknwn*NNOL*npuknwn+jnol*npuknwn+jpuknwn;
          element_matrix[indx] +=
            new_beam_matrix[(inol*NDOF+idof)*NNOL*NDOF+jnol*NDOF+jdof]*dtime;
        }
        jdof = 2;
        if ( ndim==2 )
          juknwn = rot_indx;
        else
          juknwn = rot_indx + 2*nder;
        jpuknwn = juknwn/nder;
        indx = inol*npuknwn*NNOL*npuknwn+ipuknwn*NNOL*npuknwn+jnol*npuknwn+jpuknwn;
        element_matrix[indx] +=
          new_beam_matrix[(inol*NDOF+idof)*NNOL*NDOF+jnol*NDOF+jdof];
      }
    }
    idof = 2;
    if ( ndim==2 )
      iuknwn = rot_indx;
    else
      iuknwn = rot_indx + 2*nder;
    ipuknwn = iuknwn/nder;
    element_rhside[inol*npuknwn+ipuknwn] -= new_beam_moment[inol*NDOF+idof];
    for ( jnol=0; jnol<NNOL; jnol++ ) {
      for ( jdim=0; jdim<NDIM; jdim++ ) {
        jdof = jdim;
        juknwn = vel_indx + jdim*nder;
        jpuknwn = juknwn/nder;
        indx = inol*npuknwn*NNOL*npuknwn+ipuknwn*NNOL*npuknwn+jnol*npuknwn+jpuknwn;
        element_matrix[indx] +=
          new_beam_matrix[(inol*NDOF+idof)*NNOL*NDOF+jnol*NDOF+jdof]*dtime;
      }
      jdof = 2;
      if ( ndim==2 )
        juknwn = rot_indx;
      else
        juknwn = rot_indx + 2*nder;
      jpuknwn = juknwn/nder;
      indx = inol*npuknwn*NNOL*npuknwn+ipuknwn*NNOL*npuknwn+jnol*npuknwn+jpuknwn;
      element_matrix[indx] +=
        new_beam_matrix[(inol*NDOF+idof)*NNOL*NDOF+jnol*NDOF+jdof];
    }
  }
  for ( indx=0; indx<NNOL*npuknwn; indx++ )
    element_lhside[indx] += element_matrix[indx*NNOL*npuknwn+indx];

  if ( swit ) {
    pri( "element_rhside", element_rhside, NNOL, npuknwn );
    pri( "element_lhside", element_lhside, NNOL, npuknwn );
    pri( "element_matrix", element_matrix, NNOL*npuknwn, NNOL*npuknwn );
  }
 
    // element_beam_direction (manual Professional 6.414, F5 T2 V4d):
    // local x, y and z axes; the 2D beam lies in the x-y plane
  {
    double axes[3*3] = { new_beam_direction[0], new_beam_direction[1], 0.,
      -new_beam_direction[1], new_beam_direction[0], 0., 0., 0., 1. };
    length = 3*3;
    db( ELEMENT_BEAM_DIRECTION, element, idum, axes,
      length, VERSION_NEW, PUT );
  }
  length = NNOL*NDOF;
  db( ELEMENT_BEAM_MOMENT, element, idum, new_beam_moment,
    length, VERSION_NEW, PUT );

    // element_beam_force_moment (manual Professional 6.416, F5 T2 V4c):
    // the forces and moments in the local beam axes; local x = beam
    // direction, local y = in-plane perpendicular (row 1 of the rotation
    // matrix), local z = out of the beam plane. The 2D beam only has
    // force_x, force_y and moment_z per node; the other components are 0.
  {
    double force_moment[2*6];
    array_set( force_moment, 0., 2*6 );
    for ( inol=0; inol<NNOL; inol++ ) {
      force_moment[inol*6+0] =
        new_rotation_matrix[0*NNOL*NDOF+0]*new_beam_moment[inol*NDOF+0] +
        new_rotation_matrix[0*NNOL*NDOF+1]*new_beam_moment[inol*NDOF+1];
      force_moment[inol*6+1] =
        new_rotation_matrix[1*NNOL*NDOF+0]*new_beam_moment[inol*NDOF+0] +
        new_rotation_matrix[1*NNOL*NDOF+1]*new_beam_moment[inol*NDOF+1];
      force_moment[inol*6+5] = new_beam_moment[inol*NDOF+2];
    }
    length = 2*6;
    db( ELEMENT_BEAM_FORCE_MOMENT, element, idum, force_moment,
      length, VERSION_NEW, PUT );
  }

  if ( swit ) pri( "Out function BEAM" );
}


long int beam_direction_z_record( long int element, long int element_group,
  double vec[] )

// Local beam z record of a beam element (F5 T2 V4c): returns the record
// that applies and fills vec[3] with its values, or returns 0 if none.
// Priority (Professional 02-08-2026 probes): element_beam_direction_z
// (6.415) over group_beam_direction_z (6.593) over
// group_beam_direction_z_reference_point (6.594).

{
  long int ldum=0, idum[1];

  array_set( vec, 0., 3 );
  if ( db( ELEMENT_BEAM_DIRECTION_Z, element, idum, vec, ldum,
      VERSION_NORMAL, GET_IF_EXISTS ) )
    return ELEMENT_BEAM_DIRECTION_Z;
  if ( db( GROUP_BEAM_DIRECTION_Z, element_group, idum, vec, ldum,
      VERSION_NORMAL, GET_IF_EXISTS ) )
    return GROUP_BEAM_DIRECTION_Z;
  if ( db( GROUP_BEAM_DIRECTION_Z_REFERENCE_POINT, element_group, idum,
      vec, ldum, VERSION_NORMAL, GET_IF_EXISTS ) )
    return GROUP_BEAM_DIRECTION_Z_REFERENCE_POINT;
  return 0;
}

#define BEAM_LOCAL_NLOC 12  // local beam dofs: 2 nodes x (ux uy uz rx ry rz)
#define BEAM_LOCAL_PERPENDICULAR_TOL 1.e-3 // |cos| of direction_z vs beam axis
#define BEAM_LOCAL_VERTICAL_TOL 1.e-6 // beam along z: default local z = global y
#define BEAM_LOCAL_ZERO_TOL 1.e-10 // relative length of a vanishing z axis

void beam_local( long int element, long int element_group, long int record,
  double record_vec[], double coord[], double old_dof[], double new_dof[],
  double element_lhside[], double element_matrix[], double element_rhside[] )

// Beam with local beam axes (F5 T2 V4c; manual Professional 6.415, 6.416,
// 6.593, 6.594), used when a local z record applies to the beam (see
// beam_direction_z_record) and for every 3D beam (record 0: default local
// z, F5 T2 V4d). Euler-Bernoulli beam in the local axes:
// bending in the local x-y plane with E*Izz, bending in the local x-z
// plane with E*Iyy, torsion with G*J (group_beam_young E, group_beam_shear
// G, group_beam_inertia Iyy Izz J). As the classic GNU beam it has no
// axial stiffness (the Professional neither: group_beam_area is not a
// Professional record; combine with a truss for the axial force).
// Local axes: x = beam direction node 1 -> node 2, z = the part of the
// record vector perpendicular to x (direction record: the vector itself,
// which must be perpendicular to the beam; reference point: the vector
// from the beam middle point to the point), y = z cross x.
// The forces and moments are kept in ELEMENT_BEAM_FORCE_MOMENT in the
// local axes (the record is both the output of 6.416 and the state); a
// new step adds K_local * (T * incremental global dofs), which keeps the
// local forces attached to the rotating beam (as beam_2d with its
// backrotated / forwardrotated forces). In 2D only velx, vely and the
// single rotation (about global z) are assembled.

{
  long int i=0, j=0, k=0, inol=0, jnol=0, idim=0, s=0, t=0, swit=0, ldum=0,
    length=0, icontrol=0, memory=-TOTAL_LINEAR, options_convection=-NO,
    displacement_indx=0, inertia_length=0, idum[1], options_mesh[MDIM],
    present[BEAM_LOCAL_NLOC], puknwn[BEAM_LOCAL_NLOC],
    is_disp[BEAM_LOCAL_NLOC];
  double E=0., G=0., Iyy=0., Izz=0., J=0., dtime=0., L=0., Lnew=0., c=0.,
    vlen=0., zlen=0., EI=0., ddum[MDIM], beam_inertia[3]={0.,0.,0.},
    initial_coord[2*3], old_coord[2*3], new_coord[2*3], ex[3], ey[3],
    ez[3], ex0[3], v[3], R[3*3], T[BEAM_LOCAL_NLOC*BEAM_LOCAL_NLOC],
    K[BEAM_LOCAL_NLOC*BEAM_LOCAL_NLOC], KT[BEAM_LOCAL_NLOC*BEAM_LOCAL_NLOC],
    Kg[BEAM_LOCAL_NLOC*BEAM_LOCAL_NLOC], du_g[BEAM_LOCAL_NLOC],
    du_l[BEAM_LOCAL_NLOC], f_l[BEAM_LOCAL_NLOC], f_g[BEAM_LOCAL_NLOC];
  const long int nnol=2, n=BEAM_LOCAL_NLOC;

  swit = set_swit(element,-1,"beam_local");
  if ( swit ) pri( "In routine BEAM_LOCAL." );

  db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET );
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  db( GROUP_BEAM_YOUNG, element_group, idum, &E,
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( GROUP_BEAM_SHEAR, element_group, idum, &G,
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  // group_beam_inertia Iyy Izz J; the legacy GNU single-value form is
  // the bending inertia for both axes (no torsion stiffness)
  if ( db( GROUP_BEAM_INERTIA, element_group, idum, beam_inertia,
      inertia_length, VERSION_NORMAL, GET_IF_EXISTS ) ) {
    if      ( inertia_length==1 ) {
      Iyy = beam_inertia[0]; Izz = beam_inertia[0];
    }
    else if ( inertia_length==2 ) {
      Iyy = beam_inertia[0]; Izz = beam_inertia[1];
    }
    else if ( inertia_length==3 ) {
      Iyy = beam_inertia[0]; Izz = beam_inertia[1]; J = beam_inertia[2];
    }
    else db_error( GROUP_BEAM_INERTIA, element_group );
  }
  db( GROUP_BEAM_MEMORY, element_group, &memory, ddum,
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  if ( memory==-TOTAL_LINEAR ) memory = -UPDATED_WITHOUT_ROTATION;
  if ( memory!=-UPDATED && memory!=-UPDATED_WITHOUT_ROTATION )
    db_error( GROUP_BEAM_MEMORY, element_group );

  if ( materi_velocity_integrated ) {
    options_convection = convection_apply_switch();
    db( OPTIONS_MESH, 0, options_mesh, ddum, ldum, VERSION_NORMAL, GET );
    for ( idim=0; idim<ndim; idim++ ) {
      if ( options_mesh[idim]==-FIXED_IN_SPACE ) {
        if ( options_convection==-YES ) {
          pri( "Error: set OPTIONS_CONVECTION to -NO in analysis with beams." );
          exit_tn_on_error();
        }
      }
    }
    displacement_indx = veli_indx;
  }
  else {
    assert( materi_displacement );
    displacement_indx = dis_indx;
  }

    // geometry (3 components; z = 0 in 2D): initial, old and new
    // coordinates of the two nodes
  array_set( initial_coord, 0., 2*3 );
  for ( inol=0; inol<nnol; inol++ ) {
    for ( idim=0; idim<ndim; idim++ ) {
      initial_coord[inol*3+idim] = coord[inol*ndim+idim];
      if ( materi_velocity_integrated &&
           options_mesh[idim]==-FOLLOW_MATERIAL )
        initial_coord[inol*3+idim] -=
          old_dof[inol*nuknwn+veli_indx+idim*nder];
    }
  }
  array_move( initial_coord, old_coord, 2*3 );
  array_move( initial_coord, new_coord, 2*3 );
  if ( memory==-UPDATED ) {
    for ( inol=0; inol<nnol; inol++ ) {
      for ( idim=0; idim<ndim; idim++ ) {
        old_coord[inol*3+idim] +=
          old_dof[inol*nuknwn+displacement_indx+idim*nder];
        new_coord[inol*3+idim] = old_coord[inol*3+idim] +
          new_dof[inol*nuknwn+vel_indx+idim*nder]*dtime;
      }
    }
  }
  array_subtract( &initial_coord[3], &initial_coord[0], ex0, 3 );
  L = array_size( ex0, 3 );
  array_subtract( &new_coord[3], &new_coord[0], ex, 3 );
  Lnew = array_size( ex, 3 );
  if ( L<=0. || Lnew<=0. ) {
    pri( "Error: zero length detected in beam ", element );
    exit_tn_on_error();
  }
  array_multiply( ex0, ex0, 1./L, 3 );
  array_multiply( ex, ex, 1./Lnew, 3 );

    // local z axis
  if ( !record ) {
      // no local z record (3D, F5 T2 V4d): the normal of the GNU-only
      // group_beam_plane (first axis x second axis) if given, else the
      // Professional default, global z (global y for a beam along z)
    long int group_beam_plane[2], index_plane[2];
    array_set( v, 0., 3 );
    if ( db( GROUP_BEAM_PLANE, element_group, group_beam_plane, ddum,
        ldum, VERSION_NORMAL, GET_IF_EXISTS ) ) {
      double a[3], b[3];
      for ( idim=0; idim<2; idim++ ) {
        if      ( group_beam_plane[idim]==-X ) index_plane[idim] = 0;
        else if ( group_beam_plane[idim]==-Y ) index_plane[idim] = 1;
        else if ( group_beam_plane[idim]==-Z ) index_plane[idim] = 2;
        else db_error( GROUP_BEAM_PLANE, element_group );
      }
      if ( index_plane[0]>=index_plane[1] )
        db_error( GROUP_BEAM_PLANE, element_group );
      array_set( a, 0., 3 ); a[index_plane[0]] = 1.;
      array_set( b, 0., 3 ); b[index_plane[1]] = 1.;
      array_outproduct_3D( a, b, v );
      if ( scalar_dabs( array_inproduct( v, ex0, 3 ) ) >
           1.-BEAM_LOCAL_PERPENDICULAR_TOL ) {
        pri( "Error, the beam is perpendicular to its GROUP_BEAM_PLANE. "
          "Element", element );
        exit_tn_on_error();
      }
    }
    else if ( sqrt( ex0[0]*ex0[0] + ex0[1]*ex0[1] ) <=
              BEAM_LOCAL_VERTICAL_TOL )
      v[1] = 1.;
    else
      v[2] = 1.;
  }
  else if ( record==GROUP_BEAM_DIRECTION_Z_REFERENCE_POINT ) {
    for ( idim=0; idim<3; idim++ )
      v[idim] = record_vec[idim] -
        0.5*( new_coord[idim] + new_coord[3+idim] );
  }
  else {
    array_move( record_vec, v, 3 );
    vlen = array_size( v, 3 );
    if ( vlen<=0. ||
         scalar_dabs( array_inproduct( v, ex0, 3 ) ) >
         BEAM_LOCAL_PERPENDICULAR_TOL*vlen ) {
      if ( record==ELEMENT_BEAM_DIRECTION_Z )
        pri( "Error, ELEMENT_BEAM_DIRECTION_Z should be perpendicular "
          "to the beam length direction. Element", element );
      else
        pri( "Error, GROUP_BEAM_DIRECTION_Z should be perpendicular "
          "to the beam length direction. Element", element );
      exit_tn_on_error();
    }
  }
  vlen = array_size( v, 3 );
  c = array_inproduct( v, ex, 3 );
  for ( idim=0; idim<3; idim++ ) ez[idim] = v[idim] - c*ex[idim];
  zlen = array_size( ez, 3 );
  if ( zlen<=BEAM_LOCAL_ZERO_TOL*( vlen+L ) ) {
    pri( "Error, GROUP_BEAM_DIRECTION_Z_REFERENCE_POINT lies on the "
      "beam length axis. Element", element );
    exit_tn_on_error();
  }
  array_multiply( ez, ez, 1./zlen, 3 );
  array_outproduct_3D( ez, ex, ey );
  for ( idim=0; idim<3; idim++ ) {
    R[0*3+idim] = ex[idim];
    R[1*3+idim] = ey[idim];
    R[2*3+idim] = ez[idim];
  }
  if ( swit ) {
    pri( "initial_coord", initial_coord, 2, 3 );
    pri( "new_coord", new_coord, 2, 3 );
    pri( "local axes", R, 3, 3 );
  }

    // transformation local = T global (4 blocks of R)
  array_set( T, 0., n*n );
  for ( k=0; k<4; k++ ) {
    for ( i=0; i<3; i++ ) {
      for ( j=0; j<3; j++ )
        T[(k*3+i)*n+k*3+j] = R[i*3+j];
    }
  }

    // local stiffness; local dofs per node: ux uy uz rx ry rz
  array_set( K, 0., n*n );
  {
    long int dz[4] = { 1, 5, 7, 11 }; // uy1 rz1 uy2 rz2 (x-y plane, Izz)
    long int dy[4] = { 2, 4, 8, 10 }; // uz1 ry1 uz2 ry2 (x-z plane, Iyy)
    double kz[16] = { 12., 6.*L, -12., 6.*L,
                      6.*L, 4.*L*L, -6.*L, 2.*L*L,
                      -12., -6.*L, 12., -6.*L,
                      6.*L, 2.*L*L, -6.*L, 4.*L*L };
    double ky[16] = { 12., -6.*L, -12., -6.*L,
                      -6.*L, 4.*L*L, 6.*L, 2.*L*L,
                      -12., 6.*L, 12., 6.*L,
                      -6.*L, 2.*L*L, 6.*L, 4.*L*L };
    EI = E*Izz/(L*L*L);
    for ( i=0; i<4; i++ ) {
      for ( j=0; j<4; j++ ) K[dz[i]*n+dz[j]] += EI*kz[i*4+j];
    }
    EI = E*Iyy/(L*L*L);
    for ( i=0; i<4; i++ ) {
      for ( j=0; j<4; j++ ) K[dy[i]*n+dy[j]] += EI*ky[i*4+j];
    }
    K[3*n+3] += G*J/L; K[3*n+9] -= G*J/L;
    K[9*n+3] -= G*J/L; K[9*n+9] += G*J/L;
  }

    // global beam matrix Kg = T^t K T
  for ( i=0; i<n; i++ ) {
    for ( j=0; j<n; j++ ) {
      KT[i*n+j] = 0.;
      for ( k=0; k<n; k++ ) KT[i*n+j] += K[i*n+k]*T[k*n+j];
    }
  }
  for ( i=0; i<n; i++ ) {
    for ( j=0; j<n; j++ ) {
      Kg[i*n+j] = 0.;
      for ( k=0; k<n; k++ ) Kg[i*n+j] += T[k*n+i]*KT[k*n+j];
    }
  }

    // global dofs of the 12 slots (velocity for the translations, the
    // rotation dof itself for the rotations); 2D: velx, vely, rotz
  for ( inol=0; inol<nnol; inol++ ) {
    for ( s=0; s<6; s++ ) {
      i = inol*6 + s;
      present[i] = 0; puknwn[i] = 0; is_disp[i] = 0; du_g[i] = 0.;
      if ( s<3 ) {
        if ( s<ndim ) {
          present[i] = 1; is_disp[i] = 1;
          puknwn[i] = ( vel_indx + s*nder )/nder;
          du_g[i] = new_dof[inol*nuknwn+vel_indx+s*nder]*dtime;
        }
      }
      else if ( ndim==3 || s==5 ) {
        t = ( ndim==3 ) ? rot_indx + (s-3)*nder : rot_indx;
        present[i] = 1;
        puknwn[i] = t/nder;
        du_g[i] = new_dof[inol*nuknwn+t] - old_dof[inol*nuknwn+t];
      }
    }
  }

    // new local forces = old local forces + K * T * incremental dofs
  matrix_ab( T, du_g, du_l, n, n, 1 );
  array_set( f_l, 0., n );
  db( ELEMENT_BEAM_FORCE_MOMENT, element, idum, f_l, length,
    VERSION_NORMAL, GET_IF_EXISTS );
  for ( i=0; i<n; i++ ) {
    for ( j=0; j<n; j++ ) f_l[i] += K[i*n+j]*du_l[j];
  }
  matrix_atb( T, f_l, f_g, n, n, 1 );
  if ( swit ) {
    pri( "du_g", du_g, nnol, 6 );
    pri( "f_l", f_l, nnol, 6 );
  }

    // fill in force vector and matrix
  for ( i=0; i<n; i++ ) {
    if ( !present[i] ) continue;
    inol = i/6;
    element_rhside[inol*npuknwn+puknwn[i]] -= f_g[i];
    for ( j=0; j<n; j++ ) {
      if ( !present[j] ) continue;
      jnol = j/6;
      element_matrix[(inol*npuknwn+puknwn[i])*nnol*npuknwn+
        jnol*npuknwn+puknwn[j]] += Kg[i*n+j]*( is_disp[j] ? dtime : 1. );
    }
  }
  for ( i=0; i<nnol*npuknwn; i++ )
    element_lhside[i] += element_matrix[i*nnol*npuknwn+i];

  length = n;
  db( ELEMENT_BEAM_FORCE_MOMENT, element, idum, f_l, length,
    VERSION_NEW, PUT );
  length = 3*3;
  db( ELEMENT_BEAM_DIRECTION, element, idum, R, length, VERSION_NEW, PUT );

  if ( swit ) pri( "Out routine BEAM_LOCAL." );
}
