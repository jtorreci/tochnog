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

// 3d beam 
// mechanical action of beam in a 2d plane (x-y plane, or x-z plane, or y-z plane)

{
  long int n=0, indx=0, indx_2d=0, indx_3d=0, jndx_2d=0, jndx_3d=0, 
    idim=0, jdim=0, kdim=0, ldim=0, inol=0, jnol=0, iuknwn_2d=0, 
    iuknwn_3d=0, juknwn_2d=0, juknwn_3d=0, ldum=0, displacement_indx=0, 
    swit=0, group_beam_plane[NDIM], index_plane[MDIM];
  double ddum[1], *dbl_array=NULL, *coord_2d=NULL, *old_dof_2d=NULL, 
    *new_dof_2d=NULL, *element_lhside_2d=NULL, *element_matrix_2d=NULL, 
    *element_rhside_2d=NULL;

  if ( ndim==2 ) {
    beam_2d( element, element_group, coord_3d, old_dof_3d, new_dof_3d, 
      element_lhside_3d, element_matrix_3d, element_rhside_3d );
    return;
  }

  swit = set_swit(element,-1,"beam_3d");
  if ( swit ) pri( "In routine BEAM_3D." );

    // allocate

  n = NNOL*NDIM + NNOL*nuknwn + NNOL*nuknwn +
    NNOL*npuknwn + NNOL*npuknwn*NNOL*npuknwn + NNOL*npuknwn;
  dbl_array = get_new_dbl(n);

  indx = 0;
  coord_2d = &dbl_array[indx]; indx += NNOL*NDIM;
  old_dof_2d = &dbl_array[indx]; indx += NNOL*nuknwn;
  new_dof_2d = &dbl_array[indx]; indx += NNOL*nuknwn;
  element_lhside_2d = &dbl_array[indx]; indx += NNOL*npuknwn;
  element_matrix_2d = &dbl_array[indx]; indx += NNOL*npuknwn*NNOL*npuknwn;
  element_rhside_2d = &dbl_array[indx]; indx += NNOL*npuknwn;

    // initialise

  array_set( element_lhside_2d, 0., NNOL*nuknwn );
  array_set( element_rhside_2d, 0., NNOL*nuknwn );
  array_set( element_matrix_2d, 0., NNOL*nuknwn*NNOL*nuknwn );

  if ( materi_velocity_integrated )
    displacement_indx = veli_indx;
  else {
    assert( materi_displacement );
    displacement_indx = dis_indx;
  }

    // beam plane

  group_beam_plane[0] = -X;
  group_beam_plane[1] = -Y;
  if ( db( GROUP_BEAM_PLANE, element_group, group_beam_plane, ddum, 
      ldum, VERSION_NORMAL, GET_IF_EXISTS ) ) {
    if ( labs(group_beam_plane[0])>=
         labs(group_beam_plane[1]) ) {
      db_error( GROUP_BEAM_PLANE, element_group );
    }
  }
  for ( idim=0; idim<NDIM; idim++ ) {
    if      ( group_beam_plane[idim] == -X ) 
      index_plane[idim] = 0;
    else if ( group_beam_plane[idim] == -Y ) 
      index_plane[idim] = 1;
    else if ( group_beam_plane[idim] == -Z )
      index_plane[idim] = 2;
    else
      db_error( GROUP_BEAM_PLANE, element_group );
  }
  if      ( index_plane[0]==0 && index_plane[1]==1 ) 
    index_plane[2] = 2;
  else if ( index_plane[0]==0 && index_plane[1]==2 )
    index_plane[2] = 1;
  else if ( index_plane[0]==1 && index_plane[1]==2 )
    index_plane[2] = 0;
  else
    db_error( GROUP_BEAM_PLANE, element_group );
  if ( swit ) pri( "index_plane", index_plane, MDIM );

    // map to 2d

  for ( inol=0; inol<NNOL; inol++ ) {
    for ( idim=0; idim<NDIM+1; idim++ ) {
      kdim = index_plane[idim];
      if ( idim<NDIM ) {
        coord_2d[inol*NDIM+idim] = coord_3d[inol*ndim+kdim];
        iuknwn_2d = displacement_indx + idim*nder;
        iuknwn_3d = displacement_indx + kdim*nder;
        indx_2d = inol*nuknwn + iuknwn_2d;
        indx_3d = inol*nuknwn + iuknwn_3d;
        old_dof_2d[indx_2d] = old_dof_3d[indx_3d];
        new_dof_2d[indx_2d] = new_dof_3d[indx_3d];
        iuknwn_2d = vel_indx + idim*nder;
        iuknwn_3d = vel_indx + kdim*nder;
        indx_2d = inol*nuknwn + iuknwn_2d;
        indx_3d = inol*nuknwn + iuknwn_3d;
      }
      else {
        iuknwn_2d = rot_indx + 2*nder;
        iuknwn_3d = rot_indx + kdim*nder;
        indx_2d = inol*nuknwn + iuknwn_2d;
        indx_3d = inol*nuknwn + iuknwn_3d;
      }
      old_dof_2d[indx_2d] = old_dof_3d[indx_3d];
      new_dof_2d[indx_2d] = new_dof_3d[indx_3d];
      element_lhside_2d[indx_2d] = element_lhside_3d[indx_3d];
      element_rhside_2d[indx_2d] = element_rhside_3d[indx_3d];
      for ( jnol=0; jnol<NNOL; jnol++ ) {
        for ( jdim=0; jdim<NDIM+1; jdim++ ) {
          ldim = index_plane[jdim];
          if ( jdim<NDIM ) {
            juknwn_2d = vel_indx + jdim*nder;
            juknwn_3d = vel_indx + ldim*nder;
            jndx_2d = jnol*nuknwn + juknwn_2d;
            jndx_3d = jnol*nuknwn + juknwn_3d;
          }
          else {
            juknwn_2d = rot_indx + 2*nder;
            juknwn_3d = rot_indx + ldim*nder;
            jndx_2d = jnol*nuknwn + juknwn_2d;
            jndx_3d = jnol*nuknwn + juknwn_3d;
          }
          element_matrix_2d[indx_2d*NNOL*nuknwn+jndx_2d] = 
            element_matrix_3d[indx_3d*NNOL*nuknwn+jndx_3d];
        }
      }
    }
  }
  if ( swit ) {
    pri( "coord_3d", coord_3d, NNOL, ndim );
    pri( "coord_2d", coord_2d, NNOL, NDIM );
    pri( "old_dof_3d", old_dof_3d, NNOL, nuknwn );
    pri( "old_dof_2d", old_dof_2d, NNOL, nuknwn );
    pri( "new_dof_3d", new_dof_3d, NNOL, nuknwn );
    pri( "new_dof_2d", new_dof_2d, NNOL, nuknwn );
  }

    // call 2d beam

  beam_2d( element, element_group, coord_2d, old_dof_2d, new_dof_2d, 
    element_lhside_2d, element_matrix_2d, element_rhside_2d );
  if ( swit ) {
    pri( "element_lhside_2d", element_lhside_2d, NNOL, nuknwn );
    pri( "element_rhside_2d", element_rhside_2d, NNOL, nuknwn );
    pri( "element_matrix_2d", element_matrix_2d, NNOL*nuknwn, NNOL*nuknwn );
    pri( "element_lhside_3d", element_lhside_3d, NNOL, nuknwn );
    pri( "element_rhside_3d", element_rhside_3d, NNOL, nuknwn );
    pri( "element_matrix_3d", element_matrix_3d, NNOL*nuknwn, NNOL*nuknwn );
  }

    // map to 3d

  for ( inol=0; inol<NNOL; inol++ ) {
    for ( idim=0; idim<NDIM+1; idim++ ) {
      kdim = index_plane[idim];
      if ( idim<NDIM ) {
        iuknwn_2d = vel_indx + idim*nder;
        iuknwn_3d = vel_indx + kdim*nder;
        indx_2d = inol*nuknwn + iuknwn_2d;
        indx_3d = inol*nuknwn + iuknwn_3d;
      }
      else {
        iuknwn_2d = rot_indx + 2*nder;
        iuknwn_3d = rot_indx + kdim*nder;
        indx_2d = inol*nuknwn + iuknwn_2d;
        indx_3d = inol*nuknwn + iuknwn_3d;
      }
      element_lhside_3d[indx_3d] = element_lhside_2d[indx_2d];
      element_rhside_3d[indx_3d] = element_rhside_2d[indx_2d];
      for ( jnol=0; jnol<NNOL; jnol++ ) {
        for ( jdim=0; jdim<NDIM+1; jdim++ ) {
          ldim = index_plane[jdim];
          if ( jdim<NDIM ) {
            juknwn_2d = vel_indx + jdim*nder;
            juknwn_3d = vel_indx + ldim*nder;
            jndx_2d = jnol*nuknwn + juknwn_2d;
            jndx_3d = jnol*nuknwn + juknwn_3d;
          }
          else {
            juknwn_2d = rot_indx + 2*nder;
            juknwn_3d = rot_indx + ldim*nder;
            jndx_2d = jnol*nuknwn + juknwn_2d;
            jndx_3d = jnol*nuknwn + juknwn_3d;
          }
          element_matrix_3d[indx_3d*NNOL*nuknwn+jndx_3d] = 
            element_matrix_2d[indx_2d*NNOL*nuknwn+jndx_2d];
        }
      }
    }
  }

  delete[] dbl_array;
  if ( swit ) pri( "Out routine BEAM_3D." );
}


void beam_2d( long int element, long int element_group, double coord[], 
  double old_dof[], double new_dof[], double element_lhside[], 
  double element_matrix[], double element_rhside[] )

// 2d beam in x-y plane

{
  long int idim=0, jdim=0, inol=0, jnol=0, nnol=2, iuknwn=0, juknwn=0,
    ipuknwn=0, jpuknwn=0, idof=0, jdof=0, indx=0, swit=0, ldum=0, 
    displacement_indx=0, memory=-UPDATED, options_convection=-YES,
    length=0, icontrol=0, idum[1], options_mesh[MDIM];
  double E=0., I=0., initial_L=0., dtime=0., ddum[NDIM], 
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
  db( GROUP_BEAM_INERTIA, element_group, idum, &I, 
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  if ( db( GROUP_BEAM_MEMORY, element_group, &memory, ddum,
      ldum, VERSION_NORMAL, GET_IF_EXISTS ) ) {
    if ( memory!=-UPDATED && memory!=-UPDATED_WITHOUT_ROTATION )
      db_error( GROUP_BEAM_MEMORY, element_group );
  }

  if ( materi_velocity_integrated ) {
    db( OPTIONS_CONVECTION, 0, &options_convection, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( CONTROL_OPTIONS_CONVECTION, icontrol, &options_convection, 
      ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
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
 
  length = NDIM;
  db( ELEMENT_BEAM_DIRECTION, element, idum, new_beam_direction,
    length, VERSION_NEW, PUT );
  length = NNOL*NDOF;
  db( ELEMENT_BEAM_MOMENT, element, idum, new_beam_moment,
    length, VERSION_NEW, PUT );

  if ( swit ) pri( "Out function BEAM" );
}
