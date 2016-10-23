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

#define EPS_MATRIX 1.e-10

long int *node_nel;
double   *node_lhside;
double   *node_rhside;
double   *node_dof_tmp;

void element_loop( void )
{
  long int swit=0, nthread=0, length=0, inod=0, max_node=0, 
    ithread=0, indx=0, idim=0, ipuknwn=0, iuknwn=0, ldum=0, 
    *node_nel_ptr=NULL;
  double ddum[1], *node_lhside_ptr=NULL, *node_rhside_ptr=NULL, 
    *node_dof_tmp_ptr=NULL, *node_dof_new_ptr=NULL;

  swit = set_swit(-1,-1,"element_loop");
  if ( swit ) pri( "In routine ELEMENT_LOOP" );

  db_max_index( NODE, max_node, VERSION_NORMAL, GET );
  db( OPTIONS_PROCESSORS, 0, &nthread, ddum, ldum, VERSION_NORMAL, GET );

  length = (1+max_node)*nthread;
  node_nel = get_new_int( length );
  array_set( node_nel, 0, length );

  length = (1+max_node)*npuknwn*nthread;
  node_lhside = get_new_dbl( length );
  array_set( node_lhside, 0., length );

  length = (1+max_node)*npuknwn*nthread;
  node_rhside = get_new_dbl( length );
  array_set( node_rhside, 0., length );

  length = (1+max_node)*nuknwn*nthread;
  node_dof_tmp = get_new_dbl( length );
  array_set( node_dof_tmp, 0., length );

  parallel_sys_routine( &parallel_element_loop );

  db_set_dbl( NODE_DOF_TMP, VERSION_NORMAL );
  db_set_int( NODE_NEL, VERSION_NORMAL );

    // merge results for different threads
  for ( inod=0; inod<=max_node; inod++ ) {
    if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
      node_nel_ptr = db_int( NODE_NEL, inod, VERSION_NORMAL );
      node_lhside_ptr = db_dbl( NODE_LHSIDE, inod, VERSION_NORMAL );
      node_rhside_ptr = db_dbl( NODE_RHSIDE, inod, VERSION_NORMAL );
      node_dof_new_ptr = db_dbl( NODE_DOF, inod, VERSION_NEW );
      node_dof_tmp_ptr = db_dbl( NODE_DOF_TMP, inod, VERSION_NORMAL );
      for ( ithread=0; ithread<nthread; ithread++ ) {
        indx = ithread*(1+max_node)+inod;
        node_nel_ptr[0] += node_nel[indx];
        indx = ithread*(1+max_node)*npuknwn+inod*npuknwn;
        array_add( node_lhside_ptr, &node_lhside[indx], 
          node_lhside_ptr, npuknwn );
        array_add( node_rhside_ptr, &node_rhside[indx], 
          node_rhside_ptr, npuknwn );
        indx = ithread*(1+max_node)*nuknwn+inod*nuknwn;
        array_add( node_dof_tmp_ptr, &node_dof_tmp[indx], 
          node_dof_tmp_ptr, nuknwn );
      }
      if ( derivatives && node_nel_ptr[0]>0 ) {
        for ( idim=0; idim<ndim; idim++ ) {
          for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
            iuknwn = ipuknwn*nder + idim + 1;
            node_dof_new_ptr[iuknwn] = 
              node_dof_tmp_ptr[iuknwn]/((double)node_nel_ptr[0]);
          }
        }
      }
    }
  }

  db_delete( NODE_DOF_TMP, VERSION_NORMAL );
  db_delete( NODE_NEL, VERSION_NORMAL );

  delete[] node_nel;
  delete[] node_lhside;
  delete[] node_rhside;
  delete[] node_dof_tmp;

  if ( swit ) pri( "Out routine ELEMENT_LOOP" );
}

void parallel_element_loop( void )

{
  long int element=0, max_element=0, iloop=0, nloop=0, swit=0,
    ithread=0, *next_of_loop=NULL;

  swit = set_swit(-1,-1,"parallel_element_loop");
  if ( swit ) pri( "In routine PARALLEL_ELEMENT_LOOP" );

    // loop over elements
  db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
  if ( max_element>=0 ) {
    next_of_loop = get_new_int(1+max_element);
    parallel_sys_next_of_loop( next_of_loop, max_element, nloop, ithread );
    for ( iloop=0; iloop<nloop; iloop++ ) {
      element = next_of_loop[iloop];
      if ( element>max_element )
        break;
      else if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) )
        elem( element, ithread );
    }
    delete[] next_of_loop;
  }

  if ( swit ) pri( "Out routine PARALLEL_ELEMENT_LOOP" );
}

void elem( long int element, long int ithread )

  /* 
     Left hand-side and right hand-side in an element.

     old_unknowns -> unknown values at time t
     new_unknowns -> unknown values at time t+dt
     old_grad -> gradient over shape at time 0 (total) or at time t (updated)
     new_grad -> gradient over shape at time t+dt 
     coord -> nodal coordinates at time 0 (total) or at time t (updated)
     new_coord -> nodal coordinates at time t+dt
  */

{
  long int inol=0, inod=0, ipuknwn=0, iuknwn=0, jpuknwn=0, juknwn=0, 
    iprinc=0, jnol=0, jnod=0, len_type=0, indxi=0, indxj=0,
    idim=0, type=0, swit=0, itype=0, ipoint=0, length=0, npoint=0,
    nnol=0, name=0, icontrol=0, element_group=0,  element_empty=0,
    indx=0, max_node=0, memory=-UPDATED, options_solver=-MATRIX_ITERATIVE_BICG,
    options_matrix_group=-NO, options_matrix_length=0, mnolnuknwn=npointmax*nuknwn,
    length_element_matrix_values=0, filled=0, membrane=-NO, test1=0, test2=0,
    ldum=0, n=0, axisymmetric=0, plasti_on_boundary=0, idum[1], *types=NULL, 
    *el=NULL, *nodes=NULL, *dof_type=NULL, element_dof_initialised=0,
    *dof_principal=NULL, *element_matrix_unknowns=NULL, 
    *group_matrix_unknowns=NULL, *int_array=NULL, one=1,
    length_nei=1+npointmax*ndim+npointmax+2;
  double dtime=0., volfac=0., volumeip=0., tmp=0., tmp1=0., tmp2=0.,
    element_delete_factor=0., dens=0., vel=0., materi_dens=0.,
    materi_density_minimum=0., materi_diffusion_minimum=0.,
    element_mass=0., res=0., radius=0., birth=0., death=0.,
    control_relaxation_materi_velocity=0., control_relaxation_condif_temperature=0.,
    control_relaxation_groundflow_pressure=0., control_relaxation_wave_fscalar=0.,
    control_relaxation_maxwell_e=0., element_strainenergy=0.,
    element_volume=0., time_current=0., time_total=0., 
    ddum[1], group_time[2], *coord=NULL, *new_coord=NULL, 
    *h=NULL, *coord_ip=NULL, *old_dof=NULL, *new_dof=NULL, 
    *old_unknowns=NULL, *new_unknowns=NULL, 
    *old_grad_old_unknowns=NULL, *old_grad_new_unknowns=NULL, 
    *new_grad_old_unknowns=NULL, *new_grad_new_unknowns=NULL, 
    *volume=NULL, *force_element_volume=NULL, 
    *tendon_element_rhside=NULL, 
    *residue_factor=NULL, *element_residue=NULL, 
    *element_dof_tmp_factor=NULL, *element_dof_tmp=NULL, 
    *element_softvar_tmp=NULL, 
    *element_lhside=NULL, *element_rhside=NULL, 
    *element_matrix_delete=NULL, *element_rhside_delete=NULL, 
    *condif_flow=NULL, *element_matrix_values=NULL, *element_matrix=NULL, 
    *element_matrix_second_values=NULL, 
    *element_matrix_second=NULL, *new_b=NULL, 
    *old_d=NULL, *new_d=NULL, 
    *massflow=NULL, *grad_massflow=NULL, 
    *dbl_array=NULL, *element_dof_npoint=NULL, *element_dof_npoint_new=NULL,
    *nonlocal_element_info=NULL;

  db( ELEMENT_GROUP, element, &element_group, ddum, ldum, 
    VERSION_NORMAL, GET_IF_EXISTS );
  db( GROUP_MATERI_MEMORY, element_group, &memory, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );

  db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET );

  if ( db( GROUP_TIME, element_group, idum, group_time, ldum, 
      VERSION_NORMAL, GET_IF_EXISTS ) ) {
    birth = group_time[0];
    death = group_time[1];
    db( TIME_CURRENT, 0, idum, &time_current, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    time_total = time_current + dtime;
    if      ( time_total==0. && birth==0. ) 
      ;
    else if ( time_total>birth && time_total<=death ) 
      ;
    else
      return;
  }

  swit = set_swit(element,-1,"elem");
  if ( swit ) pri( "In routine ELEM" );

  length = db_len( ELEMENT, element, VERSION_NORMAL );
  nnol = length - 1;

  n = MTYPE + nnol+1 + nnol + nuknwn + nuknwn +
    2*nnol*nprinc*nnol*nprinc + 4*nnol*nprinc*nnol*nprinc;
  int_array = get_new_int(n);
                                           
  indx = 0;
  types = &int_array[indx]; indx += MTYPE;
  el = &int_array[indx]; indx += nnol+1;
  nodes = &int_array[indx]; indx += nnol;
  dof_type = &int_array[indx]; indx += nuknwn;
  dof_principal = &int_array[indx]; indx += nuknwn;
  element_matrix_unknowns = &int_array[indx]; indx += 2*nnol*nprinc*nnol*nprinc;
  group_matrix_unknowns = &int_array[indx]; indx += 4*nnol*nprinc*nnol*nprinc;
  assert( indx<=n );

  db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
  name = el[0]; nnol = length - 1; array_move( &el[1], nodes, nnol );
  if ( swit ) {
    pri( "element", element );
    pri( "name", name );
    pri( "nodes", nodes, nnol );
  }

  n = nnol*ndim + nnol*ndim + MPOINT*nnol + 
    ndim + nnol*nuknwn +
    nnol*nuknwn + nuknwn + nuknwn + 
    ndim*nuknwn + ndim*nuknwn +
    ndim*nuknwn + ndim*nuknwn + MPOINT + nprinc + nnol*npuknwn +
    nprinc + nnol*npuknwn + nnol + nnol*nuknwn + nnol*npuknwn +
    nnol*npuknwn + nnol*npuknwn + nnol*npuknwn + nnol*npuknwn + nnol*npuknwn +
    ndim + nnol*nprinc*nnol*nprinc + 
    nnol*npuknwn*nnol*npuknwn + MPOINT*MSTRAIN*nnol*ndim +
    MPOINT*ndim*nnol + MPOINT*ndim*nnol +
    nnol*ndim + ndim*ndim +
    nnol*nprinc*nnol*nprinc + nnol*npuknwn*nnol*npuknwn  + nnol +
    mnolnuknwn + mnolnuknwn + length_nei;
  dbl_array = get_new_dbl(n);

  indx = 0;
  coord = &dbl_array[indx]; indx += nnol*ndim;
  new_coord = &dbl_array[indx]; indx += nnol*ndim;
  h = &dbl_array[indx]; indx += MPOINT*nnol;
  coord_ip = &dbl_array[indx]; indx += ndim;
  old_dof = &dbl_array[indx]; indx += nnol*nuknwn;
  new_dof = &dbl_array[indx]; indx += nnol*nuknwn;
  old_unknowns = &dbl_array[indx]; indx += nuknwn;
  new_unknowns = &dbl_array[indx]; indx += nuknwn;
  old_grad_old_unknowns = &dbl_array[indx]; indx += ndim*nuknwn;
  old_grad_new_unknowns = &dbl_array[indx]; indx += ndim*nuknwn;
  new_grad_old_unknowns = &dbl_array[indx]; indx += ndim*nuknwn;
  new_grad_new_unknowns = &dbl_array[indx]; indx += ndim*nuknwn;
  volume = &dbl_array[indx]; indx += MPOINT;
  force_element_volume = &dbl_array[indx]; indx += nprinc;
  tendon_element_rhside = &dbl_array[indx]; indx += nnol*npuknwn;
  residue_factor = &dbl_array[indx]; indx += nprinc;
  element_residue = &dbl_array[indx]; indx += nnol*npuknwn;
  element_dof_tmp_factor = &dbl_array[indx]; indx += nnol;
  element_dof_tmp = &dbl_array[indx]; indx += nnol*nuknwn;
  element_lhside = &dbl_array[indx]; indx += nnol*npuknwn;
  element_rhside = &dbl_array[indx]; indx += nnol*npuknwn;
  element_matrix_delete = &dbl_array[indx];  indx += nnol*npuknwn;
  element_rhside_delete = &dbl_array[indx]; indx += nnol*npuknwn;
  condif_flow = &dbl_array[indx]; indx += ndim;
  element_matrix_values = &dbl_array[indx]; indx += nnol*nprinc*nnol*nprinc;
  element_matrix = &dbl_array[indx]; indx += nnol*npuknwn*nnol*npuknwn;
  new_b = &dbl_array[indx]; indx += MPOINT*MSTRAIN*nnol*ndim;
  old_d = &dbl_array[indx]; indx += MPOINT*ndim*nnol;
  new_d = &dbl_array[indx]; indx += MPOINT*ndim*nnol;
  massflow = &dbl_array[indx]; indx += nnol*ndim;
  grad_massflow = &dbl_array[indx]; indx += ndim*ndim;
  element_matrix_second_values = &dbl_array[indx]; indx += nnol*nprinc*nnol*nprinc;
  element_matrix_second = &dbl_array[indx]; indx += nnol*npuknwn*nnol*npuknwn;
  element_softvar_tmp = &dbl_array[indx]; indx += nnol;
  element_dof_npoint = &dbl_array[indx]; indx += mnolnuknwn;
  element_dof_npoint_new = &dbl_array[indx]; indx += mnolnuknwn;
  nonlocal_element_info = &dbl_array[indx]; indx += length_nei;
  assert( indx<=n );

    // initialize
  array_set( tendon_element_rhside, 0., nnol*npuknwn );
  array_set( element_rhside, 0., nnol*npuknwn );
  array_set( element_matrix_delete, 0., nnol*npuknwn );
  array_set( element_rhside_delete, 0., nnol*npuknwn );
  array_set( element_lhside, 0., nnol*npuknwn );
  array_set( element_residue, 0., nnol*npuknwn );
  array_set( element_dof_tmp, 0., nnol*nuknwn );
  array_set( element_dof_tmp_factor, 0., nnol );
  array_set( element_softvar_tmp, 0., nnol );
  array_set( condif_flow, 0., ndim );
  array_set( massflow, 0., nnol*ndim );
  array_set( grad_massflow, 0., ndim*ndim );
  array_set( element_matrix, 0., nnol*npuknwn*nnol*npuknwn );
  array_set( element_matrix_second, 0., nnol*npuknwn*nnol*npuknwn );
  array_set( element_dof_npoint, 0., mnolnuknwn );
  array_set( element_dof_npoint_new, 0., mnolnuknwn );
  array_set( nonlocal_element_info, 0., length_nei );

  db( GROUP_AXISYMMETRIC, element_group, &axisymmetric, ddum, ldum, 
    VERSION_NORMAL, GET_IF_EXISTS );
  db_max_index( NODE, max_node, VERSION_NORMAL, GET );
  db( DOF_TYPE, 0, dof_type, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( DOF_PRINCIPAL, 0, dof_principal, ddum, ldum, VERSION_NORMAL, GET );
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  db( OPTIONS_SOLVER, 0, &options_solver, ddum, 
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_OPTIONS_SOLVER, icontrol, &options_solver, ddum, 
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( GROUP_MATERI_MEMBRANE, element_group, &membrane, ddum, ldum, 
    VERSION_NORMAL, GET_IF_EXISTS );
  db( GROUP_TYPE, element_group, types, ddum, len_type, 
    VERSION_NORMAL, GET_IF_EXISTS );
  if ( len_type<=0 ) goto skip_element;

  options_matrix_length = nnol*nprinc*nnol*nprinc;
  db( OPTIONS_MATRIX_LENGTH, 0, &options_matrix_length, ddum, ldum, 
    VERSION_NORMAL, GET_IF_EXISTS );

  if ( db_active_index( ELEMENT_RHSIDE_DELETE, element, VERSION_NORMAL ) ) {
    db( ELEMENT_MATRIX_DELETE, element, idum, element_matrix_delete, ldum, 
      VERSION_NORMAL, GET );
    db( ELEMENT_RHSIDE_DELETE, element, idum, element_rhside_delete, ldum, 
      VERSION_NORMAL, GET );
    db( ELEMENT_DELETE_FACTOR, element, idum, &element_delete_factor, ldum, 
      VERSION_NORMAL, GET );
    for ( inol=0; inol<nnol; inol++ ) {
      for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
        iuknwn = ipuknwn*nder;
        if ( dof_principal[iuknwn]>=0 ) {
          indx = inol*npuknwn+ipuknwn;
          element_matrix[indx*nnol*npuknwn+indx] = 
            element_delete_factor * element_matrix_delete[indx];
          element_rhside[indx] = 
            element_delete_factor * element_rhside_delete[indx];
	  
        }
      }
    }
    if(find_nonlocal_weights && scalar_dabs(options_nonlocal_softvar)>TINY) return;
    goto add_to_system_vectors;
  }
    // get dof of element
  for ( inol=0; inol<nnol; inol++ ) {
    inod = nodes[inol];
    indx = inod*npuknwn;
    array_move( db_dbl( NODE_DOF, inod, VERSION_NORMAL ), 
      &old_dof[inol*nuknwn], nuknwn );
    array_move( db_dbl( NODE_DOF, inod, VERSION_NEW ), 
      &new_dof[inol*nuknwn], nuknwn );
    if ( materi_displacement ) {
      for ( idim=0; idim<ndim; idim++ ) {
        new_dof[inol*nuknwn+dis_indx+idim*nder] = 
          old_dof[inol*nuknwn+dis_indx+idim*nder] +
          new_dof[inol*nuknwn+vel_indx+idim*nder] * dtime;   
      }
    }
  }
  	// added for options_element_dof
  if(options_element_dof==-YES) {
    db( ELEMENT_DOF, element, idum, element_dof_npoint, mnolnuknwn, VERSION_NORMAL, GET );		
    if(materi_plasti_softvar_nonlocal && !find_local_softvar)
    	db( ELEMENT_DOF, element, idum, element_dof_npoint_new, mnolnuknwn, VERSION_NEW, GET );		
    db( ELEMENT_DOF_INITIALISED, element, &element_dof_initialised, 
    		ddum, one, VERSION_NORMAL, GET );
  }
  if ( swit ) {
    pri( "old_dof", old_dof, nnol, nuknwn );
    pri( "new_dof", new_dof, nnol, nuknwn );
  }

  if ( materi_diffusion ) {
    materi_diffusion_minimum = EPS_MATERI_DIFFUSION_MINIMUM;
    db( MATERI_DIFFUSION_MINIMUM, 0, idum,
      &materi_diffusion_minimum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    for ( inol=0; inol<nnol; inol++ ) {
      if ( new_dof[inol*nuknwn+diff_indx] >= materi_diffusion_minimum )
        filled++;
    }
    if (        filled == 0    ) {
      element_empty = -YES;
    } else if ( filled == nnol ) {
      element_empty = -NO;
    } else {
      element_empty = -FRONT;
    }                                          
    length = 1;
    db( ELEMENT_EMPTY, element, &element_empty, ddum,
      length, VERSION_NEW, PUT );
    if ( element_empty==-YES ) goto skip_element;
  }

  if ( materi_density ) {
    materi_density_minimum = EPS_MATERI_DENSITY_MINIMUM;
    db( MATERI_DENSITY_MINIMUM, 0, idum,
      &materi_density_minimum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    for ( inol=0; inol<nnol; inol++ ) {
      if ( new_dof[inol*nuknwn+dens_indx] >= materi_density_minimum )
        filled++;
    }
    if (        filled == 0    ) {
      element_empty = -YES;
    } else if ( filled == nnol ) {
      element_empty = -NO;
    } else {
      element_empty = -FRONT;
    }                                          
    length = 1;
    db( ELEMENT_EMPTY, element, &element_empty, ddum,
      length, VERSION_NEW, PUT );
    if ( element_empty==-YES ) goto skip_element;
  }

  array_set( new_unknowns, 0., nuknwn );
  for ( inol=0; inol<nnol; inol++ ) {
    array_add( &new_dof[inol*nuknwn], new_unknowns, new_unknowns, nuknwn );
  }
  array_multiply( new_unknowns, new_unknowns, 1./nnol, nuknwn );
  materi_dens = get_materi_density( element, element_group, nnol, nodes, new_unknowns );

    // get element node coordinates
  for ( inol=0; inol<nnol; inol++ ) {
    inod = nodes[inol];
    db( NODE, inod, idum, &coord[inol*ndim], ldum, VERSION_NORMAL, GET );
    db( NODE, inod, idum, &new_coord[inol*ndim], ldum, VERSION_NEW, GET );
  }

  if ( db_active_index( OPTIONS_RESIDUEFACTOR, 0, VERSION_NORMAL ) )
    db( OPTIONS_RESIDUEFACTOR, 0, idum, residue_factor, ldum, VERSION_NORMAL, GET );
  else array_set( residue_factor, 1., nprinc );

    // lumped relaxation
  db( CONTROL_RELAXATION_MATERI_VELOCITY, icontrol, idum, 
    &control_relaxation_materi_velocity, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_RELAXATION_CONDIF_TEMPERATURE, icontrol, idum, 
    &control_relaxation_condif_temperature, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_RELAXATION_GROUNDFLOW_PRESSURE, icontrol, idum, 
    &control_relaxation_groundflow_pressure, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_RELAXATION_WAVE_FSCALAR, icontrol, idum, 
    &control_relaxation_wave_fscalar, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_RELAXATION_MAXWELL_E, icontrol, idum, 
    &control_relaxation_maxwell_e, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  for ( inol=0; inol<nnol; inol++ ) {
    for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
      iuknwn = ipuknwn*nder;
      indx = inol*npuknwn + ipuknwn;
      if      ( dof_type[iuknwn]==-MATERI_VELOCITY ) {
        element_lhside[indx] = control_relaxation_materi_velocity * dtime;
        element_matrix[indx*nnol*npuknwn+indx] = 
          control_relaxation_materi_velocity * dtime;
      }
      else if ( dof_type[iuknwn]==-WAVE_FSCALAR ) {
        element_lhside[indx] = control_relaxation_wave_fscalar * dtime;
        element_matrix[indx*nnol*npuknwn+indx] = 
          control_relaxation_wave_fscalar * dtime;
      }
      else if ( dof_type[iuknwn]==-CONDIF_TEMPERATURE ) {
        element_lhside[indx] = control_relaxation_condif_temperature;
        element_matrix[indx*nnol*npuknwn+indx] = 
          control_relaxation_condif_temperature;
      }
      else if ( dof_type[iuknwn]==-GROUNDFLOW_PRESSURE ) {
        element_lhside[indx] = control_relaxation_groundflow_pressure;
        element_matrix[indx*nnol*npuknwn+indx] = 
          control_relaxation_groundflow_pressure;
      }
      else if ( dof_type[iuknwn]==-MAXWELL_EI ) {
        element_lhside[indx] = control_relaxation_maxwell_e;
        element_matrix[indx*nnol*npuknwn+indx] = 
          control_relaxation_maxwell_e;
      }
      else if ( dof_type[iuknwn]==-MAXWELL_ER ) {
        element_lhside[indx] = control_relaxation_maxwell_e;
        element_matrix[indx*nnol*npuknwn+indx] = 
          control_relaxation_maxwell_e;
      }
    }
  }

    // structural elements
  if      ( name==-BEAM ) {
    beam_3d( element, element_group, coord, old_dof, 
      new_dof, element_lhside, element_matrix, element_rhside );
    if(find_nonlocal_weights && scalar_dabs(options_nonlocal_softvar)>TINY) return;
    goto add_to_system_vectors;
  }
  else if ( name==-CONTACTSPRING ) {
    contactspring( element, name, element_group, 
      nnol, nodes, coord, old_dof,
      new_dof, element_lhside, element_matrix, element_rhside );
    if(find_nonlocal_weights && scalar_dabs(options_nonlocal_softvar)>TINY) return;
    goto add_to_system_vectors;
  }
  else if ( name==-SPRING1 || name==-SPRING2 ) {
    spring( element, name, element_group, coord, old_dof,
      new_dof, element_lhside, element_matrix, element_rhside );
    if(find_nonlocal_weights && scalar_dabs(options_nonlocal_softvar)>TINY) return;
    goto add_to_system_vectors;
  }
  else if ( name==-TRUSS ) {
    truss( element, element_group, coord, old_dof, 
      new_dof, element_lhside, element_matrix, element_rhside );
    if(find_nonlocal_weights && scalar_dabs(options_nonlocal_softvar)>TINY) return;
    goto add_to_system_vectors;
  }
  else if ( name==-TRUSSBEAM ) {
    beam_3d( element, element_group, coord, old_dof, 
      new_dof, element_lhside, element_matrix, element_rhside );
    truss( element, element_group, coord, old_dof, 
      new_dof, element_lhside, element_matrix, element_rhside );
    if(find_nonlocal_weights && scalar_dabs(options_nonlocal_softvar)>TINY) return;
    goto add_to_system_vectors;
  }

    // update node coordinates for total lagrange for isoparametric elements
  for ( inol=0; inol<nnol; inol++ ) {
    inod = nodes[inol];
    if ( materi_displacement && memory!=-TOTAL_LINEAR ) {
      for ( idim=0; idim<ndim; idim++ ) new_coord[inol*ndim+idim] +=
        new_dof[inol*nuknwn+dis_indx+idim*nder];
    }
  }
  if ( swit ) {
    pri( "coord", coord, nnol, ndim );
    pri( "new_coord", new_coord, nnol, ndim );
  }

    // mass flow
  if ( materi_velocity && materi_density ) {
    for ( inol=0; inol<nnol; inol++ ) {
      dens = new_dof[inol*nuknwn+dens_indx];
      for ( idim=0; idim<ndim; idim++ ) {
        vel = new_dof[inol*nuknwn+vel_indx+idim*nder];
        indx = inol*ndim+idim;
        massflow[indx] = dens * vel;
        if ( axisymmetric==-YES ) {
          radius = new_coord[inol*ndim];
          massflow[indx] *= 2. * PIRAD * radius;
        }
      }
    }
    if ( swit ) pri( "massflow", massflow, nnol, ndim );
  }

  group_materi_plasti_boundary_evaluate( nodes, nnol, 
    element_group, plasti_on_boundary );

    // polynomials and integration point volumes
  pol( element, element_group, name, nnol, coord, new_coord, 
    npoint, h, old_d, new_d, new_b, volume );
  if ( swit ) {
    pri( "volume", volume, npoint );
    pri( "h", h, npoint, nnol );
    pri( "old_d", old_d, npoint*ndim, nnol );
    pri( "new_d", new_d, npoint*ndim, nnol );
  }

  //added for nonlocal_element_info
  if(find_nonlocal_weights && scalar_dabs(options_nonlocal_softvar)>TINY) {
    array_set(nonlocal_element_info, 0, length_nei);
    db( NONLOCAL_ELEMENT_INFO, element, idum, nonlocal_element_info, 
  	  length_nei, VERSION_NORMAL, GET );		
    nonlocal_element_info[0]=npoint;
    for ( ipoint=0; ipoint<npoint; ipoint++ ) {
        // coordinate of integration point
      matrix_ab( &h[ipoint*nnol], new_coord, coord_ip, 1, nnol, ndim );
        // volume factor
      volume_factor( element_group, coord_ip, volfac );
        // volume integration point
      volumeip = volfac*volume[ipoint];
      for(int idim=0; idim<ndim; idim++) 
	      nonlocal_element_info[1+ipoint*ndim+idim]=coord_ip[idim];
      nonlocal_element_info[1+npoint*ndim+ipoint]=volumeip;
    }
    db( NONLOCAL_ELEMENT_INFO, element, idum, nonlocal_element_info, 
  	  length_nei, VERSION_NORMAL, PUT );		

    delete[] int_array;
    delete[] dbl_array;

    return;
  }
  
  if(options_element_dof==-YES) {
    if(npointmax<npoint) {
       cout<<"set 'number_of_intagration_points' at least "<<npoint<<" in inicialization part'"<<endl;
       exit_tn_on_error(); 
    }
  }
  
    // loop over integration points
  for ( ipoint=0; ipoint<npoint; ipoint++ ) {
    if ( swit ) pri( "ipoint", ipoint );

      // old_unknowns
    matrix_ab( &h[ipoint*nnol], old_dof, old_unknowns, 1, nnol, nuknwn );
    if ( swit ) pri( "old_unknowns", old_unknowns, nuknwn );

      // new_unknowns
    matrix_ab( &h[ipoint*nnol], new_dof, new_unknowns, 1, nnol, nuknwn );
    if ( swit ) pri( "new_unknowns", new_unknowns, nuknwn );
    
	// added for options_element_dof
	// uses remembered values of dofs in integration point
	// whenever mesh changes, dofs in i. points are calculated from nodes
    if(options_element_dof==-YES) {
    	long int startindx=nuknwn*ipoint;
	if(element_dof_initialised) {
          for ( int i=0; i<nuknwn; i++ ) 
	    if((i>=hisv_indx && i<(hisv_indx+materi_history_variables))||
	      (i>=stres_indx && i<(stres_indx+MDIM*MDIM))||
	      (i>=epe_indx && i<(epe_indx+MDIM*MDIM))||
   	      (i>=epp_indx && i<(epp_indx+MDIM*MDIM))||
   	      (i>=ept_indx && i<(ept_indx+MDIM*MDIM))||
   	      (i>=epi_indx && i<(epi_indx+MDIM*MDIM))||
	      i==svloc_indx || i==svnonloc_indx) 
	       old_unknowns[i]=element_dof_npoint[i+startindx];
	}
	if( materi_plasti_softvar_nonlocal && materi_plasti_softvar_local &&
			!find_local_softvar ) {  
		new_unknowns[svnonloc_indx]=element_dof_npoint_new[svnonloc_indx+startindx];
		new_unknowns[svloc_indx]=element_dof_npoint_new[svloc_indx+startindx];
	}
	if( !(scalar_dabs(options_nonlocal_softvar)>TINY) ) {
		if( materi_plasti_softvar_nonlocal) 
			new_unknowns[svnonloc_indx]=old_unknowns[svnonloc_indx]; //No change
		if( materi_plasti_softvar_local) 
			new_unknowns[svloc_indx]=old_unknowns[svloc_indx]; 
	}
    }
        
      // gradient of unknowns
    matrix_ab( &old_d[ipoint*ndim*nnol], old_dof, old_grad_old_unknowns,
      ndim, nnol, nuknwn );
    matrix_ab( &old_d[ipoint*ndim*nnol], new_dof, old_grad_new_unknowns,
      ndim, nnol, nuknwn );
    matrix_ab( &new_d[ipoint*ndim*nnol], old_dof, new_grad_old_unknowns,
      ndim, nnol, nuknwn );
    matrix_ab( &new_d[ipoint*ndim*nnol], new_dof, new_grad_new_unknowns,
      ndim, nnol, nuknwn );
    if ( swit ) {
      pri( "old_grad_old_unknowns", old_grad_old_unknowns, ndim, nuknwn );
      pri( "old_grad_new_unknowns", old_grad_new_unknowns, ndim, nuknwn );
      pri( "new_grad_old_unknowns", new_grad_old_unknowns, ndim, nuknwn );
      pri( "new_grad_new_unknowns", new_grad_new_unknowns, ndim, nuknwn );
    }

      // gradient of conservation variables
    if ( materi_velocity && materi_density ) {
      matrix_ab( &new_d[ipoint*ndim*nnol], massflow, 
        grad_massflow, ndim, nnol, ndim );
      if ( swit ) pri( "grad_massflow", grad_massflow, ndim, ndim );
    }

      // coordinate of integration point
    matrix_ab( &h[ipoint*nnol], new_coord, coord_ip, 1, nnol, ndim );
    if ( swit ) pri( "coord_ip", coord_ip, ndim );

      // volume factor
    volume_factor( element_group, coord_ip, volfac );
    if ( swit ) pri( "volfac", volfac );

      // volume integration point
    if ( materi_diffusion )
      volumeip = volfac*volume[ipoint]*new_unknowns[diff_indx];
    else
      volumeip = volfac*volume[ipoint];

    if ( membrane==-YES && materi_strain_total ) {
      if ( ndim==2 )
        volumeip *= ( 1. + new_unknowns[ept_indx+2*nder] );
      else {
        assert( ndim==1 );
        volumeip *= ( 1. + new_unknowns[ept_indx+1*nder] ) *
          ( 1. + new_unknowns[ept_indx+2*nder] );
      }
    }

      // total element volume
    element_volume += volumeip;

      // element force
    force_element_volume_set( element, nnol, nodes, coord_ip, force_element_volume );
    for ( inol=0; inol<nnol; inol++ ) {
      iprinc=0;
      for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
        iuknwn = ipuknwn*nder;
        if ( dof_principal[iuknwn]>=0 ) {
          indx = inol*npuknwn + ipuknwn;
          tmp = h[ipoint*nnol+inol] * force_element_volume[iprinc];
          element_rhside[indx] += volumeip * tmp;
          if ( residue ) element_residue[indx] -= 
            h[ipoint*nnol+inol] * force_element_volume[iprinc];
          iprinc++;
        }
      }
    }

      // loop over differential operators
    for ( itype=0; itype<len_type; itype++ ) {
      type = types[itype];
      general( element, name, nnol, element_group, type, 
        nodes, coord_ip, old_dof, new_dof, 
        old_unknowns, new_unknowns, 
        new_grad_new_unknowns, &h[ipoint*nnol], 
        &new_d[ipoint*ndim*nnol], volumeip,
        grad_massflow, element_rhside, element_residue, 
        element_lhside, element_matrix );
      if      ( type==-CONDIF )
        condif( element, element_group, nnol, 
          &h[ipoint*nnol], volumeip,
          new_unknowns, element_lhside, element_matrix, element_rhside, 
          element_residue );
      else if ( type==-GROUNDFLOW )
        groundflow( element, element_group, nnol, 
          coord_ip, &h[ipoint*nnol], 
          &new_d[ipoint*ndim*nnol], volumeip, old_unknowns, 
          new_unknowns, new_grad_new_unknowns,
          element_matrix, element_rhside, element_residue );
      else if ( type==-MATERI ) 
        materi( element, element_group, nnol, npoint, 
          nodes, plasti_on_boundary, coord_ip, coord, 
          &h[ipoint*nnol], &new_d[ipoint*ndim*nnol], 
          &new_b[ipoint*MSTRAIN*nnol*ndim], 
          volumeip, old_unknowns, new_unknowns, old_grad_old_unknowns, 
          old_grad_new_unknowns, new_grad_new_unknowns,
          element_lhside, element_matrix,
          element_rhside, element_residue, tendon_element_rhside );
      else if ( type==-MAXWELL_FREQUENCY ||  type==-MAXWELL_TIME )
        maxwell( type, element, element_group, nnol, volumeip,
          new_unknowns, old_dof, new_dof,
          &h[ipoint*nnol], &new_d[ipoint*ndim*nnol], 
          element_lhside, element_matrix, element_matrix_second,
          element_rhside );
      else if ( type==-WAVE )
          wave( element, element_group, nnol, 
          &h[ipoint*nnol], &new_d[ipoint*ndim*nnol],
          volumeip, new_unknowns, new_grad_new_unknowns, 
          element_lhside, element_matrix, 
          element_rhside, element_residue );
      else if ( type!=-EMPTY && type!=-NONE ) {
        pri( "Error detected for element ", element );
        pri( "Illegal group_type  ", type );
        exit_tn_on_error();
      }
    }

    if(!find_local_softvar) {//Not used when searching for local values of softening variable 
    
      // contribution to unknown gradients in nodes
    if ( derivatives ) {
      for ( inol=0; inol<nnol; inol++ ) {
        for ( idim=0; idim<ndim; idim++ ) {
          for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
            iuknwn = ipuknwn*nder + idim + 1;
            element_dof_tmp[inol*nuknwn+iuknwn] += h[ipoint*nnol+inol] * 
              new_grad_new_unknowns[idim*nuknwn+ipuknwn*nder];
          }
        }
        element_dof_tmp_factor[inol] += h[ipoint*nnol+inol];
      }
    }

      // contribution to element strain energy
    if ( materi_strainenergy ) {
      element_strainenergy += volumeip * new_unknowns[ener_indx];
    }
	// added for options_element_dof
    if(options_element_dof==-YES) {
    	long int startindx=nuknwn*ipoint;
        for ( int i=0; i<nuknwn; i++ ) 
 	  if((i>=hisv_indx && i<(hisv_indx+materi_history_variables))||
	    (i>=stres_indx && i<(stres_indx+MDIM*MDIM))||
	    (i>=epe_indx && i<(epe_indx+MDIM*MDIM))||
   	    (i>=epp_indx && i<(epp_indx+MDIM*MDIM))||
   	    (i>=ept_indx && i<(ept_indx+MDIM*MDIM))||
   	    (i>=epi_indx && i<(epi_indx+MDIM*MDIM))||
	           i==svloc_indx || i==svnonloc_indx)
	    element_dof_npoint[i+startindx]=new_unknowns[i];//only for unknowns filled in materi.cc
    }
    } //Not used when searching for local values of softening variable - end
    if(find_local_softvar) {
    	long int startindx=nuknwn*ipoint;
	element_dof_npoint[svloc_indx+startindx]=new_unknowns[svloc_indx];
	db( ELEMENT_DOF, element, idum, element_dof_npoint, mnolnuknwn, VERSION_NEW, PUT );		
    }
  }
  //end loop over integration points

  if(!find_local_softvar) {//Not used when searching for local values of softening variable 

  //added for options_element_dof 
  if(options_element_dof==-YES) 
    db( ELEMENT_DOF, element, idum, element_dof_npoint, mnolnuknwn, VERSION_NEW, PUT );		
  
  array_multiply( tendon_element_rhside, tendon_element_rhside,
    1./npoint, nnol*npuknwn );
  array_add( tendon_element_rhside, element_rhside, 
    element_rhside, nnol*npuknwn );
  if ( derivatives ) {
    for ( inol=0; inol<nnol; inol++ ) {
      indx = inol * nuknwn;
      if ( element_dof_tmp_factor[inol]!=0. ) 
        array_multiply( &element_dof_tmp[indx], &element_dof_tmp[indx],
        1./element_dof_tmp_factor[inol], nuknwn );
    }
  }

    // area integrals
  area( element, name, element_group, nnol, nodes, new_coord, new_dof, 
    element_lhside, element_matrix, element_rhside );

    // residue
  if ( residue ) {
    for ( inol=0; inol<nnol; inol++ ) {
      ipuknwn = res_indx/nder;
      res = element_residue[inol*npuknwn+ipuknwn];
      iprinc = 0;
      for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
        iuknwn = ipuknwn * nder;
        if ( dof_principal[iuknwn]>=0 ) {
  	      res += residue_factor[iprinc] * 
            scalar_dabs(element_residue[inol*npuknwn+ipuknwn]);
          iprinc++;
        }
      }
      ipuknwn = res_indx/nder;
      iuknwn = ipuknwn * nder;
      indx = inol*npuknwn + ipuknwn;
      element_rhside[indx] += ( res - new_dof[inol*nuknwn+iuknwn] ) / dtime;
      element_lhside[indx] += 1. / dtime;
    }
  }

    // store elements results
  length = 1;
  db( ELEMENT_VOLUME, element, idum, &element_volume, length, VERSION_NORMAL, PUT );
  if ( materi_velocity ) {
    element_mass = materi_dens * element_volume;
    db( ELEMENT_MASS, element, idum, &element_mass, length, VERSION_NORMAL, PUT );
  }
  if ( materi_strainenergy ) {
    db( ELEMENT_STRAINENERGY, element, idum, &element_strainenergy, 
      length, VERSION_NORMAL, PUT );
  }
  }//Not used when searching for local values of softening variable - end
  
  add_to_system_vectors:
  
  if(!find_local_softvar) {//Not used when searching for local values of softening variable 

    // store element matrices sparse
  if ( options_solver!=-DIAGONAL && options_solver!=-NONE ) {

    length_element_matrix_values = 0;
    for ( inol=0; inol<nnol; inol++ ) {
      inod = nodes[inol];
      for ( jnol=0; jnol<nnol; jnol++ ) {
        jnod = nodes[jnol];
        for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
          iuknwn = ipuknwn*nder;
          for ( jpuknwn=0; jpuknwn<npuknwn; jpuknwn++ ) {
            juknwn = jpuknwn*nder;
            indxi = inol*npuknwn + ipuknwn;
            indxj = jnol*npuknwn + jpuknwn;
            if ( dof_principal[iuknwn]>=0 && dof_principal[juknwn]>=0 ) {
              tmp1 = element_matrix[indxi*nnol*npuknwn+indxj];
              tmp2 = element_matrix_second[indxi*nnol*npuknwn+indxj];
              test1 = scalar_dabs(tmp1)>EPS_MATRIX;
              test2 = scalar_dabs(tmp2)>EPS_MATRIX;
              if ( test1 || test2 ) {
                element_matrix_values[length_element_matrix_values] = tmp1;
                element_matrix_second_values[length_element_matrix_values] = tmp2;
                element_matrix_unknowns[2*length_element_matrix_values+0] = 
                  inod*npuknwn + ipuknwn;
                element_matrix_unknowns[2*length_element_matrix_values+1] = 
                  jnod*npuknwn + jpuknwn;
                group_matrix_unknowns[4*length_element_matrix_values+0] = inol;
                group_matrix_unknowns[4*length_element_matrix_values+1] = ipuknwn;
                group_matrix_unknowns[4*length_element_matrix_values+2] = jnol;
                group_matrix_unknowns[4*length_element_matrix_values+3] = jpuknwn;
                length_element_matrix_values++;
              }
            } 
          }
        }
      }
    }
    if ( length_element_matrix_values>0 ) {
      if ( length_element_matrix_values>options_matrix_length ) {
        pri( "Error: too small value for OPTIONS_MATRIX_LENGTH for element", element );
        pri( "Set OPTIONS_MATRIX_LENGTH at least to", length_element_matrix_values );
        exit(TN_EXIT_STATUS);
      }
      db( OPTIONS_MATRIX_GROUP, 0, &options_matrix_group, ddum, ldum, 
        VERSION_NORMAL, GET_IF_EXISTS );
      length = length_element_matrix_values;
      if ( options_matrix_group==-YES ) {
        if ( axisymmetric==-YES ) {
          pri( "Error: OPTIONS_MATRIX_GROUP cannot be used icw axisymmetry." );
          exit_tn_on_error();
        }
        parallel_sys_lock();
        if ( !db_active_index( GROUP_MATRIX_VALUES, element_group, VERSION_NORMAL) ) {
          db( GROUP_MATRIX_VALUES, element_group, idum, element_matrix_values, 
            length, VERSION_NORMAL, PUT );
          if ( eigen_active ) {
            db( GROUP_MATRIX_SECOND_VALUES, element_group, idum, 
              element_matrix_second_values, length, VERSION_NORMAL, PUT );
          }
          length = 4*length_element_matrix_values;
          db( GROUP_MATRIX_UNKNOWNS, element_group, group_matrix_unknowns, ddum,
          length, VERSION_NORMAL, PUT );
        }
        parallel_sys_unlock();
      }
      else {
        db( ELEMENT_MATRIX_VALUES, element, idum, element_matrix_values, 
          length, VERSION_NORMAL, PUT );
        if ( eigen_active ) {
          db( ELEMENT_MATRIX_SECOND_VALUES, element, 
            idum, element_matrix_second_values, length, VERSION_NORMAL, PUT );
        }
        length = 2*length_element_matrix_values;
        db( ELEMENT_MATRIX_UNKNOWNS, element, element_matrix_unknowns, ddum,
          length, VERSION_NORMAL, PUT );
      }
    }
  }

    // add
  for ( inol=0; inol<nnol; inol++ ) {
    inod = nodes[inol];
    indx = ithread*(1+max_node)+inod;
    node_nel[indx]++;
    indx = ithread*(1+max_node)*npuknwn+inod*npuknwn;
    array_add( &node_rhside[indx], &element_rhside[inol*npuknwn], 
      &node_rhside[indx], npuknwn );
    for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
      iuknwn = ipuknwn*nder;
      if ( options_solver==-DIAGONAL || dof_principal[iuknwn]<0 )
       	node_lhside[indx+ipuknwn] += element_lhside[inol*npuknwn+ipuknwn];
    }
    indx = ithread*(1+max_node)*nuknwn+inod*nuknwn;
    array_add( &node_dof_tmp[indx], &element_dof_tmp[inol*nuknwn], 
      &node_dof_tmp[indx], nuknwn );
  }

    // store for failure/delete options
  if ( db_active_index( ELEMENT_DELETE_FACTOR, element, VERSION_NORMAL ) && 
       !db_active_index( ELEMENT_RHSIDE_DELETE, element, VERSION_NORMAL ) ) {
    for ( inol=0; inol<nnol; inol++ ) {
      for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
        iuknwn = ipuknwn*nder;
        if ( dof_principal[iuknwn]>=0 ) {
          indx = inol*npuknwn+ipuknwn;
          element_rhside_delete[indx] = element_rhside[indx];
          element_matrix_delete[indx] = element_matrix[indx*nnol*npuknwn+indx];
        }
      }
    }
    length = nnol*npuknwn;
    db( ELEMENT_RHSIDE_DELETE, element, idum, element_rhside_delete, 
      length, VERSION_NORMAL, PUT );

    db( ELEMENT_MATRIX_DELETE, element, idum, element_matrix_delete, 
      length, VERSION_NORMAL, PUT );
    if ( swit ) {
      pri( "element_rhside_delete", element_rhside_delete, nnol, npuknwn );
      pri( "element_matrix_delete", element_matrix_delete, nnol, npuknwn );
    }
  }

    // some final printing
  if ( swit ) {
    if ( options_solver!=-DIAGONAL ) {
      pri( "element_matrix", element_matrix, nnol*npuknwn, nnol*npuknwn );
      if ( materi_stress ) pri( "element_matrix_second", 
        element_matrix_second, nnol*npuknwn, nnol*npuknwn );
    }
    pri( "element_lhside", element_lhside, nnol, npuknwn );
    pri( "element_rhside", element_rhside, nnol, npuknwn );
    pri( "element_residue", element_residue, nnol, npuknwn );
  }
  if(find_local_softvar>0) cout<<find_local_softvar<<endl;

  }//Not used when searching for local values of softening variable - end

  skip_element:

  delete[] int_array;
  delete[] dbl_array;

  if ( swit ) pri( "Out routine ELEM" );
}
