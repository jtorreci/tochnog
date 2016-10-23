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

  // this routine is only meant to initialize static variables

long int 
  echo=0, ndim=-1, derivatives=0, 
  beam_rotation=0, condif_temperature=0, 
  groundflow_velocity=0, groundflow_pressure=0, 
  materi_history_variables=0, materi_damage=0, materi_density=0, 
  materi_displacement=0, materi_diffusion=0, materi_maxwell_stress=0, 
  materi_plasti_incremental_substeps=0,
  materi_plasti_kappa=0, materi_strain_intergranular=0,
  materi_plasti_rho=0, materi_plasti_f=0, materi_plasti_f_nonlocal=0,
  materi_plasti_softvar_local=0,
  materi_plasti_softvar_nonlocal=0,
  materi_strainenergy=0, materi_strain_elasti=0, 
  materi_strain_plasti=0, materi_strain_total=0,
  materi_stress=0, materi_velocity=0, materi_velocity_integrated=0,
  materi_void_fraction=0, materi_work=0, 
  maxwell_e=0, maxwell_fe=0, maxwell_er=0, maxwell_ei=0,
  residue=0, wave_scalar=0, wave_fscalar=0,
  find_local_softvar=0, find_nonlocal_weights=0, nonlocal_first_set=0;
long int 
  any_runtime=0, nder=1, npuknwn=0, nuknwn=0, npointmax=6, nprinc=0, dam_indx=-1,
  dens_indx=-1, diff_indx=-1, dis_indx=-1, ener_indx=-1,
  epe_indx=-1, epi_indx=-1, epp_indx=-1, ept_indx=-1,
  maxfe_indx=-1, maxe_indx=-1, maxer_indx=-1, maxei_indx=0,
  gvel_indx=-1, hisv_indx=-1, kap_indx=-1, f_indx=-1, fn_indx=-1,
  mstres_indx=-1, pres_indx=-1, res_indx=-1, rho_indx=-1,
  rot_indx=-1, scal_indx=-1, stres_indx=-1, substeps_indx=-1,
  svloc_indx=-1, svnonloc_indx=-1,
  temp_indx=-1, fscal_indx=-1, vel_indx=-1,
  veli_indx=-1, void_indx=-1, work_indx=-1;
char 
  data_file[MCHAR], 
  data_file_base[MCHAR], 
  post_calcul_names[DATA_ITEM_SIZE][MCHAR],
  post_calcul_names_without_extension[DATA_ITEM_SIZE][MCHAR],
  initialization_names[DATA_ITEM_SIZE][MCHAR];
long int 
  map_version_from=0, map_version_to=0, map_always=0;
long int 
  post_found=0, post_node[4], post_node_length=0, npost_node=0;
double 
  post_point[MDIM], post_point_dof[MUKNWN], post_node_result[DATA_ITEM_SIZE];
long int 
  calcul_matrix=0, calcul_vector=0, calcul_ecomplex=0, calcul_scalar_indx=0,
  calcul_operat=0, calcul_mat_indx=0, calcul_vec_indx=0;
long int 
  geometry_ent[2], *nodes_in_geometry=NULL;
long int
  parallel_active=0;
long int
  eigen_active=0;
long int
  hyper_active=0;
long int
  timestep_decrease=0;
long int 
  *solve_global_local, solve_nlocal, solve_options_matrix_group;
double 
  *solve_b, *solve_x;
long int 
  swit_element_stack=-1, swit_node_stack=-1;
char 
  swit_routine_stack[MSTACK][MCHAR];
	// added for options_element_dof
long int options_element_dof=-YES;
double options_nonlocal_softvar=0;


void initialize( void )

{
  long int i=0;

  strcpy( data_file, "" );
  strcpy( data_file_base, "" );
  for ( i=0; i<DATA_ITEM_SIZE; i++ ) {
    strcpy( post_calcul_names[i], "" );
    strcpy( post_calcul_names_without_extension[i], "" );
    strcpy( initialization_names[i], "" );
  }
  for ( i=0; i<MSTACK; i++ ) {
    strcpy( swit_routine_stack[i], "" );
  }

}
