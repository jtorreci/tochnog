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

#define ITERATION_WISHED 4
#define EPS_DTIME 1.e-4

void top( void )

{
  long int inod=0, ielem=0, max=0, max_node=0, max_control=0, 
    start_control=0, length=0, iinc=0, ninc=0, 
    timestep_iterations=1, swit=0, swit_timestep=0, idat=0, 
    inverse_iterations=0, inverse_iter=0, ipar=0, npar=0, 
    ipar_i=0, ipar_n=0, max_ipar=0, max_ipar_i=0, icontrol=0, iteration=0, 
    number_processors=1, number_iterations=0,
    iteration_min=0, iteration_max=0,
    use_control_timestep_iterations_automatic=0, 
    use_control_timestep_size_automatic_decrease=0, 
    control_timestep_iterations_automatic_stop=-YES, 
    converged_once=0, time_at_start=0,
    options_stabilization=-STATIC, options_elementloop=-YES,
    print_control=-NO, options_solver=-MATRIX_ITERATIVE_BICG, print_where=-NO,
    ldum=0, idum[1], options_mesh[MDIM], mnolnuknwn=npointmax*nuknwn, 
    length_nei=1+npointmax*ndim+npointmax+2,
    max_elem=0, one=1, zero=0;
  double tmp=0., time_increment=0., time_current=0., dtime_initial=0.,
    dtime=0., time_old=0., time_new=0., ratio_criterium=0., 
    ratio_max=0., decrease_factor=0., min_timestep=0.,
    maximum_timestep=0., post_node_rhside_ratio=0.,
    multiplier=1., ddum[1], control_timestep_iterations_automatic[2], 
    control_timestep_size_automatic_decrease[3], 
    dwork[MUKNWN], *timestep=NULL, dworkmnol[MPOINT*MUKNWN], *dworknei=NULL;

  set_environment();

  swit = set_swit(-1,-1,"top" );
  if ( swit ) pri( "In routine TOP" );

  if ( swit ) {
    pri( "ndim", ndim );
    pri( "npuknwn", npuknwn );
    pri( "nuknwn", nuknwn );
  }

  time_at_start = (long int) time(NULL);
  length=1; db( TIME_AT_START, 0, &time_at_start, ddum, length, VERSION_NORMAL, PUT );

  if ( materi_velocity ) {
    array_set( options_mesh, -FOLLOW_MATERIAL, ndim );
    db( OPTIONS_MESH, 0, options_mesh, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    if ( materi_displacement ) array_set( options_mesh, -FIXED_IN_SPACE, ndim );
    db( OPTIONS_MESH, 0, options_mesh, ddum, ndim, VERSION_NORMAL, PUT );
  }

  db( PRINT_CONTROL, 0, &print_control, ddum, 
    ldum, VERSION_NORMAL, GET_IF_EXISTS );

  db( PRINT_WHERE, 0, &print_where, ddum, 
    ldum, VERSION_NORMAL, GET_IF_EXISTS );

  db( OPTIONS_ELEMENTLOOP, 0, &options_elementloop, ddum, 
    ldum, VERSION_NORMAL, GET_IF_EXISTS );

  db( OPTIONS_PROCESSORS, 0, &number_processors, ddum, 
    length, VERSION_NORMAL, GET_IF_EXISTS );
  if ( number_processors<1 ) db_error( OPTIONS_PROCESSORS, 0 );
  if ( number_processors>MTHREAD ) number_processors = MTHREAD;
  if ( db_max_index( GROUP_USER_UMAT, max, VERSION_NORMAL, GET ) >=0 )
    number_processors = 1;
  length=1; db( OPTIONS_PROCESSORS, 0, &number_processors, ddum, 
    length, VERSION_NORMAL, PUT );

  if ( db_active_index( OPTIONS_STABILIZATION, 0, VERSION_NORMAL ) ) {
    db( OPTIONS_STABILIZATION, 0, &options_stabilization, ddum, 
      length, VERSION_NORMAL, GET );
    if ( options_stabilization==-DYNAMIC ) {
      pri( "Error: sorry -DYNAMIC not implemented yet for OPTIONS_STABILIZATION." );
      exit(TN_EXIT_STATUS);
    }
  }
  else {
    db( OPTIONS_STABILIZATION, 0, &options_stabilization, ddum, 
      length, VERSION_NORMAL, PUT );
  }
  db( OPTIONS_NONLOCAL_SOFTVAR, 0, idum, &options_nonlocal_softvar, ldum, 
      VERSION_NORMAL, GET_IF_EXISTS );

  parallel_sys_initialize();

    // initialize node_dof, if not specified
  if ( nuknwn>0 ) {
    array_set( dwork, 0., nuknwn );
    db_max_index( NODE, max_node, VERSION_NORMAL, GET );
    for ( inod=0; inod<=max_node; inod++ ) {
      if ( db_active_index( NODE, inod, VERSION_NORMAL ) && 
           !db_active_index( NODE_DOF, inod, VERSION_NORMAL ) )
        db( NODE_DOF, inod, idum, dwork, nuknwn, VERSION_NORMAL, PUT );
    }
  }

    // added for options_element_dof
    // initialize element_dof, if not specified
  db( OPTIONS_ELEMENT_DOF, 0, &options_element_dof, ddum, ldum, 
    VERSION_NORMAL, GET_IF_EXISTS );
  if ( nuknwn>0 ) {
    array_set( dworkmnol, 0., mnolnuknwn );
    db_max_index( ELEMENT, max_elem, VERSION_NORMAL, GET );
    dworknei = get_new_dbl(length_nei);
    array_set( dworknei, 0, length_nei );
    dworknei[length_nei-1]=1;

    for ( ielem=0; ielem<=max_elem; ielem++ ) {
      if ( db_active_index( ELEMENT, ielem, VERSION_NORMAL ) && 
           !db_active_index( ELEMENT_DOF, ielem, VERSION_NORMAL ) ) 
        db( ELEMENT_DOF, ielem, idum, dworkmnol, mnolnuknwn, VERSION_NORMAL, PUT );
      if ( db_active_index( ELEMENT, ielem, VERSION_NORMAL ) && 
           !db_active_index( ELEMENT_DOF_INITIALISED, ielem, VERSION_NORMAL ) ) 
        db( ELEMENT_DOF_INITIALISED, ielem, &zero, ddum, one, VERSION_NORMAL, PUT );
      if ( db_active_index( ELEMENT, ielem, VERSION_NORMAL ) && 
           !db_active_index( NONLOCAL_ELEMENT_INFO, ielem, VERSION_NORMAL ) ) 
        db( NONLOCAL_ELEMENT_INFO, ielem, idum, dworknei, length_nei, VERSION_NORMAL, PUT );		
    }
  }

    // initialize time, if not specified
  if ( !db_active_index( TIME_CURRENT, 0, VERSION_NORMAL ) ) {
    length = 1;
    db( TIME_CURRENT, 0, idum, &time_current, length, VERSION_NORMAL, PUT );
  }

    // determine the highest index of timestep, print, etc..
  for ( idat=0; idat<MDAT; idat++ ) {
    if ( db_data_class(idat)==CONTROL ) {
      db_highest_index( idat, max, VERSION_NORMAL );
      if ( max>max_control ) max_control = max;
    }
  }

    // number of inverse parameters
  db( INVERSE_ITERATIONS, 0, &inverse_iterations, ddum, ldum, 
    VERSION_NORMAL, GET_IF_EXISTS );
  db_max_index( INVERSE_PARAMETER, max, VERSION_NORMAL, GET ); 
  if ( max>=0 ) {
    ipar_n = 3;
    for ( ipar=0; ipar<=max; ipar++ ) {
      if ( db_active_index( INVERSE_PARAMETER, ipar, VERSION_NORMAL ) ) 
        npar = ipar + 1;
    }
  }
  max_ipar = scalar_imax(1,npar);
  max_ipar_i = scalar_imax(1,ipar_n);

    // store initial data, if not specified
  db_version_copy( VERSION_NORMAL, VERSION_START );
  db_max_index( NODE_START_REFINED, max, VERSION_NORMAL, GET ); 
  if ( max<0 ) db_copy( NODE, NODE_START_REFINED, VERSION_NORMAL );
  db_max_index( NODE_DOF_START_REFINED, max, VERSION_NORMAL, GET ); 
  if ( max<0 ) db_copy( NODE_DOF, NODE_DOF_START_REFINED, VERSION_NORMAL );

    // initialize mesh
  mesh_has_changed( VERSION_NORMAL );
  
    // inverse iterations loop
  for ( inverse_iter=0; inverse_iter<=inverse_iterations; inverse_iter++ ) {
    if ( swit ) pri( "inverse_iter", inverse_iter );
    if ( npar>0 ) {
      length = 1;
      db( INVERSE_ITERATION_NUMBER, 0, &inverse_iter, ddum, length, 
        VERSION_NORMAL, PUT );
    }
      // inverse parameters loop
    for ( ipar=0; ipar<max_ipar; ipar++ ) {
      if ( npar==0 || db_active_index( INVERSE_PARAMETER, ipar, VERSION_NORMAL ) ) {
          // inverse central differences loop
        for ( ipar_i=0; ipar_i<max_ipar_i; ipar_i++ ) {
          inverse_calculation( ipar, npar, ipar_i, ipar_n, max_control, 
            INVERSE_PARAMETER_SET );
          start_control=0;
          repeat_point:
          if ( max_control>=0 ) {
              // time loop
            for ( icontrol=start_control; icontrol<=max_control; icontrol++ ) {
              if ( swit ) pri( "icontrol", icontrol );
              if ( db_partialname_any_index( "control_", icontrol ) ||
                   icontrol==max_control ) {
		   
                if ( print_control==-YES ) pri( "control index", icontrol );
                db( TIME_CURRENT, 0, idum, &time_current, ldum, VERSION_NORMAL, GET );
                swit_timestep = swit; swit = swit && set_swit(-1,-1,"");
                length=1; db( ICONTROL, 0, &icontrol, ddum, length, 
                  VERSION_NORMAL, PUT );
		  
                if ( db_active_index( CONTROL_TIMESTEP, icontrol, VERSION_NORMAL ) ) {
                  maximum_timestep = DBL_MAX;
                  ratio_criterium = DBL_MAX;
                  use_control_timestep_iterations_automatic = 0;
                  if ( db_active_index(  CONTROL_TIMESTEP_ITERATIONS_AUTOMATIC, 
                      icontrol, VERSION_NORMAL ) ) {
                    db( CONTROL_TIMESTEP_ITERATIONS_AUTOMATIC, icontrol, 
                      idum, control_timestep_iterations_automatic, ldum, 
                      VERSION_NORMAL, GET );
                    db( CONTROL_TIMESTEP_ITERATIONS_AUTOMATIC_STOP, icontrol, 
                      &control_timestep_iterations_automatic_stop, ddum, ldum, 
                      VERSION_NORMAL, GET_IF_EXISTS );
                    iteration_min = 2;
                    iteration_max = 20;
                    use_control_timestep_iterations_automatic = 1;
                    ratio_criterium = control_timestep_iterations_automatic[0];
                    maximum_timestep = control_timestep_iterations_automatic[1];
                  }
                  else if ( db_active_index( CONTROL_TIMESTEP_ITERATIONS, 
                      icontrol, VERSION_NORMAL ) ) {
                    db( CONTROL_TIMESTEP_ITERATIONS, icontrol, &timestep_iterations, 
                      ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
                    iteration_min = timestep_iterations;
                    iteration_max = timestep_iterations;
                    if ( iteration_max<2 && residue ) {
                      pri( "Error: use  at least 2 iterations with residue." );
                      exit(TN_EXIT_STATUS);
                    }
                  }
                  else {
                    iteration_min = 2;
                    iteration_max = 2;
                  }
                  if ( db_active_index(  CONTROL_TIMESTEP_SIZE_AUTOMATIC_DECREASE, 
                      icontrol, VERSION_NORMAL ) ) {
                    db( CONTROL_TIMESTEP_SIZE_AUTOMATIC_DECREASE, icontrol, 
                      idum, control_timestep_size_automatic_decrease, ldum, 
                      VERSION_NORMAL, GET );
                    use_control_timestep_size_automatic_decrease = 1;
                    ratio_max = control_timestep_size_automatic_decrease[0];
                    decrease_factor = control_timestep_size_automatic_decrease[1];
                    min_timestep = control_timestep_size_automatic_decrease[2];
                  }
                  multiplier = 1.;
                  db( CONTROL_TIMESTEP_MULTIPLIER, icontrol, idum, &multiplier, 
                    ldum, VERSION_NORMAL, GET_IF_EXISTS );
                  timestep = db_dbl( CONTROL_TIMESTEP, icontrol, VERSION_NORMAL );
                  ninc = db_len( CONTROL_TIMESTEP, icontrol, VERSION_NORMAL )/2;
                  if ( ninc==0 ) db_error( CONTROL_TIMESTEP, icontrol );
		  
                    // time increments loop
                  for ( iinc=0; iinc<ninc; iinc++ ) {
                    time_old = time_current;
                    dtime_initial = dtime = timestep[iinc*2];
                    time_increment = timestep[iinc*2+1];
                    time_new = time_old + time_increment;
                    length = 1;
                    db( TIME_OLD, 0, idum, &time_old, length, VERSION_NORMAL, PUT );
                    db( TIME_NEW, 0, idum, &time_new, length, VERSION_NORMAL, PUT );
                    if ( swit ) {
                      pri( "dtime", dtime );
                      pri( "time_old", time_old );
                      pri( "time_increment", time_increment );
                      pri( "time_new", time_new );
                      pri( "timestep_iterations", timestep_iterations );
                    }
		      
                    while ( time_current<time_new && (dtime>EPS_DTIME*dtime_initial) ) {
                      start_of_step:
                      db( TIME_CURRENT, 0, idum, &time_current, length,
                        VERSION_NORMAL, GET );
                      dtime *= multiplier;
                      time_current += dtime;
                      tmp = time_new - time_current;
                      if ( tmp>0. && tmp<0.5*dtime ) {
                        dtime += tmp;
                        time_current = time_new;
                      }
                      if ( time_current>time_new ) {
                        dtime -= time_current-time_new;
                        time_current = time_new;
                      }
                      options_solver = -MATRIX_ITERATIVE_BICG;
		         //options_solver = -MATRIX_SUPERLU;
                      db( OPTIONS_SOLVER, 0, &options_solver, ddum, length,
                        VERSION_NORMAL, GET_IF_EXISTS );
                      db( CONTROL_OPTIONS_SOLVER, icontrol, &options_solver, ddum, length,
                        VERSION_NORMAL, GET_IF_EXISTS );
						
                      step_start( YES, &options_solver, dtime, time_current );
                      db_delete( CONTROL_EIGEN_VALUES, VERSION_NORMAL );
                      db_delete( NODE_EIGEN, VERSION_NORMAL );
                      db_version_copy( VERSION_NORMAL, VERSION_NEW );
                      length=1;
                      db( DTIME, 0, idum, &dtime, length, VERSION_NEW, PUT );
                      db( TIME_CURRENT, 0, idum, &time_current, length,
                        VERSION_NEW, PUT );
                        // equilibrium loop

                      for ( iteration=1; ( iteration<=iteration_min || 
                          post_node_rhside_ratio>=ratio_criterium ) && 
                          iteration<=iteration_max; iteration++ ) {
                        if ( swit ) pri( "iteration", iteration );
                        number_iterations = iteration;
                        db( NUMBER_ITERATIONS, 0, &number_iterations, 
                          ddum, length, VERSION_NEW, PUT );
                        if ( nuknwn>0 ) {
                          iteration_start( );
                          if ( print_where==-YES ) pri( "Where: before boundary conditions." ); 
                          bounda(); 
                          if ( print_where==-YES ) pri( "Where: after boundary conditions." ); 
                          unknown_freeze(); 
                          if ( options_elementloop==-YES ) {
			    if (scalar_dabs(options_nonlocal_softvar)>TINY) {
				  if(!nonlocal_first_set) nonlocal_set();
				  nonlocal_first_set=1;
				  find_local_softvar=1;
	                          element_loop(); 
				  find_local_softvar=0;
				  calc_nonlocal_softvar();
			    }
                            if ( print_where==-YES ) pri( "Where: before element_loop." ); 
                            element_loop(); 
                            if ( print_where==-YES ) pri( "Where: after element_loop." ); 
                          }
                          parallel_sys_routine( &parallel_contact );
                          parallel_sys_routine( &parallel_new_dof_before );
                          slide(); 
                          if ( print_where==-YES ) pri( "Where: before solver." ); 
                          solve( options_solver );
                          if ( print_where==-YES ) pri( "Where: after solver." ); 
                          parallel_sys_routine( &parallel_new_dof_diagonal );
                          locate();
                          nonlocal_apply();
                          post_node_rhside_fixed_free();
                          db( POST_NODE_RHSIDE_RATIO, 0, idum, &post_node_rhside_ratio, 
                            ldum, VERSION_NORMAL, GET );
                        }
                      }

                      if ( use_control_timestep_iterations_automatic ) {
                        if ( number_iterations==iteration_max ) {
                          if ( converged_once ) {
                            dtime /= 10.;
                            if ( dtime>EPS_DTIME*dtime_initial ) goto start_of_step;
                          }
                        }
                        else {
                          converged_once = 1;
                          dtime *= 1. + (ITERATION_WISHED-number_iterations)/20.;
                          if ( dtime>maximum_timestep ) dtime = maximum_timestep;
                        }
                      }
                      if ( use_control_timestep_size_automatic_decrease ) {
			  if ((post_node_rhside_ratio > ratio_max) && !(dtime<min_timestep))
			  	dtime /= decrease_factor;
                      }
                      if ( dtime>EPS_DTIME*dtime_initial ) {
				// added for options_element_dof

 		        db_max_index( ELEMENT, max_elem, VERSION_NORMAL, GET );
 		        for ( ielem=0; ielem<=max_elem; ielem++ ) {
			  if ( db_active_index( ELEMENT, ielem, VERSION_NORMAL )) {
			    if(options_element_dof==-YES) {	
     	                      if ( nuknwn>0 ) {
			    	  db( ELEMENT_DOF_INITIALISED, ielem, &one, ddum, one, VERSION_NORMAL, PUT );
				  db( ELEMENT_DOF_INITIALISED, ielem, &one, ddum, one, VERSION_NEW, PUT );
			      }
			    }	
  		            if (scalar_dabs(options_nonlocal_softvar)>TINY) {
			      array_set(dworknei, 0., length_nei);

			      db( NONLOCAL_ELEMENT_INFO, ielem, idum, dworknei, length_nei, 
			  	  VERSION_NORMAL, GET );
			      if(dworknei[length_nei-2] != dworknei[length_nei-1]) {
				  nonlocal_first_set=0;
				  dworknei[length_nei-1] = dworknei[length_nei-2];
			      }	
			      dworknei[length_nei-2]=0.;	//set to 1 if element used for
								    	//nonlocal calulation
			      db( NONLOCAL_ELEMENT_INFO, ielem, idum, dworknei, length_nei, 
			  	  VERSION_NORMAL, PUT );		
			    }
			  }
			}

                        db_version_copy( VERSION_NEW, VERSION_NORMAL );
                        db_version_delete( VERSION_NEW );
                        step_close( YES, ipar, npar, ipar_i, ipar_n );
                        if ( repeat(start_control) ) goto repeat_point;
                      }
                      else if ( time_current<time_new ) {
                        if ( control_timestep_iterations_automatic_stop==-YES ) {
                          pri( "\nWarning: the minimal time step size is reached." );
                          pri( "The calculation is terminated.\n" );
                          exit_tn( -YES );
                        }
                      }
                    }
                  }
                  swit = swit_timestep;
                }
                else {
                  step_start( NO, &options_solver, dtime, time_current );
                  step_close( NO, ipar, npar, ipar_i, ipar_n );
                  if ( repeat(start_control) ) goto repeat_point;
                }
                inverse_calculation( ipar, npar, ipar_i, ipar_n,
                  max_control, INVERSE_DETERMINE_SENSITIVITY );
              }
            }
          }
          if ( npar>0 && !( inverse_iter==inverse_iterations && 
              ipar==npar-1 && ipar_i==ipar_n-1 ) ) {
            db_version_copy( VERSION_START, VERSION_NORMAL );
            db_copy( NODE, NODE_START_REFINED, VERSION_NORMAL );
            db_copy( NODE_DOF, NODE_DOF_START_REFINED, VERSION_NORMAL );
            db_delete( ELEMENT_TENDON_NUMBER, VERSION_NORMAL );
            db_delete( ELEMENT_TENDON_VOLUME, VERSION_NORMAL );
            db_delete( ELEMENT_TENDON_STRAIN, VERSION_NORMAL );
            db_delete( ELEMENT_TENDON_STRESS, VERSION_NORMAL );
            db_delete( ELEMENT_TENDON_DIRECTION, VERSION_NORMAL );
            mesh_has_changed( VERSION_NORMAL );
          }
        }
      }
    }
    inverse_calculation( ipar, npar, ipar_i, ipar_n, 
      max_control, INVERSE_DETERMINE_NEW_ESTIMATES );
  }

  delete[] dworknei;
  if ( swit ) pri( "Out routine TOP" );

}

void step_start( long int task, long int options_solver[], double dtime, double time_current )

{
  long int icontrol=0, control_mesh_remesh=0, control_materi_diffusion=0,
    element=0, max_element=0, length=0, nnol=0, mnol=0,
    name=0, any_beam=0, any_truss=0, any_spring=0, any_contactspring=0,
    ldum=0, options_matrix_group=-NO, options_matrix_length=0, 
    element_group=0, max_group=0, exit_tochnog=0,
    el[1+MNOL], control_adjust_geometry[4];
  double ddum[1];

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );

  input_runtime();

  if ( db_active_index( EXIT_TOCHNOG, 0, VERSION_NORMAL )  ) {
    db( EXIT_TOCHNOG, 0, &exit_tochnog, ddum, ldum, VERSION_NORMAL, GET );
    if      ( exit_tochnog==-YES ) 
      exit_tn( -YES );
    else if ( exit_tochnog==-RESTART ) {
      exit_tn( -RESTART );
    }
  }

  macro();

  if ( db_active_index( CONTROL_MESH_ADJUST_GEOMETRY, icontrol, VERSION_NORMAL )  ) {
    db( CONTROL_MESH_ADJUST_GEOMETRY, icontrol, control_adjust_geometry, 
      ddum, ldum, VERSION_NORMAL, GET );
    adjust_geom( &control_adjust_geometry[0], &control_adjust_geometry[2] ); 
  }

  change_geometry( task, dtime, time_current );

  data( task, dtime, time_current ); 

  merge(); 

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  if ( db_active_index( CONTROL_MESH_REMESH, icontrol, VERSION_NORMAL ) ) {
    db( CONTROL_MESH_REMESH, icontrol, &control_mesh_remesh, ddum, 
      ldum, VERSION_NORMAL, GET );
    if ( control_mesh_remesh==-YES  ) {
      area_node_dataitem();
      remesh( VERSION_NORMAL );
    }
  }

  db_delete( NODE_REMESH_VELOCITY, VERSION_NORMAL );

  delete_geom( time_current ); 

  extrude();

  failure( time_current ); 

  distribute();

  if ( task==YES ) area_element_group_sequence( );

  crack();

  db( CONTROL_MATERI_DIFFUSION, icontrol, &control_materi_diffusion, ddum, 
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  if ( control_materi_diffusion==-INITIALIZE ) materi_diffusion_calculate( INITIALIZE );

  if ( task==YES ) {

    db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
    for ( element=0; element<=max_element; element++ ) {
      if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) ) {
        db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
        name = el[0]; nnol = length - 1;
        if ( nnol>mnol ) mnol = nnol;
        if ( name==-BEAM || name==-TRUSSBEAM ) any_beam = 1;
        if ( name==-TRUSS || name==-TRUSSBEAM ) any_truss = 1;
        if ( name==-SPRING1 || name==-SPRING2 ) any_spring = 1;
        if ( name==-CONTACTSPRING ) any_contactspring = 1;
        element_group = 0;
        db( ELEMENT_GROUP, element, &element_group, ddum, ldum,
          VERSION_NORMAL, GET_IF_EXISTS );
        if ( element_group>max_group ) max_group = element_group;
      }
    }               

    db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
    if ( options_solver[0]!=-NONE && options_solver[0]!=-DIAGONAL ) { 

      for ( element=0; element<=max_element; element++ ) {
        if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) ) {
          length = db_len( ELEMENT, element, VERSION_NORMAL );
          nnol = length - 1;
          if ( nnol>mnol ) mnol = nnol;
        }
      }
      length = mnol*nprinc*mnol*nprinc;
      if ( db( OPTIONS_MATRIX_LENGTH, 0, &options_matrix_length, ddum, ldum,
        VERSION_NORMAL, GET_IF_EXISTS ) ) length = options_matrix_length;
      db( OPTIONS_MATRIX_GROUP, 0, &options_matrix_group, ddum, ldum,
        VERSION_NORMAL, GET_IF_EXISTS );
      eigen_active = 0;
      if ( options_matrix_group==-YES ) {
        db_data_length_put( GROUP_MATRIX_VALUES, length );
        db_allocate( GROUP_MATRIX_VALUES, max_group, VERSION_NORMAL, MINIMAL );
        db_data_length_put( GROUP_MATRIX_UNKNOWNS, 4*length );
        db_allocate( GROUP_MATRIX_UNKNOWNS, max_group, VERSION_NORMAL, MINIMAL );
        if ( db_partialname_any("control_eigen") ) {
          eigen_active = 1;
          db_data_length_put( GROUP_MATRIX_SECOND_VALUES, length );
          db_allocate( GROUP_MATRIX_SECOND_VALUES, max_group, VERSION_NORMAL, MINIMAL );
        }                  
      }
      else {
        db_data_length_put( ELEMENT_MATRIX_VALUES, length );
        db_allocate( ELEMENT_MATRIX_VALUES, max_element, VERSION_NORMAL, MINIMAL );
        db_data_length_put( ELEMENT_MATRIX_UNKNOWNS, 2*length );
        db_allocate( ELEMENT_MATRIX_UNKNOWNS, max_element, VERSION_NORMAL, MINIMAL );
        if ( db_partialname_any("control_eigen") ) {
          eigen_active = 1;
          db_data_length_put( ELEMENT_MATRIX_SECOND_VALUES, length );
          db_allocate( ELEMENT_MATRIX_SECOND_VALUES, max_element, VERSION_NORMAL, MINIMAL );
        }
      }

    }

    if ( materi_velocity ) {
      db_allocate( ELEMENT_MASS, max_element, VERSION_NORMAL, MINIMAL );
    }

    if ( db_partialname_any("group_materi_failure") || 
         db_partialname_any("control_mesh_delete") ) {
      db_allocate( ELEMENT_MATRIX_DELETE, max_element, VERSION_NORMAL, MINIMAL );
      db_allocate( ELEMENT_RHSIDE_DELETE, max_element, VERSION_NORMAL, MINIMAL );
    }

    if ( materi_diffusion ) {
      db_allocate( ELEMENT_EMPTY, max_element, VERSION_NEW, MINIMAL );
    }
    if ( materi_density ) {
      db_allocate( ELEMENT_EMPTY, max_element, VERSION_NEW, MINIMAL );
    }
    if ( materi_strainenergy ) {
      db_allocate( ELEMENT_STRAINENERGY, max_element, VERSION_NORMAL, MINIMAL );
    }
    if ( any_beam ) {
      db_allocate( ELEMENT_BEAM_DIRECTION, max_element, VERSION_NEW, MINIMAL );
      db_allocate( ELEMENT_BEAM_MOMENT, max_element, VERSION_NEW, MINIMAL );
    }
    if ( any_contactspring ) {
      db_allocate( ELEMENT_CONTACTSPRING_DIRECTION, max_element, VERSION_NEW, MINIMAL );
      db_allocate( ELEMENT_CONTACTSPRING_FORCE, max_element, VERSION_NEW, MINIMAL );
    }
    if ( any_spring ) {
      db_allocate( ELEMENT_SPRING_DIRECTION, max_element, VERSION_NEW, MINIMAL );
      db_allocate( ELEMENT_SPRING_FORCE, max_element, VERSION_NEW, MINIMAL );
    }
    if ( any_truss ) {
      db_allocate( ELEMENT_TRUSS_DIRECTION, max_element, VERSION_NEW, MINIMAL );
      db_allocate( ELEMENT_TRUSS_FORCE, max_element, VERSION_NEW, MINIMAL );
    }

  }

}

void iteration_start( void )

  // allocate before a parallel element loop is entered
{
  db_set_dbl( NODE_LHSIDE, VERSION_NORMAL );
  db_set_dbl( NODE_RHSIDE, VERSION_NORMAL );
  db_set_dbl( ELEMENT_VOLUME, VERSION_NORMAL );
}

void step_close( long int task, long int ipar, long int npar, long int ipar_i, long int ipar_n )

{
  long int i=0, nval=0, data_item=0, icontrol=0, ldum=0, control_split=0,
    use_control_refine_globally_geometry=0, length_control_refine_globally=0,
    length=0, time_of_calculation=0, time_at_start=0, time_at_end=0,
    print_lastdatabase=-NO, control_refine_globally[4], control_refine_globally_geometry[2], 
    idum[1], renumber[2], *ival=NULL;
  double ddum[MDIM];

  ival = get_new_int(DATA_ITEM_SIZE);

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );

  if ( db_active_index( CONTROL_MESH_RENUMBER, icontrol, VERSION_NORMAL ) ) {
    db( CONTROL_MESH_RENUMBER, icontrol, renumber, ddum, ldum, VERSION_NORMAL, GET );
    if ( renumber[0]<0 || renumber[1]<0 ) db_error( CONTROL_MESH_RENUMBER, icontrol );
    renumbering( VERSION_NORMAL, NO, renumber[0], renumber[1], idum, idum );
  }

  if ( db_active_index( CONTROL_EIGEN, icontrol, VERSION_NORMAL ) ) solve( -CONTROL_EIGEN );

  if ( db_active_index( CONTROL_MESH_REFINE_GLOBALLY, icontrol, VERSION_NORMAL ) ) {
    error( PUT );
    db( CONTROL_MESH_REFINE_GLOBALLY, icontrol, control_refine_globally, 
      ddum, length_control_refine_globally, VERSION_NORMAL, GET );
    use_control_refine_globally_geometry = db( CONTROL_MESH_REFINE_GLOBALLY_GEOMETRY, 
      icontrol, control_refine_globally_geometry, ddum, ldum, 
      VERSION_NORMAL, GET_IF_EXISTS );
    refine_globally( control_refine_globally, 
      length_control_refine_globally, 
      use_control_refine_globally_geometry, 
      control_refine_globally_geometry, 
      PROJECT_EXACT, VERSION_NORMAL );
  }

  if ( db_active_index( CONTROL_MESH_SPLIT, icontrol, VERSION_NORMAL ) ) {
    db( CONTROL_MESH_SPLIT, icontrol, &control_split, ddum, ldum, VERSION_NORMAL, GET );
    if ( control_split!=-NO ) {
      mesh_split( VERSION_NORMAL ); 
      renumbering( VERSION_NORMAL, NO, 1, 1, idum, idum );
    }
  }

  new_mesh();

  mesh_delete_small( VERSION_NORMAL );

  generate_beam_truss( icontrol, BEAM ); 
  generate_beam_truss( icontrol, TRUSS ); 
  generate_beam_truss( icontrol, TRUSSBEAM ); 
  generate_spring( icontrol );

  if ( task==YES ) maxwell_scatter();

  refine_locally();

  restart();

  unknown_reset();

  if ( materi_diffusion && task==YES ) {
    materi_diffusion_calculate( CLOSE );
    materi_diffusion_fill();
    materi_diffusion_temperature();
  }

  post( task ); 

  if ( task==YES ) error( GET );

  if ( task==YES ) {
    db( PRINT_LASTDATABASE, 0, &print_lastdatabase, ddum, 
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
    if ( print_lastdatabase==-YES ) {
      print_database( -1, VERSION_NORMAL, -PRINT_LASTDATABASE );
    }
  }

  calculate();

  db( TIME_AT_START, 0, &time_at_start, ddum, ldum, VERSION_NORMAL, GET );
  time_at_end = (long int) time(NULL);
  time_of_calculation = time_at_end - time_at_start;
  length=1; db( TIME_CALCULATION, 0, &time_of_calculation, 
    ddum, length, VERSION_NORMAL, PUT );

  if ( npar==0 || ( ipar==npar-1 && ipar_i==ipar_n-1 ) ) {
    if ( db_active_index( CONTROL_PRINT, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT, icontrol, ival, ddum, nval, VERSION_NORMAL, GET );
      for ( i=0; i<nval; i++ ) {
        data_item = ival[i];
        if ( data_item>=0 ) db_error( CONTROL_PRINT, icontrol );
        print_database( icontrol, VERSION_NORMAL, data_item );
      }
      cout << "\n\n";
    }
    if ( db_active_index( CONTROL_PRINT_DATABASE, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_DATABASE, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      if ( ival[0]==-RESTART ) print_restart( icontrol );
      else print_database( icontrol, VERSION_NORMAL, ival[0] );
    }
    if ( db_active_index( CONTROL_PRINT_DATA_VERSUS_DATA, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_DATA_VERSUS_DATA, icontrol, ival, ddum, length, VERSION_NORMAL, GET );
      print_data_versus_data( ival, length );
    }
    if ( db_active_index( CONTROL_PRINT_DX, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_DX, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      if ( ival[0]==-YES ) print_dx( -NO );
      else if ( ival[0]!=-NO ) db_error( CONTROL_PRINT_DX, icontrol );
    }
    if ( db_active_index( CONTROL_PRINT_ELEMENT, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_ELEMENT, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      print_element( ival[0] );
    }
    if ( db_active_index( CONTROL_PRINT_HISTORY, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_HISTORY, icontrol, ival, ddum, nval, VERSION_NORMAL, GET );
      print_history( ival, nval );
    }
    if ( db_active_index( CONTROL_PRINT_PLOTMTV, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_PLOTMTV, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      if ( ival[0]==-YES || ival[0]>=0 ) print_plotmtv( icontrol, ival );
    }
    if ( db_active_index( CONTROL_PRINT_GID, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_GID, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      if ( ival[0]!=-NO ) print_gid( ival[0] );
    }
    if ( db_active_index( CONTROL_PRINT_GMV, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_GMV, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      if ( ival[0]==-YES || ival[0]>=0 ) print_gmv( icontrol, ival );
    }
    if ( db_active_index( CONTROL_PRINT_MATLAB, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_MATLAB, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      if ( ival[0]==-YES || ival[0]>=0 ) print_matlab( );
    }
    if ( db_active_index( CONTROL_PRINT_TECPLOT, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_TECPLOT, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      if ( ival[0]!=-NO  ) print_tecplot( ival );
    }
    if ( db_active_index( CONTROL_PRINT_UNKNOWNS, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_UNKNOWNS, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      if ( ival[0]==-YES ) print_unknowns( );
    }
    if ( db_active_index( CONTROL_PRINT_UNKNOWNSRHSIDE, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_UNKNOWNSRHSIDE, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      if ( ival[0]==-YES ) print_unknownsrhside( );
    }
    if ( db_active_index( CONTROL_PRINT_VTK, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_VTK, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      if ( ival[0]==-YES ) print_vtk( icontrol );
    }
  }
  cout << flush;

  delete[] ival;

}

