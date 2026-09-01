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
    timestep_iterations_automatic_apply=-YES, 
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
    multiplier=1., ddum[1], dzero=0., control_timestep_iterations_automatic[2], 
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

  // tochnog_version (manual Professional 6.1091): the build date as a
  // queryable record (day, month, year) - parsed from the compiler's
  // __DATE__ string so the record always matches the binary
  {
    long int version[3] = { 0, 0, 0 };
    char month_str[16] = "";
    long int day = 0, year = 0, imonth = 0;
    char *months[12] = { (char*)"Jan", (char*)"Feb", (char*)"Mar",
      (char*)"Apr", (char*)"May", (char*)"Jun", (char*)"Jul",
      (char*)"Aug", (char*)"Sep", (char*)"Oct", (char*)"Nov",
      (char*)"Dec" };
    if ( sscanf( __DATE__, "%15s %ld %ld", month_str, &day, &year )==3 ) {
      for ( imonth=0; imonth<12; imonth++ )
        if ( !strcmp( month_str, months[imonth] ) ) break;
      if ( imonth<12 ) {
        version[0] = day;
        version[1] = imonth + 1;
        version[2] = year;
      }
    }
    length = 3;
    db( TOCHNOG_VERSION, 0, version, ddum, length, VERSION_NORMAL, PUT );
  }

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

  db( CHECK_USED, 0, &check_used, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CHECK_TARGET, 0, &check_target, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CHECK_ERROR, 0, &check_error, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CHECK_WARNING, 0, &check_warning, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CHECK_NAN, 0, &check_nan, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CHECK_DATA, 0, &check_data, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CHECK_SOLVER, 0, idum, &check_solver_eps, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CHECK_ELEMENT_SHAPE, 0, idum, &check_element_shape_factor, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CHECK_MEMORY, 0, &check_memory, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CHECK_MEMORY_USAGE, 0, &check_memory_usage, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );

  // check_data -yes: verify data base integrity (required items present)
  if ( check_data==-YES ) check_data_integrity();

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
      // element_intpnt_materi_plasti_hardsoil_gammap_initial (manual
      // Professional 6.440): pre-allocate the per-element record before
      // the time loop (the first-timestep PUT in set_stress runs inside
      // the parallel element loop, where allocation is not allowed).
      // Sentinel -1: the first-timestep initialization in set_stress
      // overwrites it with gamma_p_extra >= 0.
      if ( db_active_index( CONTROL_MATERI_PLASTI_HARDSOIL_GAMMAP_INITIAL,
           0, VERSION_NORMAL ) && db_active_index( ELEMENT, ielem, VERSION_NORMAL ) && 
           !db_active_index( ELEMENT_INTPNT_MATERI_PLASTI_HARDSOIL_GAMMAP_INITIAL,
             ielem, VERSION_NORMAL ) ) {
        dzero = -1.;
        db( ELEMENT_INTPNT_MATERI_PLASTI_HARDSOIL_GAMMAP_INITIAL, ielem, idum,
          &dzero, one, VERSION_NORMAL, PUT );
        dzero = 0.;
      }
    }
  }

    // initialize time, if not specified
  if ( !db_active_index( TIME_CURRENT, 0, VERSION_NORMAL ) ) {
    length = 1;
    db( TIME_CURRENT, 0, idum, &time_current, length, VERSION_NORMAL, PUT );
  }

  // debug: dump node 1 stress dofs after the input parse
  if ( db_active_index( NODE, 1, VERSION_NORMAL ) ) {
    double *ndof1 = db_dbl( NODE_DOF, 1, VERSION_NORMAL );
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
                  // timestep_iterations_automatic_apply -no (manual
                  // Professional 6.1090): neglect every
                  // control_timestep_iterations_automatic record
                  timestep_iterations_automatic_apply = -YES;
                  db( TIMESTEP_ITERATIONS_AUTOMATIC_APPLY, 0,
                    &timestep_iterations_automatic_apply, ddum, ldum,
                    VERSION_NORMAL, GET_IF_EXISTS );
                  if ( timestep_iterations_automatic_apply!=-NO &&
                       db_active_index(  CONTROL_TIMESTEP_ITERATIONS_AUTOMATIC, 
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
                    // manual Professional 6.386: [1]=minimal_timestep,
                    // [2]=maximum_timestep
                    maximum_timestep = control_timestep_iterations_automatic[2];
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
                    // materi_displacement_relative: the reference point is
                    // re-synchronized (dis_rel reset to 0) when the timestep
                    // changes (manual 4.13).
                    if ( materi_displacement_relative ) {
                      double dtime_ref = -1.;
                      db( MATERI_DISPLACEMENT_RELATIVE_REF, 0, idum, &dtime_ref,
                        ldum, VERSION_NORMAL, GET_IF_EXISTS );
                      if ( dtime_ref<0. || dtime_ref!=dtime_initial ) {
                        db_max_index( NODE, max_node, VERSION_NORMAL, GET );
                        for ( inod=0; inod<=max_node; inod++ ) {
                          if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
                            double *ndof = db_dbl( NODE_DOF, inod, VERSION_NORMAL );
                            for ( long int idim=0; idim<ndim; idim++ )
                              ndof[dis_rel_indx+idim*nder] = 0.;
                          }
                        }
                        db( MATERI_DISPLACEMENT_RELATIVE_REF, 0, idum,
                          &dtime_initial, length, VERSION_NORMAL, PUT );
                      }
                    }
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
                      // solver (manual Professional 6.1047): the GLOBAL
                      // solver type - overwrites every control_solver
                      // (the opposite precedence of options_solver,
                      // which the per-control record overwrites)
                      db( SOLVER, 0, &options_solver, ddum, length,
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

                      {
                        // control_print_number_iterations (manual
                        // Professional 6.336): console monitor - when the
                        // switch is -yes the iteration counter is printed
                        // on stdout during the equilibrium iterations of
                        // the time step. NOT the data item
                        // -inverse_iteration_number (iterations of the
                        // finished step).
                        long int print_number_iterations=-NO;
                        db( CONTROL_PRINT_NUMBER_ITERATIONS, icontrol,
                          &print_number_iterations, ddum, ldum,
                          VERSION_NORMAL, GET_IF_EXISTS );
                        for ( iteration=1; ( iteration<=iteration_min ||
                            post_node_rhside_ratio>=ratio_criterium ) &&
                            iteration<=iteration_max; iteration++ ) {
                          if ( print_number_iterations==-YES )
                            cout << "control_print_number_iterations: time "
                                 << time_current << " iteration " << iteration
                                 << "\n";
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
                          groundflow_total_pressure_limit_apply();
                          if ( print_where==-YES ) pri( "Where: after solver." ); 
                          parallel_sys_routine( &parallel_new_dof_diagonal );
                          locate();
                          nonlocal_apply();
                          post_node_rhside_fixed_free();
                          db( POST_NODE_RHSIDE_RATIO, 0, idum, &post_node_rhside_ratio, 
                            ldum, VERSION_NORMAL, GET );
                        }
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
                        // post_element_force (manual Professional
                        // 6.927): the cross-section forces/moments
                        // from the element internal nodal forces of
                        // the converged state
                        post_element_force_calculate();
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
    any_interface=0,
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

  if ( db_active_index( CONTROL_MESH_SWITCH, icontrol, VERSION_NORMAL ) ) {
    long int control_switch[DATA_ITEM_SIZE], switch_length=0;
    db( CONTROL_MESH_SWITCH, icontrol, control_switch, ddum, 
      switch_length, VERSION_NORMAL, GET );
    mesh_switch( control_switch, switch_length );
  }

  if ( db_active_index( CONTROL_MESH_MOVE, icontrol, VERSION_NORMAL ) ) {
    long int move_length=0, idum_m=0;
    double control_move[DATA_ITEM_SIZE];
    db( CONTROL_MESH_MOVE, icontrol, &idum_m, control_move, move_length,
      VERSION_NORMAL, GET );
    mesh_move( control_move, move_length );
  }

  if ( db_active_index( CONTROL_MESH_MIRROR, icontrol, VERSION_NORMAL ) ) {
    long int mirror_axis=0;
    db( CONTROL_MESH_MIRROR, icontrol, &mirror_axis, ddum, ldum,
      VERSION_NORMAL, GET );
    mesh_mirror( mirror_axis );
  }

  if ( db_active_index( CONTROL_MESH_COPY, icontrol, VERSION_NORMAL ) ) {
    long int copy_length=0, idum_c=0;
    double control_copy[MDIM];
    db( CONTROL_MESH_COPY, icontrol, &idum_c, control_copy, copy_length,
      VERSION_NORMAL, GET );
    mesh_copy( control_copy );
  }

  // control_mesh_rotate_angle: rotate the 2D mesh around z-axis
  if ( db_active_index( CONTROL_MESH_ROTATE_ANGLE, icontrol, VERSION_NORMAL ) ) {
    long int rot_length=0, idum_r=0;
    double rot_angle=0.;
    db( CONTROL_MESH_ROTATE_ANGLE, icontrol, &idum_r, &rot_angle, rot_length,
      VERSION_NORMAL, GET );
    mesh_rotate_2d( rot_angle );
  }

  // control_mesh_rotate: rotate 2D mesh to 3D (tria3->prism6, quad4->hex8)
  if ( db_active_index( CONTROL_MESH_ROTATE, icontrol, VERSION_NORMAL ) ) {
    long int nrot=1;
    db( CONTROL_MESH_ROTATE, icontrol, &nrot, ddum, ldum,
      VERSION_NORMAL, GET );
    mesh_rotate_3d( nrot );
  }

  // delete/keep elements, keep nodes, change element group
  if ( db_active_index( CONTROL_MESH_DELETE_ELEMENT, icontrol, VERSION_NORMAL ) ||
       db_active_index( CONTROL_MESH_KEEP_ELEMENT, icontrol, VERSION_NORMAL ) ||
       db_active_index( CONTROL_MESH_KEEP_ELEMENT_GROUP, icontrol, VERSION_NORMAL ) ||
       db_active_index( CONTROL_MESH_KEEP_NODE, icontrol, VERSION_NORMAL ) ||
       db_active_index( CONTROL_MESH_CHANGE_ELEMENT_GROUP, icontrol, VERSION_NORMAL ) ) {
    mesh_delete_keep( icontrol );
  }

  // control_mesh_remove: remove elements by method
  if ( db_active_index( CONTROL_MESH_REMOVE, icontrol, VERSION_NORMAL ) ) {
    mesh_remove( icontrol );
  }

  // control_mesh_convert: convert interface elements (bar2->quad4 in 2D,
  // tria3->prism6 / quad4->hex8 in 3D). Must run BEFORE the first
  // element_loop of the step so the converted elements are assembled.
  // EXTRUDE FIRST, CONVERT SECOND (2026-08-31, interface_bar2_hex8):
  // the bar2 interface is extruded along z to a quad4 in the XZ plane,
  // and only then the convert lifts it to a hex8 interface whose sides
  // connect in Y (the extrusion direction of the solids). The opposite
  // order (convert then extrude) produced a hex8 interface connecting in
  // Z - the push in Y was never transmitted (sigyy=0).
  // Both only on the first step (task==YES): re-running them on later
  // steps duplicates the mesh. They must run BEFORE data()/merge() so
  // those phases see the final 3D elements.
  if ( task==YES ) extrude();
  if ( task==YES ) interface_convert( icontrol );

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

  failure( time_current ); 

  distribute();

  // area_element_group_time -yes: re-evaluate the area_element_group
  // records at all times, not only when the mesh changes (manual 6.6)
  if ( task==YES && area_element_group_time_active() )
    area_element_group( VERSION_NORMAL );
  if ( task==YES ) area_element_group_sequence( );

  crack();

  db( CONTROL_MATERI_DIFFUSION, icontrol, &control_materi_diffusion, ddum, 
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  if ( control_materi_diffusion==-INITIALIZE ) materi_diffusion_calculate( INITIALIZE );

  // generate interface elements BEFORE scanning any_interface so the
  // interface histories are allocated for the generated elements too.
  generate_interface( icontrol );

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
        if ( db_active_index( GROUP_INTERFACE, element_group, VERSION_NORMAL ) )
          any_interface = 1;
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
    if ( any_interface ) {
      db_allocate( ELEMENT_INTERFACE_STRAIN_NORMAL, max_element, VERSION_NEW, MINIMAL );
      db_allocate( ELEMENT_INTERFACE_FORCE_TANG, max_element, VERSION_NEW, MINIMAL );
      db_allocate( ELEMENT_INTERFACE_FORCE_TANG2, max_element, VERSION_NEW, MINIMAL );
      // output records (Professional compatibility): allocated BEFORE the
      // parallel element loop (db_allocate cannot run inside it), both
      // versions so print_database / target checker (VERSION_NORMAL) and
      // the next step's loop (VERSION_NEW) can read them.
      db_allocate( ELEMENT_INTERFACE_INTPNT_STRESS, max_element, VERSION_NORMAL, MINIMAL );
      db_allocate( ELEMENT_INTERFACE_INTPNT_STRESS, max_element, VERSION_NEW, MINIMAL );
      db_allocate( ELEMENT_INTERFACE_INTPNT_STRAIN, max_element, VERSION_NORMAL, MINIMAL );
      db_allocate( ELEMENT_INTERFACE_INTPNT_STRAIN, max_element, VERSION_NEW, MINIMAL );
      db_allocate( ELEMENT_INTERFACE_STRESS_AVERAGE, max_element, VERSION_NORMAL, MINIMAL );
      db_allocate( ELEMENT_INTERFACE_STRESS_AVERAGE, max_element, VERSION_NEW, MINIMAL );
      db_allocate( ELEMENT_INTERFACE_STRAIN_AVERAGE, max_element, VERSION_NORMAL, MINIMAL );
      db_allocate( ELEMENT_INTERFACE_STRAIN_AVERAGE, max_element, VERSION_NEW, MINIMAL );
      db_allocate( ELEMENT_INTERFACE_INTPNT_MATERI_TENSION_STATUS, max_element, VERSION_NORMAL, MINIMAL );
      db_allocate( ELEMENT_INTERFACE_INTPNT_MATERI_TENSION_STATUS, max_element, VERSION_NEW, MINIMAL );
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

// control_print_frequency gate (manual Professional 6.291/6.292).
// Returns 1 when the control_print_* records of this icontrol may print
// in the current step_close, 0 when the print is gated by a
// control_print_frequency_timeinterval / control_print_frequency_timestep
// record. The 3 exceptions (control_print, control_print_history and
// control_print_data_versus_data) are NOT gated by the caller.
//
// Semantics (manual): the gated prints run each time after a time
// interval has passed (timeinterval) or after N time steps (timestep),
// and ALWAYS also at the end of the time increment. The end of the
// increment is detected with the SAME condition the timestep loop in
// top() uses: TIME_NEW is set at the start of each control_timestep
// increment and the last step is clamped exactly to it, so
// time_current>=time_new means the increment of the control_timestep
// record finished (step_close is called once per step).
//
// State (last print time / steps since last print) lives in dedicated
// VERSION_NORMAL records per icontrol (pattern of control_print_gid_time)
// so it survives restarts; on restart the cadence is re-anchored at the
// start of the current increment. If both frequency records exist for the
// same icontrol the timeinterval record wins (documented decision).
long int control_print_frequency_allowed( long int icontrol )

{
  long int ldum=0, idum[1], frequency_timestep=0, frequency_count=0,
    zero=0, one=1;
  double ddum[1], frequency_interval=0., last_print_time=0.,
    time_current=0., time_new=0., time_old=0.;
  long int end_of_increment=0;

  // no frequency record for this icontrol -> normal behaviour (allowed)
  if ( !db_active_index( CONTROL_PRINT_FREQUENCY_TIMEINTERVAL, icontrol,
       VERSION_NORMAL ) &&
       !db_active_index( CONTROL_PRINT_FREQUENCY_TIMESTEP, icontrol,
       VERSION_NORMAL ) )
    return 1;

  db( TIME_CURRENT, 0, idum, &time_current, ldum, VERSION_NORMAL, GET );

  if ( db_active_index( CONTROL_TIMESTEP, icontrol, VERSION_NORMAL ) &&
       db( TIME_NEW, 0, idum, &time_new, ldum, VERSION_NORMAL,
       GET_IF_EXISTS ) && time_current>=time_new-1.e-9*scalar_dabs(time_new) )
    end_of_increment = 1;

  if ( db_active_index( CONTROL_PRINT_FREQUENCY_TIMEINTERVAL, icontrol,
       VERSION_NORMAL ) ) {
    db( CONTROL_PRINT_FREQUENCY_TIMEINTERVAL, icontrol, idum,
      &frequency_interval, ldum, VERSION_NORMAL, GET );
    if ( frequency_interval<=0. )
      db_error( CONTROL_PRINT_FREQUENCY_TIMEINTERVAL, icontrol );
    if ( db( CONTROL_PRINT_FREQUENCY_TIMEINTERVAL_TIME, icontrol, idum,
         &last_print_time, ldum, VERSION_NORMAL, GET_IF_EXISTS ) ) {
      // restart: the calculation was rewound, re-anchor the interval at
      // the start of the current increment
      if ( last_print_time>time_current ) {
        if ( db( TIME_OLD, 0, idum, &time_old, ldum, VERSION_NORMAL,
             GET_IF_EXISTS ) && time_old<=time_current )
          last_print_time = time_old;
        else
          last_print_time = time_current;
      }
    }
    else {
      // first use: anchor the interval at the START of the current
      // increment (TIME_OLD). The manual example (interval 0.15, dt 0.04,
      // increment 0.41) prints at 0.16, 0.32, 0.41 - anchoring at the
      // first step (0.04) would print at 0.20 instead of 0.16.
      if ( db( TIME_OLD, 0, idum, &time_old, ldum, VERSION_NORMAL,
           GET_IF_EXISTS ) && time_old<=time_current )
        last_print_time = time_old;
      else
        last_print_time = time_current;
    }
    if ( end_of_increment || time_current>=last_print_time+frequency_interval ) {
      // the print is done in this step_close: re-anchor the cadence
      db( CONTROL_PRINT_FREQUENCY_TIMEINTERVAL_TIME, icontrol, idum,
        &time_current, one, VERSION_NORMAL, PUT );
      return 1;
    }
    return 0;
  }
  else {
    // control_print_frequency_timestep: print after N time steps
    // (and always at the end of the time increment)
    db( CONTROL_PRINT_FREQUENCY_TIMESTEP, icontrol, &frequency_timestep,
      ddum, ldum, VERSION_NORMAL, GET );
    if ( frequency_timestep<=0 )
      db_error( CONTROL_PRINT_FREQUENCY_TIMESTEP, icontrol );
    db( CONTROL_PRINT_FREQUENCY_TIMESTEP_COUNT, icontrol, &frequency_count,
      ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    frequency_count++;
    if ( end_of_increment || frequency_count>=frequency_timestep ) {
      db( CONTROL_PRINT_FREQUENCY_TIMESTEP_COUNT, icontrol, &zero,
        ddum, one, VERSION_NORMAL, PUT );
      return 1;
    }
    db( CONTROL_PRINT_FREQUENCY_TIMESTEP_COUNT, icontrol, &frequency_count,
      ddum, one, VERSION_NORMAL, PUT );
    return 0;
  }
}

void step_close( long int task, long int ipar, long int npar, long int ipar_i, long int ipar_n )

{
  long int i=0, nval=0, data_item=0, icontrol=0, ldum=0, control_split=0,
    use_control_refine_globally_geometry=0, length_control_refine_globally=0,
    length=0,     time_of_calculation=0, time_at_start=0, time_at_end=0,
    print_lastdatabase=-NO, control_refine_globally[4], control_refine_globally_geometry[2], 
    idum[1], renumber[2], *ival=NULL, frequency_allowed=1;
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
    // control_print_frequency_timeinterval / control_print_frequency_timestep
    // (manual Professional 6.291/6.292): gate every control_print_* of this
    // icontrol EXCEPT control_print, control_print_history and
    // control_print_data_versus_data. The gate is per icontrol: when it
    // returns 0 none of the gated prints runs in this step_close.
    frequency_allowed = control_print_frequency_allowed( icontrol );
    // print_apply -no (manual Professional 6.968): ALL control_print_*
    // records are neglected (a global gate for the whole run)
    {
      long int print_apply = -YES;
      db( PRINT_APPLY, 0, &print_apply, ddum, ldum,
        VERSION_NORMAL, GET_IF_EXISTS );
      if ( print_apply==-NO ) frequency_allowed = 0;
    }
    if ( db_active_index( CONTROL_PRINT, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT, icontrol, ival, ddum, nval, VERSION_NORMAL, GET );
      for ( i=0; i<nval; i++ ) {
        data_item = ival[i];
        if ( data_item>=0 ) db_error( CONTROL_PRINT, icontrol );
        print_database( icontrol, VERSION_NORMAL, data_item );
      }
      cout << "\n\n";
    }
    if ( frequency_allowed && db_active_index( CONTROL_PRINT_DATABASE, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_DATABASE, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      if ( ival[0]==-RESTART ) print_restart( icontrol );
      else print_database( icontrol, VERSION_NORMAL, ival[0] );
    }
    if ( frequency_allowed && db_active_index( CONTROL_PRINT_DATABASE_METHOD, icontrol, VERSION_NORMAL ) ) {
      // control_print_database_method (manual Professional 6.266): method
      // -all (default) prints all base records, -size_tot prints the size
      // of all base records (plus the system matrix), -size_tot_large only
      // the records larger than 1 Mb (plus the system matrix).
      db( CONTROL_PRINT_DATABASE_METHOD, icontrol, ival, ddum, ldum,
        VERSION_NORMAL, GET );
      if      ( ival[0]==-ALL ) print_database( icontrol, VERSION_NORMAL, -EVERYTHING );
      else if ( ival[0]==-SIZETOT ) print_database( icontrol, VERSION_NORMAL, -SIZETOT );
      else if ( ival[0]==-SIZE_TOT_LARGE ) print_database( icontrol, VERSION_NORMAL, -SIZE_TOT_LARGE );
      else db_error( CONTROL_PRINT_DATABASE_METHOD, icontrol );
    }
    if ( frequency_allowed && db_active_index( CONTROL_PRINT_PARTIALNAME, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_PARTIALNAME, icontrol, ival, ddum, nval, VERSION_NORMAL, GET );
      print_partialname( icontrol, VERSION_NORMAL, ival, nval );
    }
    if ( db_active_index( CONTROL_PRINT_DATA_VERSUS_DATA, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_DATA_VERSUS_DATA, icontrol, ival, ddum, length, VERSION_NORMAL, GET );
      print_data_versus_data( ival, length );
    }
    if ( frequency_allowed && db_active_index( CONTROL_PRINT_DX, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_DX, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      if ( ival[0]==-YES ) print_dx( -NO );
      else if ( ival[0]!=-NO ) db_error( CONTROL_PRINT_DX, icontrol );
    }
    if ( frequency_allowed && db_active_index( CONTROL_PRINT_ELEMENT, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_ELEMENT, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      print_element( ival[0] );
    }
    if ( db_active_index( CONTROL_PRINT_HISTORY, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_HISTORY, icontrol, ival, ddum, nval, VERSION_NORMAL, GET );
      print_history( ival, nval );
      if ( db_active_index( CONTROL_PRINT_HISTORY_SMOOTH, icontrol, VERSION_NORMAL ) )
        print_history_smooth( ival, nval );
    }
    if ( frequency_allowed && db_active_index( CONTROL_PRINT_PLOTMTV, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_PLOTMTV, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      if ( ival[0]==-YES || ival[0]>=0 ) print_plotmtv( icontrol, ival );
    }
    if ( frequency_allowed && db_active_index( CONTROL_PRINT_GID, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_GID, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      if ( ival[0]!=-NO ) print_gid( ival[0] );
    }
    if ( frequency_allowed && db_active_index( CONTROL_PRINT_GMV, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_GMV, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      if ( ival[0]==-YES || ival[0]>=0 ) print_gmv( icontrol, ival );
    }
    if ( frequency_allowed && db_active_index( CONTROL_PRINT_MATLAB, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_MATLAB, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      if ( ival[0]==-YES || ival[0]>=0 ) print_matlab( );
    }
    if ( frequency_allowed && db_active_index( CONTROL_PRINT_TECPLOT, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_TECPLOT, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      if ( ival[0]!=-NO  ) print_tecplot( ival );
    }
    if ( frequency_allowed && db_active_index( CONTROL_PRINT_UNKNOWNS, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_UNKNOWNS, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      if ( ival[0]==-YES ) print_unknowns( );
    }
    if ( frequency_allowed && db_active_index( CONTROL_PRINT_UNKNOWNSRHSIDE, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_UNKNOWNSRHSIDE, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      if ( ival[0]==-YES ) print_unknownsrhside( );
    }
    if ( frequency_allowed && db_active_index( CONTROL_PRINT_VTK, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_VTK, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      if ( ival[0]==-YES ) print_vtk( icontrol );
    }
    if ( frequency_allowed && db_active_index( CONTROL_PRINT_TABULAR, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_TABULAR, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      if ( ival[0]==-YES ) print_tabular( icontrol );
    }
    if ( frequency_allowed && db_active_index( CONTROL_PRINT_GMSH, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_GMSH, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      if ( ival[0]!=-NO ) print_gmsh( icontrol, ival[0] );
    }
    if ( frequency_allowed && db_active_index( CONTROL_PRINT_FRD, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_FRD, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      if ( ival[0]!=-NO ) print_frd( icontrol, ival[0] );
    }
    if ( frequency_allowed && db_active_index( CONTROL_PRINT_INTERFACE_STRESS, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_INTERFACE_STRESS, icontrol, ival, ddum, ldum,
        VERSION_NORMAL, GET );
      if ( ival[0]!=-NO ) print_interface_stress( icontrol, ival[0] );
    }
    if ( frequency_allowed && db_active_index( CONTROL_PRINT_DOF, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_DOF, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      if ( ival[0]!=-NO ) print_dof( icontrol, ival[0] );
    }
    if ( frequency_allowed && db_active_index( CONTROL_PRINT_DOF_LINE, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_DOF_LINE, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      if ( ival[0]!=-NO ) print_dof_line( icontrol, ival[0] );
    }
    if ( frequency_allowed && db_active_index( CONTROL_PRINT_DOF_POINT, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_DOF_POINT, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      if ( ival[0]!=-NO ) print_dof_point( icontrol, ival[0] );
    }
    if ( frequency_allowed && db_active_index( CONTROL_PRINT_NODE, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_NODE, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
      print_node( icontrol, ival, ldum );
    }
    if ( frequency_allowed && db_active_index( CONTROL_PRINT_BEAM_FORCE_MOMENT, icontrol, VERSION_NORMAL ) ) {
      db( CONTROL_PRINT_BEAM_FORCE_MOMENT, icontrol, ival, ddum, ldum,
        VERSION_NORMAL, GET );
      print_beam_force_moment( icontrol, ival[0] );
    }
    if ( frequency_allowed && db_active_index( CONTROL_PRINT_MATERI_STRESS_FORCE, icontrol, VERSION_NORMAL ) ) {
      // control_print_materi_stress_force (manual Professional 6.328):
      // prints the post_calcul -materi_stress -force results to
      // materi_stress_force.<icontrol> (the manual "index" is the
      // record index, like control_print_beam_force_moment); the record
      // value is the method (-all / -primary).
      db( CONTROL_PRINT_MATERI_STRESS_FORCE, icontrol, ival, ddum, ldum,
        VERSION_NORMAL, GET );
      print_materi_stress_force( icontrol, ival[0] );
    }
  }
  cout << flush;

  // remember the reaction forces of this step for bounda_time_until_force
  if ( db_active_index( NODE_RHSIDE, 0, VERSION_NORMAL ) ||
       db_max_index( NODE_RHSIDE, ldum, VERSION_NORMAL, GET )>=0 ) {
    long int max_node_prev=0, inod_prev=0, idum_p[1];
    double *rhs_p=NULL;
    db_max_index( NODE, max_node_prev, VERSION_NORMAL, GET );
    for ( inod_prev=0; inod_prev<=max_node_prev; inod_prev++ ) {
      if ( db_active_index( NODE_RHSIDE, inod_prev, VERSION_NORMAL ) ) {
        rhs_p = db_dbl( NODE_RHSIDE, inod_prev, VERSION_NORMAL );
        db( NODE_RHSIDE_PREVIOUS, inod_prev, idum_p, rhs_p, npuknwn,
          VERSION_NORMAL, PUT );
      }
    }
  }

  // check node_dof for NAN values (check_nan)
  if ( check_nan==-YES ) check_nan_results( check_nan );

  delete[] ival;

}

