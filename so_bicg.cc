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

double *d1, *d2, *p, *Ad1_thread, *Ad2_thread, *p_thread, *residue_thread;
#define EPS_dAd 1.e-16
#define EPS_TMP 1.e-16
#define EPS_P 1.e-10

long int solve_iterative_bicg( void )

{
  long int icontrol=0, iter=0, max_iter=0, ilocal=0, nthread=0, ready=0, 
    max_node=0, succesful=1, print_solver=-NO, length=0, ldum=0, 
    swit=0, idum[1];
  double error=0., alpha=0., beta=0., dAd=0., last_error=0., check_error=0.,
   bicg_error_minimum=1.e-12, size_r1r2=0., initial_error=0., 
   bicg_error=1.e-10, ddum[1], 
   *r1=NULL, *r2=NULL, *Ad1=NULL, *Ad2=NULL, *p_tmp=NULL, *residue=NULL;

  swit = set_swit(-1,-1,"solve_iterative_bicg");
  if ( swit ) pri( "In routine SOLVE_ITERATIVE_BICG" );

  db_max_index( NODE, max_node, VERSION_NORMAL, GET );
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  db( OPTIONS_PROCESSORS, 0, &nthread, ddum, ldum, VERSION_NORMAL, GET );
  db( PRINT_SOLVER, 0, &print_solver, ddum, ldum, 
    VERSION_NORMAL, GET_IF_EXISTS );

    // r1, r2 vectors
  r1 = get_new_dbl( solve_nlocal );
  r2 = get_new_dbl( solve_nlocal );
    // residue vector
  residue = get_new_dbl ( solve_nlocal );
    // search directions
  d1 = get_new_dbl( solve_nlocal );
  d2 = get_new_dbl( solve_nlocal );
  array_set( d1, 0., solve_nlocal );
  array_set( d2, 0., solve_nlocal );
    // preconditioner
  p = get_new_dbl( solve_nlocal );
  array_set( p, 1., solve_nlocal );
    // work array for inverse of preconditioner
  p_tmp = get_new_dbl( solve_nlocal );
    // A*d
  Ad1 = get_new_dbl( solve_nlocal );
  Ad2 = get_new_dbl( solve_nlocal );
    // A*d vectors for all threads
  Ad1_thread = get_new_dbl( nthread*solve_nlocal );
  Ad2_thread = get_new_dbl( nthread*solve_nlocal );
    // preconditioner for all threads
  p_thread = get_new_dbl( nthread*solve_nlocal );
    // residue vector for all threads
  residue_thread = get_new_dbl ( nthread*solve_nlocal );

    // set diagonal preconditioner
  solve_iterative_bicg_sys( Ad1, Ad2, p_tmp, residue );
  for ( ilocal=0; ilocal<solve_nlocal; ilocal++ ) 
    p[ilocal] = 1./(sqrt(scalar_dabs(p_tmp[ilocal])));

    // initial error
  solve_iterative_bicg_sys( Ad1, Ad2, p_tmp, residue );
  initial_error = array_inproduct( residue, residue, solve_nlocal );

    // relative check error
  db( CONTROL_OPTIONS_SOLVER_BICG_ERROR, icontrol, idum, &bicg_error,
	  ldum, VERSION_NORMAL, GET_IF_EXISTS );
  check_error = bicg_error * initial_error;
  if ( swit ) pri( "check_error", check_error );

    // minimum check error
  db( CONTROL_OPTIONS_SOLVER_BICG_ERROR_MINIMUM, icontrol, idum, &bicg_error_minimum,
	  ldum, VERSION_NORMAL, GET_IF_EXISTS );
  if ( check_error<bicg_error_minimum ) check_error = bicg_error_minimum;

    // start values for iterative loop
  for ( ilocal=0; ilocal<solve_nlocal; ilocal++ ) {
    r1[ilocal] = solve_b[ilocal] * p[ilocal];
    r2[ilocal] = solve_b[ilocal] * p[ilocal];
    d1[ilocal] = r1[ilocal];
    d2[ilocal] = r2[ilocal];
  }

    // iterative loop, preconditioned biconjugate gradients
  ready=0; max_iter = 10*solve_nlocal; error = 1.e10;
  for ( iter=0; !ready; iter++ ) {
    if ( swit ) pri( "iterative solve iteration", iter );
    solve_iterative_bicg_sys( Ad1, Ad2, p_tmp, residue );
    dAd = array_inproduct( d2, Ad1, solve_nlocal );
    size_r1r2 = array_inproduct( r1, r2, solve_nlocal );
    if ( swit ) {
      pri( "r1", r1, solve_nlocal );
      pri( "r2", r2, solve_nlocal );
      pri( "dAd", dAd );
      pri( "size_r1r2", size_r1r2 );
    }
    last_error = error;
    error = array_inproduct( residue, residue, solve_nlocal );
    if ( print_solver==-YES ) {
      cout << iter << " " << error << "\n";
      cout << flush;
    }
    if ( scalar_dabs(size_r1r2)<EPS_TMP || scalar_dabs(dAd)<EPS_dAd )
      ready = 1;
    else {
      alpha =  array_inproduct( r1, r2, solve_nlocal ) / dAd;
      for ( ilocal=0; ilocal<solve_nlocal; ilocal++ ) {
        solve_x[ilocal] += alpha * d1[ilocal];
        r1[ilocal] -= alpha * Ad1[ilocal];
        r2[ilocal] -= alpha * Ad2[ilocal];
      }
      beta = array_inproduct( r1, r2, solve_nlocal ) / size_r1r2;
      for ( ilocal=0; ilocal<solve_nlocal; ilocal++ ) {
        d1[ilocal] = r1[ilocal] + beta * d1[ilocal];
        d2[ilocal] = r2[ilocal] + beta * d2[ilocal];
      }
      if ( swit ) {
        pri( "alpha", alpha );
        pri( "solve_x", solve_x, solve_nlocal );
        pri( "error", error );
        pri( "beta", beta );
        pri( "d1", d1, solve_nlocal );
        pri( "d2", d2, solve_nlocal );
      }
    }
    if ( error<check_error )
      ready = 1;
    else if ( scalar_dabs(error-last_error)<1.e-1*check_error )
      ready = 1;
    else if ( iter==max_iter ) {
      pri( "" );
      pri( "Error: the, default, BI-CG solver did not converge." );
      pri( "The initial error in the linear equations is", initial_error );
      pri( "The final error in the linear equations is", error );
      pri( "That is not good enough.\n\n" );
      pri( "Probably you have a wrong input file." );
      pri( "- Not enough boundary conditions?" );
      pri( "- A diverged calculation?\n" );
      exit_tn_on_error();
    }
  }
  length = 1;
  db( GLOBAL_SOLVER_ITERATIONS, 0, &iter, ddum, length, VERSION_NORMAL, PUT );
  db( GLOBAL_SOLVER_ERROR, 0, idum, &error, length, VERSION_NORMAL, PUT );
  if ( print_solver==-YES ) {
    pri( "initial error", initial_error );
    pri( "check error", check_error );
    pri( "final error", error );
  }

    // fill solution vector
  for ( ilocal=0; ilocal<solve_nlocal; ilocal++ ) solve_b[ilocal] = solve_x[ilocal] * p[ilocal];
  if ( swit ) pri( "solution vector", solve_b, solve_nlocal );

  delete[] r1;
  delete[] r2;
  delete[] d1;
  delete[] d2;
  delete[] p;
  delete[] p_tmp;
  delete[] Ad1;
  delete[] Ad2;
  delete[] Ad1_thread;
  delete[] Ad2_thread;
  delete[] p_thread;
  delete[] residue;
  delete[] residue_thread;

  if ( swit ) pri( "Out routine SOLVE_ITERATIVE_BICG" );
  return succesful;
}

void solve_iterative_bicg_sys( double *Ad1, double *Ad2, double *p_tmp, double *residue )

{
  long int inod=0, max_node=0, ilocal=0, iglobal=0, ipuknwn=0,
    ithread=0, nthread=0, swit=0, ldum=0;
  double ddum[1], *node_lhside=NULL;

  swit = set_swit(-1,-1,"solve_iterative_bicg_sys");
  if ( swit ) pri( "In routine SOLVE_ITERATIVE_BICG" );

  db( OPTIONS_PROCESSORS, 0, &nthread, ddum, ldum, VERSION_NORMAL, GET );
  array_set( Ad1_thread, 0., nthread*solve_nlocal );
  array_set( Ad2_thread, 0., nthread*solve_nlocal );
  array_set( p_thread, 0., nthread*solve_nlocal );
  array_set( residue_thread, 0., nthread*solve_nlocal );

  parallel_sys_routine( &parallel_solve_iterative_bicg_element );

  array_set( Ad1, 0., solve_nlocal );
  array_set( Ad2, 0., solve_nlocal );
  array_set( p_tmp, 0., solve_nlocal );
  array_set( residue, 0., solve_nlocal );
  db_highest_index( NODE_LHSIDE, max_node, VERSION_NORMAL );
  for ( inod=0; inod<=max_node; inod++ ) {
    if ( db_active_index(  NODE_LHSIDE, inod, VERSION_NORMAL ) ) {
      node_lhside = db_dbl( NODE_LHSIDE, inod, VERSION_NORMAL );
      for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
        iglobal = inod*npuknwn + ipuknwn;
        ilocal = solve_global_local[iglobal];
        if ( ilocal!=-NO ) {
          Ad1[ilocal] += node_lhside[ipuknwn] * d1[ilocal] * p[ilocal] * p[ilocal];
          Ad2[ilocal] += node_lhside[ipuknwn] * d2[ilocal] * p[ilocal] * p[ilocal];
          residue[ilocal] += solve_b[ilocal] * p[ilocal] - node_lhside[ipuknwn] * 
            solve_x[ilocal] * p[ilocal] * p[ilocal];
          p_tmp[ilocal] += node_lhside[ipuknwn];
        }
      }
    }
  }
  for ( ilocal=0; ilocal<solve_nlocal; ilocal++ ) {
    for ( ithread=0; ithread<nthread; ithread++ ) {
      Ad1[ilocal] += Ad1_thread[ithread*solve_nlocal+ilocal];
      Ad2[ilocal] += Ad2_thread[ithread*solve_nlocal+ilocal];
      p_tmp[ilocal] += p_thread[ithread*solve_nlocal+ilocal];
      residue[ilocal] += residue_thread[ithread*solve_nlocal+ilocal];
    }
    if ( scalar_dabs(p_tmp[ilocal])<EPS_P ) p_tmp[ilocal] = 1.;
  }
  if ( swit ) {
    pri( "Ad1", Ad1, solve_nlocal );
    pri( "Ad2", Ad2, solve_nlocal );
    pri( "p_tmp", p_tmp, solve_nlocal );
    pri( "p_thread", p_thread, nthread, solve_nlocal );
    pri( "residue_thread", residue_thread, nthread, solve_nlocal );
  }

  if ( swit ) pri( "Out routine SOLVE_ITERATIVE_BICG_SYS" );
}

void parallel_solve_iterative_bicg_element( void )

{
  long int element=0, iloop=0, nloop=0, ithread=0, max_element=0,
    *next_of_loop=NULL;

    // loop over elements
  db_highest_index( ELEMENT, max_element, VERSION_NORMAL );
  if ( max_element>=0 ) {
    next_of_loop = get_new_int(1+max_element);
    parallel_sys_next_of_loop( next_of_loop, max_element, nloop, ithread );
    for ( iloop=0; iloop<nloop; iloop++ ) {
      element = next_of_loop[iloop];
      if ( element>max_element )
        break;
      else
        solve_iterative_bicg_element( element, ithread );
    }
    delete[] next_of_loop;
  }

}

void solve_iterative_bicg_element( long int element, long int ithread )

{
  register long int i=0, iglobal=0, ilocal=0, jglobal=0, jlocal=0,
     indx1=0, indx2=0, element_group=0, ldum=0, element_matrix_values_length=0, 
     use_element_matrix=0, use_group_matrix=0;
  long int *element_matrix_unknowns=NULL;
  double ddum[1], *element_matrix_values=NULL;

  if ( db_active_index( ELEMENT_MATRIX_VALUES, element, VERSION_NORMAL ) )
    use_element_matrix = 1;
  else {
    if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) ) {
      db( ELEMENT_GROUP, element, &element_group, ddum, ldum, 
        VERSION_NORMAL, GET_IF_EXISTS );
      if ( db_active_index( GROUP_MATRIX_VALUES, element_group, VERSION_NORMAL ) )
        use_group_matrix = 1;
    }
  }
  if ( use_element_matrix || use_group_matrix ) {
    if ( use_element_matrix ) {
      element_matrix_unknowns = 
        db_int( ELEMENT_MATRIX_UNKNOWNS, element, VERSION_NORMAL );
      element_matrix_values = 
        db_dbl( ELEMENT_MATRIX_VALUES, element, VERSION_NORMAL );
      element_matrix_values_length =
        db_len( ELEMENT_MATRIX_VALUES, element, VERSION_NORMAL );
    }
    else {
      assert( use_group_matrix );
      element_matrix_unknowns = get_new_int( 2*MNOL*nprinc*MNOL*nprinc );
      get_element_matrix_unknowns( element, element_matrix_unknowns );
      element_matrix_values = 
        db_dbl( GROUP_MATRIX_VALUES, element_group, VERSION_NORMAL );
      element_matrix_values_length =
        db_len( GROUP_MATRIX_VALUES, element_group, VERSION_NORMAL );
    }
    for ( i=0; i<element_matrix_values_length; i++ ) {
      iglobal = element_matrix_unknowns[i*2+0];
      jglobal = element_matrix_unknowns[i*2+1];
      ilocal = solve_global_local[iglobal];
      jlocal = solve_global_local[jglobal];
      if ( ilocal!=-NO && jlocal!=-NO ) {
        indx1 = ithread*solve_nlocal + ilocal;
        indx2 = ithread*solve_nlocal + jlocal;
        Ad1_thread[indx1] += element_matrix_values[i] * d1[jlocal] * 
          ( p[ilocal] * p[jlocal] ) ;  
        Ad2_thread[indx2] += element_matrix_values[i] * d2[ilocal] * 
          ( p[ilocal] * p[jlocal] ) ;  
        residue_thread[indx2] -= element_matrix_values[i] * solve_x[ilocal] * 
          ( p[ilocal] * p[jlocal] ) ;  
        if ( ilocal==jlocal ) p_thread[indx1] += element_matrix_values[i];
      }
    }
    if ( use_group_matrix ) delete[] element_matrix_unknowns;
  }

}

