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
long int solve_iterative_bicg_use_cg = 0;
#define EPS_dAd 1.e-16
#define EPS_TMP 1.e-16
#define EPS_P 1.e-10
#define EPS_SYMMETRY 1.e-8
#define MAX_PUKNWN_SYMMETRY 32

long int solve_iterative_bicg_symmetric( void );
double solve_iterative_bicg_real_residual( double *residue, double norm_b );

long int solve_iterative_bicg( void )

{
  long int icontrol=0, iter=0, max_iter=0, ilocal=0, nthread=0, ready=0, 
    max_node=0, succesful=1, print_solver=-NO, length=0, ldum=0, 
    swit=0, idum[1], bicg_stop=-YES, icontrol_bicg=0, ldum_bs=0;
  double error=0., alpha=0., beta=0., dAd=0., check_error=0.,
   bicg_error_minimum=1.e-12, size_r1r2=0., initial_error=0., 
   bicg_error=1.e-10, ddum[1], ddum_bs[1], norm_b=0., real_residual=0.,
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

    // norm of the REAL right-hand side. solve_b holds the real RHS on
    // entry; it is overwritten with the solution at the end of this
    // routine, so the norm is saved here for the honest relative-residual
    // reporting (|b-Ax|/|b|).
  norm_b = sqrt( array_inproduct( solve_b, solve_b, solve_nlocal ) );

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
  // solver_bicg_error (manual Professional 6.1048): the global record
  // overwrites the per-control one
  db( SOLVER_BICG_ERROR, 0, idum, &bicg_error,
	  ldum, VERSION_NORMAL, GET_IF_EXISTS );
  check_error = bicg_error * initial_error;
  if ( swit ) pri( "check_error", check_error );

    // minimum check error
  db( CONTROL_OPTIONS_SOLVER_BICG_ERROR_MINIMUM, icontrol, idum, &bicg_error_minimum,
	  ldum, VERSION_NORMAL, GET_IF_EXISTS );
  if ( check_error<bicg_error_minimum ) check_error = bicg_error_minimum;

    // control_solver_bicg_stop -no (manual Professional 6.376):
    // do NOT stop the calculation when the iterative solver does not
    // converge (continue with the current solution); -yes (and the
    // default) stops as before
  if ( db_active_index( CONTROL_SOLVER_BICG_STOP, 0,
       VERSION_NORMAL ) ) {
    db( ICONTROL, 0, &icontrol_bicg, ddum_bs, ldum_bs,
      VERSION_NORMAL, GET_IF_EXISTS );
    db( CONTROL_SOLVER_BICG_STOP, icontrol_bicg, &bicg_stop,
      ddum_bs, ldum_bs, VERSION_NORMAL, GET_IF_EXISTS );
  }
  // solver_bicg_stop (manual Professional 6.1050): the global record
  // overwrites the per-control one
  db( SOLVER_BICG_STOP, 0, &bicg_stop,
    ddum_bs, ldum_bs, VERSION_NORMAL, GET_IF_EXISTS );

    // the velocity/temperature/pressure system assembled by the GNU is
    // symmetric when every element matrix is symmetric (the diagonal
    // NODE_LHSIDE block is diagonal). Measured per solve; symmetric ->
    // plain CG, non-symmetric (e.g. some coupled terms) -> honest Bi-CG.
  {
    long int matrix_symmetric = -NO;
    long int measured_symmetric = solve_iterative_bicg_symmetric( );
    // solver_matrix_symmetric -yes (manual Professional 6.1052): the
    // user asserts the system is symmetric - the per-solve measurement
    // is bypassed and the symmetric solver (CG) is used. The
    // Professional symmetrizes the matrix if needed; the GNU runs CG
    // on the as-assembled matrix, so a warning is printed when the
    // measurement disagrees with the user's assertion.
    db( SOLVER_MATRIX_SYMMETRIC, 0, &matrix_symmetric, ddum_bs, ldum_bs,
      VERSION_NORMAL, GET_IF_EXISTS );
    if ( matrix_symmetric==-YES ) {
      // solver_matrix_symmetric -yes (manual Professional 6.1053): the
      // matrices are symmetrized IF NEEDED so that a symmetric equation
      // solver can be used. The GNU does not re-symmetrize the assembled
      // matrix: when the measurement disagrees with the user's assertion
      // the honest Bi-CG runs on the as-assembled system (it converges
      // where plain CG on a non-symmetric system diverges).
      solve_iterative_bicg_use_cg = measured_symmetric;
      if ( swit ) pri( "solver_matrix_symmetric -yes: using CG" );
      if ( !measured_symmetric ) {
        pri( "Warning: solver_matrix_symmetric -yes but the measured "
             "system is NOT symmetric - running Bi-CG on the "
             "as-assembled matrix instead of CG (the Professional "
             "symmetrizes the matrix in this case)" );
      }
    }
    else {
      solve_iterative_bicg_use_cg = measured_symmetric;
    }
  }
  if ( swit ) {
    if ( solve_iterative_bicg_use_cg ) pri( "system is symmetric: using CG" );
    else pri( "system is NOT symmetric: using Bi-CG" );
  }

    // start values for iterative loop
  for ( ilocal=0; ilocal<solve_nlocal; ilocal++ ) {
    r1[ilocal] = solve_b[ilocal] * p[ilocal];
    r2[ilocal] = solve_b[ilocal] * p[ilocal];
    d1[ilocal] = r1[ilocal];
    d2[ilocal] = r2[ilocal];
  }

    // iterative loop, preconditioned conjugate gradients (symmetric) or
    // preconditioned biconjugate gradients (non-symmetric).
    //
    // Honest stopping criteria (DIAG-SOLVE-MIXTO fix A+B):
    // - success is reported ONLY when the residual of the current iterate
    //   is below the tolerance: error = |b - A x|^2 (preconditioned),
    //   recomputed from the iterate on every pass, compared against
    //   check_error = max(bicg_error*|r0|^2, bicg_error_minimum). This is
    //   a RELATIVE residual test (|r|/|r0| < sqrt(bicg_error), with an
    //   absolute floor) - it is not the "error decrease" heuristic.
    // - the old breakdown exits (dAd ~ 0, r1.r2 ~ 0) and the old
    //   stagnation exit (|error-last_error| < 0.1*check_error) reported
    //   SUCCESS without convergence (x = 0 or a partial solution with
    //   error >> check_error). They are removed: a breakdown is a real
    //   singularity, a stagnation just keeps iterating, and a solve that
    //   does not converge within max_iter fails honestly (RC != 0).
  ready=0; max_iter = 10*solve_nlocal; error = 1.e10;
  for ( iter=0; !ready; iter++ ) {
    if ( swit ) pri( "iterative solve iteration", iter );
    solve_iterative_bicg_sys( Ad1, Ad2, p_tmp, residue );
    error = array_inproduct( residue, residue, solve_nlocal );
    if ( print_solver==-YES ) {
      cout << iter << " " << error << "\n";
      cout << flush;
    }
    if ( error<check_error || error<EPS_TMP ) {
        // genuine convergence: residual below the (relative) tolerance,
        // or at machine-zero level
      ready = 1;
      break;
    }
    if ( iter==max_iter ) {
        // honest failure: no convergence within the iteration limit
      real_residual = solve_iterative_bicg_real_residual( residue, norm_b );
      if ( bicg_stop==-NO ) {
        pri( "Warning: the iterative solver did not converge (continuing)." );
        pri( "Number of iterations", iter );
        pri( "Relative residual |b-Ax|/|b|", real_residual );
        ready = 1;
        break;
      }
      pri( "" );
      pri( "Error: the iterative solver did not converge." );
      pri( "The initial error in the linear equations is", initial_error );
      pri( "The final error in the linear equations is", error );
      pri( "The relative residual |b-Ax|/|b| is", real_residual );
      pri( "The number of iterations is", iter );
      pri( "That is not good enough.\n\n" );
      pri( "Probably you have a wrong input file." );
      pri( "- Not enough boundary conditions?" );
      pri( "- A singular or badly scaled system?" );
      pri( "- A diverged calculation?\n" );
      succesful = 0;
      break;
    }
    if ( solve_iterative_bicg_use_cg ) {
        // CG for the symmetric system: d = r, one direction vector.
        // alpha = |r|^2 / (d^T A d) with |r|^2 = error (recomputed);
        // the residual is updated recursively and recomputed from the
        // iterate on the next pass, which corrects rounding drift.
      dAd = array_inproduct( d1, Ad1, solve_nlocal );
      if ( swit ) pri( "dAd", dAd );
      if ( scalar_dabs(dAd)<EPS_dAd ) {
          // genuine breakdown: d^T A d ~ 0 with a non-converged residual.
          // For a symmetric matrix this means the system is (numerically)
          // singular or not positive definite - a real failure, not a
          // "success with x = 0" as in the old code.
        real_residual = solve_iterative_bicg_real_residual( residue, norm_b );
        if ( bicg_stop==-NO ) {
          pri( "Warning: the iterative solver broke down (d^T A d ~ 0, continuing)." );
          pri( "This usually means a singular or nearly singular system." );
          pri( "Relative residual |b-Ax|/|b|", real_residual );
          ready = 1;
          break;
        }
        pri( "" );
        pri( "Error: the iterative solver broke down." );
        pri( "d^T A d ~ 0 with a non-converged residual." );
        pri( "The system is (numerically) singular or not positive definite." );
        pri( "The initial error in the linear equations is", initial_error );
        pri( "The relative residual |b-Ax|/|b| is", real_residual );
        pri( "That is not good enough.\n\n" );
        pri( "Probably you have a wrong input file." );
        pri( "- Not enough boundary conditions?" );
        pri( "- An invalid element connectivity (negative Jacobians)?" );
        pri( "- A load in a zero-energy mode of the matrix?\n" );
        succesful = 0;
        break;
      }
      alpha = error / dAd;
      for ( ilocal=0; ilocal<solve_nlocal; ilocal++ ) {
        solve_x[ilocal] += alpha * d1[ilocal];
        r1[ilocal] -= alpha * Ad1[ilocal];
      }
      beta = array_inproduct( r1, r1, solve_nlocal ) / error;
      for ( ilocal=0; ilocal<solve_nlocal; ilocal++ ) {
        d1[ilocal] = r1[ilocal] + beta * d1[ilocal];
      }
      if ( swit ) {
        pri( "alpha", alpha );
        pri( "solve_x", solve_x, solve_nlocal );
        pri( "error", error );
        pri( "beta", beta );
        pri( "d1", d1, solve_nlocal );
      }
    }
    else {
        // Bi-CG for the non-symmetric system: two directions (A and A^T).
      size_r1r2 = array_inproduct( r1, r2, solve_nlocal );
      dAd = array_inproduct( d2, Ad1, solve_nlocal );
      if ( swit ) {
        pri( "dAd", dAd );
        pri( "size_r1r2", size_r1r2 );
      }
      if ( scalar_dabs(size_r1r2)<EPS_TMP || scalar_dabs(dAd)<EPS_dAd ) {
          // breakdown of the bi-Lanczos recurrence: the system is
          // (numerically) singular or the biorthogonalisation collapsed.
          // Reported as a real failure, not as "success".
        real_residual = solve_iterative_bicg_real_residual( residue, norm_b );
        if ( bicg_stop==-NO ) {
          pri( "Warning: the Bi-CG iteration broke down (r1.r2 or d^T A d ~ 0, continuing)." );
          pri( "This usually means a singular or nearly singular system." );
          pri( "Relative residual |b-Ax|/|b|", real_residual );
          ready = 1;
          break;
        }
        pri( "" );
        pri( "Error: the Bi-CG iteration broke down." );
        pri( "The bi-Lanczos recurrence collapsed (r1.r2 or d^T A d ~ 0)" );
        pri( "with a non-converged residual." );
        pri( "The system is (numerically) singular." );
        pri( "The initial error in the linear equations is", initial_error );
        pri( "The relative residual |b-Ax|/|b| is", real_residual );
        pri( "That is not good enough.\n\n" );
        pri( "Probably you have a wrong input file." );
        pri( "- Not enough boundary conditions?" );
        pri( "- An invalid element connectivity (negative Jacobians)?" );
        pri( "- A load in a zero-energy mode of the matrix?\n" );
        succesful = 0;
        break;
      }
      alpha = size_r1r2 / dAd;
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
  }
  length = 1;
  db( GLOBAL_SOLVER_ITERATIONS, 0, &iter, ddum, length, VERSION_NORMAL, PUT );
  db( GLOBAL_SOLVER_ERROR, 0, idum, &error, length, VERSION_NORMAL, PUT );
  if ( print_solver==-YES ) {
    pri( "initial error", initial_error );
    pri( "check error", check_error );
    pri( "final error", error );
    pri( "relative residual |b-Ax|/|b|",
      solve_iterative_bicg_real_residual( residue, norm_b ) );
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

// Real (unpreconditioned) relative residual |b - A x| / |b| of the current
// iterate. residue holds the PRECONDITIONED residual P (b - A x) with
// P = diag(p); the unpreconditioned residual is residue[i]/p[i]. Used only
// for honest diagnostics (the convergence test itself works in the
// preconditioned norm the algorithm iterates in).
double solve_iterative_bicg_real_residual( double *residue, double norm_b )

{
  long int ilocal=0;
  double sum=0., tmp=0.;

  if ( norm_b<=0. ) return 0.;
  for ( ilocal=0; ilocal<solve_nlocal; ilocal++ ) {
    tmp = residue[ilocal] / p[ilocal];
    sum += tmp*tmp;
  }
  return sqrt(sum)/norm_b;
}

// One serial pass over the assembled element matrices: measures whether the
// system the iterative solver iterates on is symmetric. The preconditioned
// matrix P A P is symmetric iff A is (P is diagonal). The global matrix is
// symmetric iff every element matrix is symmetric (the diagonal NODE_LHSIDE
// block is diagonal) - and the element matrices are stored dense (both
// triangles), so the symmetry error is measured per element per pair.
// Returns 1 (use CG) when for every pair of free dofs (i,j)
//   |Aij - Aji| <= EPS_SYMMETRY * max(|Aij|,|Aji|)
// and 0 (use Bi-CG) otherwise.
long int solve_iterative_bicg_symmetric( void )

{
  long int element=0, max_element=0, max_node=0, i=0, iglobal=0, jglobal=0,
    ilocal=0, jlocal=0, element_group=0, ldum=0, element_matrix_values_length=0,
    use_element_matrix=0, use_group_matrix=0, inol=0, inod=0, jnod=0,
    ipuknwn=0, jpuknwn=0, m=0, indx=0, length_el=0, nnol=0, ii=0, jj=0,
    symmetric=1, pos_i=0, pos_j=0;
  long int *element_matrix_unknowns=NULL, *el=NULL, *node_pos=NULL;
  double ddum[1], *element_matrix_values=NULL, *dense=NULL,
    maxabs=0., asym=0., rel=0., scale=0., val_i=0., val_j=0.;

  db_highest_index( ELEMENT, max_element, VERSION_NORMAL );
  db_highest_index( NODE, max_node, VERSION_NORMAL );
  if ( max_element<0 || max_node<0 ) return 1;
  if ( npuknwn>MAX_PUKNWN_SYMMETRY ) return 0;

  dense = get_new_dbl( MNOL*npuknwn*MNOL*npuknwn );
  node_pos = get_new_int( 1+max_node );
  el = get_new_int( 1+MNOL );

  for ( element=0; element<=max_element; element++ ) {
    if ( !db_active_index( ELEMENT, element, VERSION_NORMAL ) ) continue;
    element_group = 0;
    db( ELEMENT_GROUP, element, &element_group, ddum, ldum, 
      VERSION_NORMAL, GET_IF_EXISTS );
    use_element_matrix = 0; use_group_matrix = 0;
    if ( db_active_index( ELEMENT_MATRIX_VALUES, element, VERSION_NORMAL ) )
      use_element_matrix = 1;
    else if ( db_active_index( GROUP_MATRIX_VALUES, element_group, VERSION_NORMAL ) )
      use_group_matrix = 1;
    if ( !use_element_matrix && !use_group_matrix ) continue;

    db( ELEMENT, element, el, ddum, length_el, VERSION_NORMAL, GET );
    nnol = length_el - 1;
    for ( inol=0; inol<nnol; inol++ )
      node_pos[el[1+inol]] = inol;
    m = nnol * npuknwn;
    array_set( dense, 0., m*m );

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
        inod = iglobal/npuknwn; ipuknwn = iglobal%npuknwn;
        jnod = jglobal/npuknwn; jpuknwn = jglobal%npuknwn;
        pos_i = node_pos[inod]; pos_j = node_pos[jnod];
        indx = (pos_i*npuknwn + ipuknwn) * m + (pos_j*npuknwn + jpuknwn);
        dense[indx] += element_matrix_values[i];
      }
    }
    if ( use_group_matrix ) delete[] element_matrix_unknowns;
    for ( inol=0; inol<nnol; inol++ )
      node_pos[el[1+inol]] = -1;

    for ( ii=0; ii<m; ii++ ) {
      for ( jj=ii+1; jj<m; jj++ ) {
        val_i = dense[ii*m+jj];
        val_j = dense[jj*m+ii];
        scale = scalar_dabs(val_i);
        if ( scalar_dabs(val_j)>scale ) scale = scalar_dabs(val_j);
        if ( scale>EPS_P ) {
          asym = scalar_dabs(val_i-val_j);
          rel = asym/scale;
          if ( rel>maxabs ) {
            maxabs = rel;
            if ( rel>EPS_SYMMETRY ) symmetric = 0;
          }
        }
      }
    }
  }

  delete[] dense;
  delete[] node_pos;
  delete[] el;

  return symmetric;
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
          if ( !solve_iterative_bicg_use_cg )
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
      if ( !solve_iterative_bicg_use_cg )
        Ad2[ilocal] += Ad2_thread[ithread*solve_nlocal+ilocal];
      p_tmp[ilocal] += p_thread[ithread*solve_nlocal+ilocal];
      residue[ilocal] += residue_thread[ithread*solve_nlocal+ilocal];
    }
    if ( scalar_dabs(p_tmp[ilocal])<EPS_P ) p_tmp[ilocal] = 1.;
  }
    // the residue of the current iterate, in the PRIMAL system:
    // residue = b_hat - A_hat x_hat (preconditioned). For a symmetric
    // matrix this equals the transpose residual the original code
    // accumulated; for a non-symmetric matrix it is the honest monitor
    // of the system actually being solved (the old accumulation
    // measured |b_hat - A_hat^T x_hat|, which never vanishes at the
    // solution of a non-symmetric system - e.g. the beam formulation).
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
        if ( !solve_iterative_bicg_use_cg )
          Ad2_thread[indx2] += element_matrix_values[i] * d2[ilocal] * 
            ( p[ilocal] * p[jlocal] ) ;  
        residue_thread[indx1] -= element_matrix_values[i] * solve_x[jlocal] * 
          ( p[ilocal] * p[jlocal] ) ;  
        if ( ilocal==jlocal ) p_thread[indx1] += element_matrix_values[i];
      }
    }
    if ( use_group_matrix ) delete[] element_matrix_unknowns;
  }

}
