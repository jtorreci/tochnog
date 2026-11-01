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
#include <algorithm>

double *d1, *d2, *p, *Ad1_thread, *Ad2_thread, *p_thread, *residue_thread;
long int solve_iterative_bicg_use_cg = 0;
// F4 T6f C74: set by solve_iterative_bicg() when the Bi-CG recurrence broke
// down (returned 0 from the breakdown exit, solver_bicg_stop -yes); the
// diagnostics are kept for solve_iterative_bicg_breakdown_report().
long int solve_iterative_bicg_breakdown = 0;
static double bicg_breakdown_initial_error = 0., bicg_breakdown_real_residual = 0.;
#define EPS_dAd 1.e-16
// FIX A (2026-09-13): scale-relative breakdown test. |d^T A d| is compared
// against the CURRENT residual energy error = ||residue||^2 instead of the
// absolute 1e-16. The CG search direction is built from the residual, so
// |d^T A d|/error is the Rayleigh quotient of the unit-diagonal
// preconditioned operator: O(1) for a well-conditioned system and -> 0 only
// when the residual lies in a numerical null space. The old absolute test
// fired near convergence (error itself ~1e-16) on perfectly solvable
// systems (dynamic8 second DEV-B solve: dAd ~ 5e-17, error ~ 2e-16).
#define EPS_dAd_REL 1.e-14
#define EPS_TMP 1.e-16
// F4 T6f C77 relative stop (see solve_iterative_bicg)
#define BICG_FLOOR_FACTOR 100.
#define BICG_BREAKDOWN_REL 1.e-14
#define BICG_FLOOR_ACCEPT 1.e-6
#define EPS_P 1.e-10
#define EPS_SYMMETRY 1.e-8
#define MAX_PUKNWN_SYMMETRY 32

long int solve_iterative_bicg_symmetric( void );
double solve_iterative_bicg_real_residual( double *residue, double norm_b );

// ---------------------------------------------------------------------------
// Assembled CSR operator of the iterative solver (F4 T9, 2026-10-10).
//
// The historic Bi-CG/CG product is matrix free, element by element: every
// matrix-vector product walks all element matrices, maps every entry
// through solve_global_local and scatters into per-thread buffers. On the
// Professional force11 (one 3D solve, ~10^5 iterations) that costs >1 h.
// Here the element matrices are assembled once per solve into a compressed
// sparse row matrix (rows and columns = local solver dofs, the NODE_LHSIDE
// diagonal kept apart), and the products of solve_iterative_bicg_sys() run
// on it. The pattern is built without a triplet buffer: the node set of
// every element (the nodes of its matrix unknowns that are solver dofs),
// the node -> element incidence, and per node the union of the node sets
// of its elements give the column list shared by all dofs of the node.
// The products are the same sums in another order (round-off only).
//
// The element-by-element product stays the fallback: with the mpc tie
// elimination active (DIAG-SOLVE-MIXTO §16: small systems only), when an
// entry is not found in the pattern, and with the diagnostic environment
// variable TN_BICG_EBE set.
static long int bicg_csr_on = 0;
static long int *bicg_csr_rowptr = NULL;
static int *bicg_csr_col = NULL;
static double *bicg_csr_val = NULL;
static double *bicg_csr_lhside = NULL;
// 1 for the dofs of a node with NODE_LHSIDE: the element-by-element product
// adds the right-hand side to the residue in that node loop only
static double *bicg_csr_bfac = NULL;
static double *bicg_csr_Ad1 = NULL, *bicg_csr_Ad2 = NULL,
  *bicg_csr_p_tmp = NULL, *bicg_csr_residue = NULL;
static long int bicg_csr_nblock = 1;

// Visit the matrix entries of one element that map to solver dofs, with the
// same element/group matrix selection as solve_iterative_bicg_element():
// f( ilocal, jlocal, iglobal, jglobal, value ).
template <class F> static void bicg_csr_visit_element( long int element,
  long int *unk_buf, F f )

{
  long int i=0, element_group=0, ldum=0, length=0, use_element_matrix=0,
    use_group_matrix=0;
  long int *unk=NULL;
  double ddum[1], *val=NULL;

  if ( db_active_index( ELEMENT_MATRIX_VALUES, element, VERSION_NORMAL ) )
    use_element_matrix = 1;
  else if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) ) {
    db( ELEMENT_GROUP, element, &element_group, ddum, ldum,
      VERSION_NORMAL, GET_IF_EXISTS );
    if ( db_active_index( GROUP_MATRIX_VALUES, element_group, VERSION_NORMAL ) )
      use_group_matrix = 1;
  }
  if ( !use_element_matrix && !use_group_matrix ) return;
  if ( use_element_matrix ) {
    unk = db_int( ELEMENT_MATRIX_UNKNOWNS, element, VERSION_NORMAL );
    val = db_dbl( ELEMENT_MATRIX_VALUES, element, VERSION_NORMAL );
    length = db_len( ELEMENT_MATRIX_VALUES, element, VERSION_NORMAL );
  }
  else {
    unk = unk_buf;
    get_element_matrix_unknowns( element, unk );
    val = db_dbl( GROUP_MATRIX_VALUES, element_group, VERSION_NORMAL );
    length = db_len( GROUP_MATRIX_VALUES, element_group, VERSION_NORMAL );
  }
  for ( i=0; i<length; i++ ) {
    long int iglobal = unk[i*2+0], jglobal = unk[i*2+1];
    long int ilocal = solve_global_local[iglobal];
    long int jlocal = solve_global_local[jglobal];
    if ( ilocal!=-NO && jlocal!=-NO )
      f( ilocal, jlocal, iglobal, jglobal, val[i] );
  }
}

static void bicg_csr_free( void )

{
  delete[] bicg_csr_rowptr; bicg_csr_rowptr = NULL;
  delete[] bicg_csr_col; bicg_csr_col = NULL;
  delete[] bicg_csr_val; bicg_csr_val = NULL;
  delete[] bicg_csr_lhside; bicg_csr_lhside = NULL;
  delete[] bicg_csr_bfac; bicg_csr_bfac = NULL;
  bicg_csr_on = 0;
}

// Build the CSR operator of the current solve; returns 1 on success (and
// sets bicg_csr_on), 0 when the element-by-element product must be used.
static long int bicg_csr_build( void )

{
  long int n = solve_nlocal, max_element=-1, max_node=-1, element=0,
    inod=0, ipuknwn=0, k=0, nthread=1, ldum=0, ok=1;
  double ddum[1];

  bicg_csr_on = 0;
  if ( n<=0 ) return 0;
  if ( getenv( "TN_BICG_EBE" )!=NULL ) return 0;
  if ( mpc_tie_elimination_active() ) return 0;
  db_highest_index( ELEMENT, max_element, VERSION_NORMAL );
  db_highest_index( NODE, max_node, VERSION_NORMAL );
  if ( max_element<0 || max_node<0 ) return 0;
  db( OPTIONS_PROCESSORS, 0, &nthread, ddum, ldum, VERSION_NORMAL, GET );
  bicg_csr_nblock = ( nthread>1 ) ? nthread : 1;

  long int nnode = 1 + max_node;
  long int *unk_buf = get_new_int( 2*MNOL*nprinc*MNOL*nprinc );
  long int *mark = get_new_int( nnode );
  for ( inod=0; inod<nnode; inod++ ) mark[inod] = -1;

  // 1. node set of every element (nodes of entries that are solver dofs)
  long int *el_ptr = get_new_int( max_element+2 );
  std::vector<long int> el_nodes;
  el_ptr[0] = 0;
  for ( element=0; element<=max_element; element++ ) {
    bicg_csr_visit_element( element, unk_buf,
      [&]( long int, long int, long int ig, long int jg, double ) {
        long int a = ig/npuknwn, b = jg/npuknwn;
        if ( mark[a]!=element ) { mark[a] = element; el_nodes.push_back( a ); }
        if ( mark[b]!=element ) { mark[b] = element; el_nodes.push_back( b ); }
      } );
    el_ptr[element+1] = (long int) el_nodes.size();
  }

  // 2. node -> element incidence
  long int *nd_ptr = get_new_int( nnode+1 );
  for ( inod=0; inod<=nnode; inod++ ) nd_ptr[inod] = 0;
  for ( k=0; k<(long int)el_nodes.size(); k++ ) nd_ptr[el_nodes[k]+1]++;
  for ( inod=0; inod<nnode; inod++ ) nd_ptr[inod+1] += nd_ptr[inod];
  long int *nd_el = get_new_int( nd_ptr[nnode]>0 ? nd_ptr[nnode] : 1 );
  {
    long int *fill = get_new_int( nnode );
    for ( inod=0; inod<nnode; inod++ ) fill[inod] = nd_ptr[inod];
    for ( element=0; element<=max_element; element++ )
      for ( k=el_ptr[element]; k<el_ptr[element+1]; k++ )
        nd_el[fill[el_nodes[k]]++] = element;
    delete[] fill;
  }

  // 3. row lengths: a dof row of node a holds the solver dofs of every node
  //    sharing an element with a
  long int *ndof_node = get_new_int( nnode );
  for ( inod=0; inod<nnode; inod++ ) {
    ndof_node[inod] = 0;
    mark[inod] = -1;
    for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ )
      if ( solve_global_local[inod*npuknwn+ipuknwn]!=-NO ) ndof_node[inod]++;
  }
  bicg_csr_rowptr = new long int[n+1];
  for ( k=0; k<=n; k++ ) bicg_csr_rowptr[k] = 0;
  for ( inod=0; inod<nnode; inod++ ) {
    if ( ndof_node[inod]==0 ) continue;
    long int len = 0;
    for ( long int ie=nd_ptr[inod]; ie<nd_ptr[inod+1]; ie++ ) {
      element = nd_el[ie];
      for ( k=el_ptr[element]; k<el_ptr[element+1]; k++ ) {
        long int b = el_nodes[k];
        if ( mark[b]!=inod ) { mark[b] = inod; len += ndof_node[b]; }
      }
    }
    for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
      long int il = solve_global_local[inod*npuknwn+ipuknwn];
      if ( il!=-NO ) bicg_csr_rowptr[il+1] = len;
    }
  }
  for ( k=0; k<n; k++ ) bicg_csr_rowptr[k+1] += bicg_csr_rowptr[k];
  long int nnz = bicg_csr_rowptr[n];

  // 4. sorted column lists (one per node, copied to each of its dof rows)
  bicg_csr_col = new int[nnz>0 ? nnz : 1];
  bicg_csr_val = new double[nnz>0 ? nnz : 1];
  for ( k=0; k<nnz; k++ ) bicg_csr_val[k] = 0.;
  for ( inod=0; inod<nnode; inod++ ) mark[inod] = -1;
  std::vector<int> cols;
  for ( inod=0; inod<nnode; inod++ ) {
    if ( ndof_node[inod]==0 ) continue;
    cols.clear();
    for ( long int ie=nd_ptr[inod]; ie<nd_ptr[inod+1]; ie++ ) {
      element = nd_el[ie];
      for ( k=el_ptr[element]; k<el_ptr[element+1]; k++ ) {
        long int b = el_nodes[k];
        if ( mark[b]==inod ) continue;
        mark[b] = inod;
        for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
          long int jl = solve_global_local[b*npuknwn+ipuknwn];
          if ( jl!=-NO ) cols.push_back( (int) jl );
        }
      }
    }
    std::sort( cols.begin(), cols.end() );
    for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
      long int il = solve_global_local[inod*npuknwn+ipuknwn];
      if ( il==-NO ) continue;
      if ( bicg_csr_rowptr[il+1]-bicg_csr_rowptr[il]!=(long int)cols.size() )
        ok = 0;
      else
        for ( k=0; k<(long int)cols.size(); k++ )
          bicg_csr_col[bicg_csr_rowptr[il]+k] = cols[k];
    }
  }
  delete[] nd_el;
  delete[] nd_ptr;
  delete[] el_ptr;
  std::vector<long int>().swap( el_nodes );

  // 5. values (element entries summed in element order)
  if ( ok ) {
    for ( element=0; element<=max_element && ok; element++ ) {
      bicg_csr_visit_element( element, unk_buf,
        [&]( long int il, long int jl, long int, long int, double v ) {
          int *first = bicg_csr_col + bicg_csr_rowptr[il];
          int *last = bicg_csr_col + bicg_csr_rowptr[il+1];
          int *pos = std::lower_bound( first, last, (int) jl );
          if ( pos==last || *pos!=(int) jl ) { ok = 0; return; }
          bicg_csr_val[pos-bicg_csr_col] += v;
        } );
    }
  }

  // 6. NODE_LHSIDE diagonal per solver dof
  bicg_csr_lhside = new double[n];
  bicg_csr_bfac = new double[n];
  for ( k=0; k<n; k++ ) bicg_csr_lhside[k] = bicg_csr_bfac[k] = 0.;
  db_highest_index( NODE_LHSIDE, max_node, VERSION_NORMAL );
  for ( inod=0; inod<=max_node; inod++ ) {
    if ( !db_active_index( NODE_LHSIDE, inod, VERSION_NORMAL ) ) continue;
    double *node_lhside = db_dbl( NODE_LHSIDE, inod, VERSION_NORMAL );
    for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
      long int il = solve_global_local[inod*npuknwn+ipuknwn];
      if ( il!=-NO ) {
        bicg_csr_lhside[il] += node_lhside[ipuknwn];
        bicg_csr_bfac[il] = 1.;
      }
    }
  }

  delete[] ndof_node;
  delete[] mark;
  delete[] unk_buf;
  if ( !ok ) {
    bicg_csr_free();
    return 0;
  }
  bicg_csr_on = 1;
  return 1;
}

// One block of rows of the CSR products (see solve_iterative_bicg_sys):
// Ad1 = P A P d1, Ad2 = P A^T P d2 (Bi-CG only, scattered into the thread
// buffer), residue = P b - P A P x, p_tmp = diag(A) (+ NODE_LHSIDE).
static void bicg_csr_rows( long int iblock, long int ithread )

{
  long int n = solve_nlocal, i=0, k=0;
  long int chunk = ( n + bicg_csr_nblock - 1 ) / bicg_csr_nblock;
  long int i0 = iblock*chunk, i1 = i0 + chunk;
  double *Ad2t = Ad2_thread + ithread*n;

  if ( i1>n ) i1 = n;
  for ( i=i0; i<i1; i++ ) {
    double pi = p[i], a1 = 0., rx = 0., dg = bicg_csr_lhside[i],
      d2pi = d2[i]*pi;
    for ( k=bicg_csr_rowptr[i]; k<bicg_csr_rowptr[i+1]; k++ ) {
      long int j = bicg_csr_col[k];
      double v = bicg_csr_val[k], vpj = v*p[j];
      a1 += vpj*d1[j];
      rx += vpj*solve_x[j];
      if ( j==i ) dg += v;
      if ( !solve_iterative_bicg_use_cg ) Ad2t[j] += vpj*d2pi;
    }
    double lpp = bicg_csr_lhside[i]*pi*pi;
    bicg_csr_Ad1[i] = a1*pi + lpp*d1[i];
    bicg_csr_residue[i] = bicg_csr_bfac[i]*solve_b[i]*pi - rx*pi -
      lpp*solve_x[i];
    bicg_csr_p_tmp[i] = dg;
  }
}

static void parallel_solve_iterative_bicg_csr( void )

{
  long int iloop=0, nloop=0, ithread=0, *next_of_loop=NULL;

  next_of_loop = get_new_int( bicg_csr_nblock );
  parallel_sys_next_of_loop( next_of_loop, bicg_csr_nblock-1, nloop, ithread );
  for ( iloop=0; iloop<nloop; iloop++ ) {
    if ( next_of_loop[iloop]>bicg_csr_nblock-1 ) break;
    bicg_csr_rows( next_of_loop[iloop], ithread );
  }
  delete[] next_of_loop;
}

static void bicg_csr_sys( double *Ad1, double *Ad2, double *p_tmp,
  double *residue )

{
  long int n = solve_nlocal, i=0, ithread=0, nthread=1, ldum=0;
  double ddum[1];

  db( OPTIONS_PROCESSORS, 0, &nthread, ddum, ldum, VERSION_NORMAL, GET );
  if ( !solve_iterative_bicg_use_cg )
    array_set( Ad2_thread, 0., nthread*n );
  bicg_csr_Ad1 = Ad1; bicg_csr_Ad2 = Ad2;
  bicg_csr_p_tmp = p_tmp; bicg_csr_residue = residue;
  if ( bicg_csr_nblock>1 )
    parallel_sys_routine( &parallel_solve_iterative_bicg_csr );
  else
    bicg_csr_rows( 0, 0 );
  for ( i=0; i<n; i++ ) {
    if ( !solve_iterative_bicg_use_cg ) {
      double a2 = bicg_csr_lhside[i]*d2[i]*p[i]*p[i];
      for ( ithread=0; ithread<nthread; ithread++ )
        a2 += Ad2_thread[ithread*n+i];
      Ad2[i] = a2;
    }
    else
      Ad2[i] = 0.;
    if ( scalar_dabs(p_tmp[i])<EPS_P ) p_tmp[i] = 1.;
  }
}
// ---------------------------------------------------------------------------

long int solve_iterative_bicg( void )

{
  long int icontrol=0, iter=0, max_iter=0, ilocal=0, nthread=0, ready=0, 
    max_node=0, succesful=1, print_solver=-NO, length=0, ldum=0, 
    swit=0, idum[1], bicg_stop=-YES, ldum_bs=0;
  double error=0., alpha=0., beta=0., dAd=0., check_error=0.,
   bicg_error_minimum=1.e-16, size_r1r2=0., initial_error=0., 
   bicg_error=1.e-14, ddum[1], ddum_bs[1], norm_b=0., real_residual=0.,
   *r1=NULL, *r2=NULL, *Ad1=NULL, *Ad2=NULL, *p_tmp=NULL, *residue=NULL;

  swit = set_swit(-1,-1,"solve_iterative_bicg");
  solve_iterative_bicg_breakdown = 0;
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

    // assembled CSR operator (F4 T9); element-by-element when unavailable
  bicg_csr_build();
  if ( print_solver==-YES ) {
    if ( bicg_csr_on )
      cout << "iterative solver: assembled CSR operator, rows " <<
        solve_nlocal << " nonzeros " << bicg_csr_rowptr[solve_nlocal] << "\n";
    else
      cout << "iterative solver: element-by-element operator, rows " <<
        solve_nlocal << "\n";
  }

    // set diagonal preconditioner
  solve_iterative_bicg_sys( Ad1, Ad2, p_tmp, residue );
  for ( ilocal=0; ilocal<solve_nlocal; ilocal++ ) 
    p[ilocal] = 1./(sqrt(scalar_dabs(p_tmp[ilocal])));

    // initial error
  solve_iterative_bicg_sys( Ad1, Ad2, p_tmp, residue );
  initial_error = array_inproduct( residue, residue, solve_nlocal );

    // relative check error
  long int mpc_elim_tighten = 0;
  double eps_stop = EPS_TMP;
  if ( mpc_tie_elimination_active() && solve_nlocal<=MPC_TIE_TIGHT_NMAX &&
       !db_active_index( CONTROL_OPTIONS_SOLVER_BICG_ERROR, icontrol,
         VERSION_NORMAL ) &&
       !db_active_index( SOLVER_BICG_ERROR, 0, VERSION_NORMAL ) &&
       !db_active_index( CONTROL_OPTIONS_SOLVER_BICG_ERROR_MINIMUM, icontrol,
         VERSION_NORMAL ) &&
       !db_active_index( OPTIONS_SOLVER_BICG_ERROR_MINIMUM, 0,
         VERSION_NORMAL ) )
    mpc_elim_tighten = 1;
  db( CONTROL_OPTIONS_SOLVER_BICG_ERROR, icontrol, idum, &bicg_error,
	  ldum, VERSION_NORMAL, GET_IF_EXISTS );
  // solver_bicg_error (manual Professional 6.1048): the global record
  // overwrites the per-control one
  db( SOLVER_BICG_ERROR, 0, idum, &bicg_error,
	  ldum, VERSION_NORMAL, GET_IF_EXISTS );

  // DIAG-SOLVE-MIXTO §21 (dev/solve-map): the §16 energy-consistent tie
  // elimination turns the staggered momentum iteration into a map whose
  // fixed point carries EXACTLY the error of the first inner solve: the
  // outer step v += A^-1*(b(v)) contracts the previous error and then
  // STALLS at the absolute floor bicg_error_minimum (the next correction
  // rhs is already below it and the solver returns x = 0). The default
  // bicg_error = 1e-14 is a 1e-7 RELATIVE residual: measured on mpc3
  // (probe, print_solver -yes), the first solve stopped at
  // |b-Ax|/|b| = 4.28546e-08 -> ~1e-8 velocity noise, the second solve
  // was floor-limited at 9.3% of its (tiny) correction, and the map
  // stalled with disy = 1/3 - 2.2e-9. That is the intrinsic
  // ~1e-8/iteration inconsistency of the eliminated staggered map
  // (DIAG-SOLVE-MIXTO §18.3). For ELIMINATED systems both defaults are
  // therefore tightened: bicg_error = 1e-24 (a 1e-12 relative residual)
  // and bicg_error_minimum = 1e-26 (the floor only binds on the LATER
  // correction solves whose rhs is already ~1e-12 scaled; it stays far
  // above the recomputed-residual recursion noise (eps*|A|*|x|)^2 for any
  // realistic conditioning, and a system that truly cannot reach it still
  // has the honest max_iter exit icw the direct-LU retry of so.cc). The
  // genuine-convergence floor of the loop drops from EPS_TMP (1e-16, the
  // effective stop for every other solve) to the same 1e-26. The
  // eliminated map then converges to the precision of its first inner
  // solve (~1e-12 relative), which is Professional-grade for every corpus
  // target. The tightening is capped at MPC_TIE_TIGHT_NMAX equations: the
  // corpus eliminated systems (mpc3/4/5/6/7/9, force16) are small, while a
  // large eliminated system does not need the deep precision (its targets
  // are orders looser) and must not pay the breakdown-then-direct-LU-retry
  // cost on a big band matrix (measured: ground19 3min6s uncapped vs 14.6s
  // baseline). Non-eliminated paths are untouched: the gate is the active
  // tie map, and explicit user solver_bicg_error records always win.
  if ( mpc_elim_tighten ) {
    bicg_error = 1.e-24;
    bicg_error_minimum = 1.e-26;
    eps_stop = 1.e-26;
  }

  check_error = bicg_error * initial_error;
  if ( swit ) pri( "check_error", check_error );

    // minimum check error (F4 T6f C77). bicg_error_minimum is an ABSOLUTE
    // floor on the squared preconditioned residual (GNU record, not in the
    // Professional manual). Before C77 it was always active with the
    // default 1e-16, together with eps_stop = 1e-16: a solve with a small
    // right-hand side (a Newton correction, a flow near its steady state)
    // then stopped far from the relative tolerance (true |b-Ax|/|b| up to
    // 0.13 measured). Now the floor applies only when the user gives the
    // record (options_solver_bicg_error_minimum, overwritten by
    // control_options_solver_bicg_error_minimum) or for the tightened
    // eliminated tie systems; otherwise the stop is purely relative, with
    // the round-off floor eps*|A|*|x| below.
  if ( !mpc_elim_tighten ) {
    long int has_min = 0;
    if ( db( OPTIONS_SOLVER_BICG_ERROR_MINIMUM, 0, idum, &bicg_error_minimum,
           ldum, VERSION_NORMAL, GET_IF_EXISTS ) ) has_min = 1;
    if ( db( CONTROL_OPTIONS_SOLVER_BICG_ERROR_MINIMUM, icontrol, idum,
           &bicg_error_minimum, ldum, VERSION_NORMAL, GET_IF_EXISTS ) )
      has_min = 1;
    if ( !has_min ) bicg_error_minimum = 0.;
    eps_stop = 0.;
  }
  if ( check_error<bicg_error_minimum ) check_error = bicg_error_minimum;

    // solver_bicg_stop (manual Professional 6.1051) and
    // control_solver_bicg_stop (6.376): -no does NOT stop the
    // calculation when the iterative solver does not converge (continue
    // with the current solution); -yes (and the default) stops.
    // Precedence (project decision 2026-09-29, see
    // manual-developer/solve_iterative_bicg.md): the global record, then
    // the control record with index 0, then the control record of the
    // current timestep index (last read wins).
    // Measured on Professional: the control record overrules the global
    // one (as the manual says), but ONLY the index-0 control record is
    // read, for every timestep. The manual documents a per-timestep
    // index, so a record with the timestep's own index is honoured here
    // (deliberate divergence); the index-0 fallback keeps Professional
    // inputs that rely on `control_solver_bicg_stop 0 ...` unchanged.
  db( SOLVER_BICG_STOP, 0, &bicg_stop,
    ddum_bs, ldum_bs, VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_SOLVER_BICG_STOP, 0, &bicg_stop,
    ddum_bs, ldum_bs, VERSION_NORMAL, GET_IF_EXISTS );
  if ( icontrol!=0 )
    db( CONTROL_SOLVER_BICG_STOP, icontrol, &bicg_stop,
      ddum_bs, ldum_bs, VERSION_NORMAL, GET_IF_EXISTS );

    // the velocity/temperature/pressure system assembled by the GNU is
    // symmetric when every element matrix is symmetric (the diagonal
    // NODE_LHSIDE block is diagonal). Measured per solve; symmetric ->
    // plain CG, non-symmetric (e.g. some coupled terms) -> honest Bi-CG.
  {
    long int matrix_symmetric = -NO;
    long int measured_symmetric = solve_iterative_bicg_symmetric( );
    db( SOLVER_MATRIX_SYMMETRIC, 0, &matrix_symmetric, ddum_bs, ldum_bs,
      VERSION_NORMAL, GET_IF_EXISTS );
    if ( matrix_symmetric==-YES ) {
      // solver_matrix_symmetric -yes (manual Professional 6.1052): the
      // matrices are symmetrized IF NEEDED so that a symmetric equation
      // solver can be used. The GNU does not re-symmetrize the assembled
      // matrix: when the measurement disagrees with the user's assertion
      // the honest Bi-CG runs on the as-assembled system (it converges
      // where plain CG on a non-symmetric system diverges). Deliberate
      // divergence K24 (DIVERGENCES-PROFESSIONAL.md, class A).
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
    //
    // F4 T6f C77 (owner decision 2026-10-15): the stop is RELATIVE.
    // - convergence: |r|^2 < bicg_error*|r0|^2 (x0 = 0, so |r|/|b| in the
    //   preconditioned norm), plus the user's absolute bicg_error_minimum
    //   when given;
    // - round-off floor: the recomputed residual cannot drop below about
    //   eps*|A|*|x|; when |r| < BICG_FLOOR_FACTOR*eps*|A|*|x| (|A| = the
    //   largest |A d|/|d| seen) the iteration stops there;
    // - breakdown: relative, |r1.r2| < BICG_BREAKDOWN_REL*|r1||r2| or
    //   |d2.Ad1| < BICG_BREAKDOWN_REL*|d2||Ad1| (Bi-CG), |dAd| <
    //   EPS_dAd_REL*|r|^2 (CG).
    //   A floor or breakdown stop is accepted when the true residual
    //   |b-Ax|/|b| <= BICG_FLOOR_ACCEPT; otherwise it is a breakdown and
    //   so.cc retries with a direct solver.
  ready=0; max_iter = 10*solve_nlocal; error = 1.e10;
  long int stop_kind = 0; // 1 round-off floor, 2 breakdown
  double norm_A = 0.;
  for ( iter=0; !ready; iter++ ) {
    if ( swit ) pri( "iterative solve iteration", iter );
    solve_iterative_bicg_sys( Ad1, Ad2, p_tmp, residue );
    error = array_inproduct( residue, residue, solve_nlocal );
    if ( print_solver==-YES ) {
      cout << iter << " " << error << "\n";
      cout << flush;
    }
    if ( error<check_error || error<eps_stop || error==0. ) {
        // genuine convergence: residual below the relative tolerance (or
        // the user's absolute floor; eps_stop only for the tightened
        // eliminated tie systems)
      ready = 1;
      break;
    }
    {
      double nd1 = array_inproduct( d1, d1, solve_nlocal );
      double nAd1 = array_inproduct( Ad1, Ad1, solve_nlocal );
      if ( nd1>0. && nAd1>norm_A*norm_A*nd1 ) norm_A = sqrt( nAd1/nd1 );
      double fl = BICG_FLOOR_FACTOR*DBL_EPSILON*norm_A;
      if ( iter>0 &&
           error<fl*fl*array_inproduct( solve_x, solve_x, solve_nlocal ) ) {
        stop_kind = 1;
        break;
      }
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
      if ( scalar_dabs(dAd)<EPS_dAd_REL*error ) {
          // breakdown: d^T A d ~ 0 relative to the current residual
          // energy with a non-converged residual. For a symmetric matrix
          // this means the system is (numerically) singular or not
          // positive definite.
        stop_kind = 2;
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
      double nr1 = sqrt( array_inproduct( r1, r1, solve_nlocal ) );
      double nr2 = sqrt( array_inproduct( r2, r2, solve_nlocal ) );
      double nd2 = sqrt( array_inproduct( d2, d2, solve_nlocal ) );
      double nAd1 = sqrt( array_inproduct( Ad1, Ad1, solve_nlocal ) );
      if ( scalar_dabs(size_r1r2)<=BICG_BREAKDOWN_REL*nr1*nr2 ||
           scalar_dabs(dAd)<=BICG_BREAKDOWN_REL*nd2*nAd1 ) {
          // breakdown of the bi-Lanczos recurrence (scale-relative test,
          // C77): the system is (numerically) singular or the
          // biorthogonalisation collapsed.
        stop_kind = 2;
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
  if ( stop_kind ) {
      // round-off floor or breakdown: accept the iterate when its true
      // residual is small enough, otherwise report a breakdown (so.cc
      // retries with a direct solver; solver_bicg_stop -no keeps it)
    real_residual = solve_iterative_bicg_real_residual( residue, norm_b );
    if ( real_residual<=BICG_FLOOR_ACCEPT ) {
      if ( swit ) pri( "stop at the round-off floor, accepted", real_residual );
    }
    else if ( bicg_stop==-NO ) {
      pri( "Warning: the iterative solver broke down (continuing)." );
      pri( "This usually means a singular or nearly singular system." );
      pri( "Relative residual |b-Ax|/|b|", real_residual );
    }
    else {
        // F4 T6f C74: the error block is NOT printed here. The caller
        // (so.cc) retries the solve with a direct solver; only when that
        // retry fails too (or when there is no retry) does it print the
        // block, with solve_iterative_bicg_breakdown_report(). A rescued
        // breakdown gets one warning line from the caller.
      solve_iterative_bicg_breakdown = 1;
      bicg_breakdown_initial_error = initial_error;
      bicg_breakdown_real_residual = real_residual;
      succesful = 0;
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
  bicg_csr_free();

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

  if ( bicg_csr_on ) {
    bicg_csr_sys( Ad1, Ad2, p_tmp, residue );
    if ( swit ) pri( "Out routine SOLVE_ITERATIVE_BICG_SYS (CSR)" );
    return;
  }

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
  long int i=0, iglobal=0, ilocal=0, jglobal=0, jlocal=0,
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
    // energy-consistent slave elimination: hoisted per element so the
    // (inactive) hot path costs one local test per entry, not a call
    long int mpc_tie_on = mpc_tie_elimination_active();
    for ( i=0; i<element_matrix_values_length; i++ ) {
      iglobal = element_matrix_unknowns[i*2+0];
      jglobal = element_matrix_unknowns[i*2+1];
      ilocal = solve_global_local[iglobal];
      jlocal = solve_global_local[jglobal];
      if ( mpc_tie_on &&
           ( mpc_tie_of_dof(iglobal)>=0 || mpc_tie_of_dof(jglobal)>=0 ) ) {
          // energy-consistent slave elimination (DIAG-SOLVE-MIXTO §16):
          // the operator the iterative solver runs on is the REDUCED
          // system - every slave row/column of the element matrix is
          // redistributed to its tie masters with the tie weights
          // (K_red = K_FF + K_FS*C^T + C*K_SF + C*K_SS*C^T). The slave
          // dofs are not unknowns of the solve (value-constrained), so
          // their matrix couplings cannot be dropped: dropping them would
          // lose the slave consistent force share of the non-conforming
          // interface and the fixed point would be value-constrained
          // (0.454 on mpc3) instead of the homogeneous field (1/3).
        long int row_t[MPC_TIE_MAX_TARGETS], col_t[MPC_TIE_MAX_TARGETS],
          it2=0;
        double fac_t[MPC_TIE_MAX_TARGETS];
        long int nt = mpc_tie_targets( iglobal, jglobal, row_t, col_t,
          fac_t );
        for ( it2=0; it2<nt; it2++ ) {
          long int il = solve_global_local[row_t[it2]];
          long int jl = solve_global_local[col_t[it2]];
          if ( il!=-NO && jl!=-NO ) {
            double v = element_matrix_values[i] * fac_t[it2];
            indx1 = ithread*solve_nlocal + il;
            indx2 = ithread*solve_nlocal + jl;
            Ad1_thread[indx1] += v * d1[jl] * ( p[il] * p[jl] );
            if ( !solve_iterative_bicg_use_cg )
              Ad2_thread[indx2] += v * d2[il] * ( p[il] * p[jl] );
            residue_thread[indx1] -= v * solve_x[jl] * ( p[il] * p[jl] );
            if ( il==jl ) p_thread[indx1] += v;
          }
        }
        continue;
      }
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

// F4 T6f C74: full diagnostic of a Bi-CG breakdown (r1.r2 or d^T A d ~ 0
// with a non-converged residual), printed by the caller only when the
// breakdown is NOT rescued by the direct band LU retry (so.cc). Before C74
// it was printed inside solve_iterative_bicg() for every breakdown, also
// for the ones the LU retry solved exactly (599 times on validation_1).
void solve_iterative_bicg_breakdown_report( void )

{
  pri( "" );
  pri( "Error: the iterative solver broke down." );
  pri( "The recurrence collapsed (a scalar product ~ 0 relative to its" );
  pri( "vectors) with a non-converged residual." );
  pri( "The system is (numerically) singular." );
  pri( "The initial error in the linear equations is", bicg_breakdown_initial_error );
  pri( "The relative residual |b-Ax|/|b| is", bicg_breakdown_real_residual );
  pri( "That is not good enough.\n\n" );
  pri( "Probably you have a wrong input file." );
  pri( "- Not enough boundary conditions?" );
  pri( "- An invalid element connectivity (negative Jacobians)?" );
  pri( "- A load in a zero-energy mode of the matrix?\n" );
}
