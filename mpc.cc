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

// ---------------------------------------------------------------------------
// Multi point constraints (MPC) - manual Professional 6.856-6.875.
//
// mpc_node_number / mpc_node_factor (6.874/6.875): the dof dof_0 of node
// node_0 (the slave) is constrained to a linear combination of dofs of the
// master nodes:
//
//     dof_0(node_0) = sum_i factor_i * dof_i(master_i)
//
// with the factors from mpc_node_factor (default 1). Base semantics
// verified against the Professional binary (25-10-2023, .dbs): the slave
// dof is treated as a KNOWN (bounded) dof whose value is recomputed every
// equilibrium iteration from the current master values. Boundary
// conditions must not be given on slave nodes (manual 6.875).
//
// MEASURED REFINEMENT (dev/fmpc sprint, 2026-09-08): the Professional does
// NOT drop the slave equilibrium row for ties whose slave AND all master
// dofs are free unknowns - it STATICALLY CONDENSES the tie into the master
// row/column (row_m += c*row_s with the constraint substituted, the
// DIAG-SOLVE-MIXTO §16 energy-consistent elimination). Measured on mpc9
// (explicit user record, condensed diagonal K22+K33 = 11/6 = 1.833333,
// v2 = 5/11) and on force16 (the mpc_element_group-generated tie v3=v2:
// value-constrained reproduces the GNU field 3.93667, condensed reproduces
// the Pro field 6.0342857 = 1056/175 in closed form). Ties that touch a
// DIRECTLY bounded dof stay value-constrained (mpc1: Dirichlet master,
// mpc8: bounda'd slave - both byte-identical before/after).
//
// dev/solve-map (DIAG §21): the eg-generated ties join the condensation
// (they were excluded in dev/fmpc pending the map-precision fix) and the
// generated FACTORS are polished to machine precision - point_el's
// finite-difference descent left ~2.4e-8 of factor error (the Professional
// writes exact 0.5 on the same input), which was the true intrinsic
// inconsistency of the eliminated staggered map.
//
// mpc_linear_quadratic (6.873): when a quadratic element (quad9, hex27,
// ...) borders a linear element (quad4, hex8, ...) at a common interface,
// the extra nodes of the quadratic element (mid-edge, mid-face) dangle:
// they are not attached to the linear element, producing non-compatible
// solution fields. The option generates mpc_node_number/mpc_node_factor
// records that tie each dangling node to the linear element with the shape
// functions of the linear element evaluated at the node position. The
// generated records are identical to the ones the Professional writes to
// the .dbs (verified on mpc3/mpc4/mpc5: mid-edge slaves get 2 masters with
// factor 0.5, mid-face slaves 4 masters with factor 0.25).
//
// The generated records are re-created whenever the mesh changes
// (refinement, deletion, splitting, ...): the mesh fingerprint is stored
// in the internal record MPC_LINEAR_QUADRATIC_MESH_FINGERPRINT together
// with the index range of the generated records.
// ---------------------------------------------------------------------------

// tolerance on the isoparametric coordinates below which a dangling node
// is considered to be located in the linear element (same default as
// mpc_element_group_eps_iso, manual 6.863)
#define MPC_EPS_ISO 1.e-4

// weights below this value are not written to the generated records
// (the Professional lists only the non-zero masters)
#define MPC_WEIGHT_ZERO 1.e-6

// ---------------------------------------------------------------------------
// Energy-consistent slave elimination inside the staggered solve
// (DIAG-SOLVE-MIXTO §16; the mpc3/mpc4 closure). The tie map below is built
// by mpc_node_apply() every equilibrium iteration from the SAME
// mpc_node_number/mpc_node_factor records it consumes for the value
// constraint, and is read by solve() (so.cc / so_bicg.cc): the slave
// row/column of the assembled momentum system and the slave right-hand side
// are redistributed to the tie masters with the tie weights
// (K_red = K_FF + K_FS*C^T + C*K_SF + C*K_SS*C^T), so the reduced system
// contains the slave constraint IMPLICITLY in the same solve that produces
// the masters. With the slave value that feeds the constitutive stress
// re-synced to C^T*masters of the current iterate (mpc_node_apply runs
// pre-element-loop and post-solve), the staggered momentum iteration
// collapses to one pass at the energy-consistent fixed point.
//
// The elimination is applied to the ties GENERATED by
// mpc_linear_quadratic (the non-conforming quadratic/linear interface, the
// mpc3/mpc4 failure), to the ties GENERATED by mpc_element_group (measured:
// the Professional condenses them; on the conforming mpc7 patch the
// condensed and value-constrained solutions coincide so the tight 1e-8
// target keeps its values) and to EXPLICIT user mpc_node_number records
// whose slave dof and ALL master dofs are free unknowns (measured: mpc9).
// Explicit records touching a directly bounded dof (mpc1 master, mpc8
// slave) stay value-constrained - measured byte-identical.
// Active ONLY when such ties exist; otherwise every accessor reports "no
// ties" and the solver code paths are byte-identical.
// ---------------------------------------------------------------------------

static long int mpc_tie_count = 0;
static long int mpc_tie_cap = 0;            // allocated tie slots
static long int *mpc_tie_slave_ig = NULL;   // [cap] slave global dof index
static long int *mpc_tie_off = NULL;        // [cap] offset in the flat arrays
static long int *mpc_tie_nm = NULL;         // [cap] number of masters
static long int mpc_tie_master_cap = 0;     // allocated master slots
static long int *mpc_tie_master_ig = NULL;  // [master_cap] master dof indices
static double  *mpc_tie_master_fac = NULL;  // [master_cap] tie factors
static long int mpc_tie_mark_size = 0;      // (1+max_node)*npuknwn of the build
static long int *mpc_tie_mark = NULL;       // [mark_size] 1+tie index per global
                                            //   dof, 0 = not a tied slave
static long int mpc_tie_elim_lo = -1;       // index range of the records
static long int mpc_tie_elim_hi = -1;       //   generated by mpc_linear_quadratic

// global dof indices of the slaves that the current bounda() pass applied
// through an ELIMINATED tie (the -veln generated records and the eligible
// explicit user records). The post-solve re-sync (top.cc, second
// mpc_node_apply of the iteration) clears their bounded flags first so it
// can re-derive which slaves the direct bounda records own (those are
// never applied and keep their direct value) and refresh the free ones
// with the freshly solved masters; the per-pass clearing also keeps the
// slave-free eligibility test of the NEXT pass valid (without it the
// bounded flag set by the previous pass would disqualify the tie and the
// elimination would be lost after the first iteration).
static long int mpc_applied_n = 0;
static long int mpc_applied_cap = 0;
static long int *mpc_applied = NULL;

static void mpc_tie_map_clear( void )

{
  delete[] mpc_tie_slave_ig;   mpc_tie_slave_ig = NULL;
  delete[] mpc_tie_off;        mpc_tie_off = NULL;
  delete[] mpc_tie_nm;         mpc_tie_nm = NULL;
  delete[] mpc_tie_master_ig;  mpc_tie_master_ig = NULL;
  delete[] mpc_tie_master_fac; mpc_tie_master_fac = NULL;
  delete[] mpc_tie_mark;       mpc_tie_mark = NULL;
  mpc_tie_count = 0;
  mpc_tie_cap = 0;
  mpc_tie_master_cap = 0;
  mpc_tie_mark_size = 0;
  mpc_tie_elim_lo = -1;
  mpc_tie_elim_hi = -1;
}

// register a slave dof applied through an eliminated tie in the list the
// post-solve re-sync uses to clear the mpc-owned bounded flags
static void mpc_applied_add( long int iglobal )

{
  long int i = 0;
  for ( i=0; i<mpc_applied_n; i++ )
    if ( mpc_applied[i]==iglobal ) return;
  if ( mpc_applied_n>=mpc_applied_cap ) {
    long int new_cap = ( mpc_applied_cap>0 ) ?
      2*mpc_applied_cap : 32;
    long int *tmp = get_new_int( new_cap );
    for ( i=0; i<mpc_applied_n; i++ ) tmp[i] = mpc_applied[i];
    delete[] mpc_applied;
    mpc_applied = tmp;
    mpc_applied_cap = new_cap;
  }
  mpc_applied[mpc_applied_n++] = iglobal;
}

// 1 = the energy-consistent slave elimination is active for this solve
// (there are mpc_linear_quadratic ties to eliminate)
// DIAG \u00a721: solve-size cap. On a large eliminated system the folded
// rows destabilize the iterative solve (measured: ground19, 6479 nodes, CG
// breakdown at 4.5e-6 relative + direct-LU retries, 13x runtime) while the
// condensed field is not needed there (targets orders looser than the small
// eliminated systems). solve() reports its size once per solve; beyond
// MPC_TIE_TIGHT_NMAX the elimination switches off and the ties keep the
// value-constrained semantics (the map bookkeeping below is unaffected).
static long int mpc_tie_solve_nlocal = 0;

void mpc_tie_solve_size( long int solve_nlocal )

{
  mpc_tie_solve_nlocal = solve_nlocal;
}

long int mpc_tie_elimination_active( void )

{
  return ( mpc_tie_count>0 && mpc_tie_mark!=NULL &&
           mpc_tie_solve_nlocal<=MPC_TIE_TIGHT_NMAX );
}

long int mpc_tie_count_get( void )

{
  return mpc_tie_count;
}

long int mpc_tie_slave_dof( long int tie )

{
  if ( tie<0 || tie>=mpc_tie_count ) return -1;
  return mpc_tie_slave_ig[tie];
}

// tie index (0-based) of a global dof, -1 when the dof is not a tied slave
long int mpc_tie_of_dof( long int iglobal )

{
  if ( mpc_tie_mark==NULL || iglobal<0 || iglobal>=mpc_tie_mark_size ) return -1;
  return mpc_tie_mark[iglobal] - 1;
}

long int mpc_tie_nmaster( long int tie )

{
  if ( tie<0 || tie>=mpc_tie_count ) return 0;
  return mpc_tie_nm[tie];
}

long int mpc_tie_master_dof( long int tie, long int im )

{
  if ( tie<0 || tie>=mpc_tie_count || im<0 || im>=mpc_tie_nm[tie] ) return -1;
  return mpc_tie_master_ig[mpc_tie_off[tie]+im];
}

double mpc_tie_master_factor( long int tie, long int im )

{
  if ( tie<0 || tie>=mpc_tie_count || im<0 || im>=mpc_tie_nm[tie] ) return 0.;
  return mpc_tie_master_fac[mpc_tie_off[tie]+im];
}

// Enumerate the image of the matrix entry (iglobal, jglobal) under the slave
// elimination u_slave = sum_m factor_m*u_master: a free row/column maps to
// itself, a tied slave row/column maps to its masters with the tie factors,
// and an entry with both dofs tied maps to all master pairs. Returns the
// number of targets (each target carries its combined factor). The callers
// (so.cc, so_bicg.cc) accumulate the entry value times the factor into their
// matrix-vector products / fills, skipping targets whose local index is -NO.
long int mpc_tie_targets( long int iglobal, long int jglobal,
  long int row_t[], long int col_t[], double fac_t[] )

{
  long int tie_i = mpc_tie_of_dof( iglobal );
  long int tie_j = mpc_tie_of_dof( jglobal );
  long int ni = ( tie_i>=0 ) ? mpc_tie_nm[tie_i] : 1;
  long int nj = ( tie_j>=0 ) ? mpc_tie_nm[tie_j] : 1;
  long int n = 0;
  for ( long int ii=0; ii<ni; ii++ ) {
    long int ig = ( tie_i>=0 ) ? mpc_tie_master_dof( tie_i, ii ) : iglobal;
    double fg = ( tie_i>=0 ) ? mpc_tie_master_factor( tie_i, ii ) : 1.;
    for ( long int jj=0; jj<nj; jj++ ) {
      long int jg = ( tie_j>=0 ) ? mpc_tie_master_dof( tie_j, jj ) : jglobal;
      double fj = ( tie_j>=0 ) ? mpc_tie_master_factor( tie_j, jj ) : 1.;
      if ( n>=MPC_TIE_MAX_TARGETS ) break;
      row_t[n] = ig; col_t[n] = jg; fac_t[n] = fg*fj;
      n++;
    }
  }
  return n;
}

// ---------------------------------------------------------------------------
// Machine-precision tie weights (DIAG-SOLVE-MIXTO §21, dev/solve-map).
// point_el() locates a point inside a tensor-product element with a
// finite-difference descent whose final accepted step is ~DELTA/2^12 =
// 2.44e-8 iso units: on the refined mpc3 mesh every generated factor comes
// out 0.499999975586 where the Professional writes 0.5 (measured, both
// .dbs of the corpus input). The factor error enters the tie CONSTRAINT
// itself and is the intrinsic ~1e-8 inconsistency of the eliminated (and
// equally of the value-constrained) staggered map: the momentum solves then
// converge exactly to a fixed point that is ~1e-8 away from the true
// condensed solution, no matter how tightly they converge. The Gauss-Newton
// polish below recomputes the tensor-product weights with the analytic
// derivatives of the same interpolation_polynomial bases point_el uses, so
// the generated records carry machine-precision factors. Non-tensor-product
// masters stay untouched: tria/tet factors are closed-form barycentric in
// point_el (exact), bar2 masters reduce to the single-node weight 1 (exact).
// ---------------------------------------------------------------------------

static void mpc_polish_dims( long int name, long int *nnx, long int *nny,
  long int *nnz )

{
  *nnz = 1;
  if      ( name==-QUAD4 )  { *nnx=2; *nny=2; }
  else if ( name==-QUAD9 )  { *nnx=3; *nny=3; }
  else if ( name==-QUAD16 ) { *nnx=4; *nny=4; }
  else if ( name==-HEX8 )   { *nnx=2; *nny=2; *nnz=2; }
  else if ( name==-HEX27 )  { *nnx=3; *nny=3; *nnz=3; }
  else if ( name==-HEX64 )  { *nnx=4; *nny=4; *nnz=4; }
  else                      { *nnx=0; *nny=0; *nnz=0; }
}

// solve the (small, symmetric positive) normal equations a*x = b by Gaussian
// elimination with partial pivoting; returns 0 on a singular system
static long int mpc_polish_solve( long int n, double a[], double b[] )

{
  for ( long int k=0; k<n; k++ ) {
    long int piv = k;
    for ( long int i=k+1; i<n; i++ )
      if ( scalar_dabs(a[i*n+k])>scalar_dabs(a[piv*n+k]) ) piv = i;
    if ( scalar_dabs(a[piv*n+k])<1.e-300 ) return 0;
    if ( piv!=k ) {
      for ( long int j=0; j<n; j++ ) {
        double t = a[k*n+j]; a[k*n+j] = a[piv*n+j]; a[piv*n+j] = t;
      }
      double t = b[k]; b[k] = b[piv]; b[piv] = t;
    }
    for ( long int i=k+1; i<n; i++ ) {
      double f = a[i*n+k]/a[k*n+k];
      for ( long int j=k; j<n; j++ ) a[i*n+j] -= f*a[k*n+j];
      b[i] -= f*b[k];
    }
  }
  for ( long int i=n-1; i>=0; i-- ) {
    double s = b[i];
    for ( long int j=i+1; j<n; j++ ) s -= a[i*n+j]*b[j];
    b[i] = s/a[i*n+i];
  }
  return 1;
}

// rebuild the weights of a tensor-product master at the isoparametric
// coordinates iso[] (same xi-inner/eta/zeta-outer node order as point_el)
static void mpc_polish_weights( double iso[], long int nnx, long int nny,
  long int nnz, double weight[] )

{
  double h_xi[MPOINT], p_xi[MPOINT], h_eta[MPOINT], p_eta[MPOINT],
    h_zeta[MPOINT], p_zeta[MPOINT];
  interpolation_polynomial( iso[0], nnx, h_xi, p_xi );
  interpolation_polynomial( iso[1], nny, h_eta, p_eta );
  interpolation_polynomial( iso[2], nnz, h_zeta, p_zeta );
  long int inol = 0;
  for ( long int iz=0; iz<nnz; iz++ )
    for ( long int iy=0; iy<nny; iy++ )
      for ( long int ix=0; ix<nnx; ix++ )
        weight[inol++] = h_xi[ix]*h_eta[iy]*h_zeta[iz];
}

// Gauss-Newton polish of the isoparametric inversion x(iso) = point for a
// tensor-product master; on success weight[] holds the machine-precision
// shape functions, on any failure the point_el weights are restored
static void mpc_weights_polish( double point[], double coords[],
  double weight[], long int name, long int nnol )

{
  long int nnx=0, nny=0, nnz=0;
  mpc_polish_dims( name, &nnx, &nny, &nnz );
  long int niso = ( nnx>1 ) + ( nny>1 ) + ( nnz>1 );
  if ( niso<1 || niso>ndim ) return;

  double h_xi[MPOINT], p_xi[MPOINT], h_eta[MPOINT], p_eta[MPOINT],
    h_zeta[MPOINT], p_zeta[MPOINT], iso[3], iso_best[3], weight_in[MPOINT],
    resid[MDIM], jac[MDIM*3], normal[9], rhside[3], delta[3];
  for ( long int i=0; i<nnol; i++ ) weight_in[i] = weight[i];

  // element scale for the convergence test
  double scale = 0.;
  for ( long int inol=0; inol<nnol; inol++ ) {
    double d2 = 0.;
    for ( long int idim=0; idim<ndim; idim++ ) {
      double d = coords[inol*ndim+idim];
      d2 += d*d;
    }
    if ( d2>scale ) scale = d2;
  }
  scale = sqrt( scale );
  if ( scale<=0. ) return;
  double tol = 1.e-14*scale;

  iso[0] = iso[1] = iso[2] = 0.;
  iso_best[0] = iso_best[1] = iso_best[2] = 0.;
  double err_best = 0.;
  for ( long int iter=0; iter<60; iter++ ) {
    interpolation_polynomial( iso[0], nnx, h_xi, p_xi );
    interpolation_polynomial( iso[1], nny, h_eta, p_eta );
    interpolation_polynomial( iso[2], nnz, h_zeta, p_zeta );
    for ( long int idim=0; idim<ndim; idim++ ) resid[idim] = 0.;
    for ( long int j=0; j<niso*ndim; j++ ) jac[j] = 0.;
    long int inol = 0;
    for ( long int iz=0; iz<nnz; iz++ ) {
      for ( long int iy=0; iy<nny; iy++ ) {
        for ( long int ix=0; ix<nnx; ix++ ) {
          double w    = h_xi[ix]*h_eta[iy]*h_zeta[iz];
          double der[3];
          der[0] = p_xi[ix]*h_eta[iy]*h_zeta[iz];
          der[1] = h_xi[ix]*p_eta[iy]*h_zeta[iz];
          der[2] = h_xi[ix]*h_eta[iy]*p_zeta[iz];
          for ( long int idim=0; idim<ndim; idim++ ) {
            double x = coords[inol*ndim+idim];
            resid[idim] += w*x;
            for ( long int j=0; j<niso; j++ )
              jac[idim*niso+j] += der[j]*x;
          }
          inol++;
        }
      }
    }
    double err = 0.;
    for ( long int idim=0; idim<ndim; idim++ ) {
      resid[idim] -= point[idim];
      err += resid[idim]*resid[idim];
    }
    err = sqrt( err );
    if ( iter==0 || err<err_best ) {
      err_best = err;
      for ( long int j=0; j<3; j++ ) iso_best[j] = iso[j];
    }
    if ( err<=tol ) break;
    // normal equations  (J^T J) delta = -J^T resid
    for ( long int j=0; j<niso*niso; j++ ) normal[j] = 0.;
    for ( long int j=0; j<niso; j++ ) rhside[j] = 0.;
    for ( long int idim=0; idim<ndim; idim++ ) {
      for ( long int j=0; j<niso; j++ ) {
        rhside[j] -= jac[idim*niso+j]*resid[idim];
        for ( long int k=0; k<niso; k++ )
          normal[j*niso+k] += jac[idim*niso+j]*jac[idim*niso+k];
      }
    }
    if ( !mpc_polish_solve( niso, normal, rhside ) ) {
      array_move( weight_in, weight, nnol );
      return;
    }
    double dnorm = 0.;
    for ( long int j=0; j<niso; j++ ) {
      delta[j] = rhside[j];
      dnorm += delta[j]*delta[j];
    }
    if ( sqrt( dnorm )*scale<1.e-16 && err>tol ) {
      array_move( weight_in, weight, nnol );
      return;   // GN stalled short of machine precision: keep point_el
    }
    // damped update: halve the step while it does not reduce the residual
    double step = 1.;
    long int ih = 0;
    for ( ; ih<30; ih++ ) {
      double iso_new[3], weight_new[MPOINT];
      for ( long int j=0; j<3; j++ )
        iso_new[j] = iso[j] + step*delta[j];
      mpc_polish_weights( iso_new, nnx, nny, nnz, weight_new );
      double err_new = 0.;
      for ( long int idim=0; idim<ndim; idim++ ) {
        double x = 0.;
        for ( long int inol2=0; inol2<nnol; inol2++ )
          x += weight_new[inol2]*coords[inol2*ndim+idim];
        double d = x - point[idim];
        err_new += d*d;
      }
      err_new = sqrt( err_new );
      if ( err_new<=err ) {
        for ( long int j=0; j<3; j++ ) iso[j] = iso_new[j];
        break;
      }
      step *= 0.5;
    }
    if ( ih==30 ) {
      array_move( weight_in, weight, nnol );
      return;   // no progress: keep the point_el weights
    }
  }
  if ( err_best>tol ) {
    array_move( weight_in, weight, nnol );
    return;
  }
  mpc_polish_weights( iso_best, nnx, nny, nnz, weight );
}

static long int mpc_element_is_linear( long int name )
{
  return ( name==-BAR2 || name==-TRIA3 || name==-QUAD4 ||
           name==-TET4 || name==-HEX8 );
}

static long int mpc_element_is_quadratic( long int name )
{
  return ( name==-BAR3 || name==-TRIA6 || name==-QUAD9 ||
           name==-TET10 || name==-HEX27 );
}

// mesh fingerprint: a cheap hash over the element topology. Any mesh
// change (refine/delete/split) alters the element set and thus the hash.
static long int mpc_mesh_fingerprint( void )

{
  long int element=0, max_element=0, length=0, nnol=0, inol=0, inod=0,
    max_node=0, hash=0;
  double ddum[1];
  long int el[MNOL+1];

  db_max_index( NODE, max_node, VERSION_NORMAL, GET );
  db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
  hash = 100003 * ( 1 + max_node ) + 1009 * ( 1 + max_element );
  for ( element=0; element<=max_element; element++ ) {
    if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) ) {
      db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
      nnol = length - 1;
      hash ^= element * ( nnol + 7 );
      for ( inol=1; inol<=nnol; inol++ ) {
        inod = el[inol];
        hash ^= inod * ( inol + 3 );
      }
    }
  }
  return hash;

}

// exported for the bounda -veln generator (bounda.cc): the generated
// mpc records must be re-created whenever the mesh changes, with the
// same fingerprint the mpc_linear_quadratic/mpc_element_group
// generators use.
long int mpc_mesh_fingerprint_get( void )

{
  return mpc_mesh_fingerprint();
}

// mpc_linear_quadratic (manual Professional 6.873): generate the tie
// records for the dangling quadratic nodes. Called from mpc_node_apply()
// (every equilibrium iteration, but the generation runs only when the
// mesh fingerprint changed).
static void mpc_linear_quadratic_generate( void )

{
  long int lq=-NO, ldum=0, idum[1], length=0;
  double ddum[1];

  db( MPC_LINEAR_QUADRATIC, 0, &lq, ddum, ldum, VERSION_NORMAL,
    GET_IF_EXISTS );
  if ( lq != -YES ) return;

  long int fingerprint = mpc_mesh_fingerprint();
  long int stored[3] = { -1, -1, 0 };
  long int stored_len = 0;
  db( MPC_LINEAR_QUADRATIC_MESH_FINGERPRINT, 0, stored, ddum, stored_len,
    VERSION_NORMAL, GET_IF_EXISTS );
  if ( stored_len>=3 && stored[0]==fingerprint ) return;

  // delete the previously generated records
  if ( stored_len>=3 && stored[2]>0 ) {
    long int ig=0;
    for ( ig=0; ig<stored[2]; ig++ ) {
      db_delete_index( MPC_NODE_NUMBER, stored[1]+ig, VERSION_NORMAL );
      db_delete_index( MPC_NODE_FACTOR, stored[1]+ig, VERSION_NORMAL );
    }
  }

  long int max_node=0, max_element=0, element=0, ile=0, inol=0, jnol=0,
    inod=0, nnol=0, nnol_lin=0, name=0, name_lin=0, length_lin=0,
    nlinear=0, nmaster=0, i=0, ipuknwn=0, iuknwn=0, start_index=0,
    count=0, km=0, length_put=0;
  long int el[MNOL+1], el_lin[MNOL+1];
  double coord[MDIM], coords_el[MNOL*MDIM], weight[MNOL];

  db_max_index( NODE, max_node, VERSION_NORMAL, GET );
  db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
  if ( max_node<0 || max_element<0 ) {
    // nothing to tie; remember the fingerprint so the (empty) state is
    // not regenerated on every iteration
    stored[0] = fingerprint; stored[1] = -1; stored[2] = 0;
    length_put = 3;
    db( MPC_LINEAR_QUADRATIC_MESH_FINGERPRINT, 0, stored, ddum, length_put,
      VERSION_NORMAL, PUT );
    return;
  }

  // marker: nodes that belong to at least one LINEAR element
  long int *in_linear = get_new_int( 1 + max_node );
  array_set( in_linear, 0, 1 + max_node );
  // list of linear elements (master candidates)
  long int *linear_elems = get_new_int( 1 + max_element );

  for ( element=0; element<=max_element; element++ ) {
    if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) ) {
      db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
      name = el[0];
      if ( mpc_element_is_linear( name ) ) {
        linear_elems[nlinear++] = element;
        for ( inol=1; inol<length; inol++ )
          in_linear[el[inol]] = 1;
      }
    }
  }

  long int dof_label[MUKNWN], dof_principal[MUKNWN];
  db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( DOF_PRINCIPAL, 0, dof_principal, ddum, ldum, VERSION_NORMAL,
    GET_IF_EXISTS );

  // start index of the generated records: after the highest ACTIVE user
  // mpc_node_number record
  db_highest_index( MPC_NODE_NUMBER, start_index, VERSION_NORMAL );
  start_index++;

  long int vals[DATA_ITEM_SIZE];

  for ( element=0; element<=max_element; element++ ) {
    if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) ) {
      db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
      name = el[0];
      if ( !mpc_element_is_quadratic( name ) ) continue;
      nnol = length - 1;
      for ( inol=1; inol<=nnol; inol++ ) {
        inod = el[inol];
        // only dangling nodes: not attached to any linear element
        if ( in_linear[inod] ) continue;
        db( NODE_START_REFINED, inod, idum, coord, ldum, VERSION_NORMAL,
          GET );
        // find the linear element that contains the node
        for ( ile=0; ile<nlinear; ile++ ) {
          db( ELEMENT, linear_elems[ile], el_lin, ddum, length_lin,
            VERSION_NORMAL, GET );
          name_lin = el_lin[0];
          nnol_lin = length_lin - 1;
          for ( jnol=0; jnol<nnol_lin; jnol++ ) {
            db( NODE_START_REFINED, el_lin[1+jnol], idum,
              &coords_el[jnol*ndim], ldum, VERSION_NORMAL, GET );
          }
          if ( point_el( coord, coords_el, weight, name_lin, nnol_lin,
              MPC_EPS_ISO ) ) {
            // machine-precision factors (DIAG §21): point_el's
            // finite-difference descent stops ~2.4e-8 short
            mpc_weights_polish( coord, coords_el, weight, name_lin,
              nnol_lin );
            // tie the node to the linear element: one mpc record per
            // principal dof, master weights = the linear shape functions
            // evaluated at the node position
            nmaster = 0;
            for ( jnol=0; jnol<nnol_lin; jnol++ ) {
              if ( scalar_dabs( weight[jnol] ) > MPC_WEIGHT_ZERO ) {
                el_lin[1+nmaster] = el_lin[1+jnol];
                weight[nmaster] = weight[jnol];
                nmaster++;
              }
            }
            if ( nmaster==0 ) continue;
            for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
              iuknwn = ipuknwn*nder;
              if ( dof_principal[iuknwn]>=0 ) {
                // record: [node_0 dof_0 node_1 dof_0 node_2 dof_0 ...]
                i = 0;
                vals[i++] = inod;
                vals[i++] = dof_label[iuknwn];
                for ( km=0; km<nmaster; km++ ) {
                  vals[i++] = el_lin[1+km];
                  vals[i++] = dof_label[iuknwn];
                }
                length_put = i;
                db( MPC_NODE_NUMBER, start_index+count, vals, ddum,
                  length_put, VERSION_NORMAL, PUT );
                length_put = nmaster;
                db( MPC_NODE_FACTOR, start_index+count, idum, weight,
                  length_put, VERSION_NORMAL, PUT );
                count++;
              }
            }
            break;
          }
        }
      }
    }
  }

  delete[] in_linear;
  delete[] linear_elems;

  stored[0] = fingerprint; stored[1] = start_index; stored[2] = count;
  length_put = 3;
  db( MPC_LINEAR_QUADRATIC_MESH_FINGERPRINT, 0, stored, ddum, length_put,
    VERSION_NORMAL, PUT );

}

// ---------------------------------------------------------------------------
// mpc_element_group (manual Professional 6.858): each node of an element of
// group element_group_0 that is also located inside an element of group
// element_group_1 is tied to that element with multi point constraints
// consistent with the shape functions at the node's isoparametric
// location. mpc_element_group_dof (6.862) selects the dofs that are set
// equal (default: all principal dofs); mpc_element_group_geometry (6.864)
// restricts the nodes of element_group_0 to those lying on a geometry.
// The generated mpc_node_number/mpc_node_factor records are the same ones
// the Professional writes to the .dbs, and are re-created when the mesh
// changes (fingerprint bookkeeping, mirroring mpc_linear_quadratic).
//
// Master elements that already contain the node in their node list are
// skipped (the tie would reduce to the node itself and freeze its
// equation); the search continues with the other master elements
// (measured on the Professional 02-08-2026: a node shared with one
// element_group_1 element and located at a node of another one is tied to
// the other one).
//
// Options (F5 T4 U6, manual 6.256 and 6.859-6.865):
// - mpc_element_group_always -no: nodes that are a member of the node
//   list of any element of element_group_1 (the groups are connected
//   there) get no mpc. Default -yes: the behaviour above. The
//   Professional 02-08-2026 ignores -no (DIVERGENCES
//   MPC-ELEMENT-GROUP-ALWAYS).
// - mpc_element_group_closest -yes: a node that lies in no element of
//   element_group_1 is tied to the closest one: the master with the
//   smallest distance between the node and its closest location, with
//   the shape functions at that location (mpc_closest_weights; measured
//   on the Professional 02-08-2026 with quad4 / hex8 / tria3 masters).
// - mpc_element_group_coord_geometry -yes: only nodes whose last
//   coordinate (y in 2D, z in 3D) is not below the lowest node of
//   element_group_1 are tied (manual; the Professional takes only the
//   element_group_1 nodes on mpc_element_group_geometry, DIVERGENCES
//   MPC-ELEMENT-GROUP-COORD-GEOMETRY).
// - mpc_element_group_eps_iso: the tolerance on the isoparametric
//   coordinates (default MPC_EPS_ISO = 1.e-4, manual).
// - mpc_element_group_keep -yes: the records generated for that index are
//   kept for the remainder of the calculation (not re-created on a mesh
//   change or by control_mpc_element_group). They live outside the
//   regenerated range; their range is stored in the internal record
//   mpc_element_group_kept.
// - control_mpc_element_group -yes: the records are re-created at every
//   timestep of that control index (the fingerprint then also carries the
//   current time); default -no: only when the mesh changed.
// ---------------------------------------------------------------------------

// tensor-product dimensions of a master for the closest point search
// (bars included); returns 0 for an unsupported element
static long int mpc_closest_dims( long int name, long int *nnx, long int *nny,
  long int *nnz )

{
  *nny = 1; *nnz = 1;
  if      ( name==-BAR2 ) *nnx = 2;
  else if ( name==-BAR3 ) *nnx = 3;
  else if ( name==-BAR4 ) *nnx = 4;
  else mpc_polish_dims( name, nnx, nny, nnz );
  return ( *nnx>1 );
}

// closest point of a tensor-product master: every face of the
// isoparametric box (each coordinate free, -1 or +1) is solved by
// Gauss-Newton over its free coordinates; the face solution with free
// coordinates inside [-1,1] and the smallest distance wins (the
// Euclidean closest point of the element; on rectangular quad4 / hex8
// masters equal to the Professional 02-08-2026, which approximates it on
// distorted ones)
static long int mpc_closest_tensor( double point[], double coords[],
  long int name, long int nnol, double weight[], double *dist )

{
  long int nnx=0, nny=0, nnz=0;
  if ( !mpc_closest_dims( name, &nnx, &nny, &nnz ) ) return 0;
  long int nn[3] = { nnx, nny, nnz };
  long int niso = ( nnx>1 ) + ( nny>1 ) + ( nnz>1 );
  long int ncombo = 1;
  for ( long int j=0; j<niso; j++ ) ncombo *= 3;

  double h[3][MPOINT], p[3][MPOINT], iso[3], iso_best[3], x[MDIM],
    jac[MDIM*3], normal[9], rhside[3];
  double best = -1.;
  for ( long int ic=0; ic<ncombo; ic++ ) {
    long int state[3] = { 0, 0, 0 }, free_iso[3], nfree = 0, c = ic;
    for ( long int j=0; j<3; j++ ) {
      iso[j] = 0.;
      if ( j<niso ) { state[j] = c % 3; c /= 3; }
      if ( state[j]==1 ) iso[j] = -1.;
      else if ( state[j]==2 ) iso[j] = 1.;
      else if ( j<niso ) free_iso[nfree++] = j;
    }
    for ( long int iter=0; iter<=50; iter++ ) {
      for ( long int j=0; j<3; j++ )
        interpolation_polynomial( iso[j], nn[j], h[j], p[j] );
      array_set( x, 0., ndim );
      array_set( jac, 0., MDIM*3 );
      long int inol = 0;
      for ( long int iz=0; iz<nnz; iz++ ) {
        for ( long int iy=0; iy<nny; iy++ ) {
          for ( long int ix=0; ix<nnx; ix++ ) {
            double w = h[0][ix]*h[1][iy]*h[2][iz];
            double der[3] = { p[0][ix]*h[1][iy]*h[2][iz],
              h[0][ix]*p[1][iy]*h[2][iz], h[0][ix]*h[1][iy]*p[2][iz] };
            for ( long int idim=0; idim<ndim; idim++ ) {
              x[idim] += w*coords[inol*ndim+idim];
              for ( long int k=0; k<nfree; k++ )
                jac[idim*3+k] += der[free_iso[k]]*coords[inol*ndim+idim];
            }
            inol++;
          }
        }
      }
      if ( nfree==0 || iter==50 ) break;
      for ( long int k=0; k<9; k++ ) normal[k] = 0.;
      for ( long int k=0; k<3; k++ ) rhside[k] = 0.;
      for ( long int idim=0; idim<ndim; idim++ ) {
        for ( long int k=0; k<nfree; k++ ) {
          rhside[k] -= jac[idim*3+k]*( x[idim]-point[idim] );
          for ( long int l=0; l<nfree; l++ )
            normal[k*nfree+l] += jac[idim*3+k]*jac[idim*3+l];
        }
      }
      if ( !mpc_polish_solve( nfree, normal, rhside ) ) break;
      double dnorm = 0.;
      for ( long int k=0; k<nfree; k++ ) {
        iso[free_iso[k]] += rhside[k];
        dnorm += rhside[k]*rhside[k];
      }
      if ( sqrt( dnorm )<1.e-15 ) break;
    }
    long int inside = 1;
    for ( long int k=0; k<nfree; k++ )
      if ( scalar_dabs( iso[free_iso[k]] )>1.+1.e-10 ) inside = 0;
    if ( !inside ) continue;
    for ( long int k=0; k<nfree; k++ ) {
      if ( iso[free_iso[k]]>1. ) iso[free_iso[k]] = 1.;
      if ( iso[free_iso[k]]<-1. ) iso[free_iso[k]] = -1.;
    }
    mpc_polish_weights( iso, nnx, nny, nnz, weight );
    array_set( x, 0., ndim );
    for ( long int inol=0; inol<nnol; inol++ )
      for ( long int idim=0; idim<ndim; idim++ )
        x[idim] += weight[inol]*coords[inol*ndim+idim];
    double d = 0.;
    for ( long int idim=0; idim<ndim; idim++ )
      d += ( x[idim]-point[idim] )*( x[idim]-point[idim] );
    d = sqrt( d );
    if ( best<0. || d<best ) {
      best = d;
      for ( long int j=0; j<3; j++ ) iso_best[j] = iso[j];
    }
  }
  if ( best<0. ) return 0;
  mpc_polish_weights( iso_best, nnx, nny, nnz, weight );
  *dist = best;
  return 1;
}

// closest point of a simplex master (tria3/tria6/tet4/tet10) with the
// corner nodes: the barycentric coordinates of the point (first corner
// dependent) have the other corners' coordinates clamped to [0,1]; a
// negative first-corner coordinate is set to zero and the coordinates
// are normalized (measured on the Professional 02-08-2026 with tria3
// masters: (1.1,0.3) beside the tria3 (0,0) (1,0) (1,1) gets 8/11 3/11,
// (2,0.3) gets 10/13 3/13; its tet4 masters ignore the fourth corner,
// which the GNU does not copy: DIVERGENCES MPC-ELEMENT-GROUP-CLOSEST).
// The weights are the element shape functions at that location (the
// point_el formulas).
static long int mpc_closest_simplex( double point[], double coords[],
  long int name, double weight[], double *dist )

{
  long int corner[4], ncorner = 0;
  if      ( name==-TRIA3 ) { corner[0]=0; corner[1]=1; corner[2]=2; ncorner=3; }
  else if ( name==-TRIA6 ) { corner[0]=0; corner[1]=2; corner[2]=5; ncorner=3; }
  else if ( name==-TET4 ) {
    corner[0]=0; corner[1]=1; corner[2]=2; corner[3]=3; ncorner=4; }
  else if ( name==-TET10 ) {
    corner[0]=0; corner[1]=2; corner[2]=5; corner[3]=9; ncorner=4; }
  else return 0;
  if ( ncorner-1!=ndim ) return 0;

  // barycentric coordinates: point = v0 + sum_k L_k (v_k - v0)
  double *v0 = &coords[corner[0]*ndim];
  double a[9], b[3], L[4] = { 0., 0., 0., 0. };
  long int n = ncorner - 1;
  for ( long int idim=0; idim<n; idim++ ) {
    b[idim] = point[idim] - v0[idim];
    for ( long int k=0; k<n; k++ )
      a[idim*n+k] = coords[corner[k+1]*ndim+idim] - v0[idim];
  }
  if ( !mpc_polish_solve( n, a, b ) ) return 0;
  double sum = 0.;
  for ( long int k=0; k<n; k++ ) {
    L[k+1] = b[k];
    if ( L[k+1]<0. ) L[k+1] = 0.;
    if ( L[k+1]>1. ) L[k+1] = 1.;
    sum += L[k+1];
  }
  L[0] = 1. - sum;
  if ( L[0]<0. ) {
    L[0] = 0.;
    for ( long int k=1; k<ncorner; k++ ) L[k] /= sum;
  }
  double d = 0.;
  for ( long int idim=0; idim<ndim; idim++ ) {
    double xk = 0.;
    for ( long int k=0; k<ncorner; k++ )
      xk += L[k]*coords[corner[k]*ndim+idim];
    d += ( xk-point[idim] )*( xk-point[idim] );
  }
  double L1 = L[0], L2 = L[1], L3 = L[2], L4 = L[3];
  if      ( name==-TRIA3 ) {
    weight[0] = L1; weight[1] = L2; weight[2] = L3;
  }
  else if ( name==-TRIA6 ) {
    weight[0] = (2.*L1-1.)*L1;
    weight[1] = 4.*L1*L2;
    weight[2] = (2.*L2-1.)*L2;
    weight[3] = 4.*L3*L1;
    weight[4] = 4.*L2*L3;
    weight[5] = (2.*L3-1.)*L3;
  }
  else if ( name==-TET4 ) {
    weight[0] = L1; weight[1] = L2; weight[2] = L3; weight[3] = L4;
  }
  else {
    weight[0] = (2.*L1-1.)*L1;
    weight[1] = 4.*L1*L2;
    weight[2] = (2.*L2-1.)*L2;
    weight[3] = 4.*L3*L1;
    weight[4] = 4.*L2*L3;
    weight[5] = (2.*L3-1.)*L3;
    weight[6] = 4.*L1*L4;
    weight[7] = 4.*L2*L4;
    weight[8] = 4.*L3*L4;
    weight[9] = (2.*L4-1.)*L4;
  }
  *dist = sqrt( d );
  return 1;
}

// closest point of a master element (mpc_element_group_closest): the
// weights are the shape functions at the point of the element closest to
// point[]; returns 0 for an unsupported element (prisms, trusses, ...)
static long int mpc_closest_weights( double point[], double coords[],
  long int name, long int nnol, double weight[], double *dist )

{
  if ( name==-TRIA3 || name==-TRIA6 || name==-TET4 || name==-TET10 )
    return mpc_closest_simplex( point, coords, name, weight, dist );
  return mpc_closest_tensor( point, coords, name, nnol, weight, dist );
}

// generate the records of mpc_element_group index ieg at
// start_index+count..., incrementing count
static void mpc_element_group_index( long int ieg, long int start_index,
  long int &count, long int max_element, long int max_node,
  long int dof_label[], long int dof_principal[] )

{
  long int ldum=0, idum[1], length=0, i=0;
  double ddum[MDIM]; // geometry() writes normal/projection[MDIM]
  long int vals[DATA_ITEM_SIZE];

  long int *eg_vals = db_int( MPC_ELEMENT_GROUP, ieg, VERSION_NORMAL );
  long int eg_length = db_len( MPC_ELEMENT_GROUP, ieg, VERSION_NORMAL );
  if ( eg_length<2 ) db_error( MPC_ELEMENT_GROUP, ieg );
  long int gr0 = eg_vals[0];

  // master element groups: eg_vals[1..] or -all (all other groups)
  long int master_groups[DATA_ITEM_SIZE];
  long int nmaster_groups = 0, all_groups = 0;
  for ( i=1; i<eg_length; i++ ) {
    if ( eg_vals[i]==-ALL ) all_groups = 1;
    else master_groups[nmaster_groups++] = eg_vals[i];
  }

  // selected dofs (negative dof labels); empty -> all principal dofs
  long int dof_sel[MUKNWN];
  long int ndof_sel = 0;
  if ( db_active_index( MPC_ELEMENT_GROUP_DOF, ieg, VERSION_NORMAL ) ) {
    long int *dof_vals = db_int( MPC_ELEMENT_GROUP_DOF, ieg,
      VERSION_NORMAL );
    long int dof_length = db_len( MPC_ELEMENT_GROUP_DOF, ieg,
      VERSION_NORMAL );
    for ( i=0; i<dof_length; i++ ) dof_sel[ndof_sel++] = dof_vals[i];
  }

  // restricting geometry pairs (geometry entity + index per pair)
  long int geom_vals[DATA_ITEM_SIZE];
  long int geom_length = 0;
  if ( db_active_index( MPC_ELEMENT_GROUP_GEOMETRY, ieg,
      VERSION_NORMAL ) ) {
    long int *gvals = db_int( MPC_ELEMENT_GROUP_GEOMETRY, ieg,
      VERSION_NORMAL );
    geom_length = db_len( MPC_ELEMENT_GROUP_GEOMETRY, ieg,
      VERSION_NORMAL );
    for ( i=0; i<geom_length; i++ ) geom_vals[i] = gvals[i];
  }

  // options (manual 6.859-6.863)
  long int always = -YES, closest = -NO, coord_geometry = -NO;
  db( MPC_ELEMENT_GROUP_ALWAYS, ieg, &always, ddum, ldum, VERSION_NORMAL,
    GET_IF_EXISTS );
  db( MPC_ELEMENT_GROUP_CLOSEST, ieg, &closest, ddum, ldum, VERSION_NORMAL,
    GET_IF_EXISTS );
  db( MPC_ELEMENT_GROUP_COORD_GEOMETRY, ieg, &coord_geometry, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  double eps_iso = MPC_EPS_ISO;
  db( MPC_ELEMENT_GROUP_EPS_ISO, ieg, idum, &eps_iso, ldum, VERSION_NORMAL,
    GET_IF_EXISTS );
  if ( eps_iso<0. ) db_error( MPC_ELEMENT_GROUP_EPS_ISO, ieg );

  long int element=0, inod=0, inol=0, nnol=0, jelem=0, jnol=0,
    nnol_m=0, name_m=0, length_m=0, nmaster=0, k=0,
    km=0, ipuknwn=0, iuknwn=0, length_put=0, iuse=0, found=0,
    in_geom=0;
  long int el[MNOL+1], el_m[MNOL+1], el_grp=0;
  double coord[MDIM], coords_m[MNOL*MDIM], weight[MNOL], rdum=0.;
  long int geometry_entity[2];

  // the element_group_1 nodes (mpc_element_group_always -no) and their
  // lowest last coordinate (mpc_element_group_coord_geometry -yes)
  long int *in_master = NULL;
  double lowest = DBL_MAX;
  if ( ( always==-NO || coord_geometry==-YES ) && max_node>=0 ) {
    in_master = get_new_int( 1+max_node );
    array_set( in_master, 0, 1+max_node );
    for ( jelem=0; jelem<=max_element; jelem++ ) {
      if ( !db_active_index( ELEMENT, jelem, VERSION_NORMAL ) ) continue;
      el_grp = 0;
      db( ELEMENT_GROUP, jelem, &el_grp, ddum, length, VERSION_NORMAL,
        GET_IF_EXISTS );
      iuse = 0;
      if ( all_groups ) { if ( el_grp!=gr0 ) iuse = 1; }
      else {
        for ( k=0; k<nmaster_groups; k++ )
          if ( master_groups[k]==el_grp ) iuse = 1;
      }
      if ( !iuse ) continue;
      db( ELEMENT, jelem, el_m, ddum, length_m, VERSION_NORMAL, GET );
      for ( jnol=1; jnol<length_m; jnol++ ) {
        if ( in_master[el_m[jnol]] ) continue;
        in_master[el_m[jnol]] = 1;
        db( NODE_START_REFINED, el_m[jnol], idum, coord, ldum,
          VERSION_NORMAL, GET );
        if ( coord[ndim-1]<lowest ) lowest = coord[ndim-1];
      }
    }
  }

  // visited marker: each node is tied at most once per record
  long int *visited = NULL;
  if ( max_node>=0 ) {
    visited = get_new_int( 1+max_node );
    array_set( visited, 0, 1+max_node );
  }

  for ( element=0; element<=max_element; element++ ) {
    if ( !db_active_index( ELEMENT, element, VERSION_NORMAL ) ) continue;
    el_grp = 0;
    db( ELEMENT_GROUP, element, &el_grp, ddum, length, VERSION_NORMAL,
      GET_IF_EXISTS );
    if ( el_grp!=gr0 ) continue;
    db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
    nnol = length - 1;
    for ( inol=1; inol<=nnol; inol++ ) {
      inod = el[inol];
      if ( visited && visited[inod] ) continue;
      // mpc_element_group_always -no: the node is connected to
      // element_group_1
      if ( always==-NO && in_master && in_master[inod] ) continue;
      // geometry restriction
      if ( geom_length>0 ) {
        in_geom = 0;
        for ( k=0; k+1<geom_length; k+=2 ) {
          geometry_entity[0] = geom_vals[k];
          geometry_entity[1] = geom_vals[k+1];
          geometry( inod, ddum, geometry_entity, in_geom, rdum, ddum,
            rdum, ddum, NODE_START_REFINED, PROJECT_EXACT,
            VERSION_NORMAL );
          if ( in_geom ) break;
        }
        if ( !in_geom ) continue;
      }
      db( NODE_START_REFINED, inod, idum, coord, ldum, VERSION_NORMAL,
        GET );
      // mpc_element_group_coord_geometry -yes: not below the lowest
      // element_group_1 node
      if ( coord_geometry==-YES && coord[ndim-1]<lowest ) continue;
      // find a master element containing the node (first match); with
      // mpc_element_group_closest -yes, else the closest master element
      found = 0;
      double dist_best = -1.;
      long int name_best = 0, nnol_best = 0, el_best[MNOL+1];
      double weight_best[MNOL];
      for ( long int pass=0; pass<2 && !found; pass++ ) {
        if ( pass==1 && closest!=-YES ) break;
        for ( jelem=0; jelem<=max_element && !found; jelem++ ) {
          if ( !db_active_index( ELEMENT, jelem, VERSION_NORMAL ) )
            continue;
          el_grp = 0;
          db( ELEMENT_GROUP, jelem, &el_grp, ddum, length, VERSION_NORMAL,
            GET_IF_EXISTS );
          iuse = 0;
          if ( all_groups ) { if ( el_grp!=gr0 ) iuse = 1; }
          else {
            for ( k=0; k<nmaster_groups; k++ )
              if ( master_groups[k]==el_grp ) iuse = 1;
          }
          if ( !iuse ) continue;
          db( ELEMENT, jelem, el_m, ddum, length_m, VERSION_NORMAL, GET );
          name_m = el_m[0]; nnol_m = length_m - 1;
          // skip master elements that already contain the node as a
          // member (the meshes are connected there; a tie would reduce
          // to the node itself)
          for ( jnol=1; jnol<=nnol_m; jnol++ )
            if ( el_m[jnol]==inod ) iuse = 0;
          if ( !iuse ) continue;
          for ( jnol=0; jnol<nnol_m; jnol++ ) {
            db( NODE_START_REFINED, el_m[1+jnol], idum,
              &coords_m[jnol*ndim], ldum, VERSION_NORMAL, GET );
          }
          if ( pass==0 ) {
            if ( point_el( coord, coords_m, weight, name_m, nnol_m,
                eps_iso ) ) {
              // machine-precision factors (DIAG §21): point_el's
              // finite-difference descent stops ~2.4e-8 short
              mpc_weights_polish( coord, coords_m, weight, name_m,
                nnol_m );
              found = 1;
            }
          }
          else {
            // lower bound of the distance: the bounding box distance
            double lower = 0.;
            for ( long int idim=0; idim<ndim; idim++ ) {
              double lo = coords_m[idim], hi = coords_m[idim];
              for ( jnol=1; jnol<nnol_m; jnol++ ) {
                lo = scalar_dmin( lo, coords_m[jnol*ndim+idim] );
                hi = scalar_dmax( hi, coords_m[jnol*ndim+idim] );
              }
              double out = 0.;
              if ( coord[idim]<lo ) out = lo - coord[idim];
              if ( coord[idim]>hi ) out = coord[idim] - hi;
              lower += out*out;
            }
            lower = sqrt( lower );
            if ( dist_best>=0. && lower>=dist_best ) continue;
            double dist = 0.;
            if ( mpc_closest_weights( coord, coords_m, name_m, nnol_m,
                weight, &dist ) && ( dist_best<0. || dist<dist_best ) ) {
              dist_best = dist;
              name_best = name_m; nnol_best = nnol_m;
              for ( jnol=0; jnol<=nnol_m; jnol++ ) el_best[jnol] = el_m[jnol];
              for ( jnol=0; jnol<nnol_m; jnol++ )
                weight_best[jnol] = weight[jnol];
            }
          }
        }
        if ( pass==1 && dist_best>=0. ) {
          name_m = name_best; nnol_m = nnol_best;
          for ( jnol=0; jnol<=nnol_m; jnol++ ) el_m[jnol] = el_best[jnol];
          for ( jnol=0; jnol<nnol_m; jnol++ ) weight[jnol] = weight_best[jnol];
          found = 1;
        }
      }
      if ( !found ) continue;
      // masters: the master element nodes with non-zero shape function
      // weight at the node position
      nmaster = 0;
      for ( jnol=0; jnol<nnol_m; jnol++ ) {
        if ( scalar_dabs( weight[jnol] ) > MPC_WEIGHT_ZERO ) {
          el_m[1+nmaster] = el_m[1+jnol];
          weight[nmaster] = weight[jnol];
          nmaster++;
        }
      }
      if ( nmaster==0 ) continue;
      // generate one mpc record per selected dof
      for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
        iuknwn = ipuknwn*nder;
        if ( dof_principal[iuknwn]<0 ) continue;
        if ( ndof_sel>0 ) {
          long int iuse_dof = 0;
          for ( k=0; k<ndof_sel; k++ )
            if ( dof_sel[k]==dof_label[iuknwn] ) iuse_dof = 1;
          if ( !iuse_dof ) continue;
        }
        // record: [node_0 dof_0 node_1 dof_0 node_2 dof_0 ...]
        long int iv = 0;
        vals[iv++] = inod;
        vals[iv++] = dof_label[iuknwn];
        for ( km=0; km<nmaster; km++ ) {
          vals[iv++] = el_m[1+km];
          vals[iv++] = dof_label[iuknwn];
        }
        length_put = iv;
        db( MPC_NODE_NUMBER, start_index+count, vals, ddum,
          length_put, VERSION_NORMAL, PUT );
        length_put = nmaster;
        db( MPC_NODE_FACTOR, start_index+count, idum, weight,
          length_put, VERSION_NORMAL, PUT );
        count++;
      }
      if ( visited ) visited[inod] = 1;
    }
  }
  delete[] visited;
  delete[] in_master;
}

static void mpc_element_group_generate( void )

{
  long int ldum=0, idum[1], length=0, max_meg=0, ieg=0, i=0;
  double ddum[1];

  db_max_index( MPC_ELEMENT_GROUP, max_meg, VERSION_NORMAL, GET );
  if ( max_meg<0 ) return;

  long int fingerprint = mpc_mesh_fingerprint();
  // control_mpc_element_group -yes (manual 6.256): re-create at every
  // timestep of the current control index; the fingerprint then also
  // carries the current time (constant within the timestep)
  long int icontrol = 0, control_eg = -NO;
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_MPC_ELEMENT_GROUP, icontrol, &control_eg, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  if ( control_eg==-YES ) {
    double time_current = 0.;
    db( TIME_CURRENT, 0, idum, &time_current, ldum, VERSION_NORMAL,
      GET_IF_EXISTS );
    unsigned long long bits = 0;
    memcpy( &bits, &time_current, sizeof(bits) );
    fingerprint ^= (long int) ( ( bits ^ ( bits>>29 ) ) | 1ULL );
  }
  long int stored[3] = { -1, -1, 0 };
  long int stored_len = 0;
  db( MPC_ELEMENT_GROUP_MESH_FINGERPRINT, 0, stored, ddum, stored_len,
    VERSION_NORMAL, GET_IF_EXISTS );
  if ( stored_len>=3 && stored[0]==fingerprint ) return;

  // delete the previously generated records
  if ( stored_len>=3 && stored[2]>0 ) {
    for ( i=0; i<stored[2]; i++ ) {
      db_delete_index( MPC_NODE_NUMBER, stored[1]+i, VERSION_NORMAL );
      db_delete_index( MPC_NODE_FACTOR, stored[1]+i, VERSION_NORMAL );
    }
  }

  // records kept by mpc_element_group_keep -yes: drop the ones whose
  // nodes no longer exist after a mesh change
  for ( ieg=0; ieg<=max_meg; ieg++ ) {
    if ( !db_active_index( MPC_ELEMENT_GROUP_KEPT, ieg, VERSION_NORMAL ) )
      continue;
    long int kept[2] = { 0, 0 };
    db( MPC_ELEMENT_GROUP_KEPT, ieg, kept, ddum, ldum, VERSION_NORMAL,
      GET );
    for ( i=kept[0]; i<kept[0]+kept[1]; i++ ) {
      if ( !db_active_index( MPC_NODE_NUMBER, i, VERSION_NORMAL ) )
        continue;
      long int *rec = db_int( MPC_NODE_NUMBER, i, VERSION_NORMAL );
      long int rec_length = db_len( MPC_NODE_NUMBER, i, VERSION_NORMAL );
      long int alive = 1;
      for ( long int j=0; j<rec_length; j+=2 )
        if ( !db_active_index( NODE, rec[j], VERSION_NORMAL ) ) alive = 0;
      if ( !alive ) {
        db_delete_index( MPC_NODE_NUMBER, i, VERSION_NORMAL );
        db_delete_index( MPC_NODE_FACTOR, i, VERSION_NORMAL );
      }
    }
  }

  long int dof_label[MUKNWN], dof_principal[MUKNWN];
  db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( DOF_PRINCIPAL, 0, dof_principal, ddum, ldum, VERSION_NORMAL,
    GET_IF_EXISTS );

  // start index of the generated records: after the highest ACTIVE
  // mpc_node_number record (user records, any mpc_linear_quadratic
  // generated records and the kept records)
  long int start_index = 0;
  db_highest_index( MPC_NODE_NUMBER, start_index, VERSION_NORMAL );
  start_index++;
  long int count = 0;

  long int max_element=0, max_node=0;
  db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
  db_max_index( NODE, max_node, VERSION_NORMAL, GET );

  // first the indices with mpc_element_group_keep -yes that have no kept
  // records yet (their records stay outside the regenerated range), then
  // the other indices
  for ( long int pass=0; pass<2; pass++ ) {
    for ( ieg=0; ieg<=max_meg; ieg++ ) {
      if ( !db_active_index( MPC_ELEMENT_GROUP, ieg, VERSION_NORMAL ) )
        continue;
      long int keep = -NO;
      db( MPC_ELEMENT_GROUP_KEEP, ieg, &keep, ddum, ldum, VERSION_NORMAL,
        GET_IF_EXISTS );
      if ( keep==-YES ) {
        if ( pass==1 ||
             db_active_index( MPC_ELEMENT_GROUP_KEPT, ieg, VERSION_NORMAL ) )
          continue;
        long int count_before = count;
        mpc_element_group_index( ieg, start_index, count, max_element,
          max_node, dof_label, dof_principal );
        if ( count>count_before ) {
          long int kept[2] = { start_index+count_before, count-count_before };
          length = 2;
          db( MPC_ELEMENT_GROUP_KEPT, ieg, kept, ddum, length,
            VERSION_NORMAL, PUT );
        }
      }
      else if ( pass==1 )
        mpc_element_group_index( ieg, start_index, count, max_element,
          max_node, dof_label, dof_principal );
    }
    if ( pass==0 ) {
      // the regenerated range starts after the kept records
      start_index += count;
      count = 0;
    }
  }

  stored[0] = fingerprint; stored[1] = start_index; stored[2] = count;
  length = 3;
  db( MPC_ELEMENT_GROUP_MESH_FINGERPRINT, 0, stored, ddum, length,
    VERSION_NORMAL, PUT );

}

// consume the mpc_node_number / mpc_node_factor records. Called from
// bounda() AFTER the bounda records, once per equilibrium iteration: the
// slave dofs are bounded and get the value sum(factor*master dof) from
// the master values of the current iteration (VERSION_NEW).
void mpc_node_apply( void )

{
  long int ldum=0, idum[1], max_mpc=0, impc=0, max_node=0;
  double ddum[1];

  // mpc_apply (manual 6.859, no index) / control_mpc_apply (6.255,
  // current control index): switch the whole mpc machinery off
  long int mpc_apply = -YES, icontrol = 0;
  db( MPC_APPLY, 0, &mpc_apply, ddum, ldum, VERSION_NORMAL,
    GET_IF_EXISTS );
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_MPC_APPLY, icontrol, &mpc_apply, ddum, ldum, VERSION_NORMAL,
    GET_IF_EXISTS );
  if ( mpc_apply==-NO ) {
    mpc_tie_map_clear();
    return;
  }

  // mpc_linear_quadratic and mpc_element_group generation (mesh-change
  // detection inside)
  mpc_linear_quadratic_generate();
  mpc_element_group_generate();

  db_max_index( MPC_NODE_NUMBER, max_mpc, VERSION_NORMAL, GET );
  if ( max_mpc<0 ) {
    mpc_tie_map_clear();
    return;
  }

  // index range of the records generated by mpc_linear_quadratic (stored
  // by the generator in the fingerprint record): only those ties take the
  // energy-consistent elimination inside the solve
  {
    long int elim_lo = -1, elim_hi = -1;
    long int stored_mpc[3] = { -1, -1, 0 }, stored_mpc_len = 0;
    db( MPC_LINEAR_QUADRATIC_MESH_FINGERPRINT, 0, stored_mpc, ddum,
      stored_mpc_len, VERSION_NORMAL, GET_IF_EXISTS );
    if ( stored_mpc_len>=3 && stored_mpc[1]>=0 && stored_mpc[2]>0 ) {
      elim_lo = stored_mpc[1];
      elim_hi = stored_mpc[1] + stored_mpc[2] - 1;
    }
    mpc_tie_map_clear();
    mpc_tie_elim_lo = elim_lo;
    mpc_tie_elim_hi = elim_hi;
  }

  // NOTE (dev/fmpc, measured): the mpc_element_group GENERATED records are
  // NOT eligible for the elimination even though the Professional condenses
  // them too (measured on force16: the Pro field is the condensed one,
  // 6.0342857 = 1056/175, vs the GNU value-constrained 3.93667; on mpc7 the
  // Professional field is sigma = 1 uniform = the condensed state). The
  // staggered map of the GNU carries an intrinsic ~1e-8/iteration
  // inconsistency with the elimination active (pre-existing DIAG-SOLVE-MIXTO
  // §16 behaviour: mpc3 reaches 3e-9 after the 2 equilibrium iterations);
  // the mpc7 corpus target sits at +-1e-8 with the same 2-iteration cap, so
  // activating the elimination there misses the target by 1.9e-9. The
  // closure of that map inconsistency lives in the solve/element interplay,
  // not in the tie bookkeeping: the mpc_element_group ties keep the
  // value-constrained semantics (mpc7/mpc6 byte-identical), and force16
  // stays a documented diagnosis.
  long int mpc_eg_lo = -1, mpc_eg_hi = -1;
  {
    long int stored_eg[3] = { -1, -1, 0 }, stored_eg_len = 0;
    db( MPC_ELEMENT_GROUP_MESH_FINGERPRINT, 0, stored_eg, ddum,
      stored_eg_len, VERSION_NORMAL, GET_IF_EXISTS );
    if ( stored_eg_len>=3 && stored_eg[1]>=0 && stored_eg[2]>0 ) {
      mpc_eg_lo = stored_eg[1];
      mpc_eg_hi = stored_eg[1] + stored_eg[2] - 1;
    }
  }

  long int dof_label[MUKNWN];
  db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );

  // post-solve re-sync of the same iteration (top.cc calls
  // mpc_node_apply again after the solve): the bounded flags of the
  // slaves the bounda() pass applied belong to the mpc machinery, not
  // to a direct bounda record - clear them so the record loop below
  // re-derives which slaves a direct prescription owns (those keep
  // their direct value, measured on the Professional binary) and
  // refreshes the free ones with the freshly solved masters. The list
  // is re-collected below on every apply.
  if ( mpc_applied_n>0 ) {
    for ( long int ia=0; ia<mpc_applied_n; ia++ ) {
      long int ig = mpc_applied[ia];
      long int ig_nod = ig / npuknwn, ig_puk = ig % npuknwn;
      if ( db_active_index( NODE, ig_nod, VERSION_NORMAL ) ) {
        long int *nb = db_int( NODE_BOUNDED, ig_nod, VERSION_NORMAL );
        nb[ig_puk] = 0;
      }
    }
    mpc_applied_n = 0;
  }

  for ( impc=0; impc<=max_mpc; impc++ ) {
    if ( db_active_index( MPC_NODE_NUMBER, impc, VERSION_NORMAL ) ) {
      long int length_mpc = db_len( MPC_NODE_NUMBER, impc, VERSION_NORMAL );
      long int *mpc = db_int( MPC_NODE_NUMBER, impc, VERSION_NORMAL );
      if ( length_mpc<2 ) db_error( MPC_NODE_NUMBER, impc );

      // slave: node_0 dof_0
      long int node_0 = mpc[0];
      long int iuknwn_0 = -1;
      array_member( dof_label, mpc[1], nuknwn, iuknwn_0 );
      if ( iuknwn_0<0 ) db_error( MPC_NODE_NUMBER, impc );
      long int ipuknwn_0 = iuknwn_0 / nder;

      // masters: the record layout is
      // [node_1 dof_1a dof_1b ... node_2 dof_2a ...] with the node
      // numbers positive and the dof keywords negative. A record may
      // have NO master terms (length 2): the slave is then bounded to
      // the constant zero (the -veln bounda generates such records for
      // axis-aligned walls, where the normal-velocity constraint
      // reduces to a single dof being zero).
      long int nterms = ( length_mpc - 2 ) / 2 + 1;
      long int *master_nodes = get_new_int( nterms );
      long int *master_dofs  = get_new_int( nterms );
      long int nmaster = 0, i = 2;
      while ( i<length_mpc ) {
        long int mnode = mpc[i++];
        if ( mnode<=0 ) db_error( MPC_NODE_NUMBER, impc );
        while ( i<length_mpc && mpc[i]<0 ) {
          long int iuknwn_m = -1;
          array_member( dof_label, mpc[i++], nuknwn, iuknwn_m );
          if ( iuknwn_m<0 ) db_error( MPC_NODE_NUMBER, impc );
          master_nodes[nmaster] = mnode;
          master_dofs[nmaster]  = iuknwn_m;
          nmaster++;
        }
      }

      // factors (manual 6.874): same order as the master dofs; missing
      // factors default to 1
      double *master_factors = get_new_dbl( nmaster );
      for ( i=0; i<nmaster; i++ ) master_factors[i] = 1.;
      long int length_fac = 0;
      db( MPC_NODE_FACTOR, impc, idum, master_factors, length_fac,
        VERSION_NORMAL, GET_IF_EXISTS );
      if ( length_fac>nmaster ) length_fac = nmaster;

      // records generated by a bounda_dof ... -veln record carry the
      // mpc_from_bounda marker. Their semantics (measured against the
      // Professional binary 25-10-2023) differ from explicit user
      // records on ONE point: a tie whose slave dof an explicit bounda
      // record of the current pass already bounded is INERT - the
      // direct prescription wins (at the inflow corner of the -veln
      // wall the Professional keeps the direct velx while the tie is
      // recorded but never applied). The free ties are treated like the
      // mpc_linear_quadratic ones: energy-consistent elimination inside
      // the solve (the slave row/column is redistributed to the masters,
      // see DIAG-SOLVE-MIXTO §16), so the wall constraint enters the
      // same solve that produces the masters instead of lagging one
      // iteration behind.
      long int from_bounda = -NO;
      db( MPC_FROM_BOUNDA, impc, &from_bounda, ddum, ldum, VERSION_NORMAL,
        GET_IF_EXISTS );

      // the eligibility inputs, evaluated per pass: the slave dof is free
      // unless a DIRECT bounda record of the current pass owns it (the
      // per-pass clearing above removed the flags the mpc machinery itself
      // set on earlier passes), and all master dofs are free unknowns.
      // Ties touching a directly prescribed dof keep the value-constrained
      // semantics (measured byte-identical: mpc1 Dirichlet master, mpc8
      // bounda'd slave).
      long int slave_free = 1;
      if ( db_active_index( NODE, node_0, VERSION_NORMAL ) ) {
        long int *nb = db_int( NODE_BOUNDED, node_0, VERSION_NORMAL );
        if ( nb[ipuknwn_0] ) slave_free = 0;
      }
      long int masters_free = 1;
      for ( i=0; i<nmaster; i++ ) {
        if ( !db_active_index( NODE, master_nodes[i], VERSION_NORMAL ) )
          masters_free = 0;
        else {
          long int *nbm = db_int( NODE_BOUNDED, master_nodes[i],
            VERSION_NORMAL );
          if ( nbm[master_dofs[i]/nder] ) masters_free = 0;
        }
      }

      // energy-consistent elimination: register the ELIGIBLE ties in the
      // slave map consumed by solve() (only when the slave node still
      // exists). A dangling quadratic node shared by several quadratic
      // elements on the non-conforming face is generated once PER element
      // by mpc_linear_quadratic_generate (identical records, measured on
      // the refined mpc4 hex27/hex8 interface); the map keeps the FIRST
      // registration so the slave is eliminated exactly once (a duplicated
      // slave row would be added twice to the master right-hand side).
      // Eligible classes (measured on the Professional binary 25-10-2023):
      // - the mpc_linear_quadratic generated range (DIAG §16),
      // - the mpc_element_group generated range with a free slave and free
      //   masters (DIAG §21 step 2; force16 1056/175, mpc6),
      // - the -veln generated free ties (from_bounda, validation_8),
      // - explicit user records with a free slave and free masters
      //   (mpc9; condensed diagonal K22+K33 = 11/6, v2 = 5/11).
      // Ties touching a DIRECTLY bounded dof (mpc7's bounded master, mpc1's
      // master, mpc8's slave) stay value-constrained in every class.
      {
        long int already = 0, ig0 = node_0*npuknwn + ipuknwn_0;
        for ( long int it=0; it<mpc_tie_count; it++ )
          if ( mpc_tie_slave_ig[it]==ig0 ) already = 1;
        long int in_lq_range = ( impc>=mpc_tie_elim_lo &&
          impc<=mpc_tie_elim_hi );
        long int in_eg_range = ( impc>=mpc_eg_lo && impc<=mpc_eg_hi );
        // DIAG §21 step 2: the mpc_element_group range is activated
        // (the exclusion of the pre-§21 state is removed): the generated
        // factors are machine-precision now and the momentum solves converge
        // to the condensed fixed point, so the eg ties condense exactly like
        // the explicit user ties (the Professional behaviour measured on
        // force16/mpc7).
        long int eligible = in_lq_range ||
          ( from_bounda==-YES && slave_free && nmaster>0 ) ||
          ( ( in_eg_range || from_bounda!=-YES ) && !in_lq_range &&
            slave_free && masters_free && nmaster>0 );
        if ( eligible && !already &&
             db_active_index( NODE, node_0, VERSION_NORMAL ) ) {
        // grow the tie arrays by one (capacity doubling: the map is
        // rebuilt every equilibrium iteration, keep it linear)
        long int n0 = mpc_tie_count;
        if ( n0>=mpc_tie_cap ) {
          long int new_cap = ( mpc_tie_cap>0 ) ? 2*mpc_tie_cap : 16;
          long int *slave_ig = get_new_int( new_cap );
          long int *off = get_new_int( new_cap );
          long int *nm = get_new_int( new_cap );
          for ( long int ic=0; ic<n0; ic++ ) {
            slave_ig[ic] = mpc_tie_slave_ig[ic];
            off[ic] = mpc_tie_off[ic];
            nm[ic] = mpc_tie_nm[ic];
          }
          delete[] mpc_tie_slave_ig;
          delete[] mpc_tie_off;
          delete[] mpc_tie_nm;
          mpc_tie_slave_ig = slave_ig;
          mpc_tie_off = off;
          mpc_tie_nm = nm;
          mpc_tie_cap = new_cap;
        }
        long int total = ( n0>0 ) ? mpc_tie_off[n0-1] + mpc_tie_nm[n0-1] : 0;
        if ( total+nmaster>mpc_tie_master_cap ) {
          long int new_cap = ( mpc_tie_master_cap>0 ) ?
            2*mpc_tie_master_cap : 16;
          while ( new_cap<total+nmaster ) new_cap *= 2;
          long int *master_ig = get_new_int( new_cap );
          double *master_fac = get_new_dbl( new_cap );
          for ( long int ic=0; ic<total; ic++ ) {
            master_ig[ic] = mpc_tie_master_ig[ic];
            master_fac[ic] = mpc_tie_master_fac[ic];
          }
          delete[] mpc_tie_master_ig;
          delete[] mpc_tie_master_fac;
          mpc_tie_master_ig = master_ig;
          mpc_tie_master_fac = master_fac;
          mpc_tie_master_cap = new_cap;
        }
        mpc_tie_slave_ig[n0] = node_0*npuknwn + ipuknwn_0;
        mpc_tie_off[n0] = total;
        mpc_tie_nm[n0] = nmaster;
        for ( i=0; i<nmaster; i++ ) {
          double factor = 1.;
          if ( i<length_fac ) factor = master_factors[i];
          mpc_tie_master_ig[total+i] =
            master_nodes[i]*npuknwn + master_dofs[i]/nder;
          mpc_tie_master_fac[total+i] = factor;
        }
        mpc_tie_count = n0 + 1;
        }
        // an eliminated tie OWNS the bounded flag it sets on the slave:
        // register the slave so the next pass clears the flag again before
        // re-deriving the eligibility (without this the flag set by pass k
        // would disqualify the tie at pass k+1 and the elimination would
        // be lost after the first iteration)
        if ( eligible && !already &&
             db_active_index( NODE, node_0, VERSION_NORMAL ) )
          mpc_applied_add( ig0 );
      }

      // apply: bound the slave dof and set its value from the current
      // master values (VERSION_NEW). Nodes that no longer exist (mesh
      // deletion) are skipped. A -veln generated tie whose slave dof an
      // explicit bounda record of the current pass already bounded is
      // INERT: the direct prescription wins and the tie value is never
      // forced on the dof (measured on the Professional binary at the
      // inflow corner of the -veln wall).
      if ( db_active_index( NODE, node_0, VERSION_NORMAL ) ) {
        long int *node_bounded = db_int( NODE_BOUNDED, node_0,
          VERSION_NORMAL );
        double *new_node_dof = db_dbl( NODE_DOF, node_0, VERSION_NEW );
        if ( !( from_bounda==-YES && node_bounded[ipuknwn_0] ) ) {
          node_bounded[ipuknwn_0] = 1;
          double value = 0.;
          for ( i=0; i<nmaster; i++ ) {
            if ( db_active_index( NODE, master_nodes[i], VERSION_NORMAL ) ) {
              double *mnd = db_dbl( NODE_DOF, master_nodes[i], VERSION_NEW );
              double factor = 1.;
              if ( i<length_fac ) factor = master_factors[i];
              value += factor * mnd[master_dofs[i]];
            }
          }
          new_node_dof[iuknwn_0] = value;
          // A tied VELOCITY slave also carries its integrated
          // displacement: parallel_new_dof_diagonal integrates
          // dis_new = dis_old + vel_new*dtime BEFORE the post-solve
          // re-sync above refreshes the slave velocity, so without this
          // the slave displacement lags one solve behind its velocity and
          // a -total_linear memory (strain from the total displacement)
          // sees unequal displacements at tied nodes (C41,
          // mohr_coul_direct7: disy -0.001808 / -0.001692 at the tied
          // nodes 3/4, Professional 02-08-2026 -0.0018458 at both).
          if ( ( materi_displacement || materi_displacement_relative ) &&
               iuknwn_0>=vel_indx && iuknwn_0<vel_indx+ndim*nder &&
               ( iuknwn_0-vel_indx )%nder==0 ) {
            long int idim_v = ( iuknwn_0-vel_indx ) / nder, ldum_v=0,
              idum_v[1];
            double dtime_v = 0.;
            double *old_node_dof = db_dbl( NODE_DOF, node_0, VERSION_NORMAL );
            db( DTIME, 0, idum_v, &dtime_v, ldum_v, VERSION_NEW,
              GET_IF_EXISTS );
            if ( materi_displacement )
              new_node_dof[dis_indx+idim_v*nder] =
                old_node_dof[dis_indx+idim_v*nder] + value * dtime_v;
            if ( materi_displacement_relative )
              new_node_dof[dis_rel_indx+idim_v*nder] =
                old_node_dof[dis_rel_indx+idim_v*nder] + value * dtime_v;
          }
        }
      }

      delete[] master_nodes;
      delete[] master_dofs;
      delete[] master_factors;
    }
  }

  // build the per-dof mark array (lookup for the solver-side elimination)
  db_max_index( NODE, max_node, VERSION_NORMAL, GET );
  mpc_tie_mark_size = ( 1 + max_node ) * npuknwn;
  delete[] mpc_tie_mark;
  mpc_tie_mark = get_new_int( mpc_tie_mark_size );
  array_set( mpc_tie_mark, 0, mpc_tie_mark_size );
  for ( long int it=0; it<mpc_tie_count; it++ ) {
    if ( mpc_tie_slave_ig[it]>=0 && mpc_tie_slave_ig[it]<mpc_tie_mark_size )
      mpc_tie_mark[mpc_tie_slave_ig[it]] = it + 1;
  }

}
