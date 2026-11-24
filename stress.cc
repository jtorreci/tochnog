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

extern "C" 
  int umat_(double *stress, double *statev, double *ddsdde, 
    double *sse, double *spd, double *scd, double *rpl, double *ddsddt,
    double *drplde, double *drpldt, double *stran, double *dstran,     double *time, double *dtime, double *temp, double *dtemp, double *predef, 
    double *dpred, char *cmname, long int *ndi, 
    long int *nshr, long int *ntens, long int *nstatv, 
    double *props, long int *nprops, double *coords, double *drot, 
    double *pnewdt, double *celent, double *dfgrd0, double *dfgrd1, 
    long int *noel, long int *npt, long int *layer, long int *kspt, 
    long int *kstep, long int *kinc, short cmname_len);

#define MAX_ITER 1000
// default cap of the plastic loop per integration point (manual
// Professional 6.804; owner 2026-11-06, BACKLOG C112)
#define MAX_PLASTI_ITER 20000
#define EPS_F 1.e-6
#define EPS_EPS_F 1.e-3
#define EPS_SIZE 1.e-3
#define EPS_DAMAGE 1.e-5
#define EPS_TMP 1.e-1
#define EPS_LAMBDA 1.e-3
#define EPS_VISCO 3.

// materi_compression_cutoff - group_materi_plasti_compression_direct
// (+_visco) (manual Professional 6.694/6.695): principal stresses lower
// than sigy are cut off (direct stress cut-off, no plastic strains).
// Spectral: sigma = V diag(d) V^T; each eigenvalue d_i < sigy is pulled
// up towards sigy with the visco factor 1-exp(-dtime/tm) (full cut when
// no visco record: factor 1).
void materi_compression_cutoff( long int element, long int gr,
  double dtime, double new_sig[], double &f_direct )

{
  long int ldum=0, idim=0, jdim=0, kdim=0, nrot=0;
  double sigy=0., tm=0., factor=1., sig_work[MDIM*MDIM],
    d[MDIM], v[MDIM*MDIM], plasti_data[DATA_ITEM_SIZE];

  if ( !get_group_data( GROUP_MATERI_PLASTI_COMPRESSION_DIRECT, gr,
      element, new_sig, plasti_data, ldum, GET_IF_EXISTS ) )
    return;
  sigy = plasti_data[0];
  if ( get_group_data( GROUP_MATERI_PLASTI_COMPRESSION_DIRECT_VISCO, gr,
      element, new_sig, plasti_data, ldum, GET_IF_EXISTS ) )
    tm = plasti_data[0];
  // factor 1 = full cut; with the _visco record the cut relaxes with
  // 1-exp(-dtime/tm)
  factor = 1.;
  if ( get_group_data( GROUP_MATERI_PLASTI_COMPRESSION_DIRECT_VISCO, gr,
      element, new_sig, plasti_data, ldum, GET_IF_EXISTS ) ) {
    tm = plasti_data[0];
    factor = 1. - exp( -dtime / ( (tm>0.) ? tm : 1. ) );
  }

  array_move( new_sig, sig_work, MDIM*MDIM );
  matrix_jacobi( sig_work, 3, d, v, &nrot );
  {
    long int any_cut = 0, i=0;
    double d_min = DBL_MAX;
    for ( i=0; i<3; i++ )
      if ( d[i]<sigy ) { d[i] += factor*(sigy-d[i]); any_cut = 1; }
    // yield value sigy - sig_min on the final stress (K13, Professional
    // 02-08-2026: -0.76e-3 for sigy -1e-3, sig_min -0.24e-3)
    for ( i=0; i<3; i++ ) if ( d[i]<d_min ) d_min = d[i];
    if ( sigy-d_min > f_direct ) f_direct = sigy - d_min;
    if ( !any_cut ) return;
  }
  // sigma = V diag(d) V^T  (V columns are eigenvectors)
  for ( idim=0; idim<MDIM; idim++ )
    for ( jdim=0; jdim<MDIM; jdim++ ) {
      double tmp = 0.;
      for ( kdim=0; kdim<MDIM; kdim++ )
        tmp += v[idim*MDIM+kdim] * d[kdim] * v[jdim*MDIM+kdim];
      new_sig[idim*MDIM+jdim] = tmp;
    }
}

// direct_cap_principal_values - step 1 of the direct cut-offs: principal
// stresses higher than sigy are cut off (manual Professional 6.738).
// Shared by the full principal cut-off (materi_direct_full_mc) and by the
// plane mode when the tension law has no plane normal of its own
// (materi_direct_principal_cap). Returns 1 when a value was capped.
// factor is the _visco relaxation 1-exp(-dtime/tm) (1 = full cut; C33,
// verified against the Professional 02-08-2026: tension_direct 1e-3 with
// _visco tm 2, one step dt 1, sigxx 2.4e-3 -> 1.84914e-3).
long int direct_cap_principal_values( double w[], double sigy, double factor )

{
  long int i=0, any_cut=0;

  for ( i=0; i<3; i++ ) {
    if ( w[i]>sigy ) { w[i] -= factor*( w[i]-sigy ); any_cut = 1; }
  }
  return any_cut;
}

// materi_direct_principal_cap - plane-mode tension step when the tension
// law has no own normal: cap the principal stresses at sigy (verified
// against the Professional 02-08-2026: mohr_coul_direct_normal without a
// tension normal cuts the principal tension first, then the MC plane).
void materi_direct_principal_cap( double new_sig[], double sigy,
  double factor )

{
  long int idim=0, jdim=0, kdim=0, nrot=0;
  double sig_work[MDIM*MDIM], d[MDIM], v[MDIM*MDIM];

  array_move( new_sig, sig_work, MDIM*MDIM );
  matrix_jacobi( sig_work, 3, d, v, &nrot );
  if ( !direct_cap_principal_values( d, sigy, factor ) ) return;
  // sigma = V diag(d) V^T  (V columns are eigenvectors)
  for ( idim=0; idim<MDIM; idim++ )
    for ( jdim=0; jdim<MDIM; jdim++ ) {
      double tmp = 0.;
      for ( kdim=0; kdim<MDIM; kdim++ )
        tmp += v[idim*MDIM+kdim] * d[kdim] * v[jdim*MDIM+kdim];
      new_sig[idim*MDIM+jdim] = tmp;
    }
}

// materi_direct_full_mc - group_materi_plasti_mohr_coul_direct +
// group_materi_plasti_tension_direct WITHOUT a plane normal (the manual
// Professional 6.726/6.738). This is the FULL Mohr-Coulomb "direct"
// stress cut-off (the alternative programming of the MC law, "very
// stable"): principal stresses higher than sigy are cut off (spectral
// tension cap, 6.738) and principal stress differences higher than
// allowed by the mohr-coulomb criterion are cut off (6.726):
//   f = 0.5(sig1-sig3) + 0.5(sig1+sig3) sin(phi) - c cos(phi) <= 0
// with sig1 the LARGEST and sig3 the SMALLEST principal stress
// (tension-positive). The cut is a one-shot return along the
// non-associative flow direction (phi_flow): for phi_flow = 0 it reduces
// to the mean-preserving difference cut, for phi_flow > 0 the pair mean
// shifts into compression (dilatant flow). When the middle principal sits
// at the sig1 level (the sigma1=sigma2 corner of the surface, e.g. a
// shear-only state after the tension cap) it is pulled down with sig1 so
// the state ends on the surface edge (verified against the Professional
// 25-10-2023 binary: the shear-only trial maps to
// (-a(1-sin), -a(1-sin), -a(1+sin)) in the principal frame with a the
// capped shear). When sigy is not given but mohr_coul_direct is
// available, sigy = 0 (6.738). _visco tm relaxes the cut with the factor
// 1-exp(-dtime/tm); _wall provides the parameters used when the element
// is attached to a wall. tension_cap = 0 skips the principal tension cap
// (step 1): used when the tension law has its own plane normal and was
// already applied on that plane (materi_direct_cutoff), while the
// Mohr-Coulomb law has no normal (verified against the Professional
// 02-08-2026: the MC law then stays the full principal cut-off).
// Wall reduction factor of the direct Mohr-Coulomb laws for an element on
// a wall: group_materi_plasti_mpc_factor for a wall found through mpc
// records (plasti_on_boundary 2, group.cc), else
// group_materi_plasti_bounda_factor; default 2/3 (manual 6.687, 6.734).
double direct_wall_factor( long int gr, long int plasti_on_boundary )
{
  long int ldum=0, idum[1];
  double factor=BOUNDARY_REDUCTION_FACTOR;
  if ( plasti_on_boundary==2 )
    db( GROUP_MATERI_PLASTI_MPC_FACTOR, gr, idum, &factor,
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
  else
    db( GROUP_MATERI_PLASTI_BOUNDARY_FACTOR, gr, idum, &factor,
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
  return factor;
}

void materi_direct_full_mc( long int element, long int gr,
  long int plasti_on_boundary, double dtime,
  double new_sig[], double /*ddsdde*/[], long int tension_cap, double &f_direct )

{
  long int idim=0, jdim=0, kdim=0, nrot=0, mc_present=0, ten_present=0,
    ldum=0;
  double phi=0., c=0., phi_flow=0., sigy=0., factor=1., factor_ten=1.,
    tm=0., young=0., poisson=0., lambda_lame=0., gmod=0.,
    sig_work[MDIM*MDIM], d[MDIM], v[MDIM*MDIM],
    plasti_data[DATA_ITEM_SIZE], ten_data[DATA_ITEM_SIZE];

  // elastic constants for the one-shot return (the C projection on the
  // principal flow direction); the direct cut runs on the elastic stress
  // so the elastic C is the active tangent
  get_group_data( GROUP_MATERI_ELASTI_YOUNG, gr, element,
    new_sig, &young, ldum, GET_IF_EXISTS );
  get_group_data( GROUP_MATERI_ELASTI_POISSON, gr, element,
    new_sig, &poisson, ldum, GET_IF_EXISTS );

  mc_present = get_group_data( GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT, gr,
    element, new_sig, plasti_data, ldum, GET_IF_EXISTS );
  if ( mc_present ) {
    phi = plasti_data[0]; c = plasti_data[1]; phi_flow = plasti_data[2];
    if ( plasti_on_boundary ) {
      if ( get_group_data( GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_WALL, gr,
          element, new_sig, plasti_data, ldum, GET_IF_EXISTS ) ) {
        phi = plasti_data[0]; c = plasti_data[1]; phi_flow = plasti_data[2];
      }
      else {
        // no _wall record: phi, c and phi_flow times the wall factor
        // (manual 6.687, default 2/3; Professional 02-08-2026, C114)
        double wall_factor = direct_wall_factor( gr, plasti_on_boundary );
        phi *= wall_factor; c *= wall_factor; phi_flow *= wall_factor;
      }
    }
  }
  // tension_direct: gate materi_plasti_tension_apply (6.804) /
  // control_materi_plasti_tension_apply (6.150)
  ten_present = !materi_apply_gate_off( MATERI_PLASTI_TENSION_APPLY,
    CONTROL_MATERI_PLASTI_TENSION_APPLY ) &&
    get_group_data( GROUP_MATERI_PLASTI_TENSION_DIRECT, gr,
    element, new_sig, ten_data, ldum, GET_IF_EXISTS );
  // 6.738: if tension_direct is not specified, sigy is set to 0 when the
  // mohr_coul_direct is available for the group.
  if ( ten_present ) sigy = ten_data[0];
  else if ( mc_present ) sigy = 0.;
  else return;   // nothing to do (no MC and no tension record)
  if ( ten_present && plasti_on_boundary )
    if ( get_group_data( GROUP_MATERI_PLASTI_TENSION_DIRECT_WALL, gr,
        element, new_sig, ten_data, ldum, GET_IF_EXISTS ) )
      sigy = ten_data[0];
  // visco relaxation: factor 1-exp(-dtime/tm); f->0 elastic, f->1 full
  // cut. Each law relaxes with its OWN _visco record: the tension cap with
  // tension_direct_visco, the Mohr-Coulomb cut with mohr_coul_direct_visco
  // (C33, verified against the Professional 02-08-2026 on a sheared quad4,
  // c 1e-3, phi 0: MC visco only -> full tension cap + relaxed MC; tension
  // visco only -> relaxed cap + full MC; both with tm 4 / 2 -> sigxy
  // 1.6881e-3).
  if ( get_group_data( GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_VISCO, gr,
      element, new_sig, plasti_data, ldum, GET_IF_EXISTS ) ) {
    tm = plasti_data[0];
    factor = 1. - exp( -dtime / ( (tm>0.) ? tm : 1. ) );
  }
  if ( ten_present && get_group_data( GROUP_MATERI_PLASTI_TENSION_DIRECT_VISCO,
      gr, element, new_sig, plasti_data, ldum, GET_IF_EXISTS ) ) {
    tm = plasti_data[0];
    factor_ten = 1. - exp( -dtime / ( (tm>0.) ? tm : 1. ) );
  }

  array_move( new_sig, sig_work, MDIM*MDIM );
  matrix_jacobi( sig_work, 3, d, v, &nrot );
  {
    double w[3], w_new[3], tmp=0., f=0.;
    long int i=0, j=0, any_cut=0, order[3];
    // sort ascending into w[0] <= w[1] <= w[2], tracking the eigenvector
    // columns of v (w[k] belongs to the eigenvector v[.][order[k]])
    for ( i=0; i<3; i++ ) { w[i] = d[i]; order[i] = i; }
    for ( i=0; i<3; i++ )
      for ( j=i+1; j<3; j++ )
        if ( w[i]>w[j] ) {
          tmp = w[i]; w[i] = w[j]; w[j] = tmp;
          tmp = order[i]; order[i] = order[j]; order[j] = (long int)tmp;
        }
    // 1) tension cap: principal stresses higher than sigy are cut off
    if ( tension_cap && direct_cap_principal_values( w, sigy, factor_ten ) )
      any_cut = 1;
    // re-sort after the cap (capping the maximum can change the order)
    for ( i=0; i<3; i++ )
      for ( j=i+1; j<3; j++ )
        if ( w[i]>w[j] ) {
          tmp = w[i]; w[i] = w[j]; w[j] = tmp;
          tmp = order[i]; order[i] = order[j]; order[j] = (long int)tmp;
        }
    // 2) Mohr-Coulomb principal-stress-difference cut
    if ( mc_present ) {
      // w[2] largest, w[0] smallest (tension-positive)
      f = 0.5*(w[2]-w[0]) + 0.5*(w[2]+w[0])*sin(phi) - c*cos(phi);
      if ( f>0. ) {
        // One-shot return along the non-associative flow direction of the
        // plastic potential g = 0.5(s1-s3) + 0.5(s1+s3) sin(phi_flow):
        //   deps = lambda * (0.5(1+sin psi), 0, -0.5(1-sin psi))  (princ.)
        // with lambda solved so the corrected stress sits on the MC
        // surface. The stress correction -C:deps couples the principal
        // directions through the isotropic C (the intermediate principal
        // moves when the flow is dilatant, sin psi != 0). For psi = 0
        // the correction is the mean-preserving difference cut.
        double sf = sin(phi_flow);
        if ( young>0. && poisson<0.5 ) {
          lambda_lame = young*poisson/((1.+poisson)*(1.-2.*poisson));
          gmod = young/(2.*(1.+poisson));
        }
        else {   // no elastic data: decoupled correction (E-normalized)
          lambda_lame = 0.;
          gmod = 0.5;
        }
        // (C:d) components (principal frame, isotropic C)
        double cd1 = lambda_lame*sf + 2.*gmod*0.5*(1.+sf);
        double cd2 = lambda_lame*sf;
        double cd3 = lambda_lame*sf - 2.*gmod*0.5*(1.-sf);
        double denom = 0.5*(cd1-cd3) + 0.5*(cd1+cd3)*sin(phi);
        if ( denom<1.e-12 ) denom = 1.e-12;
        double dlam = f/denom;
        w_new[2] = w[2] - dlam*cd1;
        w_new[0] = w[0] - dlam*cd3;
        // middle principal: the isotropic C couples the flow to the
        // intermediate direction (dilatant flow moves it, lambda_lame*sf);
        // when it then sits at (or above) the new maximum (the sigma1 =
        // sigma2 corner, e.g. the shear-only state after the tension cap)
        // it is pulled down to the new maximum so the state lands on the
        // surface edge (Professional behaviour).
        w_new[1] = w[1] - dlam*cd2;
        if ( w_new[1] > w_new[2] ) w_new[1] = w_new[2];
        // apply with the visco factor (partial relaxation)
        for ( i=0; i<3; i++ ) {
          double corr = (w_new[i]-w[i]);
          if ( factor<1. ) corr *= factor;
          if ( corr!=0. ) any_cut = 1;
          w_new[i] = w[i] + corr;
        }
        for ( i=0; i<3; i++ ) w[i] = w_new[i];
      }
    }
    // yield value on the final principal stresses (K13): the tension
    // cap sig_max - sigy and the MC function, the largest one
    {
      double w_max = w[0], w_min = w[0], fv = 0.;
      for ( i=1; i<3; i++ ) {
        if ( w[i]>w_max ) w_max = w[i];
        if ( w[i]<w_min ) w_min = w[i];
      }
      if ( tension_cap && w_max-sigy > f_direct ) f_direct = w_max - sigy;
      if ( mc_present ) {
        fv = 0.5*(w_max-w_min) + 0.5*(w_max+w_min)*sin(phi) - c*cos(phi);
        if ( fv > f_direct ) f_direct = fv;
      }
    }
    if ( any_cut ) {
      // rebuild sigma = V diag(w) V^T (V columns are eigenvectors;
      // w[k] belongs to the column order[k] of v)
      double vv[MDIM*MDIM];
      for ( i=0; i<3; i++ )
        for ( j=0; j<3; j++ )
          vv[i*MDIM+j] = v[i*MDIM+order[j]];
      for ( idim=0; idim<MDIM; idim++ )
        for ( jdim=0; jdim<MDIM; jdim++ ) {
          double tmp2 = 0.;
          for ( kdim=0; kdim<MDIM; kdim++ )
            tmp2 += vv[idim*MDIM+kdim] * w[kdim] * vv[jdim*MDIM+kdim];
          new_sig[idim*MDIM+jdim] = tmp2;
        }
    }
  }
}

// direct_flow_factor - M/G for the plane-mode phi_flow return of
// materi_direct_cutoff: M = n_i n_j C_ijkl n_k n_l (constrained modulus
// normal to the plane) and G = 0.5 t_i n_j C_ijkl (t_k n_l + n_k t_l)
// (shear modulus for a slip along t = tau/|tau| on the plane), from the
// 4th-order elastic stiffness C (MDIM^4, row-major, as Cmem in
// set_stress). Returns 2 (the Professional constant, M/G at poisson = 0)
// when C gives no positive M and G.
double direct_flow_factor( double C[], double n[], double tau[],
  double tau_norm )

{
  long int i=0, j=0, k=0, l=0;
  double M=0., G=0., t[MDIM], cijkl=0.;

  for ( i=0; i<MDIM; i++ ) t[i] = ( i<ndim ) ? tau[i]/tau_norm : 0.;
  for ( i=0; i<MDIM; i++ )
    for ( j=0; j<MDIM; j++ )
      for ( k=0; k<MDIM; k++ )
        for ( l=0; l<MDIM; l++ ) {
          cijkl = C[((i*MDIM+j)*MDIM+k)*MDIM+l];
          if ( cijkl==0. ) continue;
          M += n[i]*n[j]*cijkl*n[k]*n[l];
          G += 0.5*t[i]*n[j]*cijkl*( t[k]*n[l] + n[k]*t[l] );
        }
  if ( M>0. && G>0. ) return M/G;
  return 2.;
}

// materi_direct_cutoff - group_materi_plasti_mohr_coul_direct(_normal[_
// automatic]) + group_materi_plasti_tension_direct(_normal[_automatic]),
// plane mode (at least one of the two laws has a plane normal; manual
// Professional 6.725-6.740, verified against the Professional 02-08-2026).
//
// The "_direct" plastic laws are DIRECT STRESS CUT-OFFS (no plastic
// strains; the manual: "cut off by Tochnog", tension_direct "does not use
// plastic strains"). set_stress stores the cut afterwards as the plastic
// strain C^-1 (sig_before - sig_after) and the yield value f (K13). Each law has its OWN plane (mc_normal, ten_normal;
// zero when the law has no normal). On a plane with normal n:
//   traction    t  = sig . n
//   normal      sig_n = n . t
//   tangential  tau  = t - sig_n n
// (tochnog stress convention: tension POSITIVE).
//
// 1) Tension first. sigy = tension_direct, or 0 when only
//    mohr_coul_direct is given (manual 6.738). With its own normal the
//    tension law caps sig_n <= sigy on that plane; without it the
//    principal stresses are capped at sigy.
// 2) Then Mohr-Coulomb on the MC plane with the updated stress (sn, tau
//    = |tau|). When tau > max(c - sn tan(phi), 0) the state returns along
//    the non-associative flow direction (phi_flow = psi), changing ONLY
//    sn and tau:
//      k       = (M/G) tan(psi)
//      tau_new = max( (c - sn tan(phi) + k tan(phi) tau) / (1 + k tan(phi)), 0 )
//      sn_new  = sn - k (tau - tau_new)
//    M = n.C.n.n.n is the constrained (oedometric) modulus and G the
//    shear modulus on the plane, both from the elastic stiffness C of
//    the element (elastic_stiffness, built in set_stress); the ratio
//    M/G is the consistent elastic return of the slip plastic strain
//    (shear slip gives dtau = -G dlambda, dilatancy gives dsn = -M
//    tan(psi) dlambda). Without usable elastic data (M or G <= 0, e.g.
//    a umat law) k = 2 tan(psi), the Professional constant (equal to
//    M/G only for poisson = 0). Documented divergence: the Professional
//    always uses k = 2 tan(psi). The floor tau_new = 0 keeps sn_new =
//    sn - k tau (verified against the Professional beyond the apex).
//    Without an MC normal the MC law is the full principal cut-off
//    (materi_direct_full_mc without its tension step).
void materi_direct_cutoff( long int element, long int gr,
  long int plasti_on_boundary, double dtime,
  double new_sig[], double ddsdde[], double mc_normal[],
  double ten_normal[], double elastic_stiffness[], double &f_direct )

{
  long int idim=0, jdim=0, mc_active=0, ten_record=0, ten_active=0,
    mc_plane=0, ten_plane=0, ldum=0;
  double phi=0., c=0., phi_flow=0., sigy=0., sig_n=0., tau_norm=0.,
    max_fric=0., scale=0., factor_mc=1., factor_ten=1., n_mc[MDIM],
    n_ten[MDIM], tau[MDIM], t[MDIM], plasti_data[DATA_ITEM_SIZE];

  array_move( mc_normal, n_mc, MDIM );
  array_move( ten_normal, n_ten, MDIM );
  mc_plane = array_normalize( n_mc, MDIM ) ? 1 : 0;
  ten_plane = array_normalize( n_ten, MDIM ) ? 1 : 0;

  mc_active = get_group_data( GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT, gr,
    element, new_sig, plasti_data, ldum, GET_IF_EXISTS );
  if ( mc_active ) {
    phi = plasti_data[0]; c = plasti_data[1]; phi_flow = plasti_data[2];
  }
  // wall values (element attached to a wall): the _wall record, else
  // phi, c and phi_flow times the wall factor (manual 6.687, default 2/3;
  // Professional 02-08-2026: no factor on top of the _wall values, C114)
  if ( mc_active && plasti_on_boundary ) {
    if ( get_group_data( GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_WALL, gr,
        element, new_sig, plasti_data, ldum, GET_IF_EXISTS ) ) {
      phi = plasti_data[0]; c = plasti_data[1]; phi_flow = plasti_data[2];
    }
    else {
      double wall_factor = direct_wall_factor( gr, plasti_on_boundary );
      phi *= wall_factor; c *= wall_factor; phi_flow *= wall_factor;
    }
  }
  // visco factor: f = 1 - exp(-dt/tm); f->0 elastic (dt<<tm), f->1 full cut
  if ( mc_active ) {
    if ( get_group_data( GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_VISCO, gr,
        element, new_sig, plasti_data, ldum, GET_IF_EXISTS ) ) {
      double tm = plasti_data[0];
      factor_mc = 1. - exp( -dtime / ( (tm>0.) ? tm : 1. ) );
    }
  }

  // materi_plasti_tension_apply / control_materi_plasti_tension_apply
  // -no (manual Professional 6.804 / 6.150): ignore any tension-plasticity
  // data (the Mohr-Coulomb direct cutoff stays active, with its default
  // sigy = 0 as in materi_direct_full_mc)
  ten_record = !materi_apply_gate_off( MATERI_PLASTI_TENSION_APPLY,
    CONTROL_MATERI_PLASTI_TENSION_APPLY ) &&
    get_group_data( GROUP_MATERI_PLASTI_TENSION_DIRECT, gr,
    element, new_sig, plasti_data, ldum, GET_IF_EXISTS );
  if ( ten_record ) {
    sigy = plasti_data[0];
    if ( plasti_on_boundary ) {
      if ( get_group_data( GROUP_MATERI_PLASTI_TENSION_DIRECT_WALL, gr,
          element, new_sig, plasti_data, ldum, GET_IF_EXISTS ) ) {
        sigy = plasti_data[0];
      }
    }
    if ( get_group_data( GROUP_MATERI_PLASTI_TENSION_DIRECT_VISCO, gr,
        element, new_sig, plasti_data, ldum, GET_IF_EXISTS ) ) {
      double tm = plasti_data[0];
      factor_ten = 1. - exp( -dtime / ( (tm>0.) ? tm : 1. ) );
    }
  }
  // 6.738: if tension_direct is not specified, sigy is set to 0 when the
  // mohr_coul_direct is available for the group
  ten_active = ten_record || mc_active;

  // 1) tension step
  if ( ten_active ) {
    if ( ten_plane ) {
      // cap the normal traction on the tension plane
      sig_n = 0.;
      for ( idim=0; idim<ndim; idim++ )
        for ( jdim=0; jdim<ndim; jdim++ )
          sig_n += n_ten[idim] * new_sig[idim*MDIM+jdim] * n_ten[jdim];
      if ( sig_n > sigy ) {
        double corr = ( sig_n - sigy ) * factor_ten;
        for ( idim=0; idim<ndim; idim++ )
          for ( jdim=0; jdim<ndim; jdim++ )
            new_sig[idim*MDIM+jdim] -= corr * n_ten[idim] * n_ten[jdim];
      }
    }
    else
      materi_direct_principal_cap( new_sig, sigy, factor_ten );
  }

  // 2) Mohr-Coulomb step
  if ( mc_active && !mc_plane )
    materi_direct_full_mc( element, gr, plasti_on_boundary, dtime,
      new_sig, ddsdde, 0, f_direct );
  if ( mc_active && mc_plane ) {
    // traction on the MC plane, from the stress after the tension step
    for ( idim=0; idim<ndim; idim++ ) {
      t[idim] = 0.;
      for ( jdim=0; jdim<ndim; jdim++ )
        t[idim] += new_sig[idim*MDIM+jdim] * n_mc[jdim];
    }
    sig_n = 0.;
    for ( idim=0; idim<ndim; idim++ ) sig_n += n_mc[idim] * t[idim];
    for ( idim=0; idim<ndim; idim++ ) tau[idim] = t[idim] - sig_n*n_mc[idim];
    tau_norm = 0.;
    for ( idim=0; idim<ndim; idim++ ) tau_norm += tau[idim]*tau[idim];
    tau_norm = sqrt( tau_norm );

    max_fric = c - sig_n * tan( phi );
    if ( max_fric < 0. ) max_fric = 0.;
    if ( tau_norm > max_fric && tau_norm > 0. ) {
      double tan_phi = tan( phi ), k_flow = 0., tau_new = 0., dsig_n = 0.;
      if ( phi_flow!=0. )
        k_flow = direct_flow_factor( elastic_stiffness, n_mc, tau,
          tau_norm ) * tan( phi_flow );
      if ( 1.+k_flow*tan_phi < 1.e-12 ) k_flow = 0.;   // degenerate data
      tau_new = ( c - sig_n*tan_phi + k_flow*tan_phi*tau_norm ) /
        ( 1. + k_flow*tan_phi );
      if ( tau_new < 0. ) tau_new = 0.;
      dsig_n = - k_flow * ( tau_norm - tau_new ) * factor_mc;
      scale = tau_new / tau_norm;
      scale = 1. - factor_mc*(1.-scale);   // visco: partial
      // new_sig -= (1-scale)*(tau x n + n x tau)  (symmetric correction)
      // new_sig += dsig_n n x n                   (dilatant normal shift)
      for ( idim=0; idim<ndim; idim++ ) {
        for ( jdim=0; jdim<ndim; jdim++ ) {
          double corr = (1.-scale) * ( tau[idim]*n_mc[jdim] +
            n_mc[idim]*tau[jdim] );
          new_sig[idim*MDIM+jdim] -= corr;
          new_sig[idim*MDIM+jdim] += dsig_n * n_mc[idim] * n_mc[jdim];
        }
      }
    }
  }

  // yield value of the direct laws on the final stress (K13, measured on
  // the Professional 02-08-2026): tension sig_n - sigy on its plane (or
  // sig_max - sigy without a tension normal), Mohr-Coulomb |tau| - (c -
  // sig_n tan(phi)) on its plane; the full MC value comes from
  // materi_direct_full_mc. The largest one is kept.
  if ( ten_active ) {
    double fv = 0.;
    if ( ten_plane ) {
      for ( idim=0; idim<ndim; idim++ )
        for ( jdim=0; jdim<ndim; jdim++ )
          fv += n_ten[idim] * new_sig[idim*MDIM+jdim] * n_ten[jdim];
    }
    else {
      long int nrot=0;
      double sig_work[MDIM*MDIM], d[MDIM], v[MDIM*MDIM];
      array_move( new_sig, sig_work, MDIM*MDIM );
      matrix_jacobi( sig_work, 3, d, v, &nrot );
      fv = d[0];
      for ( idim=1; idim<3; idim++ ) if ( d[idim]>fv ) fv = d[idim];
    }
    if ( fv-sigy > f_direct ) f_direct = fv - sigy;
  }
  if ( mc_active && mc_plane ) {
    double fv = 0.;
    for ( idim=0; idim<ndim; idim++ ) {
      t[idim] = 0.;
      for ( jdim=0; jdim<ndim; jdim++ )
        t[idim] += new_sig[idim*MDIM+jdim] * n_mc[jdim];
    }
    sig_n = 0.;
    for ( idim=0; idim<ndim; idim++ ) sig_n += n_mc[idim] * t[idim];
    tau_norm = 0.;
    for ( idim=0; idim<ndim; idim++ )
      tau_norm += ( t[idim] - sig_n*n_mc[idim] ) *
        ( t[idim] - sig_n*n_mc[idim] );
    fv = sqrt( tau_norm ) - ( c - sig_n * tan( phi ) );
    if ( fv > f_direct ) f_direct = fv;
  }
}

// direct_plastic_strain - K13: the tensor strain deps with C : deps =
// dsig for the elastic stiffness C = Cmem + Cuser + Chyper + Chypo (+
// 1/compressibility on the volumetric part), solved in the 6-component
// tensor form (shear columns 2 C_ijkl). Returns 0 when C is singular.
long int direct_plastic_strain( double dsig[],
  double Cmem[MDIM][MDIM][MDIM][MDIM], double Cuser[MDIM][MDIM][MDIM][MDIM],
  double Chyper[MDIM][MDIM][MDIM][MDIM], double Chypo[MDIM][MDIM][MDIM][MDIM],
  double compressibility, double deps[] )

{
  const long int ii[6] = { 0, 1, 2, 0, 0, 1 }, jj[6] = { 0, 1, 2, 1, 2, 2 };
  long int r=0, q=0, k=0, piv=0;
  double a[6][7], cijkl=0., tmp=0., amax=0.;

  for ( r=0; r<6; r++ ) {
    for ( q=0; q<6; q++ ) {
      long int i=ii[r], j=jj[r], kk=ii[q], l=jj[q];
      cijkl = Cmem[i][j][kk][l] + Cuser[i][j][kk][l] + Chyper[i][j][kk][l] +
        Chypo[i][j][kk][l];
      if ( compressibility!=0. && i==j && kk==l ) cijkl += 1./compressibility;
      if ( kk!=l ) {
        cijkl += Cmem[i][j][l][kk] + Cuser[i][j][l][kk] + Chyper[i][j][l][kk] +
          Chypo[i][j][l][kk];
      }
      a[r][q] = cijkl;
      if ( scalar_dabs(cijkl)>amax ) amax = scalar_dabs(cijkl);
    }
    a[r][6] = dsig[ii[r]*MDIM+jj[r]];
  }
  if ( amax<=0. ) return 0;
  // Gauss-Jordan elimination with partial pivoting; a singular C (no
  // usable elastic stiffness) stores no plastic strain
  for ( k=0; k<6; k++ ) {
    piv = k;
    for ( r=k+1; r<6; r++ )
      if ( scalar_dabs(a[r][k])>scalar_dabs(a[piv][k]) ) piv = r;
    if ( scalar_dabs(a[piv][k])<=1.e-12*amax ) return 0;
    if ( piv!=k )
      for ( q=0; q<7; q++ ) { tmp = a[k][q]; a[k][q] = a[piv][q]; a[piv][q] = tmp; }
    for ( r=0; r<6; r++ ) {
      if ( r==k || a[r][k]==0. ) continue;
      tmp = a[r][k] / a[k][k];
      for ( q=k; q<7; q++ ) a[r][q] -= tmp*a[k][q];
    }
  }
  array_set( deps, 0., MDIM*MDIM );
  for ( r=0; r<6; r++ ) {
    tmp = a[r][6] / a[r][r];
    deps[ii[r]*MDIM+jj[r]] = tmp;
    deps[jj[r]*MDIM+ii[r]] = tmp;
  }
  return 1;
}

// group_materi_plasti_visco_power_name / _value and
// group_materi_plasti_visco_exponential_name / _values: the two visco
// constants {eta p, or gamma alpha} of the active plasticity model
// plasti_type {sub-surfaces map to their record}.
// Returns 1 and fills visco_pair[0..1] when the model is named.
static long int visco_named_select( long int plasti_type, long int n_names,
  long int names[], double values[], double visco_pair[] )
{
  long int k=0, record=plasti_type;
  if ( record==GROUP_MATERI_PLASTI_DRUCKPRAG_APEX )
    record = GROUP_MATERI_PLASTI_DRUCKPRAG;
  else if ( record==GROUP_MATERI_PLASTI_MATSUOKANAKAI_APEX )
    record = GROUP_MATERI_PLASTI_MATSUOKANAKAI;
  else if ( record==GROUP_MATERI_PLASTI_MOHRCOUL_APEX ||
            record==GROUP_MATERI_PLASTI_MOHRCOUL_01 ||
            record==GROUP_MATERI_PLASTI_MOHRCOUL_12 ||
            record==GROUP_MATERI_PLASTI_MOHRCOUL_20 )
    record = GROUP_MATERI_PLASTI_MOHRCOUL;
  for ( k=0; k<n_names; k++ ) {
    if ( labs(names[k])==record ) {
      visco_pair[0] = values[2*k];
      visco_pair[1] = values[2*k+1];
      return 1;
    }
  }
  return 0;
}

void set_stress( long int element, long int gr, 
  long int plasti_on_boundary, double coord_ip[],
  double old_unknowns[], double new_unknowns[], 
  double old_grad_old_unknowns[], double new_grad_new_unknowns[], 
  double rotated_old_sig[], double new_sig[], double new_sig_eff[],
  double rotated_old_msig[], double new_msig[], 
  double inc_ept[], double new_ept[], 
  double old_epe[], double inc_epe[], 
  double old_epp[], double inc_epp[], 
  double old_rho[], double new_rho[],
  double old_epi[], double new_epi[],
  double old_hisv[], double new_hisv[], 
  double old_damage, double &new_damage, 
  double old_kappa, double &new_kappa, 
  double old_kapsh, double &new_kapsh, 
  double old_cap1pc, double &new_cap1pc, 
  double &new_f, double &new_substeps, double old_deften[], double new_deften[],
  double inc_rot[], double ddsdde[],
   double &viscosity, double &viscosity_heat_generation, double &softvar_nonl,
   double &softvar_l, double mc_normal[], double ten_normal[] )

  // Solid materials.

{
  long int i=0, j=0, ind_ddsdde=0, plasti_found=0, plasti_off=0, 
    plasti_iter=0, membrane_found=0, membrane_iter=0,
    swit=0, length=0, membrane=-NO, viscoplasti=0, 
    viscoplasti_always=-NO, plasti_type=-NONE, volumetric_young_order=0,
    visco_power_length=0,
    memory=-UPDATED, max_plasti_iter=0, total_plasti_iter=0, 
    nuser_data=0, idim=0, jdim=0, kdim=0, ldim=0, k0_active=0, 
    k0_control_swit=0,
    formulation=INCREMENTAL, ldum=0, idum[1], task[2],
    hs_gp_swit=-NO, icontrol_hs=0, hs_plasti_present=0, hs_len=1,
    damage_done=0, direct_on=0, plastic_ip=0, visco_named=0,
    n_visco_names=0, visco_names[DATA_ITEM_SIZE], max_iter_visco=0,
    max_iter_ep=0, visco_idum[DATA_ITEM_SIZE];
  double lambda=0., deps_size=0., lambda_new=0., lambda_previous=0., 
    tmp=0., tmp_old=0., tmp_inc=0., tmp_new=0., 
    f=0., f_previous=0, f_ref=0., materi_expansion_linear=0., 
    visco_limit=EPS_VISCO, visco_pair[2], 
    eta=0., pressure=0., strain_size=0., straindev_size=0., linear_ept=0.,zero_ept=0.,
	 meanstrain=0.,
    plasti_kinematic_hardening=0., young = 0., young_linear_ept=0., poisson=0., 
    compressibility=0., fac=0., kappa=0., g=0., k=0., e=0.,
    lade_1=0., lade_2=0, lade_3=0., p=0., p0=0., p1=0., young0=0., young1=0.,
    young2=0., young3=0.,
    nu0=0., nu1=0., nu2=0., nu50=0., nuur=0.,
    alpha=0., gamma=0., dtime=0., rdum=0., camclay[1], tskh[DATA_ITEM_SIZE],
    smallstrain[6],
    group_materi_elasti_lade[3], sig_dev[MDIM*MDIM], ddum[MDIM*MDIM], ddumarray[MDIM][MDIM], 
    sig_princ[MDIM],
    hs_elasti[7], hs_plasti[4],
    direct_f=NO_YIELD_F, direct_sig_before[MDIM*MDIM],
    hs_sig3=0., hs_base=0., hs_E50=0., hs_Eur=0., hs_Eref50=0., hs_sigref50=0.,
    hs_Erefur=0., hs_sigrefur=0., hs_m=0., hs_ccot=0., hs_qf=0.,
    hs_qa=0., hs_gp_extra=0., hs_phi=0., hs_c=0.,
    memmat[MDIM][MDIM], elasti_transverse_isotropy[DATA_ITEM_SIZE],
    inc_temperature_strain[MDIM*MDIM], new_temperature_strain[MDIM*MDIM], 
    plasti_dir[MDIM*MDIM], young_power[6], young_polynomial[DATA_ITEM_SIZE],
    poisson_power[5], shear_factor=0., k0_elasti=0., sph_factor=0.,
    young_strainstress[DATA_ITEM_SIZE], ept_dev[MDIM*MDIM],
    plasti_visco_exponential[2], plasti_visco_power[3], 
    inc_rho[MDIM*MDIM], test_sig[MDIM*MDIM], total_inc_epp[MDIM*MDIM], 
    work_inc_epp[MDIM*MDIM], new_epe[MDIM*MDIM], new_epp[MDIM*MDIM],
    C[MDIM][MDIM][MDIM][MDIM], Cmem[MDIM][MDIM][MDIM][MDIM], 
    Cuser[MDIM][MDIM][MDIM][MDIM], Chyper[MDIM][MDIM][MDIM][MDIM],
    Chypo[MDIM][MDIM][MDIM][MDIM],
    volumetric_young_values[DATA_ITEM_SIZE],
    user_data[DATA_ITEM_SIZE], work[DATA_ITEM_SIZE], 
    old_work[MDIM*MDIM], new_work[MDIM*MDIM], visco_values[DATA_ITEM_SIZE];

  swit = set_swit(element,-1, "set_stress");
  if ( swit ) pri( "In routine SET_STRESS" );

  plasti_visco_exponential[0] = 0.;
  plasti_visco_power[0] = 0.;
  if ( get_group_data( GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL, gr, element, 
      new_unknowns, plasti_visco_exponential, ldum, GET_IF_EXISTS ) ) {
    viscoplasti = 1;
    gamma = plasti_visco_exponential[0]; alpha = plasti_visco_exponential[1];
    pressure = ( new_sig[0] + new_sig[4] + new_sig[8] ) / 3.;
    if ( alpha<=0. || gamma<=0. ) 
      db_error(  GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL, gr );
  }
  else if ( get_group_data( GROUP_MATERI_PLASTI_VISCO_POWER, gr, element, 
      new_unknowns, plasti_visco_power, visco_power_length, GET_IF_EXISTS ) ) {
    viscoplasti = 1;
    eta = plasti_visco_power[0]; p = plasti_visco_power[1];
    // Professional layout (eta p, manual 6.748): the power law reads
    // eps_dot_pl = eta * f^p (no reference stress). The legacy GNU
    // layout (eta p f_ref) is still accepted and keeps its f/f_ref role.
    if ( visco_power_length>=3 ) f_ref = plasti_visco_power[2];
    else f_ref = 1.;
    assert( f_ref!=0. );
  }
  // group_materi_plasti_visco_power_name / _value (manual Professional
  // 6.748, 6.749): eta and p per plasticity model; a model not in the
  // names stays elasto-plastic. The choice is made per integration point
  // when the yield test selects the active model {visco_power_select}.
  // group_materi_plasti_visco_power wins when both are given, as the
  // Professional 02-08-2026.
  else if ( db( GROUP_MATERI_PLASTI_VISCO_POWER_NAME, gr, visco_names,
      work, n_visco_names, VERSION_NORMAL, GET_IF_EXISTS ) ) {
    long int n_values=0;
    visco_named = 1;
    if ( !db( GROUP_MATERI_PLASTI_VISCO_POWER_VALUE, gr, visco_idum,
          visco_values, n_values, VERSION_NORMAL, GET_IF_EXISTS ) ||
         n_values!=2*n_visco_names ) {
      pri( "Error: group_materi_plasti_visco_power_value needs two values "
        "{eta p} for each name of group_materi_plasti_visco_power_name." );
      exit(1);
    }
  }
  // group_materi_plasti_visco_exponential_name / _values (manual
  // Professional 6.745, 6.746): the same per-model selection for the
  // exponential law, with gamma alpha per name. The Professional
  // 02-08-2026 reads _name as one integer and refuses model names, so
  // these records are implemented after the manual.
  else if ( db( GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL_NAME, gr, visco_names,
      work, n_visco_names, VERSION_NORMAL, GET_IF_EXISTS ) ) {
    long int n_values=0, k=0;
    visco_named = 2;
    if ( !db( GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL_VALUES, gr, visco_idum,
          visco_values, n_values, VERSION_NORMAL, GET_IF_EXISTS ) ||
         n_values!=2*n_visco_names ) {
      pri( "Error: group_materi_plasti_visco_exponential_values needs two "
        "values {gamma alpha} for each name of "
        "group_materi_plasti_visco_exponential_name." );
      exit(1);
    }
    for ( k=0; k<n_values; k++ ) {
      if ( visco_values[k]<=0. )
        db_error( GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL_VALUES, gr );
    }
    pressure = ( new_sig[0] + new_sig[4] + new_sig[8] ) / 3.;
  }
  // group_materi_plasti_visco_exponential_limit (manual Professional
  // 6.744): the cap of the exponent argument alpha*f, default 3.
  get_group_data( GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL_LIMIT, gr, element,
    new_unknowns, &visco_limit, ldum, GET_IF_EXISTS );

  // time step of the CURRENT step (VERSION_NEW, as every other module):
  // the _visco direct cut-offs relax with 1-exp(-dt/tm) over the step from
  // t to t+dt (manual Professional 6.728) and the viscoplastic laws
  // integrate lambda*dt over the same step. VERSION_NORMAL held the
  // previous step's dt (0 in the first step: no visco relaxation at all;
  // C33, verified against the Professional 02-08-2026).
  db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET_IF_EXISTS );
  // materi_plasti_visco_apply / control_materi_plasti_visco_apply -no
  // (manual Professional 6.805 / 6.151): ignore any visco-plasticity data
  if ( materi_apply_gate_off( MATERI_PLASTI_VISCO_APPLY,
       CONTROL_MATERI_PLASTI_VISCO_APPLY ) ) {
    viscoplasti = 0;
    visco_named = 0;
  }
  db( GROUP_MATERI_PLASTI_VISCO_ALWAYS, gr, &viscoplasti_always,
    ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  memory = materi_memory_default();
  db( GROUP_MATERI_MEMORY, gr, &memory, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  // control_materi_updated_apply (manual Professional 6.152): the -yes
  // direction (unspecified memory becomes -updated) lives in
  // materi_memory_default(); the -no direction below sets any -updated
  // material memory to -updated_linear for these timesteps.
  if ( control_materi_gate_off( CONTROL_MATERI_UPDATED_APPLY ) ) {
    if ( memory==-UPDATED || memory==-UPDATED_WITHOUT_ROTATION )
      memory = -UPDATED_LINEAR;
  }
  array_set( &C[0][0][0][0], 0., MDIM*MDIM*MDIM*MDIM );
  array_set( &Cmem[0][0][0][0], 0., MDIM*MDIM*MDIM*MDIM );
  get_group_data( GROUP_MATERI_ELASTI_COMPRESSIBILITY, gr, element, new_unknowns, 
    &compressibility, ldum, GET_IF_EXISTS );
  db( GROUP_MATERI_MEMBRANE, gr, &membrane, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  task[1] = membrane;
  if ( get_group_data( GROUP_MATERI_ELASTI_TRANSVERSE_ISOTROPY, gr, 
      element, new_unknowns, elasti_transverse_isotropy, ldum, GET_IF_EXISTS ) ) {
    task[0] = GROUP_MATERI_ELASTI_TRANSVERSE_ISOTROPY;
    task[1] = -NO;
    C_matrix( young, poisson, elasti_transverse_isotropy, C, task );
    task[1] = membrane;
    C_matrix( young, poisson, elasti_transverse_isotropy, Cmem, task );
  }

  task[0] = GROUP_MATERI_ISOTROPY; //always if not changed by TRANSVERSE_ISOTROPY_GRAHOUL
  if ( get_group_data( GROUP_MATERI_ELASTI_TRANSVERSE_ISOTROPY_GRAHOUL, gr, 
      element, new_unknowns, elasti_transverse_isotropy, ldum, GET_IF_EXISTS ) ) {
    task[0] = GROUP_MATERI_ELASTI_TRANSVERSE_ISOTROPY_GRAHOUL;
  }

  // group_materi_elasti_k0 (manual Professional 6.650): when this record
  // is specified AND control_materi_elasti_k0 is set to -yes, the poisson
  // coefficient consistent with K0, nu = K0/(1+K0) (from K0 = nu/(1-nu)),
  // is used in the elastic stress law with group_materi_elasti_young or
  // group_materi_elasti_young_power (replacing group_materi_elasti_poisson,
  // also when that record is absent) and with group_materi_elasti_hardsoil
  // (replacing both nu50 and nuur). For K0 > 0.95 Tochnog takes 0.95.
  // All three combinations measured on Professional (sigxx = K0*sigyy
  // under gravity, see manual-developer).
  // The switch is read for the index of the current timestep (ICONTROL);
  // a previous guard on index 0 skipped it for any other index (elasti7).
  k0_active = 0;
  k0_control_swit = 0;
  idum[0] = 0;
  db( ICONTROL, 0, idum, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_MATERI_ELASTI_K0, idum[0], &k0_control_swit, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  if ( k0_control_swit==-YES &&
       get_group_data( GROUP_MATERI_ELASTI_K0, gr, element, new_unknowns,
         &k0_elasti, ldum, GET_IF_EXISTS ) ) {
    if ( k0_elasti>0.95 ) k0_elasti = 0.95;
    k0_active = 1;
  }

  // group_materi_elasti_young_user -yes (manual Professional 6.662):
  // user_young gives the young modulus {input: the young record or 0}
  long int young_user = -NO;
  db( GROUP_MATERI_ELASTI_YOUNG_USER, gr, &young_user, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  if ( get_group_data( GROUP_MATERI_ELASTI_YOUNG, gr, element, new_unknowns, 
      &young, ldum, GET_IF_EXISTS ) || young_user==-YES ) {
    if ( young_user==-YES )
      user_young( new_unknowns, old_hisv, rotated_old_sig, young );
    get_group_data( GROUP_MATERI_ELASTI_POISSON, gr, element,
      new_unknowns, &poisson, ldum, GET_IF_EXISTS );
    if ( k0_active ) poisson = k0_elasti/(1.+k0_elasti);
    task[1] = -NO;
    C_matrix( young, poisson, elasti_transverse_isotropy, C, task );
    task[1] = membrane;
    C_matrix( young, poisson, elasti_transverse_isotropy, Cmem, task );
  }
  if ( get_group_data( GROUP_MATERI_ELASTI_YOUNG_POLYNOMIAL, gr, element, new_unknowns, 
      young_polynomial, length, GET_IF_EXISTS ) ) {
    strain_size = array_size( new_ept, MDIM*MDIM );
    young = 0.;
    for ( i=0; i<length; i++ )
      young += young_polynomial[i] * scalar_power(strain_size,i);
    get_group_data( GROUP_MATERI_ELASTI_POISSON, gr, element,
      new_unknowns, &poisson, ldum, GET_IF_EXISTS ); 
    task[1] = -NO;
    C_matrix( young, poisson, elasti_transverse_isotropy, C, task );
    task[1] = membrane;
    C_matrix( young, poisson, elasti_transverse_isotropy, Cmem, task );
  }
  if ( get_group_data( GROUP_MATERI_ELASTI_YOUNG_POWER, gr, element, new_unknowns, 
      young_power, ldum, GET_IF_EXISTS ) ) {
    // C_matrix ACCUMULATES into its target (array_add), so any elastic C
    // built before this block (e.g. group_materi_elasti_young) must be
    // cleared: the power law IS the Young modulus (Professional 6.662),
    // it does not add to the constant young. (The GNU 2014 code let both
    // records accumulate -> the stiffness was doubled; fixed here, see
    // manual-developer.)
    array_set( &C[0][0][0][0], 0., MDIM*MDIM*MDIM*MDIM );
    array_set( &Cmem[0][0][0][0], 0., MDIM*MDIM*MDIM*MDIM );
    get_group_data( GROUP_MATERI_ELASTI_POISSON, gr, element,
      new_unknowns, &poisson, ldum, GET_IF_EXISTS );
    if ( k0_active ) poisson = k0_elasti/(1.+k0_elasti);
    // manual Professional 6.662 / theory 2.2.2: E = E0 + E1*(p/p1)^alpha
    // with the conditions E >= E2 and E <= E3, where p is the pressure
    // (p = -(sig11+sig22+sig33)/3, positive in compression; compression
    // raises E). Parameters E0 E1 E2 E3 p1 alpha. NOTE: this is the
    // Professional 6-parameter convention — the GNU 3-parameter form
    // (young = young0*|p/p0|^alpha) is no longer accepted (see
    // manual-developer).
    young0 = young_power[0];
    young1 = young_power[1];
    young2 = young_power[2];
    young3 = young_power[3];
    p1 = young_power[4];
    alpha = young_power[5];
    if ( p1<=0. ) db_error( GROUP_MATERI_ELASTI_YOUNG_POWER, gr );
    p = - ( new_sig[0] + new_sig[4] + new_sig[8] ) / 3.;
    // materi_elasti_young_power_apply -no (manual Professional 6.801):
    // the power-law nonlinearity is ignored and the constant young from
    // the record (E0) is applied at all times.
    if ( control_materi_gate_off( CONTROL_MATERI_ELASTI_YOUNG_POWER_APPLY ) )
      young = young0;
    else {
      young = young0 + young1 * scalar_power( p/p1, alpha );
      if ( young<young2 ) young = young2;
      if ( young>young3 ) young = young3;
    }
    task[1] = -NO;
    C_matrix( young, poisson, elasti_transverse_isotropy, C, task );
    task[1] = membrane;
    C_matrix( young, poisson, elasti_transverse_isotropy, Cmem, task );
  }
  // group_materi_elasti_hardsoil (manual Professional 6.649 + theory
  // "Hardening-Soil model"): 7 parameters Eref_50 sigmaref_50 nu50 m
  // Eref_ur sigmaref_ur nuur. The Young modulus follows the power law
  // E = Eref * ((sig3 + c*cot(phi))/(sigmaref + c*cot(phi)))^m where
  // sig3 is the MINOR principal stress (the LEAST compressive one; the
  // manual orders sig3 > sig2 > sig1 with sig1 the largest compressive
  // stress). In this code (tension positive, compression negative) the
  // largest compressive stress is the SMALLEST algebraic eigenvalue, so
  // sig3_code = LARGEST algebraic eigenvalue and sig3_manual = -sig3_code
  // + c*cot(phi). NOTE: taking the SMALLEST eigenvalue would give the
  // AXIAL stress (sig1_manual) and the stiffness would depend on the
  // axial load instead of the confinement -- that is NOT the HS model
  // (verified analytically in the mhardsoil tests, see manual-developer).
  // First loading uses E50/nu50, unloading/reloading Eur/nuur. The
  // switch reads the maximum |p| history dof (initia
  // materi_plasti_hardsoil_history, manual 4.22 -- the shared sph dof,
  // updated in dof.cc): if the CURRENT estimated pressure of the step
  // (p_old + dp, dp the elastic pressure increment evaluated with the
  // first-loading tangent) is SMALLER than the maximum at step start,
  // the material is unloading/reloading -> Eur; otherwise first loading
  // -> E50. Same decision logic and one-step dp estimate as
  // group_materi_elasti_stress_pressure_history_factor (lot 6).
  if ( get_group_data( GROUP_MATERI_ELASTI_HARDSOIL, gr, element, new_unknowns,
      hs_elasti, ldum, GET_IF_EXISTS ) ) {
    // C_matrix ACCUMULATES into its target: the hardsoil law IS the
    // Young modulus, any elastic C built before must be cleared
    array_set( &C[0][0][0][0], 0., MDIM*MDIM*MDIM*MDIM );
    array_set( &Cmem[0][0][0][0], 0., MDIM*MDIM*MDIM*MDIM );
    hs_Eref50 = hs_elasti[0];
    hs_sigref50 = hs_elasti[1];
    nu50 = hs_elasti[2];
    hs_m = hs_elasti[3];
    hs_Erefur = hs_elasti[4];
    hs_sigrefur = hs_elasti[5];
    nuur = hs_elasti[6];
    if ( k0_active ) {
      // group_materi_elasti_k0 with control_materi_elasti_k0 -yes: the K0
      // consistent poisson replaces nu50 and nuur (manual 6.650)
      nu50 = k0_elasti/(1.+k0_elasti);
      nuur = nu50;
    }
    if ( hs_Eref50<=0. || hs_Erefur<=0. )
      db_error( GROUP_MATERI_ELASTI_HARDSOIL, gr );
    // cohesion term c*cot(phi) comes from the plastic group; without it
    // c = 0 (base = sig3_manual). cot(phi) is never evaluated for c = 0
    // (the term vanishes identically).
    hs_ccot = 0.;
    hs_plasti_present = 0;
    if ( get_group_data( GROUP_MATERI_PLASTI_HARDSOIL, gr, element,
        new_unknowns, hs_plasti, ldum, GET_IF_EXISTS ) ) {
      hs_plasti_present = 1;
      hs_phi = hs_plasti[0];
      hs_c = hs_plasti[1];
      if ( sin(hs_phi)<=0. || cos(hs_phi)<=0. )
        db_error( GROUP_MATERI_PLASTI_HARDSOIL, gr );
      if ( hs_c!=0. ) hs_ccot = hs_c * cos(hs_phi) / sin(hs_phi);
    }
    if ( hs_sigref50+hs_ccot<=0. || hs_sigrefur+hs_ccot<=0. )
      db_error( GROUP_MATERI_ELASTI_HARDSOIL, gr );
    // control_materi_plasti_hardsoil_gammap_initial (manual 6.146,
    // theory HS): with -yes tochnog creates an extra initial
    // contribution to gamma_p exactly such that the yield function is
    // zero-valued at the INITIAL stress state (convenient to start a
    // calculation with deviatoric stresses which would be outside the
    // yield surface without this contribution). Done in the first
    // timestep (the element record does not exist yet); the value is
    // stored in element_intpnt_materi_plasti_hardsoil_gammap_initial
    // (one value per element in this port; per integration point
    // pending, see manual-developer) and ADDED to gamma_p inside the
    // hardsoil yield function (plasti_rule reads the record). The
    // extra value is gamma_p = f(initial stress, gamma_p = 0) =
    // q/(E50*(1-q/qa)) - 2*q/Eur evaluated at the initial state.
    hs_gp_swit = -NO;
    db( ICONTROL, 0, &icontrol_hs, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( CONTROL_MATERI_PLASTI_HARDSOIL_GAMMAP_INITIAL, icontrol_hs, &hs_gp_swit,
      ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    if ( hs_gp_swit==-YES && hs_plasti_present ) {
      // first timestep detection: the record was pre-allocated (top.cc)
      // with the sentinel -1; the initialization overwrites it with
      // gamma_p_extra >= 0 (values < 0 mean "not initialized yet")
      hs_gp_extra = -1.;
      db( ELEMENT_INTPNT_MATERI_PLASTI_HARDSOIL_GAMMAP_INITIAL, element, idum,
        &hs_gp_extra, ldum, VERSION_NORMAL, GET_IF_EXISTS );
      if ( hs_gp_extra<0. ) {
      double hs_q0 = 0.;
      double sig0[MDIM*MDIM];
      array_move( new_sig, sig0, MDIM*MDIM );
      matrix_eigenvalues( sig0, sig_princ );
      hs_sig3 = sig_princ[0];
      for ( i=1; i<MDIM; i++ )
        if ( sig_princ[i]>hs_sig3 ) hs_sig3 = sig_princ[i];
      hs_base = -hs_sig3 + hs_ccot;
      if ( hs_base<=0. ) {
        hs_E50 = hs_Eref50; hs_Eur = hs_Erefur;   // clamp, documented
      }
      else {
        hs_E50 = hs_Eref50 * scalar_power( hs_base/(hs_sigref50+hs_ccot), hs_m );
        hs_Eur = hs_Erefur * scalar_power( hs_base/(hs_sigrefur+hs_ccot), hs_m );
      }
      if ( hs_m==0. ) { hs_E50 = hs_Eref50; hs_Eur = hs_Erefur; }
      hs_q0 = sqrt(
        0.5*( scalar_square(sig0[0]-sig0[4]) +
              scalar_square(sig0[4]-sig0[8]) +
              scalar_square(sig0[0]-sig0[8]) ) +
        3.*( scalar_square(sig0[1]) + scalar_square(sig0[2]) +
             scalar_square(sig0[5]) ) );
      hs_qf = 2. * sin(hs_phi) * hs_base / ( 1. - sin(hs_phi) );
      hs_qa = hs_qf / hs_plasti[3];
      if ( hs_qa<=0. ) db_error( GROUP_MATERI_PLASTI_HARDSOIL, gr );
      if ( hs_q0>=hs_qa ) {
        // initial deviatoric stress already beyond the asymptote: no
        // finite gamma_p brings f to zero (clamp, documented)
        hs_gp_extra = 1.e10;
      }
      else {
        hs_gp_extra = hs_q0 / ( hs_E50 * ( 1. - hs_q0/hs_qa ) )
          - 2. * hs_q0 / hs_Eur;
        if ( hs_gp_extra<0. ) hs_gp_extra = 0.;
      }
      db( ELEMENT_INTPNT_MATERI_PLASTI_HARDSOIL_GAMMAP_INITIAL, element, idum,
        &hs_gp_extra, hs_len, VERSION_NORMAL, PUT );
      }
    }
    // minor principal stress (sig3_manual = least compressive =
    // largest algebraic eigenvalue) at the step-start stress
    matrix_eigenvalues( new_sig, sig_princ );
    hs_sig3 = sig_princ[0];
    for ( i=1; i<MDIM; i++ )
      if ( sig_princ[i]>hs_sig3 ) hs_sig3 = sig_princ[i];
    hs_base = -hs_sig3 + hs_ccot;
    // base <= 0 (sig3 very tensile beyond the cohesion): clamp to
    // E = Eref (the power law is not evaluated on a non-positive base;
    // avoids zero/negative stiffness, documented in manual-developer)
    if ( hs_base<=0. ) {
      hs_E50 = hs_Eref50; hs_Eur = hs_Erefur;
    }
    else {
      hs_E50 = hs_Eref50 * scalar_power( hs_base/(hs_sigref50+hs_ccot), hs_m );
      hs_Eur = hs_Erefur * scalar_power( hs_base/(hs_sigrefur+hs_ccot), hs_m );
    }
    if ( hs_m==0. ) { hs_E50 = hs_Eref50; hs_Eur = hs_Erefur; }
    // trial first-loading C for the pressure-increment estimate of the
    // loading/unloading decision (matrix_a4b is NOT in-place safe:
    // work is a separate scratch, see lot 6)
    task[1] = -NO;
    C_matrix( hs_E50, nu50, elasti_transverse_isotropy, C, task );
    p = - ( new_sig[0] + new_sig[4] + new_sig[8] ) / 3.;
    {
      double dp_hs = 0.;
      matrix_a4b( C, inc_ept, work );
      dp_hs = -( work[0] + work[4] + work[8] ) / 3.;
      p += dp_hs;
    }
    // normalize IEEE -0.0 (scalar_dabs(-0.0) returns -0.0 and
    // -0.0 < sph is TRUE, which would take the Eur branch from the
    // very first load step, see lot 6)
    if ( p==0. ) p = 0.;
    if ( materi_plasti_hardsoil_history &&
         scalar_dabs(p) < old_unknowns[sph_indx] ) {
      // unloading/reloading
      young = hs_Eur;
      poisson = nuur;
    }
    else {
      // first loading
      young = hs_E50;
      poisson = nu50;
    }
    array_set( &C[0][0][0][0], 0., MDIM*MDIM*MDIM*MDIM );
    array_set( &Cmem[0][0][0][0], 0., MDIM*MDIM*MDIM*MDIM );
    task[1] = -NO;
    C_matrix( young, poisson, elasti_transverse_isotropy, C, task );
    task[1] = membrane;
    C_matrix( young, poisson, elasti_transverse_isotropy, Cmem, task );
  }
  if ( get_group_data( GROUP_MATERI_ELASTI_POISSON_POWER, gr, element, new_unknowns,
      poisson_power, ldum, GET_IF_EXISTS ) ) {
    get_group_data( GROUP_MATERI_ELASTI_POISSON, gr, element,
      new_unknowns, &poisson, ldum, GET_IF_EXISTS );
    // manual Professional 6.653 / theory 2.2.2: nu = nu0 + nu1*(p/p1)^alpha
    // with the condition nu <= nu2, where p is the pressure
    // (p = -(sig11+sig22+sig33)/3, positive in compression; same sign
    // convention as GROUP_MATERI_ELASTI_YOUNG_POWER).
    nu0 = poisson_power[0];
    nu1 = poisson_power[1];
    nu2 = poisson_power[2];
    p1 = poisson_power[3];
    alpha = poisson_power[4];
    if ( p1<=0. ) db_error( GROUP_MATERI_ELASTI_POISSON_POWER, gr );
    p = - ( new_sig[0] + new_sig[4] + new_sig[8] ) / 3.;
    poisson = nu0 + nu1 * scalar_power(scalar_dabs(p/p1),alpha);
    if ( poisson>nu2 ) poisson = nu2;
    // The power law REPLACES the poisson ratio of the elastic C built
    // above (manual 6.652: "the poisson ratio is made a function of the
    // average stress state"; Professional 02-08-2026 mpower = one
    // C(young, nu(p))). C_matrix ACCUMULATES into its target, so the C
    // of group_materi_elasti_young (or _young_power) is cleared first;
    // `young` keeps the modulus of that block. Without the clear the
    // stiffness was C(young,poisson) + C(young,nu(p)), about twice the
    // Professional one (audit K06, mpower sigyy -3.29 vs -1.41).
    array_set( &C[0][0][0][0], 0., MDIM*MDIM*MDIM*MDIM );
    array_set( &Cmem[0][0][0][0], 0., MDIM*MDIM*MDIM*MDIM );
    task[1] = -NO;
    C_matrix( young, poisson, elasti_transverse_isotropy, C, task );
    task[1] = membrane;
    C_matrix( young, poisson, elasti_transverse_isotropy, Cmem, task );
  }
  // group_materi_elasti_shear_factor (manual Professional 6.654): the
  // shear stiffness following from a specified young and poisson is
  // multiplied with factor (convenient to test the effect of low shear
  // stiffness). Scales ONLY the shear entries of C and Cmem — the terms
  // C[i][j][k][l] that connect shear strain (k!=l) to shear stress
  // (i!=j); in C_matrix those are the (0,1),(0,2),(1,2) shear blocks
  // (diagonal (3,3),(4,4),(5,5) entries of the Voigt matrix in 3D, and
  // the (2,2) shear entry in 2D plane strain/plane stress).
  if ( get_group_data( GROUP_MATERI_ELASTI_SHEAR_FACTOR, gr, element, 
      new_unknowns, &shear_factor, ldum, GET_IF_EXISTS ) ) {
    for ( idim=0; idim<MDIM; idim++ ) {
      for ( jdim=0; jdim<MDIM; jdim++ ) {
        if ( idim==jdim ) continue;
        for ( kdim=0; kdim<MDIM; kdim++ ) {
          for ( ldim=0; ldim<MDIM; ldim++ ) {
            if ( kdim==ldim ) continue;
            C[idim][jdim][kdim][ldim] *= shear_factor;
            Cmem[idim][jdim][kdim][ldim] *= shear_factor;
          }
        }
      }
    }
  }
  // group_materi_elasti_stress_pressure_history_factor (manual
  // Professional 6.655): models a different soil stiffness on first
  // loading versus unloading/reloading. Requires the initia
  // materi_stress_pressure_history (manual 4.50), which stores the
  // maximum of the absolute value of the pressure over time in the
  // node_dof records (updated in dof.cc, parallel_new_dof_diagonal).
  // If the CURRENT pressure |p| (p = -sig_mean, positive in compression)
  // is SMALLER than the largest pressure in history (the interpolated
  // dof new_unknowns[sph_indx]), the material is unloading/reloading and
  // the elastic stiffness is multiplied with factor; if the current
  // pressure is the new maximum, it becomes the maximum history pressure
  // (in dof.cc) and the stiffness is NOT multiplied. Point of
  // application: the final elastic C/Cmem (after shear_factor), so the
  // factor scales the complete stiffness built from group_materi_elasti_young
  // or group_materi_elasti_young_power and combines with
  // group_materi_elasti_poisson_power/shear_factor (see manual-developer).
  if ( materi_stress_pressure_history &&
       get_group_data( GROUP_MATERI_ELASTI_STRESS_PRESSURE_HISTORY_FACTOR,
         gr, element, new_unknowns, &sph_factor, ldum, GET_IF_EXISTS ) ) {
    // The elastic stiffness is evaluated with the PREVIOUS step's stress
    // (new_sig is initialized from the old unknowns), so the pressure for
    // the unloading/reloading decision is the pressure the CURRENT step
    // will reach: p_new = p_old + dp, dp = -mean(C:inc_ept) (the elastic
    // pressure increment of this step). This catches the FIRST unloading
    // step (with p_old only, the decision lags one step: p_old == sph at
    // the peak, so the factor would start one step late). The history
    // maximum is read from the OLD unknowns (VERSION_NORMAL, the value at
    // the START of the step): during the step the sph dof is raised by
    // dof.cc (parallel_new_dof_diagonal) with the running |p|, and the
    // decision must compare against the history EXCLUDING the current
    // step — otherwise p_est == sph exactly at the peak and any rounding
    // decides loading vs unloading (spurious factor, runaway history).
    p = - ( new_sig[0] + new_sig[4] + new_sig[8] ) / 3.;
    {
      double dp=0.;
      // work = C:inc_ept (elastic stress increment of this step); the
      // source inc_ept is untouched (matrix_a4b must not be in-place)
      matrix_a4b( C, inc_ept, work );
      dp = -( work[0] + work[4] + work[8] ) / 3.;
      // K04 stress_pressure_factor (2026-10-06), rule measured on the
      // Professional 02-08-2026 (Pro test stress_pressure_factor at dt
      // 0.1 / 0.15 / 0.3 / 0.4): a step that STARTS below the history
      // estimates its end pressure with the FACTORED (unloading/
      // reloading) stiffness, the stiffness it uses while it stays
      // below the history; it stays unloading/reloading only when that
      // estimate is strictly below the history (a step that would end
      // exactly at the history, up to round-off, is first loading).
      // Before: the first-loading dp, so a reloading step switched to
      // first loading one step early (dt 0.1: final sigxx -3.10 vs
      // -3.05; dt 0.3: -2.45 vs -2.3 at disx -2.2).
      if ( scalar_dabs(p) < old_unknowns[sph_indx] ) dp *= sph_factor;
      p += dp;
    }
    // normalize IEEE -0.0: scalar_dabs(-0.0) returns -0.0 (the branch
    // a<0. is false for -0.0), and -0.0 < sph is TRUE for any sph,
    // which would apply the factor from the very first load step.
    if ( p==0. ) p = 0.;
    if ( scalar_dabs(p) < old_unknowns[sph_indx] *
         ( 1. - 1.e-10 ) ) {
      for ( idim=0; idim<MDIM; idim++ )
        for ( jdim=0; jdim<MDIM; jdim++ )
          for ( kdim=0; kdim<MDIM; kdim++ )
            for ( ldim=0; ldim<MDIM; ldim++ ) {
              C[idim][jdim][kdim][ldim] *= sph_factor;
              Cmem[idim][jdim][kdim][ldim] *= sph_factor;
            }
    }
  }
  if ( get_group_data( GROUP_MATERI_ELASTI_YOUNG_STRAINSTRESS, gr, element, new_unknowns, 
      young_strainstress, length, GET_IF_EXISTS ) ) {
    
    meanstrain = ( new_ept[0] + new_ept[4] + new_ept[8] ) / 3.;
    array_move( new_ept, ept_dev, MDIM*MDIM );
    for ( idim=0; idim<MDIM; idim++ ) ept_dev[idim*MDIM+idim] -= meanstrain;
    straindev_size = array_size( ept_dev, MDIM*MDIM ); 

    young = 0.;

    linear_ept=young_strainstress[2];
    zero_ept=young_strainstress[3]; 	
    
    for ( i=5; i<length; i++ )
    		young_linear_ept += young_strainstress[i] * scalar_power(linear_ept,i-5);	    

    if (straindev_size <= linear_ept)
    	for ( i=5; i<length; i++ )
    		young += young_strainstress[i] * scalar_power(straindev_size,i-5);
    else if (straindev_size >= linear_ept && straindev_size <= zero_ept) 
	young = young_linear_ept*(log10(zero_ept) - 
	log10(straindev_size))/(log10(zero_ept) - log10(linear_ept));
    else young = 0;	
	
    alpha = young_strainstress[0];
    p0 = young_strainstress[1];

    if ( p0<=0. ) db_error( GROUP_MATERI_ELASTI_YOUNG_STRAINSTRESS, gr );
    p = - ( new_sig[0] + new_sig[4] + new_sig[8] ) / 3.;
    young = young * scalar_power(scalar_dabs(p/p0),alpha);

    if (new_kappa>0) 
        young=young_strainstress[4]* scalar_power(scalar_dabs(p/p0),alpha);

    get_group_data( GROUP_MATERI_ELASTI_POISSON, gr, element,
      new_unknowns, &poisson, ldum, GET_IF_EXISTS ); 
    task[1] = -NO;
    C_matrix( young, poisson, elasti_transverse_isotropy, C, task );
    task[1] = membrane;
    C_matrix( young, poisson, elasti_transverse_isotropy, Cmem, task );
  }
  if (
      db( GROUP_MATERI_ELASTI_VOLUMETRIC_YOUNG_ORDER, gr,
        &volumetric_young_order, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS ) ||
      db( GROUP_MATERI_ELASTI_VOLUMETRIC_YOUNG_VALUES, gr, idum,
        volumetric_young_values, ldum, VERSION_NORMAL, GET_IF_EXISTS )) {
    get_group_data( GROUP_MATERI_ELASTI_VOLUMETRIC_POISSON, gr, element,
      new_unknowns, &poisson, ldum, GET_IF_EXISTS );
    get_group_data( GROUP_MATERI_ELASTI_VOLUMETRIC_YOUNG_VALUES, gr, element,
      new_unknowns, volumetric_young_values, length, GET );
    if ( length<4 ) db_error( GROUP_MATERI_ELASTI_VOLUMETRIC_YOUNG_VALUES, gr );
    if ( !fit_polynomial( volumetric_young_values, length/2, young_polynomial,
        volumetric_young_order ) ) {
      pri( "Error detected for GROUP_MATERI_ELASTI_VOLUMETRIC_YOUNG_*" );
      pri( "Cannot determine polynomial" );
      exit(TN_EXIT_STATUS);
    }
    strain_size = array_size( new_ept, MDIM*MDIM );
    if ( strain_size>volumetric_young_values[length-2] )
      strain_size = volumetric_young_values[length-2];
    young = 0.;
    for ( i=1; i<volumetric_young_order; i++ ) {
      young += i * young_polynomial[i] * scalar_power(strain_size,i-1) *
       (1.+poisson)*(1.-2.*poisson) / ( 1.-poisson);
    }
    task[1] = -NO;
    C_matrix( young, poisson, elasti_transverse_isotropy, C, task );
    task[1] = membrane;
    C_matrix( young, poisson, elasti_transverse_isotropy, Cmem, task );
  }                                  
  if ( get_group_data( GROUP_MATERI_ELASTI_CAMCLAY_G, gr, element, new_unknowns, 
      camclay, ldum, GET_IF_EXISTS ) || 
      get_group_data( GROUP_MATERI_ELASTI_CAMCLAY_POISSON, 
      gr, element, new_unknowns, camclay, ldum, GET_IF_EXISTS ) ) {
    if(!get_group_data( GROUP_MATERI_PLASTI_CAMCLAY, gr, element, 
      new_unknowns, work, ldum, GET_IF_EXISTS )){
      get_group_data( GROUP_MATERI_PLASTI_CAMCLAY_INCREMENTAL, gr, element, 
      new_unknowns, work, ldum, GET);
    }  
    kappa = work[1];
    e = old_hisv[0];
    pressure = - ( new_sig[0] + new_sig[4] + new_sig[8] ) / 3.;
    // group_materi_elasti_camclay_pressure_min (manual Professional
    // 6.647): minimal allowed value for the pressure in the camclay
    // bulk modulus. Pressures below pressure_min are set to
    // pressure_min, preventing numerical problems for very low (or
    // negative/tensile) bulk modulus values K = (1+e)*p/kappa.
    if ( get_group_data( GROUP_MATERI_ELASTI_CAMCLAY_PRESSURE_MIN, gr,
        element, new_unknowns, camclay, ldum, GET_IF_EXISTS ) ) {
      if ( pressure<camclay[0] ) pressure = camclay[0];
    }
    k = (1.+e)*pressure/kappa;
    if      ( get_group_data( GROUP_MATERI_ELASTI_CAMCLAY_G, gr, element, 
        new_unknowns, camclay, ldum, GET_IF_EXISTS ) ) {
      g = camclay[0];
      poisson = (3.*k-2.*g)/(2.*g+6.*k);
    }
    else if ( get_group_data( GROUP_MATERI_ELASTI_CAMCLAY_POISSON, gr, element, 
        new_unknowns, camclay, ldum, GET_IF_EXISTS ) ) {
      poisson = camclay[0];
      g = (3./2.) * k * ( 1 - 2.*poisson ) / ( 1. + poisson );
    }
    if ( scalar_dabs(poisson)==1. ) {
      pri( "\nError detected for camclay plasticity." );
      pri( "Be sure to  specify initial stresses in the NODE_DOF records." );
      pri( "Specify legal camclay data. (Or maybe your calculations diverged)." );
      exit(1);
    }
    young = 2.*g*(1.+poisson);
    task[1] = -NO;
    C_matrix( young, poisson, elasti_transverse_isotropy, C, task );
    task[1] = membrane;
    C_matrix( young, poisson, elasti_transverse_isotropy, Cmem, task );
  }
  if ( get_group_data( GROUP_MATERI_ELASTI_LADE, gr, element,
      new_unknowns, group_materi_elasti_lade, ldum, GET_IF_EXISTS ) ) {
    lade_1 = group_materi_elasti_lade[0];
    lade_2 = group_materi_elasti_lade[1];
    lade_3 = group_materi_elasti_lade[2];
    pressure = ( new_sig[0] + new_sig[4] + new_sig[8] ) / 3.;
    array_move( new_sig, sig_dev, MDIM*MDIM );
    for ( idim=0; idim<MDIM; idim++ ) sig_dev[idim*MDIM+idim] -= pressure;
    C_matrix_lade( lade_1, lade_2, lade_3, pressure, sig_dev, C );
    array_move( &C[0][0][0][0], &Cmem[0][0][0][0], MDIM*MDIM*MDIM*MDIM );
    if ( membrane==-YES ) {
      cout << "\nError: GROUP_MATERI_ELASTI_LADE not available for membrane stress state.\n";
      exit(TN_EXIT_STATUS);
    }
    if ( get_group_data( GROUP_MATERI_ELASTI_POISSON, gr, element,
        new_unknowns, &poisson, ldum, GET_IF_EXISTS ) ) {
      pri( "GROUP_MATERI_ELASTI_POISSON cannot be used with GROUP_MATERI_ELASTI_LADE." );
      exit(TN_EXIT_STATUS);
    }
  }

/***********************elastic models added**********************/

  if ( get_group_data( GROUP_MATERI_ELASTI_TSKH, gr, element, new_unknowns, 
      tskh, length, GET_IF_EXISTS ) ) {
      
    if(!get_group_data( GROUP_MATERI_PLASTI_TSKH, gr, element, 
      new_unknowns, work, ldum, GET_IF_EXISTS )) {
      get_group_data( GROUP_MATERI_PLASTI_AITSKH, gr, element, 
      new_unknowns, work, ldum, GET );
    }

    double A = tskh[0];
    double n= tskh[1];
    double m_elasti = tskh[2];
    kappa = work[1];
    e = old_hisv[0]-1;
    double p0 = old_hisv[1];
    if(p0<=0) {
    	cout<<"error in tskh_elsti"<<endl;
	exit(1);
    }
    
    double geo_sigm= - ( new_sig[0] + new_sig[4] + new_sig[8] ) / 3.;
    if(geo_sigm<=0) geo_sigm=TINY;
    double R0=2*p0 / geo_sigm;
    double g = A*scalar_power(geo_sigm, n)*scalar_power(R0, m_elasti);
    double k = geo_sigm/kappa;

    if ( task[0]==GROUP_MATERI_ISOTROPY ) {
      poisson = (3*k - 2*g)/(2*g + 6*k);
      if(poisson>0.4999) poisson=0.4999;
      else if(poisson<0.0001) poisson=0.0001;
      young = 2.*g*(1.+poisson);
    }
    else if ( task[0]==GROUP_MATERI_ELASTI_TRANSVERSE_ISOTROPY_GRAHOUL ) {
      double alpha=elasti_transverse_isotropy[0];
      poisson = (k*(3+3*alpha*alpha/2)-g*(1+2*alpha*alpha))/(k*(3+6*alpha)+g*(4*alpha-1));
      if(poisson>0.4999) poisson=0.4999;
      else if(poisson<0.0001) poisson=0.0001;
      young = 3*g*(1+poisson)*(1-2*poisson)/(1-poisson-2*alpha*poisson+alpha*alpha/2);
    }
    task[1] = -NO;
    C_matrix( young, poisson, elasti_transverse_isotropy, C, task );
    task[1] = membrane;
    C_matrix( young, poisson, elasti_transverse_isotropy, Cmem, task );
  }
  if ( get_group_data( GROUP_MATERI_ELASTI_SMALLSTRAIN, gr, element, new_unknowns, 
      smallstrain, length, GET_IF_EXISTS ) ) {
    double A=smallstrain[0];	
    double n=smallstrain[1];	
    double kappa=smallstrain[2];	
    double R=smallstrain[3];	
    double gamma=smallstrain[4];	
    double delta=smallstrain[5];	
	
    double geo_sigm= - ( new_sig[0] + new_sig[4] + new_sig[8] ) / 3.;
    if(geo_sigm<=0) geo_sigm=TINY;

    double geoept[MDIM*MDIM];
    array_set(geoept, 0, MDIM*MDIM);
    array_multiply( new_ept, geoept, -1, MDIM*MDIM );
    double geoeptsJ=0, geoepts=0, geoeptv=0;
    calc_IJlode(geoept, geoeptv, geoeptsJ, rdum, false, ddumarray, ddumarray, ddumarray);
    geoepts=geoeptsJ*2/sqrt(3.);

    double G0 = A*scalar_power(geo_sigm, n);
    double K0 = geo_sigm/kappa;

    if(geoepts<R) geoepts = R;
    if(scalar_dabs(geoeptv) < R) geoeptv = R;	
    double g=G0*scalar_power((R/geoepts),gamma);
    double k=K0*scalar_power((R/scalar_dabs(geoeptv)),delta);
	
    double ct=100;	
    if ( task[0]==GROUP_MATERI_ISOTROPY ) {
      poisson = (3*k - 2*g)/(2*g + 6*k);
      if(poisson>0.4999) poisson=0.4999;
      else if(poisson<0.0001) poisson=0.0001;
      young = 2.*g*(1.+poisson);

      double poisson0 = (3*K0 - 2*G0)/(2*G0 + 6*K0);
      double young0 = 2.*G0*(1.+poisson0);
      if(young<young0/ct) young=young0/ct;
    }
    else if ( task[0]==GROUP_MATERI_ELASTI_TRANSVERSE_ISOTROPY_GRAHOUL ) {
      double alpha=elasti_transverse_isotropy[0];
      poisson = (k*(3+3*alpha*alpha/2)-g*(1+2*alpha*alpha))/(k*(3+6*alpha)+g*(4*alpha-1));
      if(poisson>0.4999) poisson=0.4999;
      else if(poisson<0.0001) poisson=0.0001;
      young = 3*g*(1+poisson)*(1-2*poisson)/(1-poisson-2*alpha*poisson+alpha*alpha/2);

      double poisson0 = (K0*(3+3*alpha*alpha/2)-G0*(1+2*alpha*alpha))/(K0*(3+6*alpha)+G0*(4*alpha-1));
      double young0 = 3*G0*(1+poisson0)*(1-2*poisson0)/(1-poisson0-2*alpha*poisson0+alpha*alpha/2);
      if(young<young0/ct) young=young0/ct;
    }
    task[1] = -NO;
    C_matrix( young, poisson, elasti_transverse_isotropy, C, task );
    task[1] = membrane;
    C_matrix( young, poisson, elasti_transverse_isotropy, Cmem, task );
    
    if ( materi_history_variables<2 ) {
      pri( "Error: materi_history_variables should be 2 for GROUP_MATERI_ELASTI_SMALLSTRAIN (G and K for postprocessing)" );
      exit(TN_EXIT_STATUS);
    }
    new_hisv[0]=g;	//store for post-processing
    new_hisv[1]=k;
  }

/*****************************************************************************/

  // group_materi_elasti_c (manual Professional 6.642): the 81 values of
  // C_ijkl row by row {xx, xy, xz, yx, yy, yz, zx, zy, zz}, tensor strains;
  // added to the other elastic laws of the group, as the Professional
  // 02-08-2026. group_materi_elasti_c_direction (6.643): the local axes
  // dir_0 dir_1 dir_2 {global components} in which C is given, so
  // C_ijkl = dir_pi dir_qj dir_rk dir_sl C'_pqrs {unit dir}. Added after
  // every other elastic law, because several of them clear C first. No plane stress
  // condensation: refused with a membrane.
  {
    double c_local[MDIM*MDIM*MDIM*MDIM], c_dir[MDIM*MDIM];
    if ( get_group_data( GROUP_MATERI_ELASTI_C, gr, element, new_unknowns,
        c_local, ldum, GET_IF_EXISTS ) ) {
      if ( membrane==-YES ) {
        pri( "Error: group_materi_elasti_c is not available with "
          "group_materi_membrane -yes." );
        exit(1);
      }
      array_set( c_dir, 0., MDIM*MDIM );
      c_dir[0] = c_dir[4] = c_dir[8] = 1.;
      if ( get_group_data( GROUP_MATERI_ELASTI_C_DIRECTION, gr, element,
          new_unknowns, c_dir, ldum, GET_IF_EXISTS ) ) {
        // each axis is normalized, as the Professional 02-08-2026
        for ( idim=0; idim<MDIM; idim++ ) {
          if ( array_null( &c_dir[idim*MDIM], MDIM ) ) {
            pri( "Error: group_materi_elasti_c_direction has a zero "
              "direction." );
            exit(1);
          }
          array_normalize( &c_dir[idim*MDIM], MDIM );
        }
      }
      for ( idim=0; idim<MDIM; idim++ ) {
        for ( jdim=0; jdim<MDIM; jdim++ ) {
          for ( kdim=0; kdim<MDIM; kdim++ ) {
            for ( ldim=0; ldim<MDIM; ldim++ ) {
              double c_ijkl = 0.;
              long int ip=0, iq=0, ir=0, is=0;
              for ( ip=0; ip<MDIM; ip++ )
                for ( iq=0; iq<MDIM; iq++ )
                  for ( ir=0; ir<MDIM; ir++ )
                    for ( is=0; is<MDIM; is++ )
                      c_ijkl += c_dir[ip*MDIM+idim] * c_dir[iq*MDIM+jdim] *
                        c_dir[ir*MDIM+kdim] * c_dir[is*MDIM+ldim] *
                        c_local[((ip*MDIM+iq)*MDIM+ir)*MDIM+is];
              C[idim][jdim][kdim][ldim] += c_ijkl;
              Cmem[idim][jdim][kdim][ldim] += c_ijkl;
            }
          }
        }
      }
    }
  }

  if ( memory==-TOTAL || memory==-TOTAL_PIOLA ) {
    if ( !materi_strain_total  ) formulation = TOTAL;
  }
  if ( memory==-TOTAL || memory==-TOTAL_PIOLA  || memory==-TOTAL_LINEAR ) {
    // -total_linear (manual Professional 2.2.1): the stress at the new time
    // point is the previous stress plus the extra stress from the
    // INCREMENTAL strain, so formulation keeps its INCREMENTAL default.
    // Forcing TOTAL here recomputed the stress from the total strain and
    // discarded a stress set by control_reset_dof (delete3/mpc5). The
    // incremental strain increment is computed in set_deften_etc from the
    // old/new total deformation gradients (materi.cc).
    check_unknown( "materi_velocity", YES, CHECK_USAGE_AND_ERROR );
    check_unknown( "materi_displacement", YES, CHECK_USAGE_AND_ERROR );
    check_unknown( "materi_stress", YES, CHECK_USAGE_AND_ERROR );
  }
  if ( memory==-UPDATED || memory==-UPDATED_WITHOUT_ROTATION ) {
    check_unknown( "materi_velocity", YES, CHECK_USAGE_AND_ERROR );
    // The displacement restriction applies only when the updated
    // formulation is EXPLICITLY requested: with the default memory
    // (record absent) the Professional accepts materi_displacement
    // together with materi_velocity (e.g. the Masin clay corpus tests
    // hypo12/13 run without group_materi_memory).
    if ( db_active_index( GROUP_MATERI_MEMORY, gr, VERSION_NORMAL ) )
      check_unknown( "materi_displacement", NO, CHECK_USAGE_AND_ERROR );
    check_unknown( "materi_stress", YES, CHECK_USAGE_AND_ERROR );
  }


  db( GROUP_MATERI_PLASTI_KINEMATIC_HARDENING, gr, idum, &plasti_kinematic_hardening, 
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  get_group_data( GROUP_USER_DATA, gr, element, new_unknowns, user_data, 
    nuser_data, GET_IF_EXISTS );

  if ( condif_temperature ) {
    array_set( inc_temperature_strain, 0., MDIM*MDIM );
    array_set( new_temperature_strain, 0., MDIM*MDIM );
    get_group_data( GROUP_MATERI_EXPANSION_LINEAR, gr, element, new_unknowns, 
      &materi_expansion_linear, ldum, GET_IF_EXISTS );
    for ( idim=0; idim<MDIM; idim++ ) {
      inc_temperature_strain[idim*MDIM+idim] = -materi_expansion_linear * 
        ( new_unknowns[temp_indx] - old_unknowns[temp_indx] ) ;
      new_temperature_strain[idim*MDIM+idim] = -materi_expansion_linear * 
        new_unknowns[temp_indx];
    }
    if ( swit ) {
      pri( "inc_temperature_strain", inc_temperature_strain, MDIM, MDIM );
      pri( "new_temperature_strain", new_temperature_strain, MDIM, MDIM );
    }
  }

  array_set( ddsdde, 0., MSTRAIN*MSTRAIN );
  array_set( &Cuser[0][0][0][0], 0., MDIM*MDIM*MDIM*MDIM );
  array_set( &Chyper[0][0][0][0], 0., MDIM*MDIM*MDIM*MDIM );
  array_set( &Chypo[0][0][0][0], 0., MDIM*MDIM*MDIM*MDIM );
  array_set( plasti_dir, 0., MDIM*MDIM );
  array_set( total_inc_epp, 0., MDIM*MDIM );
  deps_size = array_size( inc_ept, MDIM*MDIM );
  if ( deps_size<EPS_SIZE ) deps_size = EPS_SIZE;

  if ( materi_plasti_f_nonlocal )
    max_iter_visco = 2;
  else
    max_iter_visco = 10;
  max_iter_ep = MAX_PLASTI_ITER;

  // user-specified maximum number of plastic iterations on
  // integration point level: group_materi_plasti_maximum_iterations,
  // else the global materi_plasti_maximum_iterations (manual 6.804)
  {
    long int max_plasti_iter_user=0, length_it=0;
    if ( db( GROUP_MATERI_PLASTI_MAXIMUM_ITERATIONS, gr, &max_plasti_iter_user,
        ddum, length_it, VERSION_NORMAL, GET_IF_EXISTS ) ||
         db( MATERI_PLASTI_MAXIMUM_ITERATIONS, 0, &max_plasti_iter_user,
        ddum, length_it, VERSION_NORMAL, GET_IF_EXISTS ) ) {
      if ( max_plasti_iter_user<1 ) {
        pri( "Error: group_materi_plasti_maximum_iterations / "
          "materi_plasti_maximum_iterations must be at least 1." );
        exit(1);
      }
      max_iter_visco = max_iter_ep = max_plasti_iter_user;
    }
  }
  max_plasti_iter = viscoplasti ? max_iter_visco : max_iter_ep;

  // plasti_apply (manual Professional 6.899, global) and
  // control_plasti_apply (6.257, per timestep index) -no: any plasticity
  // data is ignored for these timesteps (the control record of the
  // current timestep overrules the global one). Scope measured on the
  // Professional 02-08-2026 (see manual-developer/plasti_apply.md): the
  // yield-surface models evaluated by plasti_rule (mohr_coul, vonmises,
  // druck_prag, tension, camclay, hardsoil, ...), the direct cut-offs
  // (mohr_coul_direct, tension_direct, compression_direct) and the
  // incremental models (tskh, aitskh, camclay_incremental) fall back to
  // the elastic response; plastic strains of earlier timesteps are kept.
  plasti_off = materi_apply_gate_off( PLASTI_APPLY, CONTROL_PLASTI_APPLY );
  // group_plasti_apply -no (manual Professional 6.762): any plasticity
  // data of this group is neglected; a group -yes does not overrule a
  // global -no (Professional 02-08-2026)
  {
    long int group_plasti_apply=-YES;
    db( GROUP_PLASTI_APPLY, gr, &group_plasti_apply, ddum, ldum,
      VERSION_NORMAL, GET_IF_EXISTS );
    if ( group_plasti_apply==-NO ) plasti_off = 1;
  }

  /*********added for explicit time integration*******************************************/

  array_move(old_epi, new_epi, MDIM*MDIM);

  plasti_type = -NONE;
  //double tmpdoub=0;		
  bool plasti_incremental=false;
  double plasti_dt[DATA_ITEM_SIZE];
  for(i=0; i<DATA_ITEM_SIZE; i++) plasti_dt[i]=0;
  long int length_pl=0, i_points=0, tmp2=0;

  if(get_group_data( GROUP_MATERI_PLASTI_TSKH, 
  	gr, element, new_unknowns, plasti_dt, length_pl, GET_IF_EXISTS)) {
  	plasti_type=GROUP_MATERI_PLASTI_TSKH;
	plasti_incremental=true;
  }	
  else if(get_group_data( GROUP_MATERI_PLASTI_AITSKH, 
  	gr, element, new_unknowns, plasti_dt, length_pl, GET_IF_EXISTS)) {
  	plasti_type=GROUP_MATERI_PLASTI_AITSKH;
	plasti_incremental=true;
  }	
  else if(get_group_data( GROUP_MATERI_PLASTI_CAMCLAY_INCREMENTAL, 
  	gr, element, new_unknowns, plasti_dt, length_pl, GET_IF_EXISTS)) {
  	plasti_type=GROUP_MATERI_PLASTI_CAMCLAY_INCREMENTAL;
	plasti_incremental=true;
  }	
  if ( plasti_off ) plasti_incremental = false;
  if(plasti_incremental) {
        if (db_active_index( GROUP_INTEGRATION_POINTS, gr, VERSION_NORMAL) ) 
	    db( GROUP_INTEGRATION_POINTS, gr, &i_points, &tmp, tmp2, VERSION_NORMAL, GET );
	if((i_points!=-MAXIMAL) && (options_element_dof != -YES)) {
	    cout<<"You must set 'options_element_dof -yes'"<<endl; 
	    cout<<"for explicit time integration"<<endl<<endl;
	    exit(1);
 	}

  	//interface routine to single element program triax (D.Masin@city.ac.uk)
        plasti_incr(
		/*input*/ rotated_old_sig, inc_ept, old_epp,
			new_ept, C, old_hisv, plasti_type,
		 	plasti_dt, length_pl, gr, element, softvar_nonl, softvar_l,
     		/*output*/ new_sig, inc_epp, new_hisv, new_f, new_substeps, Cmem
	);
	array_add( inc_epp, old_epp, new_epp, MDIM*MDIM );
	array_subtract( inc_ept, inc_epp, inc_epe, MDIM*MDIM );
	array_add( inc_epe, old_epe, new_epe, MDIM*MDIM );
	
        if ( materi_plasti_kappa ) {
	  tmp = array_inproduct( inc_epp, inc_epp, MDIM*MDIM );
	  new_kappa = old_kappa + sqrt(0.5*tmp);     
	  if ( swit ) pri( "new_kappa", new_kappa );
	}
        // materi_plasti_kappa_shear (manual 4.25): the size of the SHEAR
        // plastic strain rate kappa_shear = int sqrt(0.5*dev(eps_p):dev(eps_p)).
        if ( materi_plasti_kappa_shear ) {
          double epp_dev[MDIM*MDIM];
          array_move( inc_epp, epp_dev, MDIM*MDIM );
          {
            double mean_epp = 0.;
            for ( int idim=0; idim<MDIM; idim++ ) mean_epp += epp_dev[idim*MDIM+idim];
            mean_epp /= MDIM;
            for ( int idim=0; idim<MDIM; idim++ ) epp_dev[idim*MDIM+idim] -= mean_epp;
          }
          tmp = array_inproduct( epp_dev, epp_dev, MDIM*MDIM );
          new_kapsh = old_kapsh + sqrt(0.5*tmp);
          if ( swit ) pri( "new_kapsh", new_kapsh );
        }
        else
          new_kapsh = old_kapsh;
        // cap1 hardening (materi_plasti_cap1_history, manual theory
        // cap1): pc hardens with the cap plastic volume strain
        // eps_p_cv_dot = (lambda*/kappa* - 1)/K_ref (p_ref/p*c)^m pc_dot,
        // inverted and integrated explicitly with the OLD pc:
        // new_pc = old_pc + deps_p_cv*K_ref/(lambda*/kappa* - 1)
        //                     *((old_pc + c*cot(phi))/p_ref)^m
        // with deps_p_cv = -trace(inc_epp) (positive in compression,
        // clamped >= 0: unloading never decreases pc).
        if ( materi_plasti_cap1_history ) {
          double cap1_plasti_data[8], ccotphi = 0., deps_p_cv = 0.,
            lambda_star = 0., kappa_star = 0., K_ref = 0., p_ref = 0.,
            m_cap1 = 0.;
          long int cap1_length = 0;
          new_cap1pc = old_cap1pc;
          if ( get_group_data( GROUP_MATERI_PLASTI_CAP1, gr, element,
              new_unknowns, cap1_plasti_data, cap1_length,
              GET_IF_EXISTS ) ) {
            lambda_star = cap1_plasti_data[3];
            kappa_star = cap1_plasti_data[4];
            K_ref = cap1_plasti_data[5];
            p_ref = cap1_plasti_data[6];
            m_cap1 = cap1_plasti_data[7];
            ccotphi = cap1_plasti_data[0] * cos(cap1_plasti_data[1]) /
              sin(cap1_plasti_data[1]);
            deps_p_cv = - ( inc_epp[0] + inc_epp[4] + inc_epp[8] );
            if ( deps_p_cv<0. ) deps_p_cv = 0.;
            if ( lambda_star>kappa_star && K_ref>0. && p_ref>0. ) {
              tmp = old_cap1pc + ccotphi;
              new_cap1pc = old_cap1pc + deps_p_cv * K_ref /
                ( lambda_star/kappa_star - 1. ) *
                scalar_power( tmp/p_ref, m_cap1 );
            }
          }
          if ( swit ) pri( "new_cap1pc", new_cap1pc );
        }
  }
  else {

  /*********************************************************************************/

    // plastic iterations
    // History dofs seen by the yield rules (F4 T5 / C35 follow-up,
    // 2026-10-10): plasti_rule() reads the hardening variables kappa and
    // kappa_shear from its unknowns array (hardsoil gamma_p, the kappa
    // hardened mohr_coul variants, von Mises Nadai). new_unknowns holds
    // the NODAL value written by the previous material evaluation of the
    // same step, so the result depended on the number of evaluations
    // (period-2 oscillation of mhardsoil_plast). The return mapping must
    // use the step-start value plus its own increment: pl_unknowns is
    // new_unknowns with those two slots replaced by old + increment.
  double pl_unknowns[MUKNWN];
  array_move( new_unknowns, pl_unknowns, nuknwn );
  if ( materi_plasti_kappa ) pl_unknowns[kap_indx] = old_kappa;
  if ( materi_plasti_kappa_shear ) pl_unknowns[kapsh_indx] = old_kapsh;
  while ( !plasti_found && total_plasti_iter<max_plasti_iter ) {
    plasti_iter++;
    total_plasti_iter++;
    if ( swit ) {
      pri( "plasti_iter", plasti_iter );
      pri( "total_plasti_iter", total_plasti_iter );
      pri( "lambda", lambda );
    }
      // plastic strain part
    if ( viscoplasti )
      array_multiply( plasti_dir, work_inc_epp, lambda*dtime, MDIM*MDIM );
    else
      array_multiply( plasti_dir, work_inc_epp, lambda*deps_size, MDIM*MDIM );
    array_add( work_inc_epp, total_inc_epp, inc_epp, MDIM*MDIM );
    if ( swit ) pri( "plastic strain increment", inc_epp, MDIM, MDIM );

      // plasti kappa
    if ( materi_plasti_kappa ) {
      tmp = array_inproduct( inc_epp, inc_epp, MDIM*MDIM );
      new_kappa = old_kappa + sqrt(0.5*tmp);     
      if ( swit ) pri( "new_kappa", new_kappa );
      pl_unknowns[kap_indx] = new_kappa;
    }
    // materi_plasti_kappa_shear (manual 4.25): the size of the SHEAR
    // plastic strain rate (deviatoric part), independent of the total
    // kappa dof.
    if ( materi_plasti_kappa_shear ) {
      double epp_dev2[MDIM*MDIM];
      array_move( inc_epp, epp_dev2, MDIM*MDIM );
      {
        double mean_epp2 = 0.;
        for ( int idim=0; idim<MDIM; idim++ ) mean_epp2 += epp_dev2[idim*MDIM+idim];
        mean_epp2 /= MDIM;
        for ( int idim=0; idim<MDIM; idim++ ) epp_dev2[idim*MDIM+idim] -= mean_epp2;
      }
      tmp = array_inproduct( epp_dev2, epp_dev2, MDIM*MDIM );
      new_kapsh = old_kapsh + sqrt(0.5*tmp);
      if ( swit ) pri( "new_kapsh", new_kapsh );
      pl_unknowns[kapsh_indx] = new_kapsh;
    }
    else
      new_kapsh = old_kapsh;

      // cap1 hardening (materi_plasti_cap1_history): see the incremental
      // branch above for the law; pc grows with the converged cap
      // plastic volume strain of this step.
    if ( materi_plasti_cap1_history ) {
      double cap1_plasti_data[8], ccotphi = 0., deps_p_cv = 0.,
        lambda_star = 0., kappa_star = 0., K_ref = 0., p_ref = 0.,
        m_cap1 = 0.;
      long int cap1_length = 0;
      new_cap1pc = old_cap1pc;
      if ( get_group_data( GROUP_MATERI_PLASTI_CAP1, gr, element,
          new_unknowns, cap1_plasti_data, cap1_length,
          GET_IF_EXISTS ) ) {
        lambda_star = cap1_plasti_data[3];
        kappa_star = cap1_plasti_data[4];
        K_ref = cap1_plasti_data[5];
        p_ref = cap1_plasti_data[6];
        m_cap1 = cap1_plasti_data[7];
        ccotphi = cap1_plasti_data[0] * cos(cap1_plasti_data[1]) /
          sin(cap1_plasti_data[1]);
        deps_p_cv = - ( inc_epp[0] + inc_epp[4] + inc_epp[8] );
        if ( deps_p_cv<0. ) deps_p_cv = 0.;
        if ( lambda_star>kappa_star && K_ref>0. && p_ref>0. ) {
          tmp = old_cap1pc + ccotphi;
          new_cap1pc = old_cap1pc + deps_p_cv * K_ref /
            ( lambda_star/kappa_star - 1. ) *
            scalar_power( tmp/p_ref, m_cap1 );
        }
      }
      if ( swit ) pri( "new_cap1pc", new_cap1pc );
    }

      // plasti rho
    if ( materi_plasti_rho ) {
      array_multiply( inc_epp, inc_rho, 
        plasti_kinematic_hardening, MDIM*MDIM );
      array_add( old_rho, inc_rho, new_rho, MDIM*MDIM );
      if ( swit ) pri( "new_rho", new_rho, MDIM, MDIM );
    }

      // new plastic strain
    array_add( inc_epp, old_epp, new_epp, MDIM*MDIM );
    if ( swit ) pri( "plastic strain", new_epp, MDIM, MDIM );

      // elastic strain part in elastic-plastic strain part
    array_subtract( inc_ept, inc_epp, inc_epe, MDIM*MDIM );
    array_subtract( new_ept, new_epp, new_epe, MDIM*MDIM );
    if ( condif_temperature ) {
      array_add( inc_epe, inc_temperature_strain, inc_epe, MDIM*MDIM );
      array_add( new_epe, new_temperature_strain, new_epe, MDIM*MDIM );
    }

      // membrane iterations
    membrane_found = membrane_iter = 0;
    while ( !membrane_found && membrane_iter<MAX_ITER ) {
      membrane_iter++;
      if ( swit ) {
        pri( "membrane_iter", membrane_iter );
        pri( "inc_epe", inc_epe, MDIM, MDIM );
        pri( "new_epe", new_epe, MDIM, MDIM );
      }

        // elasticity
      if ( formulation==TOTAL ) 
         matrix_a4b( C, new_epe, new_sig );
      else {
         matrix_a4b( C, inc_epe, work );
         array_add( work, rotated_old_sig, new_sig, MDIM*MDIM );
      }
      if ( swit ) pri( "stress after elasticity", new_sig, MDIM, MDIM );

        // membrane stiffness
      memmat[0][0] = C[0][0][0][0];
      memmat[0][1] = C[0][0][1][1];
      memmat[0][2] = C[0][0][2][2];
      memmat[1][0] = C[1][1][0][0];
      memmat[1][1] = C[1][1][1][1];
      memmat[1][2] = C[1][1][2][2];
      memmat[2][0] = C[2][2][0][0];
      memmat[2][1] = C[2][2][1][1];
      memmat[2][2] = C[2][2][2][2];

        // user supplied
      user_sigma( user_data, new_unknowns, inc_epe, 
        old_hisv, new_hisv, rotated_old_sig, new_sig, Cuser );
      stress_umat( element, gr, formulation, nuser_data, user_data, coord_ip,
        old_hisv, new_hisv, old_unknowns, new_unknowns,
        inc_ept, new_ept, rotated_old_sig, new_sig, 
        old_deften, new_deften, inc_rot, ddsdde );
      if ( swit ) pri( "stress after user supplied", new_sig, MDIM, MDIM );

        // hypoplasticity (ignored with plasti_apply -no: measured on the
        // Professional 02-08-2026 the hypoplastic stress stays frozen)
      if ( !plasti_off ) hypoplasticity( element, gr, formulation,
        old_hisv, new_hisv, old_unknowns, new_unknowns,
        inc_ept, old_epi, new_epi, rotated_old_sig, 
        new_sig, &Chypo[0][0][0][0], softvar_nonl, softvar_l );
      if ( swit ) pri( "stress after hypoplasticity", new_sig, MDIM, MDIM );

        // compressibility
      if ( compressibility!=0. ) {
        if      ( memory==-TOTAL_LINEAR || memory==-UPDATED_LINEAR ||
                  memory==-UPDATED_AREA || memory==-UPDATED_JAUMANN ) {
          // linear volumetric stress: the increment comes from inc_epe,
          // which is valid whether or not materi_strain_elasti is
          // initialized (old_epe is 0 without it, so new_epe-old_epe
          // would add the total strain every step)
          tmp_new = ( new_epe[0] + new_epe[4] + new_epe[8] ) / compressibility;
          tmp_inc = ( inc_epe[0] + inc_epe[4] + inc_epe[8] ) / compressibility;
        }
        else if ( memory==-TOTAL ) {
          array_move( old_epe, work, MDIM*MDIM );
          for ( idim=0; idim<MDIM; idim++ ) work[idim*MDIM+idim] += 1.;
          tmp_old = (matrix_determinant(work,MDIM)-1.) / compressibility;
          array_move( new_epe, work, MDIM*MDIM );
          for ( idim=0; idim<MDIM; idim++ ) work[idim*MDIM+idim] += 1.;
          tmp_new = (matrix_determinant(work,MDIM)-1.) / compressibility;
          tmp_inc = tmp_new - tmp_old;
        }
        else if ( memory==-TOTAL_PIOLA ) {
          array_move( old_epe, work, MDIM*MDIM );
          array_multiply( work, work, 2., MDIM*MDIM );
          for ( idim=0; idim<MDIM; idim++ ) work[idim*MDIM+idim] += 1.;
          tmp_old = (sqrt(matrix_determinant(work,MDIM))-1.) / compressibility;
          array_move( new_epe, work, MDIM*MDIM );
          array_multiply( work, work, 2., MDIM*MDIM );
          for ( idim=0; idim<MDIM; idim++ ) work[idim*MDIM+idim] += 1.;
          tmp_new = (sqrt(matrix_determinant(work,MDIM))-1.) / compressibility;
          tmp_inc = tmp_new - tmp_old;
        }
        else {
          tmp_inc = ( inc_epe[0] + inc_epe[4] + inc_epe[8] ) /  compressibility;
        }
        for ( idim=0; idim<MDIM; idim++ ) {
          if ( formulation==TOTAL ) new_sig[idim*MDIM+idim] += tmp_new;
          else new_sig[idim*MDIM+idim] += tmp_inc;
        }
        if ( swit ) pri( "stress after compressibility", new_sig, MDIM, MDIM );
      }

        // hyperelasticity
      hyperelasticity( gr, element, memory, old_unknowns, old_epe, 
        old_work, Chyper );
      hyperelasticity( gr, element, memory, new_unknowns, new_epe, 
        new_work, Chyper );
      if ( formulation==TOTAL )
        array_move( new_work, work, MDIM*MDIM );
      else
        array_subtract( new_work, old_work, work, MDIM*MDIM );
      array_add( work, new_sig, new_sig, MDIM*MDIM );
      memmat[0][0] += Chyper[0][0][0][0];
      memmat[0][1] += Chyper[0][0][1][1];
      memmat[0][2] += Chyper[0][0][2][2];
      memmat[1][0] += Chyper[1][1][0][0];
      memmat[1][1] += Chyper[1][1][1][1];
      memmat[1][2] += Chyper[1][1][2][2];
      memmat[2][0] += Chyper[2][2][0][0];
      memmat[2][1] += Chyper[2][2][1][1];
      memmat[2][2] += Chyper[2][2][2][2];
      if ( swit ) pri( "stress after hyperelasticity", new_sig, MDIM, MDIM );

        // viscosity
      viscous_stress( element, gr, user_data, old_unknowns,
        old_grad_old_unknowns, old_work, 
        viscosity, viscosity_heat_generation );
      viscous_stress( element, gr, user_data, new_unknowns,
        new_grad_new_unknowns, new_work, 
        viscosity, viscosity_heat_generation );
      if ( formulation==TOTAL )
        array_move( new_work, work, MDIM*MDIM );
      else
        array_subtract( new_work, old_work, work, MDIM*MDIM );
      array_add( work, new_sig, new_sig, MDIM*MDIM );
      if ( swit ) pri( "stress after viscosity", new_sig, MDIM, MDIM );

        // viscoelasticity
      visco_elasticity( element, gr, formulation, new_unknowns, inc_epe, 
        rotated_old_msig, new_sig, new_msig, memmat );
      if ( swit ) pri( "stress after visco elasticity", new_sig, MDIM, MDIM );

        // damage (manual Professional 2.2.5: sigma_damaged =
        // (1-d) sigma_undamaged). With materi_damage, materi() passes the
        // EFFECTIVE (undamaged) stress of the previous step in
        // rotated_old_sig (record element_intpnt_materi_stress_effective),
        // so new_sig is here the effective stress: keep it in new_sig_eff
        // (materi() stores it for the next step), evaluate the damage with
        // it and derive the damaged stress. Before 2026-09-29 the old
        // DAMAGED stress was carried and multiplied again by (1-d) every
        // step (sigma_n = (1-d_n)(sigma_n-1 + C deps)), so the damage
        // compounded over the steps (legacy defect, also in the 2014
        // source). The membrane stiffness is scaled like the stress. The
        // plastic yield test below keeps acting on the DAMAGED stress
        // (measured on the Professional 02-08-2026: von Mises yield
        // stress above the damaged peak stress gives no plastic strain).
        // materi_damage_apply / control_materi_damage_apply -no: the
        // damage data is ignored for the step (manual Professional
        // 6.799): d stays, no (1-d) scaling, and materi() passes the
        // damaged stress itself (measured on the Professional 02-08-2026:
        // the increments of a gated step are added undamaged).
      if ( materi_damage && materi_apply_gate_off( MATERI_DAMAGE_APPLY,
           CONTROL_MATERI_DAMAGE_APPLY ) )
        new_damage = old_damage;
      else if ( materi_damage ) {
        array_move( new_sig, new_sig_eff, MDIM*MDIM );
        damage( gr, new_epe, new_sig_eff, old_damage, new_damage );
        array_multiply( new_sig_eff, new_sig, 1.-new_damage, MDIM*MDIM );
        array_multiply( &memmat[0][0], &memmat[0][0], 1.-new_damage, MDIM*MDIM );
        damage_done = 1;
        if ( swit ) pri( "stress after damage", new_sig, MDIM, MDIM );
      }

          // membrane iterations ready?
      membrane_found = membrane_apply( element, gr, memmat, inc_ept,
        inc_epe, new_ept, new_epe, new_sig );
	//I have always membrane_found==TRUE
    }

      // test stresses for plastic yield functions
    // direct stress cut-off on a plane (group_materi_plasti_mohr_coul_
    // direct[_normal[_automatic]] / tension_direct[_normal[_automatic]]):
    // applied on the elastic stress BEFORE the plastic-yield test.
    // group_materi_plasti_pressure_limit / _coord_limit (manual
    // Professional 6.688/6.689): neglect the direct plasticity laws when
    // the pressure exceeds pressure_limit (free-surface problems; the
    // pressure is taken positive in compression) or when the vertical
    // coordinate exceeds coord_limit.
    {
      double pressure_now = -(
        new_sig[0] + new_sig[4] + new_sig[8] ) / 3.;
      double limit_gate = 0.;
      long int gate_off = plasti_off;
      if ( db( GROUP_MATERI_PLASTI_PRESSURE_LIMIT, gr, idum,
          &limit_gate, ldum, VERSION_NORMAL, GET_IF_EXISTS ) )
        if ( pressure_now>limit_gate ) gate_off = 1;
      if ( !gate_off &&
           db( GROUP_MATERI_PLASTI_COORD_LIMIT, gr, idum,
          &limit_gate, ldum, VERSION_NORMAL, GET_IF_EXISTS ) )
        if ( coord_ip[ndim-1]>limit_gate ) gate_off = 1;
      // K13: keep the stress before the direct cut-offs (their plastic
      // strain is C^-1 (sig_before - sig_after), stored after the loop)
      // and the yield value of the direct laws on the final stress
      direct_on = 0;
      direct_f = NO_YIELD_F;
      array_move( new_sig, direct_sig_before, MDIM*MDIM );
      if ( !gate_off ) {
        if ( db_active_index( GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT, gr,
            VERSION_NORMAL ) ||
             db_active_index( GROUP_MATERI_PLASTI_TENSION_DIRECT, gr,
            VERSION_NORMAL ) ) {
          // Without a plane normal (_normal / _normal_automatic) the
          // direct records are the FULL Mohr-Coulomb / tension spectral
          // cut-offs on the principal stresses (manual Professional
          // 6.726/6.738); when either law has a plane normal, each law
          // with a normal limits the traction on its own plane
          // (6.727/6.739, handled by materi_direct_cutoff).
          double mc_test[MDIM], ten_test[MDIM];
          array_move( mc_normal, mc_test, MDIM );
          array_move( ten_normal, ten_test, MDIM );
          if ( !array_normalize( mc_test, MDIM ) &&
               !array_normalize( ten_test, MDIM ) )
            materi_direct_full_mc( element, gr, plasti_on_boundary, dtime,
              new_sig, ddsdde, 1, direct_f );
          else {
            // elastic stiffness of the element for the phi_flow return
            // (the same linear stiffness that fills ddsdde below)
            double Cel[MDIM][MDIM][MDIM][MDIM];
            for ( idim=0; idim<MDIM; idim++ )
              for ( jdim=0; jdim<MDIM; jdim++ )
                for ( kdim=0; kdim<MDIM; kdim++ )
                  for ( ldim=0; ldim<MDIM; ldim++ ) {
                    Cel[idim][jdim][kdim][ldim] =
                      Cmem[idim][jdim][kdim][ldim] +
                      Cuser[idim][jdim][kdim][ldim] +
                      Chyper[idim][jdim][kdim][ldim] +
                      Chypo[idim][jdim][kdim][ldim];
                    if ( compressibility!=0. && idim==jdim && kdim==ldim )
                      Cel[idim][jdim][kdim][ldim] += 1./compressibility;
                  }
            materi_direct_cutoff( element, gr, plasti_on_boundary, dtime,
              new_sig, ddsdde, mc_normal, ten_normal, &Cel[0][0][0][0],
              direct_f );
          }
        }
        materi_compression_cutoff( element, gr, dtime, new_sig, direct_f );
        direct_on = ( direct_f!=NO_YIELD_F );
      }
    }
    array_move( new_sig, test_sig, MDIM*MDIM );
    if ( materi_plasti_rho ) 
      array_subtract( test_sig, new_rho, test_sig, MDIM*MDIM );
      // plasti_apply / control_plasti_apply -no: no yield test, the
      // elastic trial state is the solution
    if ( plasti_off ) {
      plasti_found = 1;
      continue;
    }
      // new lambda
    if ( plasti_iter==1 ) {
      plasti_type = -NONE; plasti_rule( element, gr, plasti_on_boundary, user_data, 
        pl_unknowns, new_grad_new_unknowns, old_hisv, new_hisv,
        old_epp, inc_epp, inc_ept, GET_YIELD_RULE, 
        plasti_type, test_sig, f, new_f, ddum );
      if ( visco_named ) {
        viscoplasti = visco_named_select( plasti_type, n_visco_names,
          visco_names, visco_values, visco_pair );
        if ( viscoplasti && visco_named==2 ) {
          gamma = plasti_visco_exponential[0] = visco_pair[0];
          alpha = plasti_visco_exponential[1] = visco_pair[1];
        }
        else if ( viscoplasti ) {
          eta = plasti_visco_power[0] = visco_pair[0];
          p = plasti_visco_power[1] = visco_pair[1];
          f_ref = 1.;
        }
        max_plasti_iter = viscoplasti ? max_iter_visco : max_iter_ep;
      }

      if ( materi_plasti_f_nonlocal && viscoplasti ) f = new_unknowns[fn_indx];
      if ( swit ) {
        pri( "plasti_iter", plasti_iter );
        pri( "plasti_type", -plasti_type );
        pri( "f", f );
      }
      //if f is negative, solution is elastic (first iteration->lambda=0)
      //otherwise plastic iterations start
      if ( ( f<EPS_F && viscoplasti_always==-NO ) || f==NO_YIELD_F ) 
        plasti_found = 1;
      else {
        plastic_ip = 1;
        if ( viscoplasti ) {
          lambda_previous = 0.; f_previous = f;
          if ( plasti_visco_exponential[0]>0. ) {
            tmp = alpha*f;
            if ( tmp>visco_limit ) tmp = visco_limit;
            lambda = gamma * (-pressure) * exp(tmp);
          }
          else {
            assert( plasti_visco_power[0]>0. );
            lambda = eta * scalar_power(f/f_ref,p);
          }
        }
        else {
          lambda_previous = 0.; f_previous = f;
          lambda = EPS_LAMBDA;
        }
        plasti_rule( element, gr, plasti_on_boundary, user_data, 
          pl_unknowns, new_grad_new_unknowns,
          old_hisv, new_hisv, old_epp, inc_epp, inc_ept,
          GET_FLOW_RULE_GRAD, plasti_type, test_sig, rdum, rdum, plasti_dir );
        if ( !viscoplasti ) array_normalize( plasti_dir, MDIM*MDIM );
        if ( swit ) pri( "plasti_dir", plasti_dir, MDIM*MDIM );
      }
    }
    else {
      plasti_rule( element, gr, plasti_on_boundary, user_data, pl_unknowns, new_grad_new_unknowns,
        old_hisv, new_hisv, old_epp, inc_epp, inc_ept, GET_YIELD_RULE, plasti_type, 
        test_sig, f, new_f, ddum );
      if ( materi_plasti_f_nonlocal && viscoplasti ) f = new_unknowns[fn_indx];
      if ( f>DBL_MAX/1.e6 ) {
        pri( "Error detected in plasticity." );
        pri( "Maybe too large time steps. Try smaller time steps." );
        exit(TN_EXIT_STATUS);
      }
      if ( swit ) {
        pri( "plasti_iter", plasti_iter );
        pri( "plasti_type", plasti_type );
        pri( "lambda", lambda );
        pri( "f", f );
      }
      if ( scalar_dabs(f)<EPS_F || scalar_dabs(f-f_previous)<EPS_EPS_F*EPS_F ) {
        array_add( work_inc_epp, total_inc_epp, total_inc_epp, MDIM*MDIM );
        plasti_type = -NONE; plasti_rule( element, gr, plasti_on_boundary, user_data, 
          pl_unknowns, new_grad_new_unknowns, old_hisv, new_hisv,
          old_epp, inc_epp, inc_ept, GET_YIELD_RULE, 
          plasti_type, test_sig, f, new_f, ddum );
        if ( f<EPS_F ) 
          plasti_found = 1;
        else {
          plasti_iter = 0; lambda = 0.;
        }
      }
      else if ( viscoplasti ) {
        lambda_previous = lambda; f_previous = f;
        if ( plasti_visco_exponential[0]>0. ) {
          tmp = alpha*f;
          if ( tmp>visco_limit ) tmp = visco_limit;
          lambda = gamma * (-pressure) * exp(tmp);
        }
        else {
          assert( plasti_visco_power[0]>0. );
          lambda = eta * scalar_power(f/f_ref,p);
        }
      }
      else {
        if      ( f>f_previous || f<0. ) {
            // adjust steps if f becomes larger
          lambda = ( 0.5 * lambda + 0.5 * lambda_previous );
        }
        else {
          tmp = - (lambda-lambda_previous)*f_previous/(f-f_previous);
          if ( deps_size>EPS_SIZE ) {
            if ( tmp>EPS_TMP ) tmp = EPS_TMP;
            if ( tmp<-EPS_TMP ) tmp = -EPS_TMP;
          }
          lambda_new = lambda_previous + tmp;
          lambda_previous = lambda; f_previous = f;
          lambda = lambda_new;
        }
      }
    }
  }
  }

  // K13: plastic strain and yield value of the direct cut-offs. The
  // direct laws cut the stress without an incremental plastic strain;
  // the cut is stored as the plastic strain that the elastic law needs
  // for sig = C (eps - eps_p): deps_p = C^-1 (sig_before - sig_after),
  // with C the elastic stiffness (Cmem + Cuser + Chyper + Chypo +
  // compressibility, as the phi_flow return). The Professional 02-08-2026
  // stores the same normal components (tension cut: eppxx 1.4e-3,
  // eppyy = eppzz -0.35e-3 for E 1, nu 0.25, trial 2.4e-3 -> 1e-3) but
  // twice the tensor value in the shear components (Delta sig_xy / G);
  // the GNU keeps the tensor convention of every other plastic strain
  // (deliberate divergence, BACKLOG). The yield value is the largest
  // direct-law value on the final stress, 0 when a law cut (also a
  // partial visco cut, as the Professional).
  if ( direct_on ) {
    double dsig[MDIM*MDIM], deps[MDIM*MDIM];
    array_subtract( direct_sig_before, new_sig, dsig, MDIM*MDIM );
    if ( array_size( dsig, MDIM*MDIM )>0. &&
         direct_plastic_strain( dsig, Cmem, Cuser, Chyper, Chypo,
           compressibility, deps ) ) {
      plastic_ip = 1;
      array_add( inc_epp, deps, inc_epp, MDIM*MDIM );
      array_subtract( inc_epe, deps, inc_epe, MDIM*MDIM );
      array_subtract( new_epe, deps, new_epe, MDIM*MDIM );
    }
    if ( direct_f>0. ) direct_f = 0.;
    if ( new_f==NO_YIELD_F || direct_f>new_f ) new_f = direct_f;
  }

  // group_materi_plasti_residual_stiffness (manual Professional 6.735):
  // fraction of the original shear stiffness used in the material matrix
  // of a plastic integration point {a yield test that started plastic
  // iterations, or a direct cut, in this stress update}. Only the
  // matrix changes, the stresses do not. The volumetric part
  // delta_ij delta_kl S/9, S = C_mmnn, is kept; the rest of C is scaled.
  // The Professional 02-08-2026 accepts the record without a measurable
  // effect {DIVERGENCES PLASTI-RESIDUAL-STIFFNESS}.
  double residual_factor = 1., c_vol = 0.;
  if ( get_group_data( GROUP_MATERI_PLASTI_RESIDUAL_STIFFNESS, gr, element,
      new_unknowns, &residual_factor, ldum, GET_IF_EXISTS ) ) {
    if ( residual_factor<=0. ) {
      pri( "Error: group_materi_plasti_residual_stiffness must be "
        "larger than 0." );
      exit(1);
    }
  }
  if ( !plastic_ip ) residual_factor = 1.;
  if ( residual_factor!=1. ) {
    for ( idim=0; idim<MDIM; idim++ ) {
      for ( kdim=0; kdim<MDIM; kdim++ ) {
        c_vol += Cmem[idim][idim][kdim][kdim] + Cuser[idim][idim][kdim][kdim]
          + Chyper[idim][idim][kdim][kdim] + Chypo[idim][idim][kdim][kdim];
      }
    }
    c_vol /= 9.;
  }

    // fill material stiffness, we use the linear stiffness
  for ( idim=0; idim<MDIM; idim++ ) {
    for ( jdim=idim; jdim<MDIM; jdim++ ) {
      for ( kdim=0; kdim<MDIM; kdim++ ) {
        for ( ldim=kdim; ldim<MDIM; ldim++ ) {
          if ( kdim==ldim ) 
            fac = 1.0;
          else 
            fac = 0.5;
          i = stress_indx(idim,jdim);
          j = stress_indx(kdim,ldim);
          ind_ddsdde = i*MSTRAIN+j;
          double c_ijkl = Cmem[idim][jdim][kdim][ldim]
            + Cuser[idim][jdim][kdim][ldim] + Chyper[idim][jdim][kdim][ldim]
            + Chypo[idim][jdim][kdim][ldim];
          if ( residual_factor!=1. ) {
            double vol = ( idim==jdim && kdim==ldim ) ? c_vol : 0.;
            c_ijkl = vol + residual_factor * ( c_ijkl - vol );
          }
          ddsdde[ind_ddsdde] += fac * c_ijkl;
          if ( compressibility!=0. && idim==jdim && kdim==ldim )
            ddsdde[ind_ddsdde] += 1./compressibility;
        }
      }
    }
  }

  // effective stress without damage (no damage model, damage switched
  // off, or the incremental plasticity branch that never applies damage):
  // the stress itself
  if ( !damage_done ) array_move( new_sig, new_sig_eff, MDIM*MDIM );

  // group_materi_factor (manual Professional 6.666): multiplication
  // factor for the material stresses AND stiffness (unit conversion of
  // a stress law specified in other units than the calculation).
  {
    double materi_factor = 1.;
    if ( get_group_data( GROUP_MATERI_FACTOR, gr, element,
        new_unknowns, &materi_factor, ldum, GET_IF_EXISTS ) ) {
      if ( scalar_dabs(materi_factor-1.)>TINY ) {
        array_multiply( new_sig, new_sig, materi_factor, MDIM*MDIM );
        array_multiply( new_sig_eff, new_sig_eff, materi_factor, MDIM*MDIM );
        array_multiply( ddsdde, ddsdde, materi_factor, MSTRAIN*MSTRAIN );
      }
    }
  }
  if ( swit ) {
    pri( "plasti_found", plasti_found );
    pri( "membrane_found", membrane_found );
    pri( "stress", new_sig, MDIM, MDIM );
    pri( "ddsdde", ddsdde, MSTRAIN, MSTRAIN );
    pri( "Out routine SET_STRESS" );
  }
}

void stress_umat( long int element, long int gr, long int formulation,
  long int nuser_data, double user_data[], double coord_ip[],
  double old_hisv[], double new_hisv[], 
  double old_unknowns[], double new_unknowns[], 
  double inc_ept[], double new_ept[], 
  double rotated_old_sig[], double new_sig[], 
  double old_deften[], double new_deften[], 
  double inc_rot[], double ddsdde[] )

  /* Interface routine to umat routine. */

{
  double stress[MSTRAIN], statev[DATA_ITEM_SIZE], 
   sse[1], spd[1], scd[1], rpl[1], ddsddt[1], drplde[1], drpldt[1], 
   stran[MSTRAIN], dstran[MSTRAIN], time[1], dtime[1], temp[1], dtemp[1], 
   predef[1], dpred[1];
  char cmname[1];
  long int ndi[1], nshr[1], ntens[1], nstatv[1];
  double props[DATA_ITEM_SIZE];
  long int nprops[1];
  double coords[MDIM], drot[MDIM*MDIM], pnewdt[1], celent[1], 
    dfgrd0[MDIM*MDIM], dfgrd1[MDIM*MDIM];
  long int noel[1], npt[1], layer[1], kspt[1], kstep[1], kinc[1];
  short cmname_len;

  long int icontrol=0, number_of_iterations=0, ldum=0, user_umat=-NO, idum[1];
  long int istrain=0, jstrain=0, idim=0, jdim=0, kdim=0, ldim=0, i=0, j=0, 
    ind_ddsdde=0, abaqus_indx[MSTRAIN][2];
  double ddum[1], ddsdde_tmp[MSTRAIN*MSTRAIN];

  db( NUMBER_ITERATIONS, 0, &number_of_iterations, ddum, ldum, VERSION_NEW, GET );
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  db( GROUP_USER_UMAT, gr, &user_umat, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  if ( user_umat==-YES ) {

      // dummies

    pnewdt[0] = 0.;
    celent[0] = 0.;
    sse[0] = 0.;
    spd[0] = 0.;
    scd[0] = 0.;
    rpl[0] = 0.;
    ddsddt[0] = 0.;
    drplde[0]= 0.;
    drpldt[0] = 0.;
    predef[0] = 0.;
    dpred[0] = 0.;
    npt[0] = 0;
    layer[0] = 0;
    kspt[0] = 0;
    kstep[0] = icontrol;
    kinc[0] = number_of_iterations;
    cmname_len = 0;

      // real data

    array_move( coord_ip, coords, ndim );

    noel[0] = element;

    db( DTIME, 0, idum, dtime, ldum, VERSION_NEW, GET );
    db( TIME_CURRENT, 0, idum, time, ldum, VERSION_NORMAL, GET );

    if ( condif_temperature ) {
      temp[0] = old_unknowns[temp_indx];
      dtemp[0] = new_unknowns[temp_indx] - old_unknowns[temp_indx];
    }
    else {
      temp[0] = 0.;
      dtemp[0] = 0.;
    }

    nprops[0] = nuser_data;
    nstatv[0] = materi_history_variables;
    ndi[0] = 3;
    nshr[0] = 3;
    ntens[0] = 6;

    array_set( ddsdde_tmp, 0., MSTRAIN*MSTRAIN );

    array_move( user_data, props, nprops[0] );

    array_move( old_hisv, statev, nstatv[0] );

    stress[0] = rotated_old_sig[0];
    stress[1] = rotated_old_sig[4];
    stress[2] = rotated_old_sig[8];
    stress[3] = rotated_old_sig[1];
    stress[4] = rotated_old_sig[2];
    stress[5] = rotated_old_sig[5];

    stran[0] = new_ept[0] - inc_ept[0];
    stran[1] = new_ept[4] - inc_ept[4];
    stran[2] = new_ept[8] - inc_ept[8];
    stran[3] = 2. * ( new_ept[1] - inc_ept[1] );
    stran[4] = 2. * ( new_ept[2] - inc_ept[2] );
    stran[5] = 2. * ( new_ept[5] - inc_ept[5] );

    dstran[0] = inc_ept[0];
    dstran[1] = inc_ept[4];
    dstran[2] = inc_ept[8];
    dstran[3] = 2. * inc_ept[1];
    dstran[4] = 2. * inc_ept[2];
    dstran[5] = 2. * inc_ept[5];

    for ( i=0; i<MDIM*MDIM; i++ ) {
      drot[i] = inc_rot[i];
      dfgrd0[i] = old_deften[i];
      dfgrd1[i] = new_deften[i];
    }    

      // stress contribution by umat
    umat_( stress, statev, ddsdde_tmp, sse, spd, scd, rpl, ddsddt,
      drplde, drpldt, stran, dstran, time, dtime, temp, dtemp, predef, 
      dpred, cmname, ndi, nshr, ntens, nstatv, props, nprops, coords, drot, 
      pnewdt, celent, dfgrd0, dfgrd1, noel, npt, layer, 
      kspt, kstep, kinc, cmname_len );

      // from abaqus fortran to tochnog c++
    abaqus_indx[0][0] = 0;
    abaqus_indx[0][1] = 0;
    abaqus_indx[1][0] = 1;
    abaqus_indx[1][1] = 1;
    abaqus_indx[2][0] = 2;
    abaqus_indx[2][1] = 2;
    abaqus_indx[3][0] = 0;
    abaqus_indx[3][1] = 1;
    abaqus_indx[4][0] = 0;
    abaqus_indx[4][1] = 2;
    abaqus_indx[5][0] = 1;
    abaqus_indx[5][1] = 2;
    for ( istrain=0; istrain<MSTRAIN; istrain++ ) {
      idim = abaqus_indx[istrain][0];
      jdim = abaqus_indx[istrain][1];
      i = stress_indx(idim,jdim); 
      for ( jstrain=0; jstrain<MSTRAIN; jstrain++ ) {
        kdim = abaqus_indx[jstrain][0];
        ldim = abaqus_indx[jstrain][1];
        j = stress_indx(kdim,ldim); 
        ind_ddsdde = i*MSTRAIN+j;
        ddsdde[ind_ddsdde] += ddsdde_tmp[jstrain*MSTRAIN+istrain];
      }
    }

    array_move( statev, new_hisv, nstatv[0] );

    if ( formulation==TOTAL ) {
      new_sig[0] += stress[0];
      new_sig[4] += stress[1];
      new_sig[8] += stress[2];
      new_sig[1] += stress[3];
      new_sig[3] += stress[3];
      new_sig[2] += stress[4];
      new_sig[6] += stress[4];
      new_sig[5] += stress[5];
      new_sig[7] += stress[5];
    }
    else {
      assert( formulation==INCREMENTAL );
      new_sig[0] += stress[0] - rotated_old_sig[0];
      new_sig[4] += stress[1] - rotated_old_sig[4];
      new_sig[8] += stress[2] - rotated_old_sig[8];
      new_sig[1] += stress[3] - rotated_old_sig[1];
      new_sig[3] += stress[3] - rotated_old_sig[3];
      new_sig[2] += stress[4] - rotated_old_sig[2];
      new_sig[6] += stress[4] - rotated_old_sig[6];
      new_sig[5] += stress[5] - rotated_old_sig[5];
      new_sig[7] += stress[5] - rotated_old_sig[7];
    }

  }
  
}
