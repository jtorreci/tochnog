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

#define EPS_SIG 1.e-5
#define EPS_DIR 1.e-10
#define EPS_PRISCO_1 1.e-8
#define EPS_PRISCO_2 1.e-5

void plasti_rule( long int element, long int gr, 
  long int plasti_on_boundary, double user_data[],
  double new_unknowns[], double new_grad_new_unknowns[],
  double old_hisv[], double new_hisv[], double old_epp[],
  double inc_epp[], double inc_ept[], long int task, long int &plasti_type, 
  double sig[], double &f, double &new_f, double dir[] )

{

  long int i=0, swit=0, idim=0, jdim=0, ind=0, tmp_plasti_type=-1, in_apex=0, 
    length=0, group_materi_plasti_user=-NO, 
    test1=0, test2=0, test3=0, test4=0, test5=0,
    test6=0, test7=0, test8=0, test9=0, ldum=0, icontrol=0,
    tension_cutoff=-NO, tension_limit=-NO,
    options_skip_plasticity=-NO, idum[1];
  double sig_yield=0., sigm=0., signorm=0., K=0., K_flow=0., alpha=0.,
    dev=0., devp=0., q=0., p=0., C=0., kappa_0=0., n=0., options_nonlocal=0.,
    f_flow_right=0., f_flow_left=0., phi=0., c=0., sig0=0.,
    phi_flow=0., alpha_flow=0., f_yield=NO_YIELD_F, f_flow=NO_YIELD_F, 
    sig_eq=0., tmp_sig_eq=0., void_fraction=0., tmp=0., tmp1=0., tmp2=0.,
    q1=0., q2=0., q3=0., m1=0., m20=0., m21=0., m22=0., m23=0.,
    prisco_gamma=0., prisco_j2eta=0., prisco_j3eta=0.,
    prisco_r=0., prisco_betaf=0., prisco_rc=0., prisco_cp=0., 
    prisco_kinc=0., prisco_bp=0., prisco_ksi=0., prisco_ksi_c=0.,
    prisco_ksi_e=0., prisco_thetahat=1.,
    prisco_theta=0., prisco_thetahat_c=0., prisco_thetahat_e=0.,
    prisco_tp=0., prisco_betafhat=0., 
    prisco_r0=0., prisco_st0=0., prisco_b4=0., prisco_b5=0., prisco_b7=0.,
    prisco_beta=0., prisco_betaf0=0., prisco_r1tmp=0.,
    prisco_r2tmp=0., prisco_rt=0.,
    pa=0., pb=0., R=0., t=0., epp_vol=0., 
    phi0=0., c0=0., phi0_flow=0., phi1=0., c1=0., 
    phi1_flow=0., kappa=0., kappa_crit=0., I1=0., I2=0., I3=0.,
    plasti_on_boundary_factor=BOUNDARY_REDUCTION_FACTOR,
    m=0., lambda=0., e=0., N=0., de=0., p0=0., dp0=0., rdum=0.,
    M=0., pc=0., p_star=0., pc_star=0.,
    hs_elasti[7], hs_sig3=0., hs_base=0., hs_E50=0., hs_Eur=0., hs_qf=0.,
    hs_qa=0., hs_gp_extra=0., hs_ccot=0., hs_Rf=0., gammap=0.,
    prisco_rv[MDIM][MDIM], prisco_st[MDIM][MDIM], prisco_rr[MDIM][MDIM],
    prisco_st1[MDIM][MDIM], prisco_st2[MDIM][MDIM], prisco_sv0[MDIM][MDIM],
    prisco_chi[MDIM*MDIM], prisco_chihat[MDIM*MDIM],
    prisco_ss[MDIM*MDIM], prisco_etas[MDIM*MDIM], prisco_epinc[MDIM*MDIM],
    prisco_epp_dev[MDIM*MDIM], prisco_incepp_dev[MDIM*MDIM],
    sig_princ_apex[MDIM], sig_dev[MDIM*MDIM], sig_tmp[MDIM*MDIM], 
    plasti_data[DATA_ITEM_SIZE], sig_princ[MDIM], ddum[MDIM*MDIM],
    e_vec[MDIM], f_vec[MDIM], work[MDIM*MDIM],
    geo_sig[MDIM*MDIM], geo_inc_epp[MDIM*MDIM], 
    geo_old_epp[MDIM*MDIM], geo_inc_ept[MDIM*MDIM];

  db( OPTIONS_SKIP_PLASTICITY, 0, &options_skip_plasticity, ddum, ldum, 
    VERSION_NORMAL, GET_IF_EXISTS );
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_OPTIONS_SKIP_PLASTICITY, icontrol, &options_skip_plasticity, ddum, ldum, 
    VERSION_NORMAL, GET_IF_EXISTS );
  if ( options_skip_plasticity==-YES ) return;

  swit = set_swit(element,-1,"plasti_rule");
  if ( task==GET_FLOW_RULE ) swit = 0;
  if ( swit ) pri( "In PLASTI_RULE" );

  if ( swit ) {
    if      ( task==GET_YIELD_RULE ) 
      pri( "GET_YIELD_RULE" );
    else if ( task==GET_FLOW_RULE ) {
      pri( "GET_FLOW_RULE" );
      pri( "plasti_type", plasti_type );
    }
    else if ( task==GET_FLOW_RULE_GRAD ) {
      pri( "GET_FLOW_RULE_GRAD" );
      pri( "plasti_type", plasti_type );
    }
    pri( "sig", sig, MDIM, MDIM );
  }

  db( OPTIONS_NONLOCAL, 0, idum, &options_nonlocal, ldum, 
    VERSION_NORMAL, GET_IF_EXISTS );
  db( GROUP_MATERI_PLASTI_USER, 0, &group_materi_plasti_user, ddum, ldum, 
    VERSION_NORMAL, GET_IF_EXISTS );
  db( GROUP_MATERI_PLASTI_BOUNDARY_FACTOR, gr, idum, &plasti_on_boundary_factor, 
    ldum, VERSION_NORMAL, GET_IF_EXISTS );

    // mean stress
  sigm = ( sig[0] + sig[4] + sig[8] ) / 3.;
  if ( swit ) pri( "sigm", sigm );

    // deviatoric strains and stresses
  array_move( sig, sig_dev, MDIM*MDIM );
  for ( idim=0; idim<MDIM; idim++ ) sig_dev[idim*MDIM+idim] -= sigm;

    // equivalent stress
  sig_eq = sqrt(scalar_dabs(array_inproduct(sig_dev,sig_dev,MDIM*MDIM))/2.);
  if ( swit ) pri( "sig_eq", sig_eq );

    // geoetechnical convention
  array_multiply( sig, geo_sig, -1., MDIM*MDIM );
  array_multiply( inc_ept, geo_inc_ept, -1., MDIM*MDIM );
  array_multiply( inc_epp, geo_inc_epp, -1., MDIM*MDIM );
  if ( materi_strain_plasti ) array_multiply( old_epp, geo_old_epp, -1., MDIM*MDIM );

  f = NO_YIELD_F;
  if ( get_group_data( GROUP_MATERI_PLASTI_CAMCLAY, gr, element, new_unknowns,
      plasti_data, ldum, GET_IF_EXISTS ) ) {
    test1 = task==GET_YIELD_RULE&&plasti_type==-NONE;
    test2 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_CAMCLAY;
    test3 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_CAMCLAY;
    if ( test1 || test2 || test3 ) {
      if ( swit ) pri( "plasti_camclay" );
      if ( materi_history_variables<2 ) {
        pri( "Error: GROUP_MATERI_PLASTI_CAMCLAY needs 'materi_history_variables 2'" );
        exit(1);
      }
      m = plasti_data[0];
      kappa = plasti_data[1];
      lambda = plasti_data[2];
      if ( scalar_dabs(lambda-kappa)==0. ) db_error( GROUP_MATERI_PLASTI_CAMCLAY, gr );
      e = old_hisv[0];
      // p0 is the stored preconsolidation pressure (cchis1). The GNU
      // legacy closure derived p0 from the void ratio e and the current
      // pressure p through the 4th material parameter N
      // (p0 = exp((N-kappa*ln(p)-(1+e))/(lambda-kappa))); the
      // Professional layout (manual 2.2.x, group_materi_plasti_camclay
      // 6.689 = M kappa lambda) tracks p0 EXPLICITLY as a history
      // variable (materi_plasti_camclay_history / node_dof -cchis1) and
      // derives N from the initial state. Fall back to the N-based
      // closure only for legacy inputs that still give 4 values AND have
      // no initialized p0.
      p = -sigm;
      p0 = old_hisv[1];
      if ( p0<=0. ) {
        if ( ldum>=4 ) {
          N = plasti_data[3];
          p0 = exp((N-kappa*log(p)-(1+e))/(lambda-kappa));
        }
        else {
          pri( "\nError: GROUP_MATERI_PLASTI_CAMCLAY needs initial p0."
               "\nSpecify initia 'materi_plasti_camclay_history' and initialize"
               " -cchis1 (p0) with control_reset_dof." );
          exit(1);
        }
      }
      if ( swit ) {
        pri( "e", e );
        pri( "p0", p0 );
      }
      dev = (geo_inc_ept[0]+geo_inc_ept[4]+geo_inc_ept[8]);
      de = -dev*(1.+e);
      devp = (geo_inc_epp[0]+geo_inc_epp[4]+geo_inc_epp[8]);
      dp0 = devp * p0 * (1.+e)/(lambda-kappa);
      if ( swit ) {
        pri( "dev", dev );
        pri( "devp", devp );
        pri( "de", de );
        pri( "dp0", dp0 );
      }
      e += de;
      p0 += dp0;
      q = sqrt( 
        0.5*( scalar_square(geo_sig[0]-geo_sig[4])+
              scalar_square(geo_sig[4]-geo_sig[8])+
              scalar_square(geo_sig[0]-geo_sig[8]) ) + 
        3.*( scalar_square(geo_sig[1]) +
             scalar_square(geo_sig[2]) +
             scalar_square(geo_sig[5]) ) );
      f_yield = q*q - m*m*(p*(p0-p));
      f_flow = f_yield;
      if      ( task==GET_YIELD_RULE ) {
        if (  f_yield>f ) {
          f = f_yield;
          tmp_plasti_type = GROUP_MATERI_PLASTI_CAMCLAY;
        }
        if ( swit ) pri( "f_yield", f_yield );

        new_hisv[0] = e;
        new_hisv[1] = p0;
      }
      else {
        assert( task==GET_FLOW_RULE );
        f = f_flow;
        if ( swit ) pri( "f_flow", f_flow );
      }
    }
  }
  if ( get_group_data( GROUP_MATERI_PLASTI_CAP, gr, element, new_unknowns, 
    plasti_data, length, GET_IF_EXISTS ) ) {
    test1 = task==GET_YIELD_RULE&&plasti_type==-NONE;
    test2 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_CAP;
    test3 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_CAP;
    if ( test1 || test2 || test3 ) {
      if ( swit ) pri( "check plasti_cap" );
        /* get data */
      c = plasti_data[0];
      phi = plasti_data[1];
      alpha = plasti_data[2];
      R = plasti_data[3];
      epp_vol = scalar_dabs( new_unknowns[epp_indx+0] + new_unknowns[epp_indx+3] + 
        new_unknowns[epp_indx+5] );
        /* epp_vol-pb table */
      length -= 4;
      if ( !table_xy( &plasti_data[4], "GROUP_MATERI_PLASTI_CAP",
          length, epp_vol, pb ) ) db_error( GROUP_MATERI_PLASTI_CAP, gr );
        /* invariants */
      p = -sigm;
      t = sqrt(3.) * sig_eq; 
        /* pa */
      tmp = (1.+R*tan(phi));
      if ( tmp==0. ) db_error( GROUP_MATERI_PLASTI_CAP, gr );
      pa = (pb-R*c)/tmp;
        /* yield rule and flow rule */
      if ( cos(phi)==0. ) db_error( GROUP_MATERI_PLASTI_CAP, gr );
      tmp = 1. + alpha - (alpha/cos(phi));
      if ( tmp==0. ) db_error( GROUP_MATERI_PLASTI_CAP, gr );
      f_yield = sqrt( scalar_square(p-pa) + scalar_square(R*t/tmp) ) - 
        R*(c+pa*tan(phi));
      f_flow = f_yield;
      if      ( task==GET_YIELD_RULE ) {
        if (  f_yield>f ) {
          f = f_yield;
          tmp_plasti_type = GROUP_MATERI_PLASTI_CAP;
        }
        if ( swit ) pri( "f_yield", f_yield );
      }
      else {
        assert( task==GET_FLOW_RULE );
        f = f_flow;
        if ( swit ) pri( "f_flow", f_flow );
      }
    }
  }
  if ( get_group_data( GROUP_MATERI_PLASTI_CAP1, gr, element, new_unknowns, 
    plasti_data, length, GET_IF_EXISTS ) ) {
    test1 = task==GET_YIELD_RULE&&plasti_type==-NONE;
    test2 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_CAP1;
    test3 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_CAP1;
    if ( test1 || test2 || test3 ) {
      if ( swit ) pri( "check plasti_cap1" );
        // manual Professional 6.691 (theory cap1): 8 parameters
        // phi c M lambda_star kappa_star K_ref p_ref m. M is READ from
        // the record (typically 6*sin(phi)/(3-sin(phi))). Surface:
        // f = q^2/M^2 + p*(p* - p*c) with p* = p + c*cot(phi),
        // p*c = pc + c*cot(phi), pc the materi_plasti_cap1_history
        // node dof (initial value given via node_dof). The hardening
        // (stress.cc) grows pc with the cap plastic volume strain, so
        // the surface translates along the p axis. Combined with shear
        // laws (DP, MC, ...) the standard driver picks the largest f.
      c = plasti_data[0];
      phi = plasti_data[1];
      M = plasti_data[2];
        // cot(phi) diverges at phi=0: the manual expects phi>0; the
        // tension pole at phi=pi/2 (cot = 0) is also rejected
      if ( cos(phi)==0. || sin(phi)==0. ) db_error( GROUP_MATERI_PLASTI_CAP1, gr );
      if ( M==0. ) db_error( GROUP_MATERI_PLASTI_CAP1, gr );
      if ( !materi_plasti_cap1_history ) {
        pri( "Error: GROUP_MATERI_PLASTI_CAP1 needs initia 'materi_plasti_cap1_history'" );
        exit(TN_EXIT_STATUS);
      }
      pc = new_unknowns[cap1_indx];
        // invariants: p compression-positive; q invariant to the sign
      p = -( sig[0] + sig[4] + sig[8] ) / 3.;
      q = sqrt( 
        0.5*( scalar_square(sig[0]-sig[4])+
              scalar_square(sig[4]-sig[8])+
              scalar_square(sig[0]-sig[8]) ) + 
        3.*( scalar_square(sig[1]) +
             scalar_square(sig[2]) +
             scalar_square(sig[5]) ) );
      tmp = c * cos(phi) / sin(phi);
      p_star = p + tmp;
      pc_star = pc + tmp;
      f_yield = q*q/(M*M) + p_star*(p_star - pc_star);
      f_flow = f_yield;
      if      ( task==GET_YIELD_RULE ) {
        if (  f_yield>f ) {
          f = f_yield;
          tmp_plasti_type = GROUP_MATERI_PLASTI_CAP1;
        }
        if ( swit ) pri( "f_yield", f_yield );
      }
      else {
        assert( task==GET_FLOW_RULE );
        f = f_flow;
        if ( swit ) pri( "f_flow", f_flow );
      }
    }
  }
  if ( get_group_data( GROUP_MATERI_PLASTI_HARDSOIL, gr, element, new_unknowns,
    plasti_data, length, GET_IF_EXISTS ) ) {
    test1 = task==GET_YIELD_RULE&&plasti_type==-NONE;
    test2 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_HARDSOIL;
    test3 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_HARDSOIL;
    if ( test1 || test2 || test3 ) {
      if ( swit ) pri( "check plasti_hardsoil" );
        // manual Professional 6.703 + theory "Hardening-Soil model":
        // 4 parameters phi c psi Rf. The yield function reads
        // f = q/(E50*(1 - q/qa)) - 2*q/Eur - gamma_p with q the
        // equivalent shear stress, qa = qf/Rf the asymptotic shear
        // stress and gamma_p the equivalent plastic shear strain
        // (here the materi_plasti_kappa dof, see below). E50 and Eur
        // are the power-law moduli of group_materi_elasti_hardsoil
        // evaluated at the CURRENT minor principal stress sig3
        // (coupled elasticity-plasticity, one-step lag documented).
        // The manual does not give the explicit qf; standard Schanz:
        // at Mohr-Coulomb failure in triaxial (q = sig1 - sig3,
        // sig1 = sig3*(1+sin(phi))/(1-sin(phi)) +
        // 2c*cos(phi)/(1-sin(phi))) qf = sig1 - sig3 =
        // 2*sin(phi)*(sig3 + c*cot(phi))/(1-sin(phi)) with sig3 the
        // minor principal stress in the manual (compression-positive)
        // convention: sig3_manual = -sig3_code + c*cot(phi),
        // sig3_code = largest algebraic eigenvalue (least compressive).
        // Flow: associative (the flow direction is the numerical
        // gradient of f, computed by the driver with central
        // differences: the law only defines f; psi is accepted for
        // input compatibility and documented as pending for the
        // non-associative flow potential).
      phi = plasti_data[0];
      c = plasti_data[1];
      hs_Rf = plasti_data[3];
      if ( sin(phi)<=0. || cos(phi)<=0. ) db_error( GROUP_MATERI_PLASTI_HARDSOIL, gr );
      if ( hs_Rf<=0. ) db_error( GROUP_MATERI_PLASTI_HARDSOIL, gr );
      if ( !get_group_data( GROUP_MATERI_ELASTI_HARDSOIL, gr, element,
          new_unknowns, hs_elasti, ldum, GET_IF_EXISTS ) ) {
        pri( "Error: GROUP_MATERI_PLASTI_HARDSOIL needs group_materi_elasti_hardsoil" );
        exit(TN_EXIT_STATUS);
      }
      if ( hs_elasti[0]<=0. || hs_elasti[4]<=0. )
        db_error( GROUP_MATERI_ELASTI_HARDSOIL, gr );
      hs_ccot = 0.;
      if ( c!=0. ) hs_ccot = c * cos(phi) / sin(phi);
        // minor principal stress: largest algebraic eigenvalue (the
        // manual's sig3, least compressive; see the elastic block)
      matrix_eigenvalues( sig, sig_princ );
      hs_sig3 = sig_princ[0];
      for ( i=1; i<MDIM; i++ )
        if ( sig_princ[i]>hs_sig3 ) hs_sig3 = sig_princ[i];
      hs_base = -hs_sig3 + hs_ccot;
      hs_E50 = hs_elasti[0];
      hs_Eur = hs_elasti[4];
      if ( hs_base>0. ) {
        if ( hs_elasti[1]+hs_ccot<=0. || hs_elasti[5]+hs_ccot<=0. )
          db_error( GROUP_MATERI_ELASTI_HARDSOIL, gr );
        hs_E50 = hs_elasti[0] *
          scalar_power( hs_base/(hs_elasti[1]+hs_ccot), hs_elasti[3] );
        hs_Eur = hs_elasti[4] *
          scalar_power( hs_base/(hs_elasti[5]+hs_ccot), hs_elasti[3] );
      }
      if ( hs_elasti[3]==0. ) { hs_E50 = hs_elasti[0]; hs_Eur = hs_elasti[4]; }
        // base <= 0 (sig3 very tensile beyond the cohesion): E = Eref
        // (same clamp as the elastic block) and the yield surface does
        // not exist -> no yielding from the hardsoil law (documented;
        // use a tension cut-off for those states)
      if ( hs_base<=0. ) {
        f_yield = NO_YIELD_F;
        f_flow = NO_YIELD_F;
      }
      else {
          // gamma_p: the accumulated plastic strain size kappa
          // (materi_plasti_kappa dof, kappa = int sqrt(0.5*deps_p:deps_p),
          // the same measure the manual calls "equivalent plastic shear
          // strain"; the yield function is calibrated to it, see
          // manual-developer) plus the optional extra initial
          // contribution from control_materi_plasti_hardsoil_gammap_initial
          // (stored per element by set_stress on the first timestep)
        gammap = 0.;
        if ( materi_plasti_kappa ) gammap = new_unknowns[kap_indx];
        hs_gp_extra = 0.;
        ldum = 0;
        db( ELEMENT_INTPNT_MATERI_PLASTI_HARDSOIL_GAMMAP_INITIAL, element, idum,
          &hs_gp_extra, ldum, VERSION_NORMAL, GET_IF_EXISTS );
        gammap += hs_gp_extra;
          // equivalent shear stress (same invariant as cap1)
        q = sqrt(
          0.5*( scalar_square(sig[0]-sig[4])+
                scalar_square(sig[4]-sig[8])+
                scalar_square(sig[0]-sig[8]) ) +
          3.*( scalar_square(sig[1]) +
               scalar_square(sig[2]) +
               scalar_square(sig[5]) ) );
        hs_qf = 2. * sin(phi) * hs_base / ( 1. - sin(phi) );
        hs_qa = hs_qf / hs_Rf;
          // q >= qa is the asymptote of the surface: the stress can
          // never reach qa, f -> +infinity approaching it. Set f to a
          // very large positive value (the driver's own "too large
          // time steps" exit threshold) so the return is driven hard;
          // if the plastic iterations cannot bring q below qa, the
          // driver exits with "too large time steps" (documented).
        if ( q>=hs_qa ) f_yield = DBL_MAX/1.e6;
        else f_yield = q / ( hs_E50 * ( 1. - q/hs_qa ) ) - 2.*q/hs_Eur - gammap;
        f_flow = f_yield;
      }
      if ( task==GET_YIELD_RULE ) {
        if (  f_yield>f ) {
          f = f_yield;
          tmp_plasti_type = GROUP_MATERI_PLASTI_HARDSOIL;
        }
        if ( swit ) pri( "f_yield", f_yield );
      }
      else {
        assert( task==GET_FLOW_RULE );
        f = f_flow;
        if ( swit ) pri( "f_flow", f_flow );
      }
    }
  }

  if ( get_group_data( GROUP_MATERI_PLASTI_COMPRESSION, gr, element, new_unknowns, 
    plasti_data, ldum, GET_IF_EXISTS ) ) {
    test1 = task==GET_YIELD_RULE&&plasti_type==-NONE;
    test2 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_COMPRESSION;
    test3 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_COMPRESSION;
    if ( test1 || test2 || test3 ) {
      if ( swit ) pri( "check plasti_compression" );
      sig_yield = plasti_data[0];
      if ( swit ) pri( "sig_yield", sig_yield );
      matrix_eigenvalues( sig, sig_princ );
      if ( swit ) pri( "sig_princ", sig_princ, MDIM );
      tmp_sig_eq = 0.;
      for ( idim=0; idim<MDIM; idim++ ) {
        if ( sig_princ[idim]<0. ) {
          tmp_sig_eq += sig_princ[idim] * sig_princ[idim];
        }
      }
      tmp_sig_eq = sqrt( scalar_dabs(tmp_sig_eq) ); 
      f_yield = tmp_sig_eq - sig_yield;
      f_flow = tmp_sig_eq - sig_yield;
      if      ( task==GET_YIELD_RULE ) {
        if (  f_yield>f ) {
          f = f_yield;
          tmp_plasti_type = GROUP_MATERI_PLASTI_COMPRESSION;
        }
        if ( swit ) pri( "f_yield", f_yield );
      }
      else {
        assert( task==GET_FLOW_RULE );
        f = f_flow;
        if ( swit ) pri( "f_flow", f_flow );
      }
    }
  }
  if ( get_group_data( GROUP_MATERI_PLASTI_DIPRISCO, gr, element, new_unknowns, 
    plasti_data, ldum, GET_IF_EXISTS ) ) {
    test1 = task==GET_YIELD_RULE&&plasti_type==-NONE;
    test2 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_DIPRISCO;
    test3 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_DIPRISCO;
    if ( test1 || test2 || test3 ) {
      get_group_data( GROUP_MATERI_PLASTI_DIPRISCO_RT, gr, element, new_unknowns, 
        &prisco_rt, ldum, GET_IF_EXISTS );
        // geotechnical notation
      array_multiply( sig, geo_sig, -1., MDIM*MDIM );
      array_multiply( inc_epp, geo_inc_epp, -1., MDIM*MDIM );
      array_multiply( old_epp, geo_old_epp, -1., MDIM*MDIM );
      if ( swit ) pri( "plasti_diprisco" );
      assert( materi_history_variables>=11 );
      sigm = ( geo_sig[0] + geo_sig[4] + geo_sig[8] ) / 3.;
      array_move( geo_sig, sig_dev, MDIM*MDIM );
      for ( idim=0; idim<MDIM; idim++ ) sig_dev[idim*MDIM+idim] -= sigm;
      if ( swit ) pri( "sig_dev", sig_dev, MDIM, MDIM );
        // history variables
      array_move( old_hisv, prisco_chi, MDIM*MDIM );
      tmp = array_size( prisco_chi, MDIM*MDIM );
      assert( tmp!=0. );
      array_multiply( prisco_chi, prisco_chi, 1./tmp, MDIM*MDIM );
      prisco_beta = old_hisv[MDIM*MDIM];
      prisco_rc = - old_hisv[MDIM*MDIM+1];
        // user data
      prisco_gamma = plasti_data[0];
      prisco_betafhat = plasti_data[1];
      prisco_bp = plasti_data[2];
      prisco_cp = plasti_data[3];
      prisco_tp = plasti_data[4];
      prisco_thetahat_c = plasti_data[5];
      prisco_thetahat_e = plasti_data[6];
      prisco_ksi_c = plasti_data[7];
      prisco_ksi_e = plasti_data[8];
      prisco_betaf0 = plasti_data[9];
        // fixed data
      prisco_r = array_inproduct( geo_sig, prisco_chi, MDIM*MDIM );
      if ( prisco_bp<=0. ) db_error( GROUP_MATERI_PLASTI_DIPRISCO, gr );
      if ( swit ) {
        pri( "prisco_chi", prisco_chi, MDIM*MDIM );
        pri( "prisco_betaf", prisco_betaf );
        pri( "prisco_rc", prisco_rc );
        pri( "prisco_gamma", prisco_gamma );
        pri( "prisco_betafhat", prisco_betafhat );
        pri( "prisco_bp", prisco_bp );
        pri( "prisco_cp", prisco_cp );
        pri( "prisco_tp", prisco_tp );
        pri( "prisco_thetahat_c", prisco_thetahat_c );
        pri( "prisco_thetahat_e", prisco_thetahat_e );
        pri( "prisco_ksi_c", prisco_ksi_c );
        pri( "prisco_ksi_e", prisco_ksi_e );
        pri( "prisco_r", prisco_r );
      }
      for ( idim=0; idim<MDIM; idim++ ) {
        for ( jdim=0; jdim<MDIM; jdim++ ) {
          ind = idim*MDIM+jdim;
          prisco_ss[ind] = geo_sig[ind] - prisco_r * prisco_chi[ind];
          if ( scalar_dabs(prisco_r)<EPS_PRISCO_1 ) {
            prisco_etas[ind] = 0.;
          }
          else
            prisco_etas[ind] = prisco_ss[ind] * sqrt(3.) / (prisco_r+prisco_rt); 
        }
      }
      if ( swit ) {
        pri( "prisco_ss", prisco_ss, MDIM, MDIM );
        pri( "prisco_etas", prisco_etas, MDIM, MDIM );
      }
      prisco_j2eta = array_inproduct( prisco_etas, prisco_etas, MDIM*MDIM );
      matrix_ab( prisco_etas, prisco_etas, work, MDIM, MDIM, MDIM );
      prisco_j3eta = array_inproduct( work, prisco_etas, MDIM*MDIM );
      if ( swit ) {
        pri( "prisco_j2eta", prisco_j2eta );
        pri( "prisco_j3eta", prisco_j3eta );
      }
        // betaf calculation
      prisco_betaf = prisco_betafhat + ( prisco_betaf0 - prisco_betafhat ) *
        exp( - prisco_beta * prisco_tp );
        // second yield function
      prisco_r1tmp = (prisco_rc/sqrt(3.))*EPS_PRISCO_2;
      prisco_r2tmp = 0.;
      for ( idim=0; idim<MDIM; idim++ ) {
        for ( jdim=0; jdim<MDIM; jdim++ ) {
          ind = idim*MDIM+jdim;
          if ( prisco_etas[ind]>prisco_r2tmp ) prisco_r2tmp = prisco_etas[ind];
        }
      }
      if ( prisco_r<prisco_r1tmp && task==GET_YIELD_RULE && prisco_rt==0. ) {
        f_yield = NO_YIELD_F;
        f_flow = NO_YIELD_F;
      }
      else if ( prisco_r2tmp >= 2.7 && task==GET_YIELD_RULE ) {
        f_yield = NO_YIELD_F;
        f_flow = NO_YIELD_F;
      }
      else {
          // normal yield rule and flow rule
        tmp = (prisco_r+prisco_rt)/(prisco_rc+prisco_rt);
        f_yield = 3. * prisco_betaf * ( prisco_gamma - 3. ) * 
          log ( scalar_dabs(tmp) ) -
          prisco_gamma * prisco_j3eta + (9./4.) * ( prisco_gamma - 1. ) * prisco_j2eta;
        f_flow  = 9. * ( prisco_gamma - 3. ) * 
          log ( scalar_dabs(tmp) ) -
          prisco_gamma * prisco_j3eta + (9./4.) * ( prisco_gamma - 1. ) * prisco_j2eta;
      }
      tmp = ( prisco_etas[0] + prisco_etas[4] + prisco_etas[8] ) ;
        // load angle
      array_set( &prisco_st1[0][0], 0., MDIM*MDIM );
      array_set( &prisco_st2[0][0], 0., MDIM*MDIM );
      prisco_st1[0][0] = 0.8164965809;
      prisco_st1[1][1] = - prisco_st1[0][0] / 2.;
      prisco_st1[2][2] = - prisco_st1[0][0] / 2.;
      tmp = array_size( geo_sig, MDIM*MDIM );
      if ( tmp!=0. )
        array_multiply( geo_sig, &prisco_sv0[0][0], 1./tmp, MDIM*MDIM );
      else 
        array_set( &prisco_sv0[0][0], 0., MDIM*MDIM );

      prisco_r0 = prisco_chi[0]+prisco_chi[4]+prisco_chi[8];
      prisco_st0 = prisco_sv0[0][0]+prisco_sv0[1][1]+
                   prisco_sv0[2][2];

      for ( idim=0; idim<MDIM; idim++ ) {
        for ( jdim=0; jdim<MDIM; jdim++ ) {
          ind = idim*MDIM+jdim;
          prisco_rv[idim][jdim] = prisco_chi[ind];
          prisco_st[idim][jdim] = prisco_sv0[idim][jdim];
          if ( idim==jdim ) {
            prisco_rv[idim][jdim] -= prisco_r0/3.;
            prisco_st[idim][jdim] -= prisco_st0/3.;
          }
          prisco_rr[idim][jdim] = prisco_st[idim][jdim] - prisco_rv[idim][jdim];
        }
      }
      if ( swit ) {
        pri( "prisco_betaf", prisco_betaf );
        pri( "prisco_r0", prisco_r0 );
        pri( "prisco_st0", prisco_st0 );
        pri( "prisco_sv0", &prisco_sv0[0][0], MDIM*MDIM );
        pri( "prisco_rv", &prisco_rv[0][0], MDIM*MDIM );
        pri( "prisco_st", &prisco_st[0][0], MDIM*MDIM );
        pri( "prisco_rr", &prisco_rr[0][0], MDIM*MDIM );
      }
      tmp = array_size( &prisco_rr[0][0], MDIM*MDIM );
      if ( tmp<EPS_PRISCO_1 ) tmp = EPS_PRISCO_1;
      array_multiply( &prisco_rr[0][0], &prisco_rr[0][0], 1./tmp, MDIM*MDIM );
//        alternative method
//      matrix_eigenvalues( &prisco_rr[0][0], work );
//      prisco_st2[0][0] = work[0];
//      prisco_st2[1][1] = work[1];
//      prisco_st2[2][2] = work[2];
      tmp1 = (prisco_rr[0][0]+prisco_rr[1][1])/2.;
      tmp2 = (prisco_rr[0][0]-prisco_rr[1][1])/2.;
      tmp2 = tmp2*tmp2+prisco_rr[0][1]*prisco_rr[0][1];
      tmp2 = sqrt(tmp2);
      prisco_st2[0][0] = tmp1+tmp2;
      prisco_st2[1][1] = tmp1-tmp2;
      prisco_st2[2][2] = prisco_rr[2][2];

      prisco_b4 = array_inproduct( &prisco_st2[0][0], &prisco_st1[0][0], MDIM*MDIM );
      prisco_b7 = scalar_dabs(prisco_b4);
      if ( prisco_b7>1. ) prisco_b7 = 1.;
      if ( prisco_b4>1. ) prisco_b4 = 1.;
      if ( prisco_b4<0. ) prisco_b4 = -prisco_b7;
      prisco_b5 = acos(prisco_b4); 
      if ( prisco_b5<(1./3)*PIRAD ) {
        prisco_b5 = scalar_dabs(prisco_b5);
      }
      else {
        if ( prisco_b5<(3./3.)*PIRAD ) {
          prisco_b5 = prisco_b5 - (2./3.)*PIRAD;
          prisco_b5 = scalar_dabs(prisco_b5);
        }
        else {
          if ( prisco_b5<(5./3.)*PIRAD ) {
            prisco_b5 = prisco_b5 - (4./3.)*PIRAD;
            prisco_b5 = scalar_dabs(prisco_b5);
          }
        }
      }
        // ksi, theta
      prisco_ksi = 0.5*(prisco_ksi_c+prisco_ksi_e) +
        0.5*(prisco_ksi_c-prisco_ksi_e)*(2.*prisco_b5*prisco_b5-4.*prisco_b5+1.);
      prisco_thetahat = 0.5*(prisco_thetahat_c+prisco_thetahat_e) +
        0.5*(prisco_thetahat_c-prisco_thetahat_e)*(2.*prisco_b5*prisco_b5-4.*prisco_b5+1.);
      prisco_theta = prisco_thetahat;
      if ( swit ) {
        pri( "prisco_b4", prisco_b4 );
        pri( "prisco_b5", prisco_b5 );
        pri( "prisco_b7", prisco_b7 );
        pri( "prisco_kinc", prisco_kinc );
        pri( "prisco_ksi", prisco_ksi );
        pri( "prisco_thetahat", prisco_thetahat );
        pri( "prisco_theta", prisco_theta );
      }
      tmp1 = array_inproduct( geo_inc_epp, prisco_chi, MDIM*MDIM );
      for ( idim=0; idim<MDIM; idim++ ) {
        for ( jdim=0; jdim<MDIM; jdim++ ) {
          ind = idim*MDIM+jdim;
          prisco_epinc[ind] = geo_inc_epp[ind] - tmp1 * prisco_chi[ind];
        }
      }
      tmp2 = prisco_ksi * array_size( prisco_epinc, MDIM*MDIM );
      prisco_rc += (prisco_rc/prisco_bp) * ( tmp1 + tmp2 );
      if (prisco_rc <= EPS_PRISCO_2) prisco_rc = EPS_PRISCO_2;
      if ( swit ) {
        pri( "prisco_epinc", prisco_epinc, MDIM, MDIM );
        pri( "tmp1", tmp1 );
        pri( "tmp2", tmp2 );
        pri( "prisco_rc", prisco_rc );
      }
      for ( idim=0; idim<MDIM; idim++ ) {
        for ( jdim=0; jdim<MDIM; jdim++ ) {
          ind = idim*MDIM+jdim;
          prisco_chihat[ind] = prisco_rr[idim][jdim] * sin( prisco_theta );
          if ( idim==jdim ) prisco_chihat[ind] += cos( prisco_theta ) / sqrt( 3. );
        }
      }
      prisco_kinc = prisco_cp * array_size( geo_inc_epp, MDIM*MDIM );
      for ( idim=0; idim<MDIM; idim++ ) {
        for ( jdim=0; jdim<MDIM; jdim++ ) {
          ind = idim*MDIM+jdim;
          prisco_chi[ind] += prisco_kinc * ( prisco_chihat[ind] - prisco_chi[ind] );
        }
      }
      array_normalize( prisco_chi, MDIM*MDIM );
      if ( swit ) {
        pri( "prisco_chihat", prisco_chihat, MDIM, MDIM );
        pri( "prisco_chi", prisco_chi, MDIM, MDIM );
      }
      array_move( geo_old_epp, work, MDIM*MDIM );
      array_add( work, geo_inc_epp, work, MDIM*MDIM );
      tmp1 = (work[0] + work[4] + work[8])/3. ;
      tmp2 = (geo_inc_epp[0] + geo_inc_epp[4] + geo_inc_epp[8])/3. ;
      array_move( work, prisco_epp_dev, MDIM*MDIM );
      array_move( geo_inc_epp, prisco_incepp_dev, MDIM*MDIM );
      for ( idim=0; idim<MDIM; idim++ ) {
        ind = idim*MDIM + idim;
        prisco_epp_dev[ind] -= tmp1;
        prisco_incepp_dev[ind] -= tmp2;
      }

      if ( swit ) {
        pri( "prisco_epp_dev", prisco_epp_dev, MDIM*MDIM );
        pri( "prisco_incepp_dev", prisco_incepp_dev, MDIM*MDIM );
      }
      tmp1 = array_inproduct( prisco_epp_dev, prisco_incepp_dev, MDIM*MDIM );
      tmp2 = array_inproduct( prisco_epp_dev, prisco_epp_dev, MDIM*MDIM );
      tmp2 = sqrt( scalar_dabs( tmp2 ) ); 
      if ( tmp2<EPS_PRISCO_1 ) tmp2 = EPS_PRISCO_1;
      tmp = tmp1 / tmp2;
      prisco_beta += tmp;
      if ( swit ) pri( "prisco_beta new", prisco_beta );
      array_move( prisco_chi, new_hisv, MDIM*MDIM );
      new_hisv[MDIM*MDIM] = prisco_beta;
      new_hisv[MDIM*MDIM+1] = - prisco_rc;
      if      ( task==GET_YIELD_RULE ) {
        if ( f_yield>f ) {
          f = f_yield;
          tmp_plasti_type = GROUP_MATERI_PLASTI_DIPRISCO;
        }
        if ( swit ) {
          pri( "f_yield", f_yield );
          pri( "new_hisv", new_hisv, materi_history_variables );
        }
      }
      else {
        assert( task==GET_FLOW_RULE );
        f = f_flow;
        if ( swit ) pri( "f_flow", f_flow );
      }
    }
  }
  long int mc_hs_direct = 0;
  if ( get_group_data( GROUP_MATERI_PLASTI_MOHR_COUL_HARDENING_SOFTENING, gr, element, new_unknowns,
      plasti_data, ldum, GET_IF_EXISTS ) ||
       ( ( mc_hs_direct = get_group_data( GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_HARDENING_SOFTENING,
         gr, element, new_unknowns, plasti_data, ldum, GET_IF_EXISTS ) ) ) ) {
    test1 = task==GET_YIELD_RULE&&plasti_type==-NONE;
    test2 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_MOHR_COUL_HARDENING_SOFTENING;
    test3 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_MOHR_COUL_HARDENING_SOFTENING;
    if ( test1 || test2 || test3 ) {
      if ( swit ) pri( "plasti_mohr_coul_hardening_softening" );
      if ( mc_hs_direct ) {
        // direct variant: the angles are in degrees (like every
        // group_materi_plasti_mohr_coul_direct record)
        plasti_data[0] *= PIRAD/180.;
        plasti_data[2] *= PIRAD/180.;
        plasti_data[3] *= PIRAD/180.;
        plasti_data[5] *= PIRAD/180.;
      }
      // manual Professional 6.731: same surface as mohr_coul, but c and
      // phi (both yield and flow) vary LINEARLY with the effective
      // plastic strain kappa_shear from the _0 values at kappa=0 up to
      // the _1 values at kappa_shear_crit, constant afterwards. The
      // accumulated kappa is the materi_plasti_kappa dof
      // (new_unknowns[kap_indx], a node dof filled by the standard
      // driver; kappa = int sqrt(0.5*depsp:depsp)).
      {
        double phi_0 = plasti_data[0], c_0 = plasti_data[1],
               phi_flow_0 = plasti_data[2], phi_1 = plasti_data[3],
               c_1 = plasti_data[4], phi_flow_1 = plasti_data[5],
               kappa_crit = plasti_data[6], kappa_now = 0., ratio;
        // manual Professional 6.731: the hardening variable is the
        // SHEAR plastic strain kappa_shear (materi_plasti_kappa_shear
        // dof kapsh, manual 4.25) - the size of the deviatoric plastic
        // strain. The total kappa (materi_plasti_kappa, 4.24) is used
        // only when the shear dof is not initialized.
        if ( materi_plasti_kappa_shear )
          kappa_now = new_unknowns[kapsh_indx];
        else if ( materi_plasti_kappa )
          kappa_now = new_unknowns[kap_indx];
        if ( kappa_crit>0. ) {
          ratio = kappa_now/kappa_crit;
          if ( ratio>1. ) ratio = 1.;
          if ( ratio<0. ) ratio = 0.;
        }
        else
          ratio = 1.;
        phi = phi_0 + ratio*(phi_1-phi_0);
        c = c_0 + ratio*(c_1-c_0);
        phi_flow = phi_flow_0 + ratio*(phi_flow_1-phi_flow_0);
      }
      if ( plasti_on_boundary ) {
        phi *= plasti_on_boundary_factor;
        phi_flow *= plasti_on_boundary_factor;
      }
      matrix_eigenvalues( sig, sig_princ );
      {
        double sig_max = sig_princ[0], sig_min = sig_princ[0];
        long int ipr;
        for ( ipr=1; ipr<MDIM; ipr++ ) {
          if ( sig_princ[ipr]>sig_max ) sig_max = sig_princ[ipr];
          if ( sig_princ[ipr]<sig_min ) sig_min = sig_princ[ipr];
        }
        f_yield = 0.5*(sig_max-sig_min) + 0.5*(sig_max+sig_min)*sin(phi)
                 - c*cos(phi);
        f_flow  = 0.5*(sig_max-sig_min) + 0.5*(sig_max+sig_min)*sin(phi_flow)
                 - c*cos(phi);
      }
      if      ( task==GET_YIELD_RULE ) {
        if ( f_yield>f ) {
          f = f_yield;
          tmp_plasti_type = GROUP_MATERI_PLASTI_MOHR_COUL_HARDENING_SOFTENING;
        }
        if ( swit ) pri( "f_yield", f_yield );
      }
      else {
        assert( task==GET_FLOW_RULE );
        f = f_flow;
        if ( swit ) pri( "f_flow", f_flow );
      }
    }
  }
  if ( get_group_data( GROUP_MATERI_PLASTI_MOHR_COUL, gr, element, new_unknowns,
      plasti_data, ldum, GET_IF_EXISTS ) ) {
    test1 = task==GET_YIELD_RULE&&plasti_type==-NONE;
    test2 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_MOHR_COUL;
    test3 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_MOHR_COUL;
    if ( test1 || test2 || test3 ) {
      if ( swit ) pri( "plasti_mohr_coul" );
      phi = plasti_data[0];
      c = plasti_data[1];
      phi_flow = plasti_data[2];
      if ( plasti_on_boundary ) {
        phi *= plasti_on_boundary_factor;
        phi_flow *= plasti_on_boundary_factor;
      }
      // Mohr-Coulomb (manual Professional: 0.5(sig1-sig3) +
      // 0.5(sig1+sig3) sin(phi) - c cos(phi) = 0, sig1 the LARGEST and
      // sig3 the SMALLEST principal stress, tension-positive). The flow
      // rule uses phi_flow (non-associative when different).
      matrix_eigenvalues( sig, sig_princ );
      {
        double sig_max = sig_princ[0], sig_min = sig_princ[0];
        long int ipr;
        for ( ipr=1; ipr<MDIM; ipr++ ) {
          if ( sig_princ[ipr]>sig_max ) sig_max = sig_princ[ipr];
          if ( sig_princ[ipr]<sig_min ) sig_min = sig_princ[ipr];
        }
        f_yield = 0.5*(sig_max-sig_min) + 0.5*(sig_max+sig_min)*sin(phi)
                 - c*cos(phi);
        f_flow  = 0.5*(sig_max-sig_min) + 0.5*(sig_max+sig_min)*sin(phi_flow)
                 - c*cos(phi);
      }
      if      ( task==GET_YIELD_RULE ) {
        if ( f_yield>f ) {
          f = f_yield;
          tmp_plasti_type = GROUP_MATERI_PLASTI_MOHR_COUL;
        }
        if ( swit ) pri( "f_yield", f_yield );
      }
      else {
        assert( task==GET_FLOW_RULE );
        f = f_flow;
        if ( swit ) pri( "f_flow", f_flow );
      }
    }
  }
  if ( get_group_data( GROUP_MATERI_PLASTI_DRUCKPRAG, gr, element, new_unknowns,
      plasti_data, ldum, GET_IF_EXISTS ) ) {
    tension_cutoff = -NO;
    db( GROUP_MATERI_PLASTI_DRUCKPRAG_TENSIONCUTOFF, 0, &tension_cutoff, ddum, ldum, 
      VERSION_NORMAL, GET_IF_EXISTS );
    tension_limit = -NO;
    db( GROUP_MATERI_PLASTI_DRUCKPRAG_TENSIONLIMIT, 0, &tension_limit, ddum, ldum, 
      VERSION_NORMAL, GET_IF_EXISTS );
    test1 = task==GET_YIELD_RULE&&plasti_type==-NONE;
    test2 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_DRUCKPRAG;
    test3 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_DRUCKPRAG_APEX;
    test4 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_DRUCKPRAG;
    test5 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_DRUCKPRAG_APEX;
    if ( test1 || test2 || test3 || test4 || test5 ) {
      if ( swit ) pri( "check plasti_druckprag" );
      phi = plasti_data[0];
      phi_flow = plasti_data[2];
      if ( plasti_on_boundary ) {
        phi *= plasti_on_boundary_factor;
        phi_flow *= plasti_on_boundary_factor;
      }
      c = plasti_data[1];
      alpha = 2.*sin(phi)/((3.-sin(phi))*sqrt(3.));
      alpha_flow = 2.*sin(phi_flow)/((3.-sin(phi_flow))*sqrt(3.));
      K = 6.*c*cos(phi)/((3.-sin(phi))*sqrt(3.));
      K_flow = 6.*c*cos(phi_flow)/((3.-sin(phi_flow))*sqrt(3.));
      test1 = task==GET_YIELD_RULE&&plasti_type==-NONE && phi>0.;
      test2 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_DRUCKPRAG_APEX;
      test3 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_DRUCKPRAG_APEX;
      if ( test1 || test2 || test3 ) {
          // separate yield and flow rule for apex
        matrix_eigenvalues( sig, sig_princ );
        array_set( sig_princ_apex, c*cos(phi)/sin(phi), MDIM );
        array_set( e_vec, 1./sqrt(3.), MDIM );
        array_subtract( sig_princ, sig_princ_apex, f_vec, MDIM );
        in_apex = array_inproduct(f_vec,e_vec,MDIM)>0. && tension_cutoff==-YES;
        if ( in_apex || plasti_type==GROUP_MATERI_PLASTI_DRUCKPRAG_APEX ) {
          f_yield = array_size( f_vec, MDIM );
          f_flow = array_size( f_vec, MDIM );
        }
        else {
          f_yield = NO_YIELD_F;
          f_flow = NO_YIELD_F;
        }
        if      ( task==GET_YIELD_RULE ) {
          if (  f_yield>f ) {
            f = f_yield;
            tmp_plasti_type = GROUP_MATERI_PLASTI_DRUCKPRAG_APEX;
          }
          if ( swit ) pri( "f_yield", f_yield );
        }
        else {
          assert( task==GET_FLOW_RULE );
          f = f_flow;
          if ( swit ) pri( "f_flow", f_flow );
        }
      }
      if ( !in_apex ) {
        test1 = task==GET_YIELD_RULE&&plasti_type==-NONE;
        test2 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_DRUCKPRAG;
        test3 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_DRUCKPRAG;
        if ( test1 || test2 || test3 ) {
          if ( tension_limit==-YES && sigm>0. ) {
            f_yield = sig_eq - K;
            f_flow = sig_eq - K_flow;
          }
          else {
            f_yield = 3.*alpha*sigm + sig_eq - K;
            f_flow = 3.*alpha_flow*sigm + sig_eq - K_flow;
          }
          if      ( task==GET_YIELD_RULE ) {
            if (  f_yield>f ) {
              f = f_yield;
              tmp_plasti_type = GROUP_MATERI_PLASTI_DRUCKPRAG;
            }
            if ( swit ) pri( "f_yield", f_yield );
          }
          else {
            assert( task==GET_FLOW_RULE );
            f = f_flow;
            if ( swit ) pri( "f_flow", f_flow );
          }
        }
      }
    }
  }
  if ( get_group_data( GROUP_MATERI_PLASTI_GURSON, gr, element, new_unknowns,
      plasti_data, ldum, GET_IF_EXISTS ) ) {
    test1 = task==GET_YIELD_RULE&&plasti_type==-NONE;
    test2 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_GURSON;
    test3 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_GURSON;
    if ( test1 || test2 || test3 ) {
      if ( swit ) pri( "plasti_gurson" );
      void_fraction = new_unknowns[void_indx];
      sig_yield = plasti_data[0];
      q1 = plasti_data[1];
      q2 = plasti_data[2];
      q3 = plasti_data[3];
      f_yield = (3.*sig_eq*sig_eq)/(sig_yield*sig_yield) +
        2.*q1*void_fraction*cosh(q2*3.*sigm/(2.*sig_yield)) -
        (1.+(q3*void_fraction)*(q3*void_fraction) );
      f_flow = f_yield;
      if      ( task==GET_YIELD_RULE ) {
        if ( f_yield>f ) {
          f = f_yield;
          tmp_plasti_type = GROUP_MATERI_PLASTI_GURSON;
        }
        if ( swit ) pri( "f_yield", f_yield );
      }
      else {
        assert( task==GET_FLOW_RULE );
        f = f_flow;
        if ( swit ) pri( "f_flow", f_flow );
      }
    }
  }
  if ( get_group_data( GROUP_MATERI_PLASTI_HLC, gr, element, new_unknowns,
      plasti_data, ldum, GET_IF_EXISTS ) ) {
    test1 = task==GET_YIELD_RULE&&plasti_type==-NONE;
    test2 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_HLC;
    test3 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_HLC;
    if ( test1 || test2 || test3 ) {
      if ( swit ) pri( "plasti_hlc" );
      void_fraction = new_unknowns[void_indx];
      sig_yield = plasti_data[0];
      m21 = plasti_data[1];
      m22 = plasti_data[2];
      m23 = plasti_data[3];
      m1  = plasti_data[4];
      m20 = (m21 - m22 * void_fraction) * exp(m23*sigm/sig_yield);
      f_yield = (3.*sig_eq*sig_eq)/(sig_yield*sig_yield) +
        (1+1/m20)*void_fraction*m1*exp(3.*sigm/(2.*sig_yield)) - 1.;
      f_flow = f_yield;
      if      ( task==GET_YIELD_RULE ) {
        if ( f_yield>f ) {
          f = f_yield;
          tmp_plasti_type = GROUP_MATERI_PLASTI_HLC;
        }
        if ( swit ) pri( "f_yield", f_yield );
      }
      else {
        assert( task==GET_FLOW_RULE );
        f = f_flow;
        if ( swit ) pri( "f_flow", f_flow );
      }
    }
  }
  if ( get_group_data( GROUP_MATERI_PLASTI_MATSUOKANAKAI, gr, element, new_unknowns,
      plasti_data, ldum, GET_IF_EXISTS ) ) {
    test1 = task==GET_YIELD_RULE&&plasti_type==-NONE;
    test2 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_MATSUOKANAKAI;
    test3 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_MATSUOKANAKAI_APEX;
    test4 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_MATSUOKANAKAI;
    test5 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_MATSUOKANAKAI_APEX;
    if ( test1 || test2 || test3 || test4 || test5 ) {
      if ( swit ) pri( "check plasti_matsuokanakai" );
      phi = plasti_data[0];
      c = plasti_data[1];
      phi_flow = plasti_data[2];
      if ( plasti_on_boundary ) {
        phi *= plasti_on_boundary_factor;
        phi_flow *= plasti_on_boundary_factor;
      }
      alpha_flow = 2.*sin(phi_flow)/((3.-sin(phi_flow))*sqrt(3.));
      if ( phi>0. )
         sig0 = c * cos(phi) / sin(phi);
      else
         sig0 = 0.;
      array_move( sig, sig_tmp, MDIM*MDIM );
      for ( idim=0; idim<MDIM; idim++ ) sig_tmp[idim*MDIM+idim] -= sig0;
      matrix_eigenvalues( sig_tmp, sig_princ );
      I1 = sig_princ[0] + sig_princ[1] + sig_princ[2];
      I2 = -sig_princ[0]*sig_princ[1] - sig_princ[1]*sig_princ[2] - 
        sig_princ[2]*sig_princ[0];
      I3 = sig_princ[0] * sig_princ[1] * sig_princ[2];
      test1 = task==GET_YIELD_RULE&&plasti_type==-NONE && phi>0.;
      test2 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_MATSUOKANAKAI_APEX;
      test3 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_MATSUOKANAKAI_APEX;
      if ( test1 || test2 || test3 ) {
          // separate yield and flow rule for apex
        array_set( sig_princ_apex, 0., MDIM );
        array_set( e_vec, 1./sqrt(3.), MDIM );
        array_subtract( sig_princ, sig_princ_apex, f_vec, MDIM );
        in_apex = array_inproduct(f_vec,e_vec,MDIM)>0. && tension_cutoff==-YES;
        if ( in_apex || plasti_type==GROUP_MATERI_PLASTI_MATSUOKANAKAI_APEX ) {
          f_yield = array_size( f_vec, MDIM );
          f_flow = array_size( f_vec, MDIM );
        }
        else {
          f_yield = NO_YIELD_F;
          f_flow = NO_YIELD_F;
        }
        if      ( task==GET_YIELD_RULE ) {
          if (  f_yield>f ) {
            f = f_yield;
            tmp_plasti_type = GROUP_MATERI_PLASTI_MATSUOKANAKAI_APEX;
          }
          if ( swit ) pri( "f_yield", f_yield );
        }
        else {
          assert( task==GET_FLOW_RULE );
          f = f_flow;
          if ( swit ) pri( "f_flow", f_flow );
        }
      }
      if ( !in_apex ) {
        test1 = task==GET_YIELD_RULE&&plasti_type==-NONE;
        test2 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_MATSUOKANAKAI;
        test3 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_MATSUOKANAKAI;
        if ( test1 || test2 || test3 ) {
          f_yield = I3 + scalar_square(cos(phi))/(9.-scalar_square(sin(phi))) * I1*I2;
          f_flow = 3.*alpha_flow*sigm + sig_eq;
          if      ( task==GET_YIELD_RULE ) {
            if (  f_yield>f ) {
              f = f_yield;
              tmp_plasti_type = GROUP_MATERI_PLASTI_MATSUOKANAKAI;
            }
            if ( swit ) pri( "f_yield", f_yield );
          }
          else {
            assert( task==GET_FLOW_RULE );
            f = f_flow;
            if ( swit ) pri( "f_flow", f_flow );
          }
        }
      }
    }
  }
  if ( get_group_data( GROUP_MATERI_PLASTI_MOHRCOUL, gr, element, new_unknowns,
      plasti_data, ldum, GET_IF_EXISTS ) ||
       get_group_data( GROUP_MATERI_PLASTI_MOHRCOUL_SOFTENING, gr, element, new_unknowns,
      plasti_data, ldum, GET_IF_EXISTS ) ) {
    db( GROUP_MATERI_PLASTI_MOHRCOUL_TENSIONCUTOFF, 0, &tension_cutoff, ddum, ldum, 
      VERSION_NORMAL, GET_IF_EXISTS );
    test1 = task==GET_YIELD_RULE&&plasti_type==-NONE;
    test2 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_MOHRCOUL_01;
    test3 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_MOHRCOUL_12;
    test4 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_MOHRCOUL_20;
    test5 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_MOHRCOUL_APEX;
    test6 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_MOHRCOUL_01;
    test7 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_MOHRCOUL_12;
    test8 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_MOHRCOUL_20;
    test9 = GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_MOHRCOUL_APEX;
    if ( test1 || test2 || test3 || test4 || test5 || test6 || test7 || test8 || test9 ) {
      if ( swit ) pri( "check plasti_mohrcoul" );
      if ( get_group_data( GROUP_MATERI_PLASTI_MOHRCOUL_SOFTENING, gr, element, 
          new_unknowns, plasti_data, ldum, GET_IF_EXISTS ) ) {
        phi0 = plasti_data[0];
        c0 = plasti_data[1];
        phi0_flow = plasti_data[2];
        phi1 = plasti_data[3];
        c1 = plasti_data[4];
        phi1_flow = plasti_data[5];
        kappa_crit = plasti_data[6];
        kappa = new_unknowns[kap_indx];
        if ( kappa>=kappa_crit ) {
          phi = phi1;
          c = c1;
          phi_flow = phi1_flow;
        }
        else {
          tmp = kappa/kappa_crit;
          phi = phi0 + tmp*(phi1-phi0);
          c = c0 + tmp*(c1-c0);
          phi_flow = phi0_flow + tmp*(phi1_flow-phi0_flow);
        }
      }
      else {
        phi = plasti_data[0];
        c = plasti_data[1];
        phi_flow = plasti_data[2];
      }
      if ( plasti_on_boundary ) {
        phi *= plasti_on_boundary_factor;
        phi_flow *= plasti_on_boundary_factor;
      }
      matrix_eigenvalues( sig, sig_princ );
      test1 = task==GET_YIELD_RULE&&plasti_type==-NONE && phi>0.;
      test2 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_MOHRCOUL_APEX;
      test3 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_MOHRCOUL_APEX;
      if ( test1 || test2 || test3 ) {
          // separate yield and flow rule for apex
        array_set( sig_princ_apex, c*cos(phi)/sin(phi), MDIM );
        array_set( e_vec, 1./sqrt(3.), MDIM );
        array_subtract( sig_princ, sig_princ_apex, f_vec, MDIM );
        in_apex = array_inproduct(f_vec,e_vec,MDIM)>0. && tension_cutoff==-YES;
        if ( in_apex || plasti_type==GROUP_MATERI_PLASTI_MOHRCOUL_APEX ) {
          f_yield = array_size( f_vec, MDIM );
          f_flow = array_size( f_vec, MDIM );
        }
        else {
          f_yield = NO_YIELD_F;
          f_flow = NO_YIELD_F;
        }
        if      ( task==GET_YIELD_RULE ) {
          if (  f_yield>f ) {
            f = f_yield;
            tmp_plasti_type = GROUP_MATERI_PLASTI_MOHRCOUL_APEX;
          }
          if ( swit ) pri( "f_yield", f_yield );
        }
        else {
          assert( task==GET_FLOW_RULE );
          f = f_flow;
          if ( swit ) pri( "f_flow", f_flow );
        }
      }
      if ( !in_apex ) {
        test1 = task==GET_YIELD_RULE&&plasti_type==-NONE;
        test2 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_MOHRCOUL_01;
        test3 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_MOHRCOUL_01;
        if ( test1 || test2 || test3 ) {
          f_yield = 0.5*scalar_dabs(sig_princ[0]-sig_princ[1]) +
           0.5*(sig_princ[0]+sig_princ[1])*sin(phi) - c*cos(phi);
          f_flow = 0.5*scalar_dabs(sig_princ[0]-sig_princ[1]) + 
            0.5*(sig_princ[0]+sig_princ[1])*sin(phi_flow) - c*cos(phi_flow);
          if      ( task==GET_YIELD_RULE ) {
            if ( f_yield>f ) {
              f = f_yield;
              tmp_plasti_type = GROUP_MATERI_PLASTI_MOHRCOUL_01;
              if ( swit ) pri( "f_yield", f_yield );
            }
          }
          else {
            assert( task==GET_FLOW_RULE );
            f = f_flow;
            if ( swit ) pri( "f_flow", f_flow );
          }
        }
        test1 = task==GET_YIELD_RULE&&plasti_type==-NONE;
        test2 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_MOHRCOUL_12;
        test3 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_MOHRCOUL_12;
        if ( test1 || test2 || test3 ) {
          f_yield = 0.5*scalar_dabs(sig_princ[1]-sig_princ[2]) +
           0.5*(sig_princ[1]+sig_princ[2])*sin(phi) - c*cos(phi);
          f_flow = 0.5*scalar_dabs(sig_princ[1]-sig_princ[2]) + 
            0.5*(sig_princ[1]+sig_princ[2])*sin(phi_flow) - c*cos(phi_flow);
          if      ( task==GET_YIELD_RULE ) {
            if ( f_yield>f ) {
              f = f_yield;
              tmp_plasti_type = GROUP_MATERI_PLASTI_MOHRCOUL_12;
              if ( swit ) pri( "f_yield", f_yield );
            }
          }
          else {
            assert( task==GET_FLOW_RULE );
            f = f_flow;
            if ( swit ) pri( "f_flow", f_flow );
          }
        }
        test1 = task==GET_YIELD_RULE&&plasti_type==-NONE;
        test2 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_MOHRCOUL_20;
        test3 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_MOHRCOUL_20;
        if ( test1 || test2 || test3 ) {
          f_yield = 0.5*scalar_dabs(sig_princ[2]-sig_princ[0]) +
           0.5*(sig_princ[2]+sig_princ[0])*sin(phi) - c*cos(phi);
          f_flow = 0.5*scalar_dabs(sig_princ[2]-sig_princ[0]) + 
            0.5*(sig_princ[2]+sig_princ[0])*sin(phi_flow) - c*cos(phi_flow);
          if      ( task==GET_YIELD_RULE ) {
            if ( f_yield>f ) {
              f = f_yield;
              tmp_plasti_type = GROUP_MATERI_PLASTI_MOHRCOUL_20;
              if ( swit ) pri( "f_yield", f_yield );
            }
          }
          else {
            assert( task==GET_FLOW_RULE );
            f = f_flow;
            if ( swit ) pri( "f_flow", f_flow );
          }
        }
      }
    }
  }
  if ( get_group_data( GROUP_MATERI_PLASTI_STRESS, gr, element, new_unknowns, 
      plasti_data, ldum, GET_IF_EXISTS ) ) {
    test1 = task==GET_YIELD_RULE&&plasti_type==-NONE;
    test2 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_STRESS;
    test3 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_STRESS;
    if ( test1 || test2 || test3 ) {
      if ( swit ) pri( "check plasti_stress" );
      sig_yield = plasti_data[0];
      if ( swit ) pri( "sig_yield", sig_yield );
      matrix_eigenvalues( sig, sig_princ );
      if ( swit ) pri( "sig_princ", sig_princ, MDIM );
      tmp_sig_eq = 0.;
      for ( idim=0; idim<MDIM; idim++ ) {
        tmp_sig_eq += sig_princ[idim] * sig_princ[idim];
      }
      tmp_sig_eq = sqrt( scalar_dabs(tmp_sig_eq) ); 
      f_yield = tmp_sig_eq - sig_yield;
      f_flow = tmp_sig_eq - sig_yield;
      if      ( task==GET_YIELD_RULE ) {
        if (  f_yield>f ) {
          f = f_yield;
          tmp_plasti_type = GROUP_MATERI_PLASTI_STRESS;
        }
        if ( swit ) pri( "f_yield", f_yield );
      }
      else {
        assert( task==GET_FLOW_RULE );
        f = f_flow;
        if ( swit ) pri( "f_flow", f_flow );
      }
    }
  }
  if ( get_group_data( GROUP_MATERI_PLASTI_TENSION, gr, element, new_unknowns, 
      plasti_data, ldum, GET_IF_EXISTS ) ) {
    test1 = task==GET_YIELD_RULE&&plasti_type==-NONE;
    test2 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_TENSION;
    test3 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_TENSION;
    if ( test1 || test2 || test3 ) {
      if ( swit ) pri( "check plasti_tension" );
      sig_yield = plasti_data[0];
      if ( swit ) pri( "sig_yield", sig_yield );
      matrix_eigenvalues( sig, sig_princ );
      if ( swit ) pri( "sig_princ", sig_princ, MDIM );
      tmp_sig_eq = 0.;
      for ( idim=0; idim<MDIM; idim++ ) {
        if ( sig_princ[idim]>0. ) {
          tmp_sig_eq += sig_princ[idim] * sig_princ[idim];
        }
      }
      tmp_sig_eq = sqrt( scalar_dabs(tmp_sig_eq) ); 
      f_yield = tmp_sig_eq - sig_yield;
      f_flow = tmp_sig_eq - sig_yield;
      if      ( task==GET_YIELD_RULE ) {
        if (  f_yield>f ) {
          f = f_yield;
          tmp_plasti_type = GROUP_MATERI_PLASTI_TENSION;
        }
        if ( swit ) pri( "f_yield", f_yield );
      }
      else {
        assert( task==GET_FLOW_RULE );
        f = f_flow;
        if ( swit ) pri( "f_flow", f_flow );
      }
    }
  }
  if ( get_group_data( GROUP_MATERI_PLASTI_VONMISES, gr, element, new_unknowns,
      plasti_data, ldum, GET_IF_EXISTS ) ) {
    test1 = task==GET_YIELD_RULE&&plasti_type==-NONE;
    test2 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_VONMISES;
    test3 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_VONMISES;
    if ( test1 || test2 || test3 ) {
      if ( swit ) pri( "plasti_vonmises" );
      if ( get_group_data( GROUP_MATERI_PLASTI_VONMISES_NADAI, gr, element, 
        new_unknowns, work, ldum, GET_IF_EXISTS ) ) {
        C = work[0];
        kappa_0 = work[1];
        n = work[2];
        sig_yield = plasti_data[0] + C * scalar_power(kappa_0+new_unknowns[kap_indx],n);
      }
      else
        sig_yield = plasti_data[0];
      f_yield = sig_eq*sqrt(3.) - sig_yield;
      f_flow = f_yield;
      if      ( task==GET_YIELD_RULE ) {
        if ( f_yield>f ) {
          f = f_yield;
          tmp_plasti_type = GROUP_MATERI_PLASTI_VONMISES;
        }
        if ( swit ) pri( "f_yield", f_yield );
      }
      else {
        assert( task==GET_FLOW_RULE );
        f = f_flow;
        if ( swit ) pri( "f_flow", f_flow );
      }
    }
  }
  test1 = task==GET_YIELD_RULE&&plasti_type==-NONE;
  test2 = task==GET_YIELD_RULE&&plasti_type==GROUP_MATERI_PLASTI_USER;
  test3 = task==GET_FLOW_RULE&&plasti_type==GROUP_MATERI_PLASTI_USER;
  if ( test1 || test2 || test3 ) {
    if ( swit ) pri( "check plasti_user" );
    if      ( task==GET_YIELD_RULE ) {
      f_yield = NO_YIELD_F;
      if ( group_materi_plasti_user==-YES )
        user_plasti( 1, user_data, new_unknowns, old_hisv, new_hisv, sig, f_yield );
      if ( f_yield>f ) {
        f = f_yield;
        tmp_plasti_type = GROUP_MATERI_PLASTI_USER;
      }
      if ( swit ) pri( "f_yield", f_yield );
    }
    else {
      assert( task==GET_FLOW_RULE );
      if ( group_materi_plasti_user==-YES )
        user_plasti( 2, user_data, new_unknowns, old_hisv, new_hisv, sig, f_flow );
      f = f_flow;
      if ( swit ) pri( "f_flow", f_flow );
    }
  }
  new_f = f;
  if ( materi_plasti_f && derivatives ) {
    for ( idim=0; idim<ndim; idim++ ) {
      tmp = new_grad_new_unknowns[idim*nuknwn+f_indx+idim*nder];
      f += options_nonlocal * tmp;
    }
  }
  if ( task==GET_YIELD_RULE && plasti_type==-NONE ) plasti_type = tmp_plasti_type;
  if ( task==GET_YIELD_RULE || task==GET_FLOW_RULE ) {
    if ( swit ) {
      pri( "new_f", new_f );
      pri( "f", f );
      pri( "plasti_type", plasti_type );
    }
  }

      // use central differences to determine flow direction
  signorm = array_size( sig, MDIM*MDIM );
  if ( task==GET_FLOW_RULE_GRAD ) {
    array_move( sig, sig_tmp, MDIM*MDIM );
    for ( i=0; i<MDIM*MDIM; i++ ) {
      sig_tmp[i] += EPS_SIG * signorm;
      plasti_rule( element, gr, plasti_on_boundary, 
        user_data, new_unknowns, new_grad_new_unknowns,
        old_hisv, new_hisv, old_epp, inc_epp, inc_ept,
        GET_FLOW_RULE, plasti_type, sig_tmp, f_flow_right, rdum, ddum );
      sig_tmp[i] -= 2. * EPS_SIG * signorm;
      plasti_rule( element, gr, plasti_on_boundary, 
        user_data, new_unknowns, new_grad_new_unknowns,
        old_hisv, new_hisv, old_epp, inc_epp, inc_ept,
        GET_FLOW_RULE, plasti_type, sig_tmp, f_flow_left, rdum, ddum );
      dir[i] = (f_flow_right-f_flow_left)/(2.*EPS_SIG*signorm);
      sig_tmp[i] += EPS_SIG * signorm;
    }
  }

  if ( swit ) pri( "Out PLASTI_RULE" );
}
