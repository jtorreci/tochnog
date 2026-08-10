/*
    Copyright (C) 2000  Dennis Roddeman
    FEAT, Finite Element Application Technology
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

#define LENGTH_LOWANGLES 10
#define LENGTH_WOLFERSDORFF 8
#define MAX_DATA_LENGTH 10
#define LENGTH_INTERGRANULARSTRAIN 5

extern "C" 
  int hypo_( double *stress, double *Mmat, double *new_hisv,
    double *inc_ept, double *time, double *dtime, long int *nhis, 
    double *data, long int *ndata, double *cohesion,
    double *epi_R, double *epi_mr, double *epi_mt, 
    double *epi_betar, double *epi_chi,
    double *old_epi, double *new_epi,
    long int *use_pres, long int *use_epi, long int *hypo_type,
    double *softvar_nonloc, double *softvar_loc,
    int *find_local_sv, int *options_nonlocal );

extern "C" 
  void masin_umat( double *stress, double *statev, double *ddsdde,
    double *dstran, double dtime, double *props, int nprops, int testing,
    int *error );

extern "C" 
  void masin_visco_umat( double *stress, double *statev, double *ddsdde,
    double *dstran, double dtime, double *props, int nprops, int testing,
    int *error );

extern "C" 
  void masin_niemunis_visco_umat( double *stress, double *statev,
    double *ddsdde, double *dstran, double dtime, double *props, int nprops,
    int testing, int *error );

extern "C" 
  void sanisand_umat( double *stress, double *statev, double *ddsdde,
    double *dstran, double dtime, double *props, int nprops, int testing,
    int *error );

void hypoplasticity( long int element, long int gr,
  long int formulation, double old_hisv[], double new_hisv[], 
  double old_unknowns[], double new_unknowns[], 
  double inc_ept[], double old_epi[], double new_epi[], 
  double rotated_old_sig[], double new_sig[], 
  double *Chypo, double softvar_nonl, double &softvar_l )

  /* Interface routine to hypoplasticity routine. */

{

  long int i=0, j=0, k=0, l=0, ldum=0, hypo_wolfersdorff=0, hypo_lowangles=0,
    length_lowangles=0, length_wolfersdorff=0, length_intergranularstrain=0, 
    pressure_dependent_void_ratio=0, idum[1],
    ndata[1], nhis[1], use_epi[1], use_pres[1], hypo_type[1];
  int find_local_sv[1], options_nonlocal[1];
  double ddum[1], cohesion[1], epi_R[1], epi_mr[1], 
    epi_mt[1], epi_betar[1], epi_chi[1], time[1], dtime[1],
    softvar_nonloc[1], softvar_loc[1],
    stress[MDIM*MDIM], data[MAX_DATA_LENGTH], 
    hypo_intergranularstrain[LENGTH_INTERGRANULARSTRAIN],
    Mmat[MDIM*MDIM*MDIM*MDIM];

    if(db_active_index( GROUP_MATERI_PLASTI_HYPO_WOLFERSDORFF, gr, VERSION_NORMAL )) hypo_wolfersdorff=1; 
    else if(db_active_index( GROUP_MATERI_PLASTI_HYPO_LOWANGLES, gr, VERSION_NORMAL )) hypo_lowangles=1;

  if ( hypo_wolfersdorff || hypo_lowangles || 
       db_active_index( GROUP_MATERI_PLASTI_HYPO_MASIN, gr, VERSION_NORMAL ) ) {

#if !HYPO_USE
    pri( "Error: HYPO_USE is not set to 1 in tnhypo.h" );
    pri( "Look in tochnog/src/makefile how to compile." );
    exit(TN_EXIT_STATUS);
#endif

    if ( formulation==TOTAL ) {
      pri( "Error: hypoplasticity not available for this group_materi_memory.");
      exit(TN_EXIT_STATUS);
    }

    if ( hypo_wolfersdorff || hypo_lowangles ) {

    length_lowangles = LENGTH_LOWANGLES;
    length_wolfersdorff = LENGTH_WOLFERSDORFF;
    length_intergranularstrain = LENGTH_INTERGRANULARSTRAIN;

    if( hypo_wolfersdorff) {
      db( GROUP_MATERI_PLASTI_HYPO_WOLFERSDORFF, gr, idum, data, length_wolfersdorff, 
        VERSION_NORMAL, GET_AND_CHECK );
      ndata[0] = length_wolfersdorff;
      hypo_type[0] = 0;
    }
    else if( hypo_lowangles) {
      db( GROUP_MATERI_PLASTI_HYPO_LOWANGLES, gr, idum, data, length_lowangles, 
        VERSION_NORMAL, GET_AND_CHECK );
      ndata[0] = length_lowangles;
      hypo_type[0] = 1;
    }
    if ( materi_history_variables<4 ) {
      pri( "Error: materi_history_variables should be 4 for GROUP_MATERI_PLASTI_HYPO_WOLFERSDORFF" );
      exit(TN_EXIT_STATUS);
    }

    use_epi[0] = 0;
    epi_R[0] = epi_mr[0] = epi_mt[0] = epi_betar[0] = epi_chi[0] = 0.;
    if ( materi_strain_intergranular ) {
      db( GROUP_MATERI_PLASTI_HYPO_INTERGRANULARSTRAIN, gr, idum,
        hypo_intergranularstrain, length_intergranularstrain, 
        VERSION_NORMAL, GET_AND_CHECK );
      use_epi[0] = 1;
      epi_R[0] = hypo_intergranularstrain[0];
      epi_mr[0] = hypo_intergranularstrain[1];
      epi_mt[0] = hypo_intergranularstrain[2];
      epi_betar[0] = hypo_intergranularstrain[3];
      epi_chi[0] = hypo_intergranularstrain[4];
      if ( epi_R[0]<=0. ) {
        pri( "Error: R in GROUP_MATERI_PLASTI_HYPO_INTERGRANULARSTRAIN should be positive." );
        exit(TN_EXIT_STATUS);
      }
    }

    use_pres[0] = 0;
    if ( db_active_index( GROUP_MATERI_PLASTI_HYPO_PRESSUREDEPENDENTVOIDRATIO, 
         gr, VERSION_NORMAL ) ) {
      db( GROUP_MATERI_PLASTI_HYPO_PRESSUREDEPENDENTVOIDRATIO, gr,
        &pressure_dependent_void_ratio, ddum, ldum, 
        VERSION_NORMAL, GET );
      if ( pressure_dependent_void_ratio==-YES ) use_pres[0] = 1;
    }

    cohesion[0] = 0.;
    db( GROUP_MATERI_PLASTI_HYPO_COHESION, gr, idum, cohesion, ldum, 
      VERSION_NORMAL, GET_IF_EXISTS );

      // data

    db( DTIME, 0, idum, dtime, ldum, VERSION_NEW, GET );
    db( TIME_CURRENT, 0, idum, time, ldum, VERSION_NORMAL, GET );

    nhis[0] = materi_history_variables;
    array_move( old_hisv, new_hisv, nhis[0] );
    array_move( rotated_old_sig, stress, MDIM*MDIM );
    array_set( &Mmat[0], 0., MDIM*MDIM*MDIM*MDIM );

    softvar_nonloc[0]=softvar_nonl;
    find_local_sv[0]=0;
    if(find_local_softvar) {
    	find_local_sv[0]=1;
        long int length_nei=1+npointmax*ndim+npointmax+2;
        double nonloc_info[length_nei-2];
        array_set(nonloc_info, 0., length_nei);
        db( NONLOCAL_ELEMENT_INFO, element, idum, nonloc_info, length_nei, VERSION_NORMAL, GET );		
	nonloc_info[1+npointmax*ndim+npointmax]=1.;
        db( NONLOCAL_ELEMENT_INFO, element, idum, nonloc_info, length_nei, VERSION_NORMAL, PUT );		
    }
    options_nonlocal[0]=0;	
    if (scalar_dabs(options_nonlocal_softvar)>TINY) options_nonlocal[0]=1;

      // stress contribution by hypoplasticity
      
#if HYPO_USE
    hypo_( stress, Mmat, new_hisv, inc_ept, time, dtime,
      nhis, data, ndata, cohesion, epi_R, epi_mr, epi_mt, 
      epi_betar, epi_chi, old_epi, new_epi, use_pres, use_epi, hypo_type,
      softvar_nonloc, softvar_loc, find_local_sv, options_nonlocal);
#endif
    if(find_local_softvar) softvar_l=softvar_loc[0];

      // from fortran to c
    for ( i=0; i<MDIM; i++ ) {
       for ( j=0; j<MDIM; j++ ) {
          for ( k=0; k<MDIM; k++ ) {
             for ( l=0; l<MDIM; l++ ) {
                Chypo[l*MDIM*MDIM*MDIM+k*MDIM*MDIM+j*MDIM+i] = 
                   Mmat[i*MDIM*MDIM*MDIM+j*MDIM*MDIM+k*MDIM+l];
             }
          }
       }
    }
    array_add( new_sig, stress, new_sig, MDIM*MDIM );
    array_subtract( new_sig, rotated_old_sig, new_sig, MDIM*MDIM );

    }   /* end wolfersdorff / lowangles block */

  }   /* end hypoplasticity dispatch */

  if ( db_active_index( GROUP_MATERI_PLASTI_HYPO_MASIN, gr, VERSION_NORMAL ) ||
       db_active_index( GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY, gr, VERSION_NORMAL ) ) {

      // Masin clay hypoplasticity (masin.c, port of umat_hcea.for)
      // ------------------------------------------------------------------
      // Conventions:
      //   tochnog stores 3x3 tensors row-major (sig[i*MDIM+j]); masin uses
      //   Voigt6 [11,22,33,12,13,23]. Compression negative in both.
      //   tochnog history hisv[] (materi_history_variables, >= 8):
      //     hisv[0..5] = intergranular strain delta (continuum -> Voigt)
      //     hisv[6]    = void ratio e
      //     hisv[7]    = sensitivity s
      //   masin statev[16] maps: [0..5]=delta, 6=e, 12=dtsub, 13=sensitivity.
      // ------------------------------------------------------------------
    long int mhis[8];
    double mstress[6], mdstran[6], mstatev[16], mddsdde[36], mprops[29];
    double ocr=0., e0=0., mdt=0.;
    int merror=0, mtesting=0, i2, j2;
    long int ocr_apply=-NO, hypo_masin_clay=0, hypo_masin_visco=0, diri=0;
    long int hypo_masin_visco_jm=0;

    if ( db_active_index( GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY, gr, VERSION_NORMAL ) )
      hypo_masin_clay = 1;
    if ( db_active_index( GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_VISCO, gr, VERSION_NORMAL ) )
      hypo_masin_visco = 1;
    if ( db_active_index( GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_VISCO_JM, gr, VERSION_NORMAL ) )
      hypo_masin_visco_jm = 1;

    if ( materi_history_variables<8 ) {
      pri( "Error: materi_history_variables should be at least 8 for GROUP_MATERI_PLASTI_HYPO_MASIN(_CLAY)." );
      pri( "   hisv[0..5] = intergranular strain, hisv[6] = void ratio e, hisv[7] = sensitivity." );
      exit(TN_EXIT_STATUS);
    }

    if ( formulation==TOTAL ) {
      pri( "Error: hypoplasticity not available for this group_materi_memory.");
      exit(TN_EXIT_STATUS);
    }

      // material parameters (29 props, raw layout of umat_hcea.for)
      // group_materi_plasti_hypo_masin        = phi_c lambda* kappa* N r
      // group_materi_plasti_hypo_masin_clay   = phi_c lambda* kappa* N nu_pp
      //   -> props[0]=phi_c, props[2]=lambda*, props[3]=kappa*,
      //      props[4]=N, props[5]=r/nu_pp
      //   props[1]=p_t is kept 0 (no cohesion shift)
    for ( i=0; i<29; i++ ) mprops[i] = 0.;
    length_wolfersdorff = 5;
    {
      double mpar[5];
      if ( hypo_masin_clay )
        db( GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY, gr, idum, mpar, length_wolfersdorff, VERSION_NORMAL, GET_AND_CHECK );
      else
        db( GROUP_MATERI_PLASTI_HYPO_MASIN, gr, idum, mpar, length_wolfersdorff, VERSION_NORMAL, GET_AND_CHECK );
      mprops[0] = mpar[0];   // phi_c [deg]
      mprops[2] = mpar[1];   // lambda*
      mprops[3] = mpar[2];   // kappa*
      mprops[4] = mpar[3];   // N
      mprops[5] = mpar[4];   // r / nu_pp
    }
    mprops[1] = 0.;                       // p_t (shift due to cohesion)
    mprops[6] = 1.;                       // alpha_G (isotropic by default)
    mprops[9] = 1.;                       // s_f (1 => no structure effect)
    mprops[13] = 0.;                      // A_g (0 => intergranular strain off)
    mprops[17] = 3.;                      // vertical direction (z in 3D)
    mprops[22] = 0.;                      // ay (0 => 0.30 default in kernel)
    mprops[23] = 0.;                      // oc (0 => 2.0 default in kernel)

    if ( hypo_masin_clay ) {
        // advanced parameters: alpha_G alpha_f ay oc
      if ( db_active_index( GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_ADVANCED_PARAMETERS, gr, VERSION_NORMAL ) ) {
        length_intergranularstrain = 4;
        db( GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_ADVANCED_PARAMETERS, gr, idum, &mprops[6], length_intergranularstrain, VERSION_NORMAL, GET_AND_CHECK );
          // mprops[6]=alpha_G, [20]=alpha_f, [22]=ay, [23]=oc
      }
        // direction diri: 0=1D(x), 1=2D(y), 2=3D(z) -> props[17]=diri+1
      if ( db_active_index( GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_AVANCED_DIRECTION, gr, VERSION_NORMAL ) ) {
        db( GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_AVANCED_DIRECTION, gr, &diri, ddum, ldum, VERSION_NORMAL, GET_AND_CHECK );
        if ( diri>=0 && diri<=2 ) mprops[17] = diri + 1.;
      }
        // structure: k A s_f
      if ( db_active_index( GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_STRUCTURE, gr, VERSION_NORMAL ) ) {
        length_intergranularstrain = 3;
        db( GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_STRUCTURE, gr, idum, &mprops[7], length_intergranularstrain, VERSION_NORMAL, GET_AND_CHECK );
      }
    }
    else {
        // basic model: optional structure
      if ( db_active_index( GROUP_MATERI_PLASTI_HYPO_MASIN_STRUCTURE, gr, VERSION_NORMAL ) ) {
        length_intergranularstrain = 3;
        db( GROUP_MATERI_PLASTI_HYPO_MASIN_STRUCTURE, gr, idum, &mprops[7], length_intergranularstrain, VERSION_NORMAL, GET_AND_CHECK );
          // mprops[7,8,9] = k, A, s_f
      }
    }
      // intergranular strain masin clay: R Ag ng mrat beta_r chi [theta]
      //   -> props[10]=R, [13]=A_g(G0), [14]=n_g, [15]=m_rat,
      //      [11]=beta_r, [12]=chi. theta has no direct slot (kernel uses chi).
    if ( db_active_index( GROUP_MATERI_PLASTI_HYPO_STRAIN_INTERGRANULAR_MASIN_CLAY, gr, VERSION_NORMAL ) ) {
      double mgr[7];
      length_wolfersdorff = 7;
      db( GROUP_MATERI_PLASTI_HYPO_STRAIN_INTERGRANULAR_MASIN_CLAY, gr, idum, mgr, length_wolfersdorff, VERSION_NORMAL, GET_AND_CHECK );
      mprops[10] = mgr[0];   // R
      mprops[13] = mgr[1];   // A_g
      mprops[14] = mgr[2];   // n_g
      mprops[15] = mgr[3];   // m_rat
      mprops[11] = mgr[4];   // beta_r
      mprops[12] = mgr[5];   // chi
    }
      // initial void ratio / OCR: props[21] = e0, or OCR+10 if > 10
    e0 = new_hisv[6];
    if ( e0>0.001 ) mprops[21] = e0;
    if ( hypo_masin_clay ) {
      db( GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_OCR, gr, idum, &ocr, ldum, VERSION_NORMAL, GET_IF_EXISTS );
      db( CONTROL_MATERI_PLASTI_HYPO_MASIN_CLAY_OCR_APPLY, gr, &ocr_apply, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    }
    else {
      db( GROUP_MATERI_PLASTI_HYPO_MASIN_OCR, gr, idum, &ocr, ldum, VERSION_NORMAL, GET_IF_EXISTS );
      db( CONTROL_MATERI_PLASTI_HYPO_MASIN_OCR_APPLY, gr, &ocr_apply, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    }
    if ( ocr_apply==-YES && ocr>0. ) mprops[21] = ocr + 10.;

      // visco parameters:
      //   _clay_visco (Dr Iv)     -> props[25]=Dr, props[26]=Iv  (Niemunis law)
      //   _clay_visco_jm (Dref)   -> props[21]=ocparam, [22]=beta_deg, [23]=ksi,
      //                              [24]=gama_deg, [25]=Dref (Jerman-Masin)
      //   e0/OCR goes to props[27] for the visco kernels.
    if ( hypo_masin_visco ) {
      double mvisco[2];
      length_wolfersdorff = 2;
      db( GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_VISCO, gr, idum, mvisco,
        length_wolfersdorff, VERSION_NORMAL, GET_AND_CHECK );
      mprops[25] = mvisco[0];   // Dr
      mprops[26] = mvisco[1];   // Iv
      mprops[27] = e0;
      if ( ocr_apply==-YES && ocr>0. ) mprops[27] = ocr + 10.;
    }
    if ( hypo_masin_visco_jm ) {
      double mvisco[5];
      length_wolfersdorff = 5;
      db( GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_VISCO_JM, gr, idum, mvisco,
        length_wolfersdorff, VERSION_NORMAL, GET_AND_CHECK );
      mprops[21] = mvisco[0];   // ocparam
      mprops[22] = mvisco[1];   // beta_deg
      mprops[23] = mvisco[2];   // ksi
      mprops[24] = mvisco[3];   // gama_deg
      mprops[25] = mvisco[4];   // Dref
      mprops[27] = e0;
      if ( ocr_apply==-YES && ocr>0. ) mprops[27] = ocr + 10.;
    }

      // strain increment: 3x3 (row-major) -> Voigt6
    mdstran[0] = inc_ept[0*MDIM+0];
    mdstran[1] = inc_ept[1*MDIM+1];
    mdstran[2] = inc_ept[2*MDIM+2];
    mdstran[3] = inc_ept[0*MDIM+1];
    mdstran[4] = inc_ept[0*MDIM+2];
    mdstran[5] = inc_ept[1*MDIM+2];

      // stress: 3x3 (row-major) -> Voigt6
    mstress[0] = rotated_old_sig[0*MDIM+0];
    mstress[1] = rotated_old_sig[1*MDIM+1];
    mstress[2] = rotated_old_sig[2*MDIM+2];
    mstress[3] = rotated_old_sig[0*MDIM+1];
    mstress[4] = rotated_old_sig[0*MDIM+2];
    mstress[5] = rotated_old_sig[1*MDIM+2];

      // history: hisv -> statev (layout of umat_hcea.for)
      // NOTE: the Fortran reference defines move_asv_hcea (which negates the
      // intergranular strain) but NEVER calls it in the integration path:
      // iniy_hcea copies asv -> y(6+i) directly. So no sign flip here.
    for ( i=0; i<16; i++ ) mstatev[i] = 0.;
    for ( i=0; i<6; i++ ) mstatev[i] = new_hisv[i];   // intergranular strain
    mstatev[6]  = new_hisv[6];   // void ratio
    mstatev[7]  = 0.;            // excess pore pressure
    mstatev[12] = 0.;            // dtsub (suggested substep, recomputed)
    mstatev[13] = new_hisv[7];   // sensitivity
    mstatev[15] = 0.;

    db( DTIME, 0, idum, &mdt, ldum, VERSION_NEW, GET );

      // stress contribution by Masin hypoplasticity
    if ( hypo_masin_visco )
      masin_niemunis_visco_umat( mstress, mstatev, mddsdde, mdstran, mdt,
        mprops, 29, mtesting, &merror );
    else if ( hypo_masin_visco_jm )
      masin_visco_umat( mstress, mstatev, mddsdde, mdstran, mdt, mprops, 29,
        mtesting, &merror );
    else
      masin_umat( mstress, mstatev, mddsdde, mdstran, mdt, mprops, 29,
        mtesting, &merror );
    if ( merror==10 ) {
      pri( "Error: severe error in Masin hypoplasticity." );
      exit(TN_EXIT_STATUS);
    }

      // statev -> hisv
    for ( i=0; i<6; i++ ) new_hisv[i] = mstatev[i];
    new_hisv[6] = mstatev[6];
    new_hisv[7] = mstatev[13];

      // Voigt6 -> 3x3 (row-major), symmetric
    for ( i2=0; i2<3; i2++ )
      for ( j2=0; j2<3; j2++ ) {
        if      ( i2==0 && j2==0 ) stress[0] = mstress[0];
        else if ( i2==1 && j2==1 ) stress[4] = mstress[1];
        else if ( i2==2 && j2==2 ) stress[8] = mstress[2];
        else if ( (i2==0&&j2==1)||(i2==1&&j2==0) ) stress[i2*MDIM+j2] = mstress[3];
        else if ( (i2==0&&j2==2)||(i2==2&&j2==0) ) stress[i2*MDIM+j2] = mstress[4];
        else if ( (i2==1&&j2==2)||(i2==2&&j2==1) ) stress[i2*MDIM+j2] = mstress[5];
      }

      // tangent ddsdde (Voigt6) -> Chypo (3x3x3x3, row-major)
    {
      int vi[3][3];
      vi[0][0]=0; vi[1][1]=1; vi[2][2]=2;
      vi[0][1]=3; vi[1][0]=3; vi[0][2]=4; vi[2][0]=4; vi[1][2]=5; vi[2][1]=5;
      for ( i=0; i<MDIM; i++ )
        for ( j=0; j<MDIM; j++ )
          for ( k=0; k<MDIM; k++ )
            for ( l=0; l<MDIM; l++ )
              Chypo[i*MDIM*MDIM*MDIM + j*MDIM*MDIM + k*MDIM + l] =
                mddsdde[ vi[i][j]*6 + vi[k][l] ];
    }

    array_add( new_sig, stress, new_sig, MDIM*MDIM );
    array_subtract( new_sig, rotated_old_sig, new_sig, MDIM*MDIM );



    (void)mhis;

    (void)mhis;
  }

  if ( db_active_index( GROUP_MATERI_PLASTI_SANISAND, gr, VERSION_NORMAL ) ) {

      // SANISAND (Dafalias & Manzari 2004), port of the reference UMAT.
      // ------------------------------------------------------------------
      // Conventions: SANISAND uses SOIL mechanics (compression positive),
      // tochnog uses Abaqus-like (tension positive). The kernel does its own
      // sign conversion (move_sig/move_eps negate), so we pass the raw
      // tensors. History: statev[0..35] is stored in hisv (36 slots).
      //   hisv[0..5]  = back stress alpha
      //   hisv[6]     = void ratio e
      //   hisv[7..12] = fabric tensor z
      //   hisv[14..19]= alpha at stress reversal
      //   hisv[28..33]= pore, p', q, cos3t, dtsub, nfev
      //   (the kernel reads/writes all 36 through the statev array)
      // ------------------------------------------------------------------
    double sstress[6], sstatev[36], sddsdde[36], sdstran[6], sprops[19];
    double mdt2=0.;
    int serror=0, stesting=0, i2b, j2b;

    if ( materi_history_variables<36 ) {
      pri( "Error: materi_history_variables should be at least 36 for GROUP_MATERI_PLASTI_SANISAND." );
      pri( "   hisv[0..5]=alpha, hisv[6]=e, hisv[7..12]=z, hisv[14..19]=alpha_sr." );
      exit(TN_EXIT_STATUS);
    }
    if ( formulation==TOTAL ) {
      pri( "Error: hypoplasticity not available for this group_materi_memory.");
      exit(TN_EXIT_STATUS);
    }

      // material parameters (19)
    for ( i=0; i<19; i++ ) sprops[i] = 0.;
    length_wolfersdorff = 19;
    db( GROUP_MATERI_PLASTI_SANISAND, gr, idum, sprops, length_wolfersdorff,
      VERSION_NORMAL, GET_AND_CHECK );

      // state: hisv -> statev (36)
    for ( i=0; i<36; i++ ) sstatev[i] = new_hisv[i];

      // stress/strain: pass raw tensors (kernel converts sign)
    sstress[0] = rotated_old_sig[0*MDIM+0];
    sstress[1] = rotated_old_sig[1*MDIM+1];
    sstress[2] = rotated_old_sig[2*MDIM+2];
    sstress[3] = rotated_old_sig[0*MDIM+1];
    sstress[4] = rotated_old_sig[0*MDIM+2];
    sstress[5] = rotated_old_sig[1*MDIM+2];
    sdstran[0] = inc_ept[0*MDIM+0];
    sdstran[1] = inc_ept[1*MDIM+1];
    sdstran[2] = inc_ept[2*MDIM+2];
    sdstran[3] = inc_ept[0*MDIM+1];
    sdstran[4] = inc_ept[0*MDIM+2];
    sdstran[5] = inc_ept[1*MDIM+2];

    db( DTIME, 0, idum, &mdt2, ldum, VERSION_NEW, GET );

    sanisand_umat( sstress, sstatev, sddsdde, sdstran, mdt2, sprops, 19,
      stesting, &serror );
    if ( serror==10 ) {
      pri( "Error: severe error in SANISAND." );
      exit(TN_EXIT_STATUS);
    }

      // statev -> hisv
    for ( i=0; i<36; i++ ) new_hisv[i] = sstatev[i];

      // stress: kernel returns Abaqus-convention (tension positive) stress
    stress[0*MDIM+0] = sstress[0];
    stress[1*MDIM+1] = sstress[1];
    stress[2*MDIM+2] = sstress[2];
    stress[0*MDIM+1] = sstress[3];  stress[1*MDIM+0] = sstress[3];
    stress[0*MDIM+2] = sstress[4];  stress[2*MDIM+0] = sstress[4];
    stress[1*MDIM+2] = sstress[5];  stress[2*MDIM+1] = sstress[5];

      // tangent: Voigt6 -> Chypo
    {
      int vi[3][3];
      vi[0][0]=0; vi[1][1]=1; vi[2][2]=2;
      vi[0][1]=3; vi[1][0]=3; vi[0][2]=4; vi[2][0]=4; vi[1][2]=5; vi[2][1]=5;
      for ( i=0; i<MDIM; i++ )
        for ( j=0; j<MDIM; j++ )
          for ( k=0; k<MDIM; k++ )
            for ( l=0; l<MDIM; l++ )
              Chypo[i*MDIM*MDIM*MDIM + j*MDIM*MDIM + k*MDIM + l] =
                sddsdde[ vi[i][j]*6 + vi[k][l] ];
    }

    array_add( new_sig, stress, new_sig, MDIM*MDIM );
    array_subtract( new_sig, rotated_old_sig, new_sig, MDIM*MDIM );

    (void)i2b; (void)j2b;
  }

}
