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
    double *drplde, double *drpldt, double *stran, double *dstran, 
    double *time, double *dtime, double *temp, double *dtemp, double *predef, 
    double *dpred, char *cmname, long int *ndi, 
    long int *nshr, long int *ntens, long int *nstatv, 
    double *props, long int *nprops, double *coords, double *drot, 
    double *pnewdt, double *celent, double *dfgrd0, double *dfgrd1, 
    long int *noel, long int *npt, long int *layer, long int *kspt, 
    long int *kstep, long int *kinc, short cmname_len);

#define MAX_ITER 1000
#define EPS_F 1.e-6
#define EPS_EPS_F 1.e-3
#define EPS_SIZE 1.e-3
#define EPS_DAMAGE 1.e-5
#define EPS_TMP 1.e-1
#define EPS_LAMBDA 1.e-3
#define EPS_VISCO 3.

void set_stress( long int element, long int gr, 
  long int plasti_on_boundary, double coord_ip[],
  double old_unknowns[], double new_unknowns[], 
  double old_grad_old_unknowns[], double new_grad_new_unknowns[], 
  double rotated_old_sig[], double new_sig[],
  double rotated_old_msig[], double new_msig[], 
  double inc_ept[], double new_ept[], 
  double old_epe[], double inc_epe[], 
  double old_epp[], double inc_epp[], 
  double old_rho[], double new_rho[],
  double old_epi[], double new_epi[],
  double old_hisv[], double new_hisv[], 
  double old_damage, double &new_damage, 
  double old_kappa, double &new_kappa, 
  double &new_f, double &new_substeps, double old_deften[], double new_deften[],
  double inc_rot[], double ddsdde[],
  double &viscosity, double &viscosity_heat_generation, double &softvar_nonl,
  double &softvar_l )

  // Solid materials.

{
  long int i=0, j=0, ind_ddsdde=0, plasti_found=0, 
    plasti_iter=0, membrane_found=0, membrane_iter=0,
    swit=0, length=0, membrane=-NO, viscoplasti=0, 
    viscoplasti_always=-NO, plasti_type=-NONE, volumetric_young_order=0,
    memory=-UPDATED, max_plasti_iter=0, total_plasti_iter=0, 
    nuser_data=0, idim=0, jdim=0, kdim=0, ldim=0, 
    formulation=INCREMENTAL, ldum=0, idum[1], task[2];
  double lambda=0., deps_size=0., lambda_new=0., lambda_previous=0., 
    tmp=0., tmp_old=0., tmp_inc=0., tmp_new=0., 
    f=0., f_previous=0, f_ref=0., materi_expansion_linear=0., 
    eta=0., pressure=0., strain_size=0., straindev_size=0., linear_ept=0.,zero_ept=0.,
	 meanstrain=0.,
    plasti_kinematic_hardening=0., young = 0., young_linear_ept=0., poisson=0., 
    compressibility=0., fac=0., kappa=0., g=0., k=0., e=0.,
    lade_1=0., lade_2=0, lade_3=0., p=0., p0=0., young0=0.,
    alpha=0., gamma=0., dtime=0., rdum=0., camclay[1], tskh[DATA_ITEM_SIZE],
    smallstrain[6],
    group_materi_elasti_lade[3], sig_dev[MDIM*MDIM], ddum[MDIM*MDIM], ddumarray[MDIM][MDIM], 
    memmat[MDIM][MDIM], elasti_transverse_isotropy[DATA_ITEM_SIZE],
    inc_temperature_strain[MDIM*MDIM], new_temperature_strain[MDIM*MDIM], 
    plasti_dir[MDIM*MDIM], young_power[3], young_polynomial[DATA_ITEM_SIZE],
    young_strainstress[DATA_ITEM_SIZE], ept_dev[MDIM*MDIM],
    plasti_visco_exponential[2], plasti_visco_power[3], 
    inc_rho[MDIM*MDIM], test_sig[MDIM*MDIM], total_inc_epp[MDIM*MDIM], 
    work_inc_epp[MDIM*MDIM], new_epe[MDIM*MDIM], new_epp[MDIM*MDIM],
    C[MDIM][MDIM][MDIM][MDIM], Cmem[MDIM][MDIM][MDIM][MDIM], 
    Cuser[MDIM][MDIM][MDIM][MDIM], Chyper[MDIM][MDIM][MDIM][MDIM],
    Chypo[MDIM][MDIM][MDIM][MDIM],
    volumetric_young_values[DATA_ITEM_SIZE],
    user_data[DATA_ITEM_SIZE], work[DATA_ITEM_SIZE], 
    old_work[MDIM*MDIM], new_work[MDIM*MDIM];

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
      new_unknowns, plasti_visco_power, ldum, GET_IF_EXISTS ) ) {
    viscoplasti = 1;
    eta = plasti_visco_power[0]; p = plasti_visco_power[1];
    f_ref = plasti_visco_power[2]; assert( f_ref!=0. );
  }

  db( DTIME, 0, idum, &dtime, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( GROUP_MATERI_PLASTI_VISCO_ALWAYS, gr, &viscoplasti_always, 
    ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( GROUP_MATERI_MEMORY, gr, &memory, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );

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

  if ( get_group_data( GROUP_MATERI_ELASTI_YOUNG, gr, element, new_unknowns, 
      &young, ldum, GET_IF_EXISTS ) ) {
    get_group_data( GROUP_MATERI_ELASTI_POISSON, gr, element,
      new_unknowns, &poisson, ldum, GET_IF_EXISTS ); 
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
    get_group_data( GROUP_MATERI_ELASTI_POISSON, gr, element,
      new_unknowns, &poisson, ldum, GET_IF_EXISTS ); 
    p0 = young_power[0];
    young0 = young_power[1];
    alpha = young_power[2];
    if ( p0<=0. ) db_error( GROUP_MATERI_ELASTI_YOUNG_POWER, gr );
    p = - ( new_sig[0] + new_sig[4] + new_sig[8] ) / 3.;
    young = young0 * scalar_power(scalar_dabs(p/p0),alpha);
    task[1] = -NO;
    C_matrix( young, poisson, elasti_transverse_isotropy, C, task );
    task[1] = membrane;
    C_matrix( young, poisson, elasti_transverse_isotropy, Cmem, task );
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

  if ( memory==-TOTAL || memory==-TOTAL_PIOLA  || memory==-TOTAL_LINEAR ) {
    if ( !materi_strain_total  ) formulation = TOTAL;
    check_unknown( "materi_velocity", YES, CHECK_USAGE_AND_ERROR );
    check_unknown( "materi_displacement", YES, CHECK_USAGE_AND_ERROR );
    check_unknown( "materi_stress", YES, CHECK_USAGE_AND_ERROR );
  }
  if ( memory==-UPDATED || memory==-UPDATED_WITHOUT_ROTATION ) {
    check_unknown( "materi_velocity", YES, CHECK_USAGE_AND_ERROR );
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

  if      ( materi_plasti_f_nonlocal && viscoplasti ) 
    max_plasti_iter = 2;
  else if ( viscoplasti )
    max_plasti_iter = 10;
  else
    max_plasti_iter = MAX_ITER;

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
  }
  else {

  /*********************************************************************************/

    // plastic iterations
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

        // hypoplasticity
      hypoplasticity( element, gr, formulation,
        old_hisv, new_hisv, old_unknowns, new_unknowns,
        inc_ept, old_epi, new_epi, rotated_old_sig, 
        new_sig, &Chypo[0][0][0][0], softvar_nonl, softvar_l );
      if ( swit ) pri( "stress after hypoplasticity", new_sig, MDIM, MDIM );

        // compressibility
      if ( compressibility!=0. ) {
        if      ( memory==-TOTAL_LINEAR ) {
          tmp_old = ( old_epe[0] + old_epe[4] + old_epe[8] ) / compressibility;
          tmp_new = ( new_epe[0] + new_epe[4] + new_epe[8] ) / compressibility;
          tmp_inc = tmp_new - tmp_old;
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

        // damage
      if ( materi_damage ) {
        damage( gr, new_epe, new_sig, old_damage, new_damage );
        array_multiply( new_sig, new_sig, 1.-new_damage, MDIM*MDIM );
        array_multiply( &memmat[0][0], &memmat[0][0], 1.-new_damage, MDIM*MDIM );
        if ( swit ) pri( "stress after damage", new_sig, MDIM, MDIM );
      }

          // membrane iterations ready?
      membrane_found = membrane_apply( element, gr, memmat, inc_ept,
        inc_epe, new_ept, new_epe, new_sig );
	//I have always membrane_found==TRUE
    }

      // test stresses for plastic yield functions
    array_move( new_sig, test_sig, MDIM*MDIM );
    if ( materi_plasti_rho ) 
      array_subtract( test_sig, new_rho, test_sig, MDIM*MDIM );
      // new lambda
    if ( plasti_iter==1 ) {
      plasti_type = -NONE; plasti_rule( element, gr, plasti_on_boundary, user_data, 
        new_unknowns, new_grad_new_unknowns, old_hisv, new_hisv,
        old_epp, inc_epp, inc_ept, GET_YIELD_RULE, 
        plasti_type, test_sig, f, new_f, ddum );

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
        if ( viscoplasti ) {
          lambda_previous = 0.; f_previous = f;
          if ( plasti_visco_exponential[0]>0. ) {
            tmp = alpha*f;
            if ( tmp>EPS_VISCO ) tmp = EPS_VISCO;
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
          new_unknowns, new_grad_new_unknowns,
          old_hisv, new_hisv, old_epp, inc_epp, inc_ept,
          GET_FLOW_RULE_GRAD, plasti_type, test_sig, rdum, rdum, plasti_dir );
        if ( !viscoplasti ) array_normalize( plasti_dir, MDIM*MDIM );
        if ( swit ) pri( "plasti_dir", plasti_dir, MDIM*MDIM );
      }
    }
    else {
      plasti_rule( element, gr, plasti_on_boundary, user_data, new_unknowns, new_grad_new_unknowns,
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
          new_unknowns, new_grad_new_unknowns, old_hisv, new_hisv,
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
          if ( tmp>EPS_VISCO ) tmp = EPS_VISCO;
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
          ddsdde[ind_ddsdde] += fac * ( Cmem[idim][jdim][kdim][ldim]
            + Cuser[idim][jdim][kdim][ldim] + Chyper[idim][jdim][kdim][ldim] 
            + Chypo[idim][jdim][kdim][ldim] );
          if ( compressibility!=0. && idim==jdim && kdim==ldim )
            ddsdde[ind_ddsdde] += 1./compressibility;
        }
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
