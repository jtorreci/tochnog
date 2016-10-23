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

  if ( hypo_wolfersdorff || hypo_lowangles ) {

#if !HYPO_USE
    pri( "Error: HYPO_USE is not set to 1 in tnhypo.h" );
    pri( "Look in tochnog/src/makefile how to compile." );
    exit(TN_EXIT_STATUS);
#endif

    if ( formulation==TOTAL ) {
      pri( "Error: hypoplasticity not available for this group_materi_memory.");
      exit(TN_EXIT_STATUS);
    }

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

  }
  
}
