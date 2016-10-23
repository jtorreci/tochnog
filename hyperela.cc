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

#define DELTA 1.e-2
#define EPS_VARIATION 1.e-3

void hyperelasticity( long int gr, long int element,
  long int memory, double unknowns[], double epe[], double sig[],
  double Chyper[MDIM][MDIM][MDIM][MDIM] )

{
  long int hyper_stiffness=-YES, swit=0, ldum=0;
  double ddum[1], new_stress[MDIM*MDIM];

  swit = set_swit(element,-1,"hyperelasticity");
  array_set( sig, 0., MDIM*MDIM );

  if ( hyper_stress( gr, element, memory, unknowns, epe, new_stress ) ) {
    array_move( new_stress, sig, MDIM*MDIM );
    db( GROUP_MATERI_HYPER_STIFFNESS, gr, &hyper_stiffness, ddum, 
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
    if ( hyper_stiffness==-YES ) 
      hyper_Cmat( gr, element, memory, unknowns, epe, Chyper );
  }
}

void hyper_Cmat( long int gr, long int element, long int memory, 
  double unknowns[], double epe[], double Chyper[MDIM][MDIM][MDIM][MDIM] )

{
  long int idim=0, jdim=0, kdim=0, ldim=0, indx=0;
  double tmp=0., variation=0, work[MDIM*MDIM], stress_left[MDIM*MDIM], 
    stress_right[MDIM*MDIM];

        // central differences to determine stresses
  variation = DELTA * array_size( epe, MDIM*MDIM );
  if ( scalar_dabs(variation)<EPS_VARIATION ) variation = EPS_VARIATION;

  array_move( epe, work, MDIM*MDIM );

  for ( kdim=0; kdim<MDIM; kdim++ ) {
    for ( ldim=0; ldim<MDIM; ldim++ ) {
      indx = kdim*MDIM + ldim;
      work[indx] += variation;
      hyper_stress( gr, element, memory, unknowns, work, stress_right );
      work[indx] -= 2. * variation;
      hyper_stress( gr, element, memory, unknowns, work, stress_left );
      for ( idim=0; idim<MDIM; idim++ ) {
        for ( jdim=0; jdim<MDIM; jdim++ ) {
          tmp = ( stress_right[idim*MDIM+jdim] -
            stress_left[idim*MDIM+jdim] ) / ( 2.*variation );
          if ( idim==kdim && jdim==ldim ) {
            if ( tmp<0. ) tmp = 0.;
          }
          Chyper[idim][jdim][kdim][ldim] = tmp;
        }
      }
      work[indx] += variation;
    }
  }
}

long int hyper_stress( long int gr, long int element, long int memory, 
  double unknowns[], double epe[], double stress[] )

{
  long int idim=0, jdim=0, indx=0;
  double variation=0, W_left=0., W_right=0., C[MDIM*MDIM], work[MDIM*MDIM];

  if ( db_active_index( GROUP_MATERI_HYPER_BESSELING, gr, VERSION_NORMAL ) ||
       db_active_index( GROUP_MATERI_HYPER_BLATZ_KO, gr, VERSION_NORMAL ) ||
       db_active_index( GROUP_MATERI_HYPER_VOLUMETRIC_MURNAGHAN, gr, VERSION_NORMAL ) ||
       db_active_index( GROUP_MATERI_HYPER_MOONEY_RIVLIN, gr, VERSION_NORMAL ) ||
       db_active_index( GROUP_MATERI_HYPER_NEOHOOKEAN, gr, VERSION_NORMAL ) ||
       db_active_index( GROUP_MATERI_HYPER_REDUCEDPOLYNOMIAL, gr, VERSION_NORMAL ) ||
       db_active_index( GROUP_MATERI_HYPER_VOLUMETRIC_LINEAR, gr, VERSION_NORMAL ) ||
       db_active_index( GROUP_MATERI_HYPER_VOLUMETRIC_SIMOTAYLOR, gr, VERSION_NORMAL ) ||
       db_active_index( GROUP_MATERI_HYPER_VOLUMETRIC_OGDEN, gr, VERSION_NORMAL ) ||
       db_active_index( GROUP_MATERI_HYPER_VOLUMETRIC_POLYNOMIAL, gr, VERSION_NORMAL ) )
    ;
  else
    return 0;

  if      ( memory==-TOTAL_PIOLA ) {
    array_move( epe, C, MDIM*MDIM );
    array_multiply( C, C, 2., MDIM*MDIM );
    for ( idim=0; idim<MDIM; idim++ ) C[idim*MDIM+idim] += 1.;
  }
  else if ( memory==-TOTAL || memory==-TOTAL_LINEAR || materi_strain_elasti ) {
    array_move( epe, work, MDIM*MDIM );
    for ( idim=0; idim<MDIM; idim++ ) work[idim*MDIM+idim] += 1.;
    matrix_atb( work, work, C, MDIM, MDIM, MDIM );
  }
  else {
    return 0;
  }

        // central differences to determine stresses
  variation = DELTA * array_size( epe, MDIM*MDIM );
  if ( scalar_dabs(variation)<EPS_VARIATION ) variation = EPS_VARIATION;

  for ( idim=0; idim<MDIM; idim++ ) {
    for ( jdim=idim; jdim<MDIM; jdim++ ) {
      indx = idim*MDIM + jdim;
      C[indx] += variation;
      hyper_law( gr, element, unknowns, C, W_right );
      C[indx] -= 2. * variation;
      hyper_law( gr, element, unknowns, C, W_left );
      stress[idim*MDIM+jdim] = stress[jdim*MDIM+idim] = 
        2.*(W_right-W_left)/(2.*variation);
      C[indx] += variation;
    }
  }

  return 1;
}

void hyper_law( long int gr, long int element, double unknowns[], double C[], double &W )

{
  long int i=0, length=0, ldum=0;
  double I1=0., I2=0., I3=0., J=0., J1=0., J2=0., K=0., K_1=0., K_2=0., 
    K_i=0., alpha=0., beta=0., tmp=0., inv_C[MDIM], hyper_data[DATA_ITEM_SIZE],G=0.0;

  W = 0.;

  matrix_invariants( C, inv_C ); I1 = inv_C[0]; I2 = inv_C[1]; I3 = inv_C[2]; 
  if ( I3==0. )
    W = 0.;
  else {
    J1 = I1 / scalar_power(I3,1./3.); J2 = I2 / scalar_power(I3,2./3.);
    J = sqrt(I3);
    if ( J==0. ) {
      pri( "Error detected in hyperelasticity." );
      pri( "Probably too large element distortion. Try smaller timesteps or so." );
      exit(TN_EXIT_STATUS);
    }
    if ( get_group_data( GROUP_MATERI_HYPER_BESSELING, gr, element, unknowns,
        hyper_data, ldum, GET_IF_EXISTS ) ) {
      K_1 = hyper_data[0];
      K_2 = hyper_data[1];
      alpha = hyper_data[2];
      W += K_1*scalar_power((J1-3.),alpha) + K_2*(J2-3.);
    }
    if ( get_group_data( GROUP_MATERI_HYPER_BLATZ_KO, gr, element, unknowns,
        hyper_data, ldum, GET_IF_EXISTS ) ) {
      G = hyper_data[0];         /*  G = shear modulus */
      beta = hyper_data[1];  /*  beta = 2*nu/(1-2*nu)  */
      if ( beta==0. ) db_error( GROUP_MATERI_HYPER_BLATZ_KO, gr );
      if ( G==0. ) db_error( GROUP_MATERI_HYPER_BLATZ_KO, gr );
      W += (G*0.5)*(I1-3.0+(2.0/beta)*(scalar_power(J,-beta)-1.));
    }
    if ( get_group_data( GROUP_MATERI_HYPER_VOLUMETRIC_MURNAGHAN, gr, element, unknowns,
        hyper_data, ldum, GET_IF_EXISTS ) ) {
      K = hyper_data[0];
      beta = hyper_data[1];
      if ( beta==1. ) db_error( GROUP_MATERI_HYPER_VOLUMETRIC_OGDEN, gr );
      W += (K/beta)*((1./(beta-1.)*scalar_power(J,-beta)+1.)*J);
    }
    if ( get_group_data( GROUP_MATERI_HYPER_MOONEY_RIVLIN, gr, element, unknowns,
        hyper_data, ldum, GET_IF_EXISTS ) ) {
      K_1 = hyper_data[0];
      K_2 = hyper_data[1];
      W += K_1*(J1-3.) + K_2*(J2-3.);
    }
    if ( get_group_data( GROUP_MATERI_HYPER_NEOHOOKEAN, gr, element, unknowns,
        hyper_data, ldum, GET_IF_EXISTS ) ) {
      K_1 = hyper_data[0];
      W += K_1*(J1-3.);
    }
    if ( get_group_data( GROUP_MATERI_HYPER_REDUCEDPOLYNOMIAL, gr, element, unknowns,
        hyper_data, length, GET_IF_EXISTS ) ) {
      for ( i=1; i<=length; i++ ) {
        K_i = hyper_data[i-1];
        W += K_i*scalar_power(J1-3.,i);
      }
    }
    if ( get_group_data( GROUP_MATERI_HYPER_VOLUMETRIC_LINEAR, gr, element, unknowns,
        hyper_data, ldum, GET_IF_EXISTS ) ) {
      K = hyper_data[0];
      W += (K/2.)*scalar_square(J-1.);
    }
    if ( get_group_data( GROUP_MATERI_HYPER_VOLUMETRIC_SIMOTAYLOR, gr, element, unknowns,
        hyper_data, ldum, GET_IF_EXISTS ) ) {
      K = hyper_data[0];
      W += (K/2.)*(scalar_square(J-1.)+scalar_square(log(J)));
    }
    if ( get_group_data( GROUP_MATERI_HYPER_VOLUMETRIC_OGDEN, gr, element, unknowns,
        hyper_data, ldum, GET_IF_EXISTS ) ) {
      K = hyper_data[0];
      beta = hyper_data[1];
      if ( beta==0. ) db_error( GROUP_MATERI_HYPER_VOLUMETRIC_OGDEN, gr );
      W += (K/beta)*((1./beta)*(scalar_power(J,-beta)-1.)+log(J));
    }
    if ( get_group_data( GROUP_MATERI_HYPER_VOLUMETRIC_POLYNOMIAL, gr, element, unknowns,
        hyper_data, length, GET_IF_EXISTS ) ) {
      tmp = 1.;
      for ( i=0; i<length; i++ ) {
        tmp *= scalar_square(J-1.);
        K_i = hyper_data[i];
        W += (K_i/2.)*tmp;
      }
    }
  }

  if ( W<0. ) W = 0.;

}
