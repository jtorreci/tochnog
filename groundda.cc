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

void groundflow_data( long int element, long int gr, long int nodes[],
  double old_unknowns[], double new_unknowns[], double coord_ip[],
  double pe[], double &C, double /*h*/[], long int nnol )

{
  long int ldum=0, idum[1], inol=0, idim=0, vert=0, inod=0,
    icontrol=0, nonsaturated_apply=-YES, permeability_length=0;
  double ddum[1], pvs[5], vg[5], sigv=0., tmp=0., sig=0., pres=0., vertmax=0.,
    *node_dof=NULL, por=0., head=0., S=0., Se=0., krel=0., dS=0., dens=0.,
    gravity=0., eps_permeability=0.;
  double force_gravity[MDIM];


  get_group_data( GROUP_GROUNDFLOW_PERMEABILITY, gr, element, new_unknowns,
    pe, permeability_length, GET_IF_EXISTS );
  get_group_data( GROUP_GROUNDFLOW_CAPACITY, gr, element, new_unknowns, 
    &C, ldum, GET_IF_EXISTS );
  // group_groundflow_permeability (manual Professional 6.618): one value
  // is used in each space direction (isotropic); otherwise one value per
  // direction (pex pey pez).
  if ( permeability_length==1 && ndim>1 ) {
    for ( idim=1; idim<ndim; idim++ ) pe[idim] = pe[0];
  }

  // group_groundflow_permeability_vertical_stress: kp = a / (sigv/sig0)^b,
  // clamped to [minimum, maximum]. sigv = vertical EFFECTIVE stress.
  // Only used in 2D/3D; falls back to group_groundflow_permeability when the
  // vertical stress is extremely small (avoid division by zero).
  if ( db_active_index( GROUP_GROUNDFLOW_PERMEABILITY_VERTICAL_STRESS, gr, VERSION_NORMAL ) ) {
    ldum = 5;
    get_group_data( GROUP_GROUNDFLOW_PERMEABILITY_VERTICAL_STRESS, gr, element,
      new_unknowns, pvs, ldum, GET_AND_CHECK );
    // pvs[0]=a, pvs[1]=b, pvs[2]=sig0, pvs[3]=minimum, pvs[4]=maximum

    force_gravity_calculate( force_gravity );

    // vertical direction = component of gravity with largest magnitude
    vert = 0; vertmax = fabs(force_gravity[0]);
    for ( idim=1; idim<ndim; idim++ ) {
      if ( fabs(force_gravity[idim]) > vertmax ) { vertmax = fabs(force_gravity[idim]); vert = idim; }
    }

    // average the effective vertical stress over the element nodes
    sigv = 0.;
    for ( inol=0; inol<nnol; inol++ ) {
      inod = nodes[inol];
      node_dof = db_dbl( NODE_DOF, inod, VERSION_NEW );
      sig = node_dof[ stres_indx + stress_indx(vert,vert)*nder ];
      pres = node_dof[ pres_indx ];
      sigv += ( sig - pres );   // effective vertical stress
    }
    sigv /= (double)nnol;
    if ( fabs(sigv) > 1.e-10 ) {
      tmp = pvs[0] / pow( fabs(sigv)/pvs[2], pvs[1] );
      if ( tmp < pvs[3] ) tmp = pvs[3];
      if ( tmp > pvs[4] ) tmp = pvs[4];
      for ( idim=0; idim<ndim; idim++ ) pe[idim] = tmp;
    }
    // else: keep pe from group_groundflow_permeability (tiny vertical stress)
  }

  // van Genuchten non-saturated ground water flow:
  //   S(phi_p) = Sres + (Ssat-Sres) * (1 + (ga*|phi_p|)^gn) ^ ((1-gn)/gn)
  //   phi_p = -pres/(dens*|gravity|)  (pore-pressure head)
  //   c = csat + por * dS/dphi_p
  //   ki = krel(S) * ksat,i
  //   Se = (S-Sres)/(Ssat-Sres)
  //   krel = Se^gl * (1 - (1 - Se^(gn/(gn-1)))^((gn-1)/gn))^2
  // Group switch group_groundflow_nonsaturated_apply (global and per timestep)
  // gates the law: -no keeps the saturated behaviour only.
  if ( db_active_index( GROUP_GROUNDFLOW_NONSATURATED_VANGENUCHTEN, gr, VERSION_NORMAL ) ) {
    db( GROUNDFLOW_NONSATURATED_APPLY, 0, &nonsaturated_apply, ddum, ldum,
      VERSION_NORMAL, GET_IF_EXISTS );
    db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( CONTROL_GROUNDFLOW_NONSATURATED_APPLY, icontrol, &nonsaturated_apply,
      ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    if ( nonsaturated_apply==-YES ) {
      ldum = 5;
      get_group_data( GROUP_GROUNDFLOW_NONSATURATED_VANGENUCHTEN, gr, element,
        new_unknowns, vg, ldum, GET_AND_CHECK );
      // vg[0]=Sresidu, vg[1]=Ssat, vg[2]=ga, vg[3]=gl, vg[4]=gn
      ldum = 1;
      get_group_data( GROUP_GROUNDFLOW_POROSITY, gr, element, new_unknowns,
        &por, ldum, GET_IF_EXISTS );
      if ( db_active_index( GROUP_GROUNDFLOW_NONSATURATED_EPS_PERMEABILITY,
          gr, VERSION_NORMAL ) )
        get_group_data( GROUP_GROUNDFLOW_NONSATURATED_EPS_PERMEABILITY, gr, element,
          new_unknowns, &eps_permeability, ldum, GET_AND_CHECK );
      db( GROUNDFLOW_DENSITY, 0, idum, &dens, ldum, VERSION_NORMAL, GET_IF_EXISTS );
      force_gravity_calculate( force_gravity );
      gravity = 0.;
      for ( idim=0; idim<ndim; idim++ )
        gravity += fabs( force_gravity[idim] );

      // pore-pressure head at integration point
      if ( dens*gravity > 0. ) {
        head = -new_unknowns[pres_indx] / ( dens * gravity );
        tmp = 1. + pow( vg[2]*fabs(head), vg[4] );
        S = vg[0] + (vg[1]-vg[0]) * pow( tmp, (1.-vg[4])/vg[4] );
        if ( S < 0. ) S = 0.;
        if ( S > 1. ) S = 1.;
        Se = (S - vg[0]) / (vg[1]-vg[0]);
        if ( Se < 0. ) Se = 0.;
        if ( Se > 1. ) Se = 1.;
        // relative permeability (Mualem): krel = Se^gl * [1-(1-Se^m2)^m1]^2
        tmp = 1. - pow( Se, vg[4]/(vg[4]-1.) );
        krel = pow( Se, vg[3] ) * pow( 1. - pow( tmp, (vg[4]-1.)/vg[4] ), 2. );
        if ( krel < 0. ) krel = 0.;
        if ( krel < eps_permeability ) krel = eps_permeability;
        for ( idim=0; idim<ndim; idim++ ) pe[idim] *= krel;
        // non-saturated capacity: dS/dphi_p, analytic derivative
        tmp = 1. + pow( vg[2]*fabs(head), vg[4] );
        dS = (vg[1]-vg[0]) * ((1.-vg[4])/vg[4]) * pow( tmp, (1.-vg[4])/vg[4]-1. )
             * vg[4] * pow( vg[2], vg[4] ) * pow( fabs(head), vg[4]-1. );
        if ( head<0. ) dS = -dS;
        if ( fabs(head)<1.e-12 ) dS = 0.;
        C += por * dS;
      }
      // store the saturation in the groundflow_saturation dof (stored, not solved)
      if ( groundflow_saturation ) {
        for ( inol=0; inol<nnol; inol++ ) {
          inod = nodes[inol];
          node_dof = db_dbl( NODE_DOF, inod, VERSION_NEW );
          if ( dens*gravity > 0. ) {
            head = -node_dof[pres_indx] / ( dens * gravity );
            tmp = 1. + pow( vg[2]*fabs(head), vg[4] );
            S = vg[0] + (vg[1]-vg[0]) * pow( tmp, (1.-vg[4])/vg[4] );
            if ( S < 0. ) S = 0.;
            if ( S > 1. ) S = 1.;
            node_dof[gsat_indx] = S;
          }
        }
      }
    }
  }

  (void)old_unknowns; (void)coord_ip; (void)ddum;
}
