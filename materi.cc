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

#define EPS_H 1.e-8

void materi( long int element, long int gr, long int name, long int nnol, 
  long int npoint, long int nodes[], long int plasti_on_boundary,
  double coord_ip[], double old_coord[], long int ipoint,
  double new_dof[], double h[], double new_d[], double new_b[],
  double volume, double old_unknowns[], 
  double new_unknowns[], double old_grad_old_unknowns[], 
  double old_grad_new_unknowns[], double new_grad_new_unknowns[],
  double element_lhside[], double element_matrix[],
  double element_rhside[], double element_residue[], 
  double tendon_element_rhside[])

{
  long int i=0, j=0, idim=0, jdim=0, kdim=0,
    inol=0, jnol=0, m=0, n=0, indx=0, ipuknwn=0, iuknwn=0, jpuknwn=0, 
    swit=0, indxi=0, indxj=0, indx1=0, indx2=0, memory=-UPDATED, 
    ind_ddsdde=0, ldum=0, idum[1],
    sri=-NO, sri_on=0,
    undr_active=0, undr_apply=-YES, undr_icontrol=0, undr_len=0, undr_len2=0, undr_k=0;
  double rdum=0., dens=0., dtime=0., materi_expansion_linear=0., 
    materi_expansion_volume=0., temp=0., tmp=0., damping=0., fac=0, 
    plasti_heatgeneration=0., viscosity_heatgeneration=0.,
    viscosity=0., old_damage=0., new_damage=0., old_kappa=0., new_kappa=0., 
    old_kapsh=0., new_kapsh=0., 
    old_cap1pc=0., new_cap1pc=0.,
    old_f=0., new_f=0., void_fraction=0., new_pres=0., old_substeps=0., new_substeps=0.,
    softvar_nonl=0, softvar_l=0, 
    md_factor=1., // materi_dynamic/control_materi_dynamic momentum factor
    static_pressure=0., total_pressure=0., location=0.,
    J=0., ddum[1], direct_normal[MDIM], *force_gravity=NULL, 
    undr_p=0., undr_C=0., undr_buf[MPOINT],
    activation_factor=1., activation_stiff=1.,
    sri_g=0., sri_g2=0., sri_g3=0., sri_volfac=1., sri_detj_weight=0.,
    sri_coord_center[MDIM],
    *old_deften=NULL, *new_deften=NULL, *inv_deften=NULL,
    *old_epe=NULL, *inc_epe=NULL, *new_epe=NULL,
    *old_epp=NULL, *inc_epp=NULL, 
    *old_ept=NULL, *inc_ept=NULL, *new_ept=NULL,
    *old_rot=NULL, *inc_rot=NULL, *new_rot=NULL, 
    *inv_rot=NULL, *old_sig=NULL, *rot=NULL, 
    *new_sig=NULL, *total_new_sig=NULL, 
    *rotated_old_sig=NULL, *rotated_new_sig=NULL, 
    *rotated_old_epi=NULL, *rotated_new_epi=NULL, 
    *rotated_old_rho=NULL, *rotated_new_rho=NULL, 
    *rotated_old_msig=NULL, *rotated_new_msig=NULL, 
    *old_msig=NULL, *new_msig=NULL,
    *old_rho=NULL, *new_rho=NULL, 
    *old_epi=NULL, *new_epi=NULL, 
    *old_hisv=NULL, *new_hisv=NULL,
    *ddsdde=NULL, *ddsdde_tendon=NULL,
    *ddsdde_total=NULL, *sigvec=NULL,
    *work=NULL, *stiffness=NULL, *force=NULL,
    *stiffness_shear=NULL, *b_shear=NULL,
    *dbl_array=NULL,
    *new_sig_nonrot=NULL;

  static long int sri_warning_done = 0;

  swit = set_swit(element,-1,"materi");
  if ( swit ) pri( "In routine MATERI." );

  n = nnol*MDIM*nnol*MDIM +
    nnol*ndim*nnol*ndim + nnol*ndim +
    MDIM + MDIM*MDIM + MDIM*MDIM + MDIM*MDIM + MDIM*MDIM +
    MDIM*MDIM + MDIM*MDIM + MDIM*MDIM + MDIM*MDIM + MDIM*MDIM +
    MDIM*MDIM + MDIM*MDIM + MDIM*MDIM + MDIM*MDIM + MDIM*MDIM +
    MDIM*MDIM + MDIM*MDIM + MDIM*MDIM + MDIM*MDIM + MDIM*MDIM +
    MDIM*MDIM + MDIM*MDIM + MDIM*MDIM +
    MDIM*MDIM + MDIM*MDIM + MDIM*MDIM +
    MDIM*MDIM + MDIM*MDIM +
    materi_maxwell_stress*MDIM*MDIM + materi_maxwell_stress*MDIM*MDIM +
    materi_maxwell_stress*MDIM*MDIM + materi_maxwell_stress*MDIM*MDIM + MDIM*MDIM + MDIM*MDIM +
    nuknwn + nuknwn + MSTRAIN*MSTRAIN + MSTRAIN*MSTRAIN +
    MSTRAIN*MSTRAIN + MSTRAIN*MSTRAIN +
    nnol*MDIM*nnol*MDIM + nnol*ndim*nnol*ndim + nnol*ndim + MDIM*MDIM +
    nnol*ndim*nnol*ndim + nnol*ndim;
  dbl_array = get_new_dbl(n);
  assert( indx<=n );
 
  indx = 0;
  work = &dbl_array[indx]; indx += nnol*MDIM*nnol*MDIM;
  stiffness_shear = &dbl_array[indx]; indx += nnol*ndim*nnol*ndim;
  b_shear = &dbl_array[indx]; indx += nnol*ndim;
  force = &dbl_array[indx]; indx += nnol*ndim;
  force_gravity = &dbl_array[indx]; indx += MDIM;
  old_deften = &dbl_array[indx]; indx += MDIM*MDIM;
  new_deften = &dbl_array[indx]; indx += MDIM*MDIM;
  inv_deften = &dbl_array[indx]; indx += MDIM*MDIM;
  old_epe = &dbl_array[indx]; indx += MDIM*MDIM;
  inc_epe = &dbl_array[indx]; indx += MDIM*MDIM;
  new_epe = &dbl_array[indx]; indx += MDIM*MDIM;
  old_epp = &dbl_array[indx]; indx += MDIM*MDIM;
  inc_epp = &dbl_array[indx]; indx += MDIM*MDIM;
  old_ept = &dbl_array[indx]; indx += MDIM*MDIM;
  inc_ept = &dbl_array[indx]; indx += MDIM*MDIM;
  new_ept = &dbl_array[indx]; indx += MDIM*MDIM;
  old_rot = &dbl_array[indx]; indx += MDIM*MDIM;
  inc_rot = &dbl_array[indx]; indx += MDIM*MDIM;
  new_rot = &dbl_array[indx]; indx += MDIM*MDIM;
  inv_rot = &dbl_array[indx]; indx += MDIM*MDIM;
  old_sig = &dbl_array[indx]; indx += MDIM*MDIM;
  rot = &dbl_array[indx]; indx += MDIM*MDIM;
  new_sig = &dbl_array[indx]; indx += MDIM*MDIM;
  total_new_sig = &dbl_array[indx]; indx += MDIM*MDIM;
  rotated_old_sig = &dbl_array[indx]; indx += MDIM*MDIM;
  rotated_new_sig = &dbl_array[indx]; indx += MDIM*MDIM;
  rotated_old_epi = &dbl_array[indx]; indx += MDIM*MDIM;
  rotated_new_epi = &dbl_array[indx]; indx += MDIM*MDIM;
  rotated_old_rho = &dbl_array[indx]; indx += MDIM*MDIM;
  rotated_new_rho = &dbl_array[indx]; indx += MDIM*MDIM;
  rotated_old_msig = &dbl_array[indx]; indx += materi_maxwell_stress*MDIM*MDIM;
  rotated_new_msig = &dbl_array[indx]; indx += materi_maxwell_stress*MDIM*MDIM;
  old_msig = &dbl_array[indx]; indx += materi_maxwell_stress*MDIM*MDIM;
  new_msig = &dbl_array[indx]; indx += materi_maxwell_stress*MDIM*MDIM;
  old_rho = &dbl_array[indx]; indx += MDIM*MDIM;
  new_rho = &dbl_array[indx]; indx += MDIM*MDIM;
  old_epi = &dbl_array[indx]; indx += MDIM*MDIM;
  new_epi = &dbl_array[indx]; indx += MDIM*MDIM;
  old_hisv = &dbl_array[indx]; indx += nuknwn;
  new_hisv = &dbl_array[indx]; indx += nuknwn;
  ddsdde = &dbl_array[indx]; indx += MSTRAIN*MSTRAIN;
  ddsdde_tendon = &dbl_array[indx]; indx += MSTRAIN*MSTRAIN;
  ddsdde_total = &dbl_array[indx]; indx += MSTRAIN*MSTRAIN;
  sigvec = &dbl_array[indx]; indx += MSTRAIN*MSTRAIN;
  work = &dbl_array[indx]; indx += nnol*MDIM*nnol*MDIM;
  stiffness = &dbl_array[indx]; indx += nnol*ndim*nnol*ndim;
  force = &dbl_array[indx]; indx += nnol*ndim;
  new_sig_nonrot = &dbl_array[indx]; indx += MDIM*MDIM;
  assert( indx<=n );

  array_set( inc_epe, 0., MDIM*MDIM );
  array_set( old_ept, 0., MDIM*MDIM );
  array_set( old_epp, 0., MDIM*MDIM );
  array_set( old_epe, 0., MDIM*MDIM );
  array_set( old_epi, 0., MDIM*MDIM );
  array_set( old_rot, 0., MDIM*MDIM );
  array_set( old_sig, 0., MDIM*MDIM );
  array_set( new_sig, 0., MDIM*MDIM );
  array_set( new_sig_nonrot, 0., MDIM*MDIM );
  array_set( old_rho, 0., MDIM*MDIM );
  array_set( old_deften, 0., MDIM*MDIM );
  array_set( old_hisv, 0., nuknwn );
  array_set( new_hisv, 0., nuknwn );
  array_set( old_msig, 0., materi_maxwell_stress*MDIM*MDIM );

  db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET );
  get_group_data( GROUP_MATERI_DAMPING, gr, element, new_unknowns, 
    &damping, ldum, GET_IF_EXISTS );
  get_group_data( GROUP_MATERI_PLASTI_HEATGENERATION, gr, element,
    new_unknowns, &plasti_heatgeneration, ldum, GET_IF_EXISTS );
  // group_materi_plasti_heat_generation (manual Professional 6.703):
  // the Professional name (with underscores) of the legacy GNU keyword
  // group_materi_plasti_heatgeneration. Dual read at the consumption
  // point (get_group_data is allocation-free here): the Professional
  // name wins when both exist.
  get_group_data( GROUP_MATERI_PLASTI_HEAT_GENERATION, gr, element,
    new_unknowns, &plasti_heatgeneration, ldum, GET_IF_EXISTS );
  db( GROUP_MATERI_MEMORY, gr, &memory, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  dens = get_materi_density( element, gr, nnol, nodes, new_unknowns );
  force_gravity_calculate( force_gravity );

  // group_element_selective_reduced_integration (SRI, Hughes):
  // opt-in fix for the shear locking of the bilinear quad4 (2D) and
  // the trilinear hex8 (3D) in bending. When -yes, the shear part of
  // the constitutive matrix (gamma_xy for the quad4; gamma_xy, gamma_xz
  // and gamma_yz for the hex8) is integrated with 1 Gauss point at the
  // element centroid while the normal/volumetric part keeps the full
  // rule - which for the SRI group becomes the classic Gauss rule
  // (2x2 / 2x2x2; pol() switches it; the codebase default Lobatto
  // corner rule would cap the benefit at ~0.35x of the exact section
  // moment, measured). The scope decision is in sri_active()
  // (miscel.cc): LINEAR ELASTICITY only, because the split
  // D = D_norm + D_shear is exact only when the tangent is constant
  // over the element (the reduced point has no material state of its
  // own); plasticity, damage, maxwell, large displacement and
  // axisymmetric groups ignore the keyword with a one-time warning.
  // WARNING (measured 2026-08-29): the hex8 SRI retains zero-energy
  // modes (3 twist modes of the isolated element; section-warping
  // modes of a mesh), so LOADED hex8 configurations assemble a
  // singular momentum matrix - see sri_active() and the developer
  // manual. The patch tests (constant states) and the rigid modes
  // stay exact.
  db( GROUP_ELEMENT_SELECTIVE_REDUCED_INTEGRATION, gr, &sri, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  sri_on = sri_active( element, gr, name, nnol );
  if ( sri==-YES && materi_stress && !sri_on && !sri_warning_done ) {
    pri( "Warning: group_element_selective_reduced_integration ignored: "
      "only linear elastic 2D quad4 / 3D hex8 groups are supported (no "
      "axisymmetry, no materi_displacement, no plasticity/damage/maxwell)." );
    sri_warning_done = 1;
  }

  // mesh_activate_gravity_time: the gravity is gradually activated for the
  // element (bottom-to-top interpolation). With method 2 the element stays
  // active (reduced stiffness) but without gravity until activation.
  activation_factor = mesh_activate_gravity_factor( element, gr, nnol, nodes,
    &activation_stiff );
  for ( idim=0; idim<ndim; idim++ )
    force_gravity[idim] *= activation_factor;

  if ( condif_temperature ) {
    get_group_data( GROUP_MATERI_EXPANSION_VOLUME, gr, element, new_unknowns, 
      &materi_expansion_volume, ldum, GET_IF_EXISTS );
    get_group_data( GROUP_MATERI_EXPANSION_LINEAR, gr, element, new_unknowns, 
      &materi_expansion_linear, ldum, GET_IF_EXISTS );
    temp = new_unknowns[temp_indx];
    dens = (1.-materi_expansion_volume*temp) * dens;
  }

    // get old stresses, etc.
  for ( idim=0; idim<MDIM; idim++ ) {
    for ( jdim=0; jdim<MDIM; jdim++ ) {
      indx = idim*MDIM + jdim;
      if ( materi_stress ) {
       tmp = old_unknowns[stres_indx+stress_indx(idim,jdim)*nder];
       old_sig[indx] = tmp;
       new_sig[indx] = tmp;
      }
      if ( materi_strain_elasti ) old_epe[indx] =
        old_unknowns[epe_indx+stress_indx(idim,jdim)*nder];
      if ( materi_strain_total ) old_ept[indx] = 
        old_unknowns[ept_indx+stress_indx(idim,jdim)*nder];
      for ( m=0; m<materi_maxwell_stress; m++ ) {
        indx = m*MDIM*MDIM+idim*MDIM+jdim;
        old_msig[indx] = old_unknowns[mstres_indx+
          (m*6+stress_indx(idim,jdim))*nder];
      }
      if ( materi_plasti_rho ) old_rho[indx] =
        old_unknowns[rho_indx+stress_indx(idim,jdim)*nder];
      if ( materi_strain_intergranular ) old_epi[indx] =
        old_unknowns[epi_indx+stress_indx(idim,jdim)*nder];
      if ( materi_strain_plasti ) old_epp[indx] =
        old_unknowns[epp_indx+stress_indx(idim,jdim)*nder];
    }
  }

  if ( materi_history_variables ) {
    for ( i=0; i<materi_history_variables; i++ ) {
      iuknwn = hisv_indx + i;
      old_hisv[i] = old_unknowns[iuknwn];
      new_hisv[i] = new_unknowns[iuknwn];
    }
  }
  if ( materi_plasti_kappa ) {
    iuknwn = kap_indx;
    old_kappa = old_unknowns[iuknwn];
    new_kappa = new_unknowns[iuknwn];
  }
  if ( materi_plasti_kappa_shear ) {
    iuknwn = kapsh_indx;
    old_kapsh = old_unknowns[iuknwn];
    new_kapsh = new_unknowns[iuknwn];
  }
  if ( materi_plasti_cap1_history ) {
    iuknwn = cap1_indx;
    old_cap1pc = old_unknowns[iuknwn];
    new_cap1pc = new_unknowns[iuknwn];
  }
  if ( materi_damage ) {
    iuknwn = dam_indx;
    old_damage = old_unknowns[iuknwn];
    new_damage = new_unknowns[iuknwn];
  }
  if ( materi_plasti_f ) {
    iuknwn = f_indx;
    old_f = old_unknowns[iuknwn];
  }
  if ( materi_plasti_incremental_substeps ) {
    iuknwn = substeps_indx;
    old_substeps = old_unknowns[iuknwn];
  }
  if ( materi_plasti_softvar_nonlocal ) {
    iuknwn = svnonloc_indx;
    softvar_nonl = new_unknowns[iuknwn];
  }
  if ( materi_plasti_softvar_local ) {
    iuknwn = svloc_indx;
    softvar_l = new_unknowns[iuknwn];
  }
  
  set_deften_etc( element, gr, nnol, h, old_coord, old_unknowns, 
    new_unknowns, old_grad_old_unknowns, old_grad_new_unknowns, 
    old_deften, new_deften, old_ept, inc_ept, new_ept, 
    old_rot, inc_rot, new_rot );

  // group_materi_undrained_capacity (manual Professional 6.760 + theory
  // 2.2.7): UNDRAINED groundwater analysis WITHOUT the groundwater dof
  // in the system matrix. When the element group carries the capacity C
  // (and control_materi_undrained_apply, 6.153, is not -no), the total
  // groundwater pressure change of the element follows from the
  // groundflow storage equation without permeability, solved on the
  // element level: C * p_dot = div(v_material). Per integration point
  // the pressure increment of this step is the volumetric strain
  // increment divided by C, accumulated over the steps in the record
  // element_intpnt_materi_undrained_pressure (6.441, one value per
  // integration point; the record lives in both db versions like
  // ELEMENT_DOF, so VERSION_NORMAL holds the converged value of the
  // previous step and VERSION_NEW the current iterate - the step-end
  // version copy promotes it). The pressure acts isotropically on the
  // skeleton: the TOTAL stress for the momentum equilibrium is the
  // effective constitutive stress plus (groundflow total pressure +
  // undrained pressure)*I (manual 2.2.7: "the fixed total pressure from
  // the hydraulic pressure heads plus the excessive undrained pressure
  // ... as the full total pressure"), while the stress DOFs keep the
  // EFFECTIVE value. The momentum stiffness gains the volumetric term
  // dt/C*(div w)*(div v) so the linear system carries the undrained
  // stiffness (drained + 1/C, measured on undrained2 of the corpus:
  // sigma = E*eps + eps/C balances the applied force).
  undr_active = 0;
  if ( !find_local_softvar &&
       db_active_index( GROUP_MATERI_UNDRAINED_CAPACITY, gr,
         VERSION_NORMAL ) ) {
    undr_C = db_dbl( GROUP_MATERI_UNDRAINED_CAPACITY, gr,
      VERSION_NORMAL )[0];
    if ( undr_C>0. && materi_stress ) {
      db( ICONTROL, 0, &undr_icontrol, ddum, ldum, VERSION_NORMAL, GET );
      db( CONTROL_MATERI_UNDRAINED_APPLY, undr_icontrol, &undr_apply,
        ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
      if ( undr_apply!=-NO ) {
        undr_active = 1;
        undr_p = 0.;
        undr_len = 0;
        // p_old basis of THIS integration point: the value at the step
        // start (VERSION_NORMAL; VERSION_NEW already holds the current
        // iterate of the other integration points of this element pass)
        if ( db( ELEMENT_INTPNT_MATERI_UNDRAINED_PRESSURE, element, idum,
            undr_buf, undr_len, VERSION_NORMAL, GET_IF_EXISTS ) &&
            ipoint<undr_len )
          undr_p = undr_buf[ipoint];
        undr_p += ( inc_ept[0*MDIM+0] + inc_ept[1*MDIM+1] +
          inc_ept[2*MDIM+2] ) / undr_C;
        // store the iterate: the record is read-modify-written per
        // integration point, so the buffer must preserve the slots
        // written by the other integration points of this element pass
        // (VERSION_NEW) and only the current slot is replaced (the
        // step-end version copy promotes the converged values of the
        // last iteration to VERSION_NORMAL).
        undr_len2 = 0;
        if ( db( ELEMENT_INTPNT_MATERI_UNDRAINED_PRESSURE, element, idum,
            undr_buf, undr_len2, VERSION_NEW, GET_IF_EXISTS ) ) {
          for ( undr_k=0; undr_k<npoint && undr_k<MPOINT; undr_k++ ) {
            if ( undr_k==ipoint ) undr_buf[undr_k] = undr_p;
            else if ( undr_k>=(long int)undr_len2 ) undr_buf[undr_k] = 0.;
          }
          db( ELEMENT_INTPNT_MATERI_UNDRAINED_PRESSURE, element, idum,
            undr_buf, npoint, VERSION_NEW, PUT );
        }
      }
    }
  }

    // back rotate to old configuration
  if      ( memory==-TOTAL || memory==-TOTAL_PIOLA ) {
    if ( !matrix_inverse( old_rot, inv_rot, rdum, MDIM ) ) {
      pri ("Error detected for element ", element );
      pri ("Probably too large distortions." );
      exit_tn_on_error();
    }                     
    if ( materi_stress ) {
      if ( memory==-TOTAL_PIOLA ) {
        if ( !matrix_inverse( old_deften, inv_deften, rdum, MDIM ) ) {
          pri ("Error detected for element ", element );
          pri ("Probably too large distortions." );
          exit_tn_on_error();
        }                     
        J = matrix_determinant( old_deften, MDIM );
        matrix_abat( inv_deften, old_sig, rotated_old_sig, work, MDIM);
        array_multiply( rotated_old_sig, rotated_old_sig, J, MDIM*MDIM );
      }
      else
        matrix_abat( inv_rot, old_sig, rotated_old_sig, work, MDIM );
    }
    if ( materi_plasti_rho ) {
      matrix_abat( inv_rot, old_rho, rotated_old_rho, work, MDIM );
      if ( swit ) pri( "rotated_old_rho", rotated_old_rho, MDIM, MDIM );
    }
    if ( materi_strain_intergranular ) {
      matrix_abat( inv_rot, old_epi, rotated_old_epi, work, MDIM );
      if ( swit ) pri( "rotated_old_epi", rotated_old_epi, MDIM, MDIM );
    }
    if ( materi_maxwell_stress ) {
      for ( m=0; m<materi_maxwell_stress; m++ ) {
        indx = m*MDIM*MDIM;
        matrix_abat( inv_rot, &old_msig[indx], &rotated_old_msig[indx], work, MDIM );
      }
    }
  }
  else if ( memory==-UPDATED || memory==-TOTAL_LINEAR ||
      memory==-UPDATED_WITHOUT_ROTATION || memory==-UPDATED_LINEAR ||
      memory==-UPDATED_AREA ) {
    if ( materi_stress )
      array_move( old_sig, rotated_old_sig, MDIM*MDIM );
    if ( materi_plasti_rho ) 
      array_move( old_rho, rotated_old_rho, MDIM*MDIM );
    if ( materi_strain_intergranular ) 
      array_move( old_epi, rotated_old_epi, MDIM*MDIM );
    if ( materi_maxwell_stress )
      array_move( old_msig, rotated_old_msig, materi_maxwell_stress*MDIM*MDIM );
  }
  else {
    array_set( rotated_old_sig, 0., MDIM*MDIM );
  }
  if ( swit ) {
    if ( materi_stress ) pri( "rotated_old_sig", rotated_old_sig, MDIM, MDIM );
  }

    // set stress, strain
  if ( materi_stress ) {
    // direct-normal plane (group_materi_plasti_*_direct_normal[_automatic]):
    // explicit normal from the keyword, or the element normal computed from
    // the element geometry (cross product of two side-1 edges).
    array_set( direct_normal, 0., MDIM );
    if ( db_active_index( GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_NORMAL, gr,
        VERSION_NORMAL ) ) {
      db( GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_NORMAL, gr, idum,
        direct_normal, ldum, VERSION_NORMAL, GET );
    }
    else if ( db_active_index( GROUP_MATERI_PLASTI_TENSION_DIRECT_NORMAL, gr,
        VERSION_NORMAL ) ) {
      db( GROUP_MATERI_PLASTI_TENSION_DIRECT_NORMAL, gr, idum,
        direct_normal, ldum, VERSION_NORMAL, GET );
    }
    else if ( db_active_index( GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_NORMAL_AUTOMATIC,
        gr, VERSION_NORMAL ) ||
              db_active_index( GROUP_MATERI_PLASTI_TENSION_DIRECT_NORMAL_AUTOMATIC,
        gr, VERSION_NORMAL ) ) {
      // element normal: cross product of the first two element edges
      double e1[MDIM], e2[MDIM];
      double *c0, *c1, *c2;
      c0 = db_dbl( NODE, nodes[0], VERSION_NORMAL );
      c1 = db_dbl( NODE, nodes[1], VERSION_NORMAL );
      c2 = db_dbl( NODE, nodes[2], VERSION_NORMAL );
      for ( i=0; i<3; i++ ) {
        e1[i] = c1[i] - c0[i];
        e2[i] = c2[i] - c0[i];
      }
      direct_normal[0] = e1[1]*e2[2] - e1[2]*e2[1];
      direct_normal[1] = e1[2]*e2[0] - e1[0]*e2[2];
      direct_normal[2] = e1[0]*e2[1] - e1[1]*e2[0];
      array_normalize( direct_normal, 3 );
    }
    set_stress( element, gr, plasti_on_boundary, coord_ip,
      old_unknowns, new_unknowns,
      old_grad_old_unknowns, new_grad_new_unknowns, 
      rotated_old_sig, new_sig, 
      rotated_old_msig, new_msig, inc_ept, new_ept,
      old_epe, inc_epe, old_epp, inc_epp, old_rho, new_rho, 
      old_epi, new_epi, old_hisv, new_hisv, 
      old_damage, new_damage, old_kappa, new_kappa, old_kapsh, new_kapsh,
      old_cap1pc, new_cap1pc,
      new_f, new_substeps,
      old_deften, new_deften, inc_rot,
      ddsdde, viscosity, viscosity_heatgeneration, softvar_nonl, softvar_l,
      direct_normal);
    tendons( element, gr, nnol, npoint, volume, new_d, old_unknowns, new_unknowns,
      new_rot, inc_ept, tendon_element_rhside, ddsdde_tendon );
    array_add( ddsdde, ddsdde_tendon, ddsdde_total, MSTRAIN*MSTRAIN );
    if ( swit ) pri( "ddsdde_total", ddsdde_total, MSTRAIN, MSTRAIN );
  }

  if(find_local_softvar) {
  for ( inol=0; inol<nnol; inol++ ) {
    if ( scalar_dabs(h[inol])>EPS_H ) {
     if ( materi_plasti_softvar_nonlocal && materi_plasti_softvar_local ) {
       		//added for options_element_dof
       if(options_element_dof==-YES) new_unknowns[svloc_indx/nder] = softvar_l;
     }
    }
   }
  }

  //Not used when searching for local values of softening variable 
  if(!find_local_softvar) {

  array_move( new_sig, new_sig_nonrot, MDIM*MDIM );
    // rotate to new configuration
  if ( memory==-UPDATED || memory==-UPDATED_LINEAR || memory==-UPDATED_AREA ||
       memory==-TOTAL || memory==-TOTAL_PIOLA ) {
    if      ( memory==-UPDATED || memory==-UPDATED_LINEAR ||
              memory==-UPDATED_AREA )
      array_move( inc_rot, rot, MDIM*MDIM );
    else {
      assert( memory==-TOTAL || memory==-TOTAL_PIOLA );
      array_move( new_rot, rot, MDIM*MDIM );
    }
    if ( materi_stress ) {
      if (  memory==-TOTAL_PIOLA ) {
        matrix_abat( new_deften, new_sig, rotated_new_sig, work, MDIM );
        J = matrix_determinant( new_deften, MDIM );
        if ( J<=0. ) {
          pri ("Error detected for element ", element );
          pri ("Non-positive jacobian." );
          pri ("Probably too large distortions." );
          exit_tn_on_error();
        }                     
        array_multiply( rotated_new_sig, rotated_new_sig, 1./J, MDIM*MDIM );
      }
      else
        matrix_abat( rot, new_sig, rotated_new_sig, work, MDIM );
      array_move( rotated_new_sig, new_sig, MDIM*MDIM );
      if ( swit ) pri( "new_sig", new_sig, MDIM, MDIM );
    }
    if ( materi_plasti_rho ) {
      matrix_abat( rot, new_rho, rotated_new_rho, work, MDIM );
      array_move( rotated_new_rho, new_rho, MDIM*MDIM );
    }
    if ( materi_strain_intergranular ) {
      matrix_abat( rot, new_epi, rotated_new_epi, work, MDIM );
      array_move( rotated_new_epi, new_epi, MDIM*MDIM );
    }
    if ( materi_maxwell_stress ) {
      for ( m=0; m<materi_maxwell_stress; m++ ) {
        indx = m*MDIM*MDIM;
        matrix_abat( rot, &new_msig[indx], &rotated_new_msig[indx], work, MDIM );
      }
      array_move( rotated_new_msig, new_msig, materi_maxwell_stress*MDIM*MDIM );
    }
  }

    // forces on nodes
  if ( materi_stress ) {
    array_move( new_sig, total_new_sig, MDIM*MDIM );
    if ( groundflow_pressure ) {
      double gpf = 1.;
      db( GROUNDFLOW_PRESSURE_FACTOR, 0, idum, &gpf, ldum,
        VERSION_NORMAL, GET_IF_EXISTS );
      new_pres = new_unknowns[pres_indx];
      // groundflow_phreatic_coord total-pressure substitution: pore
      // pressure acting on the skeleton becomes the TOTAL pressure when a
      // single groundflow_phreatic_level covers the integration point
      // (calibrated on ground15/16). NOT applied when
      // groundflow_phreatic_level_multiple records exist: their domains
      // are often combined with head-prescribing bounda_dof -pres loads
      // (the reservoir statics of ground19 of the corpus), where adding
      // the level static to the solved head double counts rho*g*z_level
      // and destabilizes the coupled system (Bi-CG breakdown). The
      // multiple-level mechanics coupling is pending calibration (see
      // SEGUIMIENTO ground8 row).
      if ( !groundflow_phreatic_level_multiple_active() &&
           groundflow_phreatic_coord( -1, coord_ip, new_unknowns, 
             total_pressure, static_pressure, location, NULL ) ) new_pres = total_pressure;
      // group_groundflow_total_pressure_tension: if the largest eigenvalue of
      // materi_strain_plastic_tension exceeds plastic_tension_minimum, use the
      // static water pore pressure determined from water_height (when it is
      // larger in absolute value than the pore pressure from the groundflow
      // equation). Takes care that in cracks in concrete the largest water
      // pressure from an environment is used.
      if ( db_active_index( GROUP_GROUNDFLOW_TOTAL_PRESSURE_TENSION, gr,
          VERSION_NORMAL ) ) {
        double gptt[2], epp_princ[3], epp_t[MDIM*MDIM], dens_water=0.;
        db( GROUP_GROUNDFLOW_TOTAL_PRESSURE_TENSION, gr, idum, gptt, ldum,
          VERSION_NORMAL, GET );
        db( GROUNDFLOW_DENSITY, 0, idum, &dens_water, ldum, VERSION_NORMAL,
          GET_IF_EXISTS );
        for ( idim=0; idim<MDIM; idim++ ) {
          for ( jdim=0; jdim<MDIM; jdim++ ) {
            if ( materi_strain_plasti )
              epp_t[idim*MDIM+jdim] =
                new_unknowns[epp_indx+stress_indx(idim,jdim)*nder];
            else
              epp_t[idim*MDIM+jdim] = 0.;
          }
        }
        matrix_eigenvalues( epp_t, epp_princ );
        tmp = epp_princ[0];
        if ( epp_princ[1] > tmp ) tmp = epp_princ[1];
        if ( epp_princ[2] > tmp ) tmp = epp_princ[2];
        if ( tmp > gptt[0] ) {
          double static_pres =
            force_gravity[ndim-1] * dens_water * ( gptt[1] - coord_ip[ndim-1] );
          if ( scalar_dabs(static_pres) > scalar_dabs(new_pres) ) new_pres = static_pres;
        }
      }
      new_pres *= gpf;
      for ( idim=0; idim<MDIM; idim++ ) total_new_sig[idim*MDIM+idim] += new_pres;
    }
    if ( undr_active ) {
      // undrained capacity: the excessive undrained pressure of this
      // step joins the fixed groundflow pressure as the full total
      // pressure acting on the skeleton (manual 2.2.7); the stress DOFs
      // keep the effective value (they are written from new_sig below).
      for ( idim=0; idim<MDIM; idim++ )
        total_new_sig[idim*MDIM+idim] += undr_p;
    }
    for ( idim=0; idim<MDIM; idim++ ) {
      for ( jdim=idim; jdim<MDIM; jdim++ ) {
        indx = stress_indx(idim,jdim);
        sigvec[indx] = total_new_sig[idim*MDIM+jdim];
      }
    }
    if ( sri_on ) {
      // SRI: D = D_norm + D_shear. The shear diagonal entries of D
      // (the engineering-shear moduli G: index stress_indx(0,1)=1 for
      // gamma_xy, (0,2)=2 for gamma_xz and (1,2)=4 for gamma_yz in 3D;
      // index 1 only in 2D) are excluded from the full integration and
      // added back with 1 point at the centroid below. For linear
      // elasticity the tangent is constant over the element, so the
      // shear moduli captured here are exact for the reduced point
      // (isotropic: the shear rows of D are decoupled, so zeroing the
      // diagonal is the exact D_norm).
      // Element-consistent momentum feedback (DIAG lot C/D, fix D-c):
      // the shear part of the feedback stress must use the SAME reduced
      // integration as the momentum matrix (1 point at the centroid),
      // otherwise the full-rule shear re-enters the right-hand side and
      // the staggered fixed point converges to K_full*u = P (the locked
      // state) instead of K_elem*u = P (the element solution) - the SRI
      // benefit is cancelled at equilibrium. The current-iterate shear
      // increments of the feedback stress are zeroed here (the old
      // shear prestresses are kept) and the reduced 1-point shear
      // internal force -dt*K_shear*v is added to the momentum RHS in
      // the velocity block below. NOTE the strain indexing: inc_ept is
      // the MDIMxMDIM TENSOR strain (inc_ept[idim*MDIM+jdim]), while
      // sigvec is Voigt-indexed (stress_indx). In 2D the indices
      // coincide; in 3D the yz entry differs (tensor 1*MDIM+2 = 5 vs
      // Voigt stress_indx(1,2) = 4), so the strain MUST use the tensor
      // index idim*MDIM+jdim.
      if ( name==-QUAD4 ) {
        sri_g = ddsdde_total[1*MSTRAIN+1];
        sigvec[stress_indx(0,1)] -= 2. * sri_g * inc_ept[stress_indx(0,1)];
        ddsdde_total[1*MSTRAIN+1] = 0.;
      }
      else {
        assert( name==-HEX8 );
        sri_g  = ddsdde_total[1*MSTRAIN+1];   // G for gamma_xy
        sri_g2 = ddsdde_total[2*MSTRAIN+2];   // G for gamma_xz
        sri_g3 = ddsdde_total[4*MSTRAIN+4];   // G for gamma_yz
        sigvec[stress_indx(0,1)] -= 2.*sri_g *inc_ept[0*MDIM+1];
        sigvec[stress_indx(0,2)] -= 2.*sri_g2*inc_ept[0*MDIM+2];
        sigvec[stress_indx(1,2)] -= 2.*sri_g3*inc_ept[1*MDIM+2];
        ddsdde_total[1*MSTRAIN+1] = 0.;
        ddsdde_total[2*MSTRAIN+2] = 0.;
        ddsdde_total[4*MSTRAIN+4] = 0.;
      }
    }
    // materi_dynamic (manual Professional 6.800) / control_materi_dynamic
    // (6.141): blend the stress at time t (old_sig, previous converged
    // step at this integration point) with the stress at time t+dt
    // (sigvec of the current iterate): sigma = (1-factor)*sigma_t +
    // factor*sigma_{t+dt}. Default factor = 1 (fully implicit, the
    // historic GNU scheme). A factor < 1 makes the scheme less
    // implicit and thus reduces numerical damping (dynamics); factor 0
    // freezes the internal force at the old stress (forward-Euler-like
    // momentum). The momentum stiffness is scaled by the same factor
    // (the tangent of factor*sigma_{t+dt}).
    {
      double materi_dynamic_factor = 1.;
      db( MATERI_DYNAMIC, 0, idum, &materi_dynamic_factor, ldum,
        VERSION_NORMAL, GET_IF_EXISTS );
      long int icontrol_md = 0;
      db( ICONTROL, 0, &icontrol_md, ddum, ldum, VERSION_NORMAL,
        GET_IF_EXISTS );
      if ( db_active_index( CONTROL_MATERI_DYNAMIC, icontrol_md,
          VERSION_NORMAL ) )
        db( CONTROL_MATERI_DYNAMIC, icontrol_md, idum,
          &materi_dynamic_factor, ldum, VERSION_NORMAL, GET );
      md_factor = materi_dynamic_factor;
      if ( materi_dynamic_factor<1. ) {
        if ( materi_dynamic_factor<0. || materi_dynamic_factor>1. )
          db_error( MATERI_DYNAMIC, 0 );
        for ( idim=0; idim<MDIM; idim++ ) {
          for ( jdim=idim; jdim<MDIM; jdim++ ) {
            long int indx_md = stress_indx(idim,jdim);
            sigvec[indx_md] = (1.-materi_dynamic_factor) *
              old_sig[idim*MDIM+jdim] +
              materi_dynamic_factor * sigvec[indx_md];
          }
        }
      }
    }
    matrix_atb( new_b, sigvec, force, MSTRAIN, nnol*ndim, 1 );
    matrix_atba( new_b, ddsdde_total, stiffness, work, MSTRAIN, nnol*ndim );
    if ( swit ) {
      pri( "force", force, nnol*ndim );
      pri( "stiffness", stiffness, nnol*ndim, nnol*ndim );
      pri( "sigvec", sigvec, MSTRAIN );
    }
    if ( sri_on ) {
      // SRI: stiffness of the shear term integrated with 1 Gauss point
      // at the element centroid. For the bilinear quad4 the local
      // derivatives at iso (0,0) are dN/dxi = xi_i/4 and
      // dN/deta = eta_i/4 with the LOCAL node coordinates (xi_i,
      // eta_i) = ((2*(i%2)-1), (2*(i/2)-1)) - the row-major-from-
      // bottom convention that polynom.cc uses (node 0 = (-1,-1),
      // 1 = (+1,-1), 2 = (-1,+1), 3 = (+1,+1); NOT the textbook
      // counter-clockwise order). For the trilinear hex8 the local
      // derivatives at iso (0,0,0) are dN/dxi = xi_i/8 etc. with
      // (xi_i, eta_i, zeta_i) = ((2*(i%2)-1), (2*((i/2)%2)-1),
      // (2*(i/4)-1)). Deriving them from the index (instead of
      // hardcoding) keeps the centroid B consistent with the B
      // matrices of pol() for ANY physical node numbering; the
      // Jacobian at the centroid handles distorted elements exactly.
      // The 1x1 / 1x1x1 Gauss weight on [-1,1]^2 / [-1,1]^3 is 4 / 8.
      if ( name==-QUAD4 ) {
        double jac[4], invjac[4], dnxi[4], dnet[4], detj_center=0.;
        for ( inol=0; inol<nnol; inol++ ) {
          dnxi[inol] = 0.25 * ( (inol%2) ? 1. : -1. );
          dnet[inol] = 0.25 * ( (inol/2) ? 1. : -1. );
        }
        for ( idim=0; idim<ndim; idim++ ) sri_coord_center[idim] = 0.;
        array_set( jac, 0., 4 );
        for ( inol=0; inol<nnol; inol++ ) {
          jac[0] += dnxi[inol] * old_coord[inol*ndim+0];
          jac[1] += dnxi[inol] * old_coord[inol*ndim+1];
          jac[2] += dnet[inol] * old_coord[inol*ndim+0];
          jac[3] += dnet[inol] * old_coord[inol*ndim+1];
          for ( idim=0; idim<ndim; idim++ )
            sri_coord_center[idim] += 0.25 * old_coord[inol*ndim+idim];
        }
        detj_center = jac[0]*jac[3] - jac[1]*jac[2];
        if ( scalar_dabs(detj_center)<TINY ) {
          array_set( stiffness_shear, 0., nnol*ndim*nnol*ndim );
        }
        else {
          invjac[0] =  jac[3]/detj_center; invjac[1] = -jac[1]/detj_center;
          invjac[2] = -jac[2]/detj_center; invjac[3] =  jac[0]/detj_center;
          for ( inol=0; inol<nnol; inol++ ) {
            // dN/dy and dN/dx at the centroid
            b_shear[inol*ndim+0] = invjac[2]*dnxi[inol] + invjac[3]*dnet[inol];
            b_shear[inol*ndim+1] = invjac[0]*dnxi[inol] + invjac[1]*dnet[inol];
          }
          sri_volfac = 1.;
          volume_factor( gr, sri_coord_center, sri_volfac );
          sri_detj_weight = 4. * detj_center;
          for ( i=0; i<nnol*ndim; i++ )
            for ( j=0; j<nnol*ndim; j++ )
              stiffness_shear[i*nnol*ndim+j] = sri_g * sri_volfac *
                sri_detj_weight * b_shear[i] * b_shear[j];
        }
        if ( swit ) pri( "stiffness_shear", stiffness_shear, nnol*ndim, nnol*ndim );
      }
      else {
        assert( name==-HEX8 );
        // 3D: the three engineering shear rows of B at the centroid:
        //   gamma_xy: [dN/dy, dN/dx, 0]
        //   gamma_xz: [dN/dz, 0, dN/dx]
        //   gamma_yz: [0, dN/dz, dN/dy]
        // (the polynom.cc:549-599 B convention; the local derivatives
        // from the index, the physical dN/dx_idim = sum_jdim
        // invjac[idim][jdim]*dN/dxi_jdim with the transpose Jacobian
        // convention of pol()/calcul_force.cc).
        double jac3[9], invjac3[9], dnxi[8], dnet[8], dnzt[8],
          b_xy[24], b_xz[24], b_yz[24], detj_center=0.;
        for ( inol=0; inol<nnol; inol++ ) {
          dnxi[inol] = 0.125 * ( 2.*(inol%2) - 1. );
          dnet[inol] = 0.125 * ( 2.*((inol/2)%2) - 1. );
          dnzt[inol] = 0.125 * ( 2.*(inol/4) - 1. );
        }
        for ( idim=0; idim<ndim; idim++ ) sri_coord_center[idim] = 0.;
        array_set( jac3, 0., 9 );
        for ( inol=0; inol<nnol; inol++ ) {
          jac3[0] += dnxi[inol]*old_coord[inol*ndim+0];
          jac3[1] += dnxi[inol]*old_coord[inol*ndim+1];
          jac3[2] += dnxi[inol]*old_coord[inol*ndim+2];
          jac3[3] += dnet[inol]*old_coord[inol*ndim+0];
          jac3[4] += dnet[inol]*old_coord[inol*ndim+1];
          jac3[5] += dnet[inol]*old_coord[inol*ndim+2];
          jac3[6] += dnzt[inol]*old_coord[inol*ndim+0];
          jac3[7] += dnzt[inol]*old_coord[inol*ndim+1];
          jac3[8] += dnzt[inol]*old_coord[inol*ndim+2];
          for ( idim=0; idim<ndim; idim++ )
            sri_coord_center[idim] += 0.125 * old_coord[inol*ndim+idim];
        }
        detj_center = jac3[0]*( jac3[4]*jac3[8] - jac3[5]*jac3[7] )
                    - jac3[1]*( jac3[3]*jac3[8] - jac3[5]*jac3[6] )
                    + jac3[2]*( jac3[3]*jac3[7] - jac3[4]*jac3[6] );
        if ( scalar_dabs(detj_center)<TINY ||
             !matrix_inverse( jac3, invjac3, detj_center, 3 ) ) {
          array_set( stiffness_shear, 0., nnol*ndim*nnol*ndim );
        }
        else {
          for ( inol=0; inol<nnol; inol++ ) {
            // physical derivatives dN/dx_idim = sum_jdim
            // invjac3[idim][jdim] * (dN/dxi_jdim at the centroid)
            double p3[3], dnx=0., dny=0., dnz=0.;
            p3[0] = dnxi[inol]; p3[1] = dnet[inol]; p3[2] = dnzt[inol];
            for ( jdim=0; jdim<3; jdim++ ) {
              dnx += invjac3[0*3+jdim]*p3[jdim];
              dny += invjac3[1*3+jdim]*p3[jdim];
              dnz += invjac3[2*3+jdim]*p3[jdim];
            }
            b_xy[inol*3+0] = dny;  b_xy[inol*3+1] = dnx;  b_xy[inol*3+2] = 0.;
            b_xz[inol*3+0] = dnz;  b_xz[inol*3+1] = 0.;   b_xz[inol*3+2] = dnx;
            b_yz[inol*3+0] = 0.;   b_yz[inol*3+1] = dnz;  b_yz[inol*3+2] = dny;
          }
          sri_volfac = 1.;
          volume_factor( gr, sri_coord_center, sri_volfac );
          sri_detj_weight = 8. * detj_center;
          for ( i=0; i<nnol*ndim; i++ )
            for ( j=0; j<nnol*ndim; j++ )
              stiffness_shear[i*nnol*ndim+j] = sri_volfac *
                sri_detj_weight * ( sri_g *b_xy[i]*b_xy[j] +
                  sri_g2*b_xz[i]*b_xz[j] + sri_g3*b_yz[i]*b_yz[j] );
        }
        if ( swit ) pri( "stiffness_shear", stiffness_shear, nnol*ndim, nnol*ndim );
      }
    }
  }

    // new elastic strains
  array_add( old_epe, inc_epe, new_epe, MDIM*MDIM );

    // add to right hand side and left hand side
  fac = ((double)nnol)/2.;
  for ( inol=0; inol<nnol; inol++ ) {

      // velocity
    for ( idim=0; idim<ndim; idim++ ) {
      ipuknwn = vel_indx/nder+idim;
      indx = inol*npuknwn + ipuknwn;
      indxi = inol*npuknwn + ipuknwn;
      iuknwn = vel_indx + idim*nder;
        // damping
      tmp = - h[inol] * damping * new_grad_new_unknowns[idim*nuknwn+iuknwn];
      element_rhside[indx] += volume * tmp;
      if ( residue ) element_residue[indx] -= tmp;
      for ( jnol=0; jnol<nnol; jnol++ ) {
        indxj = jnol*npuknwn + ipuknwn;
        tmp = volume * h[inol] * damping * new_d[idim*nnol+jnol];
        element_matrix[indxi*nnol*npuknwn+indxj] += tmp;
        if ( jnol==inol ) element_lhside[indx] += tmp;
      }
        // force_gravity
      tmp = h[inol] * dens * force_gravity[idim];
      element_rhside[indx] += volume * tmp;
      if ( residue ) element_residue[indx] -= tmp;
      if ( materi_stress ) {
          // stress gradient (rhside with green partial integration)
        tmp = force[inol*ndim+idim];
        element_rhside[indx] -= volume * tmp;
        if ( sri_on ) {
          // Element-consistent momentum feedback (DIAG lot C/D, fix D-c,
          // second half): the reduced 1-point shear internal force
          // -dt*K_shear*v of the current iterate. Together with the
          // zeroed full-rule shear increment of the feedback stress
          // (see the sigvec block above) the momentum right-hand side
          // carries the ELEMENT internal force B^T*sigma_old +
          // dt*K_SRI*v, so the staggered fixed point solves
          // K_SRI*u = P - B^T*sigma_old (the element solution) and the
          // iteration converges in one pass for linear elasticity.
          double v_shear = 0.;
          for ( jnol=0; jnol<nnol; jnol++ ) {
            for ( jdim=0; jdim<ndim; jdim++ ) {
              v_shear += stiffness_shear[(inol*ndim+idim)*nnol*ndim +
                jnol*ndim+jdim] *
                new_dof[jnol*nuknwn + vel_indx + jdim*nder];
            }
          }
          element_rhside[indx] -= dtime * v_shear / npoint * md_factor;
        }
        for ( jdim=0; jdim<ndim; jdim++ ) {
          iuknwn = stres_indx+stress_indx(idim,jdim)*nder;
          if ( residue ) element_residue[indx] += h[inol] *
            new_grad_new_unknowns[jdim*nuknwn+iuknwn];
        }
        for ( jnol=0; jnol<nnol; jnol++ ) {
            // groundflow
          if ( groundflow_pressure ) {
            jpuknwn = pres_indx/nder;
            indxj = jnol*npuknwn + jpuknwn;
            tmp = volume * new_d[idim*nnol+inol] * h[jnol];
            element_matrix[indxi*nnol*npuknwn+indxj] += tmp;
          }
              // temperature stiffness
          if ( condif_temperature ) {
            for ( kdim=0; kdim<ndim; kdim++ ) {
              jpuknwn = temp_indx/nder;
              indxj = jnol*npuknwn + jpuknwn;
              i = stress_indx(idim,idim);
              j = stress_indx(kdim,kdim);
              ind_ddsdde = i*MSTRAIN+j;
              tmp = volume * materi_expansion_linear * new_d[idim*nnol+inol] *
                 -h[jnol] * ddsdde[ind_ddsdde];
              element_matrix[indxi*nnol*npuknwn+indxj] += tmp;
            }
          }
          for ( jdim=0; jdim<ndim; jdim++ ) {
              // stiffness
            indx1 = inol*ndim + idim;
            indx2 = jnol*ndim + jdim;
            jpuknwn = vel_indx/nder + jdim;
            indxj = jnol*npuknwn + jpuknwn;
            tmp = volume * dtime * stiffness[indx1*nnol*ndim+indx2] * md_factor;
            element_matrix[indxi*nnol*npuknwn+indxj] += tmp;
            if ( indxi==indxj ) element_lhside[indx] += fac * tmp;
            // undrained capacity volumetric stiffness: the excessive
            // pressure p = p_old + (div u)/C depends on the displacement
            // of the iterate, so the momentum tangent gains
            // dt/C * (dN_i/dx_idim)*(dN_j/dx_jdim) (measured: the total
            // response is drained stiffness + 1/C on the volumetric
            // part, undrained1/undrained2 of the corpus).
            if ( undr_active ) {
              tmp = volume * dtime * (1./undr_C) *
                new_d[idim*nnol+inol] * new_d[jdim*nnol+jnol] * md_factor;
              element_matrix[indxi*nnol*npuknwn+indxj] += tmp;
              if ( indxi==indxj ) element_lhside[indx] += fac * tmp;
            }
            if ( sri_on ) {
              // SRI: the reduced-integrated shear stiffness (1x1 at the
              // centroid) is added scaled by 1/npoint because materi()
              // is called once per integration point; the npoint calls
              // sum exactly to the full reduced integral.
              tmp = dtime * stiffness_shear[indx1*nnol*ndim+indx2] / npoint
                * md_factor;
              element_matrix[indxi*nnol*npuknwn+indxj] += tmp;
              if ( indxi==indxj ) element_lhside[indx] += fac * tmp;
            }
              // viscosity
            jpuknwn = vel_indx/nder + jdim;
            indxj = jnol*npuknwn + jpuknwn;
            tmp = volume * new_d[jdim*nnol+inol] * viscosity * new_d[idim*nnol+jnol];
            element_matrix[indxi*nnol*npuknwn+indxj] += tmp;
            if ( idim==jdim && inol==jnol ) element_lhside[indx] += fac * tmp;
            jpuknwn = vel_indx/nder + idim;
            indxj = jnol*npuknwn + jpuknwn;
            tmp = volume * new_d[jdim*nnol+inol] * viscosity * new_d[jdim*nnol+jnol];
            element_matrix[indxi*nnol*npuknwn+indxj] += tmp;
            if ( idim==jdim && inol==jnol ) element_lhside[indx] += fac * tmp;
          }
        }
      }
    }
    if ( scalar_dabs(h[inol])>EPS_H ) {

      if ( materi_history_variables ) {
        for ( i=0; i<materi_history_variables; i++ ) {
          ipuknwn = hisv_indx/nder + i;
          indx = inol*npuknwn + ipuknwn;
          tmp = volume * h[inol] * ( new_hisv[i] - old_hisv[i] ) / dtime;
          element_rhside[indx] += tmp;
		//added for options_element_dof
          if(options_element_dof==-YES) new_unknowns[hisv_indx + i*nder] = new_hisv[i];
          ipuknwn++;
        }
      }

      if ( materi_plasti_kappa_shear ) {
        ipuknwn = kapsh_indx/nder;
        indx = inol*npuknwn + ipuknwn;
        tmp = volume * h[inol] * ( new_kapsh - old_kapsh ) / dtime;
        element_rhside[indx] += tmp;
        ipuknwn++;
      }

      if ( materi_damage ) {
        ipuknwn = dam_indx/nder;
        indx = inol*npuknwn + ipuknwn;
        tmp = volume * h[inol] * ( new_damage - old_damage ) / dtime;
        element_rhside[indx] += tmp;
      }
      
      if ( materi_stress ) {
        ipuknwn = stres_indx/nder;
        // stress dof recovery weight (DIAG lot C/D, fix D-b): with the
        // SRI quad4 (2x2 Gauss) the lumped h-weighted average dilutes
        // the corner nodal stresses; the consistent recovery is the
        // Lagrange extrapolation of the Gauss-point values to the nodes
        // (the "same B at the node"), which for every other quadrature
        // (node-containing rules) reduces to h. The momentum feedback
        // does NOT read these dofs (it uses the fresh constitutive
        // stress), so this changes only the OUTPUT stress field used by
        // the section-force integration and the prints. The NORMAL
        // stresses of the bilinear are superconvergent at the Gauss
        // points (the extrapolation is the exact nodal recovery); the
        // SHEAR stress is not (the interpolation error dominates), so
        // its dofs keep the h-weighting (the centroid-biased average,
        // the same estimate the SRI's reduced point provides).
        double weight = h[inol];
        if ( sri_on )
          weight = sri_stress_recovery_weight( nnol, inol, npoint,
            ipoint, h[inol], 1 );
        for ( idim=0; idim<MDIM; idim++ ) {
          for ( jdim=idim; jdim<MDIM; jdim++ ) {
            indx = inol*npuknwn + ipuknwn;
            tmp = volume * ( (idim==jdim) ? weight : h[inol] ) *
              ( new_sig[idim*MDIM+jdim] - old_sig[idim*MDIM+jdim] ) / dtime;
            element_rhside[indx] += tmp;
            if ( jdim==1 && idim==1 ) {
            }
            if ( jdim==0 && idim==0 ) {
            }
		//added for options_element_dof
            // NOTE (measured 2026-08-29, the ELEMENT_DOF zero-stress
            // mystery of the 3D `derivatives` models): the write index
            // must be the SLOT stres_indx + j*nder (the value slot of
            // the j-th stress component), NOT the unknown-number index
            // stres_indx/nder + j. For nder=1 both coincide (the 2D
            // default); for nder>1 (3D with `derivatives`, nder=5) the
            // old index landed the constitutive stress inside the
            // displacement block, so the elem.cc ELEMENT_DOF write
            // (which reads the stres slots) never saw it and the
            // element integration-point stresses stayed zero. Same
            // fix for the materi_history_variables write above.
            if(options_element_dof==-YES) new_unknowns[stres_indx + (ipuknwn - stres_indx/nder)*nder] = new_sig[idim*MDIM+jdim];// new_sig_nonrot[idim*MDIM+jdim];
            ipuknwn++;
          }
        }
      }
      if ( materi_plasti_f ) {
        ipuknwn = f_indx/nder;
        indx = inol*npuknwn + ipuknwn;
        tmp = volume * h[inol] * ( new_f - old_f ) / dtime;
        element_rhside[indx] += tmp;
      }

      if ( materi_plasti_incremental_substeps ) {
        ipuknwn = substeps_indx/nder;
        indx = inol*npuknwn + ipuknwn;
        tmp = volume * h[inol] * ( new_substeps - old_substeps ) / dtime;
        element_rhside[indx] += tmp;
      }

      if ( materi_plasti_kappa ) {
        ipuknwn = kap_indx/nder;
        indx = inol*npuknwn + ipuknwn;
        tmp = volume * h[inol] * ( new_kappa - old_kappa ) / dtime;
        element_rhside[indx] += tmp;
      }

      if ( materi_plasti_cap1_history ) {
        ipuknwn = cap1_indx/nder;
        indx = inol*npuknwn + ipuknwn;
        tmp = volume * h[inol] * ( new_cap1pc - old_cap1pc ) / dtime;
        element_rhside[indx] += tmp;
      }

      if ( materi_plasti_rho ) {
        ipuknwn = rho_indx/nder;
        for ( idim=0; idim<MDIM; idim++ ) {
          for ( jdim=idim; jdim<MDIM; jdim++ ) {
            indx = inol*npuknwn + ipuknwn;
            tmp = volume * h[inol] * ( new_rho[idim*MDIM+jdim] - 
              old_rho[idim*MDIM+jdim] ) / dtime;
            element_rhside[indx] += tmp;
            ipuknwn++;
          }
        }
      }
      if ( materi_plasti_softvar_local ) {
        ipuknwn = svloc_indx/nder;
        indx = inol*npuknwn + ipuknwn;
        tmp = volume * h[inol] * ( new_unknowns[ipuknwn] - old_unknowns[ipuknwn] ) / dtime;
        element_rhside[indx] += tmp;
      }
      if ( materi_plasti_softvar_nonlocal ) {
        ipuknwn = svnonloc_indx/nder;
        indx = inol*npuknwn + ipuknwn;
        tmp = volume * h[inol] * ( new_unknowns[ipuknwn] - old_unknowns[ipuknwn] ) / dtime;
        element_rhside[indx] += tmp;
      }
      if ( materi_strainenergy ) {
        ipuknwn = ener_indx/nder;
        indx = inol*npuknwn + ipuknwn;
        tmp = 0.5 * array_inproduct( new_sig, new_epe, MDIM*MDIM );
        tmp = volume * h[inol] * ( tmp - old_unknowns[ener_indx] ) / dtime;
        element_rhside[indx] += tmp;
      }

      if ( materi_strain_intergranular ) {
        ipuknwn = epi_indx/nder;
        for ( idim=0; idim<MDIM; idim++ ) {
          for ( jdim=idim; jdim<MDIM; jdim++ ) {
            indx = inol*npuknwn + ipuknwn;
            tmp = volume * h[inol] * ( new_epi[idim*MDIM+jdim] - 
              old_epi[idim*MDIM+jdim] ) / dtime;
            element_rhside[indx] += tmp;
		//added for options_element_dof
            if(options_element_dof==-YES) 
	      new_unknowns[ipuknwn] = new_epi[idim*MDIM+jdim];
            ipuknwn++;
          }
        }
      }

      if ( materi_strain_elasti ) {
        ipuknwn = epe_indx/nder;
        for ( idim=0; idim<MDIM; idim++ ) {
          for ( jdim=idim; jdim<MDIM; jdim++ ) {
            indx = inol*npuknwn + ipuknwn;
            tmp = volume * h[inol] * inc_epe[idim*MDIM+jdim] / dtime;
            element_rhside[indx] += tmp;
		//added for options_element_dof
            if(options_element_dof==-YES) 
	      new_unknowns[ipuknwn] = old_epe[idim*MDIM+jdim] + inc_epe[idim*MDIM+jdim];
            ipuknwn++;
          }
        }
      }

      if ( materi_strain_plasti ) {
        ipuknwn = epp_indx/nder;
        for ( idim=0; idim<MDIM; idim++ ) {
          for ( jdim=idim; jdim<MDIM; jdim++ ) {
            indx = inol*npuknwn + ipuknwn;
            tmp = volume * h[inol] * inc_epp[idim*MDIM+jdim] / dtime;
            element_rhside[indx] += tmp;
	      // added for options_element_dof
            if(options_element_dof==-YES) 
	      new_unknowns[ipuknwn] = old_epp[idim*MDIM+jdim] + inc_epp[idim*MDIM+jdim];
            ipuknwn++;
          }
        }
        if ( condif_temperature ) {
          ipuknwn = temp_indx/nder;
          indx = inol*npuknwn + ipuknwn;
          tmp = plasti_heatgeneration * h[inol] *
            array_inproduct( new_sig, inc_epp, MDIM*MDIM ) / dtime;
          element_rhside[indx] += tmp;
          element_residue[indx] -= tmp;
        }
      }

      if ( materi_strain_plasti_hardsoil || materi_strain_plasti_cap ||
           materi_strain_plasti_compression || materi_strain_plasti_diprisco ||
           materi_strain_plasti_druckprag ) {
        // materi_strain_plasti_<model> (manual Professional 4.35-4.40):
        // the plastic strain specifically for each model, added to the
        // node_dof records. All dedicated dofs accumulate the SAME
        // plastic strain increment (inc_epp, computed by the stress
        // law driver) with the SAME integration as materi_strain_plasti
        // -- the per-model initias are registration aliases pointing to
        // this single consolidated mechanism (dedicated dof per model:
        // hsepp/capepp/cepp/depp/dpepp, basenames epphs*/eppcap*/...).
        long int epp_model_active[5], epp_model_indx[5], epp_model_ii=0;
        epp_model_active[0] = materi_strain_plasti_hardsoil;
        epp_model_indx[0] = hsepp_indx;
        epp_model_active[1] = materi_strain_plasti_cap;
        epp_model_indx[1] = capepp_indx;
        epp_model_active[2] = materi_strain_plasti_compression;
        epp_model_indx[2] = cepp_indx;
        epp_model_active[3] = materi_strain_plasti_diprisco;
        epp_model_indx[3] = depp_indx;
        epp_model_active[4] = materi_strain_plasti_druckprag;
        epp_model_indx[4] = dpepp_indx;
        for ( epp_model_ii=0; epp_model_ii<5; epp_model_ii++ ) {
          if ( epp_model_active[epp_model_ii] ) {
            ipuknwn = epp_model_indx[epp_model_ii]/nder;
            for ( idim=0; idim<MDIM; idim++ ) {
              for ( jdim=idim; jdim<MDIM; jdim++ ) {
                indx = inol*npuknwn + ipuknwn;
                tmp = volume * h[inol] * inc_epp[idim*MDIM+jdim] / dtime;
                element_rhside[indx] += tmp;
                // added for options_element_dof
                if(options_element_dof==-YES)
                  new_unknowns[ipuknwn] = old_epp[idim*MDIM+jdim] +
                    inc_epp[idim*MDIM+jdim];
                ipuknwn++;
              }
            }
          }
        }
      }

      if ( condif_temperature ) {
        ipuknwn = temp_indx/nder;
        indx = inol*npuknwn + ipuknwn;
        tmp = viscosity_heatgeneration * h[inol];
        element_rhside[indx] += tmp;
        element_residue[indx] -= tmp;
      }                    

      if ( materi_strain_total ) {
        ipuknwn = ept_indx/nder;
        for ( idim=0; idim<MDIM; idim++ ) {
          for ( jdim=idim; jdim<MDIM; jdim++ ) {
            indx = inol*npuknwn + ipuknwn;
            tmp = volume * h[inol] * inc_ept[idim*MDIM+jdim] / dtime;
            element_rhside[indx] += tmp;
	      //added for options_element_dof
            if(options_element_dof==-YES) 
	      new_unknowns[ipuknwn] = old_ept[idim*MDIM+jdim] + inc_ept[idim*MDIM+jdim];
            ipuknwn++;
          }
        }
      }

      if ( materi_maxwell_stress ) {
        ipuknwn = mstres_indx/nder;
        for ( m=0; m<materi_maxwell_stress; m++ ) {
          for ( idim=0; idim<MDIM; idim++ ) {
            for ( jdim=idim; jdim<MDIM; jdim++ ) {
              indx = inol*npuknwn + ipuknwn;
              tmp = volume * h[inol] * ( new_msig[m*MDIM*MDIM+idim*MDIM+jdim] -
                  old_msig[m*MDIM*MDIM+idim*MDIM+jdim] ) / dtime;
              element_rhside[indx] += tmp;
              ipuknwn++;
            }
          }
        }
      }

      if ( materi_void_fraction && materi_strain_plasti ) {
        ipuknwn = void_indx/nder;
        indx = inol*npuknwn + ipuknwn;
        void_fraction = new_unknowns[void_indx];
        tmp = volume * h[inol] * ( 1 - void_fraction ) * void_fraction * 
          ( inc_epp[0*MDIM+0] + inc_epp[1*MDIM+1] + inc_epp[2*MDIM+2] );
        element_rhside[indx] += tmp;
      }

      if ( materi_work ) {
        ipuknwn = work_indx/nder;
        indx = inol*npuknwn + ipuknwn;
        array_subtract( new_sig, old_sig, work, MDIM*MDIM );
        tmp = array_inproduct( work, inc_ept, MDIM*MDIM );
        tmp = volume * h[inol] * ( tmp - old_unknowns[work_indx] ) / dtime;
        element_rhside[indx] += tmp;
      }

    }
    }//Not used when searching for local values of softening variable -- end
  }

  // mesh_activate_gravity method 2: the element stays active but with a
  // REDUCED STIFFNESS before/during activation (the reduced matrix keeps the
  // not-yet-activated material from artificially stiffening the model).
  if ( activation_stiff!=1. ) {
    long int nmat = nnol*npuknwn*nnol*npuknwn;
    for ( long int im=0; im<nmat; im++ )
      element_matrix[im] *= activation_stiff;
    for ( long int im=0; im<nnol*npuknwn; im++ )
      element_lhside[im] *= activation_stiff;
  }

  if ( swit ) {
    pri( "tendon_element_rhside", tendon_element_rhside, nnol, npuknwn );
    pri( "element_lhside", element_lhside, nnol, npuknwn );
    pri( "element_rhside", element_rhside, nnol, npuknwn );
    if ( residue ) pri( "element_residue", element_residue, nnol, npuknwn );
  }

  delete[] dbl_array;

  if ( swit ) pri( "Out function MATERI" );

}

// strain_settlement_creep - strain_settlement_parameters (Carril B T2).
//
// Extra vertical settlement creep strain (soil dumping). The vertical creep
// strain of a soil particle is assumed to be (saturating power law, decided
// 2026-08-18; the manual OCR is ambiguous):
//
//   eps_zz(t) = Ar * (t/t_tr)^n / (t_plus + (t/t_tr)^n)
//
// with Ar = reference_creep_strain_rate, t_plus = time_plus, t_tr =
// reference_time, n = power_n, and t = the time elapsed after the material
// has become active (dumping). The horizontal creep strains are
// eps_xx = eps_yy = lateral_factor * eps_zz. The creep starts when the
// global time reaches time_global_start. Applied in set_deften_etc by adding
// the creep strain increment to inc_ept.
void strain_settlement_creep( long int element, long int gr, long int nnol,
  double inc_ept[], double dtime )

{
  long int i=0, inod=0, length=0, ldum=0, in_geometry=0, apply=0,
    idum[1], *el=NULL, *nodes=NULL, *gr_list=NULL, time_global_start_idx=0;
  double time_global_start=0., time_plus=0., ar=0., t_ref=0., n=0.,
    lateral=0., t_total=0., t_active=0., t_creep=0., t_creep_old=0.,
    eps=0., eps_old=0., deps=0., ddum[1], *par=NULL, time_current=0.;

  if ( !db_active_index( STRAIN_SETTLEMENT_PARAMETERS, 0, VERSION_NORMAL ) )
    return;
  // element group selection (strain_settlement_element_group or -all)
  length = 0;
  apply = 0;
  if ( db_active_index( STRAIN_SETTLEMENT_ELEMENT_GROUP, 0, VERSION_NORMAL ) ) {
    gr_list = db_int( STRAIN_SETTLEMENT_ELEMENT_GROUP, 0, VERSION_NORMAL );
    length = db_len( STRAIN_SETTLEMENT_ELEMENT_GROUP, 0, VERSION_NORMAL );
    for ( i=0; i<length; i++ )
      if ( gr_list[i]==gr ) { apply = 1; break; }
  }
  else apply = 1;

  if ( !apply ) return;

  // parameters: time_global_start time_plus Ar t_ref n lateral
  par = db_dbl( STRAIN_SETTLEMENT_PARAMETERS, 0, VERSION_NORMAL );
  time_global_start = par[0];
  time_plus = par[1];
  ar = par[2];
  t_ref = par[3];
  n = par[4];
  lateral = par[5];
  if ( ar==0. || t_ref<=0. ) return;

  // strain_settlement_diagram: a parameter (1=time_plus, 2=Ar, 3=t_ref,
  // 4=n, 5=lateral) can depend on a dof value (strain_settlement_diagram_dof)
  // via the strain_settlement_diagram table. The dof is read at the first
  // node of the element.
  if ( db_active_index( STRAIN_SETTLEMENT_DIAGRAM, 0, VERSION_NORMAL ) ) {
    long int idof = 0, number = 0;
    db( STRAIN_SETTLEMENT_DIAGRAM_DOF, 0, &idof, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( STRAIN_SETTLEMENT_DIAGRAM_NUMBER, 0, &number, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    if ( idof!=0 && number>0 ) {
      el = get_new_int(MAXIMUM_NODE+1);
      nodes = get_new_int(MAXIMUM_NODE);
      db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
      long int nnol2 = length - 1;
      array_move( &el[1], nodes, nnol2 );
      long int inod0 = nodes[0];
      double *ndof = db_dbl( NODE_DOF, inod0, VERSION_NORMAL );
      long int ndof_len = db_len( NODE_DOF, inod0, VERSION_NORMAL );
      long int idx = idof;
      if ( idx<0 ) {
        long int *dof_label = get_new_int(MUKNWN);
        db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
        array_member( dof_label, idx, nuknwn, idx );
        if ( ndof_len==npuknwn ) idx /= nder;
        delete[] dof_label;
      }
      if ( idx>=0 && idx<ndof_len ) {
        double val = ndof[idx];
        double *diagram = db_dbl( STRAIN_SETTLEMENT_DIAGRAM, 0, VERSION_NORMAL );
        long int ndiag = db_len( STRAIN_SETTLEMENT_DIAGRAM, 0, VERSION_NORMAL );
        double dval = 0.;
        table_xy( diagram, "STRAIN_SETTLEMENT_DIAGRAM", ndiag, val, dval );
        if      ( number==1 ) time_plus = dval;
        else if ( number==2 ) ar = dval;
        else if ( number==3 ) t_ref = dval;
        else if ( number==4 ) n = dval;
        else if ( number==5 ) lateral = dval;
      }
      delete[] el;
      delete[] nodes;
    }
  }

  db( TIME_CURRENT, 0, idum, &time_current, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  t_total = time_current + dtime;
  if ( t_total < time_global_start ) return;

  // time of activation (dumping): the element activation time from
  // mesh_activate_gravity_time, or the global start if not used
  t_active = time_global_start;
  if ( db_active_index( MESH_ACTIVATE_GRAVITY_TIME, 0, VERSION_NORMAL ) ) {
    // element lowest-coordinate interpolation (reuse the activation factor
    // machinery: recompute the element start time here)
    el = get_new_int(MAXIMUM_NODE+1);
    nodes = get_new_int(MAXIMUM_NODE);
    db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
    nnol = length - 1;
    array_move( &el[1], nodes, nnol );
    double ts=0., te=0.;
    db( MESH_ACTIVATE_GRAVITY_TIME, 0, idum, ddum, ldum, VERSION_NORMAL, GET );
    ts = ddum[0]; te = ddum[1];
    // activation starts at ts for the lowest element; here we use ts as the
    // element activation time (single-element simplification)
    t_active = ts;
    delete[] el;
    delete[] nodes;
  }

  t_creep = t_total - t_active;
  t_creep_old = t_creep - dtime;
  if ( t_creep <= 0. ) return;
  if ( t_creep_old < 0. ) t_creep_old = 0.;

  double tn = pow( t_creep / t_ref, n );
  double tn_old = pow( t_creep_old / t_ref, n );
  eps = ar * tn / ( time_plus + tn );
  eps_old = ar * tn_old / ( time_plus + tn_old );
  deps = eps - eps_old;
  if ( deps==0. ) return;

  // add the creep increment to the vertical component (y in 2D, z in 3D)
  // and the lateral components
  long int vc = ( ndim==3 ) ? 8 : 4;   // vertical index in 3x3 row-major
  inc_ept[vc] += deps;
  if ( ndim>=1 ) inc_ept[0] += lateral * deps;          // xx
  if ( ndim>=2 ) inc_ept[4] += lateral * deps;          // yy (2D: both lateral)
  if ( ndim==3 ) inc_ept[8] += deps;                    // zz (already vertical)
}

void set_deften_etc( long int element, long int gr, long int nnol, double h[], 
  double old_coord[], double old_unknowns[], double new_unknowns[], 
  double old_grad_old_unknowns[], double old_grad_new_unknowns[], 
  double old_deften[], double new_deften[], 
  double old_ept[], double inc_ept[], double new_ept[], 
  double old_rot[], double inc_rot[], double new_rot[] )

{
  long int idim=0, jdim=0, indx=0, swit=0, memory=-UPDATED, 
    axisymmetric=-NO, ind=0, ldum=0, idum[1];
  double dtime=0., radius=0., ddum[1], 
    uTu[MDIM*MDIM], old_u[MDIM*MDIM], new_u[MDIM*MDIM], 
    inc_u[MDIM*MDIM], inc_deften[MDIM*MDIM], coord_ip[MDIM];

  swit = set_swit(element,-1,"set_deften_etc");
  if ( swit ) pri( "In routine SET_DEFTEN_ETC." );

    // initialize
  db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET );
  db( GROUP_MATERI_MEMORY, gr, &memory, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( GROUP_AXISYMMETRIC, gr, &axisymmetric, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  array_set( old_deften, 0., MDIM*MDIM ); 
  array_set( new_deften, 0., MDIM*MDIM ); 
  array_set( inc_deften, 0., MDIM*MDIM ); 
  array_set( old_rot, 0., MDIM*MDIM ); 
  array_set( new_rot, 0., MDIM*MDIM ); 
  array_set( inc_rot, 0., MDIM*MDIM ); 
  array_set( inc_ept, 0., MDIM*MDIM ); 
  array_set( new_ept, 0., MDIM*MDIM ); 

    // deformation tensors
  if ( materi_displacement || materi_velocity_integrated ) {
    old_deften[0] = old_deften[4] = old_deften[8] = 1.;
    new_deften[0] = new_deften[4] = new_deften[8] = 1.;
    for ( idim=0; idim<ndim; idim++ ) {
      for ( jdim=0; jdim<ndim; jdim++ ) {
        indx = idim*MDIM+jdim;
        old_deften[indx] +=
          old_grad_old_unknowns[jdim*nuknwn+dis_indx+idim*nder];
        new_deften[indx] += 
          old_grad_new_unknowns[jdim*nuknwn+dis_indx+idim*nder];
      }
    }
  }
  inc_deften[0] = inc_deften[4] = inc_deften[8] = 1.;
  for ( idim=0; idim<ndim; idim++ ) {
    for ( jdim=0; jdim<ndim; jdim++ ) {
      indx = idim*MDIM+jdim;
      ind = jdim*nuknwn+vel_indx+idim*nder;
      inc_deften[indx] += old_grad_new_unknowns[ind]*dtime;
    }
  }
  if ( axisymmetric==-YES ) {
    matrix_ab( h, old_coord, coord_ip, 1, nnol, ndim );
    radius = scalar_dabs(coord_ip[0]);
    if ( radius!=0. ) {
      if      ( materi_displacement ) {
        ind = dis_indx;
        old_deften[8] += old_unknowns[ind]/radius;
        new_deften[8] += new_unknowns[ind]/radius;
      }
      else if ( materi_velocity_integrated ) {
        ind = veli_indx;
        old_deften[8] += old_unknowns[ind]/radius;
        new_deften[8] += new_unknowns[ind]/radius;
      }
      ind = vel_indx;
      inc_deften[8] += new_unknowns[ind]*dtime/radius;              
    }
  }

    // rotation matrices
  if      ( memory==-UPDATED_WITHOUT_ROTATION || memory==-UPDATED_LINEAR ||
            memory==-UPDATED_AREA || memory==-TOTAL_LINEAR ) {
    for ( idim=0; idim<MDIM; idim++ ) {
      old_rot[idim*MDIM+idim] = 1.;
      new_rot[idim*MDIM+idim] = 1.;
      inc_rot[idim*MDIM+idim] = 1.;
    }
  }
  else if ( materi_displacement || materi_velocity_integrated ) {
    set_deften_u_rot( old_deften, old_u, old_rot );
    set_deften_u_rot( new_deften, new_u, new_rot );
  }
  set_deften_u_rot( inc_deften, inc_u, inc_rot );

      // strain matrices
  if ( memory==-UPDATED_WITHOUT_ROTATION || memory==-UPDATED_LINEAR ||
       memory==-UPDATED_AREA ) {
      // linear engineering strains
    for ( idim=0; idim<MDIM; idim++ ) {
      for ( jdim=0; jdim<MDIM; jdim++ ) inc_ept[idim*MDIM+jdim] = 
        0.5*(inc_deften[idim*MDIM+jdim]+inc_deften[jdim*MDIM+idim]);
    }
    for ( idim=0; idim<MDIM; idim++ ) inc_ept[idim*MDIM+idim] -= 1.;
    array_add( old_ept, inc_ept, new_ept, MDIM*MDIM );
  }
  else if ( memory==-UPDATED ) {
      // U from incremental polar decomposition
    for ( idim=0; idim<MDIM; idim++ ) {
      for ( jdim=0; jdim<MDIM; jdim++ ) inc_ept[idim*MDIM+jdim] = 
        0.5*(inc_u[idim*MDIM+jdim]+inc_u[jdim*MDIM+idim]);
    }
    for ( idim=0; idim<MDIM; idim++ ) inc_ept[idim*MDIM+idim] -= 1.;
    array_add( old_ept, inc_ept, new_ept, MDIM*MDIM );
  }
  else if ( memory==-TOTAL ) {
      // U from total polar decomposition
    for ( idim=0; idim<MDIM; idim++ ) {
      for ( jdim=0; jdim<MDIM; jdim++ ) new_ept[idim*MDIM+jdim] = 
        0.5*(new_u[idim*MDIM+jdim]+new_u[jdim*MDIM+idim]);
    }
    for ( idim=0; idim<MDIM; idim++ ) new_ept[idim*MDIM+idim] -= 1.;
    array_subtract( new_ept, old_ept, inc_ept, MDIM*MDIM );
  }
  else if ( memory==-TOTAL_LINEAR ) {
      // linear engineering strains
    for ( idim=0; idim<MDIM; idim++ ) {
      for ( jdim=0; jdim<MDIM; jdim++ ) new_ept[idim*MDIM+jdim] = 
        0.5*(new_deften[idim*MDIM+jdim]+new_deften[jdim*MDIM+idim]);
    }
    for ( idim=0; idim<MDIM; idim++ ) new_ept[idim*MDIM+idim] -= 1.;
    array_subtract( new_ept, old_ept, inc_ept, MDIM*MDIM );
  }
  else if ( memory==-TOTAL_PIOLA ) {
      // Green-Lagrange strains
    matrix_atb( new_deften, new_deften, uTu, MDIM, MDIM, MDIM );
    for ( idim=0; idim<MDIM; idim++ ) {
      for ( jdim=0; jdim<MDIM; jdim++ ) {
        new_ept[idim*MDIM+jdim] += 0.5*uTu[idim*MDIM+jdim];
        if ( idim==jdim ) new_ept[idim*MDIM+jdim] -= 0.5;
      }
    }
    array_subtract( new_ept, old_ept, inc_ept, MDIM*MDIM );
  }
  else
    db_error( GROUP_MATERI_MEMORY, gr );

  // strain_settlement_parameters: extra vertical settlement creep strain
  // (soil dumping). Adds the creep strain increment to inc_ept and updates
  // new_ept accordingly.
  strain_settlement_creep( element, gr, nnol, inc_ept, dtime );
  array_add( old_ept, inc_ept, new_ept, MDIM*MDIM );

  if ( swit ) {
    pri( "old_deften", old_deften, MDIM, MDIM );
    pri( "new_deften", new_deften, MDIM, MDIM );
    pri( "inc_deften", inc_deften, MDIM, MDIM );
    pri( "old_rot", old_rot, MDIM, MDIM );
    pri( "inc_rot", inc_rot, MDIM, MDIM );
    pri( "new_rot", new_rot, MDIM, MDIM );
    pri( "inc_ept", inc_ept, MDIM, MDIM );
    pri( "new_ept", new_ept, MDIM, MDIM );
  }

  if ( swit ) pri( "Out function SET_DEFTEN_ETC" );
}

void set_deften_u_rot( double deften[], double u[], double rot[] )

{
  long int idim=0, jdim=0, kdim=0, indx=0, idum[1];
  double rdum=0., uTu[MDIM*MDIM], val[MDIM], dir[MDIM*MDIM], work[MDIM*MDIM];

  matrix_atb( deften, deften, uTu, MDIM, MDIM, MDIM );
  matrix_jacobi( uTu, MDIM, val, dir, idum );
  for ( idim=0; idim<MDIM; idim++ ) {
    for ( jdim=0; jdim<MDIM; jdim++ ) {
      indx = idim*MDIM+jdim;
      u[indx] = 0.;
      for ( kdim=0; kdim<MDIM; kdim++ ) u[indx] += 
        sqrt(scalar_dabs(val[kdim]))*dir[idim*MDIM+kdim]*dir[jdim*MDIM+kdim];
    }
  }
  if ( matrix_inverse( u, work, rdum, MDIM ) ) {
    matrix_ab( deften, work, rot, MDIM, MDIM, MDIM );
  }
  else {
    for ( idim=0; idim<MDIM; idim++ ) {
      for ( jdim=0; jdim<MDIM; jdim++ ) {
        u[idim*MDIM+jdim] = 0.5*(deften[idim*MDIM+jdim]+deften[jdim*MDIM+idim]);
        if ( jdim==idim )
          rot[idim*MDIM+idim] = 1.;
        else
          rot[idim*MDIM+idim] = 0.;
      }
    }
  }

}

double get_materi_density( long int element, long int element_group, long int nnol, long int nodes[],
  double new_unknowns[] )

{

  long int ldum=0, all_below=0, node_phreaticlevel=0, inol=0, inod=0,
    idum[1];
  double materi_dens = 0., ddum[1], group_materi_density_groundflow[2];

  if      ( materi_density ) {
    materi_dens = new_unknowns[dens_indx];
    if ( materi_dens<0. ) materi_dens = 0.;
  }             
  else if ( db_active_index( GROUP_MATERI_DENSITY_GROUNDFLOW, element_group, VERSION_NORMAL ) ) {
    db( GROUP_MATERI_DENSITY_GROUNDFLOW, element_group, idum, 
      group_materi_density_groundflow, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    all_below = 1;
    for ( inol=0; inol<nnol; inol++ ) {
      inod = nodes[inol];
      node_phreaticlevel = -BELOW;
      db( NODE_PHREATICLEVEL, inod, &node_phreaticlevel, ddum, 
        ldum, VERSION_NORMAL, GET_IF_EXISTS );
      if ( node_phreaticlevel==-ABOVE ) all_below = 0;
    }
    if ( all_below ) 
      materi_dens = group_materi_density_groundflow[0]; // wet
    else 
      materi_dens = group_materi_density_groundflow[1]; // dry
  }
  else {
    get_group_data( GROUP_MATERI_DENSITY, element_group, element, new_unknowns, 
      &materi_dens, ldum, GET_IF_EXISTS );
  }
  return materi_dens;

}
