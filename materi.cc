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

void materi( long int element, long int gr, long int nnol, 
  long int npoint, long int nodes[], long int plasti_on_boundary,
  double coord_ip[], double old_coord[],
  double h[], double new_d[], double new_b[],
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
    ind_ddsdde=0, ldum=0, idum[1];
  double rdum=0., dens=0., dtime=0., materi_expansion_linear=0., 
    materi_expansion_volume=0., temp=0., tmp=0., damping=0., fac=0, 
    plasti_heatgeneration=0., viscosity_heatgeneration=0.,
    viscosity=0., old_damage=0., new_damage=0., old_kappa=0., new_kappa=0., 
    old_f=0., new_f=0., void_fraction=0., new_pres=0., old_substeps=0., new_substeps=0.,
    softvar_nonl=0, softvar_l=0, 
    static_pressure=0., total_pressure=0., location=0.,
    J=0., ddum[1], *force_gravity=NULL, 
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
    *dbl_array=NULL,
    *new_sig_nonrot=NULL;

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
    nnol*MDIM*nnol*MDIM + nnol*ndim*nnol*ndim + nnol*ndim + MDIM*MDIM;
  dbl_array = get_new_dbl(n);
  assert( indx<=n );
 
  indx = 0;
  work = &dbl_array[indx]; indx += nnol*MDIM*nnol*MDIM;
  stiffness = &dbl_array[indx]; indx += nnol*ndim*nnol*ndim;
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
  db( GROUP_MATERI_MEMORY, gr, &memory, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  dens = get_materi_density( element, gr, nnol, nodes, new_unknowns );
  force_gravity_calculate( force_gravity );

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
      memory==-UPDATED_WITHOUT_ROTATION ) {
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
    set_stress( element, gr, plasti_on_boundary, coord_ip,
      old_unknowns, new_unknowns,
      old_grad_old_unknowns, new_grad_new_unknowns, 
      rotated_old_sig, new_sig, 
      rotated_old_msig, new_msig, inc_ept, new_ept,
      old_epe, inc_epe, old_epp, inc_epp, old_rho, new_rho, 
      old_epi, new_epi, old_hisv, new_hisv, 
      old_damage, new_damage, old_kappa, new_kappa, new_f, new_substeps,
      old_deften, new_deften, inc_rot,
      ddsdde, viscosity, viscosity_heatgeneration, softvar_nonl, softvar_l);
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
  if ( memory==-UPDATED || memory==-TOTAL || memory==-TOTAL_PIOLA ) {
    if      ( memory==-UPDATED ) 
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
      new_pres = new_unknowns[pres_indx];
      if ( groundflow_phreatic_coord( -1, coord_ip, new_unknowns, 
        total_pressure, static_pressure, location ) ) new_pres = total_pressure;
      for ( idim=0; idim<MDIM; idim++ ) total_new_sig[idim*MDIM+idim] += new_pres;
    }
    for ( idim=0; idim<MDIM; idim++ ) {
      for ( jdim=idim; jdim<MDIM; jdim++ ) {
        indx = stress_indx(idim,jdim);
        sigvec[indx] = total_new_sig[idim*MDIM+jdim];
      }
    }
    matrix_atb( new_b, sigvec, force, MSTRAIN, nnol*ndim, 1 );
    matrix_atba( new_b, ddsdde_total, stiffness, work, MSTRAIN, nnol*ndim );
    if ( swit ) {
      pri( "force", force, nnol*ndim );
      pri( "stiffness", stiffness, nnol*ndim, nnol*ndim );
      pri( "sigvec", sigvec, MSTRAIN );
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
            tmp = volume * dtime * stiffness[indx1*nnol*ndim+indx2];
            element_matrix[indxi*nnol*npuknwn+indxj] += tmp;
            if ( indxi==indxj ) element_lhside[indx] += fac * tmp;
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
          if(options_element_dof==-YES) new_unknowns[ipuknwn] = new_hisv[i];
          ipuknwn++;
        }
      }

      if ( materi_damage ) {
        ipuknwn = dam_indx/nder;
        indx = inol*npuknwn + ipuknwn;
        tmp = volume * h[inol] * ( new_damage - old_damage ) / dtime;
        element_rhside[indx] += tmp;
      }
      
      if ( materi_stress ) {
        ipuknwn = stres_indx/nder;
        for ( idim=0; idim<MDIM; idim++ ) {
          for ( jdim=idim; jdim<MDIM; jdim++ ) {
            indx = inol*npuknwn + ipuknwn;
            tmp = volume * h[inol] * ( new_sig[idim*MDIM+jdim] - 
              old_sig[idim*MDIM+jdim] ) / dtime;
            element_rhside[indx] += tmp;
		//added for options_element_dof
            if(options_element_dof==-YES) new_unknowns[ipuknwn] = new_sig[idim*MDIM+jdim];// new_sig_nonrot[idim*MDIM+jdim];
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

  if ( swit ) {
    pri( "tendon_element_rhside", tendon_element_rhside, nnol, npuknwn );
    pri( "element_lhside", element_lhside, nnol, npuknwn );
    pri( "element_rhside", element_rhside, nnol, npuknwn );
    if ( residue ) pri( "element_residue", element_residue, nnol, npuknwn );
  }

  delete[] dbl_array;

  if ( swit ) pri( "Out function MATERI" );

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
  if      ( memory==-UPDATED_WITHOUT_ROTATION || memory==-TOTAL_LINEAR ) {
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
  if ( memory==-UPDATED_WITHOUT_ROTATION ) {
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
