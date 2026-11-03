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

#define EPS_D 1.e-8
#define EPS_LHSIDE 1.e-8
#define EPS_PECLET 1.e-10

// Set only while solve() reassembles the full-stiffness auxiliary system of
// the displacement-primary transient. Element code also recognizes f<1 on
// the normal pass; this marker preserves that reference configuration and
// makes materi_dynamic_effective() return 1 for the auxiliary pass.
long int disp_primary_reference = -NO;

void general( long int element, long int name, long int nnol, long int element_group, 
  long int type, long int nodes[], double coord_ip[], double old_dof[], double new_dof[], 
  double old_unknowns[], double new_unknowns[], 
  double new_grad_new_unknowns[], double h[], double new_d[], 
  double volume, long int npoint, long int ipoint, double grad_massflow[], 
  double element_rhside[], double element_residue[], double element_lhside[], 
  double element_matrix[] )

{
  long int inol=0, inod=0, jnol=0, jdim=0, ipuknwn=0, iuknwn=0,
    swit=0, icontrol=0, unknown_belongs_to_type=0,
    indx=0, indxi=0, indxj=0, options_stabilization=0,
    // inertia_apply (manual Professional 6.785): the inertia terms of the
    // principal dofs (material mass, heat capacity, water storage) are
    // included only when the record is set to -yes. Default (record
    // absent) is -no, matching the Professional ("Default, if
    // inertia_apply is not specified, then each of switch_0, switch_1
    // etc. is set to -no"). The GNU legacy default -yes made every
    // analysis with a material density pseudo-dynamic; the corpus
    // safety/groundflow tests (ground14/15/16) run static in the
    // Professional within the same 1 s window.
    options_inertia=-NO, options_convection=-NO, options_axisymmetric=0,
    lagrange=0, stokes=-YES, principal_unknown=0, ldum=0, 
    idum[1], options_mesh[MDIM], dof_type[MUKNWN], dof_principal[MUKNWN];
  double condif_conductivity=0., dens=0.,
    condif_capacity=0., visc=0., vel=0.,
    C=0., dtime=0., inertia=0., radius=0.,
    conv_part=0., diff_part=0., peclet=0., peclet_factor=0., 
    val=0., val_max=0., val_min=0., val_new=0.,
    artificial_diffusion=0., tmp=0., D=0., diffusion=1.,
    element_lhside_add = 0., element_rhside_add=0.,
    ddum[1], pe[MDIM], node_remesh_velocity[MDIM], condif_flow[MDIM];
  long int element_dof_initial_active=0, edi_idum[1], edi_n=0;
  double edi_vals[DATA_ITEM_SIZE], element_dof_initial_values[MUKNWN];
  long int edis_active=0, edis_label=0, edis_iuknwn=-1, edis_n=0,
    edis_idum[1], dof_label[MUKNWN];
  double edis_vals[1+MDIM], edis_node_values[MNOL];

  if ( type==-MAXWELL_FREQUENCY || type==-MAXWELL_TIME ) return;

  swit = set_swit(element,-1,"general");
  if ( swit ) pri( "In GENERAL" );

  db( DOF_TYPE, 0, dof_type, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( DOF_PRINCIPAL, 0, dof_principal, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET );
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );

  db( OPTIONS_MESH, 0, options_mesh, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  options_convection = convection_apply_switch();
  // Nodal stress recovery (2026-09-15): the absolute weighted average is
  // only valid when the stress dofs are a PURE recovered field. One case
  // keeps the historic incremental lumped update: an INDEPENDENT strain
  // unknown exists (materi_strain_total and the plasti/elasti strain
  // variants): the stress dof evolves with that unknown (measured:
  // examp9/10/23, where the Pro cannot validate the model - unsupported
  // bounda_unknown/options_mesh/contact_geometry - and the absolute form
  // changes the adaptive-remesh path).
  // C72 (F4 T6f, 2026-10-13): a CONVECTED stress dof (convection_apply
  // -yes on a mesh fixed in space) is recovered absolutely too. The
  // Professional 02-08-2026 never transports the nodal stress: it is the
  // recovery of the integration-point stress, which is not convected
  // either (a sigzz pulse stays bit-frozen). The former incremental
  // convected update was also inconsistent: its lumped mass used the
  // per-element recovery weight F/volume while the convective term used
  // h, so the nodal stress moved at about v*volume_ip instead of v.
  stress_recovery_absolute = !materi_strain_total;
  db( OPTIONS_INERTIA, 0, &options_inertia, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  {
    // multi-valued record (Professional: one switch per principal dof);
    // GNU applies the first value to every principal dof
    long int inertia_switches[DATA_ITEM_SIZE];
    if ( db( CONTROL_OPTIONS_INERTIA, icontrol, inertia_switches, ddum, ldum,
        VERSION_NORMAL, GET_IF_EXISTS ) )
      options_inertia = inertia_switches[0];
  }
  db( OPTIONS_STABILIZATION, 0, &options_stabilization, ddum, ldum, VERSION_NORMAL, GET );
  db( GROUP_AXISYMMETRIC, element_group, 
    &options_axisymmetric, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  if ( options_axisymmetric==-YES ) radius = coord_ip[0];

  array_set( pe, 0., ndim );
  array_set( condif_flow, 0., ndim );

  if      ( type==-CONDIF ) {
    if ( materi_density ) {
      dens = new_unknowns[dens_indx];
      if ( dens<0. ) dens = 0.;
    }
    get_group_data( GROUP_CONDIF_DENSITY, element_group, element, new_unknowns, 
      &dens, ldum, GET_IF_EXISTS );
    get_group_data( GROUP_CONDIF_FLOW, element_group, element, new_unknowns, 
      condif_flow, ldum, GET_IF_EXISTS );
    get_group_data( GROUP_CONDIF_CAPACITY, element_group, element, new_unknowns, 
      &condif_capacity, ldum, GET_IF_EXISTS );
    get_group_data( GROUP_CONDIF_CONDUCTIVITY, element_group, element, new_unknowns, 
      &condif_conductivity, ldum, GET_IF_EXISTS );
  }
  else if ( type==-GROUNDFLOW ) {
    groundflow_data( element, element_group, nodes, old_unknowns, new_unknowns, coord_ip, pe, C, h, nnol );
  }
  else if ( type==-MATERI ) {
    dens = get_materi_density( element, element_group, nnol, nodes, new_unknowns );
    // group_materi_expansion_volume (manual Professional 6.665): the
    // density follows rho = rho0*(1 - beta*T) in the momentum INERTIA
    // too, as materi.cc already does for the gravity force (K04
    // incnav3, 2026-10-06: tvely -2.4976 with the mass of rho0, the
    // Professional 02-08-2026 and the test's analytic value -5 with
    // rho = 0.5)
    if ( condif_temperature ) {
      double expansion_volume = 0.;
      get_group_data( GROUP_MATERI_EXPANSION_VOLUME, element_group, element,
        new_unknowns, &expansion_volume, ldum, GET_IF_EXISTS );
      dens *= 1. - expansion_volume * new_unknowns[temp_indx];
    }
    get_group_data( GROUP_MATERI_VISCOSITY, element_group, element, new_unknowns, 
      &visc, ldum, GET_IF_EXISTS );
    // group_materi_stokes (manual Professional 6.752): per element
    // group, default -yes (record absent), i.e. the momentum convection
    // term rho v.grad(v) is dropped unless the group says -no (K21 WU2,
    // 2026-10-12: the Professional 02-08-2026 without the record = -yes,
    // validation_1 velx 0.1465713). DELIBERATE DIVERGENCE (owner
    // 2026-10-12, DIVERGENCES STOKES-IDX, class A): the
    // Professional reads index 0 for every group and ignores the other
    // indices; GNU follows the manual.
    db( GROUP_MATERI_STOKES, element_group, &stokes, ddum, ldum, 
      VERSION_NORMAL, GET_IF_EXISTS );
  }

  // element_dof_initial (manual Professional 6.422): when an element
  // comes the first time to live it assumes it had in the past the dofs
  // of this record; the transient inertia of the birth step then
  // integrates from these initial dofs instead of the step-old ones
  // (the phased-analysis initial field, e.g. an initial temperature).
  // One value per element dof (the same for all element nodes); a
  // single value is used for all dofs. Applied only during the birth
  // step: the marker ELEMENT_DOF_INITIAL_APPLIED is written by
  // step_close (top.cc) once the step converged, so every iteration of
  // the birth step sees the same initial field.
  element_dof_initial_active = 0;
  if ( db( ELEMENT_DOF_INITIAL, element, edi_idum, edi_vals, edi_n,
      VERSION_NORMAL, GET_IF_EXISTS ) && edi_n>0 ) {
    if ( !db_active_index( ELEMENT_DOF_INITIAL_APPLIED, element,
        VERSION_NORMAL ) ) {
      element_dof_initial_active = 1;
      for ( iuknwn=0; iuknwn<nuknwn; iuknwn++ )
        element_dof_initial_values[iuknwn] =
          edi_vals[ ( iuknwn<edi_n ) ? iuknwn : edi_n-1 ];
    }
  }

  // element_dof_initial_specific_number / _value (manual Professional
  // 6.423 / 6.424): initial past value of one dof only, the linear field
  // value_0 + grad . x evaluated at the element nodes (Professional
  // 02-08-2026: bar 0..1 with 0.8 x gives node 2 temp 0.9, 0.1 + 0.8 x
  // gives 0.95). Same birth step and marker as element_dof_initial,
  // which takes precedence when both are given (Professional: 0.7).
  // value_0 alone is a constant field (manual; the Professional gives
  // NaN there, DIVERGENCES DOF-INITIAL-SPECIFIC-CONSTANT).
  edis_active = 0;
  if ( !element_dof_initial_active &&
       db_active_index( ELEMENT_DOF_INITIAL_SPECIFIC_NUMBER, element,
         VERSION_NORMAL ) &&
       !db_active_index( ELEMENT_DOF_INITIAL_APPLIED, element,
         VERSION_NORMAL ) ) {
    if ( !db_active_index( ELEMENT_DOF_INITIAL_SPECIFIC_VALUE, element,
        VERSION_NORMAL ) ) {
      pri( "Error: element_dof_initial_specific_value is missing for element ",
        element );
      db_error( ELEMENT_DOF_INITIAL_SPECIFIC_NUMBER, element );
    }
    db( ELEMENT_DOF_INITIAL_SPECIFIC_NUMBER, element, &edis_label, ddum,
      ldum, VERSION_NORMAL, GET );
    db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET );
    if ( !array_member( dof_label, edis_label, nuknwn, edis_iuknwn ) )
      db_error( ELEMENT_DOF_INITIAL_SPECIFIC_NUMBER, element );
    db( ELEMENT_DOF_INITIAL_SPECIFIC_VALUE, element, edis_idum, edis_vals,
      edis_n, VERSION_NORMAL, GET );
    if ( edis_n!=1 && edis_n<1+ndim )
      db_error( ELEMENT_DOF_INITIAL_SPECIFIC_VALUE, element );
    for ( inol=0; inol<nnol; inol++ ) {
      double *node_coord = db_dbl( NODE, nodes[inol], VERSION_NORMAL );
      edis_node_values[inol] = edis_vals[0];
      if ( edis_n>1 ) {
        for ( jdim=0; jdim<ndim; jdim++ )
          edis_node_values[inol] += edis_vals[1+jdim] * node_coord[jdim];
      }
    }
    edis_active = 1;
  }

  for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
    iuknwn = ipuknwn * nder;

    inertia = 0.;
    unknown_belongs_to_type = 0;
    principal_unknown = 0;
    if      ( type==-CONDIF ) {
      if      ( dof_type[iuknwn]==-CONDIF_TEMPERATURE ) {
        unknown_belongs_to_type = 1;
        inertia = dens * condif_capacity;
        principal_unknown = 1;
      }
    }
    else if ( type==-GROUNDFLOW ) {
      if      ( dof_type[iuknwn]==-GROUNDFLOW_VELOCITY ) {
        unknown_belongs_to_type = 1;
        inertia = 1.;
      }
      else if ( dof_type[iuknwn]==-GROUNDFLOW_PRESSURE ) {
        unknown_belongs_to_type = 1;
        inertia = C;
        principal_unknown = 1;
      }
    }
    else if ( type==-MATERI ) {
      if      ( dof_type[iuknwn]==-MATERI_DAMAGE ) {
        unknown_belongs_to_type = 1;
        inertia = 1.;
      }
      else if ( dof_type[iuknwn]==-MATERI_DENSITY ) {
        unknown_belongs_to_type = 1;
        inertia = 1.;
        if ( options_axisymmetric==-YES ) inertia *= 2. * PIRAD * radius;
        principal_unknown = 1;
      }
      else if ( dof_type[iuknwn]==-MATERI_HISTORY_VARIABLES ) {
        unknown_belongs_to_type = 1;
        inertia = 1.;
      }
      else if ( dof_type[iuknwn]==-MATERI_PLASTI_HYPO_HISTORY ) {
        // manual Professional 4.23: the 8 hypoplasticity history vars
        // (hyhis0..7) are integrated like the generic history variables
        unknown_belongs_to_type = 1;
        inertia = 1.;
      }
      else if ( dof_type[iuknwn]==-MATERI_PLASTI_CAMCLAY_HISTORY ) {
        // manual Professional 4.16: the 2 camclay history vars
        // (cchis0 = void ratio, cchis1 = preconsolidation pressure)
        // are integrated like the generic history variables
        unknown_belongs_to_type = 1;
        inertia = 1.;
      }
      else if ( dof_type[iuknwn]==-MATERI_PLASTI_F ) {
        unknown_belongs_to_type = 1;
        inertia = 1.;
      }
      else if ( dof_type[iuknwn]==-MATERI_PLASTI_INCREMENTAL_SUBSTEPS ) {
        unknown_belongs_to_type = 1;
        inertia = 1.;
      }
      else if ( dof_type[iuknwn]==-MATERI_PLASTI_KAPPA ) {
        unknown_belongs_to_type = 1;
        inertia = 1.;
      }
      else if ( dof_type[iuknwn]==-MATERI_PLASTI_KAPPA_SHEAR ) {
        unknown_belongs_to_type = 1;
        inertia = 1.;
      }
      else if ( dof_type[iuknwn]==-MATERI_PLASTI_CAP1_HISTORY ) {
        unknown_belongs_to_type = 1;
        inertia = 1.;
      }
      else if ( dof_type[iuknwn]==-MATERI_PLASTI_RHO ) {
        unknown_belongs_to_type = 1;
        inertia = 1.;
      }
      else if ( dof_type[iuknwn]==-MATERI_PLASTI_SOFTVAR_LOCAL ) {
        unknown_belongs_to_type = 1;
        inertia = 1.;
      }
      else if ( dof_type[iuknwn]==-MATERI_PLASTI_SOFTVAR_NONLOCAL ) {
        unknown_belongs_to_type = 1;
        inertia = 1.;
      }
      else if ( dof_type[iuknwn]==-MATERI_MAXWELL_STRESS ) {
        unknown_belongs_to_type = 1;
        inertia = 1.;
      }
      else if ( dof_type[iuknwn]==-MATERI_STRAINENERGY ) {
        unknown_belongs_to_type = 1;
        inertia = 1.;
      }
      else if ( dof_type[iuknwn]==-MATERI_STRAIN_ELASTI ) {
        unknown_belongs_to_type = 1;
        inertia = 1.;
      }
      else if ( dof_type[iuknwn]==-MATERI_STRAIN_INTERGRANULAR ) {
        unknown_belongs_to_type = 1;
        inertia = 1.;
      }
      else if ( dof_type[iuknwn]==-MATERI_STRAIN_PLASTI ) {
        unknown_belongs_to_type = 1;
        inertia = 1.;
      }
      else if ( dof_type[iuknwn]==-MATERI_STRAIN_PLASTI_CAP ) {
        unknown_belongs_to_type = 1;
        inertia = 1.;
      }
      else if ( dof_type[iuknwn]==-MATERI_STRAIN_PLASTI_COMPRESSION ) {
        unknown_belongs_to_type = 1;
        inertia = 1.;
      }
      else if ( dof_type[iuknwn]==-MATERI_STRAIN_PLASTI_DIPRISCO ) {
        unknown_belongs_to_type = 1;
        inertia = 1.;
      }
      else if ( dof_type[iuknwn]==-MATERI_STRAIN_PLASTI_DRUCKPRAG ) {
        unknown_belongs_to_type = 1;
        inertia = 1.;
      }
      else if ( dof_type[iuknwn]==-MATERI_STRAIN_PLASTI_HARDSOIL ) {
        unknown_belongs_to_type = 1;
        inertia = 1.;
      }
      else if ( dof_type[iuknwn]==-MATERI_STRAIN_TOTAL ) {
        unknown_belongs_to_type = 1;
        inertia = 1.;
      }
      else if ( dof_type[iuknwn]==-MATERI_STRESS ) {
        unknown_belongs_to_type = 1;
        inertia = 1.;
      }
      else if ( dof_type[iuknwn]==-MATERI_VELOCITY ) {
        unknown_belongs_to_type = 1;
        inertia = dens;
        principal_unknown = 1;
      }
      else if ( dof_type[iuknwn]==-MATERI_VOID_FRACTION ) {
        unknown_belongs_to_type = 1;
        inertia = 1.;
      }
      else if ( dof_type[iuknwn]==-MATERI_WORK ) {
        unknown_belongs_to_type = 1;
        inertia = 1.;
      }
    }
    else if ( type==-WAVE ) {
      if ( dof_type[iuknwn]==-WAVE_FSCALAR ) {
        unknown_belongs_to_type = 1;
        inertia = 1.;
        principal_unknown = 1;
      }
    }

    if ( unknown_belongs_to_type ) {

      for ( inol=0; inol<nnol; inol++ ) {

        indx = inol*npuknwn + ipuknwn;
        inod = nodes[inol];
        element_lhside_add = element_rhside_add = 0.;

        if ( db_active_index( NODE_REMESH_VELOCITY, inod, VERSION_NORMAL ) )
          db( NODE_REMESH_VELOCITY, inod, idum, node_remesh_velocity, ldum, 
            VERSION_NORMAL, GET );
        else
          array_set( node_remesh_velocity, 0., ndim );

          // inertia (lumped)
          //
          // Nodal stress recovery (2026-09-15, Professional alignment):
          // the MATERI_STRESS dofs are handled ENTIRELY by materi.cc as an
          // absolute weighted average
          // (sum(volume*weight*new_sig) / sum(volume*weight)), so the
          // historic incremental inertia term
          // (-volume*weight*(new_dof-past_dof)/dtime) is skipped here. The
          // 2026-09-13 least-squares weight (sri_stress_recovery_weight,
          // F/volume) is still used by materi.cc for the numerator and the
          // denominator, so the node remains the per-element arithmetic
          // mean. Every other lumped dof keeps the incremental form.
        if ( ( options_inertia==-YES || !principal_unknown ) &&
             !( dof_type[iuknwn]==-MATERI_STRESS &&
                stress_recovery_absolute ) ) {
          double weight = h[inol];
          // C31 (F4 T4, 2026-10-09): the tet10 5-point rule has h < 0 at
          // the corners (consistent corner share -V/20), so the row-sum
          // lumped inertia/capacity would be negative there (unstable
          // dynamics). The principal dofs (velocity, temperature,
          // pressure, ...) use the HRZ masses instead: the diagonal of the
          // consistent mass scaled to the element total, corner V/36 and
          // mid-edge 4V/27 (sum over the points of volume*weight = V).
          // DELIBERATE DIVERGENCE (project owner 2026-10-09): the
          // Professional keeps the negative row sum. The recovered
          // (non-principal) dofs keep h, consistent with their materi.cc
          // numerators.
          if ( name==-TET10 && npoint==5 && principal_unknown )
            weight = ( inol==0 || inol==2 || inol==5 || inol==9 ) ?
              1./36. : 4./27.;
          // C66 phase 2 (F4 T6f, 2026-10-16, owner decision; DELIBERATE
          // DIVERGENCE TRIA6-HRZ): the tria6 lumps EVERY h-weighted dof
          // (principal dofs and the h-weighted state dofs, whose materi.cc
          // numerators use the same weight) with HRZ masses (pol():
          // plane corner A/19, mid-side 16A/57; axisymmetric r-weighted
          // diagonal scaled to 2 pi A r_c). The Professional row sum gives
          // zero (plane) / negative (axisymmetric) corner masses. The
          // recovered dofs below override the weight with their fit. The
          // groundflow velocity numerator (groundfl.cc) and the condif
          // absorption (condif.cc) use the same weight (C85, 2026-10-17).
          if ( lumped_hrz_active )
            weight = lumped_hrz_weight[inol];
          // K11 (2026-10-06): the strain tensors use the recovery weight
          // of materi.cc (strain_recovery_weight), so numerator and
          // denominator of the lumped update match; kappa too (F4 T6f
          // WU2, 2026-10-11)
          if ( dof_type[iuknwn]==-MATERI_STRAIN_TOTAL ||
               dof_type[iuknwn]==-MATERI_STRAIN_ELASTI ||
               dof_type[iuknwn]==-MATERI_STRAIN_PLASTI ||
               dof_type[iuknwn]==-MATERI_PLASTI_KAPPA ||
               // C68 (F4 T6f): the other integration-point states
               dof_type[iuknwn]==-MATERI_PLASTI_KAPPA_SHEAR ||
               dof_type[iuknwn]==-MATERI_PLASTI_CAP1_HISTORY ||
               dof_type[iuknwn]==-MATERI_DAMAGE ||
               dof_type[iuknwn]==-MATERI_PLASTI_RHO ||
               dof_type[iuknwn]==-MATERI_MAXWELL_STRESS )
            weight = strain_recovery_weight(
              sri_active( element, element_group, name, nnol ), name, nnol,
              inol, npoint, ipoint, h[inol], volume, element_group,
              nodes );
          if ( dof_type[iuknwn]==-MATERI_STRESS ) {
            // Incremental stress dof (stress_recovery_absolute==0, an
            // independent strain unknown exists): keep the 2026-09-13
            // least-squares recovery weight for the incremental lumped
            // update.
            long int kcomp = iuknwn - stres_indx;
            if ( sri_active( element, element_group, name, nnol ) ) {
              if ( kcomp==0 || kcomp==3*nder || kcomp==5*nder )
                weight = sri_stress_recovery_weight( 1, name, nnol, inol,
                  npoint, ipoint, h[inol], volume, element_group,
                  nodes );
            }
            else
              weight = sri_stress_recovery_weight( 0, name, nnol, inol,
                npoint, ipoint, h[inol], volume, element_group,
                nodes );
          }
          double past_dof = old_dof[inol*nuknwn+ipuknwn*nder];
          if ( element_dof_initial_active && ipuknwn<nuknwn )
            past_dof = element_dof_initial_values[ipuknwn];
          else if ( edis_active && iuknwn==edis_iuknwn )
            past_dof = edis_node_values[inol];
          tmp = weight * inertia *
            ( new_dof[inol*nuknwn+ipuknwn*nder] - past_dof ) / dtime;
          element_rhside_add -= volume * tmp;
          element_lhside_add += volume * weight * inertia / dtime;
          if ( residue && dof_principal[iuknwn]>=0 ) element_residue[indx] -= tmp;
          element_matrix[indx*nnol*npuknwn+indx] += volume * weight * inertia / dtime;
        }
          
        for ( jdim=0; jdim<ndim; jdim++ ) {

            // C72: an absolutely recovered stress dof is a pure recovery
            // of the integration-point stress (materi.cc assembles the
            // numerator in the same element_rhside slot); no convective or
            // artificial-diffusion terms may be added to it.
          if ( dof_type[iuknwn]==-MATERI_STRESS && stress_recovery_absolute )
            continue;

            // diffusion: real dissusion
            // diff_part: diffusion term used for determination peclet number
            // conv_part: real convection
          diffusion = conv_part = diff_part = 0.;
          if      ( dof_type[iuknwn]==-CONDIF_TEMPERATURE && type==-CONDIF ) {
            conv_part = dens*condif_capacity;
            diff_part = condif_conductivity;
            diffusion = condif_conductivity;
          }
          else if ( dof_type[iuknwn]==-GROUNDFLOW_PRESSURE && type==-GROUNDFLOW ) {
            conv_part = C;
            diff_part = pe[jdim];
            diffusion = pe[jdim];
          }
          else if ( dof_type[iuknwn]==-MATERI_DENSITY && type==-MATERI ) {
            conv_part = 1.;
            if ( options_axisymmetric==-YES ) conv_part *= 2. * PIRAD * radius;
          }
          else if ( dof_type[iuknwn]==-MATERI_HISTORY_VARIABLES && type==-MATERI )
            conv_part = 1.;
          else if ( dof_type[iuknwn]==-MATERI_PLASTI_HYPO_HISTORY && type==-MATERI )
            conv_part = 1.;
          else if ( dof_type[iuknwn]==-MATERI_PLASTI_CAMCLAY_HISTORY && type==-MATERI )
            conv_part = 1.;
          else if ( dof_type[iuknwn]==-MATERI_MAXWELL_STRESS && type==-MATERI )
            conv_part = 1.;
          else if ( dof_type[iuknwn]==-MATERI_PLASTI_KAPPA && type==-MATERI )
            conv_part = 1.;
          else if ( dof_type[iuknwn]==-MATERI_PLASTI_KAPPA_SHEAR && type==-MATERI )
            conv_part = 1.;
          else if ( dof_type[iuknwn]==-MATERI_PLASTI_CAP1_HISTORY && type==-MATERI )
            conv_part = 1.;
          else if ( dof_type[iuknwn]==-MATERI_PLASTI_RHO && type==-MATERI )
            conv_part = 1.;
          else if ( dof_type[iuknwn]==-MATERI_STRAIN_ELASTI && type==-MATERI )
            conv_part = 1.;
          else if ( dof_type[iuknwn]==-MATERI_STRAIN_INTERGRANULAR && type==-MATERI )
            conv_part = 1.;
          else if ( dof_type[iuknwn]==-MATERI_STRAIN_PLASTI && type==-MATERI )
            conv_part = 1.;
          else if ( dof_type[iuknwn]==-MATERI_STRAIN_PLASTI_CAP && type==-MATERI )
            conv_part = 1.;
          else if ( dof_type[iuknwn]==-MATERI_STRAIN_PLASTI_COMPRESSION && type==-MATERI )
            conv_part = 1.;
          else if ( dof_type[iuknwn]==-MATERI_STRAIN_PLASTI_DIPRISCO && type==-MATERI )
            conv_part = 1.;
          else if ( dof_type[iuknwn]==-MATERI_STRAIN_PLASTI_DRUCKPRAG && type==-MATERI )
            conv_part = 1.;
          else if ( dof_type[iuknwn]==-MATERI_STRAIN_PLASTI_HARDSOIL && type==-MATERI )
            conv_part = 1.;
          else if ( dof_type[iuknwn]==-MATERI_STRAIN_TOTAL && type==-MATERI )
            conv_part = 1.;
          else if ( dof_type[iuknwn]==-MATERI_STRESS && type==-MATERI )
            conv_part = 1.;
          else if ( dof_type[iuknwn]==-MATERI_DAMAGE && type==-MATERI )
            conv_part = 1.;
          else if ( dof_type[iuknwn]==-MATERI_VELOCITY && type==-MATERI ) {
            if ( stokes!=-YES ) conv_part = dens;
            diff_part = visc;
          }
          else if ( dof_type[iuknwn]==-MATERI_VOID_FRACTION && type==-MATERI )
            conv_part = 1.;

          if ( options_convection==-NO || materi_displacement ||
               ( materi_velocity && options_mesh[jdim]==-FOLLOW_MATERIAL ) )
             lagrange = 1;
          else
             lagrange = 0;

// convection velocity with respect to mesh-Need to check logic below.
          // K03: convection_apply only governs the convection of the
          // MATERIAL with respect to the mesh (manual 6.395). The heat
          // transport by group_condif_flow and by the groundflow velocity
          // is applied whatever convection_apply is (measured on the
          // Professional 02-08-2026: sfnet condif1 / condif5 give the
          // convected temp = 1 with convection_apply -no, -yes or absent).
          long int lagrange_flow = ( materi_displacement ||
            ( materi_velocity && options_mesh[jdim]==-FOLLOW_MATERIAL ) );
          vel = 0.;
          if ( !lagrange && materi_velocity )
            vel += new_dof[inol*nuknwn+vel_indx+jdim*nder];
          if ( !lagrange_flow ) {
            if ( dof_type[iuknwn]==-CONDIF_TEMPERATURE ) {
              		vel += condif_flow[jdim];
              		if ( groundflow_velocity )
                	vel += new_dof[inol*nuknwn+gvel_indx+jdim*nder];
            											 }
          					 }
// Need to check lines above and logic- Need to consider option  condif_temperature && groundflow

          vel -= node_remesh_velocity[jdim];

          if      ( dof_type[iuknwn]==-MATERI_DENSITY ) {
            tmp = h[inol] * grad_massflow[jdim*ndim+jdim];
            element_rhside_add -= volume * tmp;
          }
          else if ( (name!=-TRIA3 && name!=-TET4) ||
              scalar_dabs(new_d[jdim*nnol+inol])>EPS_D ) {
              // convective terms
            tmp = h[inol] * conv_part * vel *
               new_grad_new_unknowns[jdim*nuknwn+iuknwn];
            element_rhside_add -= volume * tmp;
            if ( residue && dof_principal[iuknwn]>=0 ) element_residue[indx] -= tmp;
            indxi = inol*npuknwn+ipuknwn;
            for ( jnol=0; jnol<nnol; jnol++ ) {
              indxj = jnol*npuknwn+ipuknwn;
              element_matrix[indxi*nnol*npuknwn+indxj] += 
                volume * h[inol] * vel * conv_part * new_d[jdim*nnol+jnol];
            }
          }

            // artificial diffusion and physical diffusion
          if ( scalar_dabs(new_d[jdim*nnol+inol])>EPS_D ) {
            if      ( options_stabilization==-NO )
              peclet_factor = 0.;
            else {
              if (scalar_dabs(conv_part*h[inol]*vel)<EPS_PECLET )
                peclet_factor = 0.;
              else if ( 10.*scalar_dabs(diff_part*new_d[jdim*nnol+inol])<=
                scalar_dabs(conv_part*h[inol]*vel) ) peclet_factor = 1.;
              else {
                peclet= scalar_dabs(conv_part*h[inol]*vel/
                  (diff_part*new_d[jdim*nnol+inol]));
                peclet_factor = 
                  ((cosh(scalar_dabs(peclet)))/(sinh(scalar_dabs(peclet)))-
                  1/(scalar_dabs(peclet)));
              }
            }
            artificial_diffusion = peclet_factor *
              scalar_dabs(h[inol]*vel*conv_part/new_d[jdim*nnol+inol]+diff_part);

              // diffusive terms (with green partial integration)
            D = diffusion + artificial_diffusion;
            tmp = new_d[jdim*nnol+inol] * D *
              new_grad_new_unknowns[jdim*nuknwn+iuknwn];
            element_rhside_add -= volume * tmp;
            element_lhside_add += volume * new_d[jdim*nnol+inol] * 
              D * new_d[jdim*nnol+inol];
            if ( residue && dof_principal[iuknwn]>=0 ) element_residue[indx] += h[inol] *
              diffusion * new_grad_new_unknowns[jdim*nuknwn+iuknwn+jdim+1];
            indxi = inol*npuknwn+ipuknwn;
            for ( jnol=0; jnol<nnol; jnol++ ) {
              indxj = jnol*npuknwn+ipuknwn;
              element_matrix[indxi*nnol*npuknwn+indxj] += 
                volume * new_d[jdim*nnol+inol] * D * new_d[jdim*nnol+jnol];
            }           
          }
        }
        if ( options_stabilization==-MAXIMAL &&
              ( dof_type[iuknwn]==-CONDIF_TEMPERATURE && type==-CONDIF ) &&
            scalar_dabs(element_lhside_add)>EPS_LHSIDE ) {
          val_max = -1.e10;
          val_min = +1.e10;
          for ( jnol=0; jnol<nnol; jnol++ ) {
            val = new_dof[jnol*nuknwn+iuknwn];
            if ( val>val_max ) val_max = val;
            if ( val<val_min ) val_min = val;
          } 
          val = new_dof[inol*nuknwn+iuknwn];
          val_new = val + element_rhside_add/element_lhside_add;
          if ( val_new>val_max ) val_new = val_max;
          if ( val_new<val_min ) val_new = val_min;
          element_rhside_add = ( val_new - val ) * element_lhside_add;
        }
        element_lhside[indx] += element_lhside_add;
        element_rhside[indx] += element_rhside_add;
      }
    }
  }

  if ( swit ) {
    pri( "element_rhside", element_rhside, nnol, npuknwn );
    pri( "element_residue", element_residue, nnol, npuknwn );
    pri( "element_lhside", element_lhside, nnol, npuknwn );
  }

  if ( swit ) pri( "Out function GENERAL" );

}

// control_materi_*_apply switch helper: returns 1 when the per-timestep
// switch record (indexed by the current ICONTROL) is -no, i.e. the
// feature should be IGNORED for these timesteps (manual Professional
// 6.141-6.154). Default (record absent at the current index): 0
// (feature active).
long int control_materi_gate_off( long int control_item )
{
  long int icontrol=0, swit=-YES, ldum=0, max_index=-1;
  double ddum[1];
  db_max_index( control_item, max_index, VERSION_NORMAL, GET );
  if ( max_index>=0 ) {
    db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( control_item, icontrol, &swit, ddum, ldum, VERSION_NORMAL,
      GET_IF_EXISTS );
    if ( swit==-NO ) return 1;
  }
  return 0;
}

// materi_*_apply switch helper with a global (no index) record and its
// per-timestep control_* record, e.g. materi_damage_apply (manual
// Professional 6.798) and control_materi_damage_apply (6.140). Returns 1
// when the feature should be IGNORED for the current timestep.
// Precedence (measured on the Professional 02-08-2026): the control_*
// record with the index of the current timestep (ICONTROL) overrules the
// global record; default (no record): feature active.
long int materi_apply_gate_off( long int global_item, long int control_item )
{
  long int icontrol=0, swit=-YES, ldum=0;
  double ddum[1];
  db( global_item, 0, &swit, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( control_item, icontrol, &swit, ddum, ldum, VERSION_NORMAL,
    GET_IF_EXISTS );
  return swit==-NO;
}

// Effective materi_dynamic factor (C11). control_materi_dynamic with the
// index of the current timestep (ICONTROL) overrules the global
// materi_dynamic (manual 6.141: same as materi_dynamic, but only for the
// timesteps of that control index); default 1 (no record). While solve()
// reassembles the f=1 auxiliary system of the displacement-primary scheme
// (disp_primary_reference == -YES) the factor is 1. If specified is not
// NULL it returns 1 when one of the two records applies, else 0.
double materi_dynamic_effective( long int *specified )
{
  long int icontrol=0, ldum=0, idum[1], spec=0;
  double f=1., ddum[1];
  if ( db( MATERI_DYNAMIC, 0, idum, &f, ldum, VERSION_NORMAL,
      GET_IF_EXISTS ) ) spec = 1;
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  if ( db( CONTROL_MATERI_DYNAMIC, icontrol, idum, &f, ldum,
      VERSION_NORMAL, GET_IF_EXISTS ) ) spec = 1;
  if ( specified!=NULL ) *specified = spec;
  if ( disp_primary_reference==-YES ) return 1.;
  return f;
}

// convection_apply switch (K03). Single source of truth for every
// OPTIONS_CONVECTION reader (general, beam, truss, spring, conspr):
// convection_apply (manual Professional 6.395, GNU options_convection)
// overruled by control_convection_apply (6.113, GNU
// control_options_convection) with the index of the current timestep
// (ICONTROL). Default, if neither record is set: -no (manual 6.395
// "Default switch is set to -no"; measured on the Professional
// 02-08-2026: condif2 T = 0.5 without the record, 1 with -yes). The GNU
// legacy default was -yes.
long int convection_apply_switch( void )
{
  long int icontrol=0, swit=-NO, ldum=0;
  double ddum[1];
  db( OPTIONS_CONVECTION, 0, &swit, ddum, ldum, VERSION_NORMAL,
    GET_IF_EXISTS );
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_OPTIONS_CONVECTION, icontrol, &swit, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  return swit;
}

// Default group_materi_memory when the record is absent (C16, C18). Single
// source of truth for every GROUP_MATERI_MEMORY reader (materi, stress,
// elem, tendon). Manual Professional 6.797: the default is -total_linear
// when materi_displacement is initialised, else -updated_linear (linear in
// simple shear, measured on Professional 02-08-2026).
// control_materi_updated_apply -yes with the index of the current timestep
// (ICONTROL) makes an unspecified memory -updated (manual Professional
// 6.152) and takes precedence. An explicit group_materi_memory always
// overrides this default.
long int materi_memory_default( void )
{
  long int icontrol=0, swit=0, ldum=0, max_index=-1;
  double ddum[1];
  db_max_index( CONTROL_MATERI_UPDATED_APPLY, max_index, VERSION_NORMAL, GET );
  if ( max_index>=0 ) {
    db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( CONTROL_MATERI_UPDATED_APPLY, icontrol, &swit, ddum, ldum,
      VERSION_NORMAL, GET_IF_EXISTS );
    if ( swit==-YES ) return -UPDATED;
  }
  if ( materi_displacement ) return -TOTAL_LINEAR;
  return -UPDATED_LINEAR;
}

// ---------------------------------------------------------------------
// groundflow_apply (manual Professional 6.555) and linear_calculation_apply
// (manual Professional 6.797): database rewrites at the start of top().

static long int groundflow_apply_skip = 0;

// 1 when groundflow_apply -no was applied: the groundflow element
// equation (elem.cc) is skipped for all timesteps.
long int groundflow_apply_off( void )
{
  return groundflow_apply_skip;
}

static long int general_is_groundflow_dof( long int dof_type_value )
{
  return dof_type_value==-GROUNDFLOW_PRESSURE ||
    dof_type_value==-GROUNDFLOW_VELOCITY ||
    dof_type_value==-GROUNDFLOW_SATURATION ||
    dof_type_value==-GROUNDFLOW_PRESSURE_GRADIENT;
}

// Deletes every record whose name starts with prefix (keep: one item
// that is not deleted, -1 for none).
static void general_delete_prefix( const char *prefix, long int keep )
{
  long int idat=0, max=-1, n=strlen(prefix);
  for ( idat=0; idat<MDAT; idat++ ) {
    if ( idat==keep || strncmp( db_name(idat), prefix, n ) ) continue;
    if ( db_max_index( idat, max, VERSION_NORMAL, GET )>=0 )
      db_delete( idat, VERSION_NORMAL );
  }
}

// Removes the groundflow dof labels from the dof lists of the records
// of item idat (bounda_dof, bounda_unknown, control_reset_dof); a record
// that keeps no dof label is deleted, with the record of item
// companion (-1 none) at the same index.
static void general_remove_groundflow_labels( long int idat,
  long int companion, long int dof_label[], long int dof_type[] )
{
  long int index=0, max=-1, i=0, n=0, k=0, len=0, ndof=0, ndof_kept=0,
    iuknwn=0, *val=NULL, *kept=NULL;
  double ddum[1];
  db_max_index( idat, max, VERSION_NORMAL, GET );
  if ( max<0 ) return;
  kept = get_new_int( db_data_length(idat) );
  for ( index=0; index<=max; index++ ) {
    if ( !db_active_index( idat, index, VERSION_NORMAL ) ) continue;
    val = db_int( idat, index, VERSION_NORMAL );
    len = db_len( idat, index, VERSION_NORMAL );
    ndof = ndof_kept = n = 0;
    for ( i=0; i<len; i++ ) {
      iuknwn = -1;
      if ( val[i]<0 && array_member( dof_label, val[i], nuknwn, iuknwn ) ) {
        ndof++;
        if ( general_is_groundflow_dof( dof_type[iuknwn] ) ) continue;
        ndof_kept++;
      }
      kept[n++] = val[i];
    }
    if ( ndof==ndof_kept ) continue;
    if ( ndof_kept==0 ) {
      db_delete_index( idat, index, VERSION_NORMAL );
      if ( companion>=0 ) db_delete_index( companion, index, VERSION_NORMAL );
    }
    else {
      k = n;
      db( idat, index, kept, ddum, k, VERSION_NORMAL, PUT );
    }
  }
  delete[] kept;
}

// groundflow_apply -no: "the groundflow equation is skipped, and all
// groundflow data is ignored. This is done for all timesteps." Measured
// on the Professional 02-08-2026: the groundflow dofs end at 0 at every
// node (initial node_dof values, prescribed pressures and pressure
// resets ignored) and the solid sees no pore pressure (drained answer).
// GNU: the groundflow dofs are zeroed, removed from the bounda_dof /
// bounda_unknown / control_reset_dof dof lists, the groundflow_*,
// group_groundflow_*, group_interface_groundflow_* and
// control_groundflow_* records are deleted, and elem.cc skips the
// groundflow element equation.
void groundflow_apply_transform( void )
{
  long int swit=-YES, ldum=0, inod=0, max=-1, i=0, length=0, iuknwn=0,
    dof_label[MUKNWN], dof_type[MUKNWN];
  double ddum[1], *node_dof=NULL;

  db( GROUNDFLOW_APPLY, 0, &swit, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  if ( swit!=-NO ) return;
  groundflow_apply_skip = 1;
  if ( nuknwn>0 ) {
    db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET );
    db( DOF_TYPE, 0, dof_type, ddum, ldum, VERSION_NORMAL, GET );
    db_max_index( NODE_DOF, max, VERSION_NORMAL, GET );
    for ( inod=0; inod<=max; inod++ ) {
      if ( !db_active_index( NODE_DOF, inod, VERSION_NORMAL ) ) continue;
      node_dof = db_dbl( NODE_DOF, inod, VERSION_NORMAL );
      length = db_len( NODE_DOF, inod, VERSION_NORMAL );
      for ( i=0; i<length; i++ ) {
        iuknwn = ( length==npuknwn ) ? i*nder : i;
        if ( iuknwn<nuknwn && general_is_groundflow_dof( dof_type[iuknwn] ) )
          node_dof[i] = 0.;
      }
    }
    general_remove_groundflow_labels( BOUNDA_DOF, BOUNDA_TIME,
      dof_label, dof_type );
    general_remove_groundflow_labels( BOUNDA_UNKNOWN, BOUNDA_TIME,
      dof_label, dof_type );
    general_remove_groundflow_labels( CONTROL_RESET_DOF, -1,
      dof_label, dof_type );
  }
  general_delete_prefix( "groundflow_", GROUNDFLOW_APPLY );
  general_delete_prefix( "group_groundflow_", -1 );
  general_delete_prefix( "group_interface_groundflow_", -1 );
  general_delete_prefix( "control_groundflow_", -1 );
}

// linear_calculation_apply -yes (manual Professional 6.797): rewrites
// the database to a linear calculation. Every bullet of the manual whose
// records exist in the GNU is applied; group_materi_elasti_polynomial,
// group_groundflow_nonsaturated, group_interface_materi_hardening and
// group_interface_materi_elasti_stiffness_tangential_diagram do not
// exist in the GNU. Where the manual is silent the Professional
// 02-08-2026 decides: it also deletes the group_materi_plasti_* and
// group_interface_materi_plasti_* records and adds contact_apply -no.
void linear_calculation_apply_transform( void )
{
  long int swit=-NO, ldum=0, index=0, max=-1, idat=0, n=0, len=0,
    found=0, iuknwn=0, memory=0, one=1, no=-NO, idum[1],
    options_mesh[MDIM], dof_label[MUKNWN], *val=NULL;
  double ddum[1], young=0., poisson=0.2, k=0., *dval=NULL, *table=NULL;
  const char *item_name=NULL;

  db( LINEAR_CALCULATION_APPLY, 0, &swit, ddum, ldum, VERSION_NORMAL,
    GET_IF_EXISTS );
  if ( swit!=-YES ) return;

    // substitutions by linear elasticity, per element group
  db_max_index( GROUP_MATERI_ELASTI_HARDSOIL, max, VERSION_NORMAL, GET );
  for ( index=0; index<=max; index++ ) {
    if ( !db_active_index( GROUP_MATERI_ELASTI_HARDSOIL, index,
      VERSION_NORMAL ) ) continue;
    young = db_dbl( GROUP_MATERI_ELASTI_HARDSOIL, index, VERSION_NORMAL )[0];
    db( GROUP_MATERI_ELASTI_YOUNG, index, idum, &young, one,
      VERSION_NORMAL, PUT );
  }
  db_delete( GROUP_MATERI_ELASTI_HARDSOIL, VERSION_NORMAL );
  db_max_index( GROUP_MATERI_ELASTI_YOUNG_POWER, max, VERSION_NORMAL, GET );
  for ( index=0; index<=max; index++ ) {
    if ( !db_active_index( GROUP_MATERI_ELASTI_YOUNG_POWER, index,
      VERSION_NORMAL ) ) continue;
    young = db_dbl( GROUP_MATERI_ELASTI_YOUNG_POWER, index,
      VERSION_NORMAL )[0];
    db( GROUP_MATERI_ELASTI_YOUNG, index, idum, &young, one,
      VERSION_NORMAL, PUT );
  }
  db_delete( GROUP_MATERI_ELASTI_YOUNG_POWER, VERSION_NORMAL );
  db_max_index( GROUP_MATERI_PLASTI_HYPO_WOLFERSDORFF, max, VERSION_NORMAL,
    GET );
  for ( index=0; index<=max; index++ ) {
    if ( !db_active_index( GROUP_MATERI_PLASTI_HYPO_WOLFERSDORFF, index,
      VERSION_NORMAL ) ) continue;
      // hs (second value) as the manual says; the Professional 02-08-2026
      // takes the first value (phi), DIVERGENCES LINCALC-WOLFF
    young = db_dbl( GROUP_MATERI_PLASTI_HYPO_WOLFERSDORFF, index,
      VERSION_NORMAL )[1];
    poisson = 0.2;
    db( GROUP_MATERI_ELASTI_YOUNG, index, idum, &young, one,
      VERSION_NORMAL, PUT );
    db( GROUP_MATERI_ELASTI_POISSON, index, idum, &poisson, one,
      VERSION_NORMAL, PUT );
  }
  db_delete( GROUP_MATERI_PLASTI_HYPO_WOLFERSDORFF, VERSION_NORMAL );
  db_max_index( GROUP_SPRING_STIFFNESS_NONLINEAR, max, VERSION_NORMAL, GET );
  for ( index=0; index<=max; index++ ) {
    if ( !db_active_index( GROUP_SPRING_STIFFNESS_NONLINEAR, index,
      VERSION_NORMAL ) ) continue;
      // the stiffness at strain 0 (epsilon0 k0 epsilon1 k1 ...); outside
      // the table the nearest end value
    dval = db_dbl( GROUP_SPRING_STIFFNESS_NONLINEAR, index, VERSION_NORMAL );
    len = db_len( GROUP_SPRING_STIFFNESS_NONLINEAR, index, VERSION_NORMAL );
    table = get_new_dbl( len );
    array_move( dval, table, len );
    found = table_xy( table, "group_spring_stiffness_nonlinear", len, 0., k );
    if ( !found && len>=2 ) k = ( table[0]>0. ) ? table[1] : table[len-1];
    delete[] table;
    db( GROUP_SPRING_STIFFNESS, index, idum, &k, one, VERSION_NORMAL, PUT );
  }
  db_delete( GROUP_SPRING_STIFFNESS_NONLINEAR, VERSION_NORMAL );

    // plasticity
  db_delete( CONTROL_PLASTI_APPLY, VERSION_NORMAL );
  db( PLASTI_APPLY, 0, &no, ddum, one, VERSION_NORMAL, PUT );
  general_delete_prefix( "group_materi_plasti_", -1 );
  general_delete_prefix( "group_interface_materi_plasti_", -1 );

    // mesh -fixed_in_space
  array_set( options_mesh, -FIXED_IN_SPACE, ndim );
  db( OPTIONS_MESH, 0, options_mesh, ddum, ndim, VERSION_NORMAL, PUT );

    // every group_*_memory record: -total_linear with materi_displacement,
    // else -updated_linear
  memory = materi_displacement ? -TOTAL_LINEAR : -UPDATED_LINEAR;
  for ( idat=0; idat<MDAT; idat++ ) {
    item_name = db_name( idat );
    n = strlen( item_name );
    if ( strncmp( item_name, "group_", 6 ) || n<7 ||
         strcmp( item_name+n-7, "_memory" ) || db_type(idat)!=INTEGER )
      continue;
    db_max_index( idat, max, VERSION_NORMAL, GET );
    for ( index=0; index<=max; index++ ) {
      if ( db_active_index( idat, index, VERSION_NORMAL ) )
        db_int( idat, index, VERSION_NORMAL )[0] = memory;
    }
  }

    // dependencies of group_* data on a dof
  db_max_index( DEPENDENCY_ITEM, max, VERSION_NORMAL, GET );
  if ( max>=0 && nuknwn>0 ) {
    db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET );
    for ( index=0; index<=max; index++ ) {
      if ( !db_active_index( DEPENDENCY_ITEM, index, VERSION_NORMAL ) )
        continue;
      val = db_int( DEPENDENCY_ITEM, index, VERSION_NORMAL );
      iuknwn = -1;
      if ( !strncmp( db_name(labs(val[0])), "group_", 6 ) && val[2]<0 &&
           array_member( dof_label, val[2], nuknwn, iuknwn ) ) {
        db_delete_index( DEPENDENCY_ITEM, index, VERSION_NORMAL );
        db_delete_index( DEPENDENCY_DIAGRAM, index, VERSION_NORMAL );
      }
    }
  }

    // deleted nonlinear features
  db_delete( GROUP_INTERFACE_GAP, VERSION_NORMAL );
  db_delete( GROUP_TRUSS_ROPE, VERSION_NORMAL );
  general_delete_prefix( "group_materi_damage", -1 );
  general_delete_prefix( "group_materi_failure", -1 );
  general_delete_prefix( "groundflow_seepage_", -1 );
  general_delete_prefix( "contact_", -1 );
  db( CONTACT_APPLY, 0, &no, ddum, one, VERSION_NORMAL, PUT );
}
