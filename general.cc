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
    options_inertia=-YES, options_convection=-YES, options_axisymmetric=0,
    lagrange=0, stokes=-NO, principal_unknown=0, ldum=0, 
    idum[1], options_mesh[MDIM], dof_type[MUKNWN], dof_principal[MUKNWN];
  double condif_conductivity=0., dens=0.,
    condif_capacity=0., visc=0., vel=0.,
    C=0., dtime=0., inertia=0., radius=0., mat_diff=1.,
    conv_part=0., diff_part=0., peclet=0., peclet_factor=0., 
    val=0., val_max=0., val_min=0., val_new=0.,
    artificial_diffusion=0., tmp=0., D=0., diffusion=1.,
    element_lhside_add = 0., element_rhside_add=0.,
    ddum[1], pe[MDIM], node_remesh_velocity[MDIM], condif_flow[MDIM];
  long int element_dof_initial_active=0, edi_idum[1], edi_n=0;
  double edi_vals[DATA_ITEM_SIZE], element_dof_initial_values[MUKNWN];

  if ( type==-MAXWELL_FREQUENCY || type==-MAXWELL_TIME ) return;

  swit = set_swit(element,-1,"general");
  if ( swit ) pri( "In GENERAL" );

  db( DOF_TYPE, 0, dof_type, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( DOF_PRINCIPAL, 0, dof_principal, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET );
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );

  db( OPTIONS_MESH, 0, options_mesh, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( OPTIONS_CONVECTION, 0, &options_convection, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_OPTIONS_CONVECTION, icontrol, &options_convection, 
    ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( OPTIONS_INERTIA, 0, &options_inertia, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_OPTIONS_INERTIA, icontrol, &options_inertia, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
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
    get_group_data( GROUP_MATERI_VISCOSITY, element_group, element, new_unknowns, 
      &visc, ldum, GET_IF_EXISTS );
    db( GROUP_MATERI_STOKES, 0, &stokes, ddum, ldum, 
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
        if ( options_inertia==-YES || !principal_unknown ) {
          // stress dof recovery weight (DIAG lot C/D, fix D-b): for the
          // SRI quad4 (2x2 Gauss points) the h-weighted lumped average
          // dilutes the corner nodal stresses; the consistent recovery
          // is the Lagrange extrapolation of the Gauss-point values to
          // the nodes (the "same B at the node"). Applied to the NORMAL
          // stress components only (superconvergent at the Gauss
          // points); the SHEAR components keep h (the Q4 shear is not
          // superconvergent - the interpolation error dominates). For
          // every other quadrature the weight reduces to h
          // (node-containing rules).
          double weight = h[inol];
          if ( dof_type[iuknwn]==-MATERI_STRESS && sri_active(
              element, element_group, name, nnol ) ) {
            // normal components: stress_indx(0,0)=0, (1,1)=3, (2,2)=5
            // (times nder for the derivative slots)
            long int kcomp = iuknwn - stres_indx;
            if ( kcomp==0 || kcomp==3*nder || kcomp==5*nder )
              weight = sri_stress_recovery_weight( nnol, inol, npoint,
                ipoint, h[inol], 1 );
          }
          double past_dof = old_dof[inol*nuknwn+ipuknwn*nder];
          if ( element_dof_initial_active && ipuknwn<nuknwn )
            past_dof = element_dof_initial_values[ipuknwn];
          tmp = weight * inertia *
            ( new_dof[inol*nuknwn+ipuknwn*nder] - past_dof ) / dtime;
          element_rhside_add -= volume * tmp;
          element_lhside_add += volume * weight * inertia / dtime;
          if ( residue && dof_principal[iuknwn]>=0 ) element_residue[indx] -= tmp;
          element_matrix[indx*nnol*npuknwn+indx] += volume * weight * inertia / dtime;
        }
          
        for ( jdim=0; jdim<ndim; jdim++ ) {

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
          vel = 0.;
          if ( !lagrange ) {
            if ( materi_velocity ) 
              vel += new_dof[inol*nuknwn+vel_indx+jdim*nder];
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
  long int icontrol=0, swit=-YES, ldum=0, idum[1], max_index=-1;
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
