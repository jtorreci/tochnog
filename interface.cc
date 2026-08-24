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

// interface_element - interface elements (Carril A, Fases 1 y 3).
//
// An interface element models a joint/discontinuity between two blocks
// of material. Its strains are the DISPLACEMENT DIFFERENCES between the
// two opposite sides of the element (not field gradients). In 2D the
// element is a quadrilateral with 4 nodes: nodes {0,1} form side 1 and
// nodes {2,3} form side 2.
//
// Constitutive law (group_interface_*):
//   - elastic stiffness (Fase 1):
//       group_interface_materi_elasti_stiffness kn kt,first kt,second
//       stress_normal = kn * strain_normal
//       stress_shear  = kt * 2 * strain_shear
//   - gap (Fase 3, RF-3): the interface is OPEN when the accumulated
//       normal strain <= gap (only residual stiffness acts), CLOSED when
//       strain_normal > gap. Compression (strain_normal > 0) always closes;
//       a physical gap is a NEGATIVE gap value. Default gap = 1.e20
//       (always closed).
//   - tension limit (Fase 3, RF-2): the interface opens in traction when
//       the TOTAL accumulated normal force |kn*strain_normal| exceeds the
//       limit: group_interface_materi_plasti_tension_direct tension_limit
//   - Mohr-Coulomb (Fase 3, RF-1): cumulative. The friction limit applies
//       to the TOTAL tangential force F_t (history ELEMENT_INTERFACE_FORCE_TANG):
//       group_interface_materi_plasti_mohr_coul_direct phi c phi_flow
//       trial = F_t,old + kt*2*du_tang, clamped to +/- max_fric with
//       max_fric = max(c + Fn*tan(phi), 0), Fn = kn*strain_normal (total).
//       Active by the PRESENCE of the record (phi=0,c=0 -> max_fric=0 ->
//       free sliding). The assembled rhs increment is F_t - F_t,old and the
//       tangential stiffness is 0 while plastic.
//   - dilatancy (Fase 3, RF-4): when the tangential force plastifies,
//       strain_normal += -|du_tang|*tan(phi_flow) (plastic normal opening).
//   - residual stiffness (Fase 3):
//       group_interface_materi_residual_stiffness factor
//       (fraction of the original stiffness used in opened interfaces)
//
// Strategy: implicit penalty + control_timestep_iterations (the tochnog
// Newton scheme corrects interpenetration within the step). The stiffness
// is updated within the iterations as the interface opens/closes.
void interface_element( long int element, long int name,
  long int element_group, double coord[], double old_dof[], double new_dof[], 
  double element_lhside[], double element_matrix[], double element_rhside[] )

{
  long int idim=0, jdim=0, inol=0, jnol=0, indx=0, swit=0, ldum=0, 
    nnol=4, mc_active=0, plastified=0, memory=-UPDATED_LINEAR, idum[1];
  double dtime=0., kn=0., kt1=0., kt2=0., tmp=0., ddum[1],
    normal[MDIM], tangent[MDIM], tangent2[MDIM], du[MDIM],
    du_norm=0., du_tang=0., du_tang2=0., stress_normal=0., stress_shear=0.,
    stress_shear2=0., strain_normal=0., strain_eff=0., force_norm=0., gap=0., tension_limit=0.,
    residual_factor=0.01, phi=0., c=0., phi_flow=0., max_fric=0.,
    stiff_normal=0., stiff_tang=0., stiff_tang2=0., ddum3[3],
    f_t_old=0., f_t=0., f_t2_old=0., f_t2=0., trial=0., trial2=0.,
    ft_mag=0., fn_total=0., force_gravity[MDIM];
  long int *nodes=NULL;

  swit = set_swit(element,-1,"interface_element");
  if ( swit ) pri( "In routine INTERFACE_ELEMENT." );

  if ( !db_active_index( GROUP_INTERFACE, element_group, VERSION_NORMAL ) )
    return;

  db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET );

  // group parameters
  // ddum3 MUST be zeroed: with GET_IF_EXISTS on a missing record db() does
  // not write dval, and an uninitialized buffer made kn/kt garbage (NaN in
  // the assembled matrix with gcc -O1; exposed 2026-08-24 by the clean
  // rebuild on a newer compiler).
  array_set( ddum3, 0., 3 );
  db( GROUP_INTERFACE_MATERI_ELASTI_STIFFNESS, element_group, idum, ddum3,
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  kn = ddum3[0]; kt1 = ddum3[1]; kt2 = ddum3[2];
  db( GROUP_INTERFACE_MATERI_PLASTI_TENSION_DIRECT, element_group, idum,
    &tension_limit, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( GROUP_INTERFACE_MATERI_RESIDUAL_STIFFNESS, element_group, idum,
    &residual_factor, ldum, VERSION_NORMAL, GET_IF_EXISTS );

  // memory model (Fase 3, group_interface_materi_memory): -updated_linear
  // (default) recomputes the interface normal/tangent from the current
  // (deformed) configuration each step; -total_linear uses the time-0
  // reference geometry (NODE_START_REFINED), so the interface keeps its
  // original orientation (elem.cc pattern: TOTAL_LINEAR -> fixed reference).
  db( GROUP_INTERFACE_MATERI_MEMORY, element_group, &memory, ddum,
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  if ( memory!=-UPDATED_LINEAR && memory!=-TOTAL_LINEAR )
    db_error( GROUP_INTERFACE_MATERI_MEMORY, element_group );

  // node numbers (needed to read the time-0 reference geometry)
  long int length_el=0, *el=NULL;
  el = get_new_int(MAXIMUM_NODE+1);
  nodes = get_new_int(MAXIMUM_NODE);
  db( ELEMENT, element, el, ddum, length_el, VERSION_NORMAL, GET );
  array_move( &el[1], nodes, length_el-1 );

  // element type: 2D interface is a quadrilateral (converted bar2);
  // 3D interfaces are a prism6 (converted tria3) or hex8 (converted quad4).
  // A -bar2 reaching here means it was not converted by control_mesh_convert;
  // treat it as a degenerate 2D interface (side 1 = node 0, side 2 = node 1).
  if ( name==-BAR2 ) {
    nnol = 2;
  }
  else if ( name==-QUAD4 ) {
    nnol = 4;
  }
  else if ( name==-PRISM6 ) {
    nnol = 6;
  }
  else {
    assert( name==-HEX8 );
    nnol = 8;
  }

  // interface frame: normal perpendicular to the interface surface and two
  // in-plane tangents. 2D: tangent along side 1 (nodes 0,1), normal
  // perpendicular (in-plane). 3D: normal = cross product of two side-1
  // edges (surface normal), tangent2 completes the orthonormal frame.
  if ( ndim==2 ) {
    double *ca, *cb;
    if ( memory==-TOTAL_LINEAR ) {
      // time-0 reference geometry (NODE_START_REFINED): the interface keeps
      // its original orientation even when the mesh deforms.
      ca = db_dbl( NODE_START_REFINED, nodes[0], VERSION_NORMAL );
      cb = db_dbl( NODE_START_REFINED, nodes[1], VERSION_NORMAL );
    }
    else {
      ca = &coord[0*ndim];
      cb = &coord[1*ndim];
    }
    tangent[0] = cb[0] - ca[0];
    tangent[1] = cb[1] - ca[1];
    array_normalize( tangent, ndim );
    normal[0] = -tangent[1];
    normal[1] =  tangent[0];
    array_set( tangent2, 0., MDIM );
  }
  else {
    // 3D: the two side-1 edges define the surface plane.
    double *c0, *c1, *c2;
    if ( memory==-TOTAL_LINEAR ) {
      c0 = db_dbl( NODE_START_REFINED, nodes[0], VERSION_NORMAL );
      c1 = db_dbl( NODE_START_REFINED, nodes[1], VERSION_NORMAL );
      c2 = db_dbl( NODE_START_REFINED, nodes[2], VERSION_NORMAL );
    }
    else {
      c0 = &coord[0*ndim];
      c1 = &coord[1*ndim];
      c2 = &coord[2*ndim];
    }
    double e1[MDIM], e2[MDIM];
    for ( idim=0; idim<3; idim++ ) {
      e1[idim] = c1[idim] - c0[idim];
      e2[idim] = c2[idim] - c0[idim];
    }
    // normal = e1 x e2 (surface normal), tangent = e1, tangent2 = normal x e1
    normal[0] = e1[1]*e2[2] - e1[2]*e2[1];
    normal[1] = e1[2]*e2[0] - e1[0]*e2[2];
    normal[2] = e1[0]*e2[1] - e1[1]*e2[0];
    array_normalize( normal, 3 );
    // orientation: the 2D convention is normal pointing FROM side 2 TOWARDS
    // side 1 (so compression, side 2 moving towards side 1, gives a POSITIVE
    // normal strain). Flip the cross product if it points the wrong way.
    {
      double *cm1, *cm2;
      long int ns1 = nnol/2;
      if ( memory==-TOTAL_LINEAR ) {
        cm1 = db_dbl( NODE_START_REFINED, nodes[0], VERSION_NORMAL );
        cm2 = db_dbl( NODE_START_REFINED, nodes[ns1], VERSION_NORMAL );
      }
      else {
        cm1 = &coord[0*ndim];
        cm2 = &coord[ns1*ndim];
      }
      double dir = 0.;
      for ( idim=0; idim<3; idim++ ) dir += normal[idim]*(cm1[idim]-cm2[idim]);
      if ( dir<0. ) {
        normal[0] = -normal[0]; normal[1] = -normal[1]; normal[2] = -normal[2];
      }
    }
    tangent[0] = e1[0]; tangent[1] = e1[1]; tangent[2] = e1[2];
    array_normalize( tangent, 3 );
    tangent2[0] = normal[1]*tangent[2] - normal[2]*tangent[1];
    tangent2[1] = normal[2]*tangent[0] - normal[0]*tangent[2];
    tangent2[2] = normal[0]*tangent[1] - normal[1]*tangent[0];

    // group_interface_tangential_reference_point (manual Professional
    // 6.636): defines the first tangential direction from a reference
    // point. t1 = the part of (ref_point - element_centroid)
    // perpendicular to the normal; t2 = normal x t1. 3D only; falls back
    // to the geometric tangent if the reference point is on the normal
    // line through the centroid.
    if ( db( GROUP_INTERFACE_TANGENTIAL_REFERENCE_POINT, element_group,
        idum, ddum3, ldum, VERSION_NORMAL, GET_IF_EXISTS ) ) {
      double centroid[MDIM], vref[MDIM], dotn = 0.;
      array_set( centroid, 0., MDIM );
      if ( memory==-TOTAL_LINEAR ) {
        for ( inol=0; inol<nnol; inol++ ) {
          double *cn = db_dbl( NODE_START_REFINED, nodes[inol], VERSION_NORMAL );
          for ( idim=0; idim<3; idim++ ) centroid[idim] += cn[idim]/nnol;
        }
      }
      else {
        for ( inol=0; inol<nnol; inol++ )
          for ( idim=0; idim<3; idim++ )
            centroid[idim] += coord[inol*ndim+idim]/nnol;
      }
      for ( idim=0; idim<3; idim++ ) vref[idim] = ddum3[idim] - centroid[idim];
      for ( idim=0; idim<3; idim++ ) dotn += vref[idim]*normal[idim];
      for ( idim=0; idim<3; idim++ ) vref[idim] -= dotn*normal[idim];
      if ( array_size( vref, 3 ) > 1.e-12 ) {
        for ( idim=0; idim<3; idim++ ) tangent[idim] = vref[idim];
        array_normalize( tangent, 3 );
        tangent2[0] = normal[1]*tangent[2] - normal[2]*tangent[1];
        tangent2[1] = normal[2]*tangent[0] - normal[0]*tangent[2];
        tangent2[2] = normal[0]*tangent[1] - normal[1]*tangent[0];
      }
    }
  }

  // velocity difference between the sides (side2 - side1) on the
  // velocity dof -> incremental displacement difference (spring2 pattern)
  // For a -bar2 (not converted): node 0 = side 1, node 1 = side 2.
  // 2D quad4: sides {0,1} and {2,3}. 3D prism6: {0,1,2}/{3,4,5},
  // hex8: {0,1,2,3}/{4,5,6,7} (side 1 = first half, side 2 = second half).
  if ( name==-BAR2 ) {
    for ( idim=0; idim<ndim; idim++ ) {
      double v1 = new_dof[0*nuknwn+vel_indx+idim*nder];
      double v2 = new_dof[1*nuknwn+vel_indx+idim*nder];
      du[idim] = ( v2 - v1 ) * dtime;
    }
  }
  else {
    long int ns1 = nnol/2;   // nodes on side 1
    for ( idim=0; idim<ndim; idim++ ) {
      double v_side1 = 0., v_side2 = 0.;
      for ( inol=0; inol<ns1; inol++ ) {
        v_side1 += new_dof[inol*nuknwn+vel_indx+idim*nder];
        v_side2 += new_dof[(inol+ns1)*nuknwn+vel_indx+idim*nder];
      }
      v_side1 /= ns1; v_side2 /= ns1;
      du[idim] = ( v_side2 - v_side1 ) * dtime;
    }
  }

  du_norm  = array_inproduct( du, normal, ndim );
  du_tang  = array_inproduct( du, tangent, ndim );
  du_tang2 = ( ndim==3 ) ? array_inproduct( du, tangent2, ndim ) : 0.;

  // accumulated normal strain (history, VERSION_NORMAL). Sign convention:
  // compression is POSITIVE (verified empirically). The history is read
  // BEFORE adding the current step's increment so that gap / tension /
  // Mohr-Coulomb decisions see the accumulated total.
  db( ELEMENT_INTERFACE_STRAIN_NORMAL, element, idum, &strain_normal,
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  strain_normal += du_norm;

  // group_interface_materi_expansion_normal (manual Professional 6.630):
  // thermal strain expansion in interface thickness direction per unit
  // temperature; the temperature is the average of both sides. The
  // MECHANICAL strain seen by gap / tension / Mohr-Coulomb is the
  // accumulated strain minus the thermal expansion; the stored history
  // stays purely mechanical (no thermal re-counting per step). Only
  // meaningful with condif_temperature.
  strain_eff = strain_normal;
  if ( condif_temperature ) {
    double alpha_n = 0., t_side1 = 0., t_side2 = 0.;
    if ( db( GROUP_INTERFACE_MATERI_EXPANSION_NORMAL, element_group, idum,
        &alpha_n, ldum, VERSION_NORMAL, GET_IF_EXISTS ) ) {
      long int ns1 = nnol/2;
      double t1_old = 0., t2_old = 0.;
      for ( inol=0; inol<ns1; inol++ ) {
        t_side1 += new_dof[inol*nuknwn+temp_indx];
        t_side2 += new_dof[(inol+ns1)*nuknwn+temp_indx];
        t1_old   += old_dof[inol*nuknwn+temp_indx];
        t2_old   += old_dof[(inol+ns1)*nuknwn+temp_indx];
      }
      t_side1 /= ns1; t_side2 /= ns1; t1_old /= ns1; t2_old /= ns1;
      // total thermal expansion (for gap / tension / Mohr-Coulomb state)
      strain_eff = strain_normal - alpha_n * 0.5*(t_side1+t_side2);
      // incremental thermal expansion acts as a pseudo-load this step:
      // the normal force increment is stiff*(du_norm - d_alpha*T), the
      // same incremental pattern as stress.cc thermal strains
      du_norm -= alpha_n * ( 0.5*(t_side1+t_side2) - 0.5*(t1_old+t2_old) );
    }
  }

  // accumulated total tangential force (history, VERSION_NORMAL). Default 0:
  // without the Mohr-Coulomb record the interface stays purely elastic
  // (Fase 1 behavior). 3D: two tangential components (tangent and tangent2).
  f_t_old = 0.;
  db( ELEMENT_INTERFACE_FORCE_TANG, element, idum, &f_t_old, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  f_t2_old = 0.;
  if ( ndim==3 )
    db( ELEMENT_INTERFACE_FORCE_TANG2, element, idum, &f_t2_old, ldum,
      VERSION_NORMAL, GET_IF_EXISTS );

  // gap (FIX 3, RF-3): the interface is OPEN when strain_normal <= gap
  // (only residual stiffness acts), CLOSED when strain_normal > gap.
  // Compression (strain_normal > 0) always closes the interface. A physical
  // gap is a NEGATIVE gap value: the interface stays open until compression
  // exceeds |gap|. If no gap is specified the interface is always closed:
  // default gap = -1e20 (the OLD code used +1e20, valid only for the old
  // condition strain >= gap; with the inverted condition strain <= gap a
  // positive default would leave the interface ALWAYS open).
  if ( !db( GROUP_INTERFACE_GAP, element_group, idum, &gap, ldum,
      VERSION_NORMAL, GET_IF_EXISTS ) )
    gap = -1.e20;
  stiff_normal = kn;
  if ( strain_eff <= gap ) {
    stiff_normal = kn * residual_factor;
  }
  force_norm = stiff_normal * du_norm;
  // tension limit (FIX 2, RF-2): the interface opens in traction when the
  // TOTAL accumulated normal force |Fn_total| = |kn*strain_normal| exceeds
  // the limit, and only if it was still closed (stiff_normal==kn). On
  // opening: residual stiffness and normal force capped at the limit.
  fn_total = kn * strain_eff;
  if ( tension_limit>0. && strain_normal<0. && fabs(fn_total)>tension_limit
      && stiff_normal==kn ) {
    stiff_normal = kn * residual_factor;
    force_norm = ( du_norm>=0. ) ? tension_limit : -tension_limit;
  }
  stress_normal = force_norm;

  // cumulative Mohr-Coulomb (FIX 1, RF-1): the friction limit applies to
  // the TOTAL tangential force, not the per-step force. The MC law is
  // active by the PRESENCE of the record (D2): phi=0,c=0 gives max_fric=0
  // -> free sliding; without the record the interface stays purely elastic
  // (Fase 1 behavior). 3D: the two tangential components are clamped
  // together on the magnitude (max_fric), preserving their ratio.
  mc_active = db( GROUP_INTERFACE_MATERI_PLASTI_MOHR_COUL_DIRECT, element_group,
    idum, ddum3, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  phi = ddum3[0]; c = ddum3[1]; phi_flow = ddum3[2];
  trial  = f_t_old  + kt1 * 2. * du_tang;
  trial2 = f_t2_old + kt2 * 2. * du_tang2;
  plastified = 0;
  if ( mc_active ) {
    max_fric = c + kn * strain_eff * tan( phi );
    if ( max_fric < 0. ) max_fric = 0.;   // D5: floor at 0
    if ( ndim==3 ) {
      // clamp the magnitude of (trial, trial2) at max_fric
      ft_mag = sqrt( trial*trial + trial2*trial2 );
      if ( ft_mag > max_fric && ft_mag>0. ) {
        double scale = max_fric / ft_mag;
        trial *= scale; trial2 *= scale; plastified = 1;
      }
    }
    else {
      if      ( trial >  max_fric ) { trial =  max_fric; plastified = 1; }
      else if ( trial < -max_fric ) { trial = -max_fric; plastified = 1; }
    }
  }
  f_t  = trial;
  f_t2 = trial2;
  stiff_tang  = ( mc_active && plastified ) ? 0. : kt1 * 2.;
  stiff_tang2 = ( mc_active && plastified ) ? 0. : kt2 * 2.;

  // dilatancy (FIX 4, RF-4): plastic slip opens the interface by
  // du_n^p = -|du_tang|*tan(phi_flow) (always opening; compression positive).
  // 3D: magnitude of the tangential slip in both in-plane directions.
  if ( plastified && phi_flow>0. ) {
    double du_tang_mag = sqrt( du_tang*du_tang + du_tang2*du_tang2 );
    strain_normal += -fabs( du_tang_mag ) * tan( phi_flow );
  }

  // the rhs carries the INCREMENT F_t - F_t,old (D3): == kt*2*du_tang when
  // elastic (Fase 1 backward compatible), == clamped increment when plastic.
  stress_shear  = f_t  - f_t_old;
  stress_shear2 = f_t2 - f_t2_old;

  if ( swit ) {
    pri( "du_norm", du_norm );
    pri( "du_tang", du_tang );
    pri( "du_tang2", du_tang2 );
    pri( "strain_normal", strain_normal );
    pri( "stress_normal", stress_normal );
    pri( "f_t_old", f_t_old );
    pri( "f_t", f_t );
    pri( "stress_shear", stress_shear );
  }

  // assembly: nodal force -sign*(stress*dir) and stiffness matrix on the
  // velocity dofs (pattern spring.cc). side 1 = first half of the nodes,
  // side 2 = second half (2D quad4: {0,1}/{2,3}; prism6 {0,1,2}/{3,4,5};
  // hex8 {0,1,2,3}/{4,5,6,7}).
  long int ns1 = nnol/2;
  for ( idim=0; idim<ndim; idim++ ) {
    double dirn = normal[idim], dirt = tangent[idim], dirt2 = tangent2[idim];
    for ( inol=0; inol<nnol; inol++ ) {
      double sign = ( inol>=ns1 ) ? +1. : -1.;
      indx = inol*npuknwn + (vel_indx+idim*nder)/nder;
      tmp = -sign*( stress_normal*dirn + stress_shear*dirt +
        stress_shear2*dirt2 );
      element_rhside[indx] += tmp;
      for ( jnol=0; jnol<nnol; jnol++ ) {
        double jsign = ( jnol>=ns1 ) ? +1. : -1.;
        for ( jdim=0; jdim<ndim; jdim++ ) {
          double jdirn = normal[jdim], jdirt = tangent[jdim],
            jdirt2 = tangent2[jdim];
          double kkk = sign*jsign*( stiff_normal*dirn*jdirn +
            stiff_tang*dirt*jdirt + stiff_tang2*dirt2*jdirt2 );
          long int jndx = inol*npuknwn*nnol*npuknwn +
            ((vel_indx+idim*nder)/nder)*nnol*npuknwn +
            jnol*npuknwn + (vel_indx+jdim*nder)/nder;
          element_matrix[jndx] += kkk * dtime;
          if ( jnol==inol && jdim==idim )
            element_lhside[inol*npuknwn+(vel_indx+idim*nder)/nder] +=
              kkk * dtime;
        }
      }
    }
  }

  // store accumulated histories (used by gap / tension / Mohr-Coulomb and
  // the next step's cumulative trial)
  ldum = 1;
  db( ELEMENT_INTERFACE_STRAIN_NORMAL, element, idum, &strain_normal,
    ldum, VERSION_NEW, PUT );
  db( ELEMENT_INTERFACE_FORCE_TANG, element, idum, &f_t, ldum,
    VERSION_NEW, PUT );
  if ( ndim==3 )
    db( ELEMENT_INTERFACE_FORCE_TANG2, element, idum, &f_t2, ldum,
      VERSION_NEW, PUT );

  // group_interface_groundflow_permeability: ground flow THROUGH the
  // interface. The interface connects the pore pressures on both sides;
  // the flux across it is q = pe * (pres_side1 - pres_side2) per unit
  // length (2D) or area (3D), assembled on the pressure dofs of the
  // facing node pairs. group_interface_groundflow_capacity adds a storage
  // term on the pressure dofs of the interface nodes.
  // group_interface_groundflow_total_pressure_tension: when the accumulated
  // normal strain exceeds strain_normal_minimum (crack open), the static
  // water pressure from water_height is used instead of the pore pressure
  // from the groundflow equation when it is larger in absolute value.
  if ( groundflow_pressure ) {
    double pe_iface = 0., C_iface = 0.;
    long int has_pe = 0, has_C = 0;
    has_pe = db( GROUP_INTERFACE_GROUNDFLOW_PERMEABILITY, element_group, idum,
      &pe_iface, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    has_C  = db( GROUP_INTERFACE_GROUNDFLOW_CAPACITY, element_group, idum,
      &C_iface, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    if ( has_pe || has_C || db_active_index(
        GROUP_INTERFACE_GROUNDFLOW_TOTAL_PRESSURE_TENSION, element_group,
        VERSION_NORMAL ) ) {
      // storage term on the pressure dofs of the interface nodes (lumped)
      if ( has_C ) {
        for ( inol=0; inol<nnol; inol++ ) {
          long int jndx = inol*npuknwn + pres_indx/nder;
          double dpres = ( new_dof[inol*nuknwn+pres_indx] -
            old_dof[inol*nuknwn+pres_indx] ) / dtime;
          element_rhside[jndx] -= C_iface * dpres;
        }
      }
      // through-interface flux coupling the facing node pairs
      if ( has_pe ) {
        long int ns1 = nnol/2;
        for ( inol=0; inol<ns1; inol++ ) {
          long int jn1 = inol*npuknwn + pres_indx/nder;
          long int jn2 = (inol+ns1)*npuknwn + pres_indx/nder;
          double pres1 = new_dof[inol*nuknwn+pres_indx];
          double pres2 = new_dof[(inol+ns1)*nuknwn+pres_indx];
          double q = pe_iface * ( pres1 - pres2 );
          element_rhside[jn1] -= q;
          element_rhside[jn2] += q;
          element_matrix[jn1*nnol*npuknwn+jn1] += pe_iface;
          element_matrix[jn1*nnol*npuknwn+jn2] -= pe_iface;
          element_matrix[jn2*nnol*npuknwn+jn1] -= pe_iface;
          element_matrix[jn2*nnol*npuknwn+jn2] += pe_iface;
        }
      }
      // crack pressure: static water pressure on an opened interface
      if ( db_active_index( GROUP_INTERFACE_GROUNDFLOW_TOTAL_PRESSURE_TENSION,
          element_group, VERSION_NORMAL ) ) {
        double gitpt[2], dens=0.;
        db( GROUP_INTERFACE_GROUNDFLOW_TOTAL_PRESSURE_TENSION, element_group,
          idum, gitpt, ldum, VERSION_NORMAL, GET );
        db( GROUNDFLOW_DENSITY, 0, idum, &dens, ldum, VERSION_NORMAL,
          GET_IF_EXISTS );
        force_gravity_calculate( force_gravity );
        if ( strain_eff > gitpt[0] && dens>0. ) {
          for ( inol=0; inol<nnol; inol++ ) {
            long int jndx = inol*npuknwn + pres_indx/nder;
            double static_pres =
              force_gravity[ndim-1] * dens * gitpt[1];
            double pres_n = new_dof[inol*nuknwn+pres_indx];
            if ( scalar_dabs(static_pres) > scalar_dabs(pres_n) )
              element_rhside[jndx] += static_pres - pres_n;
          }
        }
      }
    }
  }

  // group_interface_condif_conductivity (manual Professional 6.625): heat
  // flow through the interface per unit temperature difference between the
  // facing sides, q = k * (T_side1 - T_side2), assembled on the temperature
  // dofs of the facing node pairs (thermal analog of
  // group_interface_groundflow_permeability). k is the conductivity of the
  // layer simulated by the interface (thermal thickness included), not the
  // material conductivity.
  if ( condif_temperature ) {
    double k_iface = 0.;
    if ( db( GROUP_INTERFACE_CONDIF_CONDUCTIVITY, element_group, idum,
        &k_iface, ldum, VERSION_NORMAL, GET_IF_EXISTS ) ) {
      long int ns1 = nnol/2;
      for ( inol=0; inol<ns1; inol++ ) {
        long int jn1 = inol*npuknwn + temp_indx/nder;
        long int jn2 = (inol+ns1)*npuknwn + temp_indx/nder;
        double temp1 = new_dof[inol*nuknwn+temp_indx];
        double temp2 = new_dof[(inol+ns1)*nuknwn+temp_indx];
        double q = k_iface * ( temp1 - temp2 );
        element_rhside[jn1] -= q;
        element_rhside[jn2] += q;
        element_matrix[jn1*nnol*npuknwn+jn1] += k_iface;
        element_matrix[jn1*nnol*npuknwn+jn2] -= k_iface;
        element_matrix[jn2*nnol*npuknwn+jn1] -= k_iface;
        element_matrix[jn2*nnol*npuknwn+jn2] += k_iface;
      }
    }
  }

  delete[] el;
  delete[] nodes;

  if ( swit ) pri( "Out function INTERFACE_ELEMENT" );
}

// interface_convert - control_mesh_convert (Carril A, Fase 2).
//
// Converts low-dimensional interface elements to their isoparametric
// equivalent, creating the nodes of the opposite interface side:
//   -bar2 -> -quad4  (2D): the bar2 nodes {a,b} form side 1; two new
//   nodes {a',b'} are created by copying {a,b} shifted along the
//   interface normal (side 2). The element is rewritten as -quad4
//   {a,b,a',b'}. Neighbouring isoparametric elements on the OTHER side
//   of the interface (i.e. not in the groups listed by
//   control_mesh_convert_element_group) are reconnected to {a',b'}.
//
// See ProjectDocs/DESIGN-INTERFACES.md for the full algorithm.
void interface_convert( long int icontrol )

{
  long int element=0, max_element=0, i=0, j=0,
    name=0, length=0, ldum=0, swit=0, element_group=0,
    max_node=0, max_node_old=0, length_convert_groups=0, found=0, nconv=0,
    idum[1], *el=NULL, *convert_groups=NULL, *node_element=NULL;
  double ddum[1], coord[MDIM], normal[MDIM], tangent[MDIM], shift=0.,
    *ca=NULL, *cb=NULL;

  swit = set_swit(-1,-1,"interface_convert");
  if ( swit ) pri( "In routine INTERFACE_CONVERT." );

  if ( !db_active_index( CONTROL_MESH_CONVERT, icontrol, VERSION_NORMAL ) )
    return;

  long int convert_switch = -YES;
  db( CONTROL_MESH_CONVERT, icontrol, &convert_switch, ddum, ldum,
    VERSION_NORMAL, GET );
  if ( convert_switch!= -YES ) return;

  el = get_new_int(MAXIMUM_NODE+1);
  db_highest_index( ELEMENT, max_element, VERSION_NORMAL );
  db_highest_index( NODE, max_node, VERSION_NORMAL );
  max_node_old = max_node;

  // element groups located on ONE side of the interfaces; neighbours in
  // these groups keep the original nodes, the others get the new ones
  convert_groups = get_new_int(DATA_ITEM_SIZE);
  db( CONTROL_MESH_CONVERT_ELEMENT_GROUP, icontrol, convert_groups,
    ddum, length_convert_groups, VERSION_NORMAL, GET_IF_EXISTS );

  nconv = 0;
  for ( element=0; element<=max_element; element++ ) {
    if ( !db_active_index( ELEMENT, element, VERSION_NORMAL ) ) continue;
    db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
    name = el[0];
    element_group = 0;
    db( ELEMENT_GROUP, element, &element_group, ddum, ldum,
      VERSION_NORMAL, GET_IF_EXISTS );
    if ( name!=-BAR2 && name!=-TRIA3 && name!=-QUAD4 ) continue;
    if ( !db_active_index( GROUP_INTERFACE, element_group, VERSION_NORMAL ) )
      continue;

    // side 1 nodes: for 2D bar2 = {a,b}; for 3D tria3 = 3 nodes,
    // quad4 = 4 nodes. They form the base of the interface element.
    long int ns1 = 0;
    if      ( name==-BAR2  ) ns1 = 2;
    else if ( name==-TRIA3 ) ns1 = 3;
    else                     ns1 = 4;

    // tangent along the first edge, normal perpendicular to the surface
    array_set( tangent, 0., MDIM );
    if ( ndim==2 ) {
      ca = db_dbl( NODE, el[1], VERSION_NORMAL );
      cb = db_dbl( NODE, el[2], VERSION_NORMAL );
      tangent[0] = cb[0] - ca[0];
      tangent[1] = cb[1] - ca[1];
      array_normalize( tangent, ndim );
      normal[0] = -tangent[1];
      normal[1] =  tangent[0];
    }
    else {
      // 3D: normal = cross product of two side-1 edges
      ca = db_dbl( NODE, el[1], VERSION_NORMAL );
      cb = db_dbl( NODE, el[2], VERSION_NORMAL );
      double *cc = db_dbl( NODE, el[3], VERSION_NORMAL );
      double e1[MDIM], e2[MDIM];
      for ( i=0; i<3; i++ ) {
        e1[i] = cb[i] - ca[i];
        e2[i] = cc[i] - ca[i];
      }
      normal[0] = e1[1]*e2[2] - e1[2]*e2[1];
      normal[1] = e1[2]*e2[0] - e1[0]*e2[2];
      normal[2] = e1[0]*e2[1] - e1[1]*e2[0];
      array_normalize( normal, 3 );
      tangent[0] = e1[0]; tangent[1] = e1[1]; tangent[2] = e1[2];
      array_normalize( tangent, 3 );
    }

    // shift: small fraction of the first side length (interface thickness)
    double len = 0.;
    for ( i=0; i<ndim; i++ ) len += tangent[i]*tangent[i];
    shift = 0.01 * sqrt( len );

    // create ns1 new nodes: copies of the side-1 nodes shifted along n
    for ( j=0; j<ns1; j++ ) {
      long int src = el[1+j];
      long int dst = ++max_node;
      db( NODE, src, idum, coord, ldum, VERSION_NORMAL, GET );
      for ( i=0; i<ndim; i++ ) coord[i] += shift*normal[i];
      db( NODE, dst, idum, coord, ldum, VERSION_NORMAL, PUT );
      db( NODE_START_REFINED, src, idum, coord, ldum, VERSION_NORMAL, GET );
      for ( i=0; i<ndim; i++ ) coord[i] += shift*normal[i];
      db( NODE_START_REFINED, dst, idum, coord, ldum, VERSION_NORMAL, PUT );
      double *ndof = db_dbl( NODE_DOF, src, VERSION_NORMAL );
      long int ln = db_len( NODE_DOF, src, VERSION_NORMAL );
      db( NODE_DOF, dst, idum, ndof, ln, VERSION_NORMAL, PUT );
      double *ndof_sr = db_dbl( NODE_DOF_START_REFINED, src, VERSION_NORMAL );
      long int ln_sr = db_len( NODE_DOF_START_REFINED, src, VERSION_NORMAL );
      db( NODE_DOF_START_REFINED, dst, idum, ndof_sr, ln_sr, VERSION_NORMAL, PUT );
      length = 1;
      db( NODE_MACRO_GENERATE, dst, &icontrol, ddum, length, VERSION_NORMAL, PUT );
      // side 2 node (same ordering as side 1)
      el[1+ns1+j] = dst;
    }

    // rewrite the element: bar2->quad4, tria3->prism6, quad4->hex8.
    // el[0]=name, el[1..ns1]=side1, el[ns1+1..2*ns1]=side2 (already filled).
    if ( name==-BAR2 )      el[0] = -QUAD4;
    else if ( name==-TRIA3 ) el[0] = -PRISM6;
    else                     el[0] = -HEX8;
    length = 1 + 2*ns1;
    db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, PUT );

    // reconnect neighbours on the OTHER side: elements sharing a side-1
    // node that are NOT in convert_groups get that node replaced by its
    // new duplicate.
    for ( j=0; j<ns1; j++ ) {
      long int src = el[1+j];
      long int dst = el[1+ns1+j];
      node_element = db_int( NODE_ELEMENT, src, VERSION_NORMAL );
      long int nel = db_len( NODE_ELEMENT, src, VERSION_NORMAL );
      for ( long int iel=0; iel<nel; iel++ ) {
        long int elnum = node_element[iel];
        if ( elnum==element ) continue;
        long int gr = 0;
        db( ELEMENT_GROUP, elnum, &gr, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
        // if the neighbour is in a convert group (one side), keep it
        found = 0;
        for ( long int ig=0; ig<length_convert_groups; ig++ )
          if ( convert_groups[ig]==gr ) { found = 1; break; }
        if ( found ) continue;
        // reconnect: replace src by dst in the neighbour connectivity
        long int ln = 0;
        long int *nel2 = get_new_int(MAXIMUM_NODE+1);
        db( ELEMENT, elnum, nel2, ddum, ln, VERSION_NORMAL, GET );
        for ( long int k=1; k<ln; k++ )
          if ( nel2[k]==src ) nel2[k] = dst;
        db( ELEMENT, elnum, nel2, ddum, ln, VERSION_NORMAL, PUT );
        delete[] nel2;
      }
    }
    nconv++;
  }

  delete[] el;
  delete[] convert_groups;

  if ( nconv>0 )
    mesh_has_changed( VERSION_NORMAL );

  if ( swit ) pri( "Out function INTERFACE_CONVERT" );
}
