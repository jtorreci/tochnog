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
 
    On June 30, 2013 made change to correct proper calculation of velocity 
    as the case of velocity calculated from permeability x (dh/dx) 
    had wrong sign. F. Lorenzo.
*/

#include "tochnog.h"

// relative tolerance of the "total pressure reached the limit" test of
// groundflow_total_pressure_limit (interpolated capped values)
#define EPS_TOTAL_PRESSURE_LIMIT 1.e-10
#define GROUNDFLOW_TOTAL_PRESSURE_LIMIT_DEFAULT 0.

void groundflow( long int element, long int gr, long int nnol, long int nodes[],
  double coord_ip[], double h[], double d[], 
  double volume, double old_unknowns[], 
  double new_unknowns[], double grad_new_unknowns[], 
  double element_matrix[], double element_rhside[],
  double element_residue[] )

{
  long int swit=0, inol=0, jnol=0, jdim=0, ipuknwn=0, iuknwn=0, jpuknwn=0,
    indx=0, indxi=0, indxj=0, icontrol=0, options_skip_groundflow_materidivergence=-NO,
    // groundflow_consolidation_apply family (manual Professional 6.556,
    // 6.613): the material divergence (consolidation coupling) part in the
    // groundflow equation is included only when a switch is set to -yes.
    // Default (record absent) is -no, matching the Professional ("Default
    // switch is -no"). The GNU legacy default -yes coupled every materi +
    // groundflow model into a consolidation transient; the corpus safety
    // tests (ground14/15/16, no consolidation requested) reach the drained
    // steady state within their 1 s window only without the coupling,
    // exactly like the Professional binary does.
    // Measured on the Professional (ground14 A/B): only the GLOBAL or the
    // per-timestep CONTROL record can ACTIVATE the coupling; a group-level
    // -yes alone does not (safety stays at the drained value). The
    // group-level records can only EXCLUDE (-no) the elements of the group
    // from a globally/control-activated coupling (manual 6.613).
    groundflow_consolidation_apply=-NO, control_groundflow_consolidation_apply=-NO,
    group_groundflow_consolidation_apply=-NO,
    materidivergence=-NO, group_materidivergence=-NO,
    dry=0, group_consolidation_found=0,
    global_consolidation_found=0, control_consolidation_found=0,
    group_materidivergence_found=0, ldum=0, idum[1];
  double tmp=0., C=0., dtime=0., divergence=0., dens=0., limit=0., ddum[1], pe[MDIM];

  swit = set_swit(element,-1,"groundflow");
  if ( swit ) pri( "In routine GROUNDFLOW." );

  if ( !db_active_index( GROUP_GROUNDFLOW_PERMEABILITY, gr, VERSION_NORMAL ) &&
       !db_active_index( GROUP_GROUNDFLOW_PERMEABILITY_NONLINEAR_METHOD, gr, VERSION_NORMAL ) &&
       !db_active_index( GROUP_GROUNDFLOW_CAPACITY, gr, VERSION_NORMAL ) &&
       !db_active_index( GROUP_GROUNDFLOW_CAPACITY_NONLINEAR_METHOD, gr, VERSION_NORMAL ) )
    return;
         
  db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET );
  db( GROUNDFLOW_DENSITY, 0, idum, &dens, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  group_materidivergence_found = db( GROUP_GROUNDFLOW_MATERIDIVERGENCE, gr,
    &group_materidivergence, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  group_consolidation_found = db( GROUP_GROUNDFLOW_CONSOLIDATION_APPLY, gr,
    &group_groundflow_consolidation_apply, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  global_consolidation_found = db( GROUNDFLOW_CONSOLIDATION_APPLY, 0,
    &groundflow_consolidation_apply, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( OPTIONS_SKIP_GROUNDFLOW_MATERIDIVERGENCE, 0,
    &options_skip_groundflow_materidivergence, ddum, ldum, VERSION_NORMAL,
    GET_IF_EXISTS );
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  control_consolidation_found = db( CONTROL_GROUNDFLOW_CONSOLIDATION_APPLY,
    icontrol, &control_groundflow_consolidation_apply, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_OPTIONS_SKIP_GROUNDFLOW_MATERIDIVERGENCE, icontrol, 
    &options_skip_groundflow_materidivergence, ddum, ldum, 
    VERSION_NORMAL, GET_IF_EXISTS );
  // Activation: the per-timestep CONTROL record (if present) wins, then the
  // global record; otherwise the coupling stays off (-no default).
  if ( control_consolidation_found )
    materidivergence = control_groundflow_consolidation_apply;
  else if ( global_consolidation_found )
    materidivergence = groundflow_consolidation_apply;
  else
    materidivergence = -NO;
  // Group-level exclusion: the legacy GROUP_GROUNDFLOW_MATERIDIVERGENCE and
  // the new GROUP_GROUNDFLOW_CONSOLIDATION_APPLY records only switch the
  // coupling OFF for the elements of the group (a group -yes cannot activate
  // it - measured on the Professional).
  if ( ( group_consolidation_found &&
         group_groundflow_consolidation_apply==-NO ) ||
       ( group_materidivergence_found && group_materidivergence==-NO ) )
    materidivergence = -NO;
  if ( options_skip_groundflow_materidivergence==-YES ) materidivergence = -NO;

  // groundflow_total_pressure_limit (manual Professional 6.587, default 0):
  // where the calculated total pressure at the integration point reaches
  // the limit the soil is dry and the groundflow velocity is 0 there
  // (dry rule of manual 2.4.3). Measured on the Professional 02-08-2026
  // (F4 T1 probes): the integration-point total pressure is interpolated
  // from the nodal total pressures (already capped to the limit), the
  // velocity is 0 where it is >= the limit, and neither the pres dof nor
  // the consolidation (material divergence) term is affected: with the
  // default limit 0 the free-node pres still develops (1.6667 in the
  // probe), so the former GNU 'limit 0 + pres 0 -> skip the consolidation
  // term' rule is not applied.
  if ( groundflow_velocity && groundflow_pressure ) {
    long int inod_tp=0, found_tp=0;
    double tp_node=0., sp_node=0., loc_node=0., tp_ip=0.,
      *coord_tp=NULL, *dof_tp=NULL;
    limit = groundflow_total_pressure_limit();
    found_tp = 1;
    for ( inol=0; inol<nnol && found_tp; inol++ ) {
      inod_tp = nodes[inol];
      coord_tp = db_dbl( NODE, inod_tp, VERSION_NORMAL );
      dof_tp = db_dbl( NODE_DOF, inod_tp, VERSION_NEW );
      if ( !coord_tp || !dof_tp ) { found_tp = 0; break; }
      tp_node = sp_node = 0.;
      groundflow_phreatic_coord( inod_tp, coord_tp, dof_tp, tp_node,
        sp_node, loc_node, NULL );
      tp_ip += h[inol] * tp_node;
    }
    if ( found_tp &&
         tp_ip >= limit - EPS_TOTAL_PRESSURE_LIMIT*(1.+scalar_dabs(limit)) )
      dry = 1;
  }

  groundflow_data( element, gr, nodes, old_unknowns, new_unknowns, coord_ip, pe, C, h, nnol );

  if ( materi_velocity ) {
    for ( jdim=0; jdim<ndim; jdim++ ) {
      iuknwn = vel_indx + jdim * nder;
      divergence += grad_new_unknowns[jdim*nuknwn+iuknwn];
    }
  }
  if ( swit ) pri( "divergence", divergence );

  for ( inol=0; inol<nnol; inol++ ) {
    if ( groundflow_velocity ) {
      for ( jdim=0; jdim<ndim; jdim++ ) {
        iuknwn = (gvel_indx + jdim*nder);
        ipuknwn = iuknwn/nder;
        indx = inol*npuknwn + ipuknwn;
// Groundflow velocity sign (Professional manual 8.2, "Groundflow
// velocities"): v_i^g = k_ip * dh/dx_i, i.e. the velocity follows from
// the permeability times the POSITIVE head gradient. The 2013 F. Lorenzo
// change negated this term, which inverts the convection direction in the
// coupled condif+groundflow tests (condif5: the middle node was driven to
// the downstream temperature 0 instead of the Professional's upwind 1).
// Restored to the Professional sign so that the heat convection follows
// the same upwind direction as the Professional.
// C85 (F4 T6f, 2026-10-17, owner decision; DELIBERATE DIVERGENCE
// TRIA6-HRZ): the tria6 lumps the velocity with the HRZ weight of pol(),
// the same weight as its general.cc denominator (with h the tria6 corner
// denominator sum volume h is 0 and the corner velocity is 0/0). Every
// other element keeps h.
        double h_lump = lumped_hrz_active ? lumped_hrz_weight[inol] : h[inol];
        if ( dry )
          tmp = h_lump * ( 0. - old_unknowns[iuknwn] ) / dtime;
        else
          tmp = h_lump * ( pe[jdim]*grad_new_unknowns[jdim*nuknwn+pres_indx] -
            old_unknowns[iuknwn] ) / dtime;
        element_rhside[indx] += volume * tmp;
      }
    }
    if ( materidivergence==-YES ) {
      if ( materi_velocity ) {
        ipuknwn = pres_indx/nder;
        indx = inol*npuknwn + ipuknwn;
        tmp = h[inol] * divergence;
        element_rhside[indx] += volume * tmp;
        if ( residue ) element_residue[indx] += tmp;
        indxi = inol*npuknwn + ipuknwn;
        for ( jnol=0; jnol<nnol; jnol++ ) {
          for ( jdim=0; jdim<ndim; jdim++ ) {
            jpuknwn = vel_indx/nder + jdim;
            indxj = jnol*npuknwn + jpuknwn;
            element_matrix[indxi*nnol*npuknwn+indxj] -=
              volume * h[inol] * d[jdim*nnol+jnol];
          }
        }
      }
    }
  }
  if ( swit ) {
    pri( "element_rhside", element_rhside, nnol, npuknwn );
  }

  if ( swit ) pri( "Out function GROUNDFLOW" );

}

long int groundflow_phreatic_level_multiple_active( void )

// True when at least one groundflow_phreatic_level_multiple record exists
// (any index). Each groundwater level has its own index (ground8/ground19
// of the corpus store them at 10/20/30), so testing a single index is not
// enough: db_active_index(...,0,...) only sees a record stored at index 0
// (the internal groundflow_phreatic_multiple suite test stores its two
// levels at 0 and 1).

{
  long int max_multiple=0;
  db_max_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE, max_multiple,
    VERSION_NORMAL, GET );
  return max_multiple>=0;
}

long int groundflow_phreatic_level_multiple_find_element( long int elnum )

// Returns the index (imult) of the groundflow_phreatic_level_multiple record
// whose domain owns the given ELEMENT, or -1 if none. The domain of each
// record is selected by exactly one of _element, _element_group,
// _element_geometry or _node (they are not combinable, per the manual). A
// _node domain cannot own an element. Levels are scanned in ascending
// index and the FIRST (lowest) matching level is returned.

{
  long int imult=0, max_multiple=-1, found=-1, length=0, ldum=0,
    gr=0, j=0,
    length_list=0, use_element=0, use_group=0, use_geometry=0, use_node=0,
    itmp=0, all=0, any=0;
  double ddum[MDIM], rdum=0.; // geometry() writes normal/projection[MDIM]
  long int *sel=NULL, *nodes=NULL, *el=NULL;

  db_max_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE, max_multiple, VERSION_NORMAL, GET );
  if ( max_multiple<0 ) return -1;

  nodes = get_new_int(MNOL);
  el = get_new_int(MNOL+1);
  for ( imult=0; imult<=max_multiple && found<0; imult++ ) {
    if ( !db_active_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE, imult, VERSION_NORMAL ) )
      continue;
    use_element = use_group = use_geometry = use_node = 0;
    if ( db_active_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT,
        imult, VERSION_NORMAL ) ) {
      use_element = 1;
      sel = db_int( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT, imult, VERSION_NORMAL );
      length_list = db_len( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT, imult, VERSION_NORMAL );
    }
    else if ( db_active_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GROUP,
        imult, VERSION_NORMAL ) ) {
      use_group = 1;
      sel = db_int( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GROUP, imult, VERSION_NORMAL );
      length_list = db_len( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GROUP,
        imult, VERSION_NORMAL );
    }
    else if ( db_active_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GEOMETRY,
        imult, VERSION_NORMAL ) ) {
      use_geometry = 1;
      sel = db_int( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GEOMETRY,
        imult, VERSION_NORMAL );
      length_list = db_len( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GEOMETRY,
        imult, VERSION_NORMAL );
    }
    else if ( db_active_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_NODE,
        imult, VERSION_NORMAL ) ) {
      use_node = 1;
      sel = db_int( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_NODE, imult, VERSION_NORMAL );
      length_list = db_len( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_NODE,
        imult, VERSION_NORMAL );
    }
    if ( use_node ) continue;
    if ( use_element ) {
      if ( array_member( sel, elnum, length_list, ldum ) ) found = imult;
    }
    else if ( use_group ) {
      gr = 0;
      db( ELEMENT_GROUP, elnum, &gr, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
      if ( array_member( sel, gr, length_list, ldum ) ) found = imult;
    }
    else if ( use_geometry ) {
      db( ELEMENT, elnum, el, ddum, length, VERSION_NORMAL, GET );
      long int nnol_el = length - 1;
      for ( j=1; j<=nnol_el; j++ ) nodes[j-1] = el[j];
      all = 1; any = 0;
      for ( j=0; j<nnol_el; j++ ) {
        geometry( nodes[j], ddum, sel, itmp, rdum, ddum, rdum,
          ddum, NODE_START_REFINED, PROJECT_EXACT, VERSION_NORMAL );
        if ( !itmp ) all = 0;
        if ( itmp ) any = 1;
      }
      if ( all || any ) found = imult;
    }
  }
  delete[] nodes;
  delete[] el;

  return found;
}

long int groundflow_phreatic_level_multiple_find( long int inod )

// Returns the index (imult) of the groundflow_phreatic_level_multiple record
// that owns the given node, or -1 if none. The node belongs to the domain of
// a record if it is listed in its _node list, or if one of its elements is
// in the _element/_element_group/_element_geometry domain (see
// groundflow_phreatic_level_multiple_find_element). When several records
// match, the LOWEST index wins (the level scan is ascending).

{
  long int imult=0, max_multiple=-1, found=-1, nel=0, iel=0, elnum=0, fe=0;
  long int *node_element=NULL;

  db_max_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE, max_multiple, VERSION_NORMAL, GET );
  if ( max_multiple<0 ) return -1;

  for ( imult=0; imult<=max_multiple; imult++ ) {
    if ( !db_active_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE, imult, VERSION_NORMAL ) )
      continue;
    if ( db_active_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_NODE, imult,
        VERSION_NORMAL ) ) {
      long int ldum=0;
      long int length_list = db_len( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_NODE,
        imult, VERSION_NORMAL );
      long int *sel = db_int( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_NODE,
        imult, VERSION_NORMAL );
      if ( array_member( sel, inod, length_list, ldum ) ) return imult;
    }
  }

  // nodes without a NODE_ELEMENT record (e.g. macro-generated meshes where
  // the record is only stored for nodes attached to elements) cannot
  // belong to any element-based domain
  if ( !db_active_index( NODE_ELEMENT, inod, VERSION_NORMAL ) ) return -1;
  node_element = db_int( NODE_ELEMENT, inod, VERSION_NORMAL );
  nel = db_len( NODE_ELEMENT, inod, VERSION_NORMAL );
  for ( iel=0; iel<nel; iel++ ) {
    elnum = node_element[iel];
    fe = groundflow_phreatic_level_multiple_find_element( elnum );
    if ( fe>=0 && ( found<0 || fe<found ) ) found = fe;
  }

  return found;
}

long int groundflow_phreatic_level_multiple_find_coord( double coord[] )

// Returns the index (imult) of the multiple level owning the mesh element
// that contains the given coordinate, or -1 when the coordinate lies
// outside the mesh or the containing element belongs to no level domain.
// The element is located like the post-point machinery does (point_el over
// all active elements with the NODE_START_REFINED frame). Needed for the
// inod<0 evaluations of groundflow_phreatic_coord (post points), where no
// node membership exists.

{
  long int element=0, max_element=0, length=0, inol=0, nnol=0, inod=0,
    ldum=0, idum[1], el[1+MNOL], found_el=-1;
  double ddum[1], coords[MNOL*MDIM], weight[MNOL];

  db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
  for ( element=0; element<=max_element && found_el<0; element++ ) {
    if ( !db_active_index( ELEMENT, element, VERSION_NORMAL ) ) continue;
    db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
    nnol = length - 1;
    for ( inol=0; inol<nnol; inol++ ) {
      inod = el[inol+1];
      db( NODE_START_REFINED, inod, idum, &coords[inol*ndim], ldum,
        VERSION_NORMAL, GET );
    }
    if ( point_el( coord, coords, weight, el[0], nnol ) )
      found_el = element;
  }
  if ( found_el<0 ) return -1;

  return groundflow_phreatic_level_multiple_find_element( found_el );
}

long int groundflow_phreatic_level_reference( long int inod, double coord[],
  double &z_ref )

// Head reference of a groundflow phreatic level (F5 T5 U3b/U3c).
// The pressure dof is the head minus rho*g*z_ref, with z_ref the height
// of the first point of the level table (the level itself when it is
// one value). Under a horizontal level z_ref = z_L, so the dof is the
// dynamic pressure p_total - rho*g*(z_L - z) of the GNU split. Under an
// inclined level a single constant keeps the dof a shifted head, so the
// Darcy flow follows the head gradient: measured on the Professional
// 02-08-2026, the nodes at and above an inclined level get the head
// rho*g*z_L(x) and water seeps through the domain, for a single
// groundflow_phreatic_level and for a groundflow_phreatic_level_multiple
// alike. For multiple levels the record that owns the node (inod>=0) or
// the coordinate (inod<0) gives z_ref. Returns 0 (z_ref untouched) when
// no level applies.

{
  long int length=0, imult=-1, idat=GROUNDFLOW_PHREATICLEVEL, index=0;
  double *level=NULL;

  if ( groundflow_phreatic_level_multiple_active() ) {
    if ( inod>=0 )
      imult = groundflow_phreatic_level_multiple_find( inod );
    else
      imult = groundflow_phreatic_level_multiple_find_coord( coord );
    if ( imult<0 ) return 0;
    // _static -yes records prescribe the local static pressure instead
    // (groundflow_phreatic_apply); their reference stays the local level
    if ( db_active_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_STATIC, imult,
         VERSION_NORMAL ) &&
         db_int( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_STATIC, imult,
           VERSION_NORMAL )[0]==-YES ) return 0;
    idat = GROUNDFLOW_PHREATICLEVEL_MULTIPLE;
    index = imult;
  }
  if ( !db_active_index( idat, index, VERSION_NORMAL ) ) return 0;
  length = db_len( idat, index, VERSION_NORMAL );
  level = db_dbl( idat, index, VERSION_NORMAL );
  if      ( length==1 ) z_ref = level[0];
  else if ( ndim==2 && length>=2 ) z_ref = level[1];
  else if ( ndim==3 && length>=3 ) z_ref = level[2];
  else db_error( idat, index );
  return 1;
}

long int groundflow_phreatic_coord( long int inod, double coord[], double dof[],
  double &total_pressure, double &static_pressure, double &location,
  long int *level_source )

// inod only in arguments for test printing. When level_source is not
// NULL it receives the source of the static pressure at the node:
// 0 = no level found, 1 = groundflow_phreatic_level(_multiple),
// 2 = post_calcul_static_pressure_height region. Callers that convert
// a prescribed TOTAL pressure into the pres dof (bounda_dof -topres)
// need the distinction: the head-to-dof conversion differs between the
// phreatic level (dynamic pressure uniform below the level) and the
// static-pressure-height reference (per-node static subtraction).

{
  long int length=0, found=0, ldum=0, idum[1], number[2], imult=0,
    static_switch=-NO;
  double water_level=0., dens=0., pressure_atmospheric=0., addtopressure=0.,
    total_pressure_limit=0.,
    ddum[1], force_gravity[MDIM], *groundflow_phreatic=NULL;

  if ( level_source ) *level_source = 0;
  force_gravity_calculate( force_gravity );
  db( GROUNDFLOW_DENSITY, 0, idum, &dens, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( GROUNDFLOW_PRESSURE_ATMOSPHERIC, 0, idum, &pressure_atmospheric, 
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  // The calculated total pressure is capped by groundflow_total_pressure_limit
  // (manual Professional 6.587, default 0) at every node, prescribed ones
  // included; the GNU-only groundflow_pressure_atmospheric caps the static
  // pressure only.
  total_pressure_limit = groundflow_total_pressure_limit();
  db( GROUNDFLOW_ADDTOPRESSURE, 0, idum, &addtopressure, ldum, VERSION_NORMAL, GET_IF_EXISTS );

  if ( groundflow_pressure ) {
    total_pressure = dof[pres_indx] - dens * force_gravity[ndim-1] * coord[ndim-1];
  }

  // groundflow_phreatic_level_multiple: several groundwater levels, each
  // owning a part of the domain (selected by _element/_element_group/
  // _element_geometry/_node). The node uses the level of its owning record.
  // For inod<0 (post-point evaluations, which carry no node membership) the
  // owning level is resolved by locating the element that contains the
  // coordinate (groundflow_phreatic_level_multiple_find_coord).
  if ( groundflow_phreatic_level_multiple_active() ) {
    if ( inod>=0 )
      imult = groundflow_phreatic_level_multiple_find( inod );
    else
      imult = groundflow_phreatic_level_multiple_find_coord( coord );
    if ( imult>=0 ) {
      static_switch = -NO;
      if ( db_active_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_STATIC, imult,
          VERSION_NORMAL ) ) {
        db( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_STATIC, imult, &static_switch,
          ddum, ldum, VERSION_NORMAL, GET );
      }
      length = db_len( GROUNDFLOW_PHREATICLEVEL_MULTIPLE, imult, VERSION_NORMAL );
      groundflow_phreatic = db_dbl( GROUNDFLOW_PHREATICLEVEL_MULTIPLE, imult, VERSION_NORMAL );
      if      ( ndim==1 ) {
        if ( length!=1 ) db_error( GROUNDFLOW_PHREATICLEVEL_MULTIPLE, imult );
        water_level = groundflow_phreatic[0];
        found = 1;
      }
      else if ( ndim==2 ) {
        if ( length==1 ) {
          found = 1;
          water_level = groundflow_phreatic[0];
        }
        else {
          found = table_xy( groundflow_phreatic, "GROUNDFLOW_PHREATICLEVEL_MULTIPLE",
            length, coord[0], water_level );
        }
      }
      else {
        assert( ndim==3 );
        if ( length==1 ) {
          found = 1;
          water_level = groundflow_phreatic[0];
        }
        else {
          db( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_N, imult, number, ddum, ldum,
            VERSION_NORMAL, GET );
          if ( number[0]*number[1]*3 != length )
            db_error( GROUNDFLOW_PHREATICLEVEL_MULTIPLE, imult );
          found = table_xyz( groundflow_phreatic, number, coord, water_level );
        }
      }
      if ( found ) {
        location = water_level;
        if ( level_source ) *level_source = 1;
        if ( groundflow_pressure ) {
          double z_ref = water_level;
          groundflow_phreatic_level_reference( inod, coord, z_ref );
          static_pressure =
            force_gravity[ndim-1] * dens * ( water_level - coord[ndim-1] );
          // the dof is the head minus rho*g*z_ref of the record (F5 T5 U3c)
          total_pressure = dof[pres_indx] +
            force_gravity[ndim-1] * dens * ( z_ref - coord[ndim-1] );
        }
        else if ( static_switch==-YES ) {
          // groundflow_phreatic_level_multiple_static without a
          // pressure unknown (mechanics-only static groundflow, manual
          // Professional 6.574 + 6.573): the total (pore) pressure of
          // the level is the static pressure itself
          // rho*g*(z_level - z), clamped to the atmospheric pressure
          // above the level by the common cap below the branches.
          static_pressure =
            force_gravity[ndim-1] * dens * ( water_level - coord[ndim-1] );
          total_pressure = static_pressure;
        }
      }
    }
  }
  else if ( db_active_index( GROUNDFLOW_PHREATICLEVEL, 0, VERSION_NORMAL ) ) {

    length = db_len( GROUNDFLOW_PHREATICLEVEL, 0, VERSION_NORMAL );
    groundflow_phreatic = db_dbl( GROUNDFLOW_PHREATICLEVEL, 0, VERSION_NORMAL ); 
    static_switch = -NO;
    if ( db_active_index( GROUNDFLOW_PHREATICLEVEL_STATIC, 0, VERSION_NORMAL ) ) {
      db( GROUNDFLOW_PHREATICLEVEL_STATIC, 0, &static_switch, ddum, ldum,
        VERSION_NORMAL, GET );
    }

    if      ( ndim==1 ) {
      if ( length!=1 ) db_error( GROUNDFLOW_PHREATICLEVEL, 0 );
      water_level = groundflow_phreatic[0];
      found = 1;
    }
    else if ( ndim==2 ) {
      if ( length==1 ) {
        found = 1;
        water_level = groundflow_phreatic[0];
      }
      else {
        found = table_xy( groundflow_phreatic, "GROUNDFLOW_PHREATICLEVEL",
          length, coord[0], water_level );
      }
    }
    else {
      assert( ndim==3 );
      if ( length==1 ) {
        found = 1;
        water_level = groundflow_phreatic[0];
      }
      else {
        db( GROUNDFLOW_PHREATICLEVEL_N, 0, number, ddum, ldum, VERSION_NORMAL, GET );
        if ( number[0]*number[1]*3 != length ) db_error( GROUNDFLOW_PHREATICLEVEL, 0 );
        found = table_xyz( groundflow_phreatic, number, coord, water_level );
      }
    }
    if ( found ) {
      location = water_level;
      if ( level_source ) *level_source = 1;
      if ( groundflow_pressure ) {
        double z_ref = water_level;
        groundflow_phreatic_level_reference( inod, coord, z_ref );
        static_pressure =
          force_gravity[ndim-1] * dens * ( water_level - coord[ndim-1] );
        // the dof is the head minus rho*g*z_ref (inclined level: z_ref
        // is not the local level, see groundflow_phreatic_level_reference)
        total_pressure = dof[pres_indx] +
          force_gravity[ndim-1] * dens * ( z_ref - coord[ndim-1] );
      }
      else if ( static_switch==-YES ) {
        // groundflow_phreatic_level_static without a pressure unknown
        // (corpus ground18/excavate1, manual Professional 6.573: the
        // total pressures of the nodes for which the level holds are
        // set equal to the static pressure; the phreatic line may lie
        // above the mesh part and the hydraulic heads are not solved).
        static_pressure =
          force_gravity[ndim-1] * dens * ( water_level - coord[ndim-1] );
        total_pressure = static_pressure;
      }
    }
  }

  if ( static_pressure>=pressure_atmospheric ) static_pressure = pressure_atmospheric;
  if ( total_pressure>=total_pressure_limit ) total_pressure = total_pressure_limit;
  total_pressure += addtopressure;

  // post_calcul_static_pressure_height (manual Professional 6.921/6.922):
  // when no groundwater level applies to the node, the static pressure
  // is determined relative to a reference height instead: regions of
  // the vertical coordinate (coord_min, coord_max) each with their own
  // height_ref; p_static = rho*g*(height_ref - coord_vertical). The
  // post_calcul_static_pressure_height_element_group record restricts
  // each region to an element group (-all = every group). total =
  // pres_dof + static like in the phreatic-level branch.
  if ( !found && groundflow_pressure &&
       db_active_index( POST_CALCUL_STATIC_PRESSURE_HEIGHT, 0,
         VERSION_NORMAL ) ) {
    long int iregion=0, height_length=0, ngroups=0, igroup_ok=0, k=0,
      node_groups[DATA_ITEM_SIZE], group_allowed=-ALL, nnode_groups=0,
      *group_ival=NULL;
    double *height_rec = NULL;
    height_rec = db_dbl( POST_CALCUL_STATIC_PRESSURE_HEIGHT, 0,
      VERSION_NORMAL );
    height_length = db_len( POST_CALCUL_STATIC_PRESSURE_HEIGHT, 0,
      VERSION_NORMAL );
    if ( height_length%3!=0 )
      db_error( POST_CALCUL_STATIC_PRESSURE_HEIGHT, 0 );
    if ( db_active_index( POST_CALCUL_STATIC_PRESSURE_HEIGHT_ELEMENT_GROUP,
         0, VERSION_NORMAL ) ) {
      group_ival = db_int( POST_CALCUL_STATIC_PRESSURE_HEIGHT_ELEMENT_GROUP,
        0, VERSION_NORMAL );
      ngroups = db_len( POST_CALCUL_STATIC_PRESSURE_HEIGHT_ELEMENT_GROUP,
        0, VERSION_NORMAL );
    }
    for ( iregion=0; iregion*3+2<height_length && !found; iregion++ ) {
      if ( coord[ndim-1]>=height_rec[iregion*3+0]-EPS_COORD &&
           coord[ndim-1]<=height_rec[iregion*3+1]+EPS_COORD ) {
        group_allowed = -ALL;
        if ( ngroups>0 ) {
          if ( iregion>=ngroups )
            db_error( POST_CALCUL_STATIC_PRESSURE_HEIGHT_ELEMENT_GROUP, 0 );
          group_allowed = group_ival[iregion];
        }
        igroup_ok = ( group_allowed==-ALL );
        if ( !igroup_ok ) {
          nnode_groups = 0;
          node_attached_element_groups( inod, node_groups, nnode_groups );
          for ( k=0; k<nnode_groups; k++ ) {
            if ( node_groups[k]==group_allowed ) { igroup_ok = 1; break; }
          }
        }
        if ( igroup_ok ) {
          location = height_rec[iregion*3+2];
          if ( level_source ) *level_source = 2;
          static_pressure = force_gravity[ndim-1] * dens *
            ( location - coord[ndim-1] );
          total_pressure = dof[pres_indx] + static_pressure;
          found = 1;
        }
      }
    }
  }
  if ( static_pressure>=pressure_atmospheric ) static_pressure = pressure_atmospheric;
  if ( total_pressure>=total_pressure_limit ) total_pressure = total_pressure_limit;
  total_pressure += addtopressure;

  return found;
}

long int groundflow_static_pressure_group( long int element_group )

// Mechanics-only static groundflow (manual Professional 6.574
// groundflow_phreatic_level_static): without a groundflow_pressure
// unknown, the elements of a group whose group_type holds -groundflow
// get the static pore pressure of the phreatic level. Returns 1 when
// that applies to element_group: no pressure unknown, a static level
// switch set to -yes, and -groundflow in the group_type. Both the pore
// pressure on the skeleton (materi.cc) and the wet/dry density choice
// (get_materi_density) use it. Measured on the Professional 02-08-2026
// (ground18 of the corpus with group_type -materi only): without
// -groundflow the group keeps the wet density everywhere and gets no
// pore pressure (bottom sigyy -1995, like the uncorrected GNU).

{
  long int itype=0, len_type=0, types[DATA_ITEM_SIZE];
  double ddum[1];

  if ( groundflow_pressure ) return 0;
  if ( !groundflow_static_pressure_active() ) return 0;
  db( GROUP_TYPE, element_group, types, ddum, len_type,
    VERSION_NORMAL, GET_IF_EXISTS );
  for ( itype=0; itype<len_type; itype++ ) {
    if ( types[itype]==-GROUNDFLOW ) return 1;
  }
  return 0;
}

static double groundflow_phreatic_project_height( double level[],
  long int length, double coord[] )

// groundflow_phreatic_project -yes (manual 6.582): height of the point of
// the 2D level line (x z pairs) nearest to coord, each segment clamped to
// its end points. Measured on the Professional 02-08-2026: node (0,3)
// under the line (0,1)-(2,2.5) gets the head of (0.96,1.72).

{
  long int i=0;
  double best_dist=-1., best_z=level[1], dx=0., dz=0., len2=0., t=0.,
    px=0., pz=0., dist=0.;

  for ( i=0; i+3<length; i+=2 ) {
    dx = level[i+2] - level[i];
    dz = level[i+3] - level[i+1];
    len2 = dx*dx + dz*dz;
    t = 0.;
    if ( len2>0. )
      t = ( ( coord[0]-level[i] )*dx + ( coord[1]-level[i+1] )*dz ) / len2;
    if ( t<0. ) t = 0.;
    if ( t>1. ) t = 1.;
    px = level[i] + t*dx;
    pz = level[i+1] + t*dz;
    dist = ( coord[0]-px )*( coord[0]-px ) + ( coord[1]-pz )*( coord[1]-pz );
    if ( best_dist<0. || dist<best_dist ) {
      best_dist = dist;
      best_z = pz;
    }
  }
  return best_z;
}

void groundflow_phreatic_apply( void )

{
  long int inod=0, max_node=0, iuknwn=0, ipuknwn=0, 
    groundflow_phreaticlevel_bounda=-YES, length=0, 
    node_phreaticlevel=0, ldum=0, idum[1], *node_bounded=NULL,
    imult2=0, static_switch=-NO;
  double total_pressure=0, static_pressure=0., dens=0., location=0.,
    ddum[1], force_gravity[MDIM], *coord=NULL, *node_dof=NULL;
 
  // groundflow_phreatic_bounda (manual 6.570) is a switch: -yes (the
  // Professional default, also when the record is absent) is the
  // free-surface condition of the single level below; -no switches it
  // off. -method1 / -method2 are the older GNU variants (F5 T5 U1b).
  if ( db_active_index( GROUNDFLOW_PHREATICLEVEL_BOUNDA, 0, VERSION_NORMAL ) )
    db( GROUNDFLOW_PHREATICLEVEL_BOUNDA, 0, &groundflow_phreaticlevel_bounda,
      ddum, ldum, VERSION_NORMAL, GET );
  if ( groundflow_phreaticlevel_bounda==-METHOD1 ||
       groundflow_phreaticlevel_bounda==-METHOD2 ) {
    db( GROUNDFLOW_DENSITY, 0, idum, &dens, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    force_gravity_calculate( force_gravity );
    db_max_index( NODE, max_node, VERSION_NORMAL, GET );
    for ( inod=0; inod<=max_node; inod++ ) {
      if ( db_active_index( NODE_START_REFINED, inod, VERSION_NORMAL ) ) {
        coord = db_dbl( NODE_START_REFINED, inod, VERSION_NORMAL );
        node_dof = db_dbl( NODE_DOF, inod, VERSION_NEW );
        if ( groundflow_phreatic_coord( inod, coord, node_dof, total_pressure,
            static_pressure, location, NULL ) ) {
          iuknwn = pres_indx;
          ipuknwn = iuknwn / nder;
          if ( groundflow_phreaticlevel_bounda==-METHOD1 ) {
            if      ( scalar_dabs(coord[ndim-1]-location)<EPS_COORD ) {
              node_bounded = db_int( NODE_BOUNDED, inod, VERSION_NORMAL );
              node_dof[iuknwn] = dens * force_gravity[ndim-1] * location;
              node_bounded[ipuknwn] = 1;
            }
          }
          else if ( groundflow_phreaticlevel_bounda==-METHOD2 ) {
            if ( coord[ndim-1]>(location-EPS_COORD) ) {
              node_bounded = db_int( NODE_BOUNDED, inod, VERSION_NORMAL );
              node_dof[iuknwn] = dens * force_gravity[ndim-1] * location;
              node_bounded[ipuknwn] = 1;
            }
          }
          else {
            db_error( GROUNDFLOW_PHREATICLEVEL_BOUNDA, -1 );
          }
          if ( coord[ndim-1]>location+EPS_COORD ) {
            node_phreaticlevel = -ABOVE;
          }
          else {
            node_phreaticlevel = -BELOW;
          }
          length = 1;
          db( NODE_PHREATICLEVEL, inod, &node_phreaticlevel, ddum, 
            length, VERSION_NORMAL, PUT );
        }
      }
    }
  }

  // Single groundflow_phreatic_level: free-surface condition of the
  // saturated zone. Nodes at or above the phreatic line are dry and
  // carry no water pressure (p_total = 0). In the GNU pressure
  // convention (p_total = pres_dof + static_pressure, with the static
  // part clamped to the atmospheric pressure 0 above the level) that
  // bounds the pres dof to 0 there, which confines the saturated flow
  // domain below the level. Measured against the Professional
  // (ground15/16 of the corpus): its hydraulic head on/above the
  // phreatic line is h = rho*g*z_L (p_total = 0), i.e. p_dynamic = 0 in
  // the GNU split. Explicit bounda_dof records win over this default
  // (bounda() applies them afterwards).
  // Inclined level (F5 T5 U3b, measured on the Professional 02-08-2026):
  // the head is rho*g*z_L(x) of the level above the node, so the dof is
  // rho*g*(z_L(x) - z_ref) (groundflow_phreatic_level_reference). With
  // groundflow_phreatic_project -yes (manual 6.582) z_L is the height of
  // the nearest point of the 2D level line, clamped to its end points.
  if ( groundflow_pressure && groundflow_phreaticlevel_bounda!=-NO &&
       db_active_index( GROUNDFLOW_PHREATICLEVEL, 0, VERSION_NORMAL ) &&
       !groundflow_phreatic_level_multiple_active() ) {
    long int level_len=0, gfound=0, number2[2], project=-NO;
    double water_level2=0., z_ref=0., z_head=0., *groundflow_phreatic2=NULL;
    level_len = db_len( GROUNDFLOW_PHREATICLEVEL, 0, VERSION_NORMAL );
    groundflow_phreatic2 = db_dbl( GROUNDFLOW_PHREATICLEVEL, 0,
      VERSION_NORMAL );
    db( GROUNDFLOW_DENSITY, 0, idum, &dens, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    force_gravity_calculate( force_gravity );
    if ( db_active_index( GROUNDFLOW_PHREATIC_PROJECT, 0, VERSION_NORMAL ) ) {
      db( GROUNDFLOW_PHREATIC_PROJECT, 0, &project, ddum, ldum,
        VERSION_NORMAL, GET );
      if ( project!=-YES && project!=-NO )
        db_error( GROUNDFLOW_PHREATIC_PROJECT, 0 );
      if ( project==-YES && ndim==3 && level_len>1 ) {
        pri( "Error: groundflow_phreatic_project -yes needs a 2D level line." );
        exit(TN_EXIT_STATUS);
      }
    }
    db_max_index( NODE, max_node, VERSION_NORMAL, GET );
    for ( inod=0; inod<=max_node; inod++ ) {
      if ( db_active_index( NODE_START_REFINED, inod, VERSION_NORMAL ) ) {
        coord = db_dbl( NODE_START_REFINED, inod, VERSION_NORMAL );
        if      ( ndim==1 ) {
          if ( level_len!=1 ) db_error( GROUNDFLOW_PHREATICLEVEL, 0 );
          water_level2 = groundflow_phreatic2[0];
          gfound = 1;
        }
        else if ( ndim==2 ) {
          if ( level_len==1 ) {
            water_level2 = groundflow_phreatic2[0];
            gfound = 1;
          }
          else {
            gfound = table_xy( groundflow_phreatic2,
              "GROUNDFLOW_PHREATICLEVEL", level_len, coord[0],
              water_level2 );
          }
        }
        else {
          assert( ndim==3 );
          if ( level_len==1 ) {
            water_level2 = groundflow_phreatic2[0];
            gfound = 1;
          }
          else {
            db( GROUNDFLOW_PHREATICLEVEL_N, 0, number2, ddum, ldum,
              VERSION_NORMAL, GET );
            if ( number2[0]*number2[1]*3 != level_len )
              db_error( GROUNDFLOW_PHREATICLEVEL, 0 );
            gfound = table_xyz( groundflow_phreatic2, number2, coord,
              water_level2 );
          }
        }
        if ( gfound && coord[ndim-1]>=water_level2-EPS_COORD ) {
          z_head = water_level2;
          if ( project==-YES && ndim==2 && level_len>=4 )
            z_head = groundflow_phreatic_project_height( groundflow_phreatic2,
              level_len, coord );
          z_ref = water_level2;
          groundflow_phreatic_level_reference( inod, coord, z_ref );
          node_dof = db_dbl( NODE_DOF, inod, VERSION_NEW );
          iuknwn = pres_indx;
          ipuknwn = iuknwn / nder;
            // + 0.: a zero height difference times negative gravity must
            // print 0, not -0 (output unchanged for a horizontal level)
          node_dof[iuknwn] = dens * force_gravity[ndim-1] * ( z_head - z_ref ) + 0.;
          node_bounded = db_int( NODE_BOUNDED, inod, VERSION_NORMAL );
          node_bounded[ipuknwn] = 1;
        }
      }
    }
  }

  // groundflow_phreatic_level_multiple without _static -yes: same
  // free-surface condition as the single groundflow_phreatic_level, but
  // per level and restricted to the nodes of its own domain (selected by
  // _element/_element_group/_element_geometry/_node). Nodes of the domain
  // at or above the phreatic line are dry: their pres dof is bounded to 0
  // (p_total = 0 there once the static part is clamped to the atmospheric
  // pressure), which confines the saturated flow domain below the level.
  // Measured against the Professional binary (ground8 of the corpus, .dbs
  // 25-10-2023): below each level the hydraulic head is uniform
  // h = rho*g*z_L, i.e. p_dynamic = 0 in the GNU split, and the total
  // pressure is the hydrostatic profile rho*g*(z_L - z) capped to 0 above
  // the level. Levels with _static -yes are excluded: the block below
  // prescribes the static pressure over the whole domain instead.
  // Inclined multiple level (F5 T5 U3c, measured on the Professional
  // 02-08-2026, same result as the single level): the dry node gets the
  // head rho*g*z_L(x), i.e. pres_dof = rho*g*(z_L(x) - z_ref) with z_ref
  // the first height of its record (groundflow_phreatic_level_reference).
  if ( groundflow_pressure && groundflow_phreatic_level_multiple_active() ) {
    long int static_switch3=-NO;
    double z_ref3=0.;
    db( GROUNDFLOW_DENSITY, 0, idum, &dens, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    force_gravity_calculate( force_gravity );
    db_max_index( NODE, max_node, VERSION_NORMAL, GET );
    for ( inod=0; inod<=max_node; inod++ ) {
      if ( !db_active_index( NODE_START_REFINED, inod, VERSION_NORMAL ) )
        continue;
      imult2 = groundflow_phreatic_level_multiple_find( inod );
      if ( imult2<0 ) continue;
      static_switch3 = -NO;
      if ( db_active_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_STATIC, imult2,
          VERSION_NORMAL ) ) {
        db( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_STATIC, imult2, &static_switch3,
          ddum, ldum, VERSION_NORMAL, GET );
      }
      if ( static_switch3==-YES ) continue;
      coord = db_dbl( NODE_START_REFINED, inod, VERSION_NORMAL );
      node_dof = db_dbl( NODE_DOF, inod, VERSION_NEW );
      if ( groundflow_phreatic_coord( inod, coord, node_dof, total_pressure,
          static_pressure, location, NULL ) &&
           coord[ndim-1]>=location-EPS_COORD ) {
        z_ref3 = location;
        groundflow_phreatic_level_reference( inod, coord, z_ref3 );
        iuknwn = pres_indx;
        ipuknwn = iuknwn / nder;
          // + 0.: no -0 for a zero height difference (see above)
        node_dof[iuknwn] = dens * force_gravity[ndim-1] * ( location - z_ref3 ) + 0.;
        node_bounded = db_int( NODE_BOUNDED, inod, VERSION_NORMAL );
        node_bounded[ipuknwn] = 1;
      }
    }
  }

  // groundflow_phreatic_level_multiple_static: for the nodes of a multiple
  // phreatic level with _static -yes, set the total pressure (pore pressure)
  // equal to the static pressure. Convenient when the phreatic line is located
  // above the mesh part to which it belongs (no boundary condition can be
  // imposed), and to avoid solving the hydraulic heads (saves memory and CPU).
  if ( groundflow_phreatic_level_multiple_active() ) {
    long int max_multiple=0, iuknwn=0, ipuknwn=0;
    db_max_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE, max_multiple, VERSION_NORMAL, GET );
    db( GROUNDFLOW_DENSITY, 0, idum, &dens, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    force_gravity_calculate( force_gravity );
    for ( imult2=0; imult2<=max_multiple; imult2++ ) {
      if ( !db_active_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE, imult2, VERSION_NORMAL ) )
        continue;
      if ( !db_active_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_STATIC, imult2,
          VERSION_NORMAL ) )
        continue;
      db( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_STATIC, imult2, &static_switch,
        ddum, ldum, VERSION_NORMAL, GET );
      if ( static_switch==-YES ) {
        db_max_index( NODE, max_node, VERSION_NORMAL, GET );
        for ( inod=0; inod<=max_node; inod++ ) {
          if ( db_active_index( NODE_START_REFINED, inod, VERSION_NORMAL ) ) {
            if ( groundflow_phreatic_level_multiple_find( inod )==imult2 ) {
              coord = db_dbl( NODE_START_REFINED, inod, VERSION_NORMAL );
              node_dof = db_dbl( NODE_DOF, inod, VERSION_NEW );
              if ( groundflow_phreatic_coord( inod, coord, node_dof, total_pressure,
                  static_pressure, location, NULL ) ) {
                iuknwn = pres_indx;
                ipuknwn = iuknwn / nder;
                node_dof[iuknwn] = static_pressure;
                node_bounded = db_int( NODE_BOUNDED, inod, VERSION_NORMAL );
                node_bounded[ipuknwn] = 1;
              }
            }
          }
        }
      }
    }
  }
  
}

double groundflow_total_pressure_limit( void )

// groundflow_total_pressure_limit (manual Professional 6.587): the maximum
// allowed total pressure. Any higher calculated total pressure is cut off
// to this value; "Default the limit is set to 0". Single source of truth of
// the limit and its default: the cap of the calculated total pressure in
// groundflow_phreatic_coord() and the dry rule (groundflow velocity 0) in
// groundflow() both read it here. The pres dof itself is never cut
// (measured on the Professional 02-08-2026, F4 T1).

{
  long int ldum=0, idum[1];
  double limit=GROUNDFLOW_TOTAL_PRESSURE_LIMIT_DEFAULT;

  db( GROUNDFLOW_TOTAL_PRESSURE_LIMIT, 0, idum, &limit, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  return limit;
}

// groundflow_phreatic_only -yes (manual 6.581, F5 T5 U3a): the groundflow
// data is removed for the element groups that are not listed in a
// groundflow_phreatic_level_multiple_element_group record: -groundflow
// leaves their group_type and their group_groundflow_* records are deleted.
// Measured on the Professional 02-08-2026, which edits the database the
// same way (a group whose only type was -groundflow is an input error there
// too; the GNU names the reason).
void groundflow_phreatic_only_apply( void )

{
  long int only=-NO, ldum=0, igroup=0, max_group=0, imult=0, max_mult=0,
    i=0, n=0, nnew=0, listed=0, has_groundflow=0, idat=0,
    types[DATA_ITEM_SIZE], *groups=NULL;
  double ddum[1];

  if ( !db_active_index( GROUNDFLOW_PHREATIC_ONLY, 0, VERSION_NORMAL ) ) return;
  db( GROUNDFLOW_PHREATIC_ONLY, 0, &only, ddum, ldum, VERSION_NORMAL, GET );
  if ( only==-NO ) return;
  if ( only!=-YES ) db_error( GROUNDFLOW_PHREATIC_ONLY, 0 );

  db_max_index( GROUP_TYPE, max_group, VERSION_NORMAL, GET );
  db_max_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GROUP, max_mult,
    VERSION_NORMAL, GET );
  for ( igroup=0; igroup<=max_group; igroup++ ) {
    if ( !db_active_index( GROUP_TYPE, igroup, VERSION_NORMAL ) ) continue;
    listed = 0;
    for ( imult=0; imult<=max_mult && !listed; imult++ ) {
      if ( !db_active_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GROUP,
          imult, VERSION_NORMAL ) ) continue;
      n = db_len( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GROUP, imult,
        VERSION_NORMAL );
      groups = db_int( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GROUP, imult,
        VERSION_NORMAL );
      for ( i=0; i<n; i++ ) if ( groups[i]==igroup ) listed = 1;
    }
    if ( listed ) continue;
    db( GROUP_TYPE, igroup, types, ddum, n, VERSION_NORMAL, GET );
    has_groundflow = 0;
    for ( i=0, nnew=0; i<n; i++ ) {
      if ( types[i]==-GROUNDFLOW ) has_groundflow = 1;
      else types[nnew++] = types[i];
    }
    if ( !has_groundflow ) continue;
    if ( nnew==0 ) {
      pri( "Error: group ", igroup );
      pri( "has only -groundflow and is not in "
        "groundflow_phreatic_level_multiple_element_group "
        "(groundflow_phreatic_only -yes)." );
      exit(TN_EXIT_STATUS);
    }
    db_delete_index( GROUP_TYPE, igroup, VERSION_NORMAL );
    db( GROUP_TYPE, igroup, types, ddum, nnew, VERSION_NORMAL, PUT );
    for ( idat=0; idat<LAST_DUMMY; idat++ ) {
      if ( strncmp( db_name(idat), "group_groundflow_", 17 ) ) continue;
      if ( db_active_index( idat, igroup, VERSION_NORMAL ) )
        db_delete_index( idat, igroup, VERSION_NORMAL );
    }
  }
}
