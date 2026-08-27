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

// post_calcul -materi_stress -force (manual Professional 6.913): normal
// force, shear force and moment(s) of isoparametric elements (-quad4,
// -quad9, -hex8, -hex27) with a single element over the structure
// thickness (sheet piles, tunnel shells, ...). The result is a set of
// per-node items written to NODE_DOF_CALCUL, one slot per item:
//   2D (9 items):  norx_sig nory_sig nors_sig shex_sig shey_sig shes_sig
//                  momx_sig momy_sig moms_sig
//   3D (16 items): norx_sig nory_sig norz_sig nors_sig shex_sig shey_sig
//                  shez_sig shes_sig mom1x_sig mom1y_sig mom1z_sig
//                  mom1s_sig mom2x_sig mom2y_sig mom2z_sig mom2s_sig
// The x/y/z components are GLOBAL PLOT components (the vector is drawn
// in the structure thickness direction; only the vector SIZE
// (nors_sig, shes_sig, ...) is the physical value). The items are
// generated in calculate() (calcul.cc) next to the PHIMOB/TOTAL
// branches; the configuration records below are validated here.
//
// This file is the LOT 1 (infrastructure) of the materi_stress_force
// sub-sprint: registration + dispatch + validation + print structure.
// The NUMERICAL integration (integration of the stresses over the
// element sides, normals by cross product, reference-point orientation,
// quad9/hex27 averaging) is NOT implemented yet and lands in L2/L3.
// Until then every result value is 0.

// Number of result items of the -force family: 3 vector groups of
// (ndim components + size) in 2D (nor, she, mom) and 4 in 3D (nor, she,
// mom1, mom2). Both fit in MCALCUL=20 (tochnog.h): 9 <= 20 and 16 <= 20.
// Documented limitation: a 3D -force block leaves only 4 slots of the
// per-node NODE_DOF_CALCUL record for other post_calcul items; combining
// it with e.g. a 6-value -total stress block (6+16=22 > 20) aborts with
// the "MCALCUL too small" message of parallel_calcul_node (calcul.cc).
long int post_calcul_materi_stress_force_items( void )

{
  return ( ndim==2 ? 9 : 16 );
}

// Validation of the configuration records of the -force family (runs
// once per post_calcul -materi_stress -force record in calculate(),
// BEFORE the per-node loop). Errors exit with a clear message; 2D
// cases without reference_point warn and fall back to the documented
// default (reference point at the origin).
//
// Rules (manual Professional 6.908-6.917):
//   - element_group is mandatory (the target groups)
//   - direction_exclude XOR direction_include; in 3D one of them is
//     mandatory (which element sides produce forces/moments); in 2D
//     they are meaningless (warn and ignore)
//   - reference_point: one point per element group (ndim values each)
//   - thickness_switch: one switch per element group
//   - plot_switch: one switch per item (3 in 2D, 4 in 3D)
//   - average (default -yes, quad9/hex27) and outer (default -no) are
//     single switches consumed by the L2/L3 integration
void post_calcul_materi_stress_force_validate( void )

{
  long int idum[1], ldum=0, ngroups=0, nvalues=0, i=0, value=0,
    has_exclude=0, has_include=0, nforce_items=0, *ival=NULL;
  double ddum[1], direction[DATA_ITEM_SIZE];
  static long int not_implemented_notice=0;

  // element_group (mandatory)
  if ( !db_active_index( POST_CALCUL_MATERI_STRESS_FORCE_ELEMENT_GROUP, 0,
       VERSION_NORMAL ) ) {
    pri( "Error: post_calcul -materi_stress -force requires "
         "post_calcul_materi_stress_force_element_group (the element "
         "groups for which the forces and moments are determined)" );
    exit(TN_EXIT_STATUS);
  }
  db( POST_CALCUL_MATERI_STRESS_FORCE_ELEMENT_GROUP, 0, idum, ddum,
    ngroups, VERSION_NORMAL, GET );
  if ( ngroups<1 ) {
    pri( "Error: post_calcul_materi_stress_force_element_group needs at "
         "least one element group" );
    exit(TN_EXIT_STATUS);
  }

  // direction: exclude XOR include (manual 6.913: "Only one of ... and
  // ... should be specified, not both"); 3D requires one of them
  has_exclude = db_active_index(
    POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_EXCLUDE, 0, VERSION_NORMAL );
  has_include = db_active_index(
    POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_INCLUDE, 0, VERSION_NORMAL );
  if ( has_exclude && has_include ) {
    pri( "Error: post_calcul_materi_stress_force_direction_exclude and "
         "post_calcul_materi_stress_force_direction_include are mutually "
         "exclusive (manual Professional 6.913)" );
    exit(TN_EXIT_STATUS);
  }
  if ( ndim==3 && !has_exclude && !has_include ) {
    pri( "Error: in 3D, post_calcul -materi_stress -force requires either "
         "post_calcul_materi_stress_force_direction_exclude or "
         "post_calcul_materi_stress_force_direction_include (to select the "
         "element sides, manual Professional 6.909/6.911)" );
    exit(TN_EXIT_STATUS);
  }
  if ( ndim==2 && ( has_exclude || has_include ) ) {
    pri( "Warning: post_calcul_materi_stress_force_direction_* is a 3D "
         "concept, ignored in this 2D calculation" );
  }
  if ( has_exclude ) {
    db( POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_EXCLUDE, 0, idum,
      direction, nvalues, VERSION_NORMAL, GET );
    if ( nvalues!=ndim )
      db_error( POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_EXCLUDE, 0 );
  }
  if ( has_include ) {
    db( POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_INCLUDE, 0, idum,
      direction, nvalues, VERSION_NORMAL, GET );
    if ( nvalues!=ndim )
      db_error( POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_INCLUDE, 0 );
  }
  // the direction_*_epsilon records are single doubles (default 1.e-8,
  // manual 6.910/6.912); consumed by the L2/L3 integration only

  // reference_point: one point (ndim values) per element group
  if ( db_active_index( POST_CALCUL_MATERI_STRESS_FORCE_REFERENCE_POINT,
       0, VERSION_NORMAL ) ) {
    db( POST_CALCUL_MATERI_STRESS_FORCE_REFERENCE_POINT, 0, idum,
      direction, nvalues, VERSION_NORMAL, GET );
    if ( nvalues!=ngroups*ndim ) {
      char str[256], str2[32];
      strcpy( str, "Error: post_calcul_materi_stress_force_reference_point needs one point (" );
      long_to_a( ndim, str2 );
      strcat( str, str2 );
      strcat( str, " values) per element group: " );
      long_to_a( ngroups*ndim, str2 );
      strcat( str, str2 );
      strcat( str, " values expected, " );
      long_to_a( nvalues, str2 );
      strcat( str, str2 );
      strcat( str, " given" );
      pri( str );
      exit(TN_EXIT_STATUS);
    }
  }
  else {
    if ( ndim==3 ) {
      pri( "Error: in 3D, post_calcul -materi_stress -force requires "
           "post_calcul_materi_stress_force_reference_point (one point per "
           "element group, manual Professional 6.914)" );
      exit(TN_EXIT_STATUS);
    }
    else {
      pri( "Warning: post_calcul -materi_stress -force without "
           "post_calcul_materi_stress_force_reference_point in 2D: the "
           "default reference point (0,0) is used (documented decision)" );
    }
  }

  // thickness_switch: one switch per element group (manual 6.917)
  if ( db_active_index( POST_CALCUL_MATERI_STRESS_FORCE_THICKNESS_SWITCH,
       0, VERSION_NORMAL ) ) {
    ival = get_new_int(DATA_ITEM_SIZE);
    db( POST_CALCUL_MATERI_STRESS_FORCE_THICKNESS_SWITCH, 0, ival, ddum,
      nvalues, VERSION_NORMAL, GET );
    if ( nvalues!=ngroups ) {
      char str[256], str2[32];
      strcpy( str, "Error: post_calcul_materi_stress_force_thickness_switch needs one switch per element group (" );
      long_to_a( ngroups, str2 );
      strcat( str, str2 );
      strcat( str, " expected, " );
      long_to_a( nvalues, str2 );
      strcat( str, str2 );
      strcat( str, " given)" );
      pri( str );
      delete[] ival;
      exit(TN_EXIT_STATUS);
    }
    for ( i=0; i<nvalues; i++ ) {
      value = ival[i];
      if ( value!=-YES && value!=-NO ) {
        delete[] ival;
        db_error( POST_CALCUL_MATERI_STRESS_FORCE_THICKNESS_SWITCH, 0 );
      }
    }
    delete[] ival;
  }

  // plot_switch: one switch per VECTOR item - 3 in 2D (normal force,
  // shear force, moment), 4 in 3D (normal force, shear force, two
  // moments) - manual 6.916 ("In 2D you need to specify a switch for
  // the normal force, shear force and moment. In 3D ... for the normal
  // force, shear force and two moments")
  if ( db_active_index( POST_CALCUL_MATERI_STRESS_FORCE_PLOT_SWITCH, 0,
       VERSION_NORMAL ) ) {
    ival = get_new_int(DATA_ITEM_SIZE);
    db( POST_CALCUL_MATERI_STRESS_FORCE_PLOT_SWITCH, 0, ival, ddum,
      nvalues, VERSION_NORMAL, GET );
    nforce_items = ( ndim==2 ? 3 : 4 );
    if ( nvalues!=nforce_items ) {
      char str[256], str2[32];
      strcpy( str, "Error: post_calcul_materi_stress_force_plot_switch needs one switch per item (" );
      long_to_a( nforce_items, str2 );
      strcat( str, str2 );
      strcat( str, " in this " );
      long_to_a( ndim, str2 );
      strcat( str, str2 );
      strcat( str, "D calculation, " );
      long_to_a( nvalues, str2 );
      strcat( str, str2 );
      strcat( str, " given)" );
      pri( str );
      delete[] ival;
      exit(TN_EXIT_STATUS);
    }
    for ( i=0; i<nvalues; i++ ) {
      value = ival[i];
      if ( value!=-YES && value!=-NO ) {
        delete[] ival;
        db_error( POST_CALCUL_MATERI_STRESS_FORCE_PLOT_SWITCH, 0 );
      }
    }
    delete[] ival;
  }

  // average (default -yes, quad9/hex27 only) and outer (default -no):
  // single switches consumed by the L2/L3 integration
  if ( db_active_index( POST_CALCUL_MATERI_STRESS_FORCE_AVERAGE, 0,
       VERSION_NORMAL ) ) {
    db( POST_CALCUL_MATERI_STRESS_FORCE_AVERAGE, 0, &value, ddum, ldum,
      VERSION_NORMAL, GET );
    if ( value!=-YES && value!=-NO )
      db_error( POST_CALCUL_MATERI_STRESS_FORCE_AVERAGE, 0 );
  }
  if ( db_active_index( POST_CALCUL_MATERI_STRESS_FORCE_OUTER, 0,
       VERSION_NORMAL ) ) {
    db( POST_CALCUL_MATERI_STRESS_FORCE_OUTER, 0, &value, ddum, ldum,
      VERSION_NORMAL, GET );
    if ( value!=-YES && value!=-NO )
      db_error( POST_CALCUL_MATERI_STRESS_FORCE_OUTER, 0 );
  }

  if ( !not_implemented_notice ) {
    pri( "post_calcul -materi_stress -force: the numerical integration is "
         "not yet implemented (lot 1 = registration + dispatch + print "
         "structure) - the calculated values are set to 0" );
    not_implemented_notice = 1;
  }

}

// Per-node calculation of the -force family (called from
// calculate_operat in calcul.cc). LOT 1: the NODAL layout is produced
// (length_result = 9 in 2D / 16 in 3D, all values 0) so that the
// dispatch never dies in db_error; the numerical values arrive with the
// L2/L3 integration. The family is NODAL: inod<0 marks the
// POST_LINE_DOF / POST_POINT_DOF / POST_QUADRILATERAL_DOF branch of
// calculate(), which is rejected with a clear error (there is no
// element behind such records to integrate over).
void post_calcul_materi_stress_force( double unknown_values[],
  long int inod, double coord[], double dof[], double result[],
  long int &length_result )

{
  long int nitems=0, i=0;

  if ( inod<0 ) {
    pri( "Error: post_calcul -materi_stress -force is a NODAL calculation; "
         "POST_LINE_DOF/POST_POINT_DOF/POST_QUADRILATERAL_DOF records are "
         "not supported (there is no element behind them)" );
    exit(TN_EXIT_STATUS);
  }

  nitems = post_calcul_materi_stress_force_items();
  for ( i=0; i<nitems; i++ ) result[i] = 0.;
  length_result = nitems;
}
