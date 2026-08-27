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
// LOT 1 (infrastructure, commit cc0aba1): registration + dispatch +
// validation + print structure. LOT 2 (this lot): the NUMERICAL
// integration for 2D (quad4/quad9) - end-face selection from the
// reference point, stress integration over the end faces, quad9
// middle-plane averaging, outer/plot_switch, and the per-node averaged
// flag consumed by the -primary print. The 3D integration (hex8/hex27)
// remains pending (lot 3): 3D values stay 0 with a notice.

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
//   - thickness_switch: one switch per element group (3D concept:
//     shortest/longest element direction; consumed by L3)
//   - plot_switch: one switch per item (3 in 2D, 4 in 3D)
//   - average (default -yes, quad9/hex27) and outer (default -no) are
//     single switches consumed by the L2/L3 integration
//   - materi_stress must be a solved dof (the integration reads the
//     nodal stresses from NODE_DOF)
//   - 2D: the target groups may only contain quad4/quad9 elements
//     (the manual's -hex8/-hex27 are 3D; other 2D element types have
//     no cross-section faces)
void post_calcul_materi_stress_force_validate( void )

{
  long int idum[1], ldum=0, ngroups=0, nvalues=0, i=0, value=0,
    has_exclude=0, has_include=0, nforce_items=0, *ival=NULL,
    element=0, max_element=0, element_group=0, length_el=0, name=0;
  double ddum[1], direction[DATA_ITEM_SIZE];
  static long int not_implemented_notice=0;

  if ( ndim!=2 && ndim!=3 ) {
    pri( "Error: post_calcul -materi_stress -force requires a 2D or 3D "
         "calculation" );
    exit(TN_EXIT_STATUS);
  }

  // the integration reads the solved NODAL stresses
  if ( stres_indx<0 ) {
    pri( "Error: post_calcul -materi_stress -force requires materi_stress "
         "in the initia section (the stresses are read from node_dof)" );
    exit(TN_EXIT_STATUS);
  }

  // element_group (mandatory)
  if ( !db_active_index( POST_CALCUL_MATERI_STRESS_FORCE_ELEMENT_GROUP, 0,
       VERSION_NORMAL ) ) {
    pri( "Error: post_calcul -materi_stress -force requires "
         "post_calcul_materi_stress_force_element_group (the element "
         "groups for which the forces and moments are determined)" );
    exit(TN_EXIT_STATUS);
  }
  ival = get_new_int(DATA_ITEM_SIZE);
  db( POST_CALCUL_MATERI_STRESS_FORCE_ELEMENT_GROUP, 0, ival, ddum,
    ngroups, VERSION_NORMAL, GET );
  if ( ngroups<1 ) {
    pri( "Error: post_calcul_materi_stress_force_element_group needs at "
         "least one element group" );
    exit(TN_EXIT_STATUS);
  }

  // 2D: the target groups may only contain quad4/quad9 elements
  if ( ndim==2 ) {
    db_highest_index( ELEMENT, max_element, VERSION_NORMAL );
    for ( element=0; element<=max_element; element++ ) {
      if ( !db_active_index( ELEMENT, element, VERSION_NORMAL ) ) continue;
      db( ELEMENT_GROUP, element, &element_group, ddum, ldum,
        VERSION_NORMAL, GET );
      for ( i=0; i<ngroups; i++ ) if ( ival[i]==element_group ) break;
      if ( i>=ngroups ) continue;
      db( ELEMENT, element, ival, ddum, length_el, VERSION_NORMAL, GET );
      name = ival[0];
      if ( name!=-QUAD4 && name!=-QUAD9 ) {
        char str[256], str2[32];
        strcpy( str, "Error: post_calcul -materi_stress -force (2D) "
                     "supports only -quad4/-quad9 elements; element " );
        long_to_a( element, str2 );
        strcat( str, str2 );
        strcat( str, " of the target groups is another type" );
        pri( str );
        exit(TN_EXIT_STATUS);
      }
    }
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
  // manual 6.910/6.912); consumed by the L3 integration only

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

  if ( ndim==3 && !not_implemented_notice ) {
    pri( "post_calcul -materi_stress -force: the 3D integration is not "
         "yet implemented (lot 3) - the calculated values are set to 0" );
    not_implemented_notice = 1;
  }

}

// ---------------------------------------------------------------------
// LOT 2: 2D integration (quad4/quad9). Design decisions (documented in
// ProjectDocs/manual-developer/post_calcul_materi_stress_force.md):
//   1. END FACES: the thickness direction of the section is
//      t = (element centroid - reference_point) normalized IN-PLANE
//      (manual 6.914: the reference point defines outwards/inwards in
//      thickness direction; in 2D a point "at a large perpendicular
//      distance away" from a sheet pile). The two end faces = the 2
//      sides whose EXTERIOR normals are most PERPENDICULAR to t
//      (smallest |n*t|): the cross-section faces where the section
//      forces act. The other 2 sides (the structure surfaces, |n*t|~1)
//      produce no primary values. Ambiguous selection (distorted
//      elements, e.g. |n*t| equal for 3+ sides) warns and skips the
//      element.
//   2. STRESS SOURCE: the NODAL stresses (NODE_DOF, the solved
//      unknowns). Direct along the side (exact for the linear/quadratic
//      stress interpolation) and available in VERSION_NORMAL; the
//      element IP stresses (ELEMENT_DOF) would need an IP-to-side
//      extrapolation and were not used (documented decision).
//   3. QUADRATURE: 1D along the side with npol points (polynom.cc:
//      quad4 npol=2, quad9 npol=3): integration_gauss(2) for npol=2
//      (exact for the QUADRATIC moment integrand sigma_nn*(s-s_mid);
//      the element-area Lobatto(2) trapezoid would give 3/2 x the
//      moment) and integration_lobatto(3) (Simpson, exact for the
//      CUBIC moment integrand) for npol=3. Weights sum to 1 (Tochnog
//      convention: integral = side_length * sum(w*f)).
//   4. VALUES (per unit length l; plane 2D l=1, axisymmetric
//      l=2*PI*radial coordinate of the element centroid):
//      nor = int sigma_nn ds (SIGNED: positive = tension), she = the
//      absolute value of int sigma_nt ds (only the size, manual
//      6.913), mom = int sigma_nn*(s-s_mid) ds (signed; s along the
//      side = thickness direction of the section, s_mid = side
//      midpoint).
//   5. PLOT COMPONENTS: norx/nory = nor*t, nors = |nor|;
//      shex/shey = |she|*t, shes = |she|; momx/momy = mom*t,
//      moms = |mom|. The tension/compression sign of nor and mom is
//      carried by the VECTOR DIRECTION along t (outward/inward, per
//      the reference point); the s component is the physical SIZE
//      (manual 6.913: "the size of the vector formed by these
//      components indeed is the real physical size"). plot_switch -yes
//      inverts the x/y components of the item (drawing direction).
//   6. NODE ASSIGNMENT: the nodes of the two end faces are PRIMARY
//      (they receive their face's values). quad9 middle-plane nodes
//      (the 3 nodes NOT on the end faces) receive the AVERAGE of both
//      faces with average -yes (default), nothing (0) with -no. A node
//      shared by several elements receives the mean of its element
//      contributions (identical for a conforming mesh).
//   7. outer -yes: only the PRIMARY nodes at the maximum distance from
//      the reference point receive values (manual 6.915: "the nodes
//      which have the furthest distance relative to the reference
//      point"); the averaged nodes get nothing (documented decision).
//   8. AVERAGED FLAG: the per-node record
//      POST_CALCUL_MATERI_STRESS_FORCE_AVERAGED_NODE (-YES when the
//      node received middle-plane averaged values) is written in
//      VERSION_NORMAL and consumed by msf_node_is_averaged() of the
//      -primary print (print_materi_stress_force.cc).
// ---------------------------------------------------------------------

// local-node side tables (2D): side s has npol nodes (2 for quad4,
// 3 for quad9), ordered from the FIRST side node (iso = -1) to the
// LAST (iso = +1). Same numbering as the STATIC border_nodes_* tables
// of area.cc (quad4 corners 0,1,2,3; quad9 corners 0,2,8,6 with the
// mid-edge nodes 1,3,5,7 and the center 4).
static long int msf_border_nodes_quad4[] = {
  0, 1,
  1, 3,
  3, 2,
  2, 0 };

static long int msf_border_nodes_quad9[] = {
  0, 1, 2,
  2, 5, 8,
  8, 7, 6,
  6, 3, 0 };

// Integrates the section forces over ONE element side (2D): the side
// is the CROSS-SECTION face of the structure. The stress is
// interpolated from the NODAL values (node_sig: per node the 2D
// components sxx, syy, sxy) with the element polynomial (npol) along
// the side. nrm = exterior normal, tng = side tangent (unit vectors).
// Results per unit length l (see the design notes above):
//   nor = int sigma_nn ds      (signed, tension positive)
//   she = |int sigma_nt ds|    (always positive)
//   mom = int sigma_nn*dt ds   (signed), dt = the signed distance of
//         the integration point from the element centroid IN THE
//         THICKNESS DIRECTION t (manual 6.913: "a distance in
//         thickness direction dt relative to the middle of the
//         element"). dt = (x_q - centroid)*t makes the moment of BOTH
//         end faces of a section consistent (the side parametrization
//         s alone would give opposite signs for the two faces, since
//         one runs bottom->top and the other top->bottom).
static void msf_integrate_side_2d( long int npol, long int side_nodes[],
  double node_coord[], double node_sig[], double nrm[], double tng[],
  double centroid[], double thick[], double l, double &nor, double &she,
  double &mom )

{
  double iso[MPOINT], weight[MPOINT], h_pol[MPOINT], p_pol[MPOINT],
    xa[MDIM], xb[MDIM], xq[MDIM], side_len=0., sxx=0., syy=0., sxy=0.,
    snn=0., snt=0., dt=0., in_nn=0., in_nt=0., in_mom=0.;
  long int iq=0, i=0, idim=0;

  for ( i=0; i<ndim; i++ ) {
    xa[i] = node_coord[side_nodes[0]*MDIM + i];
    xb[i] = node_coord[side_nodes[npol-1]*MDIM + i];
  }
  array_subtract( xb, xa, xb, ndim );
  side_len = array_size( xb, ndim );
  if ( side_len<1.e-12 ) { // degenerate side: no section to integrate
    nor = she = mom = 0.;
    return;
  }

  if ( npol==2 ) integration_gauss( 2, iso, weight );
  else           integration_lobatto( 3, iso, weight );

  for ( iq=0; iq<npol; iq++ ) {
    interpolation_polynomial( iso[iq], npol, h_pol, p_pol );
    sxx = syy = sxy = 0.;
    for ( i=0; i<ndim; i++ ) xq[i] = 0.;
    for ( i=0; i<npol; i++ ) {
      sxx += h_pol[i]*node_sig[side_nodes[i]*3+0];
      syy += h_pol[i]*node_sig[side_nodes[i]*3+1];
      sxy += h_pol[i]*node_sig[side_nodes[i]*3+2];
      for ( idim=0; idim<ndim; idim++ )
        xq[idim] += h_pol[i]*node_coord[side_nodes[i]*MDIM+idim];
    }
    snn = sxx*nrm[0]*nrm[0] + syy*nrm[1]*nrm[1] + 2.*sxy*nrm[0]*nrm[1];
    snt = sxx*nrm[0]*tng[0] + syy*nrm[1]*tng[1]
        + sxy*(nrm[0]*tng[1] + nrm[1]*tng[0]);
    dt = 0.;
    for ( i=0; i<ndim; i++ ) dt += ( xq[i]-centroid[i] )*thick[i];
    in_nn  += weight[iq]*snn;
    in_nt  += weight[iq]*snt;
    in_mom += weight[iq]*snn*dt;
  }
  // sum(w)=1 (Tochnog convention) with ds = L/2*d(iso): the physical
  // integrals pick up the side length L.
  in_nn  *= side_len;
  in_nt  *= side_len;
  in_mom *= side_len;
  nor = in_nn/l;
  she = ( in_nt<0. ? -in_nt : in_nt )/l;
  mom = in_mom/l;
}

// Per-element contribution for ONE node of a quad4/quad9 element (2D).
// Determines the two end faces (design note 1), integrates them (notes
// 2-4), assigns the node (notes 5-7) and returns the 9 plot components
// in node_values[] (0 when the node receives nothing). got_value = 1
// when the node received a value, is_averaged = 1 when the value is
// the quad9 middle-plane average (the -primary flag, note 8).
static void msf_element_contribution_2d( long int element, long int name,
  long int element_group, long int inod, double reference_point[],
  long int average, long int outer, long int plot_switch[],
  double node_values[], long int &got_value, long int &is_averaged )

{
  long int ldum=0, length_el=0, i=0, j=0, inol=0, iside=0,
    nside=4, npol=0, nnol=0, node=0, axisym=-NO, iface=0,
    is_face_node=0, order[4], corners[4], itmp=0,
    icorner=0, ncorner=0, inod_pos=-1, *el=NULL;
  double ddum[1], *coord=NULL, *node_dof=NULL,
    coords[MDIM*MNOL], sigmas[3*MNOL], nrm[MDIM], tng[MDIM],
    centroid[MDIM], thick[MDIM], side_nrm[4][MDIM], side_tng[4][MDIM],
    side_score[4], face_nor[2], face_she[2], face_mom[2],
    max_dist=0., d=0., nor=0., she=0., mom=0., l=1., tol=0.;
  long int side_nodes[4][3], face_nodes[2][3];
  static long int warned_centroid=0, warned_ambiguous=0, warned_axisym=0;

  got_value = 0;
  is_averaged = 0;
  for ( i=0; i<9; i++ ) node_values[i] = 0.;

  npol = ( name==-QUAD4 ? 2 : 3 );
  if ( name==-QUAD4 ) {
    ncorner = 4;
    corners[0]=0; corners[1]=1; corners[2]=2; corners[3]=3;
  }
  else {
    ncorner = 4;
    corners[0]=0; corners[1]=2; corners[2]=8; corners[3]=6;
  }

  // element nodes: coordinates + nodal stresses (the solved unknowns)
  el = get_new_int(DATA_ITEM_SIZE);
  db( ELEMENT, element, el, ddum, length_el, VERSION_NORMAL, GET );
  nnol = length_el - 1;
  for ( inol=0; inol<nnol; inol++ ) {
    node = el[1+inol];
    coord = db_dbl( NODE, node, VERSION_NORMAL );
    for ( i=0; i<ndim; i++ ) coords[inol*MDIM+i] = coord[i];
    node_dof = db_dbl( NODE_DOF, node, VERSION_NORMAL );
    sigmas[inol*3+0] = node_dof[stres_indx + stress_indx(0,0)*nder];
    sigmas[inol*3+1] = node_dof[stres_indx + stress_indx(1,1)*nder];
    sigmas[inol*3+2] = node_dof[stres_indx + stress_indx(0,1)*nder];
    if ( node==inod ) inod_pos = inol;
  }
  delete[] el;
  if ( inod_pos<0 ) return; // node not in this element (should not happen)

  // element centroid (mean of the corners)
  array_set( centroid, 0., MDIM );
  for ( icorner=0; icorner<ncorner; icorner++ ) {
    for ( i=0; i<ndim; i++ )
      centroid[i] += coords[corners[icorner]*MDIM+i]/ncorner;
  }

  // thickness direction t = (centroid - reference_point), in-plane
  for ( i=0; i<ndim; i++ ) thick[i] = centroid[i] - reference_point[i];
  d = array_size( thick, ndim );
  if ( d<1.e-12 ) {
    // the reference point coincides with the element centroid: no
    // thickness direction (documented edge case; the msf_parse L1
    // tests use exactly this configuration to stay structure-only)
    if ( !warned_centroid ) {
      pri( "Warning: post_calcul -materi_stress -force: the reference "
           "point coincides with the centroid of an element - no forces "
           "calculated for it" );
      warned_centroid = 1;
    }
    return;
  }
  array_multiply( thick, thick, 1./d, ndim );

  // sides: tangent (first -> last side node), exterior normal
  // (perpendicular, away from the centroid; pattern interface.cc)
  // and the |n*t| score
  for ( iside=0; iside<nside; iside++ ) {
    for ( j=0; j<npol; j++ ) {
      side_nodes[iside][j] = ( name==-QUAD4
        ? msf_border_nodes_quad4[iside*npol+j]
        : msf_border_nodes_quad9[iside*npol+j] );
    }
    for ( i=0; i<ndim; i++ ) {
      tng[i] = coords[side_nodes[iside][npol-1]*MDIM+i]
             - coords[side_nodes[iside][0]*MDIM+i];
    }
    array_set( nrm, 0., MDIM );
    array_normalize( tng, ndim );
    nrm[0] = -tng[1];
    nrm[1] =  tng[0];
    d = 0.;
    for ( i=0; i<ndim; i++ ) {
      double mid = 0.5*( coords[side_nodes[iside][0]*MDIM+i]
                       + coords[side_nodes[iside][npol-1]*MDIM+i] );
      d += nrm[i]*( mid - centroid[i] );
    }
    if ( d<0. ) { nrm[0] = -nrm[0]; nrm[1] = -nrm[1]; }
    for ( i=0; i<MDIM; i++ ) {
      side_nrm[iside][i] = nrm[i];
      side_tng[iside][i] = tng[i];
    }
    side_score[iside] = scalar_dabs(
      nrm[0]*thick[0] + nrm[1]*thick[1] );
  }

  // end faces = the 2 sides with the SMALLEST |n*t| (the cross-section
  // faces, perpendicular to the thickness direction)
  for ( i=0; i<nside; i++ ) order[i] = i;
  for ( i=0; i<nside; i++ ) {
    for ( j=i+1; j<nside; j++ ) {
      if ( side_score[order[j]]<side_score[order[i]] ) {
        itmp = order[i]; order[i] = order[j]; order[j] = itmp;
      }
    }
  }
  // ambiguity: the 2nd end face competes with the 3rd side (distorted
  // element or a reference point on the element diagonal)
  tol = 1.e-6*scalar_dabs( side_score[order[1]] );
  if ( tol<1.e-12 ) tol = 1.e-12;
  if ( side_score[order[2]]-side_score[order[1]]<tol ) {
    if ( !warned_ambiguous ) {
      pri( "Warning: post_calcul -materi_stress -force: ambiguous "
           "end-face selection for an element (distorted geometry or a "
           "reference point on the element diagonal) - no forces "
           "calculated for it" );
      warned_ambiguous = 1;
    }
    return;
  }

  // per unit length l (manual 6.913): plane 2D l = 1; axisymmetric
  // l = 2*PI*radius (the radial coordinate x of the element centroid,
  // area.cc convention; "values cannot be calculated at the symmetry
  // axis with zero radius")
  db( GROUP_AXISYMMETRIC, element_group, &axisym, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  if ( axisym==-YES ) {
    l = 2.*PIRAD*centroid[0];
    if ( l<1.e-12 ) {
      if ( !warned_axisym ) {
        pri( "Warning: post_calcul -materi_stress -force: axisymmetric "
             "element at the symmetry axis (zero radius) - no forces "
             "calculated for it" );
        warned_axisym = 1;
      }
      return;
    }
  }

  // integrate the two end faces
  for ( j=0; j<2; j++ ) {
    iside = order[j];
    for ( i=0; i<npol; i++ ) face_nodes[j][i] = side_nodes[iside][i];
    msf_integrate_side_2d( npol, face_nodes[j], coords, sigmas,
      side_nrm[iside], side_tng[iside], centroid, thick, l,
      face_nor[j], face_she[j], face_mom[j] );
  }

  // role of the current node
  is_face_node = 0;
  for ( j=0; j<2 && !is_face_node; j++ ) {
    for ( i=0; i<npol; i++ ) {
      if ( side_nodes[order[j]][i]==inod_pos ) {
        is_face_node = 1;
        iface = j;
      }
    }
  }

  if ( is_face_node ) {
    // primary node of an end face; with outer -yes only the nodes at
    // the MAXIMUM distance from the reference point receive values
    if ( outer==-YES ) {
      max_dist = 0.;
      for ( inol=0; inol<nnol; inol++ ) {
        d = 0.;
        for ( i=0; i<ndim; i++ ) {
          d += ( coords[inol*MDIM+i] - reference_point[i] )
             * ( coords[inol*MDIM+i] - reference_point[i] );
        }
        if ( d>max_dist ) max_dist = d;
      }
      d = 0.;
      for ( i=0; i<ndim; i++ )
        d += ( coords[inod_pos*MDIM+i] - reference_point[i] )
           * ( coords[inod_pos*MDIM+i] - reference_point[i] );
      // distance tolerance: 1.e-8 relative to the element size
      tol = 1.e-8*array_size( coords, MDIM*nnol );
      if ( max_dist-d>tol ) return; // not an outer node
    }
    nor = face_nor[iface];
    she = face_she[iface];
    mom = face_mom[iface];
    got_value = 1;
  }
  else if ( npol==3 && average==-YES && outer!=-YES ) {
    // quad9 middle-plane node (NOT on any end face): the average of
    // the two end faces (manual 6.908; with outer -yes the averaged
    // nodes receive nothing - documented decision)
    nor = 0.5*( face_nor[0] + face_nor[1] );
    she = 0.5*( face_she[0] + face_she[1] );
    mom = 0.5*( face_mom[0] + face_mom[1] );
    got_value = 1;
    is_averaged = 1;
  }
  else {
    return; // surface node (quad4) or average -no: no value
  }

  // plot components (global x/y in the thickness direction t; the s
  // component is the physical SIZE; nor/mom keep their sign in the
  // vector direction, she is always positive - manual 6.913)
  node_values[0] = nor*thick[0];
  node_values[1] = nor*thick[1];
  node_values[2] = ( nor<0. ? -nor : nor );
  node_values[3] = she*thick[0];
  node_values[4] = she*thick[1];
  node_values[5] = she;
  node_values[6] = mom*thick[0];
  node_values[7] = mom*thick[1];
  node_values[8] = ( mom<0. ? -mom : mom );
  // plot_switch -yes: invert the drawing direction of the item vector
  for ( j=0; j<3; j++ ) {
    if ( plot_switch[j]==-YES ) {
      node_values[j*3+0] = -node_values[j*3+0];
      node_values[j*3+1] = -node_values[j*3+1];
    }
  }
}

// Per-node 2D calculation (called from post_calcul_materi_stress_force
// inside the parallel node loop of parallel_calcul_node, calcul.cc):
// scans the target element groups, collects the contributions of every
// element that contains the node (averaging identical contributions of
// a conforming mesh), writes the 9 items to result[] and the per-node
// averaged flag (POST_CALCUL_MATERI_STRESS_FORCE_AVERAGED_NODE).
static void msf_calculate_node_2d( long int inod, double result[] )

{
  long int idum[1], ldum=0, ngroups=0, nvalues=0, i=0, ig=0,
    element=0, max_element=0, length_el=0, element_group=0, name=0,
    average=-YES, outer=-NO, plot_switch[3], ncontrib=0, has_avg=0,
    has_prim=0, got=0, is_avg=0, value=-NO, length=0, *groups=NULL,
    *el=NULL;
  double ddum[1], reference[DATA_ITEM_SIZE], node_values[9], sum[9];

  array_set( sum, 0., 9 );
  ncontrib = has_avg = has_prim = 0;

  groups = get_new_int(DATA_ITEM_SIZE);
  db( POST_CALCUL_MATERI_STRESS_FORCE_ELEMENT_GROUP, 0, groups, ddum,
    ngroups, VERSION_NORMAL, GET );
  array_set( reference, 0., DATA_ITEM_SIZE );
  if ( db_active_index( POST_CALCUL_MATERI_STRESS_FORCE_REFERENCE_POINT,
       0, VERSION_NORMAL ) ) {
    db( POST_CALCUL_MATERI_STRESS_FORCE_REFERENCE_POINT, 0, idum,
      reference, nvalues, VERSION_NORMAL, GET );
  }
  // no reference_point: the (0,0) default of the validate() warning
  db( POST_CALCUL_MATERI_STRESS_FORCE_AVERAGE, 0, &average, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  db( POST_CALCUL_MATERI_STRESS_FORCE_OUTER, 0, &outer, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  plot_switch[0] = plot_switch[1] = plot_switch[2] = -NO;
  db( POST_CALCUL_MATERI_STRESS_FORCE_PLOT_SWITCH, 0, plot_switch, ddum,
    nvalues, VERSION_NORMAL, GET_IF_EXISTS );

  db_highest_index( ELEMENT, max_element, VERSION_NORMAL );
  el = get_new_int(DATA_ITEM_SIZE);
  for ( element=0; element<=max_element; element++ ) {
    if ( !db_active_index( ELEMENT, element, VERSION_NORMAL ) ) continue;
    db( ELEMENT_GROUP, element, &element_group, ddum, ldum,
      VERSION_NORMAL, GET );
    for ( ig=0; ig<ngroups; ig++ ) if ( groups[ig]==element_group ) break;
    if ( ig>=ngroups ) continue;
    db( ELEMENT, element, el, ddum, length_el, VERSION_NORMAL, GET );
    name = el[0];
    if ( name!=-QUAD4 && name!=-QUAD9 ) continue; // validated earlier
    for ( i=1; i<length_el; i++ ) if ( el[i]==inod ) break;
    if ( i>=length_el ) continue; // the node is not in this element
    msf_element_contribution_2d( element, name, element_group, inod,
      reference+ig*ndim, average, outer, plot_switch, node_values,
      got, is_avg );
    if ( got ) {
      for ( i=0; i<9; i++ ) sum[i] += node_values[i];
      ncontrib++;
      if ( is_avg ) has_avg = 1;
      else          has_prim = 1;
    }
  }
  delete[] el;
  delete[] groups;

  if ( ncontrib>0 ) {
    for ( i=0; i<9; i++ ) result[i] = sum[i]/ncontrib;
    // averaged flag: -YES when the node received (also) middle-plane
    // averages; -primary (the print) skips those nodes
    value = ( has_avg ? -YES : -NO );
  }
  else {
    for ( i=0; i<9; i++ ) result[i] = 0.;
    value = -NO;
  }
  length = 1;
  db( POST_CALCUL_MATERI_STRESS_FORCE_AVERAGED_NODE, inod, &value, ddum,
    length, VERSION_NORMAL, PUT );
}

// Per-node calculation of the -force family (called from
// calculate_operat in calcul.cc). The family is NODAL: inod<0 marks the
// POST_LINE_DOF / POST_POINT_DOF / POST_QUADRILATERAL_DOF branch of
// calculate(), which is rejected with a clear error (there is no
// element behind such records to integrate over).
// 2D (quad4/quad9): the LOT 2 integration (msf_calculate_node_2d).
// 3D (hex8/hex27): NOT implemented yet (lot 3) - the layout is
// produced (16 values, all 0) so that the dispatch never dies in
// db_error.
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
  if ( ndim==2 ) {
    msf_calculate_node_2d( inod, result );
  }
  else {
    // 3D: the numerical integration lands in lot 3
    for ( i=0; i<nitems; i++ ) result[i] = 0.;
  }
  length_result = nitems;
}
