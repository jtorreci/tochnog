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
// validation + print structure. LOT 2 (commit 1babf2d): the NUMERICAL
// integration for 2D (quad4/quad9) - end-face selection from the
// reference point, stress integration over the end faces, quad9
// middle-plane averaging, outer/plot_switch, and the per-node averaged
// flag consumed by the -primary print. LOT 3 (this lot): the 3D
// integration (hex8/hex27) - face selection from direction_exclude/
// direction_include, per-face thickness/length directions, 2D face
// quadrature, mom1/mom2, hex27 middle-plane averaging. The 3D design
// decisions are documented in
// ProjectDocs/manual-developer/post_calcul_materi_stress_force.md.

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
  //   - the target groups may only contain quad4/quad9 elements (2D) or
  //     hex8/hex27 elements (3D): the manual's isoparametric elements
  //     with a cross-section; other element types have no section faces
  //     (validated below)
void post_calcul_materi_stress_force_validate( void )

{
  long int idum[1], ldum=0, ngroups=0, nvalues=0, i=0, value=0,
    has_exclude=0, has_include=0, nforce_items=0, *ival=NULL,
    element=0, max_element=0, element_group=0, length_el=0, name=0;
  double ddum[1], direction[DATA_ITEM_SIZE];
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

  // the target groups may only contain the isoparametric section
  // elements of the manual 6.913: quad4/quad9 in 2D, hex8/hex27 in 3D
  {
    db_highest_index( ELEMENT, max_element, VERSION_NORMAL );
    for ( element=0; element<=max_element; element++ ) {
      if ( !db_active_index( ELEMENT, element, VERSION_NORMAL ) ) continue;
      db( ELEMENT_GROUP, element, &element_group, ddum, ldum,
        VERSION_NORMAL, GET );
      for ( i=0; i<ngroups; i++ ) if ( ival[i]==element_group ) break;
      if ( i>=ngroups ) continue;
      db( ELEMENT, element, ival, ddum, length_el, VERSION_NORMAL, GET );
      name = ival[0];
      if ( ndim==2 && name!=-QUAD4 && name!=-QUAD9 ) {
        char str[256], str2[32];
        strcpy( str, "Error: post_calcul -materi_stress -force (2D) "
                     "supports only -quad4/-quad9 elements; element " );
        long_to_a( element, str2 );
        strcat( str, str2 );
        strcat( str, " of the target groups is another type" );
        pri( str );
        exit(TN_EXIT_STATUS);
      }
      if ( ndim==3 && name!=-HEX8 && name!=-HEX27 ) {
        char str[256], str2[32];
        strcpy( str, "Error: post_calcul -materi_stress -force (3D) "
                     "supports only -hex8/-hex27 elements; element " );
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

// ---------------------------------------------------------------------
// LOT 3: 3D integration (hex8/hex27). Design decisions (documented in
// ProjectDocs/manual-developer/post_calcul_materi_stress_force.md; the
// analytical tests msf_sheet3d / msf_tunnel3d are the arbiter):
//   1. FACE SELECTION: for each of the 6 faces the physical normal n is
//      built from the corner cross product (interface.cc pattern),
//      oriented away from the element centroid. With
//      direction_exclude, faces with |n*dir| > 1-eps are neglected
//      (the structure end caps, e.g. the tunnel end faces); with
//      direction_include, faces with |n*dir| < eps are neglected (the
//      structure surface faces). The manual 6.913 sanity check "exactly
//      4 sides should be consistent with the specified direction" =
//      exactly 4 faces with |n*dir| < eps (the section element has 4
//      faces perpendicular to the direction and 2 parallel); otherwise
//      the element is skipped with a warning. The CANDIDATES (the faces
//      not neglected) are: exclude -> the 4 perpendicular faces (outer/
//      inner/circumferential of the tunnel), include -> the 2 parallel
//      faces (the +-height sections of the sheet pile).
//   2. END FACES: the 2 opposing candidate faces where the forces and
//      moments are primarily calculated (manual 6.908). Selection rule
//      (2D decision 1 generalized): the 2 candidates with the SMALLEST
//      |n*t_hat|, t_hat = normalize(element centroid - reference
//      point) (the thickness direction of the structure). Tunnel:
//      among outer/inner/circumferential the pair perpendicular to the
//      radial direction = the 2 circumferential faces (where the hoop
//      force acts); sheet pile: the only 2 candidates. The other
//      candidates (the tunnel outer/inner surfaces) are structure
//      surfaces and produce no primary values (2D pattern). Guards:
//      the 2 end faces must point to OPPOSITE sides of the element
//      (negative normal dot product; for curved ring elements the
//      +-circumferential faces are not parallel - tolerance 0.1) and
//      the selection must be clearly separated from the 3rd-best
//      candidate (distorted element -> warning + skip).
//   3. PER-FACE FRAME: n (outward normal), t = the in-face direction
//      with the SHORTEST physical extent (manual 6.917; thickness
//      switch -yes -> the LONGEST), oriented away from the reference
//      point (plot direction, 2D decision), l = n x t (unit, the
//      length direction). Per unit length l = the element size in the
//      length direction = max(|edge*t| projected on l over the face
//      corners) (manual 6.913: "the length of an element is determined
//      from the nodal coordinates differences in length direction").
//   4. FACE QUADRATURE: 2D tensor product on the face, consistent with
//      the 1D 2D decision (lot 2 note 3): integration_gauss(2) x
//      integration_gauss(2) for hex8 (4 points, exact for the
//      quadratic moment integrand on parallelogram faces) and
//      integration_lobatto(3) x integration_lobatto(3) for hex27 (9
//      points, Simpson, exact for the cubic moment integrand on
//      regular faces; the manual 6.913 regularity condition "elements
//      should be regular shaped in length direction" keeps the face
//      Jacobian constant so the rule stays exact). The stress is
//      interpolated from the NODAL values (2D decision 2) with the
//      face shape functions.
//   5. VALUES (per unit length l): nor = int int sigma_nn dA / l
//      (signed, tension positive), she = |int int sigma_nt dA| / l
//      (always positive), mom1 = int int sigma_nn*dt dA / l with
//      dt = (x - face_mid)*t (distance in THICKNESS direction from the
//      MIDDLE OF THE FACE, manual 6.913 "a distance in thickness
//      direction dt relative to the middle of the element"; for the
//      straight elements the face middle coincides with the element
//      centroid projection - the 2D lot 2 decision - while for the
//      CURVED ring elements the corner-mean centroid lies on the chord
//      and would add a spurious moment, verified in msf_tunnel3d),
//      mom2 = int int sigma_nn*dl dA / l with dl = (x - face_mid)*l
//      (distance in LENGTH direction; manual 6.913: "moment
//      contributions of normal stresses sigmann with a distance in
//      length direction dl relative to the middle of the element,
//      integrated over thickness direction").
//   6. NODE ASSIGNMENT: the nodes of the two end faces are PRIMARY
//      (their face's values). hex27 middle-plane nodes (the 9 nodes
//      NOT on the end faces, i.e. the face between the two end faces,
//      manual 6.908 "the nodes in the plane between the two end
//      faces") receive the AVERAGE of both end faces with average -yes
//      (default). outer -yes: only the PRIMARY nodes at the maximum
//      distance from the reference point (2D decision 7).
//   7. PLOT COMPONENTS (16 items): norx/nory/norz = nor*t (the vector
//      drawn in the thickness direction of the structure, manual
//      6.913; t = the face's t for primary nodes, the element t_hat
//      for averaged nodes), nors = |nor|; shex/y/z = she*t, shes =
//      she; mom1x/y/z = mom1*t, mom1s = |mom1|; mom2x/y/z = mom2*t,
//      mom2s = |mom2|. plot_switch inverts the x/y/z components of
//      each of the 4 items (manual 6.916).
//   8. KNOWN LIMITATION (documented, fix in a separate lot): the hex8
//      element with 1 element over the thickness suffers shear locking
//      in bending exactly like the quad4 in 2D (measured in
//      msf_sheet3d_hex8: mom1 ~ 0.23 x the static value). The
//      validation tests use hex27 (msf_sheet3d, msf_tunnel3d,
//      msf_hex27_avg).
// ---------------------------------------------------------------------

// local-node face tables (3D): face f has npol*npol nodes (4 for hex8,
// 9 for hex27), ordered as the in-face tensor product position
// p = iu*npol + iv (iv fastest), i.e. the same layout as the quad4/
// quad9 side tables of the 2D integration. Same numbering as the
// STATIC border_nodes_hex8/border_nodes_hex27 tables of area.cc
// (S Y N C - keep in sync with area.cc:157-186; hex8 corners 0..7,
// hex27 corners 0..7, mid-edges 8..19, mid-faces 20..25, center 26).
static long int msf_border_nodes_hex8[] = {
  0, 1, 2, 3,
  4, 5, 6, 7,
  0, 1, 4, 5,
  1, 3, 5, 7,
  2, 3, 6, 7,
  0, 2, 4, 6 };

static long int msf_border_nodes_hex27[] = {
   0,  1,  2,  3,  4,  5,  6,  7,  8,
  18, 19, 20, 21, 22, 23, 24, 25, 26,
   0,  1,  2,  9, 10, 11, 18, 19, 20,
   2,  5,  8, 11, 14, 17, 20, 23, 26,
   6,  7,  8, 15, 16, 17, 24, 25, 26,
   0,  3,  6,  9, 12, 15, 18, 21, 24 };

// the 8 CORNER positions of the hex27 in the tensor ordering
// (i_xi, i_eta, i_zeta) in {0,2}^3 -> iz*9 + ieta*3 + ixi: the same
// order as the hex8 corners (used for the element centroid).
static long int msf_corner_nodes_hex27[] = {
   0,  2,  6,  8, 18, 20, 24, 26 };

// 3D stress components: the compacted symmetric layout read from
// node_dof[stres_indx + stress_indx(i,j)*nder] (miscel.cc stress_indx:
// sxx=0, sxy=1, sxz=2, syy=3, syz=4, szz=5).
static const long int msf_sig_sxx=0, msf_sig_sxy=1, msf_sig_sxz=2,
  msf_sig_syy=3, msf_sig_syz=4, msf_sig_szz=5;

// Integrates the section quantities over ONE face of a hex8/hex27
// element (3D): the face is one of the two END faces of the section.
// The stress is interpolated from the NODAL values (node_sig: per node
// the 6 components sxx sxy sxz syy syz szz) with the face polynomial
// (npol per in-face direction). nrm = outward normal, thick = the
// face thickness direction t (unit), leng = the face length direction
// l (unit), face_mid = the middle of the face (the mean of its nodes;
// "the middle of the element" of the manual 6.913 in the thickness/
// length directions), l = the element size in the length direction.
// Results per unit length l:
//   nor  = int int sigma_nn dA / l   (signed, tension positive)
//   she  = |int int sigma_nt dA| / l (always positive)
//   mom1 = int int sigma_nn*dt dA / l, dt = (x-face_mid)*t
//   mom2 = int int sigma_nn*dl dA / l, dl = (x-face_mid)*l
static void msf_integrate_face_3d( long int npol, long int face_nodes[],
  double node_coord[], double node_sig[], double nrm[], double thick[],
  double leng[], double face_mid[], double l, double &nor, double &she,
  double &mom1, double &mom2 )

{
  double iso[MPOINT], weight[MPOINT], h_u[MPOINT], p_u[MPOINT],
    h_v[MPOINT], p_v[MPOINT], xq[MDIM], xu[MDIM], xv[MDIM],
    sxx=0., sxy=0., sxz=0., syy=0., syz=0., szz=0., snn=0., snt=0.,
    dt=0., dl=0., in_nn=0., in_nt=0., in_m1=0., in_m2=0., jac=0.,
    nx=0., ny=0., nz=0.;
  long int iu=0, iv=0, iu_node=0, iv_node=0, i=0, idim=0, nface=0,
    inol=0;
  double hu=0., hv=0., dhu=0., dhv=0., x=0.;

  nface = npol*npol;
  if ( npol==2 ) integration_gauss( 2, iso, weight );
  else           integration_lobatto( 3, iso, weight );

  for ( iu=0; iu<npol; iu++ ) {
    interpolation_polynomial( iso[iu], npol, h_u, p_u );
    for ( iv=0; iv<npol; iv++ ) {
      interpolation_polynomial( iso[iv], npol, h_v, p_v );
      sxx = sxy = sxz = syy = syz = szz = 0.;
      array_set( xq, 0., MDIM );
      array_set( xu, 0., MDIM );
      array_set( xv, 0., MDIM );
      for ( i=0; i<nface; i++ ) {
        iu_node = i/npol;
        iv_node = i%npol;
        hu = h_u[iu_node];
        hv = h_v[iv_node];
        dhu = p_u[iu_node];
        dhv = p_v[iv_node];
        inol = face_nodes[i];
        sxx += hu*hv*node_sig[inol*6+msf_sig_sxx];
        sxy += hu*hv*node_sig[inol*6+msf_sig_sxy];
        sxz += hu*hv*node_sig[inol*6+msf_sig_sxz];
        syy += hu*hv*node_sig[inol*6+msf_sig_syy];
        syz += hu*hv*node_sig[inol*6+msf_sig_syz];
        szz += hu*hv*node_sig[inol*6+msf_sig_szz];
        for ( idim=0; idim<ndim; idim++ ) {
          x = node_coord[inol*MDIM+idim];
          xq[idim] += hu*hv*x;
          xu[idim] += dhu*hv*x;
          xv[idim] += hu*dhv*x;
        }
      }
      // face Jacobian |du x dv|
      nx = xu[1]*xv[2] - xu[2]*xv[1];
      ny = xu[2]*xv[0] - xu[0]*xv[2];
      nz = xu[0]*xv[1] - xu[1]*xv[0];
      jac = sqrt( nx*nx + ny*ny + nz*nz );
      if ( jac<1.e-20 ) continue; // degenerate face region
      snn = nrm[0]*nrm[0]*sxx + nrm[1]*nrm[1]*syy + nrm[2]*nrm[2]*szz
          + 2.*( nrm[0]*nrm[1]*sxy + nrm[0]*nrm[2]*sxz
               + nrm[1]*nrm[2]*syz );
      snt = nrm[0]*thick[0]*sxx + nrm[1]*thick[1]*syy
          + nrm[2]*thick[2]*szz
          + ( nrm[0]*thick[1] + nrm[1]*thick[0] )*sxy
          + ( nrm[0]*thick[2] + nrm[2]*thick[0] )*sxz
          + ( nrm[1]*thick[2] + nrm[2]*thick[1] )*syz;
      dt = dl = 0.;
      for ( idim=0; idim<ndim; idim++ ) {
        dt += ( xq[idim]-face_mid[idim] )*thick[idim];
        dl += ( xq[idim]-face_mid[idim] )*leng[idim];
      }
      in_nn  += 4.*weight[iu]*weight[iv]*snn*jac;
      in_nt  += 4.*weight[iu]*weight[iv]*snt*jac;
      in_m1  += 4.*weight[iu]*weight[iv]*snn*dt*jac;
      in_m2  += 4.*weight[iu]*weight[iv]*snn*dl*jac;
    }
  }
  // the weights sum to 1 per direction (Tochnog convention), the
  // integration intervals are [-1,1] each: the factor 4 = the
  // d(xi)*d(eta) = 2*2 of the two in-face coordinates (the 2D analog of
  // the implicit 2*side_len/2 = side_len of the 1D side integration of
  // the lot 2 code). Verified: the tunnel ring nor = p*R and the sheet
  // bending mom1 = E*kappa/12 come out EXACT with the factor; without
  // it they are 4x too small.
  nor = in_nn/l;
  she = ( in_nt<0. ? -in_nt : in_nt )/l;
  mom1 = in_m1/l;
  mom2 = in_m2/l;
}

// Per-element contribution for ONE node of a hex8/hex27 element (3D).
// Face selection (design notes 1-2), per-face frame (note 3),
// integration (notes 4-5), node assignment (notes 6-7). Returns the 16
// plot components in node_values[] (0 when the node receives nothing).
// got_value = 1 when the node received a value, is_averaged = 1 when
// the value is the hex27 middle-plane average (the -primary flag).
static void msf_element_contribution_3d( long int element, long int name,
  long int element_group, long int inod, double reference_point[],
  long int thickness_switch, long int average, long int outer,
  long int plot_switch[], double node_values[], long int &got_value,
  long int &is_averaged )

{
  long int ldum=0, length_el=0, i=0, j=0, inol=0, iside=0, nside=6,
    npol=0, nnol=0, node=0, iface=0, is_face_node=0, order[6], itmp=0,
    ncorner=0, icorner=0, inod_pos=-1, *el=NULL, ncand=0, nper=0,
    iend=0, jend=0, cand[6], dir_has=0, exclude=0, nface=0;
  double ddum[1], *coord=NULL, *node_dof=NULL,
    coords[MDIM*MNOL], sigmas[6*MNOL], nrm[MDIM], tng[MDIM],
    centroid[MDIM], t_global[MDIM], side_nrm[6][MDIM],
    face_t[2][MDIM], face_l[2][MDIM], face_len[2], face_nor[2],
    face_she[2], face_m1[2], face_m2[2], max_dist=0., d=0., nor=0.,
    she=0., mom1=0., mom2=0., tol=0., dir[MDIM], eps=1.e-8,
    e1[MDIM], e2[MDIM], mid[MDIM], s1=0., s2=0., len1=0., len2=0.,
    cross=0.;
  long int side_nodes[6][9], face_nodes[2][9], corners[6][4];
  static long int warned_4sides=0, warned_centroid=0,
    warned_ambiguous=0, warned_opposing=0, warned_length=0,
    warned_orient=0;

  got_value = 0;
  is_averaged = 0;
  for ( i=0; i<16; i++ ) node_values[i] = 0.;

  npol = ( name==-HEX8 ? 2 : 3 );
  nface = npol*npol;
  ncorner = 8;

  // element nodes: coordinates + nodal stresses (the solved unknowns)
  el = get_new_int(DATA_ITEM_SIZE);
  db( ELEMENT, element, el, ddum, length_el, VERSION_NORMAL, GET );
  nnol = length_el - 1;
  for ( inol=0; inol<nnol; inol++ ) {
    node = el[1+inol];
    coord = db_dbl( NODE, node, VERSION_NORMAL );
    for ( i=0; i<ndim; i++ ) coords[inol*MDIM+i] = coord[i];
    node_dof = db_dbl( NODE_DOF, node, VERSION_NORMAL );
    sigmas[inol*6+0] = node_dof[stres_indx + stress_indx(0,0)*nder];
    sigmas[inol*6+1] = node_dof[stres_indx + stress_indx(0,1)*nder];
    sigmas[inol*6+2] = node_dof[stres_indx + stress_indx(0,2)*nder];
    sigmas[inol*6+3] = node_dof[stres_indx + stress_indx(1,1)*nder];
    sigmas[inol*6+4] = node_dof[stres_indx + stress_indx(1,2)*nder];
    sigmas[inol*6+5] = node_dof[stres_indx + stress_indx(2,2)*nder];
    if ( node==inod ) inod_pos = inol;
  }
  delete[] el;
  if ( inod_pos<0 ) return; // node not in this element (should not happen)

  // element centroid (mean of the 8 CORNERS: local nodes 0-7 for the
  // hex8, the corner positions 0,2,6,8,18,20,24,26 for the hex27 -
  // local nodes 0-7 of the hex27 are the z=-1 face, NOT the corners)
  array_set( centroid, 0., MDIM );
  for ( icorner=0; icorner<ncorner; icorner++ ) {
    long int cn = ( name==-HEX8 ? icorner
                                 : msf_corner_nodes_hex27[icorner] );
    for ( i=0; i<ndim; i++ )
      centroid[i] += coords[cn*MDIM+i]/ncorner;
  }

  // thickness direction of the structure t_hat = (centroid -
  // reference_point) normalized (2D decision 1; also the plot
  // direction of the averaged nodes)
  for ( i=0; i<ndim; i++ ) t_global[i] = centroid[i] - reference_point[i];
  d = array_size( t_global, ndim );
  if ( d<1.e-12 ) {
    // the reference point coincides with the element centroid: no
    // thickness direction (documented edge case)
    if ( !warned_centroid ) {
      pri( "Warning: post_calcul -materi_stress -force: the reference "
           "point coincides with the centroid of an element - no forces "
           "calculated for it" );
      warned_centroid = 1;
    }
    return;
  }
  array_multiply( t_global, t_global, 1./d, ndim );

  // the direction record (single per run; exclude XOR include
  // validated in post_calcul_materi_stress_force_validate)
  dir_has = db_active_index(
    POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_EXCLUDE, 0,
    VERSION_NORMAL );
  if ( dir_has ) {
    exclude = 1;
    db( POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_EXCLUDE, 0, &i, dir,
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_EXCLUDE_EPSILON, 0,
      &i, &eps, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  }
  else {
    db( POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_INCLUDE, 0, &i, dir,
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_INCLUDE_EPSILON, 0,
      &i, &eps, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  }
  d = array_size( dir, ndim );
  if ( d>1.e-12 ) array_multiply( dir, dir, 1./d, ndim );

  // faces: border table, outward normal (corner cross product,
  // oriented away from the centroid; interface.cc pattern), in-face
  // edge vectors e1 (corner0->corner1) and e2 (corner0->corner3)
  for ( iside=0; iside<nside; iside++ ) {
    for ( j=0; j<nface; j++ ) {
      side_nodes[iside][j] = ( name==-HEX8
        ? msf_border_nodes_hex8[iside*nface+j]
        : msf_border_nodes_hex27[iside*nface+j] );
    }
    corners[iside][0] = side_nodes[iside][0];
    corners[iside][1] = side_nodes[iside][npol-1];
    corners[iside][2] = side_nodes[iside][npol*(npol-1)];
    corners[iside][3] = side_nodes[iside][nface-1];
    for ( i=0; i<ndim; i++ ) {
      e1[i] = coords[corners[iside][1]*MDIM+i]
            - coords[corners[iside][0]*MDIM+i];
      e2[i] = coords[corners[iside][3]*MDIM+i]
            - coords[corners[iside][0]*MDIM+i];
    }
    nrm[0] = e1[1]*e2[2] - e1[2]*e2[1];
    nrm[1] = e1[2]*e2[0] - e1[0]*e2[2];
    nrm[2] = e1[0]*e2[1] - e1[1]*e2[0];
    cross = sqrt( nrm[0]*nrm[0] + nrm[1]*nrm[1] + nrm[2]*nrm[2] );
    if ( cross<1.e-20 ) {
      // degenerate face (zero-area): the element cannot produce a
      // section (distorted or collapsed element)
      if ( !warned_length ) {
        pri( "Warning: post_calcul -materi_stress -force: an element "
             "has a degenerate face - no forces calculated for it" );
        warned_length = 1;
      }
      return;
    }
    array_multiply( nrm, nrm, 1./cross, 3 );
    // orientation: away from the centroid (face mid node average)
    array_set( mid, 0., MDIM );
    for ( j=0; j<nface; j++ ) {
      for ( i=0; i<ndim; i++ )
        mid[i] += coords[side_nodes[iside][j]*MDIM+i]/nface;
    }
    d = 0.;
    for ( i=0; i<ndim; i++ ) d += nrm[i]*( mid[i]-centroid[i] );
    if ( d<0. ) {
      nrm[0] = -nrm[0]; nrm[1] = -nrm[1]; nrm[2] = -nrm[2];
    }
    for ( i=0; i<MDIM; i++ ) side_nrm[iside][i] = nrm[i];
  }

  // the "exactly 4 sides consistent with the specified direction"
  // sanity check (manual 6.913): exactly 4 faces with |n*dir| < eps
  nper = 0;
  for ( iside=0; iside<nside; iside++ ) {
    d = scalar_dabs( side_nrm[iside][0]*dir[0]
                   + side_nrm[iside][1]*dir[1]
                   + side_nrm[iside][2]*dir[2] );
    if ( d<eps ) nper++;
  }
  if ( nper!=4 ) {
    if ( !warned_4sides ) {
      pri( "Warning: post_calcul -materi_stress -force: an element has "
           "not exactly 4 sides consistent with the specified direction "
           "(manual Professional 6.913) - no forces calculated for it" );
      warned_4sides = 1;
    }
    return;
  }

  // candidates = the faces not neglected by the direction filter
  // (exclude: |n*dir| <= 1-eps kept; include: |n*dir| >= eps kept)
  ncand = 0;
  for ( iside=0; iside<nside; iside++ ) {
    d = scalar_dabs( side_nrm[iside][0]*dir[0]
                   + side_nrm[iside][1]*dir[1]
                   + side_nrm[iside][2]*dir[2] );
    if ( exclude ) {
      if ( d<=1.-eps ) cand[ncand++] = iside;
    }
    else {
      if ( d>=eps ) cand[ncand++] = iside;
    }
  }

  // end faces = the 2 opposing candidates with the SMALLEST |n*t_hat|
  // (the cross-section faces where the section forces act)
  if ( ncand==2 ) {
    iend = cand[0];
    jend = cand[1];
  }
  else if ( ncand>2 ) {
    for ( i=0; i<ncand; i++ ) order[i] = i;
    for ( i=0; i<ncand; i++ ) {
      for ( j=i+1; j<ncand; j++ ) {
        s1 = scalar_dabs( side_nrm[cand[order[j]]][0]*t_global[0]
                        + side_nrm[cand[order[j]]][1]*t_global[1]
                        + side_nrm[cand[order[j]]][2]*t_global[2] );
        s2 = scalar_dabs( side_nrm[cand[order[i]]][0]*t_global[0]
                        + side_nrm[cand[order[i]]][1]*t_global[1]
                        + side_nrm[cand[order[i]]][2]*t_global[2] );
        if ( s1<s2 ) { itmp = order[i]; order[i] = order[j]; order[j] = itmp; }
      }
    }
    iend = cand[order[0]];
    jend = cand[order[1]];
    // ambiguity: the 2nd end face competes with the 3rd candidate
    if ( ncand>2 ) {
      s1 = scalar_dabs( side_nrm[jend][0]*t_global[0]
                      + side_nrm[jend][1]*t_global[1]
                      + side_nrm[jend][2]*t_global[2] );
      s2 = scalar_dabs( side_nrm[cand[order[2]]][0]*t_global[0]
                      + side_nrm[cand[order[2]]][1]*t_global[1]
                      + side_nrm[cand[order[2]]][2]*t_global[2] );
      tol = 1.e-6*scalar_dabs( s1 );
      if ( tol<1.e-12 ) tol = 1.e-12;
      if ( s2-s1<tol ) {
        if ( !warned_ambiguous ) {
          pri( "Warning: post_calcul -materi_stress -force: ambiguous "
               "end-face selection for an element (distorted geometry or "
               "a reference point on the element diagonal) - no forces "
               "calculated for it" );
          warned_ambiguous = 1;
        }
        return;
      }
    }
  }
  else {
    // no candidates: the direction filter dropped everything
    if ( !warned_ambiguous ) {
      pri( "Warning: post_calcul -materi_stress -force: no element sides "
           "consistent with the specified direction - no forces "
           "calculated for this element" );
      warned_ambiguous = 1;
    }
    return;
  }
  // the two end faces must point to OPPOSITE sides of the element
  // (their outward normals must have a negative dot product). For the
  // straight elements (sheet pile: the two +-height sections) the
  // normals are anti-parallel (dot = -1); for the CURVED ring elements
  // (tunnel: the two +-circumferential faces of a 45-degree element)
  // they are not parallel but still point to opposite sides (dot =
  // -cos(22.5) = -0.92 for the test mesh). The tolerance 0.1 accepts
  // the curved pairs and rejects pairs of adjacent faces (dot > 0).
  d = side_nrm[iend][0]*side_nrm[jend][0]
    + side_nrm[iend][1]*side_nrm[jend][1]
    + side_nrm[iend][2]*side_nrm[jend][2];
  if ( d>0.1 ) {
    if ( !warned_opposing ) {
      pri( "Warning: post_calcul -materi_stress -force: the two end "
           "faces of an element do not point to opposite sides (manual "
           "Professional 6.908) - no forces calculated for it" );
      warned_opposing = 1;
    }
    return;
  }

  // per-face frame: t = the in-face edge direction with the shortest
  // physical extent (thickness_switch -yes -> the longest), oriented
  // away from the reference point; l = n x t (unit); l = the element
  // size in the length direction (manual 6.913: the nodal coordinate
  // difference in the length direction)
  for ( iside=0; iside<2; iside++ ) {
    long int fs = ( iside==0 ? iend : jend );
    for ( i=0; i<ndim; i++ ) {
      e1[i] = coords[corners[fs][1]*MDIM+i]
            - coords[corners[fs][0]*MDIM+i];
      e2[i] = coords[corners[fs][3]*MDIM+i]
            - coords[corners[fs][0]*MDIM+i];
    }
    len1 = array_size( e1, ndim );
    len2 = array_size( e2, ndim );
    if ( len1<1.e-12 || len2<1.e-12 ) {
      if ( !warned_length ) {
        pri( "Warning: post_calcul -materi_stress -force: degenerate "
             "end face of an element - no forces calculated for it" );
        warned_length = 1;
      }
      return;
    }
    if ( ( len1<len2 && thickness_switch!=-YES )
      || ( len1>=len2 && thickness_switch==-YES ) ) {
      array_multiply( e1, e1, 1./len1, ndim );
      for ( i=0; i<MDIM; i++ ) face_t[iside][i] = e1[i];
      array_multiply( e2, e2, 1./len2, ndim );
    }
    else {
      array_multiply( e2, e2, 1./len2, ndim );
      for ( i=0; i<MDIM; i++ ) face_t[iside][i] = e2[i];
      array_multiply( e1, e1, 1./len1, ndim );
    }
    // orient t away from the reference point (plot direction)
    d = 0.;
    for ( i=0; i<ndim; i++ )
      d += face_t[iside][i]*( centroid[i]-reference_point[i] );
    if ( d<0. ) {
      for ( i=0; i<MDIM; i++ ) face_t[iside][i] = -face_t[iside][i];
    }
    else if ( d<1.e-12 ) {
      // reference point on the face plane normal line through the
      // centroid: keep the raw direction (documented edge case)
      if ( !warned_orient ) {
        pri( "Warning: post_calcul -materi_stress -force: the reference "
             "point lies on the thickness line through the centroid of "
             "an element - the plot direction of its values is kept "
             "geometric (documented edge case)" );
        warned_orient = 1;
      }
    }
    // l = n x t (the length direction in the face)
    face_l[iside][0] = side_nrm[fs][1]*face_t[iside][2]
                     - side_nrm[fs][2]*face_t[iside][1];
    face_l[iside][1] = side_nrm[fs][2]*face_t[iside][0]
                     - side_nrm[fs][0]*face_t[iside][2];
    face_l[iside][2] = side_nrm[fs][0]*face_t[iside][1]
                     - side_nrm[fs][1]*face_t[iside][0];
    array_normalize( face_l[iside], 3 );
    // l = the element size in the length direction (the projected
    // extent of the PHYSICAL edges on l; e1/e2 are now unit vectors
    // after the t selection above, so the projection gives the
    // direction alignment; the physical size comes from the original
    // extents len1/len2 of the edge that lies along l)
    s1 = scalar_dabs( e1[0]*face_l[iside][0] + e1[1]*face_l[iside][1]
                    + e1[2]*face_l[iside][2] );
    s2 = scalar_dabs( e2[0]*face_l[iside][0] + e2[1]*face_l[iside][1]
                    + e2[2]*face_l[iside][2] );
    face_len[iside] = ( s1>s2 ? len1*s1 : len2*s2 );
    if ( face_len[iside]<1.e-12 ) {
      if ( !warned_length ) {
        pri( "Warning: post_calcul -materi_stress -force: the length of "
             "an element in length direction is zero - no forces "
             "calculated for it" );
        warned_length = 1;
      }
      return;
    }
    // the middle of the face (the mean of its nodes): the reference
    // for the moment arms dt and dl ("the middle of the element" of
    // the manual 6.913; for the curved ring elements the corner-mean
    // element centroid lies on the chord and would add a spurious
    // moment - verified in msf_tunnel3d)
    for ( j=0; j<nface; j++ ) face_nodes[iside][j] = side_nodes[fs][j];
    array_set( mid, 0., MDIM );
    for ( j=0; j<nface; j++ ) {
      for ( i=0; i<ndim; i++ )
        mid[i] += coords[face_nodes[iside][j]*MDIM+i]/nface;
    }
    msf_integrate_face_3d( npol, face_nodes[iside], coords, sigmas,
      side_nrm[fs], face_t[iside], face_l[iside], mid,
      face_len[iside], face_nor[iside], face_she[iside], face_m1[iside],
      face_m2[iside] );
  }

  // role of the current node
  is_face_node = 0;
  for ( j=0; j<2 && !is_face_node; j++ ) {
    for ( i=0; i<nface; i++ ) {
      if ( face_nodes[j][i]==inod_pos ) {
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
      tol = 1.e-8*array_size( coords, MDIM*nnol );
      if ( max_dist-d>tol ) return; // not an outer node
    }
    nor = face_nor[iface];
    she = face_she[iface];
    mom1 = face_m1[iface];
    mom2 = face_m2[iface];
    got_value = 1;
  }
  else if ( npol==3 && average==-YES && outer!=-YES ) {
    // hex27 middle-plane node (NOT on any end face): the average of
    // the two end faces (manual 6.908; with outer -yes the averaged
    // nodes receive nothing - documented decision)
    nor = 0.5*( face_nor[0] + face_nor[1] );
    she = 0.5*( face_she[0] + face_she[1] );
    mom1 = 0.5*( face_m1[0] + face_m1[1] );
    mom2 = 0.5*( face_m2[0] + face_m2[1] );
    got_value = 1;
    is_averaged = 1;
  }
  else {
    return; // hex8 (all nodes on the end faces) or average -no: done
  }

  // plot components (global x/y/z in the thickness direction; the s
  // component is the physical SIZE; nor/mom keep their sign in the
  // vector direction, she is always positive - manual 6.913). Primary
  // nodes: the face's t; averaged nodes: the element t_hat.
  if ( is_face_node ) {
    for ( i=0; i<MDIM; i++ ) tng[i] = face_t[iface][i];
  }
  else {
    for ( i=0; i<MDIM; i++ ) tng[i] = t_global[i];
  }
  node_values[0]  = nor*tng[0];
  node_values[1]  = nor*tng[1];
  node_values[2]  = nor*tng[2];
  node_values[3]  = ( nor<0. ? -nor : nor );
  node_values[4]  = she*tng[0];
  node_values[5]  = she*tng[1];
  node_values[6]  = she*tng[2];
  node_values[7]  = she;
  node_values[8]  = mom1*tng[0];
  node_values[9]  = mom1*tng[1];
  node_values[10] = mom1*tng[2];
  node_values[11] = ( mom1<0. ? -mom1 : mom1 );
  node_values[12] = mom2*tng[0];
  node_values[13] = mom2*tng[1];
  node_values[14] = mom2*tng[2];
  node_values[15] = ( mom2<0. ? -mom2 : mom2 );
  // plot_switch -yes: invert the drawing direction of the item vector
  for ( j=0; j<4; j++ ) {
    if ( plot_switch[j]==-YES ) {
      node_values[j*4+0] = -node_values[j*4+0];
      node_values[j*4+1] = -node_values[j*4+1];
      node_values[j*4+2] = -node_values[j*4+2];
    }
  }
}

// Per-node 3D calculation (called from post_calcul_materi_stress_force
// inside the parallel node loop of parallel_calcul_node, calcul.cc):
// scans the target element groups, collects the contributions of every
// element that contains the node (averaging identical contributions of
// a conforming mesh), writes the 16 items to result[] and the per-node
// averaged flag (POST_CALCUL_MATERI_STRESS_FORCE_AVERAGED_NODE).
static void msf_calculate_node_3d( long int inod, double result[] )

{
  long int idum[1], ldum=0, ngroups=0, nvalues=0, i=0, ig=0,
    element=0, max_element=0, length_el=0, element_group=0, name=0,
    average=-YES, outer=-NO, plot_switch[4], *thick_switch=NULL,
    ncontrib=0, has_avg=0, has_prim=0, got=0, is_avg=0, value=-NO,
    length=0, *groups=NULL, *el=NULL;
  double ddum[1], reference[DATA_ITEM_SIZE], node_values[16], sum[16];

  array_set( sum, 0., 16 );
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
  db( POST_CALCUL_MATERI_STRESS_FORCE_AVERAGE, 0, &average, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  db( POST_CALCUL_MATERI_STRESS_FORCE_OUTER, 0, &outer, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  plot_switch[0] = plot_switch[1] = plot_switch[2] = plot_switch[3] = -NO;
  db( POST_CALCUL_MATERI_STRESS_FORCE_PLOT_SWITCH, 0, plot_switch, ddum,
    nvalues, VERSION_NORMAL, GET_IF_EXISTS );
  // thickness_switch: one switch per element group (manual 6.917)
  thick_switch = NULL;
  if ( db_active_index(
       POST_CALCUL_MATERI_STRESS_FORCE_THICKNESS_SWITCH, 0,
       VERSION_NORMAL ) ) {
    thick_switch = get_new_int( ngroups<1 ? 1 : ngroups );
    db( POST_CALCUL_MATERI_STRESS_FORCE_THICKNESS_SWITCH, 0,
      thick_switch, ddum, nvalues, VERSION_NORMAL, GET );
  }

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
    if ( name!=-HEX8 && name!=-HEX27 ) continue; // validated earlier
    for ( i=1; i<length_el; i++ ) if ( el[i]==inod ) break;
    if ( i>=length_el ) continue; // the node is not in this element
    msf_element_contribution_3d( element, name, element_group, inod,
      reference+ig*ndim,
      ( thick_switch ? thick_switch[ig] : -NO ),
      average, outer, plot_switch, node_values, got, is_avg );
    if ( got ) {
      for ( i=0; i<16; i++ ) sum[i] += node_values[i];
      ncontrib++;
      if ( is_avg ) has_avg = 1;
      else          has_prim = 1;
    }
  }
  delete[] el;
  if ( thick_switch ) delete[] thick_switch;
  delete[] groups;

  if ( ncontrib>0 ) {
    for ( i=0; i<16; i++ ) result[i] = sum[i]/ncontrib;
    // averaged flag: -YES when the node received (also) middle-plane
    // averages; -primary (the print) skips those nodes
    value = ( has_avg ? -YES : -NO );
  }
  else {
    for ( i=0; i<16; i++ ) result[i] = 0.;
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
// 3D (hex8/hex27): the LOT 3 integration (msf_calculate_node_3d).
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
    msf_calculate_node_3d( inod, result );
  }
  length_result = nitems;
}
