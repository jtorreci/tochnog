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
// flag consumed by the -primary print. LOT 3 (commit fa7650b): the 3D
// integration (hex8/hex27) - face selection from direction_exclude/
// direction_include, per-face thickness/length directions, 2D face
// quadrature, mom1/mom2, hex27 middle-plane averaging. The 3D design
// decisions are documented in
// ProjectDocs/manual-developer/post_calcul_materi_stress_force.md.
// LOT 4 (commit 631258e): the section stress SOURCE changed from the
// recovered nodal stresses (NODE_DOF, the L2 decision) to the ELEMENT
// integration-point stresses (ELEMENT_DOF, the constitutive stresses
// the element actually used - "the element forces needed for this
// option are setup in a timestep", manual Professional 6.913). The
// stress at a section quadrature point is reconstructed from the
// element IP field with the Lagrange polynomials of the element's own
// quadrature grid: for the node-containing Lobatto rules (quad4
// corners, quad9/hex27 nodes, hex8 corners) the section points
// COINCIDE with the element IPs on the face (direct read); for the
// interior Gauss rules (the SRI quad4 2x2 Gauss, MINIMAL rules) the
// section points are interpolated/extrapolated from the IP field.
// Axisymmetric: the section integrand now carries the physical
// circumference 2*PI*r at the section point (the manual's "integrated
// over the thickness" of a ring), so the per-unit-length values match
// the plane-2D dimension (a force per unit circumference).
// LOT 5 (this lot): the section forces are NO LONGER an integral of
// ANY raw stress field over the section (the measured evidence: the
// N/V pollution of the mixed u-sigma scheme lives IN the sigma field
// itself - the nodal values are the exact IP averages, both polluted;
// the Professional is exact because its statics do not come from a
// raw sigma field). The section forces are computed from the ELEMENT
// INTERNAL FORCES f_elem = int B^T*sigma dV (the consistent nodal
// forces of the IP stress field - the equilibrium/equilibrium-based
// section method of the Professional, manual 6.913). The section
// resultant over an end face = the sum of the internal forces of the
// face nodes = int sigma*n_hat dA of the face: for the converged
// solve the internal forces are in equilibrium with the applied loads
// BY CONSTRUCTION (sum_elements f_elem = -P at the free dofs), so the
// section forces are the EXACT statics of the loads (N/V/M ~ 1.000x
// in the arness), not the polluted field integral. The face
// quadrature and the Lagrange stress reconstruction of LOT 4 are
// retired (the stress is read at the element IPs only, the same
// ELEMENT_DOF source); the B matrix and the IP volumes replicate
// pol()/materi() (documented in the helpers).

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
  //   - materi_stress must be a solved dof AND options_element_dof -yes
  //     (the default): the integration reads the element integration-
  //     point stresses (ELEMENT_DOF, LOT 4)
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

  // the integration reads the solved stresses: since LOT 4 the source
  // is the ELEMENT integration-point record (ELEMENT_DOF), which
  // requires materi_stress solved AND options_element_dof -yes (the
  // default; elem.cc maintains the record only then)
  if ( stres_indx<0 ) {
    pri( "Error: post_calcul -materi_stress -force requires materi_stress "
         "in the initia section (the stresses are read from the element "
         "integration points)" );
    exit(TN_EXIT_STATUS);
  }
  if ( options_element_dof!=-YES ) {
    pri( "Error: post_calcul -materi_stress -force requires "
         "options_element_dof -yes (the default) in the initia section: "
         "the section forces are integrated from the element "
         "integration-point stresses (ELEMENT_DOF)" );
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
//   2. STRESS SOURCE (LOT 4+5): the ELEMENT integration-point stresses
//      (ELEMENT_DOF, element_dof_npoint[ipoint*nuknwn+i] with the
//      stress at stres_indx, the constitutive stresses of the last
//      element_loop - the "element forces ... setup in a timestep" of
//      the manual 6.913). LOT 4 read the stress at section quadrature
//      points via the Lagrange reconstruction from the element IP
//      grid; LOT 5 (this lot) reads the stress at the element IPs
//      ONLY and builds the element internal forces f_elem = int B^T*sigma
//      dV (see note 9) - the section resultants are the equilibrium
//      sums over the face nodes, not a stress field integral. The OLD
//      source (the recovered nodal stresses, NODE_DOF - the L2
//      decision) is retired except as the documented FALLBACK for the
//      all-zero ELEMENT_DOF blocks (measured: the SOLVED 3D models
//      with the `derivatives` keyword - gforce10/gforce13). The
//      measured evidence of LOT 4: the nodal values are the exact
//      element-IP averages at the shared nodes, so the N/V pollution
//      lives in the sigma FIELD itself; the element IP stresses are
//      the ones CONSISTENT with the stiffness matrix (the fixed point
//      of the staggered scheme, DIAG 12.2). The validation requires
//      options_element_dof -yes (the default; the record is
//      maintained by elem.cc).
//   3. QUADRATURE: retired from the section integration (LOT 5); the
//      element's own rule (msf_element_rule) is used ONLY to build the
//      internal-force integral (note 9). Historical (LOT 2-4): 1D
//      Gauss(2)/Lobatto(3) along the side/face.
//   4. VALUES (per unit length l; plane 2D l=1, axisymmetric
//      l = 2*PI*radial coordinate of the element centroid, manual
//      6.911: "In a axi-symmetric 2D calculation, the length of the
//      elements is set to 2*PI*radius by Tochnog"):
//      nor = n_hat . R_face / l (SIGNED: positive = tension), she =
//      |t_hat . R_face| / l (only the size, manual 6.913), mom =
//      sum_inod (n_hat . f_inod) * arm_inod / l with arm_inod =
//      (x_inod - x_mid) . t_hat measured from the MIDDLE OF THE FACE
//      (the LOT 5 equivalents of the LOT 2-4 integrals; R_face = the
//      sum of the internal forces of the face nodes - see note 9).
//      The axisymmetric f_elem carries the physical circumference
//      2*PI*r in the IP volumes (the ring section area 2*PI*r*ds),
//      divided by l = 2*PI*r_centroid: the per-unit-length values are
//      the PHYSICAL section forces per unit circumference.
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
//      contributions (identical for a conforming mesh: the face
//      resultants of the two elements at a shared section are the SAME
//      equilibrium value - the face normals and the internal-force
//      resultants flip together, verified in LOT 5).
//   7. outer -yes: only the PRIMARY nodes at the maximum distance from
//      the reference point receive values (manual 6.915: "the nodes
//      which have the furthest distance relative to the reference
//      point"); the averaged nodes get nothing (documented decision).
//   8. AVERAGED FLAG: the per-node record
//      POST_CALCUL_MATERI_STRESS_FORCE_AVERAGED_NODE (-YES when the
//      node received middle-plane averaged values) is written in
//      VERSION_NORMAL and consumed by msf_node_is_averaged() of the
//      -primary print (print_materi_stress_force.cc).
//   9. ELEMENT INTERNAL FORCES (LOT 5): f_elem[inod, idim] =
//      sum_ip volume[ip] * (B^T*sigma)[ip, inod, idim] - the
//      consistent nodal forces of the IP stress field, the same
//      quantity materi() accumulates in element_rhside with the
//      OPPOSITE sign (element_rhside -= volume*force). The kinematics
//      replicate pol()/materi(): shape functions + derivatives at the
//      element IPs (interpolation_polynomial), the physical
//      derivatives dn = invJ*p, the IP volume w*4*|detJ| (2D) /
//      w*8*|detJ| (3D) with the 2*PI*r axisymmetric factor, and the B
//      matrix of polynom.cc:549-599. PHYSICS: the internal forces are
//      in equilibrium with the applied loads BY CONSTRUCTION for the
//      converged solve (sum over all elements = -P at the free dofs;
//      the sum over ONE element's nodes is the rigid-translation
//      identity, identically zero). The section resultant over an end
//      face = sum of the internal forces of the face nodes =
//      int sigma*n_hat dA of the face (the consistent nodal force
//      identity) = the applied-load resultant of the free body: EXACT
//      statics, the same family the Professional uses ("the element
//      forces needed for this option are setup in a timestep"). The
//      moment of the consistent nodal forces about the face middle is
//      the EXACT discrete equivalent of int sigma_nn*dt dA (the shear
//      components contribute only to the in-face component of the
//      moment vector, which is not used). SRI quad4 (opt-in): the
//      full-rule shear part of f_elem is replaced by the reduced
//      1-point shear internal force (the momentum feedback of the
//      fixed point K_SRI*u = P): 4*detJ_c*b_shear*mean(IP sigma_xy)
//      (the shear modulus cancels: sigma_xy = G*gamma, the 2x2 Gauss
//      mean of the bilinear gamma = the centroid value exactly).
//      Without the correction the SRI section forces would read the
//      K_full (locked) equilibrium instead of the K_SRI one.
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

// ---------------------------------------------------------------------
// LOT 4+5: the element integration-point stress field + the element
// internal forces. Helpers.
// ---------------------------------------------------------------------

// The element's per-direction integration rule, replicated from pol()
// (polynom.cc:388-459) for the isoparametric section elements so the
// section forces use EXACTLY the quadrature the element integrated
// with:
//   - npol from the element name; integration_points default -MAXIMAL
//     (npol points per direction); -MINIMAL override = npol-1 points
//     (Gauss, interior);
//   - integration_method default -LOBATTO (the node-containing rules);
//     -GAUSS override;
//   - the SRI quad4/hex8 switches the FULL rule to GAUSS (2x2 / 2x2x2;
//     polynom.cc:421; the same single source sri_active);
//   - axisymmetric + materi_velocity forces GAUSS + MINIMAL (polynom
//     cc:401-407: the 1-point rule at the centroid).
// nper[d] = the number of points in direction d, iso[d][0..nper-1] =
// the iso coordinates and weight[d][0..nper-1] the 1D weights (the
// same arrays pol() uses for the element integration points:
// ipoint = izeta*nper[0]*nper[1] + ieta*nper[0] + ixi in 3D,
// ipoint = ieta*nper[0] + ixi in 2D; the product of the 1D weights
// is the point weight).
static void msf_element_rule( long int element, long int element_group,
  long int name, long int npol, long int nnol, long int nper[],
  double iso[][MPOINT], double weight[][MPOINT] )

{
  long int integration_method=-LOBATTO, integration_points=-MAXIMAL,
    axisymmetric=-NO, ldum=0, idim=0;
  double ddum[1];

  db( GROUP_AXISYMMETRIC, element_group, &axisymmetric, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  if ( name!=-BAR2 ) integration_points = -MAXIMAL;
  if ( axisymmetric==-YES && materi_velocity ) {
    integration_method = -GAUSS;
    integration_points = -MINIMAL;
  }
  db( GROUP_INTEGRATION_POINTS, element_group, &integration_points, ddum,
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( GROUP_INTEGRATION_METHOD, element_group, &integration_method, ddum,
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  if ( integration_points==-NORMAL ) integration_points = -MAXIMAL;
  if ( sri_active( element, element_group, name, nnol ) )
    integration_method = -GAUSS;
  for ( idim=0; idim<ndim; idim++ ) {
    nper[idim] = ( integration_points==-MINIMAL ? npol-1 : npol );
    if      ( integration_method==-GAUSS )
      integration_gauss( nper[idim], iso[idim], weight[idim] );
    else if ( integration_method==-LOBATTO )
      integration_lobatto( nper[idim], iso[idim], weight[idim] );
    else if ( nper[idim]<npol )
      integration_gauss( nper[idim], iso[idim], weight[idim] );
    else
      integration_lobatto( nper[idim], iso[idim], weight[idim] );
  }
}

// ---------------------------------------------------------------------
// LOT 5: the element internal forces f_elem = int B^T*sigma dV (the
// consistent nodal forces of the stress field) and the section
// resultants over an end face. Physics (validated in the arness and
// documented in manual-developer/post_calcul_materi_stress_force.md):
// the internal forces are in equilibrium with the applied loads BY
// CONSTRUCTION for the converged solve (sum_elements f_elem = -P at
// the free dofs), so the section resultant over a face =
// sum of the internal forces of the face nodes = the traction
// resultant int sigma*n_hat dA of the face (exact up to the solver
// tolerance). The Professional's exact statics (1e-10) are consistent
// with this equilibrium/integration ("the element forces needed for
// this option are setup in a timestep", manual 6.913), while the
// section integration over ANY raw sigma field (nodal or IP, lots 2-4)
// carries the pollution of the mixed u-sigma scheme (gforce7 N 1.24x,
// V 2.7x - the pollution lives in the sigma FIELD itself).
// ---------------------------------------------------------------------

// The element internal forces f_elem[inod*ndim+idim] =
// sum_ip volume[ip]*(B^T*sigma)[ip, inod, idim] - the same quantity
// materi() accumulates in element_rhside with the OPPOSITE sign
// (element_rhside -= volume*force, materi.cc:618-619): the internal
// force is the force the element exerts ON its nodes (K*u, opposing
// the applied loads). The kinematics replicate pol() (polynom.cc):
//   - the quadrature rule (nper/iso/weight) comes from msf_element_rule
//     (the SAME rule the element integrated with);
//   - the physical derivatives dn[inod*ndim+idim] = dN_inod/dx_idim =
//     invJ * p (polynom.cc:500-513, the Jacobian J[idim][jdim] =
//     sum_inol p[inoln jdim]*coord[inoln idim]);
//   - the IP volume = w*4*|detJ| (2D) / w*8*|detJ| (3D) with the
//     axisymmetric 2*PI*r factor (polynom.cc:514-526);
//   - the B matrix contracted with sigma (polynom.cc:549-599): in 2D
//     f[inod,0] += vol*(sxx*dN/dx + sxy*dN/dy),
//     f[inod,1] += vol*(sxy*dN/dx + syy*dN/dy); in 3D the 6-component
//     counterpart. The coordinates used are the NORMAL (reference)
//     ones - identical to the new_coord of pol() for the small-strain
//     analyses of the section family (materi_displacement-free or
//     total_linear; documented).
// SRI quad4 (opt-in): the section internal forces are the
// ELEMENT-CONSISTENT ones of the momentum feedback (DIAG-SOLVE-MIXTO
// fix D-c: the fixed point is K_SRI*u = P): the full-rule shear part
// (from the raw Gauss sigma_xy, polluted by the parasitic shear) is
// replaced by the reduced 1-point shear internal force
// sri_g*4*detJ_c*b_shear*(b_shear^T*u) = 4*detJ_c*b_shear*mean(IP
// sigma_xy) (the shear modulus cancels: sigma_xy = G*gamma and the
// 2x2 Gauss mean of the bilinear gamma = the centroid value EXACTLY).
static void msf_element_internal_forces_2d( long int npol, long int nnol,
  double coords[], long int nper[], double iso[][MPOINT],
  double weight[][MPOINT], double sig_ip[], double f_elem[], long int sri,
  long int axisym )

{
  double hx[MPOINT], px[MPOINT], hy[MPOINT], py[MPOINT], xq[MDIM],
    jac[4], invjac[4], detj=0., w=0., vol=0., sxx=0., sxy=0., syy=0.,
    dn[2*MNOL], b_shear[2*MNOL], detj_c=0., mean_xy=0.;
  long int ixi=0, ieta=0, ip=0, inol=0, idim=0;

  for ( inol=0; inol<nnol*ndim; inol++ ) f_elem[inol] = 0.;
  if ( sri ) {
    // the reduced 1-point shear internal force: the centroid B (dN/dy,
    // dN/dx at iso 0,0 - materi.cc:551-553), the centroid detJ and the
    // mean sigma_xy of the 2x2 Gauss points. The centroid Jacobian is
    // computed below at ip=0 (iso 0,0 of the GAUSS 2x2 rule).
    for ( inol=0; inol<nnol; inol++ ) {
      dn[inol*2+0] = 0.25 * ( (inol%2) ? 1. : -1. );
      dn[inol*2+1] = 0.25 * ( (inol/2) ? 1. : -1. );
    }
    array_set( jac, 0., 4 );
    for ( inol=0; inol<nnol; inol++ ) {
      jac[0] += dn[inol*2+0]*coords[inol*MDIM+0];
      jac[1] += dn[inol*2+0]*coords[inol*MDIM+1];
      jac[2] += dn[inol*2+1]*coords[inol*MDIM+0];
      jac[3] += dn[inol*2+1]*coords[inol*MDIM+1];
    }
    detj_c = jac[0]*jac[3] - jac[1]*jac[2];
    if ( detj_c<0. ) detj_c = -detj_c;
    if ( detj_c<1.e-20 ) detj_c = 0.;
    // b_shear[inod,0] = dN/dy, b_shear[inod,1] = dN/dx at the centroid
    // (the invJ pattern of materi.cc:570-576: dN/dy = invjac[2]*dnxi +
    // invjac[3]*dnet, dN/dx = invjac[0]*dnxi + invjac[1]*dnet)
    for ( inol=0; inol<nnol; inol++ ) {
      double detj_s = jac[0]*jac[3]-jac[1]*jac[2];
      b_shear[inol*2+0] = ( -jac[2]*dn[inol*2+0] + jac[0]*dn[inol*2+1] )
                        / detj_s;
      b_shear[inol*2+1] = (  jac[3]*dn[inol*2+0] - jac[1]*dn[inol*2+1] )
                        / detj_s;
    }
    mean_xy = 0.;
    for ( ip=0; ip<nper[0]*nper[1]; ip++ ) mean_xy += sig_ip[ip*3+2];
    mean_xy /= (double)( nper[0]*nper[1] );
  }

  for ( ieta=0; ieta<nper[1]; ieta++ ) {
    interpolation_polynomial( iso[1][ieta], npol, hy, py );
    for ( ixi=0; ixi<nper[0]; ixi++ ) {
      interpolation_polynomial( iso[0][ixi], npol, hx, px );
      ip = ieta*nper[0] + ixi;
      array_set( jac, 0., 4 );
      for ( inol=0; inol<nnol; inol++ ) {
        jac[0] += px[inol%npol]*hy[inol/npol]*coords[inol*MDIM+0];
        jac[1] += px[inol%npol]*hy[inol/npol]*coords[inol*MDIM+1];
        jac[2] += hx[inol%npol]*py[inol/npol]*coords[inol*MDIM+0];
        jac[3] += hx[inol%npol]*py[inol/npol]*coords[inol*MDIM+1];
      }
      detj = jac[0]*jac[3] - jac[1]*jac[2];
      if ( detj<0. ) detj = -detj;
      if ( detj<1.e-20 ) continue;
      invjac[0] =  jac[3]/( jac[0]*jac[3]-jac[1]*jac[2] );
      invjac[1] = -jac[1]/( jac[0]*jac[3]-jac[1]*jac[2] );
      invjac[2] = -jac[2]/( jac[0]*jac[3]-jac[1]*jac[2] );
      invjac[3] =  jac[0]/( jac[0]*jac[3]-jac[1]*jac[2] );
      for ( inol=0; inol<nnol; inol++ ) {
        dn[inol*2+0] = invjac[0]*px[inol%npol]*hy[inol/npol]
                     + invjac[1]*hx[inol%npol]*py[inol/npol];
        dn[inol*2+1] = invjac[2]*px[inol%npol]*hy[inol/npol]
                     + invjac[3]*hx[inol%npol]*py[inol/npol];
      }
      w = weight[0][ixi]*weight[1][ieta];
      vol = w*4.*detj;
      if ( axisym==-YES ) {
        xq[0] = 0.;
        for ( inol=0; inol<nnol; inol++ )
          xq[0] += hx[inol%npol]*hy[inol/npol]*coords[inol*MDIM+0];
        vol *= 2.*PIRAD*scalar_dabs( xq[0] );
      }
      sxx = sig_ip[ip*3+0];
      syy = sig_ip[ip*3+1];
      sxy = ( sri ? 0. : sig_ip[ip*3+2] ); // SRI: reduced 1-point shear
      for ( inol=0; inol<nnol; inol++ ) {
        f_elem[inol*2+0] += vol*( sxx*dn[inol*2+0] + sxy*dn[inol*2+1] );
        f_elem[inol*2+1] += vol*( sxy*dn[inol*2+0] + syy*dn[inol*2+1] );
      }
    }
  }
  if ( sri && detj_c>=1.e-20 ) {
    // the reduced 1-point shear internal force (the momentum feedback
    // of the fixed point K_SRI*u = P; the modulus cancels)
    for ( inol=0; inol<nnol; inol++ ) {
      f_elem[inol*2+0] += 4.*detj_c*b_shear[inol*2+0]*mean_xy;
      f_elem[inol*2+1] += 4.*detj_c*b_shear[inol*2+1]*mean_xy;
    }
  }
  (void)idim;
  (void)idim;
}

// The 3D counterpart of msf_element_internal_forces_2d: the 6 stress
// components (sxx sxy sxz syy syz szz, the msf order) contracted with
// the B rows of polynom.cc:549-599.
// SRI hex8 (opt-in): same correction as the 2D quad4 - the full-rule
// shear part of f_elem (from the raw Gauss sigma_xy/xz/yz) is replaced
// by the reduced 1-point shear internal force of the momentum feedback
// (DIAG-SOLVE-MIXTO fix D-c: the fixed point is K_elem*u = P):
// 8*detJ_c * ( b_xy*mean_xy + b_xz*mean_xz + b_yz*mean_yz ) with the
// centroid B rows and the means of the IP shear stresses (the shear
// moduli cancel: sigma = G*gamma and the 2x2x2 Gauss mean of the
// trilinear gamma = the centroid value). NOTE (measured 2026-08-29):
// the SRI hex8 is singular for loaded configurations (section-warping
// zero-energy modes), so this correction only applies to the stable
// part of the solution - documented in the developer manual.
static void msf_element_internal_forces_3d( long int npol, long int nnol,
  double coords[], long int nper[], double iso[][MPOINT],
  double weight[][MPOINT], double sig_ip[], double f_elem[], long int sri )

{
  double hx[MPOINT], px[MPOINT], hy[MPOINT], py[MPOINT], hz[MPOINT],
    pz[MPOINT], jac[9], invjac[9], detj=0., w=0., vol=0., sig[6],
    dn[3*MNOL], b_xy[3*MNOL], b_xz[3*MNOL], b_yz[3*MNOL],
    jacc[9], invjacc[9], detj_c=0., mean_xy=0., mean_xz=0., mean_yz=0.;
  long int ixi=0, ieta=0, izeta=0, ip=0, inol=0, idim=0, jdim=0;

  for ( inol=0; inol<nnol*ndim; inol++ ) f_elem[inol] = 0.;
  if ( sri ) {
    // the reduced 1-point shear internal force: the centroid B rows
    // (materi.cc SRI hex8 branch), the centroid detJ and the means of
    // the IP shear stresses over the 2x2x2 Gauss points.
    array_set( jacc, 0., 9 );
    for ( inol=0; inol<nnol; inol++ ) {
      double p3[3];
      p3[0] = 0.125 * ( 2.*(inol%2) - 1. );
      p3[1] = 0.125 * ( 2.*((inol/2)%2) - 1. );
      p3[2] = 0.125 * ( 2.*(inol/4) - 1. );
      for ( idim=0; idim<3; idim++ )
        for ( jdim=0; jdim<3; jdim++ )
          jacc[idim*3+jdim] += p3[idim]*coords[inol*MDIM+jdim];
    }
    detj_c = jacc[0]*( jacc[4]*jacc[8] - jacc[5]*jacc[7] )
           - jacc[1]*( jacc[3]*jacc[8] - jacc[5]*jacc[6] )
           + jacc[2]*( jacc[3]*jacc[7] - jacc[4]*jacc[6] );
    if ( detj_c<0. ) detj_c = -detj_c;
    if ( detj_c<1.e-20 ) detj_c = 0.;
    if ( detj_c>=1.e-20 && matrix_inverse( jacc, invjacc, detj_c, 3 ) ) {
      for ( inol=0; inol<nnol; inol++ ) {
        double p3[3], dnx=0., dny=0., dnz=0.;
        p3[0] = 0.125 * ( 2.*(inol%2) - 1. );
        p3[1] = 0.125 * ( 2.*((inol/2)%2) - 1. );
        p3[2] = 0.125 * ( 2.*(inol/4) - 1. );
        for ( jdim=0; jdim<3; jdim++ ) {
          dnx += invjacc[0*3+jdim]*p3[jdim];
          dny += invjacc[1*3+jdim]*p3[jdim];
          dnz += invjacc[2*3+jdim]*p3[jdim];
        }
        b_xy[inol*3+0] = dny;  b_xy[inol*3+1] = dnx;  b_xy[inol*3+2] = 0.;
        b_xz[inol*3+0] = dnz;  b_xz[inol*3+1] = 0.;   b_xz[inol*3+2] = dnx;
        b_yz[inol*3+0] = 0.;   b_yz[inol*3+1] = dnz;  b_yz[inol*3+2] = dny;
      }
    }
    else {
      array_set( b_xy, 0., 3*MNOL );
      array_set( b_xz, 0., 3*MNOL );
      array_set( b_yz, 0., 3*MNOL );
    }
    for ( ip=0; ip<nper[0]*nper[1]*nper[2]; ip++ ) {
      mean_xy += sig_ip[ip*6+1];
      mean_xz += sig_ip[ip*6+2];
      mean_yz += sig_ip[ip*6+4];
    }
    mean_xy /= (double)( nper[0]*nper[1]*nper[2] );
    mean_xz /= (double)( nper[0]*nper[1]*nper[2] );
    mean_yz /= (double)( nper[0]*nper[1]*nper[2] );
  }
  for ( izeta=0; izeta<nper[2]; izeta++ ) {
    interpolation_polynomial( iso[2][izeta], npol, hz, pz );
    for ( ieta=0; ieta<nper[1]; ieta++ ) {
      interpolation_polynomial( iso[1][ieta], npol, hy, py );
      for ( ixi=0; ixi<nper[0]; ixi++ ) {
        interpolation_polynomial( iso[0][ixi], npol, hx, px );
        ip = izeta*nper[0]*nper[1] + ieta*nper[0] + ixi;
        array_set( jac, 0., 9 );
        for ( inol=0; inol<nnol; inol++ ) {
          double p3[3];
          p3[0] = px[inol%npol]*hy[(inol/npol)%npol]*hz[inol/(npol*npol)];
          p3[1] = hx[inol%npol]*py[(inol/npol)%npol]*hz[inol/(npol*npol)];
          p3[2] = hx[inol%npol]*hy[(inol/npol)%npol]*pz[inol/(npol*npol)];
          for ( idim=0; idim<3; idim++ )
            for ( jdim=0; jdim<3; jdim++ )
              jac[idim*3+jdim] += p3[idim]*coords[inol*MDIM+jdim];
        }
        detj = jac[0]*( jac[4]*jac[8] - jac[5]*jac[7] )
             - jac[1]*( jac[3]*jac[8] - jac[5]*jac[6] )
             + jac[2]*( jac[3]*jac[7] - jac[4]*jac[6] );
        if ( !matrix_inverse( jac, invjac, detj, 3 ) ) continue;
        if ( detj<0. ) detj = -detj; // pol() takes |detJ| (line 514)
        if ( detj<1.e-20 ) continue;
        for ( inol=0; inol<nnol; inol++ ) {
          double p3[3];
          p3[0] = px[inol%npol]*hy[(inol/npol)%npol]*hz[inol/(npol*npol)];
          p3[1] = hx[inol%npol]*py[(inol/npol)%npol]*hz[inol/(npol*npol)];
          p3[2] = hx[inol%npol]*hy[(inol/npol)%npol]*pz[inol/(npol*npol)];
          for ( idim=0; idim<3; idim++ ) {
            dn[inol*3+idim] = 0.;
            for ( jdim=0; jdim<3; jdim++ )
              dn[inol*3+idim] += invjac[idim*3+jdim]*p3[jdim];
          }
        }
        w = weight[0][ixi]*weight[1][ieta]*weight[2][izeta];
        vol = w*8.*detj;
        for ( idim=0; idim<6; idim++ ) sig[idim] = sig_ip[ip*6+idim];
        if ( sri ) sig[1] = sig[2] = sig[4] = 0.; // reduced 1-point shear
        for ( inol=0; inol<nnol; inol++ ) {
          f_elem[inol*3+0] += vol*( sig[0]*dn[inol*3+0]
                                  + sig[1]*dn[inol*3+1]
                                  + sig[2]*dn[inol*3+2] );
          f_elem[inol*3+1] += vol*( sig[1]*dn[inol*3+0]
                                  + sig[3]*dn[inol*3+1]
                                  + sig[4]*dn[inol*3+2] );
          f_elem[inol*3+2] += vol*( sig[2]*dn[inol*3+0]
                                  + sig[4]*dn[inol*3+1]
                                  + sig[5]*dn[inol*3+2] );
        }
      }
    }
  }
  if ( sri && detj_c>=1.e-20 ) {
    // the reduced 1-point shear internal force (the momentum feedback
    // of the fixed point K_elem*u = P; the moduli cancel)
    for ( inol=0; inol<nnol; inol++ ) {
      f_elem[inol*3+0] += 8.*detj_c*( b_xy[inol*3+0]*mean_xy
                                    + b_xz[inol*3+0]*mean_xz );
      f_elem[inol*3+1] += 8.*detj_c*( b_xy[inol*3+1]*mean_xy
                                    + b_yz[inol*3+1]*mean_yz );
      f_elem[inol*3+2] += 8.*detj_c*( b_xz[inol*3+2]*mean_xz
                                    + b_yz[inol*3+2]*mean_yz );
    }
  }
}

// The section resultants over ONE element side (2D): the side is the
// CROSS-SECTION face of the structure and the resultants come from the
// ELEMENT INTERNAL FORCES f_elem (LOT 5), not from a stress field
// integration: R_face = sum of the internal forces of the face nodes =
// int sigma*n_hat dA of the face (the consistent nodal force identity;
// for the converged solve = the applied-load resultant of the free
// body, EXACT up to the solver tolerance). Values per unit length l
// (see the design notes above):
//   nor = n_hat . R_face / l   (signed, tension positive)
//   she = |t_hat . R_face| / l (always positive)
//   mom = sum_inod (n_hat . f_inod) * arm_inod / l (signed), with the
//         arm_inod = (x_inod - x_mid) . t_hat measured from the MIDDLE
//         OF THE FACE (the L3 decision; the discrete equivalent of
//         int sigma_nn*dt dA - the moment of the consistent nodal
//         forces about the face middle is EXACT, the shear components
//         contribute only to the in-face component of the moment
//         vector, which is not used - manual 6.913).
static void msf_face_resultants_2d( long int npol, long int face_nodes[],
  double coords[], double f_elem[], double nrm[], double tng[],
  double thick[], double l, double &nor, double &she, double &mom )

{
  double rx=0., ry=0., mid[MDIM], arm=0.;
  long int i=0, inol=0, idim=0;

  array_set( mid, 0., MDIM );
  for ( i=0; i<npol; i++ )
    for ( idim=0; idim<ndim; idim++ )
      mid[idim] += coords[face_nodes[i]*MDIM+idim]/npol;
  rx = ry = 0.;
  for ( i=0; i<npol; i++ ) {
    inol = face_nodes[i];
    rx += f_elem[inol*ndim+0];
    ry += f_elem[inol*ndim+1];
  }
  nor = ( nrm[0]*rx + nrm[1]*ry )/l;
  she = ( tng[0]*rx + tng[1]*ry );
  if ( she<0. ) she = -she;
  she /= l;
  mom = 0.;
  for ( i=0; i<npol; i++ ) {
    inol = face_nodes[i];
    arm = 0.;
    for ( idim=0; idim<ndim; idim++ )
      arm += ( coords[inol*MDIM+idim] - mid[idim] )*thick[idim];
    mom += ( nrm[0]*f_elem[inol*ndim+0]
           + nrm[1]*f_elem[inol*ndim+1] )*arm;
  }
  mom /= l;
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
    icorner=0, ncorner=0, inod_pos=-1, *el=NULL, npoint_ip=0,
    nper[2], ip=0, sri=0;
  double ddum[1], *coord=NULL, *node_dof=NULL, *edof=NULL,
    coords[MDIM*MNOL], sig_ip[3*MPOINT],
    sig_n[3*MNOL], iso[2][MPOINT], wrule[2][MPOINT],
    f_elem[2*MNOL], nrm[MDIM], tng[MDIM],
    centroid[MDIM], thick[MDIM], side_nrm[4][MDIM], side_tng[4][MDIM],
    side_score[4], face_nor[2], face_she[2], face_mom[2],
    max_dist=0., d=0., nor=0., she=0., mom=0., l=1., tol=0., sig_sum=0.;
  long int side_nodes[4][3], face_nodes[2][3];
  static long int warned_centroid=0, warned_ambiguous=0, warned_axisym=0,
    warned_edof=0;
  long int c=0;

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

  // element nodes: coordinates + the recovered nodal stresses (the
  // fallback source, see below)
  el = get_new_int(DATA_ITEM_SIZE);
  db( ELEMENT, element, el, ddum, length_el, VERSION_NORMAL, GET );
  nnol = length_el - 1;
  for ( inol=0; inol<nnol; inol++ ) {
    node = el[1+inol];
    coord = db_dbl( NODE, node, VERSION_NORMAL );
    for ( i=0; i<ndim; i++ ) coords[inol*MDIM+i] = coord[i];
    node_dof = db_dbl( NODE_DOF, node, VERSION_NORMAL );
    sig_n[inol*3+0] = node_dof[stres_indx + stress_indx(0,0)*nder];
    sig_n[inol*3+1] = node_dof[stres_indx + stress_indx(1,1)*nder];
    sig_n[inol*3+2] = node_dof[stres_indx + stress_indx(0,1)*nder];
    if ( node==inod ) inod_pos = inol;
  }
  delete[] el;
  if ( inod_pos<0 ) return; // node not in this element (should not happen)

  // the section stress source (LOT 4+5): the ELEMENT integration-point
  // stresses (ELEMENT_DOF, the constitutive stresses the element used;
  // the recovered nodal values of the L2 source are retired). The
  // element's own quadrature rule is replicated (msf_element_rule) and
  // the IP stresses are read in the pol() IP layout
  // ip = ieta*nper[0]+ixi (nuknwn values per IP, the stress at
  // stres_indx). A missing record (options_element_dof -no is rejected
  // in the validation) skips the element with a warning.
  msf_element_rule( element, element_group, name, npol, nnol, nper, iso,
    wrule );
  sri = sri_active( element, element_group, name, nnol );
  npoint_ip = nper[0]*nper[1];
  if ( options_element_dof==-YES &&
       db_active_index( ELEMENT_DOF, element, VERSION_NORMAL ) ) {
    edof = db_dbl( ELEMENT_DOF, element, VERSION_NORMAL );
    for ( ip=0; ip<npoint_ip; ip++ ) {
      sig_ip[ip*3+0] = edof[ip*nuknwn + stres_indx + stress_indx(0,0)*nder];
      sig_ip[ip*3+1] = edof[ip*nuknwn + stres_indx + stress_indx(1,1)*nder];
      sig_ip[ip*3+2] = edof[ip*nuknwn + stres_indx + stress_indx(0,1)*nder];
    }
    // FALLBACK (measured, 2026-08-28): for the SOLVED 3D models with
    // the `derivatives` keyword the staggered element loop does not
    // propagate the strain of the converged velocity into the element
    // integration-point stresses (gforce10/gforce13: the ELEMENT_DOF
    // stress block is all zero while the recovered NODE_DOF carries
    // the correct values). In that case the section reads the
    // recovered nodal stresses evaluated at the element IP positions
    // with the shape functions (identical to the pre-LOT4 source for
    // the node-containing rules; the same values as the element IPs
    // for the single-element faces).
    sig_sum = 0.;
    for ( ip=0; ip<npoint_ip; ip++ )
      sig_sum += scalar_dabs(sig_ip[ip*3+0]) + scalar_dabs(sig_ip[ip*3+1])
               + scalar_dabs(sig_ip[ip*3+2]);
    if ( sig_sum<1.e-12 ) {
      if ( !warned_edof ) {
        pri( "Warning: post_calcul -materi_stress -force: the element "
             "integration-point stresses of an element are all zero (the "
             "staggered element loop has not propagated the strain in "
             "this configuration) - the recovered nodal stresses are "
             "used for it (documented fallback, LOT 4)" );
        warned_edof = 1;
      }
      for ( ip=0; ip<npoint_ip; ip++ ) {
        double hx[MPOINT], px[MPOINT], hy[MPOINT], py[MPOINT];
        long int ixi = ip%nper[0];
        long int ieta = ip/nper[0];
        interpolation_polynomial( iso[0][ixi], npol, hx, px );
        interpolation_polynomial( iso[1][ieta], npol, hy, py );
        for ( c=0; c<3; c++ ) {
          sig_ip[ip*3+c] = 0.;
          for ( inol=0; inol<nnol; inol++ )
            sig_ip[ip*3+c] += hx[inol%npol]*hy[inol/npol]*sig_n[inol*3+c];
        }
      }
    }
  }
  else {
    if ( !warned_edof ) {
      pri( "Warning: post_calcul -materi_stress -force: no ELEMENT_DOF "
           "record for an element (options_element_dof -no or an element "
           "the solver skipped) - no forces calculated for it" );
      warned_edof = 1;
    }
    return;
  }

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

  // the element internal forces (LOT 5): f_elem = int B^T*sigma dV,
  // the consistent nodal forces of the IP stress field - in equilibrium
  // with the applied loads for the converged solve. The face resultants
  // (nor/she/mom) follow from the free body of the face nodes.
  msf_element_internal_forces_2d( npol, nnol, coords, nper, iso, wrule,
    sig_ip, f_elem, sri, axisym );
  for ( j=0; j<2; j++ ) {
    iside = order[j];
    for ( i=0; i<npol; i++ ) face_nodes[j][i] = side_nodes[iside][i];
    msf_face_resultants_2d( npol, face_nodes[j], coords, f_elem,
      side_nrm[iside], side_tng[iside], thick, l,
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
//      reconstructed from the ELEMENT integration-point field at the
//      face quadrature points (LOT 4, the 2D decision 2: the
//      reference coordinates of the quadrature point are interpolated
//      from the local-node reference positions with the face shape
//      functions, then the element IP field is evaluated with the
//      Lagrange weights of the element's own rule - the direct read
//      for the node-containing Lobatto rules of the hex8/hex27).
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

// The section resultants over ONE face of a hex8/hex27 element (3D):
// the face is one of the two END faces of the section and the
// resultants come from the ELEMENT INTERNAL FORCES f_elem (LOT 5, the
// 2D decision): R_face = sum of the internal forces of the face nodes
// = int sigma*n_hat dA of the face (for the converged solve = the
// applied-load resultant of the free body, EXACT up to the solver
// tolerance). nrm = outward normal, thick = the face thickness
// direction t (unit), leng = the face length direction l (unit),
// face_mid = the middle of the face (the mean of its nodes; "the
// middle of the element" of the manual 6.913 in the thickness/length
// directions - the L3 arm decision). Results per unit length l:
//   nor  = n_hat . R_face / l   (signed, tension positive)
//   she  = |t_hat . R_face| / l (always positive)
//   mom1 = sum (n_hat . f_inod) * ((x_inod-face_mid) . t_hat) / l
//   mom2 = sum (n_hat . f_inod) * ((x_inod-face_mid) . l_hat) / l
static void msf_face_resultants_3d( long int npol, long int face_nodes[],
  double coords[], double f_elem[], double nrm[], double thick[],
  double leng[], double l, double &nor, double &she, double &mom1,
  double &mom2 )

{
  double rx=0., ry=0., rz=0., mid[MDIM], arm=0.;
  long int i=0, inol=0, idim=0, nface=npol*npol;

  array_set( mid, 0., MDIM );
  for ( i=0; i<nface; i++ )
    for ( idim=0; idim<ndim; idim++ )
      mid[idim] += coords[face_nodes[i]*MDIM+idim]/nface;
  rx = ry = rz = 0.;
  for ( i=0; i<nface; i++ ) {
    inol = face_nodes[i];
    rx += f_elem[inol*ndim+0];
    ry += f_elem[inol*ndim+1];
    rz += f_elem[inol*ndim+2];
  }
  nor = ( nrm[0]*rx + nrm[1]*ry + nrm[2]*rz )/l;
  she = ( thick[0]*rx + thick[1]*ry + thick[2]*rz );
  if ( she<0. ) she = -she;
  she /= l;
  mom1 = mom2 = 0.;
  for ( i=0; i<nface; i++ ) {
    inol = face_nodes[i];
    arm = 0.;
    for ( idim=0; idim<ndim; idim++ )
      arm += ( coords[inol*MDIM+idim] - mid[idim] )*thick[idim];
    mom1 += ( nrm[0]*f_elem[inol*ndim+0]
            + nrm[1]*f_elem[inol*ndim+1]
            + nrm[2]*f_elem[inol*ndim+2] )*arm;
    arm = 0.;
    for ( idim=0; idim<ndim; idim++ )
      arm += ( coords[inol*MDIM+idim] - mid[idim] )*leng[idim];
    mom2 += ( nrm[0]*f_elem[inol*ndim+0]
            + nrm[1]*f_elem[inol*ndim+1]
            + nrm[2]*f_elem[inol*ndim+2] )*arm;
  }
  mom1 /= l;
  mom2 /= l;
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
    iend=0, jend=0, cand[6], dir_has=0, exclude=0, nface=0,
    npoint_ip=0, ip=0, c=0, nper3[3], sri=0;
  double ddum[1], *coord=NULL, *node_dof=NULL, *edof=NULL,
    coords[MDIM*MNOL], sig_ip[6*MPOINT],
    sig_n[6*MNOL], iso3[3][MPOINT], wrule3[3][MPOINT],
    f_elem[3*MNOL], nrm[MDIM], tng[MDIM],
    centroid[MDIM], t_global[MDIM], side_nrm[6][MDIM],
    face_t[2][MDIM], face_l[2][MDIM], face_len[2], face_nor[2],
    face_she[2], face_m1[2], face_m2[2], max_dist=0., d=0., nor=0.,
    she=0., mom1=0., mom2=0., tol=0., dir[MDIM], eps=1.e-8,
    e1[MDIM], e2[MDIM], mid[MDIM], s1=0., s2=0., len1=0., len2=0.,
    cross=0., sig_sum=0.;
  long int side_nodes[6][9], face_nodes[2][9], corners[6][4];
  static long int warned_4sides=0, warned_centroid=0,
    warned_ambiguous=0, warned_opposing=0, warned_length=0,
    warned_orient=0, warned_edof=0;

  got_value = 0;
  is_averaged = 0;
  for ( i=0; i<16; i++ ) node_values[i] = 0.;

  npol = ( name==-HEX8 ? 2 : 3 );
  nface = npol*npol;
  ncorner = 8;

  // element nodes: coordinates + the recovered nodal stresses (the
  // fallback source, see below)
  el = get_new_int(DATA_ITEM_SIZE);
  db( ELEMENT, element, el, ddum, length_el, VERSION_NORMAL, GET );
  nnol = length_el - 1;
  for ( inol=0; inol<nnol; inol++ ) {
    node = el[1+inol];
    coord = db_dbl( NODE, node, VERSION_NORMAL );
    for ( i=0; i<ndim; i++ ) coords[inol*MDIM+i] = coord[i];
    node_dof = db_dbl( NODE_DOF, node, VERSION_NORMAL );
    for ( c=0; c<6; c++ )
      sig_n[inol*6+c] = node_dof[stres_indx + c*nder];
    if ( node==inod ) inod_pos = inol;
  }
  delete[] el;
  if ( inod_pos<0 ) return; // node not in this element (should not happen)

  // the section stress source (LOT 4): the ELEMENT integration-point
  // stresses (ELEMENT_DOF; the recovered nodal values of the L3 source
  // are retired). The element's own quadrature rule is replicated
  // (msf_element_rule) and the IP stresses are read in the pol() IP
  // layout ip = izeta*nper[0]*nper[1] + ieta*nper[0] + ixi. A missing
  // record (options_element_dof -no is rejected in the validation)
  // skips the element with a warning.
  msf_element_rule( element, element_group, name, npol, nnol, nper3, iso3,
    wrule3 );
  npoint_ip = nper3[0]*nper3[1]*nper3[2];
  if ( options_element_dof==-YES &&
       db_active_index( ELEMENT_DOF, element, VERSION_NORMAL ) ) {
    edof = db_dbl( ELEMENT_DOF, element, VERSION_NORMAL );
    for ( ip=0; ip<npoint_ip; ip++ ) {
      for ( c=0; c<6; c++ )
        sig_ip[ip*6+c] = edof[ip*nuknwn + stres_indx + c*nder];
    }
    // FALLBACK (measured, 2026-08-28): for the SOLVED 3D models with
    // the `derivatives` keyword the staggered element loop does not
    // propagate the strain of the converged velocity into the element
    // integration-point stresses (gforce10/gforce13: the ELEMENT_DOF
    // stress block is all zero while the recovered NODE_DOF carries
    // the correct values). The section then reads the recovered nodal
    // stresses evaluated at the element IP positions with the shape
    // functions (identical to the pre-LOT4 source for the
    // node-containing rules; the same values as the element IPs for
    // the single-element faces).
    sig_sum = 0.;
    for ( ip=0; ip<npoint_ip; ip++ )
      for ( c=0; c<6; c++ ) sig_sum += scalar_dabs(sig_ip[ip*6+c]);
    if ( sig_sum<1.e-12 ) {
      if ( !warned_edof ) {
        pri( "Warning: post_calcul -materi_stress -force: the element "
             "integration-point stresses of an element are all zero (the "
             "staggered element loop has not propagated the strain in "
             "this configuration) - the recovered nodal stresses are "
             "used for it (documented fallback, LOT 4)" );
        warned_edof = 1;
      }
      for ( ip=0; ip<npoint_ip; ip++ ) {
        double hx[MPOINT], px[MPOINT], hy[MPOINT], py[MPOINT],
          hz[MPOINT], pz[MPOINT];
        long int ixi = ip%nper3[0];
        long int ieta = (ip/nper3[0])%nper3[1];
        long int izeta = ip/(nper3[0]*nper3[1]);
        interpolation_polynomial( iso3[0][ixi], npol, hx, px );
        interpolation_polynomial( iso3[1][ieta], npol, hy, py );
        interpolation_polynomial( iso3[2][izeta], npol, hz, pz );
        for ( c=0; c<6; c++ ) {
          sig_ip[ip*6+c] = 0.;
          for ( inol=0; inol<nnol; inol++ )
            sig_ip[ip*6+c] += hx[inol%npol]*hy[(inol/npol)%npol]
              *hz[inol/(npol*npol)]*sig_n[inol*6+c];
        }
      }
    }
  }
  else {
    if ( !warned_edof ) {
      pri( "Warning: post_calcul -materi_stress -force: no ELEMENT_DOF "
           "record for an element (options_element_dof -no or an element "
           "the solver skipped) - no forces calculated for it" );
      warned_edof = 1;
    }
    return;
  }

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
    // moment - verified in msf_tunnel3d). The section resultants come
    // from the ELEMENT INTERNAL FORCES (LOT 5, see below).
    for ( j=0; j<nface; j++ ) face_nodes[iside][j] = side_nodes[fs][j];
    array_set( mid, 0., MDIM );
    for ( j=0; j<nface; j++ ) {
      for ( i=0; i<ndim; i++ )
        mid[i] += coords[face_nodes[iside][j]*MDIM+i]/nface;
    }
  }
  // the element internal forces (LOT 5): f_elem = int B^T*sigma dV,
  // the consistent nodal forces of the IP stress field - in equilibrium
  // with the applied loads for the converged solve. The face resultants
  // (nor/she/mom1/mom2) follow from the free body of the face nodes.
  // SRI hex8 (opt-in): the full-rule shear part is replaced by the
  // reduced 1-point shear internal force (the momentum feedback of the
  // fixed point K_elem*u = P, the same correction as the 2D quad4).
  sri = sri_active( element, element_group, name, nnol );
  msf_element_internal_forces_3d( npol, nnol, coords, nper3, iso3, wrule3,
    sig_ip, f_elem, sri );
  for ( iside=0; iside<2; iside++ ) {
    long int fs = ( iside==0 ? iend : jend );
    msf_face_resultants_3d( npol, face_nodes[iside], coords, f_elem,
      side_nrm[fs], face_t[iside], face_l[iside], face_len[iside],
      face_nor[iside], face_she[iside], face_m1[iside], face_m2[iside] );
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
