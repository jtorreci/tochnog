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
//   - gap (Fase 3, RF-3, CONVERGENCE 2026-09-04): the interface is
//       CLOSED when the accumulated normal strain <= gap, OPEN when
//       strain_normal > gap (manual Professional 6.625: "Only when the
//       sides displacements are such that the normal strain becomes lower
//       then the specified gap value the interface will be closed and
//       start to generate stresses"; an opened interface "does not have
//       stresses"). A physical gap is a NEGATIVE gap value: the interface
//       stays open (no stress) until compression brings the accumulated
//       strain below |gap|. Default without the record = +1.e20 (always
//       closed - "if you want to allow always tension stresses set gap
//       to, by example, 1.e20"). The stress of the closed phase
//       accumulates ONLY the normal strain of the steps that end closed
//       (ELEMENT_INTERFACE_FORCE_NORM history): the free travel of the
//       open phase never enters the stress (interface2 of the corpus:
//       gap 0.1, 200 steps of -1e-3 -> final stress -101 = kn*(-0.101),
//       the 101 closed steps only, NOT kn*(-0.2)). While open the stress
//       is zero and only the residual stiffness keeps the matrix
//       regularized.
//   - tension limit (Fase 3, RF-2): the interface opens in traction when
//       the TOTAL accumulated normal force |kn*strain_normal| exceeds the
//       limit: group_interface_materi_plasti_tension_direct tension_limit
//   - Mohr-Coulomb (Fase 3, RF-1): cumulative. The friction limit applies
//       to the TOTAL tangential force F_t (history ELEMENT_INTERFACE_FORCE_TANG):
//       group_interface_materi_plasti_mohr_coul_direct phi c phi_flow
//       trial = F_t,old + kt*du_tang (the history stores the ELASTIC trial
//       kt*gamma_total, which keeps accumulating across plastic steps),
//       clamped to +/- max_fric with
//       max_fric = |c + Fn*tan(phi)|, Fn = -kn*strain_normal (POSITIVE
//       under compression: verified against the Professional plateau of
//       mohr_coul_direct4: c + |Fn|*tan(phi) = 1.20271). Active by the
//       PRESENCE of the record (phi=0,c=0 -> max_fric=0 -> free sliding).
//       The assembled rhs carries the FULL accumulated forces (spring.cc
//       pattern: node_rhside = w*sigma with sigma = kn*eps_acc for the
//       normal and sigma = clamped trial for the shear - direct3/4 and
//       interface9 verified against the Professional .dbs).
//   - dilatancy (Fase 3, RF-4): when the tangential force plastifies,
//       strain_normal += -dgamma_inc*tan(phi_flow) (plastic normal
//       opening), where dgamma_inc is the plastic slip INCREMENT of the
//       step (the accumulated-trial return multiplier minus the slip
//       already accumulated in the past - otherwise the dilatancy
//       ratchet double-counts, direct3 sigma_n -250 vs -200).
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
  long int idim=0, jdim=0, inol=0, jnol=0, knol=0, indx=0, swit=0, ldum=0, 
    nnol=4, ns1=2, mc_active=0, memory=-UPDATED_LINEAR, idum[1],
    axisymmetric=-NO;
  double damping_iface=0.;
  double dtime=0., kn=0., kt1=0., kt2=0., tmp=0., ddum[1],
    normal[MDIM], tangent[MDIM], tangent2[MDIM], du[MDIM],
    du_norm=0., du_tang=0., du_tang2=0., stress_normal=0., stress_shear=0.,
    stress_shear2=0., strain_normal=0., strain_eff=0.,
    gap=0., tension_limit=0.,
    residual_factor=0.01, phi=0., c=0., phi_flow=0.,
    ddum3[3],
    f_t_old=0., f_t=0., f_t2_old=0., f_t2=0.,
    force_gravity[MDIM];
  long int *nodes=NULL;
  swit = set_swit(element,-1,"interface_element");
  if ( swit ) pri( "In routine INTERFACE_ELEMENT." );

  if ( !db_active_index( GROUP_INTERFACE, element_group, VERSION_NORMAL ) )
    return;

  // axisymmetric (CONVERGENCE 2026-08-30, verified against Professional
  // interface9): the interface is a ring of radius r; the nodal force and
  // stiffness carry the physical weight 2*pi*r (the same pattern as
  // area.cc / elem.cc). The radius is the radial coordinate of the pair.
  db( GROUP_AXISYMMETRIC, element_group, &axisymmetric, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );

  db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET );

  // group parameters
  // ddum3 MUST be zeroed: with GET_IF_EXISTS on a missing record db() does
  // not write dval, and an uninitialized buffer made kn/kt garbage (NaN in
  // the assembled matrix with gcc -O1; exposed 2026-08-24 by the clean
  // rebuild on a newer compiler).
  array_set( ddum3, 0., 3 );
  {
    // variable-length: 2D tests give 2 values (kn, kt), 3D give 3
    // (kn, kt1, kt2) - read only what the record carries (manual
    // Professional: "kn kt,first kt,second" with kt,second omitted
    // in 2D)
    long int length_stiff = db_len(
      GROUP_INTERFACE_MATERI_ELASTI_STIFFNESS, element_group,
      VERSION_NORMAL );
    if ( length_stiff>3 ) length_stiff = 3;
    if ( length_stiff<2 ) length_stiff = 3; // default: read 3 (zeroed)
    db( GROUP_INTERFACE_MATERI_ELASTI_STIFFNESS, element_group, idum,
      ddum3, length_stiff, VERSION_NORMAL, GET_IF_EXISTS );
  }
  kn = ddum3[0]; kt1 = ddum3[1]; kt2 = ddum3[2];
  // group_interface_damping (interface12 of the corpus): viscous damping
  // on the RELATIVE velocity between the sides, d*(v_side2 - v_side1).
  // The RHS carries -d*(v2-v1) and the matrix adds d to the pair blocks
  // (the damping force is proportional to velocity, not displacement, so
  // it does NOT carry the dtime factor of the stiffness). The corpus test
  // interface12: kn=1, kt=0.5, d=1, loads -1/6,-2/3,-1/6 -> all three
  // top nodes reach disy=-0.5 exactly (the damping equalizes the nodal
  // response of the quadratic interface).
  db( GROUP_INTERFACE_DAMPING, element_group, idum, &damping_iface, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
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
  else if ( name==-QUAD6 ) {
    // the quadratic 2D interface (a converted bar3/quad8/quad9:
    // 3 nodes per side, manual Professional control_mesh_convert)
    nnol = 6;
  }
  else if ( name==-PRISM6 ) {
    nnol = 6;
  }
  else if ( name==-HEX18 ) {
    // the quadratic 3D interface (a converted -quad8 facial element of
    // the interface_quad8_hex20 family, or a native -hex18 of the suite
    // like interface3): two -quad9 sides of 9 nodes each (ns1 = 9)
    nnol = 18;
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
    // orientation (CONVERGENCE 2026-08-30, verified against Professional
    // .dbs of interface1/8/13/14/15): the Professional orients the 2D
    // interface normal so that compression gives a NEGATIVE normal strain.
    // The default normal (-t.y, t.x) matches when the first node of side 2
    // has a HIGHER number than the first node of side 1 (interface1:
    // nodes 1,2,3 | 4,5,6 -> normal (0,1); interface13/14/15: 1,2 | 3,4 ->
    // (-0.447,0.894)/(0,1)). When the side-2 nodes are numbered LOWER
    // (interface8 element 3 = quad4 5 6 3 4: side 1 = 5,6, side 2 = 3,4)
    // the Professional flips the normal to (0,-1) - the strain of the
    // compressed interface stays negative. This mirrors the mesh
    // conversion orientation: the copied side keeps the numbering of the
    // source block, so the relative numbering encodes the side order.
    if ( nodes[ns1] < nodes[0] ) {
      normal[0] = -normal[0];
      normal[1] = -normal[1];
    }
    array_set( tangent2, 0., MDIM );
  }
  else {
    // 3D: the two side-1 edges define the surface plane. For the hex18
    // the side-1 nodes are the 9 nodes of a -quad9 face in tensor order
    // (BL,BM,BR,LM,C,RM,TL,TM,TR): the surface corners are the side-1
    // nodes 0, 2 and 6 (BL, BR, TL) - nodes 1 (BM) is a mid-edge node
    // collinear with BL/BR and would give a degenerate cross product.
    long int ic1 = 1, ic2 = 2;
    if ( name==-HEX18 ) { ic1 = 2; ic2 = 6; }
    double *c0, *c1, *c2;
    if ( memory==-TOTAL_LINEAR ) {
      c0 = db_dbl( NODE_START_REFINED, nodes[0], VERSION_NORMAL );
      c1 = db_dbl( NODE_START_REFINED, nodes[ic1], VERSION_NORMAL );
      c2 = db_dbl( NODE_START_REFINED, nodes[ic2], VERSION_NORMAL );
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
  //
  // CONVERGENCE (2026-08-30): the Professional assembles the interface
  // PER INTEGRATION POINT - each facing pair of nodes (i, i+ns1) is an
  // independent spring with its own relative displacement du_i and its
  // own Lobatto weight (quad6: 1/6, 4/6, 1/6; quad4: 1/2, 1/2; bar2: 1;
  // prism6: 1/6,4/6,1/6; hex8: 1/12,5/12,5/12,1/12). The old code
  // averaged the whole side and assembled every node against every node,
  // which produced a rank-1 matrix per side -> SINGULAR for ns1>1
  // (interface1: only element 1, 3 free nodes -> band solver info!=0).
  // Verified against the Professional .dbs of interface1: the applied
  // loads -1,-4,-1 on nodes 4,5,6 ARE the Lobatto weights times the
  // total 6 (the interface nodal force = weight_i * sigma_i).
  ns1 = nnol/2;
  double *w_ip = get_new_dbl( ns1 );
  // Integration weights of the interface side:
  //  2D: Lobatto along the side (quad6: 1/6,4/6,1/6; quad4: 1/2,1/2) -
  //  verified against the Professional nodal loads of interface1
  //  (-1,-4,-1 = weights x total 6).
  //  3D: UNIFORM 1/ns1 (the converted quad4/prism6/hex8 interface is a
  //  surface; the Professional distributes the face traction evenly -
  //  interface_quad4_hex8: uniform 1/4 gives sigzz=-1.0 exactly, while
  //  the 1D Lobatto weights would give -0.2475).
  //  3D hex18 (ns1=9): the -quad9 face is quadratic, and the
  //  Professional distributes the face traction with the 2D product of
  //  the 1D Lobatto rule (interface3 of the corpus loads the top side
  //  with -1 on the corners, -4 on the mid-edge nodes and -16 on the
  //  centre = (1,4,16)/36 per tensor slot; the uniform weights would
  //  give a NON-uniform per-intpnt stress field, e.g. sigma = -20.25 on
  //  the corner intpnt instead of -9).
  if ( ndim==3 && ns1==9 ) {
    static const double lobatto_2d[9] = { 1., 4., 1., 4., 16., 4.,
      1., 4., 1. };
    for ( inol=0; inol<ns1; inol++ )
      w_ip[inol] = lobatto_2d[inol]/36.;
  }
  else if ( ndim==3 ) {
    for ( inol=0; inol<ns1; inol++ ) w_ip[inol] = 1./(double)ns1;
  }
  else if ( ns1==1 ) { w_ip[0] = 1.; }
  else if ( ns1==2 ) { w_ip[0] = 0.5; w_ip[1] = 0.5; }
  else if ( ns1==3 ) { w_ip[0] = 1./6.; w_ip[1] = 4./6.; w_ip[2] = 1./6.; }
  else               { w_ip[0] = 1./12.; w_ip[1] = 5./12.;
                       w_ip[2] = 5./12.; w_ip[3] = 1./12.; }
  // CONVERGENCE (2026-09-04, corpus patch1): the assembled per-pair spring
  // force and stiffness must carry the element MEASURE - the physical
  // integral over the interface is the Lobatto/even sum times the element
  // length (2D line) or face area (3D surface). The historical assembly
  // used the bare weights (sum = 1 = the UNIT-measure element), which is
  // exact only for unit-length/unit-area interfaces (interface1 of the
  // suite: length 1) and silently under-integrates every other element.
  // patch1 exposed it: the inclined interface (quad6 elements of length
  // 1.677 and 0.559 between two loaded quad9 blocks) converges to a
  // NON-UNIFORM equilibrium traction (sigma_n = 1610/1073/536 per intpnt,
  // mean 1431) instead of the Professional's uniform 960: without the
  // length the discrete force system loses the load-path moment arm, so a
  // uniform traction cannot balance the applied edge load (the sum
  // w_i*sigma_i balances 2400 = the load only with the non-uniform field).
  // With the measure, the pair stiffness kn*w_i*L is the spring constant
  // of the tributary length w_i*L and the uniform traction 960 equilibrates
  // (verified: the Pro's nodal loads of interface1 are w_i*6 = w_i*L*sigma
  // with L = 1, so the unit-length validations are untouched).
  double iface_measure = 1.;
  if ( ns1>1 ) {
    if ( ndim==2 ) {
      // side length between the first and the last side-1 node (reference
      // geometry: same memory branch as the interface frame above)
      double *ca, *cb;
      if ( memory==-TOTAL_LINEAR ) {
        ca = db_dbl( NODE_START_REFINED, nodes[0], VERSION_NORMAL );
        cb = db_dbl( NODE_START_REFINED, nodes[ns1-1], VERSION_NORMAL );
      }
      else {
        ca = &coord[0*ndim];
        cb = &coord[(ns1-1)*ndim];
      }
      double dx = cb[0]-ca[0], dy = cb[1]-ca[1];
      iface_measure = sqrt( dx*dx + dy*dy );
    }
    else {
      // face area of the (flat) side-1 polygon: triangle = |e1xe2|/2,
      // quad = |d1xd2|/2 with the diagonals. NOTE the quad4-in-3D
      // (nnol=4, ns1=2: the "glue" quad4 between two solids, e.g.
      // interface_quad4_hex8 of the corpus) has its FOUR nodes as the
      // face corners (0,1,2,3), not ns1 nodes per side.
      // The hex18 side-1 is a -quad9 face in tensor order: its 4 face
      // corners are the side-1 nodes 0, 2, 6 and 8 (BL, BR, TL, TR; the
      // nodes 1,3,5,7 are mid-edge nodes and 4 the centre - the corner
      // polygon read in record order would be a bowtie).
      long int inol_c[4];
      if ( nnol==18 ) { inol_c[0]=0; inol_c[1]=2; inol_c[2]=6; inol_c[3]=8; }
      else {
        inol_c[0]=0; inol_c[1]=1; inol_c[2]=2;
        long int i3 = ns1-1;
        if ( ns1==2 ) i3 = 3;
        inol_c[3]=i3;
      }
      double x0[MDIM], x1[MDIM], x2[MDIM], x3[MDIM];
      for ( idim=0; idim<3; idim++ ) {
        double *cn;
        cn = ( memory==-TOTAL_LINEAR ) ?
          db_dbl( NODE_START_REFINED, nodes[inol_c[0]], VERSION_NORMAL ) :
          &coord[inol_c[0]*ndim];
        x0[idim] = cn[idim];
        cn = ( memory==-TOTAL_LINEAR ) ?
          db_dbl( NODE_START_REFINED, nodes[inol_c[1]], VERSION_NORMAL ) :
          &coord[inol_c[1]*ndim];
        x1[idim] = cn[idim];
        cn = ( memory==-TOTAL_LINEAR ) ?
          db_dbl( NODE_START_REFINED, nodes[inol_c[2]], VERSION_NORMAL ) :
          &coord[inol_c[2]*ndim];
        x2[idim] = cn[idim];
        cn = ( memory==-TOTAL_LINEAR ) ?
          db_dbl( NODE_START_REFINED, nodes[inol_c[3]], VERSION_NORMAL ) :
          &coord[inol_c[3]*ndim];
        x3[idim] = cn[idim];
      }
      double e1[MDIM], e2[MDIM], nrm[MDIM];
      // fan of triangles from the node 0: the triangle (0,1,2) plus, for
      // the quad, the triangle (0,2,3). The crossed-diagonals formula
      // needs the CYCLIC corner order, which the tochnog hex8/quad4 face
      // numbering does not follow (the face 5,6,7,8 of the unit hex8 is a
      // bowtie order: the diagonals coincide and the area would be 0);
      // the triangle fan is order-robust and exact for planar faces.
      for ( idim=0; idim<3; idim++ ) {
        e1[idim] = x1[idim]-x0[idim];
        e2[idim] = x2[idim]-x0[idim];
      }
      nrm[0] = e1[1]*e2[2] - e1[2]*e2[1];
      nrm[1] = e1[2]*e2[0] - e1[0]*e2[2];
      nrm[2] = e1[0]*e2[1] - e1[1]*e2[0];
      iface_measure = 0.5 * sqrt( nrm[0]*nrm[0] + nrm[1]*nrm[1] +
        nrm[2]*nrm[2] );
      if ( ns1!=3 ) {
        for ( idim=0; idim<3; idim++ ) {
          e1[idim] = x2[idim]-x0[idim];
          e2[idim] = x3[idim]-x0[idim];
        }
        nrm[0] = e1[1]*e2[2] - e1[2]*e2[1];
        nrm[1] = e1[2]*e2[0] - e1[0]*e2[2];
        nrm[2] = e1[0]*e2[1] - e1[1]*e2[0];
        iface_measure += 0.5 * sqrt( nrm[0]*nrm[0] + nrm[1]*nrm[1] +
          nrm[2]*nrm[2] );
      }
    }
  }
  double *du_ip = get_new_dbl( ns1 * MDIM );
  array_set( du_ip, 0., ns1 * MDIM );
  if ( name==-BAR2 ) {
    for ( idim=0; idim<ndim; idim++ ) {
      double v1 = new_dof[0*nuknwn+vel_indx+idim*nder];
      double v2 = new_dof[1*nuknwn+vel_indx+idim*nder];
      du_ip[0*MDIM+idim] = ( v2 - v1 ) * dtime;
    }
  }
  else {
    for ( inol=0; inol<ns1; inol++ ) {
      for ( idim=0; idim<ndim; idim++ ) {
        double v1 = new_dof[inol*nuknwn+vel_indx+idim*nder];
        double v2 = new_dof[(inol+ns1)*nuknwn+vel_indx+idim*nder];
        du_ip[inol*MDIM+idim] = ( v2 - v1 ) * dtime;
      }
    }
  }
  // mean over the pairs (the value the whole element sees for the
  // records: for uniform loading all pairs carry the same du)
  ns1 = nnol/2;
  du[0] = 0.; du[1] = 0.; du[2] = 0.;
  for ( inol=0; inol<ns1; inol++ )
    for ( idim=0; idim<ndim; idim++ )
      du[idim] += w_ip[inol] * du_ip[inol*MDIM+idim];
  du_norm  = array_inproduct( du, normal, ndim );
  du_tang  = array_inproduct( du, tangent, ndim );
  du_tang2 = ( ndim==3 ) ? array_inproduct( du, tangent2, ndim ) : 0.;

  // accumulated normal strain per integration point (history arrays).
  // Sign convention: compression gives a NEGATIVE normal strain (the
  // normal is oriented so the compressed interface strain stays negative
  // - verified against the Professional .dbs of interface1/8/13/14/15).
  // The history is read BEFORE adding the current step's increment so that
  // gap / tension / Mohr-Coulomb decisions see the accumulated total.
  // (2026-08-30: the Professional stores one value per integration point
  //  - the .dbs of interface1 shows intpnt_strain with 3 entries for the
  //  quad6 - so the histories are arrays of ns1 now.)
  double *strain_normal_ip = get_new_dbl( ns1 );
  double *strain_eff_ip    = get_new_dbl( ns1 );
  double *strain_normal_old_ip = get_new_dbl( ns1 ); // pre-step value
  array_set( strain_normal_ip, 0., ns1 );
  array_set( strain_eff_ip, 0., ns1 );
  array_set( strain_normal_old_ip, 0., ns1 );
  {
    long int ln = ns1;
    double *hist = get_new_dbl( ns1 );
    array_set( hist, 0., ns1 );
    if ( db( ELEMENT_INTERFACE_STRAIN_NORMAL, element, idum, hist,
        ln, VERSION_NORMAL, GET_IF_EXISTS ) && ln>0 ) {
      for ( inol=0; inol<ns1; inol++ ) {
        strain_normal_ip[inol] = hist[inol];
        strain_normal_old_ip[inol] = hist[inol];
      }
    }
    delete[] hist;
    for ( inol=0; inol<ns1; inol++ )
      strain_normal_ip[inol] += array_inproduct( &du_ip[inol*MDIM], normal, ndim );
  }

  // accumulated contact normal STRESS per integration point (history
  // ELEMENT_INTERFACE_FORCE_NORM). CONVERGENCE (2026-09-04, interface2 +
  // control_reset_interface_strain): the stress the interface records and
  // assembles is NOT kn * (total accumulated strain) but kn times the
  // normal strain accumulated over the steps that end CLOSED (the free
  // travel of an open gap never builds stress; while open the stress is
  // exactly zero). control_reset_interface_strain (manual 6.355) resets
  // the strains but REMEMBERS this stress ("the interface stresses at
  // this moment of resetting will be remembered... the new interface
  // stresses are calculated from the interface stresses at this moment of
  // resetting plus stress due to additional deformation").
  double *force_norm_ip = get_new_dbl( ns1 );
  array_set( force_norm_ip, 0., ns1 );
  {
    long int ln = ns1;
    double *hist = get_new_dbl( ns1 );
    array_set( hist, 0., ns1 );
    if ( db( ELEMENT_INTERFACE_FORCE_NORM, element, idum, hist,
        ln, VERSION_NORMAL, GET_IF_EXISTS ) && ln>0 ) {
      for ( inol=0; inol<ns1; inol++ ) force_norm_ip[inol] = hist[inol];
    }
    delete[] hist;
  }

  // group_interface_materi_expansion_normal (manual Professional 6.630):
  // thermal strain expansion in interface thickness direction per unit
  // temperature; the temperature is the average of both sides. The
  // MECHANICAL strain seen by gap / tension / Mohr-Coulomb is the
  // accumulated strain minus the thermal expansion; the stored history
  // stays purely mechanical (no thermal re-counting per step). Only
  // meaningful with condif_temperature. The thermal strain INCREMENT of
  // this step (alpha * dT) is kept in thermal_inc_step: it belongs to the
  // stress accumulation (the Professional stress = kn * strain_eff, so a
  // heated constrained interface carries kn*(-alpha*T) - expans3 of the
  // corpus: interface sigma_n = -1 = 1*(-1) with alpha=1, T=1, u=0).
  double thermal_inc_step = 0.;
  for ( inol=0; inol<ns1; inol++ ) strain_eff_ip[inol] = strain_normal_ip[inol];
  if ( condif_temperature ) {
    double alpha_n = 0., t_side1 = 0., t_side2 = 0.;
    if ( db( GROUP_INTERFACE_MATERI_EXPANSION_NORMAL, element_group, idum,
        &alpha_n, ldum, VERSION_NORMAL, GET_IF_EXISTS ) ) {
      double t1_old = 0., t2_old = 0.;
      for ( inol=0; inol<ns1; inol++ ) {
        t_side1 += new_dof[inol*nuknwn+temp_indx];
        t_side2 += new_dof[(inol+ns1)*nuknwn+temp_indx];
        t1_old   += old_dof[inol*nuknwn+temp_indx];
        t2_old   += old_dof[(inol+ns1)*nuknwn+temp_indx];
      }
      t_side1 /= ns1; t_side2 /= ns1; t1_old /= ns1; t2_old /= ns1;
      thermal_inc_step = alpha_n *
        ( 0.5*(t_side1+t_side2) - 0.5*(t1_old+t2_old) );
      for ( inol=0; inol<ns1; inol++ ) {
        // total thermal expansion (for gap / tension / Mohr-Coulomb state)
        strain_eff_ip[inol] = strain_normal_ip[inol] -
          alpha_n * 0.5*(t_side1+t_side2);
        // incremental thermal expansion acts as a pseudo-load this step.
        // The contraction is along the interface NORMAL (the expansion is
        // in the interface thickness direction, manual 6.629): subtracting
        // the increment from the x component only (pre-2026-09-04) gave a
        // spurious tangential slip of -alpha*dT*normal_x on inclined
        // interfaces (expans3: interface at 45 degrees got a shear stress
        // -0.707 = kt * alpha*dT/sqrt(2) instead of 0).
        for ( idim=0; idim<ndim; idim++ )
          du_ip[inol*MDIM+idim] -= thermal_inc_step * normal[idim];
      }
    }
  }

  // accumulated total tangential force per IP (history arrays).
  // Default 0: without the Mohr-Coulomb record the interface stays purely
  // elastic (Fase 1 behavior). 3D: two tangential components.
  double *f_t_old_ip = get_new_dbl( ns1 );
  double *f_t_ip     = get_new_dbl( ns1 );
  double *f_t2_old_ip = get_new_dbl( ns1 );
  double *f_t2_ip     = get_new_dbl( ns1 );
  double *f_t_el_ip   = get_new_dbl( ns1 );
  double *f_t2_el_ip  = get_new_dbl( ns1 );
  double *rhs_shear_ip  = get_new_dbl( ns1 );
  double *rhs_shear2_ip = get_new_dbl( ns1 );
  array_set( f_t_old_ip, 0., ns1 );
  array_set( f_t_ip, 0., ns1 );
  array_set( f_t2_old_ip, 0., ns1 );
  array_set( f_t2_ip, 0., ns1 );
  {
    long int ln = ns1;
    double *hist = get_new_dbl( ns1 );
    array_set( hist, 0., ns1 );
    if ( db( ELEMENT_INTERFACE_FORCE_TANG, element, idum, hist,
        ln, VERSION_NORMAL, GET_IF_EXISTS ) && ln>0 ) {
      for ( inol=0; inol<ns1; inol++ ) f_t_old_ip[inol] = hist[inol];
    }
    delete[] hist;
    if ( ndim==3 ) {
      hist = get_new_dbl( ns1 );
      array_set( hist, 0., ns1 );
      if ( db( ELEMENT_INTERFACE_FORCE_TANG2, element, idum, hist,
          ln, VERSION_NORMAL, GET_IF_EXISTS ) && ln>0 ) {
        for ( inol=0; inol<ns1; inol++ ) f_t2_old_ip[inol] = hist[inol];
      }
      delete[] hist;
    }
  }

  // gap (CONVERGENCE 2026-09-04, manual Professional 6.625): the
  // interface is CLOSED when the accumulated normal strain <= gap and
  // OPEN when strain_normal > gap ("Only when the sides displacements are
  // such that the normal strain becomes lower then the specified gap
  // value the interface will be closed and start to generate stresses").
  // An opened interface "does not have stresses" (6.628) and keeps only
  // the residual stiffness for the matrix. A physical gap is a NEGATIVE
  // gap value (open until compression exceeds |gap|); the DEFAULT without
  // the record is +1e20 = always closed (allows tension stresses, 6.625:
  // "If you want to allow always tension stresses in an interface set gap
  // to, by example, 1.e20"). NOTE: the OLD code (before 2026-09-04)
  // inverted the condition (open when strain <= gap) with default -1e20:
  // it worked for the tests without a gap record but left interfaces with
  // an explicit positive gap (patch1: gap 1.e20) ALWAYS OPEN and let the
  // closed phase of a physical gap (interface2) carry only the residual
  // stress.
  if ( !db( GROUP_INTERFACE_GAP, element_group, idum, &gap, ldum,
      VERSION_NORMAL, GET_IF_EXISTS ) )
    gap = 1.e20;

  // per-IP constitutive state
  // Mohr-Coulomb (FIX 1, RF-1): friction limit on the TOTAL tangential
  // force. The law is active by the PRESENCE of the record (D2): phi=0,c=0
  // gives max_fric=0 -> free sliding; without the record the interface
  // stays purely elastic (Fase 1 behavior). The manual 6.631: the maximum
  // friction force is c + Fn*tan(phi) where Fn = kn*strain_eff is the
  // normal FORCE (negative under compression).
  mc_active = db( GROUP_INTERFACE_MATERI_PLASTI_MOHR_COUL_DIRECT, element_group,
    idum, ddum3, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  phi = ddum3[0]; c = ddum3[1]; phi_flow = ddum3[2];
  double *stiff_normal_ip = get_new_dbl( ns1 );
  double *stress_normal_ip = get_new_dbl( ns1 );
  double *stress_shear_ip  = get_new_dbl( ns1 );
  double *stress_shear2_ip = get_new_dbl( ns1 );
  array_set( f_t_el_ip, 0., ns1 );
  array_set( f_t2_el_ip, 0., ns1 );
  array_set( rhs_shear_ip, 0., ns1 );
  array_set( rhs_shear2_ip, 0., ns1 );
  double *stiff_tang_ip    = get_new_dbl( ns1 );
  double *stiff_tang2_ip   = get_new_dbl( ns1 );
  long int *plastified_ip  = get_new_int( ns1 );
  for ( inol=0; inol<ns1; inol++ ) {
    double du_tang_i  = array_inproduct( &du_ip[inol*MDIM], tangent, ndim );
    double du_tang2_i = ( ndim==3 ) ?
      array_inproduct( &du_ip[inol*MDIM], tangent2, ndim ) : 0.;
    double stiff_normal_i = kn;
    // state after this step's increment: CLOSED when the accumulated
    // normal strain <= gap, OPEN otherwise (manual 6.625; the strain_eff
    // history was accumulated with this step's du above)
    long int iface_open_i = ( strain_eff_ip[inol] > gap );
    if ( iface_open_i ) stiff_normal_i = kn * residual_factor;
    // tension limit (FIX 2, RF-2): opens in traction when the TOTAL
    // accumulated normal force |Fn_total| exceeds the limit, and only if
    // it was still closed. An interface opened by the tension limit has
    // no normal stress either (same handling as the gap-open state).
    double fn_total_i = kn * strain_eff_ip[inol];
    if ( tension_limit>0. && strain_normal_ip[inol]<0. &&
        fabs(fn_total_i)>tension_limit && !iface_open_i ) {
      stiff_normal_i = kn * residual_factor;
      iface_open_i = 1;
    }
    stiff_normal_ip[inol] = stiff_normal_i;

    // cumulative Mohr-Coulomb (FIX 1, RF-1): friction limit on the TOTAL
    // tangential force. CONVERGENCE (2026-08-30): F_t = kt*du_tang (not
    // kt*2*du_tang), and the plastic return is IMPLICIT with dilatancy
    // (see the block comments above; verified EXACT against interface15).
    // CONVERGENCE (interface_patch): the Professional's trial is the
    // ACCUMULATED tangential strain times the stiffness (kt*gamma_total),
    // NOT the last-step increment f_old + kt*du_paso. With the incremental
    // trial the plastic slip exhausts the relative displacement and the
    // friction force stalls below the cohesion (f_t=-5.55 vs target -10).
    // The stored history ELEMENT_INTERFACE_FORCE_TANG holds the ELASTIC
    // trial (grows every step, even when plastic: it is kt*gamma_total
    // with gamma_total = sum(du_tang) over all steps) and the resulting
    // force is f_t = clamp(trial) to the yield limit.
    // CONVERGENCE (2026-09-03, mohr_coul_direct3/4): the interface
    // assembles the FULL accumulated forces into the rhs (spring.cc
    // pattern, like every other element), NOT the step increment. With
    // an incremental rhs the equilibrium of a multi-step run only sees
    // the last increment: the Professional's node_rhside = w*kn*eps_acc
    // (direct3: -100 per node = 0.5*(-200)) while the old code reported
    // the last-step increment (0.5*(-100) = -50). The incremental rhs
    // also made the interface creep one increment per step under a
    // constant load (interface9: 10 equal steps vs the Professional's
    // single-step static equilibrium).
    double trial_el_i  = f_t_old_ip[inol]  + kt1 * du_tang_i;
    double trial_el2_i = f_t2_old_ip[inol] + kt2 * du_tang2_i;
    long int plast_i = 0;
    double f_t_clamped_i = trial_el_i;
    double f_t2_clamped_i = trial_el2_i;
    if ( mc_active ) {
      double strain_eff_mc_i = strain_eff_ip[inol];
      double trial_mag = sqrt( trial_el_i*trial_el_i + trial_el2_i*trial_el2_i );
      // Normal force of the yield limit: the Professional uses the normal
      // force POSITIVE under compression (direct4 plateau 1.20271 =
      // c + |kn*eps|*tan(phi) with c=1, eps=-1e-6, phi=0.2). The stored
      // accumulated strain is NEGATIVE under compression, hence the
      // minus sign: max_fric = |c - kn*eps*tan(phi)| = |c + Fn*tan(phi)|.
      double max_fric_abs = fabs( c - kn * strain_eff_mc_i * tan( phi ) );
      if ( trial_mag > max_fric_abs && trial_mag>0. ) {
        double dgamma = ( trial_mag - max_fric_abs ) /
          ( kt1 + kn * tan( phi ) * tan( phi_flow ) );
        if ( dgamma<0. ) dgamma = 0.;
        double scale = ( trial_mag - kt1 * dgamma ) / trial_mag;
        f_t_clamped_i  = trial_el_i  * scale;
        f_t2_clamped_i = trial_el2_i * scale;
        plast_i = 1;
        // dilatancy: plastic slip OPENS the interface (reduces the
        // accumulated compressive strain). The return multiplier dgamma
        // is the TOTAL plastic slip (the accumulated-trial return maps
        // the whole history back to the surface); the opening of THIS
        // step must use only the INCREMENTAL plastic slip, otherwise the
        // slip already accumulated in the past is re-counted every step
        // (direct3 sigma_n = -200 = kn*(-1e-6 - 1e-6*tan(pi/4)) for a
        // total slip of 1e-6, NOT -250). The past plastic slip follows
        // from the stored elastic trial: gamma_pl_old =
        // (|trial_old| - |clamp(trial_old)|)/kt.
        double f_el_mag_old = sqrt( f_t_old_ip[inol]*f_t_old_ip[inol] +
                                    f_t2_old_ip[inol]*f_t2_old_ip[inol] );
        double f_old_surf_mag = f_el_mag_old;
        if ( f_old_surf_mag > max_fric_abs ) f_old_surf_mag = max_fric_abs;
        double gamma_pl_old = ( kt1>0. ) ? ( f_el_mag_old - f_old_surf_mag )/kt1 : 0.;
        double dgamma_inc = dgamma - gamma_pl_old;
        if ( dgamma_inc<0. ) dgamma_inc = 0.;
        strain_normal_ip[inol] += -dgamma_inc * tan( phi_flow );
      }
    }
    // f_t: the clamped force of the current state (record + assembly)
    f_t_ip[inol]  = f_t_clamped_i;
    f_t2_ip[inol] = f_t2_clamped_i;
    // the stored history is the ELASTIC trial (keeps accumulating)
    // so the next step's trial continues from kt*gamma_total
    f_t_el_ip[inol]  = trial_el_i;
    f_t2_el_ip[inol] = trial_el2_i;
    plastified_ip[inol] = plast_i;
    stiff_tang_ip[inol]  = ( mc_active && plast_i ) ? 0. : kt1;
    stiff_tang2_ip[inol] = ( mc_active && plast_i ) ? 0. : kt2;
    // stress,normal (spring.cc pattern, CONVERGENCE 2026-09-04): the
    // accumulated contact stress of the CLOSED phase only (history
    // ELEMENT_INTERFACE_FORCE_NORM, read before the step). Each closed
    // step adds kn * (normal strain increment of this step - the relative
    // displacement du plus the dilatancy opening added above), so after n
    // closed steps the stress = kn * (sum of their strain increments).
    // Steps ending OPEN zero the history: an opened interface does not
    // have stresses (manual 6.628) and a later re-closure rebuilds the
    // stress from the penetration of the closing step (interface2 of the
    // corpus: gap 0.1 closes at step 100 of 200 and the final stress is
    // -101 = kn * (-101 * 1e-3), NOT kn * (-0.2) which would count the
    // 0.1 free gap travel; verified step-by-step against the Professional
    // per-step prints). Without a gap record the interface is always
    // closed and the history equals kn * strain,normal (the behavior of
    // all non-gap corpus tests, bit-identical up to FP round-off).
    // control_reset_interface_strain (manual 6.355) keeps this history
    // untouched: "the interface stresses at this moment of resetting will
    // be remembered... the new interface stresses are calculated from the
    // interface stresses at this moment of resetting plus stress due to
    // additional deformation" (interface10 of the corpus). stress_normal_ip
    // is set AFTER the plastic block so the record/assembly carry the
    // post-dilatancy value.
    double force_norm_new_i = force_norm_ip[inol];
    if ( iface_open_i ) {
      force_norm_new_i = 0.;
    }
    else {
      // the closed-phase stress increment: kn times the effective normal
      // strain increment of this step (relative displacement du + plastic
      // dilatancy opening - thermal contraction alpha*dT, expans3)
      force_norm_new_i += kn *
        ( strain_normal_ip[inol] - strain_normal_old_ip[inol] -
          thermal_inc_step );
    }
    force_norm_ip[inol] = force_norm_new_i;
    stress_normal_ip[inol] = force_norm_new_i;
    // stress,shear = the accumulated tangential force clamped to the
    // yield limit (WITH Mohr-Coulomb) or the accumulated elastic trial
    // kt*gamma_total (WITHOUT: elastic, single- and multi-step alike -
    // interface9 -0.159 = kt*du_total, interface14 5e3 = kt*1).
    stress_shear_ip[inol]  = f_t_ip[inol];
    stress_shear2_ip[inol] = f_t2_ip[inol];
    // the assembled rhs carries the FULL current forces (see the block
    // comment above); the names rhs_shear_* are kept for the assembly.
    rhs_shear_ip[inol]  = f_t_ip[inol];
    rhs_shear2_ip[inol] = f_t2_ip[inol];
  }

  // mean over the pairs (for the swit debug and legacy scalar view)
  du_norm = 0.; du_tang = 0.; du_tang2 = 0.; strain_normal = 0.;
  stress_normal = 0.; stress_shear = 0.; stress_shear2 = 0.; f_t = 0.;
  f_t_old = 0.; f_t2 = 0.; f_t2_old = 0.;
  for ( inol=0; inol<ns1; inol++ ) {
    du_norm  += w_ip[inol] * array_inproduct( &du_ip[inol*MDIM], normal, ndim );
    du_tang  += w_ip[inol] * array_inproduct( &du_ip[inol*MDIM], tangent, ndim );
    du_tang2 += w_ip[inol] * ( ( ndim==3 ) ?
      array_inproduct( &du_ip[inol*MDIM], tangent2, ndim ) : 0. );
    strain_normal += w_ip[inol] * strain_normal_ip[inol];
    stress_normal += w_ip[inol] * stress_normal_ip[inol];
    stress_shear  += w_ip[inol] * stress_shear_ip[inol];
    stress_shear2 += w_ip[inol] * stress_shear2_ip[inol];
    f_t  += w_ip[inol] * f_t_ip[inol];
    f_t_old += w_ip[inol] * f_t_old_ip[inol];
    f_t2 += w_ip[inol] * f_t2_ip[inol];
    f_t2_old += w_ip[inol] * f_t2_old_ip[inol];
  }

  if ( swit ) {
    pri( "du_norm", du_norm );
    pri( "du_tang", du_tang );
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
  // CONVERGENCE (2026-08-30): each facing pair (i, i+ns1) is assembled as
  // an independent spring with its own Lobatto weight w_ip[i] and its own
  // per-IP stress/stiffness. The old code assembled every node against
  // every node with the full stiffness, which produced a rank-1 matrix
  // per side (SINGULAR for ns1>1) and applied ns1 times the force.
  // The Professional nodal force is weight_i * sigma_i (interface1:
  // loads -1,-4,-1 = (1/6,4/6,1/6)*(-6)).
  for ( inol=0; inol<ns1; inol++ ) {
    double w_i = w_ip[inol] * iface_measure;
    double s_normal_i = stress_normal_ip[inol];
    double s_shear_i  = rhs_shear_ip[inol];
    double s_shear2_i = rhs_shear2_ip[inol];
    double stiff_n_i = stiff_normal_ip[inol];
    double stiff_t_i = stiff_tang_ip[inol];
    double stiff_t2_i= stiff_tang2_ip[inol];
    // axisymmetric: the interface is a ring; the pair force/stiffness
    // carry the circumference weight 2*pi*r (radius = radial coord of
    // the side-1 node, pattern area.cc)
    double w_ax = w_i;
    if ( axisymmetric==-YES ) {
      double *cr = NULL;
      if ( memory==-TOTAL_LINEAR )
        cr = db_dbl( NODE_START_REFINED, nodes[inol], VERSION_NORMAL );
      else
        cr = &coord[inol*ndim];
      w_ax *= 2. * PIRAD * cr[0];
    }
    for ( idim=0; idim<ndim; idim++ ) {
      double dirn = normal[idim], dirt = tangent[idim], dirt2 = tangent2[idim];
      // side 1 node (negative sign), side 2 node (positive sign)
      for ( jnol=0; jnol<2; jnol++ ) {
        long int inod = inol + ( jnol ? ns1 : 0 );
        double sign = ( jnol ) ? +1. : -1.;
        indx = inod*npuknwn + (vel_indx+idim*nder)/nder;
        tmp = -sign * w_ax * ( s_normal_i*dirn + s_shear_i*dirt +
          s_shear2_i*dirt2 );
        // damping: -d*(v_side2 - v_side1) on the relative velocity
        // (v = du/dtime), assembled like the stiffness but WITHOUT the
        // dtime factor (the damping force is proportional to velocity)
        if ( damping_iface>0. ) {
          double v_rel = du_ip[inol*MDIM+idim] / dtime;
          tmp += -sign * w_ax * damping_iface * v_rel;
        }
        element_rhside[indx] += tmp;
        for ( jdim=0; jdim<ndim; jdim++ ) {
          double jdirn = normal[jdim], jdirt = tangent[jdim],
            jdirt2 = tangent2[jdim];
          double kkk = sign * w_ax * ( stiff_n_i*dirn*jdirn +
            stiff_t_i*dirt*jdirt + stiff_t2_i*dirt2*jdirt2 );
          double kkk_damp = 0.;
          if ( damping_iface>0. && jdim==idim )
            kkk_damp = sign * w_ax * damping_iface;
          // jnode: same pair's side-1 or side-2 node
          for ( knol=0; knol<2; knol++ ) {
            long int jnod = inol + ( knol ? ns1 : 0 );
            double jsign = ( knol ) ? +1. : -1.;
            long int jndx = inod*npuknwn*nnol*npuknwn +
              ((vel_indx+idim*nder)/nder)*nnol*npuknwn +
              jnod*npuknwn + (vel_indx+jdim*nder)/nder;
            element_matrix[jndx] += kkk * jsign * dtime;
            element_matrix[jndx] += kkk_damp * jsign;
            if ( jnod==inod && jdim==idim )
              element_lhside[inod*npuknwn+(vel_indx+idim*nder)/nder] +=
                kkk * jsign * dtime + kkk_damp * jsign;
          }
        }
      }
    }
  }

  // store accumulated histories per IP (used by gap / tension /
  // Mohr-Coulomb and the next step's cumulative trial)
  {
    long int ln = ns1;
    db( ELEMENT_INTERFACE_STRAIN_NORMAL, element, idum, strain_normal_ip,
      ln, VERSION_NEW, PUT );
    // the accumulated contact normal stress of the closed phase (see the
    // per-IP block); zeroed while open, remembered by
    // control_reset_interface_strain (data.cc)
    db( ELEMENT_INTERFACE_FORCE_NORM, element, idum, force_norm_ip,
      ln, VERSION_NEW, PUT );
    // the history holds the ELASTIC trial (kt*gamma_total) so the
    // accumulated tangential strain keeps growing after plastic slip
    // (the clamped force is rebuilt from it every step)
    db( ELEMENT_INTERFACE_FORCE_TANG, element, idum, f_t_el_ip,
      ln, VERSION_NEW, PUT );
    if ( ndim==3 )
      db( ELEMENT_INTERFACE_FORCE_TANG2, element, idum, f_t2_el_ip,
        ln, VERSION_NEW, PUT );
  }

  // output records (Professional compatibility, verified against .dbs of
  // interface1/13/14/15): the Professional fills element_interface_intpnt_stress,
  // _intpnt_strain, _stress_average, _strain_average and the per-intpnt
  // statuses at the end of the calculation. Semantics (deduced from the
  // .dbs + manual 6.431-6.436):
  //   - n_intpnt = nnol/2 (nodes per side): bar2=1, quad4=2, quad6=3,
  //     prism6=3, hex8=4. Each intpnt carries (normal, shear_first,
  //     [shear_second 3D]).
  //   - strain,normal = du_norm of the CURRENT step (incremental, not the
  //     accumulated history: interface7 with 2 steps of vely=-1 targets
  //     -1e10 = kn*(-1) = kn * du of the last step);
  //   - strain,shear = du_tang/2 (interface13: du=(2,-1), tangent=(2,1)/|..|
  //     -> du_tang=3/sqrt(5)=1.34164, gamma=0.67082 EXACT; the manual:
  //     shear,gamma = 2*strain,shear so strain,shear = du_tang/2);
  //   - stress,normal = the accumulated contact stress (FORCE_NORM
  //     history): kn times the strain accumulated over CLOSED steps only,
  //     dilatancy opening INCLUDED (interface15 sigma_n = -1730.48 =
  //     1e4 * (-0.17305)); zero while the interface is open (interface2:
  //     the gap closes at step 100/200 and the stress grows from there to
  //     -101 = kn * (-0.101));
  //   - stress,shear = F_t - F_t,old (the incremental tangential force,
  //     == kt*du_tang elastic, == clamped increment plastic);
  //   - stress_average/strain_average = mean over the intpnts.
  // The intpnts of a uniform interface coincide with the nodes (the
  // interpolation is linear/quadratic; for uniform loading all intpnts
  // carry the same value).
  {
    long int n_intpnt = nnol/2;
    long int nval = ( ndim==3 ) ? 3 : 2;
    long int i_intpnt = 0, iv = 0;
    double *rec_stress = get_new_dbl( n_intpnt * nval );
    double *rec_strain = get_new_dbl( n_intpnt * nval );
    double *rec_status = (double*)get_new_int( n_intpnt );
    double *avg_stress = get_new_dbl( nval );
    double *avg_strain = get_new_dbl( nval );
    array_set( rec_stress, 0., n_intpnt * nval );
    array_set( rec_strain, 0., n_intpnt * nval );
    array_set( avg_stress, 0., nval );
    array_set( avg_strain, 0., nval );
    for ( i_intpnt=0; i_intpnt<n_intpnt; i_intpnt++ ) {
      // the current-step values per integration point
      rec_stress[i_intpnt*nval+0] = stress_normal_ip[i_intpnt];
      rec_stress[i_intpnt*nval+1] = stress_shear_ip[i_intpnt];
      rec_strain[i_intpnt*nval+0] =
        array_inproduct( &du_ip[i_intpnt*MDIM], normal, ndim );
      rec_strain[i_intpnt*nval+1] =
        array_inproduct( &du_ip[i_intpnt*MDIM], tangent, ndim ) / 2.;
      if ( ndim==3 ) {
        rec_stress[i_intpnt*nval+2] = stress_shear2_ip[i_intpnt];
        rec_strain[i_intpnt*nval+2] =
          array_inproduct( &du_ip[i_intpnt*MDIM], tangent2, ndim ) / 2.;
      }
      rec_status[i_intpnt] =
        ( strain_eff_ip[i_intpnt] <= gap ) ? CLOSED : OPENED;
      for ( iv=0; iv<nval; iv++ ) {
        avg_stress[iv] += rec_stress[i_intpnt*nval+iv] / n_intpnt;
        avg_strain[iv] += rec_strain[i_intpnt*nval+iv] / n_intpnt;
      }
    }
    ldum = n_intpnt * nval;
    // VERSION_NORMAL (t) is what print_database and the target checker read;
    // VERSION_NEW (t+dt) keeps the record current for the next step's loop.
    db( ELEMENT_INTERFACE_INTPNT_STRESS, element, idum, rec_stress,
      ldum, VERSION_NORMAL, PUT );
    db( ELEMENT_INTERFACE_INTPNT_STRESS, element, idum, rec_stress,
      ldum, VERSION_NEW, PUT );
    db( ELEMENT_INTERFACE_INTPNT_STRAIN, element, idum, rec_strain,
      ldum, VERSION_NORMAL, PUT );
    db( ELEMENT_INTERFACE_INTPNT_STRAIN, element, idum, rec_strain,
      ldum, VERSION_NEW, PUT );
    ldum = n_intpnt;   // tension_status: one value per integration point
    db( ELEMENT_INTERFACE_INTPNT_MATERI_TENSION_STATUS, element, idum,
      rec_status, ldum, VERSION_NORMAL, PUT );
    db( ELEMENT_INTERFACE_INTPNT_MATERI_TENSION_STATUS, element, idum,
      rec_status, ldum, VERSION_NEW, PUT );
    ldum = nval;
    db( ELEMENT_INTERFACE_STRESS_AVERAGE, element, idum, avg_stress,
      ldum, VERSION_NORMAL, PUT );
    db( ELEMENT_INTERFACE_STRESS_AVERAGE, element, idum, avg_stress,
      ldum, VERSION_NEW, PUT );
    db( ELEMENT_INTERFACE_STRAIN_AVERAGE, element, idum, avg_strain,
      ldum, VERSION_NORMAL, PUT );
    db( ELEMENT_INTERFACE_STRAIN_AVERAGE, element, idum, avg_strain,
      ldum, VERSION_NEW, PUT );
    free( rec_stress ); free( rec_strain ); free( rec_status );
    free( avg_stress ); free( avg_strain );
  }

  // release per-IP work arrays
  delete[] w_ip; delete[] du_ip;
  delete[] strain_normal_ip; delete[] strain_eff_ip;
  delete[] strain_normal_old_ip; delete[] force_norm_ip;
  delete[] f_t_old_ip; delete[] f_t_ip;
  delete[] f_t2_old_ip; delete[] f_t2_ip;
  delete[] stiff_normal_ip; delete[] stress_normal_ip;
  delete[] stress_shear_ip; delete[] stress_shear2_ip;
  delete[] stiff_tang_ip; delete[] stiff_tang2_ip;
  delete[] plastified_ip;

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
//   -bar3 -> -quad6 (2D quadratic): the same split with 3 nodes per
//   side (the Professional auto-converts the quadratic interfaces;
//   interface_bar3_quad8.dat of the corpus states it textually). The
//   interface_element() routine already supports the -quad6 (3+3
//   nodes, Lobatto weights 1/6,4/6,1/6).
//
// See ProjectDocs/DESIGN-INTERFACES.md for the full algorithm.
void interface_convert( long int icontrol )

{
  long int element=0, max_element=0, max_element_c=0, i=0, j=0, jnod=0,
    name=0, length=0, ldum=0, swit=0, element_group=0,
    max_node=0, length_convert_groups=0, found=0, nconv=0,
    idum[1], *el=NULL, *convert_groups=NULL;
  double ddum[1], coord[MDIM], normal[MDIM], tangent[MDIM], shift=0.,
    *ca=NULL, *cb=NULL;

  swit = set_swit(-1,-1,"interface_convert");
  if ( swit ) pri( "In routine INTERFACE_CONVERT." );

  // control_mesh_convert may be declared with ANY control index (the
  // corpus uses e.g. control_mesh_convert 20 -yes with
  // control_timestep 30), so look it up over the whole control range
  // instead of only the current icontrol. The first active -yes record
  // found drives the conversion.
  {
    long int ic_max = 0, ic_found = -1;
    long int conv_sw = -NO;
    db_max_index( CONTROL_MESH_CONVERT, ic_max, VERSION_NORMAL, GET );
    for ( long int ic2=0; ic2<=ic_max; ic2++ ) {
      if ( db_active_index( CONTROL_MESH_CONVERT, ic2, VERSION_NORMAL ) ) {
        db( CONTROL_MESH_CONVERT, ic2, &conv_sw, ddum, ldum,
          VERSION_NORMAL, GET );
        if ( conv_sw==-YES ) { ic_found = ic2; break; }
      }
    }
    if ( ic_found<0 ) return;
    icontrol = ic_found;
  }

  el = get_new_int(MAXIMUM_NODE+1);
  db_highest_index( ELEMENT, max_element, VERSION_NORMAL );
  max_element_c = max_element;
  db_highest_index( NODE, max_node, VERSION_NORMAL );

  // element groups located on ONE side of the interfaces; neighbours in
  // these groups keep the original nodes, the others get the new ones
  convert_groups = get_new_int(DATA_ITEM_SIZE);
  db( CONTROL_MESH_CONVERT_ELEMENT_GROUP, icontrol, convert_groups,
    ddum, length_convert_groups, VERSION_NORMAL, GET_IF_EXISTS );
  nconv = 0;

  // PRE-PASS (2026-09-07, interface_tria3_prism6): per interface element
  // to convert, record on which side of its normal the element-record
  // side 1 sits. Professional layout of the converted family (verified
  // against its .dbs of interface_tria3_prism6 / interface_quad4_hex8 /
  // interface_bar2_quad4 and variants with swapped element numbers):
  // side 1 of the record is the block of the LOWEST-numbered volume
  // element that shares the whole interface face. side1_minus[e] = 1 when
  // that block lies on the -normal side of the element, 0 when it lies on
  // the +normal side (a clockwise triangle of a triangulated quad keeps
  // its exclusive corner on the opposite block than the counter-clockwise
  // one: interface_tria3_prism6 - bottom hex8 1 2 3 4 5 6 7 22 with the
  // COPY of node 8, top hex8 19 20 21 8 ... with the ORIGINAL 8).
  // Only the 3D surface conversions need it (the 2D orientation is the
  // numbering heuristic of interface_element(); there the record keeps
  // side 1 = the originals = the -normal side, like the old code).
  long int *side1_minus = get_new_int( max_element+1 );
  long int nconv_total = 0;
  long int *nel_pre = get_new_int( MAXIMUM_NODE+1 );
  for ( long int e2=0; e2<=max_element; e2++ ) {
    side1_minus[e2] = 1;
    if ( !db_active_index( ELEMENT, e2, VERSION_NORMAL ) ) continue;
    db( ELEMENT, e2, el, ddum, length, VERSION_NORMAL, GET );
    name = el[0];
    long int gr2 = 0;
    db( ELEMENT_GROUP, e2, &gr2, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    if ( name!=-BAR2 && name!=-BAR3 && name!=-TRIA3 && name!=-QUAD4 &&
         !( name==-QUAD8 && ndim==3 ) ) continue;
    if ( !db_active_index( GROUP_INTERFACE, gr2, VERSION_NORMAL ) ) continue;
    nconv_total++;
    if ( !( ndim==3 && ( name==-TRIA3 || name==-QUAD4 ) ) ) continue;
    long int ns1p = ( name==-TRIA3 ) ? 3 : 4;
    // element normal (same formulas as the conversion below)
    double np[MDIM], e1p[MDIM], e2p[MDIM];
    double *cap = db_dbl( NODE, el[1], VERSION_NORMAL );
    double *cbp = db_dbl( NODE, el[2], VERSION_NORMAL );
    double *ccp = db_dbl( NODE, el[3], VERSION_NORMAL );
    for ( i=0; i<3; i++ ) { e1p[i] = cbp[i]-cap[i]; e2p[i] = ccp[i]-cap[i]; }
    np[0] = e1p[1]*e2p[2] - e1p[2]*e2p[1];
    np[1] = e1p[2]*e2p[0] - e1p[0]*e2p[2];
    np[2] = e1p[0]*e2p[1] - e1p[1]*e2p[0];
    array_normalize( np, 3 );
    // face centroid of the interface element
    double fcp[MDIM]; array_set( fcp, 0., MDIM );
    for ( j=0; j<ns1p; j++ ) {
      double *cnp = db_dbl( NODE, el[1+j], VERSION_NORMAL );
      for ( i=0; i<3; i++ ) fcp[i] += cnp[i]/(double)ns1p;
    }
    long int min_minus = -1, min_plus = -1;
    for ( long int iel2=0; iel2<=max_element; iel2++ ) {
      if ( iel2==e2 ) continue;
      if ( !db_active_index( ELEMENT, iel2, VERSION_NORMAL ) ) continue;
      long int lnp = 0, grp = 0;
      db( ELEMENT_GROUP, iel2, &grp, ddum, ldum, VERSION_NORMAL,
        GET_IF_EXISTS );
      if ( db_active_index( GROUP_INTERFACE, grp, VERSION_NORMAL ) ) continue;
      db( ELEMENT, iel2, nel_pre, ddum, lnp, VERSION_NORMAL, GET );
      long int shared = 0;
      for ( long int jj=0; jj<ns1p; jj++ ) {
        long int s2p = el[1+jj];
        for ( long int k2=1; k2<lnp; k2++ )
          if ( nel_pre[k2]==s2p ) { shared++; break; }
      }
      if ( shared!=ns1p ) continue;   // not sharing the whole face
      double cx = 0., cy = 0., cz = 0.;
      long int nn2 = ( lnp>1 ) ? lnp-1 : 1;
      for ( long int k2=1; k2<lnp; k2++ ) {
        double *cn2 = db_dbl( NODE, nel_pre[k2], VERSION_NORMAL );
        cx += cn2[0]; cy += cn2[1]; cz += cn2[2];
      }
      cx /= nn2; cy /= nn2; cz /= nn2;
      double dot = (cx-fcp[0])*np[0] + (cy-fcp[1])*np[1] +
        (cz-fcp[2])*np[2];
      if ( dot<0. && ( min_minus<0 || iel2<min_minus ) ) min_minus = iel2;
      if ( dot>0. && ( min_plus<0 || iel2<min_plus ) ) min_plus = iel2;
    }
    if ( min_minus<0 && min_plus>=0 ) side1_minus[e2] = 0;
    else if ( min_minus>=0 && min_plus>=0 && min_plus<min_minus )
      side1_minus[e2] = 0;
    // else: side 1 = the -normal side (default; also when only the
    // -normal side has volumes or none was found)
  }
  delete[] nel_pre;
  // per-node pair state of the conversion: every interface node is
  // duplicated at most once. node_s1/s2 hold the record side-1/side-2
  // node of the pair (original or duplicate) once the node is classified.
  long int node_buf = max_node + 10*nconv_total + 8;
  long int *node_done = get_new_int( node_buf+1 );
  long int *node_dup  = get_new_int( node_buf+1 );
  long int *node_s1   = get_new_int( node_buf+1 );
  long int *node_s2   = get_new_int( node_buf+1 );
  // get_new_int does NOT zero the memory (plain new[]): the pair state
  // must start all-zero (node 0 is never a valid duplicate target).
  for ( long int kz=0; kz<=node_buf; kz++ ) {
    node_done[kz] = 0; node_dup[kz] = 0;
    node_s1[kz] = 0;   node_s2[kz] = 0;
  }

  for ( element=0; element<=max_element; element++ ) {
    if ( !db_active_index( ELEMENT, element, VERSION_NORMAL ) ) continue;
    db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
    name = el[0];
    element_group = 0;
    db( ELEMENT_GROUP, element, &element_group, ddum, ldum,
      VERSION_NORMAL, GET_IF_EXISTS );
    if ( name!=-BAR2 && name!=-BAR3 && name!=-TRIA3 && name!=-QUAD4 &&
         !( name==-QUAD8 && ndim==3 ) )
      continue;
    if ( !db_active_index( GROUP_INTERFACE, element_group, VERSION_NORMAL ) )
      continue;

    // -quad8 (3D): quadratic FACIAL interface element (the
    // interface_quad8_hex20 family: a -quad8 face between two -hex20 /
    // -hex27 volumes, converted to the -hex18 interface, 9+9 nodes).
    // The quad8 record carries only 8 nodes (corners BL,BR,TL,TR then
    // mid-edge BM,LM,RM,TM - no centre): the 9th side-1 node of the
    // hex18 is the FACE CENTRE, the node the mesh_convert_hex20
    // auto-conversion created (deduplicated) on the shared face of the
    // hex27 volumes - looked up by coordinates, created if the mesh has
    // no hex20. The side-1 ordering is rewritten to the GNU -quad9
    // tensor order [BL,BM,BR,LM,C,RM,TL,TM,TR] (the Professional .dbs
    // of the corpus interface_quad8_hex20 shows exactly this layout).
    if ( name==-QUAD8 ) {
      if ( length!=1+8 ) db_error( ELEMENT, element );
      long int q8[1+8], ns1q = 9;
      for ( j=0; j<=8; j++ ) q8[j] = el[j];
      // face centre = average of the 4 corners (q8[1..4])
      double centre[MDIM];
      array_set( centre, 0., MDIM );
      for ( j=0; j<4; j++ ) {
        double *cq = db_dbl( NODE, q8[1+j], VERSION_NORMAL );
        for ( i=0; i<3; i++ ) centre[i] += cq[i];
      }
      for ( i=0; i<3; i++ ) centre[i] /= 4.;
      long int cn = -1;
      for ( jnod=0; jnod<=max_node; jnod++ ) {
        if ( !db_active_index( NODE, jnod, VERSION_NORMAL ) ) continue;
        double *cq = db_dbl( NODE, jnod, VERSION_NORMAL );
        long int ok = 1;
        for ( i=0; i<3 && ok; i++ )
          if ( fabs( cq[i]-centre[i] )>1.e-10 ) ok = 0;
        if ( ok ) { cn = jnod; break; }
      }
      if ( cn<0 ) {
        cn = ++max_node;
        db( NODE, cn, idum, centre, ndim, VERSION_NORMAL, PUT );
        db( NODE_START_REFINED, cn, idum, centre, ndim,
          VERSION_NORMAL, PUT );
        double *ndof = db_dbl( NODE_DOF, q8[1], VERSION_NORMAL );
        long int ln = db_len( NODE_DOF, q8[1], VERSION_NORMAL );
        db( NODE_DOF, cn, idum, ndof, ln, VERSION_NORMAL, PUT );
        double *ndof_sr = db_dbl( NODE_DOF_START_REFINED, q8[1],
          VERSION_NORMAL );
        long int ln_sr = db_len( NODE_DOF_START_REFINED, q8[1],
          VERSION_NORMAL );
        db( NODE_DOF_START_REFINED, cn, idum, ndof_sr, ln_sr,
          VERSION_NORMAL, PUT );
        length = 1;
        db( NODE_MACRO_GENERATE, cn, &icontrol, ddum, length,
          VERSION_NORMAL, PUT );
      }
      // side-1 nodes in the hex18 (quad9 tensor) order
      long int s1[9];
      s1[0]=q8[1]; s1[1]=q8[5]; s1[2]=q8[2]; s1[3]=q8[6]; s1[4]=cn;
      s1[5]=q8[7]; s1[6]=q8[3]; s1[7]=q8[8]; s1[8]=q8[4];
      // normal = cross product of two side-1 face edges (corners
      // BL,BR,TL of the quad8), tangent = the first edge
      {
        ca = db_dbl( NODE, q8[1], VERSION_NORMAL );
        cb = db_dbl( NODE, q8[2], VERSION_NORMAL );
        double *cc = db_dbl( NODE, q8[3], VERSION_NORMAL );
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
      // CONVERGENCE (2026-09-07): zero-thickness side-2 copies. The
      // Professional .dbs of interface_quad8_hex20 duplicates the 9
      // side-1 nodes AT THE SAME COORDINATES; the old 0.01 shift put
      // side 2 beyond side 1 along +n, interface_element() then flipped
      // the normal (dir = n.(cm1-cm2) < 0) and the interface reported
      // COMPRESSION WITH POSITIVE STRESS (+1.005 vs Professional -1.0).
      shift = 0.;
      // create the 9 side-2 nodes: coincident copies of side 1
      for ( j=0; j<ns1q; j++ ) {
        long int src = s1[j];
        long int dst = ++max_node;
        db( NODE, src, idum, coord, ldum, VERSION_NORMAL, GET );
        for ( i=0; i<3; i++ ) coord[i] += shift*normal[i];
        db( NODE, dst, idum, coord, ldum, VERSION_NORMAL, PUT );
        db( NODE_START_REFINED, src, idum, coord, ldum,
          VERSION_NORMAL, GET );
        for ( i=0; i<3; i++ ) coord[i] += shift*normal[i];
        db( NODE_START_REFINED, dst, idum, coord, ldum,
          VERSION_NORMAL, PUT );
        double *ndof = db_dbl( NODE_DOF, src, VERSION_NORMAL );
        long int ln = db_len( NODE_DOF, src, VERSION_NORMAL );
        db( NODE_DOF, dst, idum, ndof, ln, VERSION_NORMAL, PUT );
        double *ndof_sr = db_dbl( NODE_DOF_START_REFINED, src,
          VERSION_NORMAL );
        long int ln_sr = db_len( NODE_DOF_START_REFINED, src,
          VERSION_NORMAL );
        db( NODE_DOF_START_REFINED, dst, idum, ndof_sr, ln_sr,
          VERSION_NORMAL, PUT );
        length = 1;
        db( NODE_MACRO_GENERATE, dst, &icontrol, ddum, length,
          VERSION_NORMAL, PUT );
        el[1+ns1q+j] = dst;
      }
      // rewrite: -hex18, side 1 = tensor order, side 2 = the dups
      el[0] = -HEX18;
      for ( j=0; j<ns1q; j++ ) el[1+j] = s1[j];
      length = 1 + 2*ns1q;
      db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, PUT );
      // reconnect the neighbours on the other side (same logic as the
      // linear conversions below: all 9 side-1 nodes shared + centroid
      // on the +normal side)
      {
        long int *nel_neigh = get_new_int(MAXIMUM_NODE+1);
        for ( long int iel=0; iel<=max_element_c; iel++ ) {
          if ( !db_active_index( ELEMENT, iel, VERSION_NORMAL ) ) continue;
          long int gr = 0;
          db( ELEMENT_GROUP, iel, &gr, ddum, ldum,
            VERSION_NORMAL, GET_IF_EXISTS );
          found = 0;
          for ( long int ig=0; ig<length_convert_groups; ig++ )
            if ( convert_groups[ig]==gr ) { found = 1; break; }
          if ( found ) continue;
          if ( iel==element ) continue;
          long int ln_n = 0;
          db( ELEMENT, iel, nel_neigh, ddum, ln_n, VERSION_NORMAL, GET );
          long int shared = 0;
          for ( long int jj2=0; jj2<ns1q; jj2++ ) {
            long int s2 = el[1+jj2];
            for ( long int k2=1; k2<ln_n; k2++ )
              if ( nel_neigh[k2]==s2 ) { shared++; break; }
          }
          if ( shared!=ns1q ) continue;   // NOT on the other side
          {
            double icx = 0., icy = 0., icz = 0., ncx = 0., ncy = 0.,
              ncz = 0.;
            long int nshared_nodes = 0;
            for ( long int k2=1; k2<ln_n; k2++ ) {
              double *cn2 = db_dbl( NODE, nel_neigh[k2], VERSION_NORMAL );
              ncx += cn2[0]; ncy += cn2[1]; ncz += cn2[2];
              nshared_nodes++;
            }
            if ( nshared_nodes>0 ) {
              ncx /= nshared_nodes; ncy /= nshared_nodes;
              ncz /= nshared_nodes;
              icx = 0.; icy = 0.; icz = 0.;
              for ( long int jj2=0; jj2<ns1q; jj2++ ) {
                double *cn1 = db_dbl( NODE, el[1+jj2],
                  VERSION_NORMAL );
                icx += cn1[0]; icy += cn1[1]; icz += cn1[2];
              }
              icx /= ns1q; icy /= ns1q; icz /= ns1q;
              double ddx = ncx-icx, ddy = ncy-icy, ddz = ncz-icz;
              double dot = ddx*normal[0] + ddy*normal[1] +
                ddz*normal[2];
              if ( dot<=0. ) continue;
            }
          }
          for ( long int jj2=0; jj2<ns1q; jj2++ ) {
            long int src2 = el[1+jj2];
            long int dst2 = el[1+ns1q+jj2];
            for ( long int k2=1; k2<ln_n; k2++ )
              if ( nel_neigh[k2]==src2 ) nel_neigh[k2] = dst2;
          }
          db( ELEMENT, iel, nel_neigh, ddum, ln_n, VERSION_NORMAL, PUT );
        }
        delete[] nel_neigh;
      }
      nconv++;
      continue;
    }

    // side 1 nodes: for 2D bar2 = {a,b} / bar3 = {a,b,c}; for 3D
    // tria3 = 3 nodes, quad4 = 4 nodes. They form the base of the
    // interface element.
    long int ns1 = 0;
    if      ( name==-BAR2 ) ns1 = 2;
    else if ( name==-BAR3 ) ns1 = 3;
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
    else if ( name==-BAR2 || name==-BAR3 ) {
      // 3D bar2/bar3: a LINE (only 2/3 collinear nodes) - the normal
      // cannot come from a cross product of two edges. Verified against
      // the Professional .dbs of interface_bar2_hex8: the extruded
      // quad4 interface has normal (0,1,0) for a bar2 along +x in the
      // xy-plane, which is z_hat x tangent (the extrusion direction is
      // z). The converted quad4 lives in the xy-plane; later
      // control_mesh_convert lifts it to the hex8 interface.
      ca = db_dbl( NODE, el[1], VERSION_NORMAL );
      cb = db_dbl( NODE, el[2], VERSION_NORMAL );
      for ( i=0; i<3; i++ ) tangent[i] = cb[i] - ca[i];
      array_normalize( tangent, 3 );
      normal[0] = 0.*tangent[2] - tangent[1]*1.;
      normal[1] = tangent[0]*1. - 0.*tangent[2];
      normal[2] = 0.;
      array_normalize( normal, 3 );
      if ( array_size( normal, 3 )<1.e-12 ) {
        // degenerate (bar2 along z): pick +x
        normal[0] = 1.; normal[1] = 0.; normal[2] = 0.;
      }
    }
    else {
      // 3D tria3/quad4: normal = cross product of two side-1 edges
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

    // CONVERGENCE (2026-09-07, interface_tria3_prism6): the side-2
    // copies are created COINCIDENT with side 1 (zero-thickness
    // interface) and ONCE PER INTERFACE NODE (per-corner pairs shared by
    // every interface element that contains the node). The Professional
    // .dbs of the whole converted family (bar2_quad4, bar3_quad8,
    // tria3_prism6, quad4_hex8, quad8_hex20) shows coincident copies; the
    // old code shifted them 0.01 along +n and duplicated the diagonal
    // nodes of a triangulated quad a second time, reconnecting them to
    // the opposite block - the originals ended up attached to nothing but
    // the interface elements, and the converged interface stress was
    // +1.51 (compression POSITIVE, flipped normal) with a 1/6-1/3-1/3-1/6
    // nodal split instead of the Professional -1.0 with 1/4 per node.
    // Per-corner rule (Professional layout): the ORIGINAL node stays with
    // the block on the -normal side of the FIRST interface element that
    // contains it (lowest element number); the DUPLICATE goes to the
    // +normal side. Both blocks share the node at the input; the
    // reconnection below moves the duplicate into every volume element on
    // the +normal side that still references the original.
    shift = 0.;

    // face centroid of the element (side-1 nodes)
    double fc[MDIM];
    array_set( fc, 0., MDIM );
    for ( j=0; j<ns1; j++ ) {
      double *ccf = db_dbl( NODE, el[1+j], VERSION_NORMAL );
      for ( i=0; i<ndim; i++ ) fc[i] += ccf[i]/(double)ns1;
    }
    long int new_name = ( name==-BAR2 ) ? -QUAD4 :
      ( name==-BAR3 ) ? -QUAD6 : ( name==-TRIA3 ) ? -PRISM6 : -HEX8;
    long int *nel_neigh = get_new_int(MAXIMUM_NODE+1);
    for ( j=0; j<ns1; j++ ) {
      long int c = el[1+j];
      if ( !node_done[c] ) {
        node_done[c] = 1;
        long int dst = ++max_node;
        node_dup[c] = dst;
        db( NODE, c, idum, coord, ldum, VERSION_NORMAL, GET );
        db( NODE, dst, idum, coord, ldum, VERSION_NORMAL, PUT );
        db( NODE_START_REFINED, c, idum, coord, ldum, VERSION_NORMAL, GET );
        db( NODE_START_REFINED, dst, idum, coord, ldum, VERSION_NORMAL, PUT );
        double *ndof = db_dbl( NODE_DOF, c, VERSION_NORMAL );
        long int ln = db_len( NODE_DOF, c, VERSION_NORMAL );
        db( NODE_DOF, dst, idum, ndof, ln, VERSION_NORMAL, PUT );
        double *ndof_sr = db_dbl( NODE_DOF_START_REFINED, c, VERSION_NORMAL );
        long int ln_sr = db_len( NODE_DOF_START_REFINED, c, VERSION_NORMAL );
        db( NODE_DOF_START_REFINED, dst, idum, ndof_sr, ln_sr,
          VERSION_NORMAL, PUT );
        length = 1;
        db( NODE_MACRO_GENERATE, dst, &icontrol, ddum, length,
          VERSION_NORMAL, PUT );
        // record sides of the pair (side 1 = lowest-numbered volume side)
        if ( side1_minus[element] ) { node_s1[c] = c;  node_s2[c] = dst; }
        else                        { node_s1[c] = dst; node_s2[c] = c;  }
        // reconnect the volumes on the +normal side of THIS element that
        // still contain the original node (other interface elements keep
        // their original records - they are rewritten by their own
        // conversion with the same per-node pairs)
        for ( long int iel=0; iel<=max_element_c; iel++ ) {
          if ( !db_active_index( ELEMENT, iel, VERSION_NORMAL ) ) continue;
          long int gr = 0;
          db( ELEMENT_GROUP, iel, &gr, ddum, ldum,
            VERSION_NORMAL, GET_IF_EXISTS );
          if ( db_active_index( GROUP_INTERFACE, gr, VERSION_NORMAL ) )
            continue;
          found = 0;
          for ( long int ig=0; ig<length_convert_groups; ig++ )
            if ( convert_groups[ig]==gr ) { found = 1; break; }
          if ( found ) continue;      // keep-side group of the control
          long int ln_n = 0;
          db( ELEMENT, iel, nel_neigh, ddum, ln_n, VERSION_NORMAL, GET );
          long int has_c = 0;
          for ( long int k2=1; k2<ln_n; k2++ )
            if ( nel_neigh[k2]==c ) { has_c = 1; break; }
          if ( !has_c ) continue;
          // centroid side test: only the +normal-side volumes move to the
          // duplicate (the -normal side keeps the original node)
          double ncx = 0., ncy = 0., ncz = 0.;
          long int nsh = ( ln_n>1 ) ? ln_n-1 : 1;
          for ( long int k2=1; k2<ln_n; k2++ ) {
            double *cn2 = db_dbl( NODE, nel_neigh[k2], VERSION_NORMAL );
            ncx += cn2[0]; ncy += cn2[1]; ncz += cn2[2];
          }
          ncx /= nsh; ncy /= nsh; ncz /= nsh;
          double dot = (ncx-fc[0])*normal[0] + (ncy-fc[1])*normal[1] +
            (ncz-fc[2])*normal[2];
          if ( dot<=0. ) continue;
          for ( long int k2=1; k2<ln_n; k2++ )
            if ( nel_neigh[k2]==c ) nel_neigh[k2] = dst;
          db( ELEMENT, iel, nel_neigh, ddum, ln_n, VERSION_NORMAL, PUT );
        }
      }
      // record slot of the corner: (side-1 node, side-2 node) of the pair
      el[1+j]     = node_s1[c];
      el[1+ns1+j] = node_s2[c];
    }
    delete[] nel_neigh;

    // 3D surface conversions only: order side 1 counter-clockwise around
    // the interface normal pointing from the side-1 block to the side-2
    // block. interface_element() derives the normal of a zero-thickness
    // interface as the cross product of the first three side-1 nodes (the
    // geometric dir-flip of a thick interface never triggers: dir=0), so
    // the record must encode the orientation. A transposition of slots 0
    // and 1 reverses the cross and preserves the per-slot pairs (the
    // integration weights of the side are uniform 1/ns1 in 3D).
    if ( ndim==3 && ( name==-TRIA3 || name==-QUAD4 ) ) {
      double *c0 = db_dbl( NODE, el[1], VERSION_NORMAL );
      double *c1 = db_dbl( NODE, el[2], VERSION_NORMAL );
      double *c2 = db_dbl( NODE, el[3], VERSION_NORMAL );
      double e1c[MDIM], e2c[MDIM], ncross[MDIM];
      for ( i=0; i<3; i++ ) { e1c[i] = c1[i]-c0[i]; e2c[i] = c2[i]-c0[i]; }
      ncross[0] = e1c[1]*e2c[2] - e1c[2]*e2c[1];
      ncross[1] = e1c[2]*e2c[0] - e1c[0]*e2c[2];
      ncross[2] = e1c[0]*e2c[1] - e1c[1]*e2c[0];
      // side 1 = the -normal side (side1_minus=1): the normal points
      // along +n; otherwise it points along -n
      double sgn = side1_minus[element] ? 1. : -1.;
      double dot = sgn * ( ncross[0]*normal[0] + ncross[1]*normal[1] +
        ncross[2]*normal[2] );
      if ( dot<0. ) {
        long int t;
        t = el[1]; el[1] = el[2]; el[2] = t;
        t = el[1+ns1]; el[1+ns1] = el[2+ns1]; el[2+ns1] = t;
      }
    }

    // rewrite the element: bar2->quad4, bar3->quad6, tria3->prism6,
    // quad4->hex8. el[0]=name, el[1..ns1]=side1,
    // el[ns1+1..2*ns1]=side2.
    el[0] = new_name;
    length = 1 + 2*ns1;
    db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, PUT );
    nconv++;
  }
  delete[] el;
  delete[] convert_groups;
  delete[] side1_minus;
  delete[] node_done;
  delete[] node_dup;
  delete[] node_s1;
  delete[] node_s2;
  // mesh_has_changed: ALWAYS called after a conversion (the 3D extrude
  // now runs BEFORE the convert, so by the time the convert lifts the
  // extruded quad4 interface to hex8 the solids are already 3D and
  // area_element_group processes valid hex8 - the old order (convert
  // first) left 2D quad4/bar2 that area_element_group could not handle).
  if ( nconv>0 )
    mesh_has_changed( VERSION_NORMAL );

  if ( swit ) pri( "Out function INTERFACE_CONVERT" );
}
