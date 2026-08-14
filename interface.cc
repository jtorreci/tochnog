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
//   - gap (Fase 3): the interface only generates stresses when the
//       accumulated normal strain < gap (closed). Otherwise the stiffness
//       is the residual one.
//   - tension limit (Fase 3):
//       group_interface_materi_plasti_tension_direct tension_limit
//   - Mohr-Coulomb (Fase 3):
//       group_interface_materi_plasti_mohr_coul_direct phi c phi_flow
//       max friction force = c + Fn * tan(phi)
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
    nnol=4, idum[1];
  double dtime=0., kn=0., kt1=0., kt2=0., tmp=0., ddum[1],
    normal[MDIM], tangent[MDIM], du[MDIM],
    du_norm=0., du_tang=0., stress_normal=0., stress_shear=0.,
    strain_normal=0., force_norm=0., force_tang=0., gap=0., tension_limit=0.,
    residual_factor=0.01, phi=0., c=0., phi_flow=0., max_fric=0.,
    stiff_normal=0., stiff_tang=0., ddum3[3];

  swit = set_swit(element,-1,"interface_element");
  if ( swit ) pri( "In routine INTERFACE_ELEMENT." );

  if ( !db_active_index( GROUP_INTERFACE, element_group, VERSION_NORMAL ) )
    return;

  db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET );

  // group parameters
  db( GROUP_INTERFACE_MATERI_ELASTI_STIFFNESS, element_group, idum, ddum3,
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  kn = ddum3[0]; kt1 = ddum3[1]; kt2 = ddum3[2];
  db( GROUP_INTERFACE_MATERI_PLASTI_TENSION_DIRECT, element_group, idum,
    &tension_limit, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( GROUP_INTERFACE_MATERI_RESIDUAL_STIFFNESS, element_group, idum,
    &residual_factor, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( GROUP_INTERFACE_MATERI_PLASTI_MOHR_COUL_DIRECT, element_group, idum,
    ddum3, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  phi = ddum3[0]; c = ddum3[1]; phi_flow = ddum3[2];

  // 2D interface: sides are node pairs {0,1} and {2,3}. A -bar2 element
  // reaching here means it was not converted by control_mesh_convert;
  // treat it as a degenerate interface (side 1 = node 0, side 2 = node 1)
  // rather than aborting.
  if ( name==-BAR2 ) {
    nnol = 2;
  }
  else {
    assert( name==-QUAD4 );
    nnol = 4;
  }

  // normal to the interface: perpendicular to side 1 (nodes 0,1)
  if ( ndim==2 ) {
    tangent[0] = coord[1*ndim+0] - coord[0*ndim+0];
    tangent[1] = coord[1*ndim+1] - coord[0*ndim+1];
    array_normalize( tangent, ndim );
    normal[0] = -tangent[1];
    normal[1] =  tangent[0];
  }
  else {
    array_set( normal, 0., MDIM );
    normal[2] = 1.;
    tangent[0] = 1.; tangent[1] = 0.; tangent[2] = 0.;
  }

  // velocity difference between the sides (side2 - side1) on the
  // velocity dof -> incremental displacement difference (spring2 pattern)
  // For a -bar2 (not converted): node 0 = side 1, node 1 = side 2.
  if ( name==-BAR2 ) {
    for ( idim=0; idim<ndim; idim++ ) {
      double v1 = new_dof[0*nuknwn+vel_indx+idim*nder];
      double v2 = new_dof[1*nuknwn+vel_indx+idim*nder];
      du[idim] = ( v2 - v1 ) * dtime;
    }
  }
  else {
    for ( idim=0; idim<ndim; idim++ ) {
      double v_side1 = 0.5*( new_dof[0*nuknwn+vel_indx+idim*nder] +
                             new_dof[1*nuknwn+vel_indx+idim*nder] );
      double v_side2 = 0.5*( new_dof[2*nuknwn+vel_indx+idim*nder] +
                             new_dof[3*nuknwn+vel_indx+idim*nder] );
      du[idim] = ( v_side2 - v_side1 ) * dtime;
    }
  }

  du_norm  = array_inproduct( du, normal, ndim );
  du_tang  = array_inproduct( du, tangent, ndim );

  // accumulated normal strain (history) - used only to decide closed/open
  db( ELEMENT_INTERFACE_STRAIN_NORMAL, element, idum, &strain_normal,
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  strain_normal += du_norm;

  // normal force: INCREMENTAL (validated Fase 1 approach, spring2 pattern:
  // force = kn*du). The stiffness depends on the closed/open state: if
  // the interface is opened (accumulated strain >= gap), only the residual
  // stiffness acts. If gap is not specified the interface is always
  // closed (gap = +1e20 -> never opens).
  if ( !db( GROUP_INTERFACE_GAP, element_group, idum, &gap, ldum,
      VERSION_NORMAL, GET_IF_EXISTS ) )
    gap = 1.e20;
  stiff_normal = kn;
  if ( strain_normal >= gap ) {
    stiff_normal = kn * residual_factor;
  }
  force_norm = stiff_normal * du_norm;
  // tension limit: if the interface opens in tension, switch to residual
  if ( tension_limit>0. && strain_normal>0. && force_norm>tension_limit ) {
    stiff_normal = kn * residual_factor;
    force_norm = tension_limit;
  }
  stress_normal = force_norm;

  // tangential force: incremental elastic, limited by Mohr-Coulomb.
  // max friction = c + Fn*tan(phi), where Fn = kn*strain_normal is the
  // total normal force (compression negative).
  force_tang = kt1 * 2. * du_tang;
  stiff_tang = kt1 * 2.;
  if ( phi>0. || c>0. ) {
    double fn = kn * strain_normal;
    max_fric = c + fn * tan(phi);
    if      ( force_tang>  max_fric ) { force_tang =  max_fric; stiff_tang = 0.; }
    else if ( force_tang< -max_fric ) { force_tang = -max_fric; stiff_tang = 0.; }
  }
  stress_shear = force_tang;

  if ( swit ) {
    pri( "du_norm", du_norm );
    pri( "du_tang", du_tang );
    pri( "strain_normal", strain_normal );
    pri( "stress_normal", stress_normal );
    pri( "force_tang", force_tang );
  }

  // assembly: nodal force -sign*(stress*dir) and stiffness matrix on the
  // velocity dofs (pattern spring.cc)
  for ( idim=0; idim<ndim; idim++ ) {
    double dirn = normal[idim], dirt = tangent[idim];
    for ( inol=0; inol<nnol; inol++ ) {
      double sign = ( inol>=2 ) ? +1. : -1.;
      indx = inol*npuknwn + (vel_indx+idim*nder)/nder;
      tmp = -sign*( stress_normal*dirn + stress_shear*dirt );
      element_rhside[indx] += tmp;
      for ( jnol=0; jnol<nnol; jnol++ ) {
        double jsign = ( jnol>=2 ) ? +1. : -1.;
        for ( jdim=0; jdim<ndim; jdim++ ) {
          double jdirn = normal[jdim], jdirt = tangent[jdim];
          double kkk = sign*jsign*( stiff_normal*dirn*jdirn +
            stiff_tang*dirt*jdirt );
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

  // store accumulated normal strain history (used by gap / Mohr-Coulomb)
  ldum = 1;
  db( ELEMENT_INTERFACE_STRAIN_NORMAL, element, idum, &strain_normal,
    ldum, VERSION_NEW, PUT );

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
    if ( name!=-BAR2 ) continue;
    if ( !db_active_index( GROUP_INTERFACE, element_group, VERSION_NORMAL ) )
      continue;

    // side 1 nodes {a,b} = el[1], el[2]
    long int a = el[1], b = el[2];

    // tangent along a-b, normal perpendicular (2D)
    array_set( tangent, 0., MDIM );
    if ( ndim==2 ) {
      ca = db_dbl( NODE, a, VERSION_NORMAL );
      cb = db_dbl( NODE, b, VERSION_NORMAL );
      tangent[0] = cb[0] - ca[0];
      tangent[1] = cb[1] - ca[1];
      array_normalize( tangent, ndim );
      normal[0] = -tangent[1];
      normal[1] =  tangent[0];
    }
    else {
      array_set( normal, 0., MDIM );
      normal[2] = 1.;
    }

    // shift: small fraction of the side length (interface thickness)
    double len = 0.;
    for ( i=0; i<ndim; i++ ) len += tangent[i]*tangent[i];
    shift = 0.01 * sqrt( len );

    // create the two new nodes {a',b'} = copies of {a,b} shifted along n
    long int an = ++max_node, bn = ++max_node;
    for ( j=0; j<2; j++ ) {
      long int src = ( j==0 ) ? a : b;
      long int dst = ( j==0 ) ? an : bn;
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
    }

    // rewrite element as -quad4 {a,b,a',b'}
    el[0] = -QUAD4;
    el[1] = a; el[2] = b; el[3] = an; el[4] = bn;
    length = 5;
    db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, PUT );

    // reconnect neighbours on the OTHER side: elements sharing {a,b}
    // that are NOT in convert_groups get {a,b} replaced by {a',b'}
    for ( j=0; j<2; j++ ) {
      long int src = ( j==0 ) ? a : b;
      long int dst = ( j==0 ) ? an : bn;
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
