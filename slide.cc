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

void slide_spring( long int inod, long int islide, double *new_node_dof,
  double *node_lhside, double *node_rhside, double normal[MDIM],
  double dtime, long int icontrol, long int swit );
void slide_penalty_law( long int inod, long int islide, double *new_node_dof,
  double *node_lhside, double *node_rhside, double normal[MDIM],
  double dtime, long int swit );

// slide() implements the slide_geometry family (manual Professional
// 6.1041-6.1047 + 6.370-6.372):
//
// LEGACY path (no slide_stiffness record for the slide index): velocity
// penalty constraint - the node is not allowed to move normal to the
// slide plane (slide_penalty), with an optional Coulomb friction
// (slide_friction, slide_axisymmetric scales it by 2*pi*r). Kept for
// the classic inputs (e.g. examp14) that only use slide_geometry.
//
// PROFESSIONAL path (slide_stiffness kn kt present): elastic-plastic
// slide springs on the TOTAL displacements, per node, following the
// support_edge_normal pattern (area.cc):
//   - the normal spring reacts on the total normal displacement
//     u_n = u.n  (the node started ON the slide plane, so the plane is
//     the fixed reference):  F_n = -k_n*u_n  on the node;
//   - the tangential spring reacts on the total tangential displacement
//     u_t = u - u_n n: elastic predictor F_t = -k_t*u_t capped by the
//     Mohr-Coulomb law  |F_t| <= c + F_n*tan(phi)  (slide_plasti_friction)
//     and by the tension cut-off of the connection
//     (slide_plasti_tension sig_t: maximum pull-off force; absent =
//     compression-only/gap);
//   - slide_plasti_residual_stiffness (knr ktr): FRACTION of the elastic
//     stiffness added to the matrix only (not the RHS) once the
//     corresponding yield branch is active - it only steadies the
//     iterations (the Professional writes a default 1e-2 1e-2);
//   - in axisymmetric problems (problem group_axisymmetric -yes or the
//     slide_axisymmetric record) every force/stiffness acts on the whole
//     ring: multiplied by 2*pi*r (r = radial node coordinate), like the
//     element assembly.
// The matrix terms are the per-direction diagonal (|n_i|/|t_i| weights)
// scaled by dtime, the standard node-based assembly (contact.cc).
// Output records per slide node (Professional .dbs): node_slide_direction
// = (n, t) with t the friction direction ON the node, node_slide_f = the
// plastic yield function value and node_slide_force = (force on the slide
// geometry along n, force on the slide geometry along t) = -force on the
// node in the local frame.

void slide( void )

{
  long int inod=0, max_node=0, idim=0, swit=0, islide=0, max_slide=0,
    in_geometry=0, ldum=0, icontrol=0, idum[1], slide_geometry[2];
  double dtime=0., tmp=0., rdum=0., radius=0., ddum[MDIM], 
    *new_node_dof=NULL, *node_lhside=NULL, *node_rhside=NULL;

  if ( db_max_index( SLIDE_GEOMETRY, max_slide, VERSION_NORMAL, GET ) >= 0 ) {
    swit = set_swit(-1,-1,"slide");
    if ( swit ) pri( "In routine SLIDE" );
    db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET_IF_EXISTS );
    db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    for ( islide=0; islide<=max_slide; islide++ ) {
      if ( db_active_index( SLIDE_GEOMETRY, islide, VERSION_NORMAL ) ) {
        long int has_spring = db_active_index( SLIDE_STIFFNESS, islide,
          VERSION_NORMAL );
        db( SLIDE_GEOMETRY, islide, slide_geometry, ddum, 
          ldum, VERSION_NORMAL, GET );
        db_max_index( NODE, max_node, VERSION_NORMAL, GET );
        for ( inod=0; inod<=max_node; inod++ ) {
          if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
            // node_slide: ADDITIONAL membership by explicit node record
            // (manual Professional 6.893): the node belongs to the sliding
            // geometry with index islide when its node_slide record
            // carries slide_number islide, EVEN IF the node is not inside
            // the geometry itself. The plane normal still comes from
            // slide_geometry (for all members).
            long int node_slide_member = 0, node_slide_number = -1;
            double normal[MDIM];
            if ( db( NODE_SLIDE, inod, &node_slide_number, ddum, ldum,
                 VERSION_NORMAL, GET_IF_EXISTS ) )
              if ( node_slide_number==islide ) node_slide_member = 1;
            // The membership test uses the START (undeformed) node
            // coordinates: the nodes belong to the plane where they were
            // created even after they slid far along it.
            geometry( inod, ddum, slide_geometry, in_geometry, rdum, normal,
              rdum, ddum, NODE_START_REFINED, PROJECT_EXACT, VERSION_NORMAL );
            if ( node_slide_member ) in_geometry = 1;
            if ( in_geometry ) {
              new_node_dof = db_dbl( NODE_DOF, inod, VERSION_NEW );
              node_lhside = db_dbl( NODE_LHSIDE, inod, VERSION_NORMAL );
              node_rhside = db_dbl( NODE_RHSIDE, inod, VERSION_NORMAL );
              if ( has_spring )
                slide_spring( inod, islide, new_node_dof, node_lhside,
                  node_rhside, normal, dtime, icontrol, swit );
              else
                slide_penalty_law( inod, islide, new_node_dof, node_lhside,
                  node_rhside, normal, dtime, swit );
            }
          }
        }
      }
    }
    if ( swit ) pri( "Out routine SLIDE" );
  }

}

// Professional elastic-plastic slide spring law (per node).
void slide_spring( long int inod, long int islide, double *new_node_dof,
  double *node_lhside, double *node_rhside, double normal[MDIM],
  double dtime, long int icontrol, long int swit )

{
  long int idim=0, ldum=0, idum[1];
  double ddum[MDIM], values[DATA_ITEM_SIZE];
  double kn=0., kt=0., un=0., ut_mag=0., S=0., cap_t=0., c=0., mu=0.,
    sig_t=0., res_t=1.e-2, ax=1., radius=0.,
    tangent[MDIM], fn_node[MDIM], ft_node[MDIM], *coord=NULL;
  long int plasti_apply=-YES, stiffness_apply=-YES, has_friction=0,
    axisym=-NO;

  if ( !materi_velocity ) return;
  // control gates (manual Professional 6.371/6.372): -no disables the
  // plastic caps respectively the slide_stiffness springs.
  db( CONTROL_SLIDE_PLASTI_APPLY, icontrol, &plasti_apply, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_SLIDE_STIFFNESS_APPLY, icontrol, &stiffness_apply, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  if ( stiffness_apply==-NO ) return;

  db( SLIDE_STIFFNESS, islide, idum, values, ldum, VERSION_NORMAL, GET );
  kn = values[0];
  kt = values[1];
  // plastic friction phi c (manual 6.1042) and tension cut-off sig_t
  // (manual 6.1043): the records of the slide with the same index.
  if ( db( SLIDE_PLASTI_FRICTION, islide, idum, values, ldum,
      VERSION_NORMAL, GET_IF_EXISTS ) ) {
    c = values[1];
    mu = tan( values[0] );
    has_friction = 1;
  }
  if ( db( SLIDE_PLASTI_TENSION, islide, idum, values, ldum,
      VERSION_NORMAL, GET_IF_EXISTS ) )
    sig_t = values[0];
  // slide_plasti_residual_stiffness (tangential fraction of the elastic
  // stiffness, matrix only, active while the node slips); the
  // Professional defaults to 1e-2 when the record is absent.
  if ( db( SLIDE_PLASTI_RESIDUAL_STIFFNESS, islide, idum, values, ldum,
      VERSION_NORMAL, GET_IF_EXISTS ) )
    res_t = values[1];
  // axisymmetric: the slide acts on the ring of circumference 2*pi*r.
  // Triggered by the slide_axisymmetric record (legacy, index 0) or by a
  // problem group with group_axisymmetric -yes (the corpus slide2/slide3).
  db( SLIDE_AXISYMMETRIC, 0, &axisym, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  if ( axisym!=-YES ) {
    long int max_grp=0, gr=0;
    db_max_index( GROUP_TYPE, max_grp, VERSION_NORMAL, GET );
    for ( gr=0; gr<=max_grp && axisym!=-YES; gr++ )
      if ( db_active_index( GROUP_TYPE, gr, VERSION_NORMAL ) )
        db( GROUP_AXISYMMETRIC, gr, &axisym, ddum, ldum,
          VERSION_NORMAL, GET_IF_EXISTS );
  }
  if ( axisym==-YES ) {
    coord = db_dbl( NODE, inod, VERSION_NORMAL );
    radius = coord[0];
    ax = 2. * PIRAD * radius;
  }

  // total displacement of the node: the dis dofs (materi_displacement)
  // or the integrated velocities (materi_velocity_integrated). Without
  // either there is no displacement to spring on - nothing to do.
  long int displ_indx = -1;
  if ( materi_displacement ) displ_indx = dis_indx;
  else if ( materi_velocity_integrated ) displ_indx = veli_indx;
  if ( displ_indx<0 ) return;

  array_set( tangent, 0., MDIM );
  array_set( fn_node, 0., MDIM );
  array_set( ft_node, 0., MDIM );
  un = 0.;
  for ( idim=0; idim<ndim; idim++ )
    un += new_node_dof[displ_indx+idim*nder] * normal[idim];
  for ( idim=0; idim<ndim; idim++ )
    tangent[idim] = new_node_dof[displ_indx+idim*nder] - un*normal[idim];
  ut_mag = array_size( tangent, ndim );
  // the tangential FRICTION direction on the node opposes the tangential
  // displacement (the Professional node_slide_direction record)
  for ( idim=0; idim<ndim; idim++ ) tangent[idim] *= -1.;
  if ( ut_mag>1.e-12 ) array_normalize( tangent, ndim );

  // normal force (compression positive): S = -k_n*u_n on the node. The
  // normal spring is TWO-SIDED (the node is elastically attached to the
  // slide plane); slide_plasti_tension (manual 6.1044) only caps the
  // maximum TENSILE force the connection can take (S >= -sig_t) - absent
  // means no tension limit. S enters the friction cap and the applied
  // normal force.
  S = -kn*un;
  if ( S<-sig_t ) S = -sig_t;
  // friction cap c + Fn*tan(phi), zero when the cap turns negative
  // (contact.cc convention); without the plasti_friction record, or with
  // control_slide_plasti_apply -no, the tangential spring stays elastic
  cap_t = 1.e300;
  if ( has_friction && plasti_apply!=-NO ) {
    cap_t = c + mu*S;
    if ( cap_t<0. ) cap_t = 0.;
  }

  // the normal force on the node (out of the plane when compressed) and
  // the tangential force on the node: the elastic predictor magnitude
  // k_t*|u_t| capped by the Mohr-Coulomb limit. In the SLIP branch the
  // tangential force is capped (constant force) and the matrix keeps only
  // the residual fraction of the tangential stiffness (the elastic
  // tangential stiffness in the matrix biases the 2-iteration steps of
  // the velocity solver - measured on slide1). The residual fraction is
  // added to the matrix only (slide_plasti_residual_stiffness, default
  // 1e-2 1e-2 like the Professional) for stable iterations.
  for ( idim=0; idim<ndim; idim++ ) {
    fn_node[idim] = ax*S*normal[idim];
    if ( ut_mag>1.e-12 ) {
      double ft_el = ax*kt*ut_mag;          // elastic predictor magnitude
      double ft_cap = ax*cap_t;             // cap on the ring force
      double ft_mag = ft_el;
      if ( ft_el>ft_cap ) {
        ft_mag = ft_cap;
        node_lhside[vel_indx/nder+idim] +=
          dtime*ax*res_t*kt*scalar_dabs( tangent[idim] );
      }
      else
        node_lhside[vel_indx/nder+idim] +=
          dtime*ax*kt*scalar_dabs( tangent[idim] );
      ft_node[idim] = ft_mag*tangent[idim];
    }
    else
      // node sits at the origin of the tangent frame: elastic spring
      node_lhside[vel_indx/nder+idim] +=
        dtime*ax*kt*( 1. - scalar_dabs( normal[idim] ) );
    // normal spring in the matrix (diagonal |n_i| weights, contact.cc
    // pattern)
    node_lhside[vel_indx/nder+idim] +=
      dtime*ax*kn*scalar_dabs( normal[idim] );
    node_rhside[vel_indx/nder+idim] += fn_node[idim] + ft_node[idim];
  }
  {
    double fout[MDIM], dirs[6], fval = 0., fric_mag = 0.;
    long int lfout = ndim, ldirs = 6, lfval = 1;
    array_set( fout, 0., MDIM );
    array_set( dirs, 0., 6 );
    for ( idim=0; idim<ndim; idim++ ) fric_mag += ft_node[idim]*ft_node[idim];
    fric_mag = sqrt( fric_mag );
    // node_slide_force = the force the material applies ON the slide
    // geometry (= -force on the node) in the local (n, t) frame:
    // slot 0 = k_n*u_n (compression negative), slot 1 = - tangential
    // force magnitude on the node.
    fout[0] = -S*ax;
    if ( ut_mag>1.e-12 ) fout[1] = -fric_mag;
    if ( has_friction && plasti_apply!=-NO )
      fval = fric_mag/ax - ( c + mu*S );
    for ( idim=0; idim<ndim; idim++ ) {
      dirs[idim] = normal[idim];
      dirs[3+idim] = tangent[idim];
    }
    db( NODE_SLIDE_DIRECTION, inod, idum, dirs, ldirs, VERSION_NEW, PUT );
    db( NODE_SLIDE_F, inod, idum, &fval, lfval, VERSION_NEW, PUT );
    db( NODE_SLIDE_FORCE, inod, idum, fout, lfout, VERSION_NEW, PUT );
    if ( swit ) pri( "node_slide_force", fout, ndim );
  }
}

// Legacy slide penalty law (velocity constraint + friction) for slide
// geometries WITHOUT the Professional slide_stiffness record.
void slide_penalty_law( long int inod, long int islide, double *new_node_dof,
  double *node_lhside, double *node_rhside, double normal[MDIM],
  double dtime, long int swit )

{
  long int idim=0, ldum=0, idum[1], slide_axisymmetric=-NO;
  double ddum[MDIM], velocity[MDIM], slide_velocity[MDIM];
  double normal_velocity=0., slide_force=0., normal_force=0.,
    slide_friction=0., slide_penalty=1.e15, tmp=0., radius=0.,
    *coord=NULL;

  db( SLIDE_FRICTION, inod, idum, &slide_friction,
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( SLIDE_PENALTY, islide, idum, &slide_penalty, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  db( SLIDE_AXISYMMETRIC, 0, &slide_axisymmetric, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  for ( idim=0; idim<ndim; idim++ ) velocity[idim] =
    new_node_dof[vel_indx+idim*nder];
  normal_velocity = array_inproduct( normal, velocity, ndim );
  for ( idim=0; idim<ndim; idim++ ) {
    tmp = velocity[idim] - normal_velocity * normal[idim];
    slide_velocity[idim] = tmp;
  }
  normal_force = array_inproduct( node_rhside, normal, ndim );
  slide_force = slide_friction * normal_force;
  // slide_axisymmetric -yes: the slide friction acts on the
  // whole ring of circumference 2*pi*r (r = radial distance of
  // the node to the axis, the x coordinate in axisymmetric).
  if ( slide_axisymmetric==-YES ) {
    coord = db_dbl( NODE, inod, VERSION_NORMAL );
    radius = coord[0];
    slide_force *= 2. * PIRAD * radius;
  }
  if ( array_normalize( slide_velocity, ndim ) ) {
    for ( idim=0; idim<ndim; idim++ ) {
      node_lhside[vel_indx+idim*nder] += slide_penalty * dtime;
      node_rhside[vel_indx+idim*nder] +=
        - slide_penalty * normal_velocity * dtime * normal[idim]
        - slide_force * slide_velocity[idim];
    }
  }
}
