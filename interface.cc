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

// interface_element - interface elements (Carril A, Fase 1).
//
// An interface element models a joint/discontinuity between two blocks
// of material. Its strains are the DISPLACEMENT DIFFERENCES between the
// two opposite sides of the element (not field gradients). In 2D the
// element is a quadrilateral with 4 nodes: nodes {0,1} form side 1 and
// nodes {2,3} form side 2.
//
// Elastic law (group_interface_materi_elasti_stiffness kn kt,first
// kt,second):
//   stress_normal = kn * strain_normal
//   stress_shear1 = kt,first * shear_gamma1   (shear_gamma = 2*strain)
//   stress_shear2 = kt,second * shear_gamma2
// where strain = displacement difference between the sides divided by
// the interface thickness (taken as 1, so kn/kt absorb the thickness).
void interface_element( long int element, long int name,
  long int element_group, double coord[], double old_dof[], double new_dof[], 
  double element_lhside[], double element_matrix[], double element_rhside[] )

{
  long int idim=0, jdim=0, inol=0, jnol=0, indx=0, swit=0, ldum=0, 
    nnol=4, idum[1];
  double dtime=0., kn=0., kt1=0., kt2=0., tmp=0., ddum[1],
    normal[MDIM], tangent[MDIM], mid[MDIM], du[MDIM],
    du_norm=0., du_tang=0., stress_normal=0., stress_shear1=0.,
    stress_shear2=0., stiff[MDIM];

  swit = set_swit(element,-1,"interface_element");
  if ( swit ) pri( "In routine INTERFACE_ELEMENT." );

  if ( !db_active_index( GROUP_INTERFACE, element_group, VERSION_NORMAL ) )
    return;

  db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET );

  // stiffness kn, kt,first, kt,second (defaults to 0 if not given)
  db( GROUP_INTERFACE_MATERI_ELASTI_STIFFNESS, element_group, idum, ddum,
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  // note: db with GET_IF_EXISTS and a scalar target fills ddum[0..]
  kn = ddum[0]; kt1 = ddum[1]; kt2 = ddum[2];

  // 2D interface: sides are node pairs {0,1} and {2,3}
  assert( name==-QUAD4 );

  // normal to the interface: perpendicular to the segment connecting the
  // midpoints of the two sides
  for ( idim=0; idim<ndim; idim++ ) {
    mid[idim] = 0.5*( coord[0*ndim+idim] + coord[1*ndim+idim]
      + coord[2*ndim+idim] + coord[3*ndim+idim] )/2.;
  }
  if ( ndim==2 ) {
    // tangent along side 1
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

  // displacement difference between the sides (side2 - side1).
  // old_dof is the previous displacement, new_dof the current one; the
  // incremental displacement difference drives the stress increment.
  for ( idim=0; idim<ndim; idim++ ) {
    double u_side1_old = 0.5*( old_dof[0*nuknwn+dis_indx+idim*nder] +
                               old_dof[1*nuknwn+dis_indx+idim*nder] );
    double u_side2_old = 0.5*( old_dof[2*nuknwn+dis_indx+idim*nder] +
                               old_dof[3*nuknwn+dis_indx+idim*nder] );
    double u_side1_new = 0.5*( new_dof[0*nuknwn+dis_indx+idim*nder] +
                               new_dof[1*nuknwn+dis_indx+idim*nder] );
    double u_side2_new = 0.5*( new_dof[2*nuknwn+dis_indx+idim*nder] +
                               new_dof[3*nuknwn+dis_indx+idim*nder] );
    du[idim] = ( u_side2_new - u_side1_new ) -
               ( u_side2_old - u_side1_old );
  }

  du_norm  = array_inproduct( du, normal, ndim );
  du_tang  = array_inproduct( du, tangent, ndim );

  // strains = du / thickness (thickness = 1 -> strains = du)
  // stresses from the elastic interface law
  stress_normal  = kn * du_norm;
  stress_shear1  = kt1 * 2. * du_tang;
  stress_shear2  = kt2 * 0.;
  if ( swit ) {
    pri( "du_norm", du_norm );
    pri( "du_tang", du_tang );
    pri( "stress_normal", stress_normal );
    pri( "stress_shear1", stress_shear1 );
  }

  // assembly: nodal force and stiffness matrix on the VELOCITY dofs
  // (same pattern as spring.cc). The stiffness K = [K -K; -K K] acts on
  // the velocity difference between the two sides, with K = kn (normal)
  // and kt*2 (tangential).
  stiff[0] = kn;            // normal
  stiff[1] = kt1 * 2.;      // tangential (shear)
  stiff[2] = kt2 * 2.;      // second tangential (unused in 2D)

  for ( idim=0; idim<ndim; idim++ ) {
    double dirn = normal[idim], dirt = tangent[idim];
    for ( inol=0; inol<nnol; inol++ ) {
      // sign: + on side2 nodes (2,3), - on side1 nodes (0,1)
      double sign = ( inol>=2 ) ? +1. : -1.;
      indx = inol*npuknwn + (vel_indx+idim*nder)/nder;
      tmp = sign*( stress_normal*dirn + stress_shear1*dirt );
      element_rhside[indx] += tmp;
      for ( jnol=0; jnol<nnol; jnol++ ) {
        double jsign = ( jnol>=2 ) ? +1. : -1.;
        for ( jdim=0; jdim<ndim; jdim++ ) {
          double jdirn = normal[jdim], jdirt = tangent[jdim];
          double kkk = sign*jsign*( stiff[0]*dirn*jdirn +
            stiff[1]*dirt*jdirt );
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

  if ( swit ) pri( "Out function INTERFACE_ELEMENT" );
}
