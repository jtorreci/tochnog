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

#define EPS_ISOP 1.e-3
#define EPS_SIZE 1.e-6
#define MITER 10000
#define DELTA 1.e-4

long int point_el( double point[], double coords[], double weight[],
  long int name, long int nnol )

{
  long int found=1, idim=0, inol=0, jnol=0, knol=0, lnol=0, i=0,
    inol_xi=0, inol_eta=0, inol_zeta=0, nnol_xi=1, nnol_eta=1, 
    nnol_zeta=1, iiso=0, iter=0, converged=0;
  double dist_high=0., dist_low=0., dist_old=0., dist_tmp=0.,
    element_largest_size=0.,
    tmp=0., step=1., point_distance=0., total=0.,
    L1=0., L2=0., L3=0., L4=0., h_xi[MNOL], h_eta[MNOL], 
    h_zeta[MNOL], element_coord[MDIM], 
    ddum[MPOINT*MDIM*MNOL], iso[MDIM], iso_old[MDIM], 
    iso_new[MDIM], iso_tmp[MDIM], work[MNOL];

  if ( name==-TRUSS || name==-BEAM || name==-TRUSSBEAM ||
    name==-SPRING1 || name==-SPRING2 || name==-CONTACTSPRING ) return 0;

  if      ( name==-TRIA3 || name==-TRIA6 ) {
    if ( name==-TRIA3 ) {
      inol = 0;
      jnol = 1;
      knol = 2;
    }
    else {
      assert( name==-TRIA6 );
      inol = 0;
      jnol = 2;
      knol = 5;
    }
    project_point_on_triangle( point, &coords[inol*ndim], &coords[jnol*ndim],
      &coords[knol*ndim], work );
    L1 = work[0];
    L2 = work[1];
    L3 = work[2];
    if      ( L1<-EPS_ISOP ) 
      found = 0;
    else if ( L2<-EPS_ISOP ) 
      found = 0;
    else if ( L3<-EPS_ISOP ) 
      found = 0;
    else if ( L1>(1.+EPS_ISOP) ) 
      found = 0;
    else if ( L2>(1.+EPS_ISOP) ) 
      found = 0;
    else if ( L3>(1.+EPS_ISOP) ) 
      found = 0;
    else {
      if ( name==-TRIA3 ) {
        weight[0] = L1;
        weight[1] = L2;
        weight[2] = L3;
      }
      else {
        assert( name==-TRIA6 );
        weight[0] = (2.*L1-1.)*L1;
        weight[1] = 4.*L1*L2;
        weight[2] = (2.*L2-1.)*L2;
        weight[3] = 4.*L3*L1;
        weight[4] = 4.*L2*L3;
        weight[5] = (2.*L3-1.)*L3;
      }
    }
  }
  else if ( name==-TET4 || name==-TET10 ) {
    if ( name==-TET4 ) {
      inol = 0;
      jnol = 1;
      knol = 2;
      lnol = 3;
    }
    else {
      assert( name==-TET10 );
      inol = 0;
      jnol = 2;
      knol = 5;
      lnol = 9;
    }
    total = tetrahedron_volume( &coords[inol*ndim], &coords[jnol*ndim], 
      &coords[knol*ndim], &coords[lnol*ndim] ) ;
    if ( total<1.e-10 ) {
      found = 0;
    }
    else {
      tmp = tetrahedron_volume( &coords[jnol*ndim], &coords[knol*ndim],
        &coords[lnol*ndim], point );
      L1 = tmp / total;
      tmp = tetrahedron_volume( &coords[inol*ndim], &coords[knol*ndim],
        &coords[lnol*ndim], point );
      L2 = tmp / total;
      tmp = tetrahedron_volume( &coords[inol*ndim], &coords[jnol*ndim],
        &coords[lnol*ndim], point );
      L3 = tmp / total;
      L4 = 1. - L1 - L2 -L3;
      if ( name==-TET4 ) {
        weight[0] = L1;
        weight[1] = L2;
        weight[2] = L3;
        weight[3] = 1. - L1 - L2 - L3;
      }
      else {
        assert( name==-TET10 );
        weight[0] = (2.*L1-1.)*L1;
        weight[1] = 4.*L1*L2;
        weight[2] = (2.*L2-1.)*L2;
        weight[3] = 4.*L3*L1;
        weight[4] = 4.*L2*L3;
        weight[5] = (2.*L3-1.)*L3;
        weight[6] = 4.*L1*L4;
        weight[7] = 4.*L2*L4;
        weight[8] = 4.*L3*L4;
        weight[9] = (2.*L4-1.)*L4;
      }
      for ( idim=0; idim<ndim; idim++ ) {
        work[idim] = 0.;
        for ( i=0; i<nnol; i++ )
          work[idim] += weight[i] * coords[i*ndim+idim];
      }
      if ( L1<-EPS_ISOP || L1>(1.+EPS_ISOP) ||
           L2<-EPS_ISOP || L2>(1.+EPS_ISOP) ||
           L3<-EPS_ISOP || L3>(1.+EPS_ISOP) ||
           L4<-EPS_ISOP || L4>(1.+EPS_ISOP) ||
           ( scalar_dabs(point[0]-work[0]) > EPS_SIZE ) ||
           ( scalar_dabs(point[1]-work[1]) > EPS_SIZE ) ||
           ( scalar_dabs(point[2]-work[2]) > EPS_SIZE )
         ) found = 0;
    }
  }
  else {
    if      ( name==-BAR2 ) 
      nnol_xi = 2;
    else if ( name==-BAR3 ) 
      nnol_xi = 3;
    else if ( name==-BAR4 ) 
      nnol_xi = 4;
    else if ( name==-QUAD4 ) {
      nnol_xi = 2;
      nnol_eta = 2;
    }
    else if ( name==-QUAD9) {
      nnol_xi = 3;
      nnol_eta = 3;
    }
    else if ( name==-QUAD16 ) {
      nnol_xi = 4;
      nnol_eta = 4;
    }
    else if ( name==-HEX8 ) {
      nnol_xi = 2;
      nnol_eta = 2;
      nnol_zeta = 2;
    }
    else if ( name==-HEX27) {
      nnol_xi = 3;
      nnol_eta = 3;
      nnol_zeta = 3;
    }
    else {
      assert( name==-HEX64 );
      nnol_xi = 4;
      nnol_eta = 4;
      nnol_zeta = 4;
    }
      // average element coordinate
    array_set( element_coord, 0., ndim );
    for ( inol=0; inol<nnol; inol++ ) {
      for ( idim=0; idim<ndim; idim++ )
        element_coord[idim] += coords[inol*ndim+idim]/nnol;
    }
      // maximum element size
    element_largest_size = -1.;
    for ( inol=0; inol<nnol; inol++ ) {
      tmp = array_distance( element_coord, &coords[inol*ndim], work, 
        ndim );
      if ( tmp>element_largest_size ) element_largest_size = tmp;
    }
    point_distance = array_distance( element_coord, point, 
      work, ndim );
    if ( point_distance>(1.+EPS_SIZE)*element_largest_size )
      found = 0;
    else {
        /* iterate with the iso parametric coordinates such that these
           minimize the distance between the point and the element */
      converged=0; 
      array_set( iso_old, 0., MDIM ); 
      array_set( iso_new, 0., MDIM );
      for ( iter=0; iter<=MITER; iter++ ) {
        interpolation_polynomial( iso_old[0], nnol_xi, h_xi, ddum );
        interpolation_polynomial( iso_old[1], nnol_eta, h_eta, ddum );
        interpolation_polynomial( iso_old[2], nnol_zeta, h_zeta, ddum );
        inol = 0; array_set( element_coord, 0., ndim );
        for ( inol_zeta=0; inol_zeta<nnol_zeta; inol_zeta++ ) {
          for ( inol_eta=0; inol_eta<nnol_eta; inol_eta++ ) {
            for ( inol_xi=0; inol_xi<nnol_xi; inol_xi++ ) {
              weight[inol] = h_xi[inol_xi]*h_eta[inol_eta]*
                h_zeta[inol_zeta];
              for ( idim=0; idim<ndim; idim++ ) element_coord[idim] += 
                weight[inol] * coords[inol*ndim+idim];
              inol++;
            }
          }
        }
        if ( converged )
          break;
        else {
          dist_old = array_distance( element_coord, point, work, ndim );
              // use new iso-parametric coordinates if distance decreases
          if ( iter>0 && dist_old>dist_tmp ) {
            step = step / 2.;
            dist_old = dist_tmp;
            array_move( iso_tmp, iso_old, ndim );
          }
          else
            step = 1.;
          converged = 1;
          for ( iiso=0; iiso<ndim; iiso++ ) {
              /* central differences to determine derivative of
                 distance with respect to iso parametric coordinates */
            array_move( iso_old, iso, MDIM ); 
            iso[iiso] = iso_old[iiso] + DELTA;
            interpolation_polynomial( iso[0], nnol_xi, h_xi, ddum );
            interpolation_polynomial( iso[1], nnol_eta, h_eta, ddum );
            interpolation_polynomial( iso[2], nnol_zeta, h_zeta, ddum );
            inol = 0; array_set( element_coord, 0., ndim );
            for ( inol_zeta=0; inol_zeta<nnol_zeta; inol_zeta++ ) {
              for ( inol_eta=0; inol_eta<nnol_eta; inol_eta++ ) {
                for ( inol_xi=0; inol_xi<nnol_xi; inol_xi++ ) {
                  weight[inol] = h_xi[inol_xi]*h_eta[inol_eta]*
                    h_zeta[inol_zeta];
                  for ( idim=0; idim<ndim; idim++ ) 
                    element_coord[idim] += weight[inol] * 
                    coords[inol*ndim+idim];
                  inol++;
                }
              }
            }
            dist_high = array_distance( element_coord, point, work, ndim );
            array_move( iso_old, iso, MDIM ); 
            iso[iiso] = iso_old[iiso] - DELTA;
            interpolation_polynomial( iso[0], nnol_xi, h_xi, ddum );
            interpolation_polynomial( iso[1], nnol_eta, h_eta, ddum );
            interpolation_polynomial( iso[2], nnol_zeta, h_zeta, ddum );
            inol = 0; array_set( element_coord, 0., ndim );
            for ( inol_zeta=0; inol_zeta<nnol_zeta; inol_zeta++ ) {
              for ( inol_eta=0; inol_eta<nnol_eta; inol_eta++ ) {
                for ( inol_xi=0; inol_xi<nnol_xi; inol_xi++ ) {
                  weight[inol] = h_xi[inol_xi]*h_eta[inol_eta]*
                    h_zeta[inol_zeta];
                  for ( idim=0; idim<ndim; idim++ ) 
                    element_coord[idim] += 
                    weight[inol] * coords[inol*ndim+idim];
                  inol++;
                }
              }
            }
            dist_low = array_distance( element_coord, point, work, ndim );
            iso_new[iiso] = iso_old[iiso];
            if ( scalar_dabs(dist_high-dist_low)>
                EPS_SIZE*element_largest_size ) {
              tmp = step*dist_old*2.*DELTA/(ndim*(dist_low-dist_high));
              iso_new[iiso] += tmp;
              if ( scalar_dabs(tmp)>EPS_SIZE ) converged = 0;
            }
          }
          dist_tmp = dist_old;
          array_move( iso_old, iso_tmp, MDIM );
          array_move( iso_new, iso_old, MDIM );
        }
      }
      if ( iso_old[0]<(-(1.+EPS_ISOP)) || iso_old[0]>(1.+EPS_ISOP) ||
           iso_old[1]<(-(1.+EPS_ISOP)) || iso_old[1]>(1.+EPS_ISOP) ||
           iso_old[2]<(-(1.+EPS_ISOP)) || iso_old[2]>(1.+EPS_ISOP) ||
           dist_old>element_largest_size*EPS_ISOP || !converged ) found = 0;
    }
  }

  return found;
}
