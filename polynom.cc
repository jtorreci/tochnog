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

#define EPS_TRIANGLE 1.e-6

void pol( long int element, long int element_group,
  long int name, long int nnol, double old_coord[], 
  double new_coord[], long int &npoint, double h[], 
  double old_d[], double new_d[], double new_b[], 
  double volume[] )

{
  long int idim=0, jdim=0, ixi=0, ieta=0, izeta=0, nxi=1, neta=1, nzeta=1,
    ipoint=0, npol=0, length=0, axisymmetric=-NO, npoint_per_dir=0,
    inol_xi=0, inol_eta=0, inol_zeta=0, nnol_xi=1, nnol_eta=1, nnol_zeta=1, 
    inol=0, istrain=0, indx=0, ldum=0,
    integration_method=-LOBATTO, integration_points=-MINIMAL;
  double detj=0., radius=0., fac=0., L1=0., L2=0., L3=0., L4=0.,
    ddum[1], iso[MPOINT], xi[MPOINT], 
    eta[MPOINT], zeta[MPOINT], weight_iso[MPOINT], 
    weight_xi[MPOINT], weight_eta[MPOINT], weight_zeta[MPOINT], 
    xj[MDIM*MDIM], h_xi[MNOL], p_xi[MNOL], 
    h_eta[MNOL], p_eta[MNOL], h_zeta[MNOL], p_zeta[MNOL], 
    p[MPOINT*MDIM*MNOL], 
    xj_inv[MPOINT*MDIM*MDIM], weight[MPOINT], coord_ip[MDIM];

  db( GROUP_AXISYMMETRIC, element_group, &axisymmetric, ddum, ldum, 
    VERSION_NORMAL, GET_IF_EXISTS );

  if ( name==-TRIA3 ) {
    nnol = 3; 
    db( GROUP_INTEGRATION_POINTS, element_group, &integration_points, ddum, 
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
    if ( axisymmetric==-YES ) integration_points = -MINIMAL;
    if ( integration_points==-NORMAL ) integration_points = -MINIMAL;
    if      ( integration_points==-MINIMAL ) {
      npoint = 1;
      weight[0] = 1.;
      array_set( h, 1./3., 3 );
      p[0] =  1.;
      p[1] =  0.;
      p[2] = -1.;
      p[3] =  0.;
      p[4] =  1.;
      p[5] = -1.;
    }
    else if ( integration_points==-MAXIMAL ) {
      npoint = nnol;
      for ( ipoint=0; ipoint<npoint; ipoint++ ) {
        weight[ipoint] = 1./npoint;
        p[ipoint*ndim*nnol+0] =  1.;
        p[ipoint*ndim*nnol+1] =  0.;
        p[ipoint*ndim*nnol+2] = -1.;
        p[ipoint*ndim*nnol+3] =  0.;
        p[ipoint*ndim*nnol+4] =  1.;
        p[ipoint*ndim*nnol+5] = -1.;
        for ( inol=0; inol<nnol; inol++ ) {
          if ( inol==ipoint )
            h[ipoint*nnol+inol] = 1.;
          else
            h[ipoint*nnol+inol] = 0.;
        }
      }
    }
    else
      db_error( GROUP_INTEGRATION_POINTS, element_group );
  }
  else if ( name==-TRIA6 ) {
    if ( axisymmetric==-YES ) {
      pri( "Error: not available for axisymmetric analysis: ", name );
      exit_tn_on_error();
    }
    db( GROUP_INTEGRATION_POINTS, element_group, &integration_points, ddum, 
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
    npoint = nnol;
    for ( ipoint=0; ipoint<npoint; ipoint++ ) {
      if      ( ipoint==0 ) {
        L1 = 1.;
        L2 = 0.;
        weight[ipoint] = EPS_TRIANGLE;
      }
      else if ( ipoint==1 ) {
        L1 = 0.;
        L2 = 1.;
        weight[ipoint] = EPS_TRIANGLE;
      }
      else if ( ipoint==2 ) {
        L1 = 0.;
        L2 = 0.;
        weight[ipoint] = EPS_TRIANGLE;
      }
      else if ( ipoint==3 ) {
        L1 = 0.5;
        L2 = 0.5;
        weight[ipoint] = 1./3.-EPS_TRIANGLE;
      }
      else if ( ipoint==4 ) {
        L1 = 0.0;
        L2 = 0.5;
        weight[ipoint] = 1./3.-EPS_TRIANGLE;
      }
      else if ( ipoint==5 ) {
        L1 = 0.5;
        L2 = 0.0;
        weight[ipoint] = 1./3.-EPS_TRIANGLE;
      }
      L3 = 1. - L1 - L2;
      h[ipoint*nnol+0] = (2.*L1-1.)*L1;
      h[ipoint*nnol+1] = 4.*L1*L2;
      h[ipoint*nnol+2] = (2.*L2-1.)*L2;
      h[ipoint*nnol+3] = 4.*L3*L1;
      h[ipoint*nnol+4] = 4.*L2*L3;
      h[ipoint*nnol+5] = (2.*L3-1.)*L3;
      p[ipoint*ndim*nnol+0*nnol+0] =  2.*L1 + (2.*L1-1.);
      p[ipoint*ndim*nnol+0*nnol+1] =  4.*L2;
      p[ipoint*ndim*nnol+0*nnol+2] =  0.;
      p[ipoint*ndim*nnol+0*nnol+3] =  4.*-1.*L1 + 4*L3;
      p[ipoint*ndim*nnol+0*nnol+4] =  4.*L2*-1.;
      p[ipoint*ndim*nnol+0*nnol+5] =  -2.*L3 + (2.*L3-1.)*-1.;
      p[ipoint*ndim*nnol+1*nnol+0] =  0.;
      p[ipoint*ndim*nnol+1*nnol+1] =  4.*L1;
      p[ipoint*ndim*nnol+1*nnol+2] =  2.*L2 + (2.*L2-1.);
      p[ipoint*ndim*nnol+1*nnol+3] =  4.*-1.*L1;
      p[ipoint*ndim*nnol+1*nnol+4] =  4.*L3 + 4.*L2*-1.;
      p[ipoint*ndim*nnol+1*nnol+5] =  -2.*L3 + (2.*L3-1.)*-1.;
    }
  }
  else if ( name==-PRISM6 ) {
    pri( "Error: PRISM6 is not available yet." );
    exit_tn_on_error();
  }
  else if ( name==-TET4 ) {
    db( GROUP_INTEGRATION_POINTS, element_group, &integration_points, ddum, 
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
    if ( axisymmetric==-YES ) integration_points = -MINIMAL;
    if ( integration_points==-NORMAL ) integration_points = -MINIMAL;
    nnol = 4; 
    if ( integration_points==-MINIMAL ) {
      npoint = 1; weight[0] = 1.;
      array_set( h, 1./4., 4 );
      p[0]  =  1.;
      p[1]  =  0.;
      p[2]  =  0.;
      p[3]  = -1.;
      p[4]  =  0.;
      p[5]  =  0.;
      p[6]  =  1.;
      p[7]  = -1.;
      p[8]  =  0.;
      p[9]  =  1.;
      p[10] =  0.;
      p[11] = -1.;
    }
    else if ( integration_points==-MAXIMAL ) {
      npoint = nnol;
      for ( ipoint=0; ipoint<npoint; ipoint++ ) {
        weight[ipoint] = 1./npoint;
        p[ipoint*ndim*nnol+0]  =  1.;
        p[ipoint*ndim*nnol+1]  =  0.;
        p[ipoint*ndim*nnol+2]  =  0.;
        p[ipoint*ndim*nnol+3]  = -1.;
        p[ipoint*ndim*nnol+4]  =  0.;
        p[ipoint*ndim*nnol+5]  =  0.;
        p[ipoint*ndim*nnol+6]  =  1.;
        p[ipoint*ndim*nnol+7]  = -1.;
        p[ipoint*ndim*nnol+8]  =  0.;
        p[ipoint*ndim*nnol+9]  =  1.;
        p[ipoint*ndim*nnol+10] =  0.;
        p[ipoint*ndim*nnol+11] = -1.;
        for ( inol=0; inol<nnol; inol++ ) {
          if ( inol==ipoint )
            h[ipoint*nnol+inol] = 1.;
          else
            h[ipoint*nnol+inol] = 0.;
        }
      }
    }
    else
      db_error( GROUP_INTEGRATION_POINTS, element_group );
  }
  else if ( name==-TET10 ) {
    db( GROUP_INTEGRATION_POINTS, element_group, &integration_points, ddum, 
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
    npoint = nnol;
    for ( ipoint=0; ipoint<npoint; ipoint++ ) {
      if      ( ipoint==0 ) {
        L1 = 1.;
        L2 = 0.;
        L3 = 0.;
        weight[ipoint] = EPS_TRIANGLE;
      }
      else if ( ipoint==1 ) {
        L1 = 0.;
        L2 = 1.;
        L3 = 0.;
        weight[ipoint] = EPS_TRIANGLE;
      }
      else if ( ipoint==2 ) {
        L1 = 0.;
        L2 = 0.;
        L3 = 1.;
        weight[ipoint] = EPS_TRIANGLE;
      }
      else if ( ipoint==3 ) {
        L1 = 0.;
        L2 = 0.;
        L3 = 0.;
        weight[ipoint] = EPS_TRIANGLE;
      }
      else if ( ipoint==4 ) {
        L1 = 0.5;
        L2 = 0.5;
        L3 = 0.0;
        weight[ipoint] = 1./6.-(4./6.)*EPS_TRIANGLE;
      }
      else if ( ipoint==5 ) {
        L1 = 0.5;
        L2 = 0.0;
        L3 = 0.5;
        weight[ipoint] = 1./6.-(4./6.)*EPS_TRIANGLE;
      }
      else if ( ipoint==6 ) {
        L1 = 0.0;
        L2 = 0.5;
        L3 = 0.5;
        weight[ipoint] = 1./6.-(4./6.)*EPS_TRIANGLE;
      }
      else if ( ipoint==7 ) {
        L1 = 0.5;
        L2 = 0.0;
        L3 = 0.0;
        weight[ipoint] = 1./6.-(4./6.)*EPS_TRIANGLE;
      }
      else if ( ipoint==8 ) {
        L1 = 0.0;
        L2 = 0.5;
        L3 = 0.0;
        weight[ipoint] = 1./6.-(4./6.)*EPS_TRIANGLE;
      }
      else if ( ipoint==9 ) {
        L1 = 0.0;
        L2 = 0.0;
        L3 = 0.5;
        weight[ipoint] = 1./6.-(4./6.)*EPS_TRIANGLE;
      }
      L4 = 1. - L1 - L2 - L3;

      h[ipoint*nnol+0] = (2.*L1-1.)*L1;
      h[ipoint*nnol+1] = 4.*L1*L2;
      h[ipoint*nnol+2] = (2.*L2-1.)*L2;
      h[ipoint*nnol+3] = 4.*L3*L1;
      h[ipoint*nnol+4] = 4.*L2*L3;
      h[ipoint*nnol+5] = (2.*L3-1.)*L3;
      h[ipoint*nnol+6] = 4.*L1*L4;
      h[ipoint*nnol+7] = 4.*L2*L4;
      h[ipoint*nnol+8] = 4.*L3*L4;
      h[ipoint*nnol+9] = (2.*L4-1.)*L4;

      p[ipoint*ndim*nnol+0*nnol+0] = 2.*L1 + (2.*L1-1.);
      p[ipoint*ndim*nnol+0*nnol+1] = 4.*L2;  
      p[ipoint*ndim*nnol+0*nnol+2] = 0.;
      p[ipoint*ndim*nnol+0*nnol+3] = 4.*L3;
      p[ipoint*ndim*nnol+0*nnol+4] = 0.;
      p[ipoint*ndim*nnol+0*nnol+5] = 0.;
      p[ipoint*ndim*nnol+0*nnol+6] = 4.*L4 + 4.*L1*-1.;
      p[ipoint*ndim*nnol+0*nnol+7] = 4.*L2*-1.;
      p[ipoint*ndim*nnol+0*nnol+8] = 4.*L3*-1.;
      p[ipoint*ndim*nnol+0*nnol+9] = 2.*-1.*L4 + (2.*L4-1.)*-1.;

      p[ipoint*ndim*nnol+1*nnol+0] = 0.;
      p[ipoint*ndim*nnol+1*nnol+1] = 4.*L1;
      p[ipoint*ndim*nnol+1*nnol+2] = 2.*L2 + (2.*L2-1.);
      p[ipoint*ndim*nnol+1*nnol+3] = 0.;
      p[ipoint*ndim*nnol+1*nnol+4] = 4.*L3;
      p[ipoint*ndim*nnol+1*nnol+5] = 0.;
      p[ipoint*ndim*nnol+1*nnol+6] = 4.*L1*-1.;
      p[ipoint*ndim*nnol+1*nnol+7] = 4.*L4 + 4.*L2*-1.;
      p[ipoint*ndim*nnol+1*nnol+8] = 4.*L3*-1.;
      p[ipoint*ndim*nnol+1*nnol+9] = 2.*-1.*L4 + (2.*L4-1.)*-1.;

      p[ipoint*ndim*nnol+2*nnol+0] = 0.;
      p[ipoint*ndim*nnol+2*nnol+1] = 0.;
      p[ipoint*ndim*nnol+2*nnol+2] = 0.;
      p[ipoint*ndim*nnol+2*nnol+3] = 4.*L1;
      p[ipoint*ndim*nnol+2*nnol+4] = 4.*L2;
      p[ipoint*ndim*nnol+2*nnol+5] = 2.*L3 + (2.*L3-1.);
      p[ipoint*ndim*nnol+2*nnol+6] = 4.*L1*-1.;
      p[ipoint*ndim*nnol+2*nnol+7] = 4.*L2*-1.;
      p[ipoint*ndim*nnol+2*nnol+8] = 4.*L4 + 4.*L3*-1.;
      p[ipoint*ndim*nnol+2*nnol+9] = 2.*-1.*L4 + (2.*L4-1.)*-1.;

    }
  }
  else {
    if      ( name==-BAR2 || name==-QUAD4  || name==-HEX8   ) npol = 2;
    else if ( name==-BAR3 || name==-QUAD9  || name==-HEX27  ) npol = 3;
    else if ( name==-BAR4 || name==-QUAD16 || name==-HEX64  ) npol = 4;
    else db_error( ELEMENT, element );
    if ( ndim>=1 ) nnol_xi = npol;
    if ( ndim>=2 ) nnol_eta = npol;
    if ( ndim==3 ) nnol_zeta = npol;
    nnol = nnol_xi * nnol_eta * nnol_zeta;
    if ( name==-BAR2 )
      integration_points = -MINIMAL;
    else
      integration_points = -MAXIMAL;
    if ( axisymmetric==-YES && materi_velocity ) {
      if ( npol>2 ) {
        pri( "Error, not available for axisymmetric analysis: ", name );
        exit_tn_on_error();
      }
      integration_method = -GAUSS;
      integration_points = -MINIMAL;
    }
    db( GROUP_INTEGRATION_POINTS, element_group, &integration_points, ddum, 
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( GROUP_INTEGRATION_METHOD, element_group, &integration_method, ddum, 
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
    if ( integration_points==-NORMAL ) integration_points = -MAXIMAL;
    npoint = 1;
    for ( idim=0; idim<MDIM; idim++ ) {
      if ( idim>ndim-1 )
        npoint_per_dir = 1;
      else if ( integration_points==-MINIMAL ) 
        npoint_per_dir = npol - 1;
      else if ( integration_points==-MAXIMAL ) 
        npoint_per_dir = npol;
      else
        db_error( GROUP_INTEGRATION_POINTS, element_group );
      npoint *= npoint_per_dir;
      if      ( integration_method==-GAUSS )
        integration_gauss( npoint_per_dir, iso, weight_iso );
      else if ( integration_method==-LOBATTO )
        integration_lobatto( npoint_per_dir, iso, weight_iso );
      else if ( npoint_per_dir<npol )
        integration_gauss( npoint_per_dir, iso, weight_iso );
      else
        integration_lobatto( npoint_per_dir, iso, weight_iso );
      if ( idim==0 ) {
        nxi = npoint_per_dir;
        array_move( iso, xi, nxi );
        array_move( weight_iso, weight_xi, nxi );
      }
      else if ( idim==1 ) {
        neta = npoint_per_dir;
        array_move( iso, eta, neta );
        array_move( weight_iso, weight_eta, neta );
      }
      else {
        assert( idim==2 );
        nzeta = npoint_per_dir;
        array_move( iso, zeta, nzeta );
        array_move( weight_iso, weight_zeta, nzeta );
      }
    }
    if ( npoint>MPOINT ) {
      cout << "\nError: too many int. points in element " << element << ".\n";
      exit(TN_EXIT_STATUS);
    }
    for ( izeta=0; izeta<nzeta; izeta++ ) {
      interpolation_polynomial( zeta[izeta], nnol_zeta, h_zeta, p_zeta );
      for ( ieta=0; ieta<neta; ieta++ ) {
        interpolation_polynomial( eta[ieta], nnol_eta, h_eta, p_eta );
        for ( ixi=0; ixi<nxi; ixi++ ) {
          interpolation_polynomial( xi[ixi], nnol_xi, h_xi, p_xi );
          inol = 0;
          ipoint = izeta*neta*nxi + ieta*nxi + ixi;
          weight[ipoint] = weight_xi[ixi]*weight_eta[ieta]*weight_zeta[izeta];
          for ( inol_zeta=0; inol_zeta<nnol_zeta; inol_zeta++ ) {
            for ( inol_eta=0; inol_eta<nnol_eta; inol_eta++ ) {
              for ( inol_xi=0; inol_xi<nnol_xi; inol_xi++ ) {
                h[ipoint*nnol+inol] = h_xi[inol_xi] * h_eta[inol_eta] *
                  h_zeta[inol_zeta];
                if ( ndim>=1 ) {
                  p[ipoint*ndim*nnol+0*nnol+inol] =
                    p_xi[inol_xi]*h_eta[inol_eta]*h_zeta[inol_zeta];
                }
                if ( ndim>=2 ) {
                  p[ipoint*ndim*nnol+1*nnol+inol] =
                    h_xi[inol_xi]*p_eta[inol_eta]*h_zeta[inol_zeta];
                }
                if ( ndim==3 ) {
                  p[ipoint*ndim*nnol+2*nnol+inol] =
                    h_xi[inol_xi]*h_eta[inol_eta]*p_zeta[inol_zeta];
                }
                inol++;
              }
            }
          }
        }
      }
    }
  }

  for ( ipoint=0; ipoint<npoint; ipoint++ ) {
    matrix_ab( &p[ipoint*ndim*nnol], old_coord, xj, ndim, nnol, ndim );
    if ( !matrix_inverse( xj, &xj_inv[ipoint*ndim*ndim], detj, ndim ) ) {
      detj = 0.;
      array_set( &xj_inv[ipoint*ndim*ndim], 0., ndim*ndim );
    }
    matrix_ab( &xj_inv[ipoint*ndim*ndim], &p[ipoint*ndim*nnol], 
      &old_d[ipoint*ndim*nnol], ndim, ndim, nnol );
    matrix_ab( &p[ipoint*ndim*nnol], new_coord, xj, ndim, nnol, ndim );
    if ( !matrix_inverse( xj, &xj_inv[ipoint*ndim*ndim], detj, ndim ) ) {
      detj = 0.;
      array_set( &xj_inv[ipoint*ndim*ndim], 0., ndim*ndim );
    }
    matrix_ab( &xj_inv[ipoint*ndim*ndim], &p[ipoint*ndim*nnol], 
      &new_d[ipoint*ndim*nnol], ndim, ndim, nnol );
    detj = scalar_dabs( detj );
    if      ( name==-TRIA3 ) volume[ipoint] = weight[ipoint]*detj/2.;
    else if ( name==-TRIA6 ) volume[ipoint] = weight[ipoint]*detj/2.;
    else if ( name==-TET4  ) volume[ipoint] = weight[ipoint]*detj/6.;
    else if ( name==-TET10 ) volume[ipoint] = weight[ipoint]*detj/6.;
    else if ( ndim==1  )     volume[ipoint] = weight[ipoint]*2.*detj;
    else if ( ndim==2 )      volume[ipoint] = weight[ipoint]*4.*detj;
    else                     volume[ipoint] = weight[ipoint]*8.*detj;
    if ( axisymmetric==-YES ) {
      matrix_ab( &h[ipoint*nnol], new_coord, coord_ip, 1, nnol, ndim );
      radius = scalar_dabs(coord_ip[0]);
      volume[ipoint] *= 2. * PIRAD * radius;
    }
  }

  if ( materi_velocity ) {
    array_set( new_b, 0., npoint*MSTRAIN*nnol*ndim );
    for ( ipoint=0; ipoint<npoint; ipoint++ ) {
      for ( inol=0; inol<nnol; inol++ ) {
        for ( istrain=0; istrain<MSTRAIN; istrain++ ) {
          if      ( istrain==0 ) {
            idim = 0;
            jdim = 0;
            fac  = 1.;
          }
          else if ( istrain==1 ) {
            idim = 0;
            jdim = 1;
            fac  = 2.;
          }
          else if ( istrain==2 ) {
            idim = 0;
            jdim = 2;
            fac  = 2.;
          }
          else if ( istrain==3 ) {
            idim = 1;
            jdim = 1;
            fac  = 1.;
          }
          else if ( istrain==4 ) {
            idim = 1;
            jdim = 2;
            fac  = 2.;
          }
          else {
            assert( istrain==5 );
            idim = 2;
            jdim = 2;
            fac  = 1.;
          }
          if ( jdim<ndim ) {
            indx = ipoint*MSTRAIN*nnol*ndim + istrain*nnol*ndim + inol*ndim + idim;
            new_b[indx] += 0.5 * fac * new_d[ipoint*ndim*nnol+jdim*nnol+inol];
            indx = ipoint*MSTRAIN*nnol*ndim + istrain*nnol*ndim + inol*ndim + jdim;
            new_b[indx] += 0.5 * fac * new_d[ipoint*ndim*nnol+idim*nnol+inol];
          }
          if ( axisymmetric==-YES && istrain==5 ) {
            matrix_ab( &h[ipoint*nnol], new_coord, coord_ip, 1, nnol, ndim );
            radius = scalar_dabs(coord_ip[0]);
            indx = ipoint*MSTRAIN*nnol*ndim + istrain*nnol*ndim + inol*ndim + 0;
            new_b[indx] += h[ipoint*nnol+inol] / radius;
          }
        }
      }
    }
  }

  length = db_len( ELEMENT, element, VERSION_NORMAL );
  if ( length!=nnol+1 ) db_error( ELEMENT, element );

}

void interpolation_polynomial( double iso, long int npol, double h_pol[],
  double p_pol[] )

{
  long int i=0, j=0, k=0;
  double tmp=0., iso_points[MPOINT], ddum[MPOINT];

  if      ( npol==1 ) {
    h_pol[0] = 1.;
    p_pol[0] = 0.;
  }
  else if ( npol==2 ) {
    h_pol[0] = 0.5*(1.-iso);
    h_pol[1] = 0.5*(1.+iso);
    p_pol[0] = -0.5;
    p_pol[1] = +0.5;
  }
  else {
    assert( integration_lobatto( npol, iso_points, ddum ) );
    for ( k=0; k<npol; k++ ) {
      h_pol[k] = 1.;
      p_pol[k] = 0.;
      for ( i=0; i<npol; i++ ) {
        if ( i!=k ) {
          h_pol[k] *= (iso-iso_points[i])/(iso_points[k]-iso_points[i]);
          tmp = 1.;
          for ( j=0; j<npol; j++ ) {
            if ( j==i )
              tmp *= 1./(iso_points[k]-iso_points[j]);
            else if ( j!=k )
              tmp *= (iso-iso_points[j])/(iso_points[k]-iso_points[j]);
          }
          p_pol[k] += tmp;
        }
      }
    }
  }
}
