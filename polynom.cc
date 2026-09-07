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
    // 6-noded prism (wedge): triangular base (nodes 1-3) at z=0 and
    // triangular top (nodes 4-6) at z=1 of the REFERENCE cell: the
    // standard triangle (area coords L1,L2,L3, area 1/2) times
    // zeta in [0,1] (volume 1/2). Shape functions: N_i = (1-z)*L_i
    // for the base, z*L_i for the top.
    //
    // CONVERGENCE (2026-09-07, interface-vs-wedge/tet4 bug, sprint
    // dev/iface-tri3d): the OLD rule of this branch was broken in two
    // ways. (1) The 6 integration weights were 0.25 each (sum 1.5 = 3x
    // the reference volume 0.5) AND the volume[] assembly below fell
    // through to the hex8 branch (weight*8*detj, the cube volume 8), so
    // every volume integral of the wedge came out 24x the physical one
    // (8*1.5/0.5): the stiffness, the mass, the gravity loads and the
    // nodal forces on the triangular faces (measured: a single wedge
    // with the exact field u=-z, E=1, sigma=-1 gives nodal reactions
    // +-4 instead of the consistent +-A/3 = +-1/6). Zero-thickness
    // prism6 interfaces against -prism6/-tet4 neighbours therefore
    // converged to sigma_iface = -24 instead of -1 (minimal model
    // zw.dat; interface11 el4 = -13.74 vs Pro -1). Hex8 neighbours
    // were exact because the hex8 branch is correct. (2) The old
    // z-levels (z=0.5 and z=0.75, both in the upper half) do NOT
    // integrate the zeta direction: the wedge of interface11 is an
    // OBLIQUE prism (the cut pieces connect a z=0.6 face to the z=1
    // face) whose Jacobian varies with zeta, and the mis-integrated
    // stiffness bends the field away from u=-z (interface11 node_dof
    // of the cut plane: -0.626..-0.758 instead of -0.6, per-pair
    // sigma_iface -0.32..-0.88 instead of uniform -1).
    //
    // New rule: 3-point degree-2 triangle rule (weights 1/6, sum 1/2
    // = the triangle area) x 2-point Gauss in zeta in [0,1] (weights
    // 1/2), i.e. 6 points of weight 1/12 (sum 1/2 = the reference
    // volume). This is the layout the Professional integrates (its
    // element_intpnt_coord of the zw.dbs prism elements shows zeta =
    // 0.2113248654 / 0.7886751346 = 1/2 +- 1/(2*sqrt(3))). With
    // sum(weight) = 1/2 the volume[] branch below integrates
    // weight*detj (no extra factor), exactly like the PRISM15 branch.
    db( GROUP_INTEGRATION_POINTS, element_group, &integration_points, ddum,
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
    nnol = 6;
    npoint = 6;   // 2 Gauss zeta-levels x 3 area points
    {
      // Gauss-Legendre 2 points mapped to zeta in [0,1]
      double zeta_pt[2], zeta_w[2];
      zeta_pt[0] = 0.5 - 0.5/sqrt(3.);   // 0.211324865405187
      zeta_pt[1] = 0.5 + 0.5/sqrt(3.);   // 0.788675134594813
      zeta_w[0] = 0.5; zeta_w[1] = 0.5;
      // degree-2 triangle rule: (2/3,1/6,1/6) permutations, weights 1/6
      double tri_L1[3] = { 2./3., 1./6., 1./6. };
      double tri_L2[3] = { 1./6., 2./3., 1./6. };
      double tri_w[3]  = { 1./6., 1./6., 1./6. };
      for ( ipoint=0; ipoint<npoint; ipoint++ ) {
        long int iz = ipoint/3, it = ipoint%3;
        double zz = zeta_pt[iz];
        L1 = tri_L1[it]; L2 = tri_L2[it]; L3 = 1. - L1 - L2;
        weight[ipoint] = zeta_w[iz]*tri_w[it]; // 1/12 each, sum 1/2
        double zbase = 1.-zz, ztop = zz;
        // shape functions at zz: (1-z)*L on the base, z*L on the top
        h[ipoint*nnol+0] = zbase*L1;
        h[ipoint*nnol+1] = zbase*L2;
        h[ipoint*nnol+2] = zbase*L3;
        h[ipoint*nnol+3] = ztop*L1;
        h[ipoint*nnol+4] = ztop*L2;
        h[ipoint*nnol+5] = ztop*L3;
        // derivatives d/d(L1) (x-direction of the reference triangle):
        // base d/dxi of (1-z)*L1 = (1-z), etc.; top = z
        p[ipoint*ndim*nnol+0*nnol+0] = zbase;
        p[ipoint*ndim*nnol+0*nnol+1] = 0.;
        p[ipoint*ndim*nnol+0*nnol+2] = -zbase;
        p[ipoint*ndim*nnol+0*nnol+3] = ztop;
        p[ipoint*ndim*nnol+0*nnol+4] = 0.;
        p[ipoint*ndim*nnol+0*nnol+5] = -ztop;
        // derivatives d/d(L2) (y-direction of the reference triangle)
        p[ipoint*ndim*nnol+1*nnol+0] = 0.;
        p[ipoint*ndim*nnol+1*nnol+1] = zbase;
        p[ipoint*ndim*nnol+1*nnol+2] = -zbase;
        p[ipoint*ndim*nnol+1*nnol+3] = 0.;
        p[ipoint*ndim*nnol+1*nnol+4] = ztop;
        p[ipoint*ndim*nnol+1*nnol+5] = -ztop;
        // derivatives d/d(zeta): d/dz of (1-z)*L = -L (base), z*L = +L
        p[ipoint*ndim*nnol+2*nnol+0] = -L1;
        p[ipoint*ndim*nnol+2*nnol+1] = -L2;
        p[ipoint*ndim*nnol+2*nnol+2] = -L3;
        p[ipoint*ndim*nnol+2*nnol+3] =  L1;
        p[ipoint*ndim*nnol+2*nnol+4] =  L2;
        p[ipoint*ndim*nnol+2*nnol+5] =  L3;
      }
    }
  }
  else if ( name==-PRISM15 ) {
    // 15-node serendipity quadratic prism (wedge). Node order (verified
    // against the Professional .dbs of the corpus prism15.dat):
    //   1-3  base triangle corners (area coords L1,L2,L3), zeta=-1
    //   4-6  top triangle corners (L1,L2,L3), zeta=+1
    //   7-9  mid nodes of the vertical edges above corners 1-3, zeta=0
    //   10-12 base triangle edge mids (1-2),(2-3),(3-1), zeta=-1
    //   13-15 top triangle edge mids (1-2),(2-3),(3-1), zeta=+1
    // with L3 = 1-L1-L2 and zeta in [-1,1] (z_phys = (1+zeta)/2 for a
    // unit-height prism). In-plane interpolation is the quadratic
    // 6-node triangle; along z the interpolation is quadratic Lagrange
    // on the vertical edges (serendipity: the zeta^2 vertical-mid term
    // corrects the corner/mid functions, there are no mid nodes on the
    // rectangular side faces - that would be the 18-node prism).
    //   corner base/top:  N = L_i (1-/+zeta) (2L_i - 2 -/+ zeta)/2
    //   vertical mid:     N = L_i (1-zeta^2)
    //   base edge mid:    N = 2 L_i L_j (1-zeta)   (i,j the edge corners)
    //   top edge mid:     N = 2 L_i L_j (1+zeta)
    // Integration rule (matches the Professional, measured from its
    // element_intpnt_coord): 3-point Gauss in zeta x the 7-point
    // degree-5 triangle rule (Dunavant) = 21 points. The shape
    // functions and rule reproduce the Professional element_intpnt_h
    // of prism15.dat exactly (dev-checked). The natural-domain measure
    // (triangle area 1/2 x zeta in [-1,1]) is 1, so volume[] uses the
    // plain weight*detj factor (see the volume chain below).
    nnol = 15;
    npoint = 21;
    {
      double sqrt15 = sqrt( 15. );
      double zeta_pt[3] = { -sqrt(0.6), 0., sqrt(0.6) };
      double zeta_w[3]  = { 5./9., 8./9., 5./9. };
      double a = (6.+sqrt15)/21., b = (9.-2.*sqrt15)/21.;
      double c = (6.-sqrt15)/21., d = (9.+2.*sqrt15)/21.;
      double tri_pt[7][2] = { {1./3.,1./3.}, {a,b}, {b,a}, {b,b},
        {c,d}, {d,c}, {d,d} };
      double tri_w[7] = { 9./80., (155.+sqrt15)/2400., (155.+sqrt15)/2400.,
        (155.+sqrt15)/2400., (155.-sqrt15)/2400., (155.-sqrt15)/2400.,
        (155.-sqrt15)/2400. };
      for ( ipoint=0; ipoint<npoint; ipoint++ ) {
        double zz = zeta_pt[ipoint/7];
        L1 = tri_pt[ipoint%7][0];
        L2 = tri_pt[ipoint%7][1];
        L3 = 1. - L1 - L2;
        weight[ipoint] = tri_w[ipoint%7]*zeta_w[ipoint/7];
        for ( inol=0; inol<nnol; inol++ )
          h[ipoint*nnol+inol] = 0.;
        array_set( &p[ipoint*ndim*nnol], 0., ndim*nnol );
        // base corners (inol 0-2): N = L_i (1-zz) (2L_i-2-zz) / 2
        if ( true ) {
          double Lk[3] = { L1, L2, L3 };
          for ( long int k=0; k<3; k++ ) {
            long int i = k;
            double f = 0.5*(1.-zz)*(2.*Lk[k]-2.-zz);
            double g = 0.5*(1.-zz)*(4.*Lk[k]-2.-zz);
            double dh = 0.5*Lk[k]*(1.-2.*Lk[k]+2.*zz);
            h[ipoint*nnol+i] = Lk[k]*f;
            if ( k==0 ) {
              p[ipoint*ndim*nnol+0*nnol+i] = g;
              p[ipoint*ndim*nnol+1*nnol+i] = 0.;
            }
            else if ( k==1 ) {
              p[ipoint*ndim*nnol+0*nnol+i] = 0.;
              p[ipoint*ndim*nnol+1*nnol+i] = g;
            }
            else {
              p[ipoint*ndim*nnol+0*nnol+i] = -g;
              p[ipoint*ndim*nnol+1*nnol+i] = -g;
            }
            p[ipoint*ndim*nnol+2*nnol+i] = dh;
          }
        }
        // top corners (inol 3-5): N = L_i (1+zz) (2L_i-2+zz) / 2
        if ( true ) {
          double Lk[3] = { L1, L2, L3 };
          for ( long int k=0; k<3; k++ ) {
            long int i = 3+k;
            double f = 0.5*(1.+zz)*(2.*Lk[k]-2.+zz);
            double g = 0.5*(1.+zz)*(4.*Lk[k]-2.+zz);
            double dh = 0.5*Lk[k]*(2.*Lk[k]-1.+2.*zz);
            h[ipoint*nnol+i] = Lk[k]*f;
            if ( k==0 ) {
              p[ipoint*ndim*nnol+0*nnol+i] = g;
              p[ipoint*ndim*nnol+1*nnol+i] = 0.;
            }
            else if ( k==1 ) {
              p[ipoint*ndim*nnol+0*nnol+i] = 0.;
              p[ipoint*ndim*nnol+1*nnol+i] = g;
            }
            else {
              p[ipoint*ndim*nnol+0*nnol+i] = -g;
              p[ipoint*ndim*nnol+1*nnol+i] = -g;
            }
            p[ipoint*ndim*nnol+2*nnol+i] = dh;
          }
        }
        // vertical edge mids (inol 6-8): N = L_i (1-zz^2)
        if ( true ) {
          double Lk[3] = { L1, L2, L3 };
          for ( long int k=0; k<3; k++ ) {
            long int i = 6+k;
            h[ipoint*nnol+i] = Lk[k]*(1.-zz*zz);
            if ( k==0 ) {
              p[ipoint*ndim*nnol+0*nnol+i] = 1.-zz*zz;
              p[ipoint*ndim*nnol+1*nnol+i] = 0.;
            }
            else if ( k==1 ) {
              p[ipoint*ndim*nnol+0*nnol+i] = 0.;
              p[ipoint*ndim*nnol+1*nnol+i] = 1.-zz*zz;
            }
            else {
              p[ipoint*ndim*nnol+0*nnol+i] = -(1.-zz*zz);
              p[ipoint*ndim*nnol+1*nnol+i] = -(1.-zz*zz);
            }
            p[ipoint*ndim*nnol+2*nnol+i] = -2.*zz*Lk[k];
          }
        }
        // base triangle edge mids (inol 9-11): edges (1,2),(2,3),(3,1)
        h[ipoint*nnol+9] = 2.*L1*L2*(1.-zz);
        p[ipoint*ndim*nnol+0*nnol+9] = 2.*L2*(1.-zz);
        p[ipoint*ndim*nnol+1*nnol+9] = 2.*L1*(1.-zz);
        p[ipoint*ndim*nnol+2*nnol+9] = -2.*L1*L2;
        h[ipoint*nnol+10] = 2.*L2*L3*(1.-zz);
        p[ipoint*ndim*nnol+0*nnol+10] = -2.*L2*(1.-zz);
        p[ipoint*ndim*nnol+1*nnol+10] = 2.*(L3-L2)*(1.-zz);
        p[ipoint*ndim*nnol+2*nnol+10] = -2.*L2*L3;
        h[ipoint*nnol+11] = 2.*L3*L1*(1.-zz);
        p[ipoint*ndim*nnol+0*nnol+11] = 2.*(L3-L1)*(1.-zz);
        p[ipoint*ndim*nnol+1*nnol+11] = -2.*L1*(1.-zz);
        p[ipoint*ndim*nnol+2*nnol+11] = -2.*L3*L1;
        // top triangle edge mids (inol 12-14): edges (1,2),(2,3),(3,1)
        h[ipoint*nnol+12] = 2.*L1*L2*(1.+zz);
        p[ipoint*ndim*nnol+0*nnol+12] = 2.*L2*(1.+zz);
        p[ipoint*ndim*nnol+1*nnol+12] = 2.*L1*(1.+zz);
        p[ipoint*ndim*nnol+2*nnol+12] = 2.*L1*L2;
        h[ipoint*nnol+13] = 2.*L2*L3*(1.+zz);
        p[ipoint*ndim*nnol+0*nnol+13] = -2.*L2*(1.+zz);
        p[ipoint*ndim*nnol+1*nnol+13] = 2.*(L3-L2)*(1.+zz);
        p[ipoint*ndim*nnol+2*nnol+13] = 2.*L2*L3;
        h[ipoint*nnol+14] = 2.*L3*L1*(1.+zz);
        p[ipoint*ndim*nnol+0*nnol+14] = 2.*(L3-L1)*(1.+zz);
        p[ipoint*ndim*nnol+1*nnol+14] = -2.*L1*(1.+zz);
        p[ipoint*ndim*nnol+2*nnol+14] = 2.*L3*L1;
      }
    }
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
    // group_element_selective_reduced_integration (SRI, Hughes): for the
    // bilinear quad4 and the trilinear hex8 with the keyword active the
    // FULL rule becomes the classic Gauss rule (2x2 / 2x2x2; the
    // codebase default is the Lobatto rule at the element corners,
    // which over-integrates the quadratic bending energy and caps a
    // shear-only SRI at ~0.35x of the exact section moment; measured).
    // The shear terms are handled in materi(): D is split into D_norm
    // (integrated here with this full rule) and D_shear (integrated
    // with 1 Gauss point at the centroid).
    if ( sri_active( element, element_group, name, nnol ) )
      integration_method = -GAUSS;
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
    else if ( name==-PRISM15 ) volume[ipoint] = weight[ipoint]*detj;
    else if ( name==-PRISM6 ) {
      // CONVERGENCE (2026-09-07, interface-vs-wedge/tet4 bug, sprint
      // dev/iface-tri3d): the wedge FELL THROUGH to the hex8 branch
      // below (weight*8*detj), but its reference cell is the standard
      // triangle (area 1/2) x zeta in [0,1] (volume 1/2), NOT the hex8
      // cube of volume 8; the old PRISM6 branch weights also summed to
      // 1.5 = 3x the reference volume. The assembled element volume
      // (and with it the stiffness, gravity loads, mass and every
      // nodal force derived from volume[]) came out 24x the physical
      // one (8*1.5/0.5 = 24), so the nodal force a compressed wedge
      // exerts on its triangular faces was 24x the consistent load
      // (measured: u=-z, E=1, sigma=-1 -> reactions +-4 per node
      // instead of +-A/3 = +-1/6) - the exact factor seen in the
      // zero-thickness prism6 interface minimal model (sigma_iface
      // -24 vs Professional -1) and in interface11 (sigma_el4 -13.74
      // vs -1). Hex8-neighbour interfaces were exact because the hex8
      // integration is correct. The rewritten PRISM6 branch above uses
      // weights summing to 1/2 (the reference volume), so here the
      // integral is the plain weight*detj (PRISM15 pattern).
      volume[ipoint] = weight[ipoint]*detj;
    }
    else if ( ndim==1  )     volume[ipoint] = weight[ipoint]*2.*detj;
    else if ( ndim==2 )      volume[ipoint] = weight[ipoint]*4.*detj;
    else                     volume[ipoint] = weight[ipoint]*8.*detj;
    if ( axisymmetric==-YES ) {
      matrix_ab( &h[ipoint*nnol], new_coord, coord_ip, 1, nnol, ndim );
      radius = scalar_dabs(coord_ip[0]);
      volume[ipoint] *= 2. * PIRAD * radius;
    }
  }

  // check_element_shape: warn if element is too distorted
  if ( check_element_shape_factor>0. ) {
    double vol_avg=0., distortion=0.;
    long int ip2=0;
    for ( ip2=0; ip2<npoint; ip2++ ) vol_avg += scalar_dabs(volume[ip2]);
    vol_avg /= npoint;
    if ( vol_avg>0. ) {
      for ( ip2=0; ip2<npoint; ip2++ )
        distortion += scalar_dabs(volume[ip2]-vol_avg)/vol_avg;
      distortion /= npoint;
      if ( distortion>check_element_shape_factor ) {
        char warn_str[MCHAR];
        snprintf( warn_str, MCHAR,
          "Warning: element with distorted shape, distortion=%g > factor=%g",
          distortion, check_element_shape_factor );
        pri( warn_str );
      }
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
