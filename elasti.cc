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

#define EPS_LADE 1.e-12

void C_matrix( double young, double poisson,
  double transverse_isotropy[], double Ctot[MDIM][MDIM][MDIM][MDIM], long int task[] )

{
  long int idim=0, jdim=0, kdim=0, ldim=0, pdim=0, qdim=0, rdim=0, sdim=0,
    indx_Ctmp=0;
  double fac=0., Caaaa=0., Cbbbb=0., Caabb=0., Cabab=0., Cbcbc=0.,
   Cacac=0, Cbbcc=0, Ccccc=0, e[MDIM*MDIM], dir[MDIM], 
   dirtmp[MDIM*MDIM], C[MDIM][MDIM][MDIM][MDIM],
   Ctmp[MDIM*MDIM*MDIM*MDIM];

  array_set( &C[0][0][0][0], 0., MDIM*MDIM*MDIM*MDIM );
  if ( task[0]==GROUP_MATERI_ISOTROPY ) {
    if ( task[1]==-YES ) {
      if      ( ndim==1 )
        C[0][0][0][0] += young;
      else {
        assert( ndim==2 );
        fac = young/(1.-poisson*poisson);
        C[0][0][0][0] = C[1][1][1][1] += fac;
        C[0][0][1][1] = C[1][1][0][0] += fac*poisson;
        C[0][1][0][1] = C[1][0][1][0] += fac*(1.-poisson);
      }
    }
    else {
      fac = (young*(1.-poisson))/((1.+poisson)*(1.-2*poisson));
      C[0][0][0][0] = C[1][1][1][1] = C[2][2][2][2] += 
        fac;
      C[0][0][1][1] = C[0][0][2][2] = C[1][1][2][2] =
      C[1][1][0][0] = C[2][2][0][0] = C[2][2][1][1] += 
        (poisson/(1.-poisson))*fac;
      C[0][1][0][1] = C[0][2][0][2] = C[1][2][1][2] =
      C[1][0][1][0] = C[2][0][2][0] = C[2][1][2][1] +=
        ((1.-2.*poisson)/(1.-poisson))*fac;
    }
  }
  else if ( task[0]==GROUP_MATERI_ELASTI_TRANSVERSE_ISOTROPY_GRAHOUL ) {

      double factris = young/((1.+poisson)*(1.-2*poisson));	
      double alpha=transverse_isotropy[0];
      if ( ndim==2 ) {
	C[1][1][1][1] += factris*(1-poisson);
	C[0][0][0][0] = C[2][2][2][2] += factris*alpha*alpha*(1-poisson);
	C[1][1][0][0] = C[1][1][2][2] = C[0][0][1][1] = C[2][2][1][1] += alpha*poisson*factris;
	C[2][2][0][0] = C[0][0][2][2] += alpha*alpha*poisson*factris;
	C[1][0][1][0] = C[1][2][1][2] = C[0][1][0][1] = C[2][1][2][1] += alpha*(1.-2.*poisson)*factris;
	C[0][2][0][2] = C[2][0][2][0] += alpha*alpha*(1.-2.*poisson)*factris;
      }
      else if ( ndim==3 ) {
	C[2][2][2][2] += factris*(1-poisson);
	C[1][1][1][1] = C[0][0][0][0] += factris*alpha*alpha*(1-poisson);
	C[2][2][1][1] = C[2][2][0][0] = C[1][1][2][2] = C[0][0][2][2] += alpha*poisson*factris;
	C[0][0][1][1] = C[1][1][0][0] += alpha*alpha*poisson*factris;
	C[2][1][2][1] = C[2][0][2][0] = C[1][2][1][2] = C[0][2][0][2] += alpha*(1.-2.*poisson)*factris;
	C[1][0][1][0] = C[0][1][0][1] += alpha*alpha*(1.-2.*poisson)*factris;
      }
      else {
      	cout<<"Transverse isotropy grahoul not implemented for 1D"<<endl;
	exit_tn(-NO);
      }
  }
  else {
    assert( task[0]==GROUP_MATERI_ELASTI_TRANSVERSE_ISOTROPY );
    array_set( Ctmp, 0., MDIM*MDIM*MDIM*MDIM );
    array_move( transverse_isotropy, dir, MDIM ); 
    array_normalize( dir, MDIM );
    array_move( dir, dirtmp, MDIM );
    array_set( e, 0., MDIM*MDIM ); e[0] = 1.; e[4] = 1.; e[8] = 1.;
    array_outproduct_3D( &dirtmp[0], &e[0], &dirtmp[3] );
    if ( array_null( &dirtmp[3], MDIM ) ) 
      array_outproduct_3D( &dirtmp[0], &e[3], &dirtmp[3] );
    array_normalize( &dirtmp[3], MDIM );
    array_outproduct_3D( &dirtmp[0], &dirtmp[3], &dirtmp[6] );
    Caaaa = transverse_isotropy[MDIM+0];
    Cbbbb = transverse_isotropy[MDIM+1];
    Caabb = transverse_isotropy[MDIM+2];
    Cabab = transverse_isotropy[MDIM+3];
    Cbcbc = transverse_isotropy[MDIM+4];
    Cacac = Cabab;
    Cbbcc = Cbbbb - 2.*Cbcbc;
    Ccccc = Cbbbb;
    indx_Ctmp = 0*MDIM*MDIM*MDIM + 0*MDIM*MDIM + 0*MDIM + 0; 
    Ctmp[indx_Ctmp] = Caaaa;
    indx_Ctmp = 0*MDIM*MDIM*MDIM + 0*MDIM*MDIM + 1*MDIM + 1; 
    Ctmp[indx_Ctmp] = Caabb;
    indx_Ctmp = 1*MDIM*MDIM*MDIM + 1*MDIM*MDIM + 0*MDIM + 0; 
    Ctmp[indx_Ctmp] = Caabb;
    indx_Ctmp = 0*MDIM*MDIM*MDIM + 1*MDIM*MDIM + 0*MDIM + 1; 
    Ctmp[indx_Ctmp] = Cabab;
    indx_Ctmp = 1*MDIM*MDIM*MDIM + 0*MDIM*MDIM + 0*MDIM + 1; 
    Ctmp[indx_Ctmp] = Cabab;
    indx_Ctmp = 0*MDIM*MDIM*MDIM + 1*MDIM*MDIM + 1*MDIM + 0; 
    Ctmp[indx_Ctmp] = Cabab;
    indx_Ctmp = 1*MDIM*MDIM*MDIM + 0*MDIM*MDIM + 1*MDIM + 0; 
    Ctmp[indx_Ctmp] = Cabab;
    indx_Ctmp = 0*MDIM*MDIM*MDIM + 2*MDIM*MDIM + 0*MDIM + 2; 
    Ctmp[indx_Ctmp] = Cacac;
    indx_Ctmp = 2*MDIM*MDIM*MDIM + 0*MDIM*MDIM + 0*MDIM + 2; 
    Ctmp[indx_Ctmp] = Cacac;
    indx_Ctmp = 0*MDIM*MDIM*MDIM + 2*MDIM*MDIM + 2*MDIM + 0; 
    Ctmp[indx_Ctmp] = Cacac;
    indx_Ctmp = 2*MDIM*MDIM*MDIM + 0*MDIM*MDIM + 2*MDIM + 0; 
    Ctmp[indx_Ctmp] = Cacac;
    indx_Ctmp = 1*MDIM*MDIM*MDIM + 1*MDIM*MDIM + 1*MDIM + 1; 
    Ctmp[indx_Ctmp] = Cbbbb;
    indx_Ctmp = 1*MDIM*MDIM*MDIM + 1*MDIM*MDIM + 2*MDIM + 2; 
    Ctmp[indx_Ctmp] = Cbbcc;
    indx_Ctmp = 2*MDIM*MDIM*MDIM + 2*MDIM*MDIM + 1*MDIM + 1; 
    Ctmp[indx_Ctmp] = Cbbcc;
    indx_Ctmp = 1*MDIM*MDIM*MDIM + 2*MDIM*MDIM + 1*MDIM + 2; 
    Ctmp[indx_Ctmp] = Cbcbc;
    indx_Ctmp = 2*MDIM*MDIM*MDIM + 1*MDIM*MDIM + 1*MDIM + 2; 
    Ctmp[indx_Ctmp] = Cbcbc;
    indx_Ctmp = 1*MDIM*MDIM*MDIM + 2*MDIM*MDIM + 2*MDIM + 1; 
    Ctmp[indx_Ctmp] = Cbcbc;
    indx_Ctmp = 2*MDIM*MDIM*MDIM + 1*MDIM*MDIM + 2*MDIM + 1; 
    Ctmp[indx_Ctmp] = Cbcbc;
    indx_Ctmp = 2*MDIM*MDIM*MDIM + 2*MDIM*MDIM + 2*MDIM + 2; 
    Ctmp[indx_Ctmp] = Ccccc;
    for ( idim=0; idim<MDIM; idim++ ) {
      for ( jdim=0; jdim<MDIM; jdim++ ) {
        for ( kdim=0; kdim<MDIM; kdim++ ) {
          for ( ldim=0; ldim<MDIM; ldim++ ) {
            for ( pdim=0; pdim<MDIM; pdim++ ) {
              for ( qdim=0; qdim<MDIM; qdim++ ) {
                for ( rdim=0; rdim<MDIM; rdim++ ) {
                  for ( sdim=0; sdim<MDIM; sdim++ ) {
                    indx_Ctmp = pdim*MDIM*MDIM*MDIM + qdim*MDIM*MDIM +
                      rdim*MDIM + sdim;
                    C[idim][jdim][kdim][ldim] += dirtmp[pdim*MDIM+idim] * 
                      dirtmp[qdim*MDIM+jdim] * dirtmp[rdim*MDIM+kdim] * 
                      dirtmp[sdim*MDIM+ldim] * Ctmp[indx_Ctmp];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  array_add( &C[0][0][0][0], &Ctot[0][0][0][0], &Ctot[0][0][0][0], MDIM*MDIM*MDIM*MDIM );

}

void C_matrix_lade( double lade_1, double lade_2, double lade_3, 
  double I1, double sig_dev[], double Ctot[MDIM][MDIM][MDIM][MDIM] )

{
  long int i=0, j=0, k=0, l=0, m=0, indx1=0, indx2=0, nrot=0;
  double X=0., J2=0., B=0., B1=0., delta[MDIM][MDIM],
    compliance[MDIM][MDIM][MDIM][MDIM], C[MDIM][MDIM][MDIM][MDIM], val[MDIM*MDIM],
    dir[MDIM*MDIM*MDIM*MDIM], work[MDIM*MDIM][MDIM*MDIM];

  array_set( &C[0][0][0][0], 0., MDIM*MDIM*MDIM*MDIM );
  J2 = array_inproduct( sig_dev, sig_dev, MDIM*MDIM );
  X = (-I1)*(-I1) + lade_2*scalar_dabs(J2);
  if ( X<EPS_LADE  ) {
    return;
  }
  else {
    B = lade_1 * scalar_power( X, lade_3 );
    B1 = lade_1 * lade_3 * scalar_power( X, lade_3 -1 );
  }

  for ( l=0; l<MDIM; l++ ) {
    for ( m=0; m<MDIM; m++ ) {
      if ( l==m ) delta[l][m] = 1.;
      else delta[l][m] = 0.;
    }
  }

  for ( i=0; i<MDIM; i++ ) {
    for ( j=0; j<MDIM; j++ ) {
      for ( k=0; k<MDIM; k++ ) {
        for ( l=0; l<MDIM; l++ ) {
          compliance[i][j][k][l] =
            ( (1./3. - lade_2 ) * delta[i][j]*delta[k][l] / 3. + 
               lade_2 * delta[i][k] * delta[j][l] ) / B -
            2 * ( I1 * delta[i][j] / 3. + lade_2 * sig_dev[i*MDIM+j] ) *
                ( I1 * delta[k][l] / 3. + lade_2 * sig_dev[k*MDIM+l] ) * B1 / ( B*B );
        }
      }
    }
  }

    // invert compliance matrix to material stiffness matrix
  for ( i=0; i<MDIM; i++ ) {
    for ( j=0; j<MDIM; j++ ) {
      for ( k=0; k<MDIM; k++ ) {
        for ( l=0; l<MDIM; l++ ) {
          work[i*MDIM+j][k*MDIM+l] = compliance[i][j][k][l];
        }
      }
    }
  }
  matrix_jacobi( &work[0][0], MDIM*MDIM, val, dir, &nrot );
  for ( i=0; i<MDIM; i++ ) {
    for ( j=0; j<MDIM; j++ ) {
      for ( k=0; k<MDIM; k++ ) {
        for ( l=0; l<MDIM; l++ ) {
          for ( m=0; m<MDIM*MDIM; m++ ) {
            if ( val[m]!=0. ) {
              indx1 = i*MDIM +j;
              indx2 = k*MDIM +l;
              C[i][j][k][l] += 
                (1./val[m])*dir[indx1*MDIM*MDIM+m]*dir[indx2*MDIM*MDIM+m];
            }
          }
        }
      }
    }
  }
  array_add( &C[0][0][0][0], &Ctot[0][0][0][0], &Ctot[0][0][0][0], MDIM*MDIM*MDIM*MDIM );

}
