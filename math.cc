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

#define EPS_Q 1.e-10
#define EPS_DET 1.e-10

void array_add( double a[], double b[], double c[], long int n )

{
  register long int i=0;

  for ( i=0; i<n; i++ ) c[i] = a[i] + b[i];

}

double array_distance( double a[], double b[], double work[], long int n )

{
  array_subtract( a, b, work, n );
  return array_size( work, n );
}

double array_inproduct( double a[], double b[], long int n )

{
  register long int i=0;
  double result=0.;

  for ( i=0; i<n; i++ ) result += a[i]*b[i];
  return result;

}

long int array_member( long int list[], long int i, long int n, long int &indx )

{
  register long int j=0, found=0;

  indx = -1;

  while ( j<n && !found ) {
    if ( list[j]==i ) {
      indx = j;
      found = 1;
    }
    j++;
  }
  return found;

}

void array_move( long int from[], long int to[], long int n )

{
  long int i=0;

  for ( i=0; i<n; i++ ) to[i] = from[i];

}


void array_move( double from[], double to[], long int n )

{
  register long int i=0;

  for ( i=0; i<n; i++ ) to[i] = from[i];

}

void array_set( double *ptr, double value, long int n )

{
 register long int i=0;

 for ( i=0; i<n; i++ ) *(ptr+i) = value;
}

void array_set( long int *ptr, long int value, long int n )

{
 register long int i=0;

 for ( i=0; i<n; i++ ) *(ptr+i) = value;

}


void array_multiply( double a[], double b[], double c, long int n )

{
  register long int i=0;

  for ( i=0; i<n; i++ ) b[i] = c * a[i];

}

long int array_normalize( double a[], long int n )

{
  long int i=0;
  double l=0;

  l = array_size( a, n );
  if ( l<1.e-10 ) return 0;
  for ( i=0; i<n; i++ ) a[i] = a[i] / l;

  return 1;
}

long int array_null( double dval[], long int n )

{
  long int i=0;

  for ( i=0; i<n; i++ ) {
    if ( dval[i]!=0. ) return 0;
  }
  return 1;

}

void array_outproduct_2D( double a[], double b[] )
 
{
  double vec_tmp[MDIM], three[MDIM], normal[MDIM];

  array_move( a, vec_tmp, 2 ); vec_tmp[2] = 0.;
  array_set( three, 0., MDIM ); three[2] = 1.;
  array_outproduct_3D( three, vec_tmp, normal );          
  array_move( normal, b, 2 );
}


void array_outproduct_3D( double a[], double b[], double c[] )

{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

void array_subtract( double a[], double b[], double c[], long int n )

{
  register long int i=0;

  for ( i=0; i<n; i++ ) c[i] = a[i] - b[i];

}

double array_size( double a[], long int n )

{
  double size=0.;

  size = sqrt( scalar_dabs( array_inproduct(a,a,n) ) );

  return size;
}

long int equations_factor( long int n, double **a, long int *p, long int *f )

/* factorize matrix

   From Computers and Structures Vol. 52, No. 4, pp. 743-748, 1994
   O. Hededal and S. Krenk, A Profile Solver in C for Finite Element
   Equations.

   Adapted to catch large ratio of pivots (Ton van den Boogaard)

   Factors profile matrix [A] defined in the lower triangular
   part into off-diagonal factors [L] and the diagonal [D]:
                [A] = [L][D][L]T
   Rows and columns identified by f=1 are left untouched.

   input:
      n         = Dimension of the matrix a[][].
      a[i][j]   = Lower profile matrix to be factored.
                  i=0 to n-1, j=p[i] to i+1.
      p[i]      = Number of leading zeroes in row i.
      f[i]      = 0: unknown dof, known load.
                = 1: known dof, unknown load.

   output:
      a[i][j]   = i>j: Profile factor [L] below diagonal.
                = i=j: Diagonal factor [D] on diagonal.
      status    = No. of TINY diagonal elements
*/

{

   #define TINY 1E-20 /* to avoid numerical instability */
   #define max(x,y) (x>y ? x:y)

   long int i,j,k, status=0;
   double u=0., maxpiv=0., minpiv=0., abspiv=0.;

   for ( i=0; i<n; i++ ) if (!f[i])
   {
      for ( j=p[i]; j<i; j++ ) if (!f[j])
         for ( k=max(p[i],p[j]); k<j; k++ ) if (!f[k])
            a[i][j-p[i]] -= a[j][k-p[j]]*a[i][k-p[i]];

      for ( j=p[i]; j<i; j++ ) if (!f[j])
      {
         u = a[i][j-p[i]];
         a[i][j-p[i]] /= a[j][j-p[j]];
         a[i][i-p[i]] -= a[i][j-p[i]]*u;
      }

      abspiv = fabs(a[i][i-p[i]]);
      if ( abspiv < TINY )
      {
         a[i][i-p[i]] = TINY;
         status++;
      }

      if ( i==0 ) maxpiv = minpiv = abspiv;
      else if ( abspiv > maxpiv ) maxpiv = abspiv;
      else if ( abspiv < minpiv ) minpiv = abspiv;
   }

   if ( minpiv / maxpiv < 1.e-12 ) status = max( status, 1 );

   #undef max
   #undef TINY

   if ( status==0 ) status = 1;
   else status = 0;
   return status;
}

void equations_solve( long int n, double **a, long int *p, double *x,
  double *b, long int *f )

/*  solve:

   From Computers and Structures Vol. 52, No. 4, pp. 743-748, 1994
   O. Hededal and S. Krenk, A Profile Solver in C for Finite Element
   Equations.

   Adapted to catch large ratio of pivots (Ton van den Boogaard)

   Solves a symmetric system of equations, when the profile
   matrix [A] has been factored in the form
             [L][D][L]T {x} = {b}
   [L] is a lower profile matrix and [D] is a diagonal matrix.

   input:
     n        = Dimension of equation system.
     a[i][j]  = i>j: Off-diagonal factor [L].
              = i=j: Diagonal factor [D].
                i=0 to n-1, j=p[i]+1 to i.
     p[i]     = Number of leading zeroes in a row.
     x[i]     = System degrees of freedom.
     b[i]     = Load vector.
     f[i]     = 0: unknown dof, known load.
              = 1: known dof, unknown load.

   output:
     x[i]     = Degrees of freedom.
     b[i]     = Loads.
*/
{
   long int i,j;

   for ( i=0; i<n; i++ )
   {
      if (!f[i])
         x[i] = b[i];
      else
         b[i] = 0;
   }

   for ( i=0; i<n; i++ )
      for ( j=p[i]; j<i; j++ )
      {
         if (f[i] && !f[j])
            x[j] -= a[i][j-p[i]]*x[i];
         else if (f[j] && !f[i] )
            x[i] -= a[i][j-p[i]]*x[j];
      }

   for ( i=1; i<n; i++ ) if (!f[i])
      for ( j=p[i]; j<i; j++ ) if (!f[j])
         x[i] -= a[i][j-p[i]]*x[j];

   for ( i=0; i<n; i++ ) if (!f[i])
      x[i] /= a[i][i-p[i]];

   for ( i=n-1; i>0; i-- ) if (!f[i])
      for ( j=p[i]; j<i; j++ ) if (!f[j])
         x[j] -= a[i][j-p[i]]*x[i];

   for ( i=0; i<n; i++ )
      for ( j=p[i]; j<=i; j++ )
      {
         if (f[i])
            b[i] += a[i][j-p[i]]*x[j];
         if (f[j] && (i!=j) )
            b[j] += a[i][j-p[i]]*x[i];
      }

}

long int fit_polynomial( double points[], long int npoint, 
  double coefficients[], long int ncoefficient )

{
  long int icoeff=0, jcoeff=0, ipoint=0, return_value=0, length=0, n=0,
    *p=NULL, *f=NULL;
  double x=0., y=0., **a=NULL, *b=NULL;

  if ( npoint>0 && ncoefficient<=npoint ) {
    length = ncoefficient;
    a = new double * [length];
    for ( icoeff=0; icoeff<ncoefficient; icoeff++ ) {
      a[icoeff] = get_new_dbl(length);
      array_set( a[icoeff], 0., length );
    }
    b = get_new_dbl(length);
    array_set( b, 0., length );
    f = get_new_int(length);
    array_set( f, 0, length );
    p = get_new_int(length);
    array_set( p, 0, length );
    for ( ipoint=0; ipoint<npoint; ipoint++ ) {
      x = points[ipoint*2+0];
      y = points[ipoint*2+1];
      for ( icoeff=0; icoeff<ncoefficient; icoeff++ ) {
        b[icoeff] += scalar_power(x,icoeff)*y;
        for ( jcoeff=0; jcoeff<ncoefficient; jcoeff++ )
          a[icoeff][jcoeff] += scalar_power(x,icoeff+jcoeff);
      }
    }
    n = ncoefficient;
    if ( equations_factor( n, a, p, f ) ) {
      equations_solve( n, a, p, coefficients, b, f );    
      return_value = 1;
    }
    for ( icoeff=0; icoeff<ncoefficient; icoeff++ ) delete[] a[icoeff];
    delete[] a;
    delete[] b;
    delete[] f;
    delete[] p;
  }

  return return_value;
}

void itoa( int n, char str[] )

{
  int i=0, sign=0;

  if ( (sign=n)<0 ) n = -n;

  do {
    str[i++] = n % 10 + '0';
  } while ( (n/=10)>0 );

  if ( sign<0 ) str[i++] = '-';

  str[i] = '\0';

  string_reverse(str);
}

void matrix_ab( double *a, double *b, double *c, long int n, long int m,
  long int k )

  // c[n][k] = a[n][m] * b[m][k]

{
    register long int i=0, j=0, l=0;

  for ( i=0; i<n; i++ ) {
    for ( j=0; j<k; j++ ) {
      *( c + i*k + j ) = 0;
      for ( l=0; l<m; l++ ) {
        *( c + i*k + j ) += ( *( a + i*m + l ) ) *
                             ( *( b + l*k + j ) );
       }
    }                             

  }

}

void matrix_abat( double a[], double b[], double c[], 
  double work[], long int n )

{

  matrix_ab( a, b, work, n, n, n );
  matrix_abt( work, a, c, n, n, n );

}



void matrix_abt( double *a, double *b, double *c, long int n, long int m,
  long int k )

  // c[n][k] = a[n][m] * b[k][m]Transposed

{
  register long int i=0, j=0, l=0;

  for ( i=0; i<n; i++ ) {
    for ( j=0; j<k; j++ ) {
      *( c + i*k + j) = 0;
      for ( l=0; l<m; l++ ) {
        *( c + i*k + j) += ( *( a + i*m + l ) ) *
                           ( *( b + j*k + l ) );        
      }
    }
  }

}

void matrix_atb( double *a, double *b, double *c, long int n, long int m,
  long int k )

  // c[m][k] = a[n][m]Transposed * b[n][k]

{
  register long int i=0, j=0, l=0;

  for ( i=0; i<m; i++ ) {
    for ( j=0; j<k; j++ ) {
      *( c + i*k + j ) = 0;
      for ( l=0; l<n; l++ ) {
        *( c + i*k + j ) += ( *( a + l*m + i ) ) *
                            ( *( b + l*k + j ) );
      }
    }                                        
  }

}

void matrix_atba( double a[], double b[], double c[], 
  double work[], long int n, long int m )

  // c[m][m] = a[n][m]Transposed * b[n][n] * a[n][m]
{

  matrix_atb( a, b, work, n, m, n );
  matrix_ab( work, a, c, m, n, m );
}

void matrix_a4b( double a[3][3][3][3], double b[], double c[] )

{
  register long int i=0, j=0, k=0, l=0;

  for ( i=0; i<3; i++ ) {
    for ( j=0; j<3; j++ ) {
      c[i*3+j] = 0.;
      for ( k=0; k<3; k++ ) {
        for ( l=0; l<3; l++ ) {
          c[i*3+j] += a[i][j][k][l] * b[k*3+l];
        }
      }
    }
  }

}

double matrix_determinant( double a[], long int n )

{
  double result=0.;

  if ( n==1 )
    result = a[0];
  else if ( n==2 )
    result = a[0]*a[3] - a[1]*a[2];
  else {
    assert( n==3 );
    result = a[0]*(a[4]*a[8]-a[7]*a[5]) -
      a[1]*(a[3]*a[8]-a[6]*a[5]) + a[2]*(a[3]*a[7]-a[6]*a[4]);
  }
  return result;
}

void matrix_eigenvalues( double mat[], double eigenvalues[] )

{
  double I1=0., I2=0., I3=0., r=0., s=0., t=0., p=0., q=0.,
    bigR=0., phi=0., y0=0., y1=0., y2=0., tmp=0., inv[3];

  matrix_invariants( mat, inv );
  I1 = inv[0];
  I2 = inv[1];
  I3 = inv[2];
  r = -I1;
  s = +I2;
  t = -I3;
  p = (3.*s-r*r)/3.;
  q = 2.*r*r*r/27. - r*s/3. + t;
  if ( scalar_dabs(q)<EPS_Q ) {
    y0 = -sqrt(scalar_dabs(p));
    y1 = +sqrt(scalar_dabs(p));
    y2 = 0.;
  }
  else {
    bigR = sqrt(scalar_dabs(p)/3.); if ( q<0. ) bigR = -bigR;
    tmp = q/(2.*bigR*bigR*bigR);
    if ( tmp<-1. ) tmp = -1.;
    if ( tmp>+1. ) tmp = +1.;
    phi = acos(tmp);
    y0 = -2.*bigR*cos(phi/3.);
    y1 = -2.*bigR*cos(phi/3.+2.*PIRAD/3.);
    y2 = -2.*bigR*cos(phi/3.+4.*PIRAD/3.);
  }
  eigenvalues[0] = y0 - r/3.;
  eigenvalues[1] = y1 - r/3.;
  eigenvalues[2] = y2 - r/3.;
}

void matrix_insert( double a[], long int n, long int m,
  double b[], long int k, long int l, long int p )

  // add matrix a[n,m] into matrix b[*,p] at location k, l

{
  long int i=0, j=0, indxi=0, indxj=0;

  for ( i=0; i<n; i++ ) {
    for ( j=0; j<m; j++ ) {
      indxi = k + i;
      indxj = l + j;
      b[indxi*p+indxj] += a[i*m+j];
    }
  }
}

void matrix_invariants( double *mat, double *inv )

{

  inv[0] = mat[0] + mat[4] + mat[8];
  inv[1] = mat[0]*mat[4] + mat[4]*mat[8] + mat[8]*mat[0] -
           mat[1]*mat[3] - mat[5]*mat[7] - mat[6]*mat[2];
  inv[2] = matrix_determinant( mat, 3 );

}

long int matrix_inverse( double *mat, double *inv_mat, double &det, long int n )

{
  long int i=0;
  double inv_det=0., a1=0., a2=0., a3=0., norm=0.;

  for ( i=0; i<n*n; i++ ) norm += EPS_DET*mat[i]*mat[i];

  if ( n==1 ) {
    det = mat[0];
    if ( scalar_dabs(det)<=norm ) return 0;
    inv_mat[0] = 1. / mat[0];
  }
  else if ( n==2 ) {
    det = mat[0]*mat[3] - mat[1]*mat[2];
    if ( scalar_dabs(det)<=norm ) return 0;
    inv_det=1./det;
    inv_mat[0] =  mat[3]*inv_det;
    inv_mat[1] = -mat[1]*inv_det;
    inv_mat[2] = -mat[2]*inv_det;
    inv_mat[3] =  mat[0]*inv_det;
  }
  else {
    assert( n==3 );
    a1 = mat[4]*mat[8] - mat[7]*mat[5];
    a2 = mat[7]*mat[2] - mat[1]*mat[8];
    a3 = mat[1]*mat[5] - mat[4]*mat[2];
    det = mat[0]*a1+mat[3]*a2+mat[6]*a3;
    if ( scalar_dabs(det)<=norm ) return 0;
    inv_det = 1./det;
    inv_mat[0] = inv_det*a1;
    inv_mat[1] = inv_det*a2;
    inv_mat[2] = inv_det*a3;
    inv_mat[3] = inv_det*(mat[5]*mat[6]-mat[3]*mat[8]);
    inv_mat[4] = inv_det*(mat[0]*mat[8]-mat[6]*mat[2]);
    inv_mat[5] = inv_det*(mat[3]*mat[2]-mat[0]*mat[5]);
    inv_mat[6] = inv_det*(mat[3]*mat[7]-mat[6]*mat[4]);
    inv_mat[7] = inv_det*(mat[6]*mat[1]-mat[0]*mat[7]);
    inv_mat[8] = inv_det*(mat[0]*mat[4]-mat[3]*mat[1]);
  }
  return 1;

}

#define ROTATE(a,i,j,k,l) g=a[i*n+j];h=a[k*n+l];\
  a[i*n+j]=g-s*(h+g*tau);\
  a[k*n+l]=h+s*(g-h*tau);

void matrix_jacobi(double *a, long int n, double d[], double *v, long int *nrot)

   // destroys *a!

{
	long int j,iq,ip,i;
	double tresh,theta,tau,t,sm,s,h,g,c,b[100],z[100];

	if ( n>10 ) {
	  pri( "Program error: maximum of n exceeded in jacobi." );
	  exit(TN_EXIT_STATUS);
	}

	for (ip=0;ip<n;ip++) {
		for (iq=0;iq<n;iq++) v[ip*n+iq]=0.0;
		v[ip*n+ip]=1.0;
	}
	for (ip=0;ip<n;ip++) {
		b[ip]=d[ip]=a[ip*n+ip];
		z[ip]=0.0;
	}
	*nrot=0;
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++)
				sm += scalar_dabs(a[ip*n+iq]);
		}
		if (sm == 0.0) {
			return;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++) {
				g=100.0*scalar_dabs(a[ip*n+iq]);
				if (i > 4 && (double)(scalar_dabs(d[ip])+g) == (double)scalar_dabs(d[ip])
					&& (double)(scalar_dabs(d[iq])+g) == (double)scalar_dabs(d[iq]))
					a[ip*n+iq]=0.0;
				else if (scalar_dabs(a[ip*n+iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((double)(scalar_dabs(h)+g) == (double)scalar_dabs(h))
						t=(a[ip*n+iq])/h;
					else {
						theta=0.5*h/(a[ip*n+iq]);
						t=1.0/(scalar_dabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip*n+iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip*n+iq]=0.0;
					for (j=0;j<=ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=0;j<n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}
		for (ip=0;ip<n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
/*
	pri( "Error: max. number of iterations in JACOBI exceeded." );
	pri( "Use smaller time steps." );
	exit(TN_EXIT_STATUS);
*/
}

void matrix_symmetric( double a[], long int n )

{
  long int i=0, j=0, indx1=0, indx2=0;

  for ( i=0; i<n; i++ ) {
    for ( j=0; j<n; j++ ) {
      indx1 = i*n+j;
      indx2 = j*n+i;
      a[indx2] = a[indx1];
    }
  }
}

double scalar_dabs( double a ) 

{
  double result=0.;
  if ( a < 0. ) result = -a;
  else result = a;
  return result;
};

double scalar_dmax( double a, double b ) 

{
  double result=0;
  if ( a > b ) result = a;
  else result = b;
  return result;
};

double scalar_dmin( double a, double b ) 

{
  double result=0;
  if ( a < b ) result = a;
  else result = b;
  return result;
};

long int scalar_imax( long int a, long int b ) 

{  
  long int result=0;
  if ( a > b ) result = a;
  else result = b;
  return result;
};                                                 

double scalar_power( double a, double b )

{
  double result=0.;

  if ( b==0. )
    result = 1.;
  else
    result = pow( a, b );

  return result;
}

double scalar_ran_normal( int &idum )
  /* Return a normal distributed random number 
     with zero mean and unit variance.
     From: numerical recipes in c gasdev.
     Set idum to any negative value to initialize or
     reinitialize the sequence.*/
{
  static int iset=0;
  static double gset;
  double fac,r,v1,v2;

  if ( iset==0 ) {
    do {
      v1 = 2.0*scalar_ran_uniform(idum) - 1.0;
      v2 = 2.0*scalar_ran_uniform(idum) - 1.0;
      r = v1*v1+v2*v2;
    } while ( r>=1.0 || r==0.0 );  
    fac = sqrt( -2.0*log(r)/r );
    gset = v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return (double) gset;
  }
}


double scalar_ran_uniform( int &idum )
  /* Return a uniform random number between 0 and 1.
     From: numerical recipes in c */

#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243000
#define IA3 4561
#define IC3 51349

{
  static long ix1,ix2,ix3;
  static double r[98];
  double temp;
  static int iff=0;
  int j;
  
  if ( idum<0 || iff==0 ) {
    iff=1;
    ix1 = (IC1-(idum)) % M1;
    ix1 = (IA1*ix1+IC1) % M1;
    ix2 = ix1 % M2;
    ix1 = (IA1*ix1+IC1) % M1;
    ix3 = ix1 % M3;
    for ( j=1; j<=97; j++ ) {
      ix1 = (IA1*ix1+IC1) % M1;
      ix2 = (IA2*ix2+IC2) % M2;
      r[j] = (ix1+ix2*RM2)*RM1;
    }
    idum = 1;
  }
  ix1 = (IA1*ix1+IC1) % M1;
  ix2 = (IA2*ix2+IC2) % M2;
  ix3 = (IA3*ix3+IC3) % M3;
  j = (int) 1 + ((97*ix3)/M3);
  if ( j>97 || j<1 ) {
    pri( "Error in random generator." );
    exit(1);
  }
  temp = r[j];
  r[j] = (ix1+ix2*RM2)*RM1;
  return temp;
}  


double scalar_sign( double a )

{
  if ( a>0. )
    return 1.;
  else
    return -1.;
}

double scalar_square( double a )

{
  return a * a;
}

void sort( double val[], double vec[] )

{
  long int idim=0, min_indx=0, middle_indx=0, max_indx=0;
  double min_val=1.e20, max_val=-1.e20, work_val[MDIM], work_vec[MDIM*MDIM];

  array_move( val, work_val, MDIM );
  array_move( vec, work_vec, MDIM*MDIM );

  for ( idim=0; idim<MDIM; idim++ ) {
    if ( work_val[idim]<min_val ) {
      min_indx = idim;
      min_val = work_val[idim];
    }
    if ( work_val[idim]>max_val ) {
      max_indx = idim;
      max_val = work_val[idim];
    }
  } 
  assert( max_indx>=0 && max_indx<MDIM ); 
  assert( min_indx>=0 && min_indx<MDIM ); 
  if ( min_indx!=0 && max_indx!=0 )
    middle_indx = 0;
  else if ( min_indx!=1 && max_indx!=1 )
    middle_indx = 1;
  else
    middle_indx = 2;

  val[0] = work_val[max_indx];
  val[1] = work_val[middle_indx];
  val[2] = work_val[min_indx];
  array_move( &work_vec[max_indx*MDIM], &vec[0*MDIM], MDIM );
  array_move( &work_vec[middle_indx*MDIM], &vec[1*MDIM], MDIM );
  array_move( &work_vec[min_indx*MDIM], &vec[2*MDIM], MDIM );

}

void string_convert_to_lower_case( char str[] )

  // convert upper case to lower case

{
  int i=0, length=0;

  length = strlen(str);
  for ( i=0; i<length; i++ ) {
      // convert to lower case
    str[i] = tolower(str[i]);
  }
}

long int string_isinteger( char name[] )

  // test if string is an integer

{
  long int i=0, length=0, result=1;

  length = strlen(name);
  for ( i=0; i<length; i++ ) {
    if ( !isdigit(name[i]) ) {
      if ( name[i]!='-' && name[i]!='+' && name[i]!='e' ) result = 0;
    }
  }

  return result;
}


long int string_isdouble( char name[] )

  // test if string is a double

{
  long int i=0, length=0, result=1;

  length = strlen(name);
  for ( i=0; i<length; i++ ) {
    if ( !isdigit(name[i]) ) {
      if ( name[i]!='-' && name[i]!='+' &&
           name[i]!='e' && name[i]!='.' ) result = 0;
    }
  }

  return result;
}

void string_replace( char s[], char from, char to )

  // replace character in string

{
  long int i=0, l=0;

  l = strlen(s);
  for ( i=0; i<l; i++ ) {
    if ( s[i]==from ) s[i] = to;
  }

}


void string_reverse( char s[] )

  // reverse string

{
  long int i=0, l=0;
  char str[MCHAR];

  l = strlen(s);
  assert(l<MCHAR);

  for ( i=0; i<l; i++ ) str[i] = s[l-1-i];
  str[l] = '\0';
  strcpy(s,str);                   
}

void string_shorten( char s[], long int length )

  // shorten string

{
  long int l=0;

  l = strlen(s);

  if ( l>length ) {
    s[length] = '\0';
  }
}

long int table_xy( double table[], const char* table_name,
  long int length, double x, double &y )

{
  long int found=0, i=0, n=0, i2=0,i22=0,i21=0,i221=0;
  double x0=0., x1=0., y0=0., y1=0., dummy=0.;

  y = 0.;
  if      ( length==1 ) {
    y = table[0];
    return 1;
  }              
  else if ( length==2 ) {
    y = table[1];
    return 1;
  }              
  else if ( length>=4 ) {
    n = length / 2; found = 0; y = 0.;
      // order the table using simple sort  ofb
    for (i = 0;i < n-1; i++ ) {
      i2=i*2;i21=i2+1;i22=i21+1;i221=i22+1;
      if ( table[i2] > table[i22] ) { // swap
        dummy=table[i2];table[i2]=table[i22];table[i22]=dummy;  // swap xs
        dummy=table[i21];table[i21]=table[i221];table[i221]=dummy; // swap ys
       }
    }
      // end of sorting
    for ( i=0; !found && i<n-1; i++ ) {
      x0 = table[i*2+0];
      y0 = table[i*2+1];
      x1 = table[i*2+2];
      y1 = table[i*2+3];
      if ( x1<=x0 ) {
        pri( "Error found in ", table_name );
        pri( "The value  ", x1 );
        pri( "is not larger than ", x0 );
        pri( "Total table ", table, length );
        exit_tn_on_error();
      }
      if ( x>=(x0-1.e-10) && x<=x1 ) {
        found = 1;
        if ( x0==x1 )
          y = y0;
        else
          y = y0 + (y1-y0)*(x-x0)/(x1-x0);
      }
    }
  }

  return found;
}

long int table_xyz( double table[], long int number[], double coord[], double &z )

{
  long int inol=0, nnol=3, found=0, ix=0, iy=0, i0=0, i1=0, i2=0, i3=0, nx=0, ny=0;
  double coord0[MDIM], coord1[MDIM], coord2[MDIM], coord3[MDIM], zcoord[3], weight[3];

  if ( number[0]==1 && number[1]==1 ) {
    z = table[2];
    return 1;
  }
  
  nx = number[0];
  ny = number[1];
  for ( ix=0; !found && ix<nx-1; ix++ ) {
    for ( iy=0; !found && iy<ny-1; iy++ ) {
        // quad in grid
      i0 = (iy+0)*3*nx+((ix+0)*3);
      i1 = (iy+0)*3*nx+((ix+1)*3);
      i2 = (iy+1)*3*nx+((ix+0)*3);
      i3 = (iy+1)*3*nx+((ix+1)*3);
      coord0[0] = table[i0+0];
      coord0[1] = table[i0+1];
      coord0[2] = 0.;
      coord1[0] = table[i1+0];
      coord1[1] = table[i1+1];
      coord1[2] = 0.;
      coord2[0] = table[i2+0];
      coord2[1] = table[i2+1];
      coord2[2] = 0.;
      coord3[0] = table[i3+0];
      coord3[1] = table[i3+1];
      coord3[2] = 0.;
        // first triangle in quad
      project_point_on_triangle( coord, coord0, coord1, coord2, weight );
      found = 1;
      for ( inol=0; inol<nnol; inol++ ) {
        if ( weight[inol]<-EPS_ISO || weight[inol]>(1.+EPS_ISO) ) found = 0;
      }
      if ( found ) {
        z = 0.;
        zcoord[0] = table[i0+2];
        zcoord[1] = table[i1+2];
        zcoord[2] = table[i2+2];
        for ( inol=0; inol<nnol; inol++ ) z += weight[inol] * zcoord[inol];
      }
        // second triangle in quad
      if ( !found ) {
        project_point_on_triangle( coord, coord1, coord2, coord3, weight );
        found = 1;
        for ( inol=0; inol<nnol; inol++ ) {
          if ( weight[inol]<-EPS_ISO || weight[inol]>(1.+EPS_ISO) ) found = 0;
        }
        if ( found ) {
          z = 0.;
          zcoord[0] = table[i1+2];
          zcoord[1] = table[i2+2];
          zcoord[2] = table[i3+2];
          for ( inol=0; inol<3; inol++ ) z += weight[inol] * zcoord[inol];
        }
      }
    }
  }

  return found;
}

double triangle_area( double c0[], double c1[], double c2[] )

{
  double area=0., vec0[MDIM], vec1[MDIM], vec2[MDIM];

  array_set( vec0, 0., MDIM );
  array_set( vec1, 0., MDIM );
  array_set( vec2, 0., MDIM );

  array_subtract( c1, c0, vec0, ndim );
  array_subtract( c2, c0, vec1, ndim );
  array_outproduct_3D( vec0, vec1, vec2 ); 
  area = array_size( vec2, MDIM ) / 2.;
  return area;
}

double tetrahedron_volume( double c0[], double c1[], double c2[], double c3[] )

{
  long int nnol=4;
  double volume=0., detj=0., p[MDIM*MNOL], d[MDIM*MNOL], 
    coord[MNOL*MDIM], xj[MDIM*MDIM], xj_inv[MDIM*MDIM];

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

  array_move( c0, &coord[0*MDIM], MDIM );
  array_move( c1, &coord[1*MDIM], MDIM );
  array_move( c2, &coord[2*MDIM], MDIM );
  array_move( c3, &coord[3*MDIM], MDIM );

  matrix_ab( p, coord, xj, MDIM, nnol, MDIM );
  if ( matrix_inverse( xj, xj_inv, detj, MDIM ) ) {
    matrix_ab( xj_inv, p, d, MDIM, MDIM, nnol );
    detj = scalar_dabs( detj );
    volume = detj / 6.;
  }
  else volume = 0.;

  return volume;

}

/****************************************************************************************
	From Triax (D.Masin@city.ac.uk)
****************************************************************************************/

void matrix_inverse_general(double *matr, double *inv, int P) {
	int N=int(sqrt((double)P));
        register long int i=0, j=0 ;
	Matrix matice(N,N);
	Matrix img(N,N);
	Matrix imginv(N,N);
	Matrix invmatrix(N,N);
	for(i=1; i<=N; i++) {
		for(j=1; j<=N; j++) {
			matice(i,j)=matr[N*(i-1)+j-1];
		}
	}
	cinv( matice, img, invmatrix, imginv);
	for(i=1; i<=N; i++) {
		for(j=1; j<=N; j++) {
			inv[N*(i-1)+j-1]=invmatrix(i,j);
		}
	}
}

void make_dev(double tnz[9], double dev[9]) {
	array_set(dev, 0, 9);
	double tnzm = ( tnz[0] + tnz[4] + tnz[8] ) / 3.;    		
	array_move( tnz, dev, MDIM*MDIM );				
	for ( int idim=0; idim<MDIM; idim++ ) dev[idim*MDIM+idim] -= tnzm;
}


void matrix4_ab( double a[], double b[], double c[3][3][3][3] ) {
  register long int i=0, j=0, k=0, l=0;

  for ( i=0; i<3; i++ ) {
    for ( j=0; j<3; j++ ) {
      for ( k=0; k<3; k++ ) {
        for ( l=0; l<3; l++ ) {
          c[i][j][k][l] = a[i*3+j] * b[k*3+l];
        }
      }
    }
  }

}

void matrix_a_contr_b( double a[], double b[], double &c ) {
  register long int i=0, j=0 ;
  c=0;	

  for ( i=0; i<3; i++ ) {
    for ( j=0; j<3; j++ ) {
      c += a[i*3+j]*b[i*3+j];
    }
  }

}

void matrix_ab4( double a[], double b[3][3][3][3], double c[] ) {
  register long int i=0, j=0, k=0, l=0;

  for ( i=0; i<3; i++ ) {
    for ( j=0; j<3; j++ ) {
      c[i*3+j] = 0.;
      for ( k=0; k<3; k++ ) {
        for ( l=0; l<3; l++ ) {
          c[i*3+j] += a[k*3+l] * b[k][l][i][j];
        }
      }
    }
  }
}

void calc_IJlode(double geosigma[9], double &I, double &J, double &lode,
	bool calcderiv, double dIdsig[3][3], double dJdsig[3][3], double dlodedsig[3][3]) {
	
	I=geosigma[0]+geosigma[4]+geosigma[8];		//geosigma==positive in compression
	double sij[3][3];
	array_set(*sij, 0, 9);
	make_dev(geosigma, *sij);
	J=sqrt(array_inproduct(*sij, *sij, 9)/2);
	double kron_delta[MDIM][MDIM];
	array_set(*kron_delta, 0., MDIM*MDIM);
	kron_delta[0][0]=kron_delta[1][1]=kron_delta[2][2]=1;
	double S=0;
	double sijsjkski=0;
	for (int i=0; i<3; i++ ) {
		for (int j=0; j<3; j++ ) {
			for (int k=0; k<3; k++ ) {
		          sijsjkski += sij[i][j]*sij[j][k]*sij[k][i];
			}
		}
	}
	if(sijsjkski>=0) S=scalar_power(sijsjkski/3, 0.333333333333333333333333333333333333333333);
	else S=-scalar_power(-sijsjkski/3, 0.333333333333333333333333333333333333333333);
	if(J==0) J=TINY;
	double inlodebracket=3.0*sqrt(3.0)*scalar_power(S,3)/(2*scalar_power(J,3));
	if(inlodebracket>=1) inlodebracket=1;
	if(inlodebracket<=-1) inlodebracket=-1;
	lode=asin(inlodebracket)/3;
	if (calcderiv){
		array_move(*kron_delta, *dIdsig, MDIM*MDIM);
		for (int i=0; i<3; i++ ) {
			for (int j=0; j<3; j++ ) {
		          	dJdsig[i][j]= sij[i][j]/(2*J);
			}
		}
		double sikskj[3][3];
		array_set(*sikskj,0,9);
		for (int i=0; i<3; i++ ) {
			for (int j=0; j<3; j++ ) {
				for (int k=0; k<3; k++ ) {
		        	  sikskj[i][j] += sij[i][k]*sij[k][j];
				}
			}
		}
		double Sij[3][3];
		array_set(*Sij,0,9);
		array_set(*dlodedsig,0,9);
		if(lode!=PIRAD/6 && lode!=-PIRAD/6) {		//at these points grad lode == 0
			for (int i=0; i<3; i++ ) {
				for (int j=0; j<3; j++ ) {
			          	Sij[i][j]= sqrt(2.0)*(sqrt(3.0)*sikskj[i][j]/(2*J*J)-
						sin(3*lode)*sij[i][j]/(2*J)-
						kron_delta[i][j]/sqrt(3.0))/cos(3*lode);
				}
			}
			for (int i=0; i<3; i++ ) {
				for (int j=0; j<3; j++ ) {
			          	dlodedsig[i][j]=Sij[i][j]/(J*sqrt(2.0));
				}
			}
		}
	}
}

/****************************************************************************************
*********from http://www.algarcia.org/nummeth/Programs2E.html******************************
****************************************************************************************/

// Compute inverse of complex matrix
void cinv( Matrix RealA, Matrix ImagA, 
			 Matrix& RealAinv, Matrix& ImagAinv ) 
// Inputs
//   RealA  -    Real part of matrix A (N by N)
//   ImagA  -    Imaginary part of matrix A (N by N)
// Outputs
//   RealAinv  -    Real part of inverse of matrix A (N by N)
//   ImagAinv  -    Imaginary part of A inverse (N by N)
{

  int N = RealA.nRow();
  assert( N == RealA.nCol() && N == ImagA.nRow() 
	                        && N == ImagA.nCol());
    RealAinv = RealA; // Copy matrices to ensure they are same size
  ImagAinv = ImagA;
  
  int i, j, k;
  Matrix scale(N);	 // Scale factor
  int *index;  index = new int [N+1];

  //* Matrix B is initialized to the identity matrix
  Matrix RealB(N,N), ImagB(N,N);
  RealB.set(0.0);  ImagB.set(0.0);
  for( i=1; i<=N; i++ )
    RealB(i,i) = 1.0;

  //* Set scale factor, scale(i) = max( |a(i,j)| ), for each row
  for( i=1; i<=N; i++ ) {
    index[i] = i;			  // Initialize row index list
    double scaleMax = 0.;
    for( j=1; j<=N; j++ ) {
	  double MagA = RealA(i,j)*RealA(i,j) + ImagA(i,j)*ImagA(i,j);
      scaleMax = (scaleMax > MagA) ? scaleMax : MagA;
    }
	scale(i) = scaleMax;
  }

  //* Loop over rows k = 1, ..., (N-1)
  for( k=1; k<=N-1; k++ ) {
	//* Select pivot row from max( |a(j,k)/s(j)| )
    double ratiomax = 0.0;
	int jPivot = k;
    for( i=k; i<=N; i++ ) {
	  double MagA = RealA(index[i],k)*RealA(index[i],k) + 
		            ImagA(index[i],k)*ImagA(index[i],k);
      if(scale(index[i])==0) {
	cout<<endl<<" Error in matrix inversion ";
	exit(0);
      }
      double ratio = MagA/scale(index[i]);
      if( ratio > ratiomax ) {
        jPivot=i;
        ratiomax = ratio;
      }
    }
	//* Perform pivoting using row index list
	int indexJ = index[k];
	if( jPivot != k ) {	          // Pivot
      indexJ = index[jPivot];
      index[jPivot] = index[k];   // Swap index jPivot and k
      index[k] = indexJ;
	}
	//* Perform forward elimination
    for( i=k+1; i<=N; i++ ) {
	  double denom = RealA(indexJ,k)*RealA(indexJ,k) 
		           + ImagA(indexJ,k)*ImagA(indexJ,k);
	  if(denom==0) {
		cout<<endl<<" Error in matrix inversion ";
		exit(0);
	  }
      double RealCoeff = (RealA(index[i],k)*RealA(indexJ,k)
		               + ImagA(index[i],k)*ImagA(indexJ,k))/denom;
      double ImagCoeff = (ImagA(index[i],k)*RealA(indexJ,k)
		               - RealA(index[i],k)*ImagA(indexJ,k))/denom;
      for( j=k+1; j<=N; j++ ) {
        RealA(index[i],j) -= RealCoeff*RealA(indexJ,j)
		                   - ImagCoeff*ImagA(indexJ,j);
        ImagA(index[i],j) -= RealCoeff*ImagA(indexJ,j)
		                   + ImagCoeff*RealA(indexJ,j);
      }
	  RealA(index[i],k) = RealCoeff;
	  ImagA(index[i],k) = ImagCoeff;
      for( j=1; j<=N; j++ ) {
        RealB(index[i],j) -= RealA(index[i],k)*RealB(indexJ,j)
			               - ImagA(index[i],k)*ImagB(indexJ,j);
        ImagB(index[i],j) -= RealA(index[i],k)*ImagB(indexJ,j)
			               + ImagA(index[i],k)*RealB(indexJ,j);
	  }
    }
  }
  //* Perform backsubstitution
  for( k=1; k<=N; k++ ) {
	double denom = RealA(index[N],N)*RealA(index[N],N) 
		         + ImagA(index[N],N)*ImagA(index[N],N);
	if(denom==0) {
		cout<<endl<<" Error in matrix inversion ";
		exit(0);
	}
    RealAinv(N,k) = (RealB(index[N],k)*RealA(index[N],N) 
		          + ImagB(index[N],k)*ImagA(index[N],N))/denom;
    ImagAinv(N,k) = (ImagB(index[N],k)*RealA(index[N],N) 
		          - RealB(index[N],k)*ImagA(index[N],N))/denom;
    for( i=N-1; i>=1; i--) {
      double RealSum = RealB(index[i],k);
      double ImagSum = ImagB(index[i],k);
      for( j=i+1; j<=N; j++ ) {
        RealSum -= RealA(index[i],j)*RealAinv(j,k)
			     - ImagA(index[i],j)*ImagAinv(j,k);
        ImagSum -= RealA(index[i],j)*ImagAinv(j,k)
			     + ImagA(index[i],j)*RealAinv(j,k);
	  }
	  double denom = RealA(index[i],i)*RealA(index[i],i) 
		           + ImagA(index[i],i)*ImagA(index[i],i);
          if(denom==0) {
		cout<<endl<<" Error in matrix inversion ";
		exit(0);
	  }
      RealAinv(i,k) = (RealSum*RealA(index[i],i) 
		            + ImagSum*ImagA(index[i],i))/denom;
      ImagAinv(i,k) = (ImagSum*RealA(index[i],i) 
		            - RealSum*ImagA(index[i],i))/denom;
    }
  }

  delete [] index;	// Release allocated memory
}

