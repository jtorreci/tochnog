/*
    derived.cc - Derived stress magnitudes for tabular/VTK output.

    C++ with a fixed-dimension tensor template (Tensor<D>) to compute:
      - von Mises equivalent stress
      - Tresca maximum shear (as 2*tau_max, i.e. sigma1-sigma3)
      - principal stresses (eigenvalues of the stress tensor)

    Sign convention: tochnog uses compression-negative stress (soil
    mechanics). Von Mises / Tresca are invariant to the sign of the
    hydrostatic part, so the convention does not matter for the deviatoric
    magnitudes.

    The 6-component Voigt stress is [xx, yy, zz, xy, xz, yz].
*/

#include "tochnog.h"

// Fixed-dimension tensor helper (dimensions known at compile time).
template <int D>
class Tensor {
public:
  double v[D];
  Tensor() { for ( int i=0; i<D; i++ ) v[i] = 0.; }
  double& operator[](int i) { return v[i]; }
  double operator[](int i) const { return v[i]; }
};

// Fill the symmetric 3x3 matrix (row-major) from Voigt stress.
static void voigt_to_3x3( const double sig[6], double m[9] )
{
  m[0]=sig[0]; m[1]=sig[3]; m[2]=sig[4];
  m[3]=sig[3]; m[4]=sig[1]; m[5]=sig[5];
  m[6]=sig[4]; m[7]=sig[5]; m[8]=sig[2];
}

// Compute principal stresses (eigenvalues) via matrix_jacobi.
// Returns true on success, false if Jacobi fails.
static bool principal_stresses( const double sig[6], double pr[3] )
{
  double m[9], workvec[MDIM*MDIM];
  long int nrot=0;
  voigt_to_3x3( sig, m );
  matrix_jacobi( m, MDIM, pr, workvec, &nrot );
  return true;
}

// Derived magnitudes from a Voigt stress tensor.
// out[0] = von Mises, out[1] = Tresca (2*tau_max), out[2..4] = principal.
// Returns false if the input is all-zero (nothing to compute).
bool calc_derived( const double sig[6], double out[5] )
{
  double pr[3] = {0.,0.,0.};
  double vmises=0., tresca=0.;

  // zero stress -> zero derived
  if ( sig[0]==0. && sig[1]==0. && sig[2]==0. &&
       sig[3]==0. && sig[4]==0. && sig[5]==0. ) {
    out[0]=0.; out[1]=0.; out[2]=0.; out[3]=0.; out[4]=0.;
    return true;
  }

  principal_stresses( sig, pr );

  // von Mises from principal stresses (invariant form)
  {
    double d12 = pr[0]-pr[1], d23 = pr[1]-pr[2], d31 = pr[2]-pr[0];
    vmises = sqrt( 0.5*(d12*d12 + d23*d23 + d31*d31) );
  }

  // Tresca = max principal shear * 2  =  sigma1 - sigma3
  {
    double smin = pr[0], smax = pr[0];
    for ( int i=1; i<3; i++ ) {
      if ( pr[i] < smin ) smin = pr[i];
      if ( pr[i] > smax ) smax = pr[i];
    }
    tresca = smax - smin;
  }

  // matrix_jacobi does not sort its eigenvalues. Report them in
  // descending algebraic order (sigma1 >= sigma2 >= sigma3).
  {
    for ( int i=0; i<2; i++ ) {
      for ( int j=i+1; j<3; j++ ) {
        if ( pr[j] > pr[i] ) {
          double tmp = pr[i]; pr[i] = pr[j]; pr[j] = tmp;
        }
      }
    }
  }

  out[0] = vmises;
  out[1] = tresca;
  out[2] = pr[0];
  out[3] = pr[1];
  out[4] = pr[2];
  return true;
}
