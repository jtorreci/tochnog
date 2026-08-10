/*
    masin_visco.c - C port of umat_visco.f (and umat_cleanvisco.f, identical)
    Visco-hypoplasticity of clays (Jerman & Masin 2020), the rate-dependent
    extension of the Masin clay model. Built on the same RKF23 integration
    skeleton as masin.c but with the LD (loading-direction) formulation:
    parameters ocparam, beta, ksi, gama and Dref add the visco-plastic
    scaling and the shear-band rotation of the flow direction.

    Original Fortran: Jerman & Masin (GPL). This file is a faithful C port of
    the reference implementation umat_visco.f (3181 lines). Function names map
    1:1 to the Fortran subroutines; indices are converted from Fortran 1-based
    to C 0-based.

    Conventions (identical to the Fortran):
      - Voigt notation: [11,22,33,12,13,23]
      - tension/extension positive
      - dot_vect flag: 1=stress-like (shear weight 2), 2=strain-like (shear
        weight 0.5), 3=ordinary
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NASV    10
#define NYDIM   (6 + NASV)

static const double PI = 3.14159265358979323846;

/* ------------------------------------------------------------------ */
/* utility: is the number NaN / huge?                                  */
/* ------------------------------------------------------------------ */
static void umatisnan(double chcknum, int *testnan)
{
  if (!(chcknum >= 0.0 || chcknum < 0.0)) *testnan = 1;
  if (chcknum >  1.0e30) *testnan = 1;
  if (chcknum < -1.0e30) *testnan = 1;
  if (chcknum != chcknum) *testnan = 1;
}

/* ------------------------------------------------------------------ */
/* dot product of Voigt vectors                                        */
/* ------------------------------------------------------------------ */
static double dot_vect(int flag, const double *a, const double *b, int n)
{
  int i;
  double coeff, res = 0.0;
  if (flag == 1) coeff = 2.0;        /* stress-like */
  else if (flag == 2) coeff = 0.5;   /* strain-like */
  else coeff = 1.0;                  /* ordinary */
  for (i = 0; i < n; i++) {
    if (i < 3) res += a[i] * b[i];
    else res += coeff * a[i] * b[i];
  }
  return res;
}

/* ------------------------------------------------------------------ */
/* matrix multiplication c = a * b  (a lxm, b mxn, c lxn)              */
/* row-major, flat arrays                                              */
/* ------------------------------------------------------------------ */
static void matmul(const double *a, const double *b, double *c,
  int l, int m, int n)
{
  int i, j, k;
  for (i = 0; i < l; i++)
    for (j = 0; j < n; j++) {
      c[i * n + j] = 0.0;
      for (k = 0; k < m; k++)
        c[i * n + j] += a[i * m + k] * b[k * n + j];
    }
}

/* ------------------------------------------------------------------ */
/* invariants of a strain tensor (Voigt)                               */
/* ------------------------------------------------------------------ */
static void inv_eps(double *eps, double *eps_v, double *eps_s, double *sin3t)
{
  int i;
  double edev[6], edev2[6], ev3, norm2, numer, denom, tredev3;
  double onethird = 1.0 / 3.0, twothirds = 2.0 / 3.0, sqrt6 = sqrt(6.0);

  *eps_v = eps[0] + eps[1] + eps[2];
  ev3 = onethird * (*eps_v);

  edev[0] = eps[0] - ev3;
  edev[1] = eps[1] - ev3;
  edev[2] = eps[2] - ev3;
  edev[3] = eps[3] / 2.0;
  edev[4] = eps[4] / 2.0;
  edev[5] = eps[5] / 2.0;

  norm2 = edev[0]*edev[0] + edev[1]*edev[1] + edev[2]*edev[2] +
    2.0 * (edev[3]*edev[3] + edev[4]*edev[4] + edev[5]*edev[5]);
  *eps_s = sqrt(twothirds * norm2);

  edev2[0] = edev[0]*edev[0] + edev[3]*edev[3] + edev[4]*edev[4];
  edev2[1] = edev[3]*edev[3] + edev[1]*edev[1] + edev[5]*edev[5];
  edev2[2] = edev[5]*edev[5] + edev[4]*edev[4] + edev[2]*edev[2];
  edev2[3] = 2.0*(edev[0]*edev[3] + edev[3]*edev[1] + edev[5]*edev[4]);
  edev2[4] = 2.0*(edev[4]*edev[0] + edev[5]*edev[3] + edev[2]*edev[4]);
  edev2[5] = 2.0*(edev[3]*edev[4] + edev[1]*edev[5] + edev[5]*edev[2]);

  if (*eps_s == 0.0) {
    *sin3t = -1.0;
  } else {
    tredev3 = 0.0;
    for (i = 0; i < 6; i++) tredev3 += edev[i] * edev2[i];
    numer = sqrt6 * tredev3;
    denom = pow(sqrt(norm2), 3.0);
    *sin3t = numer / denom;
    if (fabs(*sin3t) > 1.0) *sin3t = *sin3t / fabs(*sin3t);
  }
}

/* ------------------------------------------------------------------ */
/* invariants of a stress tensor (Voigt)                               */
/* ------------------------------------------------------------------ */
static void inv_sig(double *sig, double *pp, double *qq, double *cos3t,
  double *cos3t_rot, double *I1, double *I2, double *I3,
  double *I1rot, double *I2rot, double *I3rot, double *sig_rot,
  double *parms)
{
  double sdev[6], eta[6], eta_d[6], eta_d2[6];
  double xmin1, xmin2, xmin3, tretadev3;
  double norm2, norm2sig, norm2eta, numer, denom;
  double half = 0.5, one = 1.0, three = 3.0;
  double onethird = 1.0 / 3.0, threehalves = 3.0 / 2.0, sqrt6 = sqrt(6.0);
  double tiny = 1.0e-18;
  double vertical, beta, beta_deg, trace_sig, fullbeta[6];
  double numer_rot, denom_rot, norm2sig_rot, norm2eta_rot;
  double eta_rot[6], eta_d_rot[6], eta_d2_rot[6];
  double tretadev3_rot;
  double xmin1_rot, xmin2_rot, xmin3_rot;
  int i;

  *I1 = sig[0] + sig[1] + sig[2];
  *pp = onethird * (*I1);

  sdev[0] = sig[0] - *pp;
  sdev[1] = sig[1] - *pp;
  sdev[2] = sig[2] - *pp;
  sdev[3] = sig[3];
  sdev[4] = sig[4];
  sdev[5] = sig[5];

  if (*I1 != 0.0) {
    for (i = 0; i < 6; i++) eta[i] = sig[i] / (*I1);
  } else {
    for (i = 0; i < 6; i++) eta[i] = sig[i] / tiny;
  }
  eta_d[0] = eta[0] - onethird;
  eta_d[1] = eta[1] - onethird;
  eta_d[2] = eta[2] - onethird;
  eta_d[3] = eta[3];
  eta_d[4] = eta[4];
  eta_d[5] = eta[5];

  norm2 = dot_vect(1, sdev, sdev, 6);
  norm2sig = dot_vect(1, sig, sig, 6);
  norm2eta = dot_vect(1, eta_d, eta_d, 6);

  *qq = sqrt(threehalves * norm2);
  *I2 = half * (norm2sig - (*I1) * (*I1));

  eta_d2[0] = eta_d[0]*eta_d[0] + eta_d[3]*eta_d[3] + eta_d[4]*eta_d[4];
  eta_d2[1] = eta_d[3]*eta_d[3] + eta_d[1]*eta_d[1] + eta_d[5]*eta_d[5];
  eta_d2[2] = eta_d[5]*eta_d[5] + eta_d[4]*eta_d[4] + eta_d[2]*eta_d[2];
  eta_d2[3] = eta_d[0]*eta_d[3] + eta_d[3]*eta_d[1] + eta_d[5]*eta_d[4];
  eta_d2[4] = eta_d[4]*eta_d[0] + eta_d[5]*eta_d[3] + eta_d[2]*eta_d[4];
  eta_d2[5] = eta_d[3]*eta_d[4] + eta_d[1]*eta_d[5] + eta_d[5]*eta_d[2];

  if (norm2eta < tiny) {
    *cos3t = -one;
  } else {
    tretadev3 = dot_vect(1, eta_d, eta_d2, 6);
    numer = -sqrt6 * tretadev3;
    denom = pow(sqrt(norm2eta), 3.0);
    *cos3t = numer / denom;
    if (fabs(*cos3t) > one) *cos3t = *cos3t / fabs(*cos3t);
  }

  xmin1 = sig[1]*sig[2] - sig[5]*sig[5];
  xmin2 = sig[3]*sig[2] - sig[5]*sig[4];
  xmin3 = sig[3]*sig[5] - sig[4]*sig[1];
  *I3 = sig[0]*xmin1 - sig[3]*xmin2 + sig[4]*xmin3;

  /* rotated stress tensor (beta shear-band rotation) */
  vertical = parms[17];
  beta_deg = parms[22];
  for (i = 3; i < 6; i++) fullbeta[i] = 0.0;
  beta = beta_deg * 3.14159265358979323846 / 180.0;
  trace_sig = sig[0] + sig[1] + sig[2];

  if (vertical == 1) {
    fullbeta[0] = -2.0 * tan(beta) * trace_sig / 9.0;
    fullbeta[1] = tan(beta) * trace_sig / 9.0;
    fullbeta[2] = tan(beta) * trace_sig / 9.0;
  } else if (vertical == 2) {
    fullbeta[1] = -2.0 * tan(beta) * trace_sig / 9.0;
    fullbeta[0] = tan(beta) * trace_sig / 9.0;
    fullbeta[2] = tan(beta) * trace_sig / 9.0;
  } else {
    fullbeta[2] = -2.0 * tan(beta) * trace_sig / 9.0;
    fullbeta[0] = tan(beta) * trace_sig / 9.0;
    fullbeta[1] = tan(beta) * trace_sig / 9.0;
  }
  for (i = 0; i < 6; i++) sig_rot[i] = fullbeta[i] + sig[i];

  *I1rot = sig_rot[0] + sig_rot[1] + sig_rot[2];
  norm2sig_rot = dot_vect(1, sig_rot, sig_rot, 6);
  *I2rot = half * (norm2sig_rot - (*I1rot)*(*I1rot));

  xmin1_rot = sig_rot[1]*sig_rot[2] - sig_rot[5]*sig_rot[5];
  xmin2_rot = sig_rot[3]*sig_rot[2] - sig_rot[5]*sig_rot[4];
  xmin3_rot = sig_rot[3]*sig_rot[5] - sig_rot[4]*sig_rot[1];
  *I3rot = sig_rot[0]*xmin1_rot - sig_rot[3]*xmin2_rot + sig_rot[4]*xmin3_rot;

  if (*I1rot != 0.0) {
    for (i = 0; i < 6; i++) eta_rot[i] = sig_rot[i] / (*I1);
  } else {
    for (i = 0; i < 6; i++) eta_rot[i] = sig_rot[i] / tiny;
  }
  for (i = 0; i < 3; i++) eta_d_rot[i] = eta_rot[i] - onethird;
  for (i = 3; i < 6; i++) eta_d_rot[i] = eta[i];

  norm2eta_rot = dot_vect(1, eta_d_rot, eta_d_rot, 6);

  eta_d2_rot[0] = eta_d_rot[0]*eta_d_rot[0] + eta_d_rot[3]*eta_d_rot[3] +
    eta_d_rot[4]*eta_d_rot[4];
  eta_d2_rot[1] = eta_d_rot[3]*eta_d_rot[3] + eta_d_rot[1]*eta_d_rot[1] +
    eta_d_rot[5]*eta_d_rot[5];
  eta_d2_rot[2] = eta_d_rot[5]*eta_d_rot[5] + eta_d_rot[4]*eta_d_rot[4] +
    eta_d_rot[2]*eta_d_rot[2];
  eta_d2_rot[3] = eta_d_rot[0]*eta_d_rot[3] + eta_d_rot[3]*eta_d_rot[1] +
    eta_d_rot[5]*eta_d_rot[4];
  eta_d2_rot[4] = eta_d_rot[4]*eta_d_rot[0] + eta_d_rot[5]*eta_d_rot[3] +
    eta_d_rot[2]*eta_d_rot[4];
  eta_d2_rot[5] = eta_d_rot[3]*eta_d_rot[4] + eta_d_rot[1]*eta_d_rot[5] +
    eta_d_rot[5]*eta_d_rot[2];

  if (norm2eta_rot < tiny) {
    *cos3t_rot = -one;
  } else {
    tretadev3_rot = dot_vect(1, eta_d_rot, eta_d2_rot, 6);
    numer_rot = -sqrt6 * tretadev3_rot;
    denom_rot = pow(sqrt(norm2eta_rot), 3.0);
    *cos3t_rot = numer_rot / denom_rot;
    if (fabs(*cos3t) > one) *cos3t = *cos3t / fabs(*cos3t);
  }

  (void)three;
}

/* ------------------------------------------------------------------ */
/* tangent operators MM, HH, LL, NN for the Masin model                */
/* istrain=1: full model with intergranular strains                    */
/* istrain=0: basic model without intergranular strains                */
/* ------------------------------------------------------------------ */
static void get_tan(double *deps, double *sig, double *q, int nasv,
  double *parms, int nparms, double *MM, double *HH, double *LL,
  double *LL_unl, double *hypo_Dsom_ld, double *NN, int istrain, int *error)
{
  double eta[6], eta_dev[6], del[6], evoid, sig_star[6], sensit;
  double H_s[6], eta_del[6], eta_delta[6], eta_eps[6];
  double norm_del, norm_del2, norm_deps, norm_deps2;
  double pp, qq, cos3t, cos3t_rot, I1, I2, I3, I1rot, I2rot, I3rot;
  double a, a2, alpha, fd, fs, fdsbs, fddivfdA;
  double H_del[36], H_e[6], IU[36], AAhce[36], NNvec[6];
  double krondelta[6], hypo_Dsom_rot[6];
  double load, rho, N_par, Stf, kparam, Aparam, sfparam;
  double kap_par, lambda_struct;
  double zero=0.0, one=1.0, two=2.0, three=3.0;
  double tiny=1.0e-17, half=0.5;
  double onethird, sqrt3, twosqrt2, ln2m1;
  double temp1, temp2, alpha_power;
  double phi, lam_star, kap_star, N_star, nuvh;
  double m_T, beta_r, chi, p_t, sinphi, sinphi2;
  double nparam, m_Trat, G0, Gvh, alphanu;
  double alphaG, alphaE, Am, nuhh, pmeangt1;
  double npow, cos2phic, ocparam, peast;
  double pmean, kpow, sinphickpow, sinphimkpow, Amult;
  double nvect[3], pmat[9], kron_delta[9];
  double kck[81], kdk[81], pdk[81], kdp[81], pck[81], pdp[81], LLfour[81];
  double an1, an2, an3, an4, an5, norm_m, norm_m2;
  double eta_rot[6], eta_dev_rot[6], sig_rot[6];
  double Fmfactor_rot, Dref;
  double trace_sig, gama, gama_deg, tr_Dsom_rot, mintwooneone[6];
  double Dsominproduct, tangama, gamafact, fullgama[6];
  double norm_Dsom, norm_Dsom2, norm_tr_hypo_Dsom_rot;
  double wyfact, acorr_fact, acorrwy, wy, ksi;
  int i, j, k, l, softmodel, vertical;

  onethird = one/three;
  sqrt3 = sqrt(three);
  twosqrt2 = two*sqrt(two);
  ln2m1 = one/log(two);

  for (i = 0; i < 6; i++)
    for (j = 0; j < 6; j++) {
      MM[i*6+j] = zero; LL[i*6+j] = zero; LL_unl[i*6+j] = zero;
      IU[i*6+j] = zero; H_del[i*6+j] = zero;
    }
  for (i = 0; i < 6; i++) { eta_del[i] = zero; eta_delta[i] = zero; eta_eps[i] = zero; }
  for (i = 0; i < nasv; i++)
    for (j = 0; j < 6; j++) HH[i*6+j] = zero;

  IU[0*6+0]=one; IU[1*6+1]=one; IU[2*6+2]=one;
  IU[3*6+3]=one; IU[4*6+4]=one; IU[5*6+5]=one;

  /* recover material parameters */
  phi = parms[0];
  p_t = parms[1];
  lam_star = parms[2];
  kap_par = parms[3];
  N_par = parms[4];
  nuhh = parms[5];
  alphaG = parms[6];
  if (alphaG < 0.01) alphaG = 1.0;
  kparam = parms[7];
  Aparam = parms[8];
  sfparam = parms[9];
  beta_r = parms[11];
  chi = parms[12];
  G0 = parms[13];
  nparam = parms[14];
  m_Trat = parms[15];
  vertical = (int)parms[17];
  alphaE = parms[18];
  alphanu = parms[19];
  if (alphaE < 0.01) alphaE = pow(alphaG, 1.25);
  if (alphanu < 0.01) alphanu = alphaG;
  ocparam = parms[21];
  beta_r = parms[11];
  ksi = parms[23];
  gama_deg = parms[24];
  Dref = parms[25];

  sinphi = sin(phi);
  sinphi2 = sinphi*sinphi;
  nuvh = nuhh/alphanu;

  a = sqrt3*(three - sin(phi))/(twosqrt2*sin(phi));
  a2 = a*a;
  temp1 = (lam_star - kap_par)/(lam_star + kap_par);
  temp2 = (three + a2)/(sqrt3*a);
  alpha_power = ln2m1*log(temp1*temp2);
  if (nparms >= 21 && parms[20] >= 0.00001) alpha_power = parms[20];

  /* recover internal state variables */
  for (i = 0; i < 6; i++) del[i] = q[i];
  evoid = q[6];
  sensit = q[7];

  /* soft clay model */
  softmodel = 0;
  Stf = 1;
  N_star = N_par;
  kap_star = kap_par;
  lambda_struct = lam_star;
  if (sensit >= 1) {
    softmodel = 1;
    N_star = N_par + lam_star*log(sensit);
    Stf = (sensit - kparam*(sensit - sfparam))/sensit;
    lambda_struct = lam_star/Stf;
  }

  /* axis translation due to cohesion */
  sig_star[0]=sig[0]-p_t; sig_star[1]=sig[1]-p_t; sig_star[2]=sig[2]-p_t;
  sig_star[3]=sig[3]; sig_star[4]=sig[4]; sig_star[5]=sig[5];

  /* strain increment and intergranular strain directions */
  norm_deps2 = dot_vect(2, deps, deps, 6);
  norm_del2 = dot_vect(2, del, del, 6);
  norm_deps = sqrt(norm_deps2);
  norm_del = sqrt(norm_del2);

  if (norm_del >= tiny)
    for (i = 0; i < 6; i++) eta_del[i] = del[i]/norm_del;

  eta_delta[0]=eta_del[0]; eta_delta[1]=eta_del[1]; eta_delta[2]=eta_del[2];
  eta_delta[3]=half*eta_del[3]; eta_delta[4]=half*eta_del[4]; eta_delta[5]=half*eta_del[5];

  if (norm_deps >= tiny)
    for (i = 0; i < 6; i++) eta_eps[i] = deps[i]/norm_deps;

  /* auxiliary stress tensors (with rotation) */
  inv_sig(sig_star, &pp, &qq, &cos3t, &cos3t_rot, &I1, &I2, &I3,
    &I1rot, &I2rot, &I3rot, sig_rot, parms);

  for (i = 0; i < 6; i++) eta[i] = sig_star[i]/I1;

  eta_dev[0]=eta[0]-onethird; eta_dev[1]=eta[1]-onethird; eta_dev[2]=eta[2]-onethird;
  eta_dev[3]=eta[3]; eta_dev[4]=eta[4]; eta_dev[5]=eta[5];

  krondelta[0]=one; krondelta[1]=one; krondelta[2]=one;
  krondelta[3]=zero; krondelta[4]=zero; krondelta[5]=zero;

  for (i = 0; i < 6; i++) eta_rot[i] = sig_rot[i]/I1rot;
  for (i = 0; i < 3; i++) eta_dev_rot[i] = eta_rot[i] - onethird;
  for (i = 3; i < 6; i++) eta_dev_rot[i] = eta_rot[i];

  /* explicit clay hypoplasticity specific (rotated invariants) */
  peast = exp((N_star - log(one+evoid))/lam_star);

  if ((I3rot + I1rot*I2rot) != 0) Fmfactor_rot = (9*I3rot + I1rot*I2rot)/(I3rot + I1rot*I2rot);
  else Fmfactor_rot = 1;
  if (Fmfactor_rot > 1) Fmfactor_rot = 1;
  if (Fmfactor_rot < 0) Fmfactor_rot = 0;

  cos2phic = 1 - sin(phi)*sin(phi);
  npow = -log(cos2phic)/log(ocparam) + 0.30*(Fmfactor_rot - sin(phi)*sin(phi));
  fdsbs = ocparam*pow(1 - Fmfactor_rot, 1/npow);
  pmean = -I1/3.0;
  fd = pow((ocparam*pmean)/peast, alpha_power);
  fdsbs = pow(fdsbs, alpha_power);
  fddivfdA = fd/fdsbs;

  if (Fmfactor_rot < 1.e-10) { Fmfactor_rot = 0; cos3t = -1; }

  Am = nuvh*nuvh*(4*alphaE*alphanu - 2*alphaE*alphaE*alphanu*alphanu +
        2*alphaE*alphaE - alphanu*alphanu) +
    nuvh*(4*alphaE + 2*alphaE*alphanu) + 1 + 2*alphaE;

  fs = 9*pmean/2*(1/kap_par + 1/lambda_struct)/Am;

  for (i = 0; i < 3; i++) nvect[i] = 0;
  if (vertical >= 1 && vertical <= 3) nvect[vertical-1] = 1;
  else { *error = 10; return; }

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++) {
      pmat[i*3+j] = nvect[i]*nvect[j];
      kron_delta[i*3+j] = 0;
    }
  kron_delta[0*3+0]=1; kron_delta[1*3+1]=1; kron_delta[2*3+2]=1;

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      for (k = 0; k < 3; k++)
        for (l = 0; l < 3; l++) {
          int idx = (i*3+j)*9 + (k*3+l);
          kck[idx] = (kron_delta[i*3+k]*kron_delta[j*3+l] +
            kron_delta[i*3+l]*kron_delta[j*3+k] +
            kron_delta[j*3+l]*kron_delta[i*3+k] +
            kron_delta[j*3+k]*kron_delta[i*3+l])/2;
          kdk[idx] = kron_delta[i*3+j]*kron_delta[k*3+l];
          pdk[idx] = pmat[i*3+j]*kron_delta[k*3+l];
          kdp[idx] = kron_delta[i*3+j]*pmat[k*3+l];
          pck[idx] = (pmat[i*3+k]*kron_delta[j*3+l] +
            pmat[i*3+l]*kron_delta[j*3+k] +
            pmat[j*3+l]*kron_delta[i*3+k] +
            pmat[j*3+k]*kron_delta[i*3+l])/2;
          pdp[idx] = pmat[i*3+j]*pmat[k*3+l];
        }

  an1 = alphaE*(1 - alphanu*nuvh - 2*alphaE*nuvh*nuvh);
  an2 = alphaE*nuvh*(alphanu + alphaE*nuvh);
  an3 = alphaE*nuvh*(1 + alphanu*nuvh - alphanu - alphaE*nuvh);
  an4 = (1 - alphanu*nuvh - 2*alphaE*nuvh*nuvh)*(alphaE*(1 - alphaG))/alphaG;
  an5 = alphaE*(1 - alphaE*nuvh*nuvh) + 1 - alphanu*alphanu*nuvh*nuvh -
    2*alphaE*nuvh*(1 + alphanu*nuvh) -
    2*alphaE*(1 - alphanu*nuvh - 2*alphaE*nuvh*nuvh)/alphaG;

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      for (k = 0; k < 3; k++)
        for (l = 0; l < 3; l++) {
          int idx = (i*3+j)*9 + (k*3+l);
          LLfour[idx] = an1*kck[idx]/2 + an2*kdk[idx] +
            an3*(pdk[idx] + kdp[idx]) + an4*pck[idx] + an5*pdp[idx];
        }

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++) LL[i*6+j] = LLfour[(i*3+i)*9 + (j*3+j)];
  LL[3*6+3] = LLfour[(0*3+1)*9 + (0*3+1)];
  LL[4*6+4] = LLfour[(0*3+2)*9 + (0*3+2)];
  LL[5*6+5] = LLfour[(1*3+2)*9 + (1*3+2)];

  for (i = 0; i < 6; i++) hypo_Dsom_rot[i] = 0.0;

  kpow = 1.70 + 3.90*sin(phi)*sin(phi);
  sinphickpow = pow(sin(phi), kpow);
  sinphimkpow = pow(sqrt(Fmfactor_rot), kpow);

  Amult = 2.0/3.0 - sqrt(sqrt(Fmfactor_rot))*(cos3t_rot + 1.0)/4.0;

  for (i = 0; i < 6; i++)
    hypo_Dsom_rot[i] = -eta_dev_rot[i] + krondelta[i]*
      (sinphimkpow - sinphickpow)/(1 - sinphickpow)*Amult;
  for (i = 3; i < 6; i++) hypo_Dsom_rot[i] = 2.0*hypo_Dsom_rot[i];

  norm_m2 = dot_vect(2, hypo_Dsom_rot, hypo_Dsom_rot, 6);
  norm_m = sqrt(norm_m2);
  for (i = 0; i < 6; i++) hypo_Dsom_rot[i] = hypo_Dsom_rot[i]/norm_m;

  for (i = 0; i < 6; i++)
    for (j = 0; j < 6; j++)
      AAhce[i*6+j] = fs*LL[i*6+j] + sig_star[i]*krondelta[j]/lambda_struct;

  matmul(AAhce, hypo_Dsom_rot, NNvec, 6, 6, 1);
  for (i = 0; i < 6; i++) NN[i] = -NNvec[i]*fddivfdA/(fd*fs);

  /* end basic model; void ratio evolution (tension positive) */
  for (i = 0; i < 6; i++) H_e[i] = (i < 3) ? (one+evoid) : zero;

  for (i = 0; i < 6; i++) {
    H_s[i] = zero;
    if (softmodel == 1) {
      if (i < 3) H_s[i] = -kparam*(sensit - sfparam)/lam_star;
    }
  }

  for (i = 0; i < nasv; i++) {
    if (i < 6) {
      for (j = 0; j < 6; j++) HH[i*6+j] = 0;
    } else if (i == 6) {
      for (j = 0; j < 6; j++) HH[i*6+j] = H_e[j];
    } else if (i == 7) {
      for (j = 0; j < 6; j++) HH[i*6+j] = H_s[j];
    }
  }

  for (i = 0; i < 6; i++)
    for (j = 0; j < 6; j++) LL[i*6+j] = LL[i*6+j]*fs;
  for (i = 0; i < 6; i++) NN[i] = NN[i]*fs*fd;

  /* LD approach */
  trace_sig = sig[0]+sig[1]+sig[2];
  gama = gama_deg*3.14159265358979323846/180.0;

  tr_Dsom_rot = hypo_Dsom_rot[0]+hypo_Dsom_rot[1]+hypo_Dsom_rot[2];

  if (vertical == 1) { mintwooneone[0]=-2.0; mintwooneone[1]=1.0; mintwooneone[2]=1.0; }
  else if (vertical == 2) { mintwooneone[0]=1.0; mintwooneone[1]=-2.0; mintwooneone[2]=1.0; }
  else { mintwooneone[0]=1.0; mintwooneone[1]=1.0; mintwooneone[2]=-2.0; }
  for (i = 3; i < 6; i++) mintwooneone[i] = 0.0;

  Dsominproduct = dot_vect(3, hypo_Dsom_rot, mintwooneone, 6);
  tangama = -Dsominproduct/(2*tr_Dsom_rot);

  if (gama < -1.570796) gama = atan(tangama);

  gamafact = 3*tan(gama)*tr_Dsom_rot/9.0;

  for (i = 0; i < 6; i++) fullgama[i] = 0.0;
  if (vertical == 1) { fullgama[0]=-2.0*gamafact; fullgama[1]=gamafact; fullgama[2]=gamafact; }
  else if (vertical == 2) { fullgama[0]=gamafact; fullgama[1]=-2.0*gamafact; fullgama[2]=gamafact; }
  else { fullgama[0]=gamafact; fullgama[1]=gamafact; fullgama[2]=-2.0*gamafact; }

  for (i = 0; i < 6; i++) hypo_Dsom_ld[i] = hypo_Dsom_rot[i] + fullgama[i];

  norm_Dsom2 = dot_vect(2, hypo_Dsom_rot, hypo_Dsom_rot, 6);
  norm_Dsom = sqrt(norm_Dsom2);

  if (fabs(norm_Dsom) > 0.00000001)
    norm_tr_hypo_Dsom_rot = (hypo_Dsom_rot[0]+hypo_Dsom_rot[1]+hypo_Dsom_rot[2])/norm_Dsom;

  wyfact = sqrt(one/three)*norm_tr_hypo_Dsom_rot + 1;
  if (wyfact > 1) wyfact = 1;

  wy = pow(wyfact, ksi);
  if (ksi > 9.5) wy = 0;

  acorr_fact = dot_vect(2, hypo_Dsom_rot, hypo_Dsom_ld, 6);
  acorrwy = 1 - wy;
  if (fabs(acorr_fact) > 0.00000001) acorrwy = (1 - wy)/acorr_fact;

  if (Dref < 0.000000000001) {
    for (i = 0; i < 6; i++)
      for (j = 0; j < 3; j++)
        LL_unl[i*6+j] = LL[i*6+j] - acorrwy*NN[i]*hypo_Dsom_ld[j];
    for (i = 0; i < 6; i++)
      for (j = 3; j < 6; j++)
        LL_unl[i*6+j] = LL[i*6+j] - acorrwy*NN[i]*hypo_Dsom_ld[j]/2.0;
    for (i = 0; i < 6; i++)
      for (j = 0; j < 3; j++)
        LL[i*6+j] = LL[i*6+j] + acorrwy*NN[i]*hypo_Dsom_ld[j];
    for (i = 0; i < 6; i++)
      for (j = 3; j < 6; j++)
        LL[i*6+j] = LL[i*6+j] + acorrwy*NN[i]*hypo_Dsom_ld[j]/2.0;
  } else {
    for (i = 0; i < 6; i++)
      for (j = 0; j < 3; j++)
        LL_unl[i*6+j] = acorrwy*NN[i]*hypo_Dsom_ld[j];
    for (i = 0; i < 6; i++)
      for (j = 3; j < 6; j++)
        LL_unl[i*6+j] = acorrwy*NN[i]*hypo_Dsom_ld[j]/2.0;
  }

  for (i = 0; i < 6; i++) NN[i] = NN[i]*wy;

  (void)alpha; (void)load; (void)rho; (void)m_T; (void)Gvh; (void)pmeangt1;
  (void)nparam; (void)m_Trat; (void)beta_r; (void)nparms; (void)chi;
  (void)sinphi2;
}

/* F_sig = MM*deps (istrain=1) or LL*deps+NN*|deps| (istrain=0)        */
/* F_q   = HH*deps with soft-model correction                          */
/* ------------------------------------------------------------------ */
static void get_F_sig_q(double *sig, double *q, int nasv, double *parms,
  int nparms, double *deps, double *F_sig, double *F_q, double dtime,
  int *error)
{
  double MM[36], HH[NASV*6], LL[36], LL_unl[36], hypo_Dsom_ld[6], NN[6];
  double Ltoadd[36], norm_D, norm_D2, norm2;
  double trD, depsh, depd, Aparam, edev[6];
  double Dref, depsDdir;
  int istrain, ii, jj;

  istrain = 0;
  Dref = parms[25];

  get_tan(deps, sig, q, nasv, parms, nparms, MM, HH, LL, LL_unl,
    hypo_Dsom_ld, NN, istrain, error);

  if (Dref < 0.000000000001) {
    depsDdir = dot_vect(2, deps, hypo_Dsom_ld, 6);
    if (depsDdir > 0) {
      /* LL stays */
    } else {
      for (ii = 0; ii < 6; ii++)
        for (jj = 0; jj < 6; jj++) LL[ii*6+jj] = LL_unl[ii*6+jj];
    }
    matmul(LL, deps, F_sig, 6, 6, 1);
    norm_D2 = dot_vect(2, deps, deps, 6);
    norm_D = sqrt(norm_D2);
    for (ii = 0; ii < 6; ii++) F_sig[ii] = F_sig[ii] + NN[ii]*norm_D;
  } else {
    for (ii = 0; ii < 6; ii++) deps[ii] = deps[ii]/dtime;
    norm_D2 = dot_vect(2, deps, deps, 6);
    norm_D = sqrt(norm_D2);
    if (norm_D > 0.00000001) {
      for (ii = 0; ii < 6; ii++)
        for (jj = 0; jj < 6; jj++) Ltoadd[ii*6+jj] = LL_unl[ii*6+jj]*Dref/norm_D;
    }
    if (dot_vect(2, deps, hypo_Dsom_ld, 6) < 0.0) {
      for (ii = 0; ii < 6; ii++)
        for (jj = 0; jj < 6; jj++) Ltoadd[ii*6+jj] = -Ltoadd[ii*6+jj];
    }
    for (ii = 0; ii < 6; ii++)
      for (jj = 0; jj < 6; jj++) LL[ii*6+jj] = LL[ii*6+jj] + Ltoadd[ii*6+jj];
    matmul(LL, deps, F_sig, 6, 6, 1);
    for (ii = 0; ii < 6; ii++) F_sig[ii] = F_sig[ii] + NN[ii]*Dref;
    for (ii = 0; ii < 6; ii++) deps[ii] = deps[ii]*dtime;
    for (ii = 0; ii < 6; ii++) F_sig[ii] = F_sig[ii]*dtime;
  }

  matmul(HH, deps, F_q, nasv, 6, 1);

  trD = deps[0]+deps[1]+deps[2];
  edev[0]=deps[0]-trD/3.; edev[1]=deps[1]-trD/3.; edev[2]=deps[2]-trD/3.;
  edev[3]=deps[3]/2.; edev[4]=deps[4]/2.; edev[5]=deps[5]/2.;
  norm2 = edev[0]*edev[0]+edev[1]*edev[1]+edev[2]*edev[2] +
    2.*(edev[3]*edev[3]+edev[4]*edev[4]+edev[5]*edev[5]);
  depsh = sqrt(norm2*2./3.);
  Aparam = parms[8];
  depd = sqrt(trD*trD + Aparam/(1-Aparam)*depsh*depsh);

  F_q[nasv-1] = HH[7*6+0]*depd;
}


static void rhs(double *y, int ny, int nasv, double *parms, int nparms,
  double *deps, double *kRK, int *nfev, double dtime, int *error)
{
  double sig[6], q[NASV], F_sig[6], F_q[NASV];
  int i;

  (*nfev)++;
  for (i = 0; i < ny; i++) kRK[i] = 0.0;

  for (i = 0; i < 6; i++) sig[i] = y[i];
  for (i = 0; i < nasv; i++) q[i] = y[6+i];

  get_F_sig_q(sig, q, nasv, parms, nparms, deps, F_sig, F_q, dtime, error);
  if (*error == 10) return;

  for (i = 0; i < 6; i++) kRK[i] = F_sig[i];
  for (i = 0; i < nasv; i++) kRK[6+i] = F_q[i];
}

/* ------------------------------------------------------------------ */
/* relative error norm for the RKF23 step                              */
/* ------------------------------------------------------------------ */
static void norm_res(double *y_til, double *y_hat, int ny, int nasv,
  double *norm_R)
{
  int i;
  double sig_hat[6], sig_til[6], del_sig[6];
  double q_hat[8], q_til[8], del_q[8];
  double void_hat, void_til, del_void;
  double sensit_hat, sensit_til, del_sensit;
  double err[NYDIM];
  double norm_R2, norm_sig2, norm_q2, norm_sig, norm_q;
  double zero = 0.0;
  int testnan = 0;

  for (i = 0; i < ny; i++) err[i] = zero;

  for (i = 0; i < 6; i++) {
    sig_hat[i] = y_hat[i]; sig_til[i] = y_til[i];
    del_sig[i] = fabs(sig_hat[i] - sig_til[i]);
  }
  for (i = 0; i < 6; i++) {
    q_hat[i] = y_hat[6+i]; q_til[i] = y_til[6+i];
    del_q[i] = fabs(q_hat[i] - q_til[i]);
  }
  void_hat = y_hat[6+6]; void_til = y_til[6+6];
  del_void = fabs(void_hat - void_til);

  sensit_hat = y_hat[6+7]; sensit_til = y_til[6+7];
  del_sensit = fabs(sensit_hat - sensit_til);

  norm_sig2 = dot_vect(1, sig_hat, sig_hat, 6);
  norm_q2 = dot_vect(2, q_hat, q_hat, 6);
  norm_sig = sqrt(norm_sig2);
  norm_q = sqrt(norm_q2);

  if (norm_sig > zero)
    for (i = 0; i < 6; i++) err[i] = del_sig[i]/norm_sig;
  if (norm_q > zero)
    for (i = 0; i < 6; i++) err[6+i] = del_q[i]/norm_q;

  /* Fortran: err(6+nasv-1)=del_void/void_hat  -> 1-based idx 13 -> 0-based 12 */
  err[12] = del_void/void_hat;
  /* Fortran: err(6+8)=del_sensit/sensit_hat   -> 1-based idx 14 -> 0-based 13 */
  err[13] = 0;
  if (sensit_hat > 0) err[13] = del_sensit/sensit_hat;

  norm_R2 = dot_vect(3, err, err, ny);
  *norm_R = sqrt(norm_R2);

  testnan = 0;
  umatisnan(norm_sig, &testnan);
  umatisnan(norm_q, &testnan);
  umatisnan(void_hat, &testnan);
  umatisnan(sensit_hat, &testnan);
  if (testnan == 1) *norm_R = 1.0e20;
}

/* ------------------------------------------------------------------ */
/* check whether the RKF23 solution vector is admissible               */
/* ------------------------------------------------------------------ */
static void check_RKF(int *error_RKF, double *y, int ny, int nasv,
  double *parms, int nparms)
{
  double sig[6], pmean, sig_star[6];
  double I1, I2, I3, pp, qq, cos3t, cos3t_rot;
  double I1rot, I2rot, I3rot, sig_rot[6];
  double p_t, minstress, tolerance;
  double OCR, omega, fSBS, sensit, cos2phic;
  double ocparam, Fmfactor_rot;
  double ashape;
  int i, testnan;

  p_t = parms[1];
  minstress = p_t/100.0;
  for (i = 0; i < 6; i++) sig[i] = y[i];

  sig_star[0]=sig[0]-p_t; sig_star[1]=sig[1]-p_t; sig_star[2]=sig[2]-p_t;
  sig_star[3]=sig[3]; sig_star[4]=sig[4]; sig_star[5]=sig[5];

  inv_sig(sig_star, &pp, &qq, &cos3t, &cos3t_rot, &I1, &I2, &I3,
    &I1rot, &I2rot, &I3rot, sig_rot, parms);
  pmean = -I1/3;
  if (pmean < minstress) *error_RKF = 1;

  if ((I3rot + I1rot*I2rot) != 0) Fmfactor_rot = (9*I3rot + I1rot*I2rot)/(I3rot + I1rot*I2rot);
  else Fmfactor_rot = 1;

  sensit = 1;
  if (y[6+7] >= 1) sensit = y[6+7];
  OCR = -sensit*exp((parms[4] - log(1 + y[6+6]))/parms[2])/pp;
  cos2phic = 1 - sin(parms[0])*sin(parms[0]);

  ashape = 0.3;
  ocparam = parms[21];
  omega = -log(cos2phic)/log(ocparam) + ashape*(Fmfactor_rot - sin(parms[0])*sin(parms[0]));
  fSBS = Fmfactor_rot + pow(1/OCR, omega) - 1;

  if (Fmfactor_rot >= 1) *error_RKF = 1;

  tolerance = 0.1;
  if (OCR <= ocparam) {
    if (Fmfactor_rot >= (1 - cos2phic + tolerance)) *error_RKF = 1;
  }
  if (OCR > ocparam) {
    if (fSBS >= tolerance) *error_RKF = 1;
    if (fSBS >= tolerance) *error_RKF = 1;
  }

  testnan = 0;
  for (i = 0; i < ny; i++) umatisnan(y[i], &testnan);
  umatisnan(0.0, &testnan);   /* sin2phim unused in visco */
  umatisnan(OCR, &testnan);
  umatisnan(fSBS, &testnan);
  if (testnan == 1) *error_RKF = 1;

  (void)nasv; (void)nparms;
}


static void rkf23_update(double *y, int n, int nasv, double *dtsub,
  double err_tol, int maxnint, double DTmin, double *deps_np1,
  double *parms, int nparms, int *nfev, double dtime, int *error)
{
  int i, ksubst, kreject, error_RKF;
  double y_k[NYDIM], y_2[NYDIM], y_3[NYDIM], y_til[NYDIM], y_hat[NYDIM];
  double kRK_1[NYDIM], kRK_2[NYDIM], kRK_3[NYDIM];
  double T_k, DT_k, norm_R, S_hull;
  double one=1.0, two=2.0, three=3.0, four=4.0, six=6.0;
  double half=0.5, ptnine=0.9;
  double onesixth, onethird, twothirds, temp;

  for (i = 0; i < n; i++) y_k[i] = 0.0;
  onesixth = one/six;
  onethird = one/three;
  twothirds = two/three;

  error_RKF = 0;
  T_k = 0.0;
  DT_k = (*dtsub)/dtime;
  ksubst = 0;
  kreject = 0;
  *nfev = 0;
  for (i = 0; i < n; i++) y_k[i] = y[i];

  while (T_k < one) {
    ksubst = ksubst + 1;
    if (ksubst > maxnint) { *error = 3; return; }

    error_RKF = 0;
    check_RKF(&error_RKF, y_k, n, nasv, parms, nparms);
    if (error_RKF == 1) { *error = 3; return; }
    else rhs(y_k, n, nasv, parms, nparms, deps_np1, kRK_1, nfev, dtime, error);
    if (*error == 10) return;

    temp = half*DT_k;
    for (i = 0; i < n; i++) y_2[i] = y_k[i] + temp*kRK_1[i];

    error_RKF = 0;
    check_RKF(&error_RKF, y_2, n, nasv, parms, nparms);
    if (error_RKF == 1) { *error = 3; return; }
    else rhs(y_2, n, nasv, parms, nparms, deps_np1, kRK_2, nfev, dtime, error);
    if (*error == 10) return;

    for (i = 0; i < n; i++)
      y_3[i] = y_k[i] - DT_k*kRK_1[i] + two*DT_k*kRK_2[i];

    error_RKF = 0;
    check_RKF(&error_RKF, y_3, n, nasv, parms, nparms);
    if (error_RKF == 1) { *error = 3; return; }
    else rhs(y_3, n, nasv, parms, nparms, deps_np1, kRK_3, nfev, dtime, error);
    if (*error == 10) return;

    for (i = 0; i < n; i++) y_til[i] = y_k[i] + DT_k*kRK_2[i];
    for (i = 0; i < n; i++)
      y_hat[i] = y_k[i] + DT_k*
        (onesixth*kRK_1[i] + twothirds*kRK_2[i] + onesixth*kRK_3[i]);

    norm_res(y_til, y_hat, n, nasv, &norm_R);

    error_RKF = 0;
    check_RKF(&error_RKF, y_hat, n, nasv, parms, nparms);
    if (error_RKF != 0) { *error = 3; return; }

    if (norm_R != 0) S_hull = ptnine*DT_k*pow(err_tol/norm_R, onethird);
    else S_hull = 1;

    if (norm_R < err_tol) {
      for (i = 0; i < n; i++) y_k[i] = y_hat[i];
      T_k = T_k + DT_k;
      DT_k = fmin(four*DT_k, S_hull);
      *dtsub = DT_k*dtime;
      DT_k = fmin((one - T_k), DT_k);
    } else {
      DT_k = fmax(DT_k/four, S_hull);
      if (DT_k < DTmin) { *error = 3; return; }
    }
  }

  for (i = 0; i < n; i++) y[i] = y_k[i];

  (void)kreject; (void)three; (void)half;
}

/* ------------------------------------------------------------------ */
/* compute numerically consistent tangent (elastic-like, scaled L)     */
/* ------------------------------------------------------------------ */
static void perturbate(double *y_n, int nasv, double *parms, int nparms,
  double *deps_np1, double *DD, int *error)
{
  double sig[6], q[NASV], deps[6];
  double MM[36], HH[NASV*6], LL[36], LL_unl[36], hypo_Dsom_ld[6], NN[6];
  double p_t;
  int istrain, jj, kk;

  p_t = parms[1];
  for (kk = 0; kk < 6; kk++) sig[kk] = y_n[kk];

  if (-(sig[0]+sig[1]+sig[2])/3 < p_t) {
    sig[0]=-p_t; sig[1]=-p_t; sig[2]=-p_t;
    sig[3]=0; sig[4]=0; sig[5]=0;
  }

  for (kk = 0; kk < nasv; kk++) q[kk] = y_n[6+kk];
  for (kk = 0; kk < 6; kk++) deps[kk] = 0.0;

  istrain = 0;

  for (kk = 0; kk < 6; kk++)
    for (jj = 0; jj < 6; jj++) DD[kk*6+jj] = 0.0;

  if (*error != 10) {
    get_tan(deps, sig, q, nasv, parms, nparms, MM, HH, LL, LL_unl,
      hypo_Dsom_ld, NN, istrain, error);
  }
  for (kk = 0; kk < 6; kk++)
    for (jj = 0; jj < 6; jj++) DD[kk*6+jj] = LL[kk*6+jj];

  (void)nparms; (void)deps_np1;
}


static void move_sig(double *stress, double pore, double *sig)
{
  int i;
  for (i = 0; i < 6; i++) sig[i] = 0.0;
  for (i = 0; i < 3; i++) sig[i] = stress[i] + pore;
  sig[3] = stress[3]; sig[4] = stress[4]; sig[5] = stress[5];
}

/* y = [sig, qq] */
static void iniy(double *y, double *sig, double *qq)
{
  int i;
  for (i = 0; i < NYDIM; i++) y[i] = 0.0;
  for (i = 0; i < 6; i++) y[i] = sig[i];
  for (i = 0; i < NASV; i++) y[6+i] = qq[i];
}

/* output: stress, asv, ddsdde from solution y and tangent DD */
static void solout(double *stress, double *asv, int nasv, double *y,
  double pore, double depsv_np1, double *parms, double *DD, double *ddsdde)
{
  int i, j;
  double bulk_w = parms[16];

  pore = pore - bulk_w*depsv_np1;

  for (i = 0; i < 3; i++) stress[i] = y[i] - pore;
  for (i = 3; i < 6; i++) stress[i] = y[i];

  for (i = 0; i < nasv; i++) asv[i] = y[6+i];

  for (j = 0; j < 6; j++)
    for (i = 0; i < 6; i++) ddsdde[i*6+j] = DD[i*6+j];
  for (j = 0; j < 3; j++)
    for (i = 0; i < 3; i++) ddsdde[i*6+j] = ddsdde[i*6+j] + bulk_w;
}

/* additional state variables for postprocessing: statev(11),(12),(15) */
static void calc_statev(double *stress, double *statev, double *parms,
  int nparms, int nasv, double *deps)
{
  double sig_star[6], I1, I2, I3, cos3t, pp, qq;
  double sin2phi, sinphi, p_t;
  double norm_del, norm_del2, del[6], sensit;
  double cos3t_rot, I1rot, I2rot, I3rot, sig_rot[6];
  int i;

  p_t = parms[1];
  for (i = 0; i < 3; i++) sig_star[i] = stress[i]-p_t;
  for (i = 3; i < 6; i++) sig_star[i] = stress[i];

  inv_sig(sig_star, &pp, &qq, &cos3t, &cos3t_rot, &I1, &I2, &I3,
    &I1rot, &I2rot, &I3rot, sig_rot, parms);
  if (I3 != 0) sin2phi = (9. + I1*I2/I3)/(1. + I1*I2/I3);
  else sin2phi = 0;
  if (sin2phi < 0) sin2phi = 0;
  if (sin2phi > 1) sin2phi = 1;
  sinphi = sqrt(sin2phi);

  statev[10] = asin(sinphi)*180.0/3.141592;
  /* statev(11) 1-based -> statev[10] 0-based; rho set to 0 (visco) */

  statev[11] = 0;

  sensit = 1;
  if (statev[13] >= 1) sensit = statev[13];
  statev[14] = -sensit*exp((parms[4] - log(1 + statev[6]))/parms[2])/pp;

  (void)deps; (void)nasv; (void)nparms;
  (void)norm_del; (void)norm_del2; (void)del;
}

/* ------------------------------------------------------------------ */
/* main entry: single integration-point update (Abaqus UMAT interface) */
/* props[0..28] raw material parameters (degrees for phi_c), props[21] */
/* = e0 or OCR (>10 => OCR = props[21]-10); props[28] = s0             */
/* statev[] has 16 entries, same layout as the Fortran UMAT.           */
/* testing: 0 = normal Abaqus, 1 = PLAXIS first call (loose tol),      */
/*          2 = stiffness-only (norm_D == 0)                           */
/* ------------------------------------------------------------------ */
void masin_visco_umat(double *stress, double *statev, double *ddsdde,
  double *dstran, double dtime, double *props, int nprops, int testing,
  int *error)
{
  double parms[40];
  double sig_n[6], deps_np1[6], depsv_np1;
  double pore, asv[NASV], y[NYDIM], y_n[NYDIM];
  double DDtan[36], F_sig[6];
  double dtsub, tolintT, tolintTtest, DTmin;
  double ameanstress, avoid, aOCR;
  double sensit;
  double theta = 0.0;
  double phi_corig_deg, beta_deg, M, Mnew, sinphiMnew, sinphiM;
  double phi_deg;
  int nparms, nasv, nyact, nfev, error_RKF, inittension;
  int i;
  double pp;

  nparms = nprops;
  for (i = 0; i < nprops && i < 40; i++) parms[i] = props[i];

  /* check_parms: reduce phi_c by the shear-band angle beta */
  phi_corig_deg = parms[0];
  beta_deg = parms[22];
  M = 6.0*sin(phi_corig_deg*PI/180.0)/(3.0 - sin(phi_corig_deg*PI/180.0));
  Mnew = M - tan(beta_deg*PI/180.0);
  sinphiMnew = (3.0*Mnew)/(6.0 + Mnew);
  sinphiM = (3.0*M)/(6.0 + M);
  phi_deg = phi_corig_deg -
    (asin(sinphiM) - asin(sinphiMnew))*180.0/PI;
  parms[0] = phi_deg*PI/180.0; /* phi deg->rad (reduced) */

  tolintT = 1.0e-3;
  tolintTtest = 1.0e-1;
  DTmin = 1.0e-17;

  nasv = NASV; /* 10 */
  nyact = 6 + nasv;

  dtsub = statev[12];
  pore = -statev[7];

  sensit = 1;
  if (statev[13] >= 1) sensit = statev[13];

  if (statev[6] < 0.001) {
    ameanstress = -(stress[0]+stress[1]+stress[2])/3;
    avoid = 0;
    if (props[27] <= 10.0) {
      avoid = props[27];
    } else if (props[27] > 10.0) {
      aOCR = props[27] - 10.0;
      avoid = exp(props[4] - props[2]*log(ameanstress + props[1]) -
        props[2]*log(aOCR/sensit)) - 1;
    }
    statev[6] = avoid;
    statev[15] = 0;
  }

  for (i = 0; i < nasv - 1; i++) asv[i] = statev[i];
  asv[7] = statev[13];  /* sensitivity */

  for (i = 0; i < 6; i++) { sig_n[i] = 0; deps_np1[i] = 0; }
  move_sig(stress, pore, sig_n);
  for (i = 0; i < 6; i++) deps_np1[i] = dstran[i];
  depsv_np1 = deps_np1[0]+deps_np1[1]+deps_np1[2];

  {
    double norm_D2 = dot_vect(2, deps_np1, deps_np1, 6);
    double norm_D = sqrt(norm_D2);
    int testnan = 0;
    umatisnan(norm_D, &testnan);
    if (testnan == 1) { *error = 3; return; }
  }

  iniy(y, sig_n, asv);
  for (i = 0; i < NYDIM; i++) y_n[i] = y[i];

  /* check initial state (inittension) */
  {
    int tmp = 0;
    check_RKF(&tmp, y, nyact, nasv, parms, nparms);
    inittension = tmp;
  }

  if ((dtsub <= 0.0) || (dtsub > dtime)) dtsub = dtime;

  nfev = 0;

  if (inittension == 0) {
    if (testing == 1) {
      rkf23_update(y, nyact, nasv, &dtsub, tolintTtest, 1000, DTmin,
        deps_np1, parms, nparms, &nfev, dtime, error);
      if (*error == 3) {
        for (i = 0; i < nyact; i++) y[i] = y_n[i];
        *error = 0;
      }
    } else if (testing == 2) {
      for (i = 0; i < nyact; i++) y[i] = y_n[i];
    } else {
      rkf23_update(y, nyact, nasv, &dtsub, tolintT, 10000, DTmin,
        deps_np1, parms, nparms, &nfev, dtime, error);
    }

    if (*error == 3) {
      for (i = 0; i < nyact; i++) y[i] = y_n[i];
      statev[15] = 1;
    } else if (*error == 10) {
      *error = 10;
      return;
    }

    perturbate(y_n, nasv, parms, nparms, deps_np1, DDtan, error);

    error_RKF = 0;
    check_RKF(&error_RKF, y, nyact, nasv, parms, nparms);
    if (error_RKF == 1) {
      matmul(DDtan, deps_np1, F_sig, 6, 6, 1);
      for (i = 0; i < 6; i++) y[i] = y_n[i] + F_sig[i];
    }
    error_RKF = 0;
    check_RKF(&error_RKF, y, nyact, nasv, parms, nparms);
    if (error_RKF == 1)
      for (i = 0; i < 6; i++) y[i] = y_n[i];
  } else {
    perturbate(y_n, nasv, parms, nparms, deps_np1, DDtan, error);
    matmul(DDtan, deps_np1, F_sig, 6, 6, 1);
    for (i = 0; i < 6; i++) y[i] = y_n[i] + F_sig[i];
    error_RKF = 0;
    check_RKF(&error_RKF, y, nyact, nasv, parms, nparms);
    if (error_RKF == 1)
      for (i = 0; i < 6; i++) y[i] = y_n[i];
    statev[15] = 1;
  }

  if (dtsub <= 0.0) dtsub = 0;
  else if (dtsub >= dtime) dtsub = dtime;
  statev[12] = dtsub;
  statev[9] = (double)nfev;

  solout(stress, asv, nasv, y, pore, depsv_np1, parms, DDtan, ddsdde);

  for (i = 0; i < nasv - 1; i++) statev[i] = asv[i];
  statev[13] = asv[7];

  for (i = 0; i < 6; i++) sig_n[i] = y[i];
  pp = -(sig_n[0]+sig_n[1]+sig_n[2])/3;

  statev[7] = -pore;
  statev[8] = pp;

  if (inittension == 0) calc_statev(sig_n, statev, parms, nparms, nasv,
    deps_np1);

  (void)theta; (void)ameanstress;
}

/* ------------------------------------------------------------------ */
/* Niemunis visco law (Dr, Iv) — theory from the professional manual   */
/* (section "Niemunis visco law", the extension used by                */
/*  group_materi_plasti_hypo_masin_clay_visco Dr Iv).                  */
/*                                                                     */
/*   stress rate:  sig_dot = M:D - L:eps_vis                          */
/*   L = fb*Lhat,  fb = -sigkk / ((1+a^2/3)*kappa)                    */
/*   Lhat = F^2 I + a^2 shat s_hat + b^2 (I - 1/3 II)                */
/*   eps_vis = Dr * mhat * (1/OCR)^(1/Iv)                              */
/*   m = -(F^2/a^2 (shat + shat*) + (shat:shat) shat*                 */
/*        - shat (shat:shat*))                                         */
/*   OCR = pe / pe+                                                    */
/*                                                                     */
/* Defaults derived from the clay parameters (the manual exposes only  */
/* Dr Iv): kappa=lambda*, ee0=e_initial, pe0=p_initial, betaR=1.       */
/* props: [0]=phi_c, [2]=lambda*, [3]=kappa*(unused, see default),     */
/*        [5]=nu_pp, [17]=vertical, [25]=Dr, [26]=Iv, [27]=e0          */
/* ------------------------------------------------------------------ */
void masin_niemunis_visco_umat(double *stress, double *statev,
  double *ddsdde, double *dstran, double dtime, double *props, int nprops,
  int testing, int *error)
{
  double phi, p_t, lam_star, nu_pp, vertical;
  double Dr, Iv, e0, evoid;
  double sig[6], sig_star[6], sig_rot_tmp[6], deps[6];
  double I1, I2, I3, pp, qq, cos3t;
  double F, F2, a, a2, b2, fb;
  double sig_hat[6], sig_dev[6], pmean, s2;
  double m[6], mnorm, mhat[6];
  double pe, pe0, peplus, eta, M, betaR, kappa, lambda;
  double OCR, creep_rate, Lhat[36], L[36], Ldeps[6];
  int i, j;
  double vnu;

  (void)testing; (void)nprops;

  phi = props[0]*PI/180.0;
  p_t = props[1];
  lam_star = props[2];
  nu_pp = props[5];
  vertical = props[17];
  Dr = props[25];
  Iv = props[26];
  e0 = props[27];
  evoid = statev[6];

  /* defaults from the clay parameters */
  lambda = lam_star;
  kappa = lam_star;   /* kappa not exposed by the clay; use lambda* */
  betaR = 1.0;
  pe0 = -(stress[0]+stress[1]+stress[2])/3.0;
  if (pe0 <= 0) pe0 = 1.0;
  if (evoid <= 0) evoid = e0;

  /* current stress (effective), axis-shifted */
  sig[0]=stress[0]; sig[1]=stress[1]; sig[2]=stress[2];
  sig[3]=stress[3]; sig[4]=stress[4]; sig[5]=stress[5];
  sig_star[0]=sig[0]-p_t; sig_star[1]=sig[1]-p_t; sig_star[2]=sig[2]-p_t;
  sig_star[3]=sig[3]; sig_star[4]=sig[4]; sig_star[5]=sig[5];

  inv_sig(sig_star, &pp, &qq, &cos3t, &cos3t, &I1, &I2, &I3,
    &I1, &I2, &I3, sig_rot_tmp, props);
  /* NOTE: inv_sig is the rotated version; for the Niemunis law we
     use the non-rotated invariants I1,I2,I3 (rotated outputs ignored). */

  pmean = -I1/3.0;
  if (pmean < 1.0e-10) { *error = 10; return; }

  /* Matsuoka-Nakai F */
  if ((I3 + I1*I2) != 0) {
    double sin2phim = (9*I3 + I1*I2)/(I3 + I1*I2);
    if (sin2phim > 1) sin2phim = 1;
    if (sin2phim < 0) sin2phim = 0;
    F = sqrt(sin2phim);
  } else {
    F = 1.0;
  }
  F2 = F*F;

  a = sqrt(3.0)*(3.0 - sin(phi))/(2.0*sqrt(2.0)*sin(phi));
  a2 = a*a;
  b2 = (1.0 + a2/3.0)*(1.0 - 2.0*nu_pp)/(1.0 + nu_pp) - 1.0;
  if (b2 < 0) b2 = 0.0;

  fb = -I1/((1.0 + a2/3.0)*kappa);

  /* normalised stress shat = sig_star / I1 */
  for (i = 0; i < 6; i++) sig_hat[i] = sig_star[i]/I1;
  sig_dev[0]=sig_hat[0]-1.0/3.0; sig_dev[1]=sig_hat[1]-1.0/3.0;
  sig_dev[2]=sig_hat[2]-1.0/3.0;
  sig_dev[3]=sig_hat[3]; sig_dev[4]=sig_hat[4]; sig_dev[5]=sig_hat[5];

  /* Lhat = F^2 I + a^2 shat shat + b^2 (I - 1/3 II)  in Voigt */
  for (i = 0; i < 36; i++) Lhat[i] = 0.0;
  Lhat[0*6+0]=1; Lhat[1*6+1]=1; Lhat[2*6+2]=1;
  Lhat[3*6+3]=1; Lhat[4*6+4]=1; Lhat[5*6+5]=1;
  for (i = 0; i < 6; i++)
    for (j = 0; j < 6; j++) {
      double Kron[6];
      Kron[0]=1; Kron[1]=1; Kron[2]=1; Kron[3]=0; Kron[4]=0; Kron[5]=0;
      Lhat[i*6+j] = F2*Lhat[i*6+j] +
        a2*sig_hat[i]*sig_hat[j] +
        b2*( (i==j?1.0:0.0) - (1.0/3.0)*Kron[i]*Kron[j] );
    }

  /* L = fb * Lhat */
  for (i = 0; i < 36; i++) L[i] = fb*Lhat[i];

  /* flow rule m */
  s2 = 0.0;
  for (i = 0; i < 6; i++) s2 += sig_dev[i]*sig_dev[i];
  for (i = 0; i < 6; i++)
    m[i] = -(F2/a2)*(sig_hat[i] + sig_dev[i]) + s2*sig_dev[i] -
      sig_hat[i]*s2;

  mnorm = 0.0;
  for (i = 0; i < 6; i++) mnorm += m[i]*m[i];
  if (mnorm > 1.0e-12) mnorm = sqrt(mnorm);
  else mnorm = 1.0;
  for (i = 0; i < 6; i++) mhat[i] = m[i]/mnorm;

  /* OCR = pe / pe+  (Niemunis) */
  pe = pe0*exp((1.0/lambda)*log((1.0+e0)/(1.0+evoid)));
  M = 6.0*F*sin(phi)/(3.0 - sin(phi));
  eta = qq/(M*pmean);
  if (eta < 1.0) {
    peplus = pmean/(betaR - 1.0)*
      (betaR*sqrt(1.0 + eta*eta*(betaR*betaR - 1.0)) - 1.0);
  } else {
    peplus = pmean*pow(1.0 + eta*eta, (1.0 + betaR)/2.0);
  }
  if (peplus > 1.0e-12) OCR = pe/peplus;
  else OCR = 1.0;
  if (OCR <= 0) OCR = 1.0e-6;

  /* creep rate; clamp the exponent to keep the model stable. In the
     Niemunis theory OCR>=1 in normal use (1/OCR <= 1); the clamp only
     protects against pathological states from the simplified defaults. */
  creep_rate = Dr*pow(1.0/OCR, 1.0/Iv);
  if (creep_rate > 1.0e3*Dr) creep_rate = 1.0e3*Dr;

  /* sig_dot = L:D - L:eps_vis ; eps_vis = creep_rate * mhat */
  for (i = 0; i < 6; i++) {
    deps[i] = dstran[i]/dtime;   /* strain rate */
    Ldeps[i] = 0.0;
    for (j = 0; j < 6; j++) Ldeps[i] += L[i*6+j]*deps[j];
  }
  for (i = 0; i < 6; i++) {
    double Lm = 0.0;
    for (j = 0; j < 6; j++) Lm += L[i*6+j]*mhat[j];
    stress[i] += (Ldeps[i] - Lm*creep_rate)*dtime;
  }

  /* void ratio evolution: de = (1+e)*tr(D) - visco volumetric part */
  {
    double trD = (deps[0]+deps[1]+deps[2])*dtime;
    double trm_vis = (mhat[0]+mhat[1]+mhat[2])*creep_rate*dtime;
    evoid += (1.0+evoid)*(trD - trm_vis);
  }
  statev[6] = evoid;

  /* tangent: use L (elastic-like), scaled by (1 - creep sensitivity) */
  for (i = 0; i < 36; i++) ddsdde[i] = L[i];

  (void)vertical;
}
