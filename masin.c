/*
    masin.c - C port of umat_hcea.for
    Clay hypoplasticity model of Masin (Masin 2014, Geotechnique 64(3):232-238)
    including stiffness anisotropy, intergranular strain (Niemunis-Herle),
    structure (sensitivity) and adaptive RKF23 substepping.

    Original Fortran: C. Tamagnini, E. Sellari, D. Masin, P.A. von Wolffersdorff
    (GPL). This file is a faithful C port of the reference implementation
    umat_hcea.for (2885 lines). Function names map 1:1 to the Fortran
    subroutines; indices are converted from Fortran 1-based to C 0-based.

    Conventions (identical to the Fortran):
      - Voigt notation: [11,22,33,12,13,23]
      - tension/extension positive
      - dot_vect flag: 1=stress-like (shear weight 2), 2=strain-like (shear
        weight 0.5), 3=ordinary
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NASV    8
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
/* invariants of a stress tensor (Voigt)                               */
/* ------------------------------------------------------------------ */
static void inv_sig(double *sig, double *pp, double *qq, double *cos3t,
  double *I1, double *I2, double *I3)
{
  double sdev[6], eta[6], eta_d[6], eta_d2[6];
  double xmin1, xmin2, xmin3, tretadev3;
  double norm2, norm2sig, norm2eta, numer, denom;
  double half = 0.5, one = 1.0, three = 3.0;
  double onethird = 1.0 / 3.0, threehalves = 3.0 / 2.0, sqrt6 = sqrt(6.0);
  double tiny = 1.0e-18;
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

  (void)three;
}

/* ------------------------------------------------------------------ */
/* tangent operators MM, HH, LL, NN for the Masin model                */
/* istrain=1: full model with intergranular strains                    */
/* istrain=0: basic model without intergranular strains                */
/* ------------------------------------------------------------------ */
static void get_tan(double *deps, double *sig, double *q, int nasv,
  double *parms, int nparms, double *MM, double *HH, double *LL,
  double *NN, int istrain, int *error)
{
  double eta[6], eta_dev[6], del[6], evoid, sig_star[6], sensit;
  double H_s[6], eta_del[6], eta_delta[6], eta_eps[6];
  double norm_del, norm_del2, norm_deps, norm_deps2;
  double pp, qq, cos3t, I1, I2, I3;
  double a, a2, fd, fs, fdsbs, fddivfdA;
  double AA[36], H_del[36], H_e[6], IU[36], Leta[6], hypo_Dsom[6];
  double krondelta[6], AAhce[36], NNvec[6];
  double load, rho, N_par, Stf, kparam, Aparam, sfparam;
  double kap_par, lambda_struct;
  double zero=0.0, one=1.0, two=2.0, three=3.0;
  double tiny=1.0e-17, half=0.5;
  double onethird, sqrt3, twosqrt2, sqrt2, ln2m1;
  double temp1, temp2, temp3, temp4, gamma;
  double phi, lam_star, kap_star, N_star, nuvh, r_uc;
  double m_R, m_T, beta_r, chi, p_t, sinphi, sinphi2;
  double nparam, m_Trat, G0, Gvh, alphanu, ocrcs;
  double alphaG, alphaE, Am, nuhh, pmeangt1;
  double sin2phim, ashape, npow, cos2phic, alpha_power;
  double pmean, kpow, sinphickpow, sinphimkpow, Amult, peast;
  double nvect[3], pmat[9], kron_delta[9];
  double kck[81], kdk[81], pdk[81], kdp[81], pck[81], pdp[81], LLfour[81];
  double an1, an2, an3, an4, an5, norm_m, norm_m2;
  int i, j, k, l, softmodel, vert;
  double t1;

  onethird = one/three;
  sqrt3 = sqrt(three);
  twosqrt2 = two*sqrt(two);
  sqrt2 = sqrt(two);
  ln2m1 = one/log(two);

  for (i = 0; i < 6; i++)
    for (j = 0; j < 6; j++) {
      MM[i*6+j] = zero; LL[i*6+j] = zero; IU[i*6+j] = zero; H_del[i*6+j] = zero;
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
  (void)Aparam;
  sfparam = parms[9];
  r_uc = parms[10];
  beta_r = parms[11];
  chi = parms[12];
  gamma = chi;
  if (gamma == 0) gamma = chi;
  G0 = parms[13];
  nparam = parms[14];
  m_Trat = parms[15];
  vert = (int)parms[17];
  alphaE = parms[18];
  alphanu = parms[19];
  if (alphaE < 0.01) alphaE = pow(alphaG, 1.25);
  if (alphanu < 0.01) alphanu = alphaG;

  sinphi = sin(phi);
  sinphi2 = sinphi*sinphi;
  (void)sinphi2;
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
  sig_star[0] = sig[0]-p_t; sig_star[1] = sig[1]-p_t; sig_star[2] = sig[2]-p_t;
  sig_star[3] = sig[3]; sig_star[4] = sig[4]; sig_star[5] = sig[5];

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

  /* auxiliary stress tensors */
  inv_sig(sig_star, &pp, &qq, &cos3t, &I1, &I2, &I3);

  for (i = 0; i < 6; i++) eta[i] = sig_star[i]/I1;

  eta_dev[0]=eta[0]-onethird; eta_dev[1]=eta[1]-onethird; eta_dev[2]=eta[2]-onethird;
  eta_dev[3]=eta[3]; eta_dev[4]=eta[4]; eta_dev[5]=eta[5];

  krondelta[0]=one; krondelta[1]=one; krondelta[2]=one;
  krondelta[3]=zero; krondelta[4]=zero; krondelta[5]=zero;

  /* explicit clay hypoplasticity specific */
  peast = exp((N_star - log(one+evoid))/lam_star);

  if ((I3 + I1*I2) != 0) sin2phim = (9*I3 + I1*I2)/(I3 + I1*I2);
  else sin2phim = 1;
  if (sin2phim > 1) sin2phim = 1;
  if (sin2phim < 0) sin2phim = 0;

  cos2phic = 1 - sin(phi)*sin(phi);
  ashape = parms[22];          /* ay, default 0.30 */
  if (ashape < 1.0e-6) ashape = 0.30;
  ocrcs = parms[23];           /* oc, default 2.0 */
  if (ocrcs < 1.0e-6) ocrcs = 2.0;
  npow = -log(cos2phic)/log(ocrcs) + ashape*(sin2phim - sin(phi)*sin(phi));
  fdsbs = ocrcs*pow(1 - sin2phim, 1/npow);

  pmean = -I1/3.0;
  fd = pow((ocrcs*pmean)/peast, alpha_power);
  fdsbs = pow(fdsbs, alpha_power);
  fddivfdA = fd/fdsbs;

  if (sin2phim < 1.e-10) { sin2phim = 0; cos3t = -1; }

  Am = nuvh*nuvh*(4*alphaE*alphanu - 2*alphaE*alphaE*alphanu*alphanu +
        2*alphaE*alphaE - alphanu*alphanu) +
    nuvh*(4*alphaE + 2*alphaE*alphanu) + 1 + 2*alphaE;

  fs = 9*pmean/2*(1/kap_par + 1/lambda_struct)/Am;

  for (i = 0; i < 3; i++) nvect[i] = 0;
  if (vert >= 1 && vert <= 3) nvect[vert-1] = 1;
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

  for (i = 0; i < 6; i++) hypo_Dsom[i] = 0;

  kpow = 1.70 + 3.90*sin(phi)*sin(phi);
  sinphickpow = pow(sin(phi), kpow);
  sinphimkpow = pow(sqrt(sin2phim), kpow);

  Amult = 2.0/3.0 - sqrt(sqrt(sin2phim))*(cos3t + 1.0)/4.0;

  for (i = 0; i < 6; i++)
    hypo_Dsom[i] = -eta_dev[i] + krondelta[i]*
      (sinphimkpow - sinphickpow)/(1 - sinphickpow)*Amult;
  for (i = 3; i < 6; i++) hypo_Dsom[i] = 2.0*hypo_Dsom[i];

  norm_m2 = dot_vect(2, hypo_Dsom, hypo_Dsom, 6);
  norm_m = sqrt(norm_m2);
  for (i = 0; i < 6; i++) hypo_Dsom[i] = hypo_Dsom[i]/norm_m;

  for (i = 0; i < 6; i++)
    for (j = 0; j < 6; j++)
      AAhce[i*6+j] = fs*LL[i*6+j] + sig_star[i]*krondelta[j]/lambda_struct;

  matmul(AAhce, hypo_Dsom, NNvec, 6, 6, 1);
  for (i = 0; i < 6; i++) NN[i] = -NNvec[i]*fddivfdA/(fd*fs);

  /* end basic model */
  if (istrain == 1) {
    load = dot_vect(2, eta_del, eta_eps, 6);

    rho = norm_del/r_uc;
    if (rho > one) rho = one;

    matmul(LL, eta_del, Leta, 6, 6, 1);

    pmeangt1 = p_t/4.0;
    if (pmean > pmeangt1) pmeangt1 = pmean;

    Gvh = G0*pow(pmeangt1, nparam);
    m_R = Gvh*4*Am*alphaG/(9*pmeangt1*alphaE)*lambda_struct*kap_star/
      (lambda_struct + kap_star)/(1 - alphanu*nuvh - 2*alphaE*nuvh*nuvh);
    if (m_R > 40) m_R = 40;
    m_T = m_R*m_Trat;
    temp1 = ((pow(rho,chi)*m_T + (1 - pow(rho,chi))*m_R)*fs);

    if (load > zero) {
      temp2 = pow(rho,chi)*(one - m_T)*fs;
      temp3 = pow(rho,gamma)*fs*fd;
      for (i = 0; i < 6; i++)
        for (j = 0; j < 6; j++) {
          AA[i*6+j] = temp2*Leta[i]*eta_delta[j] + temp3*NN[i]*eta_delta[j];
          MM[i*6+j] = temp1*LL[i*6+j] + AA[i*6+j];
        }
    } else {
      temp4 = pow(rho,chi)*(m_R - m_T)*fs;
      for (i = 0; i < 6; i++)
        for (j = 0; j < 6; j++) {
          AA[i*6+j] = temp4*Leta[i]*eta_delta[j];
          MM[i*6+j] = temp1*LL[i*6+j] + AA[i*6+j];
        }
    }

    if (load > zero) {
      for (i = 0; i < 6; i++)
        for (j = 0; j < 6; j++)
          H_del[i*6+j] = IU[i*6+j] - pow(rho,beta_r)*eta_del[i]*eta_delta[j];
    } else {
      for (i = 0; i < 6; i++) H_del[i*6+i] = one;
    }
  }

  for (i = 0; i < 6; i++) H_e[i] = (i < 3) ? (one+evoid) : zero;

  for (i = 0; i < 6; i++) {
    H_s[i] = zero;
    if (softmodel == 1) {
      if (i < 3) H_s[i] = -kparam*(sensit - sfparam)/lam_star;
    }
  }

  for (i = 0; i < nasv; i++) {
    if (i < 6) {
      if (istrain == 1)
        for (j = 0; j < 6; j++) HH[i*6+j] = H_del[i*6+j];
    } else if (i == 6) {
      for (j = 0; j < 6; j++) HH[i*6+j] = H_e[j];
    } else if (i == 7) {
      for (j = 0; j < 6; j++) HH[i*6+j] = H_s[j];
    }
  }

  for (i = 0; i < 6; i++)
    for (j = 0; j < 6; j++) LL[i*6+j] = LL[i*6+j]*fs;
  for (i = 0; i < 6; i++) NN[i] = NN[i]*fs*fd;

  (void)t1; (void)sqrt2; (void)pmean;
}

/* ------------------------------------------------------------------ */
/* F_sig = MM*deps (istrain=1) or LL*deps+NN*|deps| (istrain=0)        */
/* F_q   = HH*deps with soft-model correction                          */
/* ------------------------------------------------------------------ */
static void get_F_sig_q(double *sig, double *q, int nasv, double *parms,
  int nparms, double *deps, double *F_sig, double *F_q, int *error)
{
  double MM[36], HH[NASV*6], LL[36], NN[6], norm_D, norm_D2, norm2;
  double trD, depsh, depd, Aparam, edev[6];
  int istrain, ii;

  if (parms[13] <= 0.5) istrain = 0;
  else istrain = 1;

  get_tan(deps, sig, q, nasv, parms, nparms, MM, HH, LL, NN, istrain, error);

  if (istrain == 1) {
    matmul(MM, deps, F_sig, 6, 6, 1);
  } else {
    matmul(LL, deps, F_sig, 6, 6, 1);
    norm_D2 = dot_vect(2, deps, deps, 6);
    norm_D = sqrt(norm_D2);
    for (ii = 0; ii < 6; ii++) F_sig[ii] = F_sig[ii] + NN[ii]*norm_D;
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

/* ------------------------------------------------------------------ */
/* RHS of the ODE system y' = f(y)                                     */
/* ------------------------------------------------------------------ */
static void rhs(double *y, int ny, int nasv, double *parms, int nparms,
  double *deps, double *kRK, int *nfev, int *error)
{
  double sig[6], q[NASV], F_sig[6], F_q[NASV];
  int i;

  (*nfev)++;
  for (i = 0; i < ny; i++) kRK[i] = 0.0;

  for (i = 0; i < 6; i++) sig[i] = y[i];
  for (i = 0; i < nasv; i++) q[i] = y[6+i];

  get_F_sig_q(sig, q, nasv, parms, nparms, deps, F_sig, F_q, error);
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
  (void)nasv;
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
  double I1, I2, I3, pp, qq, cos3t;
  double p_t, minstress, sin2phim, tolerance;
  double OCR, omega, fSBS, sensit, cos2phic;
  double ashape, ocrcs;
  int i, testnan;

  p_t = parms[1];
  minstress = p_t/100.0;
  for (i = 0; i < 6; i++) sig[i] = y[i];

  sig_star[0]=sig[0]-p_t; sig_star[1]=sig[1]-p_t; sig_star[2]=sig[2]-p_t;
  sig_star[3]=sig[3]; sig_star[4]=sig[4]; sig_star[5]=sig[5];

  inv_sig(sig_star, &pp, &qq, &cos3t, &I1, &I2, &I3);
  pmean = -I1/3;
  if (pmean < minstress) *error_RKF = 1;

  if ((I3 + I1*I2) != 0) sin2phim = (9*I3 + I1*I2)/(I3 + I1*I2);
  else sin2phim = 1;

  sensit = 1;
  if (y[6+7] >= 1) sensit = y[6+7];
  OCR = -sensit*exp((parms[4] - log(1 + y[6+6]))/parms[2])/pp;
  cos2phic = 1 - sin(parms[0])*sin(parms[0]);

  ashape = parms[22];
  if (ashape < 1.0e-6) ashape = 0.3;
  ocrcs = parms[23];
  if (ocrcs < 1.0e-6) ocrcs = 2.;
  omega = -log(cos2phic)/log(ocrcs) + ashape*(sin2phim - sin(parms[0])*sin(parms[0]));
  fSBS = sin2phim + pow(1/OCR, omega) - 1;

  if (sin2phim >= 1) *error_RKF = 1;

  tolerance = 0.1;
  if (parms[13] >= 0.5) {  /* istrain active */
    if (OCR <= ocrcs) {
      if (sin2phim >= (1 - cos2phic + tolerance)) *error_RKF = 1;
    }
    if (OCR > ocrcs) {
      if (fSBS >= tolerance) *error_RKF = 1;
    }
  } else {
    if (fSBS >= tolerance) *error_RKF = 1;
  }

  testnan = 0;
  for (i = 0; i < ny; i++) umatisnan(y[i], &testnan);
  umatisnan(sin2phim, &testnan);
  umatisnan(OCR, &testnan);
  umatisnan(fSBS, &testnan);
  if (testnan == 1) *error_RKF = 1;

  (void)nasv; (void)nparms;
}

/* ------------------------------------------------------------------ */
/* adaptive RKF23 integration with substepping                         */
/* ------------------------------------------------------------------ */
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
    else rhs(y_k, n, nasv, parms, nparms, deps_np1, kRK_1, nfev, error);
    if (*error == 10) return;

    temp = half*DT_k;
    for (i = 0; i < n; i++) y_2[i] = y_k[i] + temp*kRK_1[i];

    error_RKF = 0;
    check_RKF(&error_RKF, y_2, n, nasv, parms, nparms);
    if (error_RKF == 1) { *error = 3; return; }
    else rhs(y_2, n, nasv, parms, nparms, deps_np1, kRK_2, nfev, error);
    if (*error == 10) return;

    for (i = 0; i < n; i++)
      y_3[i] = y_k[i] - DT_k*kRK_1[i] + two*DT_k*kRK_2[i];

    error_RKF = 0;
    check_RKF(&error_RKF, y_3, n, nasv, parms, nparms);
    if (error_RKF == 1) { *error = 3; return; }
    else rhs(y_3, n, nasv, parms, nparms, deps_np1, kRK_3, nfev, error);
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
  double sig[6], q[NASV];
  double HHtmp[NASV*6], LL[36], NN[6];
  double m_R, Gvh, pmean, G0, alphaG, alphaE, alphanu, nuvh, nuhh;
  double lambda, kappa, Am, pmeangt1, p_t;
  int istrain, jj, kk;

  G0 = parms[13];
  p_t = parms[1];
  for (kk = 0; kk < 6; kk++) sig[kk] = y_n[kk];

  if (-(sig[0]+sig[1]+sig[2])/3 < p_t) {
    sig[0]=-p_t; sig[1]=-p_t; sig[2]=-p_t;
    sig[3]=0; sig[4]=0; sig[5]=0;
  }

  for (kk = 0; kk < nasv; kk++) q[kk] = y_n[6+kk];

  pmean = -(sig[0]+sig[1]+sig[2])/3;

  nuhh = parms[5];
  lambda = parms[2];
  kappa = parms[3];
  alphaG = parms[6];
  if (alphaG < 0.01) alphaG = 1.0;
  alphaE = parms[18];
  alphanu = parms[19];
  if (alphaE < 0.01) alphaE = pow(alphaG, 1.25);
  if (alphanu < 0.01) alphanu = alphaG;
  nuvh = nuhh/alphanu;

  Am = nuvh*nuvh*(4*alphaE*alphanu - 2*alphaE*alphaE*alphanu*alphanu +
        2*alphaE*alphaE - alphanu*alphanu) +
    nuvh*(4*alphaE + 2*alphaE*alphanu) + 1 + 2*alphaE;
  pmeangt1 = p_t/4.0;
  if (pmean > pmeangt1) pmeangt1 = pmean;
  Gvh = G0*pow(pmeangt1, parms[14]);
  m_R = Gvh*4*Am*alphaG/(9*pmeangt1*alphaE)*lambda*kappa/(lambda+kappa)/
    (1 - alphanu*nuvh - 2*alphaE*nuvh*nuvh);
  if (m_R > 40) m_R = 40;

  if (G0 <= 0.5) istrain = 0;
  else istrain = 1;

  for (kk = 0; kk < 6; kk++)
    for (jj = 0; jj < 6; jj++) DD[kk*6+jj] = 0.0;

  if (*error != 10) {
    get_tan(deps_np1, sig, q, nasv, parms, nparms, DD, HHtmp, LL, NN,
      istrain, error);
  }
  if (istrain == 0) {
    for (kk = 0; kk < 6; kk++)
      for (jj = 0; jj < 6; jj++) DD[kk*6+jj] = LL[kk*6+jj];
  } else {
    for (kk = 0; kk < 6; kk++)
      for (jj = 0; jj < 6; jj++) DD[kk*6+jj] = m_R*LL[kk*6+jj];
  }

  (void)nparms;
}

/* ------------------------------------------------------------------ */
/* helpers used by masin_umat (mapped from Fortran)                    */
/* ------------------------------------------------------------------ */

/* sig = stress + pore (effective from total) */
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
  int i;

  p_t = parms[1];
  for (i = 0; i < 3; i++) sig_star[i] = stress[i]-p_t;
  for (i = 3; i < 6; i++) sig_star[i] = stress[i];

  inv_sig(sig_star, &pp, &qq, &cos3t, &I1, &I2, &I3);
  if (I3 != 0) sin2phi = (9. + I1*I2/I3)/(1. + I1*I2/I3);
  else sin2phi = 0;
  if (sin2phi < 0) sin2phi = 0;
  if (sin2phi > 1) sin2phi = 1;
  sinphi = sqrt(sin2phi);

  statev[10] = asin(sinphi)*180.0/3.141592;
  /* statev(11) 1-based -> statev[10] 0-based */

  if (parms[13] > 0.5) {  /* istrain active */
    for (i = 0; i < 6; i++) del[i] = statev[i];
    norm_del2 = dot_vect(2, del, del, 6);
    norm_del = sqrt(norm_del2);
    statev[11] = norm_del/parms[10];
  } else {
    statev[11] = 0;
  }

  sensit = 1;
  if (statev[13] >= 1) sensit = statev[13];
  statev[14] = -sensit*exp((parms[4] - log(1 + statev[6]))/parms[2])/pp;

  (void)deps; (void)nasv; (void)nparms;
}

/* ------------------------------------------------------------------ */
/* main entry: single integration-point update (Abaqus UMAT interface) */
/* props[0..28] raw material parameters (degrees for phi_c), props[21] */
/* = e0 or OCR (>10 => OCR = props[21]-10); props[28] = s0             */
/* statev[] has 16 entries, same layout as the Fortran UMAT.           */
/* testing: 0 = normal Abaqus, 1 = PLAXIS first call (loose tol),      */
/*          2 = stiffness-only (norm_D == 0)                           */
/* ------------------------------------------------------------------ */
void masin_umat(double *stress, double *statev, double *ddsdde,
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
  int nparms, nasv, nyact, nfev, error_RKF, inittension;
  int i;
  double pp;

  nparms = nprops;
  for (i = 0; i < nprops && i < 40; i++) parms[i] = props[i];
  parms[0] = parms[0]*PI/180.0; /* phi deg->rad */

  tolintT = 1.0e-3;
  tolintTtest = 1.0e-1;
  DTmin = 1.0e-17;

  nasv = NASV; /* 8 */
  nyact = 6 + nasv;

  dtsub = statev[12];
  pore = -statev[7];

  sensit = 1;
  if (statev[13] >= 1) sensit = statev[13];

  if (statev[6] < 0.001) {
    ameanstress = -(stress[0]+stress[1]+stress[2])/3;
    avoid = 0;
    if (props[21] <= 10.0) {
      avoid = props[21];
    } else if (props[21] > 10.0) {
      aOCR = props[21] - 10.0;
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
