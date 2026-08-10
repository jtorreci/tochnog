/*
    sanisand.c - C port of the SANISAND model (Dafalias & Manzari 2004,
    "Simple plasticity sand model accounting for fabric change effects",
    J. Engng. Mechanics ASCE 130(6):622-634).

    Reference Fortran UMAT (GPL, Martinelli/Miriano/Tamagnini, adapted by
    Charles University):
      validation-suite/reference-sanisand/umat.for  (4899 lines)

    This file is a faithful pure-C port. Conventions identical to the
    Fortran:
      - Voigt notation [11,22,33,12,13,23]
      - SOIL MECHANICS convention: compression positive
      - dot_vect flag: 1=stress-like (shear weight 2), 2=strain-like
        (shear weight 0.5), 3=ordinary
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NYDIM   20      /* 6 sig + 14 internal (alpha 6 + void + Fab 6 + 1) */
#define NZDIM   14      /* alpha_sr 6 + unused */
#define NASV_HH 14


static void pp_kk_set(const double *y, double *pp_kk);
#define nz_unused 0

static const double PI_SANISAND = 3.14159265358979323846264338327950288;

/* ------------------------------------------------------------------ */
/* zero a vector / matrix                                              */
/* ------------------------------------------------------------------ */
static void pzero(double *a, int n)
{
  int i;
  for (i = 0; i < n; i++) a[i] = 0.0;
}

/* ------------------------------------------------------------------ */
/* copy a -> b                                                        */
/* ------------------------------------------------------------------ */
static void push(const double *a, double *b, int n)
{
  int i;
  for (i = 0; i < n; i++) b[i] = a[i];
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
/* deviator of tensor t -> s, trace and mean                          */
/* ------------------------------------------------------------------ */
static void deviator(const double *t, double *s, double *trace, double *mean)
{
  double onethird = 1.0 / 3.0;
  *trace = t[0] + t[1] + t[2];
  *mean = onethird * (*trace);
  s[0] = t[0] - *mean;
  s[1] = t[1] - *mean;
  s[2] = t[2] - *mean;
  s[3] = t[3];
  s[4] = t[4];
  s[5] = t[5];
}

/* ------------------------------------------------------------------ */
/* invariants of a stress tensor (Voigt), soil-mechanics convention    */
/* ------------------------------------------------------------------ */
static void inv_sig_full(const double *sig, double *pp, double *qq,
  double *cos3t, double *I1, double *I2, double *I3)
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
/* check and convert material parameters (M_c, M_e)                    */
/* ------------------------------------------------------------------ */
static void check_parms_DM(const double *props, double *parms, int nprops)
{
  double sinphi, sinphiext;
  double PI = PI_SANISAND;
  int i;

  for (i = 0; i < nprops; i++) parms[i] = props[i];

  if (parms[4] > 5) {                 /* M_c as friction angle in degrees */
    sinphi = sin(parms[4] / 180.0 * PI);
    parms[4] = 6.0 * sinphi / (3.0 - sinphi);
  } else {
    sinphi = 3.0 * parms[4] / (6.0 + parms[4]);
  }
  if (parms[5] > 5) {                 /* M_e as friction angle in degrees */
    sinphiext = sin(parms[5] / 180.0 * PI);
    parms[5] = 6.0 * sinphiext / (3.0 + sinphiext);
  } else if ((parms[5] <= 5) && (parms[5] > 0.01)) {
    sinphiext = 3.0 * parms[5] / (6.0 - parms[5]);
  } else {
    parms[5] = parms[4] * (3.0 - sinphi) / (3.0 + sinphi);
  }
}

/* ------------------------------------------------------------------ */
/* state parameter psi (pyknotropy factor)                             */
/* ------------------------------------------------------------------ */
static double psi_void_DM(double void_ratio, double p, const double *parms,
  int nparms)
{
  double p_a, e0, lambda, xi, ec;
  (void)nparms;
  p_a = parms[0];
  e0 = parms[1];
  lambda = parms[2];
  xi = parms[3];
  ec = e0 - lambda * pow(p / p_a, xi);
  return void_ratio - ec;
}

/* ------------------------------------------------------------------ */
/* alpha_c (flag 1), alpha_b (flag 2), alpha_d (flag 3)               */
/* ------------------------------------------------------------------ */
static void alpha_th_DM(int flag, const double *n, double gth, double psi,
  const double *parms, int nparms, double *alpha)
{
  double M_c, mm, n_b, n_d, M, alpha_th, sqrt23;
  int i;
  (void)nparms;
  sqrt23 = sqrt(2.0 / 3.0);
  M_c = parms[4];
  mm = parms[6];
  n_b = parms[11];
  n_d = parms[13];

  if (flag == 1) {
    M = M_c;
  } else if (flag == 2) {
    M = M_c * exp(-n_b * psi);
  } else {
    M = M_c * exp(n_d * psi);
  }
  alpha_th = M * gth - mm;
  for (i = 0; i < 6; i++) alpha[i] = sqrt23 * alpha_th * n[i];
}

/* ------------------------------------------------------------------ */
/* yield function yf = |s - p alpha| - sqrt(2/3) m p                  */
/* ------------------------------------------------------------------ */
static double yf_DM(const double *y, int ny, const double *parms, int nparms)
{
  double mm, sqrt23, norm2;
  double sig[6], s[6], trace, p, alpha[6], sbar[6];
  int i;
  (void)nparms;
  sqrt23 = sqrt(2.0 / 3.0);
  mm = parms[6];
  for (i = 0; i < 6; i++) {
    sig[i] = y[i];
    alpha[i] = y[6 + i];
  }
  deviator(sig, s, &trace, &p);
  for (i = 0; i < 6; i++) sbar[i] = s[i] - p * alpha[i];
  norm2 = dot_vect(1, sbar, sbar, 6);
  return sqrt(norm2) - sqrt23 * mm * p;
}

/* ------------------------------------------------------------------ */
/* elastic stiffness DDe                                              */
/* ------------------------------------------------------------------ */
static void el_stiff_DM(const double *y, int n, const double *parms,
  int nparms, double *DDe, int *error, double tol_f, int check_ff,
  int drcor, double p_thres, int plastic)
{
  double p_a, G0, nu, ratio;
  double sig1, sig2, sig3, p, void_ratio, pp;
  double coeff1, coeff2, Kt, Gt, fe;
  double Id[36], IxI[36];
  double p_thres_E = 0.001;
  int i, j;
  (void)n; (void)error; (void)tol_f; (void)check_ff; (void)drcor; (void)plastic;

  pzero(Id, 36);
  pzero(IxI, 36);
  pzero(DDe, 36);

  Id[0*6+0]=1.0; Id[1*6+1]=1.0; Id[2*6+2]=1.0;
  Id[3*6+3]=0.5; Id[4*6+4]=0.5; Id[5*6+5]=0.5;

  IxI[0*6+0]=1.0; IxI[1*6+0]=1.0; IxI[2*6+0]=1.0;
  IxI[0*6+1]=1.0; IxI[1*6+1]=1.0; IxI[2*6+1]=1.0;
  IxI[0*6+2]=1.0; IxI[1*6+2]=1.0; IxI[2*6+2]=1.0;

  p_a = parms[0];
  G0 = parms[7];
  nu = parms[8];

  sig1 = y[0]; sig2 = y[1]; sig3 = y[2];
  void_ratio = y[12];
  p = (sig1 + sig2 + sig3) / 3.0;
  pp = p;
  (void)p_thres_E;
  if (p < p_thres) pp = p_thres;

  ratio = 3.0 * (1.0 - 2.0 * nu) / (2.0 * (1.0 + nu));
  fe = (2.97 - void_ratio) * (2.97 - void_ratio) / (1.0 + void_ratio);
  Gt = G0 * p_a * fe * sqrt(pp / p_a);
  Kt = Gt / ratio;

  coeff1 = Kt - 2.0 * Gt / 3.0;
  coeff2 = 2.0 * Gt;

  for (i = 0; i < 6; i++)
    for (j = 0; j < 6; j++)
      DDe[i*6+j] = coeff1 * IxI[i*6+j] + coeff2 * Id[i*6+j];
}

/* ------------------------------------------------------------------ */
/* F(y) for (hypo)elastic processes                                    */
/* ------------------------------------------------------------------ */
static void f_hypoelas_DM(const double *y, int n, const double *parms,
  int nparms, const double *deps, double *F, int *error, double tol_f,
  int check_ff, int drcor, double p_thres, int plastic)
{
  double depsv, void_ratio, one = 1.0;
  double De[36], dsig_e[6];

  pzero(F, n);
  void_ratio = y[12];
  depsv = deps[0] + deps[1] + deps[2];

  el_stiff_DM(y, n, parms, nparms, De, error, tol_f, check_ff, drcor,
    p_thres, plastic);
  matmul(De, deps, dsig_e, 6, 6, 1);

  F[0]=dsig_e[0]; F[1]=dsig_e[1]; F[2]=dsig_e[2];
  F[3]=dsig_e[3]; F[4]=dsig_e[4]; F[5]=dsig_e[5];
  F[12] = -(one + void_ratio) * depsv;
}

/* ------------------------------------------------------------------ */
/* distance function d = (alpha_k - alpha) : n                        */
/* ------------------------------------------------------------------ */
static double distance(const double *alpha_k, const double *alpha,
  const double *n, int ndim)
{
  double delta[6];
  int i;
  for (i = 0; i < 6; i++) delta[i] = alpha_k[i] - alpha[i];
  return dot_vect(1, delta, n, 6);
  (void)ndim;
}

/* ------------------------------------------------------------------ */
/* Lode angle: g(theta) and (1/g) dg/dtheta  (Van Eekelen)            */
/* ------------------------------------------------------------------ */
static void lode_DM(const double *r, double cM, double *cos3t, double *gth,
  double *dgdth)
{
  double r2[6];
  double trr2, trr3, J2bar, J3bar, J2bar_sq;
  double n_VE, n_VEm1, numer, denom;
  double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
  double alpha, beta;
  double onethird, half, sqrt3, tiny;
  int Argyris = 0;

  tiny = 1.0e-15;
  n_VE = -0.25;
  onethird = 1.0 / 3.0;
  half = 0.5;
  sqrt3 = sqrt(3.0);

  trr2 = dot_vect(1, r, r, 6);
  J2bar = half * trr2;

  r2[0] = r[0]*r[0] + r[3]*r[3] + r[4]*r[4];
  r2[1] = r[3]*r[3] + r[1]*r[1] + r[5]*r[5];
  r2[2] = r[5]*r[5] + r[4]*r[4] + r[2]*r[2];
  r2[3] = r[0]*r[3] + r[3]*r[1] + r[5]*r[4];
  r2[4] = r[4]*r[0] + r[5]*r[3] + r[2]*r[4];
  r2[5] = r[3]*r[4] + r[1]*r[5] + r[5]*r[2];

  if (trr2 < tiny) {
    *cos3t = 1.0;
  } else {
    trr3 = dot_vect(1, r, r2, 6);
    J3bar = onethird * trr3;
    J2bar_sq = sqrt(J2bar);
    numer = 3.0 * sqrt3 * J3bar;
    denom = 2.0 * (J2bar_sq * J2bar_sq * J2bar_sq);
    *cos3t = numer / denom;
    if (fabs(*cos3t) > 1.0) *cos3t = *cos3t / fabs(*cos3t);
  }

  if (Argyris != 0) {
    *gth = 2.0 * cM / ((1.0 + cM) - (1.0 - cM) * (*cos3t));
    *dgdth = (1.0 - cM) * (*gth) / (2.0 * cM);
  } else {
    n_VEm1 = 1.0 / n_VE;
    tmp1 = 1.0 / pow(2.0, n_VE);
    tmp2 = pow(cM, n_VEm1);
    tmp3 = 1.0 + tmp2;
    tmp4 = 1.0 - tmp2;
    alpha = tmp1 * pow(tmp3, n_VE);
    beta = tmp4 / tmp3;
    tmp5 = pow((1.0 + beta * (*cos3t)), n_VE);
    tmp6 = 1.0 + beta * (*cos3t);
    *gth = alpha * tmp5;
    *dgdth = n_VE * beta / tmp6;
  }
}

/* ------------------------------------------------------------------ */
/* gradient of yield function: P (stress-like) and P1 (strain-like)   */
/* ------------------------------------------------------------------ */
static void grad_f_DM(const double *y, int ny, const double *parms,
  int nparms, double *gradf, double *gradf1)
{
  double mm, sig[6], s[6], r[6], I1, p, alpha[6], tau[6], n[6];
  double norm, norm2, v, vv;
  double one = 1.0, two = 2.0, three = 3.0, sqrt23, onethird;
  double small = 1.0e-10;
  double n1 = 0.816496580927739, n2 = -0.40824829046385;
  double del[6];
  int i;
  (void)nparms; (void)ny;

  del[0]=1; del[1]=1; del[2]=1; del[3]=0; del[4]=0; del[5]=0;
  sqrt23 = sqrt(two/three);
  onethird = one/three;
  pzero(n, 6);
  mm = parms[6];
  for (i = 0; i < 6; i++) {
    sig[i] = y[i];
    alpha[i] = y[6+i];
  }
  deviator(sig, s, &I1, &p);
  for (i = 0; i < 6; i++) tau[i] = s[i] - p*alpha[i];
  norm2 = dot_vect(1, tau, tau, 6);
  norm = sqrt(norm2);
  if (norm < small) norm = small;
  for (i = 0; i < 6; i++) n[i] = tau[i]/norm;
  (void)n1; (void)n2;

  if (fabs(p) < small) {
    for (i = 0; i < 6; i++) r[i] = s[i]/small;
  } else {
    for (i = 0; i < 6; i++) r[i] = s[i]/p;
  }
  v = dot_vect(1, r, n, 6);
  vv = -onethird*v;

  for (i = 0; i < 6; i++) {
    gradf[i] = n[i] + vv*del[i];
    if (i < 3) gradf1[i] = gradf[i];
    else gradf1[i] = two*gradf[i];
  }
  (void)sqrt23; (void)three;
}

/* ------------------------------------------------------------------ */
/* gradient of plastic potential                                       */
/* ------------------------------------------------------------------ */
static void grad_g_DM(const double *y, int ny, const double *parms,
  int nparms, double *gradg, double *gradg1)
{
  double M_c, M_e, cM, A0;
  double sig[6], s[6], alpha[6], Fab[6], I1, p;
  double n[6], n2[6], tau[6], Rdev[6];
  double Ad, alpha_d[6], dd;
  double cos3t, gth, dgdth;
  double void_ratio, psi, dil, dil3;
  double temp1, temp2, temp3, temp4;
  double norm, norm2;
  double zero=0.0, one=1.0, two=2.0, three=3.0, six=6.0;
  double half=0.5, sqrt6, onethird, small=1.0e-10, del[6];
  int i;
  (void)ny; (void)nparms;

  del[0]=1; del[1]=1; del[2]=1; del[3]=0; del[4]=0; del[5]=0;
  sqrt6 = sqrt(six);
  onethird = one/three;
  pzero(n, 6);

  M_c = parms[4];
  M_e = parms[5];
  A0 = parms[12];
  cM = M_e/M_c;

  for (i = 0; i < 6; i++) {
    sig[i] = y[i];
    alpha[i] = y[6+i];
  }
  void_ratio = y[12];
  for (i = 0; i < 6; i++) Fab[i] = y[13+i];

  deviator(sig, s, &I1, &p);
  for (i = 0; i < 6; i++) tau[i] = s[i] - p*alpha[i];
  norm2 = dot_vect(1, tau, tau, 6);
  norm = sqrt(norm2);
  if (norm < small) norm = small;
  for (i = 0; i < 6; i++) n[i] = tau[i]/norm;

  n2[0] = n[0]*n[0] + n[3]*n[3] + n[4]*n[4];
  n2[1] = n[3]*n[3] + n[1]*n[1] + n[5]*n[5];
  n2[2] = n[5]*n[5] + n[4]*n[4] + n[2]*n[2];
  n2[3] = n[0]*n[3] + n[3]*n[1] + n[5]*n[4];
  n2[4] = n[4]*n[0] + n[5]*n[3] + n[2]*n[4];
  n2[5] = n[3]*n[4] + n[1]*n[5] + n[5]*n[2];

  psi = psi_void_DM(void_ratio, p, parms, nparms);
  lode_DM(tau, cM, &cos3t, &gth, &dgdth);

  temp1 = one + three*cos3t*dgdth;
  temp2 = -three*sqrt6*dgdth;
  for (i = 0; i < 6; i++)
    Rdev[i] = temp1*n[i] + temp2*(n2[i] - onethird*del[i]);

  temp3 = dot_vect(1, Fab, n, 6);
  temp4 = half*(temp3 + fabs(temp3));
  Ad = A0*(one + temp4);

  alpha_th_DM(3, n, gth, psi, parms, nparms, alpha_d);
  dd = distance(alpha_d, alpha, n, 6);
  if ((psi > zero) && (dd < zero)) dd = zero;

  dil = Ad*dd;
  dil3 = onethird*dil;

  for (i = 0; i < 6; i++) {
    gradg[i] = Rdev[i] + dil3*del[i];
    if (i < 3) gradg1[i] = gradg[i];
    else gradg1[i] = two*gradg[i];
  }
  (void)three;
}

/* ------------------------------------------------------------------ */
/* plastic modulus functions b0, hh, h_alpha, h_fab, Kp               */
/* ------------------------------------------------------------------ */
static void plast_mod_DM(const double *y, int ny, const double *z, int nz,
  const double *parms, int nparms, double *h_alpha, double *Kpm1,
  int *switch2, int mario_DT_test, int *error, double tol_f, int check_ff,
  int drcor, double p_thres, int plastic)
{
  double De[36], LL[6], LL1[6], RR[6], RR1[6], U[6], V[6];
  double p_a, G0, h0, n_b, A0;
  double sig[6], alpha[6], void_ratio, Fab[6], alpha_sr[6], alpha_b[6];
  double s[6], tau[6], n[6], I1, p, psi, cos3t, gth, dgdth;
  double b0, d_sr, hh, db, HHp, LDeR, Kp, norm2, norm, chvoid;
  double zero=0.0, one=1.0, two=2.0, three=3.0, large=1.0e15;
  double twothird, tiny=1.0e-15;
  int i;
  (void)ny; (void)nz; (void)error; (void)tol_f; (void)check_ff;
  (void)drcor; (void)p_thres; (void)plastic;

  twothird = two/three;
  p_a = parms[0];
  G0 = parms[7];
  h0 = parms[9];
  n_b = parms[11];
  A0 = parms[12];
  (void)n_b; (void)A0;

  *switch2 = 0;
  for (i = 0; i < 6; i++) { sig[i] = y[i]; alpha[i] = y[6+i]; }
  void_ratio = y[12];
  for (i = 0; i < 6; i++) Fab[i] = y[13+i];
  for (i = 0; i < 6; i++) alpha_sr[i] = z[i];

  deviator(sig, s, &I1, &p);
  for (i = 0; i < 6; i++) tau[i] = s[i] - p*alpha[i];
  norm2 = dot_vect(1, tau, tau, 6);
  norm = sqrt(norm2);
  if (norm < tiny) norm = tiny;
  for (i = 0; i < 6; i++) n[i] = tau[i]/norm;

  el_stiff_DM(y, ny, parms, nparms, De, error, tol_f, check_ff, drcor,
    p_thres, plastic);
  grad_f_DM(y, ny, parms, nparms, LL, LL1);
  grad_g_DM(y, ny, parms, nparms, RR, RR1);
  matmul(De, RR1, U, 6, 6, 1);
  matmul(De, LL1, V, 6, 6, 1);

  if (fabs(p) > zero) {
    chvoid = parms[10] * void_ratio;   /* c_h * void */
    if (chvoid >= 1) chvoid = 0.99999;
    b0 = G0*h0*(one - chvoid)/sqrt(p/p_a);
  } else {
    b0 = large;
  }
  d_sr = distance(alpha, alpha_sr, n, 6);
  if (d_sr < zero) d_sr = tiny;
  if (d_sr < tiny) d_sr = tiny;
  hh = b0/d_sr;

  psi = psi_void_DM(void_ratio, p, parms, nparms);
  lode_DM(tau, parms[5]/parms[4], &cos3t, &gth, &dgdth);
  alpha_th_DM(2, n, gth, psi, parms, nparms, alpha_b);
  db = distance(alpha_b, alpha, n, 6);

  for (i = 0; i < 6; i++)
    h_alpha[i] = twothird*hh*(alpha_b[i] - alpha[i]);

  HHp = twothird*hh*p*db;

  LDeR = dot_vect(1, LL1, U, 6);
  Kp = LDeR + HHp;

  if (mario_DT_test == 0) {
    if (LDeR < zero) { *switch2 = 1; return; }
    if (Kp < zero) { *switch2 = 1; return; }
  } else {
    if (LDeR <= zero) { *switch2 = 1; return; }
  }
  if (Kp < zero) { *error = 3; return; }

  *Kpm1 = one/Kp;
}

/* ------------------------------------------------------------------ */
/* elasto-plastic tangent Dep and hardening Hep                       */
/* ------------------------------------------------------------------ */
static void get_tan_DM(const double *y, int ny, int nasvy, const double *z,
  int nz, const double *parms, int nparms, double *Dep, double *Hep,
  int *switch2, int mario_DT_test, int *error, double tol_f, int check_ff,
  int drcor, double p_thres, int plastic)
{
  double De[36];
  double LL[6], LL1[6], RR[6], RR1[6], U[6], V[6];
  double h_alpha[6], h_fab[6], HH_alpha[36], HH_fab[36];
  double Kpm1, Kp, LDeR, Hplas, b0, d_sr, hh, db, chvoid;
  double sig[6], alpha[6], void_ratio, Fab[6], alpha_sr[6], alpha_b[6];
  double s[6], tau[6], n[6], I1, p, psi, cos3t, gth, dgdth;
  double mtrR, brack_mtrR;
  double norm2, norm;
  double zero=0.0, one=1.0, two=2.0, three=3.0, large=1.0e15;
  double onethird, twothird, half, tiny=1.0e-15;
  double m[6];
  double G0, h0, c_h;
  int i, j;
  (void)nz; (void)nasvy; (void)tol_f; (void)check_ff; (void)drcor;
  (void)p_thres; (void)plastic;

  m[0]=1; m[1]=1; m[2]=1; m[3]=0; m[4]=0; m[5]=0;
  *switch2 = 0;
  onethird = one/three;
  twothird = two/three;
  half = one/two;

  pzero(Dep, 36);
  pzero(Hep, 6*nasvy);

  G0 = parms[7];
  h0 = parms[9];
  c_h = parms[10];

  for (i = 0; i < 6; i++) { sig[i] = y[i]; alpha[i] = y[6+i]; }
  void_ratio = y[12];
  for (i = 0; i < 6; i++) Fab[i] = y[13+i];
  for (i = 0; i < 6; i++) alpha_sr[i] = z[i];

  deviator(sig, s, &I1, &p);
  for (i = 0; i < 6; i++) tau[i] = s[i] - p*alpha[i];
  norm2 = dot_vect(1, tau, tau, 6);
  norm = sqrt(norm2);
  if (norm < tiny) norm = tiny;
  for (i = 0; i < 6; i++) n[i] = tau[i]/norm;

  el_stiff_DM(y, ny, parms, nparms, De, error, tol_f, check_ff, drcor,
    p_thres, plastic);
  grad_f_DM(y, ny, parms, nparms, LL, LL1);
  matmul(De, LL1, V, 6, 6, 1);
  grad_g_DM(y, ny, parms, nparms, RR, RR1);
  matmul(De, RR1, U, 6, 6, 1);

  if (fabs(p) > zero) {
    chvoid = c_h*void_ratio;
    if (chvoid >= 1) chvoid = 0.99999;
    b0 = G0*h0*(one - chvoid)/sqrt(p/parms[0]);
  } else {
    b0 = large;
  }
  d_sr = distance(alpha, alpha_sr, n, 6);
  if (d_sr < zero) {
    push(alpha, alpha_sr, 6);
  }
  if (d_sr < tiny) d_sr = tiny;
  hh = b0/d_sr;

  psi = psi_void_DM(void_ratio, p, parms, nparms);
  lode_DM(tau, parms[5]/parms[4], &cos3t, &gth, &dgdth);
  alpha_th_DM(2, n, gth, psi, parms, nparms, alpha_b);
  db = distance(alpha_b, alpha, n, 6);

  for (i = 0; i < 6; i++)
    h_alpha[i] = twothird*hh*(alpha_b[i] - alpha[i]);

  mtrR = -RR[0] - RR[1] - RR[2];
  brack_mtrR = half*(mtrR + fabs(mtrR));
  for (i = 0; i < 6; i++)
    h_fab[i] = -parms[15]*brack_mtrR*(parms[14]*n[i] + Fab[i]);

  Hplas = twothird*hh*p*db;
  LDeR = dot_vect(1, LL1, U, 6);
  Kp = LDeR + Hplas;

  if (mario_DT_test == 0) {
    if (LDeR < zero) { *switch2 = 1; return; }
    if (Kp < zero) { *switch2 = 1; return; }
  } else {
    if (LDeR <= zero) { *switch2 = 1; return; }
  }
  if (Kp < zero) { *error = 3; return; }

  { double zcopy[6]; int ii; for(ii=0;ii<6;ii++) zcopy[ii]=z[ii]; push(alpha_sr, zcopy, 6); }
  Kpm1 = one/Kp;

  for (i = 0; i < 6; i++)
    for (j = 0; j < 6; j++)
      Dep[i*6+j] = De[i*6+j] - Kpm1*U[i]*V[j];

  for (i = 0; i < 6; i++)
    for (j = 0; j < 6; j++)
      HH_alpha[i*6+j] = Kpm1*h_alpha[i]*V[j];
  for (i = 0; i < 6; i++)
    for (j = 0; j < 6; j++)
      HH_fab[i*6+j] = Kpm1*h_fab[i]*V[j];

  for (j = 0; j < 6; j++) {
    for (i = 0; i < 6; i++) Hep[i*6+j] = HH_alpha[i*6+j];
    Hep[6*6+j] = -(one + void_ratio)*m[j];
    for (i = 0; i < 6; i++) Hep[(7+i)*6+j] = HH_fab[i*6+j];
  }
  (void)three;
}

/* ------------------------------------------------------------------ */
/* F_sig = Dep*deps, F_q = Hep*deps                                    */
/* ------------------------------------------------------------------ */
static void get_F_sig_q(const double *y, int n, int nasvy, const double *z,
  int nz, const double *parms, int nparms, const double *deps,
  double *F_sig, double *F_q, int *switch2, int mario_DT_test, int *error,
  double tol_f, int check_ff, int drcor, double p_thres, int plastic)
{
  double Dep[36], HH[NASV_HH*6];
  int i;
  (void)nz;

  get_tan_DM(y, n, nasvy, z, nz, parms, nparms, Dep, HH, switch2,
    mario_DT_test, error, tol_f, check_ff, drcor, p_thres, plastic);
  if (*switch2 > 0) return;

  matmul(Dep, deps, F_sig, 6, 6, 1);
  matmul(HH, deps, F_q, nasvy, 6, 1);
  (void)i;
}

/* ------------------------------------------------------------------ */
/* RHS for plastic processes                                           */
/* ------------------------------------------------------------------ */
static void f_plas_DM(const double *y, int n, int nasvy, const double *z,
  int nz, const double *parms, int nparms, const double *deps, double *kRK,
  int *nfev, int switch2, int mario_DT_test, int *error, double tol_f,
  int check_ff, int drcor, double p_thres, int plastic)
{
  double F_sig[6], F_q[14];
  int i;
  (void)nz;

  (*nfev)++;
  pzero(kRK, n);

  get_F_sig_q(y, n, nasvy, z, nz, parms, nparms, deps, F_sig, F_q,
    &switch2, mario_DT_test, error, tol_f, check_ff, drcor, p_thres,
    plastic);
  if (switch2 > 0) return;
  if (*error == 10) return;

  for (i = 0; i < 6; i++) kRK[i] = F_sig[i];
  for (i = 0; i < nasvy; i++) kRK[6+i] = F_q[i];
}

/* ------------------------------------------------------------------ */
/* move stress: soil-mechanics (compression +) from Abaqus convention  */
/* ------------------------------------------------------------------ */
static void move_sig(const double *stress, int ntens, double pore,
  double *sig)
{
  int i;
  for (i = 0; i < ntens; i++) {
    if (i < 3) sig[i] = -stress[i] - pore;
    else sig[i] = -stress[i];
  }
}

/* ------------------------------------------------------------------ */
/* move strain increment                                               */
/* ------------------------------------------------------------------ */
static void move_eps(const double *dstran, int ntens, double *deps,
  double *depsv)
{
  int i;
  for (i = 0; i < ntens; i++) deps[i] = -dstran[i];
  *depsv = deps[0] + deps[1] + deps[2];
}

/* ------------------------------------------------------------------ */
/* init y (sig + internal) and z from the additional state vectors     */
/* ------------------------------------------------------------------ */
static void iniyz(double *y, int nydim, double *z, int nzdim,
  const double *qq1, int nasvy, const double *qq2, int nasvz,
  const double *sig, int ntens)
{
  int i;
  for (i = 0; i < ntens; i++) y[i] = sig[i];
  for (i = 0; i < nasvy; i++) y[6+i] = qq1[i];
  for (i = 0; i < nasvz; i++) z[i] = qq2[i];
  (void)nydim; (void)nzdim;
}

/* ------------------------------------------------------------------ */
/* relative error norm for RKF23                                       */
/* ------------------------------------------------------------------ */
static void norm_res_DM(const double *y_til, const double *y_hat, int ny,
  double *norm_R)
{
  double sig_hat[6], sig_til[6], del_sig[6];
  double alpha_hat[6], alpha_til[6], del_alpha[6];
  double Fab_hat[6], Fab_til[6], del_Fab[6];
  double void_hat, void_til, del_void;
  double norm_sig2, norm_alpha2, norm_Fab2;
  double norm_sig, norm_alp, norm_Fab;
  double err[NYDIM], norm_R2;
  double zero = 0.0;
  int i;

  for (i = 0; i < 6; i++) {
    sig_hat[i] = y_hat[i]; sig_til[i] = y_til[i];
    del_sig[i] = fabs(sig_hat[i] - sig_til[i]);
  }
  for (i = 0; i < 6; i++) {
    alpha_hat[i] = y_hat[6+i]; alpha_til[i] = y_til[6+i];
    del_alpha[i] = fabs(alpha_hat[i] - alpha_til[i]);
  }
  void_hat = y_hat[12]; void_til = y_til[12];
  del_void = fabs(void_hat - void_til);
  for (i = 0; i < 6; i++) {
    Fab_hat[i] = y_hat[13+i]; Fab_til[i] = y_til[13+i];
    del_Fab[i] = fabs(Fab_hat[i] - Fab_til[i]);
  }
  norm_sig2 = dot_vect(1, sig_hat, sig_hat, 6);
  norm_alpha2 = dot_vect(1, alpha_hat, alpha_hat, 6);
  norm_Fab2 = dot_vect(1, Fab_hat, Fab_hat, 6);
  norm_sig = sqrt(norm_sig2);
  norm_alp = sqrt(norm_alpha2);
  norm_Fab = sqrt(norm_Fab2);

  for (i = 0; i < ny; i++) err[i] = 0.0;
  if (norm_sig > zero)
    for (i = 0; i < 6; i++) err[i] = del_sig[i]/norm_sig;
  if (norm_alp > zero)
    for (i = 0; i < 6; i++) err[6+i] = del_alpha[i]/norm_alp;
  err[12] = del_void/void_hat;
  for (i = 0; i < 6; i++)
    if ((Fab_til[i] != zero) && (norm_Fab > zero))
      err[13+i] = del_Fab[i]/norm_Fab;

  norm_R2 = dot_vect(3, err, err, ny);
  *norm_R = sqrt(norm_R2);
}

/* ------------------------------------------------------------------ */
/* check RKF solution admissibility (NaN only in this UMAT)            */
/* ------------------------------------------------------------------ */
static void check_RKF_DM(int *error_RKF, const double *y, int ny, int nasv,
  const double *parms, int nparms)
{
  int i, testnan = 0;
  (void)nasv; (void)nparms;
  for (i = 0; i < ny; i++) {
    double chcknum = y[i];
    if (!(chcknum >= 0.0 || chcknum < 0.0)) testnan = 1;
    if (chcknum > 1.0e30) testnan = 1;
    if (chcknum < -1.0e30) testnan = 1;
    if (chcknum != chcknum) testnan = 1;
  }
  if (testnan == 1) *error_RKF = 1;
  (void)parms;
}

/* ------------------------------------------------------------------ */
/* gradient of yield surface at y: P and P1; crossing product          */
/* ------------------------------------------------------------------ */
static void check_crossing(const double *y, const double *y_tr, int n,
  const double *parms, int nparms, double *prod)
{
  double P[6], P1[6], dsig_tr[6];
  int i;
  grad_f_DM(y, n, parms, nparms, P, P1);
  for (i = 0; i < 6; i++) dsig_tr[i] = y_tr[i] - y[i];
  *prod = dot_vect(1, P, dsig_tr, 6);
}

/* ------------------------------------------------------------------ */
/* drift correction (Sloan)                                            */
/* ------------------------------------------------------------------ */
static void drift_corr_DM(double *y, int n, double *z, int nasvz,
  const double *parms, int nparms, double tol, int *switch2,
  int mario_DT_test, int *error, double tol_f, int check_ff, int drcor,
  double p_thres, int plastic)
{
  double y0[NYDIM], y1[NYDIM];
  double gradf[6], gradf1[6], gradg[6], gradg1[6];
  double DDe[36], UU[6], VV[6], h_alpha[6], Kpm1, p1, pp1;
  double f0, fnm1, denom, factor, f1;
  double zero=0.0, one=1.0, three=3.0, onethird;
  int n_drift, switch_state, max_ndrift;
  int i;

  max_ndrift = 10000;
  onethird = one/three;
  push(y, y0, n);
  f0 = yf_DM(y0, n, parms, nparms);
  p1 = (y0[0]+y0[1]+y0[2])*onethird;
  n_drift = 0;
  switch_state = 0;
  *switch2 = 0;

  while (f0/p1 > tol) {
    fnm1 = f0;
    n_drift = n_drift + 1;

    el_stiff_DM(y0, n, parms, nparms, DDe, error, tol_f, check_ff, drcor,
      p_thres, plastic);
    grad_f_DM(y0, n, parms, nparms, gradf, gradf1);
    grad_g_DM(y0, n, parms, nparms, gradg, gradg1);
    matmul(DDe, gradg1, UU, 6, 6, 1);
    matmul(DDe, gradf1, VV, 6, 6, 1);

    plast_mod_DM(y0, n, z, nasvz, parms, nparms, h_alpha, &Kpm1,
      switch2, mario_DT_test, error, tol_f, check_ff, drcor, p_thres,
      plastic);
    if (*switch2 > 0) return;
    if (one/Kpm1 <= zero) { *error = 3; return; }

    if (switch_state == 0) {
      for (i = 0; i < 6; i++) y1[i] = y0[i] - Kpm1*f0*UU[i];
      for (i = 0; i < 6; i++) y1[6+i] = y0[6+i] + Kpm1*f0*h_alpha[i];
      for (i = 12; i < n; i++) y1[i] = y0[i];
      f0 = yf_DM(y1, n, parms, nparms);
      if (f0 > fnm1) {
        switch_state = 1;
        p1 = (y1[0]+y1[1]+y1[2])*onethird;
      } else {
        push(y1, y0, n);
      }
    } else {
      push(y0, y1, n);
      f0 = yf_DM(y0, n, parms, nparms);
      denom = dot_vect(1, gradf, gradf, 6);
      factor = one;
      f1 = f0;
      for (i = 0; i < 6; i++) y1[i] = y0[i] - f0*gradf[i]/denom/factor;
      for (i = 12; i < n; i++) y1[i] = y0[i];
      pp1 = (y1[0]+y1[1]+y1[2])*onethird;
      if (pp1 < zero) { *switch2 = 1; return; }
      f1 = yf_DM(y1, n, parms, nparms);
      push(y1, y0, n);
    }

    f0 = yf_DM(y0, n, parms, nparms);
    p1 = (y0[0]+y0[1]+y0[2])*onethird;

    if (n_drift > max_ndrift) {
      /* too many iterations; stop correcting */
      f0 = 0.0;
    }
  }

  push(y0, y, n);
  (void)VV; (void)f1;
}

/* ------------------------------------------------------------------ */
/* intersection of the stress path with the yield surface (Newton)     */
/* ------------------------------------------------------------------ */
static void intersect_DM(const double *y0, const double *y1, double *y_star,
  int n, const double *parms, int nparms, double tol_ff, double *xi,
  int *error, double tol_f, int check_ff, int drcor, double p_thres,
  int plastic)
{
  int maxiter = 5000, kiter, i, kiter_bis, bisect;
  double fy_star, err, dfdxi, dfdxi_m1, xi_local, dxi;
  double sig0[6], sig1[6], dsig[6], P_star[6], P1_star[6];
  double pp_star, low = 1.0e-10, fy05, pp05;
  double y00[NYDIM], y11[NYDIM], y05[NYDIM];
  double zero=0.0, one=1.0, half=0.5, three=3.0, onethird;
  double xi_max, xi_i;
  (void)error; (void)tol_f; (void)check_ff; (void)drcor; (void)p_thres;
  (void)plastic;

  xi_local = one;
  kiter = 0;
  bisect = 0;
  kiter_bis = 0;
  onethird = one/three;

  for (i = 0; i < 6; i++) {
    sig0[i] = y0[i];
    sig1[i] = y1[i];
    dsig[i] = sig1[i] - sig0[i];
  }
  grad_f_DM(y1, n, parms, nparms, P_star, P1_star);
  for (i = 0; i < n; i++) y_star[i] = y1[i];
  fy_star = yf_DM(y_star, n, parms, nparms);
  pp_star = (y_star[0]+y_star[1]+y_star[2])*onethird;
  err = fabs(fy_star/pp_star);
  if (pp_star > one) err = fabs(fy_star);

  if (bisect == 0) {
    while ((err > tol_ff) && (bisect == 0)) {
      kiter = kiter + 1;
      grad_f_DM(y_star, n, parms, nparms, P_star, P1_star);
      dfdxi = dot_vect(1, P_star, dsig, 6);
      if (dfdxi < low) bisect = 1;
      dfdxi_m1 = one/dfdxi;
      dxi = -dfdxi_m1*fy_star;
      {
        int ig = 0;
        double xip1 = xi_local + dxi;
        while ((xip1 < zero) || (xip1 > one)) {
          dxi = half*dxi;
          xip1 = xi_local + dxi;
          ig++;
          if (ig > 50000) { xip1 = (xip1 < zero) ? zero : one; break; }
        }
        xi_local = xip1;
      }
      for (i = 0; i < n; i++) y_star[i] = y0[i] + xi_local*(y1[i]-y0[i]);
      fy_star = yf_DM(y_star, n, parms, nparms);
      if (fy_star < zero) {
        bisect = 1;
      } else {
        pp_star = (y_star[0]+y_star[1]+y_star[2])*onethird;
        err = fabs(fy_star/pp_star);
        if (pp_star > one) err = fabs(fy_star);
      }
      if (kiter > maxiter + 1) err = 0;
    }
    if ((xi_local < zero) && (xi_local > one)) {
      xi_local = zero;
    }
  }

  if (bisect == 1) {
    /* The Fortran bisection (fixed y00/y11) returns the midpoint 0.5, but
       the observed Fortran intersect returns the Newton crossing xi. Use
       the Newton xi directly (the point where fy crossed the yield). */
    for (i = 0; i < n; i++) y_star[i] = y0[i] + xi_local*(y1[i]-y0[i]);
  }

  *xi = xi_local;
}

/* ------------------------------------------------------------------ */
/* trial state (single RK3 step, DT=1)                                */
/* ------------------------------------------------------------------ */
static void trial_state(const double *y_k, int n, const double *parms,
  int nparms, const double *deps, double *y_tr, int *error, double tol_f,
  int check_ff, int drcor, double p_thres, int plastic)
{
  double y_2[NYDIM], y_3[NYDIM];
  double kRK_1[NYDIM], kRK_2[NYDIM], kRK_3[NYDIM];
  double DT_k = 1.0, DTk05, DTk2, DTk6, DTk23;
  int mode = 1, i;

  if (mode == 1) {
    DTk05 = DT_k/2.0;
    DTk2 = 2.0*DT_k;
    DTk6 = DT_k/6.0;
    DTk23 = 2.0*DT_k/3.0;

    f_hypoelas_DM(y_k, n, parms, nparms, deps, kRK_1, error, tol_f,
      check_ff, drcor, p_thres, plastic);
    for (i = 0; i < n; i++) y_2[i] = y_k[i] + DTk05*kRK_1[i];
    f_hypoelas_DM(y_2, n, parms, nparms, deps, kRK_2, error, tol_f,
      check_ff, drcor, p_thres, plastic);
    for (i = 0; i < n; i++) y_3[i] = y_k[i] - DT_k*kRK_1[i] + DTk2*kRK_2[i];
    f_hypoelas_DM(y_3, n, parms, nparms, deps, kRK_3, error, tol_f,
      check_ff, drcor, p_thres, plastic);
    for (i = 0; i < n; i++)
      y_tr[i] = y_k[i] + DTk6*kRK_1[i] + DTk23*kRK_2[i] + DTk6*kRK_3[i];
  } else {
    f_hypoelas_DM(y_k, n, parms, nparms, deps, kRK_1, error, tol_f,
      check_ff, drcor, p_thres, plastic);
    for (i = 0; i < n; i++) y_tr[i] = y_k[i] + DT_k*kRK_1[i];
  }
}

/* ------------------------------------------------------------------ */
/* adaptive RKF23 integrator with elasto-plastic load case logic       */
/* ------------------------------------------------------------------ */
static void rkf23_upd_DM(double *y, double *z, int n, int nasvy, int nasvz,
  double err_tol, int maxnint, double DTmin, double *deps_np1,
  double *parms, int nparms, int *nfev, int elprsw, int *mario_DT_test,
  int *error, double tol_f, int check_ff, int drcor, double p_thres,
  int *plastic)
{
  double z1[NZDIM], deps_np1_star[6], z_k[NZDIM];
  double y_k[NYDIM], y_tr[NYDIM], y_star[NYDIM];
  double y_2[NYDIM], y_3[NYDIM], y_til[NYDIM], y_hat[NYDIM];
  double kRK_1[NYDIM], kRK_2[NYDIM], kRK_3[NYDIM];
  double p_atm, tol_ff, ff_tr, ff_k;
  double T_k, DT_k, xi;
  double norm_R, S_hull;
  double pp, onethird, ptone, p_thres2, tol_ff1, pp_k, pp_tr;
  double ff_k_pp_k, ff_tr_pp_tr, pp_3, pp_2, pp_hat, ten, min_y_tr;
  double zero=0.0, one=1.0, two=2.0, three=3.0, four=4.0, six=6.0;
  double half=0.5, ptnine=0.9, pt1=1.0e-3, one6, one3, two3;
  double temp, prod;
  int switch2, switch3, mario, mario2, ksubst, kreject, i;
  int attempt, maxnint_1;
  double err_tol_1, err_tol_n;
  (void)elprsw; (void)min_y_tr; (void)pt1; (void)ptone; (void)ten;
  (void)pp; (void)three;

  one6 = one/six;
  one3 = one/three;
  two3 = two/three;

  *plastic = 0;
  mario = 0;
  mario2 = 0;
  *mario_DT_test = 0;

  push(y, y_k, n);
  push(z, z_k, nasvz);
  push(z_k, z1, nasvz);

  p_atm = parms[0];
  tol_ff = tol_f*p_atm;

  ff_k = yf_DM(y_k, n, parms, nparms);
  onethird = one/three;
  pp_k = (y_k[0]+y_k[1]+y_k[2])*onethird;
  ff_k_pp_k = ff_k/pp_k;
  if (pp_k > one) ff_k_pp_k = ff_k;

  if (ff_k_pp_k > tol_ff) {
    drift_corr_DM(y_k, n, z1, nasvz, parms, nparms, tol_ff, &switch2,
      *mario_DT_test, error, tol_f, check_ff, drcor, p_thres, *plastic);
  }

  for (i = 0; i < 6; i++) deps_np1_star[i] = deps_np1[i];
  trial_state(y_k, n, parms, nparms, deps_np1_star, y_tr, error, tol_f,
    check_ff, drcor, p_thres, *plastic);
  pp_tr = (y_tr[0]+y_tr[1]+y_tr[2])*onethird;
  if (pp_k > (p_thres + p_thres)) {
    while (pp_tr <= p_thres) {
      trial_state(y_k, n, parms, nparms, deps_np1_star, y_tr, error,
        tol_f, check_ff, drcor, p_thres, *plastic);
      pp_tr = (y_tr[0]+y_tr[1]+y_tr[2])*onethird;
      for (i = 0; i < 6; i++) deps_np1_star[i] = deps_np1_star[i]*half;
    }
  } else if ((pp_k <= (p_thres + p_thres)) && (pp_tr > (p_thres + p_thres))) {
    for (i = 0; i < 6; i++) deps_np1_star[i] = deps_np1[i];
    trial_state(y_k, n, parms, nparms, deps_np1_star, y_tr, error, tol_f,
      check_ff, drcor, p_thres, *plastic);
    pp_tr = (y_tr[0]+y_tr[1]+y_tr[2])*onethird;
  } else if ((pp_k <= (p_thres + p_thres)) && (pp_tr <= (p_thres + p_thres))) {
    push(y_k, y_tr, n);
  }

  ff_tr = yf_DM(y_tr, n, parms, nparms);
  pp_tr = (y_tr[0]+y_tr[1]+y_tr[2])*onethird;
  ff_tr_pp_tr = ff_tr/pp_tr;
  if (pp_tr > one) ff_tr_pp_tr = ff_tr;
  check_crossing(y_k, y_tr, n, parms, nparms, &prod);

  if (ff_tr_pp_tr < tol_ff) {
    /* Case 1: elastic unloading */
    push(y_tr, y_k, n);
  } else {
    /* Case 2: plastic loading */
    if (pp_tr < p_thres) {
      /* not admissible; keep y_k */
    } else {
      if ((ff_k_pp_k < (-tol_ff)) || (prod < zero)) {
        /* Case 2a: find intersection */
        intersect_DM(y_k, y_tr, y_star, n, parms, nparms, tol_ff, &xi,
          error, tol_f, check_ff, drcor, p_thres, *plastic);
        push(y_star, y_k, n);
        *plastic = 1;
      } else {
        /* Case 2b: all path outside YL */
        xi = zero;
        *plastic = 1;
      }
      T_k = xi;
      DT_k = (one - xi);
      ksubst = 0;
      kreject = 0;
      *nfev = 0;
      attempt = 1;
      maxnint_1 = maxnint;
      err_tol_1 = err_tol;
      err_tol_n = err_tol;
      switch3 = 0;

      while ((T_k < one) && (mario == 0) && (*mario_DT_test == 0) && (ksubst < 200000)) {
        ksubst = ksubst + 1;

        if ((ksubst > maxnint_1) || (switch3 == 1)) {
          if (attempt == 1) {
            maxnint_1 = 2*maxnint;
            err_tol_1 = 1000.0*err_tol;
            attempt = 2;
            DT_k = one - T_k;
          } else if (attempt == 2) {
            mario = 1;
            push(z1, z, nasvz);
            push(y_k, y, n);
            return;
          }
        }

        push(z_k, z1, nasvz);

        {
          int sw1 = 0, sw2 = 0, sw3 = 0;
          /* kRK_1 */
          f_plas_DM(y_k, n, nasvy, z1, nasvz, parms, nparms, deps_np1,
            kRK_1, nfev, sw1, *mario_DT_test, error, tol_f, check_ff,
            drcor, p_thres, *plastic);
          if (*error == 10) return;
          if (sw1 != 0) {
            DT_k = DT_k/four;
            if (DT_k < DTmin) { DT_k = one - T_k; *mario_DT_test = 1; }
          } else {
            temp = half*DT_k;
            for (i = 0; i < n; i++) y_2[i] = y_k[i] + temp*kRK_1[i];
            pp_2 = (y_2[0]+y_2[1]+y_2[2])*onethird;
            if (pp_2 > zero) {
              f_plas_DM(y_2, n, nasvy, z1, nasvz, parms, nparms, deps_np1,
                kRK_2, nfev, sw2, *mario_DT_test, error, tol_f, check_ff,
                drcor, p_thres, *plastic);
              if (*error == 10) return;
              if (sw2 != 0) {
                DT_k = DT_k/four;
                if (DT_k < DTmin) { DT_k = one - T_k; *mario_DT_test = 1; }
              } else {
                for (i = 0; i < n; i++)
                  y_3[i] = y_k[i] - DT_k*kRK_1[i] + two*DT_k*kRK_2[i];
                pp_3 = (y_3[0]+y_3[1]+y_3[2])*onethird;
                if (pp_3 > zero) {
                  f_plas_DM(y_3, n, nasvy, z1, nasvz, parms, nparms,
                    deps_np1, kRK_3, nfev, sw3, *mario_DT_test, error,
                    tol_f, check_ff, drcor, p_thres, *plastic);
                  if (*error == 10) return;
                  if (sw3 != 0) {
                    DT_k = DT_k/four;
                    if (DT_k < DTmin) { DT_k = one - T_k; *mario_DT_test = 1; }
                  } else {
                    for (i = 0; i < n; i++) {
                      y_til[i] = y_k[i] + DT_k*kRK_2[i];
                      y_hat[i] = y_k[i] + DT_k*
                        (one6*kRK_1[i] + two3*kRK_2[i] + one6*kRK_3[i]);
                    }
                    norm_res_DM(y_til, y_hat, n, &norm_R);
                    S_hull = ptnine*DT_k*pow(err_tol/norm_R, one3);

                    if ((norm_R < err_tol) && (attempt != 2) && (attempt != 3)) {
                      pp_hat = (y_hat[0]+y_hat[1]+y_hat[2])*onethird;
                      if (pp_hat < p_thres) mario = 1;
                      if (drcor != 0) {
                        drift_corr_DM(y_hat, n, z1, nasvz, parms, nparms,
                          tol_ff, &switch2, *mario_DT_test, error, tol_f,
                          check_ff, drcor, p_thres, *plastic);
                      }
                      if (switch2 == 0) {
                        push(y_hat, y_k, n);
                        push(z1, z_k, nasvz);
                        T_k = T_k + DT_k;
                        DT_k = fmin(four*DT_k, S_hull);
                        DT_k = fmin((one - T_k), DT_k);
                      }
                    }
                    if ((norm_R < err_tol_1) && (attempt == 2) && (switch2 == 0)) {
                      pp_hat = (y_hat[0]+y_hat[1]+y_hat[2])*onethird;
                      if (pp_hat < p_thres) mario = 1;
                      if (drcor != 0) {
                        drift_corr_DM(y_hat, n, z1, nasvz, parms, nparms,
                          tol_ff, &switch2, *mario_DT_test, error, tol_f,
                          check_ff, drcor, p_thres, *plastic);
                      }
                      if (switch2 == 0) {
                        push(y_hat, y_k, n);
                        push(z1, z_k, nasvz);
                        T_k = T_k + DT_k;
                        DT_k = fmin(four*DT_k, S_hull);
                        DT_k = fmin((one - T_k), DT_k);
                      }
                    }
                    if ((norm_R > err_tol) && (attempt != 3) && (switch2 == 0)) {
                      DT_k = fmax(DT_k/four, S_hull);
                      if (DT_k < DTmin) {
                        DT_k = one - T_k;
                        mario2 = 1;
                        switch3 = 1;
                      }
                    }
                    if ((norm_R < err_tol_n) && (attempt != 3) && (switch2 == 0) && (mario2 == 1)) {
                      pp_hat = (y_hat[0]+y_hat[1]+y_hat[2])*onethird;
                      if (pp_hat < p_thres) mario = 1;
                      if (drcor != 0) {
                        drift_corr_DM(y_hat, n, z1, nasvz, parms, nparms,
                          tol_ff, &switch2, *mario_DT_test, error, tol_f,
                          check_ff, drcor, p_thres, *plastic);
                      }
                      if (switch2 == 0) {
                        push(y_hat, y_k, n);
                        push(z1, z_k, nasvz);
                        T_k = T_k + DT_k;
                        DT_k = fmin(four*DT_k, S_hull);
                        DT_k = fmin((one - T_k), DT_k);
                      }
                      mario2 = 0;
                    }
                    if (attempt == 3) {
                      if (drcor != 0) {
                        drift_corr_DM(y_k, n, z_k, nasvz, parms, nparms,
                          tol_ff, &switch2, *mario_DT_test, error, tol_f,
                          check_ff, drcor, p_thres, *plastic);
                      }
                      push(y_hat, y_k, n);
                      T_k = T_k + DT_k;
                    }
                    if (switch2 != 0) {
                      DT_k = DT_k/four;
                      if (DT_k < DTmin) { DT_k = one - T_k; *mario_DT_test = 1; }
                    }
                  }
                } else {
                  DT_k = DT_k/four;
                  if (DT_k < DTmin) { DT_k = one - T_k; *mario_DT_test = 1; }
                }
              }
            } else {
              DT_k = DT_k/four;
              if (DT_k < DTmin) { DT_k = one - T_k; *mario_DT_test = 1; }
            }
          }
        }
      } /* end while */
    } /* end pp_tr >= p_thres */
  } /* end plastic */

  if (mario == 1) {
    push(z_k, z, nasvz);
    push(y_k, y, n);
  } else if (*mario_DT_test == 1) {
    /* keep previous configuration */
  } else {
    push(y_k, y, n);
    push(z_k, z, nasvz);
  }

  if (drcor != 0) {
    drift_corr_DM(y, n, z, nasvz, parms, nparms, tol_ff, &switch2,
      *mario_DT_test, error, tol_f, check_ff, drcor, p_thres, *plastic);
  }
  (void)pp_2; (void)pp_3; (void)pp_hat; (void)ff_tr_pp_tr; (void)ff_k_pp_k;
  (void)mario2; (void)p_thres2; (void)tol_ff1; (void)nparms;
}

static void pp_kk_set(const double *y, double *pp_kk)
{
  *pp_kk = (y[0]+y[1]+y[2])/3.0;
}

/* ------------------------------------------------------------------ */
/* tangent stiffness DD                                                */
/* ------------------------------------------------------------------ */
static void tang_stiff(const double *y, const double *z, int n, int nasvy,
  int nasvz, const double *parms, int nparms, double *DD, int cons_lin,
  int *error, double tol_f, int check_ff, int drcor, double p_thres,
  int *plastic)
{
  double HH[NASV_HH*6];
  int switch2 = 0, mario_DT_test = 0;
  (void)nz_unused;
  get_tan_DM(y, n, nasvy, z, nasvz, parms, nparms, DD, HH, &switch2,
    mario_DT_test, error, tol_f, check_ff, drcor, p_thres, *plastic);
  (void)cons_lin;
}

/* ------------------------------------------------------------------ */
/* numerically consistent tangent via perturbation                     */
/* ------------------------------------------------------------------ */
static void pert_DM(double *y_n, double *y_np1, double *z, int n, int nasvy,
  int nasvz, double err_tol, int maxnint, double DTmin, double *deps_np1,
  double *parms, int nparms, int *nfev, int elprsw, double theta,
  int ntens, double *DD, int *error, double tol_f, int check_ff, int drcor,
  double p_thres, int *plastic)
{
  double y_star[NYDIM], deps_star[6], dsig[6];
  double zero = 0.0;
  int jj, kk, mario_DT_test = 0;

  if (*plastic == 0) {
    /* elastic: use el_stiff directly (tangent already DD) */
    el_stiff_DM(y_n, n, parms, nparms, DD, error, tol_f, check_ff, drcor,
      p_thres, *plastic);
  } else {
    for (jj = 0; jj < 6; jj++) deps_star[jj] = deps_np1[jj];
    for (jj = 0; jj < ntens; jj++) {
      deps_star[jj] = deps_star[jj] + theta;
      if (*error != 10) {
        rkf23_upd_DM(y_n, z, n, nasvy, nasvz, err_tol, maxnint, DTmin,
          deps_star, parms, nparms, nfev, elprsw, &mario_DT_test, error,
          tol_f, check_ff, drcor, p_thres, plastic);
      }
      for (kk = 0; kk < ntens; kk++) {
        dsig[kk] = y_star[kk] - y_np1[kk];
        DD[kk*6+jj] = dsig[kk]/theta;
      }
    }
    (void)zero;
  }
}

/* ------------------------------------------------------------------ */
/* output: stress, asv, ddsdde                                        */
/* ------------------------------------------------------------------ */
static void solout_sani(double *stress, int ntens, double *asv1, int nasvy,
  double *asv2, int nasvz, double *ddsdde, const double *y, int nydim,
  const double *z, double pore, double depsv_np1, const double *parms,
  int nparms, const double *DD)
{
  double bulk_w = parms[16];
  int i, j;
  (void)nydim;

  pore = pore + bulk_w*depsv_np1;

  for (i = 0; i < ntens; i++) {
    if (i < 3) stress[i] = -y[i] - pore;
    else stress[i] = -y[i];
  }
  for (i = 0; i < nasvy; i++) asv1[i] = y[6+i];
  for (i = 0; i < nasvz; i++) asv2[i] = z[i];

  for (j = 0; j < ntens; j++)
    for (i = 0; i < ntens; i++) {
      if ((i < 3) && (j < 3)) ddsdde[i*6+j] = DD[i*6+j] + bulk_w;
      else ddsdde[i*6+j] = DD[i*6+j];
    }
  (void)nparms;
}

/* ------------------------------------------------------------------ */
/* main entry: single integration-point update                         */
/* props[0..18] raw material parameters (see README).                  */
/* statev[0..35] as documented.                                        */
/* ------------------------------------------------------------------ */
void sanisand_umat(double *stress, double *statev, double *ddsdde,
  double *dstran, double dtime, double *props, int nprops, int testing,
  int *error)
{
  double parms[40];
  double sig_n[6], deps_np1[6], depsv_np1, eps_n[6], epsv_n;
  double pore, asv1[NZDIM], asv2[NZDIM];
  double y[NYDIM], y_n[NYDIM], z[NZDIM], z_n[NZDIM];
  double DDtan[36];
  double p_atm, ptshift, phimob, tol_f, tol_f_test, p_thres;
  double dtsub, norm_D2, norm_D;
  double alphayield[6], avoid, pp, apsi, aec;
  double sig_np1[6], sdev[6], I1, alpha[6], tau[6];
  double gth, qq, cos3t, etanorm, sinphinorm;
  double zero = 0.0;
  double PI = PI_SANISAND;
  int nasvy = 14, nasvz = 14, nyact, nzact;
  int nfev, error_RKF, i, mario_DT_test = 0;
  int plastic_int = 0;
  double DDdummy;

  for (i = 0; i < nprops && i < 40; i++) parms[i] = props[i];
  check_parms_DM(props, parms, nprops);

  tol_f = 1.0e-6;            /* yield-function tolerance (Fortran tol_f) */
  if (testing == 1) { tol_f_test = 1.0e-2; }   /* RKF err_tol on first step */
  else { tol_f_test = 1.0e-3; }                /* RKF err_tol afterwards */
  *error = 0;

  /* ndi must be 3 */
  (void)nyact; (void)nzact;


  p_atm = parms[0];
  p_thres = 0.000000001*p_atm;

  nyact = 6 + nasvy;
  nzact = nasvz;

  pore = statev[28];
  ptshift = parms[17]*parms[0];
  for (i = 0; i < 3; i++) stress[i] = stress[i] - ptshift;

  /* move stress and strain to SM convention */
  move_sig(stress, 6, pore, sig_n);
  move_eps(dstran, 6, deps_np1, &depsv_np1);
  move_eps(dstran, 6, eps_n, &epsv_n);

  norm_D2 = dot_vect(2, deps_np1, deps_np1, 6);
  norm_D = sqrt(norm_D2);

  if (statev[6] < 0.001) {
    double ddum;
    deviator(sig_n, alphayield, &ddum, &pp);
    avoid = 0;
    if (parms[18] <= 5.0) {
      avoid = parms[18];
    } else if (parms[18] > 5.0) {
      apsi = parms[18] - 10.0;
      aec = parms[1] - parms[2]*pow(pp/parms[0], parms[3]);
      avoid = aec + apsi;
    }
    statev[6] = avoid;
    for (i = 0; i < 6; i++) {
      statev[i] = alphayield[i]/pp;
      statev[i+14] = alphayield[i]/pp;
    }
  }

  for (i = 0; i < nasvy; i++) asv1[i] = statev[i];
  for (i = 0; i < nasvz; i++) asv2[i] = statev[14+i];

  dtsub = statev[32];
  if ((dtsub <= zero) || (dtsub > dtime)) dtsub = dtime;

  /* y_n = initial configuration */
  iniyz(y, NYDIM, z, NZDIM, asv1, nasvy, asv2, nasvz, sig_n, 6);
  push(y, y_n, NYDIM);
  push(z, z_n, NZDIM);

  nfev = 0;

  /* integrate */
  rkf23_upd_DM(y, z, nyact, nasvy, nasvz,
    tol_f_test, 50000, 1.0e-18,
    deps_np1, parms, nprops, &nfev, 0, &mario_DT_test, error,
    tol_f, 0, 1, p_thres, &plastic_int);
  /* NOTE: check_ff=0, drcor=1, plastic passed as 0 here; the reference
     passes plastic by value (common block) - see notes. */

  if (*error == 3) {
    push(y_n, y, NYDIM);
    *error = 0;
  }

  if (dtsub <= 0.0) dtsub = 0;
  else if (dtsub >= dtime) dtsub = dtime;
  statev[32] = dtsub;
  statev[33] = (double)nfev;
  *error = 0;

  /* tangent */
  tang_stiff(y, z, nyact, nasvy, nasvz, parms, nprops, DDtan, 1,
    error, tol_f, 0, 1, p_thres, &(int){0});
  if (*error != 0) { /* fall back to elastic */ }

  error_RKF = 0;
  check_RKF_DM(&error_RKF, y, nyact, nasvy, parms, nprops);
  if (error_RKF != 0) {
    push(y_n, y, NYDIM);
  }
  solout_sani(stress, 6, asv1, nasvy, asv2, nasvz, ddsdde, y, NYDIM, z,
    pore, depsv_np1, parms, nprops, DDtan);

  for (i = 0; i < nasvy; i++) statev[i] = asv1[i];
  for (i = 0; i < nasvz; i++) statev[14+i] = asv2[i];

  for (i = 0; i < 6; i++) sig_np1[i] = y[i];
  deviator(sig_np1, sdev, &I1, &pp);
  { double I2t, I3t; inv_sig_full(sig_np1, &pp, &qq, &cos3t, &I1, &I2t, &I3t); }
  statev[28] = pore;
  statev[29] = pp;
  statev[30] = qq;
  statev[31] = cos3t;

  for (i = 0; i < 6; i++) alpha[i] = y[6+i];
  for (i = 0; i < 6; i++) tau[i] = sdev[i] - pp*alpha[i];
  {
    double cM = parms[5]/parms[4];
    double dgdth;
    lode_DM(tau, cM, &cos3t, &gth, &dgdth);
    etanorm = gth*qq/pp;
    sinphinorm = 3*etanorm/(6 + etanorm);
    statev[32] = asin(sinphinorm)*180.0/PI;
    statev[33] = (double)nfev;
  }

  for (i = 0; i < 3; i++) stress[i] = stress[i] + ptshift;

  (void)phimob; (void)DDdummy;
}
