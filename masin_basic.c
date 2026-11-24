/*
 * masin_basic.c - Mašín (2005) five-parameter clay hypoplasticity.
 *
 * This is deliberately separate from masin.c.  masin.c is the HCEA
 * implementation used by group_materi_plasti_hypo_masin_clay; the basic
 * group uses Appendix C, equations (73)--(85), of Mašín (2005).
 */

#include <math.h>

#define MB_PI 3.14159265358979323846

static double mb_dot(const double *a, const double *b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] +
    2.0*(a[3]*b[3] + a[4]*b[4] + a[5]*b[5]);
}

static double mb_det(const double *s)
{
  return s[0]*(s[1]*s[2] - s[5]*s[5]) -
    s[3]*(s[3]*s[2] - s[5]*s[4]) +
    s[4]*(s[3]*s[5] - s[1]*s[4]);
}

/* Appendix C rate equation.  Stress and strain are tensorial Voigt values:
 * [11,22,33,12,13,23], with compression negative.
 *
 * The optional structure (sensitivity) extension of the Professional manual
 * "Masin law" section is applied here.  props[7..9] = k, A, sf and statev[13]
 * = s (hyhis4).  When s < 1 the structure is inactive (Si = 1, ṡ = 0), which
 * is the basic five-parameter model and leaves the pre-existing path intact. */
static int mb_rate(const double *stress, double e, double s, const double *d,
  const double *p, double *dstress, double *dsens)
{
  double phi, lambda, kappa, nvoid, r;
  double k, A, sf, si, structure_on;
  double tr, a, alpha, c1, c2, den, fs, fd;
  double that[6], tdev[6], m[6], larg[6], nn[6], ddev[6];
  double i1, i2, i3, y, tanpsi, cos3, f, mnorm;
  double ts2, ts3, farg, thatdot, that2, dnorm, trd, depd;
  int i;

  phi = p[0]*MB_PI/180.0;
  lambda = p[2]; kappa = p[3]; nvoid = p[4]; r = p[5];
  k = p[7]; A = p[8]; sf = p[9];
  tr = stress[0] + stress[1] + stress[2];
  if (lambda <= 0.0 || kappa <= 0.0 || r <= 0.0 || tr >= -1.e-12 || e <= -0.999999)
    return 0;

  /* Eq. (13): sensitivity factor Si = (s - k(s - sf))/s.  Inactive structure
   * (s < 1, e.g. no structure record and statev[13] = 0) falls back to Si=1. */
  structure_on = (s >= 1.0);
  si = 1.0;
  if (structure_on) {
    si = (s - k*(s - sf))/s;
    if (si < 1.e-12) si = 1.e-12;
  }

  a = sqrt(3.0)*(3.0-sin(phi))/(2.0*sqrt(2.0)*sin(phi));
  /* Eq. (18): alpha uses kappa*Si (equivalently lambda_star/Si). */
  if (structure_on)
    alpha = log((lambda-kappa*si)/(lambda+kappa*si)*(3.0+a*a)/(a*sqrt(3.0)))/log(2.0);
  else
    alpha = log((lambda-kappa)/(lambda+kappa)*(3.0+a*a)/(a*sqrt(3.0)))/log(2.0);
  den = 3.0+a*a-pow(2.0, alpha)*a*sqrt(3.0);
  if (fabs(den) < 1.e-14) return 0;
  /* Eq. (19): c1 has an extra 1/Si with respect to Eq. (85). */
  c1 = 2.0*den/(9.0*r*si);
  c2 = 1.0 + (1.0-c1)*3.0/(a*a);
  for (i=0; i<6; ++i) that[i] = stress[i]/tr;
  for (i=0; i<6; ++i) tdev[i] = that[i];
  tdev[0] -= 1.0/3.0; tdev[1] -= 1.0/3.0; tdev[2] -= 1.0/3.0;

  i1 = tr;
  i2 = 0.5*(mb_dot(stress, stress)-i1*i1);
  i3 = mb_det(stress);
  if (i3 >= -1.e-18) return 0;
  y = (sqrt(3.0)*a/(3.0+a*a)-1.0) *
    ((i1*i2+9.0*i3)*(1.0-sin(phi)*sin(phi))/(8.0*i3*sin(phi)*sin(phi))) +
    sqrt(3.0)*a/(3.0+a*a);

  ts2 = mb_dot(tdev, tdev);
  tanpsi = sqrt(3.0*ts2);
  /* tr(T_hat^3) expansion: sum(aii^3) + 3*a12^2*(a11+a22)
   * + 3*a13^2*(a11+a33) + 3*a23^2*(a22+a33) + 6*a12*a13*a23. */
  ts3 = tdev[0]*(tdev[0]*tdev[0]+3.0*tdev[3]*tdev[3]+3.0*tdev[4]*tdev[4]) +
    tdev[1]*(tdev[1]*tdev[1]+3.0*tdev[3]*tdev[3]+3.0*tdev[5]*tdev[5]) +
    tdev[2]*(tdev[2]*tdev[2]+3.0*tdev[4]*tdev[4]+3.0*tdev[5]*tdev[5]) +
    6.0*tdev[3]*tdev[4]*tdev[5];
  cos3 = ts2 < 1.e-24 ? 1.0 : -sqrt(6.0)*ts3/pow(ts2, 1.5);
  if (cos3 > 1.0) cos3 = 1.0;
  if (cos3 < -1.0) cos3 = -1.0;
  farg = 0.125*tanpsi*tanpsi +
    (2.0-tanpsi*tanpsi)/(2.0+sqrt(2.0)*tanpsi*cos3);
  if (farg <= 0.0) return 0;
  f = sqrt(farg) - tanpsi/(2.0*sqrt(2.0));
  if (f <= 1.e-14) return 0;

  that2 = mb_dot(that, that);
  for (i=0; i<6; ++i)
    m[i] = -(a/f)*(that[i] + tdev[i] - that[i]/3.0 *
      ((6.0*that2-1.0)/((f/a)*(f/a)+that2)));
  mnorm = sqrt(mb_dot(m, m));
  if (mnorm <= 1.e-14) return 0;
  for (i=0; i<6; ++i) larg[i] = -y*m[i]/mnorm;

  /* N = L : (-Y m / |m|), Eq. (74). */
  thatdot = mb_dot(that, larg);
  for (i=0; i<6; ++i) nn[i] = 3.0*(c1*larg[i] + c2*a*a*that[i]*thatdot);
  /* Eq. (12): fs is multiplied by Si. */
  fs = -si*tr/lambda/den;
  /* Eq. (8) of Masin (2007): the pyknotropy factor replaces p*_e by s*p*_e,
   * i.e. fd = (2p/(s p*_e))^alpha.  The manual Eq. (12) prints the reference
   * (s = 1) form; the structure factor divides the bracket. */
  if (structure_on)
    fd = pow(-2.0*tr/3.0*exp((log(1.0+e)-nvoid)/lambda)/s, alpha);
  else
    fd = pow(-2.0*tr/3.0*exp((log(1.0+e)-nvoid)/lambda), alpha);
  thatdot = mb_dot(that, d);
  dnorm = sqrt(mb_dot(d, d));
  for (i=0; i<6; ++i)
    dstress[i] = fs*(3.0*(c1*d[i] + c2*a*a*that[i]*thatdot) + fd*nn[i]*dnorm);

  /* Eq. (21): sdot = -(k/lambda*)(s-sf) sqrt((trD)^2 +
   * (A/(1-A))(eps_s_dot)^2), with eps_s_dot = sqrt(2/3)|devD|. */
  *dsens = 0.0;
  if (structure_on && A < 1.0) {
    trd = d[0]+d[1]+d[2];
    for (i=0; i<3; ++i) ddev[i] = d[i] - trd/3.0;
    ddev[3] = d[3]; ddev[4] = d[4]; ddev[5] = d[5];
    depd = sqrt(trd*trd + A/(1.0-A)*(2.0/3.0)*mb_dot(ddev, ddev));
    *dsens = -(k/lambda)*(s - sf)*depd;
  }
  return 1;
}

static int mb_integrate(const double *stress0, double e0, double s0,
  const double *d, const double *p, double *stress, double *e, double *s)
{
  /* Third-order Runge--Kutta gives a finite-increment integration while
   * retaining the Appendix-C rate law as the sole constitutive equation. */
  const int nsub = 32;
  double st[6], s2[6], s3[6], k1[6], k2[6], k3[6], ds[6];
  double ee, e2, e3, de1, de2, de3, trd, h;
  double ss, ss2, ss3, g1, g2, g3, dg1, dg2, dg3;
  int i, k;
  for (i=0; i<6; ++i) { st[i]=stress0[i]; ds[i]=d[i]/nsub; }
  ee=e0; ss=s0; h=1.0/(double)nsub; trd=d[0]+d[1]+d[2];
  for (k=0; k<nsub; ++k) {
    if (!mb_rate(st, ee, ss, d, p, k1, &g1)) return 0;
    de1=(1.0+ee)*trd; dg1=g1;
    for (i=0; i<6; ++i) s2[i]=st[i]+0.5*h*k1[i];
    e2=ee+0.5*h*de1; ss2=ss+0.5*h*dg1;
    if (!mb_rate(s2, e2, ss2, d, p, k2, &g2)) return 0;
    de2=(1.0+e2)*trd; dg2=g2;
    for (i=0; i<6; ++i) s3[i]=st[i]-h*k1[i]+2.0*h*k2[i];
    e3=ee-h*de1+2.0*h*de2; ss3=ss-h*dg1+2.0*h*dg2;
    if (!mb_rate(s3, e3, ss3, d, p, k3, &g3)) return 0;
    de3=(1.0+e3)*trd; dg3=g3;
    for (i=0; i<6; ++i) st[i]+=h*(k1[i]+4.0*k2[i]+k3[i])/6.0;
    ee += h*(de1+4.0*de2+de3)/6.0;
    ss += h*(dg1+4.0*dg2+dg3)/6.0;
  }
  for (i=0; i<6; ++i) stress[i]=st[i];
  *e=ee; *s=ss;
  return 1;
}

void masin_basic_umat(double *stress, double *statev, double *ddsdde,
  double *dstran, double dtime, double *props, int nprops, int testing,
  int *error)
{
  double s0[6], sp[6], sm[6], ep, em, dp[6], dm[6], eps, e0, sens0;
  double sensp, sensm;
  int i, j;
  (void)dtime; (void)nprops; (void)testing;
  *error=0;
  for (i=0; i<6; ++i) s0[i]=stress[i];
  e0=statev[6];
  sens0=statev[13];
  if (!mb_integrate(s0, e0, sens0, dstran, props, stress, &statev[6],
      &statev[13])) {
    *error=10;
    return;
  }
  /* A numerical algorithmic tangent is required by the global Newton loop.
   * It differentiates the same finite update, not the HCEA tangent. */
  eps=1.e-7;
  for (i=0; i<36; ++i) ddsdde[i]=0.0;
  for (j=0; j<6; ++j) {
    for (i=0; i<6; ++i) { dp[i]=dstran[i]; dm[i]=dstran[i]; }
    dp[j]+=eps; dm[j]-=eps;
    if (!mb_integrate(s0, e0, sens0, dp, props, sp, &ep, &sensp) ||
        !mb_integrate(s0, e0, sens0, dm, props, sm, &em, &sensm)) continue;
    for (i=0; i<6; ++i) ddsdde[i*6+j]=(sp[i]-sm[i])/(2.0*eps);
  }
}
