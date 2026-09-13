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
 * [11,22,33,12,13,23], with compression negative. */
static int mb_rate(const double *stress, double e, const double *d,
  const double *p, double *dstress)
{
  double phi, lambda, kappa, nvoid, r;
  double tr, a, alpha, c1, c2, den, fs, fd;
  double that[6], tdev[6], m[6], larg[6], nn[6];
  double i1, i2, i3, y, tanpsi, cos3, f, mnorm;
  double ts2, ts3, farg, thatdot, that2, dnorm;
  int i;

  phi = p[0]*MB_PI/180.0;
  lambda = p[2]; kappa = p[3]; nvoid = p[4]; r = p[5];
  tr = stress[0] + stress[1] + stress[2];
  if (lambda <= 0.0 || kappa <= 0.0 || r <= 0.0 || tr >= -1.e-12 || e <= -0.999999)
    return 0;

  a = sqrt(3.0)*(3.0-sin(phi))/(2.0*sqrt(2.0)*sin(phi));
  alpha = log((lambda-kappa)/(lambda+kappa)*(3.0+a*a)/(a*sqrt(3.0)))/log(2.0);
  den = 3.0+a*a-pow(2.0, alpha)*a*sqrt(3.0);
  if (fabs(den) < 1.e-14) return 0;
  c1 = 2.0*den/(9.0*r);
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
  fs = -tr/lambda/den;
  fd = pow(-2.0*tr/3.0*exp((log(1.0+e)-nvoid)/lambda), alpha);
  thatdot = mb_dot(that, d);
  dnorm = sqrt(mb_dot(d, d));
  for (i=0; i<6; ++i)
    dstress[i] = fs*(3.0*(c1*d[i] + c2*a*a*that[i]*thatdot) + fd*nn[i]*dnorm);
  return 1;
}

static int mb_integrate(const double *stress0, double e0, const double *d,
  const double *p, double *stress, double *e)
{
  /* Third-order Runge--Kutta gives a finite-increment integration while
   * retaining the Appendix-C rate law as the sole constitutive equation. */
  const int nsub = 32;
  double s[6], s2[6], s3[6], k1[6], k2[6], k3[6], ds[6];
  double ee, e2, e3, de1, de2, de3, trd, h;
  int i, k;
  for (i=0; i<6; ++i) { s[i]=stress0[i]; ds[i]=d[i]/nsub; }
  ee=e0; h=1.0/(double)nsub; trd=d[0]+d[1]+d[2];
  for (k=0; k<nsub; ++k) {
    if (!mb_rate(s, ee, d, p, k1)) return 0;
    de1=(1.0+ee)*trd;
    for (i=0; i<6; ++i) s2[i]=s[i]+0.5*h*k1[i];
    e2=ee+0.5*h*de1;
    if (!mb_rate(s2, e2, d, p, k2)) return 0;
    de2=(1.0+e2)*trd;
    for (i=0; i<6; ++i) s3[i]=s[i]-h*k1[i]+2.0*h*k2[i];
    e3=ee-h*de1+2.0*h*de2;
    if (!mb_rate(s3, e3, d, p, k3)) return 0;
    de3=(1.0+e3)*trd;
    for (i=0; i<6; ++i) s[i]+=h*(k1[i]+4.0*k2[i]+k3[i])/6.0;
    ee += h*(de1+4.0*de2+de3)/6.0;
  }
  for (i=0; i<6; ++i) stress[i]=s[i];
  *e=ee;
  return 1;
}

void masin_basic_umat(double *stress, double *statev, double *ddsdde,
  double *dstran, double dtime, double *props, int nprops, int testing,
  int *error)
{
  double s0[6], sp[6], sm[6], ep, em, dp[6], dm[6], eps, e0;
  int i, j;
  (void)dtime; (void)nprops; (void)testing;
  *error=0;
  for (i=0; i<6; ++i) s0[i]=stress[i];
  e0=statev[6];
  if (!mb_integrate(s0, e0, dstran, props, stress, &statev[6])) {
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
    if (!mb_integrate(s0, e0, dp, props, sp, &ep) ||
        !mb_integrate(s0, e0, dm, props, sm, &em)) continue;
    for (i=0; i<6; ++i) ddsdde[i*6+j]=(sp[i]-sm[i])/(2.0*eps);
  }
}
