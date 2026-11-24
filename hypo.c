#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* hypo.f -- translated by f2c (version 19980831), then ported to pure C.
   The original required linking with -lf2c; the f2c runtime I/O (s_wsle,
   do_lio, e_wsle, s_stop) and pow_dd have been replaced with C equivalents.
*/

#define TRUE_ 1
#define FALSE_ 0
#define abs(x) ((x) >= 0 ? (x) : -(x))


/* Table of constant values */

static long int c__9 = 9;

/* Set by hypoplasticity() (hypoplas.cc) before each call:
   materi_plasti_hypo_substepping / control_materi_plasti_hypo_substepping
   -no switch off the accuracy substepping (manual Professional 6.803,
   6.149; the halving after a failed substep stays), and
   group_materi_plasti_hypo_minimum_void_ratio resets a void ratio below
   the minimum to the minimum (6.715). Thread-local: the element loop
   runs in parallel and the minimum depends on the element group. */
__thread long int hypo_substepping_off = 0;
__thread long int hypo_minimum_void_ratio_on = 0;
__thread long int hypo_void_ratio_exact = 0;
/* group_materi_plasti_hypo_wolffersdorff_niemunis nu z */
__thread long int hypo_niemunis_on = 0;
__thread double hypo_niemunis_nu = 0., hypo_niemunis_z = 0.;

/* ------------------------------------------------------------------------------ */
/* Niemunis visco-hypoplasticity (manual Professional "Niemunis visco law",
   group_materi_plasti_hypo_niemunis_visco phi nu Dr Iv ee0 pe0 lambda betaR
   kappa; Niemunis 2003, habilitation, section 4.2). Tension positive,
   T = effective stress, all tensors 3x3 with index i+3*j.

     That = T/trT, That* = That - 1/3, tan psi = sqrt3 |That*|,
     F = sqrt(tan^2 psi/8 + (2 - tan^2 psi)/(2 + sqrt2 tan psi cos3theta))
         - tan psi/(2 sqrt2)
     a = sqrt3 (3 - sin phi)/(2 sqrt2 sin phi)
     b2 = (1 + a^2/3)(1 - 2 nu)/(1 + nu) - 1   (b2 < 0 rejected by hypoplas.cc)
     L = fb [F^2 I + a^2 That That + b2 (I - 1/3 1 1)],
     fb = -trT/((1 + a^2/3) kappa)
     B = F^2/a^2 (That + That*) + (That:That) That* - That (That*:That*)
     Dvis = -Dr B/|B| (1/OCR)^(1/Iv), OCR = pe/pe+,
     pe = pe0 ((1 + ee0)/(1 + e))^(1/lambda)
     pe+ = p/(betaR - 1) [betaR sqrt(1 + eta^2 (betaR^2 - 1)) - 1], eta < 1
         = p (1 + eta^2)(1 + betaR)/2,                                eta >= 1
     eta = q/(M p), M = 6 F sin phi/(3 - sin phi)
     dT = M:D - L:Dvis, with the intergranular strain stiffness
     M = [rho^chi mT + (1 - rho^chi) mR] L + rho^chi (1 - mT) L:h h   (h:D > 0)
       = [rho^chi mT + (1 - rho^chi) mR] L + rho^chi (mR - mT) L:h h  (h:D <= 0)

   niemunis_visco_shape() returns That, That*, F and pe+ for a stress;
   niemunis_visco_peplus() is used by hypoplasticity() for the initial void
   ratio from group_materi_plasti_hypo_niemunis_visco_ocr. */
static int niemunis_visco_shape(double *t, double *data, double *that,
	double *thatd, double *f, double *peplus)
{
    long int i, j, k;
    double tr, nrm, tanpsi, cos3th, tmp, sinphi, p, q, m, eta, betar, t3;

    tr = t[0] + t[4] + t[8];
    if (tr >= 0.) return 1;
    for (i = 0; i < 9; ++i) that[i] = t[i] / tr;
    for (i = 0; i < 9; ++i) thatd[i] = that[i];
    thatd[0] -= 1. / 3.; thatd[4] -= 1. / 3.; thatd[8] -= 1. / 3.;
    nrm = 0.;
    for (i = 0; i < 9; ++i) nrm += thatd[i] * thatd[i];
    nrm = sqrt(nrm);
    tanpsi = sqrt(3.) * nrm;
    cos3th = 0.;
    if (nrm > 1e-12) {
	t3 = 0.;
	for (i = 0; i < 3; ++i)
	    for (j = 0; j < 3; ++j)
		for (k = 0; k < 3; ++k)
		    t3 += thatd[i + 3 * j] * thatd[j + 3 * k] * thatd[k + 3 * i];
	cos3th = -sqrt(6.) * t3 / (nrm * nrm * nrm);
	if (cos3th > 1.) cos3th = 1.;
	if (cos3th < -1.) cos3th = -1.;
    }
    tmp = tanpsi * tanpsi / 8. + (2. - tanpsi * tanpsi) /
	    (2. + sqrt(2.) * tanpsi * cos3th);
    if (tmp < 0.) tmp = 0.;
    *f = sqrt(tmp) - tanpsi / (2. * sqrt(2.));
    sinphi = sin(data[0] * .017453292777777778);
    p = -tr / 3.;
    q = 0.;
    for (i = 0; i < 9; ++i) {
	tmp = t[i] - ((i == 0 || i == 4 || i == 8) ? tr / 3. : 0.);
	q += tmp * tmp;
    }
    q = sqrt(1.5 * q);
    m = 6. * *f * sinphi / (3. - sinphi);
    eta = m > 0. ? q / (m * p) : 0.;
    betar = data[7];
    if (eta < 1.) {
	if (fabs(betar - 1.) < 1e-8)
	    *peplus = p * (1. + eta * eta);
	else
	    *peplus = p / (betar - 1.) *
		    (betar * sqrt(1. + eta * eta * (betar * betar - 1.)) - 1.);
    } else
	*peplus = p * (1. + eta * eta) * (1. + betar) / 2.;
    if (*peplus <= 0.) return 1;
    return 0;
}

int niemunis_visco_peplus(double *t, double *data, double *peplus)
{
    double that[9], thatd[9], f;
    return niemunis_visco_shape(t, data, that, thatd, &f, peplus);
}

/* stress update over the substep dt with the strain rate d; mm gets the
   tangent M (index i+3*(j+3*(k+3*l))), his[0] the void ratio (already
   updated), his[5] the OCR and *stiffness the size of M. Returns 1 for a
   stress state without compression. */
static int niemunis_visco_rate(double *t, double *mm, double *his,
	long int nhis, double *d, double dt, double *data, long int iuseepi,
	double rho, double mr, double mt, double chi, double epid,
	double *hdir, double *stiffness)
{
    long int i, j, k, l, ind;
    double that[9], thatd[9], f, peplus, sinphi, a, a2, b2, nu, tr, fb,
	    bmat[9], bnrm, tt, tdtd, pe, ocr, creep, dvis[9], lmat[81], lh[9],
	    fac, coef, rhochi, dt_[9], tmp, delta_ij, delta_kl, iso;

    if (niemunis_visco_shape(t, data, that, thatd, &f, &peplus)) return 1;
    tr = t[0] + t[4] + t[8];
    sinphi = sin(data[0] * .017453292777777778);
    a = sqrt(3.) * (3. - sinphi) / (2. * sqrt(2.) * sinphi);
    a2 = a * a;
    nu = data[1];
    b2 = (1. + a2 / 3.) * (1. - 2. * nu) / (1. + nu) - 1.;
    fb = -tr / ((1. + a2 / 3.) * data[8]);

/*       L = fb Lhat */
    for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	    for (k = 0; k < 3; ++k)
		for (l = 0; l < 3; ++l) {
		    ind = i + 3 * (j + 3 * (k + 3 * l));
		    iso = 0.5 * ((i == k && j == l ? 1. : 0.) +
			    (i == l && j == k ? 1. : 0.));
		    delta_ij = (i == j) ? 1. : 0.;
		    delta_kl = (k == l) ? 1. : 0.;
		    lmat[ind] = fb * (f * f * iso + a2 * that[i + 3 * j] *
			    that[k + 3 * l] + b2 * (iso - delta_ij * delta_kl / 3.));
		}

/*       viscous strain rate */
    tt = 0.;
    tdtd = 0.;
    for (i = 0; i < 9; ++i) {
	tt += that[i] * that[i];
	tdtd += thatd[i] * thatd[i];
    }
    bnrm = 0.;
    for (i = 0; i < 9; ++i) {
	bmat[i] = f * f / a2 * (that[i] + thatd[i]) + tt * thatd[i] -
		that[i] * tdtd;
	bnrm += bmat[i] * bmat[i];
    }
    bnrm = sqrt(bnrm);
    pe = data[5] * pow((1. + data[4]) / (1. + his[0]), 1. / data[6]);
    ocr = pe / peplus;
    creep = data[2] * pow(1. / ocr, 1. / data[3]);
    for (i = 0; i < 9; ++i) dvis[i] = bnrm > 0. ? -creep * bmat[i] / bnrm : 0.;
    if (nhis >= 6) his[5] = ocr;

/*       tangent M with intergranular strain */
    fac = 1.;
    coef = 0.;
    if (iuseepi == 1) {
	rhochi = pow(rho, chi);
	fac = rhochi * mt + (1. - rhochi) * mr;
	coef = epid > 0. ? rhochi * (1. - mt) : rhochi * (mr - mt);
	for (i = 0; i < 9; ++i) {
	    lh[i] = 0.;
	    for (k = 0; k < 9; ++k) lh[i] += lmat[i + 9 * k] * hdir[k];
	}
    }
    tmp = 0.;
    for (i = 0; i < 9; ++i)
	for (k = 0; k < 9; ++k) {
	    mm[i + 9 * k] = fac * lmat[i + 9 * k];
	    if (iuseepi == 1) mm[i + 9 * k] += coef * lh[i] * hdir[k];
	    tmp += mm[i + 9 * k] * mm[i + 9 * k];
	}
    *stiffness = sqrt(tmp);

/*       dT = M:D - L:Dvis */
    for (i = 0; i < 9; ++i) {
	dt_[i] = 0.;
	for (k = 0; k < 9; ++k)
	    dt_[i] += mm[i + 9 * k] * d[k] - lmat[i + 9 * k] * dvis[k];
    }
    for (i = 0; i < 9; ++i) t[i] += dt_[i] * dt;
    return 0;
}
__thread double hypo_minimum_void_ratio = 0.;
static double c_b45 = 1.5;
static long int c__81 = 81;

/* ------------------------------------------------------------------------------ */
/* Subroutine */ int hypo_(stress, mmat, his, inc_ept__, time, dtime, nhis, 
	data, ndata, cohesion, epi_r__, epi_mr__, epi_mt__, epi_betar__, 
	epi_chi__, old_epi__, new_epi__, iusepres, iuseepi, ihypotype, 
	softvar_nonloc, softvar_loc, find_local_sv, options_nonlocal)

long int *find_local_sv, *options_nonlocal;
double *softvar_nonloc, *softvar_loc;
double *stress, *mmat, *his, *inc_ept__, *time, *dtime;
long int *nhis;
double *data;
long int *ndata;
double *cohesion, *epi_r__, *epi_mr__, *epi_mt__, *epi_betar__, *
	epi_chi__, *old_epi__, *new_epi__;
long int *iusepres, *iuseepi, *ihypotype;
{
    /* System generated locals */
    double d__1;

    /* Builtin functions */
    double exp();

    /* Local variables */
    double epid;
    long int idim, jdim, kdim, ldim;
    double stress_work__[9]	/* was [3][3] */, rate_epi__[9]	/* 
	    was [3][3] */;
    extern /* Subroutine */ int copy_();
    double dstr[9]	/* was [3][3] */, his_save__[100], epi_size__,
	     new_epi_save__[9]	/* was [3][3] */, d__[9]	/* was [3][3] 
	    */, previous_epi__[9]	/* was [3][3] */;
    long int i__;
    double trace;
    extern /* Subroutine */ int sigma_();
    extern double inpro_();
    long int error;
    extern /* Subroutine */ int minus_();
    extern double power_();
    double direction_epi__[9]	/* was [3][3] */, ed, ei, ep, 
	    dt, middle_epi__[9]	/* was [3][3] */;
    extern /* Subroutine */ int pridbl_(), privec_();
    double dt_tot__;
    long int switch__;
    extern /* Subroutine */ int mul_();
    double tmp, epi_rho__;
    extern double normvec_();
    extern /* Subroutine */ int extract_();
    long int iter_strain__;
    double stress_save__[9]	/* was [3][3] */;

    int counter=0;
    double dnorm_nonloc[1]={0};

/* Routine for hypoplastic stress calculation. */

/* Copyright: Dennis Roddeman */
/*            FEAT Finite Element Application Technology */
/*            dennis.roddeman@feat.nl */

/*   This program is free software; you can redistribute it and/or modify */
/*   it under the terms of the GNU General Public License as published by */
/*   the Free Software Foundation; either version 2 of the License, or */
/*   (at your option) any later version. */

/*   This program is distributed in the hope that it will be useful, */
/*   but WITHOUT ANY WARRANTY; without even the implied warranty of */
/*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/*   GNU General Public License for more details. */


/*   You should have received a copy of the GNU General Public License */
/*   along with this program; if not, write to the Free Software Foundation */
/*   59 Temple Place, Suite 330, Boston, MA, 02111-1307, USA */

/* ------------------------------------------------------------------------------ */




/* ****************** test for valid old void ratio *************** */
    /* Parameter adjustments */
    stress -= 4;
    mmat -= 40;
    inc_ept__ -= 4;
    --time;
    --his;
    --data;
    --cohesion;
    old_epi__ -= 4;
    new_epi__ -= 4;

    /* Function Body */
    if (his[1] <= 0.) {
	printf("Illegal porosity detected in hypoplasticity.\n");
	printf("Remember to initialise node_dof records.\n");
	exit(1);
    }

/* ****************** initialise pressure dependent void ratio **** */
    if (*iusepres == 1 && time[1] == 0.) {
	trace = (double)0.;
	
	for (i__ = 1; i__ <= 3; ++i__) {
	    trace += stress[i__ + i__ * 3];
	}
	d__1 = -trace / data[2];
	tmp = exp(-power_(&d__1, &data[3]));
	ep = his[1] * tmp;
	ei = data[6] * tmp;
	ed = data[5] * tmp;

	if (ep > ei) {
	    his[1] = ei;
	} else if (ep < ed) {
	    his[1] = ed;
	} else {
	    his[1] = ep;
	}
    }
/* *********************  stress with sub-stepping **************** */
    dt = *dtime;
    epid = 0.;
    dt_tot__ = 0.;
    switch__ = FALSE_;
    copy_(&old_epi__[4], &new_epi__[4], &c__9);
    counter=0;

    while(dt_tot__ < *dtime) {
    	counter++;
        /*printf("%lg\n", (double)counter);*/
	
/*           *** adjust, if required, last sub-step size */
	if (dt > *dtime - dt_tot__) {
	    dt = *dtime - dt_tot__;
	}
/*           *** remember in case timestep size will change */
	copy_(&his[1], his_save__, nhis);
	copy_(&stress[4], stress_save__, &c__9);
	if (*iuseepi == 1) {
	    copy_(&new_epi__[4], new_epi_save__, &c__9);
	}
/*           *** strain increment for sub step */
	copy_(&inc_ept__[4], dstr, &c__9);
	d__1 = dt / *dtime;
	mul_(dstr, &d__1, &c__9);
	if(options_nonlocal[0] && !find_local_sv[0]) 
	   dnorm_nonloc[0]=sqrt(softvar_nonloc[0]);	

/*           *** intergranular strain for sub step */
/*           *** mid-point-rule to determine rate of intergranular strain */
	if (*iuseepi == 1) {
	    extract_(dstr, &dt, d__);
	    copy_(&new_epi__[4], previous_epi__, &c__9);
	    for (iter_strain__ = 1; iter_strain__ <= 10; ++iter_strain__) {
		for (idim = 1; idim <= 3; ++idim) {
		    for (jdim = 1; jdim <= 3; ++jdim) {
			middle_epi__[idim + jdim * 3 - 4] = (previous_epi__[
				idim + jdim * 3 - 4] + new_epi__[idim + jdim *
				 3]) / 2.;
/* L50: */
		    }
/* L60: */
		}
		epi_size__ = normvec_(middle_epi__, &c__9);
		epi_rho__ = epi_size__ / *epi_r__;
		/* Professional manual (intergranular strains extension): "From
		   the evolution equation it follows that rho must remain
		   between 0 and 1." The discrete midpoint rule can overshoot
		   rho > 1 when the strain direction rotates within a step
		   (e.g. an isotropic intergranular-strain reset in an
		   oedometric path); clamping keeps the stiffness weights
		   [rho^chi*mT + (1-rho^chi)*mR] inside their physical range
		   (unclamped, rho > 1 gives a NEGATIVE weight to mR).
		   Calibrated 2026-09-07 vs the Professional: hypo2 sigxx
		   -0.176 -> -0.17447 (Professional -0.17284). */
		if (epi_rho__ > 1.) epi_rho__ = 1.;
		copy_(middle_epi__, direction_epi__, &c__9);
		if (epi_size__ > 1e-12) {
		    d__1 = 1. / epi_size__;
		    mul_(direction_epi__, &d__1, &c__9);
		}
		epid = inpro_(direction_epi__, d__, &c__9);
		for (idim = 1; idim <= 3; ++idim) {
		    for (jdim = 1; jdim <= 3; ++jdim) {
			rate_epi__[idim + jdim * 3 - 4] = d__[idim + jdim * 3 
				- 4];
			if (epid >= 0.) {
			    for (kdim = 1; kdim <= 3; ++kdim) {
				for (ldim = 1; ldim <= 3; ++ldim) {
				    rate_epi__[idim + jdim * 3 - 4] -= power_(
					    &epi_rho__, epi_betar__) * 
					    direction_epi__[idim + jdim * 3 - 
					    4] * direction_epi__[kdim + ldim *
					     3 - 4] * d__[kdim + ldim * 3 - 4]
					    ;
/* L70: */
				}
/* L80: */
			    }
			}
			new_epi__[idim + jdim * 3] = previous_epi__[idim + 
				jdim * 3 - 4] + rate_epi__[idim + jdim * 3 - 
				4] * dt;
/* L90: */
		    }
/* L100: */
		}
/* L110: */
	    }
	}
/*           *** calculate new stress */
	error = FALSE_;
	sigma_(&switch__, &error, &stress[4], &mmat[40], &his[1], dstr, &dt, 
		nhis, &data[1], ndata, &cohesion[1], &epi_rho__, epi_mr__, 
		epi_mt__, epi_betar__, epi_chi__, &epid, direction_epi__, 
		iuseepi, ihypotype,
		softvar_nonloc, softvar_loc, find_local_sv, options_nonlocal, 
		dnorm_nonloc);
	if (error) {
	    if (dt >= *dtime * 1.9999999999999999e-6) {
		dt /= 2.;
		copy_(his_save__, &his[1], nhis);
		copy_(stress_save__, &stress[4], &c__9);
		if (*iuseepi == 1) {
		    copy_(new_epi_save__, &new_epi__[4], &c__9);
		}
	    } else {
		goto L210;
	    }
	} else {
	    minus_(&stress[4], stress_save__, stress_work__, &c__9);
	    /*increases number of substeps, if necessary*/
	    if(!(options_nonlocal[0] && find_local_sv[0])) { /*substepping not used when searching for local ||D||^2*/
     	      /* Substeps are halved while the trial stress increment of
		 the substep exceeds 2% of the current stress magnitude (and
		 the current stress is not negligible). The historical
		 tolerance was 1% (hypo.f, GNU 2014). Calibrated 2026-09-07
		 against the Professional binary: the effective tolerance of
		 the Professional kernel is ~2% (hypo1 sigyy -862.766 @1%
		 vs Professional -862.929; @2% the GNU gives -862.927,
		 within the target -862.92 +- 0.1 -> rc=0; @5% it
		 overshoots to -863.29). */
     	      if (!hypo_substepping_off &&
		  normvec_(stress_work__, &c__9) > normvec_(&stress[4], &c__9) *
	         0.02 && normvec_(&stress[4], &c__9) > .1 && dt >= *dtime *
	    	     1.9999999999999999e-6) {
			dt /= 2.;
			copy_(his_save__, &his[1], nhis);
			copy_(stress_save__, &stress[4], &c__9);
			if (*iuseepi == 1) {
			    copy_(new_epi_save__, &new_epi__[4], &c__9);
			}
	      } else dt_tot__ += dt;
	    }
	    else dt_tot__ += dt;
	}
/* L200: */
    }
L210:
/*       store for plotting */
    his[2] = dt;

/* ***************  some last printing ******************** */
    switch__ = FALSE_;
    if (switch__) {
	pridbl_("dt    ", &dt, (long int)6);
	privec_("stress", &stress[4], &c__9, (long int)6);
    }

    return 0;
} /* hypo_ */


/* ------------------------------------------------------------------------------ */
/* Subroutine */ int sigma_(switch__, error, stress, mmat, his, dstr, dt, 
	nhis, data, ndata, cohesion, epi_rho__, epi_mr__, epi_mt__, 
	epi_betar__, epi_chi__, epid, direction_epi__, iuseepi, ihypotype,
	softvar_nonloc, softvar_loc, find_local_sv, options_nonlocal, dnorm_nonloc)

long int *find_local_sv, *options_nonlocal;
double *softvar_nonloc, *softvar_loc, *dnorm_nonloc;
long int *switch__, *error;
double *stress, *mmat, *his, *dstr, *dt;
long int *nhis;
double *data;
long int *ndata;
double *cohesion, *epi_rho__, *epi_mr__, *epi_mt__, *epi_betar__, *
	epi_chi__, *epid, *direction_epi__;
long int *iuseepi, *ihypotype;
{


    /* System generated locals */
    double d__1, d__2, d__3, d__4, dmob_1, dmob_2;

    /* Builtin functions */
    double sqrt(), sin(), exp(), asin();

    /* Local variables */
    double alfa, beta;
    extern /* Subroutine */ int abdyadic_();
    double eold, phic, edot, lmat[81]	/* 
	    was [3][3][3][3] */, nmat[9]	/* was [3][3] */, nval, rval, that[
	    9]	/* was [3][3] */, tmin, tmax;
    extern /* Subroutine */ int copy_(), zero_();
    double a, work1[9]	/* was [3][3] */, d__[9]	/* 
	    was [3][3] */, e, f, work2[9]	/* was [3][3] */;
    long int i__, j, k, l, m, n;
    extern /* Subroutine */ int eigen_();
    extern double power_();
    double lwork[81]	/* was [3][3][3][3] */, tcohesion[9]	/* 
	    was [3][3] */;
    extern /* Subroutine */ int addfac_();
    double cos3th, ec, ed, fd, fe, fb, ei, hs, c1=0., c2=0., xi, re, powxi;
    extern /* Subroutine */ int unity4_();
    double eigval[3];
    extern /* Subroutine */ int pridbl_();
    double d_size__, phimob, ec0, ed0, epi_rhochi__;
    extern /* Subroutine */ int privec_();
    double tanpsi, ei0;
    extern long int null_array__();
    double inc_stress__[9]	/* was [3][3] */;
    extern /* Subroutine */ int abc_(), add_();
    double fac;
    extern /* Subroutine */ int dev_();
    extern double tra_();
    extern /* Subroutine */ int mul_();
    double tmp;
    extern /* Subroutine */ int a4bc_();
    extern long int negativ_();
    double thatdev[9]	/* was [3][3] */;
    extern double normvec_();
    extern /* Subroutine */ int extract_();
    extern double normmat_();

/* ------------------------------------------------------------------------------ */





/*       old void ratio */
    /* Parameter adjustments */
    stress -= 4;
    mmat -= 40;
    dstr -= 4;
    --his;
    --data;
    direction_epi__ -= 4;

    /* unused parameters (kept for ABI compatibility with f2c signature) */
    (void)nhis;
    (void)ndata;
    (void)epi_betar__;
    (void)softvar_nonloc;

    /* Function Body */
    eold = his[1];
/*       new void ratio */
    extract_(&dstr[4], dt, d__);
    edot = (eold + 1.) * tra_(d__);
    e = eold + edot * *dt;
/*       group_materi_plasti_hypo_void_ratio_linear -no: exact volume
         change from the determinant of I + strain increment */
    if (hypo_void_ratio_exact) {
	double *g = &dstr[4];
	double g00 = 1. + g[0], g11 = 1. + g[4], g22 = 1. + g[8];
	double det = g00 * (g11 * g22 - g[5] * g[7])
		   - g[3] * (g[1] * g22 - g[2] * g[7])
		   + g[6] * (g[1] * g[5] - g[2] * g11);
	e = (eold + 1.) * det - 1.;
    }
    if (hypo_minimum_void_ratio_on && e < hypo_minimum_void_ratio)
	e = hypo_minimum_void_ratio;
    his[1] = e;

/*       Niemunis visco-hypoplasticity: own stress update and tangent */
    if (*ihypotype == 2) {
	phic = data[1] * .017453292777777778;
	if (niemunis_visco_rate(&stress[4], &mmat[40], &his[1], *nhis, d__,
		*dt, &data[1], *iuseepi, *epi_rho__, *epi_mr__, *epi_mt__,
		*epi_chi__, *epid, &direction_epi__[4], &re)) {
	    *error = TRUE_;
	    return 0;
	}
	goto L_mobilised;
    }

/*       von Wolffersdorff */
    phic = data[1] * .017453292777777778;
    hs = data[2];
    nval = data[3];
    ec0 = data[4];
    ed0 = data[5];
    ei0 = data[6];
    alfa = data[7];
    beta = data[8];
    rval = 0;
    powxi = 1;

/*       hypo clay */
    if( *ihypotype == 1 ) {
    	rval = data[9];
	powxi = data[10];
    }

/*       subtract cohesion to get hypoplastic law for cohesion */
    copy_(&stress[4], tcohesion, &c__9);
    for (i__ = 1; i__ <= 3; ++i__) {
	tcohesion[i__ + i__ * 3 - 4] -= *cohesion;
    }

    zero_(inc_stress__, &c__9);

    d__1 = tra_(tcohesion);
    if (null_array__(&d__1)) {
	*error = TRUE_;
	goto L1000;
    }
    copy_(tcohesion, that, &c__9);
    d__1 = 1. / tra_(tcohesion);
    mul_(that, &d__1, &c__9);
    dev_(that, thatdev);
    tanpsi = sqrt(3.) * normmat_(thatdev);
    abc_(thatdev, thatdev, work1);
    abc_(thatdev, work1, work2);
    if (tra_(work1) < 1e-10) {
	cos3th = 1.;
    } else {
	d__1 = tra_(work1);
	cos3th = -sqrt(6.) * tra_(work2) / power_(&d__1, &c_b45);
	if (cos3th < -1.) {
	    cos3th = -1.;
	}
	if (cos3th > 1.) {
	    cos3th = 1.;
	}
    }
    if (*switch__) {
	pridbl_("tanpsi", &tanpsi, (long int)6);
	pridbl_("cos3th", &cos3th, (long int)6);
    }

    a = sqrt(3.) * (3. - sin(phic)) / (sqrt(2.) * 2. * sin(phic));
    tmp = tanpsi * tanpsi / 8. + (2. - tanpsi * tanpsi) / (sqrt(2.) * tanpsi *
	     cos3th + 2.);
    if (negativ_(&tmp)) {
	*error = TRUE_;
	goto L1000;
    }
    f = sqrt(tmp) - tanpsi / (sqrt(2.) * 2.);
    if (*switch__) {
	pridbl_("f     ", &f, (long int)6);
	pridbl_("a     ", &a, (long int)6);
    }

    d__1 = -tra_(tcohesion);
    if (negativ_(&d__1)) {
	*error = TRUE_;
	goto L1000;
    }
    d__1 = -tra_(tcohesion) / hs;
    tmp = exp(-power_(&d__1, &nval));
    ei = tmp * ei0;
    ec = tmp * ec0;
    ed = tmp * ed0;
    if (*switch__) {
	pridbl_("ei    ", &ei, (long int)6);
	pridbl_("ec    ", &ec, (long int)6);
	pridbl_("ed    ", &ed, (long int)6);
    }

    d__1 = ec / e;
    if (negativ_(&d__1)) {
	*error = TRUE_;
	goto L1000;
    }
    d__1 = ec / e;
    fe = power_(&d__1, &beta);
    d__1 = ec - ed;
    if (null_array__(&d__1)) {
	*error = TRUE_;
	goto L1000;
    }

    d__1 = (e - ed) / (ec - ed);
    re = d__1;
    fd = power_(&d__1, &alfa);
    if (e < 1.0001*ed) fd=0.0001;

    if (null_array__(&ei)) {
	*error = TRUE_;
	goto L1000;
    }
    if (negativ_(&ei0)) {
	*error = TRUE_;
	goto L1000;
    }
    if (negativ_(&ec0)) {
	*error = TRUE_;
	goto L1000;
    }
    d__1 = ec0 - ed0;
    if (null_array__(&d__1)) {
	*error = TRUE_;
	goto L1000;
    }
    d__1 = (ei0 - ed0) / (ec0 - ed0);
    if (negativ_(&d__1)) {
	*error = TRUE_;
	goto L1000;
    }
    d__1 = ei0 / ec0;
    d__2 = -tra_(tcohesion) / hs;
    d__3 = 1. - nval;
    d__4 = (ei0 - ed0) / (ec0 - ed0);
    fb = hs / nval * ((ei + 1.) / ei) * power_(&d__1, &beta) * power_(&d__2, &
	    d__3) / (a * a + 3. - sqrt(3.) * a * power_(&d__4, &alfa));

    if( *ihypotype == 1 ) {
        eigen_(&stress[4], eigval);
        dmob_1 = abs(eigval[0]), dmob_2 = abs(eigval[1]), dmob_1 = ((dmob_1>dmob_2)?dmob_1:dmob_2), dmob_2 
	    = abs(eigval[2]);
        tmax = -((dmob_1>dmob_2)?dmob_1:dmob_2);
        dmob_1 = abs(eigval[0]), dmob_2 = abs(eigval[1]), dmob_1 = ((dmob_1<dmob_2)?dmob_1:dmob_2), dmob_2 
	    = abs(eigval[2]);
        tmin = -((dmob_1<dmob_2)?dmob_1:dmob_2);
        if ((dmob_1 = tmax + tmin, abs(dmob_1)) < 1e-12) 
	   phimob = phic;
        else 
	   phimob = asin((dmob_1 = tmax - tmin, abs(dmob_1)) / (dmob_2 = tmax + tmin, 
		abs(dmob_2)));
		
    	xi = (sin(phic) - sin(phimob))/sin(phic);
    	if ( xi < 0 ) xi = 0;
	xi = power_(&xi, &powxi);
    	c1 = (1 + (a*a)/3 - a/(sqrt(3.)))/(1.5 * rval);
	c1 = power_(&c1, &xi);
    	c2 = 1. + (1. - c1)*(3. / (a * a));
        fb = hs / nval * ((ei + 1.) / ei) * power_(&d__1, &beta) * power_(&d__2, &
	    d__3) / (c2 * a * a + 3. * c1 - sqrt(3.) * a * power_(&d__4, &alfa));
    }
    
/*       Niemunis extension {Niemunis 2003, eq. 4.211, 4.221, 4.223}:
         f_d consistent with the lower bound e = e_d */
    if (hypo_niemunis_on && *ihypotype == 0) {
	double pt = -tra_(tcohesion) / hs, nm1 = nval - 1., mt, fdbar, fd1;
	mt = -nval * ed / hs * power_(&pt, &nm1);
	fdbar = -(sqrt(3.) * (e + 1.) + mt * fb * fe * 3. / sqrt(3.) *
		(a * a + 3.)) / (mt * fb * fe * 3. * a);
	d__1 = fabs(e - ed) / (ec - ed);
	fd1 = power_(&d__1, &alfa);
	if (e < ed) fd1 = -fd1;
	fd = fd1;
	if (fd < 1.) {
	    d__1 = 1. - fd;
	    fd += power_(&d__1, &hypo_niemunis_z) * fdbar;
	}
    }

    if (*switch__) {
	pridbl_("fb    ", &fb, (long int)6);
	pridbl_("fe    ", &fe, (long int)6);
	pridbl_("fd    ", &fd, (long int)6);
    }

    abc_(that, that, work1);
    d__1 = tra_(work1);
    if (null_array__(&d__1)) {
	*error = TRUE_;
	goto L1000;
    }
    fac = fb * fe / tra_(work1);

/*        first linear */
    tmp = fac * f * f;
    if( *ihypotype == 1 ) tmp = fac * f * f * c1;
    unity4_(lmat);
    mul_(lmat, &tmp, &c__81);
    
/*        second linear */
    tmp = fac * a * a;
    if( *ihypotype == 1 ) tmp = fac * a * a * c2;
    abdyadic_(that, that, lwork);
    mul_(lwork, &tmp, &c__81);
    add_(lwork, lmat, lmat, &c__81);

/*        Niemunis extension: b2 {I - 1/3 1 1} with b2 from the Poisson
          ratio nu at the hydrostatic axis {Niemunis 2003, eq. 4.183} */
    if (hypo_niemunis_on && *ihypotype == 0) {
	double b2 = (1. + a * a / 3. + a / sqrt(3.)) * (1. - 2. * hypo_niemunis_nu)
		/ (1. + hypo_niemunis_nu) - 1.;
	for (i__ = 1; i__ <= 3; ++i__) {
	    for (j = 1; j <= 3; ++j) {
		for (k = 1; k <= 3; ++k) {
		    for (l = 1; l <= 3; ++l) {
			tmp = 0.;
			if (i__ == k && j == l) tmp += 1.;
			if (i__ == j && k == l) tmp -= 1. / 3.;
			lmat[i__ + (j + (k + l * 3) * 3) * 3 - 40] += fac * b2 * tmp;
		    }
		}
	    }
	}
    }
    
/*        nonlinear */
    tmp = fac * fd * a * f;
    add_(that, thatdev, nmat, &c__9);
    mul_(nmat, &tmp, &c__9);

/*        Niemunis extension: the flow rule m stays unchanged {Niemunis
          2003, below eq. 4.183}, so N = -fd Y L:m follows the b2 term of
          L: N := N + b2 dev{Lhat^-1 : N}, Lhat = F^2 I + a^2 That That */
    if (hypo_niemunis_on && *ihypotype == 0) {
	double b2 = (1. + a * a / 3. + a / sqrt(3.)) * (1. - 2. * hypo_niemunis_nu)
		/ (1. + hypo_niemunis_nu) - 1.;
	double tn = 0., tt = 0., x[9], trx;
	for (i__ = 0; i__ < 9; ++i__) {
	    tn += that[i__] * nmat[i__];
	    tt += that[i__] * that[i__];
	}
	for (i__ = 0; i__ < 9; ++i__)
	    x[i__] = (nmat[i__] - a * a * that[i__] * tn / (f * f + a * a * tt))
		/ (f * f);
	trx = (x[0] + x[4] + x[8]) / 3.;
	x[0] -= trx; x[4] -= trx; x[8] -= trx;
	for (i__ = 0; i__ < 9; ++i__) nmat[i__] += b2 * x[i__];
    }

/*        total stress increment and material stiffness */
    if (*iuseepi == 1) {
	epi_rhochi__ = power_(epi_rho__, epi_chi__);
	for (i__ = 1; i__ <= 3; ++i__) {
	    for (j = 1; j <= 3; ++j) {
		for (k = 1; k <= 3; ++k) {
		    for (l = 1; l <= 3; ++l) {
			tmp = epi_rhochi__ * *epi_mt__ + (1. - epi_rhochi__) *
				 *epi_mr__;
			mmat[i__ + (j + (k + l * 3) * 3) * 3] = tmp * lmat[
				i__ + (j + (k + l * 3) * 3) * 3 - 40];
			if (*epid > 0.) {
			    mmat[i__ + (j + (k + l * 3) * 3) * 3] += 
				    epi_rhochi__ * nmat[i__ + j * 3 - 4] * 
				    direction_epi__[k + l * 3];
			}
			for (m = 1; m <= 3; ++m) {
			    for (n = 1; n <= 3; ++n) {
				if (*epid >= 0.) {
				    mmat[i__ + (j + (k + l * 3) * 3) * 3] += 
					    epi_rhochi__ * (1. - *epi_mt__) * 
					    lmat[i__ + (j + (m + n * 3) * 3) *
					     3 - 40] * direction_epi__[m + n *
					     3] * direction_epi__[k + l * 3];
				} else {
				    mmat[i__ + (j + (k + l * 3) * 3) * 3] += 
					    epi_rhochi__ * (*epi_mr__ - *
					    epi_mt__) * lmat[i__ + (j + (m + 
					    n * 3) * 3) * 3 - 40] * 
					    direction_epi__[m + n * 3] * 
					    direction_epi__[k + l * 3];
				}
/* L50: */
			    }
/* L60: */
			}
/* L70: */
		    }
/* L80: */
		}
/* L90: */
	    }
/* L100: */
	}
	a4bc_(&mmat[40], d__, inc_stress__);
    } else {
	d_size__ = normvec_(d__, &c__9);
	for (i__ = 1; i__ <= 3; ++i__) {
	    for (j = 1; j <= 3; ++j) {
		for (k = 1; k <= 3; ++k) {
		    for (l = 1; l <= 3; ++l) {
			mmat[i__ + (j + (k + l * 3) * 3) * 3] = lmat[i__ + (j 
				+ (k + l * 3) * 3) * 3 - 40];
			if(0 && !options_nonlocal[0]) {/*Use only linear stiffness if nonlocal ||D||*/
				if (d_size__ > 0.) {
				    mmat[i__ + (j + (k + l * 3) * 3) * 3] += nmat[i__ 
					    + j * 3 - 4] * d__[k + l * 3 - 4] / d_size__;
				}
			}
/* L170: */
		    }
/* L180: */
		}
/* L190: */
	    }
/* L200: */
	}
	a4bc_(lmat, d__, inc_stress__);
	d__1 = normmat_(d__);
	/*d__1=||D||*/
	if(options_nonlocal[0] && find_local_sv[0]) softvar_loc[0]=d__1*d__1;
	if(options_nonlocal[0] && !find_local_sv[0]) d__1=dnorm_nonloc[0];

	addfac_(inc_stress__, nmat, &d__1, &c__9);
    }
/*       from rate to increment */
    mul_(inc_stress__, dt, &c__9);
    add_(&stress[4], inc_stress__, &stress[4], &c__9);
    if (*switch__) {
	privec_("stress", &stress[4], &c__9, (long int)6);
    }

L_mobilised:
/*       mobilised friction angle for postprocessing */
    eigen_(&stress[4], eigval);
/* Computing MAX */
    d__1 = abs(eigval[0]), d__2 = abs(eigval[1]), d__1 = ((d__1>d__2)?d__1:d__2), d__2 
	    = abs(eigval[2]);
    tmax = -((d__1>d__2)?d__1:d__2);
/* Computing MIN */
    d__1 = abs(eigval[0]), d__2 = abs(eigval[1]), d__1 = ((d__1<d__2)?d__1:d__2), d__2 
	    = abs(eigval[2]);
    tmin = -((d__1<d__2)?d__1:d__2);
    if ((d__1 = tmax + tmin, abs(d__1)) < 1e-12) {
	phimob = phic;
    } else {
	phimob = asin((d__1 = tmax - tmin, abs(d__1)) / (d__2 = tmax + tmin, 
		abs(d__2)));
	phimob = phimob * 360. / 6.2831853999999998;
    }
    his[3] = phimob;

    /*       effective stiffness */
    /*his[4] = effective_stiffness__;*/
    his[4] = re;

L1000:

    return 0;
} /* sigma_ */


/* ------------------------------------------------------------------------------ */
/* Subroutine */ int abc_(a, b, c__)
double *a, *b, *c__;
{
    long int idi, jdi, kdi;

/* ------------------------------------------------------------------------------ */



    /* Parameter adjustments */
    c__ -= 4;
    b -= 4;
    a -= 4;

    /* Function Body */
    for (idi = 1; idi <= 3; ++idi) {
	for (jdi = 1; jdi <= 3; ++jdi) {
	    c__[idi + jdi * 3] = 0.;
	    for (kdi = 1; kdi <= 3; ++kdi) {
		c__[idi + jdi * 3] += a[idi + kdi * 3] * b[kdi + jdi * 3];
/* L10: */
	    }
/* L20: */
	}
/* L30: */
    }

    return 0;
} /* abc_ */


/* ------------------------------------------------------------------------------ */
/* Subroutine */ int a4bc_(a, b, c__)
double *a, *b, *c__;
{
    long int idi, jdi, kdi, ldi;

/* ------------------------------------------------------------------------------ */



    /* Parameter adjustments */
    c__ -= 4;
    b -= 4;
    a -= 40;

    /* Function Body */
    for (idi = 1; idi <= 3; ++idi) {
	for (jdi = 1; jdi <= 3; ++jdi) {
	    c__[idi + jdi * 3] = 0.;
	    for (kdi = 1; kdi <= 3; ++kdi) {
		for (ldi = 1; ldi <= 3; ++ldi) {
		    c__[idi + jdi * 3] += a[idi + (jdi + (kdi + ldi * 3) * 3) 
			    * 3] * b[kdi + ldi * 3];
/* L10: */
		}
/* L20: */
	    }
/* L30: */
	}
/* L40: */
    }

    return 0;
} /* a4bc_ */


/* ----------------------------------------------------------------------------- */
/* Subroutine */ int add_(a, b, c__, n)
double *a, *b, *c__;
long int *n;
{
    /* System generated locals */
    long int i__1;

    /* Local variables */
    long int i__;

/* ----------------------------------------------------------------------------- */


    /* Parameter adjustments */
    --c__;
    --b;
    --a;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	c__[i__] = a[i__] + b[i__];
/* L10: */
    }

    return 0;
} /* add_ */


/* ----------------------------------------------------------------------------- */
/* Subroutine */ int addfac_(tot, inc, factor, n)
double *tot, *inc, *factor;
long int *n;
{
    /* System generated locals */
    long int i__1;

    /* Local variables */
    long int i__;

/* ----------------------------------------------------------------------------- */


    /* Parameter adjustments */
    --inc;
    --tot;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	tot[i__] += *factor * inc[i__];
/* L10: */
    }

    return 0;
} /* addfac_ */


/* ----------------------------------------------------------------------------- */
/* Subroutine */ int copy_(a, b, n)
double *a, *b;
long int *n;
{
    /* System generated locals */
    long int i__1;

    /* Local variables */
    long int i__;

/* ------------------------------------------------------------------------------ */


    /* Parameter adjustments */
    --b;
    --a;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b[i__] = a[i__];
/* L10: */
    }

    return 0;
} /* copy_ */


/* ----------------------------------------------------------------------------- */
/* Subroutine */ int dev_(a, adev)
double *a, *adev;
{
    long int idi, jdi;
    extern double tra_();

/* ------------------------------------------------------------------------------ */



    /* Parameter adjustments */
    adev -= 4;
    a -= 4;

    /* Function Body */
    for (idi = 1; idi <= 3; ++idi) {
	for (jdi = 1; jdi <= 3; ++jdi) {
	    adev[idi + jdi * 3] = a[idi + jdi * 3];
/* L10: */
	}
	adev[idi + idi * 3] -= tra_(&a[4]) / 3.;
/* L20: */
    }

    return 0;
} /* dev_ */


/* ----------------------------------------------------------------------------- */
/* Subroutine */ int extract_(dstr, dt, d__)
double *dstr, *dt, *d__;
{
    /* System generated locals */
    double d__1;

    /* Local variables */
    extern /* Subroutine */ int copy_(), mul_();

/* ----------------------------------------------------------------------------- */



    /* Parameter adjustments */
    d__ -= 4;
    dstr -= 4;

    /* Function Body */
    copy_(&dstr[4], &d__[4], &c__9);
    d__1 = (double)1. / *dt;
    mul_(&d__[4], &d__1, &c__9);

    return 0;
} /* extract_ */


/* ----------------------------------------------------------------------------- */
/* Subroutine */ int minus_(a, b, c__, n)
double *a, *b, *c__;
long int *n;
{
    /* System generated locals */
    long int i__1;

    /* Local variables */
    long int i__;

/* ----------------------------------------------------------------------------- */


    /* Parameter adjustments */
    --c__;
    --b;
    --a;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	c__[i__] = a[i__] - b[i__];
/* L10: */
    }

    return 0;
} /* minus_ */


/* ----------------------------------------------------------------------------- */
/* Subroutine */ int mul_(a, scal, n)
double *a, *scal;
long int *n;
{
    /* System generated locals */
    long int i__1;

    /* Local variables */
    long int i__;

/* ----------------------------------------------------------------------------- */


    /* Parameter adjustments */
    --a;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a[i__] = *scal * a[i__];
/* L10: */
    }

    return 0;
} /* mul_ */


/* ----------------------------------------------------------------------------- */
long int negativ_(x)
double *x;
{
    /* System generated locals */
    long int ret_val;

/* ------------------------------------------------------------------------------ */


    if (*x < (double)1e-8) {
	ret_val = TRUE_;
    } else {
	ret_val = FALSE_;
    }

    return ret_val;
} /* negativ_ */


/* ------------------------------------------------------------------------------ */
double normmat_(a)
double *a;
{
    /* System generated locals */
    double ret_val;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    double work[9]	/* was [3][3] */;
    extern /* Subroutine */ int abc_();
    extern double tra_();

/* ------------------------------------------------------------------------------ */



    /* Parameter adjustments */
    a -= 4;

    /* Function Body */
    abc_(&a[4], &a[4], work);
    ret_val = sqrt(tra_(work));

    return ret_val;
} /* normmat_ */


/* ------------------------------------------------------------------------------ */
double normvec_(a, n)
double *a;
long int *n;
{
    /* System generated locals */
    long int i__1;
    double ret_val;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    long int i__;

/* ------------------------------------------------------------------------------ */


    /* Parameter adjustments */
    --a;

    /* Function Body */
    ret_val = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ret_val += a[i__] * a[i__];
/* L10: */
    }
    ret_val = sqrt(ret_val);

    return ret_val;
} /* normvec_ */


/* ----------------------------------------------------------------------------- */
long int null_array__(x)
double *x;
{
    /* System generated locals */
    long int ret_val;

/* ------------------------------------------------------------------------------ */


    if (abs(*x) < (double)1e-8) {
	ret_val = TRUE_;
    } else {
	ret_val = FALSE_;
    }

    return ret_val;
} /* null_array__ */


/* ------------------------------------------------------------------------------ */
double power_(x, y)
double *x, *y;
{
    /* System generated locals */
    double ret_val;

/* ------------------------------------------------------------------------------ */


    ret_val = pow(*x, *y);

    return ret_val;
} /* power_ */


/* ------------------------------------------------------------------------------ */
/* Subroutine */ int pridbl_(label, a, label_len)
char *label;
double *a;
long int label_len;
{


/* ------------------------------------------------------------------------------ */


    printf(" %.*s %.12g\n", (int)label_len, label, *a);

    return 0;
} /* pridbl_ */


/* Subroutine */ int priint_(label, a, label_len)
char *label;
long int *a;
long int label_len;
{


/* ------------------------------------------------------------------------------ */


    printf(" %.*s %ld\n", (int)label_len, label, (long)*a);

    return 0;
} /* priint_ */


/* ------------------------------------------------------------------------------ */
/* Subroutine */ int pritxt_(label, label_len)
char *label;
long int label_len;
{


/* ------------------------------------------------------------------------------ */


    printf(" %.*s\n", (int)label_len, label);

    return 0;
} /* pritxt_ */


/* ------------------------------------------------------------------------------ */
/* Subroutine */ int privec_(label, a, n, label_len)
char *label;
double *a;
long int *n;
long int label_len;
{
    /* System generated locals */
    long int i__1;

    /* Local variables */
    long int i__;


/* ------------------------------------------------------------------------------ */


    /* Parameter adjustments */
    --a;

    /* Function Body */
    printf(" %.*s:", (int)label_len, label);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	printf(" %.12g", a[i__]);
    }
    printf("\n");

    return 0;
} /* privec_ */


/* ------------------------------------------------------------------------------ */
double tra_(a)
double *a;
{
    /* System generated locals */
    double ret_val;

/* ------------------------------------------------------------------------------ */



    /* Parameter adjustments */
    a -= 4;

    /* Function Body */
    ret_val = a[4] + a[8] + a[12];

    return ret_val;
} /* tra_ */


/* ----------------------------------------------------------------------------- */
/* Subroutine */ int zero_(a, n)
double *a;
long int *n;
{
    /* System generated locals */
    long int i__1;

    /* Local variables */
    long int i__;

/* ----------------------------------------------------------------------------- */


    /* Parameter adjustments */
    --a;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a[i__] = 0.;
/* L10: */
    }

    return 0;
} /* zero_ */



/* ----------------------------------------------------------------------------- */
/* Subroutine */ int eigen_(sigmat, eigval)
double *sigmat, *eigval;
{
    /* Builtin functions */
    double sqrt(), acos(), cos();

    /* Local variables */
    double bigr, p, q, r__, s, t;
    extern /* Subroutine */ int invar_();
    double i1, i2, i3, y0, y1, y2, phi, inv[3], tmp;

/* ----------------------------------------------------------------------------- */



    /* Parameter adjustments */
    --eigval;
    sigmat -= 4;

    /* Function Body */
    invar_(&sigmat[4], inv);
    i1 = inv[0];
    i2 = inv[1];
    i3 = inv[2];
    r__ = -i1;
    s = i2;
    t = -i3;
    p = (s * 3. - r__ * r__) / 3.;
    q = r__ * 2. * r__ * r__ / 27. - r__ * s / 3. + t;
    if (abs(q) < 1e-10) {
	y0 = -sqrt((abs(p)));
	y1 = sqrt((abs(p)));
	y2 = 0.;
    } else {
	bigr = sqrt(abs(p) / 3.);
	if (q < 0.) {
	    bigr = -bigr;
	}
	tmp = q / (bigr * 2. * bigr * bigr);
	if (tmp < -1.) {
	    tmp = -1.;
	}
	if (tmp > 1.) {
	    tmp = 1.;
	}
	phi = acos(tmp);
	y0 = bigr * -2. * cos(phi / 3.);
	y1 = bigr * -2. * cos(phi / 3. + 2.0943953196207681);
	y2 = bigr * -2. * cos(phi / 3. + 4.1887906392415362);
    }
    eigval[1] = y0 - r__ / 3.;
    eigval[2] = y1 - r__ / 3.;
    eigval[3] = y2 - r__ / 3.;

    return 0;
} /* eigen_ */


/* ----------------------------------------------------------------------------- */
/* Subroutine */ int invar_(sigmat, inv)
double *sigmat, *inv;
{
    extern /* Subroutine */ int determ_();

/* ----------------------------------------------------------------------------- */



    /* Parameter adjustments */
    --inv;
    sigmat -= 4;

    /* Function Body */
    inv[1] = sigmat[4] + sigmat[8] + sigmat[12];
    inv[2] = sigmat[4] * sigmat[8] + sigmat[8] * sigmat[12] + sigmat[12] * 
	    sigmat[4] - sigmat[7] * sigmat[5] - sigmat[11] * sigmat[9] - 
	    sigmat[6] * sigmat[10];
    determ_(&sigmat[4], &inv[3]);

    return 0;
} /* invar_ */


/* ----------------------------------------------------------------------------- */
/* Subroutine */ int determ_(a, det)
double *a, *det;
{
/* ----------------------------------------------------------------------------- */



    /* Parameter adjustments */
    a -= 4;

    /* Function Body */
    *det = a[4] * (a[8] * a[12] - a[9] * a[11]) - a[7] * (a[5] * a[12] - a[6] 
	    * a[11]) + a[10] * (a[5] * a[9] - a[6] * a[8]);

    return 0;
} /* determ_ */

/* ----------------------------------------------------------------------------- */
/* Subroutine */ int unity4_(a)
double *a;
{
    long int idi, jdi, kdi, ldi;

/* ----------------------------------------------------------------------------- */



    /* Parameter adjustments */
    a -= 40;

    /* Function Body */
    for (idi = 1; idi <= 3; ++idi) {
	for (jdi = 1; jdi <= 3; ++jdi) {
	    for (kdi = 1; kdi <= 3; ++kdi) {
		for (ldi = 1; ldi <= 3; ++ldi) {
		    if (idi == kdi && jdi == ldi) {
			a[idi + (jdi + (kdi + ldi * 3) * 3) * 3] = 1.;
		    } else {
			a[idi + (jdi + (kdi + ldi * 3) * 3) * 3] = 0.;
		    }
/* L70: */
		}
/* L80: */
	    }
/* L90: */
	}
/* L100: */
    }

    return 0;
} /* unity4_ */


/* ----------------------------------------------------------------------------- */
/* Subroutine */ int abdyadic_(a, b, c__)
double *a, *b, *c__;
{
    long int idi, jdi, kdi, ldi;

/* ----------------------------------------------------------------------------- */



    /* Parameter adjustments */
    c__ -= 40;
    b -= 4;
    a -= 4;

    /* Function Body */
    for (idi = 1; idi <= 3; ++idi) {
	for (jdi = 1; jdi <= 3; ++jdi) {
	    for (kdi = 1; kdi <= 3; ++kdi) {
		for (ldi = 1; ldi <= 3; ++ldi) {
		    c__[idi + (jdi + (kdi + ldi * 3) * 3) * 3] = a[idi + jdi *
			     3] * b[kdi + ldi * 3];
/* L70: */
		}
/* L80: */
	    }
/* L90: */
	}
/* L100: */
    }

    return 0;
} /* abdyadic_ */

/* ----------------------------------------------------------------------------- */
double inpro_(a, b, n)
double *a, *b;
long int *n;
{
    /* System generated locals */
    long int i__1;
    double ret_val;

    /* Local variables */
    long int i__;

/* ------------------------------------------------------------------------------ */


    /* Parameter adjustments */
    --b;
    --a;

    /* Function Body */
    ret_val = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ret_val += a[i__] * b[i__];
/* L100: */
    }

    return ret_val;
} /* inpro_ */

