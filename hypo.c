#include <stdio.h>
#include <math.h>
/* hypo.f -- translated by f2c (version 19980831).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static doublereal c_b45 = 1.5;
static integer c__81 = 81;
static integer c__5 = 5;
static integer c__3 = 3;

/* ------------------------------------------------------------------------------ */
/* Subroutine */ int hypo_(stress, mmat, his, inc_ept__, time, dtime, nhis, 
	data, ndata, cohesion, epi_r__, epi_mr__, epi_mt__, epi_betar__, 
	epi_chi__, old_epi__, new_epi__, iusepres, iuseepi, ihypotype, 
	softvar_nonloc, softvar_loc, find_local_sv, options_nonlocal)

integer *find_local_sv, *options_nonlocal;
doublereal *softvar_nonloc, *softvar_loc;
doublereal *stress, *mmat, *his, *inc_ept__, *time, *dtime;
integer *nhis;
doublereal *data;
integer *ndata;
doublereal *cohesion, *epi_r__, *epi_mr__, *epi_mt__, *epi_betar__, *
	epi_chi__, *old_epi__, *new_epi__;
integer *iusepres, *iuseepi, *ihypotype;
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    integer s_wsle(), do_lio(), e_wsle();
    /* Subroutine */ int s_stop();
    double exp();

    /* Local variables */
    static doublereal epid;
    static integer idim, jdim, kdim, ldim;
    static doublereal stress_work__[9]	/* was [3][3] */, rate_epi__[9]	/* 
	    was [3][3] */;
    extern /* Subroutine */ int copy_();
    static doublereal dstr[9]	/* was [3][3] */, his_save__[100], epi_size__,
	     new_epi_save__[9]	/* was [3][3] */, d__[9]	/* was [3][3] 
	    */, previous_epi__[9]	/* was [3][3] */;
    static integer i__;
    static doublereal trace;
    extern /* Subroutine */ int sigma_();
    extern doublereal inpro_();
    static logical error;
    extern /* Subroutine */ int minus_();
    extern doublereal power_();
    static doublereal direction_epi__[9]	/* was [3][3] */, ed, ei, ep, 
	    dt, middle_epi__[9]	/* was [3][3] */;
    extern /* Subroutine */ int pridbl_(), privec_();
    static doublereal dt_tot__;
    static logical switch__;
    extern /* Subroutine */ int mul_();
    static doublereal tmp, epi_rho__;
    extern doublereal normvec_();
    extern /* Subroutine */ int extract_();
    static integer iter_strain__;
    static doublereal stress_save__[9]	/* was [3][3] */;

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, 0, 0 };
    static cilist io___2 = { 0, 6, 0, 0, 0 };

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
	s_wsle(&io___1);
	do_lio(&c__9, &c__1, "Illegal porosity detected in hypoplasticity.", (
		ftnlen)44);
	e_wsle();
	s_wsle(&io___2);
	do_lio(&c__9, &c__1, "Remember to initialise node_dof records.", (
		ftnlen)40);
	e_wsle();
	s_stop("", (ftnlen)0);
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
     	      if (normvec_(stress_work__, &c__9) > normvec_(&stress[4], &c__9) *
	         .01 && normvec_(&stress[4], &c__9) > .1 && dt >= *dtime *
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
	pridbl_("dt    ", &dt, (ftnlen)6);
	privec_("stress", &stress[4], &c__9, (ftnlen)6);
    }

} /* hypo_ */


/* ------------------------------------------------------------------------------ */
/* Subroutine */ int sigma_(switch__, error, stress, mmat, his, dstr, dt, 
	nhis, data, ndata, cohesion, epi_rho__, epi_mr__, epi_mt__, 
	epi_betar__, epi_chi__, epid, direction_epi__, iuseepi, ihypotype,
	softvar_nonloc, softvar_loc, find_local_sv, options_nonlocal, dnorm_nonloc)

integer *find_local_sv, *options_nonlocal;
doublereal *softvar_nonloc, *softvar_loc, *dnorm_nonloc;
logical *switch__, *error;
doublereal *stress, *mmat, *his, *dstr, *dt;
integer *nhis;
doublereal *data;
integer *ndata;
doublereal *cohesion, *epi_rho__, *epi_mr__, *epi_mt__, *epi_betar__, *
	epi_chi__, *epid, *direction_epi__;
integer *iuseepi, *ihypotype;
{


    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4, dmob_1, dmob_2;

    /* Builtin functions */
    double sqrt(), sin(), exp(), asin();

    /* Local variables */
    static doublereal alfa, beta;
    extern /* Subroutine */ int abdyadic_();
    static doublereal effective_stiffness__, eold, phic, edot, lmat[81]	/* 
	    was [3][3][3][3] */, nmat[9]	/* was [3][3] */, nval, rval, that[
	    9]	/* was [3][3] */, tmin, tmax;
    extern /* Subroutine */ int copy_(), zero_();
    static doublereal a, work1[9]	/* was [3][3] */, d__[9]	/* 
	    was [3][3] */, e, f, work2[9]	/* was [3][3] */;
    static integer i__, j, k, l, m, n;
    extern /* Subroutine */ int eigen_();
    extern doublereal power_();
    static doublereal lwork[81]	/* was [3][3][3][3] */, tcohesion[9]	/* 
	    was [3][3] */;
    extern /* Subroutine */ int addfac_();
    static doublereal cos3th, ec, ed, fd, fe, fb, ei, hs, c1, c2, xi, re, powxi;
    extern /* Subroutine */ int unity4_();
    static doublereal eigval[3];
    extern /* Subroutine */ int pridbl_();
    static doublereal d_size__, phimob, ec0, ed0, epi_rhochi__;
    extern /* Subroutine */ int privec_();
    static doublereal tanpsi, ei0;
    extern logical null_array__();
    static doublereal inc_stress__[9]	/* was [3][3] */;
    extern /* Subroutine */ int abc_(), add_();
    static doublereal fac;
    extern /* Subroutine */ int dev_();
    extern doublereal tra_();
    extern /* Subroutine */ int mul_();
    static doublereal tmp;
    extern /* Subroutine */ int a4bc_();
    extern logical negativ_();
    static doublereal thatdev[9]	/* was [3][3] */;
    extern doublereal normvec_();
    extern /* Subroutine */ int extract_();
    extern doublereal normmat_();

/* ------------------------------------------------------------------------------ */





/*       old void ratio */
    /* Parameter adjustments */
    stress -= 4;
    mmat -= 40;
    dstr -= 4;
    --his;
    --data;
    direction_epi__ -= 4;

    /* Function Body */
    eold = his[1];
/*       new void ratio */
    extract_(&dstr[4], dt, d__);
    edot = (eold + 1.) * tra_(d__);
    e = eold + edot * *dt;
    his[1] = e;

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
	pridbl_("tanpsi", &tanpsi, (ftnlen)6);
	pridbl_("cos3th", &cos3th, (ftnlen)6);
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
	pridbl_("f     ", &f, (ftnlen)6);
	pridbl_("a     ", &a, (ftnlen)6);
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
	pridbl_("ei    ", &ei, (ftnlen)6);
	pridbl_("ec    ", &ec, (ftnlen)6);
	pridbl_("ed    ", &ed, (ftnlen)6);
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
        dmob_1 = abs(eigval[0]), dmob_2 = abs(eigval[1]), dmob_1 = max(dmob_1,dmob_2), dmob_2 
	    = abs(eigval[2]);
        tmax = -max(dmob_1,dmob_2);
        dmob_1 = abs(eigval[0]), dmob_2 = abs(eigval[1]), dmob_1 = min(dmob_1,dmob_2), dmob_2 
	    = abs(eigval[2]);
        tmin = -min(dmob_1,dmob_2);
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
    
    if (*switch__) {
	pridbl_("fb    ", &fb, (ftnlen)6);
	pridbl_("fe    ", &fe, (ftnlen)6);
	pridbl_("fd    ", &fd, (ftnlen)6);
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
    
/*        nonlinear */
    tmp = fac * fd * a * f;
    add_(that, thatdev, nmat, &c__9);
    mul_(nmat, &tmp, &c__9);

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
	privec_("stress", &stress[4], &c__9, (ftnlen)6);
    }

/*       mobilised friction angle for postprocessing */
    eigen_(&stress[4], eigval);
/* Computing MAX */
    d__1 = abs(eigval[0]), d__2 = abs(eigval[1]), d__1 = max(d__1,d__2), d__2 
	    = abs(eigval[2]);
    tmax = -max(d__1,d__2);
/* Computing MIN */
    d__1 = abs(eigval[0]), d__2 = abs(eigval[1]), d__1 = min(d__1,d__2), d__2 
	    = abs(eigval[2]);
    tmin = -min(d__1,d__2);
    if ((d__1 = tmax + tmin, abs(d__1)) < 1e-12) {
	phimob = phic;
    } else {
	phimob = asin((d__1 = tmax - tmin, abs(d__1)) / (d__2 = tmax + tmin, 
		abs(d__2)));
	phimob = phimob * 360. / 6.2831853999999998;
    }
    his[3] = phimob;

/*       effective stiffness */
    effective_stiffness__ = normvec_(&mmat[40], &c__81);
    /*his[4] = effective_stiffness__;*/
    his[4] = re;

L1000:

    ;
} /* sigma_ */


/* ------------------------------------------------------------------------------ */
/* Subroutine */ int abc_(a, b, c__)
doublereal *a, *b, *c__;
{
    static integer idi, jdi, kdi;

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

} /* abc_ */


/* ------------------------------------------------------------------------------ */
/* Subroutine */ int a4bc_(a, b, c__)
doublereal *a, *b, *c__;
{
    static integer idi, jdi, kdi, ldi;

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

} /* a4bc_ */


/* ----------------------------------------------------------------------------- */
/* Subroutine */ int add_(a, b, c__, n)
doublereal *a, *b, *c__;
integer *n;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

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

} /* add_ */


/* ----------------------------------------------------------------------------- */
/* Subroutine */ int addfac_(tot, inc, factor, n)
doublereal *tot, *inc, *factor;
integer *n;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

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

} /* addfac_ */


/* ----------------------------------------------------------------------------- */
/* Subroutine */ int copy_(a, b, n)
doublereal *a, *b;
integer *n;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

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

} /* copy_ */


/* ----------------------------------------------------------------------------- */
/* Subroutine */ int dev_(a, adev)
doublereal *a, *adev;
{
    static integer idi, jdi;
    extern doublereal tra_();
    static doublereal tmp;

/* ------------------------------------------------------------------------------ */



    /* Parameter adjustments */
    adev -= 4;
    a -= 4;

    /* Function Body */
    tmp = tra_(&a[4]);
    for (idi = 1; idi <= 3; ++idi) {
	for (jdi = 1; jdi <= 3; ++jdi) {
	    adev[idi + jdi * 3] = a[idi + jdi * 3];
/* L10: */
	}
	adev[idi + idi * 3] -= tra_(&a[4]) / 3.;
/* L20: */
    }

} /* dev_ */


/* ----------------------------------------------------------------------------- */
/* Subroutine */ int extract_(dstr, dt, d__)
doublereal *dstr, *dt, *d__;
{
    /* System generated locals */
    doublereal d__1;

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

} /* extract_ */


/* ----------------------------------------------------------------------------- */
/* Subroutine */ int minus_(a, b, c__, n)
doublereal *a, *b, *c__;
integer *n;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

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

} /* minus_ */


/* ----------------------------------------------------------------------------- */
/* Subroutine */ int mul_(a, scal, n)
doublereal *a, *scal;
integer *n;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/* ----------------------------------------------------------------------------- */


    /* Parameter adjustments */
    --a;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a[i__] = *scal * a[i__];
/* L10: */
    }

} /* mul_ */


/* ----------------------------------------------------------------------------- */
logical negativ_(x)
doublereal *x;
{
    /* System generated locals */
    logical ret_val;

/* ------------------------------------------------------------------------------ */


    if (*x < (double)1e-8) {
	ret_val = TRUE_;
    } else {
	ret_val = FALSE_;
    }

    return ret_val;
} /* negativ_ */


/* ------------------------------------------------------------------------------ */
doublereal normmat_(a)
doublereal *a;
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal work[9]	/* was [3][3] */;
    extern /* Subroutine */ int abc_();
    extern doublereal tra_();

/* ------------------------------------------------------------------------------ */



    /* Parameter adjustments */
    a -= 4;

    /* Function Body */
    abc_(&a[4], &a[4], work);
    ret_val = sqrt(tra_(work));

    return ret_val;
} /* normmat_ */


/* ------------------------------------------------------------------------------ */
doublereal normvec_(a, n)
doublereal *a;
integer *n;
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static integer i__;

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
logical null_array__(x)
doublereal *x;
{
    /* System generated locals */
    logical ret_val;

/* ------------------------------------------------------------------------------ */


    if (abs(*x) < (double)1e-8) {
	ret_val = TRUE_;
    } else {
	ret_val = FALSE_;
    }

    return ret_val;
} /* null_array__ */


/* ------------------------------------------------------------------------------ */
doublereal power_(x, y)
doublereal *x, *y;
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double pow_dd();

/* ------------------------------------------------------------------------------ */


    ret_val = pow_dd(x, y);

    return ret_val;
} /* power_ */


/* ------------------------------------------------------------------------------ */
/* Subroutine */ int pridbl_(label, a, label_len)
char *label;
doublereal *a;
ftnlen label_len;
{
    /* Builtin functions */
    integer s_wsle(), do_lio(), e_wsle();

    /* Fortran I/O blocks */
    static cilist io___94 = { 0, 6, 0, 0, 0 };


/* ------------------------------------------------------------------------------ */


    s_wsle(&io___94);
    do_lio(&c__9, &c__1, label, (ftnlen)6);
    do_lio(&c__5, &c__1, (char *)&(*a), (ftnlen)sizeof(doublereal));
    e_wsle();

} /* pridbl_ */


/* Subroutine */ int priint_(label, a, label_len)
char *label;
integer *a;
ftnlen label_len;
{
    /* Builtin functions */
    integer s_wsle(), do_lio(), e_wsle();

    /* Fortran I/O blocks */
    static cilist io___95 = { 0, 6, 0, 0, 0 };


/* ------------------------------------------------------------------------------ */


    s_wsle(&io___95);
    do_lio(&c__9, &c__1, label, (ftnlen)6);
    do_lio(&c__3, &c__1, (char *)&(*a), (ftnlen)sizeof(integer));
    e_wsle();

} /* priint_ */


/* ------------------------------------------------------------------------------ */
/* Subroutine */ int pritxt_(label, label_len)
char *label;
ftnlen label_len;
{
    /* Builtin functions */
    integer s_wsle(), do_lio(), e_wsle();

    /* Fortran I/O blocks */
    static cilist io___96 = { 0, 6, 0, 0, 0 };


/* ------------------------------------------------------------------------------ */


    s_wsle(&io___96);
    do_lio(&c__9, &c__1, label, (ftnlen)6);
    e_wsle();

} /* pritxt_ */


/* ------------------------------------------------------------------------------ */
/* Subroutine */ int privec_(label, a, n, label_len)
char *label;
doublereal *a;
integer *n;
ftnlen label_len;
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsle(), do_lio(), e_wsle();

    /* Local variables */
    static integer i__;

    /* Fortran I/O blocks */
    static cilist io___97 = { 0, 6, 0, 0, 0 };
    static cilist io___99 = { 0, 6, 0, 0, 0 };


/* ------------------------------------------------------------------------------ */


    /* Parameter adjustments */
    --a;

    /* Function Body */
    s_wsle(&io___97);
    do_lio(&c__9, &c__1, label, (ftnlen)6);
    e_wsle();
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_wsle(&io___99);
	do_lio(&c__5, &c__1, (char *)&a[i__], (ftnlen)sizeof(doublereal));
	e_wsle();
/* L10: */
    }

} /* privec_ */


/* ------------------------------------------------------------------------------ */
doublereal tra_(a)
doublereal *a;
{
    /* System generated locals */
    doublereal ret_val;

/* ------------------------------------------------------------------------------ */



    /* Parameter adjustments */
    a -= 4;

    /* Function Body */
    ret_val = a[4] + a[8] + a[12];

    return ret_val;
} /* tra_ */


/* ----------------------------------------------------------------------------- */
/* Subroutine */ int zero_(a, n)
doublereal *a;
integer *n;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/* ----------------------------------------------------------------------------- */


    /* Parameter adjustments */
    --a;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a[i__] = 0.;
/* L10: */
    }

} /* zero_ */



/* ----------------------------------------------------------------------------- */
/* Subroutine */ int eigen_(sigmat, eigval)
doublereal *sigmat, *eigval;
{
    /* Builtin functions */
    double sqrt(), acos(), cos();

    /* Local variables */
    static doublereal bigr, p, q, r__, s, t;
    extern /* Subroutine */ int invar_();
    static doublereal i1, i2, i3, y0, y1, y2, phi, inv[3], tmp;

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

} /* eigen_ */


/* ----------------------------------------------------------------------------- */
/* Subroutine */ int invar_(sigmat, inv)
doublereal *sigmat, *inv;
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

} /* invar_ */


/* ----------------------------------------------------------------------------- */
/* Subroutine */ int determ_(a, det)
doublereal *a, *det;
{
/* ----------------------------------------------------------------------------- */



    /* Parameter adjustments */
    a -= 4;

    /* Function Body */
    *det = a[4] * (a[8] * a[12] - a[9] * a[11]) - a[7] * (a[5] * a[12] - a[6] 
	    * a[11]) + a[10] * (a[5] * a[9] - a[6] * a[8]);

} /* determ_ */

/* ----------------------------------------------------------------------------- */
/* Subroutine */ int unity4_(a)
doublereal *a;
{
    static integer idi, jdi, kdi, ldi;

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

} /* unity4_ */


/* ----------------------------------------------------------------------------- */
/* Subroutine */ int abdyadic_(a, b, c__)
doublereal *a, *b, *c__;
{
    static integer idi, jdi, kdi, ldi;

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

} /* abdyadic_ */

/* ----------------------------------------------------------------------------- */
doublereal inpro_(a, b, n)
doublereal *a, *b;
integer *n;
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer i__;

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

