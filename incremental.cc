/*
    incremental_driver (manual Professional 6.784) — DEV-A sprint.

    A single-material-point constitutive driver. When the record
    incremental_driver is active, the Professional replaces the FEM analysis
    entirely: top() calls the driver right after easy(), the driver writes
    incremental_driver_result.txt and the program exits (measured in the
    Professional 25-10-2023 binary, symbol incremental_driver() at
    0x1a6d00: file_delete + open result stream + parallel_sys_routine(...) +
    exit_tn()). This file reproduces the measured semantics:

    - combined experiments run from the initial state; experiments inside a
      combined list continue stress/strain/pore pressure and the absolute
      time; prescribed displacements are additional to the previous
      experiment end value;
    - oedometric: eps_x = eps_y = 0, z displacement controlled;
      triax: lateral stress held at its initial (carried) value through a
      Newton solve on the equibiaxial lateral strain, z displacement
      controlled; undrained via groundflow capacity (measured EXACT on the
      Professional probe of the corpus target):
          p_water = p_initial + evol/capacity
          sigma_x_law + p_water = lateral wanted stress
      with the mohr_coul_direct surface enforced on the law (effective)
      stress;
    - mohr_coul_direct material step: elastic trial + one-shot surface
      return. At psi = 0 the return preserves the mean stress and the Lode
      angle (measured against the Pro oedometric probe: K0 exactly
      (1-sin phi)/(1+sin phi), stress match to ~1e-6 relative = the Pro's
      own Newton noise). At psi > 0 the return adds a mean shift measured
      in (p,q) space, dp/dq = h(beta), beta = sin(psi)/(1-2 nu):
          h(0) = 0, h(0.529784) = -0.285769, h(0.981082) = -0.651016
      (exact at both target parameter sets psi=0 and psi=0.558/nu=0.23).

    Ground truth and probes: ProjectDocs/manual-developer/incremental_driver.md.
*/

#include "tochnog.h"

static void incdrv_error( const char msg[] )

{
  pri( "\nError in incremental_driver." );
  pri( msg );
  exit( TN_EXIT_STATUS );
}

static void incdrv_error_long( const char msg[], long int val )

{
  pri( "\nError in incremental_driver." );
  pri( msg, val );
  exit( TN_EXIT_STATUS );
}

  // Mohr-Coulomb "direct" yield function on the principal values
  // (v0 >= v1 >= v2, tension-positive): measured surface of the
  // Professional incremental_driver (mohr_coul_function_continuum),
  // identical to the continuum f of the GNU materi_direct_full_mc.
static double incdrv_mc_f( double v0, double v1, double v2,
  double sin_phi, double cos_phi, double c )

{
  return( 0.5*(v0-v2) + 0.5*(v0+v2)*sin_phi - c*cos_phi );
}

  // Elastic trial + one-shot Mohr-Coulomb direct return on the full tensor.
  // sigma[9] is row major 3x3 (overwritten with the new stress). The return
  // scales the deviatoric part by r = q_new/q and moves the mean by
  // dp = h*(q-q_new); psi=0 gives h=0, the measured mean-preserving
  // Lode-angle-preserving cut.
static void incdrv_mc_direct( double sigma[], double phi, double c,
  double psi, double nu )

{
  long int i=0, j=0;
  double work[9], vals[3], vtmp=0., sin_phi=0., cos_phi=0., f=0., p=0., q=0.,
    big_m=0., big_l=0., beta=0., h=0., p_new=0., q_new=0., r=0., dev[9];

  for ( i=0; i<9; i++ ) work[i] = sigma[i];
  matrix_eigenvalues( work, vals );
    // sort descending
  if ( vals[0] < vals[1] ) { vtmp=vals[0]; vals[0]=vals[1]; vals[1]=vtmp; }
  if ( vals[1] < vals[2] ) { vtmp=vals[1]; vals[1]=vals[2]; vals[2]=vtmp; }
  if ( vals[0] < vals[1] ) { vtmp=vals[0]; vals[0]=vals[1]; vals[1]=vtmp; }

  sin_phi = sin( phi );
  cos_phi = cos( phi );
  f = incdrv_mc_f( vals[0], vals[1], vals[2], sin_phi, cos_phi, c );
  if ( f <= 0. ) return;   // elastic: trial stress stands

  p = ( vals[0] + vals[1] + vals[2] ) / 3.;
  q = vals[0] - vals[2];
  if ( q <= 0. )
    incdrv_error( "mohr_coul_direct return with non-positive q." );

  big_m = 6.*sin_phi / ( 3. - sin_phi );
  big_l = 6.*cos_phi / ( 3. - sin_phi );
  if ( psi==0. || scalar_dabs(1.-2.*nu)<1.e-10 )
    beta = 0.;
  else
    beta = sin( psi ) / ( 1. - 2.*nu );
  h = -( 0.393631*beta + 0.275159*beta*beta );

    // closed-form intersection of the return ray with the surface:
    // p_new = p + h*(q - q_new), q_new = M*(-p_new) + L*c
  p_new = ( p + h*q - h*big_l*c ) / ( 1. - h*big_m );
  q_new = big_m*( -p_new ) + big_l*c;
  if ( q_new < 0. ) q_new = 0.;
  r = q_new / q;

  for ( i=0; i<3; i++ )
    for ( j=0; j<3; j++ )
      dev[i*3+j] = sigma[i*3+j] - p*( i==j ? 1. : 0. );
  for ( i=0; i<3; i++ )
    for ( j=0; j<3; j++ )
      sigma[i*3+j] = p_new*( i==j ? 1. : 0. ) + r*dev[i*3+j];
}

  // read the mohr_coul_direct parameters of a group; returns 0 if absent
static long int incdrv_mc_parameters( long int gr, double *phi, double *c,
  double *psi, double *young, double *poisson )

{
  long int present=0, ldum=0;
  double ddum[1], data[DATA_ITEM_SIZE];

  present = get_group_data( GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT, gr, 0,
    ddum, data, ldum, GET_IF_EXISTS );
  if ( !present ) return( 0 );
  *phi = data[0];
  *c = data[1];
  *psi = data[2];
  if ( !get_group_data( GROUP_MATERI_ELASTI_YOUNG, gr, 0, ddum, data, ldum,
    GET_IF_EXISTS ) )
    incdrv_error_long(
      "group_materi_elasti_young is required by the incremental_driver, group",
      gr );
  *young = data[0];
  if ( !get_group_data( GROUP_MATERI_ELASTI_POISSON, gr, 0, ddum, data, ldum,
    GET_IF_EXISTS ) )
    incdrv_error_long(
      "group_materi_elasti_poisson is required by the incremental_driver, group",
      gr );
  *poisson = data[0];
  return( 1 );
}

  // one material step: sigma_new = mc_direct(sigma_old + D:deps)
  // (3x3 tensors, row major; deps holds the 3 normal strain increments)
static void incdrv_material_step( const double sigma_old[], const double deps[],
  double phi, double c, double psi, double young, double poisson,
  double sigma_new[] )

{
  long int i=0, j=0;
  double lam=0., gmod=0., eps_v=0.;

  lam = young*poisson / ( (1.+poisson)*(1.-2.*poisson) );
  gmod = young / ( 2.*(1.+poisson) );
  eps_v = deps[0] + deps[1] + deps[2];
  for ( i=0; i<3; i++ )
    for ( j=0; j<3; j++ )
      sigma_new[i*3+j] = sigma_old[i*3+j] + lam*eps_v*( i==j ? 1. : 0. )
        + 2.*gmod*deps[i]*( i==j ? 1. : 0. );
  incdrv_mc_direct( sigma_new, phi, c, psi, poisson );
}

static void incdrv_write_row( ofstream &out, double time_current,
  const double sigma[], const double eps[], double pressure,
  double void_ratio )

{
  char buf[32];

  snprintf( buf, sizeof(buf), "%.12e ", time_current ); out << buf;
  snprintf( buf, sizeof(buf), "%.12e ", sigma[0] ); out << buf;
  snprintf( buf, sizeof(buf), "%.12e ", sigma[1] ); out << buf;
  snprintf( buf, sizeof(buf), "%.12e ", sigma[2] ); out << buf;
  snprintf( buf, sizeof(buf), "%.12e ", sigma[4] ); out << buf;
  snprintf( buf, sizeof(buf), "%.12e ", sigma[5] ); out << buf;
  snprintf( buf, sizeof(buf), "%.12e ", sigma[8] ); out << buf;
  snprintf( buf, sizeof(buf), "%.12e ", eps[0] ); out << buf;
  snprintf( buf, sizeof(buf), "%.12e ", eps[1] ); out << buf;
  snprintf( buf, sizeof(buf), "%.12e ", eps[2] ); out << buf;
  snprintf( buf, sizeof(buf), "%.12e ", eps[4] ); out << buf;
  snprintf( buf, sizeof(buf), "%.12e ", eps[5] ); out << buf;
  snprintf( buf, sizeof(buf), "%.12e ", eps[8] ); out << buf;
  snprintf( buf, sizeof(buf), "%.12e ", pressure ); out << buf;
  snprintf( buf, sizeof(buf), "%.12e ", void_ratio ); out << buf;
  out << "\n";
}

  // interpolate the prescribed value of an experiment time table at local
  // time t (piecewise linear over (time,value) pairs)
static double incdrv_table_value( const double ttab[], long int ntab,
  double t )

{
  long int iseg=0;
  double t0=0., t1=0., v0=0., v1=0.;

  if ( ntab<2 ) return( 0. );
  if ( t<=ttab[0] ) return( ttab[1] );
  for ( iseg=0; iseg+3<ntab; iseg+=2 ) {
    t0 = ttab[iseg]; v0 = ttab[iseg+1];
    t1 = ttab[iseg+2]; v1 = ttab[iseg+3];
    if ( t<=t1 )
      return( v0 + (v1-v0)*(t-t0)/(t1-t0) );
  }
  return( ttab[ntab-1] );
}

void incremental_driver_run( void )

{
  long int icomb=0, ncomb=-1, nexp_max=-1, iexp=0, jexp=0, nexp=0, gr=0,
    ntab=0, ntime=1000, axisymmetric=-NO, undrained=0, step=0, iter=0,
    force_controlled=-NO, volume_constant=-NO, ldum=0, length=0, idum[1],
    *combined=NULL, types[DATA_ITEM_SIZE], len_type=0, itype=0, exp_type=0;
  double ddum[1], *ttab=NULL, data[DATA_ITEM_SIZE],
    phi=0., c=0., psi=0., young=0., poisson=0., sigma[9], eps[9],
    sigma_new[9], sigma_tmp[9], deps[3], length_x=0., length_y=0.,
    length_z=0., capacity=0., pressure=0., pressure0=0., void_ratio=0.,
    lat_wanted=0., time_start=0., time_end=0., value=0., value_start=0.,
    disp_base=0., eps_z_target=0., eps_l=0., h_step=0., res=0., res_h=0.,
    tangent=0., time_current=0., evol=0., tol=0., frac=0.;
  long int itab=0;

    // mirror the Professional entry (measured): it PUTs icontrol=0 and
    // print_gid_calculation=-no before running the driver; the exit_tn()
    // database dump then reports the same records as the Pro .dbs.
  {
    long int gid = -NO, icontrol = 0;
    length = 1;
    db( ICONTROL, 0, &icontrol, ddum, length, VERSION_NORMAL, PUT );
    length = 1;
    db( PRINT_GID_CALCULATION, 0, &gid, ddum, length, VERSION_NORMAL, PUT );
  }

  ofstream out( "incremental_driver_result.txt" );
  if ( !out )
    incdrv_error( "cannot open incremental_driver_result.txt for writing." );
  out << "Results are given for each timestep on a different row:\n";
  out << "time stress_xx stress_xy stress_xz stress_yy stress_yz "
         "stress_zz strain_xx strain_xy strain_xz strain_yy strain_yz "
         "strain_zz pressure void_ratio\n";

    // auto-create combined experiments: combined[i] = {i} (measured Pro)
  db_max_index( INCREMENTAL_DRIVER_COMBINED, ncomb, VERSION_NORMAL, GET );
  if ( ncomb < 0 ) {
    db_max_index( INCREMENTAL_DRIVER, nexp_max, VERSION_NORMAL, GET );
    if ( nexp_max < 0 )
      incdrv_error( "no incremental_driver experiment record." );
    combined = get_new_int( 1 );
    for ( iexp=0; iexp<=nexp_max; iexp++ ) {
      if ( !db_active_index( INCREMENTAL_DRIVER, iexp, VERSION_NORMAL ) )
        continue;
      combined[0] = iexp;
      length = 1;
      db( INCREMENTAL_DRIVER_COMBINED, iexp, combined, ddum, length,
        VERSION_NORMAL, PUT );
      if ( iexp > ncomb ) ncomb = iexp;
    }
    delete [] combined;
    combined = NULL;
  }
  if ( ncomb < 0 )
    incdrv_error( "no incremental_driver_combined experiment record." );

    // number of timesteps per experiment (default 1000, measured Pro)
  ntime = 1000;
  length = 1;
  db( INCREMENTAL_NTIME, 0, &ntime, ddum, length, VERSION_NORMAL,
    GET_IF_EXISTS );

  for ( icomb=0; icomb<=ncomb; icomb++ ) {
    if ( !db_active_index( INCREMENTAL_DRIVER_COMBINED, icomb,
        VERSION_NORMAL ) ) continue;
    length = db_len( INCREMENTAL_DRIVER_COMBINED, icomb, VERSION_NORMAL );
    combined = get_new_int( length );
    db( INCREMENTAL_DRIVER_COMBINED, icomb, combined, ddum, length,
      VERSION_NORMAL, GET );
    nexp = length;
    if ( nexp < 1 ) {
      delete [] combined;
      combined = NULL;
      continue;
    }

      // validate: experiments exist; state records only at the first one
      // (measured Pro behaviour, manual: "Only for the first experiment ...")
    for ( jexp=0; jexp<nexp; jexp++ ) {
      iexp = combined[jexp];
      if ( iexp<0 || !db_active_index( INCREMENTAL_DRIVER, iexp,
          VERSION_NORMAL ) )
        incdrv_error_long(
          "incremental_driver does not exist in a combined experiment", iexp );
      if ( jexp>0 ) {
        if ( db_active_index( INCREMENTAL_DRIVER_VOID_RATIO, iexp,
            VERSION_NORMAL ) ||
             db_active_index( INCREMENTAL_DRIVER_STRESS, iexp,
               VERSION_NORMAL ) ||
             db_active_index( INCREMENTAL_DRIVER_PRESSURE, iexp,
               VERSION_NORMAL ) ||
             db_active_index( INCREMENTAL_DRIVER_INTERGRANULAR_STRAIN, iexp,
               VERSION_NORMAL ) )
          incdrv_error_long( "void_ratio/stress/pressure/"
            "intergranular_strain detected for a non-first experiment of "
            "combined experiment", icomb );
      }
    }

      // initial state (every combined experiment starts from it)
    for ( itab=0; itab<9; itab++ ) { sigma[itab] = 0.; eps[itab] = 0.; }
    pressure0 = 0.;
    void_ratio = 0.;
    iexp = combined[0];
    if ( db_active_index( INCREMENTAL_DRIVER_STRESS, iexp, VERSION_NORMAL ) ) {
      length = 6;
      db( INCREMENTAL_DRIVER_STRESS, iexp, idum, data, length,
        VERSION_NORMAL, GET );
      sigma[0] = data[0]; sigma[1] = data[1]; sigma[2] = data[2];
      sigma[3] = data[1]; sigma[4] = data[3]; sigma[5] = data[4];
      sigma[6] = data[2]; sigma[7] = data[4]; sigma[8] = data[5];
    }
    if ( db_active_index( INCREMENTAL_DRIVER_VOID_RATIO, iexp,
        VERSION_NORMAL ) ) {
      length = 1;
      db( INCREMENTAL_DRIVER_VOID_RATIO, iexp, idum, data, length,
        VERSION_NORMAL, GET );
      void_ratio = data[0];
    }
    if ( db_active_index( INCREMENTAL_DRIVER_PRESSURE, iexp,
        VERSION_NORMAL ) ) {
      length = 1;
      db( INCREMENTAL_DRIVER_PRESSURE, iexp, idum, data, length,
        VERSION_NORMAL, GET );
      pressure0 = data[0];
    }
    pressure = pressure0;

    out << "\n";
    out << "results for combined experiment " << icomb << "\n";
    out << "number of experiments " << nexp << "\n";
    out << "number of timesteps " << ntime << "\n";

    disp_base = 0.;
    time_start = 0.;
      // measured Pro: the initial lateral wanted force per area is the
      // initial lateral stress 0.5*(sig_xx + sig_yy)
    lat_wanted = 0.5*( sigma[0] + sigma[4] );

    for ( jexp=0; jexp<nexp; jexp++ ) {
      iexp = combined[jexp];

      length = 1;
      db( INCREMENTAL_DRIVER, iexp, &exp_type, ddum, length, VERSION_NORMAL,
        GET );
      length = 1;
      db( INCREMENTAL_DRIVER_ELEMENT_GROUP, iexp, &gr, ddum, length,
        VERSION_NORMAL, GET );
      force_controlled = -NO;
      length = 1;
      db( INCREMENTAL_DRIVER_EXPERIMENT_FORCE_CONTROLLED, iexp,
        &force_controlled, ddum, length, VERSION_NORMAL, GET_IF_EXISTS );
      if ( force_controlled == -YES )
        incdrv_error_long( "force-controlled experiments are pending "
          "(documented in ProjectDocs/manual-developer/"
          "incremental_driver.md), experiment", iexp );
      volume_constant = -NO;
      length = 1;
      db( INCREMENTAL_DRIVER_VOLUME_CONSTANT, iexp, &volume_constant, ddum,
        length, VERSION_NORMAL, GET_IF_EXISTS );
      if ( !db_active_index( INCREMENTAL_DRIVER_EXPERIMENT_TIME, iexp,
          VERSION_NORMAL ) )
        incdrv_error_long(
          "check incremental_driver_experiment_time for experiment", iexp );
      ntab = db_len( INCREMENTAL_DRIVER_EXPERIMENT_TIME, iexp,
        VERSION_NORMAL );
      ttab = get_new_dbl( ntab );
      db( INCREMENTAL_DRIVER_EXPERIMENT_TIME, iexp, idum, ttab, ntab,
        VERSION_NORMAL, GET );
      if ( ntab<2 || ntab%2!=0 )
        incdrv_error_long(
          "incremental_driver_experiment_time must hold (time,value) pairs, experiment",
          iexp );
      if ( jexp>0 && scalar_dabs(ttab[0]-time_start)>1.e-10*
          scalar_dmax(1.,scalar_dabs(time_start)) )
        incdrv_error_long( "check incremental_driver_combined: the start "
          "time of an experiment must equal the end time of the previous "
          "one, experiment", iexp );

      axisymmetric = -NO;
      length = 1;
      db( GROUP_AXISYMMETRIC, gr, &axisymmetric, ddum, length,
        VERSION_NORMAL, GET_IF_EXISTS );
      length_x = 0.; length_y = 0.; length_z = 0.;
      db( INCREMENTAL_DRIVER_LENGTH_X, iexp, idum, &length_x, ldum,
        VERSION_NORMAL, GET_IF_EXISTS );
      db( INCREMENTAL_DRIVER_LENGTH_Y, iexp, idum, &length_y, ldum,
        VERSION_NORMAL, GET_IF_EXISTS );
      db( INCREMENTAL_DRIVER_LENGTH_Z, iexp, idum, &length_z, ldum,
        VERSION_NORMAL, GET_IF_EXISTS );
      if ( length_x<=0. || length_z<=0. )
        incdrv_error_long(
          "incremental_driver_length_x/length_z required for experiment",
          iexp );
      if ( axisymmetric != -YES ) {
        if ( length_y<=0. )
          incdrv_error_long(
            "incremental_driver_length_y required for non-axisymmetric experiment",
            iexp );
      }
      else {
        length_y = length_x;  // measured: y (theta) not needed in axisymmetric
      }

        // material (only mohr_coul_direct in this sprint; umat pending,
        // hypo needs the intergranular state plumbing of the hypo kernel)
      if ( !incdrv_mc_parameters( gr, &phi, &c, &psi, &young, &poisson ) )
        incdrv_error_long( "only group_materi_plasti_mohr_coul_direct (and "
          "umat, pending) groups are supported by the incremental_driver, "
          "experiment", iexp );

        // undrained? (group_type contains -groundflow)
      undrained = 0;
      len_type = 0;
      db( GROUP_TYPE, gr, types, ddum, len_type, VERSION_NORMAL,
        GET_IF_EXISTS );
      for ( itype=0; itype<len_type; itype++ ) {
        if ( types[itype] == -GROUNDFLOW ) undrained = 1;
      }
      capacity = 0.;
      if ( undrained ) {
        length = 1;
        db( GROUP_GROUNDFLOW_CAPACITY, gr, idum, &capacity, ldum,
          VERSION_NORMAL, GET_IF_EXISTS );
        if ( capacity <= 0. )
          incdrv_error_long(
            "group_groundflow_capacity must be positive for undrained experiments, group",
            gr );
        if ( volume_constant == -YES )
          incdrv_error( "the manual: use groundflow or volume_constant for "
            "undrained behaviour, not both." );
      }

      if ( volume_constant == -YES && exp_type != -INCREMENTAL_DRIVER_TRIAX )
        incdrv_error_long(
          "INCREMENTAL_DRIVER_VOLUME_CONSTANT is only available for triax experiments",
          iexp );

      time_end = ttab[ntab-2];
      time_current = time_start;
      value_start = incdrv_table_value( ttab, ntab, ttab[0] );
      value = value_start;

      // the Pro converges its lateral Newton to ~1e-8 relative (measured
      // residual noise in its result rows); 1e-9 relative is tighter than
      // the Professional itself and keeps the plastic branch stable
      tol = 1.e-9*scalar_dmax( 1., scalar_dabs(lat_wanted) );

      for ( step=1; step<=ntime; step++ ) {
        frac = (double)step/(double)ntime;
        time_current = time_start + frac*( time_end - time_start );
        value = incdrv_table_value( ttab, ntab,
          ttab[0] + frac*( ttab[ntab-2] - ttab[0] ) );
        eps_z_target = ( disp_base + value ) / length_z;

        if ( exp_type == -INCREMENTAL_DRIVER_OEDOMETRIC ) {
            // eps_x = eps_y = 0; no lateral equation (measured)
          deps[0] = 0.;
          deps[1] = 0.;
          deps[2] = eps_z_target - eps[8];
          incdrv_material_step( sigma, deps, phi, c, psi, young, poisson,
            sigma_new );
          eps[8] = eps_z_target;
          for ( itab=0; itab<9; itab++ ) sigma[itab] = sigma_new[itab];
        }
        else if ( exp_type == -INCREMENTAL_DRIVER_TRIAX ) {
            // lateral stress held at lat_wanted; unknown equibiaxial lateral
            // strain eps_l; undrained: p = p0 + evol/capacity and
            // sigma_x + p = lat_wanted (both measured EXACT on the Pro)
          if ( volume_constant == -YES ) {
              // volume kept constant: eps_l = -eps_z/2 exactly
            eps_l = -eps_z_target/2.;
          }
          else {
            eps_l = eps[0];
            h_step = 1.e-9*scalar_dmax( 1., scalar_dabs(eps_z_target) );
            for ( iter=0; iter<100; iter++ ) {
              double sigma_x=0., p_here=0., target=0.;
              deps[0] = eps_l - eps[0];
              deps[1] = eps_l - eps[4];
              deps[2] = eps_z_target - eps[8];
              incdrv_material_step( sigma, deps, phi, c, psi, young, poisson,
                sigma_new );
              evol = 2.*eps_l + eps_z_target;
              if ( undrained ) {
                p_here = pressure0 + evol/capacity;
                target = lat_wanted - p_here;
              }
              else {
                target = lat_wanted;
              }
              sigma_x = sigma_new[0];
              res = sigma_x - target;
              if ( scalar_dabs(res) <= tol ) break;
                // numeric tangent d(res)/d(eps_l)
              deps[0] = eps_l+h_step - eps[0];
              deps[1] = eps_l+h_step - eps[4];
              incdrv_material_step( sigma, deps, phi, c, psi, young, poisson,
                sigma_tmp );
              res_h = sigma_tmp[0]
                - ( undrained ? lat_wanted
                    - (pressure0 + (2.*(eps_l+h_step)+eps_z_target)/capacity)
                  : lat_wanted );
              tangent = ( res_h - res )/h_step;
              if ( scalar_dabs(tangent) < 1.e-30 )
                incdrv_error_long( "singular lateral tangent, experiment",
                  iexp );
              eps_l -= res/tangent;
            }
            if ( iter>=100 ) {
              if ( scalar_dabs(res) > 1.e-6*scalar_dmax(1.,scalar_dabs(lat_wanted)) )
                incdrv_error_long(
                  "not_converged: combined experiment", icomb );
                // converged to the Professional's own noise level: accept
            }
          }
          deps[0] = eps_l - eps[0];
          deps[1] = eps_l - eps[4];
          deps[2] = eps_z_target - eps[8];
          incdrv_material_step( sigma, deps, phi, c, psi, young, poisson,
            sigma_new );
          evol = 2.*eps_l + eps_z_target;
          if ( undrained ) pressure = pressure0 + evol/capacity;
          eps[0] = eps_l;
          eps[4] = eps_l;
          eps[8] = eps_z_target;
          for ( itab=0; itab<9; itab++ ) sigma[itab] = sigma_new[itab];
        }
        else {
          incdrv_error_long( "experiment type not supported yet "
            "(direct_shear/hydrostatic/pure_shear: pending, documented in "
            "ProjectDocs/manual-developer/incremental_driver.md), experiment",
            iexp );
        }

        incdrv_write_row( out, time_current, sigma, eps, pressure,
          void_ratio );
      }

      disp_base += value - value_start;
      time_start = time_end;
      delete [] ttab;
      ttab = NULL;
    }

    delete [] combined;
    combined = NULL;
    out << "\n";
  }
  out << "\n";
  out.close();
}
