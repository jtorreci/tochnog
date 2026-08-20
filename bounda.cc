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

#define EPS_ATAN 1.e-10

void bounda( )

{
  long int in=0, ready=0, iuknwn=0, inod=0, iboun=0, found=0,
    max_bounda=0, max_bounda_unknown=0, max_bounda_dof=0, max_bounda_force=0,
    bounda_time_user=0, ind1=0, ipuknwn=0, iu=0, iu_start=0, iu_end=0,
    inc=0, ninc=0, length=0, range_length=0, unknown=0, force=0,
    idim=0, time=0, sine=0, user=0, rotate=0, rotate_axis=0,
    use_geom=0, swit=0, swit_node=0, bounda_length=0,
    use_range=0, use_all=0, use_node_set=0, ifreq=0, nfreq=0, indx=0, max_node=0,
    bounda_time_file=0, length_bounda_time=0, ldum=0, idum[1], 
    *val=NULL, *integer_range=NULL, 
    *dof_label=NULL, *dof_type=NULL, *node_bounded=NULL;
  double load=0., time0=0., time1=0., load0=0., load1=0., factor=0.,
    time_current=0., time_total=0., dtime=0., amplitude=0., frequency=0.,
    time_start=0., radius=0., angle_start=0., on_time=0., off_time=0.,
    until_force=0., until_factor=1., reaction=0.,
    angle_total=0., rdum=0., ddum[MDIM], coord_start[MDIM], coord_total[MDIM],
    *bounda_time=NULL, *new_node_dof=NULL, 
    *node_dof=NULL, *bounda_sine=NULL, *node_rhside=NULL;
  long int bounda_on_off=0, bounda_until_force=0, bounda_constant=0,
    bounda_geometry_method=0, bounda_alternate_list[DATA_ITEM_SIZE],
    bounda_alternate_n=0, iteration=0, bounda_water=0;
  double bounda_normal_vec[3], bounda_dof_radial[3],
    bounda_dof_cylindrical[6];
  double bounda_time_increment=0., bounda_time_offset=0.;
  double bounda_factor[4], bounda_factor_px[3], bounda_time_units[2];

  swit = set_swit(-1,-1,"bounda");
  if ( swit ) pri( "In routine BOUNDA" );

  val = get_new_int(MBOUNDA);
  integer_range = get_new_int(MRANGE);
  dof_label = get_new_int(MUKNWN);
  dof_type = get_new_int(MUKNWN);
  bounda_time = get_new_dbl(DATA_ITEM_SIZE);

  db_set_int( NODE_BOUNDED, VERSION_NORMAL );

  groundflow_phreatic_apply();

  db_max_index( BOUNDA_UNKNOWN, max_bounda_unknown, VERSION_NORMAL, GET );
  db_max_index( BOUNDA_DOF, max_bounda_dof, VERSION_NORMAL, GET );
  db_max_index( BOUNDA_FORCE, max_bounda_force, VERSION_NORMAL, GET );
  if ( max_bounda_unknown<0 && max_bounda_dof<0 && max_bounda_force<0 )
    goto end_of_bounda;

  db_max_index( NODE, max_node, VERSION_NORMAL, GET );
  if ( max_node<0 ) goto end_of_bounda;

  nodes_in_geometry = get_new_int( 1+max_node );

  db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( DOF_TYPE, 0, dof_type, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET_IF_EXISTS );
  db( TIME_CURRENT, 0, idum, &time_current, ldum, 
    VERSION_NORMAL, GET_IF_EXISTS );
  time_total = time_current + dtime;
  if ( swit ) pri( "time_total", time_total );
  db( NUMBER_ITERATIONS, 0, &iteration, ddum, ldum, VERSION_NEW, GET_IF_EXISTS );

  if ( max_bounda_unknown>max_bounda ) max_bounda = max_bounda_unknown;
  if ( max_bounda_dof>max_bounda ) max_bounda = max_bounda_dof;
  if ( max_bounda_force>max_bounda ) max_bounda = max_bounda_force;
  for ( iboun=0; iboun<=max_bounda; iboun++ ) {
    unknown = db_active_index( BOUNDA_UNKNOWN, iboun, VERSION_NORMAL ) ||
              db_active_index( BOUNDA_DOF, iboun, VERSION_NORMAL );
    force   = db_active_index( BOUNDA_FORCE, iboun, VERSION_NORMAL );
    bounda_time_user = -NO; db( BOUNDA_TIME_USER, iboun, &bounda_time_user, ddum, 
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
    bounda_on_off = 0;
    if ( unknown || force ) {
      if ( swit ) pri( "iboun", iboun );
      // bounda_alternate: in successive iterations, omit one of the
      // listed bounda_dof indices (rotating). Useful for very large runs.
      bounda_alternate_n = 0;
      db( BOUNDA_ALTERNATE, 0, bounda_alternate_list, ddum,
        bounda_alternate_n, VERSION_NORMAL, GET_IF_EXISTS );
      if ( bounda_alternate_n>0 ) {
        long int jalt=0, skip_bounda=0;
        for ( jalt=0; jalt<bounda_alternate_n; jalt++ ) {
          if ( bounda_alternate_list[jalt]==iboun &&
               (iteration % bounda_alternate_n)==jalt ) {
            skip_bounda = 1;
            break;
          }
        }
        if ( skip_bounda ) continue;
      }
      time =  sine = user = 0; ninc = 2;
      if      ( db_active_index( BOUNDA_SINE, iboun, VERSION_NORMAL ) ) {
        bounda_sine = db_dbl( BOUNDA_SINE, iboun, VERSION_NORMAL );
        time_start = bounda_sine[0];
        nfreq = ( db_len( BOUNDA_SINE, iboun, VERSION_NORMAL ) - 1 ) / 2;
        ninc = 2;
        sine = 1;
      }
      else if ( bounda_time_user==-YES ) {
        ninc = 2;
        user = 1;
      }          
      else if ( db_active_index( BOUNDA_TIME_FILE, iboun, VERSION_NORMAL ) ) {
        db( BOUNDA_TIME_FILE, iboun, &bounda_time_file, ddum, ldum, 
          VERSION_NORMAL, GET_IF_EXISTS );
        if ( bounda_time_file==-YES ) {
          bounda_time_file_apply( iboun, time_total, bounda_time, ninc );
          time = 1;
        }
        else 
          db_error( BOUNDA_TIME_FILE, iboun );
      }
      else if ( db_active_index( BOUNDA_TIME, iboun, VERSION_NORMAL ) ) {
        db( BOUNDA_TIME, iboun, idum, bounda_time, length_bounda_time, 
          VERSION_NORMAL, GET );
        time = 1;
        // bounda_time_units: convert time and length units in bounda_time
        bounda_time_units[0] = 1.; bounda_time_units[1] = 1.;
        db( BOUNDA_TIME_UNITS, iboun, idum, bounda_time_units, ldum,
          VERSION_NORMAL, GET_IF_EXISTS );
        if ( bounda_time_units[0]!=1. || bounda_time_units[1]!=1. ) {
          long int iu=0;
          for ( iu=0; iu<length_bounda_time; iu++ ) {
            if ( iu%2==0 ) bounda_time[iu] *= bounda_time_units[0];
            else bounda_time[iu] *= bounda_time_units[1];
          }
        }
        db( BOUNDA_TIME_INCREMENT, iboun, idum, &bounda_time_increment, ldum,
          VERSION_NORMAL, GET_IF_EXISTS );
        if ( bounda_time_increment>0. )
          ninc = length_bounda_time;
        else if ( length_bounda_time==1 )
          ninc = 2;
        else {
          if ( length_bounda_time<4 ) db_error( BOUNDA_TIME, iboun );
          ninc = length_bounda_time / 2;
        }
      }
      else {
        ninc = 2;
        time = 1;
        length_bounda_time = 0;
      }
      // only periodically use the bounda_time values
      if ( db_active_index( BOUNDA_TIME_ON_OFF, iboun, VERSION_NORMAL ) ) {
        double on_off_tmp[2];
        db( BOUNDA_TIME_ON_OFF, iboun, idum, on_off_tmp, ldum, 
          VERSION_NORMAL, GET );
        on_time = on_off_tmp[0];
        off_time = on_off_tmp[1];
        if ( on_time<=0. || off_time<0. ) db_error( BOUNDA_TIME_ON_OFF, iboun );
        bounda_on_off = 1;
      }
      // limit the force response by reducing the prescribed velocity
      if ( db_active_index( BOUNDA_TIME_UNTIL_FORCE, iboun, VERSION_NORMAL ) ) {
        double until_tmp[2];
        db( BOUNDA_TIME_UNTIL_FORCE, iboun, idum, until_tmp, ldum, 
          VERSION_NORMAL, GET );
        until_force = until_tmp[0];
        until_factor = until_tmp[1];
        if ( until_factor<=0. || until_factor>1. )
          db_error( BOUNDA_TIME_UNTIL_FORCE, iboun );
        bounda_until_force = 1;
      }
      // keep the prescribed dofs constant (bounda_constant)
      db( BOUNDA_CONSTANT, iboun, &bounda_constant, ddum, ldum,
        VERSION_NORMAL, GET_IF_EXISTS );
      db( BOUNDA_TIME_OFFSET, iboun, idum, &bounda_time_offset, ldum,
        VERSION_NORMAL, GET_IF_EXISTS );
      array_set( bounda_factor, 0., 4 );
      db( BOUNDA_FACTOR, iboun, idum, bounda_factor, ldum,
        VERSION_NORMAL, GET_IF_EXISTS );
      array_set( bounda_factor_px, 0., 3 );
      db( BOUNDA_FACTOR_PARABOLIC_X, iboun, idum, bounda_factor_px, ldum,
        VERSION_NORMAL, GET_IF_EXISTS );
      db( BOUNDA_GEOMETRY_METHOD, iboun, &bounda_geometry_method, ddum, ldum,
        VERSION_NORMAL, GET_IF_EXISTS );
      array_set( bounda_normal_vec, 0., 3 );
      db( BOUNDA_NORMAL, iboun, idum, bounda_normal_vec, ldum,
        VERSION_NORMAL, GET_IF_EXISTS );
      array_set( bounda_dof_radial, 0., 3 );
      db( BOUNDA_DOF_RADIAL, iboun, idum, bounda_dof_radial, ldum,
        VERSION_NORMAL, GET_IF_EXISTS );
      array_set( bounda_dof_cylindrical, 0., 6 );
      db( BOUNDA_DOF_CYLINDRICAL, iboun, idum, bounda_dof_cylindrical, ldum,
        VERSION_NORMAL, GET_IF_EXISTS );
      bounda_water = 0;
      db( BOUNDA_WATER, iboun, &bounda_water, ddum, ldum,
        VERSION_NORMAL, GET_IF_EXISTS );

      if ( unknown ) {
        if ( db_active_index( BOUNDA_DOF, iboun, VERSION_NORMAL ) )
          db( BOUNDA_DOF, iboun, val, ddum, bounda_length, 
            VERSION_NORMAL, GET );
        else
          db( BOUNDA_UNKNOWN, iboun, val, ddum, bounda_length, 
            VERSION_NORMAL, GET );
        if ( bounda_length<2 ) db_error( BOUNDA_UNKNOWN, iboun );
        rotate = 0; rotate_axis = val[bounda_length-1];
        if      ( rotate_axis==-ROTATION_X_AXIS ) {
          rotate = 1;
          val[bounda_length-1] = dof_label[vel_indx+1*nder];
          val[bounda_length-0] = dof_label[vel_indx+2*nder];
          bounda_length++;
        }
        else if ( rotate_axis==-ROTATION_Y_AXIS ) {
          rotate = 1;
          val[bounda_length-1] = dof_label[vel_indx+0*nder];
          val[bounda_length-0] = dof_label[vel_indx+2*nder];
          bounda_length++;
        }
        else if ( rotate_axis==-ROTATION_Z_AXIS ) {
          rotate = 1;
          val[bounda_length-1] = dof_label[vel_indx+0*nder];
          val[bounda_length-0] = dof_label[vel_indx+1*nder];
          bounda_length++;
        }
      }
      else {
        assert( force );
        db( BOUNDA_FORCE, iboun, val, ddum, bounda_length, 
          VERSION_NORMAL, GET );
        if ( bounda_length<2 ) db_error( BOUNDA_FORCE, iboun );
      }
      if ( val[0]<0 && db_data_class(val[0])==GEOMETRY ) {
        array_move( val, geometry_ent, 2 );
        array_set( nodes_in_geometry, 0, 1+max_node );
        parallel_sys_routine( &parallel_geometry );
      }

      found = 0;
      for ( inc=0; !found && inc<ninc-1; inc++ ) {

        if ( time ) {
          if      ( length_bounda_time==0 ) {
            load = 0.;
            found = 1;
          }
          else if ( length_bounda_time==1 ) {
            load =  bounda_time[0];
            found = 1;
          }
          else {
            if ( bounda_time_increment>0. ) {
              // bounda_time holds load-only values; times are offset + k*increment
              long int nload = length_bounda_time;
              if ( time_total>=bounda_time_offset ) {
                long int k = (long int) floor( (time_total-bounda_time_offset)
                  / bounda_time_increment + 1.e-9 );
                if ( k<0 ) k = 0;
                if ( k>nload-1 ) k = nload-1;
                load = bounda_time[k];
                found = 1;
              }
            }
            else {
              time0 = bounda_time[inc*2+0];
              load0 = bounda_time[inc*2+1];
              time1 = bounda_time[inc*2+2];
              load1 = bounda_time[inc*2+3];
              if ( time0>=time1 ) db_error( BOUNDA_TIME, iboun );
              if ( swit ) {
                pri( "time_total", time_total );
                pri( "time0", time0 );
                pri( "time1", time1 );
                pri( "load0", load0 );
                pri( "load1", load1 );
              }
              if ( time_total>=(time0-1.e-10) && time_total<=time1 ) {
                found = 1;
                if ( time0==time1 ) load = load0;
                else load = load0 + (load1-load0)*(time_total-time0)/(time1-time0);
              }
            }
            // only periodically use the bounda_time values
            if ( found && bounda_on_off ) {
              double phase = fmod( time_total, on_time+off_time );
              if ( phase>=on_time ) { load = 0.; found = 0; }
            }
          }
        }
        else if ( user ) {
          load = 0.;
          found = 1;
          user_bounda_time( iboun, time_total, load );
        }
        else if ( sine ) 
          found = time_current>=time_start;
        else {
          load = 0.;
          found = 1;
        }

        if ( found ) {

          iu_end = bounda_length - 1;
          use_geom = use_range = use_all = use_node_set = 0;
          if      ( val[0]==-RA ) {
            use_range = 1;
            range_expand( val, integer_range, length, range_length );
            iu_start = length;
          }
          else if ( val[0]==-ALL ) {
            use_all  = 1;
            iu_start = 1;
          }
          else if ( val[0]==-NODE_SET ) {
            use_node_set = 1;
            iu_start = 1;
          }
          else if ( val[0]<0 && db_data_class(val[0])==GEOMETRY ) {
            use_geom = 1;
            iu_start = 2;
          }
          else
            iu_start = 1;

          for ( in=0, ready=0; !ready; in++ ) {
            found = 0;
            if      ( use_range ) {
              inod = integer_range[in];
              factor = 1.;
              found = 1;
            }
            else if ( use_all ) {
              if ( db_active_index( NODE, in, VERSION_NORMAL ) ) {
                inod = in;
                factor = 1.;
                found = 1;
              }
            }
            else if ( use_geom  ) {
              if ( nodes_in_geometry[in] ) {
                found = 1;
                inod = in;
                long int node_type = NODE_START_REFINED;
                if ( bounda_geometry_method!=0 )
                  node_type = bounda_geometry_method;
                geometry( in, ddum, val, found, factor, ddum, rdum,
                  ddum, node_type, PROJECT_EXACT, VERSION_NORMAL );
              }
            }
            else if ( use_node_set  ) {
              inod = in;
              if ( db_active_index( NODE_SET, inod, VERSION_NORMAL ) ) {
                found = 1;
                factor = 1.;
              }
            }
            else {
              inod = val[0];
              factor = 1.;
              found = 1;
            }
            if ( found ) {
              if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
                swit_node = swit;
                swit = set_swit(-1,inod,"bounda");
                if ( swit ) pri( "inod", inod );
                for ( iu=iu_start; iu<=iu_end; iu++ ) {
                  array_member( dof_label, val[iu], nuknwn, iuknwn );
                  if ( iuknwn<0 ) {
                    if ( unknown ) 
                      db_error( BOUNDA_UNKNOWN, iboun );
                    else 
                      db_error( BOUNDA_FORCE, iboun );
                  }
                  ipuknwn = iuknwn / nder;
                  if ( unknown ) {
                    node_bounded = db_int( NODE_BOUNDED, inod, VERSION_NORMAL );
                    new_node_dof = db_dbl( NODE_DOF, inod, VERSION_NEW );
                    node_dof = db_dbl( NODE_DOF, inod, VERSION_NORMAL );
                    if ( swit ) pri( "node_dof", node_dof, nuknwn );
                    ind1 = ipuknwn*nder + nder - 1;
                    node_bounded[ipuknwn] = 1;
                    // groundflow seepage: on a seepage edge the pore pressure
                    // is prescribed (typically 0, free air) ONLY when water
                    // flows OUT of the domain. When water would flow IN, the
                    // edge is closed (no boundary condition). The normal of the
                    // seepage geometry must point outwards the material.
                    if ( iuknwn==pres_indx && groundflow_pressure &&
                         db_active_index( GROUNDFLOW_SEEPAGE_EPS, 0,
                           VERSION_NORMAL ) ) {
                      long int in_seep=0, inod_tmp=0;
                      double coords_tmp[MDIM], seep_normal[MDIM], flow_dot=0.;
                      double seep_eps=0.1;
                      db( GROUNDFLOW_SEEPAGE_EPS, 0, idum, &seep_eps, ldum,
                        VERSION_NORMAL, GET_IF_EXISTS );
                      // node-based seepage list
                      if ( db_active_index( GROUNDFLOW_SEEPAGE_NODE, 0,
                          VERSION_NORMAL ) ) {
                        long int snode[DATA_ITEM_SIZE], lsn=0;
                        db( GROUNDFLOW_SEEPAGE_NODE, 0, snode, ddum, lsn,
                          VERSION_NORMAL, GET );
                        if ( array_member( snode, inod, lsn, ldum ) ) in_seep = 1;
                      }
                      // geometry-based seepage list
                      long int max_seep=0, iseep=0;
                      db_max_index( GROUNDFLOW_SEEPAGE_GEOMETRY, max_seep,
                        VERSION_NORMAL, GET );
                      for ( iseep=0; iseep<=max_seep && !in_seep; iseep++ ) {
                        if ( db_active_index( GROUNDFLOW_SEEPAGE_GEOMETRY,
                            iseep, VERSION_NORMAL ) ) {
                          long int gent[2];
                          db( GROUNDFLOW_SEEPAGE_GEOMETRY, iseep, gent, ddum,
                            ldum, VERSION_NORMAL, GET );
                          db( NODE, inod, idum, coords_tmp, ldum,
                            VERSION_NORMAL, GET );
                          long int gfound=0;
                          double gfac=0., gnor[MDIM], gpen=0., gproj[MDIM];
                          geometry( inod, coords_tmp, gent, gfound, gfac,
                            gnor, gpen, gproj, NODE_START_REFINED,
                            PROJECT_EXACT, VERSION_NORMAL );
                          if ( gfound ) {
                            in_seep = 1;
                            array_move( gnor, seep_normal, ndim );
                            array_normalize( seep_normal, ndim );
                            // flow = -k*grad(pres) ~ -gvel (Darcy velocity);
                            // check the component along the outward normal.
                            flow_dot = 0.;
                            if ( groundflow_velocity ) {
                              double *gvel = db_dbl( NODE_DOF, inod,
                                VERSION_NORMAL );
                              for ( inod_tmp=0; inod_tmp<ndim; inod_tmp++ )
                                flow_dot += gvel[gvel_indx+inod_tmp*nder] *
                                  seep_normal[inod_tmp];
                            }
                          }
                        }
                      }
                      // water flows OUT when flow_dot>0 (Darcy velocity points
                      // outwards). If it flows IN, close the edge: do not
                      // impose the pressure boundary condition.
                      if ( in_seep && flow_dot < -seep_eps ) {
                        node_bounded[ipuknwn] = 0;
                      }
                    }
                    if ( sine ) {
                      new_node_dof[iuknwn] = new_node_dof[ind1] = 0.;
                      for ( ifreq=0; ifreq<nfreq; ifreq++ ) {
                        frequency = bounda_sine[1+ifreq*2+0];
                        amplitude = bounda_sine[1+ifreq*2+1];
                        new_node_dof[iuknwn] += factor * amplitude *
                          sin( 2. * PIRAD * frequency * time_total );
                        if ( derivatives ) new_node_dof[ind1] += 
                          factor * amplitude * 2. * PIRAD * frequency *
                          cos( 2. * PIRAD * frequency * time_total );
                      }
                    }
                    else if ( rotate ) {
                      db( NODE, inod, idum, coord_start, ldum, 
                        VERSION_NORMAL, GET );
                      if ( materi_displacement ) {
                        for ( idim=0; idim<ndim; idim++ ) 
                          coord_start[idim] += node_dof[dis_indx+idim*nder];
                      }
                      if      ( rotate_axis==-ROTATION_X_AXIS ) {
                        radius = sqrt( coord_start[1]*coord_start[1] +
                          coord_start[2]*coord_start[2] );
                        if ( scalar_dabs(coord_start[1])<EPS_ATAN ) {
                          if ( coord_start[2]>=0. )
                            angle_start = 1.*PIRAD/2.;
                          else
                            angle_start = 3.*PIRAD/2.;
                        }
                        else {
                          angle_start = 
                            atan(scalar_dabs(coord_start[2]/coord_start[1]));
                          if      ( coord_start[1]<0. && coord_start[2]>=0. )
                            angle_start = 1.*PIRAD - angle_start;
                          else if ( coord_start[1]<0. && coord_start[2]<0. )
                            angle_start = 1.*PIRAD + angle_start;
                          else if ( coord_start[1]>0. && coord_start[2]<0. )
                            angle_start = 2.*PIRAD - angle_start;
                        }
                        angle_total = angle_start + 
                          factor * load * dtime * PIRAD / 180.;
                        coord_total[0] = coord_start[0];                        
                        coord_total[1] = radius * cos(angle_total);                        
                        coord_total[2] = radius * sin(angle_total);                        
                      }
                      else if ( rotate_axis==-ROTATION_Y_AXIS ) {
                        radius = sqrt( coord_start[0]*coord_start[0] +
                          coord_start[2]*coord_start[2] );
                        if ( scalar_dabs(coord_start[0])<EPS_ATAN ) {
                          if ( coord_start[2]>=0. )
                            angle_start = 1.*PIRAD/2.;
                          else
                            angle_start = 3.*PIRAD/2.;
                        }
                        else {
                          angle_start = 
                            atan(scalar_dabs(coord_start[2]/coord_start[0]));
                          if      ( coord_start[0]<0. && coord_start[2]>=0. )
                            angle_start = 1.*PIRAD - angle_start;
                          else if ( coord_start[0]<0. && coord_start[2]<0. )
                            angle_start = 1.*PIRAD + angle_start;
                          else if ( coord_start[0]>0. && coord_start[2]<0. )
                            angle_start = 2.*PIRAD - angle_start;
                        }
                        angle_total = angle_start - 
                          factor * load * dtime * PIRAD / 180.;
                        coord_total[1] = coord_start[1];                        
                        coord_total[0] = radius * cos(angle_total);                        
                        coord_total[2] = radius * sin(angle_total);                        
                      }
                      else {
                        assert( rotate_axis==-ROTATION_Z_AXIS );
                        radius = sqrt( coord_start[0]*coord_start[0] +
                          coord_start[1]*coord_start[1] );
                        if ( scalar_dabs(coord_start[0])<EPS_ATAN ) {
                          if ( coord_start[1]>=0. )
                            angle_start = 1.*PIRAD/2.;
                          else
                            angle_start = 3.*PIRAD/2.;
                        }
                        else {
                          angle_start = 
                            atan(scalar_dabs(coord_start[1]/coord_start[0]));
                          if      ( coord_start[0]<0. && coord_start[1]>=0. )
                            angle_start = 1.*PIRAD - angle_start;
                          else if ( coord_start[0]<0. && coord_start[1]<0. )
                            angle_start = 1.*PIRAD + angle_start;
                          else if ( coord_start[0]>0. && coord_start[1]<0. )
                            angle_start = 2.*PIRAD - angle_start;
                        }
                        angle_total = angle_start + 
                          factor * load * dtime * PIRAD / 180.;
                        coord_total[2] = coord_start[2];                        
                        coord_total[0] = radius * cos(angle_total);                        
                        coord_total[1] = radius * sin(angle_total);                        
                      }
                      idim = ( iuknwn - vel_indx ) / nder;
                      new_node_dof[iuknwn] = 
                        ( coord_total[idim] - coord_start[idim] ) / dtime;
                      if ( derivatives ) new_node_dof[ind1] = 
                        ( new_node_dof[iuknwn] - node_dof[iuknwn] ) / dtime;
                    }
                    else {
                      if ( bounda_constant==-YES && iuknwn>=vel_indx &&
                           node_dof!=NULL && !force &&
                           node_dof[iuknwn]!=0. )
                        new_node_dof[iuknwn] = node_dof[iuknwn];
                      else {
                        // coordinate-dependent factor on the load
                        double load_factor = 1.;
                        if ( bounda_factor[0]!=0. || bounda_factor[1]!=0. ||
                             bounda_factor[2]!=0. || bounda_factor[3]!=0. ) {
                          db( NODE, inod, idum, coord_start, ldum,
                            VERSION_NORMAL, GET );
                          load_factor = bounda_factor[0];
                          if ( ndim>=1 ) load_factor += bounda_factor[1]*coord_start[0];
                          if ( ndim>=2 ) load_factor += bounda_factor[2]*coord_start[1];
                          if ( ndim==3 ) load_factor += bounda_factor[3]*coord_start[2];
                        }
                        if ( bounda_factor_px[0]!=0. || bounda_factor_px[1]!=0. ||
                             bounda_factor_px[2]!=0. ) {
                          db( NODE, inod, idum, coord_start, ldum,
                            VERSION_NORMAL, GET );
                          load_factor = bounda_factor_px[0]
                            + bounda_factor_px[1]*coord_start[0]
                            + bounda_factor_px[2]*coord_start[0]*coord_start[0];
                        }
                        if ( bounda_water==-YES && iuknwn==pres_indx ) {
                          // pore pressure from the water column height:
                          // density_water * g * (water_level - y)
                          double coords_w[MDIM], sp=0., wl=0., dens=0.,
                            fg[MDIM];
                          db( NODE, inod, idum, coords_w, ndim,
                            VERSION_NORMAL, GET );
                          force_gravity_calculate( fg );
                          if ( db_active_index( GROUNDFLOW_DENSITY, 0,
                              VERSION_NORMAL ) )
                            dens = db_dbl( GROUNDFLOW_DENSITY, 0,
                              VERSION_NORMAL )[0];
                          if ( db_active_index( GROUNDFLOW_PHREATICLEVEL, 0,
                              VERSION_NORMAL ) ) {
                            long int plen=0, idum2[1];
                            db( GROUNDFLOW_PHREATICLEVEL, 0, idum2, &wl, plen,
                              VERSION_NORMAL, GET );
                            if ( plen>=1 ) {
                              double *gpv = db_dbl(
                                GROUNDFLOW_PHREATICLEVEL, 0, VERSION_NORMAL );
                              wl = gpv[0];
                            }
                          }
                          sp = fg[ndim-1] * dens * ( wl - coords_w[ndim-1] );
                          new_node_dof[iuknwn] = factor * sp;
                        }
                        else
                          new_node_dof[iuknwn] = factor * load * load_factor;
                        // bounda_dof_radial: prescribe velocity radial to a point
                        if ( (bounda_dof_radial[0]!=0. || bounda_dof_radial[1]!=0.
                              || bounda_dof_radial[2]!=0.) &&
                             iuknwn>=vel_indx && iuknwn<vel_indx+ndim*nder ) {
                          long int idim_r = ( iuknwn - vel_indx ) / nder;
                          double coords_r[MDIM], r=0.;
                          db( NODE, inod, idum, coords_r, ndim,
                            VERSION_NORMAL, GET );
                          for ( long int k=0; k<ndim; k++ ) {
                            double dk = coords_r[k] - bounda_dof_radial[k];
                            r += dk*dk;
                          }
                          r = sqrt(r);
                          if ( r>0. && idim_r<ndim ) {
                            double dr = coords_r[idim_r]
                              - bounda_dof_radial[idim_r];
                            new_node_dof[iuknwn] *= dr / r;
                          }
                        }
                        // bounda_dof_cylindrical: prescribe velocity cylindrical
                        // to a line defined by two points (radial to the line)
                        if ( (bounda_dof_cylindrical[0]!=0. ||
                              bounda_dof_cylindrical[1]!=0. ||
                              bounda_dof_cylindrical[2]!=0.) &&
                             iuknwn>=vel_indx && iuknwn<vel_indx+ndim*nder ) {
                          long int idim_c = ( iuknwn - vel_indx ) / nder;
                          double coords_c[MDIM], p1[MDIM], p2[MDIM],
                            axis[MDIM], rc[MDIM], r=0.;
                          for ( long int k=0; k<MDIM; k++ ) {
                            p1[k] = bounda_dof_cylindrical[k];
                            p2[k] = bounda_dof_cylindrical[k+3];
                            axis[k] = p2[k] - p1[k];
                          }
                          db( NODE, inod, idum, coords_c, ndim,
                            VERSION_NORMAL, GET );
                          // project the node onto the line: rc = coord - p1
                          for ( long int k=0; k<ndim; k++ )
                            rc[k] = coords_c[k] - p1[k];
                          // distance from the node to the line
                          double t=0., aa=0.;
                          for ( long int k=0; k<ndim; k++ ) {
                            t   += rc[k]*axis[k];
                            aa  += axis[k]*axis[k];
                          }
                          if ( aa>0. ) {
                            t /= aa;
                            double proj[MDIM];
                            for ( long int k=0; k<ndim; k++ ) {
                              proj[k] = p1[k] + t*axis[k];
                              double dk = coords_c[k] - proj[k];
                              r += dk*dk;
                            }
                            r = sqrt(r);
                            if ( r>0. && idim_c<ndim ) {
                              double dr = coords_c[idim_c]
                                - proj[idim_c];
                              new_node_dof[iuknwn] *= dr / r;
                            }
                          }
                        }
                      }
                      if ( derivatives ) new_node_dof[ind1] = 
                        ( new_node_dof[iuknwn] - node_dof[iuknwn] ) / dtime;
                      // limit the force response by reducing the velocity
                      if ( bounda_until_force && force==0 ) {
                        long int ireac = ( iuknwn - vel_indx ) / nder;
                        if ( db_active_index( NODE_RHSIDE_PREVIOUS, inod,
                            VERSION_NORMAL ) ) {
                          node_rhside = db_dbl( NODE_RHSIDE_PREVIOUS, inod,
                            VERSION_NORMAL );
                          reaction = scalar_dabs( node_rhside[ireac] );
                          if ( reaction>until_force )
                            new_node_dof[iuknwn] *= 
                              until_factor*(until_force/reaction);
                        }
                      }
                    }
                    if ( dof_type[iuknwn]==-MATERI_DISPLACEMENT ) {
                      indx = vel_indx + ( iuknwn - dis_indx );
                      new_node_dof[indx] = 
                      ( new_node_dof[iuknwn] - node_dof[iuknwn] ) / dtime;
                      node_bounded[indx/nder] = 1;
                    }
                    if ( dof_type[iuknwn]==-MATERI_VELOCITY_INTEGRATED ) {
                      indx = vel_indx + ( iuknwn - veli_indx );
                      new_node_dof[indx] = 
                      ( new_node_dof[iuknwn] - node_dof[iuknwn] ) / dtime;
                      node_bounded[indx/nder] = 1;
                    }
                    if ( dof_type[iuknwn]==-MAXWELL_E ) {
                      indx = maxfe_indx + ( iuknwn - maxe_indx );
                      new_node_dof[indx] = 
                      ( new_node_dof[iuknwn] - node_dof[iuknwn] ) / dtime;
                      node_bounded[indx/nder] = 1;
                    }
                    if ( dof_type[iuknwn]==-WAVE_SCALAR ) {
                      new_node_dof[fscal_indx] = 
                      ( new_node_dof[iuknwn] - node_dof[iuknwn] ) / dtime;
                      node_bounded[fscal_indx/nder] = 1;
                    }
                    if ( swit ) {
                      pri( "node_bounded", node_bounded, npuknwn );
                      pri( "new_node_dof", new_node_dof, nuknwn );
                    }
                    // bounda_normal: restrict the node to slide on a plane
                    // (velocity component normal to the plane is set to zero)
                    if ( (bounda_normal_vec[0]!=0. || bounda_normal_vec[1]!=0. ||
                          bounda_normal_vec[2]!=0.) && materi_velocity ) {
                      double nn = bounda_normal_vec[0]*bounda_normal_vec[0]
                        + bounda_normal_vec[1]*bounda_normal_vec[1]
                        + bounda_normal_vec[2]*bounda_normal_vec[2];
                      if ( nn>0. ) {
                        double vn=0.;
                        long int idim_n=0;
                        for ( idim_n=0; idim_n<ndim; idim_n++ ) {
                          long int iv = vel_indx + idim_n*nder;
                          vn += new_node_dof[iv]*bounda_normal_vec[idim_n];
                        }
                        for ( idim_n=0; idim_n<ndim; idim_n++ ) {
                          long int iv = vel_indx + idim_n*nder;
                          new_node_dof[iv] -= vn*bounda_normal_vec[idim_n]/nn;
                        }
                      }
                    }
                  }
                  else if ( force ) {
                    if      ( dof_type[iuknwn]==-MATERI_DISPLACEMENT )
                      indx = vel_indx + ( iuknwn - dis_indx );
                    else if ( dof_type[iuknwn]==-MATERI_VELOCITY_INTEGRATED )
                      indx = vel_indx + ( iuknwn - veli_indx );
                    else if ( dof_type[iuknwn]==-MATERI_DISPLACEMENT )
                      indx = fscal_indx + ( iuknwn - scal_indx );
                    else
                      indx = iuknwn;
                    indx /= nder;
                    node_rhside = db_dbl( NODE_RHSIDE, inod, VERSION_NORMAL );
                    if ( sine ) {
                      node_rhside[indx] = 0.;
                      for ( ifreq=0; ifreq<nfreq; ifreq++ ) { 
                        frequency = bounda_sine[1+ifreq*2+0];
                        amplitude = bounda_sine[1+ifreq*2+1];
                        node_rhside[indx] += factor*amplitude*
                          sin(2.*PIRAD*frequency*time_total);
                      }
                    }
                    else {
                      node_rhside[indx] = factor * load;
                    }
                  }
                }
                swit = swit_node;
              }
            }
            if      ( use_range )
              ready = ( (in+1)==range_length );
            else if ( use_all )
              ready = ( in==max_node );
            else if ( use_geom )
              ready = ( in==max_node );
            else if ( use_node_set )
              ready = ( in==max_node );
            else
              ready = 1;
          }
        }
      }
    }
  }

  delete[] nodes_in_geometry;

  end_of_bounda:

  delete[] val;
  delete[] integer_range;
  delete[] dof_label;
  delete[] dof_type;
  delete[] bounda_time;

  if ( swit ) pri( "Out routine BOUNDA" );

}

void bounda_time_file_apply( long int iboun, double time_total,
  double bounda_time[], long int &ninc )

{
  long int n=0, found=0;
  char str[MCHAR], filename[MCHAR];
  double old_load=0., old_time=0., new_load=0., new_time=0.;

  strcpy( filename, long_to_a(iboun,str) );
  strcat( filename, ".bounda" );
  ifstream in( filename );
  if ( !in ) db_error( BOUNDA_TIME_FILE, iboun );

  ninc = 0;
  for ( ; !in.eof() && !found ; ) {
    n++;
    old_load = new_load;
    old_time = new_time;
    if ( !(in>>new_time) || !(in>>new_load) ) 
      db_error( BOUNDA_TIME_FILE, iboun );
    if ( n>1 && time_total>=old_time && time_total<=new_time ) {
      found = 1;
      ninc = 2;
      bounda_time[0] = old_time;
      bounda_time[1] = old_load;
      bounda_time[2] = new_time;
      bounda_time[3] = new_load;
    }
  }
  if ( ninc<2 ) db_error( BOUNDA_TIME_FILE, iboun );
  
  in.close();
}
