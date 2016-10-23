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
    max_bounda=0, max_bounda_unknown=0, max_bounda_force=0,
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
    time_start=0., radius=0., angle_start=0.,
    angle_total=0., rdum=0., ddum[MDIM], coord_start[MDIM], coord_total[MDIM],
    *bounda_time=NULL, *new_node_dof=NULL, 
    *node_dof=NULL, *bounda_sine=NULL, *node_rhside=NULL;

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
  db_max_index( BOUNDA_FORCE, max_bounda_force, VERSION_NORMAL, GET );
  if ( max_bounda_unknown<0 && max_bounda_force<0 ) goto end_of_bounda;

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

  if ( max_bounda_unknown>max_bounda ) max_bounda = max_bounda_unknown;
  if ( max_bounda_force>max_bounda ) max_bounda = max_bounda_force;
  for ( iboun=0; iboun<=max_bounda; iboun++ ) {
    unknown = db_active_index( BOUNDA_UNKNOWN, iboun, VERSION_NORMAL );
    force   = db_active_index( BOUNDA_FORCE, iboun, VERSION_NORMAL );
    bounda_time_user = -NO; db( BOUNDA_TIME_USER, iboun, &bounda_time_user, ddum, 
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
    if ( unknown || force ) {
      if ( swit ) pri( "iboun", iboun );
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
        if ( length_bounda_time==1 )
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

      if ( unknown ) {
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
                geometry( in, ddum, val, found, factor, ddum, rdum,
                  ddum, NODE_START_REFINED, PROJECT_EXACT, VERSION_NORMAL );
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
                      new_node_dof[iuknwn] = factor * load;
                      if ( derivatives ) new_node_dof[ind1] = 
                        ( new_node_dof[iuknwn] - node_dof[iuknwn] ) / dtime;
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
