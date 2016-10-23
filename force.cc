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

void force_element_volume_set( long int element, long int nnol, 
  long int nodes[], double coord[], double force_element_volume[] )

{
  long int length=0, iforce=0, max_force=0, use_geom=0, nel=0,
    inol=0, inod=0, ok=0, ifreq=0, nfreq=0,
    ldum=0, idum[1], *force_element_volume_geometry=NULL;
  double factor=1., forcefac=0., frequency=0., amplitude=0., rdum=0., 
    dtime=0., time_current=0., time_total=0., time_start=0., ddum[MDIM], 
    force_work[DATA_ITEM_SIZE], *force_element_volume_time=NULL,
     *force_element_volume_sine=NULL;

  db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET_IF_EXISTS );
  db( TIME_CURRENT, 0, idum, &time_current, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  time_total = time_current + dtime;

  array_set( force_element_volume, 0., nprinc );
  db_max_index( FORCE_ELEMENT_VOLUME, max_force, VERSION_NORMAL, GET );
  for ( iforce=0; iforce<=max_force; iforce++ ) {
    if ( db_active_index( FORCE_ELEMENT_VOLUME, iforce, VERSION_NORMAL ) ) {
      db( FORCE_ELEMENT_VOLUME, iforce, idum, force_work, ldum, VERSION_NORMAL, GET );
      ok = 1;
      if ( db_active_index( FORCE_ELEMENT_VOLUME_GEOMETRY, iforce, VERSION_NORMAL ) ) {
        force_element_volume_geometry = 
          db_int( FORCE_ELEMENT_VOLUME_GEOMETRY, iforce, VERSION_NORMAL );
        if ( force_element_volume_geometry[0]>0 ) {
          nel = db_len( FORCE_ELEMENT_VOLUME_GEOMETRY, iforce, VERSION_NORMAL );
          ok = array_member( force_element_volume_geometry, element, nel, ldum );
        }
        else {
          use_geom = force_element_volume_geometry[0]<0 &&
            db_data_class(force_element_volume_geometry[0])==GEOMETRY;
          if ( !use_geom ) db_error( FORCE_ELEMENT_VOLUME_GEOMETRY, iforce );
          for ( inol=0; inol<nnol && ok; inol++ ) {
            inod = nodes[inol];
            geometry( inod, ddum, force_element_volume_geometry, ok, rdum, ddum, rdum,
              ddum, NODE_START_REFINED, PROJECT_EXACT, VERSION_NORMAL );
          }
        }
      }
      if ( ok ) {
        force_factor( FORCE_ELEMENT_VOLUME_FACTOR, iforce, coord, factor );
        forcefac = factor;
        if      ( db_active_index( FORCE_ELEMENT_VOLUME_SINE, iforce, VERSION_NORMAL ) ) {
          force_element_volume_sine = db_dbl( FORCE_ELEMENT_VOLUME_SINE, iforce, VERSION_NORMAL );
          nfreq = ( db_len( FORCE_ELEMENT_VOLUME_SINE, iforce, VERSION_NORMAL ) - 1 ) / 2;
          time_start = force_element_volume_sine[0];
          factor = 0.;
          if ( time_total>time_start ) {
            for ( ifreq=0; ifreq<nfreq; ifreq++ ) {
              frequency = force_element_volume_sine[1+ifreq*2+0];
              amplitude = force_element_volume_sine[1+ifreq*2+1];
              factor += amplitude * sin( 2. * PIRAD * frequency * time_total );
            }
          }
          forcefac *= factor;
        }
        else if ( db_active_index( FORCE_ELEMENT_VOLUME_TIME, iforce, VERSION_NORMAL ) ) {
          length = db_len( FORCE_ELEMENT_VOLUME_TIME, iforce, VERSION_NORMAL );
          force_element_volume_time = 
            db_dbl( FORCE_ELEMENT_VOLUME_TIME, iforce, VERSION_NORMAL );
          force_time( force_element_volume_time, "FORCE_ELEMENT_VOLUME_TIME",
            length, factor );
          forcefac *= factor;
        }
        array_multiply( force_work, force_work, forcefac, nprinc );
        array_add( force_work, force_element_volume, force_element_volume, nprinc );
      }
    }
  }

}

long int force_time( double time_table[], const char* table_name, 
  long int length, double &load )

{
  long int found=0, ldum=0, idum[1];
  double dtime=0., time_current=0., time_total=0.;

  db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET_IF_EXISTS );
  db( TIME_CURRENT, 0, idum, &time_current, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  time_total = time_current + dtime;

  found = table_xy( time_table, table_name, length, time_total, load );
  if ( !found ) load = 0.;

  return found;
}

void force_time_file_apply( long int iforce, long int table_name, double &load )

{
  long int n=0, found=0, ninc=0, ldum=0, idum[1];
  char str[MCHAR], filename[MCHAR];
  double old_load=0., old_time=0., new_load=0., new_time=0.;
  double force_time_file_table[4];

  double dtime=0., time_current=0., time_total=0.;
  db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET_IF_EXISTS );
  db( TIME_CURRENT, 0, idum, &time_current, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  time_total = time_current + dtime;

  strcpy( filename, long_to_a(iforce,str) );
  strcat( filename, ".force" );
  ifstream in( filename );
  if ( !in ) db_error( table_name, iforce );

  ninc = 0;
  for ( ; !in.eof() && !found ; ) {
    n++;
    old_load = new_load;
    old_time = new_time;
    if ( !(in>>new_time) || !(in>>new_load) ) 
      db_error( table_name, iforce );
    if ( n>1 && time_total>=old_time && time_total<=new_time ) {
      found = 1;
      ninc = 2;
      force_time_file_table[0] = old_time;
      force_time_file_table[1] = old_load;
      force_time_file_table[2] = new_time;
      force_time_file_table[3] = new_load;
    }
  }
  if ( ninc<2 ) db_error( table_name, iforce );
  
  long int length=2*ninc;
  found = table_xy( force_time_file_table, "FORCE_TIME_FILE", length, time_total, load );
  if ( !found ) load = 0.;

  in.close();
}

void force_factor( long int factor_name, long int iforce, 
  double coord[], double &factor )

{
  long int i=0, j=0, level=0, length=0;
  double x=0., y=0., z=0., *factor_name_values=NULL;

  factor = 1;
  if ( db_active_index( factor_name, iforce, VERSION_NORMAL ) ) {
    length = db_len( factor_name, iforce, VERSION_NORMAL );
    factor_name_values = db_dbl( factor_name, iforce, VERSION_NORMAL );
    if ( ndim==1 ) {
      x = coord[0];
      factor = factor_name_values[0];
      for ( j=level=1; j<length; level++ ) {
        factor += factor_name_values[j] * scalar_power(x,level); j++;
      }
    }
    else if( ndim==2 ) {
      x = coord[0];
      y = coord[1];
      j = 0;
      factor = factor_name_values[j]; j++;
      for ( level=1; j<length; level++ ) {
        if ( j<length ) {
          factor += factor_name_values[j] * scalar_power(x,level); j++;
        }
        for ( i=level-1; i>0 && j<length; i-- ) {
          factor += factor_name_values[j] * 
            scalar_power(x,i) * scalar_power(y,level-i); j++;
        }
        if ( j<length ) {
          factor += factor_name_values[j] * scalar_power(y,level); j++;
        }
      }
    }
    else {
      assert( ndim==3 );
      z = coord[2];
      factor = factor_name_values[0];
      for ( j=0; j<length; j++ ) {
        factor += factor_name_values[j] * scalar_power(z,j); 
      }
    }
  }
}

void force_gravity_calculate( double force_gravity[] ) 

{
  long int idim=0, length=0, icontrol=0, options_skip_gravity=-NO, 
    ldum=0, idum[1];
  double factor=0., *force_gravity_time=NULL, ddum[1];

  array_set( force_gravity, 0., ndim );                    
  if ( db_active_index( FORCE_GRAVITY, 0, VERSION_NORMAL ) ) {
    db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
    db( OPTIONS_SKIP_GRAVITY, 0, &options_skip_gravity, 
      ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( CONTROL_OPTIONS_SKIP_GRAVITY, icontrol, &options_skip_gravity, 
      ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    if ( options_skip_gravity==-NO ) {
      db( FORCE_GRAVITY, 0, idum, force_gravity, ldum, VERSION_NORMAL, GET );
      if ( db_active_index( FORCE_GRAVITY_TIME, 0, VERSION_NORMAL ) ) {
        length = db_len( FORCE_GRAVITY_TIME, 0, VERSION_NORMAL );
        force_gravity_time = db_dbl( FORCE_GRAVITY_TIME, 0, VERSION_NORMAL );
        force_time( force_gravity_time, "FORCE_GRAVITY_TIME", length, factor );
        for ( idim=0; idim<ndim; idim++ ) force_gravity[idim] *= factor;
      }
    }
  }

}
