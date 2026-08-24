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
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#include "tochnog.h"

void condif( long int element, long int gr, long int nnol, double h[],
  double coord_ip[], double volume, double new_unknowns[],
  double element_lhside[],
  double element_matrix[], double element_rhside[],
  double element_residue[] )

{

  long int swit=0, inol=0, ipuknwn=0, indx=0, ldum=0;
  double condif_absorption=0.;

  if ( get_group_data( GROUP_CONDIF_ABSORPTION, gr, element, new_unknowns,
      &condif_absorption, ldum, GET_IF_EXISTS ) ) {
    swit = set_swit(element,-1,"condif");
    if ( swit ) pri( "In routine CONDIF" );
      // absorption
    for ( inol=0; inol<nnol; inol++ ) {
      ipuknwn = temp_indx/nder;
      indx = inol*npuknwn + ipuknwn;
      element_rhside[indx] -= volume * h[inol] * condif_absorption *
        new_unknowns[temp_indx];
      if ( residue ) element_residue[indx] -= h[inol] * condif_absorption *
        new_unknowns[temp_indx];
      element_lhside[indx] += volume * h[inol] * condif_absorption;
      element_matrix[indx*nnol*npuknwn+indx] +=
        volume * h[inol] * condif_absorption;
    }
    if ( swit ) {
      pri( "element_lhside", element_lhside, nnol, npuknwn );
      pri( "element_rhside", element_rhside, nnol, npuknwn );
      if ( residue ) pri( "element_residue", element_residue, nnol, npuknwn );
      pri( "Out routine CONDIF" );
    }
  }

  // condif_heat_volume (manual Professional 6.83-6.91): distributed volume
  // heat source S. Contribution to the temperature equation:
  //   rhs_i += volume * h_i * S * load * factor
  // with load from _user/_sine/_time, factor the spatial polynomial of
  // _factor at the integration point, and restrictions by _element,
  // _element_group and _geometry (all element nodes in the geometry).
  {
    long int max_heat=0, iheat=0, length=0, idum[1], length_el=0,
      heat_volume_user=0, nfreq=0, ifreq=0, all_in=0, in_geometry=0,
      inod=0, i=0, el[MNOL+1], geometry_entity[2];
    double ddum[1], heat=0., load=1., factor=1., time_current=0., dtime=0.,
      time_total=0., time_start=0., frequency=0., amplitude=0.,
      values[DATA_ITEM_SIZE], coord[MDIM], rdum=0.;
    long int *heat_element=NULL, *heat_group=NULL;
    double *heat_sine=NULL, *heat_time=NULL;

    db_max_index( CONDIF_HEAT_VOLUME, max_heat, VERSION_NORMAL, GET );
    if ( max_heat>=0 ) {
      swit = set_swit(element,-1,"condif");
      if ( swit ) pri( "In routine CONDIF (heat_volume)" );
      db( TIME_CURRENT, 0, idum, &time_current, ldum,
        VERSION_NORMAL, GET_IF_EXISTS );
      db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET_IF_EXISTS );
      time_total = time_current + dtime;
      for ( iheat=0; iheat<=max_heat; iheat++ ) {
        if ( !db_active_index( CONDIF_HEAT_VOLUME, iheat, VERSION_NORMAL ) )
          continue;

        // restrictions
        if ( db_active_index( CONDIF_HEAT_VOLUME_ELEMENT, iheat,
            VERSION_NORMAL ) ) {
          heat_element = db_int( CONDIF_HEAT_VOLUME_ELEMENT, iheat,
            VERSION_NORMAL );
          length = db_len( CONDIF_HEAT_VOLUME_ELEMENT, iheat,
            VERSION_NORMAL );
          if ( !array_member( heat_element, element, length, ldum ) )
            continue;
        }
        if ( db_active_index( CONDIF_HEAT_VOLUME_ELEMENT_GROUP, iheat,
            VERSION_NORMAL ) ) {
          heat_group = db_int( CONDIF_HEAT_VOLUME_ELEMENT_GROUP, iheat,
            VERSION_NORMAL );
          length = db_len( CONDIF_HEAT_VOLUME_ELEMENT_GROUP, iheat,
            VERSION_NORMAL );
          if ( !array_member( heat_group, gr, length, ldum ) ) continue;
        }
        if ( db_active_index( CONDIF_HEAT_VOLUME_GEOMETRY, iheat,
            VERSION_NORMAL ) ) {
          db( CONDIF_HEAT_VOLUME_GEOMETRY, iheat, geometry_entity, ddum,
            ldum, VERSION_NORMAL, GET );
          db( ELEMENT, element, el, ddum, length_el, VERSION_NORMAL, GET );
          all_in = 1;
          for ( i=1; i<length_el; i++ ) {
            inod = el[i];
            geometry( inod, ddum, geometry_entity, in_geometry, rdum,
              ddum, rdum, ddum, NODE_START_REFINED,
              CONDIF_HEAT_VOLUME_GEOMETRY, VERSION_NORMAL );
            if ( !in_geometry ) all_in = 0;
          }
          if ( !all_in ) continue;
        }

        // heat value: the record, or the user function when _user -yes
        db( CONDIF_HEAT_VOLUME, iheat, idum, values, ldum,
          VERSION_NORMAL, GET );
        heat = values[0];
        heat_volume_user = -NO;
        db( CONDIF_HEAT_VOLUME_USER, iheat, &heat_volume_user, ddum, ldum,
          VERSION_NORMAL, GET_IF_EXISTS );
        if ( heat_volume_user==-YES )
          user_condif_heat_volume( iheat, time_total, coord_ip, heat );

        // temporal load
        if ( db_active_index( CONDIF_HEAT_VOLUME_SINE, iheat,
            VERSION_NORMAL ) ) {
          heat_sine = db_dbl( CONDIF_HEAT_VOLUME_SINE, iheat,
            VERSION_NORMAL );
          nfreq = ( db_len( CONDIF_HEAT_VOLUME_SINE, iheat,
            VERSION_NORMAL ) - 1 ) / 2;
          time_start = heat_sine[0];
          load = 0.;
          if ( time_total>time_start ) {
            for ( ifreq=0; ifreq<nfreq; ifreq++ ) {
              frequency = heat_sine[1+ifreq*2+0];
              amplitude = heat_sine[1+ifreq*2+1];
              load += amplitude * sin( 2. * PIRAD * frequency * time_total );
            }
          }
        }
        else if ( db_active_index( CONDIF_HEAT_VOLUME_TIME, iheat,
            VERSION_NORMAL ) ) {
          heat_time = db_dbl( CONDIF_HEAT_VOLUME_TIME, iheat,
            VERSION_NORMAL );
          length = db_len( CONDIF_HEAT_VOLUME_TIME, iheat,
            VERSION_NORMAL );
          force_time( heat_time, "CONDIF_HEAT_VOLUME_TIME", length, load );
        }
        else
          load = 1.;

        // spatial factor at the integration point
        for ( i=0; i<ndim; i++ ) coord[i] = coord_ip[i];
        force_factor( CONDIF_HEAT_VOLUME_FACTOR, iheat, coord, factor );

        if ( swit ) {
          pri( "heat_volume heat", heat );
          pri( "heat_volume load", load );
          pri( "heat_volume factor", factor );
        }

        for ( inol=0; inol<nnol; inol++ ) {
          ipuknwn = temp_indx/nder;
          indx = inol*npuknwn + ipuknwn;
          element_rhside[indx] +=
            volume * h[inol] * heat * load * factor;
        }
      }
      if ( swit ) pri( "Out routine CONDIF (heat_volume)" );
    }
  }

}
