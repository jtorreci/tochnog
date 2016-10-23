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
 
    On June 30, 2013 made change to correct proper calculation of velocity 
    as the case of velocity calculated from permeability x (dh/dx) 
    had wrong sign. F. Lorenzo.
*/

#include "tochnog.h"

void groundflow( long int element, long int gr, long int nnol,
  double coord_ip[], double h[], double d[], 
  double volume, double old_unknowns[], 
  double new_unknowns[], double grad_new_unknowns[], 
  double element_matrix[], double element_rhside[],
  double element_residue[] )

{
  long int swit=0, inol=0, jnol=0, jdim=0, ipuknwn=0, iuknwn=0, jpuknwn=0,
    indx=0, indxi=0, indxj=0, icontrol=0, options_skip_groundflow_materidivergence=-NO,
    materidivergence=-YES, ldum=0, idum[1];
  double tmp=0., C=0., dtime=0., divergence=0., dens=0., ddum[1], pe[MDIM];

  swit = set_swit(element,-1,"groundflow");
  if ( swit ) pri( "In routine GROUNDFLOW." );

  if ( !db_active_index( GROUP_GROUNDFLOW_PERMEABILITY, gr, VERSION_NORMAL ) &&
       !db_active_index( GROUP_GROUNDFLOW_PERMEABILITY_NONLINEAR_METHOD, gr, VERSION_NORMAL ) &&
       !db_active_index( GROUP_GROUNDFLOW_CAPACITY, gr, VERSION_NORMAL ) &&
       !db_active_index( GROUP_GROUNDFLOW_CAPACITY_NONLINEAR_METHOD, gr, VERSION_NORMAL ) )
    return;
         
  db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET );
  db( GROUNDFLOW_DENSITY, 0, idum, &dens, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( GROUP_GROUNDFLOW_MATERIDIVERGENCE, gr, &materidivergence, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  db( OPTIONS_SKIP_GROUNDFLOW_MATERIDIVERGENCE, 0, &options_skip_groundflow_materidivergence, 
    ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_OPTIONS_SKIP_GROUNDFLOW_MATERIDIVERGENCE, icontrol, 
    &options_skip_groundflow_materidivergence, ddum, ldum, 
    VERSION_NORMAL, GET_IF_EXISTS );
  if ( options_skip_groundflow_materidivergence==-YES ) materidivergence = -NO;

  groundflow_data( element, gr, old_unknowns, new_unknowns, coord_ip, pe, C );

  if ( materi_velocity ) {
    for ( jdim=0; jdim<ndim; jdim++ ) {
      iuknwn = vel_indx + jdim * nder;
      divergence += grad_new_unknowns[jdim*nuknwn+iuknwn];
    }
  }
  if ( swit ) pri( "divergence", divergence );

  for ( inol=0; inol<nnol; inol++ ) {
    if ( groundflow_velocity ) {
      for ( jdim=0; jdim<ndim; jdim++ ) {
        iuknwn = (gvel_indx + jdim*nder);
        ipuknwn = iuknwn/nder;
        indx = inol*npuknwn + ipuknwn;
// made change in line above to correct for right sign of velocity
        tmp = h[inol] * (-pe[jdim]*grad_new_unknowns[jdim*nuknwn+pres_indx] -
          old_unknowns[iuknwn] ) / dtime;
        element_rhside[indx] += volume * tmp;
      }
    }
    if ( materidivergence==-YES ) {
      if ( materi_velocity ) {
        ipuknwn = pres_indx/nder;
        indx = inol*npuknwn + ipuknwn;
        tmp = h[inol] * divergence;
        element_rhside[indx] += volume * tmp;
        if ( residue ) element_residue[indx] += tmp;
        indxi = inol*npuknwn + ipuknwn;
        for ( jnol=0; jnol<nnol; jnol++ ) {
          for ( jdim=0; jdim<ndim; jdim++ ) {
            jpuknwn = vel_indx/nder + jdim;
            indxj = jnol*npuknwn + jpuknwn;
            element_matrix[indxi*nnol*npuknwn+indxj] -=
              volume * h[inol] * d[jdim*nnol+jnol];
          }
        }
      }
    }
  }
  if ( swit ) {
    pri( "element_rhside", element_rhside, nnol, npuknwn );
  }

  if ( swit ) pri( "Out function GROUNDFLOW" );

}

long int groundflow_phreatic_coord( long int inod, double coord[], double dof[],
  double &total_pressure, double &static_pressure, double &location )

// inod only in arguments for test printing

{

  long int length=0, found=0, ldum=0, idum[1], number[2];
  double water_level=0., dens=0., pressure_atmospheric=0., addtopressure=0.,
    ddum[1], force_gravity[MDIM], *groundflow_phreatic=NULL;

  force_gravity_calculate( force_gravity );
  db( GROUNDFLOW_DENSITY, 0, idum, &dens, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( GROUNDFLOW_PRESSURE_ATMOSPHERIC, 0, idum, &pressure_atmospheric, 
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( GROUNDFLOW_ADDTOPRESSURE, 0, idum, &addtopressure, ldum, VERSION_NORMAL, GET_IF_EXISTS );

  if ( groundflow_pressure ) {
    total_pressure = dof[pres_indx] - dens * force_gravity[ndim-1] * coord[ndim-1];
  }

  if ( db_active_index( GROUNDFLOW_PHREATICLEVEL, 0, VERSION_NORMAL ) ) {

    length = db_len( GROUNDFLOW_PHREATICLEVEL, 0, VERSION_NORMAL );
    groundflow_phreatic = db_dbl( GROUNDFLOW_PHREATICLEVEL, 0, VERSION_NORMAL ); 

    if      ( ndim==1 ) {
      if ( length!=1 ) db_error( GROUNDFLOW_PHREATICLEVEL, 0 );
      water_level = groundflow_phreatic[0];
      found = 1;
    }
    else if ( ndim==2 ) {
      if ( length==1 ) {
        found = 1;
        water_level = groundflow_phreatic[0];
      }
      else {
        found = table_xy( groundflow_phreatic, "GROUNDFLOW_PHREATICLEVEL",
          length, coord[0], water_level );
      }
    }
    else {
      assert( ndim==3 );
      if ( length==1 ) {
        found = 1;
        water_level = groundflow_phreatic[0];
      }
      else {
        db( GROUNDFLOW_PHREATICLEVEL_N, 0, number, ddum, ldum, VERSION_NORMAL, GET );
        if ( number[0]*number[1]*3 != length ) db_error( GROUNDFLOW_PHREATICLEVEL, 0 );
        found = table_xyz( groundflow_phreatic, number, coord, water_level );
      }
    }
    if ( found ) {
      location = water_level;
      if ( groundflow_pressure ) {
      	static_pressure = 
        force_gravity[ndim-1] * dens * ( water_level - coord[ndim-1] );
  	total_pressure = dof[pres_indx] + static_pressure; // added -> bug 
      }
    }
  }

  if ( static_pressure>=pressure_atmospheric ) static_pressure = pressure_atmospheric;
  if ( total_pressure>=pressure_atmospheric ) total_pressure = pressure_atmospheric;
  total_pressure += addtopressure;

  return found;
}

void groundflow_phreatic_apply( void )
 
{
  long int inod=0, max_node=0, iuknwn=0, ipuknwn=0, 
    groundflow_phreaticlevel_bounda=-NO, length=0, 
    node_phreaticlevel=0, ldum=0, idum[1], *node_bounded=NULL;
  double total_pressure=0, static_pressure=0., dens=0., location=0.,
    ddum[1], force_gravity[MDIM], *coord=NULL, *node_dof=NULL;
 
  if ( db_active_index( GROUNDFLOW_PHREATICLEVEL_BOUNDA, 0, VERSION_NORMAL ) ) {
    db( GROUNDFLOW_PHREATICLEVEL_BOUNDA, 0, &groundflow_phreaticlevel_bounda, ddum, 
      ldum, VERSION_NORMAL, GET );
    db( GROUNDFLOW_DENSITY, 0, idum, &dens, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    force_gravity_calculate( force_gravity );
    db_max_index( NODE, max_node, VERSION_NORMAL, GET );
    for ( inod=0; inod<=max_node; inod++ ) {
      if ( db_active_index( NODE_START_REFINED, inod, VERSION_NORMAL ) ) {
        coord = db_dbl( NODE_START_REFINED, inod, VERSION_NORMAL );
        node_dof = db_dbl( NODE_DOF, inod, VERSION_NEW );
        if ( groundflow_phreatic_coord( inod, coord, node_dof, total_pressure,
            static_pressure, location ) ) {
          iuknwn = pres_indx;
          ipuknwn = iuknwn / nder;
          if ( groundflow_phreaticlevel_bounda==-METHOD1 ) {
            if      ( scalar_dabs(coord[ndim-1]-location)<EPS_COORD ) {
              node_bounded = db_int( NODE_BOUNDED, inod, VERSION_NORMAL );
              node_dof[iuknwn] = dens * force_gravity[ndim-1] * location;
              node_bounded[ipuknwn] = 1;
            }
          }
          else if ( groundflow_phreaticlevel_bounda==-METHOD2 ) {
            if ( coord[ndim-1]>(location-EPS_COORD) ) {
              node_bounded = db_int( NODE_BOUNDED, inod, VERSION_NORMAL );
              node_dof[iuknwn] = dens * force_gravity[ndim-1] * location;
              node_bounded[ipuknwn] = 1;
            }
          }
          else {
            db_error( GROUNDFLOW_PHREATICLEVEL_BOUNDA, -1 );
          }
          if ( coord[ndim-1]>location+EPS_COORD ) {
            node_phreaticlevel = -ABOVE;
          }
          else {
            node_phreaticlevel = -BELOW;
          }
          length = 1;
          db( NODE_PHREATICLEVEL, inod, &node_phreaticlevel, ddum, 
            length, VERSION_NORMAL, PUT );
        }
      }
    }
  }
 
}
