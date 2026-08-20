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

void groundflow( long int element, long int gr, long int nnol, long int nodes[],
  double coord_ip[], double h[], double d[], 
  double volume, double old_unknowns[], 
  double new_unknowns[], double grad_new_unknowns[], 
  double element_matrix[], double element_rhside[],
  double element_residue[] )

{
  long int swit=0, inol=0, jnol=0, jdim=0, ipuknwn=0, iuknwn=0, jpuknwn=0,
    indx=0, indxi=0, indxj=0, icontrol=0, options_skip_groundflow_materidivergence=-NO,
    groundflow_consolidation_apply=-YES, control_groundflow_consolidation_apply=-YES,
    group_groundflow_consolidation_apply=-YES,
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
  db( GROUP_GROUNDFLOW_CONSOLIDATION_APPLY, gr, &group_groundflow_consolidation_apply,
    ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( GROUNDFLOW_CONSOLIDATION_APPLY, 0, &groundflow_consolidation_apply,
    ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( OPTIONS_SKIP_GROUNDFLOW_MATERIDIVERGENCE, 0, &options_skip_groundflow_materidivergence, 
    ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_GROUNDFLOW_CONSOLIDATION_APPLY, icontrol,
    &control_groundflow_consolidation_apply, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_OPTIONS_SKIP_GROUNDFLOW_MATERIDIVERGENCE, icontrol, 
    &options_skip_groundflow_materidivergence, ddum, ldum, 
    VERSION_NORMAL, GET_IF_EXISTS );
  if ( group_groundflow_consolidation_apply==-NO ) materidivergence = -NO;
  if ( groundflow_consolidation_apply==-NO ) materidivergence = -NO;
  if ( control_groundflow_consolidation_apply==-NO ) materidivergence = -NO;
  if ( options_skip_groundflow_materidivergence==-YES ) materidivergence = -NO;

  groundflow_data( element, gr, nodes, old_unknowns, new_unknowns, coord_ip, pe, C, h, nnol );

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

long int groundflow_phreatic_level_multiple_find( long int inod )

// Returns the index (imult) of the groundflow_phreatic_level_multiple record
// that owns the given node, or -1 if none. The domain of each record is
// selected by exactly one of _element, _element_group, _element_geometry or
// _node (they are not combinable, per the manual). The node belongs to the
// domain if it is listed in _node, or if one of its elements is in the
// _element/_element_group/_element_geometry list.

{
  long int imult=0, max_multiple=-1, found=-1, length=0, ldum=0,
    idum[1], nel=0, iel=0, elnum=0, gr=0, j=0,
    length_list=0, use_element=0, use_group=0, use_geometry=0, use_node=0,
    itmp=0, all=0, any=0;
  double ddum[1], rdum=0.;
  long int *node_element=NULL, *sel=NULL, *nodes=NULL, *el=NULL;

  db_max_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE, max_multiple, VERSION_NORMAL, GET );
  if ( max_multiple<0 ) return -1;

  nodes = get_new_int(MNOL);
  el = get_new_int(MNOL+1);
  for ( imult=0; imult<=max_multiple && found<0; imult++ ) {
    if ( !db_active_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE, imult, VERSION_NORMAL ) )
      continue;
    use_element = use_group = use_geometry = use_node = 0;
    if ( db_active_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT,
        imult, VERSION_NORMAL ) ) {
      use_element = 1;
      sel = db_int( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT, imult, VERSION_NORMAL );
      length_list = db_len( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT, imult, VERSION_NORMAL );
    }
    else if ( db_active_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GROUP,
        imult, VERSION_NORMAL ) ) {
      use_group = 1;
      sel = db_int( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GROUP, imult, VERSION_NORMAL );
      length_list = db_len( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GROUP,
        imult, VERSION_NORMAL );
    }
    else if ( db_active_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GEOMETRY,
        imult, VERSION_NORMAL ) ) {
      use_geometry = 1;
      sel = db_int( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GEOMETRY,
        imult, VERSION_NORMAL );
      length_list = db_len( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GEOMETRY,
        imult, VERSION_NORMAL );
    }
    else if ( db_active_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_NODE,
        imult, VERSION_NORMAL ) ) {
      use_node = 1;
      sel = db_int( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_NODE, imult, VERSION_NORMAL );
      length_list = db_len( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_NODE,
        imult, VERSION_NORMAL );
    }
    if ( use_node ) {
      if ( array_member( sel, inod, length_list, ldum ) ) found = imult;
    }
    else if ( use_element || use_group || use_geometry ) {
      node_element = db_int( NODE_ELEMENT, inod, VERSION_NORMAL );
      nel = db_len( NODE_ELEMENT, inod, VERSION_NORMAL );
      for ( iel=0; iel<nel && found<0; iel++ ) {
        elnum = node_element[iel];
        if ( use_element ) {
          if ( array_member( sel, elnum, length_list, ldum ) ) found = imult;
        }
        else if ( use_group ) {
          gr = 0;
          db( ELEMENT_GROUP, elnum, &gr, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
          if ( array_member( sel, gr, length_list, ldum ) ) found = imult;
        }
        else if ( use_geometry ) {
          db( ELEMENT, elnum, el, ddum, length, VERSION_NORMAL, GET );
          long int nnol_el = length - 1;
          for ( j=1; j<=nnol_el; j++ ) nodes[j-1] = el[j];
          all = 1; any = 0;
          for ( j=0; j<nnol_el; j++ ) {
            geometry( nodes[j], ddum, sel, itmp, rdum, ddum, rdum,
              ddum, NODE_START_REFINED, PROJECT_EXACT, VERSION_NORMAL );
            if ( !itmp ) all = 0;
            if ( itmp ) any = 1;
          }
          if ( all || any ) found = imult;
        }
      }
    }
  }
  delete[] nodes;
  delete[] el;

  return found;
}

long int groundflow_phreatic_coord( long int inod, double coord[], double dof[],
  double &total_pressure, double &static_pressure, double &location )

// inod only in arguments for test printing

{

  long int length=0, found=0, ldum=0, idum[1], number[2], imult=0;
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

  // groundflow_phreatic_level_multiple: several groundwater levels, each
  // owning a part of the domain (selected by _element/_element_group/
  // _element_geometry/_node). The node uses the level of its owning record.
  if ( db_active_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE, 0, VERSION_NORMAL ) ) {
    imult = groundflow_phreatic_level_multiple_find( inod );
    if ( imult>=0 ) {
      length = db_len( GROUNDFLOW_PHREATICLEVEL_MULTIPLE, imult, VERSION_NORMAL );
      groundflow_phreatic = db_dbl( GROUNDFLOW_PHREATICLEVEL_MULTIPLE, imult, VERSION_NORMAL );
      if      ( ndim==1 ) {
        if ( length!=1 ) db_error( GROUNDFLOW_PHREATICLEVEL_MULTIPLE, imult );
        water_level = groundflow_phreatic[0];
        found = 1;
      }
      else if ( ndim==2 ) {
        if ( length==1 ) {
          found = 1;
          water_level = groundflow_phreatic[0];
        }
        else {
          found = table_xy( groundflow_phreatic, "GROUNDFLOW_PHREATICLEVEL_MULTIPLE",
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
          db( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_N, imult, number, ddum, ldum,
            VERSION_NORMAL, GET );
          if ( number[0]*number[1]*3 != length )
            db_error( GROUNDFLOW_PHREATICLEVEL_MULTIPLE, imult );
          found = table_xyz( groundflow_phreatic, number, coord, water_level );
        }
      }
      if ( found ) {
        location = water_level;
        if ( groundflow_pressure ) {
          static_pressure = 
            force_gravity[ndim-1] * dens * ( water_level - coord[ndim-1] );
          total_pressure = dof[pres_indx] + static_pressure;
        }
      }
    }
  }
  else if ( db_active_index( GROUNDFLOW_PHREATICLEVEL, 0, VERSION_NORMAL ) ) {

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
    node_phreaticlevel=0, ldum=0, idum[1], *node_bounded=NULL,
    imult2=0, static_switch=-NO;
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

  // groundflow_phreatic_level_multiple_static: for the nodes of a multiple
  // phreatic level with _static -yes, set the total pressure (pore pressure)
  // equal to the static pressure. Convenient when the phreatic line is located
  // above the mesh part to which it belongs (no boundary condition can be
  // imposed), and to avoid solving the hydraulic heads (saves memory and CPU).
  if ( db_active_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE, 0, VERSION_NORMAL ) ) {
    long int max_multiple=0, iuknwn=0, ipuknwn=0;
    db_max_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE, max_multiple, VERSION_NORMAL, GET );
    db( GROUNDFLOW_DENSITY, 0, idum, &dens, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    force_gravity_calculate( force_gravity );
    for ( imult2=0; imult2<=max_multiple; imult2++ ) {
      if ( !db_active_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE, imult2, VERSION_NORMAL ) )
        continue;
      if ( !db_active_index( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_STATIC, imult2,
          VERSION_NORMAL ) )
        continue;
      db( GROUNDFLOW_PHREATICLEVEL_MULTIPLE_STATIC, imult2, &static_switch,
        ddum, ldum, VERSION_NORMAL, GET );
      if ( static_switch==-YES ) {
        db_max_index( NODE, max_node, VERSION_NORMAL, GET );
        for ( inod=0; inod<=max_node; inod++ ) {
          if ( db_active_index( NODE_START_REFINED, inod, VERSION_NORMAL ) ) {
            if ( groundflow_phreatic_level_multiple_find( inod )==imult2 ) {
              coord = db_dbl( NODE_START_REFINED, inod, VERSION_NORMAL );
              node_dof = db_dbl( NODE_DOF, inod, VERSION_NEW );
              if ( groundflow_phreatic_coord( inod, coord, node_dof, total_pressure,
                  static_pressure, location ) ) {
                iuknwn = pres_indx;
                ipuknwn = iuknwn / nder;
                node_dof[iuknwn] = static_pressure;
                node_bounded = db_int( NODE_BOUNDED, inod, VERSION_NORMAL );
                node_bounded[ipuknwn] = 1;
              }
            }
          }
        }
      }
    }
  }
 
}
