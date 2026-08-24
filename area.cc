/*
    Copyright (C) 1998  Dennis Roddeman
    email: dennis.roddeman@feat.nl

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is element_edge in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.


    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software Foundation 
    59 Temple Place, Suite 330, Boston, MA, 02111-1307, USA
*/

#include "tochnog.h"

#define MTYPES 9

// companion items for the convection/radiation edge families
// (condif_convection_edge_normal / condif_radiation_edge_normal, the
// Professional names of the legacy condif_convection / condif_radiation).
// which: 0=_ELEMENT 1=_ELEMENT_GROUP 2=_ELEMENT_SIDE 3=_NODE
// 4=_ELEMENT_NODE
static long int conv_rad_companion( long int master, long int which )
{
  static long int convection[] = {
    CONDIF_CONVECTION_EDGE_NORMAL_ELEMENT,
    CONDIF_CONVECTION_EDGE_NORMAL_ELEMENT_GROUP,
    CONDIF_CONVECTION_EDGE_NORMAL_ELEMENT_SIDE,
    CONDIF_CONVECTION_EDGE_NORMAL_NODE,
    CONDIF_CONVECTION_EDGE_NORMAL_ELEMENT_NODE };
  static long int radiation[] = {
    CONDIF_RADIATION_EDGE_NORMAL_ELEMENT,
    CONDIF_RADIATION_EDGE_NORMAL_ELEMENT_GROUP,
    CONDIF_RADIATION_EDGE_NORMAL_ELEMENT_SIDE,
    CONDIF_RADIATION_EDGE_NORMAL_NODE,
    CONDIF_RADIATION_EDGE_NORMAL_ELEMENT_NODE };
  return ( master==CONDIF_CONVECTION_EDGE_NORMAL ?
    convection[which] : radiation[which] );
}

static long int conv_rad_is_master( long int item )
{
  return item==CONDIF_CONVECTION_EDGE_NORMAL ||
    item==CONDIF_RADIATION_EDGE_NORMAL;
}

// companion items for the edge-normal flux families: the groundflow water
// flux (groundflow_flux_edge_normal) and the condif heat flux
// (condif_heat_edge_normal) share the same machinery. which selects:
// 0=_ELEMENT 1=_ELEMENT_GROUP 2=_ELEMENT_SIDE 3=_SINE 4=_TIME 5=_FACTOR
// 6=_NODE 7=_ELEMENT_NODE 8=_ELEMENT_NODE_FACTOR
static long int flux_edge_companion( long int master, long int which )
{
  static long int groundflow[] = {
    GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT,
    GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_GROUP,
    GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_SIDE,
    GROUNDFLOW_FLUX_EDGE_NORMAL_SINE,
    GROUNDFLOW_FLUX_EDGE_NORMAL_TIME,
    GROUNDFLOW_FLUX_EDGE_NORMAL_FACTOR,
    GROUNDFLOW_FLUX_EDGE_NORMAL_NODE,
    GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_NODE,
    GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_NODE_FACTOR };
  static long int condif[] = {
    CONDIF_HEAT_EDGE_NORMAL_ELEMENT,
    CONDIF_HEAT_EDGE_NORMAL_ELEMENT_GROUP,
    CONDIF_HEAT_EDGE_NORMAL_ELEMENT_SIDE,
    CONDIF_HEAT_EDGE_NORMAL_SINE,
    CONDIF_HEAT_EDGE_NORMAL_TIME,
    CONDIF_HEAT_EDGE_NORMAL_FACTOR,
    CONDIF_HEAT_EDGE_NORMAL_NODE,
    CONDIF_HEAT_EDGE_NORMAL_ELEMENT_NODE,
    CONDIF_HEAT_EDGE_NORMAL_ELEMENT_NODE_FACTOR };
  return ( master==GROUNDFLOW_FLUX_EDGE_NORMAL ? groundflow[which] : condif[which] );
}

static long int flux_edge_is_master( long int item )
{
  return item==GROUNDFLOW_FLUX_EDGE_NORMAL || item==CONDIF_HEAT_EDGE_NORMAL;
}

// companion items for the force_element_edge families (Professional
// force_edge_*): _ELEMENT _ELEMENT_GROUP _ELEMENT_SIDE _NODE
// _ELEMENT_NODE _NODE_FACTOR(DOUBLE, normal family only)
static long int force_edge_companion( long int master, long int which )
{
  static long int edge[] = {
    FORCE_ELEMENT_EDGE_ELEMENT,
    FORCE_ELEMENT_EDGE_ELEMENT_GROUP,
    FORCE_ELEMENT_EDGE_ELEMENT_SIDE,
    FORCE_ELEMENT_EDGE_NODE,
    FORCE_ELEMENT_EDGE_ELEMENT_NODE,
    -1 };
  static long int normal[] = {
    FORCE_ELEMENT_EDGE_NORMAL_ELEMENT,
    FORCE_ELEMENT_EDGE_NORMAL_ELEMENT_GROUP,
    FORCE_ELEMENT_EDGE_NORMAL_ELEMENT_SIDE,
    FORCE_ELEMENT_EDGE_NORMAL_NODE,
    FORCE_ELEMENT_EDGE_NORMAL_ELEMENT_NODE,
    FORCE_ELEMENT_EDGE_NORMAL_NODE_FACTOR };
  static long int water[] = {
    FORCE_ELEMENT_EDGE_WATER_ELEMENT,
    FORCE_ELEMENT_EDGE_WATER_ELEMENT_GROUP,
    FORCE_ELEMENT_EDGE_WATER_ELEMENT_SIDE,
    FORCE_ELEMENT_EDGE_WATER_NODE,
    FORCE_ELEMENT_EDGE_WATER_ELEMENT_NODE,
    -1 };
  if ( master==FORCE_ELEMENT_EDGE ) return edge[which];
  if ( master==FORCE_ELEMENT_EDGE_NORMAL ) return normal[which];
  if ( master==FORCE_ELEMENT_EDGE_WATER ) return water[which];
  return -1;
}

static long int force_edge_is_master( long int item )
{
  return item==FORCE_ELEMENT_EDGE || item==FORCE_ELEMENT_EDGE_NORMAL ||
    item==FORCE_ELEMENT_EDGE_WATER;
}

static long int border_nodes_tria3[] = {
    0, 1,
    1, 2,
    2, 0};

static long int border_nodes_tria6[] =  {
    0, 1, 2,
    2, 4, 5,
    5, 3, 0};

static long int border_nodes_quad4[] =  {
    0, 1,
    1, 3,
    3, 2,
    2, 0};

static long int border_nodes_quad9[] =  {
    0, 1, 2,
    2, 5, 8,
    8, 7, 6,
    6, 3, 0};

static long int border_nodes_quad16[] =  {
     0,  1,  2,  3,
     3,  7, 11, 15,
    15, 14, 13, 12,
    12,  8,  4,  0};

static long int border_nodes_tet4[] =  {
    0, 1, 2,
    0, 1, 3,
    0, 2, 3,
    1, 2, 3};

static long int border_nodes_hex8[] =  {
    0, 1, 2, 3,
    4, 5, 6, 7,
    0, 1, 4, 5,
    1, 3, 5, 7,
    2, 3, 6, 7,
    0, 2, 4, 6};

static long int border_nodes_hex27[] =  {
   0,  1,  2,  3,  4,  5,  6,  7,  8,
  18, 19, 20, 21, 22, 23, 24, 25, 26,
   0,  1,  2,  9, 10, 11, 18, 19, 20,
   2,  5,  8, 11, 14, 17, 20, 23, 26,
   6,  7,  8, 15, 16, 17, 24, 25, 26,
   0,  3,  6,  9, 12, 15, 18, 21, 24};

  // things with area integral, e.g. convection, radiation, element_edge force

void area( long int element, long int name, 
  long int gr, long int nnol, long int nodes[], 
  double new_coord[], double new_dof[], double element_lhside[], 
  double element_matrix[], double element_rhside[] )

{
  long int inol=0, jnol=0, knol=0, inol_side=0, max=0, ind=0, nnod=0, ipuknwn=0,
    iside=0, nside=0, ok=0, ok_tmp=0, itype=0, use_geom=0, swit=0, inod=0, idim=0,
    ind1=0, ind2=0, i=0, j=0, n=0, axisymmetric=-NO, any_area_integral=0,
    length=0, indx=0, nnol_side=0, iprinc=0, iuknwn=0, 
    all_under_phreatic_level=0, ifreq=0, nfreq=0, ldum=0, idum[1], 
    *sides, type[MTYPES], type_area[MTYPES], 
    geometry_entity[DATA_ITEM_SIZE], dof_principal[MUKNWN], *area=NULL;
  double alpha=0., temp=0., env_temp=0., ar=0., 
    frequency=0., amplitude=0., area_size=0., heat_flux=0.,
    a=0., b=0., c=0., radius=0., pressure=0.,
    dtime=0., time_start=0., time_current=0., time_total=0.,
    heat_flux_stiffness=0., load=0., delta_z=0., rdum=0., 
    factor=0., water_level=0., tmp=0., ddum[MDIM+MUKNWN],
    iso[MNOL], average_element_coord[MDIM], average_side_coord[MDIM],
    weight[MNOL], weight_tmp[MNOL], vec01[MDIM], vec02[MDIM],
    values[DATA_ITEM_SIZE], vec[MDIM], normal[MDIM], normal_tmp[MDIM],
    values_fac[DATA_ITEM_SIZE],
    *force_element_edge_time=NULL, *force_element_edge_normal_time=NULL,
    *force_element_edge_water_time=NULL, *force_element_edge_sine=NULL,
    *force_element_edge_normal_sine=NULL,
    *groundflow_flux_edge_normal_time=NULL, *groundflow_flux_edge_normal_sine=NULL;


  type[0] = CONDIF_RADIATION;
  type[1] = CONDIF_CONVECTION;
  type[2] = FORCE_ELEMENT_EDGE;
  type[3] = FORCE_ELEMENT_EDGE_NORMAL;
  type[4] = FORCE_ELEMENT_EDGE_WATER;
  type[5] = GROUNDFLOW_FLUX_EDGE_NORMAL;
  type[6] = CONDIF_HEAT_EDGE_NORMAL;
  type[7] = CONDIF_CONVECTION_EDGE_NORMAL;
  type[8] = CONDIF_RADIATION_EDGE_NORMAL;
  type_area[0] = CONDIF_RADIATION_GEOMETRY;
  type_area[1] = CONDIF_CONVECTION_GEOMETRY;
  type_area[2] = FORCE_ELEMENT_EDGE_GEOMETRY;
  type_area[3] = FORCE_ELEMENT_EDGE_NORMAL_GEOMETRY;
  type_area[4] = FORCE_ELEMENT_EDGE_WATER_GEOMETRY;
  type_area[5] = GROUNDFLOW_FLUX_EDGE_NORMAL_GEOMETRY;
  type_area[6] = CONDIF_HEAT_EDGE_NORMAL_GEOMETRY;
  type_area[7] = CONDIF_CONVECTION_EDGE_NORMAL_GEOMETRY;
  type_area[8] = CONDIF_RADIATION_EDGE_NORMAL_GEOMETRY;
  db( DOF_PRINCIPAL, 0, dof_principal, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );

  db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET_IF_EXISTS );
  db( TIME_CURRENT, 0, idum, &time_current, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( GROUP_AXISYMMETRIC, gr, &axisymmetric, ddum, ldum, 
    VERSION_NORMAL, GET_IF_EXISTS );
  time_total = time_current + dtime;

  for ( itype=0; itype<MTYPES; itype++ ) {
    if ( db_max_index( type[itype], max, VERSION_NORMAL, GET ) >=0 ) 
      any_area_integral = 1;
  }

  if ( any_area_integral ) {

    swit = set_swit(-1,-1,"area");
    if ( swit ) pri( "In routine AREA" );

    array_set( ddum, 0., MDIM+MUKNWN );

    if      ( name==-TRIA3 ) {
      nside = 3;
      nnol_side = 2;
      sides = border_nodes_tria3;
    }
    else if ( name==-TRIA6 ) {
      nside = 3;
      nnol_side = 3;
      sides = border_nodes_tria6;
    }
    else if ( name==-QUAD4 ) {
      nside = 4;
      nnol_side = 2;
      sides = border_nodes_quad4;
    }
    else if ( name==-QUAD9 ) {
      nside = 4;
      nnol_side = 3;
      sides = border_nodes_quad9;
    }
    else if ( name==-QUAD16 ) {
      nside = 4;
      nnol_side = 4;
      sides = border_nodes_quad16;
    }
    else if ( name==-TET4 ) {
      nside = 4;
      nnol_side = 3;
      sides = border_nodes_tet4;
    }
    else if ( name==-HEX8 ) {
      nside = 6;
      nnol_side = 4;
      sides = border_nodes_hex8;
    }
    else if ( name==-HEX27 ) {
      nside = 6;
      nnol_side = 9;
      sides = border_nodes_hex27;
    }
    else {
      nside = 0;
      sides = NULL;
    }
    array_set( average_element_coord, 0., ndim );
    for ( inol=0; inol<nnol; inol++ )
      array_add( &new_coord[inol*ndim], average_element_coord, average_element_coord, ndim );
    array_multiply( average_element_coord, average_element_coord, 1./nnol, ndim );

    for ( itype=0; itype<MTYPES; itype++ ) {
      if ( swit ) pri( "itype", (int)itype );

      db_max_index( type[itype], max, VERSION_NORMAL, GET );
      for ( ind=0; ind<=max; ind++ ) {
        if ( db_active_index( type[itype], ind, VERSION_NORMAL ) ) {
          if ( nside==0 ) {
            cout << "\nError: " << db_name( name ) << " not available for "; 
            cout << db_name( type[itype] ) << ".\n";
            exit(TN_EXIT_STATUS);
          }

          use_geom = nnod = 0;
          area = db_int( type_area[itype], ind, VERSION_NORMAL );
          if ( area[0]>0 ) {
            nnod = db_len( type_area[itype], ind, VERSION_NORMAL );
            renumbering_check( type_area[itype] );
          }
          else {
            db( type_area[itype], ind, geometry_entity, ddum, 
              ldum, VERSION_NORMAL, GET );
            use_geom = ( geometry_entity[0]<0 && 
              db_data_class(geometry_entity[0])==GEOMETRY );
            if ( !use_geom ) db_error( type[itype], ind );
          }
          if ( conv_rad_is_master(type[itype]) ) {
            // restriction variants for the convection/radiation families
            if ( db_active_index( conv_rad_companion(type[itype],0),
                ind, VERSION_NORMAL ) ) {
              long int elt[DATA_ITEM_SIZE], length_elt=0;
              db( conv_rad_companion(type[itype],0), ind, elt, ddum,
                length_elt, VERSION_NORMAL, GET );
              if ( !array_member( elt, element, length_elt, ldum ) ) continue;
            }
            if ( db_active_index( conv_rad_companion(type[itype],1),
                ind, VERSION_NORMAL ) ) {
              long int grp[DATA_ITEM_SIZE], length_grp=0;
              db( conv_rad_companion(type[itype],1), ind, grp, ddum,
                length_grp, VERSION_NORMAL, GET );
              if ( !array_member( grp, gr, length_grp, ldum ) ) continue;
            }
            if ( db_active_index( conv_rad_companion(type[itype],2),
                ind, VERSION_NORMAL ) ) {
              long int side_sel[DATA_ITEM_SIZE], length_side=0;
              db( conv_rad_companion(type[itype],2), ind, side_sel, ddum,
                length_side, VERSION_NORMAL, GET );
              long int ok_side = 0;
              for ( i=0; i+1<length_side; i+=2 )
                if ( side_sel[i]==element ) ok_side = 1;
              if ( !ok_side ) continue;
            }
          }
          if ( flux_edge_is_master(type[itype]) ) {
            // restriction variants for the edge-normal flux families
            if ( db_active_index( flux_edge_companion(type[itype],0),
                ind, VERSION_NORMAL ) ) {
              long int elt[DATA_ITEM_SIZE], length_elt=0;
              db( flux_edge_companion(type[itype],0), ind, elt, ddum,
                length_elt, VERSION_NORMAL, GET );
              if ( !array_member( elt, element, length_elt, ldum ) ) continue;
            }
            if ( db_active_index( flux_edge_companion(type[itype],1),
                ind, VERSION_NORMAL ) ) {
              long int grp[DATA_ITEM_SIZE], length_grp=0;
              db( flux_edge_companion(type[itype],1), ind, grp, ddum,
                length_grp, VERSION_NORMAL, GET );
              if ( !array_member( grp, gr, length_grp, ldum ) ) continue;
            }
            if ( db_active_index( flux_edge_companion(type[itype],2),
                ind, VERSION_NORMAL ) ) {
              long int side_sel[DATA_ITEM_SIZE], length_side=0;
              db( flux_edge_companion(type[itype],2), ind, side_sel, ddum,
                length_side, VERSION_NORMAL, GET );
              // pairs (element, side); skip this element if not listed
              long int ok_side = 0;
              for ( i=0; i+1<length_side; i+=2 )
                if ( side_sel[i]==element ) ok_side = 1;
              if ( !ok_side ) continue;
            }
          }

          if      ( type[itype]==FORCE_ELEMENT_EDGE ) {
            if      ( db_active_index( FORCE_ELEMENT_EDGE_SINE, ind, VERSION_NORMAL ) ) {
              force_element_edge_sine = db_dbl( FORCE_ELEMENT_EDGE_SINE, ind, VERSION_NORMAL );
              nfreq = ( db_len( FORCE_ELEMENT_EDGE_SINE, ind, VERSION_NORMAL ) - 1 ) / 2;
              time_start = force_element_edge_sine[0];
              load = 0.;
              if ( time_total>time_start ) {
                for ( ifreq=0; ifreq<nfreq; ifreq++ ) {
                  frequency = force_element_edge_sine[1+ifreq*2+0];
                  amplitude = force_element_edge_sine[1+ifreq*2+1];
                  load += amplitude * sin( 2. * PIRAD * frequency * time_total );
                }
              }
            }
            else if ( db_active_index( FORCE_ELEMENT_EDGE_TIME, ind, VERSION_NORMAL ) ) {
              force_element_edge_time = db_dbl( FORCE_ELEMENT_EDGE_TIME, 
                ind, VERSION_NORMAL );
              length = db_len( FORCE_ELEMENT_EDGE_TIME, ind, VERSION_NORMAL );
              force_time( force_element_edge_time, "FORCE_ELEMENT_EDGE_TIME",
                length, load );
            }
            else if ( db_active_index( FORCE_ELEMENT_EDGE_TIME_FILE, ind, VERSION_NORMAL ) ) {
	      long int force_time_file=0, ninc=0;
              db( FORCE_ELEMENT_EDGE_TIME_FILE, ind, &force_time_file, ddum, ldum, 
                 VERSION_NORMAL, GET_IF_EXISTS );
              if ( force_time_file==-YES ) force_time_file_apply( ind, FORCE_ELEMENT_EDGE_TIME_FILE, load );
              else db_error( FORCE_ELEMENT_EDGE_TIME_FILE, ind );
            }
            else
              load = 1.;
          }
          else if ( type[itype]==FORCE_ELEMENT_EDGE_NORMAL ) {
            if      ( db_active_index( FORCE_ELEMENT_EDGE_NORMAL_SINE, ind, VERSION_NORMAL ) ) {
              force_element_edge_normal_sine = db_dbl( FORCE_ELEMENT_EDGE_NORMAL_SINE, ind, VERSION_NORMAL );
              nfreq = ( db_len( FORCE_ELEMENT_EDGE_NORMAL_SINE, ind, VERSION_NORMAL ) - 1 ) / 2;
              time_start = force_element_edge_normal_sine[0];
              load = 0.;
              if ( time_total>time_start ) {
                for ( ifreq=0; ifreq<nfreq; ifreq++ ) {
                  frequency = force_element_edge_normal_sine[1+ifreq*2+0];
                  amplitude = force_element_edge_normal_sine[1+ifreq*2+1];
                  load += amplitude * sin( 2. * PIRAD * frequency * time_total );
                }
              }
            }
            else if ( db_active_index( FORCE_ELEMENT_EDGE_NORMAL_TIME, 
                ind, VERSION_NORMAL ) ) {
              force_element_edge_normal_time = db_dbl( FORCE_ELEMENT_EDGE_NORMAL_TIME, 
                ind, VERSION_NORMAL );
              length = db_len( FORCE_ELEMENT_EDGE_NORMAL_TIME, ind, VERSION_NORMAL );
              force_time( force_element_edge_normal_time,
                "FORCE_ELEMENT_EDGE_NORMAL_TIME", length, load );
            }
            else
              load = 1.;
          }
          else if ( type[itype]==FORCE_ELEMENT_EDGE_WATER ) {
            if ( db_active_index( FORCE_ELEMENT_EDGE_WATER_TIME, 
                ind, VERSION_NORMAL ) ) {
              force_element_edge_water_time = db_dbl( FORCE_ELEMENT_EDGE_WATER_TIME, 
                ind, VERSION_NORMAL );
              length = db_len( FORCE_ELEMENT_EDGE_WATER_TIME, ind, VERSION_NORMAL );
              force_time( force_element_edge_water_time,
                "FORCE_ELEMENT_EDGE_WATER_TIME", length, load );
            }
            else
              load = 1.;
          }
          else if ( flux_edge_is_master(type[itype]) ) {
            if ( db_active_index( flux_edge_companion(type[itype],3),
                ind, VERSION_NORMAL ) ) {
              groundflow_flux_edge_normal_sine = db_dbl(
                flux_edge_companion(type[itype],3), ind, VERSION_NORMAL );
              nfreq = ( db_len( flux_edge_companion(type[itype],3), ind,
                VERSION_NORMAL ) - 1 ) / 2;
              time_start = groundflow_flux_edge_normal_sine[0];
              load = 0.;
              if ( time_total>time_start ) {
                for ( ifreq=0; ifreq<nfreq; ifreq++ ) {
                  frequency = groundflow_flux_edge_normal_sine[1+ifreq*2+0];
                  amplitude = groundflow_flux_edge_normal_sine[1+ifreq*2+1];
                  load += amplitude * sin( 2. * PIRAD * frequency * time_total );
                }
              }
            }
            else if ( db_active_index( flux_edge_companion(type[itype],4),
                ind, VERSION_NORMAL ) ) {
              groundflow_flux_edge_normal_time = db_dbl(
                flux_edge_companion(type[itype],4), ind, VERSION_NORMAL );
              length = db_len( flux_edge_companion(type[itype],4), ind,
                VERSION_NORMAL );
              force_time( groundflow_flux_edge_normal_time,
                "EDGE_NORMAL_TIME", length, load );
            }
            else
              load = 1.;
          }
          if ( swit ) pri( "load", load );

          for ( iside=0; iside<nside; iside++ ) {
            all_under_phreatic_level = 1;
            ok = 1;
            array_set( normal, 0., ndim );
            array_set( average_side_coord, 0., ndim );
            for ( inol_side=0; inol_side<nnol_side && ok; inol_side++ ) {
              inol = sides[iside*nnol_side + inol_side]; 
              inod = nodes[inol];
              array_add( &new_coord[inol*ndim], average_side_coord,
                average_side_coord, ndim );
              if      ( use_geom ) {
                geometry( inod, ddum, geometry_entity, ok_tmp, rdum, normal_tmp, rdum,
                  ddum, NODE_START_REFINED, PROJECT_EXACT, VERSION_NORMAL );
                array_add( normal_tmp, normal, normal, ndim );
              }
              else {
                assert( nnod>0 );
                ok_tmp = array_member( area, inod, nnod, ldum );
              }
              if ( !ok_tmp ) ok = 0;
              if ( type[itype]==FORCE_ELEMENT_EDGE_WATER ) {
                if ( groundflow_phreatic_coord( inod, &new_coord[inol*ndim], 
                    ddum, ddum[0], ddum[0], water_level ) ) {
                  db( FORCE_ELEMENT_EDGE_WATER, ind, idum, values, 
                    ldum, VERSION_NORMAL, GET );
                  if ( new_coord[inol*ndim+ndim-1]>=water_level+EPS_COORD ) 
                    all_under_phreatic_level = 0;
                }
                else {
                  pri( "Error: no valid phreatic level for FORCE_ELEMENT_EDGE_WATER." );
                  exit(TN_EXIT_STATUS);
                }
              }
              else
                all_under_phreatic_level = 0;
            }
            array_multiply( average_side_coord, average_side_coord, 1./nnol_side, ndim );
            array_subtract( average_side_coord, average_element_coord, vec, ndim );
            if ( ok ) {
              if ( !use_geom ) {
                if ( ndim==2 ) {
                  inol = sides[iside * nnol_side + 0];
                  jnol = sides[iside * nnol_side + 1];
                  array_subtract( &new_coord[inol*ndim], &new_coord[jnol*ndim],
                    vec01, ndim );
                  array_outproduct_2D( vec01, normal );
                }
                else {
                  inol = sides[iside * nnol_side + 0];
                  jnol = sides[iside * nnol_side + 1];
                  knol = sides[iside * nnol_side + 2];
                  array_subtract( &new_coord[inol*ndim], &new_coord[jnol*ndim],
                    vec01, ndim );
                  array_subtract( &new_coord[inol*ndim], &new_coord[knol*ndim],
                    vec02, ndim );
                  array_outproduct_3D( vec01, vec02, normal );
                }
              }
              if ( array_inproduct( normal, vec, ndim ) < 0. )
                array_multiply( normal, normal, -1., ndim );
              array_normalize( normal, ndim );
              if ( ndim==2 ) {
                array_subtract( &new_coord[sides[iside * nnol_side + 0]*ndim], 
                    &new_coord[sides[iside * nnol_side + nnol_side-1]*ndim], vec, ndim );
                ar = array_size( vec, ndim );
                integration_lobatto( nnol_side, iso, weight );
              }
              else {
                assert( ndim==3 );
                if ( name==-TET4 ) {
                  inol = sides[iside * nnol_side + 0];
                  jnol = sides[iside * nnol_side + 1];
                  knol = sides[iside * nnol_side + 2];
                  array_subtract( &new_coord[inol*ndim], &new_coord[jnol*ndim], 
                    vec, ndim );
                  a = array_size( vec, ndim );
                  array_subtract( &new_coord[inol*ndim], &new_coord[knol*ndim], 
                    vec, ndim );
                  b = array_size( vec, ndim );
                  array_subtract( &new_coord[jnol*ndim], &new_coord[knol*ndim], 
                    vec, ndim );
                  c = array_size( vec, ndim );
                  ar = sqrt( (a+b+c)*(a+b-c)*(a-b+c)*(-a+b+c)/16 );
                  weight[0] = 1./3.;
                  weight[1] = 1./3.;
                  weight[2] = 1./3.;
                }
                else {
                  if ( name==-HEX8 ) {
                    ind1 = 1;
                    ind2 = 2;
                    n = 2;
                  }
                  else {
                    assert( name==-HEX27 );
                    ind1 = 2;
                    ind2 = 6;
                    n = 3;
                  }
                  assert( n*n==nnol_side );
                  ar = 
                    triangle_area( &new_coord[sides[iside * nnol_side + 0]*ndim],
                      &new_coord[sides[iside * nnol_side + ind1]*ndim],
                      &new_coord[sides[iside * nnol_side + ind2]*ndim] );
                  ar +=  triangle_area( &new_coord[sides[iside * nnol_side + ind1]*ndim],
                      &new_coord[sides[iside * nnol_side + ind2]*ndim],
                      &new_coord[sides[iside * nnol_side + nnol_side-1]*ndim] );
                  integration_lobatto( n, iso, weight_tmp );
                  for ( j=0; j<n; j++ ) {
                    for ( i=0; i<n; i++ ) {
                      weight[j*n+i] = weight_tmp[i]*weight_tmp[j];
                    }
                  }
                }
              }
              if ( swit ) {
                pri( "ar", ar );
                pri( "weight", weight, nnol_side );
              }
              for ( inol_side=0; inol_side<nnol_side; inol_side++ ) {
                inol = sides[iside * nnol_side + inol_side]; 
                inod = nodes[inol];
                if ( axisymmetric==-YES ) {
                  radius = new_coord[inol*ndim];
                  area_size = ar * 2. * PIRAD * radius;
                }
                else
                  area_size = ar;
                if ( type[itype]==CONDIF_RADIATION ||
                     type[itype]==CONDIF_CONVECTION ||
                     type[itype]==CONDIF_RADIATION_EDGE_NORMAL ||
                     type[itype]==CONDIF_CONVECTION_EDGE_NORMAL ) {
                  // node restrictions for the Professional-name masters
                  long int use_it = 1;
          if ( force_edge_is_master(type[itype]) ) {
            // element-level restrictions for the force_element_edge
            // families (Professional force_edge_*)
            if ( force_edge_companion(type[itype],0)>=0 &&
                 db_active_index( force_edge_companion(type[itype],0),
                     ind, VERSION_NORMAL ) ) {
              long int elt[DATA_ITEM_SIZE], length_elt=0;
              db( force_edge_companion(type[itype],0), ind, elt, ddum,
                length_elt, VERSION_NORMAL, GET );
              if ( !array_member( elt, element, length_elt, ldum ) ) continue;
            }
            if ( force_edge_companion(type[itype],1)>=0 &&
                 db_active_index( force_edge_companion(type[itype],1),
                     ind, VERSION_NORMAL ) ) {
              long int grp[DATA_ITEM_SIZE], length_grp=0;
              db( force_edge_companion(type[itype],1), ind, grp, ddum,
                length_grp, VERSION_NORMAL, GET );
              if ( !array_member( grp, gr, length_grp, ldum ) ) continue;
            }
            if ( force_edge_companion(type[itype],2)>=0 &&
                 db_active_index( force_edge_companion(type[itype],2),
                     ind, VERSION_NORMAL ) ) {
              long int side_sel[DATA_ITEM_SIZE], length_side=0;
              db( force_edge_companion(type[itype],2), ind, side_sel, ddum,
                length_side, VERSION_NORMAL, GET );
              long int ok_side = 0;
              for ( i=0; i+1<length_side; i+=2 )
                if ( side_sel[i]==element ) ok_side = 1;
              if ( !ok_side ) continue;
            }
          }
          if ( conv_rad_is_master(type[itype]) ) {
                    if ( db_active_index( conv_rad_companion(type[itype],3),
                        ind, VERSION_NORMAL ) ) {
                      long int nds[DATA_ITEM_SIZE], length_nds=0;
                      db( conv_rad_companion(type[itype],3), ind, nds, ddum,
                        length_nds, VERSION_NORMAL, GET );
                      if ( !array_member( nds, inod, length_nds, ldum ) )
                        use_it = 0;
                    }
                    if ( db_active_index( conv_rad_companion(type[itype],4),
                        ind, VERSION_NORMAL ) ) {
                      long int en[DATA_ITEM_SIZE], length_en=0;
                      db( conv_rad_companion(type[itype],4), ind, en, ddum,
                        length_en, VERSION_NORMAL, GET );
                      // en[0]=element, en[1..]=local node numbers
                      if ( en[0]!=element || !array_member( &en[1], inol,
                          length_en-1, ldum ) ) use_it = 0;
                    }
                  }
                  if ( !use_it ) continue;
                  db( type[itype], ind, idum, values,
                    ldum, VERSION_NORMAL, GET );
                  alpha = values[0]; env_temp = values[1];
                  if ( swit ) {
                    pri( "alpha", alpha );
                    pri( "env_temp", env_temp );
                  }
                  temp = new_dof[inol*nuknwn+temp_indx];
                  if ( swit ) pri( "temp", temp );
                  if ( type[itype]==CONDIF_RADIATION ||
                       type[itype]==CONDIF_RADIATION_EDGE_NORMAL ) {
                    heat_flux = alpha * weight[inol_side] * area_size *
                      (scalar_power(env_temp,4)-scalar_power(temp,4));
                    if ( swit ) pri( "heat_flux", heat_flux );
                    heat_flux_stiffness = alpha * weight[inol_side] * area_size *
                      4.*scalar_power(temp,3);
                    if ( swit ) pri( "heat_flux_stiffness", heat_flux_stiffness );
                  }
                  else {
                    heat_flux = alpha * weight[inol_side] * area_size * (env_temp-temp);
                    heat_flux_stiffness = alpha * weight[inol_side] * area_size;
                  }
                  if ( swit ) {
                    pri( "inol", inol );
                    pri( "temp", temp );
                    pri( "heat_flux_stiffness", heat_flux_stiffness );
                    pri( "heat_flux", heat_flux );
                  }
                  indx = inol*npuknwn + temp_indx/nder;
                  element_lhside[indx] += heat_flux_stiffness;
                  element_matrix[indx*nnol*npuknwn+indx] += heat_flux_stiffness;
                  element_rhside[inol*npuknwn+temp_indx/nder] +=
                    heat_flux;
                }
                else if ( type[itype]==FORCE_ELEMENT_EDGE ) {
                  // node restrictions (Professional force_edge_*)
                  long int use_it = 1;
                  if ( force_edge_companion(type[itype],3)>=0 &&
                       db_active_index( force_edge_companion(type[itype],3),
                           ind, VERSION_NORMAL ) ) {
                    long int nds[DATA_ITEM_SIZE], length_nds=0;
                    db( force_edge_companion(type[itype],3), ind, nds, ddum,
                      length_nds, VERSION_NORMAL, GET );
                    if ( !array_member( nds, inod, length_nds, ldum ) )
                      use_it = 0;
                  }
                  if ( force_edge_companion(type[itype],4)>=0 &&
                       db_active_index( force_edge_companion(type[itype],4),
                           ind, VERSION_NORMAL ) ) {
                    long int en[DATA_ITEM_SIZE], length_en=0;
                    db( force_edge_companion(type[itype],4), ind, en, ddum,
                      length_en, VERSION_NORMAL, GET );
                    if ( en[0]!=element || !array_member( &en[1], inol,
                        length_en-1, ldum ) ) use_it = 0;
                  }
                  double node_factor = 1.;
                  if ( db_active_index( FORCE_ELEMENT_EDGE_NODE_FACTOR,
                      ind, VERSION_NORMAL ) ) {
                    long int length_nf=0;
                    double values_nf[DATA_ITEM_SIZE];
                    db( FORCE_ELEMENT_EDGE_NODE_FACTOR, ind,
                      idum, values_nf, length_nf, VERSION_NORMAL, GET );
                    long int nf_el = (long int)values_nf[0], jnf;
                    if ( nf_el==element )
                      for ( jnf=0; jnf+1<length_nf; jnf++ )
                        if ( jnf==inol ) node_factor = values_nf[jnf+1];
                  }
                  if ( !use_it ) continue;
                  force_factor( FORCE_ELEMENT_EDGE_FACTOR, ind,
                    &new_coord[inol*ndim], factor );
                  db( FORCE_ELEMENT_EDGE, ind, idum, values, 
                    ldum, VERSION_NORMAL, GET );
                  iprinc = 0;
                  for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
                    iuknwn = ipuknwn*nder;
                    if ( dof_principal[iuknwn]>=0 ) {
                      element_rhside[inol*npuknwn+ipuknwn] += factor * node_factor *
                        load * weight[inol_side] * area_size * values[iprinc];
                      iprinc++;
                    }
                  }
                }
                else if ( type[itype]==FORCE_ELEMENT_EDGE_NORMAL ) {
                  // node restrictions + per-node factor (Professional
                  // force_edge_normal_node_factor)
                  long int use_it = 1;
                  double node_factor = 1.;
                  if ( db_active_index( force_edge_companion(type[itype],3),
                      ind, VERSION_NORMAL ) ) {
                    long int nds[DATA_ITEM_SIZE], length_nds=0;
                    db( force_edge_companion(type[itype],3), ind, nds, ddum,
                      length_nds, VERSION_NORMAL, GET );
                    if ( !array_member( nds, inod, length_nds, ldum ) )
                      use_it = 0;
                  }
                  if ( db_active_index( force_edge_companion(type[itype],4),
                      ind, VERSION_NORMAL ) ) {
                    long int en[DATA_ITEM_SIZE], length_en=0;
                    db( force_edge_companion(type[itype],4), ind, en, ddum,
                      length_en, VERSION_NORMAL, GET );
                    if ( en[0]!=element || !array_member( &en[1], inol,
                        length_en-1, ldum ) ) use_it = 0;
                  }
                  if ( db_active_index( FORCE_ELEMENT_EDGE_NORMAL_NODE_FACTOR,
                      ind, VERSION_NORMAL ) ) {
                    long int length_nf=0;
                    double values_nf[DATA_ITEM_SIZE];
                    db( FORCE_ELEMENT_EDGE_NORMAL_NODE_FACTOR, ind,
                      idum, values_nf, length_nf, VERSION_NORMAL, GET );
                    long int nf_el = (long int)values_nf[0], jnf;
                    if ( nf_el==element )
                      for ( jnf=0; jnf+1<length_nf; jnf++ )
                        if ( jnf==inol ) node_factor = values_nf[jnf+1];
                  }
                  if ( !use_it ) continue;
                  force_factor( FORCE_ELEMENT_EDGE_NORMAL_FACTOR, ind,
                    &new_coord[inol*ndim], factor );
                  db( FORCE_ELEMENT_EDGE_NORMAL, ind, idum, values,
                    ldum, VERSION_NORMAL, GET );
                  for ( idim=0; idim<ndim; idim++ ) {
                    ipuknwn = vel_indx/nder + idim;
                    tmp = factor * node_factor *
                      load * weight[inol_side] * area_size *
                      values[0] * normal[idim];
                    element_rhside[inol*npuknwn+ipuknwn] += tmp;
                  }
                }
                else if ( type[itype]==FORCE_ELEMENT_EDGE_WATER ) {
                  // node restrictions (Professional force_edge_water_*)
                  long int use_it = 1;
                  if ( db_active_index( force_edge_companion(type[itype],3),
                      ind, VERSION_NORMAL ) ) {
                    long int nds[DATA_ITEM_SIZE], length_nds=0;
                    db( force_edge_companion(type[itype],3), ind, nds, ddum,
                      length_nds, VERSION_NORMAL, GET );
                    if ( !array_member( nds, inod, length_nds, ldum ) )
                      use_it = 0;
                  }
                  if ( db_active_index( force_edge_companion(type[itype],4),
                      ind, VERSION_NORMAL ) ) {
                    long int en[DATA_ITEM_SIZE], length_en=0;
                    db( force_edge_companion(type[itype],4), ind, en, ddum,
                      length_en, VERSION_NORMAL, GET );
                    if ( en[0]!=element || !array_member( &en[1], inol,
                        length_en-1, ldum ) ) use_it = 0;
                  }
                  if ( !use_it ) continue;
                  if ( all_under_phreatic_level ) {
                    double water_factor = 1.;
                    force_factor( FORCE_ELEMENT_EDGE_WATER_FACTOR, ind,
                      &new_coord[inol*ndim], water_factor );
                    groundflow_phreatic_coord( inod, &new_coord[inol*ndim],
                      ddum, ddum[0], ddum[0], water_level );
                    db( FORCE_ELEMENT_EDGE_WATER, ind, idum, values,
                      ldum, VERSION_NORMAL, GET );
                    array_normalize( &values[2], ndim );
                    delta_z = water_level - new_coord[inol*ndim+ndim-1];
                    pressure = values[0] * values[1] * delta_z;
                    for ( idim=0; idim<ndim; idim++ ) {
                      ipuknwn = vel_indx/nder + idim;
                      tmp = load * weight[inol_side] * area_size *
                        pressure * values[2+idim] * water_factor;
                      element_rhside[inol*npuknwn+ipuknwn] += tmp;
                    }
                  }
                }
                else if ( flux_edge_is_master(type[itype]) ) {
                  force_factor( flux_edge_companion(type[itype],5), ind,
                    &new_coord[inol*ndim], factor );
                  db( type[itype], ind, idum, values,
                    ldum, VERSION_NORMAL, GET );
                  // restriction by node
                  long int use_it = 1;
                  if ( db_active_index( flux_edge_companion(type[itype],6),
                      ind, VERSION_NORMAL ) ) {
                    long int nds[DATA_ITEM_SIZE], length_nds=0;
                    db( flux_edge_companion(type[itype],6), ind, nds, ddum,
                      length_nds, VERSION_NORMAL, GET );
                    if ( !array_member( nds, inod, length_nds, ldum ) ) use_it = 0;
                  }
                  if ( db_active_index( flux_edge_companion(type[itype],7),
                      ind, VERSION_NORMAL ) ) {
                    long int en[DATA_ITEM_SIZE], length_en=0;
                    db( flux_edge_companion(type[itype],7), ind, en, ddum,
                      length_en, VERSION_NORMAL, GET );
                    // en[0]=element, en[1..]=local node numbers
                    if ( en[0]!=element || !array_member( &en[1], inol,
                        length_en-1, ldum ) ) use_it = 0;
                  }
                  double node_factor = 1.;
                  if ( db_active_index( flux_edge_companion(type[itype],8),
                      ind, VERSION_NORMAL ) ) {
                    long int length_enf=0;
                    db( flux_edge_companion(type[itype],8), ind,
                      idum, values_fac, length_enf, VERSION_NORMAL, GET );
                    // values_fac[0]=element, values_fac[1..]=factors for the
                    // local nodes of that element
                    long int enf_el = (long int)values_fac[0];
                    long int j;
                    if ( enf_el==element ) {
                      for ( j=0; j+1<length_enf; j++ )
                        if ( j==inol ) node_factor = values_fac[j+1];
                    }
                  }
                  if ( use_it ) {
                    if ( type[itype]==GROUNDFLOW_FLUX_EDGE_NORMAL )
                      ipuknwn = pres_indx/nder;
                    else
                      ipuknwn = temp_indx/nder;
                    tmp = factor * node_factor *
                      load * weight[inol_side] * area_size * values[0];
                    element_rhside[inol*npuknwn+ipuknwn] += tmp;
                  }
                }
                else
                  db_error( type[itype], ind );
              }
            }
          }
        }
      }
    }

    if ( swit ) pri( "Out routine AREA" );
  }

}

void area_node_dataitem( )

{
  long int inod=0, max_node=0, found=0, iarea=0, max_area=0,
    length=0, ldum=0, data_item_name=0, idum[1], 
    geometry_entity[2], area_node_dataitem[3],
    *area_node_dataitem_integer=NULL;
  double rdum=0., ddum[MDIM], *area_node_dataitem_double=NULL;

  area_node_dataitem_integer = get_new_int(DATA_ITEM_SIZE);
  area_node_dataitem_double = get_new_dbl(DATA_ITEM_SIZE);
  db_max_index( AREA_NODE_DATAITEM, max_area, VERSION_NORMAL, GET );
  for ( iarea=0; iarea<=max_area; iarea++ ) {
    if ( db_active_index( AREA_NODE_DATAITEM, iarea, VERSION_NORMAL ) ) {
      db( AREA_NODE_DATAITEM, iarea, area_node_dataitem, ddum, ldum, 
        VERSION_NORMAL, GET );
      array_move( area_node_dataitem, geometry_entity, 2 );
      data_item_name = area_node_dataitem[2];
      if ( db_type(data_item_name)==INTEGER )
        db( AREA_NODE_DATAITEM_INTEGER, iarea, 
          area_node_dataitem_integer, ddum, length, VERSION_NORMAL, GET );
      else
        db( AREA_NODE_DATAITEM_DOUBLE, iarea, idum, 
          area_node_dataitem_double, length, VERSION_NORMAL, GET );
      db_max_index( NODE, max_node, VERSION_NORMAL, GET );
      for ( inod=0; inod<=max_node; inod++ ) {
        if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
          geometry( inod, ddum, geometry_entity, found, rdum, ddum, rdum,
            ddum, NODE_START_REFINED, PROJECT_EXACT, VERSION_NORMAL );
          if ( found ) {
            db( data_item_name, inod, area_node_dataitem_integer, 
              area_node_dataitem_double, length, VERSION_NORMAL, PUT );
          }
        }
      }
    }
  }
  delete[] area_node_dataitem_integer;
  delete[] area_node_dataitem_double;

}
