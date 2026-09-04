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

#define MTYPES 11

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
  static long int projected[] = {
    FORCE_ELEMENT_EDGE_PROJECTED_ELEMENT,
    FORCE_ELEMENT_EDGE_PROJECTED_ELEMENT_GROUP,
    FORCE_ELEMENT_EDGE_PROJECTED_ELEMENT_SIDE,
    FORCE_ELEMENT_EDGE_PROJECTED_NODE,
    FORCE_ELEMENT_EDGE_PROJECTED_ELEMENT_NODE,
    FORCE_ELEMENT_EDGE_PROJECTED_NODE_FACTOR };
  if ( master==FORCE_ELEMENT_EDGE ) return edge[which];
  if ( master==FORCE_ELEMENT_EDGE_NORMAL ) return normal[which];
  if ( master==FORCE_ELEMENT_EDGE_WATER ) return water[which];
  if ( master==FORCE_ELEMENT_EDGE_PROJECTED ) return projected[which];
  return -1;
}

// Sprint 13 lot 1: the support_edge_normal family companions (same
// index): element / element_group / element_side (element-level) and
// node / element_node (node-level)
static long int support_edge_companion( long int which )
{
  static long int sup[] = {
    SUPPORT_EDGE_NORMAL_ELEMENT,
    SUPPORT_EDGE_NORMAL_ELEMENT_GROUP,
    SUPPORT_EDGE_NORMAL_ELEMENT_SIDE,
    SUPPORT_EDGE_NORMAL_NODE,
    SUPPORT_EDGE_NORMAL_ELEMENT_NODE,
    -1 };
  return sup[which];
}

static long int force_edge_is_master( long int item )
{
  return item==FORCE_ELEMENT_EDGE || item==FORCE_ELEMENT_EDGE_NORMAL ||
    item==FORCE_ELEMENT_EDGE_WATER || item==FORCE_ELEMENT_EDGE_PROJECTED;
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
  type[9] = FORCE_ELEMENT_EDGE_PROJECTED;
  type[10] = SUPPORT_EDGE_NORMAL;
  type_area[0] = CONDIF_RADIATION_GEOMETRY;
  type_area[1] = CONDIF_CONVECTION_GEOMETRY;
  type_area[2] = FORCE_ELEMENT_EDGE_GEOMETRY;
  type_area[3] = FORCE_ELEMENT_EDGE_NORMAL_GEOMETRY;
  type_area[4] = FORCE_ELEMENT_EDGE_WATER_GEOMETRY;
  type_area[5] = GROUNDFLOW_FLUX_EDGE_NORMAL_GEOMETRY;
  type_area[6] = CONDIF_HEAT_EDGE_NORMAL_GEOMETRY;
  type_area[7] = CONDIF_CONVECTION_EDGE_NORMAL_GEOMETRY;
  type_area[8] = CONDIF_RADIATION_EDGE_NORMAL_GEOMETRY;
  type_area[9] = FORCE_ELEMENT_EDGE_PROJECTED_GEOMETRY;
  type_area[10] = SUPPORT_EDGE_NORMAL_GEOMETRY;
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

    // Sprint 13 lot 1: re-zero the node_support_edge_normal_force
    // records ONCE per area() call at the start of each assembly sweep
    // (the element loop runs ascending and area() is called once per
    // element per sweep: a non-increasing element number means a new
    // sweep started). The records accumulate over the elements sharing
    // a node within the sweep. Single-threaded only: with
    // OPTIONS_PROCESSORS > 1 the interleaving breaks the accumulation
    // (the support forces themselves stay correct) - warned once.
    {
      static long int s_last_element = -1;
      static long int s_warned_thread = 0;
      if ( db_max_index( SUPPORT_EDGE_NORMAL, max, VERSION_NORMAL, GET )
           >=0 ) {
        long int nthread = 1;
        db( OPTIONS_PROCESSORS, 0, &nthread, ddum, ldum,
          VERSION_NORMAL, GET_IF_EXISTS );
        if ( nthread>1 && !s_warned_thread ) {
          pri( "Warning: node_support_edge_normal_force accumulation "
               "requires OPTIONS_PROCESSORS 1 (the support forces "
               "themselves are correct)" );
          s_warned_thread = 1;
        }
        if ( element<=s_last_element ) {
          long int max_nd = -1, ind2;
          db_max_index( NODE_SUPPORT_EDGE_NORMAL_FORCE, max_nd,
            VERSION_NORMAL, GET );
          for ( ind2=0; ind2<=max_nd; ind2++ )
            if ( db_active_index( NODE_SUPPORT_EDGE_NORMAL_FORCE, ind2,
                 VERSION_NORMAL ) ) {
              double z[MDIM];
              array_set( z, 0., MDIM );
              db( NODE_SUPPORT_EDGE_NORMAL_FORCE, ind2, idum, z, ndim,
                VERSION_NEW, PUT );
            }
        }
        s_last_element = element;
      }
    }

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
          else if ( type[itype]==FORCE_ELEMENT_EDGE_PROJECTED ) {
            if ( db_active_index( FORCE_ELEMENT_EDGE_PROJECTED_SINE,
                ind, VERSION_NORMAL ) ) {
              groundflow_flux_edge_normal_sine = db_dbl(
                FORCE_ELEMENT_EDGE_PROJECTED_SINE, ind, VERSION_NORMAL );
              nfreq = ( db_len( FORCE_ELEMENT_EDGE_PROJECTED_SINE, ind,
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
            else if ( db_active_index( FORCE_ELEMENT_EDGE_PROJECTED_TIME,
                ind, VERSION_NORMAL ) ) {
              groundflow_flux_edge_normal_time = db_dbl(
                FORCE_ELEMENT_EDGE_PROJECTED_TIME, ind, VERSION_NORMAL );
              length = db_len( FORCE_ELEMENT_EDGE_PROJECTED_TIME, ind,
                VERSION_NORMAL );
              force_time( groundflow_flux_edge_normal_time,
                "FORCE_ELEMENT_EDGE_PROJECTED_TIME", length, load );
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
          if ( type[itype]==SUPPORT_EDGE_NORMAL ) {
            // support_edge_normal_time (manual Professional 6.1084):
            // a multiplication factor for the support edge force
            // (linear interpolation between the diagram points)
            if ( db_active_index( SUPPORT_EDGE_NORMAL_TIME, ind,
                 VERSION_NORMAL ) ) {
              double *sup_time = db_dbl( SUPPORT_EDGE_NORMAL_TIME, ind,
                VERSION_NORMAL );
              length = db_len( SUPPORT_EDGE_NORMAL_TIME, ind,
                VERSION_NORMAL );
              force_time( sup_time, "SUPPORT_EDGE_NORMAL_TIME",
                length, load );
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
              // The load/flux direction of the edge families is the
              // PHYSICAL side normal, computed from the side node
              // coordinates (not from the selector geometry: a
              // geometry_point selector returns the RADIAL normal of
              // each side node, which for a side not aligned with the
              // point tilts the accumulated normal away from the true
              // face normal - corpus force12 loads the face between
              // group 1 and group 2 via geometry_point with a huge
              // tolerance).
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
              // element_side i  side_0 element_1 side_1...:
              // BOTH the element AND the (1-based local) side must
              // match for the side to be selected
              for ( i=0; i+1<length_side; i+=2 )
                if ( side_sel[i]==element &&
                     side_sel[i+1]==iside+1 ) ok_side = 1;
              if ( !ok_side ) continue;
            }
            // element_node (Professional manual 6.1072): node-level
            // restriction by (element, local node indices)
            if ( force_edge_companion(type[itype],4)>=0 &&
                 db_active_index( force_edge_companion(type[itype],4),
                     ind, VERSION_NORMAL ) ) {
              long int en[DATA_ITEM_SIZE], length_en=0;
              db( force_edge_companion(type[itype],4), ind, en, ddum,
                length_en, VERSION_NORMAL, GET );
              long int ok_en = 0;
              // element_node i en_0 en_1...:
              for ( i=0; i+1<length_en; i+=2 )
                if ( en[i]==element &&
                     array_member( &en[i+1], inol, length_en-i-1, ldum ) )
                  ok_en = 1;
              if ( !ok_en ) use_it = 0;
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
                  // one value per space direction (Professional
                  // force_edge, manual 6.454); the number of stored
                  // values (ldum) bounds the principal dof mapping
                  iprinc = 0;
                  for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
                    iuknwn = ipuknwn*nder;
                    if ( dof_principal[iuknwn]>=0 ) {
                      if ( iprinc>=ldum ) break;
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
                else if ( type[itype]==SUPPORT_EDGE_NORMAL ) {
                  // Sprint 13 lot 1: the distributed Winkler support
                  // (manual Professional 6.1067). The support force is
                  // computed from the TOTAL DISPLACEMENTS of the side
                  // nodes (the Professional: "supports calculate forces
                  // directly from total displacements") and applied as
                  // a consistent nodal force of the side quadrature
                  // (Lobatto: the nodes themselves), per unit length in
                  // 2D / unit area in 3D:
                  //   f = -( k_n * u_n ) * n - k_t * ( u - u_n * n )
                  // with n the OUTWARD side normal (into the support:
                  // u_n > 0 compresses the support and it pushes back
                  // along -n). Element-level restrictions (element /
                  // element_group / element_side) and node-level
                  // restrictions (node / element_node), same index.
                  {
                    long int use_it = 1;
                    if ( db_active_index( support_edge_companion(0),
                        ind, VERSION_NORMAL ) ) {
                      long int elt[DATA_ITEM_SIZE], length_elt=0;
                      db( support_edge_companion(0), ind, elt, ddum,
                        length_elt, VERSION_NORMAL, GET );
                      if ( !array_member( elt, element, length_elt, ldum ) )
                        continue;
                    }
                    if ( db_active_index( support_edge_companion(1),
                        ind, VERSION_NORMAL ) ) {
                      long int grp[DATA_ITEM_SIZE], length_grp=0;
                      db( support_edge_companion(1), ind, grp, ddum,
                        length_grp, VERSION_NORMAL, GET );
                      if ( !array_member( grp, gr, length_grp, ldum ) )
                        continue;
                    }
                    if ( db_active_index( support_edge_companion(2),
                        ind, VERSION_NORMAL ) ) {
                      long int side_sel[DATA_ITEM_SIZE], length_side=0;
                      db( support_edge_companion(2), ind, side_sel, ddum,
                        length_side, VERSION_NORMAL, GET );
                      long int ok_side = 0;
                      for ( i=0; i+1<length_side; i+=2 )
                        if ( side_sel[i]==element &&
                             side_sel[i+1]==iside+1 ) ok_side = 1;
                      if ( !ok_side ) continue;
                    }
                    if ( db_active_index( support_edge_companion(3),
                        ind, VERSION_NORMAL ) ) {
                      long int nds[DATA_ITEM_SIZE], length_nds=0;
                      db( support_edge_companion(3), ind, nds, ddum,
                        length_nds, VERSION_NORMAL, GET );
                      if ( !array_member( nds, inod, length_nds, ldum ) )
                        use_it = 0;
                    }
                    if ( db_active_index( support_edge_companion(4),
                        ind, VERSION_NORMAL ) ) {
                      long int en[DATA_ITEM_SIZE], length_en=0;
                      db( support_edge_companion(4), ind, en, ddum,
                        length_en, VERSION_NORMAL, GET );
                      // en[0]=element, en[1..]=local node numbers
                      if ( en[0]!=element || !array_member( &en[1], inol,
                          length_en-1, ldum ) ) use_it = 0;
                    }
                    if ( !use_it ) continue;
                    if ( !materi_displacement ) {
                      pri( "Error: support_edge_normal needs "
                        "materi_displacement (the support force is "
                        "computed from total displacements)" );
                      exit(TN_EXIT_STATUS);
                    }
                    db( SUPPORT_EDGE_NORMAL, ind, idum, values,
                      ldum, VERSION_NORMAL, GET );
                    {
                      double un = 0., vn = 0., an = 0., f0 = 0.,
                        sup_nodal[MDIM], old_dof_node[MDIM+MUKNWN];
                      double cn = 0., ct = 0., fac_k = 1.;
                      long int idim2 = 0;
                      // the spatial factor scales the STIFFNESSES only
                      // (manual 6.1075: "for the support stiffnesses
                      // and not the force")
                      force_factor( SUPPORT_EDGE_NORMAL_FACTOR, ind,
                        &new_coord[inol*ndim], fac_k );
                      // the velocity and acceleration of the node
                      // (acceleration: backward difference of the
                      // velocity dofs; the old values from NODE_DOF
                      // VERSION_NORMAL)
                      array_move( db_dbl( NODE_DOF, inod,
                        VERSION_NORMAL ), old_dof_node, nuknwn );
                      for ( idim=0; idim<ndim; idim++ ) {
                        un += new_dof[inol*nuknwn+dis_indx+idim*nder]
                             *normal[idim];
                        vn += new_dof[inol*nuknwn+vel_indx+idim*nder]
                             *normal[idim];
                        an += ( new_dof[inol*nuknwn+vel_indx+idim*nder]
                              - old_dof_node[vel_indx+idim*nder] )
                              /dtime*normal[idim];
                      }
                      // the initial normal force (manual 6.1076):
                      // a0 + a1*depth, a compression preload of the
                      // support (the reaction pushes the element even
                      // at zero displacement)
                      if ( db_active_index( SUPPORT_EDGE_NORMAL_FORCE_INITIAL,
                           ind, VERSION_NORMAL ) ) {
                        double fi[2];
                        db( SUPPORT_EDGE_NORMAL_FORCE_INITIAL, ind,
                          idum, fi, ldum, VERSION_NORMAL, GET );
                        f0 = fi[0] + fi[1]
                          *new_coord[inol*ndim+ndim-1];
                      }
                      // the damping coefficients (manual 6.1068):
                      // viscous dampers on the support. Gated by
                      // control_support_edge_normal_damping_apply -no
                      // (manual 6.380). The automatic variants
                      // (6.1069/6.1070) compute them from the
                      // attached element group: c_n = Cn*rho*Vn with
                      // Cn = 1, Vn = sqrt(Eoed/rho), Eoed =
                      // (1-nu)E/((1+nu)(1-2nu)); c_t = Ct*rho*Vt with
                      // Ct = 0.25, Vt = sqrt(G/rho), G = E/(2(1+nu)).
                      // The _apparent variant uses the apparent
                      // moduli from the CURRENT nodal state (for
                      // elastic behavior identical to the nominal
                      // values; guards fall back to nominal).
                      {
                        long int damping_apply = -YES;
                        long int icontrol_d = 0;
                        db( ICONTROL, 0, &icontrol_d, ddum, ldum,
                          VERSION_NORMAL, GET_IF_EXISTS );
                        db( CONTROL_SUPPORT_EDGE_NORMAL_DAMPING_APPLY,
                          icontrol_d, &damping_apply, ddum, ldum,
                          VERSION_NORMAL, GET_IF_EXISTS );
                        if ( damping_apply!=-NO ) {
                          if ( db_active_index(
                               SUPPORT_EDGE_NORMAL_DAMPING, ind,
                               VERSION_NORMAL ) ) {
                            double cdamp[2];
                            db( SUPPORT_EDGE_NORMAL_DAMPING, ind,
                              idum, cdamp, ldum, VERSION_NORMAL, GET );
                            cn = cdamp[0];
                            ct = cdamp[1];
                          }
                          else if ( db_active_index(
                               SUPPORT_EDGE_NORMAL_DAMPING_AUTOMATIC,
                               ind, VERSION_NORMAL ) ||
                                    db_active_index(
                               SUPPORT_EDGE_NORMAL_DAMPING_AUTOMATIC_APPARENT,
                               ind, VERSION_NORMAL ) ) {
                            double e_mod = 0., nu = 0., rho = 0.;
                            db( GROUP_MATERI_ELASTI_YOUNG, gr, idum,
                              &e_mod, ldum, VERSION_NORMAL,
                              GET_IF_EXISTS );
                            db( GROUP_MATERI_ELASTI_POISSON, gr, idum,
                              &nu, ldum, VERSION_NORMAL,
                              GET_IF_EXISTS );
                            db( GROUP_MATERI_DENSITY, gr, idum, &rho,
                              ldum, VERSION_NORMAL, GET_IF_EXISTS );
                            if ( rho<=0. ) {
                              pri( "Warning: support_edge_normal_"
                                "damping_automatic needs "
                                "group_materi_density > 0 on the "
                                "supported group - no damping "
                                "applied" );
                            }
                            else if ( e_mod>0. && nu<0.5 ) {
                              if ( db_active_index(
                                   SUPPORT_EDGE_NORMAL_DAMPING_AUTOMATIC_APPARENT,
                                   ind, VERSION_NORMAL ) ) {
                                // the apparent moduli: from the
                                // current nodal stress/strain along
                                // the loading axis, with guards
                                double eps =
                                  new_dof[inol*nuknwn+ept_indx
                                    +(ndim-1)*nder];
                                double sig =
                                  new_dof[inol*nuknwn+stres_indx
                                    +(ndim-1)*nder];
                                if ( scalar_dabs(eps)>1.e-12 &&
                                     sig/eps>0. && sig/eps<e_mod )
                                  e_mod = sig/eps;
                              }
                              double eoed = (1.-nu)*e_mod
                                /((1.+nu)*(1.-2.*nu));
                              double g_mod = e_mod/(2.*(1.+nu));
                              cn = 1.    *sqrt( rho*eoed );
                              ct = 0.25 *sqrt( rho*g_mod );
                            }
                          }
                        }
                      }
                      // the distributed support density (manual
                      // 6.1071): inertia of the support mass
                      double dn = 0., dt_ = 0.;
                      if ( db_active_index( SUPPORT_EDGE_NORMAL_DENSITY,
                           ind, VERSION_NORMAL ) ) {
                        double densv[2];
                        db( SUPPORT_EDGE_NORMAL_DENSITY, ind, idum,
                          densv, ldum, VERSION_NORMAL, GET );
                        dn = densv[0];
                        dt_ = densv[1];
                      }
                      // the consistent nodal support force (Lobatto:
                      // the side nodes carry the line/area load),
                      // all terms times the time factor
                      // -- pre-loop: compute the RAW support reaction
                      // and tangential vector and apply the PLASTIC
                      // caps before the per-direction loop. My sign:
                      // S = k*u_n + c*v_n + d*a_n + f0 (POSITIVE =
                      // compression, manual convention 6.1079); fn
                      // = -S is the NORMAL REACTION (positive = pushes
                      // the element OUT of the support).
                      double fn_raw = -( fac_k*values[0]*un
                        + cn*vn + dn*an + f0 );
                      double ft_raw_local[MDIM];
                      double ft_mag2_local = 0.;
                      // the TANGENTIAL components use the tangential
                      // displacement u_t = u - (u.n)n, velocity
                      // v_t = v - (v.n)n and acceleration a_t, NOT the
                      // full vectors (the normal part is carried by fn)
                      for ( long int idim_p = 0; idim_p<ndim; idim_p++ ) {
                        double ut_p =
                          new_dof[inol*nuknwn+dis_indx+idim_p*nder]
                          - un*normal[idim_p];
                        double vt_p =
                          new_dof[inol*nuknwn+vel_indx+idim_p*nder]
                          - vn*normal[idim_p];
                        double at_p =
                          ( new_dof[inol*nuknwn+vel_indx+idim_p*nder]
                          - old_dof_node[vel_indx+idim_p*nder] )
                          /dtime - an*normal[idim_p];
                        ft_raw_local[idim_p] = -( fac_k*values[1]*ut_p
                          + ct*vt_p + dt_*at_p );
                        ft_mag2_local += ft_raw_local[idim_p]
                          *ft_raw_local[idim_p];
                      }
                      // the PLASTICITY records (manual 6.1079-6.1082)
                      long int pt_gap = 0, has_tension_max = 0,
                        has_comp_min = 0, has_tang_fac = 0,
                        has_friction = 0;
                      double comp_min_mag = 1.e300, tang_factor = 1.e300,
                        tension_max = -1.e300, cohesion = 0.,
                        friction_mu = 0.;
                      if ( db_active_index(
                           SUPPORT_EDGE_NORMAL_PLASTI_TENSION, ind,
                           VERSION_NORMAL ) ) {
                        long int sw_pt = 0;
                        long int *ptr2 = db_int(
                          SUPPORT_EDGE_NORMAL_PLASTI_TENSION, ind,
                          VERSION_NORMAL );
                        sw_pt = *ptr2;
                        {
                        const char *sw_str = (const char *)ptr2;
                        // the parser stores -yes/-no as ival = -YES (the
                        // NEGATIVE of the YES enum ordinal); check the
                        // matching pattern
                        if ( sw_pt == -YES ) pt_gap = 1;
                        }
                      }
                      if ( db_active_index(
                           SUPPORT_EDGE_NORMAL_PLASTI_TENSION_DOUBLE, ind,
                           VERSION_NORMAL ) ) {
                        long int need = 1;
                        double ptmax_v = 0.;
                        if ( db_active_index(
                             SUPPORT_EDGE_NORMAL_PLASTI_TENSION_DOUBLE,
                             ind, VERSION_NEW ) ) {
                          double *ptr2 = db_dbl(
                            SUPPORT_EDGE_NORMAL_PLASTI_TENSION_DOUBLE,
                            ind, VERSION_NEW );
                          ptmax_v = *ptr2;
                        }
                        else {
                          db( SUPPORT_EDGE_NORMAL_PLASTI_TENSION_DOUBLE,
                            ind, idum, &ptmax_v, ldum,
                            VERSION_NORMAL, GET );
                        }
                        (void)need;
                        tension_max = ptmax_v;
                        has_tension_max = 1;
                      }
                      if ( db_active_index(
                           SUPPORT_EDGE_NORMAL_PLASTI_COMPRESSION, ind,
                           VERSION_NORMAL ) ) {
                        double pc[2];
                        db( SUPPORT_EDGE_NORMAL_PLASTI_COMPRESSION, ind,
                          idum, pc, ldum, VERSION_NORMAL, GET );
                        comp_min_mag = scalar_dabs( pc[0] );
                        has_comp_min = 1;
                        tang_factor = pc[1];
                        has_tang_fac = 1;
                      }
                      if ( db_active_index(
                           SUPPORT_EDGE_NORMAL_PLASTI_FRICTION, ind,
                           VERSION_NORMAL ) ) {
                        double pf[2];
                        db( SUPPORT_EDGE_NORMAL_PLASTI_FRICTION, ind,
                          idum, pf, ldum, VERSION_NORMAL, GET );
                        cohesion = pf[0];
                        friction_mu = pf[1];
                        has_friction = 1;
                      }
                      // apply the PLASTICITY caps (manual 6.1079-6.1082)
                      double S = -fn_raw;
                      long int side_opened = 0;
                      if ( has_tension_max ) {
                        // tension cap (manual 6.1082): when the
                        // tension exceeded, S is clipped and the
                        // tangential force is zeroed
                        if ( S<-tension_max ) S = -tension_max;
                      }
                      if ( has_comp_min ) {
                        // compression floor: the compression cannot
                        // exceed the floor (manual 6.1079)
                        if ( S>comp_min_mag ) S = comp_min_mag;
                      }
                      if ( pt_gap && S<0. ) {
                        // tension gap (manual 6.1081): when the support
                        // is in tension, all forces are zero
                        S = 0.; ft_mag2_local = 0.;
                        for ( long int idim_q = 0;
                          idim_q<ndim; idim_q++ )
                          ft_raw_local[idim_q] = 0.;
                        side_opened = 1;
                      }
                      if ( has_tension_max && S==-tension_max ) {
                        // tension cap zeroes the tangential force
                        for ( long int idim_r = 0; idim_r<ndim;
                          idim_r++ )
                          ft_raw_local[idim_r] = 0.;
                        ft_mag2_local = 0.;
                      }
                      // the tangential magnitude caps (manual 6.1079
                      // compression tangential factor and 6.1080
                      // Coulomb friction): the smallest of the two
                      // governs
                      if ( has_tang_fac || has_friction ) {
                        double ft_cap = 1.e300;
                        if ( has_tang_fac )
                          ft_cap = tang_factor*scalar_dabs(S);
                        if ( has_friction ) {
                          double f_lim = cohesion
                            + friction_mu*scalar_dabs(S);
                          if ( f_lim<ft_cap ) ft_cap = f_lim;
                        }
                        if ( ft_mag2_local>ft_cap*ft_cap &&
                             ft_mag2_local>0. ) {
                          double scale = ft_cap/sqrt(ft_mag2_local);
                          for ( long int idim_s = 0; idim_s<ndim;
                            idim_s++ )
                            ft_raw_local[idim_s] *= scale;
                        }
                      }
                      // the final capped values: fn (from S) and
                      // ft_raw_local[idim] for each direction
                      double fn = -S;
                      // now the per-direction loop builds fidim from
                      // the CAPPED fn and the CAPPED ft_raw_local
                      // element-level restrictions (support_edge_normal_
                      // element_side + _element_node, same index)
                      {
                        long int ok_side_l = 1, ok_node_l = 1;
                        if ( support_edge_companion(2)>=0 &&
                             db_active_index(
                                 support_edge_companion(2),
                                 ind, VERSION_NORMAL ) ) {
                          long int side_sel[DATA_ITEM_SIZE],
                            length_side=0;
                          db( support_edge_companion(2), ind, side_sel, ddum,
                            length_side, VERSION_NORMAL, GET );
                          ok_side_l = 0;
                          for ( long int ii=0; ii+1<length_side; ii+=2 )
                            if ( side_sel[ii]==element &&
                                 side_sel[ii+1]==iside+1 ) ok_side_l = 1;
                        }
                        if ( support_edge_companion(4)>=0 &&
                             db_active_index(
                                 support_edge_companion(4),
                                 ind, VERSION_NORMAL ) ) {
                          long int en[DATA_ITEM_SIZE], length_en=0;
                          db( support_edge_companion(4), ind, en, ddum,
                            length_en, VERSION_NORMAL, GET );
                          ok_node_l = 0;
                          for ( long int ii=0; ii+1<length_en; ii+=2 )
                            if ( en[ii]==element &&
                                 array_member( &en[ii+1], inol,
                                   length_en-ii-1, ldum ) )
                              ok_node_l = 1;
                        }
                        if ( !ok_side_l || !ok_node_l ) continue;
                      }
                      for ( idim=0; idim<ndim; idim++ ) {
                        double ut_idim =
                          new_dof[inol*nuknwn+dis_indx+idim*nder]
                          - un*normal[idim];
                        double vt_idim =
                          new_dof[inol*nuknwn+vel_indx+idim*nder]
                          - vn*normal[idim];
                        double at_idim =
                          ( new_dof[inol*nuknwn+vel_indx+idim*nder]
                          - old_dof_node[vel_indx+idim*nder] )
                          /dtime - an*normal[idim];
                        // fidim uses the CAPPED fn and ft_raw_local
                        // (the gap/tension-cap branch zeros ft_raw_local
                        // and the per-side prefactor ft_raw_local is
                        // the only thing the tang force depends on)

                        double fidim = fn*normal[idim]
                          + ft_raw_local[idim];

                        sup_nodal[idim] =
                          load * weight[inol_side] * area_size * fidim;
                        ipuknwn = vel_indx/nder + idim;
                        element_rhside[inol*npuknwn+ipuknwn]
                          += sup_nodal[idim];
                      }
                      // the output record
                      // node_support_edge_normal_plasti_tension_status
                      // (manual 6.895): 0 = closed (no gap), 1 =
                      // opened (tension gap fired at this side).
                      // Same per-sweep re-zeroing as the force record
                      // (accumulated per node with OR semantics).
                      if ( db_active_index(
                           NODE_SUPPORT_EDGE_NORMAL_PLASTI_TENSION_STATUS,
                           inod, VERSION_NORMAL ) ) {
                        long int acc_t = 0;
                        if ( db_active_index(
                             NODE_SUPPORT_EDGE_NORMAL_PLASTI_TENSION_STATUS,
                             inod, VERSION_NEW ) ) {
                          long int *ptr = db_int(
                            NODE_SUPPORT_EDGE_NORMAL_PLASTI_TENSION_STATUS,
                            inod, VERSION_NEW );
                          acc_t = *ptr;
                        }
                        if ( side_opened>0.5 ) acc_t = 1;
                        {
                          long int *ptr = db_int(
                            NODE_SUPPORT_EDGE_NORMAL_PLASTI_TENSION_STATUS,
                            inod, VERSION_NEW );
                          *ptr = acc_t;
                        }
                      }
                      // the consistent support stiffness in the
                      // MATRIX (and the diagonal element_lhside):
                      // dt * int ( k_n n(x)n + k_t (I - n(x)n) ) N_i
                      // N_j dA over the side. Without it a body
                      // resting on the support keeps a zero-energy
                      // rigid mode in the velocity matrix and the
                      // solver breaks down (measured); with it the
                      // fixed point is unchanged (the equilibrium is
                      // in the RHS force; the matrix term only
                      // carries the iteration - the Professional's
                      // plasti_residual_stiffness documentation states
                      // the same for the plastic case). Same dtime
                      // scaling as the element stiffness (materi.cc:
                      // volume*dtime*stiffness).
                      // support_edge_normal_plasti_residual_stiffness
                      // (manual 6.1083): a fraction factor in [0,1] of the
                      // original elastic stiffness is added to the matrix
                      // for stability. Default factor=1 means NO extra
                      // stiffness beyond the elastic baseline.
                      {
                        double fac_res = 0.;
                        if ( db_active_index(
                             SUPPORT_EDGE_NORMAL_PLASTI_RESIDUAL_STIFFNESS,
                             ind, VERSION_NORMAL ) ) {
                          double res_v = 0.;
                          db( SUPPORT_EDGE_NORMAL_PLASTI_RESIDUAL_STIFFNESS,
                            ind, idum, &res_v, ldum,
                            VERSION_NORMAL, GET );
                          fac_res = res_v;
                        }
                        {
                          double gxi[3], gw[3];
                        long int ngs = ( nnol_side<3 ? 2 : 3 ), igs=0,
                          jgs=0, inol_i=0, inol_j=0, ia=0, ib=0;
                        long int n1 = ( ndim==2 ? nnol_side
                          : (long int) sqrt((double)nnol_side) );
                        // 2D: the side is a 1D interval [-1,1],
                        // ds = ar/2. 3D: the face is the tensor
                        // [-1,1]^2, dA = ar/4; the border-table face
                        // node order is tensor with the first
                        // direction fastest (the same order the 3D
                        // Lobatto weights above use).
                        double measfac = ( ndim==2 ? ar/2. : ar/4. );
                        double cij[MNOL*MNOL];
                        for ( inol_i=0; inol_i<nnol_side; inol_i++ )
                          for ( inol_j=0; inol_j<nnol_side; inol_j++ )
                            cij[inol_i*nnol_side+inol_j] = 0.;
                        integration_gauss( ngs, gxi, gw );
                        for ( igs=0; igs<ngs; igs++ ) {
                          double hi1[MNOL], pdum[MPOINT];
                          interpolation_polynomial( gxi[igs], n1,
                            hi1, pdum );
                          for ( jgs=0; jgs<( ndim==2 ? 1 : ngs );
                            jgs++ ) {
                            double hj1[MNOL];
                            double wt = gw[igs]*measfac;
                            interpolation_polynomial(
                              ( ndim==2 ? 0. : gxi[jgs] ), n1,
                              hj1, pdum );
                            if ( ndim==3 ) wt *= gw[jgs];
                            for ( inol_i=0; inol_i<nnol_side;
                              inol_i++ ) {
                              double hii = ( ndim==2 ? hi1[inol_i] :
                                hi1[inol_i%n1]*hj1[inol_i/n1] );
                              for ( inol_j=0; inol_j<nnol_side;
                                inol_j++ ) {
                                double hjj = ( ndim==2 ?
                                  hi1[inol_j] :
                                  hi1[inol_j%n1]*hj1[inol_j/n1] );
                                cij[inol_i*nnol_side+inol_j]
                                  += wt*hii*hjj;
                              }
                            }
                          }
                        }
                        for ( inol_i=0; inol_i<nnol_side; inol_i++ ) {
                          long int node_i =
                            sides[iside*nnol_side+inol_i];
                          for ( inol_j=0; inol_j<nnol_side; inol_j++ ) {
                            long int node_j =
                              sides[iside*nnol_side+inol_j];
                            for ( ia=0; ia<ndim; ia++ ) {
                              for ( ib=0; ib<ndim; ib++ ) {
                                double kterm = dtime *
                                  ( values[0]*normal[ia]*normal[ib]
                                    + values[1]*( (ia==ib?1.:0.)
                                      - normal[ia]*normal[ib] ) ) *
                                  cij[inol_i*nnol_side+inol_j];
                                long int indxi = node_i*npuknwn
                                  + vel_indx/nder + ia;
                                long int indxj = node_j*npuknwn
                                  + vel_indx/nder + ib;
                                element_matrix[indxi*nnol*npuknwn+indxj]
                                  += kterm;
                                  // support_edge_normal_plasti_residual_
                                  // stiffness (manual 6.1083): a fraction
                                  // factor in [0,1] of the original
                                  // elastic stiffness is added to the matrix
                                  // for stability; default 1 (no extra
                                  // stiffness beyond the elastic baseline).
                                  if ( fac_res>0. ) {
                                    double rterm = dtime*fac_res*
                                      ( values[0]*normal[ia]*normal[ib]
                                       + values[1]*((ia==ib?1.:0.)
                                         - normal[ia]*normal[ib]) )*
                                      cij[inol_i*nnol_side+inol_j];
                                    element_matrix[indxi*nnol*npuknwn+indxj]
                                      += rterm;
                                    if ( inol_i==inol_j && ia==ib )
                                      element_lhside[indxi] += rterm;
                                  }
                                if ( inol_i==inol_j && ia==ib )
                                  element_lhside[indxi] += kterm;
                              }
                            }
                          }
                        }
                      }
                      // the output record
                      // node_support_edge_normal_force: the consistent
                      // nodal support force, accumulated over the
                      // elements sharing the node. Re-zeroed at EVERY
                      // assembly sweep (detected by the element number
                      // NOT increasing - one call per element per
                      // sweep, the loop runs ascending). Single
                      // threaded only: with OPTIONS_PROCESSORS > 1 the
                      // interleaving breaks the accumulation (the
                      // support forces themselves stay correct) -
                      // warned once.
                      if ( db_active_index( NODE_SUPPORT_EDGE_NORMAL_FORCE,
                           inod, VERSION_NORMAL ) ) {
                        double acc[MDIM];
                        array_set( acc, 0., MDIM );
                        db( NODE_SUPPORT_EDGE_NORMAL_FORCE, inod,
                          idum, acc, ldum, VERSION_NEW, GET_IF_EXISTS );
                        for ( idim=0; idim<ndim; idim++ )
                          acc[idim] += sup_nodal[idim];
                        db( NODE_SUPPORT_EDGE_NORMAL_FORCE, inod,
                          idum, acc, ndim, VERSION_NEW, PUT );
                      }
                      }
                    }
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
                else if ( type[itype]==FORCE_ELEMENT_EDGE_PROJECTED ) {
                  // force_edge_projected (manual Professional 6.478):
                  // Terzaghi tunnel load. A linear ground stress field
                  //   ph = ph0 + ph_grad . x   (horizontal, ⊥ tunnel axis)
                  //   pv = pv0 + pv_grad . x   (vertical, along v_dir)
                  // is projected on the edge: with the edge outward
                  // normal n and tangent t,
                  //   sig_radial     = ph (n.hd)^2 + pv (n.vd)^2
                  //   sig_tangential = ph (t.hd)(n.hd) + pv (t.vd)(n.vd)
                  // applied as factor_normal*sig_radial*n +
                  // factor_tangential*sig_tangential*t (the release load
                  // pushes the excavation boundary outward into the void).
                  // hd = tunnel_dir x v_dir (2D: out-of-plane x v_dir).
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
                  if ( db_active_index( FORCE_ELEMENT_EDGE_PROJECTED_NODE_FACTOR,
                      ind, VERSION_NORMAL ) ) {
                    long int length_nf=0;
                    double values_nf[DATA_ITEM_SIZE];
                    db( FORCE_ELEMENT_EDGE_PROJECTED_NODE_FACTOR, ind,
                      idum, values_nf, length_nf, VERSION_NORMAL, GET );
                    long int nf_el = (long int)values_nf[0], jnf;
                    if ( nf_el==element )
                      for ( jnf=0; jnf+1<length_nf; jnf++ )
                        if ( jnf==inol ) node_factor = values_nf[jnf+1];
                  }
                  if ( !use_it ) continue;
                  force_factor( FORCE_ELEMENT_EDGE_PROJECTED_FACTOR, ind,
                    &new_coord[inol*ndim], factor );
                  {
                    double pdata[16], ph=0., pv=0., fn=1., ft=0.;
                    long int idir;
                    double vd[MDIM], hd[MDIM], tang[MDIM], n2=0., t2=0.;
                    db( FORCE_ELEMENT_EDGE_PROJECTED, ind, idum, pdata,
                      ldum, VERSION_NORMAL, GET );
                    // stress field at the node
                    if ( ndim==2 ) {
                      ph = pdata[0] + pdata[1]*new_coord[inol*ndim+0]
                                   + pdata[2]*new_coord[inol*ndim+1];
                      pv = pdata[3] + pdata[4]*new_coord[inol*ndim+0]
                                   + pdata[5]*new_coord[inol*ndim+1];
                      fn = pdata[6]; ft = pdata[7];
                      vd[0] = pdata[8]; vd[1] = pdata[9]; vd[2] = 0.;
                      // 2D: tunnel axis out-of-plane; hd = z_axis x vd
                      hd[0] = -vd[1]; hd[1] = vd[0]; hd[2] = 0.;
                    }
                    else {
                      ph = pdata[0]  + pdata[1]*new_coord[inol*ndim+0]
                                    + pdata[2]*new_coord[inol*ndim+1]
                                    + pdata[3]*new_coord[inol*ndim+2];
                      pv = pdata[4]  + pdata[5]*new_coord[inol*ndim+0]
                                    + pdata[6]*new_coord[inol*ndim+1]
                                    + pdata[7]*new_coord[inol*ndim+2];
                      fn = pdata[8]; ft = pdata[9];
                      vd[0] = pdata[10]; vd[1] = pdata[11]; vd[2] = pdata[12];
                      // hd = tunnel_dir x vd (tunnel pdata[13..15])
                      hd[0] = pdata[14]*vd[2] - pdata[15]*vd[1];
                      hd[1] = pdata[15]*vd[0] - pdata[13]*vd[2];
                      hd[2] = pdata[13]*vd[1] - pdata[14]*vd[0];
                    }
                    if ( !array_normalize( vd, 3 ) ) {
                      vd[0] = 0.; vd[1] = -1.; vd[2] = 0.;
                    }
                    if ( !array_normalize( hd, 3 ) ) {
                      pri( "Error: force_edge_projected tunnel/vertical directions degenerate." );
                      exit(TN_EXIT_STATUS);
                    }
                    // edge tangent from the outward normal (2D rotation)
                    tang[0] = -normal[1]; tang[1] = normal[0]; tang[2] = 0.;
                    {
                      double ndothd=0., ndotvd=0., tdothd=0., tdotvd=0.;
                      long int kk;
                      for ( kk=0; kk<3; kk++ ) {
                        ndothd += hd[kk]*normal[kk];
                        ndotvd += vd[kk]*normal[kk];
                        tdothd += hd[kk]*tang[kk];
                        tdotvd += vd[kk]*tang[kk];
                      }
                      n2 = ph*ndothd*ndothd + pv*ndotvd*ndotvd;
                      t2 = ph*ndothd*tdothd + pv*ndotvd*tdotvd;
                    }
                    for ( idim=0; idim<ndim; idim++ ) {
                      ipuknwn = vel_indx/nder + idim;
                      tmp = factor * node_factor * load *
                        weight[inol_side] * area_size *
                        ( fn * n2 * normal[idim] + ft * t2 * tang[idim] );
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
