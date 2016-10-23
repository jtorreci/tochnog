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

#define EPS_DELETE 1.e-12
#define EPS_ELEMENT_DELETE_FACTOR 1.e-6

void delete_element( long int element, long int version )

{
  long int idat=0, data_class=0, swit=0;

  swit = set_swit(-1,-1,"delete_element");
  if ( swit ) pri( "In routine DELETE_ELEMENT" );

  for ( idat=0; idat<MDAT; idat++ ) {
    data_class = db_data_class( idat );
    if ( data_class==ELEMENT && db_version(idat,version) ) 
      db_delete_index( idat, element, version );
  }

  if ( swit ) pri( "Out routine DELETE_ELEMENT" );

}

void delete_node( long int inod, long int version )

{
  long int idat=0, data_class=0, swit=0;

  swit = set_swit(-1,-1,"delete_node");
  if ( swit ) pri( "In routine DELETE_NODE" );

  for ( idat=0; idat<MDAT; idat++ ) {
    data_class = db_data_class( idat );
    if ( data_class==NODE && db_version(idat,version) ) 
      db_delete_index( idat, inod, version );
  }

  if ( swit ) pri( "Out routine DELETE_NODE" );

}

void delete_geom( double time_current )

{
  long int element=0, max_element=0, inol=0, nnol=0, inod=0, length=0,
    all_in_geometry=0, in_geometry=0, element_group=0,
    max_node=0, any_element_deleted=0,
    node_boundary=-YES, swit=0, icontrol=0, 
    control_mesh_delete_geometry_movenodes=-YES, 
    length_control_mesh_delete_geometry_elementgroup=0,
    length_control_mesh_delete_geometry_element=0,
    name=0, ldum=0, idum[1], control_mesh_delete_geometry[2], 
    control_mesh_delete_geometry_elementgroup[DATA_ITEM_SIZE],
    control_mesh_delete_geometry_element[DATA_ITEM_SIZE], 
    nodes[MNOL], el[MNOL+1], *node_in_geometry=NULL;
  double element_volume=0., old_element_volume=0., 
    element_delete_factor=0., time_old=0., 
    time_new=0., rdum=0., f0=0., f1=0., 
    ddum[MDIM], coord[MDIM], diff_coord[MDIM], 
    node_start_refined[MDIM], projection[MDIM],
    control_mesh_delete_geometry_factor[2]; 

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  if ( db_active_index( CONTROL_MESH_DELETE_GEOMETRY, 
      icontrol, VERSION_NORMAL )  ) {

    swit = set_swit(-1,-1,"delete_geom");
    if ( swit ) pri( "In routine DELETE_GEOM" );

    db_max_index( NODE, max_node, VERSION_NORMAL, GET );
    length = 1+max_node;
    node_in_geometry = get_new_int( length );
    array_set( node_in_geometry, 0, length );

    db( CONTROL_MESH_DELETE_GEOMETRY, icontrol, control_mesh_delete_geometry, 
      ddum, ldum, VERSION_NORMAL, GET );

    control_mesh_delete_geometry_element[0] = -ALL;
    length_control_mesh_delete_geometry_element = 1;
    db( CONTROL_MESH_DELETE_GEOMETRY_ELEMENT, icontrol, 
      control_mesh_delete_geometry_element, ddum, 
      length_control_mesh_delete_geometry_element, VERSION_NORMAL,
      GET_IF_EXISTS );

    control_mesh_delete_geometry_factor[0] = 0.;
    control_mesh_delete_geometry_factor[1] = 1.;
    db( CONTROL_MESH_DELETE_GEOMETRY_FACTOR, icontrol, idum,
      control_mesh_delete_geometry_factor, ldum, VERSION_NORMAL,
      GET_IF_EXISTS );
    f0 = 1. - control_mesh_delete_geometry_factor[0];
    f1 = 1. - control_mesh_delete_geometry_factor[1];
    if ( swit ) {
      pri( "f0", f0 );
      pri( "f1", f1 );
    }

    control_mesh_delete_geometry_elementgroup[0] = -ALL;
    length_control_mesh_delete_geometry_elementgroup = 1;
    db( CONTROL_MESH_DELETE_GEOMETRY_ELEMENTGROUP, icontrol, 
      control_mesh_delete_geometry_elementgroup, ddum, 
      length_control_mesh_delete_geometry_elementgroup, VERSION_NORMAL,
      GET_IF_EXISTS );

    db( CONTROL_MESH_DELETE_GEOMETRY_MOVENODES, icontrol, 
      &control_mesh_delete_geometry_movenodes, ddum, ldum, VERSION_NORMAL,
      GET_IF_EXISTS );

    db_version_copy( VERSION_NORMAL, VERSION_TMP );
    db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
    if ( db_active_index( TIME_OLD, 0, VERSION_NORMAL ) ) {
      db( TIME_OLD, 0, idum, &time_old, ldum, VERSION_NORMAL, GET );
      db( TIME_NEW, 0, idum, &time_new, ldum, VERSION_NORMAL, GET );
      if ( scalar_dabs(time_new-time_old)>EPS_DELETE ) {
        element_delete_factor = f0 - ( f0 - f1 ) *
          (time_current-time_old)/(time_new-time_old);
      }
    }

      // determine which nodes are in the geometry
    for ( inod=0; inod<=max_node; inod++ ) {
      if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
         geometry( inod, ddum, control_mesh_delete_geometry, 
           in_geometry, rdum, ddum, rdum, ddum, NODE_START_REFINED, 
           CONTROL_MESH_DELETE_GEOMETRY, VERSION_NORMAL );
         if ( in_geometry ) node_in_geometry[inod] = 1;
      }
    }

      // delete elements which are totally in geometry
    for ( element=0; element<=max_element; element++ ) {
      if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) ) {
        db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
        name = el[0];
        element_group = 0;
        db( ELEMENT_GROUP, element, &element_group, ddum, 
          ldum, VERSION_NORMAL, GET_IF_EXISTS );
        if ( control_mesh_delete_geometry_element[0]==-ALL || 
             array_member(control_mesh_delete_geometry_element,name,
             length_control_mesh_delete_geometry_element,ldum) ) {
          if ( control_mesh_delete_geometry_elementgroup[0]==-ALL || 
               array_member(control_mesh_delete_geometry_elementgroup,element_group,
               length_control_mesh_delete_geometry_elementgroup,ldum) ) {   
            nnol = length - 1; array_move( &el[1], nodes, nnol );
            all_in_geometry = 1;
            for ( inol=0; inol<nnol; inol++ ) {
              inod = nodes[inol];
              if ( !node_in_geometry[inod] ) all_in_geometry = 0;
            }
            if ( all_in_geometry ) {
              if ( element_delete_factor>EPS_ELEMENT_DELETE_FACTOR ) {
                length = 1;
                db( ELEMENT_DELETE_FACTOR, element, idum, &element_delete_factor, 
                  length, VERSION_NORMAL, PUT );
              }
              else {
                delete_element( element, VERSION_NORMAL );
                any_element_deleted = 1;
              }
            }
          }
        }
      }
    }

    if ( any_element_deleted ) {
         // move remaining nodes in geometry to edges
      if ( control_mesh_delete_geometry_movenodes==-YES ) {
        for ( inod=0; inod<=max_node; inod++ ) {
          if ( node_in_geometry[inod] ) {
            geometry( inod, ddum, control_mesh_delete_geometry, 
              in_geometry, rdum, ddum, rdum, projection, 
              NODE_START_REFINED, PROJECT_ON_EDGE, VERSION_NORMAL );
            db( NODE_START_REFINED, inod, idum, node_start_refined, 
              ldum, VERSION_NORMAL, GET );
            db( NODE_START_REFINED, inod, idum, projection, 
              ndim, VERSION_NORMAL, PUT );
            db( NODE, inod, idum, coord, ldum, VERSION_NORMAL, GET );
            array_subtract( projection, node_start_refined, diff_coord, ndim );
            array_add( coord, diff_coord, coord, ndim );
            db( NODE, inod, idum, coord, ndim, VERSION_NORMAL, PUT );
            length=1; db( NODE_BOUNDARY, inod, &node_boundary, ddum, 
              length, VERSION_NORMAL, PUT );
          }
        }
         // delete too small elements (collapsed because of moving nodes to edge)
        for ( element=0; element<=max_element; element++ ) {
          if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) ) {
            db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
            name = el[0];
            if ( name==-BAR2 || name==-TRIA3 || name==-TET4 ) {
              nnol = length - 1; array_move( &el[1], nodes, nnol );
              element_volume_set( name, nodes, VERSION_NORMAL, element_volume );
              element_volume_set( name, nodes, VERSION_TMP, old_element_volume );
              if ( element_volume<EPS_VOLUME*old_element_volume ) 
                delete_element( element, VERSION_NORMAL );
            }
          }
        }
      }
      mesh_has_changed( VERSION_NORMAL );
    }
    db_version_delete( VERSION_TMP );

    delete[] node_in_geometry;

    if ( swit ) pri( "Out routine DELETE_GEOM" );

  }

}

void mesh_delete_small( long int version )

{

  long int icontrol=0, element=0, max_element=0, swit=0, 
    any_element_deleted=0, ldum=0, idum[1];
  double small=0., element_volume=0., ddum[1];

  swit = set_swit(-1,-1,"mesh_delete_small");
  if ( swit ) pri( "In routine MESH_DELETE_SMALL" );

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  if ( db_active_index( CONTROL_MESH_DELETE_SMALL, icontrol, VERSION_NORMAL ) ) {
    db( CONTROL_MESH_DELETE_SMALL, icontrol, idum, &small, ldum, VERSION_NORMAL, GET );
    db_max_index( ELEMENT, max_element, version, GET );
    for ( element=0; element<=max_element; element++ ) {
      if ( db_active_index( ELEMENT_VOLUME, element, version ) ) {
        db( ELEMENT_VOLUME, element, idum, &element_volume, ldum, version, GET );
        if ( element_volume<small ) {
          if ( swit ) pri( "Deleting element ", element );
          any_element_deleted = 1;
          delete_element( element, version );
        }
      }
    }
  }

  if ( any_element_deleted ) mesh_has_changed( version );

  if ( swit ) pri( "Out routine MESH_DELETE_SMALL" );
}
