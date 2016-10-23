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

void renumbering( long int version, long int fill_old_numbers, long int lowest_element, 
  long int lowest_node, long int old_node_numbers[], 
  long int old_element_numbers[] )

{
  long int element=0, max_element=0, new_element_number=0,
    inol=0, inod=0, max_node=0, tmp_node_number=0, 
    swit=0, length=0, *ival=NULL, *new_node_numbers=NULL;
  double ddum[1], *tmp_node=NULL, *tmp_node_dof=NULL, 
    *tmp_node_start_refined=NULL, *tmp_node_dof_start_refined=NULL;

  swit = set_swit(-1,-1,"renumbering");
  if ( swit ) pri( "In routine RENUMBERING" );

  db_max_index( ELEMENT, max_element, version, GET );
  db_max_index( NODE, max_node, version, GET );

  ival = get_new_int( DATA_ITEM_SIZE );
  new_node_numbers = get_new_int( max_node+1 );
  array_set( new_node_numbers, -1, max_node+1 );

  tmp_node_number = lowest_node;
  for ( inod=0; inod<=max_node; inod++ ) {
    if ( db_active_index( NODE, inod, version ) ) {
      ival[0] = inod;
      tmp_node = db_dbl( NODE, inod, version );
      tmp_node_start_refined = db_dbl( NODE_START_REFINED, inod, version );
      if ( nuknwn>0 ) {
        tmp_node_dof = db_dbl( NODE_DOF, inod, version );
        tmp_node_dof_start_refined = db_dbl( NODE_DOF_START_REFINED, inod, version );
      }
      else {
        tmp_node_dof = ddum;
        tmp_node_dof_start_refined = ddum;
      }
      create_node( ival, 1, tmp_node_number, tmp_node, tmp_node_dof,
        tmp_node_start_refined, tmp_node_dof_start_refined,
        version, VERSION_TMP );
      new_node_numbers[inod] = tmp_node_number;
      if ( fill_old_numbers==YES ) old_node_numbers[tmp_node_number]=inod;
      tmp_node_number++;
    }   
  }

  new_element_number = lowest_element;
  for ( element=0; element<=max_element; element++ ) {
    if ( db_active_index( ELEMENT, element, version ) ) {
      db( ELEMENT, element, ival, ddum, length, version, GET );
      for ( inol=0; inol<length-1; inol++ ) {
        inod = ival[1+inol];
        tmp_node_number = new_node_numbers[inod];
        if ( tmp_node_number<0 ) {
          pri( "Error: non-existing node detected for element ", element );
          exit(TN_EXIT_STATUS);
        }
        ival[1+inol] = tmp_node_number;
      }
      create_element( element, new_element_number, ival, length,
        version, VERSION_TMP );
      if ( fill_old_numbers==YES ) old_element_numbers[new_element_number]=element;
      new_element_number++;
    }
  }

  db_version_copy( VERSION_TMP, version );
  db_version_delete( VERSION_TMP );

  if ( version==VERSION_NORMAL ) mesh_has_changed( version );

  delete[] ival;
  delete[] new_node_numbers;

  if ( swit ) pri( "Out routine RENUMBERING" );

}

void renumbering_check( long int idat )

{
  long int type=NONE, max=0;

  if      ( db_max_index( CONTROL_MESH_REFINE_LOCALLY, max, VERSION_NORMAL, GET ) >= 0 ) 
    type = CONTROL_MESH_REFINE_LOCALLY;
  else if ( db_max_index( CONTROL_MESH_REFINE_GLOBALLY, max, VERSION_NORMAL, GET ) >= 0 ) 
    type = CONTROL_MESH_REFINE_GLOBALLY;
  else if ( db_max_index( CONTROL_MESH_RENUMBER, max, VERSION_NORMAL, GET ) >= 0 ) 
    type = CONTROL_MESH_RENUMBER;

  if ( type!=NONE ) {
    cout << "\nError: direct node numbers in data item " << db_name(idat);
    cout << "\ncannot be combined with " << db_name(type) << ".\n";
    cout << "\nUse a geometrical entity instead.\n";
    exit(TN_EXIT_STATUS);
  }

}
