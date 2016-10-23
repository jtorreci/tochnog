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

void create_element( long int old_element, long int new_element,
  long int el[], long int length_el, long int version_from, long int version_to )

{
  long int idat=0, data_class=0, length=0, idum[1], *ival=NULL,
  	mnolnuknwn=npointmax*nuknwn, zero=0, one=1,
        length_nei=1+npointmax*ndim+npointmax+2;
  double ddum[1], *dval=NULL, *tmp_element_dof=NULL, *dworknei=NULL;

  ival = get_new_int(DATA_ITEM_SIZE);
  dval = get_new_dbl(DATA_ITEM_SIZE);
  // added for options_element_dof
  tmp_element_dof = get_new_dbl(mnolnuknwn);
  dworknei = get_new_dbl(length_nei);
  array_set(tmp_element_dof, 0, mnolnuknwn);
  array_set(dworknei, 0, length_nei);
  dworknei[length_nei-1]=1;

  for ( idat=0; idat<MDAT; idat++ ) {
    data_class = db_data_class( idat );
    if ( data_class==ELEMENT && db_version(idat,version_from) &&
         db_version(idat,version_to) ) {
      if ( idat==ELEMENT )
        db( ELEMENT, new_element, el, ddum, length_el, version_to, PUT );
      else if ( db_active_index( idat, old_element, version_from ) 
      	&& idat!=ELEMENT_DOF && idat!=ELEMENT_DOF_INITIALISED
	&& idat!=NONLOCAL_ELEMENT_INFO) {
        if ( db_type(idat)==DOUBLE_PRECISION ) {
          db( idat, old_element, idum, dval, length, version_from, GET );
          db( idat, new_element, idum, dval, length, version_to, PUT );
        }
        else {
          db( idat, old_element, ival, ddum, length, version_from, GET );
          db( idat, new_element, ival, ddum, length, version_to, PUT );
        }
      }
    }
    if ( data_class==ELEMENT && db_version(idat,version_to)) {
      //element_dof will be re--initialised anyway->just version_to
      // added for options_element_dof
      if ( idat==ELEMENT_DOF ) 
        db( idat, new_element, idum, tmp_element_dof, mnolnuknwn, version_to, PUT );
      else if ( idat==ELEMENT_DOF_INITIALISED ) 
        db( idat, new_element, &zero, ddum, one, version_to, PUT );
    }

    	//Problem with renumbering in print_gid. Updates NEI only if new mesh is generated
    if( version_from == VERSION_MACRO || version_from == VERSION_NORMAL || 
    	version_from == VERSION_NEW_MESH_TMP || version_from == VERSION_NEW_MESH_GENERATED ) {
         if ( data_class==ELEMENT && idat==NONLOCAL_ELEMENT_INFO ) {
           db( idat, new_element, idum, dworknei, length_nei, VERSION_NORMAL, PUT );		
	   nonlocal_first_set=0;
	 }
     }
  }
  delete[] ival;
  delete[] dval;
  delete[] tmp_element_dof;
  delete[] dworknei;
}

void create_node( long int old_nodes[], long int nnod, long int tmp_node_number,
  double tmp_node[], double tmp_node_dof[], double tmp_node_start_refined[],
  double tmp_node_dof_start_refined[], long int version_from, long int version_to )

{
  long int inod=0, idat=0, data_class=0, length=0, 
    present_in_all_old_nodes=0, idum[1], *ival=NULL;
  double ddum[1], *dval=NULL;

  for ( idat=0; idat<MDAT; idat++ ) {
    data_class = db_data_class( idat );
    if ( data_class==NODE && 
         db_version(idat,version_from) &&
         db_version(idat,version_to) ) {
      present_in_all_old_nodes = 1;
      for ( inod=0; inod<nnod; inod++ ) if ( !db_active_index( idat, 
        old_nodes[inod], version_from ) ) present_in_all_old_nodes = 0;
      if ( present_in_all_old_nodes ) {
        length = db_len( idat, old_nodes[0], version_from );
        if      ( idat==NODE )
          db( idat, tmp_node_number, idum, tmp_node, length, version_to, PUT );
        else if ( idat==NODE_DOF )
          db( idat, tmp_node_number, idum, tmp_node_dof, length, version_to, PUT );
        else if ( idat==NODE_START_REFINED )
          db( idat, tmp_node_number, idum, tmp_node_start_refined, 
            length, version_to, PUT );
        else if ( idat==NODE_DOF_START_REFINED )
          db( idat, tmp_node_number, idum, tmp_node_dof_start_refined, 
            length, version_to, PUT );
        else if ( db_type(idat)==DOUBLE_PRECISION ) {
          dval = db_dbl( idat, old_nodes[0], version_from );
          db( idat, tmp_node_number, idum, dval, length, version_to, PUT );
        }
        else {
          ival = db_int( idat, old_nodes[0], version_from );
          db( idat, tmp_node_number, ival, ddum, length, version_to, PUT );
        }
      }
    }
  }

}
