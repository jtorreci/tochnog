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

#define EPS_DELETE_FACTOR 1.e-6

void failure( double time_current )

{
  long int element=0, max_element=0, any_element_deleted=0, length=0, 
    inol=0, idim=0, jdim=0, indx=0, nnol=0, inod=0, swit=0, swit_element=0, 
    element_group=0, max_dependency_item=0, any_failure=0, idep=0, 
    print_failure=-NO, ldum=0, idum[1], el[MNOL+1], nodes[MNOL], *dependency_item=NULL;
  double tmp=0., threshold=0., delete_time=0., element_delete_factor=0., 
    ddum[1], element_delete_times[2], dval[2], 
    strain[MDIM*MDIM], workval[MDIM], workvec[MDIM*MDIM],
    average_node_dof[MUKNWN], *node_dof=NULL;

  if ( materi_strain_plasti || materi_strain_total || materi_damage || materi_void_fraction ) {
    
    if ( db_partialname_any("group_materi_failure") ) any_failure = 1;

    db_max_index( DEPENDENCY_ITEM, max_dependency_item, VERSION_NORMAL, GET );
    for ( idep=0; idep<=max_dependency_item; idep++ ) {
      if ( db_active_index( DEPENDENCY_ITEM, idep, VERSION_NORMAL ) ) {
        dependency_item = db_int( DEPENDENCY_ITEM, idep, VERSION_NORMAL );
        if ( labs(dependency_item[0])==GROUP_MATERI_FAILURE_PLASTI_KAPPA ||
             labs(dependency_item[0])==GROUP_MATERI_FAILURE_VOIDFRACTION ||
             labs(dependency_item[0])==GROUP_MATERI_FAILURE_CRUCHING ||
             labs(dependency_item[0])==GROUP_MATERI_FAILURE_RUPTURE ||
             labs(dependency_item[0])==GROUP_MATERI_FAILURE_VOIDFRACTION ||
             labs(dependency_item[0])==GROUP_MATERI_FAILURE_DAMAGE ) 
          any_failure = 1;
      }
    }

    if ( any_failure ) {
      swit = set_swit(-1,-1,"failure");
      if ( swit ) pri( "In routine FAILURE" );

      db( PRINT_FAILURE, 0, &print_failure, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );

        // loop over elements
      db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
      for ( element=0; element<=max_element; element++ ) {
        if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) ) {
          swit_element = swit; swit = swit && set_swit(element,-1,"");
          db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
          nnol = length - 1; array_move( &el[1], nodes, nnol );
          array_set( average_node_dof, 0., nuknwn );
          for ( inol=0; inol<nnol; inol++ ) {
            inod = nodes[inol];
            node_dof = db_dbl( NODE_DOF, inod, VERSION_NORMAL );
            array_add( average_node_dof, node_dof, average_node_dof, nuknwn );
          }
          array_multiply( average_node_dof, average_node_dof, 1./((double)nnol), nuknwn );
          if ( db_active_index( ELEMENT_GROUP, element, VERSION_NORMAL ) )
            db( ELEMENT_GROUP, element, &element_group, ddum, ldum, VERSION_NORMAL, GET );
          else
            element_group = 0;
          if ( get_group_data( GROUP_MATERI_FAILURE_PLASTI_KAPPA, 
              element_group, element,
              average_node_dof, dval, ldum, GET_IF_EXISTS ) ) {
            tmp = scalar_dabs( average_node_dof[kap_indx] );
          }
          if ( get_group_data( GROUP_MATERI_FAILURE_DAMAGE, element_group, 
              element, average_node_dof, dval, ldum, GET_IF_EXISTS ) ) {
            tmp = scalar_dabs( average_node_dof[dam_indx] );
          }
          if ( get_group_data( GROUP_MATERI_FAILURE_VOIDFRACTION, element_group, 
              element, average_node_dof, dval, ldum, GET_IF_EXISTS ) ) {
            tmp = scalar_dabs( average_node_dof[void_indx] );
          }
          if ( get_group_data( GROUP_MATERI_FAILURE_CRUCHING, element_group, 
              element, average_node_dof, dval, ldum, GET_IF_EXISTS ) ) {
            for ( idim=0; idim<MDIM; idim++ ) {
              for ( jdim=0; jdim<MDIM; jdim++ ) {
                indx = idim*MDIM + jdim;
                strain[indx] = average_node_dof[ept_indx+stress_indx(idim,jdim)*nder];
              }
            }
            matrix_jacobi( strain, MDIM, workval, workvec, idum );
            sort( workval, workvec );
            tmp = workval[2];
            if ( tmp>0. ) tmp = 0.;
          }
          if ( get_group_data( GROUP_MATERI_FAILURE_RUPTURE, element_group, 
              element, average_node_dof, dval, ldum, GET_IF_EXISTS ) ) {
            for ( idim=0; idim<MDIM; idim++ ) {
              for ( jdim=0; jdim<MDIM; jdim++ ) {
                indx = idim*MDIM + jdim;
                strain[indx] = average_node_dof[ept_indx+stress_indx(idim,jdim)*nder];
              }
            }
            matrix_jacobi( strain, MDIM, workval, workvec, idum );
            sort( workval, workvec );
            tmp = workval[0];
            if ( tmp<0. ) tmp = 0.;
          }
          threshold = dval[0];
          delete_time = dval[1];
          if ( tmp>threshold && 
               !db_active_index( ELEMENT_DELETE_TIMES, element, VERSION_NORMAL ) ) {
            if ( print_failure==-YES ) pri( "Failure for element ", element );
            element_delete_times[0] = time_current;
            element_delete_times[1] = time_current + delete_time;
            length = 2; db( ELEMENT_DELETE_TIMES, element, idum, element_delete_times, 
              length, VERSION_NORMAL, PUT );
          }
          if ( db_active_index( ELEMENT_DELETE_TIMES, element, VERSION_NORMAL ) ) {
            db( ELEMENT_DELETE_TIMES, element, idum, element_delete_times, 
              ldum, VERSION_NORMAL, GET );
            if ( time_current>element_delete_times[1] ) {
              delete_element( element, VERSION_NORMAL );
              any_element_deleted = 1;
            }
            else {
              element_delete_factor = 1. - ( time_current - element_delete_times[0] ) / 
                ( element_delete_times[1] - element_delete_times[0] );
              if ( element_delete_factor<EPS_DELETE_FACTOR ) 
                element_delete_factor = EPS_DELETE_FACTOR;
              length = 1; db( ELEMENT_DELETE_FACTOR, element, idum, &element_delete_factor, 
                length, VERSION_NORMAL, PUT );
            }
          }
          swit = swit_element;
        }
      }
      if ( any_element_deleted ) mesh_has_changed( VERSION_NORMAL );

      if ( swit ) pri( "Out routine FAILURE" );
    }

  }

}
