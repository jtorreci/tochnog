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

#define EPS_ERROR 1.e-10

void error( long int task )

   /*  Calculate error estimates.

       POST_ERROR_MESH2 is on the oldest mesh.
       POST_ERROR_MESH1 is on the previous mesh.
       value is on the current mesh.
   */
{
  long int data_item_name=0, data_item_index=0, number=0, ipost=0,
    max_post_error_item=0, length=0, ldum=0, 
    idum[1], post_error_item[3], *dof_label=NULL;
  double post_error_result=0., value=0., post_error_mesh1=0., 
    post_error_mesh2=0., a=0., b=0., estimated_value=0., 
    ddum[1], *ddat=NULL;

  dof_label = get_new_int(MUKNWN);
  ddat = get_new_dbl(DATA_ITEM_SIZE);
  db_max_index( POST_ERROR_ITEM, max_post_error_item, VERSION_NORMAL, GET );
  db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );

  for ( ipost=0; ipost<=max_post_error_item; ipost++ ) {
    if ( db_active_index( POST_ERROR_ITEM, ipost, VERSION_NORMAL ) ) {
      db( POST_ERROR_ITEM, ipost, post_error_item, ddum, ldum, 
        VERSION_NORMAL, GET );
      data_item_name = post_error_item[0];
      data_item_index = post_error_item[1];
      if ( db_type(data_item_name)!=DOUBLE_PRECISION ) 
        db_error( POST_ERROR_ITEM, ipost );
      if ( db( data_item_name, data_item_index, idum, ddat, length, 
          VERSION_NORMAL, GET_IF_EXISTS ) ) {
        if ( post_error_item[2]<0 ) {
          array_member(dof_label,post_error_item[2],nuknwn,number);
          if ( length==npuknwn ) number /= nder;
        }
        else
          number = post_error_item[2];
        if ( number<0 || number>length ) 
          db_error( POST_ERROR_ITEM, ipost );
        else
          value = ddat[number];
        if ( task==GET && 
             db_active_index( POST_ERROR_MESH1, ipost, VERSION_NORMAL ) &&
             db_active_index( POST_ERROR_MESH2, ipost, VERSION_NORMAL ) ) {
          db( POST_ERROR_MESH1, ipost, idum, &post_error_mesh1, length, 
            VERSION_NORMAL, GET );
          db( POST_ERROR_MESH2, ipost, idum, &post_error_mesh2, length, 
            VERSION_NORMAL, GET );
          if ( scalar_dabs(post_error_mesh2)>EPS_ERROR ) {
            a = scalar_dabs( post_error_mesh2 - post_error_mesh1 );
            b = scalar_dabs( post_error_mesh1 - value );
            if ( scalar_dabs(a)>EPS_ERROR )
              estimated_value = value - (b/a)*( value - post_error_mesh1 );
            else
              estimated_value = value - 0.5 * ( value - post_error_mesh1 );
            if ( scalar_dabs(estimated_value)>EPS_ERROR ) {
              post_error_result = 
                100. * scalar_dabs(estimated_value-value)/estimated_value;
              db( POST_ERROR_RESULT, ipost, idum, &post_error_result, length, 
                VERSION_NORMAL, PUT );
            }
          }
        }
        else if ( task==PUT ) {
          if ( db_active_index( POST_ERROR_MESH1, ipost, VERSION_NORMAL ) ) {
            db( POST_ERROR_MESH1, ipost, idum, &post_error_mesh1, length, 
              VERSION_NORMAL, GET );
            db( POST_ERROR_MESH2, ipost, idum, &post_error_mesh1, length, 
              VERSION_NORMAL, PUT );
            db( POST_ERROR_MESH1, ipost, idum, &value, length, 
              VERSION_NORMAL, PUT );
          }
          else {
            length = 1;
            db( POST_ERROR_MESH1, ipost, idum, &value, length, VERSION_NORMAL, PUT );
          }
        }
      }
    }
  }
  delete[] dof_label;
  delete[] ddat;
}
