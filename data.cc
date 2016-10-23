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

#define EPS_TIME 1.e-12

void data( long int task, double dtime, double time_current )

{
  long int idat=0, in=0, iv=0, index=0, range_length=0, icontrol=0, length=0, 
    swit=0, max_index=0, inod=0, max_node=0, found=0, 
    ichange=0, max_change=0, idim=0, operat=0, ldum=0, 
    data_item_name=0, data_item_index=0, data_item_number=0,
    change_dataitem_time_discrete=-NO, change_dataitem_time_user=0,
    idum[1], change_dataitem[4], *dof_label=NULL, *integer_range=NULL, 
    *data_delete=NULL, *data_put=NULL;
  double rdum=0., val=0., ddum[MDIM], *change_dataitem_time=NULL, 
    *dval=NULL, *coord=NULL;

  dof_label = get_new_int(MUKNWN);
  integer_range = get_new_int(MRANGE);
  change_dataitem_time = get_new_dbl(DATA_ITEM_SIZE);
  dval = get_new_dbl(DATA_ITEM_SIZE);
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  db_max_index( NODE, max_node, VERSION_NORMAL, GET );

  if ( db_active_index( CONTROL_DATA_DELETE, icontrol, VERSION_NORMAL )  ) {
    swit = set_swit(-1,-1,"data");
    if ( swit ) pri( "In routine DATA" );
    length = db_len( CONTROL_DATA_DELETE, icontrol, VERSION_NORMAL );
    data_delete = db_int( CONTROL_DATA_DELETE, icontrol, VERSION_NORMAL );
    idat = data_delete[0];
    if ( data_delete[1]>=0 ) {
      index = data_delete[1];
      if      ( idat==-ELEMENT )
        delete_element( index, VERSION_NORMAL );
      else if ( idat==-NODE )
        delete_node( index, VERSION_NORMAL );
      else
        db_delete_index( idat, index, VERSION_NORMAL );
    }
    else if ( data_delete[1]==-RA ) {
      range_expand( &data_delete[1], integer_range, length, range_length );
      for ( in=0; in<range_length; in++ ) {
        index = integer_range[in];
        if      ( idat==-ELEMENT )
          delete_element( index, VERSION_NORMAL );
        else if ( idat==-NODE )
          delete_node( index, VERSION_NORMAL );
        else
          db_delete_index( idat, index, VERSION_NORMAL );
      }
    }
    else if ( data_delete[1]==-ALL ) {
      db_max_index( idat, max_index, VERSION_NORMAL, GET );
      for ( index=0; index<max_index; index++ ) {
        if      ( idat==-ELEMENT )
          delete_element( index, VERSION_NORMAL );
        else if ( idat==-NODE )
          delete_node( index, VERSION_NORMAL );
        else
          db_delete_index( idat, index, VERSION_NORMAL );
      }
    }
    else if ( db_data_class(data_delete[1])==GEOMETRY &&
              db_data_class(data_delete[0])==NODE ) {
      for ( inod=0; inod<=max_node; inod++ ) {
        if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
          geometry( inod, ddum, &data_delete[1], found, rdum, ddum, rdum,
            ddum, NODE_START_REFINED, PROJECT_EXACT, VERSION_NORMAL );
          if ( found ) {
            if ( idat==-NODE )
              delete_node( inod, VERSION_NORMAL );
            else
              db_delete_index( idat, inod, VERSION_NORMAL );
          }
        }
      }
    }
    else
      db_error( CONTROL_DATA_DELETE, icontrol );
    if ( swit ) pri( "Out routine DATA" );
  }

  if ( db_active_index( CONTROL_DATA_PUT, icontrol, VERSION_NORMAL )  ) {
    swit = set_swit(-1,-1,"data");
    if ( swit ) pri( "In routine DATA" );
    length = db_len( CONTROL_DATA_PUT, icontrol, VERSION_NORMAL );
    data_put = db_int( CONTROL_DATA_PUT, icontrol, VERSION_NORMAL );
    idat = data_put[0];
    if ( db_data_class(idat)==TENDON ) {
      pri( "Error: CONTROL_DATA_PUT cannot be used for tendon data. " );
      exit(TN_EXIT_STATUS);
    }
    if ( data_put[1]>=0 ) {
      index = data_put[1];
      if ( db_type(idat)==INTEGER ) {
        length = db_len( CONTROL_DATA_PUT_INTEGER, icontrol, VERSION_NORMAL );
        db( idat, index, db_int(CONTROL_DATA_PUT_INTEGER, icontrol, 
          VERSION_NORMAL), ddum, length, VERSION_NORMAL, PUT );
      }
      else if ( db_data_class(idat)==NODE && 
          db_active_index( CONTROL_DATA_PUT_DOUBLE_NODE, icontrol, VERSION_NORMAL ) ) {
        length = db_len( idat, index, VERSION_NORMAL );
        ldum = (1+ndim)*length;
        db( CONTROL_DATA_PUT_DOUBLE_NODE, icontrol, idum, dval, 
          ldum, VERSION_NORMAL, GET_AND_CHECK );
        coord = db_dbl( NODE, index, VERSION_NORMAL );
        for ( idim=0; idim<ndim; idim++ ) {
          for ( iv=0; iv<length; iv++ ) {
            dval[iv] += coord[idim] * dval[(1+idim)*length+iv];
          }
        }
        if ( length > db_data_length(idat) ) {
          db_error( CONTROL_DATA_PUT_DOUBLE_NODE, icontrol );
        }
        db( idat, index, idum, dval, length, VERSION_NORMAL, PUT );
        if ( labs(idat)==NODE_DOF && time_current==0. )
          db( NODE_DOF_START_REFINED, index, idum, dval, length, VERSION_NORMAL, PUT );
      }
      else {
        length = db_len( CONTROL_DATA_PUT_DOUBLE, icontrol, VERSION_NORMAL );
        if ( length > db_data_length(idat) ) {
          db_error( CONTROL_DATA_PUT_DOUBLE, icontrol );
        }
        db( idat, index, idum, db_dbl(CONTROL_DATA_PUT_DOUBLE, icontrol, 
          VERSION_NORMAL), length, VERSION_NORMAL, PUT );
        if ( labs(idat)==NODE_DOF && time_current==0. )
          db( NODE_DOF_START_REFINED, index, idum, db_dbl(CONTROL_DATA_PUT_DOUBLE, icontrol,
            VERSION_NORMAL), length, VERSION_NORMAL, PUT );
      }
    }
    else if ( data_put[1]==-RA ) {
      range_expand( &data_put[1], integer_range, length, range_length );
      for ( in=0; in<range_length; in++ ) {
        index = integer_range[in];
        if ( db_type(idat)==INTEGER ) {
          length = db_len( CONTROL_DATA_PUT_INTEGER, icontrol, VERSION_NORMAL );
          db( idat, index, db_int(CONTROL_DATA_PUT_INTEGER, icontrol, 
            VERSION_NORMAL), ddum, length, VERSION_NORMAL, PUT );
        }
        else if ( db_data_class(idat)==NODE && 
            db_active_index( CONTROL_DATA_PUT_DOUBLE_NODE, icontrol, VERSION_NORMAL ) ) {
          length = db_len( idat, index, VERSION_NORMAL );
          ldum = (1+ndim)*length;
          db( CONTROL_DATA_PUT_DOUBLE_NODE, icontrol, idum, dval, 
            ldum, VERSION_NORMAL, GET_AND_CHECK );
          coord = db_dbl( NODE, index, VERSION_NORMAL );
          for ( idim=0; idim<ndim; idim++ ) {
            for ( iv=0; iv<length; iv++ ) {
              dval[iv] += coord[idim] * dval[(1+idim)*length+iv];
            }
          }
          db( idat, index, idum, dval, length, VERSION_NORMAL, PUT );
          if ( labs(idat)==NODE_DOF && time_current==0. )
            db( NODE_DOF_START_REFINED, index, idum, dval, length, VERSION_NORMAL, PUT );
        }
        else {
          length = db_len( CONTROL_DATA_PUT_DOUBLE, icontrol, VERSION_NORMAL );
          db( idat, index, idum, db_dbl(CONTROL_DATA_PUT_DOUBLE, icontrol, 
            VERSION_NORMAL), length, VERSION_NORMAL, PUT );
          if ( labs(idat)==NODE_DOF && time_current==0. )
            db( NODE_DOF_START_REFINED, index, idum, db_dbl(CONTROL_DATA_PUT_DOUBLE, icontrol,
              VERSION_NORMAL), length, VERSION_NORMAL, PUT );
        }
      }
    }
    else if ( data_put[1]==-ALL && db_data_class(idat)==NODE ) {
      db_max_index( NODE, max_index, VERSION_NORMAL, GET );
      for ( index=0; index<=max_index; index++ ) {
        if ( db_active_index( NODE, index, VERSION_NORMAL ) ) {
          if ( db_type(idat)==INTEGER ) {
            length = db_len( CONTROL_DATA_PUT_INTEGER, icontrol, VERSION_NORMAL );
            db( idat, index, db_int(CONTROL_DATA_PUT_INTEGER, icontrol, 
              VERSION_NORMAL), ddum, length, VERSION_NORMAL, PUT );
          }
          else if ( db_data_class(idat)==NODE && 
              db_active_index( CONTROL_DATA_PUT_DOUBLE_NODE, 
                icontrol, VERSION_NORMAL ) ) {
            length = db_len( idat, index, VERSION_NORMAL );
            ldum = (1+ndim)*length;
            db( CONTROL_DATA_PUT_DOUBLE_NODE, icontrol, idum, dval, 
              ldum, VERSION_NORMAL, GET_AND_CHECK );
            coord = db_dbl( NODE, index, VERSION_NORMAL );
            for ( idim=0; idim<ndim; idim++ ) {
              for ( iv=0; iv<length; iv++ ) {
                dval[iv] += coord[idim] * dval[(1+idim)*length+iv];
              }
            }
            if ( length > db_data_length(idat) ) {
              db_error( CONTROL_DATA_PUT_DOUBLE_NODE, icontrol );
            }
            db( idat, index, idum, dval, length, VERSION_NORMAL, PUT );
            if ( labs(idat)==NODE_DOF && time_current==0. ) {
              db( NODE_DOF_START_REFINED, index, idum, dval, length, VERSION_NORMAL, PUT );
            }
          }
          else {
            length = db_len( CONTROL_DATA_PUT_DOUBLE, icontrol, VERSION_NORMAL );
            if ( length > db_data_length(idat) ) {
              db_error( CONTROL_DATA_PUT_DOUBLE_NODE, icontrol );
            }
            db( idat, index, idum, db_dbl(CONTROL_DATA_PUT_DOUBLE, icontrol, 
              VERSION_NORMAL), length, VERSION_NORMAL, PUT );
            if ( labs(idat)==NODE_DOF && time_current==0. )
              db( NODE_DOF_START_REFINED, index, idum, db_dbl(CONTROL_DATA_PUT_DOUBLE, icontrol,
                VERSION_NORMAL), length, VERSION_NORMAL, PUT );
          }
        }
      }
    }
    else if ( db_data_class(data_put[1])==GEOMETRY &&
              db_data_class(data_put[0])==NODE ) {
      for ( inod=0; inod<=max_node; inod++ ) {
        if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
          geometry( inod, ddum, &data_put[1], found, rdum, ddum, rdum,
            ddum, NODE_START_REFINED, PROJECT_EXACT, VERSION_NORMAL );
          if ( found ) {
            if ( db_type(idat)==INTEGER ) {
              length = db_len( CONTROL_DATA_PUT_INTEGER, icontrol, VERSION_NORMAL );
              db( idat, inod, db_int(CONTROL_DATA_PUT_INTEGER, icontrol, 
                VERSION_NORMAL), ddum, length, VERSION_NORMAL, PUT );
            }
            else if ( db_data_class(idat)==NODE && 
                db_active_index( CONTROL_DATA_PUT_DOUBLE_NODE, 
                  icontrol, VERSION_NORMAL ) ) {
              length = db_len( idat, inod, VERSION_NORMAL );
              ldum = (1+ndim)*length;
              db( CONTROL_DATA_PUT_DOUBLE_NODE, icontrol, idum, dval, 
                ldum, VERSION_NORMAL, GET_AND_CHECK );
              coord = db_dbl( NODE, inod, VERSION_NORMAL );
              for ( idim=0; idim<ndim; idim++ ) {
                for ( iv=0; iv<length; iv++ ) {
                  dval[iv] += coord[idim] * dval[(1+idim)*length+iv];
                }
              }
              db( idat, inod, idum, dval, length, VERSION_NORMAL, PUT );
              if ( labs(idat)==NODE_DOF && time_current==0. )
                db( NODE_DOF_START_REFINED, inod, idum, dval, length, VERSION_NORMAL, PUT );
            }
            else {
              length = db_len( CONTROL_DATA_PUT_DOUBLE, icontrol, VERSION_NORMAL );
              db( idat, inod, idum, db_dbl(CONTROL_DATA_PUT_DOUBLE, icontrol, 
                VERSION_NORMAL), length, VERSION_NORMAL, PUT );
            }
          }
        }
      }
    }
    else
      db_error( CONTROL_DATA_PUT, icontrol );
    if ( swit ) pri( "Out routine DATA" );
  }
  else {
    if ( db_active_index( CONTROL_DATA_PUT_INTEGER, icontrol, VERSION_NORMAL ) )
      db_error( CONTROL_DATA_PUT, icontrol );
    if ( db_active_index( CONTROL_DATA_PUT_DOUBLE, icontrol, VERSION_NORMAL ) )
      db_error( CONTROL_DATA_PUT, icontrol );
    if ( db_active_index( CONTROL_DATA_PUT_DOUBLE_NODE, icontrol, VERSION_NORMAL ) )
      db_error( CONTROL_DATA_PUT, icontrol );
  }

  if ( db_active_index( CONTROL_DATA_INITELDOF_GEOMETRY, icontrol, VERSION_NORMAL )  ) {

    long int in_geometry=0, max_element=0, *node_in_geometry=NULL,
    	control_data_initeldof_geometry[2], nodes[MNOL], el[MNOL+1], zero=0, one=1,
	length_nodes = 1+max_node;
    node_in_geometry = get_new_int( length_nodes );
    array_set(el, 0, MNOL+1);	
    array_set( node_in_geometry, 0, length_nodes );
    array_set( nodes, 0, MNOL );
    array_set( control_data_initeldof_geometry, 0, 2 );

    db( CONTROL_DATA_INITELDOF_GEOMETRY, icontrol, control_data_initeldof_geometry, 
      ddum, ldum, VERSION_NORMAL, GET );

    db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );

    // determine which nodes are in the geometry
    for ( inod=0; inod<=max_node; inod++ ) {
      if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
         geometry( inod, ddum, control_data_initeldof_geometry, 
           in_geometry, rdum, ddum, rdum, ddum, NODE_START_REFINED, 
           CONTROL_DATA_INITELDOF_GEOMETRY, VERSION_NORMAL );
         if ( in_geometry ) node_in_geometry[inod] = 1;
      }
    }
	
      // set element_dof_initialised = 0 for elements which are totally in geometry
    for ( long int element=0; element<=max_element; element++ ) {
      if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) ) {
        db( ELEMENT, element, el, ddum, length_nodes, VERSION_NORMAL, GET );
        long int nnol = length_nodes - 1; 
	array_move( &el[1], nodes, nnol );
        long int all_in_geometry = 1;
        for ( long inol=0; inol<nnol; inol++ ) {
          inod = nodes[inol];
          if ( !node_in_geometry[inod] ) all_in_geometry = 0;
        }
        if ( all_in_geometry ) 
	  db( ELEMENT_DOF_INITIALISED, element, &zero, ddum, one, VERSION_NORMAL, PUT );        
      }
    }
    delete[] node_in_geometry;
  }

  if ( idat==-NODE || idat==-ELEMENT ) mesh_has_changed( VERSION_NORMAL );

  db_max_index( CHANGE_DATAITEM, max_change, VERSION_NORMAL, GET );
  if ( max_change>=0 ) {
    swit = set_swit(-1,-1,"data");
    for ( ichange=0; ichange<=max_change; ichange++ ) {
      if ( db_active_index( CHANGE_DATAITEM, ichange, VERSION_NORMAL ) ) {
        db( CHANGE_DATAITEM, ichange, change_dataitem, ddum, ldum, VERSION_NORMAL, GET );
        data_item_name = change_dataitem[0];
        data_item_index = change_dataitem[1];
        data_item_number = change_dataitem[2];
        operat = change_dataitem[3];
        change_dataitem_time_user = -NO;
        found = 0;
        db( CHANGE_DATAITEM_TIME_USER, ichange, &change_dataitem_time_user, ddum, ldum, 
          VERSION_NORMAL, GET_IF_EXISTS );
        if ( change_dataitem_time_user==-YES ) {
          user_change_dataitem_time( ichange, time_current, val );
          found = 1;
        }
        else {
          db( CHANGE_DATAITEM_TIME, ichange, idum, change_dataitem_time, length, 
            VERSION_NORMAL, GET );
          change_dataitem_time_discrete = -NO;
          db( CHANGE_DATAITEM_TIME_DISCRETE, ichange, &change_dataitem_time_discrete, 
            ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
          if ( change_dataitem_time_discrete==-YES ) {
            for ( iv=0; iv<length/2; iv++ ) {
              if ( time_current>(change_dataitem_time[iv*2+0]+EPS_TIME) ) {
                found = 1;
                val = change_dataitem_time[iv*2+1];
              }
            }
          }
          else {
            found = table_xy( change_dataitem_time, "CHANGE_DATAITEM_TIME",
              length, time_current, val );
          }
        }
        if ( found && db_active_index( data_item_name, data_item_index, VERSION_NORMAL ) ) {
          length = db_len( data_item_name, data_item_index, VERSION_NORMAL );
          if ( data_item_number<0 ) {
            db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET );
            array_member( dof_label, data_item_number, nuknwn, data_item_number );
            if ( length==npuknwn ) data_item_number /= nder;
          }
          if ( data_item_number<0 || data_item_number>length ) 
            db_error( CHANGE_DATAITEM, ichange );
          if ( db_type(data_item_name)==INTEGER ) db_error( CHANGE_DATAITEM, ichange );
          db( data_item_name, data_item_index, idum, dval, ldum, VERSION_NORMAL, GET );
          if      ( operat==-USE )
            dval[data_item_number] = val;
          else if ( operat==-ADD && task==-YES )
            dval[data_item_number] += val*dtime;
          else
            db_error( CHANGE_DATAITEM, ichange );
          db( data_item_name, data_item_index, idum, dval, ldum, VERSION_NORMAL, PUT );
        }
      }
    }
  }

  delete[] dof_label;
  delete[] integer_range;
  delete[] change_dataitem_time;
  delete[] dval;

}
