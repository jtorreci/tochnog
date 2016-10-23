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

void area_element_group( long int version )

{
  long int element=0, max_element=0, inol=0, nnol=0, inod=0, 
    any=0, all=0, length=0, itmp=0, iarea=0, 
    max_area_element_group=0, area_element_group[3], 
    element_group=0, method=0, ldum=0, *el=NULL, *nodes=NULL;
  double rdum=0., ddum[MDIM];

  db_max_index( AREA_ELEMENT_GROUP, max_area_element_group, VERSION_NORMAL, GET );
  if ( max_area_element_group>=0 ) {
    el = get_new_int(MNOL+1);
    nodes = get_new_int(MNOL);
    db_max_index( ELEMENT, max_element, version, GET );
    for ( iarea=0; iarea<=max_area_element_group; iarea++ ) {
      if ( db_active_index( AREA_ELEMENT_GROUP, iarea, VERSION_NORMAL ) ) {
        db( AREA_ELEMENT_GROUP, iarea, area_element_group, ddum, 
          ldum, VERSION_NORMAL, GET );
        method = -ALL;
        db( AREA_ELEMENT_GROUP_METHOD, iarea, &method, ddum, 
          ldum, VERSION_NORMAL, GET_IF_EXISTS );
        for ( element=0; element<=max_element; element++ ) {
          if ( db_active_index( ELEMENT, element, version ) ) {
            db( ELEMENT, element, el, ddum, length, version, GET );
            nnol = length - 1; array_move( &el[1], nodes, nnol );
            all = 1;
            any = 0;
            for ( inol=0; inol<nnol; inol++ ) {
              inod = nodes[inol];
              geometry( inod, ddum, area_element_group, itmp, rdum, ddum, rdum,
                ddum, NODE_START_REFINED, PROJECT_EXACT, version );
              if ( !itmp ) all = 0;
              if ( itmp ) any = 1;
            }
            if ( ( method==-ALL && all ) || ( method==-ANY && any ) ) {
              element_group = area_element_group[2];
              length = 1; db( ELEMENT_GROUP, element, &element_group, ddum, 
                length, version, PUT );
              length = 1; db( ELEMENT_GROUP_AREA_ELEMENT_GROUP, 
                element, &iarea, ddum, length, VERSION_NORMAL, PUT );
            }
          }
        }
      }
    }
    delete[] el;
    delete[] nodes;
  }

}

void area_element_group_sequence( void )

{
  long int element=0, max_element=0, itime=0, inol=0, nnol=0, 
    inod=0, ok=0, length=0, name=0, length_elementgroup=0,
    itmp=0, all=0, any=0, method=0, found=0, iarea=0, max_area_element_group=0, 
    length_area_element_group_sequence=0,
    area_element_group_sequence_element[1], 
    area_element_group_sequence_elementgroup[DATA_ITEM_SIZE], 
    area_element_group_sequence_geometry[2], 
    element_group=0, use_geometry=0, use_element=0,
    ldum=0, idum[1], *el=NULL, *nodes=NULL, 
    *area_element_group_sequence=NULL;
  double time=0., time_total=0., rdum=0., 
    ddum[MDIM], area_element_group_sequence_time[DATA_ITEM_SIZE];

  db_max_index( AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP, 
    max_area_element_group, VERSION_NORMAL, GET );
  if ( max_area_element_group>=0 ) {
    el = get_new_int(MNOL+1);
    nodes = get_new_int(MNOL);
    db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
    db( TIME_CURRENT, 0, idum, &time_total, ldum, 
      VERSION_NORMAL, GET );
    for ( iarea=0; iarea<=max_area_element_group; iarea++ ) {
      if ( db_active_index( AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP, 
          iarea, VERSION_NORMAL ) ) {
        use_geometry = 0;
        if ( db_active_index( AREA_ELEMENT_GROUP_SEQUENCE_GEOMETRY, 
            iarea, VERSION_NORMAL ) ) {
          use_geometry = 1;
          db( AREA_ELEMENT_GROUP_SEQUENCE_GEOMETRY, iarea, 
            area_element_group_sequence_geometry, ddum, 
            ldum, VERSION_NORMAL, GET );
        }
        use_element = 0;
        if ( db_active_index( AREA_ELEMENT_GROUP_SEQUENCE, 
            iarea, VERSION_NORMAL ) ) {
          use_element = 1;
          area_element_group_sequence = db_int( AREA_ELEMENT_GROUP_SEQUENCE, 
            iarea, VERSION_NORMAL );
          length_area_element_group_sequence = db_len( AREA_ELEMENT_GROUP_SEQUENCE,
            iarea, VERSION_NORMAL );
        }
        if ( !use_geometry && !use_element ) {
          pri( "Error: AREA_ELEMENT_GROUP_SEQUENCE_GEOMETRY or AREA_ELEMENT_GROUP_SEQUENCE should be specified." );
          exit(TN_EXIT_STATUS);
        }
        area_element_group_sequence_element[0] = -ALL;
        db( AREA_ELEMENT_GROUP_SEQUENCE_ELEMENT, iarea, 
          area_element_group_sequence_element, ddum, 
          ldum, VERSION_NORMAL, GET_IF_EXISTS );
        method = -ALL;
        db( AREA_ELEMENT_GROUP_SEQUENCE_METHOD, iarea, 
          &method, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
        db( AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP, iarea, 
          area_element_group_sequence_elementgroup, ddum, 
          length_elementgroup, VERSION_NORMAL, GET );
        db( AREA_ELEMENT_GROUP_SEQUENCE_TIME, iarea, 
          idum, area_element_group_sequence_time,
          length_elementgroup, VERSION_NORMAL, GET_AND_CHECK );
        for ( element=0; element<=max_element; element++ ) {
          if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) ) {
            ok = 0;
            if ( use_geometry ) {
              db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
              nnol = length - 1; array_move( &el[1], nodes, nnol );
              name = el[0];
              all = 1;
              any = 0;
              for ( inol=0; inol<nnol; inol++ ) {
                inod = nodes[inol];
                geometry( inod, ddum, area_element_group_sequence_geometry, 
                  itmp, rdum, ddum, rdum, ddum, NODE_START_REFINED, 
                  PROJECT_EXACT, VERSION_NORMAL );
                if ( !itmp ) all = 0;
                if ( itmp ) any = 1;
              }
              if ( ( method==-ALL && all ) || ( method==-ANY && any ) ) {
                ok = 1;
              }
              if ( area_element_group_sequence_element[0]!=-ALL ) {
                if ( name!=area_element_group_sequence_element[0] ) ok = 0;
              }
            }
            if ( use_element ) {
              if ( array_member( area_element_group_sequence, element,
                length_area_element_group_sequence, ldum ) ) ok = 1;
            }
            if ( ok ) {
              found = 0;
              for ( itime=0; itime<length_elementgroup; itime++ ) {
                time = area_element_group_sequence_time[itime];
                if ( time_total>=(time-EPS_SMALL) ) {
                  element_group = area_element_group_sequence_elementgroup[itime];
                  found = 1;
                }
              }
              if ( found ) {
                length = 1; db( ELEMENT_GROUP, element, &element_group, ddum, 
                  length, VERSION_NORMAL, PUT );
                length = 1; db( ELEMENT_GROUP_AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP, 
                  element, &iarea, ddum, length, VERSION_NORMAL, PUT );
              }
            }
          }
        }
      }
    }
    delete[] el;
    delete[] nodes;
  }

}

long int get_group_data( long int idat, long int gr, long int element,
  double new_unknowns[], double values[], long int &nvalue, long int task )

{
  long int i=0, n=0, idep=0, max_dep=0, found=0, go_ahead=0,
    iuknwn=0, ival=0, nval=0, length=0,
    data_item_name=0, data_item_number=0, 
    idistribute=0, ndistribute=0, ldum=0, 
    idum[1], element_distribute[DATA_ITEM_SIZE], dof_label[MUKNWN], *dependency_item=NULL;
  double tmp=0., time_current=0., dtime=0., time_left=0., time_right=0,
    val_left=0., val_right=0., ddum[1], element_distribute_values[DATA_ITEM_SIZE], 
    *dependency_diagram=NULL;

  db_max_index( DEPENDENCY_ITEM, max_dep, VERSION_NORMAL, GET );
  if ( max_dep>=0 ) {
    db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    for ( idep=0; idep<=max_dep && !found; idep++ ) {
      if ( db_active_index( DEPENDENCY_ITEM, idep, VERSION_NORMAL ) ) {
        dependency_item = db_int( DEPENDENCY_ITEM, idep, VERSION_NORMAL );
        n = dependency_item[3];
        if ( n<2 ) db_error( DEPENDENCY_ITEM, idep );
        if ( labs(dependency_item[0])==idat && dependency_item[1]==gr ) {
          if ( dependency_item[2]==-TIME_CURRENT ) {
            db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET );
            db( TIME_CURRENT, 0, idum, &time_current, 
              ldum, VERSION_NORMAL, GET_IF_EXISTS );
            go_ahead = 1;
            tmp = time_current + dtime;
          }
          else {
            array_member(dof_label,dependency_item[2],nuknwn,iuknwn);
            if ( iuknwn>=0 && iuknwn<=nuknwn ){
              go_ahead = 1;
              tmp = new_unknowns[iuknwn];
            }
            else
              db_error( DEPENDENCY_ITEM, idep );
          }
          if ( go_ahead ) {
            found = 1;
            dependency_diagram = db_dbl( DEPENDENCY_DIAGRAM, idep, VERSION_NORMAL );
            length = db_len( DEPENDENCY_DIAGRAM, idep, VERSION_NORMAL );
            if ( length%n!=0 ) db_error( DEPENDENCY_DIAGRAM, idep );
            if ( db_fixed_length( idat ) )
              nval = db_data_length( idat );
            else
              nval = ( length -  n ) / n;
            for ( ival=0; ival<nval; ival++ ) {
              if      ( tmp<dependency_diagram[0] ) 
                values[ival] = dependency_diagram[n+ival*n+0];
              else if ( tmp>dependency_diagram[n-1] ) 
                values[ival] = dependency_diagram[n+ival*n+n-1];
              else {
                for( i=0; i<n-1; i++ ) {
                  time_left = dependency_diagram[i];
                  time_right = dependency_diagram[i+1];
                  val_left = dependency_diagram[n+ival*n+i];
                  val_right = dependency_diagram[n+ival*n+i+1];
                  if ( time_right<=time_left ) db_error( DEPENDENCY_DIAGRAM, idep );
                  if ( tmp>=time_left && tmp<=time_right ) {
                    values[ival] = val_left + (tmp-dependency_diagram[i])*
                      (val_right-val_left) / (time_right-time_left);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if ( found ) 
    nvalue = nval;
  else
    found = db( idat, gr, idum, values, nvalue, VERSION_NORMAL, task );

  if ( found ) {
    if ( db_active_index( ELEMENT_DISTRIBUTE, element, VERSION_NORMAL ) ) {
      db( ELEMENT_DISTRIBUTE, element, element_distribute, ddum, 
        ldum, VERSION_NORMAL, GET );
      db( ELEMENT_DISTRIBUTE_VALUES, element, idum, element_distribute_values,
        ndistribute, VERSION_NORMAL, GET );
      for ( idistribute=0; idistribute<ndistribute; idistribute++ ) {
        data_item_name = element_distribute[idistribute*2+0];
        data_item_number = element_distribute[idistribute*2+1];
        if ( labs(idat)==labs(data_item_name) ) {
          if ( data_item_number<0 || data_item_number>nvalue-1 ) {
            pri( "Error detected in CONTROL_DISTRIBUTE." );
            exit(1);
          }
          values[data_item_number] += element_distribute_values[idistribute];
        }
      }
    }
  }

  return found;
}

void group_materi_plasti_boundary_evaluate( long int nodes[], long int nnol,
  long int element_group, long int &plasti_on_boundary )

{
   long int inol=0, inod=0, length=0, iel=0, nel=0, elnum=0, gr=0, ldum=0,
     group_materi_plasti_boundary[DATA_ITEM_SIZE], *node_element=NULL;
   double ddum[1];

   plasti_on_boundary = 0;

   if ( db( GROUP_MATERI_PLASTI_BOUNDARY, element_group, group_materi_plasti_boundary, 
       ddum, length, VERSION_NORMAL, GET_IF_EXISTS ) ) {
     for ( inol=0; inol<nnol && !plasti_on_boundary; inol++ ) {
       inod = nodes[inol];
       node_element = db_int( NODE_ELEMENT, inod, VERSION_NORMAL );
       nel = db_len( NODE_ELEMENT, inod, VERSION_NORMAL );
       for ( iel=0; iel<nel && !plasti_on_boundary; iel++ ) {
         elnum = node_element[iel];
         gr = 0;
         db( ELEMENT_GROUP, elnum, &gr, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
         if ( array_member( group_materi_plasti_boundary, gr, length, ldum ) )
           plasti_on_boundary = 1;
       }
     }
   }

}
