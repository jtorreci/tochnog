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
    any=0, all=0, length=0, itmp=0, iarea=0, count_in=0,
    max_area_element_group=0, area_element_group[3],
    element_group=0, method=0, ldum=0, el_name=0, use_node_list=0,
    length_nodes_list=0, *el=NULL, *nodes=NULL,
    *nodes_list=NULL;
  double rdum=0., ddum[MDIM];

  db_max_index( AREA_ELEMENT_GROUP, max_area_element_group, VERSION_NORMAL, GET );
  if ( max_area_element_group>=0 ) {
    el = get_new_int(MNOL+1);
    nodes = get_new_int(MNOL);
    nodes_list = get_new_int(DATA_ITEM_SIZE);
    db_max_index( ELEMENT, max_element, version, GET );
    for ( iarea=0; iarea<=max_area_element_group; iarea++ ) {
      if ( db_active_index( AREA_ELEMENT_GROUP, iarea, VERSION_NORMAL ) ) {
        db( AREA_ELEMENT_GROUP, iarea, area_element_group, ddum,
          ldum, VERSION_NORMAL, GET );
        method = -ALL;
        db( AREA_ELEMENT_GROUP_METHOD, iarea, &method, ddum,
          ldum, VERSION_NORMAL, GET_IF_EXISTS );
        // area_element_group_element: only elements with this element
        // name (manual 6.2); default all elements
        el_name = -ALL;
        db( AREA_ELEMENT_GROUP_ELEMENT, iarea, &el_name, ddum,
          ldum, VERSION_NORMAL, GET_IF_EXISTS );
        // area_element_group_node: direct global node numbers instead of
        // a geometry (manual 6.5); the last value of the master record is
        // the element group
        use_node_list = 0;
        if ( db_active_index( AREA_ELEMENT_GROUP_NODE, iarea,
            VERSION_NORMAL ) ) {
          db( AREA_ELEMENT_GROUP_NODE, iarea, nodes_list, ddum,
            length_nodes_list, VERSION_NORMAL, GET );
          use_node_list = 1;
        }
        // area_element_group_time -yes: evaluate at all times, not only
        // at the start (manual 6.6). The caller re-runs this routine per
        // step when any record asks for it (top.cc step_start).
        for ( element=0; element<=max_element; element++ ) {
          if ( db_active_index( ELEMENT, element, version ) ) {
            db( ELEMENT, element, el, ddum, length, version, GET );
            nnol = length - 1; array_move( &el[1], nodes, nnol );
            if ( el_name!=-ALL && el[0]!=el_name ) continue;
            // area_element_group_interface -no (the only legal value,
            // manual 6.3): interface elements are excluded from the
            // re-grouping (their group keeps GROUP_INTERFACE parameters)
            {
              long int grp_now = -1, interface_switch = -NO;
              db( ELEMENT_GROUP, element, &grp_now, ddum, ldum, version,
                GET_IF_EXISTS );
              if ( grp_now>=0 &&
                   db_active_index( GROUP_INTERFACE, grp_now, VERSION_NORMAL ) ) {
                db( AREA_ELEMENT_GROUP_INTERFACE, iarea, &interface_switch,
                  ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
                if ( interface_switch!=-YES ) continue;
              }
            }
            all = 1;
            any = 0;
            count_in = 0;
            for ( inol=0; inol<nnol; inol++ ) {
              inod = nodes[inol];
              if ( use_node_list ) {
                itmp = array_member( nodes_list, inod,
                  length_nodes_list, ldum );
              }
              else {
                geometry( inod, ddum, area_element_group, itmp, rdum, ddum, rdum,
                  ddum, NODE_START_REFINED, PROJECT_EXACT, version );
              }
              if ( !itmp ) all = 0;
              if ( itmp ) { any = 1; count_in++; }
            }
            if ( ( method==-ALL && all ) ||
                 ( method==-ANY && any ) ||
                 ( method==-ANY_BUT_NOT_ALL && any && !all ) ||
                 ( method>=0 && count_in>=method ) ) {
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
    delete[] nodes_list;
  }

}

// does any area_element_group record ask for evaluation at all times?
// (area_element_group_time -yes, manual 6.6)
long int area_element_group_time_active( void )
{
  long int iarea=0, max_area=0, swit=0, ldum=0;
  double ddum[1];
  db_max_index( AREA_ELEMENT_GROUP_TIME, max_area, VERSION_NORMAL, GET );
  for ( iarea=0; iarea<=max_area; iarea++ ) {
    if ( db_active_index( AREA_ELEMENT_GROUP_TIME, iarea, VERSION_NORMAL ) ) {
      db( AREA_ELEMENT_GROUP_TIME, iarea, &swit, ddum, ldum,
        VERSION_NORMAL, GET );
      if ( swit==-YES ) return 1;
    }
  }
  return 0;
}

void area_element_group_sequence( void )

{
  long int element=0, max_element=0, itime=0, inol=0, nnol=0, count_in=0,
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

  // area_element_group_sequence_element_group (with underscores) is the
  // Professional manual name of the legacy GNU keyword
  // area_element_group_sequence_elementgroup (without). Copy the record
  // when only the Professional name is used (safe PUT: this routine runs
  // in step_close, not in a parallel loop). This runs BEFORE the
  // max_index query of the legacy item.
  {
    long int ialias=0, max_alias=0, length_alias=0,
      alias_groups[DATA_ITEM_SIZE];
    db_max_index( AREA_ELEMENT_GROUP_SEQUENCE_ELEMENT_GROUP, max_alias,
      VERSION_NORMAL, GET );
    for ( ialias=0; ialias<=max_alias; ialias++ ) {
      if ( db_active_index( AREA_ELEMENT_GROUP_SEQUENCE_ELEMENT_GROUP,
           ialias, VERSION_NORMAL ) &&
           !db_active_index( AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP,
           ialias, VERSION_NORMAL ) ) {
        db( AREA_ELEMENT_GROUP_SEQUENCE_ELEMENT_GROUP, ialias,
          alias_groups, ddum, length_alias, VERSION_NORMAL, GET );
        db( AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP, ialias,
          alias_groups, ddum, length_alias, VERSION_NORMAL, PUT );
      }
    }
  }
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
        // manual Professional 6.9: as a completely separate option, NEITHER
        // geometry NOR the element list is used - then the elements of the
        // PREVIOUS group number group_(i-1) get the new group number
        // group_i at time_i (the previous group selects the elements).
        // The old GNU code required one of the selectors and errored out
        // (dam_building only uses _element_group + _time).
        area_element_group_sequence_element[0] = -ALL;
        db( AREA_ELEMENT_GROUP_SEQUENCE_ELEMENT, iarea, 
          area_element_group_sequence_element, ddum, 
          ldum, VERSION_NORMAL, GET_IF_EXISTS );
        method = -ALL;
        db( AREA_ELEMENT_GROUP_SEQUENCE_METHOD, iarea,
          &method, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
        // area_element_group_sequence_geometry_method (manual 6.11): the
        // -all/-any selection of the geometry, kept separate from
        // AREA_ELEMENT_GROUP_SEQUENCE_METHOD for Professional parity;
        // it wins for the geometry test when both are given
        {
          long int geometry_method = 0;
          if ( db( AREA_ELEMENT_GROUP_SEQUENCE_GEOMETRY_METHOD, iarea,
              &geometry_method, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS ) )
            method = geometry_method;
        }
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
              count_in = 0;
              for ( inol=0; inol<nnol; inol++ ) {
                inod = nodes[inol];
                geometry( inod, ddum, area_element_group_sequence_geometry,
                  itmp, rdum, ddum, rdum, ddum, NODE_START_REFINED,
                  PROJECT_EXACT, VERSION_NORMAL );
                if ( !itmp ) all = 0;
                if ( itmp ) { any = 1; count_in++; }
              }
              if ( ( method==-ALL && all ) ||
                   ( method==-ANY && any ) ||
                   ( method==-ANY_BUT_NOT_ALL && any && !all ) ||
                   ( method>=0 && count_in>=method ) ) {
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
            // area_element_group_sequence_interface -no (the only legal
            // value, manual 6.12): interface elements are excluded
            if ( ok ) {
              long int grp_now = -1, interface_switch = -NO;
              db( ELEMENT_GROUP, element, &grp_now, ddum, ldum,
                VERSION_NORMAL, GET_IF_EXISTS );
              if ( grp_now>=0 &&
                   db_active_index( GROUP_INTERFACE, grp_now, VERSION_NORMAL ) ) {
                db( AREA_ELEMENT_GROUP_SEQUENCE_INTERFACE, iarea,
                  &interface_switch, ddum, ldum, VERSION_NORMAL,
                  GET_IF_EXISTS );
                if ( interface_switch!=-YES ) ok = 0;
              }
            }
            if ( ok ) {
              // mode without geometry/element selectors: the element must
              // currently have the PREVIOUS group in the sequence (the
              // i-th time point maps group_(i-1) -> group_i)
              if ( !use_geometry && !use_element ) {
                long int grp_now2 = -1;
                db( ELEMENT_GROUP, element, &grp_now2, ddum, ldum,
                  VERSION_NORMAL, GET_IF_EXISTS );
                if ( grp_now2>=0 ) {
                  long int prev_grp = -1;
                  for ( long int it2=0; it2<length_elementgroup; it2++ ) {
                    double tt2 = area_element_group_sequence_time[it2];
                    if ( time_total>=(tt2-EPS_SMALL) )
                      prev_grp = area_element_group_sequence_elementgroup[it2];
                  }
                  // find the group BEFORE the current time window
                  long int prev2 = -1;
                  for ( long int it2=0; it2<length_elementgroup; it2++ ) {
                    double tt2 = area_element_group_sequence_time[it2];
                    if ( time_total < (tt2-EPS_SMALL) ) break;
                    prev2 = area_element_group_sequence_elementgroup[it2];
                  }
                  if ( prev2<0 ) prev2 = prev_grp;
                  if ( grp_now2!=prev2 ) ok = 0;
                }
              }
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

  // dependency_apply (manual Professional 6.125, global) and
  // control_dependency_apply (manual Professional 6.126, per timestep
  // index): when the switch is -no the dependency_item/dependency_diagram
  // machinery is disabled. Precedence: control_* (same index as the
  // current timestep) overrides the global record; default -yes.
  {
    long int icontrol=0, dependency_apply=-YES;
    db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
    db( DEPENDENCY_APPLY, 0, &dependency_apply, ddum, ldum,
      VERSION_NORMAL, GET_IF_EXISTS );
    db( CONTROL_DEPENDENCY_APPLY, icontrol, &dependency_apply, ddum, ldum,
      VERSION_NORMAL, GET_IF_EXISTS );
    if ( dependency_apply==-NO ) {
      return db( idat, gr, idum, values, nvalue, VERSION_NORMAL, task );
    }
  }

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
       // group_materi_plasti_bounda (manual Professional 6.231): the
       // listed values are the indices of the bounda_dof records that put
       // the element on a wall ("an element is on a wall when at least
       // one of the velocities of the elements is prescribed via
       // bounda_dof"). When a listed value matches an ACTIVE bounda_dof
       // record the element is on the wall when one of its nodes is
       // bounded (node_bounded) on the velocity or displacement parts.
       if ( !plasti_on_boundary ) {
         long int ib=0;
         for ( ib=0; ib<length && !plasti_on_boundary; ib++ ) {
           long int iboun = group_materi_plasti_boundary[ib];
           if ( iboun>=0 && db_active_index( BOUNDA_DOF, iboun,
               VERSION_NORMAL ) ) {
             long int ip=0, *node_bounded=NULL;
             if ( !db_active_index( NODE_BOUNDED, inod, VERSION_NORMAL ) )
               continue;
             node_bounded = db_int( NODE_BOUNDED, inod, VERSION_NORMAL );
             // velocity/displacement parts (nder = 1 in the current
             // formulation: one value per dof part)
             for ( ip=0; ip<npuknwn && !plasti_on_boundary; ip++ ) {
               if ( node_bounded[ip] &&
                    ( ( vel_indx>=0 && ip>=vel_indx && ip<vel_indx+ndim ) ||
                      ( dis_indx>=0 && ip>=dis_indx && ip<dis_indx+ndim ) ) )
                 plasti_on_boundary = 1;
             }
           }
         }
       }
     }
   }

}

long int node_attached_element_groups( long int inod, long int groups[],
  long int &n )

// Fills groups[0..n) with the distinct element groups of the elements
// the node belongs to (unsorted). Returns 1 when the node is a node of
// at least one active element, 0 otherwise. Used by the
// geometry_element_group filter (geometry.cc) and by the region/group
// restriction of post_calcul_static_pressure_height (groundfl.cc).

{
  long int ielem=0, max_element=0, inol=0, length=0, igroup=0, k=0,
    ldum=0, idum[1];
  double ddum[1];
  long int *el=NULL;

  n = 0;
  db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
  for ( ielem=0; ielem<=max_element; ielem++ ) {
    if ( !db_active_index( ELEMENT, ielem, VERSION_NORMAL ) ) continue;
    el = db_int( ELEMENT, ielem, VERSION_NORMAL );
    length = db_len( ELEMENT, ielem, VERSION_NORMAL );
    for ( inol=0; inol+1<length; inol++ ) {
      if ( el[inol+1]==inod ) {
        igroup = 0;
        db( ELEMENT_GROUP, ielem, &igroup, ddum, ldum,
          VERSION_NORMAL, GET_IF_EXISTS );
        for ( k=0; k<n; k++ )
          if ( groups[k]==igroup ) break;
        if ( k==n ) groups[n++] = igroup;
        break;
      }
    }
  }
  return ( n>0 );
}
