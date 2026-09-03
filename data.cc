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

// helper for control_data_copy / control_data_copy_index (manual
// Professional 6.117/6.119): copy one record index from idat_from to
// idat_to with an optional multiplication factor. Integer copies require
// factor 1; the manual's special node_inertia -> node_force case is a
// double copy with factor -1 (d'alembert), which this covers naturally.
static void data_copy_apply( long int idat_from, long int index_from,
  long int idat_to, long int index_to, double factor,
  long int idat_control, long int icontrol )
{
  long int ldum=0, idum[1], length_from=0, i=0;
  double ddum[1], *dval_copy=NULL;
  long int *ival_copy=NULL;

  ival_copy = get_new_int(DATA_ITEM_SIZE);
  dval_copy = get_new_dbl(DATA_ITEM_SIZE);

  length_from = db_len( idat_from, index_from, VERSION_NORMAL );
  if ( db_type(idat_from)==INTEGER && db_type(idat_to)==INTEGER ) {
    if ( scalar_dabs(factor-1.)>TINY ) db_error( idat_control, icontrol );
    db( idat_from, index_from, ival_copy, ddum, length_from,
      VERSION_NORMAL, GET );
    db( idat_to, index_to, ival_copy, ddum, length_from,
      VERSION_NORMAL, PUT );
  }
  else if ( db_type(idat_from)==DOUBLE_PRECISION &&
            db_type(idat_to)==DOUBLE_PRECISION ) {
    db( idat_from, index_from, idum, dval_copy, length_from,
      VERSION_NORMAL, GET );
    if ( scalar_dabs(factor-1.)>TINY )
      for ( i=0; i<length_from; i++ ) dval_copy[i] *= factor;
    db( idat_to, index_to, idum, dval_copy, length_from,
      VERSION_NORMAL, PUT );
  }
  else
    db_error( idat_control, icontrol );
}

void data( long int task, double dtime, double time_current )

{
  long int idat=0, in=0, iv=0, index=0, range_length=0, icontrol=0, length=0, 
    swit=0, max_index=0, inod=0, max_node=0, found=0, 
    ichange=0, max_change=0, idim=0, operat=0, ldum=0, 
    ireset=0, max_reset=0, idof_reset=0, idof_value=0, ireset_val=0,
    length_diagram=0, change_dataitem_apply=-YES,
    reset_method=-USE,
    data_item_name=0, data_item_index=0, data_item_number=0,
    change_dataitem_time_discrete=-NO, change_dataitem_time_user=0,
    change_dataitem_time_method=0, change_dataitem_geometry[2]={0,0},
    idum[1], change_dataitem[4], *dof_label=NULL, *integer_range=NULL, 
    *data_delete=NULL, *data_put=NULL, *reset_dof=NULL, *reset_value_dof=NULL;
  double rdum=0., val=0., ddum[MDIM], *change_dataitem_time=NULL, 
    *dval=NULL, *coord=NULL, *node_dof=NULL, *reset_value_diagram=NULL,
    reset_value_constant=0., coords[MDIM];

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

  // data_activate / data_delete (manual Professional 6.397-6.400):
  // time-gated variants of control_data_activate/control_data_delete
  // (not per control-timestep). Evaluated when time_current >= the
  // data_*_time point (default: at the start of the calculation).
  {
    long int idata_recs = 0, max_idata = 0, iv2 = 0, idat2 = 0,
      max_del2 = 0, *data_list = NULL, switch_val = -YES, length_list = 0;
    double time_gate = 0.;

    db_max_index( DATA_ACTIVATE, max_idata, VERSION_NORMAL, GET );
    for ( idata_recs=0; idata_recs<=max_idata; idata_recs++ ) {
      if ( !db_active_index( DATA_ACTIVATE, idata_recs, VERSION_NORMAL ) )
        continue;
      time_gate = 0.;
      db( DATA_ACTIVATE_TIME, idata_recs, idum, &time_gate, ldum,
        VERSION_NORMAL, GET_IF_EXISTS );
      if ( time_current<(time_gate-EPS_TIME) ) continue;
      data_list = db_int( DATA_ACTIVATE, idata_recs, VERSION_NORMAL );
      length_list = db_len( DATA_ACTIVATE, idata_recs, VERSION_NORMAL );
      switch_val = data_list[length_list-1];
      if ( switch_val!=-NO ) continue;
      for ( iv2=0; iv2<length_list-1; iv2++ ) {
        idat2 = data_list[iv2];
        if ( idat2>=0 ) continue;
        db_max_index( idat2, max_del2, VERSION_NORMAL, GET );
        for ( index=0; index<=max_del2; index++ )
          if ( db_active_index( idat2, index, VERSION_NORMAL ) )
            db_delete_index( idat2, index, VERSION_NORMAL );
      }
    }

    db_max_index( DATA_DELETE, max_idata, VERSION_NORMAL, GET );
    for ( idata_recs=0; idata_recs<=max_idata; idata_recs++ ) {
      if ( !db_active_index( DATA_DELETE, idata_recs, VERSION_NORMAL ) )
        continue;
      time_gate = 0.;
      db( DATA_DELETE_TIME, idata_recs, idum, &time_gate, ldum,
        VERSION_NORMAL, GET_IF_EXISTS );
      if ( time_current<(time_gate-EPS_TIME) ) continue;
      data_list = db_int( DATA_DELETE, idata_recs, VERSION_NORMAL );
      length_list = db_len( DATA_DELETE, idata_recs, VERSION_NORMAL );
      idat2 = data_list[0];
      if ( idat2>=0 ) continue;
      if ( data_list[1]>=0 ) {
        if      ( idat2==-ELEMENT ) delete_element( data_list[1], VERSION_NORMAL );
        else if ( idat2==-NODE ) delete_node( data_list[1], VERSION_NORMAL );
        else
          db_delete_index( idat2, data_list[1], VERSION_NORMAL );
      }
      else if ( data_list[1]==-RA ) {
        range_expand( &data_list[1], integer_range, length_list, range_length );
        for ( in=0; in<range_length; in++ ) {
          index = integer_range[in];
          if      ( idat2==-ELEMENT ) delete_element( index, VERSION_NORMAL );
          else if ( idat2==-NODE ) delete_node( index, VERSION_NORMAL );
          else
            db_delete_index( idat2, index, VERSION_NORMAL );
        }
      }
      else if ( data_list[1]==-ALL ) {
        db_max_index( idat2, max_del2, VERSION_NORMAL, GET );
        for ( index=0; index<=max_del2; index++ ) {
          if      ( idat2==-ELEMENT ) delete_element( index, VERSION_NORMAL );
          else if ( idat2==-NODE ) delete_node( index, VERSION_NORMAL );
          else
            db_delete_index( idat2, index, VERSION_NORMAL );
        }
      }
    }
  }

  // print_mesh_dof (manual Professional 6.31-6.33 via bounda_print_mesh_dof):
  // one-shot dump of node coordinates and the listed dof values (all
  // dofs when none are listed) at the first evaluation; nodes can be
  // restricted to a geometry. Written to print_mesh_dof.dat.
  {
    static long int print_mesh_dof_done = 0;
    long int *pmd_dofs = NULL, length_pmd = 0, pmd_geometry[2] = {0,0},
      in_geometry_pmd = 0, idof_pmd = 0, indx_pmd = 0;
    double *coord_pmd = NULL, *node_dof_pmd = NULL, rdum_pmd = 0.;

    if ( !print_mesh_dof_done &&
         db_active_index( PRINT_MESH_DOF, 0, VERSION_NORMAL ) ) {
      print_mesh_dof_done = 1;
      pmd_dofs = db_int( PRINT_MESH_DOF, 0, VERSION_NORMAL );
      length_pmd = db_len( PRINT_MESH_DOF, 0, VERSION_NORMAL );
      if ( db_active_index( PRINT_MESH_DOF_GEOMETRY, 0, VERSION_NORMAL ) )
        db( PRINT_MESH_DOF_GEOMETRY, 0, pmd_geometry, ddum, ldum,
          VERSION_NORMAL, GET );
      ofstream pmd_out( "print_mesh_dof.dat" );
      for ( inod=0; inod<=max_node; inod++ ) {
        if ( !db_active_index( NODE, inod, VERSION_NORMAL ) ) continue;
        if ( pmd_geometry[0] ) {
          geometry( inod, ddum, pmd_geometry, in_geometry_pmd, rdum_pmd,
            ddum, rdum_pmd, ddum, NODE_START_REFINED, PROJECT_EXACT,
            VERSION_NORMAL );
          if ( !in_geometry_pmd ) continue;
        }
        coord_pmd = db_dbl( NODE_START_REFINED, inod, VERSION_NORMAL );
        node_dof_pmd = db_dbl( NODE_DOF, inod, VERSION_NORMAL );
        length = db_len( NODE_DOF, inod, VERSION_NORMAL );
        pmd_out << inod;
        for ( idim=0; idim<ndim; idim++ )
          pmd_out << " " << coord_pmd[idim];
        for ( iv=0; iv<length_pmd; iv++ ) {
          indx_pmd = pmd_dofs[iv];
          if ( indx_pmd<0 ) {
            array_member( dof_label, indx_pmd, nuknwn, indx_pmd );
            if ( length==npuknwn ) indx_pmd /= nder;
          }
          if ( indx_pmd>=0 && indx_pmd<length )
            pmd_out << " " << node_dof_pmd[indx_pmd];
        }
        pmd_out << "\n";
      }
      pmd_out.close();
    }
  }

  // control_data_activate (manual Professional 6.114): activate/deactivate
  // data items. -no deletes all records of the listed items (the records
  // stop being used by the solvers); -yes is a no-op: input records are
  // active by default, and the GNU deletion is destructive, so
  // re-activation would require re-putting the records (documented
  // difference with Professional).
  if ( db_active_index( CONTROL_DATA_ACTIVATE, icontrol, VERSION_NORMAL ) ) {
    long int *data_activate=NULL, activate_switch=0, idat_activate=0,
      max_activate=0;
    data_activate = db_int( CONTROL_DATA_ACTIVATE, icontrol, VERSION_NORMAL );
    length = db_len( CONTROL_DATA_ACTIVATE, icontrol, VERSION_NORMAL );
    activate_switch = data_activate[length-1];
    for ( in=0; in<length-1; in++ ) {
      idat_activate = data_activate[in];
      if ( idat_activate>=0 )
        db_error( CONTROL_DATA_ACTIVATE, icontrol );
      if ( activate_switch==-NO ) {
        db_max_index( idat_activate, max_activate, VERSION_NORMAL, GET );
        for ( index=0; index<=max_activate; index++ )
          if ( db_active_index( idat_activate, index, VERSION_NORMAL ) )
            db_delete_index( idat_activate, index, VERSION_NORMAL );
      }
    }
  }

  // control_data_arithmetic (manual Professional 6.115/6.116): change the
  // selected data item with the value of control_data_arithmetic_double
  // (same index) using -plus/-minus/-multiply/-divide. Record layout:
  //   [ item, index | -RA range..., number | -ALL, operat ]
  // number and operat are the LAST two slots; -ALL applies the operation
  // to all numbers of the record.
  if ( db_active_index( CONTROL_DATA_ARITHMETIC, icontrol, VERSION_NORMAL ) ) {
    long int *data_arith=NULL, arith_operat=0, arith_number=0, idat_arith=0,
      iindx=0, length_record=0, number_indx=0, *arith_doflabel=NULL;
    double arith_val=0.;
    db( CONTROL_DATA_ARITHMETIC_DOUBLE, icontrol, idum, &arith_val, ldum,
      VERSION_NORMAL, GET );
    data_arith = db_int( CONTROL_DATA_ARITHMETIC, icontrol, VERSION_NORMAL );
    length = db_len( CONTROL_DATA_ARITHMETIC, icontrol, VERSION_NORMAL );
    idat_arith = data_arith[0];
    arith_number = data_arith[length-2];
    arith_operat = data_arith[length-1];
    if ( idat_arith>=0 || db_type(idat_arith)==INTEGER )
      db_error( CONTROL_DATA_ARITHMETIC, icontrol );
    if ( arith_operat==-DIVIDE && scalar_dabs(arith_val)<TINY )
      db_error( CONTROL_DATA_ARITHMETIC, icontrol );
    if ( data_arith[1]==-RA )
      range_expand( &data_arith[1], integer_range, length, range_length );
    else {
      integer_range[0] = data_arith[1];
      range_length = 1;
    }
    arith_doflabel = get_new_int(MUKNWN);
    for ( iindx=0; iindx<range_length; iindx++ ) {
      index = integer_range[iindx];
      if ( db_active_index( idat_arith, index, VERSION_NORMAL ) ) {
        length_record = db_len( idat_arith, index, VERSION_NORMAL );
        db( idat_arith, index, idum, dval, ldum, VERSION_NORMAL, GET );
        if ( arith_number==-ALL ) {
          for ( in=0; in<length_record; in++ ) {
            if      ( arith_operat==-PLUS )    dval[in] += arith_val;
            else if ( arith_operat==-MINUS )   dval[in] -= arith_val;
            else if ( arith_operat==-MULTIPLY ) dval[in] *= arith_val;
            else if ( arith_operat==-DIVIDE )  dval[in] /= arith_val;
          }
        }
        else {
          number_indx = arith_number;
          if ( number_indx<0 ) {
            db( DOF_LABEL, 0, arith_doflabel, ddum, ldum,
              VERSION_NORMAL, GET );
            array_member( arith_doflabel, number_indx, nuknwn, number_indx );
            if ( length_record==npuknwn ) number_indx /= nder;
          }
          if ( number_indx<0 || number_indx>length_record-1 )
            db_error( CONTROL_DATA_ARITHMETIC, icontrol );
          if      ( arith_operat==-PLUS )    dval[number_indx] += arith_val;
          else if ( arith_operat==-MINUS )   dval[number_indx] -= arith_val;
          else if ( arith_operat==-MULTIPLY ) dval[number_indx] *= arith_val;
          else if ( arith_operat==-DIVIDE )  dval[number_indx] /= arith_val;
          else
            db_error( CONTROL_DATA_ARITHMETIC, icontrol );
        }
        db( idat_arith, index, idum, dval, length_record,
          VERSION_NORMAL, PUT );
      }
    }
  }

  // control_data_copy (manual Professional 6.117): copy ALL indices of
  // data_item_from to data_item_to, with optional multiplication factor
  // from control_data_copy_factor (same index).
  if ( db_active_index( CONTROL_DATA_COPY, icontrol, VERSION_NORMAL ) ) {
    long int *data_copy=NULL, idat_from=0, idat_to=0, max_copy=0;
    double copy_factor=1.;
    db( CONTROL_DATA_COPY_FACTOR, icontrol, idum, &copy_factor, ldum,
      VERSION_NORMAL, GET_IF_EXISTS );
    data_copy = db_int( CONTROL_DATA_COPY, icontrol, VERSION_NORMAL );
    idat_from = data_copy[0];
    idat_to = data_copy[1];
    if ( idat_from>=0 || idat_to>=0 ) db_error( CONTROL_DATA_COPY, icontrol );
    db_max_index( idat_from, max_copy, VERSION_NORMAL, GET );
    for ( index=0; index<=max_copy; index++ )
      if ( db_active_index( idat_from, index, VERSION_NORMAL ) )
        data_copy_apply( idat_from, index, idat_to, index, copy_factor,
          CONTROL_DATA_COPY, icontrol );
  }

  // control_data_copy_index (manual Professional 6.119): copy a single
  // record index_from of data_item_from to index_to of data_item_to, with
  // optional multiplication factor from control_data_copy_index_factor.
  if ( db_active_index( CONTROL_DATA_COPY_INDEX, icontrol, VERSION_NORMAL ) ) {
    long int *data_copy=NULL, idat_from=0, idat_to=0;
    double copy_factor=1.;
    db( CONTROL_DATA_COPY_INDEX_FACTOR, icontrol, idum, &copy_factor, ldum,
      VERSION_NORMAL, GET_IF_EXISTS );
    data_copy = db_int( CONTROL_DATA_COPY_INDEX, icontrol, VERSION_NORMAL );
    idat_from = data_copy[0];
    idat_to = data_copy[2];
    if ( idat_from>=0 || idat_to>=0 )
      db_error( CONTROL_DATA_COPY_INDEX, icontrol );
    if ( db_active_index( idat_from, data_copy[1], VERSION_NORMAL ) )
      data_copy_apply( idat_from, data_copy[1], idat_to, data_copy[3],
        copy_factor, CONTROL_DATA_COPY_INDEX, icontrol );
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
    db( CONTROL_CHANGE_DATAITEM_APPLY, icontrol, &change_dataitem_apply,
      ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    if ( change_dataitem_apply!=-NO ) {
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
        // change_dataitem_time_method (manual Professional 6.51): the time
        // table contains cosinus, sinus or tangent values; the inverse
        // trigonometric angle is stored instead of the value itself.
        // Typically used for phi-c reduction: the table gives tan(phi),
        // the stored parameter is phi = atan(val), so cohesion and tangent
        // of the friction angle can be decreased at the same ratio.
        change_dataitem_time_method = 0;
        db( CHANGE_DATAITEM_TIME_METHOD, ichange, &change_dataitem_time_method,
          ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
        if ( found ) {
          if      ( change_dataitem_time_method==-COSINUS ) val = acos(val);
          else if ( change_dataitem_time_method==-SINUS )   val = asin(val);
          else if ( change_dataitem_time_method==-TANGENT ) val = atan(val);
        }

        // change_dataitem_geometry (manual Professional 6.48): restrict the
        // change of a group_* record to the elements of the group that are
        // inside the geometry. Materialized by splitting the group: on the
        // first application a clone group with identical group records is
        // created and the elements fully inside the geometry move to it;
        // the value change is then applied to the clone only, so elements
        // outside the geometry (and the original group record) keep the old
        // value. Element internal state (NODE_DOF history) is untouched.
        if ( db( CHANGE_DATAITEM_GEOMETRY, ichange, change_dataitem_geometry,
              ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS ) ) {
          if ( strncmp( db_name(data_item_name), "group_", 6 ) )
            db_error( CHANGE_DATAITEM_GEOMETRY, ichange );
          {
            static long int *clone_map = NULL;
            static long int clone_map_len = 0;
            long int idat2=0, ielem2=0, max_element2=0, in_geometry2=0,
              all_in=0, new_group=0, max_group=0, length_clone=0,
              el_group2=0, length_nodes2=0, *node_in_geometry2=NULL,
              *ival_clone=NULL, el2[MNOL+1];
            double rdum2=0.;
            if ( found ) {
              db_max_index( GROUP_TYPE, max_group, VERSION_NORMAL, GET );
              if ( !clone_map ) {
                clone_map_len = 2*max_group + 64;
                clone_map = new long int[clone_map_len];
                for ( long int k=0; k<clone_map_len; k++ ) clone_map[k] = -1;
              }
              if ( data_item_index>=clone_map_len ) {
                long int new_len = 2*data_item_index + 64, *new_map =
                  new long int[new_len];
                for ( long int k=0; k<new_len; k++ )
                  new_map[k] = ( k<clone_map_len ? clone_map[k] : -1 );
                delete[] clone_map; clone_map = new_map;
                clone_map_len = new_len;
              }
              new_group = clone_map[data_item_index];
              if ( new_group<0 ) {
                new_group = max_group + 1;
                clone_map[data_item_index] = new_group;
                // clone all per-group records ("group_*") of the original
                ival_clone = get_new_int(DATA_ITEM_SIZE);
                for ( idat2=0; idat2<MDAT; idat2++ ) {
                  if ( !strncmp(db_name(idat2),"group_",6) &&
                       db_active_index(idat2, data_item_index,
                         VERSION_NORMAL) ) {
                    length_clone = db_len( idat2, data_item_index,
                      VERSION_NORMAL );
                    if ( db_type(idat2)==INTEGER ) {
                      db( idat2, data_item_index, ival_clone, ddum,
                        length_clone, VERSION_NORMAL, GET );
                      db( idat2, new_group, ival_clone, ddum,
                        length_clone, VERSION_NORMAL, PUT );
                    }
                    else {
                      db( idat2, data_item_index, idum, dval,
                        length_clone, VERSION_NORMAL, GET );
                      db( idat2, new_group, idum, dval,
                        length_clone, VERSION_NORMAL, PUT );
                    }
                  }
                }
                // move the elements fully inside the geometry to the clone
                node_in_geometry2 = get_new_int( 1+max_node );
                array_set( node_in_geometry2, 0, 1+max_node );
                for ( inod=0; inod<=max_node; inod++ ) {
                  if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
                    geometry( inod, ddum, change_dataitem_geometry,
                      in_geometry2, rdum2, ddum, rdum2, ddum,
                      NODE_START_REFINED, CHANGE_DATAITEM_GEOMETRY,
                      VERSION_NORMAL );
                    if ( in_geometry2 ) node_in_geometry2[inod] = 1;
                  }
                }
                db_max_index( ELEMENT, max_element2, VERSION_NORMAL, GET );
                for ( ielem2=0; ielem2<=max_element2; ielem2++ ) {
                  if ( db_active_index( ELEMENT, ielem2, VERSION_NORMAL ) ) {
                    el_group2 = -1;
                    db( ELEMENT_GROUP, ielem2, &el_group2, ddum, ldum,
                      VERSION_NORMAL, GET_IF_EXISTS );
                    if ( el_group2==data_item_index ) {
                      db( ELEMENT, ielem2, el2, ddum, length_nodes2,
                        VERSION_NORMAL, GET );
                      all_in = 1;
                      for ( long int k=1; k<length_nodes2; k++ )
                        if ( !node_in_geometry2[el2[k]] ) all_in = 0;
                      if ( all_in )
                        db( ELEMENT_GROUP, ielem2, &new_group, ddum, ldum,
                          VERSION_NORMAL, PUT );
                    }
                  }
                }
                delete[] node_in_geometry2;
              }
              // from now on the change applies to the clone group
              data_item_index = new_group;
            }
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
  }

  // control_reset_dof: reset the listed node dofs to a value that is
  // constant, or depends on another dof (via control_reset_value_dof +
  // control_reset_value_dof_diagram). Method -add or -multiply relative
  // to the current value, or -use to set the value.
  reset_dof = get_new_int(DATA_ITEM_SIZE);
  reset_value_dof = get_new_int(DATA_ITEM_SIZE);
  reset_value_diagram = get_new_dbl(DATA_ITEM_SIZE);
  db_max_index( CONTROL_RESET_DOF, max_reset, VERSION_NORMAL, GET );
  // control_reset_interface / control_reset_interface_strain (manual
  // Professional 6.354/6.355) may exist WITHOUT any control_reset_dof
  // (interface7 of the corpus). FIX (2026-09-03): the interface history
  // reset is scanned with db_active_index over the control range (its
  // db_max_index returns -1 although the records are active - measured
  // on interface7, record 15 active but max index -1), and it is no
  // longer nested inside the control_reset_dof gate, so it fires even
  // when only the interface reset record is present.
  if ( max_reset>=0 ) {
    swit = set_swit(-1,-1,"data");
    if ( swit ) pri( "In routine DATA (control_reset)" );
    db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET );
  }

  // control_reset_interface / control_reset_interface_strain: reset the
  // accumulated histories of the interface elements located in the
  // geometry. _interface resets ALL interface data (strains + tangential
  // forces); _interface_strain resets only the normal strain, keeping
  // the tangential force history (the stresses are "remembered": new
  // strains start at 0 and new stresses grow from the remembered ones
  // through the stiffness). Control records fire in the control step of
  // their own index (same rule as control_reset_dof, manual Professional
  // 6.350): interface7 resets at index 15 between the control 10/20
  // phases.
  {
    long int ireset_i = 0, iel_i = 0, max_element_i = 0,
      inol_i = 0, length_el_i = 0, all_in_i = 0, in_geom_i = 0,
      geometry_i[2], zero_one = 0;
    double rdum_i = 0.;
    db_max_index( ELEMENT, max_element_i, VERSION_NORMAL, GET );
    for ( ireset_i=0; ireset_i<=1000; ireset_i++ ) {
      long int full_reset = 0, strain_reset = 0;
      if ( db_active_index( CONTROL_RESET_INTERFACE, ireset_i,
           VERSION_NORMAL ) ) full_reset = 1;
      if ( db_active_index( CONTROL_RESET_INTERFACE_STRAIN, ireset_i,
           VERSION_NORMAL ) ) strain_reset = 1;
      if ( !full_reset && !strain_reset ) continue;
      if ( ireset_i!=icontrol ) continue;
      {
        long int geometry_src[2];
        if ( full_reset )
          db( CONTROL_RESET_INTERFACE, ireset_i, geometry_src, ddum,
            ldum, VERSION_NORMAL, GET );
        else
          db( CONTROL_RESET_INTERFACE_STRAIN, ireset_i, geometry_src, ddum,
            ldum, VERSION_NORMAL, GET );
        geometry_i[0] = geometry_src[0];
        geometry_i[1] = geometry_src[1];
      }
      for ( iel_i=0; iel_i<=max_element_i; iel_i++ ) {
        long int el_i[MNOL+1];
        if ( !db_active_index( ELEMENT, iel_i, VERSION_NORMAL ) )
          continue;
        db( ELEMENT, iel_i, el_i, ddum, length_el_i, VERSION_NORMAL, GET );
        all_in_i = 1;
        for ( inol_i=1; inol_i<length_el_i; inol_i++ ) {
          geometry( el_i[inol_i], ddum, geometry_i, in_geom_i, rdum_i,
            ddum, rdum_i, ddum, NODE_START_REFINED,
            CONTROL_RESET_INTERFACE, VERSION_NORMAL );
          if ( !in_geom_i ) all_in_i = 0;
        }
        if ( all_in_i ) {
          zero_one = 0;
          // per-integration-point histories (one value per facing pair,
          // ns1 = nnol/2); zero ALL slots in BOTH versions (FIX
          // 2026-09-03: the old PUT passed a leftover length and a
          // single value, so the reset never cleared the whole record
          // and the accumulated strain survived the reset).
          if ( strain_reset || full_reset ) {
            long int ns1_r = ( length_el_i-1 )/2;
            double zero_arr[4];
            for ( inol_i=0; inol_i<ns1_r && inol_i<4; inol_i++ )
              zero_arr[inol_i] = 0.;
            db( ELEMENT_INTERFACE_STRAIN_NORMAL, iel_i, idum, zero_arr,
              ns1_r, VERSION_NORMAL, PUT );
            db( ELEMENT_INTERFACE_STRAIN_NORMAL, iel_i, idum, zero_arr,
              ns1_r, VERSION_NEW, PUT );
          }
          if ( full_reset ) {
            long int ns1_r = ( length_el_i-1 )/2;
            double zero_arr[4];
            for ( inol_i=0; inol_i<ns1_r && inol_i<4; inol_i++ )
              zero_arr[inol_i] = 0.;
            db( ELEMENT_INTERFACE_FORCE_TANG, iel_i, idum, zero_arr,
              ns1_r, VERSION_NORMAL, PUT );
            db( ELEMENT_INTERFACE_FORCE_TANG, iel_i, idum, zero_arr,
              ns1_r, VERSION_NEW, PUT );
            db( ELEMENT_INTERFACE_FORCE_TANG2, iel_i, idum, zero_arr,
              ns1_r, VERSION_NORMAL, PUT );
            db( ELEMENT_INTERFACE_FORCE_TANG2, iel_i, idum, zero_arr,
              ns1_r, VERSION_NEW, PUT );
          }
        }
      }
    }
  }

    long int *reset_dof_node_filter = NULL, reset_dof_length = 0, idof_list = 0;
    for ( ireset=0; ireset<=max_reset; ireset++ ) {
      if ( db_active_index( CONTROL_RESET_DOF, ireset, VERSION_NORMAL ) ) {
        // manual Professional 6.350: control_reset_dof is indexed like
        // every control; it is applied only in the control steps of ITS
        // own index (icontrol==ireset). The corpus tests (hypo2/3,
        // mohrcou4, ...) reset -hyhis0/-epixx/-sigxx with indices 1..3
        // while the timestep runs at index 20. The Professional REJECTS
        // a reset sharing the timestep index ("Error detected for data
        // item : control_reset_dof, record : 20" - measured), so the
        // reset runs ONCE in its own (timeless) control step.
        if ( ireset!=icontrol ) continue;
        db( CONTROL_RESET_DOF, ireset, reset_dof, ddum, ldum, VERSION_NORMAL, GET );
        // the record may list SEVERAL dofs (-sigxx -sigyy -sigzz ...);
        // every listed dof gets the reset value.
        reset_dof_length = ldum;
        if ( reset_dof_length<1 ) db_error( CONTROL_RESET_DOF, ireset );
        for ( idof_list=0; idof_list<reset_dof_length; idof_list++ ) {
        idof_reset = reset_dof[idof_list];
        // materi_displacement_relative: a displacement reset re-synchronizes
        // the relative displacement reference (manual 4.13).
        if ( materi_displacement_relative ) {
          long int idof_reset_idx = idof_reset;
          if ( idof_reset_idx<0 ) {
            array_member( dof_label, idof_reset_idx, nuknwn, idof_reset_idx );
            if ( db_len( NODE_DOF, 1, VERSION_NORMAL )==npuknwn )
              idof_reset_idx /= nder;
          }
          if ( idof_reset_idx==dis_indx ) {
            for ( inod=0; inod<=max_node; inod++ ) {
              if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
                node_dof = db_dbl( NODE_DOF, inod, VERSION_NORMAL );
                for ( long int idim=0; idim<ndim; idim++ )
                  node_dof[dis_rel_indx+idim*nder] = 0.;
              }
            }
          }
        }
        reset_method = -USE;
        db( CONTROL_RESET_VALUE_METHOD, ireset, &reset_method, ddum, ldum,
          VERSION_NORMAL, GET_IF_EXISTS );

        // node selection filters (manual Professional 6.352/6.353/6.356):
        // control_reset_geometry (nodes of elements completely inside the
        // geometry), control_reset_node (nodes of elements with all their
        // nodes listed) and control_reset_element_group (restrict to
        // elements of the listed groups). Without any filter all nodes are
        // treated (previous behaviour).
        {
          long int use_filter = 0, length_rn = 0, length_reg = 0,
            iel2 = 0, inol2 = 0, max_element2 = 0, length_el2 = 0,
            all_in = 0, in_geometry2 = 0, all_nodes_listed = 0,
            el_group2 = 0, el2[MNOL+1];
          long int *reset_nodes_list = NULL;
          double rdum2 = 0.;
          reset_nodes_list = get_new_int(1+max_node);
          array_set( reset_nodes_list, 0, 1+max_node );
          if ( db_active_index( CONTROL_RESET_GEOMETRY, ireset,
              VERSION_NORMAL ) ||
               db_active_index( CONTROL_RESET_NODE, ireset,
              VERSION_NORMAL ) ||
               db_active_index( CONTROL_RESET_ELEMENT_GROUP, ireset,
              VERSION_NORMAL ) ) {
            long int reset_geometry[2] = {0,0};
            use_filter = 1;
            if ( db_active_index( CONTROL_RESET_GEOMETRY, ireset,
                VERSION_NORMAL ) )
              db( CONTROL_RESET_GEOMETRY, ireset, reset_geometry, ddum,
                ldum, VERSION_NORMAL, GET );
            if ( db_active_index( CONTROL_RESET_NODE, ireset,
                VERSION_NORMAL ) ) {
              db( CONTROL_RESET_NODE, ireset, reset_nodes_list, ddum,
                length_rn, VERSION_NORMAL, GET );
              // convert the listed nodes to a marker set (avoid clashing
              // with the output array usage below)
              for ( long int k=0; k<length_rn; k++ )
                reset_nodes_list[reset_nodes_list[k]] = 1;
            }
            db_max_index( ELEMENT, max_element2, VERSION_NORMAL, GET );
            for ( iel2=0; iel2<=max_element2; iel2++ ) {
              if ( !db_active_index( ELEMENT, iel2, VERSION_NORMAL ) )
                continue;
              if ( db_active_index( CONTROL_RESET_ELEMENT_GROUP, ireset,
                  VERSION_NORMAL ) ) {
                long int go = 0;
                long int *reg = NULL;
                reg = db_int( CONTROL_RESET_ELEMENT_GROUP, ireset,
                  VERSION_NORMAL );
                length_reg = db_len( CONTROL_RESET_ELEMENT_GROUP, ireset,
                  VERSION_NORMAL );
                el_group2 = -1;
                db( ELEMENT_GROUP, iel2, &el_group2, ddum, ldum,
                  VERSION_NORMAL, GET_IF_EXISTS );
                if ( array_member( reg, el_group2, length_reg, ldum ) )
                  go = 1;
                if ( !go ) continue;
              }
              db( ELEMENT, iel2, el2, ddum, length_el2, VERSION_NORMAL, GET );
              all_in = 1;
              all_nodes_listed = 1;
              for ( inol2=1; inol2<length_el2; inol2++ ) {
                if ( db_active_index( CONTROL_RESET_GEOMETRY, ireset,
                    VERSION_NORMAL ) ) {
                  geometry( el2[inol2], ddum, reset_geometry, in_geometry2,
                    rdum2, ddum, rdum2, ddum, NODE_START_REFINED,
                    CONTROL_RESET_GEOMETRY, VERSION_NORMAL );
                  if ( !in_geometry2 ) all_in = 0;
                }
                if ( db_active_index( CONTROL_RESET_NODE, ireset,
                    VERSION_NORMAL ) ) {
                  if ( !reset_nodes_list[el2[inol2]] )
                    all_nodes_listed = 0;
                }
              }
              if ( all_in && all_nodes_listed ) {
                for ( inol2=1; inol2<length_el2; inol2++ )
                  reset_nodes_list[el2[inol2]] = 2;
              }
            }
          }
          if ( use_filter ) {
            // marker 2 = selected; collapse to 0/1
            for ( inod=0; inod<=max_node; inod++ )
              reset_nodes_list[inod] = ( reset_nodes_list[inod]==2 );
          }
          else {
            for ( inod=0; inod<=max_node; inod++ )
              reset_nodes_list[inod] = 1;
          }
          reset_dof_node_filter = reset_nodes_list;
        }

        if ( db_active_index( CONTROL_RESET_VALUE_CONSTANT, ireset, VERSION_NORMAL ) ) {
          db( CONTROL_RESET_VALUE_CONSTANT, ireset, idum, &reset_value_constant,
            ldum, VERSION_NORMAL, GET );
          for ( inod=0; inod<=max_node; inod++ ) {
            if ( db_active_index( NODE, inod, VERSION_NORMAL ) &&
                 ( !reset_dof_node_filter || reset_dof_node_filter[inod] ) ) {
              node_dof = db_dbl( NODE_DOF, inod, VERSION_NORMAL );
              length = db_len( NODE_DOF, inod, VERSION_NORMAL );
              long int indx = idof_reset;
              if ( indx<0 ) {
                array_member( dof_label, indx, nuknwn, indx );
                if ( length==npuknwn ) indx /= nder;
              }
              if ( indx<0 || indx>length-1 )
                db_error( CONTROL_RESET_DOF, ireset );
              if      ( reset_method==-ADD ) node_dof[indx] += reset_value_constant;
              else if ( reset_method==-MULTIPLY ) node_dof[indx] *= reset_value_constant;
              else                            node_dof[indx] = reset_value_constant;
            }
          }
          // CONVERGENCE (2026-09-03, mohr_coul_direct4): a reset of a
          // normal stress dof (-sigxx/-sigyy/-sigzz) also initialises the
          // accumulated normal strain of the interface elements whose
          // normal aligns with the reset axis: sigma_n := reset_value ->
          // ELEMENT_INTERFACE_STRAIN_NORMAL := reset_value/kn. Verified
          // against the Professional .dbs: mohr_coul_direct4 resets
          // -sigyy to -1 and the horizontal interface carries
          // sigma_n = -1 (eps_n = -1e-6 with kn = 1e6) through the whole
          // run, so the yield limit sees c + |Fn|*tan(phi) = 1.20271
          // instead of the unconfined cohesion c.
          if ( reset_method!=-ADD && reset_method!=-MULTIPLY ) {
            long int reset_axis = -1;
            const char *reset_name = db_name( labs( idof_reset ) );
            if ( reset_name && !strcmp( reset_name, "sigxx" ) ) reset_axis = 0;
            else if ( reset_name && !strcmp( reset_name, "sigyy" ) ) reset_axis = 1;
            else if ( reset_name && !strcmp( reset_name, "sigzz" ) ) reset_axis = 2;
            if ( reset_axis>=0 && reset_axis<ndim ) {
              long int max_element_r = 0, iel_r = 0, length_el_r = 0,
                element_group_r = 0, inol_r = 0;
              double ddum_r[1], kn_r = 0., normal_r[MDIM];
              db_max_index( ELEMENT, max_element_r, VERSION_NORMAL, GET );
              for ( iel_r=0; iel_r<=max_element_r; iel_r++ ) {
                if ( !db_active_index( ELEMENT, iel_r, VERSION_NORMAL ) ) continue;
                element_group_r = 0;
                db( ELEMENT_GROUP, iel_r, &element_group_r, ddum_r, ldum,
                  VERSION_NORMAL, GET_IF_EXISTS );
                if ( !db_active_index( GROUP_INTERFACE, element_group_r,
                    VERSION_NORMAL ) ) continue;
                long int el_r[MNOL+1];
                db( ELEMENT, iel_r, el_r, ddum_r, length_el_r, VERSION_NORMAL, GET );
                if ( length_el_r<3 ) continue;
                // normal of the interface plane from the reference
                // geometry (NODE_START_REFINED, total_linear frame of
                // interface_element()): 2D normal perpendicular to the
                // side-1 edge; 3D normal = cross product of the side-1
                // edges.
                {
                  double *ca = db_dbl( NODE_START_REFINED, el_r[1], VERSION_NORMAL );
                  double *cb = db_dbl( NODE_START_REFINED, el_r[2], VERSION_NORMAL );
                  if ( ndim==2 ) {
                    normal_r[0] = -( cb[1] - ca[1] );
                    normal_r[1] =    cb[0] - ca[0];
                  }
                  else if ( length_el_r>=4 ) {
                    double *cc = db_dbl( NODE_START_REFINED, el_r[3], VERSION_NORMAL );
                    double e1[MDIM], e2[MDIM];
                    for ( long int d=0; d<3; d++ ) {
                      e1[d] = cb[d] - ca[d];
                      e2[d] = cc[d] - ca[d];
                    }
                    normal_r[0] = e1[1]*e2[2] - e1[2]*e2[1];
                    normal_r[1] = e1[2]*e2[0] - e1[0]*e2[2];
                    normal_r[2] = e1[0]*e2[1] - e1[1]*e2[0];
                  }
                  else continue;
                }
                {
                  double dot = 0.;
                  for ( long int d=0; d<ndim; d++ ) dot += normal_r[d]*normal_r[d];
                  if ( dot<1.e-24 ) continue;
                  dot = normal_r[reset_axis]*normal_r[reset_axis]/dot;
                  if ( dot<0.99 ) continue;
                }
                {
                  double ddum3_r[3];
                  array_set( ddum3_r, 0., 3 );
                  db( GROUP_INTERFACE_MATERI_ELASTI_STIFFNESS, element_group_r,
                    idum, ddum3_r, ldum, VERSION_NORMAL, GET_IF_EXISTS );
                  kn_r = ddum3_r[0];
                }
                if ( kn_r<=0. ) continue;
                // all nodes of the element inside the reset filter (when
                // a geometry/node/element_group filter is active)
                if ( reset_dof_node_filter ) {
                  long int all_in_r = 1;
                  for ( inol_r=1; inol_r<length_el_r; inol_r++ ) {
                    if ( !reset_dof_node_filter[el_r[inol_r]] ) all_in_r = 0;
                  }
                  if ( !all_in_r ) continue;
                }
                {
                  long int ns1_r = ( length_el_r-1 )/2;
                  double strain_r[4];
                  for ( inol_r=0; inol_r<ns1_r; inol_r++ )
                    strain_r[inol_r] = reset_value_constant / kn_r;
                  ldum = ns1_r;
                  db( ELEMENT_INTERFACE_STRAIN_NORMAL, iel_r, idum, strain_r,
                    ldum, VERSION_NORMAL, PUT );
                  db( ELEMENT_INTERFACE_STRAIN_NORMAL, iel_r, idum, strain_r,
                    ldum, VERSION_NEW, PUT );
                }
              }
            }
          }
        }
        else if ( db_active_index( CONTROL_RESET_VALUE_DOF, ireset, VERSION_NORMAL ) ) {
          db( CONTROL_RESET_VALUE_DOF, ireset, &idof_value, ddum, ldum,
            VERSION_NORMAL, GET );
          db( CONTROL_RESET_VALUE_DOF_DIAGRAM, ireset, idum, reset_value_diagram,
            length_diagram, VERSION_NORMAL, GET );
          for ( inod=0; inod<=max_node; inod++ ) {
            if ( db_active_index( NODE, inod, VERSION_NORMAL ) &&
                 ( !reset_dof_node_filter || reset_dof_node_filter[inod] ) ) {
              node_dof = db_dbl( NODE_DOF, inod, VERSION_NORMAL );
              length = db_len( NODE_DOF, inod, VERSION_NORMAL );
              long int indx_reset = idof_reset, indx_val = idof_value;
              if ( indx_reset<0 ) {
                array_member( dof_label, indx_reset, nuknwn, indx_reset );
                if ( length==npuknwn ) indx_reset /= nder;
              }
              if ( indx_val<0 ) {
                array_member( dof_label, indx_val, nuknwn, indx_val );
                if ( length==npuknwn ) indx_val /= nder;
              }
              if ( indx_reset<0 || indx_reset>length-1 )
                db_error( CONTROL_RESET_DOF, ireset );
              if ( indx_val<0 || indx_val>length-1 )
                db_error( CONTROL_RESET_DOF, ireset );
               table_xy( reset_value_diagram, "CONTROL_RESET_VALUE_DOF_DIAGRAM",
                 length_diagram, node_dof[indx_val], val );
               if      ( reset_method==-ADD ) node_dof[indx_reset] += val;
               else if ( reset_method==-MULTIPLY ) node_dof[indx_reset] *= val;
               else                             node_dof[indx_reset] = val;
             }
           }
         }
         // spatial distributions: the reset value depends on the node
         // coordinates (x, y, z). Variants: _linear, _exponent, _power,
         // _square_root, _logarithmic, _logarithmic_second, _multi_linear.
         else if ( db_active_index( CONTROL_RESET_VALUE_LINEAR, ireset, VERSION_NORMAL ) ||
             db_active_index( CONTROL_RESET_VALUE_EXPONENT, ireset, VERSION_NORMAL ) ||
             db_active_index( CONTROL_RESET_VALUE_POWER, ireset, VERSION_NORMAL ) ||
             db_active_index( CONTROL_RESET_VALUE_SQUARE_ROOT, ireset, VERSION_NORMAL ) ||
             db_active_index( CONTROL_RESET_VALUE_LOGARITHMIC, ireset, VERSION_NORMAL ) ||
             db_active_index( CONTROL_RESET_VALUE_LOGARITHMIC_SECOND, ireset, VERSION_NORMAL ) ||
             db_active_index( CONTROL_RESET_VALUE_MULTI_LINEAR, ireset, VERSION_NORMAL ) ) {
           long int spatial_data[7] = { CONTROL_RESET_VALUE_LINEAR,
             CONTROL_RESET_VALUE_EXPONENT, CONTROL_RESET_VALUE_POWER,
             CONTROL_RESET_VALUE_SQUARE_ROOT, CONTROL_RESET_VALUE_LOGARITHMIC,
             CONTROL_RESET_VALUE_LOGARITHMIC_SECOND, CONTROL_RESET_VALUE_MULTI_LINEAR };
           long int sdat = 0, spatial_active = -1;
           for ( sdat=0; sdat<7; sdat++ )
             if ( db_active_index( spatial_data[sdat], ireset, VERSION_NORMAL ) )
               spatial_active = spatial_data[sdat];
           long int nl = 0;
           if ( spatial_active==CONTROL_RESET_VALUE_LINEAR ) nl = ndim;
           else if ( spatial_active==CONTROL_RESET_VALUE_POWER ) nl = 2*ndim;
           else if ( spatial_active==CONTROL_RESET_VALUE_SQUARE_ROOT ) nl = 3*ndim;
           else if ( spatial_active==CONTROL_RESET_VALUE_EXPONENT ) nl = 5*ndim;
           else if ( spatial_active==CONTROL_RESET_VALUE_LOGARITHMIC ) nl = 5*ndim;
           else if ( spatial_active==CONTROL_RESET_VALUE_LOGARITHMIC_SECOND ) nl = 7*ndim;
           for ( inod=0; inod<=max_node; inod++ ) {
             if ( db_active_index( NODE, inod, VERSION_NORMAL ) &&
                  ( !reset_dof_node_filter || reset_dof_node_filter[inod] ) ) {
               node_dof = db_dbl( NODE_DOF, inod, VERSION_NORMAL );
               length = db_len( NODE_DOF, inod, VERSION_NORMAL );
               long int indx_reset = idof_reset;
               if ( indx_reset<0 ) {
                 array_member( dof_label, indx_reset, nuknwn, indx_reset );
                 if ( length==npuknwn ) indx_reset /= nder;
               }
               if ( indx_reset<0 || indx_reset>length-1 )
                 db_error( CONTROL_RESET_DOF, ireset );
               coord = db_dbl( NODE, inod, VERSION_NORMAL );
               for ( idim=0; idim<ndim; idim++ ) coords[idim] = coord[idim];
               db( spatial_active, ireset, idum, reset_value_diagram, ldum,
                 VERSION_NORMAL, GET );
               val = 0.;
               if ( spatial_active==CONTROL_RESET_VALUE_LINEAR ) {
                 // ax x + ay y + az z
                 for ( idim=0; idim<ndim; idim++ )
                   val += reset_value_diagram[idim] * coords[idim];
               }
               else if ( spatial_active==CONTROL_RESET_VALUE_POWER ) {
                 // ax x^bx + ay y^by + az z^bz
                 for ( idim=0; idim<ndim; idim++ )
                   val += reset_value_diagram[2*idim] *
                     scalar_power( coords[idim], reset_value_diagram[2*idim+1] );
               }
               else if ( spatial_active==CONTROL_RESET_VALUE_SQUARE_ROOT ) {
                 // ax sqrt( bx x + cx x^2 ) ... (3 coefs per dim)
                 for ( idim=0; idim<ndim; idim++ ) {
                   double b = reset_value_diagram[3*idim+1];
                   double c = reset_value_diagram[3*idim+2];
                   val += reset_value_diagram[3*idim] *
                     sqrt( scalar_dabs( b*coords[idim] + c*coords[idim]*coords[idim] ) );
                 }
               }
               else if ( spatial_active==CONTROL_RESET_VALUE_EXPONENT ) {
                 // ax e^(bx + cx x dx + ex x) (5 coefs per dim)
                 for ( idim=0; idim<ndim; idim++ ) {
                   double a = reset_value_diagram[5*idim+0];
                   double b = reset_value_diagram[5*idim+1];
                   double c = reset_value_diagram[5*idim+2];
                   double d = reset_value_diagram[5*idim+3];
                   double e = reset_value_diagram[5*idim+4];
                   val += a * exp( b + c * coords[idim] * d + e * coords[idim] );
                 }
               }
               else if ( spatial_active==CONTROL_RESET_VALUE_LOGARITHMIC ) {
                 // ax ln( bx + cx x dx + ex x ) (5 coefs per dim)
                 for ( idim=0; idim<ndim; idim++ ) {
                   double a = reset_value_diagram[5*idim+0];
                   double b = reset_value_diagram[5*idim+1];
                   double c = reset_value_diagram[5*idim+2];
                   double d = reset_value_diagram[5*idim+3];
                   double e = reset_value_diagram[5*idim+4];
                   val += a * log( scalar_dabs( b + c * coords[idim] * d + e * coords[idim] ) );
                 }
               }
               else if ( spatial_active==CONTROL_RESET_VALUE_LOGARITHMIC_SECOND ) {
                 // (ax + bx) e^(cx ln(dx (x+ex)/fx)) + gx  (7 coefs per dim)
                 for ( idim=0; idim<ndim; idim++ ) {
                   double a = reset_value_diagram[7*idim+0];
                   double b = reset_value_diagram[7*idim+1];
                   double c = reset_value_diagram[7*idim+2];
                   double d = reset_value_diagram[7*idim+3];
                   double e = reset_value_diagram[7*idim+4];
                   double f = reset_value_diagram[7*idim+5];
                   double g = reset_value_diagram[7*idim+6];
                   val += (a+b) * exp( c * log( scalar_dabs( d*(coords[idim]+e)/f ) ) ) + g;
                 }
               }
               else if ( spatial_active==CONTROL_RESET_VALUE_MULTI_LINEAR ) {
                 // table (z0 value0 z1 value1 ...); vertical coordinate:
                 // 1D -> x, 2D -> y, 3D -> z
                 long int vdim = ( ndim==1 ) ? 0 : ( ndim==2 ) ? 1 : 2;
                 table_xy( reset_value_diagram, "CONTROL_RESET_VALUE_MULTI_LINEAR",
                   nl/2, coords[vdim], val );
               }
                if      ( reset_method==-ADD ) node_dof[indx_reset] += val;
                else if ( reset_method==-MULTIPLY ) node_dof[indx_reset] *= val;
                else                             node_dof[indx_reset] = val;
              }
            }
          }
        }
        }   // end dof list loop
      }
      delete[] reset_dof;
      delete[] reset_value_dof;
      delete[] reset_value_diagram;

  delete[] dof_label;
  delete[] integer_range;
  delete[] change_dataitem_time;
  delete[] dval;

}
