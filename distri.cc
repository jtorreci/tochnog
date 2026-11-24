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

void distribute( void )

{
  int distribute_ran=-1;
  long int distribution_type=0, data_item_name=0, data_item_number=0, length=0,
    icontrol=0, index=0, max_index=0, group_data=0, name=0,
    length_control_distribute=0,
    length_control_distribute_values=0,
    length_element_distribute=0,
    length_element_distribute_values=0,
    ldum=0, idistribute=0, ndistribute=0,
    idum[1], element_distribute[2],
    *control_distribute=NULL, *dof_label=NULL;
  double range=0., variance=0., ran=0., tmp=0., ddum[1],
    *element_distribute_values=NULL, *control_distribute_values=NULL,
    *dval=NULL;
  char str[MCHAR];

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  if ( db_active_index( CONTROL_DISTRIBUTE, icontrol, VERSION_NORMAL ) ) {
    control_distribute = get_new_int(DATA_ITEM_SIZE);
    dof_label = get_new_int(MUKNWN);
    element_distribute_values = get_new_dbl(DATA_ITEM_SIZE);
    control_distribute_values = get_new_dbl(DATA_ITEM_SIZE);
    db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    if ( db_active_index( CONTROL_DISTRIBUTE_VALUES, icontrol, VERSION_NORMAL ) )
      db( CONTROL_DISTRIBUTE_VALUES, icontrol, idum, control_distribute_values,
        length_control_distribute_values, VERSION_NORMAL, GET );
    else
      length_control_distribute_values = 0;
    db( CONTROL_DISTRIBUTE, icontrol, control_distribute, ddum,
      length_control_distribute, VERSION_NORMAL, GET );

    // Professional layout (manual 6.127): one distribution per record
    //   [ distribution_type, data_item_name, data_item_index(-all|n),
    //     data_item_number ]
    // with mean/std from control_distribute_parameters, optional seed,
    // spatial correlation and minimum/maximum clamping. The GNU legacy
    // layout is triplets (type, name, number) x N with per-record deltas
    // in control_distribute_values; it is kept unchanged below.
    if ( length_control_distribute==4 ) {

      long int data_item_index=0, use_corr=0, use_minmax=0, nent=0, ient=0,
        jent=0, max_element=0, nnol=0, element_group=0, length_corr=0,
        length_nodes=0, el[MNOL+1], constant_field=0, have_coords=0;
      double mean=0., std=0., mu_ln=0., sigma_ln=0., z=0.,
        corr_length[MDIM], corr_distance=0., minmax[2], coord[MDIM],
        seed_value=0.;
      long int *entities=NULL;
      double *zval=NULL, *ecoord=NULL, *wgt=NULL, *val=NULL;

      distribution_type = control_distribute[0];
      data_item_name    = control_distribute[1];
      data_item_index   = control_distribute[2];
      data_item_number  = control_distribute[3];
      if ( distribution_type!=-NORMAL && distribution_type!=-LOGNORMAL )
        db_error( CONTROL_DISTRIBUTE, icontrol );

      // control_distribute_seed: negative idum reinitializes the random
      // sequence deterministically (numerical recipes ran1), so equal
      // seeds give reproducible identical fields.
      if ( db( CONTROL_DISTRIBUTE_SEED, icontrol, idum, &seed_value,
          ldum, VERSION_NORMAL, GET_IF_EXISTS ) )
        distribute_ran = (int)( -(labs((long int)seed_value)) - 1 );

      // control_distribute_parameters: mean value and standard deviation
      // of the distributed value itself.
      db( CONTROL_DISTRIBUTE_PARAMETERS, icontrol, idum, minmax,
        ldum, VERSION_NORMAL, GET );
      mean = minmax[0];
      std  = minmax[1];
      if ( std<0. ) db_error( CONTROL_DISTRIBUTE_PARAMETERS, icontrol );
      if ( distribution_type==-LOGNORMAL ) {
        if ( mean<=0. ) db_error( CONTROL_DISTRIBUTE_PARAMETERS, icontrol );
        // lognormal with mean/std of the variable itself (not of the
        // underlying normal): mu_ln and sigma_ln of ln(X)
        sigma_ln = sqrt( log( 1. + (std/mean)*(std/mean) ) );
        mu_ln    = log(mean) - 0.5*sigma_ln*sigma_ln;
      }

      // spatial correlation (manual 6.128/6.129): correlation length per
      // direction (1..ndim values); > 1.e12 means a constant field.
      // Data are only correlated below the correlation distance
      // (default 4 times the correlation length).
      array_set( corr_length, 0., MDIM );
      length_corr = 0;
      use_corr = db( CONTROL_DISTRIBUTE_CORRELATION_LENGTH, icontrol,
        idum, corr_length, length_corr, VERSION_NORMAL, GET_IF_EXISTS );
      if ( use_corr ) {
        long int idim=0;
        double lmax=0.;
        if ( length_corr<1 || length_corr>ndim )
          db_error( CONTROL_DISTRIBUTE_CORRELATION_LENGTH, icontrol );
        for ( idim=length_corr; idim<ndim; idim++ )
          corr_length[idim] = corr_length[0];
        for ( idim=0; idim<ndim; idim++ )
          if ( corr_length[idim]>lmax ) lmax = corr_length[idim];
        corr_distance = 4. * lmax;
        db( CONTROL_DISTRIBUTE_CORRELATION_DISTANCE, icontrol, idum,
          &corr_distance, ldum, VERSION_NORMAL, GET_IF_EXISTS );
        if ( corr_length[0]>1.e12 ) constant_field = 1;
      }

      // minimum/maximum clamp on the drawn value
      minmax[0] = -DBL_MAX; minmax[1] = DBL_MAX;
      use_minmax = db( CONTROL_DISTRIBUTE_MINIMUM_MAXIMUM, icontrol,
        idum, minmax, ldum, VERSION_NORMAL, GET_IF_EXISTS );

      strcpy( str, db_name( data_item_name ) );
      group_data = ( str[0]=='g' && str[1]=='r' && str[2]=='o' && str[3]=='u' &&
        str[4]=='p' );

      // collect the target entities
      db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
      if ( max_element>max_index ) max_index = max_element;
      entities = get_new_int(max_index+2);
      if ( group_data ) {
        // elements which use the group record: with a specific index the
        // elements of that element group, with -all the elements whose
        // group has the record active
        db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
        for ( index=0; index<=max_element; index++ ) {
          if ( !db_active_index( ELEMENT, index, VERSION_NORMAL ) ) continue;
          element_group = -1;
          db( ELEMENT_GROUP, index, &element_group, ddum, ldum,
            VERSION_NORMAL, GET_IF_EXISTS );
          if ( data_item_index==-ALL ) {
            if ( !db_active_index( data_item_name, element_group,
                VERSION_NORMAL ) ) continue;
          }
          else if ( element_group!=data_item_index ) continue;
          entities[nent++] = index;
        }
      }
      else {
        db_max_index( data_item_name, max_index, VERSION_NORMAL, GET );
        for ( index=0; index<=max_index; index++ )
          if ( db_active_index( data_item_name, index, VERSION_NORMAL ) )
            entities[nent++] = index;
        // de-activate on boundary nodes is a GNU-legacy behaviour of the
        // node path; the Professional path does not mention it.
      }

      // draw standard normals per entity (correlated in space when asked)
      zval    = get_new_dbl(nent+1);
      val     = get_new_dbl(nent+1);
      ecoord  = get_new_dbl((nent+1)*ndim);
      have_coords = 0;
      if ( group_data ) {
        // element centroids
        for ( ient=0; ient<nent; ient++ ) {
          db( ELEMENT, entities[ient], el, ddum, length_nodes,
            VERSION_NORMAL, GET );
          nnol = length_nodes - 1;
          array_set( coord, 0., MDIM );
          for ( index=1; index<length_nodes; index++ ) {
            dval = db_dbl( NODE, el[index], VERSION_NORMAL );
            for ( long int idim=0; idim<ndim; idim++ )
              coord[idim] += dval[idim]/nnol;
          }
          for ( long int idim=0; idim<ndim; idim++ )
            ecoord[ient*ndim+idim] = coord[idim];
        }
        have_coords = 1;
      }
      else if ( db_data_class(data_item_name)==NODE ) {
        for ( ient=0; ient<nent; ient++ ) {
          dval = db_dbl( NODE, entities[ient], VERSION_NORMAL );
          for ( long int idim=0; idim<ndim; idim++ )
            ecoord[ient*ndim+idim] = dval[idim];
        }
        have_coords = 1;
      }

      for ( ient=0; ient<nent; ient++ ) {
        zval[ient] = scalar_ran_normal( distribute_ran );
        if ( constant_field && ient>0 ) zval[ient] = zval[0];
      }

      if ( use_corr && !constant_field && have_coords && nent>1 ) {
        // correlated field: v_i = sum_j w_ij z_j / sqrt(sum_j w_ij^2)
        // with w_ij = exp(-d_eff), d_eff the distance scaled per direction
        // by the correlation lengths; pairs beyond the correlation
        // distance (physical) get zero weight. The normalization keeps
        // unit variance.
        wgt = get_new_dbl(nent+1);
        for ( ient=0; ient<nent; ient++ ) {
          double wsum=0., w2=0.;
          array_set( wgt, 0., nent+1 );
          for ( jent=0; jent<nent; jent++ ) {
            double d2=0., dphys=0., deff=0.;
            for ( long int idim=0; idim<ndim; idim++ ) {
              double d = ecoord[ient*ndim+idim]-ecoord[jent*ndim+idim];
              dphys += d*d;
              if ( corr_length[idim]>TINY )
                d2 += (d/corr_length[idim])*(d/corr_length[idim]);
            }
            dphys = sqrt(dphys);
            if ( dphys>corr_distance ) continue;
            deff = ( d2>TINY ? sqrt(d2) : 0. );
            wgt[jent] = exp(-deff);
          }
          for ( jent=0; jent<nent; jent++ ) { wsum += wgt[jent]; w2 += wgt[jent]*wgt[jent]; }
          if ( w2>TINY ) {
            double v = 0.;
            for ( jent=0; jent<nent; jent++ ) v += wgt[jent]*zval[jent];
            zval[ient] = v / sqrt(w2);
          }
        }
        delete[] wgt;
      }

      // transform to the distributed value and clamp
      for ( ient=0; ient<nent; ient++ ) {
        z = zval[ient];
        if      ( distribution_type==-NORMAL )
          val[ient] = mean + std*z;
        else
          val[ient] = exp( mu_ln + sigma_ln*z );
        if ( use_minmax ) {
          if ( val[ient]<minmax[0] ) val[ient] = minmax[0];
          if ( val[ient]>minmax[1] ) val[ient] = minmax[1];
        }
      }

      if ( group_data ) {
        // per-element delta so that the += in get_group_data yields
        // exactly the drawn value (manual: the record itself is not
        // changed, the item is changed for the elements using it).
        // distribute() runs at every step start: refresh the per-element
        // records (redraw per step, like the GNU legacy path).
        db_delete( ELEMENT_DISTRIBUTE, VERSION_NORMAL );
        db_delete( ELEMENT_DISTRIBUTE_VALUES, VERSION_NORMAL );
        double record_value = 0.;
        for ( ient=0; ient<nent; ient++ ) {
          element_group = -1;
          db( ELEMENT_GROUP, entities[ient], &element_group, ddum, ldum,
            VERSION_NORMAL, GET_IF_EXISTS );
          length = db_len( data_item_name, element_group, VERSION_NORMAL );
          if ( data_item_number<0 || data_item_number>length-1 )
            db_error( CONTROL_DISTRIBUTE, icontrol );
          dval = db_dbl( data_item_name, element_group, VERSION_NORMAL );
          record_value = dval[data_item_number];
          tmp = val[ient] - record_value;
          if ( db_active_index( ELEMENT_DISTRIBUTE, entities[ient],
              VERSION_NORMAL ) ) {
            db( ELEMENT_DISTRIBUTE, entities[ient], element_distribute, ddum,
              length_element_distribute, VERSION_NORMAL, GET );
            db( ELEMENT_DISTRIBUTE_VALUES, entities[ient], idum,
              element_distribute_values, length_element_distribute_values,
              VERSION_NORMAL, GET );
          }
          else {
            length_element_distribute = 0;
            length_element_distribute_values = 0;
          }
          element_distribute[length_element_distribute*2+0] = data_item_name;
          element_distribute[length_element_distribute*2+1] = data_item_number;
          element_distribute_values[length_element_distribute_values] = tmp;
          length_element_distribute += 2;
          length_element_distribute_values += 1;
          db( ELEMENT_DISTRIBUTE, entities[ient], element_distribute, ddum,
            length_element_distribute, VERSION_NORMAL, PUT );
          db( ELEMENT_DISTRIBUTE_VALUES, entities[ient], idum,
            element_distribute_values, length_element_distribute_values,
            VERSION_NORMAL, PUT );
        }
      }
      else {
        // node-like items: the drawn value replaces the selected number
        // of the record (manual examples: nodal temperatures, y
        // coordinates)
        for ( ient=0; ient<nent; ient++ ) {
          dval = db_dbl( data_item_name, entities[ient], VERSION_NORMAL );
          length = db_len( data_item_name, entities[ient], VERSION_NORMAL );
          if ( data_item_number<0 ) {
            array_member( dof_label, data_item_number, nuknwn,
              data_item_number );
            if ( length==npuknwn ) data_item_number /= nder;
          }
          if ( data_item_number<0 || data_item_number>length-1 )
            db_error( CONTROL_DISTRIBUTE, icontrol );
          dval[data_item_number] = val[ient];
        }
      }

      delete[] entities; delete[] zval; delete[] val; delete[] ecoord;

    }
    else {

      // GNU legacy layout: triplets (type, name, number) x N with deltas
      if ( length_control_distribute!=3*length_control_distribute_values ) {
        pri( "Error: lengths of CONTROL_DISTRIBUTE and CONTROL_DISTRIBUTE_VALUES do not match." );
        exit(TN_EXIT_STATUS);
      }
      db_delete( ELEMENT_DISTRIBUTE, VERSION_NORMAL );
      db_delete( ELEMENT_DISTRIBUTE_VALUES, VERSION_NORMAL );
      ndistribute = length_control_distribute_values;
      for ( idistribute=0; idistribute<ndistribute; idistribute++ ) {
        distribution_type = control_distribute[idistribute*3+0];
        data_item_name = control_distribute[idistribute*3+1];
        data_item_number = control_distribute[idistribute*3+2];
        strcpy( str, db_name( data_item_name ) );
        group_data = ( str[0]=='g' && str[1]=='r' && str[2]=='o' && str[3]=='u' &&
          str[4]=='p' );
        if ( group_data )
          name = -ELEMENT;
        else
          name = data_item_name;
        if      ( distribution_type==-UNIFORM )
          range = control_distribute_values[idistribute];
        else if ( distribution_type==-NORMAL )
          variance = control_distribute_values[idistribute];
        else
          db_error(  CONTROL_DISTRIBUTE, icontrol );
        distribute_ran = -1;
        db_max_index( name, max_index, VERSION_NORMAL, GET );
        for ( index=0; index<=max_index; index++ ) {
          if ( db_active_index( name, index, VERSION_NORMAL ) ) {
            if ( distribution_type==-UNIFORM ) {
              ran = scalar_ran_uniform( distribute_ran );
              tmp = range*ran - 0.5*range;
            }
            else if ( distribution_type==-NORMAL ) {
              ran = scalar_ran_normal( distribute_ran );
              tmp = variance*ran;
            }
            else
              db_error(  CONTROL_DISTRIBUTE, icontrol );
            if ( name==-ELEMENT ) {
              if ( db_active_index( ELEMENT_DISTRIBUTE, index, VERSION_NORMAL ) ) {
                db( ELEMENT_DISTRIBUTE, index, element_distribute, ddum,
                  length_element_distribute, VERSION_NORMAL, GET );
                db( ELEMENT_DISTRIBUTE_VALUES, index, idum, element_distribute_values,
                  length_element_distribute_values, VERSION_NORMAL, GET );
              }
              else {
                length_element_distribute = 0;
                length_element_distribute_values = 0;
              }
              element_distribute[length_element_distribute*2+0] = data_item_name;
              element_distribute[length_element_distribute*2+1] = data_item_number;
              element_distribute_values[length_element_distribute_values] = tmp;
              length_element_distribute += 2;
              length_element_distribute_values += 1;
              db( ELEMENT_DISTRIBUTE, index, element_distribute, ddum,
                length_element_distribute, VERSION_NORMAL, PUT );
              db( ELEMENT_DISTRIBUTE_VALUES, index, idum, element_distribute_values,
                length_element_distribute_values, VERSION_NORMAL, PUT );
            }
            else {
              if ( data_item_name==-NODE ) {
                if ( db_active_index( NODE_BOUNDARY, index, VERSION_NORMAL ) ) tmp = 0.;
              }
              dval = db_dbl( data_item_name, index, VERSION_NORMAL );
              length = db_len( data_item_name, index, VERSION_NORMAL );
              data_item_number = control_distribute[idistribute*3+2];
              if ( data_item_number<0 ) {
                array_member(dof_label,data_item_number,nuknwn,data_item_number);
                if ( db_len(data_item_name,index,VERSION_NORMAL)==npuknwn )
                  data_item_number /= nder;
              }
              if ( data_item_number<0 || data_item_number>length-1 )
                db_error(  CONTROL_DISTRIBUTE, icontrol );
              dval[data_item_number] += tmp;
            }
          }
        }
      }
    }
    delete[] control_distribute;
    delete[] dof_label;
    delete[] element_distribute_values;
    delete[] control_distribute_values;
  }
}
