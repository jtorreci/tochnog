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

#define DEFAULT_VARIATION 1.e-2
#define EPS_PARAMETER 1.e-10

void inverse_calculation( long int ipar, long int npar, long int ipar_i,
  long int ipar_n, long int max_time, long int task )

  /* Numerical central differences are used to determine the sensitivities
     of the targets with respect to the parameters 
  */

{

  long int itar=0, mtar=0, data_item_name=0, swit=0,
    data_item_index=0, data_item_number=0, length=0,
    inverse_target_timestep=0, length_inverse_target_data=0,
    icontrol=0, ldum=0, idum[1], *inverse_parameter=NULL, 
    *inverse_target=NULL, dof_label[MUKNWN];
  double value=0., penalty=0., difference=0.,
    lower_value=0., higher_value=0., inverse_parameter_step=0.,
    inverse_parameter_variation=0., delta=0., tmp=0., dfdx=0., d2fdx2=0.,
    inverse_history_step=1., inverse_history_sum=0.,
    ddum[1], inverse_history[2], inverse_parameter_limits[2], 
    inverse_parameter_sensitivity[3], *actual_values=NULL, 
    *inverse_target_data=NULL, *parameter=NULL;

  if ( npar>0 ) {
    db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    swit = set_swit(-1,-1,"inverse_calculation");
    if ( task==INVERSE_PARAMETER_SET ) {
      if ( db_active_index( INVERSE_PARAMETER_VARIATION, ipar, VERSION_NORMAL ) )
        db( INVERSE_PARAMETER_VARIATION, ipar, idum, 
          &inverse_parameter_variation, ldum, VERSION_NORMAL, GET );
      else
        inverse_parameter_variation = DEFAULT_VARIATION;
      inverse_parameter = db_int( INVERSE_PARAMETER, ipar, VERSION_NORMAL );
      data_item_name = inverse_parameter[0];
      data_item_index = inverse_parameter[1];
      data_item_number = inverse_parameter[2];
      if ( data_item_number<0 ) {
        array_member(dof_label,data_item_number,nuknwn,data_item_number);
        if ( db_len(data_item_name,data_item_index,VERSION_NORMAL)==npuknwn ) 
          data_item_number /= nder;
      }
      length = db_len(data_item_name,data_item_index, VERSION_NORMAL);
      if ( data_item_number<0 || data_item_number>=length ) 
        db_error( INVERSE_PARAMETER, ipar );
      parameter = db_dbl( data_item_name, data_item_index, VERSION_NORMAL );
        // Set parameter for central differences.
      if ( scalar_dabs( parameter[data_item_number] ) > EPS_PARAMETER ) {
        tmp = inverse_parameter_variation;
        if      ( ipar_i==0 ) {
          parameter[data_item_number] *= (1.-tmp);
        }
        else if ( ipar_i==1 ) {
          parameter[data_item_number] /= (1.-tmp);
          parameter[data_item_number] *= (1.+tmp);
        }
        else {
          assert( ipar_i==ipar_n-1 );
          parameter[data_item_number] /= (1.+tmp);
        }
      }
      if ( db_data_class(data_item_name)==TENDON ) tendon_distribute();
    }
    if ( task==INVERSE_DETERMINE_SENSITIVITY ) {
      db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
      if ( db_active_index( INVERSE_PARAMETER_SENSITIVITY, ipar, VERSION_NORMAL ) )
        db( INVERSE_PARAMETER_SENSITIVITY, ipar, idum, 
          inverse_parameter_sensitivity, ldum, VERSION_NORMAL, GET_IF_EXISTS );
      else
        array_set( inverse_parameter_sensitivity, 0., ipar_n );
      db_max_index( INVERSE_TARGET, mtar, VERSION_NORMAL, GET );
      for ( itar=0; itar<=mtar; itar++ ) {
        if ( db_active_index( INVERSE_TARGET, itar, VERSION_NORMAL ) ) {
          if ( db_active_index( INVERSE_TARGET_TIMESTEP, itar, VERSION_NORMAL ) )
            db( INVERSE_TARGET_TIMESTEP, itar, &inverse_target_timestep, 
              ddum, ldum, VERSION_NORMAL, GET );
          else
            inverse_target_timestep = max_time;
          if ( icontrol==inverse_target_timestep ) {
            inverse_target_data = db_dbl( INVERSE_TARGET_DATA, itar, VERSION_NORMAL );
            inverse_target = db_int( INVERSE_TARGET, itar, VERSION_NORMAL );
            data_item_name = inverse_target[0];
            data_item_index = inverse_target[1];
            data_item_number = inverse_target[2];
            actual_values = db_dbl( data_item_name, data_item_index, VERSION_NORMAL );
            if ( data_item_number<0 ) {
              array_member( dof_label, data_item_number, nuknwn, data_item_number );
              if ( db_len(data_item_name,data_item_index,VERSION_NORMAL)==npuknwn ) 
                data_item_number /= nder;
            }
            if ( data_item_number<0 || data_item_number>=
              db_len(data_item_name,data_item_index, VERSION_NORMAL) ) 
              db_error( INVERSE_TARGET, data_item_index );
            length_inverse_target_data = db_len( INVERSE_TARGET_DATA, 
              itar, VERSION_NORMAL );
            if ( length_inverse_target_data==2 ) {
              value = inverse_target_data[0];
              penalty = inverse_target_data[1];
              difference = actual_values[data_item_number] - value;
            }
            else {
              assert( length_inverse_target_data==3 );
              lower_value = inverse_target_data[0];
              higher_value = inverse_target_data[1];
              penalty = inverse_target_data[2];
              if ( actual_values[data_item_number]<lower_value )
                difference = actual_values[data_item_number] - lower_value;
              else if ( actual_values[data_item_number]>higher_value )
                difference = actual_values[data_item_number] - higher_value;
              else
                difference = 0.;
            }
            inverse_parameter_sensitivity[ipar_i] += penalty * 
              difference * difference;
          }
        }
      }
      db( INVERSE_PARAMETER_SENSITIVITY, ipar, idum, 
        inverse_parameter_sensitivity, ipar_n, VERSION_NORMAL, PUT );
    }
    if ( task==INVERSE_DETERMINE_NEW_ESTIMATES ) {
      if ( swit ) {
        pri( "In routine INVERSE_CALCULATION" );
        pri( "determine new estimates" );
      }
      for ( ipar=0; ipar<npar; ipar++ ) {
        if ( db_active_index( INVERSE_PARAMETER, ipar, VERSION_NORMAL ) ) {
          inverse_parameter = db_int( INVERSE_PARAMETER, ipar, VERSION_NORMAL );
          data_item_name = inverse_parameter[0];
          data_item_index = inverse_parameter[1];
          data_item_number = inverse_parameter[2];
          if ( data_item_number<0 ) {
            array_member(dof_label,data_item_number,nuknwn,data_item_number);
            if ( db_len(data_item_name,data_item_index,VERSION_NORMAL)==npuknwn ) 
              data_item_number /= nder;
          }
          parameter = db_dbl( data_item_name, data_item_index, VERSION_NORMAL );
          db( INVERSE_PARAMETER_SENSITIVITY, ipar, idum, 
            inverse_parameter_sensitivity, ldum, VERSION_NORMAL, GET );
          if ( db_active_index( INVERSE_PARAMETER_VARIATION, ipar, VERSION_NORMAL ) )
            db( INVERSE_PARAMETER_VARIATION, ipar, idum, 
              &inverse_parameter_variation, ldum, VERSION_NORMAL, GET );
          else
            inverse_parameter_variation = DEFAULT_VARIATION;
          if ( db_active_index( INVERSE_PARAMETER_STEP, ipar, VERSION_NORMAL ) )
            db( INVERSE_PARAMETER_STEP, ipar, idum,
              &inverse_parameter_step, ldum, VERSION_NORMAL, GET );
          else
            inverse_parameter_step = 1.e10;
          if ( db_active_index( INVERSE_HISTORY, 0, VERSION_NORMAL ) ) {
            db( INVERSE_HISTORY, 0, idum, inverse_history, ldum, VERSION_NORMAL, GET );
            inverse_history_step = inverse_history[0];
          }
          if ( scalar_dabs( parameter[data_item_number] ) > EPS_PARAMETER ) {
            delta = parameter[data_item_number] * inverse_parameter_variation *
              inverse_history_step; 
            if ( delta!=0. ) {
              dfdx = 
                ( inverse_parameter_sensitivity[1] - 
                   inverse_parameter_sensitivity[0] ) / ( 2.*delta );
              d2fdx2 = 
                ( ( inverse_parameter_sensitivity[1] -
                    inverse_parameter_sensitivity[2] ) -
                  ( inverse_parameter_sensitivity[2] -
                    inverse_parameter_sensitivity[0] ) ) / ( delta*delta );
              if ( d2fdx2>0. )
                    /* quadratic approximation to squared sum */
                tmp = -dfdx/d2fdx2;
              else if ( dfdx!=0. )
                  /* linear approximation to squared sum,
                     because second derivative of the squared sum is 
                     not positive (thus no minimum) */
                tmp = -inverse_parameter_sensitivity[2]/dfdx;
              else
                tmp = 0.;
              if      ( tmp>(+inverse_parameter_step) ) 
                tmp = +inverse_parameter_step;
              else if ( tmp<(-inverse_parameter_step) ) 
                tmp = -inverse_parameter_step;
              parameter[data_item_number] += inverse_history_step * tmp;
              if ( db_active_index( INVERSE_PARAMETER_LIMITS, ipar, VERSION_NORMAL ) ) {
                db( INVERSE_PARAMETER_LIMITS, ipar, idum, 
                  inverse_parameter_limits, ldum, VERSION_NORMAL, GET );
                if ( parameter[data_item_number]<inverse_parameter_limits[0] )
                  parameter[data_item_number] = inverse_parameter_limits[0];
                else if ( parameter[data_item_number]>inverse_parameter_limits[1])
                  parameter[data_item_number] = inverse_parameter_limits[1];
              }
              if ( swit ) {
                pri( "ipar", ipar );
                pri( "inverse_parameter_sensitivity", 
                  inverse_parameter_sensitivity, ipar_n );
                pri( "dfdx", dfdx );
                pri( "d2fdx2", d2fdx2 );
                pri( "inverse_parameter_variation", inverse_parameter_variation );
                pri( "parameter", parameter[data_item_number] );
              }
            }
          }
        }
      }

        // store history for next inverse iteration 
      if ( db_active_index( INVERSE_HISTORY, 0, VERSION_NORMAL ) ) {
        db( INVERSE_HISTORY, 0, idum, inverse_history, ldum, VERSION_NORMAL, GET );
        inverse_history_step = inverse_history[0];
        inverse_history_sum = inverse_history[1];
        if ( inverse_parameter_sensitivity[2] > inverse_history_sum )
          inverse_history_step /= 2.;
      }
      inverse_history_sum = inverse_parameter_sensitivity[2];
      inverse_history[0] = inverse_history_step;
      inverse_history[1] = inverse_history_sum;
      length = 2; db( INVERSE_HISTORY, 0, idum, inverse_history, 
        length, VERSION_NORMAL, PUT );
      if ( swit ) pri( "inverse_history", inverse_history, 2 );

      db_delete( INVERSE_PARAMETER_SENSITIVITY, VERSION_NORMAL );

      if ( swit ) cout << "Out routine INVERSE_CALCULATION\n";
    }
  }

}
