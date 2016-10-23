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
/*  Modified on April 1st 2011 by Fernando Lorenzo to get the Von Mises stresses 
	to print in the results file and for post processing
*/

#include "tochnog.h"

#define NTYPE 3
#define EPS_PRI 1.e-12

void calculate( void )

{
  long int idim=0, jdim=0, indx=0, length=0, length_result=0,
    icalcul=0, nmat=0, unknown=0, max_type_post=0, ipost=0, itype=0, 
    ncalcul=0, max_node=0, ldum=0, idum[1], dof_amount[MUKNWN],
    calcul[DATA_ITEM_SIZE], type_post_dof[NTYPE], 
    type_post_dof_calcul[NTYPE], post_calcul_scal_vec_mat[DATA_ITEM_SIZE],
    post_calcul_unknown_operat[DATA_ITEM_SIZE];
  double ddum[1], unknown_values[MUKNWN], dof_calcul[DATA_ITEM_SIZE], 
    result[DATA_ITEM_SIZE], *coord=NULL, *post_dof=NULL;
  char str[MCHAR], outname[MCHAR], outname_without_extension[MCHAR], unknown_name[MCHAR];

  if ( db_active_index( POST_CALCUL, 0, VERSION_NORMAL ) ) {

    db( DOF_AMOUNT, 0, dof_amount, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );

    db_delete( NODE_DOF_CALCUL, VERSION_NORMAL );
    db_delete( POST_POINT_DOF_CALCUL, VERSION_NORMAL );
    db_delete( POST_LINE_DOF_CALCUL, VERSION_NORMAL );
    db_delete( POST_QUADRILATERAL_DOF_CALCUL, VERSION_NORMAL );

    type_post_dof[0] = POST_LINE_DOF;
    type_post_dof_calcul[0] = POST_LINE_DOF_CALCUL;
    type_post_dof[1] = POST_POINT_DOF;
    type_post_dof_calcul[1] = POST_POINT_DOF_CALCUL;
    type_post_dof[2] = POST_QUADRILATERAL_DOF;
    type_post_dof_calcul[2] = POST_QUADRILATERAL_DOF_CALCUL;

    db( POST_CALCUL, 0, calcul, ddum, length, VERSION_NORMAL, GET );
    if ( (length%2)!=0 ) db_error( POST_CALCUL, 0 );
    nmat = length / 2; 
    for ( icalcul=0; icalcul<nmat; icalcul++ ) {
      unknown = calcul[icalcul*2+0];
      calcul_matrix = calcul_vector = calcul_ecomplex = 0;
      if       ( unknown==-CONDIF_TEMPERATURE ) {
        calcul_scalar_indx = temp_indx;
        strcpy( unknown_name, "temp" );
      }      
      else if ( unknown==-GROUNDFLOW_PRESSURE ) {
        calcul_scalar_indx = pres_indx;
        strcpy( unknown_name, "pres" );
      }      
      else if ( unknown==-MATERI_STRESS ) {
        calcul_matrix = 1;
        calcul_mat_indx = stres_indx;
        strcpy( unknown_name, "sig" );
      }
      else if ( unknown==-MATERI_STRAIN_ELASTI ) {
        calcul_matrix = 1;
        calcul_mat_indx = epe_indx;
        strcpy( unknown_name, "epe" );
      }
      else if ( unknown==-MATERI_STRAIN_PLASTI ) {
        calcul_matrix = 1;
        calcul_mat_indx = epp_indx;
        strcpy( unknown_name, "epp" );
      }
      else if ( unknown==-MATERI_STRAIN_TOTAL ) {
        calcul_matrix = 1;
        calcul_mat_indx = ept_indx;
        strcpy( unknown_name, "ept" );
      }
      else if ( unknown==-MATERI_VELOCITY ) {
        calcul_vector = 1;
        calcul_vec_indx = vel_indx;
        strcpy( unknown_name, "vel" );
      }
      else if ( unknown==-MATERI_DISPLACEMENT ) {
        calcul_vector = 1;
        calcul_vec_indx = dis_indx;
        strcpy( unknown_name, "dis" );
      }
      else if ( unknown==-MAXWELL_ECOMPLEX ) {
        calcul_ecomplex = 1;
        strcpy( unknown_name, "ec" );
      }
      else if ( unknown==-MAXWELL_E ) {
        calcul_vector = 1;
        calcul_vec_indx = maxe_indx;
        strcpy( unknown_name, "e" );
      }
      else if ( unknown==-MAXWELL_EI ) {
        calcul_vector = 1;
        calcul_vec_indx = maxei_indx;
        strcpy( unknown_name, "ei" );
      }
      else if ( unknown==-MAXWELL_ER ) {
        calcul_vector = 1;
        calcul_vec_indx = maxer_indx;
        strcpy( unknown_name, "er" );
      }
      else if ( unknown==-MAXWELL_FE ) {
        calcul_vector = 1;
        calcul_vec_indx = maxfe_indx;
        strcpy( unknown_name, "fe" );
      }
      else {
        db_error( POST_CALCUL, 0 );
      }
      calcul_operat = calcul[icalcul*2+1];

      if ( calcul_matrix ) {
        if ( dof_amount[calcul_mat_indx]==MDIM*MDIM ) {
          pri( "Error: POST_CALCUL is not available for non-symmetric matrices" );
          exit(TN_EXIT_STATUS);
        }
      }

        // to prevent memory problems in parallel comp.
      db_max_index( NODE, max_node, VERSION_NORMAL, GET );
      db_allocate( -NODE_DOF_CALCUL, max_node, VERSION_NORMAL, MINIMAL );
      parallel_sys_routine( &parallel_calcul_node );

      for ( itype=0; itype<NTYPE; itype++ ) {
        db_max_index( type_post_dof[itype], max_type_post, VERSION_NORMAL, GET );
        for ( ipost=0; ipost<=max_type_post; ipost++ ) {
          if ( db_active_index( type_post_dof[itype], ipost, VERSION_NORMAL ) ) {
            post_dof = db_dbl( type_post_dof[itype], ipost, VERSION_NORMAL );
            if ( groundflow_pressure ) {
              if ( type_post_dof[itype]==POST_POINT_DOF )
                coord = db_dbl( POST_POINT, ipost, VERSION_NORMAL );
              else {
                pri( "Not available for groundflow analysis: ", -type_post_dof[itype] );
                exit(TN_EXIT_STATUS);
              }
            }
            else
              coord = ddum;
            if      ( calcul_ecomplex ) {
              array_move( post_dof, unknown_values, nuknwn );
            }
            else if ( calcul_matrix ) {
              for ( idim=0; idim<MDIM; idim++ ) {
                for ( jdim=0; jdim<MDIM; jdim++ ) {
                  indx = idim*MDIM + jdim;
                  unknown_values[indx] = 
                    post_dof[calcul_mat_indx+stress_indx(idim,jdim)*nder];
                }
              }
            }
            else if ( calcul_vector ) {
              for ( idim=0; idim<ndim; idim++ ) {
                indx = idim;
                unknown_values[indx] = post_dof[calcul_vec_indx+idim*nder];
              }
            }
            else
              unknown_values[0] = post_dof[calcul_scalar_indx];
            calculate_operat( unknown_values, -1, coord, post_dof, result, length_result );
            db( type_post_dof_calcul[itype], ipost, idum, dof_calcul, 
              length, VERSION_NORMAL, GET_IF_EXISTS );
            if ( length+length_result>DATA_ITEM_SIZE ) 
              db_error( POST_CALCUL, 0 );
            array_move( result, &dof_calcul[length], length_result );
            length += length_result;
            db( type_post_dof_calcul[itype], ipost, idum, dof_calcul, 
              length, VERSION_NORMAL, PUT );
          }
        }
      }

      if      ( labs(calcul_operat)==ABSOL ) {
        strcpy( outname, "abs" );
        strcat( outname, unknown_name );
        ncalcul++;
        strcpy( post_calcul_names[ncalcul-1], outname );
        strcpy( post_calcul_names_without_extension[ncalcul-1], outname );
        post_calcul_scal_vec_mat[ncalcul-1] = -SCALAR;
        post_calcul_unknown_operat[(ncalcul-1)*2+0] = unknown;
        post_calcul_unknown_operat[(ncalcul-1)*2+1] = calcul_operat;
      }            
      else if ( labs(calcul_operat)==AVERAGE ) {
        strcpy( outname, "a" );
        strcat( outname, unknown_name );
        ncalcul++;
        strcpy( post_calcul_names[ncalcul-1], outname );
        strcpy( post_calcul_names_without_extension[ncalcul-1], outname );
        post_calcul_scal_vec_mat[ncalcul-1] = -SCALAR;
        post_calcul_unknown_operat[(ncalcul-1)*2+0] = unknown;
        post_calcul_unknown_operat[(ncalcul-1)*2+1] = calcul_operat;
      }
      else if ( labs(calcul_operat)==NEGATIVE ) {
        strcpy( outname, "n" );
        strcat( outname, unknown_name );
        ncalcul++;
        strcpy( post_calcul_names[ncalcul-1], outname );
        strcpy( post_calcul_names_without_extension[ncalcul-1], outname );
        post_calcul_scal_vec_mat[ncalcul-1] = -SCALAR;
        post_calcul_unknown_operat[(ncalcul-1)*2+0] = unknown;
        post_calcul_unknown_operat[(ncalcul-1)*2+1] = calcul_operat;
      }
      else if ( labs(calcul_operat)==POSITIVE ) {
        strcpy( outname, "p" );
        strcat( outname, unknown_name );
        ncalcul++;
        strcpy( post_calcul_names[ncalcul-1], outname );
        strcpy( post_calcul_names_without_extension[ncalcul-1], outname );
        post_calcul_scal_vec_mat[ncalcul-1] = -SCALAR;
        post_calcul_unknown_operat[(ncalcul-1)*2+0] = unknown;
        post_calcul_unknown_operat[(ncalcul-1)*2+1] = calcul_operat;
      }
      else if ( labs(calcul_operat)==PRIVAL ) {
        strcpy( outname, "va" );
        for ( idim=0; idim<MDIM; idim++ ) {
          strcpy( outname, "va" );
          long_to_a( idim, str );
          strcat( outname, str );
          strcat( outname, unknown_name );
          ncalcul++;
          strcpy( post_calcul_names[ncalcul-1], outname );
          strcpy( post_calcul_names_without_extension[ncalcul-1], outname );
          post_calcul_scal_vec_mat[ncalcul-1] = -SCALAR;
          post_calcul_unknown_operat[(ncalcul-1)*2+0] = unknown;
          post_calcul_unknown_operat[(ncalcul-1)*2+1] = calcul_operat;
        }
      }
      else if ( labs(calcul_operat)==PRIVEC ) {
        for ( idim=0; idim<MDIM; idim++ ) {
          strcpy( str, "ve" );
          long_to_a( idim, outname );
          strcat( str, outname );
          for ( jdim=0; jdim<MDIM; jdim++ ) {
            strcpy( outname, str );
            strcpy( outname_without_extension, str );
            if      ( jdim==0 )
              strcat( outname, "x" );
            else if ( jdim==1 )
              strcat( outname, "y" );
            else {
              assert( jdim==2 );
              strcat( outname, "z" );
            }
            strcat( outname, unknown_name );
            strcat( outname_without_extension, unknown_name );
            ncalcul++;
            strcpy( post_calcul_names[ncalcul-1], outname );
            strcpy( post_calcul_names_without_extension[ncalcul-1], 
              outname_without_extension );
            post_calcul_scal_vec_mat[ncalcul-1] = -VECTOR;
            post_calcul_unknown_operat[(ncalcul-1)*2+0] = unknown;
            post_calcul_unknown_operat[(ncalcul-1)*2+1] = calcul_operat;
          }
        }
      }
      else if ( labs(calcul_operat)==SIZETOT ) {
        strcpy( outname, "s" );
        strcat( outname, unknown_name );
        ncalcul++;
        strcpy( post_calcul_names[ncalcul-1], outname );
        strcpy( post_calcul_names_without_extension[ncalcul-1], outname );
        post_calcul_scal_vec_mat[ncalcul-1] = -SCALAR;
        post_calcul_unknown_operat[(ncalcul-1)*2+0] = unknown;
        post_calcul_unknown_operat[(ncalcul-1)*2+1] = calcul_operat;
      }
      else if ( labs(calcul_operat)==SIZEDEV ) {
        strcpy( outname, "sd" );
        strcat( outname, unknown_name );
        ncalcul++;
        strcpy( post_calcul_names[ncalcul-1], outname );
        strcpy( post_calcul_names_without_extension[ncalcul-1], outname );
        post_calcul_scal_vec_mat[ncalcul-1] = -SCALAR;
        post_calcul_unknown_operat[(ncalcul-1)*2+0] = unknown;
        post_calcul_unknown_operat[(ncalcul-1)*2+1] = calcul_operat;
	  }
      else if ( labs(calcul_operat)==MISES ) {
		  strcpy( outname, "mises-" );
		  strcat( outname, unknown_name );
		  ncalcul++;
		  strcpy( post_calcul_names[ncalcul-1], outname );
		  strcpy( post_calcul_names_without_extension[ncalcul-1], outname );
		  post_calcul_scal_vec_mat[ncalcul-1] = -SCALAR;
		  post_calcul_unknown_operat[(ncalcul-1)*2+0] = unknown;
		  post_calcul_unknown_operat[(ncalcul-1)*2+1] = calcul_operat;
      }
      else if ( unknown==-MATERI_STRESS && labs(calcul_operat)==PHIMOB ) {
        strcpy( outname, "pm" );
        strcat( outname, unknown_name );
        ncalcul++;
        strcpy( post_calcul_names[ncalcul-1], outname );
        strcpy( post_calcul_names_without_extension[ncalcul-1], outname );
        post_calcul_scal_vec_mat[ncalcul-1] = -SCALAR;
        post_calcul_unknown_operat[(ncalcul-1)*2+0] = unknown;
        post_calcul_unknown_operat[(ncalcul-1)*2+1] = calcul_operat;
      }
      else if ( unknown==-MATERI_STRESS && labs(calcul_operat)==TOTAL &&
                groundflow_pressure ) {
        strcpy( str, "to" );
        for ( idim=0; idim<MDIM; idim++ ) {
          for ( jdim=idim; jdim<MDIM; jdim++ ) {
            strcpy( outname, str );
            strcpy( outname_without_extension, str );
            strcat( outname, unknown_name );
            strcat( outname_without_extension, unknown_name );
            if      ( idim==0 && jdim==0 )
              strcat( outname, "xx" );
            else if ( idim==0 && jdim==1 )
              strcat( outname, "xy" );
            else if ( idim==0 && jdim==2 )
              strcat( outname, "xz" );
            else if ( idim==1 && jdim==1 )
              strcat( outname, "yy" );
            else if ( idim==1 && jdim==2 )
              strcat( outname, "yz" );
            else if ( idim==2 && jdim==2 )
              strcat( outname, "zz" );
            ncalcul++;
            strcpy( post_calcul_names[ncalcul-1], outname );
            strcpy( post_calcul_names_without_extension[ncalcul-1], 
              outname_without_extension );
            post_calcul_scal_vec_mat[ncalcul-1] = -MATRIX;
            post_calcul_unknown_operat[(ncalcul-1)*2+0] = unknown;
            post_calcul_unknown_operat[(ncalcul-1)*2+1] = calcul_operat;
          }
        }
      }
      else if ( unknown==-GROUNDFLOW_PRESSURE && labs(calcul_operat)==TOTAL &&
                groundflow_pressure ) {
        strcpy( outname, "to" );
        strcat( outname, unknown_name );
        ncalcul++;
        strcpy( post_calcul_names[ncalcul-1], outname );
        strcpy( post_calcul_names_without_extension[ncalcul-1], outname );
        post_calcul_scal_vec_mat[ncalcul-1] = -SCALAR;
        post_calcul_unknown_operat[(ncalcul-1)*2+0] = unknown;
        post_calcul_unknown_operat[(ncalcul-1)*2+1] = calcul_operat;
      }             
      else if ( unknown==-GROUNDFLOW_PRESSURE && labs(calcul_operat)==STATIC &&
                groundflow_pressure ) {
        strcpy( outname, "st" );
        strcat( outname, unknown_name );
        ncalcul++;
        strcpy( post_calcul_names[ncalcul-1], outname );
        strcpy( post_calcul_names_without_extension[ncalcul-1], outname );
        post_calcul_scal_vec_mat[ncalcul-1] = -SCALAR;
        post_calcul_unknown_operat[(ncalcul-1)*2+0] = unknown;
        post_calcul_unknown_operat[(ncalcul-1)*2+1] = calcul_operat;
      }                
      else if ( unknown==-GROUNDFLOW_PRESSURE && labs(calcul_operat)==DYNAMIC &&
                groundflow_pressure ) {
        strcpy( outname, "dy" );
        strcat( outname, unknown_name );
        ncalcul++;
        strcpy( post_calcul_names[ncalcul-1], outname );
        strcpy( post_calcul_names_without_extension[ncalcul-1], outname );
        post_calcul_scal_vec_mat[ncalcul-1] = -SCALAR;
        post_calcul_unknown_operat[(ncalcul-1)*2+0] = unknown;
        post_calcul_unknown_operat[(ncalcul-1)*2+1] = calcul_operat;           
      }
      else
        db_error( POST_CALCUL, 0 );
      assert( ncalcul<db_data_length(POST_CALCUL_SCAL_VEC_MAT)-1 );

    }
    if ( ncalcul>0 ) {
      length = ncalcul;
      db( POST_CALCUL_SCAL_VEC_MAT, 0, post_calcul_scal_vec_mat, 
        ddum, length, VERSION_NORMAL, PUT );
      length = 2*ncalcul;
      db( POST_CALCUL_UNKNOWN_OPERAT, 0, post_calcul_unknown_operat, 
        ddum, length, VERSION_NORMAL, PUT );
    }

  }

}

void parallel_calcul_node( void )

{
  long int length=0, idim=0, jdim=0, inod=0, max_node=0,
    length_result=0, indx=0, iloop=0, nloop=0, ithread=0, idum[1], 
    *next_of_loop=NULL;
  double unknown_values[MUKNWN], result[DATA_ITEM_SIZE],
    dof_calcul[MCALCUL], *coord=NULL, *node_dof=NULL;

  db_max_index( NODE_DOF, max_node, VERSION_NORMAL, GET );
  if ( max_node>=0 ) {
    next_of_loop = get_new_int(1+max_node);
    parallel_sys_next_of_loop( next_of_loop, max_node, nloop, ithread );
    for ( iloop=0; iloop<nloop; iloop++ ) {
      inod = next_of_loop[iloop];
      if ( inod>max_node )
        break;
      else if ( db_active_index( NODE_DOF, inod, VERSION_NORMAL ) ) {
        coord = db_dbl( NODE, inod, VERSION_NORMAL );
        node_dof = db_dbl( NODE_DOF, inod, VERSION_NORMAL );
        if ( calcul_ecomplex )
          array_move( node_dof, unknown_values, nuknwn );
        else if ( calcul_matrix ) {
          for ( idim=0; idim<MDIM; idim++ ) {
            for ( jdim=0; jdim<MDIM; jdim++ ) {
              indx = idim*MDIM + jdim;
              unknown_values[indx] = 
                node_dof[calcul_mat_indx+stress_indx(idim,jdim)*nder];
            }
          }
        }
        else if ( calcul_vector ) {
          for ( idim=0; idim<ndim; idim++ ) {
            indx = idim;
            unknown_values[indx] = node_dof[calcul_vec_indx+idim*nder];
          }
        }
        else
           unknown_values[0] = node_dof[calcul_scalar_indx];
        calculate_operat( unknown_values, inod, coord, node_dof, result, length_result );
        db( NODE_DOF_CALCUL, inod, idum, dof_calcul, length, 
          VERSION_NORMAL, GET_IF_EXISTS );
        if ( (length+length_result)>MCALCUL ) {
          cout << "\nMCALCUL too small. Increase it in tochnog.h and recompile.\n";
          exit(TN_EXIT_STATUS);
        }
        array_move( result, &dof_calcul[length], length_result );
        length += length_result;
        db( NODE_DOF_CALCUL, inod, idum, dof_calcul, length, 
          VERSION_NORMAL, PUT );
      }
    }
    delete[] next_of_loop;
  }

}

void calculate_operat( double unknown_values[], long int inod, 
  double coord[], double dof[], double result[], long int &length_result )

{

  long int idim=0, jdim=0, indx=0, n=0, idum[1];
  double average=0., phim=0., d__1=0., d__2=0., tmax=0., tmin=0., 
    cohesion=0., tmp1=0., tmp2=0., a=0., b=0., c=0., x1=0., x2=0.,
    pres=0., static_pres=0., total_pres=0., location=0.,
    prival[MDIM], privec[MDIM*MDIM], tmp_mat[MDIM*MDIM], 
    dev_mat[MDIM*MDIM], workval[MDIM], workvec[MDIM*MDIM];

  if ( calcul_matrix )
    average = ( unknown_values[0] + unknown_values[4] + 
      unknown_values[8] ) / 3.;
  else if ( calcul_vector ) {
    for ( idim=0; idim<ndim; idim++ )
      average += unknown_values[idim] / ( (double) ndim );
  }

  if ( calcul_matrix ) {
    array_move( unknown_values, tmp_mat, MDIM*MDIM );
    matrix_jacobi( tmp_mat, MDIM, workval, workvec, idum );
    sort( workval, workvec );
    array_move( workval, prival, MDIM );
    array_move( workvec, privec, MDIM*MDIM );
  }
  
  if      ( calcul_ecomplex ) {
    result[0] = 0.;
    for ( idim=0; idim<MDIM; idim++ ) {
      result[0] += scalar_square(unknown_values[maxer_indx+idim*nder]) + 
                   scalar_square(unknown_values[maxei_indx+idim*nder] );
    }
    result[0] = sqrt( result[0] );
    length_result = 1;
  }
  else if ( labs(calcul_operat)==ABSOL ) {
    result[0] = scalar_dabs( unknown_values[0] );
    length_result = 1;
  }
  else if ( labs(calcul_operat)==AVERAGE ) {
    result[0] = average;
    length_result = 1;
  }
  else if ( labs(calcul_operat)==NEGATIVE ) {
    if ( !calcul_matrix ) db_error( POST_CALCUL, 0 );
    result[0] = 0.;
    for ( idim=0; idim<MDIM; idim++ ) {
      if ( prival[idim]<0. ) {
        n++;
        result[0] += prival[idim];
      }
    }
    if ( n>0 ) result[0] = result[0] / n;
    length_result = 1;
  }
  else if ( labs(calcul_operat)==POSITIVE ) {
    if ( !calcul_matrix ) db_error( POST_CALCUL, 0 );
    result[0] = 0.;
    for ( idim=0; idim<MDIM; idim++ ) {
      if ( prival[idim]>0. ) {
        n++;
        result[0] += prival[idim];
      }
    }
    if ( n>0 ) result[0] = result[0] / n;
    length_result = 1;
  }
  else if ( labs(calcul_operat)==PRIVAL ) {
    if ( !calcul_matrix ) db_error( POST_CALCUL, 0 );
    array_move( prival, result, MDIM );
    length_result = MDIM;
  }
  else if ( labs(calcul_operat)==PRIVEC ) {
    if ( !calcul_matrix ) db_error( POST_CALCUL, 0 );
    for ( idim=0; idim<MDIM; idim++ ) {
      indx = idim * MDIM;
      array_normalize( &privec[indx], MDIM );
      array_multiply( &privec[indx], &result[indx], prival[idim], MDIM );
    }
    length_result = MDIM * MDIM;
  }
  else if ( labs(calcul_operat)==SIZETOT ) {
    if ( calcul_matrix )
      result[0] = sqrt( 0.5 * array_inproduct( unknown_values, 
        unknown_values, MDIM*MDIM ) );
    else
      result[0] = sqrt( 0.5 * array_inproduct( unknown_values, 
        unknown_values, ndim ) );
    length_result = 1;
  }
  else if ( labs(calcul_operat)==SIZEDEV ) {
    if ( !calcul_matrix ) db_error( POST_CALCUL, 0 );
    for ( idim=0; idim<MDIM; idim++ ) {
      for ( jdim=0; jdim<MDIM; jdim++ ) {
        indx = idim*MDIM + jdim;
        dev_mat[indx] = unknown_values[indx];
        if ( idim==jdim ) dev_mat[indx] -= average;
      }
    }
    result[0] = sqrt( .5 * array_inproduct( dev_mat, 
      dev_mat, MDIM*MDIM ) );
    length_result = 1;
  }
  else if ( labs(calcul_operat)==MISES ) {
	  if ( !calcul_matrix ) db_error( POST_CALCUL, 0 );
	  for ( idim=0; idim<MDIM; idim++ ) {
		  for ( jdim=0; jdim<MDIM; jdim++ ) {
			  indx = idim*MDIM + jdim;
			  dev_mat[indx] = unknown_values[indx];
			  if ( idim==jdim ) dev_mat[indx] -= average;
		  }
	  }
	  result[0] = sqrt( 1.5 * array_inproduct( dev_mat, 
											  dev_mat, MDIM*MDIM ) );
	  length_result = 1;
  }
  else if ( labs(calcul_operat)==PHIMOB ) {
    matrix_eigenvalues( unknown_values, prival );
    d__1 = scalar_dabs(prival[0]), 
    d__2 = scalar_dabs(prival[1]), 
    d__1 = scalar_dmax(d__1,d__2), 
    d__2 = scalar_dabs(prival[2]);
    tmax = -scalar_dmax(d__1,d__2);
    d__1 = scalar_dabs(prival[0]), 
    d__2 = scalar_dabs(prival[1]), 
    d__1 = scalar_dmin(d__1,d__2), 
    d__2 = scalar_dabs(prival[2]);
    tmin = -scalar_dmin(d__1,d__2);
    a = 0.5 * scalar_dabs(tmax-tmin);
    b = 0.5 * scalar_dabs(tmax+tmin);
    c = - cohesion;
    tmp1 = -a*a + b*b + c*c;
    tmp2 =  b*b + c*c;
    if ( tmp1>=0 && tmp2!=0. ) {
      x1 = acos( scalar_dabs( (-2.*a*c+2.*b*sqrt(tmp1))/(2.*tmp2) ) );
      x2 = acos( scalar_dabs( (-2.*a*c-2.*b*sqrt(tmp1))/(2.*tmp2) ) );
      if ( x1>x2 ) 
        phim = x1;
      else 
        phim = x2;
      if ( phim>PIRAD/2. ) phim = PIRAD/2.;
    }
    else
      phim = PIRAD/2.;
    result[0] = phim * 360. / ( 2.*PIRAD );
    length_result = 1;
  }
  else if ( labs(calcul_operat)==TOTAL ) {
    groundflow_phreatic_coord( inod, coord, dof, total_pres, 
       static_pres, location );
    pres = total_pres;
    if ( calcul_matrix ) {
      result[0] = unknown_values[0] + pres;
      result[1] = unknown_values[1];
      result[2] = unknown_values[2];
      result[3] = unknown_values[4] + pres;
      result[4] = unknown_values[5];
      result[5] = unknown_values[8] + pres;
      length_result = 6;
    }
    else {
      result[0] = pres;
      length_result = 1;                 
    }
  }                        
  else if ( labs(calcul_operat)==STATIC ) {
    pres = 0.;
    if ( groundflow_phreatic_coord( inod, coord, dof, total_pres, 
        static_pres, location ) )
      pres = static_pres;
    result[0] = pres;
    length_result = 1;
  }
  else if ( labs(calcul_operat)==DYNAMIC ) {
    pres = dof[pres_indx];
    if ( groundflow_phreatic_coord( inod, coord, dof, total_pres, 
        static_pres, location ) ) 
      pres = total_pres - static_pres;
    result[0] = pres;
    length_result = 1;                 
  }
  else
    db_error( POST_CALCUL, 0 );

}

