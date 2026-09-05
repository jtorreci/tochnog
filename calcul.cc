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
    // per-node averaged flag of post_calcul -materi_stress -force
    // (calcul_force.cc, consumed by the -primary print): refreshed
    // every calculate() like NODE_DOF_CALCUL.
    db_delete( POST_CALCUL_MATERI_STRESS_FORCE_AVERAGED_NODE, VERSION_NORMAL );

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

      // post_calcul -materi_stress -force (manual Professional 6.913):
      // the family is NODAL - reject the POST_LINE_DOF/POST_POINT_DOF/
      // POST_QUADRILATERAL_DOF records (there is no element behind them
      // to integrate over) and validate the configuration records once
      // per record (calcul_force.cc).
      if ( unknown==-MATERI_STRESS && labs(calcul_operat)==FORCE ) {
        for ( itype=0; itype<NTYPE; itype++ ) {
          db_max_index( type_post_dof[itype], max_type_post,
            VERSION_NORMAL, GET );
          for ( ipost=0; ipost<=max_type_post; ipost++ ) {
            if ( db_active_index( type_post_dof[itype], ipost,
                 VERSION_NORMAL ) ) {
              pri( "Error: post_calcul -materi_stress -force is a NODAL "
                   "calculation; POST_LINE_DOF/POST_POINT_DOF/"
                   "POST_QUADRILATERAL_DOF records are not supported" );
              exit(TN_EXIT_STATUS);
            }
          }
        }
        post_calcul_materi_stress_force_validate();
      }

      if ( calcul_matrix ) {
        if ( dof_amount[calcul_mat_indx]==MDIM*MDIM ) {
          pri( "Error: POST_CALCUL is not available for non-symmetric matrices" );
          exit(TN_EXIT_STATUS);
        }
      }

        // to prevent memory problems in parallel comp.
      db_max_index( NODE, max_node, VERSION_NORMAL, GET );
      db_allocate( -NODE_DOF_CALCUL, max_node, VERSION_NORMAL, MINIMAL );
      db_allocate( -POST_CALCUL_MATERI_STRESS_FORCE_AVERAGED_NODE,
        max_node, VERSION_NORMAL, MINIMAL );
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
        // manual Professional 6.913 area: the item label is -to_pres
        // (post_calcul_label -to_pres in the .dbs of the Professional)
        strcpy( outname, "to_" );
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
        // item label -st_pres (Professional post_calcul_label naming)
        strcpy( outname, "st_" );
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
        // item label -dy_pres (Professional post_calcul_label naming)
        strcpy( outname, "dy_" );
        strcat( outname, unknown_name );
        ncalcul++;
        strcpy( post_calcul_names[ncalcul-1], outname );
        strcpy( post_calcul_names_without_extension[ncalcul-1], outname );
        post_calcul_scal_vec_mat[ncalcul-1] = -SCALAR;
        post_calcul_unknown_operat[(ncalcul-1)*2+0] = unknown;
        post_calcul_unknown_operat[(ncalcul-1)*2+1] = calcul_operat;           
      }
      else if ( unknown==-MATERI_STRESS && groundflow_pressure &&
                ( labs(calcul_operat)==SAFETY_PIPING ||
                  labs(calcul_operat)==SAFETY_LIFTING ) ) {
        // post_calcul -materi_stress -safety_piping/-safety_lifting
        // (manual Professional 6.919): hydraulic piping/lifting safety
        // factors. The generated item names follow the
        // post_calcul_safety_method: -vertical (default) one value with
        // the plain name safety_piping/safety_lifting, -prival three
        // values safety_*_prival_0..2 and -global three values
        // safety_*_global_x/y/z (post_calcul_label naming measured on
        // the Professional .dbs of ground15/16).
        long int safety_method=-VERTICAL, sm_idum[1];
        double sm_ddum[1];
        long int sm_ldum=0;
        db( POST_CALCUL_SAFETY_METHOD, 0, &safety_method, sm_ddum,
          sm_ldum, VERSION_NORMAL, GET_IF_EXISTS );
        char safety_stem[MCHAR];
        strcpy( safety_stem,
          ( labs(calcul_operat)==SAFETY_PIPING ) ?
          "safety_piping" : "safety_lifting" );
        if ( safety_method==-PRIVAL ) {
          for ( idim=0; idim<MDIM; idim++ ) {
            strcpy( outname, safety_stem );
            strcat( outname, "_prival_" );
            long_to_a( idim, str );
            strcat( outname, str );
            ncalcul++;
            strcpy( post_calcul_names[ncalcul-1], outname );
            strcpy( post_calcul_names_without_extension[ncalcul-1], outname );
            post_calcul_scal_vec_mat[ncalcul-1] = -SCALAR;
            post_calcul_unknown_operat[(ncalcul-1)*2+0] = unknown;
            post_calcul_unknown_operat[(ncalcul-1)*2+1] = calcul_operat;
          }
        }
        else if ( safety_method==-GLOBAL ) {
          for ( idim=0; idim<MDIM; idim++ ) {
            strcpy( outname, safety_stem );
            strcat( outname, "_global_" );
            if      ( idim==0 ) strcat( outname, "x" );
            else if ( idim==1 ) strcat( outname, "y" );
            else                strcat( outname, "z" );
            ncalcul++;
            strcpy( post_calcul_names[ncalcul-1], outname );
            strcpy( post_calcul_names_without_extension[ncalcul-1], outname );
            post_calcul_scal_vec_mat[ncalcul-1] = -SCALAR;
            post_calcul_unknown_operat[(ncalcul-1)*2+0] = unknown;
            post_calcul_unknown_operat[(ncalcul-1)*2+1] = calcul_operat;
          }
        }
        else {
          assert( safety_method==-VERTICAL );
          ncalcul++;
          strcpy( post_calcul_names[ncalcul-1], safety_stem );
          strcpy( post_calcul_names_without_extension[ncalcul-1],
            safety_stem );
          post_calcul_scal_vec_mat[ncalcul-1] = -SCALAR;
          post_calcul_unknown_operat[(ncalcul-1)*2+0] = unknown;
          post_calcul_unknown_operat[(ncalcul-1)*2+1] = calcul_operat;
        }
      }
      else if ( unknown==-MATERI_STRESS &&
                ( labs(calcul_operat)==YOUNG_APPARENT ||
                  labs(calcul_operat)==POISSON_APPARENT ) ) {
        // post_calcul -materi_stress -young_apparent/-poisson_apparent
        // (manual Professional 6.903): apparent Young modulus and
        // Poisson ratio from the INCREMENTAL strains and INCREMENTAL
        // stresses of the last time step. One scalar value per node with
        // the plain operator name (the .dbs post_calcul_label of the
        // Professional is "-young_apparent -poisson_apparent", measured
        // on ground17). The value itself is computed in
        // calculate_operat() from the node dofs and the internal
        // node_dof_previous_step snapshot.
        ncalcul++;
        strcpy( post_calcul_names[ncalcul-1],
          ( labs(calcul_operat)==YOUNG_APPARENT ) ?
          "young_apparent" : "poisson_apparent" );
        strcpy( post_calcul_names_without_extension[ncalcul-1],
          post_calcul_names[ncalcul-1] );
        post_calcul_scal_vec_mat[ncalcul-1] = -SCALAR;
        post_calcul_unknown_operat[(ncalcul-1)*2+0] = unknown;
        post_calcul_unknown_operat[(ncalcul-1)*2+1] = calcul_operat;
      }
      else if ( unknown==-MATERI_STRESS && labs(calcul_operat)==FORCE ) {
        long int iforce=0, icomp=0, nforce_stems=0, nforce_comp=0;
        char force_stem[MCHAR], force_comp[MCHAR];
        // post_calcul -materi_stress -force (manual Professional 6.913):
        // normal force, shear force and moment(s) per node. The items are
        // GLOBAL PLOT vectors (components in the structure thickness
        // direction; only the SIZE is the physical value):
        //   2D (9 items):  norx nory nors | shex shey shes | momx momy moms
        //   3D (16 items): norx nory norz nors | shex shey shez shes |
        //                  mom1x mom1y mom1z mom1s | mom2x mom2y mom2z mom2s
        // The flat NODE_DOF_CALCUL layout is ONE slot per item (the
        // SCALARS section of print_vt.cc reads them one by one); the
        // -VECTOR marks group the item with its ndim component siblings
        // (print_vt.cc advances nval=MDIM for -VECTOR entries). Both
        // layouts fit in MCALCUL=20 (9 and 16); the 3D block leaves only
        // 4 slots for other post_calcul items (documented limitation).
        if ( ndim==2 ) { nforce_stems = 3; nforce_comp = 3; }
        else           { nforce_stems = 4; nforce_comp = 4; }
        for ( iforce=0; iforce<nforce_stems; iforce++ ) {
          if      ( iforce==0 ) strcpy( force_stem, "nor" );
          else if ( iforce==1 ) strcpy( force_stem, "she" );
          else if ( ndim==2 )   strcpy( force_stem, "mom" );
          else if ( iforce==2 ) strcpy( force_stem, "mom1" );
          else                  strcpy( force_stem, "mom2" );
          for ( icomp=0; icomp<nforce_comp; icomp++ ) {
            if      ( icomp==0 ) strcpy( force_comp, "x" );
            else if ( icomp==1 ) strcpy( force_comp, "y" );
            else if ( ndim==2 )  strcpy( force_comp, "s" );
            else if ( icomp==2 ) strcpy( force_comp, "z" );
            else                 strcpy( force_comp, "s" );
            strcpy( outname, force_stem );
            strcat( outname, force_comp );
            strcat( outname, "_" );
            strcat( outname, unknown_name );
            ncalcul++;
            strcpy( post_calcul_names[ncalcul-1], outname );
            strcpy( post_calcul_names_without_extension[ncalcul-1], outname );
            post_calcul_scal_vec_mat[ncalcul-1] =
              ( icomp==nforce_comp-1 ) ? -SCALAR : -VECTOR;
            post_calcul_unknown_operat[(ncalcul-1)*2+0] = unknown;
            post_calcul_unknown_operat[(ncalcul-1)*2+1] = calcul_operat;
          }
        }
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

  long int idim=0, jdim=0, indx=0, n=0, ldum=0, idum[1];
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
       static_pres, location, NULL );
    pres = total_pres;
    // node_total_pressure: user override of the calculated total pressure
    // (manual Professional 6.898)
    db( NODE_TOTAL_PRESSURE, inod, idum, &pres, ldum, VERSION_NORMAL,
      GET_IF_EXISTS );
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
        static_pres, location, NULL ) )
      pres = static_pres;
    // node_static_pressure: user override (manual Professional 6.894)
    db( NODE_STATIC_PRESSURE, inod, idum, &pres, ldum, VERSION_NORMAL,
      GET_IF_EXISTS );
    result[0] = pres;
    length_result = 1;
  }
  else if ( labs(calcul_operat)==DYNAMIC ) {
    pres = dof[pres_indx];
    if ( groundflow_phreatic_coord( inod, coord, dof, total_pres,
        static_pres, location, NULL ) )
      pres = total_pres - static_pres;
    // node_dynamic_pressure: user override (manual Professional 6.886)
    db( NODE_DYNAMIC_PRESSURE, inod, idum, &pres, ldum, VERSION_NORMAL,
      GET_IF_EXISTS );
    result[0] = pres;
    length_result = 1;
  }
  else if ( groundflow_pressure && calcul_matrix &&
            ( labs(calcul_operat)==SAFETY_PIPING ||
              labs(calcul_operat)==SAFETY_LIFTING ) ) {
    // post_calcul -materi_stress -safety_piping/-safety_lifting
    // (manual Professional 6.919): hydraulic safety factors
    //   safety_piping  = (sigma_i + p_dynamic)/p_dynamic
    //   safety_lifting = (sigma_i + p_total)/p_total
    // evaluated for the stress of the method of
    // post_calcul_safety_method (-vertical: the vertical normal stress,
    // one value; -prival: the three principal stresses; -global: the
    // three global normal stresses). p_dynamic = p_total - p_static with
    // p_total/p_static from groundflow_phreatic_coord(). Measured
    // against the Professional on ground15/16: the principal values are
    // listed from the MOST compressive (prival_0 = the minimum), while
    // the GNU prival sort is descending, hence the reversed indexing;
    // when the pressure denominator is zero the factor is set to 0
    // (post_calcul_safety_default eps very small, value 0).
    long int safety_method=-VERTICAL, sm_idum[1], ipipe=0, nres=1;
    double sm_ddum[1], safety_pipe=0., p_div=0., safety_max=0.;
    long int sm_ldum=0, safety_maximum_set=0, res_idum[1];
    double res_ddum[1];
    db( POST_CALCUL_SAFETY_METHOD, 0, &safety_method, sm_ddum,
      sm_ldum, VERSION_NORMAL, GET_IF_EXISTS );
    safety_maximum_set = db( POST_CALCUL_SAFETY_MAXIMUM, 0, res_idum,
      &safety_max, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    groundflow_phreatic_coord( inod, coord, dof, total_pres,
      static_pres, location, NULL );
    if ( labs(calcul_operat)==SAFETY_PIPING )
      p_div = total_pres - static_pres;   // p_dynamic
    else
      p_div = total_pres;                  // p_total
    if ( safety_method==-PRIVAL ) nres = MDIM;
    else if ( safety_method==-GLOBAL ) nres = MDIM;
    else assert( safety_method==-VERTICAL );
    for ( ipipe=0; ipipe<nres; ipipe++ ) {
      if ( safety_method==-VERTICAL ) {
        // vertical normal stress: 1D sigma_xx, 2D sigma_yy, 3D sigma_zz
        indx = (ndim-1)*MDIM + (ndim-1);
        safety_pipe = unknown_values[indx];
      }
      else if ( safety_method==-PRIVAL ) {
        // principal stresses, most compressive first (the GNU prival
        // is sorted descending, reversed here)
        safety_pipe = prival[MDIM-1-ipipe];
      }
      else {
        // global normal stresses sigma_xx/sigma_yy/sigma_zz
        assert( safety_method==-GLOBAL );
        indx = ipipe*MDIM + ipipe;
        safety_pipe = unknown_values[indx];
      }
      if ( scalar_dabs(p_div)<TINY )
        result[ipipe] = 0.;
      else {
        result[ipipe] = ( safety_pipe + p_div ) / p_div;
        if ( safety_maximum_set && result[ipipe]>safety_max )
          result[ipipe] = safety_max;
      }
    }
    length_result = nres;
  }
  else if ( labs(calcul_operat)==FORCE ) {
    // post_calcul -materi_stress -force (manual Professional 6.913):
    // per-node normal/shear force and moment(s); the values come from
    // the stub in calcul_force.cc (lot 1: zero-filled layout, the
    // numerical integration lands in L2/L3).
    post_calcul_materi_stress_force( unknown_values, inod, coord, dof,
      result, length_result );
  }
  else if ( labs(calcul_operat)==YOUNG_APPARENT ||
            labs(calcul_operat)==POISSON_APPARENT ) {
    // post_calcul -materi_stress -young_apparent/-poisson_apparent
    // (manual Professional 6.903): apparent Young modulus and Poisson
    // ratio determined from the INCREMENTAL strains and INCREMENTAL
    // stresses of the last time step (0 when the determination is not
    // possible, e.g. almost zero incremental strains). The increment is
    // the difference between the converged node dofs and the internal
    // node_dof_previous_step snapshot (captured by top() at the start
    // of the step). K = dp/deps_v and G = dq/(3*deps_q) from the mean
    // and deviatoric parts; E = 9KG/(3K+G), nu = (3K-2G)/(2(3K+G)).
    // Verified against the Professional .dbs of ground17 (E = 1e7,
    // nu = 0 EXACT from the incremental uniaxial compression).
    if ( calcul_matrix ) {
      double d_sig[MDIM*MDIM], d_ept[MDIM*MDIM];
      double sig_prev[MDIM*MDIM], ept_prev[MDIM*MDIM];
      double p_now=0., p_prev=0., d_p=0., d_vol=0.;
      double dev_s[MDIM*MDIM], dev_e[MDIM*MDIM];
      double q=0., e_q=0., K_app=0., G_app=0.;
      long int have_prev=0, prev_idum[1], prev_len=0;
      double prev_dof[MUKNWN];
      long int ept_present = ( ept_indx>=0 );
      array_set( d_sig, 0., MDIM*MDIM );
      array_set( d_ept, 0., MDIM*MDIM );
      if ( db( NODE_DOF_PREVIOUS_STEP, inod, prev_idum, prev_dof, prev_len,
          VERSION_NORMAL, GET_IF_EXISTS ) && prev_len>0 )
        have_prev = 1;
      if ( have_prev && ept_present ) {
        for ( idim=0; idim<MDIM; idim++ ) {
          for ( jdim=0; jdim<MDIM; jdim++ ) {
            indx = idim*MDIM + jdim;
            long int sidx = stress_indx(idim,jdim);
            d_sig[indx] = unknown_values[indx] -
              prev_dof[stres_indx + sidx*nder];
            d_ept[indx] = dof[ept_indx + sidx*nder] -
              prev_dof[ept_indx + sidx*nder];
          }
        }
        d_p = ( d_sig[0] + d_sig[4] + d_sig[8] ) / 3.;
        d_vol = d_ept[0] + d_ept[4] + d_ept[8];
        for ( idim=0; idim<MDIM; idim++ ) {
          for ( jdim=0; jdim<MDIM; jdim++ ) {
            indx = idim*MDIM + jdim;
            dev_s[indx] = d_sig[indx];
            dev_e[indx] = d_ept[indx];
            if ( idim==jdim ) {
              dev_s[indx] -= d_p;
              dev_e[indx] -= d_vol/3.;
            }
          }
        }
        // incremental bulk and shear modulus (q = sqrt(1.5 s:s),
        // eps_q = sqrt(2/3 e:e), G = dq/(3*deps_q))
        q = sqrt( 1.5 * array_inproduct( dev_s, dev_s, MDIM*MDIM ) );
        e_q = sqrt( (2./3.) * array_inproduct( dev_e, dev_e,
          MDIM*MDIM ) );
        K_app = ( scalar_dabs(d_vol)>TINY ) ? d_p/d_vol : 0.;
        G_app = ( q>1.e-20 && e_q>1.e-20 ) ? q/(3.*e_q) : 0.;
        if ( K_app==0. || G_app==0. ) {
          result[0] = 0.;
        }
        else if ( labs(calcul_operat)==YOUNG_APPARENT ) {
          // E = 9KG/(3K+G)
          result[0] = 9.*K_app*G_app/( 3.*K_app + G_app );
        }
        else {
          // nu = (3K-2G)/(2(3K+G))
          result[0] = ( 3.*K_app - 2.*G_app ) /
            ( 2.*( 3.*K_app + G_app ) );
        }
      }
      else {
        result[0] = 0.;
      }
      length_result = 1;
    }
    else
      db_error( POST_CALCUL, 0 );
  }
  else
    db_error( POST_CALCUL, 0 );

}

