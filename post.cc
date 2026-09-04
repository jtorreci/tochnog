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
/*  +++++++With contribution in post_global routine contributed by
    +++++++Roman Putanowicz dated September 5, 2013  +++++++++++*/

#include "tochnog.h"

#define EPS_RHSIDE_FIXED 1.e-6
#define EPS_LARGE 1.e20
#define MGROUP 1000

void post( long int task )

{
  long int i=0, j=0, n=0, ipost=0, max_post=0,
    post_line_moment=0, post_line_operat=0, ldum=0, idum[1],
    post_point_move=0;
  double h0=0., h1=0., h2=0., h3=0., total_weight=0.,
    segment_size=0., line_size=0., tmp=0., dtime=0.,
    ddum[1], *xi=NULL, *eta=NULL, *weight_xi=NULL, *weight_eta=NULL, 
    vec0[MDIM], vec1[MDIM], vec2[MDIM], vec3[MDIM], 
    post_line[2*MDIM], post_quadrilateral[4*MDIM],
    post_line_dof[MUKNWN], post_quadrilateral_dof[MUKNWN], 
    line_middle[MDIM];

  if ( nuknwn>0 ) {

    // -post_force_edge_summed (manual Professional 6.936): total force
    // following from the -force_edge records, integrated over the edges.
    // Computed on demand (target_item / control_print) at every
    // step_close, BEFORE the control_print section prints the record.
    post_force_edge_summed_calculate();

    db( POST_POINT_MOVE, 0, &post_point_move, ddum, ldum, 
      VERSION_NORMAL, GET_IF_EXISTS );

    db_max_index( POST_POINT, max_post, VERSION_NORMAL, GET );
    for( ipost=0; ipost<=max_post; ipost++ ) {
      if ( db_active_index( POST_POINT, ipost, VERSION_NORMAL ) ) {
        db( POST_POINT, ipost, idum, post_point, ldum, VERSION_NORMAL, GET );
        array_set( post_point_dof, 0., nuknwn );
        post_found = 0;
        parallel_sys_routine( &parallel_post_point );
        if ( post_found ) {
          db( POST_POINT_DOF, ipost, idum, post_point_dof, 
            nuknwn, VERSION_NORMAL, PUT );
          // let the post point follow the material particle,
          // using the interpolated velocity field
          if ( post_point_move==-YES && materi_velocity ) {
            if ( db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET_IF_EXISTS ) ||
                 db( DTIME, 0, idum, &dtime, ldum, VERSION_NORMAL, GET_IF_EXISTS ) ) {
              for ( i=0; i<ndim; i++ )
                post_point[i] += post_point_dof[vel_indx+i*nder]*dtime;
              db( POST_POINT, ipost, idum, post_point, ndim, VERSION_NORMAL, PUT );
            }
          }
        }
      }
    }

    db_max_index( POST_LINE, max_post, VERSION_NORMAL, GET );
    for( ipost=0; ipost<=max_post; ipost++ ) {
      if ( db_active_index( POST_LINE, ipost, VERSION_NORMAL ) ) {
        db( POST_LINE, ipost, idum, post_line, ldum, VERSION_NORMAL, GET );
        if ( db_active_index( POST_LINE_MOMENT, ipost, VERSION_NORMAL ) )
          db( POST_LINE_MOMENT, ipost, &post_line_moment, ddum, ldum, VERSION_NORMAL, GET );
        else
          post_line_moment = -NO;
        array_add( &post_line[0], &post_line[ndim], line_middle, ndim );
        array_multiply( line_middle, line_middle, 0.5, ndim );
        array_subtract( &post_line[0], &post_line[ndim], vec0, ndim );
        line_size = array_size( vec0, ndim );
        if ( db_active_index( POST_LINE_N, ipost, VERSION_NORMAL ) )
          db( POST_LINE_N, ipost, &n, ddum, ldum, VERSION_NORMAL, GET );
        else
          n = 5;
        if ( db_active_index( POST_LINE_OPERAT, ipost, VERSION_NORMAL ) )
          db( POST_LINE_OPERAT, ipost, &post_line_operat, ddum, 
            ldum, VERSION_NORMAL, GET );
        else
          post_line_operat = -AVERAGE;
        if      ( post_line_operat==-AVERAGE ) 
          segment_size = 1./n;
        else if ( post_line_operat==-SUM && n>1 )
          segment_size = line_size/(n-1);
        else
          db_error( POST_LINE_OPERAT, ipost );
        if ( n>0 ) {
          xi = get_new_dbl(n);
          weight_xi = get_new_dbl(n);
          if ( !integration_gauss(n,xi,weight_xi) ) {
            for ( i=0; i<n; i++ ) {
              xi[i] = -1. + i*2./(n-1);
              weight_xi[i] = 1./n;
            }
          }
          array_set( post_line_dof, 0., nuknwn );
          total_weight = 0.;
          for ( i=0; i<n; i++ ) {
            h0 = (1.-xi[i])/2.; h1 = (1.+xi[i])/2.;
            array_set( post_point, 0., ndim );
            array_multiply( &post_line[0*ndim], vec0, h0, ndim );
            array_multiply( &post_line[1*ndim], vec1, h1, ndim );
            array_add( vec0, vec1, post_point, ndim );
            post_found = 0;
            parallel_sys_routine( &parallel_post_point );
            if ( post_found ) {
              if ( post_line_moment==-YES )
                array_multiply( post_point_dof, post_point_dof,
                  (weight_xi[i]*line_size)*(xi[i]*line_size/2.), 
                  nuknwn );
              else {
                if      ( post_line_operat==-AVERAGE )
                  tmp = segment_size;
                else {
                  assert( post_line_operat==-SUM );
                  if ( i==0 || i==(n-1) )
                    tmp = 0.5 * segment_size;
                  else
                    tmp = segment_size;
                }
                array_multiply( post_point_dof, post_point_dof, tmp, nuknwn );
              }
              array_add( post_point_dof, post_line_dof, post_line_dof, nuknwn);
            }
          }
          db( POST_LINE_DOF, ipost, idum, post_line_dof, 
            nuknwn, VERSION_NORMAL, PUT );
          delete[] xi;
          delete[] weight_xi;
        }
      }
    }

    db_max_index( POST_QUADRILATERAL, max_post, VERSION_NORMAL, GET );
    for( ipost=0; ipost<=max_post; ipost++ ) {
      if ( db_active_index( POST_QUADRILATERAL, ipost, VERSION_NORMAL ) ) {
        db( POST_QUADRILATERAL, ipost, idum, post_quadrilateral, ldum, VERSION_NORMAL, GET );
        if ( db_active_index( POST_QUADRILATERAL_N, ipost, VERSION_NORMAL ) )
          db( POST_QUADRILATERAL_N, ipost, &n, ddum, ldum, VERSION_NORMAL, GET );
        else
          n = 5;
        if ( n>0 ) {
          xi = get_new_dbl(n);
          weight_xi = get_new_dbl(n);
          eta = get_new_dbl(n);
          weight_eta = get_new_dbl(n);
          if ( !integration_gauss(n,xi,weight_xi) ) {
            for ( i=0; i<n; i++ ) {
              xi[i] = -1. + i*2./(n-1);
              weight_xi[i] = 1./n;
            }
          }
          if ( !integration_gauss(n,eta,weight_eta) ) {
            for ( i=0; i<n; i++ ) {
              eta[i] = -1. + i*2./(n-1);
              weight_eta[i] = 1./n;
            }
          }
          array_set( post_quadrilateral_dof, 0., nuknwn ); total_weight = 0.;
          for ( i=0; i<n; i++ ) {
            for ( j=0; j<n; j++ ) {
              h0 = (1.-xi[i])*(1.-eta[j])/4.; 
              h1 = (1.+xi[i])*(1.-eta[j])/4.;
              h2 = (1.-xi[i])*(1.+eta[j])/4.;
              h3 = (1.+xi[i])*(1.+eta[j])/4.;
              array_set( post_point, 0., ndim );
              array_multiply( &post_quadrilateral[0*ndim], vec0, h0, ndim );
              array_multiply( &post_quadrilateral[1*ndim], vec1, h1, ndim );
              array_multiply( &post_quadrilateral[2*ndim], vec2, h2, ndim );
              array_multiply( &post_quadrilateral[3*ndim], vec3, h3, ndim );
              array_add( vec0, post_point, post_point, ndim );
              array_add( vec1, post_point, post_point, ndim );
              array_add( vec2, post_point, post_point, ndim );
              array_add( vec3, post_point, post_point, ndim );
              post_found = 0;
              parallel_sys_routine( &parallel_post_point );
              if ( post_found ) {
                total_weight += weight_xi[i] * weight_eta[j];
                array_multiply( post_point_dof, post_point_dof, 
                  weight_xi[i]*weight_eta[j], nuknwn );
                array_add( post_point_dof, post_quadrilateral_dof,
                  post_quadrilateral_dof, nuknwn );
              }
            }
          }
          delete[] xi;
          delete[] weight_xi;
          delete[] eta;
          delete[] weight_eta;
          if ( total_weight ) array_multiply( post_quadrilateral_dof, 
            post_quadrilateral_dof, 1./total_weight, nuknwn );
          db( POST_QUADRILATERAL_DOF, ipost, idum, post_quadrilateral_dof, 
            nuknwn, VERSION_NORMAL, PUT );
        }
      }
    }

  }

  db_max_index( POST_NODE, max_post, VERSION_NORMAL, GET );
  for( ipost=0; ipost<=max_post; ipost++ ) {
    if ( db_active_index( POST_NODE, ipost, VERSION_NORMAL ) ) {
      db( POST_NODE, ipost, post_node, ddum, ldum, VERSION_NORMAL, GET );
      npost_node = 0;
      array_set( post_node_result, 0., DATA_ITEM_SIZE );
      parallel_sys_routine( &parallel_post_node );
      if ( npost_node>0 ) {
        if      ( post_node[1]==-AVERAGE )
          array_multiply( post_node_result, post_node_result, 1./npost_node, post_node_length );
        else if ( post_node[1]!=-SUM && post_node[1]!=-MOMENT )
          db_error( POST_NODE, ipost );
        db( POST_NODE_RESULT, ipost, idum, post_node_result, post_node_length, VERSION_NORMAL, PUT );
      }
    }
  }

  post_global();

  if ( task==YES ) post_integrate();

}

void parallel_post_point( void )

{
  long int element=0, length=0, max_element=0, name=0, inol=0, nnol=0, 
    inod=0, iloop=0, nloop=0, ithread=0, ldum=0, idum[1], el[1+MNOL], nodes[MNOL],
    *next_of_loop=NULL;
  double ddum[1], coords[MNOL*MDIM], tmp_node_dof[MUKNWN], weight[MNOL];

  db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
  if ( max_element>=0 ) {
    next_of_loop = get_new_int(1+max_element);
    parallel_sys_next_of_loop( next_of_loop, max_element, nloop, ithread );
    for ( iloop=0; iloop<nloop; iloop++ ) {
      element = next_of_loop[iloop];
      if ( element>max_element )
        break;
      else if ( !post_found && db_active_index( ELEMENT, element, 
          VERSION_NORMAL ) ) {
        db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
        name = el[0]; nnol = length - 1; array_move( &el[1], nodes, nnol );
        for ( inol=0; inol<nnol; inol++ ) {
          inod = nodes[inol];
          db( NODE_START_REFINED, inod, idum, &coords[inol*ndim], ldum, 
            VERSION_NORMAL, GET );
        }
        if ( point_el( post_point, coords, weight, name, nnol ) ) {
          parallel_sys_lock();
          post_found = 1;
          array_set( post_point_dof, 0., nuknwn );
          for ( inol=0; inol<nnol; inol++ ) {
            inod = nodes[inol];
            if ( db_active_index( NODE_DOF, inod, VERSION_NEW ) )
              db( NODE_DOF, inod, idum, tmp_node_dof, ldum, VERSION_NEW, GET );
            else
              db( NODE_DOF, inod, idum, tmp_node_dof, ldum, VERSION_NORMAL, GET );
            array_multiply( tmp_node_dof, tmp_node_dof, weight[inol], nuknwn );
            array_add( tmp_node_dof, post_point_dof, post_point_dof, nuknwn );
          }
          parallel_sys_unlock();
        }
      }
    }
    delete[] next_of_loop;
  }

}

void parallel_post_node( void )

{
  long int inod=0, max_node=0, iloop=0, nloop=0, found=0,
    ithread=0, length=0, idum[1], *next_of_loop=NULL;
  double rdum=0., ddum[MDIM], work[DATA_ITEM_SIZE], *coord=NULL;

  db_max_index( NODE, max_node, VERSION_NORMAL, GET );
  next_of_loop = get_new_int(1+max_node);
  if ( max_node>=0 ) {
    for ( ;; ) {
      parallel_sys_next_of_loop( next_of_loop, max_node, nloop, ithread );
      for ( iloop=0; iloop<nloop; iloop++ ) {
        inod = next_of_loop[iloop];
        if ( inod>max_node )
          goto after_loop;
        else if ( db_active_index( post_node[0], inod, VERSION_NORMAL ) ) {
          if      ( post_node[2]>0 )
            found = inod==post_node[2];
          else if ( post_node[2]==-ALL )
            found = 1;
          else
            geometry( inod, ddum, &post_node[2], found, rdum, ddum, rdum,
              ddum, NODE_START_REFINED, PROJECT_EXACT, VERSION_NORMAL );
          if ( found ) {
            db( post_node[0], inod, idum, work, length, VERSION_NORMAL, GET );
            parallel_sys_lock();
            npost_node++;
            if ( post_node[0]==-NODE_RHSIDE && post_node[1]==-MOMENT ) {
              post_node_length = MDIM;
              coord = db_dbl( NODE, inod, VERSION_NORMAL );
              if ( ndim>=3 ) {
                post_node_result[1] -= work[vel_indx+2*nder] * coord[0];
                post_node_result[0] += work[vel_indx+2*nder] * coord[1];
                post_node_result[1] += work[vel_indx+0*nder] * coord[2];
                post_node_result[0] -= work[vel_indx+1*nder] * coord[2];
              }
              if ( ndim>=2 ) {
                post_node_result[2] += work[vel_indx+1*nder] * coord[0];
                post_node_result[2] -= work[vel_indx+0*nder] * coord[1];
              }
            }
            else {
              post_node_length = length;
              array_add( work, post_node_result, post_node_result, post_node_length );
            }
            parallel_sys_unlock();
          }
        }
      }
    }
  }
  after_loop:
  delete[] next_of_loop;
  return;
}

void post_node_rhside_fixed_free( void )

{
  long int i=0, n=0, ipuknwn=0, iuknwn=0, inod=0, max_node=0, length=0, 
    ready=0, length_unknowntypes=0, use_this_unknown=0,
    number=0, ldum=0, idum[1], node_bounded[MPUKNWN], 
    dof_principal[MUKNWN], dof_amount[MUKNWN], dof_type[MUKNWN],
    unknowntypes[DATA_ITEM_SIZE];
  double tmp=0., tmp1=0., tmp2=0., ratio=0., post_node_rhside_ratio=0.,
    post_node_rhside_free[MPUKNWN], post_node_rhside_fixed[MPUKNWN], 
    ddum[1], *node_rhside=NULL;

  if ( npuknwn>0 ) {
    array_set( post_node_rhside_fixed, 0., npuknwn );
    array_set( post_node_rhside_free, 0., npuknwn );
    db_highest_index( NODE_RHSIDE, max_node, VERSION_NORMAL );
    db( DOF_PRINCIPAL, 0, dof_principal, ddum, ldum, VERSION_NORMAL, GET );
    db( DOF_AMOUNT, 0, dof_amount, ddum, ldum, VERSION_NORMAL, GET );
    db( DOF_TYPE, 0, dof_type, ddum, ldum, VERSION_NORMAL, GET );
    db( POST_NODE_RHSIDE_RATIO_UNKNOWNTYPES, 0, unknowntypes, ddum, 
      length_unknowntypes, VERSION_NORMAL, GET_IF_EXISTS );
    for ( inod=0; inod<=max_node; inod++ ) {
	    if ( db_active_index( NODE_RHSIDE, inod, VERSION_NORMAL ) ) {
        node_rhside = db_dbl( NODE_RHSIDE, inod, VERSION_NORMAL );
        array_set( node_bounded, 0, npuknwn );
        db( NODE_BOUNDED, inod, node_bounded, ddum, ldum, 
          VERSION_NORMAL, GET_IF_EXISTS );
        for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
          iuknwn = ipuknwn*nder;
          if ( dof_principal[iuknwn]>=0 ) {
            tmp = scalar_dabs(node_rhside[ipuknwn]);
            if ( node_bounded[ipuknwn] ) {
              if ( tmp>post_node_rhside_fixed[ipuknwn] ) 
                post_node_rhside_fixed[ipuknwn] = tmp;
            }
            else {
              if ( tmp>post_node_rhside_free[ipuknwn] ) 
                post_node_rhside_free[ipuknwn] = tmp;
            }
          }
        }
      }
    }
    db( POST_NODE_RHSIDE_FIXED, 0, idum, post_node_rhside_fixed, 
      npuknwn, VERSION_NORMAL, PUT );
    db( POST_NODE_RHSIDE_FREE, 0, idum, post_node_rhside_free, 
      npuknwn, VERSION_NORMAL, PUT );
    ipuknwn = iuknwn = ready = 0;
    while ( !ready ) {
      iuknwn = ipuknwn * nder;
      n = dof_amount[iuknwn];
      tmp1 = tmp2 = 0.;
      for ( i=0; i<n; i++ ) {
        use_this_unknown = 1;
        if ( length_unknowntypes>0 ) {
          if ( !array_member(unknowntypes,dof_type[iuknwn],
            length_unknowntypes,number) ) use_this_unknown = 0;
        }
        if ( use_this_unknown ) {
          tmp1 += post_node_rhside_free[ipuknwn];
          tmp2 += post_node_rhside_fixed[ipuknwn];
        }
        ipuknwn++;
      }
      if ( tmp2>EPS_RHSIDE_FIXED ) {
        ratio = tmp1/tmp2;
        if ( ratio>post_node_rhside_ratio ) post_node_rhside_ratio = ratio;
      }
      ready = (ipuknwn>=npuknwn);
    }
    length = 1;
    db( POST_NODE_RHSIDE_RATIO, 0, idum, &post_node_rhside_ratio, 
      length, VERSION_NORMAL, PUT );
  }

}

void post_global( void )

{
    long int ipost=0, npost=0, post_type=0, element=0, max_element=0,
    max_node=0, length=0, global_elements=0, global_nodes=0,
    global_unknown_number=0, element_empty=0, inod=0, ipuknwn=0, iuknwn=0,
    ldum=0, idum[1], post_global[DATA_ITEM_SIZE],
    node_bounded[MUKNWN], dof_principal[MUKNWN];
    double element_mass=0., element_strainenergy=0.,
    global_mass=0., global_strainenergy=0.,
    element_volume=0., global_volume=0., ddum[1],
    global_unknown_sum[MUKNWN], global_unknown_max[MUKNWN],
    global_unknown_average[MUKNWN], global_unknown_min[MUKNWN], node_dof[MUKNWN];
    
    if ( db_active_index( POST_GLOBAL, 0, VERSION_NORMAL ) ) {
        db( POST_GLOBAL, 0, post_global, ddum, npost, VERSION_NORMAL, GET );
        db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
        db_max_index( NODE, max_node, VERSION_NORMAL, GET );
        if ( max_element>=0 ) {
            for ( element=0; element<=max_element; element++ ) {
                element_empty = -NO;
                db( ELEMENT_EMPTY, element, &element_empty, ddum,
                   ldum, VERSION_NORMAL, GET_IF_EXISTS );
                if ( element_empty==-NO || element_empty==-FRONT ) {
                    if ( db_active_index( ELEMENT_MASS, element, VERSION_NORMAL ) ) {
                        db( ELEMENT_MASS, element, idum, &element_mass,
                           ldum, VERSION_NORMAL, GET );
                        global_mass += element_mass;
                    }
                    if ( db_active_index( ELEMENT_STRAINENERGY, element, VERSION_NORMAL ) ) {
                        db( ELEMENT_STRAINENERGY, element, idum, &element_strainenergy,
                           ldum, VERSION_NORMAL, GET );
                        global_strainenergy += element_strainenergy;
                    }
                    if ( db_active_index( ELEMENT_VOLUME, element, VERSION_NORMAL ) ) {
                        db( ELEMENT_VOLUME, element, idum, &element_volume,
                           ldum, VERSION_NORMAL, GET );
                        global_volume += element_volume;
                    }
                    if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) ) {
                        global_elements++;
                    }
                }
            }
        }
        if ( max_node>=0 ) {
            db( DOF_PRINCIPAL, 0, dof_principal, ddum, ldum,
               VERSION_NORMAL, GET_IF_EXISTS );
            array_set( global_unknown_sum, 0., nuknwn );
            array_set( global_unknown_max, -EPS_LARGE, nuknwn );
            array_set( global_unknown_min, +EPS_LARGE, nuknwn );
            for ( inod=0; inod<=max_node; inod++ ) {
                if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
                    global_nodes++;
                    array_set( node_bounded, 0, npuknwn );
                    db( NODE_BOUNDED, inod, node_bounded, ddum, ldum,
                       VERSION_NORMAL, GET_IF_EXISTS );
                    db( NODE_DOF, inod, idum, node_dof, ldum,
                       VERSION_NORMAL, GET );
                    for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
                        iuknwn = ipuknwn*nder;
                        if ( dof_principal[iuknwn]>=0 ) {
                            if ( !node_bounded[ipuknwn] ) global_unknown_number++;
                        }
                        global_unknown_sum[iuknwn] += node_dof[iuknwn];
                        if ( node_dof[iuknwn]>global_unknown_max[iuknwn] )
                            global_unknown_max[iuknwn] = node_dof[iuknwn];
                        if ( node_dof[iuknwn]<global_unknown_min[iuknwn] )
                            global_unknown_min[iuknwn] = node_dof[iuknwn];
                        /* calculate extreme values for derivatives */
                        { /* LOCAL_BLOCK */
                            int j;
                            for(j=0; j<nder; j++) {
                                if ( node_dof[iuknwn+j]>global_unknown_max[iuknwn+j] )
                                    global_unknown_max[iuknwn+j] = node_dof[iuknwn+j];
                                if ( node_dof[iuknwn+j]<global_unknown_min[iuknwn+j] )
                                    global_unknown_min[iuknwn+j] = node_dof[iuknwn+j];
                            }
                        }  /* END_LOCAL_BLOCK */
                    }
                }
            }
        }
        for ( ipost=0; ipost<npost; ipost++ ) {
            post_type = post_global[ipost];
            length = nuknwn;
            if ( post_type==-GLOBAL_UNKNOWN_AVERAGE && global_nodes>0 ) {
                array_multiply( global_unknown_sum, global_unknown_average,
                               1./global_nodes, length );
                db( -GLOBAL_UNKNOWN_AVERAGE, 0, idum, global_unknown_average,
                   length, VERSION_NORMAL, PUT );
            }
            if ( post_type==-GLOBAL_UNKNOWN_SUM ) db( -GLOBAL_UNKNOWN_SUM, 0,
                                                     idum, global_unknown_sum, length, VERSION_NORMAL, PUT );
            if ( post_type==-GLOBAL_UNKNOWN_MAX ) db( -GLOBAL_UNKNOWN_MAX, 0,
                                                     idum, global_unknown_max, length, VERSION_NORMAL, PUT );
            if ( post_type==-GLOBAL_UNKNOWN_MIN ) db( -GLOBAL_UNKNOWN_MIN, 0,
                                                     idum, global_unknown_min, length, VERSION_NORMAL, PUT );
            length = 1;
            if ( post_type==-GLOBAL_MASS ) db( -GLOBAL_MASS, 0, idum,
                                              &global_mass, length, VERSION_NORMAL, PUT );
            if ( post_type==-GLOBAL_VOLUME ) db( -GLOBAL_VOLUME, 0, idum,
                                                &global_volume, length, VERSION_NORMAL, PUT );
            if ( post_type==-GLOBAL_ELEMENTS ) db( -GLOBAL_ELEMENTS, 0,
                                                  &global_elements, ddum, length, VERSION_NORMAL, PUT );
            if ( post_type==-GLOBAL_NODES ) db( -GLOBAL_NODES, 0, &global_nodes,
                                               ddum, length, VERSION_NORMAL, PUT );
            if ( post_type==-GLOBAL_UNKNOWN_NUMBER ) db( -GLOBAL_UNKNOWN_NUMBER, 0,
                                                        &global_unknown_number, ddum, length, VERSION_NORMAL, PUT );
            if ( post_type==-GLOBAL_STRAINENERGY ) db( -GLOBAL_STRAINENERGY, 0, idum,
                                                      &global_strainenergy, length, VERSION_NORMAL, PUT );
        }
    }
}

// post_force_edge_summed_calculate - the Professional post_global item
// -post_force_edge_summed (manual 6.936): "total force following from
// -force_edge integrated over edges in x,y,z directions",
// number_of_space_dimensions values, stored in the flat record
// post_force_edge_summed (no index). The Professional computes it when
// post_global is -yes (the default); the GNU computes it ON DEMAND -
// when a target_item or a control_print asks for it (the corpus
// elasti6 checks target_item -post_force_edge_summed 0 1 = 10: the
// total vertical force 1.0*10 on the top edge).
//
// The integration mirrors the load assembly of area()
// (force_element_edge): for every element side whose nodes ALL lie on
// the record's geometry entity (or in its node list) and that passes
// the element/group/side/node restrictions, the total force vector
// accumulates, per side node,
//     w_lobatto * ar * load * factor * node_factor * values[dir]
// with the SAME quadrature that distributes the traction to the nodes
// (2D: side length x Lobatto; 3D: face area via triangle fan x tensor
// Lobatto) - so the summed value is exactly the total force the
// record applies. time/load factor: force_element_edge_time /
// _sine / _time_file, default 1; spatial factor:
// force_element_edge_factor (multi_linear_factor_x included via
// force_factor()).
static long int post_fee_border_tria3[] = { 0, 1,  1, 2,  2, 0 };
static long int post_fee_border_tria6[] = { 0, 1, 2,  2, 4, 5,  5, 3, 0 };
static long int post_fee_border_quad4[] = { 0, 1,  1, 3,  3, 2,  2, 0 };
static long int post_fee_border_quad9[] = { 0, 1, 2,  2, 5, 8,  8, 7, 6,  6, 3, 0 };
static long int post_fee_border_quad16[] =
  { 0, 1, 2, 3,  3, 7, 11, 15,  15, 14, 13, 12,  12, 8, 4, 0 };
static long int post_fee_border_tet4[] =
  { 0, 1, 2,  0, 1, 3,  0, 2, 3,  1, 2, 3 };
static long int post_fee_border_hex8[] =
  { 0, 1, 2, 3,  4, 5, 6, 7,  0, 1, 4, 5,  1, 3, 5, 7,
    2, 3, 6, 7,  0, 2, 4, 6 };
static long int post_fee_border_hex27[] =
  { 0, 1, 2, 3, 4, 5, 6, 7, 8,  18, 19, 20, 21, 22, 23, 24, 25, 26,
    0, 1, 2, 9, 10, 11, 18, 19, 20,  2, 5, 8, 11, 14, 17, 20, 23, 26,
    6, 7, 8, 15, 16, 17, 24, 25, 26,  0, 3, 6, 9, 12, 15, 18, 21, 24 };

void post_force_edge_summed_calculate( void )

{
  long int i=0, j=0, inod=0, inol=0, iside=0, nside=0, nnol_side=0,
    inol_side=0, ind=0, imax=0, element=0, max_element=0, length=0,
    ldum=0, nnod=0, nnol=0, name=0, requested=0, found=0, swit=0,
    nfreq=0, ifreq=0, idum[1], el[1+MNOL], nodes[MNOL],
    *area_int=NULL, *sides=NULL, *target=NULL, *cpr=NULL,
    dof_principal[MUKNWN];
  double ddum[1], values[DATA_ITEM_SIZE], total[MDIM], coord[MDIM],
    load=0., factor=0., node_factor=0., ar=0., time_start=0.,
    time_current=0., dtime=0., time_total=0., frequency=0., amplitude=0.,
    w[MNOL], normal_tmp[MDIM], rdum=0., ddum2[1], geom_work[MDIM], tmp=0.,
    iso_l[MNOL],
    wt_l[MNOL], *force_time_tab=NULL, *sine_tab=NULL, *force_vals=NULL,
    vec[MDIM];
  long int geometry_entity[DATA_ITEM_SIZE];

  swit = set_swit(-1,-1,"post_force_edge_summed_calculate");
  if ( swit ) pri( "In routine POST_FORCE_EDGE_SUMMED_CALCULATE" );

  // on demand: any target_item or control_print asking for the item
  db_max_index( TARGET_ITEM, imax, VERSION_NORMAL, GET );
  for ( ind=0; ind<=imax && !requested; ind++ ) {
    if ( !db_active_index( TARGET_ITEM, ind, VERSION_NORMAL ) ) continue;
    target = db_int( TARGET_ITEM, ind, VERSION_NORMAL );
    if ( labs(target[0])==POST_FORCE_EDGE_SUMMED ) requested = 1;
  }
  if ( !requested ) {
    db_max_index( CONTROL_PRINT, imax, VERSION_NORMAL, GET );
    for ( ind=0; ind<=imax && !requested; ind++ ) {
      if ( !db_active_index( CONTROL_PRINT, ind, VERSION_NORMAL ) ) continue;
      length = db_len( CONTROL_PRINT, ind, VERSION_NORMAL );
      cpr = db_int( CONTROL_PRINT, ind, VERSION_NORMAL );
      for ( i=0; i<length; i++ )
        if ( cpr[i]==-POST_FORCE_EDGE_SUMMED ) requested = 1;
    }
  }
  if ( !requested ) return;

  array_set( total, 0., MDIM );
  db( DOF_PRINCIPAL, 0, dof_principal, ddum, ldum, VERSION_NORMAL,
    GET_IF_EXISTS );
  db( TIME_CURRENT, 0, idum, &time_current, ldum, VERSION_NORMAL,
    GET_IF_EXISTS );
  db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET_IF_EXISTS );
  time_total = time_current + dtime;

  db_highest_index( ELEMENT, max_element, VERSION_NORMAL );

  db_max_index( FORCE_ELEMENT_EDGE, imax, VERSION_NORMAL, GET );
  for ( ind=0; ind<=imax; ind++ ) {
    if ( !db_active_index( FORCE_ELEMENT_EDGE, ind, VERSION_NORMAL ) )
      continue;
    force_vals = db_dbl( FORCE_ELEMENT_EDGE, ind, VERSION_NORMAL );
    length = db_len( FORCE_ELEMENT_EDGE, ind, VERSION_NORMAL );
    if ( length>DATA_ITEM_SIZE ) length = DATA_ITEM_SIZE;
    for ( i=0; i<length; i++ ) values[i] = force_vals[i];

    // the "area" record: node list (area[0]>0) or geometry entity
    area_int = db_int( FORCE_ELEMENT_EDGE_GEOMETRY, ind, VERSION_NORMAL );
    nnod = db_len( FORCE_ELEMENT_EDGE_GEOMETRY, ind, VERSION_NORMAL );
    if ( area_int[0]>0 ) {
      // node list: the load applies to the element sides whose nodes
      // are all members of the list
      geometry_entity[0] = 0;
    }
    else {
      geometry_entity[0] = area_int[0];
      geometry_entity[1] = area_int[1];
      if ( geometry_entity[0]>=0 ||
           db_data_class(geometry_entity[0])!=GEOMETRY ) continue;
    }

    // load factor over time (same resolution order as area(): sine,
    // time table, time file, default 1)
    load = 1.;
    if ( db_active_index( FORCE_ELEMENT_EDGE_SINE, ind, VERSION_NORMAL ) ) {
      sine_tab = db_dbl( FORCE_ELEMENT_EDGE_SINE, ind, VERSION_NORMAL );
      nfreq = ( db_len( FORCE_ELEMENT_EDGE_SINE, ind, VERSION_NORMAL ) - 1 ) / 2;
      time_start = sine_tab[0];
      load = 0.;
      if ( time_total>time_start ) {
        for ( ifreq=0; ifreq<nfreq; ifreq++ ) {
          frequency = sine_tab[1+ifreq*2+0];
          amplitude  = sine_tab[1+ifreq*2+1];
          load += amplitude * sin( 2. * PIRAD * frequency * time_total );
        }
      }
    }
    else if ( db_active_index( FORCE_ELEMENT_EDGE_TIME, ind, VERSION_NORMAL ) ) {
      force_time_tab = db_dbl( FORCE_ELEMENT_EDGE_TIME, ind, VERSION_NORMAL );
      length = db_len( FORCE_ELEMENT_EDGE_TIME, ind, VERSION_NORMAL );
      load = 0.;
      force_time( force_time_tab, "FORCE_ELEMENT_EDGE_TIME", length, load );
    }
    else if ( db_active_index( FORCE_ELEMENT_EDGE_TIME_FILE, ind,
        VERSION_NORMAL ) ) {
      long int force_time_file = 0;
      db( FORCE_ELEMENT_EDGE_TIME_FILE, ind, &force_time_file, ddum,
        ldum, VERSION_NORMAL, GET_IF_EXISTS );
      if ( force_time_file==-YES ) {
        load = 0.;
        force_time_file_apply( ind, FORCE_ELEMENT_EDGE_TIME_FILE, load );
      }
    }

    for ( element=0; element<=max_element; element++ ) {
      if ( !db_active_index( ELEMENT, element, VERSION_NORMAL ) ) continue;
      db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
      name = el[0]; nnol = length - 1;
      for ( inol=0; inol<nnol; inol++ ) nodes[inol] = el[1+inol];
      nside = 0; sides = NULL;
      if      ( name==-TRIA3 ) { nside = 3; nnol_side = 2; sides = post_fee_border_tria3; }
      else if ( name==-TRIA6 ) { nside = 3; nnol_side = 3; sides = post_fee_border_tria6; }
      else if ( name==-QUAD4 ) { nside = 4; nnol_side = 2; sides = post_fee_border_quad4; }
      else if ( name==-QUAD9 ) { nside = 4; nnol_side = 3; sides = post_fee_border_quad9; }
      else if ( name==-QUAD16 ){ nside = 4; nnol_side = 4; sides = post_fee_border_quad16; }
      else if ( name==-TET4  ) { nside = 4; nnol_side = 3; sides = post_fee_border_tet4; }
      else if ( name==-HEX8  ) { nside = 6; nnol_side = 4; sides = post_fee_border_hex8; }
      else if ( name==-HEX27 ) { nside = 6; nnol_side = 9; sides = post_fee_border_hex27; }
      if ( !sides ) continue;

      // element-level restrictions (same companion order as area.cc)
      long int grp = 0;
      db( ELEMENT_GROUP, element, &grp, ddum, ldum, VERSION_NORMAL,
        GET_IF_EXISTS );
      if ( db_active_index( FORCE_ELEMENT_EDGE_ELEMENT, ind, VERSION_NORMAL ) ) {
        long int elt[DATA_ITEM_SIZE], length_elt = 0;
        db( FORCE_ELEMENT_EDGE_ELEMENT, ind, elt, ddum, length_elt,
          VERSION_NORMAL, GET );
        if ( !array_member( elt, element, length_elt, ldum ) ) continue;
      }
      if ( db_active_index( FORCE_ELEMENT_EDGE_ELEMENT_GROUP, ind,
          VERSION_NORMAL ) ) {
        long int gl[DATA_ITEM_SIZE], length_gl = 0;
        db( FORCE_ELEMENT_EDGE_ELEMENT_GROUP, ind, gl, ddum, length_gl,
          VERSION_NORMAL, GET );
        if ( !array_member( gl, grp, length_gl, ldum ) ) continue;
      }

      for ( iside=0; iside<nside; iside++ ) {
        // side in the geometry / node list: EVERY side node must match
        found = 1;
        for ( inol_side=0; inol_side<nnol_side && found; inol_side++ ) {
          inol = sides[iside*nnol_side + inol_side];
          inod = nodes[inol];
          if ( geometry_entity[0]==0 ) {
            if ( !array_member( area_int, inod, nnod, ldum ) ) found = 0;
          }
          else {
            geometry( inod, geom_work, geometry_entity, found, rdum,
              normal_tmp, rdum, geom_work, NODE_START_REFINED, PROJECT_EXACT,
              VERSION_NORMAL );
          }
        }
        if ( !found ) continue;
        if ( db_active_index( FORCE_ELEMENT_EDGE_ELEMENT_SIDE, ind,
            VERSION_NORMAL ) ) {
          long int ss[DATA_ITEM_SIZE], length_ss = 0, ok_side = 0;
          db( FORCE_ELEMENT_EDGE_ELEMENT_SIDE, ind, ss, ddum, length_ss,
            VERSION_NORMAL, GET );
          for ( i=0; i+1<length_ss; i+=2 )
            if ( ss[i]==element && ss[i+1]==iside+1 ) ok_side = 1;
          if ( !ok_side ) continue;
        }

        // edge/face measure of the side (same quadrature as area())
        ar = 0.;
        array_set( w, 0., MNOL );
        if ( ndim==2 ) {
          double *c0 = db_dbl( NODE, nodes[sides[iside*nnol_side+0]],
            VERSION_NORMAL );
          double *c1 = db_dbl( NODE, nodes[sides[iside*nnol_side+nnol_side-1]],
            VERSION_NORMAL );
          for ( i=0; i<ndim; i++ ) vec[i] = c1[i] - c0[i];
          ar = array_size( vec, ndim );
          integration_lobatto( nnol_side, iso_l, w );
        }
        else {
          assert( ndim==3 );
          if ( name==-TET4 ) {
            double *c0 = db_dbl( NODE, nodes[sides[iside*3+0]], VERSION_NORMAL );
            double *c1 = db_dbl( NODE, nodes[sides[iside*3+1]], VERSION_NORMAL );
            double *c2 = db_dbl( NODE, nodes[sides[iside*3+2]], VERSION_NORMAL );
            double a=0., b=0., c=0.;
            for ( i=0; i<3; i++ ) {
              vec[i] = c1[i]-c0[i]; a += vec[i]*vec[i];
              vec[i] = c2[i]-c0[i]; b += vec[i]*vec[i];
              vec[i] = c2[i]-c1[i]; c += vec[i]*vec[i];
            }
            a = sqrt(a); b = sqrt(b); c = sqrt(c);
            ar = sqrt( (a+b+c)*(a+b-c)*(a-b+c)*(-a+b+c)/16 );
            w[0] = w[1] = w[2] = 1./3.;
          }
          else {
            long int ind1 = ( name==-HEX8 ) ? 1 : 2;
            long int ind2 = ( name==-HEX8 ) ? 2 : 6;
            long int nq  = ( name==-HEX8 ) ? 2 : 3;
            double *f0 = db_dbl( NODE, nodes[sides[iside*nnol_side+0]], VERSION_NORMAL );
            double *f1 = db_dbl( NODE, nodes[sides[iside*nnol_side+ind1]], VERSION_NORMAL );
            double *f2 = db_dbl( NODE, nodes[sides[iside*nnol_side+ind2]], VERSION_NORMAL );
            double *f3 = db_dbl( NODE, nodes[sides[iside*nnol_side+nnol_side-1]], VERSION_NORMAL );
            ar = triangle_area( f0, f1, f2 ) + triangle_area( f1, f2, f3 );
            integration_lobatto( nq, iso_l, wt_l );
            for ( j=0; j<nq; j++ )
              for ( i=0; i<nq; i++ ) w[j*nq+i] = wt_l[i]*wt_l[j];
          }
        }

        for ( inol_side=0; inol_side<nnol_side; inol_side++ ) {
          inol = sides[iside*nnol_side + inol_side];
          inod = nodes[inol];
          // node-level restriction + factor + per-node factor
          if ( db_active_index( FORCE_ELEMENT_EDGE_NODE, ind,
              VERSION_NORMAL ) ) {
            long int nds[DATA_ITEM_SIZE], length_nds = 0;
            db( FORCE_ELEMENT_EDGE_NODE, ind, nds, ddum, length_nds,
              VERSION_NORMAL, GET );
            if ( !array_member( nds, inod, length_nds, ldum ) ) continue;
          }
          if ( db_active_index( FORCE_ELEMENT_EDGE_ELEMENT_NODE, ind,
              VERSION_NORMAL ) ) {
            long int en[DATA_ITEM_SIZE], length_en = 0, ok_en = 0;
            db( FORCE_ELEMENT_EDGE_ELEMENT_NODE, ind, en, ddum, length_en,
              VERSION_NORMAL, GET );
            if ( en[0]==element )
              for ( i=1; i<length_en; i++ )
                if ( en[i]==inol ) ok_en = 1;
            if ( !ok_en ) continue;
          }
          db( NODE, inod, idum, coord, ndim, VERSION_NORMAL, GET );
          node_factor = 1.;
          if ( db_active_index( FORCE_ELEMENT_EDGE_NODE_FACTOR, ind,
              VERSION_NORMAL ) ) {
            long int length_nf = 0;
            double nf_vals[DATA_ITEM_SIZE];
            db( FORCE_ELEMENT_EDGE_NODE_FACTOR, ind, idum, nf_vals,
              length_nf, VERSION_NORMAL, GET );
            if ( (long int)nf_vals[0]==element )
              for ( j=0; j+1<length_nf; j++ )
                if ( j==inol ) node_factor = nf_vals[j+1];
          }
          force_factor( FORCE_ELEMENT_EDGE_FACTOR, ind, coord, factor );
          // map the record values onto the space directions via the
          // principal unknowns (the same loop as the area() assembly)
          long int iprinc = 0;
          double nf = w[inol_side] * ar * load * factor * node_factor;
          for ( long int ipuknwn=0; ipuknwn<npuknwn && iprinc<length;
              ipuknwn++ ) {
            long int iuknwn = ipuknwn*nder;
            if ( dof_principal[iuknwn]>=0 ) {
              if ( iprinc<ndim ) total[iprinc] += nf * values[iprinc];
              iprinc++;
            }
          }
        }
      }
    }
  }
  length = ndim;
  db( POST_FORCE_EDGE_SUMMED, 0, idum, total, length, VERSION_NORMAL, PUT );

  if ( swit ) pri( "Out function POST_FORCE_EDGE_SUMMED_CALCULATE" );
}

void post_integrate( void )
 
{
  long int ipost=0, max_post=0, data_item_name=0, data_item_index=0, 
    data_item_number=0, length=0, number=0, ldum=0, idum[1], 
    post_integrate[3], dof_label[MUKNWN];
  double dtime=0., post_integrate_result=0., ddum[1], dval[DATA_ITEM_SIZE];

  db( DTIME, 0, idum, &dtime, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db_max_index( POST_INTEGRATE, max_post, VERSION_NORMAL, GET );
  for ( ipost=0; ipost<=max_post; ipost++ ) {
    if ( db_active_index( POST_INTEGRATE, ipost, VERSION_NORMAL ) ) {
      db( POST_INTEGRATE, ipost, post_integrate, ddum, ldum, VERSION_NORMAL, GET );
      data_item_name = post_integrate[0];
      data_item_index = post_integrate[1];
      data_item_number = post_integrate[2];
      if ( db_active_index( data_item_name, data_item_index, VERSION_NORMAL ) ) {
        if ( data_item_number<0 ) {
          array_member(dof_label,data_item_number,nuknwn,number);
          if ( db_len(data_item_name,data_item_index,VERSION_NORMAL)==npuknwn )
            data_item_number = number / nder;
        }
        post_integrate_result = 0.;
        db( POST_INTEGRATE_RESULT, ipost, idum, &post_integrate_result, ldum, 
          VERSION_NORMAL, GET_IF_EXISTS );
        db( data_item_name, data_item_index, idum, dval, length, 
          VERSION_NORMAL, GET );
        post_integrate_result += dval[data_item_number] * dtime;
        length = 1;
        db( POST_INTEGRATE_RESULT, ipost, idum, &post_integrate_result, length, 
          VERSION_NORMAL, PUT );
      }
    }
  }
}                          
