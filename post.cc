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
    post_line_moment=0, post_line_operat=0, ldum=0, idum[1];
  double h0=0., h1=0., h2=0., h3=0., total_weight=0.,
    segment_size=0., line_size=0., tmp=0.,
    ddum[1], *xi=NULL, *eta=NULL, *weight_xi=NULL, *weight_eta=NULL, 
    vec0[MDIM], vec1[MDIM], vec2[MDIM], vec3[MDIM], 
    post_line[2*MDIM], post_quadrilateral[4*MDIM],
    post_line_dof[MUKNWN], post_quadrilateral_dof[MUKNWN], 
    line_middle[MDIM];

  if ( nuknwn>0 ) {

    db_max_index( POST_POINT, max_post, VERSION_NORMAL, GET );
    for( ipost=0; ipost<=max_post; ipost++ ) {
      if ( db_active_index( POST_POINT, ipost, VERSION_NORMAL ) ) {
        db( POST_POINT, ipost, idum, post_point, ldum, VERSION_NORMAL, GET );
        array_set( post_point_dof, 0., nuknwn );
        post_found = 0;
        parallel_sys_routine( &parallel_post_point );
        if ( post_found ) db( POST_POINT_DOF, ipost, idum, post_point_dof, 
          nuknwn, VERSION_NORMAL, PUT );
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
