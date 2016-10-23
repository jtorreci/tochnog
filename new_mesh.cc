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

#define MNOL_SEGMENT 3
#define MNOL_ELEMENT 8
#define FAC 1.e-1

void new_mesh( void )

{
  long int control_mesh_new_mesh_region=0, 
    element=0, inod=0, max_element=0, length=0, element_group=0, 
    inol=0, nnol=0, swit=0, icontrol=0, ldum=0, idum[1], el[MNOL+1], nodes[MNOL];
  double control_new_mesh=0., delta=0., ddum[1], 
    *tmp_node=NULL, *tmp_node_dof=NULL, *tmp_node_start_refined=NULL, 
    *tmp_node_dof_start_refined=NULL;

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  if ( db_active_index( CONTROL_MESH_NEW_MESH, icontrol, VERSION_NORMAL ) ) {

    db( CONTROL_MESH_NEW_MESH, icontrol, idum, &control_new_mesh, 
      ldum, VERSION_NORMAL, GET );

    swit = set_swit(-1,-1,"new_mesh");
    if ( swit ) pri( "In routine NEW_MESH" );

        /* Initialize */
    delta = control_new_mesh;
    db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
    if ( db_active_index( CONTROL_MESH_NEW_MESH_REGION, 
        icontrol, VERSION_NORMAL ) )
      db( CONTROL_MESH_NEW_MESH_REGION, icontrol, &control_mesh_new_mesh_region, 
        ddum, ldum, VERSION_NORMAL, GET );
    else
      control_mesh_new_mesh_region = -ALL;

      /* Put elements and nodes of the 
         element group on VERSION_NEW_MESH_GENERATED */
    db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
    for ( element=0; element<=max_element; element++ ) {
      if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) ) {
        if ( db_active_index( ELEMENT_GROUP, element, VERSION_NORMAL ) )
          db( ELEMENT_GROUP, element, &element_group, ddum, ldum, 
            VERSION_NORMAL, GET );
        else
          element_group = 0;
        if ( control_mesh_new_mesh_region==-ALL ||
             element_group==control_mesh_new_mesh_region ) {
          db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
          nnol = length - 1; array_move( &el[1], nodes, nnol );
          create_element( element, element, el, length,
            VERSION_NORMAL, VERSION_NEW_MESH_GENERATED );
          for ( inol=0; inol<nnol; inol++ ) {
            inod = nodes[inol];
            tmp_node = db_dbl( NODE, inod, VERSION_NORMAL );
            tmp_node_start_refined = 
              db_dbl( NODE_START_REFINED, inod, VERSION_NORMAL );
            if ( nuknwn>0 ) {
              tmp_node_dof = 
                db_dbl( NODE_DOF, inod, VERSION_NORMAL );
              tmp_node_dof_start_refined = 
                db_dbl( NODE_DOF_START_REFINED, inod, VERSION_NORMAL );
            }
            create_node( &inod, 1, inod, tmp_node, tmp_node_dof,
              tmp_node_start_refined, tmp_node_dof_start_refined,
              VERSION_NORMAL, VERSION_NEW_MESH_GENERATED );
          }
          delete_element( element, VERSION_NORMAL );
        }
      }
    }
    mesh_has_changed( VERSION_NORMAL );
    mesh_has_changed( VERSION_NEW_MESH_GENERATED );

      /* Build a new mesh on VERSION_NEW_MESH_GENERATED */
    new_mesh_version( VERSION_NEW_MESH_GENERATED, delta );

      /* Merge VERSION_NEW_MESH_GENERATED with VERSION_NORMAL */
    mesh_add( VERSION_NEW_MESH_GENERATED, VERSION_NORMAL );
    renumbering( VERSION_NORMAL, NO, 1, 1, idum, idum );

    if ( swit ) pri( "Out routine NEW_MESH" );

  }

}

void new_mesh_version( long int version, double delta )

{
  long int inol=0, jnol=0, knol=0, nnol=0, isegment=0,
    inod=0, max_node=0, element=0, 
    max_element=0, length=0, nallocate=0,
    any_node_inside_old_mesh=0, in=0, jn=0, kn=0, idim=0, 
    found_in=0, found_jn=0, found_kn=0, swit=0,
    number_of_boundary_segments=0, closest_segment=0,
    max_node_boundary=0, icontrol=0, control_new_mesh_element=0, 
    node_boundary_in=0, node_boundary_jn=0, node_boundary_kn=0,
    test1=0, test2=0, test3=0, 
    ldum=0, idum[1], number_of_nodes[MDIM], number_of_elements[MDIM],
    nodes[MNOL], el[MNOL+1], *boundary_segments_node_numbers=NULL, 
    *map_node_inside=NULL;
  double distance=0., closest_distance=0.,
    ddum[MDIM], min_coord[MDIM], max_coord[MDIM], coord[MDIM], 
    projection[MDIM], closest_projection[MDIM],
    work[MDIM], weight[MDIM], coord0[MDIM], coord1[MDIM], coord2[MDIM],
    node_dof[MUKNWN], delta_adjusted[MDIM], *boundary_segments_node=NULL;

  swit = set_swit(-1,-1,"new_mesh_version");
  if ( swit ) pri( "In routine NEW_MESH_VERSION" );

  array_set( coord, 0., MDIM );
  array_set( closest_projection, 0., MDIM );

  if ( db_max_index( NODE_BOUNDARY, max_node_boundary, version, GET ) < 0 ) {
    pri( "Error: NODE_BOUNDARY should be specified if you use control_new_mesh." );
    exit(TN_EXIT_STATUS);
  }

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  if ( db_active_index( CONTROL_MESH_NEW_MESH_ELEMENT, icontrol, VERSION_NORMAL ) )
    db( CONTROL_MESH_NEW_MESH_ELEMENT, icontrol, &control_new_mesh_element, 
      ddum, ldum, VERSION_NORMAL, GET );
  else {
    if      ( ndim==1 ) 
      control_new_mesh_element = -BAR2;
    else if ( ndim==2 ) 
      control_new_mesh_element = -TRIA3;
    else {
      assert( ndim==3 );
      control_new_mesh_element = -TET4;
    }
  }

  if ( swit ) pri( "Determine minimal and maximal coordinates." );
  array_set( min_coord, +1.e10, ndim );
  array_set( max_coord, -1.e10, ndim );
  db_max_index( NODE, max_node, version, GET );
  for ( inod=0; inod<=max_node; inod++ ) {
    if ( db_active_index( NODE, inod, version ) ) {
      db( NODE, inod, idum, coord, ldum, version, GET );
      for ( idim=0; idim<ndim; idim++ ) {
        if ( coord[idim]<min_coord[idim] ) min_coord[idim] = coord[idim];
        if ( coord[idim]>max_coord[idim] ) max_coord[idim] = coord[idim];
      }
    }
  }

    // set number of nodes in each direction
  array_set( number_of_nodes, 1, MDIM );
  if      ( ndim==1 ) {
    number_of_nodes[0] = 1 + (long int) ((max_coord[0]-min_coord[0])/delta);
  }
  else if ( ndim==2 ) {
    number_of_nodes[0] = 1 + (long int) ((max_coord[0]-min_coord[0])/delta);
    number_of_nodes[1] = 1 + (long int) ((max_coord[1]-min_coord[1])/delta);
  }
  else {
    assert( ndim==3 );
    number_of_nodes[0] = 1 + (long int) ((max_coord[0]-min_coord[0])/delta);
    number_of_nodes[1] = 1 + (long int) ((max_coord[1]-min_coord[1])/delta);
    number_of_nodes[2] = 1 + (long int) ((max_coord[2]-min_coord[2])/delta);
  }
  for ( idim=0; idim<ndim; idim++ ) {
    if ( number_of_nodes[idim]%2!=0 ) number_of_nodes[idim]++;
    if ( number_of_nodes[idim]<4 ) number_of_nodes[idim] = 4;
    delta_adjusted[idim] = 
      (max_coord[idim]-min_coord[idim]+2.*FAC*delta)/(number_of_nodes[idim]-1);
  }
  if ( swit ) pri( "number_of_nodes", number_of_nodes, MDIM );

  nallocate = MNOL_ELEMENT*(1+number_of_nodes[0])*
    (1+number_of_nodes[1])*(1+number_of_nodes[2]);
  length = nallocate*MNOL_SEGMENT*ndim*ndim;
  boundary_segments_node = get_new_dbl( length );
  array_set( boundary_segments_node, 0., length );
  length = nallocate*MNOL_SEGMENT*ndim;
  boundary_segments_node_numbers = get_new_int( length );
  array_set( boundary_segments_node_numbers, -1, length );

  if ( swit ) pri( "Determine old boundary segments." );
  db_max_index( ELEMENT, max_element, version, GET );
  for ( element=0; element<=max_element; element++ ) {
    if ( db_active_index( ELEMENT, element, version ) ) {
      db( ELEMENT, element, el, ddum, length, version, GET );
      nnol = length - 1; array_move( &el[1], nodes, nnol );
      if      ( ndim==1 ) {
        for ( inol=0; inol<nnol; inol++ ) {
          in = nodes[inol];
          if ( db_active_index( NODE_BOUNDARY, in, version ) )
            found_in = 1;
          else
            found_in = 0;
          if ( found_in ) {
            db( NODE, in, idum, coord, ldum, version, GET );
            array_move(coord, 
              &boundary_segments_node[number_of_boundary_segments*
              ndim*ndim+0*ndim],ndim);
            boundary_segments_node_numbers[number_of_boundary_segments*
              ndim+0] = in;
            number_of_boundary_segments++;
          }
        }
      }
      else if ( ndim==2 ) {
        for ( inol=0; inol<nnol; inol++ ) {
          in = nodes[inol];
          if ( db_active_index( NODE_BOUNDARY, in, version ) )
            found_in = 1;
          else
            found_in = 0;
          if ( found_in ) {
            for ( jnol=inol+1; jnol<nnol; jnol++ ) {
              jn = nodes[jnol];
              if ( db_active_index( NODE_BOUNDARY, jn, version ) )
                found_jn = 1;
              else
                found_jn = 0;
              if ( found_jn ) {
                db( NODE, in, idum, coord, ldum, version, GET );
                array_move(coord, 
                  &boundary_segments_node[number_of_boundary_segments*
                  ndim*ndim+0*ndim],ndim);
                boundary_segments_node_numbers[number_of_boundary_segments*
                  ndim+0] = in;
                db( NODE, jn, idum, coord, ldum, version, GET );
                array_move(coord, 
                  &boundary_segments_node[number_of_boundary_segments*
                  ndim*ndim+1*ndim],ndim);
                boundary_segments_node_numbers[number_of_boundary_segments*
                  ndim+1] = jn;
                number_of_boundary_segments++;
              }
            }
          }
        }
      }
      else {
        assert( ndim==3 );
        for ( inol=0; inol<nnol; inol++ ) {
          in = nodes[inol];
          if ( db_active_index( NODE_BOUNDARY, in, version ) ) 
            found_in = 1;
          else
            found_in = 0;
          if ( found_in ) {
            for ( jnol=inol+1; jnol<nnol; jnol++ ) {
              jn = nodes[jnol];
              if ( db_active_index( NODE_BOUNDARY, jn, version ) ) 
                found_jn = 1;
              else
                found_jn = 0;
              if ( found_jn ) {
                for ( knol=jnol+1; knol<nnol; knol++ ) {
                  kn = nodes[knol];
                  if ( db_active_index( NODE_BOUNDARY, kn, version ) ) 
                    found_kn = 1;
                  else
                    found_kn = 0;
                  if ( found_kn ) {
                    db( NODE, in, idum, coord, ldum, version, GET );
                    array_move(coord, 
                      &boundary_segments_node[number_of_boundary_segments*
                      ndim*ndim+0*ndim],ndim);
                    boundary_segments_node_numbers[number_of_boundary_segments*
                      ndim+0] = in;
                    db( NODE, jn, idum, coord, ldum, version, GET );
                    array_move(coord,
                      &boundary_segments_node[number_of_boundary_segments*
                      ndim*ndim+1*ndim],ndim);
                    boundary_segments_node_numbers[number_of_boundary_segments*
                      ndim+1] = jn;
                    db( NODE, kn, idum, coord, ldum, version, GET );
                    array_move(coord,
                      &boundary_segments_node[number_of_boundary_segments*
                      ndim*ndim+2*ndim],ndim);
                    boundary_segments_node_numbers[number_of_boundary_segments*
                      ndim+2] = kn;
                    number_of_boundary_segments++;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if ( swit ) pri( "number_of_boundary_segments", number_of_boundary_segments );

  if ( swit ) pri( "Generate new nodes." );
  inod = 0;
  array_set( node_dof, 0., nuknwn );
  for ( kn=0; kn<number_of_nodes[2]; kn++ ) {
    for ( jn=0; jn<number_of_nodes[1]; jn++ ) {
      for ( in=0; in<number_of_nodes[0]; in++ ) {
        if ( ndim>0 ) coord[0] = min_coord[0]-FAC*delta + in*delta_adjusted[0]; 
        if ( ndim>1 ) coord[1] = min_coord[1]-FAC*delta + jn*delta_adjusted[1]; 
        if ( ndim>2 ) coord[2] = min_coord[2]-FAC*delta + kn*delta_adjusted[2];
        db( NODE, inod, idum, coord, ndim, VERSION_NEW_MESH_TMP, PUT );
        if ( nuknwn>0 ) db( NODE_DOF, inod, idum, node_dof, nuknwn, 
          VERSION_NEW_MESH_TMP, PUT );
        inod++;
      }
    }
  }

  if ( swit ) pri( "Generate new elements." );
  element = 0;
  for ( idim=0; idim<MDIM; idim++ )
    number_of_elements[idim] = scalar_imax(number_of_nodes[idim]-1,1);
  if ( swit ) pri( "number_of_elements", number_of_elements, MDIM );
  for ( kn=0; kn<number_of_elements[2]; kn++ ) {
    for ( jn=0; jn<number_of_elements[1]; jn++ ) {
      for ( in=0; in<number_of_elements[0]; in++ ) {
        if      ( ndim==1 ) {
          el[0] = -BAR2; nnol=2;
          el[1] = in+0; 
          el[2] = in+1;
        }
        else if ( ndim==2 ) {
          el[0] = -QUAD4; nnol=4;
          el[1] = (jn+0)*number_of_nodes[0]+in+0; 
          el[2] = (jn+0)*number_of_nodes[0]+in+1; 
          el[3] = (jn+1)*number_of_nodes[0]+in+0; 
          el[4] = (jn+1)*number_of_nodes[0]+in+1; 
        }
        else {
          assert( ndim==3 );
          el[0] = -HEX8; nnol=8;
          el[1] = (kn+0)*number_of_nodes[1]*number_of_nodes[0]+
                  (jn+0)*number_of_nodes[0]+in+0; 
          el[2] = (kn+0)*number_of_nodes[1]*number_of_nodes[0]+
                  (jn+0)*number_of_nodes[0]+in+1; 
          el[3] = (kn+0)*number_of_nodes[1]*number_of_nodes[0]+
                  (jn+1)*number_of_nodes[0]+in+0; 
          el[4] = (kn+0)*number_of_nodes[1]*number_of_nodes[0]+
                  (jn+1)*number_of_nodes[0]+in+1; 
          el[5] = (kn+1)*number_of_nodes[1]*number_of_nodes[0]+
                  (jn+0)*number_of_nodes[0]+in+0; 
          el[6] = (kn+1)*number_of_nodes[1]*number_of_nodes[0]+
                  (jn+0)*number_of_nodes[0]+in+1; 
          el[7] = (kn+1)*number_of_nodes[1]*number_of_nodes[0]+
                  (jn+1)*number_of_nodes[0]+in+0; 
          el[8] = (kn+1)*number_of_nodes[1]*number_of_nodes[0]+
                  (jn+1)*number_of_nodes[0]+in+1; 
        }
        length=1+nnol; db( ELEMENT, element, el, ddum, length, 
          VERSION_NEW_MESH_TMP, PUT );
        element++; 
      }
    }
  }
  mesh_has_changed( VERSION_NEW_MESH_TMP );

  if ( ndim==2 && control_new_mesh_element==-TRIA3 ||
       ndim==3 && control_new_mesh_element==-TET4 ) {
    if ( swit ) pri( "Split into simplex elements." );
    mesh_split( VERSION_NEW_MESH_TMP ); 
  }

  if ( swit ) pri( "Map." );
  db_max_index( ELEMENT, max_element, VERSION_NEW_MESH_TMP, GET );
  db_max_index( NODE, max_node, VERSION_NEW_MESH_TMP, GET );
  db_allocate_class( ELEMENT, max_element, VERSION_NEW_MESH_TMP );
  db_allocate_class( NODE, max_node, VERSION_NEW_MESH_TMP );
  map_always = 0;
  map_version_from = version;
  map_version_to = VERSION_NEW_MESH_TMP;
  parallel_sys_routine( &parallel_map_element );
  parallel_sys_routine( &parallel_map_node );

  if ( swit ) pri( "Which new nodes are inside the old mesh?" );
  map_node_inside = get_new_int( 1+max_node );
  array_set( map_node_inside, 0, 1+max_node );
  for ( inod=0; inod<=max_node; inod++ ) {
    if ( db_active_index( NODE_START_REFINED, inod, VERSION_NEW_MESH_TMP ) )
      map_node_inside[inod] = 1;
  }
  if ( swit ) pri( "map_node_inside", map_node_inside, 1+max_node );

  for ( inod=0; inod<=max_node; inod++ ) {
    if ( db_active_index( NODE, inod, VERSION_NEW_MESH_TMP ) ) {
      db( NODE, inod, idum, coord, ldum, VERSION_NEW_MESH_TMP, GET );
      db( NODE_START_REFINED, inod, idum, coord, 
        ldum, VERSION_NEW_MESH_TMP, PUT );
      if ( nuknwn>0 ) {
        db( NODE_DOF, inod, idum, node_dof, nuknwn, 
          VERSION_NEW_MESH_TMP, GET );
        db( NODE_DOF_START_REFINED, inod, idum, node_dof, nuknwn, 
          VERSION_NEW_MESH_TMP, PUT );
      }
    }
  }

  if ( swit ) pri( "Delete elements totally outside old mesh." );
  db_max_index( ELEMENT, max_element, VERSION_NEW_MESH_TMP, GET );
  for ( element=0; element<=max_element; element++ ) {
    if ( db_active_index( ELEMENT, element, VERSION_NEW_MESH_TMP ) ) {
      db( ELEMENT, element, el, ddum, length, VERSION_NEW_MESH_TMP, GET );
      nnol = length - 1; array_move( &el[1], nodes, nnol );
      any_node_inside_old_mesh = 0;
      for ( inol=0; inol<nnol; inol++ ) {
        inod = nodes[inol];
        if ( map_node_inside[inod] ) any_node_inside_old_mesh = 1;
      }
      if ( !any_node_inside_old_mesh ) {
        delete_element( element, VERSION_NEW_MESH_TMP );
      }
    }
  }

  if ( swit ) pri( "Move nodes to boundary." );
  if ( number_of_boundary_segments>0 ) {
    db_max_index( NODE, max_node, VERSION_NEW_MESH_TMP, GET );
    for ( inod=0; inod<=max_node; inod++ ) {
      if ( db_active_index( NODE, inod, VERSION_NEW_MESH_TMP ) ) {
        db( NODE, inod, idum, coord, ldum, VERSION_NEW_MESH_TMP, GET );
        if ( !map_node_inside[inod] ) {
          closest_distance = 1.e10;
          for ( isegment=0; isegment<number_of_boundary_segments; isegment++ ) {
            if      ( ndim==1 ) {
              weight[0] = 1.;
            }
            else if ( ndim==2 ) {
              in = isegment*ndim*ndim+0*ndim;
              jn = isegment*ndim*ndim+1*ndim;
              array_move( &boundary_segments_node[in], coord0, ndim );
              array_move( &boundary_segments_node[jn], coord1, ndim );
              project_point_exactly_on_line( coord, coord0, coord1, weight );
            }
            else {
              assert( ndim==3 );
              in = isegment*ndim*ndim+0*ndim;
              jn = isegment*ndim*ndim+1*ndim;
              kn = isegment*ndim*ndim+2*ndim;
              array_move( &boundary_segments_node[in], coord0, ndim );
              array_move( &boundary_segments_node[jn], coord1, ndim );
              array_move( &boundary_segments_node[kn], coord2, ndim );
              project_point_exactly_on_triangle( coord, coord0, 
                coord1, coord2, weight );
            }
            array_set( projection, 0., ndim );
            for ( inol=0; inol<ndim; inol++ ) {
              for ( idim=0; idim<ndim; idim++ ) {
                if ( scalar_dabs(weight[inol])>EPS_ISO ) {
                  projection[idim] += weight[inol]*
                    boundary_segments_node[isegment*ndim*ndim+inol*ndim+idim];
                }
              }
            }
            distance = array_distance( coord, projection, work, ndim );
            if ( distance<closest_distance ) {
              closest_segment = isegment;
              closest_distance = distance;
              array_move( projection, closest_projection, ndim );
            }
          }
          db( NODE, inod, idum, closest_projection, ndim,
             VERSION_NEW_MESH_TMP, PUT );
          db( NODE_START_REFINED, inod, idum, closest_projection, ndim,
             VERSION_NEW_MESH_TMP, PUT );
          in = boundary_segments_node_numbers[closest_segment*ndim+0];
          db( NODE_BOUNDARY, in, &node_boundary_in,
             ddum, length, version, GET );
          if ( ndim>1 ) {
            jn = boundary_segments_node_numbers[closest_segment*ndim+1];
            db( NODE_BOUNDARY, jn, &node_boundary_jn,
              ddum, length, version, GET );
          }
          if ( ndim>2 ) {
            kn = boundary_segments_node_numbers[closest_segment*ndim+2];
            db( NODE_BOUNDARY, kn, &node_boundary_kn,
              ddum, length, version, GET );
          }
          test1 = ndim==1;
          test2 = ndim==2 && node_boundary_in==node_boundary_jn;
          test3 = ndim==3 && node_boundary_in==node_boundary_jn &&
            node_boundary_in==node_boundary_kn;
          if ( test1 || test2 || test3 ) {
            db( NODE_BOUNDARY, inod, &node_boundary_in,
              ddum, length, VERSION_NEW_MESH_TMP, PUT );
          }
        }
      }
    }
  }
  db_delete( NODE_BOUNDARY, version );

    // map because of changed node positions
  map_always = 1;
  map_version_from = version;
  map_version_to = VERSION_NEW_MESH_TMP;
  parallel_sys_routine( &parallel_map_element );
  parallel_sys_routine( &parallel_map_node );

    // get rid of allocations
  delete[] boundary_segments_node;
  delete[] boundary_segments_node_numbers;
  delete[] map_node_inside;

    // store the new version
  db_version_copy( VERSION_NEW_MESH_TMP, version );
  mesh_has_changed( version );

    // clean database
  db_version_delete( VERSION_NEW_MESH_TMP );

  if ( swit ) pri( "Out routine NEW_MESH_VERSION" );

}
