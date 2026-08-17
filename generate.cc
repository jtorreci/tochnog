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

#define EPS_COORD 1.e-10

// interface_face_subdivide - subdivide a shared face into linear sub-faces.
//
// The interface elements implemented are LINEAR (quad4/prism6/hex8). When
// the two contacting elements are QUADRATIC (quad9/tet10/hex27/bar3), the
// shared face has mid-side nodes (2D: 3 nodes; 3D tria6: 6 nodes; 3D
// quad9: 9 nodes). Using a single linear interface over only the corner
// nodes would leave the mid-side nodes kinematically uncoupled, creating
// gaps. Instead we SUBDIVIDE the quadratic face into linear sub-faces
// (pattern: mesh the face with triangles/quadrilaterals) and generate one
// linear interface per sub-face, so ALL face nodes are coupled.
//
//   face_nnodes: number of shared nodes (2,3 linear; 3/6/9 quadratic)
//   node_xyz[][]: coordinates of the shared nodes (side A)
//   sub_faces[][]: filled with the node indices (into node_xyz) of each
//                  linear sub-face
//   sub_nnodes[]: number of nodes of each sub-face (2 or 3 or 4)
//   returns: number of sub-faces
//
// The shared-node ordering is arbitrary (depends on the element scan), so
// the geometry is classified by COORDINATES: a node is a mid-side node if
// it is the average of two other nodes; the centre of a quad9 is the
// average of the four corners. Sub-division rules:
//   2D, 3 nodes (quadratic edge): 2 linear segments [c0 m] [m c1]
//   3D, 6 nodes (tria6 from tet10): 4 linear triangles
//   3D, 9 nodes (quad9 from hex27): 4 linear quads
long int interface_face_subdivide( long int face_nnodes, double node_xyz[][MDIM],
  long int sub_faces[][4], long int sub_nnodes[] )

{
  long int i=0, j=0, k=0, d=0;

  // 2D linear edge (2 nodes): single sub-face
  if ( face_nnodes==2 && ndim==2 ) {
    sub_nnodes[0] = 2; sub_faces[0][0] = 0; sub_faces[0][1] = 1;
    return 1;
  }
  // 3D linear tria (3 nodes): single sub-face
  if ( face_nnodes==3 && ndim==3 ) {
    sub_nnodes[0] = 3; sub_faces[0][0] = 0; sub_faces[0][1] = 1; sub_faces[0][2] = 2;
    return 1;
  }
  // 3D linear quad (4 nodes): single sub-face
  if ( face_nnodes==4 && ndim==3 ) {
    sub_nnodes[0] = 4; sub_faces[0][0] = 0; sub_faces[0][1] = 1; sub_faces[0][2] = 2; sub_faces[0][3] = 3;
    return 1;
  }

  // 2D quadratic edge (3 nodes): the mid node is the average of the two ends.
  if ( face_nnodes==3 ) {
    long int m = -1;
    for ( i=0; i<3; i++ ) {
      j = (i+1)%3; k = (i+2)%3;
      long int ok = 1;
      for ( d=0; d<ndim && ok; d++ )
        if ( fabs( node_xyz[i][d] - 0.5*(node_xyz[j][d]+node_xyz[k][d]) ) > EPS_COORD )
          ok = 0;
      if ( ok ) { m = i; break; }
    }
    if ( m<0 ) return 0;
    long int c0 = (m+1)%3, c1 = (m+2)%3;
    // order each linear segment so its first node has the smaller
    // coordinate along the edge direction -> tangent points +, normal
    // points from side 2 towards side 1 (compression = positive strain).
    sub_nnodes[0]=2; sub_nnodes[1]=2;
    long int seg0a = c0, seg0b = m, seg1a = m, seg1b = c1;
    // compare along the dominant edge direction
    long int d0 = 0;
    double span = -1.;
    for ( d=0; d<ndim; d++ ) {
      double s = fabs( node_xyz[c1][d] - node_xyz[c0][d] );
      if ( s>span ) { span = s; d0 = d; }
    }
    if ( node_xyz[seg0b][d0] < node_xyz[seg0a][d0] ) { long int t=seg0a; seg0a=seg0b; seg0b=t; }
    if ( node_xyz[seg1b][d0] < node_xyz[seg1a][d0] ) { long int t=seg1a; seg1a=seg1b; seg1b=t; }
    sub_faces[0][0]=seg0a; sub_faces[0][1]=seg0b;
    sub_faces[1][0]=seg1a; sub_faces[1][1]=seg1b;
    return 2;
  }

  // 3D quadratic tria (6 nodes, from tet10). Identify the 3 corners
  // (not an average of another pair) and the 3 mid-edge nodes.
  if ( face_nnodes==6 ) {
    long int corner[3], ncorner=0, mid[3], nmid=0;
    for ( i=0; i<6; i++ ) {
      long int is_mid = 0;
      for ( j=0; j<6 && !is_mid; j++ ) {
        if ( j==i ) continue;
        for ( k=j+1; k<6 && !is_mid; k++ ) {
          if ( k==i ) continue;
          long int ok = 1;
          for ( d=0; d<3 && ok; d++ )
            if ( fabs( node_xyz[i][d] - 0.5*(node_xyz[j][d]+node_xyz[k][d]) ) > EPS_COORD )
              ok = 0;
          if ( ok ) is_mid = 1;
        }
      }
      if ( is_mid ) mid[nmid++] = i; else corner[ncorner++] = i;
    }
    if ( ncorner!=3 || nmid!=3 ) return 0;
    // corner i is opposite mid[i] (the mid of the two edges from corner i
    // to the other two corners). Build the corner triangles so the
    // mid-edges pair up: triangle i = (corner_i, mid_i, mid_k) with
    // k such that mid_k connects corner_i's other two neighbours.
    long int mid_opp[3][2];   // for corner i, the two mid nodes adjacent
    for ( i=0; i<3; i++ ) {
      long int a = (i+1)%3, b = (i+2)%3;
      mid_opp[i][0] = -1; mid_opp[i][1] = -1;
      for ( j=0; j<3; j++ ) {
        // mid[j] is adjacent to corner i if it is the average of
        // corner[i] and one of the other corners
        long int ok = 1;
        for ( d=0; d<3 && ok; d++ )
          if ( fabs( node_xyz[mid[j]][d] - 0.5*(node_xyz[corner[i]][d]+node_xyz[corner[a]][d]) ) > EPS_COORD )
            ok = 0;
        if ( ok ) { mid_opp[i][0] = mid[j]; break; }
      }
      for ( j=0; j<3; j++ ) {
        long int ok = 1;
        for ( d=0; d<3 && ok; d++ )
          if ( fabs( node_xyz[mid[j]][d] - 0.5*(node_xyz[corner[i]][d]+node_xyz[corner[b]][d]) ) > EPS_COORD )
            ok = 0;
        if ( ok ) { mid_opp[i][1] = mid[j]; break; }
      }
      if ( mid_opp[i][0]<0 || mid_opp[i][1]<0 ) return 0;
    }
    // central triangle = the 3 mid nodes
    sub_nnodes[0]=3; sub_faces[0][0]=mid[0]; sub_faces[0][1]=mid[1]; sub_faces[0][2]=mid[2];
    for ( i=0; i<3; i++ ) {
      sub_nnodes[1+i]=3;
      sub_faces[1+i][0]=corner[i];
      sub_faces[1+i][1]=mid_opp[i][0];
      sub_faces[1+i][2]=mid_opp[i][1];
    }
    return 4;
  }

  // 3D quadratic quad (9 nodes, from hex27). 4 corners, 4 mid-edge, 1 centre.
  if ( face_nnodes==9 ) {
    long int corner[4], ncorner=0, mid[4], nmid=0, centre=-1;
    for ( i=0; i<9; i++ ) {
      long int is_mid = 0;
      for ( j=0; j<9 && !is_mid; j++ ) {
        if ( j==i ) continue;
        for ( k=j+1; k<9 && !is_mid; k++ ) {
          if ( k==i ) continue;
          long int ok = 1;
          for ( d=0; d<3 && ok; d++ )
            if ( fabs( node_xyz[i][d] - 0.5*(node_xyz[j][d]+node_xyz[k][d]) ) > EPS_COORD )
              ok = 0;
          if ( ok ) is_mid = 1;
        }
      }
      if ( is_mid ) mid[nmid++] = i;
      else corner[ncorner++] = i;
    }
    if ( ncorner!=4 || nmid!=4 ) return 0;
    // centre = the node that is the average of the 4 corners
    double cx=0., cy=0., cz=0.;
    for ( i=0; i<4; i++ ) { cx+=node_xyz[corner[i]][0]; cy+=node_xyz[corner[i]][1]; cz+=node_xyz[corner[i]][2]; }
    cx/=4.; cy/=4.; cz/=4.;
    for ( i=0; i<9; i++ ) {
      long int is_corner=0;
      for ( j=0; j<4 && !is_corner; j++ ) if ( i==corner[j] ) is_corner=1;
      if ( !is_corner ) {
        long int is_mid=0;
        for ( j=0; j<4 && !is_mid; j++ ) if ( i==mid[j] ) is_mid=1;
        if ( !is_mid ) { centre=i; break; }
      }
    }
    if ( centre<0 ) return 0;
    // order the corners around the centre by angle
    long int ordered[4]; ordered[0]=corner[0];
    // mid edge between corner[i] and corner[j]
    long int mid_ij[4][4]; for ( i=0;i<4;i++) for (j=0;j<4;j++) mid_ij[i][j]=-1;
    for ( i=0; i<4; i++ ) {
      for ( j=i+1; j<4; j++ ) {
        for ( k=0; k<4; k++ ) {
          long int ok = 1;
          for ( d=0; d<3 && ok; d++ )
            if ( fabs( node_xyz[mid[k]][d] - 0.5*(node_xyz[corner[i]][d]+node_xyz[corner[j]][d]) ) > EPS_COORD )
              ok = 0;
          if ( ok ) { mid_ij[i][j]=mid_ij[j][i]=mid[k]; break; }
        }
      }
    }
    // walk the perimeter: from corner 0, the two adjacent mids connect to
    // two of the other corners; pick the next corner that shares a mid with
    // corner 0, then continue.
    for ( i=1; i<4; i++ ) {
      long int prev = ordered[i-1];
      long int nxt = -1;
      for ( j=0; j<4 && nxt<0; j++ ) {
        if ( j==prev ) continue;
        long int used = 0;
        for ( k=0; k<i; k++ ) if ( ordered[k]==j ) used=1;
        if ( used ) continue;
        if ( mid_ij[prev][j]>=0 ) nxt = j;
      }
      if ( nxt<0 ) return 0;
      ordered[i] = nxt;
    }
    // sub-faces: (corner_i, mid_i->i+1, centre, mid_i-1->i)
    for ( i=0; i<4; i++ ) {
      long int a = ordered[i], b = ordered[(i+1)%4], c = ordered[(i+3)%4];
      sub_nnodes[i]=4;
      sub_faces[i][0]=a;
      sub_faces[i][1]=mid_ij[a][b];
      sub_faces[i][2]=centre;
      sub_faces[i][3]=mid_ij[c][a];
    }
    return 4;
  }

  return 0;
}

void generate_spring( long int icontrol )

{
  long int j=0, n=0, inod=0, jnod=0, max_node=0, max_element=0, 
    element_group=0, in_geometry=0, length=0, swit=0, ldum=0, 
    correct_elements=0, length_node_element_inod=0, length_node_element_jnod=0, 
    control_mesh_generate_contactspring_element_specified=0, 
    iel=0, element=0, name=0, element_name0=0, element_name1=0,
    element0_in_node_element_inod=0, element0_in_node_element_jnod=0,
    element1_in_node_element_inod=0, element1_in_node_element_jnod=0,
    control_mesh_generate_spring[3], el[1+MNOL],
    control_mesh_generate_contactspring_element[2], 
    geometry_entity[2], *in_geometry_list=NULL, 
    *node_element_inod=NULL, *node_element_jnod=NULL;
  double distance=0., rdum=0., ddum[MDIM], *coordi=NULL, *coordj=NULL;
  long int zero=0, mnolnuknwn=npointmax*nuknwn, idum[1]={0}, 
  	length_nei=1+npointmax*ndim+npointmax+2;
  double *tmp_element_dof=NULL, *dworknei=NULL;
  tmp_element_dof = get_new_dbl(mnolnuknwn);
  dworknei = get_new_dbl(length_nei);

  array_set( control_mesh_generate_contactspring_element, -ALL, 2 );
  array_set( dworknei, 0, length_nei );
  db_highest_index( ELEMENT, max_element, VERSION_NORMAL );
  db_max_index( NODE_START_REFINED, max_node, VERSION_NORMAL, GET );

  if ( db_active_index(CONTROL_MESH_GENERATE_SPRING1,icontrol,VERSION_NORMAL) ) {
    swit = set_swit(-1,-1,"generate_spring");
    if ( swit ) pri( "In routine GENERATE_SPRING." );
    db( CONTROL_MESH_GENERATE_SPRING1, icontrol, control_mesh_generate_spring, 
      ddum, ldum, VERSION_NORMAL, GET );
    element_group = control_mesh_generate_spring[0];
    array_move( &control_mesh_generate_spring[1], geometry_entity, 2 );
    for ( inod=0; inod<=max_node; inod++ ) {
      if ( db_active_index( NODE_START_REFINED, inod, VERSION_NORMAL ) ) {
        geometry( inod, ddum, geometry_entity, in_geometry, rdum, ddum, rdum,
          ddum, NODE_START_REFINED, PROJECT_EXACT, VERSION_NORMAL );
        if ( in_geometry ) {
          max_element++;
          el[0] = -SPRING1;
          el[1] = inod;
          length = 2;
          db( ELEMENT, max_element, el, ddum, length, 
            VERSION_NORMAL, PUT );
          length = 1;
          db( ELEMENT_GROUP, max_element, &element_group, ddum, length, 
            VERSION_NORMAL, PUT );
          length = 1;
          db( ELEMENT_MACRO_GENERATE, max_element, &icontrol, 
            ddum, length, VERSION_NORMAL, PUT );
          db( ELEMENT_DOF, max_element, idum, tmp_element_dof, mnolnuknwn, VERSION_NORMAL, PUT );
          db( ELEMENT_DOF_INITIALISED, max_element, &zero, ddum, length, VERSION_NORMAL, PUT );
          db( NONLOCAL_ELEMENT_INFO, max_element, idum, dworknei, length_nei, VERSION_NORMAL, PUT );		
        }
      }
    }

    mesh_has_changed( VERSION_NORMAL );
    if ( swit ) pri( "Out routine GENERATE_SPRING." );
  }


  if ( db_active_index(CONTROL_MESH_GENERATE_SPRING2,icontrol,VERSION_NORMAL) ) {
    swit = set_swit(-1,-1,"generate_spring");
    if ( swit ) pri( "In routine GENERATE_SPRING." );
    db( CONTROL_MESH_GENERATE_SPRING2, icontrol, control_mesh_generate_spring, ddum, ldum, 
      VERSION_NORMAL, GET );
    in_geometry_list = get_new_int(1+max_node);
    element_group = control_mesh_generate_spring[0];
    array_move( &control_mesh_generate_spring[1], geometry_entity, 2 );
    for ( inod=0; inod<=max_node; inod++ ) {
      if ( db_active_index( NODE_START_REFINED, inod, VERSION_NORMAL ) ) {
        geometry( inod, ddum, geometry_entity, in_geometry, rdum, ddum, rdum,
          ddum, NODE_START_REFINED, PROJECT_EXACT, VERSION_NORMAL );
        if ( in_geometry ) {
          in_geometry_list[n] = inod;     
          coordi = db_dbl( NODE_START_REFINED, inod, VERSION_NORMAL );
          for ( j=0; j<n; j++ ) {             
            jnod = in_geometry_list[j];
            coordj = db_dbl( NODE_START_REFINED, jnod, VERSION_NORMAL );
            distance = array_distance( coordi, coordj, ddum, ndim );
            if ( distance < EPS_COORD ) {
              max_element++;
              el[0] = -SPRING2;
              el[1] = inod;
              el[2] = jnod;
              length = 3;
              db( ELEMENT, max_element, el, ddum, length, 
                VERSION_NORMAL, PUT );
              length = 1;
              db( ELEMENT_GROUP, max_element, &element_group, ddum, length, 
                VERSION_NORMAL, PUT );
              length = 1;
              db( ELEMENT_MACRO_GENERATE, max_element, &icontrol, 
                ddum, length, VERSION_NORMAL, PUT );
              db( ELEMENT_DOF, max_element, idum, tmp_element_dof, mnolnuknwn, VERSION_NORMAL, PUT );
              db( ELEMENT_DOF_INITIALISED, max_element, &zero, ddum, length, VERSION_NORMAL, PUT );
              db( NONLOCAL_ELEMENT_INFO, max_element, idum, dworknei, length_nei, VERSION_NORMAL, PUT );		
            }
          }
          n++;
        }
      }
    }
    delete[] in_geometry_list;

    mesh_has_changed( VERSION_NORMAL );
    if ( swit ) pri( "Out routine GENERATE_SPRING." );
  }

  if ( db_active_index(CONTROL_MESH_GENERATE_CONTACTSPRING,icontrol,VERSION_NORMAL) ) {
    swit = set_swit(-1,-1,"generate_spring");
    if ( swit ) pri( "In routine GENERATE_SPRING." );
    db( CONTROL_MESH_GENERATE_CONTACTSPRING, icontrol, 
      control_mesh_generate_spring, ddum, ldum, VERSION_NORMAL, GET );
    if ( db( CONTROL_MESH_GENERATE_CONTACTSPRING_ELEMENT, icontrol, 
        control_mesh_generate_contactspring_element, ddum, ldum, 
        VERSION_NORMAL, GET_IF_EXISTS ) ) {
      control_mesh_generate_contactspring_element_specified = 1;
      element_name0 = control_mesh_generate_contactspring_element[0];
      element_name1 = control_mesh_generate_contactspring_element[1];
      length = 1+max_element;
      node_element_inod = get_new_int(1+max_element);
      node_element_jnod = get_new_int(1+max_element);
    }
    in_geometry_list = get_new_int(1+max_node);
    element_group = control_mesh_generate_spring[0];
    array_move( &control_mesh_generate_spring[1], geometry_entity, 2 );
    for ( inod=0; inod<=max_node; inod++ ) {
      if ( db_active_index( NODE_START_REFINED, inod, VERSION_NORMAL ) ) {
        geometry( inod, ddum, geometry_entity, in_geometry, rdum, ddum, rdum,
          ddum, NODE_START_REFINED, PROJECT_EXACT, VERSION_NORMAL );
        if ( in_geometry ) {
          in_geometry_list[n] = inod;     
          coordi = db_dbl( NODE_START_REFINED, inod, VERSION_NORMAL );
          for ( j=0; j<n; j++ ) {             
            jnod = in_geometry_list[j];
            coordj = db_dbl( NODE_START_REFINED, jnod, VERSION_NORMAL );
            distance = array_distance( coordi, coordj, ddum, ndim );
            correct_elements = 1;
            if ( control_mesh_generate_contactspring_element_specified ) {
              db( NODE_ELEMENT, inod, node_element_inod, ddum, 
                length_node_element_inod, VERSION_NORMAL, GET );
              db( NODE_ELEMENT, jnod, node_element_jnod, ddum, 
                length_node_element_jnod, VERSION_NORMAL, GET );
              element0_in_node_element_inod = 0;
              element1_in_node_element_inod = 0;
              for ( iel=0; iel<length_node_element_inod; iel++ ) {
                element = node_element_inod[iel];
                db( ELEMENT, element, el, ddum, ldum, VERSION_NORMAL, GET );
                name = el[0];
                if ( name==element_name0 ) element0_in_node_element_inod = 1;
                if ( name==element_name1 ) element1_in_node_element_inod = 1;
              }
              element0_in_node_element_jnod = 0;
              element1_in_node_element_jnod = 0;
              for ( iel=0; iel<length_node_element_jnod; iel++ ) {
                element = node_element_jnod[iel];
                db( ELEMENT, element, el, ddum, ldum, VERSION_NORMAL, GET );
                name = el[0];
                if ( name==element_name0 ) element0_in_node_element_jnod = 1;
                if ( name==element_name1 ) element1_in_node_element_jnod = 1;
              }
              if      ( element0_in_node_element_inod && element1_in_node_element_jnod )
                correct_elements = 1;
              else if ( element0_in_node_element_jnod && element1_in_node_element_inod )
                correct_elements = 1;
              else
                correct_elements = 0;
            }
            if ( distance<EPS_COORD && correct_elements ) {
              max_element++;
              el[0] = -CONTACTSPRING;
              el[1] = inod;
              el[2] = jnod;
              length = 3;
              db( ELEMENT, max_element, el, ddum, length, 
                VERSION_NORMAL, PUT );
              length = 1;
              db( ELEMENT_GROUP, max_element, &element_group, ddum, length, 
                VERSION_NORMAL, PUT );
              length = 1;
              db( ELEMENT_MACRO_GENERATE, max_element, &icontrol, 
                ddum, length, VERSION_NORMAL, PUT );
              db( ELEMENT_DOF, max_element, idum, tmp_element_dof, mnolnuknwn, VERSION_NORMAL, PUT );
              db( ELEMENT_DOF_INITIALISED, max_element, &zero, ddum, length, VERSION_NORMAL, PUT );
              db( NONLOCAL_ELEMENT_INFO, max_element, idum, dworknei, length_nei, VERSION_NORMAL, PUT );		
            }
          }
          n++;
        }
      }
    }
    delete[] tmp_element_dof;
    delete[] dworknei;
    delete[] in_geometry_list;
    if ( control_mesh_generate_contactspring_element_specified ) {
      delete[] node_element_inod;
      delete[] node_element_jnod;
    }

    mesh_has_changed( VERSION_NORMAL );
    if ( swit ) pri( "Out routine GENERATE_SPRING." );
  }

}

void generate_beam_truss( long int icontrol, long int task )

{
  long int jn=0, inod=0, jnod=0, max_node_old=0, max_node=0, max_element=0, 
    element_group=0, in_geometry=0, igenerated=0, ngenerated=0, mgenerated=0, 
    length_node_node=0, already_generated=0, swit=0, length=0, loose=-NO, ldum=0,
    node_macro_generate=0, length_macro=0, idum[1], control_mesh_generate[3], 
    geometry_entity[2], el[1+MNOL], macro[DATA_ITEM_SIZE],
    *node_node=NULL, *in_geometry_list=0, 
    *generated_list=NULL, *new_node_list=NULL;
  double rdum=0., ddum[MDIM], coord[MDIM], node_dof[MUKNWN];
  long int zero=0, mnolnuknwn=npointmax*nuknwn, length_nei=1+npointmax*ndim+npointmax+2;
  double *tmp_element_dof=NULL, *dworknei=NULL;
  tmp_element_dof = get_new_dbl(mnolnuknwn);
  dworknei = get_new_dbl(length_nei);
  array_set(dworknei, 0, length_nei);

  if     ( task==TRUSS ) {
    if ( db_active_index(CONTROL_MESH_GENERATE_TRUSS,icontrol,VERSION_NORMAL) )
      db( CONTROL_MESH_GENERATE_TRUSS, icontrol, control_mesh_generate, ddum, 
        ldum, VERSION_NORMAL, GET );
    else
      return;
  }
  else if ( task==TRUSSBEAM ) {
    if ( db_active_index(CONTROL_MESH_GENERATE_TRUSSBEAM,icontrol,VERSION_NORMAL) )
      db( CONTROL_MESH_GENERATE_TRUSSBEAM, icontrol, control_mesh_generate, ddum, 
        ldum, VERSION_NORMAL, GET );
    else
      return;
  }
  else {
    assert( task==BEAM );
    if ( db_active_index(CONTROL_MESH_GENERATE_BEAM,icontrol,VERSION_NORMAL) )
      db( CONTROL_MESH_GENERATE_BEAM, icontrol, control_mesh_generate, ddum, 
        ldum, VERSION_NORMAL, GET );
    else
      return;
  }

  swit = set_swit(-1,-1,"generate_beam_truss");
  if ( swit ) pri( "In routine GENERATE_BEAM_TRUSS." );

  db( CONTROL_MESH_GENERATE_TRUSS_BEAM_LOOSE, icontrol, &loose, 
    ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_MESH_GENERATE_TRUSS_BEAM_MACRO, icontrol, macro, 
    ddum, length_macro, VERSION_NORMAL, GET_IF_EXISTS );

  db_highest_index( ELEMENT, max_element, VERSION_NORMAL );
  db_max_index( NODE_START_REFINED, max_node_old, VERSION_NORMAL, GET );

  element_group = control_mesh_generate[0];
  array_move( &control_mesh_generate[1], geometry_entity, 2 );

  length = db_data_length(NODE_NODE);
  node_node = get_new_int(length);

    // list for nodes in geometry
  in_geometry_list = get_new_int(1+max_node_old);
  array_set( in_geometry_list, 0, (1+max_node_old) );

    // list for generated beams/trusses, 
  mgenerated = 5*ndim*(1+max_node_old);
  generated_list = get_new_int(mgenerated*2);
  array_set( generated_list, -1, (mgenerated*2) );

    // list for new node numbers (for generate with contact spring)
  new_node_list = get_new_int(1+max_node_old);
  array_set( new_node_list, -1, (1+max_node_old) );

  for ( inod=0; inod<=max_node_old; inod++ ) {
    if ( db_active_index( NODE_START_REFINED, inod, VERSION_NORMAL ) ) {
      geometry( inod, ddum, geometry_entity, in_geometry, rdum, ddum, rdum,
        ddum, NODE_START_REFINED, PROJECT_EXACT, VERSION_NORMAL );
      if ( in_geometry ) in_geometry_list[inod] = 1;     
    }
  }

  if ( db_active_index( NODE, 0, VERSION_NORMAL ) ) {
    pri( "Error: node number 0 not allowed if you generate trusses, beams, or so." );
    exit(TN_EXIT_STATUS);
  }

  max_node = max_node_old;
  for ( inod=0; inod<=max_node_old; inod++ ) {
    if ( db_active_index( NODE_NODE, inod, VERSION_NORMAL ) ) {
      if ( in_geometry_list[inod] ) {
        node_macro_generate = -ALL;
        db( NODE_MACRO_GENERATE, inod, &node_macro_generate, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
        if ( length_macro==0 || array_member(macro,node_macro_generate,length_macro,ldum) ) {
          db( NODE_NODE, inod, node_node, ddum, length_node_node, VERSION_NORMAL, GET );
          for ( jn=0; jn<length_node_node; jn++ ) {
            jnod = node_node[jn];
            if ( jnod>=0 ) {
              if ( in_geometry_list[jnod] ) {
                node_macro_generate = -ALL;
                db( NODE_MACRO_GENERATE, jnod, &node_macro_generate, ddum, ldum, 
                  VERSION_NORMAL, GET_IF_EXISTS );
                if ( length_macro==0 || array_member(macro,node_macro_generate,length_macro,ldum) ) {
                  already_generated = 0;
                  for ( igenerated=0; igenerated<ngenerated; igenerated++ ) {
                    if ( ( generated_list[igenerated*2+0]==inod && 
                           generated_list[igenerated*2+1]==jnod ) ||
                         ( generated_list[igenerated*2+0]==jnod && 
                           generated_list[igenerated*2+1]==inod ) )
                      already_generated = 1;
                  }         
                  if ( !already_generated ) {
                    ngenerated++;
                    if ( ngenerated>mgenerated ) {
                      pri( "Error: mgenerated too small in routine generate_beam_truss. ");
                      exit(TN_EXIT_STATUS );
                    }
                    generated_list[(ngenerated-1)*2+0]=inod;
                    generated_list[(ngenerated-1)*2+1]=jnod;
                    if      ( task==TRUSS ) 
                      el[0] = -TRUSS;
                    else if ( task==TRUSSBEAM ) 
                      el[0] = -TRUSSBEAM;
                    else {
                      assert( task==BEAM );
                      el[0] = -BEAM;
                    }
                    if ( loose==-YES ) {
                         // generate beam / truss
                      if ( new_node_list[inod]>=0 ) {
                        el[1] = new_node_list[inod];
                      }
                      else {
                        max_node++;
                        el[1] = max_node;
                        new_node_list[inod] = max_node;
                        db( NODE, inod, idum, coord, ldum, VERSION_NORMAL, GET );
                        db( NODE, max_node, idum, coord, ldum, VERSION_NORMAL, PUT );
                        db( NODE_START_REFINED, inod, idum, coord, ldum, VERSION_NORMAL, GET );
                        db( NODE_START_REFINED, max_node, idum, coord, ldum, VERSION_NORMAL, PUT );
                        length = 1; db( NODE_MACRO_GENERATE, max_node, &icontrol, ddum, length, VERSION_NORMAL, PUT );
                        db( NODE_DOF, inod, idum, node_dof, ldum, VERSION_NORMAL, GET );
                        db( NODE_DOF, max_node, idum, node_dof, ldum, VERSION_NORMAL, PUT );
                        db( NODE_DOF_START_REFINED, inod, idum, node_dof, ldum, VERSION_NORMAL, GET );
                        db( NODE_DOF_START_REFINED, max_node, idum, node_dof, ldum, VERSION_NORMAL, PUT );
                      }
                      if ( new_node_list[jnod]>=0 ) {
                        el[2] = new_node_list[jnod];
                      }
                      else {
                        max_node++;
                        el[2] = max_node;
                        new_node_list[jnod] = max_node;
                        db( NODE, jnod, idum, coord, ldum, VERSION_NORMAL, GET );
                        db( NODE, max_node, idum, coord, ldum, VERSION_NORMAL, PUT );
                        db( NODE_START_REFINED, jnod, idum, coord, ldum, VERSION_NORMAL, GET );
                        db( NODE_START_REFINED, max_node, idum, coord, ldum, VERSION_NORMAL, PUT );
                        length = 1; db( NODE_MACRO_GENERATE, max_node, &icontrol, ddum, length, VERSION_NORMAL, PUT );
                        db( NODE_DOF, jnod, idum, node_dof, ldum, VERSION_NORMAL, GET );
                        db( NODE_DOF, max_node, idum, node_dof, ldum, VERSION_NORMAL, PUT );
                        db( NODE_DOF_START_REFINED, jnod, idum, node_dof, ldum, VERSION_NORMAL, GET );
                        db( NODE_DOF_START_REFINED, max_node, idum, node_dof, ldum, VERSION_NORMAL, PUT );
                      }
                      length = 3;
                      max_element++;
                      db( ELEMENT, max_element, el, ddum, length, VERSION_NORMAL, PUT );
                      length = 1;
                      db( ELEMENT_GROUP, max_element, &element_group, ddum, length, 
                        VERSION_NORMAL, PUT );
                      length = 1;
                      db( ELEMENT_MACRO_GENERATE, max_element, &icontrol, ddum, length, 
                        VERSION_NORMAL, PUT );
    	              db( ELEMENT_DOF, max_element, idum, tmp_element_dof, mnolnuknwn, VERSION_NORMAL, PUT );
                      db( ELEMENT_DOF_INITIALISED, max_element, &zero, ddum, length, VERSION_NORMAL, PUT );
                      db( NONLOCAL_ELEMENT_INFO, max_element, idum, dworknei, length_nei, VERSION_NORMAL, PUT );		
                    }
                    else {
                      if      ( task==TRUSS ) 
                        el[0] = -TRUSS;
                      else if ( task==TRUSSBEAM ) 
                        el[0] = -TRUSSBEAM;
                      else {
                        assert( task==BEAM );
                        el[0] = -BEAM;
                      }
                      el[1] = inod;
                      el[2] = jnod;
                      length = 3;
                      max_element++;
                      db( ELEMENT, max_element, el, ddum, length, VERSION_NORMAL, PUT );
                      length = 1;
                      db( ELEMENT_GROUP, max_element, &element_group, ddum, length, 
                        VERSION_NORMAL, PUT );
                      length = 1;
                      db( ELEMENT_MACRO_GENERATE, max_element, &icontrol, ddum, length, 
                        VERSION_NORMAL, PUT );
    	              db( ELEMENT_DOF, max_element, idum, tmp_element_dof, mnolnuknwn, VERSION_NORMAL, PUT );
                      db( ELEMENT_DOF_INITIALISED, max_element, &zero, ddum, length, VERSION_NORMAL, PUT );
                      db( NONLOCAL_ELEMENT_INFO, max_element, idum, dworknei, length_nei, VERSION_NORMAL, PUT );		
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  delete[] in_geometry_list;
  delete[] generated_list;
  delete[] new_node_list;
  delete[] node_node;
  delete[] tmp_element_dof;
  delete[] dworknei;

  mesh_has_changed( VERSION_NORMAL );
  if ( swit ) pri( "Out routine GENERATE_BEAM_TRUSS." );

}

// generate_interface - control_mesh_generate_interface (Carril B).
//
// Generates interface elements between two element groups that share a
// common face (spatially coincident nodes, e.g. duplicated nodes between
// two blocks). The record syntax is:
//
//   control_mesh_generate_interface index
//     eg0 eg00 eg01 eg1 eg10 eg11 ...
//
// For each triple (eg_i, eg_a, eg_b): an interface element is generated
// for every element pair (one in group eg_a, one in group eg_b) that
// shares a complete face. The interface element is assigned to group
// eg_i. The generated element type follows the shared face:
//   - 2D, 2 shared nodes  : -quad4  {nA0 nA1 nB0 nB1}
//   - 3D, 3 shared nodes  : -prism6 {nA0 nA1 nA2 nB0 nB1 nB2}
//   - 3D, 4 shared nodes  : -hex8   {nA0..nA3 nB0..nB3}
// control_mesh_generate_interface_geometry restricts generation to the
// given geometry (all shared nodes must be inside it).
void generate_interface( long int icontrol )

{
  long int i=0, k=0, iel=0, jel=0, inol=0, jnol=0, isub=0, nsub=0,
    max_element=0, max_node=0, element_group=0, in_geometry=0, ldum=0,
    swit=0, length=0, nshared=0, sharedA[MNOL], sharedB[MNOL],
    length_geometry=0, length_gen=0,
    new_name=0, sub_faces[8][4], sub_nnodes[8],
    idum[1], *elA=NULL, *elB=NULL, *nodesA=NULL, *nodesB=NULL,
    *geometry_entity=NULL, *gen=NULL;
  double rdum=0., ddum[MDIM], *cA=NULL, *cB=NULL, node_xyz[MNOL][MDIM];
  long int zero=0, mnolnuknwn=npointmax*nuknwn,
    length_nei=1+npointmax*ndim+npointmax+2;
  double *tmp_element_dof=NULL, *dworknei=NULL;

  swit = set_swit(-1,-1,"generate_interface");
  if ( swit ) pri( "In routine GENERATE_INTERFACE." );

  if ( !db_active_index( CONTROL_MESH_GENERATE_INTERFACE, icontrol, VERSION_NORMAL ) )
    return;

  tmp_element_dof = get_new_dbl(mnolnuknwn);
  dworknei = get_new_dbl(length_nei);
  array_set( dworknei, 0, length_nei );

  gen = get_new_int(DATA_ITEM_SIZE);
  db( CONTROL_MESH_GENERATE_INTERFACE, icontrol, gen, ddum, length_gen,
    VERSION_NORMAL, GET );

  // method: method_select / method_generate. Default = element_group.
  // -element_geometry selects by element_geometry and/or generates an
  // element_geometry record for the interface element (manual 6.192).
  long int method_select=0, method_generate=0;
  {
    long int method[2] = {0, 0};
    if ( db( CONTROL_MESH_GENERATE_INTERFACE_METHOD, icontrol, method,
        ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS ) ) {
      method_select  = method[0];
      method_generate = method[1];
    }
  }

  geometry_entity = get_new_int(DATA_ITEM_SIZE);
  if ( db_active_index( CONTROL_MESH_GENERATE_INTERFACE_GEOMETRY, icontrol,
      VERSION_NORMAL ) ) {
    db( CONTROL_MESH_GENERATE_INTERFACE_GEOMETRY, icontrol, geometry_entity,
      ddum, length_geometry, VERSION_NORMAL, GET );
  }

  elA = get_new_int(MAXIMUM_NODE+1);
  elB = get_new_int(MAXIMUM_NODE+1);
  nodesA = get_new_int(MAXIMUM_NODE);
  nodesB = get_new_int(MAXIMUM_NODE);

  db_highest_index( ELEMENT, max_element, VERSION_NORMAL );
  db_max_index( NODE_START_REFINED, max_node, VERSION_NORMAL, GET );
  long int max_element_old = max_element;

  // the record is a list of triples (eg_i, eg_a, eg_b)
  for ( i=0; i+2<length_gen; i+=3 ) {
    long int eg_iface = gen[i];
    long int eg_a = gen[i+1];
    long int eg_b = gen[i+2];

    for ( iel=0; iel<=max_element_old; iel++ ) {
      if ( !db_active_index( ELEMENT, iel, VERSION_NORMAL ) ) continue;
      // select by element_group (default) or by element_geometry (method)
      long int grA = 0;
      if ( method_select==-ELEMENT_GEOMETRY )
        db( ELEMENT_GEOMETRY, iel, &grA, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
      else
        db( ELEMENT_GROUP, iel, &grA, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
      if ( grA!=eg_a ) continue;
      // already generated an interface for this element in a previous step
      long int iface_done = -1;
      db( ELEMENT_MACRO_GENERATE, iel, &iface_done, ddum, ldum,
        VERSION_NORMAL, GET_IF_EXISTS );
      if ( iface_done==icontrol ) continue;
      db( ELEMENT, iel, elA, ddum, length, VERSION_NORMAL, GET );
      long int nnolA = length - 1;
      array_move( &elA[1], nodesA, nnolA );

      for ( jel=0; jel<=max_element_old; jel++ ) {
        if ( jel==iel ) continue;
        if ( !db_active_index( ELEMENT, jel, VERSION_NORMAL ) ) continue;
        long int grB = 0;
        if ( method_select==-ELEMENT_GEOMETRY )
          db( ELEMENT_GEOMETRY, jel, &grB, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
        else
          db( ELEMENT_GROUP, jel, &grB, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
        if ( grB!=eg_b ) continue;
        db( ELEMENT, jel, elB, ddum, length, VERSION_NORMAL, GET );
        long int nnolB = length - 1;
        array_move( &elB[1], nodesB, nnolB );

        // shared face = node pairs of A and B with coincident coordinates
        nshared = 0;
        for ( inol=0; inol<nnolA && nshared<MNOL; inol++ ) {
          cA = db_dbl( NODE_START_REFINED, nodesA[inol], VERSION_NORMAL );
          for ( jnol=0; jnol<nnolB; jnol++ ) {
            cB = db_dbl( NODE_START_REFINED, nodesB[jnol], VERSION_NORMAL );
            if ( array_distance( cA, cB, ddum, ndim ) < EPS_COORD ) {
              sharedA[nshared] = nodesA[inol];
              sharedB[nshared] = nodesB[jnol];
              nshared++;
              break;
            }
          }
        }
        if ( nshared<2 ) continue;

        // geometry restriction: every shared node must be inside
        if ( length_geometry>0 ) {
          for ( k=0; k<nshared; k++ ) {
            in_geometry = 0;
            geometry( sharedA[k], ddum, geometry_entity, in_geometry, rdum,
              ddum, rdum, ddum, NODE_START_REFINED, PROJECT_EXACT,
              VERSION_NORMAL );
            if ( !in_geometry ) break;
          }
          if ( !in_geometry ) continue;
        }

        // avoid duplicate generation of the symmetric pair (jel,iel)
        if ( iel>jel ) continue;

        // collect the shared-node coordinates (side A) for the subdivision
        for ( k=0; k<nshared; k++ ) {
          cA = db_dbl( NODE_START_REFINED, sharedA[k], VERSION_NORMAL );
          for ( long int d=0; d<ndim; d++ ) node_xyz[k][d] = cA[d];
        }

        // subdivide the (possibly quadratic) face into linear sub-faces.
        // A quadratic face (3 nodes in 2D, 6/9 in 3D) is split so ALL
        // face nodes (including mid-side) get coupled by the interface.
        nsub = interface_face_subdivide( nshared, node_xyz, sub_faces,
          sub_nnodes );
        if ( nsub<=0 ) continue;

        // generate one linear interface element per sub-face
        for ( isub=0; isub<nsub; isub++ ) {
          long int nnsub = sub_nnodes[isub];
          if      ( nnsub==2 ) new_name = -QUAD4;
          else if ( nnsub==3 ) new_name = -PRISM6;
          else                 new_name = -HEX8;
          max_element++;
          elB[0] = new_name;
          length = 1 + 2*nnsub;
          for ( k=0; k<nnsub; k++ ) {
            elB[1+k]           = sharedA[ sub_faces[isub][k] ];
            elB[1+nnsub+k]     = sharedB[ sub_faces[isub][k] ];
          }
          db( ELEMENT, max_element, elB, ddum, length, VERSION_NORMAL, PUT );
          length = 1;
          if ( method_generate==-ELEMENT_GEOMETRY ) {
            // generate an element_geometry record instead of element_group
            element_group = eg_iface;
            db( ELEMENT_GEOMETRY, max_element, &element_group, ddum, length,
              VERSION_NORMAL, PUT );
          }
          else {
            element_group = eg_iface;
            db( ELEMENT_GROUP, max_element, &element_group, ddum, length,
              VERSION_NORMAL, PUT );
          }
          db( ELEMENT_MACRO_GENERATE, max_element, &icontrol, ddum, length,
            VERSION_NORMAL, PUT );
          db( ELEMENT_DOF, max_element, idum, tmp_element_dof, mnolnuknwn,
            VERSION_NORMAL, PUT );
          db( ELEMENT_DOF_INITIALISED, max_element, &zero, ddum, length,
            VERSION_NORMAL, PUT );
          db( NONLOCAL_ELEMENT_INFO, max_element, idum, dworknei, length_nei,
            VERSION_NORMAL, PUT );
        }
        // mark the source elements so the interface is generated only once
        db( ELEMENT_MACRO_GENERATE, iel, &icontrol, ddum, length,
          VERSION_NORMAL, PUT );
        db( ELEMENT_MACRO_GENERATE, jel, &icontrol, ddum, length,
          VERSION_NORMAL, PUT );
      }
    }
  }

  delete[] gen;
  delete[] geometry_entity;
  delete[] elA;
  delete[] elB;
  delete[] nodesA;
  delete[] nodesB;
  delete[] tmp_element_dof;
  delete[] dworknei;

  mesh_has_changed( VERSION_NORMAL );
  if ( swit ) pri( "Out routine GENERATE_INTERFACE." );
}
