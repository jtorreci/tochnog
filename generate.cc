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
