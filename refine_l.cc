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

#define MSIDE 25
#define MEL 8
#define NNOL_SIDE 10
#define EPS_NORM 1.e-6

long int refine_get_node_number( long int max_node_refined, long int &tmp_max_node,
  long int version, long int edge_nodes[2], long int ordered_nodes[] );

void refine_locally( void )

  /* 
     First, elements are refined.
     Second, elements which got an extra node on a side
     (from a refined neighbour) will also be refined (=corrected)
     in order to attach correctly to that neighbour.
  */

{
  long int element=0, inod=0, inol=0, correct_this_one=0,
    name=0, max_element=0, max_node=0, nnol=0, swit=0, max_node_refined=0,
    length=0, tmp_node_number_1=0, tmp_node_number_2=0,
    tmp_node_number_extra_1=0, tmp_node_number_extra_2=0,
    tmp_max_element=0, tmp_max_node=-1,
    iside=0, nside=0, inol_side=0, iel=0, nel=0, still_being_refined=1,
    refine_this_one=0, still_being_corrected=0,
    all_not=0, any_only=0, itmp=0, correcting_has_been_done=0,
    do_the_refining=1, do_the_correcting=0,
    nothing=0, control_refine_locally_passes=1,
    ipas=0, icontrol=0, version=0, ldum=0, nnol_side=0,
    control_refine_locally_unknown=0,
    use_control_refine_locally_geometry=0, 
    use_control_refine_locally_not=0, 
    use_control_refine_locally_only=0,
    idum[1], ilist[NNOL_SIDE], control_refine_locally_not[2],
    control_refine_locally_only[2], control_refine_locally_geometry[2],
    tmp_extra_side_nodes[MSIDE], tmp_off_side_nodes[MSIDE][NNOL_SIDE],
    side_nodes[MSIDE][NNOL_SIDE], tmp_side_nodes[MSIDE][NNOL_SIDE], 
    el_nodes[MEL][MNOL], el[MNOL+1], nodes[MNOL], tmp_nodes[MNOL], 
    edge_nodes[2], *ordered_nodes=NULL;
  double tmp=0., h=0., element_residue_norm=0., average_element_residue_norm=0.,
    total=0., larger=0., eps_coord=EPS_COORD, rdum=0., control_refine_locally=0.,
    ddum[MDIM], work[MUKNWN], tmp_node_dof[MUKNWN], tmp_node[MDIM],
    test_coord[MDIM], coord0[MDIM], coord1[MDIM],  
    tmp_node_start_refined[MDIM], tmp_node_dof_start_refined[MUKNWN], 
    test_coord_start_refined[MDIM], *all_element_residue_norm=NULL;

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  if ( db_active_index( CONTROL_MESH_REFINE_LOCALLY, icontrol, VERSION_NORMAL ) ) {

    error( PUT );
    db( CONTROL_MESH_REFINE_LOCALLY, icontrol, idum, &control_refine_locally, 
      ldum, VERSION_NORMAL, GET );
    if ( residue )
      control_refine_locally_unknown = -RESIDUE;
    else
      control_refine_locally_unknown = -NOTHING;
    db( CONTROL_MESH_REFINE_LOCALLY_UNKNOWN, icontrol, 
      &control_refine_locally_unknown, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    use_control_refine_locally_not = db( CONTROL_MESH_REFINE_LOCALLY_NOT, icontrol, 
      control_refine_locally_not, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    use_control_refine_locally_only = db( CONTROL_MESH_REFINE_LOCALLY_ONLY, icontrol, 
      control_refine_locally_only, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    use_control_refine_locally_geometry =  db( CONTROL_MESH_REFINE_LOCALLY_GEOMETRY, 
      icontrol, control_refine_locally_geometry, ddum, ldum, 
      VERSION_NORMAL, GET_IF_EXISTS );

    swit = set_swit(-1,-1,"refine_locally");
    if ( swit ) pri( "In routine REFINE_LOCALLY" );

    db_max_index( NODE, max_node, version, GET ); 
    max_node_refined = 10 * max_node; // heuristic
    ordered_nodes = get_new_int(1+max_node_refined);

    while ( do_the_refining || do_the_correcting ) {
      nel = 0; 
      db_max_index( ELEMENT, max_element, version, GET ); 
      db_max_index( NODE, max_node, version, GET ); 
      tmp_max_element = still_being_refined = still_being_corrected = 0;
      if ( do_the_refining ) {
        ipas++;
        nel = 0;
        average_element_residue_norm = 0.;
        all_element_residue_norm = get_new_dbl( 1+max_element );
        for ( element=0; element<=max_element; element++ ) {
          if ( db_active_index( ELEMENT, element, version ) ) {
            element_residue_norm_set( icontrol, control_refine_locally_unknown,
              element, tmp, version );
            all_element_residue_norm[element] = tmp;
            average_element_residue_norm += tmp;
            nel++;
          }
        }
        if ( nel>0 ) average_element_residue_norm /= nel;
      }
      for ( element=0; element<=max_element; element++ ) {
        if ( db_active_index( ELEMENT, element, version ) ) {
          db( ELEMENT, element, el, ddum, length, version, GET );
          name = el[0];
          nnol = length - 1; array_move( &el[1], nodes, nnol );
          if ( swit ) {
            pri( "element", element );
            pri( "nodes", nodes, nnol );
          }
          all_not = 0; any_only = 0;
          if ( use_control_refine_locally_not ) all_not = 1;
          if ( !use_control_refine_locally_only ) any_only = 1;
          if ( do_the_correcting ) 
            array_move( nodes, tmp_nodes, nnol );
          else {
            assert( do_the_refining );
            for ( inol=0; inol<nnol; inol++ ) {
              inod = nodes[inol];
              if ( use_control_refine_locally_not ) {
                geometry( inod, ddum, control_refine_locally_not, itmp, rdum, 
                  ddum, rdum, ddum, NODE_START_REFINED, PROJECT_EXACT, version );
                if ( !itmp ) all_not = 0;
              }
              if ( use_control_refine_locally_only ) {
                geometry( inod, ddum, control_refine_locally_only, itmp, rdum, 
                  ddum, rdum, ddum, NODE_START_REFINED, PROJECT_EXACT, version );
                if ( itmp ) any_only = 1;
              }
              db( NODE, inod, idum, tmp_node, ldum, version, GET );
              db( NODE_START_REFINED, inod, idum, 
                tmp_node_start_refined, ldum, version, GET );
              if ( nuknwn>0 ) {
                db( NODE_DOF_START_REFINED, inod, idum, 
                  tmp_node_dof_start_refined, ldum, version, GET );
                db( NODE_DOF, inod, idum, tmp_node_dof, ldum, version, GET );
              }
              ordered_list_apply( inod, ordered_nodes, tmp_max_node, tmp_node, eps_coord, 
                tmp_node_number_1, ADD );
              if ( tmp_node_number_1<0 ) {
                tmp_max_node++; 
                if ( tmp_max_node>max_node_refined ) {
                  cout << "max_node_refined should become " << 
                    max_node << " in refine_l.cc\n";
                  exit(TN_EXIT_STATUS);
                }
                tmp_node_number_1 = tmp_max_node;
                create_node( &inod, 1, tmp_node_number_1, tmp_node,
                  tmp_node_dof, tmp_node_start_refined,
                  tmp_node_dof_start_refined,
                  version, VERSION_TMP );
              }
              tmp_nodes[inol] = tmp_node_number_1;
            }
          }
          if ( swit ) pri( "tmp_nodes", tmp_nodes, nnol );
          refine_this_one = 0; correct_this_one = -1;
          if ( name!=-BAR2 && name!=-BAR3 && name!=-TRIA3 && name!=-TRIA6 && 
               name!=-TET4 && name!=-TET10 ) {
            pri( "Error: illegal element for CONTROL_MESH_REFINE_LOCALLY detected." );
            exit(TN_EXIT_STATUS);
          }
          else {
            refine_this_one = 0;
            if ( do_the_refining && !all_not && any_only ) {
              nothing = !element_residue_norm_set( icontrol,
                control_refine_locally_unknown, element, 
                element_residue_norm, version );
              larger = total = 0.;
              db_max_index( ELEMENT, max_element, version, GET ); 
              for ( iel=0; iel<=max_element; iel++ ) {
                if ( db_active_index( ELEMENT, iel, version ) ) {
                  total += 1.;
                  if ( element_residue_norm>= (all_element_residue_norm[iel]+
                    EPS_NORM*average_element_residue_norm) ) larger += 1.;
                }
              }
              refine_this_one = nothing || 
                (larger/total)>=(1.-control_refine_locally/100.);
              if ( refine_this_one ) still_being_refined = 1;
            }
            if ( name==-BAR2 ) {
              nside = 1;
              nnol_side = 2;
            }
            if ( name==-BAR3 ) {
              if ( do_the_refining ) {
                nside = 2;
                nnol_side = 2;
              }
              else {
                assert( do_the_correcting );
                nside = 1;
                nnol_side = 3;
              }
            }
            else if ( name==-TRIA3 ) {
              nside = 3;
              nnol_side = 2;
            }
            else if ( name==-TRIA6 ) {
              if ( do_the_refining ) {
                nside = 9;
                nnol_side = 2;
              }
              else {
                assert( do_the_correcting );
                nside = 3;
                nnol_side = 3;
              }
            }
            else if ( name==-TET4 ) {
              nside = 6;
              nnol_side = 2;
            }
            else {
              assert( name==-TET10 );
              if ( do_the_refining ) {
                nside = 25;
                nnol_side = 2;
              }
              else {
                assert( do_the_correcting );
                nside = 6;
                nnol_side = 3;
              }
            }
            for ( iside=0; iside<nside && correct_this_one<0; iside++ ) {
              if      ( name==-BAR2 ) {
                ilist[0] = 0; ilist[1] = 1;
              }
              else if ( name==-BAR2 ) {
                if ( do_the_refining ) {
                  if ( iside==0 ) {
                    ilist[0] = 0; ilist[1] = 1;
                  }
                  else {
                    ilist[0] = 1; ilist[1] = 2;
                  }
                }
                else {
                  assert( do_the_correcting );
                  ilist[0] = 0; ilist[1] = 1; ilist[2] = 2;
                }
              }
              else if ( name==-TRIA3 ) {
                if      ( iside==0 ) {
                  ilist[0] = 0; ilist[1] = 1;
                  tmp_off_side_nodes[iside][0] = tmp_nodes[2];
                }
                else if ( iside==1 ) {
                  ilist[0] = 1; ilist[1] = 2;
                  tmp_off_side_nodes[iside][0] = tmp_nodes[0];
                }
                else {
                  assert( iside==2 );
                  ilist[0] = 2; ilist[1] = 0;
                  tmp_off_side_nodes[iside][0] = tmp_nodes[1];
                }
              }
              else if ( name==-TRIA6 ) {
                if ( do_the_refining ) {
                  if      ( iside==0 ) {
                    ilist[0] = 1; ilist[1] = 0;
                  }
                  else if ( iside==1 ) {
                    ilist[0] = 3; ilist[1] = 0;
                  }
                  else if ( iside==2 ) {
                    ilist[0] = 1; ilist[1] = 3;
                  }
                  else if ( iside==3 ) {
                    ilist[0] = 1; ilist[1] = 2;
                  }
                  else if ( iside==4 ) {
                    ilist[0] = 4; ilist[1] = 2;
                  }
                  else if ( iside==5 ) {
                    ilist[0] = 1; ilist[1] = 4;
                  }
                  else if ( iside==6 ) {
                    ilist[0] = 4; ilist[1] = 5;
                  }
                  else if ( iside==7 ) {
                    ilist[0] = 3; ilist[1] = 5;
                  }
                  else if ( iside==8 ) {
                    ilist[0] = 3; ilist[1] = 4;
                  }                        
                }
                else {
                  assert( do_the_correcting );
                  if      ( iside==0 ) {
                    ilist[0] = 0; ilist[1] = 1; ilist[2] = 2;
                    tmp_off_side_nodes[iside][0] = tmp_nodes[3];
                    tmp_off_side_nodes[iside][1] = tmp_nodes[4];
                    tmp_off_side_nodes[iside][2] = tmp_nodes[5];
                  }
                  else if ( iside==1 ) {
                    ilist[0] = 2; ilist[1] = 4; ilist[2] = 5;
                    tmp_off_side_nodes[iside][0] = tmp_nodes[1];
                    tmp_off_side_nodes[iside][1] = tmp_nodes[3];
                    tmp_off_side_nodes[iside][2] = tmp_nodes[0];
                  }
                  else {
                    assert( iside==2 );
                    ilist[0] = 5; ilist[1] = 3; ilist[2] = 0;
                    tmp_off_side_nodes[iside][0] = tmp_nodes[4];
                    tmp_off_side_nodes[iside][1] = tmp_nodes[1];
                    tmp_off_side_nodes[iside][2] = tmp_nodes[2];
                  }
                }
              }
              else if ( name==-TET4 ) {
                if      ( iside==0 ) {
                  ilist[0] = 0; ilist[1] = 1;
                  tmp_off_side_nodes[iside][0] = tmp_nodes[2];
                  tmp_off_side_nodes[iside][1] = tmp_nodes[3];
                }
                else if ( iside==1 ) {
                  ilist[0] = 1; ilist[1] = 2;
                  tmp_off_side_nodes[iside][0] = tmp_nodes[0];
                  tmp_off_side_nodes[iside][1] = tmp_nodes[3];
                }
                else if ( iside==2 ) {
                  ilist[0] = 2; ilist[1] = 0;
                  tmp_off_side_nodes[iside][0] = tmp_nodes[1];
                  tmp_off_side_nodes[iside][1] = tmp_nodes[3];
                }
                else if ( iside==3 ) {
                  ilist[0] = 0; ilist[1] = 3;
                  tmp_off_side_nodes[iside][0] = tmp_nodes[1];
                  tmp_off_side_nodes[iside][1] = tmp_nodes[2];
                }
                else if ( iside==4 ) {
                  ilist[0] = 1; ilist[1] = 3;
                  tmp_off_side_nodes[iside][0] = tmp_nodes[0];
                  tmp_off_side_nodes[iside][1] = tmp_nodes[2];
                }
                else {
                  assert ( iside==5 );
                  ilist[0] = 2; ilist[1] = 3;
                  tmp_off_side_nodes[iside][0] = tmp_nodes[0];
                  tmp_off_side_nodes[iside][1] = tmp_nodes[1];
                }
              }
              else {
                assert( name==-TET10 );
                if ( do_the_refining ) {
                  if      ( iside==0 ) {
                    ilist[0] = 0; ilist[1] = 1;
                  }
                  else if ( iside==1 ) {
                    ilist[0] = 0; ilist[1] = 3;
                  }
                  else if ( iside==2 ) {
                    ilist[0] = 0; ilist[1] = 6;
                  }
                  else if ( iside==3 ) {
                    ilist[0] = 1; ilist[1] = 3;
                  }
                  else if ( iside==4 ) {
                    ilist[0] = 1; ilist[1] = 6;
                  }
                  else if ( iside==5 ) {
                    ilist[0] = 3; ilist[1] = 6;
                  }
                  else if ( iside==6 ) {
                    ilist[0] = 2; ilist[1] = 4;
                  }
                  else if ( iside==7 ) {
                    ilist[0] = 2; ilist[1] = 1;
                  }
                  else if ( iside==8 ) {
                    ilist[0] = 2; ilist[1] = 7;
                  }
                  else if ( iside==9 ) {
                    ilist[0] = 4; ilist[1] = 1;
                  }
                  else if ( iside==10 ) {
                    ilist[0] = 4; ilist[1] = 7;
                  }
                  else if ( iside==11 ) {
                    ilist[0] = 7; ilist[1] = 1;
                  }
                  else if ( iside==12 ) {
                    ilist[0] = 5; ilist[1] = 3;
                  }
                  else if ( iside==13 ) {
                    ilist[0] = 5; ilist[1] = 4;
                  }
                  else if ( iside==14 ) {
                    ilist[0] = 5; ilist[1] = 8;
                  }
                  else if ( iside==15 ) {
                    ilist[0] = 3; ilist[1] = 4;
                  }
                  else if ( iside==16 ) {
                    ilist[0] = 3; ilist[1] = 8;
                  }
                  else if ( iside==17 ) {
                    ilist[0] = 4; ilist[1] = 8;
                  }
                  else if ( iside==18 ) {
                    ilist[0] = 9; ilist[1] = 7;
                  }
                  else if ( iside==19 ) {
                    ilist[0] = 9; ilist[1] = 6;
                  }
                  else if ( iside==20 ) {
                    ilist[0] = 9; ilist[1] = 8;
                  }
                  else if ( iside==21 ) {
                    ilist[0] = 7; ilist[1] = 6;
                  }
                  else if ( iside==22 ) {
                    ilist[0] = 7; ilist[1] = 8;
                  }
                  else if ( iside==23 ) {
                    ilist[0] = 6; ilist[1] = 8;
                  }
                  else if ( iside==24 ) {
                    ilist[0] = 7; ilist[1] = 3;
                  }
                }
                else {
                  assert( do_the_correcting );
                  if      ( iside==0 ) {
                    ilist[0] = 0; ilist[1] = 1; ilist[2] = 2;
                    tmp_off_side_nodes[iside][0] = tmp_nodes[3];
                    tmp_off_side_nodes[iside][1] = tmp_nodes[4];
                    tmp_off_side_nodes[iside][2] = tmp_nodes[6];
                    tmp_off_side_nodes[iside][3] = tmp_nodes[7];
                    tmp_off_side_nodes[iside][4] = tmp_nodes[5];
                    tmp_off_side_nodes[iside][5] = tmp_nodes[8];
                    tmp_off_side_nodes[iside][6] = tmp_nodes[9];
                  }
                  else if ( iside==1 ) {
                    ilist[0] = 0; ilist[1] = 3; ilist[2] = 5;
                    tmp_off_side_nodes[iside][0] = tmp_nodes[1];
                    tmp_off_side_nodes[iside][1] = tmp_nodes[4];
                    tmp_off_side_nodes[iside][2] = tmp_nodes[6];
                    tmp_off_side_nodes[iside][3] = tmp_nodes[8];
                    tmp_off_side_nodes[iside][4] = tmp_nodes[2];
                    tmp_off_side_nodes[iside][5] = tmp_nodes[7];
                    tmp_off_side_nodes[iside][6] = tmp_nodes[9];
                  }
                  else if ( iside==2 ) {
                    ilist[0] = 2; ilist[1] = 4; ilist[2] = 5;
                    tmp_off_side_nodes[iside][0] = tmp_nodes[1];
                    tmp_off_side_nodes[iside][1] = tmp_nodes[3];
                    tmp_off_side_nodes[iside][2] = tmp_nodes[7];
                    tmp_off_side_nodes[iside][3] = tmp_nodes[8];
                    tmp_off_side_nodes[iside][4] = tmp_nodes[0];
                    tmp_off_side_nodes[iside][5] = tmp_nodes[6];
                    tmp_off_side_nodes[iside][6] = tmp_nodes[9];
                  }
                  else if ( iside==3 ) {
                    ilist[0] = 0; ilist[1] = 6; ilist[2] = 9;
                    tmp_off_side_nodes[iside][0] = tmp_nodes[1];
                    tmp_off_side_nodes[iside][1] = tmp_nodes[7];
                    tmp_off_side_nodes[iside][2] = tmp_nodes[3];
                    tmp_off_side_nodes[iside][3] = tmp_nodes[8];
                    tmp_off_side_nodes[iside][4] = tmp_nodes[2];
                    tmp_off_side_nodes[iside][5] = tmp_nodes[4];
                    tmp_off_side_nodes[iside][6] = tmp_nodes[5];
                  }
                  else if ( iside==4 ) {
                    ilist[0] = 2; ilist[1] = 7; ilist[2] = 9;
                    tmp_off_side_nodes[iside][0] = tmp_nodes[1];
                    tmp_off_side_nodes[iside][1] = tmp_nodes[6];
                    tmp_off_side_nodes[iside][2] = tmp_nodes[4];
                    tmp_off_side_nodes[iside][3] = tmp_nodes[8];
                    tmp_off_side_nodes[iside][4] = tmp_nodes[0];
                    tmp_off_side_nodes[iside][5] = tmp_nodes[3];
                    tmp_off_side_nodes[iside][6] = tmp_nodes[5];
                  }
                  else {
                    assert( iside==5 );
                    ilist[0] = 5; ilist[1] = 8; ilist[2] = 9;
                    tmp_off_side_nodes[iside][0] = tmp_nodes[4];
                    tmp_off_side_nodes[iside][1] = tmp_nodes[7];
                    tmp_off_side_nodes[iside][2] = tmp_nodes[3];
                    tmp_off_side_nodes[iside][3] = tmp_nodes[6];
                    tmp_off_side_nodes[iside][4] = tmp_nodes[2];
                    tmp_off_side_nodes[iside][5] = tmp_nodes[1];
                    tmp_off_side_nodes[iside][6] = tmp_nodes[0];
                  }
                }
              }
              for ( inol=0; inol<nnol_side; inol++ ) {
                side_nodes[iside][inol] = nodes[ilist[inol]];
                tmp_side_nodes[iside][inol] = tmp_nodes[ilist[inol]];
              }
              db( NODE, side_nodes[iside][0], idum, coord0, ldum, version, GET );
              db( NODE, side_nodes[iside][1], idum, coord1, ldum, version, GET );
              tmp = array_distance( coord0, coord1, work, ndim );
            }
            for ( iside=0; iside<nside && correct_this_one<0; iside++ ) {
              array_set( test_coord, 0., ndim );
              array_set( test_coord_start_refined, 0., ndim );
              if ( nuknwn>0 ) {
                array_set( tmp_node_dof_start_refined, 0., nuknwn );
                array_set( tmp_node_dof, 0., nuknwn );
              }
              for ( inol_side=0; inol_side<2; inol_side++ ) {
                inod = side_nodes[iside][inol_side];
                h = 1./2.;
                db( NODE, inod, idum, work, ldum, version, GET );
                array_multiply( work, work, h, ndim );
                array_add( test_coord, work, test_coord, ndim );
                db( NODE_START_REFINED, inod, idum, work, ldum, version, GET );
                array_multiply( work, work, h, ndim );
                array_add( test_coord_start_refined, work, 
                  test_coord_start_refined, ndim );
                if ( nuknwn>0 ) {
                  db( NODE_DOF_START_REFINED, inod, idum, work, ldum, version, GET );
                  array_multiply( work, work, h, nuknwn );
                  array_add( tmp_node_dof_start_refined, work, 
                    tmp_node_dof_start_refined, nuknwn );
                  db( NODE_DOF, inod, idum, work, ldum, version, GET );
                  array_multiply( work, work, h, nuknwn );
                  array_add( tmp_node_dof, work, tmp_node_dof, nuknwn );
                }
              }
              if ( do_the_correcting ) {
                array_move( test_coord, tmp_node, ndim );
                array_move( test_coord_start_refined, 
                  tmp_node_start_refined, ndim );
                ordered_list_apply( inod, ordered_nodes, tmp_max_node, tmp_node, eps_coord, 
                  tmp_node_number_1, CHECK );
                if ( tmp_node_number_1>0 ) {
                  correct_this_one = iside;
                  if ( name==-TRIA6 ) {
                    edge_nodes[0] = side_nodes[iside][1];
                    edge_nodes[1] = tmp_off_side_nodes[iside][2];
                    tmp_node_number_extra_1 = 
                      refine_get_node_number( max_node_refined, tmp_max_node,
                      version, edge_nodes, ordered_nodes ); 
                    edge_nodes[0] = side_nodes[iside][1];
                    edge_nodes[1] = side_nodes[iside][2];
                    tmp_node_number_2 = 
                      refine_get_node_number( max_node_refined, tmp_max_node,
                      version, edge_nodes, ordered_nodes ); 
                  }
                  else if ( name==-TET10 ) {
                    edge_nodes[0] = side_nodes[iside][1];
                    edge_nodes[1] = tmp_off_side_nodes[iside][4];
                    tmp_node_number_extra_1 = 
                      refine_get_node_number( max_node_refined, tmp_max_node,
                      version, edge_nodes, ordered_nodes ); 
                    edge_nodes[0] = side_nodes[iside][1];
                    edge_nodes[1] = tmp_off_side_nodes[iside][6];
                    tmp_node_number_extra_2 = 
                      refine_get_node_number( max_node_refined, tmp_max_node,
                      version, edge_nodes, ordered_nodes ); 
                    edge_nodes[0] = side_nodes[iside][1];
                    edge_nodes[1] = side_nodes[iside][2];
                    tmp_node_number_2 = 
                      refine_get_node_number( max_node_refined, tmp_max_node,
                      version, edge_nodes, ordered_nodes ); 
                  }
                }
              }
              else {
                if ( refine_this_one ) {
                  array_move( test_coord, tmp_node, ndim );
                  array_move( test_coord_start_refined, tmp_node_start_refined, ndim );
                  if ( use_control_refine_locally_geometry ) 
                    interpolate_geometry( control_refine_locally_geometry,
                      side_nodes[iside], 2, test_coord, tmp_node, 
                      test_coord_start_refined, tmp_node_start_refined, 
                      PROJECT_EXACT, version );
                  ordered_list_apply( inod, ordered_nodes, tmp_max_node, tmp_node, eps_coord, 
                    tmp_node_number_1, ADD );
                  if ( tmp_node_number_1<0 ) {
                    tmp_max_node++;
                    if ( tmp_max_node>max_node_refined ) {
                      cout << "max_node_refined should become " << 
                        max_node << " in refine_l.cc\n";
                      exit(TN_EXIT_STATUS);
                    }
                    tmp_node_number_1 = tmp_max_node;
                    create_node( side_nodes[iside], 2,
                      tmp_node_number_1, tmp_node, tmp_node_dof,
                      tmp_node_start_refined, tmp_node_dof_start_refined,
                      version, VERSION_TMP );
                  }
                  tmp_extra_side_nodes[iside] = tmp_node_number_1;
                }
              }
            }
            if ( swit ) {
              pri( "refine_this_one", refine_this_one );
              pri( "correct_this_one", correct_this_one );
            }
            if ( correct_this_one>=0 ) {
              still_being_corrected = 1;
              iside = correct_this_one;
              if      ( name==-TRIA3 ) {
                nel = 2;
                el_nodes[0][0] = tmp_side_nodes[iside][0];
                el_nodes[0][1] = tmp_node_number_1;
                el_nodes[0][2] = tmp_off_side_nodes[iside][0];
                el_nodes[1][0] = tmp_node_number_1;
                el_nodes[1][1] = tmp_side_nodes[iside][1];
                el_nodes[1][2] = tmp_off_side_nodes[iside][0];
              }
              else if ( name==-TRIA6 ) {
                nel = 2;
                el_nodes[0][0] = tmp_side_nodes[iside][0];
                el_nodes[0][1] = tmp_node_number_1;
                el_nodes[0][2] = tmp_side_nodes[iside][1];
                el_nodes[0][3] = tmp_off_side_nodes[iside][0];
                el_nodes[0][4] = tmp_node_number_extra_1;
                el_nodes[0][5] = tmp_off_side_nodes[iside][2];
                el_nodes[1][0] = tmp_side_nodes[iside][2];
                el_nodes[1][1] = tmp_node_number_2;
                el_nodes[1][2] = tmp_side_nodes[iside][1];
                el_nodes[1][3] = tmp_off_side_nodes[iside][1];
                el_nodes[1][4] = tmp_node_number_extra_1;
                el_nodes[1][5] = tmp_off_side_nodes[iside][2];
              }
              else if ( name==-TET4 ) {
                nel = 2;
                el_nodes[0][0] = tmp_side_nodes[iside][0];
                el_nodes[0][1] = tmp_node_number_1;
                el_nodes[0][2] = tmp_off_side_nodes[iside][0];
                el_nodes[0][3] = tmp_off_side_nodes[iside][1];
                el_nodes[1][0] = tmp_node_number_1;
                el_nodes[1][1] = tmp_side_nodes[iside][1];
                el_nodes[1][2] = tmp_off_side_nodes[iside][0];
                el_nodes[1][3] = tmp_off_side_nodes[iside][1];
              }
              else {
                assert( name==-TET10 );
                nel = 2;
                el_nodes[0][0] = tmp_side_nodes[iside][0];
                el_nodes[0][1] = tmp_node_number_1;
                el_nodes[0][2] = tmp_side_nodes[iside][1];
                el_nodes[0][3] = tmp_side_nodes[iside][3];
                el_nodes[0][4] = tmp_node_number_extra_1;
                el_nodes[0][5] = tmp_side_nodes[iside][5];
                el_nodes[0][6] = tmp_side_nodes[iside][6];
                el_nodes[0][7] = tmp_node_number_extra_2;
                el_nodes[0][8] = tmp_side_nodes[iside][8];
                el_nodes[0][9] = tmp_side_nodes[iside][9];
                el_nodes[1][0] = tmp_side_nodes[iside][2];
                el_nodes[1][1] = tmp_node_number_2;
                el_nodes[1][2] = tmp_side_nodes[iside][1];
                el_nodes[1][3] = tmp_side_nodes[iside][4];
                el_nodes[1][4] = tmp_node_number_extra_1;
                el_nodes[1][5] = tmp_side_nodes[iside][5];
                el_nodes[1][6] = tmp_side_nodes[iside][7];
                el_nodes[1][7] = tmp_node_number_extra_2;
                el_nodes[1][8] = tmp_side_nodes[iside][8];
                el_nodes[1][9] = tmp_side_nodes[iside][9];
              }
            }
            else if ( name==-BAR2 ) {
              if ( refine_this_one ) {
                nel = 2;
                el_nodes[0][0] = tmp_nodes[0];
                el_nodes[0][1] = tmp_extra_side_nodes[0];
                el_nodes[1][0] = tmp_extra_side_nodes[0];
                el_nodes[1][1] = tmp_nodes[1];
              }
              else {
                nel = 1;
                array_move( tmp_nodes, el_nodes[0], nnol );
              }
            }
            else if ( name==-BAR3 ) {
              if ( refine_this_one ) {
                nel = 2;
                el_nodes[0][0] = tmp_nodes[0];
                el_nodes[0][1] = tmp_extra_side_nodes[0];
                el_nodes[1][0] = tmp_extra_side_nodes[0];
                el_nodes[1][1] = tmp_nodes[1];
              }
              else {
                nel = 1;
                array_move( tmp_nodes, el_nodes[0], nnol );
              }
            }
            else if ( name==-TRIA3 ) {
              if ( refine_this_one ) {
                nel = 4;
                el_nodes[0][0] = tmp_nodes[0];
                el_nodes[0][1] = tmp_extra_side_nodes[2];
                el_nodes[0][2] = tmp_extra_side_nodes[0];
                el_nodes[1][0] = tmp_nodes[1];
                el_nodes[1][1] = tmp_extra_side_nodes[0];
                el_nodes[1][2] = tmp_extra_side_nodes[1];
                el_nodes[2][0] = tmp_nodes[2];
                el_nodes[2][1] = tmp_extra_side_nodes[1];
                el_nodes[2][2] = tmp_extra_side_nodes[2];
                el_nodes[3][0] = tmp_extra_side_nodes[0];
                el_nodes[3][1] = tmp_extra_side_nodes[1];
                el_nodes[3][2] = tmp_extra_side_nodes[2];
              }
              else {
                nel = 1;
                array_move( tmp_nodes, el_nodes[0], nnol );
              }
            }
            else if ( name==-TRIA6 ) {
              if ( refine_this_one ) {
                nel = 4;
                el_nodes[0][0] = tmp_nodes[0];
                el_nodes[0][1] = tmp_extra_side_nodes[0];
                el_nodes[0][2] = tmp_nodes[1];
                el_nodes[0][3] = tmp_extra_side_nodes[1];
                el_nodes[0][4] = tmp_extra_side_nodes[2];
                el_nodes[0][5] = tmp_nodes[3];
                el_nodes[1][0] = tmp_nodes[1];
                el_nodes[1][1] = tmp_extra_side_nodes[3];
                el_nodes[1][2] = tmp_nodes[2];
                el_nodes[1][3] = tmp_extra_side_nodes[5];
                el_nodes[1][4] = tmp_extra_side_nodes[4];
                el_nodes[1][5] = tmp_nodes[4];
                el_nodes[2][0] = tmp_nodes[3];
                el_nodes[2][1] = tmp_extra_side_nodes[8];
                el_nodes[2][2] = tmp_nodes[4];
                el_nodes[2][3] = tmp_extra_side_nodes[2];
                el_nodes[2][4] = tmp_extra_side_nodes[5];
                el_nodes[2][5] = tmp_nodes[1];
                el_nodes[3][0] = tmp_nodes[3];
                el_nodes[3][1] = tmp_extra_side_nodes[8];
                el_nodes[3][2] = tmp_nodes[4];
                el_nodes[3][3] = tmp_extra_side_nodes[7];
                el_nodes[3][4] = tmp_extra_side_nodes[6];
                el_nodes[3][5] = tmp_nodes[5];
              }
              else {
                nel = 1;
                array_move( tmp_nodes, el_nodes[0], nnol );
              }
            }
            else if ( name==-TET4 ) {
              if ( refine_this_one ) {
                nel = 8;
                el_nodes[0][0] = tmp_nodes[0];
                el_nodes[0][1] = tmp_extra_side_nodes[2];
                el_nodes[0][2] = tmp_extra_side_nodes[3];
                el_nodes[0][3] = tmp_extra_side_nodes[0];
                el_nodes[1][0] = tmp_nodes[1];
                el_nodes[1][1] = tmp_extra_side_nodes[0];
                el_nodes[1][2] = tmp_extra_side_nodes[4];
                el_nodes[1][3] = tmp_extra_side_nodes[1];
                el_nodes[2][0] = tmp_nodes[2];
                el_nodes[2][1] = tmp_extra_side_nodes[1];
                el_nodes[2][2] = tmp_extra_side_nodes[5];
                el_nodes[2][3] = tmp_extra_side_nodes[2];
                el_nodes[3][0] = tmp_nodes[3];
                el_nodes[3][1] = tmp_extra_side_nodes[5];
                el_nodes[3][2] = tmp_extra_side_nodes[4];
                el_nodes[3][3] = tmp_extra_side_nodes[3];
                el_nodes[4][0] = tmp_extra_side_nodes[0];
                el_nodes[4][1] = tmp_extra_side_nodes[2];
                el_nodes[4][2] = tmp_extra_side_nodes[3];
                el_nodes[4][3] = tmp_extra_side_nodes[4];
                el_nodes[5][0] = tmp_extra_side_nodes[0];
                el_nodes[5][1] = tmp_extra_side_nodes[1];
                el_nodes[5][2] = tmp_extra_side_nodes[2];
                el_nodes[5][3] = tmp_extra_side_nodes[4];
                el_nodes[6][0] = tmp_extra_side_nodes[1];
                el_nodes[6][1] = tmp_extra_side_nodes[2];
                el_nodes[6][2] = tmp_extra_side_nodes[4];
                el_nodes[6][3] = tmp_extra_side_nodes[5];
                el_nodes[7][0] = tmp_extra_side_nodes[2];
                el_nodes[7][1] = tmp_extra_side_nodes[3];
                el_nodes[7][2] = tmp_extra_side_nodes[4];
                el_nodes[7][3] = tmp_extra_side_nodes[5];
              }
              else {
                nel = 1;
                array_move( tmp_nodes, el_nodes[0], nnol );
              }
            }
            else {
              assert( name==-TET10 );
              if ( refine_this_one ) {
                nel = 8;
                  // first corner element
                el_nodes[0][0] = tmp_nodes[0];
                el_nodes[0][1] = tmp_extra_side_nodes[0];
                el_nodes[0][2] = tmp_nodes[1];
                el_nodes[0][3] = tmp_extra_side_nodes[1];
                el_nodes[0][4] = tmp_extra_side_nodes[3];
                el_nodes[0][5] = tmp_nodes[3];
                el_nodes[0][6] = tmp_extra_side_nodes[2];
                el_nodes[0][7] = tmp_extra_side_nodes[4];
                el_nodes[0][8] = tmp_extra_side_nodes[5];
                el_nodes[0][9] = tmp_nodes[6];
                  // second corner element
                el_nodes[1][0] = tmp_nodes[2];
                el_nodes[1][1] = tmp_extra_side_nodes[6];
                el_nodes[1][2] = tmp_nodes[4];
                el_nodes[1][3] = tmp_extra_side_nodes[7];
                el_nodes[1][4] = tmp_extra_side_nodes[9];
                el_nodes[1][5] = tmp_nodes[1];
                el_nodes[1][6] = tmp_extra_side_nodes[8];
                el_nodes[1][7] = tmp_extra_side_nodes[10];
                el_nodes[1][8] = tmp_extra_side_nodes[11];
                el_nodes[1][9] = tmp_nodes[7];
                  // third corner element
                el_nodes[2][0] = tmp_nodes[5];
                el_nodes[2][1] = tmp_extra_side_nodes[12];
                el_nodes[2][2] = tmp_nodes[3];
                el_nodes[2][3] = tmp_extra_side_nodes[13];
                el_nodes[2][4] = tmp_extra_side_nodes[15];
                el_nodes[2][5] = tmp_nodes[4];
                el_nodes[2][6] = tmp_extra_side_nodes[14];
                el_nodes[2][7] = tmp_extra_side_nodes[16];
                el_nodes[2][8] = tmp_extra_side_nodes[17];
                el_nodes[2][9] = tmp_nodes[8];
                  // fourth corner element
                el_nodes[3][0] = tmp_nodes[9];
                el_nodes[3][1] = tmp_extra_side_nodes[18];
                el_nodes[3][2] = tmp_nodes[7];
                el_nodes[3][3] = tmp_extra_side_nodes[19];
                el_nodes[3][4] = tmp_extra_side_nodes[21];
                el_nodes[3][5] = tmp_nodes[6];
                el_nodes[3][6] = tmp_extra_side_nodes[20];
                el_nodes[3][7] = tmp_extra_side_nodes[22];
                el_nodes[3][8] = tmp_extra_side_nodes[23];
                el_nodes[3][9] = tmp_nodes[8];
                  // first middle element on the left
                el_nodes[4][0] = tmp_nodes[1];
                el_nodes[4][1] = tmp_extra_side_nodes[9];
                el_nodes[4][2] = tmp_nodes[4];
                el_nodes[4][3] = tmp_extra_side_nodes[3];
                el_nodes[4][4] = tmp_extra_side_nodes[15];
                el_nodes[4][5] = tmp_nodes[3];
                el_nodes[4][6] = tmp_extra_side_nodes[11];
                el_nodes[4][7] = tmp_extra_side_nodes[10];
                el_nodes[4][8] = tmp_extra_side_nodes[24];
                el_nodes[4][9] = tmp_nodes[7];
                  // second middle element on the left
                el_nodes[5][0] = tmp_nodes[1];
                el_nodes[5][1] = tmp_extra_side_nodes[3];
                el_nodes[5][2] = tmp_nodes[3];
                el_nodes[5][3] = tmp_extra_side_nodes[11];
                el_nodes[5][4] = tmp_extra_side_nodes[24];
                el_nodes[5][5] = tmp_nodes[7];
                el_nodes[5][6] = tmp_extra_side_nodes[4];
                el_nodes[5][7] = tmp_extra_side_nodes[5];
                el_nodes[5][8] = tmp_extra_side_nodes[21];
                el_nodes[5][9] = tmp_nodes[6];
                  // first middle element on the right
                el_nodes[6][0] = tmp_nodes[6];
                el_nodes[6][1] = tmp_extra_side_nodes[23];
                el_nodes[6][2] = tmp_nodes[8];
                el_nodes[6][3] = tmp_extra_side_nodes[21];
                el_nodes[6][4] = tmp_extra_side_nodes[22];
                el_nodes[6][5] = tmp_nodes[7];
                el_nodes[6][6] = tmp_extra_side_nodes[5];
                el_nodes[6][7] = tmp_extra_side_nodes[16];
                el_nodes[6][8] = tmp_extra_side_nodes[24];
                el_nodes[6][9] = tmp_nodes[3];
                  // second middle element on the right
                el_nodes[7][0] = tmp_nodes[3];
                el_nodes[7][1] = tmp_extra_side_nodes[15];
                el_nodes[7][2] = tmp_nodes[4];
                el_nodes[7][3] = tmp_extra_side_nodes[16];
                el_nodes[7][4] = tmp_extra_side_nodes[17];
                el_nodes[7][5] = tmp_nodes[8];
                el_nodes[7][6] = tmp_extra_side_nodes[24];
                el_nodes[7][7] = tmp_extra_side_nodes[10];
                el_nodes[7][8] = tmp_extra_side_nodes[22];
                el_nodes[7][9] = tmp_nodes[7];
              }
              else {
                nel = 1;
                array_move( tmp_nodes, el_nodes[0], nnol );
              }
            }
            assert( nel<=MEL );
            for ( iel=0; iel<nel; iel++ ) {
              el[0] = name; array_move( el_nodes[iel], &el[1], nnol );
              tmp_max_element++; length = nnol + 1;
              create_element( element, tmp_max_element, el, length,
                version, VERSION_TMP );
              if ( swit ) {
                pri( "tmp_max_element", tmp_max_element );
                pri( "el", el, length );
              }
            }
          }
        }
      }
      if ( do_the_refining ) delete[] all_element_residue_norm;
      db_version_copy( VERSION_TMP, version );
      if ( swit ) {
        pri( "still_being_refined", still_being_refined );
        pri( "still_being_corrected", still_being_corrected );
      }
      if ( nothing || !still_being_refined ||
        ipas>=control_refine_locally_passes ) do_the_refining = 0;
      if ( !do_the_refining ) {
        if ( !correcting_has_been_done ) {
          do_the_correcting = 1;
        }
        else
          do_the_correcting = still_being_corrected;
        correcting_has_been_done++;
      }
      if ( swit ) {
        pri( "do_the_refining", do_the_refining );
        pri( "do_the_correcting", do_the_correcting );
      }
      if ( tmp_max_element>10*tmp_max_node ) {
        pri( "Error detected in refine_locally." );
        pri( "Is the mesh not correct before the local refinement starts?" );
        pri( "Else the mesh size is too small." );
        exit(TN_EXIT_STATUS);
      }
    }
    delete[] ordered_nodes;

    db_version_delete( VERSION_TMP );
    renumbering( version, NO, 1, 1, idum, idum ); 
    mesh_has_changed( version );

    if ( swit ) pri( "Out routine REFINE_LOCALLY" );

  }

}

long int element_residue_norm_set( long int icontrol, 
  long int control_refine_locally_unknown,
  long int element, double &element_residue_norm, long int version )

{
  long int inol=0, nnol=0, inod=0, idim=0, jdim=0, indx=0, length=0, 
    name=0, ldum=0, idum[1], el[MNOL+1], nodes[MNOL];
  double res=0., element_volume=0., ddum[1], result[1], node_dof[MUKNWN], 
    unknown_values[MDIM*MDIM];

  element_residue_norm = 0.;

  calcul_matrix = 0;
  if       ( control_refine_locally_unknown==-MATERI_DAMAGE )
    calcul_scalar_indx = dam_indx;
  else if ( control_refine_locally_unknown==-MATERI_DISPLACEMENT ) {
    calcul_vector = 1;
    calcul_vec_indx = dis_indx;
  }
  else if ( control_refine_locally_unknown==-MATERI_PLASTI_KAPPA ) {
    calcul_scalar_indx = kap_indx;
  }
  else if ( control_refine_locally_unknown==-MATERI_STRAIN_ELASTI ) {
    calcul_matrix = 1;
    calcul_mat_indx = epe_indx;
  }
  else if ( control_refine_locally_unknown==-MATERI_STRAIN_PLASTI ) {
    calcul_matrix = 1;
    calcul_mat_indx = epp_indx;
  }
  else if ( control_refine_locally_unknown==-MATERI_STRAIN_TOTAL ) {
    calcul_matrix = 1;
    calcul_mat_indx = ept_indx;
  }
  else if ( control_refine_locally_unknown==-MATERI_STRESS ) {
    calcul_matrix = 1;
    calcul_mat_indx = stres_indx;
  }
  else if ( control_refine_locally_unknown==-MATERI_VELOCITY ) {
    calcul_vector = 1;
    calcul_vec_indx = vel_indx;
  }
  else if ( control_refine_locally_unknown==-MATERI_VOID_FRACTION ) {
    calcul_scalar_indx = void_indx;
  }
  else if ( control_refine_locally_unknown==-NOTHING )
    return 0;
  else if ( control_refine_locally_unknown==-RESIDUE &&
    residue ) {
    calcul_scalar_indx = res_indx;
  }
  else
    db_error( CONTROL_MESH_REFINE_LOCALLY_UNKNOWN, icontrol );
  calcul_operat = -SIZETOT;

  db( ELEMENT, element, el, ddum, length, version, GET );
  name = el[0];
  nnol = length - 1; array_move( &el[1], nodes, nnol );
  for ( inol=0; inol<nnol; inol++ ) {
    inod = nodes[inol];
    db( NODE_DOF, inod, idum, node_dof, ldum, version, GET );
    if ( calcul_matrix ) {
      for ( idim=0; idim<MDIM; idim++ ) {
        for ( jdim=0; jdim<MDIM; jdim++ ) {
          indx = idim*MDIM + jdim;
          unknown_values[indx] = 
            node_dof[calcul_mat_indx+stress_indx(idim,jdim)*nder];
        }
      }
      calculate_operat( unknown_values, inod, ddum, ddum, result, ldum );
    }
    else if ( calcul_vector ) {
      for ( idim=0; idim<ndim; idim++ ) {
        indx = idim;
        unknown_values[indx] = node_dof[calcul_vec_indx+idim*nder];
      }
      calculate_operat( unknown_values, inod, ddum, ddum, result, ldum );
    }
    else
      result[0] = node_dof[calcul_scalar_indx];
    res += result[0];
  }
  res = res/nnol;

  element_volume_set( name, nodes, version, element_volume );
  element_residue_norm = element_volume * res;
  return 1;

}

long int refine_get_node_number( long int max_node_refined, long int &tmp_max_node,
  long int version, long int edge_nodes[2], long int ordered_nodes[] )

{
  long int inol=0, inod=0, nnol_edge=2, tmp_node_number=0, 
    ldum=0, idum[1];
  double h=0., eps_coord=EPS_COORD, work[MUKNWN], 
    tmp_node[MDIM], tmp_node_start_refined[MDIM], 
    tmp_node_dof[MUKNWN], tmp_node_dof_start_refined[MUKNWN];

  array_set( tmp_node, 0., ndim );
  array_set( tmp_node_start_refined, 0., ndim );
  if ( nuknwn>0 ) {
    array_set( tmp_node_dof_start_refined, 0., nuknwn );
    array_set( tmp_node_dof, 0., nuknwn );
  }
  for ( inol=0; inol<nnol_edge; inol++ ) {
    inod = edge_nodes[inol];
    h = 1./2.;
    db( NODE, inod, idum, work, ldum, version, GET );
    array_multiply( work, work, h, ndim );
    array_add( tmp_node, work, tmp_node, ndim );
    db( NODE_START_REFINED, inod, idum, work, ldum, version, GET );
    array_multiply( work, work, h, ndim );
    array_add( tmp_node_start_refined, work, tmp_node_start_refined, ndim );
    if ( nuknwn>0 ) {
      db( NODE_DOF_START_REFINED, inod, idum, work, ldum, version, GET );
      array_multiply( work, work, h, nuknwn );
      array_add( tmp_node_dof_start_refined, work,
        tmp_node_dof_start_refined, nuknwn );
      db( NODE_DOF, inod, idum, work, ldum, version, GET );
      array_multiply( work, work, h, nuknwn );
      array_add( tmp_node_dof, work, tmp_node_dof, nuknwn );
    }                                       
  }
  ordered_list_apply( inod, ordered_nodes, tmp_max_node, tmp_node, eps_coord,
    tmp_node_number, ADD );
  if ( tmp_node_number<0 ) {
    tmp_max_node++;
    assert( tmp_max_node<=max_node_refined );
    tmp_node_number = tmp_max_node;
    create_node( edge_nodes, nnol_edge, tmp_node_number, tmp_node, tmp_node_dof,
      tmp_node_start_refined, tmp_node_dof_start_refined, version, VERSION_TMP );
  }
  return tmp_node_number;


}
