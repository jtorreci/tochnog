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

#define MSIDE 6
#define MUL 20
#define EPS_DISTANCE 1.e-8

void refine_globally( long int control_refine_globally[4], 
  long int length_control_refine_globally, 
  long int use_control_refine_globally_geometry, 
  long int control_refine_globally_geometry[2],
  long int project_type, long int version )

{
  long int n=0, element=0, swit_element=0, inod=0, inol=0,
    name=0, max_element=0, max_node=0, nnol=0, swit=0, 
    length=0, tmp_node_number=0, tmp_max_element=0, tmp_max_node=-1,
    iside=0, name_new=0, npol=0, npol_new=0, iel_xi_new=0,
    iel_eta_new=0, iel_zeta_new=0, nel_xi_new=1, nel_eta_new=1,
    nel_zeta_new=1, inol_xi_new=0, inol_eta_new=0, inol_zeta_new=0,
    inol_new=0, nnol_xi_new=1, nnol_eta_new=1, nnol_zeta_new=1, nnol_new=0,
    inol_xi=0, inol_eta=0, inol_zeta=0, nnol_xi=1, nnol_eta=1, nnol_zeta=1,
    first_xi=0, last_xi=0, first_eta=0, last_eta=0, first_zeta=0, last_zeta=0,
    icontrol=0, indx=0, part_of_all_sides=0,
    ldum=0, idum[1], sides_member[MSIDE], sides[MSIDE][MNOL], el[MNOL+1], 
    nodes[MNOL], tmp_nodes[MNOL], node_numbers[MNOL], sides_counter[MSIDE], 
    *ordered_nodes=NULL;
  double xi=0., eta=0., zeta=0., h=0., eps_coord=EPS_COORD,
    tmp_node_dof[MUKNWN], tmp_node[MDIM], 
    ddum[MPOINT*MDIM*MNOL], h_xi[MNOL], h_eta[MNOL], h_zeta[MNOL], 
    iso_new[MPOINT], work[MUKNWN], test_coord[MDIM], 
    tmp_node_start_refined[MDIM], test_coord_start_refined[MDIM],
    tmp_node_dof_start_refined[MUKNWN];

  swit = set_swit(-1,-1,"refine_globally");
  if ( swit ) pri( "In routine MESH_REFINE_GLOBALLY" );

  if ( control_refine_globally[0]==-H_REFINEMENT ) {
    if ( length_control_refine_globally==1+ndim ) {
      if ( ndim>=1 && control_refine_globally[1]==-YES ) nel_xi_new = 2;
      if ( ndim>=2 && control_refine_globally[2]==-YES ) nel_eta_new = 2;
      if ( ndim>=3 && control_refine_globally[3]==-YES ) nel_zeta_new = 2;
    }
    else {
      if ( ndim>=1 ) nel_xi_new   = 2;
      if ( ndim>=2 ) nel_eta_new  = 2;
      if ( ndim==3 ) nel_zeta_new = 2;
    }
  }

  db_max_index( ELEMENT, max_element, version, GET ); 
  db_max_index( NODE, max_node, version, GET ); 

  ordered_nodes = get_new_int(MUL*(1+max_node));

  for ( element=0; element<=max_element; element++ ) {
    if ( db_active_index( ELEMENT, element, version ) ) {
      swit_element = swit; swit = swit && set_swit( element, -1, "" );
      db( ELEMENT, element, el, ddum, length, version, GET );
      name = el[0];
      nnol = length - 1; array_move( &el[1], nodes, nnol );
      if ( swit ) {
        pri( "element", element );
        pri( "nodes", nodes, nnol );
      }
      if ( name==-TRIA3 || name==-TET4 ) {
        pri( "Error: -TRIA3 and -TET4 not available for CONTROL_MESH_REFINE_GLOBALLY.");
        exit(TN_EXIT_STATUS);
      }
      else {
        if      ( name==-BAR2 || name==-QUAD4  || name==-HEX8   ) npol = 2;
        else if ( name==-BAR3 || name==-QUAD9  || name==-HEX27  ) npol = 3;
        else if ( name==-BAR4 || name==-QUAD16 || name==-HEX64  ) npol = 4;
        else db_error( ELEMENT, element );
        if      ( control_refine_globally[0]==-P_REFINEMENT ) {
          if      ( name==-BAR2 )   name_new = -BAR3;
          else if ( name==-BAR3 )   name_new = -BAR4;
          else if ( name==-QUAD4 )  name_new = -QUAD9;
          else if ( name==-QUAD9 )  name_new = -QUAD16;
          else if ( name==-HEX8 )   name_new = -HEX27;
          else if ( name==-HEX27 )  name_new = -HEX64;
          npol_new = npol + 1;
        }
        else if ( control_refine_globally[0]==-P_COARSEN ) {
          if      ( name==-BAR2 )   name_new = -BAR2;
          else if ( name==-BAR3 )   name_new = -BAR2;
          else if ( name==-BAR4 )   name_new = -BAR3;
          else if ( name==-QUAD4 )  name_new = -QUAD4;
          else if ( name==-QUAD9 )  name_new = -QUAD4;
          else if ( name==-QUAD16 ) name_new = -QUAD9;
          else if ( name==-HEX8 )   name_new = -HEX8;
          else if ( name==-HEX27 )  name_new = -HEX8;
          else if ( name==-HEX64 )  name_new = -HEX27;
          if ( name_new==name ) npol_new = npol;
          else npol_new = npol - 1;
        }
        else if ( control_refine_globally[0]==-H_REFINEMENT ) {
          name_new = name;
          npol_new = npol;
        }
        else
          db_error( CONTROL_MESH_REFINE_GLOBALLY, icontrol );
        if ( npol_new>=1 && npol_new<=5 ) {
          if ( ndim>=1 ) {
            nnol_xi = npol;
            nnol_xi_new = npol_new;
          }
          if ( ndim>=2 ) {
            nnol_eta = npol;
            nnol_eta_new = npol_new;
          }
          if ( ndim>=3 ) {
            nnol_zeta = npol;
            nnol_zeta_new = npol_new;
          }
          nnol_new = nnol_xi_new*nnol_eta_new*nnol_zeta_new;
          assert( nnol_new<=MNOL );
          integration_lobatto( npol_new, iso_new, ddum );
          if ( swit ) {
            pri( "npol", npol );
            pri( "nnol_xi", nnol_xi );
            pri( "nnol_eta", nnol_eta );
            pri( "nnol_zeta", nnol_zeta );
            pri( "nel_xi_new", nel_xi_new );
            pri( "nel_eta_new", nel_eta_new );
            pri( "nel_zeta_new", nel_zeta_new );
            pri( "npol_new", npol_new );
            pri( "nnol_xi_new", nnol_xi_new );
            pri( "nnol_eta_new", nnol_eta_new );
            pri( "nnol_zeta_new", nnol_zeta_new );
            pri( "nnol_new", nnol_new );
            pri( "iso_new", iso_new, npol_new );
          }
          inol = 0; array_set( sides_counter, 0, MSIDE ); 
          array_set( &sides[0][0], 0, MSIDE*MNOL );
          for ( inol_zeta=0; inol_zeta<nnol_zeta; inol_zeta++ ) {
            for ( inol_eta=0; inol_eta<nnol_eta; inol_eta++ ) {
              for ( inol_xi=0; inol_xi<nnol_xi; inol_xi++ ) {
                if ( ndim>=1 && inol_xi==0 ) {
                  iside = 0;
                  sides[iside][sides_counter[iside]] = nodes[inol];
                  sides_counter[iside]++;
                }
                if ( ndim>=1 && inol_xi==(nnol_xi-1) ) {
                  iside = 1;
                  sides[iside][sides_counter[iside]] = nodes[inol];
                  sides_counter[iside]++;
                }
                if ( ndim>=2 && inol_eta==0 ) {
                  iside = 2;
                  sides[iside][sides_counter[iside]] = nodes[inol];
                  sides_counter[iside]++;
                }
                else if ( ndim>=2 && inol_eta==(nnol_eta-1) ) {
                  iside = 3;
                  sides[iside][sides_counter[iside]] = nodes[inol];
                  sides_counter[iside]++;
                }
                if ( ndim==3 && inol_zeta==0 ) {
                  iside = 4;
                  sides[iside][sides_counter[iside]] = nodes[inol];
                  sides_counter[iside]++;
                }
                if ( ndim==3 && inol_zeta==(nnol_zeta-1) ) {
                  iside = 5;
                  sides[iside][sides_counter[iside]] = nodes[inol];
                  sides_counter[iside]++;
                }
                inol++;
              }
            }
          }
          if ( swit ) pri( "sides_counter", sides_counter, MSIDE );
          for ( iel_zeta_new=0; iel_zeta_new<nel_zeta_new; iel_zeta_new++ ) {
            for ( iel_eta_new=0; iel_eta_new<nel_eta_new; iel_eta_new++ ) {
              for ( iel_xi_new=0; iel_xi_new<nel_xi_new; iel_xi_new++ ) {
                inol_new = 0;
                for (inol_zeta_new=0;inol_zeta_new<nnol_zeta_new;
                     inol_zeta_new++){
                  if ( nel_zeta_new==1 )
                    zeta = iso_new[inol_zeta_new];
                  else {
                    if ( iel_zeta_new==0 )
                      zeta = -0.5 + iso_new[inol_zeta_new]/2.;
                    else
                      zeta = +0.5 + iso_new[inol_zeta_new]/2.;
                  }
                  interpolation_polynomial( zeta, nnol_zeta, h_zeta, ddum );
                  for ( inol_eta_new=0; inol_eta_new<nnol_eta_new;
                    inol_eta_new++){
                    if ( nel_eta_new==1 )
                      eta = iso_new[inol_eta_new];
                    else {
                      if ( iel_eta_new==0 )
                        eta = -0.5 + iso_new[inol_eta_new]/2.;
                      else
                        eta = +0.5 + iso_new[inol_eta_new]/2.;
                    }
                    interpolation_polynomial( eta, nnol_eta, h_eta, ddum );
                    for ( inol_xi_new=0; inol_xi_new<nnol_xi_new;
                          inol_xi_new++ ) {
                      if ( nel_xi_new==1 )
                        xi = iso_new[inol_xi_new];
                      else {
                        if ( iel_xi_new==0 )
                          xi = -0.5 + iso_new[inol_xi_new]/2.;
                        else
                          xi = +0.5 + iso_new[inol_xi_new]/2.;
                      }
                      interpolation_polynomial( xi, nnol_xi, h_xi, ddum );
                      inol = 0; 
                      array_set( test_coord, 0., ndim );
                      array_set( test_coord_start_refined, 0., ndim );
                      if ( nuknwn>0 ) {
                        array_set( tmp_node_dof_start_refined, 0., nuknwn );
                        array_set( tmp_node_dof, 0., nuknwn );
                      }
                      for ( inol_zeta=0; inol_zeta<nnol_zeta; inol_zeta++ ) {
                        for ( inol_eta=0; inol_eta<nnol_eta; inol_eta++ ) {
                          for ( inol_xi=0; inol_xi<nnol_xi; inol_xi++ ) {
                            inod = nodes[inol];
                            h = h_xi[inol_xi] * h_eta[inol_eta] * h_zeta[inol_zeta];
                            db( NODE, inod, idum, work, ldum, version, GET );
                            array_multiply( work, work, h, ndim );
                            array_add( test_coord, work, test_coord, ndim );
                            db( NODE_START_REFINED, inod, idum, work, 
                              ldum, version, GET );
                            array_multiply( work, work, h, ndim );
                            array_add( test_coord_start_refined, work, 
                              test_coord_start_refined, ndim );
                            if ( nuknwn>0 ) {
                              db( NODE_DOF_START_REFINED, inod, idum, work, 
                                ldum, version, GET );
                              array_multiply( work, work, h, nuknwn );
                              array_add( tmp_node_dof_start_refined, work, 
                                tmp_node_dof_start_refined, nuknwn);
                              db( NODE_DOF, inod, idum, work, ldum, version, GET );
                              array_multiply( work, work, h, nuknwn );
                              array_add( tmp_node_dof, work, 
                                tmp_node_dof, nuknwn);
                            }
                            inol++;
                          }
                        }
                      }
                      first_xi = ( iel_xi_new==0 && inol_xi_new==0 );
                      last_xi = ( iel_xi_new==(nel_xi_new-1) &&
                        inol_xi_new==(nnol_xi_new-1) );
                      first_eta = ( iel_eta_new==0 && inol_eta_new==0 );
                      last_eta = ( iel_eta_new==(nel_eta_new-1) &&
                        inol_eta_new==(nnol_eta_new-1) );
                      first_zeta = ( iel_zeta_new==0 && inol_zeta_new==0 );
                      last_zeta = ( iel_zeta_new==(nel_zeta_new-1) &&
                        inol_zeta_new==(nnol_zeta_new-1) );
                      array_set( sides_member, 0, MSIDE );
                      if      ( ndim==1 ) {
                        if ( first_xi )   sides_member[0] = 1;
                        if ( last_xi )    sides_member[1] = 1;
                      }
                      else if ( ndim==2 ) {
                        if ( first_xi )   sides_member[0] = 1;
                        if ( last_xi )    sides_member[1] = 1;
                        if ( first_eta )  sides_member[2] = 1;
                        if ( last_eta )   sides_member[3] = 1;
                      }
                      else {
                        assert( ndim==3 );
                        if ( first_xi )   sides_member[0] = 1;
                        if ( last_xi )    sides_member[1] = 1;
                        if ( first_eta )  sides_member[2] = 1;
                        if ( last_eta )   sides_member[3] = 1;
                        if ( first_zeta ) sides_member[4] = 1;
                        if ( last_zeta )  sides_member[5] = 1;
                      }
                      n = 0;
                      for ( inol=0; inol<nnol; inol++ ) {
                        inod = nodes[inol];
                        if ( !array_member(node_numbers,inod,n,indx) ) {
                          part_of_all_sides = 1;
                          for ( iside=0; iside<MSIDE; iside++ ) {
                            if ( sides_member[iside] ) {
                              if ( !array_member(&sides[iside][0],inod,
                                sides_counter[iside],indx) ) part_of_all_sides = 0;
                            }
                          }
                          if ( part_of_all_sides ) {
                            node_numbers[n] = inod;
                            n++;
                          }
                        }
                      }
                      if ( swit ) pri( "node_numbers", node_numbers, n );
                      if ( swit ) {
                        pri( "test_coord", test_coord, ndim );
                        pri( "test_coord_start_refined", 
                          test_coord_start_refined, ndim );
                      }
                      array_move( test_coord, tmp_node, ndim );
                      array_move( test_coord_start_refined, 
                          tmp_node_start_refined, ndim );
                      if ( use_control_refine_globally_geometry ) {
                        interpolate_geometry( control_refine_globally_geometry,
                          node_numbers, n, test_coord, 
                          tmp_node, test_coord_start_refined, 
                          tmp_node_start_refined, project_type, version );
                      }
                      ordered_list_apply( inod, ordered_nodes, tmp_max_node, tmp_node, 
                        eps_coord, tmp_node_number, ADD );
                      if ( tmp_node_number<0 ) {
                        tmp_max_node++; assert( tmp_max_node<MUL*max_node );
                        tmp_node_number = tmp_max_node;
                        if ( n>=0 ) create_node( node_numbers, n, 
                          tmp_node_number, tmp_node, tmp_node_dof,
                          tmp_node_start_refined, tmp_node_dof_start_refined,
                          version, VERSION_TMP );
                        else create_node( nodes, nnol, tmp_node_number, 
                          tmp_node, tmp_node_dof, tmp_node_start_refined,
                          tmp_node_dof_start_refined,
                          version, VERSION_TMP );
                        if ( swit ) {
                          pri( "tmp_node_number", tmp_node_number );
                          pri( "coord", tmp_node, ndim );
                        }
                      }
                      tmp_nodes[inol_new] = tmp_node_number; inol_new++;
                    }
                  }
                }
                tmp_max_element++;
                el[0] = name_new; array_move( tmp_nodes, &el[1], nnol_new );
                length = nnol_new + 1;
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
      }
      swit = swit_element;
    }
  }

  delete[] ordered_nodes;

  db_version_copy( VERSION_TMP, version );
  db_version_delete( VERSION_TMP );

  renumbering( version, NO, 1, 1, idum, idum );
  mesh_has_changed( version );

  if ( swit ) pri( "Out routine REFINE_GLOBALLY" );

}
