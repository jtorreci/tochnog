/*
    copyright (c) 1998  dennis roddeman
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

long int print_gid_5( long int task )

{
  long int inol=0, inod=0, element=0, max_node=0, max_element=0, name=0,
    element_group=0, length=0, idim=0, jdim=0, ipuknwn=0, iuknwn=0,
    icalcul=0, nval=0, ready=0, indx=0, eigen_analysis=0, neigen=0,
    isoparametric_name=0, length_post_calcul_scal_vec_mat=0,
    load_type=0, data_type=0, data_loc=0, desc_comp=0,
    n=0, swit=0, ldum=0, nnol=0, type=0, ieigen=0,
    mesh_elements=0, bound_elements=0, bound_nodes=0, 
    intersection=0, itmp=0, itendon=0, ntendon=0, 
    nmesh=0, nbound=0, first_time_gid=0,
    print_deformation=0, element_empty=-NO, icontrol=0, 
    control_print_gid_mesh=0, idum[1], 
    *dof_amount=NULL, *dof_label=NULL, *dof_type=NULL,
    *post_calcul_scal_vec_mat=NULL, *dof_scal_vec_mat=NULL, 
    *nodes=NULL, *el=NULL;
  double time_current=0., control_print_gid_time=0.,
     step_val=0., tmp=0., ddum[1], coord[MDIM],
     coord_start_refined[MDIM], 
     *element_tendon_intersections=NULL, *control_eigen_values=NULL, 
     *node_eigen=NULL, *node_dof=NULL, 
     *node_dof_calcul=NULL;
  char str[MCHAR], descr_menu[MCHAR], filename[MCHAR], 
     filename_with_number[MCHAR], outname[MCHAR];

  if ( ndim==1 ) return 1;

  swit = set_swit(-1,-1,"print_gid_5");
  if ( swit ) pri( "In routine PRINT_GID_5" );

  if ( task!=-YES && task!=-SEPARATE ) {
    pri( "Error detected in gid 5 printing." );
    pri( "Gid 5 printing doesn't understand ", task );
    exit_tn_on_error();
  }

  dof_amount = get_new_int(MUKNWN);
  dof_label = get_new_int(MUKNWN);
  dof_type = get_new_int(MUKNWN);
  dof_scal_vec_mat = get_new_int(MUKNWN);
  post_calcul_scal_vec_mat = get_new_int(DATA_ITEM_SIZE);
  nodes = get_new_int(MAXIMUM_NODE);
  el = get_new_int(MAXIMUM_NODE+1);
  element_tendon_intersections = get_new_dbl(DATA_ITEM_SIZE);
  control_eigen_values = get_new_dbl(DATA_ITEM_SIZE);

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( TIME_CURRENT, 0, idum, &time_current, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  if ( db( CONTROL_PRINT_GID_TIME, 0, idum, &control_print_gid_time, 
       ldum, VERSION_NORMAL, GET_IF_EXISTS ) ) {
    if ( time_current<=control_print_gid_time && task!=-SEPARATE ) return 1;
  }
  else {
    first_time_gid = 1;
  }
  length = 1;
  db( CONTROL_PRINT_GID_TIME, 0, idum, &time_current, length, VERSION_NORMAL, PUT );

  db_copy( NODE_RHSIDE, NODE_RHSIDE_PRINT, VERSION_NORMAL );
  db_version_copy( VERSION_NORMAL, VERSION_PRINT );
  renumbering( VERSION_PRINT, NO, 1, 1, idum, idum );

  db_highest_index( NODE, max_node, VERSION_PRINT );
  if ( max_node<0 ) return 1;
  db_highest_index( ELEMENT, max_element, VERSION_PRINT );
  if ( max_element<0 ) return 1;

  if ( db( CONTROL_EIGEN_VALUES, 0, idum, control_eigen_values,
    neigen, VERSION_NORMAL, GET_IF_EXISTS ) ) eigen_analysis = 1;
  else {
    eigen_analysis = 0;
    neigen = 1;
  }

  strcpy( filename_with_number, data_file_base );   
  if ( task==-SEPARATE ) {
    db( CONTROL_PRINT_GID_MESH, 0, &control_print_gid_mesh, ddum, 
      length, VERSION_NORMAL, GET_IF_EXISTS );
    if      ( control_print_gid_mesh<10 )
      strcat( filename_with_number, "00" );
    else if ( control_print_gid_mesh<100 )
      strcat( filename_with_number, "0" );
    long_to_a( control_print_gid_mesh, str );
    strcat( filename_with_number, str );
    control_print_gid_mesh++; length=1;
    db( CONTROL_PRINT_GID_MESH, 0, &control_print_gid_mesh, ddum, 
      length, VERSION_NORMAL, PUT );
  }                    

  strcpy( filename, filename_with_number );
  if ( ndim==2 ) 
    strcat( filename, "5.flavia.dat" );
  else 
    strcat( filename, "5.flavia.msh" );
  if ( first_time_gid || task==-SEPARATE ) {
    ofstream outmesh( filename );
    outmesh.close();
  }
  ofstream outmesh( filename, ios::app );
  outmesh.precision(TN_PRECISION);

  strcpy( filename, filename_with_number );
  strcat( filename, "5.flavia.res" );
  if ( first_time_gid || task==-SEPARATE ) {
    ofstream outres( filename );
    outres.close();
  }
  ofstream outres( filename, ios::app );
  outres.precision(TN_PRECISION);
  outres.setf(ios::showpoint);

    // determine the type of the isoparametric elements in the mesh
  for ( element=1; element<=max_element; element++ ) {
    if ( db_active_index(  ELEMENT, element, VERSION_PRINT ) ) {
      db( ELEMENT, element, el, ddum, length, VERSION_PRINT, GET );
      name = el[0];
      if      ( name==-TRIA3 ) itmp = 3;
      else if ( name==-TRIA6 ) itmp = 6;
      else if ( name==-QUAD4 ) itmp = 4;
      else if ( name==-QUAD9 ) itmp = 9;
      else if ( name==-TET4  ) itmp = 3;
      else if ( name==-TET10 ) itmp = 3;
      else if ( name==-HEX8  ) itmp = 1;
      else if ( name==-HEX27 ) itmp = 1;
      else itmp = -1;
      if ( itmp>0 ) {
        type = itmp;
        isoparametric_name = name;
      }
    }
  }

    // determine the number of elements
  mesh_elements = bound_elements = 0;
  bound_nodes = max_node;
  for ( element=1; element<=max_element; element++ ) {
    if ( db_active_index(  ELEMENT, element, VERSION_PRINT ) ) {
      db( ELEMENT, element, el, ddum, length, VERSION_PRINT, GET );
      name = el[0];
      element_empty=-NO;
      db( ELEMENT_EMPTY, element, &element_empty, ddum,
        ldum, VERSION_PRINT, GET_IF_EXISTS );
      if ( element_empty==-NO || element_empty==-FRONT ) {       
        if ( name==isoparametric_name ) mesh_elements++;
        if ( name==-TRUSS || name==-BEAM ) bound_elements++;
        if ( db_active_index( ELEMENT_TENDON_NUMBER, element, VERSION_PRINT ) ) {
          ntendon = db_len( ELEMENT_TENDON_NUMBER, element, VERSION_PRINT );
          bound_elements += ntendon;
          bound_nodes += ntendon * 2;
        }
      }
    }
  }

  if ( bound_elements )
    strcpy( filename, "lines.bon" );
  else
    strcpy( filename, "tn.tmp" );
  if ( first_time_gid ) {
    ofstream outbound( filename );
    outbound.close();
  }
  ofstream outbound( filename, ios::app );

    // free initial lines for mesh file
  outmesh << " \n";
  outmesh << " \n";
  outmesh << " \n";
  outmesh << " \n";
  outmesh << " \n";
  outmesh << "n_mesh_element n_mesh_points n_element_type:\n";
  outmesh << mesh_elements << " " << max_node << " " << type << "\n";

    // free initial lines for boundary file
  if ( bound_elements ) {
    outbound << " \n";
    outbound << " \n";
    outbound << " \n";
    outbound << " \n";
    outbound << " \n";
    outbound << "n_bound_element n_bound_points n_element_type:\n";
    outbound << bound_elements << " " << bound_nodes << " 11" << "\n";
  }

  outmesh << "Coordinates:\n";
  if ( bound_elements ) outbound << "Coordinates:\n";
  for ( inod=1; inod<=max_node; inod++ ) {
    if ( db_active_index( NODE_START_REFINED, inod, VERSION_PRINT ) ) {
      outmesh << inod << " ";
      if ( bound_elements ) outbound << inod << " ";
      db( NODE_START_REFINED, inod, idum, coord_start_refined, ldum, VERSION_PRINT, GET );
      for ( idim=0; idim<ndim; idim++ ) {
        outmesh << coord_start_refined[idim] << " ";
        if ( bound_elements ) outbound << coord_start_refined[idim] << " ";
      }
      outmesh << "\n";
      if ( bound_elements ) outbound << "\n";
    }
  }
    // tendon coordinates
  for ( element=1; element<=max_element; element++ ) {
    if ( db_active_index(  ELEMENT, element, VERSION_PRINT ) ) {
      db( ELEMENT, element, el, ddum, length, VERSION_PRINT, GET );
      if ( db_active_index( ELEMENT_TENDON_NUMBER, element, VERSION_PRINT ) ) {
        ntendon = db_len( ELEMENT_TENDON_NUMBER, element, VERSION_PRINT );
        ldum =ntendon * 2 * ndim;
        db( ELEMENT_TENDON_INTERSECTIONS, element, idum,
          element_tendon_intersections, ldum, VERSION_PRINT, GET_AND_CHECK );
        for ( itendon=0; itendon<ntendon; itendon++ ) {
          for ( intersection=0; intersection<2; intersection++ ) {
            outbound << inod << " ";
            for ( idim=0; idim<ndim; idim++ ) {
              outbound << 
                element_tendon_intersections[itendon*2*ndim+intersection*ndim+idim] << " ";
            }
            outbound << "\n";
            inod++;
          }
        }
      }
    }
  }

    // connectivity
  nmesh = nbound = 0;
  outmesh << "Connectivities:\n";
  if ( bound_elements ) outbound << "Connectivities:\n";
  for ( element=1; element<=max_element; element++ ) {
    if ( db_active_index( ELEMENT, element, VERSION_PRINT ) ) {
      if ( db_active_index( ELEMENT_GROUP, element, VERSION_PRINT ) )
        db( ELEMENT_GROUP, element, &element_group, ddum, ldum, VERSION_PRINT, GET );
      else
        element_group = 0;
      element_empty = -NO;
      db( ELEMENT_EMPTY, element, &element_empty, ddum, 
        ldum, VERSION_PRINT, GET_IF_EXISTS );
      if ( element_empty==-NO || element_empty==-FRONT ) {
        db( ELEMENT, element, el, ddum, length, VERSION_PRINT, GET );
        name = el[0];
        nnol = length - 1; array_move( &el[1], nodes, nnol );
        if      ( name==-TRIA3 ) {
          nmesh++;
          outmesh << nmesh << " ";
          outmesh << nodes[0] << " " << nodes[1] << " " << nodes[2] << " ";
          outmesh << element_group << "\n";
        }
        else if ( name==-TRIA6 ) {
          nmesh++;
          outmesh << nmesh << " ";
          outmesh << nodes[0] << " " << nodes[1] << " " << nodes[2] << " " <<
                     nodes[4] << " " << nodes[5] << " " << nodes[3] << " ";
          outmesh << element_group << "\n";
        }
        else if ( name==-QUAD4 ) {
          nmesh++;
          outmesh << nmesh << " ";
          outmesh << nodes[0] << " " << nodes[1] << " " << nodes[3] << " " << nodes[2] << " ";
          outmesh << element_group << "\n";
        }
        else if ( name==-QUAD9 ) {
          nmesh++;
          outmesh << nmesh << " ";
          outmesh << nodes[0] << " " << nodes[1] << " " << nodes[2] << " " << 
                     nodes[5] << " " << nodes[8] << " " << nodes[7] << " " << 
                     nodes[6] << " " << nodes[3] << " " << nodes[4] << " ";
          outmesh << element_group << "\n";
        }
        else if ( name==-TET4 ) {
          nmesh++;
          outmesh << nmesh << " ";
          outmesh << nodes[0] << " " << nodes[1] << " " << nodes[2] << " " << nodes[3] << " ";
          outmesh << element_group << "\n";
        }
        else if ( name==-TET10 ) {
          nmesh++;
          outmesh << nmesh << " ";
          outmesh << nodes[0] << " " << nodes[2] << " " << nodes[5] << " " << nodes[9] << " ";
          outmesh << element_group << "\n";
        }
        else if ( name==-HEX8 ) {
          nmesh++;
          outmesh << nmesh << " ";
          outmesh << nodes[0] << " " << nodes[1] << " " << nodes[3] << " " << nodes[2] << " " <<
                     nodes[4] << " " << nodes[5] << " " << nodes[7] << " " << nodes[6] << " ";
          outmesh << element_group << "\n";
        }
        else if ( name==-HEX27 ) {
          nmesh++;
          outmesh << nmesh << " ";
          outmesh << nodes[0] << " " << nodes[2] << " " << nodes[8] << " " << nodes[6] << " " <<
                     nodes[18] << " " << nodes[20] << " " << nodes[26] << " " << nodes[24] << " ";
          outmesh << element_group << "\n";
        }
        else if ( name==-TRUSS ) {
          nbound++;
          outbound << nbound << " ";
          outbound << nodes[0] << " " << nodes[1] << " ";
          outbound << element_group << "\n";
        }
        else if ( name==-BEAM ) {
          nbound++;
          outbound << nbound << " ";
          outbound << nodes[0] << " " << nodes[1] << " ";
          outbound << element_group << "\n";
        }
      }
    }
  }
    // tendon connectivity
  inod = max_node;
  for ( element=1; element<=max_element; element++ ) {
    if ( db_active_index(  ELEMENT, element, VERSION_PRINT ) ) {
      db( ELEMENT, element, el, ddum, length, VERSION_PRINT, GET );
      if ( db_active_index( ELEMENT_TENDON_NUMBER, element, VERSION_PRINT ) ) {
        ntendon = db_len( ELEMENT_TENDON_NUMBER, element, VERSION_PRINT );
        for ( itendon=0; itendon<ntendon; itendon++ ) {
          nbound++;
          outbound << nbound << " ";
          for ( inol=0; inol<2; inol++ ) {
            inod++;
            outbound << inod << " ";
          }
          outbound << "100" << "\n";
        }
      }
    }
  }

    // results
  load_type = 1;
  data_loc = 1;
  for ( ieigen=0; ieigen<neigen; ieigen++ ) {
    if ( eigen_analysis ) step_val = ieigen;
    else step_val = time_current;

    if ( npuknwn>0 ) {

      db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET );
      db( DOF_TYPE, 0, dof_type, ddum, ldum, VERSION_NORMAL, GET );
      db( DOF_AMOUNT, 0, dof_amount, ddum, ldum, VERSION_NORMAL, GET );
      db( DOF_SCAL_VEC_MAT, 0, dof_scal_vec_mat, ddum, ldum, VERSION_NORMAL, GET );

        // write scalars, vectors and tensors for primary unknowns
      ipuknwn = 0; ready = 0;
      while ( !ready ) {
        iuknwn = ipuknwn*nder;
        if      ( print_deformation ) {
          nval = ndim;
          data_type = 2;
          desc_comp = 1;
          strcpy( descr_menu, "mesh_deform" );
          outres << descr_menu << "                    " << load_type << " " << step_val;
          outres << " " << data_type << " " << data_loc << " " << desc_comp << "\n";
          for ( idim=0; idim<ndim; idim++ ) {
            if      ( idim==0 )
              outres << "mesh_deform-x" << "\n";
            else if ( idim==1 )
              outres << "mesh_deform-y" << "\n";
            else {
              assert( idim==2 );
              outres << "mesh_deform-z" << "\n";
            }
          }
        }
        else if ( dof_scal_vec_mat[iuknwn]==-SCALAR ) {
          nval = 1;
          data_type = 1;
          desc_comp = 1;
          strcpy( descr_menu, db_name(dof_label[iuknwn]) );
          string_shorten( descr_menu, 15 );
          outres << descr_menu << "                    " << load_type << " " << step_val;
          outres << " " << data_type << " " << data_loc << " " << desc_comp << "\n";
        }
        else if ( dof_scal_vec_mat[iuknwn]==-VECTOR ) {
          nval = dof_amount[iuknwn];
          data_type = 2;
          desc_comp = 1;
          if ( dof_type[iuknwn]==-MATERI_VELOCITY_INTEGRATED )
            strcpy( descr_menu, "materi_velint" );
          else
            strcpy( descr_menu, db_name(dof_type[iuknwn]) );
          string_shorten( descr_menu, 15 );
          outres << descr_menu << "                    " << load_type << " " << step_val;
          outres << " " << data_type << " " << data_loc << " " << desc_comp << "\n";
          for ( idim=0; idim<nval; idim++ ) {
            indx = iuknwn + idim*nder;
            outres << db_name(dof_label[indx]) << "\n";
          }
        }
        else {
          assert( dof_scal_vec_mat[iuknwn]==-MATRIX );
          nval = 6;
          data_type = 3;
          desc_comp = 1;
          strcpy( descr_menu, db_name(dof_type[iuknwn]) );
          string_shorten( descr_menu, 15 );
          outres << descr_menu << "                    " << load_type << " " << step_val;
          outres << " " << data_type << " " << data_loc << " " << desc_comp << "\n";
          for ( idim=0; idim<MDIM; idim++ ) {
            for ( jdim=idim; jdim<MDIM; jdim++ ) {
              indx = iuknwn + stress_indx(idim,jdim)*nder;
              outres << db_name(dof_label[indx]) << "\n";
            }
          }
        }
        for ( inod=1; inod<=max_node; inod++ ) {
          outres << inod << " ";
          if ( eigen_analysis )
            node_eigen = db_dbl( NODE_EIGEN, inod, VERSION_PRINT );
          else
            node_dof = db_dbl( NODE_DOF, inod, VERSION_PRINT );
          if      ( print_deformation ) {
            db( NODE, inod, idum, 
              coord, ldum, VERSION_PRINT, GET );
            db( NODE_START_REFINED, inod, idum, 
              coord_start_refined, ldum, VERSION_PRINT, GET );
            for ( idim=0; idim<ndim; idim++ ) {
              if ( materi_displacement )
                tmp = coord[idim] + node_dof[dis_indx+idim*nder];
              else
                tmp = coord[idim];
              outres << tmp-coord_start_refined[idim] << " ";
            }
          }
          else if ( dof_scal_vec_mat[iuknwn]==-SCALAR ) {
            if ( eigen_analysis )
              tmp = node_eigen[ieigen*nuknwn+iuknwn];
            else
              tmp = node_dof[iuknwn];
            outres << tmp;
          }
          else if ( dof_scal_vec_mat[iuknwn]==-VECTOR ) {
            n = dof_amount[iuknwn];
            for ( idim=0; idim<n; idim++ ) {
              indx = iuknwn+idim*nder;
              if ( eigen_analysis )
                tmp = node_eigen[ieigen*nuknwn+indx];
              else
                tmp = node_dof[indx];
              outres << tmp << " ";
            }
          }
          else {
            assert( dof_scal_vec_mat[iuknwn]==-MATRIX );
            for ( idim=0; idim<MDIM; idim++ ) {
              for ( jdim=idim; jdim<MDIM; jdim++ ) {
                indx = iuknwn + stress_indx(idim,jdim)*nder;
                if ( eigen_analysis )
                  tmp = node_eigen[ieigen*nuknwn+indx];
                else
                  tmp = node_dof[indx];
                outres << tmp << " ";
              }
            }
          }
          outres << "\n";
        }
        ipuknwn += nval;
        if ( ipuknwn>=npuknwn ) {
          if      ( materi_velocity && !print_deformation ) {
            print_deformation = 1;
          }
          else
            ready = 1;
        }
      }

      if ( !eigen_analysis && db_active_index( POST_CALCUL, 0, VERSION_NORMAL ) ) {

        db( POST_CALCUL_SCAL_VEC_MAT, 0, post_calcul_scal_vec_mat, 
          ddum, length_post_calcul_scal_vec_mat, VERSION_NORMAL, GET );

          // write scalars, vectors and tensors for calculated data
        icalcul = idim = ready = 0;
        while ( !ready ) {
          strcpy( outname, post_calcul_names_without_extension[icalcul] );
          if      ( post_calcul_scal_vec_mat[icalcul]==-SCALAR ) {
            nval = 1;
            data_type = 1;
            desc_comp = 0;
            strcpy( descr_menu, outname );
            string_shorten( descr_menu, 15 );
            outres << descr_menu << "                    " << load_type << " " << step_val;
            outres << " " << data_type << " " << data_loc << " " << desc_comp << "\n";
          }
          else if ( post_calcul_scal_vec_mat[icalcul]==-VECTOR ) {
            nval = MDIM;
            data_type = 2;
            desc_comp = 1;
            strcpy( descr_menu, outname );
            string_shorten( descr_menu, 15 );
            outres << descr_menu << "                    " << load_type << " " << step_val;
            outres << " " << data_type << " " << data_loc << " " << desc_comp << "\n";
            for ( idim=0; idim<MDIM; idim++ ) {
              outres << post_calcul_names[icalcul+idim] << "\n";
            }
          }
          else {
            assert( post_calcul_scal_vec_mat[icalcul]==-MATRIX );
            nval = 6;
            data_type = 3;
            desc_comp = 1;
            strcpy( descr_menu, outname );
            string_shorten( descr_menu, 15 );
            outres << descr_menu << "                    " << load_type << " " << step_val;
            outres << " " << data_type << " " << data_loc << " " << desc_comp << "\n";
            for ( idim=0; idim<MDIM; idim++ ) {
              for ( jdim=idim; jdim<MDIM; jdim++ ) {
                indx = icalcul + stress_indx(idim,jdim)*nder;
                outres << post_calcul_names[icalcul+indx] << "\n";
              }
            }
          }
          for ( inod=1; inod<=max_node; inod++ ) {
            outres << inod << " ";
            node_dof_calcul = db_dbl( NODE_DOF_CALCUL, inod, VERSION_PRINT );
            if      ( post_calcul_scal_vec_mat[icalcul]==-SCALAR ) {
              tmp = node_dof_calcul[icalcul];
              outres << tmp;
            }
            else if ( post_calcul_scal_vec_mat[icalcul]==-VECTOR ) {
              for ( idim=0; idim<MDIM; idim++ ) {
                indx = icalcul+idim;
                tmp = node_dof_calcul[indx];
                outres << tmp << " ";
              }
            }
            else {
              assert( post_calcul_scal_vec_mat[icalcul]==-MATRIX );
              for ( idim=0; idim<MDIM; idim++ ) {
                for ( jdim=idim; jdim<MDIM; jdim++ ) {
                  indx = icalcul + stress_indx(idim,jdim)*nder;
                  tmp = node_dof_calcul[indx];
                  outres << tmp << " ";
                }
              }
            }
            outres  << "\n";
          }
          icalcul += nval;
          ready = (icalcul>=length_post_calcul_scal_vec_mat);
        }

      }

    }

  }

  delete[] dof_amount;
  delete[] dof_label;
  delete[] dof_type;
  delete[] post_calcul_scal_vec_mat;
  delete[] dof_scal_vec_mat;
  delete[] nodes;
  delete[] el;
  delete[] element_tendon_intersections;
  delete[] control_eigen_values;

  outmesh.close();
  outbound.close();
  outres.close();
  db_version_delete( VERSION_PRINT );

  if ( swit ) pri( "Out routine PRINT_GID_5" );

  return 1;
}

