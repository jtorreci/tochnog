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

void print_gmv( long int icontrol, long int ival[] )

{
  long int inod=0, element=0, max_node=0, max_element=0, nnol=0, name=0, max=0, ngroup=0,
    element_group=0, length=0, idim=0, jdim=0, ipuknwn=0, iuknwn=0,
    icalcul=0, ncalcul=0, ieigen=0, neigen=0,
    igroup=0, control_print_gmv_mesh=0, itrace=0, ntrace=0, 
    swit=0, ldum=0, idum[1], *groups=NULL, *dof_label=NULL, 
    *nodes=NULL, *el=NULL, *dof_principal=NULL;
  double tmp=0., time_current=0., ddum[1], post_point[MDIM], coord[MDIM], 
    *post_point_dof=NULL, *node_eigen=NULL, 
    *node_dof=NULL, *node_dof_calcul=NULL;
  char filename[MCHAR], str[MCHAR];

  db_version_copy( VERSION_NORMAL, VERSION_PRINT );
  renumbering( VERSION_PRINT, NO, 1, 1, idum, idum );

  db_highest_index( NODE, max_node, VERSION_PRINT );
  db_highest_index( ELEMENT, max_element, VERSION_PRINT );
  if ( ndim==1 || max_element<0 ) return;

  swit = set_swit(-1,-1,"print_gmv");
  if ( swit ) pri( "In routine PRINT_GMV" );

  groups = get_new_int(DATA_ITEM_SIZE);
  dof_label = get_new_int(MUKNWN);
  nodes = get_new_int(MAXIMUM_NODE);
  el = get_new_int(MAXIMUM_NODE+1);
  dof_principal = get_new_int(MUKNWN);
  post_point_dof = get_new_dbl(MUKNWN);
  node_eigen = get_new_dbl(DATA_ITEM_SIZE);

  db( CONTROL_PRINT_GMV_MESH, 0, &control_print_gmv_mesh, ddum, 
    length, VERSION_NORMAL, GET_IF_EXISTS );
  control_print_gmv_mesh++; length=1;
  db( CONTROL_PRINT_GMV_MESH, 0, &control_print_gmv_mesh, ddum, 
    length, VERSION_NORMAL, PUT );
  db( DOF_PRINCIPAL, 0, dof_principal, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );

  if ( db_max_index( POST_CALCUL_SCAL_VEC_MAT, ldum, VERSION_NORMAL, GET ) >=0 )
    ncalcul = db_len( POST_CALCUL_SCAL_VEC_MAT, 0, VERSION_NORMAL );
  if ( db_active_index( CONTROL_EIGEN_VALUES, 0, VERSION_NORMAL ) )
    neigen = db_len( CONTROL_EIGEN_VALUES, 0, VERSION_NORMAL );

  ieigen = ival[0];
  if ( ieigen>=0 && ieigen>(neigen-1) ) db_error( CONTROL_PRINT_PLOTMTV, icontrol );

  strcpy( filename, "gmv" );
  if ( icontrol>=0 ) {
    if      ( control_print_gmv_mesh<10 )
      strcat( filename, "00" );
    else if ( control_print_gmv_mesh<100 )
      strcat( filename, "0" );
    long_to_a( control_print_gmv_mesh, str );
    strcat( filename, str );
  }
  strcat( filename, ".inp" );
  ofstream out( filename );
  out.precision(TN_PRECISION);

  out << "gmvinput ascii\n\n";

  out << "nodes " << max_node << "\n";
  for ( idim=0; idim<MDIM; idim++ ) {
    for ( inod=1; inod<=max_node; inod++ ) {
      db( NODE, inod, idum, coord, ldum, VERSION_PRINT, GET );
      if ( ieigen>=0 ) {
        db( NODE_EIGEN, inod, idum, node_eigen, ldum, VERSION_PRINT, GET );
        for ( jdim=0; jdim<ndim; jdim++ )
          coord[jdim] += node_eigen[ieigen*nuknwn+vel_indx+jdim*nder];
      }
      if      ( ndim==2 && idim==2 )
        out << "0." << " ";
      else if ( materi_displacement ) {
        node_dof = db_dbl( NODE_DOF, inod, VERSION_PRINT );
        out << coord[idim]+node_dof[dis_indx+idim*nder] << " ";
      }
      else
        out << coord[idim] << " ";
    }
    out << "\n";
  } 
  out << "\n";

  out << "cells " << max_element << "\n";
  for ( element=1; element<=max_element; element++ ) {
    db( ELEMENT, element, el, ddum, length, VERSION_PRINT, GET );
    name = el[0];
    nnol = length - 1; array_move( &el[1], nodes, nnol );
    if      ( name==-TRIA3 ) {
      out << "tri 3\n";
      out << nodes[0] << " " << nodes[1] << " " << nodes[2];
    }
    else if ( name==-QUAD4 ) {
      out << "quad 4\n";
      out << nodes[0] << " " << nodes[1] << " " << nodes[3] << " " << nodes[2];
    }
    else if ( name==-QUAD9 ) {
      out << "quad 4\n";
      out << nodes[0] << " " << nodes[2] << " " << nodes[8] << " " << nodes[6];
    }
    else if ( name==-TET4 ) {
      out << "tet 4\n";
      out << nodes[0] << " " << nodes[1] << " " << nodes[2] << " " << nodes[3];
    }
    else if ( name==-HEX8 ) {
      out << "hex 8\n";
      out << nodes[0] << " " << nodes[1] << " " << nodes[3] << " " << nodes[2] << " " <<
             nodes[4] << " " << nodes[5] << " " << nodes[7] << " " << nodes[6];
    }
    else if ( name==-HEX27 ) {
      out << "hex 8\n";
      out << nodes[0] << " " << nodes[2] << " " << nodes[8] << " " << nodes[6] << " " <<
             nodes[18] << " " << nodes[20] << " " << nodes[26] << " " << nodes[24];
    }
    else {
      pri( "Error: illegal element type detected for CONTROL_PRINT_GMV.\n");
      exit(TN_EXIT_STATUS);
    }
    out << "\n";
  }
  out << "\n";

  if ( materi_velocity ) {
    out << "velocity 1\n";
    for ( idim=0; idim<MDIM; idim++ ) {
      for ( inod=1; inod<=max_node; inod++ ) {
        node_dof = db_dbl( NODE_DOF, inod, VERSION_PRINT );
        if      ( ndim==2 && idim==2 )
          out << "0." << " ";
        else {
          node_dof = db_dbl( NODE_DOF, inod, VERSION_PRINT );
          out << node_dof[vel_indx+idim*nder] << " ";
        }
      }
      out << "\n";
    } 
    out << "\n";
  }

  if ( npuknwn>0 ) {
    db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET );
    out << "variable\n";
    for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
      iuknwn = ipuknwn*nder;
      out << db_name(dof_label[iuknwn]) << " 1\n";
      for ( inod=1; inod<=max_node; inod++ ) {
        node_dof = db_dbl( NODE_DOF, inod, VERSION_PRINT );
        out << node_dof[iuknwn] << " ";
      }
      if ( neigen>0 ) {
        if ( dof_principal[iuknwn]>=0 ) {
          for ( ieigen=0; ieigen<neigen; ieigen++ ) {
            long_to_a( ieigen, str );
            out << "E" << str << "_" << db_name(dof_label[iuknwn]) << " 1\n";
            for ( inod=1; inod<=max_node; inod++ ) {
              db( NODE_EIGEN, inod, idum, node_eigen, ldum, VERSION_PRINT, GET );
              out << node_eigen[ieigen*nuknwn+iuknwn] << " ";
            }
          }
        }
      }
      out << "\n";
    }
    for ( icalcul=0; icalcul<ncalcul; icalcul++ ) {
      out << post_calcul_names[icalcul] << " 1\n";
      for ( inod=1; inod<=max_node; inod++ ) {
        node_dof_calcul = db_dbl( NODE_DOF_CALCUL, inod, VERSION_PRINT );
        out << node_dof_calcul[icalcul] << " ";
      }
      out << "\n";
    }
    out << "endvars\n";
    out << "\n";
  }

  if ( db_max_index( ELEMENT_GROUP, max, VERSION_PRINT, GET ) >= 0 ) {
    array_set( groups, 0, DATA_ITEM_SIZE );
    for ( element=1; element<=max_element; element++ ) {
      if ( db_active_index( ELEMENT_GROUP, element, VERSION_PRINT ) )
        db( ELEMENT_GROUP, element, &element_group, 
          ddum, ldum, VERSION_PRINT, GET );
      else element_group = 0;
      groups[element_group] = 1;
    }
    for ( igroup=0; igroup<DATA_ITEM_SIZE; igroup++ ) {
      if ( groups[igroup] ) {
        ngroup++;
        groups[igroup] = ngroup;
      }
    }
    out << "material " << ngroup << " 0\n";
    for ( igroup=0; igroup<DATA_ITEM_SIZE; igroup++ ) {
      if ( groups[igroup] ) out << "group_" << igroup << " ";
    }
    out << "\n";
    for ( element=1; element<=max_element; element++ ) {
      if ( db_active_index( ELEMENT_GROUP, element, VERSION_PRINT ) )
        db( ELEMENT_GROUP, element, &element_group, 
          ddum, ldum, VERSION_PRINT, GET );
      else element_group = 0;
      out << groups[element_group] << " ";
    }
    out << "\n\n";
  }

  if ( db_max_index( POST_POINT, max, VERSION_NORMAL, GET ) >= 0 ) {
    for ( itrace=0; itrace<=max; itrace++ ) {
      if ( db_active_index( POST_POINT, itrace, VERSION_NORMAL ) ) ntrace++;
    }
    out << "tracers " << ntrace << "\n";
    for ( idim=0; idim<MDIM; idim++ ) {
      for ( itrace=0; itrace<=max; itrace++ ) {
        if ( db_active_index( POST_POINT, itrace, VERSION_NORMAL ) ) {
          if ( idim>ndim-1 )
            tmp = 0.;
          else {
            db( POST_POINT, itrace, idum, post_point, 
              ldum, VERSION_NORMAL, GET );
            tmp = post_point[idim];
          }
          out << tmp << " ";
          ntrace++;
        }
      }
      out << "\n";
    }
    for ( ipuknwn=0; ipuknwn<npuknwn && ipuknwn<19; ipuknwn++ ) {
      iuknwn = ipuknwn*nder;
      out << db_name(dof_label[iuknwn]) << "\n";
      for ( itrace=0; itrace<=max; itrace++ ) {
        if ( db_active_index( POST_POINT_DOF, itrace, VERSION_NORMAL ) ) {
          db( POST_POINT_DOF, itrace, idum, post_point_dof, 
            ldum, VERSION_NORMAL, GET );
          out << post_point_dof[iuknwn] << " ";
        }
      }
      out << "\n";
    }
    out << "endtrace\n\n";
  }

  db( TIME_CURRENT, 0, idum, &time_current, ldum, VERSION_NORMAL, GET );
  out << "probtime " << time_current << "\n";
  out << "\n";

  out << "endgmv";
  out.close();

  db_version_delete( VERSION_PRINT );
  delete[] groups;
  delete[] dof_label;
  delete[] nodes;
  delete[] el;
  delete[] dof_principal;
  delete[] post_point_dof;
  delete[] node_eigen;

  if ( swit ) pri( "Out routine PRINT_GMV" );
}

// print_gmsh - control_print_gmsh: Gmsh 2.2 ASCII output (.msh).
// Writes the mesh once (nodes + elements, plus dummy point elements for
// vector plots) and appends NodeData / ElementData per time step.
// switch task: -yes (single <base>.msh, mesh written only the first time),
// -separate_index (<base><icontrol>.msh), -separate_sequential
// (<base><n>.msh with increasing n).
void print_gmsh( long int icontrol, long int task )

{
  long int inod=0, element=0, max_node=0, max_element=0, nnol=0, name=0,
    length=0, idim=0, jdim=0, kdim=0, ldim=0, ipuknwn=0, iuknwn=0,
    nder_=0, nuknwn_=0, element_group=0, swit=0, ldum=0, first=1,
    dummy=-YES, element_data=-YES, node_method=-NODE;
  long int idum[1], *dof_label=NULL, *dof_scal_vec_mat=NULL, *nodes=NULL,
    *el=NULL;
  double ddum[1], time_current=0., coord[MDIM], *node_dof=NULL;
  char filename[MCHAR], str[MCHAR];

  swit = set_swit(-1,-1,"print_gmsh");
  if ( swit ) pri( "In routine PRINT_GMSH" );

  // options (per control record): control_print_gmsh_dummy (default -yes),
  // gmsh_element_data (default -yes -> ElementData; -no -> ElementNodeData),
  // gmsh_node_method (default -node).
  db( CONTROL_PRINT_GMSH_DUMMY, icontrol, &dummy, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_PRINT_GMSH_ELEMENT_DATA, icontrol, &element_data, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_PRINT_GMSH_NODE_METHOD, icontrol, &node_method, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );

  db_version_copy( VERSION_NORMAL, VERSION_PRINT );
  renumbering( VERSION_PRINT, NO, 0, 0, idum, idum );
  db_highest_index( NODE, max_node, VERSION_PRINT );
  db_highest_index( ELEMENT, max_element, VERSION_PRINT );
  if ( max_node<0 || max_element<0 ) return;
  db( TIME_CURRENT, 0, idum, &time_current, ldum, VERSION_NORMAL, GET );
  nder_ = nder;
  nuknwn_ = nuknwn;

  dof_label = get_new_int(MUKNWN);
  dof_scal_vec_mat = get_new_int(MUKNWN);
  nodes = get_new_int(MAXIMUM_NODE);
  el = get_new_int(MAXIMUM_NODE+1);
  db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( DOF_SCAL_VEC_MAT, 0, dof_scal_vec_mat, ddum, ldum, VERSION_NORMAL,
    GET_IF_EXISTS );

  // file name: -yes -> <base>.msh (mesh once); -separate_index ->
  // <base><icontrol>.msh; -separate_sequential -> <base><n>.msh.
  strcpy( filename, data_file_base );
  if      ( task==-SEPARATE_INDEX && icontrol>=0 ) {
    long_to_a( icontrol, str );
    strcat( filename, str );
  }
  else if ( task==-SEPARATE_SEQUENTIAL ) {
    static long int gmsh_seq=0;
    long_to_a( gmsh_seq++, str );
    strcat( filename, str );
  }
  strcat( filename, ".msh" );

  {
    std::ifstream fexists( filename );
    first = !fexists.is_open();
    fexists.close();
  }

  std::ofstream out( filename, std::ios::app );
  out.precision(TN_PRECISION);

  if ( first ) {
    out << "$MeshFormat\n";
    out << "2.2 0 8\n";
    out << "$EndMeshFormat\n";

    out << "$Nodes\n" << max_node+1 << "\n";
    for ( inod=0; inod<=max_node; inod++ ) {
      db( NODE, inod, idum, coord, ldum, VERSION_PRINT, GET );
      if ( node_method==-NODE_DEFORMED_MESH && materi_displacement ) {
        node_dof = db_dbl( NODE_DOF, inod, VERSION_PRINT );
        out << inod+1;
        for ( idim=0; idim<ndim; idim++ )
          out << " " << coord[idim]+node_dof[dis_indx+idim*nder_];
        for ( idim=ndim; idim<MDIM; idim++ ) out << " 0";
      }
      else if ( node_method==-NODE_START_REFINED &&
                db_active_index( NODE_START_REFINED, inod, VERSION_PRINT ) ) {
        db( NODE_START_REFINED, inod, idum, coord, ldum, VERSION_PRINT, GET );
        out << inod+1;
        for ( idim=0; idim<MDIM; idim++ ) out << " " << coord[idim];
      }
      else {
        out << inod+1;
        for ( idim=0; idim<ndim; idim++ ) out << " " << coord[idim];
        for ( idim=ndim; idim<MDIM; idim++ ) out << " 0";
      }
      out << "\n";
    }
    out << "$EndNodes\n";

    // element connectivity (Gmsh 2.2 types; node order matches print_vtk)
    out << "$Elements\n";
    // dummy point element in each node for vector-field plots (group 1234)
    long int ntotal = max_element+1;
    if ( dummy!=-NO ) ntotal += max_node+1;
    out << ntotal << "\n";
    long int nelem=0;
    for ( element=0; element<=max_element; element++ ) {
      if ( !db_active_index( ELEMENT, element, VERSION_PRINT ) ) continue;
      db( ELEMENT, element, el, ddum, length, VERSION_PRINT, GET );
      name = el[0];
      nnol = length - 1; array_move( &el[1], nodes, nnol );
      if ( db_active_index( ELEMENT_GROUP, element, VERSION_PRINT ) )
        db( ELEMENT_GROUP, element, &element_group, ddum, ldum,
          VERSION_PRINT, GET );
      else
        element_group = 0;
      if      ( name==-BAR2 ) {
        nelem++;
        out << nelem << " 1 2 " << element_group << " " << element_group;
        out << " " << nodes[0]+1 << " " << nodes[1]+1 << "\n";
      }
      else if ( name==-TRIA3 ) {
        nelem++;
        out << nelem << " 2 2 " << element_group << " " << element_group;
        out << " " << nodes[0]+1 << " " << nodes[1]+1 << " " << nodes[2]+1 << "\n";
      }
      else if ( name==-QUAD4 ) {
        nelem++;
        out << nelem << " 3 2 " << element_group << " " << element_group;
        out << " " << nodes[0]+1 << " " << nodes[1]+1 << " " << nodes[3]+1
            << " " << nodes[2]+1 << "\n";
      }
      else if ( name==-TET4 ) {
        nelem++;
        out << nelem << " 4 2 " << element_group << " " << element_group;
        out << " " << nodes[0]+1 << " " << nodes[1]+1 << " " << nodes[2]+1
            << " " << nodes[3]+1 << "\n";
      }
      else if ( name==-HEX8 ) {
        nelem++;
        out << nelem << " 5 2 " << element_group << " " << element_group;
        out << " " << nodes[0]+1 << " " << nodes[1]+1 << " " << nodes[3]+1
            << " " << nodes[2]+1 << " " << nodes[4]+1 << " " << nodes[5]+1
            << " " << nodes[7]+1 << " " << nodes[6]+1 << "\n";
      }
    }
    if ( dummy!=-NO ) {
      for ( inod=0; inod<=max_node; inod++ ) {
        nelem++;
        out << nelem << " 15 2 1234 1234 " << inod+1 << "\n";
      }
    }
    out << "$EndElements\n";
  }

  // node / element data for each exported dof (append per time step)
  if ( nuknwn_>0 ) {
    for ( ipuknwn=0; ipuknwn<nuknwn_; ipuknwn++ ) {
      if ( dof_scal_vec_mat[ipuknwn]!=-SCALAR &&
           dof_scal_vec_mat[ipuknwn]!=-VECTOR &&
           dof_scal_vec_mat[ipuknwn]!=-MATRIX ) continue;
      long int nval = 1;
      if      ( dof_scal_vec_mat[ipuknwn]==-VECTOR ) nval = ndim;
      else if ( dof_scal_vec_mat[ipuknwn]==-MATRIX ) nval = 6;
      long int base_indx = ipuknwn*nder_;
      for ( idim=0; idim<nval; idim++ ) {
        char label[MCHAR];
        if      ( dof_scal_vec_mat[ipuknwn]==-SCALAR )
          strcpy( label, db_name(dof_label[ipuknwn]) );
        else if ( dof_scal_vec_mat[ipuknwn]==-VECTOR ) {
          sprintf( label, "%s_%ld", db_name(dof_label[ipuknwn]), idim );
        }
        else {
          // matrix: idim 0..5 = xx,yy,zz,xy,xz,yz (Voigt)
          const char* comp[6] = { "xx","yy","zz","xy","xz","yz" };
          sprintf( label, "%s_%s", db_name(dof_label[ipuknwn]), comp[idim] );
        }
        // node_* data
        out << "$NodeData\n";
        out << "1\n\"node_" << label << "\"\n1\n" << time_current
            << "\n3\n0\n0\n1\n" << max_node+1 << "\n";
        for ( inod=0; inod<=max_node; inod++ ) {
          node_dof = db_dbl( NODE_DOF, inod, VERSION_PRINT );
          long int indx = base_indx;
          if      ( dof_scal_vec_mat[ipuknwn]==-VECTOR )
            indx = base_indx + idim*nder_;
          else if ( dof_scal_vec_mat[ipuknwn]==-MATRIX ) {
            long int kk, ll;
            if      ( idim==0 ) { kk=0; ll=0; }
            else if ( idim==1 ) { kk=1; ll=1; }
            else if ( idim==2 ) { kk=2; ll=2; }
            else if ( idim==3 ) { kk=0; ll=1; }
            else if ( idim==4 ) { kk=0; ll=2; }
            else                 { kk=1; ll=2; }
            indx = base_indx + stress_indx(kk,ll)*nder_;
          }
          out << inod+1 << " " << node_dof[indx] << "\n";
        }
        out << "$EndNodeData\n";

        // element_* data: averaged over the element (ElementData) or
        // per element node (ElementNodeData), per gmsh_element_data.
        long int nact=0;
        for ( element=0; element<=max_element; element++ )
          if ( db_active_index( ELEMENT, element, VERSION_PRINT ) ) nact++;
        out << "$" << ( element_data==-NO ? "ElementNodeData" : "ElementData" )
            << "\n";
        out << "1\n\"element_" << label << "\"\n1\n" << time_current
            << "\n3\n0\n0\n1\n" << nact << "\n";
        for ( element=0; element<=max_element; element++ ) {
          if ( !db_active_index( ELEMENT, element, VERSION_PRINT ) ) continue;
          db( ELEMENT, element, el, ddum, length, VERSION_PRINT, GET );
          name = el[0];
          nnol = length - 1; array_move( &el[1], nodes, nnol );
          double sum=0.;
          for ( inod=0; inod<nnol; inod++ ) {
            node_dof = db_dbl( NODE_DOF, nodes[inod], VERSION_PRINT );
            long int indx = base_indx;
            if      ( dof_scal_vec_mat[ipuknwn]==-VECTOR )
              indx = base_indx + idim*nder_;
            else if ( dof_scal_vec_mat[ipuknwn]==-MATRIX ) {
              long int kk, ll;
              if      ( idim==0 ) { kk=0; ll=0; }
              else if ( idim==1 ) { kk=1; ll=1; }
              else if ( idim==2 ) { kk=2; ll=2; }
              else if ( idim==3 ) { kk=0; ll=1; }
              else if ( idim==4 ) { kk=0; ll=2; }
              else                 { kk=1; ll=2; }
              indx = base_indx + stress_indx(kk,ll)*nder_;
            }
            sum += node_dof[indx];
          }
          double avg = ( nnol>0 ) ? sum/((double)nnol) : 0.;
          if ( element_data==-NO ) {
            out << element+1 << " " << nnol;
            for ( inod=0; inod<nnol; inod++ ) {
              node_dof = db_dbl( NODE_DOF, nodes[inod], VERSION_PRINT );
              long int indx = base_indx;
              if      ( dof_scal_vec_mat[ipuknwn]==-VECTOR )
                indx = base_indx + idim*nder_;
              else if ( dof_scal_vec_mat[ipuknwn]==-MATRIX ) {
                long int kk, ll;
                if      ( idim==0 ) { kk=0; ll=0; }
                else if ( idim==1 ) { kk=1; ll=1; }
                else if ( idim==2 ) { kk=2; ll=2; }
                else if ( idim==3 ) { kk=0; ll=1; }
                else if ( idim==4 ) { kk=0; ll=2; }
                else                 { kk=1; ll=2; }
                indx = base_indx + stress_indx(kk,ll)*nder_;
              }
              out << " " << nodes[inod]+1 << " " << node_dof[indx];
            }
            out << "\n";
          }
          else {
            out << element+1 << " " << avg << "\n";
          }
        }
        out << "$End" << ( element_data==-NO ? "ElementNodeData" : "ElementData" )
            << "\n";
      }
    }
  }

  out.close();

  db_version_delete( VERSION_PRINT );
  delete[] dof_label;
  delete[] dof_scal_vec_mat;
  delete[] nodes;
  delete[] el;

  if ( swit ) pri( "Out routine PRINT_GMSH" );
}
