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
