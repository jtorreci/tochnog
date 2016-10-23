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

void print_tecplot( long int ival[] )

{
  long int i=0, element=0, inol=0, idim=0, nplot=0, nnol=0, name=0, 
    length=0, inod=0, max_node=0, max_element=0, new_file=0, 
    print_tecplot_new_mesh=-NO, number_of_elements=0, icalcul=0,
    number_of_nodes=0, ipuknwn=0, iuknwn=0,
    ncalcul=0, ieigen=0, neigen=0, icontrol=0,
    swit=0, ldum=0, idum[1], *dof_label=NULL, *dof_principal=NULL,
    *el=NULL, *nodes=NULL, *plot_nodes=NULL;
  double time_current=0., ddum[1], node[MDIM], 
    *control_eigen_values=NULL, *node_eigen=NULL, 
    *node_dof=NULL, *node_dof_calcul=NULL;
  char outname[MCHAR];

  if ( ival[0]==-MESH || !db_active_index( CONTROL_PRINT_TECPLOT_MESH, 0, VERSION_NORMAL ) ) {
    print_tecplot_new_mesh = -YES; length=1;
    db( CONTROL_PRINT_TECPLOT_MESH, 0, &print_tecplot_new_mesh, ddum, 
      length, VERSION_NORMAL, PUT );
  }

    // tecplot wants sequential numbering starting from 1
  db_version_copy( VERSION_NORMAL, VERSION_PRINT );
  renumbering( VERSION_PRINT, NO, 1, 1, idum, idum );

  db_max_index( NODE, max_node, VERSION_PRINT, GET );
  db_max_index( ELEMENT, max_element, VERSION_PRINT, GET );

  if ( max_element<0 ) return;

  for ( inod=0; inod<=max_node; inod++ )
    if ( db_active_index( NODE, inod, VERSION_PRINT ) ) number_of_nodes++;

  for ( element=0; element<=max_element; element++ )
    if ( db_active_index( ELEMENT, element, VERSION_PRINT ) ) number_of_elements++;

  swit = set_swit(-1,-1,"print_tecplot");
  if ( swit ) pri( "In routine PRINT_TECPLOT" );

  dof_label = get_new_int(MUKNWN);
  dof_principal = get_new_int(MUKNWN);
  el = get_new_int(MAXIMUM_NODE+1);
  nodes = get_new_int(MAXIMUM_NODE);
  plot_nodes = get_new_int(MAXIMUM_NODE);
  control_eigen_values = get_new_dbl(DATA_ITEM_SIZE);
  node_eigen = get_new_dbl(DATA_ITEM_SIZE);

  if ( db_max_index( POST_CALCUL_SCAL_VEC_MAT, ldum, VERSION_NORMAL, GET ) >=0 )
    ncalcul = db_len( POST_CALCUL_SCAL_VEC_MAT, 0, VERSION_NORMAL );

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  db( DOF_PRINCIPAL, 0, dof_principal, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_EIGEN_VALUES, 0, idum, control_eigen_values,
    neigen, VERSION_NORMAL, GET_IF_EXISTS );

  ieigen = ival[0];
  if ( ieigen>=0 && ieigen>(neigen-1) ) db_error( CONTROL_PRINT_TECPLOT, icontrol );

  ifstream in( "tecplot.plt" );
  if ( !in ) new_file = 1;
  in.close();

  ofstream out( "tecplot.plt", ios::app );
  out.precision(TN_PRECISION);

    // general header
  if ( new_file ) {
    out << "TITLE = \"results with " << data_file << "\"\n";
    out << "VARIABLES = ";
    for ( idim=0; idim<ndim; idim++ ) {
      if      ( idim==0 ) 
        out << "\"x\"";
      else if ( idim==1 ) 
        out << ", \"y\"";
      else {
        assert( idim==2 );
        out << ", \"z\"";
      }
    }
    if ( npuknwn>0 ) {
      db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
      out << ", \"" << "mesh" << "\" ";
      for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
        iuknwn = ipuknwn*nder;
        out << ", \"" << db_name(dof_label[iuknwn]) << "\" ";
        for ( i=0; i<neigen; i++ ) {
          out << ", \"" << "Eigen mode " << control_eigen_values[i] <<
            " for " << db_name(dof_label[iuknwn]) << "\" ";
        }
      }
      for ( icalcul=0; icalcul<ncalcul; icalcul++ ) {
        strcpy( outname, post_calcul_names[icalcul] );
        out << ", \"" << outname << "\" ";
      }
    }
    out << "\n\n";
  }

    // zone header
  db( TIME_CURRENT, 0, idum, &time_current, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  out << "ZONE T=\"" << "t=" << time_current << "\"";
  if ( ndim==1 )
    out << ", F=POINT, I=" << number_of_nodes << "\n";
  else {
    out << ", F=FEPOINT, N=" << number_of_nodes << ", E=" << number_of_elements;
    if ( ndim==2 )
      out << ", ET=QUADRILATERAL\n";
    else {
      assert( ndim==3 );
      out << ", ET=BRICK\n";
    }
  }
  if ( ndim>1 && print_tecplot_new_mesh!=-YES ) {
    if ( ndim==2 )
      out << "D=(1,2,FECONNECT)\n";
    else {
      assert( ndim==3 );
      out << "D=(1,2,3,FECONNECT)\n";
    }
  }
  out << "\n";

    // nodal locations and primary unknowns
  for ( inod=0; inod<=max_node; inod++ ) {
    if ( db_active_index( NODE, inod, VERSION_PRINT ) ) {
      if ( db_active_index( NODE_DOF, inod, VERSION_PRINT ) )
        node_dof = db_dbl( NODE_DOF, inod, VERSION_PRINT );
      if ( db_active_index( NODE_DOF_CALCUL, inod, VERSION_PRINT ) )
        node_dof_calcul = db_dbl( NODE_DOF_CALCUL, inod, VERSION_PRINT );
      db( NODE_EIGEN, inod, idum, node_eigen, ldum, VERSION_PRINT, GET_IF_EXISTS );
      if ( print_tecplot_new_mesh==-YES ) {
        db( NODE, inod, idum, node, ldum, VERSION_PRINT, GET );
        if ( ieigen>=0 ) {
          db( NODE_EIGEN, inod, idum, node_eigen, ldum, VERSION_PRINT, GET );
          for ( idim=0; idim<ndim; idim++ )
            node[idim] += node_eigen[ieigen*nuknwn+vel_indx+idim*nder];
        }
        for ( idim=0; idim<ndim; idim++ ) {
          if ( materi_displacement )
            out << node[idim]+node_dof[dis_indx+idim*nder] << " ";
          else
            out << node[idim] << " ";
        }
      }
      out << "0.0" << " ";
      if ( npuknwn>0 ) {
        for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
          iuknwn = ipuknwn*nder;
          out << node_dof[iuknwn] << " ";
          for ( i=0; i<neigen; i++ )
            out << node_eigen[i*nuknwn+iuknwn] << " ";
        }
        if ( db_active_index( NODE_DOF_CALCUL, inod, VERSION_PRINT ) ) {
          length = db_len( NODE_DOF_CALCUL, inod, VERSION_PRINT );
          for ( i=0; i<length; i++ )
            out << node_dof_calcul[i] << " ";
        }
      }
      out << "\n";
    }
  }
  out << "\n";

    // connectivity
  if ( ndim>1 && print_tecplot_new_mesh==-YES ) {
    for ( element=0; element<=max_element; element++ ) {
      if ( db_active_index( ELEMENT, element, VERSION_PRINT ) ) {
        db( ELEMENT, element, el, ddum, length, VERSION_PRINT, GET );
        name = el[0];
        nnol = length - 1; array_move( &el[1], nodes, nnol );
        if      ( name==-BAR2 ) {
          nplot = 2;
          plot_nodes[0] = nodes[0];
          plot_nodes[1] = nodes[1];
        }
        else if ( name==-BAR3 ) {
          nplot = 2;
          plot_nodes[0] = nodes[0];
          plot_nodes[1] = nodes[2];
        }
        else if ( name==-BAR4 ) {
          nplot = 2;
          plot_nodes[0] = nodes[0];
          plot_nodes[1] = nodes[3];
        }
        else if ( name==-TRIA3 ) {
          nplot = 4;
          plot_nodes[0] = nodes[0];
          plot_nodes[1] = nodes[0];
          plot_nodes[2] = nodes[1];
          plot_nodes[3] = nodes[2];
        }
        else if ( name==-QUAD4 ) {
          nplot = 4;
          plot_nodes[0] = nodes[0];
          plot_nodes[1] = nodes[1];
          plot_nodes[2] = nodes[3];
          plot_nodes[3] = nodes[2];
        }
        else if ( name==-QUAD9 ) {
          nplot = 4;
          plot_nodes[0] = nodes[0];
          plot_nodes[1] = nodes[2];
          plot_nodes[2] = nodes[8];
          plot_nodes[3] = nodes[6];
        }
        else if ( name==-QUAD16 ) {
          nplot = 4;
          plot_nodes[0] = nodes[0];
          plot_nodes[1] = nodes[3];
          plot_nodes[2] = nodes[15];
          plot_nodes[3] = nodes[12];
        }
        else if ( name==-TET4 ) {
          nplot = 8;
          plot_nodes[0] = nodes[0];
          plot_nodes[1] = nodes[0];
          plot_nodes[2] = nodes[0];
          plot_nodes[3] = nodes[0];
          plot_nodes[4] = nodes[1];
          plot_nodes[5] = nodes[1];
          plot_nodes[6] = nodes[2];
          plot_nodes[7] = nodes[3];
        }
        else if ( name==-HEX8 ) {
          nplot = 8;
          plot_nodes[0] = nodes[0];
          plot_nodes[1] = nodes[1];
          plot_nodes[2] = nodes[3];
          plot_nodes[3] = nodes[2];
          plot_nodes[4] = nodes[4];
          plot_nodes[5] = nodes[5];
          plot_nodes[6] = nodes[7];
          plot_nodes[7] = nodes[6];
        }
        else if ( name==-HEX27 ) {
          nplot = 8;
          plot_nodes[0] = nodes[0];
          plot_nodes[1] = nodes[2];
          plot_nodes[2] = nodes[8];
          plot_nodes[3] = nodes[6];
          plot_nodes[4] = nodes[18];
          plot_nodes[5] = nodes[20];
          plot_nodes[6] = nodes[26];
          plot_nodes[7] = nodes[24];
        }
        else if ( name==-HEX64 ) {
          nplot = 8;
          plot_nodes[0] = nodes[0];
          plot_nodes[1] = nodes[3];
          plot_nodes[2] = nodes[15];
          plot_nodes[3] = nodes[12];
          plot_nodes[4] = nodes[48];
          plot_nodes[5] = nodes[51];
          plot_nodes[6] = nodes[63];
          plot_nodes[7] = nodes[60];
        }
        else {
          pri( "Error: CONTROL_PRINT_TECPLOT not available for ", name );
          exit(TN_EXIT_STATUS );
        }      
        for ( inol=0; inol<nplot; inol++ )
          out << plot_nodes[inol] << " ";
        out << "\n";
      }
    }
    out << "\n";
  }

  out.close();
  db_version_delete( VERSION_PRINT );

  delete[] dof_label;
  delete[] dof_principal;
  delete[] el;
  delete[] nodes;
  delete[] plot_nodes;
  delete[] control_eigen_values;
  delete[] node_eigen;

  if ( swit ) pri( "Out routine PRINT_TECPLOT" );
}
