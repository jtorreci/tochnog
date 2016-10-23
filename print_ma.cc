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

void print_matlab( void )

{
  long int nnol=0, inod=0, ipuknwn=0, iuknwn=0, icontrol=0, 
    icalcul=0, ncalcul=0, swit=0, element=0, max_element=0, 
    itask=0, ntask=0, length=0, name=0, iplot=0, nplot=0,
    print_coord=0, print_unknown=0, print_calcul=0, ready=0, ldum=0, 
    *plot_nodes=NULL, *dof_label=NULL, *nodes=NULL, *el=NULL;
  double ddum[1], *coord=NULL, *node_dof=NULL, *node_dof_calcul=NULL;
  char str[MCHAR], x_name[MCHAR], y_name[MCHAR], z_name[MCHAR], unknown_name[MCHAR];

  swit = set_swit(-1,-1,"print_matlab");
  if ( swit ) pri( "In routine PRINT_MATLAB" );

  plot_nodes = get_new_int(MAXIMUM_NODE);
  dof_label = get_new_int(MUKNWN);
  nodes = get_new_int(MAXIMUM_NODE);
  el = get_new_int(1+MNOL);

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
  if ( db_max_index( POST_CALCUL_SCAL_VEC_MAT, ldum, VERSION_NORMAL, GET ) >=0 )
    ncalcul = db_len( POST_CALCUL_SCAL_VEC_MAT, 0, VERSION_NORMAL );
  long_to_a( icontrol, str );

  ntask = 1; // coordinates
  ntask += npuknwn; // primary unknowns
  ntask += ncalcul; // calcul results

    // empty the files
  for ( itask=0; itask<ntask; itask++ ) {
    if      ( itask>=0 && itask<1 ) {
      ;
    }
    else if ( itask>=1 && itask<1+npuknwn ) {
      ipuknwn = itask - 1;
      iuknwn = ipuknwn*nder;
    }
    else {
      icalcul = itask - (1+npuknwn);
    }
    strcpy( x_name, "x." );
    strcat( x_name, str );
    ofstream x_out( x_name );
    x_out.precision(TN_PRECISION);
    strcpy( y_name, "y." );
    strcat( y_name, str );
    ofstream y_out( y_name );
    y_out.precision(TN_PRECISION);
    strcpy( z_name, "z." );
    strcat( z_name, str );
    ofstream z_out( z_name );
    z_out.precision(TN_PRECISION);
    strcpy( unknown_name, db_name(dof_label[iuknwn]) );
    strcat( unknown_name, "." );
    strcat( unknown_name, str );
    ofstream unknown_out( unknown_name );
    unknown_out.precision(TN_PRECISION);
    strcpy( unknown_name, post_calcul_names[icalcul] );
    strcat( unknown_name, "." );
    strcat( unknown_name, str );
    ofstream calcul_out( unknown_name );
    calcul_out.precision(TN_PRECISION);
    x_out.close();
    y_out.close();
    z_out.close();
    unknown_out.close();
    calcul_out.close();
  }

    // print
  for ( itask=0; itask<ntask; itask++ ) {
    if ( swit ) pri( "itask", itask );
    print_coord = print_unknown = print_calcul = 0;
    if      ( itask>=0 && itask<1 ) {
      print_coord = 1;
      if ( swit ) pri( "print_coord" );
    }
    else if ( itask>=1 && itask<1+npuknwn ) {
      ipuknwn = itask - 1;
      iuknwn = ipuknwn*nder;
      print_unknown = 1;
      if ( swit ) pri( "print_unknown" );
    }
    else {
      icalcul = itask - (1+npuknwn);
      print_calcul = 1;
      if ( swit ) pri( "print_calcul" );
    }
    strcpy( x_name, "x." );
    strcat( x_name, str );
    ofstream x_out( x_name, ios::app );
    x_out.precision(TN_PRECISION);
    strcpy( y_name, "y." );
    strcat( y_name, str );
    ofstream y_out( y_name, ios::app );
    y_out.precision(TN_PRECISION);
    strcpy( z_name, "z." );
    strcat( z_name, str );
    ofstream z_out( z_name, ios::app );
    z_out.precision(TN_PRECISION);
    strcpy( unknown_name, db_name(dof_label[iuknwn]) );
    strcat( unknown_name, "." );
    strcat( unknown_name, str );
    ofstream unknown_out( unknown_name, ios::app );
    unknown_out.precision(TN_PRECISION);
    strcpy( unknown_name, post_calcul_names[icalcul] );
    strcat( unknown_name, "." );
    strcat( unknown_name, str );
    ofstream calcul_out( unknown_name, ios::app );
    calcul_out.precision(TN_PRECISION);
    ready = 0;
    for ( iplot=0; !ready; iplot++ ) {
      if ( swit ) pri( "iplot", iplot );
      for ( element=0; element<=max_element; element++ ) {
        if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) ) {
          if ( swit ) pri( "element", element );
          db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
          name = el[0]; nnol = length - 1; array_move( &el[1], nodes, nnol );
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
            pri( "Error: CONTROL_PRINT_MATLAB not available for ", name );
            exit(TN_EXIT_STATUS );
          }
          if ( swit ) pri( "plot_nodes", plot_nodes, nplot );
          inod = plot_nodes[iplot];
          if ( swit ) pri( "inod", inod );
          if ( print_coord ) {
            coord = db_dbl( NODE, inod, VERSION_NORMAL );
            if ( db_active_index( NODE_DOF, inod, VERSION_NORMAL ) )
              node_dof = db_dbl( NODE_DOF, inod, VERSION_NORMAL );
            if ( ndim>0 ) {
              if ( materi_displacement )
                x_out << coord[0]+node_dof[dis_indx+0*nder] << " ";
              else
                x_out << coord[0] << " ";
            }
            if ( ndim>1 ) {
              if ( materi_displacement )
                y_out << coord[1]+node_dof[dis_indx+1*nder] << " ";
              else
                y_out << coord[1] << " ";
            }
            if ( ndim>2 ) {
              if ( materi_displacement )
                z_out << coord[2]+node_dof[dis_indx+2*nder] << " ";
              else
                z_out << coord[2] << " ";
            }
          }
          else if ( print_unknown ) {
            node_dof = db_dbl( NODE_DOF, inod, VERSION_NORMAL );
            unknown_out << node_dof[iuknwn] << " ";
          }
          else {
            assert( print_calcul );
            node_dof_calcul = db_dbl( NODE_DOF_CALCUL, inod, VERSION_NORMAL );
            calcul_out << node_dof_calcul[icalcul] << " ";
          }
        }
      }
      ready = ( iplot==(nplot-1) );
      if      ( print_coord ) {
        x_out << "\n";
        y_out << "\n";
        z_out << "\n";
      }
      else if ( print_unknown ) {
        unknown_out << "\n";
      }
      else {
        assert( print_calcul );
        calcul_out << "\n";
      }
    }
    x_out.close();
    y_out.close();
    z_out.close();
    unknown_out.close();
    calcul_out.close();
  }

  delete[] plot_nodes;
  delete[] dof_label;
  delete[] nodes;
  delete[] el;

  if ( swit ) pri( "Out routine PRINT_MATLAB" );
}
