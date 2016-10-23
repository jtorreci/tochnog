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

#define MCURVE 8
#define SOLID 1
#define DOTTED 3

void print_plotmtv( long int icontrol, long iarg[] )

{
  long int ival=0, nval=0, element=0, inol=0, nnol=0, inod=0, print_mesh=0,
    print_displacement=0, print_velocity=0, print_unknown=0,
    print_mesh_deformed=0, print_calcul_unknown=0, print_calcul_vector=0,
    print_eigen_unknown=0, idim=0, length=0, indx=0,
    name=0, max_element=0, ipuknwn=0, iuknwn=0, icalcul=0, 
    icurve=0, ncurve=0, ncurve_nodes=0, ncalcul=0,
    ieigen=0, neigen=0, eigen_mode=0,
    control_print_plotmtv_mesh=0, swit=0, ldum=0, idum[1], 
    *dof_principal=NULL, *dof_label=NULL, *el=NULL, 
    *nodes=NULL, *curve_nodes=NULL, *post_calcul_scal_vec_mat=NULL;
  double time_current=0., ddum[1], work[MDIM], node[MDIM], 
    *control_eigen_values=NULL, *node_eigen=NULL, 
    *node_dof=NULL, *node_dof_calcul=NULL;
  char filename[MCHAR], str[MCHAR];

  db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
  if ( max_element<0 ) return;

  swit = set_swit(-1,-1,"print_plotmtv");
  if ( swit ) pri( "In routine PRINT_PLOTMTV" );

  dof_principal = get_new_int(MUKNWN);
  dof_label = get_new_int(MUKNWN);
  el = get_new_int(MNOL+1);
  nodes = get_new_int(MNOL);
  curve_nodes = get_new_int(MCURVE*MNOL);
  post_calcul_scal_vec_mat = get_new_int(MUKNWN);
  control_eigen_values = get_new_dbl(DATA_ITEM_SIZE);
  node_eigen = get_new_dbl(DATA_ITEM_SIZE);

  db( POST_CALCUL_SCAL_VEC_MAT, 0, post_calcul_scal_vec_mat, ddum, ncalcul, 
     VERSION_NORMAL, GET_IF_EXISTS );
  db( TIME_CURRENT, 0, idum, &time_current, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( DOF_PRINCIPAL, 0, dof_principal, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );

  db( CONTROL_PRINT_PLOTMTV_MESH, 0, &control_print_plotmtv_mesh, ddum, 
    length, VERSION_NORMAL, GET_IF_EXISTS );
  control_print_plotmtv_mesh++; length=1;
  db( CONTROL_PRINT_PLOTMTV_MESH, 0, &control_print_plotmtv_mesh, ddum, 
    length, VERSION_NORMAL, PUT );
  db( CONTROL_EIGEN_VALUES, 0, idum, control_eigen_values,
    neigen, VERSION_NORMAL, GET_IF_EXISTS );

  ieigen = iarg[0];
  if ( ieigen>=0 && ieigen>(neigen-1) ) db_error( CONTROL_PRINT_PLOTMTV, icontrol );

  strcpy( filename, "plotmtv" );
  if ( icontrol>=0 ) {
    long_to_a( control_print_plotmtv_mesh, str );
    strcat( filename, str );
  }
  strcat( filename, ".dat" );
  ofstream out( filename );
  out.precision(TN_PRECISION);

  nval += 1; // mesh
  nval += 1; // deformed mesh
  nval += 1; // displacement vector
  nval += 1; // velocity vector
  nval += npuknwn; // unknowns
  nval += ncalcul; // calcul
  nval += neigen; // deformed mesh eigen
  nval += neigen*npuknwn; // eigen modes of principal unknowns

  for ( ival=0; ival<nval; ival++ ) {

    print_mesh = print_mesh_deformed = print_displacement =
      print_velocity = print_unknown = print_calcul_unknown = 
      print_calcul_vector = print_eigen_unknown = 0;
    indx = 0;
    if ( ival==indx ) {
      print_mesh = 1;
    }
    indx += 1;
    if ( ival==indx && materi_displacement ) {
      print_mesh_deformed = 1;
    }
    indx += 1;
    if ( ival==indx && materi_displacement ) {
      print_displacement = 1;
    }
    indx += 1;
    if ( ival==indx && materi_velocity ) {
      print_velocity = 1;
    }
    indx += 1;
    if ( ival>=indx && ival<indx+npuknwn && npuknwn>0 ) {
      ipuknwn = ival - indx;
      iuknwn = ipuknwn * nder;
      if ( ndim<MDIM ) print_unknown = 1;
    }
    indx += npuknwn;
    if ( ival>=indx && ival<indx+ncalcul && ncalcul>0 ) {
      icalcul = ival - indx;
      if      ( post_calcul_scal_vec_mat[icalcul]==-SCALAR ) {
        if ( ndim<MDIM ) print_calcul_unknown = 1;
      }
      else if ( post_calcul_scal_vec_mat[icalcul]==-VECTOR )
        print_calcul_vector = 1;
    }
    indx += ncalcul;
    if ( ival>=indx && ival<indx+neigen*npuknwn && neigen>0 ) {
      ipuknwn = ival - indx - eigen_mode*npuknwn;
      if ( ival==indx )
        eigen_mode = 0;
      else if ( ipuknwn==npuknwn ) {
        ipuknwn = 0;
        eigen_mode++;
      }
      iuknwn = ipuknwn * nder;
      if ( ndim<MDIM && eigen_mode<neigen && dof_principal[iuknwn]>=0 )
        print_eigen_unknown = 1;
    }
    indx += neigen*npuknwn;

    if ( print_mesh || print_mesh_deformed || print_displacement ||
         print_velocity || print_unknown || print_calcul_unknown ||
         print_calcul_vector || print_eigen_unknown ) {

      if      ( print_mesh ) {
        out << "$ DATA=CURVE3D\n";
        out << "% linetype=" << SOLID << "\n";
        out << "% equalscale=true\n";
        if ( ndim!=3 ) out << "% applyfill=false\n";
      }
      else if ( print_mesh_deformed ) {
        out << "$ DATA=CURVE3D\n";
        out << "% linetype=" << SOLID << "\n";
        out << "% equalscale=true\n";
        if ( ndim!=3 ) out << "% applyfill=false\n";
      }
      else if ( print_displacement || print_velocity ) {
        out << "$ DATA=VECTOR\n";
        out << "% equalscale=true\n";
        out << "% vlog=false\n";
      }
      else if ( print_calcul_vector ) {
        out << "$ DATA=VECTOR\n";
        out << "% equalscale=true\n";
        out << "% vlog=false\n";
        out << "% vhead=false\n";
      }
      else if ( print_unknown || print_calcul_unknown ||
                print_eigen_unknown ) {
        if ( ndim==1 )
          out << "$ DATA=CURVE2D\n"; 
        else if ( ndim==2 ) {
          out << "$ DATA=CONTCURVE\n";
          out << "% contstyle=2\n";
          out << "% contfill=true\n";
          out << "% equalscale=true\n";
        }
      }

      out << "% toplabel=\"";
      if      ( print_unknown )
        out << db_name(dof_label[iuknwn]);
      else if ( print_calcul_unknown )
        out << "calcul " << post_calcul_names[icalcul];
      else if ( print_eigen_unknown )
        out << "eigen mode " << control_eigen_values[eigen_mode] << 
          " for " << db_name(dof_label[iuknwn]);
      else if ( print_calcul_vector ) {
        strcpy( str, post_calcul_names[icalcul] );
        length = strlen( str );
        str[length-1] = '\0';
        out << str;
      }
      else if ( print_displacement )
        out << "displacement vectors";
      else if ( print_velocity )
        out << "velocity vectors";
      else if ( print_mesh )
        out << "mesh";
      else if ( print_mesh_deformed )
        out << "deformed mesh";
      out << " at time_current=" << time_current << "\"\n";

      out << "% xlabel=\" \"\n";
      out << "% ylabel=\" \"\n";
      out << "% zlabel=\" \"\n";
      out << "\n";

      for ( element=0; element<=max_element; element++ ) {
        if ( db_active_index(ELEMENT,element,VERSION_NORMAL) ) {
          db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
          name = el[0];
          nnol = length - 1; array_move( &el[1], nodes, nnol );
          ncurve = 0;
          if      ( ndim==1 ) {
            ncurve = 1; ncurve_nodes = nnol;
            array_move( nodes, curve_nodes, nnol );
          }
          else if ( name==-TRIA3 ) {
            if ( print_unknown || print_calcul_unknown || print_eigen_unknown ) {
              ncurve = 1; ncurve_nodes = 3;
              curve_nodes[0] = nodes[0];
              curve_nodes[1] = nodes[1];
              curve_nodes[2] = nodes[2];
            }
            else {
              ncurve = 1; ncurve_nodes = 4;
              curve_nodes[0] = nodes[0];
              curve_nodes[1] = nodes[1];
              curve_nodes[2] = nodes[2];
              curve_nodes[3] = nodes[0];

            }
          }
          else if ( name==-TRIA6 ) {
            if ( print_unknown || print_calcul_unknown || print_eigen_unknown ) {
              ncurve = 1; ncurve_nodes = 3;
              curve_nodes[0] = nodes[0];
              curve_nodes[1] = nodes[2];
              curve_nodes[2] = nodes[5];
            }
            else {
              ncurve = 1; ncurve_nodes = 4;
              curve_nodes[0] = nodes[0];
              curve_nodes[1] = nodes[2];
              curve_nodes[2] = nodes[5];
              curve_nodes[3] = nodes[0];

            }
          }
          else if ( name==-QUAD4 ) {
            if ( print_unknown || print_calcul_unknown || print_eigen_unknown ) {
              ncurve = 4; ncurve_nodes = 3;

              curve_nodes[0+0] = nodes[0];
              curve_nodes[0+1] = nodes[1];
              curve_nodes[0+2] = -1;

              curve_nodes[3+0] = nodes[1];
              curve_nodes[3+1] = nodes[3];
              curve_nodes[3+2] = -1;

              curve_nodes[6+0] = nodes[3];
              curve_nodes[6+1] = nodes[2];
              curve_nodes[6+2] = -1;

              curve_nodes[9+0] = nodes[2];
              curve_nodes[9+1] = nodes[0];
              curve_nodes[9+2] = -1;

            }
            else {
              ncurve = 1; ncurve_nodes = 5;

              curve_nodes[0] = nodes[0];
              curve_nodes[1] = nodes[1];
              curve_nodes[2] = nodes[3];
              curve_nodes[3] = nodes[2];
              curve_nodes[4] = nodes[0];

            }

          }
          else if ( name==-QUAD9 ) {
            if ( print_unknown || print_calcul_unknown || print_eigen_unknown ) {
              ncurve = 8; ncurve_nodes = 3;

              curve_nodes[0+0] = nodes[0];
              curve_nodes[0+1] = nodes[1];
              curve_nodes[0+2] = nodes[4];

              curve_nodes[3+0] = nodes[0];
              curve_nodes[3+1] = nodes[4];
              curve_nodes[3+2] = nodes[3];

              curve_nodes[6+0] = nodes[1];
              curve_nodes[6+1] = nodes[2];
              curve_nodes[6+2] = nodes[5];

              curve_nodes[9+0] = nodes[1];
              curve_nodes[9+1] = nodes[5];
              curve_nodes[9+2] = nodes[4];

              curve_nodes[12+0] = nodes[3];
              curve_nodes[12+1] = nodes[4];
              curve_nodes[12+2] = nodes[7];

              curve_nodes[15+0] = nodes[3];
              curve_nodes[15+1] = nodes[7];
              curve_nodes[15+2] = nodes[6];

              curve_nodes[18+0] = nodes[4];
              curve_nodes[18+1] = nodes[5];
              curve_nodes[18+2] = nodes[8];

              curve_nodes[21+0] = nodes[4];
              curve_nodes[21+1] = nodes[8];
              curve_nodes[21+2] = nodes[7];

            }
            else {
              ncurve = 1; ncurve_nodes = 9;

              curve_nodes[0] = nodes[0];
              curve_nodes[1] = nodes[1];
              curve_nodes[2] = nodes[2];
              curve_nodes[3] = nodes[5];
              curve_nodes[4] = nodes[8];
              curve_nodes[5] = nodes[7];
              curve_nodes[6] = nodes[6];
              curve_nodes[7] = nodes[3];
              curve_nodes[8] = nodes[0];

            }
          }
          else if ( name==-TET4 ) {
            ncurve = 4; ncurve_nodes = 4;

            curve_nodes[0+0] = nodes[0];
            curve_nodes[0+1] = nodes[1];
            curve_nodes[0+2] = nodes[2];
            curve_nodes[0+3] = nodes[0];

            curve_nodes[4+0] = nodes[0];
            curve_nodes[4+1] = nodes[1];
            curve_nodes[4+2] = nodes[3];
            curve_nodes[4+3] = nodes[0];

            curve_nodes[8+0] = nodes[1];
            curve_nodes[8+1] = nodes[2];
            curve_nodes[8+2] = nodes[3];
            curve_nodes[8+3] = nodes[1];

            curve_nodes[12+0] = nodes[2];
            curve_nodes[12+1] = nodes[3];
            curve_nodes[12+2] = nodes[0];
            curve_nodes[12+3] = nodes[2];

          }
          else if ( name==-HEX8 ) {
            ncurve = 6; ncurve_nodes = 5;

            curve_nodes[0+0] = nodes[0];
            curve_nodes[0+1] = nodes[1];
            curve_nodes[0+2] = nodes[3];
            curve_nodes[0+3] = nodes[2];
            curve_nodes[0+4] = nodes[0];

            curve_nodes[5+0] = nodes[4];
            curve_nodes[5+1] = nodes[5];
            curve_nodes[5+2] = nodes[7];
            curve_nodes[5+3] = nodes[6];
            curve_nodes[5+4] = nodes[4];

            curve_nodes[10+0] = nodes[0];
            curve_nodes[10+1] = nodes[1];
            curve_nodes[10+2] = nodes[5];
            curve_nodes[10+3] = nodes[4];
            curve_nodes[10+4] = nodes[0];

            curve_nodes[15+0] = nodes[1];
            curve_nodes[15+1] = nodes[3];
            curve_nodes[15+2] = nodes[7];
            curve_nodes[15+3] = nodes[5];
            curve_nodes[15+4] = nodes[1];

            curve_nodes[20+0] = nodes[2];
            curve_nodes[20+1] = nodes[3];
            curve_nodes[20+2] = nodes[7];
            curve_nodes[20+3] = nodes[6];
            curve_nodes[20+4] = nodes[2];

            curve_nodes[25+0] = nodes[0];
            curve_nodes[25+1] = nodes[2];
            curve_nodes[25+2] = nodes[6];
            curve_nodes[25+3] = nodes[4];
            curve_nodes[25+4] = nodes[0];

          }
          else {
            pri( "Error: CONTROL_PRINT_PLOTMTV not available for ", name );
            exit(TN_EXIT_STATUS );
          }      

          if ( ncurve==0 ) db_error( CONTROL_PRINT_PLOTMTV, icontrol );
          assert( ncurve<=MCURVE );
          for ( icurve=0; icurve<ncurve; icurve++ ) {
            for ( inol=0; inol<ncurve_nodes; inol++ ) {
              inod = curve_nodes[icurve*ncurve_nodes+inol];
              if ( inod==-1 ) {
                array_set( node, 0., ndim );
                for ( inol=0; inol<nnol; inol++ ) {
                  inod = nodes[inol];
                  db( NODE, inod, idum, work, ldum, VERSION_NORMAL, GET );
                  if ( ieigen>=0 ) {
                    db( NODE_EIGEN, inod, idum, node_eigen, ldum, VERSION_NORMAL, GET );
                    for ( idim=0; idim<ndim; idim++ )
                      node[idim] += node_eigen[ieigen*nuknwn+vel_indx+idim*nder];
                  }
                  array_add( node, work, node, ndim );
                }
                array_multiply( node, node, 1./nnol, ndim );
              }
              else {
                db( NODE, inod, idum, node, ldum, VERSION_NORMAL, GET );
                if ( ieigen>=0 ) {
                  db( NODE_EIGEN, inod, idum, node_eigen, ldum, VERSION_NORMAL, GET );
                  for ( idim=0; idim<ndim; idim++ )
                    node[idim] += node_eigen[ieigen*nuknwn+vel_indx+idim*nder];
                }
              }
              db( NODE_EIGEN, inod, idum, node_eigen, ldum, VERSION_NORMAL, GET_IF_EXISTS );
              if ( db_active_index( NODE_DOF, inod, VERSION_NORMAL ) )
                node_dof = db_dbl( NODE_DOF, inod, VERSION_NORMAL );
              if ( db_active_index( NODE_DOF_CALCUL, inod, VERSION_NORMAL ) )
                node_dof_calcul = db_dbl( NODE_DOF_CALCUL, inod, VERSION_NORMAL );
              for ( idim=0; idim<MDIM; idim++ ) {
                if     ( print_unknown || print_calcul_unknown || print_eigen_unknown ) {
                  if ( idim<=ndim-1 ) {
                    if ( materi_displacement )
                      out << node[idim]+node_dof[dis_indx+idim*nder] << " ";
                    else
                      out << node[idim] << " ";
                  }
                }
                else if ( print_mesh ) {
                  if ( idim>ndim-1 ) 
                    out << "0 ";
                  else
                    out << node[idim] << " ";
                }
                else if ( print_mesh_deformed ) {
                  if ( idim>ndim-1 ) 
                    out << "0 ";
                  else
                    out << node[idim]+node_dof[dis_indx+idim*nder] << " ";
                }
                else {
                  assert( print_displacement || print_velocity ||
                    print_calcul_vector );
                  if ( idim>ndim-1 ) 
                    out << "0 ";
                  else
                    out << node[idim] << " ";
                }
              }
              if      ( print_displacement ) {
                for ( idim=0; idim<MDIM; idim++ ) {
                  if ( idim>ndim-1 )
                    out << "0" << " ";
                  else
                    out << node_dof[dis_indx+idim*nder] << " ";
                }
              }
              else if ( print_velocity ) {
                for ( idim=0; idim<MDIM; idim++ ) {
                  if ( idim>ndim-1 )
                    out << "0" << " ";
                  else
                    out << node_dof[vel_indx+idim*nder] << " ";
                }
              }
              else if ( print_calcul_vector ) {
                for ( idim=0; idim<MDIM; idim++ )
                  out << node_dof_calcul[icalcul+idim] << " ";
              }
              else if ( print_unknown )
                out << node_dof[iuknwn] << " ";
              else if ( print_calcul_unknown )
                out << node_dof_calcul[icalcul] << " ";
              else if ( print_eigen_unknown )
                out << node_eigen[eigen_mode*nuknwn+iuknwn] << " ";
              out << "\n";
            }
            out << "\n";
          }
        }
      }
    }
  }

  out << "$ END\n";
  out << "$ QUIT\n";
  out.close();

  delete[] dof_principal;
  delete[] dof_label;
  delete[] el;
  delete[] nodes;
  delete[] curve_nodes;
  delete[] post_calcul_scal_vec_mat;
  delete[] control_eigen_values;
  delete[] node_eigen;

  if ( swit ) pri( "Out routine PRINT_PLOTMTV" );
}
