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

void print_vtk( long int icontrol )

{
  long int inod=0, element=0, max_node=0, max_element=0, nnol=0, 
    name=0, length=0, i=0, idim=0, jdim=0, ipuknwn=0, iuknwn=0, 
    icalcul=0, ready=0, indx=0, nval=0, length_cells=0,
    length_post_calcul_scal_vec_mat=0,
    calcul_unknown=0, calcul_operat=0, swit=0, ldum=0, 
    idum[1], *dof_label=NULL, *dof_type=NULL, *dof_scal_vec_mat=NULL, 
    *post_calcul_scal_vec_mat=NULL, *post_calcul_unknown_operat=NULL, 
    *nodes=NULL, *el=NULL;
  double ddum[1], coord[MDIM], *node_dof=NULL, *node_dof_calcul=NULL;
  char str[MCHAR], outputname[MCHAR], filename[MCHAR];

  swit = set_swit(-1,-1,"print_vtk");
  if ( swit ) pri( "In routine PRINT_VTK" );

  dof_label = get_new_int(MUKNWN);
  dof_type = get_new_int(MUKNWN);
  dof_scal_vec_mat = get_new_int(MUKNWN);
  post_calcul_scal_vec_mat = get_new_int(DATA_ITEM_SIZE);
  post_calcul_unknown_operat = get_new_int(DATA_ITEM_SIZE);
  nodes = get_new_int(MAXIMUM_NODE);
  el = get_new_int(MAXIMUM_NODE+1);

  db_version_copy( VERSION_NORMAL, VERSION_PRINT );
  renumbering( VERSION_PRINT, NO, 0, 0, idum, idum );
  db_highest_index( NODE, max_node, VERSION_PRINT );
  db_highest_index( ELEMENT, max_element, VERSION_PRINT );
  if ( max_element<0 ) return;

  strcpy( filename, "tn" );
  if ( icontrol>=0 ) {
    long_to_a( icontrol, str );
    strcat( filename, str );
  }
  strcat( filename, ".vtk" );
  ofstream outvtk( filename );
  outvtk.precision(TN_PRECISION);

  outvtk << "# vtk DataFile Version 2.0\n";
  outvtk << "Calculation " << data_file_base << "\n";
  outvtk << "ASCII\n\n";

  outvtk << "DATASET UNSTRUCTURED_GRID\n\n";

  outvtk << "POINTS " << max_node+1 << " double\n";
  for ( inod=0; inod<=max_node; inod++ ) {
    db( NODE, inod, idum, coord, ldum, VERSION_PRINT, GET );
    for ( idim=0; idim<MDIM; idim++ ) {
      if      ( idim>ndim-1 )
        outvtk << "0.0" << " ";
      else if ( materi_displacement ) {
        node_dof = db_dbl( NODE_DOF, inod, VERSION_PRINT );
        if ( coord[idim]+node_dof[dis_indx+idim*nder]==0.0 )
          outvtk << "0.0" << " ";
        else
          outvtk << coord[idim]+node_dof[dis_indx+idim*nder] << " ";
      }
      else {
        if ( coord[idim]==0.0 )
          outvtk << "0.0" << " ";
        else
          outvtk << coord[idim] << " ";
      }
    }
    outvtk << "\n";
  } 
  outvtk << "\n";

  length_cells = 0;
  for ( element=0; element<=max_element; element++ ) {
    db( ELEMENT, element, el, ddum, length, VERSION_PRINT, GET );
    name = el[0];
    if      ( name==-BAR2 )  length_cells += 3;
    else if ( name==-BAR3 )  length_cells += 3;
    else if ( name==-TRIA3 ) length_cells += 4;
    else if ( name==-QUAD4 ) length_cells += 5;
    else if ( name==-QUAD9 ) length_cells += 5;
    else if ( name==-TET10)  length_cells += 5;
    else if ( name==-HEX8 )  length_cells += 9;
    else if ( name==-HEX27 ) length_cells += 9;
    else {
      pri( "Error: illegal element type detected for CONTROL_PRINT_VTK.\n");
      exit(TN_EXIT_STATUS);
    }
  }

  outvtk << "CELLS " << max_element+1 << " " << length_cells << " \n";
  for ( element=0; element<=max_element; element++ ) {
    db( ELEMENT, element, el, ddum, length, VERSION_PRINT, GET );
    name = el[0];
    nnol = length - 1; array_move( &el[1], nodes, nnol );
    if      ( name==-BAR2 ) {
      outvtk << "2 ";
      outvtk << nodes[0] << " " << nodes[1];
    }
    else if ( name==-BAR3 ) {
      outvtk << "2 ";
      outvtk << nodes[0] << " " << nodes[2];
    }
    else if ( name==-TRIA3 ) {
      outvtk << "3 ";
      outvtk << nodes[0] << " " << nodes[1] << " " << nodes[2];
    }
    else if ( name==-TRIA6 ) {
      outvtk << "3 ";
      outvtk << nodes[0] << " " << nodes[2] << " " << nodes[5];
    }
    else if ( name==-QUAD4 ) {
      outvtk << "4 ";
      outvtk << nodes[0] << " " << nodes[1] << " " << nodes[3] << " " << nodes[2];
    }
    else if ( name==-QUAD9 ) {
      outvtk << "4 ";
      outvtk << nodes[0] << " " << nodes[2] << " " << nodes[8] << " " << nodes[6];
    }
    else if ( name==-TET4 ) {
      outvtk << "4 ";
      outvtk << nodes[0] << " " << nodes[1] << " " << nodes[2] << " " << nodes[3];
    }
    else if ( name==-TET10 ) {
      outvtk << "4 ";
      outvtk << nodes[0] << " " << nodes[2] << " " << nodes[5] << " " << nodes[9];
    }
    else if ( name==-HEX8 ) {
      outvtk << "8 ";
      outvtk << nodes[0] << " " << nodes[1] << " " << nodes[3] << " " << nodes[2] << " " <<
             nodes[4] << " " << nodes[5] << " " << nodes[7] << " " << nodes[6];
    }
    else {
      assert( name==-HEX27 );
      outvtk << "8 ";
      outvtk << nodes[0] << " " << nodes[2] << " " << nodes[8] << " " << nodes[6] << " " <<
             nodes[18] << " " << nodes[20] << " " << nodes[26] << " " << nodes[24];
    }
    outvtk << "\n";
  }
  outvtk << "\n";

  outvtk << "CELL_TYPES " << max_element+1 << "\n";
  for ( element=0; element<=max_element; element++ ) {
    db( ELEMENT, element, el, ddum, length, VERSION_PRINT, GET );
    name = el[0];
    if      ( name==-BAR2 )  outvtk << "3 ";
    else if ( name==-BAR3 )  outvtk << "3 ";
    else if ( name==-TRIA3 ) outvtk << "5 ";
    else if ( name==-TRIA6 ) outvtk << "5 ";
    else if ( name==-QUAD4 ) outvtk << "9 ";
    else if ( name==-QUAD9 ) outvtk << "9 ";
    else if ( name==-TET4  ) outvtk << "10 ";
    else if ( name==-TET10 ) outvtk << "10 ";
    else if ( name==-HEX8 )  outvtk << "12 ";
    else if ( name==-HEX27 ) outvtk << "12 ";
    else {
      pri( "Error: illegal element type detected for CONTROL_PRINT_VTK.\n");
      exit(TN_EXIT_STATUS);
    }
    outvtk << "\n";
  }
  outvtk << "\n";

  if ( npuknwn>0 ) {
    db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET );
    db( DOF_TYPE, 0, dof_type, ddum, ldum, VERSION_NORMAL, GET );
    db( DOF_SCAL_VEC_MAT, 0, dof_scal_vec_mat, ddum, ldum, VERSION_NORMAL, GET );
    outvtk << "POINT_DATA " << max_node+1 << "\n\n";

      // write scalars, vectors and tensors for primary unknowns
    ipuknwn = 0; ready = 0;
    while ( !ready ) {
      iuknwn = ipuknwn*nder;
      if      ( dof_scal_vec_mat[iuknwn]==-SCALAR ) {
        nval = 1;
        outvtk << "SCALARS " << db_name(dof_label[iuknwn]) << " double\n";
        outvtk << "LOOKUP_TABLE default\n";
      }
      else if ( dof_scal_vec_mat[iuknwn]==-VECTOR ) {
        nval = ndim;
        outvtk << "VECTORS " << db_name(dof_type[iuknwn]) << " double\n";
      }
      else {
        assert( dof_scal_vec_mat[iuknwn]==-MATRIX );
        nval = 6;
        outvtk << "TENSORS " << db_name(dof_type[iuknwn]) << " double\n";
      }
      for ( inod=0; inod<=max_node; inod++ ) {
        node_dof = db_dbl( NODE_DOF, inod, VERSION_PRINT );
        if      ( dof_scal_vec_mat[iuknwn]==-SCALAR ) {
          if ( node_dof[iuknwn]==0.0 )
            outvtk << "0.0";
          else
            outvtk << node_dof[iuknwn];
        }
        else if ( dof_scal_vec_mat[iuknwn]==-VECTOR ) {
          for ( idim=0; idim<MDIM; idim++ ) {
            if      ( idim>ndim-1 )
              outvtk << "0.0" << " ";
            else {
              indx = iuknwn+idim*nder;
              if ( node_dof[indx]==0.0 )
                outvtk << "0.0" << " ";
              else
                outvtk << node_dof[indx] << " ";
            }
          }
        }
        else {
          assert( dof_scal_vec_mat[iuknwn]==-MATRIX );
          for ( idim=0; idim<MDIM; idim++ ) {
            for ( jdim=0; jdim<MDIM; jdim++ ) {
              indx = iuknwn + stress_indx(idim,jdim)*nder;
              if ( node_dof[indx]==0.0 )
                outvtk << "0.0" << " ";
              else
                outvtk << node_dof[indx] << " ";
            }
            if ( idim!=MDIM-1 ) outvtk << "\n";
          }
        }
        outvtk << "\n";
      }
      outvtk << "\n";
      ipuknwn += nval;
      ready = (ipuknwn>=npuknwn);
    }

      // write vector components and tensors components for primary unknowns
    ipuknwn = 0; ready = 0;
    while ( !ready ) {
      iuknwn = ipuknwn*nder;
      if ( dof_scal_vec_mat[iuknwn]!=-SCALAR ) {
        outvtk << "SCALARS " << db_name(dof_label[iuknwn]) << " double\n";
        outvtk << "LOOKUP_TABLE default\n";
        for ( inod=0; inod<=max_node; inod++ ) {
          node_dof = db_dbl( NODE_DOF, inod, VERSION_PRINT );
          if ( node_dof[iuknwn]==0.0 )
            outvtk << "0.0";
          else
          outvtk << node_dof[iuknwn];
          outvtk << "\n";
        }
        outvtk << "\n";
      }
      ipuknwn += 1;
      ready = (ipuknwn>=npuknwn);
    }

    if ( db_active_index( POST_CALCUL, 0,  VERSION_NORMAL ) ) {
      db( POST_CALCUL_UNKNOWN_OPERAT, 0, post_calcul_unknown_operat, 
        ddum, ldum, VERSION_NORMAL, GET );
      db( POST_CALCUL_SCAL_VEC_MAT, 0, post_calcul_scal_vec_mat, 
        ddum, length_post_calcul_scal_vec_mat, VERSION_NORMAL, GET );

        // write scalars, vectors and tensors for calculated data
      icalcul = idim = ready = 0;
      while ( !ready ) {
        calcul_unknown = post_calcul_unknown_operat[icalcul*2+0];
        calcul_operat = post_calcul_unknown_operat[icalcul*2+1];
        strcpy( outputname, db_name(calcul_unknown) );
        strcat( outputname, "_" );
        strcat( outputname, db_name(calcul_operat) );
        if      ( post_calcul_scal_vec_mat[icalcul]==-SCALAR ) {
          nval = 1;
          outvtk << "SCALARS " << outputname << " double\n";
          outvtk << "LOOKUP_TABLE default\n";
        }
        else {
          assert( post_calcul_scal_vec_mat[icalcul]==-VECTOR );
          nval = MDIM;
          strcat( outputname, "_" );
          long_to_a( idim, str); 
          strcat( outputname, str );
          idim++;
          if ( idim==MDIM ) {
            idim = 0;
          }
          outvtk << "VECTORS " << outputname << " double\n";
        }
        for ( inod=0; inod<=max_node; inod++ ) {
          node_dof_calcul = db_dbl( NODE_DOF_CALCUL, inod, VERSION_PRINT );
          if      ( post_calcul_scal_vec_mat[icalcul]==-SCALAR ) {
            if ( node_dof_calcul[icalcul]==0.0 )
              outvtk << "0.0";
            else
              outvtk << node_dof_calcul[icalcul];
          }
          else {
            assert( post_calcul_scal_vec_mat[icalcul]==-VECTOR );
            for ( i=0; i<MDIM; i++ ) {
              indx = icalcul+i;
              if ( node_dof_calcul[indx]==0.0 )
                outvtk << "0.0" << " ";
              else
                outvtk << node_dof_calcul[indx] << " ";
            }
          }
          outvtk  << "\n";
        }
        outvtk << "\n";
        icalcul += nval;
        ready = (icalcul>=length_post_calcul_scal_vec_mat);
      }

        // write vector components and tensor components for calculated data
      icalcul = 0; ready=0;
      while ( !ready ) {
        if ( post_calcul_scal_vec_mat[icalcul]!=-SCALAR ) {
          outvtk << "SCALARS " << post_calcul_names[icalcul] << " double\n";
          outvtk << "LOOKUP_TABLE default\n";
          for ( inod=0; inod<=max_node; inod++ ) {
            node_dof_calcul = db_dbl( NODE_DOF_CALCUL, inod, VERSION_PRINT );
            if ( node_dof_calcul[icalcul]==0.0 )
              outvtk << "0.0";
            else
              outvtk << node_dof_calcul[icalcul];
            outvtk  << "\n";
          }
          outvtk << "\n";
        }
        icalcul++;
        ready = (icalcul>=length_post_calcul_scal_vec_mat);
      }

    }
    outvtk << "\n";

  }

  outvtk.close();

  db_version_delete( VERSION_PRINT );

  delete[] dof_label;
  delete[] dof_type;
  delete[] dof_scal_vec_mat;
  delete[] post_calcul_scal_vec_mat;
  delete[] post_calcul_unknown_operat;
  delete[] nodes;
  delete[] el;

  if ( swit ) pri( "Out routine PRINT_VTK" );
}
