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

// an element is "empty" when its ELEMENT_EMPTY record is -YES (computed by
// the solver for materi_diffusion / materi_density); the print_g5 pattern
// counts an element as present when the record is -NO or -FRONT.
static long int vtk_element_is_empty( long int element )

{
  long int element_empty = -NO, ldum=0;
  double ddum[1];

  db( ELEMENT_EMPTY, element, &element_empty, ddum, ldum,
    VERSION_PRINT, GET_IF_EXISTS );
  return ( element_empty==-YES );
}

void print_vtk( long int icontrol )

{
  long int inod=0, element=0, max_node=0, max_element=0, nnol=0, 
    name=0, length=0, i=0, idim=0, jdim=0, ipuknwn=0, iuknwn=0, 
    icalcul=0, ready=0, indx=0, nval=0, length_cells=0,
    length_post_calcul_scal_vec_mat=0,
    calcul_unknown=0, calcul_operat=0, swit=0, ldum=0, 
    nvtk_dof=0, nvtk_dof_calcul=0, ifilter=0, print_field=1,
    vtk_coord=-YES, vtk_empty=-YES, vtk_node_method=-NODE_START_REFINED,
    vtk_other=-YES, ncell=0,
    idum[1], *dof_label=NULL, *dof_type=NULL, *dof_scal_vec_mat=NULL, 
    *post_calcul_scal_vec_mat=NULL, *post_calcul_unknown_operat=NULL, 
    *nodes=NULL, *el=NULL, *vtk_dof=NULL, *vtk_dof_calcul=NULL,
    *node_bounded=NULL, *print_post_field=NULL;
  double ddum[1], coord[MDIM], *node_dof=NULL, *node_dof_calcul=NULL;
  char str[MCHAR], outputname[MCHAR], filename[MCHAR];

  swit = set_swit(-1,-1,"print_vtk");
  if ( swit ) pri( "In routine PRINT_VTK" );

  // control_print_vtk_dof: only the listed solution fields are written
  // (e.g. -condif_temperature, -materi_velocity, ...); -none -> nothing.
  vtk_dof = get_new_int(DATA_ITEM_SIZE);
  db( CONTROL_PRINT_VTK_DOF, icontrol, vtk_dof, ddum, nvtk_dof,
    VERSION_NORMAL, GET_IF_EXISTS );

  // control_print_vtk_coord (6.340): -yes (default) writes the node
  // coordinates (the POINTS block); -no omits them.
  db( CONTROL_PRINT_VTK_COORD, icontrol, &vtk_coord, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  // control_print_vtk_empty (6.343): -yes (default) includes empty
  // elements; -no skips the elements whose ELEMENT_EMPTY record is -YES
  // (pattern print_g5.cc: element_empty==-NO || -FRONT -> not empty).
  db( CONTROL_PRINT_VTK_EMPTY, icontrol, &vtk_empty, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  // control_print_vtk_node_method (6.345): which node coordinates are
  // written - -node (stored coordinates), -node_start_refined (default:
  // NODE_START_REFINED when available, stored coordinates otherwise) or
  // -node_deformed_mesh (stored coordinates + nodal displacement).
  db( CONTROL_PRINT_VTK_NODE_METHOD, icontrol, &vtk_node_method, ddum,
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  // control_print_vtk_other (6.346): -yes (default) also writes the
  // "other" fields (boundary conditions, mesh deformation - see the
  // partial subset documented in manual-developer); -no omits them.
  db( CONTROL_PRINT_VTK_OTHER, icontrol, &vtk_other, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  // control_print_vtk_dof_calcul (6.342): only the listed post fields are
  // written; -none -> no post field; without the record all post fields.
  vtk_dof_calcul = get_new_int(DATA_ITEM_SIZE);
  db( CONTROL_PRINT_VTK_DOF_CALCUL, icontrol, vtk_dof_calcul, ddum,
    nvtk_dof_calcul, VERSION_NORMAL, GET_IF_EXISTS );

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

  // control_print_vtk_coord: -no omits the POINTS block (the only place
  // where coordinates appear in the vtk output). A file without POINTS is
  // not a valid visualization dataset; the switch exists for smaller debug
  // dumps (limitation documented in manual-developer).
  if ( vtk_coord!=-NO ) {
    outvtk << "POINTS " << max_node+1 << " double\n";
    for ( inod=0; inod<=max_node; inod++ ) {
      db( NODE, inod, idum, coord, ldum, VERSION_PRINT, GET );
      if ( vtk_node_method==-NODE_DEFORMED_MESH && materi_displacement ) {
        node_dof = db_dbl( NODE_DOF, inod, VERSION_PRINT );
        for ( idim=0; idim<MDIM; idim++ ) {
          if      ( idim>ndim-1 )
            outvtk << "0.0" << " ";
          else {
            ddum[0] = coord[idim]+node_dof[dis_indx+idim*nder];
            if ( ddum[0]==0.0 )
              outvtk << "0.0" << " ";
            else
              outvtk << ddum[0] << " ";
          }
        }
      }
      else if ( vtk_node_method==-NODE_START_REFINED &&
                db_active_index( NODE_START_REFINED, inod, VERSION_PRINT ) ) {
        db( NODE_START_REFINED, inod, idum, coord, ldum, VERSION_PRINT, GET );
        for ( idim=0; idim<MDIM; idim++ ) {
          if ( coord[idim]==0.0 )
            outvtk << "0.0" << " ";
          else
            outvtk << coord[idim] << " ";
        }
      }
      else {
        for ( idim=0; idim<MDIM; idim++ ) {
          if ( coord[idim]==0.0 )
            outvtk << "0.0" << " ";
          else
            outvtk << coord[idim] << " ";
        }
      }
      outvtk << "\n";
    } 
    outvtk << "\n";
  }

  // control_print_vtk_empty: -no excludes the empty elements
  // (ELEMENT_EMPTY == -YES) from CELLS and CELL_TYPES.
  length_cells = 0; ncell = 0;
  for ( element=0; element<=max_element; element++ ) {
    if ( vtk_empty==-NO && vtk_element_is_empty( element ) ) continue;
    ncell++;
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

  outvtk << "CELLS " << ncell << " " << length_cells << " \n";
  for ( element=0; element<=max_element; element++ ) {
    if ( vtk_empty==-NO && vtk_element_is_empty( element ) ) continue;
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

  outvtk << "CELL_TYPES " << ncell << "\n";
  for ( element=0; element<=max_element; element++ ) {
    if ( vtk_empty==-NO && vtk_element_is_empty( element ) ) continue;
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
      // filter by control_print_vtk_dof (dof_type names, or -none)
      print_field = 1;
      if ( nvtk_dof>0 && vtk_dof ) {
        if ( vtk_dof[0]==-NONE ) print_field = 0;
        else {
          print_field = 0;
          for ( ifilter=0; ifilter<nvtk_dof; ifilter++ )
            if ( dof_type[iuknwn]==vtk_dof[ifilter] ) { print_field = 1; break; }
        }
      }
      if      ( dof_scal_vec_mat[iuknwn]==-SCALAR )
        nval = 1;
      else if ( dof_scal_vec_mat[iuknwn]==-VECTOR )
        nval = ndim;
      else {
        assert( dof_scal_vec_mat[iuknwn]==-MATRIX );
        nval = 6;
      }
      if ( print_field ) {
        if      ( dof_scal_vec_mat[iuknwn]==-SCALAR ) {
          outvtk << "SCALARS " << db_name(dof_label[iuknwn]) << " double\n";
          outvtk << "LOOKUP_TABLE default\n";
        }
        else if ( dof_scal_vec_mat[iuknwn]==-VECTOR ) {
          outvtk << "VECTORS " << db_name(dof_type[iuknwn]) << " double\n";
        }
        else {
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
      }
      ipuknwn += nval;
      ready = (ipuknwn>=npuknwn);
    }

      // write vector components and tensors components for primary unknowns
    ipuknwn = 0; ready = 0;
    while ( !ready ) {
      iuknwn = ipuknwn*nder;
      print_field = 1;
      if ( nvtk_dof>0 && vtk_dof ) {
        if ( vtk_dof[0]==-NONE ) print_field = 0;
        else {
          print_field = 0;
          for ( ifilter=0; ifilter<nvtk_dof; ifilter++ )
            if ( dof_type[iuknwn]==vtk_dof[ifilter] ) { print_field = 1; break; }
        }
      }
      if ( print_field && dof_scal_vec_mat[iuknwn]!=-SCALAR ) {
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

    // derived magnitudes from the first stress tensor dof (von Mises,
    // Tresca, principal stresses), reusing calc_derived() from derived.cc
    {
      long int sig_indx = -1;
      for ( i=0; i<nuknwn; i++ )
        if ( dof_scal_vec_mat[i]==-MATRIX ) { sig_indx = i; break; }
      // respect control_print_vtk_dof: only write derived magnitudes when
      // the filter (if any) includes a stress field
      long int print_derived = ( nvtk_dof<=0 || !vtk_dof );
      if ( !print_derived && vtk_dof ) {
        for ( ifilter=0; ifilter<nvtk_dof; ifilter++ )
          if ( vtk_dof[ifilter]==-MATERI_STRESS ) { print_derived = 1; break; }
      }
      if ( sig_indx>=0 && print_derived ) {
        const char* derived_names[5] =
          { "vmises", "tresca", "sig1", "sig2", "sig3" };
        for ( idim=0; idim<5; idim++ ) {
          outvtk << "SCALARS " << derived_names[idim] << " double\n";
          outvtk << "LOOKUP_TABLE default\n";
          for ( inod=0; inod<=max_node; inod++ ) {
            node_dof = db_dbl( NODE_DOF, inod, VERSION_PRINT );
            double sig[6], dout[5];
            sig[0]=node_dof[sig_indx+stress_indx(0,0)*nder];
            sig[1]=node_dof[sig_indx+stress_indx(1,1)*nder];
            sig[2]=node_dof[sig_indx+stress_indx(2,2)*nder];
            sig[3]=node_dof[sig_indx+stress_indx(0,1)*nder];
            sig[4]=node_dof[sig_indx+stress_indx(0,2)*nder];
            sig[5]=node_dof[sig_indx+stress_indx(1,2)*nder];
            calc_derived( sig, dout );
            outvtk << dout[idim] << "\n";
          }
          outvtk << "\n";
        }
      }
    }

    if ( db_active_index( POST_CALCUL, 0,  VERSION_NORMAL ) ) {
      db( POST_CALCUL_UNKNOWN_OPERAT, 0, post_calcul_unknown_operat, 
        ddum, ldum, VERSION_NORMAL, GET );
      db( POST_CALCUL_SCAL_VEC_MAT, 0, post_calcul_scal_vec_mat, 
        ddum, length_post_calcul_scal_vec_mat, VERSION_NORMAL, GET );

      // control_print_vtk_dof_calcul: one flag per post field. Without the
      // record every field is written (default). -none -> no post field.
      // Otherwise a field is written when any listed name matches its
      // underlying unknown exactly (db_name of the initia name, e.g.
      // -materi_stress matches every operator of the stress tensor) or is
      // a substring of its label (post_calcul_names[icalcul], e.g. -sigyy
      // matches sigyy and tosigyy). GNU adaptation: post_calcul_label does
      // not exist in the GNU, the labels are post_calcul_names.
      print_post_field = get_new_int( length_post_calcul_scal_vec_mat );
      for ( icalcul=0; icalcul<length_post_calcul_scal_vec_mat; icalcul++ )
        print_post_field[icalcul] = 1;
      if ( nvtk_dof_calcul>0 && vtk_dof_calcul ) {
        if ( vtk_dof_calcul[0]==-NONE ) {
          for ( icalcul=0; icalcul<length_post_calcul_scal_vec_mat; icalcul++ )
            print_post_field[icalcul] = 0;
        }
        else {
          for ( icalcul=0; icalcul<length_post_calcul_scal_vec_mat; icalcul++ ) {
            print_post_field[icalcul] = 0;
            for ( ifilter=0; ifilter<nvtk_dof_calcul; ifilter++ ) {
              if ( post_calcul_unknown_operat[icalcul*2+0]==
                   vtk_dof_calcul[ifilter] ) {
                print_post_field[icalcul] = 1;
                break;
              }
              if ( strstr( post_calcul_names[icalcul],
                           db_name(vtk_dof_calcul[ifilter]) ) ) {
                print_post_field[icalcul] = 1;
                break;
              }
            }
          }
        }
      }

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
          if ( print_post_field[icalcul] ) {
            outvtk << "SCALARS " << outputname << " double\n";
            outvtk << "LOOKUP_TABLE default\n";
          }
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
          if ( print_post_field[icalcul] ) {
            outvtk << "VECTORS " << outputname << " double\n";
          }
        }
        if ( print_post_field[icalcul] ) {
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
        }
        icalcul += nval;
        ready = (icalcul>=length_post_calcul_scal_vec_mat);
      }

        // write vector components and tensor components for calculated data
      icalcul = 0; ready=0;
      while ( !ready ) {
        if ( post_calcul_scal_vec_mat[icalcul]!=-SCALAR &&
             print_post_field[icalcul] ) {
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

    // control_print_vtk_other (6.346, default -yes): "other things" like
    // boundary conditions and mesh deformation. Partial subset (the
    // Professional manual does not detail the list): 1) boundary_condition
    // scalar = 1 when any primary dof of the node is bounded
    // (node_bounded), 0 otherwise; 2) mesh_deformation vector = the nodal
    // displacement, written only when materi_displacement is active.
    // node_bounded has no version_all -> read from VERSION_NORMAL, where
    // node indices are 1-based: the print index inod (0-based VERSION_PRINT
    // after renumbering with lowest_node=0) maps to inod+1 as long as no
    // node is deleted (limitation documented in manual-developer).
    if ( vtk_other!=-NO ) {
      long int node_bounded_max = -1;
      db_max_index( NODE_BOUNDED, node_bounded_max, VERSION_NORMAL, GET );
      if ( node_bounded_max>=0 ) {
        outvtk << "SCALARS boundary_condition double\n";
        outvtk << "LOOKUP_TABLE default\n";
        for ( inod=0; inod<=max_node; inod++ ) {
          long int bounded = 0;
          if ( db_active_index( NODE_BOUNDED, inod+1, VERSION_NORMAL ) ) {
            node_bounded = db_int( NODE_BOUNDED, inod+1, VERSION_NORMAL );
            for ( i=0; i<npuknwn && !bounded; i++ )
              if ( node_bounded[i] ) bounded = 1;
          }
          if ( bounded )
            outvtk << "1.0" << "\n";
          else
            outvtk << "0.0" << "\n";
        }
        outvtk << "\n";
      }
      if ( materi_displacement ) {
        outvtk << "VECTORS mesh_deformation double\n";
        for ( inod=0; inod<=max_node; inod++ ) {
          node_dof = db_dbl( NODE_DOF, inod, VERSION_PRINT );
          for ( idim=0; idim<MDIM; idim++ ) {
            if      ( idim>ndim-1 )
              outvtk << "0.0" << " ";
            else {
              ddum[0] = node_dof[dis_indx+idim*nder];
              if ( ddum[0]==0.0 )
                outvtk << "0.0" << " ";
              else
                outvtk << ddum[0] << " ";
            }
          }
          outvtk << "\n";
        }
        outvtk << "\n";
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
  delete[] vtk_dof;
  delete[] vtk_dof_calcul;
  if ( print_post_field ) delete[] print_post_field;

  if ( swit ) pri( "Out routine PRINT_VTK" );
}
