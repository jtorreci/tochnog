/*
    copyright (c) 2000  dennis roddeman
    email: d.g.roddeman@wb.utwente.nl

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

#define MAX_ELEMENT_TYPE 100
#define MAX_FIELD 1000
#define EPS_SMALL 1.e-12

void print_dx( long int final_call )

{
  long int n=0, indx=0, nnol=0, inod=0, max_node=0, idim=0, jdim=0, 
    nval=0, ready=0, icalcul=0, element=0, max_element=0, 
    swit=0, itype=0, ntype=0, ldum=0, length=0, name=0, 
    element_name=0, ipuknwn=0, iuknwn=0, old_name=0,
    rank=0, ifield=0, nfield=0, length_post_calcul_scal_vec_mat=0, 
    itime=0, length_control_print_dx_time=0,
    element_group=0, length_max=0, idum[1], type_nnol[MAX_ELEMENT_TYPE], 
    type_group[MAX_ELEMENT_TYPE], type_name[MAX_ELEMENT_TYPE], 
    type_max[MAX_ELEMENT_TYPE], el[1+MAXIMUM_NODE], nodes[MAXIMUM_NODE], 
    dof_label[MUKNWN], dof_amount[MUKNWN], 
    dof_type[MUKNWN], dof_scal_vec_mat[MUKNWN],
    post_calcul_scal_vec_mat[DATA_ITEM_SIZE],
    *old_node_numbers=NULL, *old_element_numbers=NULL;
 double tmp=0., time_current=0., ddum[MDIM], average_coord[MDIM],
    coord[MDIM], *node_dof=NULL, 
    *node_dof_calcul=NULL, *control_print_dx_time=NULL;
  char dof_type_name[MCHAR], outname[MCHAR], filename[MCHAR],
    field_names[MCHAR][MAX_FIELD];

    // initialize
  db( TIME_CURRENT, 0, idum, &time_current, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( DOF_AMOUNT, 0, dof_amount, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( DOF_TYPE, 0, dof_type, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( DOF_SCAL_VEC_MAT, 0, dof_scal_vec_mat, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( POST_CALCUL_SCAL_VEC_MAT, 0, post_calcul_scal_vec_mat, 
    ddum, length_post_calcul_scal_vec_mat, VERSION_NORMAL, GET_IF_EXISTS );

  db_highest_index( NODE, max_node, VERSION_NORMAL );
  db_highest_index( ELEMENT, max_element, VERSION_NORMAL );
  if ( max_node<0 || max_element<0 ) return;

    // allocate
  old_node_numbers = get_new_int(1+max_node);
  old_element_numbers = get_new_int(1+max_element);

    // create print version
  element_middle_radius_set();
  db_version_copy( VERSION_NORMAL, VERSION_PRINT );
  renumbering( VERSION_PRINT, YES, 0, 0, old_node_numbers, old_element_numbers );
  db_highest_index( NODE, max_node, VERSION_PRINT );
  db_highest_index( ELEMENT, max_element, VERSION_PRINT );

    // empty the dx file at the start of calculation
  strcpy( filename, data_file_base );
  strcat( filename, ".dx" );           
  if ( !db_active_index( CONTROL_PRINT_DX_TIME, 0, VERSION_NORMAL ) ) {
    ofstream outdx( filename );
    outdx.close();
  }

  swit = set_swit(-1,-1,"print_dx");
  if ( swit ) pri( "In routine PRINT_DX" );

    // open dx file
  ofstream outdx( filename, ios::app );
  outdx.precision(TN_PRECISION);
  outdx.setf(ios::showpoint);

    // allocate space
  if ( db_active_index( CONTROL_PRINT_DX_TIME, 0, VERSION_NORMAL ) )
    length_control_print_dx_time = db_len( CONTROL_PRINT_DX_TIME, 0, VERSION_NORMAL );
  control_print_dx_time = get_new_dbl(length_control_print_dx_time+1);

    // get previous time points
  if ( db_active_index( CONTROL_PRINT_DX_TIME, 0, VERSION_NORMAL ) )
    db( CONTROL_PRINT_DX_TIME, 0, idum, control_print_dx_time, 
      length_control_print_dx_time, VERSION_NORMAL, GET );

    // increase length of time points record if necessary
  length_max = db_data_length( CONTROL_PRINT_DX_TIME );
  if ( length_control_print_dx_time+1>length_max ) {
    length_max += 100;
    db_data_length_put( CONTROL_PRINT_DX_TIME, length_max );
    db( CONTROL_PRINT_DX_TIME, 0, idum, control_print_dx_time, 
      length_control_print_dx_time, VERSION_NORMAL, PUT );
  }

    // determine the element types in the mesh
  array_set( type_group, 0, MAX_ELEMENT_TYPE );
  array_set( type_name, 0, MAX_ELEMENT_TYPE );
  array_set( type_nnol, 0, MAX_ELEMENT_TYPE );
  array_set( type_max, -1, MAX_ELEMENT_TYPE );
  for ( element=0; element<=max_element; element++ ) {
    db( ELEMENT, element, el, ddum, ldum, VERSION_PRINT, GET );
    name = el[0];
    if ( element==0 ) old_name = name;
    element_group = 0; db( ELEMENT_GROUP, element, &element_group, ddum, 
      ldum, VERSION_PRINT, GET_IF_EXISTS );
    if ( !array_member(type_group,element_group,ntype,itype) ) {
      assert(ntype<MAX_ELEMENT_TYPE-1);
      nnol = 0;
      if      ( name==-BAR2   )        nnol = 2;
      else if ( name==-BAR3   )        nnol = 2;
      else if ( name==-BAR4   )        nnol = 2;
      else if ( name==-BEAM   )        nnol = 2;
      else if ( name==-CONTACTSPRING ) nnol = 2;
      else if ( name==-TRIA3  )        nnol = 3;
      else if ( name==-TRIA6  )        nnol = 3;
      else if ( name==-QUAD4  )        nnol = 4;
      else if ( name==-QUAD9  )        nnol = 4;
      else if ( name==-QUAD16 )        nnol = 4;
      else if ( name==-TET4   )        nnol = 4;
      else if ( name==-TET10  )        nnol = 4;
      else if ( name==-HEX8   )        nnol = 8;
      else if ( name==-HEX27  )        nnol = 8;
      else if ( name==-HEX64  )        nnol = 8;
      else if ( name==-SPRING1 )       nnol = 1;
      else if ( name==-SPRING2 )       nnol = 2;
      else if ( name==-TRUSS  )        nnol = 2;
      else if ( name==-TRUSSBEAM )     nnol = 2;
      if ( nnol>0 ) {
        itype = ntype;
        type_group[itype] = element_group;
        type_name[itype] = name;
        type_nnol[itype] = nnol;
        type_max[itype]++;
        ntype++;
      }
    }
    else {
      type_max[itype]++;
    }
  }

    // field names for scalar unknowns
  for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
    iuknwn = ipuknwn*nder;
    strcpy( dof_type_name, db_name(dof_label[iuknwn]) );
    assert( nfield<MAX_FIELD-1 );
    strcpy( field_names[nfield], dof_type_name ); nfield++;
  }
    // field names for vector and matrix unknowns
  ipuknwn = 0;
  if ( nuknwn==0 )
    ready = 1;
  else
    ready = 0;
  while ( !ready ) {
    iuknwn = ipuknwn*nder;
    if      ( dof_scal_vec_mat[iuknwn]==-SCALAR ) {
      rank = 0;
      nval = dof_amount[iuknwn];
    }
    else if ( dof_scal_vec_mat[iuknwn]==-VECTOR ) {
      rank = 1;
      nval = dof_amount[iuknwn];
    }
    else {
      assert( dof_scal_vec_mat[iuknwn]==-MATRIX );
      rank = 1;
      nval = dof_amount[iuknwn];
    }                                  
    if ( rank>0 ) {
      strcpy( dof_type_name, db_name(dof_type[iuknwn]) );
      assert( nfield<MAX_FIELD-1 );
      strcpy( field_names[nfield], dof_type_name ); nfield++;
    }
    ipuknwn += nval;
    if ( ipuknwn>=npuknwn ) ready = 1;
  }
    // field names for post_calcul scalars
  if ( db_active_index( POST_CALCUL, 0, VERSION_NORMAL ) ) {
      // field names for post_calcul scalars
    for ( icalcul=0; icalcul<length_post_calcul_scal_vec_mat; icalcul++ ) {
      if ( post_calcul_scal_vec_mat[icalcul]==-SCALAR ) {
        strcpy( outname, post_calcul_names_without_extension[icalcul] );
        assert( nfield<MAX_FIELD-1 );
        strcpy( field_names[nfield], outname ); nfield++;
      }
    }
      // field names for post_calcul vectors and matrices
    icalcul = 0;
    if ( length_post_calcul_scal_vec_mat<=0 )
      ready = 1;
    else
      ready = 0;
    while ( !ready ) {
      if      ( post_calcul_scal_vec_mat[icalcul]==-SCALAR ) {
        rank = 0;
        nval = 1;
      }
      else if ( post_calcul_scal_vec_mat[icalcul]==-VECTOR ) {
        rank = 1;
        nval = MDIM;
      }
      else {
        assert( post_calcul_scal_vec_mat[icalcul]==-MATRIX );
        rank = 2;
        nval = 6;
      }
      if ( rank>0 ) {
        strcpy( outname, post_calcul_names_without_extension[icalcul] );
        assert( nfield<MAX_FIELD-1 );
        strcpy( field_names[nfield], outname ); nfield++;
      }
      icalcul += nval;
      ready = (icalcul>=length_post_calcul_scal_vec_mat);
    }
  }

    // check if time_current newer than the last time the dx file was written
  if ( db_active_index( CONTROL_PRINT_DX_TIME, 0, VERSION_NORMAL ) &&
       time_current==control_print_dx_time[length_control_print_dx_time-1] ) 
    goto close_dx_file;
  else {
    control_print_dx_time[length_control_print_dx_time] = time_current;
    length_control_print_dx_time++;
    db( CONTROL_PRINT_DX_TIME, 0, idum, control_print_dx_time, 
      length_control_print_dx_time, VERSION_NORMAL, PUT );
  }

    // print nodes coordinates
  outdx << "#\n";
  outdx << "object \"nodes time=" << time_current << "\" ";
  outdx << "class array type double rank 1 shape "
    << ndim << " items " << 1+max_node << " data follows\n";
  for ( inod=0; inod<=max_node; inod++ ) {
    db( NODE, inod, idum, coord, ldum, VERSION_PRINT, GET );
    if ( materi_displacement )
      node_dof = db_dbl( NODE_DOF, inod, VERSION_PRINT );
    for ( idim=0; idim<ndim; idim++ ) {
      if ( materi_displacement )
        tmp = coord[idim] + node_dof[dis_indx+idim*nder];
      else
        tmp = coord[idim];
      if ( scalar_dabs(tmp)<EPS_SMALL ) tmp = 0.;
      outdx << tmp << " ";
    }
    outdx << "\n";
  }

    // print nodes numbers
  outdx << "#\n";
  outdx << "object \"node_numbers"
    << " time=" << time_current << "\" ";
  outdx << "class array type string rank 1 shape 5"
    " items " << 1+max_node << " data follows\n";
  for ( inod=0; inod<=max_node; inod++ ) {
    outdx << "\"" << old_node_numbers[inod] << "\"";
    outdx << "\n";
  }
  outdx << "attribute \"dep\" string \"positions\"\n";   
  outdx << "#\n";
  outdx << "object \"" << "node_numbers" << " field"
    << " time=" << time_current << "\" ";
  outdx << "class field\n";
  outdx << "component \"positions\" value \"nodes" 
    << " time=" << time_current << "\"\n";
  outdx << "component \"data\" value \"" << "node_numbers" 
    << " time=" << time_current << "\"\n";

    // print elements for each type as a separate object
  for ( itype=0; itype<ntype; itype++ ) {
    name = type_name[itype];
    nnol = type_nnol[itype];
    outdx << "#\n";
    outdx << "object \"elements group " << type_group[itype] << " time=" << time_current << "\" ";;
    outdx  << "class array type int rank 1 shape "
      << nnol << " items " << 1+type_max[itype] << " data follows\n";
    for ( element=0; element<=max_element; element++ ) {
      element_group = 0;
      db( ELEMENT_GROUP, element, &element_group, ddum, ldum, 
        VERSION_PRINT, GET_IF_EXISTS );
      if ( element_group==type_group[itype] ) {
        db( ELEMENT, element, el, ddum, length, VERSION_PRINT, GET );
        element_name = el[0];
        if ( element_name!=name ) {
          pri( "Warning: the elements should be all the same in an group for dx files.");
          pri( "No valid dx files generated.");
          goto end_of_dx;
        }
        array_move( &el[1], nodes, length-1 );
        if      ( name==-BAR2 )
          outdx << nodes[0] << " " << nodes[1];
        else if ( name==-BAR3 )
          outdx << nodes[0] << " " << nodes[2];
        else if ( name==-BAR4 )
          outdx << nodes[0] << " " << nodes[3];
        else if ( name==-BEAM )
          outdx << nodes[0] << " " << nodes[1];
        else if ( name==-CONTACTSPRING )
          outdx << nodes[0] << " " << nodes[1];
        else if ( name==-TRIA3 )
          outdx << nodes[0] << " " << nodes[1] << " " 
                << nodes[2];
        else if ( name==-TRIA6 )
          outdx << nodes[0] << " " << nodes[2] << " " 
                << nodes[5];
        else if ( name==-QUAD4 )
          outdx << nodes[0] << " " << nodes[1] << " " 
                << nodes[2] << " " << nodes[3];
        else if ( name==-QUAD9 )
          outdx << nodes[0] << " " << nodes[2] << " " 
                << nodes[6] << " " << nodes[8];
        else if ( name==-QUAD16 )
          outdx << nodes[0] << " " << nodes[3] << " " 
                << nodes[12] << " " << nodes[15];
        else if ( name==-TET4 )
          outdx << nodes[0] << " " << nodes[1] << " " 
                << nodes[2] << " " << nodes[3];
        else if ( name==-TET10 )
          outdx << nodes[0] << " " << nodes[2] << " " 
                << nodes[5] << " " << nodes[9];
        else if ( name==-HEX8 )
          outdx << nodes[0] << " " << nodes[1] << " " 
                << nodes[2] << " " << nodes[3] << " "
                << nodes[4] << " " << nodes[5] << " "
                << nodes[6] << " " << nodes[7];
        else if ( name==-HEX27 )
          outdx << nodes[0] << " " << nodes[2] << " " 
                << nodes[6] << " " << nodes[8] << " "
                << nodes[18] << " " << nodes[20] << " "
                << nodes[24] << " " << nodes[26];
        else if ( name==-HEX64 )
          outdx << nodes[0] << " " << nodes[3] << " " 
                << nodes[12] << " " << nodes[15] << " "
                << nodes[48] << " " << nodes[51] << " "
                << nodes[60] << " " << nodes[63];
        else if ( name==-SPRING1 )
          outdx << nodes[0] ;
        else if ( name==-SPRING2 )
          outdx << nodes[0] << " " << nodes[1];
        else if ( name==-TRUSS )
          outdx << nodes[0] << " " << nodes[1];
        else {
          assert ( name==-TRUSSBEAM );
          outdx << nodes[0] << " " << nodes[1];
        }
        outdx << "\n";
      }
    }
    outdx << "attribute \"element type\" string \"";
    if      ( name==-BEAM || name==-TRUSS || name==-TRUSSBEAM ) 
      outdx << "lines";   
    else if ( name==-BAR2 || name==-BAR3 || name==-BAR4 ) 
      outdx << "lines";   
    else if ( name==-SPRING2 || name==-CONTACTSPRING ) 
      outdx << "lines";   
    else if ( name==-TRIA3 || name==-TRIA6 ) 
      outdx << "triangles";   
    else if ( name==-QUAD4 || name==-QUAD9 || name==-QUAD16 ) 
      outdx << "quads";   
    else if ( name==-TET4 || name==-TET10 ) 
      outdx << "tetrahedra";   
    else if ( name==-HEX8 || name==-HEX27 || name==-HEX64 ) 
      outdx << "cubes";   
    outdx <<"\"\n";
    outdx << "attribute \"ref\" string \"positions\"\n";   
  }

    // print element coordinates
  outdx << "#\n";
  outdx << "object \"element_middle"
    << " time=" << time_current << "\" ";
  outdx << "class array type double rank 1 shape "
    << ndim << " items " << 1+max_element << " data follows\n";
  for ( element=0; element<=max_element; element++ ) {
    db( ELEMENT_MIDDLE, element, idum, average_coord, 
      ldum, VERSION_PRINT, GET );
    for ( idim=0; idim<ndim; idim++ ) {
      tmp = average_coord[idim];
      if ( scalar_dabs(tmp)<EPS_SMALL ) tmp = 0.;
      outdx << tmp << " ";
    }
    outdx << "\n";
  }

    // print element numbers
  outdx << "#\n";
  outdx << "object \"element_numbers"
    << " time=" << time_current << "\" ";
  outdx << "class array type string rank 1 shape 5"
    " items " << 1+max_element << " data follows\n";
  for ( element=0; element<=max_element; element++ ) {
    outdx << "\"" << old_element_numbers[element] << "\"";
    outdx << "\n";
  }
  outdx << "attribute \"dep\" string \"positions\"\n";   
  outdx << "#\n";
  outdx << "object \"" << "element_numbers" << " field"
    << " time=" << time_current << "\" ";
  outdx << "class field\n";
  outdx << "component \"positions\" value \"element_middle" 
    << " time=" << time_current << "\"\n";
  outdx << "component \"data\" value \"" << "element_numbers" 
    << " time=" << time_current << "\"\n";

    // write scalars fields
  for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
    for ( itype=0; itype<ntype; itype++ ) {
      rank = 0;
      iuknwn = ipuknwn*nder;
      strcpy( dof_type_name, db_name(dof_label[iuknwn]) );
      outdx << "#\n";
      outdx << "object \"" << dof_type_name << 
        " group " << type_group[itype] << " time=" << time_current << "\" ";
      outdx << "class array type double rank " << rank << " ";
      outdx << "items " << 1+max_node << " data follows\n";
      for ( inod=0; inod<=max_node; inod++ ) {
        node_dof = db_dbl( NODE_DOF, inod, VERSION_PRINT );
        tmp = node_dof[iuknwn];
        if ( scalar_dabs(tmp)<EPS_SMALL ) tmp = 0.;
        outdx << tmp;
        outdx << "\n";
      }
      outdx << "attribute \"dep\" string \"positions\"\n";   
      outdx << "#\n";
      outdx << "object \"" << dof_type_name << " field"
        << " group " << type_group[itype] << " time=" << time_current << "\" ";
      outdx << "class field\n";
      outdx << "component \"positions\" value \"nodes" << " time=" << time_current << "\"\n";
      outdx << "component \"connections\" value \"elements group " <<
        type_group[itype] << " time=" << time_current << "\"\n";
      outdx << "component \"data\" value \"" << dof_type_name 
        << " group " << type_group[itype] << " time=" << time_current << "\"\n";
      outdx << "attribute \"name\" string \"" << dof_type_name 
        << " group " << type_group[itype] << " time=" << time_current << "\"\n";
    }
  }

    // write "mesh" field (just for viewing the shape of the numerical domain)
  for ( itype=0; itype<ntype; itype++ ) {
    rank = 0;
    strcpy( dof_type_name, "mesh" );
    outdx << "#\n";
    outdx << "object \"" << dof_type_name << 
      " group " << type_group[itype] << " time=" << time_current << "\" ";
    outdx << "class array type double rank " << rank << " ";
    outdx << "items " << 1+max_node << " data follows\n";
    for ( inod=0; inod<=max_node; inod++ ) {
      tmp = 1.; // dummy value
      outdx << tmp;
      outdx << "\n";
    }
    outdx << "attribute \"dep\" string \"positions\"\n";   
    outdx << "#\n";
    outdx << "object \"" << dof_type_name << " field"
      << " group " << type_group[itype] << " time=" << time_current << "\" ";
    outdx << "class field\n";
    outdx << "component \"positions\" value \"nodes" << " time=" << time_current << "\"\n";
    outdx << "component \"connections\" value \"elements group " <<
      type_group[itype] << " time=" << time_current << "\"\n";
    outdx << "component \"data\" value \"" << dof_type_name 
        << " group " << type_group[itype] << " time=" << time_current << "\"\n";
    outdx << "attribute \"name\" string \"" << dof_type_name 
        << " group " << type_group[itype] << " time=" << time_current << "\"\n";
  }

    // write vectors and tensors for primary unknowns
  ipuknwn = 0;
  if ( nuknwn==0 ) 
    ready = 1;
  else
    ready = 0;
  while ( !ready ) {
    iuknwn = ipuknwn*nder;
    if      ( dof_scal_vec_mat[iuknwn]==-SCALAR ) {
      rank = 0;
      nval = dof_amount[iuknwn];
    }
    else if ( dof_scal_vec_mat[iuknwn]==-VECTOR ) {
      rank = 1;
      nval = dof_amount[iuknwn];
    }
    else {
      assert( dof_scal_vec_mat[iuknwn]==-MATRIX );
      rank = 2;
      nval = dof_amount[iuknwn];
    }
    if ( rank>0 ) {
      strcpy( dof_type_name, db_name(dof_type[iuknwn]) );
      outdx << "#\n";
      for ( itype=0; itype<ntype; itype++ ) {
        outdx << "object \"" << dof_type_name << 
          " group " << type_group[itype] << " time=" << time_current << "\" ";
        outdx << "class array type double rank " << rank << " ";
        if      ( dof_scal_vec_mat[iuknwn]==-VECTOR ) 
          outdx << "shape " << nval << " ";
        else if ( dof_scal_vec_mat[iuknwn]==-MATRIX ) 
          outdx << "shape 3 3 ";
        outdx << "items " << 1+max_node << " data follows\n";
        for ( inod=0; inod<=max_node; inod++ ) {
          node_dof = db_dbl( NODE_DOF, inod, VERSION_PRINT );
          if ( dof_scal_vec_mat[iuknwn]==-VECTOR ) {
            n = dof_amount[iuknwn];
            for ( idim=0; idim<n; idim++ ) {
              indx = iuknwn+idim*nder;
              tmp = node_dof[indx];
              if ( scalar_dabs(tmp)<EPS_SMALL ) tmp = 0.;
              outdx << tmp << " ";
            }
          }
          else {
            assert( dof_scal_vec_mat[iuknwn]==-MATRIX );
            for ( idim=0; idim<MDIM; idim++ ) {
              for ( jdim=0; jdim<MDIM; jdim++ ) {
                indx = iuknwn + stress_indx(idim,jdim)*nder;
                tmp = node_dof[indx];
                if ( scalar_dabs(tmp)<EPS_SMALL ) tmp = 0.;
                outdx << tmp << " ";
              }
            }
          }
          outdx << "\n";
        }
        outdx << "attribute \"dep\" string \"positions\"\n";   
        outdx << "#\n";
        outdx << "object \"" << dof_type_name << " field"
          << " group " << type_group[itype] << " time=" << time_current << "\" ";
        outdx << "class field\n";
        outdx << "component \"positions\" value \"nodes" << " time=" << time_current << "\"\n";
        outdx << "component \"connections\" value \"elements group " <<
          type_group[itype] << " time=" << time_current << "\"\n";
        outdx << "component \"data\" value \"" << dof_type_name 
        << " group " << type_group[itype] << " time=" << time_current << "\"\n";
        outdx << "attribute \"name\" string \"" << dof_type_name 
        << " group " << type_group[itype] << " time=" << time_current << "\"\n";
      }
    }
    ipuknwn += nval;
    if ( ipuknwn>=npuknwn ) ready = 1;
  }


    // post_calcul 
  if ( db_active_index( POST_CALCUL, 0, VERSION_NORMAL ) ) {

      // post_calcul scalars
    for ( icalcul=0; icalcul<length_post_calcul_scal_vec_mat; icalcul++ ) {
      if ( post_calcul_scal_vec_mat[icalcul]==-SCALAR ) {
        rank = 0;
        strcpy( outname, post_calcul_names_without_extension[icalcul] );
        outdx << "#\n";
        for ( itype=0; itype<ntype; itype++ ) {
          outdx << "object \"" << outname << 
            " group " << type_group[itype] << " time=" << time_current << "\" ";
          outdx << "class array type double rank " << rank << " ";
          outdx << "items " << 1+max_node << " data follows\n";
          for ( inod=0; inod<=max_node; inod++ ) {
            node_dof_calcul = db_dbl( NODE_DOF_CALCUL, inod, VERSION_PRINT );
            tmp = node_dof_calcul[icalcul];
            if ( scalar_dabs(tmp)<EPS_SMALL ) tmp = 0.;
            outdx << tmp;
            outdx  << "\n";
          }
          outdx << "attribute \"dep\" string \"positions\"\n";   
          outdx << "#\n";
          outdx << "object \"" << outname << " field"
            << " group " << type_group[itype] << " time=" << time_current << "\" ";
          outdx << "class field\n";
          outdx << "component \"positions\" value \"nodes" << " time=" << time_current << "\"\n";
          outdx << "component \"connections\" value \"elements group " <<
            type_group[itype] << " time=" << time_current << "\"\n";
          outdx << "component \"data\" value \"" << outname 
            << " group " << type_group[itype] << " time=" << time_current << "\"\n";
          outdx << "attribute \"name\" string \"" << outname 
            << " group " << type_group[itype] << " time=" << time_current << "\"\n";
        }
      }
    }

      // post_calcul vectors and tensors
    icalcul = 0;
    if ( length_post_calcul_scal_vec_mat<=0 )
      ready = 1;
    else
      ready = 0;
    while ( !ready ) {
      strcpy( outname, post_calcul_names_without_extension[icalcul] );
      if      ( post_calcul_scal_vec_mat[icalcul]==-SCALAR ) {
        rank = 0;
        nval = 1;
      }
      else if ( post_calcul_scal_vec_mat[icalcul]==-VECTOR ) {
        rank = 1;
        nval = MDIM;
      }
      else {
        assert( post_calcul_scal_vec_mat[icalcul]==-MATRIX );
        rank = 2;
        nval = 6;
      }
      if ( rank>0 ) {
        outdx << "#\n";
        for ( itype=0; itype<ntype; itype++ ) {
          outdx << "object \"" << outname << 
            " group " << type_group[itype] << " time=" << time_current << "\" ";
          outdx << "class array type double rank " << rank << " ";
          if      ( post_calcul_scal_vec_mat[icalcul]==-VECTOR ) 
            outdx << "shape " << nval << " ";
          else if ( post_calcul_scal_vec_mat[icalcul]==-MATRIX ) 
            outdx << "shape 3 3 ";
          outdx << "items " << 1+max_node << " data follows\n";
          for ( inod=0; inod<=max_node; inod++ ) {
            node_dof_calcul = db_dbl( NODE_DOF_CALCUL, inod, VERSION_PRINT );
            if      ( post_calcul_scal_vec_mat[icalcul]==-SCALAR ) {
              tmp = node_dof_calcul[icalcul];
              if ( scalar_dabs(tmp)<EPS_SMALL ) tmp = 0.;
              outdx << tmp;
            }
            else if ( post_calcul_scal_vec_mat[icalcul]==-VECTOR ) {
              for ( idim=0; idim<MDIM; idim++ ) {
                indx = icalcul+idim;
                tmp = node_dof_calcul[indx];
                if ( scalar_dabs(tmp)<EPS_SMALL ) tmp = 0.;
                outdx << tmp << " ";
              }
            }
            else {
              assert( post_calcul_scal_vec_mat[icalcul]==-MATRIX );
              for ( idim=0; idim<MDIM; idim++ ) {
                for ( jdim=0; jdim<MDIM; jdim++ ) {
                  indx = icalcul + stress_indx(idim,jdim)*nder;
                  tmp = node_dof_calcul[indx];
                  if ( scalar_dabs(tmp)<EPS_SMALL ) tmp = 0.;
                  outdx << tmp << " ";
                }
              }
            }
            outdx  << "\n";
          }
          outdx << "attribute \"dep\" string \"positions\"\n";   
          outdx << "#\n";
          outdx << "object \"" << outname << " field" 
            << " group " << type_group[itype] << " time=" << time_current << "\" ";
          outdx << "class field\n";
          outdx << "component \"positions\" value \"nodes" << " time=" << time_current << "\"\n";
          outdx << "component \"connections\" value \"elements " 
            << " group " << type_group[itype] << " time=" << time_current << "\"\n";
          outdx << "component \"data\" value \"" << outname 
            << " group " << type_group[itype] << " time=" << time_current << "\"\n";
          outdx << "attribute \"name\" string \"" << outname 
            << " group " << type_group[itype] << " time=" << time_current << "\"\n";
        }
      }
      icalcul += nval;
      ready = (icalcul>=length_post_calcul_scal_vec_mat);
    }
  }

  db_version_delete( VERSION_PRINT );

    // close dx file
  close_dx_file:

    // series
  if ( final_call==-YES ) {

      // elements series
    for ( itype=0; itype<ntype; itype++ ) {
      outdx << "#\n";
      outdx << "object \"elements" 
        << " group " << type_group[itype] << " series\" class series\n";
      for ( itime=0; itime<length_control_print_dx_time; itime++ ) {
        outdx << "member " << itime << " position " << control_print_dx_time[itime] << " ";
        outdx << "value \"elements"
          << " group " << type_group[itype] << " time=" << control_print_dx_time[itime] << "\"\n";
      }
    }

      // field variables series
    if (nfield>0 ) {
      for ( ifield=0; ifield<nfield; ifield++ ) {
        for ( itype=0; itype<ntype; itype++ ) {
          outdx << "#\n";
          outdx << "object \""  << field_names[ifield] 
            << " group " << type_group[itype] << " series\" class series\n";
          for ( itime=0; itime<length_control_print_dx_time; itime++ ) {
            outdx << "member " << itime << " position " << control_print_dx_time[itime] << " ";
            outdx << "value \"" << field_names[ifield] << " field"
              << " group " << type_group[itype] << " time=" << control_print_dx_time[itime] << "\"\n";
          }
        }
      }
    }

      // mesh series
    for ( itype=0; itype<ntype; itype++ ) {
      outdx << "#\n";
      outdx << "object \""  << "mesh" 
        << " group " << type_group[itype] << " series\" class series\n";
      for ( itime=0; itime<length_control_print_dx_time; itime++ ) {
        outdx << "member " << itime << " position " << control_print_dx_time[itime] << " ";
        outdx << "value \"" << "mesh" << " field"
          << " group " << type_group[itype] << " time=" << control_print_dx_time[itime] << "\"\n";
      }
   }

      // node_numbers series
    outdx << "#\n";
    outdx << "object \""  << "node_numbers" 
      << " series\" class series\n";
    for ( itime=0; itime<length_control_print_dx_time; itime++ ) {
      outdx << "member " << itime << " position " << control_print_dx_time[itime] << " ";
      outdx << "value \"" << "node_numbers" << " field"
        << " time=" << control_print_dx_time[itime] << "\"\n";
    }

      // element_numbers series
    outdx << "#\n";
    outdx << "object \""  << "element_numbers" 
      << " series\" class series\n";
    for ( itime=0; itime<length_control_print_dx_time; itime++ ) {
      outdx << "member " << itime << " position " << control_print_dx_time[itime] << " ";
      outdx << "value \"" << "element_numbers" << " field"
        << " time=" << control_print_dx_time[itime] << "\"\n";
    }     

      // everything grouped

    outdx << "#\n";
    outdx << "object \"mesh grouped" << "\" class group\n";
    for ( itype=0; itype<ntype; itype++ ) {
      outdx << "member \"" << "mesh" 
        << " group " << type_group[itype]
        << "\" \"" << "mesh" 
        << " group " << type_group[itype] << " series\"\n";
    }

    outdx << "#\n";
    outdx << "object \"variables grouped" << "\" class group\n";
    for ( itype=0; itype<ntype; itype++ ) {
      outdx << "member \"" << "no_result" 
        << " group " << type_group[itype]
        << "\" \"" << "mesh" 
        << " group " << type_group[itype] << " series\"\n";
      for ( ifield=0; ifield<nfield; ifield++ ) {
        outdx << "member \"" << field_names[ifield] 
          << " group " << type_group[itype]
          << "\" \"" << field_names[ifield] 
          << " group " << type_group[itype] << " series\"\n";
      }
    }

    outdx << "#\n";
    outdx << "end\n";

  }
  end_of_dx:

    // clean
  outdx.close();
  delete[] control_print_dx_time;
  delete[] old_node_numbers;
  delete[] old_element_numbers;

  if ( swit ) pri( "Out routine PRINT_DX" );
}
