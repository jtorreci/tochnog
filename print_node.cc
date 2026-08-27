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

// print_node - control_print_node (manual Professional 6.330) with the
// companion records control_print_node_angular (6.331),
// control_print_node_angular_middle (6.332), control_print_node_geometry
// (6.333), control_print_node_sort (6.334) and control_print_node_zero
// (6.335): print NODAL data records (node_dof, node_dof_calcul, or ANY
// record whose name starts with "node") to plain ASCII files with one
// "x y z <value>" line per node (in 1D only x). One file per selected
// part (manual example: control_print_node index -node_dof -velx -vely
// produces velx.index and vely.index with x y velx / x y vely columns).
//
// Parts (number_0 number_1 ...) select the values to print:
//   - for node_dof: dof labels (-velx, -vely, ...) or numbers 0, 1, ...
//   - for node_dof_calcul: post_calcul labels or numbers
//   - for any other node record: numbers only
// No parts -> ALL parts of the record (documented decision).
//
// File naming (documented decision): label.<index> for dof labels and
// post_calcul items (consistent with control_print_dof_line), and
// <record_name>_<n>.<index> for numeric parts (e.g. node_dof_0.10).

// Angle in DEGREES of the position (coord) around the middle point
// (manual 6.331): axis 0 = from +x directed to +y (atan2(y-ym, x-xm)),
// axis 1 = from +y directed to +z, axis 2 = from +x directed to +z.
static double node_angle_degrees( double coord[], double middle[],
  long int axis )

{
  double dx=0., dy=0.;

  if      ( axis==0 ) { dx = coord[0]-middle[0]; dy = coord[1]-middle[1]; }
  else if ( axis==1 ) { dx = coord[1]-middle[1]; dy = coord[2]-middle[2]; }
  else                { dx = coord[0]-middle[0]; dy = coord[2]-middle[2]; }
  return 180./PIRAD * atan2( dy, dx );
}

void print_node( long int icontrol, long int ival[], long int nval )

{
  long int data_item=0, inod=0, idim=0, max_node=0, swit=0, ldum=0,
    ipart=0, ifile=0, nfiles=0, i=0, nlist=0, sort_axis=-1,
    angular=0, ang_axis=0, zero=-YES, ngeom=0, in_geometry=0,
    geom_node_type=NODE, ncalcul=0, icalcul=0, len=0, have=0,
    matched=0, is_match=0, idum[1], *dof_label=NULL,
    *geom_entity=NULL, *post_calcul_scal_vec_mat=NULL,
    *post_calcul_unknown_operat=NULL, *node_list=NULL, *order=NULL,
    *ival2=NULL;
  double ddum[1], coord[MDIM], middle[MDIM], *key=NULL, *dval=NULL,
    factor=0., normal[MDIM], penetration=0., projection[MDIM], value=0.;
  char filename[MCHAR], str[MCHAR];
  static char base[DATA_ITEM_SIZE][MCHAR]; // file base name per part
  static long int vind[DATA_ITEM_SIZE];    // value index per part

  swit = set_swit(-1,-1,"print_node");
  if ( swit ) pri( "In routine PRINT_NODE" );

  // data item: a nodal record whose name starts with "node"
  if ( nval<1 || ival[0]>=0 ) db_error( CONTROL_PRINT_NODE, icontrol );
  data_item = labs( ival[0] );
  if ( strncmp( db_name(data_item), "node", 4 ) ||
       ( db_type(data_item)!=INTEGER &&
         db_type(data_item)!=DOUBLE_PRECISION ) )
    db_error( CONTROL_PRINT_NODE, icontrol );

  db_highest_index( NODE, max_node, VERSION_NORMAL );
  if ( max_node<0 ) return;

  // control_print_node_angular (6.331): angle instead of coordinates.
  // 1D cannot use it; 2D only -yes -yes (and switch_z NOT given); 3D the
  // three combinations -yes -yes -no / -no -yes -yes / -yes -no -yes.
  if ( db_active_index( CONTROL_PRINT_NODE_ANGULAR, icontrol,
       VERSION_NORMAL ) ) {
    long int ang[DATA_ITEM_SIZE], ang_len=0;
    if ( ndim==1 ) db_error( CONTROL_PRINT_NODE_ANGULAR, icontrol );
    db( CONTROL_PRINT_NODE_ANGULAR, icontrol, ang, ddum, ang_len,
      VERSION_NORMAL, GET );
    if ( ndim==2 ) {
      if ( ang_len!=2 || ang[0]!=-YES || ang[1]!=-YES )
        db_error( CONTROL_PRINT_NODE_ANGULAR, icontrol );
      angular = 1; ang_axis = 0;
    }
    else {
      if ( ang_len!=3 ) db_error( CONTROL_PRINT_NODE_ANGULAR, icontrol );
      if      ( ang[0]==-YES && ang[1]==-YES && ang[2]==-NO )
        { angular = 1; ang_axis = 0; }
      else if ( ang[0]==-NO && ang[1]==-YES && ang[2]==-YES )
        { angular = 1; ang_axis = 1; }
      else if ( ang[0]==-YES && ang[1]==-NO && ang[2]==-YES )
        { angular = 1; ang_axis = 2; }
      else db_error( CONTROL_PRINT_NODE_ANGULAR, icontrol );
    }
  }
  array_set( middle, 0., MDIM ); // default middle point (0,0,0)
  if ( db_active_index( CONTROL_PRINT_NODE_ANGULAR_MIDDLE, icontrol,
       VERSION_NORMAL ) ) {
    double md[MDIM];
    long int md_len=0;
    db( CONTROL_PRINT_NODE_ANGULAR_MIDDLE, icontrol, idum, md, md_len,
      VERSION_NORMAL, GET );
    if ( ( ndim==2 && md_len!=2 ) || ( ndim==3 && md_len!=3 ) )
      db_error( CONTROL_PRINT_NODE_ANGULAR_MIDDLE, icontrol );
    array_move( md, middle, md_len );
  }

  // control_print_node_geometry (6.333): only nodes in the geometry
  geom_entity = get_new_int(DATA_ITEM_SIZE);
  if ( db_active_index( CONTROL_PRINT_NODE_GEOMETRY, icontrol,
       VERSION_NORMAL ) ) {
    db( CONTROL_PRINT_NODE_GEOMETRY, icontrol, geom_entity, ddum, ngeom,
      VERSION_NORMAL, GET );
    if ( ngeom!=2 || geom_entity[0]>=0 ||
         db_data_class(geom_entity[0])!=GEOMETRY )
      db_error( CONTROL_PRINT_NODE_GEOMETRY, icontrol );
    long int max_nsr=-1;
    db_max_index( NODE_START_REFINED, max_nsr, VERSION_NORMAL, GET );
    if ( max_nsr>=0 ) geom_node_type = NODE_START_REFINED;
  }

  // control_print_node_sort (6.334): -angle with angular (the manual's
  // "control_print_node_method" is a typo), -x / -y (2D/3D) / -z (3D)
  // otherwise. Ascending order.
  if ( db_active_index( CONTROL_PRINT_NODE_SORT, icontrol,
       VERSION_NORMAL ) ) {
    long int sm=0;
    db( CONTROL_PRINT_NODE_SORT, icontrol, &sm, ddum, ldum,
      VERSION_NORMAL, GET );
    if ( angular ) {
      if ( sm!=-ANGLE ) db_error( CONTROL_PRINT_NODE_SORT, icontrol );
      sort_axis = 3;
    }
    else {
      if      ( sm==-X ) sort_axis = 0;
      else if ( sm==-Y && ndim>=2 ) sort_axis = 1;
      else if ( sm==-Z && ndim==3 ) sort_axis = 2;
      else db_error( CONTROL_PRINT_NODE_SORT, icontrol );
    }
  }

  // control_print_node_zero (6.335): default -yes prints zero valued
  // results; -no suppresses them (exact zero comparison, decision)
  db( CONTROL_PRINT_NODE_ZERO, icontrol, &zero, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  if ( zero!=-YES && zero!=-NO )
    db_error( CONTROL_PRINT_NODE_ZERO, icontrol );

  // post_calcul administration for node_dof_calcul
  if ( data_item==NODE_DOF_CALCUL ) {
    ncalcul = db_len( POST_CALCUL_SCAL_VEC_MAT, 0, VERSION_NORMAL );
    if ( ncalcul<=0 ) db_error( CONTROL_PRINT_NODE, icontrol );
    post_calcul_scal_vec_mat = get_new_int(DATA_ITEM_SIZE);
    post_calcul_unknown_operat = get_new_int(DATA_ITEM_SIZE);
    db( POST_CALCUL_SCAL_VEC_MAT, 0, post_calcul_scal_vec_mat, ddum, ldum,
      VERSION_NORMAL, GET );
    db( POST_CALCUL_UNKNOWN_OPERAT, 0, post_calcul_unknown_operat, ddum,
      ldum, VERSION_NORMAL, GET );
  }

  dof_label = get_new_int(MUKNWN);
  db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );

  // expand the requested parts into (file base name, value index) pairs.
  // No parts -> ALL parts (nuknwn for node_dof, ncalcul for
  // node_dof_calcul, the record length otherwise).
  nfiles = 0;
  if ( nval<2 ) {
    long int np=0, p=0;
    if      ( data_item==NODE_DOF ) np = nuknwn;
    else if ( data_item==NODE_DOF_CALCUL ) np = ncalcul;
    else np = db_data_length(data_item);
    for ( p=0; p<np && nfiles<DATA_ITEM_SIZE; p++ ) {
      snprintf( base[nfiles], MCHAR, "%s_%ld", db_name(data_item), p );
      vind[nfiles] = p;
      nfiles++;
    }
  }
  else {
    for ( ipart=1; ipart<nval && nfiles<DATA_ITEM_SIZE; ipart++ ) {
      long int part = ival[ipart];
      if ( part<0 ) {
        if ( data_item==NODE_DOF ) {
          long int idx=0;
          if ( !array_member( dof_label, part, nuknwn, idx ) )
            db_error( CONTROL_PRINT_NODE, icontrol );
          strcpy( base[nfiles], db_name(labs(part)) );
          vind[nfiles] = idx;
          nfiles++;
        }
        else if ( data_item==NODE_DOF_CALCUL ) {
          // GNU adaptation: post_calcul_label does not exist in the GNU;
          // a name matches the underlying unknown exactly (e.g.
          // -materi_velocity selects every operator of that unknown) or
          // is a substring of its label (post_calcul_names, e.g. -sigyy
          // matches asigyy). One file per matched calcul item.
          matched = 0;
          for ( icalcul=0; icalcul<ncalcul && nfiles<DATA_ITEM_SIZE;
               icalcul++ ) {
            is_match = 0;
            if ( post_calcul_unknown_operat[icalcul*2+0]==part )
              is_match = 1;
            else if ( strstr( post_calcul_names[icalcul],
                              db_name(part) ) )
              is_match = 1;
            if ( is_match ) {
              strcpy( base[nfiles],
                post_calcul_names_without_extension[icalcul] );
              vind[nfiles] = icalcul;
              nfiles++;
              matched = 1;
            }
          }
          if ( !matched ) db_error( CONTROL_PRINT_NODE, icontrol );
        }
        else db_error( CONTROL_PRINT_NODE, icontrol );
      }
      else {
        snprintf( base[nfiles], MCHAR, "%s_%ld", db_name(data_item),
          part );
        vind[nfiles] = part;
        nfiles++;
      }
    }
  }
  if ( nfiles==0 ) {
    delete[] geom_entity; delete[] dof_label;
    if ( post_calcul_scal_vec_mat ) delete[] post_calcul_scal_vec_mat;
    if ( post_calcul_unknown_operat ) delete[] post_calcul_unknown_operat;
    if ( swit ) pri( "Out routine PRINT_NODE" );
    return;
  }

  // collect the printable nodes (geometry filter) with their sort key
  node_list = get_new_int(max_node+1);
  key = get_new_dbl(max_node+1);
  for ( inod=0; inod<=max_node; inod++ ) {
    if ( !db_active_index( NODE, inod, VERSION_NORMAL ) ) continue;
    if ( ngeom==2 ) {
      geometry( inod, ddum, geom_entity, in_geometry, factor, normal,
        penetration, projection, geom_node_type, PROJECT_EXACT,
        VERSION_NORMAL );
      if ( !in_geometry ) continue;
    }
    node_list[nlist] = inod;
    key[nlist] = 0.;
    if ( sort_axis>=0 ) {
      db( NODE, inod, idum, coord, ldum, VERSION_NORMAL, GET );
      if ( angular ) key[nlist] = node_angle_degrees( coord, middle,
        ang_axis );
      else           key[nlist] = coord[sort_axis];
    }
    nlist++;
  }

  // sort ascending (collect+sort pattern of
  // print_interface_stress.cc:72-89)
  order = get_new_int(nlist);
  for ( i=0; i<nlist; i++ ) order[i] = i;
  if ( sort_axis>=0 ) {
    for ( i=1; i<nlist; i++ ) {
      long int k = i;
      while ( k>0 ) {
        long int a = order[k-1], b = order[k];
        if ( key[b]>=key[a] ) break;
        order[k] = a; order[k-1] = b; k--;
      }
    }
  }

  // write one file per part
  ival2 = get_new_int(DATA_ITEM_SIZE);
  dval = get_new_dbl(DATA_ITEM_SIZE);
  for ( ifile=0; ifile<nfiles; ifile++ ) {
    strcpy( filename, base[ifile] );
    strcat( filename, "." );
    long_to_a( icontrol, str );
    strcat( filename, str );
    ofstream out( filename, ios::app );
    out.precision(TN_PRECISION);
    for ( i=0; i<nlist; i++ ) {
      inod = node_list[ ( sort_axis>=0 ) ? order[i] : i ];
      if ( !db_active_index( data_item, inod, VERSION_NORMAL ) ) continue;
      have = 0;
      if ( db_type(data_item)==INTEGER ) {
        db( data_item, inod, ival2, ddum, len, VERSION_NORMAL, GET );
        if ( vind[ifile]<len ) { value = ival2[vind[ifile]]; have = 1; }
      }
      else {
        db( data_item, inod, idum, dval, len, VERSION_NORMAL, GET );
        if ( vind[ifile]<len ) { value = dval[vind[ifile]]; have = 1; }
      }
      if ( !have ) continue; // part beyond the record length: omitted
      if ( zero==-NO && value==0. ) continue;
      db( NODE, inod, idum, coord, ldum, VERSION_NORMAL, GET );
      if ( angular )
        out << node_angle_degrees( coord, middle, ang_axis ) << " ";
      else
        for ( idim=0; idim<ndim; idim++ ) out << coord[idim] << " ";
      if ( db_type(data_item)==INTEGER )
        out << (long int) value << "\n";
      else
        out << value << "\n";
    }
    out.close();
  }

  delete[] geom_entity;
  delete[] dof_label;
  if ( post_calcul_scal_vec_mat ) delete[] post_calcul_scal_vec_mat;
  if ( post_calcul_unknown_operat ) delete[] post_calcul_unknown_operat;
  delete[] node_list;
  delete[] key;
  delete[] order;
  delete[] ival2;
  delete[] dval;

  if ( swit ) pri( "Out routine PRINT_NODE" );
}
