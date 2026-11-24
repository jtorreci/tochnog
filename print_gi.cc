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

// icontrol: index of the control_print_gid that prints, -1 for the
// end-of-run print (only the global print_gid_* records apply then).
void print_gid( long int task, long int icontrol )

{
  if      ( print_gid_6( task, icontrol ) ) {
    return;
  }
  else if ( print_gid_5( task ) ) {
    return;
  }
  else {
    pri( "Internal error: no valid GID printing found." );
    pri( "Please report this error with your input file." );
    exit(TN_EXIT_STATUS);
  }

}

/*
  GiD output options (manual Professional 6.295-6.315 control_print_gid_*,
  6.982-6.989 print_gid_*). print_gid_options_get reads them once per GiD
  print (print_gid_6) for the control_print_gid index that prints; the
  global print_gid_X is the default of control_print_gid_X (print_option).
  Without any of these records the GiD files are the ones the GNU always
  wrote: the records that add results (coord, group, spring2, vectors) only
  add them when they are set; see ProjectDocs/manual-developer/
  print_gid_options.md.
*/

static long int gid_coord=0, gid_group=0, gid_other=-YES,
  gid_node_type=-NODE_START_REFINED, gid_spring2=0, gid_beam_vectors=-NO,
  gid_truss_vector=-NO, gid_element_mpc=-YES, gid_smooth_n=10,
  ngid_dof=0, ngid_dof_calcul=0, ngid_element_group=0, ngid_smooth_dof=0,
  gid_max_node=-1, gid_mesh_activate_gravity=-NO, gid_nold_element=0,
  *gid_old_element_numbers=NULL,
  gid_dof[DATA_ITEM_SIZE], gid_dof_calcul[DATA_ITEM_SIZE],
  gid_element_group[DATA_ITEM_SIZE], gid_smooth_dof[DATA_ITEM_SIZE],
  *gid_mpc_node=NULL;
static double gid_beam_normal[MDIM], gid_truss_normal[MDIM],
  gid_time_current=0.;

static void print_gid_option_error( const char *names )

{
  pri( "Error detected for the GiD output option ", names );
  exit_tn_on_error();
}

// value of a control_print_gid_* record that has no global print_gid_* form
static long int print_gid_control_value( long int control_name,
  long int icontrol, long int default_value )

{
  long int value=default_value, ldum=0;
  double ddum[1];

  if ( icontrol>=0 )
    db( control_name, icontrol, &value, ddum, ldum, VERSION_NORMAL,
      GET_IF_EXISTS );
  return value;
}

static void print_gid_switch_check( long int value, const char *names )

{
  if ( value!=-YES && value!=-NO ) print_gid_option_error( names );
}

// plotting-plane normal of the beam/truss vectors, default 0 0 1
static void print_gid_normal_get( long int control_name, long int icontrol,
  double normal[] )

{
  long int ldum=0, idum[1];

  normal[0] = 0.; normal[1] = 0.; normal[2] = 1.;
  if ( icontrol>=0 )
    db( control_name, icontrol, idum, normal, ldum, VERSION_NORMAL,
      GET_IF_EXISTS );
}

void print_gid_options_get( long int icontrol, long int max_node,
  long int old_node_numbers[], long int nold_element,
  long int old_element_numbers[] )

{
  long int ldum=0, inod=0, max_mpc=-1, impc=0, max_old=-1, length=0,
    idum[1], *mpc=NULL, *mpc_old=NULL;
  double ddum[1];

  // switches and numbers: global default, overruled by the control record
  // defaults of the Professional 02-08-2026 (owner decision, odd/tasks/gid.md
  // L13): coord and element_group are written unless set to -no
  gid_coord = print_option( PRINT_GID_COORD, CONTROL_PRINT_GID_COORD,
    icontrol, -YES );
  gid_group = print_option( PRINT_GID_GROUP, CONTROL_PRINT_GID_GROUP,
    icontrol, -YES );
  gid_node_type = print_option( PRINT_GID_NODE_TYPE,
    CONTROL_PRINT_GID_NODE_TYPE, icontrol, -NODE_START_REFINED );
  gid_spring2 = print_option( PRINT_GID_SPRING2, CONTROL_PRINT_GID_SPRING2,
    icontrol, 1 );
  gid_smooth_n = print_option( PRINT_GID_SMOOTH_N, CONTROL_PRINT_GID_SMOOTH_N,
    icontrol, 10 );
  // records without a global form
  gid_other = print_gid_control_value( CONTROL_PRINT_GID_OTHER, icontrol, -YES );
  gid_beam_vectors = print_gid_control_value( CONTROL_PRINT_GID_BEAM_VECTORS,
    icontrol, -NO );
  gid_truss_vector = print_gid_control_value( CONTROL_PRINT_GID_TRUSS_VECTOR,
    icontrol, -NO );
  gid_element_mpc = print_gid_control_value( CONTROL_PRINT_GID_ELEMENT_MPC,
    icontrol, -YES );
  print_gid_normal_get( CONTROL_PRINT_GID_BEAM_VECTORS_NORMAL, icontrol,
    gid_beam_normal );
  print_gid_normal_get( CONTROL_PRINT_GID_TRUSS_VECTOR_NORMAL, icontrol,
    gid_truss_normal );
  // name lists
  print_option_list( -1, CONTROL_PRINT_GID_DOF, icontrol, gid_dof, ngid_dof );
  print_option_list( -1, CONTROL_PRINT_GID_DOF_CALCUL, icontrol,
    gid_dof_calcul, ngid_dof_calcul );
  print_option_list( -1, CONTROL_PRINT_GID_ELEMENT_GROUP, icontrol,
    gid_element_group, ngid_element_group );
  print_option_list( PRINT_GID_SMOOTH_DOF, CONTROL_PRINT_GID_SMOOTH_DOF,
    icontrol, gid_smooth_dof, ngid_smooth_dof );

  if ( gid_coord!=0 )
    print_gid_switch_check( gid_coord, "control_print_gid_coord/print_gid_coord" );
  if ( gid_group!=0 )
    print_gid_switch_check( gid_group, "control_print_gid_group/print_gid_group" );
  print_gid_switch_check( gid_other, "control_print_gid_other" );
  print_gid_switch_check( gid_beam_vectors, "control_print_gid_beam_vectors" );
  print_gid_switch_check( gid_truss_vector, "control_print_gid_truss_vector" );
  print_gid_switch_check( gid_element_mpc, "control_print_gid_element_mpc" );
  if ( gid_node_type!=-NODE && gid_node_type!=-NODE_START_REFINED &&
       gid_node_type!=-NODE_DEFORMED_MESH )
    print_gid_option_error( "control_print_gid_node_type/print_gid_node_type" );
  if ( gid_spring2!=1 && gid_spring2!=2 )
    print_gid_option_error( "control_print_gid_spring2/print_gid_spring2" );
  if ( gid_smooth_n<0 )
    print_gid_option_error( "control_print_gid_smooth_n/print_gid_smooth_n" );
  // mesh_activate_gravity (6.305/6.985): -no (the Professional default,
  // odd/tasks/gid.md L14) writes only the elements active at this print;
  // print_gid_define (6.983): accepted, no effect (see the user page)
  gid_mesh_activate_gravity = print_option( PRINT_GID_MESH_ACTIVATE_GRAVITY,
    CONTROL_PRINT_GID_MESH_ACTIVATE_GRAVITY, icontrol, -NO );
  print_gid_switch_check( gid_mesh_activate_gravity,
    "control_print_gid_mesh_activate_gravity/print_gid_mesh_activate_gravity" );
  gid_time_current = 0.;
  db( TIME_CURRENT, 0, idum, &gid_time_current, ldum, VERSION_NORMAL,
    GET_IF_EXISTS );
  gid_nold_element = nold_element;
  gid_old_element_numbers = old_element_numbers;
  print_gid_switch_check( print_option( PRINT_GID_DEFINE, -1, -1, -YES ),
    "print_gid_define" );

  // control_print_gid_element_mpc -no: flag the mpc'd nodes, the first
  // node of each mpc_node_number (the slave; = Professional 02-08-2026).
  // The print version is renumbered: old_node_numbers maps back.
  gid_max_node = max_node;
  if ( gid_element_mpc==-NO && max_node>=0 ) {
    gid_mpc_node = get_new_int( max_node+1 );
    array_set( gid_mpc_node, 0, max_node+1 );
    db_max_index( NODE, max_old, VERSION_NORMAL, GET );
    db_max_index( MPC_NODE_NUMBER, max_mpc, VERSION_NORMAL, GET );
    if ( max_old>=0 && max_mpc>=0 ) {
      mpc_old = get_new_int( max_old+1 );
      array_set( mpc_old, 0, max_old+1 );
      mpc = get_new_int( DATA_ITEM_SIZE );
      for ( impc=0; impc<=max_mpc; impc++ ) {
        if ( db( MPC_NODE_NUMBER, impc, mpc, ddum, length, VERSION_NORMAL,
             GET_IF_EXISTS ) && mpc[0]>=0 && mpc[0]<=max_old )
          mpc_old[mpc[0]] = 1;
      }
      for ( inod=0; inod<=max_node; inod++ ) {
        ldum = old_node_numbers[inod];
        if ( ldum>=0 && ldum<=max_old ) gid_mpc_node[inod] = mpc_old[ldum];
      }
      delete[] mpc_old;
      delete[] mpc;
    }
  }
}

void print_gid_options_free( void )

{
  if ( gid_mpc_node ) delete[] gid_mpc_node;
  gid_mpc_node = NULL;
  gid_old_element_numbers = NULL;
  gid_nold_element = 0;
}

// control_print_gid_element_group (6.301), control_print_gid_element_mpc
// (6.302) and control_print_gid_mesh_activate_gravity (6.305): is the
// element (print version numbering) written to the GiD files? spring2
// elements only as 1 node (Point) by default, or 2 nodes with
// control_print_gid_spring2 2.
long int print_gid_element_selected( long int element, long int name,
  long int element_group, long int nnol, long int nodes[] )

{
  long int i=0, found=0, old=element;

  if ( ngid_element_group>0 ) {
    for ( i=0; i<ngid_element_group; i++ )
      if ( gid_element_group[i]==element_group ) found = 1;
    if ( !found ) return 0;
  }
  if ( gid_mpc_node ) {
    for ( i=0; i<nnol; i++ )
      if ( nodes[i]>=0 && nodes[i]<=gid_max_node && gid_mpc_node[nodes[i]] )
        return 0;
  }
  // -no: only the elements active in mesh_activate_gravity_time at this
  // print (the .flavia.msh is written at the first print only)
  if ( gid_mesh_activate_gravity==-NO ) {
    if ( gid_old_element_numbers && element>=0 && element<gid_nold_element )
      old = gid_old_element_numbers[element];
    if ( old>=0 && !mesh_activate_gravity_active( old, gid_time_current ) )
      return 0;
  }
  return 1;
}

// control_print_gid_dof (6.299): dof field of initialisation name dof_type
long int print_gid_dof_selected( long int dof_type )

{
  long int i=0;

  if ( ngid_dof==0 ) return 1;
  for ( i=0; i<ngid_dof; i++ ) {
    if ( gid_dof[i]==-NONE ) return 0;
    if ( gid_dof[i]==dof_type ) return 1;
  }
  return 0;
}

// control_print_gid_dof_calcul (6.300): post field icalcul. A name selects
// the field when it is its unknown (-materi_stress), its operator
// (-mises) or a part of its label (as control_print_vtk_dof_calcul).
// -none selects no field. Manual behaviour (owner decision L8).
long int print_gid_post_selected( long int icalcul )

{
  long int i=0, ldum=0, idum[1], *unknown_operat=NULL, selected=0;
  double ddum[1];

  if ( ngid_dof_calcul==0 ) return 1;
  if ( gid_dof_calcul[0]==-NONE ) return 0;
  unknown_operat = get_new_int( DATA_ITEM_SIZE );
  if ( db( POST_CALCUL_UNKNOWN_OPERAT, 0, unknown_operat, ddum, ldum,
       VERSION_NORMAL, GET_IF_EXISTS ) ) {
    for ( i=0; i<ngid_dof_calcul && !selected; i++ ) {
      if ( unknown_operat[icalcul*2+0]==gid_dof_calcul[i] ||
           unknown_operat[icalcul*2+1]==gid_dof_calcul[i] ||
           strstr( post_calcul_names[icalcul],
             db_name(-gid_dof_calcul[i]) ) ) selected = 1;
    }
  }
  (void) idum;
  delete[] unknown_operat;
  return selected;
}

// control_print_gid_other (6.308): -no drops mesh deform and element_group
long int print_gid_other( void )

{
  return gid_other!=-NO;
}

// control_print_gid_spring2/print_gid_spring2 (6.313/6.987): 1 (Point mesh,
// the default) or 2 (Linear mesh)
long int print_gid_spring2_nodes( void )

{
  return gid_spring2;
}

// control_print_gid_node_type/print_gid_node_type (6.307/6.986): the node
// coordinates of the .flavia.msh. -node: stored coordinates;
// -node_start_refined (default): node_start_refined; -node_deformed_mesh:
// stored coordinates + displacement (as control_print_vtk_node_method).
void print_gid_node_coord( long int inod, double coord[] )

{
  long int idim=0, ldum=0, idum[1];
  double *node_dof=NULL;

  if ( gid_node_type==-NODE_START_REFINED &&
       db_active_index( NODE_START_REFINED, inod, VERSION_PRINT ) ) {
    db( NODE_START_REFINED, inod, idum, coord, ldum, VERSION_PRINT, GET );
    return;
  }
  db( NODE, inod, idum, coord, ldum, VERSION_PRINT, GET );
  if ( gid_node_type==-NODE_DEFORMED_MESH && materi_displacement ) {
    node_dof = db_dbl( NODE_DOF, inod, VERSION_PRINT );
    for ( idim=0; idim<ndim; idim++ )
      coord[idim] += node_dof[dis_indx+idim*nder];
  }
}

// control_print_gid_smooth_dof/print_gid_smooth_dof (6.311/6.988) with
// control_print_gid_smooth_n/print_gid_smooth_n (6.312/6.989, default 10):
// smooth the selected nodal dofs of the print version. Each pass (Jacobi)
// replaces a value by 0.75 itself + 0.25 the mean of the other nodes of the
// elements around the node; a dof bounded at the node keeps its value
// (= Professional 02-08-2026 with -all). -all smooths every dof, names
// (initialisation names or dof labels) only those dofs (manual; the
// Professional smooths nothing with names: owner decision L7).
void print_gid_smooth( long int max_node, long int old_node_numbers[] )

{
  long int i=0, inod=0, jnod=0, ipuknwn=0, iuknwn=0, ismooth=0, inol=0,
    jnol=0, element=0, max_element=-1, max_bounded=-1, length=0, nnol=0,
    old=0, all=0, nsel=0, *smooth=NULL, *el=NULL, *neighbours=NULL,
    *nneighbour=NULL, *dof_label=NULL, *dof_type=NULL, *node_bounded=NULL,
    *fixed=NULL;
  long int const mneighbour = 4*MAXIMUM_NODE;
  double ddum[1], *value=NULL, *value_new=NULL, *node_dof=NULL, sum=0.;

  if ( ngid_smooth_dof==0 || gid_smooth_n==0 || max_node<0 ||
       npuknwn==0 ) return;

  dof_label = get_new_int( MUKNWN );
  dof_type = get_new_int( MUKNWN );
  db( DOF_LABEL, 0, dof_label, ddum, length, VERSION_NORMAL, GET );
  db( DOF_TYPE, 0, dof_type, ddum, length, VERSION_NORMAL, GET );
  smooth = get_new_int( npuknwn );
  for ( i=0; i<ngid_smooth_dof; i++ )
    if ( gid_smooth_dof[i]==-ALL ) all = 1;
  for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
    iuknwn = ipuknwn*nder;
    smooth[ipuknwn] = all;
    for ( i=0; i<ngid_smooth_dof; i++ )
      if ( gid_smooth_dof[i]==dof_label[iuknwn] ||
           gid_smooth_dof[i]==dof_type[iuknwn] ) smooth[ipuknwn] = 1;
    if ( smooth[ipuknwn] ) nsel++;
  }

  if ( nsel>0 ) {
    // neighbour nodes: the other nodes of the elements around each node
    neighbours = get_new_int( (max_node+1)*mneighbour );
    nneighbour = get_new_int( max_node+1 );
    array_set( nneighbour, 0, max_node+1 );
    el = get_new_int( MAXIMUM_NODE+1 );
    db_max_index( ELEMENT, max_element, VERSION_PRINT, GET );
    for ( element=0; element<=max_element; element++ ) {
      if ( !db_active_index( ELEMENT, element, VERSION_PRINT ) ) continue;
      db( ELEMENT, element, el, ddum, length, VERSION_PRINT, GET );
      nnol = length - 1;
      for ( inol=0; inol<nnol; inol++ ) {
        inod = el[1+inol];
        if ( inod<0 || inod>max_node ) continue;
        for ( jnol=0; jnol<nnol; jnol++ ) {
          jnod = el[1+jnol];
          if ( jnod==inod || jnod<0 || jnod>max_node ) continue;
          for ( i=0; i<nneighbour[inod]; i++ )
            if ( neighbours[inod*mneighbour+i]==jnod ) break;
          if ( i==nneighbour[inod] && nneighbour[inod]<mneighbour ) {
            neighbours[inod*mneighbour+i] = jnod;
            nneighbour[inod]++;
          }
        }
      }
    }

    // bounded dofs (node_bounded is kept in the normal version, old
    // node numbers). materi_displacement and materi_velocity_integrated
    // are fixed where materi_velocity is bounded (= Professional).
    fixed = get_new_int( (max_node+1)*npuknwn );
    array_set( fixed, 0, (max_node+1)*npuknwn );
    db_max_index( NODE_BOUNDED, max_bounded, VERSION_NORMAL, GET );
    for ( inod=0; inod<=max_node; inod++ ) {
      old = old_node_numbers[inod];
      if ( old<0 || old>max_bounded ||
           !db_active_index( NODE_BOUNDED, old, VERSION_NORMAL ) ) continue;
      node_bounded = db_int( NODE_BOUNDED, old, VERSION_NORMAL );
      for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
        fixed[inod*npuknwn+ipuknwn] = node_bounded[ipuknwn];
        iuknwn = ipuknwn*nder;
        if ( materi_velocity && materi_displacement && iuknwn>=dis_indx &&
             iuknwn<dis_indx+ndim*nder )
          fixed[inod*npuknwn+ipuknwn] |=
            node_bounded[(vel_indx+iuknwn-dis_indx)/nder];
        if ( materi_velocity && materi_velocity_integrated &&
             iuknwn>=veli_indx && iuknwn<veli_indx+ndim*nder )
          fixed[inod*npuknwn+ipuknwn] |=
            node_bounded[(vel_indx+iuknwn-veli_indx)/nder];
      }
    }

    value = get_new_dbl( max_node+1 );
    value_new = get_new_dbl( max_node+1 );
    for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
      if ( !smooth[ipuknwn] ) continue;
      iuknwn = ipuknwn*nder;
      for ( inod=0; inod<=max_node; inod++ ) {
        value[inod] = 0.;
        if ( db_active_index( NODE_DOF, inod, VERSION_PRINT ) )
          value[inod] = db_dbl( NODE_DOF, inod, VERSION_PRINT )[iuknwn];
      }
      for ( ismooth=0; ismooth<gid_smooth_n; ismooth++ ) {
        for ( inod=0; inod<=max_node; inod++ ) {
          value_new[inod] = value[inod];
          if ( nneighbour[inod]==0 || fixed[inod*npuknwn+ipuknwn] ) continue;
          sum = 0.;
          for ( i=0; i<nneighbour[inod]; i++ )
            sum += value[neighbours[inod*mneighbour+i]];
          value_new[inod] = 0.75*value[inod] + 0.25*sum/nneighbour[inod];
        }
        array_move( value_new, value, max_node+1 );
      }
      for ( inod=0; inod<=max_node; inod++ ) {
        if ( !db_active_index( NODE_DOF, inod, VERSION_PRINT ) ) continue;
        node_dof = db_dbl( NODE_DOF, inod, VERSION_PRINT );
        node_dof[iuknwn] = value[inod];
      }
    }
    delete[] neighbours;
    delete[] nneighbour;
    delete[] el;
    delete[] fixed;
    delete[] value;
    delete[] value_new;
  }
  delete[] smooth;
  delete[] dof_label;
  delete[] dof_type;
}

// element types of the GiD Gauss point results (the mesh types of
// print_gid_6 plus spring2)
static const long int gid_ntype = 12;
static long int print_gid_type( long int itype )

{
  const long int types[] = { -TRIA3, -TRIA6, -QUAD4, -QUAD9, -TET4, -TET10,
    -HEX8, -HEX27, -TRUSS, -BEAM, -TRUSSBEAM, -SPRING2 };
  return types[itype];
}

static long int print_gid_type_linear( long int name )

{
  return name==-TRUSS || name==-BEAM || name==-TRUSSBEAM || name==-SPRING2;
}

static const char *print_gid_type_gid_name( long int name )

{
  if      ( name==-TRIA3 || name==-TRIA6 ) return "Triangle";
  else if ( name==-QUAD4 || name==-QUAD9 ) return "Quadrilateral";
  else if ( name==-TET4 || name==-TET10 ) return "Tetrahedra";
  else if ( name==-HEX8 || name==-HEX27 ) return "Hexahedra";
  return "Linear";
}

// is any element result (control_print_gid_group, _beam_vectors,
// _truss_vector, or print_gid_group) requested for any control index?
// The Gauss points are declared once, when the .flavia.res is created.
static long int print_gid_gauss_needed( void )

{
  long int i=0, icontrol=0, max_index=-1, value=0, ldum=0;
  double ddum[1];
  const long int names[] = { CONTROL_PRINT_GID_GROUP,
    CONTROL_PRINT_GID_BEAM_VECTORS, CONTROL_PRINT_GID_TRUSS_VECTOR };

  // element_group is written by default (control_print_gid_group -yes)
  if ( print_option( PRINT_GID_GROUP, -1, -1, -YES )==-YES ) return 1;
  for ( i=0; i<3; i++ ) {
    db_max_index( names[i], max_index, VERSION_NORMAL, GET );
    for ( icontrol=0; icontrol<=max_index; icontrol++ ) {
      value = 0;
      db( names[i], icontrol, &value, ddum, ldum, VERSION_NORMAL,
        GET_IF_EXISTS );
      if ( value==-YES ) return 1;
    }
  }
  return 0;
}

// which element types have elements to write (types[itype] set)
static void print_gid_types_present( long int max_element, long int present[] )

{
  long int itype=0, element=0, element_group=0, length=0, ldum=0,
    *el=NULL;
  double ddum[1];

  array_set( present, 0, gid_ntype );
  el = get_new_int( MAXIMUM_NODE+1 );
  for ( element=0; element<=max_element; element++ ) {
    if ( !db_active_index( ELEMENT, element, VERSION_PRINT ) ) continue;
    db( ELEMENT, element, el, ddum, length, VERSION_PRINT, GET );
    element_group = 0;
    db( ELEMENT_GROUP, element, &element_group, ddum, ldum, VERSION_PRINT,
      GET_IF_EXISTS );
    if ( !print_gid_element_selected( element, el[0], element_group,
         length-1, &el[1] ) ) continue;
    for ( itype=0; itype<gid_ntype; itype++ )
      if ( el[0]==print_gid_type(itype) ) present[itype] = 1;
  }
  // spring2 as a Point mesh has no Gauss points
  if ( gid_spring2!=2 ) present[gid_ntype-1] = 0;
  delete[] el;
}

// Gauss points of the element results: "element_<type>", 1 internal point
// for the isoparametric types, 2 points with the nodes included for the
// line types (as the Professional 02-08-2026)
void print_gid_gauss_points( ofstream &outres, long int max_element )

{
  long int itype=0, name=0, present[gid_ntype];

  if ( !print_gid_gauss_needed() ) return;
  print_gid_types_present( max_element, present );
  for ( itype=0; itype<gid_ntype; itype++ ) {
    if ( !present[itype] ) continue;
    name = print_gid_type( itype );
    outres << "GaussPoints \"element_" << db_name(-name) << "\" ElemType "
           << print_gid_type_gid_name(name) << endl;
    if ( print_gid_type_linear(name) ) {
      outres << "Nodes included" << endl;
      outres << "Number Of Gauss Points: 2" << endl;
    }
    else
      outres << "Number Of Gauss Points: 1" << endl;
    outres << "Natural Coordinates: Internal" << endl;
    outres << "End GaussPoints" << endl;
  }
}

// control_print_gid_coord/print_gid_coord -yes (6.298/6.982): nodal result
// "coord" with the stored node coordinates
void print_gid_coord_result( ofstream &outres, double step_val,
  long int max_node )

{
  long int inod=0, idim=0, ldum=0, idum[1];
  double coord[MDIM];

  if ( gid_coord!=-YES ) return;
  outres << "Result  \"coord\"     \"Load Analysis\" \t" << step_val
         << " Vector OnNodes" << endl;
  outres << "Values" << endl;
  for ( inod=1; inod<=max_node; inod++ ) {
    if ( !db_active_index( NODE, inod, VERSION_PRINT ) ) continue;
    db( NODE, inod, idum, coord, ldum, VERSION_PRINT, GET );
    outres << inod << "\t";
    for ( idim=0; idim<ndim; idim++ ) outres << coord[idim] << "\t";
    outres << endl;
  }
  outres << "End Values" << endl;
}

// unit element axis t (node 1 to node 2) cross the plotting normal n
static void print_gid_axis_cross_normal( long int nodes[], double normal[],
  double result[] )

{
  long int idim=0, ldum=0, idum[1];
  double c0[MDIM], c1[MDIM], t[MDIM], n[MDIM], size=0.;

  array_set( c0, 0., MDIM ); array_set( c1, 0., MDIM );
  db( NODE, nodes[0], idum, c0, ldum, VERSION_PRINT, GET );
  db( NODE, nodes[1], idum, c1, ldum, VERSION_PRINT, GET );
  for ( idim=0; idim<MDIM; idim++ ) {
    t[idim] = idim<ndim ? c1[idim]-c0[idim] : 0.;
    n[idim] = normal[idim];
  }
  size = array_size( t, MDIM );
  if ( size>0. ) for ( idim=0; idim<MDIM; idim++ ) t[idim] /= size;
  size = array_size( n, MDIM );
  if ( size>0. ) for ( idim=0; idim<MDIM; idim++ ) n[idim] /= size;
  result[0] = t[1]*n[2] - t[2]*n[1];
  result[1] = t[2]*n[0] - t[0]*n[2];
  result[2] = t[0]*n[1] - t[1]*n[0];
}

// one element result of type itype on its Gauss points. kind 0:
// element_group; 1: truss_vector; 2..7: beam_vectors_0..5
static void print_gid_element_result( ofstream &outres, double step_val,
  long int max_element, long int itype, long int kind )

{
  long int element=0, element_group=0, length=0, ldum=0, idim=0, ipoint=0,
    npoint=0, name=print_gid_type(itype), any=0, idum[1], *el=NULL;
  double ddum[1], axis[MDIM], value[2], force=0., force_moment[12];
  char result_name[MCHAR];

  if      ( kind==0 ) strcpy( result_name, "element_group" );
  else if ( kind==1 ) strcpy( result_name, "truss_vector" );
  else snprintf( result_name, MCHAR, "beam_vectors_%ld", kind-2 );
  npoint = print_gid_type_linear(name) ? 2 : 1;
  el = get_new_int( MAXIMUM_NODE+1 );
  for ( element=0; element<=max_element; element++ ) {
    if ( !db_active_index( ELEMENT, element, VERSION_PRINT ) ) continue;
    db( ELEMENT, element, el, ddum, length, VERSION_PRINT, GET );
    if ( el[0]!=name ) continue;
    element_group = 0;
    db( ELEMENT_GROUP, element, &element_group, ddum, ldum, VERSION_PRINT,
      GET_IF_EXISTS );
    if ( !print_gid_element_selected( element, el[0], element_group,
         length-1, &el[1] ) ) continue;
    if ( !any ) {
      outres << "Result  \"" << result_name << "\"     \"Load Analysis\" \t"
             << step_val << ( kind==0 ? " Scalar" : " Vector" )
             << " OnGaussPoints \"element_" << db_name(-name) << "\"" << endl;
      outres << "Values" << endl;
      any = 1;
    }
    if ( kind==0 ) {
      for ( ipoint=0; ipoint<npoint; ipoint++ ) {
        if ( ipoint==0 ) outres << element << "\t";
        outres << element_group << endl;
      }
      continue;
    }
    if ( kind==1 ) {
      force = 0.;
      db( ELEMENT_TRUSS_FORCE, element, idum, &force, ldum, VERSION_PRINT,
        GET_IF_EXISTS );
      value[0] = value[1] = force;
      print_gid_axis_cross_normal( &el[1], gid_truss_normal, axis );
    }
    else {
      array_set( force_moment, 0., 12 );
      db( ELEMENT_BEAM_FORCE_MOMENT, element, idum, force_moment, ldum,
        VERSION_PRINT, GET_IF_EXISTS );
      value[0] = -force_moment[kind-2];
      value[1] = force_moment[6+kind-2];
      print_gid_axis_cross_normal( &el[1], gid_beam_normal, axis );
    }
    for ( ipoint=0; ipoint<npoint; ipoint++ ) {
      if ( ipoint==0 ) outres << element << "\t";
      for ( idim=0; idim<ndim; idim++ )
        outres << value[ipoint]*axis[idim] << "\t";
      outres << endl;
    }
  }
  if ( any ) outres << "End Values" << endl;
  delete[] el;
}

// element results on the Gauss points: element_group (control_print_gid_group
// /print_gid_group -yes, dropped by control_print_gid_other -no),
// truss_vector (control_print_gid_truss_vector -yes; -truss and
// -truss_beam) and beam_vectors_0..5 (control_print_gid_beam_vectors -yes;
// -beam and -truss_beam)
void print_gid_element_results( ofstream &outres, double step_val,
  long int max_element )

{
  long int itype=0, name=0, kind=0, present[gid_ntype];

  if ( !( gid_group==-YES && gid_other!=-NO ) &&
       gid_truss_vector!=-YES && gid_beam_vectors!=-YES ) return;
  print_gid_types_present( max_element, present );
  for ( itype=0; itype<gid_ntype; itype++ ) {
    if ( !present[itype] ) continue;
    name = print_gid_type( itype );
    if ( gid_group==-YES && gid_other!=-NO )
      print_gid_element_result( outres, step_val, max_element, itype, 0 );
    if ( gid_truss_vector==-YES && ( name==-TRUSS || name==-TRUSSBEAM ) )
      print_gid_element_result( outres, step_val, max_element, itype, 1 );
    if ( gid_beam_vectors==-YES && ( name==-BEAM || name==-TRUSSBEAM ) ) {
      for ( kind=2; kind<8; kind++ )
        print_gid_element_result( outres, step_val, max_element, itype,
          kind );
    }
  }
}
