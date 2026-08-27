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

// print_materi_stress_force - control_print_materi_stress_force
// (manual Professional 6.328): prints the forces and moments calculated
// by post_calcul -materi_stress -force to the special purpose ASCII
// file materi_stress_force.<index>. The "index" of the manual is the
// RECORD INDEX (icontrol), the same convention as
// control_print_beam_force_moment: "the name of the file will be
// materi_stress_force.100 if the index is 100" means the record
// control_print_materi_stress_force 100 -all -> materi_stress_force.100.
//
// Record value (INTEGER, fixed 1): the method - -all prints every
// result, -primary only the primarily calculated results (i.e. WITHOUT
// the averaged results of post_calcul_materi_stress_force_average -yes
// on the middle-plane nodes of quad9/hex27 elements).
//
// File structure (the header comments below explain it, as the manual
// requires: "the files themselves will contain comments explaining the
// detailed structure"): one line per node; the FIRST column is the node
// (position in the compacted print version, 0-based - the same layout
// the vtk post_calcul block of print_vt.cc reads), followed by the
// post_calcul -materi_stress -force items in the order of
// post_calcul_names:
//   2D (9 columns):  norx_sig nory_sig nors_sig shex_sig shey_sig
//                    shes_sig momx_sig momy_sig moms_sig
//   3D (16 columns): norx_sig nory_sig norz_sig nors_sig shex_sig
//                    shey_sig shez_sig shes_sig mom1x_sig mom1y_sig
//                    mom1z_sig mom1s_sig mom2x_sig mom2y_sig mom2z_sig
//                    mom2s_sig
// The x/y/z components are GLOBAL PLOT components (the vector is drawn
// in the structure thickness direction); the s component is the PHYSICAL
// vector size (the design value). The values are read from
// NODE_DOF_CALCUL (VERSION_PRINT, one slot per item - the first slot of
// the -force block is located by scanning POST_CALCUL_UNKNOWN_OPERAT).
//
// The file is opened in append mode (one block per print, pattern of
// print_beam_force_moment). LOT 1 (infrastructure): the numerical
// integration of the -force family is not implemented yet, so the values
// are 0; -primary and -all write the same lines until the averaging of
// L2/L3 marks non-primary nodes (the filter hook msf_node_is_averaged
// returns 0 for every node). If no -force post_calcul block exists no
// file is written (decision, like print_beam_force_moment without
// crossed elements).

// Filter hook for the -primary method: returns 1 when the node values
// are the AVERAGED results (quad9/hex27 middle-plane nodes with
// post_calcul_materi_stress_force_average -yes) that -primary must
// skip. LOT 1: no averaging exists yet, so every node is primary.
// The L2/L3 integration will flag the averaged nodes (e.g. a per-node
// record) and this hook will read that flag.
static long int msf_node_is_averaged( long int inod )

{
  return 0;
}

void print_materi_stress_force( long int icontrol, long int method )

{
  long int swit=0, idum[1], max_node=0,
    inod=0, i=0, nitems=0, offset=0, ncalcul=0, j=0, nforce=0,
    *post_calcul_unknown_operat=NULL;
  double ddum[1], *node_dof_calcul=NULL;
  char filename[MCHAR], str[MCHAR];

  swit = set_swit(-1,-1,"print_materi_stress_force");
  if ( swit ) pri( "In routine PRINT_MATERI_STRESS_FORCE" );

  // record value: the method (-all / -primary)
  if ( method!=-ALL && method!=-PRIMARY )
    db_error( CONTROL_PRINT_MATERI_STRESS_FORCE, icontrol );

  // locate the -materi_stress -force block in post_calcul: the first
  // FORCE item of POST_CALCUL_UNKNOWN_OPERAT; its slot in
  // NODE_DOF_CALCUL equals its item index (one slot per item, the flat
  // layout of print_vt.cc). No -force block -> no file (decision).
  if ( !db_active_index( POST_CALCUL, 0, VERSION_NORMAL ) ) {
    if ( swit ) pri( "Out function PRINT_MATERI_STRESS_FORCE (no post_calcul)" );
    return;
  }
  post_calcul_unknown_operat = get_new_int(DATA_ITEM_SIZE);
  db( POST_CALCUL_UNKNOWN_OPERAT, 0, post_calcul_unknown_operat, ddum,
    ncalcul, VERSION_NORMAL, GET );
  for ( i=0; i+1<ncalcul; i+=2 ) {
    if ( labs(post_calcul_unknown_operat[i+1])==FORCE ) break;
  }
  if ( i+1>=ncalcul ) {
    delete[] post_calcul_unknown_operat;
    if ( swit )
      pri( "Out function PRINT_MATERI_STRESS_FORCE (no -force post_calcul)" );
    return;
  }
  offset = i/2; // item index == slot index of the first FORCE item
  nitems = post_calcul_materi_stress_force_items();
  for ( j=i; j+1<ncalcul; j+=2 )
    if ( labs(post_calcul_unknown_operat[j+1])==FORCE ) nforce++;
  delete[] post_calcul_unknown_operat;
  if ( nforce<nitems ) {
    // incomplete -force block (e.g. cut by the MCALCUL=20 limit):
    // nothing reliable to print
    if ( swit )
      pri( "Out function PRINT_MATERI_STRESS_FORCE (incomplete -force block)" );
    return;
  }

  // compacted print version of the database (pattern of print_vt.cc:
  // db_version_copy + renumbering + highest_index + read + delete)
  db_version_copy( VERSION_NORMAL, VERSION_PRINT );
  renumbering( VERSION_PRINT, NO, 0, 0, idum, idum );
  db_highest_index( NODE, max_node, VERSION_PRINT );
  if ( max_node<0 ) {
    db_version_delete( VERSION_PRINT );
    if ( swit ) pri( "Out function PRINT_MATERI_STRESS_FORCE (no nodes)" );
    return;
  }

  strcpy( filename, "materi_stress_force." );
  long_to_a( icontrol, str ); // the manual "index" = record index
  strcat( filename, str );
  ofstream out( filename, ios::app );
  out.precision(TN_PRECISION);
  out << "# materi_stress_force." << icontrol
      << " - normal force, shear force and moment(s) per node\n";
  out << "# (post_calcul -materi_stress -force; manual Professional "
         "6.328/6.913)\n";
  out << "# One line per node. Columns: node";
  for ( i=0; i<nitems; i++ ) out << " " << post_calcul_names[offset+i];
  out << "\n";
  out << "# x/y/z components: global plot vector (drawn in the structure "
         "thickness direction);\n";
  out << "# s component: physical vector size (the design value). Method: ";
  if ( method==-ALL ) out << "-all (all results)";
  else                out << "-primary (only the primarily calculated "
                             "results, without the quad9/hex27 averaged ones)";
  out << "\n";

  for ( inod=0; inod<=max_node; inod++ ) {
    if ( method==-PRIMARY && msf_node_is_averaged( inod ) ) continue;
    node_dof_calcul = db_dbl( NODE_DOF_CALCUL, inod, VERSION_PRINT );
    out << inod;
    for ( i=0; i<nitems; i++ ) out << " " << node_dof_calcul[offset+i];
    out << "\n";
  }
  out.close();

  db_version_delete( VERSION_PRINT );

  if ( swit ) pri( "Out function PRINT_MATERI_STRESS_FORCE" );
}
