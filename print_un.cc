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

void print_unknowns( void )

{
  long int icontrol=0, ipuknwn=0, iuknwn=0, inod=0, max_node=0, ldum=0, 
    idim=0, ready=0, icalcul=0, ncalcul=0, test=0, length_print_filter_index=0,
    print_filter_index[DATA_ITEM_SIZE], *dof_label=NULL;
  double ddum[1], *node_dof=NULL, *node_dof_calcul=NULL, *coord=NULL;
  char filename[MCHAR], str[MCHAR];

  set_swit(-1,-1,"print_unknowns");

  dof_label = get_new_int(DATA_ITEM_SIZE);

  db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET );
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );

  if ( !db( CONTROL_PRINT_FILTER, icontrol, print_filter_index, ddum, 
      length_print_filter_index, VERSION_NORMAL, GET_IF_EXISTS ) ) {
    print_filter_index[0] = -ALL; length_print_filter_index = 1;
  }

  db_max_index( NODE, max_node, VERSION_NORMAL, GET );
  for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
    iuknwn = ipuknwn*nder;
    test = filter( print_filter_index, length_print_filter_index,
      -NODE_DOF, inod, iuknwn, CHECK_NUMBER );
    if ( test ) {
      strcpy( filename, db_name(dof_label[iuknwn]) );
      long_to_a( icontrol, str );
      strcat( filename, "." );
      strcat( filename, str );
      ofstream out( filename );
      out.precision(TN_PRECISION);
      for ( inod=0; inod<=max_node; inod++ ) {
        if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
          test = filter( print_filter_index, length_print_filter_index,
            -NODE_DOF, inod, iuknwn, CHECK_INDEX );
          if ( test ) {
            coord = db_dbl( NODE, inod, VERSION_NORMAL );     
            node_dof = db_dbl( NODE_DOF, inod, VERSION_NORMAL );     
            for ( idim=0; idim<ndim; idim++ ) out << coord[idim] << " ";
            out << node_dof[iuknwn] << "\n";
          }
        }
      }
      out.close();
    }
  }

  if ( db_active_index( POST_CALCUL, 0, VERSION_NORMAL ) ) {
    icalcul = ready = 0;
    ncalcul = db_len( POST_CALCUL_SCAL_VEC_MAT, 0, VERSION_NORMAL );
    while ( !ready ) {
      test = filter( print_filter_index, length_print_filter_index,
        -NODE_DOF_CALCUL, inod, icalcul, CHECK_NUMBER );
      if ( test ) {
        strcpy( filename, post_calcul_names_without_extension[icalcul] );
        long_to_a( icontrol, str );
        strcat( filename, "." );
        strcat( filename, str );
        ofstream out( filename );
        out.precision(TN_PRECISION);
        for ( inod=0; inod<=max_node; inod++ ) {
          if ( db_active_index( NODE_DOF_CALCUL, inod, VERSION_NORMAL ) ) {
            test = filter( print_filter_index, length_print_filter_index,
              -NODE_DOF_CALCUL, inod, icalcul, CHECK_INDEX );
            if ( test ) {
              node_dof_calcul = db_dbl( NODE_DOF_CALCUL, inod, VERSION_NORMAL );
              coord = db_dbl( NODE, inod, VERSION_NORMAL );     
              for ( idim=0; idim<ndim; idim++ ) out << coord[idim] << " ";
              out << node_dof_calcul[icalcul] << "\n";
            }
          }
        }
        out.close();
      }
      icalcul++;
      ready = (icalcul>=ncalcul);
    }
  }

  delete[] dof_label;

}

void print_unknownsrhside( void )

{
  long int icontrol=0, ipuknwn=0, iuknwn=0, inod=0, max_node=0, ldum=0, idim=0, 
    test=0, length_print_filter_index=0, print_filter_index[DATA_ITEM_SIZE], 
    *dof_label=NULL;
  double ddum[1], *node_rhside=NULL, *coord=NULL;
  char filename[MCHAR], str[MCHAR];

  set_swit(-1,-1,"print_unknownsrhside");

  dof_label = get_new_int(DATA_ITEM_SIZE);

  db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET );
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );

  if ( !db( CONTROL_PRINT_FILTER, icontrol, print_filter_index, ddum, 
      length_print_filter_index, VERSION_NORMAL, GET_IF_EXISTS ) ) {
    print_filter_index[0] = -ALL; length_print_filter_index = 1;
  }

  db_max_index( NODE, max_node, VERSION_NORMAL, GET );
  for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
    iuknwn = ipuknwn*nder;
    test = filter( print_filter_index, length_print_filter_index,
      -NODE_RHSIDE, inod, iuknwn, CHECK_NUMBER );
    if ( test ) {
      strcpy( filename, db_name(dof_label[iuknwn]) );
      strcat( filename, "_rhside" );
      long_to_a( icontrol, str );
      strcat( filename, "." );
      strcat( filename, str );
      ofstream out( filename );
      out.precision(TN_PRECISION);
      for ( inod=0; inod<=max_node; inod++ ) {
        if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
          test = filter( print_filter_index, length_print_filter_index,
            -NODE_RHSIDE, inod, iuknwn, CHECK_INDEX );
          if ( test ) {
            coord = db_dbl( NODE, inod, VERSION_NORMAL );     
            node_rhside = db_dbl( NODE_RHSIDE, inod, VERSION_NORMAL );     
            for ( idim=0; idim<ndim; idim++ ) out << coord[idim] << " ";
            out << node_rhside[iuknwn] << "\n";
          }
        }
      }
      out.close();
    }
  }

  delete[] dof_label;

}
