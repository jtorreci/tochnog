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

void unknown_freeze( void )

{
  long int i=0, iuknwn=0, icontrol=0, length=0, swit=0, ldum=0,
    inod=0, max_node=0, control_unknown_freeze[MUKNWN], 
    dof_label[MUKNWN], *node_bounded=NULL;
  double ddum[1];

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  if ( db_active_index( CONTROL_UNKNOWN_FREEZE, icontrol, VERSION_NORMAL ) ) {
    swit = set_swit(-1,-1,"unknown_freeze");
    if ( swit ) pri( "In routine UNKNOWN_FREEZE" );
    db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( CONTROL_UNKNOWN_FREEZE, icontrol, control_unknown_freeze, 
      ddum, length, VERSION_NORMAL, GET );
    db_max_index( NODE, max_node, VERSION_NORMAL, GET );
    for ( inod=0; inod<=max_node; inod++ ) {
      if ( db_active_index( NODE_BOUNDED, inod, VERSION_NORMAL ) ) {         
        node_bounded = db_int( NODE_BOUNDED, inod, VERSION_NORMAL );
        for ( i=0; i<length; i++ ) {
          if ( control_unknown_freeze[i]<0 )
            array_member(dof_label,control_unknown_freeze[i],nuknwn,iuknwn);
          else
            iuknwn = control_unknown_freeze[i];             
          node_bounded[iuknwn] = 1;
        }
      }
    }
    if ( swit ) pri( "Out routine UNKNOWN_FREEZE" );
  }

}

void unknown_reset( void )

{
  long int i=0, iuknwn=0, icontrol=0, 
    length_control_unknown_reset_unknown=0, swit=0, ldum=0,
    inod=0, max_node=0, use_geometry=0, in_geometry=0, idum[1],
    control_unknown_reset_unknown[MUKNWN], control_unknown_reset_geometry[2], 
    dof_label[MUKNWN];
  double rdum=0., control_unknown_reset_value=0.,
    ddum[1], normal[MDIM], *node_dof=NULL;

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  if ( db_active_index( CONTROL_UNKNOWN_RESET_UNKNOWN, icontrol, VERSION_NORMAL ) ) {
    swit = set_swit(-1,-1,"unknown_reset");
    if ( swit ) pri( "In routine UNKNOWN_RESET" );
    db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( CONTROL_UNKNOWN_RESET_UNKNOWN, icontrol, control_unknown_reset_unknown, 
      ddum, length_control_unknown_reset_unknown, VERSION_NORMAL, GET );
    db( CONTROL_UNKNOWN_RESET_VALUE, icontrol, idum, &control_unknown_reset_value, 
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
    use_geometry = db( CONTROL_UNKNOWN_RESET_GEOMETRY, icontrol, 
      control_unknown_reset_geometry, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db_max_index( NODE, max_node, VERSION_NORMAL, GET );
    for ( inod=0; inod<=max_node; inod++ ) {
      if ( db_active_index( NODE_DOF, inod, VERSION_NORMAL ) ) {         
        if ( use_geometry )
          geometry( inod, ddum, control_unknown_reset_geometry, in_geometry, rdum, 
            normal, rdum, ddum, NODE_START_REFINED, PROJECT_EXACT, VERSION_NORMAL );
        else
          in_geometry = 1;
        if ( in_geometry ) {        
          node_dof = db_dbl( NODE_DOF, inod, VERSION_NORMAL );
          for ( i=0; i<length_control_unknown_reset_unknown; i++ ) {
            if ( control_unknown_reset_unknown[i] == -ALL ) {
              for ( iuknwn=0; iuknwn<nuknwn; iuknwn++ ) node_dof[iuknwn] = 0.;
            }
            else {
              if ( control_unknown_reset_unknown[i]<0 )
                array_member(dof_label,control_unknown_reset_unknown[i],nuknwn,iuknwn);
              else
                iuknwn = control_unknown_reset_unknown[i];             
              node_dof[iuknwn] = control_unknown_reset_value;
            }
          }
        }
      }
    }
    if ( swit ) pri( "Out routine UNKNOWN_RESET" );
  }

}
