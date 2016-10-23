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

void restart( void )

{
  long int icontrol=0, control_restart=0, length=0, swit=0,
    ldum=0, idum[1];
  double start_time_current=0., ddum[1];

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  if ( db_active_index( CONTROL_RESTART, icontrol, VERSION_NORMAL ) ) {
    swit = set_swit(-1,-1,"restart");
    if ( swit ) pri( "In routine RESTART" );
    db( CONTROL_RESTART, icontrol, &control_restart, ddum, ldum, VERSION_NORMAL, GET );
    if ( control_restart==-YES ) {
      db( TIME_CURRENT, 0, idum, &start_time_current, length, VERSION_START, GET );
      db( TIME_CURRENT, 0, idum, &start_time_current, length, VERSION_NORMAL, PUT );
      db_copy( NODE_START_REFINED, NODE, VERSION_NORMAL );
      db_copy( NODE_DOF_START_REFINED, NODE_DOF, VERSION_NORMAL );
    }
    if ( swit ) pri( "Out routine RESTART" );
  }

}
