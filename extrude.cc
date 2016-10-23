/*
    Copyright (C) 2000  Dennis Roddeman
    email: d.g.roddeman@wb.utwente.nl

*/

#include "tochnog.h"

void extrude( void )

{
  long int icontrol=0, ldum=0;
  double ddum[1];

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  if ( db_active_index( CONTROL_MESH_EXTRUDE, icontrol, VERSION_NORMAL ) ) {
    pri( "Error: extrude is not available." );
    exit(TN_EXIT_STATUS);
  }
}
