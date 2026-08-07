/*
    Copyright (C) 2000  Dennis Roddeman
    email: d.g.roddeman@wb.utwente.nl

*/

#include "tochnog.h"

void extrude( void )

{
  long int icontrol=0, ldum=0, ext_length=0, n_layer=0, idum_e[1];
  double ddum[1], ext_z[DATA_ITEM_SIZE];

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  if ( db_active_index( CONTROL_MESH_EXTRUDE, icontrol, VERSION_NORMAL ) ) {
    db( CONTROL_MESH_EXTRUDE, icontrol, idum_e, ext_z, ext_length,
      VERSION_NORMAL, GET );
    n_layer = ext_length;
    if ( n_layer>0 ) mesh_extrude( ext_z, n_layer );
  }
}
