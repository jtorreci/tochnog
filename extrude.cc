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
  // control_mesh_extrude may be declared with ANY control index (the
  // corpus uses control_mesh_extrude 10 with control_timestep 30), so
  // look it up over the whole control range like the convert.
  {
    long int ic_max = 0, ic_found = -1;
    db_max_index( CONTROL_MESH_EXTRUDE, ic_max, VERSION_NORMAL, GET );
    for ( long int ic2=0; ic2<=ic_max; ic2++ ) {
      if ( db_active_index( CONTROL_MESH_EXTRUDE, ic2, VERSION_NORMAL ) ) {
        ic_found = ic2; break;
      }
    }
    if ( ic_found>=0 ) icontrol = ic_found;
  }
  if ( db_active_index( CONTROL_MESH_EXTRUDE, icontrol, VERSION_NORMAL ) ) {
    db( CONTROL_MESH_EXTRUDE, icontrol, idum_e, ext_z, ext_length,
      VERSION_NORMAL, GET );
    n_layer = ext_length;
    // control_mesh_extrude_n: the FIRST value is the number of layers
    // of the quadratic extrusion (quad9 -> hex27, the Professional's
    // force11 ring: control_mesh_extrude 10 0. 10. + _n 10 2 1 -> 2
    // layers of 5 over the 0..10 extent). The linear branches
    // (tria3/bar2/quad4) keep the layer-boundary convention of the
    // record values themselves. The _n record is read into a SEPARATE
    // buffer (ext_z still holds the extrude record for mesh_extrude).
    long int quad9_layers = 0, ext_length_n = 0, idum_n[DATA_ITEM_SIZE];
    if ( db_active_index( CONTROL_MESH_EXTRUDE_N, icontrol,
         VERSION_NORMAL ) ) {
      db( CONTROL_MESH_EXTRUDE_N, icontrol, idum_n, ddum,
        ext_length_n, VERSION_NORMAL, GET );
      if ( ext_length_n>0 ) quad9_layers = idum_n[0];
    }
    if ( n_layer>0 ) mesh_extrude( ext_z, n_layer, quad9_layers );
  }
}
