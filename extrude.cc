/*
    Copyright (C) 2000  Dennis Roddeman
    email: d.g.roddeman@wb.utwente.nl

*/

#include "tochnog.h"

void extrude( void )

{
  // control_mesh_extrude runs once, when its own control index is
  // processed, also an index without control_timestep (Professional
  // 02-08-2026 other/extrude1 has no timestep at all). Before F5 T4 U5a
  // the first record of any index was applied at the first timestep.
  static long int extruded = 0;
  long int icontrol=0, ldum=0;
  double ddum[1];

  if ( extruded ) return;
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  if ( !db_active_index( CONTROL_MESH_EXTRUDE, icontrol, VERSION_NORMAL ) )
    return;
  extruded = 1;
  // mesh.cc mesh_extrude reads control_mesh_extrude and its options
  // (_n, _direction, _element_group_new, _contact_spring_element_group
  // and _new) at that control index
  mesh_extrude( icontrol );
}
