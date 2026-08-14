// print_interface_stress - control_print_interface_stress (Carril A, Fase 4).
// 2D: prints the interface stresses (interface_sign, interface_sigt)
// through a set of interfaces starting at (xstart,ystart) and ending at
// (xend,yend) as specified by control_print_interface_stress_2d_coordinates.
// The switch must be -separate_index or -separate_sequential. The stresses
// are written to the file interface_stress.<index>. The first column is

#include "tochnog.h"
// the distance from the start point; the following columns are
// interface_sign and interface_sigt. A line is written for each node of
// each interface element.
void print_interface_stress( long int icontrol, long int task )

{
  long int element=0, max_element=0, inod=0, inol=0, max_node=0,
    length=0, ldum=0, swit=0, element_group=0, idum[1], *el=NULL,
    *nodes=NULL;
  double ddum[1], xstart=0., ystart=0., xend=0., yend=0.,
    dx=0., dy=0., len=0., strain_normal=0., sign=0., sigt=0., kn=0.,
    f_t=0., ddum3[3];
  char filename[MCHAR], str[MCHAR];

  swit = set_swit(-1,-1,"print_interface_stress");
  if ( swit ) pri( "In routine PRINT_INTERFACE_STRESS." );

  // file name: interface_stress.<index> (task -separate_index /
  // -separate_sequential)
  strcpy( filename, "interface_stress." );
  if ( task==-SEPARATE_SEQUENTIAL ) {
    static long int iseq=0;
    long_to_a( iseq++, str );
    strcat( filename, str );
  }
  else {
    long_to_a( icontrol, str );
    strcat( filename, str );
  }
  ofstream out( filename, ios::app );
  out.precision(TN_PRECISION);

  // 2D cut line
  if ( db_active_index( CONTROL_PRINT_INTERFACE_STRESS_2D_COORDINATES,
      icontrol, VERSION_NORMAL ) ) {
    db( CONTROL_PRINT_INTERFACE_STRESS_2D_COORDINATES, icontrol, idum,
      ddum, length, VERSION_NORMAL, GET );
    xstart = ddum[0]; ystart = ddum[1]; xend = ddum[2]; yend = ddum[3];
  }
  else {
    // default: cut through the whole mesh along x at y=0
    db_highest_index( NODE, max_node, VERSION_NORMAL );
    xstart = 0.; ystart = 0.; xend = 1.; yend = 0.;
  }
  dx = xend - xstart;
  dy = yend - ystart;
  len = sqrt( dx*dx + dy*dy );

  el = get_new_int(MAXIMUM_NODE+1);
  nodes = get_new_int(MAXIMUM_NODE);
  db_highest_index( ELEMENT, max_element, VERSION_NORMAL );

  for ( element=0; element<=max_element; element++ ) {
    if ( !db_active_index( ELEMENT, element, VERSION_NORMAL ) ) continue;
    db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
    element_group = 0;
    db( ELEMENT_GROUP, element, &element_group, ddum, ldum,
      VERSION_NORMAL, GET_IF_EXISTS );
    if ( !db_active_index( GROUP_INTERFACE, element_group, VERSION_NORMAL ) )
      continue;

    long int nnol = length - 1;
    array_move( &el[1], nodes, nnol );

    // stiffness of this interface group
    db( GROUP_INTERFACE_MATERI_ELASTI_STIFFNESS, element_group, idum,
      ddum3, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    kn = ddum3[0];

    // accumulated normal strain -> normal stress (total normal force)
    strain_normal = 0.;
    db( ELEMENT_INTERFACE_STRAIN_NORMAL, element, idum, &strain_normal,
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
    sign = kn * strain_normal;
    // accumulated total tangential force (history ELEMENT_INTERFACE_FORCE_TANG,
    // Fase 3): cumulative Mohr-Coulomb or elastic total. Same accumulated
    // semantics as sign -> consistent sigt.
    f_t = 0.;
    db( ELEMENT_INTERFACE_FORCE_TANG, element, idum, &f_t, ldum,
      VERSION_NORMAL, GET_IF_EXISTS );
    sigt = f_t;

    for ( inol=0; inol<nnol; inol++ ) {
      inod = nodes[inol];
      double *c = db_dbl( NODE, inod, VERSION_NORMAL );
      // distance from the start point, projected on the cut direction
      double proj = ( (c[0]-xstart)*dx + (c[1]-ystart)*dy ) / len;
      out << proj << " " << sign << " " << sigt << "\n";
    }
  }

  out.close();
  delete[] el;
  delete[] nodes;

  if ( swit ) pri( "Out function PRINT_INTERFACE_STRESS" );
}
