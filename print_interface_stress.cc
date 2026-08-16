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
    *nodes=NULL, in_geometry=0;
  double ddum[1], xstart=0., ystart=0., xend=0., yend=0.,
    dx=0., dy=0., len=0., strain_normal=0., sign=0., sigt=0., kn=0.,
    f_t=0., f_t2=0., ddum3[3], centroid[MDIM];
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
  if ( ndim==2 ) {
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
  }

  el = get_new_int(MAXIMUM_NODE+1);
  nodes = get_new_int(MAXIMUM_NODE);
  db_highest_index( ELEMENT, max_element, VERSION_NORMAL );

  // 3D filter: only interface elements whose middle lies in the geometry
  // (control_print_interface_stress_3d_geometry name index) are printed.
  long int geometry_entity[DATA_ITEM_SIZE], geometry_length=0;
  array_set( geometry_entity, 0., DATA_ITEM_SIZE );
  if ( ndim==3 && db_active_index( CONTROL_PRINT_INTERFACE_STRESS_3D_GEOMETRY,
      icontrol, VERSION_NORMAL ) ) {
    db( CONTROL_PRINT_INTERFACE_STRESS_3D_GEOMETRY, icontrol, geometry_entity,
      ddum, geometry_length, VERSION_NORMAL, GET );
  }
  long int order=0, order_axis=-1;
  if ( ndim==3 ) {
    db( CONTROL_PRINT_INTERFACE_STRESS_3D_ORDER, icontrol, &order,
      ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    if ( order==-X ) order_axis = 0;
    else if ( order==-Y ) order_axis = 1;
    else if ( order==-Z ) order_axis = 2;
  }

  // 3D output is collected first so it can be sorted (control_print_interface_
  // stress_3d_order) or written in element order. Each entry is the sort key
  // (element number or the requested coordinate) + the line.
  const long int MAX_3D_LINES = 4096;
  long int n3d=0;
  double *key3d = ( order_axis>=0 ) ? new double[MAX_3D_LINES] : NULL;
  long int *elem3d = ( order_axis<0 ) ? new long int[MAX_3D_LINES] : NULL;
  char **line3d = new char*[MAX_3D_LINES];
  for ( long int i=0; i<MAX_3D_LINES; i++ ) line3d[i] = NULL;

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
    // semantics as sign -> consistent sigt. 3D: two tangential components.
    f_t = 0.; f_t2 = 0.;
    db( ELEMENT_INTERFACE_FORCE_TANG, element, idum, &f_t, ldum,
      VERSION_NORMAL, GET_IF_EXISTS );
    if ( ndim==3 )
      db( ELEMENT_INTERFACE_FORCE_TANG2, element, idum, &f_t2, ldum,
        VERSION_NORMAL, GET_IF_EXISTS );
    sigt = f_t;

    if ( ndim==2 ) {
      for ( inol=0; inol<nnol; inol++ ) {
        inod = nodes[inol];
        double *c = db_dbl( NODE, inod, VERSION_NORMAL );
        // distance from the start point, projected on the cut direction
        double proj = ( (c[0]-xstart)*dx + (c[1]-ystart)*dy ) / len;
        out << proj << " " << sign << " " << sigt << "\n";
      }
    }
    else {
      // 3D: element middle + average stresses (interface_sign, interface_sigt1,
      // interface_sigt2). A line is written for each interface element.
      array_set( centroid, 0., MDIM );
      for ( inol=0; inol<nnol; inol++ ) {
        double *c = db_dbl( NODE, nodes[inol], VERSION_NORMAL );
        for ( long int idim=0; idim<3; idim++ ) centroid[idim] += c[idim];
      }
      for ( long int idim=0; idim<3; idim++ ) centroid[idim] /= nnol;
      // geometry filter: at least one interface-element node must lie in the
      // geometry (control_print_interface_stress_3d_geometry)
      if ( geometry_length>0 ) {
        in_geometry = 0;
        for ( inol=0; inol<nnol && !in_geometry; inol++ ) {
          geometry( nodes[inol], ddum, geometry_entity, in_geometry, ddum[0],
            ddum, ddum[0], ddum, NODE_START_REFINED, PROJECT_EXACT,
            VERSION_NORMAL );
        }
        if ( !in_geometry ) continue;
      }
      if ( n3d>=MAX_3D_LINES ) continue;
      char *line = new char[MCHAR];
      snprintf( line, MCHAR, "%g %g %g %g %g %g",
        centroid[0], centroid[1], centroid[2], sign, sigt, f_t2 );
      if ( order_axis>=0 ) key3d[n3d] = centroid[order_axis];
      else                 elem3d[n3d] = element;
      line3d[n3d] = line;
      n3d++;
    }
  }

  // write 3D lines: sorted (control_print_interface_stress_3d_order) or in
  // element order (default)
  if ( ndim==3 ) {
    long int idx[MAX_3D_LINES];
    for ( long int i=0; i<n3d; i++ ) idx[i] = i;
    for ( long int i=1; i<n3d; i++ ) {
      long int j = i;
      while ( j>0 ) {
        long int a = idx[j-1], b = idx[j];
        long int swap = 0;
        if ( order_axis>=0 ) {
          if ( key3d[b]<key3d[a] ) swap = 1;
        }
        else {
          if ( elem3d[b]<elem3d[a] ) swap = 1;
        }
        if ( !swap ) break;
        idx[j] = a; idx[j-1] = b; j--;
      }
    }
    for ( long int i=0; i<n3d; i++ ) out << line3d[idx[i]] << "\n";
  }
  for ( long int i=0; i<MAX_3D_LINES; i++ ) delete[] line3d[i];
  delete[] line3d;
  delete[] key3d;
  delete[] elem3d;

  out.close();
  delete[] el;
  delete[] nodes;

  if ( swit ) pri( "Out function PRINT_INTERFACE_STRESS" );
}
