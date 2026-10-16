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
#include "tochnog_exceptions.h"
#include "string_utils.h"
#include <sys/resource.h>

void element_middle_radius_set( )

{
  long int element=0, max_element=0, idim=0, inol=0, inod=0, 
    nnol=0, length=0, ldum=0, idum[1], el[1+MNOL], nodes[MNOL];
  double distance=0., radius=0., middle[MDIM], coord[MDIM], 
    work[MDIM], ddum[1], *node_dof=NULL;

  db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
  for ( element=0; element<=max_element; element++ ) {
    if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) ) {
      db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
      nnol = length - 1; array_move( &el[1], nodes, nnol );
      array_set( middle, 0., MDIM );
      for ( inol=0; inol<nnol; inol++ ) {
        inod = nodes[inol];
        db( NODE, inod, idum, coord, ldum, VERSION_NORMAL, GET );
        if ( materi_displacement ) {
          node_dof = db_dbl( NODE_DOF, inod, VERSION_NORMAL );
          for ( idim=0; idim<ndim; idim++ ) {
            coord[idim] += node_dof[dis_indx+idim*nder];
          }
        }
        array_add( coord, middle, middle, ndim );
      }
      array_multiply( middle, middle, (double) 1/nnol, ndim );
      db( ELEMENT_MIDDLE, element, idum, middle, ndim, VERSION_NORMAL, PUT );
      radius = 0.;
      for ( inol=0; inol<nnol; inol++ ) {
        inod = nodes[inol];
        db( NODE, inod, idum, coord, ldum, VERSION_NORMAL, GET );
        distance = array_distance( coord, middle, work, ndim );
        if ( distance>radius ) radius = distance;
      }
      length = 1;
      db( ELEMENT_RADIUS, element, idum, &radius, length, VERSION_NORMAL, PUT );
    }
  }
}

void element_volume_set( long int name, long int nodes[], long int version, 
  double &element_volume )

{

  long int inod=0, jnod=0, knod=0, lnod=0, ldum=0, idum[1];
  double coord0[MDIM], coord1[MDIM], coord2[MDIM], coord3[MDIM], work[MDIM];

  if      ( ndim==1 ) {
    inod = nodes[0];
    jnod = nodes[1];
    db( NODE, inod, idum, coord0, ldum, version, GET );
    db( NODE, jnod, idum, coord1, ldum, version, GET );
    element_volume = array_distance( coord0, coord1, work, ndim );
  }
  else if ( ndim==2 ) {
    if ( name==-TRIA3 ) {
      inod = nodes[0];
      jnod = nodes[1];
      knod = nodes[2];
    }
    else {
      assert( name==-TRIA6 );
      inod = nodes[0];
      jnod = nodes[2];
      knod = nodes[5];
    }
    db( NODE, inod, idum, coord0, ldum, version, GET );
    db( NODE, jnod, idum, coord1, ldum, version, GET );
    db( NODE, knod, idum, coord2, ldum, version, GET );
    element_volume = triangle_area( coord0, coord1, coord2 );
  }
  else {
    assert( ndim==3 );
    if ( name==-TET4 ) {
      inod = nodes[0];
      jnod = nodes[1];
      knod = nodes[2];
      lnod = nodes[3];
    }
    else {
      assert( name==-TET10 );
      inod = nodes[0];
      jnod = nodes[2];
      knod = nodes[5];
      lnod = nodes[9];
    }
    db( NODE, inod, idum, coord0, ldum, version, GET );
    db( NODE, jnod, idum, coord1, ldum, version, GET );
    db( NODE, knod, idum, coord2, ldum, version, GET );
    db( NODE, lnod, idum, coord3, ldum, version, GET );
    element_volume = tetrahedron_volume( coord0, coord1, coord2, coord3 );
  }

}

void get_element_matrix_unknowns( long int element,
  long int element_matrix_unknowns[] )

{
  long int i=0, n=0, inol=0, jnol=0, inod=0, jnod=0,
    ipuknwn=0, jpuknwn=0, length=0, nnol=0, element_group=0, ldum=0,
    el[1+MNOL], nodes[MNOL], *group_matrix_unknowns=NULL;
  double ddum[1];

  if ( db_active_index( ELEMENT_MATRIX_UNKNOWNS, element, VERSION_NORMAL ) ) 
    db( ELEMENT_MATRIX_UNKNOWNS, element, element_matrix_unknowns, ddum, ldum,
      VERSION_NORMAL, GET );
  else {
    db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
    nnol = length - 1; array_move( &el[1], nodes, nnol );
    db( ELEMENT_GROUP, element, &element_group, ddum, ldum, 
      VERSION_NORMAL, GET_IF_EXISTS );
    group_matrix_unknowns = 
      db_int( GROUP_MATRIX_UNKNOWNS, element_group, VERSION_NORMAL ); 
    length = db_len( GROUP_MATRIX_UNKNOWNS, element_group, VERSION_NORMAL );
    n = length / 4;
    for ( i=0; i<n; i++ ) {
      inol = group_matrix_unknowns[4*i+0]; inod = nodes[inol];
      ipuknwn = group_matrix_unknowns[4*i+1];
      jnol = group_matrix_unknowns[4*i+2]; jnod = nodes[jnol];
      jpuknwn = group_matrix_unknowns[4*i+3];
      element_matrix_unknowns[2*i+0] = inod*npuknwn + ipuknwn;
      element_matrix_unknowns[2*i+1] = jnod*npuknwn + jpuknwn;
    }
  }

}

char *get_new_char( long int n )

{
  char *ptr=NULL;

  if ( n<=0 ) n = 1;
  if ( !(ptr = new char[n] ) ) {
    pri( "Error: cannot allocate enough memory." );
#ifdef USE_EXCEPTIONS
    throw OutOfMemoryException();
#else
    exit(TN_EXIT_STATUS);
#endif
  }
  return ptr;
}

double *get_new_dbl( long int n )

{
  double *ptr=NULL;

  if ( n<=0 ) n = 1;
  if ( !(ptr = new double[n] ) ) {
    pri( "Error: cannot allocate enough memory." );
#ifdef USE_EXCEPTIONS
    throw OutOfMemoryException();
#else
    exit(TN_EXIT_STATUS);
#endif
  }
  return ptr;
}

long int *get_new_int( long int n )

{
  long int *ptr=NULL;

  if ( n<=0 ) n = 1;
  if ( !(ptr = new long int[n] ) ) {
    pri( "Error: cannot allocate enough memory." );
#ifdef USE_EXCEPTIONS
    throw OutOfMemoryException();
#else
    exit(TN_EXIT_STATUS);
#endif
  }
  return ptr;
}

int *get_new_int_short( long int n )
 
{
  int *ptr=NULL;
 
  if ( n<=0 ) n = 1;
  if ( !(ptr = new int[n] ) ) {
    pri( "Error: cannot allocate enough memory." );
    exit(TN_EXIT_STATUS);
  }
  return ptr;
}
                   

void set_environment( void )

{
  long int length=0, options_processors=1,
    options_solver=-MATRIX_ITERATIVE_BICG;
  double ddum[1];
  char *str=NULL;

  if ( !db_active_index( OPTIONS_PROCESSORS, 0, VERSION_NORMAL ) ) {
    str = getenv("TOCHNOG_OPTIONS_PROCESSORS");
    if ( str!=NULL ) {
      options_processors = atoi( str );
      length = 1;
      db( OPTIONS_PROCESSORS, 0, &options_processors, ddum, 
        length, VERSION_NORMAL, PUT );
    }
  }
  if ( !db_active_index( OPTIONS_SOLVER, 0, VERSION_NORMAL ) ) {
    str = getenv("TOCHNOG_OPTIONS_SOLVER");
    if ( str!=NULL ) {
      options_solver = -db_number( str );
      length = 1;
      db( OPTIONS_SOLVER, 0, &options_solver, ddum, 
        length, VERSION_NORMAL, PUT );
    }
  }

}

long int set_swit( long int element, long int inod, const char* routine )

{
  long int i=0, result=0, iteration=0, icontrol=0, ldum=0;
  double ddum[1];
  char *str=NULL;

  if ( strcmp(routine,swit_routine_stack[0]) ) {
    for ( i=MSTACK-2; i>=0; i-- ) {
      strcpy( swit_routine_stack[i+1], swit_routine_stack[i] );
    }
    strcpy( swit_routine_stack[0], routine );
  }
  if ( element>=0 )
    swit_element_stack = element;
  else
    swit_element_stack = -1;
  if ( inod>=0 ) 
    swit_node_stack = inod;
  else
    swit_node_stack = -1;

  if ( getenv("TOCHNOG_REPORT_ROUTINE")!=NULL ) {
    pri( "In routine ", routine );
  }

  if ( getenv("TOCHNOG_DEBUG")!=NULL ) {
    str = getenv("TOCHNOG_DEBUG");
    if ( !strcmp("yes",str ) ) result = 1;
  }
  if ( result==0 ) return 0;

  if ( element>=0 ) {
    if ( getenv("TOCHNOG_ELEMENT")!=NULL ) {
      str = getenv("TOCHNOG_ELEMENT");
      if ( element!=atol(str) ) result=0;
    }
  }

  if ( db_active_index( ICONTROL, 0, VERSION_NORMAL ) ) {
    if ( getenv("TOCHNOG_ICONTROL")!=NULL ) {
      db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
      str = getenv("TOCHNOG_ICONTROL");
      if ( icontrol!=atol(str) ) result=0;
    }
  }

  if ( db_active_index( NUMBER_ITERATIONS, 0, VERSION_NORMAL ) ) {
    if ( getenv("TOCHNOG_ITERATION")!=NULL ) {
      db( NUMBER_ITERATIONS, 0, &iteration, ddum, ldum, VERSION_NORMAL, GET );
      str = getenv("TOCHNOG_ITERATION");
      if ( iteration!=atol(str) ) result=0;
    }
  }

  if ( inod>=0 ) {
    if ( getenv("TOCHNOG_NODE")!=NULL ) {
      str = getenv("TOCHNOG_NODE");
      if ( inod!=atol(str) ) result=0;
    }
  }

  if ( strcmp(routine,"") ) {
    if ( getenv("TOCHNOG_ROUTINE")!=NULL ) {
      str = getenv("TOCHNOG_ROUTINE");
      if ( strcmp(routine,str) ) result=0;
    }
  }

  return result;

}

long int stress_indx( long int idim, long int jdim )

{
  long int kdim=0, ldim=0, indx=0;

  if ( idim<jdim ) {
    kdim = idim;
    ldim = jdim;
  }
  else {
    kdim = jdim;
    ldim = idim;
  }

  if      ( kdim==0 ) indx = 0 + ldim;
  else if ( kdim==1 ) indx = 2 + ldim;
  else if ( kdim==2 ) indx = 3 + ldim;

  return indx;
}

long int sri_active( long int /*element*/, long int element_group,
  long int name, long int nnol )

  // group_element_selective_reduced_integration (SRI, Hughes): opt-in
  // fix for the shear locking of the bilinear elements in bending.
  // Returns the mode mask (SRI_SHEAR for -yes, SRI_VOLUMETRIC for
  // -volumetric, C80.1; SRI_LAMBDA for -lambda, C80.2) when the full SRI
  // applies to the element, 0 otherwise:
  //   - 2D bilinear quad4 (nnol==4, ndim==2)
  //   - 3D trilinear hex8 (nnol==8, ndim==3) — EXTENSION 2026-08-29
  //     (same mechanism: the shear strains gamma_xy/gamma_xz/gamma_yz
  //     integrated at 1 Gauss point at the centroid)
  //   - NOT axisymmetric
  //   - NOT large displacement (materi_displacement): the reduced shear
  //     point is built from the reference coordinates passed to materi()
  //   - LINEAR ELASTICITY only: the split D = D_norm + D_shear is exact
  //     only when the tangent is constant over the element (the reduced
  //     point has no material state of its own)
  // Used by both pol() (which switches the full rule to Gauss 2x2 /
  // 2x2x2 for the normal terms) and materi() (which splits D and adds
  // the shear terms with 1 Gauss point at the centroid). Kept in ONE
  // place so the quadrature and the stiffness split can never disagree.
  //
  // WARNING (measured 2026-08-29, documented in the developer manual):
  // the 2D quad4 SRI is hourglass-free (the 2D "section warping" mode
  // has a normal strain and is stabilized by the full normal rule),
  // but the 3D hex8 SRI retains ZERO-ENERGY modes: the isolated
  // element has 9 (6 rigid + 3 twist) and a mesh has the additional
  // "section warping" modes (u_y = A(x)*(2z-1), u_z = A(x)*(2y-1):
  // all normal strains zero, the shear strains vanish at the section
  // centroid, so the 1-point rule misses them). The loaded
  // configurations (cantilever) therefore assemble a SINGULAR matrix.
  // The patch tests (constant states) and the rigid modes remain exact.
  // This is the known limitation of the shear-only selective
  // integration of the 8-node brick (the literature moved to the
  // B-bar/assumed-strain formulations for this reason).

{
  long int sri=-NO, axisymmetric=-NO;
  long int ldum=0;
  double ddum[1];

  if ( !( ( name==-QUAD4 && nnol==4 && ndim==2 ) ||
          ( name==-HEX8 && nnol==8 && ndim==3 ) ) ) return 0;
  db( GROUP_ELEMENT_SELECTIVE_REDUCED_INTEGRATION, element_group, &sri,
    ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  if ( sri!=-NO && sri!=-YES && sri!=-VOLUMETRIC && sri!=-LAMBDA )
    db_error( GROUP_ELEMENT_SELECTIVE_REDUCED_INTEGRATION, element_group );
  if ( sri==-NO ) return 0;
  if ( !materi_stress ) return 0;
  db( GROUP_AXISYMMETRIC, element_group, &axisymmetric, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  if ( axisymmetric==-YES ) return 0;
  if ( materi_displacement ) return 0;
  if ( materi_plasti_kappa || materi_strain_plasti || materi_plasti_f ||
    materi_plasti_softvar_local || materi_plasti_softvar_nonlocal ||
    materi_damage || materi_maxwell_stress ) return 0;
  // C80.1 (2026-10-16): mode mask. -yes = the shear terms only
  // (SRI_SHEAR, unchanged); -volumetric = the volumetric part kappa*m*m^T
  // only (SRI_VOLUMETRIC, B-bar / selective volumetric integration).
  // C80.2 (2026-10-16): -lambda = the lambda*m*m^T part only
  // (SRI_LAMBDA, Hughes 1977 lambda/mu split; the former -both value
  // was removed and is rejected by the input reader). The volumetric
  // and lambda parts need a 3D or plane-strain kinematics: a membrane
  // (plane stress) group has a dependent eps_zz, an effective lambda
  // 2*lambda*mu/(lambda+2*mu) that stays bounded and hence no
  // volumetric locking, so the bit is dropped there (full Gauss; the
  // caller warns once).
  long int mask = 0, membrane=-NO;
  if ( sri==-YES ) mask |= SRI_SHEAR;
  if ( sri==-VOLUMETRIC || sri==-LAMBDA ) {
    db( GROUP_MATERI_MEMBRANE, element_group, &membrane, ddum, ldum,
      VERSION_NORMAL, GET_IF_EXISTS );
    if ( membrane!=-YES )
      mask |= ( sri==-LAMBDA ? SRI_LAMBDA : SRI_VOLUMETRIC );
  }
  return mask;
}

// strain_recovery_weight - K11 (2026-10-06): nodal recovery weight of the
// strain tensors (materi_strain_total / _elasti / _plasti): the stress
// recovery weight (Professional per-element least-squares fit) except for
// SRI elements, which keep h (their stress recovery mixes Lagrange and h).
// Used for the numerator (materi.cc) and the lumped denominator and
// inertia term (general.cc) of the incremental strain update.
double strain_recovery_weight( long int sri_on, long int name, long int nnol,
  long int inol, long int npoint, long int ipoint, double h_inol, double volume,
  long int element_group, const long int nodes[] )

{
  if ( sri_on ) return h_inol;
  return sri_stress_recovery_weight( 0, name, nnol, inol, npoint, ipoint,
    h_inol, volume, element_group, nodes );
}

// Geometry-dependent per-element recovery operator F(inol,ipoint),
// measured on the Professional 02-08-2026 (element_intpnt_dof ->
// element_dof, residual <= 1e-12 of the field scale on distorted and
// curved elements; K21 2026-10-12 for the 2D quad rules, C70 2026-10-13
// for the 3D and 3x3 rules):
//   kind 1  linear least squares of 1,x,y (2D) / 1,x,y,z (3D) in PHYSICAL
//           coordinates over the integration points, evaluated at the node:
//           F = [1,x_n] (A^T A)^-1 [1,x_g]^T. Rules: quad4 2x2, quad9 2x2 and
//           3x3, hex8 2x2x2, hex27 2x2x2 and 3x3x3 Gauss, tet10 5-point
//           Keast, prism6 6-point (natural-coordinate fits coincide with it
//           only on affine elements: parallelograms, parallelepipeds,
//           straight-edged tet10);
//   kind 2  quad4 / hex8 1 point: the point value (F = 1);
//   axis    axisymmetric quad4 / quad9 (2x2 and 3x3) with a node EXACTLY at
//           r=0 (r=1e-6 does not trigger it): F = h_n(g) / sum_g h_n(g);
//   kind 3  Lobatto rules with the points at the nodes (quad4 2x2, quad9
//           3x3, hex8 2x2x2, hex27 3x3x3, C79 2026-10-13): the value of the
//           point at the node (F = Kronecker delta). The Lobatto -minimal
//           quad9 2x2 / hex27 2x2x2 (corner points) are kind 1 with the
//           Lobatto point positions.
//   C67 (2026-10-16): the same three rules for the 1D bars (bar2, bar3,
//           bar4: kind 1 = linear fit 1,x) and the quad16 (4x4 nodes):
//           Gauss (every point count), -reduced and Lobatto -minimal are
//           kind 1, Lobatto -maximal kind 3 (the bar4/quad16 interior nodes
//           sit at the Lobatto positions iso +-1/sqrt(5), interpolation_
//           polynomial, so the points are at the nodes); axis rule for the
//           on-axis quad16 too. hex64 is rejected by both codes.
// The operator depends on the element geometry, so it is built from the
// start-of-step nodal coordinates (NODE, VERSION_NORMAL) and cached per
// thread for the current element (the element loop may run in parallel
// threads). The integration points and shape functions below repeat the
// layout of pol() (polynom.cc): tensor point/node index xi fastest, tet10
// and prism6 point order of their pol() branches.
struct RecoveryFitCache {
  long int valid, name, nnol, npoint, group, axis, lobatto;
  double riso;
  long int nodes[MNOL];
  double coord[MNOL*MDIM];
  double F[MNOL][MPOINT];
  long int ok;
};

// Nodes per direction of the tensor Lagrange elements recovered by the
// fit (pol() layout: bars in 1D, quads in 2D, hexes in 3D), 0 otherwise.
static long int recovery_tensor_npol( long int name )

{
  if ( ndim==1 ) {
    if ( name==-BAR2 ) return 2;
    if ( name==-BAR3 ) return 3;
    if ( name==-BAR4 ) return 4;
  }
  else if ( ndim==2 ) {
    if ( name==-QUAD4 ) return 2;
    if ( name==-QUAD9 ) return 3;
    if ( name==-QUAD16 ) return 4;
  }
  else if ( ndim==3 ) {
    if ( name==-HEX8 ) return 2;
    if ( name==-HEX27 ) return 3;
  }
  return 0;
}

static long int recovery_fit_kind( long int name, long int nnol,
  long int npoint )

{
  long int npol=0, nfull=1, nminimal=1, d=0;

  if ( npoint==1 && ( ( ndim==2 && name==-QUAD4 && nnol==4 ) ||
                      ( ndim==3 && name==-HEX8  && nnol==8 ) ) ) return 2;
  if ( ndim==3 && name==-TET10  && nnol==10 && npoint==5 ) return 1;
  if ( ndim==3 && name==-PRISM6 && nnol==6  && npoint==6 ) return 1;
  npol = recovery_tensor_npol( name );
  if ( !npol ) return 0;
  for ( d=0; d<ndim; d++ ) {
    nfull *= npol;
    nminimal *= npol-1;
  }
  if ( nnol!=nfull ) return 0;
  // tensor elements. Lobatto rule with the points at the nodes (C79,
  // F4 T6f 2026-10-13; bars and quad16 C67 2026-10-16): the Professional
  // 02-08-2026 takes the point value per element and the arithmetic mean
  // over the elements (measured exact on distorted meshes of unequal
  // element size, plane and axisymmetric off-axis). Before, the h weight
  // times the point volume made the node a volume-weighted average.
  if ( !stress_recovery_gauss ) {
    if ( stress_recovery_lobatto && npoint==nnol ) return 3;
    return 0;
  }
  // Gauss rules (npol or npol-1 points per direction), the -reduced rules
  // (C80/C67, stress_recovery_gauss set by pol()) and the Lobatto -minimal
  // rules (npol-1 Lobatto points per direction; Professional: the same
  // physical linear fit, C79/C67). The 1-point rules (bar2 default, quad4/
  // hex8 -minimal) are the point value (kind 2 above, bar2 in
  // sri_stress_recovery_weight).
  if ( npoint>1 && ( npoint==nfull || npoint==nminimal ) ) return 1;
  return 0;
}

// Shape functions hp[ipoint][inol] of the element at its integration points.
static void recovery_fit_shape( long int name, long int nnol, long int npoint,
  long int lobatto, double riso, double hp[MPOINT][MNOL] )

{
  long int ipoint=0, inol=0, k=0, l=0, d=0, npol=0, m=1, mpow=1, stride=1;
  double niso[MPOINT], wn[MPOINT], xg[MPOINT], wg[MPOINT], lag[MDIM][MPOINT],
    iso[MDIM];

  if ( name==-TET10 ) {
    for ( ipoint=0; ipoint<npoint; ipoint++ ) {
      double L1 = ( ipoint==0 ) ? 0.25 : ( ( ipoint==1 ) ? 0.5 : 1./6. );
      double L2 = ( ipoint==0 ) ? 0.25 : ( ( ipoint==2 ) ? 0.5 : 1./6. );
      double L3 = ( ipoint==0 ) ? 0.25 : ( ( ipoint==3 ) ? 0.5 : 1./6. );
      double L4 = 1. - L1 - L2 - L3;
      hp[ipoint][0] = (2.*L1-1.)*L1;
      hp[ipoint][1] = 4.*L1*L2;
      hp[ipoint][2] = (2.*L2-1.)*L2;
      hp[ipoint][3] = 4.*L3*L1;
      hp[ipoint][4] = 4.*L2*L3;
      hp[ipoint][5] = (2.*L3-1.)*L3;
      hp[ipoint][6] = 4.*L1*L4;
      hp[ipoint][7] = 4.*L2*L4;
      hp[ipoint][8] = 4.*L3*L4;
      hp[ipoint][9] = (2.*L4-1.)*L4;
    }
    return;
  }
  if ( name==-PRISM6 ) {
    static const double tri_L1[3] = { 2./3., 1./6., 1./6. };
    static const double tri_L2[3] = { 1./6., 2./3., 1./6. };
    for ( ipoint=0; ipoint<npoint; ipoint++ ) {
      double zz = ( ipoint/3==0 ) ? 0.5-0.5/sqrt(3.) : 0.5+0.5/sqrt(3.);
      double L1 = tri_L1[ipoint%3], L2 = tri_L2[ipoint%3], L3 = 1.-L1-L2;
      hp[ipoint][0] = (1.-zz)*L1;
      hp[ipoint][1] = (1.-zz)*L2;
      hp[ipoint][2] = (1.-zz)*L3;
      hp[ipoint][3] = zz*L1;
      hp[ipoint][4] = zz*L2;
      hp[ipoint][5] = zz*L3;
    }
    return;
  }

  // tensor Lagrange elements (bar2/3/4, quad4/9/16, hex8/27) with m^ndim
  // Gauss points (Lobatto points for the Lobatto rules, C79). The nodes
  // sit at the Lobatto positions of npol points (interpolation_polynomial:
  // -1,1 / -1,0,1 / -1,-1/sqrt5,1/sqrt5,1).
  npol = recovery_tensor_npol( name );
  integration_lobatto( npol, niso, wn );
  while ( mpow<npoint ) {
    m++;
    mpow = 1;
    for ( d=0; d<ndim; d++ ) mpow *= m;
  }
  if ( lobatto ) integration_lobatto( m, xg, wg );
  else integration_gauss( m, xg, wg );
  // C80: the quad4/hex8 -reduced points at iso +-1/f
  if ( riso>0. && m==2 ) {
    xg[0] = -riso;
    xg[1] =  riso;
  }
  for ( ipoint=0; ipoint<npoint; ipoint++ ) {
    stride = 1;
    for ( d=0; d<ndim; d++ ) {
      iso[d] = xg[(ipoint/stride)%m];
      stride *= m;
    }
    for ( d=0; d<ndim; d++ ) {
      for ( k=0; k<npol; k++ ) {
        lag[d][k] = 1.;
        for ( l=0; l<npol; l++ ) {
          if ( l==k ) continue;
          lag[d][k] *= ( iso[d]-niso[l] ) / ( niso[k]-niso[l] );
        }
      }
    }
    for ( inol=0; inol<nnol; inol++ ) {
      hp[ipoint][inol] = 1.;
      stride = 1;
      for ( d=0; d<ndim; d++ ) {
        hp[ipoint][inol] *= lag[d][(inol/stride)%npol];
        stride *= npol;
      }
    }
  }
}

static void recovery_fit_build( RecoveryFitCache &c, long int kind )

{
  long int nt=ndim+1, inol=0, ipoint=0, k=0, l=0, r=0, d=0, pmax=0;
  double hp[MPOINT][MNOL], xp[MPOINT][MDIM], xc[MDIM], scale=0.,
    aug[MDIM+1][2*(MDIM+1)], piv=0., fac=0., phn[MDIM+1], php[MDIM+1];

  c.ok = 0;
  if ( kind==2 ) {
    for ( inol=0; inol<c.nnol; inol++ ) c.F[inol][0] = 1.;
    c.ok = 1;
    return;
  }
  if ( kind==3 ) {
    // Lobatto points at the nodes: the point value (same tensor index)
    for ( inol=0; inol<c.nnol; inol++ )
      for ( ipoint=0; ipoint<c.npoint; ipoint++ )
        c.F[inol][ipoint] = ( inol==ipoint ) ? 1. : 0.;
    c.ok = 1;
    return;
  }

  // physical coordinates of the integration points
  recovery_fit_shape( c.name, c.nnol, c.npoint, c.lobatto, c.riso, hp );
  for ( ipoint=0; ipoint<c.npoint; ipoint++ ) {
    for ( d=0; d<ndim; d++ ) {
      xp[ipoint][d] = 0.;
      for ( inol=0; inol<c.nnol; inol++ )
        xp[ipoint][d] += hp[ipoint][inol] * c.coord[inol*ndim+d];
    }
  }

  if ( c.axis ) {
    for ( inol=0; inol<c.nnol; inol++ ) {
      double sum = 0.;
      for ( ipoint=0; ipoint<c.npoint; ipoint++ ) sum += hp[ipoint][inol];
      if ( scalar_dabs(sum)<TINY ) return;
      for ( ipoint=0; ipoint<c.npoint; ipoint++ )
        c.F[inol][ipoint] = hp[ipoint][inol] / sum;
    }
    c.ok = 1;
    return;
  }

  // linear least squares 1,x(,y,z) on centred and scaled coordinates (the
  // fit is affine invariant; this only keeps the normal matrix well scaled)
  for ( d=0; d<ndim; d++ ) {
    xc[d] = 0.;
    for ( ipoint=0; ipoint<c.npoint; ipoint++ )
      xc[d] += xp[ipoint][d] / c.npoint;
  }
  for ( ipoint=0; ipoint<c.npoint; ipoint++ )
    for ( d=0; d<ndim; d++ )
      scale = scalar_dmax( scale, scalar_dabs(xp[ipoint][d]-xc[d]) );
  if ( scale<TINY ) return;
  for ( k=0; k<nt; k++ )
    for ( l=0; l<2*nt; l++ ) aug[k][l] = ( l==nt+k ) ? 1. : 0.;
  for ( ipoint=0; ipoint<c.npoint; ipoint++ ) {
    php[0] = 1.;
    for ( d=0; d<ndim; d++ ) php[1+d] = ( xp[ipoint][d]-xc[d] ) / scale;
    for ( k=0; k<nt; k++ )
      for ( l=0; l<nt; l++ ) aug[k][l] += php[k]*php[l];
  }
  // invert the normal matrix by Gauss-Jordan with partial pivoting
  for ( k=0; k<nt; k++ ) {
    pmax = k;
    for ( l=k+1; l<nt; l++ )
      if ( scalar_dabs(aug[l][k])>scalar_dabs(aug[pmax][k]) ) pmax = l;
    if ( pmax!=k ) {
      for ( l=0; l<2*nt; l++ ) {
        double tmp = aug[k][l]; aug[k][l] = aug[pmax][l]; aug[pmax][l] = tmp;
      }
    }
    piv = aug[k][k];
    if ( scalar_dabs(piv)<1.e-12 ) return;
    for ( l=0; l<2*nt; l++ ) aug[k][l] /= piv;
    for ( l=0; l<nt; l++ ) {
      if ( l!=k ) {
        fac = aug[l][k];
        for ( r=0; r<2*nt; r++ ) aug[l][r] -= fac*aug[k][r];
      }
    }
  }
  for ( inol=0; inol<c.nnol; inol++ ) {
    phn[0] = 1.;
    for ( d=0; d<ndim; d++ )
      phn[1+d] = ( c.coord[inol*ndim+d]-xc[d] ) / scale;
    for ( ipoint=0; ipoint<c.npoint; ipoint++ ) {
      php[0] = 1.;
      for ( d=0; d<ndim; d++ ) php[1+d] = ( xp[ipoint][d]-xc[d] ) / scale;
      double F = 0.;
      for ( k=0; k<nt; k++ )
        for ( l=0; l<nt; l++ ) F += phn[k]*aug[k][nt+l]*php[l];
      c.F[inol][ipoint] = F;
    }
  }
  c.ok = 1;
}

// Returns 1 and F(inol,ipoint) of the geometry-dependent recovery of the
// current element, 0 when it does not apply (or the fit is singular).
static long int recovery_fit_factor( long int name, long int nnol,
  long int inol, long int npoint, long int ipoint, long int element_group,
  const long int nodes[], double &F )

{
  static thread_local RecoveryFitCache cache = { 0, 0, 0, 0, 0, 0, 0, 0., {0}, {0.}, {{0.}}, 0 };
  long int kind=0, axis=0, axisymmetric=-NO, same=0, k=0, d=0, ldum=0,
    idum[1], lobatto=stress_recovery_lobatto;
  double riso=stress_recovery_reduced_iso;
  double xn[MDIM], coord[MNOL*MDIM], ddum[1];

  // Applies to every stress dof, convected ones included (K21 WU1b,
  // 2026-10-12: the convection gate was removed; the old test-2014
  // cd_method target was met by chance by an ill-posed model, and the
  // convected nodal stress is a known class B divergence, BACKLOG C72).
  if ( nodes==NULL ) return 0;
  kind = recovery_fit_kind( name, nnol, npoint );
  if ( !kind ) return 0;

  // Fast path: same element (type, group and node numbers) as the cached
  // one. The coordinates and the axis flag are re-read on the (inol 0,
  // ipoint 0) call, which is the first call of every element visit
  // (general() and materi() loop ipoint, then inol, from 0), so a geometry
  // change between passes (updated mesh, refinement) is always caught.
  same = cache.valid && cache.name==name && cache.nnol==nnol &&
    cache.npoint==npoint && cache.group==element_group &&
    cache.lobatto==lobatto && cache.riso==riso;
  for ( k=0; same && k<nnol; k++ )
    if ( cache.nodes[k]!=nodes[k] ) same = 0;
  if ( same && !( inol==0 && ipoint==0 ) ) {
    if ( !cache.ok ) return 0;
    F = cache.F[inol][ipoint];
    return 1;
  }

  for ( k=0; k<nnol; k++ ) {
    db( NODE, nodes[k], idum, xn, ldum, VERSION_NORMAL, GET );
    for ( d=0; d<ndim; d++ ) coord[k*ndim+d] = xn[d];
  }
  if ( kind==1 && ndim==2 ) {
    db( GROUP_AXISYMMETRIC, element_group, &axisymmetric, ddum, ldum,
      VERSION_NORMAL, GET_IF_EXISTS );
    if ( axisymmetric==-YES )
      for ( k=0; k<nnol; k++ ) if ( coord[k*ndim+0]==0. ) axis = 1;
  }

  same = same && cache.axis==axis;
  for ( k=0; same && k<nnol*ndim; k++ )
    if ( cache.coord[k]!=coord[k] ) same = 0;
  if ( !same ) {
    cache.name = name; cache.nnol = nnol; cache.npoint = npoint;
    cache.group = element_group; cache.axis = axis; cache.lobatto = lobatto;
    cache.riso = riso;
    for ( k=0; k<nnol; k++ ) cache.nodes[k] = nodes[k];
    for ( k=0; k<nnol*ndim; k++ ) cache.coord[k] = coord[k];
    recovery_fit_build( cache, kind );
    cache.valid = 1;
  }
  if ( !cache.ok ) return 0;
  F = cache.F[inol][ipoint];
  return 1;
}

double sri_stress_recovery_weight( long int sri, long int name, long int nnol,
  long int inol, long int npoint, long int ipoint, double h_inol, double volume,
  long int element_group, const long int nodes[] )

  // Nodal stress recovery weight of the lumped sigma-dof update.
  //
  // The Professional (binary 2-8-2026) recovers the nodal stress from the
  // integration-point stresses with a per-element least-squares fit and
  // assembles a node by the ARITHMETIC MEAN of the adjacent elements'
  // per-element nodal values (its node_lhside for a stress dof equals the
  // number of adjacent elements). The per-element fit was identified from
  // the Professional's own .dbs dumps (element_intpnt_dof -> element_dof)
  // and is, per element type:
  //   quad4 2x2, quad9 2x2 / 3x3, hex8 2x2x2, hex27 2x2x2 / 3x3x3 Gauss,
  //   tet10 5-point Keast, prism6 6 points: linear least squares 1,x,y(,z)
  //                  in PHYSICAL coordinates (K21 2026-10-12 quad4/quad9 2x2,
  //                  C70 2026-10-13 the others; before: natural-coordinate
  //                  fits, equal only on affine elements; prism6 had h);
  //                  axisymmetric quad4/quad9 with a node at r=0:
  //                  sum h sig / sum h
  //   quad4 / hex8   1 pt  the point value (K21 / C70; before: volume
  //                  weighted h)
  //   tria3 1 pt / 3 pt nodal, tria6 6 pt nodal, tet4 1 pt, bar2 1 pt:
  //            the point value (F=1) / the value at the node (F=h)
  //            (F4 T6f, 2026-10-11; before, these were area-weighted)
  //   quad4 2x2 / hex8 2x2x2 -reduced (points at iso +-1/f) and quad9 3x3
  //                  / hex27 3x3x3 -reduced (Lobatto points): the physical
  //                  linear fit, also with the points at the nodes (C80,
  //                  2026-10-16)
  //   quad4 / quad9 / hex8 / hex27 Lobatto with the points at the nodes:
  //                  the point value; quad9 / hex27 Lobatto -minimal
  //                  (corner points): the physical linear fit (C79,
  //                  2026-10-13; before: h, i.e. volume weighted, and
  //                  no contribution at the quad9/hex27 non-corner nodes
  //                  of the -minimal rules)
  //   bar2 2 points, bar3, bar4, quad16: the same rules per rule family
  //                  (Gauss / -reduced / Lobatto -minimal: the physical
  //                  linear fit 1,x(,y); Lobatto -maximal: the point value;
  //                  C67, 2026-10-16; before: h, volume weighted)
  // Every other element/rule returns h_inol: the unmeasured rules (tet4
  // -maximal, prism15, ...; hex64 is rejected by both codes).
  //
  // The lumped sigma-dof equation (general.cc/materi.cc) multiplies the
  // weight by the integration-point volume, so the assembled average is
  // sum(volume_g*w_g*sigma_g)/sum(volume_g*w_g). The Professional's
  // per-element value is the plain interpolation sum(F_g*sigma_g) with
  // sum_g F_g = 1, and its assembly is unweighted. Returning F_g/volume_g
  // makes volume_g*weight_g = F_g, so each element contributes weight 1
  // and the assembled node is the arithmetic mean, exactly as the
  // Professional. Without the volume cancellation the quadrature weights
  // bias the average (measured on a controlled quad9 linear patch: the
  // exact sigxx = 1000 came out 1333 = 4/3).

{
  double F=0.;

  // SRI path (GNU-specific feature, unchanged 2026-09-13): the classic
  // selective reduced integration of the quad4/hex8 uses the Lagrange
  // extrapolation of the Gauss-point values for the NORMAL stress
  // components (the caller only invokes this for those; the shear keeps
  // h). The volume-weighted assembly is kept (no volume cancellation)
  // so the SRI results are byte-identical to the pre-alignment build.
  if ( sri ) {
    if ( !( ( nnol==4 && npoint==4 ) || ( nnol==8 && npoint==8 ) ) )
      return h_inol;
    double xi_n  = ( inol%2   == 0 ? -1. : 1. );
    double eta_n = ( (inol/2)%2 == 0 ? -1. : 1. );
    double xi_g  = ( ipoint%2 == 0 ? -1. : 1. ) / sqrt(3.);
    double eta_g = ( (ipoint/2)%2 == 0 ? -1. : 1. ) / sqrt(3.);
    double w_xi  = ( xi_n  + xi_g  ) / ( 2.*xi_g  );
    double w_eta = ( eta_n + eta_g ) / ( 2.*eta_g  );
    if ( nnol==4 ) return w_xi * w_eta;
    double zeta_n = ( inol/4 == 0 ? -1. : 1. );
    double zeta_g = ( ipoint/4 == 0 ? -1. : 1. ) / sqrt(3.);
    double w_zeta = ( zeta_n + zeta_g ) / ( 2.*zeta_g  );
    return w_xi * w_eta * w_zeta;
  }

  // Simplex elements and the 1-point bar2 (F4 T6f, 2026-10-11): the
  // Professional 02-08-2026 assembles their nodes by the arithmetic mean
  // of the adjacent elements' values too (measured exact on two-element
  // probes of different size: tria3 1-point and 3-point nodal rule,
  // tria6 6-point nodal rule, tet4 1-point, bar2 1-point; GNU was
  // area/volume weighted). Per element the nodal value is the single
  // point value (1-point rules, F=1) or the value of the point at the
  // node (nodal rules, F=h, the Kronecker delta), so the weight is F/volume
  // and each element contributes weight 1, exactly as for the table
  // elements below. Other simplex rules (tet4 -maximal, bar2 2 points,
  // ...) are unmeasured and keep h_inol.
  if ( ( name==-TRIA3 && ( npoint==1 || npoint==3 ) ) ||
       ( name==-TRIA6 && npoint==6 ) ||
       ( name==-TET4  && npoint==1 ) ||
       ( name==-BAR2  && npoint==1 ) ) {
    if ( scalar_dabs(volume)>TINY ) return ( npoint==1 ? 1. : h_inol )/volume;
    return h_inol;
  }

  // Geometry-dependent per-element fit (K21 + C70, see recovery_fit_factor):
  // the quad4/quad9/hex8/hex27 Gauss rules, the 1-point quad4/hex8, the
  // tet10 Keast and the prism6 rules, including the axisymmetric on-axis
  // rule, and the quad4/quad9/hex8/hex27 Lobatto rules (C79). Each element
  // contributes weight 1 (F/volume).
  if ( recovery_fit_factor( name, nnol, inol, npoint, ipoint, element_group,
         nodes, F ) ) {
    if ( scalar_dabs(volume)>TINY ) return F/volume;
  }
  return h_inol;
}

char *long_to_a( long int n, char s[] )

{

  int i=0, sign=0, m=0;
  char *ptr=NULL;

  m = (int) n;

  if ( ( sign = m ) < 0 ) 
    m = -m;

  i = 0;
  do {
    s[i++] = m % 10 + '0';
  } while ( ( m/= 10 ) > 0 );

  if ( sign<0 ) 
    s[i++] = '-';
  s[i] = '\0';
  string_reverse( s );

  ptr = s;
  return ptr;
  
}

// Modern C++ version of long_to_a using std::string
std::string long_to_string_modern(long int n) {
    return std::to_string(n);
}

void exit_tn( long int print_database_type )

{

  long int itarget=0, ntarget=0, data_item_name=0, 
    data_item_index=0, number=0, correct=0, length=0, ldum=0,
    dof_label[MUKNWN], *target_item=NULL, *int_data=NULL;
  double value=0., tolerance=0., ddum[1], *dbl_data=NULL, *target_value=NULL;

  if ( any_runtime ) {
    pri( "*** Runtime file used in calculation. ***" );
  }

  if ( print_database_type==-RESTART ) {
    print_restart( -1 );
  }
  else {
    assert( print_database_type==-YES );
    // print_database_calculation (manual Professional 6.972): switch of
    // the final .dbs output at the end of the run; -no skips it (the
    // large1 performance test keeps the huge 3D dump away). Default -yes.
    long int print_database_calculation = -YES;
    db( PRINT_DATABASE_CALCULATION, 0, &print_database_calculation, ddum,
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
    if ( print_database_calculation!=-NO )
      print_database( -1, VERSION_NORMAL, -YES );
  }

  // print_gid_calculation: switch of the final GiD (.flavia) output.
  long int print_gid_calculation = -YES;
  db( PRINT_GID_CALCULATION, 0, &print_gid_calculation, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  if ( print_gid_calculation!=-NO )
    print_gid( -YES );
  if ( db_partialname_any("control_print_dx") ) print_dx( -YES );

  db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db_max_index( TARGET_ITEM, ntarget, VERSION_NORMAL, GET );
  if ( ntarget>=0 ) {
    for ( itarget=0; itarget<=ntarget; itarget++ ) {
      if ( db_active_index( TARGET_ITEM, itarget, VERSION_NORMAL ) &&
           db_active_index( TARGET_VALUE, itarget, VERSION_NORMAL ) ) {
        target_item     = db_int( TARGET_ITEM, itarget, VERSION_NORMAL );
        data_item_name  = labs(target_item[0]);
        data_item_index = target_item[1];
        if ( !db_active_index( data_item_name, data_item_index, VERSION_NORMAL ) ) {
          ofstream out( "tn.log", ios::app );
          out << "\nError in calculation with data file " << data_file << ".";
          out.close();
          exit(TN_EXIT_STATUS);
        }
        else {
          length = db_len( data_item_name, data_item_index, VERSION_NORMAL );
          if ( target_item[2]<0 ) {
            array_member(dof_label,target_item[2],nuknwn,number);
            if ( number<0 &&
                 ( data_item_name==NODE_DOF_CALCUL ||
                   data_item_name==POST_POINT_DOF_CALCUL ||
                   data_item_name==POST_LINE_DOF_CALCUL ||
                   data_item_name==POST_QUADRILATERAL_DOF_CALCUL ) ) {
              // Item names of the post_calcul output slots: -norx_sig,
              // -nory_sig, ... (post_calcul -materi_stress -force) and
              // -to_pres/-st_pres/-dy_pres (post_calcul
              // -groundflow_pressure -total_pressure/-static_pressure/
              // -dynamic_pressure). The slot is the position of the
              // generated item name in post_calcul_names (calcul.cc
              // calculate()); not found -> number stays -1 ->
              // db_error below.
              long int icalc=0;
              const char *item_name = db_name(labs(target_item[2]));
              for ( icalc=0; icalc<DATA_ITEM_SIZE; icalc++ ) {
                if ( post_calcul_names[icalc][0] &&
                     !strcmp( post_calcul_names[icalc], item_name ) ) {
                  number = icalc;
                  break;
                }
              }
            }
            if ( number>=0 &&
                 db_len(data_item_name,data_item_index,VERSION_NORMAL)==npuknwn ) 
              number /= nder;
          }
          else
            number = target_item[2];
          if ( number<0 || number>length-1 ) db_error( TARGET_ITEM, itarget );
          target_value = db_dbl( TARGET_VALUE, itarget, VERSION_NORMAL );
          value        = target_value[0];
          tolerance    = target_value[1];
          correct      = 1;
          if ( db_type(data_item_name)==INTEGER ) {
            int_data = db_int( data_item_name, data_item_index, VERSION_NORMAL );
            if ( int_data[number]<(value-tolerance) ||
                 int_data[number]>(value+tolerance) ) correct = 0;
          }
          else {
            dbl_data = db_dbl( data_item_name, data_item_index, VERSION_NORMAL );
            if ( dbl_data[number]<(value-tolerance) ||
                 dbl_data[number]>(value+tolerance) ) correct = 0;
          }
          if ( !correct ) {
            if ( check_target==-YES ) {
              ofstream out( "tn.log", ios::app );
              out << "\nError in calculation with data file " << data_file << ".";
              out << "\nTarget value for: ";
              out << "data item " << db_name(data_item_name) << " ";
              out << "with index " << data_item_index << " ";
              if ( target_item[2]<0 ) {
                out << "for " << db_name(target_item[2]) << " ";
              }
              else
                out << "and value number " << number << " ";
              out << "is " << value << ".";
              if ( db_type(data_item_name)==INTEGER )
                out << "\nThe actual value is " << int_data[number] << ".";
              else
                out << "\nThe actual value is " << dbl_data[number] << ".";
              out << "\n";
              out.close();
              exit(TN_EXIT_STATUS);
            }
            else {
              ofstream out( "tn.log", ios::app );
              out << "\nNote in calculation with data file " << data_file << ".";
              out << "\nTarget value (check_target -no) not met for: ";
              out << db_name(data_item_name) << " index " << data_item_index;
              out << " (wanted " << value << ").";
              out << "\n";
              out.close();
            }
          }
        }
      }
    }
  }

  if ( check_used==-YES ) check_used_report();

  // check_memory / check_memory_usage: report peak memory usage (GB)
  if ( check_memory==-YES || check_memory_usage==-YES ) {
    struct rusage ru;
    long int idum_m[1];
    if ( getrusage( RUSAGE_SELF, &ru )==0 ) {
      check_memory_usage_result = (double) ru.ru_maxrss / 1048576.;
      cout << "Peak memory usage: " << check_memory_usage_result << " GB.\n";
      cout << flush;
      long int mem_len=1;
      db( CHECK_MEMORY_USAGE_RESULT, 0, idum_m, &check_memory_usage_result,
        mem_len, VERSION_NORMAL, PUT );
    }
  }

  // zip (manual Professional 6.1095): zip all *flavia*, *msh, vtk,
  // *.plt and *.dbs files with the gzip program at the end of the
  // calculation (the gzip program must be installed). Read BEFORE
  // db_close (the record lives in the database).
  {
    long int zip = -NO;
    db( ZIP, 0, &zip, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    if ( zip==-YES ) {
      // only existing files are gzipped (a glob without matches must
      // not count as a failure)
      int rc = system( "for f in *flavia* *msh vtk* *.plt *.dbs ; do "
                       "[ -f \"$f\" ] && gzip -f \"$f\" ; done "
                       ">/dev/null 2>&1" );
      if ( rc!=0 )
        pri( "Warning: zip -yes but gzipping the output files failed "
             "(is the gzip program installed?)" );
      else
        pri( "zip -yes: output files zipped with gzip" );
    }
  }

  db_close();

  ofstream out( "tn.log", ios::app );
  out << "\nCalculation with data file " << data_file << " ready.\n";
  out.close();
  exit(0);

}


void exit_tn_on_error( void )

{
  long int i=0, element_group=0, iarea=0, ldum=0;
  double ddum[1];

  if ( strcmp(swit_routine_stack[0],"") ) {
    pri( "\n\nLast routines:" );
    for ( i=0; i<MSTACK; i++ ) {
      if ( strcmp( swit_routine_stack[i], "" ) ) pri( swit_routine_stack[i] );
    }
  }
  if ( swit_node_stack>=0 ) pri( "\n\nLast node:" , swit_node_stack );
  if ( swit_element_stack>=0 ) {
    pri( "\n\nLast element:" , swit_element_stack );
    element_group = 0;
    db( ELEMENT_GROUP, swit_element_stack, &element_group, 
      ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    pri( "Last element has ELEMENT_GROUP", element_group );
    if ( db_active_index( ELEMENT_GROUP_AREA_ELEMENT_GROUP, swit_element_stack, VERSION_NORMAL ) ) {
      db( ELEMENT_GROUP_AREA_ELEMENT_GROUP, swit_element_stack, &iarea, 
        ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
      pri( "Last element got ELEMENT_GROUP from AREA_ELEMENT_GROUP with index", iarea );
    }
    if ( db_active_index( ELEMENT_GROUP_AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP, 
        swit_element_stack, VERSION_NORMAL ) ) {
      db( ELEMENT_GROUP_AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP, swit_element_stack, &iarea, 
        ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
      pri( "Last element got ELEMENT_GROUP from AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP with index", iarea );
    }
  }

  parallel_sys_lock();
  parallel_active = 0;
  print_database( -1, VERSION_NORMAL, -EVERYTHING );
  print_gid( -YES );
  cout << flush;
  exit(TN_EXIT_STATUS);
  parallel_sys_unlock();
}
