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
  // Returns 1 when the keyword is -yes AND the full SRI applies to the
  // element:
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
  if ( sri!=-YES ) return 0;
  if ( !materi_stress ) return 0;
  db( GROUP_AXISYMMETRIC, element_group, &axisymmetric, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  if ( axisymmetric==-YES ) return 0;
  if ( materi_displacement ) return 0;
  if ( materi_plasti_kappa || materi_strain_plasti || materi_plasti_f ||
    materi_plasti_softvar_local || materi_plasti_softvar_nonlocal ||
    materi_damage || materi_maxwell_stress ) return 0;
  return 1;
}

// strain_recovery_weight - K11 (2026-10-06): nodal recovery weight of the
// strain tensors (materi_strain_total / _elasti / _plasti): the stress
// recovery weight (Professional per-element least-squares fit) except for
// SRI elements, which keep h (their stress recovery mixes Lagrange and h).
// Used for the numerator (materi.cc) and the lumped denominator and
// inertia term (general.cc) of the incremental strain update.
double strain_recovery_weight( long int sri_on, long int name, long int nnol,
  long int inol, long int npoint, long int ipoint, double h_inol, double volume )

{
  if ( sri_on ) return h_inol;
  return sri_stress_recovery_weight( 0, name, nnol, inol, npoint, ipoint,
    h_inol, volume );
}

double sri_stress_recovery_weight( long int sri, long int name, long int nnol,
  long int inol, long int npoint, long int ipoint, double h_inol, double volume )

  // Nodal stress recovery weight of the lumped sigma-dof update.
  //
  // The Professional (binary 2-8-2026) recovers the nodal stress from the
  // integration-point stresses with a per-element least-squares fit and
  // assembles a node by the ARITHMETIC MEAN of the adjacent elements'
  // per-element nodal values (its node_lhside for a stress dof equals the
  // number of adjacent elements). The per-element fit was identified from
  // the Professional's own .dbs dumps (element_intpnt_dof -> element_dof)
  // and is, per element type:
  //   quad4   4 pts  bilinear interpolation  1,x,y,xy  (4 terms)
  //   quad9   9 pts  linear least squares    1,x,y     (3 terms)
  //   hex8    8 pts  linear least squares    1,x,y,z   (4 terms)
  //   hex27  27 pts  linear least squares    1,x,y,z   (4 terms)
  // Every other element/rule returns h_inol. The node-containing Lobatto
  // rules already have h as the exact weight (h is the Kronecker delta at
  // the nodes), and the non-tensor elements (bar/tria/tet) keep the
  // historic h weighting.
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
  long int nterm=0, npol=0, nxi=1, neta=1, nzeta=1, ixi=0, ieta=0, izeta=0,
    inol_xi=0, inol_eta=0, inol_zeta=0, k=0, l=0, d=0, p=0, m=0, pmax=0;
  long int expo[4][MDIM];
  double xg[MPOINT], wg[MPOINT], node[MDIM], ipc[MDIM], pc[MDIM],
    phn[4], phi[4], M[4][4], aug[4][8], piv=0., fac=0., F=0.;

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

  if ( !stress_recovery_gauss ) return h_inol;

  for ( k=0; k<4; k++ )
    for ( d=0; d<MDIM; d++ ) expo[k][d] = 0;

  if ( name==-QUAD4 && nnol==4 && npoint==4 && ndim==2 ) {
    npol = 2; nterm = 4;               // bilinear 1,x,y,xy
    expo[1][0]=1; expo[2][1]=1; expo[3][0]=1; expo[3][1]=1;
  }
  else if ( name==-QUAD9 && nnol==9 && npoint==9 && ndim==2 ) {
    npol = 3; nterm = 3;               // linear 1,x,y
    expo[1][0]=1; expo[2][1]=1;
  }
  else if ( name==-HEX8 && nnol==8 && npoint==8 && ndim==3 ) {
    npol = 2; nterm = 4;               // linear 1,x,y,z
    expo[1][0]=1; expo[2][1]=1; expo[3][2]=1;
  }
  else if ( name==-HEX27 && nnol==27 && npoint==27 && ndim==3 ) {
    npol = 3; nterm = 4;               // linear 1,x,y,z
    expo[1][0]=1; expo[2][1]=1; expo[3][2]=1;
  }
  else return h_inol;

  if ( ndim>=1 ) nxi = npol;
  if ( ndim>=2 ) neta = npol;
  if ( ndim==3 ) nzeta = npol;

  inol_xi = inol % nxi; inol_eta = (inol/nxi) % neta;
  inol_zeta = inol/(nxi*neta);
  node[0] = ( nxi>1 ) ? -1. + 2.*inol_xi/(nxi-1) : 0.;
  if ( ndim>=2 ) node[1] = ( neta>1 ) ? -1. + 2.*inol_eta/(neta-1) : 0.;
  if ( ndim==3 ) node[2] = ( nzeta>1 ) ? -1. + 2.*inol_zeta/(nzeta-1) : 0.;

  ixi = ipoint % nxi; ieta = (ipoint/nxi) % neta; izeta = ipoint/(nxi*neta);
  integration_gauss( npol, xg, wg );
  ipc[0] = xg[ixi];
  if ( ndim>=2 ) ipc[1] = xg[ieta];
  if ( ndim==3 ) ipc[2] = xg[izeta];

  // normal matrix M = sum_p phi(p) phi(p)^T (UNIFORM weights, as measured)
  for ( k=0; k<nterm; k++ )
    for ( l=0; l<nterm; l++ ) M[k][l] = 0.;
  for ( p=0; p<npoint; p++ ) {
    long int pxi = p % nxi, peta = (p/nxi) % neta, pzeta = p/(nxi*neta);
    pc[0] = xg[pxi];
    if ( ndim>=2 ) pc[1] = xg[peta];
    if ( ndim==3 ) pc[2] = xg[pzeta];
    for ( k=0; k<nterm; k++ ) {
      phi[k] = 1.;
      for ( d=0; d<ndim; d++ )
        for ( long int e=0; e<expo[k][d]; e++ ) phi[k] *= pc[d];
    }
    for ( k=0; k<nterm; k++ )
      for ( l=0; l<nterm; l++ ) M[k][l] += phi[k]*phi[l];
  }

  // invert M (nterm <= 4) by Gauss-Jordan with partial pivoting
  for ( k=0; k<nterm; k++ ) {
    for ( l=0; l<nterm; l++ ) aug[k][l] = M[k][l];
    for ( l=0; l<nterm; l++ ) aug[k][nterm+l] = ( k==l ) ? 1. : 0.;
  }
  for ( k=0; k<nterm; k++ ) {
    pmax = k;
    for ( l=k+1; l<nterm; l++ )
      if ( scalar_dabs(aug[l][k])>scalar_dabs(aug[pmax][k]) ) pmax = l;
    if ( pmax!=k ) {
      for ( l=0; l<2*nterm; l++ ) {
        double tmp = aug[k][l]; aug[k][l] = aug[pmax][l]; aug[pmax][l] = tmp;
      }
    }
    piv = aug[k][k];
    if ( scalar_dabs(piv)<TINY ) return h_inol;
    for ( l=0; l<2*nterm; l++ ) aug[k][l] /= piv;
    for ( l=0; l<nterm; l++ ) {
      if ( l!=k ) {
        fac = aug[l][k];
        for ( m=0; m<2*nterm; m++ ) aug[l][m] -= fac*aug[k][m];
      }
    }
  }

  // F = phi(node)^T M^{-1} phi(ipoint)
  for ( k=0; k<nterm; k++ ) {
    phn[k] = 1.;
    for ( d=0; d<ndim; d++ )
      for ( long int e=0; e<expo[k][d]; e++ ) phn[k] *= node[d];
  }
  for ( k=0; k<nterm; k++ ) {
    phi[k] = 1.;
    for ( d=0; d<ndim; d++ )
      for ( long int e=0; e<expo[k][d]; e++ ) phi[k] *= ipc[d];
  }
  F = 0.;
  for ( k=0; k<nterm; k++ )
    for ( l=0; l<nterm; l++ ) F += phn[k]*aug[k][nterm+l]*phi[l];

  if ( scalar_dabs(volume)>TINY ) return F/volume;
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
