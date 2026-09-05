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

long int sri_active( long int element, long int element_group,
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

double sri_stress_recovery_weight( long int nnol, long int inol,
  long int npoint, long int ipoint, double h_inol, long int sri_active )

  // Nodal stress recovery weight of the lumped sigma-dof update
  // (DIAG-SOLVE-MIXTO lot C/D, fix D-b). The sigma dofs are advanced by
  // the lumped "inertia" equation (general.cc) with the shape function
  // h as the weight: sigma_node = sum_gp h*sigma_gp / sum_gp h. With
  // the NODE-CONTAINING quadratures (the default 2x2 Lobatto corners of
  // the quad4, the 2x2x2 Lobatto corners of the hex8, the quad9/hex27
  // Lobatto rules) h is the Kronecker delta and the recovery is exact.
  // With the interior GAUSS rules (2x2 for the SRI quad4, 2x2x2 for the
  // SRI hex8) the h-weighted average DILUTES the nodal values at the
  // corners (measured: 0.577x of the exact value for the bilinear), so
  // the section moments of the materi_stress_force integration read
  // systematically low values even when the displacement field is
  // correct (SRI quad4: mom = 0.29x instead of the element's 93.75%).
  //
  // The consistent recovery is the evaluation of the element's stress
  // field at the nodes - the bilinear/trilinear Lagrange extrapolation
  // of the Gauss-point values (the "same B at the node"). For the
  // 2-point 1D rules the extrapolation weight of the node iso
  // coordinate xi_n in {+1,-1} for the Gauss point xi_g in
  // {+1/sqrt(3), -1/sqrt(3)} is
  //   w(xi_n, xi_g) = prod_{g' != g} (xi_n - xi_g')/(xi_g - xi_g')
  // which for the corner Lobatto rule (xi_g = +-1) reduces to the
  // Kronecker delta, i.e. w = h (the node IS the integration point).
  // Returns the 2D tensor product w_xi * w_eta for the bilinear quad4
  // with the 2x2 rule and the 3D tensor product w_xi*w_eta*w_zeta for
  // the trilinear hex8 with the 2x2x2 rule (row-major-from-bottom
  // ordering of nodes and points, the same convention as pol() and the
  // SRI centroid B), and h_inol for every other case (quad9/hex27/
  // 1-point rules: unchanged).

{
  if ( !( sri_active && ( ( nnol==4 && npoint==4 ) ||
                          ( nnol==8 && npoint==8 ) ) ) )
    return h_inol;

  // node iso coordinates (+-1, row-major from bottom): quad4 inol 0 =
  // (-1,-1), 1 = (+1,-1), 2 = (-1,+1), 3 = (+1,+1); hex8 inol 0 =
  // (-1,-1,-1), 1 = (+1,-1,-1), ..., 7 = (+1,+1,+1). Point iso
  // coordinates (+-1/sqrt(3)): ipoint = izeta*nxi*neta + ieta*nxi +
  // ixi (the pol() izeta->ieta->ixi loop with nxi = neta = nzeta = 2).
  double xi_n  = ( inol%2   == 0 ? -1. : 1. );
  double eta_n = ( (inol/2)%2 == 0 ? -1. : 1. );
  double xi_g  = ( ipoint%2 == 0 ? -1. : 1. ) / sqrt(3.);
  double eta_g = ( (ipoint/2)%2 == 0 ? -1. : 1. ) / sqrt(3.);
  double w_xi  = ( xi_n  + xi_g  ) / ( 2.*xi_g  );
  double w_eta = ( eta_n + eta_g ) / ( 2.*eta_g );
  if ( nnol==4 ) return w_xi * w_eta;

  double zeta_n = ( inol/4 == 0 ? -1. : 1. );
  double zeta_g = ( ipoint/4 == 0 ? -1. : 1. ) / sqrt(3.);
  double w_zeta = ( zeta_n + zeta_g ) / ( 2.*zeta_g );
  return w_xi * w_eta * w_zeta;
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
