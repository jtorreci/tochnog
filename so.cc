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

/* Modified on September 28, 2010 by Fernando Lorenzo to add the LAPACK solver
   capability to tochnog
*/

#include "tochnog.h"

  // c routines
extern "C" 
  int dsbev_(char *jobz, char *uplo, integer *n, integer *kd, 
    doublereal *ab, integer *ldab, doublereal *w, doublereal *z, integer *
    ldz, doublereal *work, integer *info);
extern "C" 
  int dsbgv_(char *jobz, char *uplo, integer *n, integer *ka, 
    integer *kb, doublereal *ab, integer *ldab, doublereal *bb, integer *
    ldbb, doublereal *w, doublereal *z, integer *ldz, doublereal *work, 
    integer *info);
// Lapack band Solver
extern "C"
  int dgbsv_(integer *n, integer *kl, integer *ku, integer *
		   nrhs, double *mat, integer *ldmat, integer *ipiv, double *solve_b, 
		   integer *ldb, integer *info);
//Lapack Band Solver
extern "C" 
  long int solve_iterative_petsc( double **solve_A, double solve_b[],
    long int solve_nlocal, int nnz[], int **inz, 
    char *ksptype, char *pctype, long int mg, long int &iter );

extern "C" 
  long int solve_superlu( double *superlu_A, int *superlu_asub, int *superlu_xa, 
    int superlu_nnz, double solve_b[], long int solve_nlocal, long int nthread );

double *solve_solution;

#define EPS_MAT 1.e-8
#define EPS_EIGEN 1.e-10
#define EPS_SOL 1.e15
#define INITIAL_SIZE 100

void solve( long int task )

{
  integer n=0, kl=0, ku=0, ldmat=0, ldmatlin=0, 
    ldmatss=0, lmat=0, lmatlin=0, lmatss=0, lz=0, 
    lw=0, lwork=0, info=0, kd=0, ldz=0, nrhs=0;
  long int i=0, a=0, b=0, imat=0, ii=0, jj=0, index=0, inod=0, jnod=0, nnod=0,
    max_node=0, ipuknwn=0, iuknwn=0, element=0, max_element=0,
    element_group=0, swit=0, length=0, band=0, itmp=0,
    iprinc=0, ilocal=0, jlocal=0, iglobal=0, jglobal=0, indx=0,
    lglobal_local=0, ready=0, iter=0, anything_ordered=0, 
    icontrol=0, ieigen=0, neigen=0, neigen_positive=0,
    new_node=0, ldum=0, succesful=0, test1=0, test2=0,
    mg=0, nthread=1, band_solver=0, bicg_solver=0,lapack_solver=0,
    petsc_solver=0, superlu_solver=0, use_node_lhside=0, 
    control_options_solver_petsc_ksptype=-1, control_options_solver_petsc_pctype=-1,
    idum[1], control_eigen[2], *node_bounded=NULL, *dof_principal=NULL, 
    *node_node=NULL, *dof_label=NULL, *dof_type=NULL,
    *element_matrix_unknowns=NULL, *ipiv=NULL, 
    *ordered_global=NULL, *global_ordered=NULL, *ordering_has_been_done=NULL;
  int superlu_nnz=0, *nnz=NULL, **inz=NULL,  *length_inz=NULL, *int_ptr=NULL,
    *superlu_asub=NULL, *superlu_xa=NULL;
  double max=0., tmp=0., control_eigen_scale=0., ddum[1],
    *node_eigen=NULL, *options_relaxation=NULL,
    *mat=NULL, *matlin=NULL, *matss=NULL, *w=NULL, *z=NULL, *work=NULL, 
    *node_lhside=NULL, *node_rhside=NULL, *element_matrix_values=NULL,  
    *element_matrix_second_values=NULL, *node_dof_new=NULL,
    **solve_A=NULL, *superlu_A=NULL, *dbl_ptr=NULL, *solve_b_temp=NULL;
  char jobz[10], uplo[10], *ksptype=NULL, *pctype=NULL, cdum[MCHAR];

  db_highest_index( ELEMENT, max_element, VERSION_NORMAL );
  if ( max_element<0 ) return;

  db_highest_index( NODE, max_node, VERSION_NORMAL );
  if ( max_node<0 ) return;

  swit = set_swit(-1,-1,"solve");
  if ( swit ) pri( "In routine SOLVE" );
  if ( swit ) pri( "task", task );

  node_bounded = get_new_int(MPUKNWN);
  dof_principal = get_new_int(MUKNWN);
  length = db_data_length( NODE_NODE );
  node_node = get_new_int(length);     
  dof_label = get_new_int(MUKNWN);
  dof_type = get_new_int(MUKNWN);
  node_eigen = get_new_dbl(DATA_ITEM_SIZE);
  options_relaxation = get_new_dbl(MPRINC);

  db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( DOF_TYPE, 0, dof_type, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( DOF_PRINCIPAL, 0, dof_principal, ddum, ldum, VERSION_NORMAL, GET );
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  db( OPTIONS_PROCESSORS, 0, &nthread, ddum, ldum, VERSION_NORMAL, GET );
  array_set( options_relaxation, 1., nprinc );
  db( OPTIONS_RELAXATION, 0, idum, options_relaxation, ldum, 
    VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_OPTIONS_RELAXATION, icontrol, idum, options_relaxation, ldum, 
    VERSION_NORMAL, GET_IF_EXISTS );

  if ( SUPERLU_USE && PETSC_USE ) {
    pri( "Error: do not set both of SUPERLU_USE and PETSC_USE to 1" );
    exit(1);
  }

  if ( task==-MATRIX_ITERATIVE_BICG ) {
    bicg_solver = 1;
  }
  
  if ( task==-MATRIX_ITERATIVE_PETSC ) {
    if ( !PETSC_USE ) {
      pri( "Error: set PETSC_USE to 1 in tnpetsc.h" );
      exit(TN_EXIT_STATUS);
    }
    petsc_solver = 1;
    strcpy(cdum,"none");
    ksptype = cdum;
    pctype = cdum;
  }

  if ( task==-MATRIX_SUPERLU ) {
    if ( !SUPERLU_USE && !SUPERLU_MT_USE && !SUPERLU_DIST_USE ) {
      pri( "Error: set SUPERLU_USE to 1 in tnsuplu.h and recompile." );
      exit(TN_EXIT_STATUS);
    }
    superlu_solver = 1;
  }
  if ( task==-MATRIX_SUPERLU_DIST ) {
    if ( !SUPERLU_DIST_USE ) {
      pri( "Error: set SUPERLU_DIST_USE to 1 in tnsuplu.h and recompile." );
      exit(TN_EXIT_STATUS);
    }
    superlu_solver = 1;
  }
  if ( task==-MATRIX_SUPERLU_MT ) 
  {
    if ( !SUPERLU_MT_USE ) 
	{
      pri( "Error: set SUPERLU_MT_USE to 1 in tnsuplu.h and recompile." );
      exit(TN_EXIT_STATUS);
    }
    superlu_solver = 1;
  }

	/*----- lapack  ----*/  
	if ( task==-MATRIX_LAPACK ) 
	{   cout<< "Solving using Lapack Band Solver ..."<<endl;
		lapack_solver=1;
		band_solver = 1;
	}
	/*----- lapack  ----*/  

  array_set( control_eigen, 0, 2 );
  if ( task==-CONTROL_EIGEN )
  {
    band_solver = 1;
    db( CONTROL_EIGEN, icontrol, control_eigen, ddum, length, VERSION_NORMAL, GET );
    if ( swit ) pri( "control_eigen", control_eigen, 2 );
    if ( !db( CONTROL_EIGEN_SCALE, icontrol, idum, &control_eigen_scale,
      ldum, VERSION_NORMAL, GET_IF_EXISTS ) ) control_eigen_scale = 1;
    db_highest_index( NODE_LHSIDE, max_node, VERSION_NORMAL );
    if ( max_node<0 )
    {
      cout << "\nError: CONTROL_EIGEN is only possible after at least ";
      cout << "once CONTROL_TIMESTEP\n";
    }
  }

  solve_options_matrix_group = -NO;
  db( OPTIONS_MATRIX_GROUP, 0, &solve_options_matrix_group, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );

  length = 2*MNOL*nprinc*MNOL*nprinc;
  element_matrix_unknowns = get_new_int( length );

  length = (1+max_node)*nprinc;
  solve_b = get_new_dbl( length );
  solve_b_temp = get_new_dbl( length );
  array_set( solve_b, 0., length );
  array_set( solve_b_temp, 0., length );

  length = (1+max_node)*nprinc;
  solve_x = get_new_dbl( length );
  array_set( solve_x, 0., length );

  lglobal_local = (1+max_node)*npuknwn*2;
  solve_global_local = get_new_int( lglobal_local );
  array_set( solve_global_local, -NO, lglobal_local );

  length = 1+max_node;
  global_ordered = get_new_int( length );
  ordered_global = get_new_int( length );
  ordering_has_been_done = get_new_int( length );
  array_set( global_ordered, -NO, length );
  array_set( ordered_global, -NO, length );
  array_set( ordering_has_been_done, 0, length );

    // order nodes numbers to minimize band width
  if ( bicg_solver )
  {
    for ( inod=0; inod<=max_node; inod++ )
    {
      if ( db_active_index( NODE, inod, VERSION_NORMAL ) )
      {
        global_ordered[inod] = inod;
        ordered_global[inod] = inod;
      }
    }
  }
  else
  {
    new_node = 0;
    global_ordered[max_node] = new_node;
    ordered_global[new_node] = max_node;
    while ( !ready )
    {
      anything_ordered = 0;
      for ( inod=0; inod<=max_node; inod++ )
      {
        if ( global_ordered[inod]!=-NO && !ordering_has_been_done[inod] )
        {
          anything_ordered = 1;
          ordering_has_been_done[inod] = 1;
          db( NODE_NODE, inod, node_node, ddum, nnod, VERSION_NORMAL, GET );
          for ( jj=0; jj<nnod; jj++ )
          {
            jnod = labs(node_node[jj]);
            if ( global_ordered[jnod]==-NO )
            {
              new_node++;
              global_ordered[jnod] = new_node;
              ordered_global[new_node] = jnod;
              itmp = ( 1 + labs( global_ordered[inod] -
                global_ordered[jnod] ) ) * nprinc;
              if ( itmp>band ) band = itmp;
            }
          }
        }
      }
      if ( anything_ordered )
        ready = 0;
      else
      {
        ready = 1;
          // if needed, take new node if front has ended
        for ( inod=0; ready && inod<=max_node; inod++ )
        {
          if ( db_active_index( NODE, inod, VERSION_NORMAL ) &&
	       global_ordered[inod]==-NO )
          {
            new_node++;
	        global_ordered[inod] = new_node;
            ordered_global[new_node] = inod;
	        ready = 0;
          }
        }
      }
    }
  }
  if ( swit )
  {
    pri( "global_ordered", global_ordered, 1+max_node );
    pri( "ordered_global", ordered_global, 1+max_node );
  }

    // fill global = TOCHNOG -> local = SOLVER = TOCHNOG - BOUNDED info
  solve_nlocal = 0;
  for ( ii=0; ii<=max_node; ii++ )
  {
    inod = ordered_global[ii];
    if ( inod>=0 ) {
      array_set( node_bounded, 0, npuknwn );
      db( NODE_BOUNDED, inod, node_bounded, ddum, ldum, 
        VERSION_NORMAL, GET_IF_EXISTS );
      node_rhside = db_dbl( NODE_RHSIDE, inod, VERSION_NORMAL );
      for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ )
      {
        iuknwn = ipuknwn*nder;
        iglobal = inod*npuknwn + ipuknwn;
        test1 = !node_bounded[ipuknwn];
        test2 = dof_principal[iuknwn]>=0;
        if ( test1 && test2 )
        {
          solve_global_local[iglobal] = ilocal;
          solve_b[ilocal] = node_rhside[ipuknwn];
		  solve_b_temp[ilocal+1] = node_rhside[ipuknwn];
          ilocal++;
          solve_nlocal++;
        }
      }
    }
  }
  if ( band>solve_nlocal ) band = solve_nlocal;
  if ( swit )
  {
    pri( "solve_nlocal", solve_nlocal );
    pri( "band", band );
    pri( "right hand side", solve_b, solve_nlocal );
    pri( "solve_global_local", solve_global_local, lglobal_local );
  }

  if ( solve_nlocal!=0 )
  {
    if ( task==-CONTROL_EIGEN )
    {
      neigen = control_eigen[1];
      if ( swit ) pri( "neigen", neigen );
      if ( neigen>solve_nlocal ) neigen = solve_nlocal;
      if ( neigen>0 ) db_set_dbl( NODE_EIGEN, VERSION_NORMAL );
      if ( neigen*nuknwn>DATA_ITEM_SIZE )
      {
        pri( "Error: DATA_ITEM_SIZE is too small for NODE_EIGEN." );
        pri( "It should become at least", neigen*nuknwn );
        pri( "Change it in tochnog.h and recompile." );
        exit(TN_EXIT_STATUS);
      }
    }
    if ( petsc_solver )
    {
      solve_solution = get_new_dbl( solve_nlocal );
      array_set( solve_solution, 0., solve_nlocal );
    }
    if ( petsc_solver || superlu_solver )
    {
      nnz = get_new_int_short( solve_nlocal );
      length_inz = get_new_int_short( solve_nlocal );
      solve_A = new double * [solve_nlocal];
      inz = new int * [solve_nlocal];
      for ( ilocal=0; ilocal<solve_nlocal; ilocal++ )
      {
        nnz[ilocal] = 0;
        length_inz[ilocal] = INITIAL_SIZE;
        inz[ilocal] = get_new_int_short(INITIAL_SIZE);
        solve_A[ilocal] = get_new_dbl(INITIAL_SIZE);
        array_set( solve_A[ilocal], 0., INITIAL_SIZE );
      }
    }

    if ( band_solver )
    {
      ipiv = get_new_int( solve_nlocal );
      array_set( ipiv, 0, solve_nlocal );
      n = solve_nlocal;
      kl = band;
      ku = band;
      ldmat = 2*kl+ku+1;
      lmat = ldmat*solve_nlocal;
      mat = get_new_dbl( lmat );
      array_set( mat, 0., lmat );
      if ( swit )
      {
        pri( "ldmat", ldmat );
        pri( "kl", kl );
        pri( "ku", ku );
      }
    }

    if ( neigen>0 ) {
      lw = solve_nlocal;
      w = get_new_dbl( lw );
      array_set( w, 0., lw );
      ldz = solve_nlocal;
      lz = ldz*solve_nlocal;
      z = get_new_dbl( lz );
      array_set( z, 0., lz );
      lwork = 3*solve_nlocal;
      work = get_new_dbl( lwork );
      array_set( work, 0., lwork );
      ldmatlin = ku+1;
      lmatlin = ldmatlin*solve_nlocal;
      matlin = get_new_dbl( lmatlin );
      array_set( matlin, 0., lmatlin );
      if ( swit ) pri( "ldmatlin", ldmatlin );
      if ( control_eigen[0]==-GENERALIZED ) {
        ldmatss = ku+1;
        lmatss = ldmatss*solve_nlocal;
        matss = get_new_dbl( lmatss );
        array_set( matss, 0., lmatss );
        if ( swit ) pri( "ldmatss", ldmatss );
      }
    }

// fill matrices in band format or profile format
    if ( band_solver || petsc_solver || superlu_solver ) {
      for ( element=0; element<=max_element; element++ ) {
        if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) ) {
          element_group = 0;
          db( ELEMENT_GROUP, element, &element_group, ddum, ldum, 
            VERSION_NORMAL, GET_IF_EXISTS );
          test1 = db_active_index( ELEMENT_MATRIX_VALUES, element, VERSION_NORMAL );
          test2 = db_active_index( GROUP_MATRIX_VALUES, element_group, VERSION_NORMAL );
          if ( test1 || test2 ) {
            if ( solve_options_matrix_group==-YES ) {
              length = db_len( GROUP_MATRIX_VALUES, element_group, VERSION_NORMAL );
              element_matrix_values = 
                db_dbl( GROUP_MATRIX_VALUES, element_group, VERSION_NORMAL );
              if ( db_active_index( GROUP_MATRIX_SECOND_VALUES, 
                element_group, VERSION_NORMAL ) ) {
                element_matrix_second_values = 
                  db_dbl( GROUP_MATRIX_SECOND_VALUES, element_group, VERSION_NORMAL );
              }
            }
            else {
              length = db_len( ELEMENT_MATRIX_VALUES, element, VERSION_NORMAL );
              element_matrix_values = 
                db_dbl( ELEMENT_MATRIX_VALUES, element, VERSION_NORMAL );
              if ( db_active_index( ELEMENT_MATRIX_SECOND_VALUES, 
                element, VERSION_NORMAL ) ) {
                element_matrix_second_values = 
                  db_dbl( ELEMENT_MATRIX_SECOND_VALUES, element, VERSION_NORMAL );
              }
            }
            get_element_matrix_unknowns( element, element_matrix_unknowns );
            for ( imat=0; imat<length; imat++ ) {
              iglobal = element_matrix_unknowns[imat*2+0];
              jglobal = element_matrix_unknowns[imat*2+1];
              ilocal = solve_global_local[iglobal];
              jlocal = solve_global_local[jglobal];
              if ( ilocal!=-NO && jlocal!=-NO ) {
                if ( band_solver ) {
                  ii = kl + ku + 1 + (ilocal+1) - (jlocal+1);
                  jj = jlocal + 1;
                  indx = (jj-1)*ldmat + ii;
                  mat[indx-1] += element_matrix_values[imat];
                  if ( neigen>0 && jlocal>=ilocal ) {
                    ii = ku + 1 + (ilocal+1) - (jlocal+1);
                    jj = jlocal + 1;
                    indx = (jj-1)*ldmatlin + ii;
                    matlin[indx-1] += element_matrix_values[imat];
                    if ( control_eigen[0]==-GENERALIZED ) {
                      matss[indx-1] -= element_matrix_second_values[imat];
// prevent small negative diagonals due to limited numerical accuracy
                      if ( jlocal==ilocal ) matss[indx-1] += EPS_EIGEN; 
                    }
                  }
                }
                else {
                  assert( petsc_solver || superlu_solver );
                  if ( element_matrix_values[imat]!=0. ) {
                    if ( petsc_solver ) {
                      a = ilocal;
                      b = jlocal;
                    }
                    else {
                      assert( superlu_solver );
                      a = jlocal;
                      b = ilocal;
                    }
                    array_member( (long int *) inz[a], b, nnz[a], index );
                    if ( index<0 ) {
                      index = nnz[a];
                      nnz[a]++;
                      if ( nnz[a]>length_inz[a] ) {
                        int_ptr = get_new_int_short( length_inz[a] );
                        dbl_ptr = get_new_dbl( length_inz[a] );
                        for ( i=0; i<length_inz[a]; i++ ) {
                          int_ptr[i] = inz[a][i];
                          dbl_ptr[i] = solve_A[a][i];
                        }
                        delete[] inz[a];
                        delete[] solve_A[a];
                        inz[a] = get_new_int_short( 2*length_inz[a] );
                        solve_A[a] = get_new_dbl( 2*length_inz[a] );
                        for ( i=0; i<2*length_inz[a]; i++ ) {
                          inz[a][i] = 0;
                          solve_A[a][i] = 0.;
                        }
                        for ( i=0; i<length_inz[a]; i++ ) {
                          inz[a][i] = int_ptr[i];
                          solve_A[a][i] = dbl_ptr[i];
                        }
                        length_inz[a] *= 2;
                        delete[] int_ptr;
                        delete[] dbl_ptr;
                      }
                      inz[a][index] = b;
                    }
                    solve_A[a][index] += element_matrix_values[imat];
                  }
                }
              }
            }
          }
        }
      }
      for ( inod=0; inod<=max_node; inod++ ) {
        if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
          use_node_lhside = 0;
          if ( db_active_index( NODE_LHSIDE, inod, VERSION_NORMAL ) ) {
            use_node_lhside = 1;
            node_lhside = db_dbl( NODE_LHSIDE, inod, VERSION_NORMAL );
          }
          for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ )
          {
            iglobal = inod*npuknwn + ipuknwn;
            ilocal = solve_global_local[iglobal];
            jlocal = ilocal;
            if ( ilocal!=-NO ) {
              if ( band_solver ) {
                ii = kl + ku + 1 + (ilocal+1) - (jlocal+1);
                jj = jlocal + 1;
                indx = (jj-1)*ldmat + ii;
                if ( use_node_lhside ) mat[indx-1] += node_lhside[ipuknwn];
                if ( ilocal==jlocal && scalar_dabs(mat[indx-1])<EPS_MAT ) {
                  mat[indx-1] = EPS_MAT;
                }
                if ( neigen>0 ) {
                  ii = ku + 1 + (ilocal+1) - (jlocal+1);
                  jj = jlocal + 1;
                  indx = (jj-1)*ldmatlin + ii;
                  if ( use_node_lhside ) matlin[indx-1] += node_lhside[ipuknwn];
                }
              }
              else {
                assert( petsc_solver || superlu_solver );
                if ( petsc_solver ) {
                  a = ilocal;
                  b = jlocal;
                }
                else {
                  assert( superlu_solver );
                  a = jlocal;
                  b = ilocal;
                }              
                array_member( (long int *) inz[a], b, nnz[a], index );
                if ( index<0 ) {
                  index = nnz[a];
                  nnz[a]++;
                  if ( nnz[a]>length_inz[a] ) {
                    int_ptr = get_new_int_short( length_inz[a] );
                    dbl_ptr = get_new_dbl( length_inz[a] );
                    for ( i=0; i<length_inz[a]; i++ ) {
                      int_ptr[i] = inz[a][i];
                      dbl_ptr[i] = solve_A[a][i];
                    }
                    delete[] inz[a];
                    delete[] solve_A[a];
                    inz[a] = get_new_int_short( 2*length_inz[a] );
                    solve_A[a] = get_new_dbl( 2*length_inz[a] );
                    for ( i=0; i<2*length_inz[a]; i++ ) {
                      inz[a][i] = 0;
                      solve_A[a][i] = 0.;
                    }
                    for ( i=0; i<length_inz[a]; i++ ) {
                      inz[a][i] = int_ptr[i];
                      solve_A[a][i] = dbl_ptr[i];
                    }
                    length_inz[a] *= 2;
                    delete[] int_ptr;
                    delete[] dbl_ptr;
                  }
                  inz[a][index] = b;
                }
                if ( use_node_lhside ) {
                  solve_A[a][index] += node_lhside[ipuknwn];
                }
                if ( a==b && scalar_dabs(solve_A[a][index])<EPS_MAT ) {
                  solve_A[a][index] = EPS_MAT;
                }
              }
            }
          }
        }
      }
    }

// call solver
    if      ( bicg_solver ) 
      succesful = solve_iterative_bicg( );
    else if ( petsc_solver ) {
      if ( db( CONTROL_OPTIONS_SOLVER_PETSC_KSPTYPE, icontrol, &control_options_solver_petsc_ksptype, 
          ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS ) ) {
        ksptype = db_name( control_options_solver_petsc_ksptype );
      }
      if ( db( CONTROL_OPTIONS_SOLVER_PETSC_PCTYPE, icontrol, &control_options_solver_petsc_pctype, 
          ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS ) ) {
        pctype = db_name( control_options_solver_petsc_pctype );
      }
      if ( db( CONTROL_OPTIONS_SOLVER_PETSC_MG, icontrol, &mg, 
          ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS ) );
#if PETSC_USE
      succesful = solve_iterative_petsc( solve_A, solve_b, solve_nlocal, 
        nnz, inz, ksptype, pctype, mg, iter );
#endif
      array_move( solve_solution, solve_b, solve_nlocal );
      length = 1;
      db( GLOBAL_SOLVER_ITERATIONS, 0, &iter, ddum, length, VERSION_NORMAL, PUT );    
    }
    else if ( superlu_solver  ) {
      superlu_nnz = 0;
      for ( jlocal=0; jlocal<solve_nlocal; jlocal++ ) superlu_nnz += nnz[jlocal];
      superlu_A = get_new_dbl(superlu_nnz);
      superlu_asub = get_new_int_short(superlu_nnz);
      superlu_xa = get_new_int_short(solve_nlocal+1);
      for ( i=0; i<superlu_nnz; i++ ) superlu_A[i] = 0.;
      for ( i=0; i<superlu_nnz; i++ ) superlu_asub[i] = 0;
      for ( i=0; i<solve_nlocal+1; i++ ) superlu_xa[i] = 0;
      index = 0;
      for ( jlocal=0; jlocal<solve_nlocal; jlocal++ ) {
        for ( ilocal=0; ilocal<nnz[jlocal]; ilocal++ ) {
          assert(index<superlu_nnz);
          if ( ilocal==0 ) superlu_xa[jlocal] = index;
          superlu_A[index] = solve_A[jlocal][ilocal];
          superlu_asub[index] = inz[jlocal][ilocal];
          index++;
        }
        superlu_xa[solve_nlocal] = superlu_nnz;
      }
#if ( SUPERLU_USE || SUPERLU_MT_USE || SUPERLU_DIST_USE )
      succesful = solve_superlu( superlu_A, superlu_asub, superlu_xa, 
        superlu_nnz, solve_b, solve_nlocal, nthread );
#endif
    }
    else
      succesful = 1;
    if ( !succesful ) {
      pri( "The solver did not find a solution." );
      exit_tn_on_error();
    }
//+++++Lapack Solver Section
	  
if (band_solver) 
{
	n = solve_nlocal;
	cout << "number of equations/size of matrix= "<<n<<endl;
	nrhs= 1;
	dgbsv_( &n, &kl, &ku, &nrhs, mat, &ldmat, ipiv, solve_b, &n, &info );
	
	if ( info !=0 ) 
	{
		pri( "Error with Lapack band solver." );
		pri( "Try one of the other solver option:" );
		exit(TN_EXIT_STATUS);
	}
/*	cout<<"Value of info is "<<info<<"\n";	*/
}
//+++++Lapack Solver Section	  
    if ( neigen>0 ) {
      jobz[0] = 'V';
      uplo[0] = 'U';
      n = solve_nlocal;
      kd = ku;
      if ( solve_nlocal==1 ) {
        if ( control_eigen[0]==-MATRIX ) {
          w[0] = matlin[lmatlin-1];
          z[0] = 1;
        }
        else {
          assert( control_eigen[0]==-GENERALIZED );
          w[0] = matlin[lmatlin-1]/matss[lmatss-1];
          z[0] = 1;
        }
      }
      else {
        if ( control_eigen[0]==-MATRIX )
          dsbev_( jobz, uplo, &n, &kd, matlin, &ldmatlin, 
            w, z, &ldz, work, &info);
        else {
          assert( control_eigen[0]==-GENERALIZED );
          dsbgv_(jobz, uplo, &n, &kd, &kd, matlin, &ldmatlin, matss, 
            &ldmatss, w, z, &ldz, work, &info );
        }
        if ( info ) {
          pri( "Error with eigen solver." );
          pri( "Try one of the following:" );
          pri( "  Linear elements." );
          pri( "  Less elements." );
          pri( "  -matrix_direct for CONTROL_OPTIONS_SOLVER." );
          exit(TN_EXIT_STATUS);
        }
        for ( jlocal=0; jlocal<solve_nlocal; jlocal++ ) {
          max = 0.;
          for ( ilocal=0; ilocal<solve_nlocal; ilocal++ ) {
            tmp = scalar_dabs(z[jlocal*solve_nlocal+ilocal]);
            if ( tmp>max ) max = tmp;
          }
          if ( control_eigen[0]==-MATRIX || max>0. ) {
            for ( ilocal=0; ilocal<solve_nlocal; ilocal++ )
              z[jlocal*solve_nlocal+ilocal] *= control_eigen_scale/max;
          }
        }
      }
      if ( swit ) pri( "w", w, solve_nlocal );
      neigen_positive = 0;
      for ( ieigen=0; ieigen<neigen; ieigen++ ) {
        if ( control_eigen[0]==-MATRIX || w[ieigen]>0. ) {
          work[neigen_positive] = w[ieigen];
          neigen_positive++;
        }
      }
      if ( neigen_positive>0 ) db( CONTROL_EIGEN_VALUES, 0, idum, work, 
        neigen_positive, VERSION_NORMAL, PUT );
    }

// fill dofs
    for ( inod=0; inod<=max_node; inod++ ) {
      if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
        if ( bicg_solver || petsc_solver || superlu_solver || lapack_solver) 
          node_dof_new = db_dbl( NODE_DOF, inod, VERSION_NEW );
        for ( i=0; i<neigen_positive*nuknwn; i++ ) node_eigen[i] = 0.;
        iprinc=0;
        for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
          iuknwn = ipuknwn*nder;
          iglobal = inod*npuknwn + ipuknwn;
          if ( solve_global_local[iglobal]>=0 ) {
            ilocal = solve_global_local[iglobal];
            if ( ilocal!=-NO ) {
              if ( bicg_solver || petsc_solver || superlu_solver || lapack_solver) {
                tmp = solve_b[ilocal];
                if ( dof_principal[iuknwn]>=0 ) tmp *= options_relaxation[iprinc];
                node_dof_new[iuknwn] += tmp;
                if ( node_dof_new[iuknwn]<-EPS_SOL || node_dof_new[iuknwn]>EPS_SOL ) {
                  pri( "Error: the solution seems to diverge for unknown ", 
                    db_name(dof_label[iuknwn]) );
                  pri( "value", node_dof_new[iuknwn] );
                  exit_tn_on_error();
                }
              }
              neigen_positive = 0;
              for ( ieigen=0; ieigen<neigen; ieigen++ ) {
                if ( w[ieigen]>0. ) {
                  node_eigen[neigen_positive*nuknwn+iuknwn] = z[ieigen*solve_nlocal+ilocal];
                  neigen_positive++;
                }
              }
            }
          }
          if ( dof_principal[iuknwn]>=0 ) iprinc++;
        }
        if ( neigen_positive>0 ) {
          length = neigen_positive*nuknwn;
          db( NODE_EIGEN, inod, idum, node_eigen, length, VERSION_NORMAL, PUT );
        }
      }
    }

    if   ( band_solver ) 
	{
      delete[] mat;
      delete[] ipiv;
    }
    if ( petsc_solver ) 
	{
      delete[] solve_solution;
    }
    if ( petsc_solver || superlu_solver ) 
	{
      delete[] nnz;
      delete[] length_inz;
      for ( ilocal=0; ilocal<solve_nlocal; ilocal++ ) 
	  {
        delete[] solve_A[ilocal];
        delete[] inz[ilocal];
      }
      delete[] solve_A;
      delete[] inz;
    }
    if ( neigen>0 ) 
	{
      delete[] w;
      delete[] z;
      delete[] work;
      delete[] matlin;
      if ( control_eigen[0]==-GENERALIZED ) delete[] matss;
    }

  }

  delete[] element_matrix_unknowns;
  delete[] solve_b;
  delete[] solve_x;
  delete[] ordered_global;
  delete[] global_ordered;
  delete[] solve_global_local;
  delete[] ordering_has_been_done;

  delete[] node_bounded;
  delete[] dof_principal;
  delete[] node_node;
  delete[] dof_label;
  delete[] dof_type;
  delete[] node_eigen;
  delete[] options_relaxation;

  if ( swit ) pri( "Out routine SOLVE" );

}
