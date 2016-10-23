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

/* #include <petscksp.h> */

#if PETSC_USE 

ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);

extern double *solve_solution;

long int solve_iterative_petsc( double *solve_petsc_A[], double solve_b[],
  long int solve_nlocal, int *nnz, int *inz[], 
  char *ksptype, char *pctype, long int mg, long int *iter ) 

{
  Vec     x, b, xtot;   /* approx solution, RHS, exact solution */
  Mat     A;            /* linear system matrix */
  PetscErrorCode ierr;  /* Petsc error code */
  PC      pc;           /* preconditioner context */
  KSP     ksp;          /* Krylov subspace method context */
  VecScatter scatter;   /* scatter context */ 
  IS      from, to;     /* index sets that define the scatter */ 
  int     its=0, i=0, j=0;
  int     *cols, *idx_from, *idx_to ,*loc_size; 
  int     nlocal;
  
/* 
  initialize (PetscInitialize and PetscFinalize are called in main)
*/
  
/* 
  Create vectors 
*/
  nlocal = solve_nlocal;
  ierr = VecCreateMPI( PETSC_COMM_WORLD, 
    PETSC_DECIDE, nlocal, &x ); CHKERRA(ierr);
  loc_size=(int *)PetscMalloc(sizeof(int));  
  ierr = VecGetLocalSize( x, loc_size); CHKERRA(ierr);
  ierr = VecCreateMPIWithArray( PETSC_COMM_WORLD, 
    *loc_size, nlocal, solve_b, &b ); CHKERRA(ierr);
  ierr = VecAssemblyBegin(b); CHKERRA(ierr);
  ierr = VecAssemblyEnd(b); CHKERRA(ierr);        

/*
  Create matrix 
*/
  cols = (int *)PetscMalloc(nlocal*sizeof(int)); CHKPTRA(cols);
  ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD, *loc_size, *loc_size,
    nlocal, nlocal, 0, PETSC_NULL, 0, nnz, &A ); CHKERRA(ierr);
  for (i=0; i<nlocal; i++ ) {
    for ( j=0; j<nnz[i]; j++ ) cols[j] = inz[i][j];
    ierr = MatSetValues(A,1,&i,nnz[i],cols,
      solve_petsc_A[i],INSERT_VALUES); CHKERRA(ierr);
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRA(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRA(ierr);
  PetscFree(cols);

/*
  Create linear solver
*/
  ierr = SLESCreate(PETSC_COMM_WORLD,&sles); CHKERRA(ierr);
  ierr = SLESSetOperators(sles,A,A,DIFFERENT_NONZERO_PATTERN); CHKERRA(ierr);

/*
  solver options
*/
  ierr = SLESGetKSP(sles,&ksp); CHKERRA(ierr);
  ierr = SLESGetPC(sles,&pc); CHKERRA(ierr);
  if ( strcmp(ksptype,"none" ) ) {
    KSPSetType(ksp,ksptype); CHKERRA(ierr);
  }
  if ( strcmp(pctype,"none" ) ) {
    PCSetType(pc,pctype); CHKERRA(ierr);
  }
  if ( mg>0 ) {
    ierr = PCSetType(pc,PCMG); CHKERRA(ierr);
/* ierr = MGSetLevels( pc, nlevels, mg ); */
    ierr = MGSetType( pc, mg );
  }
  ierr = KSPSetTolerances(ksp,1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,
    10*nlocal ); CHKERRA(ierr);
  ierr = SLESSetFromOptions(sles); CHKERRA(ierr);

/*
  Solve
*/
  ierr = SLESSolve(sles,b,x,&its); CHKERRA(ierr); 
  

/*
  Get solution vector
*/
  idx_from = (int *)PetscMalloc(nlocal*sizeof(int)); CHKPTRA(idx_from);
  idx_to = (int *)PetscMalloc(nlocal*sizeof(int)); CHKPTRA(idx_to);
  for (i=0; i<nlocal; i++ ) {
    idx_from[i] = i;
    idx_to[i] = i;
  }
  VecCreateSeq(PETSC_COMM_SELF,nlocal,&xtot);
  ISCreateGeneral(PETSC_COMM_SELF,nlocal,idx_from,&from); 
  ISCreateGeneral(PETSC_COMM_SELF,nlocal,idx_to,&to); 
  VecScatterCreate(x,from,xtot,to,&scatter); 
  VecScatterBegin(x,xtot,INSERT_VALUES,SCATTER_FORWARD,scatter); 
  VecScatterEnd(x,xtot,INSERT_VALUES,SCATTER_FORWARD,scatter); 
  VecGetArray(xtot,&solve_b); 
  for (i=0; i<nlocal; i++ ) solve_solution[i] = solve_b[i];
  ISDestroy(from); 
  ISDestroy(to);  
  VecScatterDestroy(scatter); 
  PetscFree(idx_from);
  PetscFree(idx_to);
  PetscFree(loc_size);
  VecDestroy(xtot);

/*
  Clean
*/
  ierr = VecDestroy(x); CHKERRA(ierr); 
  ierr = VecDestroy(b); CHKERRA(ierr); 
  ierr = MatDestroy(A); CHKERRA(ierr);
  ierr = SLESDestroy(sles); CHKERRA(ierr);

  if ( its>=0 ) {
    *iter = its;
    return 1;
  }
  else {
    *iter = -its;
    return 0;
  }
}

#endif
