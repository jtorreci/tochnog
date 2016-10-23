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

#if MPI_USE
#include <mpi.h>
#endif


#if PETSC_USE
 static char help[] = "";
 static char out[] = "";
 extern "C"
  int PetscInitialize(int *argc,char ***argv,char *file_name,char *help_message);
 extern "C"
  int PetscFinalize(void);
#endif

int main( int argc, char* argv[] )

{
    //added for CPU time
  Time CPU;
  CPU.firsttime=CPU.taketime();

  long int i=0, l=0, any_point=0;

    // initialise static variables
  initialize();

    // test arguments
#if MPI_USE
  strcpy( data_file, "tn.dat" );
#else
  if ( argc==1 ) strcpy( data_file, "tn.dat" );
  else if ( argc==2 ) strcpy( data_file, argv[1] );
  else {
    pri( "Usage: tochnog" );
    pri( "or:    tochnog file.dat" );
    pri( "or:    tochnog file.dat > file.out" );
    exit(TN_EXIT_STATUS);
  }
#endif

    // append .dat to input file name
  l = strlen( data_file );
  strcpy( data_file_base, "" );
  for ( i=0; i<l && !any_point; i++ ) {
    if ( data_file[i]=='.' ) any_point = 1;
    else strncat( data_file_base, &data_file[i], 1 );
  }
  if ( !any_point ) strcat( data_file, ".dat" );

    // empty the tn.dvd file at the start of calculation
  ofstream outdvd( "tn.dvd" );
  outdvd.close();

    // read input file
  input();

    // initialize petsc solver
#if PETSC_USE
  PetscInitialize(&argc,&argv,out,help);
#endif
    // initialize MPI
#if MPI_USE
  MPI_Init(&argc, &argv );
#endif

    // fe analysis
  top();

    // finalize petsc solver
#if PETSC_USE
  PetscFinalize();
#endif

    // finalize MPI
#if MPI_USE
  MPI_Finalize();
#endif

    //added for CPU time
  CPU.secondtime=CPU.taketime();
  CPU.CPUtime();

    // round up things
  exit_tn( -YES );

    // and exit
  return 0;
}
