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
#include "pthread.h"

pthread_mutex_t mutex;
long int inext=0;
pthread_t threads[MTHREAD];

void parallel_sys_initialize( void )

{
  pthread_mutex_init( &mutex, NULL );
}

void parallel_sys_lock( void )

{
  pthread_mutex_lock( &mutex );
}

void parallel_sys_next_of_loop( long int next_of_loop[], long int max_loop,
  long int &nloop, long int &ithread )

{
  long int i=0, nthread=0, ldum=0;
  double ddum[1];

  db( OPTIONS_PROCESSORS, 0, &nthread, ddum, ldum, VERSION_NORMAL, GET );

  parallel_sys_lock();

  inext++; 

  nloop = 1 + max_loop/nthread;
  for ( i=0; i<nloop; i++ ) {
    assert( i<=max_loop );
    next_of_loop[i] = inext*nloop + i;
  }

  ithread = inext;

  parallel_sys_unlock();

}

void parallel_sys_routine( void (*routine)(void) )

{
  long int ithread=0, nthread=0, icheck=0, ldum=0;
  double ddum[1];

  parallel_active = 1;
  inext = -1;

  db( OPTIONS_PROCESSORS, 0, &nthread, ddum, ldum, VERSION_NORMAL, GET );

  for ( ithread=0; ithread<nthread; ithread++ ) {
    icheck = pthread_create( &threads[ithread], NULL, 
      (void * (*)(void *)) routine, 0 );
    if ( icheck!=0 ) {
      cout << "Error in creating threads.\n";
      exit(TN_EXIT_STATUS);
    }
  }

  for ( ithread=0; ithread<nthread; ithread++ ) {
    icheck = pthread_join( threads[ithread], NULL );
    if ( icheck!=0 ) {
      cout << "Error in joining threads.\n";
      exit(TN_EXIT_STATUS);
    }
  }

  parallel_active = 0;

}

void parallel_sys_unlock( void )

{
  pthread_mutex_unlock( &mutex );
}
