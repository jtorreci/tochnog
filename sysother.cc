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

long int inext=0;

void parallel_sys_initialize( void )

{
  long int length=1, number_of_processors=1;
  double ddum[1];

  db( OPTIONS_PROCESSORS, 0, &number_of_processors, ddum, 
    length, VERSION_NORMAL, PUT );
}

void parallel_sys_lock( void )

{
}

void parallel_sys_next_of_loop( long int next_of_loop[], long int max_loop,
  long int &nloop, long int &ithread )

{
  long int i=0;

  inext++; 

  nloop = 1+max_loop;
  for ( i=0; i<nloop; i++ ) next_of_loop[i] = inext*nloop+i;

  ithread = inext;
}

void parallel_sys_routine( void (*routine)(void) )

{
  inext = -1;

  parallel_active = 1;
  routine();
  parallel_active = 0;
}

void parallel_sys_unlock( void )

{
}
