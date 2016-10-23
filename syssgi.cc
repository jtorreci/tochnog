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
#include "task.h"

void parallel_sys_initialize( void )

{
  long int nthread=0, ldum=0;
  double ddum[1];

  db( OPTIONS_PROCESSORS, 0, &nthread, ddum, ldum, VERSION_NORMAL, GET );

  m_set_procs( nthread );
}

void parallel_sys_lock( void )

{
  m_lock();
}

void parallel_sys_next_of_loop( long int next_of_loop[], long int max_loop,
  long int &nloop, long int &ithread )

{
  long int i=0, inext=0, ldum=0, nthread=0;
  double ddum[1];

  inext = m_next();

  db( OPTIONS_PROCESSORS, 0, &nthread, ddum, ldum, VERSION_NORMAL, GET );

  nloop = 1 + max_loop/nthread;
  for ( i=0; i<nloop; i++ ) {
    next_of_loop[i] = inext*nloop+i;
  }

  ithread = inext;
}

void parallel_sys_routine( void (*routine)(void) )

{
  parallel_active = 1;
  m_fork( routine );
  parallel_active = 0;
}

void parallel_sys_unlock( void )

{
  m_unlock();
}
