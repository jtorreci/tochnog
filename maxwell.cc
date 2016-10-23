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

void maxwell( long int type, long int element, 
  long int gr, long int nnol, double volume,
  double new_unknowns[], 
  double old_dof[], double new_dof[],
  double h[], double d[], double element_lhside[],
  double element_matrix[], double element_matrix_second[],
  double element_rhside[] )         

{
  pri( "maxwell is not available." );
  exit(TN_EXIT_STATUS);
}

void maxwell_scatter( void )

{
}
