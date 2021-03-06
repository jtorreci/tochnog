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

void groundflow_data( long int element, long int gr, double old_unknowns[], 
  double new_unknowns[], double coord_ip[], double pe[], double &C )

{
  long int ldum=0;

  get_group_data( GROUP_GROUNDFLOW_PERMEABILITY, gr, element, new_unknowns,
    pe, ldum, GET_IF_EXISTS );
  get_group_data( GROUP_GROUNDFLOW_CAPACITY, gr, element, new_unknowns, 
    &C, ldum, GET_IF_EXISTS );

}
