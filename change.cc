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


void change_geometry( long int task, double dtime, double time_current )

{
  long int i=0, ichange=0, max_change=0, geometry_name=0, 
    geometry_index=0, igeom=0, ngeom=0, operat=0, 
    index=0, nindex=0, length=0, ldum=0, change_geometry_time_user=0,
    idum[1], change_geometry[3], *indices=NULL, *geometry_entities=NULL;
  double val=0., ddum[1], *dval=NULL;

  db_max_index( CHANGE_GEOMETRY, max_change, VERSION_NORMAL, GET );
  if ( max_change>=0 ) {
    set_swit(-1,-1,"change_geometry");
    indices = get_new_int(DATA_ITEM_SIZE);
    geometry_entities = get_new_int(DATA_ITEM_SIZE);
    dval = get_new_dbl(DATA_ITEM_SIZE);
    for ( ichange=0; ichange<=max_change; ichange++ ) {
      if ( db_active_index( CHANGE_GEOMETRY, ichange, VERSION_NORMAL ) ) {
        db( CHANGE_GEOMETRY, ichange, change_geometry, 
          ddum, ldum, VERSION_NORMAL, GET );
        if ( change_geometry[0]==-GEOMETRY_SET ) {
          db( GEOMETRY_SET, change_geometry[1], geometry_entities, 
            ddum, ngeom, VERSION_NORMAL, GET );
          ngeom /= 2;
        }
        else {
          geometry_entities[0] = change_geometry[0];
          geometry_entities[1] = change_geometry[1];
          ngeom = 1;
        }
        operat = change_geometry[2];
        db( CHANGE_GEOMETRY_TIME_USER, ichange, &change_geometry_time_user, 
          ddum, ldum, VERSION_NORMAL, GET );
        if ( change_geometry_time_user==-YES )
          user_change_geometry_time( ichange, time_current, val );
        else {
          pri( "Error: CHANGE_GEOMETRY_TIME_USER should be set to -YES" );
          exit(TN_EXIT_STATUS);
        }
        for ( igeom=0; igeom<ngeom; igeom++ ) {
          geometry_name = geometry_entities[igeom*2+0];
          geometry_index = geometry_entities[igeom*2+1];
          db( geometry_name, geometry_index, idum, dval, length, VERSION_NORMAL, GET );
          if      ( geometry_name==-GEOMETRY_BRICK ) {
            nindex = 1;
            indices[0] = 1*ndim-1;
          }
          else if ( geometry_name==-GEOMETRY_CIRCLE ) {
            nindex = 1;
            indices[0] = 1*ndim-1;
          }
          else if ( geometry_name==-GEOMETRY_CIRCLE_SEGMENT ) {
            nindex = 1;
            indices[0] = 1;
          }
          else if ( geometry_name==-GEOMETRY_CYLINDER ) {
            nindex = 2;
            indices[0] = 1*ndim-1;
            indices[1] = 2*ndim-1;
          }
          else if ( geometry_name==-GEOMETRY_CYLINDER_SEGMENT ) {
            nindex = 2;
            indices[0] = 1*ndim-1;
            indices[1] = 2*ndim-1;
          }
          else if ( geometry_name==-GEOMETRY_ELLIPSE ) {
            nindex = 1;
            indices[0] = ndim-1;
          }
          else if ( geometry_name==-GEOMETRY_LINE ) {
            nindex = 2;
            indices[0] = 1*ndim-1;
            indices[1] = 2*ndim-1;
          }
          else if ( geometry_name==-GEOMETRY_POINT ) {
            nindex = 1;
            indices[0] = ndim-1;
          }
          else if ( geometry_name==-GEOMETRY_QUADRILATERAL ) {
            nindex = 4;
            indices[0] = 1*ndim-1;
            indices[1] = 2*ndim-1;
            indices[2] = 3*ndim-1;
            indices[3] = 4*ndim-1;
          }
          else if ( geometry_name==-GEOMETRY_SPHERE ) {
            nindex = 1;
            indices[0] = ndim-1;
          }
          else if ( geometry_name==-GEOMETRY_SPHERE_SEGMENT ) {
            nindex = 1;
            indices[0] = ndim-1;
          }
          else if ( geometry_name==-GEOMETRY_TRIANGLE ) {
            nindex = 3;
            indices[0] = 1*ndim-1;
            indices[1] = 2*ndim-1;
            indices[2] = 3*ndim-1;
          }
          else {
            pri( "Sorry: CHANGE_GEOMETRY not available for ", geometry_name );
            exit(TN_EXIT_STATUS);
          }
          for ( i=0; i<nindex; i++ ) {
            index = indices[i];
            assert(index<length);
            if      ( operat==-USE )
              dval[index] = val;
            else if ( operat==-ADD ) {
              if ( task==YES ) {
                dval[index] += val*dtime;
              }
            }
          }
          db( geometry_name, geometry_index, idum, dval, length, VERSION_NORMAL, PUT );
        }
      }
    }
    delete[] indices;
    delete[] geometry_entities;
    delete[] dval;
  }

}
