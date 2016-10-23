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
    along with this program; if not, write to the Free Software
    59 Temple Place, Suite 330, Boston, MA, 02111-1307, USA
*/

#include "tochnog.h"

#define MNOD_MACRO 27
#define MEL_MACRO 8
#define EPSPHI 1.e-5

void macro( void )

{
  long int i=0, n=0, nthick=0, ncirc=0, ny=0, nz=0, nx=0, icontrol=0,
    ix=0, ithick=0, icirc=0, iy=0, iz=0, inol=0, inod=0, macro_item=0, 
    element=0, idim=0,
    length_mesh_macro=0, length_mesh_macro_parameters=0, isnot360=0,
    use_control_refine_globally_geometry=0, nnod_macro=0, nel_macro=0,
    length_control_refine_globally=0, max_geometry=0,ncircx=0,
    node_boundary=-YES, length=0, nnol=0, element_group=0,
    ok=0, geometry_macro_name=0, swit=0,
    test1=0, test2=0, test3=0,
    control_mesh_macro_element=0,
    control_mesh_macro_set_node_boundary=-YES,
    ldum=0, idum[1], node_remesh_allowed[MDIM], nodes[MNOL], el[1+MNOL], 
    control_refine_globally[4], control_refine_globally_geometry[2],
    nodes_macro[MNOL][MEL_MACRO], *control_mesh_macro=NULL;
  double radius=0., thickness=0., lenx=0., leny=0., lenz=0.,
    r=0., phi=0., dx=0., dy=0., dz=0., phi0=0., phi1=2.0*PIRAD,
    xuij=0., yuij=0., zuij=0., ij1=0., alpha=0., beta=0., ONE=0., TWO=0., THREE=0.,
    ddum[1], work[MDIM], coord[MDIM], coordloc[MDIM], rot[MDIM*MDIM],
    xc[MDIM], yc[MDIM], zc[MDIM], xc0[MDIM], yc0[MDIM], zc0[MDIM], 
    node_dof[MUKNWN], coord_macro[MDIM][MNOD_MACRO], 
    geometry_macro[MDIM+2], *control_mesh_macro_parameters=NULL;

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  if ( db_active_index( CONTROL_MESH_MACRO, icontrol, VERSION_NORMAL ) ) {
    swit = set_swit(-1,-1,"macro");
    if ( swit ) pri( "In routine MACRO" );
    db( CONTROL_MESH_MACRO_SET_NODE_BOUNDARY, icontrol, 
      &control_mesh_macro_set_node_boundary, ddum, ldum, VERSION_NORMAL,
      GET_IF_EXISTS );
    control_mesh_macro = 
      db_int( CONTROL_MESH_MACRO, icontrol, VERSION_NORMAL );
    length_mesh_macro = db_len( CONTROL_MESH_MACRO, icontrol, VERSION_NORMAL );
    length_mesh_macro_parameters = db_len( CONTROL_MESH_MACRO_PARAMETERS, 
      icontrol, VERSION_NORMAL );
    macro_item = control_mesh_macro[0];
    element_group = control_mesh_macro[1];

    ok = 0;
    if      ( macro_item==-BAR && length_mesh_macro==3 ) ok = 1;
    else if ( macro_item==-RECTANGLE && length_mesh_macro==4 ) ok = 1;
    else if ( macro_item==-CIRCLE && length_mesh_macro==3 ) ok = 1;
    else if ( macro_item==-CIRCLE_HOLLOW && length_mesh_macro==4 ) ok = 1;
    else if ( macro_item==-BRICK && length_mesh_macro==5 ) ok = 1;
    else if ( macro_item==-CYLINDER_HOLLOW && length_mesh_macro==5 ) ok = 1;
    else if ( macro_item==-SPHERE && length_mesh_macro==3 ) ok = 1;
    if ( !ok ) db_error( CONTROL_MESH_MACRO, icontrol );

    ok = 0;
    if      ( macro_item==-BAR && length_mesh_macro_parameters==2 ) ok = 1;
    else if ( macro_item==-RECTANGLE && length_mesh_macro_parameters==4 ) ok = 1;
    else if ( macro_item==-CIRCLE && length_mesh_macro_parameters==3 ) ok = 1;
    else if ( macro_item==-CIRCLE_HOLLOW && length_mesh_macro_parameters==6 ) ok = 1;
    else if ( macro_item==-BRICK && length_mesh_macro_parameters==6 ) ok = 1;
    else if ( macro_item==-CYLINDER_HOLLOW && length_mesh_macro_parameters==10 ) ok = 1;
    else if ( macro_item==-SPHERE && length_mesh_macro_parameters==4 ) ok = 1;
    if ( !ok ) db_error( CONTROL_MESH_MACRO_PARAMETERS, icontrol );

    control_mesh_macro_parameters = 
      db_dbl( CONTROL_MESH_MACRO_PARAMETERS, icontrol, VERSION_NORMAL );
    if ( macro_item==-CYLINDER_HOLLOW ) {
      lenx = array_distance( &control_mesh_macro_parameters[0], 
        &control_mesh_macro_parameters[ndim], work, ndim );
      xc[0] = control_mesh_macro_parameters[0];
      yc[0] = control_mesh_macro_parameters[1];
      zc[0] = control_mesh_macro_parameters[2];
      xc[1] = control_mesh_macro_parameters[3];
      yc[1] = control_mesh_macro_parameters[4];
      zc[1] = control_mesh_macro_parameters[5];
      array_set( xc0, 0., 2 );
      array_set( yc0, 0., 2 );
      array_set( zc0, 0., 2 );
      xc0[1] = lenx;
      xuij = (xc[1]-xc0[1]) - (xc[0]-xc0[0]);
      yuij = (yc[1]-yc0[1]) - (yc[0]-yc0[0]);
      zuij = (zc[1]-zc0[1]) - (zc[0]-zc0[0]);
      ij1 = sqrt( scalar_square(lenx+xuij) + scalar_square(zuij) );
      if ( scalar_dabs(ij1)<1.e-8 )
        alpha = PIRAD/2.;
      else
        alpha = acos( (lenx+xuij)/ij1 );
      beta = asin( yuij/lenx );
      rot[0] = cos(alpha)*cos(beta);
      rot[1] = sin(beta);
      rot[2] = sin(alpha)*cos(beta);
      rot[3] = -cos(alpha)*sin(beta);
      rot[4] = cos(beta);
      rot[5] = -sin(alpha)*sin(beta);
      rot[6] = -sin(alpha);
      rot[7] = 0.;
      rot[8] = cos(alpha);
    }
    if      ( macro_item==-CIRCLE_HOLLOW || macro_item==-CYLINDER_HOLLOW ) {
      if ( macro_item==-CIRCLE_HOLLOW ) {
        nx = 1;
        phi0=control_mesh_macro_parameters[4]*PIRAD/180.0;
        phi1=control_mesh_macro_parameters[5]*PIRAD/180.0;
        nthick = control_mesh_macro[2];
        if ( nthick<2 ) db_error( CONTROL_MESH_MACRO, icontrol );
        ncirc = control_mesh_macro[3];
        if ( ncirc<2 ) db_error( CONTROL_MESH_MACRO, icontrol );
        if ( scalar_dabs(phi1-phi0-2.0*PIRAD)<EPSPHI )
          isnot360 = 0; 
        else 
          isnot360 = 1;
        ncircx = ncirc + isnot360;
        lenx = 0.;
        radius = control_mesh_macro_parameters[2];
        thickness = control_mesh_macro_parameters[3];
      }
      else {
        assert( macro_item==-CYLINDER_HOLLOW );
        nx = control_mesh_macro[2];
        if ( nx<2 ) db_error( CONTROL_MESH_MACRO, icontrol );
        nthick = control_mesh_macro[3];
        if ( nthick<2 ) db_error( CONTROL_MESH_MACRO, icontrol );
        ncirc = control_mesh_macro[4];
        if ( ncirc<2 ) db_error( CONTROL_MESH_MACRO, icontrol );
        radius = control_mesh_macro_parameters[6];
        thickness = control_mesh_macro_parameters[7];
        phi0=control_mesh_macro_parameters[8]*PIRAD/180.0;
        phi1=control_mesh_macro_parameters[9]*PIRAD/180.0;
        if ( scalar_dabs(phi1-phi0-2.0*PIRAD)<EPSPHI ) 
          isnot360 = 0; 
        else 
          isnot360 = 1;
        ncircx = ncirc + isnot360;
      }
      for ( ix=0; ix<nx; ix++ ) {
        dx = ix*lenx/(nx-1);
        for ( icirc=0; icirc<ncircx; icirc++ ) {
          phi = icirc*(phi1-phi0)/(ncirc)+phi0;
          for ( ithick=0; ithick<nthick; ithick++ ) {
            r = radius - thickness/2. + ithick*thickness/(nthick-1);
            if      ( ndim==2 ) {
              coordloc[0] = r*cos(phi);
              coordloc[1] = r*sin(phi);
              array_move( coordloc, coord, ndim );
            }
            else {
              assert( ndim==3 );
              coordloc[0] = dx;
              coordloc[1] = r*cos(phi);
              coordloc[2] = r*sin(phi);
              matrix_atb( rot, coordloc, coord, ndim, ndim, 1 );
            }              
            array_add( coord, control_mesh_macro_parameters, coord, ndim );
            db( NODE, inod, idum, coord, ndim, VERSION_MACRO, PUT );
            db( NODE_START_REFINED, inod, idum, coord, ndim, VERSION_MACRO, PUT );
            length = 1; db( NODE_MACRO_GENERATE, inod, &icontrol, ddum, length, VERSION_MACRO, PUT );
            if ( nuknwn>0 ) {
              array_set( node_dof, 0., nuknwn );
              db( NODE_DOF, inod, idum, node_dof, nuknwn, VERSION_MACRO, PUT );
              db( NODE_DOF_START_REFINED, inod, idum, node_dof, nuknwn, VERSION_MACRO, PUT );
            }
            if ( ithick==0 || ithick==nthick-1 || (isnot360==1 && (icirc==0 || icirc==ncirc)) ||
                ( ndim==3 && (ix==0 || ix==nx-1) ) ) {
              if (  control_mesh_macro_set_node_boundary==-YES) {
                length = 1; db( NODE_BOUNDARY, inod, &node_boundary, ddum, 
                  length, VERSION_MACRO, PUT );
              }
            }
            inod++;
          }
        }
      }
      for ( ix=0; ix<scalar_imax(1,nx-1); ix++ ) {
        for ( icirc=0; icirc<ncirc; icirc++ ) {
          for ( ithick=0; ithick<nthick-1; ithick++ ) {
            if ( ndim==2 ) {
              el[0] = -QUAD4;
              nnol = 4;
              nodes[0] = icirc*nthick+ithick;
              nodes[1] = icirc*nthick+ithick + 1;
              if ( icirc==ncirc-1 && isnot360==0 ) {
                nodes[2] = 0*nthick+ithick;
                nodes[3] = 0*nthick+ithick + 1;
              }
              else {
                nodes[2] = (icirc+1)*nthick+ithick;
                nodes[3] = (icirc+1)*nthick+ithick + 1;
              }
            }
            else {
              assert( ndim==3 );
              el[0] = -HEX8;
              nnol = 8;
              nodes[0] = ix*ncircx*nthick+icirc*nthick+ithick;
              nodes[1] = ix*ncircx*nthick+icirc*nthick+ithick + 1;
              nodes[4] = (ix+1)*ncircx*nthick+icirc*nthick+ithick;
              nodes[5] = (ix+1)*ncircx*nthick+icirc*nthick+ithick + 1;
              if ( icirc==ncirc-1 && isnot360==0 ) {
                nodes[2] = ix*ncircx*nthick+0*nthick+ithick;
                nodes[3] = ix*ncircx*nthick+0*nthick+ithick + 1;
                nodes[6] = (ix+1)*ncircx*nthick+0*nthick+ithick;
                nodes[7] = (ix+1)*ncircx*nthick+0*nthick+ithick + 1;
              }
              else {
                nodes[2] = ix*ncircx*nthick+(icirc+1)*nthick+ithick;
                nodes[3] = ix*ncircx*nthick+(icirc+1)*nthick+ithick + 1;
                nodes[6] = (ix+1)*ncircx*nthick+(icirc+1)*nthick+ithick;
                nodes[7] = (ix+1)*ncircx*nthick+(icirc+1)*nthick+ithick + 1;
              }
            }
            
            array_move( nodes, &el[1], nnol ); 
            length = 1 + nnol;
            db( ELEMENT, element, el, ddum, length, VERSION_MACRO, PUT );
            length = 1;
            db( ELEMENT_GROUP, element, &element_group, ddum, length, VERSION_MACRO, PUT );
            length = 1;
            db( ELEMENT_MACRO_GENERATE, element, &icontrol, ddum, length, VERSION_MACRO, PUT );
            element++;
          }
        }
      }
    }
    else if ( macro_item==-SPHERE || macro_item==-CIRCLE ) {
      ONE = 1.;
      TWO = 1./sqrt(2.);
      THREE = 1./sqrt(3.);
      n = control_mesh_macro[2];
      if ( n<0 ) db_error( CONTROL_MESH_MACRO, icontrol );
      radius = control_mesh_macro_parameters[ndim];
      array_set( &coord_macro[0][0], 0., MDIM*MNOD_MACRO );
      if ( macro_item==-SPHERE ) {
        nnod_macro = 27;
        nel_macro = 8;
        coord_macro[0][0]   = -THREE; coord_macro[1][0]   = -THREE; coord_macro[2][0]   = -THREE;
        coord_macro[0][1]   =  0.;    coord_macro[1][1]   = -TWO;   coord_macro[2][1]   = -TWO;
        coord_macro[0][2]   = +THREE; coord_macro[1][2]   = -THREE; coord_macro[2][2]   = -THREE;
        coord_macro[0][3]   = -TWO;   coord_macro[1][3]   =  0.;    coord_macro[2][3]   = -TWO;
        coord_macro[0][4]   =  0.;    coord_macro[1][4]   = -0.;    coord_macro[2][4]   = -ONE;
        coord_macro[0][5]   = +TWO;   coord_macro[1][5]   =  0.;    coord_macro[2][5]   = -TWO;
        coord_macro[0][6]   = -THREE; coord_macro[1][6]   = +THREE; coord_macro[2][6]   = -THREE;
        coord_macro[0][7]   =  0.;    coord_macro[1][7]   = +TWO;   coord_macro[2][7]   = -TWO;
        coord_macro[0][8]   = +THREE; coord_macro[1][8]   = +THREE; coord_macro[2][8]   = -THREE;

        coord_macro[0][9]   = -TWO;   coord_macro[1][9]   = -TWO;   coord_macro[2][9]   =  0.;
        coord_macro[0][10]  =  0.;    coord_macro[1][10]  = -ONE;   coord_macro[2][10]  =  0.;
        coord_macro[0][11]  = +TWO;   coord_macro[1][11]  = -TWO;   coord_macro[2][11]  =  0.;
        coord_macro[0][12]  = -ONE;   coord_macro[1][12]  =  0.;    coord_macro[2][12]  =  0.;
        coord_macro[0][13]  =  0.;    coord_macro[1][13]  =  0.;    coord_macro[2][13]  =  0.;
        coord_macro[0][14]  = +ONE;   coord_macro[1][14]  =  0.;    coord_macro[2][14]  =  0.;
        coord_macro[0][15]  = -TWO;   coord_macro[1][15]  = +TWO;   coord_macro[2][15]  =  0.;
        coord_macro[0][16]  =  0.;    coord_macro[1][16]  = +ONE;   coord_macro[2][16]  =  0.;
        coord_macro[0][17]  = +TWO;   coord_macro[1][17]  = +TWO;   coord_macro[2][17]  =  0.;

        coord_macro[0][18]  = -THREE; coord_macro[1][18]  = -THREE; coord_macro[2][18]  = +THREE;
        coord_macro[0][19]  =  0.;    coord_macro[1][19]  = -TWO;   coord_macro[2][19]  = +TWO;
        coord_macro[0][20]  = +THREE; coord_macro[1][20]  = -THREE; coord_macro[2][20]  = +THREE;
        coord_macro[0][21]  = -TWO;   coord_macro[1][21]  =  0.;    coord_macro[2][21]  = +TWO;
        coord_macro[0][22]  =  0.;    coord_macro[1][22]  =  0.;    coord_macro[2][22]  = +ONE;
        coord_macro[0][23]  = +TWO;   coord_macro[1][23]  =  0.;    coord_macro[2][23]  = +TWO;
        coord_macro[0][24]  = -THREE; coord_macro[1][24]  = +THREE; coord_macro[2][24]  = +THREE;
        coord_macro[0][25]  =  0.;    coord_macro[1][25]  = +TWO;   coord_macro[2][25]  = +TWO;
        coord_macro[0][26]  = +THREE; coord_macro[1][26]  = +THREE; coord_macro[2][26]  = +THREE;

        nodes_macro[0][0]  = 0;  nodes_macro[1][0]  = 3;  nodes_macro[2][0]  = 1;  nodes_macro[3][0]  = 4;
        nodes_macro[4][0]  = 9;  nodes_macro[5][0]  = 12; nodes_macro[6][0]  = 10; nodes_macro[7][0]  = 13;

        nodes_macro[0][1]  = 1;  nodes_macro[1][1]  = 4;  nodes_macro[2][1]  = 2;  nodes_macro[3][1]  = 5;
        nodes_macro[4][1]  = 10; nodes_macro[5][1]  = 13; nodes_macro[6][1]  = 11; nodes_macro[7][1]  = 14;

        nodes_macro[0][2]  = 3;  nodes_macro[1][2]  = 6;  nodes_macro[2][2]  = 4;  nodes_macro[3][2]  = 7;
        nodes_macro[4][2]  = 12; nodes_macro[5][2]  = 15; nodes_macro[6][2]  = 13; nodes_macro[7][2]  = 16;

        nodes_macro[0][3]  = 4;  nodes_macro[1][3]  = 7;  nodes_macro[2][3]  = 5;  nodes_macro[3][3]  = 8;
        nodes_macro[4][3]  = 13; nodes_macro[5][3]  = 16; nodes_macro[6][3]  = 14; nodes_macro[7][3]  = 17;

        nodes_macro[0][4]  = 9;  nodes_macro[1][4]  = 12; nodes_macro[2][4]  = 10; nodes_macro[3][4]  = 13;
        nodes_macro[4][4]  = 18; nodes_macro[5][4]  = 21; nodes_macro[6][4]  = 19; nodes_macro[7][4]  = 22;

        nodes_macro[0][5]  = 10; nodes_macro[1][5]  = 13; nodes_macro[2][5]  = 11; nodes_macro[3][5]  = 14;
        nodes_macro[4][5]  = 19; nodes_macro[5][5]  = 22; nodes_macro[6][5]  = 20; nodes_macro[7][5]  = 23;

        nodes_macro[0][6]  = 12; nodes_macro[1][6]  = 15; nodes_macro[2][6]  = 13; nodes_macro[3][6]  = 16;
        nodes_macro[4][6]  = 21; nodes_macro[5][6]  = 24; nodes_macro[6][6]  = 22; nodes_macro[7][6]  = 25;

        nodes_macro[0][7]  = 13; nodes_macro[1][7]  = 16; nodes_macro[2][7]  = 14; nodes_macro[3][7]  = 17;
        nodes_macro[4][7]  = 22; nodes_macro[5][7]  = 25; nodes_macro[6][7]  = 23; nodes_macro[7][7]  = 26;

      }
      else {
        assert( macro_item==-CIRCLE );
        nnod_macro = 9;
        nel_macro = 4;

        coord_macro[0][0]   = -TWO; coord_macro[1][0]   = -TWO;
        coord_macro[0][1]   =   0.; coord_macro[1][1]   = -ONE;
        coord_macro[0][2]   = +TWO; coord_macro[1][2]   = -TWO;

        coord_macro[0][3]   = -ONE; coord_macro[1][3]   =    0.;
        coord_macro[0][4]   =   0.; coord_macro[1][4]   =    0.;
        coord_macro[0][5]   = +ONE; coord_macro[1][5]   =    0.;

        coord_macro[0][6]   = -TWO; coord_macro[1][6]   = +TWO;
        coord_macro[0][7]   =   0.; coord_macro[1][7]   = +ONE;
        coord_macro[0][8]   = +TWO; coord_macro[1][8]   = +TWO;

        nodes_macro[0][0]  = 0;  nodes_macro[1][0]  = 1;  nodes_macro[2][0]  = 3;  nodes_macro[3][0]  = 4;

        nodes_macro[0][1]  = 1;  nodes_macro[1][1]  = 2;  nodes_macro[2][1]  = 4;  nodes_macro[3][1]  = 5;

        nodes_macro[0][2]  = 3;  nodes_macro[1][2]  = 4;  nodes_macro[2][2]  = 6;  nodes_macro[3][2]  = 7;

        nodes_macro[0][3]  = 4;  nodes_macro[1][3]  = 5;  nodes_macro[2][3]  = 7;  nodes_macro[3][3]  = 8;
      }
      array_multiply( &coord_macro[0][0], &coord_macro[0][0], radius, MDIM*MNOD_MACRO );

      for ( inod=0; inod<nnod_macro; inod++ ) {
        for ( idim=0; idim<ndim; idim++ ) coordloc[idim] = coord_macro[idim][inod];
        array_move( coordloc, coord, ndim );
        array_add( coord, control_mesh_macro_parameters, coord, ndim );
        db( NODE, inod, idum, coord, ndim, VERSION_MACRO, PUT );
        db( NODE_START_REFINED, inod, idum, coord, ndim, VERSION_MACRO, PUT );
        length = 1; db( NODE_MACRO_GENERATE, inod, &icontrol, ddum, length, VERSION_MACRO, PUT );
        if ( nuknwn>0 ) {
          array_set( node_dof, 0., nuknwn );
          db( NODE_DOF, inod, idum, node_dof, nuknwn, VERSION_MACRO, PUT );
          db( NODE_DOF_START_REFINED, inod, idum, node_dof, nuknwn, VERSION_MACRO, PUT );
        }
        if ( control_mesh_macro_set_node_boundary==-YES && inod!=(nnod_macro-1)/2 ) {
          length = 1; db( NODE_BOUNDARY, inod, &node_boundary, ddum, 
            length, VERSION_MACRO, PUT );
        }
      }
      for ( element=0; element<nel_macro; element++ ) {
        if ( ndim==3 ) {
          el[0] = -HEX8; nnol = 8;
          geometry_macro_name = -GEOMETRY_SPHERE;
        }
        else {
          assert( ndim==2 );
          el[0] = -QUAD4; nnol = 4;
          geometry_macro_name = -GEOMETRY_CIRCLE;
        }
        for ( inol=0; inol<nnol; inol++ ) el[1+inol] = nodes_macro[inol][element];
        length = 1 + nnol;
        db( ELEMENT, element, el, ddum, length, VERSION_MACRO, PUT );
        length = 1;
        db( ELEMENT_GROUP, element, &element_group, ddum, length, VERSION_MACRO, PUT );
        length = 1;
        db( ELEMENT_MACRO_GENERATE, element, &icontrol, ddum, length, VERSION_MACRO, PUT );
      }

      array_move( control_mesh_macro_parameters, geometry_macro, ndim+1 );
      geometry_macro[ndim+1] = EPS_COORD;
      db_max_index( geometry_macro_name, max_geometry, VERSION_NORMAL, GET );
      max_geometry++; length = ndim + 2;
      db( geometry_macro_name, max_geometry, idum, geometry_macro, length, VERSION_NORMAL, PUT );
      use_control_refine_globally_geometry = 1;
      control_refine_globally_geometry[0] = geometry_macro_name;
      control_refine_globally_geometry[1] = max_geometry;
      for ( i=0; i<n; i++ ) {
        length_control_refine_globally = 1;
        control_refine_globally[0] = -H_REFINEMENT;
        refine_globally( control_refine_globally, length_control_refine_globally, 
          use_control_refine_globally_geometry, control_refine_globally_geometry, 
          PROJECT_EXACT, VERSION_MACRO );
      }
      db_delete_index( geometry_macro_name, max_geometry, VERSION_NORMAL );

    }
    else if ( macro_item==-BAR || macro_item==-RECTANGLE || macro_item==-BRICK ) {
      if      ( macro_item==-BAR ) {
        nx = 1;
        ny = 1;
        nz = control_mesh_macro[2];
        if ( nz<2 ) db_error( CONTROL_MESH_MACRO, icontrol );
        lenz = control_mesh_macro_parameters[ndim];
      }
      else if ( macro_item==-RECTANGLE ) {
        nx = 1;
        ny = control_mesh_macro[2];
        if ( ny<2 ) db_error( CONTROL_MESH_MACRO, icontrol );
        nz = control_mesh_macro[3];
        if ( nz<2 ) db_error( CONTROL_MESH_MACRO, icontrol );
        leny = control_mesh_macro_parameters[ndim];
        lenz = control_mesh_macro_parameters[ndim+1];
      }
      else {
        assert( macro_item==-BRICK );
        nx = control_mesh_macro[2];
        if ( nx<2 ) db_error( CONTROL_MESH_MACRO, icontrol );
        ny = control_mesh_macro[3];
        if ( ny<2 ) db_error( CONTROL_MESH_MACRO, icontrol );
        nz = control_mesh_macro[4];
        if ( nz<2 ) db_error( CONTROL_MESH_MACRO, icontrol );
        lenx = control_mesh_macro_parameters[ndim];
        leny = control_mesh_macro_parameters[ndim+1];
        lenz = control_mesh_macro_parameters[ndim+2];
      }
      for ( ix=0; ix<nx; ix++ ) {
        if ( ndim>=3 ) dx = ix*lenx/(nx-1);
        for ( iy=0; iy<ny; iy++ ) {
          if ( ndim>=2 ) dy = iy*leny/(ny-1);
          for ( iz=0; iz<nz; iz++ ) {
            if ( ndim>=1 ) dz = iz*lenz/(nz-1);
            if      ( ndim==1 ) {
              coord[0] = -0.5*lenz + dz;
            }
            else if ( ndim==2 ) {
              coord[0] = -0.5*leny + dy;
              coord[1] = -0.5*lenz + dz;
            }
            else {
              assert( ndim==3 );
              coord[0] = -0.5*lenx + dx;
              coord[1] = -0.5*leny + dy;
              coord[2] = -0.5*lenz + dz;
            }
            array_add( coord, control_mesh_macro_parameters, coord, ndim );
            db( NODE, inod, idum, coord, ndim, VERSION_MACRO, PUT );
            db( NODE_START_REFINED, inod, idum, coord, ndim, VERSION_MACRO, PUT );
            length = 1; db( NODE_MACRO_GENERATE, inod, &icontrol, ddum, length, VERSION_MACRO, PUT );
            if ( nuknwn>0 ) {
              array_set( node_dof, 0., nuknwn );
              db( NODE_DOF, inod, idum, node_dof, nuknwn, VERSION_MACRO, PUT );
              db( NODE_DOF_START_REFINED, inod, idum, node_dof, nuknwn, VERSION_MACRO, PUT );
            }
            test1 = ndim>=3 && (ix==0 || ix==nx-1);
            test2 = ndim>=2 && (iy==0 || iy==ny-1);
            test3 = ndim>=1 && (iz==0 || iz==nz-1);
            if ( test1 || test2 || test3 ) {
              if (  control_mesh_macro_set_node_boundary==-YES ) {
                length = 1; db( NODE_BOUNDARY, inod, &node_boundary, ddum, 
                  length, VERSION_MACRO, PUT );
              }
            }
            array_set( node_remesh_allowed, -YES, ndim );
            if ( test1 ) node_remesh_allowed[0] = -NO;
            if ( test2 ) node_remesh_allowed[1] = -NO;
            if ( test3 ) node_remesh_allowed[2] = -NO;
            length = ndim; db( NODE_REMESH_ALLOWED, inod, node_remesh_allowed, ddum, 
               length, VERSION_MACRO, PUT );
            inod++;
          }
        }
      }
      for ( ix=0; ix<scalar_imax(1,nx-1); ix++ ) {
        for ( iy=0; iy<scalar_imax(1,ny-1); iy++ ) {
          for ( iz=0; iz<scalar_imax(1,nz-1); iz++ ) {
            if      ( ndim==1 ) {
              el[0] = -BAR2;
              nnol = 2;
              nodes[0] = iz;
              nodes[1] = iz + 1;
            }
            else if ( ndim==2 ) {
              el[0] = -QUAD4;
              nnol = 4;
              nodes[0] = iy*nz+iz;
              nodes[1] = iy*nz+iz + 1;
              nodes[2] = (iy+1)*nz+iz;
              nodes[3] = (iy+1)*nz+iz + 1;
            }
            else {
              assert( ndim==3 );
              el[0] = -HEX8;
              nnol = 8;
              nodes[0] = ix*nz*ny+iy*nz+iz;
              nodes[1] = ix*nz*ny+iy*nz+iz + 1;
              nodes[2] = ix*nz*ny+(iy+1)*nz+iz;
              nodes[3] = ix*nz*ny+(iy+1)*nz+iz + 1;
              nodes[4] = (ix+1)*nz*ny+iy*nz+iz;
              nodes[5] = (ix+1)*nz*ny+iy*nz+iz + 1;
              nodes[6] = (ix+1)*nz*ny+(iy+1)*nz+iz;
              nodes[7] = (ix+1)*nz*ny+(iy+1)*nz+iz + 1;
            }
            array_move( nodes, &el[1], nnol ); 
            length = 1 + nnol;
            db( ELEMENT, element, el, ddum, length, VERSION_MACRO, PUT );
            length = 1;
            db( ELEMENT_GROUP, element, &element_group, ddum, length, VERSION_MACRO, PUT );
            length = 1;
            db( ELEMENT_MACRO_GENERATE, element, &icontrol, ddum, length, VERSION_MACRO, PUT );
            element++;
          }
        }
      }
    }
    else
      db_error( CONTROL_MESH_MACRO, icontrol );

      // split or p_refine mesh
    if ( db( CONTROL_MESH_MACRO_ELEMENT, icontrol, &control_mesh_macro_element, 
        ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS ) ) {
      if      ( control_mesh_macro_element==-TRIA3 )
        mesh_split( VERSION_MACRO );
      else if ( control_mesh_macro_element==-TET4 )
        mesh_split( VERSION_MACRO );
      else if ( control_mesh_macro_element==-BAR3  || 
                control_mesh_macro_element==-QUAD9 ||
                control_mesh_macro_element==-HEX27 ) {
        use_control_refine_globally_geometry = 0;
        length_control_refine_globally = 1;
        control_refine_globally[0] = -P_REFINEMENT;
        refine_globally( control_refine_globally, length_control_refine_globally, 
          use_control_refine_globally_geometry, control_refine_globally_geometry, 
          PROJECT_EXACT, VERSION_MACRO );
      }
      else if ( control_mesh_macro_element==-TRIA6  || 
                control_mesh_macro_element==-TET10 ) {
        use_control_refine_globally_geometry = 0;
        length_control_refine_globally = 1;
        control_refine_globally[0] = -P_REFINEMENT;
        refine_globally( control_refine_globally, length_control_refine_globally, 
          use_control_refine_globally_geometry, control_refine_globally_geometry, 
          PROJECT_EXACT, VERSION_MACRO );
        mesh_split( VERSION_MACRO );
      }
    }

      // renumber
    renumbering( VERSION_MACRO, NO, 1, 1, idum, idum );

      // add to existing mesh
    mesh_add( VERSION_MACRO, VERSION_NORMAL );

    if ( swit ) pri( "Out routine MACRO" );

  }
  
}

