 /*
    Copyright (C) 2000  Dennis Roddeman
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

// osman buyukisik  osman@fuse.net
// crack stress intensity calcs using displacement correlation method.
// works for 2D only
// only needs "crack_nodes" and "crack_elementgroup" commands to specify
// and "control_crack -calculate_stressintensityfactor " to calculate SIF
// print using "control_print  -crack_stressintensityfactor "
// due to problems, quarter point singular displacement correlation method is
// not working with tochnog. Use corner points only. i.e.
// crack_nodes  tip  cn1 cn1 cn2 cn2
//                            --*-
//                       --     ^-- corner node top
//                  --
//              -
//    --------   <---- tip
//              -
//                   --
//                        --    v-- corner node bottom
//                            --*-  
// Not as efficient as the singular point method, but at least it works :-)
// (needs a lot more elements  to get the same accuracy)

#include "tochnog.h"

#define MAX_END_NODE 5

void crack( void )

{
  long int k=0, idim=0, symm_crack=0, end_node=MAX_END_NODE,
    icontrol=0, length=0, ldum=0, membrane=-NO, crack_elementgroup=0, 
    swit=0, inode=0, control_crack=0, idum[1], crack_nodes[MAX_END_NODE];
  double vfactor=1., dummy=0., young=0.,g=0., poisson=0., kappa=0., factor=1., 
    coeff1=0.,coeff2=0., crack_length=0., ddum[1], crack_stressintensityfactor[2], 
    node_dof[MUKNWN], coords1[MDIM], coords2[MDIM], crack_direction[MDIM], 
    crack_disp[MAX_END_NODE][MUKNWN];

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  db( CONTROL_CRACK, icontrol, &control_crack, ddum, ldum, 
      VERSION_NORMAL, GET_IF_EXISTS);

  if ( control_crack==-CALCULATE_STRESSINTENSITYFACTOR ) {

    swit = set_swit(-1,-1,"crack");
    if ( swit ) pri( "In routine CRACK" );

    array_set(crack_direction, 0., ndim );
    for( k=0; k<MAX_END_NODE; k++ ){ 
      crack_disp[k][0]=0.0;
      crack_disp[k][1]=0.0;
    }

      //
      // crack_nodes are set of nodes defining the crack_face
      // 1st is the crack tip
      // 2nd lower face 1/4 node
      // 3rd lower face corner node
      // 4th upper face 1/4 node    :  if it exists otherwise 0
      // 5th upper face corner node :  ditto
      // for a symmetric model, one of the faces does not exist.
      // just put 0 for nodes 4,5.
      // crack_length is the distance from CTIP to one of the corner nodes.
      // 
      //
      // get displacements (velocities ) at the crack nodes
    db( CRACK_NODES, 0, crack_nodes, ddum, ldum, VERSION_NORMAL, GET );
    length = MUKNWN; end_node = MAX_END_NODE; symm_crack=0;
    if ( crack_nodes[3]==0 ) { 
      symm_crack = 1; 
      end_node = 3; 
    }
    for ( inode=0; inode<end_node; inode++) {
      db( NODE_DOF, crack_nodes[inode], idum, node_dof, length, VERSION_NORMAL, GET );
      for( idim=0; idim<ndim; idim++ ) crack_disp[inode][idim] = node_dof[vel_indx+idim];
    }


    if (symm_crack==0 ) {   
        //  non symmetric crack model nodes 4,5 exists.
        // get the nodal coords of CT and  2 CORNER nodes
        // define a vector from CT to avg of 2 CORNER NODES
        // normalize it. this is the direction vector
      db( NODE, crack_nodes[2], idum, coords1, ldum, VERSION_NORMAL, GET );
      db( NODE, crack_nodes[4], idum, coords2, ldum, VERSION_NORMAL, GET );
      array_add( coords1 , coords2, coords1 , ndim );
      array_multiply(coords1, coords1 , 0.5, ndim );
      db( NODE, crack_nodes[0], idum, coords2, ldum, VERSION_NORMAL, GET );
    }
    else {   
        // symmetric crack model. nodes 4,5 are 0.
        // direction is simply line from CT to CORNER node
      db( NODE, crack_nodes[2], idum, coords1, ldum, VERSION_NORMAL, GET );
      db( NODE, crack_nodes[0], idum, coords2, ldum, VERSION_NORMAL, GET );
    }

    array_subtract(coords2, coords1, crack_direction, ndim);
    dummy = 0.;
    for( idim=0; idim<ndim; idim++ ) dummy += crack_direction[idim]*crack_direction[idim];
    crack_length=sqrt(dummy);
    dummy = 1.0/crack_length;
    array_multiply( crack_direction, crack_direction, dummy, ndim );

      // now rotate the displacements into the crack coord system
      // using crack_direction vector.
    for (k=0; k<end_node; k++ ) {
      dummy = crack_direction[0]*crack_disp[k][0] + 
        crack_direction[1]*crack_disp[k][1];
      crack_disp[k][1]= -crack_direction[1]*crack_disp[k][0] +
        crack_direction[0]*crack_disp[k][1];
      crack_disp[k][0]=dummy;
    }

     // get properties for the crack region
    db( CRACK_ELEMENTGROUP,0,&crack_elementgroup,ddum,ldum,
      VERSION_NORMAL, GET_IF_EXISTS );
    db( GROUP_VOLUME_FACTOR, crack_elementgroup, idum, &vfactor,
      ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( GROUP_MATERI_ELASTI_YOUNG, crack_elementgroup,
      idum, &young, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( GROUP_MATERI_ELASTI_POISSON, crack_elementgroup,
      idum, &poisson, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( GROUP_MATERI_MEMBRANE, crack_elementgroup, 
      &membrane, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );   
    if ( membrane==-YES ) {
      kappa=(3.0-poisson)/(1.0+poisson);
    }
    else {     
        // plane strain or axisymmetric . 
      kappa=3.0-(4.0*poisson);
    }
    g = young*0.5/(1.0+poisson);
    factor=(g/(kappa+1.0))*sqrt((2.0*PIRAD)/crack_length);
    //factor=1.27*factor;  // tochnog only correction, see below
      // for singular T6 elements (midside nodes at 1/4 )
      // coeff1 = 4.0
      // coeff2 = -1.0
      // for regular (non-singular elements)
      // coeff1 = 0.0
      // coeff2 = 1.0
    coeff1 = 4.0; // this is the only option that works due to 
    coeff2 = -1.0; // differences in displacements (velocities)
                  // between tochnog and other FE programs.
    if (symm_crack==0 ) {
      crack_stressintensityfactor[0] = factor*(
        (coeff1*(crack_disp[3][1]-crack_disp[1][1]))+
        (coeff2*(crack_disp[4][1]-crack_disp[2][1])));
      crack_stressintensityfactor[1] = factor*(
        (coeff1*(crack_disp[3][0]-crack_disp[1][0]))+
        (coeff2*(crack_disp[4][0]-crack_disp[2][0])));
    }
    else {
      crack_stressintensityfactor[0] = 2.0*factor*(
        (coeff1*crack_disp[1][1])+ coeff2*crack_disp[2][1] );
      crack_stressintensityfactor[1] = 2.0*factor*(
        (coeff1*crack_disp[1][1])+ coeff2*crack_disp[2][0] );
    }

/*  these restrict the output to positive numbers. 
    if (crack_stressintensityfactor[0]<0.0)
      crack_stressintensityfactor[0]=-crack_stressintensityfactor[0];
    if (crack_stressintensityfactor[1]<0.0)
      crack_stressintensityfactor[1]=-crack_stressintensityfactor[1];
*/
    length = 2;
    db( CRACK_STRESSINTENSITYFACTOR, 0, idum, crack_stressintensityfactor, 
        length, VERSION_NORMAL, PUT );

    if ( swit ) {
      pri( "crack length", crack_length );
      pri( "kappa", kappa );
      pri( "factor", factor );
      pri( "crack_stressintensityfactor1", crack_stressintensityfactor[0] );
      pri( "crack_stressintensityfactor2", crack_stressintensityfactor[1] );
    }

    if ( swit ) pri( "Out routine CRACK" );
  }
}
