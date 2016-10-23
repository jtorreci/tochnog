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

#define MAX_SIDE 12

void tendon_distribute( )

{
 long int it=0, element=0, max_element=0, itendon=0, max_tendon=0,
   ldum=0, len=0, nnol=0, iel=0, inol=0, jnol=0, knol=0, 
   length_tendon_split_element=0, indx=0, tmp_indx=0,
   length_element_tendon_number=0, 
   length_element_tendon_strain=0,
   length_element_tendon_stress=0,
   tmp_length_element_tendon_strain=0,
   tmp_length_element_tendon_stress=0,
   length_element_tendon_volume=0, 
   length_element_tendon_intersections=0, 
   length_element_tendon_direction=0,
   inod=0, jnod=0, knod=0, name=0, iside=0, nside=0, nnol_side=0, 
   nintersect=0, swit=0, swit_element=0, already_used=0,
   itendon_split=0, max_tendon_split=0, test_if_new_intersection=1, 
   tmp_length_element_tendon_number=0,
   yes_it_is_a_new_intersection=0, idum[1], 
   sides[MAX_SIDE*MNOL], el[MNOL+1], nodes[MNOL], 
   element_tendon_number[DATA_ITEM_SIZE], 
   tmp_element_tendon_number[DATA_ITEM_SIZE], 
   tendon_split_element[DATA_ITEM_SIZE];
 double length=0., volume=0., tendon_area=0., tendon_length=0., 
   iso_tendon=0., tendon_elasti=0.,
   tendon_stress=0., iso_2d=0.,
   ddum[1], iso_3d[2], work[MDIM], tendon_coord[MDIM], 
   nodei[MDIM], nodej[MDIM], nodek[MDIM], 
   tendon_vector[MDIM], tmp_vector[MDIM], tendon_split[DATA_ITEM_SIZE],
   intersect_coord[2*MDIM], intersect_iso_tendon[2], 
   tendon[DATA_ITEM_SIZE], 
   element_tendon_volume[DATA_ITEM_SIZE], 
   element_tendon_intersections[DATA_ITEM_SIZE], 
   element_tendon_strain[DATA_ITEM_SIZE],
   element_tendon_stress[DATA_ITEM_SIZE],
   tmp_element_tendon_strain[DATA_ITEM_SIZE],
   tmp_element_tendon_stress[DATA_ITEM_SIZE],
   element_tendon_direction[DATA_ITEM_SIZE];

 db_max_index( TENDON, max_tendon, VERSION_NORMAL, GET );
 if ( max_tendon>=0 ) {

   swit = set_swit(-1,-1,"tendon_distribute");
   if ( swit ) pri( "In routine TENDON_DISTRIBUTE" );

   db_version_copy( VERSION_NORMAL, VERSION_TMP );

   db_delete( ELEMENT_TENDON_DIRECTION, VERSION_NORMAL );
   db_delete( ELEMENT_TENDON_VOLUME, VERSION_NORMAL );
   db_delete( ELEMENT_TENDON_INTERSECTIONS, VERSION_NORMAL );
   db_delete( ELEMENT_TENDON_NUMBER, VERSION_NORMAL );
   db_delete( ELEMENT_TENDON_STRAIN, VERSION_NORMAL );
   db_delete( ELEMENT_TENDON_STRESS, VERSION_NORMAL );

   db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
   
   for ( itendon=0; itendon<=max_tendon; itendon++ ) {
     if ( db_active_index( TENDON, itendon, VERSION_NORMAL ) ) {
       db( TENDON, itendon, idum, tendon, ldum, VERSION_NORMAL, GET );
       tendon_area = tendon[2*ndim];
       array_move( tendon, tendon_coord, ndim );
       array_move( &tendon[ndim], tmp_vector, ndim ); 
       array_subtract( tmp_vector, tendon_coord, tendon_vector, ndim );
       tendon_length = array_size( tendon_vector, ndim );
       tendon_elasti = 0.; 
       tendon_stress = 0.; 
       db( TENDON_ELASTI, itendon, idum, 
         &tendon_elasti, ldum, VERSION_NORMAL, GET_IF_EXISTS );
       db( TENDON_STRESS, itendon, idum, 
         &tendon_stress, ldum, VERSION_NORMAL, GET_IF_EXISTS );
       if ( tendon_length==0. ) db_error( TENDON, itendon );
       if ( swit ) {
         pri( "tendon_coord", tendon_coord, ndim );
         pri( "tendon_vector", tendon_vector, ndim );
         pri( "tendon_area", tendon_area );
       }
       for ( element=0; element<=max_element; element++ ) {
         if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) ) {
           swit_element = swit;
           swit = swit && set_swit(element,-1,"tendon_distribute");
           if ( swit ) pri( "element", element );
           db( ELEMENT, element, el, ddum, len, VERSION_NORMAL, GET );
           name = el[0]; nnol = len - 1; array_move( &el[1], nodes, nnol );
           if      ( name==-BAR2 ) {
             nside = 2; nnol_side = 1;
             sides[0*nnol_side+0] = 0; sides[0*nnol_side+1] = 1;
           }
           else if ( name==-BAR3 ) {
             nside = 2; nnol_side = 1;
             sides[0*nnol_side+0] = 0; sides[0*nnol_side+1] = 2;
           }
           else if ( name==-TRIA3 ) {
             nside = 3; nnol_side = 2;
             sides[0*nnol_side+0] = 0; sides[0*nnol_side+1] = 1;
             sides[1*nnol_side+0] = 1; sides[1*nnol_side+1] = 2;
             sides[2*nnol_side+0] = 2; sides[2*nnol_side+1] = 0;
           }
           else if ( name==-QUAD4 ) {
             nside = 4; nnol_side = 2;
             sides[0*nnol_side+0] = 0; sides[0*nnol_side+1] = 1;
             sides[1*nnol_side+0] = 1; sides[1*nnol_side+1] = 3;
             sides[2*nnol_side+0] = 3; sides[2*nnol_side+1] = 2;
             sides[3*nnol_side+0] = 2; sides[3*nnol_side+1] = 0;
           }
           else if ( name==-QUAD9 ) {
             nside = 4; nnol_side = 2;
             sides[0*nnol_side+0] = 0; sides[0*nnol_side+1] = 2;
             sides[1*nnol_side+0] = 2; sides[1*nnol_side+1] = 8;
             sides[2*nnol_side+0] = 8; sides[2*nnol_side+1] = 6;
             sides[3*nnol_side+0] = 6; sides[3*nnol_side+1] = 0;
           }
           else if ( name==-TET4 ) {
             nside = 4; nnol_side = 3;
             sides[0*nnol_side+0] = 0; sides[0*nnol_side+1] = 1;
             sides[0*nnol_side+2] = 2;
             sides[1*nnol_side+0] = 1; sides[1*nnol_side+1] = 2;
             sides[1*nnol_side+2] = 3;
             sides[2*nnol_side+0] = 0; sides[2*nnol_side+1] = 2;
             sides[2*nnol_side+2] = 3;
             sides[3*nnol_side+0] = 0; sides[3*nnol_side+1] = 1;
             sides[3*nnol_side+2] = 3;
           }
           else if ( name==-HEX8 ) {
             nside = 12; nnol_side = 3;
             sides[0*nnol_side+0] = 0; sides[0*nnol_side+1] = 1;
             sides[0*nnol_side+2] = 2; 
             sides[1*nnol_side+0] = 1; sides[1*nnol_side+1] = 3;
             sides[1*nnol_side+2] = 2; 
             sides[2*nnol_side+0] = 4; sides[2*nnol_side+1] = 5;
             sides[2*nnol_side+2] = 6; 
             sides[3*nnol_side+0] = 5; sides[3*nnol_side+1] = 7;
             sides[3*nnol_side+2] = 6; 
             sides[4*nnol_side+0] = 0; sides[4*nnol_side+1] = 1;
             sides[4*nnol_side+2] = 4; 
             sides[5*nnol_side+0] = 1; sides[5*nnol_side+1] = 5;
             sides[5*nnol_side+2] = 4; 
             sides[6*nnol_side+0] = 2; sides[6*nnol_side+1] = 3;
             sides[6*nnol_side+2] = 6; 
             sides[7*nnol_side+0] = 3; sides[7*nnol_side+1] = 6;
             sides[7*nnol_side+2] = 7; 
             sides[8*nnol_side+0] = 1; sides[8*nnol_side+1] = 3;
             sides[8*nnol_side+2] = 5; 
             sides[9*nnol_side+0] = 3; sides[9*nnol_side+1] = 7;
             sides[9*nnol_side+2] = 5; 
             sides[10*nnol_side+0] = 0; sides[10*nnol_side+1] = 2;
             sides[10*nnol_side+2] = 4; 
             sides[11*nnol_side+0] = 2; sides[11*nnol_side+1] = 6;
             sides[11*nnol_side+2] = 4; 
             if ( swit ) pri( "sides", sides, nside, nnol_side );
           }
           else if ( name==-HEX27 ) {
             nside = 12; nnol_side = 3;
             sides[0*nnol_side+0] = 0; sides[0*nnol_side+1] = 2;
             sides[0*nnol_side+2] = 6; 
             sides[1*nnol_side+0] = 2; sides[1*nnol_side+1] = 8;
             sides[1*nnol_side+2] = 6; 
             sides[2*nnol_side+0] = 18; sides[2*nnol_side+1] = 20;
             sides[2*nnol_side+2] = 24; 
             sides[3*nnol_side+0] = 20; sides[3*nnol_side+1] = 26;
             sides[3*nnol_side+2] = 24; 
             sides[4*nnol_side+0] = 0; sides[4*nnol_side+1] = 2;
             sides[4*nnol_side+2] = 18; 
             sides[5*nnol_side+0] = 2; sides[5*nnol_side+1] = 20;
             sides[5*nnol_side+2] = 18; 
             sides[6*nnol_side+0] = 6; sides[6*nnol_side+1] = 8;
             sides[6*nnol_side+2] = 24; 
             sides[7*nnol_side+0] = 8; sides[7*nnol_side+1] = 24;
             sides[7*nnol_side+2] = 26; 
             sides[8*nnol_side+0] = 2; sides[8*nnol_side+1] = 8;
             sides[8*nnol_side+2] = 20; 
             sides[9*nnol_side+0] = 8; sides[9*nnol_side+1] = 26;
             sides[9*nnol_side+2] = 20; 
             sides[10*nnol_side+0] = 0; sides[10*nnol_side+1] = 6;
             sides[10*nnol_side+2] = 18; 
             sides[11*nnol_side+0] = 6; sides[11*nnol_side+1] = 24;
             sides[11*nnol_side+2] = 18; 
             if ( swit ) pri( "sides", sides, nside, nnol_side );
           }
           else {
             pri( "Tendons not available for ", name );
             exit(TN_EXIT_STATUS);
           }
           if ( nside ) {
             nintersect = 0; 
             for ( iside=0; iside<nside && nintersect<2; iside++ ) {
               if ( swit ) pri( "iside", iside );
               test_if_new_intersection = 0;
               if ( ndim==1 ) {
                 inod = nodes[sides[iside*nnol_side+0]];
                 db( NODE, inod, idum, nodei, ldum, VERSION_NORMAL, GET );
                 array_subtract( nodei, tendon_coord, tmp_vector, ndim );
                 iso_tendon = array_inproduct( tmp_vector, tendon_vector, 
                   ndim ) / ( tendon_length*tendon_length );
                 if ( iso_tendon>=0. && iso_tendon<=1. ) {
                   array_move( nodei, tmp_vector, ndim );
                   test_if_new_intersection = 1;
                 }
               }
               else if ( ndim==2 ) {
                 inol = sides[iside*nnol_side+0]; inod = nodes[inol];
                 jnol = sides[iside*nnol_side+1]; jnod = nodes[jnol];
                 db( NODE, inod, idum, nodei, ldum, VERSION_NORMAL, GET );
                 db( NODE, jnod, idum, nodej, ldum, VERSION_NORMAL, GET );
                 if ( intersect_line_with_line( &tendon[0], &tendon[ndim], 
                   nodei, nodej, iso_tendon, iso_2d ) ) {
                   if ( iso_tendon>=-EPS_ISO && iso_tendon<=1.+EPS_ISO &&
                        iso_2d>=-EPS_ISO && iso_2d<=1.+EPS_ISO ) {
                     array_multiply( tendon_vector, tmp_vector, 
                       iso_tendon, ndim );
                     array_add( tendon_coord, tmp_vector, tmp_vector, ndim );
                     test_if_new_intersection = 1;
                   }
                 }
               }
               else {
                 assert( ndim==3 );
                 inol = sides[iside*nnol_side+0]; inod = nodes[inol];
                 jnol = sides[iside*nnol_side+1]; jnod = nodes[jnol];
                 knol = sides[iside*nnol_side+2]; knod = nodes[knol];
                 db( NODE, inod, idum, nodei, ldum, VERSION_NORMAL, GET );
                 db( NODE, jnod, idum, nodej, ldum, VERSION_NORMAL, GET );
                 db( NODE, knod, idum, nodek, ldum, VERSION_NORMAL, GET );
                 if ( intersect_line_with_triangle( &tendon[0], &tendon[ndim], 
                   nodei, nodej, nodek, iso_tendon, iso_3d ) ) {
                   if ( iso_tendon>=-EPS_ISO && iso_tendon<=1.+EPS_ISO &&
                        iso_3d[0]>=-EPS_ISO && iso_3d[0]<=1.+EPS_ISO &&
                        iso_3d[1]>=-EPS_ISO && iso_3d[1]<=1.+EPS_ISO ) {
                     array_multiply( tendon_vector, tmp_vector, 
                       iso_tendon, ndim );
                     array_add( tendon_coord, tmp_vector, tmp_vector, ndim );
                     test_if_new_intersection = 1;
                   }
                 }
               }
               if ( test_if_new_intersection ) {
                 if ( nintersect==0 ) yes_it_is_a_new_intersection = 1;
                 else {
                   if ( array_distance(tmp_vector,
                        intersect_coord,work,ndim)>1.e-8 )
                     yes_it_is_a_new_intersection = 1;
                   else
                     yes_it_is_a_new_intersection = 0;
                 }
               }
               else yes_it_is_a_new_intersection = 0;
               if ( yes_it_is_a_new_intersection ) {
                 array_move( tmp_vector, &intersect_coord[nintersect*ndim], 
                   ndim );
                 intersect_iso_tendon[nintersect] = iso_tendon;
                 nintersect++;
                 if ( swit ) {
                   pri( "intersect_coord", intersect_coord, nintersect, 
                     ndim );
                   pri( "intersect_iso_tendon", intersect_iso_tendon, 
                     nintersect );
                 }
               }
             }
             if ( nintersect==2 ) {
                 // volume
               array_subtract( &intersect_coord[1*ndim], 
                 &intersect_coord[0*ndim], tmp_vector, ndim );
               length = array_size( tmp_vector, ndim );
               volume = length * tendon_area;
                 // test if this tendon part is not already used
               db_max_index( TENDON_SPLIT, max_tendon_split, VERSION_NORMAL, GET );
               iso_tendon =  ( intersect_iso_tendon[0] + 
                 intersect_iso_tendon[1] ) / 2.;
               already_used = itendon_split = 0;
               while ( itendon_split<=max_tendon_split && !already_used ) {
                 if ( db_active_index( TENDON_SPLIT, itendon_split, VERSION_NORMAL ) ) {
                   db( TENDON_SPLIT, itendon_split, idum, tendon_split, 
                     len, VERSION_NORMAL, GET );
                   if ( iso_tendon>tendon_split[0] && 
                     iso_tendon<tendon_split[1] ) already_used = 1;
                   else
                     itendon_split++;
                 }
                 else
                   itendon_split++;
               }
               if ( swit ) pri( "already_used", already_used );
               if ( already_used )
                 db( TENDON_SPLIT_ELEMENT, itendon_split, tendon_split_element, 
                   ddum, length_tendon_split_element, VERSION_NORMAL, GET );
               else {
                 if ( intersect_iso_tendon[0]<intersect_iso_tendon[1] ) {
                   tendon_split[0] = intersect_iso_tendon[0];
                   tendon_split[1] = intersect_iso_tendon[1];
                 }
                 else {
                   tendon_split[0] = intersect_iso_tendon[1];
                   tendon_split[1] = intersect_iso_tendon[0];
                 }
                 itendon_split = max_tendon_split + 1; len = 2; 
                 db( TENDON_SPLIT, itendon_split, idum, tendon_split, len, 
                   VERSION_NORMAL, PUT );
                 length_tendon_split_element = 0;
               }
               tendon_split_element[length_tendon_split_element] = element; 
               length_tendon_split_element++;
               if ( swit ) {
                 pri( "itendon_split", itendon_split );
                 pri( "tendon_split_element", tendon_split_element,
                   length_tendon_split_element );
               }
               db( TENDON_SPLIT_ELEMENT, itendon_split, tendon_split_element,
                 ddum, length_tendon_split_element, VERSION_NORMAL, PUT );
                 // (re)distribute tendon over elements
               for ( it=0; it<length_tendon_split_element; it++ ) {
                 iel = tendon_split_element[it];
                 if ( swit ) pri( "iel", iel );
                 if ( db_active_index( ELEMENT_TENDON_NUMBER, iel, VERSION_NORMAL ) ) {
                   db( ELEMENT_TENDON_NUMBER, iel, element_tendon_number, 
                     ddum, length_element_tendon_number, VERSION_NORMAL, GET );
                   db( ELEMENT_TENDON_STRAIN, iel, idum, 
                     element_tendon_strain, length_element_tendon_strain, 
                     VERSION_NORMAL, GET );
                   db( ELEMENT_TENDON_STRESS, iel, idum, 
                     element_tendon_stress, length_element_tendon_stress, 
                     VERSION_NORMAL, GET );
                   db( ELEMENT_TENDON_VOLUME, iel, idum, 
                     element_tendon_volume, length_element_tendon_volume, 
                     VERSION_NORMAL, GET );
                   db( ELEMENT_TENDON_INTERSECTIONS, iel, idum, 
                     element_tendon_intersections, 
                     length_element_tendon_intersections, 
                     VERSION_NORMAL, GET );
                   db( ELEMENT_TENDON_DIRECTION, iel, idum, 
                     element_tendon_direction, 
                     length_element_tendon_direction, VERSION_NORMAL, GET );
                   array_member(element_tendon_number,itendon,
                     length_element_tendon_number,indx);
                 }
                 else {
                   length_element_tendon_number = 0;
                   length_element_tendon_strain = 0;
                   length_element_tendon_stress = 0;
                   length_element_tendon_volume = 0;
                   length_element_tendon_intersections = 0;
                   length_element_tendon_direction = 0;
                   indx = -1;
                 }
                 if ( indx<0 ) {
                   indx = length_element_tendon_number;
                   length_element_tendon_number++; 
                   length_element_tendon_strain++; 
                   length_element_tendon_stress++; 
                   length_element_tendon_volume++; 
                   length_element_tendon_intersections += 2*ndim; 
                   length_element_tendon_direction += MDIM; 
                   if ( length_element_tendon_number>
                     db_data_length(ELEMENT_TENDON_NUMBER) ) {
                     pri( "MTENDON needs to become ", length_element_tendon_number );
                     pri( "Change it in tochnog.h and recompile.");
                     exit(TN_EXIT_STATUS);
                   }
                 }
                 if ( length_element_tendon_number>MTENDON ) {
                   cout << "Error: to many tendons in element " << iel << ".\n";
                   exit(TN_EXIT_STATUS);
                 }
                 element_tendon_number[indx] = itendon;
                 if ( db_active_index( ELEMENT_TENDON_NUMBER, iel, VERSION_TMP ) ) {
                     // use old element tendon stress if present for this tendon number
                   db( ELEMENT_TENDON_NUMBER, iel, tmp_element_tendon_number, 
                     ddum, tmp_length_element_tendon_number, VERSION_TMP, GET );
                   array_member(tmp_element_tendon_number,itendon,
                     tmp_length_element_tendon_number,tmp_indx);
                   if ( tmp_indx>=0 ) {
                     db( ELEMENT_TENDON_STRAIN, iel, idum, tmp_element_tendon_strain, 
                       tmp_length_element_tendon_strain, VERSION_TMP, GET );
                     element_tendon_strain[indx] = tmp_element_tendon_strain[tmp_indx];
                     db( ELEMENT_TENDON_STRESS, iel, idum, tmp_element_tendon_stress, 
                       tmp_length_element_tendon_stress, VERSION_TMP, GET );
                     element_tendon_stress[indx] = tmp_element_tendon_stress[tmp_indx];
                   }
                   else {
                     element_tendon_stress[indx] = tendon_stress;
                     if ( tendon_elasti==0. ) 
                       element_tendon_strain[indx] = 0.;
                     else
                       element_tendon_strain[indx] = tendon_stress / tendon_elasti;
                   }
                 }
                 else {
                   element_tendon_stress[indx] = tendon_stress;
                   if ( tendon_elasti==0. ) 
                     element_tendon_strain[indx] = 0.;
                   else
                     element_tendon_strain[indx] = tendon_stress / tendon_elasti;
                 }
                 element_tendon_volume[indx] = 
                   volume/length_tendon_split_element;
                 array_move( intersect_coord, 
                   &element_tendon_intersections[indx*2*ndim], 2*ndim );
                 array_set( &element_tendon_direction[indx*MDIM], 0., MDIM );
                 array_move( tendon_vector, 
                   &element_tendon_direction[indx*MDIM], ndim );
                 array_normalize( &element_tendon_direction[indx*MDIM], 
                   MDIM );
                 db( ELEMENT_TENDON_NUMBER, iel, element_tendon_number, 
                   ddum, length_element_tendon_number, VERSION_NORMAL, PUT );
                 db( ELEMENT_TENDON_STRAIN, iel, idum, element_tendon_strain, 
                   length_element_tendon_strain, VERSION_NORMAL, PUT );
                 db( ELEMENT_TENDON_STRESS, iel, idum, element_tendon_stress, 
                   length_element_tendon_stress, VERSION_NORMAL, PUT );
                 db( ELEMENT_TENDON_VOLUME, iel, idum, element_tendon_volume, 
                   length_element_tendon_volume, VERSION_NORMAL, PUT );
                 db( ELEMENT_TENDON_INTERSECTIONS, iel, idum,
                   element_tendon_intersections,
                   length_element_tendon_intersections, VERSION_NORMAL, PUT );
                 db( ELEMENT_TENDON_DIRECTION, iel, idum, 
                   element_tendon_direction, 
                   length_element_tendon_direction, VERSION_NORMAL, PUT );
                 if ( swit ) {
                   pri( "element_tendon_number", element_tendon_number, 
                     length_element_tendon_number );
                   pri( "element_tendon_strain", element_tendon_strain, 
                     length_element_tendon_strain );
                   pri( "element_tendon_stress", element_tendon_stress, 
                     length_element_tendon_stress );
                   pri( "element_tendon_volume", element_tendon_volume, 
                     length_element_tendon_volume );
                   pri( "element_tendon_intersections",
                     element_tendon_intersections,
                     length_element_tendon_intersections );
                   pri( "element_tendon_direction", element_tendon_direction, 
                     length_element_tendon_direction );
                 }
               }
             }
           }
           swit = swit_element;
         }
       }
       db_delete(TENDON_SPLIT, VERSION_NORMAL);
       db_delete(TENDON_SPLIT_ELEMENT, VERSION_NORMAL);
     }
   }
   db_version_delete( VERSION_TMP );

   if ( swit ) pri( "Out routine TENDON_DISTRIBUTE" );
 }

}

void tendons( long int element, long int gr, long int nnol, 
  long int npoint, double volume, double new_d[], double old_unknowns[], 
  double new_unknowns[], double new_rot[], double inc_ept[], 
  double tendon_element_rhside[], double ddsdde_tendon[] )

{
  long int i=0, j=0, swit=0, idim=0, jdim=0, kdim=0, ldim=0,
    inol=0, iten=0, itendon=0, ntendon=0, ind_ddsdde=0, ipuknwn=0, 
    indx=0, plastic=0, fixed_stress=0, length=0, ldum=0, memory=-UPDATED, 
    idum[1], *element_tendon_number=NULL;
  double inc_tendon_ept=0., tendon_expansion=0., young=0., factor=0., 
    tendon_stress=0., ddum[1], new_element_tendon_strain[MTENDON], 
    new_element_tendon_stress[MTENDON], 
    direction[MDIM], tmp_vec[MDIM], tendon_force[MDIM*MDIM], 
    *tendon_elasti=NULL, *tendon_plasti=NULL, 
    *element_tendon_volume=NULL, *element_tendon_direction=NULL,
    *tendon_stress_time=NULL;

  array_set( ddsdde_tendon, 0., MSTRAIN*MSTRAIN );

  if ( db_active_index( ELEMENT_TENDON_NUMBER, element, VERSION_NORMAL ) ) {
    swit = set_swit(-1,-1,"tendons");
    if ( swit ) pri( "In routine TENDONS" );

    ntendon = db_len( ELEMENT_TENDON_NUMBER, element, VERSION_NORMAL );

    element_tendon_number = 
      db_int( ELEMENT_TENDON_NUMBER, element, VERSION_NORMAL );
    element_tendon_volume = 
      db_dbl( ELEMENT_TENDON_VOLUME, element, VERSION_NORMAL );
    element_tendon_direction = 
      db_dbl( ELEMENT_TENDON_DIRECTION, element, VERSION_NORMAL );

    db( ELEMENT_TENDON_STRAIN, element, idum, new_element_tendon_strain, 
      ntendon, VERSION_NORMAL, GET_AND_CHECK );
    db( ELEMENT_TENDON_STRESS, element, idum, new_element_tendon_stress, 
      ntendon, VERSION_NORMAL, GET_AND_CHECK );
    db( GROUP_MATERI_MEMORY, gr, &memory, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );

    array_set( tendon_force, 0., ndim*ndim );
    for ( iten=0; iten<ntendon; iten++ ) {
      array_move( &element_tendon_direction[iten*MDIM], direction, MDIM ); 
      matrix_ab( inc_ept, direction, tmp_vec, MDIM, MDIM, 1 );
      matrix_ab( tmp_vec, direction, &inc_tendon_ept, 1, MDIM, 1 );
      new_element_tendon_strain[iten] += inc_tendon_ept;
      itendon = element_tendon_number[iten];
      if ( db_active_index( TENDON_ELASTI, itendon, VERSION_NORMAL ) ) {
        tendon_elasti = db_dbl( TENDON_ELASTI, itendon, VERSION_NORMAL );
        young = tendon_elasti[0];
        new_element_tendon_stress[iten] += young * inc_tendon_ept;
        if ( condif_temperature ) {
          db( TENDON_EXPANSION, itendon, idum, &tendon_expansion, ldum, 
            VERSION_NORMAL, GET_IF_EXISTS );
          new_element_tendon_stress[iten] += -tendon_expansion * young * 
            ( new_unknowns[temp_indx] - old_unknowns[temp_indx] );
        }
      }
      plastic = 0;
      if ( db_active_index( TENDON_PLASTI, itendon, VERSION_NORMAL ) ) {
        tendon_plasti = db_dbl( TENDON_PLASTI, itendon, VERSION_NORMAL );
        if      ( new_element_tendon_stress[iten]>tendon_plasti[0] ) {
          plastic = 1;
          new_element_tendon_stress[iten] = tendon_plasti[0];
        }
        else if ( new_element_tendon_stress[iten]<-tendon_plasti[0] ) {
          plastic = 1;
          new_element_tendon_stress[iten] = -tendon_plasti[0];
        }
      }
      if ( db_active_index( TENDON_STRESS, itendon, VERSION_NORMAL ) &&
           db_active_index( TENDON_STRESS_TIME, itendon, VERSION_NORMAL ) ) {
        db( TENDON_STRESS, itendon, idum, &tendon_stress, ldum, 
          VERSION_NORMAL, GET );
        tendon_stress_time = db_dbl( TENDON_STRESS_TIME, itendon, VERSION_NORMAL );
        length = db_len( TENDON_STRESS_TIME, itendon, VERSION_NORMAL ); 
        fixed_stress = force_time( tendon_stress_time, "TENDON_STRESS_TIME", length, factor );
        if ( fixed_stress ) {
          new_element_tendon_stress[iten] = factor * tendon_stress;
          if ( young==0. ) 
             new_element_tendon_strain[iten] =  0.;
          else
             new_element_tendon_strain[iten] = new_element_tendon_stress[iten] / young;
          young = 0.;
        }
      }
      if      ( memory==-TOTAL )
        matrix_ab( new_rot, direction, tmp_vec, MDIM, MDIM, 1 );
      else if ( memory==-TOTAL_PIOLA )
        matrix_ab( new_rot, direction, tmp_vec, MDIM, MDIM, 1 );
      else if ( memory==-UPDATED ) {
        pri( "Sorry: -UPDATED lagrange memory not available for tendons." );
        pri( "Use -UPDATED_WITHOUT_ROTATION." );
        exit(TN_EXIT_STATUS);
      }
      else {
        assert( memory==-UPDATED_WITHOUT_ROTATION );
        array_move( direction, tmp_vec, MDIM );
      }
      for ( idim=0; idim<ndim; idim++ ) {
        for ( jdim=0; jdim<ndim; jdim++ ) {
          tendon_force[idim*ndim+jdim] +=  new_element_tendon_stress[iten] * 
            element_tendon_volume[iten] * tmp_vec[idim] * tmp_vec[jdim];
          if ( !plastic && !fixed_stress ) {
            for ( kdim=0; kdim<ndim; kdim++ ) {
              for ( ldim=0; ldim<ndim; ldim++ ) {
                i = stress_indx(idim,jdim);
                j = stress_indx(kdim,ldim);
                ind_ddsdde = i*MSTRAIN+j;
                ddsdde_tendon[ind_ddsdde] += young * 
                 (element_tendon_volume[iten]/(npoint*volume)) *
                 tmp_vec[idim] * tmp_vec[jdim] * tmp_vec[kdim] * tmp_vec[ldim];
              }
            }
          }
        }
      }
    }
    parallel_sys_lock();
    db( ELEMENT_TENDON_STRAIN, element, idum, new_element_tendon_strain, 
      ntendon, VERSION_NEW, PUT );
    db( ELEMENT_TENDON_STRESS, element, idum, new_element_tendon_stress, 
      ntendon, VERSION_NEW, PUT );
    parallel_sys_unlock();
    for ( inol=0; inol<nnol; inol++ ) {
      for ( idim=0; idim<ndim; idim++ ) {
        ipuknwn = vel_indx/nder+idim;
        indx = inol*npuknwn + ipuknwn;
        for ( jdim=0; jdim<ndim; jdim++ ) {
          tendon_element_rhside[indx] -= new_d[jdim*nnol+inol] * 
            tendon_force[idim*ndim+jdim];
        }
      }
    }
    if ( swit ) {
      pri( "new_element_tendon_stress", new_element_tendon_stress, ntendon );
      pri( "tendon_force", tendon_force, ndim, ndim );
      pri( "tendon_element_rhside", tendon_element_rhside, nnol, npuknwn );
      pri( "ddsdde_tendon", ddsdde_tendon, MSTRAIN, MSTRAIN );
    }

    if ( swit ) pri( "Out routine TENDONS" );
  }

}
