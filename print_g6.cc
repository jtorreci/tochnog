/***************************************************************************
                          print_g6.cc  -  description
                             -------------------
    begin                : Thu Apr 4 2002
    copyright            : (Cl) 2002 by Jesús Torrecilla
    email                : jtorreci@unex.es
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

/*
    copyright (c) 1998  dennis roddeman
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

/* Comment: Jesús Torrecilla Pinero
            I have begun with just the code from Dennis, cleaning the code
						for Gid 6
*/
#include "tochnog.h"
#define MTYPES 11 //In print_g5.cc this is set to 10, so, no trussbeam are printed


int print_coordinates (ofstream & outmesh)
{
double coord_start_refined[MDIM];
long int ldum=0, idum[1], max_node=0;
int inod,idim;

 db_highest_index( NODE, max_node, VERSION_PRINT );
 outmesh <<"# node number   coordinate_x  coordinate_y  coordinate_z"<<endl;	      		
   for ( inod=1; inod<=max_node; inod++ )
   {
	if ( db_active_index( NODE_START_REFINED, inod, VERSION_PRINT ) )
	{
 		outmesh << inod << "\t";      		
 		db( NODE_START_REFINED, inod, idum, coord_start_refined, ldum, VERSION_PRINT, GET );
 		for ( idim=0; idim<ndim; idim++ )
   		outmesh << coord_start_refined[idim] << "\t";      						
 		outmesh << "\n";
    	} //end if
   }	//end for	

return 1;

}


long int print_gid_6( long int task )

{
  long int inod=0, element=0, max_node=0, max_element=0, max_group=0, name=0, itype,
    eltype[MTYPES],
    element_group=0, length=0, idim=0, jdim=0, ipuknwn=0, iuknwn=0,
    icalcul=0, nval=0, ready=0, indx=0, eigen_analysis=0, neigen=0,
    length_post_calcul_scal_vec_mat=0,
    load_type=0, data_type=0, data_loc=0, desc_comp=0,
    n=0, swit=0, ldum=0, nnol=0, ieigen=0,
    bound_elements=0,
    nmesh=0,first_time_gid=0,
    print_deformation=0, element_empty=-NO, icontrol=0,
    control_print_gid_mesh=0, idum[1],
    *dof_amount=NULL, *dof_label=NULL, *dof_type=NULL,
    *post_calcul_scal_vec_mat=NULL, *dof_scal_vec_mat=NULL,
    *nodes=NULL, *el=NULL;
  double time_current=0., control_print_gid_time=0.,
     step_val=0., tmp=0., ddum[1], coord[MDIM],
     coord_start_refined[MDIM],
     *element_tendon_intersections=NULL, *control_eigen_values=NULL,
     *node_eigen=NULL, *node_dof=NULL,
     *node_dof_calcul=NULL;
  char str[MCHAR], descr_menu[MCHAR], filename[MCHAR],
     filename_with_number[MCHAR], outname[MCHAR];

  if ( ndim==1 ) return 1;

  swit = set_swit(-1,-1,"print_gid_6");
  if ( swit ) pri( "In routine PRINT_GID_6" );

  if ( task!=-YES && task!=-SEPARATE ) {
    pri( "Error detected in gid 6 printing." );
    pri( "Gid 6 printing doesn't understand ", task );
    exit_tn_on_error();
  }

  dof_amount = get_new_int(MUKNWN);
  dof_label = get_new_int(MUKNWN);
  dof_type = get_new_int(MUKNWN);
  dof_scal_vec_mat = get_new_int(MUKNWN);
  post_calcul_scal_vec_mat = get_new_int(DATA_ITEM_SIZE);
  nodes = get_new_int(MAXIMUM_NODE);
  el = get_new_int(MAXIMUM_NODE+1);
  element_tendon_intersections = get_new_dbl(DATA_ITEM_SIZE);
  control_eigen_values = get_new_dbl(DATA_ITEM_SIZE);

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( TIME_CURRENT, 0, idum, &time_current, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  if ( db( CONTROL_PRINT_GID_TIME, 0, idum, &control_print_gid_time,
       ldum, VERSION_NORMAL, GET_IF_EXISTS ) ) {
    if ( time_current<=control_print_gid_time && task!=-SEPARATE ) return 1;
  }
  else {
    first_time_gid = 1;
  }
  length = 1;
  db( CONTROL_PRINT_GID_TIME, 0, idum, &time_current, length, VERSION_NORMAL, PUT );

  db_copy( NODE_RHSIDE, NODE_RHSIDE_PRINT, VERSION_NORMAL );
  db_version_copy( VERSION_NORMAL, VERSION_PRINT );

  renumbering( VERSION_PRINT, NO, 1, 1, idum, idum );

  db_highest_index( NODE, max_node, VERSION_PRINT );
  if ( max_node<0 ) return 1;
  db_highest_index( ELEMENT, max_element, VERSION_PRINT );
  if ( max_element<0 ) return 1;
  db_highest_index( ELEMENT_GROUP, max_group, VERSION_PRINT );
  //if (max_group<0) return 1;
  if (max_group<0) max_group=0;

  if ( db( CONTROL_EIGEN_VALUES, 0, idum, control_eigen_values,
    neigen, VERSION_NORMAL, GET_IF_EXISTS ) ) eigen_analysis = 1;
  else {
    eigen_analysis = 0;
    neigen = 1;
  }

  strcpy( filename_with_number, data_file_base );
  if ( task==-SEPARATE ) {
    db( CONTROL_PRINT_GID_MESH, 0, &control_print_gid_mesh, ddum,
      length, VERSION_NORMAL, GET_IF_EXISTS );
    if      ( control_print_gid_mesh<10 )
      strcat( filename_with_number, "00" );
    else if ( control_print_gid_mesh<100 )
      strcat( filename_with_number, "0" );
    long_to_a( control_print_gid_mesh, str );
    strcat( filename_with_number, str );
    control_print_gid_mesh++; length=1;
    db( CONTROL_PRINT_GID_MESH, 0, &control_print_gid_mesh, ddum,
      length, VERSION_NORMAL, PUT );
  }

  /*If number of dimensions == 2 or number of dimensions == 3, anyway, the mesh
    is stored in  .flavia.msh file. The stream for this file is outmesh.
    Next, create, or open the files for postprocessing.*/
  strcpy( filename, filename_with_number );
  strcat( filename, ".flavia.msh" );
  if ( first_time_gid || task==-SEPARATE ) {
    ofstream outmesh( filename );
    outmesh.close();
  }
  ofstream outmesh( filename, ios::app );
  outmesh.precision(TN_PRECISION);

  strcpy( filename, filename_with_number );
  strcat( filename, ".flavia.res" );
  if ( first_time_gid || task==-SEPARATE ) {
    ofstream outres( filename );
    outres << "GiD Post Results File 1.0" <<endl;
    outres.close();
  }
  ofstream outres( filename, ios::app );
  outres.precision(TN_PRECISION);
  outres.setf(ios::showpoint);

    /* In the original prin_g5.cc, here we determined the type of the isoparametric
       elements in mesh, but, in fact, there may be more than one type, so, in this
       .cc file, we first write all the .msh file, with all the elements and then we
       write the .res file */

int number_of_meshes=0;
eltype[0]=  -TRIA3;
eltype[1]=  -TRIA6;
eltype[2]=  -QUAD4;
eltype[3]=  -QUAD9;
eltype[4]=  -TET4;
eltype[5]=  -TET10;
eltype[6]=  -HEX8;
eltype[7]=  -HEX27;
eltype[8]=  -TRUSS;
eltype[9]=  -BEAM;
eltype[10]= -TRUSSBEAM;

if(first_time_gid || task==-SEPARATE) {//prints mesh only once 
for (itype=0; itype<MTYPES; itype++)
for (int igroup=0; igroup<=max_group; igroup++)
{
{
  nmesh=0;
  for (element=1; element <=max_element; element ++)
  {
      //Does the element have to be printed?
     if ( db_active_index(  ELEMENT, element, VERSION_PRINT ) )
     {
      	db( ELEMENT, element, el, ddum, length, VERSION_PRINT, GET );
      name = el[0];
      if ( db_active_index( ELEMENT_GROUP, element, VERSION_PRINT ) )
        db( ELEMENT_GROUP, element, &element_group, ddum, ldum, VERSION_PRINT, GET );
      else
        element_group = 0;
      if ((name==eltype[itype])&&(element_group==igroup))
      {
	    element_empty=-NO;
  	    db( ELEMENT_EMPTY, element, &element_empty, ddum,
    	    ldum, VERSION_PRINT, GET_IF_EXISTS );
      	if ( element_empty==-NO || element_empty==-FRONT )
      	{
         	nnol = length - 1; array_move( &el[1], nodes, nnol );         	    	
         		nmesh++;
	      	if (nmesh==1)
	      	{
	      	  //If we still haven't writen any mesh, we have to write the nodal coordinates
	      	  //just now.
           	switch (eltype[itype])
           	{
           		case -TRIA3:
           				 outmesh <<"MESH Tria-3Gr"<< igroup<<" dimension "<<ndim;
		           		 outmesh <<" ElemType Triangle Nnode 3"<<endl;
			      		 outmesh <<" Coordinates" << endl; break;
           		case -TRIA6:
              				 outmesh <<"MESH Tria-6Gr"<< igroup<<" dimension "<<ndim;
		           		 outmesh <<" ElemType Triangle Nnode 6"<<endl;
			      		 outmesh <<" Coordinates" << endl; break;
     			case -QUAD4:
           				 outmesh <<"MESH Quad-4Gr"<< igroup<<" dimension "<<ndim;
		           		 outmesh <<" ElemType Quadrilateral Nnode 4"<<endl;
			      		 outmesh <<" Coordinates" << endl;break;       						
     			case -QUAD9:
           				 outmesh <<"MESH Quad-9Gr"<< igroup<<" dimension "<<ndim;
		           		 outmesh <<" ElemType Quadrilateral Nnode 9"<<endl;
			      		 outmesh <<" Coordinates" << endl;break;
                    case -TET4:
           				 outmesh <<"MESH Tetra-3Gr"<< igroup<<" dimension "<<ndim;
		           		 outmesh <<" ElemType Tetrahedra Nnode 4"<<endl;
			      		 outmesh <<" Coordinates" << endl;break;
                    case -TET10:
           				 outmesh <<"MESH Tetra-10Gr"<< igroup<<" dimension "<<ndim;
		           		 outmesh <<" ElemType Tetrahedra Nnode 10"<<endl;
			      		 outmesh <<" Coordinates" << endl; break;       		
                    case -HEX8:
            				 outmesh <<"MESH Hexa-8Gr"<< igroup<<" dimension "<<ndim;
		           		 outmesh <<" ElemType Hexahedra Nnode 8"<<endl;
			      		 outmesh <<" Coordinates" << endl; break;
                    case -HEX27:
             				 outmesh <<"MESH Hexa-27Gr"<< igroup<<" dimension "<<ndim;
		           		 outmesh <<" ElemType Hexahedra Nnode 27"<<endl;
			      		 outmesh <<" Coordinates" << endl; break;
                    case -TRUSS:
            				 outmesh <<"MESH Truss-Gr"<<igroup<<" dimension "<<ndim;
		           		 outmesh <<" ElemType Linear Nnode 2"<<endl;
			      		 outmesh <<" Coordinates" << endl; break;
                    case -BEAM:
           				 outmesh <<"MESH BeamGr"<<igroup<<" dimension "<<ndim;
		           		 outmesh <<" ElemType Linear Nnode 2"<<endl;
			      		 outmesh <<" Coordinates" << endl; break;
                    case -TRUSSBEAM:
           				 outmesh <<"MESH TrussBeamGr"<< igroup<<" dimension "<<ndim;
		           		 outmesh <<" ElemType Linear Nnode 2"<<endl;
			      		 outmesh <<" Coordinates" << endl; break;
			      		
			     }		
	      		if (number_of_meshes==0)
	      		{	      		
	        		  number_of_meshes=1;
	      		  print_coordinates (outmesh);
	      		 }
				outmesh <<"end coordinates" <<endl << endl;
		      	outmesh<<"Elements" <<endl;
		      	outmesh<<"# element  node_1   node_2  node_3 ... material_number	"<<endl;
	      		} //end if number_of_meshes
      	outmesh << element << "\t";
      	switch (eltype[itype])
      	{
      		case -TRIA3:
      					          outmesh << nodes[0] << "\t" << nodes[1] << "\t" << nodes[2] << "\t";break;
      		case -TRIA6:
							      outmesh << nodes[0] << "\t" << nodes[2] << "\t" << nodes[5] << "\t " <<
						                     nodes[1] << "\t" << nodes[4] << "\t" << nodes[3] << "\t"; break;
			case -QUAD4:
						          outmesh << nodes[0] << "\t" << nodes[1] << "\t" << nodes[3] << "\t" << nodes[2] << "\t";break;
			case -QUAD9:
                                  outmesh << nodes[0] << "\t" << nodes[2] << "\t" << nodes[8] << "\t" <<
                                  nodes[6] << "\t" << nodes[1] << "\t" << nodes[5] << "\t" <<
                                  nodes[7] << "\t" << nodes[3] << "\t" << nodes[4] << "\t"; break;
            case -TET4:
                                  outmesh << nodes[0] << "\t" << nodes[1] << "\t" << nodes[2] << "\t" << nodes[3] << "\t"; break;
            case -TET10:
                                  outmesh << nodes[0] << "\t" << nodes[2] << "\t" << nodes[5] << "\t" << nodes[9] << "\t";
                                  outmesh << nodes[1] << "\t" << nodes[4] << "\t" << nodes[3] << "\t";
                                  outmesh << nodes[6] << "\t" << nodes[7] << "\t" << nodes[8] << "\t"; break;
 

/*                              outmesh << nodes[0] << "\t" << nodes[1] << "\t" << nodes[2] << "\t" << nodes[3] << "\t";
                                outmesh << nodes[4] << "\t" << nodes[5] << "\t" << nodes[6] << "\t";
                                outmesh << nodes[7] << "\t" << nodes[8] << "\t" << nodes[9] << "\t"; break;
 */
            case -HEX8:
						          outmesh << nodes[0] << "\t" << nodes[1] << "\t" << nodes[3] << "\t" << nodes[2] << "\t" <<
                                  nodes[4] << "\t" << nodes[5] << "\t" << nodes[7] << "\t" << nodes[6] << "\t"; break;

               case -HEX27:
    /*in g6 by jesus, wrong order
							  outmesh << nodes[0]  << "\t" << nodes[2]  << "\t" << nodes[8]  << "\t" << nodes[6]  << "\t";
							  outmesh << nodes[18] << "\t" << nodes[20] << "\t" << nodes[26] << "\t" << nodes[24] << "\t";
                              outmesh << nodes[1]  << "\t" << nodes[5]  << "\t" << nodes[7]  << "\t" << nodes[3]  << "\t";
                              outmesh << nodes[9]  << "\t" << nodes[11] << "\t" << nodes[17] << "\t" << nodes[15] << "\t";
                              outmesh << nodes[19] << "\t" << nodes[23] << "\t" << nodes[25] << "\t" << nodes[22] << "\t";
                              outmesh << nodes[4] ;
                              outmesh << nodes[10] << "\t" << nodes[14] << "\t" << nodes[16] << "\t" << nodes[12] << "\t";
                              outmesh << nodes[22] << "\t" << nodes[13] << "\t" ;break;
			*/
                    outmesh << nodes[20] << " " << nodes[26] << " " << nodes[24] << " " << nodes[18] << " " ;
                    outmesh << nodes[2]  << " " << nodes[8]  << " " << nodes[6]  << " " << nodes[0]  << " " ;
                    outmesh << nodes[23] << " " << nodes[25] << " " << nodes[21] << " " << nodes[19] << " " ;
                    outmesh << nodes[11] << " " << nodes[17] << " " << nodes[15] << " " << nodes[9]  << " " ;
                    outmesh << nodes[5]  << " " << nodes[7]  << " " << nodes[3]  << " " << nodes[1]  << " " ;
                    outmesh << nodes[22] << " " << nodes[14] << " " << nodes[16] << " " << nodes[12] << " " ;
                    outmesh << nodes[10] << " " << nodes[4]  << " " << nodes[13] << " "; break;

            case -TRUSS:
			                         outmesh << nodes[0] << "\t" << nodes[1] << "\t"; break;
            case -BEAM:
			                         outmesh << nodes[0] << "\t" << nodes[1] << "\t"; break;
			case -TRUSSBEAM:
			                         outmesh << nodes[0] << "\t" << nodes[1] << "\t"; break;
	
      	}
    	     outmesh << element_group << endl;
      	}
      }
	}
   }
 if (nmesh)  outmesh << "end elements" << endl;
}
}
}//end first time gid

  if ( bound_elements )
    strcpy( filename, "lines.bon" );
  else
    strcpy( filename, "tn.tmp" );
  if ( first_time_gid ) {
    ofstream outbound( filename );
    outbound.close();
  }
  ofstream outbound( filename, ios::app );


    // tendon connectivity. We have to rewrite this.

    // Writing of the results. Stream outres.

  load_type = 1; //For timestep analysis in old postprocess format
  data_loc = 1;  //On Nodes
  for ( ieigen=0; ieigen<neigen; ieigen++ ) {
    if ( eigen_analysis ) step_val = ieigen;
    else step_val = time_current;

    if ( npuknwn>0 ) {

      db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET );
      db( DOF_TYPE, 0, dof_type, ddum, ldum, VERSION_NORMAL, GET );
      db( DOF_AMOUNT, 0, dof_amount, ddum, ldum, VERSION_NORMAL, GET );
      db( DOF_SCAL_VEC_MAT, 0, dof_scal_vec_mat, ddum, ldum, VERSION_NORMAL, GET );

    // write scalars, vectors and tensors for primary unknowns
      ipuknwn = 0; ready = 0;
      while ( !ready ) {
        iuknwn = ipuknwn*nder;
        if      ( print_deformation )
        {
        	outres <<"Result "<<"\"mesh deform\""<<"\t"<<"\"Load Analysis\"";
        	outres <<"\t"<<step_val<<"\t Vector  OnNodes"<<endl;
        	outres << "ComponentNames ";
            nval = ndim;
            data_type = 2;
            desc_comp = 1;
        for ( idim=0; idim<ndim; idim++ )
          {
            if      ( idim==0 )
              outres << "\"mesh_deform-x\"" << ", ";
            else if ( idim==1 )
              outres << "\"mesh_deform-y\"" << ", ";
            else
                {
              assert( idim==2 );
              outres << "\"mesh_deform-z\"" << ", ";
                }
          }
          outres << endl;
        }
        else if ( dof_scal_vec_mat[iuknwn]==-SCALAR )
        {
          nval = 1;
          data_type = 1;
          desc_comp = 1;
          strcpy( descr_menu, db_name(dof_label[iuknwn]) );
          string_shorten( descr_menu, 15 );
          outres <<"Result  \""<<descr_menu<<"\"     \""<<"Load Analysis"<< "\" \t"<<step_val<<" Scalar OnNodes"<<endl;
        }
        else if ( dof_scal_vec_mat[iuknwn]==-VECTOR )
        {
          nval = dof_amount[iuknwn];
          data_type = 2;
          desc_comp = 1;
          if ( dof_type[iuknwn]==-MATERI_VELOCITY_INTEGRATED )
            strcpy( descr_menu, "materi_velint" );
          else
            strcpy( descr_menu, db_name(dof_type[iuknwn]) );
          string_shorten( descr_menu, 15 );
          outres <<"Result  \""<<descr_menu<<"\" \t  \""<<"Load Analysis"<< "\" \t"<<step_val<<" Vector OnNodes"<<endl;
          outres <<"ComponentNames ";
          for ( idim=0; idim<nval; idim++ )
          {
            indx = iuknwn + idim*nder;
            outres << "\""<<db_name(dof_label[indx]) <<"\" ";
            if (idim!=(nval-1)) outres <<",  "; else outres <<endl;
          }
        }
        else
        {
          assert( dof_scal_vec_mat[iuknwn]==-MATRIX );
          nval = 6;
          data_type = 3;
          desc_comp = 1;
          strcpy( descr_menu, db_name(dof_type[iuknwn]) );
          string_shorten( descr_menu, 15 );
          outres <<"Result  \""<<descr_menu<<"\"     \""<<"Load Analysis"<< "\" \t"<<step_val<<" Matrix OnNodes"<<endl;
          outres <<"ComponentNames " ;
          //First, main diagonal, and then the rest of diagonals going down
           for ( idim=0; idim<MDIM; idim++ ) {
            for ( jdim=idim; jdim<MDIM; jdim++ )
            {
              indx = iuknwn + stress_indx(jdim,jdim-idim)*nder;
              outres << "\""<<db_name(dof_label[indx])<<"\", \t";
            }
          }
        }
        outres << endl;
        outres <<"Values"<<endl;
        for ( inod=1; inod<=max_node; inod++ ) {
          outres << inod << "\t";
          if ( eigen_analysis )
            node_eigen = db_dbl( NODE_EIGEN, inod, VERSION_PRINT );
          else
            node_dof = db_dbl( NODE_DOF, inod, VERSION_PRINT );
          if      ( print_deformation ) {
            db( NODE, inod, idum,
              coord, ldum, VERSION_PRINT, GET );
            db( NODE_START_REFINED, inod, idum,
              coord_start_refined, ldum, VERSION_PRINT, GET );
            for ( idim=0; idim<ndim; idim++ ) {
              if ( materi_displacement )
                tmp = coord[idim] + node_dof[dis_indx+idim*nder];
              else
                tmp = coord[idim];
              outres << tmp-coord_start_refined[idim] << "\t";
            }
          }
          else if ( dof_scal_vec_mat[iuknwn]==-SCALAR ) {
            if ( eigen_analysis )
              tmp = node_eigen[ieigen*nuknwn+iuknwn];
            else
              tmp = node_dof[iuknwn];
            outres << tmp;
          }
          else if ( dof_scal_vec_mat[iuknwn]==-VECTOR ) {
            n = dof_amount[iuknwn];
            for ( idim=0; idim<n; idim++ ) {
              indx = iuknwn+idim*nder;
              if ( eigen_analysis )
                tmp = node_eigen[ieigen*nuknwn+indx];
              else
                tmp = node_dof[indx];
              outres << tmp << "\t";
            }
          }
          else {
            assert( dof_scal_vec_mat[iuknwn]==-MATRIX );
            for ( idim=0; idim<MDIM; idim++ ) {
              for ( jdim=idim; jdim<MDIM; jdim++ ) {
                indx = iuknwn + stress_indx(jdim,jdim-idim)*nder;
                if ( eigen_analysis )
                  tmp = node_eigen[ieigen*nuknwn+indx];
                else
                  tmp = node_dof[indx];
                outres << tmp << "\t";
              }
            }
          }
          outres << endl;
        }
        outres <<"End Values" <<endl;
        ipuknwn += nval;
        if ( ipuknwn>=npuknwn ) {
          if      ( materi_velocity && !print_deformation ) {
            print_deformation = 1;
          }
          else
            ready = 1;
        }
      }

      if ( !eigen_analysis && db_active_index( POST_CALCUL, 0, VERSION_NORMAL ) ) {

        db( POST_CALCUL_SCAL_VEC_MAT, 0, post_calcul_scal_vec_mat,
          ddum, length_post_calcul_scal_vec_mat, VERSION_NORMAL, GET );

    // write scalars, vectors and tensors for calculated data
        icalcul = idim = ready = 0;
        while ( !ready ) {
          strcpy( outname, post_calcul_names_without_extension[icalcul] );
          if      ( post_calcul_scal_vec_mat[icalcul]==-SCALAR ) {
            nval = 1;
            data_type = 1;
            desc_comp = 0;
            strcpy( descr_menu, outname );
            string_shorten( descr_menu, 15 );
            outres <<"Result  \""<<descr_menu<<"\"     \""<<"Load Analysis"<< "\" \t"<<step_val<<" Scalar OnNodes"<<endl;
          }
          else if ( post_calcul_scal_vec_mat[icalcul]==-VECTOR ) {
            nval = MDIM;
            data_type = 2;
            desc_comp = 1;
            strcpy( descr_menu, outname );
            string_shorten( descr_menu, 15 );
            outres <<"Result  \""<<descr_menu<<"\"     \""<<"Load Analysis"<< "\"\t"<<step_val<<" Vector OnNodes"<<endl;
            outres <<"ComponentNames " ;
            for ( idim=0; idim<MDIM; idim++ ) {
              outres << "\""<<post_calcul_names[icalcul+idim]<<"\", \t";
            }
          }
          else {
            assert( post_calcul_scal_vec_mat[icalcul]==-MATRIX );
            nval = 6;
            data_type = 3;
            desc_comp = 1;
            strcpy( descr_menu, outname );
            string_shorten( descr_menu, 15 );
            outres <<"Result  \""<<descr_menu<<"\"     \""<<"Load Analysis"<< "\"\t"<<step_val<<" Matrix OnNodes"<<endl;
            outres <<"ComponentNames " ;
            for ( idim=0; idim<MDIM; idim++ ) {
              for ( jdim=idim; jdim<MDIM; jdim++ ) {
                indx = icalcul + stress_indx(jdim,jdim-idim)*nder;
                outres << "\""<<post_calcul_names[icalcul+idim]<<"\", \t";
              }
            }
          }
          outres <<endl<<"Values "<<endl;
          for ( inod=1; inod<=max_node; inod++ ) {
            outres << inod << "\t";
            node_dof_calcul = db_dbl( NODE_DOF_CALCUL, inod, VERSION_PRINT );
            if      ( post_calcul_scal_vec_mat[icalcul]==-SCALAR ) {
              tmp = node_dof_calcul[icalcul];
              outres << tmp;
            }
            else if ( post_calcul_scal_vec_mat[icalcul]==-VECTOR ) {
              for ( idim=0; idim<MDIM; idim++ ) {
                indx = icalcul+idim;
                tmp = node_dof_calcul[indx];
                outres << tmp << "\t";
              }
            }
            else {
              assert( post_calcul_scal_vec_mat[icalcul]==-MATRIX );
              for ( idim=0; idim<MDIM; idim++ ) {
                for ( jdim=idim; jdim<MDIM; jdim++ ) {
                  indx = icalcul + stress_indx(jdim,jdim-idim)*nder;
                  tmp = node_dof_calcul[indx];
                  outres << tmp << "\t";
                }
              }
            }
            outres  << endl;
          }
          outres <<"End Values "<<endl;
          icalcul += nval;
          ready = (icalcul>=length_post_calcul_scal_vec_mat);
        }

      }

    }

  }

  delete[] dof_amount;
  delete[] dof_label;
  delete[] dof_type;
  delete[] post_calcul_scal_vec_mat;
  delete[] dof_scal_vec_mat;
  delete[] nodes;
  delete[] el;
  delete[] element_tendon_intersections;
  delete[] control_eigen_values;

  outmesh.close();
  outres.close();
  db_version_delete( VERSION_PRINT );

  if ( swit ) pri( "Out routine PRINT_GID_6" );

  return 1;
}





