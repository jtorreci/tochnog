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

void nonlocal_set( void )

{
  long int i=0, jn=0, inod=0, jnod=0, kn=0, knod=0, ldum=0, 
    length=0, length_node_node=0, max_node=0, nnonlocal=0, ready=0, 
    swit=0, idum[1], *node_node=NULL, *node_nonlocal=NULL,
    *element_nonlocal=NULL, max_element=0, ielem=0, kelem=0, ipoint_ielem=0, npoint_ielem=0,
    ipoint_kelem=0, npoint_kelem=0, *element_nonlocal_ipoint=NULL, kelem_nonloc=1, ielem_nonloc=1;
  double x=0., normal_distribution=0., distance=0., 
    options_nonlocal=0., total_weight=0., 
    ddum[1], icoord[MDIM], kcoord[MDIM], work[MDIM], 
    *node_nonlocal_weight=NULL, bell_function=0,
    *element_nonlocal_weight=NULL, 
    *nonlocal_ielem_info=NULL, *nonlocal_kelem_info=NULL, kelem_ipoint_volume=0;

  if ( materi_plasti_f_nonlocal ) {
    swit = set_swit(-1,-1,"nonlocal_set");
    if ( swit ) pri( "In routine NONLOCAL_SET" );

    length = db_data_length(NODE_NODE);
    node_node = get_new_int(length);
    node_nonlocal = get_new_int(NONLOCAL_ITEM_SIZE);
    node_nonlocal_weight = get_new_dbl(NONLOCAL_ITEM_SIZE);
    db_delete( NODE_NONLOCAL, VERSION_NORMAL );
    db_delete( NODE_NONLOCAL_WEIGHT, VERSION_NORMAL );
    db_max_index( NODE, max_node, VERSION_NORMAL, GET );
    db( OPTIONS_NONLOCAL, 0, idum, &options_nonlocal, ldum, 
      VERSION_NORMAL, GET_IF_EXISTS );
    if(scalar_dabs(options_nonlocal_softvar)>0) options_nonlocal=options_nonlocal_softvar;
    
    for ( inod=0; inod<=max_node; inod++ ) {
      if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
          // search for neighbours of node inod
        db( NODE, inod, idum, icoord, ldum, VERSION_NORMAL, GET );
        ready=0; nnonlocal = 0;
        node_nonlocal[nnonlocal] = inod;
        node_nonlocal_weight[nnonlocal] = 1.;
        nnonlocal++;
        while ( !ready ) {
          ready = 1;
            // loop over already found neighbours jnod
          for ( jn=0; jn<nnonlocal; jn++ ) {
            jnod = node_nonlocal[jn];
              // loop over nodes knod connected to neighbours jnod
            db( NODE_NODE, jnod, node_node, ddum, length_node_node, VERSION_NORMAL, GET );
            for ( kn=0; kn<length_node_node; kn++ ) {
              knod = labs(node_node[kn]);
                // add this node knod to inod list if its distance is small enough
              db( NODE, knod, idum, kcoord, ldum, VERSION_NORMAL, GET );
              distance = array_distance( icoord, kcoord, work, ndim );
              if ( distance<=options_nonlocal && 
                   !array_member(node_nonlocal,knod,nnonlocal,ldum) ) {
                ready = 0;
                node_nonlocal[nnonlocal] = knod;
                x = 2.*distance/options_nonlocal;
                normal_distribution = (1./(sqrt(2.*PIRAD))) * exp(-x*x/(2.));
                node_nonlocal_weight[nnonlocal] = normal_distribution;
                nnonlocal++;
                if ( nnonlocal==NONLOCAL_ITEM_SIZE ) {
                  pri( "NONLOCAL_ITEM_SIZE in tochnog.h too small for nonlocal calculation." );
                  pri( "Increase it and recompile all routines." );
                  exit(TN_EXIT_STATUS);
                }
              }
            }
          }
        }
        db( NODE_NONLOCAL, inod, node_nonlocal, ddum, 
          nnonlocal, VERSION_NORMAL, PUT );
        db( NODE_NONLOCAL_WEIGHT, inod, idum, node_nonlocal_weight, 
          nnonlocal, VERSION_NORMAL, PUT );
      }
    }
      // normalize the weight factors to 1
    for ( inod=0; inod<=max_node; inod++ ) {
      if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
        db( NODE_NONLOCAL_WEIGHT, inod, idum, node_nonlocal_weight, 
          nnonlocal, VERSION_NORMAL, GET );
        total_weight = 0.;
        for ( i=0; i<nnonlocal; i++ ) total_weight += node_nonlocal_weight[i];
        assert( total_weight!=0. );
        for ( i=0; i<nnonlocal; i++ ) node_nonlocal_weight[i] /= total_weight;
        db( NODE_NONLOCAL_WEIGHT, inod, idum, node_nonlocal_weight, 
          nnonlocal, VERSION_NORMAL, PUT );
      }
    }
    delete[] node_node;
    delete[] node_nonlocal;
    delete[] node_nonlocal_weight;
    if ( swit ) pri( "Out routine NONLOCAL_SET" );
  }

/****************Averages values directly in integration points************************/
/*                    Programmed by D.Masin@city.ac.uk                                */
	
  if ( materi_plasti_softvar_nonlocal && scalar_dabs(options_nonlocal_softvar)>TINY) {
    cout<<"Set nonlocal weights"<<endl;
    swit = set_swit(-1,-1,"nonlocal_set");
    if ( swit ) pri( "In routine NONLOCAL_SET" );
    if(options_element_dof != -YES) {
    	cout<<"Please set 'options_element_dof -yes' for nonlocal calculation"<<endl;
	exit_tn( -YES ); 
    }

    long int length_nei=1+npointmax*ndim+npointmax+2;
    long int length_nonloc=NONLOCAL_ITEM_SIZE*npointmax;
    element_nonlocal = get_new_int(NONLOCAL_ITEM_SIZE*npointmax);
    element_nonlocal_ipoint = get_new_int(NONLOCAL_ITEM_SIZE*npointmax);
    element_nonlocal_weight = get_new_dbl(NONLOCAL_ITEM_SIZE*npointmax);
    nonlocal_ielem_info = get_new_dbl(length_nei);
    nonlocal_kelem_info = get_new_dbl(length_nei);
    db_delete( ELEMENT_NONLOCAL, VERSION_NORMAL );
    db_delete( ELEMENT_NONLOCAL_IPOINT, VERSION_NORMAL );
    db_delete( ELEMENT_NONLOCAL_WEIGHT, VERSION_NORMAL );
    db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );

    find_nonlocal_weights=1;
    element_loop();	//fills nonlocal_element_info
    find_nonlocal_weights=0;

    for ( ielem=0; ielem<=max_element; ielem++ ) {
     if ( db_active_index( ELEMENT, ielem, VERSION_NORMAL )) {
      array_set( element_nonlocal, -10, NONLOCAL_ITEM_SIZE*npointmax );
      array_set( element_nonlocal_ipoint, -10, NONLOCAL_ITEM_SIZE*npointmax );
      array_set( element_nonlocal_weight, -10, NONLOCAL_ITEM_SIZE*npointmax );
       
      db( NONLOCAL_ELEMENT_INFO, ielem, idum, nonlocal_ielem_info, 
 	length_nei, VERSION_NORMAL, GET );		
      npoint_ielem=(long int)nonlocal_ielem_info[0];
      ielem_nonloc=(long int)nonlocal_ielem_info[length_nei-1];
      if(ielem_nonloc) {

       for ( ipoint_ielem=0; ipoint_ielem<npoint_ielem; ipoint_ielem++ ) {
        // search for neighbours of integration point ipoint

        for(int idim=0; idim<ndim; idim++) 
	    icoord[idim]=nonlocal_ielem_info[1+ipoint_ielem*ndim+idim];
	nnonlocal = 0;
        
           // loop over all elements
        for ( kelem=0; kelem<=max_element; kelem++ ) {
	  if ( db_active_index( ELEMENT, kelem, VERSION_NORMAL )) {
	    db( NONLOCAL_ELEMENT_INFO, kelem, idum, nonlocal_kelem_info, 
	     length_nei, VERSION_NORMAL, GET );		
            npoint_kelem=(long int)nonlocal_kelem_info[0];
            kelem_nonloc=(long int)nonlocal_kelem_info[length_nei-1];

	    if(kelem_nonloc) {
              for ( ipoint_kelem=0; ipoint_kelem<npoint_kelem; ipoint_kelem++ ) {
	        for(int idim=0; idim<ndim; idim++) 
		    kcoord[idim]=nonlocal_kelem_info[1+ipoint_kelem*ndim+idim];
	        kelem_ipoint_volume=nonlocal_kelem_info[1+npoint_kelem*ndim+ipoint_kelem];
                distance = array_distance( icoord, kcoord, work, ndim );
                if ( distance<=options_nonlocal_softvar ) {

		    // I use bell-shaped function with bounded support
       	          if(distance>=0 && distance<=options_nonlocal_softvar)
		    bell_function=scalar_power((1-((distance*distance)/
		    	(options_nonlocal_softvar*options_nonlocal_softvar))) ,2);
 		  else bell_function=0;
		  if(distance==0) {	// If integration points same coordinate, take only the actual one
			if (ielem==kelem && ipoint_ielem==ipoint_kelem) bell_function=1;
			else bell_function=0;
		  }

                  element_nonlocal[ipoint_ielem*(NONLOCAL_ITEM_SIZE)+nnonlocal] = kelem;
                  element_nonlocal_ipoint[ipoint_ielem*(NONLOCAL_ITEM_SIZE)+nnonlocal] = ipoint_kelem;
                  element_nonlocal_weight[ipoint_ielem*(NONLOCAL_ITEM_SIZE)+nnonlocal] = 
			bell_function*kelem_ipoint_volume;

                  nnonlocal++;
                  if ( nnonlocal==NONLOCAL_ITEM_SIZE ) {
                    pri( "NONLOCAL_ITEM_SIZE in tochnog.h too small for nonlocal calculation." );
                    pri( "Increase it and recompile all routines." );
                    exit(TN_EXIT_STATUS);
                  }
		}
              }
	    }
           }
	  }
	 }
    	}
        db( ELEMENT_NONLOCAL, ielem, element_nonlocal, ddum, 
          length_nonloc, VERSION_NORMAL, PUT );
        db( ELEMENT_NONLOCAL_IPOINT, ielem, element_nonlocal_ipoint, ddum, 
          length_nonloc, VERSION_NORMAL, PUT );
        db( ELEMENT_NONLOCAL_WEIGHT, ielem, idum, element_nonlocal_weight, 
          length_nonloc, VERSION_NORMAL, PUT );
      }
    }
    
      // normalize the weight factors to 1
    for ( ielem=0; ielem<=max_element; ielem++ ) {
     if ( db_active_index( ELEMENT, ielem, VERSION_NORMAL )) {
        db( ELEMENT_NONLOCAL, ielem, element_nonlocal, ddum, 
          length_nonloc, VERSION_NORMAL, GET );
        db( ELEMENT_NONLOCAL_WEIGHT, ielem, idum, element_nonlocal_weight, 
          length_nonloc, VERSION_NORMAL, GET );
        db( NONLOCAL_ELEMENT_INFO, ielem, idum, nonlocal_ielem_info, 
          length_nei, VERSION_NORMAL, GET );		
        npoint_ielem=(long int)nonlocal_ielem_info[0];
        ielem_nonloc=(long int)nonlocal_ielem_info[length_nei-1];
	if(ielem_nonloc) {
          for ( ipoint_ielem=0; ipoint_ielem<npoint_ielem; ipoint_ielem++ ) {
            total_weight = 0.;
            for ( i=ipoint_ielem*NONLOCAL_ITEM_SIZE; i<(ipoint_ielem+1)*NONLOCAL_ITEM_SIZE; i++ ) 
	      if(element_nonlocal[i]>-5) //-10 if no element added
		total_weight += element_nonlocal_weight[i];
		
            assert( total_weight!=0. );
            for ( i=ipoint_ielem*NONLOCAL_ITEM_SIZE; i<(ipoint_ielem+1)*NONLOCAL_ITEM_SIZE; i++ ) 
	      if(element_nonlocal[i]>-5) //-10 if no element added
		element_nonlocal_weight[i] /= total_weight;
          }
	}
	db( ELEMENT_NONLOCAL_WEIGHT, ielem, idum, element_nonlocal_weight, 
           length_nonloc, VERSION_NORMAL, PUT );
      }
    }    
    delete[] element_nonlocal;
    delete[] element_nonlocal_ipoint;
    delete[] element_nonlocal_weight;
    delete[] nonlocal_ielem_info;
    delete[] nonlocal_kelem_info;
    if ( swit ) pri( "Out routine NONLOCAL_SET" );

    cout<<"Nonlocal weights ready"<<endl<<endl;
  }
}

void calc_nonlocal_softvar() {

/****************Averages values directly in integration points************************/
/*                    Programmed by D.Masin@city.ac.uk                                */

  long int jn=0, swit=0, idum[1],
    *element_nonlocal=NULL, max_element=0, ipoint_ielem=0, npoint_ielem=0,
    ielem_nonloc=0, jelem_nonloc=0,
    *element_nonlocal_ipoint=NULL, ielem=0, jelem=0, jelem_ipoint;
  double ddum[1], *ielem_dof=NULL, *jelem_dof=NULL,
    *element_nonlocal_weight=NULL, *nonlocal_jelem_info=NULL,
    *nonlocal_ielem_info=NULL;

  if ( materi_plasti_softvar_nonlocal && scalar_dabs(options_nonlocal_softvar)>TINY) {
    long int length_nei=1+npointmax*ndim+npointmax+2;
    long int length_nonloc=NONLOCAL_ITEM_SIZE*npointmax;
    long int mnolnuknwn=npointmax*nuknwn;
    swit = set_swit(-1,-1,"calc_nonlocal_softvar");
    if ( swit ) pri( "In routine CALC_NONLOCAL_SOFTVAR" );

    element_nonlocal = get_new_int(NONLOCAL_ITEM_SIZE*npointmax);
    element_nonlocal_ipoint = get_new_int(NONLOCAL_ITEM_SIZE*npointmax);
    element_nonlocal_weight = get_new_dbl(NONLOCAL_ITEM_SIZE*npointmax);
    ielem_dof = get_new_dbl(npointmax*nuknwn);
    jelem_dof = get_new_dbl(npointmax*nuknwn);
    nonlocal_ielem_info = get_new_dbl(length_nei);
    nonlocal_jelem_info = get_new_dbl(length_nei);

    array_set( element_nonlocal, -10, NONLOCAL_ITEM_SIZE*npointmax );
    array_set( element_nonlocal_ipoint, -10, NONLOCAL_ITEM_SIZE*npointmax );
    array_set( element_nonlocal_weight, -10, NONLOCAL_ITEM_SIZE*npointmax );
    array_set( nonlocal_ielem_info, 0, length_nei );
    array_set( nonlocal_jelem_info, 0, length_nei );
    array_set( ielem_dof, 0, npointmax*nuknwn );
    array_set( jelem_dof, 0, npointmax*nuknwn );
    db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );

    for ( ielem=0; ielem<=max_element; ielem++ ) {
     if ( db_active_index( ELEMENT, ielem, VERSION_NORMAL )) {
        db( NONLOCAL_ELEMENT_INFO, ielem, idum, nonlocal_ielem_info, 
          length_nei, VERSION_NORMAL, GET );	
	ielem_nonloc=(long int)nonlocal_ielem_info[length_nei-1];

        if(ielem_nonloc) {
          db( ELEMENT_NONLOCAL, ielem, element_nonlocal, ddum, 
            length_nonloc, VERSION_NORMAL, GET );	
          db( ELEMENT_NONLOCAL_IPOINT, ielem, element_nonlocal_ipoint, ddum, 
            length_nonloc, VERSION_NORMAL, GET );
          db( ELEMENT_NONLOCAL_WEIGHT, ielem, idum, element_nonlocal_weight, 
            length_nonloc, VERSION_NORMAL, GET );
          db( ELEMENT_DOF, ielem, idum, ielem_dof, mnolnuknwn, VERSION_NEW, GET );		
          npoint_ielem=(long int)nonlocal_ielem_info[0];	
          for ( ipoint_ielem=0; ipoint_ielem<npoint_ielem; ipoint_ielem++ ) {
            ielem_dof[nuknwn*ipoint_ielem+svnonloc_indx] = 0.;
            for ( jn=ipoint_ielem*NONLOCAL_ITEM_SIZE; jn<(ipoint_ielem+1)*NONLOCAL_ITEM_SIZE; jn++ ) {
               jelem = element_nonlocal[jn];
	       if(jelem>-5) { //jelem==-10 -- no record
                 db( NONLOCAL_ELEMENT_INFO, jelem, idum, nonlocal_jelem_info, 
                   length_nei, VERSION_NORMAL, GET );	

   	         jelem_nonloc=(long int)nonlocal_jelem_info[1+npointmax*ndim+npointmax];

	         if(jelem_nonloc) {// not used if sv_loc not filled->elem. delete, non-nonlocal groups...
   	            jelem_ipoint = element_nonlocal_ipoint[jn];
        	    db( ELEMENT_DOF, jelem, idum, jelem_dof, mnolnuknwn, VERSION_NEW, GET );		
	            ielem_dof[nuknwn*ipoint_ielem+svnonloc_indx] += 
		        element_nonlocal_weight[jn] * jelem_dof[nuknwn*jelem_ipoint+svloc_indx];  
		 }
	       }
	    }	 
	  }
	}
        db( ELEMENT_DOF, ielem, idum, ielem_dof, mnolnuknwn, VERSION_NEW, PUT );		
      }
    }
    delete[] element_nonlocal;
    delete[] element_nonlocal_ipoint;
    delete[] element_nonlocal_weight;
    delete[] nonlocal_ielem_info;
    delete[] ielem_dof;
    delete[] jelem_dof;

    if ( swit ) pri( "Out routine NONLOCAL_APPLY" );
  }
}

void nonlocal_apply( void )

{
  long int jn=0, inod=0, jnod=0, nnonlocal=0, max_node=0,
    swit=0, idum[1], *node_nonlocal=NULL;
  double ddum[1], *node_nonlocal_weight=NULL, 
    *inode_dof=NULL, *jnode_dof=NULL;

  if ( materi_plasti_f_nonlocal ) {
    swit = set_swit(-1,-1,"nonlocal_apply");
    if ( swit ) pri( "In routine NONLOCAL_APPLY" );

    node_nonlocal = get_new_int(NONLOCAL_ITEM_SIZE);
    node_nonlocal_weight = get_new_dbl(NONLOCAL_ITEM_SIZE);
    db_max_index( NODE, max_node, VERSION_NORMAL, GET );
    for ( inod=0; inod<=max_node; inod++ ) {
      if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
        db( NODE_NONLOCAL, inod, node_nonlocal, ddum, 
          nnonlocal, VERSION_NORMAL, GET );
        db( NODE_NONLOCAL_WEIGHT, inod, idum, node_nonlocal_weight, 
          nnonlocal, VERSION_NORMAL, GET );
        inode_dof = db_dbl( NODE_DOF, inod, VERSION_NEW );
        inode_dof[fn_indx] = 0.;
        for ( jn=0; jn<nnonlocal; jn++ ) {
          jnod = node_nonlocal[jn];
          jnode_dof = db_dbl( NODE_DOF, jnod, VERSION_NEW );
          inode_dof[fn_indx] += node_nonlocal_weight[jn] * jnode_dof[f_indx];  
        }
      }
    }
    delete[] node_nonlocal;
    delete[] node_nonlocal_weight;

    if ( swit ) pri( "Out routine NONLOCAL_APPLY" );
  }
}
