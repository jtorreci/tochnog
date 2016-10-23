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

void element_middle_radius_set( )

{
  long int element=0, max_element=0, idim=0, inol=0, inod=0, 
    nnol=0, length=0, ldum=0, idum[1], el[1+MNOL], nodes[MNOL];
  double distance=0., radius=0., middle[MDIM], coord[MDIM], 
    work[MDIM], ddum[1], *node_dof=NULL;

  db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
  for ( element=0; element<=max_element; element++ ) {
    if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) ) {
      db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
      nnol = length - 1; array_move( &el[1], nodes, nnol );
      array_set( middle, 0., MDIM );
      for ( inol=0; inol<nnol; inol++ ) {
        inod = nodes[inol];
        db( NODE, inod, idum, coord, ldum, VERSION_NORMAL, GET );
        if ( materi_displacement ) {
          node_dof = db_dbl( NODE_DOF, inod, VERSION_NORMAL );
          for ( idim=0; idim<ndim; idim++ ) {
            coord[idim] += node_dof[dis_indx+idim*nder];
          }
        }
        array_add( coord, middle, middle, ndim );
      }
      array_multiply( middle, middle, (double) 1/nnol, ndim );
      db( ELEMENT_MIDDLE, element, idum, middle, ndim, VERSION_NORMAL, PUT );
      radius = 0.;
      for ( inol=0; inol<nnol; inol++ ) {
        inod = nodes[inol];
        db( NODE, inod, idum, coord, ldum, VERSION_NORMAL, GET );
        distance = array_distance( coord, middle, work, ndim );
        if ( distance>radius ) radius = distance;
      }
      length = 1;
      db( ELEMENT_RADIUS, element, idum, &radius, length, VERSION_NORMAL, PUT );
    }
  }
}

void element_volume_set( long int name, long int nodes[], long int version, 
  double &element_volume )

{

  long int inod=0, jnod=0, knod=0, lnod=0, ldum=0, idum[1];
  double coord0[MDIM], coord1[MDIM], coord2[MDIM], coord3[MDIM], work[MDIM];

  if      ( ndim==1 ) {
    inod = nodes[0];
    jnod = nodes[1];
    db( NODE, inod, idum, coord0, ldum, version, GET );
    db( NODE, jnod, idum, coord1, ldum, version, GET );
    element_volume = array_distance( coord0, coord1, work, ndim );
  }
  else if ( ndim==2 ) {
    if ( name==-TRIA3 ) {
      inod = nodes[0];
      jnod = nodes[1];
      knod = nodes[2];
    }
    else {
      assert( name==-TRIA6 );
      inod = nodes[0];
      jnod = nodes[2];
      knod = nodes[5];
    }
    db( NODE, inod, idum, coord0, ldum, version, GET );
    db( NODE, jnod, idum, coord1, ldum, version, GET );
    db( NODE, knod, idum, coord2, ldum, version, GET );
    element_volume = triangle_area( coord0, coord1, coord2 );
  }
  else {
    assert( ndim==3 );
    if ( name==-TET4 ) {
      inod = nodes[0];
      jnod = nodes[1];
      knod = nodes[2];
      lnod = nodes[3];
    }
    else {
      assert( name==-TET10 );
      inod = nodes[0];
      jnod = nodes[2];
      knod = nodes[5];
      lnod = nodes[9];
    }
    db( NODE, inod, idum, coord0, ldum, version, GET );
    db( NODE, jnod, idum, coord1, ldum, version, GET );
    db( NODE, knod, idum, coord2, ldum, version, GET );
    db( NODE, lnod, idum, coord3, ldum, version, GET );
    element_volume = tetrahedron_volume( coord0, coord1, coord2, coord3 );
  }

}

void get_element_matrix_unknowns( long int element,
  long int element_matrix_unknowns[] )

{
  long int i=0, n=0, inol=0, jnol=0, inod=0, jnod=0,
    ipuknwn=0, jpuknwn=0, length=0, nnol=0, element_group=0, ldum=0,
    el[1+MNOL], nodes[MNOL], *group_matrix_unknowns=NULL;
  double ddum[1];

  if ( db_active_index( ELEMENT_MATRIX_UNKNOWNS, element, VERSION_NORMAL ) ) 
    db( ELEMENT_MATRIX_UNKNOWNS, element, element_matrix_unknowns, ddum, ldum,
      VERSION_NORMAL, GET );
  else {
    db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
    nnol = length - 1; array_move( &el[1], nodes, nnol );
    db( ELEMENT_GROUP, element, &element_group, ddum, ldum, 
      VERSION_NORMAL, GET_IF_EXISTS );
    group_matrix_unknowns = 
      db_int( GROUP_MATRIX_UNKNOWNS, element_group, VERSION_NORMAL ); 
    length = db_len( GROUP_MATRIX_UNKNOWNS, element_group, VERSION_NORMAL );
    n = length / 4;
    for ( i=0; i<n; i++ ) {
      inol = group_matrix_unknowns[4*i+0]; inod = nodes[inol];
      ipuknwn = group_matrix_unknowns[4*i+1];
      jnol = group_matrix_unknowns[4*i+2]; jnod = nodes[jnol];
      jpuknwn = group_matrix_unknowns[4*i+3];
      element_matrix_unknowns[2*i+0] = inod*npuknwn + ipuknwn;
      element_matrix_unknowns[2*i+1] = jnod*npuknwn + jpuknwn;
    }
  }

}

char *get_new_char( long int n )

{
  char *ptr=NULL;

  if ( n<=0 ) n = 1;
  if ( !(ptr = new char[n] ) ) {
    pri( "Error: cannot allocate enough memory." );
    exit(TN_EXIT_STATUS);
  }
  return ptr;
}

double *get_new_dbl( long int n )

{
  double *ptr=NULL;

  if ( n<=0 ) n = 1;
  if ( !(ptr = new double[n] ) ) {
    pri( "Error: cannot allocate enough memory." );
    exit(TN_EXIT_STATUS);
  }
  return ptr;
}

long int *get_new_int( long int n )

{
  long int *ptr=NULL;

  if ( n<=0 ) n = 1;
  if ( !(ptr = new long int[n] ) ) {
    pri( "Error: cannot allocate enough memory." );
    exit(TN_EXIT_STATUS);
  }
  return ptr;
}

int *get_new_int_short( long int n )
 
{
  int *ptr=NULL;
 
  if ( n<=0 ) n = 1;
  if ( !(ptr = new int[n] ) ) {
    pri( "Error: cannot allocate enough memory." );
    exit(TN_EXIT_STATUS);
  }
  return ptr;
}
                   

void set_environment( void )

{
  long int length=0, options_processors=1,
    options_solver=-MATRIX_ITERATIVE_BICG;
  double ddum[1];
  char *str=NULL;

  if ( !db_active_index( OPTIONS_PROCESSORS, 0, VERSION_NORMAL ) ) {
    str = getenv("TOCHNOG_OPTIONS_PROCESSORS");
    if ( str!=NULL ) {
      options_processors = atoi( str );
      length = 1;
      db( OPTIONS_PROCESSORS, 0, &options_processors, ddum, 
        length, VERSION_NORMAL, PUT );
    }
  }
  if ( !db_active_index( OPTIONS_SOLVER, 0, VERSION_NORMAL ) ) {
    str = getenv("TOCHNOG_OPTIONS_SOLVER");
    if ( str!=NULL ) {
      options_solver = -db_number( str );
      length = 1;
      db( OPTIONS_SOLVER, 0, &options_solver, ddum, 
        length, VERSION_NORMAL, PUT );
    }
  }

}

long int set_swit( long int element, long int inod, const char* routine )

{
  long int i=0, result=0, iteration=0, icontrol=0, ldum=0;
  double ddum[1];
  char *str=NULL;

  if ( strcmp(routine,swit_routine_stack[0]) ) {
    for ( i=MSTACK-2; i>=0; i-- ) {
      strcpy( swit_routine_stack[i+1], swit_routine_stack[i] );
    }
    strcpy( swit_routine_stack[0], routine );
  }
  if ( element>=0 )
    swit_element_stack = element;
  else
    swit_element_stack = -1;
  if ( inod>=0 ) 
    swit_node_stack = inod;
  else
    swit_node_stack = -1;

  if ( getenv("TOCHNOG_REPORT_ROUTINE")!=NULL ) {
    pri( "In routine ", routine );
  }

  if ( getenv("TOCHNOG_DEBUG")!=NULL ) {
    str = getenv("TOCHNOG_DEBUG");
    if ( !strcmp("yes",str ) ) result = 1;
  }
  if ( result==0 ) return 0;

  if ( element>=0 ) {
    if ( getenv("TOCHNOG_ELEMENT")!=NULL ) {
      str = getenv("TOCHNOG_ELEMENT");
      if ( element!=atol(str) ) result=0;
    }
  }

  if ( db_active_index( ICONTROL, 0, VERSION_NORMAL ) ) {
    if ( getenv("TOCHNOG_ICONTROL")!=NULL ) {
      db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
      str = getenv("TOCHNOG_ICONTROL");
      if ( icontrol!=atol(str) ) result=0;
    }
  }

  if ( db_active_index( NUMBER_ITERATIONS, 0, VERSION_NORMAL ) ) {
    if ( getenv("TOCHNOG_ITERATION")!=NULL ) {
      db( NUMBER_ITERATIONS, 0, &iteration, ddum, ldum, VERSION_NORMAL, GET );
      str = getenv("TOCHNOG_ITERATION");
      if ( iteration!=atol(str) ) result=0;
    }
  }

  if ( inod>=0 ) {
    if ( getenv("TOCHNOG_NODE")!=NULL ) {
      str = getenv("TOCHNOG_NODE");
      if ( inod!=atol(str) ) result=0;
    }
  }

  if ( strcmp(routine,"") ) {
    if ( getenv("TOCHNOG_ROUTINE")!=NULL ) {
      str = getenv("TOCHNOG_ROUTINE");
      if ( strcmp(routine,str) ) result=0;
    }
  }

  return result;

}

long int stress_indx( long int idim, long int jdim )

{
  long int kdim=0, ldim=0, indx=0;

  if ( idim<jdim ) {
    kdim = idim;
    ldim = jdim;
  }
  else {
    kdim = jdim;
    ldim = idim;
  }

  if      ( kdim==0 ) indx = 0 + ldim;
  else if ( kdim==1 ) indx = 2 + ldim;
  else if ( kdim==2 ) indx = 3 + ldim;

  return indx;
}

char *long_to_a( long int n, char s[] )

{

  int i=0, sign=0, m=0;
  char *ptr=NULL;

  m = (int) n;

  if ( ( sign = m ) < 0 ) 
    m = -m;

  i = 0;
  do {
    s[i++] = m % 10 + '0';
  } while ( ( m/= 10 ) > 0 );

  if ( sign<0 ) 
    s[i++] = '-';
  s[i] = '\0';
  string_reverse( s );

  ptr = s;
  return ptr;
  
}

void exit_tn( long int print_database_type )

{

  long int itarget=0, ntarget=0, data_item_name=0, 
    data_item_index=0, number=0, correct=0, length=0, ldum=0,
    dof_label[MUKNWN], *target_item=NULL, *int_data=NULL;
  double value=0., tolerance=0., ddum[1], *dbl_data=NULL, *target_value=NULL;

  if ( any_runtime ) {
    pri( "*** Runtime file used in calculation. ***" );
  }

  if ( print_database_type==-RESTART ) {
    print_restart( -1 );
  }
  else {
    assert( print_database_type==-YES );
    print_database( -1, VERSION_NORMAL, -YES );
  }

  print_gid( -YES );
  if ( db_partialname_any("control_print_dx") ) print_dx( -YES );

  db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db_max_index( TARGET_ITEM, ntarget, VERSION_NORMAL, GET );
  if ( ntarget>=0 ) {
    for ( itarget=0; itarget<=ntarget; itarget++ ) {
      if ( db_active_index( TARGET_ITEM, itarget, VERSION_NORMAL ) &&
           db_active_index( TARGET_VALUE, itarget, VERSION_NORMAL ) ) {
        target_item     = db_int( TARGET_ITEM, itarget, VERSION_NORMAL );
        data_item_name  = labs(target_item[0]);
        data_item_index = target_item[1];
        if ( !db_active_index( data_item_name, data_item_index, VERSION_NORMAL ) ) {
          ofstream out( "tn.log", ios::app );
          out << "\nError in calculation with data file " << data_file << ".";
          out.close();
          exit(TN_EXIT_STATUS);
        }
        else {
          length = db_len( data_item_name, data_item_index, VERSION_NORMAL );
          if ( target_item[2]<0 ) {
            array_member(dof_label,target_item[2],nuknwn,number);
            if ( db_len(data_item_name,data_item_index,VERSION_NORMAL)==npuknwn ) 
              number /= nder;
          }
          else
            number = target_item[2];
          if ( number>length-1 ) db_error( TARGET_ITEM, itarget );
          target_value = db_dbl( TARGET_VALUE, itarget, VERSION_NORMAL );
          value        = target_value[0];
          tolerance    = target_value[1];
          correct      = 1;
          if ( db_type(data_item_name)==INTEGER ) {
            int_data = db_int( data_item_name, data_item_index, VERSION_NORMAL );
            if ( int_data[number]<(value-tolerance) ||
                 int_data[number]>(value+tolerance) ) correct = 0;
          }
          else {
            dbl_data = db_dbl( data_item_name, data_item_index, VERSION_NORMAL );
            if ( dbl_data[number]<(value-tolerance) ||
                 dbl_data[number]>(value+tolerance) ) correct = 0;
          }
          if ( !correct ) {
            ofstream out( "tn.log", ios::app );
            out << "\nError in calculation with data file " << data_file << ".";
            out << "\nTarget value for: ";
            out << "data item " << db_name(data_item_name) << " ";
            out << "with index " << data_item_index << " ";
            if ( target_item[2]<0 ) {
              out << "for " << db_name(target_item[2]) << " ";
            }
            else
              out << "and value number " << number << " ";
            out << "is " << value << ".";
            if ( db_type(data_item_name)==INTEGER )
              out << "\nThe actual value is " << int_data[number] << ".";
            else
              out << "\nThe actual value is " << dbl_data[number] << ".";
            out << "\n";
            out.close();
            exit(TN_EXIT_STATUS);
          }
        }
      }
    }
  }

  db_close();

  ofstream out( "tn.log", ios::app );
  out << "\nCalculation with data file " << data_file << " ready.\n";
  out.close();
  exit(0);

}


void exit_tn_on_error( void )

{
  long int i=0, element_group=0, iarea=0, ldum=0;
  double ddum[1];

  if ( strcmp(swit_routine_stack[0],"") ) {
    pri( "\n\nLast routines:" );
    for ( i=0; i<MSTACK; i++ ) {
      if ( strcmp( swit_routine_stack[i], "" ) ) pri( swit_routine_stack[i] );
    }
  }
  if ( swit_node_stack>=0 ) pri( "\n\nLast node:" , swit_node_stack );
  if ( swit_element_stack>=0 ) {
    pri( "\n\nLast element:" , swit_element_stack );
    element_group = 0;
    db( ELEMENT_GROUP, swit_element_stack, &element_group, 
      ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    pri( "Last element has ELEMENT_GROUP", element_group );
    if ( db_active_index( ELEMENT_GROUP_AREA_ELEMENT_GROUP, swit_element_stack, VERSION_NORMAL ) ) {
      db( ELEMENT_GROUP_AREA_ELEMENT_GROUP, swit_element_stack, &iarea, 
        ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
      pri( "Last element got ELEMENT_GROUP from AREA_ELEMENT_GROUP with index", iarea );
    }
    if ( db_active_index( ELEMENT_GROUP_AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP, 
        swit_element_stack, VERSION_NORMAL ) ) {
      db( ELEMENT_GROUP_AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP, swit_element_stack, &iarea, 
        ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
      pri( "Last element got ELEMENT_GROUP from AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP with index", iarea );
    }
  }

  parallel_sys_lock();
  parallel_active = 0;
  print_database( -1, VERSION_NORMAL, -EVERYTHING );
  print_gid( -YES );
  cout << flush;
  exit(TN_EXIT_STATUS);
  parallel_sys_unlock();
}
