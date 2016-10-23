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

/* Modified on April 1st 2011 by Fernando Lorenzo to get the Von Mises Stresses 
 in the .res file
*/

#include "tochnog.h" 

double   *dbl_data[MDAT][MVERSION]; // pointer to actual double data
long int *int_data[MDAT][MVERSION]; // pointers to actual integer data
long int max_index[MDAT][MVERSION]; // maximum index allocated
long int no_index[MDAT];       // 0: index in datafile, 1: no index in data file
long int type[MDAT];           // INTEGER or DOUBLE_PRECISION
char     name[MDAT][MCHAR];    // "element", "node", ...
long int external[MDAT];       // 0: internal in TOCHNOG only, 1: also for datafile
long int data_length[MDAT];    // length (maximum for actual length)
long int fixed_length[MDAT];   // 0: all records same length; 1: records diff. length
long int data_class[MDAT];     // ELEMENT or NODE or so
long int data_required[MDAT];  // data required for this data (for the same index)
long int print_only[MDAT];     // 0: for reading and printing; 1: for printing only
long int version_all[MDAT];    // 1: with versions, 0: without versions

void db_initialize( long int dof_type[], long int dof_label[] )

{
  long int iversion=0, idim=0, ipuknwn=0, iuknwn=0, idat=0, n=0, m=0;
  char basename[MCHAR], str[MCHAR], tmpname[MCHAR];

    // fill data base administration with defaults 
  for ( idat=0; idat<MDAT; idat++ ) {
    strcpy( name[idat], " " );
    for ( iversion=0; iversion<MVERSION; iversion++ )
      max_index[idat][iversion] = -1;
  }
  array_set( version_all, 0, MDAT );
  array_set( data_class, -1, MDAT );
  array_set( data_required, -1, MDAT );
  array_set( data_length, 0, MDAT );
  array_set( external, 1, MDAT );
  array_set( fixed_length, 1, MDAT );
  array_set( no_index, 0, MDAT );
  array_set( type, 0, MDAT );
  array_set( print_only, 0, MDAT );

  strcpy(name[ABOVE],"above" );

  strcpy(name[ABSOL],"absol" );

  strcpy(name[ADD],"add" );

  strcpy(name[ADD_ALWAYS],"add_always" );

  strcpy(name[ALL],"all" );

  strcpy(name[AREA],"area");

  strcpy(name[AREA_ELEMENT_GROUP],"area_element_group");
  type[AREA_ELEMENT_GROUP] = INTEGER;
  data_length[AREA_ELEMENT_GROUP] = 3;
  data_class[AREA_ELEMENT_GROUP] = AREA;

  strcpy(name[AREA_ELEMENT_GROUP_SEQUENCE],"area_element_group_sequence");
  type[AREA_ELEMENT_GROUP_SEQUENCE] = INTEGER;
  data_length[AREA_ELEMENT_GROUP_SEQUENCE] = 1000;
  fixed_length[AREA_ELEMENT_GROUP_SEQUENCE] = 0;
  data_class[AREA_ELEMENT_GROUP_SEQUENCE] = AREA;
  data_required[AREA_ELEMENT_GROUP_SEQUENCE] = AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP;

  strcpy(name[AREA_ELEMENT_GROUP_SEQUENCE_ELEMENT],"area_element_group_sequence_element");
  type[AREA_ELEMENT_GROUP_SEQUENCE_ELEMENT] = INTEGER;
  data_length[AREA_ELEMENT_GROUP_SEQUENCE_ELEMENT] = 1;
  data_class[AREA_ELEMENT_GROUP_SEQUENCE_ELEMENT] = AREA;
  data_required[AREA_ELEMENT_GROUP_SEQUENCE_ELEMENT] = AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP;

  strcpy(name[AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP],"area_element_group_sequence_elementgroup");
  type[AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP] = INTEGER;
  data_length[AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP] = DATA_ITEM_SIZE;
  data_class[AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP] = AREA;
  fixed_length[AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP] = 0;

  strcpy(name[AREA_ELEMENT_GROUP_SEQUENCE_GEOMETRY],"area_element_group_sequence_geometry");
  type[AREA_ELEMENT_GROUP_SEQUENCE_GEOMETRY] = INTEGER;
  data_length[AREA_ELEMENT_GROUP_SEQUENCE_GEOMETRY] = 2;
  data_class[AREA_ELEMENT_GROUP_SEQUENCE_GEOMETRY] = AREA;
  data_required[AREA_ELEMENT_GROUP_SEQUENCE_GEOMETRY] = AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP;

  strcpy(name[AREA_ELEMENT_GROUP_SEQUENCE_TIME],"area_element_group_sequence_time");
  type[AREA_ELEMENT_GROUP_SEQUENCE_TIME] = DOUBLE_PRECISION;
  data_length[AREA_ELEMENT_GROUP_SEQUENCE_TIME] = DATA_ITEM_SIZE;
  data_class[AREA_ELEMENT_GROUP_SEQUENCE_TIME] = AREA;
  fixed_length[AREA_ELEMENT_GROUP_SEQUENCE_TIME] = 0;
  data_required[AREA_ELEMENT_GROUP_SEQUENCE_GEOMETRY] = AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP;

  strcpy(name[AREA_NODE_DATAITEM],"area_node_dataitem");
  type[AREA_NODE_DATAITEM] = INTEGER;
  data_length[AREA_NODE_DATAITEM] = 3;
  data_class[AREA_NODE_DATAITEM] = AREA;

  strcpy(name[AREA_NODE_DATAITEM_DOUBLE],"area_node_dataitem_double");
  type[AREA_NODE_DATAITEM_DOUBLE] = DOUBLE_PRECISION;
  data_length[AREA_NODE_DATAITEM_DOUBLE] = DATA_ITEM_SIZE;
  data_class[AREA_NODE_DATAITEM_DOUBLE] = AREA;
  fixed_length[AREA_NODE_DATAITEM_DOUBLE] = 0;
  data_required[AREA_NODE_DATAITEM_DOUBLE] = AREA_NODE_DATAITEM;

  strcpy(name[AREA_NODE_DATAITEM_INTEGER],"area_node_dataitem_integer");
  type[AREA_NODE_DATAITEM_INTEGER] = INTEGER;
  data_length[AREA_NODE_DATAITEM_INTEGER] = DATA_ITEM_SIZE;
  data_class[AREA_NODE_DATAITEM_INTEGER] = AREA;
  fixed_length[AREA_NODE_DATAITEM_INTEGER] = 0;
  data_required[AREA_NODE_DATAITEM_INTEGER] = AREA_NODE_DATAITEM;

  strcpy(name[ASM],"asm");

  strcpy(name[AVERAGE],"average");

  strcpy(name[BAR],"bar" );

  strcpy(name[BAR2],"bar2");

  strcpy(name[BAR3],"bar3");

  strcpy(name[BAR4],"bar4");

  strcpy(name[BCGS],"bcgs");

  strcpy(name[BEAM],"beam");

  strcpy(name[BEAM_ROTATION],"beam_rotation");

  strcpy(name[BELOW],"below" );

  strcpy(name[BICG],"bicg");

  strcpy(name[BJACOBI],"bjacobi");

  strcpy(name[BOUNDA],"bounda");

  strcpy(name[BOUNDA_FORCE],"bounda_force");
  type[BOUNDA_FORCE] = INTEGER;
  data_length[BOUNDA_FORCE] = MBOUNDA;
  fixed_length[BOUNDA_FORCE] = 0;
  data_class[BOUNDA_FORCE] = BOUNDA;

  strcpy(name[BOUNDA_SINE],"bounda_sine");
  type[BOUNDA_SINE] = DOUBLE_PRECISION;
  data_length[BOUNDA_SINE] = DATA_ITEM_SIZE;
  fixed_length[BOUNDA_SINE] = 0;
  data_class[BOUNDA_SINE] = BOUNDA;

  strcpy(name[BOUNDA_TIME],"bounda_time");
  type[BOUNDA_TIME] = DOUBLE_PRECISION;
  data_length[BOUNDA_TIME] = DATA_ITEM_SIZE;
  fixed_length[BOUNDA_TIME] = 0;
  data_class[BOUNDA_TIME] = BOUNDA;

  strcpy(name[BOUNDA_TIME_FILE],"bounda_time_file");
  type[BOUNDA_TIME_FILE] = INTEGER;
  data_length[BOUNDA_TIME_FILE] = 1;
  data_class[BOUNDA_TIME_FILE] = BOUNDA;

  strcpy(name[BOUNDA_TIME_USER],"bounda_time_user");
  type[BOUNDA_TIME_USER] = INTEGER;
  data_length[BOUNDA_TIME_USER] = 1;
  data_class[BOUNDA_TIME_USER] = BOUNDA;

  strcpy(name[BOUNDA_UNKNOWN],"bounda_unknown");
  type[BOUNDA_UNKNOWN] = INTEGER;
  data_length[BOUNDA_UNKNOWN] = MBOUNDA;
  fixed_length[BOUNDA_UNKNOWN] = 0;
  data_class[BOUNDA_UNKNOWN] = BOUNDA;

  strcpy(name[BRICK],"brick" );

  strcpy(name[CALCULATE_STRESSINTENSITYFACTOR],"calculate_stressintensityfactor");

  strcpy(name[CG],"cg");

  strcpy(name[CGS],"cgs");

  strcpy(name[CHANGE],"change");

  strcpy(name[CHANGE_DATAITEM],"change_dataitem");
  type[CHANGE_DATAITEM] = INTEGER;
  data_length[CHANGE_DATAITEM] = 4;

  strcpy(name[CHANGE_DATAITEM_TIME],"change_dataitem_time");
  type[CHANGE_DATAITEM_TIME] = DOUBLE_PRECISION;
  data_length[CHANGE_DATAITEM_TIME] = DATA_ITEM_SIZE;
  fixed_length[CHANGE_DATAITEM_TIME] = 0;
  data_required[CHANGE_DATAITEM_TIME] = CHANGE_DATAITEM;

  strcpy(name[CHANGE_DATAITEM_TIME_DISCRETE],"change_dataitem_time_discrete");
  type[CHANGE_DATAITEM_TIME_DISCRETE] = INTEGER;
  data_length[CHANGE_DATAITEM_TIME_DISCRETE] = 1;
  data_required[CHANGE_DATAITEM_TIME_DISCRETE] = CHANGE_DATAITEM;

  strcpy(name[CHANGE_DATAITEM_TIME_USER],"change_dataitem_time_user");
  type[CHANGE_DATAITEM_TIME_USER] = INTEGER;
  data_length[CHANGE_DATAITEM_TIME_USER] = 1;

  strcpy(name[CHANGE_GEOMETRY],"change_geometry");
  type[CHANGE_GEOMETRY] = INTEGER;
  data_length[CHANGE_GEOMETRY] = 3;

  strcpy(name[CHANGE_GEOMETRY_TIME_USER],"change_geometry_time_user");
  type[CHANGE_GEOMETRY_TIME_USER] = INTEGER;
  data_length[CHANGE_GEOMETRY_TIME_USER] = 1;

  strcpy(name[CHEBYCHEV],"chebychev");

  strcpy(name[CHECK],"check" );

  strcpy(name[CHECK_COMBINATION],"check_combination");
  type[CHECK_COMBINATION] = INTEGER;
  data_length[CHECK_COMBINATION] = 1;
  no_index[CHECK_COMBINATION] = 1;
  data_class[CHECK_COMBINATION] = CHECK_COMBINATION;

  strcpy(name[CHECK_INDEX],"check_index" );

  strcpy(name[CHECK_NUMBER],"check_number" );

  strcpy(name[CIRCLE],"circle" );

  strcpy(name[CIRCLE_HOLLOW],"circle_hollow" );

  strcpy(name[COMPOSITE],"composite");

  strcpy(name[CONDIF],"condif");

  strcpy(name[CONDIF_CONVECTION],"condif_convection");
  type[CONDIF_CONVECTION] = DOUBLE_PRECISION;
  data_length[CONDIF_CONVECTION] = 2;
  data_class[CONDIF_CONVECTION] = CONDIF;

  strcpy(name[CONDIF_CONVECTION_GEOMETRY],"condif_convection_geometry");
  type[CONDIF_CONVECTION_GEOMETRY] = INTEGER;
  data_length[CONDIF_CONVECTION_GEOMETRY] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_CONVECTION_GEOMETRY] = 0;
  data_class[CONDIF_CONVECTION_GEOMETRY] = CONDIF;
  data_required[CONDIF_CONVECTION_GEOMETRY] = CONDIF_CONVECTION;

  strcpy(name[CONDIF_RADIATION],"condif_radiation");
  type[CONDIF_RADIATION] = DOUBLE_PRECISION;
  data_length[CONDIF_RADIATION] = 2;
  data_class[CONDIF_RADIATION] = CONDIF;

  strcpy(name[CONDIF_RADIATION_GEOMETRY],"condif_radiation_geometry");
  type[CONDIF_RADIATION_GEOMETRY] = INTEGER;
  data_length[CONDIF_RADIATION_GEOMETRY] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_RADIATION_GEOMETRY] = 0;
  data_class[CONDIF_RADIATION_GEOMETRY] = CONDIF;
  data_required[CONDIF_RADIATION_GEOMETRY] = CONDIF_RADIATION;

  strcpy(name[CONDIF_TEMPERATURE],"condif_temperature");

  strcpy(name[CONTACT],"contact");

  strcpy(name[CONTACTSPRING],"contactspring");

  strcpy(name[CONTACT_FRICTION],"contact_friction");
  type[CONTACT_FRICTION] = DOUBLE_PRECISION;
  data_length[CONTACT_FRICTION] = 1;
  no_index[CONTACT_FRICTION] = 1;
  data_class[CONTACT_FRICTION] = CONTACT;

  strcpy(name[CONTACT_GEOMETRY],"contact_geometry");
  type[CONTACT_GEOMETRY] = INTEGER;
  data_length[CONTACT_GEOMETRY] = 2;
  data_class[CONTACT_GEOMETRY] = CONTACT;

  strcpy(name[CONTACT_GEOMETRY_SWITCH],"contact_geometry_switch");
  type[CONTACT_GEOMETRY_SWITCH] = INTEGER;
  data_length[CONTACT_GEOMETRY_SWITCH] = 1;
  data_class[CONTACT_GEOMETRY_SWITCH] = CONTACT;
  data_required[CONTACT_GEOMETRY_SWITCH] = CONTACT_GEOMETRY;

  strcpy(name[CONTACT_HEATGENERATION],"contact_heatgeneration");
  type[CONTACT_HEATGENERATION] = DOUBLE_PRECISION;
  data_length[CONTACT_HEATGENERATION] = 1;
  no_index[CONTACT_HEATGENERATION] = 1;
  data_class[CONTACT_HEATGENERATION] = CONTACT;

  strcpy(name[CONTACT_PENALTY_PRESSURE],"contact_penalty_pressure");
  type[CONTACT_PENALTY_PRESSURE] = DOUBLE_PRECISION;
  data_length[CONTACT_PENALTY_PRESSURE] = 1;
  no_index[CONTACT_PENALTY_PRESSURE] = 1;
  data_class[CONTACT_PENALTY_PRESSURE] = CONTACT;

  strcpy(name[CONTACT_PENALTY_TEMPERATURE],"contact_penalty_temperature");
  type[CONTACT_PENALTY_TEMPERATURE] = DOUBLE_PRECISION;
  data_length[CONTACT_PENALTY_TEMPERATURE] = 1;
  no_index[CONTACT_PENALTY_TEMPERATURE] = 1;
  data_class[CONTACT_PENALTY_TEMPERATURE] = CONTACT;

  strcpy(name[CONTACT_PENALTY_VELOCITY],"contact_penalty_velocity");
  type[CONTACT_PENALTY_VELOCITY] = DOUBLE_PRECISION;
  data_length[CONTACT_PENALTY_VELOCITY] = 1;
  no_index[CONTACT_PENALTY_VELOCITY] = 1;
  data_class[CONTACT_PENALTY_VELOCITY] = CONTACT;

  strcpy(name[CONTACT_RELAXATION],"contact_relaxation");
  type[CONTACT_RELAXATION] = DOUBLE_PRECISION;
  data_length[CONTACT_RELAXATION] = 1;
  no_index[CONTACT_RELAXATION] = 1;
  data_class[CONTACT_RELAXATION] = CONTACT;

  strcpy(name[CONTACT_STICK],"contact_stick");
  type[CONTACT_STICK] = INTEGER;
  data_length[CONTACT_STICK] = 1;
  no_index[CONTACT_STICK] = 1;
  data_class[CONTACT_STICK] = CONTACT;

  strcpy(name[CONTROL_CRACK],"control_crack");
  type[CONTROL_CRACK] = INTEGER;
  data_length[CONTROL_CRACK] = 1;
  data_class[CONTROL_CRACK] = CONTROL;     

  strcpy(name[CONTROL_DATA_DELETE],"control_data_delete");
  type[CONTROL_DATA_DELETE] = INTEGER;
  data_length[CONTROL_DATA_DELETE] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_DATA_DELETE] = 0;
  data_class[CONTROL_DATA_DELETE] = CONTROL;

  strcpy(name[CONTROL_DATA_INITELDOF_GEOMETRY],"control_data_initeldof_geometry");
  type[CONTROL_DATA_INITELDOF_GEOMETRY] = INTEGER;
  data_length[CONTROL_DATA_INITELDOF_GEOMETRY] = 2;
  data_class[CONTROL_DATA_INITELDOF_GEOMETRY] = CONTROL;

  strcpy(name[CONTROL_DATA_PUT],"control_data_put");
  type[CONTROL_DATA_PUT] = INTEGER;
  data_length[CONTROL_DATA_PUT] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_DATA_PUT] = 0;
  data_class[CONTROL_DATA_PUT] = CONTROL;

  strcpy(name[CONTROL_DATA_PUT_DOUBLE],"control_data_put_double");
  type[CONTROL_DATA_PUT_DOUBLE] = DOUBLE_PRECISION;
  data_length[CONTROL_DATA_PUT_DOUBLE] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_DATA_PUT_DOUBLE] = 0;
  data_class[CONTROL_DATA_PUT_DOUBLE] = CONTROL;
  data_required[CONTROL_DATA_PUT_DOUBLE] = CONTROL_DATA_PUT;

  strcpy(name[CONTROL_DATA_PUT_DOUBLE_NODE],"control_data_put_double_node");
  type[CONTROL_DATA_PUT_DOUBLE_NODE] = DOUBLE_PRECISION;
  data_length[CONTROL_DATA_PUT_DOUBLE_NODE] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_DATA_PUT_DOUBLE_NODE] = 0;
  data_class[CONTROL_DATA_PUT_DOUBLE_NODE] = CONTROL;
  data_required[CONTROL_DATA_PUT_DOUBLE_NODE] = CONTROL_DATA_PUT;

  strcpy(name[CONTROL_DATA_PUT_INTEGER],"control_data_put_integer");
  type[CONTROL_DATA_PUT_INTEGER] = INTEGER;
  data_length[CONTROL_DATA_PUT_INTEGER] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_DATA_PUT_INTEGER] = 0;
  data_class[CONTROL_DATA_PUT_INTEGER] = CONTROL;
  data_required[CONTROL_DATA_PUT_INTEGER] = CONTROL_DATA_PUT;

  strcpy(name[CONTROL_DISTRIBUTE],"control_distribute");
  type[CONTROL_DISTRIBUTE] = INTEGER;
  data_length[CONTROL_DISTRIBUTE] = DATA_ITEM_SIZE;
  data_class[CONTROL_DISTRIBUTE] = CONTROL;
  fixed_length[CONTROL_DISTRIBUTE] = 0;

  strcpy(name[CONTROL_DISTRIBUTE_VALUES],"control_distribute_values");
  type[CONTROL_DISTRIBUTE_VALUES] = DOUBLE_PRECISION;
  data_length[CONTROL_DISTRIBUTE_VALUES] = DATA_ITEM_SIZE;
  data_class[CONTROL_DISTRIBUTE_VALUES] = CONTROL;
  fixed_length[CONTROL_DISTRIBUTE_VALUES] = 0;
  data_required[CONTROL_DISTRIBUTE_VALUES] = CONTROL_DISTRIBUTE;

  strcpy(name[CONTROL_EIGEN],"control_eigen");
  type[CONTROL_EIGEN] = INTEGER;
  data_length[CONTROL_EIGEN] = 2;
  data_class[CONTROL_EIGEN] = CONTROL;

  strcpy(name[CONTROL_EIGEN_SCALE],"control_eigen_scale");
  type[CONTROL_EIGEN_SCALE] = DOUBLE_PRECISION;
  data_length[CONTROL_EIGEN_SCALE] = 1;
  data_class[CONTROL_EIGEN_SCALE] = CONTROL;
  data_required[CONTROL_EIGEN_SCALE] = CONTROL_EIGEN;

  strcpy(name[CONTROL_EIGEN_VALUES],"control_eigen_values");
  type[CONTROL_EIGEN_VALUES] = DOUBLE_PRECISION;
  data_length[CONTROL_EIGEN_VALUES] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_EIGEN_VALUES] = 0;
  data_class[CONTROL_EIGEN_VALUES] = CONTROL;
  no_index[CONTROL_EIGEN_VALUES] = 1;

  strcpy(name[CONTROL_MATERI_DIFFUSION],"control_materi_diffusion");
  type[CONTROL_MATERI_DIFFUSION] = INTEGER;
  data_length[CONTROL_MATERI_DIFFUSION] = 1;
  data_class[CONTROL_MATERI_DIFFUSION] = CONTROL;

  strcpy(name[CONTROL_MESH_ADJUST_GEOMETRY],"control_mesh_adjust_geometry");
  type[CONTROL_MESH_ADJUST_GEOMETRY] = INTEGER;
  data_length[CONTROL_MESH_ADJUST_GEOMETRY] = 4;
  data_class[CONTROL_MESH_ADJUST_GEOMETRY] = CONTROL;

  strcpy(name[CONTROL_MESH_DELETE_GEOMETRY],"control_mesh_delete_geometry");
  type[CONTROL_MESH_DELETE_GEOMETRY] = INTEGER;
  data_length[CONTROL_MESH_DELETE_GEOMETRY] = 2;
  data_class[CONTROL_MESH_DELETE_GEOMETRY] = CONTROL;

  strcpy(name[CONTROL_MESH_DELETE_GEOMETRY_ELEMENT],"control_mesh_delete_geometry_element");
  type[CONTROL_MESH_DELETE_GEOMETRY_ELEMENT] = INTEGER;
  data_length[CONTROL_MESH_DELETE_GEOMETRY_ELEMENT] = DATA_ITEM_SIZE;
  data_class[CONTROL_MESH_DELETE_GEOMETRY_ELEMENT] = CONTROL;     
  fixed_length[CONTROL_MESH_DELETE_GEOMETRY_ELEMENT] = 0;     
  data_required[CONTROL_MESH_DELETE_GEOMETRY_ELEMENT] = CONTROL_MESH_DELETE_GEOMETRY;     

  strcpy(name[CONTROL_MESH_DELETE_GEOMETRY_ELEMENTGROUP],"control_mesh_delete_geometry_elementgroup");
  type[CONTROL_MESH_DELETE_GEOMETRY_ELEMENTGROUP] = INTEGER;
  data_length[CONTROL_MESH_DELETE_GEOMETRY_ELEMENTGROUP] = DATA_ITEM_SIZE;
  data_class[CONTROL_MESH_DELETE_GEOMETRY_ELEMENTGROUP] = CONTROL;     
  fixed_length[CONTROL_MESH_DELETE_GEOMETRY_ELEMENTGROUP] = 0;     
  data_required[CONTROL_MESH_DELETE_GEOMETRY_ELEMENTGROUP] = CONTROL_MESH_DELETE_GEOMETRY;     

  strcpy(name[CONTROL_MESH_DELETE_GEOMETRY_FACTOR],"control_mesh_delete_geometry_factor");
  type[CONTROL_MESH_DELETE_GEOMETRY_FACTOR] = DOUBLE_PRECISION;
  data_length[CONTROL_MESH_DELETE_GEOMETRY_FACTOR] = 2;
  data_class[CONTROL_MESH_DELETE_GEOMETRY_FACTOR] = CONTROL;     
  fixed_length[CONTROL_MESH_DELETE_GEOMETRY_FACTOR] = 0;     
  data_required[CONTROL_MESH_DELETE_GEOMETRY_FACTOR] = CONTROL_MESH_DELETE_GEOMETRY;     

  strcpy(name[CONTROL_MESH_DELETE_GEOMETRY_MOVENODES],"control_mesh_delete_geometry_movenodes");
  type[CONTROL_MESH_DELETE_GEOMETRY_MOVENODES] = INTEGER;
  data_length[CONTROL_MESH_DELETE_GEOMETRY_MOVENODES] = 1;
  data_class[CONTROL_MESH_DELETE_GEOMETRY_MOVENODES] = CONTROL;     
  data_required[CONTROL_MESH_DELETE_GEOMETRY_MOVENODES] = CONTROL_MESH_DELETE_GEOMETRY;     

  strcpy(name[CONTROL_MESH_DELETE_SMALL],"control_mesh_delete_small");
  type[CONTROL_MESH_DELETE_SMALL] = DOUBLE_PRECISION;
  data_length[CONTROL_MESH_DELETE_SMALL] = 1;
  data_class[CONTROL_MESH_DELETE_SMALL] = CONTROL;

  strcpy(name[CONTROL_MESH_EXTRUDE],"control_mesh_extrude");
  type[CONTROL_MESH_EXTRUDE] = DOUBLE_PRECISION;
  data_length[CONTROL_MESH_EXTRUDE] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_MESH_EXTRUDE] = 0;
  data_class[CONTROL_MESH_EXTRUDE] = CONTROL;

  strcpy(name[CONTROL_MESH_EXTRUDE_N],"control_mesh_extrude_n");
  type[CONTROL_MESH_EXTRUDE_N] = INTEGER;
  data_length[CONTROL_MESH_EXTRUDE_N] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_MESH_EXTRUDE_N] = 0;
  data_class[CONTROL_MESH_EXTRUDE_N] = CONTROL;
  data_required[CONTROL_MESH_EXTRUDE_N] = CONTROL_MESH_EXTRUDE;

  strcpy(name[CONTROL_MESH_GENERATE_BEAM],"control_mesh_generate_beam");
  type[CONTROL_MESH_GENERATE_BEAM] = INTEGER;
  data_length[CONTROL_MESH_GENERATE_BEAM] = 3;
  data_class[CONTROL_MESH_GENERATE_BEAM] = CONTROL;

  strcpy(name[CONTROL_MESH_GENERATE_CONTACTSPRING],"control_mesh_generate_contactspring");
  type[CONTROL_MESH_GENERATE_CONTACTSPRING] = INTEGER;
  data_length[CONTROL_MESH_GENERATE_CONTACTSPRING] = 3;
  data_class[CONTROL_MESH_GENERATE_CONTACTSPRING] = CONTROL;

  strcpy(name[CONTROL_MESH_GENERATE_CONTACTSPRING_ELEMENT],"control_mesh_generate_contactspring_element");
  type[CONTROL_MESH_GENERATE_CONTACTSPRING_ELEMENT] = INTEGER;
  data_length[CONTROL_MESH_GENERATE_CONTACTSPRING_ELEMENT] = 2;
  data_class[CONTROL_MESH_GENERATE_CONTACTSPRING_ELEMENT] = CONTROL;

  strcpy(name[CONTROL_MESH_GENERATE_SPRING1],"control_mesh_generate_spring1");
  type[CONTROL_MESH_GENERATE_SPRING1] = INTEGER;
  data_length[CONTROL_MESH_GENERATE_SPRING1] = 3;
  data_class[CONTROL_MESH_GENERATE_SPRING1] = CONTROL;

  strcpy(name[CONTROL_MESH_GENERATE_SPRING2],"control_mesh_generate_spring2");
  type[CONTROL_MESH_GENERATE_SPRING2] = INTEGER;
  data_length[CONTROL_MESH_GENERATE_SPRING2] = 3;
  data_class[CONTROL_MESH_GENERATE_SPRING2] = CONTROL;

  strcpy(name[CONTROL_MESH_GENERATE_TRUSS],"control_mesh_generate_truss");
  type[CONTROL_MESH_GENERATE_TRUSS] = INTEGER;
  data_length[CONTROL_MESH_GENERATE_TRUSS] = 3;
  data_class[CONTROL_MESH_GENERATE_TRUSS] = CONTROL;

  strcpy(name[CONTROL_MESH_GENERATE_TRUSSBEAM],"control_mesh_generate_trussbeam");
  type[CONTROL_MESH_GENERATE_TRUSSBEAM] = INTEGER;
  data_length[CONTROL_MESH_GENERATE_TRUSSBEAM] = 3;
  data_class[CONTROL_MESH_GENERATE_TRUSSBEAM] = CONTROL;

  strcpy(name[CONTROL_MESH_GENERATE_TRUSS_BEAM_LOOSE],"control_mesh_generate_truss_beam_loose");
  type[CONTROL_MESH_GENERATE_TRUSS_BEAM_LOOSE] = INTEGER;
  data_length[CONTROL_MESH_GENERATE_TRUSS_BEAM_LOOSE] = 1;
  data_class[CONTROL_MESH_GENERATE_TRUSS_BEAM_LOOSE] = CONTROL;

  strcpy(name[CONTROL_MESH_GENERATE_TRUSS_BEAM_MACRO],"control_mesh_generate_truss_beam_macro");
  type[CONTROL_MESH_GENERATE_TRUSS_BEAM_MACRO] = INTEGER;
  data_length[CONTROL_MESH_GENERATE_TRUSS_BEAM_MACRO] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_MESH_GENERATE_TRUSS_BEAM_MACRO] = 0;
  data_class[CONTROL_MESH_GENERATE_TRUSS_BEAM_MACRO] = CONTROL;

  strcpy(name[CONTROL_MESH_MACRO],"control_mesh_macro");
  type[CONTROL_MESH_MACRO] = INTEGER;
  data_length[CONTROL_MESH_MACRO] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_MESH_MACRO] = 0;
  data_class[CONTROL_MESH_MACRO] = CONTROL;

  strcpy(name[CONTROL_MESH_MACRO_ELEMENT],"control_mesh_macro_element");
  type[CONTROL_MESH_MACRO_ELEMENT] = INTEGER;
  data_length[CONTROL_MESH_MACRO_ELEMENT] = 1;
  data_class[CONTROL_MESH_MACRO_ELEMENT] = CONTROL;

  strcpy(name[CONTROL_MESH_MACRO_PARAMETERS],"control_mesh_macro_parameters");
  type[CONTROL_MESH_MACRO_PARAMETERS] = DOUBLE_PRECISION;
  data_length[CONTROL_MESH_MACRO_PARAMETERS] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_MESH_MACRO_PARAMETERS] = 0;
  data_class[CONTROL_MESH_MACRO_PARAMETERS] = CONTROL;
  data_required[CONTROL_MESH_MACRO_PARAMETERS] = CONTROL_MESH_MACRO;

  strcpy(name[CONTROL_MESH_MACRO_SET_NODE_BOUNDARY],"control_mesh_macro_set_node_boundary");
  type[CONTROL_MESH_MACRO_SET_NODE_BOUNDARY] = INTEGER;
  data_length[CONTROL_MESH_MACRO_SET_NODE_BOUNDARY] = 1;
  data_class[CONTROL_MESH_MACRO_SET_NODE_BOUNDARY] = CONTROL;
  data_required[CONTROL_MESH_MACRO_SET_NODE_BOUNDARY] = CONTROL_MESH_MACRO;

  strcpy(name[CONTROL_MESH_MERGE],"control_mesh_merge");
  type[CONTROL_MESH_MERGE] = INTEGER;
  data_length[CONTROL_MESH_MERGE] = 1;
  data_class[CONTROL_MESH_MERGE] = CONTROL;

  strcpy(name[CONTROL_MESH_MERGE_EPSCOORD],"control_mesh_merge_epscoord");
  type[CONTROL_MESH_MERGE_EPSCOORD] = DOUBLE_PRECISION;
  data_length[CONTROL_MESH_MERGE_EPSCOORD] = 1;
  data_class[CONTROL_MESH_MERGE_EPSCOORD] = CONTROL;

  strcpy(name[CONTROL_MESH_MERGE_MACRO_GENERATE],"control_mesh_merge_macro_generate");
  type[CONTROL_MESH_MERGE_MACRO_GENERATE] = INTEGER;
  data_length[CONTROL_MESH_MERGE_MACRO_GENERATE] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_MESH_MERGE_MACRO_GENERATE] = 0;
  data_class[CONTROL_MESH_MERGE_MACRO_GENERATE] = CONTROL;

  strcpy(name[CONTROL_MESH_MERGE_NOT],"control_mesh_merge_not");
  type[CONTROL_MESH_MERGE_NOT] = INTEGER;
  data_length[CONTROL_MESH_MERGE_NOT] = 2;
  data_class[CONTROL_MESH_MERGE_NOT] = CONTROL;

  strcpy(name[CONTROL_MESH_NEW_MESH],"control_mesh_new_mesh");
  type[CONTROL_MESH_NEW_MESH] = DOUBLE_PRECISION;
  data_length[CONTROL_MESH_NEW_MESH] = 1;
  data_class[CONTROL_MESH_NEW_MESH] = CONTROL;

  strcpy(name[CONTROL_MESH_NEW_MESH_ELEMENT],"control_mesh_new_mesh_element");
  type[CONTROL_MESH_NEW_MESH_ELEMENT] = INTEGER;
  data_length[CONTROL_MESH_NEW_MESH_ELEMENT] = 1;
  data_class[CONTROL_MESH_NEW_MESH_ELEMENT] = CONTROL;
  data_required[CONTROL_MESH_NEW_MESH_ELEMENT] = CONTROL_MESH_NEW_MESH;

  strcpy(name[CONTROL_MESH_NEW_MESH_REGION],"control_mesh_new_mesh_region");
  type[CONTROL_MESH_NEW_MESH_REGION] = INTEGER;
  data_length[CONTROL_MESH_NEW_MESH_REGION] = 6+nuknwn;
  fixed_length[CONTROL_MESH_NEW_MESH_REGION] = 0;
  data_class[CONTROL_MESH_NEW_MESH_REGION] = CONTROL;
  data_required[CONTROL_MESH_NEW_MESH_REGION] = CONTROL_MESH_NEW_MESH;

  strcpy(name[CONTROL_MESH_REFINE_GLOBALLY],"control_mesh_refine_globally");
  type[CONTROL_MESH_REFINE_GLOBALLY] = INTEGER;
  data_length[CONTROL_MESH_REFINE_GLOBALLY] = 4;
  fixed_length[CONTROL_MESH_REFINE_GLOBALLY] = 0;
  data_class[CONTROL_MESH_REFINE_GLOBALLY] = CONTROL;

  strcpy(name[CONTROL_MESH_REFINE_GLOBALLY_GEOMETRY],"control_mesh_refine_globally_geometry");
  type[CONTROL_MESH_REFINE_GLOBALLY_GEOMETRY] = INTEGER;
  data_length[CONTROL_MESH_REFINE_GLOBALLY_GEOMETRY] = 2;
  data_class[CONTROL_MESH_REFINE_GLOBALLY_GEOMETRY] = CONTROL;
  data_required[CONTROL_MESH_REFINE_GLOBALLY_GEOMETRY] = CONTROL_MESH_REFINE_GLOBALLY;

  strcpy(name[CONTROL_MESH_REFINE_LOCALLY],"control_mesh_refine_locally");
  type[CONTROL_MESH_REFINE_LOCALLY] = DOUBLE_PRECISION;
  data_length[CONTROL_MESH_REFINE_LOCALLY] = 1;
  data_class[CONTROL_MESH_REFINE_LOCALLY] = CONTROL;

  strcpy(name[CONTROL_MESH_REFINE_LOCALLY_GEOMETRY],"control_mesh_refine_locally_geometry");
  type[CONTROL_MESH_REFINE_LOCALLY_GEOMETRY] = INTEGER;
  data_length[CONTROL_MESH_REFINE_LOCALLY_GEOMETRY] = 2;
  data_class[CONTROL_MESH_REFINE_LOCALLY_GEOMETRY] = CONTROL;
  data_required[CONTROL_MESH_REFINE_LOCALLY_GEOMETRY] = CONTROL_MESH_REFINE_LOCALLY;

  strcpy(name[CONTROL_MESH_REFINE_LOCALLY_NOT],"control_mesh_refine_locally_not");
  type[CONTROL_MESH_REFINE_LOCALLY_NOT] = INTEGER;
  data_length[CONTROL_MESH_REFINE_LOCALLY_NOT] = 2;
  data_class[CONTROL_MESH_REFINE_LOCALLY_NOT] = CONTROL;
  data_required[CONTROL_MESH_REFINE_LOCALLY_NOT] = CONTROL_MESH_REFINE_LOCALLY;

  strcpy(name[CONTROL_MESH_REFINE_LOCALLY_ONLY],"control_mesh_refine_locally_only");
  type[CONTROL_MESH_REFINE_LOCALLY_ONLY] = INTEGER;
  data_length[CONTROL_MESH_REFINE_LOCALLY_ONLY] = 2;
  data_class[CONTROL_MESH_REFINE_LOCALLY_ONLY] = CONTROL;
  data_required[CONTROL_MESH_REFINE_LOCALLY_ONLY] = CONTROL_MESH_REFINE_LOCALLY;

  strcpy(name[CONTROL_MESH_REFINE_LOCALLY_UNKNOWN],"control_mesh_refine_locally_unknown");
  type[CONTROL_MESH_REFINE_LOCALLY_UNKNOWN] = INTEGER;
  data_length[CONTROL_MESH_REFINE_LOCALLY_UNKNOWN] = 1;
  data_class[CONTROL_MESH_REFINE_LOCALLY_UNKNOWN] = CONTROL;
  data_required[CONTROL_MESH_REFINE_LOCALLY_UNKNOWN] = CONTROL_MESH_REFINE_LOCALLY;

  strcpy(name[CONTROL_MESH_REMESH],"control_mesh_remesh");
  type[CONTROL_MESH_REMESH] = INTEGER;
  data_length[CONTROL_MESH_REMESH] = 1;
  data_class[CONTROL_MESH_REMESH] = CONTROL;

  strcpy(name[CONTROL_MESH_REMESH_FACTOR],"control_mesh_remesh_factor");
  type[CONTROL_MESH_REMESH_FACTOR] = DOUBLE_PRECISION;
  data_length[CONTROL_MESH_REMESH_FACTOR] = 2;
  data_class[CONTROL_MESH_REMESH_FACTOR] = CONTROL;
  data_required[CONTROL_MESH_REMESH_FACTOR] = CONTROL_MESH_REMESH;

  strcpy(name[CONTROL_MESH_RENUMBER],"control_mesh_renumber");
  type[CONTROL_MESH_RENUMBER] = INTEGER;
  data_length[CONTROL_MESH_RENUMBER] = 2;
  data_class[CONTROL_MESH_RENUMBER] = CONTROL;

  strcpy(name[CONTROL_MESH_SPLIT],"control_mesh_split");
  type[CONTROL_MESH_SPLIT] = INTEGER;
  data_length[CONTROL_MESH_SPLIT] = 1;
  data_class[CONTROL_MESH_SPLIT] = CONTROL;

  strcpy(name[CONTROL_MESH_SPLIT_ONLY],"control_mesh_split_only");
  type[CONTROL_MESH_SPLIT_ONLY] = INTEGER;
  data_length[CONTROL_MESH_SPLIT_ONLY] = 2;
  data_class[CONTROL_MESH_SPLIT_ONLY] = CONTROL;

  strcpy(name[CONTROL_OPTIONS_CONVECTION],"control_options_convection");
  type[CONTROL_OPTIONS_CONVECTION] = INTEGER;
  data_length[CONTROL_OPTIONS_CONVECTION] = 1;
  data_class[CONTROL_OPTIONS_CONVECTION] = CONTROL;

  strcpy(name[CONTROL_OPTIONS_INERTIA],"control_options_inertia");
  type[CONTROL_OPTIONS_INERTIA] = INTEGER;
  data_length[CONTROL_OPTIONS_INERTIA] = 1;
  data_class[CONTROL_OPTIONS_INERTIA] = CONTROL;

  strcpy(name[CONTROL_OPTIONS_RELAXATION],"control_options_relaxation");
  type[CONTROL_OPTIONS_RELAXATION] = DOUBLE_PRECISION;
  data_length[CONTROL_OPTIONS_RELAXATION] = 1;
  data_class[CONTROL_OPTIONS_RELAXATION] = CONTROL;

  strcpy(name[CONTROL_OPTIONS_SKIP_GRAVITY],"control_options_skip_gravity");
  type[CONTROL_OPTIONS_SKIP_GRAVITY] = INTEGER;
  data_length[CONTROL_OPTIONS_SKIP_GRAVITY] = 1;
  data_class[CONTROL_OPTIONS_SKIP_GRAVITY] = CONTROL;

  strcpy(name[CONTROL_OPTIONS_SKIP_GROUNDFLOW_NONLINEAR],"control_options_skip_groundflow_nonlinear");
  type[CONTROL_OPTIONS_SKIP_GROUNDFLOW_NONLINEAR] = INTEGER;
  data_length[CONTROL_OPTIONS_SKIP_GROUNDFLOW_NONLINEAR] = 1;
  data_class[CONTROL_OPTIONS_SKIP_GROUNDFLOW_NONLINEAR] = CONTROL;

  strcpy(name[CONTROL_OPTIONS_SKIP_PLASTICITY],"control_options_skip_plasticity");
  type[CONTROL_OPTIONS_SKIP_PLASTICITY] = INTEGER;
  data_length[CONTROL_OPTIONS_SKIP_PLASTICITY] = 1;
  data_class[CONTROL_OPTIONS_SKIP_PLASTICITY] = CONTROL;

  strcpy(name[CONTROL_OPTIONS_SKIP_GROUNDFLOW_NONLINEAR],"control_options_skip_groundflow_nonlinear");
  type[CONTROL_OPTIONS_SKIP_GROUNDFLOW_NONLINEAR] = INTEGER;
  data_length[CONTROL_OPTIONS_SKIP_GROUNDFLOW_NONLINEAR] = 1;
  data_class[CONTROL_OPTIONS_SKIP_GROUNDFLOW_NONLINEAR] = CONTROL;

  strcpy(name[CONTROL_OPTIONS_SOLVER],"control_options_solver");
  type[CONTROL_OPTIONS_SOLVER] = INTEGER;
  data_length[CONTROL_OPTIONS_SOLVER] = 1;
  data_class[CONTROL_OPTIONS_SOLVER] = CONTROL;

  strcpy(name[CONTROL_OPTIONS_SOLVER_BICG_ERROR],"control_options_solver_bicg_error");
  type[CONTROL_OPTIONS_SOLVER_BICG_ERROR] = DOUBLE_PRECISION;
  data_length[CONTROL_OPTIONS_SOLVER_BICG_ERROR] = 1;
  data_class[CONTROL_OPTIONS_SOLVER_BICG_ERROR] = CONTROL;

  strcpy(name[CONTROL_OPTIONS_SOLVER_BICG_ERROR_MINIMUM],"control_options_solver_bicg_error_minimum");
  type[CONTROL_OPTIONS_SOLVER_BICG_ERROR_MINIMUM] = DOUBLE_PRECISION;
  data_length[CONTROL_OPTIONS_SOLVER_BICG_ERROR_MINIMUM] = 1;
  data_class[CONTROL_OPTIONS_SOLVER_BICG_ERROR_MINIMUM] = CONTROL;

  strcpy(name[CONTROL_OPTIONS_SOLVER_PETSC_KSPTYPE],"control_options_solver_petsc_ksptype");
  type[CONTROL_OPTIONS_SOLVER_PETSC_KSPTYPE] = INTEGER;
  data_length[CONTROL_OPTIONS_SOLVER_PETSC_KSPTYPE] = 1;
  data_class[CONTROL_OPTIONS_SOLVER_PETSC_KSPTYPE] = CONTROL;
  data_required[CONTROL_OPTIONS_SOLVER_PETSC_KSPTYPE] = CONTROL_OPTIONS_SOLVER;

  strcpy(name[CONTROL_OPTIONS_SOLVER_PETSC_MG],"control_options_solver_petsc_mg");
  type[CONTROL_OPTIONS_SOLVER_PETSC_MG] = INTEGER;
  data_length[CONTROL_OPTIONS_SOLVER_PETSC_MG] = 1;
  data_class[CONTROL_OPTIONS_SOLVER_PETSC_MG] = CONTROL;
  data_required[CONTROL_OPTIONS_SOLVER_PETSC_MG] = CONTROL_OPTIONS_SOLVER;

  strcpy(name[CONTROL_OPTIONS_SOLVER_PETSC_PCTYPE],"control_options_solver_petsc_pctype");
  type[CONTROL_OPTIONS_SOLVER_PETSC_PCTYPE] = INTEGER;
  data_length[CONTROL_OPTIONS_SOLVER_PETSC_PCTYPE] = 1;
  data_class[CONTROL_OPTIONS_SOLVER_PETSC_PCTYPE] = CONTROL;
  data_required[CONTROL_OPTIONS_SOLVER_PETSC_PCTYPE] = CONTROL_OPTIONS_SOLVER;

  strcpy(name[CONTROL_PRINT],"control_print");
  type[CONTROL_PRINT] = INTEGER;
  data_length[CONTROL_PRINT] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_PRINT] = 0;
  data_class[CONTROL_PRINT] = CONTROL;

  strcpy(name[CONTROL_PRINT_DATABASE],"control_print_database");
  type[CONTROL_PRINT_DATABASE] = INTEGER;
  data_length[CONTROL_PRINT_DATABASE] = 1;
  data_class[CONTROL_PRINT_DATABASE] = CONTROL;

  strcpy(name[CONTROL_PRINT_ELEMENT],"control_print_element");
  type[CONTROL_PRINT_ELEMENT] = INTEGER;
  data_length[CONTROL_PRINT_ELEMENT] = 1;
  data_class[CONTROL_PRINT_ELEMENT] = CONTROL;

  strcpy(name[CONTROL_PRINT_DATA_VERSUS_DATA],"control_print_data_versus_data");
  type[CONTROL_PRINT_DATA_VERSUS_DATA] = INTEGER;
  data_length[CONTROL_PRINT_DATA_VERSUS_DATA] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_PRINT_DATA_VERSUS_DATA] = 0;
  data_class[CONTROL_PRINT_DATA_VERSUS_DATA] = CONTROL;

  strcpy(name[CONTROL_PRINT_DX],"control_print_dx");
  type[CONTROL_PRINT_DX] = INTEGER;
  data_length[CONTROL_PRINT_DX] = 1;
  data_class[CONTROL_PRINT_DX] = CONTROL;

  strcpy(name[CONTROL_PRINT_DX_TIME],"control_print_dx_time");
  type[CONTROL_PRINT_DX_TIME] = DOUBLE_PRECISION;
  data_length[CONTROL_PRINT_DX_TIME] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_PRINT_DX_TIME] = 0;
  data_class[CONTROL_PRINT_DX_TIME] = CONTROL;
  no_index[CONTROL_PRINT_DX_TIME] = 1;

  strcpy(name[CONTROL_PRINT_FILTER],"control_print_filter");
  type[CONTROL_PRINT_FILTER] = INTEGER;
  data_length[CONTROL_PRINT_FILTER] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_PRINT_FILTER] = 0;
  data_class[CONTROL_PRINT_FILTER] = CONTROL;

  strcpy(name[CONTROL_PRINT_GID],"control_print_gid");
  type[CONTROL_PRINT_GID] = INTEGER;
  data_length[CONTROL_PRINT_GID] = 1;
  data_class[CONTROL_PRINT_GID] = CONTROL;

  strcpy(name[CONTROL_PRINT_GID_EMPTY],"control_print_gid_empty");
  type[CONTROL_PRINT_GID_EMPTY] = INTEGER;
  data_length[CONTROL_PRINT_GID_EMPTY] = 1;
  data_class[CONTROL_PRINT_GID_EMPTY] = CONTROL;

  strcpy(name[CONTROL_PRINT_GID_MESH],"control_print_gid_mesh");
  type[CONTROL_PRINT_GID_MESH] = INTEGER;
  data_length[CONTROL_PRINT_GID_MESH] = 1;
  no_index[CONTROL_PRINT_GID_MESH] = 1;
  external[CONTROL_PRINT_GID_MESH] = 0;
  data_class[CONTROL_PRINT_GID_MESH] = CONTROL;

  strcpy(name[CONTROL_PRINT_GID_TIME],"control_print_gid_time");
  type[CONTROL_PRINT_GID_TIME] = DOUBLE_PRECISION;
  data_length[CONTROL_PRINT_GID_TIME] = 1;
  no_index[CONTROL_PRINT_GID_TIME] = 1;
  external[CONTROL_PRINT_GID_TIME] = 0;
  data_class[CONTROL_PRINT_GID_TIME] = CONTROL;

  strcpy(name[CONTROL_PRINT_GMV],"control_print_gmv");
  type[CONTROL_PRINT_GMV] = INTEGER;
  data_length[CONTROL_PRINT_GMV] = 1;
  data_class[CONTROL_PRINT_GMV] = CONTROL;

  strcpy(name[CONTROL_PRINT_GMV_MESH],"control_print_gmv_mesh");
  type[CONTROL_PRINT_GMV_MESH] = INTEGER;
  data_length[CONTROL_PRINT_GMV_MESH] = 1;
  no_index[CONTROL_PRINT_GMV_MESH] = 1;
  external[CONTROL_PRINT_GMV_MESH] = 0;
  data_class[CONTROL_PRINT_GMV_MESH] = CONTROL;

  strcpy(name[CONTROL_PRINT_HISTORY],"control_print_history");
  type[CONTROL_PRINT_HISTORY] = INTEGER;
  data_length[CONTROL_PRINT_HISTORY] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_PRINT_HISTORY] = 0;
  data_class[CONTROL_PRINT_HISTORY] = CONTROL;

  strcpy(name[CONTROL_PRINT_PLOTMTV],"control_print_plotmtv");
  type[CONTROL_PRINT_PLOTMTV] = INTEGER;
  data_length[CONTROL_PRINT_PLOTMTV] = 1;
  data_class[CONTROL_PRINT_PLOTMTV] = CONTROL;

  strcpy(name[CONTROL_PRINT_PLOTMTV_MESH],"control_print_plotmtv_mesh");
  type[CONTROL_PRINT_PLOTMTV_MESH] = INTEGER;
  data_length[CONTROL_PRINT_PLOTMTV_MESH] = 1;
  no_index[CONTROL_PRINT_PLOTMTV_MESH] = 1;
  external[CONTROL_PRINT_PLOTMTV_MESH] = 0;
  data_class[CONTROL_PRINT_PLOTMTV_MESH] = CONTROL;

  strcpy(name[CONTROL_PRINT_MATLAB],"control_print_matlab");
  type[CONTROL_PRINT_MATLAB] = INTEGER;
  data_length[CONTROL_PRINT_MATLAB] = 1;
  data_class[CONTROL_PRINT_MATLAB] = CONTROL;

  strcpy(name[CONTROL_PRINT_TECPLOT],"control_print_tecplot");
  type[CONTROL_PRINT_TECPLOT] = INTEGER;
  data_length[CONTROL_PRINT_TECPLOT] = 1;
  data_class[CONTROL_PRINT_TECPLOT] = CONTROL;

  strcpy(name[CONTROL_PRINT_TECPLOT_MESH],"control_print_tecplot_mesh");
  type[CONTROL_PRINT_TECPLOT_MESH] = INTEGER;
  data_length[CONTROL_PRINT_TECPLOT_MESH] = 1;
  no_index[CONTROL_PRINT_TECPLOT_MESH] = 1;
  external[CONTROL_PRINT_TECPLOT_MESH] = 0;
  data_class[CONTROL_PRINT_TECPLOT_MESH] = CONTROL;

  strcpy(name[CONTROL_PRINT_UNKNOWNS],"control_print_unknowns");
  type[CONTROL_PRINT_UNKNOWNS] = INTEGER;
  data_length[CONTROL_PRINT_UNKNOWNS] = 1;
  data_class[CONTROL_PRINT_UNKNOWNS] = CONTROL;

  strcpy(name[CONTROL_PRINT_UNKNOWNSRHSIDE],"control_print_unknownsrhside");
  type[CONTROL_PRINT_UNKNOWNSRHSIDE] = INTEGER;
  data_length[CONTROL_PRINT_UNKNOWNSRHSIDE] = 1;
  data_class[CONTROL_PRINT_UNKNOWNSRHSIDE] = CONTROL;

  strcpy(name[CONTROL_PRINT_VTK],"control_print_vtk");
  type[CONTROL_PRINT_VTK] = INTEGER;
  data_length[CONTROL_PRINT_VTK] = 1;
  data_class[CONTROL_PRINT_VTK] = CONTROL;

  strcpy(name[CONTROL_RELAXATION_CONDIF_TEMPERATURE],"control_relaxation_condif_temperature");
  type[CONTROL_RELAXATION_CONDIF_TEMPERATURE] = DOUBLE_PRECISION;
  data_length[CONTROL_RELAXATION_CONDIF_TEMPERATURE] = 1;
  data_class[CONTROL_RELAXATION_CONDIF_TEMPERATURE] = CONTROL;

  strcpy(name[CONTROL_RELAXATION_GROUNDFLOW_PRESSURE],"control_relaxation_groundflow_pressure");
  type[CONTROL_RELAXATION_GROUNDFLOW_PRESSURE] = DOUBLE_PRECISION;
  data_length[CONTROL_RELAXATION_GROUNDFLOW_PRESSURE] = 1;
  data_class[CONTROL_RELAXATION_GROUNDFLOW_PRESSURE] = CONTROL;

  strcpy(name[CONTROL_RELAXATION_MAXWELL_E],"control_relaxation_maxwell_e");
  type[CONTROL_RELAXATION_MAXWELL_E] = DOUBLE_PRECISION;
  data_length[CONTROL_RELAXATION_MAXWELL_E] = 1;
  data_class[CONTROL_RELAXATION_MAXWELL_E] = CONTROL;      

  strcpy(name[CONTROL_RELAXATION_MATERI_VELOCITY],"control_relaxation_materi_velocity");
  type[CONTROL_RELAXATION_MATERI_VELOCITY] = DOUBLE_PRECISION;
  data_length[CONTROL_RELAXATION_MATERI_VELOCITY] = 1;
  data_class[CONTROL_RELAXATION_MATERI_VELOCITY] = CONTROL;

  strcpy(name[CONTROL_RELAXATION_WAVE_FSCALAR],"control_relaxation_wave_fscalar");
  type[CONTROL_RELAXATION_WAVE_FSCALAR] = DOUBLE_PRECISION;
  data_length[CONTROL_RELAXATION_WAVE_FSCALAR] = 1;
  data_class[CONTROL_RELAXATION_WAVE_FSCALAR] = CONTROL;

  strcpy(name[CONTROL_REPEAT],"control_repeat");
  type[CONTROL_REPEAT] = INTEGER;
  data_length[CONTROL_REPEAT] = 2;
  data_class[CONTROL_REPEAT] = CONTROL;

  strcpy(name[CONTROL_REPEAT_UNTIL_ITEM],"control_repeat_until_item");
  type[CONTROL_REPEAT_UNTIL_ITEM] = INTEGER;
  data_length[CONTROL_REPEAT_UNTIL_ITEM] = 5;
  data_class[CONTROL_REPEAT_UNTIL_ITEM] = CONTROL;

  strcpy(name[CONTROL_REPEAT_UNTIL_TOLERANCE],"control_repeat_until_tolerance");
  type[CONTROL_REPEAT_UNTIL_TOLERANCE] = DOUBLE_PRECISION;
  data_length[CONTROL_REPEAT_UNTIL_TOLERANCE] = 1;
  data_class[CONTROL_REPEAT_UNTIL_TOLERANCE] = CONTROL;

  strcpy(name[CONTROL_REPEAT_UNTIL_VALUE],"control_repeat_until_value");
  type[CONTROL_REPEAT_UNTIL_VALUE] = DOUBLE_PRECISION;
  data_length[CONTROL_REPEAT_UNTIL_VALUE] = 1;
  external[CONTROL_REPEAT_UNTIL_VALUE] = 0;
  data_class[CONTROL_REPEAT_UNTIL_VALUE] = CONTROL;

  strcpy(name[CONTROL_RESTART],"control_restart");
  type[CONTROL_RESTART] = INTEGER;
  data_length[CONTROL_RESTART] = 1;
  data_class[CONTROL_RESTART] = CONTROL;

  strcpy(name[CONTROL_UNKNOWN_FREEZE],"control_unknown_freeze");
  type[CONTROL_UNKNOWN_FREEZE] = INTEGER;
  data_length[CONTROL_UNKNOWN_FREEZE] = nuknwn;
  data_class[CONTROL_UNKNOWN_FREEZE] = CONTROL;
  fixed_length[CONTROL_UNKNOWN_FREEZE] = 0;

  strcpy(name[CONTROL_UNKNOWN_RESET_GEOMETRY],"control_unknown_reset_geometry");
  type[CONTROL_UNKNOWN_RESET_GEOMETRY] = INTEGER;
  data_length[CONTROL_UNKNOWN_RESET_GEOMETRY] =2;
  data_class[CONTROL_UNKNOWN_RESET_GEOMETRY] = CONTROL;

  strcpy(name[CONTROL_UNKNOWN_RESET_UNKNOWN],"control_unknown_reset_unknown");
  type[CONTROL_UNKNOWN_RESET_UNKNOWN] = INTEGER;
  data_length[CONTROL_UNKNOWN_RESET_UNKNOWN] = nuknwn;
  data_class[CONTROL_UNKNOWN_RESET_UNKNOWN] = CONTROL;
  fixed_length[CONTROL_UNKNOWN_RESET_UNKNOWN] = 0;

  strcpy(name[CONTROL_UNKNOWN_RESET_VALUE],"control_unknown_reset_value");
  type[CONTROL_UNKNOWN_RESET_VALUE] = DOUBLE_PRECISION;
  data_length[CONTROL_UNKNOWN_RESET_VALUE] = 1;
  data_class[CONTROL_UNKNOWN_RESET_VALUE] = CONTROL;

  strcpy(name[CONTROL_TIMESTEP],"control_timestep");
  type[CONTROL_TIMESTEP] = DOUBLE_PRECISION;
  data_length[CONTROL_TIMESTEP] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_TIMESTEP] = 0;
  data_class[CONTROL_TIMESTEP] = CONTROL;

  strcpy(name[CONTROL_TIMESTEP_ITERATIONS],"control_timestep_iterations");
  type[CONTROL_TIMESTEP_ITERATIONS] = INTEGER;
  data_length[CONTROL_TIMESTEP_ITERATIONS] = 1;
  data_class[CONTROL_TIMESTEP_ITERATIONS] = CONTROL;
  data_required[CONTROL_TIMESTEP_ITERATIONS] = CONTROL_TIMESTEP;

  strcpy(name[CONTROL_TIMESTEP_ITERATIONS_AUTOMATIC],"control_timestep_iterations_automatic");
  type[CONTROL_TIMESTEP_ITERATIONS_AUTOMATIC] = DOUBLE_PRECISION;
  data_length[CONTROL_TIMESTEP_ITERATIONS_AUTOMATIC] = 2;
  data_class[CONTROL_TIMESTEP_ITERATIONS_AUTOMATIC] = CONTROL;
  data_required[CONTROL_TIMESTEP_ITERATIONS_AUTOMATIC] = CONTROL_TIMESTEP;

  strcpy(name[CONTROL_TIMESTEP_ITERATIONS_AUTOMATIC_STOP],"control_timestep_iterations_automatic_stop");
  type[CONTROL_TIMESTEP_ITERATIONS_AUTOMATIC_STOP] = INTEGER;
  data_length[CONTROL_TIMESTEP_ITERATIONS_AUTOMATIC_STOP] = 1;
  data_class[CONTROL_TIMESTEP_ITERATIONS_AUTOMATIC_STOP] = CONTROL;
  data_required[CONTROL_TIMESTEP_ITERATIONS_AUTOMATIC_STOP] = CONTROL_TIMESTEP;

  strcpy(name[CONTROL_TIMESTEP_SIZE_AUTOMATIC_DECREASE],"control_timestep_size_automatic_decrease");
  type[CONTROL_TIMESTEP_SIZE_AUTOMATIC_DECREASE] = DOUBLE_PRECISION;
  data_length[CONTROL_TIMESTEP_SIZE_AUTOMATIC_DECREASE] = 3;
  data_class[CONTROL_TIMESTEP_SIZE_AUTOMATIC_DECREASE] = CONTROL;
  data_required[CONTROL_TIMESTEP_SIZE_AUTOMATIC_DECREASE] = CONTROL_TIMESTEP;

  strcpy(name[CONTROL_TIMESTEP_MULTIPLIER],"control_timestep_multiplier");
  type[CONTROL_TIMESTEP_MULTIPLIER] = DOUBLE_PRECISION;
  data_length[CONTROL_TIMESTEP_MULTIPLIER] = 1;
  fixed_length[CONTROL_TIMESTEP_MULTIPLIER] = 0;
  data_class[CONTROL_TIMESTEP_MULTIPLIER] = CONTROL;

  strcpy(name[CR],"cr");

  strcpy(name[CRACK_DIRECTION],"crack_direction");
  type[CRACK_DIRECTION] = DOUBLE_PRECISION;
  data_length[CRACK_DIRECTION] = ndim;
  data_class[CRACK_DIRECTION] = CRACK;
  no_index[CRACK_DIRECTION] = 1;

  strcpy(name[CRACK_NODES],"crack_nodes");
  type[CRACK_NODES] = INTEGER;
  data_length[CRACK_NODES] = 5;
  data_class[CRACK_NODES] = CRACK;
  no_index[CRACK_NODES] = 1;

  strcpy(name[CRACK_ELEMENTGROUP],"crack_elementgroup");
  type[CRACK_ELEMENTGROUP] = INTEGER;
  data_length[CRACK_ELEMENTGROUP] = 1;
  data_class[CRACK_ELEMENTGROUP] = CRACK;
  no_index[CRACK_ELEMENTGROUP] = 1;

  strcpy(name[CRACK_LENGTH],"crack_length");
  type[CRACK_LENGTH] = DOUBLE_PRECISION;
  data_length[CRACK_LENGTH] = ndim;
  data_class[CRACK_LENGTH] = CRACK;
  no_index[CRACK_LENGTH] = 1;

  strcpy(name[CRACK_STRESSINTENSITYFACTOR],"crack_stressintensityfactor");
  type[CRACK_STRESSINTENSITYFACTOR] = DOUBLE_PRECISION;
  data_length[CRACK_STRESSINTENSITYFACTOR] = 2;
  data_class[CRACK_STRESSINTENSITYFACTOR] = CRACK;
  no_index[CRACK_STRESSINTENSITYFACTOR] = 1;

  strcpy(name[CRACK_TIP],"crack_tip");
  type[CRACK_TIP] = DOUBLE_PRECISION;
  data_length[CRACK_TIP] = ndim;
  data_class[CRACK_TIP] = CRACK;
  no_index[CRACK_TIP] = 1;

  strcpy(name[CYLINDER_HOLLOW],"cylinder_hollow" );

  strcpy(name[DATABASE],"database");

  strcpy(name[DEPENDENCY],"dependency");

  strcpy(name[DEPENDENCY_DIAGRAM],"dependency_diagram");
  type[DEPENDENCY_DIAGRAM] = DOUBLE_PRECISION;
  data_length[DEPENDENCY_DIAGRAM] = DATA_ITEM_SIZE;
  fixed_length[DEPENDENCY_DIAGRAM] = 0;
  data_class[DEPENDENCY_DIAGRAM] = DEPENDENCY;
  data_required[DEPENDENCY_DIAGRAM] = DEPENDENCY_ITEM;

  strcpy(name[DEPENDENCY_ITEM],"dependency_item");
  type[DEPENDENCY_ITEM] = INTEGER;
  data_length[DEPENDENCY_ITEM] = 4;
  data_class[DEPENDENCY_ITEM] = DEPENDENCY;
  data_required[DEPENDENCY_ITEM] = DEPENDENCY_DIAGRAM;

  strcpy(name[DIAGONAL],"diagonal");

  strcpy(name[EMPTY],"empty" );

  strcpy(name[MATERI_DISPLACEMENT],"materi_displacement");

  strcpy(name[DOF],"dof");

  strcpy(name[DOF_AMOUNT],"dof_amount");
  type[DOF_AMOUNT] = INTEGER;
  data_length[DOF_AMOUNT] = nuknwn;
  no_index[DOF_AMOUNT] = 1;
  external[DOF_AMOUNT] = 0;
  data_class[DOF_AMOUNT] = DOF;

  strcpy(name[DOF_LABEL],"dof_label");
  type[DOF_LABEL] = INTEGER;
  data_length[DOF_LABEL] = nuknwn;
  print_only[DOF_LABEL] = 1;
  no_index[DOF_LABEL] = 1;
  data_class[DOF_LABEL] = DOF;

  strcpy(name[DOF_PRINCIPAL],"dof_principal");
  type[DOF_PRINCIPAL] = INTEGER;
  data_length[DOF_PRINCIPAL] = nuknwn;
  no_index[DOF_PRINCIPAL] = 1;
  external[DOF_PRINCIPAL] = 0;
  data_class[DOF_PRINCIPAL] = DOF;

  strcpy(name[DOF_SCAL_VEC_MAT],"dof_scal_vec_mat");
  type[DOF_SCAL_VEC_MAT] = INTEGER;
  data_length[DOF_SCAL_VEC_MAT] = nuknwn;
  no_index[DOF_SCAL_VEC_MAT] = 1;
  external[DOF_SCAL_VEC_MAT] = 0;
  data_class[DOF_SCAL_VEC_MAT] = DOF;

  strcpy(name[DOF_TYPE],"dof_type");
  type[DOF_TYPE] = INTEGER;
  data_length[DOF_TYPE] = nuknwn;
  no_index[DOF_TYPE] = 1;
  external[DOF_TYPE] = 0;
  data_class[DOF_TYPE] = DOF;

  strcpy(name[DTIME],"dtime");
  type[DTIME] = DOUBLE_PRECISION;
  data_length[DTIME] = 1;
  no_index[DTIME] = 1;
  version_all[DTIME] = 1;
  external[DTIME] = 0;

  strcpy(name[DYNAMIC],"dynamic" );

  strcpy(name[ELEMENT],"element");
  type[ELEMENT] = INTEGER;
  if      ( ndim==1 ) data_length[ELEMENT] = 1+4;
  else if ( ndim==2 ) data_length[ELEMENT] = 1+16;
  else if ( ndim==3 ) data_length[ELEMENT] = 1+MNOL;
  fixed_length[ELEMENT] = 0;
  version_all[ELEMENT] = 1;
  data_class[ELEMENT] = ELEMENT;

  strcpy(name[ELEMENT_ADJUST],"element_adjust");
  type[ELEMENT_ADJUST] = INTEGER;
  data_length[ELEMENT_ADJUST] = 1;
  version_all[ELEMENT_ADJUST] = 1;
  data_class[ELEMENT_ADJUST] = ELEMENT;
  external[ELEMENT_ADJUST] = 0;
  data_required[ELEMENT_ADJUST] = ELEMENT;

  strcpy(name[ELEMENT_BEAM_DIRECTION],"element_beam_direction");
  type[ELEMENT_BEAM_DIRECTION] = DOUBLE_PRECISION;
  data_length[ELEMENT_BEAM_DIRECTION] = 2;
  version_all[ELEMENT_BEAM_DIRECTION] = 1;
  print_only[ELEMENT_BEAM_DIRECTION] = 1;
  data_class[ELEMENT_BEAM_DIRECTION] = ELEMENT;
  data_required[ELEMENT_BEAM_DIRECTION] = ELEMENT;

  strcpy(name[ELEMENT_BEAM_MOMENT],"element_beam_moment");
  type[ELEMENT_BEAM_MOMENT] = DOUBLE_PRECISION;
  data_length[ELEMENT_BEAM_MOMENT] = 6;
  version_all[ELEMENT_BEAM_MOMENT] = 1;
  data_class[ELEMENT_BEAM_MOMENT] = ELEMENT;
  data_required[ELEMENT_BEAM_MOMENT] = ELEMENT;

  strcpy(name[ELEMENT_CONTACTSPRING_DIRECTION],"element_contactspring_direction");
  type[ELEMENT_CONTACTSPRING_DIRECTION] = DOUBLE_PRECISION;
  data_length[ELEMENT_CONTACTSPRING_DIRECTION] = MDIM*MDIM;
  version_all[ELEMENT_CONTACTSPRING_DIRECTION] = 1;
  print_only[ELEMENT_CONTACTSPRING_DIRECTION] = 1;
  data_class[ELEMENT_CONTACTSPRING_DIRECTION] = ELEMENT;
  data_required[ELEMENT_CONTACTSPRING_DIRECTION] = ELEMENT;

  strcpy(name[ELEMENT_CONTACTSPRING_FORCE],"element_contactspring_force");
  type[ELEMENT_CONTACTSPRING_FORCE] = DOUBLE_PRECISION;
  data_length[ELEMENT_CONTACTSPRING_FORCE] = MDIM;
  version_all[ELEMENT_CONTACTSPRING_FORCE] = 1;
  print_only[ELEMENT_CONTACTSPRING_FORCE] = 1;
  data_class[ELEMENT_CONTACTSPRING_FORCE] = ELEMENT;
  data_required[ELEMENT_CONTACTSPRING_FORCE] = ELEMENT;

  strcpy(name[ELEMENT_DELETE_FACTOR],"element_delete_factor");
  type[ELEMENT_DELETE_FACTOR] = DOUBLE_PRECISION;
  data_length[ELEMENT_DELETE_FACTOR] = 1;
  external[ELEMENT_DELETE_FACTOR] = 0;
  version_all[ELEMENT_DELETE_FACTOR] = 1;
  data_class[ELEMENT_DELETE_FACTOR] = ELEMENT;
  data_required[ELEMENT_DELETE_FACTOR] = ELEMENT;

  strcpy(name[ELEMENT_DELETE_TIMES],"element_delete_times");
  type[ELEMENT_DELETE_TIMES] = DOUBLE_PRECISION;
  data_length[ELEMENT_DELETE_TIMES] = 2;
  external[ELEMENT_DELETE_TIMES] = 0;
  version_all[ELEMENT_DELETE_TIMES] = 1;
  data_class[ELEMENT_DELETE_TIMES] = ELEMENT;
  data_required[ELEMENT_DELETE_TIMES] = ELEMENT;

  strcpy(name[ELEMENT_DISTRIBUTE],"element_distribute");
  type[ELEMENT_DISTRIBUTE] = INTEGER;
  data_length[ELEMENT_DISTRIBUTE] = 2;
  external[ELEMENT_DISTRIBUTE] = 0;
  version_all[ELEMENT_DISTRIBUTE] = 1;
  data_class[ELEMENT_DISTRIBUTE] = ELEMENT;
  data_required[ELEMENT_DISTRIBUTE] = ELEMENT;

  strcpy(name[ELEMENT_DISTRIBUTE_VALUES],"element_distribute_values");
  type[ELEMENT_DISTRIBUTE_VALUES] = DOUBLE_PRECISION;
  data_length[ELEMENT_DISTRIBUTE_VALUES] = 1;
  external[ELEMENT_DISTRIBUTE_VALUES] = 0;
  version_all[ELEMENT_DISTRIBUTE_VALUES] = 1;
  data_class[ELEMENT_DISTRIBUTE_VALUES] = ELEMENT;
  data_required[ELEMENT_DISTRIBUTE_VALUES] = ELEMENT;

  strcpy(name[ELEMENT_DOF],"element_dof");
  type[ELEMENT_DOF] = DOUBLE_PRECISION;
  data_length[ELEMENT_DOF] = npointmax*nuknwn;
  version_all[ELEMENT_DOF] = 1;
  data_class[ELEMENT_DOF] = ELEMENT;
  data_required[ELEMENT_DOF] = ELEMENT;

  strcpy(name[ELEMENT_DOF_INITIALISED],"element_dof_initialised");
  type[ELEMENT_DOF_INITIALISED] = INTEGER;
  data_length[ELEMENT_DOF_INITIALISED] = 1;
  version_all[ELEMENT_DOF_INITIALISED] = 1;
  data_class[ELEMENT_DOF_INITIALISED] = ELEMENT;
  data_required[ELEMENT_DOF_INITIALISED] = ELEMENT;

  strcpy(name[ELEMENT_EMPTY],"element_empty");
  type[ELEMENT_EMPTY] = INTEGER;
  data_length[ELEMENT_EMPTY] = 1;
  version_all[ELEMENT_EMPTY] = 1;
  external[ELEMENT_EMPTY] = 0;   
  data_class[ELEMENT_EMPTY] = ELEMENT;   
  data_required[ELEMENT_EMPTY] = ELEMENT;   

  strcpy(name[ELEMENT_GROUP],"element_group");
  type[ELEMENT_GROUP] = INTEGER;
  data_length[ELEMENT_GROUP] = 1;
  version_all[ELEMENT_GROUP] = 1;
  data_class[ELEMENT_GROUP] = ELEMENT;
  data_required[ELEMENT_GROUP] = ELEMENT;

  strcpy(name[ELEMENT_GROUP_AREA_ELEMENT_GROUP],
    "element_group_area_element_group");
  type[ELEMENT_GROUP_AREA_ELEMENT_GROUP] = INTEGER;
  data_length[ELEMENT_GROUP_AREA_ELEMENT_GROUP] = 1;
  version_all[ELEMENT_GROUP_AREA_ELEMENT_GROUP] = 1;
  data_class[ELEMENT_GROUP_AREA_ELEMENT_GROUP] = ELEMENT;
  data_required[ELEMENT_GROUP_AREA_ELEMENT_GROUP] = ELEMENT;

  strcpy(name[ELEMENT_GROUP_AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP],
    "element_group_area_element_group_sequence_elementgroup");
  type[ELEMENT_GROUP_AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP] = INTEGER;
  data_length[ELEMENT_GROUP_AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP] = 1;
  version_all[ELEMENT_GROUP_AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP] = 1;
  data_class[ELEMENT_GROUP_AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP] = ELEMENT;
  data_required[ELEMENT_GROUP_AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP] = ELEMENT;

  strcpy(name[ELEMENT_MACRO_GENERATE],"element_macro_generate");
  type[ELEMENT_MACRO_GENERATE] = INTEGER;
  data_length[ELEMENT_MACRO_GENERATE] = 1;
  version_all[ELEMENT_MACRO_GENERATE] = 1;
  data_class[ELEMENT_MACRO_GENERATE] = ELEMENT;
  data_required[ELEMENT_MACRO_GENERATE] = ELEMENT;
  external[ELEMENT_MACRO_GENERATE] = 0;

  strcpy(name[ELEMENT_MASS],"element_mass");
  type[ELEMENT_MASS] = DOUBLE_PRECISION;
  data_length[ELEMENT_MASS] = 1;
  external[ELEMENT_MASS] = 0;
  data_class[ELEMENT_MASS] = ELEMENT;
  data_required[ELEMENT_MASS] = ELEMENT;

  strcpy(name[ELEMENT_MATRIX_DELETE],"element_matrix_delete");
  type[ELEMENT_MATRIX_DELETE] = DOUBLE_PRECISION;
  data_length[ELEMENT_MATRIX_DELETE] = MNOL*npuknwn;
  fixed_length[ELEMENT_MATRIX_DELETE] = 0;
  external[ELEMENT_MATRIX_DELETE] = 0;
  data_class[ELEMENT_MATRIX_DELETE] = ELEMENT;
  data_required[ELEMENT_MATRIX_DELETE] = ELEMENT;

  strcpy(name[ELEMENT_MATRIX_UNKNOWNS],"element_matrix_unknowns");
  type[ELEMENT_MATRIX_UNKNOWNS] = INTEGER;
  data_length[ELEMENT_MATRIX_UNKNOWNS] = 0; // set run-time
  fixed_length[ELEMENT_MATRIX_UNKNOWNS] = 0;
  external[ELEMENT_MATRIX_UNKNOWNS] = 0;
  data_class[ELEMENT_MATRIX_UNKNOWNS] = ELEMENT;
  data_required[ELEMENT_MATRIX_UNKNOWNS] = ELEMENT;

  strcpy(name[ELEMENT_MATRIX_VALUES],"element_matrix_values");
  type[ELEMENT_MATRIX_VALUES] = DOUBLE_PRECISION;
  data_length[ELEMENT_MATRIX_VALUES] = 0; // set run-time
  fixed_length[ELEMENT_MATRIX_VALUES] = 0;
  external[ELEMENT_MATRIX_VALUES] = 0;
  data_class[ELEMENT_MATRIX_VALUES] = ELEMENT;
  data_required[ELEMENT_MATRIX_VALUES] = ELEMENT;

  strcpy(name[ELEMENT_MATRIX_SECOND_VALUES],"element_matrix_second_values");
  type[ELEMENT_MATRIX_SECOND_VALUES] = DOUBLE_PRECISION;
  data_length[ELEMENT_MATRIX_SECOND_VALUES] = 0; // set run-time
  fixed_length[ELEMENT_MATRIX_SECOND_VALUES] = 0;
  external[ELEMENT_MATRIX_SECOND_VALUES] = 0;
  data_class[ELEMENT_MATRIX_SECOND_VALUES] = ELEMENT;
  data_required[ELEMENT_MATRIX_SECOND_VALUES] = ELEMENT;

  strcpy(name[ELEMENT_MIDDLE],"element_middle");
  type[ELEMENT_MIDDLE] = DOUBLE_PRECISION;
  data_length[ELEMENT_MIDDLE] = ndim;
  external[ELEMENT_MIDDLE] = 0;
  version_all[ELEMENT_MIDDLE] = 1;
  data_class[ELEMENT_MIDDLE] = ELEMENT;
  data_required[ELEMENT_MIDDLE] = ELEMENT;

  strcpy(name[ELEMENT_NONLOCAL],"element_nonlocal");
  type[ELEMENT_NONLOCAL] = INTEGER;
  data_length[ELEMENT_NONLOCAL] = NONLOCAL_ITEM_SIZE*npointmax;
  external[ELEMENT_NONLOCAL] = 0;
  version_all[ELEMENT_NONLOCAL] = 0;
  fixed_length[ELEMENT_NONLOCAL] = 0;
  data_class[ELEMENT_NONLOCAL] = ELEMENT;
  data_required[ELEMENT_NONLOCAL] = ELEMENT;

  strcpy(name[ELEMENT_NONLOCAL_IPOINT],"element_nonlocal_ipoint");
  type[ELEMENT_NONLOCAL_IPOINT] = INTEGER;
  data_length[ELEMENT_NONLOCAL_IPOINT] = NONLOCAL_ITEM_SIZE*npointmax;
  external[ELEMENT_NONLOCAL_IPOINT] = 0;
  version_all[ELEMENT_NONLOCAL_IPOINT] = 0;
  fixed_length[ELEMENT_NONLOCAL_IPOINT] = 0;
  data_class[ELEMENT_NONLOCAL_IPOINT] = ELEMENT;
  data_required[ELEMENT_NONLOCAL_IPOINT] = ELEMENT;

  strcpy(name[ELEMENT_NONLOCAL_WEIGHT],"element_nonlocal_weight");
  type[ELEMENT_NONLOCAL_WEIGHT] = DOUBLE_PRECISION;
  data_length[ELEMENT_NONLOCAL_WEIGHT] = NONLOCAL_ITEM_SIZE*npointmax;
  external[ELEMENT_NONLOCAL_WEIGHT] = 0;
  version_all[ELEMENT_NONLOCAL_WEIGHT] = 0;
  fixed_length[ELEMENT_NONLOCAL_WEIGHT] = 0;
  data_class[ELEMENT_NONLOCAL_WEIGHT] = ELEMENT;
  data_required[ELEMENT_NONLOCAL_WEIGHT] = ELEMENT;

  strcpy(name[ELEMENT_RADIUS],"element_radius");
  type[ELEMENT_RADIUS] = DOUBLE_PRECISION;
  data_length[ELEMENT_RADIUS] = 1;
  external[ELEMENT_RADIUS] = 0;
  version_all[ELEMENT_RADIUS] = 1;
  data_class[ELEMENT_RADIUS] = ELEMENT;
  data_required[ELEMENT_RADIUS] = ELEMENT;

  strcpy(name[ELEMENT_RHSIDE_DELETE],"element_rhside_delete");
  type[ELEMENT_RHSIDE_DELETE] = DOUBLE_PRECISION;
  data_length[ELEMENT_RHSIDE_DELETE] = MNOL*npuknwn;
  fixed_length[ELEMENT_RHSIDE_DELETE] = 0;
  external[ELEMENT_RHSIDE_DELETE] = 0;
  data_class[ELEMENT_RHSIDE_DELETE] = ELEMENT;
  data_required[ELEMENT_RHSIDE_DELETE] = ELEMENT;

  strcpy(name[ELEMENT_SPRING_DIRECTION],"element_spring_direction");
  type[ELEMENT_SPRING_DIRECTION] = DOUBLE_PRECISION;
  data_length[ELEMENT_SPRING_DIRECTION] = ndim;
  version_all[ELEMENT_SPRING_DIRECTION] = 1;
  print_only[ELEMENT_SPRING_DIRECTION] = 1;
  data_class[ELEMENT_SPRING_DIRECTION] = ELEMENT;
  data_required[ELEMENT_SPRING_DIRECTION] = ELEMENT;

  strcpy(name[ELEMENT_SPRING_FORCE],"element_spring_force");
  type[ELEMENT_SPRING_FORCE] = DOUBLE_PRECISION;
  data_length[ELEMENT_SPRING_FORCE] = 1;
  version_all[ELEMENT_SPRING_FORCE] = 1;
  print_only[ELEMENT_SPRING_FORCE] = 1;
  data_class[ELEMENT_SPRING_FORCE] = ELEMENT;
  data_required[ELEMENT_SPRING_FORCE] = ELEMENT;

  strcpy(name[ELEMENT_STRAINENERGY],"element_strainenergy");
  type[ELEMENT_STRAINENERGY] = DOUBLE_PRECISION;
  data_length[ELEMENT_STRAINENERGY] = 1;
  external[ELEMENT_STRAINENERGY] = 0;
  data_class[ELEMENT_STRAINENERGY] = ELEMENT;
  data_required[ELEMENT_STRAINENERGY] = ELEMENT;

  strcpy(name[ELEMENT_TENDON_DIRECTION],"element_tendon_direction");
  type[ELEMENT_TENDON_DIRECTION] = DOUBLE_PRECISION;
  data_length[ELEMENT_TENDON_DIRECTION] = MTENDON * MDIM;
  fixed_length[ELEMENT_TENDON_DIRECTION] = 0;
  version_all[ELEMENT_TENDON_DIRECTION] = 1;
  print_only[ELEMENT_TENDON_DIRECTION] = 1;
  data_class[ELEMENT_TENDON_DIRECTION] = ELEMENT;
  data_required[ELEMENT_TENDON_DIRECTION] = ELEMENT;

  strcpy(name[ELEMENT_TENDON_INTERSECTIONS],"element_tendon_intersections");
  type[ELEMENT_TENDON_INTERSECTIONS] = DOUBLE_PRECISION;
  data_length[ELEMENT_TENDON_INTERSECTIONS] = MTENDON * 2 * ndim;
  fixed_length[ELEMENT_TENDON_INTERSECTIONS] = 0;
  version_all[ELEMENT_TENDON_INTERSECTIONS] = 1;
  external[ELEMENT_TENDON_INTERSECTIONS] = 0;
  data_class[ELEMENT_TENDON_INTERSECTIONS] = ELEMENT;
  data_required[ELEMENT_TENDON_INTERSECTIONS] = ELEMENT;

  strcpy(name[ELEMENT_TENDON_NUMBER],"element_tendon_number");
  type[ELEMENT_TENDON_NUMBER] = INTEGER;
  data_length[ELEMENT_TENDON_NUMBER] = MTENDON;
  fixed_length[ELEMENT_TENDON_NUMBER] = 0;
  version_all[ELEMENT_TENDON_NUMBER] = 1;
  print_only[ELEMENT_TENDON_NUMBER] = 1;
  data_class[ELEMENT_TENDON_NUMBER] = ELEMENT;
  data_required[ELEMENT_TENDON_NUMBER] = ELEMENT;

  strcpy(name[ELEMENT_TENDON_STRAIN],"element_tendon_strain");
  type[ELEMENT_TENDON_STRAIN] = DOUBLE_PRECISION;
  data_length[ELEMENT_TENDON_STRAIN] = MTENDON;
  fixed_length[ELEMENT_TENDON_STRAIN] = 0;
  version_all[ELEMENT_TENDON_STRAIN] = 1;
  print_only[ELEMENT_TENDON_STRAIN] = 1;
  data_class[ELEMENT_TENDON_STRAIN] = ELEMENT;
  data_required[ELEMENT_TENDON_STRAIN] = ELEMENT;

  strcpy(name[ELEMENT_TENDON_STRESS],"element_tendon_stress");
  type[ELEMENT_TENDON_STRESS] = DOUBLE_PRECISION;
  data_length[ELEMENT_TENDON_STRESS] = MTENDON;
  fixed_length[ELEMENT_TENDON_STRESS] = 0;
  version_all[ELEMENT_TENDON_STRESS] = 1;
  print_only[ELEMENT_TENDON_STRESS] = 1;
  data_class[ELEMENT_TENDON_STRESS] = ELEMENT;
  data_required[ELEMENT_TENDON_STRESS] = ELEMENT;

  strcpy(name[ELEMENT_TENDON_VOLUME],"element_tendon_volume");
  type[ELEMENT_TENDON_VOLUME] = DOUBLE_PRECISION;
  data_length[ELEMENT_TENDON_VOLUME] = MTENDON;
  fixed_length[ELEMENT_TENDON_VOLUME] = 0;
  version_all[ELEMENT_TENDON_VOLUME] = 1;
  print_only[ELEMENT_TENDON_VOLUME] = 1;
  data_class[ELEMENT_TENDON_VOLUME] = ELEMENT;
  data_required[ELEMENT_TENDON_VOLUME] = ELEMENT;

  strcpy(name[ELEMENT_TRUSS_DIRECTION],"element_truss_direction");
  type[ELEMENT_TRUSS_DIRECTION] = DOUBLE_PRECISION;
  data_length[ELEMENT_TRUSS_DIRECTION] = ndim;
  version_all[ELEMENT_TRUSS_DIRECTION] = 1;
  print_only[ELEMENT_TRUSS_DIRECTION] = 1;
  data_class[ELEMENT_TRUSS_DIRECTION] = ELEMENT;
  data_required[ELEMENT_TRUSS_DIRECTION] = ELEMENT;

  strcpy(name[ELEMENT_TRUSS_FORCE],"element_truss_force");
  type[ELEMENT_TRUSS_FORCE] = DOUBLE_PRECISION;
  data_length[ELEMENT_TRUSS_FORCE] = 1;
  version_all[ELEMENT_TRUSS_FORCE] = 1;
  data_class[ELEMENT_TRUSS_FORCE] = ELEMENT;
  data_required[ELEMENT_TRUSS_FORCE] = ELEMENT;

  strcpy(name[ELEMENT_VOLUME],"element_volume");
  type[ELEMENT_VOLUME] = DOUBLE_PRECISION;
  data_length[ELEMENT_VOLUME] = 1;
  external[ELEMENT_VOLUME] = 0;
  data_class[ELEMENT_VOLUME] = ELEMENT;
  data_required[ELEMENT_VOLUME] = ELEMENT;

  strcpy(name[EISENSTAT],"eisenstat");

  strcpy(name[EVERYTHING],"everything" );

  strcpy(name[EXIT_TOCHNOG],"exit_tochnog");
  type[EXIT_TOCHNOG] = INTEGER;
  data_length[EXIT_TOCHNOG] = 1;
  data_class[EXIT_TOCHNOG] = EXIT_TOCHNOG;
  no_index[EXIT_TOCHNOG] = 1;

  strcpy(name[FIXED_IN_SPACE],"fixed_in_space");

  strcpy(name[FOLLOW_MATERIAL],"follow_material");

  strcpy(name[FORCE],"force");

  strcpy(name[FORCE_ELEMENT_EDGE],"force_element_edge");
  type[FORCE_ELEMENT_EDGE] = DOUBLE_PRECISION;
  data_length[FORCE_ELEMENT_EDGE] = nprinc;
  data_class[FORCE_ELEMENT_EDGE] = FORCE;

  strcpy(name[FORCE_ELEMENT_EDGE_FACTOR],"force_element_edge_factor");
  type[FORCE_ELEMENT_EDGE_FACTOR] = DOUBLE_PRECISION;
  data_length[FORCE_ELEMENT_EDGE_FACTOR] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_FACTOR] = 0;
  data_class[FORCE_ELEMENT_EDGE_FACTOR] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_FACTOR] = FORCE_ELEMENT_EDGE;

  strcpy(name[FORCE_ELEMENT_EDGE_GEOMETRY],"force_element_edge_geometry");
  type[FORCE_ELEMENT_EDGE_GEOMETRY] = INTEGER;
  data_length[FORCE_ELEMENT_EDGE_GEOMETRY] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_GEOMETRY] = 0;
  data_class[FORCE_ELEMENT_EDGE_GEOMETRY] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_GEOMETRY] = FORCE_ELEMENT_EDGE;

  strcpy(name[FORCE_ELEMENT_EDGE_SINE],"force_element_edge_sine");
  type[FORCE_ELEMENT_EDGE_SINE] = DOUBLE_PRECISION;
  data_length[FORCE_ELEMENT_EDGE_SINE] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_SINE] = 0;
  data_class[FORCE_ELEMENT_EDGE_SINE] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_SINE] = FORCE_ELEMENT_EDGE;

  strcpy(name[FORCE_ELEMENT_EDGE_TIME],"force_element_edge_time");
  type[FORCE_ELEMENT_EDGE_TIME] = DOUBLE_PRECISION;
  data_length[FORCE_ELEMENT_EDGE_TIME] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_TIME] = 0;
  data_class[FORCE_ELEMENT_EDGE_TIME] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_TIME] = FORCE_ELEMENT_EDGE;

  strcpy(name[FORCE_ELEMENT_EDGE_TIME_FILE],"force_element_edge_time_file");
  type[FORCE_ELEMENT_EDGE_TIME_FILE] = INTEGER;
  data_length[FORCE_ELEMENT_EDGE_TIME_FILE] = 1;
  data_class[FORCE_ELEMENT_EDGE_TIME_FILE] = FORCE_ELEMENT_EDGE;

  strcpy(name[FORCE_ELEMENT_EDGE_NORMAL],"force_element_edge_normal");
  type[FORCE_ELEMENT_EDGE_NORMAL] = DOUBLE_PRECISION;
  data_length[FORCE_ELEMENT_EDGE_NORMAL] = 1;
  data_class[FORCE_ELEMENT_EDGE_NORMAL] = FORCE;

  strcpy(name[FORCE_ELEMENT_EDGE_NORMAL_FACTOR],"force_element_edge_normal_factor");
  type[FORCE_ELEMENT_EDGE_NORMAL_FACTOR] = DOUBLE_PRECISION;
  data_length[FORCE_ELEMENT_EDGE_NORMAL_FACTOR] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_NORMAL_FACTOR] = 0;
  data_class[FORCE_ELEMENT_EDGE_NORMAL_FACTOR] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_NORMAL_FACTOR] = FORCE_ELEMENT_EDGE_NORMAL;

  strcpy(name[FORCE_ELEMENT_EDGE_NORMAL_GEOMETRY],"force_element_edge_normal_geometry");
  type[FORCE_ELEMENT_EDGE_NORMAL_GEOMETRY] = INTEGER;
  data_length[FORCE_ELEMENT_EDGE_NORMAL_GEOMETRY] = 2;
  fixed_length[FORCE_ELEMENT_EDGE_NORMAL_GEOMETRY] = 0;
  data_class[FORCE_ELEMENT_EDGE_NORMAL_GEOMETRY] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_NORMAL_GEOMETRY] = FORCE_ELEMENT_EDGE_NORMAL;

  strcpy(name[FORCE_ELEMENT_EDGE_NORMAL_SINE],"force_element_edge_normal_sine");
  type[FORCE_ELEMENT_EDGE_NORMAL_SINE] = DOUBLE_PRECISION;
  data_length[FORCE_ELEMENT_EDGE_NORMAL_SINE] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_NORMAL_SINE] = 0;
  data_class[FORCE_ELEMENT_EDGE_NORMAL_SINE] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_NORMAL_SINE] = FORCE_ELEMENT_EDGE_NORMAL;

  strcpy(name[FORCE_ELEMENT_EDGE_NORMAL_TIME],"force_element_edge_normal_time");
  type[FORCE_ELEMENT_EDGE_NORMAL_TIME] = DOUBLE_PRECISION;
  data_length[FORCE_ELEMENT_EDGE_NORMAL_TIME] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_NORMAL_TIME] = 0;
  data_class[FORCE_ELEMENT_EDGE_NORMAL_TIME] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_NORMAL_TIME] = FORCE_ELEMENT_EDGE_NORMAL;

  strcpy(name[FORCE_ELEMENT_EDGE_WATER],"force_element_edge_water");
  type[FORCE_ELEMENT_EDGE_WATER] = DOUBLE_PRECISION;
  if      ( ndim==2 )
    data_length[FORCE_ELEMENT_EDGE_WATER] = 4;
  else if ( ndim==3 )
    data_length[FORCE_ELEMENT_EDGE_WATER] = 5;
  data_class[FORCE_ELEMENT_EDGE_WATER] = FORCE;

  strcpy(name[FORCE_ELEMENT_EDGE_WATER_GEOMETRY],"force_element_edge_water_geometry");
  type[FORCE_ELEMENT_EDGE_WATER_GEOMETRY] = INTEGER;
  data_length[FORCE_ELEMENT_EDGE_WATER_GEOMETRY] = 2;
  fixed_length[FORCE_ELEMENT_EDGE_WATER_GEOMETRY] = 0;
  data_class[FORCE_ELEMENT_EDGE_WATER_GEOMETRY] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_WATER_GEOMETRY] = FORCE_ELEMENT_EDGE_WATER;

  strcpy(name[FORCE_ELEMENT_EDGE_WATER_TIME],"force_element_edge_water_time");
  type[FORCE_ELEMENT_EDGE_WATER_TIME] = DOUBLE_PRECISION;
  data_length[FORCE_ELEMENT_EDGE_WATER_TIME] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_WATER_TIME] = 0;
  data_class[FORCE_ELEMENT_EDGE_WATER_TIME] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_WATER_TIME] = FORCE_ELEMENT_EDGE_WATER;

  strcpy(name[FORCE_ELEMENT_VOLUME],"force_element_volume");
  type[FORCE_ELEMENT_VOLUME] = DOUBLE_PRECISION;
  data_length[FORCE_ELEMENT_VOLUME] = nprinc;
  version_all[FORCE_ELEMENT_VOLUME] = 1;
  data_class[FORCE_ELEMENT_VOLUME] = FORCE;

  strcpy(name[FORCE_ELEMENT_VOLUME_FACTOR],"force_element_volume_factor");
  type[FORCE_ELEMENT_VOLUME_FACTOR] = DOUBLE_PRECISION;
  data_length[FORCE_ELEMENT_VOLUME_FACTOR] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_VOLUME_FACTOR] = 0;
  data_class[FORCE_ELEMENT_VOLUME_FACTOR] = FORCE;
  data_required[FORCE_ELEMENT_VOLUME_FACTOR] = FORCE_ELEMENT_VOLUME;

  strcpy(name[FORCE_ELEMENT_VOLUME_GEOMETRY],"force_element_volume_geometry");
  type[FORCE_ELEMENT_VOLUME_GEOMETRY] = INTEGER;
  data_length[FORCE_ELEMENT_VOLUME_GEOMETRY] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_VOLUME_GEOMETRY] = 0;
  data_class[FORCE_ELEMENT_VOLUME_GEOMETRY] = FORCE;
  data_required[FORCE_ELEMENT_VOLUME_GEOMETRY] = FORCE_ELEMENT_VOLUME;

  strcpy(name[FORCE_ELEMENT_VOLUME_SINE],"force_element_volume_sine");
  type[FORCE_ELEMENT_VOLUME_SINE] = DOUBLE_PRECISION;
  data_length[FORCE_ELEMENT_VOLUME_SINE] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_VOLUME_SINE] = 0;
  data_class[FORCE_ELEMENT_VOLUME_SINE] = FORCE;
  data_required[FORCE_ELEMENT_VOLUME_SINE] = FORCE_ELEMENT_VOLUME;

  strcpy(name[FORCE_ELEMENT_VOLUME_TIME],"force_element_volume_time");
  type[FORCE_ELEMENT_VOLUME_TIME] = DOUBLE_PRECISION;
  data_length[FORCE_ELEMENT_VOLUME_TIME] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_VOLUME_TIME] = 0;
  data_class[FORCE_ELEMENT_VOLUME_TIME] = FORCE;
  data_required[FORCE_ELEMENT_VOLUME_TIME] = FORCE_ELEMENT_VOLUME;

  strcpy(name[FORCE_GRAVITY],"force_gravity");
  type[FORCE_GRAVITY] = DOUBLE_PRECISION;
  data_length[FORCE_GRAVITY] = ndim;
  no_index[FORCE_GRAVITY] = 1;
  data_class[FORCE_GRAVITY] = FORCE_GRAVITY;

  strcpy(name[FORCE_GRAVITY_TIME],"force_gravity_time");
  type[FORCE_GRAVITY_TIME] = DOUBLE_PRECISION;
  no_index[FORCE_GRAVITY_TIME] = 1;
  data_length[FORCE_GRAVITY_TIME] = DATA_ITEM_SIZE;
  fixed_length[FORCE_GRAVITY_TIME] = 0;
  data_class[FORCE_GRAVITY_TIME] = FORCE_GRAVITY;
  data_required[FORCE_GRAVITY_TIME] = FORCE_GRAVITY;

  strcpy(name[FROM],"from");

  strcpy(name[FRONT],"front");

  strcpy(name[GAUSS],"gauss");

  strcpy(name[GEOMETRY],"geometry");

  strcpy(name[GEOMETRY_BOUNDA_FACTOR],"geometry_bounda_factor");
  type[GEOMETRY_BOUNDA_FACTOR] = DOUBLE_PRECISION;
  data_length[GEOMETRY_BOUNDA_FACTOR] = 3;
  fixed_length[GEOMETRY_BOUNDA_FACTOR] = 0;
  data_class[GEOMETRY_BOUNDA_FACTOR] = GEOMETRY;

  strcpy(name[GEOMETRY_BOUNDA_SINE_X],"geometry_bounda_sine_x");
  type[GEOMETRY_BOUNDA_SINE_X] = DOUBLE_PRECISION;
  data_length[GEOMETRY_BOUNDA_SINE_X] = 2;
  data_class[GEOMETRY_BOUNDA_SINE_X] = GEOMETRY;

  strcpy(name[GEOMETRY_BOUNDA_SINE_Y],"geometry_bounda_sine_y");
  type[GEOMETRY_BOUNDA_SINE_Y] = DOUBLE_PRECISION;
  data_length[GEOMETRY_BOUNDA_SINE_Y] = 2;
  data_class[GEOMETRY_BOUNDA_SINE_Y] = GEOMETRY;

  strcpy(name[GEOMETRY_BOUNDA_SINE_Z],"geometry_bounda_sine_z");
  type[GEOMETRY_BOUNDA_SINE_Z] = DOUBLE_PRECISION;
  data_length[GEOMETRY_BOUNDA_SINE_Z] = 2;
  data_class[GEOMETRY_BOUNDA_SINE_Z] = GEOMETRY;

  strcpy(name[GEOMETRY_BRICK],"geometry_brick");
  type[GEOMETRY_BRICK] = DOUBLE_PRECISION;
  data_length[GEOMETRY_BRICK] = 2*MDIM+1;
  data_class[GEOMETRY_BRICK] = GEOMETRY;

  strcpy(name[GEOMETRY_CIRCLE],"geometry_circle");
  type[GEOMETRY_CIRCLE] = DOUBLE_PRECISION;
  data_length[GEOMETRY_CIRCLE] = ndim+2;
  data_class[GEOMETRY_CIRCLE] = GEOMETRY;

  strcpy(name[GEOMETRY_CIRCLE_SEGMENT],"geometry_circle_segment");
  type[GEOMETRY_CIRCLE_SEGMENT] = DOUBLE_PRECISION;
  data_length[GEOMETRY_CIRCLE_SEGMENT] = ndim+1+ndim+1;
  data_class[GEOMETRY_CIRCLE_SEGMENT] = GEOMETRY;         

  strcpy(name[GEOMETRY_CIRCLE_SMALLSEGMENT],"geometry_circle_smallsegment");
  type[GEOMETRY_CIRCLE_SMALLSEGMENT] = DOUBLE_PRECISION;
  data_length[GEOMETRY_CIRCLE_SMALLSEGMENT] = ndim+1+2*ndim+1;
  data_class[GEOMETRY_CIRCLE_SMALLSEGMENT] = GEOMETRY;         

  strcpy(name[GEOMETRY_ELLIPSE],"geometry_ellipse");
  type[GEOMETRY_ELLIPSE] = DOUBLE_PRECISION;
  data_length[GEOMETRY_ELLIPSE] = ndim+3;
  data_class[GEOMETRY_ELLIPSE] = GEOMETRY;

  strcpy(name[GEOMETRY_CYLINDER],"geometry_cylinder");
  type[GEOMETRY_CYLINDER] = DOUBLE_PRECISION;
  data_length[GEOMETRY_CYLINDER] = 8;
  data_class[GEOMETRY_CYLINDER] = GEOMETRY;

  strcpy(name[GEOMETRY_CYLINDER_SEGMENT],"geometry_cylinder_segment");
  type[GEOMETRY_CYLINDER_SEGMENT] = DOUBLE_PRECISION;
  data_length[GEOMETRY_CYLINDER_SEGMENT] = ndim+ndim+1+ndim+1;
  data_class[GEOMETRY_CYLINDER_SEGMENT] = GEOMETRY;

  strcpy(name[GEOMETRY_LINE],"geometry_line");
  type[GEOMETRY_LINE] = DOUBLE_PRECISION;
  data_length[GEOMETRY_LINE] = 2*ndim+1;
  data_class[GEOMETRY_LINE] = GEOMETRY;

  strcpy(name[GEOMETRY_NUMBER],"geometry_number");
  type[GEOMETRY_NUMBER] = INTEGER;
  data_length[GEOMETRY_NUMBER] = 1;
  data_class[GEOMETRY_NUMBER] = GEOMETRY;

  strcpy(name[GEOMETRY_POINT],"geometry_point");
  type[GEOMETRY_POINT] = DOUBLE_PRECISION;
  data_length[GEOMETRY_POINT] = ndim+1;
  data_class[GEOMETRY_POINT] = GEOMETRY;

  strcpy(name[GEOMETRY_POLYNOMIAL],"geometry_polynomial");
  type[GEOMETRY_POLYNOMIAL] = DOUBLE_PRECISION;
  data_length[GEOMETRY_POLYNOMIAL] = DATA_ITEM_SIZE;
  fixed_length[GEOMETRY_POLYNOMIAL] = 0;
  data_class[GEOMETRY_POLYNOMIAL] = GEOMETRY;

  strcpy(name[GEOMETRY_QUADRILATERAL],"geometry_quadrilateral");
  type[GEOMETRY_QUADRILATERAL] = DOUBLE_PRECISION;
  data_length[GEOMETRY_QUADRILATERAL] = 4*ndim+1;
  data_class[GEOMETRY_QUADRILATERAL] = GEOMETRY;

  strcpy(name[GEOMETRY_SET],"geometry_set");
  type[GEOMETRY_SET] = INTEGER;
  data_length[GEOMETRY_SET] = DATA_ITEM_SIZE;
  fixed_length[GEOMETRY_SET] = 0;
  data_class[GEOMETRY_SET] = GEOMETRY;

  strcpy(name[GEOMETRY_SPHERE],"geometry_sphere");
  type[GEOMETRY_SPHERE] = DOUBLE_PRECISION;
  data_length[GEOMETRY_SPHERE] = ndim+2;
  data_class[GEOMETRY_SPHERE] = GEOMETRY;

  strcpy(name[GEOMETRY_SPHERE_SEGMENT],"geometry_sphere_segment");
  type[GEOMETRY_SPHERE_SEGMENT] = DOUBLE_PRECISION;
  data_length[GEOMETRY_SPHERE_SEGMENT] = ndim+1+ndim+1;
  data_class[GEOMETRY_SPHERE_SEGMENT] = GEOMETRY;

  strcpy(name[GEOMETRY_TRIANGLE],"geometry_triangle");
  type[GEOMETRY_TRIANGLE] = DOUBLE_PRECISION;
  data_length[GEOMETRY_TRIANGLE] = 3*ndim+1;
  data_class[GEOMETRY_TRIANGLE] = GEOMETRY;

  strcpy(name[GEOMETRY_TRIANGLE_EPSISO],"geometry_triangle_epsiso");
  type[GEOMETRY_TRIANGLE_EPSISO] = DOUBLE_PRECISION;
  data_length[GEOMETRY_TRIANGLE_EPSISO] = 1;
  data_class[GEOMETRY_TRIANGLE_EPSISO] = GEOMETRY;

  strcpy(name[GENERALIZED],"generalized");

  strcpy(name[GET],"get");

  strcpy(name[GLOBAL_ELEMENTS],"global_elements");
  type[GLOBAL_ELEMENTS] = INTEGER;
  data_length[GLOBAL_ELEMENTS] = 1;
  data_class[GLOBAL_ELEMENTS] = GLOBAL_ELEMENTS;
  no_index[GLOBAL_ELEMENTS] = 1;
  print_only[GLOBAL_ELEMENTS] = 1;         

  strcpy(name[GLOBAL_MASS],"global_mass");
  type[GLOBAL_MASS] = DOUBLE_PRECISION;
  data_length[GLOBAL_MASS] = 1;
  data_class[GLOBAL_MASS] = GLOBAL_MASS;
  no_index[GLOBAL_MASS] = 1;
  print_only[GLOBAL_MASS] = 1;            

  strcpy(name[GLOBAL_NODES],"global_nodes");
  type[GLOBAL_NODES] = INTEGER;
  data_length[GLOBAL_NODES] = 1;
  data_class[GLOBAL_NODES] = GLOBAL_NODES;
  no_index[GLOBAL_NODES] = 1;
  print_only[GLOBAL_NODES] = 1;

  strcpy(name[GLOBAL_POINT_MATERI_DIFFUSION_LOST],"global_point_materi_diffusion_lost");
  type[GLOBAL_POINT_MATERI_DIFFUSION_LOST] = INTEGER;
  data_length[GLOBAL_POINT_MATERI_DIFFUSION_LOST] = 1;
  no_index[GLOBAL_POINT_MATERI_DIFFUSION_LOST] = 1;

  strcpy(name[GLOBAL_POINT_MATERI_DIFFUSION_TOTAL],"global_point_materi_diffusion_total");
  type[GLOBAL_POINT_MATERI_DIFFUSION_TOTAL] = INTEGER;
  data_length[GLOBAL_POINT_MATERI_DIFFUSION_TOTAL] = 1;
  no_index[GLOBAL_POINT_MATERI_DIFFUSION_TOTAL] = 1;

  strcpy(name[GLOBAL_SOLVER_ITERATIONS],"global_solver_iterations");
  type[GLOBAL_SOLVER_ITERATIONS] = INTEGER;
  data_length[GLOBAL_SOLVER_ITERATIONS] = 1;
  data_class[GLOBAL_SOLVER_ITERATIONS] = GLOBAL_SOLVER_ITERATIONS;
  no_index[GLOBAL_SOLVER_ITERATIONS] = 1;
  print_only[GLOBAL_SOLVER_ITERATIONS] = 1;

  strcpy(name[GLOBAL_SOLVER_ERROR],"global_solver_error");
  type[GLOBAL_SOLVER_ERROR] = DOUBLE_PRECISION;
  data_length[GLOBAL_SOLVER_ERROR] = 1;
  data_class[GLOBAL_SOLVER_ERROR] = GLOBAL_SOLVER_ERROR;
  no_index[GLOBAL_SOLVER_ERROR] = 1;
  print_only[GLOBAL_SOLVER_ERROR] = 1;

  strcpy(name[GLOBAL_STRAINENERGY],"global_strainenergy");
  type[GLOBAL_STRAINENERGY] = DOUBLE_PRECISION;
  data_length[GLOBAL_STRAINENERGY] = 1;
  data_class[GLOBAL_STRAINENERGY] = GLOBAL_STRAINENERGY;
  no_index[GLOBAL_STRAINENERGY] = 1;
  print_only[GLOBAL_STRAINENERGY] = 1;            

  strcpy(name[GLOBAL_UNKNOWN_AVERAGE],"global_unknown_average");
  type[GLOBAL_UNKNOWN_AVERAGE] = DOUBLE_PRECISION;
  data_length[GLOBAL_UNKNOWN_AVERAGE] = nuknwn;
  data_class[GLOBAL_UNKNOWN_AVERAGE] = GLOBAL_UNKNOWN_AVERAGE;
  no_index[GLOBAL_UNKNOWN_AVERAGE] = 1;
  print_only[GLOBAL_UNKNOWN_AVERAGE] = 1;

  strcpy(name[GLOBAL_UNKNOWN_MIN],"global_unknown_min");
  type[GLOBAL_UNKNOWN_MIN] = DOUBLE_PRECISION;
  data_length[GLOBAL_UNKNOWN_MIN] = nuknwn;
  data_class[GLOBAL_UNKNOWN_MIN] = GLOBAL_UNKNOWN_MIN;
  no_index[GLOBAL_UNKNOWN_MIN] = 1;
  print_only[GLOBAL_UNKNOWN_MIN] = 1;

  strcpy(name[GLOBAL_UNKNOWN_MAX],"global_unknown_max");
  type[GLOBAL_UNKNOWN_MAX] = DOUBLE_PRECISION;
  data_length[GLOBAL_UNKNOWN_MAX] = nuknwn;
  data_class[GLOBAL_UNKNOWN_MAX] = GLOBAL_UNKNOWN_MAX;
  no_index[GLOBAL_UNKNOWN_MAX] = 1;
  print_only[GLOBAL_UNKNOWN_MAX] = 1;

  strcpy(name[GLOBAL_UNKNOWN_NUMBER],"global_unknown_number");
  type[GLOBAL_UNKNOWN_NUMBER] = INTEGER;
  data_length[GLOBAL_UNKNOWN_NUMBER] = 1;
  data_class[GLOBAL_UNKNOWN_NUMBER] = GLOBAL_UNKNOWN_NUMBER;
  no_index[GLOBAL_UNKNOWN_NUMBER] = 1;
  print_only[GLOBAL_UNKNOWN_NUMBER] = 1;

  strcpy(name[GLOBAL_UNKNOWN_SUM],"global_unknown_sum");
  type[GLOBAL_UNKNOWN_SUM] = DOUBLE_PRECISION;
  data_length[GLOBAL_UNKNOWN_SUM] = nuknwn;
  data_class[GLOBAL_UNKNOWN_SUM] = GLOBAL_UNKNOWN_SUM;
  no_index[GLOBAL_UNKNOWN_SUM] = 1;
  print_only[GLOBAL_UNKNOWN_SUM] = 1;

  strcpy(name[GLOBAL_VOLUME],"global_volume");
  type[GLOBAL_VOLUME] = DOUBLE_PRECISION;
  data_length[GLOBAL_VOLUME] = 1;
  data_class[GLOBAL_VOLUME] = GLOBAL_VOLUME;
  no_index[GLOBAL_VOLUME] = 1;
  print_only[GLOBAL_VOLUME] = 1;

  strcpy(name[GMRES],"gmres");

  strcpy(name[GROUND],"ground");

  strcpy(name[GROUNDFLOW],"groundflow");

  strcpy(name[GROUNDFLOW_ADDTOPRESSURE],"groundflow_addtopressure");
  type[GROUNDFLOW_ADDTOPRESSURE] = DOUBLE_PRECISION;
  data_length[GROUNDFLOW_ADDTOPRESSURE] = 1;
  no_index[GROUNDFLOW_ADDTOPRESSURE] = 1;
  data_class[GROUNDFLOW_ADDTOPRESSURE] = GROUNDFLOW;

  strcpy(name[GROUNDFLOW_DENSITY],"groundflow_density");
  type[GROUNDFLOW_DENSITY] = DOUBLE_PRECISION;
  data_length[GROUNDFLOW_DENSITY] = 1;
  no_index[GROUNDFLOW_DENSITY] = 1;
  data_class[GROUNDFLOW_DENSITY] = GROUNDFLOW;

  strcpy(name[GROUNDFLOW_PHREATICLEVEL],"groundflow_phreaticlevel");
  type[GROUNDFLOW_PHREATICLEVEL] = DOUBLE_PRECISION;
  data_length[GROUNDFLOW_PHREATICLEVEL] = DATA_ITEM_SIZE;
  fixed_length[GROUNDFLOW_PHREATICLEVEL] = 0;
  no_index[GROUNDFLOW_PHREATICLEVEL] = 1;
  data_class[GROUNDFLOW_PHREATICLEVEL] = GROUNDFLOW;

  strcpy(name[GROUNDFLOW_PHREATICLEVEL_N],"groundflow_phreaticlevel_n");
  type[GROUNDFLOW_PHREATICLEVEL_N] = INTEGER;
  data_length[GROUNDFLOW_PHREATICLEVEL_N] = 2;
  no_index[GROUNDFLOW_PHREATICLEVEL_N] = 1;
  data_class[GROUNDFLOW_PHREATICLEVEL_N] = GROUNDFLOW;           
  data_required[GROUNDFLOW_PHREATICLEVEL_N] = GROUNDFLOW_PHREATICLEVEL;           

  strcpy(name[GROUNDFLOW_PHREATICLEVEL_BOUNDA],"groundflow_phreaticlevel_bounda");
  type[GROUNDFLOW_PHREATICLEVEL_BOUNDA] = INTEGER;
  data_length[GROUNDFLOW_PHREATICLEVEL_BOUNDA] = 1;
  no_index[GROUNDFLOW_PHREATICLEVEL_BOUNDA] = 1;
  data_class[GROUNDFLOW_PHREATICLEVEL_BOUNDA] = GROUNDFLOW;           
  data_required[GROUNDFLOW_PHREATICLEVEL_BOUNDA] = GROUNDFLOW_PHREATICLEVEL;           

  strcpy(name[GROUNDFLOW_PHREATICLEVEL_MINIMUM],"groundflow_phreaticlevel_minimum");
  type[GROUNDFLOW_PHREATICLEVEL_MINIMUM] = DOUBLE_PRECISION;
  data_length[GROUNDFLOW_PHREATICLEVEL_MINIMUM] = 1;
  no_index[GROUNDFLOW_PHREATICLEVEL_MINIMUM] = 1;
  data_class[GROUNDFLOW_PHREATICLEVEL_MINIMUM] = GROUNDFLOW;           

  strcpy(name[GROUNDFLOW_PRESSURE],"groundflow_pressure");

  strcpy(name[GROUNDFLOW_PRESSURE_ATMOSPHERIC],"groundflow_pressure_atmospheric");
  type[GROUNDFLOW_PRESSURE_ATMOSPHERIC] = DOUBLE_PRECISION;
  data_length[GROUNDFLOW_PRESSURE_ATMOSPHERIC] = 1;
  no_index[GROUNDFLOW_PRESSURE_ATMOSPHERIC] = 1;
  data_class[GROUNDFLOW_PRESSURE_ATMOSPHERIC] = GROUNDFLOW;

  strcpy(name[GROUNDFLOW_VELOCITY],"groundflow_velocity");

  strcpy(name[GROUP_AXISYMMETRIC],"group_axisymmetric");
  type[GROUP_AXISYMMETRIC] = INTEGER;
  data_length[GROUP_AXISYMMETRIC] = 1;
  data_class[GROUP_AXISYMMETRIC] = GROUP_AXISYMMETRIC;
  data_required[GROUP_AXISYMMETRIC] = GROUP_TYPE;

  strcpy(name[GROUP_BEAM_INERTIA],"group_beam_inertia");
  type[GROUP_BEAM_INERTIA] = DOUBLE_PRECISION;
  data_length[GROUP_BEAM_INERTIA] = 1;
  data_class[GROUP_BEAM_INERTIA] = BEAM;
  data_required[GROUP_BEAM_INERTIA] = GROUP_TYPE;

  strcpy(name[GROUP_BEAM_MEMORY],"group_beam_memory");
  type[GROUP_BEAM_MEMORY] = INTEGER;
  data_length[GROUP_BEAM_MEMORY] = 1;
  data_class[GROUP_BEAM_MEMORY] = BEAM;
  data_required[GROUP_BEAM_MEMORY] = GROUP_TYPE;

  strcpy(name[GROUP_BEAM_PLANE],"group_beam_plane");
  type[GROUP_BEAM_PLANE] = INTEGER;
  data_length[GROUP_BEAM_PLANE] = 2;
  data_class[GROUP_BEAM_PLANE] = BEAM;
  data_required[GROUP_BEAM_PLANE] = GROUP_TYPE;

  strcpy(name[GROUP_BEAM_YOUNG],"group_beam_young");
  type[GROUP_BEAM_YOUNG] = DOUBLE_PRECISION;
  data_length[GROUP_BEAM_YOUNG] = 1;
  data_class[GROUP_BEAM_YOUNG] = BEAM;
  data_required[GROUP_BEAM_YOUNG] = GROUP_TYPE;

  strcpy(name[GROUP_CONDIF_ABSORPTION],"group_condif_absorption");
  type[GROUP_CONDIF_ABSORPTION] = DOUBLE_PRECISION;
  data_length[GROUP_CONDIF_ABSORPTION] = 1;
  data_class[GROUP_CONDIF_ABSORPTION] = CONDIF;
  data_required[GROUP_CONDIF_ABSORPTION] = GROUP_TYPE;

  strcpy(name[GROUP_CONDIF_DENSITY],"group_condif_density");
  type[GROUP_CONDIF_DENSITY] = DOUBLE_PRECISION;
  data_length[GROUP_CONDIF_DENSITY] = 1;
  data_class[GROUP_CONDIF_DENSITY] = CONDIF;
  data_required[GROUP_CONDIF_DENSITY] = GROUP_TYPE;

  strcpy(name[GROUP_CONDIF_CAPACITY],"group_condif_capacity");
  type[GROUP_CONDIF_CAPACITY] = DOUBLE_PRECISION;
  data_length[GROUP_CONDIF_CAPACITY] = 1;
  data_class[GROUP_CONDIF_CAPACITY] = CONDIF;
  data_required[GROUP_CONDIF_CAPACITY] = GROUP_TYPE;

  strcpy(name[GROUP_CONDIF_FLOW],"group_condif_flow");
  type[GROUP_CONDIF_FLOW] = DOUBLE_PRECISION;
  data_length[GROUP_CONDIF_FLOW] = ndim;
  data_class[GROUP_CONDIF_FLOW] = CONDIF;
  data_required[GROUP_CONDIF_FLOW] = GROUP_TYPE;

  strcpy(name[GROUP_CONDIF_CONDUCTIVITY],"group_condif_conductivity");
  type[GROUP_CONDIF_CONDUCTIVITY] = DOUBLE_PRECISION;
  data_length[GROUP_CONDIF_CONDUCTIVITY] = 1;
  data_class[GROUP_CONDIF_CONDUCTIVITY] = CONDIF;
  data_required[GROUP_CONDIF_CONDUCTIVITY] = GROUP_TYPE;

  strcpy(name[GROUP_CONTACTSPRING_COHESION],"group_contactspring_cohesion");
  type[GROUP_CONTACTSPRING_COHESION] = DOUBLE_PRECISION;
  data_length[GROUP_CONTACTSPRING_COHESION] = 1;
  data_class[GROUP_CONTACTSPRING_COHESION] = CONTACTSPRING;
  data_required[GROUP_CONTACTSPRING_COHESION] = GROUP_TYPE;

  strcpy(name[GROUP_CONTACTSPRING_DIRECTION],"group_contactspring_direction");
  type[GROUP_CONTACTSPRING_DIRECTION] = DOUBLE_PRECISION;
  data_length[GROUP_CONTACTSPRING_DIRECTION] = MDIM;
  data_class[GROUP_CONTACTSPRING_DIRECTION] = CONTACTSPRING;
  data_required[GROUP_CONTACTSPRING_DIRECTION] = GROUP_TYPE;

  strcpy(name[GROUP_CONTACTSPRING_FRICTION],"group_contactspring_friction");
  type[GROUP_CONTACTSPRING_FRICTION] = DOUBLE_PRECISION;
  data_length[GROUP_CONTACTSPRING_FRICTION] = 1;
  data_class[GROUP_CONTACTSPRING_FRICTION] = CONTACTSPRING;
  data_required[GROUP_CONTACTSPRING_FRICTION] = GROUP_TYPE;

  strcpy(name[GROUP_CONTACTSPRING_FRICTION_AUTOMATIC],"group_contactspring_friction_automatic");
  type[GROUP_CONTACTSPRING_FRICTION_AUTOMATIC] = INTEGER;
  data_length[GROUP_CONTACTSPRING_FRICTION_AUTOMATIC] = 1;
  data_class[GROUP_CONTACTSPRING_FRICTION_AUTOMATIC] = CONTACTSPRING;
  data_required[GROUP_CONTACTSPRING_FRICTION_AUTOMATIC] = GROUP_TYPE;

  strcpy(name[GROUP_CONTACTSPRING_MEMORY],"group_contactspring_memory");
  type[GROUP_CONTACTSPRING_MEMORY] = INTEGER;
  data_length[GROUP_CONTACTSPRING_MEMORY] = 1;
  data_class[GROUP_CONTACTSPRING_MEMORY] = CONTACTSPRING;
  data_required[GROUP_CONTACTSPRING_MEMORY] = GROUP_TYPE;

  strcpy(name[GROUP_CONTACTSPRING_STIFFNESS],"group_contactspring_stiffness");
  type[GROUP_CONTACTSPRING_STIFFNESS] = DOUBLE_PRECISION;
  data_length[GROUP_CONTACTSPRING_STIFFNESS] = 2;
  data_class[GROUP_CONTACTSPRING_STIFFNESS] = CONTACTSPRING;
  data_required[GROUP_CONTACTSPRING_STIFFNESS] = GROUP_TYPE;

  strcpy(name[GROUP_GROUNDFLOW_CAPACITY],"group_groundflow_capacity");
  type[GROUP_GROUNDFLOW_CAPACITY] = DOUBLE_PRECISION;
  data_length[GROUP_GROUNDFLOW_CAPACITY] = 1;
  data_class[GROUP_GROUNDFLOW_CAPACITY] = GROUNDFLOW;
  data_required[GROUP_GROUNDFLOW_CAPACITY] = GROUP_TYPE;

  strcpy(name[GROUP_GROUNDFLOW_CAPACITY_NONLINEAR_METHOD],"group_groundflow_capacity_nonlinear_method");
  type[GROUP_GROUNDFLOW_CAPACITY_NONLINEAR_METHOD] = INTEGER;
  data_length[GROUP_GROUNDFLOW_CAPACITY_NONLINEAR_METHOD] = 1;
  data_class[GROUP_GROUNDFLOW_CAPACITY_NONLINEAR_METHOD] = GROUNDFLOW;
  data_required[GROUP_GROUNDFLOW_CAPACITY_NONLINEAR_METHOD] = GROUP_TYPE;

  strcpy(name[GROUP_GROUNDFLOW_CAPACITY_NONLINEAR_PARAMETERS],"group_groundflow_capacity_nonlinear_parameters");
  type[GROUP_GROUNDFLOW_CAPACITY_NONLINEAR_PARAMETERS] = DOUBLE_PRECISION;
  data_length[GROUP_GROUNDFLOW_CAPACITY_NONLINEAR_PARAMETERS] = DATA_ITEM_SIZE;
  data_class[GROUP_GROUNDFLOW_CAPACITY_NONLINEAR_PARAMETERS] = GROUNDFLOW;
  data_required[GROUP_GROUNDFLOW_CAPACITY_NONLINEAR_PARAMETERS] = GROUP_TYPE;
  fixed_length[GROUP_GROUNDFLOW_CAPACITY_NONLINEAR_PARAMETERS] = 0;

  strcpy(name[GROUP_GROUNDFLOW_MATERIDIVERGENCE],"group_groundflow_materidivergence");
  type[GROUP_GROUNDFLOW_MATERIDIVERGENCE] = INTEGER;
  data_length[GROUP_GROUNDFLOW_MATERIDIVERGENCE] = 1;
  data_class[GROUP_GROUNDFLOW_MATERIDIVERGENCE] = GROUNDFLOW;
  data_required[GROUP_GROUNDFLOW_MATERIDIVERGENCE] = GROUP_TYPE;

  strcpy(name[GROUP_GROUNDFLOW_PERMEABILITY],"group_groundflow_permeability");
  type[GROUP_GROUNDFLOW_PERMEABILITY] = DOUBLE_PRECISION;
  data_length[GROUP_GROUNDFLOW_PERMEABILITY] = ndim;
  data_class[GROUP_GROUNDFLOW_PERMEABILITY] = GROUNDFLOW;
  data_required[GROUP_GROUNDFLOW_PERMEABILITY] = GROUP_TYPE;

  strcpy(name[GROUP_GROUNDFLOW_PERMEABILITY_NONLINEAR_METHOD],"group_groundflow_permeability_nonlinear_method");
  type[GROUP_GROUNDFLOW_PERMEABILITY_NONLINEAR_METHOD] = INTEGER;
  data_length[GROUP_GROUNDFLOW_PERMEABILITY_NONLINEAR_METHOD] = 1;
  data_class[GROUP_GROUNDFLOW_PERMEABILITY_NONLINEAR_METHOD] = GROUNDFLOW;
  data_required[GROUP_GROUNDFLOW_PERMEABILITY_NONLINEAR_METHOD] = GROUP_TYPE;

  strcpy(name[GROUP_GROUNDFLOW_PERMEABILITY_NONLINEAR_PARAMETERS],"group_groundflow_permeability_nonlinear_parameters");
  type[GROUP_GROUNDFLOW_PERMEABILITY_NONLINEAR_PARAMETERS] = DOUBLE_PRECISION;
  data_length[GROUP_GROUNDFLOW_PERMEABILITY_NONLINEAR_PARAMETERS] = DATA_ITEM_SIZE;
  data_class[GROUP_GROUNDFLOW_PERMEABILITY_NONLINEAR_PARAMETERS] = GROUNDFLOW;
  data_required[GROUP_GROUNDFLOW_PERMEABILITY_NONLINEAR_PARAMETERS] = GROUP_TYPE;
  fixed_length[GROUP_GROUNDFLOW_PERMEABILITY_NONLINEAR_PARAMETERS] = 0;

  strcpy(name[GROUP_GROUNDFLOW_POROSITY],"group_groundflow_porosity");
  type[GROUP_GROUNDFLOW_POROSITY] = DOUBLE_PRECISION;
  data_length[GROUP_GROUNDFLOW_POROSITY] = 1;
  data_class[GROUP_GROUNDFLOW_POROSITY] = GROUNDFLOW;
  data_required[GROUP_GROUNDFLOW_POROSITY] = GROUP_TYPE;

  strcpy(name[GROUP_INTEGRATION_METHOD],"group_integration_method");
  type[GROUP_INTEGRATION_METHOD] = INTEGER;
  data_length[GROUP_INTEGRATION_METHOD] = 1;
  version_all[GROUP_INTEGRATION_METHOD] = 1;
  data_class[GROUP_INTEGRATION_METHOD] = GROUP_INTEGRATION_METHOD;  
  data_required[GROUP_INTEGRATION_METHOD] = GROUP_TYPE;  

  strcpy(name[GROUP_INTEGRATION_POINTS],"group_integration_points");
  type[GROUP_INTEGRATION_POINTS] = INTEGER;
  data_length[GROUP_INTEGRATION_POINTS] = 1;
  version_all[GROUP_INTEGRATION_POINTS] = 1;
  data_class[GROUP_INTEGRATION_POINTS] = GROUP_INTEGRATION_POINTS;
  data_required[GROUP_INTEGRATION_POINTS] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_ELASTI_CAMCLAY_G],"group_materi_elasti_camclay_g");
  type[GROUP_MATERI_ELASTI_CAMCLAY_G] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_ELASTI_CAMCLAY_G] = 1;
  data_class[GROUP_MATERI_ELASTI_CAMCLAY_G] = MATERI;
  data_required[GROUP_MATERI_ELASTI_CAMCLAY_G] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_ELASTI_CAMCLAY_POISSON],"group_materi_elasti_camclay_poisson");
  type[GROUP_MATERI_ELASTI_CAMCLAY_POISSON] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_ELASTI_CAMCLAY_POISSON] = 1;
  data_class[GROUP_MATERI_ELASTI_CAMCLAY_POISSON] = MATERI;
  data_required[GROUP_MATERI_ELASTI_CAMCLAY_POISSON] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_ELASTI_TSKH],"group_materi_elasti_tskh");
  type[GROUP_MATERI_ELASTI_TSKH] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_ELASTI_TSKH] = DATA_ITEM_SIZE;
  data_class[GROUP_MATERI_ELASTI_TSKH] = MATERI;
  fixed_length[GROUP_MATERI_ELASTI_TSKH] = 0;
  data_required[GROUP_MATERI_ELASTI_TSKH] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_DAMAGE_MAZARS],"group_materi_damage_mazars");
  type[GROUP_MATERI_DAMAGE_MAZARS] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_DAMAGE_MAZARS] = 6;
  data_required[GROUP_MATERI_DAMAGE_MAZARS] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_DAMPING],"group_materi_damping");
  type[GROUP_MATERI_DAMPING] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_DAMPING] = 1;
  data_required[GROUP_MATERI_DAMPING] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_DENSITY],"group_materi_density");
  type[GROUP_MATERI_DENSITY] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_DENSITY] = 1;
  data_required[GROUP_MATERI_DENSITY] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_DENSITY_GROUNDFLOW],"group_materi_density_groundflow");
  type[GROUP_MATERI_DENSITY_GROUNDFLOW] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_DENSITY_GROUNDFLOW] = 2;
  data_required[GROUP_MATERI_DENSITY_GROUNDFLOW] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_ELASTI_COMPRESSIBILITY],"group_materi_elasti_compressibility");
  type[GROUP_MATERI_ELASTI_COMPRESSIBILITY] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_ELASTI_COMPRESSIBILITY] = 1;
  data_class[GROUP_MATERI_ELASTI_COMPRESSIBILITY] = MATERI;
  data_required[GROUP_MATERI_ELASTI_COMPRESSIBILITY] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_ELASTI_LADE],"group_materi_elasti_lade");
  type[GROUP_MATERI_ELASTI_LADE] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_ELASTI_LADE] = 3;
  data_class[GROUP_MATERI_ELASTI_LADE] = MATERI;
  data_required[GROUP_MATERI_ELASTI_LADE] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_ELASTI_POISSON],"group_materi_elasti_poisson");
  type[GROUP_MATERI_ELASTI_POISSON] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_ELASTI_POISSON] = 1;
  data_class[GROUP_MATERI_ELASTI_POISSON] = MATERI;
  data_required[GROUP_MATERI_ELASTI_POISSON] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_ELASTI_SMALLSTRAIN],"group_materi_elasti_smallstrain");
  type[GROUP_MATERI_ELASTI_SMALLSTRAIN] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_ELASTI_SMALLSTRAIN] = 6;
  data_class[GROUP_MATERI_ELASTI_SMALLSTRAIN] = MATERI;
  data_required[GROUP_MATERI_ELASTI_SMALLSTRAIN] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_ELASTI_TRANSVERSE_ISOTROPY],"group_materi_elasti_transverse_isotropy");
  type[GROUP_MATERI_ELASTI_TRANSVERSE_ISOTROPY] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_ELASTI_TRANSVERSE_ISOTROPY] = 8;
  data_class[GROUP_MATERI_ELASTI_TRANSVERSE_ISOTROPY] = MATERI;
  data_required[GROUP_MATERI_ELASTI_TRANSVERSE_ISOTROPY] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_ELASTI_TRANSVERSE_ISOTROPY_GRAHOUL],"group_materi_elasti_transverse_isotropy_grahoul");
  type[GROUP_MATERI_ELASTI_TRANSVERSE_ISOTROPY_GRAHOUL] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_ELASTI_TRANSVERSE_ISOTROPY_GRAHOUL] = 1;
  data_class[GROUP_MATERI_ELASTI_TRANSVERSE_ISOTROPY_GRAHOUL] = MATERI;
  data_required[GROUP_MATERI_ELASTI_TRANSVERSE_ISOTROPY_GRAHOUL] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_ELASTI_VOLUMETRIC_YOUNG_ORDER],"group_materi_elasti_volumetric_young_order");
  type[GROUP_MATERI_ELASTI_VOLUMETRIC_YOUNG_ORDER] = INTEGER;
  data_length[GROUP_MATERI_ELASTI_VOLUMETRIC_YOUNG_ORDER] = 1;
  data_class[GROUP_MATERI_ELASTI_VOLUMETRIC_YOUNG_ORDER] = MATERI;
  data_required[GROUP_MATERI_ELASTI_VOLUMETRIC_YOUNG_ORDER] = GROUP_TYPE;
 
  strcpy(name[GROUP_MATERI_ELASTI_VOLUMETRIC_POISSON],"group_materi_elasti_volumetric_poisson");
  type[GROUP_MATERI_ELASTI_VOLUMETRIC_POISSON] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_ELASTI_VOLUMETRIC_POISSON] = 1;
  data_class[GROUP_MATERI_ELASTI_VOLUMETRIC_POISSON] = MATERI;
  data_required[GROUP_MATERI_ELASTI_VOLUMETRIC_POISSON] = GROUP_TYPE;
 
  strcpy(name[GROUP_MATERI_ELASTI_VOLUMETRIC_YOUNG_VALUES],"group_materi_elasti_volumetric_young_values");
  type[GROUP_MATERI_ELASTI_VOLUMETRIC_YOUNG_VALUES] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_ELASTI_VOLUMETRIC_YOUNG_VALUES] = DATA_ITEM_SIZE;
  data_class[GROUP_MATERI_ELASTI_VOLUMETRIC_YOUNG_VALUES] = MATERI;
  fixed_length[GROUP_MATERI_ELASTI_VOLUMETRIC_YOUNG_VALUES] = 0;
  data_required[GROUP_MATERI_ELASTI_VOLUMETRIC_YOUNG_VALUES] = GROUP_TYPE;
                                                                           
  strcpy(name[GROUP_MATERI_ELASTI_YOUNG],"group_materi_elasti_young");
  type[GROUP_MATERI_ELASTI_YOUNG] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_ELASTI_YOUNG] = 1;
  data_class[GROUP_MATERI_ELASTI_YOUNG] = MATERI;
  data_required[GROUP_MATERI_ELASTI_YOUNG] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_ELASTI_YOUNG_POLYNOMIAL],"group_materi_elasti_young_polynomial");
  type[GROUP_MATERI_ELASTI_YOUNG_POLYNOMIAL] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_ELASTI_YOUNG_POLYNOMIAL] = DATA_ITEM_SIZE;
  data_class[GROUP_MATERI_ELASTI_YOUNG_POLYNOMIAL] = MATERI;
  fixed_length[GROUP_MATERI_ELASTI_YOUNG_POLYNOMIAL] = 0;
  data_required[GROUP_MATERI_ELASTI_YOUNG_POLYNOMIAL] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_ELASTI_YOUNG_POWER],"group_materi_elasti_young_power");
  type[GROUP_MATERI_ELASTI_YOUNG_POWER] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_ELASTI_YOUNG_POWER] = 3;
  data_class[GROUP_MATERI_ELASTI_YOUNG_POWER] = MATERI;
  data_required[GROUP_MATERI_ELASTI_YOUNG_POWER] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_ELASTI_YOUNG_STRAINSTRESS],"group_materi_elasti_young_strainstress");
  type[GROUP_MATERI_ELASTI_YOUNG_STRAINSTRESS] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_ELASTI_YOUNG_STRAINSTRESS] = DATA_ITEM_SIZE;
  data_class[GROUP_MATERI_ELASTI_YOUNG_STRAINSTRESS] = MATERI;
  fixed_length[GROUP_MATERI_ELASTI_YOUNG_STRAINSTRESS] = 0;
  data_required[GROUP_MATERI_ELASTI_YOUNG_STRAINSTRESS] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_EXPANSION_LINEAR],"group_materi_expansion_linear");
  type[GROUP_MATERI_EXPANSION_LINEAR] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_EXPANSION_LINEAR] = 1;
  data_class[GROUP_MATERI_EXPANSION_LINEAR] = MATERI;
  data_required[GROUP_MATERI_EXPANSION_LINEAR] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_EXPANSION_VOLUME],"group_materi_expansion_volume");
  type[GROUP_MATERI_EXPANSION_VOLUME] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_EXPANSION_VOLUME] = 1;
  data_class[GROUP_MATERI_EXPANSION_VOLUME] = MATERI;
  data_required[GROUP_MATERI_EXPANSION_VOLUME] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_FAILURE_CRUCHING],"group_materi_failure_cruching");
  type[GROUP_MATERI_FAILURE_CRUCHING] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_FAILURE_CRUCHING] = 2;
  data_class[GROUP_MATERI_FAILURE_CRUCHING] = MATERI;
  data_required[GROUP_MATERI_FAILURE_CRUCHING] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_FAILURE_DAMAGE],"group_materi_failure_damage");
  type[GROUP_MATERI_FAILURE_DAMAGE] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_FAILURE_DAMAGE] = 2;
  data_class[GROUP_MATERI_FAILURE_DAMAGE] = MATERI;
  data_required[GROUP_MATERI_FAILURE_DAMAGE] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_FAILURE_PLASTI_KAPPA],"group_materi_failure_plasti_kappa");
  type[GROUP_MATERI_FAILURE_PLASTI_KAPPA] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_FAILURE_PLASTI_KAPPA] = 2;
  data_class[GROUP_MATERI_FAILURE_PLASTI_KAPPA] = MATERI;
  data_required[GROUP_MATERI_FAILURE_PLASTI_KAPPA] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_FAILURE_RUPTURE],"group_materi_failure_rupture");
  type[GROUP_MATERI_FAILURE_RUPTURE] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_FAILURE_RUPTURE] = 2;
  data_class[GROUP_MATERI_FAILURE_RUPTURE] = MATERI;
  data_required[GROUP_MATERI_FAILURE_RUPTURE] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_FAILURE_VOIDFRACTION],"group_materi_failure_voidfraction");
  type[GROUP_MATERI_FAILURE_VOIDFRACTION] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_FAILURE_VOIDFRACTION] = 2;
  data_class[GROUP_MATERI_FAILURE_VOIDFRACTION] = MATERI;
  data_required[GROUP_MATERI_FAILURE_VOIDFRACTION] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_HYPER_BESSELING],"group_materi_hyper_besseling");
  type[GROUP_MATERI_HYPER_BESSELING] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_HYPER_BESSELING] = 3;
  data_class[GROUP_MATERI_HYPER_BESSELING] = MATERI;
  data_required[GROUP_MATERI_HYPER_BESSELING] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_HYPER_BLATZ_KO],"group_materi_hyper_blatz_ko");
  type[GROUP_MATERI_HYPER_BLATZ_KO] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_HYPER_BLATZ_KO] = 2;
  data_class[GROUP_MATERI_HYPER_BLATZ_KO] = MATERI;
  data_required[GROUP_MATERI_HYPER_BLATZ_KO] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_HYPER_MOONEY_RIVLIN],"group_materi_hyper_mooney_rivlin");
  type[GROUP_MATERI_HYPER_MOONEY_RIVLIN] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_HYPER_MOONEY_RIVLIN] = 2;
  data_class[GROUP_MATERI_HYPER_MOONEY_RIVLIN] = MATERI;
  data_required[GROUP_MATERI_HYPER_MOONEY_RIVLIN] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_HYPER_NEOHOOKEAN],"group_materi_hyper_neohookean");
  type[GROUP_MATERI_HYPER_NEOHOOKEAN] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_HYPER_NEOHOOKEAN] = 1;
  data_class[GROUP_MATERI_HYPER_NEOHOOKEAN] = MATERI;
  data_required[GROUP_MATERI_HYPER_NEOHOOKEAN] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_HYPER_REDUCEDPOLYNOMIAL],"group_materi_hyper_reducedpolynomial");
  type[GROUP_MATERI_HYPER_REDUCEDPOLYNOMIAL] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_HYPER_REDUCEDPOLYNOMIAL] = DATA_ITEM_SIZE;
  data_class[GROUP_MATERI_HYPER_REDUCEDPOLYNOMIAL] = MATERI;
  fixed_length[GROUP_MATERI_HYPER_REDUCEDPOLYNOMIAL] = 0;
  data_required[GROUP_MATERI_HYPER_REDUCEDPOLYNOMIAL] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_HYPER_STIFFNESS],"group_materi_hyper_stiffness");
  type[GROUP_MATERI_HYPER_STIFFNESS] = INTEGER;
  data_length[GROUP_MATERI_HYPER_STIFFNESS] = 1;
  data_class[GROUP_MATERI_HYPER_STIFFNESS] = MATERI;
  data_required[GROUP_MATERI_HYPER_STIFFNESS] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_HYPER_VOLUMETRIC_LINEAR],"group_materi_hyper_volumetric_linear");
  type[GROUP_MATERI_HYPER_VOLUMETRIC_LINEAR] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_HYPER_VOLUMETRIC_LINEAR] = 1;
  data_class[GROUP_MATERI_HYPER_VOLUMETRIC_LINEAR] = MATERI;
  data_required[GROUP_MATERI_HYPER_VOLUMETRIC_LINEAR] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_HYPER_VOLUMETRIC_MURNAGHAN],"group_materi_hyper_volumetric_murnaghan");
  type[GROUP_MATERI_HYPER_VOLUMETRIC_MURNAGHAN] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_HYPER_VOLUMETRIC_MURNAGHAN] = 2;
  data_class[GROUP_MATERI_HYPER_VOLUMETRIC_MURNAGHAN] = MATERI;
  data_required[GROUP_MATERI_HYPER_VOLUMETRIC_MURNAGHAN] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_HYPER_VOLUMETRIC_OGDEN],"group_materi_hyper_volumetric_ogden");
  type[GROUP_MATERI_HYPER_VOLUMETRIC_OGDEN] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_HYPER_VOLUMETRIC_OGDEN] = 2;
  data_class[GROUP_MATERI_HYPER_VOLUMETRIC_OGDEN] = MATERI;
  data_required[GROUP_MATERI_HYPER_VOLUMETRIC_OGDEN] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_HYPER_VOLUMETRIC_POLYNOMIAL],"group_materi_hyper_volumetric_polynomial");
  type[GROUP_MATERI_HYPER_VOLUMETRIC_POLYNOMIAL] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_HYPER_VOLUMETRIC_POLYNOMIAL] = DATA_ITEM_SIZE;
  data_class[GROUP_MATERI_HYPER_VOLUMETRIC_POLYNOMIAL] = MATERI;
  fixed_length[GROUP_MATERI_HYPER_VOLUMETRIC_POLYNOMIAL] = 0;
  data_required[GROUP_MATERI_HYPER_VOLUMETRIC_POLYNOMIAL] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_HYPER_VOLUMETRIC_SIMOTAYLOR],"group_materi_hyper_volumetric_simotaylor");
  type[GROUP_MATERI_HYPER_VOLUMETRIC_SIMOTAYLOR] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_HYPER_VOLUMETRIC_SIMOTAYLOR] = 1;
  data_class[GROUP_MATERI_HYPER_VOLUMETRIC_SIMOTAYLOR] = MATERI;
  data_required[GROUP_MATERI_HYPER_VOLUMETRIC_SIMOTAYLOR] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_MAXWELL_CHAIN],"group_materi_maxwell_chain");
  type[GROUP_MATERI_MAXWELL_CHAIN] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_MAXWELL_CHAIN] = materi_maxwell_stress*2;
  data_class[GROUP_MATERI_MAXWELL_CHAIN] = MATERI;
  data_required[GROUP_MATERI_MAXWELL_CHAIN] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_MAXWELL_CHAIN_NONLINEAR],"group_materi_maxwell_chain_nonlinear");
  type[GROUP_MATERI_MAXWELL_CHAIN_NONLINEAR] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_MAXWELL_CHAIN_NONLINEAR] = materi_maxwell_stress;
  data_class[GROUP_MATERI_MAXWELL_CHAIN_NONLINEAR] = MATERI;
  data_required[GROUP_MATERI_MAXWELL_CHAIN_NONLINEAR] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_MEMBRANE],"group_materi_membrane");
  type[GROUP_MATERI_MEMBRANE] = INTEGER;
  data_length[GROUP_MATERI_MEMBRANE] = 1;
  data_class[GROUP_MATERI_MEMBRANE] = MATERI;
  data_required[GROUP_MATERI_MEMBRANE] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_MEMORY],"group_materi_memory");
  type[GROUP_MATERI_MEMORY] = INTEGER;
  data_length[GROUP_MATERI_MEMORY] = 1;
  data_class[GROUP_MATERI_MEMORY] = MATERI;
  data_required[GROUP_MATERI_MEMORY] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_AITSKH],"group_materi_plasti_aitskh");
  type[GROUP_MATERI_PLASTI_AITSKH] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_AITSKH] = DATA_ITEM_SIZE;
  data_class[GROUP_MATERI_PLASTI_AITSKH] = MATERI;
  fixed_length[GROUP_MATERI_PLASTI_AITSKH] = 0;
  data_required[GROUP_MATERI_PLASTI_AITSKH] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_BOUNDARY],"group_materi_plasti_boundary");
  type[GROUP_MATERI_PLASTI_BOUNDARY] = INTEGER;
  data_length[GROUP_MATERI_PLASTI_BOUNDARY] = DATA_ITEM_SIZE;
  data_class[GROUP_MATERI_PLASTI_BOUNDARY] = MATERI;
  fixed_length[GROUP_MATERI_PLASTI_BOUNDARY] = 0;
  data_required[GROUP_MATERI_PLASTI_BOUNDARY] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_BOUNDARY_FACTOR],"group_materi_plasti_boundary_factor");
  type[GROUP_MATERI_PLASTI_BOUNDARY_FACTOR] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_BOUNDARY_FACTOR] = 1;
  data_class[GROUP_MATERI_PLASTI_BOUNDARY_FACTOR] = MATERI;
  data_required[GROUP_MATERI_PLASTI_BOUNDARY_FACTOR] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_CAMCLAY],"group_materi_plasti_camclay");
  type[GROUP_MATERI_PLASTI_CAMCLAY] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_CAMCLAY] = 4;
  data_class[GROUP_MATERI_PLASTI_CAMCLAY] = MATERI;
  data_required[GROUP_MATERI_PLASTI_CAMCLAY] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_CAMCLAY_INCREMENTAL],"group_materi_plasti_camclay_incremental");
  type[GROUP_MATERI_PLASTI_CAMCLAY_INCREMENTAL] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_CAMCLAY_INCREMENTAL] = 4;
  data_class[GROUP_MATERI_PLASTI_CAMCLAY_INCREMENTAL] = MATERI;
  data_required[GROUP_MATERI_PLASTI_CAMCLAY_INCREMENTAL] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_CAP],"group_materi_plasti_cap");
  type[GROUP_MATERI_PLASTI_CAP] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_CAP] = DATA_ITEM_SIZE;
  data_class[GROUP_MATERI_PLASTI_CAP] = MATERI;
  fixed_length[GROUP_MATERI_PLASTI_CAP] = 0;
  data_required[GROUP_MATERI_PLASTI_CAP] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_COMPRESSION],"group_materi_plasti_compression");
  type[GROUP_MATERI_PLASTI_COMPRESSION] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_COMPRESSION] = 1;
  data_class[GROUP_MATERI_PLASTI_COMPRESSION] = MATERI;
  data_required[GROUP_MATERI_PLASTI_COMPRESSION] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_DIPRISCO],"group_materi_plasti_diprisco");
  type[GROUP_MATERI_PLASTI_DIPRISCO] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_DIPRISCO] = 10;
  data_class[GROUP_MATERI_PLASTI_DIPRISCO] = MATERI;
  data_required[GROUP_MATERI_PLASTI_DIPRISCO] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_DIPRISCO_RT],"group_materi_plasti_diprisco_rt");
  type[GROUP_MATERI_PLASTI_DIPRISCO_RT] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_DIPRISCO_RT] = 1;
  data_class[GROUP_MATERI_PLASTI_DIPRISCO_RT] = MATERI;
  data_required[GROUP_MATERI_PLASTI_DIPRISCO_RT] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_DRUCKPRAG],"group_materi_plasti_druckprag");
  type[GROUP_MATERI_PLASTI_DRUCKPRAG] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_DRUCKPRAG] = 3;
  data_class[GROUP_MATERI_PLASTI_DRUCKPRAG] = MATERI;
  data_required[GROUP_MATERI_PLASTI_DRUCKPRAG] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_DRUCKPRAG_TENSIONCUTOFF],"group_materi_plasti_druckprag_tensioncutoff");
  type[GROUP_MATERI_PLASTI_DRUCKPRAG_TENSIONCUTOFF] = INTEGER;
  data_length[GROUP_MATERI_PLASTI_DRUCKPRAG_TENSIONCUTOFF] = 1;
  data_class[GROUP_MATERI_PLASTI_DRUCKPRAG_TENSIONCUTOFF] = MATERI;
  data_required[GROUP_MATERI_PLASTI_DRUCKPRAG_TENSIONCUTOFF] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_DRUCKPRAG_TENSIONLIMIT],"group_materi_plasti_druckprag_tensionlimit");
  type[GROUP_MATERI_PLASTI_DRUCKPRAG_TENSIONLIMIT] = INTEGER;
  data_length[GROUP_MATERI_PLASTI_DRUCKPRAG_TENSIONLIMIT] = 1;
  data_class[GROUP_MATERI_PLASTI_DRUCKPRAG_TENSIONLIMIT] = MATERI;
  data_required[GROUP_MATERI_PLASTI_DRUCKPRAG_TENSIONLIMIT] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_GURSON],"group_materi_plasti_gurson");
  type[GROUP_MATERI_PLASTI_GURSON] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_GURSON] = 4;
  data_class[GROUP_MATERI_PLASTI_GURSON] = MATERI;
  data_required[GROUP_MATERI_PLASTI_GURSON] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_HLC],"group_materi_plasti_hlc");
  type[GROUP_MATERI_PLASTI_HLC] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_HLC] = 5;
  data_class[GROUP_MATERI_PLASTI_HLC] = MATERI;
  data_required[GROUP_MATERI_PLASTI_HLC] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_HEATGENERATION], "group_materi_plasti_heatgeneration");
  type[GROUP_MATERI_PLASTI_HEATGENERATION] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_HEATGENERATION] = 1;
  data_class[GROUP_MATERI_PLASTI_HEATGENERATION] = MATERI;
  data_required[GROUP_MATERI_PLASTI_HEATGENERATION] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_HYPO_LOWANGLES],"group_materi_plasti_hypo_lowangles");
  type[GROUP_MATERI_PLASTI_HYPO_LOWANGLES] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_HYPO_LOWANGLES] = 10;
  data_class[GROUP_MATERI_PLASTI_HYPO_LOWANGLES] = MATERI;
  data_required[GROUP_MATERI_PLASTI_HYPO_LOWANGLES] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_HYPO_COHESION],"group_materi_plasti_hypo_cohesion");
  type[GROUP_MATERI_PLASTI_HYPO_COHESION] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_HYPO_COHESION] = 1;
  data_class[GROUP_MATERI_PLASTI_HYPO_COHESION] = MATERI;
  data_required[GROUP_MATERI_PLASTI_HYPO_COHESION] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_HYPO_INTERGRANULARSTRAIN],"group_materi_plasti_hypo_intergranularstrain");
  type[GROUP_MATERI_PLASTI_HYPO_INTERGRANULARSTRAIN] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_HYPO_INTERGRANULARSTRAIN] = 5;
  data_class[GROUP_MATERI_PLASTI_HYPO_INTERGRANULARSTRAIN] = MATERI;
  data_required[GROUP_MATERI_PLASTI_HYPO_INTERGRANULARSTRAIN] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_HYPO_PRESSUREDEPENDENTVOIDRATIO],
    "group_materi_plasti_hypo_pressuredependentvoidratio");
  type[GROUP_MATERI_PLASTI_HYPO_PRESSUREDEPENDENTVOIDRATIO] = INTEGER;
  data_length[GROUP_MATERI_PLASTI_HYPO_PRESSUREDEPENDENTVOIDRATIO] = 1;
  data_class[GROUP_MATERI_PLASTI_HYPO_PRESSUREDEPENDENTVOIDRATIO] = MATERI;
  data_required[GROUP_MATERI_PLASTI_HYPO_PRESSUREDEPENDENTVOIDRATIO] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_HYPO_WOLFERSDORFF],"group_materi_plasti_hypo_wolfersdorff");
  type[GROUP_MATERI_PLASTI_HYPO_WOLFERSDORFF] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_HYPO_WOLFERSDORFF] = 8;
  data_class[GROUP_MATERI_PLASTI_HYPO_WOLFERSDORFF] = MATERI;
  data_required[GROUP_MATERI_PLASTI_HYPO_WOLFERSDORFF] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_INCREMENTAL_ERASERECENTHISTORY], "group_materi_plasti_incremental_eraserecenthistory");
  type[GROUP_MATERI_PLASTI_INCREMENTAL_ERASERECENTHISTORY] = INTEGER;
  data_length[GROUP_MATERI_PLASTI_INCREMENTAL_ERASERECENTHISTORY] = 1;
  data_class[GROUP_MATERI_PLASTI_INCREMENTAL_ERASERECENTHISTORY] = MATERI;
  data_required[GROUP_MATERI_PLASTI_INCREMENTAL_ERASERECENTHISTORY] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_INCREMENTAL_FEERROR],"group_materi_plasti_incremental_feerror");
  type[GROUP_MATERI_PLASTI_INCREMENTAL_FEERROR] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_INCREMENTAL_FEERROR] = 1;
  data_class[GROUP_MATERI_PLASTI_INCREMENTAL_FEERROR] = MATERI;
  data_required[GROUP_MATERI_PLASTI_INCREMENTAL_FEERROR] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_INCREMENTAL_FESUBSTEPS],"group_materi_plasti_incremental_fesubsteps");
  type[GROUP_MATERI_PLASTI_INCREMENTAL_FESUBSTEPS] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_INCREMENTAL_FESUBSTEPS] = 1;
  data_class[GROUP_MATERI_PLASTI_INCREMENTAL_FESUBSTEPS] = MATERI;
  data_required[GROUP_MATERI_PLASTI_INCREMENTAL_FESUBSTEPS] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_INCREMENTAL_MAXSUBSTEPS],"group_materi_plasti_incremental_maxsubsteps");
  type[GROUP_MATERI_PLASTI_INCREMENTAL_MAXSUBSTEPS] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_INCREMENTAL_MAXSUBSTEPS] = 1;
  data_class[GROUP_MATERI_PLASTI_INCREMENTAL_MAXSUBSTEPS] = MATERI;
  data_required[GROUP_MATERI_PLASTI_INCREMENTAL_MAXSUBSTEPS] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_INCREMENTAL_MINSUBSTEPS],"group_materi_plasti_incremental_minsubsteps");
  type[GROUP_MATERI_PLASTI_INCREMENTAL_MINSUBSTEPS] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_INCREMENTAL_MINSUBSTEPS] = 1;
  data_class[GROUP_MATERI_PLASTI_INCREMENTAL_MINSUBSTEPS] = MATERI;
  data_required[GROUP_MATERI_PLASTI_INCREMENTAL_MINSUBSTEPS] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_INCREMENTAL_PRINTF],"group_materi_plasti_incremental_printf");
  type[GROUP_MATERI_PLASTI_INCREMENTAL_PRINTF] = INTEGER;
  data_length[GROUP_MATERI_PLASTI_INCREMENTAL_PRINTF] = 1;
  data_class[GROUP_MATERI_PLASTI_INCREMENTAL_PRINTF] = MATERI;
  data_required[GROUP_MATERI_PLASTI_INCREMENTAL_PRINTF] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_INCREMENTAL_USEMATRIX],"group_materi_plasti_incremental_usematrix");
  type[GROUP_MATERI_PLASTI_INCREMENTAL_USEMATRIX] = INTEGER;
  data_length[GROUP_MATERI_PLASTI_INCREMENTAL_USEMATRIX] = 1;
  data_class[GROUP_MATERI_PLASTI_INCREMENTAL_USEMATRIX] = MATERI;
  data_required[GROUP_MATERI_PLASTI_INCREMENTAL_USEMATRIX] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_KINEMATIC_HARDENING],"group_materi_plasti_kinematic_hardening");
  type[GROUP_MATERI_PLASTI_KINEMATIC_HARDENING] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_KINEMATIC_HARDENING] = 1;
  data_class[GROUP_MATERI_PLASTI_KINEMATIC_HARDENING] = MATERI;
  data_required[GROUP_MATERI_PLASTI_KINEMATIC_HARDENING] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_MATSUOKANAKAI],"group_materi_plasti_matsuokanakai");
  type[GROUP_MATERI_PLASTI_MATSUOKANAKAI] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_MATSUOKANAKAI] = 3;
  data_class[GROUP_MATERI_PLASTI_MATSUOKANAKAI] = MATERI;
  data_required[GROUP_MATERI_PLASTI_MATSUOKANAKAI] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_MATSUOKANAKAI_TENSIONCUTOFF],"group_materi_plasti_matsuokanakai_tensioncutoff");
  type[GROUP_MATERI_PLASTI_MATSUOKANAKAI_TENSIONCUTOFF] = INTEGER;
  data_length[GROUP_MATERI_PLASTI_MATSUOKANAKAI_TENSIONCUTOFF] = 1;
  data_class[GROUP_MATERI_PLASTI_MATSUOKANAKAI_TENSIONCUTOFF] = MATERI;
  data_required[GROUP_MATERI_PLASTI_MATSUOKANAKAI_TENSIONCUTOFF] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_MOHRCOUL],"group_materi_plasti_mohrcoul");
  type[GROUP_MATERI_PLASTI_MOHRCOUL] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_MOHRCOUL] = 3;
  data_class[GROUP_MATERI_PLASTI_MOHRCOUL] = MATERI;
  data_required[GROUP_MATERI_PLASTI_MOHRCOUL] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_MOHRCOUL_SOFTENING],"group_materi_plasti_mohrcoul_softening");
  type[GROUP_MATERI_PLASTI_MOHRCOUL_SOFTENING] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_MOHRCOUL_SOFTENING] = 7;
  data_class[GROUP_MATERI_PLASTI_MOHRCOUL_SOFTENING] = MATERI;
  data_required[GROUP_MATERI_PLASTI_MOHRCOUL_SOFTENING] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_MOHRCOUL_TENSIONCUTOFF],"group_materi_plasti_mohrcoul_tensioncutoff");
  type[GROUP_MATERI_PLASTI_MOHRCOUL_TENSIONCUTOFF] = INTEGER;
  data_length[GROUP_MATERI_PLASTI_MOHRCOUL_TENSIONCUTOFF] = 1;
  data_class[GROUP_MATERI_PLASTI_MOHRCOUL_TENSIONCUTOFF] = MATERI;
  data_required[GROUP_MATERI_PLASTI_MOHRCOUL_TENSIONCUTOFF] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_STRESS],"group_materi_plasti_stress");
  type[GROUP_MATERI_PLASTI_STRESS] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_STRESS] = 1;
  data_class[GROUP_MATERI_PLASTI_STRESS] = MATERI;
  data_required[GROUP_MATERI_PLASTI_STRESS] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_TENSION],"group_materi_plasti_tension");
  type[GROUP_MATERI_PLASTI_TENSION] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_TENSION] = 1;
  data_class[GROUP_MATERI_PLASTI_TENSION] = MATERI;
  data_required[GROUP_MATERI_PLASTI_TENSION] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_TSKH],"group_materi_plasti_tskh");
  type[GROUP_MATERI_PLASTI_TSKH] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_TSKH] = DATA_ITEM_SIZE;
  data_class[GROUP_MATERI_PLASTI_TSKH] = MATERI;
  fixed_length[GROUP_MATERI_PLASTI_TSKH] = 0;
  data_required[GROUP_MATERI_PLASTI_TSKH] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_USER],"group_materi_plasti_user");
  type[GROUP_MATERI_PLASTI_USER] = INTEGER;
  data_length[GROUP_MATERI_PLASTI_USER] = 1;
  data_class[GROUP_MATERI_PLASTI_USER] = MATERI;
  data_required[GROUP_MATERI_PLASTI_USER] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_VISCO_ALWAYS],"group_materi_plasti_visco_always");
  type[GROUP_MATERI_PLASTI_VISCO_ALWAYS] = INTEGER;
  data_length[GROUP_MATERI_PLASTI_VISCO_ALWAYS] = 1;
  data_class[GROUP_MATERI_PLASTI_VISCO_ALWAYS] = MATERI;
  data_required[GROUP_MATERI_PLASTI_VISCO_ALWAYS] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL],"group_materi_plasti_visco_exponential");
  type[GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL] = 2;
  data_class[GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL] = MATERI;
  data_required[GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_VISCO_POWER],"group_materi_plasti_visco_power");
  type[GROUP_MATERI_PLASTI_VISCO_POWER] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_VISCO_POWER] = 3;
  data_class[GROUP_MATERI_PLASTI_VISCO_POWER] = MATERI;
  data_required[GROUP_MATERI_PLASTI_VISCO_POWER] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_VONMISES],"group_materi_plasti_vonmises");
  type[GROUP_MATERI_PLASTI_VONMISES] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_VONMISES] = 1;
  data_class[GROUP_MATERI_PLASTI_VONMISES] = MATERI;
  data_required[GROUP_MATERI_PLASTI_VONMISES] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_VONMISES_NADAI],"group_materi_plasti_vonmises_nadai");
  type[GROUP_MATERI_PLASTI_VONMISES_NADAI] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_VONMISES_NADAI] = 3;
  data_class[GROUP_MATERI_PLASTI_VONMISES_NADAI] = MATERI;
  data_required[GROUP_MATERI_PLASTI_VONMISES_NADAI] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_STOKES],"group_materi_stokes");
  type[GROUP_MATERI_STOKES] = INTEGER;
  data_length[GROUP_MATERI_STOKES] = 1;
  data_class[GROUP_MATERI_STOKES] = MATERI;          
  data_required[GROUP_MATERI_STOKES] = GROUP_TYPE;          

  strcpy(name[GROUP_MATERI_VISCOSITY],"group_materi_viscosity");
  type[GROUP_MATERI_VISCOSITY] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_VISCOSITY] = 1;
  data_class[GROUP_MATERI_VISCOSITY] = MATERI;
  data_required[GROUP_MATERI_VISCOSITY] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_VISCOSITY_HEATGENERATION],"group_materi_viscosity_heatgeneration");
  type[GROUP_MATERI_VISCOSITY_HEATGENERATION] = INTEGER;
  data_length[GROUP_MATERI_VISCOSITY_HEATGENERATION] = 1;
  data_class[GROUP_MATERI_VISCOSITY_HEATGENERATION] = MATERI;
  data_required[GROUP_MATERI_VISCOSITY_HEATGENERATION] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_VISCOSITY_USER],"group_materi_viscosity_user");
  type[GROUP_MATERI_VISCOSITY_USER] = INTEGER;
  data_length[GROUP_MATERI_VISCOSITY_USER] = 1;
  data_class[GROUP_MATERI_VISCOSITY_USER] = MATERI;
  data_required[GROUP_MATERI_VISCOSITY_USER] = GROUP_TYPE;

  strcpy(name[GROUP_MATRIX_SECOND_VALUES],"group_matrix_second_values");
  type[GROUP_MATRIX_SECOND_VALUES] = DOUBLE_PRECISION;
  data_length[GROUP_MATRIX_SECOND_VALUES] = 0; // set run-time
  fixed_length[GROUP_MATRIX_SECOND_VALUES] = 0;
  data_class[GROUP_MATRIX_SECOND_VALUES] = GROUP_MATRIX_SECOND_VALUES;
  external[GROUP_MATRIX_SECOND_VALUES] = 0;
  data_required[GROUP_MATRIX_SECOND_VALUES] = GROUP_TYPE;

  strcpy(name[GROUP_MATRIX_UNKNOWNS],"group_matrix_unknowns");
  type[GROUP_MATRIX_UNKNOWNS] = INTEGER;
  data_length[GROUP_MATRIX_UNKNOWNS] = 0; // set run-time
  fixed_length[GROUP_MATRIX_UNKNOWNS] = 0;
  data_class[GROUP_MATRIX_UNKNOWNS] = GROUP_MATRIX_UNKNOWNS;
  external[GROUP_MATRIX_UNKNOWNS] = 0;
  data_required[GROUP_MATRIX_UNKNOWNS] = GROUP_TYPE;

  strcpy(name[GROUP_MATRIX_VALUES],"group_matrix_values");
  type[GROUP_MATRIX_VALUES] = DOUBLE_PRECISION;
  data_length[GROUP_MATRIX_VALUES] = 0; // set run-time
  fixed_length[GROUP_MATRIX_VALUES] = 0;
  data_class[GROUP_MATRIX_VALUES] = GROUP_MATRIX_VALUES;
  external[GROUP_MATRIX_VALUES] = 0;
  data_required[GROUP_MATRIX_VALUES] = GROUP_TYPE;

  strcpy(name[GROUP_MAXWELL_FREQUENCY_EIGEN],"group_maxwell_frequency_eigen");
  type[GROUP_MAXWELL_FREQUENCY_EIGEN] = INTEGER;
  data_length[GROUP_MAXWELL_FREQUENCY_EIGEN] = 1;
  data_class[GROUP_MAXWELL_FREQUENCY_EIGEN] = MAXWELL_FREQUENCY;
  data_required[GROUP_MAXWELL_FREQUENCY_EIGEN] = GROUP_TYPE;

  strcpy(name[GROUP_MAXWELL_FREQUENCY_EPSILON],"group_maxwell_frequency_epsilon");
  type[GROUP_MAXWELL_FREQUENCY_EPSILON] = DOUBLE_PRECISION;
  data_length[GROUP_MAXWELL_FREQUENCY_EPSILON] = 2;
  data_class[GROUP_MAXWELL_FREQUENCY_EPSILON] = MAXWELL_FREQUENCY;
  data_required[GROUP_MAXWELL_FREQUENCY_EPSILON] = GROUP_TYPE;

  strcpy(name[GROUP_MAXWELL_FREQUENCY_EPSILON_ANISOTROPIC],"group_maxwell_frequency_epsilon_anisotropic");
  type[GROUP_MAXWELL_FREQUENCY_EPSILON_ANISOTROPIC] = DOUBLE_PRECISION;
  data_length[GROUP_MAXWELL_FREQUENCY_EPSILON_ANISOTROPIC] = MDIM*2;
  data_class[GROUP_MAXWELL_FREQUENCY_EPSILON_ANISOTROPIC] = MAXWELL_FREQUENCY;
  data_required[GROUP_MAXWELL_FREQUENCY_EPSILON_ANISOTROPIC] = GROUP_TYPE;
 
  strcpy(name[GROUP_MAXWELL_FREQUENCY_J],"group_maxwell_frequency_j");
  type[GROUP_MAXWELL_FREQUENCY_J] = DOUBLE_PRECISION;
  data_length[GROUP_MAXWELL_FREQUENCY_J] = 2*MDIM;
  data_class[GROUP_MAXWELL_FREQUENCY_J] = MAXWELL_FREQUENCY;
  data_required[GROUP_MAXWELL_FREQUENCY_J] = GROUP_TYPE;

  strcpy(name[GROUP_MAXWELL_FREQUENCY_MU],"group_maxwell_frequency_mu");
  type[GROUP_MAXWELL_FREQUENCY_MU] = DOUBLE_PRECISION;
  data_length[GROUP_MAXWELL_FREQUENCY_MU] = 2;
  data_class[GROUP_MAXWELL_FREQUENCY_MU] = MAXWELL_FREQUENCY;
  data_required[GROUP_MAXWELL_FREQUENCY_MU] = GROUP_TYPE;

  strcpy(name[GROUP_MAXWELL_FREQUENCY_PENALTY],"group_maxwell_frequency_penalty");
  type[GROUP_MAXWELL_FREQUENCY_PENALTY] = DOUBLE_PRECISION;
  data_length[GROUP_MAXWELL_FREQUENCY_PENALTY] = 1;
  data_class[GROUP_MAXWELL_FREQUENCY_PENALTY] = MAXWELL_FREQUENCY;
  data_required[GROUP_MAXWELL_FREQUENCY_PENALTY] = GROUP_TYPE;

  strcpy(name[GROUP_MAXWELL_FREQUENCY_PML_EPSILONANDMU],
    "group_maxwell_frequency_pml_epsilonandmu");
  type[GROUP_MAXWELL_FREQUENCY_PML_EPSILONANDMU] = DOUBLE_PRECISION;
  data_length[GROUP_MAXWELL_FREQUENCY_PML_EPSILONANDMU] = 2;
  data_class[GROUP_MAXWELL_FREQUENCY_PML_EPSILONANDMU] = MAXWELL_FREQUENCY;
  data_required[GROUP_MAXWELL_FREQUENCY_PML_EPSILONANDMU] = GROUP_TYPE;
 
  strcpy(name[GROUP_MAXWELL_FREQUENCY_PML_PLANES],
    "group_maxwell_frequency_pml_planes");
  type[GROUP_MAXWELL_FREQUENCY_PML_PLANES] = INTEGER;
  data_length[GROUP_MAXWELL_FREQUENCY_PML_PLANES] = 2;
  data_class[GROUP_MAXWELL_FREQUENCY_PML_PLANES] = MAXWELL_FREQUENCY;
  data_required[GROUP_MAXWELL_FREQUENCY_PML_PLANES] = GROUP_TYPE;

  strcpy(name[GROUP_MAXWELL_TIME_EPSILON],"group_maxwell_time_epsilon");
  type[GROUP_MAXWELL_TIME_EPSILON] = DOUBLE_PRECISION;
  data_length[GROUP_MAXWELL_TIME_EPSILON] = 1;
  data_class[GROUP_MAXWELL_TIME_EPSILON] = MAXWELL_TIME;
  data_required[GROUP_MAXWELL_TIME_EPSILON] = GROUP_TYPE;

  strcpy(name[GROUP_MAXWELL_TIME_EPSILON_ANISOTROPIC],"group_maxwell_time_epsilon_anisotropic");
  type[GROUP_MAXWELL_TIME_EPSILON_ANISOTROPIC] = DOUBLE_PRECISION;
  data_length[GROUP_MAXWELL_TIME_EPSILON_ANISOTROPIC] = MDIM;
  data_class[GROUP_MAXWELL_TIME_EPSILON_ANISOTROPIC] = MAXWELL_TIME;
  data_required[GROUP_MAXWELL_TIME_EPSILON_ANISOTROPIC] = GROUP_TYPE;
 
  strcpy(name[GROUP_MAXWELL_TIME_J],"group_maxwell_time_j");
  type[GROUP_MAXWELL_TIME_J] = DOUBLE_PRECISION;
  data_length[GROUP_MAXWELL_TIME_J] = MDIM;
  data_class[GROUP_MAXWELL_TIME_J] = MAXWELL_TIME; 
  data_required[GROUP_MAXWELL_TIME_J] = GROUP_TYPE; 

  strcpy(name[GROUP_MAXWELL_TIME_MU],"group_maxwell_time_mu");
  type[GROUP_MAXWELL_TIME_MU] = DOUBLE_PRECISION;
  data_length[GROUP_MAXWELL_TIME_MU] = 1;
  data_class[GROUP_MAXWELL_TIME_MU] = MAXWELL_TIME;
  data_required[GROUP_MAXWELL_TIME_MU] = GROUP_TYPE;

  strcpy(name[GROUP_MAXWELL_TIME_PENALTY],"group_maxwell_time_penalty");
  type[GROUP_MAXWELL_TIME_PENALTY] = DOUBLE_PRECISION;
  data_length[GROUP_MAXWELL_TIME_PENALTY] = 1;
  data_class[GROUP_MAXWELL_TIME_PENALTY] = MAXWELL_TIME;
  data_required[GROUP_MAXWELL_TIME_PENALTY] = GROUP_TYPE;

  strcpy(name[GROUP_SPRING_DIRECTION],"group_spring_direction");
  type[GROUP_SPRING_DIRECTION] = DOUBLE_PRECISION;
  data_length[GROUP_SPRING_DIRECTION] = ndim;
  data_class[GROUP_SPRING_DIRECTION] = SPRING;
  data_required[GROUP_SPRING_DIRECTION] = GROUP_TYPE;

  strcpy(name[GROUP_SPRING_STIFFNESS],"group_spring_stiffness");
  type[GROUP_SPRING_STIFFNESS] = DOUBLE_PRECISION;
  data_length[GROUP_SPRING_STIFFNESS] = 1;
  data_class[GROUP_SPRING_STIFFNESS] = SPRING;
  data_required[GROUP_SPRING_STIFFNESS] = GROUP_TYPE;

  strcpy(name[GROUP_SPRING_PLASTI],"group_spring_plasti");
  type[GROUP_SPRING_PLASTI] = DOUBLE_PRECISION;
  data_length[GROUP_SPRING_PLASTI] = 1;
  data_class[GROUP_SPRING_PLASTI] = SPRING;
  data_required[GROUP_SPRING_PLASTI] = GROUP_TYPE;

  strcpy(name[GROUP_TIME],"group_time");
  type[GROUP_TIME] = DOUBLE_PRECISION;
  data_length[GROUP_TIME] = 2;
  data_class[GROUP_TIME] = GROUP_TIME;
  data_required[GROUP_TIME] = GROUP_TYPE;

  strcpy(name[GROUP_TRUSS_AREA],"group_truss_area");
  type[GROUP_TRUSS_AREA] = DOUBLE_PRECISION;
  data_length[GROUP_TRUSS_AREA] = 1;
  data_class[GROUP_TRUSS_AREA] = TRUSS;
  data_required[GROUP_TRUSS_AREA] = GROUP_TYPE;

  strcpy(name[GROUP_TRUSS_DENSITY],"group_truss_density");
  type[GROUP_TRUSS_DENSITY] = DOUBLE_PRECISION;
  data_length[GROUP_TRUSS_DENSITY] = 1;
  data_class[GROUP_TRUSS_DENSITY] = TRUSS;
  data_required[GROUP_TRUSS_DENSITY] = GROUP_TYPE;

  strcpy(name[GROUP_TRUSS_MEMORY],"group_truss_memory");
  type[GROUP_TRUSS_MEMORY] = INTEGER;
  data_length[GROUP_TRUSS_MEMORY] = 1;
  data_class[GROUP_TRUSS_MEMORY] = TRUSS;
  data_required[GROUP_TRUSS_MEMORY] = GROUP_TYPE;

  strcpy(name[GROUP_TRUSS_ROPE],"group_truss_rope");
  type[GROUP_TRUSS_ROPE] = INTEGER;
  data_length[GROUP_TRUSS_ROPE] = 1;
  data_class[GROUP_TRUSS_ROPE] = TRUSS;
  data_required[GROUP_TRUSS_ROPE] = GROUP_TYPE;

  strcpy(name[GROUP_TRUSS_PLASTI],"group_truss_plasti");
  type[GROUP_TRUSS_PLASTI] = DOUBLE_PRECISION;
  data_length[GROUP_TRUSS_PLASTI] = 1;
  data_class[GROUP_TRUSS_PLASTI] = TRUSS;
  data_required[GROUP_TRUSS_PLASTI] = GROUP_TYPE;

  strcpy(name[GROUP_TRUSS_YOUNG],"group_truss_young");
  type[GROUP_TRUSS_YOUNG] = DOUBLE_PRECISION;
  data_length[GROUP_TRUSS_YOUNG] = 1;
  data_class[GROUP_TRUSS_YOUNG] = TRUSS;
  data_required[GROUP_TRUSS_YOUNG] = GROUP_TYPE;

  strcpy(name[GROUP_TYPE],"group_type");
  type[GROUP_TYPE] = INTEGER;
  data_length[GROUP_TYPE] = MTYPE;
  fixed_length[GROUP_TYPE] = 0;
  data_required[GROUP_TYPE] = GROUP_TYPE;

  strcpy(name[GROUP_USER_DATA],"group_user_data");
  type[GROUP_USER_DATA] = DOUBLE_PRECISION;
  data_length[GROUP_USER_DATA] = DATA_ITEM_SIZE;
  fixed_length[GROUP_USER_DATA] = 0;
  data_required[GROUP_USER_DATA] = GROUP_TYPE;

  strcpy(name[GROUP_USER_UMAT],"group_user_umat");
  type[GROUP_USER_UMAT] = INTEGER;
  data_length[GROUP_USER_UMAT] = 1;
  data_required[GROUP_USER_UMAT] = GROUP_TYPE;

  strcpy(name[GROUP_VOLUME_FACTOR],"group_volume_factor");
  type[GROUP_VOLUME_FACTOR] = DOUBLE_PRECISION;
  data_length[GROUP_VOLUME_FACTOR] = 1;
  version_all[GROUP_VOLUME_FACTOR] = 1;
  data_class[GROUP_VOLUME_FACTOR] = VOLUME;
  data_required[GROUP_VOLUME_FACTOR] = GROUP_TYPE;

  strcpy(name[GROUP_WAVE_SPEED_OF_SOUND],"group_wave_speed_of_sound");
  type[GROUP_WAVE_SPEED_OF_SOUND] = DOUBLE_PRECISION;
  data_length[GROUP_WAVE_SPEED_OF_SOUND] = 1;
  data_class[GROUP_WAVE_SPEED_OF_SOUND] = WAVE;
  data_required[GROUP_WAVE_SPEED_OF_SOUND] = GROUP_TYPE;

  strcpy(name[GROWTH],"growth");

  strcpy(name[HEX8],"hex8");

  strcpy(name[HEX27],"hex27");

  strcpy(name[HEX64],"hex64");

  strcpy(name[H_REFINEMENT],"h_refinement");

  strcpy(name[ICC],"icc");

  strcpy(name[ILU],"ilu");

  strcpy(name[ICONTROL],"icontrol");
  type[ICONTROL] = INTEGER;
  data_length[ICONTROL] = 1;
  no_index[ICONTROL] = 1;
  external[ICONTROL] = 0;

  strcpy(name[INITIALIZE],"initialize");

  strcpy(name[INITIALIZATION_VALUES],"initialization_values");
  type[INITIALIZATION_VALUES] = INTEGER;
  data_length[INITIALIZATION_VALUES] = DATA_ITEM_SIZE;
  no_index[INITIALIZATION_VALUES] = 1;
  external[INITIALIZATION_VALUES] = 0;
  fixed_length[INITIALIZATION_VALUES] = 0;

  strcpy(name[INVERSE],"inverse");

  strcpy(name[INVERSE_HISTORY],"inverse_history");
  type[INVERSE_HISTORY] = DOUBLE_PRECISION;
  data_length[INVERSE_HISTORY] = 2;
  no_index[INVERSE_HISTORY] = 1;
  external[INVERSE_HISTORY] = 0;
  data_class[INVERSE_HISTORY] = INVERSE;

  strcpy(name[INVERSE_ITERATIONS],"inverse_iterations");
  type[INVERSE_ITERATIONS] = INTEGER;
  data_length[INVERSE_ITERATIONS] = 1;
  no_index[INVERSE_ITERATIONS] = 1;
  data_class[INVERSE_ITERATIONS] = INVERSE;

  strcpy(name[INVERSE_ITERATION_NUMBER],"inverse_iteration_number");
  type[INVERSE_ITERATION_NUMBER] = INTEGER;
  data_length[INVERSE_ITERATION_NUMBER] = 1;
  no_index[INVERSE_ITERATION_NUMBER] = 1;
  data_class[INVERSE_ITERATION_NUMBER] = INVERSE;

  strcpy(name[INVERSE_PARAMETER],"inverse_parameter");
  type[INVERSE_PARAMETER] = INTEGER;
  data_length[INVERSE_PARAMETER] = 3;
  data_class[INVERSE_PARAMETER] = INVERSE;

  strcpy(name[INVERSE_PARAMETER_LIMITS],"inverse_parameter_limits");
  type[INVERSE_PARAMETER_LIMITS] = DOUBLE_PRECISION;
  data_length[INVERSE_PARAMETER_LIMITS] = 2;
  data_class[INVERSE_PARAMETER_LIMITS] = INVERSE;
  data_required[INVERSE_PARAMETER_LIMITS] = INVERSE_PARAMETER;

  strcpy(name[INVERSE_PARAMETER_SENSITIVITY],"inverse_parameter_sensitivity");
  type[INVERSE_PARAMETER_SENSITIVITY] = DOUBLE_PRECISION;
  data_length[INVERSE_PARAMETER_SENSITIVITY] = 3;
  external[INVERSE_PARAMETER_SENSITIVITY] = 0;
  data_class[INVERSE_PARAMETER_SENSITIVITY] = INVERSE;
  data_required[INVERSE_PARAMETER_SENSITIVITY] = INVERSE_PARAMETER;

  strcpy(name[INVERSE_PARAMETER_STEP],"inverse_parameter_step");
  type[INVERSE_PARAMETER_STEP] = DOUBLE_PRECISION;
  data_length[INVERSE_PARAMETER_STEP] = 1;
  data_class[INVERSE_PARAMETER_STEP] = INVERSE;
  data_required[INVERSE_PARAMETER_STEP] = INVERSE_PARAMETER;

  strcpy(name[INVERSE_PARAMETER_VARIATION],"inverse_parameter_variation");
  type[INVERSE_PARAMETER_VARIATION] = DOUBLE_PRECISION;
  data_length[INVERSE_PARAMETER_VARIATION] = 1;
  data_class[INVERSE_PARAMETER_VARIATION] = INVERSE;
  data_required[INVERSE_PARAMETER_VARIATION] = INVERSE_PARAMETER;

  strcpy(name[INVERSE_TARGET],"inverse_target");
  type[INVERSE_TARGET] = INTEGER;
  data_length[INVERSE_TARGET] = 3;
  data_class[INVERSE_TARGET] = INVERSE;

  strcpy(name[INVERSE_TARGET_DATA],"inverse_target_data");
  type[INVERSE_TARGET_DATA] = DOUBLE_PRECISION;
  data_length[INVERSE_TARGET_DATA] = 3;
  fixed_length[INVERSE_TARGET_DATA] = 0;
  data_class[INVERSE_TARGET_DATA] = INVERSE;
  data_required[INVERSE_TARGET_DATA] = INVERSE_TARGET;

  strcpy(name[INVERSE_TARGET_TIMESTEP],"inverse_target_timestep");
  type[INVERSE_TARGET_TIMESTEP] = INTEGER;
  data_length[INVERSE_TARGET_TIMESTEP] = 1;
  data_class[INVERSE_TARGET_TIMESTEP] = INVERSE;
  data_required[INVERSE_TARGET_TIMESTEP] = INVERSE_TARGET;

  strcpy(name[JACOBI],"jacobi");

  strcpy(name[LSQR],"lsqr");

  strcpy(name[LOBATTO],"lobatto");

  strcpy(name[LU],"lu");

  strcpy(name[MACRO],"macro");

  strcpy(name[MATERI],"materi");

  strcpy(name[MATERI_DAMAGE],"materi_damage");

  strcpy(name[MATERI_DENSITY],"materi_density");

  strcpy(name[MATERI_DENSITY_MINIMUM],"materi_density_minimum");
  type[MATERI_DENSITY_MINIMUM] = DOUBLE_PRECISION;
  data_length[MATERI_DENSITY_MINIMUM] = 1;
  no_index[MATERI_DENSITY_MINIMUM] = 1;

  strcpy(name[MATERI_DIFFUSION],"materi_diffusion");

  strcpy(name[MATERI_DIFFUSION_MINIMUM],"materi_diffusion_minimum");
  type[MATERI_DIFFUSION_MINIMUM] = DOUBLE_PRECISION;
  data_length[MATERI_DIFFUSION_MINIMUM] = 1;
  no_index[MATERI_DIFFUSION_MINIMUM] = 1;

  strcpy(name[MATERI_DIFFUSION_ADJUST_GEOMETRY],"materi_diffusion_adjust_geometry");
  type[MATERI_DIFFUSION_ADJUST_GEOMETRY] = INTEGER;
  data_length[MATERI_DIFFUSION_ADJUST_GEOMETRY] = 4;
  data_class[MATERI_DIFFUSION_ADJUST_GEOMETRY] = MATERI_DIFFUSION;

  strcpy(name[MATERI_DIFFUSION_CORRECT],"materi_diffusion_correct");
  type[MATERI_DIFFUSION_CORRECT] = INTEGER;
  data_length[MATERI_DIFFUSION_CORRECT] = 1;
  data_class[MATERI_DIFFUSION_CORRECT] = MATERI_DIFFUSION;
  no_index[MATERI_DIFFUSION_CORRECT] = 1;

  strcpy(name[MATERI_DIFFUSION_FILL_GEOMETRY],"materi_diffusion_fill_geometry");
  type[MATERI_DIFFUSION_FILL_GEOMETRY] = INTEGER;
  data_length[MATERI_DIFFUSION_FILL_GEOMETRY] = 2;
  data_class[MATERI_DIFFUSION_FILL_GEOMETRY] = MATERI_DIFFUSION;

  strcpy(name[MATERI_DIFFUSION_FILL_EPSVELOCITY],"materi_diffusion_fill_epsvelocity");
  type[MATERI_DIFFUSION_FILL_EPSVELOCITY] = DOUBLE_PRECISION;
  data_length[MATERI_DIFFUSION_FILL_EPSVELOCITY] = 1;
  data_class[MATERI_DIFFUSION_FILL_EPSVELOCITY] = MATERI_DIFFUSION;
  no_index[MATERI_DIFFUSION_FILL_EPSVELOCITY] = 1;

  strcpy(name[MATERI_DIFFUSION_SMOOTH],"materi_diffusion_smooth");
  type[MATERI_DIFFUSION_SMOOTH] = INTEGER;
  data_length[MATERI_DIFFUSION_SMOOTH] = 1;
  data_class[MATERI_DIFFUSION_SMOOTH] = MATERI_DIFFUSION;
  no_index[MATERI_DIFFUSION_SMOOTH] = 1;

  strcpy(name[MATERI_DIFFUSION_TEMPERATURE],"materi_diffusion_temperature");
  type[MATERI_DIFFUSION_TEMPERATURE] = DOUBLE_PRECISION;
  data_length[MATERI_DIFFUSION_TEMPERATURE] = 1;
  data_class[MATERI_DIFFUSION_TEMPERATURE] = MATERI_DIFFUSION;
  no_index[MATERI_DIFFUSION_TEMPERATURE] = 1;

  strcpy(name[MATERI_DISPLACEMENT],"materi_displacement");

  strcpy(name[MATERI_HISTORY_VARIABLES],"materi_history_variables");

  strcpy(name[MATERI_MAXWELL_STRESS],"materi_maxwell_stress");

  strcpy(name[MATERI_PLASTI_F],"materi_plasti_f");

  strcpy(name[MATERI_PLASTI_F_NONLOCAL],"materi_plasti_f_nonlocal");

  strcpy(name[MATERI_PLASTI_INCREMENTAL_SUBSTEPS],"materi_plasti_incremental_substeps");

  strcpy(name[MATERI_PLASTI_KAPPA],"materi_plasti_kappa");

  strcpy(name[MATERI_PLASTI_RHO],"materi_plasti_rho");

  strcpy(name[MATERI_PLASTI_SOFTVAR_LOCAL],"materi_plasti_softvar_local");

  strcpy(name[MATERI_PLASTI_SOFTVAR_NONLOCAL],"materi_plasti_softvar_nonlocal");

  strcpy(name[MATERI_ROTATION],"materi_rotation");

  strcpy(name[MATERI_STRAINENERGY],"materi_strainenergy");

  strcpy(name[MATERI_STRAIN_ELASTI],"materi_strain_elasti");

  strcpy(name[MATERI_STRAIN_INTERGRANULAR],"materi_strain_intergranular");

  strcpy(name[MATERI_STRAIN_PLASTI],"materi_strain_plasti");

  strcpy(name[MATERI_STRAIN_TOTAL],"materi_strain_total" );

  strcpy(name[MATERI_STRESS],"materi_stress");

  strcpy(name[MATERI_VELOCITY],"materi_velocity");

  strcpy(name[MATERI_VELOCITY_INTEGRATED],"materi_velocity_integrated");

  strcpy(name[MATERI_VOID_FRACTION],"materi_void_fraction");

  strcpy(name[MATERI_WORK],"materi_work");

  strcpy(name[MATRIX],"matrix");

  strcpy(name[MATRIX_ITERATIVE_BICG],"matrix_iterative_bicg");

  strcpy(name[MATRIX_ITERATIVE_PETSC],"matrix_iterative_petsc");

  strcpy(name[MATRIX_SUPERLU],"matrix_superlu");

  strcpy(name[MATRIX_SUPERLU_DIST],"matrix_superlu_dist");

  strcpy(name[MATRIX_SUPERLU_MT],"matrix_superlu_mt");
  
  strcpy(name[MATRIX_LAPACK],"matrix_lapack");

  strcpy(name[MAXFRE],"maxfre");

  strcpy(name[MAXIMAL],"maximal");

  strcpy(name[MAXTIM],"maxtim");

  strcpy(name[MAXWELL],"maxwell");

  strcpy(name[MAXWELL_ECOMPLEX],"maxwell_ecomplex");

  strcpy(name[MAXWELL_E],"maxwell_e");

  strcpy(name[MAXWELL_EI],"maxwell_ei");

  strcpy(name[MAXWELL_ER],"maxwell_er");

  strcpy(name[MAXWELL_FE],"maxwell_fe");

  strcpy(name[MAXWELL_FREQUENCY],"maxwell_frequency");

  strcpy(name[MAXWELL_FREQUENCY_EXCITATION],"maxwell_frequency_excitation");
  type[MAXWELL_FREQUENCY_EXCITATION] = DOUBLE_PRECISION;
  data_length[MAXWELL_FREQUENCY_EXCITATION] = 1;
  data_class[MAXWELL_FREQUENCY_EXCITATION] = MAXWELL_FREQUENCY;
  no_index[MAXWELL_FREQUENCY_EXCITATION] = 1;

  strcpy(name[MAXWELL_SCATTER_ENERGYCONSERVATION],"maxwell_scatter_energyconservation");
  type[MAXWELL_SCATTER_ENERGYCONSERVATION] = INTEGER;
  data_length[MAXWELL_SCATTER_ENERGYCONSERVATION] = 1;
  data_class[MAXWELL_SCATTER_ENERGYCONSERVATION] = MAXWELL;
  no_index[MAXWELL_SCATTER_ENERGYCONSERVATION] = 1;

  strcpy(name[MAXWELL_SCATTER_MATRIX_AMPLITUDE],"maxwell_scatter_matrix_amplitude");
  type[MAXWELL_SCATTER_MATRIX_AMPLITUDE] = DOUBLE_PRECISION;
  data_length[MAXWELL_SCATTER_MATRIX_AMPLITUDE] = DATA_ITEM_SIZE;
  data_class[MAXWELL_SCATTER_MATRIX_AMPLITUDE] = MAXWELL;
  fixed_length[MAXWELL_SCATTER_MATRIX_AMPLITUDE] = 0;
  no_index[MAXWELL_SCATTER_MATRIX_AMPLITUDE] = 1;

  strcpy(name[MAXWELL_SCATTER_MATRIX_AMPLITUDEDB],"maxwell_scatter_matrix_amplitudedb");
  type[MAXWELL_SCATTER_MATRIX_AMPLITUDEDB] = DOUBLE_PRECISION;
  data_length[MAXWELL_SCATTER_MATRIX_AMPLITUDEDB] = DATA_ITEM_SIZE;
  data_class[MAXWELL_SCATTER_MATRIX_AMPLITUDEDB] = MAXWELL;
  fixed_length[MAXWELL_SCATTER_MATRIX_AMPLITUDEDB] = 0;
  no_index[MAXWELL_SCATTER_MATRIX_AMPLITUDEDB] = 1;

  strcpy(name[MAXWELL_SCATTER_MATRIX_IMAGINARY],"maxwell_scatter_matrix_imaginary");
  type[MAXWELL_SCATTER_MATRIX_IMAGINARY] = DOUBLE_PRECISION;
  data_length[MAXWELL_SCATTER_MATRIX_IMAGINARY] = DATA_ITEM_SIZE;
  data_class[MAXWELL_SCATTER_MATRIX_IMAGINARY] = MAXWELL;
  fixed_length[MAXWELL_SCATTER_MATRIX_IMAGINARY] = 0;
  no_index[MAXWELL_SCATTER_MATRIX_IMAGINARY] = 1;

  strcpy(name[MAXWELL_SCATTER_MATRIX_PHASE],"maxwell_scatter_matrix_phase");
  type[MAXWELL_SCATTER_MATRIX_PHASE] = DOUBLE_PRECISION;
  data_length[MAXWELL_SCATTER_MATRIX_PHASE] = DATA_ITEM_SIZE;
  data_class[MAXWELL_SCATTER_MATRIX_PHASE] = MAXWELL;
  fixed_length[MAXWELL_SCATTER_MATRIX_PHASE] = 0;
  no_index[MAXWELL_SCATTER_MATRIX_PHASE] = 1;

  strcpy(name[MAXWELL_SCATTER_MATRIX_REAL],"maxwell_scatter_matrix_real");
  type[MAXWELL_SCATTER_MATRIX_REAL] = DOUBLE_PRECISION;
  data_length[MAXWELL_SCATTER_MATRIX_REAL] = DATA_ITEM_SIZE;
  data_class[MAXWELL_SCATTER_MATRIX_REAL] = MAXWELL;
  fixed_length[MAXWELL_SCATTER_MATRIX_REAL] = 0;
  no_index[MAXWELL_SCATTER_MATRIX_REAL] = 1;

  strcpy(name[MAXWELL_SCATTER_PORT_INPUT],"maxwell_scatter_port_input");
  type[MAXWELL_SCATTER_PORT_INPUT] = INTEGER;
  data_length[MAXWELL_SCATTER_PORT_INPUT] = 2;
  data_class[MAXWELL_SCATTER_PORT_INPUT] = MAXWELL;
  no_index[MAXWELL_SCATTER_PORT_INPUT] = 1;

  strcpy(name[MAXWELL_SCATTER_PORT_OUTPUT],"maxwell_scatter_port_output");
  type[MAXWELL_SCATTER_PORT_OUTPUT] = INTEGER;
  data_length[MAXWELL_SCATTER_PORT_OUTPUT] = DATA_ITEM_SIZE;
  data_class[MAXWELL_SCATTER_PORT_OUTPUT] = MAXWELL;
  no_index[MAXWELL_SCATTER_PORT_OUTPUT] = 1;
  fixed_length[MAXWELL_SCATTER_PORT_OUTPUT] = 0;

  strcpy(name[MAXWELL_SCATTER_PARAMETERS],"maxwell_scatter_parameters");
  type[MAXWELL_SCATTER_PARAMETERS] = DOUBLE_PRECISION;
  data_length[MAXWELL_SCATTER_PARAMETERS] = 2;
  data_class[MAXWELL_SCATTER_PARAMETERS] = MAXWELL;
  no_index[MAXWELL_SCATTER_PARAMETERS] = 1;

  strcpy(name[MAXWELL_TIME],"maxwell_time");

  strcpy(name[MESH],"mesh");

  strcpy(name[METHOD1],"method1");

  strcpy(name[METHOD2],"method2");

  strcpy(name[MINIMAL],"minimal");

  strcpy(name[MINUS_ONE],"minus_one" );

  strcpy(name[MOMENT],"moment");

  strcpy(name[NEGATIVE],"negative");

  strcpy(name[NO],"no");

  strcpy(name[NODE],"node");
  type[NODE] = DOUBLE_PRECISION;
  data_length[NODE] = ndim;
  version_all[NODE] = 1;
  data_class[NODE] = NODE;

  strcpy(name[NODE_ADJUST],"node_adjust");
  type[NODE_ADJUST] = INTEGER;
  data_length[NODE_ADJUST] = 1;
  version_all[NODE_ADJUST] = 1;
  data_class[NODE_ADJUST] = NODE;
  external[NODE_ADJUST] = 0;
  data_required[NODE_ADJUST] = NODE;

  strcpy(name[NODE_BOUNDARY],"node_boundary");
  type[NODE_BOUNDARY] = INTEGER;
  data_length[NODE_BOUNDARY] = 1;
  version_all[NODE_BOUNDARY] = 1;
  data_class[NODE_BOUNDARY] = NODE;
  data_required[NODE_BOUNDARY] = NODE;

  strcpy(name[NODE_BOUNDED],"node_bounded");
  type[NODE_BOUNDED] = INTEGER;
  data_length[NODE_BOUNDED] = npuknwn;
  external[NODE_BOUNDED] = 0;
  data_class[NODE_BOUNDED] = NODE;
  data_required[NODE_BOUNDED] = NODE;

  strcpy(name[NODE_DAMPING],"node_damping");
  type[NODE_DAMPING] = DOUBLE_PRECISION;
  data_length[NODE_DAMPING] = ndim;
  version_all[NODE_DAMPING] = 1;
  data_class[NODE_DAMPING] = NODE;
  data_required[NODE_DAMPING] = NODE;

  strcpy(name[NODE_DOF],"node_dof");
  type[NODE_DOF] = DOUBLE_PRECISION;
  data_length[NODE_DOF] = nuknwn;
  version_all[NODE_DOF] = 1;
  data_class[NODE_DOF] = NODE;
  data_required[NODE_DOF] = NODE;

  strcpy(name[NODE_DOF_CALCUL],"node_dof_calcul");
  type[NODE_DOF_CALCUL] = DOUBLE_PRECISION;
  data_length[NODE_DOF_CALCUL] = MCALCUL;
  fixed_length[NODE_DOF_CALCUL] = 0;
  version_all[NODE_DOF_CALCUL] = 1;
  data_class[NODE_DOF_CALCUL] = NODE;
  data_required[NODE_DOF_CALCUL] = NODE;

  strcpy(name[NODE_DOF_START_REFINED],"node_dof_start_refined");
  type[NODE_DOF_START_REFINED] = DOUBLE_PRECISION;
  data_length[NODE_DOF_START_REFINED] = nuknwn;
  version_all[NODE_DOF_START_REFINED] = 1;
  data_class[NODE_DOF_START_REFINED] = NODE;
  data_required[NODE_DOF_START_REFINED] = NODE;

  strcpy(name[NODE_DOF_TMP],"node_dof_tmp");
  type[NODE_DOF_TMP] = DOUBLE_PRECISION;
  data_length[NODE_DOF_TMP] = nuknwn;
  external[NODE_DOF_TMP] = 0;
  data_class[NODE_DOF_TMP] = NODE;
  data_required[NODE_DOF_TMP] = NODE;

  strcpy(name[NODE_EIGEN],"node_eigen");
  type[NODE_EIGEN] = DOUBLE_PRECISION;
  data_length[NODE_EIGEN] = DATA_ITEM_SIZE;
  fixed_length[NODE_EIGEN] = 0;
  version_all[NODE_EIGEN] = 1;
  data_class[NODE_EIGEN] = NODE;
  data_required[NODE_EIGEN] = NODE;

  strcpy(name[NODE_ELEMENT],"node_element");
  type[NODE_ELEMENT] = INTEGER;
  data_length[NODE_ELEMENT] = 0; // set run-time
  external[NODE_ELEMENT] = 0;
  fixed_length[NODE_ELEMENT] = 0;
  version_all[NODE_ELEMENT] = 1;
  data_class[NODE_ELEMENT] = NODE;
  data_required[NODE_ELEMENT] = NODE;

  strcpy(name[NODE_LHSIDE],"node_lhside");
  type[NODE_LHSIDE] = DOUBLE_PRECISION;
  data_length[NODE_LHSIDE] = npuknwn;
  external[NODE_LHSIDE] = 0;
  data_class[NODE_LHSIDE] = NODE;
  data_required[NODE_LHSIDE] = NODE;

  strcpy(name[NODE_MACRO_GENERATE],"node_macro_generate");
  type[NODE_MACRO_GENERATE] = INTEGER;
  data_length[NODE_MACRO_GENERATE] = 1;
  version_all[NODE_MACRO_GENERATE] = 1;
  data_class[NODE_MACRO_GENERATE] = NODE;
  data_required[NODE_MACRO_GENERATE] = NODE;
  external[NODE_MACRO_GENERATE] = 0;

  strcpy(name[NODE_MASS],"node_mass");
  type[NODE_MASS] = DOUBLE_PRECISION;
  data_length[NODE_MASS] = 1;
  version_all[NODE_MASS] = 1;
  data_class[NODE_MASS] = NODE;
  data_required[NODE_MASS] = NODE;

  strcpy(name[NODE_REMESH_ALLOWED],"node_remesh_allowed");
  strcpy(name[NODE_NEL],"node_nel");
  type[NODE_NEL] = INTEGER;
  data_length[NODE_NEL] = 1;
  external[NODE_NEL] = 0;
  version_all[NODE_NEL] = 1;
  data_class[NODE_NEL] = NODE;
  data_required[NODE_NEL] = NODE;

  strcpy(name[NODE_NODE],"node_node");
  type[NODE_NODE] = INTEGER;
  data_length[NODE_NODE] = 0; // set run-time
  external[NODE_NODE] = 0;
  fixed_length[NODE_NODE] = 0;
  version_all[NODE_NODE] = 1;
  data_class[NODE_NODE] = NODE;
  data_required[NODE_NODE] = NODE;

  strcpy(name[NODE_NONLOCAL],"node_nonlocal");
  type[NODE_NONLOCAL] = INTEGER;
  data_length[NODE_NONLOCAL] = NONLOCAL_ITEM_SIZE;
  external[NODE_NONLOCAL] = 0;
  version_all[NODE_NONLOCAL] = 1;
  fixed_length[NODE_NONLOCAL] = 0;
  data_class[NODE_NONLOCAL] = NODE;
  data_required[NODE_NONLOCAL] = NODE;

  strcpy(name[NODE_NONLOCAL_WEIGHT],"node_nonlocal_weight");
  type[NODE_NONLOCAL_WEIGHT] = DOUBLE_PRECISION;
  data_length[NODE_NONLOCAL_WEIGHT] = NONLOCAL_ITEM_SIZE;
  external[NODE_NONLOCAL_WEIGHT] = 0;
  version_all[NODE_NONLOCAL_WEIGHT] = 1;
  fixed_length[NODE_NONLOCAL_WEIGHT] = 0;
  data_class[NODE_NONLOCAL_WEIGHT] = NODE;
  data_required[NODE_NONLOCAL_WEIGHT] = NODE;

  strcpy(name[NODE_PHREATICLEVEL],"node_phreaticlevel");
  type[NODE_PHREATICLEVEL] = INTEGER;
  data_length[NODE_PHREATICLEVEL] = 1;
  version_all[NODE_PHREATICLEVEL] = 1;
  data_class[NODE_PHREATICLEVEL] = NODE;
  data_required[NODE_PHREATICLEVEL] = NODE;
  external[NODE_PHREATICLEVEL] = 0;

  strcpy(name[NODE_PRINT],"node_node_print");
  type[NODE_PRINT] = DOUBLE_PRECISION;
  data_length[NODE_PRINT] = ndim;
  version_all[NODE_PRINT] = 1;
  data_class[NODE_PRINT] = NODE;
  external[NODE_PRINT] = 0;

  strcpy(name[NODE_REMESH_ALLOWED],"node_remesh_allowed");
  type[NODE_REMESH_ALLOWED] = INTEGER;
  data_length[NODE_REMESH_ALLOWED] = ndim;
  version_all[NODE_REMESH_ALLOWED] = 1;
  external[NODE_REMESH_ALLOWED] = 0;
  data_class[NODE_REMESH_ALLOWED] = NODE;
  data_required[NODE_REMESH_ALLOWED] = NODE;

  strcpy(name[NODE_REMESH_VELOCITY],"node_remesh_velocity");
  type[NODE_REMESH_VELOCITY] = DOUBLE_PRECISION;
  data_length[NODE_REMESH_VELOCITY] = ndim;
  version_all[NODE_REMESH_VELOCITY] = 1;
  external[NODE_REMESH_VELOCITY] = 0;
  data_class[NODE_REMESH_VELOCITY] = NODE;
  data_required[NODE_REMESH_VELOCITY] = NODE;

  strcpy(name[NODE_RHSIDE],"node_rhside");
  type[NODE_RHSIDE] = DOUBLE_PRECISION;
  data_length[NODE_RHSIDE] = npuknwn;
  data_class[NODE_RHSIDE] = NODE;
  data_required[NODE_RHSIDE] = NODE;

  strcpy(name[NODE_RHSIDE_PRINT],"node_rhside_print");
  type[NODE_RHSIDE_PRINT] = DOUBLE_PRECISION;
  data_length[NODE_RHSIDE_PRINT] = npuknwn;
  version_all[NODE_RHSIDE_PRINT] = 1;
  external[NODE_RHSIDE_PRINT] = 0;
  data_required[NODE_RHSIDE_PRINT] = NODE;

  strcpy(name[NODE_START_REFINED],"node_start_refined");
  type[NODE_START_REFINED] = DOUBLE_PRECISION;
  data_length[NODE_START_REFINED] = ndim;
  version_all[NODE_START_REFINED] = 1;
  data_class[NODE_START_REFINED] = NODE;
  data_required[NODE_START_REFINED] = NODE;

  strcpy(name[NODE_SET],"node_set");
  type[NODE_SET] = INTEGER;
  data_length[NODE_SET] = 1;
  version_all[NODE_SET] = 1;
  data_class[NODE_SET] = NODE;

  strcpy(name[NODE_STIFFNESS],"node_stiffness");
  type[NODE_STIFFNESS] = DOUBLE_PRECISION;
  data_length[NODE_STIFFNESS] = ndim;
  version_all[NODE_STIFFNESS] = 1;
  data_class[NODE_STIFFNESS] = NODE;
  data_required[NODE_STIFFNESS] = NODE;

  strcpy(name[NONE],"none" );

  strcpy(name[NONLOCAL_ELEMENT_INFO],"nonlocal_element_info");
  type[NONLOCAL_ELEMENT_INFO] = DOUBLE_PRECISION;
  data_length[NONLOCAL_ELEMENT_INFO] = 1+npointmax*ndim+npointmax+2;
  version_all[NONLOCAL_ELEMENT_INFO] = 0;
  external[NONLOCAL_ELEMENT_INFO] = 0;
  data_class[NONLOCAL_ELEMENT_INFO] = ELEMENT;
  data_required[NONLOCAL_ELEMENT_INFO] = ELEMENT;

  strcpy(name[NORMAL],"normal");

  strcpy(name[NOTHING],"nothing");

  strcpy(name[NUMBER_ITERATIONS],"number_iterations");
  type[NUMBER_ITERATIONS] = INTEGER;
  version_all[NUMBER_ITERATIONS] = 1;
  data_length[NUMBER_ITERATIONS] = 1;
  no_index[NUMBER_ITERATIONS] = 1;
  external[NUMBER_ITERATIONS] = 0;
  data_class[NUMBER_ITERATIONS] = NUMBER_ITERATIONS;

  strcpy(name[OPTIONS_CONVECTION],"options_convection");
  type[OPTIONS_CONVECTION] = INTEGER;
  data_length[OPTIONS_CONVECTION] = 1;
  no_index[OPTIONS_CONVECTION] = 1;

  strcpy(name[OPTIONS_ELEMENT_DOF],"options_element_dof");
  type[OPTIONS_ELEMENT_DOF] = INTEGER;
  data_length[OPTIONS_ELEMENT_DOF] = 1;
  no_index[OPTIONS_ELEMENT_DOF] = 1;

  strcpy(name[OPTIONS_ELEMENTLOOP],"options_elementloop");
  type[OPTIONS_ELEMENTLOOP] = INTEGER;
  data_length[OPTIONS_ELEMENTLOOP] = 1;
  no_index[OPTIONS_ELEMENTLOOP] = 1;

  strcpy(name[OPTIONS_INERTIA],"options_inertia");
  type[OPTIONS_INERTIA] = INTEGER;
  data_length[OPTIONS_INERTIA] = 1;
  no_index[OPTIONS_INERTIA] = 1;

  strcpy(name[OPTIONS_MATRIX_GROUP],"options_matrix_group");
  type[OPTIONS_MATRIX_GROUP] = INTEGER;
  data_length[OPTIONS_MATRIX_GROUP] = 1;
  no_index[OPTIONS_MATRIX_GROUP] = 1;
  external[OPTIONS_MATRIX_GROUP] = 0;

  strcpy(name[OPTIONS_MATRIX_LENGTH],"options_matrix_length");
  type[OPTIONS_MATRIX_LENGTH] = INTEGER;
  data_length[OPTIONS_MATRIX_LENGTH] = 1;
  no_index[OPTIONS_MATRIX_LENGTH] = 1;

  strcpy(name[OPTIONS_MESH],"options_mesh");
  type[OPTIONS_MESH] = INTEGER;
  data_length[OPTIONS_MESH] = ndim;
  no_index[OPTIONS_MESH] = 1;

  strcpy(name[OPTIONS_NONLOCAL],"options_nonlocal");
  type[OPTIONS_NONLOCAL] = DOUBLE_PRECISION;
  data_length[OPTIONS_NONLOCAL] = 1;
  no_index[OPTIONS_NONLOCAL] = 1;

  strcpy(name[OPTIONS_NONLOCAL_SOFTVAR],"options_nonlocal_softvar");
  type[OPTIONS_NONLOCAL_SOFTVAR] = DOUBLE_PRECISION;
  data_length[OPTIONS_NONLOCAL_SOFTVAR] = 1;
  no_index[OPTIONS_NONLOCAL_SOFTVAR] = 1;

  strcpy(name[OPTIONS_PROCESSORS],"options_processors");
  type[OPTIONS_PROCESSORS] = INTEGER;
  data_length[OPTIONS_PROCESSORS] = 1;
  no_index[OPTIONS_PROCESSORS] = 1;

  strcpy(name[OPTIONS_RELAXATION],"options_relaxation");
  type[OPTIONS_RELAXATION] = DOUBLE_PRECISION;
  data_length[OPTIONS_RELAXATION] = nprinc;
  no_index[OPTIONS_RELAXATION] = 1;

  strcpy(name[OPTIONS_RESIDUEFACTOR],"options_residuefactor");
  type[OPTIONS_RESIDUEFACTOR] = DOUBLE_PRECISION;
  data_length[OPTIONS_RESIDUEFACTOR] = nprinc;
  no_index[OPTIONS_RESIDUEFACTOR] = 1;

  strcpy(name[OPTIONS_SKIP_GRAVITY],"options_skip_gravity");
  type[OPTIONS_SKIP_GRAVITY] = INTEGER;
  data_length[OPTIONS_SKIP_GRAVITY] = 1;
  no_index[OPTIONS_SKIP_GRAVITY] = 1;

  strcpy(name[OPTIONS_SKIP_GROUNDFLOW_NONLINEAR],"options_skip_groundflow_nonlinear");
  type[OPTIONS_SKIP_GROUNDFLOW_NONLINEAR] = INTEGER;
  data_length[OPTIONS_SKIP_GROUNDFLOW_NONLINEAR] = 1;
  no_index[OPTIONS_SKIP_GROUNDFLOW_NONLINEAR] = 1;

  strcpy(name[OPTIONS_SKIP_PLASTICITY],"options_skip_plasticity");
  type[OPTIONS_SKIP_PLASTICITY] = INTEGER;
  data_length[OPTIONS_SKIP_PLASTICITY] = 1;
  no_index[OPTIONS_SKIP_PLASTICITY] = 1;

  strcpy(name[OPTIONS_SOLVER],"options_solver");
  type[OPTIONS_SOLVER] = INTEGER;
  data_length[OPTIONS_SOLVER] = 1;
  no_index[OPTIONS_SOLVER] = 1;

  strcpy(name[OPTIONS_SOLVER_BICG_ERROR],"options_solver_bicg_error");
  type[OPTIONS_SOLVER_BICG_ERROR] = DOUBLE_PRECISION;
  data_length[OPTIONS_SOLVER_BICG_ERROR] = 1;
  no_index[OPTIONS_SOLVER_BICG_ERROR] = 1;

  strcpy(name[OPTIONS_SOLVER_BICG_ERROR_MINIMUM],"options_solver_bicg_error_minimum");
  type[OPTIONS_SOLVER_BICG_ERROR_MINIMUM] = DOUBLE_PRECISION;
  data_length[OPTIONS_SOLVER_BICG_ERROR_MINIMUM] = 1;
  no_index[OPTIONS_SOLVER_BICG_ERROR_MINIMUM] = 1;

  strcpy(name[OPTIONS_STABILIZATION],"options_stabilization");
  type[OPTIONS_STABILIZATION] = INTEGER;
  data_length[OPTIONS_STABILIZATION] = 1;
  no_index[OPTIONS_STABILIZATION] = 1;

  strcpy(name[PHIMOB],"phimob" );

  strcpy(name[POINT_MATERI_DIFFUSION],"point_materi_diffusion");
  type[POINT_MATERI_DIFFUSION] = DOUBLE_PRECISION;
  data_length[POINT_MATERI_DIFFUSION] = ndim;
  version_all[POINT_MATERI_DIFFUSION] = 1;
  external[POINT_MATERI_DIFFUSION] = 0;

  strcpy(name[POINT_MATERI_DIFFUSION_PREVIOUS],"point_materi_diffusion_previous");
  type[POINT_MATERI_DIFFUSION_PREVIOUS] = DOUBLE_PRECISION;
  data_length[POINT_MATERI_DIFFUSION_PREVIOUS] = ndim;
  version_all[POINT_MATERI_DIFFUSION_PREVIOUS] = 1;
  external[POINT_MATERI_DIFFUSION_PREVIOUS] = 0;

  strcpy(name[POSITIVE],"positive");

  strcpy(name[POST],"post");

  strcpy(name[POST_CALCUL],"post_calcul");
  type[POST_CALCUL] = INTEGER;
  data_length[POST_CALCUL] = DATA_ITEM_SIZE;
  no_index[POST_CALCUL] = 1;
  fixed_length[POST_CALCUL] = 0;
  data_class[POST_CALCUL] = POST;

  strcpy(name[POST_CALCUL_SCAL_VEC_MAT],"post_calcul_scal_vec_mat");
  type[POST_CALCUL_SCAL_VEC_MAT] = INTEGER;
  data_length[POST_CALCUL_SCAL_VEC_MAT] = DATA_ITEM_SIZE;
  fixed_length[POST_CALCUL_SCAL_VEC_MAT] = 0;
  no_index[POST_CALCUL_SCAL_VEC_MAT] = 1;
  external[POST_CALCUL_SCAL_VEC_MAT] = 0;
  data_class[POST_CALCUL_SCAL_VEC_MAT] = POST;

  strcpy(name[POST_CALCUL_UNKNOWN_OPERAT],"post_calcul_unknown_operat");
  type[POST_CALCUL_UNKNOWN_OPERAT] = INTEGER;
  data_length[POST_CALCUL_UNKNOWN_OPERAT] = DATA_ITEM_SIZE;
  fixed_length[POST_CALCUL_UNKNOWN_OPERAT] = 0;
  no_index[POST_CALCUL_UNKNOWN_OPERAT] = 1;
  external[POST_CALCUL_UNKNOWN_OPERAT] = 0;
  data_class[POST_CALCUL_UNKNOWN_OPERAT] = POST;

  strcpy(name[POST_ERROR_ITEM],"post_error_item");
  type[POST_ERROR_ITEM] = INTEGER;
  data_length[POST_ERROR_ITEM] = 3;
  data_class[POST_ERROR_ITEM] = POST;

  strcpy(name[POST_ERROR_MESH1],"post_error_mesh1");
  type[POST_ERROR_MESH1] = DOUBLE_PRECISION;
  data_length[POST_ERROR_MESH1] = 1;
  external[POST_ERROR_MESH1] = 0;
  data_class[POST_ERROR_MESH1] = POST;

  strcpy(name[POST_ERROR_MESH2],"post_error_mesh2");
  type[POST_ERROR_MESH2] = DOUBLE_PRECISION;
  data_length[POST_ERROR_MESH2] = 1;
  external[POST_ERROR_MESH2] = 0;
  data_class[POST_ERROR_MESH2] = POST;

  strcpy(name[POST_ERROR_RESULT],"post_error_result");
  type[POST_ERROR_RESULT] = DOUBLE_PRECISION;
  data_length[POST_ERROR_RESULT] = 1;
  data_class[POST_ERROR_RESULT] = POST;

  strcpy(name[POST_GLOBAL],"post_global");
  type[POST_GLOBAL] = INTEGER;
  data_length[POST_GLOBAL] = DATA_ITEM_SIZE;
  data_class[POST_GLOBAL] = POST;
  fixed_length[POST_GLOBAL] = 0;
  no_index[POST_GLOBAL] = 1;

  strcpy(name[POST_INTEGRATE],"post_integrate");
  type[POST_INTEGRATE] = INTEGER;
  data_length[POST_INTEGRATE] = 3;
  data_class[POST_INTEGRATE] = POST;

  strcpy(name[POST_INTEGRATE_RESULT],"post_integrate_result");
  type[POST_INTEGRATE_RESULT] = DOUBLE_PRECISION;
  data_length[POST_INTEGRATE_RESULT] = 1;
  data_class[POST_INTEGRATE_RESULT] = POST;
                                         
  strcpy(name[POST_LINE],"post_line");
  type[POST_LINE] = DOUBLE_PRECISION;
  data_length[POST_LINE] = 2*ndim;
  data_class[POST_LINE] = POST;

  strcpy(name[POST_LINE_DOF],"post_line_dof");
  type[POST_LINE_DOF] = DOUBLE_PRECISION;
  data_length[POST_LINE_DOF] = nuknwn;
  data_class[POST_LINE_DOF] = POST;

  strcpy(name[POST_LINE_DOF_CALCUL],"post_line_dof_calcul");
  type[POST_LINE_DOF_CALCUL] = DOUBLE_PRECISION;
  data_length[POST_LINE_DOF_CALCUL] = MCALCUL;
  fixed_length[POST_LINE_DOF_CALCUL] = 0;
  data_class[POST_LINE_DOF_CALCUL] = POST;

  strcpy(name[POST_LINE_MOMENT],"post_line_moment");
  type[POST_LINE_MOMENT] = INTEGER;
  data_length[POST_LINE_MOMENT] = 1;
  data_class[POST_LINE_MOMENT] = POST;

  strcpy(name[POST_LINE_N],"post_line_n");
  type[POST_LINE_N] = INTEGER;
  data_length[POST_LINE_N] = 1;
  data_class[POST_LINE_N] = POST;

  strcpy(name[POST_LINE_OPERAT],"post_line_operat");
  type[POST_LINE_OPERAT] = INTEGER;
  data_length[POST_LINE_OPERAT] = 1;
  data_class[POST_LINE_OPERAT] = POST;

  strcpy(name[POST_NODE],"post_node");
  type[POST_NODE] = INTEGER;
  data_length[POST_NODE] = 4;
  fixed_length[POST_NODE] = 0;
  data_class[POST_NODE] = POST;

  strcpy(name[POST_NODE_RESULT],"post_node_result");
  type[POST_NODE_RESULT] = DOUBLE_PRECISION;
  data_length[POST_NODE_RESULT] = DATA_ITEM_SIZE;
  fixed_length[POST_NODE_RESULT] = 0;
  data_class[POST_NODE_RESULT] = POST;

  strcpy(name[POST_NODE_RHSIDE_FIXED],"post_node_rhside_fixed");
  type[POST_NODE_RHSIDE_FIXED] = DOUBLE_PRECISION;
  data_length[POST_NODE_RHSIDE_FIXED] = npuknwn;
  data_class[POST_NODE_RHSIDE_FIXED] = POST;
  no_index[POST_NODE_RHSIDE_FIXED] = 1;

  strcpy(name[POST_NODE_RHSIDE_FREE],"post_node_rhside_free");
  type[POST_NODE_RHSIDE_FREE] = DOUBLE_PRECISION;
  data_length[POST_NODE_RHSIDE_FREE] = npuknwn;
  data_class[POST_NODE_RHSIDE_FREE] = POST;
  no_index[POST_NODE_RHSIDE_FREE] = 1;

  strcpy(name[POST_NODE_RHSIDE_RATIO],"post_node_rhside_ratio");
  type[POST_NODE_RHSIDE_RATIO] = DOUBLE_PRECISION;
  data_length[POST_NODE_RHSIDE_RATIO] = 1;
  data_class[POST_NODE_RHSIDE_RATIO] = POST;
  no_index[POST_NODE_RHSIDE_RATIO] = 1;

  strcpy(name[POST_NODE_RHSIDE_RATIO_UNKNOWNTYPES],"post_node_rhside_ratio_unknowntypes");
  type[POST_NODE_RHSIDE_RATIO_UNKNOWNTYPES] = INTEGER;
  data_length[POST_NODE_RHSIDE_RATIO_UNKNOWNTYPES] = DATA_ITEM_SIZE;
  data_class[POST_NODE_RHSIDE_RATIO_UNKNOWNTYPES] = POST;
  fixed_length[POST_NODE_RHSIDE_RATIO_UNKNOWNTYPES] = 0;
  no_index[POST_NODE_RHSIDE_RATIO_UNKNOWNTYPES] = 1;

  strcpy(name[POST_POINT],"post_point");
  type[POST_POINT] = DOUBLE_PRECISION;
  data_length[POST_POINT] = ndim;
  data_class[POST_POINT] = POST;

  strcpy(name[POST_POINT_DOF],"post_point_dof");
  type[POST_POINT_DOF] = DOUBLE_PRECISION;
  data_length[POST_POINT_DOF] = nuknwn;
  data_class[POST_POINT_DOF] = POST;

  strcpy(name[POST_POINT_DOF_CALCUL],"post_point_dof_calcul");
  type[POST_POINT_DOF_CALCUL] = DOUBLE_PRECISION;
  data_length[POST_POINT_DOF_CALCUL] = MCALCUL;
  fixed_length[POST_POINT_DOF_CALCUL] = 0;
  data_class[POST_POINT_DOF_CALCUL] = POST;

  strcpy(name[POST_QUADRILATERAL],"post_quadrilateral");
  type[POST_QUADRILATERAL] = DOUBLE_PRECISION;
  data_length[POST_QUADRILATERAL] = 4*ndim;
  data_class[POST_QUADRILATERAL] = POST;

  strcpy(name[POST_QUADRILATERAL_DOF],"post_quadrilateral_dof");
  type[POST_QUADRILATERAL_DOF] = DOUBLE_PRECISION;
  data_length[POST_QUADRILATERAL_DOF] = nuknwn;
  data_class[POST_QUADRILATERAL_DOF] = POST;

  strcpy(name[POST_QUADRILATERAL_DOF_CALCUL],"post_quadrilateral_dof_calcul");
  type[POST_QUADRILATERAL_DOF_CALCUL] = DOUBLE_PRECISION;
  data_length[POST_QUADRILATERAL_DOF_CALCUL] = DATA_ITEM_SIZE;
  fixed_length[POST_QUADRILATERAL_DOF_CALCUL] = 0;
  data_class[POST_QUADRILATERAL_DOF_CALCUL] = POST;

  strcpy(name[POST_QUADRILATERAL_N],"post_quadrilateral_n");
  type[POST_QUADRILATERAL_N] = INTEGER;
  data_length[POST_QUADRILATERAL_N] = 1;
  data_class[POST_QUADRILATERAL_N] = POST;

  strcpy(name[PREONLY],"preonly");

  strcpy(name[PRINT],"print");

  strcpy(name[PRINT_ARITHMETIC],"print_arithmetic");
  type[PRINT_ARITHMETIC] = INTEGER;
  data_length[PRINT_ARITHMETIC] = 1;
  no_index[PRINT_ARITHMETIC] = 1;

  strcpy(name[PRINT_CONTROL],"print_control");
  type[PRINT_CONTROL] = INTEGER;
  data_length[PRINT_CONTROL] = 1;
  no_index[PRINT_CONTROL] = 1;

  strcpy(name[PRINT_DEFINE],"print_define");
  type[PRINT_DEFINE] = INTEGER;
  data_length[PRINT_DEFINE] = 1;
  no_index[PRINT_DEFINE] = 1;

  strcpy(name[PRINT_FAILURE],"print_failure");
  type[PRINT_FAILURE] = INTEGER;
  data_length[PRINT_FAILURE] = 1;
  data_class[PRINT_FAILURE] = PRINT;
  no_index[PRINT_FAILURE] = 1;

  strcpy(name[PRINT_FILTER],"print_filter");
  type[PRINT_FILTER] = INTEGER;
  data_length[PRINT_FILTER] = DATA_ITEM_SIZE;
  fixed_length[PRINT_FILTER] = 0;
  data_class[PRINT_FILTER] = PRINT;

  strcpy(name[PRINT_LASTDATABASE],"print_lastdatabase");
  type[PRINT_LASTDATABASE] = INTEGER;
  data_length[PRINT_LASTDATABASE] = 1;
  data_class[PRINT_LASTDATABASE] = PRINT;
  no_index[PRINT_LASTDATABASE] = 1;

  strcpy(name[PRINT_SOLVER],"print_solver");
  type[PRINT_SOLVER] = INTEGER;
  data_length[PRINT_SOLVER] = 1;
  no_index[PRINT_SOLVER] = 1;

  strcpy(name[PRINT_WHERE],"print_where");
  type[PRINT_WHERE] = INTEGER;
  data_length[PRINT_WHERE] = 1;
  data_class[PRINT_WHERE] = PRINT;
  no_index[PRINT_WHERE] = 1;

  strcpy(name[PRISM6],"prism6");

  strcpy(name[PRIVAL],"prival");

  strcpy(name[PRIVEC],"privec");

  strcpy(name[P_COARSEN],"p_coarsen");

  strcpy(name[P_REFINEMENT],"p_refinement");

  strcpy(name[PUT],"put");

  strcpy(name[QUAD4],"quad4");

  strcpy(name[QUAD9],"quad9");

  strcpy(name[QUAD16],"quad16");

  strcpy(name[RA],"ra");

  strcpy(name[RECTANGLE],"rectangle" );

  strcpy(name[RESIDUE],"residue");

  strcpy(name[RESTART],"restart");

  strcpy(name[RICHARDSON],"richardson");

  strcpy(name[ROTATION_X_AXIS],"rotation_x_axis");

  strcpy(name[ROTATION_Y_AXIS],"rotation_y_axis");

  strcpy(name[ROTATION_Z_AXIS],"rotation_z_axis");

  strcpy(name[SCALAR],"scalar" );

  strcpy(name[SEPARATE],"separate");

  strcpy(name[SEPARATE_INDEX],"separate_index");

  strcpy(name[SEPARATE_SEQUENTIAL],"separate_sequential");

  strcpy(name[SHELL],"shell");

  strcpy(name[SIZEDEV],"sizedev");
	
  strcpy(name[MISES],"mises");

  strcpy(name[SIZETOT],"sizetot");

  strcpy(name[SLES],"sles");

  strcpy(name[SLIDE_FRICTION],"slide_friction");
  type[SLIDE_FRICTION] = DOUBLE_PRECISION;
  data_length[SLIDE_FRICTION] = 1;
  data_class[SLIDE_FRICTION] = SLIDE;

  strcpy(name[SLIDE_GEOMETRY],"slide_geometry");
  type[SLIDE_GEOMETRY] = INTEGER;
  data_length[SLIDE_GEOMETRY] = 2;
  data_class[SLIDE_GEOMETRY] = SLIDE;

  strcpy(name[SLIDE_PENALTY],"slide_penalty");
  type[SLIDE_PENALTY] = DOUBLE_PRECISION;
  data_length[SLIDE_PENALTY] = 1;
  data_class[SLIDE_PENALTY] = SLIDE;

  strcpy(name[SOR],"sor");

  strcpy(name[SPHERE],"sphere");

  strcpy(name[SPRING],"spring");

  strcpy(name[SPRING1],"spring1");

  strcpy(name[SPRING2],"spring2");

  strcpy(name[STATIC],"static" );

  strcpy(name[STEP],"step");

  strcpy(name[STRESS],"stress");

  strcpy(name[SUM],"sum");

  strcpy(name[TARGET],"target");

  strcpy(name[TARGET_ITEM],"target_item");
  type[TARGET_ITEM] = INTEGER;
  data_length[TARGET_ITEM] = 3;
  data_class[TARGET_ITEM] = TARGET;
  data_required[TARGET_ITEM] = TARGET_VALUE;

  strcpy(name[TARGET_VALUE],"target_value");
  type[TARGET_VALUE] = DOUBLE_PRECISION;
  data_length[TARGET_VALUE] = 2;
  data_class[TARGET_VALUE] = TARGET;
  data_required[TARGET_VALUE] = TARGET_ITEM;

  strcpy(name[TCQMR],"tcqmr");

  strcpy(name[TENDON],"tendon");
  type[TENDON] = DOUBLE_PRECISION;
  data_length[TENDON] = 2*ndim + 1;
  data_class[TENDON] = TENDON;

  strcpy(name[TENDON_ELASTI],"tendon_elasti");
  type[TENDON_ELASTI] = DOUBLE_PRECISION;
  data_length[TENDON_ELASTI] = 1;
  data_class[TENDON_ELASTI] = TENDON;
  data_required[TENDON_ELASTI] = TENDON;

  strcpy(name[TENDON_EXPANSION],"tendon_expansion");
  type[TENDON_EXPANSION] = DOUBLE_PRECISION;
  data_length[TENDON_EXPANSION] = 1;
  data_class[TENDON_EXPANSION] = TENDON;
  data_required[TENDON_EXPANSION] = TENDON;

  strcpy(name[TENDON_PLASTI],"tendon_plasti");
  type[TENDON_PLASTI] = DOUBLE_PRECISION;
  data_length[TENDON_PLASTI] = 1;
  data_class[TENDON_PLASTI] = TENDON;
  data_required[TENDON_PLASTI] = TENDON;

  strcpy(name[TENDON_SPLIT],"tendon_split");
  type[TENDON_SPLIT] = DOUBLE_PRECISION;
  data_length[TENDON_SPLIT] = 2;
  external[TENDON_SPLIT] = 0;
  data_class[TENDON_SPLIT] = TENDON;
  data_required[TENDON_SPLIT] = TENDON;

  strcpy(name[TENDON_SPLIT_ELEMENT],"tendon_split_element");
  type[TENDON_SPLIT_ELEMENT] = INTEGER;
  data_length[TENDON_SPLIT_ELEMENT] = MTENDON;
  fixed_length[TENDON_SPLIT_ELEMENT] = 0;
  external[TENDON_SPLIT_ELEMENT] = 0;
  data_class[TENDON_SPLIT_ELEMENT] = TENDON;
  data_required[TENDON_SPLIT_ELEMENT] = TENDON;

  strcpy(name[TENDON_STRESS],"tendon_stress");
  type[TENDON_STRESS] = DOUBLE_PRECISION;
  data_length[TENDON_STRESS] = 1;
  data_class[TENDON_STRESS] = TENDON;
  data_required[TENDON_STRESS] = TENDON;

  strcpy(name[TENDON_STRESS_TIME],"tendon_stress_time");
  type[TENDON_STRESS_TIME] = DOUBLE_PRECISION;
  data_length[TENDON_STRESS_TIME] = DATA_ITEM_SIZE;
  fixed_length[TENDON_STRESS_TIME] = 0;
  data_class[TENDON_STRESS_TIME] = TENDON;
  data_required[TENDON_STRESS_TIME] = TENDON;

  strcpy(name[TET4],"tet4");

  strcpy(name[TET10],"tet10");

  strcpy(name[TFQMR],"tfqmr");

  strcpy(name[THERMAL],"thermal");

  strcpy(name[TIME],"time");

  strcpy(name[TIME_AT_START],"time_at_start");
  type[TIME_AT_START] = INTEGER;
  data_length[TIME_AT_START] = 1;
  external[TIME_AT_START] = 0;
  data_class[TIME_AT_START] = TIME;

  strcpy(name[TIME_CURRENT],"time_current");
  type[TIME_CURRENT] = DOUBLE_PRECISION;
  data_length[TIME_CURRENT] = 1;
  no_index[TIME_CURRENT] = 1;
  version_all[TIME_CURRENT] = 1;
  data_class[TIME_CURRENT] = TIME;

  strcpy(name[TIME_CALCULATION],"time_calculation");
  type[TIME_CALCULATION] = INTEGER;
  data_length[TIME_CALCULATION] = 1;
  no_index[TIME_CALCULATION] = 1;
  data_class[TIME_CALCULATION] = TIME;

  strcpy(name[TIME_NEW],"time_new");
  type[TIME_NEW] = DOUBLE_PRECISION;
  external[TIME_NEW] = 0;
  data_length[TIME_NEW] = 1;
  no_index[TIME_NEW] = 1;
  version_all[TIME_NEW] = 1;
  data_class[TIME_NEW] = TIME;

  strcpy(name[TIME_OLD],"time_old");
  type[TIME_OLD] = DOUBLE_PRECISION;
  external[TIME_OLD] = 0;
  data_length[TIME_OLD] = 1;
  no_index[TIME_OLD] = 1;
  version_all[TIME_OLD] = 1;
  data_class[TIME_OLD] = TIME;

  strcpy(name[TO],"to");

  strcpy(name[TOTAL],"total");

  strcpy(name[TOTAL_LINEAR],"total_linear");

  strcpy(name[TOTAL_PIOLA],"total_piola");

  strcpy(name[TRIA3],"tria3");

  strcpy(name[TRIA6],"tria6");

  strcpy(name[TRUSS],"truss");

  strcpy(name[TRUSSBEAM],"trussbeam");

  strcpy(name[UNIFORM],"uniform");

  strcpy(name[UPDATED],"updated");

  strcpy(name[UPDATED_WITHOUT_ROTATION],"updated_without_rotation");

  strcpy(name[USE],"use" );

  strcpy(name[USER],"user" );

  strcpy(name[VALUE],"value");

  strcpy(name[VECTOR],"vector");

  strcpy(name[VOLUME_FACTOR],"volume_factor");
  type[VOLUME_FACTOR] = DOUBLE_PRECISION;
  data_length[VOLUME_FACTOR] = DATA_ITEM_SIZE;
  fixed_length[VOLUME_FACTOR] = 0;
  no_index[VOLUME_FACTOR] = 1;
  data_class[VOLUME_FACTOR] = VOLUME;

  strcpy(name[WAVE],"wave");

  strcpy(name[WAVE_SCALAR],"wave_scalar");

  strcpy(name[WAVE_FSCALAR],"wave_fscalar");

  strcpy(name[YES],"yes");

  strcpy(name[X],"x");

  strcpy(name[Y],"y");

  strcpy(name[Z],"z");

  for ( idat=0; idat<MDAT; idat++ ) {
    if ( data_length[idat]<1 ) data_length[idat] = 1;
  }

  idat = LAST_DUMMY;
  for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
    iuknwn = ipuknwn*nder;
    if      ( dof_type[iuknwn]==-BEAM_ROTATION ) {
      if ( iuknwn==rot_indx ) n = 0;
      n++;
      if ( ndim==2 ) {
        assert( n==1 );
        strcpy( basename, "rotz"  );
      }
      else {
        assert( ndim==3 );
        if      ( n==1 ) strcpy( basename, "rotx"  );
        else if ( n==2 ) strcpy( basename, "roty"  );
        else if ( n==3 ) strcpy( basename, "rotz"  );
      }
    }              
    else if ( dof_type[iuknwn]==-CONDIF_TEMPERATURE ) 
      strcpy( basename, "temp"  );
    else if ( dof_type[iuknwn]==-GROUNDFLOW_PRESSURE ) 
      strcpy( basename, "pres" );
    else if ( dof_type[iuknwn]==-GROUNDFLOW_VELOCITY ) {
      if ( iuknwn==gvel_indx ) n = 0;
      n++;
      if      ( n==1 ) strcpy( basename, "gvelx"  );
      else if ( n==2 ) strcpy( basename, "gvely"  );
      else if ( n==3 ) strcpy( basename, "gvelz"  );
    }
    else if ( dof_type[iuknwn]==-MATERI_DAMAGE ) 
      strcpy( basename, "dam" );
    else if ( dof_type[iuknwn]==-MATERI_DENSITY ) 
      strcpy( basename, "dens" );
    else if ( dof_type[iuknwn]==-MATERI_DIFFUSION ) 
      strcpy( basename, "diff" );
    else if ( dof_type[iuknwn]==-MATERI_DISPLACEMENT ) {
      if ( iuknwn==dis_indx ) n = 0;
      n++;
      if      ( n==1 ) strcpy( basename, "disx"  );
      else if ( n==2 ) strcpy( basename, "disy"  );
      else if ( n==3 ) strcpy( basename, "disz"  );
    }
    else if ( dof_type[iuknwn]==-MATERI_HISTORY_VARIABLES ) {
      if ( iuknwn==hisv_indx ) n = 0;
      strcpy( basename, "hisv" );
      long_to_a( n, str );
      strcat( basename, str );
      n++;
    }
    else if ( dof_type[iuknwn]==-MATERI_MAXWELL_STRESS ) {
      if ( iuknwn==mstres_indx ) { 
        m = 1; 
        n = 0; 
      }
      n++;
      if ( n>6 ) { 
        m++; 
        n = 1; 
      }
      if      ( n==1 ) strcpy( tmpname, "msigxx" );
      else if ( n==2 ) strcpy( tmpname, "msigxy" );
      else if ( n==3 ) strcpy( tmpname, "msigxz" );
      else if ( n==4 ) strcpy( tmpname, "msigyy" );
      else if ( n==5 ) strcpy( tmpname, "msigyz" );
      else if ( n==6 ) strcpy( tmpname, "msigzz" );
      long_to_a( m, basename );
      strcat( basename, tmpname );
    }
    else if ( dof_type[iuknwn]==-MATERI_PLASTI_F )
      strcpy( basename, "f" );
    else if ( dof_type[iuknwn]==-MATERI_PLASTI_F_NONLOCAL )
      strcpy( basename, "fn" );
    else if ( dof_type[iuknwn]==-MATERI_PLASTI_INCREMENTAL_SUBSTEPS )
      strcpy( basename, "substeps" );
    else if ( dof_type[iuknwn]==-MATERI_PLASTI_KAPPA )
      strcpy( basename, "kap" );
    else if ( dof_type[iuknwn]==-MATERI_PLASTI_RHO ) {
      if ( iuknwn==rho_indx ) n = 0;
      n++;
      if      ( n==1 ) strcpy( basename, "rhoxx" );
      else if ( n==2 ) strcpy( basename, "rhoxy" );
      else if ( n==3 ) strcpy( basename, "rhoxz" );
      else if ( n==4 ) strcpy( basename, "rhoyy" );
      else if ( n==5 ) strcpy( basename, "rhoyz" );
      else if ( n==6 ) strcpy( basename, "rhozz" );
    }
    else if ( dof_type[iuknwn]==-MATERI_PLASTI_SOFTVAR_LOCAL )
      strcpy( basename, "softvar_loc" );
    else if ( dof_type[iuknwn]==-MATERI_PLASTI_SOFTVAR_NONLOCAL )
      strcpy( basename, "softvar_nonloc" );
    else if ( dof_type[iuknwn]==-MATERI_STRAINENERGY ) 
      strcpy( basename, "ener" );
    else if ( dof_type[iuknwn]==-MATERI_STRAIN_ELASTI ) {
      if ( iuknwn==epe_indx ) n = 0;
      n++;
      if      ( n==1 ) strcpy( basename, "epexx" );
      else if ( n==2 ) strcpy( basename, "epexy" );
      else if ( n==3 ) strcpy( basename, "epexz" );
      else if ( n==4 ) strcpy( basename, "epeyy" );
      else if ( n==5 ) strcpy( basename, "epeyz" );
      else if ( n==6 ) strcpy( basename, "epezz" );
    }
    else if ( dof_type[iuknwn]==-MATERI_STRAIN_INTERGRANULAR ) {
      if ( iuknwn==epi_indx ) n = 0;
      n++;
      if      ( n==1 ) strcpy( basename, "epixx" );
      else if ( n==2 ) strcpy( basename, "epixy" );
      else if ( n==3 ) strcpy( basename, "epixz" );
      else if ( n==4 ) strcpy( basename, "epiyy" );
      else if ( n==5 ) strcpy( basename, "epiyz" );
      else if ( n==6 ) strcpy( basename, "epizz" );
    }
    else if ( dof_type[iuknwn]==-MATERI_STRAIN_PLASTI ) {
      if ( iuknwn==epp_indx ) n = 0;
      n++;
      if      ( n==1 ) strcpy( basename, "eppxx" );
      else if ( n==2 ) strcpy( basename, "eppxy" );
      else if ( n==3 ) strcpy( basename, "eppxz" );
      else if ( n==4 ) strcpy( basename, "eppyy" );
      else if ( n==5 ) strcpy( basename, "eppyz" );
      else if ( n==6 ) strcpy( basename, "eppzz" );
    }
    else if ( dof_type[iuknwn]==-MATERI_STRAIN_TOTAL ) {
      if ( iuknwn==ept_indx ) n = 0;
      n++;
      if      ( n==1 ) strcpy( basename, "eptxx" );
      else if ( n==2 ) strcpy( basename, "eptxy" );
      else if ( n==3 ) strcpy( basename, "eptxz" );
      else if ( n==4 ) strcpy( basename, "eptyy" );
      else if ( n==5 ) strcpy( basename, "eptyz" );
      else if ( n==6 ) strcpy( basename, "eptzz" );
    }
    else if ( dof_type[iuknwn]==-MATERI_STRESS ) {
      if ( iuknwn==stres_indx ) n = 0;
      n++;
      if      ( n==1 ) strcpy( basename, "sigxx" );
      else if ( n==2 ) strcpy( basename, "sigxy" );
      else if ( n==3 ) strcpy( basename, "sigxz" );
      else if ( n==4 ) strcpy( basename, "sigyy" );
      else if ( n==5 ) strcpy( basename, "sigyz" );
      else if ( n==6 ) strcpy( basename, "sigzz" );
    }
    else if ( dof_type[iuknwn]==-MATERI_VELOCITY ) {
      if ( iuknwn==vel_indx ) n = 0;
      n++;
      if      ( n==1 ) strcpy( basename, "velx"  );
      else if ( n==2 ) strcpy( basename, "vely"  );
      else if ( n==3 ) strcpy( basename, "velz"  );
    }
    else if ( dof_type[iuknwn]==-MATERI_VELOCITY_INTEGRATED ) {
      if ( iuknwn==veli_indx ) n = 0;
      n++;
      if      ( n==1 ) strcpy( basename, "velix"  );
      else if ( n==2 ) strcpy( basename, "veliy"  );
      else if ( n==3 ) strcpy( basename, "veliz"  );
    }
    else if ( dof_type[iuknwn]==-MATERI_VOID_FRACTION ) 
      strcpy( basename, "void" );
    else if ( dof_type[iuknwn]==-MATERI_WORK ) 
      strcpy( basename, "work" );
    else if ( dof_type[iuknwn]==-MAXWELL_E ) {
      if ( iuknwn==maxe_indx ) n = 0;
      n++;
      if      ( n==1 ) strcpy( basename, "ex"  );
      else if ( n==2 ) strcpy( basename, "ey"  );
      else if ( n==3 ) strcpy( basename, "ez"  );
    }                                                                   
    else if ( dof_type[iuknwn]==-MAXWELL_EI ) {
      if ( iuknwn==maxei_indx ) n = 0;
      n++;
      if      ( n==1 ) strcpy( basename, "eix"  );
      else if ( n==2 ) strcpy( basename, "eiy"  );
      else if ( n==3 ) strcpy( basename, "eiz"  );
    }                                                                   
    else if ( dof_type[iuknwn]==-MAXWELL_ER ) {
      if ( iuknwn==maxer_indx ) n = 0;
      n++;
      if      ( n==1 ) strcpy( basename, "erx"  );
      else if ( n==2 ) strcpy( basename, "ery"  );
      else if ( n==3 ) strcpy( basename, "erz"  );
    }                                                                   
    else if ( dof_type[iuknwn]==-MAXWELL_FE ) {
      if ( iuknwn==maxfe_indx ) n = 0;
      n++;
      if      ( n==1 ) strcpy( basename, "fex"  );
      else if ( n==2 ) strcpy( basename, "fey"  );
      else if ( n==3 ) strcpy( basename, "fez"  );
    }                                                                   
    else if ( dof_type[iuknwn]==-RESIDUE ) 
      strcpy( basename, "res" );
    else if ( dof_type[iuknwn]==-WAVE_SCALAR ) 
      strcpy( basename, "scal"  );
    else if ( dof_type[iuknwn]==-WAVE_FSCALAR ) 
      strcpy( basename, "fscal"  );

    idat++; assert( idat<MDAT );
    strcpy( name[idat], basename );
    dof_label[iuknwn] = -idat;

    if ( derivatives ) {
      for ( idim=0; idim<ndim; idim++ ) {
        iuknwn++;
        if      ( idim==0 ) strcpy( str, "x" );
        else if ( idim==1 ) strcpy( str, "y" );
        else if ( idim==2 ) strcpy( str, "z" );
        strcat( str, basename );
        idat++; strcpy( name[idat], str );
        dof_label[iuknwn] = -idat;
      }
      iuknwn++;
      strcpy( str, "t" );
      strcat( str, basename );
      idat++; strcpy( name[idat], str );
      dof_label[iuknwn] = -idat;
    }
  }

}

long int db( long int idat, long int index, long int *ival,
  double *dval, long int &length, long int version, long int task )

{
  long int i=0, l=0, data_ptr=0, data_number=0;

  data_number = labs(idat);
  if ( version<0 || version>=MVERSION ||
       (!version_all[data_number]&&version!=VERSION_NORMAL) ) {
    pri( "Version failure for version ", version );
    db_error( idat, index );
  }

  if      ( task==GET ) {
    if ( !db_active_index(data_number,index,version) ) db_error( idat, index );
    length = db_len( data_number, index, version );
  }
  else if ( task==GET_IF_EXISTS ) {
    if ( !db_active_index(data_number,index,version) ) {
      length = 0;
      return 0;
    }
    length = db_len( data_number, index, version );
  }
  else if ( task==GET_AND_CHECK ) {
    if ( !db_active_index(data_number,index,version) ) db_error( idat, index );
    l = db_len( data_number, index, version );
    if ( l!=length ) db_error( idat, index );
  }
  else if ( task==PUT ) {
    if ( index<0 )
      db_error( idat, index );
    if ( index>0 && db_no_index(idat) )
      db_error( idat, index );
    if ( length<1 )
      db_error( idat, index );
    if ( length>db_data_length(data_number) ) {
      pri( "Length too small of ", db_name(data_number) );
      pri( "The length is ", db_data_length(data_number) );
      pri( "It should become at least", length );
      pri( "Increase it and recompile." );
      exit(TN_EXIT_STATUS);
    }
    if ( db_fixed_length( idat ) && length!=db_data_length(data_number) )
      db_error( idat, index );
    db_allocate( data_number, index, version, MAXIMAL );
  }
  data_ptr = db_data_length(data_number) * index;

  if ( db_type(data_number)==INTEGER ) {
    for ( i=0; i<length; i++ ) {
      if ( task==PUT )
        int_data[data_number][version][data_ptr+i] = ival[i];
      else
        ival[i] = int_data[data_number][version][data_ptr+i];
    }
  }
  else {
    for ( i=0; i<length; i++ ) {
      if ( task==PUT )
        dbl_data[data_number][version][data_ptr+i] = dval[i];
      else
        dval[i] = dbl_data[data_number][version][data_ptr+i];
    }
  }

  if ( task==PUT ) {
    if ( length<db_data_length(data_number) ) {
      if ( db_type(data_number)==INTEGER )
        int_data[data_number][version][data_ptr+length] = LONG_MIN;
      else
        dbl_data[data_number][version][data_ptr+length] = DBL_MAX;
    }
  }

  return 1;

}


long int db_active_index( long int idat, long int index, long int version )

{
  long int result=1, max=-1, data_number=0, data_ptr=0;

  data_number = labs(idat);
  max = max_index[data_number][version];
  data_ptr = db_data_length(data_number) * index;

  if      ( index<0 || index>max )
    result = 0;
  else if ( type[data_number]==DOUBLE_PRECISION )
    result = (dbl_data[data_number][version][data_ptr]!=DBL_MAX);
  else {
    result = (int_data[data_number][version][data_ptr]!=LONG_MIN);
  }
  return result;

}

void db_allocate( long int idat, long int index, long int version, long int task )

{
  long int n=0, n_old=0, max_old=0, data_number=0,
    max=0, length=0, increase=0, *int_old=0;
  double *dbl_old=0;

  data_number = labs(idat);
  if ( version<0 || version>=MVERSION ||
       (!version_all[data_number]&&version!=VERSION_NORMAL) ) {
    pri( "Version failure for version ", version );
    db_error( idat, index );
  }

  db_max_index( data_number, max_old, version, GET );
  if ( index<=max_old ) return;

  if ( parallel_active ) {
    pri( "Program error detected for data", -idat );
    pri( "Data is allocated in a parallel loop.");
    pri( "Please report this error with the input file.");
    exit_tn_on_error();
  }

  max = index;
  if ( task==MAXIMAL && !db_no_index(data_number) ) {
    increase = index/10; // a bit extra for future use, heuristic
    if ( increase<1 ) increase = 1;
    max += increase; 
  }
  db_max_index( data_number, max, version, PUT );

  length = db_data_length( data_number );
  n     = (1+max    ) * length;
  n_old = (1+max_old) * length;
  if ( db_type(data_number)==INTEGER ) {
    if ( max_old>=0 ) int_old = int_data[data_number][version];
    int_data[data_number][version] = get_new_int(n);
    array_set( int_data[data_number][version], LONG_MIN, n );
    if ( max_old>=0 ) {
      array_move( int_old, int_data[data_number][version], n_old );
      delete[] int_old;
    }
  }
  else if ( db_type(data_number)==DOUBLE_PRECISION ) {
    if ( max_old>=0 ) dbl_old = dbl_data[data_number][version];
    dbl_data[data_number][version] = get_new_dbl(n);
    array_set( dbl_data[data_number][version], DBL_MAX, n );
    if ( max_old>=0 ) {
      array_move( dbl_old, dbl_data[data_number][version], n_old );
      delete[] dbl_old;
    }
  }
  else {
    db_error( data_number, version );
  }

}

void db_allocate_class( long int cl, long int index, long int version )

{
  long int idat=0, ic=0;

  ic = labs( cl );

  for ( idat=0; idat<MDAT; idat++ ) {
    if ( db_data_class(idat)==ic && db_version( idat, version ) )
      db_allocate( idat, index, version, -MINIMAL );
  }
}

long int db_partialname( long int idat, char *str )

{
  long int return_value=0;
  char *name;

  name = db_name( idat );
  if ( strstr(name,str)!=NULL ) return_value = 1;

  return return_value;
}

long int db_partialname_any( const char *str )

{
  long int idat=0, max=0, return_value=0;
  char *name;

  for ( idat=0; idat<MDAT; idat++ ) {
    name = db_name( idat );
    if ( strstr(name,str)!=NULL ) {
      db_highest_index( idat, max, VERSION_NORMAL );
      if ( max>=0 ) return_value = 1;
    }
  }

  return return_value;
}

long int db_partialname_any_index( const char *str, long int index )

{
  long int idat=0, return_value=0;
  char *name;

  for ( idat=0; idat<MDAT; idat++ ) {
    name = db_name( idat );
    if ( strstr(name,str)!=NULL ) {
      if ( db_active_index( idat, index, VERSION_NORMAL ) )
        return_value = 1;
    }
  }

  return return_value;
}

long int db_data_class( long int idat )

{
  long int data_number=0;

  data_number = labs(idat);
  return data_class[data_number];
}

long int db_data_required( long int idat )

{
  long int data_number=0;

  data_number = labs(idat);
  return data_required[data_number];
}

void db_close( )

{
  long int idat=0, version=0, data_number=0;

  for ( version=0; version<MVERSION; version++ ) {
    for ( idat=0; idat<MDAT; idat++ ) {
     data_number = labs(idat);
      if ( version==0 || version_all[data_number] ) db_delete( idat, version );
    }
  }
}

void db_copy( long int idat, long int jdat, long int version )

{
  long int n=0, max=0, data_numberi=0, data_numberj=0;

  data_numberi = labs(idat);
  data_numberj = labs(jdat);
  if ( version<0 || version>=MVERSION ||
       (!version_all[data_numberi]&&version!=VERSION_NORMAL) ) {
    pri( "Version failure for version ", version );
    db_error( idat, -1 );
  }
  if ( version<0 || version>=MVERSION ||
       (!version_all[data_numberj]&&version!=VERSION_NORMAL) ) {
    pri( "Version failure for version ", version );
    db_error( jdat, -1 );
  }

  db_max_index( data_numberi, max, version, GET );
  if ( max>=0 ) {
    if ( db_data_length(data_numberi)!=db_data_length(data_numberj) )
      db_error( data_numberj, -1 );
    db_delete( data_numberj, version );
    db_allocate( data_numberj, max, version, MINIMAL );
    n = (1+max)*db_data_length(data_numberi);
    if ( db_type(data_numberi)==INTEGER )
      array_move( int_data[data_numberi][version], int_data[data_numberj][version], n );
    else
      array_move( dbl_data[data_numberi][version], dbl_data[data_numberj][version], n );
  }

}

long int db_data_length( long int idat )

{

  long int data_number=0;

  data_number = labs(idat);

  return data_length[data_number];

}

void db_data_length_put( long int idat, long int length )

{

  long int data_number=0, version=0;

  data_number = labs(idat);

  for ( version=0; version<MVERSION; version++ ) {
    if ( db_version( idat, version ) ) db_delete( idat, version );
  }
  data_length[data_number] = length;

}

double *db_dbl( long int idat, long int version )

{
  long int data_number=0, ldum=0;
  double *ptr=NULL;

  data_number = labs(idat);
  if ( version<0 || version>=MVERSION ||
       (!version_all[data_number]&&version!=VERSION_NORMAL) ) {
    pri( "Version failure for version ", version );
    db_error( idat, -1 );
  }

  if ( db_max_index(data_number,ldum,version,GET)<0 ) db_error( idat, -1 );
  if ( db_type(data_number)!=DOUBLE_PRECISION ) db_error( idat, -1 );

  ptr = &dbl_data[data_number][version][0];

  return ptr;

}

double *db_dbl( long int idat, long int index, long int version )

{
  long int data_ptr=0, data_number=0;
  double *ptr=NULL;

  data_number = labs(idat);
  if ( version<0 || version>=MVERSION ||
       (!version_all[data_number]&&version!=VERSION_NORMAL) ) {
    pri( "Version failure for version ", version );
    db_error( idat, index );
  }

  if ( !db_active_index(data_number,index,version) ) db_error( idat, index );
  if ( db_type(data_number)!=DOUBLE_PRECISION ) db_error( idat, index );

  data_ptr = data_length[data_number] * index;
  ptr = &dbl_data[data_number][version][data_ptr];

  return ptr;

}

void db_delete( long int idat, long int version )

{
  long int data_number=0, max=-1, ldum=0;

  data_number = labs(idat);
  if ( version<0 || version>=MVERSION ||
       (!version_all[data_number]&&version!=VERSION_NORMAL) ) {
    pri( "Version failure for version ", version );
    db_error( idat, -1 );
  }

  if ( db_max_index(idat,ldum,version,GET)>=0 ) {
    if      ( db_type(data_number)==INTEGER )
      delete[] int_data[data_number][version];
    else
      delete[] dbl_data[data_number][version];
    db_max_index(idat,max,version,PUT);
  }

}

void db_delete_index( long int idat, long int index, long int version )

{
  long int length=0, data_number=0;

  data_number = labs(idat);
  if ( version<0 || version>=MVERSION ||
       (!version_all[data_number]&&version!=VERSION_NORMAL) ) {
    pri( "Version failure for version ", version );
    db_error( idat, index );
  }

  if ( db_active_index(idat,index,version) ) {
    length = db_data_length( idat );
    if      ( db_type(idat)==INTEGER )
      array_set( db_int(idat,index,version), LONG_MIN, length );
    else
      array_set( db_dbl(idat,index,version), DBL_MAX, length );
  }
}

void db_error( long int idat, long int index )

{
  pri( "Error detected for data item ", db_name(idat) );
  if ( index>=0 && !db_no_index(idat) )
    pri( "Error detected for record ", index );

  exit(TN_EXIT_STATUS);
}


long int db_external( long int idat )

{
  long int data_number = 0;

  data_number = labs(idat);

  return external[data_number];

}

long int db_fixed_length( long int idat )

{
  long int data_number=0;

  data_number = labs(idat);

  return fixed_length[data_number];

}

void db_highest_index( long int idat, long int &max, long int version )

{
  long int index=0, data_number=0;

  data_number = labs(idat);
  if ( version<0 || version>=MVERSION ||
       (!version_all[data_number]&&version!=VERSION_NORMAL) ) {
    pri( "Version failure for version ", version );
    db_error( idat, -1 );
  }

  db_max_index( idat, max, version, GET );

  index = max;
  while ( index>=0 && !db_active_index(idat,index,version) ) index--;
  max = index;

}

long int *db_int( long int idat, long int version )

{
  long int data_number=0, ldum=0, *ptr=NULL;

  data_number = labs(idat);
  if ( version<0 || version>=MVERSION ||
       (!version_all[data_number]&&version!=VERSION_NORMAL) ) {
    pri( "Version failure for version ", version );
    db_error( idat, -1 );
  }

  if ( db_max_index(data_number,ldum,version,GET)<0 ) db_error( idat, -1 );
  if ( db_type(data_number)!=INTEGER ) db_error( idat, -1 );

  ptr = &int_data[data_number][version][0];

  return ptr;

}

long int *db_int( long int idat, long int index, long int version )

{
  long int data_ptr=0, data_number=0, *ptr=NULL;

  data_number = labs(idat);
  if ( version<0 || version>=MVERSION ||
       (!version_all[data_number]&&version!=VERSION_NORMAL) ) {
    pri( "Version failure for version ", version );
    db_error( idat, index );
  }

  if ( !db_active_index(data_number,index,version) ) db_error( idat, index );
  if ( db_type(data_number)!=INTEGER ) db_error( idat, index );

  data_ptr = data_length[data_number] * index;
  ptr = &int_data[data_number][version][data_ptr];

  return ptr;

}

long int db_len( long int idat, long int index, long int version )

{
  long int i=0, l=0, data_ptr=0, data_number=0, length=0, left=0, right=0;

  data_number = labs( idat );
  if ( version<0 || version>=MVERSION ||
       (!version_all[data_number]&&version!=VERSION_NORMAL) ) {
    pri( "Version failure for version ", version );
    db_error( idat, index );
  }

  l = data_length[data_number];
  if ( fixed_length[data_number] )
    length = l;
  else if ( db_active_index( idat, index, version ) ) {
    data_ptr = l * index;
    length = -1;
    left = 0;
    right = l;
    do {
      i = ( left + right ) / 2;
      if      ( left==right-1 )
        length = right;
      else if ( type[data_number]==INTEGER ) {
        if ( int_data[data_number][version][data_ptr+i]==LONG_MIN ) 
          right = i;
        else
          left = i;
      }
      else {
        if ( dbl_data[data_number][version][data_ptr+i]==DBL_MAX  ) 
          right = i;
        else
          left = i;
      }
    }
    while ( length<0 );
  }
  else {
    db_error( idat, index );
  }

  return length;
}


long int db_max_index( long int idat, long int &max, long int version, long int task )

{
  long int data_number=0;

  data_number = labs(idat);
  if ( version<0 || version>=MVERSION ||
       (!version_all[data_number]&&version!=VERSION_NORMAL) ) {
    pri( "Version failure for version ", version );
    db_error( idat, -1 );
  }

  assert( data_number>=0 || data_number<MDAT );

  if ( task==GET ) 
    max = max_index[data_number][version];
  else {
    assert( task==PUT );
    max_index[data_number][version] = max;
  }

  return max;
}

char *db_name( long int idat )

{
  long int data_number=0;

  data_number = labs(idat);
  return name[data_number];
}

long int db_no_index( long int idat )
{
  long int data_number=0;

  data_number = labs(idat);
  return no_index[data_number];
}

long int db_number( char str[] )

{
  long int data_number=0, found=-1;

  for ( data_number=0; data_number<MDAT && found<0; data_number++ ) {
    if ( !strcmp(str,db_name(data_number)) ) found = data_number;
  }

  return found;

}

long int db_print_only( long int idat )

{
  long int data_number=0, return_value=0;

  data_number = labs(idat);

  return_value = print_only[data_number];
  return return_value;
}

void db_set_dbl( long int jdat, long int version )

{
  long int k=0, max=0, length=0, index=0, data_ptr=0, data_numberi=0, 
    data_numberj=0, idat=0, ldum=0;

  assert( db_type(jdat)==DOUBLE_PRECISION );

  assert( version>=0 && version<MVERSION );
  idat = db_data_class( jdat );
  if ( db_max_index(idat,ldum,version,GET)<0 ) return;

  data_numberi = labs(idat);
  data_numberj = labs(jdat);
  if ( !version_all[data_numberi] ) assert( version==VERSION_NORMAL );

  db_max_index( data_numberi, max, version, GET ); 
  db_delete( data_numberj, version );
  db_allocate( data_numberj, max, version, MINIMAL );

  length = db_data_length(data_numberj);
  for ( index=0; index<=max; index++ ) {
    data_ptr = length * index;
    if ( db_active_index(data_numberi,index,version) ) {
      for ( k=0; k<length; k++ )
        dbl_data[data_numberj][version][data_ptr+k] = 0.;
    }
    else {
      for ( k=0; k<length; k++ ) {
        dbl_data[data_numberj][version][data_ptr+k] = DBL_MAX;
      }
    }
  }

}

void db_set_int( long int jdat, long int version )

{
  long int k=0, max=0, length=0, index=0, data_ptr=0, data_numberi=0, 
    data_numberj=0, idat=0, ldum=0;

  assert( db_type(jdat)==INTEGER );

  assert( version>=0 && version<MVERSION );
  idat = db_data_class( jdat );
  if ( db_max_index(idat,ldum,version,GET)<0 ) return;

  data_numberi = labs(idat);
  data_numberj = labs(jdat);
  if ( !version_all[data_numberi] ) assert( version==VERSION_NORMAL );

  db_max_index( data_numberi, max, version, GET );
  db_delete( data_numberj, version );
  db_allocate( data_numberj, max, version, MINIMAL );

  length = db_data_length(data_numberj);
  for ( index=0; index<=max; index++ ) {
    data_ptr = length * index;
    if ( db_active_index(data_numberi,index,version) ) {
      for ( k=0; k<length; k++ )
        int_data[data_numberj][version][data_ptr+k] = 0;
    }
    else {
      for ( k=0; k<length; k++ )
        int_data[data_numberj][version][data_ptr+k] = LONG_MIN;
    }
  }

}

long int db_type( long int idat )

{
  long int data_number=0;

  data_number = labs(idat);

  return type[data_number];

}

long int db_version( long int idat, long int version )

{
  long int data_number=0, return_value=0;

  data_number = labs(idat);
  assert( version>=0 && version<MVERSION );

  if      ( version_all[data_number] )
    return_value = 1;
  else
    return_value = (version==VERSION_NORMAL);

  return return_value;
}

void db_version_copy( long int version_from, long int version_to )

{
  long int idat=0;

  for ( idat=0; idat<MDAT; idat++ ) {
    if ( version_all[idat] )
      db_version_copy_data( idat, version_from, version_to );
  }

}

void db_version_copy_data( long int idat, long int version_from, long int version_to )

{
  long int max_from=0, max_to=0, max=0, n=0, *int_ptr_from=NULL, *int_ptr_to=NULL;
  double *dbl_ptr_from=NULL, *dbl_ptr_to=NULL;

  assert( version_all[idat] );
  assert( version_from>=0 && version_from<MVERSION );
  assert( version_to>=0 && version_to<MVERSION );

  db_max_index( idat, max_from, version_from, GET );
  db_max_index( idat, max_to, version_to, GET );
  if ( max_from>max_to ) 
    max = max_from;
  else 
    max = max_to;

  if ( max_from>=0 ) {
    db_delete( idat, version_to );
    if ( db_type(idat)==INTEGER ) {
      int_ptr_from = db_int( idat, version_from );
      db_allocate( idat, max, version_to, MINIMAL );
      int_ptr_to = db_int( idat, version_to );
      n = (1+max_from)*db_data_length(idat);
      array_move( int_ptr_from, int_ptr_to, n );
    }
    else { 
      assert( db_type(idat)==DOUBLE_PRECISION );
      dbl_ptr_from = db_dbl( idat, version_from );
      db_allocate( idat, max, version_to, MINIMAL );
      dbl_ptr_to = db_dbl( idat, version_to );
      n = (1+max_from)*db_data_length(idat);
      array_move( dbl_ptr_from, dbl_ptr_to, n );
    }
    db_max_index( idat, max, version_to, PUT );
  }

}

void db_version_delete( long int version )

{
  long int idat=0;

  assert( version>=0 && version<MVERSION );
  for ( idat=0; idat<MDAT; idat++ ) {
    if ( db_version(idat,version) ) db_delete( idat, version );
  }

}

