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
#include "tochnog_exceptions.h"

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
long int db_read[MDAT];        // set to 1 when item is read (GET/GET_IF_EXISTS)

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
  array_set( db_read, 0, MDAT );
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
  strcpy(name[MULTIPLY],"multiply" );

  strcpy(name[DIVIDE],"divide" );

  strcpy(name[MINUS],"minus" );

  strcpy(name[PLUS],"plus" );

  strcpy(name[ADD_ALWAYS],"add_always" );

  strcpy(name[ALL],"all" );

  strcpy(name[ANY],"any" );

  strcpy(name[ANY_BUT_NOT_ALL],"any_but_not_all" );

  strcpy(name[AREA],"area");

  strcpy(name[AREA_ELEMENT_GROUP],"area_element_group");
  type[AREA_ELEMENT_GROUP] = INTEGER;
  data_length[AREA_ELEMENT_GROUP] = 3;
  data_class[AREA_ELEMENT_GROUP] = AREA;

  strcpy(name[AREA_ELEMENT_GROUP_ELEMENT],"area_element_group_element");
  type[AREA_ELEMENT_GROUP_ELEMENT] = INTEGER;
  data_length[AREA_ELEMENT_GROUP_ELEMENT] = 1;
  data_class[AREA_ELEMENT_GROUP_ELEMENT] = AREA;
  data_required[AREA_ELEMENT_GROUP_ELEMENT] = AREA_ELEMENT_GROUP;

  strcpy(name[AREA_ELEMENT_GROUP_INTERFACE],"area_element_group_interface");
  type[AREA_ELEMENT_GROUP_INTERFACE] = INTEGER;
  data_length[AREA_ELEMENT_GROUP_INTERFACE] = 1;
  data_class[AREA_ELEMENT_GROUP_INTERFACE] = AREA;
  data_required[AREA_ELEMENT_GROUP_INTERFACE] = AREA_ELEMENT_GROUP;

  strcpy(name[AREA_ELEMENT_GROUP_NODE],"area_element_group_node");
  type[AREA_ELEMENT_GROUP_NODE] = INTEGER;
  data_length[AREA_ELEMENT_GROUP_NODE] = DATA_ITEM_SIZE;
  fixed_length[AREA_ELEMENT_GROUP_NODE] = 0;
  data_class[AREA_ELEMENT_GROUP_NODE] = AREA;
  data_required[AREA_ELEMENT_GROUP_NODE] = AREA_ELEMENT_GROUP;

  strcpy(name[AREA_ELEMENT_GROUP_TIME],"area_element_group_time");
  type[AREA_ELEMENT_GROUP_TIME] = INTEGER;
  data_length[AREA_ELEMENT_GROUP_TIME] = 1;
  data_class[AREA_ELEMENT_GROUP_TIME] = AREA;
  data_required[AREA_ELEMENT_GROUP_TIME] = AREA_ELEMENT_GROUP;

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
  // combination checked at runtime (area_element_group_sequence accepts
  // both the legacy elementgroup and the Professional element_group alias)

  strcpy(name[AREA_ELEMENT_GROUP_SEQUENCE_TIME],"area_element_group_sequence_time");
  type[AREA_ELEMENT_GROUP_SEQUENCE_TIME] = DOUBLE_PRECISION;
  data_length[AREA_ELEMENT_GROUP_SEQUENCE_TIME] = DATA_ITEM_SIZE;
  data_class[AREA_ELEMENT_GROUP_SEQUENCE_TIME] = AREA;
  fixed_length[AREA_ELEMENT_GROUP_SEQUENCE_TIME] = 0;
  // combination checked at runtime (area_element_group_sequence accepts
  // both the legacy elementgroup and the Professional element_group alias)

  strcpy(name[AREA_ELEMENT_GROUP_METHOD],"area_element_group_method");
  type[AREA_ELEMENT_GROUP_METHOD] = INTEGER;
  data_length[AREA_ELEMENT_GROUP_METHOD] = 1;
  data_class[AREA_ELEMENT_GROUP_METHOD] = AREA;
  data_required[AREA_ELEMENT_GROUP_METHOD] = AREA_ELEMENT_GROUP;

  strcpy(name[AREA_ELEMENT_GROUP_SEQUENCE_METHOD],"area_element_group_sequence_method");
  type[AREA_ELEMENT_GROUP_SEQUENCE_METHOD] = INTEGER;
  data_length[AREA_ELEMENT_GROUP_SEQUENCE_METHOD] = 1;
  data_class[AREA_ELEMENT_GROUP_SEQUENCE_METHOD] = AREA;
  data_required[AREA_ELEMENT_GROUP_SEQUENCE_METHOD] = AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP;

  strcpy(name[AREA_ELEMENT_GROUP_SEQUENCE_ELEMENT_GROUP],"area_element_group_sequence_element_group");
  type[AREA_ELEMENT_GROUP_SEQUENCE_ELEMENT_GROUP] = INTEGER;
  data_length[AREA_ELEMENT_GROUP_SEQUENCE_ELEMENT_GROUP] = DATA_ITEM_SIZE;
  data_class[AREA_ELEMENT_GROUP_SEQUENCE_ELEMENT_GROUP] = AREA;
  fixed_length[AREA_ELEMENT_GROUP_SEQUENCE_ELEMENT_GROUP] = 0;

  strcpy(name[AREA_ELEMENT_GROUP_SEQUENCE_GEOMETRY_METHOD],"area_element_group_sequence_geometry_method");
  type[AREA_ELEMENT_GROUP_SEQUENCE_GEOMETRY_METHOD] = INTEGER;
  data_length[AREA_ELEMENT_GROUP_SEQUENCE_GEOMETRY_METHOD] = 1;
  data_class[AREA_ELEMENT_GROUP_SEQUENCE_GEOMETRY_METHOD] = AREA;
  data_required[AREA_ELEMENT_GROUP_SEQUENCE_GEOMETRY_METHOD] = AREA_ELEMENT_GROUP_SEQUENCE_GEOMETRY;

  strcpy(name[AREA_ELEMENT_GROUP_SEQUENCE_INTERFACE],"area_element_group_sequence_interface");
  type[AREA_ELEMENT_GROUP_SEQUENCE_INTERFACE] = INTEGER;
  data_length[AREA_ELEMENT_GROUP_SEQUENCE_INTERFACE] = 1;
  data_class[AREA_ELEMENT_GROUP_SEQUENCE_INTERFACE] = AREA;
  data_required[AREA_ELEMENT_GROUP_SEQUENCE_INTERFACE] = AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP;

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

  strcpy(name[BOUNDA_ALTERNATE],"bounda_alternate");
  type[BOUNDA_ALTERNATE] = INTEGER;
  data_length[BOUNDA_ALTERNATE] = DATA_ITEM_SIZE;
  fixed_length[BOUNDA_ALTERNATE] = 0;
  data_class[BOUNDA_ALTERNATE] = BOUNDA;

  strcpy(name[BOUNDA_FORCE],"bounda_force");
  type[BOUNDA_FORCE] = INTEGER;
  data_length[BOUNDA_FORCE] = MBOUNDA;
  fixed_length[BOUNDA_FORCE] = 0;
  data_class[BOUNDA_FORCE] = BOUNDA;

  strcpy(name[BOUNDA_FACTOR],"bounda_factor");
  type[BOUNDA_FACTOR] = DOUBLE_PRECISION;
  data_length[BOUNDA_FACTOR] = DATA_ITEM_SIZE;
  fixed_length[BOUNDA_FACTOR] = 0;
  data_class[BOUNDA_FACTOR] = BOUNDA;

  strcpy(name[BOUNDA_FACTOR_PARABOLIC_X],"bounda_factor_parabolic_x");
  type[BOUNDA_FACTOR_PARABOLIC_X] = DOUBLE_PRECISION;
  data_length[BOUNDA_FACTOR_PARABOLIC_X] = DATA_ITEM_SIZE;
  fixed_length[BOUNDA_FACTOR_PARABOLIC_X] = 0;
  data_class[BOUNDA_FACTOR_PARABOLIC_X] = BOUNDA;

  strcpy(name[BOUNDA_FOUND],"bounda_found");
  type[BOUNDA_FOUND] = INTEGER;
  data_length[BOUNDA_FOUND] = 1;
  data_class[BOUNDA_FOUND] = BOUNDA;

  strcpy(name[BOUNDA_GEOMETRY_METHOD],"bounda_geometry_method");
  type[BOUNDA_GEOMETRY_METHOD] = INTEGER;
  data_length[BOUNDA_GEOMETRY_METHOD] = 1;
  data_class[BOUNDA_GEOMETRY_METHOD] = BOUNDA;

  strcpy(name[BOUNDA_DOF],"bounda_dof");
  type[BOUNDA_DOF] = INTEGER;
  data_length[BOUNDA_DOF] = MBOUNDA;
  fixed_length[BOUNDA_DOF] = 0;
  data_class[BOUNDA_DOF] = BOUNDA;

  strcpy(name[BOUNDA_DOF_CYLINDRICAL],"bounda_dof_cylindrical");
  type[BOUNDA_DOF_CYLINDRICAL] = DOUBLE_PRECISION;
  data_length[BOUNDA_DOF_CYLINDRICAL] = 6;
  data_class[BOUNDA_DOF_CYLINDRICAL] = BOUNDA;
  data_required[BOUNDA_DOF_CYLINDRICAL] = BOUNDA_DOF;

  strcpy(name[BOUNDA_DOF_RADIAL],"bounda_dof_radial");
  type[BOUNDA_DOF_RADIAL] = DOUBLE_PRECISION;
  data_length[BOUNDA_DOF_RADIAL] = 3;
  data_class[BOUNDA_DOF_RADIAL] = BOUNDA;
  data_required[BOUNDA_DOF_RADIAL] = BOUNDA_DOF;

  strcpy(name[BOUNDA_SINE],"bounda_sine");
  type[BOUNDA_SINE] = DOUBLE_PRECISION;
  data_length[BOUNDA_SINE] = DATA_ITEM_SIZE;
  fixed_length[BOUNDA_SINE] = 0;
  data_class[BOUNDA_SINE] = BOUNDA;

  strcpy(name[BOUNDA_CONSTANT],"bounda_constant");
  type[BOUNDA_CONSTANT] = INTEGER;
  data_length[BOUNDA_CONSTANT] = 1;
  data_class[BOUNDA_CONSTANT] = BOUNDA;
  data_required[BOUNDA_CONSTANT] = BOUNDA_UNKNOWN;

  strcpy(name[BOUNDA_NORMAL],"bounda_normal");
  type[BOUNDA_NORMAL] = DOUBLE_PRECISION;
  data_length[BOUNDA_NORMAL] = 3;
  data_class[BOUNDA_NORMAL] = BOUNDA;
  data_required[BOUNDA_NORMAL] = BOUNDA_UNKNOWN;

  strcpy(name[BOUNDA_TIME],"bounda_time");
  type[BOUNDA_TIME] = DOUBLE_PRECISION;
  data_length[BOUNDA_TIME] = DATA_ITEM_SIZE;
  fixed_length[BOUNDA_TIME] = 0;
  data_class[BOUNDA_TIME] = BOUNDA;

  strcpy(name[BOUNDA_TIME_FILE],"bounda_time_file");
  type[BOUNDA_TIME_FILE] = INTEGER;
  data_length[BOUNDA_TIME_FILE] = 1;
  data_class[BOUNDA_TIME_FILE] = BOUNDA;

  strcpy(name[BOUNDA_TIME_INCREMENT],"bounda_time_increment");
  type[BOUNDA_TIME_INCREMENT] = DOUBLE_PRECISION;
  data_length[BOUNDA_TIME_INCREMENT] = 1;
  data_class[BOUNDA_TIME_INCREMENT] = BOUNDA;
  data_required[BOUNDA_TIME_INCREMENT] = BOUNDA_TIME;

  strcpy(name[BOUNDA_TIME_OFFSET],"bounda_time_offset");
  type[BOUNDA_TIME_OFFSET] = DOUBLE_PRECISION;
  data_length[BOUNDA_TIME_OFFSET] = 1;
  data_class[BOUNDA_TIME_OFFSET] = BOUNDA;
  data_required[BOUNDA_TIME_OFFSET] = BOUNDA_TIME;

  strcpy(name[BOUNDA_TIME_FACTOR],"bounda_time_factor");
  type[BOUNDA_TIME_FACTOR] = DOUBLE_PRECISION;
  data_length[BOUNDA_TIME_FACTOR] = 1;
  data_class[BOUNDA_TIME_FACTOR] = BOUNDA;
  data_required[BOUNDA_TIME_FACTOR] = BOUNDA_TIME;

  strcpy(name[BOUNDA_TIME_ON_OFF],"bounda_time_on_off");
  type[BOUNDA_TIME_ON_OFF] = DOUBLE_PRECISION;
  data_length[BOUNDA_TIME_ON_OFF] = 2;
  data_class[BOUNDA_TIME_ON_OFF] = BOUNDA;
  data_required[BOUNDA_TIME_ON_OFF] = BOUNDA_TIME;

  strcpy(name[BOUNDA_TIME_UNTIL_FORCE],"bounda_time_until_force");
  type[BOUNDA_TIME_UNTIL_FORCE] = DOUBLE_PRECISION;
  data_length[BOUNDA_TIME_UNTIL_FORCE] = 2;
  data_class[BOUNDA_TIME_UNTIL_FORCE] = BOUNDA;
  data_required[BOUNDA_TIME_UNTIL_FORCE] = BOUNDA_TIME;

  // manual Professional 6.40: bounda_time_until_data index
  // data_item_name data_item_index data_item_number. Monitor the data
  // item; when it falls from start (6.41) the load of the bounda_time
  // record with the same index is reduced, reaching 0 when the monitor
  // reaches wanted. The reduction factor is quadratic:
  // ((monitor/first - wanted)/(start - wanted))^2, clamped [0,1],
  // where first is the initial value of the monitor (verified against
  // the Professional binary 25-10-2023: until1.dat E=1 and E=2 runs).
  strcpy(name[BOUNDA_TIME_UNTIL_DATA],"bounda_time_until_data");
  type[BOUNDA_TIME_UNTIL_DATA] = INTEGER;
  data_length[BOUNDA_TIME_UNTIL_DATA] = 3;
  fixed_length[BOUNDA_TIME_UNTIL_DATA] = 1;
  data_class[BOUNDA_TIME_UNTIL_DATA] = BOUNDA;
  data_required[BOUNDA_TIME_UNTIL_DATA] = BOUNDA_TIME;

  strcpy(name[BOUNDA_TIME_UNTIL_VALUE_MINIMUM],"bounda_time_until_value_minimum");
  type[BOUNDA_TIME_UNTIL_VALUE_MINIMUM] = DOUBLE_PRECISION;
  data_length[BOUNDA_TIME_UNTIL_VALUE_MINIMUM] = 2;
  fixed_length[BOUNDA_TIME_UNTIL_VALUE_MINIMUM] = 1;
  data_class[BOUNDA_TIME_UNTIL_VALUE_MINIMUM] = BOUNDA;
  data_required[BOUNDA_TIME_UNTIL_VALUE_MINIMUM] = BOUNDA_TIME_UNTIL_DATA;

  // records written to the .dbs database (same output semantics as the
  // Professional): bounda_time_until_first = initial monitor value,
  // bounda_time_until_used = reduction factor applied in the step.
  strcpy(name[BOUNDA_TIME_UNTIL_FIRST],"bounda_time_until_first");
  type[BOUNDA_TIME_UNTIL_FIRST] = DOUBLE_PRECISION;
  data_length[BOUNDA_TIME_UNTIL_FIRST] = 1;
  data_class[BOUNDA_TIME_UNTIL_FIRST] = BOUNDA;

  strcpy(name[BOUNDA_TIME_UNTIL_USED],"bounda_time_until_used");
  type[BOUNDA_TIME_UNTIL_USED] = DOUBLE_PRECISION;
  data_length[BOUNDA_TIME_UNTIL_USED] = 1;
  data_class[BOUNDA_TIME_UNTIL_USED] = BOUNDA;

  // manual Professional 6.41 variant used by the 02-08-2026 corpus
  // (validation_14_mesh): bounda_time_until_value index min max start.
  // Registered as a known keyword so the parser accepts it; the
  // consumption semantics could not be verified (the 25-10-2023 binary
  // rejects the record) and is documented as PENDING.
  strcpy(name[BOUNDA_TIME_UNTIL_VALUE],"bounda_time_until_value");
  type[BOUNDA_TIME_UNTIL_VALUE] = DOUBLE_PRECISION;
  data_length[BOUNDA_TIME_UNTIL_VALUE] = 3;
  fixed_length[BOUNDA_TIME_UNTIL_VALUE] = 1;
  data_class[BOUNDA_TIME_UNTIL_VALUE] = BOUNDA;
  data_required[BOUNDA_TIME_UNTIL_VALUE] = BOUNDA_TIME_UNTIL_DATA;

  strcpy(name[BOUNDA_TIME_UNITS],"bounda_time_units");
  type[BOUNDA_TIME_UNITS] = DOUBLE_PRECISION;
  data_length[BOUNDA_TIME_UNITS] = 2;
  data_class[BOUNDA_TIME_UNITS] = BOUNDA;
  data_required[BOUNDA_TIME_UNITS] = BOUNDA_TIME;

  strcpy(name[BOUNDA_TIME_USER],"bounda_time_user");
  type[BOUNDA_TIME_USER] = INTEGER;
  data_length[BOUNDA_TIME_USER] = 1;
  data_class[BOUNDA_TIME_USER] = BOUNDA;

  // bounda_time_smc/_offset/_units (manual Professional 6.42/6.43/6.44):
  // base acceleration read from an SMC (Strong Motion CD) file. The
  // records are REGISTERED and PARSED; the SMC file reader is not yet
  // implemented (PENDING). bounda_time_smc index -yes would need the
  // file <index>.smc next to the input file.
  strcpy(name[BOUNDA_TIME_SMC],"bounda_time_smc");
  type[BOUNDA_TIME_SMC] = INTEGER;
  data_length[BOUNDA_TIME_SMC] = 1;
  data_class[BOUNDA_TIME_SMC] = BOUNDA;

  strcpy(name[BOUNDA_TIME_SMC_OFFSET],"bounda_time_smc_offset");
  type[BOUNDA_TIME_SMC_OFFSET] = DOUBLE_PRECISION;
  data_length[BOUNDA_TIME_SMC_OFFSET] = 1;
  data_class[BOUNDA_TIME_SMC_OFFSET] = BOUNDA;

  strcpy(name[BOUNDA_TIME_SMC_UNITS],"bounda_time_smc_units");
  type[BOUNDA_TIME_SMC_UNITS] = DOUBLE_PRECISION;
  data_length[BOUNDA_TIME_SMC_UNITS] = 2;
  data_class[BOUNDA_TIME_SMC_UNITS] = BOUNDA;

  strcpy(name[BOUNDA_WATER],"bounda_water");
  type[BOUNDA_WATER] = INTEGER;
  data_length[BOUNDA_WATER] = 1;
  data_class[BOUNDA_WATER] = BOUNDA;
  data_required[BOUNDA_WATER] = BOUNDA_UNKNOWN;

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

  strcpy(name[CHANGE_DATAITEM_GEOMETRY],"change_dataitem_geometry");
  type[CHANGE_DATAITEM_GEOMETRY] = INTEGER;
  data_length[CHANGE_DATAITEM_GEOMETRY] = 2;
  data_required[CHANGE_DATAITEM_GEOMETRY] = CHANGE_DATAITEM;

  strcpy(name[CHANGE_DATAITEM_TIME],"change_dataitem_time");
  type[CHANGE_DATAITEM_TIME] = DOUBLE_PRECISION;
  data_length[CHANGE_DATAITEM_TIME] = DATA_ITEM_SIZE;
  fixed_length[CHANGE_DATAITEM_TIME] = 0;
  data_required[CHANGE_DATAITEM_TIME] = CHANGE_DATAITEM;

  strcpy(name[CHANGE_DATAITEM_TIME_DISCRETE],"change_dataitem_time_discrete");
  type[CHANGE_DATAITEM_TIME_DISCRETE] = INTEGER;
  data_length[CHANGE_DATAITEM_TIME_DISCRETE] = 1;
  data_required[CHANGE_DATAITEM_TIME_DISCRETE] = CHANGE_DATAITEM;

  strcpy(name[CHANGE_DATAITEM_TIME_METHOD],"change_dataitem_time_method");
  type[CHANGE_DATAITEM_TIME_METHOD] = INTEGER;
  data_length[CHANGE_DATAITEM_TIME_METHOD] = 1;
  data_required[CHANGE_DATAITEM_TIME_METHOD] = CHANGE_DATAITEM;

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

  strcpy(name[CHECK_DATA],"check_data");
  type[CHECK_DATA] = INTEGER;
  data_length[CHECK_DATA] = 1;
  no_index[CHECK_DATA] = 1;

  strcpy(name[CHECK_ELEMENT_NODE],"check_element_node");
  type[CHECK_ELEMENT_NODE] = INTEGER;
  data_length[CHECK_ELEMENT_NODE] = 1;
  no_index[CHECK_ELEMENT_NODE] = 1;

  strcpy(name[CHECK_ELEMENT_SHAPE],"check_element_shape");
  type[CHECK_ELEMENT_SHAPE] = DOUBLE_PRECISION;
  data_length[CHECK_ELEMENT_SHAPE] = 1;
  no_index[CHECK_ELEMENT_SHAPE] = 1;

  strcpy(name[CHECK_ERROR],"check_error");
  type[CHECK_ERROR] = INTEGER;
  data_length[CHECK_ERROR] = 1;
  no_index[CHECK_ERROR] = 1;

  strcpy(name[CHECK_MEMORY],"check_memory");
  type[CHECK_MEMORY] = INTEGER;
  data_length[CHECK_MEMORY] = 1;
  no_index[CHECK_MEMORY] = 1;

  strcpy(name[CHECK_MEMORY_USAGE],"check_memory_usage");
  type[CHECK_MEMORY_USAGE] = INTEGER;
  data_length[CHECK_MEMORY_USAGE] = 1;
  no_index[CHECK_MEMORY_USAGE] = 1;

  strcpy(name[CHECK_MEMORY_USAGE_RESULT],"check_memory_usage_result");
  type[CHECK_MEMORY_USAGE_RESULT] = DOUBLE_PRECISION;
  data_length[CHECK_MEMORY_USAGE_RESULT] = 1;
  no_index[CHECK_MEMORY_USAGE_RESULT] = 1;

  strcpy(name[CHECK_NAN],"check_nan");
  type[CHECK_NAN] = INTEGER;
  data_length[CHECK_NAN] = 1;
  no_index[CHECK_NAN] = 1;

  strcpy(name[CHECK_NUMBER],"check_number" );

  strcpy(name[CHECK_SOLVER],"check_solver");
  type[CHECK_SOLVER] = DOUBLE_PRECISION;
  data_length[CHECK_SOLVER] = 1;
  no_index[CHECK_SOLVER] = 1;

  strcpy(name[CHECK_TARGET],"check_target");
  type[CHECK_TARGET] = INTEGER;
  data_length[CHECK_TARGET] = 1;
  no_index[CHECK_TARGET] = 1;

  strcpy(name[CHECK_WARNING],"check_warning");
  type[CHECK_WARNING] = INTEGER;
  data_length[CHECK_WARNING] = 1;
  no_index[CHECK_WARNING] = 1;

  strcpy(name[CHECK_USED],"check_used");
  type[CHECK_USED] = INTEGER;
  data_length[CHECK_USED] = 1;
  no_index[CHECK_USED] = 1;

  strcpy(name[CIRCLE],"circle" );

  strcpy(name[CIRCLE_HOLLOW],"circle_hollow" );

  strcpy(name[COMPOSITE],"composite");

  strcpy(name[CONDIF],"condif");

  strcpy(name[CONDIF_CONVECTION],"condif_convection");
  type[CONDIF_CONVECTION] = DOUBLE_PRECISION;
  data_length[CONDIF_CONVECTION] = 2;
  data_class[CONDIF_CONVECTION] = CONDIF;

  strcpy(name[CONDIF_CONVECTION_EDGE_NORMAL],"condif_convection_edge_normal");
  type[CONDIF_CONVECTION_EDGE_NORMAL] = DOUBLE_PRECISION;
  data_length[CONDIF_CONVECTION_EDGE_NORMAL] = 2;
  data_class[CONDIF_CONVECTION_EDGE_NORMAL] = CONDIF;

  strcpy(name[CONDIF_CONVECTION_EDGE_NORMAL_ELEMENT],"condif_convection_edge_normal_element");
  type[CONDIF_CONVECTION_EDGE_NORMAL_ELEMENT] = INTEGER;
  data_length[CONDIF_CONVECTION_EDGE_NORMAL_ELEMENT] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_CONVECTION_EDGE_NORMAL_ELEMENT] = 0;
  data_class[CONDIF_CONVECTION_EDGE_NORMAL_ELEMENT] = CONDIF;
  data_required[CONDIF_CONVECTION_EDGE_NORMAL_ELEMENT] = CONDIF_CONVECTION_EDGE_NORMAL;

  strcpy(name[CONDIF_CONVECTION_EDGE_NORMAL_ELEMENT_GROUP],"condif_convection_edge_normal_element_group");
  type[CONDIF_CONVECTION_EDGE_NORMAL_ELEMENT_GROUP] = INTEGER;
  data_length[CONDIF_CONVECTION_EDGE_NORMAL_ELEMENT_GROUP] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_CONVECTION_EDGE_NORMAL_ELEMENT_GROUP] = 0;
  data_class[CONDIF_CONVECTION_EDGE_NORMAL_ELEMENT_GROUP] = CONDIF;
  data_required[CONDIF_CONVECTION_EDGE_NORMAL_ELEMENT_GROUP] = CONDIF_CONVECTION_EDGE_NORMAL;

  strcpy(name[CONDIF_CONVECTION_EDGE_NORMAL_ELEMENT_NODE],"condif_convection_edge_normal_element_node");
  type[CONDIF_CONVECTION_EDGE_NORMAL_ELEMENT_NODE] = INTEGER;
  data_length[CONDIF_CONVECTION_EDGE_NORMAL_ELEMENT_NODE] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_CONVECTION_EDGE_NORMAL_ELEMENT_NODE] = 0;
  data_class[CONDIF_CONVECTION_EDGE_NORMAL_ELEMENT_NODE] = CONDIF;
  data_required[CONDIF_CONVECTION_EDGE_NORMAL_ELEMENT_NODE] = CONDIF_CONVECTION_EDGE_NORMAL;

  strcpy(name[CONDIF_CONVECTION_EDGE_NORMAL_ELEMENT_SIDE],"condif_convection_edge_normal_element_side");
  type[CONDIF_CONVECTION_EDGE_NORMAL_ELEMENT_SIDE] = INTEGER;
  data_length[CONDIF_CONVECTION_EDGE_NORMAL_ELEMENT_SIDE] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_CONVECTION_EDGE_NORMAL_ELEMENT_SIDE] = 0;
  data_class[CONDIF_CONVECTION_EDGE_NORMAL_ELEMENT_SIDE] = CONDIF;
  data_required[CONDIF_CONVECTION_EDGE_NORMAL_ELEMENT_SIDE] = CONDIF_CONVECTION_EDGE_NORMAL;

  strcpy(name[CONDIF_CONVECTION_EDGE_NORMAL_GEOMETRY],"condif_convection_edge_normal_geometry");
  type[CONDIF_CONVECTION_EDGE_NORMAL_GEOMETRY] = INTEGER;
  data_length[CONDIF_CONVECTION_EDGE_NORMAL_GEOMETRY] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_CONVECTION_EDGE_NORMAL_GEOMETRY] = 0;
  data_class[CONDIF_CONVECTION_EDGE_NORMAL_GEOMETRY] = CONDIF;
  data_required[CONDIF_CONVECTION_EDGE_NORMAL_GEOMETRY] = CONDIF_CONVECTION_EDGE_NORMAL;

  strcpy(name[CONDIF_CONVECTION_EDGE_NORMAL_NODE],"condif_convection_edge_normal_node");
  type[CONDIF_CONVECTION_EDGE_NORMAL_NODE] = INTEGER;
  data_length[CONDIF_CONVECTION_EDGE_NORMAL_NODE] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_CONVECTION_EDGE_NORMAL_NODE] = 0;
  data_class[CONDIF_CONVECTION_EDGE_NORMAL_NODE] = CONDIF;
  data_required[CONDIF_CONVECTION_EDGE_NORMAL_NODE] = CONDIF_CONVECTION_EDGE_NORMAL;

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

  strcpy(name[CONDIF_RADIATION_EDGE_NORMAL],"condif_radiation_edge_normal");
  type[CONDIF_RADIATION_EDGE_NORMAL] = DOUBLE_PRECISION;
  data_length[CONDIF_RADIATION_EDGE_NORMAL] = 2;
  data_class[CONDIF_RADIATION_EDGE_NORMAL] = CONDIF;

  strcpy(name[CONDIF_RADIATION_EDGE_NORMAL_ELEMENT],"condif_radiation_edge_normal_element");
  type[CONDIF_RADIATION_EDGE_NORMAL_ELEMENT] = INTEGER;
  data_length[CONDIF_RADIATION_EDGE_NORMAL_ELEMENT] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_RADIATION_EDGE_NORMAL_ELEMENT] = 0;
  data_class[CONDIF_RADIATION_EDGE_NORMAL_ELEMENT] = CONDIF;
  data_required[CONDIF_RADIATION_EDGE_NORMAL_ELEMENT] = CONDIF_RADIATION_EDGE_NORMAL;

  strcpy(name[CONDIF_RADIATION_EDGE_NORMAL_ELEMENT_GROUP],"condif_radiation_edge_normal_element_group");
  type[CONDIF_RADIATION_EDGE_NORMAL_ELEMENT_GROUP] = INTEGER;
  data_length[CONDIF_RADIATION_EDGE_NORMAL_ELEMENT_GROUP] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_RADIATION_EDGE_NORMAL_ELEMENT_GROUP] = 0;
  data_class[CONDIF_RADIATION_EDGE_NORMAL_ELEMENT_GROUP] = CONDIF;
  data_required[CONDIF_RADIATION_EDGE_NORMAL_ELEMENT_GROUP] = CONDIF_RADIATION_EDGE_NORMAL;

  strcpy(name[CONDIF_RADIATION_EDGE_NORMAL_ELEMENT_NODE],"condif_radiation_edge_normal_element_node");
  type[CONDIF_RADIATION_EDGE_NORMAL_ELEMENT_NODE] = INTEGER;
  data_length[CONDIF_RADIATION_EDGE_NORMAL_ELEMENT_NODE] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_RADIATION_EDGE_NORMAL_ELEMENT_NODE] = 0;
  data_class[CONDIF_RADIATION_EDGE_NORMAL_ELEMENT_NODE] = CONDIF;
  data_required[CONDIF_RADIATION_EDGE_NORMAL_ELEMENT_NODE] = CONDIF_RADIATION_EDGE_NORMAL;

  strcpy(name[CONDIF_RADIATION_EDGE_NORMAL_ELEMENT_SIDE],"condif_radiation_edge_normal_element_side");
  type[CONDIF_RADIATION_EDGE_NORMAL_ELEMENT_SIDE] = INTEGER;
  data_length[CONDIF_RADIATION_EDGE_NORMAL_ELEMENT_SIDE] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_RADIATION_EDGE_NORMAL_ELEMENT_SIDE] = 0;
  data_class[CONDIF_RADIATION_EDGE_NORMAL_ELEMENT_SIDE] = CONDIF;
  data_required[CONDIF_RADIATION_EDGE_NORMAL_ELEMENT_SIDE] = CONDIF_RADIATION_EDGE_NORMAL;

  strcpy(name[CONDIF_RADIATION_EDGE_NORMAL_GEOMETRY],"condif_radiation_edge_normal_geometry");
  type[CONDIF_RADIATION_EDGE_NORMAL_GEOMETRY] = INTEGER;
  data_length[CONDIF_RADIATION_EDGE_NORMAL_GEOMETRY] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_RADIATION_EDGE_NORMAL_GEOMETRY] = 0;
  data_class[CONDIF_RADIATION_EDGE_NORMAL_GEOMETRY] = CONDIF;
  data_required[CONDIF_RADIATION_EDGE_NORMAL_GEOMETRY] = CONDIF_RADIATION_EDGE_NORMAL;

  strcpy(name[CONDIF_RADIATION_EDGE_NORMAL_NODE],"condif_radiation_edge_normal_node");
  type[CONDIF_RADIATION_EDGE_NORMAL_NODE] = INTEGER;
  data_length[CONDIF_RADIATION_EDGE_NORMAL_NODE] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_RADIATION_EDGE_NORMAL_NODE] = 0;
  data_class[CONDIF_RADIATION_EDGE_NORMAL_NODE] = CONDIF;
  data_required[CONDIF_RADIATION_EDGE_NORMAL_NODE] = CONDIF_RADIATION_EDGE_NORMAL;

  strcpy(name[CONDIF_RADIATION_GEOMETRY],"condif_radiation_geometry");
  type[CONDIF_RADIATION_GEOMETRY] = INTEGER;
  data_length[CONDIF_RADIATION_GEOMETRY] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_RADIATION_GEOMETRY] = 0;
  data_class[CONDIF_RADIATION_GEOMETRY] = CONDIF;
  data_required[CONDIF_RADIATION_GEOMETRY] = CONDIF_RADIATION;

  strcpy(name[CONDIF_HEAT_EDGE_NORMAL],"condif_heat_edge_normal");
  type[CONDIF_HEAT_EDGE_NORMAL] = DOUBLE_PRECISION;
  data_length[CONDIF_HEAT_EDGE_NORMAL] = 1;
  data_class[CONDIF_HEAT_EDGE_NORMAL] = CONDIF;

  strcpy(name[CONDIF_HEAT_EDGE_NORMAL_ELEMENT],"condif_heat_edge_normal_element");
  type[CONDIF_HEAT_EDGE_NORMAL_ELEMENT] = INTEGER;
  data_length[CONDIF_HEAT_EDGE_NORMAL_ELEMENT] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_HEAT_EDGE_NORMAL_ELEMENT] = 0;
  data_class[CONDIF_HEAT_EDGE_NORMAL_ELEMENT] = CONDIF;
  data_required[CONDIF_HEAT_EDGE_NORMAL_ELEMENT] = CONDIF_HEAT_EDGE_NORMAL;

  strcpy(name[CONDIF_HEAT_EDGE_NORMAL_ELEMENT_GROUP],"condif_heat_edge_normal_element_group");
  type[CONDIF_HEAT_EDGE_NORMAL_ELEMENT_GROUP] = INTEGER;
  data_length[CONDIF_HEAT_EDGE_NORMAL_ELEMENT_GROUP] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_HEAT_EDGE_NORMAL_ELEMENT_GROUP] = 0;
  data_class[CONDIF_HEAT_EDGE_NORMAL_ELEMENT_GROUP] = CONDIF;
  data_required[CONDIF_HEAT_EDGE_NORMAL_ELEMENT_GROUP] = CONDIF_HEAT_EDGE_NORMAL;

  strcpy(name[CONDIF_HEAT_EDGE_NORMAL_ELEMENT_NODE],"condif_heat_edge_normal_element_node");
  type[CONDIF_HEAT_EDGE_NORMAL_ELEMENT_NODE] = INTEGER;
  data_length[CONDIF_HEAT_EDGE_NORMAL_ELEMENT_NODE] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_HEAT_EDGE_NORMAL_ELEMENT_NODE] = 0;
  data_class[CONDIF_HEAT_EDGE_NORMAL_ELEMENT_NODE] = CONDIF;
  data_required[CONDIF_HEAT_EDGE_NORMAL_ELEMENT_NODE] = CONDIF_HEAT_EDGE_NORMAL;

  strcpy(name[CONDIF_HEAT_EDGE_NORMAL_ELEMENT_NODE_FACTOR],"condif_heat_edge_normal_element_node_factor");
  type[CONDIF_HEAT_EDGE_NORMAL_ELEMENT_NODE_FACTOR] = DOUBLE_PRECISION;
  data_length[CONDIF_HEAT_EDGE_NORMAL_ELEMENT_NODE_FACTOR] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_HEAT_EDGE_NORMAL_ELEMENT_NODE_FACTOR] = 0;
  data_class[CONDIF_HEAT_EDGE_NORMAL_ELEMENT_NODE_FACTOR] = CONDIF;
  data_required[CONDIF_HEAT_EDGE_NORMAL_ELEMENT_NODE_FACTOR] = CONDIF_HEAT_EDGE_NORMAL;

  strcpy(name[CONDIF_HEAT_EDGE_NORMAL_ELEMENT_SIDE],"condif_heat_edge_normal_element_side");
  type[CONDIF_HEAT_EDGE_NORMAL_ELEMENT_SIDE] = INTEGER;
  data_length[CONDIF_HEAT_EDGE_NORMAL_ELEMENT_SIDE] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_HEAT_EDGE_NORMAL_ELEMENT_SIDE] = 0;
  data_class[CONDIF_HEAT_EDGE_NORMAL_ELEMENT_SIDE] = CONDIF;
  data_required[CONDIF_HEAT_EDGE_NORMAL_ELEMENT_SIDE] = CONDIF_HEAT_EDGE_NORMAL;

  strcpy(name[CONDIF_HEAT_EDGE_NORMAL_FACTOR],"condif_heat_edge_normal_factor");
  type[CONDIF_HEAT_EDGE_NORMAL_FACTOR] = DOUBLE_PRECISION;
  data_length[CONDIF_HEAT_EDGE_NORMAL_FACTOR] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_HEAT_EDGE_NORMAL_FACTOR] = 0;
  data_class[CONDIF_HEAT_EDGE_NORMAL_FACTOR] = CONDIF;
  data_required[CONDIF_HEAT_EDGE_NORMAL_FACTOR] = CONDIF_HEAT_EDGE_NORMAL;

  strcpy(name[CONDIF_HEAT_EDGE_NORMAL_GEOMETRY],"condif_heat_edge_normal_geometry");
  type[CONDIF_HEAT_EDGE_NORMAL_GEOMETRY] = INTEGER;
  data_length[CONDIF_HEAT_EDGE_NORMAL_GEOMETRY] = 2;
  data_class[CONDIF_HEAT_EDGE_NORMAL_GEOMETRY] = CONDIF;
  data_required[CONDIF_HEAT_EDGE_NORMAL_GEOMETRY] = CONDIF_HEAT_EDGE_NORMAL;

  strcpy(name[CONDIF_HEAT_EDGE_NORMAL_NODE],"condif_heat_edge_normal_node");
  type[CONDIF_HEAT_EDGE_NORMAL_NODE] = INTEGER;
  data_length[CONDIF_HEAT_EDGE_NORMAL_NODE] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_HEAT_EDGE_NORMAL_NODE] = 0;
  data_class[CONDIF_HEAT_EDGE_NORMAL_NODE] = CONDIF;
  data_required[CONDIF_HEAT_EDGE_NORMAL_NODE] = CONDIF_HEAT_EDGE_NORMAL;

  strcpy(name[CONDIF_HEAT_EDGE_NORMAL_SINE],"condif_heat_edge_normal_sine");
  type[CONDIF_HEAT_EDGE_NORMAL_SINE] = DOUBLE_PRECISION;
  data_length[CONDIF_HEAT_EDGE_NORMAL_SINE] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_HEAT_EDGE_NORMAL_SINE] = 0;
  data_class[CONDIF_HEAT_EDGE_NORMAL_SINE] = CONDIF;
  data_required[CONDIF_HEAT_EDGE_NORMAL_SINE] = CONDIF_HEAT_EDGE_NORMAL;

  strcpy(name[CONDIF_HEAT_EDGE_NORMAL_TIME],"condif_heat_edge_normal_time");
  type[CONDIF_HEAT_EDGE_NORMAL_TIME] = DOUBLE_PRECISION;
  data_length[CONDIF_HEAT_EDGE_NORMAL_TIME] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_HEAT_EDGE_NORMAL_TIME] = 0;
  data_class[CONDIF_HEAT_EDGE_NORMAL_TIME] = CONDIF;
  data_required[CONDIF_HEAT_EDGE_NORMAL_TIME] = CONDIF_HEAT_EDGE_NORMAL;

  strcpy(name[CONDIF_HEAT_VOLUME],"condif_heat_volume");
  type[CONDIF_HEAT_VOLUME] = DOUBLE_PRECISION;
  data_length[CONDIF_HEAT_VOLUME] = 1;
  data_class[CONDIF_HEAT_VOLUME] = CONDIF;

  strcpy(name[CONDIF_HEAT_VOLUME_ELEMENT],"condif_heat_volume_element");
  type[CONDIF_HEAT_VOLUME_ELEMENT] = INTEGER;
  data_length[CONDIF_HEAT_VOLUME_ELEMENT] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_HEAT_VOLUME_ELEMENT] = 0;
  data_class[CONDIF_HEAT_VOLUME_ELEMENT] = CONDIF;
  data_required[CONDIF_HEAT_VOLUME_ELEMENT] = CONDIF_HEAT_VOLUME;

  strcpy(name[CONDIF_HEAT_VOLUME_ELEMENT_GROUP],"condif_heat_volume_element_group");
  type[CONDIF_HEAT_VOLUME_ELEMENT_GROUP] = INTEGER;
  data_length[CONDIF_HEAT_VOLUME_ELEMENT_GROUP] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_HEAT_VOLUME_ELEMENT_GROUP] = 0;
  data_class[CONDIF_HEAT_VOLUME_ELEMENT_GROUP] = CONDIF;
  data_required[CONDIF_HEAT_VOLUME_ELEMENT_GROUP] = CONDIF_HEAT_VOLUME;

  strcpy(name[CONDIF_HEAT_VOLUME_FACTOR],"condif_heat_volume_factor");
  type[CONDIF_HEAT_VOLUME_FACTOR] = DOUBLE_PRECISION;
  data_length[CONDIF_HEAT_VOLUME_FACTOR] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_HEAT_VOLUME_FACTOR] = 0;
  data_class[CONDIF_HEAT_VOLUME_FACTOR] = CONDIF;
  data_required[CONDIF_HEAT_VOLUME_FACTOR] = CONDIF_HEAT_VOLUME;

  strcpy(name[CONDIF_HEAT_VOLUME_GEOMETRY],"condif_heat_volume_geometry");
  type[CONDIF_HEAT_VOLUME_GEOMETRY] = INTEGER;
  data_length[CONDIF_HEAT_VOLUME_GEOMETRY] = 2;
  data_class[CONDIF_HEAT_VOLUME_GEOMETRY] = CONDIF;
  data_required[CONDIF_HEAT_VOLUME_GEOMETRY] = CONDIF_HEAT_VOLUME;

  strcpy(name[CONDIF_HEAT_VOLUME_SINE],"condif_heat_volume_sine");
  type[CONDIF_HEAT_VOLUME_SINE] = DOUBLE_PRECISION;
  data_length[CONDIF_HEAT_VOLUME_SINE] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_HEAT_VOLUME_SINE] = 0;
  data_class[CONDIF_HEAT_VOLUME_SINE] = CONDIF;
  data_required[CONDIF_HEAT_VOLUME_SINE] = CONDIF_HEAT_VOLUME;

  strcpy(name[CONDIF_HEAT_VOLUME_TIME],"condif_heat_volume_time");
  type[CONDIF_HEAT_VOLUME_TIME] = DOUBLE_PRECISION;
  data_length[CONDIF_HEAT_VOLUME_TIME] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_HEAT_VOLUME_TIME] = 0;
  data_class[CONDIF_HEAT_VOLUME_TIME] = CONDIF;
  data_required[CONDIF_HEAT_VOLUME_TIME] = CONDIF_HEAT_VOLUME;

  strcpy(name[CONDIF_HEAT_VOLUME_USER],"condif_heat_volume_user");
  type[CONDIF_HEAT_VOLUME_USER] = INTEGER;
  data_length[CONDIF_HEAT_VOLUME_USER] = 1;
  data_class[CONDIF_HEAT_VOLUME_USER] = CONDIF;
  data_required[CONDIF_HEAT_VOLUME_USER] = CONDIF_HEAT_VOLUME;

  strcpy(name[CONDIF_HEAT_VOLUME_USER_PARAMETERS],"condif_heat_volume_user_parameters");
  type[CONDIF_HEAT_VOLUME_USER_PARAMETERS] = DOUBLE_PRECISION;
  data_length[CONDIF_HEAT_VOLUME_USER_PARAMETERS] = DATA_ITEM_SIZE;
  fixed_length[CONDIF_HEAT_VOLUME_USER_PARAMETERS] = 0;
  data_class[CONDIF_HEAT_VOLUME_USER_PARAMETERS] = CONDIF;
  data_required[CONDIF_HEAT_VOLUME_USER_PARAMETERS] = CONDIF_HEAT_VOLUME;

  strcpy(name[CONDIF_TEMPERATURE],"condif_temperature");

  strcpy(name[COSINUS],"cosinus");

  strcpy(name[CONTACT],"contact");

  strcpy(name[CONTACTSPRING],"contactspring");

  strcpy(name[CONTACT_APPLY],"contact_apply");
  type[CONTACT_APPLY] = INTEGER;
  data_length[CONTACT_APPLY] = 1;
  no_index[CONTACT_APPLY] = 1;
  data_class[CONTACT_APPLY] = CONTACT;

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

  strcpy(name[CONTACT_HEAT_GENERATION],"contact_heat_generation");
  type[CONTACT_HEAT_GENERATION] = DOUBLE_PRECISION;
  data_length[CONTACT_HEAT_GENERATION] = 1;
  no_index[CONTACT_HEAT_GENERATION] = 1;
  data_class[CONTACT_HEAT_GENERATION] = CONTACT;

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

  strcpy(name[CONTACT_PLASTI_FRICTION],"contact_plasti_friction");
  type[CONTACT_PLASTI_FRICTION] = DOUBLE_PRECISION;
  data_length[CONTACT_PLASTI_FRICTION] = 2;
  no_index[CONTACT_PLASTI_FRICTION] = 1;
  data_class[CONTACT_PLASTI_FRICTION] = CONTACT;

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

  strcpy(name[CONTACT_TARGET_ELEMENT_GROUP],"contact_target_element_group");
  type[CONTACT_TARGET_ELEMENT_GROUP] = INTEGER;
  data_length[CONTACT_TARGET_ELEMENT_GROUP] = DATA_ITEM_SIZE;
  fixed_length[CONTACT_TARGET_ELEMENT_GROUP] = 0;
  no_index[CONTACT_TARGET_ELEMENT_GROUP] = 1;
  data_class[CONTACT_TARGET_ELEMENT_GROUP] = CONTACT;

  strcpy(name[CONTACT_TARGET_GEOMETRY],"contact_target_geometry");
  type[CONTACT_TARGET_GEOMETRY] = INTEGER;
  data_length[CONTACT_TARGET_GEOMETRY] = 2;
  data_class[CONTACT_TARGET_GEOMETRY] = CONTACT;

  strcpy(name[CONTACT_TARGET_GEOMETRY_SWITCH],"contact_target_geometry_switch");
  type[CONTACT_TARGET_GEOMETRY_SWITCH] = INTEGER;
  data_length[CONTACT_TARGET_GEOMETRY_SWITCH] = 1;
  data_class[CONTACT_TARGET_GEOMETRY_SWITCH] = CONTACT;
  data_required[CONTACT_TARGET_GEOMETRY_SWITCH] = CONTACT_TARGET_GEOMETRY;

  strcpy(name[CONTROL_CHANGE_DATAITEM_APPLY],"control_change_dataitem_apply");
  type[CONTROL_CHANGE_DATAITEM_APPLY] = INTEGER;
  data_length[CONTROL_CHANGE_DATAITEM_APPLY] = 1;
  data_class[CONTROL_CHANGE_DATAITEM_APPLY] = CONTROL;

  strcpy(name[CONTROL_BOUNDA_RELAX],"control_bounda_relax");
  type[CONTROL_BOUNDA_RELAX] = INTEGER;
  data_length[CONTROL_BOUNDA_RELAX] = 1;
  data_class[CONTROL_BOUNDA_RELAX] = CONTROL;

  strcpy(name[CONTROL_BOUNDA_RELAX_GEOMETRY],"control_bounda_relax_geometry");
  type[CONTROL_BOUNDA_RELAX_GEOMETRY] = INTEGER;
  data_length[CONTROL_BOUNDA_RELAX_GEOMETRY] = 2;
  data_class[CONTROL_BOUNDA_RELAX_GEOMETRY] = CONTROL;
  data_required[CONTROL_BOUNDA_RELAX_GEOMETRY] = CONTROL_BOUNDA_RELAX;

  // SMALL-FAMILY BATCH (2026-09-05): see the enum comment block.
  // post_calcul operator -k0 (manual Professional 6.901: ratio of the
  // average horizontal stress over the vertical stress) and its output
  // item name -k0_sig (post_calcul_label of the .dbs). Pure name
  // entries (INTEGER, resolved by the input reader like -force).
  strcpy(name[K0],"k0");
  type[K0] = INTEGER;
  data_length[K0] = 1;

  strcpy(name[K0_SIG],"k0_sig");
  type[K0_SIG] = INTEGER;
  data_length[K0_SIG] = 1;

  // element_dof_initial (manual Professional 6.422): the dofs the
  // element assumes it had in the past when it comes the first time to
  // live; consumed by the inertia terms of the transient integration
  // (general.cc). Per element index, one value per element dof.
  strcpy(name[ELEMENT_DOF_INITIAL],"element_dof_initial");
  type[ELEMENT_DOF_INITIAL] = DOUBLE_PRECISION;
  data_length[ELEMENT_DOF_INITIAL] = DATA_ITEM_SIZE;
  fixed_length[ELEMENT_DOF_INITIAL] = 0;
  data_class[ELEMENT_DOF_INITIAL] = ELEMENT;

  // internal marker: set at the element's birth step so the initial
  // field is only used once (versioned like the element state records).
  strcpy(name[ELEMENT_DOF_INITIAL_APPLIED],"element_dof_initial_applied");
  type[ELEMENT_DOF_INITIAL_APPLIED] = INTEGER;
  data_length[ELEMENT_DOF_INITIAL_APPLIED] = 1;
  fixed_length[ELEMENT_DOF_INITIAL_APPLIED] = 1;
  version_all[ELEMENT_DOF_INITIAL_APPLIED] = 1;
  external[ELEMENT_DOF_INITIAL_APPLIED] = 0;
  data_class[ELEMENT_DOF_INITIAL_APPLIED] = ELEMENT;

  // post_apply (manual Professional 6.900): global switch of the post
  // processing commands (post_* records evaluated per step). Default
  // -yes; only the post_node_rhside_ratio is exempt. Registered with
  // the -no consumption pending (no corpus test uses -no).
  strcpy(name[POST_APPLY],"post_apply");
  type[POST_APPLY] = INTEGER;
  data_length[POST_APPLY] = 1;
  no_index[POST_APPLY] = 1;

  // print_database_calculation (6.972) / print_gid_calculation: global
  // switches of the final .dbs/.flavia output at the end of the run.
  // Consumed in exit_tn (miscel.cc): -no skips the final database/gid
  // dump (large1 uses it to keep the huge 3D run lean).
  strcpy(name[PRINT_DATABASE_CALCULATION],"print_database_calculation");
  type[PRINT_DATABASE_CALCULATION] = INTEGER;
  data_length[PRINT_DATABASE_CALCULATION] = 1;
  no_index[PRINT_DATABASE_CALCULATION] = 1;

  strcpy(name[PRINT_GID_CALCULATION],"print_gid_calculation");
  type[PRINT_GID_CALCULATION] = INTEGER;
  data_length[PRINT_GID_CALCULATION] = 1;
  no_index[PRINT_GID_CALCULATION] = 1;

  // print_group_data (manual Professional 6.990): plot group_* data
  // items in the GiD output for isoparametric elements (and fill the
  // element_print_group_data records). Registered parse-only in this
  // batch (the GiD group-data writing is pending; distri3 only lists
  // the young modulus distribution without a target on it).
  strcpy(name[PRINT_GROUP_DATA],"print_group_data");
  type[PRINT_GROUP_DATA] = INTEGER;
  data_length[PRINT_GROUP_DATA] = 1;
  fixed_length[PRINT_GROUP_DATA] = 1;
  no_index[PRINT_GROUP_DATA] = 1;

  // geometry_node_type (manual Professional 6.540) / geometry_
  // projection_type (6.543): per-geometry records (same index as the
  // geometry entity) overriding the coordinates used to check nodes on
  // the geometry (-node / -node_start_refined / -plus_displacement,
  // default -node_start_refined) and the projection semantics
  // (-project_inside / -project_exact, default -project_exact).
  // Consumed in geometry() (geometry.cc).
  strcpy(name[GEOMETRY_NODE_TYPE],"geometry_node_type");
  type[GEOMETRY_NODE_TYPE] = INTEGER;
  data_length[GEOMETRY_NODE_TYPE] = 1;
  fixed_length[GEOMETRY_NODE_TYPE] = 1;
  data_class[GEOMETRY_NODE_TYPE] = GEOMETRY;

  strcpy(name[GEOMETRY_PROJECTION_TYPE],"geometry_projection_type");
  type[GEOMETRY_PROJECTION_TYPE] = INTEGER;
  data_length[GEOMETRY_PROJECTION_TYPE] = 1;
  fixed_length[GEOMETRY_PROJECTION_TYPE] = 1;
  data_class[GEOMETRY_PROJECTION_TYPE] = GEOMETRY;

  // -project_inside keyword value (6.543): everything inside the
  // geometry is used (the "filled" semantics of the delete/cut family)
  strcpy(name[PROJECT_INSIDE],"project_inside");
  type[PROJECT_INSIDE] = INTEGER;
  data_length[PROJECT_INSIDE] = 1;

  // print_node_geometry_present (6.995) + _node_type (6.996): switch on
  // the filling of the per-node node_geometry_present record (6.886:
  // the list of geometries in which each node is present) and the
  // default node_type of that check. Filled in step_start (top.cc).
  strcpy(name[PRINT_NODE_GEOMETRY_PRESENT],"print_node_geometry_present");
  type[PRINT_NODE_GEOMETRY_PRESENT] = INTEGER;
  data_length[PRINT_NODE_GEOMETRY_PRESENT] = 1;
  no_index[PRINT_NODE_GEOMETRY_PRESENT] = 1;

  strcpy(name[PRINT_NODE_GEOMETRY_PRESENT_NODE_TYPE],
    "print_node_geometry_present_node_type");
  type[PRINT_NODE_GEOMETRY_PRESENT_NODE_TYPE] = INTEGER;
  data_length[PRINT_NODE_GEOMETRY_PRESENT_NODE_TYPE] = 1;
  no_index[PRINT_NODE_GEOMETRY_PRESENT_NODE_TYPE] = 1;

  // node_geometry_present (6.886): per node index, the list of
  // geometries in which the node is present, stored as pairs
  // (geometry name value, geometry index). Fill on/off with
  // print_node_geometry_present.
  strcpy(name[NODE_GEOMETRY_PRESENT],"node_geometry_present");
  type[NODE_GEOMETRY_PRESENT] = INTEGER;
  data_length[NODE_GEOMETRY_PRESENT] = DATA_ITEM_SIZE;
  fixed_length[NODE_GEOMETRY_PRESENT] = 0;
  version_all[NODE_GEOMETRY_PRESENT] = 1;
  data_class[NODE_GEOMETRY_PRESENT] = NODE;

  strcpy(name[CONTROL_CONTACT_APPLY],"control_contact_apply");
  type[CONTROL_CONTACT_APPLY] = INTEGER;
  data_length[CONTROL_CONTACT_APPLY] = 1;
  data_class[CONTROL_CONTACT_APPLY] = CONTACT;

  strcpy(name[CONTROL_CRACK],"control_crack");
  type[CONTROL_CRACK] = INTEGER;
  data_length[CONTROL_CRACK] = 1;
  data_class[CONTROL_CRACK] = CONTROL;     

  strcpy(name[CONTROL_DATA_ACTIVATE],"control_data_activate");
  type[CONTROL_DATA_ACTIVATE] = INTEGER;
  data_length[CONTROL_DATA_ACTIVATE] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_DATA_ACTIVATE] = 0;
  data_class[CONTROL_DATA_ACTIVATE] = CONTROL;

  strcpy(name[CONTROL_DATA_ARITHMETIC],"control_data_arithmetic");
  type[CONTROL_DATA_ARITHMETIC] = INTEGER;
  data_length[CONTROL_DATA_ARITHMETIC] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_DATA_ARITHMETIC] = 0;
  data_class[CONTROL_DATA_ARITHMETIC] = CONTROL;

  strcpy(name[CONTROL_DATA_ARITHMETIC_DOUBLE],"control_data_arithmetic_double");
  type[CONTROL_DATA_ARITHMETIC_DOUBLE] = DOUBLE_PRECISION;
  data_length[CONTROL_DATA_ARITHMETIC_DOUBLE] = 1;
  data_class[CONTROL_DATA_ARITHMETIC_DOUBLE] = CONTROL;
  data_required[CONTROL_DATA_ARITHMETIC_DOUBLE] = CONTROL_DATA_ARITHMETIC;

  strcpy(name[CONTROL_DATA_COPY],"control_data_copy");
  type[CONTROL_DATA_COPY] = INTEGER;
  data_length[CONTROL_DATA_COPY] = 2;
  data_class[CONTROL_DATA_COPY] = CONTROL;

  strcpy(name[CONTROL_DATA_COPY_FACTOR],"control_data_copy_factor");
  type[CONTROL_DATA_COPY_FACTOR] = DOUBLE_PRECISION;
  data_length[CONTROL_DATA_COPY_FACTOR] = 1;
  data_class[CONTROL_DATA_COPY_FACTOR] = CONTROL;
  data_required[CONTROL_DATA_COPY_FACTOR] = CONTROL_DATA_COPY;

  strcpy(name[CONTROL_DATA_COPY_INDEX],"control_data_copy_index");
  type[CONTROL_DATA_COPY_INDEX] = INTEGER;
  data_length[CONTROL_DATA_COPY_INDEX] = 4;
  data_class[CONTROL_DATA_COPY_INDEX] = CONTROL;

  strcpy(name[CONTROL_DATA_COPY_INDEX_FACTOR],"control_data_copy_index_factor");
  type[CONTROL_DATA_COPY_INDEX_FACTOR] = DOUBLE_PRECISION;
  data_length[CONTROL_DATA_COPY_INDEX_FACTOR] = 1;
  data_class[CONTROL_DATA_COPY_INDEX_FACTOR] = CONTROL;
  data_required[CONTROL_DATA_COPY_INDEX_FACTOR] = CONTROL_DATA_COPY_INDEX;

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

  strcpy(name[CONTROL_DISTRIBUTE_CORRELATION_DISTANCE],"control_distribute_correlation_distance");
  type[CONTROL_DISTRIBUTE_CORRELATION_DISTANCE] = DOUBLE_PRECISION;
  data_length[CONTROL_DISTRIBUTE_CORRELATION_DISTANCE] = 1;
  data_class[CONTROL_DISTRIBUTE_CORRELATION_DISTANCE] = CONTROL;
  data_required[CONTROL_DISTRIBUTE_CORRELATION_DISTANCE] = CONTROL_DISTRIBUTE;

  strcpy(name[CONTROL_DISTRIBUTE_CORRELATION_LENGTH],"control_distribute_correlation_length");
  type[CONTROL_DISTRIBUTE_CORRELATION_LENGTH] = DOUBLE_PRECISION;
  data_length[CONTROL_DISTRIBUTE_CORRELATION_LENGTH] = MDIM;
  fixed_length[CONTROL_DISTRIBUTE_CORRELATION_LENGTH] = 0;
  data_class[CONTROL_DISTRIBUTE_CORRELATION_LENGTH] = CONTROL;
  data_required[CONTROL_DISTRIBUTE_CORRELATION_LENGTH] = CONTROL_DISTRIBUTE;

  strcpy(name[CONTROL_DISTRIBUTE_MINIMUM_MAXIMUM],"control_distribute_minimum_maximum");
  type[CONTROL_DISTRIBUTE_MINIMUM_MAXIMUM] = DOUBLE_PRECISION;
  data_length[CONTROL_DISTRIBUTE_MINIMUM_MAXIMUM] = 2;
  data_class[CONTROL_DISTRIBUTE_MINIMUM_MAXIMUM] = CONTROL;
  data_required[CONTROL_DISTRIBUTE_MINIMUM_MAXIMUM] = CONTROL_DISTRIBUTE;

  strcpy(name[CONTROL_DISTRIBUTE_PARAMETERS],"control_distribute_parameters");
  type[CONTROL_DISTRIBUTE_PARAMETERS] = DOUBLE_PRECISION;
  data_length[CONTROL_DISTRIBUTE_PARAMETERS] = 2;
  data_class[CONTROL_DISTRIBUTE_PARAMETERS] = CONTROL;
  data_required[CONTROL_DISTRIBUTE_PARAMETERS] = CONTROL_DISTRIBUTE;

  strcpy(name[CONTROL_DISTRIBUTE_SEED],"control_distribute_seed");
  type[CONTROL_DISTRIBUTE_SEED] = DOUBLE_PRECISION;
  data_length[CONTROL_DISTRIBUTE_SEED] = 1;
  data_class[CONTROL_DISTRIBUTE_SEED] = CONTROL;
  data_required[CONTROL_DISTRIBUTE_SEED] = CONTROL_DISTRIBUTE;

  strcpy(name[CONTROL_DISTRIBUTE_VALUES],"control_distribute_values");
  type[CONTROL_DISTRIBUTE_VALUES] = DOUBLE_PRECISION;
  data_length[CONTROL_DISTRIBUTE_VALUES] = DATA_ITEM_SIZE;
  data_class[CONTROL_DISTRIBUTE_VALUES] = CONTROL;
  fixed_length[CONTROL_DISTRIBUTE_VALUES] = 0;
  data_required[CONTROL_DISTRIBUTE_VALUES] = CONTROL_DISTRIBUTE;

  strcpy(name[CONTROL_GROUNDFLOW_CONSOLIDATION_APPLY],"control_groundflow_consolidation_apply");
  type[CONTROL_GROUNDFLOW_CONSOLIDATION_APPLY] = INTEGER;
  data_length[CONTROL_GROUNDFLOW_CONSOLIDATION_APPLY] = 1;
  data_class[CONTROL_GROUNDFLOW_CONSOLIDATION_APPLY] = CONTROL;

  strcpy(name[CONTROL_GROUNDFLOW_NONSATURATED_APPLY],"control_groundflow_nonsaturated_apply");
  type[CONTROL_GROUNDFLOW_NONSATURATED_APPLY] = INTEGER;
  data_length[CONTROL_GROUNDFLOW_NONSATURATED_APPLY] = 1;
  data_class[CONTROL_GROUNDFLOW_NONSATURATED_APPLY] = CONTROL;

  // control_groundflow_seepage_apply (manual Professional 6.105): per
  // control step switch of the seepage faces (groundflow_seepage_*).
  // The seepage machinery lives in bounda.cc; the control record gates
  // it per timestep (consumption in bounda.cc seepage branch).
  strcpy(name[CONTROL_GROUNDFLOW_SEEPAGE_APPLY],"control_groundflow_seepage_apply");
  type[CONTROL_GROUNDFLOW_SEEPAGE_APPLY] = INTEGER;
  data_length[CONTROL_GROUNDFLOW_SEEPAGE_APPLY] = 1;
  data_class[CONTROL_GROUNDFLOW_SEEPAGE_APPLY] = CONTROL;

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

  strcpy(name[CONTROL_MESH_ACTIVATE_GRAVITY_APPLY],"control_mesh_activate_gravity_apply");
  type[CONTROL_MESH_ACTIVATE_GRAVITY_APPLY] = INTEGER;
  data_length[CONTROL_MESH_ACTIVATE_GRAVITY_APPLY] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_MESH_ACTIVATE_GRAVITY_APPLY] = 0;
  data_class[CONTROL_MESH_ACTIVATE_GRAVITY_APPLY] = CONTROL;

  strcpy(name[CONTROL_MESH_ADJUST_GEOMETRY],"control_mesh_adjust_geometry");
  type[CONTROL_MESH_ADJUST_GEOMETRY] = INTEGER;
  data_length[CONTROL_MESH_ADJUST_GEOMETRY] = 4;
  data_class[CONTROL_MESH_ADJUST_GEOMETRY] = CONTROL;

  strcpy(name[CONTROL_MESH_CHANGE_ELEMENT_GROUP],"control_mesh_change_element_group");
  type[CONTROL_MESH_CHANGE_ELEMENT_GROUP] = INTEGER;
  data_length[CONTROL_MESH_CHANGE_ELEMENT_GROUP] = 2;
  data_class[CONTROL_MESH_CHANGE_ELEMENT_GROUP] = CONTROL;

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

  strcpy(name[CONTROL_MESH_DELETE_ELEMENT],"control_mesh_delete_element");
  type[CONTROL_MESH_DELETE_ELEMENT] = INTEGER;
  data_length[CONTROL_MESH_DELETE_ELEMENT] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_MESH_DELETE_ELEMENT] = 0;
  data_class[CONTROL_MESH_DELETE_ELEMENT] = CONTROL;

  strcpy(name[CONTROL_MESH_DELETE_SMALL],"control_mesh_delete_small");
  type[CONTROL_MESH_DELETE_SMALL] = DOUBLE_PRECISION;
  data_length[CONTROL_MESH_DELETE_SMALL] = 1;
  data_class[CONTROL_MESH_DELETE_SMALL] = CONTROL;

  strcpy(name[CONTROL_MESH_COPY],"control_mesh_copy");
  type[CONTROL_MESH_COPY] = DOUBLE_PRECISION;
  data_length[CONTROL_MESH_COPY] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_MESH_COPY] = 0;
  data_class[CONTROL_MESH_COPY] = CONTROL;

  strcpy(name[CONTROL_MESH_CUT_GEOMETRY],"control_mesh_cut_geometry");
  type[CONTROL_MESH_CUT_GEOMETRY] = INTEGER;
  data_length[CONTROL_MESH_CUT_GEOMETRY] = 2;
  data_class[CONTROL_MESH_CUT_GEOMETRY] = CONTROL;

  strcpy(name[CONTROL_MESH_CUT_NODE_FORCE],"control_mesh_cut_node_force");
  type[CONTROL_MESH_CUT_NODE_FORCE] = INTEGER;
  data_length[CONTROL_MESH_CUT_NODE_FORCE] = MDIM;
  fixed_length[CONTROL_MESH_CUT_NODE_FORCE] = 0;
  data_class[CONTROL_MESH_CUT_NODE_FORCE] = CONTROL;
  data_required[CONTROL_MESH_CUT_NODE_FORCE] = CONTROL_MESH_CUT_GEOMETRY;

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
  strcpy(name[CONTROL_MESH_GENERATE_CONTACTSPRING_ELEMENT_GROUP],"control_mesh_generate_contactspring_element_group");
  type[CONTROL_MESH_GENERATE_CONTACTSPRING_ELEMENT_GROUP] = INTEGER;
  data_length[CONTROL_MESH_GENERATE_CONTACTSPRING_ELEMENT_GROUP] = 2;
  data_class[CONTROL_MESH_GENERATE_CONTACTSPRING_ELEMENT_GROUP] = CONTROL;
  type[CONTROL_MESH_GENERATE_CONTACTSPRING_ELEMENT] = INTEGER;
  data_length[CONTROL_MESH_GENERATE_CONTACTSPRING_ELEMENT] = 2;
  data_class[CONTROL_MESH_GENERATE_CONTACTSPRING_ELEMENT] = CONTROL;

  strcpy(name[CONTROL_MESH_GENERATE_INTERFACE],"control_mesh_generate_interface");
  type[CONTROL_MESH_GENERATE_INTERFACE] = INTEGER;
  data_length[CONTROL_MESH_GENERATE_INTERFACE] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_MESH_GENERATE_INTERFACE] = 0;
  data_class[CONTROL_MESH_GENERATE_INTERFACE] = CONTROL;

  strcpy(name[CONTROL_MESH_GENERATE_INTERFACE_GEOMETRY],"control_mesh_generate_interface_geometry");
  type[CONTROL_MESH_GENERATE_INTERFACE_GEOMETRY] = INTEGER;
  data_length[CONTROL_MESH_GENERATE_INTERFACE_GEOMETRY] = 2;
  data_class[CONTROL_MESH_GENERATE_INTERFACE_GEOMETRY] = CONTROL;
  data_required[CONTROL_MESH_GENERATE_INTERFACE_GEOMETRY] = CONTROL_MESH_GENERATE_INTERFACE;

  strcpy(name[CONTROL_MESH_GENERATE_INTERFACE_METHOD],"control_mesh_generate_interface_method");
  type[CONTROL_MESH_GENERATE_INTERFACE_METHOD] = INTEGER;
  data_length[CONTROL_MESH_GENERATE_INTERFACE_METHOD] = 2;
  data_class[CONTROL_MESH_GENERATE_INTERFACE_METHOD] = CONTROL;
  data_required[CONTROL_MESH_GENERATE_INTERFACE_METHOD] = CONTROL_MESH_GENERATE_INTERFACE;

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

  strcpy(name[CONTROL_MESH_CONVERT],"control_mesh_convert");
  type[CONTROL_MESH_CONVERT] = INTEGER;
  data_length[CONTROL_MESH_CONVERT] = 1;
  data_class[CONTROL_MESH_CONVERT] = CONTROL;

  strcpy(name[CONTROL_MESH_CONVERT_ELEMENT_GROUP],"control_mesh_convert_element_group");
  type[CONTROL_MESH_CONVERT_ELEMENT_GROUP] = INTEGER;
  data_length[CONTROL_MESH_CONVERT_ELEMENT_GROUP] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_MESH_CONVERT_ELEMENT_GROUP] = 0;
  data_class[CONTROL_MESH_CONVERT_ELEMENT_GROUP] = CONTROL;
  data_required[CONTROL_MESH_CONVERT_ELEMENT_GROUP] = CONTROL_MESH_CONVERT;

  strcpy(name[CONTROL_MESH_KEEP_ELEMENT],"control_mesh_keep_element");
  type[CONTROL_MESH_KEEP_ELEMENT] = INTEGER;
  data_length[CONTROL_MESH_KEEP_ELEMENT] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_MESH_KEEP_ELEMENT] = 0;
  data_class[CONTROL_MESH_KEEP_ELEMENT] = CONTROL;

  strcpy(name[CONTROL_MESH_KEEP_ELEMENT_GROUP],"control_mesh_keep_element_group");
  type[CONTROL_MESH_KEEP_ELEMENT_GROUP] = INTEGER;
  data_length[CONTROL_MESH_KEEP_ELEMENT_GROUP] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_MESH_KEEP_ELEMENT_GROUP] = 0;
  data_class[CONTROL_MESH_KEEP_ELEMENT_GROUP] = CONTROL;

  strcpy(name[CONTROL_MESH_KEEP_GEOMETRY],"control_mesh_keep_geometry");
  type[CONTROL_MESH_KEEP_GEOMETRY] = INTEGER;
  data_length[CONTROL_MESH_KEEP_GEOMETRY] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_MESH_KEEP_GEOMETRY] = 0;
  data_class[CONTROL_MESH_KEEP_GEOMETRY] = CONTROL;

  strcpy(name[CONTROL_MESH_KEEP_NODE],"control_mesh_keep_node");
  type[CONTROL_MESH_KEEP_NODE] = INTEGER;
  data_length[CONTROL_MESH_KEEP_NODE] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_MESH_KEEP_NODE] = 0;
  data_class[CONTROL_MESH_KEEP_NODE] = CONTROL;

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

  strcpy(name[CONTROL_MESH_MERGE_NOT],"control_mesh_merge_geometry_not");
  type[CONTROL_MESH_MERGE_NOT] = INTEGER;
  data_length[CONTROL_MESH_MERGE_NOT] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_MESH_MERGE_NOT] = 0;
  data_class[CONTROL_MESH_MERGE_NOT] = CONTROL;

  strcpy(name[CONTROL_MESH_MERGE_GEOMETRY],"control_mesh_merge_geometry");
  type[CONTROL_MESH_MERGE_GEOMETRY] = INTEGER;
  data_length[CONTROL_MESH_MERGE_GEOMETRY] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_MESH_MERGE_GEOMETRY] = 0;
  data_class[CONTROL_MESH_MERGE_GEOMETRY] = CONTROL;

  strcpy(name[CONTROL_MESH_MIRROR],"control_mesh_mirror");
  type[CONTROL_MESH_MIRROR] = INTEGER;
  data_length[CONTROL_MESH_MIRROR] = 1;
  data_class[CONTROL_MESH_MIRROR] = CONTROL;

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

  strcpy(name[CONTROL_MESH_MOVE],"control_mesh_move");
  type[CONTROL_MESH_MOVE] = DOUBLE_PRECISION;
  data_length[CONTROL_MESH_MOVE] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_MESH_MOVE] = 0;
  data_class[CONTROL_MESH_MOVE] = CONTROL;

  strcpy(name[CONTROL_MESH_MULTIPLY],"control_mesh_multiply");
  type[CONTROL_MESH_MULTIPLY] = INTEGER;
  data_length[CONTROL_MESH_MULTIPLY] = 1;
  data_class[CONTROL_MESH_MULTIPLY] = CONTROL;

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

  strcpy(name[CONTROL_MESH_REMOVE],"control_mesh_remove");
  type[CONTROL_MESH_REMOVE] = INTEGER;
  data_length[CONTROL_MESH_REMOVE] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_MESH_REMOVE] = 0;
  data_class[CONTROL_MESH_REMOVE] = CONTROL;

  strcpy(name[CONTROL_MESH_RENUMBER],"control_mesh_renumber");
  type[CONTROL_MESH_RENUMBER] = INTEGER;
  data_length[CONTROL_MESH_RENUMBER] = 2;
  data_class[CONTROL_MESH_RENUMBER] = CONTROL;

  strcpy(name[CONTROL_MESH_ROTATE],"control_mesh_rotate");
  type[CONTROL_MESH_ROTATE] = INTEGER;
  data_length[CONTROL_MESH_ROTATE] = 1;
  data_class[CONTROL_MESH_ROTATE] = CONTROL;
  data_required[CONTROL_MESH_ROTATE] = CONTROL_MESH_ROTATE;

  strcpy(name[CONTROL_MESH_ROTATE_ANGLE],"control_mesh_rotate_angle");
  type[CONTROL_MESH_ROTATE_ANGLE] = DOUBLE_PRECISION;
  data_length[CONTROL_MESH_ROTATE_ANGLE] = 1;
  data_class[CONTROL_MESH_ROTATE_ANGLE] = CONTROL;

  strcpy(name[CONTROL_MESH_SPLIT],"control_mesh_split");
  type[CONTROL_MESH_SPLIT] = INTEGER;
  data_length[CONTROL_MESH_SPLIT] = 1;
  data_class[CONTROL_MESH_SPLIT] = CONTROL;

  strcpy(name[CONTROL_MESH_SPLIT_ONLY],"control_mesh_split_only");
  type[CONTROL_MESH_SPLIT_ONLY] = INTEGER;
  data_length[CONTROL_MESH_SPLIT_ONLY] = 2;
  data_class[CONTROL_MESH_SPLIT_ONLY] = CONTROL;

  strcpy(name[CONTROL_MESH_SWITCH],"control_mesh_switch");
  type[CONTROL_MESH_SWITCH] = INTEGER;
  data_length[CONTROL_MESH_SWITCH] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_MESH_SWITCH] = 0;
  data_class[CONTROL_MESH_SWITCH] = CONTROL;

  strcpy(name[CONTROL_RESET_DOF],"control_reset_dof");
  type[CONTROL_RESET_DOF] = INTEGER;
  data_length[CONTROL_RESET_DOF] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_RESET_DOF] = 0;
  data_class[CONTROL_RESET_DOF] = CONTROL;

  strcpy(name[CONTROL_RESET_ELEMENT_DOF],"control_reset_element_dof");
  type[CONTROL_RESET_ELEMENT_DOF] = INTEGER;
  data_length[CONTROL_RESET_ELEMENT_DOF] = 1;
  data_class[CONTROL_RESET_ELEMENT_DOF] = CONTROL;
  data_required[CONTROL_RESET_ELEMENT_DOF] = CONTROL_RESET_DOF;

  strcpy(name[CONTROL_RESET_ELEMENT_GROUP],"control_reset_element_group");
  type[CONTROL_RESET_ELEMENT_GROUP] = INTEGER;
  data_length[CONTROL_RESET_ELEMENT_GROUP] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_RESET_ELEMENT_GROUP] = 0;
  data_class[CONTROL_RESET_ELEMENT_GROUP] = CONTROL;
  data_required[CONTROL_RESET_ELEMENT_GROUP] = CONTROL_RESET_DOF;

  strcpy(name[CONTROL_RESET_GEOMETRY],"control_reset_geometry");
  type[CONTROL_RESET_GEOMETRY] = INTEGER;
  data_length[CONTROL_RESET_GEOMETRY] = 2;
  data_class[CONTROL_RESET_GEOMETRY] = CONTROL;
  data_required[CONTROL_RESET_GEOMETRY] = CONTROL_RESET_DOF;

  strcpy(name[CONTROL_RESET_INTERFACE],"control_reset_interface");
  type[CONTROL_RESET_INTERFACE] = INTEGER;
  data_length[CONTROL_RESET_INTERFACE] = 2;
  data_class[CONTROL_RESET_INTERFACE] = CONTROL;

  strcpy(name[CONTROL_RESET_INTERFACE_STRAIN],"control_reset_interface_strain");
  type[CONTROL_RESET_INTERFACE_STRAIN] = INTEGER;
  data_length[CONTROL_RESET_INTERFACE_STRAIN] = 2;
  data_class[CONTROL_RESET_INTERFACE_STRAIN] = CONTROL;

  strcpy(name[CONTROL_RESET_NODE],"control_reset_node");
  type[CONTROL_RESET_NODE] = INTEGER;
  data_length[CONTROL_RESET_NODE] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_RESET_NODE] = 0;
  data_class[CONTROL_RESET_NODE] = CONTROL;
  data_required[CONTROL_RESET_NODE] = CONTROL_RESET_DOF;

  strcpy(name[CONTROL_RESET_VALUE_CONSTANT],"control_reset_value_constant");
  type[CONTROL_RESET_VALUE_CONSTANT] = DOUBLE_PRECISION;
  data_length[CONTROL_RESET_VALUE_CONSTANT] = 1;
  data_class[CONTROL_RESET_VALUE_CONSTANT] = CONTROL;
  data_required[CONTROL_RESET_VALUE_CONSTANT] = CONTROL_RESET_DOF;

  strcpy(name[CONTROL_RESET_VALUE_DOF],"control_reset_value_dof");
  type[CONTROL_RESET_VALUE_DOF] = INTEGER;
  data_length[CONTROL_RESET_VALUE_DOF] = 1;
  data_class[CONTROL_RESET_VALUE_DOF] = CONTROL;
  data_required[CONTROL_RESET_VALUE_DOF] = CONTROL_RESET_DOF;

  strcpy(name[CONTROL_RESET_VALUE_DOF_DIAGRAM],"control_reset_value_dof_diagram");
  type[CONTROL_RESET_VALUE_DOF_DIAGRAM] = DOUBLE_PRECISION;
  data_length[CONTROL_RESET_VALUE_DOF_DIAGRAM] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_RESET_VALUE_DOF_DIAGRAM] = 0;
  data_class[CONTROL_RESET_VALUE_DOF_DIAGRAM] = CONTROL;
  data_required[CONTROL_RESET_VALUE_DOF_DIAGRAM] = CONTROL_RESET_VALUE_DOF;

  strcpy(name[CONTROL_RESET_VALUE_EXPONENT],"control_reset_value_exponent");
  type[CONTROL_RESET_VALUE_EXPONENT] = DOUBLE_PRECISION;
  data_length[CONTROL_RESET_VALUE_EXPONENT] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_RESET_VALUE_EXPONENT] = 0;
  data_class[CONTROL_RESET_VALUE_EXPONENT] = CONTROL;
  data_required[CONTROL_RESET_VALUE_EXPONENT] = CONTROL_RESET_DOF;

  strcpy(name[CONTROL_RESET_VALUE_LINEAR],"control_reset_value_linear");
  type[CONTROL_RESET_VALUE_LINEAR] = DOUBLE_PRECISION;
  data_length[CONTROL_RESET_VALUE_LINEAR] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_RESET_VALUE_LINEAR] = 0;
  data_class[CONTROL_RESET_VALUE_LINEAR] = CONTROL;
  data_required[CONTROL_RESET_VALUE_LINEAR] = CONTROL_RESET_DOF;

  strcpy(name[CONTROL_RESET_VALUE_LOGARITHMIC],"control_reset_value_logarithmic");
  type[CONTROL_RESET_VALUE_LOGARITHMIC] = DOUBLE_PRECISION;
  data_length[CONTROL_RESET_VALUE_LOGARITHMIC] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_RESET_VALUE_LOGARITHMIC] = 0;
  data_class[CONTROL_RESET_VALUE_LOGARITHMIC] = CONTROL;
  data_required[CONTROL_RESET_VALUE_LOGARITHMIC] = CONTROL_RESET_DOF;

  strcpy(name[CONTROL_RESET_VALUE_LOGARITHMIC_SECOND],"control_reset_value_logarithmic_second");
  type[CONTROL_RESET_VALUE_LOGARITHMIC_SECOND] = DOUBLE_PRECISION;
  data_length[CONTROL_RESET_VALUE_LOGARITHMIC_SECOND] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_RESET_VALUE_LOGARITHMIC_SECOND] = 0;
  data_class[CONTROL_RESET_VALUE_LOGARITHMIC_SECOND] = CONTROL;
  data_required[CONTROL_RESET_VALUE_LOGARITHMIC_SECOND] = CONTROL_RESET_DOF;

  strcpy(name[CONTROL_RESET_VALUE_METHOD],"control_reset_value_method");
  type[CONTROL_RESET_VALUE_METHOD] = INTEGER;
  data_length[CONTROL_RESET_VALUE_METHOD] = 1;
  data_class[CONTROL_RESET_VALUE_METHOD] = CONTROL;
  data_required[CONTROL_RESET_VALUE_METHOD] = CONTROL_RESET_DOF;

  strcpy(name[CONTROL_RESET_VALUE_MULTI_LINEAR],"control_reset_value_multi_linear");
  type[CONTROL_RESET_VALUE_MULTI_LINEAR] = DOUBLE_PRECISION;
  data_length[CONTROL_RESET_VALUE_MULTI_LINEAR] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_RESET_VALUE_MULTI_LINEAR] = 0;
  data_class[CONTROL_RESET_VALUE_MULTI_LINEAR] = CONTROL;
  data_required[CONTROL_RESET_VALUE_MULTI_LINEAR] = CONTROL_RESET_DOF;

  strcpy(name[CONTROL_RESET_VALUE_POWER],"control_reset_value_power");
  type[CONTROL_RESET_VALUE_POWER] = DOUBLE_PRECISION;
  data_length[CONTROL_RESET_VALUE_POWER] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_RESET_VALUE_POWER] = 0;
  data_class[CONTROL_RESET_VALUE_POWER] = CONTROL;
  data_required[CONTROL_RESET_VALUE_POWER] = CONTROL_RESET_DOF;

  strcpy(name[CONTROL_RESET_VALUE_SQUARE_ROOT],"control_reset_value_square_root");
  type[CONTROL_RESET_VALUE_SQUARE_ROOT] = DOUBLE_PRECISION;
  data_length[CONTROL_RESET_VALUE_SQUARE_ROOT] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_RESET_VALUE_SQUARE_ROOT] = 0;
  data_class[CONTROL_RESET_VALUE_SQUARE_ROOT] = CONTROL;
  data_required[CONTROL_RESET_VALUE_SQUARE_ROOT] = CONTROL_RESET_DOF;

  strcpy(name[CONTROL_OPTIONS_CONVECTION],"control_options_convection");
  type[CONTROL_OPTIONS_CONVECTION] = INTEGER;
  data_length[CONTROL_OPTIONS_CONVECTION] = 1;
  data_class[CONTROL_OPTIONS_CONVECTION] = CONTROL;

  strcpy(name[CONTROL_OPTIONS_INERTIA],"control_options_inertia");
  type[CONTROL_OPTIONS_INERTIA] = INTEGER;
  data_length[CONTROL_OPTIONS_INERTIA] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_OPTIONS_INERTIA] = 0;
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

  strcpy(name[CONTROL_OPTIONS_SKIP_GROUNDFLOW_MATERIDIVERGENCE],"control_options_skip_groundflow_materidivergence");
  type[CONTROL_OPTIONS_SKIP_GROUNDFLOW_MATERIDIVERGENCE] = INTEGER;
  data_length[CONTROL_OPTIONS_SKIP_GROUNDFLOW_MATERIDIVERGENCE] = 1;
  data_class[CONTROL_OPTIONS_SKIP_GROUNDFLOW_MATERIDIVERGENCE] = CONTROL;

  strcpy(name[CONTROL_OPTIONS_SOLVER],"control_options_solver");
  type[CONTROL_OPTIONS_SOLVER] = INTEGER;
  data_length[CONTROL_OPTIONS_SOLVER] = 1;
  data_class[CONTROL_OPTIONS_SOLVER] = CONTROL;

  strcpy(name[CONTROL_SOLVER_BICG_RESTART],"control_solver_bicg_restart");
  type[CONTROL_SOLVER_BICG_RESTART] = INTEGER;
  data_length[CONTROL_SOLVER_BICG_RESTART] = 1;
  data_class[CONTROL_SOLVER_BICG_RESTART] = CONTROL;

  strcpy(name[CONTROL_SOLVER_BICG_STOP],"control_solver_bicg_stop");
  type[CONTROL_SOLVER_BICG_STOP] = INTEGER;
  data_length[CONTROL_SOLVER_BICG_STOP] = 1;
  data_class[CONTROL_SOLVER_BICG_STOP] = CONTROL;

  strcpy(name[CONTROL_SOLVER_MATRIX_SAVE],"control_solver_matrix_save");
  type[CONTROL_SOLVER_MATRIX_SAVE] = INTEGER;
  data_length[CONTROL_SOLVER_MATRIX_SAVE] = 1;
  data_class[CONTROL_SOLVER_MATRIX_SAVE] = CONTROL;

  strcpy(name[CONTROL_SOLVER_PARDISO_ORDERING],"control_solver_pardiso_ordering");
  type[CONTROL_SOLVER_PARDISO_ORDERING] = INTEGER;
  data_length[CONTROL_SOLVER_PARDISO_ORDERING] = 1;
  data_class[CONTROL_SOLVER_PARDISO_ORDERING] = CONTROL;

  strcpy(name[CONTROL_SOLVER_PARDISO_OUT_OF_CORE],"control_solver_pardiso_out_of_core");
  type[CONTROL_SOLVER_PARDISO_OUT_OF_CORE] = INTEGER;
  data_length[CONTROL_SOLVER_PARDISO_OUT_OF_CORE] = 1;
  data_class[CONTROL_SOLVER_PARDISO_OUT_OF_CORE] = CONTROL;

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

  strcpy(name[CONTROL_DATA_SAVE],"control_data_save");
  type[CONTROL_DATA_SAVE] = INTEGER;
  data_length[CONTROL_DATA_SAVE] = 1;
  data_class[CONTROL_DATA_SAVE] = CONTROL;

  strcpy(name[CONTROL_PRINT_GID_SAVE_DIFFERENCE],"control_print_gid_save_difference");
  type[CONTROL_PRINT_GID_SAVE_DIFFERENCE] = INTEGER;
  data_length[CONTROL_PRINT_GID_SAVE_DIFFERENCE] = 1;
  data_class[CONTROL_PRINT_GID_SAVE_DIFFERENCE] = CONTROL;

  strcpy(name[CONTROL_PRINT_DATABASE_METHOD],"control_print_database_method");
  type[CONTROL_PRINT_DATABASE_METHOD] = INTEGER;
  data_length[CONTROL_PRINT_DATABASE_METHOD] = 1;
  data_class[CONTROL_PRINT_DATABASE_METHOD] = CONTROL;

  strcpy(name[CONTROL_PRINT_ELEMENT],"control_print_element");
  type[CONTROL_PRINT_ELEMENT] = INTEGER;
  data_length[CONTROL_PRINT_ELEMENT] = 1;
  data_class[CONTROL_PRINT_ELEMENT] = CONTROL;

  strcpy(name[CONTROL_PRINT_ELEMENT_METHOD],"control_print_element_method");
  type[CONTROL_PRINT_ELEMENT_METHOD] = INTEGER;
  data_length[CONTROL_PRINT_ELEMENT_METHOD] = 1;
  data_class[CONTROL_PRINT_ELEMENT_METHOD] = CONTROL;
  data_required[CONTROL_PRINT_ELEMENT_METHOD] = CONTROL_PRINT_ELEMENT;

  strcpy(name[CONTROL_PRINT_DATA_VERSUS_DATA],"control_print_data_versus_data");
  type[CONTROL_PRINT_DATA_VERSUS_DATA] = INTEGER;
  data_length[CONTROL_PRINT_DATA_VERSUS_DATA] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_PRINT_DATA_VERSUS_DATA] = 0;
  data_class[CONTROL_PRINT_DATA_VERSUS_DATA] = CONTROL;

  strcpy(name[CONTROL_PRINT_DATA_VERSUS_DATA_FACTOR],"control_print_data_versus_data_factor");
  type[CONTROL_PRINT_DATA_VERSUS_DATA_FACTOR] = DOUBLE_PRECISION;
  data_length[CONTROL_PRINT_DATA_VERSUS_DATA_FACTOR] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_PRINT_DATA_VERSUS_DATA_FACTOR] = 0;
  data_class[CONTROL_PRINT_DATA_VERSUS_DATA_FACTOR] = CONTROL;
  data_required[CONTROL_PRINT_DATA_VERSUS_DATA_FACTOR] = CONTROL_PRINT_DATA_VERSUS_DATA;

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

  strcpy(name[CONTROL_PRINT_HISTORY_FACTOR],"control_print_history_factor");
  type[CONTROL_PRINT_HISTORY_FACTOR] = DOUBLE_PRECISION;
  data_length[CONTROL_PRINT_HISTORY_FACTOR] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_PRINT_HISTORY_FACTOR] = 0;
  data_class[CONTROL_PRINT_HISTORY_FACTOR] = CONTROL;
  data_required[CONTROL_PRINT_HISTORY_FACTOR] = CONTROL_PRINT_HISTORY;

  strcpy(name[CONTROL_PRINT_HISTORY_SMOOTH],"control_print_history_smooth");
  type[CONTROL_PRINT_HISTORY_SMOOTH] = INTEGER;
  data_length[CONTROL_PRINT_HISTORY_SMOOTH] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_PRINT_HISTORY_SMOOTH] = 0;
  data_class[CONTROL_PRINT_HISTORY_SMOOTH] = CONTROL;
  data_required[CONTROL_PRINT_HISTORY_SMOOTH] = CONTROL_PRINT_HISTORY;

  strcpy(name[CONTROL_PRINT_HISTORY_RELATIVE_TIME],"control_print_history_relative_time");
  type[CONTROL_PRINT_HISTORY_RELATIVE_TIME] = DOUBLE_PRECISION;
  data_length[CONTROL_PRINT_HISTORY_RELATIVE_TIME] = 1;
  data_class[CONTROL_PRINT_HISTORY_RELATIVE_TIME] = CONTROL;
  data_required[CONTROL_PRINT_HISTORY_RELATIVE_TIME] = CONTROL_PRINT_HISTORY;

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

  strcpy(name[CONTROL_PRINT_INTERFACE_STRESS],"control_print_interface_stress");
  type[CONTROL_PRINT_INTERFACE_STRESS] = INTEGER;
  data_length[CONTROL_PRINT_INTERFACE_STRESS] = 1;
  data_class[CONTROL_PRINT_INTERFACE_STRESS] = CONTROL;

  strcpy(name[CONTROL_PRINT_INTERFACE_STRESS_2D_COORDINATES],"control_print_interface_stress_2d_coordinates");
  type[CONTROL_PRINT_INTERFACE_STRESS_2D_COORDINATES] = DOUBLE_PRECISION;
  data_length[CONTROL_PRINT_INTERFACE_STRESS_2D_COORDINATES] = 4;
  data_class[CONTROL_PRINT_INTERFACE_STRESS_2D_COORDINATES] = CONTROL;
  data_required[CONTROL_PRINT_INTERFACE_STRESS_2D_COORDINATES] = CONTROL_PRINT_INTERFACE_STRESS;

  strcpy(name[CONTROL_PRINT_INTERFACE_STRESS_3D_GEOMETRY],"control_print_interface_stress_3d_geometry");
  type[CONTROL_PRINT_INTERFACE_STRESS_3D_GEOMETRY] = INTEGER;
  data_length[CONTROL_PRINT_INTERFACE_STRESS_3D_GEOMETRY] = 2;
  data_class[CONTROL_PRINT_INTERFACE_STRESS_3D_GEOMETRY] = CONTROL;
  data_required[CONTROL_PRINT_INTERFACE_STRESS_3D_GEOMETRY] = CONTROL_PRINT_INTERFACE_STRESS;

  strcpy(name[CONTROL_PRINT_INTERFACE_STRESS_3D_ORDER],"control_print_interface_stress_3d_order");
  type[CONTROL_PRINT_INTERFACE_STRESS_3D_ORDER] = INTEGER;
  data_length[CONTROL_PRINT_INTERFACE_STRESS_3D_ORDER] = 1;
  data_class[CONTROL_PRINT_INTERFACE_STRESS_3D_ORDER] = CONTROL;
  data_required[CONTROL_PRINT_INTERFACE_STRESS_3D_ORDER] = CONTROL_PRINT_INTERFACE_STRESS;

  strcpy(name[CONTROL_PRINT_MATLAB],"control_print_matlab");
  type[CONTROL_PRINT_MATLAB] = INTEGER;
  data_length[CONTROL_PRINT_MATLAB] = 1;
  data_class[CONTROL_PRINT_MATLAB] = CONTROL;

  strcpy(name[CONTROL_PRINT_NUMBER_ITERATIONS],"control_print_number_iterations");
  type[CONTROL_PRINT_NUMBER_ITERATIONS] = INTEGER;
  data_length[CONTROL_PRINT_NUMBER_ITERATIONS] = 1;
  data_class[CONTROL_PRINT_NUMBER_ITERATIONS] = CONTROL;

  strcpy(name[CONTROL_PRINT_PARTIALNAME],"control_print_partialname");
  type[CONTROL_PRINT_PARTIALNAME] = INTEGER;
  data_length[CONTROL_PRINT_PARTIALNAME] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_PRINT_PARTIALNAME] = 0;
  data_class[CONTROL_PRINT_PARTIALNAME] = CONTROL;

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

  strcpy(name[CONTROL_PRINT_VTK_DOF],"control_print_vtk_dof");
  type[CONTROL_PRINT_VTK_DOF] = INTEGER;
  data_length[CONTROL_PRINT_VTK_DOF] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_PRINT_VTK_DOF] = 0;
  data_class[CONTROL_PRINT_VTK_DOF] = CONTROL;
  data_required[CONTROL_PRINT_VTK_DOF] = CONTROL_PRINT_VTK;

  strcpy(name[CONTROL_PRINT_VTK_COORD],"control_print_vtk_coord");
  type[CONTROL_PRINT_VTK_COORD] = INTEGER;
  data_length[CONTROL_PRINT_VTK_COORD] = 1;
  data_class[CONTROL_PRINT_VTK_COORD] = CONTROL;
  data_required[CONTROL_PRINT_VTK_COORD] = CONTROL_PRINT_VTK;

  strcpy(name[CONTROL_PRINT_VTK_DOF_CALCUL],"control_print_vtk_dof_calcul");
  type[CONTROL_PRINT_VTK_DOF_CALCUL] = INTEGER;
  data_length[CONTROL_PRINT_VTK_DOF_CALCUL] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_PRINT_VTK_DOF_CALCUL] = 0;
  data_class[CONTROL_PRINT_VTK_DOF_CALCUL] = CONTROL;
  data_required[CONTROL_PRINT_VTK_DOF_CALCUL] = CONTROL_PRINT_VTK;

  strcpy(name[CONTROL_PRINT_VTK_EMPTY],"control_print_vtk_empty");
  type[CONTROL_PRINT_VTK_EMPTY] = INTEGER;
  data_length[CONTROL_PRINT_VTK_EMPTY] = 1;
  data_class[CONTROL_PRINT_VTK_EMPTY] = CONTROL;
  data_required[CONTROL_PRINT_VTK_EMPTY] = CONTROL_PRINT_VTK;

  strcpy(name[CONTROL_PRINT_VTK_NODE_METHOD],"control_print_vtk_node_method");
  type[CONTROL_PRINT_VTK_NODE_METHOD] = INTEGER;
  data_length[CONTROL_PRINT_VTK_NODE_METHOD] = 1;
  data_class[CONTROL_PRINT_VTK_NODE_METHOD] = CONTROL;
  data_required[CONTROL_PRINT_VTK_NODE_METHOD] = CONTROL_PRINT_VTK;

  strcpy(name[CONTROL_PRINT_VTK_OTHER],"control_print_vtk_other");
  type[CONTROL_PRINT_VTK_OTHER] = INTEGER;
  data_length[CONTROL_PRINT_VTK_OTHER] = 1;
  data_class[CONTROL_PRINT_VTK_OTHER] = CONTROL;
  data_required[CONTROL_PRINT_VTK_OTHER] = CONTROL_PRINT_VTK;

  strcpy(name[CONTROL_PRINT_TABULAR],"control_print_tabular");
  type[CONTROL_PRINT_TABULAR] = INTEGER;
  data_length[CONTROL_PRINT_TABULAR] = 1;
  fixed_length[CONTROL_PRINT_TABULAR] = 0;
  data_class[CONTROL_PRINT_TABULAR] = CONTROL;

  strcpy(name[CONTROL_PRINT_GMSH],"control_print_gmsh");
  type[CONTROL_PRINT_GMSH] = INTEGER;
  data_length[CONTROL_PRINT_GMSH] = 1;
  fixed_length[CONTROL_PRINT_GMSH] = 0;
  data_class[CONTROL_PRINT_GMSH] = CONTROL;

  strcpy(name[CONTROL_PRINT_GMSH_DUMMY],"control_print_gmsh_dummy");
  type[CONTROL_PRINT_GMSH_DUMMY] = INTEGER;
  data_length[CONTROL_PRINT_GMSH_DUMMY] = 1;
  fixed_length[CONTROL_PRINT_GMSH_DUMMY] = 0;
  data_class[CONTROL_PRINT_GMSH_DUMMY] = CONTROL;

  strcpy(name[CONTROL_PRINT_GMSH_ELEMENT_DATA],"control_print_gmsh_element_data");
  type[CONTROL_PRINT_GMSH_ELEMENT_DATA] = INTEGER;
  data_length[CONTROL_PRINT_GMSH_ELEMENT_DATA] = 1;
  fixed_length[CONTROL_PRINT_GMSH_ELEMENT_DATA] = 0;
  data_class[CONTROL_PRINT_GMSH_ELEMENT_DATA] = CONTROL;

  strcpy(name[CONTROL_PRINT_GMSH_NODE_METHOD],"control_print_gmsh_node_method");
  type[CONTROL_PRINT_GMSH_NODE_METHOD] = INTEGER;
  data_length[CONTROL_PRINT_GMSH_NODE_METHOD] = 1;
  fixed_length[CONTROL_PRINT_GMSH_NODE_METHOD] = 0;
  data_class[CONTROL_PRINT_GMSH_NODE_METHOD] = CONTROL;

  strcpy(name[CONTROL_PRINT_FRD],"control_print_frd");
  type[CONTROL_PRINT_FRD] = INTEGER;
  data_length[CONTROL_PRINT_FRD] = 1;
  fixed_length[CONTROL_PRINT_FRD] = 0;
  data_class[CONTROL_PRINT_FRD] = CONTROL;

  strcpy(name[CONTROL_PRINT_FRD_FREECAD],"control_print_frd_freecad");
  type[CONTROL_PRINT_FRD_FREECAD] = INTEGER;
  data_length[CONTROL_PRINT_FRD_FREECAD] = 1;
  fixed_length[CONTROL_PRINT_FRD_FREECAD] = 0;
  data_class[CONTROL_PRINT_FRD_FREECAD] = CONTROL;

  strcpy(name[CONTROL_PRINT_FRD_PREPOMAX],"control_print_frd_prepomax");
  type[CONTROL_PRINT_FRD_PREPOMAX] = INTEGER;
  data_length[CONTROL_PRINT_FRD_PREPOMAX] = 1;
  fixed_length[CONTROL_PRINT_FRD_PREPOMAX] = 0;
  data_class[CONTROL_PRINT_FRD_PREPOMAX] = CONTROL;

  strcpy(name[CONTROL_PRINT_DOF],"control_print_dof");
  type[CONTROL_PRINT_DOF] = INTEGER;
  data_length[CONTROL_PRINT_DOF] = 1;
  fixed_length[CONTROL_PRINT_DOF] = 0;
  data_class[CONTROL_PRINT_DOF] = CONTROL;

  strcpy(name[CONTROL_PRINT_DOF_ID],"control_print_dof_id");
  type[CONTROL_PRINT_DOF_ID] = INTEGER;
  data_length[CONTROL_PRINT_DOF_ID] = 1;
  data_class[CONTROL_PRINT_DOF_ID] = CONTROL;
  data_required[CONTROL_PRINT_DOF_ID] = CONTROL_PRINT_DOF;

  strcpy(name[CONTROL_PRINT_DOF_LINE],"control_print_dof_line");
  type[CONTROL_PRINT_DOF_LINE] = INTEGER;
  data_length[CONTROL_PRINT_DOF_LINE] = 1;
  fixed_length[CONTROL_PRINT_DOF_LINE] = 0;
  data_class[CONTROL_PRINT_DOF_LINE] = CONTROL;

  strcpy(name[CONTROL_PRINT_DOF_LINE_COORDINATES],"control_print_dof_line_coordinates");
  type[CONTROL_PRINT_DOF_LINE_COORDINATES] = DOUBLE_PRECISION;
  data_length[CONTROL_PRINT_DOF_LINE_COORDINATES] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_PRINT_DOF_LINE_COORDINATES] = 0;
  data_class[CONTROL_PRINT_DOF_LINE_COORDINATES] = CONTROL;
  data_required[CONTROL_PRINT_DOF_LINE_COORDINATES] = CONTROL_PRINT_DOF_LINE;

  strcpy(name[CONTROL_PRINT_DOF_LINE_ELEMENT_GROUP],"control_print_dof_line_element_group");
  type[CONTROL_PRINT_DOF_LINE_ELEMENT_GROUP] = INTEGER;
  data_length[CONTROL_PRINT_DOF_LINE_ELEMENT_GROUP] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_PRINT_DOF_LINE_ELEMENT_GROUP] = 0;
  data_class[CONTROL_PRINT_DOF_LINE_ELEMENT_GROUP] = CONTROL;
  data_required[CONTROL_PRINT_DOF_LINE_ELEMENT_GROUP] = CONTROL_PRINT_DOF_LINE;

  strcpy(name[CONTROL_PRINT_DOF_LINE_EPS_ISO],"control_print_dof_line_eps_iso");
  type[CONTROL_PRINT_DOF_LINE_EPS_ISO] = DOUBLE_PRECISION;
  data_length[CONTROL_PRINT_DOF_LINE_EPS_ISO] = 1;
  data_class[CONTROL_PRINT_DOF_LINE_EPS_ISO] = CONTROL;
  data_required[CONTROL_PRINT_DOF_LINE_EPS_ISO] = CONTROL_PRINT_DOF_LINE;

  strcpy(name[CONTROL_PRINT_DOF_LINE_METHOD],"control_print_dof_line_method");
  type[CONTROL_PRINT_DOF_LINE_METHOD] = INTEGER;
  data_length[CONTROL_PRINT_DOF_LINE_METHOD] = 1;
  data_class[CONTROL_PRINT_DOF_LINE_METHOD] = CONTROL;
  data_required[CONTROL_PRINT_DOF_LINE_METHOD] = CONTROL_PRINT_DOF_LINE;

  strcpy(name[CONTROL_PRINT_DOF_LINE_MOVE],"control_print_dof_line_move");
  type[CONTROL_PRINT_DOF_LINE_MOVE] = INTEGER;
  data_length[CONTROL_PRINT_DOF_LINE_MOVE] = 1;
  data_class[CONTROL_PRINT_DOF_LINE_MOVE] = CONTROL;
  data_required[CONTROL_PRINT_DOF_LINE_MOVE] = CONTROL_PRINT_DOF_LINE;

  strcpy(name[CONTROL_PRINT_DOF_LINE_N],"control_print_dof_line_n");
  type[CONTROL_PRINT_DOF_LINE_N] = INTEGER;
  data_length[CONTROL_PRINT_DOF_LINE_N] = 1;
  data_class[CONTROL_PRINT_DOF_LINE_N] = CONTROL;
  data_required[CONTROL_PRINT_DOF_LINE_N] = CONTROL_PRINT_DOF_LINE;

  strcpy(name[CONTROL_PRINT_DOF_LINE_TIME],"control_print_dof_line_time");
  type[CONTROL_PRINT_DOF_LINE_TIME] = INTEGER;
  data_length[CONTROL_PRINT_DOF_LINE_TIME] = 1;
  data_class[CONTROL_PRINT_DOF_LINE_TIME] = CONTROL;
  data_required[CONTROL_PRINT_DOF_LINE_TIME] = CONTROL_PRINT_DOF_LINE;

  strcpy(name[CONTROL_PRINT_DOF_POINT],"control_print_dof_point");
  type[CONTROL_PRINT_DOF_POINT] = INTEGER;
  data_length[CONTROL_PRINT_DOF_POINT] = 1;
  fixed_length[CONTROL_PRINT_DOF_POINT] = 0;
  data_class[CONTROL_PRINT_DOF_POINT] = CONTROL;

  strcpy(name[CONTROL_PRINT_DOF_POINT_COORDINATES],"control_print_dof_point_coordinates");
  type[CONTROL_PRINT_DOF_POINT_COORDINATES] = DOUBLE_PRECISION;
  data_length[CONTROL_PRINT_DOF_POINT_COORDINATES] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_PRINT_DOF_POINT_COORDINATES] = 0;
  data_class[CONTROL_PRINT_DOF_POINT_COORDINATES] = CONTROL;
  data_required[CONTROL_PRINT_DOF_POINT_COORDINATES] = CONTROL_PRINT_DOF_POINT;

  strcpy(name[CONTROL_PRINT_DOF_POINT_TIME],"control_print_dof_point_time");
  type[CONTROL_PRINT_DOF_POINT_TIME] = INTEGER;
  data_length[CONTROL_PRINT_DOF_POINT_TIME] = 1;
  data_class[CONTROL_PRINT_DOF_POINT_TIME] = CONTROL;
  data_required[CONTROL_PRINT_DOF_POINT_TIME] = CONTROL_PRINT_DOF_POINT;

  strcpy(name[CONTROL_PRINT_DOF_SMOOTH_DOF],"control_print_dof_smooth_dof");
  type[CONTROL_PRINT_DOF_SMOOTH_DOF] = INTEGER;
  data_length[CONTROL_PRINT_DOF_SMOOTH_DOF] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_PRINT_DOF_SMOOTH_DOF] = 0;
  data_class[CONTROL_PRINT_DOF_SMOOTH_DOF] = CONTROL;
  data_required[CONTROL_PRINT_DOF_SMOOTH_DOF] = CONTROL_PRINT_DOF;

  strcpy(name[CONTROL_PRINT_DOF_SMOOTH_N],"control_print_dof_smooth_n");
  type[CONTROL_PRINT_DOF_SMOOTH_N] = INTEGER;
  data_length[CONTROL_PRINT_DOF_SMOOTH_N] = 1;
  data_class[CONTROL_PRINT_DOF_SMOOTH_N] = CONTROL;
  data_required[CONTROL_PRINT_DOF_SMOOTH_N] = CONTROL_PRINT_DOF_SMOOTH_DOF;

  strcpy(name[CONTROL_PRINT_NODE],"control_print_node");
  type[CONTROL_PRINT_NODE] = INTEGER;
  data_length[CONTROL_PRINT_NODE] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_PRINT_NODE] = 0;
  data_class[CONTROL_PRINT_NODE] = CONTROL;

  strcpy(name[CONTROL_PRINT_NODE_ANGULAR],"control_print_node_angular");
  type[CONTROL_PRINT_NODE_ANGULAR] = INTEGER;
  data_length[CONTROL_PRINT_NODE_ANGULAR] = 3;
  fixed_length[CONTROL_PRINT_NODE_ANGULAR] = 0;
  data_class[CONTROL_PRINT_NODE_ANGULAR] = CONTROL;
  data_required[CONTROL_PRINT_NODE_ANGULAR] = CONTROL_PRINT_NODE;

  strcpy(name[CONTROL_PRINT_NODE_ANGULAR_MIDDLE],"control_print_node_angular_middle");
  type[CONTROL_PRINT_NODE_ANGULAR_MIDDLE] = DOUBLE_PRECISION;
  data_length[CONTROL_PRINT_NODE_ANGULAR_MIDDLE] = 3;
  fixed_length[CONTROL_PRINT_NODE_ANGULAR_MIDDLE] = 0;
  data_class[CONTROL_PRINT_NODE_ANGULAR_MIDDLE] = CONTROL;
  data_required[CONTROL_PRINT_NODE_ANGULAR_MIDDLE] = CONTROL_PRINT_NODE_ANGULAR;

  strcpy(name[CONTROL_PRINT_NODE_GEOMETRY],"control_print_node_geometry");
  type[CONTROL_PRINT_NODE_GEOMETRY] = INTEGER;
  data_length[CONTROL_PRINT_NODE_GEOMETRY] = 2;
  data_class[CONTROL_PRINT_NODE_GEOMETRY] = CONTROL;
  data_required[CONTROL_PRINT_NODE_GEOMETRY] = CONTROL_PRINT_NODE;

  strcpy(name[CONTROL_PRINT_NODE_SORT],"control_print_node_sort");
  type[CONTROL_PRINT_NODE_SORT] = INTEGER;
  data_length[CONTROL_PRINT_NODE_SORT] = 1;
  data_class[CONTROL_PRINT_NODE_SORT] = CONTROL;
  data_required[CONTROL_PRINT_NODE_SORT] = CONTROL_PRINT_NODE;

  strcpy(name[CONTROL_PRINT_NODE_ZERO],"control_print_node_zero");
  type[CONTROL_PRINT_NODE_ZERO] = INTEGER;
  data_length[CONTROL_PRINT_NODE_ZERO] = 1;
  data_class[CONTROL_PRINT_NODE_ZERO] = CONTROL;
  data_required[CONTROL_PRINT_NODE_ZERO] = CONTROL_PRINT_NODE;

  strcpy(name[CONTROL_PRINT_FREQUENCY_TIMEINTERVAL],"control_print_frequency_timeinterval");
  type[CONTROL_PRINT_FREQUENCY_TIMEINTERVAL] = DOUBLE_PRECISION;
  data_length[CONTROL_PRINT_FREQUENCY_TIMEINTERVAL] = 1;
  data_class[CONTROL_PRINT_FREQUENCY_TIMEINTERVAL] = CONTROL;
  data_required[CONTROL_PRINT_FREQUENCY_TIMEINTERVAL] = CONTROL_TIMESTEP;

  strcpy(name[CONTROL_PRINT_FREQUENCY_TIMEINTERVAL_TIME],"control_print_frequency_timeinterval_time");
  type[CONTROL_PRINT_FREQUENCY_TIMEINTERVAL_TIME] = DOUBLE_PRECISION;
  data_length[CONTROL_PRINT_FREQUENCY_TIMEINTERVAL_TIME] = 1;
  external[CONTROL_PRINT_FREQUENCY_TIMEINTERVAL_TIME] = 0;
  data_class[CONTROL_PRINT_FREQUENCY_TIMEINTERVAL_TIME] = CONTROL;

  strcpy(name[CONTROL_PRINT_FREQUENCY_TIMESTEP],"control_print_frequency_timestep");
  type[CONTROL_PRINT_FREQUENCY_TIMESTEP] = INTEGER;
  data_length[CONTROL_PRINT_FREQUENCY_TIMESTEP] = 1;
  data_class[CONTROL_PRINT_FREQUENCY_TIMESTEP] = CONTROL;
  data_required[CONTROL_PRINT_FREQUENCY_TIMESTEP] = CONTROL_TIMESTEP;

  strcpy(name[CONTROL_PRINT_FREQUENCY_TIMESTEP_COUNT],"control_print_frequency_timestep_count");
  type[CONTROL_PRINT_FREQUENCY_TIMESTEP_COUNT] = INTEGER;
  data_length[CONTROL_PRINT_FREQUENCY_TIMESTEP_COUNT] = 1;
  external[CONTROL_PRINT_FREQUENCY_TIMESTEP_COUNT] = 0;
  data_class[CONTROL_PRINT_FREQUENCY_TIMESTEP_COUNT] = CONTROL;

  strcpy(name[CONTROL_PRINT_BEAM_FORCE_MOMENT],"control_print_beam_force_moment");
  type[CONTROL_PRINT_BEAM_FORCE_MOMENT] = INTEGER;
  data_length[CONTROL_PRINT_BEAM_FORCE_MOMENT] = 1;
  data_class[CONTROL_PRINT_BEAM_FORCE_MOMENT] = CONTROL;

  strcpy(name[CONTROL_PRINT_BEAM_FORCE_MOMENT_COORDINATES],"control_print_beam_force_moment_coordinates");
  type[CONTROL_PRINT_BEAM_FORCE_MOMENT_COORDINATES] = DOUBLE_PRECISION;
  data_length[CONTROL_PRINT_BEAM_FORCE_MOMENT_COORDINATES] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_PRINT_BEAM_FORCE_MOMENT_COORDINATES] = 0;
  data_class[CONTROL_PRINT_BEAM_FORCE_MOMENT_COORDINATES] = CONTROL;
  data_required[CONTROL_PRINT_BEAM_FORCE_MOMENT_COORDINATES] = CONTROL_PRINT_BEAM_FORCE_MOMENT;

  strcpy(name[CONTROL_PRINT_BEAM_FORCE_MOMENT_SWITCH],"control_print_beam_force_moment_switch");
  type[CONTROL_PRINT_BEAM_FORCE_MOMENT_SWITCH] = INTEGER;
  data_length[CONTROL_PRINT_BEAM_FORCE_MOMENT_SWITCH] = 1;
  data_class[CONTROL_PRINT_BEAM_FORCE_MOMENT_SWITCH] = CONTROL;
  data_required[CONTROL_PRINT_BEAM_FORCE_MOMENT_SWITCH] = CONTROL_PRINT_BEAM_FORCE_MOMENT;

  strcpy(name[CONTROL_PRINT_MATERI_STRESS_FORCE],"control_print_materi_stress_force");
  type[CONTROL_PRINT_MATERI_STRESS_FORCE] = INTEGER;
  data_length[CONTROL_PRINT_MATERI_STRESS_FORCE] = 1;
  data_class[CONTROL_PRINT_MATERI_STRESS_FORCE] = CONTROL;

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

  strcpy(name[CONTROL_REPEAT_SAVE],"control_repeat_save");
  type[CONTROL_REPEAT_SAVE] = INTEGER;
  data_length[CONTROL_REPEAT_SAVE] = DATA_ITEM_SIZE;
  fixed_length[CONTROL_REPEAT_SAVE] = 0;
  data_class[CONTROL_REPEAT_SAVE] = CONTROL;
  data_required[CONTROL_REPEAT_SAVE] = CONTROL_REPEAT;

  strcpy(name[CONTROL_REPEAT_SAVE_CALCULATE],"control_repeat_save_calculate");
  type[CONTROL_REPEAT_SAVE_CALCULATE] = INTEGER;
  data_length[CONTROL_REPEAT_SAVE_CALCULATE] = 1;
  data_class[CONTROL_REPEAT_SAVE_CALCULATE] = CONTROL;
  data_required[CONTROL_REPEAT_SAVE_CALCULATE] = CONTROL_REPEAT;

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
  // manual Professional 6.386: ratio_criterium minimal_timestep
  // maximum_timestep (3 values). The GNU used 2 (ratio, maximum) -
  // the corpus slope tests give all three.
  data_length[CONTROL_TIMESTEP_ITERATIONS_AUTOMATIC] = 3;
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

  strcpy(name[DATA_ACTIVATE],"data_activate");
  type[DATA_ACTIVATE] = INTEGER;
  data_length[DATA_ACTIVATE] = DATA_ITEM_SIZE;
  fixed_length[DATA_ACTIVATE] = 0;
  data_class[DATA_ACTIVATE] = CONTROL;

  strcpy(name[DATA_ACTIVATE_TIME],"data_activate_time");
  type[DATA_ACTIVATE_TIME] = DOUBLE_PRECISION;
  data_length[DATA_ACTIVATE_TIME] = 1;
  data_class[DATA_ACTIVATE_TIME] = CONTROL;
  data_required[DATA_ACTIVATE_TIME] = DATA_ACTIVATE;

  strcpy(name[DATA_DELETE],"data_delete");
  type[DATA_DELETE] = INTEGER;
  data_length[DATA_DELETE] = DATA_ITEM_SIZE;
  fixed_length[DATA_DELETE] = 0;
  data_class[DATA_DELETE] = CONTROL;

  strcpy(name[DATA_DELETE_TIME],"data_delete_time");
  type[DATA_DELETE_TIME] = DOUBLE_PRECISION;
  data_length[DATA_DELETE_TIME] = 1;
  data_class[DATA_DELETE_TIME] = CONTROL;
  data_required[DATA_DELETE_TIME] = DATA_DELETE;

  strcpy(name[DATA_IGNORE],"data_ignore");
  type[DATA_IGNORE] = INTEGER;
  data_length[DATA_IGNORE] = 1;
  no_index[DATA_IGNORE] = 1;
  data_class[DATA_IGNORE] = CONTROL;

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
  data_class[ELEMENT_EMPTY] = ELEMENT;
  data_required[ELEMENT_EMPTY] = ELEMENT;

  strcpy(name[ELEMENT_GEOMETRY],"element_geometry");
  type[ELEMENT_GEOMETRY] = INTEGER;
  data_length[ELEMENT_GEOMETRY] = 1;
  version_all[ELEMENT_GEOMETRY] = 1;
  data_class[ELEMENT_GEOMETRY] = ELEMENT;
  data_required[ELEMENT_GEOMETRY] = ELEMENT;

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

  strcpy(name[ELEMENT_INTPNT_MATERI_PLASTI_HARDSOIL_GAMMAP_INITIAL],
    "element_intpnt_materi_plasti_hardsoil_gammap_initial");
  type[ELEMENT_INTPNT_MATERI_PLASTI_HARDSOIL_GAMMAP_INITIAL] = DOUBLE_PRECISION;
  data_length[ELEMENT_INTPNT_MATERI_PLASTI_HARDSOIL_GAMMAP_INITIAL] = 1;
  data_class[ELEMENT_INTPNT_MATERI_PLASTI_HARDSOIL_GAMMAP_INITIAL] = ELEMENT;
  data_required[ELEMENT_INTPNT_MATERI_PLASTI_HARDSOIL_GAMMAP_INITIAL] = ELEMENT;

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

  strcpy(name[ELEMENT_INTERFACE_STRAIN_NORMAL],"element_interface_strain_normal");
  type[ELEMENT_INTERFACE_STRAIN_NORMAL] = DOUBLE_PRECISION;
  data_length[ELEMENT_INTERFACE_STRAIN_NORMAL] = 4;
  fixed_length[ELEMENT_INTERFACE_STRAIN_NORMAL] = 0;
  version_all[ELEMENT_INTERFACE_STRAIN_NORMAL] = 1;
  print_only[ELEMENT_INTERFACE_STRAIN_NORMAL] = 1;
  data_class[ELEMENT_INTERFACE_STRAIN_NORMAL] = ELEMENT;
  data_required[ELEMENT_INTERFACE_STRAIN_NORMAL] = ELEMENT;

  strcpy(name[ELEMENT_INTERFACE_FORCE_TANG],"element_interface_force_tang");
  type[ELEMENT_INTERFACE_FORCE_TANG] = DOUBLE_PRECISION;
  data_length[ELEMENT_INTERFACE_FORCE_TANG] = 4;
  fixed_length[ELEMENT_INTERFACE_FORCE_TANG] = 0;
  version_all[ELEMENT_INTERFACE_FORCE_TANG] = 1;
  print_only[ELEMENT_INTERFACE_FORCE_TANG] = 1;
  data_class[ELEMENT_INTERFACE_FORCE_TANG] = ELEMENT;
  data_required[ELEMENT_INTERFACE_FORCE_TANG] = ELEMENT;

  strcpy(name[ELEMENT_INTERFACE_FORCE_TANG2],"element_interface_force_tang2");
  type[ELEMENT_INTERFACE_FORCE_TANG2] = DOUBLE_PRECISION;
  data_length[ELEMENT_INTERFACE_FORCE_TANG2] = 4;
  fixed_length[ELEMENT_INTERFACE_FORCE_TANG2] = 0;
  version_all[ELEMENT_INTERFACE_FORCE_TANG2] = 1;
  print_only[ELEMENT_INTERFACE_FORCE_TANG2] = 1;
  data_class[ELEMENT_INTERFACE_FORCE_TANG2] = ELEMENT;
  data_required[ELEMENT_INTERFACE_FORCE_TANG2] = ELEMENT;

  strcpy(name[ELEMENT_INTERFACE_FORCE_NORM],"element_interface_force_norm");
  type[ELEMENT_INTERFACE_FORCE_NORM] = DOUBLE_PRECISION;
  data_length[ELEMENT_INTERFACE_FORCE_NORM] = 4;
  fixed_length[ELEMENT_INTERFACE_FORCE_NORM] = 0;
  version_all[ELEMENT_INTERFACE_FORCE_NORM] = 1;
  print_only[ELEMENT_INTERFACE_FORCE_NORM] = 1;
  data_class[ELEMENT_INTERFACE_FORCE_NORM] = ELEMENT;
  data_required[ELEMENT_INTERFACE_FORCE_NORM] = ELEMENT;

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
  // variable length: one value per space direction (Professional
  // force_edge, manual 6.454) or one value per principal dof (legacy)
  fixed_length[FORCE_ELEMENT_EDGE] = 0;

  strcpy(name[FORCE_ELEMENT_EDGE_ELEMENT],"force_element_edge_element");
  type[FORCE_ELEMENT_EDGE_ELEMENT] = INTEGER;
  data_length[FORCE_ELEMENT_EDGE_ELEMENT] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_ELEMENT] = 0;
  data_class[FORCE_ELEMENT_EDGE_ELEMENT] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_ELEMENT] = FORCE_ELEMENT_EDGE;

  strcpy(name[FORCE_ELEMENT_EDGE_ELEMENT_GROUP],"force_element_edge_element_group");
  type[FORCE_ELEMENT_EDGE_ELEMENT_GROUP] = INTEGER;
  data_length[FORCE_ELEMENT_EDGE_ELEMENT_GROUP] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_ELEMENT_GROUP] = 0;
  data_class[FORCE_ELEMENT_EDGE_ELEMENT_GROUP] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_ELEMENT_GROUP] = FORCE_ELEMENT_EDGE;

  strcpy(name[FORCE_ELEMENT_EDGE_ELEMENT_NODE],"force_element_edge_element_node");
  type[FORCE_ELEMENT_EDGE_ELEMENT_NODE] = INTEGER;
  data_length[FORCE_ELEMENT_EDGE_ELEMENT_NODE] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_ELEMENT_NODE] = 0;
  data_class[FORCE_ELEMENT_EDGE_ELEMENT_NODE] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_ELEMENT_NODE] = FORCE_ELEMENT_EDGE;

  strcpy(name[FORCE_ELEMENT_EDGE_ELEMENT_SIDE],"force_element_edge_element_side");
  type[FORCE_ELEMENT_EDGE_ELEMENT_SIDE] = INTEGER;
  data_length[FORCE_ELEMENT_EDGE_ELEMENT_SIDE] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_ELEMENT_SIDE] = 0;
  data_class[FORCE_ELEMENT_EDGE_ELEMENT_SIDE] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_ELEMENT_SIDE] = FORCE_ELEMENT_EDGE;

  strcpy(name[FORCE_ELEMENT_EDGE_NODE],"force_element_edge_node");
  type[FORCE_ELEMENT_EDGE_NODE] = INTEGER;
  data_length[FORCE_ELEMENT_EDGE_NODE] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_NODE] = 0;
  data_class[FORCE_ELEMENT_EDGE_NODE] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_NODE] = FORCE_ELEMENT_EDGE;

  strcpy(name[FORCE_ELEMENT_EDGE_NODE_FACTOR],"force_element_edge_node_factor");
  type[FORCE_ELEMENT_EDGE_NODE_FACTOR] = DOUBLE_PRECISION;
  data_length[FORCE_ELEMENT_EDGE_NODE_FACTOR] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_NODE_FACTOR] = 0;
  data_class[FORCE_ELEMENT_EDGE_NODE_FACTOR] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_NODE_FACTOR] = FORCE_ELEMENT_EDGE;

  strcpy(name[FORCE_ELEMENT_EDGE_FACTOR],"force_element_edge_factor");
  type[FORCE_ELEMENT_EDGE_FACTOR] = DOUBLE_PRECISION;
  data_length[FORCE_ELEMENT_EDGE_FACTOR] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_FACTOR] = 0;
  data_class[FORCE_ELEMENT_EDGE_FACTOR] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_FACTOR] = FORCE_ELEMENT_EDGE;

  strcpy(name[FORCE_ELEMENT_EDGE_MULTI_LINEAR_FACTOR_X],"force_element_edge_multi_linear_factor_x");
  type[FORCE_ELEMENT_EDGE_MULTI_LINEAR_FACTOR_X] = DOUBLE_PRECISION;
  data_length[FORCE_ELEMENT_EDGE_MULTI_LINEAR_FACTOR_X] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_MULTI_LINEAR_FACTOR_X] = 0;
  data_class[FORCE_ELEMENT_EDGE_MULTI_LINEAR_FACTOR_X] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_MULTI_LINEAR_FACTOR_X] = FORCE_ELEMENT_EDGE;

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

  strcpy(name[FORCE_ELEMENT_EDGE_NORMAL_ELEMENT],"force_element_edge_normal_element");
  type[FORCE_ELEMENT_EDGE_NORMAL_ELEMENT] = INTEGER;
  data_length[FORCE_ELEMENT_EDGE_NORMAL_ELEMENT] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_NORMAL_ELEMENT] = 0;
  data_class[FORCE_ELEMENT_EDGE_NORMAL_ELEMENT] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_NORMAL_ELEMENT] = FORCE_ELEMENT_EDGE_NORMAL;

  strcpy(name[FORCE_ELEMENT_EDGE_NORMAL_ELEMENT_GROUP],"force_element_edge_normal_element_group");
  type[FORCE_ELEMENT_EDGE_NORMAL_ELEMENT_GROUP] = INTEGER;
  data_length[FORCE_ELEMENT_EDGE_NORMAL_ELEMENT_GROUP] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_NORMAL_ELEMENT_GROUP] = 0;
  data_class[FORCE_ELEMENT_EDGE_NORMAL_ELEMENT_GROUP] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_NORMAL_ELEMENT_GROUP] = FORCE_ELEMENT_EDGE_NORMAL;

  strcpy(name[FORCE_ELEMENT_EDGE_NORMAL_ELEMENT_NODE],"force_element_edge_normal_element_node");
  type[FORCE_ELEMENT_EDGE_NORMAL_ELEMENT_NODE] = INTEGER;
  data_length[FORCE_ELEMENT_EDGE_NORMAL_ELEMENT_NODE] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_NORMAL_ELEMENT_NODE] = 0;
  data_class[FORCE_ELEMENT_EDGE_NORMAL_ELEMENT_NODE] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_NORMAL_ELEMENT_NODE] = FORCE_ELEMENT_EDGE_NORMAL;

  strcpy(name[FORCE_ELEMENT_EDGE_NORMAL_ELEMENT_SIDE],"force_element_edge_normal_element_side");
  type[FORCE_ELEMENT_EDGE_NORMAL_ELEMENT_SIDE] = INTEGER;
  data_length[FORCE_ELEMENT_EDGE_NORMAL_ELEMENT_SIDE] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_NORMAL_ELEMENT_SIDE] = 0;
  data_class[FORCE_ELEMENT_EDGE_NORMAL_ELEMENT_SIDE] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_NORMAL_ELEMENT_SIDE] = FORCE_ELEMENT_EDGE_NORMAL;

  strcpy(name[FORCE_ELEMENT_EDGE_NORMAL_NODE],"force_element_edge_normal_node");
  type[FORCE_ELEMENT_EDGE_NORMAL_NODE] = INTEGER;
  data_length[FORCE_ELEMENT_EDGE_NORMAL_NODE] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_NORMAL_NODE] = 0;
  data_class[FORCE_ELEMENT_EDGE_NORMAL_NODE] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_NORMAL_NODE] = FORCE_ELEMENT_EDGE_NORMAL;

  strcpy(name[FORCE_ELEMENT_EDGE_NORMAL_NODE_FACTOR],"force_element_edge_normal_node_factor");
  type[FORCE_ELEMENT_EDGE_NORMAL_NODE_FACTOR] = DOUBLE_PRECISION;
  data_length[FORCE_ELEMENT_EDGE_NORMAL_NODE_FACTOR] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_NORMAL_NODE_FACTOR] = 0;
  data_class[FORCE_ELEMENT_EDGE_NORMAL_NODE_FACTOR] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_NORMAL_NODE_FACTOR] = FORCE_ELEMENT_EDGE_NORMAL;

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

  strcpy(name[FORCE_ELEMENT_EDGE_PROJECTED],"force_element_edge_projected");
  type[FORCE_ELEMENT_EDGE_PROJECTED] = DOUBLE_PRECISION;
  if      ( ndim==2 )
    data_length[FORCE_ELEMENT_EDGE_PROJECTED] = 10;
  else
    data_length[FORCE_ELEMENT_EDGE_PROJECTED] = 16;
  data_class[FORCE_ELEMENT_EDGE_PROJECTED] = FORCE;

  strcpy(name[FORCE_ELEMENT_EDGE_PROJECTED_ELEMENT],"force_element_edge_projected_element");
  type[FORCE_ELEMENT_EDGE_PROJECTED_ELEMENT] = INTEGER;
  data_length[FORCE_ELEMENT_EDGE_PROJECTED_ELEMENT] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_PROJECTED_ELEMENT] = 0;
  data_class[FORCE_ELEMENT_EDGE_PROJECTED_ELEMENT] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_PROJECTED_ELEMENT] = FORCE_ELEMENT_EDGE_PROJECTED;

  strcpy(name[FORCE_ELEMENT_EDGE_PROJECTED_ELEMENT_GROUP],"force_element_edge_projected_element_group");
  type[FORCE_ELEMENT_EDGE_PROJECTED_ELEMENT_GROUP] = INTEGER;
  data_length[FORCE_ELEMENT_EDGE_PROJECTED_ELEMENT_GROUP] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_PROJECTED_ELEMENT_GROUP] = 0;
  data_class[FORCE_ELEMENT_EDGE_PROJECTED_ELEMENT_GROUP] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_PROJECTED_ELEMENT_GROUP] = FORCE_ELEMENT_EDGE_PROJECTED;

  strcpy(name[FORCE_ELEMENT_EDGE_PROJECTED_ELEMENT_NODE],"force_element_edge_projected_element_node");
  type[FORCE_ELEMENT_EDGE_PROJECTED_ELEMENT_NODE] = INTEGER;
  data_length[FORCE_ELEMENT_EDGE_PROJECTED_ELEMENT_NODE] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_PROJECTED_ELEMENT_NODE] = 0;
  data_class[FORCE_ELEMENT_EDGE_PROJECTED_ELEMENT_NODE] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_PROJECTED_ELEMENT_NODE] = FORCE_ELEMENT_EDGE_PROJECTED;

  strcpy(name[FORCE_ELEMENT_EDGE_PROJECTED_ELEMENT_SIDE],"force_element_edge_projected_element_side");
  type[FORCE_ELEMENT_EDGE_PROJECTED_ELEMENT_SIDE] = INTEGER;
  data_length[FORCE_ELEMENT_EDGE_PROJECTED_ELEMENT_SIDE] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_PROJECTED_ELEMENT_SIDE] = 0;
  data_class[FORCE_ELEMENT_EDGE_PROJECTED_ELEMENT_SIDE] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_PROJECTED_ELEMENT_SIDE] = FORCE_ELEMENT_EDGE_PROJECTED;

  strcpy(name[FORCE_ELEMENT_EDGE_PROJECTED_FACTOR],"force_element_edge_projected_factor");
  type[FORCE_ELEMENT_EDGE_PROJECTED_FACTOR] = DOUBLE_PRECISION;
  data_length[FORCE_ELEMENT_EDGE_PROJECTED_FACTOR] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_PROJECTED_FACTOR] = 0;
  data_class[FORCE_ELEMENT_EDGE_PROJECTED_FACTOR] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_PROJECTED_FACTOR] = FORCE_ELEMENT_EDGE_PROJECTED;

  strcpy(name[FORCE_ELEMENT_EDGE_PROJECTED_GEOMETRY],"force_element_edge_projected_geometry");
  type[FORCE_ELEMENT_EDGE_PROJECTED_GEOMETRY] = INTEGER;
  data_length[FORCE_ELEMENT_EDGE_PROJECTED_GEOMETRY] = 2;
  data_class[FORCE_ELEMENT_EDGE_PROJECTED_GEOMETRY] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_PROJECTED_GEOMETRY] = FORCE_ELEMENT_EDGE_PROJECTED;

  strcpy(name[FORCE_ELEMENT_EDGE_PROJECTED_NODE],"force_element_edge_projected_node");
  type[FORCE_ELEMENT_EDGE_PROJECTED_NODE] = INTEGER;
  data_length[FORCE_ELEMENT_EDGE_PROJECTED_NODE] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_PROJECTED_NODE] = 0;
  data_class[FORCE_ELEMENT_EDGE_PROJECTED_NODE] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_PROJECTED_NODE] = FORCE_ELEMENT_EDGE_PROJECTED;

  strcpy(name[FORCE_ELEMENT_EDGE_PROJECTED_NODE_FACTOR],"force_element_edge_projected_node_factor");
  type[FORCE_ELEMENT_EDGE_PROJECTED_NODE_FACTOR] = DOUBLE_PRECISION;
  data_length[FORCE_ELEMENT_EDGE_PROJECTED_NODE_FACTOR] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_PROJECTED_NODE_FACTOR] = 0;
  data_class[FORCE_ELEMENT_EDGE_PROJECTED_NODE_FACTOR] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_PROJECTED_NODE_FACTOR] = FORCE_ELEMENT_EDGE_PROJECTED;

  strcpy(name[FORCE_ELEMENT_EDGE_PROJECTED_SINE],"force_element_edge_projected_sine");
  type[FORCE_ELEMENT_EDGE_PROJECTED_SINE] = DOUBLE_PRECISION;
  data_length[FORCE_ELEMENT_EDGE_PROJECTED_SINE] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_PROJECTED_SINE] = 0;
  data_class[FORCE_ELEMENT_EDGE_PROJECTED_SINE] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_PROJECTED_SINE] = FORCE_ELEMENT_EDGE_PROJECTED;

  strcpy(name[FORCE_ELEMENT_EDGE_PROJECTED_TIME],"force_element_edge_projected_time");
  type[FORCE_ELEMENT_EDGE_PROJECTED_TIME] = DOUBLE_PRECISION;
  data_length[FORCE_ELEMENT_EDGE_PROJECTED_TIME] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_PROJECTED_TIME] = 0;
  data_class[FORCE_ELEMENT_EDGE_PROJECTED_TIME] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_PROJECTED_TIME] = FORCE_ELEMENT_EDGE_PROJECTED;

  strcpy(name[FORCE_ELEMENT_EDGE_WATER],"force_element_edge_water");
  type[FORCE_ELEMENT_EDGE_WATER] = DOUBLE_PRECISION;
  if      ( ndim==2 )
    data_length[FORCE_ELEMENT_EDGE_WATER] = 4;
  else if ( ndim==3 )
    data_length[FORCE_ELEMENT_EDGE_WATER] = 5;
  else
    data_length[FORCE_ELEMENT_EDGE_WATER] = 3;
  // variable length: legacy GNU layout (rho g dirx [diry dirz]) or the
  // single Professional switch force_edge_water index -yes (auto
  // hydrostatic; manual 6.489). The parser stores -yes as the negative
  // YES enum value (see input.cc).
  fixed_length[FORCE_ELEMENT_EDGE_WATER] = 0;
  data_class[FORCE_ELEMENT_EDGE_WATER] = FORCE;

  strcpy(name[FORCE_ELEMENT_EDGE_WATER_ELEMENT],"force_element_edge_water_element");
  type[FORCE_ELEMENT_EDGE_WATER_ELEMENT] = INTEGER;
  data_length[FORCE_ELEMENT_EDGE_WATER_ELEMENT] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_WATER_ELEMENT] = 0;
  data_class[FORCE_ELEMENT_EDGE_WATER_ELEMENT] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_WATER_ELEMENT] = FORCE_ELEMENT_EDGE_WATER;

  strcpy(name[FORCE_ELEMENT_EDGE_WATER_ELEMENT_GROUP],"force_element_edge_water_element_group");
  type[FORCE_ELEMENT_EDGE_WATER_ELEMENT_GROUP] = INTEGER;
  data_length[FORCE_ELEMENT_EDGE_WATER_ELEMENT_GROUP] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_WATER_ELEMENT_GROUP] = 0;
  data_class[FORCE_ELEMENT_EDGE_WATER_ELEMENT_GROUP] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_WATER_ELEMENT_GROUP] = FORCE_ELEMENT_EDGE_WATER;

  strcpy(name[FORCE_ELEMENT_EDGE_WATER_ELEMENT_NODE],"force_element_edge_water_element_node");
  type[FORCE_ELEMENT_EDGE_WATER_ELEMENT_NODE] = INTEGER;
  data_length[FORCE_ELEMENT_EDGE_WATER_ELEMENT_NODE] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_WATER_ELEMENT_NODE] = 0;
  data_class[FORCE_ELEMENT_EDGE_WATER_ELEMENT_NODE] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_WATER_ELEMENT_NODE] = FORCE_ELEMENT_EDGE_WATER;

  strcpy(name[FORCE_ELEMENT_EDGE_WATER_ELEMENT_SIDE],"force_element_edge_water_element_side");
  type[FORCE_ELEMENT_EDGE_WATER_ELEMENT_SIDE] = INTEGER;
  data_length[FORCE_ELEMENT_EDGE_WATER_ELEMENT_SIDE] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_WATER_ELEMENT_SIDE] = 0;
  data_class[FORCE_ELEMENT_EDGE_WATER_ELEMENT_SIDE] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_WATER_ELEMENT_SIDE] = FORCE_ELEMENT_EDGE_WATER;

  strcpy(name[FORCE_ELEMENT_EDGE_WATER_FACTOR],"force_element_edge_water_factor");
  type[FORCE_ELEMENT_EDGE_WATER_FACTOR] = DOUBLE_PRECISION;
  data_length[FORCE_ELEMENT_EDGE_WATER_FACTOR] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_WATER_FACTOR] = 0;
  data_class[FORCE_ELEMENT_EDGE_WATER_FACTOR] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_WATER_FACTOR] = FORCE_ELEMENT_EDGE_WATER;

  strcpy(name[FORCE_ELEMENT_EDGE_WATER_NODE],"force_element_edge_water_node");
  type[FORCE_ELEMENT_EDGE_WATER_NODE] = INTEGER;
  data_length[FORCE_ELEMENT_EDGE_WATER_NODE] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_EDGE_WATER_NODE] = 0;
  data_class[FORCE_ELEMENT_EDGE_WATER_NODE] = FORCE;
  data_required[FORCE_ELEMENT_EDGE_WATER_NODE] = FORCE_ELEMENT_EDGE_WATER;

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

  strcpy(name[FORCE_ELEMENT_VOLUME_ELEMENT],"force_element_volume_element");
  type[FORCE_ELEMENT_VOLUME_ELEMENT] = INTEGER;
  data_length[FORCE_ELEMENT_VOLUME_ELEMENT] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_VOLUME_ELEMENT] = 0;
  data_class[FORCE_ELEMENT_VOLUME_ELEMENT] = FORCE;
  data_required[FORCE_ELEMENT_VOLUME_ELEMENT] = FORCE_ELEMENT_VOLUME;

  strcpy(name[FORCE_ELEMENT_VOLUME_ELEMENT_GROUP],"force_element_volume_element_group");
  type[FORCE_ELEMENT_VOLUME_ELEMENT_GROUP] = INTEGER;
  data_length[FORCE_ELEMENT_VOLUME_ELEMENT_GROUP] = DATA_ITEM_SIZE;
  fixed_length[FORCE_ELEMENT_VOLUME_ELEMENT_GROUP] = 0;
  data_class[FORCE_ELEMENT_VOLUME_ELEMENT_GROUP] = FORCE;
  data_required[FORCE_ELEMENT_VOLUME_ELEMENT_GROUP] = FORCE_ELEMENT_VOLUME;

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

  strcpy(name[FORCE_POINT],"force_point");
  type[FORCE_POINT] = DOUBLE_PRECISION;
  data_length[FORCE_POINT] = ndim + MUKNWN;
  fixed_length[FORCE_POINT] = 0;
  data_class[FORCE_POINT] = FORCE_POINT;
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

  strcpy(name[GROUNDFLOW_CONSOLIDATION_APPLY],"groundflow_consolidation_apply");
  type[GROUNDFLOW_CONSOLIDATION_APPLY] = INTEGER;
  data_length[GROUNDFLOW_CONSOLIDATION_APPLY] = 1;
  no_index[GROUNDFLOW_CONSOLIDATION_APPLY] = 1;
  data_class[GROUNDFLOW_CONSOLIDATION_APPLY] = GROUNDFLOW;

  strcpy(name[GROUNDFLOW_DENSITY],"groundflow_density");
  type[GROUNDFLOW_DENSITY] = DOUBLE_PRECISION;
  data_length[GROUNDFLOW_DENSITY] = 1;
  no_index[GROUNDFLOW_DENSITY] = 1;
  data_class[GROUNDFLOW_DENSITY] = GROUNDFLOW;

  strcpy(name[GROUNDFLOW_FLUX_EDGE_NORMAL],"groundflow_flux_edge_normal");
  type[GROUNDFLOW_FLUX_EDGE_NORMAL] = DOUBLE_PRECISION;
  data_length[GROUNDFLOW_FLUX_EDGE_NORMAL] = 1;
  data_class[GROUNDFLOW_FLUX_EDGE_NORMAL] = GROUNDFLOW;

  strcpy(name[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT],"groundflow_flux_edge_normal_element");
  type[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT] = INTEGER;
  data_length[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT] = DATA_ITEM_SIZE;
  fixed_length[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT] = 0;
  data_class[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT] = GROUNDFLOW;
  data_required[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT] = GROUNDFLOW_FLUX_EDGE_NORMAL;

  strcpy(name[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_GROUP],"groundflow_flux_edge_normal_element_group");
  type[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_GROUP] = INTEGER;
  data_length[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_GROUP] = DATA_ITEM_SIZE;
  fixed_length[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_GROUP] = 0;
  data_class[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_GROUP] = GROUNDFLOW;
  data_required[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_GROUP] = GROUNDFLOW_FLUX_EDGE_NORMAL;

  strcpy(name[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_NODE],"groundflow_flux_edge_normal_element_node");
  type[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_NODE] = INTEGER;
  data_length[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_NODE] = DATA_ITEM_SIZE;
  fixed_length[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_NODE] = 0;
  data_class[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_NODE] = GROUNDFLOW;
  data_required[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_NODE] = GROUNDFLOW_FLUX_EDGE_NORMAL;

  strcpy(name[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_NODE_FACTOR],"groundflow_flux_edge_normal_element_node_factor");
  type[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_NODE_FACTOR] = DOUBLE_PRECISION;
  data_length[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_NODE_FACTOR] = DATA_ITEM_SIZE;
  fixed_length[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_NODE_FACTOR] = 0;
  data_class[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_NODE_FACTOR] = GROUNDFLOW;
  data_required[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_NODE_FACTOR] = GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_NODE;

  strcpy(name[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_SIDE],"groundflow_flux_edge_normal_element_side");
  type[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_SIDE] = INTEGER;
  data_length[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_SIDE] = DATA_ITEM_SIZE;
  fixed_length[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_SIDE] = 0;
  data_class[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_SIDE] = GROUNDFLOW;
  data_required[GROUNDFLOW_FLUX_EDGE_NORMAL_ELEMENT_SIDE] = GROUNDFLOW_FLUX_EDGE_NORMAL;

  strcpy(name[GROUNDFLOW_FLUX_EDGE_NORMAL_FACTOR],"groundflow_flux_edge_normal_factor");
  type[GROUNDFLOW_FLUX_EDGE_NORMAL_FACTOR] = DOUBLE_PRECISION;
  data_length[GROUNDFLOW_FLUX_EDGE_NORMAL_FACTOR] = DATA_ITEM_SIZE;
  fixed_length[GROUNDFLOW_FLUX_EDGE_NORMAL_FACTOR] = 0;
  data_class[GROUNDFLOW_FLUX_EDGE_NORMAL_FACTOR] = GROUNDFLOW;
  data_required[GROUNDFLOW_FLUX_EDGE_NORMAL_FACTOR] = GROUNDFLOW_FLUX_EDGE_NORMAL;

  strcpy(name[GROUNDFLOW_FLUX_EDGE_NORMAL_GEOMETRY],"groundflow_flux_edge_normal_geometry");
  type[GROUNDFLOW_FLUX_EDGE_NORMAL_GEOMETRY] = INTEGER;
  data_length[GROUNDFLOW_FLUX_EDGE_NORMAL_GEOMETRY] = 2;
  data_class[GROUNDFLOW_FLUX_EDGE_NORMAL_GEOMETRY] = GROUNDFLOW;
  data_required[GROUNDFLOW_FLUX_EDGE_NORMAL_GEOMETRY] = GROUNDFLOW_FLUX_EDGE_NORMAL;

  strcpy(name[GROUNDFLOW_FLUX_EDGE_NORMAL_NODE],"groundflow_flux_edge_normal_node");
  type[GROUNDFLOW_FLUX_EDGE_NORMAL_NODE] = INTEGER;
  data_length[GROUNDFLOW_FLUX_EDGE_NORMAL_NODE] = DATA_ITEM_SIZE;
  fixed_length[GROUNDFLOW_FLUX_EDGE_NORMAL_NODE] = 0;
  data_class[GROUNDFLOW_FLUX_EDGE_NORMAL_NODE] = GROUNDFLOW;
  data_required[GROUNDFLOW_FLUX_EDGE_NORMAL_NODE] = GROUNDFLOW_FLUX_EDGE_NORMAL;

  strcpy(name[GROUNDFLOW_FLUX_EDGE_NORMAL_SINE],"groundflow_flux_edge_normal_sine");
  type[GROUNDFLOW_FLUX_EDGE_NORMAL_SINE] = DOUBLE_PRECISION;
  data_length[GROUNDFLOW_FLUX_EDGE_NORMAL_SINE] = DATA_ITEM_SIZE;
  fixed_length[GROUNDFLOW_FLUX_EDGE_NORMAL_SINE] = 0;
  data_class[GROUNDFLOW_FLUX_EDGE_NORMAL_SINE] = GROUNDFLOW;
  data_required[GROUNDFLOW_FLUX_EDGE_NORMAL_SINE] = GROUNDFLOW_FLUX_EDGE_NORMAL;

  strcpy(name[GROUNDFLOW_FLUX_EDGE_NORMAL_TIME],"groundflow_flux_edge_normal_time");
  type[GROUNDFLOW_FLUX_EDGE_NORMAL_TIME] = DOUBLE_PRECISION;
  data_length[GROUNDFLOW_FLUX_EDGE_NORMAL_TIME] = DATA_ITEM_SIZE;
  fixed_length[GROUNDFLOW_FLUX_EDGE_NORMAL_TIME] = 0;
  data_class[GROUNDFLOW_FLUX_EDGE_NORMAL_TIME] = GROUNDFLOW;
  data_required[GROUNDFLOW_FLUX_EDGE_NORMAL_TIME] = GROUNDFLOW_FLUX_EDGE_NORMAL;

  strcpy(name[GROUNDFLOW_NONSATURATED_APPLY],"groundflow_nonsaturated_apply");
  type[GROUNDFLOW_NONSATURATED_APPLY] = INTEGER;
  data_length[GROUNDFLOW_NONSATURATED_APPLY] = 1;
  no_index[GROUNDFLOW_NONSATURATED_APPLY] = 1;
  data_class[GROUNDFLOW_NONSATURATED_APPLY] = GROUNDFLOW;

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

  strcpy(name[GROUNDFLOW_PHREATICLEVEL_MULTIPLE],"groundflow_phreatic_level_multiple");
  type[GROUNDFLOW_PHREATICLEVEL_MULTIPLE] = DOUBLE_PRECISION;
  data_length[GROUNDFLOW_PHREATICLEVEL_MULTIPLE] = DATA_ITEM_SIZE;
  fixed_length[GROUNDFLOW_PHREATICLEVEL_MULTIPLE] = 0;
  data_class[GROUNDFLOW_PHREATICLEVEL_MULTIPLE] = GROUNDFLOW;

  strcpy(name[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT],"groundflow_phreatic_level_multiple_element");
  type[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT] = INTEGER;
  data_length[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT] = DATA_ITEM_SIZE;
  fixed_length[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT] = 0;
  data_class[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT] = GROUNDFLOW;
  data_required[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT] = GROUNDFLOW_PHREATICLEVEL_MULTIPLE;

  strcpy(name[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GEOMETRY],"groundflow_phreatic_level_multiple_element_geometry");
  type[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GEOMETRY] = INTEGER;
  data_length[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GEOMETRY] = DATA_ITEM_SIZE;
  fixed_length[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GEOMETRY] = 0;
  data_class[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GEOMETRY] = GROUNDFLOW;
  data_required[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GEOMETRY] = GROUNDFLOW_PHREATICLEVEL_MULTIPLE;

  strcpy(name[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GROUP],"groundflow_phreatic_level_multiple_element_group");
  type[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GROUP] = INTEGER;
  data_length[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GROUP] = DATA_ITEM_SIZE;
  fixed_length[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GROUP] = 0;
  data_class[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GROUP] = GROUNDFLOW;
  data_required[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GROUP] = GROUNDFLOW_PHREATICLEVEL_MULTIPLE;

  strcpy(name[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_N],"groundflow_phreatic_level_multiple_n");
  type[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_N] = INTEGER;
  data_length[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_N] = 2;
  data_class[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_N] = GROUNDFLOW;
  data_required[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_N] = GROUNDFLOW_PHREATICLEVEL_MULTIPLE;

  strcpy(name[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_NODE],"groundflow_phreatic_level_multiple_node");
  type[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_NODE] = INTEGER;
  data_length[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_NODE] = DATA_ITEM_SIZE;
  fixed_length[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_NODE] = 0;
  data_class[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_NODE] = GROUNDFLOW;
  data_required[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_NODE] = GROUNDFLOW_PHREATICLEVEL_MULTIPLE;

  strcpy(name[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_STATIC],"groundflow_phreatic_level_multiple_static");
  type[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_STATIC] = INTEGER;
  data_length[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_STATIC] = 1;
  data_class[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_STATIC] = GROUNDFLOW;
  data_required[GROUNDFLOW_PHREATICLEVEL_MULTIPLE_STATIC] = GROUNDFLOW_PHREATICLEVEL_MULTIPLE;

  // Professional groundflow_phreatic_level_static (manual 6.57x, single
  // level): total pressures at the phreatic-level nodes are set to the
  // static pressure (record of the corpus excavate1; consumed by the
  // groundflow family).
  strcpy(name[GROUNDFLOW_PHREATICLEVEL_STATIC],"groundflow_phreatic_level_static");
  type[GROUNDFLOW_PHREATICLEVEL_STATIC] = INTEGER;
  data_length[GROUNDFLOW_PHREATICLEVEL_STATIC] = 1;
  no_index[GROUNDFLOW_PHREATICLEVEL_STATIC] = 1;
  data_class[GROUNDFLOW_PHREATICLEVEL_STATIC] = GROUNDFLOW;
  data_required[GROUNDFLOW_PHREATICLEVEL_STATIC] = GROUNDFLOW_PHREATICLEVEL;

  strcpy(name[GROUNDFLOW_PRESSURE],"groundflow_pressure");

  strcpy(name[GROUNDFLOW_PRESSURE_GRADIENT],"groundflow_pressure_gradient");

  strcpy(name[GROUNDFLOW_PRESSURE_ATMOSPHERIC],"groundflow_pressure_atmospheric");
  type[GROUNDFLOW_PRESSURE_ATMOSPHERIC] = DOUBLE_PRECISION;
  data_length[GROUNDFLOW_PRESSURE_ATMOSPHERIC] = 1;
  no_index[GROUNDFLOW_PRESSURE_ATMOSPHERIC] = 1;
  data_class[GROUNDFLOW_PRESSURE_ATMOSPHERIC] = GROUNDFLOW;

  strcpy(name[GROUNDFLOW_PRESSURE_FACTOR],"groundflow_pressure_factor");
  type[GROUNDFLOW_PRESSURE_FACTOR] = DOUBLE_PRECISION;
  data_length[GROUNDFLOW_PRESSURE_FACTOR] = 1;
  no_index[GROUNDFLOW_PRESSURE_FACTOR] = 1;
  data_class[GROUNDFLOW_PRESSURE_FACTOR] = GROUNDFLOW;

  strcpy(name[GROUNDFLOW_SATURATION],"groundflow_saturation");
  type[GROUNDFLOW_SATURATION] = INTEGER;
  data_length[GROUNDFLOW_SATURATION] = 1;
  no_index[GROUNDFLOW_SATURATION] = 1;
  data_class[GROUNDFLOW_SATURATION] = GROUNDFLOW;

  strcpy(name[GROUNDFLOW_SEEPAGE_EPS],"groundflow_seepage_eps");
  type[GROUNDFLOW_SEEPAGE_EPS] = DOUBLE_PRECISION;
  data_length[GROUNDFLOW_SEEPAGE_EPS] = 1;
  no_index[GROUNDFLOW_SEEPAGE_EPS] = 1;
  data_class[GROUNDFLOW_SEEPAGE_EPS] = GROUNDFLOW;

  strcpy(name[GROUNDFLOW_SEEPAGE_GEOMETRY],"groundflow_seepage_geometry");
  type[GROUNDFLOW_SEEPAGE_GEOMETRY] = INTEGER;
  data_length[GROUNDFLOW_SEEPAGE_GEOMETRY] = 2;
  data_class[GROUNDFLOW_SEEPAGE_GEOMETRY] = GROUNDFLOW;

  strcpy(name[GROUNDFLOW_SEEPAGE_NODE],"groundflow_seepage_node");
  type[GROUNDFLOW_SEEPAGE_NODE] = INTEGER;
  data_length[GROUNDFLOW_SEEPAGE_NODE] = DATA_ITEM_SIZE;
  fixed_length[GROUNDFLOW_SEEPAGE_NODE] = 0;
  data_class[GROUNDFLOW_SEEPAGE_NODE] = GROUNDFLOW;

  strcpy(name[GROUNDFLOW_TOTAL_PRESSURE_LIMIT],"groundflow_total_pressure_limit");
  type[GROUNDFLOW_TOTAL_PRESSURE_LIMIT] = DOUBLE_PRECISION;
  data_length[GROUNDFLOW_TOTAL_PRESSURE_LIMIT] = 1;
  no_index[GROUNDFLOW_TOTAL_PRESSURE_LIMIT] = 1;
  data_class[GROUNDFLOW_TOTAL_PRESSURE_LIMIT] = GROUNDFLOW;

  strcpy(name[GROUNDFLOW_VELOCITY],"groundflow_velocity");

  strcpy(name[GROUP_AXISYMMETRIC],"group_axisymmetric");
  type[GROUP_AXISYMMETRIC] = INTEGER;
  data_length[GROUP_AXISYMMETRIC] = 1;
  data_class[GROUP_AXISYMMETRIC] = GROUP_AXISYMMETRIC;
  data_required[GROUP_AXISYMMETRIC] = GROUP_TYPE;

  strcpy(name[GROUP_BEAM_INERTIA],"group_beam_inertia");
  type[GROUP_BEAM_INERTIA] = DOUBLE_PRECISION;
  // manual Professional 6.590: index Iyy Izz J (3 values). The GNU beam
  // (2D x-y bending about the local z axis) consumes Izz = value 2; the
  // legacy GNU single-value form (1 value, used as Izz) is still read.
  data_length[GROUP_BEAM_INERTIA] = 3;
  fixed_length[GROUP_BEAM_INERTIA] = 0;
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

  strcpy(name[GROUP_CONTACTSPRING_DIRECTION_AUTOMATIC],"group_contactspring_direction_automatic");
  type[GROUP_CONTACTSPRING_DIRECTION_AUTOMATIC] = INTEGER;
  data_length[GROUP_CONTACTSPRING_DIRECTION_AUTOMATIC] = 1;
  data_class[GROUP_CONTACTSPRING_DIRECTION_AUTOMATIC] = CONTACTSPRING;
  data_required[GROUP_CONTACTSPRING_DIRECTION_AUTOMATIC] = GROUP_TYPE;

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

  strcpy(name[GROUP_GROUNDFLOW_CONSOLIDATION_APPLY],"group_groundflow_consolidation_apply");
  type[GROUP_GROUNDFLOW_CONSOLIDATION_APPLY] = INTEGER;
  data_length[GROUP_GROUNDFLOW_CONSOLIDATION_APPLY] = 1;
  data_class[GROUP_GROUNDFLOW_CONSOLIDATION_APPLY] = GROUNDFLOW;
  data_required[GROUP_GROUNDFLOW_CONSOLIDATION_APPLY] = GROUP_TYPE;

  strcpy(name[GROUP_GROUNDFLOW_MATERIDIVERGENCE],"group_groundflow_materidivergence");
  type[GROUP_GROUNDFLOW_MATERIDIVERGENCE] = INTEGER;
  data_length[GROUP_GROUNDFLOW_MATERIDIVERGENCE] = 1;
  data_class[GROUP_GROUNDFLOW_MATERIDIVERGENCE] = GROUNDFLOW;
  data_required[GROUP_GROUNDFLOW_MATERIDIVERGENCE] = GROUP_TYPE;

  strcpy(name[GROUP_GROUNDFLOW_NONSATURATED_EPS_PERMEABILITY],"group_groundflow_nonsaturated_eps_permeability");
  type[GROUP_GROUNDFLOW_NONSATURATED_EPS_PERMEABILITY] = DOUBLE_PRECISION;
  data_length[GROUP_GROUNDFLOW_NONSATURATED_EPS_PERMEABILITY] = 1;
  data_class[GROUP_GROUNDFLOW_NONSATURATED_EPS_PERMEABILITY] = GROUNDFLOW;
  data_required[GROUP_GROUNDFLOW_NONSATURATED_EPS_PERMEABILITY] = GROUP_TYPE;

  strcpy(name[GROUP_GROUNDFLOW_NONSATURATED_VANGENUCHTEN],"group_groundflow_nonsaturated_vangenuchten");
  type[GROUP_GROUNDFLOW_NONSATURATED_VANGENUCHTEN] = DOUBLE_PRECISION;
  data_length[GROUP_GROUNDFLOW_NONSATURATED_VANGENUCHTEN] = 5;
  data_class[GROUP_GROUNDFLOW_NONSATURATED_VANGENUCHTEN] = GROUNDFLOW;
  data_required[GROUP_GROUNDFLOW_NONSATURATED_VANGENUCHTEN] = GROUP_TYPE;

  strcpy(name[GROUP_GROUNDFLOW_PERMEABILITY],"group_groundflow_permeability");
  type[GROUP_GROUNDFLOW_PERMEABILITY] = DOUBLE_PRECISION;
  data_length[GROUP_GROUNDFLOW_PERMEABILITY] = ndim;
  data_class[GROUP_GROUNDFLOW_PERMEABILITY] = GROUNDFLOW;
  data_required[GROUP_GROUNDFLOW_PERMEABILITY] = GROUP_TYPE;
  // variable length 1..ndim: a single value is used in each space
  // direction (manual Professional 6.618)
  fixed_length[GROUP_GROUNDFLOW_PERMEABILITY] = 0;

  strcpy(name[GROUP_GROUNDFLOW_PERMEABILITY_NONLINEAR_METHOD],"group_groundflow_permeability_nonlinear_method");
  type[GROUP_GROUNDFLOW_PERMEABILITY_NONLINEAR_METHOD] = INTEGER;
  data_length[GROUP_GROUNDFLOW_PERMEABILITY_NONLINEAR_METHOD] = 1;
  data_class[GROUP_GROUNDFLOW_PERMEABILITY_NONLINEAR_METHOD] = GROUNDFLOW;
  data_required[GROUP_GROUNDFLOW_PERMEABILITY_NONLINEAR_METHOD] = GROUP_TYPE;

  strcpy(name[GROUP_GROUNDFLOW_PERMEABILITY_NONLINEAR_PARAMETERS],"group_groundflow_permeability_nonlinear_parameters");
  type[GROUP_GROUNDFLOW_PERMEABILITY_NONLINEAR_PARAMETERS] = DOUBLE_PRECISION;
  data_length[GROUP_GROUNDFLOW_PERMEABILITY_NONLINEAR_PARAMETERS] = 1;
  data_class[GROUP_GROUNDFLOW_PERMEABILITY_NONLINEAR_PARAMETERS] = GROUNDFLOW;
  data_required[GROUP_GROUNDFLOW_PERMEABILITY_NONLINEAR_PARAMETERS] = GROUP_TYPE;

  strcpy(name[GROUP_GROUNDFLOW_PERMEABILITY_VERTICAL_STRESS],"group_groundflow_permeability_vertical_stress");
  type[GROUP_GROUNDFLOW_PERMEABILITY_VERTICAL_STRESS] = DOUBLE_PRECISION;
  data_length[GROUP_GROUNDFLOW_PERMEABILITY_VERTICAL_STRESS] = 5;
  data_class[GROUP_GROUNDFLOW_PERMEABILITY_VERTICAL_STRESS] = GROUNDFLOW;
  data_required[GROUP_GROUNDFLOW_PERMEABILITY_VERTICAL_STRESS] = GROUP_TYPE;
  fixed_length[GROUP_GROUNDFLOW_PERMEABILITY_NONLINEAR_PARAMETERS] = 0;

  strcpy(name[GROUP_GROUNDFLOW_POROSITY],"group_groundflow_porosity");
  type[GROUP_GROUNDFLOW_POROSITY] = DOUBLE_PRECISION;
  data_length[GROUP_GROUNDFLOW_POROSITY] = 1;
  data_class[GROUP_GROUNDFLOW_POROSITY] = GROUNDFLOW;
  data_required[GROUP_GROUNDFLOW_POROSITY] = GROUP_TYPE;

  strcpy(name[GROUP_GROUNDFLOW_TOTAL_PRESSURE_TENSION],"group_groundflow_total_pressure_tension");
  type[GROUP_GROUNDFLOW_TOTAL_PRESSURE_TENSION] = DOUBLE_PRECISION;
  data_length[GROUP_GROUNDFLOW_TOTAL_PRESSURE_TENSION] = 2;
  data_class[GROUP_GROUNDFLOW_TOTAL_PRESSURE_TENSION] = GROUNDFLOW;
  data_required[GROUP_GROUNDFLOW_TOTAL_PRESSURE_TENSION] = GROUP_TYPE;

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

  strcpy(name[GROUP_ELEMENT_SELECTIVE_REDUCED_INTEGRATION],
    "group_element_selective_reduced_integration");
  type[GROUP_ELEMENT_SELECTIVE_REDUCED_INTEGRATION] = INTEGER;
  data_length[GROUP_ELEMENT_SELECTIVE_REDUCED_INTEGRATION] = 1;
  version_all[GROUP_ELEMENT_SELECTIVE_REDUCED_INTEGRATION] = 1;
  data_class[GROUP_ELEMENT_SELECTIVE_REDUCED_INTEGRATION] = GROUP_TYPE;
  data_required[GROUP_ELEMENT_SELECTIVE_REDUCED_INTEGRATION] = GROUP_TYPE;

  strcpy(name[GROUP_INTERFACE],"group_interface");
  type[GROUP_INTERFACE] = INTEGER;
  data_length[GROUP_INTERFACE] = 1;
  data_class[GROUP_INTERFACE] = GROUP_TYPE;
  data_required[GROUP_INTERFACE] = GROUP_TYPE;

  strcpy(name[GROUP_INTERFACE_CONDIF_CONDUCTIVITY],"group_interface_condif_conductivity");
  type[GROUP_INTERFACE_CONDIF_CONDUCTIVITY] = DOUBLE_PRECISION;
  data_length[GROUP_INTERFACE_CONDIF_CONDUCTIVITY] = 1;
  data_class[GROUP_INTERFACE_CONDIF_CONDUCTIVITY] = GROUP_INTERFACE;
  data_required[GROUP_INTERFACE_CONDIF_CONDUCTIVITY] = GROUP_INTERFACE;

  strcpy(name[GROUP_INTERFACE_MATERI_EXPANSION_NORMAL],"group_interface_materi_expansion_normal");
  type[GROUP_INTERFACE_MATERI_EXPANSION_NORMAL] = DOUBLE_PRECISION;
  data_length[GROUP_INTERFACE_MATERI_EXPANSION_NORMAL] = 1;
  data_class[GROUP_INTERFACE_MATERI_EXPANSION_NORMAL] = GROUP_INTERFACE;
  data_required[GROUP_INTERFACE_MATERI_EXPANSION_NORMAL] = GROUP_INTERFACE;

  strcpy(name[GROUP_INTERFACE_TANGENTIAL_REFERENCE_POINT],"group_interface_tangential_reference_point");
  type[GROUP_INTERFACE_TANGENTIAL_REFERENCE_POINT] = DOUBLE_PRECISION;
  data_length[GROUP_INTERFACE_TANGENTIAL_REFERENCE_POINT] = MDIM;
  data_class[GROUP_INTERFACE_TANGENTIAL_REFERENCE_POINT] = GROUP_INTERFACE;
  data_required[GROUP_INTERFACE_TANGENTIAL_REFERENCE_POINT] = GROUP_INTERFACE;


  // Sprint 13: the interface element print items + damping (the
  // Professional's interface tests)
  strcpy(name[ELEMENT_INTERFACE_STRESS_AVERAGE],"element_interface_stress_average");
  type[ELEMENT_INTERFACE_STRESS_AVERAGE] = DOUBLE_PRECISION;
  data_length[ELEMENT_INTERFACE_STRESS_AVERAGE] = 3;
  version_all[ELEMENT_INTERFACE_STRESS_AVERAGE] = 1;
  print_only[ELEMENT_INTERFACE_STRESS_AVERAGE] = 1;
  data_class[ELEMENT_INTERFACE_STRESS_AVERAGE] = ELEMENT;
  data_required[ELEMENT_INTERFACE_STRESS_AVERAGE] = ELEMENT;
  fixed_length[ELEMENT_INTERFACE_STRESS_AVERAGE] = 0;

  strcpy(name[ELEMENT_INTERFACE_INTPNT_STRESS],"element_interface_intpnt_stress");
  type[ELEMENT_INTERFACE_INTPNT_STRESS] = DOUBLE_PRECISION;
  data_length[ELEMENT_INTERFACE_INTPNT_STRESS] = 12;
  version_all[ELEMENT_INTERFACE_INTPNT_STRESS] = 1;
  print_only[ELEMENT_INTERFACE_INTPNT_STRESS] = 1;
  data_class[ELEMENT_INTERFACE_INTPNT_STRESS] = ELEMENT;
  data_required[ELEMENT_INTERFACE_INTPNT_STRESS] = ELEMENT;
  fixed_length[ELEMENT_INTERFACE_INTPNT_STRESS] = 0;

  strcpy(name[ELEMENT_INTERFACE_STRAIN_AVERAGE],"element_interface_strain_average");
  type[ELEMENT_INTERFACE_STRAIN_AVERAGE] = DOUBLE_PRECISION;
  data_length[ELEMENT_INTERFACE_STRAIN_AVERAGE] = 3;
  version_all[ELEMENT_INTERFACE_STRAIN_AVERAGE] = 1;
  print_only[ELEMENT_INTERFACE_STRAIN_AVERAGE] = 1;
  data_class[ELEMENT_INTERFACE_STRAIN_AVERAGE] = ELEMENT;
  data_required[ELEMENT_INTERFACE_STRAIN_AVERAGE] = ELEMENT;
  fixed_length[ELEMENT_INTERFACE_STRAIN_AVERAGE] = 0;

  strcpy(name[ELEMENT_INTERFACE_INTPNT_MATERI_TENSION_STATUS],"element_interface_intpnt_materi_tension_status");
  type[ELEMENT_INTERFACE_INTPNT_MATERI_TENSION_STATUS] = INTEGER;
  data_length[ELEMENT_INTERFACE_INTPNT_MATERI_TENSION_STATUS] = 4;
  version_all[ELEMENT_INTERFACE_INTPNT_MATERI_TENSION_STATUS] = 1;
  print_only[ELEMENT_INTERFACE_INTPNT_MATERI_TENSION_STATUS] = 1;
  data_class[ELEMENT_INTERFACE_INTPNT_MATERI_TENSION_STATUS] = ELEMENT;
  data_required[ELEMENT_INTERFACE_INTPNT_MATERI_TENSION_STATUS] = ELEMENT;
  fixed_length[ELEMENT_INTERFACE_INTPNT_MATERI_TENSION_STATUS] = 0;

  strcpy(name[ELEMENT_INTERFACE_INTPNT_STRAIN],"element_interface_intpnt_strain");
  type[ELEMENT_INTERFACE_INTPNT_STRAIN] = DOUBLE_PRECISION;
  data_length[ELEMENT_INTERFACE_INTPNT_STRAIN] = 12;
  version_all[ELEMENT_INTERFACE_INTPNT_STRAIN] = 1;
  print_only[ELEMENT_INTERFACE_INTPNT_STRAIN] = 1;
  data_class[ELEMENT_INTERFACE_INTPNT_STRAIN] = ELEMENT;
  data_required[ELEMENT_INTERFACE_INTPNT_STRAIN] = ELEMENT;
  fixed_length[ELEMENT_INTERFACE_INTPNT_STRAIN] = 0;

  strcpy(name[GROUP_INTERFACE_DAMPING],"group_interface_damping");
  type[GROUP_INTERFACE_DAMPING] = DOUBLE_PRECISION;
  data_length[GROUP_INTERFACE_DAMPING] = 1;
  data_class[GROUP_INTERFACE_DAMPING] = GROUP_INTERFACE;
  strcpy(name[GROUP_INTERFACE_MATERI_PLASTI_MOHR_COUL_DIRECT],"group_interface_materi_plasti_mohr_coul_direct");

  strcpy(name[GROUP_INTERFACE_GAP],"group_interface_gap");
  type[GROUP_INTERFACE_GAP] = DOUBLE_PRECISION;
  data_length[GROUP_INTERFACE_GAP] = 1;
  data_class[GROUP_INTERFACE_GAP] = GROUP_TYPE;
  data_required[GROUP_INTERFACE_GAP] = GROUP_INTERFACE;

  strcpy(name[GROUP_INTERFACE_GROUNDFLOW_CAPACITY],"group_interface_groundflow_capacity");
  type[GROUP_INTERFACE_GROUNDFLOW_CAPACITY] = DOUBLE_PRECISION;
  data_length[GROUP_INTERFACE_GROUNDFLOW_CAPACITY] = 1;
  data_class[GROUP_INTERFACE_GROUNDFLOW_CAPACITY] = GROUNDFLOW;
  data_required[GROUP_INTERFACE_GROUNDFLOW_CAPACITY] = GROUP_INTERFACE;

  strcpy(name[GROUP_INTERFACE_GROUNDFLOW_PERMEABILITY],"group_interface_groundflow_permeability");
  type[GROUP_INTERFACE_GROUNDFLOW_PERMEABILITY] = DOUBLE_PRECISION;
  data_length[GROUP_INTERFACE_GROUNDFLOW_PERMEABILITY] = 1;
  data_class[GROUP_INTERFACE_GROUNDFLOW_PERMEABILITY] = GROUNDFLOW;
  data_required[GROUP_INTERFACE_GROUNDFLOW_PERMEABILITY] = GROUP_INTERFACE;

  strcpy(name[GROUP_INTERFACE_GROUNDFLOW_TOTAL_PRESSURE_TENSION],"group_interface_groundflow_total_pressure_tension");
  type[GROUP_INTERFACE_GROUNDFLOW_TOTAL_PRESSURE_TENSION] = DOUBLE_PRECISION;
  data_length[GROUP_INTERFACE_GROUNDFLOW_TOTAL_PRESSURE_TENSION] = 2;
  data_class[GROUP_INTERFACE_GROUNDFLOW_TOTAL_PRESSURE_TENSION] = GROUNDFLOW;
  data_required[GROUP_INTERFACE_GROUNDFLOW_TOTAL_PRESSURE_TENSION] = GROUP_INTERFACE;

  strcpy(name[GROUP_INTERFACE_MATERI_ELASTI_STIFFNESS],"group_interface_materi_elasti_stiffness");
  type[GROUP_INTERFACE_MATERI_ELASTI_STIFFNESS] = DOUBLE_PRECISION;
  data_length[GROUP_INTERFACE_MATERI_ELASTI_STIFFNESS] = 3;
  // variable-length: 2 values in 2D (kn, kt), 3 in 3D (kn, kt1, kt2)
  data_length[GROUP_INTERFACE_MATERI_ELASTI_STIFFNESS] = DATA_ITEM_SIZE;
  fixed_length[GROUP_INTERFACE_MATERI_ELASTI_STIFFNESS] = 0;
  data_class[GROUP_INTERFACE_MATERI_ELASTI_STIFFNESS] = GROUP_TYPE;
  data_required[GROUP_INTERFACE_MATERI_ELASTI_STIFFNESS] = GROUP_INTERFACE;

  strcpy(name[GROUP_INTERFACE_MATERI_MEMORY],"group_interface_materi_memory");
  type[GROUP_INTERFACE_MATERI_MEMORY] = INTEGER;
  data_length[GROUP_INTERFACE_MATERI_MEMORY] = 1;
  data_class[GROUP_INTERFACE_MATERI_MEMORY] = GROUP_TYPE;
  data_required[GROUP_INTERFACE_MATERI_MEMORY] = GROUP_INTERFACE;

  strcpy(name[GROUP_INTERFACE_MATERI_PLASTI_MOHR_COUL_DIRECT],"group_interface_materi_plasti_mohr_coul_direct");
  type[GROUP_INTERFACE_MATERI_PLASTI_MOHR_COUL_DIRECT] = DOUBLE_PRECISION;
  data_length[GROUP_INTERFACE_MATERI_PLASTI_MOHR_COUL_DIRECT] = 3;
  data_class[GROUP_INTERFACE_MATERI_PLASTI_MOHR_COUL_DIRECT] = GROUP_TYPE;
  data_required[GROUP_INTERFACE_MATERI_PLASTI_MOHR_COUL_DIRECT] = GROUP_INTERFACE;

  strcpy(name[GROUP_INTERFACE_MATERI_PLASTI_TENSION_DIRECT],"group_interface_materi_plasti_tension_direct");
  type[GROUP_INTERFACE_MATERI_PLASTI_TENSION_DIRECT] = DOUBLE_PRECISION;
  data_length[GROUP_INTERFACE_MATERI_PLASTI_TENSION_DIRECT] = 1;
  data_class[GROUP_INTERFACE_MATERI_PLASTI_TENSION_DIRECT] = GROUP_TYPE;
  data_required[GROUP_INTERFACE_MATERI_PLASTI_TENSION_DIRECT] = GROUP_INTERFACE;

  strcpy(name[GROUP_INTERFACE_MATERI_RESIDUAL_STIFFNESS],"group_interface_materi_residual_stiffness");
  type[GROUP_INTERFACE_MATERI_RESIDUAL_STIFFNESS] = DOUBLE_PRECISION;
  data_length[GROUP_INTERFACE_MATERI_RESIDUAL_STIFFNESS] = 1;
  data_class[GROUP_INTERFACE_MATERI_RESIDUAL_STIFFNESS] = GROUP_TYPE;
  data_required[GROUP_INTERFACE_MATERI_RESIDUAL_STIFFNESS] = GROUP_INTERFACE;

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

  strcpy(name[GROUP_MATERI_ELASTI_CAMCLAY_PRESSURE_MIN],"group_materi_elasti_camclay_pressure_min");
  type[GROUP_MATERI_ELASTI_CAMCLAY_PRESSURE_MIN] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_ELASTI_CAMCLAY_PRESSURE_MIN] = 1;
  data_class[GROUP_MATERI_ELASTI_CAMCLAY_PRESSURE_MIN] = MATERI;
  data_required[GROUP_MATERI_ELASTI_CAMCLAY_PRESSURE_MIN] = GROUP_TYPE;

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

  strcpy(name[GROUP_MATERI_DAMPING_METHOD],"group_materi_damping_method");
  type[GROUP_MATERI_DAMPING_METHOD] = INTEGER;
  data_length[GROUP_MATERI_DAMPING_METHOD] = 1;
  data_required[GROUP_MATERI_DAMPING_METHOD] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_FACTOR],"group_materi_factor");
  type[GROUP_MATERI_FACTOR] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_FACTOR] = 1;
  data_required[GROUP_MATERI_FACTOR] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL_LIMIT],"group_materi_plasti_visco_exponential_limit");
  type[GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL_LIMIT] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL_LIMIT] = 1;
  data_required[GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL_LIMIT] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL_NAME],"group_materi_plasti_visco_exponential_name");
  type[GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL_NAME] = INTEGER;
  data_length[GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL_NAME] = 1;
  data_required[GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL_NAME] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL_VALUES],"group_materi_plasti_visco_exponential_values");
  type[GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL_VALUES] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL_VALUES] = DATA_ITEM_SIZE;
  fixed_length[GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL_VALUES] = 0;
  data_required[GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL_VALUES] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_ELASTI_COMPRESSIBILITY],"group_materi_elasti_compressibility");
  type[GROUP_MATERI_ELASTI_COMPRESSIBILITY] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_ELASTI_COMPRESSIBILITY] = 1;
  data_class[GROUP_MATERI_ELASTI_COMPRESSIBILITY] = MATERI;
  data_required[GROUP_MATERI_ELASTI_COMPRESSIBILITY] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_ELASTI_K0],"group_materi_elasti_k0");
  type[GROUP_MATERI_ELASTI_K0] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_ELASTI_K0] = 1;
  data_class[GROUP_MATERI_ELASTI_K0] = MATERI;
  data_required[GROUP_MATERI_ELASTI_K0] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_ELASTI_HARDSOIL],"group_materi_elasti_hardsoil");
  type[GROUP_MATERI_ELASTI_HARDSOIL] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_ELASTI_HARDSOIL] = 7;
  data_class[GROUP_MATERI_ELASTI_HARDSOIL] = MATERI;
  data_required[GROUP_MATERI_ELASTI_HARDSOIL] = GROUP_TYPE;

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

  strcpy(name[GROUP_MATERI_ELASTI_POISSON_POWER],"group_materi_elasti_poisson_power");
  type[GROUP_MATERI_ELASTI_POISSON_POWER] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_ELASTI_POISSON_POWER] = 5;
  data_class[GROUP_MATERI_ELASTI_POISSON_POWER] = MATERI;
  data_required[GROUP_MATERI_ELASTI_POISSON_POWER] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_ELASTI_SHEAR_FACTOR],"group_materi_elasti_shear_factor");
  type[GROUP_MATERI_ELASTI_SHEAR_FACTOR] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_ELASTI_SHEAR_FACTOR] = 1;
  data_class[GROUP_MATERI_ELASTI_SHEAR_FACTOR] = MATERI;
  data_required[GROUP_MATERI_ELASTI_SHEAR_FACTOR] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_ELASTI_SMALLSTRAIN],"group_materi_elasti_smallstrain");
  type[GROUP_MATERI_ELASTI_SMALLSTRAIN] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_ELASTI_SMALLSTRAIN] = 6;
  data_class[GROUP_MATERI_ELASTI_SMALLSTRAIN] = MATERI;
  data_required[GROUP_MATERI_ELASTI_SMALLSTRAIN] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_ELASTI_STRESS_PRESSURE_HISTORY_FACTOR],"group_materi_elasti_stress_pressure_history_factor");
  type[GROUP_MATERI_ELASTI_STRESS_PRESSURE_HISTORY_FACTOR] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_ELASTI_STRESS_PRESSURE_HISTORY_FACTOR] = 1;
  data_class[GROUP_MATERI_ELASTI_STRESS_PRESSURE_HISTORY_FACTOR] = MATERI;
  data_required[GROUP_MATERI_ELASTI_STRESS_PRESSURE_HISTORY_FACTOR] = GROUP_TYPE;

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
  data_length[GROUP_MATERI_ELASTI_YOUNG_POWER] = 6;
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

  strcpy(name[GROUP_MATERI_FAILURE_CRUCHING],"group_materi_failure_crunching");
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

  strcpy(name[GROUP_MATERI_FAILURE_VOIDFRACTION],"group_materi_failure_void_fraction");
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
  // manual Professional 6.689: index M kappa lambda (3 values). The
  // preconsolidation pressure p0 is a HISTORY variable
  // (materi_plasti_camclay_history / node_dof cchis0/cchis1), NOT a
  // fourth material parameter (the GNU legacy 4th value N is derived
  // from the initial state instead).
  data_length[GROUP_MATERI_PLASTI_CAMCLAY] = 3;
  fixed_length[GROUP_MATERI_PLASTI_CAMCLAY] = 0;
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

  strcpy(name[GROUP_MATERI_PLASTI_CAP1],"group_materi_plasti_cap1");
  type[GROUP_MATERI_PLASTI_CAP1] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_CAP1] = 8;
  data_class[GROUP_MATERI_PLASTI_CAP1] = MATERI;
  data_required[GROUP_MATERI_PLASTI_CAP1] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_HARDSOIL],"group_materi_plasti_hardsoil");
  type[GROUP_MATERI_PLASTI_HARDSOIL] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_HARDSOIL] = 4;
  data_class[GROUP_MATERI_PLASTI_HARDSOIL] = MATERI;
  data_required[GROUP_MATERI_PLASTI_HARDSOIL] = GROUP_TYPE;

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

  strcpy(name[GROUP_MATERI_PLASTI_MOHR_COUL_HARDENING_SOFTENING],"group_materi_plasti_mohr_coul_hardening_softening");
  type[GROUP_MATERI_PLASTI_MOHR_COUL_HARDENING_SOFTENING] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_MOHR_COUL_HARDENING_SOFTENING] = 7;
  data_class[GROUP_MATERI_PLASTI_MOHR_COUL_HARDENING_SOFTENING] = MATERI;
  data_required[GROUP_MATERI_PLASTI_MOHR_COUL_HARDENING_SOFTENING] = GROUP_TYPE;

  // manual Professional (dam_building): the "direct" variant of the
  // hardening-softening record - same 7 values (phi_0 c_0 phiflow_0
  // phi_1 c_1 phiflow_1 kappashear_crit) but the angles are in DEGREES
  // like every group_materi_plasti_mohr_coul_direct record.
  strcpy(name[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_HARDENING_SOFTENING],"group_materi_plasti_mohr_coul_direct_hardening_softening");
  type[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_HARDENING_SOFTENING] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_HARDENING_SOFTENING] = 7;
  data_class[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_HARDENING_SOFTENING] = MATERI;
  data_required[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_HARDENING_SOFTENING] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_MOHR_COUL],"group_materi_plasti_mohr_coul");
  type[GROUP_MATERI_PLASTI_MOHR_COUL] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_MOHR_COUL] = 3;
  data_class[GROUP_MATERI_PLASTI_MOHR_COUL] = MATERI;
  data_required[GROUP_MATERI_PLASTI_MOHR_COUL] = GROUP_TYPE;

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

  strcpy(name[GROUP_MATERI_PLASTI_HEAT_GENERATION],"group_materi_plasti_heat_generation");
  type[GROUP_MATERI_PLASTI_HEAT_GENERATION] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_HEAT_GENERATION] = 1;
  data_class[GROUP_MATERI_PLASTI_HEAT_GENERATION] = MATERI;
  data_required[GROUP_MATERI_PLASTI_HEAT_GENERATION] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_COMPRESSION_DIRECT],"group_materi_plasti_compression_direct");
  type[GROUP_MATERI_PLASTI_COMPRESSION_DIRECT] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_COMPRESSION_DIRECT] = 1;
  data_class[GROUP_MATERI_PLASTI_COMPRESSION_DIRECT] = MATERI;
  data_required[GROUP_MATERI_PLASTI_COMPRESSION_DIRECT] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_COMPRESSION_DIRECT_VISCO],"group_materi_plasti_compression_direct_visco");
  type[GROUP_MATERI_PLASTI_COMPRESSION_DIRECT_VISCO] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_COMPRESSION_DIRECT_VISCO] = 1;
  data_class[GROUP_MATERI_PLASTI_COMPRESSION_DIRECT_VISCO] = MATERI;
  data_required[GROUP_MATERI_PLASTI_COMPRESSION_DIRECT_VISCO] = GROUP_MATERI_PLASTI_COMPRESSION_DIRECT;

  strcpy(name[GROUP_MATERI_PLASTI_PRESSURE_LIMIT],"group_materi_plasti_pressure_limit");
  type[GROUP_MATERI_PLASTI_PRESSURE_LIMIT] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_PRESSURE_LIMIT] = 1;
  data_class[GROUP_MATERI_PLASTI_PRESSURE_LIMIT] = MATERI;
  data_required[GROUP_MATERI_PLASTI_PRESSURE_LIMIT] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_COORD_LIMIT],"group_materi_plasti_coord_limit");
  type[GROUP_MATERI_PLASTI_COORD_LIMIT] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_COORD_LIMIT] = 1;
  data_class[GROUP_MATERI_PLASTI_COORD_LIMIT] = MATERI;
  data_required[GROUP_MATERI_PLASTI_COORD_LIMIT] = GROUP_TYPE;

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
  // manual Professional 6.715: R mR mT beta_r chi theta (6 params).
  // theta is the exponent of the rho^theta f_d N S_hat term in the
  // stiffness; for monotonic loading the manual recommends theta=chi,
  // and the GNU kernel evaluates that term with chi, so both coincide.
  data_length[GROUP_MATERI_PLASTI_HYPO_INTERGRANULARSTRAIN] = 6;
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

  strcpy(name[GROUP_MATERI_PLASTI_HYPO_MASIN],"group_materi_plasti_hypo_masin");
  type[GROUP_MATERI_PLASTI_HYPO_MASIN] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_HYPO_MASIN] = 5;
  data_class[GROUP_MATERI_PLASTI_HYPO_MASIN] = MATERI;
  data_required[GROUP_MATERI_PLASTI_HYPO_MASIN] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_HYPO_MASIN_STRUCTURE],"group_materi_plasti_hypo_masin_structure");
  type[GROUP_MATERI_PLASTI_HYPO_MASIN_STRUCTURE] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_HYPO_MASIN_STRUCTURE] = 3;
  data_class[GROUP_MATERI_PLASTI_HYPO_MASIN_STRUCTURE] = MATERI;
  data_required[GROUP_MATERI_PLASTI_HYPO_MASIN_STRUCTURE] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_HYPO_MASIN_OCR],"group_materi_plasti_hypo_masin_ocr");
  type[GROUP_MATERI_PLASTI_HYPO_MASIN_OCR] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_HYPO_MASIN_OCR] = 1;
  data_class[GROUP_MATERI_PLASTI_HYPO_MASIN_OCR] = MATERI;
  data_required[GROUP_MATERI_PLASTI_HYPO_MASIN_OCR] = GROUP_TYPE;

  strcpy(name[CONTROL_MATERI_DAMAGE_APPLY],"control_materi_damage_apply");
  type[CONTROL_MATERI_DAMAGE_APPLY] = INTEGER;
  data_length[CONTROL_MATERI_DAMAGE_APPLY] = 1;
  data_class[CONTROL_MATERI_DAMAGE_APPLY] = CONTROL;

  strcpy(name[CONTROL_MATERI_DYNAMIC],"control_materi_dynamic");
  type[CONTROL_MATERI_DYNAMIC] = DOUBLE_PRECISION;
  data_length[CONTROL_MATERI_DYNAMIC] = 1;
  data_class[CONTROL_MATERI_DYNAMIC] = CONTROL;

  strcpy(name[MATERI_DYNAMIC],"materi_dynamic");
  type[MATERI_DYNAMIC] = DOUBLE_PRECISION;
  data_length[MATERI_DYNAMIC] = 1;
  data_class[MATERI_DYNAMIC] = CONTROL;
  no_index[MATERI_DYNAMIC] = 1;

  strcpy(name[CONTROL_MATERI_ELASTI_K0],"control_materi_elasti_k0");
  type[CONTROL_MATERI_ELASTI_K0] = INTEGER;
  data_length[CONTROL_MATERI_ELASTI_K0] = 1;
  data_class[CONTROL_MATERI_ELASTI_K0] = CONTROL;

  strcpy(name[CONTROL_MATERI_ELASTI_YOUNG_POWER_APPLY],"control_materi_elasti_young_power_apply");
  type[CONTROL_MATERI_ELASTI_YOUNG_POWER_APPLY] = INTEGER;
  data_length[CONTROL_MATERI_ELASTI_YOUNG_POWER_APPLY] = 1;
  data_class[CONTROL_MATERI_ELASTI_YOUNG_POWER_APPLY] = CONTROL;

  strcpy(name[CONTROL_MATERI_FAILURE_APPLY],"control_materi_failure_apply");
  type[CONTROL_MATERI_FAILURE_APPLY] = INTEGER;
  data_length[CONTROL_MATERI_FAILURE_APPLY] = 1;
  data_class[CONTROL_MATERI_FAILURE_APPLY] = CONTROL;

  strcpy(name[CONTROL_MATERI_PLASTI_HARDSOIL_GAMMAP_INITIAL],"control_materi_plasti_hardsoil_gammap_initial");
  type[CONTROL_MATERI_PLASTI_HARDSOIL_GAMMAP_INITIAL] = INTEGER;
  data_length[CONTROL_MATERI_PLASTI_HARDSOIL_GAMMAP_INITIAL] = 1;
  data_class[CONTROL_MATERI_PLASTI_HARDSOIL_GAMMAP_INITIAL] = CONTROL;

  strcpy(name[CONTROL_MATERI_PLASTI_HYPO_NIEMUNIS_VISCO_OCR_APPLY],"control_materi_plasti_hypo_niemunis_visco_ocr_apply");
  type[CONTROL_MATERI_PLASTI_HYPO_NIEMUNIS_VISCO_OCR_APPLY] = INTEGER;
  data_length[CONTROL_MATERI_PLASTI_HYPO_NIEMUNIS_VISCO_OCR_APPLY] = 1;
  data_class[CONTROL_MATERI_PLASTI_HYPO_NIEMUNIS_VISCO_OCR_APPLY] = CONTROL;

  strcpy(name[CONTROL_MATERI_PLASTI_HYPO_PRESSURE_DEPENDENT_VOID_RATIO],"control_materi_plasti_hypo_pressure_dependent_void_ratio");
  type[CONTROL_MATERI_PLASTI_HYPO_PRESSURE_DEPENDENT_VOID_RATIO] = INTEGER;
  data_length[CONTROL_MATERI_PLASTI_HYPO_PRESSURE_DEPENDENT_VOID_RATIO] = 1;
  data_class[CONTROL_MATERI_PLASTI_HYPO_PRESSURE_DEPENDENT_VOID_RATIO] = CONTROL;

  strcpy(name[CONTROL_MATERI_PLASTI_HYPO_SUBSTEPPING],"control_materi_plasti_hypo_substepping");
  type[CONTROL_MATERI_PLASTI_HYPO_SUBSTEPPING] = INTEGER;
  data_length[CONTROL_MATERI_PLASTI_HYPO_SUBSTEPPING] = 1;
  data_class[CONTROL_MATERI_PLASTI_HYPO_SUBSTEPPING] = CONTROL;

  strcpy(name[CONTROL_MATERI_PLASTI_TENSION_APPLY],"control_materi_plasti_tension_apply");
  type[CONTROL_MATERI_PLASTI_TENSION_APPLY] = INTEGER;
  data_length[CONTROL_MATERI_PLASTI_TENSION_APPLY] = 1;
  data_class[CONTROL_MATERI_PLASTI_TENSION_APPLY] = CONTROL;

  strcpy(name[CONTROL_MATERI_PLASTI_VISCO_APPLY],"control_materi_plasti_visco_apply");
  type[CONTROL_MATERI_PLASTI_VISCO_APPLY] = INTEGER;
  data_length[CONTROL_MATERI_PLASTI_VISCO_APPLY] = 1;
  data_class[CONTROL_MATERI_PLASTI_VISCO_APPLY] = CONTROL;

  strcpy(name[CONTROL_MATERI_UNDRAINED_APPLY],"control_materi_undrained_apply");
  type[CONTROL_MATERI_UNDRAINED_APPLY] = INTEGER;
  data_length[CONTROL_MATERI_UNDRAINED_APPLY] = 1;
  data_class[CONTROL_MATERI_UNDRAINED_APPLY] = CONTROL;

  // UNDRAINED-CAPACITY FAMILY (2026-09-05, manual Professional 6.760,
  // 6.441, 6.153 + theory 2.2.7): group_materi_undrained_capacity C
  // models an undrained groundwater analysis without adding the
  // groundwater equation to the system matrix: the total groundwater
  // pressure change of an element follows from C * p_dot =
  // div(v_material), solved on the element level. The pressure is stored
  // per element (one slot per integration point; the average record
  // holds the mean over the integration points). control_materi_
  // undrained_apply (registered above) switches the analysis on/off
  // (default -yes).
  strcpy(name[GROUP_MATERI_UNDRAINED_CAPACITY],"group_materi_undrained_capacity");
  type[GROUP_MATERI_UNDRAINED_CAPACITY] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_UNDRAINED_CAPACITY] = 1;
  data_class[GROUP_MATERI_UNDRAINED_CAPACITY] = MATERI;
  data_required[GROUP_MATERI_UNDRAINED_CAPACITY] = GROUP_TYPE;

  strcpy(name[ELEMENT_INTPNT_MATERI_UNDRAINED_PRESSURE],
    "element_intpnt_materi_undrained_pressure");
  type[ELEMENT_INTPNT_MATERI_UNDRAINED_PRESSURE] = DOUBLE_PRECISION;
  data_length[ELEMENT_INTPNT_MATERI_UNDRAINED_PRESSURE] = npointmax;
  data_class[ELEMENT_INTPNT_MATERI_UNDRAINED_PRESSURE] = ELEMENT;
  data_required[ELEMENT_INTPNT_MATERI_UNDRAINED_PRESSURE] = ELEMENT;
  // version_all like ELEMENT_DOF: the step-end version copy promotes
  // the converged per-step value to VERSION_NORMAL (materi() reads the
  // previous step value from VERSION_NORMAL and stores the current
  // iterate in VERSION_NEW). fixed_length 0: the record holds one slot
  // per integration point of the element (variable record length).
  version_all[ELEMENT_INTPNT_MATERI_UNDRAINED_PRESSURE] = 1;
  fixed_length[ELEMENT_INTPNT_MATERI_UNDRAINED_PRESSURE] = 0;

  strcpy(name[ELEMENT_INTPNT_MATERI_UNDRAINED_PRESSURE_AVERAGE],
    "element_intpnt_materi_undrained_pressure_average");
  type[ELEMENT_INTPNT_MATERI_UNDRAINED_PRESSURE_AVERAGE] = DOUBLE_PRECISION;
  data_length[ELEMENT_INTPNT_MATERI_UNDRAINED_PRESSURE_AVERAGE] = 1;
  data_class[ELEMENT_INTPNT_MATERI_UNDRAINED_PRESSURE_AVERAGE] = ELEMENT;
  data_required[ELEMENT_INTPNT_MATERI_UNDRAINED_PRESSURE_AVERAGE] = ELEMENT;

  // post_calcul operator values -young_apparent/-poisson_apparent
  // (manual Professional 6.903): the apparent Young modulus and Poisson
  // ratio from the INCREMENTAL strains and INCREMENTAL stresses of the
  // last time step (0 when the determination is not possible, e.g.
  // almost zero incremental strains). Resolved like -safety_piping as
  // INTEGER name entries of the post_calcul record.
  strcpy(name[YOUNG_APPARENT],"young_apparent");
  type[YOUNG_APPARENT] = INTEGER;
  data_length[YOUNG_APPARENT] = 1;

  strcpy(name[POISSON_APPARENT],"poisson_apparent");
  type[POISSON_APPARENT] = INTEGER;
  data_length[POISSON_APPARENT] = 1;

  strcpy(name[CONTROL_MATERI_UPDATED_APPLY],"control_materi_updated_apply");
  type[CONTROL_MATERI_UPDATED_APPLY] = INTEGER;
  data_length[CONTROL_MATERI_UPDATED_APPLY] = 1;
  data_class[CONTROL_MATERI_UPDATED_APPLY] = CONTROL;

  strcpy(name[CONTROL_MATERI_VISCOSITY_APPLY],"control_materi_viscosity_apply");
  type[CONTROL_MATERI_VISCOSITY_APPLY] = INTEGER;
  data_length[CONTROL_MATERI_VISCOSITY_APPLY] = 1;
  data_class[CONTROL_MATERI_VISCOSITY_APPLY] = CONTROL;

  strcpy(name[CONTROL_MATERI_PLASTI_HYPO_MASIN_OCR_APPLY],"control_materi_plasti_hypo_masin_ocr_apply");
  type[CONTROL_MATERI_PLASTI_HYPO_MASIN_OCR_APPLY] = INTEGER;
  data_length[CONTROL_MATERI_PLASTI_HYPO_MASIN_OCR_APPLY] = 1;
  data_class[CONTROL_MATERI_PLASTI_HYPO_MASIN_OCR_APPLY] = CONTROL;
  // per-timestep switch indexed by ICONTROL (manual 6.147): no group
  // is required at the same index.

  strcpy(name[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY],"group_materi_plasti_hypo_masin_clay");
  type[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY] = 5;
  data_class[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY] = MATERI;
  data_required[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_ADVANCED_PARAMETERS],"group_materi_plasti_hypo_masin_clay_advanced_parameters");
  type[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_ADVANCED_PARAMETERS] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_ADVANCED_PARAMETERS] = 4;
  data_class[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_ADVANCED_PARAMETERS] = MATERI;
  data_required[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_ADVANCED_PARAMETERS] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_AVANCED_DIRECTION],"group_materi_plasti_hypo_masin_clay_avanced_direction");
  type[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_AVANCED_DIRECTION] = INTEGER;
  data_length[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_AVANCED_DIRECTION] = 1;
  data_class[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_AVANCED_DIRECTION] = MATERI;
  data_required[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_AVANCED_DIRECTION] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_OCR],"group_materi_plasti_hypo_masin_clay_ocr");
  type[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_OCR] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_OCR] = 1;
  data_class[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_OCR] = MATERI;
  data_required[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_OCR] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_STRUCTURE],"group_materi_plasti_hypo_masin_clay_structure");
  type[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_STRUCTURE] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_STRUCTURE] = 3;
  data_class[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_STRUCTURE] = MATERI;
  data_required[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_STRUCTURE] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_VISCO],"group_materi_plasti_hypo_masin_clay_visco");
  type[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_VISCO] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_VISCO] = 2;
  data_class[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_VISCO] = MATERI;
  data_required[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_VISCO] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_VISCO_JM],"group_materi_plasti_hypo_masin_clay_visco_jm");
  type[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_VISCO_JM] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_VISCO_JM] = 5;
  data_class[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_VISCO_JM] = MATERI;
  data_required[GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_VISCO_JM] = GROUP_TYPE;

  strcpy(name[CONTROL_MATERI_PLASTI_HYPO_MASIN_CLAY_OCR_APPLY],"control_materi_plasti_hypo_masin_clay_ocr_apply");
  type[CONTROL_MATERI_PLASTI_HYPO_MASIN_CLAY_OCR_APPLY] = INTEGER;
  data_length[CONTROL_MATERI_PLASTI_HYPO_MASIN_CLAY_OCR_APPLY] = 1;
  data_class[CONTROL_MATERI_PLASTI_HYPO_MASIN_CLAY_OCR_APPLY] = CONTROL;
  // per-timestep switch indexed by ICONTROL (manual 6.147): no group
  // is required at the same index.

  strcpy(name[GROUP_MATERI_PLASTI_HYPO_STRAIN_INTERGRANULAR_MASIN_CLAY],"group_materi_plasti_hypo_strain_intergranular_masin_clay");
  type[GROUP_MATERI_PLASTI_HYPO_STRAIN_INTERGRANULAR_MASIN_CLAY] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_HYPO_STRAIN_INTERGRANULAR_MASIN_CLAY] = 7;
  data_class[GROUP_MATERI_PLASTI_HYPO_STRAIN_INTERGRANULAR_MASIN_CLAY] = MATERI;
  data_required[GROUP_MATERI_PLASTI_HYPO_STRAIN_INTERGRANULAR_MASIN_CLAY] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_SANISAND],"group_materi_plasti_sanisand");
  type[GROUP_MATERI_PLASTI_SANISAND] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_SANISAND] = 19;
  data_class[GROUP_MATERI_PLASTI_SANISAND] = MATERI;
  data_required[GROUP_MATERI_PLASTI_SANISAND] = GROUP_TYPE;

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

  strcpy(name[GROUP_MATERI_PLASTI_MAXIMUM_ITERATIONS],"group_materi_plasti_maximum_iterations");
  type[GROUP_MATERI_PLASTI_MAXIMUM_ITERATIONS] = INTEGER;
  data_length[GROUP_MATERI_PLASTI_MAXIMUM_ITERATIONS] = 1;
  data_class[GROUP_MATERI_PLASTI_MAXIMUM_ITERATIONS] = MATERI;
  data_required[GROUP_MATERI_PLASTI_MAXIMUM_ITERATIONS] = GROUP_TYPE;

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

  strcpy(name[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT],"group_materi_plasti_mohr_coul_direct");
  type[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT] = 3;
  data_class[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT] = MATERI;
  data_required[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_NORMAL],"group_materi_plasti_mohr_coul_direct_normal");
  type[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_NORMAL] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_NORMAL] = 3;
  data_class[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_NORMAL] = MATERI;
  data_required[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_NORMAL] = GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT;
  // the plane normal is ndim-flexible (manual Professional 6.727:
  // "In 1d only specify normal_x, etc."; in 2D nx ny; in 3D nx ny nz)
  fixed_length[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_NORMAL] = 0;

  strcpy(name[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_NORMAL_AUTOMATIC],"group_materi_plasti_mohr_coul_direct_normal_automatic");
  type[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_NORMAL_AUTOMATIC] = INTEGER;
  data_length[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_NORMAL_AUTOMATIC] = 1;
  data_class[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_NORMAL_AUTOMATIC] = MATERI;
  data_required[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_NORMAL_AUTOMATIC] = GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT;

  strcpy(name[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_VISCO],"group_materi_plasti_mohr_coul_direct_visco");
  type[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_VISCO] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_VISCO] = 1;
  data_class[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_VISCO] = MATERI;
  data_required[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_VISCO] = GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT;

  strcpy(name[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_WALL],"group_materi_plasti_mohr_coul_direct_wall");
  type[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_WALL] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_WALL] = 3;
  data_class[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_WALL] = MATERI;
  data_required[GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_WALL] = GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT;

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

  strcpy(name[GROUP_MATERI_PLASTI_TENSION_DIRECT],"group_materi_plasti_tension_direct");
  type[GROUP_MATERI_PLASTI_TENSION_DIRECT] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_TENSION_DIRECT] = 1;
  data_class[GROUP_MATERI_PLASTI_TENSION_DIRECT] = MATERI;
  data_required[GROUP_MATERI_PLASTI_TENSION_DIRECT] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_TENSION_DIRECT_NORMAL],"group_materi_plasti_tension_direct_normal");
  type[GROUP_MATERI_PLASTI_TENSION_DIRECT_NORMAL] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_TENSION_DIRECT_NORMAL] = 3;
  data_class[GROUP_MATERI_PLASTI_TENSION_DIRECT_NORMAL] = MATERI;
  data_required[GROUP_MATERI_PLASTI_TENSION_DIRECT_NORMAL] = GROUP_MATERI_PLASTI_TENSION_DIRECT;
  // ndim-flexible plane normal (manual Professional 6.739)
  fixed_length[GROUP_MATERI_PLASTI_TENSION_DIRECT_NORMAL] = 0;

  strcpy(name[GROUP_MATERI_PLASTI_TENSION_DIRECT_NORMAL_AUTOMATIC],"group_materi_plasti_tension_direct_normal_automatic");
  type[GROUP_MATERI_PLASTI_TENSION_DIRECT_NORMAL_AUTOMATIC] = INTEGER;
  data_length[GROUP_MATERI_PLASTI_TENSION_DIRECT_NORMAL_AUTOMATIC] = 1;
  data_class[GROUP_MATERI_PLASTI_TENSION_DIRECT_NORMAL_AUTOMATIC] = MATERI;
  data_required[GROUP_MATERI_PLASTI_TENSION_DIRECT_NORMAL_AUTOMATIC] = GROUP_MATERI_PLASTI_TENSION_DIRECT;

  strcpy(name[GROUP_MATERI_PLASTI_TENSION_DIRECT_VISCO],"group_materi_plasti_tension_direct_visco");
  type[GROUP_MATERI_PLASTI_TENSION_DIRECT_VISCO] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_TENSION_DIRECT_VISCO] = 1;
  data_class[GROUP_MATERI_PLASTI_TENSION_DIRECT_VISCO] = MATERI;
  data_required[GROUP_MATERI_PLASTI_TENSION_DIRECT_VISCO] = GROUP_MATERI_PLASTI_TENSION_DIRECT;

  strcpy(name[GROUP_MATERI_PLASTI_TENSION_DIRECT_WALL],"group_materi_plasti_tension_direct_wall");
  type[GROUP_MATERI_PLASTI_TENSION_DIRECT_WALL] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_TENSION_DIRECT_WALL] = 1;
  data_class[GROUP_MATERI_PLASTI_TENSION_DIRECT_WALL] = MATERI;
  data_required[GROUP_MATERI_PLASTI_TENSION_DIRECT_WALL] = GROUP_MATERI_PLASTI_TENSION_DIRECT;

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
  // variable length: Professional layout (eta p, manual 6.748) and the
  // legacy GNU layout (eta p f_ref) are both accepted
  fixed_length[GROUP_MATERI_PLASTI_VISCO_POWER] = 0;

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

  strcpy(name[GROUP_TRUSS_EXPANSION],"group_truss_expansion");
  type[GROUP_TRUSS_EXPANSION] = DOUBLE_PRECISION;
  data_length[GROUP_TRUSS_EXPANSION] = 1;
  data_class[GROUP_TRUSS_EXPANSION] = TRUSS;
  data_required[GROUP_TRUSS_EXPANSION] = GROUP_TYPE;

  strcpy(name[GROUP_TRUSS_INITIAL_FORCE],"group_truss_initial_force");
  type[GROUP_TRUSS_INITIAL_FORCE] = DOUBLE_PRECISION;
  data_length[GROUP_TRUSS_INITIAL_FORCE] = 1;
  data_class[GROUP_TRUSS_INITIAL_FORCE] = TRUSS;
  data_required[GROUP_TRUSS_INITIAL_FORCE] = GROUP_TYPE;

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

  strcpy(name[HEX18],"hex18");

  strcpy(name[HEX20],"hex20");

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

  strcpy(name[INCLUDE],"include");
  type[INCLUDE] = INTEGER;
  data_length[INCLUDE] = 1;
  no_index[INCLUDE] = 1;
  external[INCLUDE] = 0;

  strcpy(name[INITIALIZATION_VALUES],"initialization_values");
  type[INITIALIZATION_VALUES] = INTEGER;
  data_length[INITIALIZATION_VALUES] = DATA_ITEM_SIZE;
  no_index[INITIALIZATION_VALUES] = 1;
  external[INITIALIZATION_VALUES] = 0;
  fixed_length[INITIALIZATION_VALUES] = 0;

  strcpy(name[INPUT_ABAQUS],"input_abaqus");
  type[INPUT_ABAQUS] = INTEGER;
  data_length[INPUT_ABAQUS] = 1;
  no_index[INPUT_ABAQUS] = 1;
  external[INPUT_ABAQUS] = 0;

  strcpy(name[INPUT_ABAQUS_CONTINUE],"input_abaqus_continue");
  type[INPUT_ABAQUS_CONTINUE] = INTEGER;
  data_length[INPUT_ABAQUS_CONTINUE] = 1;
  no_index[INPUT_ABAQUS_CONTINUE] = 1;
  external[INPUT_ABAQUS_CONTINUE] = 0;

  strcpy(name[INPUT_ABAQUS_GROUP],"input_abaqus_group");
  type[INPUT_ABAQUS_GROUP] = INTEGER;
  data_length[INPUT_ABAQUS_GROUP] = 1;
  no_index[INPUT_ABAQUS_GROUP] = 1;
  external[INPUT_ABAQUS_GROUP] = 0;

  strcpy(name[INPUT_ABAQUS_MESH],"input_abaqus_mesh");
  type[INPUT_ABAQUS_MESH] = INTEGER;
  data_length[INPUT_ABAQUS_MESH] = 1;
  no_index[INPUT_ABAQUS_MESH] = 1;
  external[INPUT_ABAQUS_MESH] = 0;

  strcpy(name[INPUT_ABAQUS_NAME],"input_abaqus_name");
  type[INPUT_ABAQUS_NAME] = INTEGER;
  data_length[INPUT_ABAQUS_NAME] = DATA_ITEM_SIZE;
  no_index[INPUT_ABAQUS_NAME] = 1;
  external[INPUT_ABAQUS_NAME] = 0;
  fixed_length[INPUT_ABAQUS_NAME] = 0;

  strcpy(name[INPUT_ABAQUS_SET],"input_abaqus_set");
  type[INPUT_ABAQUS_SET] = INTEGER;
  data_length[INPUT_ABAQUS_SET] = DATA_ITEM_SIZE;
  no_index[INPUT_ABAQUS_SET] = 1;
  external[INPUT_ABAQUS_SET] = 0;
  fixed_length[INPUT_ABAQUS_SET] = 0;

  strcpy(name[INPUT_FEFLOW_FEM],"input_feflow_fem");
  type[INPUT_FEFLOW_FEM] = INTEGER;
  data_length[INPUT_FEFLOW_FEM] = 1;
  no_index[INPUT_FEFLOW_FEM] = 1;
  external[INPUT_FEFLOW_FEM] = 0;

  strcpy(name[INPUT_FEFLOW_MESH],"input_feflow_mesh");
  type[INPUT_FEFLOW_MESH] = INTEGER;
  data_length[INPUT_FEFLOW_MESH] = 1;
  no_index[INPUT_FEFLOW_MESH] = 1;
  external[INPUT_FEFLOW_MESH] = 0;

  strcpy(name[INPUT_FEFLOW_MESH_HYDRAULIC_HEAD],"input_feflow_mesh_hydraulic_head");
  type[INPUT_FEFLOW_MESH_HYDRAULIC_HEAD] = INTEGER;
  data_length[INPUT_FEFLOW_MESH_HYDRAULIC_HEAD] = 1;
  no_index[INPUT_FEFLOW_MESH_HYDRAULIC_HEAD] = 1;
  external[INPUT_FEFLOW_MESH_HYDRAULIC_HEAD] = 0;

  strcpy(name[INPUT_GMSH],"input_gmsh");
  type[INPUT_GMSH] = INTEGER;
  data_length[INPUT_GMSH] = 1;
  no_index[INPUT_GMSH] = 1;
  external[INPUT_GMSH] = 0;

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

  strcpy(name[MATERI_DISPLACEMENT_RELATIVE],"materi_displacement_relative");

  strcpy(name[MATERI_DISPLACEMENT_RELATIVE_REF],"materi_displacement_relative_ref");
  type[MATERI_DISPLACEMENT_RELATIVE_REF] = DOUBLE_PRECISION;
  data_length[MATERI_DISPLACEMENT_RELATIVE_REF] = 1;
  no_index[MATERI_DISPLACEMENT_RELATIVE_REF] = 1;
  version_all[MATERI_DISPLACEMENT_RELATIVE_REF] = 1;
  data_class[MATERI_DISPLACEMENT_RELATIVE_REF] = MATERI;

  strcpy(name[MATERI_HISTORY_VARIABLES],"materi_history_variables");

  strcpy(name[MATERI_MAXWELL_STRESS],"materi_maxwell_stress");

  strcpy(name[MATERI_PLASTI_F],"materi_plasti_f");

  strcpy(name[MATERI_PLASTI_F_NONLOCAL],"materi_plasti_f_nonlocal");

  strcpy(name[MATERI_PLASTI_INCREMENTAL_SUBSTEPS],"materi_plasti_incremental_substeps");

  strcpy(name[MATERI_PLASTI_KAPPA],"materi_plasti_kappa");

  strcpy(name[MATERI_PLASTI_KAPPA_SHEAR],"materi_plasti_kappa_shear");

  strcpy(name[MATERI_PLASTI_HYPO_HISTORY],"materi_plasti_hypo_history");

  strcpy(name[MATERI_PLASTI_CAP1_HISTORY],"materi_plasti_cap1_history");

  strcpy(name[MATERI_PLASTI_DIPRISCO_HISTORY],"materi_plasti_diprisco_history");

  strcpy(name[MATERI_PLASTI_HARDSOIL_HISTORY],"materi_plasti_hardsoil_history");

  strcpy(name[MATERI_PLASTI_CAMCLAY_HISTORY],"materi_plasti_camclay_history");

  strcpy(name[MATERI_PLASTI_RHO],"materi_plasti_rho");

  strcpy(name[MATERI_PLASTI_SOFTVAR_LOCAL],"materi_plasti_softvar_local");

  strcpy(name[MATERI_PLASTI_SOFTVAR_NONLOCAL],"materi_plasti_softvar_nonlocal");

  strcpy(name[MATERI_ROTATION],"materi_rotation");

  strcpy(name[MATERI_STRAINENERGY],"materi_strainenergy");

  strcpy(name[MATERI_STRAIN_ELASTI],"materi_strain_elasti");

  strcpy(name[MATERI_STRAIN_INTERGRANULAR],"materi_strain_intergranular");

  strcpy(name[MATERI_STRAIN_PLASTI],"materi_strain_plasti");

  strcpy(name[MATERI_STRAIN_PLASTI_CAP],"materi_strain_plasti_cap");

  strcpy(name[MATERI_STRAIN_PLASTI_COMPRESSION],"materi_strain_plasti_compression");

  strcpy(name[MATERI_STRAIN_PLASTI_DIPRISCO],"materi_strain_plasti_diprisco");

  strcpy(name[MATERI_STRAIN_PLASTI_DRUCKPRAG],"materi_strain_plasti_druckprag");

  strcpy(name[MATERI_STRAIN_PLASTI_HARDSOIL],"materi_strain_plasti_hardsoil");

  strcpy(name[MATERI_STRAIN_TOTAL],"materi_strain_total" );

  strcpy(name[MATERI_STRESS],"materi_stress");

  strcpy(name[MATERI_STRESS_PRESSURE_HISTORY],"materi_stress_pressure_history");

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
  strcpy(name[MATRIX_PARDISO],"matrix_pardiso");

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

  strcpy(name[MESH_ACTIVATE_GRAVITY_ELEMENT],"mesh_activate_gravity_element");
  type[MESH_ACTIVATE_GRAVITY_ELEMENT] = INTEGER;
  data_length[MESH_ACTIVATE_GRAVITY_ELEMENT] = DATA_ITEM_SIZE;
  fixed_length[MESH_ACTIVATE_GRAVITY_ELEMENT] = 0;
  data_class[MESH_ACTIVATE_GRAVITY_ELEMENT] = CONTROL;

  strcpy(name[MESH_ACTIVATE_GRAVITY_ELEMENT_GROUP],"mesh_activate_gravity_element_group");
  type[MESH_ACTIVATE_GRAVITY_ELEMENT_GROUP] = INTEGER;
  data_length[MESH_ACTIVATE_GRAVITY_ELEMENT_GROUP] = DATA_ITEM_SIZE;
  fixed_length[MESH_ACTIVATE_GRAVITY_ELEMENT_GROUP] = 0;
  data_class[MESH_ACTIVATE_GRAVITY_ELEMENT_GROUP] = CONTROL;

  strcpy(name[MESH_ACTIVATE_GRAVITY_GEOMETRY],"mesh_activate_gravity_geometry");
  type[MESH_ACTIVATE_GRAVITY_GEOMETRY] = INTEGER;
  data_length[MESH_ACTIVATE_GRAVITY_GEOMETRY] = 2;
  data_class[MESH_ACTIVATE_GRAVITY_GEOMETRY] = CONTROL;

  strcpy(name[MESH_ACTIVATE_GRAVITY_METHOD],"mesh_activate_gravity_method");
  type[MESH_ACTIVATE_GRAVITY_METHOD] = INTEGER;
  data_length[MESH_ACTIVATE_GRAVITY_METHOD] = 1;
  data_class[MESH_ACTIVATE_GRAVITY_METHOD] = CONTROL;

  strcpy(name[MESH_ACTIVATE_GRAVITY_STIFFNESS_FACTOR],"mesh_activate_gravity_stiffness_factor");
  type[MESH_ACTIVATE_GRAVITY_STIFFNESS_FACTOR] = DOUBLE_PRECISION;
  data_length[MESH_ACTIVATE_GRAVITY_STIFFNESS_FACTOR] = 1;
  data_class[MESH_ACTIVATE_GRAVITY_STIFFNESS_FACTOR] = CONTROL;

  strcpy(name[MESH_ACTIVATE_GRAVITY_TIME],"mesh_activate_gravity_time");
  type[MESH_ACTIVATE_GRAVITY_TIME] = DOUBLE_PRECISION;
  data_length[MESH_ACTIVATE_GRAVITY_TIME] = 2;
  data_class[MESH_ACTIVATE_GRAVITY_TIME] = CONTROL;

  strcpy(name[MESH_ACTIVATE_GRAVITY_TIME_INITIAL],"mesh_activate_gravity_time_initial");
  type[MESH_ACTIVATE_GRAVITY_TIME_INITIAL] = DOUBLE_PRECISION;
  data_length[MESH_ACTIVATE_GRAVITY_TIME_INITIAL] = 1;
  data_class[MESH_ACTIVATE_GRAVITY_TIME_INITIAL] = CONTROL;

  strcpy(name[MESH_ACTIVATE_GRAVITY_TIME_STRAIN_SETTLEMENT],"mesh_activate_gravity_time_strain_settlement");
  type[MESH_ACTIVATE_GRAVITY_TIME_STRAIN_SETTLEMENT] = INTEGER;
  data_length[MESH_ACTIVATE_GRAVITY_TIME_STRAIN_SETTLEMENT] = 1;
  data_class[MESH_ACTIVATE_GRAVITY_TIME_STRAIN_SETTLEMENT] = CONTROL;

  strcpy(name[METHOD1],"method1");

  strcpy(name[METHOD2],"method2");

  strcpy(name[MIDDLE],"middle");

  strcpy(name[MINIMAL],"minimal");

  strcpy(name[MINUS_ONE],"minus_one" );

  strcpy(name[MOMENT],"moment");

  strcpy(name[NEGATIVE],"negative");

  strcpy(name[NO],"no");

  strcpy(name[OPENED],"opened");
  strcpy(name[CLOSED],"closed");

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

  strcpy(name[NODE_CONVECTION_APPLY],"node_convection_apply");
  type[NODE_CONVECTION_APPLY] = INTEGER;
  data_length[NODE_CONVECTION_APPLY] = 1;
  version_all[NODE_CONVECTION_APPLY] = 1;
  data_class[NODE_CONVECTION_APPLY] = NODE;
  data_required[NODE_CONVECTION_APPLY] = NODE;

  strcpy(name[NODE_DYNAMIC_PRESSURE],"node_dynamic_pressure");
  type[NODE_DYNAMIC_PRESSURE] = DOUBLE_PRECISION;
  data_length[NODE_DYNAMIC_PRESSURE] = 1;
  version_all[NODE_DYNAMIC_PRESSURE] = 1;
  data_class[NODE_DYNAMIC_PRESSURE] = NODE;
  data_required[NODE_DYNAMIC_PRESSURE] = NODE;

  strcpy(name[NODE_FORCE],"node_force");
  type[NODE_FORCE] = DOUBLE_PRECISION;
  data_length[NODE_FORCE] = ndim;
  version_all[NODE_FORCE] = 1;
  data_class[NODE_FORCE] = NODE;
  data_required[NODE_FORCE] = NODE;

  strcpy(name[NODE_INERTIA],"node_inertia");
  type[NODE_INERTIA] = DOUBLE_PRECISION;
  data_length[NODE_INERTIA] = ndim;
  version_all[NODE_INERTIA] = 1;
  data_class[NODE_INERTIA] = NODE;
  data_required[NODE_INERTIA] = NODE;

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

  // node_dof_previous_step: INTERNAL snapshot of the node dofs at the
  // START of the current time step (captured by top() before the
  // equilibrium iterations, one copy per step). Consumed by the
  // post_calcul operators -materi_stress -young_apparent and
  // -poisson_apparent (manual Professional 6.903): the apparent E and
  // Poisson ratio are determined from the INCREMENTAL strains and
  // stresses of the last time step, i.e. the difference between the
  // converged node dofs (VERSION_NORMAL at step_close) and this
  // snapshot. Not printed to the .dbs (external 0) and not versioned.
  strcpy(name[NODE_DOF_PREVIOUS_STEP],"node_dof_previous_step");
  type[NODE_DOF_PREVIOUS_STEP] = DOUBLE_PRECISION;
  data_length[NODE_DOF_PREVIOUS_STEP] = nuknwn;
  external[NODE_DOF_PREVIOUS_STEP] = 0;
  data_class[NODE_DOF_PREVIOUS_STEP] = NODE;
  data_required[NODE_DOF_PREVIOUS_STEP] = NODE;

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
  // FIX (2026-08-24): was 1, so node_mass in 2D/3D consumed only the
  // first component and the parser choked on the rest (node_mass was
  // effectively unusable since the GNU origins; node_damping and
  // node_stiffness correctly use ndim)
  data_length[NODE_MASS] = ndim;
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

  strcpy(name[NODE_RHSIDE_PREVIOUS],"node_rhside_previous");
  type[NODE_RHSIDE_PREVIOUS] = DOUBLE_PRECISION;
  data_length[NODE_RHSIDE_PREVIOUS] = npuknwn;
  data_class[NODE_RHSIDE_PREVIOUS] = NODE;

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

  strcpy(name[NODE_DEFORMED_MESH],"node_deformed_mesh");

  strcpy(name[NODE_MESH],"node_mesh");
  type[NODE_MESH] = INTEGER;
  data_length[NODE_MESH] = 1;
  version_all[NODE_MESH] = 1;
  data_class[NODE_MESH] = NODE;
  data_required[NODE_MESH] = NODE;

  strcpy(name[NODE_SLIDE],"node_slide");
  type[NODE_SLIDE] = INTEGER;
  data_length[NODE_SLIDE] = 1;
  version_all[NODE_SLIDE] = 1;
  data_class[NODE_SLIDE] = NODE;
  data_required[NODE_SLIDE] = NODE;

  // Output records of the elastic-plastic slide law (slide.cc, manual
  // Professional 6.1042-6.1047): node_slide_direction = the local slide
  // frame (3 values normal + 3 values tangential), node_slide_f = the
  // plastic yield function value and node_slide_force = the slide forces
  // (first value = force on the slide geometry along the normal, second
  // = force on the slide geometry along the tangential direction; scaled
  // by 2*pi*r in axisymmetric problems). Computed per slide node at every
  // equilibrium iteration; the target_item tests (slide2) read them at
  // the end of the calculation.
  strcpy(name[NODE_SLIDE_DIRECTION],"node_slide_direction");
  type[NODE_SLIDE_DIRECTION] = DOUBLE_PRECISION;
  data_length[NODE_SLIDE_DIRECTION] = 6;
  version_all[NODE_SLIDE_DIRECTION] = 1;
  data_class[NODE_SLIDE_DIRECTION] = NODE;
  data_required[NODE_SLIDE_DIRECTION] = NODE;

  strcpy(name[NODE_SLIDE_F],"node_slide_f");
  type[NODE_SLIDE_F] = DOUBLE_PRECISION;
  data_length[NODE_SLIDE_F] = 1;
  version_all[NODE_SLIDE_F] = 1;
  data_class[NODE_SLIDE_F] = NODE;
  data_required[NODE_SLIDE_F] = NODE;

  strcpy(name[NODE_SLIDE_FORCE],"node_slide_force");
  type[NODE_SLIDE_FORCE] = DOUBLE_PRECISION;
  data_length[NODE_SLIDE_FORCE] = ndim;
  version_all[NODE_SLIDE_FORCE] = 1;
  data_class[NODE_SLIDE_FORCE] = NODE;
  data_required[NODE_SLIDE_FORCE] = NODE;

  strcpy(name[NODE_STATIC_PRESSURE],"node_static_pressure");
  type[NODE_STATIC_PRESSURE] = DOUBLE_PRECISION;
  data_length[NODE_STATIC_PRESSURE] = 1;
  version_all[NODE_STATIC_PRESSURE] = 1;
  data_class[NODE_STATIC_PRESSURE] = NODE;
  data_required[NODE_STATIC_PRESSURE] = NODE;

  strcpy(name[NODE_TOTAL_PRESSURE],"node_total_pressure");
  type[NODE_TOTAL_PRESSURE] = DOUBLE_PRECISION;
  data_length[NODE_TOTAL_PRESSURE] = 1;
  version_all[NODE_TOTAL_PRESSURE] = 1;
  data_class[NODE_TOTAL_PRESSURE] = NODE;
  data_required[NODE_TOTAL_PRESSURE] = NODE;
  type[NODE_DEFORMED_MESH] = DOUBLE_PRECISION;
  data_length[NODE_DEFORMED_MESH] = 1;
  version_all[NODE_DEFORMED_MESH] = 1;
  data_class[NODE_DEFORMED_MESH] = NODE;
  data_required[NODE_DEFORMED_MESH] = NODE;

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

  strcpy(name[LOGNORMAL],"lognormal");

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

  // nonlocal_name (manual Professional 6.898): the name of the
  // plasticity model treated nonlocal. The GNU applies the nonlocal
  // yield rule contribution to every model (plasti.cc), so the record
  // is accepted but has no per-model gate (partial).
  strcpy(name[NONLOCAL_NAME],"nonlocal_name");
  type[NONLOCAL_NAME] = INTEGER;
  data_length[NONLOCAL_NAME] = 1;
  no_index[NONLOCAL_NAME] = 1;

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

  strcpy(name[POST_CALCUL_MATERI_STRESS_FORCE_AVERAGE],"post_calcul_materi_stress_force_average");
  type[POST_CALCUL_MATERI_STRESS_FORCE_AVERAGE] = INTEGER;
  data_length[POST_CALCUL_MATERI_STRESS_FORCE_AVERAGE] = 1;
  no_index[POST_CALCUL_MATERI_STRESS_FORCE_AVERAGE] = 1;
  data_class[POST_CALCUL_MATERI_STRESS_FORCE_AVERAGE] = POST;
  data_required[POST_CALCUL_MATERI_STRESS_FORCE_AVERAGE] = POST_CALCUL;

  strcpy(name[POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_EXCLUDE],"post_calcul_materi_stress_force_direction_exclude");
  type[POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_EXCLUDE] = DOUBLE_PRECISION;
  data_length[POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_EXCLUDE] = DATA_ITEM_SIZE;
  fixed_length[POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_EXCLUDE] = 0;
  no_index[POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_EXCLUDE] = 1;
  data_class[POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_EXCLUDE] = POST;
  data_required[POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_EXCLUDE] = POST_CALCUL;

  strcpy(name[POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_EXCLUDE_EPSILON],"post_calcul_materi_stress_force_direction_exclude_epsilon");
  type[POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_EXCLUDE_EPSILON] = DOUBLE_PRECISION;
  data_length[POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_EXCLUDE_EPSILON] = 1;
  no_index[POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_EXCLUDE_EPSILON] = 1;
  data_class[POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_EXCLUDE_EPSILON] = POST;
  data_required[POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_EXCLUDE_EPSILON] = POST_CALCUL;

  strcpy(name[POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_INCLUDE],"post_calcul_materi_stress_force_direction_include");
  type[POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_INCLUDE] = DOUBLE_PRECISION;
  data_length[POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_INCLUDE] = DATA_ITEM_SIZE;
  fixed_length[POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_INCLUDE] = 0;
  no_index[POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_INCLUDE] = 1;
  data_class[POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_INCLUDE] = POST;
  data_required[POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_INCLUDE] = POST_CALCUL;

  strcpy(name[POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_INCLUDE_EPSILON],"post_calcul_materi_stress_force_direction_include_epsilon");
  type[POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_INCLUDE_EPSILON] = DOUBLE_PRECISION;
  data_length[POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_INCLUDE_EPSILON] = 1;
  no_index[POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_INCLUDE_EPSILON] = 1;
  data_class[POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_INCLUDE_EPSILON] = POST;
  data_required[POST_CALCUL_MATERI_STRESS_FORCE_DIRECTION_INCLUDE_EPSILON] = POST_CALCUL;

  strcpy(name[POST_CALCUL_MATERI_STRESS_FORCE_ELEMENT_GROUP],"post_calcul_materi_stress_force_element_group");
  type[POST_CALCUL_MATERI_STRESS_FORCE_ELEMENT_GROUP] = INTEGER;
  data_length[POST_CALCUL_MATERI_STRESS_FORCE_ELEMENT_GROUP] = DATA_ITEM_SIZE;
  fixed_length[POST_CALCUL_MATERI_STRESS_FORCE_ELEMENT_GROUP] = 0;
  no_index[POST_CALCUL_MATERI_STRESS_FORCE_ELEMENT_GROUP] = 1;
  data_class[POST_CALCUL_MATERI_STRESS_FORCE_ELEMENT_GROUP] = POST;
  data_required[POST_CALCUL_MATERI_STRESS_FORCE_ELEMENT_GROUP] = POST_CALCUL;

  strcpy(name[POST_CALCUL_MATERI_STRESS_FORCE_OUTER],"post_calcul_materi_stress_force_outer");
  type[POST_CALCUL_MATERI_STRESS_FORCE_OUTER] = INTEGER;
  data_length[POST_CALCUL_MATERI_STRESS_FORCE_OUTER] = 1;
  no_index[POST_CALCUL_MATERI_STRESS_FORCE_OUTER] = 1;
  data_class[POST_CALCUL_MATERI_STRESS_FORCE_OUTER] = POST;
  data_required[POST_CALCUL_MATERI_STRESS_FORCE_OUTER] = POST_CALCUL;

  strcpy(name[POST_CALCUL_MATERI_STRESS_FORCE_PLOT_SWITCH],"post_calcul_materi_stress_force_plot_switch");
  type[POST_CALCUL_MATERI_STRESS_FORCE_PLOT_SWITCH] = INTEGER;
  data_length[POST_CALCUL_MATERI_STRESS_FORCE_PLOT_SWITCH] = DATA_ITEM_SIZE;
  fixed_length[POST_CALCUL_MATERI_STRESS_FORCE_PLOT_SWITCH] = 0;
  no_index[POST_CALCUL_MATERI_STRESS_FORCE_PLOT_SWITCH] = 1;
  data_class[POST_CALCUL_MATERI_STRESS_FORCE_PLOT_SWITCH] = POST;
  data_required[POST_CALCUL_MATERI_STRESS_FORCE_PLOT_SWITCH] = POST_CALCUL;

  strcpy(name[POST_CALCUL_MATERI_STRESS_FORCE_REFERENCE_POINT],"post_calcul_materi_stress_force_reference_point");
  type[POST_CALCUL_MATERI_STRESS_FORCE_REFERENCE_POINT] = DOUBLE_PRECISION;
  data_length[POST_CALCUL_MATERI_STRESS_FORCE_REFERENCE_POINT] = DATA_ITEM_SIZE;
  fixed_length[POST_CALCUL_MATERI_STRESS_FORCE_REFERENCE_POINT] = 0;
  no_index[POST_CALCUL_MATERI_STRESS_FORCE_REFERENCE_POINT] = 1;
  data_class[POST_CALCUL_MATERI_STRESS_FORCE_REFERENCE_POINT] = POST;
  data_required[POST_CALCUL_MATERI_STRESS_FORCE_REFERENCE_POINT] = POST_CALCUL;

  strcpy(name[POST_CALCUL_MATERI_STRESS_FORCE_THICKNESS_SWITCH],"post_calcul_materi_stress_force_thickness_switch");
  type[POST_CALCUL_MATERI_STRESS_FORCE_THICKNESS_SWITCH] = INTEGER;
  data_length[POST_CALCUL_MATERI_STRESS_FORCE_THICKNESS_SWITCH] = DATA_ITEM_SIZE;
  fixed_length[POST_CALCUL_MATERI_STRESS_FORCE_THICKNESS_SWITCH] = 0;
  no_index[POST_CALCUL_MATERI_STRESS_FORCE_THICKNESS_SWITCH] = 1;
  data_class[POST_CALCUL_MATERI_STRESS_FORCE_THICKNESS_SWITCH] = POST;
  data_required[POST_CALCUL_MATERI_STRESS_FORCE_THICKNESS_SWITCH] = POST_CALCUL;

  // per-node flag record of the -force family: -YES when the node
  // received the AVERAGED results of the quad9 middle plane
  // (post_calcul_materi_stress_force_average -yes), consumed by
  // control_print_materi_stress_force -primary through the
  // msf_node_is_averaged() hook (print_materi_stress_force.cc).
  // NODE class + version_all=1 so that db_version_copy and the
  // renumbering of the print (VERSION_PRINT) carry it like
  // NODE_DOF_CALCUL.
  strcpy(name[POST_CALCUL_MATERI_STRESS_FORCE_AVERAGED_NODE],"post_calcul_materi_stress_force_averaged_node");
  type[POST_CALCUL_MATERI_STRESS_FORCE_AVERAGED_NODE] = INTEGER;
  data_length[POST_CALCUL_MATERI_STRESS_FORCE_AVERAGED_NODE] = 1;
  version_all[POST_CALCUL_MATERI_STRESS_FORCE_AVERAGED_NODE] = 1;
  data_class[POST_CALCUL_MATERI_STRESS_FORCE_AVERAGED_NODE] = NODE;
  data_required[POST_CALCUL_MATERI_STRESS_FORCE_AVERAGED_NODE] = NODE;

  // item names of post_calcul -materi_stress -force (manual
  // Professional 6.913): -norx_sig, -nory_sig, -nors_sig, ... are the
  // VALUE names of target_item records against the flat NODE_DOF_CALCUL
  // layout (target_item N -node_dof_calcul <node> -norx_sig). These are
  // PURE NAME entries (no record: type/class stay default): exit_tn
  // (miscel.cc) maps them to the slot of the item via post_calcul_names.
  strcpy(name[NORX_SIG],"norx_sig");
  strcpy(name[NORY_SIG],"nory_sig");
  strcpy(name[NORZ_SIG],"norz_sig");
  strcpy(name[NORS_SIG],"nors_sig");
  strcpy(name[SHEX_SIG],"shex_sig");
  strcpy(name[SHEY_SIG],"shey_sig");
  strcpy(name[SHEZ_SIG],"shez_sig");
  strcpy(name[SHES_SIG],"shes_sig");
  strcpy(name[MOMX_SIG],"momx_sig");
  strcpy(name[MOMY_SIG],"momy_sig");
  strcpy(name[MOMS_SIG],"moms_sig");
  strcpy(name[MOM1X_SIG],"mom1x_sig");
  strcpy(name[MOM1Y_SIG],"mom1y_sig");
  strcpy(name[MOM1Z_SIG],"mom1z_sig");
  strcpy(name[MOM1S_SIG],"mom1s_sig");
  strcpy(name[MOM2X_SIG],"mom2x_sig");
  strcpy(name[MOM2Y_SIG],"mom2y_sig");
  strcpy(name[MOM2Z_SIG],"mom2z_sig");
  strcpy(name[MOM2S_SIG],"mom2s_sig");

  // item names of the groundflow pressure split (manual Professional
  // 6.913 area): post_calcul -groundflow_pressure -total_pressure/
  // -static_pressure/-dynamic_pressure generate the item names -to_pres,
  // -st_pres, -dy_pres (target_item N -post_point_dof_calcul <post> -to_pres).
  // Pure name entries, same resolution as the *_sig family.
  strcpy(name[TO_PRES],"to_pres");
  strcpy(name[ST_PRES],"st_pres");
  strcpy(name[DY_PRES],"dy_pres");

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

  strcpy(name[POST_FORCE_EDGE_SUMMED],"post_force_edge_summed");
  type[POST_FORCE_EDGE_SUMMED] = DOUBLE_PRECISION;
  data_length[POST_FORCE_EDGE_SUMMED] = ndim;
  data_class[POST_FORCE_EDGE_SUMMED] = POST;
  no_index[POST_FORCE_EDGE_SUMMED] = 1;

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

  strcpy(name[POST_POINT_MOVE],"post_point_move");
  type[POST_POINT_MOVE] = INTEGER;
  data_length[POST_POINT_MOVE] = 1;
  no_index[POST_POINT_MOVE] = 1;
  data_class[POST_POINT_MOVE] = POST;

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

  strcpy(name[PRIMARY],"primary");

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

  // print_debug -yes/-no (Professional input of slide3): switches the
  // routine-level debug prints (set_swit) on/off. Parsed for
  // compatibility; the GNU's debug output is driven by print_where.
  strcpy(name[PRINT_DEBUG],"print_debug");
  type[PRINT_DEBUG] = INTEGER;
  data_length[PRINT_DEBUG] = 1;
  data_class[PRINT_DEBUG] = PRINT;
  no_index[PRINT_DEBUG] = 1;

  strcpy(name[PRINT_FAILURE],"print_failure");
  type[PRINT_FAILURE] = INTEGER;
  data_length[PRINT_FAILURE] = 1;
  data_class[PRINT_FAILURE] = PRINT;
  no_index[PRINT_FAILURE] = 1;

  strcpy(name[PRINT_MESH_DOF],"print_mesh_dof");
  type[PRINT_MESH_DOF] = INTEGER;
  data_length[PRINT_MESH_DOF] = DATA_ITEM_SIZE;
  fixed_length[PRINT_MESH_DOF] = 0;
  no_index[PRINT_MESH_DOF] = 1;
  data_class[PRINT_MESH_DOF] = PRINT;

  strcpy(name[PRINT_MESH_DOF_GEOMETRY],"print_mesh_dof_geometry");
  type[PRINT_MESH_DOF_GEOMETRY] = INTEGER;
  data_length[PRINT_MESH_DOF_GEOMETRY] = 2;
  no_index[PRINT_MESH_DOF_GEOMETRY] = 1;
  data_class[PRINT_MESH_DOF_GEOMETRY] = PRINT;

  strcpy(name[PRINT_MESH_DOF_VALUES],"print_mesh_dof_values");
  type[PRINT_MESH_DOF_VALUES] = DOUBLE_PRECISION;
  data_length[PRINT_MESH_DOF_VALUES] = DATA_ITEM_SIZE;
  fixed_length[PRINT_MESH_DOF_VALUES] = 0;
  no_index[PRINT_MESH_DOF_VALUES] = 1;
  data_class[PRINT_MESH_DOF_VALUES] = PRINT;

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
  strcpy(name[QUAD6],"quad6");
  strcpy(name[QUAD8],"quad8");

  strcpy(name[QUAD9],"quad9");

  strcpy(name[QUAD16],"quad16");

  strcpy(name[RA],"ra");

  strcpy(name[RECTANGLE],"rectangle" );

  strcpy(name[REPEAT_CALCULATE_RESULT],"repeat_calculate_result");
  type[REPEAT_CALCULATE_RESULT] = DOUBLE_PRECISION;
  data_length[REPEAT_CALCULATE_RESULT] = 2;
  data_class[REPEAT_CALCULATE_RESULT] = CALCUL;

  strcpy(name[REPEAT_SAVE_RESULT],"repeat_save_result");
  type[REPEAT_SAVE_RESULT] = DOUBLE_PRECISION;
  data_length[REPEAT_SAVE_RESULT] = MCALCUL;
  fixed_length[REPEAT_SAVE_RESULT] = 0;
  data_class[REPEAT_SAVE_RESULT] = CALCUL;

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

  strcpy(name[SIZE_TOT_LARGE],"size_tot_large");

  strcpy(name[SINUS],"sinus");

  strcpy(name[SLES],"sles");

  strcpy(name[SLIDE_AXISYMMETRIC],"slide_axisymmetric");
  type[SLIDE_AXISYMMETRIC] = INTEGER;
  data_length[SLIDE_AXISYMMETRIC] = 1;
  data_class[SLIDE_AXISYMMETRIC] = SLIDE;

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

  // slide_plasti_friction (manual Professional 6.1042): phi c — friction
  // angle in RADIANS + cohesion of the slide (slide_geometry) law.
  // Maximum slide friction force = c + Fn*tan(phi).
  strcpy(name[SLIDE_PLASTI_FRICTION],"slide_plasti_friction");
  type[SLIDE_PLASTI_FRICTION] = DOUBLE_PRECISION;
  data_length[SLIDE_PLASTI_FRICTION] = 2;
  data_class[SLIDE_PLASTI_FRICTION] = SLIDE;
  data_required[SLIDE_PLASTI_FRICTION] = SLIDE_GEOMETRY;

  // slide_plasti_tension (manual Professional 6.1043): maximum tensile
  // (pull-off) force of the slide connection.
  strcpy(name[SLIDE_PLASTI_TENSION],"slide_plasti_tension");
  type[SLIDE_PLASTI_TENSION] = DOUBLE_PRECISION;
  data_length[SLIDE_PLASTI_TENSION] = 1;
  data_class[SLIDE_PLASTI_TENSION] = SLIDE;
  data_required[SLIDE_PLASTI_TENSION] = SLIDE_GEOMETRY;

  // slide_stiffness (manual Professional 6.1046): stiffness_n
  // stiffness_t of the elastic slide connection (per node).
  strcpy(name[SLIDE_STIFFNESS],"slide_stiffness");
  type[SLIDE_STIFFNESS] = DOUBLE_PRECISION;
  data_length[SLIDE_STIFFNESS] = 2;
  data_class[SLIDE_STIFFNESS] = SLIDE;
  data_required[SLIDE_STIFFNESS] = SLIDE_GEOMETRY;

  // slide_plasti_residual_stiffness (manual Professional 6.1047 area):
  // residual stiffness fraction after plastification; the Professional
  // writes default 1e-2 1e-2 into its .dbs when the record is absent.
  strcpy(name[SLIDE_PLASTI_RESIDUAL_STIFFNESS],"slide_plasti_residual_stiffness");
  type[SLIDE_PLASTI_RESIDUAL_STIFFNESS] = DOUBLE_PRECISION;
  data_length[SLIDE_PLASTI_RESIDUAL_STIFFNESS] = 2;
  data_class[SLIDE_PLASTI_RESIDUAL_STIFFNESS] = SLIDE;
  data_required[SLIDE_PLASTI_RESIDUAL_STIFFNESS] = SLIDE_GEOMETRY;

  // control_slide_plasti_apply (manual Professional 6.371): -no turns
  // the slide plasti law off (keeps the elastic slide).
  strcpy(name[CONTROL_SLIDE_PLASTI_APPLY],"control_slide_plasti_apply");
  type[CONTROL_SLIDE_PLASTI_APPLY] = INTEGER;
  data_length[CONTROL_SLIDE_PLASTI_APPLY] = 1;
  data_class[CONTROL_SLIDE_PLASTI_APPLY] = CONTROL;

  // control_slide_stiffness_apply (manual Professional 6.372): -no turns
  // the elastic slide stiffness off.
  strcpy(name[CONTROL_SLIDE_STIFFNESS_APPLY],"control_slide_stiffness_apply");
  type[CONTROL_SLIDE_STIFFNESS_APPLY] = INTEGER;
  data_length[CONTROL_SLIDE_STIFFNESS_APPLY] = 1;
  data_class[CONTROL_SLIDE_STIFFNESS_APPLY] = CONTROL;

  strcpy(name[SOR],"sor");

  strcpy(name[SPHERE],"sphere");

  strcpy(name[SPRING],"spring");

  strcpy(name[SPRING1],"spring1");

  strcpy(name[SPRING2],"spring2");

  strcpy(name[STATIC],"static" );

  strcpy(name[STEP],"step");

  strcpy(name[STRAIN_SETTLEMENT_DIAGRAM],"strain_settlement_diagram");
  type[STRAIN_SETTLEMENT_DIAGRAM] = DOUBLE_PRECISION;
  data_length[STRAIN_SETTLEMENT_DIAGRAM] = DATA_ITEM_SIZE;
  fixed_length[STRAIN_SETTLEMENT_DIAGRAM] = 0;
  data_class[STRAIN_SETTLEMENT_DIAGRAM] = CONTROL;
  data_required[STRAIN_SETTLEMENT_DIAGRAM] = STRAIN_SETTLEMENT_PARAMETERS;

  strcpy(name[STRAIN_SETTLEMENT_DIAGRAM_DOF],"strain_settlement_diagram_dof");
  type[STRAIN_SETTLEMENT_DIAGRAM_DOF] = INTEGER;
  data_length[STRAIN_SETTLEMENT_DIAGRAM_DOF] = 1;
  data_class[STRAIN_SETTLEMENT_DIAGRAM_DOF] = CONTROL;
  data_required[STRAIN_SETTLEMENT_DIAGRAM_DOF] = STRAIN_SETTLEMENT_DIAGRAM;

  strcpy(name[STRAIN_SETTLEMENT_DIAGRAM_NUMBER],"strain_settlement_diagram_number");
  type[STRAIN_SETTLEMENT_DIAGRAM_NUMBER] = INTEGER;
  data_length[STRAIN_SETTLEMENT_DIAGRAM_NUMBER] = 1;
  data_class[STRAIN_SETTLEMENT_DIAGRAM_NUMBER] = CONTROL;
  data_required[STRAIN_SETTLEMENT_DIAGRAM_NUMBER] = STRAIN_SETTLEMENT_DIAGRAM;

  strcpy(name[STRAIN_SETTLEMENT_ELEMENT_GROUP],"strain_settlement_element_group");
  type[STRAIN_SETTLEMENT_ELEMENT_GROUP] = INTEGER;
  data_length[STRAIN_SETTLEMENT_ELEMENT_GROUP] = DATA_ITEM_SIZE;
  fixed_length[STRAIN_SETTLEMENT_ELEMENT_GROUP] = 0;
  data_class[STRAIN_SETTLEMENT_ELEMENT_GROUP] = CONTROL;
  data_required[STRAIN_SETTLEMENT_ELEMENT_GROUP] = STRAIN_SETTLEMENT_PARAMETERS;

  strcpy(name[STRAIN_SETTLEMENT_PARAMETERS],"strain_settlement_parameters");
  type[STRAIN_SETTLEMENT_PARAMETERS] = DOUBLE_PRECISION;
  data_length[STRAIN_SETTLEMENT_PARAMETERS] = 6;
  data_class[STRAIN_SETTLEMENT_PARAMETERS] = CONTROL;

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

  strcpy(name[TANGENT],"tangent");

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

  strcpy(name[UPDATED_AREA],"updated_area");

  strcpy(name[UPDATED_LINEAR],"updated_linear");

  strcpy(name[UPDATED_WITHOUT_ROTATION],"updated_without_rotation");

  strcpy(name[USE],"use" );

  strcpy(name[USER],"user" );

  strcpy(name[VALUE],"value");

  strcpy(name[VECTOR],"vector");


  // Sprint 13 lote 1: support_edge_normal - the distributed Winkler
  // support of an edge (manual Professional 6.1067-6.1078). Master
  // record: stiffness_normal + stiffness_tangential per unit length
  // (2D) / area (3D). Companions: the side selectors (geometry /
  // element_side / element_node / node) and the element_group
  // restriction, same-index. Output: node_support_edge_normal_force
  // (the consistent nodal support forces, filled by area()).
  strcpy(name[SUPPORT_EDGE_NORMAL],"support_edge_normal");
  type[SUPPORT_EDGE_NORMAL] = DOUBLE_PRECISION;
  data_length[SUPPORT_EDGE_NORMAL] = 2;
  data_class[SUPPORT_EDGE_NORMAL] = FORCE;
  // the _GEOMETRY record is required by the area() machinery (it can
  // carry a geometry entity OR a node list); the other companions
  // (element_side / element_node / node / element_group) are
  // additional restrictions on top of it, manual Professional 6.1067
  // "Also the record support_edge_normal_geometry should be specified"
  data_required[SUPPORT_EDGE_NORMAL] = SUPPORT_EDGE_NORMAL_GEOMETRY;

  strcpy(name[SUPPORT_EDGE_NORMAL_GEOMETRY],"support_edge_normal_geometry");
  type[SUPPORT_EDGE_NORMAL_GEOMETRY] = INTEGER;
  data_length[SUPPORT_EDGE_NORMAL_GEOMETRY] = DATA_ITEM_SIZE;
  fixed_length[SUPPORT_EDGE_NORMAL_GEOMETRY] = 0;
  data_class[SUPPORT_EDGE_NORMAL_GEOMETRY] = FORCE;
  data_required[SUPPORT_EDGE_NORMAL_GEOMETRY] = SUPPORT_EDGE_NORMAL;

  strcpy(name[SUPPORT_EDGE_NORMAL_ELEMENT],"support_edge_normal_element");
  type[SUPPORT_EDGE_NORMAL_ELEMENT] = INTEGER;
  data_length[SUPPORT_EDGE_NORMAL_ELEMENT] = DATA_ITEM_SIZE;
  fixed_length[SUPPORT_EDGE_NORMAL_ELEMENT] = 0;
  data_class[SUPPORT_EDGE_NORMAL_ELEMENT] = FORCE;
  data_required[SUPPORT_EDGE_NORMAL_ELEMENT] = SUPPORT_EDGE_NORMAL;

  strcpy(name[SUPPORT_EDGE_NORMAL_ELEMENT_GROUP],"support_edge_normal_element_group");
  type[SUPPORT_EDGE_NORMAL_ELEMENT_GROUP] = INTEGER;
  data_length[SUPPORT_EDGE_NORMAL_ELEMENT_GROUP] = DATA_ITEM_SIZE;
  fixed_length[SUPPORT_EDGE_NORMAL_ELEMENT_GROUP] = 0;
  data_class[SUPPORT_EDGE_NORMAL_ELEMENT_GROUP] = FORCE;
  data_required[SUPPORT_EDGE_NORMAL_ELEMENT_GROUP] = SUPPORT_EDGE_NORMAL;

  strcpy(name[SUPPORT_EDGE_NORMAL_ELEMENT_SIDE],"support_edge_normal_element_side");
  type[SUPPORT_EDGE_NORMAL_ELEMENT_SIDE] = INTEGER;
  data_length[SUPPORT_EDGE_NORMAL_ELEMENT_SIDE] = DATA_ITEM_SIZE;
  fixed_length[SUPPORT_EDGE_NORMAL_ELEMENT_SIDE] = 0;
  data_class[SUPPORT_EDGE_NORMAL_ELEMENT_SIDE] = FORCE;
  data_required[SUPPORT_EDGE_NORMAL_ELEMENT_SIDE] = SUPPORT_EDGE_NORMAL;

  strcpy(name[SUPPORT_EDGE_NORMAL_ELEMENT_NODE],"support_edge_normal_element_node");
  type[SUPPORT_EDGE_NORMAL_ELEMENT_NODE] = INTEGER;
  data_length[SUPPORT_EDGE_NORMAL_ELEMENT_NODE] = DATA_ITEM_SIZE;
  fixed_length[SUPPORT_EDGE_NORMAL_ELEMENT_NODE] = 0;
  data_class[SUPPORT_EDGE_NORMAL_ELEMENT_NODE] = FORCE;
  data_required[SUPPORT_EDGE_NORMAL_ELEMENT_NODE] = SUPPORT_EDGE_NORMAL;

  strcpy(name[SUPPORT_EDGE_NORMAL_NODE],"support_edge_normal_node");
  type[SUPPORT_EDGE_NORMAL_NODE] = INTEGER;
  data_length[SUPPORT_EDGE_NORMAL_NODE] = DATA_ITEM_SIZE;
  fixed_length[SUPPORT_EDGE_NORMAL_NODE] = 0;
  data_class[SUPPORT_EDGE_NORMAL_NODE] = FORCE;
  data_required[SUPPORT_EDGE_NORMAL_NODE] = SUPPORT_EDGE_NORMAL;

  strcpy(name[SUPPORT_EDGE_NORMAL_DAMPING],"support_edge_normal_damping");
  type[SUPPORT_EDGE_NORMAL_DAMPING] = DOUBLE_PRECISION;
  data_length[SUPPORT_EDGE_NORMAL_DAMPING] = 2;
  data_class[SUPPORT_EDGE_NORMAL_DAMPING] = FORCE;
  data_required[SUPPORT_EDGE_NORMAL_DAMPING] = SUPPORT_EDGE_NORMAL;

  strcpy(name[SUPPORT_EDGE_NORMAL_DAMPING_AUTOMATIC],"support_edge_normal_damping_automatic");
  type[SUPPORT_EDGE_NORMAL_DAMPING_AUTOMATIC] = INTEGER;
  data_length[SUPPORT_EDGE_NORMAL_DAMPING_AUTOMATIC] = 1;
  data_class[SUPPORT_EDGE_NORMAL_DAMPING_AUTOMATIC] = FORCE;
  data_required[SUPPORT_EDGE_NORMAL_DAMPING_AUTOMATIC] = SUPPORT_EDGE_NORMAL;

  strcpy(name[SUPPORT_EDGE_NORMAL_DAMPING_AUTOMATIC_APPARENT],"support_edge_normal_damping_automatic_apparent");
  type[SUPPORT_EDGE_NORMAL_DAMPING_AUTOMATIC_APPARENT] = INTEGER;
  data_length[SUPPORT_EDGE_NORMAL_DAMPING_AUTOMATIC_APPARENT] = 1;
  data_class[SUPPORT_EDGE_NORMAL_DAMPING_AUTOMATIC_APPARENT] = FORCE;
  data_required[SUPPORT_EDGE_NORMAL_DAMPING_AUTOMATIC_APPARENT] = SUPPORT_EDGE_NORMAL;

  strcpy(name[SUPPORT_EDGE_NORMAL_DENSITY],"support_edge_normal_density");
  type[SUPPORT_EDGE_NORMAL_DENSITY] = DOUBLE_PRECISION;
  data_length[SUPPORT_EDGE_NORMAL_DENSITY] = 2;
  data_class[SUPPORT_EDGE_NORMAL_DENSITY] = FORCE;
  data_required[SUPPORT_EDGE_NORMAL_DENSITY] = SUPPORT_EDGE_NORMAL;

  strcpy(name[SUPPORT_EDGE_NORMAL_FACTOR],"support_edge_normal_factor");
  type[SUPPORT_EDGE_NORMAL_FACTOR] = DOUBLE_PRECISION;
  data_length[SUPPORT_EDGE_NORMAL_FACTOR] = DATA_ITEM_SIZE;
  fixed_length[SUPPORT_EDGE_NORMAL_FACTOR] = 0;
  data_class[SUPPORT_EDGE_NORMAL_FACTOR] = FORCE;
  data_required[SUPPORT_EDGE_NORMAL_FACTOR] = SUPPORT_EDGE_NORMAL;

  strcpy(name[SUPPORT_EDGE_NORMAL_FORCE_INITIAL],"support_edge_normal_force_initial");
  type[SUPPORT_EDGE_NORMAL_FORCE_INITIAL] = DOUBLE_PRECISION;
  data_length[SUPPORT_EDGE_NORMAL_FORCE_INITIAL] = 2;
  data_class[SUPPORT_EDGE_NORMAL_FORCE_INITIAL] = FORCE;
  data_required[SUPPORT_EDGE_NORMAL_FORCE_INITIAL] = SUPPORT_EDGE_NORMAL;

  strcpy(name[SUPPORT_EDGE_NORMAL_TIME],"support_edge_normal_time");
  type[SUPPORT_EDGE_NORMAL_TIME] = DOUBLE_PRECISION;
  data_length[SUPPORT_EDGE_NORMAL_TIME] = DATA_ITEM_SIZE;
  fixed_length[SUPPORT_EDGE_NORMAL_TIME] = 0;
  data_class[SUPPORT_EDGE_NORMAL_TIME] = FORCE;
  data_required[SUPPORT_EDGE_NORMAL_TIME] = SUPPORT_EDGE_NORMAL;

  strcpy(name[NODE_SUPPORT_EDGE_NORMAL_FORCE],"node_support_edge_normal_force");
  type[NODE_SUPPORT_EDGE_NORMAL_FORCE] = DOUBLE_PRECISION;
  data_length[NODE_SUPPORT_EDGE_NORMAL_FORCE] = ndim;
  version_all[NODE_SUPPORT_EDGE_NORMAL_FORCE] = 1;
  data_class[NODE_SUPPORT_EDGE_NORMAL_FORCE] = NODE;
  data_required[NODE_SUPPORT_EDGE_NORMAL_FORCE] = NODE;


  strcpy(name[SUPPORT_EDGE_NORMAL_PLASTI_COMPRESSION],"support_edge_normal_plasti_compression");
  type[SUPPORT_EDGE_NORMAL_PLASTI_COMPRESSION] = DOUBLE_PRECISION;
  data_length[SUPPORT_EDGE_NORMAL_PLASTI_COMPRESSION] = 2;
  data_class[SUPPORT_EDGE_NORMAL_PLASTI_COMPRESSION] = FORCE;
  data_required[SUPPORT_EDGE_NORMAL_PLASTI_COMPRESSION] = SUPPORT_EDGE_NORMAL;

  strcpy(name[SUPPORT_EDGE_NORMAL_PLASTI_FRICTION],"support_edge_normal_plasti_friction");
  type[SUPPORT_EDGE_NORMAL_PLASTI_FRICTION] = DOUBLE_PRECISION;
  data_length[SUPPORT_EDGE_NORMAL_PLASTI_FRICTION] = 2;
  data_class[SUPPORT_EDGE_NORMAL_PLASTI_FRICTION] = FORCE;
  data_required[SUPPORT_EDGE_NORMAL_PLASTI_FRICTION] = SUPPORT_EDGE_NORMAL;

  strcpy(name[SUPPORT_EDGE_NORMAL_PLASTI_RESIDUAL_STIFFNESS],"support_edge_normal_plasti_residual_stiffness");
  type[SUPPORT_EDGE_NORMAL_PLASTI_RESIDUAL_STIFFNESS] = DOUBLE_PRECISION;
  data_length[SUPPORT_EDGE_NORMAL_PLASTI_RESIDUAL_STIFFNESS] = 1;
  data_class[SUPPORT_EDGE_NORMAL_PLASTI_RESIDUAL_STIFFNESS] = FORCE;
  data_required[SUPPORT_EDGE_NORMAL_PLASTI_RESIDUAL_STIFFNESS] = SUPPORT_EDGE_NORMAL;

  strcpy(name[SUPPORT_EDGE_NORMAL_PLASTI_TENSION],"support_edge_normal_plasti_tension");
  type[SUPPORT_EDGE_NORMAL_PLASTI_TENSION] = INTEGER;
  data_length[SUPPORT_EDGE_NORMAL_PLASTI_TENSION] = 1;
  data_class[SUPPORT_EDGE_NORMAL_PLASTI_TENSION] = FORCE;
  data_required[SUPPORT_EDGE_NORMAL_PLASTI_TENSION] = SUPPORT_EDGE_NORMAL;

  strcpy(name[SUPPORT_EDGE_NORMAL_PLASTI_TENSION_DOUBLE],"support_edge_normal_plasti_tension_double");
  type[SUPPORT_EDGE_NORMAL_PLASTI_TENSION_DOUBLE] = DOUBLE_PRECISION;
  data_length[SUPPORT_EDGE_NORMAL_PLASTI_TENSION_DOUBLE] = 1;
  data_class[SUPPORT_EDGE_NORMAL_PLASTI_TENSION_DOUBLE] = FORCE;
  data_required[SUPPORT_EDGE_NORMAL_PLASTI_TENSION_DOUBLE] = SUPPORT_EDGE_NORMAL;

  strcpy(name[NODE_SUPPORT_EDGE_NORMAL_PLASTI_TENSION_STATUS],"node_support_edge_normal_plasti_tension_status");
  type[NODE_SUPPORT_EDGE_NORMAL_PLASTI_TENSION_STATUS] = INTEGER;
  data_length[NODE_SUPPORT_EDGE_NORMAL_PLASTI_TENSION_STATUS] = 1;
  version_all[NODE_SUPPORT_EDGE_NORMAL_PLASTI_TENSION_STATUS] = 1;
  data_class[NODE_SUPPORT_EDGE_NORMAL_PLASTI_TENSION_STATUS] = NODE;
  data_required[NODE_SUPPORT_EDGE_NORMAL_PLASTI_TENSION_STATUS] = NODE;

  strcpy(name[CONTROL_SUPPORT_EDGE_NORMAL_DAMPING_APPLY],"control_support_edge_normal_damping_apply");
  type[CONTROL_SUPPORT_EDGE_NORMAL_DAMPING_APPLY] = INTEGER;
  data_length[CONTROL_SUPPORT_EDGE_NORMAL_DAMPING_APPLY] = 1;
  data_class[CONTROL_SUPPORT_EDGE_NORMAL_DAMPING_APPLY] = CONTROL;

  strcpy(name[CONTROL_SUPPORT_EDGE_NORMAL_STIFFNESS_FREEZE],"control_support_edge_normal_stiffness_freeze");
  type[CONTROL_SUPPORT_EDGE_NORMAL_STIFFNESS_FREEZE] = INTEGER;
  data_length[CONTROL_SUPPORT_EDGE_NORMAL_STIFFNESS_FREEZE] = 1;
  data_class[CONTROL_SUPPORT_EDGE_NORMAL_STIFFNESS_FREEZE] = CONTROL;

  // Sprint 13: Professional convergence backlog (CORPUS-BASELINE)
  strcpy(name[GROUP_BEAM_SHEAR],"group_beam_shear");
  type[GROUP_BEAM_SHEAR] = DOUBLE_PRECISION;
  data_length[GROUP_BEAM_SHEAR] = 1;
  data_class[GROUP_BEAM_SHEAR] = BEAM;

  // Sprint 13: post_element_force family (manual Professional 6.927-
  // 6.935): cross-section forces/moments from the element internal
  // nodal forces (the free-body statics of the L5 machinery)
  strcpy(name[POST_ELEMENT_FORCE],"post_element_force");
  type[POST_ELEMENT_FORCE] = DOUBLE_PRECISION;
  data_length[POST_ELEMENT_FORCE] = DATA_ITEM_SIZE;
  fixed_length[POST_ELEMENT_FORCE] = 0;

  strcpy(name[POST_ELEMENT_FORCE_GEOMETRY],"post_element_force_geometry");
  type[POST_ELEMENT_FORCE_GEOMETRY] = INTEGER;
  data_length[POST_ELEMENT_FORCE_GEOMETRY] = DATA_ITEM_SIZE;
  fixed_length[POST_ELEMENT_FORCE_GEOMETRY] = 0;
  data_required[POST_ELEMENT_FORCE_GEOMETRY] = POST_ELEMENT_FORCE;

  strcpy(name[POST_ELEMENT_FORCE_GROUP],"post_element_force_group");
  type[POST_ELEMENT_FORCE_GROUP] = INTEGER;
  data_length[POST_ELEMENT_FORCE_GROUP] = DATA_ITEM_SIZE;
  fixed_length[POST_ELEMENT_FORCE_GROUP] = 0;
  data_required[POST_ELEMENT_FORCE_GROUP] = POST_ELEMENT_FORCE;

  strcpy(name[POST_ELEMENT_FORCE_NUMBER],"post_element_force_number");
  type[POST_ELEMENT_FORCE_NUMBER] = INTEGER;
  data_length[POST_ELEMENT_FORCE_NUMBER] = DATA_ITEM_SIZE;
  fixed_length[POST_ELEMENT_FORCE_NUMBER] = 0;
  data_required[POST_ELEMENT_FORCE_NUMBER] = POST_ELEMENT_FORCE;

  strcpy(name[POST_ELEMENT_FORCE_NORMAL],"post_element_force_normal");
  type[POST_ELEMENT_FORCE_NORMAL] = INTEGER;
  data_length[POST_ELEMENT_FORCE_NORMAL] = 1;
  data_required[POST_ELEMENT_FORCE_NORMAL] = POST_ELEMENT_FORCE;

  strcpy(name[POST_ELEMENT_FORCE_FORCE],"post_element_force_force");
  type[POST_ELEMENT_FORCE_FORCE] = INTEGER;
  data_length[POST_ELEMENT_FORCE_FORCE] = 1;
  data_required[POST_ELEMENT_FORCE_FORCE] = POST_ELEMENT_FORCE;

  strcpy(name[POST_ELEMENT_FORCE_INERTIA],"post_element_force_inertia");
  type[POST_ELEMENT_FORCE_INERTIA] = INTEGER;
  data_length[POST_ELEMENT_FORCE_INERTIA] = 1;
  data_required[POST_ELEMENT_FORCE_INERTIA] = POST_ELEMENT_FORCE;

  strcpy(name[POST_ELEMENT_FORCE_MULTIPLY_FACTOR],"post_element_force_multiply_factor");
  type[POST_ELEMENT_FORCE_MULTIPLY_FACTOR] = DOUBLE_PRECISION;
  data_length[POST_ELEMENT_FORCE_MULTIPLY_FACTOR] = 1;
  data_required[POST_ELEMENT_FORCE_MULTIPLY_FACTOR] = POST_ELEMENT_FORCE;

  strcpy(name[POST_ELEMENT_FORCE_RESULT],"post_element_force_result");
  type[POST_ELEMENT_FORCE_RESULT] = DOUBLE_PRECISION;
  data_length[POST_ELEMENT_FORCE_RESULT] = 5;

  strcpy(name[PRINT_APPLY],"print_apply");
  type[PRINT_APPLY] = INTEGER;
  data_length[PRINT_APPLY] = 1;
  no_index[PRINT_APPLY] = 1;

  strcpy(name[VOLUME_FACTOR],"volume_factor");
  type[VOLUME_FACTOR] = DOUBLE_PRECISION;
  data_length[VOLUME_FACTOR] = DATA_ITEM_SIZE;
  fixed_length[VOLUME_FACTOR] = 0;
  no_index[VOLUME_FACTOR] = 1;
  data_class[VOLUME_FACTOR] = VOLUME;

  // Sprint 12 lot 3: the Professional's global solver_* family (manual
  // 6.1047-6.1056), timestep_predict_velocity / _iterations_automatic_
  // apply (6.1089-6.1090), tochnog_version (6.1091), volume_factor_x
  // (6.1094) and zip (6.1095). Plain (non-indexed) records, unlike
  // their control_solver_* counterparts of Sprint 9.
  strcpy(name[SOLVER],"solver");
  type[SOLVER] = INTEGER;
  data_length[SOLVER] = 1;
  no_index[SOLVER] = 1;

  strcpy(name[SOLVER_BICG_ERROR],"solver_bicg_error");
  type[SOLVER_BICG_ERROR] = DOUBLE_PRECISION;
  data_length[SOLVER_BICG_ERROR] = 1;
  no_index[SOLVER_BICG_ERROR] = 1;

  strcpy(name[SOLVER_BICG_RESTART],"solver_bicg_restart");
  type[SOLVER_BICG_RESTART] = INTEGER;
  data_length[SOLVER_BICG_RESTART] = 1;
  no_index[SOLVER_BICG_RESTART] = 1;

  strcpy(name[SOLVER_BICG_STOP],"solver_bicg_stop");
  type[SOLVER_BICG_STOP] = INTEGER;
  data_length[SOLVER_BICG_STOP] = 1;
  no_index[SOLVER_BICG_STOP] = 1;

  strcpy(name[SOLVER_MATRIX_SAVE],"solver_matrix_save");
  type[SOLVER_MATRIX_SAVE] = INTEGER;
  data_length[SOLVER_MATRIX_SAVE] = 1;
  no_index[SOLVER_MATRIX_SAVE] = 1;

  strcpy(name[SOLVER_MATRIX_SYMMETRIC],"solver_matrix_symmetric");
  type[SOLVER_MATRIX_SYMMETRIC] = INTEGER;
  data_length[SOLVER_MATRIX_SYMMETRIC] = 1;
  no_index[SOLVER_MATRIX_SYMMETRIC] = 1;

  strcpy(name[SOLVER_PARDISO_ORDERING],"solver_pardiso_ordering");
  type[SOLVER_PARDISO_ORDERING] = INTEGER;
  data_length[SOLVER_PARDISO_ORDERING] = 1;
  no_index[SOLVER_PARDISO_ORDERING] = 1;

  strcpy(name[SOLVER_PARDISO_OUT_OF_CORE],"solver_pardiso_out_of_core");
  type[SOLVER_PARDISO_OUT_OF_CORE] = INTEGER;
  data_length[SOLVER_PARDISO_OUT_OF_CORE] = 1;
  no_index[SOLVER_PARDISO_OUT_OF_CORE] = 1;

  strcpy(name[SOLVER_PARDISO_PROCESSORS],"solver_pardiso_processors");
  type[SOLVER_PARDISO_PROCESSORS] = INTEGER;
  data_length[SOLVER_PARDISO_PROCESSORS] = 1;
  no_index[SOLVER_PARDISO_PROCESSORS] = 1;

  strcpy(name[SOLVER_PARDISO_PROCESSORS_MAXIMUM],"solver_pardiso_processors_maximum");
  type[SOLVER_PARDISO_PROCESSORS_MAXIMUM] = INTEGER;
  data_length[SOLVER_PARDISO_PROCESSORS_MAXIMUM] = 1;
  no_index[SOLVER_PARDISO_PROCESSORS_MAXIMUM] = 1;

  strcpy(name[TIMESTEP_ITERATIONS_AUTOMATIC_APPLY],"timestep_iterations_automatic_apply");
  type[TIMESTEP_ITERATIONS_AUTOMATIC_APPLY] = INTEGER;
  data_length[TIMESTEP_ITERATIONS_AUTOMATIC_APPLY] = 1;
  no_index[TIMESTEP_ITERATIONS_AUTOMATIC_APPLY] = 1;

  strcpy(name[TIMESTEP_PREDICT_VELOCITY],"timestep_predict_velocity");
  type[TIMESTEP_PREDICT_VELOCITY] = INTEGER;
  data_length[TIMESTEP_PREDICT_VELOCITY] = 1;
  no_index[TIMESTEP_PREDICT_VELOCITY] = 1;

  strcpy(name[TOCHNOG_VERSION],"tochnog_version");
  type[TOCHNOG_VERSION] = INTEGER;
  data_length[TOCHNOG_VERSION] = 3;

  strcpy(name[VOLUME_FACTOR_X],"volume_factor_x");
  type[VOLUME_FACTOR_X] = DOUBLE_PRECISION;
  data_length[VOLUME_FACTOR_X] = DATA_ITEM_SIZE;
  fixed_length[VOLUME_FACTOR_X] = 0;
  no_index[VOLUME_FACTOR_X] = 1;
  data_class[VOLUME_FACTOR_X] = VOLUME;

  strcpy(name[ZIP],"zip");
  type[ZIP] = INTEGER;
  data_length[ZIP] = 1;
  no_index[ZIP] = 1;

  strcpy(name[WAVE],"wave");

  strcpy(name[WAVE_SCALAR],"wave_scalar");

  strcpy(name[WAVE_FSCALAR],"wave_fscalar");

  // --- corpus parser keywords (registered for the parser only; the
  //     consumption of these families is PENDING, see
  //     SEGUIMIENTO-CONVERGENCIA.md) ---

  // post_calcul_length (written by the Professional in the .dbs and
  // used as input by the corpus tests): length of the post_calcul
  // record. The GNU derives the length from the parsed values.
  // bounda_used (written by the Professional in the .dbs database):
  // switch whether the bounda record with the index was applied in the
  // step. Output record for .dbs parity.
  strcpy(name[BOUNDA_USED],"bounda_used");
  type[BOUNDA_USED] = INTEGER;
  data_length[BOUNDA_USED] = 1;
  data_class[BOUNDA_USED] = BOUNDA;

  strcpy(name[POST_CALCUL_LENGTH],"post_calcul_length");
  type[POST_CALCUL_LENGTH] = INTEGER;
  data_length[POST_CALCUL_LENGTH] = 1;
  no_index[POST_CALCUL_LENGTH] = 1;
  data_class[POST_CALCUL_LENGTH] = POST;

  // family mpc_* (manual Professional 6.856-6.875): multi point
  // constraints. Consumption PENDING.
  strcpy(name[MPC_ELEMENT_GROUP],"mpc_element_group");
  type[MPC_ELEMENT_GROUP] = INTEGER;
  data_length[MPC_ELEMENT_GROUP] = DATA_ITEM_SIZE;
  fixed_length[MPC_ELEMENT_GROUP] = 0;
  data_class[MPC_ELEMENT_GROUP] = CONTROL;

  strcpy(name[MPC_ELEMENT_GROUP_DOF],"mpc_element_group_dof");
  type[MPC_ELEMENT_GROUP_DOF] = INTEGER;
  data_length[MPC_ELEMENT_GROUP_DOF] = DATA_ITEM_SIZE;
  fixed_length[MPC_ELEMENT_GROUP_DOF] = 0;
  data_class[MPC_ELEMENT_GROUP_DOF] = CONTROL;

  strcpy(name[MPC_ELEMENT_GROUP_GEOMETRY],"mpc_element_group_geometry");
  type[MPC_ELEMENT_GROUP_GEOMETRY] = INTEGER;
  data_length[MPC_ELEMENT_GROUP_GEOMETRY] = DATA_ITEM_SIZE;
  fixed_length[MPC_ELEMENT_GROUP_GEOMETRY] = 0;
  data_class[MPC_ELEMENT_GROUP_GEOMETRY] = CONTROL;

  strcpy(name[MPC_APPLY],"mpc_apply");
  type[MPC_APPLY] = INTEGER;
  data_length[MPC_APPLY] = 1;
  no_index[MPC_APPLY] = 1;
  data_class[MPC_APPLY] = CONTROL;

  strcpy(name[CONTROL_MPC_APPLY],"control_mpc_apply");
  type[CONTROL_MPC_APPLY] = INTEGER;
  data_length[CONTROL_MPC_APPLY] = 1;
  data_class[CONTROL_MPC_APPLY] = CONTROL;

  strcpy(name[MPC_GEOMETRY],"mpc_geometry");
  type[MPC_GEOMETRY] = INTEGER;
  data_length[MPC_GEOMETRY] = DATA_ITEM_SIZE;
  fixed_length[MPC_GEOMETRY] = 0;
  data_class[MPC_GEOMETRY] = CONTROL;

  strcpy(name[MPC_GEOMETRY_DOF],"mpc_geometry_dof");
  type[MPC_GEOMETRY_DOF] = INTEGER;
  data_length[MPC_GEOMETRY_DOF] = DATA_ITEM_SIZE;
  fixed_length[MPC_GEOMETRY_DOF] = 0;
  data_class[MPC_GEOMETRY_DOF] = CONTROL;

  strcpy(name[MPC_GEOMETRY_SWITCH],"mpc_geometry_switch");
  type[MPC_GEOMETRY_SWITCH] = INTEGER;
  data_length[MPC_GEOMETRY_SWITCH] = DATA_ITEM_SIZE;
  fixed_length[MPC_GEOMETRY_SWITCH] = 0;
  data_class[MPC_GEOMETRY_SWITCH] = CONTROL;

  strcpy(name[MPC_LINEAR_QUADRATIC],"mpc_linear_quadratic");
  type[MPC_LINEAR_QUADRATIC] = INTEGER;
  data_length[MPC_LINEAR_QUADRATIC] = 1;
  no_index[MPC_LINEAR_QUADRATIC] = 1;
  data_class[MPC_LINEAR_QUADRATIC] = CONTROL;

  strcpy(name[MPC_NODE_FACTOR],"mpc_node_factor");
  type[MPC_NODE_FACTOR] = DOUBLE_PRECISION;
  data_length[MPC_NODE_FACTOR] = DATA_ITEM_SIZE;
  fixed_length[MPC_NODE_FACTOR] = 0;
  data_class[MPC_NODE_FACTOR] = CONTROL;

  strcpy(name[MPC_NODE_NUMBER],"mpc_node_number");
  type[MPC_NODE_NUMBER] = INTEGER;
  data_length[MPC_NODE_NUMBER] = DATA_ITEM_SIZE;
  fixed_length[MPC_NODE_NUMBER] = 0;
  data_class[MPC_NODE_NUMBER] = CONTROL;

  // mpc_linear_quadratic_mesh_fingerprint: INTERNAL record (external 0,
  // never written to the .dbs) that stores the mesh fingerprint of the
  // last mpc_linear_quadratic tie generation plus the index range of the
  // generated mpc_node_number/mpc_node_factor records:
  // [fingerprint, start_index, count].
  strcpy(name[MPC_LINEAR_QUADRATIC_MESH_FINGERPRINT],
    "mpc_linear_quadratic_mesh_fingerprint");
  type[MPC_LINEAR_QUADRATIC_MESH_FINGERPRINT] = INTEGER;
  data_length[MPC_LINEAR_QUADRATIC_MESH_FINGERPRINT] = 3;
  fixed_length[MPC_LINEAR_QUADRATIC_MESH_FINGERPRINT] = 3;
  external[MPC_LINEAR_QUADRATIC_MESH_FINGERPRINT] = 0;
  no_index[MPC_LINEAR_QUADRATIC_MESH_FINGERPRINT] = 1;
  data_class[MPC_LINEAR_QUADRATIC_MESH_FINGERPRINT] = CONTROL;

  // mpc_element_group_mesh_fingerprint: INTERNAL bookkeeping record of
  // the mpc_element_group tie generation (same layout as the
  // mpc_linear_quadratic one): [mesh fingerprint, start_index, count].
  strcpy(name[MPC_ELEMENT_GROUP_MESH_FINGERPRINT],
    "mpc_element_group_mesh_fingerprint");
  type[MPC_ELEMENT_GROUP_MESH_FINGERPRINT] = INTEGER;
  data_length[MPC_ELEMENT_GROUP_MESH_FINGERPRINT] = 3;
  fixed_length[MPC_ELEMENT_GROUP_MESH_FINGERPRINT] = 3;
  external[MPC_ELEMENT_GROUP_MESH_FINGERPRINT] = 0;
  no_index[MPC_ELEMENT_GROUP_MESH_FINGERPRINT] = 1;
  data_class[MPC_ELEMENT_GROUP_MESH_FINGERPRINT] = CONTROL;

  // -veln family (manual Professional 6.22): VELN is the keyword a
  // bounda_dof/bounda_unknown record lists as dof to prescribe the
  // zero velocity NORMAL to a plane (geometry normal or the bounda_normal
  // of the same index). MPC_FROM_BOUNDA marks the mpc_node_number/
  // mpc_node_factor records that the -veln bounda generated (index k =
  // -yes beside the generated record k; the Professional .dbs stores the
  // same marker). BOUNDA_VELN_MESH_FINGERPRINT is the INTERNAL
  // bookkeeping record of the generator (indexed by bounda record,
  // [mesh fingerprint, start_index, count] - the mpc_linear_quadratic
  // pattern). See bounda.cc::bounda_veln_mpc().
  strcpy(name[VELN],"veln");

  strcpy(name[MPC_FROM_BOUNDA],"mpc_from_bounda");
  type[MPC_FROM_BOUNDA] = INTEGER;
  data_length[MPC_FROM_BOUNDA] = 1;
  fixed_length[MPC_FROM_BOUNDA] = 0;
  data_class[MPC_FROM_BOUNDA] = CONTROL;

  strcpy(name[BOUNDA_VELN_MESH_FINGERPRINT],"bounda_veln_mesh_fingerprint");
  type[BOUNDA_VELN_MESH_FINGERPRINT] = INTEGER;
  data_length[BOUNDA_VELN_MESH_FINGERPRINT] = 3;
  fixed_length[BOUNDA_VELN_MESH_FINGERPRINT] = 3;
  external[BOUNDA_VELN_MESH_FINGERPRINT] = 0;
  data_class[BOUNDA_VELN_MESH_FINGERPRINT] = CONTROL;

  // control_mesh_truss_distribute_mpc (manual Professional 6.245) and
  // the _exact variant (6.250): distribute truss nodes over the
  // isoparametric elements with multi point constraints. Consumption
  // PENDING.
  strcpy(name[CONTROL_MESH_TRUSS_DISTRIBUTE_MPC],"control_mesh_truss_distribute_mpc");
  type[CONTROL_MESH_TRUSS_DISTRIBUTE_MPC] = INTEGER;
  data_length[CONTROL_MESH_TRUSS_DISTRIBUTE_MPC] = 1;
  fixed_length[CONTROL_MESH_TRUSS_DISTRIBUTE_MPC] = 1;
  data_class[CONTROL_MESH_TRUSS_DISTRIBUTE_MPC] = CONTROL;

  strcpy(name[CONTROL_MESH_TRUSS_DISTRIBUTE_MPC_EXACT],"control_mesh_truss_distribute_mpc_exact");
  type[CONTROL_MESH_TRUSS_DISTRIBUTE_MPC_EXACT] = INTEGER;
  data_length[CONTROL_MESH_TRUSS_DISTRIBUTE_MPC_EXACT] = 1;
  fixed_length[CONTROL_MESH_TRUSS_DISTRIBUTE_MPC_EXACT] = 1;
  data_class[CONTROL_MESH_TRUSS_DISTRIBUTE_MPC_EXACT] = CONTROL;

  // control_mesh_truss_distribute_mpc_element_group_truss (manual
  // Professional 6.248): the truss element groups whose nodes are
  // coupled to the isoparametric elements. Consumption PENDING
  // (same family as the parent record above; registered so the truss
  // distribute corpus tests parse - the target physics of truss11/12
  // is unaffected: all nodes are prescribed/fixed there).
  strcpy(name[CONTROL_MESH_TRUSS_DISTRIBUTE_MPC_ELEMENT_GROUP_TRUSS],
    "control_mesh_truss_distribute_mpc_element_group_truss");
  type[CONTROL_MESH_TRUSS_DISTRIBUTE_MPC_ELEMENT_GROUP_TRUSS] = INTEGER;
  data_length[CONTROL_MESH_TRUSS_DISTRIBUTE_MPC_ELEMENT_GROUP_TRUSS] = 1;
  fixed_length[CONTROL_MESH_TRUSS_DISTRIBUTE_MPC_ELEMENT_GROUP_TRUSS] = 1;
  data_class[CONTROL_MESH_TRUSS_DISTRIBUTE_MPC_ELEMENT_GROUP_TRUSS] =
    CONTROL;

  // strain_volume_* (manual Professional 6.96x): prescribed volume
  // strain and its element. Consumption PENDING.
  strcpy(name[STRAIN_VOLUME_ABSOLUTE_TIME],"strain_volume_absolute_time");
  type[STRAIN_VOLUME_ABSOLUTE_TIME] = DOUBLE_PRECISION;
  data_length[STRAIN_VOLUME_ABSOLUTE_TIME] = DATA_ITEM_SIZE;
  fixed_length[STRAIN_VOLUME_ABSOLUTE_TIME] = 0;
  data_class[STRAIN_VOLUME_ABSOLUTE_TIME] = CONTROL;

  strcpy(name[STRAIN_VOLUME_ELEMENT],"strain_volume_element");
  type[STRAIN_VOLUME_ELEMENT] = INTEGER;
  data_length[STRAIN_VOLUME_ELEMENT] = 1;
  fixed_length[STRAIN_VOLUME_ELEMENT] = 1;
  data_class[STRAIN_VOLUME_ELEMENT] = CONTROL;

  // strain_volume output items (manual Professional 6.966/6.968):
  // consumed by the target_item records of the corpus tests.
  // Consumption PENDING.
  strcpy(name[POST_STRAIN_VOLUME_ABSOLUTE],"post_strain_volume_absolute");
  type[POST_STRAIN_VOLUME_ABSOLUTE] = DOUBLE_PRECISION;
  data_length[POST_STRAIN_VOLUME_ABSOLUTE] = 1;
  data_class[POST_STRAIN_VOLUME_ABSOLUTE] = POST;

  strcpy(name[POST_STRAIN_VOLUME_RELATIVE],"post_strain_volume_relative");
  type[POST_STRAIN_VOLUME_RELATIVE] = DOUBLE_PRECISION;
  data_length[POST_STRAIN_VOLUME_RELATIVE] = 1;
  data_class[POST_STRAIN_VOLUME_RELATIVE] = POST;

  strcpy(name[YES],"yes");

  strcpy(name[X],"x");

  strcpy(name[Y],"y");

  strcpy(name[Z],"z");

  strcpy(name[ANGLE],"angle");

  // PRISM15 element type (DEV-B prism15 sprint): the 15-node quadratic
  // prism of the Professional (corpus prism15.dat). Registered at the
  // end like the appended enum value (see tochnog.h): only the name is
  // needed, the ELEMENT record machinery is name-generic.
  strcpy(name[PRISM15],"prism15");

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
    else if ( dof_type[iuknwn]==-GROUNDFLOW_SATURATION ) 
      strcpy( basename, "gsat" );
    else if ( dof_type[iuknwn]==-GROUNDFLOW_PRESSURE_GRADIENT ) {
      // manual Professional 4.7: the gradient of the hydraulic pressure
      // head dh/dx, dh/dy, dh/dz added to the node_dof records
      // (Professional basenames pres_gradx/pres_grady/pres_gradz)
      if ( iuknwn==pres_grad_indx ) n = 0;
      n++;
      strcpy( basename, "pres_grad" );
      if      ( n==1 ) strcat( basename, "x" );
      else if ( n==2 ) strcat( basename, "y" );
      else if ( n==3 ) strcat( basename, "z" );
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
    else if ( dof_type[iuknwn]==-MATERI_DISPLACEMENT_RELATIVE ) {
      if ( iuknwn==dis_rel_indx ) n = 0;
      n++;
      if      ( n==1 ) strcpy( basename, "disrx"  );
      else if ( n==2 ) strcpy( basename, "disry"  );
      else if ( n==3 ) strcpy( basename, "disrz"  );
    }
    else if ( dof_type[iuknwn]==-MATERI_HISTORY_VARIABLES ) {
      if ( iuknwn==hisv_indx ) n = 0;
      // manual Professional 4.18: materi_plasti_diprisco_history declares
      // n history variables named dipriscohis0..dipriscohis(n-1) in the
      // node_dof records (the generic name is hisv<n>)
      if ( materi_plasti_diprisco_history )
        strcpy( basename, "dipriscohis" );
      else
        strcpy( basename, "hisv" );
      long_to_a( n, str );
      strcat( basename, str );
      n++;
    }
    else if ( dof_type[iuknwn]==-MATERI_PLASTI_HYPO_HISTORY ) {
      // manual Professional 4.23: hyhis0..hyhis7 (void ratio e,
      // substep size, mobilized friction angle, stiffness measure,
      // structure s, OCR, density index, intergranular rho)
      if ( iuknwn==hisv_indx ) n = 0;
      strcpy( basename, "hyhis" );
      long_to_a( n, str );
      strcat( basename, str );
      n++;
    }
    else if ( dof_type[iuknwn]==-MATERI_PLASTI_CAMCLAY_HISTORY ) {
      // manual Professional 4.16: cchis0 = void ratio e0 and
      // cchis1 = preconsolidation pressure p0 of the camclay model
      if ( iuknwn==hisv_indx ) n = 0;
      strcpy( basename, "cchis" );
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
    else if ( dof_type[iuknwn]==-MATERI_PLASTI_KAPPA_SHEAR )
      strcpy( basename, "kapsh" );
    else if ( dof_type[iuknwn]==-MATERI_PLASTI_CAP1_HISTORY )
      strcpy( basename, "pc" );
    else if ( dof_type[iuknwn]==-MATERI_PLASTI_HARDSOIL_HISTORY )
      // maximum |p| history of the hardsoil model (manual Professional
      // 4.22): same concept as materi_stress_pressure_history (4.50),
      // same shared dof (sph_indx), same basename
      strcpy( basename, "sph" );
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
    else if ( dof_type[iuknwn]==-MATERI_STRAIN_PLASTI_CAP ) {
      // cap plastic strain (manual Professional 4.35): same
      // 6-component layout as materi_strain_plasti, dedicated dof
      if ( iuknwn==capepp_indx ) n = 0;
      n++;
      if      ( n==1 ) strcpy( basename, "eppcapxx" );
      else if ( n==2 ) strcpy( basename, "eppcapxy" );
      else if ( n==3 ) strcpy( basename, "eppcapxz" );
      else if ( n==4 ) strcpy( basename, "eppcappyy" );
      else if ( n==5 ) strcpy( basename, "eppcapyz" );
      else if ( n==6 ) strcpy( basename, "eppcapzz" );
    }
    else if ( dof_type[iuknwn]==-MATERI_STRAIN_PLASTI_COMPRESSION ) {
      // compression plastic strain (manual Professional 4.36)
      if ( iuknwn==cepp_indx ) n = 0;
      n++;
      if      ( n==1 ) strcpy( basename, "eppcmpxx" );
      else if ( n==2 ) strcpy( basename, "eppcmpxy" );
      else if ( n==3 ) strcpy( basename, "eppcmpxz" );
      else if ( n==4 ) strcpy( basename, "eppcmpyy" );
      else if ( n==5 ) strcpy( basename, "eppcmpyz" );
      else if ( n==6 ) strcpy( basename, "eppcmpzz" );
    }
    else if ( dof_type[iuknwn]==-MATERI_STRAIN_PLASTI_DIPRISCO ) {
      // di Prisco plastic strain (manual Professional 4.37)
      if ( iuknwn==depp_indx ) n = 0;
      n++;
      if      ( n==1 ) strcpy( basename, "eppdipxx" );
      else if ( n==2 ) strcpy( basename, "eppdipxy" );
      else if ( n==3 ) strcpy( basename, "eppdipxz" );
      else if ( n==4 ) strcpy( basename, "eppdipyy" );
      else if ( n==5 ) strcpy( basename, "eppdipyz" );
      else if ( n==6 ) strcpy( basename, "eppdipzz" );
    }
    else if ( dof_type[iuknwn]==-MATERI_STRAIN_PLASTI_DRUCKPRAG ) {
      // drucker-prager plastic strain (manual Professional 4.39)
      if ( iuknwn==dpepp_indx ) n = 0;
      n++;
      if      ( n==1 ) strcpy( basename, "eppdrpxx" );
      else if ( n==2 ) strcpy( basename, "eppdrpxy" );
      else if ( n==3 ) strcpy( basename, "eppdrpxz" );
      else if ( n==4 ) strcpy( basename, "eppdrpyy" );
      else if ( n==5 ) strcpy( basename, "eppdrpyz" );
      else if ( n==6 ) strcpy( basename, "eppdrpzz" );
    }
    else if ( dof_type[iuknwn]==-MATERI_STRAIN_PLASTI_HARDSOIL ) {
      // hardsoil plastic strain (manual Professional 4.40): same
      // 6-component layout as materi_strain_plasti, dedicated dof
      if ( iuknwn==hsepp_indx ) n = 0;
      n++;
      if      ( n==1 ) strcpy( basename, "epphsxx" );
      else if ( n==2 ) strcpy( basename, "epphsxy" );
      else if ( n==3 ) strcpy( basename, "epphsxz" );
      else if ( n==4 ) strcpy( basename, "epphsyy" );
      else if ( n==5 ) strcpy( basename, "epphsyz" );
      else if ( n==6 ) strcpy( basename, "epphszz" );
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
    else if ( dof_type[iuknwn]==-MATERI_STRESS_PRESSURE_HISTORY )
      strcpy( basename, "sph" );
    else if ( dof_type[iuknwn]==-MATERI_ACCELERATION ) {
      if ( iuknwn==acc_indx ) n = 0;
      n++;
      if      ( n==1 ) strcpy( basename, "accx"  );
      else if ( n==2 ) strcpy( basename, "accy"  );
      else if ( n==3 ) strcpy( basename, "accz"  );
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

  // ------------------------------------------------------------------
  // KEYWORDS BATCH C (2026-09-04): small keyword clusters blocking
  // corpus tests. All records appended at the END of the enum (see
  // tochnog.h); registration blocks below stay at the end of
  // db_initialize so the enumeration order is irrelevant.
  // ------------------------------------------------------------------

  // spring memory/nonlinear-stiffness records (family group_type
  // springs; spring1/spring6 of the corpus). group_spring_memory:
  // memory model -updated_linear/-total_linear/-updated (manual
  // Professional 6.765); parse-consumed by the spring element only for
  // the diagram variant, the memory model itself is accepted (the GNU
  // spring law is incremental on the current configuration, which
  // equals -updated_linear; -total_linear coincides in the 1D corpus
  // tests).
  strcpy(name[GROUP_SPRING_MEMORY],"group_spring_memory");
  type[GROUP_SPRING_MEMORY] = INTEGER;
  data_length[GROUP_SPRING_MEMORY] = 1;
  data_class[GROUP_SPRING_MEMORY] = SPRING;
  data_required[GROUP_SPRING_MEMORY] = GROUP_TYPE;

  // group_spring_stiffness_nonlinear (manual Professional 6.768):
  // diagram epsilon0 k0 epsilon1 k1 ... of the spring stiffness vs the
  // total spring elongation (strain). Consumed in spring.cc when the
  // linear group_spring_stiffness record is absent.
  strcpy(name[GROUP_SPRING_STIFFNESS_NONLINEAR],"group_spring_stiffness_nonlinear");
  type[GROUP_SPRING_STIFFNESS_NONLINEAR] = DOUBLE_PRECISION;
  data_length[GROUP_SPRING_STIFFNESS_NONLINEAR] = DATA_ITEM_SIZE;
  fixed_length[GROUP_SPRING_STIFFNESS_NONLINEAR] = 0;
  data_class[GROUP_SPRING_STIFFNESS_NONLINEAR] = SPRING;
  data_required[GROUP_SPRING_STIFFNESS_NONLINEAR] = GROUP_TYPE;

  // element_spring_strain (manual Professional spring family): total
  // spring strain (= total elongation of the spring) written per
  // element, sibling of element_spring_force.
  strcpy(name[ELEMENT_SPRING_STRAIN],"element_spring_strain");
  type[ELEMENT_SPRING_STRAIN] = DOUBLE_PRECISION;
  data_length[ELEMENT_SPRING_STRAIN] = 1;
  version_all[ELEMENT_SPRING_STRAIN] = 1;
  print_only[ELEMENT_SPRING_STRAIN] = 1;
  data_class[ELEMENT_SPRING_STRAIN] = ELEMENT;
  data_required[ELEMENT_SPRING_STRAIN] = ELEMENT;

  // dependency_apply / control_dependency_apply (manual Professional
  // 6.126/6.125): global/per-timestep switches that enable or disable
  // the dependency_item/dependency_diagram machinery (get_group_data in
  // group.cc). Precedence: control_* overrides the global record;
  // default -yes.
  strcpy(name[DEPENDENCY_APPLY],"dependency_apply");
  type[DEPENDENCY_APPLY] = INTEGER;
  data_length[DEPENDENCY_APPLY] = 1;
  no_index[DEPENDENCY_APPLY] = 1;
  data_class[DEPENDENCY_APPLY] = DEPENDENCY;

  strcpy(name[CONTROL_DEPENDENCY_APPLY],"control_dependency_apply");
  type[CONTROL_DEPENDENCY_APPLY] = INTEGER;
  data_length[CONTROL_DEPENDENCY_APPLY] = 1;
  data_class[CONTROL_DEPENDENCY_APPLY] = CONTROL;

  // geometry_element_group (manual Professional 6.524/6.525): restrict
  // the geometry record with the same index to nodes that are also a
  // node of elements of one of the specified element groups. Method
  // (geometry_element_group_method) -all/-any/-only; the GNU consumes
  // the filter inside geometry() for node tests (see geometry.cc).
  strcpy(name[GEOMETRY_ELEMENT_GROUP],"geometry_element_group");
  type[GEOMETRY_ELEMENT_GROUP] = INTEGER;
  data_length[GEOMETRY_ELEMENT_GROUP] = DATA_ITEM_SIZE;
  fixed_length[GEOMETRY_ELEMENT_GROUP] = 0;
  data_class[GEOMETRY_ELEMENT_GROUP] = GEOMETRY;

  strcpy(name[GEOMETRY_ELEMENT_GROUP_METHOD],"geometry_element_group_method");
  type[GEOMETRY_ELEMENT_GROUP_METHOD] = INTEGER;
  data_length[GEOMETRY_ELEMENT_GROUP_METHOD] = 1;
  data_class[GEOMETRY_ELEMENT_GROUP_METHOD] = GEOMETRY;
  data_required[GEOMETRY_ELEMENT_GROUP_METHOD] = GEOMETRY_ELEMENT_GROUP;

  // group_materi_plasti_element_group/_factor (manual Professional
  // 6.698/6.699): model frictional slip of granular materials on other
  // materials (concrete/steel): phi/c/phiflow (and the equivalents of
  // the other plasticity models) of the granular element_group are
  // reduced with a factor (default 2./3., group_materi_plasti_element_
  // group_factor overrides per neighbor group) for granular elements
  // which are a DIRECT NEIGHBOR of an element of one of the listed
  // groups group_0 group_1 ...
  strcpy(name[GROUP_MATERI_PLASTI_ELEMENT_GROUP],"group_materi_plasti_element_group");
  type[GROUP_MATERI_PLASTI_ELEMENT_GROUP] = INTEGER;
  data_length[GROUP_MATERI_PLASTI_ELEMENT_GROUP] = DATA_ITEM_SIZE;
  fixed_length[GROUP_MATERI_PLASTI_ELEMENT_GROUP] = 0;
  data_class[GROUP_MATERI_PLASTI_ELEMENT_GROUP] = MATERI;
  data_required[GROUP_MATERI_PLASTI_ELEMENT_GROUP] = GROUP_TYPE;

  strcpy(name[GROUP_MATERI_PLASTI_ELEMENT_GROUP_FACTOR],"group_materi_plasti_element_group_factor");
  type[GROUP_MATERI_PLASTI_ELEMENT_GROUP_FACTOR] = DOUBLE_PRECISION;
  data_length[GROUP_MATERI_PLASTI_ELEMENT_GROUP_FACTOR] = DATA_ITEM_SIZE;
  fixed_length[GROUP_MATERI_PLASTI_ELEMENT_GROUP_FACTOR] = 0;
  data_class[GROUP_MATERI_PLASTI_ELEMENT_GROUP_FACTOR] = MATERI;
  data_required[GROUP_MATERI_PLASTI_ELEMENT_GROUP_FACTOR] = GROUP_MATERI_PLASTI_ELEMENT_GROUP;

  // post_calcul_static_pressure_height (manual Professional 6.921/6.922):
  // determine the static pressure (post_calcul -groundflow_pressure
  // -static_pressure) relative to a reference height instead of a
  // groundwater level. Record values: region triples (coord_min,
  // coord_max, height_ref) along the vertical coordinate; the
  // post_calcul_static_pressure_height_element_group record restricts
  // each region to an element group (-all = any group).
  strcpy(name[POST_CALCUL_STATIC_PRESSURE_HEIGHT],"post_calcul_static_pressure_height");
  type[POST_CALCUL_STATIC_PRESSURE_HEIGHT] = DOUBLE_PRECISION;
  data_length[POST_CALCUL_STATIC_PRESSURE_HEIGHT] = DATA_ITEM_SIZE;
  fixed_length[POST_CALCUL_STATIC_PRESSURE_HEIGHT] = 0;
  no_index[POST_CALCUL_STATIC_PRESSURE_HEIGHT] = 1;
  data_class[POST_CALCUL_STATIC_PRESSURE_HEIGHT] = POST;

  strcpy(name[POST_CALCUL_STATIC_PRESSURE_HEIGHT_ELEMENT_GROUP],"post_calcul_static_pressure_height_element_group");
  type[POST_CALCUL_STATIC_PRESSURE_HEIGHT_ELEMENT_GROUP] = INTEGER;
  data_length[POST_CALCUL_STATIC_PRESSURE_HEIGHT_ELEMENT_GROUP] = DATA_ITEM_SIZE;
  fixed_length[POST_CALCUL_STATIC_PRESSURE_HEIGHT_ELEMENT_GROUP] = 0;
  no_index[POST_CALCUL_STATIC_PRESSURE_HEIGHT_ELEMENT_GROUP] = 1;
  data_class[POST_CALCUL_STATIC_PRESSURE_HEIGHT_ELEMENT_GROUP] = POST;
  data_required[POST_CALCUL_STATIC_PRESSURE_HEIGHT_ELEMENT_GROUP] =
    POST_CALCUL_STATIC_PRESSURE_HEIGHT;

  // post_calcul_safety_method (manual Professional 6.919): how the
  // hydraulic piping/lifting safety factors of post_calcul
  // -materi_stress -safety_piping/-safety_lifting are determined:
  // -vertical (default, one value) / -prival (three values, principal
  // stresses) / -global (three values, global normal stresses).
  // post_calcul_safety_maximum (6.918) caps the factor.
  strcpy(name[POST_CALCUL_SAFETY_METHOD],"post_calcul_safety_method");
  type[POST_CALCUL_SAFETY_METHOD] = INTEGER;
  data_length[POST_CALCUL_SAFETY_METHOD] = 1;
  no_index[POST_CALCUL_SAFETY_METHOD] = 1;
  data_class[POST_CALCUL_SAFETY_METHOD] = POST;

  strcpy(name[POST_CALCUL_SAFETY_MAXIMUM],"post_calcul_safety_maximum");
  type[POST_CALCUL_SAFETY_MAXIMUM] = DOUBLE_PRECISION;
  data_length[POST_CALCUL_SAFETY_MAXIMUM] = 1;
  no_index[POST_CALCUL_SAFETY_MAXIMUM] = 1;
  data_class[POST_CALCUL_SAFETY_MAXIMUM] = POST;

  // post_calcul operator values -safety_piping/-safety_lifting (manual
  // Professional 6.9xx post_calcul command): the operators resolve as
  // INTEGER data items (like -force/-total_pressure) so the input
  // parser accepts them in the post_calcul record. Pure name entries.
  strcpy(name[SAFETY_PIPING],"safety_piping");
  type[SAFETY_PIPING] = INTEGER;
  data_length[SAFETY_PIPING] = 1;

  strcpy(name[SAFETY_LIFTING],"safety_lifting");
  type[SAFETY_LIFTING] = INTEGER;
  data_length[SAFETY_LIFTING] = 1;

  // safety method selector values -vertical/-global (pure name entries,
  // the same resolution pattern as -prival; manual Professional 6.919).
  strcpy(name[VERTICAL],"vertical");
  type[VERTICAL] = INTEGER;
  data_length[VERTICAL] = 1;

  strcpy(name[GLOBAL],"global");
  type[GLOBAL] = INTEGER;
  data_length[GLOBAL] = 1;

  // item names of the hydraulic safety factors in the post_calcul_label
  // of the .dbs (manual Professional 6.919): with method -prival the
  // names are <safety>_prival_0..2 and with -global <safety>_global_x/y/z;
  // the -vertical method reuses the plain names safety_piping/
  // safety_lifting (one value). Pure name entries so target_item
  // records parse and resolve the generated slots (same resolution
  // pattern as -to_pres/-st_pres/-dy_pres).
  strcpy(name[SAFETY_PIPING_PRIVAL_0],"safety_piping_prival_0");
  strcpy(name[SAFETY_PIPING_PRIVAL_1],"safety_piping_prival_1");
  strcpy(name[SAFETY_PIPING_PRIVAL_2],"safety_piping_prival_2");
  strcpy(name[SAFETY_LIFTING_PRIVAL_0],"safety_lifting_prival_0");
  strcpy(name[SAFETY_LIFTING_PRIVAL_1],"safety_lifting_prival_1");
  strcpy(name[SAFETY_LIFTING_PRIVAL_2],"safety_lifting_prival_2");
  strcpy(name[SAFETY_PIPING_GLOBAL_X],"safety_piping_global_x");
  strcpy(name[SAFETY_PIPING_GLOBAL_Y],"safety_piping_global_y");
  strcpy(name[SAFETY_PIPING_GLOBAL_Z],"safety_piping_global_z");
  strcpy(name[SAFETY_LIFTING_GLOBAL_X],"safety_lifting_global_x");
  strcpy(name[SAFETY_LIFTING_GLOBAL_Y],"safety_lifting_global_y");
  strcpy(name[SAFETY_LIFTING_GLOBAL_Z],"safety_lifting_global_z");

  // group_groundflow_expansion (manual Professional 6.614): thermal
  // expansion coefficient of the pore fluid/soil used by the
  // heat-exchanger groundflow analyses (corpus heat_exchanger_pile_dt
  // and _load_dt). Registered parse-only in this batch (consumption
  // belongs to the coupled thermo-groundflow sprint).
  strcpy(name[GROUP_GROUNDFLOW_EXPANSION],"group_groundflow_expansion");
  type[GROUP_GROUNDFLOW_EXPANSION] = DOUBLE_PRECISION;
  data_length[GROUP_GROUNDFLOW_EXPANSION] = 1;
  data_class[GROUP_GROUNDFLOW_EXPANSION] = GROUNDFLOW;
  data_required[GROUP_GROUNDFLOW_EXPANSION] = GROUP_TYPE;

  // control_print_gid_contact_spring2 (manual Professional 6.297):
  // number of nodes (1 or 2) used to draw the contact_spring2 elements
  // in the GiD output. Pure rendering switch - registered parse-only
  // (the GNU GiD printer draws the contact springs with the element
  // nodes of the model).
  // control_print_gid_contact_spring2 (manual Professional 6.297):
  // number of nodes (1 or 2) used to draw the contact_spring2 elements
  // in the GiD output. Pure rendering switch - registered parse-only
  // (the GNU GiD printer draws the contact springs with the element
  // nodes of the model). no_index: the corpus spelling
  // "print_gid_contact_spring2 1" (tutorial_3) carries the node count
  // as the only value; the manual "index" behaves as a run-wide switch.
  strcpy(name[CONTROL_PRINT_GID_CONTACT_SPRING2],"control_print_gid_contact_spring2");
  type[CONTROL_PRINT_GID_CONTACT_SPRING2] = INTEGER;
  data_length[CONTROL_PRINT_GID_CONTACT_SPRING2] = 1;
  no_index[CONTROL_PRINT_GID_CONTACT_SPRING2] = 1;
  data_class[CONTROL_PRINT_GID_CONTACT_SPRING2] = CONTROL;

  // control_mesh_macro_concentrate (manual Professional 6.207): mesh
  // fineness concentration factors of the -rectangle macro (two factors
  // per direction: at the beginning and at the end). Registered
  // parse-only in this batch (the graded-macro generation belongs to
  // the mesh macro sprint).
  strcpy(name[CONTROL_MESH_MACRO_CONCENTRATE],"control_mesh_macro_concentrate");
  type[CONTROL_MESH_MACRO_CONCENTRATE] = DOUBLE_PRECISION;
  data_length[CONTROL_MESH_MACRO_CONCENTRATE] = 4;
  data_class[CONTROL_MESH_MACRO_CONCENTRATE] = CONTROL;

  // mesh_interface_triangle family (manual Professional 6.856/6.857/
  // 6.201, interface11 of the corpus): generate zero-thickness interface
  // elements by cutting a 3D tet4 mesh with a triangulated plane.
  //   mesh_interface_triangle_coordinate  index x0 y0 z0 x1 y1 z1 x2 y2 z2
  //     [x3 y3 z3 ...]   - the plane triangle(s), 9 coordinates each.
  //     The SINGULAR spelling "coordinate" is the one the Professional
  //     binary 25-10-2023 accepts and echoes (the manual text/TOC say
  //     "coordinates"; the corpus and the Pro .dbs use "coordinate").
  //   mesh_interface_triangle_element_group index element_group - the
  //     group attributed to the generated interface elements.
  //   control_mesh_interface_triangle index -yes/-no - switch that
  //     activates the generation at the given control index.
  // The interface elements are generated by generate_interface_triangle()
  // in generate.cc (dispatch in step_start, top.cc).
  strcpy(name[MESH_INTERFACE_TRIANGLE_COORDINATE],
    "mesh_interface_triangle_coordinate");
  type[MESH_INTERFACE_TRIANGLE_COORDINATE] = DOUBLE_PRECISION;
  data_length[MESH_INTERFACE_TRIANGLE_COORDINATE] = DATA_ITEM_SIZE;
  fixed_length[MESH_INTERFACE_TRIANGLE_COORDINATE] = 0;
  data_class[MESH_INTERFACE_TRIANGLE_COORDINATE] = CONTROL;

  strcpy(name[MESH_INTERFACE_TRIANGLE_ELEMENT_GROUP],
    "mesh_interface_triangle_element_group");
  type[MESH_INTERFACE_TRIANGLE_ELEMENT_GROUP] = INTEGER;
  data_length[MESH_INTERFACE_TRIANGLE_ELEMENT_GROUP] = 1;
  data_class[MESH_INTERFACE_TRIANGLE_ELEMENT_GROUP] = CONTROL;

  strcpy(name[CONTROL_MESH_INTERFACE_TRIANGLE],
    "control_mesh_interface_triangle");
  type[CONTROL_MESH_INTERFACE_TRIANGLE] = INTEGER;
  data_length[CONTROL_MESH_INTERFACE_TRIANGLE] = 1;
  data_class[CONTROL_MESH_INTERFACE_TRIANGLE] = CONTROL;

  // control_print_gid_method (manual Professional 6.29x): GiD output
  // method selector (-element / -all_nodes / ...). Registered parse-only:
  // the corpus test interface11 carries the record and the GNU GiD
  // printer is self-contained (it does not dispatch on the method).
  strcpy(name[CONTROL_PRINT_GID_METHOD],"control_print_gid_method");
  type[CONTROL_PRINT_GID_METHOD] = INTEGER;
  data_length[CONTROL_PRINT_GID_METHOD] = 1;
  data_class[CONTROL_PRINT_GID_METHOD] = CONTROL;

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
    db_read[data_number] = 1;
    length = db_len( data_number, index, version );
  }
  else if ( task==GET_IF_EXISTS ) {
    if ( !db_active_index(data_number,index,version) ) {
      length = 0;
      return 0;
    }
    db_read[data_number] = 1;
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
  db_read[data_number] = 1;

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
  db_read[data_number] = 1;

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

#ifdef USE_EXCEPTIONS
  throw DatabaseException("Database error for item " + std::to_string(idat) + " and index " + std::to_string(index));
#else
  exit(TN_EXIT_STATUS);
#endif
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
  db_read[data_number] = 1;

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
  db_read[data_number] = 1;

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

long int db_number( const char str[] )

{
  long int data_number=0, found=-1;

  // "mesh" is the Professional short alias of options_mesh (the GNU
  // MESH record is a dead placeholder with no type/class; routing it
  // through the exact-name loop would return MESH with no_index=0 and
  // the parser would read "-fixed_in_space" as an illegal index).
  if ( !strcmp( str, "mesh" ) )
    return OPTIONS_MESH;
  // "nonlocal" (manual Professional 6.897) is the short name of
  // options_nonlocal: the radius of the nonlocal averaging (the GNU
  // machinery reads OPTIONS_NONLOCAL in nonloc.cc).
  if ( !strcmp( str, "nonlocal" ) )
    return OPTIONS_NONLOCAL;
  // the Professional writes the visco exponential table with the
  // singular "value"; the GNU canonical name is ..._values.
  if ( !strcmp( str, "group_materi_plasti_visco_exponential_value" ) )
    return GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL_VALUES;
  // control_mesh_refine_locally_dof (manual Professional 6.222) is the
  // Professional name of the legacy control_mesh_refine_locally_unknown.
  if ( !strcmp( str, "control_mesh_refine_locally_dof" ) )
    return CONTROL_MESH_REFINE_LOCALLY_UNKNOWN;
  // control_mesh_merge_not: legacy GNU spelling of the Professional
  // control_mesh_merge_geometry_not (manual 6.215)
  if ( !strcmp( str, "control_mesh_merge_not" ) )
    return CONTROL_MESH_MERGE_NOT;
  // control_mesh_cut_force (manual Professional 6.164): the manual
  // name of the record the corpus/Professional write as
  // control_mesh_cut_node_force (6.163 companion; the .dbs of the
  // Professional 25-10-2023 stores control_mesh_cut_node_force).
  if ( !strcmp( str, "control_mesh_cut_force" ) )
    return CONTROL_MESH_CUT_NODE_FORCE;
  // control_mesh_merge_eps_coord (manual 6.212): Professional spelling
  // of the legacy control_mesh_merge_epscoord
  if ( !strcmp( str, "control_mesh_merge_eps_coord" ) )
    return CONTROL_MESH_MERGE_EPSCOORD;

  for ( data_number=0; data_number<MDAT && found<0; data_number++ ) {
    if ( !strcmp(str,db_name(data_number)) ) found = data_number;
  }
  if ( found>=0 ) return found;

  // Professional manual names of GNU keywords with another prefix:
  // force_edge_* -> force_element_edge_*, force_volume_* ->
  // force_element_volume_*. Resolve here so that BOTH the keyword
  // detection and the end-of-variable-values detection (input.cc
  // checks db_number(str)>=0 to stop reading data values) see the
  // translated item. Equivalence table in SEGUIMIENTO-CONVERGENCIA.md.
  {
    static char translated[MCHAR];
    if ( !strncmp( str, "force_edge", 10 ) ) {
      strcpy( translated, "force_element_" );
      strncat( translated, &str[6], MCHAR-20 );
      return db_number( translated );
    }
    else if ( !strncmp( str, "force_volume", 12 ) ) {
      strcpy( translated, "force_element_volume" );
      strncat( translated, &str[12], MCHAR-40 );
      return db_number( translated );
    }
    else if ( !strcmp( str, "control_solver" ) )
      return CONTROL_OPTIONS_SOLVER;
    else if ( !strcmp( str, "group_truss_elasti_young" ) )
      return GROUP_TRUSS_YOUNG;
    else if ( !strcmp( str, "control_mesh_delete_geometry_move_node" ) )
      // Professional spelling (manual 6.19x) of the GNU canonical
      // control_mesh_delete_geometry_movenodes record
      return CONTROL_MESH_DELETE_GEOMETRY_MOVENODES;
    else if ( !strcmp( str, "repeat_save_calculate_result" ) )
      // Professional name of the GNU canonical repeat_calculate_result
      // record (the average/variance analysis of control_repeat_save)
      return REPEAT_CALCULATE_RESULT;
    else if ( !strcmp( str, "group_materi_plasti_hypo_wolffersdorff" ) )
      // Professional spelling (double f); GNU canonical is wolfersdorff
      return GROUP_MATERI_PLASTI_HYPO_WOLFERSDORFF;
    else if ( !strcmp( str, "group_materi_plasti_hypo_strain_intergranular" ) )
      // Professional name of the GNU canonical intergranularstrain record
      return GROUP_MATERI_PLASTI_HYPO_INTERGRANULARSTRAIN;
    else if ( !strcmp( str, "group_materi_plasti_hypo_void_ratio_linear" ) )
      // Professional name of the GNU canonical pressuredependentvoidratio
      return GROUP_MATERI_PLASTI_HYPO_PRESSUREDEPENDENTVOIDRATIO;
    else if ( !strcmp( str, "group_materi_plasti_hypo_masin_clay_advanced_direction" ) )
      // Professional spelling; GNU canonical has the avanced typo
      return GROUP_MATERI_PLASTI_HYPO_MASIN_CLAY_AVANCED_DIRECTION;
    else if ( !strcmp( str, "truss_beam" ) )
      return TRUSSBEAM;
    else if ( !strcmp( str, "mesh" ) )
      // Professional manual 6.799? "mesh" is the short alias of
      // options_mesh: mesh -fixed_in_space -fixed_in_space sets the
      // mesh motion (manual Professional 6.10xx). The GNU parser had a
      // dead MESH record (no type/class); route it to OPTIONS_MESH.
      return OPTIONS_MESH;
    else if ( !strcmp( str, "groundflow_phreatic_level" ) )
      return GROUNDFLOW_PHREATICLEVEL;
    else if ( !strcmp( str, "group_porosity" ) )
      // Professional short name (6.637 area) of the GNU canonical
      // group_groundflow_porosity record
      return GROUP_GROUNDFLOW_POROSITY;
    else if ( !strcmp( str, "size_dev" ) )
      return SIZEDEV;
    else if ( !strcmp( str, "contact_spring" ) )
      return CONTACTSPRING;
    else if ( !strcmp( str, "contact_spring2" ) )
      return CONTACTSPRING;
    else if ( !strcmp( str, "control_mesh_generate_contact_spring" ) )
      return CONTROL_MESH_GENERATE_CONTACTSPRING;
    else if ( !strcmp( str, "control_mesh_generate_contact_spring_element" ) )
      return CONTROL_MESH_GENERATE_CONTACTSPRING_ELEMENT;
    else if ( !strcmp( str, "control_mesh_generate_contact_spring_element_group" ) )
      return CONTROL_MESH_GENERATE_CONTACTSPRING_ELEMENT_GROUP;
    else if ( !strcmp( str, "group_contact_spring_direction_automatic" ) )
      return GROUP_CONTACTSPRING_DIRECTION_AUTOMATIC;
    else if ( !strcmp( str, "element_contact_spring_force" ) )
      return ELEMENT_CONTACTSPRING_FORCE;
    else if ( !strcmp( str, "element_contact_spring_strain" ) )
      return ELEMENT_CONTACTSPRING_FORCE;
    else if ( !strcmp( str, "group_contact_spring_stiffness" ) )
      return db_number( "group_contactspring_stiffness" );
    else if ( !strcmp( str, "group_contact_spring_direction" ) )
      return db_number( "group_contactspring_direction" );
    else if ( !strcmp( str, "group_contact_spring_plasti_friction" ) )
      return db_number( "group_contactspring_friction" );
    else if ( !strcmp( str, "group_contact_spring_friction_automatic" ) )
      return db_number( "group_contactspring_friction_automatic" );
    else if ( !strcmp( str, "group_contact_spring_cohesion" ) )
      return db_number( "group_contactspring_cohesion" );
    else if ( !strcmp( str, "group_contact_spring_memory" ) )
      return db_number( "group_contactspring_memory" );
    else if ( !strcmp( str, "total_pressure" ) )
      return TOTAL;
    else if ( !strcmp( str, "static_pressure" ) )
      return STATIC;
    else if ( !strcmp( str, "dynamic_pressure" ) )
      return DYNAMIC;
    else if ( !strcmp( str, "topres" ) )
      return GROUNDFLOW_PRESSURE;
    else if ( !strcmp( str, "tpres" ) )
      // corpus spelling (undrained1/2 of the sfnet suite) of the same
      // total-pressure dof of bounda_dof (manual Professional 2.4.1)
      return GROUNDFLOW_PRESSURE;
    else if ( !strcmp( str, "inertia_apply" ) )
      return OPTIONS_INERTIA;
    else if ( !strcmp( str, "control_inertia_apply" ) )
      return CONTROL_OPTIONS_INERTIA;
    // convection_apply (manual Professional 6.395) is the input name of
    // options_convection; control_convection_apply (6.113) the indexed
    // per-timestep form; convection_stabilization (6.396) maps to
    // options_stabilization (-no/-yes/-maximal, default -yes).
    else if ( !strcmp( str, "convection_apply" ) )
      return OPTIONS_CONVECTION;
    else if ( !strcmp( str, "control_convection_apply" ) )
      return CONTROL_OPTIONS_CONVECTION;
    else if ( !strcmp( str, "convection_stabilization" ) )
      return OPTIONS_STABILIZATION;
    else if ( !strcmp( str, "materi_elasti_young_power_apply" ) )
      return CONTROL_MATERI_ELASTI_YOUNG_POWER_APPLY;
    else if ( !strcmp( str, "control_solver_bicg_error" ) )
      return CONTROL_OPTIONS_SOLVER_BICG_ERROR;
    else if ( !strcmp( str, "axisymmetric" ) )
      return GROUP_AXISYMMETRIC;
    else if ( !strcmp( str, "group_materi_plasti_druck_prag" ) )
      return GROUP_MATERI_PLASTI_DRUCKPRAG;
    else if ( !strcmp( str, "group_materi_plasti_cap2" ) )
      return GROUP_MATERI_PLASTI_CAP;
    else if ( !strcmp( str, "group_materi_failure_cruching" ) )
      return GROUP_MATERI_FAILURE_CRUCHING;
    else if ( !strcmp( str, "group_materi_failure_voidfraction" ) )
      return GROUP_MATERI_FAILURE_VOIDFRACTION;
    else if ( !strcmp( str, "group_materi_plasti_bounda" ) )
      return GROUP_MATERI_PLASTI_BOUNDARY;
    else if ( !strcmp( str, "group_materi_plasti_bounda_factor" ) )
      return GROUP_MATERI_PLASTI_BOUNDARY_FACTOR;
    else if ( !strncmp( str, "bounda_print_mesh_dof", 21 ) ) {
      strcpy( translated, "print_mesh_dof" );
      strncat( translated, &str[21], MCHAR-40 );
      return db_number( translated );
    }
    else if ( !strcmp( str, "control_print_dof_rhside" ) )
      return CONTROL_PRINT_UNKNOWNSRHSIDE;
    else if ( !strcmp( str, "print_gid_contact_spring2" ) )
      // Professional short alias (tutorial_3 of the corpus) of
      // control_print_gid_contact_spring2 (manual Professional 6.297)
      return CONTROL_PRINT_GID_CONTACT_SPRING2;
    else if ( !strcmp( str, "control_mesh_generate_truss_beam" ) )
      // Professional spelling (with underscores) of the GNU canonical
      // control_mesh_generate_trussbeam (tutorial_3 of the corpus)
      return CONTROL_MESH_GENERATE_TRUSSBEAM;
    else if ( !strcmp( str, "size_tot" ) )
      // Professional spelling (manual 6.266) of the GNU switch "sizetot".
      return SIZETOT;
    else if ( !strncmp( str, "control_print_mesh_dof", 22 ) ) {
      strcpy( translated, "print_mesh_dof" );
      strncat( translated, &str[22], MCHAR-40 );
      return db_number( translated );
    }
    else if ( !strcmp( str, "geometry_factor" ) )
      // Professional manual 6.527: spatial weighting factors of
      // boundary/force_edge loads along a geometry entity. Same record
      // as the GNU geometry_bounda_factor (read by geometry()): same
      // index as the geometry entity, 2 values = linear, 3 = quadratic.
      return GEOMETRY_BOUNDA_FACTOR;
    else if ( !strcmp( str, "processors" ) )
      // Professional short name of options_processors: number of
      // solver threads (consumed in area.cc/elem.cc).
      return OPTIONS_PROCESSORS;
    else if ( ndim==2 && ( !strcmp( str, "rotx" ) || !strcmp( str, "roty" ) ) ) {
      // 2D beam_rotation declares a SINGLE rotation unknown named rotz
      // (the in-plane rotation). The Professional beam inputs prescribe
      // all three of -rotx -roty -rotz (its 2D beam model keeps the
      // three rotation dofs); map the two out-of-plane names onto the
      // same unknown so the bounda_dof records parse.
      long int idat_tmp=0;
      for ( idat_tmp=0; idat_tmp<MDAT; idat_tmp++ )
        if ( !strcmp( name[idat_tmp], "rotz" ) ) return idat_tmp;
    }
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

void check_used_report( void )

{
  long int idat=0, nused=0;

  // report data items that were defined in the data file
  // but never read during the calculation
  for ( idat=0; idat<MDAT; idat++ ) {
    if ( external[idat] && max_index[idat][VERSION_NORMAL]>=0 &&
         !db_read[idat] ) {
      if ( nused==0 )
        cout << "\nCheck_used: data items defined in data file but not used:\n";
      cout << "  " << db_name(idat) << "\n";
      nused++;
    }
  }
  if ( nused==0 )
    cout << "\nCheck_used: all data items from data file were used.\n";
  else
    cout << "Check_used: " << nused << " data items not used.\n";
}


void check_data_integrity( void )

{
  // check_data -yes: verify that data items required by active items exist
  // (data_required refers to another data item that must be active).
  long int idat=0, ireq=0, max=0, index=0;
  long int nmissing=0;

  for ( idat=0; idat<MDAT; idat++ ) {
    if ( data_required[idat]>0 && external[idat] ) {
      ireq = data_required[idat];
      db_max_index( idat, max, VERSION_NORMAL, GET );
      for ( index=0; index<=max; index++ ) {
        if ( db_active_index( idat, index, VERSION_NORMAL ) &&
             !db_active_index( ireq, index, VERSION_NORMAL ) ) {
          cout << "Warning: data item " << db_name(idat) << " (index "
               << index << ") requires " << db_name(ireq)
               << " which is not specified.\n";
          nmissing++;
        }
      }
    }
  }
  if ( nmissing>0 )
    cout << "check_data: " << nmissing
         << " data item(s) have a missing required item.\n";
  else
    cout << "check_data: no missing required data items.\n";
}
