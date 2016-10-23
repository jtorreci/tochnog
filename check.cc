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

long int check( long int idat, long int task )

{
  long int data_number=0, ok=1;

  data_number = labs( idat );

  if ( data_number==HEX8 )
    return 0;
  if ( data_number==HEX27 )
    return 0;
  if ( data_number==TET4 )
    return 0;
  if ( data_number==TET10 )
    return 0;
  if ( data_number==QUAD4 )
    return 0;
  if ( data_number==QUAD9 )
    return 0;
  if ( data_number==QUAD16 )
    return 0;
  if ( data_number==TRIA3 )
    return 0;
  if ( data_number==TRIA6 )
    return 0;
  if ( data_number==BAR2 )
    return 0;
  if ( data_number==BAR3 )
    return 0;

  if ( data_number==NODE ) 
    ok = 1;
  if ( data_number==ELEMENT ) 
    ok = 1;
  if ( data_number==TRUSS ) {
    ok = check_unknown_atleastone( "materi_velocity_integrated", 
      "materi_displacement", task );
  }
  if ( data_number==TRUSSBEAM ) {
    ok = check_unknown_atleastone( "materi_velocity_integrated", 
      "materi_displacement", task );
    ok = ok && check_unknown( "beam_rotation", YES, task );
    ok = ok && check_ndim( 2, 3, task );
  }
  if ( data_number==BAR )
    ok = check_ndim( 1, 1, task );
  if ( data_number==BOUNDA_FORCE )
    ok = check_unknowns_are_specified( task );
  if ( data_number==BOUNDA_SINE )
    ok = check_unknowns_are_specified( task );
  if ( data_number==BOUNDA_TIME )
    ok = check_unknowns_are_specified( task );
  if ( data_number==BOUNDA_TIME_FILE )
    ok = check_unknowns_are_specified( task );
  if ( data_number==BOUNDA_UNKNOWN )
    ok = check_unknowns_are_specified( task );
  if ( data_number==BRICK )
    ok = check_ndim( 3, 3, task );
  if ( data_number==CIRCLE )
    ok = check_ndim( 2, 2, task );
  if ( data_number==CIRCLE_HOLLOW )
    ok = check_ndim( 2, 2, task );
  if ( data_number==CONDIF ) {
    ok = check_unknown( "condif_temperature", YES, task );
    ok = ok && check_unknown( "wave_scalar", NO, task );
  }
  if ( data_number==CONDIF_CONVECTION ) {
    ok = check_ndim( 2, 3, task );
    ok = ok && check_unknown( "condif_temperature", YES, task );
  }
  if ( data_number==CONDIF_CONVECTION_GEOMETRY ) {
    ok = check_ndim( 2, 3, task );
    ok = ok && check_unknown( "condif_temperature", YES, task );
  }
  if ( data_number==CONDIF_RADIATION ) {
    ok = check_ndim( 2, 3, task );
    ok = check_unknown( "condif_temperature", YES, task );
  }
  if ( data_number==CONDIF_RADIATION_GEOMETRY ) {
    ok = check_ndim( 2, 3, task );
    ok = check_unknown( "condif_temperature", YES, task );
  }
  if ( data_number==CONTACT_FRICTION )
    ok = check_unknown( "materi_velocity", YES, task );
  if ( data_number==CONTACT_GEOMETRY )
    ok = check_unknown( "materi_velocity", YES, task );
  if ( data_number==CONTACT_GEOMETRY_SWITCH )
    ok = check_unknown( "materi_velocity", YES, task );
  if ( data_number==CONTACT_HEATGENERATION ) {
    ok = check_unknown( "materi_velocity", YES, task );
    ok = ok && check_unknown( "condif_temperature", YES, task );
  }
  if ( data_number==CONTACT_PENALTY_PRESSURE ) {
    ok = check_unknown( "materi_velocity", YES, task );
    ok = ok && check_unknown( "groundflow_pressure", YES, task );
  }
  if ( data_number==CONTACT_PENALTY_TEMPERATURE ) {
    ok = check_unknown( "materi_velocity", YES, task );
    ok = ok && check_unknown( "condif_temperature", YES, task );
  }
  if ( data_number==CONTACT_PENALTY_VELOCITY )
    ok = check_unknown( "materi_velocity", YES, task );
  if ( data_number==CONTACT_RELAXATION )
    ok = check_unknown( "materi_velocity", YES, task );
  if ( data_number==CONTACTSPRING ) {
    if      ( !materi_displacement )
      ok = check_unknown( "materi_velocity_integrated", YES, task );
    else if ( !materi_velocity_integrated )
      ok = check_unknown( "materi_displacement", YES, task );
  }
  if ( data_number==CONTROL_CRACK ) {
    ok = check_unknown( "materi_velocity", YES, task );
    ok = ok && check_unknown( "materi_stress", YES, task );
  }
  if ( data_number==CONTROL_EIGEN )
    ok = check_unknowns_are_specified( task );
  if ( data_number==CONTROL_EIGEN_SCALE )
    ok = check_unknowns_are_specified( task );
  if ( data_number==CONTROL_MESH_GENERATE_BEAM ) {
    ok = check_unknown_atleastone( "materi_velocity_integrated", 
      "materi_displacement", task );
    ok = ok && check_unknown( "beam_rotation", YES, task );
    ok = ok && check_ndim( 2, 3, task );
  }
  if ( data_number==CONTROL_MESH_GENERATE_CONTACTSPRING )
    ok = check_unknown( "materi_velocity", YES, task );
  if ( data_number==CONTROL_MESH_GENERATE_CONTACTSPRING_ELEMENT )
    ok = check_unknown( "materi_velocity", YES, task );
  if ( data_number==CONTROL_MESH_GENERATE_SPRING1 )
    ok = check_unknown( "materi_velocity", YES, task );
  if ( data_number==CONTROL_MESH_GENERATE_SPRING2 )
    ok = check_unknown( "materi_velocity", YES, task );
  if ( data_number==CONTROL_MESH_GENERATE_TRUSS )
    ok = check_unknown( "materi_velocity", YES, task );
  if ( data_number==CONTROL_MESH_GENERATE_TRUSSBEAM ) {
    ok = check_unknown_atleastone( "materi_velocity_integrated", 
      "materi_displacement", task );
    ok = ok && check_unknown( "beam_rotation", YES, task );
    ok = ok && check_ndim( 2, 3, task );
  }
  if ( data_number==CONTROL_MESH_GENERATE_TRUSS_BEAM_LOOSE )
    ok = check_unknown( "materi_velocity", YES, task );
  if ( data_number==CONTROL_MESH_GENERATE_TRUSS_BEAM_MACRO )
    ok = check_unknown( "materi_velocity", YES, task );
  if ( data_number==CONTROL_MESH_MACRO_ELEMENT )
    ok = check_ndim( 2, 3, task );
  if ( data_number==CONTROL_MESH_SPLIT )
    ok = check_ndim( 2, 3, task );
  if ( data_number==CONTROL_OPTIONS_SOLVER )
    ok = check_unknowns_are_specified( task );
  if ( data_number==CONTROL_PRINT_GMV )
    ok = check_ndim( 2, 3, task );
  if ( data_number==CONTROL_PRINT_GID )
    ok = check_ndim( 2, 3, task );
  if ( data_number==CONTROL_PRINT_UNKNOWNS )
    ok = check_unknowns_are_specified( task );
  if ( data_number==CONTROL_PRINT_UNKNOWNSRHSIDE )
    ok = check_unknowns_are_specified( task );
  if ( data_number==CONTROL_RELAXATION_CONDIF_TEMPERATURE )
    ok = check_unknown( "condif_temperature", YES, task );
  if ( data_number==CONTROL_RELAXATION_MAXWELL_E ) {
    ok = check_ndim( 2, 3, task );
    ok = ok && check_unknown( "maxwell_er", YES, task );
    ok = ok && check_unknown( "maxwell_ei", YES, task );
  }
  if ( data_number==CONTROL_RELAXATION_GROUNDFLOW_PRESSURE )
    ok = check_unknown( "groundflow_pressure", YES, task );
  if ( data_number==CONTROL_RELAXATION_WAVE_FSCALAR )
    ok = check_unknown( "wave_fscalar", YES, task );
  if ( data_number==CONTROL_RELAXATION_MATERI_VELOCITY )
    ok = check_unknown( "materi_velocity", YES, task );
  if ( data_number==CONTROL_UNKNOWN_FREEZE )
    ok = check_unknowns_are_specified( task );
  if ( data_number==CONTROL_UNKNOWN_RESET_GEOMETRY )
    ok = check_unknowns_are_specified( task );
  if ( data_number==CONTROL_UNKNOWN_RESET_UNKNOWN )
    ok = check_unknowns_are_specified( task );
  if ( data_number==CONTROL_UNKNOWN_RESET_VALUE )
    ok = check_unknowns_are_specified( task );
  if ( data_number==CONTROL_TIMESTEP_ITERATIONS_AUTOMATIC )
    ok = check_unknowns_are_specified( task );
  if ( data_number==CRACK_DIRECTION ) {
    ok = check_unknown( "materi_velocity", YES, task );
    ok = ok && check_unknown( "materi_stress", YES, task );
  }
  if ( data_number==CRACK_ELEMENTGROUP ) {
    ok = check_unknown( "materi_velocity", YES, task );
    ok = ok && check_unknown( "materi_stress", YES, task );
  }
  if ( data_number==CRACK_LENGTH ) {
    ok = check_unknown( "materi_velocity", YES, task );
    ok = ok && check_unknown( "materi_stress", YES, task );
  }
  if ( data_number==CRACK_NODES ) {
    ok = check_unknown( "materi_velocity", YES, task );
    ok = ok && check_unknown( "materi_stress", YES, task );
  }
  if ( data_number==CRACK_STRESSINTENSITYFACTOR ) {
    ok = check_unknown( "materi_velocity", YES, task );
    ok = ok && check_unknown( "materi_stress", YES, task );
  }
  if ( data_number==CRACK_TIP ) {
    ok = check_unknown( "materi_velocity", YES, task );
    ok = ok && check_unknown( "materi_stress", YES, task );
  }
  if ( data_number==CYLINDER_HOLLOW )
    ok = check_ndim( 3, 3, task );
  if ( data_number==DEPENDENCY_DIAGRAM )
    ok = check_unknowns_are_specified( task );
  if ( data_number==DEPENDENCY_ITEM )
    ok = check_unknowns_are_specified( task );
  if ( data_number==DOF_LABEL )
    ok = check_unknowns_are_specified( task );
  if ( data_number==ELEMENT_BEAM_MOMENT ) {
    ok = check_unknown_atleastone( "materi_velocity_integrated", 
      "materi_displacement", task );
    ok = ok && check_unknown( "beam_rotation", YES, task );
    ok = ok && check_ndim( 2, 3, task );
  }
  if ( data_number==ELEMENT_TRUSS_FORCE )
    ok = check_unknown( "materi_velocity", YES, task );
  if ( data_number==FORCE_ELEMENT_EDGE )
    ok = check_unknowns_are_specified( task );
  if ( data_number==FORCE_ELEMENT_EDGE_NORMAL )
    ok = check_unknown( "materi_velocity", YES, task );
  if ( data_number==FORCE_ELEMENT_EDGE_NORMAL_SINE )
    ok = check_unknowns_are_specified( task );
  if ( data_number==FORCE_ELEMENT_EDGE_NORMAL_TIME )
    ok = check_unknowns_are_specified( task );
  if ( data_number==FORCE_ELEMENT_EDGE_GEOMETRY )
    ok = check_unknowns_are_specified( task );
  if ( data_number==FORCE_ELEMENT_EDGE_SINE )
    ok = check_unknowns_are_specified( task );
  if ( data_number==FORCE_ELEMENT_EDGE_TIME )
    ok = check_unknowns_are_specified( task );
  if ( data_number==FORCE_ELEMENT_VOLUME )
    ok = check_unknowns_are_specified( task );
  if ( data_number==FORCE_ELEMENT_VOLUME_FACTOR )
    ok = check_unknowns_are_specified( task );
  if ( data_number==FORCE_ELEMENT_VOLUME_GEOMETRY )
    ok = check_unknowns_are_specified( task );
  if ( data_number==FORCE_ELEMENT_VOLUME_SINE ) {
    ok = check_unknowns_are_specified( task );
    ok = ok && check_ndim( 1, 2, task );
  }
  if ( data_number==FORCE_ELEMENT_VOLUME_TIME ) {
    ok = check_unknowns_are_specified( task );
    ok = ok && check_ndim( 1, 2, task );
  }
  if ( data_number==FORCE_GRAVITY ) {
    ok = check_unknown_atleastone( "materi_stress", 
      "groundflow_pressure", task );
  }
  if ( data_number==FORCE_GRAVITY_TIME ) {
    ok = check_unknown_atleastone( "materi_stress", 
      "groundflow_pressure", task );
  }
  if ( data_number==GEOMETRY_BOUNDA_FACTOR )
    ok = check_ndim( 2, 3, task );
  if ( data_number==GEOMETRY_BOUNDA_SINE_X )
    ok = check_ndim( 1, 3, task );
  if ( data_number==GEOMETRY_BOUNDA_SINE_Y )
    ok = check_ndim( 2, 3, task );
  if ( data_number==GEOMETRY_BOUNDA_SINE_Z )
    ok = check_ndim( 3, 3, task );
  if ( data_number==GEOMETRY_CIRCLE )
    ok = check_ndim( 2, 2, task );
  if ( data_number==GEOMETRY_ELLIPSE )
    ok = check_ndim( 2, 2, task );
  if ( data_number==GEOMETRY_CYLINDER )
    ok = check_ndim( 3, 3, task );
  if ( data_number==GEOMETRY_CYLINDER_SEGMENT )
    ok = check_ndim( 3, 3, task );
  if ( data_number==GEOMETRY_POLYNOMIAL )
    ok = check_ndim( 2, 3, task );
  if ( data_number==GEOMETRY_QUADRILATERAL )
    ok = check_ndim( 2, 3, task );
  if ( data_number==GEOMETRY_SPHERE )
    ok = check_ndim( 3, 3, task );
  if ( data_number==GEOMETRY_SPHERE_SEGMENT )
    ok = check_ndim( 3, 3, task );
  if ( data_number==GEOMETRY_TRIANGLE )
    ok = check_ndim( 2, 3, task );
  if ( data_number==GROUNDFLOW ) {
    ok = check_unknown( "groundflow_pressure", YES, task );
    ok = check_unknown( "wave_scalar", NO, task );
  }
  if ( data_number==GROUNDFLOW_DENSITY )
    ok = check_unknown( "groundflow_pressure", YES, task );
  if ( data_number==GROUNDFLOW_PHREATICLEVEL )
    ok = check_unknown_atleastone( "materi_stress", 
      "groundflow_pressure", task );
  if ( data_number==GROUNDFLOW_PHREATICLEVEL_N ) {
    ok = check_unknown( "groundflow_pressure", YES, task );
    ok = ok && check_ndim( 3, 3, task );
  }
  if ( data_number==GROUNDFLOW_PHREATICLEVEL_BOUNDA )
    ok = check_unknown( "groundflow_pressure", YES, task );
  if ( data_number==GROUNDFLOW_PRESSURE )
    ok = check_unknown( "groundflow_pressure", YES, task );
  if ( data_number==GROUNDFLOW_PRESSURE_ATMOSPHERIC )
    ok = check_unknown( "groundflow_pressure", YES, task );
  if ( data_number==GROUP_BEAM_INERTIA ) {
    ok = check_unknown_atleastone( "materi_velocity_integrated", 
      "materi_displacement", task );
    ok = ok && check_unknown( "beam_rotation", YES, task );
    ok = ok && check_ndim( 2, 3, task );
  }
  if ( data_number==GROUP_BEAM_MEMORY ) {
    ok = check_unknown_atleastone( "materi_velocity_integrated", 
      "materi_displacement", task );
    ok = ok && check_unknown( "beam_rotation", YES, task );
    ok = ok && check_ndim( 2, 3, task );
  }
  if ( data_number==GROUP_BEAM_PLANE ) {
    ok = check_unknown_atleastone( "materi_velocity_integrated", 
      "materi_displacement", task );
    ok = ok && check_unknown( "beam_rotation", YES, task );
    ok = ok && check_ndim( 3, 3, task );
  }
  if ( data_number==GROUP_BEAM_YOUNG ) {
    ok = check_unknown_atleastone( "materi_velocity_integrated", 
      "materi_displacement", task );
    ok = ok && check_unknown( "beam_rotation", YES, task );
    ok = ok && check_ndim( 2, 3, task );
  }
  if ( data_number==GROUP_CONDIF_ABSORPTION )
    ok = check_unknown( "condif_temperature", YES, task );
  if ( data_number==GROUP_CONDIF_CONDUCTIVITY )
    ok = check_unknown( "condif_temperature", YES, task );
  if ( data_number==GROUP_CONDIF_DENSITY )
    ok = check_unknown( "condif_temperature", YES, task );
  if ( data_number==GROUP_CONDIF_CAPACITY )
    ok = check_unknown( "condif_temperature", YES, task );
  if ( data_number==GROUP_CONDIF_FLOW )
    ok = check_unknown( "condif_temperature", YES, task );
  if ( data_number==GROUP_CONTACTSPRING_COHESION ) {
    if      ( !materi_displacement )
      ok = check_unknown( "materi_velocity_integrated", YES, task );
    else if ( !materi_velocity_integrated )
      ok = check_unknown( "materi_displacement", YES, task );
  }
  if ( data_number==GROUP_CONTACTSPRING_DIRECTION ) {
    if      ( !materi_displacement )
      ok = check_unknown( "materi_velocity_integrated", YES, task );
    else if ( !materi_velocity_integrated )
      ok = check_unknown( "materi_displacement", YES, task );
  }
  if ( data_number==GROUP_CONTACTSPRING_FRICTION ) {
    if      ( !materi_displacement )
      ok = check_unknown( "materi_velocity_integrated", YES, task );
    else if ( !materi_velocity_integrated )
      ok = check_unknown( "materi_displacement", YES, task );
  }
  if ( data_number==GROUP_CONTACTSPRING_FRICTION_AUTOMATIC ) {
    if      ( !materi_displacement )
      ok = check_unknown( "materi_velocity_integrated", YES, task );
    else if ( !materi_velocity_integrated )
      ok = check_unknown( "materi_displacement", YES, task );
  }
  if ( data_number==GROUP_CONTACTSPRING_MEMORY ) {
    if      ( !materi_displacement )
      ok = check_unknown( "materi_velocity_integrated", YES, task );
    else if ( !materi_velocity_integrated )
      ok = check_unknown( "materi_displacement", YES, task );
  }
  if ( data_number==GROUP_CONTACTSPRING_STIFFNESS ) {
    if      ( !materi_displacement )
      ok = check_unknown( "materi_velocity_integrated", YES, task );
    else if ( !materi_velocity_integrated )
      ok = check_unknown( "materi_displacement", YES, task );
  }
  if ( data_number==GROUP_GROUNDFLOW_CAPACITY )
    ok = check_unknown( "groundflow_pressure", YES, task );
  if ( data_number==GROUP_GROUNDFLOW_MATERIDIVERGENCE ) {
    ok = check_unknown( "groundflow_pressure", YES, task );
    ok = ok && check_unknown( "materi_velocity", YES, task );
  }
  if ( data_number==GROUP_GROUNDFLOW_PERMEABILITY )
    ok = check_unknown( "groundflow_pressure", YES, task );
  if ( data_number==GROUP_INTEGRATION_POINTS )
    ok = check_unknowns_are_specified( task );
  if ( data_number==GROUP_MATERI_DAMAGE_MAZARS )
    ok = check_unknown( "materi_damage", YES, task );
  if ( data_number==GROUP_MATERI_DAMPING )
    ok = check_unknown( "materi_velocity", YES, task );
  if ( data_number==GROUP_MATERI_DENSITY )
    ok = check_unknown( "materi_velocity", YES, task );
  if ( data_number==GROUP_MATERI_DENSITY_GROUNDFLOW ) {
    ok = check_unknown( "materi_velocity", YES, task );
    ok = ok && check_unknown( "groundflow_pressure", YES, task );
  }
  if ( data_number==GROUP_MATERI_ELASTI_CAMCLAY_G ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_history_variables", YES, task );
  }
  if ( data_number==GROUP_MATERI_ELASTI_CAMCLAY_POISSON ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_history_variables", YES, task );
  }
  if ( data_number==GROUP_MATERI_ELASTI_COMPRESSIBILITY ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_velocity", YES, task );
  }
  if ( data_number==GROUP_MATERI_ELASTI_LADE ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_velocity", YES, task );
    ok = ok && check_unknown( "materi_strain_total", YES, task );
  }
  if ( data_number==GROUP_MATERI_ELASTI_POISSON ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_velocity", YES, task );
  }
  if ( data_number==GROUP_MATERI_ELASTI_TRANSVERSE_ISOTROPY ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_velocity", YES, task );
    ok = ok && check_unknown( "materi_strain_total", YES, task );
  }
  if ( data_number==GROUP_MATERI_ELASTI_YOUNG ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_velocity", YES, task );
  }
  if ( data_number==GROUP_MATERI_ELASTI_YOUNG_POLYNOMIAL ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_velocity", YES, task );
    ok = ok && check_unknown( "materi_strain_total", YES, task );
  }
  if ( data_number==GROUP_MATERI_ELASTI_YOUNG_POWER ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_velocity", YES, task );
    ok = ok && check_unknown( "materi_strain_total", YES, task );
  }
  if ( data_number==GROUP_MATERI_EXPANSION_LINEAR ) {
    ok = check_unknown( "condif_temperature", YES, task );
    ok = ok && check_unknown( "materi_velocity", YES, task );
  }
  if ( data_number==GROUP_MATERI_EXPANSION_VOLUME ) {
    ok = check_unknown( "condif_temperature", YES, task );
    ok = ok && check_unknown( "materi_velocity", YES, task );
  }
  if ( data_number==GROUP_MATERI_FAILURE_DAMAGE ) {
    ok = check_unknown( "materi_damage", YES, task );
    ok = ok && check_unknown( "materi_strain_total", YES, task );
  } 
  if ( data_number==GROUP_MATERI_FAILURE_CRUCHING )
    ok = check_unknown( "materi_strain_total", YES, task );
  if ( data_number==GROUP_MATERI_FAILURE_PLASTI_KAPPA )
    ok = check_unknown( "materi_plasti_kappa", YES, task );
  if ( data_number==GROUP_MATERI_FAILURE_RUPTURE )
    ok = check_unknown( "materi_strain_total", YES, task );
  if ( data_number==GROUP_MATERI_FAILURE_VOIDFRACTION )
    ok = check_unknown( "materi_void_fraction", YES, task );
  if ( data_number==GROUP_MATERI_HYPER_BESSELING ) {
    ok = ok && check_unknown( "materi_strain_total", YES, task );
    ok = ok && check_unknown( "materi_strain_elasti", YES, task );
    ok = ok && check_unknown( "materi_stress", YES, task );
  }
  if ( data_number==GROUP_MATERI_HYPER_MOONEY_RIVLIN ) {
    ok = ok && check_unknown( "materi_strain_total", YES, task );
    ok = ok && check_unknown( "materi_strain_elasti", YES, task );
    ok = ok && check_unknown( "materi_stress", YES, task );
  }
  if ( data_number==GROUP_MATERI_HYPER_VOLUMETRIC_LINEAR ) {
    ok = ok && check_unknown( "materi_strain_total", YES, task );
    ok = ok && check_unknown( "materi_strain_elasti", YES, task );
    ok = ok && check_unknown( "materi_stress", YES, task );
  }
  if ( data_number==GROUP_MATERI_HYPER_VOLUMETRIC_MURNAGHAN ) {
    ok = ok && check_unknown( "materi_strain_total", YES, task );
    ok = ok && check_unknown( "materi_strain_elasti", YES, task );
    ok = ok && check_unknown( "materi_stress", YES, task );
  }
  if ( data_number==GROUP_MATERI_HYPER_VOLUMETRIC_OGDEN ) {
    ok = ok && check_unknown( "materi_strain_total", YES, task );
    ok = ok && check_unknown( "materi_strain_elasti", YES, task );
    ok = ok && check_unknown( "materi_stress", YES, task );
  }
  if ( data_number==GROUP_MATERI_HYPER_VOLUMETRIC_SIMOTAYLOR ) {
    ok = ok && check_unknown( "materi_strain_total", YES, task );
    ok = ok && check_unknown( "materi_strain_elasti", YES, task );
    ok = ok && check_unknown( "materi_stress", YES, task );
  }
  if ( data_number==GROUP_MATERI_MAXWELL_CHAIN )
    ok = check_unknown( "materi_maxwell_stress", YES, task );
  if ( data_number==GROUP_MATERI_MAXWELL_CHAIN_NONLINEAR )
    ok = check_unknown( "materi_maxwell_stress", YES, task );
  if ( data_number==GROUP_MATERI_MEMBRANE ) {
    ok = check_ndim( 1, 2, task );
    ok = ok & check_unknown( "materi_stress", YES, task );
    ok = ok & check_unknown( "groundflow_pressure", NO, task );
  }
  if ( data_number==GROUP_MATERI_MEMORY )
    ok = check_unknown( "materi_stress", YES, task );
  if ( data_number==GROUP_MATERI_PLASTI_CAMCLAY ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_history_variables", YES, task );
    ok = ok && check_unknown( "materi_strain_plasti", YES, task );
  }
  if ( data_number==GROUP_MATERI_PLASTI_CAP ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_strain_plasti", YES, task );
  }
  if ( data_number==GROUP_MATERI_PLASTI_COMPRESSION ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_strain_plasti", YES, task );
  }
  if ( data_number==GROUP_MATERI_PLASTI_DIPRISCO ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_history_variables", YES, task );
    ok = ok && check_unknown( "materi_strain_plasti", YES, task );
  }
  if ( data_number==GROUP_MATERI_PLASTI_DIPRISCO_RT ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_history_variables", YES, task );
    ok = ok && check_unknown( "materi_strain_plasti", YES, task );
  }
  if ( data_number==GROUP_MATERI_PLASTI_DRUCKPRAG ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_strain_plasti", YES, task );
  }
  if ( data_number==GROUP_MATERI_PLASTI_DRUCKPRAG_TENSIONCUTOFF ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_strain_plasti", YES, task );
  }
  if ( data_number==GROUP_MATERI_PLASTI_GURSON ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_strain_plasti", YES, task );
    ok = ok && check_unknown( "materi_void_fraction", YES, task );
  }
  if ( data_number==GROUP_MATERI_PLASTI_HLC ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_strain_plasti", YES, task );
    ok = ok && check_unknown( "materi_void_fraction", YES, task );
  }
  if ( data_number==GROUP_MATERI_PLASTI_HEATGENERATION ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_strain_plasti", YES, task );
    ok = ok && check_unknown( "condif_temperature", YES, task );
  }
  if ( data_number==GROUP_MATERI_PLASTI_HYPO_INTERGRANULARSTRAIN ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_strain_intergranular", YES, task );
  }
  if ( data_number==GROUP_MATERI_PLASTI_HYPO_PRESSUREDEPENDENTVOIDRATIO ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_history_variables", YES, task );
  }
  if ( data_number==GROUP_MATERI_PLASTI_HYPO_WOLFERSDORFF ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_history_variables", YES, task );
  }
  if ( data_number==GROUP_MATERI_PLASTI_KINEMATIC_HARDENING ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_strain_plasti", YES, task );
    ok = ok && check_unknown( "materi_plasti_rho", YES, task );
  }
  if ( data_number==GROUP_MATERI_PLASTI_MATSUOKANAKAI ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_strain_plasti", YES, task );
  }
  if ( data_number==GROUP_MATERI_PLASTI_MOHRCOUL ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_strain_plasti", YES, task );
  }
  if ( data_number==GROUP_MATERI_PLASTI_MOHRCOUL_SOFTENING ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_strain_plasti", YES, task );
    ok = ok && check_unknown( "materi_plasti_kappa", YES, task );
  }
  if ( data_number==GROUP_MATERI_PLASTI_STRESS ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_strain_plasti", YES, task );
  }
  if ( data_number==GROUP_MATERI_PLASTI_TENSION ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_strain_plasti", YES, task );
  }
  if ( data_number==GROUP_MATERI_PLASTI_USER ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_strain_plasti", YES, task );
  }
  if ( data_number==GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_strain_plasti", YES, task );
  }
  if ( data_number==GROUP_MATERI_PLASTI_VISCO_POWER ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_strain_plasti", YES, task );
  }
  if ( data_number==GROUP_MATERI_PLASTI_VONMISES ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_strain_plasti", YES, task );
  }
  if ( data_number==GROUP_MATERI_PLASTI_VONMISES_NADAI ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_strain_plasti", YES, task );
    ok = ok && check_unknown( "materi_plasti_kappa", YES, task );
  }
  if ( data_number==GROUP_MAXWELL_FREQUENCY_EIGEN ) {
    ok = check_unknown( "residue", NO, task );
    ok = ok && check_unknown( "maxwell_er", YES, task );
    ok = ok && check_unknown( "maxwell_ei", YES, task );
  }
  if ( data_number==GROUP_MAXWELL_FREQUENCY_EPSILON ) {
    ok = check_unknown( "residue", NO, task );
    ok = ok && check_unknown( "maxwell_er", YES, task );
    ok = ok && check_unknown( "maxwell_ei", YES, task );
  }
  if ( data_number==GROUP_MAXWELL_FREQUENCY_PML_EPSILONANDMU ) {
    ok = check_unknown( "residue", NO, task );
    ok = ok && check_unknown( "maxwell_er", YES, task );
    ok = ok && check_unknown( "maxwell_ei", YES, task );
  }
  if ( data_number==GROUP_MAXWELL_FREQUENCY_PML_PLANES ) {
    ok = check_unknown( "residue", NO, task );
    ok = ok && check_unknown( "maxwell_er", YES, task );
    ok = ok && check_unknown( "maxwell_ei", YES, task );
  }
  if ( data_number==GROUP_MAXWELL_FREQUENCY_EPSILON_ANISOTROPIC ) {
    ok = check_unknown( "residue", NO, task );
    ok = ok && check_unknown( "maxwell_er", YES, task );
    ok = ok && check_unknown( "maxwell_ei", YES, task );
  }
  if ( data_number==GROUP_MAXWELL_FREQUENCY_J ) {
    ok = check_unknown( "residue", NO, task );
    ok = ok && check_unknown( "maxwell_er", YES, task );
    ok = ok && check_unknown( "maxwell_ei", YES, task );
  }
  if ( data_number==GROUP_MAXWELL_FREQUENCY_MU ) {
    ok = check_unknown( "residue", NO, task );
    ok = ok && check_unknown( "maxwell_er", YES, task );
    ok = ok && check_unknown( "maxwell_ei", YES, task );
  }
  if ( data_number==GROUP_MAXWELL_FREQUENCY_PENALTY ) {
    ok = check_unknown( "residue", NO, task );
    ok = ok && check_unknown( "maxwell_er", YES, task );
    ok = ok && check_unknown( "maxwell_ei", YES, task );
  }
  if ( data_number==GROUP_MAXWELL_TIME_EPSILON ) {
    ok = check_unknown( "residue", NO, task );
    ok = ok && check_unknown( "maxwell_e", YES, task );
  }
  if ( data_number==GROUP_MAXWELL_TIME_EPSILON_ANISOTROPIC ) {
    ok = check_unknown( "residue", NO, task );
    ok = ok && check_unknown( "maxwell_e", YES, task );
  }
  if ( data_number==GROUP_MAXWELL_TIME_J ) {
    ok = check_unknown( "residue", NO, task );
    ok = ok && check_unknown( "maxwell_e", YES, task );
  }
  if ( data_number==GROUP_MAXWELL_TIME_MU ) {
    ok = check_unknown( "residue", NO, task );
    ok = ok && check_unknown( "maxwell_e", YES, task );
  }
  if ( data_number==GROUP_MAXWELL_TIME_PENALTY ) {
    ok = check_unknown( "residue", NO, task );
    ok = ok && check_unknown( "maxwell_e", YES, task );
  }
  if ( data_number==GROUP_SPRING_DIRECTION ) {
    if      ( !materi_displacement )
      ok = check_unknown( "materi_velocity_integrated", YES, task );
    else if ( !materi_velocity_integrated )
      ok = check_unknown( "materi_displacement", YES, task );
  }
  if ( data_number==GROUP_SPRING_PLASTI ) {
    if      ( !materi_displacement )
      ok = check_unknown( "materi_velocity_integrated", YES, task );
    else if ( !materi_velocity_integrated )
      ok = check_unknown( "materi_displacement", YES, task );
  }           
  if ( data_number==GROUP_SPRING_STIFFNESS ) {
    if      ( !materi_displacement )
      ok = check_unknown( "materi_velocity_integrated", YES, task );
    else if ( !materi_velocity_integrated )
      ok = check_unknown( "materi_displacement", YES, task );
  }
  if ( data_number==NODE_RHSIDE )
    ok = check_unknowns_are_specified( task );
  if ( data_number==GROUP_MATERI_STOKES )
    ok = check_unknown( "materi_velocity", YES, task );
  if ( data_number==GROUP_MATERI_VISCOSITY )
    ok = check_unknown( "materi_velocity", YES, task );
  if ( data_number==GROUP_MATERI_VISCOSITY_USER )
    ok = check_unknown( "materi_velocity", YES, task );
  if ( data_number==GROUP_MATERI_VISCOSITY_HEATGENERATION ) 
  {
    ok = check_unknown( "materi_velocity", YES, task );
    ok = ok && check_unknown( "condif_temperature", YES, task );
  }
  if ( data_number==GROUP_TRUSS_AREA )
    ok = check_unknown_atleastone( "materi_velocity_integrated", 
      "materi_displacement", task );
  if ( data_number==GROUP_TRUSS_DENSITY )
    ok = check_unknown_atleastone( "materi_velocity_integrated", 
      "materi_displacement", task );
  if ( data_number==GROUP_TRUSS_MEMORY )
    ok = check_unknown_atleastone( "materi_velocity_integrated", 
      "materi_displacement", task );
  if ( data_number==GROUP_TRUSS_PLASTI )
    ok = check_unknown_atleastone( "materi_velocity_integrated", 
      "materi_displacement", task );
  if ( data_number==GROUP_TRUSS_ROPE )
    ok = check_unknown_atleastone( "materi_velocity_integrated", 
      "materi_displacement", task );
  if ( data_number==GROUP_TRUSS_YOUNG )
    ok = check_unknown_atleastone( "materi_velocity_integrated", 
      "materi_displacement", task );
  if ( data_number==GROUP_TYPE )
    ok = check_unknowns_are_specified( task );
  if ( data_number==GROUP_USER_DATA )
    ok = check_unknown( "materi_strain_total", YES, task );
  if ( data_number==GROUP_USER_UMAT )
    ok = check_unknown( "materi_strain_total", YES, task );
  if ( data_number==GROUP_WAVE_SPEED_OF_SOUND ) {
    ok = check_unknown( "wave_scalar", YES, task );
    ok = ok && check_unknown( "wave_fscalar", YES, task );
  }
  if ( data_number==MATERI )
    ok = ok && check_unknown( "wave_scalar", NO, task );
  if ( data_number==MATERI_STRAINENERGY )
    ok = check_unknown( "materi_strainenergy", YES, task );
  if ( data_number==MATERI_STRAIN_ELASTI )
    ok = check_unknown( "materi_strain_elasti", YES, task );
  if ( data_number==MATERI_STRAIN_PLASTI )
    ok = check_unknown( "materi_strain_plasti", YES, task );
  if ( data_number==MATERI_STRAIN_TOTAL )
    ok = check_unknown( "materi_strain_total", YES, task );
  if ( data_number==MATERI_STRESS )
    ok = check_unknown( "materi_stress", YES, task );
  if ( data_number==MATERI_VELOCITY )
    ok = check_unknown( "materi_velocity", YES, task );
  if ( data_number==MATERI_WORK ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = check_unknown( "materi_velocity", YES, task );
  }
  if ( data_number==MAXWELL_ECOMPLEX ) {
    ok = check_unknown( "maxwell_er", YES, task );
    ok = ok && check_unknown( "maxwell_ei", YES, task );
  }
  if ( data_number==MAXWELL_FREQUENCY_EXCITATION ) {
    ok = check_unknown( "maxwell_er", YES, task );
    ok = ok && check_unknown( "maxwell_ei", YES, task );
    ok = ok && check_unknown( "residue", NO, task );
  }
  if ( data_number==MAXWELL_SCATTER_ENERGYCONSERVATION ) {
    ok = check_unknown( "residue", NO, task );
    ok = ok && check_unknown( "maxwell_er", YES, task );
    ok = ok && check_unknown( "maxwell_ei", YES, task );
  }
  if ( data_number==MAXWELL_SCATTER_MATRIX_AMPLITUDE ) {
    ok = check_ndim( 2, 3, task );
    ok = ok && check_unknown( "maxwell_er", YES, task );
    ok = ok && check_unknown( "maxwell_ei", YES, task );
  }
  if ( data_number==MAXWELL_SCATTER_MATRIX_AMPLITUDEDB ) {
    ok = check_ndim( 2, 3, task );
    ok = ok && check_unknown( "maxwell_er", YES, task );
    ok = ok && check_unknown( "maxwell_ei", YES, task );
  }
  if ( data_number==MAXWELL_SCATTER_MATRIX_IMAGINARY ) {
    ok = check_ndim( 2, 3, task );
    ok = ok && check_unknown( "maxwell_er", YES, task );
    ok = ok && check_unknown( "maxwell_ei", YES, task );
  }
  if ( data_number==MAXWELL_SCATTER_MATRIX_PHASE ) {
    ok = check_ndim( 2, 3, task );
    ok = ok && check_unknown( "maxwell_er", YES, task );
    ok = ok && check_unknown( "maxwell_ei", YES, task );
  }
  if ( data_number==MAXWELL_SCATTER_MATRIX_REAL ) {
    ok = check_ndim( 2, 3, task );
    ok = ok && check_unknown( "maxwell_er", YES, task );
    ok = ok && check_unknown( "maxwell_ei", YES, task );
  }
  if ( data_number==MAXWELL_SCATTER_PARAMETERS ) {
    ok = check_ndim( 2, 3, task );
    ok = ok && check_unknown( "maxwell_er", YES, task );
    ok = ok && check_unknown( "maxwell_ei", YES, task );
  }
  if ( data_number==MAXWELL_SCATTER_PORT_INPUT ) {
    ok = check_ndim( 2, 3, task );
    ok = ok && check_unknown( "maxwell_er", YES, task );
    ok = ok && check_unknown( "maxwell_ei", YES, task );
  }
  if ( data_number==MAXWELL_SCATTER_PORT_OUTPUT ) {
    ok = check_ndim( 2, 3, task );
    ok = ok && check_unknown( "maxwell_er", YES, task );
    ok = ok && check_unknown( "maxwell_ei", YES, task );
  }
  if ( data_number==NODE_DAMPING )
    ok = check_unknown( "materi_velocity", YES, task );
  if ( data_number==NODE_DOF )
    ok = check_unknowns_are_specified( task );
  if ( data_number==NODE_DOF_CALCUL )
    ok = check_unknowns_are_specified( task );
  if ( data_number==NODE_EIGEN )
    ok = check_unknowns_are_specified( task );
  if ( data_number==NODE_MASS )
    ok = check_unknown( "materi_velocity", YES, task );
  if ( data_number==NODE_STIFFNESS ) {
    if      ( !materi_displacement )
      ok = check_unknown( "materi_velocity_integrated", YES, task );
    else if ( !materi_velocity_integrated )
      ok = check_unknown( "materi_displacement", YES, task );
  }
  if ( data_number==GROUP_AXISYMMETRIC ) 
    ok = check_ndim( 1, 2, task );
  if ( data_number==OPTIONS_INERTIA )
    ok = check_unknowns_are_specified( task );
  if ( data_number==OPTIONS_MESH )
    ok = check_unknown( "materi_velocity", YES, task );
  if ( data_number==OPTIONS_NONLOCAL ) {
    ok = check_unknown( "materi_plasti_f", YES, task );
    ok = ok && check_unknown( "materi_plasti_f_nonlocal", YES, task );
  }
  if ( data_number==OPTIONS_NONLOCAL_SOFTVAR ) {
    ok = check_unknown( "materi_plasti_softvar_local", YES, task );
    ok = ok && check_unknown( "materi_plasti_softvar_nonlocal", YES, task );
  }
  if ( data_number==OPTIONS_RESIDUEFACTOR )
    ok = check_unknowns_are_specified( task );
  if ( data_number==OPTIONS_STABILIZATION )
    ok = check_unknowns_are_specified( task );
  if ( data_number==POST_LINE )
    ok = check_unknowns_are_specified( task );
  if ( data_number==POST_LINE_N )
    ok = check_unknowns_are_specified( task );
  if ( data_number==POST_LINE_DOF )
    ok = check_unknowns_are_specified( task );
  if ( data_number==POST_POINT )
    ok = check_unknowns_are_specified( task );
  if ( data_number==POST_POINT_DOF )
    ok = check_unknowns_are_specified( task );
  if ( data_number==POST_QUADRILATERAL ) {
    ok = check_ndim( 2, 3, task );
    ok = ok && check_unknowns_are_specified( task );
  }
  if ( data_number==POST_QUADRILATERAL_N ) {
    ok = check_ndim( 2, 3, task );
    ok = ok && check_unknowns_are_specified( task );
  }
  if ( data_number==POST_QUADRILATERAL_DOF ) {
    ok = check_ndim( 2, 3, task  );
    ok = ok && check_unknowns_are_specified( task );
  }
  if ( data_number==PRINT_FAILURE )
    ok = check_unknown( "materi_stress", YES, task );
  if ( data_number==RECTANGLE )
    ok = check_ndim( 2, 2, task );
  if ( data_number==ROTATION_X_AXIS ) {
    ok = check_unknown( "materi_velocity", YES, task );
    ok = ok && check_ndim( 3, 3, task );
  }
  if ( data_number==ROTATION_Y_AXIS ) {
    ok = check_unknown( "materi_velocity", YES, task );
    ok = ok && check_ndim( 3, 3, task );
  }
  if ( data_number==ROTATION_Z_AXIS ) {
    ok = check_unknown( "materi_velocity", YES, task );
    ok = ok && check_ndim( 2, 3, task );
  }
  if ( data_number==SLIDE_GEOMETRY )
    ok = check_unknown( "materi_velocity", YES, task );
  if ( data_number==SLIDE_PENALTY )
    ok = check_unknown( "materi_velocity", YES, task );
  if ( data_number==SLIDE_FRICTION )
    ok = check_unknown( "materi_velocity", YES, task );
  if ( data_number==SPRING1 ||  data_number==SPRING2 ) {
    if      ( !materi_displacement )
      ok = check_unknown( "materi_velocity_integrated", YES, task );
    else if ( !materi_velocity_integrated )
      ok = check_unknown( "materi_displacement", YES, task );
  }
  if ( data_number==TENDON ) {
    ok = check_unknown( "materi_stress", YES, task );
  }
  if ( data_number==TENDON_ELASTI ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_velocity", YES, task );
  }
  if ( data_number==TENDON_EXPANSION ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "condif_temperature", YES, task );
  }
  if ( data_number==TENDON_PLASTI ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_velocity", YES, task );
  }
  if ( data_number==TENDON_STRESS ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_velocity", YES, task );
  }
  if ( data_number==TENDON_STRESS_TIME ) {
    ok = check_unknown( "materi_stress", YES, task );
    ok = ok && check_unknown( "materi_velocity", YES, task );
  }
  if ( data_number==GROUP_VOLUME_FACTOR )
    ok = check_ndim( 1, 2, task );
  if ( data_number==VOLUME_FACTOR )
    ok = check_ndim( 1, 2, task );
  if ( data_number==WAVE ) {
    ok = check_unknown( "wave_scalar", YES, task );
    ok = ok && check_unknown( "wave_fscalar", YES, task );
    ok = ok && check_unknown( "materi_velocity", NO, task );
  }

  return ok;
}


long int check_ndim( long int lower, long int higher, long int task ) 

{
  long int ok=1;

  if ( ndim<lower || ndim>higher ) {
    if ( task==CHECK_USAGE_AND_ERROR ) {
      pri( "\n\n Error: inconsistent with number_of_space_dimensions." );
      exit(TN_EXIT_STATUS);
    }
    else {
      assert( task==CHECK_USAGE );
      ok = 0;
    }
  }

  return ok;
}

long int check_unknown(const char* str, long int initialization_needs_to_exist,
  long int task ) 

{
  long int iinitia=0, initialization_exists=0, ok=1, ninitia=0;

  ninitia = db_len( INITIALIZATION_VALUES, 0, VERSION_NORMAL );

  for ( iinitia=0; iinitia<ninitia; iinitia++ ) {
    if ( !strcmp(str,initialization_names[iinitia]) ) 
      initialization_exists = 1;
  }

  if ( initialization_needs_to_exist==YES && !initialization_exists ) {
    if ( task==CHECK_USAGE_AND_ERROR ) {
      cout << "\n\nError: " << str << " should be initialized.\n";
      exit(TN_EXIT_STATUS);
    }
    else {
      assert( task==CHECK_USAGE );
      ok = 0;
    }
  }
  else if ( initialization_needs_to_exist==NO && initialization_exists ) {
    if ( task==CHECK_USAGE_AND_ERROR ) {
      cout << "\n\nError: " << str << " should not be initialized.\n";
      exit(TN_EXIT_STATUS);
    }
    else {
      assert( task==CHECK_USAGE );
      ok = 0;
    }
  }

  return ok;
}

long int check_unknown_atleastone( const char* str1, const char* str2, long int task )

{
  long int iinitia=0, initialization1_exists=0, 
    initialization2_exists=0, ok=1, ninitia=0;

  ninitia = db_len( INITIALIZATION_VALUES, 0, VERSION_NORMAL );

  for ( iinitia=0; iinitia<ninitia; iinitia++ ) {
    if ( !strcmp(str1,initialization_names[iinitia]) ) 
      initialization1_exists = 1;
    if ( !strcmp(str2,initialization_names[iinitia]) ) 
      initialization2_exists = 1;
  }

  if ( !initialization1_exists && !initialization2_exists ) {
    if ( task==CHECK_USAGE_AND_ERROR ) {
      cout << "\n\nError:  at least one of " << str1 <<  " or " << 
        str2 << " should be initialized.\n";
      exit(TN_EXIT_STATUS);
    }
    else {
      assert( task==CHECK_USAGE );
      ok = 0;
    }
  }

  return ok;
}

long int check_unknown_minimum( char str[], long int min, long int task ) 

{
  long int ok=1, iinitia=0, ninitia=0, initialization_value=0,
    initialization_values[DATA_ITEM_SIZE];
  double ddum[1];

  db( INITIALIZATION_VALUES, 0, initialization_values, ddum, 
    ninitia, VERSION_NORMAL, GET );

  for ( iinitia=0; iinitia<ninitia; iinitia++ ) {
    if ( !strcmp(str,initialization_names[iinitia]) )
      initialization_value = initialization_values[iinitia];
  }

  if ( initialization_value<min ) {
    if ( task==CHECK_USAGE_AND_ERROR ) {
      cout << "\n\nError: " << str << " should be " << min << ".\n";
      exit(TN_EXIT_STATUS);
    }
    else {
      assert( task==CHECK_USAGE );
      ok = 0;
    }
  }

  return ok;
}


long int check_unknowns_are_specified( long int task )

{
  long int ok=1;

  if ( nuknwn==0 ) {
    if ( task==CHECK_USAGE_AND_ERROR ) {
      pri( "\n\nError: no unknowns (e.g. condif_temperature, ..) initialized");
      exit(TN_EXIT_STATUS);
    }
    else {
      assert( task==CHECK_USAGE );
      ok = 0;
    }
  }

  return ok;
}
