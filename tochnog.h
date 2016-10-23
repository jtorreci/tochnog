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
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/
#include <iostream>
#include <cassert>
#include <cctype>
#include <cfloat>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <ctime>
#include <inttypes.h>
#include <string.h>
#include "tnpetsc.h"
#include "tnsuplu.h"
#include "tnlapack.h"
#include "tnhypo.h"
#include "matrix.h"
#include "time.h"
#include "f2c.h"

using namespace std;
 
typedef long int integer;
typedef double doublereal; 

  // version stuff: NODE and ELEMENT data can have more versions (i.e. meshes)
enum {
  VERSION_NORMAL = 0,   // time=t
  VERSION_START,        // time=start_time
  VERSION_NEW,          // time=t+dt
  VERSION_TMP,          // trash version used by several routines
  VERSION_NEW_MESH_TMP, // trash version used by new_mesh
  VERSION_NEW_MESH_GENERATED, // generated mesh in new_mesh
  VERSION_PRINT,        // mesh for printing
  VERSION_MACRO,        // mesh for control_macro
  VERSION_EXTRUDE,      // mesh for extrude
  MVERSION              // maximum number of versions, this must be the last item
}; 

  // constants
#define MCHAR 100  // maximum length of names
#define MDIM 3 // maximum number of space dimensions
#define MNOL 27 // maximum number of nodes in an element, set if to 64 for HEX64
#define MSTRAIN 6 // maximum number of strain components
#define MTENDON 10 // maximum number of tendons in an element
#define MTYPE 10 // maximum number of types
#define MCALCUL 20 // maximum length of calcul records
#define MRANGE 500000 // maximum range length
#define MMAXWELL 50 // maximum number of maxwell chains
#define MTHREAD 64 // maximum number of threads
#define MAXIMUM_NODE 64 // always 64
#define DATA_ITEM_SIZE 190 // maximum length of (almost all) records
#define NONLOCAL_ITEM_SIZE 160 // maximum number of integration points for nonlocal calculation
#define TN_PRECISION 12 // precision in writing output file
#define MPOINT MNOL // maximum number of integration points in an element, always MNOL
#define MUKNWN DATA_ITEM_SIZE // maximum number of unknowns, always DATA_ITEM_SIZE
#define MPUKNWN MUKNWN // maximum number of primary unknowns, always MUKNWN
#define MPRINC 10 // maximum number of principal unknowns
#define MBOUNDA 1000 // maximum length bounda_unknown and bounda_force
#define MSTACK 10 // maximum lnumber of routines in routine stack
#define TN_EXIT_STATUS 1

const double EPS_ISO=1.e-9;  // epsilon for isoparametric coordinates
const double EPS_COORD=1.e-9; // epsilon for coordinates
const double EPS_SMALL=1.e-12; // epsilon generic
const double EPS_VOLUME=1.e-6; // epsilon for element volume
const double EPS_MATERI_DIFFUSION_MINIMUM=0.5; // epsilon diffusion if an element is empty
const double EPS_MATERI_DENSITY_MINIMUM=1.e-9; // epsilon density if an element is empty
const double PIRAD = 3.14159265358979323846; // Abramowitz and M_PI in GLIBC
const double NO_YIELD_F=-1.e10; // yield function value for no yielding
const double BOUNDARY_REDUCTION_FACTOR=2./3.; // reduction factor granular material - wall
const double TINY=1.e-20;

  // database items + others
enum {
  MINUS_ONE = 1,
  ABOVE,
  ABSOL,
  ADD,
  ADD_ALWAYS,
  ALL,
  ANY,
  AREA,
  AREA_ELEMENT_GROUP,
  AREA_ELEMENT_GROUP_METHOD,
  AREA_ELEMENT_GROUP_SEQUENCE,
  AREA_ELEMENT_GROUP_SEQUENCE_ELEMENT,
  AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP,
  AREA_ELEMENT_GROUP_SEQUENCE_GEOMETRY,
  AREA_ELEMENT_GROUP_SEQUENCE_METHOD,
  AREA_ELEMENT_GROUP_SEQUENCE_TIME,
  AREA_NODE_DATAITEM,
  AREA_NODE_DATAITEM_DOUBLE,
  AREA_NODE_DATAITEM_INTEGER,
  ASM,
  AVERAGE,
  BAR,
  BAR2,
  BAR3,
  BAR4,
  BEAM,
  BEAM_ROTATION,
  BELOW,
  BCGS,
  BICG,
  BJACOBI,
  BOUNDA,
  BOUNDA_FORCE,
  BOUNDA_SINE,
  BOUNDA_TIME,
  BOUNDA_TIME_FILE,
  BOUNDA_TIME_USER,
  BOUNDA_UNKNOWN,
  BRICK,
  CALCULATE_STRESSINTENSITYFACTOR,
  CG,
  CGS,
  CHANGE,
  CHANGE_DATAITEM,
  CHANGE_DATAITEM_TIME,
  CHANGE_DATAITEM_TIME_DISCRETE,
  CHANGE_DATAITEM_TIME_USER,
  CHANGE_GEOMETRY,
  CHANGE_GEOMETRY_TIME_USER,
  CHEBYCHEV,
  CHECK,
  CHECK_COMBINATION,
  CHECK_INDEX,
  CHECK_NUMBER,
  CHECK_USAGE,
  CHECK_USAGE_AND_ERROR,
  CIRCLE,
  CIRCLE_HOLLOW,
  CLOSE,
  COMPOSITE,
  CONDIF,
  CONDIF_CONVECTION,
  CONDIF_CONVECTION_GEOMETRY,
  CONDIF_RADIATION,
  CONDIF_RADIATION_GEOMETRY,
  CONDIF_TEMPERATURE,
  CONTACT,
  CONTACTSPRING,
  CONTACT_FRICTION,
  CONTACT_GEOMETRY,
  CONTACT_GEOMETRY_SWITCH,
  CONTACT_HEATGENERATION,
  CONTACT_PENALTY_PRESSURE,
  CONTACT_PENALTY_TEMPERATURE,
  CONTACT_PENALTY_VELOCITY,
  CONTACT_RELAXATION,
  CONTACT_STICK,
  CONTROL,
  CONTROL_CRACK,
  CONTROL_DATA_DELETE,
  CONTROL_DATA_INITELDOF_GEOMETRY,
  CONTROL_DATA_PUT,
  CONTROL_DATA_PUT_DOUBLE,
  CONTROL_DATA_PUT_DOUBLE_NODE,
  CONTROL_DATA_PUT_INTEGER,
  CONTROL_DISTRIBUTE,
  CONTROL_DISTRIBUTE_VALUES,
  CONTROL_EIGEN,
  CONTROL_EIGEN_SCALE,
  CONTROL_EIGEN_VALUES,
  CONTROL_MATERI_DIFFUSION,
  CONTROL_MESH_ADJUST_GEOMETRY,
  CONTROL_MESH_DELETE_GEOMETRY,
  CONTROL_MESH_DELETE_GEOMETRY_ELEMENT,
  CONTROL_MESH_DELETE_GEOMETRY_ELEMENTGROUP,
  CONTROL_MESH_DELETE_GEOMETRY_FACTOR,
  CONTROL_MESH_DELETE_GEOMETRY_MOVENODES,
  CONTROL_MESH_DELETE_SMALL,
  CONTROL_MESH_EXTRUDE,
  CONTROL_MESH_EXTRUDE_N,
  CONTROL_MESH_GENERATE_BEAM,
  CONTROL_MESH_GENERATE_CONTACTSPRING,
  CONTROL_MESH_GENERATE_CONTACTSPRING_ELEMENT,
  CONTROL_MESH_GENERATE_SPRING1,
  CONTROL_MESH_GENERATE_SPRING2,
  CONTROL_MESH_GENERATE_TRUSS,
  CONTROL_MESH_GENERATE_TRUSSBEAM,
  CONTROL_MESH_GENERATE_TRUSS_BEAM_LOOSE,
  CONTROL_MESH_GENERATE_TRUSS_BEAM_MACRO,
  CONTROL_MESH_MACRO,
  CONTROL_MESH_MACRO_ELEMENT,
  CONTROL_MESH_MACRO_PARAMETERS,
  CONTROL_MESH_MACRO_SET_NODE_BOUNDARY,
  CONTROL_MESH_MERGE,
  CONTROL_MESH_MERGE_EPSCOORD,
  CONTROL_MESH_MERGE_MACRO_GENERATE,
  CONTROL_MESH_MERGE_NOT,
  CONTROL_MESH_NEW_MESH,
  CONTROL_MESH_NEW_MESH_ELEMENT,
  CONTROL_MESH_NEW_MESH_REGION,
  CONTROL_MESH_REFINE_GLOBALLY,
  CONTROL_MESH_REFINE_GLOBALLY_GEOMETRY,
  CONTROL_MESH_REFINE_LOCALLY,
  CONTROL_MESH_REFINE_LOCALLY_GEOMETRY,
  CONTROL_MESH_REFINE_LOCALLY_NOT,
  CONTROL_MESH_REFINE_LOCALLY_ONLY,
  CONTROL_MESH_REFINE_LOCALLY_UNKNOWN,
  CONTROL_MESH_REMESH,
  CONTROL_MESH_REMESH_FACTOR,
  CONTROL_MESH_RENUMBER,
  CONTROL_MESH_SPLIT,
  CONTROL_MESH_SPLIT_ONLY,
  CONTROL_OPTIONS_CONVECTION,
  CONTROL_OPTIONS_INERTIA,
  CONTROL_OPTIONS_RELAXATION,
  CONTROL_OPTIONS_SKIP_GRAVITY,
  CONTROL_OPTIONS_SKIP_GROUNDFLOW_MATERIDIVERGENCE,
  CONTROL_OPTIONS_SKIP_GROUNDFLOW_NONLINEAR,
  CONTROL_OPTIONS_SKIP_PLASTICITY,
  CONTROL_OPTIONS_SOLVER,
  CONTROL_OPTIONS_SOLVER_BICG_ERROR,
  CONTROL_OPTIONS_SOLVER_BICG_ERROR_MINIMUM,
  CONTROL_OPTIONS_SOLVER_PETSC_KSPTYPE,
  CONTROL_OPTIONS_SOLVER_PETSC_MG,
  CONTROL_OPTIONS_SOLVER_PETSC_PCTYPE,
  CONTROL_PRINT,
  CONTROL_PRINT_DATABASE,
  CONTROL_PRINT_DATA_VERSUS_DATA,
  CONTROL_PRINT_DX,
  CONTROL_PRINT_DX_TIME,
  CONTROL_PRINT_ELEMENT,
  CONTROL_PRINT_FILTER,
  CONTROL_PRINT_GID,
  CONTROL_PRINT_GID_EMPTY,
  CONTROL_PRINT_GID_TIME,
  CONTROL_PRINT_GID_MESH,
  CONTROL_PRINT_GMV,
  CONTROL_PRINT_GMV_MESH,
  CONTROL_PRINT_HISTORY,
  CONTROL_PRINT_MATLAB,
  CONTROL_PRINT_PLOTMTV,
  CONTROL_PRINT_PLOTMTV_MESH,
  CONTROL_PRINT_TECPLOT,
  CONTROL_PRINT_TECPLOT_MESH,
  CONTROL_PRINT_UNKNOWNS,
  CONTROL_PRINT_UNKNOWNSRHSIDE,
  CONTROL_PRINT_VTK,
  CONTROL_RELAXATION_CONDIF_TEMPERATURE,
  CONTROL_RELAXATION_GROUNDFLOW_PRESSURE,
  CONTROL_RELAXATION_MATERI_VELOCITY,
  CONTROL_RELAXATION_MAXWELL_E,
  CONTROL_RELAXATION_WAVE_FSCALAR,
  CONTROL_REPEAT,
  CONTROL_REPEAT_UNTIL_ITEM,
  CONTROL_REPEAT_UNTIL_TOLERANCE,
  CONTROL_REPEAT_UNTIL_VALUE,
  CONTROL_RESTART,
  CONTROL_UNKNOWN_FREEZE,
  CONTROL_UNKNOWN_RESET_UNKNOWN,
  CONTROL_UNKNOWN_RESET_GEOMETRY,
  CONTROL_UNKNOWN_RESET_VALUE,
  CONTROL_TIMESTEP,
  CONTROL_TIMESTEP_ITERATIONS,
  CONTROL_TIMESTEP_ITERATIONS_AUTOMATIC,
  CONTROL_TIMESTEP_ITERATIONS_AUTOMATIC_STOP,
  CONTROL_TIMESTEP_SIZE_AUTOMATIC_DECREASE,
  CONTROL_TIMESTEP_MULTIPLIER,
  CR,
  CRACK,
  CRACK_DIRECTION,
  CRACK_ELEMENTGROUP,
  CRACK_LENGTH,
  CRACK_NODES,
  CRACK_STRESSINTENSITYFACTOR,
  CRACK_TIP,
  CYLINDER_HOLLOW,
  DATABASE,
  DEPENDENCY,
  DEPENDENCY_DIAGRAM,
  DEPENDENCY_ITEM,
  DIAGONAL,
  DOF,
  DOF_AMOUNT,
  DOF_LABEL,
  DOF_PRINCIPAL,
  DOF_SCAL_VEC_MAT,
  DOF_TYPE,
  DOUBLE_PRECISION,
  DTIME,
  DYNAMIC,
  EISENSTAT,
  ELEMENT,
  ELEMENT_ADJUST,
  ELEMENT_BEAM_DIRECTION,
  ELEMENT_BEAM_MOMENT,
  ELEMENT_CONTACTSPRING_DIRECTION,
  ELEMENT_CONTACTSPRING_FORCE,
  ELEMENT_DELETE_FACTOR,
  ELEMENT_DELETE_TIMES,
  ELEMENT_DISTRIBUTE,
  ELEMENT_DISTRIBUTE_VALUES,
  ELEMENT_DOF,
  ELEMENT_DOF_INITIALISED,
  ELEMENT_EMPTY,
  ELEMENT_GROUP,
  ELEMENT_GROUP_AREA_ELEMENT_GROUP,
  ELEMENT_GROUP_AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP,
  ELEMENT_MACRO_GENERATE,
  ELEMENT_MASS,
  ELEMENT_MATRIX_DELETE,
  ELEMENT_MATRIX_SECOND_VALUES,
  ELEMENT_MATRIX_UNKNOWNS,
  ELEMENT_MATRIX_VALUES,
  ELEMENT_MIDDLE,
  ELEMENT_NONLOCAL,
  ELEMENT_NONLOCAL_IPOINT,
  ELEMENT_NONLOCAL_WEIGHT,
  ELEMENT_RADIUS,
  ELEMENT_RHSIDE_DELETE,
  ELEMENT_SPRING_DIRECTION,
  ELEMENT_SPRING_FORCE,
  ELEMENT_STRAINENERGY,
  ELEMENT_TENDON_DIRECTION,
  ELEMENT_TENDON_INTERSECTIONS,
  ELEMENT_TENDON_NUMBER,
  ELEMENT_TENDON_STRAIN,
  ELEMENT_TENDON_STRESS,
  ELEMENT_TENDON_VOLUME,
  ELEMENT_TRUSS_DIRECTION,
  ELEMENT_TRUSS_FORCE,
  ELEMENT_VOLUME,
  EMPTY,
  EVERYTHING,
  EXIT_TOCHNOG,
  FIXED_IN_SPACE,
  FOLLOW_MATERIAL,
  FORCE,
  FORCE_ELEMENT_EDGE,
  FORCE_ELEMENT_EDGE_FACTOR,
  FORCE_ELEMENT_EDGE_GEOMETRY,
  FORCE_ELEMENT_EDGE_SINE,
  FORCE_ELEMENT_EDGE_TIME,
  FORCE_ELEMENT_EDGE_TIME_FILE,
  FORCE_ELEMENT_EDGE_NORMAL,
  FORCE_ELEMENT_EDGE_NORMAL_FACTOR,
  FORCE_ELEMENT_EDGE_NORMAL_GEOMETRY,
  FORCE_ELEMENT_EDGE_NORMAL_SINE,
  FORCE_ELEMENT_EDGE_NORMAL_TIME,
  FORCE_ELEMENT_EDGE_WATER,
  FORCE_ELEMENT_EDGE_WATER_GEOMETRY,
  FORCE_ELEMENT_EDGE_WATER_TIME,
  FORCE_ELEMENT_VOLUME,
  FORCE_ELEMENT_VOLUME_FACTOR,
  FORCE_ELEMENT_VOLUME_GEOMETRY,
  FORCE_ELEMENT_VOLUME_SINE,
  FORCE_ELEMENT_VOLUME_TIME,
  FORCE_GRAVITY,
  FORCE_GRAVITY_TIME,
  FROM,
  FRONT,
  GAUSS,
  GEOMETRY,
  GEOMETRY_BOUNDA_FACTOR,
  GEOMETRY_BOUNDA_SINE_X,
  GEOMETRY_BOUNDA_SINE_Y,
  GEOMETRY_BOUNDA_SINE_Z,
  GEOMETRY_BRICK,
  GEOMETRY_CIRCLE,
  GEOMETRY_CIRCLE_SEGMENT,
  GEOMETRY_CIRCLE_SMALLSEGMENT,
  GEOMETRY_CYLINDER,
  GEOMETRY_CYLINDER_SEGMENT,
  GEOMETRY_ELLIPSE,
  GEOMETRY_LINE,
  GEOMETRY_NUMBER,
  GEOMETRY_POINT,
  GEOMETRY_POLYNOMIAL,
  GEOMETRY_QUADRILATERAL,
  GEOMETRY_SET,
  GEOMETRY_SPHERE,
  GEOMETRY_SPHERE_SEGMENT,
  GEOMETRY_TRIANGLE,
  GEOMETRY_TRIANGLE_EPSISO,
  GENERALIZED,
  GET,
  GET_AND_CHECK,
  GET_FLOW_RULE,
  GET_FLOW_RULE_GRAD,
  GET_IF_EXISTS,
  GET_YIELD_RULE,
  GLOBAL_ELEMENTS,
  GLOBAL_MASS,
  GLOBAL_NODES,
  GLOBAL_POINT_MATERI_DIFFUSION_LOST,
  GLOBAL_POINT_MATERI_DIFFUSION_TOTAL,
  GLOBAL_SOLVER_ERROR,
  GLOBAL_SOLVER_ITERATIONS,
  GLOBAL_STRAINENERGY,
  GLOBAL_UNKNOWN_AVERAGE,
  GLOBAL_UNKNOWN_MAX,
  GLOBAL_UNKNOWN_MIN,
  GLOBAL_UNKNOWN_NUMBER,
  GLOBAL_UNKNOWN_SUM,
  GLOBAL_VOLUME,
  GMRES,
  GROUND,
  GROUNDFLOW,
  GROUNDFLOW_ADDTOPRESSURE,
  GROUNDFLOW_PHREATICLEVEL,
  GROUNDFLOW_PHREATICLEVEL_N,
  GROUNDFLOW_PHREATICLEVEL_BOUNDA,
  GROUNDFLOW_PHREATICLEVEL_MINIMUM,
  GROUNDFLOW_DENSITY,
  GROUNDFLOW_PRESSURE,
  GROUNDFLOW_PRESSURE_ATMOSPHERIC,
  GROUNDFLOW_VELOCITY,
  GROUP_AXISYMMETRIC,
  GROUP_BEAM_AREA,
  GROUP_BEAM_YOUNG,
  GROUP_BEAM_INERTIA,
  GROUP_BEAM_MEMORY,
  GROUP_BEAM_PLANE,
  GROUP_CONDIF_ABSORPTION,
  GROUP_CONDIF_CAPACITY,
  GROUP_CONDIF_CONDUCTIVITY,
  GROUP_CONDIF_DENSITY,
  GROUP_CONDIF_FLOW,
  GROUP_CONTACTSPRING_COHESION,
  GROUP_CONTACTSPRING_DIRECTION,
  GROUP_CONTACTSPRING_FRICTION,
  GROUP_CONTACTSPRING_FRICTION_AUTOMATIC,
  GROUP_CONTACTSPRING_MEMORY,
  GROUP_CONTACTSPRING_STIFFNESS,
  GROUP_GROUNDFLOW_CAPACITY,
  GROUP_GROUNDFLOW_CAPACITY_NONLINEAR_METHOD,
  GROUP_GROUNDFLOW_CAPACITY_NONLINEAR_PARAMETERS,
  GROUP_GROUNDFLOW_MATERIDIVERGENCE,
  GROUP_GROUNDFLOW_PERMEABILITY,
  GROUP_GROUNDFLOW_PERMEABILITY_NONLINEAR_METHOD,
  GROUP_GROUNDFLOW_PERMEABILITY_NONLINEAR_PARAMETERS,
  GROUP_GROUNDFLOW_POROSITY,
  GROUP_INTEGRATION_METHOD,
  GROUP_INTEGRATION_POINTS,
  GROUP_MATERI_DAMAGE_MAZARS,
  GROUP_MATERI_DAMPING,
  GROUP_MATERI_DENSITY,
  GROUP_MATERI_DENSITY_GROUNDFLOW,
  GROUP_MATERI_ELASTI_CAMCLAY_G,
  GROUP_MATERI_ELASTI_CAMCLAY_POISSON,
  GROUP_MATERI_ELASTI_COMPRESSIBILITY,
  GROUP_MATERI_ELASTI_LADE,
  GROUP_MATERI_ELASTI_POISSON,
  GROUP_MATERI_ELASTI_SMALLSTRAIN,
  GROUP_MATERI_ELASTI_TRANSVERSE_ISOTROPY,
  GROUP_MATERI_ELASTI_TRANSVERSE_ISOTROPY_GRAHOUL,
  GROUP_MATERI_ELASTI_TSKH,
  GROUP_MATERI_ELASTI_VOLUMETRIC_POISSON,
  GROUP_MATERI_ELASTI_VOLUMETRIC_YOUNG_ORDER,
  GROUP_MATERI_ELASTI_VOLUMETRIC_YOUNG_VALUES,      
  GROUP_MATERI_ELASTI_YOUNG,
  GROUP_MATERI_ELASTI_YOUNG_POLYNOMIAL,
  GROUP_MATERI_ELASTI_YOUNG_POWER,
  GROUP_MATERI_ELASTI_YOUNG_STRAINSTRESS,
  GROUP_MATERI_EXPANSION_LINEAR,
  GROUP_MATERI_EXPANSION_VOLUME,
  GROUP_MATERI_FAILURE_CRUCHING,
  GROUP_MATERI_FAILURE_DAMAGE,
  GROUP_MATERI_FAILURE_PLASTI_KAPPA,
  GROUP_MATERI_FAILURE_RUPTURE,
  GROUP_MATERI_FAILURE_VOIDFRACTION,
  GROUP_MATERI_HYPER_BESSELING,
  GROUP_MATERI_HYPER_BLATZ_KO,
  GROUP_MATERI_HYPER_MOONEY_RIVLIN,
  GROUP_MATERI_HYPER_NEOHOOKEAN,
  GROUP_MATERI_HYPER_REDUCEDPOLYNOMIAL,
  GROUP_MATERI_HYPER_STIFFNESS,
  GROUP_MATERI_HYPER_VOLUMETRIC_OGDEN,
  GROUP_MATERI_HYPER_VOLUMETRIC_LINEAR,
  GROUP_MATERI_HYPER_VOLUMETRIC_MURNAGHAN,
  GROUP_MATERI_HYPER_VOLUMETRIC_POLYNOMIAL,
  GROUP_MATERI_HYPER_VOLUMETRIC_SIMOTAYLOR,
  GROUP_MATERI_ISOTROPY,
  GROUP_MATERI_MAXWELL_CHAIN,
  GROUP_MATERI_MAXWELL_CHAIN_NONLINEAR,
  GROUP_MATERI_MEMBRANE,
  GROUP_MATERI_MEMORY,
  GROUP_MATERI_PLASTI_AITSKH,
  GROUP_MATERI_PLASTI_BOUNDARY,
  GROUP_MATERI_PLASTI_BOUNDARY_FACTOR,
  GROUP_MATERI_PLASTI_CAMCLAY,
  GROUP_MATERI_PLASTI_CAMCLAY_INCREMENTAL,
  GROUP_MATERI_PLASTI_CAP,
  GROUP_MATERI_PLASTI_COMPRESSION,
  GROUP_MATERI_PLASTI_DIPRISCO,
  GROUP_MATERI_PLASTI_DIPRISCO_RT,
  GROUP_MATERI_PLASTI_DRUCKPRAG,
  GROUP_MATERI_PLASTI_DRUCKPRAG_APEX,
  GROUP_MATERI_PLASTI_DRUCKPRAG_TENSIONCUTOFF,
  GROUP_MATERI_PLASTI_DRUCKPRAG_TENSIONLIMIT,
  GROUP_MATERI_PLASTI_GURSON,
  GROUP_MATERI_PLASTI_HLC,
  GROUP_MATERI_PLASTI_HEATGENERATION,
  GROUP_MATERI_PLASTI_HYPO_LOWANGLES,
  GROUP_MATERI_PLASTI_HYPO_COHESION,
  GROUP_MATERI_PLASTI_HYPO_INTERGRANULARSTRAIN,
  GROUP_MATERI_PLASTI_HYPO_PRESSUREDEPENDENTVOIDRATIO,
  GROUP_MATERI_PLASTI_HYPO_WOLFERSDORFF,
  GROUP_MATERI_PLASTI_INCREMENTAL_ERASERECENTHISTORY,
  GROUP_MATERI_PLASTI_INCREMENTAL_FEERROR,
  GROUP_MATERI_PLASTI_INCREMENTAL_FESUBSTEPS,
  GROUP_MATERI_PLASTI_INCREMENTAL_MAXSUBSTEPS,
  GROUP_MATERI_PLASTI_INCREMENTAL_MINSUBSTEPS,
  GROUP_MATERI_PLASTI_INCREMENTAL_PRINTF,
  GROUP_MATERI_PLASTI_INCREMENTAL_USEMATRIX,
  GROUP_MATERI_PLASTI_KINEMATIC_HARDENING,
  GROUP_MATERI_PLASTI_MATSUOKANAKAI,
  GROUP_MATERI_PLASTI_MATSUOKANAKAI_APEX,
  GROUP_MATERI_PLASTI_MATSUOKANAKAI_TENSIONCUTOFF,
  GROUP_MATERI_PLASTI_MOHRCOUL,
  GROUP_MATERI_PLASTI_MOHRCOUL_01,
  GROUP_MATERI_PLASTI_MOHRCOUL_12,
  GROUP_MATERI_PLASTI_MOHRCOUL_20,
  GROUP_MATERI_PLASTI_MOHRCOUL_APEX,
  GROUP_MATERI_PLASTI_MOHRCOUL_SOFTENING,
  GROUP_MATERI_PLASTI_MOHRCOUL_TENSIONCUTOFF,
  GROUP_MATERI_PLASTI_STRESS,
  GROUP_MATERI_PLASTI_TENSION,
  GROUP_MATERI_PLASTI_TSKH,
  GROUP_MATERI_PLASTI_USER,
  GROUP_MATERI_PLASTI_VISCO_ALWAYS,
  GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL,
  GROUP_MATERI_PLASTI_VISCO_POWER,
  GROUP_MATERI_PLASTI_VONMISES,
  GROUP_MATERI_PLASTI_VONMISES_NADAI,
  GROUP_MATERI_STOKES,
  GROUP_MATERI_VISCOSITY,
  GROUP_MATERI_VISCOSITY_HEATGENERATION,
  GROUP_MATERI_VISCOSITY_USER,
  GROUP_MATRIX_SECOND_VALUES,
  GROUP_MATRIX_UNKNOWNS,
  GROUP_MATRIX_VALUES,
  GROUP_MAXWELL_FREQUENCY_PML_EPSILONANDMU,
  GROUP_MAXWELL_FREQUENCY_PML_PLANES,     
  GROUP_MAXWELL_FREQUENCY_EIGEN,
  GROUP_MAXWELL_FREQUENCY_EPSILON,
  GROUP_MAXWELL_FREQUENCY_EPSILON_ANISOTROPIC,
  GROUP_MAXWELL_FREQUENCY_J,
  GROUP_MAXWELL_FREQUENCY_MU,
  GROUP_MAXWELL_FREQUENCY_PENALTY,
  GROUP_MAXWELL_TIME_EPSILON,
  GROUP_MAXWELL_TIME_EPSILON_ANISOTROPIC,
  GROUP_MAXWELL_TIME_J,
  GROUP_MAXWELL_TIME_MU,
  GROUP_MAXWELL_TIME_PENALTY,
  GROUP_SPRING_DIRECTION,
  GROUP_SPRING_STIFFNESS,
  GROUP_SPRING_PLASTI,
  GROUP_TIME,
  GROUP_TRUSS_AREA,
  GROUP_TRUSS_ROPE,
  GROUP_TRUSS_DENSITY,
  GROUP_TRUSS_MEMORY,
  GROUP_TRUSS_PLASTI,
  GROUP_TRUSS_YOUNG,
  GROUP_TYPE,
  GROUP_USER_DATA,
  GROUP_USER_UMAT,
  GROUP_VOLUME_FACTOR,
  GROUP_WAVE_SPEED_OF_SOUND,
  GROWTH,
  HEX27,
  HEX64,
  HEX8,
  H_REFINEMENT,
  ICC,
  ICONTROL,
  ILU,
  INCREMENTAL,
  INITIALIZE,
  INITIALIZATION_VALUES,
  INTEGER,
  INVERSE,
  INVERSE_DETERMINE_NEW_ESTIMATES,
  INVERSE_DETERMINE_SENSITIVITY,
  INVERSE_HISTORY,
  INVERSE_ITERATIONS,
  INVERSE_ITERATION_NUMBER,
  INVERSE_PARAMETER,
  INVERSE_PARAMETER_LIMITS,
  INVERSE_PARAMETER_SENSITIVITY,
  INVERSE_PARAMETER_SET,
  INVERSE_PARAMETER_STEP,
  INVERSE_PARAMETER_VARIATION,
  INVERSE_TARGET,
  INVERSE_TARGET_DATA,
  INVERSE_TARGET_TIMESTEP,
  JACOBI,
  LOBATTO,
  LSQR,
  LU,
  MACRO,
  MATERI,
  MATERI_DAMAGE,
  MATERI_DENSITY,
  MATERI_DENSITY_MINIMUM,
  MATERI_DIFFUSION,
  MATERI_DIFFUSION_CORRECT,
  MATERI_DIFFUSION_MINIMUM,
  MATERI_DIFFUSION_ADJUST_GEOMETRY,
  MATERI_DIFFUSION_FILL_GEOMETRY,
  MATERI_DIFFUSION_FILL_EPSVELOCITY,
  MATERI_DIFFUSION_SMOOTH,
  MATERI_DIFFUSION_TEMPERATURE,
  MATERI_DISPLACEMENT,
  MATERI_HISTORY_VARIABLES,
  MATERI_MAXWELL_STRESS,
  MATERI_PLASTI_F,
  MATERI_PLASTI_F_NONLOCAL,
  MATERI_PLASTI_INCREMENTAL_SUBSTEPS,
  MATERI_PLASTI_KAPPA,
  MATERI_PLASTI_RHO,
  MATERI_PLASTI_SOFTVAR_LOCAL,
  MATERI_PLASTI_SOFTVAR_NONLOCAL,
  MATERI_ROTATION,
  MATERI_STRAINENERGY,
  MATERI_STRAIN_ELASTI,
  MATERI_STRAIN_INTERGRANULAR,
  MATERI_STRAIN_PLASTI,
  MATERI_STRAIN_TOTAL,
  MATERI_STRESS,
  MATERI_VELOCITY,
  MATERI_VELOCITY_INTEGRATED,
  MATERI_VOID_FRACTION,
  MATERI_WORK,
  MATRIX,
  MATRIX_ITERATIVE_BICG,
  MATRIX_ITERATIVE_PETSC,
  MATRIX_SUPERLU,
  MATRIX_SUPERLU_DIST,
  MATRIX_SUPERLU_MT,
  MATRIX_LAPACK,
  MATRIX_LAPACK_GEN,
  MAXFRE,
  MAXIMAL,
  MAXTIM,
  MAXWELL,
  MAXWELL_ECOMPLEX,
  MAXWELL_E,
  MAXWELL_FE,
  MAXWELL_EI,
  MAXWELL_ER,
  MAXWELL_FREQUENCY,
  MAXWELL_FREQUENCY_EXCITATION,
  MAXWELL_SCATTER_ENERGYCONSERVATION,
  MAXWELL_SCATTER_MATRIX_AMPLITUDE,
  MAXWELL_SCATTER_MATRIX_AMPLITUDEDB,
  MAXWELL_SCATTER_MATRIX_IMAGINARY,
  MAXWELL_SCATTER_MATRIX_PHASE,
  MAXWELL_SCATTER_MATRIX_REAL,
  MAXWELL_SCATTER_PARAMETERS,
  MAXWELL_SCATTER_PORT_INPUT,
  MAXWELL_SCATTER_PORT_OUTPUT,
  MAXWELL_TIME,
  MESH,
  METHOD1,
  METHOD2,
  MG,
  MOMENT,
  MINIMAL,
  MISES,
  NEGATIVE,
  NO,
  NODE,
  NODE_ADJUST,
  NODE_BOUNDARY,
  NODE_BOUNDED,
  NODE_DAMPING,
  NODE_DOF,
  NODE_DOF_CALCUL,
  NODE_DOF_MATERI_DIFFUSION,
  NODE_DOF_START_REFINED,
  NODE_DOF_TMP,
  NODE_EIGEN,
  NODE_ELEMENT,
  NODE_LHSIDE,
  NODE_MACRO_GENERATE,
  NODE_MASS,
  NODE_NEL,
  NODE_NODE,
  NODE_NONLOCAL,
  NODE_NONLOCAL_WEIGHT,
  NODE_PHREATICLEVEL,
  NODE_PRINT,
  NODE_REMESH_ALLOWED,
  NODE_REMESH_VELOCITY,
  NODE_RHSIDE,
  NODE_RHSIDE_PRINT,
  NODE_SET,
  NODE_START_REFINED,
  NODE_STIFFNESS,
  NONE,
  NONLOCAL_ELEMENT_INFO,
  NORMAL,
  NOTHING,
  NUMBER,
  NUMBER_ITERATIONS,
  ON,
  OPTIONS_CONVECTION,
  OPTIONS_ELEMENT_DOF,
  OPTIONS_ELEMENTLOOP,
  OPTIONS_INERTIA,
  OPTIONS_MATRIX_GROUP,
  OPTIONS_MATRIX_LENGTH,
  OPTIONS_MESH,
  OPTIONS_NONLOCAL,
  OPTIONS_NONLOCAL_SOFTVAR,
  OPTIONS_PROCESSORS,
  OPTIONS_RELAXATION,
  OPTIONS_RESIDUEFACTOR,
  OPTIONS_SKIP_GRAVITY,
  OPTIONS_SKIP_GROUNDFLOW_MATERIDIVERGENCE,
  OPTIONS_SKIP_GROUNDFLOW_NONLINEAR,
  OPTIONS_SKIP_PLASTICITY,
  OPTIONS_SOLVER,
  OPTIONS_SOLVER_BICG_ERROR,
  OPTIONS_SOLVER_BICG_ERROR_MINIMUM,
  OPTIONS_STABILIZATION,
  PCG,
  PHIMOB,
  PLUS_DISPLACEMENT,
  POINT_MATERI_DIFFUSION,
  POINT_MATERI_DIFFUSION_PREVIOUS,
  POSITIVE,
  POST,
  POST_CALCUL,
  POST_CALCUL_SCAL_VEC_MAT,
  POST_CALCUL_UNKNOWN_OPERAT,
  POST_ERROR_ITEM,
  POST_ERROR_MESH1,
  POST_ERROR_MESH2,
  POST_ERROR_RESULT,
  POST_GLOBAL,
  POST_INTEGRATE,
  POST_INTEGRATE_RESULT,
  POST_LINE,
  POST_LINE_DOF,
  POST_LINE_DOF_CALCUL,
  POST_LINE_MOMENT,
  POST_LINE_N,
  POST_LINE_OPERAT,
  POST_NODE,
  POST_NODE_RESULT,
  POST_NODE_RHSIDE_FIXED,
  POST_NODE_RHSIDE_FREE,
  POST_NODE_RHSIDE_RATIO,
  POST_NODE_RHSIDE_RATIO_UNKNOWNTYPES,
  POST_POINT,
  POST_POINT_DOF,
  POST_POINT_DOF_CALCUL,
  POST_QUADRILATERAL,
  POST_QUADRILATERAL_DOF,
  POST_QUADRILATERAL_DOF_CALCUL,
  POST_QUADRILATERAL_N,
  PREONLY,
  PRINT,
  PRINT_ARITHMETIC,
  PRINT_CONTROL,
  PRINT_DEFINE,
  PRINT_LASTDATABASE,
  PRINT_FAILURE,
  PRINT_FILTER,
  PRINT_SOLVER,
  PRINT_WHERE,
  PRISM6,
  PRIVAL,
  PRIVEC,
  PROJECT_EXACT,
  PROJECT_ON_EDGE,
  PUT,
  P_COARSEN,
  P_REFINEMENT,
  RICHARDSON,
  QUAD16,
  QUAD4,
  QUAD9,
  RA,
  RECTANGLE,
  REMOVE,
  RESIDUE,
  RESTART,
  ROTATION_X_AXIS,
  ROTATION_Y_AXIS,
  ROTATION_Z_AXIS,
  SCALAR,
  SEPARATE,
  SEPARATE_INDEX,
  SEPARATE_SEQUENTIAL,
  SHELL,
  SIZEDEV,
  SIZETOT,
  SLES,
  SLIDE,
  SLIDE_FRICTION,
  SLIDE_GEOMETRY,
  SLIDE_PENALTY,
  SOR,
  SPHERE,
  SPRING,
  SPRING1,
  SPRING2,
  START,
  STATIC,
  STEP,
  STRESS,
  SUM,
  TARGET,
  TARGET_ITEM,
  TARGET_VALUE,
  TCQMR,
  TENDON,
  TENDON_ELASTI,
  TENDON_EXPANSION,
  TENDON_PLASTI,
  TENDON_SPLIT,
  TENDON_SPLIT_ELEMENT,
  TENDON_STRESS,
  TENDON_STRESS_TIME,
  TET10,
  TET4,
  TFQMR,
  THERMAL,
  TIME,
  TIME_AT_START,
  TIME_CALCULATION,
  TIME_CURRENT,
  TIME_NEW,
  TIME_OLD,
  TO,
  TOTAL,
  TOTAL_PIOLA,
  TOTAL_LINEAR,
  TRIA3,
  TRIA6,
  TRUSS,
  TRUSSBEAM,
  UNIFORM,
  UPDATED,
  UPDATED_WITHOUT_ROTATION,
  USE,
  USER,
  VALUE,
  VECTOR,
  VOLUME,
  VOLUME_FACTOR,
  WAVE,
  WAVE_FSCALAR,
  WAVE_SCALAR,
  YES,
  X,
  Y,
  Z,
  LAST_DUMMY }; // keep LAST_DUMMY always the last one

#define MDAT LAST_DUMMY+DATA_ITEM_SIZE  // reserve space for unknowns

  // see initialization part in manual
extern long int echo, ndim, derivatives, 
  beam_rotation, condif_temperature, 
  groundflow_velocity, groundflow_pressure, materi_history_variables,
  materi_damage, materi_density, 
  materi_diffusion, materi_displacement, 
  materi_maxwell_stress, materi_plasti_kappa, 
  materi_plasti_f, materi_plasti_f_nonlocal, materi_plasti_incremental_substeps,
  materi_plasti_softvar_local,
  materi_plasti_softvar_nonlocal,
  materi_strain_intergranular, materi_plasti_rho, materi_strainenergy,
  materi_strain_elasti, materi_strain_plasti, materi_strain_total, 
  materi_stress, materi_velocity, materi_velocity_integrated,
  materi_void_fraction, materi_work, 
  maxwell_e, maxwell_fe,
  maxwell_er, maxwell_ei,
  residue, wave_scalar, wave_fscalar,
  find_local_softvar, find_nonlocal_weights, nonlocal_first_set; 

  // for internal use
extern long int
  any_runtime, // 1 if runtime file found
  nder, // number of derivatives (incl. primary unknown itself)
  npuknwn, // number of primary unknowns
  nuknwn, // number of unknowns (= npuknwn*nder)
  nprinc, // number of principal unknowns
  npoint, // number of integration points in element
  npointmax, // number of integration points in element with maximal number of i. p.
  dam_indx, // index stating start of materi_damage in node_dof
  dens_indx, // index stating start of materi_density in node_dof
  diff_indx, // index stating start of materi_diffusion in node_dof
  dis_indx, // index stating start of materi_displacement in node_dof
  ener_indx, // index stating start of materi_strainenergy in node_dof
  epe_indx, // index stating start of materi_strain_elasti in node_dof
  epi_indx, // index stating start of strain_intergranular in node_dof
  epp_indx, // index stating start of materi_strain_plasti in node_dof
  ept_indx, // index stating start of materi_strain_total in node_dof
  maxe_indx, // index stating start of maxwell_e in node_dof
  maxfe_indx, // index stating start of maxwell_fe in node_dof
  maxer_indx, // index stating start of maxwell_e in node_dof
  maxei_indx, // index stating start of maxwell_fe in node_dof
  ei_indx, // index stating start of ei in node_dof
  f_indx, // index stating start of materi_plasti_f in node_dof
  fn_indx, // index stating start of materi_plasti_f_nonlocal in node_dof
  fscal_indx, // index stating start of wave_fscalar in node_dof
  gvel_indx, // index stating start of groundflow_velocity in node_dof
  hisv_indx, // index stating start of materi_history_variables in node_dof
  mstres_indx, // index stating start of materi_maxwell stress in node_dof
  pres_indx, // index stating start of groundflow_pressure in node_dof
  res_indx, // index stating start of residue in node_dof
  kap_indx, // index stating start of materi_plasti_kappa in node_dof
  rho_indx, // index stating start of materi_plasti_rho in node_dof
  rot_indx, // index stating start of beam_rotation in node_dof
  scal_indx, // index stating start of wave_scalar in node_dof
  stres_indx, // index stating start of materi_stress in node_dof
  substeps_indx, // index stating start of materi_plasti_incremental_substeps in node_dof
  svloc_indx, // index stating start of materi_plasti_softvar_local in node_dof
  svnonloc_indx, // index stating start of materi_plasti_softvar_nonlocal in node_dof
  temp_indx, // index stating start of condif_temperature in node_dof
  vel_indx, // index stating start of materi_velocity in node_dof
  veli_indx, // index stating start of materi_velocity in node_dof
  void_indx, // index stating start of materi_void_fraction in node_dof
  work_indx; // index stating start of materi_work in node_dof
extern char
  data_file[MCHAR], // data file name
  data_file_base[MCHAR], // data file name base
  post_calcul_names[DATA_ITEM_SIZE][MCHAR], // names for calcul
  post_calcul_names_without_extension[DATA_ITEM_SIZE][MCHAR], // names for calcul
  initialization_names[DATA_ITEM_SIZE][MCHAR]; // unknown initialization records

  // for map routine
extern long int map_version_from, map_version_to, map_always;

  // for post routines
extern long int post_found, post_node[4], post_node_length, npost_node;
extern double post_point[MDIM], post_point_dof[MUKNWN], 
  post_node_result[DATA_ITEM_SIZE];

  // for calcul routine
extern long int calcul_operat, calcul_matrix, calcul_ecomplex, 
  calcul_vector, calcul_scalar_indx, calcul_mat_indx, calcul_vec_indx;

  // for geometry routine
extern long int geometry_ent[2], *nodes_in_geometry;

  // for parallel computations
extern long int parallel_active;

  // for eigenvalue computations
extern long int eigen_active;

  // for solver
extern long int *solve_global_local, solve_nlocal, solve_options_matrix_group;
extern double *solve_b, *solve_x;

  // for swit
extern long int swit_element_stack, swit_node_stack;
extern char swit_routine_stack[MSTACK][MCHAR]; // routines called

  // added for options_element_dof
extern long int options_element_dof; 	   //in initia.cc
extern double options_nonlocal_softvar;  //in initia.cc

  // routines
void      adjust_geom( long int geometry_entity[], long int geometry_entity_edge[] );
void      area( long int element, long int name, 
            long int gr, long int nnol, long int nodes[],
            double coord[], double new_dof[], double element_lhside[], 
            double element_matrix[], double element_rhside[] );
void      area_element_group( long int version );
void      area_element_group_sequence( );
void      area_node_dataitem( void );
void      array_add( double a[], double b[], double c[], long int n );
double    array_distance( double a[], double b[], double work[], long int n );
double    array_inproduct( double a[], double b[], long int n );
long int  array_member( long int i_list[], long int i, long int n, long int &indx );
void      array_move( long int from[], long int to[], long int n );
void      array_move( double from[], double to[], long int n );
void      array_multiply( double a[], double b[], double c, long int n );
long int  array_normalize( double a[], long int n );
long int  array_null( double dval[], long int n );
void      array_outproduct_2D( double a[], double b[] );
void      array_outproduct_3D( double a[], double b[], double c[] );
void      array_set( long int i[], long int ival, long int n );
void      array_set( double d[], double dval, long int n );
double    array_size( double a[], long int n );
void      array_subtract( double a[], double b[], double c[], long int n );
void      beam_2d( long int element, long int element_group, 
            double coord[], double old_dof[], double new_dof[], 
            double element_lhside[], double element_matrix[],
            double element_rhside[] );    
void      beam_3d( long int element, long int element_group, 
            double coord[], double old_dof[], double new_dof[], 
            double element_lhside[], double element_matrix[],
            double element_rhside[] );    
void      bounda( void );
void      bounda_time_file_apply( long int iboun, double total_time,
            double bounda_time[], long int &ninc );
void      calculate( void );
void      calculate_operat( double unknown_values[], long int inod, 
            double coord[], double dof[],
            double result[], long int &length_result );
long int  check( long int idat, long int task );
long int  check_ndim( long int lower, long int higher, long int task );
long int  check_unknown( const char* str, long int initialization_needs_to_exist,
            long int task );
long int  check_unknown_atleastone( const char* str1, const char* str2, long int task );
long int  check_unknown_minimum( char str[], long int min, long int task );
long int  check_unknowns_are_specified( long int task );
void      crack( void );
void      step_close( long int task, long int ipar, long int npar, long int ipar_i,
            long int ipar_n );
void      step_start( long int task, long int control_solver[],
            double dtime, double time_current );
void      calc_IJlode(double geosigma[9], double &I, double &J, double &lode, bool calcderiv, 
	    double dIdsig[3][3], double dJdsig[3][3], double dlodedsig[3][3]);
void      calc_nonlocal_softvar();
void      cinv( Matrix RealA, Matrix ImagA, 
			 Matrix& RealAinv, Matrix& ImagAinv ); 
void      condif( long int element, long int gr, long int nnol, double h[], 
            double volume, double new_unknowns[], double element_lhside[],
            double element_matrix[], double element_rhside[], 
            double element_residue[] );
void      C_matrix( double young, double poisson,
            double transverse_isotropy[], double C[MDIM][MDIM][MDIM][MDIM], 
            long int task[] );
void      C_matrix_lade( double lade_1, double lade_2, double lade_3, 
            double I1, double sig_dev[], double C[MDIM][MDIM][MDIM][MDIM] );
void      change_geometry( long int task, double dtime, double time_current );
void      contactspring( long int element, long int name, long int element_group, 
            long int nnol, long int nodes[], double coord[], double old_dof[], double new_dof[], 
            double element_lhside[], double element_matrix[],
            double element_rhside[] );    
void      create_element( long int old_element, long int new_element,
            long int el[], long int length_el,
            long int version_from, long int version_to );
void      create_node( long int old_nodes[], long int nnod, 
            long int tmp_node_number, double tmp_node[], double tmp_node_dof[], 
            double tmp_node_start_refined[], double tmp_node_dof_start_refined[],
            long int version_from, long int version_to );
void      damage( long int gr, double new_epe[], double new_sig[], 
            double old_damage, double &new_damage );
void      damage_mazars( double materi_damage_mazars[], double new_epe[], double new_sig[], 
            double old_damage, double &new_damage );
void      data( long int task, double dtime, double time_current );
void      date( void );
long int  db( long int idat, long int index, long int *ival,
            double *dval, long int &length, long int version,
            long int task );
long int  db_active_index( long int idat, long int index, long int version );
void      db_allocate( long int idat, long int index, long int version,
            long int task );
void      db_allocate_class( long int cl, long int index, long int version );
void      db_close( );
void      db_copy( long int idat, long int jdat, long int version );
long int  db_data_length( long int idat );
void      db_data_length_put( long int idat, long int length );
long int  db_data_class( long int idat );
long int  db_data_required( long int idat );
double   *db_dbl( long int idat, long int version );
double   *db_dbl( long int idat, long int index, long int version );
void      db_delete( long int idat, long int version );
void      db_delete_index( long int idat, long int index, long int version );
void      db_error( long int idat, long int index );
long int  db_external( long int idat );
long int  db_fixed_length( long int idat );
void      db_highest_index( long int idat, long int &max, long int version );
void      db_initialize( long int dof_type[], long int dof_label[] );
long int *db_int( long int idat, long int version );
long int *db_int( long int idat, long int index, long int version );
long int  db_len( long int idat, long int index, long int version );
long int  db_max_index( long int idat, long int &max, long int version,
            long int task );
char     *db_name( long int idat );
long int  db_no_index( long int idat );
long int  db_number( char name[] );
long int  db_partialname( long int idat, char *str );
long int  db_partialname_any( const char *str );
long int  db_partialname_any_index( const char *str, long int index );
long int  db_print_only( long int idat );
long int  db_version( long int idat, long int version );
void      db_version_copy( long int version_from, long int version_to );
void      db_version_copy_data( long int idat, long int version_from, long int version_to );
void      db_version_delete( long int version );
void      db_set_int( long int idat, long int version );
void      db_set_dbl( long int idat, long int version );
long int  db_type( long int idat );
void      delete_element( long int element, long int version );
void      delete_geom( double time_current );
void      delete_node( long int inod, long int version );
void      distribute( void );
void      elem( long int element, long int ithread );
void      element_loop( void );
void      element_middle_radius_set( void );
long int  element_residue_norm_set( long int icontrol,
            long int control_refine_locally_unknown,
            long int element, double &element_residue_norm, long int version );
void      element_volume_set( long int name, long int nodes[], long int version, 
            double &element_volume );
void      equal( void );
long int  equations_factor( long int n, double **a, long int *p, long int *f );
void      equations_solve( long int n, double **a, long int *p, double *x,
            double *b, long int *f );
void      error( long int task );
void      exit_tn( long int print_database_type );
void      exit_tn_on_error( void );
void      extrude( void );
void      failure( double time_current );
long int  filter( long int print_filter_index[], long int length_print_filter_index, 
            long int idat, long int index, long int number, long int task );
long int  fit_polynomial( double points[], long int npoint, 
            double coefficients[], long int ncoefficient );
void      force_element_volume_set( long int element, long int nnol, long int nodes[],
            double coord[], double force_element_volume[] );
void      force_time_file_apply( long int iforce, long int table_name, double &load );
void      force_gravity_calculate( double force_gravity[] );
long int  force_time( double time_table[], const char* table_name,
            long int length, double &load );
void      force_factor( long int factor_name, long int iforce, 
            double coord[], double &factor );
void      general( long int element, long int name, long int nnol, long int gr, 
            long int type, long int nodes[], double coord_ip[],
			      double old_dof[], double new_dof[],
            double old_unknowns[], double new_unknowns[], 
            double grad_new_unknowns[], double h[], 
            double d[], double volume, double grad_massflow[],
            double element_rhside[], double element_residue[], 
            double element_lhside[], double element_matrix[] );
void      generate_beam_truss( long int icontrol, long int task );
void      generate_spring( long int icontrol );
void      geometry( long int inod, double co[], long int geometry_entity[],
            long int &found, double &factor, double normal[],
            double &penetration, double projection[],
            long int node_type, long int projection_type, long int version );
void      get_element_matrix_unknowns( long int element,
            long int element_matrix_unknowns[] );
long int  get_group_data( long int idat, long int group, long int element,
            double new_unknowns[], double values[], long int &nvalue, long int task );
double    get_materi_density( long int element, long int element_group, 
            long int nnol, long int nodes[], double new_unknowns[] );
char     *get_new_char( long int n );
double   *get_new_dbl( long int n );
long int *get_new_int( long int n );
int      *get_new_int_short( long int n );
void      groundflow_data( long int element, long int gr, 
            double old_unknowns[], double new_unknowns[], 
            double coord_ip[], double pe[], double &C );
void      groundflow( long int element, long int group, long int nnol,
            double coord_ip[], double h[], double d[], 
            double volume, double old_unknowns[], 
            double new_unknowns[], double grad_new_unknowns[],
            double element_matrix[], double element_rhside[],
            double element_residue[] );
void      groundflow_phreatic_apply( void );
long int  groundflow_phreatic_coord( long int inod, double coord[], double dof[], 
            double &total_pressure, double &static_pressure,
            double &location );
void      group_materi_plasti_boundary_evaluate( long int nodes[], long int nnol,
            long int element_group, long int &plasti_on_boundary );
void      hypoplasticity( long int element, long int gr,
            long int formulation, double old_hisv[], double new_hisv[], 
            double old_unknowns[], double new_unknowns[], 
            double inc_ept[], double old_epi[], double new_epi[], 
            double rotated_old_sig[], double new_sig[],
            double *Chypo, double softvar_nonl, double &softvar_l );
void      hyperelasticity( long int gr, long int element, long int memory,
            double unknowns[], double new_epe[], 
            double sig[], double Chyper[MDIM][MDIM][MDIM][MDIM] );
void      hyper_Cmat( long int gr, long int element, long int memory,
            double unknowns[], double epe[], double Chyper[MDIM][MDIM][MDIM][MDIM] );
void      hyper_law( long int gr, long int element, double unknowns[], double C[], double &W );
long int  hyper_stress( long int gr, long int element, long int memory, double unknowns[],
            double epe[], double stress[] );
void      initialize( void );
void      input();
void      input_check_required();
void      input_convert_to_lower_case( char str[] );
void      input_read_string( long int echo, char str[], double &d, long int &d_is_set );
void      input_runtime( void );
void      input_skip_comment( char str[] );
long int  integration_gauss( long int niso, double iso[], double weight_iso[]);
long int  integration_lobatto( long int niso, double iso[], 
            double weight_iso[] );
void      interpolate_geometry( long int geometry_entity[],
            long int node_numbers[], long int n,
            double test_coord[], double tmp_node[],
            double test_coord_start_refined[], 
            double tmp_node_start_refined[], 
            long int project_type, long int version );
void      interpolation_polynomial( double iso, long int npol, double h_pol[],
            double p_pol[] );
long int  intersect_line_with_line( double line_a0[], double line_a1[], 
            double line_b0[], double line_b1[], double &iso_line_a,
            double &iso_line_b );
long int  intersect_line_with_point( double line0[], double line1[], 
            double point[], double &iso_line );
long int  intersect_line_with_triangle( double line0[], double line1[], 
            double triangle0[], double triangle1[], double triangle2[], 
            double &iso_line, double iso_triangle[] );
void      inverse_calculation( long int ipar, long int npar, long int ipar_i,
            long int ipar_n, long int max_time, long int task );
void      iteration_start( void );
void      itoa( int n, char str[] );
void      locate( void );
char     *long_to_a( long int i, char str[] );
void      macro( void );
void      make_dev(double tnz[9], double dev[9]);
void      map_element( long int element );
void      map_node( long int inod );
void      materi( long int element, long int group, long int nnol, 
            long int npoint, long int nodes[], 
            long int plasti_on_boundary, double coord_ip[],
            double old_coord[], double h[], double new_d[], 
            double new_b[], double volume,
            double old_unknowns[], double new_unknowns[],
            double old_grad_old_unknowns[], double old_grad_new_unknowns[], 
            double new_grad_new_unknowns[], double element_lhside[], 
            double element_matrix[], double element_rhside[], 
            double element_residue[], double tendon_element_rhside[] );
void      materi_diffusion_adjust( void );
void      materi_diffusion_adjust_geom( long int geometry_entity[], 
            long int geometry_entity_edge[] );
void      materi_diffusion_calculate( long int task );
void      materi_diffusion_fill( void );
void      materi_diffusion_smooth( void );
void      materi_diffusion_temperature( void );
void      matrix_ab( double *a, double *b, double *c, long int n, long int m,
            long int k );
void      matrix_abat( double a[], double b[], double c[], double work[],
            long int n );
void      matrix_abt( double *a, double *b, double *c, long int n, long int m,
            long int k );
void      matrix_atb( double *a, double *b, double *c, long int n, long int m,
            long int k );
void      matrix_atba( double a[], double b[], double c[], 
            double work[], long int n, long int m );
void      matrix_a4b( double a[MDIM][MDIM][MDIM][MDIM], double b[], 
            double c[] );
void      matrix4_ab( double a[], double b[], double c[3][3][3][3] );
void      matrix_a_contr_b( double a[], double b[], double &c );
void      matrix_ab4( double a[], double b[3][3][3][3], double c[] );
double    matrix_determinant( double a[], long int n );
void      matrix_eigenvalues( double a[], double eigenvalues[] );
void      matrix_insert( double *a, long int n, long int m,
            double *b, long int k, long int l, long int p );
void      matrix_invariants( double matrix[], double inv[] );
long int  matrix_inverse( double mat[], double inv_mat[], double &det, 
            long int n );
void      matrix_inverse_general(double *matr, double *inv, int P);
void      matrix_jacobi( double *a, long int n, double d[], 
            double *v, long int *nrot);
void      matrix_symmetric( double a[], long int n );
void      maxwell( long int type, long int element, 
            long int gr, long int nnol, double volume,
            double new_unknowns[], 
            double old_dof[], double new_dof[],
            double h[], double d[], double element_lhside[],
            double element_matrix[], double element_matrix_second[],
            double element_rhside[] );
void      maxwell_scatter( void );
long int  membrane_apply( long int element, long int gr, double memmat[3][3], 
            double inc_ept[], double inc_epe[], double new_ept[],
            double new_epe[], double new_sig[] );
void      merge( void );
void      mesh_add( long int version_from, long int version_to );
void      mesh_delete_small( long int version );
void      mesh_has_changed( long int version );
void      mesh_split( long int version );
void      new_mesh( void );
void      new_mesh_version( long int version, double delta );
void      nod_nod( long int version );
void      nonlocal_set( );
void      nonlocal_apply( );
void      ordered_list_apply( long int inod, long int ordered_nodes[], 
            long int new_max_node, double coord_new[], 
            double eps_coord, long int &equal, long int task );
void      parallel_calcul_node( void );
void      parallel_contact( void );
void      parallel_element_loop( void );
void      parallel_geometry( void );
void      parallel_map_element( void );
void      parallel_map_node( void );
void      parallel_new_dof_before( void );
void      parallel_new_dof_diagonal( void );
void      parallel_post_point( void );
void      parallel_post_node( void );
void      parallel_solve_iterative_bicg_element( void );
void      parallel_solve_iterative_bicgs_element( void );
void      parallel_sys_initialize( void );
void      parallel_sys_lock( void );
void      parallel_sys_next_of_loop( long int next_of_loop[], 
            long int max_loop, long int &nloop, long int &ithread );
void      parallel_sys_routine( void (*routine)(void) );
void      parallel_sys_unlock( void );
void      plasti_incr( double rotated_old_sig[], double inc_ept[], double old_epp[], 
	    double new_ept[], double C[3][3][3][3], double old_hisv[], long int plasti_type, 
	    double plasti_dt[], long int lenght_pl, long int gr, long int element,
	    double softvar_nonl, double &softvar_l,
            double new_sigv[], double inc_epp[], double new_hisv[], double &new_f, 
	    double &new_substeps, 
	    double Cep_cons[3][3][3][3]);
void      plasti_rule( long int element, long int group, 
            long int plasti_on_boundary, double user_data[],
            double new_unknowns[], double new_grad_new_unknowns[],
            double old_hisv[], double new_hisv[], 
            double old_epp[], double inc_epp[], double inc_ept[],
            long int task, long int &plasti_type, 
            double sig[], double &f, double &new_f, double dir[] );
long int  point_el( double point[], double coord[], double weight[],
            long int name, long int nnol );
void      pol( long int element, long int element_group,
            long int name, long int nnol, double old_coord[], 
            double new_coord[], long int &npoint, double h[], 
            double old_d[], double new_d[], double new_b[], 
            double volume[] );
void      post( long int task );
void      post_global( void );
void      post_integrate( void );
void      post_node_rhside_fixed_free( void );
void      pri( const char *s );
void      pri( const char *s, const char *st );
void      pri( const char *s, int i );
void      pri( const char *s, long int i);
void      pri( const char *s, double d );
void      pri( const char *s, int *i, int n );
void      pri( const char *s, long int *i, long int n );
void      pri( const char *s, double *d, long int n );
void      pri( const char *s, long int *iar, long int n, long int m );
void      pri( const char *s, double *d, long int n, long int m );
void      print_cmd( void );
void      print_data_versus_data( long int ival[], long int length );
void      print_database( long int task, long int version, long int idat );
void      print_dx( long int final_call );
void      print_element( long int data_item_name );
void      print_gid( long int task );
long int  print_gid_5( long int task );
long int  print_gid_6( long int task );
void      print_giddata( long int task );
void      print_gmv( long int icontrol, long int ival[] );
void      print_history( long int ival[], long int nval );
void      print_plotmtv( long int icontrol, long int ival[] );
void      print_matlab( void );
void      print_restart( long int icontrol );
void      print_tecplot( long int ival[] );
void      print_unknowns( void );
void      print_unknownsrhside( void );
void      print_vtk( long int icontrol );
long int  project_point_exactly_on_line( double coord[], double coord0[], 
            double coord1[], double weight[] );
long int  project_point_exactly_on_quad( double coord[], double coord0[], 
            double coord1[], double coord2[], double coord3[], double weight[] );
long int  project_point_exactly_on_triangle( double coord[], double coord0[], 
            double coord1[], double coord2[], double weight[] );
long int  project_point_on_line( double coord[], double coord0[], 
            double coord1[], double weight[] );
long int  project_point_on_triangle( double coord[], double coord0[], 
            double coord1[], double coord2[], double weight[] );
void      range_expand( long int ival[], long int integer_range[],
            long int &length, long int &range_length );
void      range_scan( long int echo, double d, long int d_is_set,
            long int ival[], long int &length );
void      range_test_expand( long int i );
void      range_test_scan( long int i );
void      refine_locally( void );
void      refine_globally( long int control_refine_globally[4], 
            long int length_control_refine_globally, 
            long int use_control_refine_globally_geometry, 
            long int control_refine_globally_geometry[2],
            long int project_type, long int version );
void      remesh( long int version );
void      renumbering( long int version, long int fill_old_numbers, long int lowest_element, 
            long int lowest_node, long int old_node_numbers[], 
            long int old_element_numbers[] );
void      renumbering_check( long int idat );
long int  repeat( long int &start_control );
void      restart( void );
double    scalar_dabs( double a );
double    scalar_dmax( double a, double b );
double    scalar_dmin( double a, double b );
long int  scalar_imax( long int a, long int b );
double    scalar_power( double a, double b );
double    scalar_ran_normal( int &idum );
double    scalar_ran_uniform( int &idum );
double    scalar_sign( double a );
double    scalar_square( double a );
void      set_deften_etc( long int element, long int gr, long int nnol, 
            double h[], double old_coord[],
            double old_unknowns[], double new_unknowns[], 
            double old_grad_old_unknowns[], double old_grad_new_unknowns[], 
            double old_deften[], double new_deften[], 
            double old_ept[], double inc_ept[], double new_ept[],
            double old_rot[], double inc_rot[], double new_rot[] );
void      set_deften_u_rot( double deften[], double u[], double rot[] );
void      set_environment( void );
void      set_stress( long int element, long int gr, 
            long int plasti_on_boundary, double coord_ip[],
            double old_unknowns[], double new_unknowns[], 
			      double old_grad_old_unknowns[], double new_grad_new_unknowns[],
            double rotated_old_sig[], double new_sig[],
            double rotated_old_msig[], double new_msig[], 
            double inc_ept[], double new_ept[], 
			      double old_epe[], double inc_epe[], 
            double old_epp[], double inc_epp[], 
            double old_rho[], double new_rho[],
            double old_epi[], double new_epi[],
            double old_hisv[], double new_hisv[], 
            double old_damage, double &new_damage,
            double old_kappa, double &new_kappa,
            double &new_f, double &new_substeps, double old_deften[], double new_deften[],
            double inc_rot[], double ddsdde[],
            double &viscosity, double &viscosity_heat_generation, 
	    double &softvar_nonl, double &softvar_l );
long int  set_swit( long int element, long int inod, const char *routine
 );
void      slide( void );
void      solve( long int task );
long int  solve_direct_symmetric_factor( long int n, double *a, 
            long int band, long int *p, long int *f );
long int  solve_direct_symmetric_substi( long int n, double *a, 
            long int band, long int *p, long int *f );
long int  solve_iterative_bicg( void );
void      solve_iterative_bicg_sys( double *Ad1, double *Ad2, double *p_tmp, double *residue );
void      solve_iterative_bicg_element( long int element, long int ithread );
void      sort( double val[], double vec[] );
void      spring( long int element, long int name, long int element_group, 
            double coord[], double old_dof[], double new_dof[], 
            double element_lhside[], double element_matrix[],
            double element_rhside[] );    
long int  stress_indx( long int idim, long int jdim );
void      stress_umat( long int element, long int gr, long int formulation,
            long int nuser_data, double user_data[], double coord_ip[],
            double old_hisv[], double new_hisv[], 
            double old_unknowns[], double new_unknowns[], 
            double inc_ept[], double new_ept[], 
            double rotated_old_sig[], double new_sig[],
            double old_deften[], double new_deften[],
            double inc_rot[], double ddsdde[] );
void      string_convert_to_lower_case( char str[] );
long int  string_isdouble( char str[] );
long int  string_isinteger( char str[] );
void      string_replace( char str[], char from, char to );
void      string_reverse( char str[] );
void      string_shorten( char str[], long int length );
long int  table_xy( double table[], const char* table_name, 
            long int length, double x, double &y );
long int  table_xyz( double table[], long int number[2], double coord[], double &z );
void      tendon_distribute();
void      tendons( long int element, long int gr, long int nnol, 
            long int npoint, double volume,
            double new_d[], double old_unknowns[], double new_unknowns[],
            double new_rot[], double inc_ept[], double tendon_element_rhside[],
            double ddsdde_tendon[] );
double    tetrahedron_volume( double c0[], double c1[], 
            double c2[], double c3[] );
void      top( void );
double    triangle_area( double c0[], double c1[], double c2[] );
void      truss( long int element, long int element_group, 
            double coord[], double old_dof[], double new_dof[], 
            double element_lhside[], double element_matrix[],
            double element_rhside[] );    
void      unknown_freeze( void );
void      unknown_reset( void );
void      user_bounda_time( long int index, double time_current, 
            double &load );
void      user_change_dataitem_time( long int index, 
            double time_current, double &val );
void      user_change_geometry_time( long int index, 
            double time_current, double &val );
void      user_plasti( long int task, double user_data[],
            double new_unknowns[], double old_hisv[],
            double new_hisv[], double sig[], double &f );
void      user_sigma( double user_data[], double new_unknowns[], 
            double inc_epe[], double old_hisv[], double new_hisv[], 
            double old_sig[], double new_sig[], 
            double Cuser[MDIM][MDIM][MDIM][MDIM] );
void      user_viscosity( double user_data[], double new_unknowns[], 
            double &visc );
void      visco_elasticity( long int element, long int group, long int formulation,
            double new_unknowns[], double inc_epe[], double rotated_old_msig[],
            double new_sig[], double new_msig[], double memmat[MDIM][MDIM] );
void      visco_elastiticity_nonlinear( long int m, long int gr, long int element, 
            double new_unknowns[], double new_sig[], double new_msig[], 
            double em, double tm );
void      viscous_stress( long int element, long int gr, double user_data[],
            double unknowns[], double grad_unknowns[],
            double sig[], double &viscosity,
            double &viscosity_heat_generation );
void      volume_factor( long int element_group, double coord[], double &volfac );
void      wave( long int element, long int gr, long int nnol,
            double h[], double d[], double volume, double new_unknowns[], 
            double grad_new_unknowns[], double element_lhside[], double element_matrix[],
            double element_rhside[], double element_residue[] );
