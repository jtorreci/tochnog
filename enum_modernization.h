/*
 * enum_modernization.h - Modern enum classes for Tochnog
 * Provides type-safe enum classes as an alternative to traditional enums
 */

#ifndef ENUM_MODERNIZATION_H
#define ENUM_MODERNIZATION_H

#include "tochnog.h"

// Modern enum class for versions - provides type safety
enum class VersionType : int {
    Normal = VERSION_NORMAL,      // time=t
    Start = VERSION_START,        // time=start_time 
    New = VERSION_NEW,            // time=t+dt
    Temporary = VERSION_TMP,      // trash version used by several routines
    NewMeshTmp = VERSION_NEW_MESH_TMP, // trash version used by new_mesh
    NewMeshGenerated = VERSION_NEW_MESH_GENERATED, // generated mesh in new_mesh
    Print = VERSION_PRINT,        // mesh for printing
    Macro = VERSION_MACRO,        // mesh for control_macro
    Extrude = VERSION_EXTRUDE,    // mesh for extrude
    Max = MVERSION                // maximum number of versions
};

// Modern enum class for data types - provides type safety
enum class DataType : int {
    // Using the same values as the original enum for compatibility
    MinusOne = MINUS_ONE,
    Above = ABOVE,
    Absol = ABSOL,
    Add = ADD,
    AddAlways = ADD_ALWAYS,
    All = ALL,
    Any = ANY,
    Area = AREA,
    AreaElementGroup = AREA_ELEMENT_GROUP,
    AreaElementGroupMethod = AREA_ELEMENT_GROUP_METHOD,
    AreaElementGroupSequence = AREA_ELEMENT_GROUP_SEQUENCE,
    AreaElementGroupSequenceElement = AREA_ELEMENT_GROUP_SEQUENCE_ELEMENT,
    AreaElementGroupSequenceElementGroup = AREA_ELEMENT_GROUP_SEQUENCE_ELEMENTGROUP,
    AreaElementGroupSequenceGeometry = AREA_ELEMENT_GROUP_SEQUENCE_GEOMETRY,
    AreaElementGroupSequenceMethod = AREA_ELEMENT_GROUP_SEQUENCE_METHOD,
    AreaElementGroupSequenceTime = AREA_ELEMENT_GROUP_SEQUENCE_TIME,
    AreaNodeDataitem = AREA_NODE_DATAITEM,
    AreaNodeDataitemDouble = AREA_NODE_DATAITEM_DOUBLE,
    AreaNodeDataitemInteger = AREA_NODE_DATAITEM_INTEGER,
    Asm = ASM,
    Average = AVERAGE,
    Bar = BAR,
    Bar2 = BAR2,
    Bar3 = BAR3,
    Bar4 = BAR4,
    Beam = BEAM,
    BeamRotation = BEAM_ROTATION,
    Below = BELOW,
    Bcgs = BCGS,
    Bicg = BICG,
    Bjacobi = BJACOBI,
    Bounda = BOUNDA,
    // ... (would continue for all values)
};

// Conversion functions to maintain compatibility
inline int to_legacy_version(VersionType v) { return static_cast<int>(v); }
inline VersionType from_legacy_version(int v) { return static_cast<VersionType>(v); }

inline int to_legacy_data_type(DataType dt) { return static_cast<int>(dt); }
inline DataType from_legacy_data_type(int dt) { return static_cast<DataType>(dt); }

// Helper function to check version validity
inline bool is_valid_version(VersionType version) {
    return (static_cast<int>(version) >= 0 && static_cast<int>(version) < MVERSION);
}

#endif // ENUM_MODERNIZATION_H