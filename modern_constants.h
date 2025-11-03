/*
 * modern_constants.h - Modern C++ constants replacing C-style macros
 * Provides type-safe, scoped alternatives to #define constants
 */

#ifndef MODERN_CONSTANTS_H
#define MODERN_CONSTANTS_H

#include "tochnog.h"

// Modern alternatives to macro constants using constexpr
namespace TochnogConstants {
    // Size constants - replacing #define macros
    constexpr long int MCHAR = 100;  ///< Maximum length of names
    constexpr long int MDIM = 3;     ///< Maximum number of space dimensions  
    constexpr long int MNOL = 27;    ///< Maximum number of nodes in an element, set if to 64 for HEX64
    constexpr long int MSTRAIN = 6;  ///< Maximum number of strain components
    constexpr long int MTENDON = 10; ///< Maximum number of tendons in an element
    constexpr long int MTYPE = 10;   ///< Maximum number of types
    constexpr long int MCALCUL = 20; ///< Maximum length of calcul records
    constexpr long int MRANGE = 500000; ///< Maximum range length
    constexpr long int MMAXWELL = 50;   ///< Maximum number of maxwell chains
    constexpr long int MTHREAD = 64;    ///< Maximum number of threads
    constexpr long int MAXIMUM_NODE = 64; ///< Always 64
    constexpr long int DATA_ITEM_SIZE = 190; ///< Maximum length of (almost all) records
    constexpr long int NONLOCAL_ITEM_SIZE = 160; ///< Maximum number of integration points for nonlocal calculation
    constexpr long int TN_PRECISION = 12; ///< Precision in writing output file
    constexpr long int MPOINT = MNOL; ///< Maximum number of integration points in an element, always MNOL
    constexpr long int MUKNWN = DATA_ITEM_SIZE; ///< Maximum number of unknowns, always DATA_ITEM_SIZE
    constexpr long int MPUKNWN = MUKNWN; ///< Maximum number of primary unknowns, always MUKNWN
    constexpr long int MPRINC = 10; ///< Maximum number of principal unknowns
    constexpr long int MBOUNDA = 1000; ///< Maximum length bounda_unknown and bounda_force
    constexpr long int MSTACK = 10; ///< Maximum number of routines in routine stack
    constexpr long int TN_EXIT_STATUS = 1; ///< Exit status code
}

// For backward compatibility, keep the original macros
// In new code, encourage use of TochnogConstants instead

// Modern alternative to macro-based min/max functions using templates
template<typename T>
constexpr T modern_min(const T& a, const T& b) {
    return (a < b) ? a : b;
}

template<typename T>
constexpr T modern_max(const T& a, const T& b) {
    return (a > b) ? a : b;
}

// Safe version of macro with type checking
template<typename T>
constexpr T safe_abs(const T& value) {
    return value < 0 ? -value : value;
}

#endif // MODERN_CONSTANTS_H