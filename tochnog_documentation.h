/*
 * tochnog_documentation.h - Comprehensive function documentation for Tochnog
 * English documentation for complex function signatures to improve maintainability
 */

#ifndef TOCHNOG_DOCUMENTATION_H
#define TOCHNOG_DOCUMENTATION_H

#include "tochnog.h"

/**
 * @brief Main database access function for Tochnog
 * 
 * This function serves as the primary interface to the Tochnog database system,
 * providing unified access for getting, putting, and managing data across various
 * categories like nodes, elements, degrees of freedom, etc.
 * 
 * @param idat        Data item identifier (e.g., NODE, ELEMENT, DOF_LABEL)
 *                    Use negative values to indicate data item types for reading
 *                    Use positive values for different data categories
 * @param index       Index within the data category (record number)
 *                    For data without index, use 0. Range typically 0 to max elements/nodes
 * @param int_arr     Integer data array for input/output
 *                    Input: When action is PUT, provides integer values to store
 *                    Output: When action is GET, receives integer values retrieved
 *                    Size must match 'length' parameter
 * @param dbl_arr     Double data array for input/output  
 *                    Input: When action is PUT, provides double values to store
 *                    Output: When action is GET, receives double values retrieved
 *                    Size must match 'length' parameter
 * @param length      Reference to array length
 *                    Input: When action is PUT, specifies length of data to store
 *                    Output: When action is GET, updated with actual length retrieved
 *                    Must be properly initialized before call
 * @param version     Temporal version of data to access
 *                    Common values: VERSION_NORMAL (current), VERSION_START (initial)
 *                                  VERSION_NEW (next time step)
 *                    Allows handling different temporal states of the mesh/model
 * @param action      Action to perform on the data
 *                    Possible values: GET (retrieve existing data)
 *                                   PUT (store new data, creates record if needed)  
 *                                   GET_IF_EXISTS (retrieve if exists, no error if missing)
 *                                   GET_AND_CHECK (retrieve and verify length matches)
 * 
 * @return long int   Status code indicating success or failure
 *                    1: Operation successful
 *                    0: Operation failed (typically for GET_IF_EXISTS when data absent)
 * 
 * @throws DatabaseException When compiled with USE_EXCEPTIONS and errors occur
 * 
 * @note This is the most frequently called function in Tochnog, used thousands of times
 * @note Handles automatic space allocation if action is PUT and record doesn't exist
 * @warning Ensure arrays have sufficient size to prevent buffer overflows
 * @warning Index bounds must be valid to prevent database corruption
 * 
 * @par Example usage for reading data:
 * @code
 *     long int int_data[10]; 
 *     double dbl_data[20];
 *     long int length = 0;
 *     
 *     // Retrieve DOF information for node 100
 *     long int result = db(NODE_DOF, 100, int_data, dbl_data, length, 
 *                         VERSION_NORMAL, GET);
 *     if (result == 1) {
 *         // Data retrieved successfully, length now contains actual data length
 *     }
 * @endcode
 * 
 * @par Example usage for writing data:
 * @code  
 *     long int new_int_data[] = {1, 2, 3};
 *     double new_dbl_data[] = {1.0, 2.0, 3.0};
 *     long int length = 3;  // Length of data to write
 *     
 *     // Store GROUP information for element 50
 *     long int result = db(ELEMENT_GROUP, 50, new_int_data, new_dbl_data, 
 *                         length, VERSION_NORMAL, PUT);
 * @endcode
 * 
 * @sa db_int - Direct access to integer database arrays
 * @sa db_dbl - Direct access to double database arrays  
 * @sa db_active_index - Check if a specific index is active
 * @sa db_max_index - Get/set maximum allocated index
 */

/**
 * @brief Safe memory allocation for character arrays
 * 
 * Wrapper function that allocates character arrays with proper error handling.
 * Replaces raw 'new char[n]' allocations with a consistent error-handling pattern.
 * 
 * @param n Number of characters to allocate (will be adjusted to minimum of 1 if <= 0)
 * @return char* Pointer to newly allocated character array, or terminates if allocation fails
 * 
 * @note Always returns valid pointer or exits program (no null returns)
 * @warning Caller is responsible for ensuring proper deallocation if needed elsewhere
 * @sa get_new_dbl, get_new_int for equivalent functions for other types
 */
char *get_new_char( long int n );

/**
 * @brief Safe memory allocation for double arrays
 * 
 * Wrapper function that allocates double arrays with proper error handling.
 * Replaces raw 'new double[n]' allocations with consistent error handling.
 * 
 * @param n Number of doubles to allocate (will be adjusted to minimum of 1 if <= 0)  
 * @return double* Pointer to newly allocated double array, or terminates if allocation fails
 * 
 * @note Always returns valid pointer or exits program (no null returns)
 * @warning Caller is responsible for ensuring proper deallocation if needed elsewhere
 * @sa get_new_char, get_new_int for equivalent functions for other types
 */
double *get_new_dbl( long int n );

/**
 * @brief Safe memory allocation for integer arrays
 * 
 * Wrapper function that allocates long integer arrays with proper error handling.
 * Replaces raw 'new long int[n]' allocations with consistent error handling.
 * 
 * @param n Number of integers to allocate (will be adjusted to minimum of 1 if <= 0)
 * @return long int* Pointer to newly allocated integer array, or terminates if allocation fails
 * 
 * @note Always returns valid pointer or exits program (no null returns) 
 * @warning Caller is responsible for ensuring proper deallocation if needed elsewhere
 * @sa get_new_char, get_new_dbl for equivalent functions for other types
 */
long int *get_new_int( long int n );

/**
 * @brief Check if a specific database index is active
 * 
 * Determines whether a particular index within a database category is currently in use.
 * Used to verify existence of data before accessing it.
 * 
 * @param idat    Data item identifier (NODE, ELEMENT, etc.)
 * @param index   Index to check within the category
 * @param version Version of data to check (VERSION_NORMAL, VERSION_NEW, etc.)
 * 
 * @return bool True if the index is active (contains data), false otherwise
 * 
 * @note Used extensively in loops to avoid invalid data access
 * @sa db - Main database access function that performs similar checks internally
 */
bool db_active_index( long int idat, long int index, long int version );

/**
 * @brief Get or set the maximum allocated index for a data category
 * 
 * Manages the highest index that has been allocated for a particular type of data.
 * Used internally by allocation functions and can be used to traverse all allocated entries.
 * 
 * @param idat    Data item identifier
 * @param max     Reference to store/receive the maximum index value  
 * @param version Version of data to query
 * @param task    Operation to perform (GET to retrieve, PUT to set)
 * 
 * @sa db_allocate - Function that increases maximum index when allocating new space
 */
void db_max_index( long int idat, long int &max, long int version, long int task );

#endif // TOCHNOG_DOCUMENTATION_H