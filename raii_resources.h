/*
 * raii_resources.h - RAII resource management for Tochnog
 * Provides safe wrappers for resource management using RAII principle
 */

#ifndef TOCHNOG_RAI_RESOURCES_H
#define TOCHNOG_RAI_RESOURCES_H

#include "tochnog.h"
#include <memory>
#include <vector>

// RAII wrapper for the get_new_* functions
// Instead of returning raw pointers, these could return safer containers

// Safe array wrapper that follows RAII
template<typename T>
class SafeArray {
private:
    std::vector<T> data;
    
public:
    explicit SafeArray(long int n) {
        if (n <= 0) n = 1;
        data.resize(n);
    }
    
    // Access operators
    T& operator[](long int index) { return data[index]; }
    const T& operator[](long int index) const { return data[index]; }
    
    // Get raw pointer when needed for compatibility (dangerous but sometimes necessary)
    T* get() { return data.data(); }
    const T* get() const { return data.data(); }
    
    long int size() const { return static_cast<long int>(data.size()); }
    
    // Fill with default value
    void fill(const T& value) { std::fill(data.begin(), data.end(), value); }
};

// RAII wrappers for the get_new_* functions
inline std::unique_ptr<char[]> make_safe_char_array(long int n) {
    if (n <= 0) n = 1;
    try {
        return std::make_unique<char[]>(n);
    } catch (const std::bad_alloc&) {
        pri("Error: cannot allocate enough memory.");
        exit(TN_EXIT_STATUS);
        return nullptr; // Unreachable but makes compiler happy
    }
}

inline std::unique_ptr<double[]> make_safe_dbl_array(long int n) {
    if (n <= 0) n = 1;
    try {
        return std::make_unique<double[]>(n);
    } catch (const std::bad_alloc&) {
        pri("Error: cannot allocate enough memory.");
        exit(TN_EXIT_STATUS);
        return nullptr; // Unreachable but makes compiler happy
    }
}

inline std::unique_ptr<long int[]> make_safe_int_array(long int n) {
    if (n <= 0) n = 1;
    try {
        return std::make_unique<long int[]>(n);
    } catch (const std::bad_alloc&) {
        pri("Error: cannot allocate enough memory.");
        exit(TN_EXIT_STATUS);
        return nullptr; // Unreachable but makes compiler happy
    }
}

// For compatibility, we keep the old functions but mark them for deprecation
// In a future version, we might replace their implementation with RAII versions

#endif // TOCHNOG_RAI_RESOURCES_H