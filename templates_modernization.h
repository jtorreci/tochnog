/*
 * templates_modernization.h - Modern C++ templates and generics
 * Replaces macros and duplicated code with generic templates
 */

#ifndef TEMPLATES_MODERNIZATION_H
#define TEMPLATES_MODERNIZATION_H

#include "tochnog.h"
#include <vector>
#include <algorithm>
#include <functional>
#include <type_traits>
#include <memory>

// Generic array operations to replace duplicated code for different types
template<typename T>
class GenericArrayOps {
public:
    // Find maximum value in array
    static T max(const T* arr, long int length) {
        if (!arr || length <= 0) return T{};
        T result = arr[0];
        for (long int i = 1; i < length; ++i) {
            if (arr[i] > result) {
                result = arr[i];
            }
        }
        return result;
    }
    
    // Find minimum value in array
    static T min(const T* arr, long int length) {
        if (!arr || length <= 0) return T{};
        T result = arr[0];
        for (long int i = 1; i < length; ++i) {
            if (arr[i] < result) {
                result = arr[i];
            }
        }
        return result;
    }
    
    // Set all values in array
    static void set(T* arr, const T& value, long int length) {
        if (!arr || length <= 0) return;
        for (long int i = 0; i < length; ++i) {
            arr[i] = value;
        }
    }
    
    // Move/copy array elements
    static void move(const T* src, T* dest, long int length) {
        if (!src || !dest || length <= 0) return;
        for (long int i = 0; i < length; ++i) {
            dest[i] = src[i];
        }
    }
    
    // Add two arrays element-wise
    static void add(const T* a, const T* b, T* result, long int length) {
        if (!a || !b || !result || length <= 0) return;
        for (long int i = 0; i < length; ++i) {
            result[i] = a[i] + b[i];
        }
    }
};

// Type-traits based utilities for numeric types
template<typename T>
struct is_arithmetic {
    static constexpr bool value = std::is_arithmetic<T>::value;
};

// Generic version of array operations that were previously duplicated for int/double
template<typename T>
inline T safe_array_max(const T* arr, long int length) {
    static_assert(is_arithmetic<T>::value, "T must be an arithmetic type");
    return GenericArrayOps<T>::max(arr, length);
}

template<typename T>
inline T safe_array_min(const T* arr, long int length) {
    static_assert(is_arithmetic<T>::value, "T must be an arithmetic type");
    return GenericArrayOps<T>::min(arr, length);
}

template<typename T>
inline void safe_array_set(T* arr, const T& value, long int length) {
    static_assert(is_arithmetic<T>::value, "T must be an arithmetic type");
    GenericArrayOps<T>::set(arr, value, length);
}

template<typename T>
inline void safe_array_move(const T* src, T* dest, long int length) {
    static_assert(is_arithmetic<T>::value, "T must be an arithmetic type");
    GenericArrayOps<T>::move(src, dest, length);
}

// Generic function to replace type-specific functions where possible
template<typename T>
inline T clamp_value(T value, T min_val, T max_val) {
    static_assert(is_arithmetic<T>::value, "T must be an arithmetic type");
    if (value < min_val) return min_val;
    if (value > max_val) return max_val;
    return value;
}

// Generic smart pointer wrapper for the get_new functions
template<typename T>
class TochnogArray {
private:
    std::unique_ptr<T[]> data;
    long int size_val;
    
public:
    explicit TochnogArray(long int n) : size_val(n > 0 ? n : 1) {
        data = std::make_unique<T[]>(size_val);
    }
    
    T& operator[](long int index) {
        if (index < 0 || index >= size_val) {
            // In debug mode, this would throw an exception
            // For compatibility, we'll just return first element
            #ifndef NDEBUG
            throw std::out_of_range("Index out of bounds");
            #endif
            return data[0];
        }
        return data[index];
    }
    
    const T& operator[](long int index) const {
        if (index < 0 || index >= size_val) {
            #ifndef NDEBUG
            throw std::out_of_range("Index out of bounds");
            #endif
            return data[0];
        }
        return data[index];
    }
    
    T* get() { return data.get(); }
    const T* get() const { return data.get(); }
    long int size() const { return size_val; }
    
    void fill(const T& value) {
        for (long int i = 0; i < size_val; ++i) {
            data[i] = value;
        }
    }
};

// Template alias for common array types
template<typename T>
using SafeArray = TochnogArray<T>;

// Generic swap function (type-safe replacement for macro)
template<typename T>
inline void safe_swap(T& a, T& b) {
    T temp = std::move(a);
    a = std::move(b);
    b = std::move(temp);
}

// Generic comparison function
template<typename T>
inline bool safe_compare(const T& a, const T& b, const std::function<bool(const T&, const T&)>& comparator) {
    return comparator(a, b);
}

#endif // TEMPLATES_MODERNIZATION_H