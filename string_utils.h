/*
 * string_utils.h - Safe string utilities for Tochnog
 * Replaces C-style string operations with std::string equivalents
 */

#ifndef TOCHNOG_STRING_UTILS_H
#define TOCHNOG_STRING_UTILS_H

#include "tochnog.h"
#include <string>
#include <vector>
#include <memory>
#include <cstdarg>
#include <algorithm>

// Safe string operations to replace strcpy/strcat/strlen etc.
inline std::string safe_string_copy(const std::string& source, size_t max_len = MCHAR-1) {
    return source.substr(0, max_len);
}

inline std::string safe_string_concat(const std::string& str1, const std::string& str2, size_t max_len = MCHAR-1) {
    std::string result = str1 + str2;
    return result.substr(0, max_len);
}

inline std::string safe_string_n_concat(const std::string& str1, const std::string& str2, size_t max_len = MCHAR-1) {
    std::string result = str1;
    if (result.length() + str2.length() <= max_len) {
        result += str2;
    } else {
        // Ensure we don't exceed max length
        size_t remaining = max_len - result.length();
        if (remaining > 0) {
            result += str2.substr(0, remaining);
        }
    }
    return result;
}

// String formatting utility (safer alternative to sprintf)
template<typename... Args>
std::string string_format(const std::string& format, Args... args) {
    // This is a simplified version - in practice would use more sophisticated formatting
    // For the actual implementation, consider using std::format (C++20) or fmt library
    int size_s = std::snprintf(nullptr, 0, format.c_str(), args...) + 1; // Extra space for '\0'
    auto size = static_cast<size_t>(size_s);
    std::unique_ptr<char[]> buf(new char[size]);
    std::snprintf(buf.get(), size, format.c_str(), args...);
    return std::string(buf.get(), buf.get() + size - 1); // We don't want the '\0' inside
}

// Safe version of long_to_a using std::string
inline std::string long_to_string_safe(long int n) {
    return std::to_string(n);
}

// Convert C-style char array to std::string with bounds checking
inline std::string char_array_to_string(const char* c_str, size_t max_len = MCHAR) {
    if (!c_str) return "";
    std::string result(c_str);
    if (result.length() > max_len) {
        return result.substr(0, max_len);
    }
    return result;
}

// Safe string operations for file paths and names
class FileNameBuilder {
private:
    std::string base_name;
    std::string extension;
    
public:
    FileNameBuilder(const std::string& name = "") : base_name(name) {}
    
    FileNameBuilder& set_base(const std::string& name) {
        base_name = name.length() <= MCHAR - 20 ? name : name.substr(0, MCHAR - 20); // Reserve space for extension
        return *this;
    }
    
    FileNameBuilder& set_extension(const std::string& ext) {
        extension = ext.length() <= 20 ? ext : ext.substr(0, 20); // Limit extension size
        return *this;
    }
    
    std::string build() const {
        std::string result = base_name;
        if (!extension.empty() && extension[0] != '.') {
            result += "." + extension;
        } else if (!extension.empty()) {
            result += extension;
        }
        return result.length() <= MCHAR ? result : result.substr(0, MCHAR);
    }
};

// Utility function to safely copy std::string to char array
inline void safe_string_copy_to_array(const std::string& src, char* dest, size_t dest_size) {
    if (!dest) return;
    size_t copy_len = src.length() < dest_size - 1 ? src.length() : dest_size - 1;
    src.copy(dest, copy_len);
    dest[copy_len] = '\0';  // Ensure null termination
}

#endif // TOCHNOG_STRING_UTILS_H