/*
 * tochnog_exceptions.h - Exception hierarchy for Tochnog
 * Provides modern C++ exception handling instead of exit() calls
 */

#ifndef TOCHNOG_EXCEPTIONS_H
#define TOCHNOG_EXCEPTIONS_H

#include <exception>
#include <string>
#include <iostream>

// Base exception class for Tochnog
class TochnogException : public std::exception {
protected:
    std::string msg;
    
public:
    explicit TochnogException(const std::string& message) : msg(message) {}
    const char* what() const noexcept override { return msg.c_str(); }
};

// Database-related exceptions
class DatabaseException : public TochnogException {
public:
    explicit DatabaseException(const std::string& message) 
        : TochnogException("Database Error: " + message) {}
};

class DatabaseIndexException : public DatabaseException {
public:
    explicit DatabaseIndexException(long int index, long int idat) 
        : DatabaseException("Invalid index " + std::to_string(index) + " for data " + std::to_string(idat)) {}
};

class DatabaseVersionException : public DatabaseException {
public:
    explicit DatabaseVersionException(long int version) 
        : DatabaseException("Invalid version " + std::to_string(version)) {}
};

// Memory-related exceptions
class MemoryException : public TochnogException {
public:
    explicit MemoryException(const std::string& message) 
        : TochnogException("Memory Error: " + message) {}
};

class OutOfMemoryException : public MemoryException {
public:
    explicit OutOfMemoryException() : MemoryException("Cannot allocate enough memory") {}
};

// Solver-related exceptions
class SolverException : public TochnogException {
public:
    explicit SolverException(const std::string& message) 
        : TochnogException("Solver Error: " + message) {}
};

// Input/Output related exceptions
class IOException : public TochnogException {
public:
    explicit IOException(const std::string& message) 
        : TochnogException("I/O Error: " + message) {}
};

// Utility function to convert exceptions to the old exit behavior (for compatibility)
inline void handle_exception_for_compatibility(const std::exception& e) {
    std::cerr << "Exception caught: " << e.what() << std::endl;
    // For now, maintain compatibility by exiting as before
    // In future versions, this could be configured
    exit(TN_EXIT_STATUS);
}

// Macro to wrap potentially throwing functions for compatibility
#define TOCHNOG_TRY_BEGIN try {
#define TOCHNOG_CATCH_END } \
    catch (const DatabaseException& e) { handle_exception_for_compatibility(e); } \
    catch (const MemoryException& e) { handle_exception_for_compatibility(e); } \
    catch (const SolverException& e) { handle_exception_for_compatibility(e); } \
    catch (const IOException& e) { handle_exception_for_compatibility(e); } \
    catch (const std::exception& e) { handle_exception_for_compatibility(e); } \
    catch (...) { \
        std::cerr << "Unknown exception caught" << std::endl; \
        exit(TN_EXIT_STATUS); \
    }

#endif // TOCHNOG_EXCEPTIONS_H