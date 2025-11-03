/*
 * Test for Exception Handling implementation
 * Tests modern C++ exception handling instead of exit() calls
 */

#include <iostream>
#include <stdexcept>
#include <cassert>

// Define custom exception hierarchy for Tochnog
class TochnogException : public std::exception {
protected:
    std::string msg;
public:
    explicit TochnogException(const std::string& message) : msg(message) {}
    const char* what() const noexcept override { return msg.c_str(); }
};

class DatabaseException : public TochnogException {
public:
    explicit DatabaseException(const std::string& message) : TochnogException("Database Error: " + message) {}
};

class MemoryException : public TochnogException {
public:
    explicit MemoryException(const std::string& message) : TochnogException("Memory Error: " + message) {}
};

// Function that might fail and throw exception instead of calling exit()
void potentially_failing_function(bool should_fail) {
    if (should_fail) {
        throw DatabaseException("Failed to access database element");
    }
    std::cout << "Function executed successfully" << std::endl;
}

// Function that demonstrates exception propagation and handling
void demonstrate_exception_handling() {
    std::cout << "\n=== Testing Exception Handling ===" << std::endl;
    
    try {
        // Test successful case
        potentially_failing_function(false);
        
        // Test failure case - this should throw
        potentially_failing_function(true);
        
        std::cout << "This line should not be reached" << std::endl;
    }
    catch (const DatabaseException& e) {
        std::cout << "Caught expected DatabaseException: " << e.what() << std::endl;
        assert(std::string(e.what()).find("Database Error") != std::string::npos);
    }
    catch (const TochnogException& e) {
        std::cout << "Caught TochnogException: " << e.what() << std::endl;
    }
    catch (...) {
        std::cout << "Caught unexpected exception" << std::endl;
        throw; // Re-throw if not expected
    }
    
    std::cout << "Exception handling test passed!" << std::endl;
}

// Test nested exception handling
void nested_function() {
    throw MemoryException("Out of memory");
}

void test_nested_exceptions() {
    std::cout << "\n=== Testing Nested Exception Handling ===" << std::endl;
    
    try {
        nested_function();
    }
    catch (const TochnogException& e) {
        std::cout << "Caught nested exception: " << e.what() << std::endl;
        assert(std::string(e.what()).find("Memory Error") != std::string::npos);
    }
    
    std::cout << "Nested exception test passed!" << std::endl;
}

int main() {
    std::cout << "Testing Exception Handling implementation..." << std::endl;
    
    demonstrate_exception_handling();
    test_nested_exceptions();
    
    std::cout << "All exception tests passed!" << std::endl;
    return 0;
}