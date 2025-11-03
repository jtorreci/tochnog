/*
 * Test for String Modernization implementation
 * Tests std::string usage instead of C-style char arrays
 */

#include <iostream>
#include <string>
#include <cassert>

// Test safe string operations replacing strcpy/strcat
void test_safe_string_operations() {
    std::cout << "\n=== Testing Safe String Operations ===" << std::endl;
    
    // Instead of: char buffer[MCHAR];
    // Use: std::string
    
    std::string str1 = "Hello";
    std::string str2 = "World";
    
    // Safe concatenation - no buffer overflows
    std::string result = str1 + " " + str2;
    assert(result == "Hello World");
    
    std::cout << "Concatenation result: " << result << std::endl;
    
    // Safe substring operations
    std::string sub = result.substr(0, 5);  // "Hello"
    assert(sub == "Hello");
    
    std::cout << "Substring result: " << sub << std::endl;
    
    // Safe append operations
    std::string append_test = "Test";
    append_test += "_suffix";
    assert(append_test == "Test_suffix");
    
    std::cout << "Append result: " << append_test << std::endl;
    
    std::cout << "Safe string operations test passed!" << std::endl;
}

// Test string formatting functions
void test_string_formatting() {
    std::cout << "\n=== Testing String Formatting ===" << std::endl;
    
    int number = 42;
    std::string formatted = "Number: " + std::to_string(number);
    assert(formatted == "Number: 42");
    
    std::cout << "Formatted string: " << formatted << std::endl;
    
    double pi = 3.14159;
    std::string pi_str = "Pi: " + std::to_string(pi);
    std::cout << "Pi string: " << pi_str << std::endl;
    
    std::cout << "String formatting test passed!" << std::endl;
}

// Test string to char* conversion for C compatibility
void test_c_compatibility() {
    std::cout << "\n=== Testing C Compatibility ===" << std::endl;
    
    std::string modern_string = "Hello C Compatibility";
    
    // Safe conversion to C-string when needed for legacy functions
    const char* c_str = modern_string.c_str();
    assert(std::string(c_str) == modern_string);
    
    std::cout << "C-string: " << c_str << std::endl;
    std::cout << "C compatibility test passed!" << std::endl;
}

int main() {
    std::cout << "Testing String Modernization implementation..." << std::endl;
    
    test_safe_string_operations();
    test_string_formatting();
    test_c_compatibility();
    
    std::cout << "\nAll string modernization tests passed!" << std::endl;
    return 0;
}