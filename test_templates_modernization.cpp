/*
 * Test for Templates and Generics implementation
 * Tests modern C++ templates replacing macros and duplicated code
 */

#include <iostream>
#include <vector>
#include <cassert>
#include <type_traits>

// Generic container wrapper using templates
template<typename T>
class GenericContainer {
private:
    std::vector<T> data;

public:
    void add(const T& item) {
        data.push_back(item);
    }
    
    const T& get(size_t index) const {
        return data.at(index);
    }
    
    size_t size() const {
        return data.size();
    }
    
    // Generic algorithm to find an element
    template<typename Predicate>
    int find_index(Predicate pred) const {
        for (size_t i = 0; i < data.size(); ++i) {
            if (pred(data[i])) {
                return static_cast<int>(i);
            }
        }
        return -1;
    }
};

// Generic function similar to what might replace macro-based code
template<typename T>
T array_max(const std::vector<T>& arr) {
    if (arr.empty()) return T{};
    T max_val = arr[0];
    for (const auto& val : arr) {
        if (val > max_val) {
            max_val = val;
        }
    }
    return max_val;
}

template<typename T>
T array_min(const std::vector<T>& arr) {
    if (arr.empty()) return T{};
    T min_val = arr[0];
    for (const auto& val : arr) {
        if (val < min_val) {
            min_val = val;
        }
    }
    return min_val;
}

// Generic swap function (replacing potential macro)
template<typename T>
void generic_swap(T& a, T& b) {
    T temp = std::move(a);
    a = std::move(b);
    b = std::move(temp);
}

void test_generic_containers() {
    std::cout << "\n=== Testing Generic Containers ===" << std::endl;
    
    GenericContainer<int> int_container;
    int_container.add(10);
    int_container.add(20);
    int_container.add(30);
    
    assert(int_container.size() == 3);
    assert(int_container.get(1) == 20);
    
    // Test generic find with lambda
    auto index = int_container.find_index([](int val) { return val == 20; });
    assert(index == 1);
    
    std::cout << "Generic container test passed!" << std::endl;
}

void test_generic_algorithms() {
    std::cout << "\n=== Testing Generic Algorithms ===" << std::endl;
    
    std::vector<int> int_vec = {5, 2, 8, 1, 9};
    std::vector<double> double_vec = {5.5, 2.2, 8.8, 1.1, 9.9};
    
    assert(array_max(int_vec) == 9);
    assert(array_min(int_vec) == 1);
    assert(array_max(double_vec) > 9.8 && array_max(double_vec) < 10.0);
    assert(array_min(double_vec) > 1.0 && array_min(double_vec) < 1.2);
    
    std::cout << "Max int: " << array_max(int_vec) << std::endl;
    std::cout << "Min int: " << array_min(int_vec) << std::endl;
    std::cout << "Max double: " << array_max(double_vec) << std::endl;
    std::cout << "Min double: " << array_min(double_vec) << std::endl;
    
    std::cout << "Generic algorithms test passed!" << std::endl;
}

void test_generic_swap() {
    std::cout << "\n=== Testing Generic Swap ===" << std::endl;
    
    int a = 5, b = 10;
    std::cout << "Before swap: a=" << a << ", b=" << b << std::endl;
    generic_swap(a, b);
    std::cout << "After swap: a=" << a << ", b=" << b << std::endl;
    
    assert(a == 10 && b == 5);
    
    double x = 3.14, y = 2.71;
    std::cout << "Before swap: x=" << x << ", y=" << y << std::endl;
    generic_swap(x, y);
    std::cout << "After swap: x=" << x << ", y=" << y << std::endl;
    
    assert(x > 2.70 && x < 2.72 && y > 3.13 && y < 3.15);
    
    std::cout << "Generic swap test passed!" << std::endl;
}

int main() {
    std::cout << "Testing Templates and Generics implementation..." << std::endl;
    
    test_generic_containers();
    test_generic_algorithms();
    test_generic_swap();
    
    std::cout << "\nAll templates and generics tests passed!" << std::endl;
    return 0;
}