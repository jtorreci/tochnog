/*
 * Test for RAII and Smart Pointers implementation
 * Tests resource management using RAII principle
 */

#include <iostream>
#include <memory>
#include <cassert>

class TestResource {
private:
    int* data;
    size_t size;
    static int resource_count;
    
public:
    explicit TestResource(size_t s) : size(s) {
        data = new int[size];
        for(size_t i = 0; i < size; ++i) {
            data[i] = i;
        }
        resource_count++;
        std::cout << "Resource allocated with size: " << size << std::endl;
    }
    
    // RAII: Destructor automatically frees resources
    ~TestResource() {
        delete[] data;
        resource_count--;
        std::cout << "Resource freed" << std::endl;
    }
    
    // Disable copy to follow RAII best practices
    TestResource(const TestResource&) = delete;
    TestResource& operator=(const TestResource&) = delete;
    
    // Enable move for efficiency
    TestResource(TestResource&& other) noexcept 
        : data(other.data), size(other.size) {
        other.data = nullptr;
        other.size = 0;
    }
    
    TestResource& operator=(TestResource&& other) noexcept {
        if (this != &other) {
            delete[] data;
            data = other.data;
            size = other.size;
            other.data = nullptr;
            other.size = 0;
        }
        return *this;
    }
    
    int* get() const { return data; }
    size_t get_size() const { return size; }
    
    static int get_resource_count() { return resource_count; }
};

int TestResource::resource_count = 0;

// Test RAII with smart pointers
void test_smart_pointers() {
    std::cout << "\n=== Testing Smart Pointers ===" << std::endl;
    
    {
        // unique_ptr automatically manages memory
        std::unique_ptr<TestResource> ptr = std::make_unique<TestResource>(10);
        assert(ptr->get_size() == 10);
        std::cout << "unique_ptr test passed" << std::endl;
    } // ptr is automatically destroyed here (RAII)
    
    std::cout << "Resource count after unique_ptr scope: " << TestResource::get_resource_count() << std::endl;
    
    {
        // shared_ptr for shared ownership
        std::shared_ptr<TestResource> ptr1 = std::make_shared<TestResource>(5);
        {
            std::shared_ptr<TestResource> ptr2 = ptr1; // Shared ownership
            std::cout << "Reference count: " << ptr1.use_count() << std::endl;
            assert(ptr1.use_count() == 2);
        } // ptr2 goes out of scope, but resource still exists
        std::cout << "Reference count after inner scope: " << ptr1.use_count() << std::endl;
        assert(ptr1.use_count() == 1);
    } // ptr1 goes out of scope, resource is destroyed
    std::cout << "Resource count after shared_ptr scope: " << TestResource::get_resource_count() << std::endl;
}

int main() {
    std::cout << "Testing RAII and Smart Pointers implementation..." << std::endl;
    
    test_smart_pointers();
    
    std::cout << "All RAII tests passed!" << std::endl;
    return 0;
}