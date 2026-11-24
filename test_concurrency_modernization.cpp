/*
 * Test for Concurrency Modernization implementation
 * Tests std::thread and RAII synchronization primitives instead of pthread
 */

#include <iostream>
#include <thread>
#include <mutex>
#include <vector>
#include <chrono>
#include <cassert>

class ModernConcurrencyTest {
private:
    int shared_counter = 0;
    std::mutex counter_mutex;
    
public:
    void increment_counter() {
        std::lock_guard<std::mutex> lock(counter_mutex);
        shared_counter++;
    }
    
    int get_counter() const {
        std::lock_guard<std::mutex> lock(counter_mutex);
        return shared_counter;
    }
    
    void worker_function(int worker_id, int iterations) {
        for (int i = 0; i < iterations; ++i) {
            increment_counter();
            // Simulate some work
            std::this_thread::sleep_for(std::chrono::microseconds(100));
        }
        std::cout << "Worker " << worker_id << " completed " << iterations << " iterations\n";
    }
};

void test_modern_concurrency() {
    std::cout << "\n=== Testing Modern Concurrency (std::thread) ===" << std::endl;
    
    ModernConcurrencyTest test_obj;
    const int num_threads = 4;
    const int iterations_per_thread = 100;
    
    std::vector<std::thread> threads;
    
    // Create multiple threads using std::thread
    for (int i = 0; i < num_threads; ++i) {
        threads.emplace_back(&ModernConcurrencyTest::worker_function, 
                           &test_obj, i, iterations_per_thread);
    }
    
    // Wait for all threads to complete
    for (auto& t : threads) {
        if (t.joinable()) {
            t.join();
        }
    }
    
    // Verify the result
    int expected = num_threads * iterations_per_thread;
    int actual = test_obj.get_counter();
    std::cout << "Expected: " << expected << ", Actual: " << actual << std::endl;
    
    assert(actual == expected);
    
    std::cout << "Modern concurrency test passed!" << std::endl;
}

void test_raii_synchronization() {
    std::cout << "\n=== Testing RAII Synchronization ===" << std::endl;
    
    std::mutex test_mutex;
    int shared_resource = 0;
    
    auto worker_lambda = [&test_mutex, &shared_resource](int id) {
        for (int i = 0; i < 10; ++i) {
            // RAII: lock_guard automatically acquires and releases the lock
            std::lock_guard<std::mutex> lock(test_mutex);
            shared_resource += id;
            // lock automatically released when going out of scope
        }
        std::cout << "Thread " << id << " completed\n";
    };
    
    std::thread t1(worker_lambda, 1);
    std::thread t2(worker_lambda, 2);
    
    if (t1.joinable()) t1.join();
    if (t2.joinable()) t2.join();
    
    std::cout << "Final shared resource value: " << shared_resource << std::endl;
    std::cout << "RAII synchronization test passed!" << std::endl;
}

int main() {
    std::cout << "Testing Concurrency Modernization implementation..." << std::endl;
    
    test_modern_concurrency();
    test_raii_synchronization();
    
    std::cout << "\nAll concurrency modernization tests passed!" << std::endl;
    return 0;
}