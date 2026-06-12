#include <iostream>
#include "BigInt.hpp"
#include <chrono>
using namespace std;
using namespace chrono;

BigInt direct_discretelog(const BigInt& a, const BigInt& b, const BigInt& p, int timeout_seconds) {
    if (p == 0) return -1;
    if (b >= p) return -1;
    if (b == 1) return 0;
    auto start_time = steady_clock::now();
    auto max_duration = seconds(timeout_seconds);
    BigInt current = 1;
    // 5 min -> uint64_t for x
    for (uint64_t x = 0; BigInt(x) < p; ++x) {
        if (x % 100000 == 0 && x > 0) {
            auto current_time = steady_clock::now();
            if (current_time - start_time >= max_duration) return -2;
        }
        if (current == b) return BigInt(x);
        current = (current * a) % p;
        if (x > 0 && current == 1) break;
    }
    return -1;
}

int main()
{
    cout << "=== Direct search ===" << endl;
    BigInt a1, b1, p1;
    a1.fromDecimalString("5");
    b1.fromDecimalString("94512");
    p1.fromDecimalString("100003");
    int timeout1 = 10; // Limit = 10 sec
    cout << "Test 1: 5^x == 94512 (mod 100003)" << endl;
    auto start1 = high_resolution_clock::now();
    BigInt result1 = direct_discretelog(a1, b1, p1, timeout1);
    auto end1 = high_resolution_clock::now();
    duration<double> elapsed1 = end1 - start1;
    if (result1 == BigInt(-1)) {
        cout << "There is no solution in this group" << endl;
    }
    else if (result1 == BigInt(-2)) {
        cout << "Time limit exceeded (" << timeout1 << " sec)." << endl;
    }
    else {
        cout << "x = " << result1.toDecimalString() << endl;
    }
    cout << "Runtime: " << elapsed1.count() << " sec." << endl;
    cout << "--------------------------------------------" << endl;

    BigInt a2, b2, p2;
    a2.fromDecimalString("2");
    b2.fromDecimalString("12345678901234567890");
    p2.fromDecimalString("1000000000000000000000000000003");
    int timeout2 = 10;
    cout << "Test 2: 2^x == 12345678901234567890 (mod 1000000000000000000000000000003)" << endl;
    auto start2 = high_resolution_clock::now();
    BigInt result2 = direct_discretelog(a2, b2, p2, timeout2);
    auto end2 = high_resolution_clock::now();
    duration<double> elapsed2 = end2 - start2;
    if (result2 == BigInt(-2)) {
        cout << "Time limit exceeded (" << timeout2 << " sec)." << endl;
    }
    else if (result2 == BigInt(-1)) {
        cout << "There is no solution in this group" << endl;
    }
    else {
        cout << "x = " << result2.toDecimalString() << endl;
    }
    cout << "Runtime: " << elapsed2.count() << " sec." << endl;
    cout << "--------------------------------------------" << endl;

    return 0;
}
