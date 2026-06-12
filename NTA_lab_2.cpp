#include <iostream>
#include "BigInt.hpp"
#include "CanonForm.hpp"
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

// p-1 = q_1^e_1 * q_2^e_2 * ... * q_n^e_n
// Struct for q_i^e_i
struct PrimeFactor {
    BigInt q;
    size_t e;
    BigInt q_e;
};

BigInt discretelog_SPH(const BigInt& a, const BigInt& b, const BigInt& p, const vector<PrimeFactor>& factors) {
    BigInt zero, one;
    zero.fromDecimalString("0");
    one.fromDecimalString("1");
    BigInt group_order = p - one;
    vector<BigInt> remainders;
    vector<BigInt> moduli;
    for (const auto& factor : factors) {
        BigInt x_val = zero; // = X_prev
        BigInt current_q_power = one;
        // gamma = a^(group_order / q) mod p
        BigInt gamma_exp = group_order / factor.q;
        BigInt gamma = a.GorMod(gamma_exp, p);
        for (size_t l = 0; l < factor.e; ++l) {
            BigInt exp_denominator = current_q_power * factor.q;
            // exp = n / q_i^l
            BigInt current_exp = group_order / exp_denominator;
            // (beta * alpha^(-x_0 - ... - x_(k-1)q_i^(k-1)))^q_i^(n/(k+1))$
            // g = a^exp mod p
            BigInt g = a.GorMod(current_exp, p);
            BigInt g_order = current_q_power * factor.q;
            // g^(-1) = g^(O-1) mod p
            BigInt g_inv = g.GorMod(g_order - one, p);
            // alpha^exp^(-X_prev) mod p
            BigInt correction = g_inv.GorMod(x_val, p);
            // beta_scaled = beta^exp mod p
            BigInt b_scaled = b.GorMod(current_exp, p);
            // delta = (beta * alpha^(-X_prev))^exp mod p
            BigInt delta = (b_scaled * correction) % p;
            // Searching x_l: gamma^x_l = delta mod p
            BigInt x_l = zero;
            BigInt gamma_pow = one;
            while (!(gamma_pow == delta)) {
                gamma_pow = (gamma_pow * gamma) % p;
                x_l = x_l + one;
            }
            //  y_i = x mod q_i^{l_i}}$
            x_val = x_val + (x_l * current_q_power);
            current_q_power = current_q_power * factor.q;
        }
        remainders.push_back(x_val);
        moduli.push_back(factor.q_e);
    }
    // CRT
    BigInt global_x = zero;
    for (size_t i = 0; i < moduli.size(); ++i) {
        BigInt M_i = group_order / moduli[i];
        // phi(q^e) = q^e - q^(e-1)
        BigInt q_minus_1 = factors[i].q - one;
        BigInt phi = q_minus_1 * (factors[i].q_e / factors[i].q);
        // N_i = M_i^(-1) mod moduli[i]
        BigInt N_i = M_i.GorMod(phi - one, moduli[i]);
        BigInt term = (remainders[i] * M_i) % group_order;
        term = (term * N_i) % group_order;
        global_x = (global_x + term) % group_order;
    }
    return global_x;
}

int main()
{
    /*cout << "=== Direct search ===" << endl;
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

    BigInt a, b, p;
    a.fromDecimalString("3");
    b.fromDecimalString("215");
    p.fromDecimalString("1013");
    cout << "Test 3: 3^x == 215 (mod 1013)" << endl;
    // p - 1 = 1012 = 2^2 * 11^1 * 23^1
    vector<PrimeFactor> factors(3);
    factors[0].q.fromDecimalString("2");
    factors[0].e = 2;
    factors[0].q_e.fromDecimalString("4"); // 2^2
    factors[1].q.fromDecimalString("11");
    factors[1].e = 1;
    factors[1].q_e.fromDecimalString("11"); // 11^1
    factors[2].q.fromDecimalString("23");
    factors[2].e = 1;
    factors[2].q_e.fromDecimalString("23"); // 23^1
    cout << "SPH..." << endl;
    auto start_sph = chrono::high_resolution_clock::now();
    BigInt res_sph = discretelog_SPH(a, b, p, factors);
    auto end_sph = chrono::high_resolution_clock::now();
    chrono::duration<double> time_sph = end_sph - start_sph;
    cout << "x = " << res_sph.toDecimalString() << endl;
    cout << "Runtime: " << time_sph.count() << " sec." << endl;
    cout << "--------------------------------------------" << endl;*/

    BigInt a, b, p;
    a.fromDecimalString("2402221");
    b.fromDecimalString("3809591");
    p.fromDecimalString("8582963");
    cout << "Test SPH" << endl;
    BigInt one;
    one.fromDecimalString("1");
    BigInt group_order = p - one;
    map<BigInt, int> raw_factors;
    cout << "Factoring p - 1..." << endl;
    findCannonicalFactorization(group_order, raw_factors);
    vector<PrimeFactor> factors;
    for (const auto& factor : raw_factors) {
        PrimeFactor pf;
        pf.q = factor.first;
        pf.e = factor.second;
        BigInt exp_bigint;
        exp_bigint.fromDecimalString(to_string(pf.e));
        pf.q_e = pf.q.GorMod(exp_bigint, p);
        factors.push_back(pf);
    }
    cout << "SPH..." << endl;
    auto start_sph = high_resolution_clock::now();
    BigInt res_sph = discretelog_SPH(a, b, p, factors);
    auto end_sph = high_resolution_clock::now();
    duration<double> time_sph = end_sph - start_sph;
    cout << "x = " << res_sph.toDecimalString() << endl;
    cout << "Runtime: " << time_sph.count() << " sec." << endl;
    cout << "--------------------------------------------" << endl;
    return 0;
}
