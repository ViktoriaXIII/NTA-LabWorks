#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include "BigInt.hpp"
#include "CanonForm.hpp"
using namespace std;

struct FactorBase {
    double B;
    vector<int> S;
};

// k = c_1*log(p_1) + c_2*log(p_2) + ... + c_t*log(p_t) (mod n)
struct Equation {
    long long k;
    vector<int> coefficients;
};

// Sieve of Eratosthenes
vector<int> sieve(double limit) {
    if (limit < 2) return {};
    int limit_int = static_cast<int>(limit);
    vector<bool> is_prime(limit_int + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (int i = 2; i * i <= limit_int; ++i) {
        if (is_prime[i]) {
            for (int j = i * i; j <= limit_int; j += i) is_prime[j] = false;
        }
    }
    vector<int> primes;
    for (int i = 2; i <= limit_int; ++i) {
        // p < B
        if (is_prime[i] && i < limit) primes.push_back(i);
    }
    return primes;
}

FactorBase build_factor_base(double n, double c = 3.38) {
    double log_n = log(n);
    double log_log_n = log(log_n);
    // 0.5 * sqrt(log(n) * log(log(n)))
    double exponent = 0.5 * sqrt(log_n * log_log_n);
    double B = c * exp(exponent);
    vector<int> S = sieve(B);
    return { B, S };
}

bool check_smoothness(long long value, const vector<int>& S, vector<int>& coefficients) {
    coefficients.assign(S.size(), 0);
    long long temp = value;
    for (size_t i = 0; i < S.size(); ++i) {
        while (temp % S[i] == 0) {
            coefficients[i]++;
            temp /= S[i];
        }
    }
    return (temp == 1);
}

int main()
{
    
}