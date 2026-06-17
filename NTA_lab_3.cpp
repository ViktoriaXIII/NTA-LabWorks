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
    BigInt k;
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

bool check_smoothness(BigInt value, const vector<int>& S, vector<int>& coefficients) {
    coefficients.assign(S.size(), 0);
    BigInt temp = value;
    BigInt zero(0);
    BigInt one(1);
    for (size_t i = 0; i < S.size(); ++i) {
        BigInt prime_big(S[i]);
        while (temp % S[i] == 0) {
            coefficients[i]++;
            temp = temp / prime_big;
        }
    }
    return (temp == one);
}

// SLE
vector<Equation> generate_equations(long long p, BigInt alpha, const vector<int>& S, int extra_equations = 10) {
    long long n = p - 1; // Z_p*
    size_t t = S.size();
    size_t required_equations = t + extra_equations; // t + c equations
    vector<Equation> equations;
    BigInt mod_p(p);
    random_device rd;
    mt19937_64 gen(rd());
    uniform_int_distribution<long long> dis(1, n - 1);
    cout << "Starting the search for smoth numbers. Creating " << required_equations << " equations\n";
    while (equations.size() < required_equations) {
        long long k = dis(gen); // k
        BigInt exp_k(k);
        BigInt val = alpha.GorMod(k, p); // Обчислюємо alpha^k mod p
        vector<int> coefficients;
        if (check_smoothness(val, S, coefficients)) {
            Equation eq;
            eq.k = k;
            eq.coefficients = coefficients;
            equations.push_back(eq);
            if (equations.size() % 5 == 0 || equations.size() == required_equations) {
                cout << "Found " << equations.size() << " / " << required_equations << " equations...\n";
            }
        }
    }

    return equations;
}

int main()
{
    
}