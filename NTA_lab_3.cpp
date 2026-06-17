#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <numeric>
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

// EEA
long long EEA(long long a, long long b, long long& x, long long& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    long long x1, y1;
    long long d = EEA(b, a % b, x1, y1);
    x = y1;
    y = x1 - y1 * (a / b);
    return d;
}

long long mod_inv(long long a, long long m) {
    long long x, y;
    long long g = EEA((a % m + m) % m, m, x, y);
    if (g != 1) return -1;
    return (x % m + m) % m;
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
vector<Equation> generateSLE(long long p, BigInt alpha, const vector<int>& S, int extra_equations = 10) {
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

// Solving SLE with reduction
vector<long long> solve_sle(const vector<Equation>& equations, long long n, size_t t) {
    size_t num_eq = equations.size();
    // [equations] x [t + 1]
    vector<vector<long long>> mat(num_eq, vector<long long>(t + 1, 0));
    for (size_t i = 0; i < num_eq; ++i) {
        for (size_t j = 0; j < t; ++j) mat[i][j] = equations[i].coefficients[j];
        mat[i][t] = equations[i].k;
    }
    size_t r = 0;
    for (size_t c = 0; c < t && r < num_eq; ++c) {
        size_t pivot_row = r;
        while (pivot_row < num_eq && mat[pivot_row][c] == 0) pivot_row++;
        if (pivot_row == num_eq) continue;
        swap(mat[r], mat[pivot_row]);
        for (size_t i = r + 1; i < num_eq; ++i) {
            while (mat[i][c] != 0) {
                long long q = mat[r][c] / mat[i][c];
                for (size_t col = c; col <= t; ++col) {
                    mat[r][col] = (mat[r][col] - q * mat[i][col]) % n;
                    if (mat[r][col] < 0) mat[r][col] += n;
                }
                swap(mat[r], mat[i]);
            }
        }
        r++;
    }
    // Back substitution
    vector<long long> x(t, -1); // logs
    for (int i = static_cast<int>(r) - 1; i >= 0; --i) {
        int lead_c = -1;
        for (size_t j = 0; j < t; ++j) {
            if (mat[i][j] != 0) {
                lead_c = j;
                break;
            }
        }
        if (lead_c == -1) continue;
        long long rhs = mat[i][t];
        for (size_t j = lead_c + 1; j < t; ++j) {
            if (x[j] != -1 && mat[i][j] != 0) {
                rhs = (rhs - mat[i][j] * x[j]) % n;
                if (rhs < 0) rhs += n;
            }
        }
        long long coeff = mat[i][lead_c];
        long long inv = mod_inv(coeff, n);
        if (inv != -1) x[lead_c] = (rhs * inv) % n;
        else {
            long long g = gcd(coeff, n);
            if (rhs % g == 0) {
                long long reduced_coeff = coeff / g;
                long long reduced_rhs = rhs / g;
                long long reduced_n = n / g;
                long long reduced_inv = mod_inv(reduced_coeff, reduced_n);
                if (reduced_inv != -1) x[lead_c] = (reduced_rhs * reduced_inv) % reduced_n;
            }
        }
    }
    return x; // x[i] = log_alpha(p_i)
}

long long log_beta(BigInt beta, BigInt alpha, long long p, const vector<int>& S, const vector<long long>& log_S) {
    long long n = p - 1;
    BigInt mod_p(p);
    random_device rd;
    mt19937_64 gen(rd());
    uniform_int_distribution<long long> dis(0, n - 1);
    while (true) {
        long long l = dis(gen);
        BigInt exp_l(l);
        BigInt alpha_l = alpha.GorMod(exp_l, mod_p); // alpha^l mod p
        BigInt target = (beta * alpha_l) % mod_p; // beta * alpha^l mod p
        vector<int> d;
        if (check_smoothness(target, S, d)) {
            cout << "l = " << l << "\n";
            // sum(d_i * log_alpha(p_i)) mod n
            long long sum_logs = 0;
            for (size_t i = 0; i < S.size(); ++i) {
                if (d[i] > 0) {
                    // d_i * log_S[i]
                    long long term = (static_cast<long long>(d[i]) * log_S[i]) % n;
                    sum_logs = (sum_logs + term) % n;
                }
            }
            long long log_beta = (sum_logs - l) % n; // log_alpha(beta) = (sum_logs - l) mod n
            if (log_beta < 0) log_beta += n;
            return log_beta;
        }
    }
}

int main()
{
    
}