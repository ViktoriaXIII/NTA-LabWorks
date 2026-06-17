#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

using Vector = vector<double>;
using Matrix = vector<Vector>;

void basis_mu(int k, const Matrix& B, Matrix& B_star, Matrix& mu) {
    int n = B[k].size();
    B_star[k] = B[k];
    for (int j = 0; j < k; ++j) {
        double num = 0.0;
        double den = 0.0;
        for (int r = 0; r < n; ++r) {
            num += B[k][r] * B_star[j][r];
            den += B_star[j][r] * B_star[j][r];
        }
        mu[k][j] = (den == 0.0) ? 0.0 : num / den; // \mu_{k,j}
        // b_k^* := b_k^* - \mu_{k,j} * b_j^*
        for (int r = 0; r < n; ++r) B_star[k][r] -= mu[k][j] * B_star[j][r];
    }
}

Matrix lll_reduction(Matrix B, double delta) {
    int m = B.size();
    if (m == 0) return B;
    int n = B[0].size();
    // Matrix for the b^*-basis
    Matrix B_star(m, Vector(n, 0.0));
    // Matrix for the \mu (m x m)
    Matrix mu(m, Vector(m, 0.0));
    // 1.Build the basis
    for (int i = 0; i < m; ++i) basis_mu(i, B, B_star, mu);
    // 2.
    int k = 1;
    while (k < m) {
        // 2.1. Size Reduction
        for (int j = k - 1; j >= 0; --j) {
            if (abs(mu[k][j]) > 0.5) {
                double q = round(mu[k][j]);
                for (int r = 0; r < n; ++r) {
                    B[k][r] -= q * B[j][r]; // b_k := b_k - q * b_j
                }
                basis_mu(k, B, B_star, mu);
            }
        }
        // ||b_k^*||^2 та ||b_{k-1}^*||^2
        double b_star_k_norm2 = 0.0;
        double b_star_k_minus_1_norm2 = 0.0;
        for (int r = 0; r < n; ++r) {
            b_star_k_norm2 += B_star[k][r] * B_star[k][r];
            b_star_k_minus_1_norm2 += B_star[k - 1][r] * B_star[k - 1][r];
        }
        // 2.2. Lovasz
        if (b_star_k_norm2 < (delta - mu[k][k - 1] * mu[k][k - 1]) * b_star_k_minus_1_norm2) {
            // 2.2.1. swap b_{k-1} та b_k
            swap(B[k], B[k - 1]);
            for (int i = k - 1; i < m; ++i) basis_mu(i, B, B_star, mu);
            k = max(k - 1, 1); // 2.2.3. k = max(k - 1, 1)
        }
        else k = k + 1;
    }

    return B;
}

int main() {
    
}