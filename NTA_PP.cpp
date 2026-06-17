#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <random>
#include <algorithm>
#include <fstream>
using namespace std;

using Vector = vector<double>;
using Matrix = vector<Vector>;

double c_norm(const Vector& v) {
    double sum = 0.0;
    for (double x : v) sum += x * x;
    return sqrt(sum);
}

double get_min_row_norm(const Matrix& M) {
    if (M.empty()) return 0.0;
    double min_n = c_norm(M[0]);
    for (size_t i = 1; i < M.size(); ++i) {
        double current_norm = c_norm(M[i]);
        if (current_norm < min_n) min_n = current_norm;
    }
    return min_n;
}

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

// Gauss method
bool is_full_rank(Matrix M) {
    int n = M.size();
    for (int i = 0; i < n; ++i) {
        // main element in the column i
        int pivot = i;
        for (int j = i + 1; j < n; ++j) {
            if (abs(M[j][i]) > abs(M[pivot][i])) {
                pivot = j;
            }
        }
        // (Max element -> 0), -> It isn't a full-rank mmatrix
        if (abs(M[pivot][i]) < 1e-9) return false;
        if (pivot != i) swap(M[i], M[pivot]);
        for (int j = i + 1; j < n; ++j) {
            double factor = M[j][i] / M[i][i];
            for (int k = i; k < n; ++k) {
                M[j][k] -= factor * M[i][k];
            }
        }
    }
    return true;
}

// Random Full-Rank matrix
Matrix generate_rfr_matrix(int n, int min_val, int max_val, mt19937& gen) {
    uniform_int_distribution<> distr(min_val, max_val);
    Matrix M(n, Vector(n));
    while (true) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) M[i][j] = distr(gen);
        }
        if (is_full_rank(M)) return M;
    }
}

void matrix_to_file(ofstream& out, const Matrix& M) {
    for (const auto& row : M) {
        for (double val : row) out << setw(4) << static_cast<int>(val) << " ";
        out << "\n";
    }
}

int main() {
    random_device rd;
    mt19937 gen(rd());
    const int num_matrices = 50;
    const int n = 30;
    const int min_element = -5;
    const int max_element = 5;
    const double delta = 0.75;
    ofstream out_file("lll_50_matrixes.txt");
    if (!out_file.is_open()) {
        cerr << "ERROR!!! Cannot open the file." << std::endl;
        return 1;
    }
    out_file << "Results for the " << num_matrices << " random matrices (" << n << "x" << n << ")\n";
    out_file << "delta = " << delta << "\n";
    for (int i = 1; i <= num_matrices; ++i) {
        Matrix B = generate_rfr_matrix(n, min_element, max_element, gen);
        double min_norm_before = get_min_row_norm(B);
        Matrix reduced_B = lll_reduction(B, delta);
        double min_norm_after = get_min_row_norm(reduced_B);
        out_file << "Matrix #" << i << "\n";
        out_file << "Min norm before hte reduction:   " << fixed << setprecision(4) << min_norm_before << "\n";
        out_file << "Min norm after the reduction: " << min_norm_after << "\n\n";
        out_file << "Basis matrix:\n";
        matrix_to_file(out_file, B);
        out_file << "\n";
        out_file << "LLL-reduced basis:\n";
        matrix_to_file(out_file, reduced_B);
        out_file << "\n\n";
    }
    out_file.close();
    cout << "\nResults are in the file lll_50_matrixes.txt" << endl;
    return 0;
}