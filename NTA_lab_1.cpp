#include "BigInt.hpp"
using namespace std;

bool MillerRabin(BigInt p, int k = 10) {
	if (p == BigInt(2) || p == BigInt(3)) return true;
	if (p < BigInt(2) || (p % BigInt(2) == BigInt(0))) return false;
	// Step 0. p - 1 = d * 2^s
	BigInt p_1 = p - BigInt(1);
	BigInt d = p_1;
	size_t s = 0;
	while (d % BigInt(2) == BigInt(0)) {
		d = d / BigInt(2);
		s++;
	}
	int counter = 0;
	vector<string> bases = { "2" , "3", "5", "7", "11", "13", "17", "19", "23", "29" };
	while (counter < k) {
		// Step 1
		BigInt x;
		if (counter < bases.size()) x.fromDecimalString(bases[counter]);
		else x = BigInt::generateRandom(1) % (p - BigInt(4)) + BigInt(2);
		if (x >= p_1) {
			counter++;
			k++;
			continue;
		}
		BigInt g = BigInt::gcd(x, p);
		if (g > BigInt(1)) return false;
		// Step 2
		bool isStornglyPseudoPrime = false;
		// x^d mod p = +-1
		BigInt val = x.GorMod(d, p);
		if (val == BigInt(1) || val == p_1) isStornglyPseudoPrime = true;
		else {
			BigInt x_r = val;
			for (size_t r = 1; r < s; r++) {
				// x_r = x^(d * 2^r) mod p
				x_r = x_r.mulMod(x_r, p);
				if (x_r == p_1) {
					isStornglyPseudoPrime = true;
					break;
				}
				if (x_r == BigInt(1)) {
					isStornglyPseudoPrime = false;
					break;
				}
			}
		}
		// Step 3
		if (!isStornglyPseudoPrime) return false;
		// Step 4
		counter++;
	}
	return true;
}

bool Pascal(const BigInt& n, uint32_t m) {
	if (m == 0) return false;
	if (m == 1) return true;
	uint64_t sum = 0;
	uint64_t r_i = 1; // r_0 = 1
	uint64_t B_mod_m = 0x100000000ULL % m;
	for (size_t i = 0; i < n.getDigitsCount(); i++) {
		uint64_t a_i = n.getDigit(i);
		// n = sum(a_i * r_i (mod m))
		uint64_t term = (a_i % m * r_i) % m;
		sum = (sum + term) % m;
		// r_{i + 1} = (r_i * B) mod m
		r_i = (r_i * B_mod_m) % m;
	}
	return static_cast<uint32_t>(sum);
}

void trialDiv(BigInt n) {
	BigInt temp = n;
	vector<uint32_t> factors;
	while (temp > BigInt(1) && (temp.getDigit(0) % 2 == 0)) {
		factors.push_back(2);
		temp = temp / BigInt(2);
	}
	uint32_t limit = 100000;
	for (uint32_t m = 3; m <= limit; m += 2) {
		if (Pascal(temp, m) == 0) {
			BigInt divisor(m);
			while (temp > BigInt(1) && (temp % divisor == BigInt(0))) {
				factors.push_back(m);
				temp = temp / divisor;
			}
		}
		if (temp == BigInt(1)) break;
	}
	if (factors.empty()) cout << "No factors found by Trial Division up to " << limit << endl;
	else {
		cout << "Prime factors found: ";
		for (size_t i = 0; i < factors.size(); i++) cout << factors[i] << (i == factors.size() - 1 ? "" : " * ");
		cout << endl;
	}
	if (temp > BigInt(1)) cout << "Remaining unfactored part: " << temp.toDecimalString() << endl;
}

BigInt rhoPollard(const BigInt& n) {
	if (n % BigInt(2) == BigInt(0)) return BigInt(2);
	// Step 1
	BigInt x(2), y(2), d(1), c(1);
	// f(x) = (x^2 + c) mod n;
	auto f = [&](const BigInt& val, const BigInt& mod, const BigInt c) {
		return (val.mulMod(val, mod) + c) % mod;
		};
	while (true) {
		// Step 2
		x = f(x, n, c);
		y = f(f(y, n, c), n, c);
		BigInt diff = (x > y) ? (x - y) : (y - x);
		d = BigInt::gcd(diff, n);
		// Step 3
		if (x == y) {
			cout << "x = y with c = " << c.toDecimalString() << "." << endl;
			c = c + BigInt(1);
			x = BigInt(2);
			y = BigInt(2);
			d = BigInt(1);
			if (c > BigInt(100)) {
				cout << "Rho-Pollard failed to find a factor (cycle detected)" << endl;
			}
		}
		// Step 4
		if (d != BigInt(1)) {
			cout << "Found non-trivial factor: " << d.toDecimalString() << endl;
			return d;
		}
	}
}

int Legendre(BigInt n, BigInt p) {
	if (n % p == BigInt(0)) return 0;
	BigInt n_mod_p = n % p;
	BigInt result = n_mod_p.GorMod(((p - BigInt(1)) / BigInt(2)), p);
	if (result == BigInt(1)) return 1;
	else return -1;
}

vector<int> factorBase(const BigInt& n) {
	vector<int> base = { -1 };
	int input;
	cout << "~~~Formation of the factor base for n = " << n.toDecimalString() << endl;
	cout << "Enter prime number (or 0 to complete):" << endl;
	while (true) {
		cout << "Current factor base B = { ";
		for (size_t i = 0; i < base.size(); i++) cout << base[i] << (i == base.size() - 1 ? "" : ", ");
		cout << " }" << endl;
		cout << "Enter p: ";
		cin >> input;
		if (input == 0) {
			if (base.size() < 2) {
				cout << "Base is too small. Add more numbers" << endl;
				continue;
			}
			cout << "Base formation completed" << endl;
			break;
		}
		BigInt p(input);
		int symbol = Legendre(n, p);
		if (symbol == 1) {
			bool exists = false;
			for (int b : base) if (b == input) exists = true;
			if (!exists) {
				base.push_back(input);
			}
			else cout << "This number has already been added. Try another one" << endl;
		}
		else if (symbol == 0) cout << "Error!!! p divides n. Try another one" << endl;
		else cout << "Error!!! (n/" << input << ") = -1. Try another one" << endl;
	}
	return base;
}

bool tryFactorize(BigInt val, const vector<int>& base, vector<int>& degrees) {
	degrees.assign(base.size(), 0);
	if (val < BigInt(0)) {
		degrees[0] = 1;
		val = val * BigInt(-1);
	}
	for (size_t i = 1; i < base.size(); i++) {
		BigInt p(base[i]);
		while (val > BigInt(0) && val % p == BigInt(0)) {
			degrees[i]++;
			val = val / p;
		}
	}
	return (val == BigInt(1));
}

bool zeroVector(const vector<int>& degrees) {
	for (int deg : degrees) if (deg % 2 != 0) return false;
	return true;
}

void solution(const BigInt& n, const vector<BigInt>& selected_b, const vector<vector<int>>& selected_degrees, const vector<int>& base) {
	BigInt X(1);
	for (const auto& b : selected_b) X = (X * b) % n;
	vector<int> total_degrees(base.size(), 0);
	for (const auto& vec : selected_degrees) {
		for (size_t j = 0; j < base.size(); j++) total_degrees[j] += vec[j];
	}
	BigInt Y(1);
	for (size_t j = 0; j < base.size(); j++) {
		if (total_degrees[j] > 0) {
			BigInt p_val = (base[j] == -1) ? (n - BigInt(1)) : BigInt(base[j]);
			BigInt degree(total_degrees[j] / 2);
			Y = (Y * p_val.GorMod(degree, n)) % n;
		}
	}
	BigInt Y_neg = n - Y;
	if (X % n != Y % n && X % n != Y_neg % n) {
		BigInt d = BigInt::gcd(X - Y, n);
		if (d > BigInt(1) && d < n) cout << "Divisor found: " << d.toDecimalString() << endl;
		else cout << "A trivial divisor was found" << endl;
	}
	else cout << "X = +-Y (mod n). Looking for another solution" << endl;
}

vector<vector<int>> buildGF2Matrix(const vector<vector<int>>& degree_vectors) {
	if (degree_vectors.empty()) return {};
	size_t rows = degree_vectors[0].size();
	size_t cols = degree_vectors.size();
	vector<vector<int>> matrix(rows, vector<int>(cols));
	for (size_t i = 0; i < rows; ++i) {
		for (size_t j = 0; j < cols; ++j) {
			matrix[i][j] = degree_vectors[j][i] % 2;
		}
	}
	return matrix;
}

vector<vector<int>> solveGaussianGF2(vector<vector<int>> matrix) {
	int rows = matrix.size();
	int cols = matrix[0].size();
	vector<int> pivot_row(cols, -1);
	int current_row = 0;
	for (int j = 0; j < cols && current_row < rows; ++j) {
		int sel = current_row;
		while (sel < rows && matrix[sel][j] == 0) sel++;
		if (sel == rows) continue;
		swap(matrix[sel], matrix[current_row]);
		pivot_row[j] = current_row;
		for (int i = 0; i < rows; ++i) {
			if (i != current_row && matrix[i][j] == 1) {
				for (int k = j; k < cols; ++k) matrix[i][k] ^= matrix[current_row][k];
			}
		}
		current_row++;
	}
	vector<vector<int>> solutions;
	for (int j = 0; j < cols; ++j) {
		if (pivot_row[j] == -1) {
			vector<int> sol(cols, 0);
			sol[j] = 1;
			for (int k = 0; k < j; ++k) {
				if (pivot_row[k] != -1 && matrix[pivot_row[k]][j] == 1) sol[k] = 1;
			}
			solutions.push_back(sol);
		}
	}
	return solutions;
}

void factorsFromSystem(const BigInt& n, const vector<BigInt>& b_values, const vector<vector<int>>& degree_vectors, const vector<int>& base) {
	vector<vector<int>> matrix = buildGF2Matrix(degree_vectors);
	vector<vector<int>> solutions = solveGaussianGF2(matrix);
	cout << "Found " << solutions.size() << " potential combinations" << endl;
	for (const auto& sol : solutions) {
		vector<BigInt> selected_b;
		vector<vector<int>> selected_degrees;
		for (size_t i = 0; i < sol.size(); ++i) {
			if (sol[i] == 1) {
				selected_b.push_back(b_values[i]);
				selected_degrees.push_back(degree_vectors[i]);
			}
		}
		solution(n, selected_b, selected_degrees, base);
	}
}

void BrillhartMorrison(const BigInt& n) {
	// Step 1
	vector<int> base = factorBase(n);
	size_t k = base.size();
	vector<BigInt> b_values;
	vector<vector<int>> degree_vectors;
	// Step 2
	BigInt v = 1;
	BigInt u = n.sqrt();
	BigInt a = u;
	BigInt b_prev2 = 0;
	BigInt b_prev1 = 1;
	while (degree_vectors.size() < k + 1) {
		BigInt b_current = (a * b_prev1 + b_prev2) % n;
		BigInt a_val = (b_current * b_current) % n;
		if (a_val > n / BigInt(2)) a_val = a_val - n;
		// Step 3
		vector<int> degrees;
		if (tryFactorize(a_val, base, degrees)) {
			if (zeroVector(degrees)) {
				// Step 4 and 5
				solution(n, { b_current }, { degrees }, base);
				return;
			}
			b_values.push_back(b_current);
			degree_vectors.push_back(degrees);
			cout << "Found B-smooth number #" << degree_vectors.size() << endl;
		}
		BigInt v_next = (n - u * v) / v;
		BigInt a_next = (n.sqrt() + u) / v_next;
		BigInt u_next = a_next * v_next - u;
		b_prev2 = b_prev1;
		b_prev1 = b_current;
		v = v_next;
		u = u_next;
		a = a_next;
	}
	factorsFromSystem(n, b_values, degree_vectors, base);
}

int main()
{
}