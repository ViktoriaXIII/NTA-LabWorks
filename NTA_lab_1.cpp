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
					isStornglyPseudoPrime == false;
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



int main()
{
}