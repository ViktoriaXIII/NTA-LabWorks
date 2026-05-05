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

int main()
{
}