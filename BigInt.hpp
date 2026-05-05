#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <random>

using namespace std;
using namespace std::chrono;

class BigInt {
private:
    // 2048 = 64 блоки по 32 біт
    vector<uint32_t> digits;
    static const uint64_t BASE = 0x100000000ULL; // 2^32

    void remove_leading_zeros() {
        while (digits.size() > 1 && digits.back() == 0) {
            digits.pop_back();
        }
    }

    uint32_t divInt(uint32_t v) {
        uint64_t rem = 0;
        for (int i = static_cast<int>(digits.size()) - 1; i >= 0; --i) {
            uint64_t cur = digits[i] + rem * BASE;
            digits[i] = static_cast<uint32_t>(cur / v);
            rem = cur % v;
        }
        remove_leading_zeros();
        return static_cast<uint32_t>(rem);
    }

    void multiplyInt(uint64_t v) {
        if (v == 0) { *this = 0; return; }
        if (v == 1) return;

        uint64_t carry = 0;
        for (size_t i = 0; i < digits.size() || carry; ++i) {
            if (i == digits.size()) digits.push_back(0);
            uint64_t cur = carry + digits[i] * 1ULL * v;
            digits[i] = static_cast<uint32_t>(cur & 0xFFFFFFFF);
            carry = cur >> 32;
        }
        remove_leading_zeros();
    }

    void addInt(uint32_t v) {
        uint64_t carry = v;
        for (size_t i = 0; i < digits.size() || carry; ++i) {
            if (i == digits.size()) digits.push_back(0);
            uint64_t cur = digits[i] * 1ULL + carry;
            digits[i] = static_cast<uint32_t>(cur & 0xFFFFFFFF);
            carry = cur >> 32;
        }
    }

public:
    BigInt() { digits.push_back(0); }
    BigInt(uint32_t n) {
        digits.clear();
        digits.push_back(n);
    }
    size_t getDigitsCount() const {
        return digits.size();
    }
    uint32_t getDigit(size_t i) const {
        if (i >= digits.size()) return 0;
        return digits[i];
    }
    BigInt& operator=(uint32_t n) {
        digits.clear();
        digits.push_back(n);
        return *this;
    }

    bool operator==(const BigInt& other) const {
        return digits == other.digits;
    }

    bool operator!=(const BigInt& other) const {
        return !(*this == other);
    }

    bool operator>(const BigInt& other) const {
        if (digits.size() != other.digits.size()) {
            return digits.size() > other.digits.size();
        }
        for (int i = (int)digits.size() - 1; i >= 0; i--) {
            if (digits[i] != other.digits[i]) {
                return digits[i] > other.digits[i];
            }
        }
        return false;
    }

    bool operator<(const BigInt& other) const {
        return other > *this;
    }

    bool operator<=(const BigInt& other) const {
        return !(*this > other);
    }

    bool operator>=(const BigInt& other) const {
        return !(*this < other);
    }

    BigInt operator+(const BigInt& other) const {
        BigInt res;
        res.digits.clear();

        uint64_t carry = 0;
        size_t n = max(digits.size(), other.digits.size());

        for (size_t i = 0; i < n || carry; ++i) {
            uint64_t sum = carry +
                (i < digits.size() ? digits[i] : 0) +
                (i < other.digits.size() ? other.digits[i] : 0);

            res.digits.push_back(static_cast<uint32_t>(sum & 0xFFFFFFFF));
            carry = sum >> 32;
        }
        return res;
    }

    BigInt operator-(const BigInt& other) const {
        if (*this < other) {
            throw runtime_error("Underflow error: Result of subtraction would be negative.");
        }

        BigInt res = *this;
        uint64_t borrow = 0;

        for (size_t i = 0; i < other.digits.size() || borrow; ++i) {
            uint64_t sub = (i < other.digits.size() ? (uint64_t)other.digits[i] : 0) + borrow;

            if ((uint64_t)res.digits[i] < sub) {
                res.digits[i] = static_cast<uint32_t>(0x100000000ULL + res.digits[i] - sub);
                borrow = 1;
            }
            else {
                res.digits[i] = static_cast<uint32_t>(res.digits[i] - sub);
                borrow = 0;
            }
        }
        res.remove_leading_zeros();
        return res;
    }

    BigInt operator*(const BigInt& other) const {
        if (*this == 0 || other == 0) return BigInt(0);

        BigInt res;
        res.digits.assign(digits.size() + other.digits.size(), 0);

        for (size_t i = 0; i < digits.size(); ++i) {
            if (digits[i] == 0) continue;
            uint64_t carry = 0;
            for (size_t j = 0; j < other.digits.size() || carry; ++j) {
                uint64_t cur = res.digits[i + j] +
                    (uint64_t)digits[i] * (j < other.digits.size() ? other.digits[j] : 0) +
                    carry;

                res.digits[i + j] = static_cast<uint32_t>(cur & 0xFFFFFFFFULL);
                carry = cur >> 32;
            }
        }
        res.remove_leading_zeros();
        return res;
    }

    static pair<BigInt, BigInt> divmod(BigInt a, BigInt b) {
        if (b == 0) throw runtime_error("Division by zero");
        if (a < b) return { BigInt(0), a };

        BigInt quotient, remainder;
        quotient.digits.assign(a.digits.size(), 0);

        for (int i = (int)a.digits.size() - 1; i >= 0; i--) {
            remainder.multiplyInt(0x100000000ULL);
            remainder.addInt(a.digits[i]);

            uint64_t low = 0, high = 0xFFFFFFFFULL;
            uint32_t q = 0;

            while (low <= high) {
                uint64_t mid = low + (high - low) / 2;
                BigInt temp = b;
                temp.multiplyInt(mid);

                if (temp <= remainder) {
                    q = (uint32_t)mid;
                    low = mid + 1;
                }
                else {
                    high = mid - 1;
                }
            }

            quotient.digits[i] = q;
            BigInt subtractor = b;
            subtractor.multiplyInt(q);
            remainder = remainder - subtractor;
        }

        quotient.remove_leading_zeros();
        remainder.remove_leading_zeros();
        return { quotient, remainder };
    }

    BigInt operator/(const BigInt& other) const {
        return divmod(*this, other).first;
    }

    BigInt operator%(const BigInt& other) const {
        return divmod(*this, other).second;
    }

    BigInt square() const {
        size_t n = digits.size();

        if (n < 4) {
            BigInt res;
            res.digits.assign(2 * n, 0);

            for (size_t i = 0; i < n; i++) {
                uint64_t carry = 0;
                for (size_t j = i + 1; j < n; j++) {
                    uint64_t prod = (uint64_t)digits[i] * digits[j];
                    uint64_t sum = (uint64_t)res.digits[i + j] + prod + carry;
                    res.digits[i + j] = static_cast<uint32_t>(sum & 0xFFFFFFFF);
                    carry = sum >> 32;
                }
                res.digits[i + n] = static_cast<uint32_t>(carry);
            }

            uint32_t shiftCarry = 0;
            for (size_t i = 0; i < res.digits.size(); i++) {
                uint64_t cur = ((uint64_t)res.digits[i] << 1) | shiftCarry;
                res.digits[i] = static_cast<uint32_t>(cur & 0xFFFFFFFF);
                shiftCarry = static_cast<uint32_t>(cur >> 32);
            }

            uint64_t addCarry = 0;
            for (size_t i = 0; i < n; i++) {
                uint64_t sq = (uint64_t)digits[i] * digits[i];
                uint64_t sum1 = (uint64_t)res.digits[2 * i] + (sq & 0xFFFFFFFF) + addCarry;
                res.digits[2 * i] = static_cast<uint32_t>(sum1 & 0xFFFFFFFF);
                uint64_t sum2 = (uint64_t)res.digits[2 * i + 1] + (sq >> 32) + (sum1 >> 32);
                res.digits[2 * i + 1] = static_cast<uint32_t>(sum2 & 0xFFFFFFFF);
                addCarry = sum2 >> 32;
            }
            res.remove_leading_zeros();
            return res;
        }

        BigInt a_copy = *this;
        if (n % 2 != 0) {
            a_copy.digits.push_back(0);
            n++;
        }
        size_t k = n / 2;

        // a = a1 + x^k * a2
        BigInt a1, a2;
        a1.digits.assign(a_copy.digits.begin(), a_copy.digits.begin() + k);
        a2.digits.assign(a_copy.digits.begin() + k, a_copy.digits.end());

        a1.remove_leading_zeros();
        a2.remove_leading_zeros();

        BigInt p1 = a1.square();
        BigInt p2 = a2.square();

        BigInt sum_a = a1 + a2;
        BigInt t = sum_a.square();

        // c = p1 + x^k * (t - p1 - p2) + x^(2k) * p2

        BigInt middle = t - p1 - p2; // (t - p1 - p2)

        // x^k (зсув на k блоків)
        auto shiftBlocks = [](BigInt& val, size_t blocks) {
            if (val.digits.size() == 1 && val.digits[0] == 0) return;
            val.digits.insert(val.digits.begin(), blocks, 0);
            };

        shiftBlocks(middle, k);
        shiftBlocks(p2, 2 * k);

        BigInt res = p1 + middle + p2;
        res.remove_leading_zeros();
        return res;
    }

    size_t bitLength() const {
        if (digits.empty() || (digits.size() == 1 && digits[0] == 0)) return 0;

        size_t bits = (digits.size() - 1) * 32;

        uint32_t lastBlock = digits.back();
        while (lastBlock > 0) {
            lastBlock >>= 1;
            bits++;
        }
        return bits;
    }

    bool testBit(size_t bitIdx) const {
        size_t blockIdx = bitIdx / 32;
        size_t bitInBlock = bitIdx % 32;
        if (blockIdx >= digits.size()) return false;
        return (digits[blockIdx] >> bitInBlock) & 1;
    }

    BigInt Gor(const BigInt& exponent) const {
        if (exponent == 0) return BigInt(1);
        if (*this == 0) return BigInt(0);
        if (*this == 1) return BigInt(1);

        BigInt res(1);
        size_t nBits = exponent.bitLength();

        for (int i = (int)nBits - 1; i >= 0; i--) {
            res = res.square();

            if (exponent.testBit(i)) {
                res = res * (*this);
            }
        }

        return res;
    }

    BigInt killLastDigits(size_t k) const { // = зсув вправо
        if (k >= digits.size()) return BigInt(0);
        BigInt res;
        res.digits.assign(digits.begin() + k, digits.end());
        return res;
    }
    static BigInt calculateMu(const BigInt& n) {
        size_t k = n.digits.size();

        BigInt beta2k;
        beta2k.digits.assign(2 * k + 1, 0);
        beta2k.digits[2 * k] = 1;

        return beta2k / n;
    }

    static BigInt barrettReduction(const BigInt& x, const BigInt& n, const BigInt& mu) {
        size_t k = n.digits.size();

        BigInt q = x.killLastDigits(k - 1);
        q = q * mu;
        q = q.killLastDigits(k + 1);

        BigInt r = x - (q * n);
        while (r > n || r == n) {
            r = r - n;
        }

        return r;
    }

    static BigInt gcd(BigInt a, BigInt b) {
        while (!(b == 0)) {
            a = a % b;
            swap(a, b);
        }
        return a;
    }

    static BigInt lcm(BigInt a, BigInt b) {
        if (a == 0 || b == 0) return BigInt(0);
        return (a / gcd(a, b)) * b;
    }

    BigInt addMod(const BigInt& other, const BigInt& n) const {
        BigInt res = *this + other;
        while (res > n || res == n) {
            res = res - n;
        }
        return res;
    }

    BigInt subMod(const BigInt& other, const BigInt& n) const {
        if (n == 0) throw runtime_error("Division by zero (modulus is 0)");

        BigInt a = *this % n;
        BigInt b = other % n;

        if (a > b || a == b) {
            return (a - b) % n;
        }
        else {
            return n - (b - a);
        }
    }

    BigInt mulMod(const BigInt& other, const BigInt& n) const {
        if (n == 0) throw runtime_error("Division by zero");
        if (n == 1) return BigInt(0);

        BigInt mu = BigInt::calculateMu(n);

        BigInt multiplicationResult = (*this) * other;

        return BigInt::barrettReduction(multiplicationResult, n, mu);
    }

    BigInt squareMod(const BigInt& n) const {
        if (n == 0) throw runtime_error("Division by zero");

        BigInt mu = BigInt::calculateMu(n);
        return BigInt::barrettReduction(this->square(), n, mu);
    }

    BigInt GorMod(const BigInt& d, const BigInt& n) const {
        if (n == BigInt(1)) return BigInt(0);
        if (d == BigInt(0)) return BigInt(1);

        BigInt mu = BigInt::calculateMu(n);

        BigInt res(1);
        BigInt a = *this % n;
        size_t bits = d.bitLength();

        for (int i = (int)bits - 1; i >= 0; i--) {
            res = BigInt::barrettReduction(res.square(), n, mu);

            if (d.testBit(i)) {
                res = BigInt::barrettReduction(res * a, n, mu);
            }
        }
        return res;
    }

    BigInt sqrt() const {
        if (*this <= BigInt(1)) return *this;
        BigInt x = *this / BigInt(2);
        BigInt y = (x + *this / x) / BigInt(2);
        while (y < x) {
            x = y;
            y = (x + *this / x) / BigInt(2);
        }
        return x;
    }

    void fromHexString(string hexStr) {
        digits.clear();
        for (int i = (int)hexStr.length(); i > 0; i -= 8) {
            int start = max(0, i - 8);
            int len = i - start;
            string part = hexStr.substr(start, len);
            digits.push_back(stoul(part, nullptr, 16));
        }
        remove_leading_zeros();
    }

    string toHexString() const {
        if (digits.empty()) return "0";
        stringstream ss;
        ss << hex << digits.back();
        for (int i = static_cast<int>(digits.size()) - 2; i >= 0; --i) {
            ss << setfill('0') << setw(8) << hex << digits[i];
        }
        return ss.str();
    }

    void fromDecimalString(string s) {
        digits.clear();
        digits.push_back(0);

        size_t start = 0;
        while (start < s.length()) {
            size_t len = min((size_t)9, s.length() - start);
            string part = s.substr(start, len);
            uint32_t v = stoul(part);

            uint32_t multiplier = 1;
            for (size_t i = 0; i < len; ++i) multiplier *= 10;

            multiplyInt(multiplier);
            addInt(v);

            start += len;
        }
        remove_leading_zeros();
    }

    string toDecimalString() const {
        if (digits.size() == 1 && digits[0] == 0) return "0";
        BigInt temp = *this;
        string res = "";
        while (!(temp.digits.size() == 1 && temp.digits[0] == 0)) {
            uint32_t rem = temp.divInt(1000000000);
            string s = to_string(rem);
            if (!(temp.digits.size() == 1 && temp.digits[0] == 0)) {
                while (s.length() < 9) s = "0" + s;
            }
            res = s + res;
        }
        return res;
    }

    friend ostream& operator<<(ostream& os, const BigInt& bi) {
        os << bi.toHexString();
        return os;
    }

    /*uint32_t getDigit(size_t i) const {
        if (i >= digits.size()) return 0;
        return digits[i];
    }*/

    static BigInt generateRandom(size_t numBlocks) {
        BigInt res;
        res.digits.resize(numBlocks);

        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<uint32_t> dis(0, 0xFFFFFFFF);

        for (size_t i = 0; i < numBlocks; ++i) {
            res.digits[i] = dis(gen);
        }

        res.remove_leading_zeros();
        return res;
    }

    void timeMeasure() {
        const int iterations = 500;
        const size_t size2048 = 64;

        vector<pair<BigInt, BigInt>> data;
        for (int i = 0; i < iterations; ++i) {
            data.push_back({ BigInt::generateRandom(size2048), BigInt::generateRandom(size2048) });
        }

        auto measure = [&](string name, auto func) {
            auto start = high_resolution_clock::now();
            for (int i = 0; i < iterations; ++i) {
                func(data[i].first, data[i].second);
            }
            auto end = high_resolution_clock::now();

            auto avg = duration_cast<nanoseconds>(end - start).count() / iterations;
            cout << left << setw(20) << name << ": " << avg << " ns" << endl;
            };

        cout << "--- Time measure (500) ---" << endl;

        measure("Addition (+)", [](BigInt& a, BigInt& b) { BigInt r = a + b; });
        measure("Multiplication (*)", [](BigInt& a, BigInt& b) { BigInt r = a * b; });
        measure("Division (/)", [](BigInt& a, BigInt& b) { BigInt r = a / b; });
        /*measure("Exponentiation (pow)", [](BigInt& a, BigInt& b) {
            BigInt exp(b.getDigit(0)); // лише перший блок, інакше занадто довго виконується код
            BigInt r = a.Gor(exp);
            });*/
    }
};

