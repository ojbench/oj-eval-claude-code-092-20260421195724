#ifndef FRACTION_HPP
#define FRACTION_HPP

#include <iostream>
#include <numeric>

class fraction {
    long long num, den;
    void simplify() {
        if (den == 0) return;
        if (den < 0) { num = -num; den = -den; }
        long long common = std::gcd(std::abs(num), den);
        num /= common;
        den /= common;
    }

public:
    fraction(long long n = 0, long long d = 1) : num(n), den(d) { simplify(); }

    fraction operator+(const fraction& other) const {
        return fraction(num * other.den + other.num * den, den * other.den);
    }
    fraction operator-(const fraction& other) const {
        return fraction(num * other.den - other.num * den, den * other.den);
    }
    fraction operator*(const fraction& other) const {
        return fraction(num * other.num, den * other.den);
    }
    fraction operator/(const fraction& other) const {
        return fraction(num * other.den, den * other.num);
    }
    fraction& operator+=(const fraction& other) { *this = *this + other; return *this; }
    fraction& operator-=(const fraction& other) { *this = *this - other; return *this; }
    fraction& operator*=(const fraction& other) { *this = *this * other; return *this; }
    fraction& operator/=(const fraction& other) { *this = *this / other; return *this; }

    bool operator==(const fraction& other) const { return num == other.num && den == other.den; }
    bool operator!=(const fraction& other) const { return !(*this == other); }

    friend std::ostream& operator<<(std::ostream& os, const fraction& f) {
        if (f.den == 1) return os << f.num;
        return os << f.num << "/" << f.den;
    }
};

#endif
