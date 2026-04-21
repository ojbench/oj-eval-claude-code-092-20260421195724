#ifndef SRC_HPP
#define SRC_HPP

#include <vector>
#include <iostream>
#include <numeric>

// Fraction class included in src.hpp to avoid missing file issues on OJ
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
    bool operator<(const fraction& other) const { return num * other.den < other.num * den; }

    friend std::ostream& operator<<(std::ostream& os, const fraction& f) {
        if (f.den == 1) return os << f.num;
        return os << f.num << "/" << f.den;
    }
};

template <typename T>
class matrix {
    int r, c;
    std::vector<std::vector<T>> data;
public:
    matrix(int rows, int cols) : r(rows), c(cols), data(rows, std::vector<T>(cols, T(0))) {}
    T& operator()(int i, int j) { return data[i][j]; }
    const T& operator()(int i, int j) const { return data[i][j]; }
    int rows() const { return r; }
    int cols() const { return c; }

    std::vector<T> solve(const std::vector<T>& b) const {
        int n = r;
        matrix A = *this;
        std::vector<T> x = b;
        for (int i = 0; i < n; ++i) {
            int pivot = i;
            for(int j = i + 1; j < n; ++j) {
                if (A(pivot, i) == T(0) && A(j, i) != T(0)) {
                    pivot = j;
                    break;
                }
            }
            if (A(pivot, i) == T(0)) continue;
            std::swap(A.data[i], A.data[pivot]);
            std::swap(x[i], x[pivot]);

            T factor = A(i, i);
            for (int j = i; j < n; ++j) A(i, j) /= factor;
            x[i] /= factor;

            for (int k = 0; k < n; ++k) {
                if (k != i) {
                    T f = A(k, i);
                    if (f == T(0)) continue;
                    for (int j = i; j < n; ++j) A(k, j) -= f * A(i, j);
                    x[k] -= f * x[i];
                }
            }
        }
        return x;
    }
};

class resistive_network {
    int n, m;
    std::vector<std::vector<fraction>> L; // Adjacency/Laplacian info
public:
    resistive_network(int n, int m, const std::vector<int>& from, const std::vector<int>& to, const std::vector<fraction>& r)
        : n(n), m(m), L(n, std::vector<fraction>(n, fraction(0))) {
        for (int i = 0; i < m; ++i) {
            int u = from[i] - 1;
            int v = to[i] - 1;
            fraction g = fraction(1) / r[i];
            L[u][u] += g;
            L[v][v] += g;
            L[u][v] -= g;
            L[v][u] -= g;
        }
    }

    fraction get_equivalent_resistance(int u, int v) {
        if (u == v) return fraction(0);
        u--; v--;
        
        int sz = n - 1;
        matrix<fraction> L_red(sz, sz);
        for(int i=0; i<sz; ++i)
            for(int j=0; j<sz; ++j)
                L_red(i, j) = L[i][j];
        
        std::vector<fraction> I_vec(sz, fraction(0));
        if (u < sz) I_vec[u] = fraction(1);
        if (v < sz) I_vec[v] = fraction(-1);
        
        std::vector<fraction> V = L_red.solve(I_vec);
        fraction vu = (u < sz) ? V[u] : fraction(0);
        fraction vv = (v < sz) ? V[v] : fraction(0);
        
        return vu - vv;
    }

    fraction get_voltage(const std::vector<fraction>& I, int k) {
        k--;
        int sz = n - 1;
        matrix<fraction> L_red(sz, sz);
        for(int i=0; i<sz; ++i)
            for(int j=0; j<sz; ++j)
                L_red(i, j) = L[i][j];
        
        std::vector<fraction> I_vec(sz);
        for(int i=0; i<sz; ++i) I_vec[i] = I[i];
        
        std::vector<fraction> V = L_red.solve(I_vec);
        if (k < sz) return V[k];
        return fraction(0);
    }

    fraction get_power(const std::vector<fraction>& U) {
        fraction power(0);
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                if (L[i][j] != fraction(0)) {
                    fraction g = fraction(0) - L[i][j];
                    fraction du = U[i] - U[j];
                    power += du * du * g;
                }
            }
        }
        return power;
    }
};

#endif
