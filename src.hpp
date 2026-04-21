#ifndef SRC_HPP
#define SRC_HPP

#include "fraction.hpp"
#include <vector>
#include <iostream>

template <typename T>
class matrix {
    int r, c;
    std::vector<std::vector<T>> data;
public:
    matrix(int rows, int cols) : r(rows), c(cols), data(rows, std::vector<T>(cols)) {}
    T& operator()(int i, int j) { return data[i][j]; }
    const T& operator()(int i, int j) const { return data[i][j]; }
    int rows() const { return r; }
    int cols() const { return c; }

    static matrix identity(int n) {
        matrix res(n, n);
        for (int i = 0; i < n; ++i) res(i, i) = T(1);
        return res;
    }

    matrix solve(const std::vector<T>& b) const {
        int n = r;
        matrix A = *this;
        std::vector<T> x = b;
        for (int i = 0; i < n; ++i) {
            int pivot = i;
            while (pivot < n && A(pivot, i) == T(0)) pivot++;
            if (pivot == n) continue;
            std::swap(A.data[i], A.data[pivot]);
            std::swap(x[i], x[pivot]);

            T factor = A(i, i);
            for (int j = i; j < n; ++j) A(i, j) /= factor;
            x[i] /= factor;

            for (int k = 0; k < n; ++k) {
                if (k != i) {
                    T f = A(k, i);
                    for (int j = i; j < n; ++j) A(k, j) -= f * A(i, j);
                    x[k] -= f * x[i];
                }
            }
        }
        matrix res(n, 1);
        for (int i = 0; i < n; ++i) res(i, 0) = x[i];
        return res;
    }
};

class resistive_network {
    int n, m;
    matrix<fraction> L; // Laplacian matrix
public:
    resistive_network(int n, int m, const std::vector<int>& from, const std::vector<int>& to, const std::vector<fraction>& r)
        : n(n), m(m), L(n, n) {
        for (int i = 0; i < m; ++i) {
            int u = from[i] - 1;
            int v = to[i] - 1;
            fraction g = fraction(1) / r[i];
            L(u, u) += g;
            L(v, v) += g;
            L(u, v) -= g;
            L(v, u) -= g;
        }
    }

    fraction get_equivalent_resistance(int u, int v) {
        u--; v--;
        // Reduced Laplacian L_n-1 by removing row/col n-1 (node n)
        // Solve L_reduced * V = I where I_u = 1, I_v = -1, others 0
        // But u_n = 0 is fixed.
        if (u == v) return fraction(0);
        
        int sz = n - 1;
        matrix<fraction> L_red(sz, sz);
        for(int i=0; i<sz; ++i)
            for(int j=0; j<sz; ++j)
                L_red(i, j) = L(i, j);
        
        std::vector<fraction> I_vec(sz, fraction(0));
        if (u < sz) I_vec[u] = fraction(1);
        if (v < sz) I_vec[v] = fraction(-1);
        
        matrix<fraction> V = L_red.solve(I_vec);
        fraction vu = (u < sz) ? V(u, 0) : fraction(0);
        fraction vv = (v < sz) ? V(v, 0) : fraction(0);
        
        return vu - vv;
    }

    fraction get_voltage(const std::vector<fraction>& I, int k) {
        k--;
        int sz = n - 1;
        matrix<fraction> L_red(sz, sz);
        for(int i=0; i<sz; ++i)
            for(int j=0; j<sz; ++j)
                L_red(i, j) = L(i, j);
        
        std::vector<fraction> I_vec(sz);
        for(int i=0; i<sz; ++i) I_vec[i] = I[i];
        
        matrix<fraction> V = L_red.solve(I_vec);
        if (k < sz) return V(k, 0);
        return fraction(0);
    }

    fraction get_power(const std::vector<fraction>& U) {
        fraction power(0);
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                if (L(i, j) != fraction(0)) {
                    fraction g = fraction(0) - L(i, j); // conductance is -L(i,j)
                    fraction du = U[i] - U[j];
                    power += du * du * g;
                }
            }
        }
        return power;
    }
};

#endif
