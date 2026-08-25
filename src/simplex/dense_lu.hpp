#pragma once

#include <cmath>
#include <vector>

namespace igaos::simplex {

struct DenseLU {
    int n = 0;
    std::vector<double> a;
    std::vector<double> kept;
    std::vector<int> perm;
    bool ok = false;

    void factor(std::vector<double> m_in) {
        n = static_cast<int>(std::sqrt((double)m_in.size()));
        kept = m_in;
        a = std::move(m_in);
        perm.resize(n);
        ok = true;
        for (int i = 0; i < n; ++i) perm[i] = i;
        for (int k = 0; k < n; ++k) {
            int piv = k;
            double best = std::fabs(a[k * n + k]);
            for (int i = k + 1; i < n; ++i) {
                double v = std::fabs(a[i * n + k]);
                if (v > best) {
                    best = v;
                    piv = i;
                }
            }
            if (!(best > 0.0)) {
                ok = false;
                return;
            }
            if (piv != k) {
                for (int j = 0; j < n; ++j)
                    std::swap(a[k * n + j], a[piv * n + j]);
                std::swap(perm[k], perm[piv]);
            }
            double d = a[k * n + k];
            for (int i = k + 1; i < n; ++i) {
                double f = a[i * n + k] / d;
                a[i * n + k] = f;
                if (f == 0.0) continue;
                for (int j = k + 1; j < n; ++j)
                    a[i * n + j] -= f * a[k * n + j];
            }
        }
    }

    void solve(const std::vector<double>& b, std::vector<double>& x) const {
        x.assign(n, 0.0);
        for (int i = 0; i < n; ++i) x[i] = b[perm[i]];
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < i; ++j) x[i] -= a[i * n + j] * x[j];
        for (int i = n - 1; i >= 0; --i) {
            for (int j = i + 1; j < n; ++j) x[i] -= a[i * n + j] * x[j];
            x[i] /= a[i * n + i];
        }
    }

    void solve_transpose(const std::vector<double>& b,
                         std::vector<double>& x) const {
        x.assign(n, 0.0);
        for (int i = 0; i < n; ++i) {
            double s = b[i];
            for (int j = 0; j < i; ++j) s -= a[j * n + i] * x[j];
            x[i] = s / a[i * n + i];
        }
        for (int i = n - 1; i >= 0; --i) {
            double s = x[i];
            for (int j = i + 1; j < n; ++j) s -= a[j * n + i] * x[j];
            x[i] = s;
        }
        std::vector<double> tmp = x;
        for (int i = 0; i < n; ++i) x[perm[i]] = tmp[i];
    }
};

}  // namespace igaos::simplex
