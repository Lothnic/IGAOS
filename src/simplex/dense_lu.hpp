#pragma once

#include <algorithm>
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
        // Near-singularity floor: a basis whose LU pivot falls below
        // 1e-8 x (median column magnitude) is numerically singular — solves
        // against it amplify round-off by ~1/floor into the duals (observed:
        // reduced costs of 1e18 driving an eternal pivot cycle on scsd1).
        // The MEDIAN column max (not the global max) calibrates the floor:
        // a single wildly-scaled basis column (wood1p carries entries of
        // 5e5) must not inflate the floor and reject usable bases. Rejecting
        // here lets the caller roll back the offending pivot instead of
        // pivoting on noise forever.
        constexpr double PIVOT_FLOOR = 1e-8;
        std::vector<double> colmax(n, 0.0);
        for (int c = 0; c < n; ++c)
            for (int r = 0; r < n; ++r)
                colmax[c] = std::max(colmax[c],
                                     std::fabs(a[(size_t)r * n + c]));
        auto mid = colmax.begin() + n / 2;
        std::nth_element(colmax.begin(), mid, colmax.end());
        double med = n > 0 ? *mid : 0.0;
        const double pivot_floor = PIVOT_FLOOR * med;
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
            if (!(best > pivot_floor)) {
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

// Product-form eta file: B_k = B_0 · E_1 ··· E_k where E_i is the
// identity with column p_i replaced by alpha_i (= B_{i-1}^{-1} a_enter).
// Refactorize from scratch only every REFACTOR_ETA pivots (design sheet
// #5); between refactorizations solves apply the eta transformations.
struct EtaFile {
    DenseLU base;
    std::vector<int> piv;
    std::vector<std::vector<double>> etas;
    int m = 0;
    bool ok = false;
    bool unstable = false;  // an eta pivot was relatively tiny
    int since_unstable = 0;  // pivots since the unstable flag was raised

    static constexpr int REFACTOR_ETA = 50;
    // Backoff: degenerate problems route pivots through Bland/last-resort
    // fallbacks whose pivots sit below the stability threshold, and
    // refactoring on every one stormed (modszk1: 430 O(m^3) factors in 15s
    // at m=687). Tolerate a few unstable etas before paying for a refactor;
    // the honesty guards (refactor before declaring optimal or unbounded)
    // bound the damage of the extra drift.
    static constexpr int UNSTABLE_BACKOFF = 8;

    void factor(std::vector<double> Bm) {
        m = static_cast<int>(std::sqrt((double)Bm.size()));
        piv.clear();
        etas.clear();
        unstable = false;
        since_unstable = 0;
        base.factor(std::move(Bm));
        ok = base.ok;
    }

    // register the pivot replacing basis column p with a_enter whose
    // tableau column (B^{-1} a_enter) is alpha
    void update(int p, const std::vector<double>& alpha) {
        piv.push_back(p);
        etas.push_back(alpha);
        ++since_unstable;
        // small relative pivot elements wreck product-form accuracy —
        // flag for refactorization at the next opportunity
        double amax = 0.0;
        for (double v : alpha) amax = std::max(amax, std::fabs(v));
        // 1e-4 relative: the 1e-2 trigger refactorized on nearly every
        // pivot of ill-conditioned Netlib rows (perold: 505 FACTs/1670
        // pivots — the O(m^3) factorization dominated the solve). Etas with
        // pivot >= 1e-4*amax lose ~4 digits worst-case, which the honesty
        // guards (refactor-before-declare, PIVOT_ABS floor) already cover.
        if (std::fabs(alpha[p]) < 1e-4 * amax) unstable = true;
        ok = base.ok;
    }

    bool needs_refactor() const {
        return (int)etas.size() >= REFACTOR_ETA ||
               (unstable && since_unstable >= UNSTABLE_BACKOFF);
    }

    // x = B^{-1} v: base solve, then etas in insertion order
    void solve(const std::vector<double>& v, std::vector<double>& x) const {
        base.solve(v, x);
        for (size_t k = 0; k < etas.size(); ++k) {
            const std::vector<double>& a = etas[k];
            int p = piv[k];
            double ap = a[p];
            double xp = x[p] / ap;
            if (xp != 0.0)
                for (int i = 0; i < m; ++i)
                    if (i != p) x[i] -= a[i] * xp;
            x[p] = xp;
        }
    }

    // x = B^{-T} v: etas in reverse order, then base transpose solve
    void solve_transpose(const std::vector<double>& v,
                         std::vector<double>& x) const {
        x = v;
        for (size_t kk = etas.size(); kk-- > 0;) {
            const std::vector<double>& a = etas[kk];
            int p = piv[kk];
            double s = x[p];
            for (int i = 0; i < m; ++i)
                if (i != p) s -= a[i] * x[i];
            x[p] = s / a[p];
        }
        std::vector<double> tmp;
        base.solve_transpose(x, tmp);
        x = tmp;
    }
};

}  // namespace igaos::simplex
