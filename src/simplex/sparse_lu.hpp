#pragma once

// Sparse LU for simplex bases: Gilbert-Peierls-style left-looking
// factorization with column ordering by original column count (static
// minimum-degree surrogate) and threshold partial pivoting on the spike
// (max |entry| among active rows). All triangular solves run in
// O(m + nnz(L,U)) via column scatter/gather — no DFS at solve time.
//
// ponytail: column ordering is static (original column count). A Markowitz
// dynamic ordering (min row-count x col-count over the active submatrix) is
// the upgrade path if fill becomes the bottleneck; measured fill on Netlib
// bases is low enough that the simpler order wins on factor time.

#include <algorithm>
#include <cmath>
#include <vector>

#include "dense_lu.hpp"

namespace igaos::simplex {

struct SparseLU {
    int n = 0;
    bool ok = false;
    bool dense_mode = false;  // dense fallback active (factor failed sparse)

    // Factors in step space (permuted). L: unit diagonal implicit, stored
    // column-major; Lcol[t] holds (row_step, val) pairs with row_step > t.
    // U: stored column-major; Ucol[j] holds (row_step <= j, val) incl. diag.
    std::vector<int> Lstart, Lrow;   // CSC over steps
    std::vector<double> Lval;
    std::vector<int> Ustart, Urow;
    std::vector<double> Uval;
    std::vector<int> prow;    // step -> original row (pivot row)
    std::vector<int> rowstep; // original row -> step (inverse)

    DenseLU dn;  // dense fallback / IGAOS_DENSE_LU escape hatch

    // CSC input: column j of the basis matrix is
    // (rowind[cp[j]..cp[j+1]), vals[same)). n x n.
    void factor_csc(int n_, const std::vector<int>& cp,
                    const std::vector<int>& rowind,
                    const std::vector<double>& vals,
                    bool allow_dense_fallback = true) {
        n = n_;
        dense_mode = false;
        ok = factor_sparse(cp, rowind, vals);
        if (!ok && allow_dense_fallback && n <= 3000) {
            // sparse rejected the basis: cross-check dense before declaring
            // singular — a disagreement here routes the simplex into basis
            // rollback, so only trust the rejection if dense agrees
            std::vector<double> Bm((size_t)n * n, 0.0);
            for (int j = 0; j < n; ++j)
                for (int p = cp[j]; p < cp[j + 1]; ++p)
                    Bm[(size_t)rowind[p] * n + j] = vals[p];
            dn.factor(std::move(Bm));
            if (dn.ok) dense_mode = true;
        }
    }

    void factor_dense(std::vector<double> Bm) {
        n = (int)std::sqrt((double)Bm.size());
        dense_mode = true;
        dn.factor(std::move(Bm));
        ok = dn.ok;
    }

    double pivot_mag(int t) const {
        // |U(t,t)| — diagnostics (FACT trace)
        if (dense_mode) return std::fabs(dn.a[(size_t)t * n + t]);
        for (int p = Ustart[t]; p < Ustart[t + 1]; ++p)
            if (Urow[p] == t) return std::fabs(Uval[p]);
        return 0.0;
    }

    void solve(const std::vector<double>& b, std::vector<double>& x) const {
        if (dense_mode) { dn.solve(b, x); return; }
        // P B Q = L U.  x = B^{-1} b: L z = P b, U w = z, x = Q w.
        x.assign(n, 0.0);
        for (int t = 0; t < n; ++t) x[t] = b[prow[t]];
        for (int t = 0; t < n; ++t) {          // L z = z: column scatter
            double zt = x[t];
            if (zt == 0.0) continue;
            for (int p = Lstart[t]; p < Lstart[t + 1]; ++p)
                x[Lrow[p]] -= Lval[p] * zt;
        }
        for (int j = n - 1; j >= 0; --j) {     // U w = z: column scatter
            double zj = x[j];
            if (zj == 0.0) continue;
            double d = 0.0;
            for (int p = Ustart[j]; p < Ustart[j + 1]; ++p)
                if (Urow[p] == j) d = Uval[p];
            double wj = zj / d;
            x[j] = wj;
            if (wj == 0.0) continue;
            for (int p = Ustart[j]; p < Ustart[j + 1]; ++p)
                if (Urow[p] != j) x[Urow[p]] -= Uval[p] * wj;
        }
        // w(j) lives in step space; the answer is indexed by original
        // column: x[qcol[j]] = w(j)
        std::vector<double> w = x;
        for (int j = 0; j < n; ++j) x[qcol[j]] = w[j];
    }

    void solve_transpose(const std::vector<double>& b,
                         std::vector<double>& x) const {
        if (dense_mode) { dn.solve_transpose(b, x); return; }
        // B^T = Q U^T L^T P.  x = B^{-T} b:
        //   U^T s = Q^T b,  L^T u = s,  x = P^T u.
        solve_transpose_impl(b, x);
    }

    // Pivot floor mirrors DenseLU: 1e-8 x median column max of the basis.
    double pivot_floor_from(const std::vector<int>& cp,
                            const std::vector<double>& vals) const {
        std::vector<double> colmax(n, 0.0);
        for (int j = 0; j < n; ++j) {
            double cm = 0.0;
            for (int p = cp[j]; p < cp[j + 1]; ++p)
                cm = std::max(cm, std::fabs(vals[p]));
            colmax[j] = cm;
        }
        auto mid = colmax.begin() + n / 2;
        std::nth_element(colmax.begin(), mid, colmax.end());
        return 1e-8 * (n > 0 ? *mid : 0.0);
    }

private:
    void solve_transpose_impl(const std::vector<double>& b,
                              std::vector<double>& x) const {
        // v(j) = b[q_j]: the column permutation q is baked into the CSC
        // input order (step j eliminates original column qcol[j]).
        std::vector<double> s(n), u(n);
        // U^T s = v (lower solve, gather over U columns)
        for (int j = 0; j < n; ++j) {
            double acc = b[qcol[j]];
            for (int p = Ustart[j]; p < Ustart[j + 1]; ++p)
                if (Urow[p] < j) acc -= Uval[p] * s[Urow[p]];
            double d = 0.0;
            for (int p = Ustart[j]; p < Ustart[j + 1]; ++p)
                if (Urow[p] == j) d = Uval[p];
            s[j] = acc / d;
        }
        // L^T u = s (unit upper solve, gather over L columns)
        for (int t = n - 1; t >= 0; --t) {
            double acc = s[t];
            for (int p = Lstart[t]; p < Lstart[t + 1]; ++p)
                acc -= Lval[p] * u[Lrow[p]];
            u[t] = acc;
        }
        x.assign(n, 0.0);
        for (int t = 0; t < n; ++t) x[prow[t]] = u[t];
    }

    std::vector<int> qcol;  // step -> original column (the column permutation)

    bool factor_sparse(const std::vector<int>& cp_in,
                       const std::vector<int>& rowind_in,
                       const std::vector<double>& vals_in) {
        // 1. column order: increasing original column count (ties: index)
        std::vector<int> cnt(n);
        for (int j = 0; j < n; ++j) cnt[j] = cp_in[j + 1] - cp_in[j];
        qcol.resize(n);
        for (int j = 0; j < n; ++j) qcol[j] = j;
        std::stable_sort(qcol.begin(), qcol.end(),
                         [&](int a, int b) { return cnt[a] < cnt[b]; });

        const double floor_piv = pivot_floor_from(cp_in, vals_in);

        Lstart.assign(n + 1, 0);
        Ustart.assign(n + 1, 0);
        Lrow.clear(); Lval.clear();
        Urow.clear(); Uval.clear();
        prow.assign(n, -1);
        rowstep.assign(n, -1);

        std::vector<double> x(n, 0.0);       // dense work array (orig rows)
        std::vector<char> touched(n, 0);     // spike candidates
        std::vector<int> cand;               // candidate rows this column
        std::vector<char> reached(n, 0);     // steps reached this column
        std::vector<int> reach;              // reached steps (unsorted)
        std::vector<int> stack;              // BFS stack over steps

        for (int j = 0; j < n; ++j) {
            int col = qcol[j];
            // reset per-column scratch
            for (int r : cand) { touched[r] = 0; x[r] = 0.0; }
            cand.clear();
            for (int t : reach) reached[t] = 0;
            reach.clear();

            // scatter original column
            for (int p = cp_in[col]; p < cp_in[col + 1]; ++p) {
                int r = rowind_in[p];
                if (touched[r] == 0) { touched[r] = 1; cand.push_back(r); }
                x[r] = vals_in[p];
            }
            // BFS over L columns from the assigned rows in the pattern.
            // Edge s -> t exists iff column s of L (original-row space)
            // contains row prow[t]; unassigned rows are terminal.
            stack.clear();
            for (int p = cp_in[col]; p < cp_in[col + 1]; ++p) {
                int t = rowstep[rowind_in[p]];
                if (t >= 0 && !reached[t]) {
                    reached[t] = 1;
                    reach.push_back(t);
                    stack.push_back(t);
                }
            }
            while (!stack.empty()) {
                int s = stack.back(); stack.pop_back();
                for (int p = Lstart[s]; p < Lstart[s + 1]; ++p) {
                    int t = rowstep[Lrow[p]];
                    if (t >= 0 && !reached[t]) {
                        reached[t] = 1;
                        reach.push_back(t);
                        stack.push_back(t);
                    }
                }
            }
            // numeric: ascending step order is topological (all L edges go
            // s -> t with t > s)
            std::sort(reach.begin(), reach.end());
            for (int s : reach) {
                double v = x[prow[s]];
                if (v == 0.0) continue;
                Urow.push_back(s); Uval.push_back(v);  // U(s, j)
                for (int p = Lstart[s]; p < Lstart[s + 1]; ++p) {
                    int r = Lrow[p];
                    x[r] -= Lval[p] * v;
                    if (touched[r] == 0) { touched[r] = 1; cand.push_back(r); }
                }
            }
            // pivot: max |x[r]| among unassigned candidate rows
            int piv_row = -1;
            double best = 0.0;
            for (int r : cand) {
                if (rowstep[r] >= 0) continue;
                double v = std::fabs(x[r]);
                if (v > best) { best = v; piv_row = r; }
            }
            if (piv_row < 0 || !(best > floor_piv)) return false;
            prow[j] = piv_row;
            rowstep[piv_row] = j;
            Urow.push_back(j); Uval.push_back(x[piv_row]);  // U(j, j)
            Ustart[j + 1] = (int)Urow.size();
            double d = x[piv_row];
            for (int r : cand) {
                if (rowstep[r] >= 0 || r == piv_row) continue;
                if (x[r] == 0.0) continue;
                Lrow.push_back(r); Lval.push_back(x[r] / d);  // orig-row space
            }
            Lstart[j + 1] = (int)Lrow.size();
        }
        // convert L rows to step space for the solves
        for (auto& r : Lrow) r = rowstep[r];
        return true;
    }
};

}  // namespace igaos::simplex
