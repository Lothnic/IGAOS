// OSQP-style ADMM for convex QP (from mathematical foundations, per
// docs/DEPENDENCIES.md — no solver library):
//
//   minimize    0.5 x'Px + q'x
//   subject to  l <= Ax <= u
//
// Scaled dual ADMM (Stellato et al. 2020, implemented from the paper):
//   [[P + sigma I, A'], [A, -1/rho I]] [x; nu]
//        = [sigma*x^k - q; z^k - y^k/rho]
//   xtil <- x^{k+1} + alpha (x^{k+1} - x^k)
//   z    <- Pi_[l,u](A xtil + y^k/rho)
//   y    <- y^k + rho (A xtil - z)
// Dense KKT factorization via the simplex's DenseLU — correct for
// Haverly-class sizes (m ~ tens); sparse KKT is the upgrade path.
#include "qp.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

#include "dense_lu.hpp"

namespace igaos::qp {

namespace {
constexpr double INF = std::numeric_limits<double>::infinity();

inline double clamp_box(double v, double l, double u) {
    if (v < l) return l;
    if (v > u) return u;
    return v;
}
}  // namespace

Solution solve(const io::Model& model, const Options& opt) {
    Solution sol;
    auto t0 = std::chrono::steady_clock::now();
    auto elapsed = [&]() {
        return std::chrono::duration<double>(
                   std::chrono::steady_clock::now() - t0)
            .count();
    };

    const int n = model.n, m = model.m;

    // Constraint system: model rows PLUS one row per variable bound
    // (variable bounds must live in [l,u] or ADMM optimizes the wrong
    // problem — Haverly without x-bounds is unbounded).
    int nb = 0;
    for (int j = 0; j < n; ++j)
        if (std::isfinite(model.cl[j]) || std::isfinite(model.cu[j])) ++nb;
    const int M = m + nb;

    // dense constraint matrix A (row-major) — fine at this scale
    std::vector<double> A((size_t)M * n, 0.0);
    for (int j = 0; j < n; ++j)
        for (int p = model.cp[j]; p < model.cp[j + 1]; ++p)
            A[(size_t)model.ci[p] * n + j] = model.acx[p];
    {
        int r = m;
        for (int j = 0; j < n; ++j) {
            if (std::isfinite(model.cl[j]) || std::isfinite(model.cu[j])) {
                A[(size_t)r * n + j] = 1.0;
                ++r;
            }
        }
    }
    std::vector<double> l(M), u(M);
    for (int i = 0; i < m; ++i) { l[i] = model.rmin[i]; u[i] = model.rmax[i]; }
    {
        int r = m;
        for (int j = 0; j < n; ++j) {
            if (std::isfinite(model.cl[j]) || std::isfinite(model.cu[j])) {
                l[r] = model.cl[j];
                u[r] = model.cu[j];
                ++r;
            }
        }
    }

    // dense P (off-diagonal QUADOBJ entries doubled)
    std::vector<double> P((size_t)n * n, 0.0);
    for (size_t k = 0; k < model.q_v.size(); ++k) {
        int i = model.q_i[k], j = model.q_j[k];
        P[(size_t)i * n + j] += model.q_v[k];
        if (i != j) P[(size_t)j * n + i] += model.q_v[k];
    }
    std::vector<double> q(n);
    for (int j = 0; j < n; ++j) q[j] = model.c[j];

    const double sigma = 1e-3, rho = 0.1, alpha = 0.0;
    // ponytail: fixed (sigma, rho), no over-relaxation — alpha>0
    // diverges on Haverly-class problems in this implementation;
    // revisit with OSQP §5 adaptive rho if large instances land here

    // KKT = [[P + sigma I, A'], [A, -1/rho I]]  (size n+m)
    const int K = n + M;
    std::vector<double> Kmat((size_t)K * K, 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            Kmat[(size_t)i * K + j] =
                P[(size_t)i * n + j] + (i == j ? sigma : 0.0);
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < n; ++j)
            Kmat[(size_t)(n + i) * K + j] = A[(size_t)i * n + j];
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < n; ++j)
            Kmat[(size_t)j * K + (n + i)] = A[(size_t)i * n + j];
    for (int i = 0; i < M; ++i) Kmat[(size_t)(n + i) * K + (n + i)] = -1.0 / rho;

    simplex::DenseLU lu;
    {
        std::vector<double> tmp = Kmat;
        lu.factor(std::move(tmp));
    }
    if (!lu.ok) {
        sol.status = Status::Error;
        sol.message = "KKT factorization failed";
        sol.solve_time_ms = elapsed() * 1000.0;
        return sol;
    }

    std::vector<double> x(n, 0.0), z(M, 0.0), yv(M, 0.0);
    std::vector<double> rhs(K), solk(K), xt(n), Ax(M), Aty(n);

    auto matvec_A = [&]() {
        for (int i = 0; i < M; ++i) {
            double s = 0.0;
            for (int j = 0; j < n; ++j)
                if (A[(size_t)i * n + j] != 0.0)
                    s += A[(size_t)i * n + j] * x[j];
            Ax[i] = s;
        }
    };
    auto objective = [&]() {
        double obj = model.obj_const;
        for (int j = 0; j < n; ++j) obj += q[j] * x[j];
        for (int k = 0; k < n; ++k)
            for (int j = 0; j < n; ++j)
                if (P[(size_t)k * n + j] != 0.0)
                    obj += 0.5 * P[(size_t)k * n + j] * x[k] * x[j];
        return obj;
    };

    // data-scale norms for the (absolute) tolerance ladder — never scaled
    // by the running objective (diverging iterates would self-pass)
    double cmax = 0.0, Pmax = 0.0, Amax = 0.0;
    for (double v : q) cmax = std::max(cmax, std::fabs(v));
    for (double v : P) Pmax = std::max(Pmax, std::fabs(v));
    for (double v : A) Amax = std::max(Amax, std::fabs(v));

    long it = 0;
    double pinf = INF, dinf = INF;
    for (; it < opt.max_iterations; ++it) {
        if (elapsed() > opt.time_limit_s) {
            sol.status = Status::TimeLimit;
            sol.message = "wall-clock budget exhausted";
            break;
        }
        std::vector<double> xprev = x;
        for (int j = 0; j < n; ++j) rhs[j] = sigma * xprev[j] - q[j];
        for (int i = 0; i < M; ++i) rhs[n + i] = z[i] - yv[i] / rho;
        lu.solve(rhs, solk);
        for (int j = 0; j < n; ++j) x[j] = solk[j];

        for (int j = 0; j < n; ++j) xt[j] = x[j] + alpha * (x[j] - xprev[j]);
        for (int i = 0; i < M; ++i) {
            double s = 0.0;
            for (int j = 0; j < n; ++j)
                if (A[(size_t)i * n + j] != 0.0)
                    s += A[(size_t)i * n + j] * xt[j];
            Ax[i] = s;
        }
        for (int i = 0; i < M; ++i)
            z[i] = clamp_box(Ax[i] + yv[i] / rho, l[i], u[i]);
        for (int i = 0; i < M; ++i) yv[i] += rho * (Ax[i] - z[i]);

        if (it % 25 == 24 || it == 0) {
            matvec_A();
            pinf = 0.0;
            for (int i = 0; i < M; ++i) {
                double viol = std::max(l[i] - Ax[i], Ax[i] - u[i]);
                pinf = std::max(pinf, viol);
            }
            for (int j = 0; j < n; ++j) Aty[j] = 0.0;
            for (int i = 0; i < M; ++i)
                if (yv[i] != 0.0)
                    for (int j = 0; j < n; ++j)
                        if (A[(size_t)i * n + j] != 0.0)
                            Aty[j] += A[(size_t)i * n + j] * yv[i];
            dinf = 0.0;
            for (int j = 0; j < n; ++j) {
                double dr = q[j] + Aty[j];
                for (int k = 0; k < n; ++k)
                    if (P[(size_t)j * n + k] != 0.0)
                        dr += P[(size_t)j * n + k] * x[k];
                dinf = std::max(dinf, std::fabs(dr));
            }
            double scale = 1.0 + cmax + Pmax + Amax;
            if (pinf <= opt.tolerance * scale && dinf <= opt.tolerance * scale) {
                sol.status = Status::Optimal;
                break;
            }
        }
    }
    if (sol.status != Status::Optimal && sol.status != Status::TimeLimit &&
        it >= opt.max_iterations) {
        sol.status = Status::IterationLimit;
        sol.message = "iteration cap reached";
    }

    sol.objective = objective();
    sol.x = x;
    sol.y.assign(yv.begin(), yv.begin() + m);
    sol.pinf = pinf;
    sol.dinf = dinf;
    sol.rel_gap = pinf + dinf;
    sol.iterations = it;
    sol.solve_time_ms = elapsed() * 1000.0;
    return sol;
}

}  // namespace igaos::qp
