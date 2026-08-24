#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "igaos/options.h"
#include "igaos/solution.h"
#include "igaos/status.h"
#include "model.hpp"

using igaos::io::Model;

namespace {

void usage() {
    std::printf(
        "igaos — sovereign optimization solver (SIH26119)\n\n"
        "USAGE:\n"
        "  igaos info   <model.mps>\n"
        "  igaos solve  <model.mps> [flags]\n\n"
        "FLAGS (solve):\n"
        "  --time-limit <sec>      wall-clock budget (default 120)\n"
        "  --tol <float>           termination tolerance ladder rung "
        "(default 1e-4)\n"
        "  --max-iter <int>        iteration cap (default 500000)\n"
        "  --presolve on|off       presolve pass toggle (default on)\n"
        "  --seed <int>            perturbation seed (default 0)\n"
        "  --out <file.json>       write solution JSON here (else stdout)\n"
        "  --verbose               engine chatter to stderr\n");
}

struct Args {
    std::string cmd, model_path, out_path;
    igaos::Options opts;
    bool ok = true;
};

Args parse_args(int argc, char** argv) {
    Args a;
    if (argc < 3) {
        a.ok = false;
        return a;
    }
    a.cmd = argv[1];
    a.model_path = argv[2];
    if (a.cmd != "info" && a.cmd != "solve") a.ok = false;
    for (int i = 3; i < argc && a.ok; ++i) {
        std::string f = argv[i];
        auto need = [&](double& dst) {
            if (i + 1 >= argc) {
                a.ok = false;
                return;
            }
            dst = std::atof(argv[++i]);
        };
        auto needl = [&](long& dst) {
            if (i + 1 >= argc) {
                a.ok = false;
                return;
            }
            dst = std::atol(argv[++i]);
        };
        if (f == "--time-limit") need(a.opts.time_limit_s);
        else if (f == "--tol") need(a.opts.tolerance);
        else if (f == "--max-iter") needl(a.opts.max_iterations);
        else if (f == "--seed") needl(reinterpret_cast<long&>(a.opts.seed));
        else if (f == "--presolve" && i + 1 < argc)
            a.opts.presolve = (std::string(argv[++i]) == "on");
        else if (f == "--verbose") a.opts.verbosity = 1;
        else if (f == "--out" && i + 1 < argc) a.out_path = argv[++i];
        else a.ok = false;
    }
    return a;
}

std::string json_escape(const std::string& s) {
    std::string r;
    for (char ch : s) {
        if (ch == '"' || ch == '\\') r += '\\';
        r += ch;
    }
    return r;
}

void emit_json(const Model& mdl, const igaos::Solution& sol,
               const std::string& instance, const std::string& out_path) {
    char buf[512];
    std::snprintf(buf, sizeof(buf),
                  "{\n  \"instance\": \"%s\",\n  \"m\": %d,\n  \"n\": %d,\n"
                  "  \"nnz\": %d,\n  \"status\": \"%s\",\n"
                  "  \"objective\": %.10g,\n  \"pinf\": %.3e,\n"
                  "  \"dinf\": %.3e,\n  \"rel_gap\": %.3e,\n"
                  "  \"iterations\": %ld,\n  \"solve_time_ms\": %.1f,\n"
                  "  \"message\": \"%s\"\n}\n",
                  json_escape(instance).c_str(), mdl.m, mdl.n, mdl.nnz(),
                  igaos::status_name(sol.status), sol.objective, sol.pinf,
                  sol.dinf, sol.rel_gap, sol.iterations, sol.solve_time_ms,
                  json_escape(sol.message).c_str());
    if (out_path.empty()) {
        std::fputs(buf, stdout);
    } else {
        std::ofstream of(out_path);
        of << buf;
    }
}

int cmd_info(const Model& mdl) {
    const double INF = std::numeric_limits<double>::infinity();
    int nE = 0, nL = 0, nG = 0, nRange = 0;
    for (int i = 0; i < mdl.m; ++i) {
        if (mdl.rmin[i] == mdl.rmax[i]) ++nE;
        else {
            if (std::isfinite(mdl.rmin[i])) ++nG;
            if (std::isfinite(mdl.rmax[i])) ++nL;
            if (std::isfinite(mdl.rmin[i]) && std::isfinite(mdl.rmax[i]) &&
                mdl.rmin[i] != -INF && mdl.rmax[i] != INF)
                ++nRange;
        }
    }
    int nfree = 0, nfixed = 0, nboxed = 0, nupper = 0;
    for (int j = 0; j < mdl.n; ++j) {
        bool fl = std::isfinite(mdl.cl[j]), fu = std::isfinite(mdl.cu[j]);
        if (!fl && !fu) ++nfree;
        else if (fl && fu && mdl.cl[j] == mdl.cu[j]) ++nfixed;
        else if (fl && fu) ++nboxed;
        else if (fl || fu) ++nupper;
    }
    double cnz_max = 0;
    for (double v : mdl.c) cnz_max = std::max(cnz_max, std::fabs(v));
    std::printf("rows           %d\n", mdl.m);
    std::printf("cols           %d\n", mdl.n);
    std::printf("nnz            %d\n", mdl.nnz());
    std::printf("row types      E=%d L=%d G=%d ranged=%d\n", nE,
                nL - nRange > 0 ? nL - nRange : 0, nG, nRange);
    std::printf("col bounds     free=%d fixed=%d boxed=%d one-sided=%d\n",
                nfree, nfixed, nboxed, nupper);
    std::printf("objective      |c|max=%.6g\n", cnz_max);
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        usage();
        return 2;
    }
    if (std::strcmp(argv[1], "--help") == 0 || std::strcmp(argv[1], "-h") == 0) {
        usage();
        return 0;
    }
    Args a = parse_args(argc, argv);
    if (!a.ok) {
        usage();
        return 2;
    }
    Model mdl;
    try {
        igaos::io::parse_mps(a.model_path, mdl);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 3;
    }
    if (a.cmd == "info") return cmd_info(mdl);

    igaos::Solution sol;
    sol.status = igaos::Status::Error;
    sol.message = "interface prototype: no engine linked in this build";
    emit_json(mdl, sol, a.model_path, a.out_path);
    return 0;
}
