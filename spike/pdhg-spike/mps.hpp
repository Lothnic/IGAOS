#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

struct LP {
    int m = 0, n = 0;
    std::vector<double> rmin, rmax;
    std::vector<double> cl, cu;
    std::vector<double> c;
    std::vector<int> ap, ai, cp, ci;
    std::vector<double> ax, acx;
};

inline void counting_sort_csr(const std::vector<std::tuple<int,int,double>>& trips,
                              int rows, int cols,
                              std::vector<int>& ap, std::vector<int>& ai,
                              std::vector<double>& ax,
                              std::vector<int>& cp, std::vector<int>& ci,
                              std::vector<double>& acx) {
    ap.assign(rows + 1, 0);
    cp.assign(cols + 1, 0);
    for (const auto& t : trips) { ap[std::get<0>(t) + 1]++; cp[std::get<1>(t) + 1]++; }
    for (int i = 0; i < rows; ++i) ap[i + 1] += ap[i];
    for (int j = 0; j < cols; ++j) cp[j + 1] += cp[j];
    std::vector<int> nr = ap, nc = cp;
    ai.resize(trips.size()); ci.resize(trips.size());
    ax.resize(trips.size()); acx.resize(trips.size());
    for (const auto& t : trips) {
        int r = std::get<0>(t), c = std::get<1>(t);
        double v = std::get<2>(t);
        ai[nr[r]] = c; ax[nr[r]] = v; nr[r]++;
        ci[nc[c]] = r; acx[nc[c]] = v; nc[c]++;
    }
}

inline void parse_mps(const std::string& path, LP& lp) {
    const double INF = std::numeric_limits<double>::infinity();
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::string line, section, objname;
    std::unordered_map<std::string, int> rowidx, colidx;
    std::vector<char> rsense;
    std::vector<double> rhs;
    std::vector<std::tuple<int, int, double>> trips;
    std::vector<std::pair<int, double>> ranges;

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '*') continue;
        bool header = !(line[0] == ' ' || line[0] == '\t');
        std::istringstream ss(line);
        std::vector<std::string> t;
        for (std::string w; ss >> w;) t.push_back(w);
        if (t.empty()) continue;
        if (header) { section = t[0]; continue; }

        if (section == "ROWS") {
            const char s = t[0][0];
            if (s == 'N') { if (objname.empty()) objname = t[1]; rowidx[t[1]] = -1; }
            else { rowidx[t[1]] = (int)rsense.size(); rsense.push_back(s); rhs.push_back(0.0); }
        } else if (section == "COLUMNS") {
            if (line.find("'MARKER'") != std::string::npos) continue;
            auto ci = colidx.find(t[0]);
            int j;
            if (ci == colidx.end()) {
                j = (int)colidx.size();
                colidx[t[0]] = j;
                lp.c.push_back(0.0);
                lp.cl.push_back(0.0);
                lp.cu.push_back(INF);
            } else j = ci->second;
            for (size_t k = 1; k + 1 < t.size(); k += 2) {
                double v = std::stod(t[k + 1]);
                if (t[k] == objname) { lp.c[j] += v; continue; }
                auto rit = rowidx.find(t[k]);
                if (rit == rowidx.end() || rit->second < 0) continue;
                trips.emplace_back(rit->second, j, v);
            }
        } else if (section == "RHS") {
            for (size_t k = 1; k + 1 < t.size(); k += 2) {
                auto rit = rowidx.find(t[k]);
                if (rit == rowidx.end() || rit->second < 0) continue;
                rhs[rit->second] = std::stod(t[k + 1]);
            }
        } else if (section == "RANGES") {
            for (size_t k = 1; k + 1 < t.size(); k += 2) {
                auto rit = rowidx.find(t[k]);
                if (rit == rowidx.end() || rit->second < 0) continue;
                ranges.emplace_back(rit->second, std::stod(t[k + 1]));
            }
        } else if (section == "BOUNDS") {
            const std::string& type = t[0];
            int j = colidx.at(t[2]);
            double v = (t.size() > 3) ? std::stod(t[3]) : 0.0;
            if (type == "UP") { lp.cu[j] = v; if (v < 0 && lp.cl[j] == 0.0) lp.cl[j] = -INF; }
            else if (type == "LO") lp.cl[j] = v;
            else if (type == "FX") { lp.cl[j] = lp.cu[j] = v; }
            else if (type == "FR") { lp.cl[j] = -INF; lp.cu[j] = INF; }
            else if (type == "MI") lp.cl[j] = -INF;
            else if (type == "PL") lp.cu[j] = INF;
            else if (type == "BV") { lp.cl[j] = 0.0; lp.cu[j] = 1.0; }
            else if (type == "LI") lp.cl[j] = v;
            else if (type == "UI") lp.cu[j] = v;
        }
    }

    lp.m = (int)rsense.size();
    lp.n = (int)colidx.size();
    lp.rmin.assign(lp.m, -INF);
    lp.rmax.assign(lp.m, INF);
    for (int i = 0; i < lp.m; ++i) {
        if (rsense[i] == 'E') lp.rmin[i] = lp.rmax[i] = rhs[i];
        else if (rsense[i] == 'L') lp.rmax[i] = rhs[i];
        else lp.rmin[i] = rhs[i];
    }
    for (const auto& pr : ranges) {
        int i = pr.first;
        double r = pr.second, b = rhs[i];
        if (rsense[i] == 'E') {
            if (r >= 0) { lp.rmin[i] = b; lp.rmax[i] = b + r; }
            else { lp.rmin[i] = b + r; lp.rmax[i] = b; }
        } else if (rsense[i] == 'G') { lp.rmin[i] = b; lp.rmax[i] = b + std::fabs(r); }
        else { lp.rmin[i] = b - std::fabs(r); lp.rmax[i] = b; }
    }
    counting_sort_csr(trips, lp.m, lp.n, lp.ap, lp.ai, lp.ax, lp.cp, lp.ci, lp.acx);
}
