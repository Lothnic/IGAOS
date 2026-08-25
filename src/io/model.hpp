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

namespace igaos::io {

struct Model {
    int m = 0, n = 0;
    std::vector<double> rmin, rmax;
    std::vector<double> cl, cu;
    std::vector<double> c;
    std::vector<int> ap, ai;
    std::vector<double> ax;
    std::vector<int> cp, ci;
    std::vector<double> acx;
    std::vector<unsigned char> integ;

    int nnz() const { return static_cast<int>(ax.size()); }
};

inline void counting_sort(const std::vector<std::tuple<int, int, double>>& trips,
                          int rows, int cols, std::vector<int>& ap,
                          std::vector<int>& ai, std::vector<double>& ax,
                          std::vector<int>& cp, std::vector<int>& ci,
                          std::vector<double>& acx) {
    ap.assign(rows + 1, 0);
    cp.assign(cols + 1, 0);
    for (const auto& t : trips) {
        ap[std::get<0>(t) + 1]++;
        cp[std::get<1>(t) + 1]++;
    }
    for (int i = 0; i < rows; ++i) ap[i + 1] += ap[i];
    for (int j = 0; j < cols; ++j) cp[j + 1] += cp[j];
    std::vector<int> nr = ap, nc = cp;
    ai.resize(trips.size());
    ax.resize(trips.size());
    ci.resize(trips.size());
    acx.resize(trips.size());
    for (const auto& t : trips) {
        int r = std::get<0>(t), c = std::get<1>(t);
        double v = std::get<2>(t);
        ai[nr[r]] = c;
        ax[nr[r]] = v;
        nr[r]++;
        ci[nc[c]] = r;
        acx[nc[c]] = v;
        nc[c]++;
    }
}

inline void parse_mps(const std::string& path, Model& model) {
    const double INF = std::numeric_limits<double>::infinity();
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::string line, section, objname;
    std::unordered_map<std::string, int> rowidx, colidx;
    std::vector<char> rsense;
    std::vector<double> rhs;
    std::vector<std::tuple<int, int, double>> trips;
    std::vector<std::pair<int, double>> ranges;
    bool pending_int = false;

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
            if (s == 'N') {
                if (objname.empty()) objname = t[1];
                rowidx[t[1]] = -1;
            } else {
                rowidx[t[1]] = static_cast<int>(rsense.size());
                rsense.push_back(s);
                rhs.push_back(0.0);
            }
        } else if (section == "COLUMNS") {
            if (line.find("'MARKER'") != std::string::npos) {
                pending_int = line.find("INTORG") != std::string::npos;
                continue;
            }
            auto ci = colidx.find(t[0]);
            int j;
            if (ci == colidx.end()) {
                j = static_cast<int>(colidx.size());
                colidx[t[0]] = j;
                model.c.push_back(0.0);
                model.cl.push_back(0.0);
                model.cu.push_back(INF);
                model.integ.push_back(
                    static_cast<unsigned char>(pending_int ? 1 : 0));
            } else {
                j = ci->second;
            }
            for (size_t k = 1; k + 1 < t.size(); k += 2) {
                double v = std::stod(t[k + 1]);
                if (t[k] == objname) {
                    model.c[j] += v;
                    continue;
                }
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
            if (type == "UP") {
                model.cu[j] = v;
                if (v < 0 && model.cl[j] == 0.0) model.cl[j] = -INF;
            } else if (type == "LO") {
                model.cl[j] = v;
            } else if (type == "FX") {
                model.cl[j] = model.cu[j] = v;
            } else if (type == "FR") {
                model.cl[j] = -INF;
                model.cu[j] = INF;
            } else if (type == "MI") {
                model.cl[j] = -INF;
            } else if (type == "PL") {
                model.cu[j] = INF;
            } else if (type == "BV") {
                model.cl[j] = 0.0;
                model.cu[j] = 1.0;
                model.integ[j] = 1;
            } else if (type == "LI") {
                model.cl[j] = v;
                model.integ[j] = 1;
            } else if (type == "UI") {
                model.cu[j] = v;
                model.integ[j] = 1;
            }
        }
    }

    model.m = static_cast<int>(rsense.size());
    model.n = static_cast<int>(colidx.size());
    model.rmin.assign(model.m, -INF);
    model.rmax.assign(model.m, INF);
    for (int i = 0; i < model.m; ++i) {
        if (rsense[i] == 'E') model.rmin[i] = model.rmax[i] = rhs[i];
        else if (rsense[i] == 'L') model.rmax[i] = rhs[i];
        else model.rmin[i] = rhs[i];
    }
    for (const auto& pr : ranges) {
        int i = pr.first;
        double r = pr.second, b = rhs[i];
        if (rsense[i] == 'E') {
            if (r >= 0) { model.rmin[i] = b; model.rmax[i] = b + r; }
            else { model.rmin[i] = b + r; model.rmax[i] = b; }
        } else if (rsense[i] == 'G') {
            model.rmin[i] = b;
            model.rmax[i] = b + std::fabs(r);
        } else {
            model.rmin[i] = b - std::fabs(r);
            model.rmax[i] = b;
        }
    }
    counting_sort(trips, model.m, model.n, model.ap, model.ai, model.ax,
                  model.cp, model.ci, model.acx);
}

}
