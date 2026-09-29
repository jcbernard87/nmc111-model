// NMC111 half-cell model: Li foil | separator | porous NMC111 cathode, with uniform particles
// or porous agglomerates of crystals (particle_model = 'uniform' | 'agglomerate').
//
// A self-contained C++17 program: edit the input file, run, and it writes Time_Voltage.txt.
// It reads the same namelist input as the Fortran program and the Python package, and
// contains its own copy of the BAND block-tridiagonal solver (adapted from bandsolver
// v0.1.2, https://github.com/jcbernard87/bandsolver, BSD-3-Clause), so it needs no library.
//
//   build:  c++ -std=c++17 -O2 -ffp-contract=off nmc.cpp -o nmc_cpp
//   usage:  nmc_cpp [input.nml]          (default input file: nmc.nml)
//
// Equations: docs/model.md. Parameters: docs/parameters.md. mode = 'faithful' reproduces
// the original research code including the defects in docs/deviations.md; mode =
// 'corrected' applies the fixes.
//
// SPDX-License-Identifier: BSD-3-Clause

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// ============================== BAND solver ==============================
// Solves A[j] dc[j-1] + B[j] dc[j] + D[j] dc[j+1] = G[j], j = 0..nj-1 (row-major n x n blocks).
namespace band {

enum class Pivot { partial, legacy };

// Solve Bm * S = R in place (Bm n x n, R n x m, row-major). Returns false if singular.
bool solve_partial(int n, int m, double* Bm, double* R) {
    const double eps = std::numeric_limits<double>::epsilon();
    double bscale = 0;
    for (int i = 0; i < n * n; ++i) bscale = std::max(bscale, std::abs(Bm[i]));
    if (bscale == 0) return false;
    for (int k = 0; k < n; ++k) {
        int p = k;
        for (int i = k + 1; i < n; ++i)
            if (std::abs(Bm[i * n + k]) > std::abs(Bm[p * n + k])) p = i;
        if (std::abs(Bm[p * n + k]) <= n * eps * bscale) return false;
        if (p != k) {
            std::swap_ranges(Bm + p * n, Bm + p * n + n, Bm + k * n);
            std::swap_ranges(R + p * m, R + p * m + m, R + k * m);
        }
        double f = 1.0 / Bm[k * n + k];
        for (int c = k; c < n; ++c) Bm[k * n + c] *= f;
        for (int c = 0; c < m; ++c) R[k * m + c] *= f;
        for (int i = 0; i < n; ++i) {
            if (i == k) continue;
            f = Bm[i * n + k];
            if (f == 0) continue;
            for (int c = k; c < n; ++c) Bm[i * n + c] -= f * Bm[k * n + c];
            for (int c = 0; c < m; ++c) R[i * m + c] -= f * R[k * m + c];
        }
    }
    return true;
}

// Pivot heuristic of the archival MATINV (Newman, Appendix C), with the operation order kept.
bool solve_legacy(int n, int m, double* Bm, double* R, bool archival = false) {
    const double eps = std::numeric_limits<double>::epsilon();
    const double bmax0 = static_cast<double>(1.1f);
    double bscale = 0;
    for (int i = 0; i < n * n; ++i) bscale = std::max(bscale, std::abs(Bm[i]));
    std::vector<char> used(n, 0);
    int irow = 0, jcol = 0, jc = 0;
    for (int nn = 0; nn < n; ++nn) {
        double bmax = bmax0;
        bool found = false;
        for (int i = 0; i < n; ++i) {
            if (used[i]) continue;
            double bnext = 0, btry = 0;
            for (int j = 0; j < n; ++j) {
                if (used[j]) continue;
                const double a = std::abs(Bm[i * n + j]);
                if (a <= bnext) continue;
                bnext = a;
                if (bnext <= btry) continue;
                bnext = btry;
                btry = a;
                jc = j;
                found = true;
            }
            if (bnext >= bmax * btry) continue;
            bmax = bnext / btry;
            irow = i;
            jcol = jc;
        }
        if (!found) return false;
        used[jcol] = 1;
        if (jcol != irow) {
            std::swap_ranges(Bm + irow * n, Bm + irow * n + n, Bm + jcol * n);
            std::swap_ranges(R + irow * m, R + irow * m + m, R + jcol * m);
        }
        // archival (bandsolver singular="exact"): like the original MATINV, a block is singular only
        // when no nonzero pivot is left
        if (!archival && std::abs(Bm[jcol * n + jcol]) <= n * eps * bscale) return false;
        double f = 1.0 / Bm[jcol * n + jcol];
        for (int j = 0; j < n; ++j) Bm[jcol * n + j] *= f;
        for (int k = 0; k < m; ++k) R[jcol * m + k] *= f;
        for (int i = 0; i < n; ++i) {
            if (i == jcol) continue;
            f = Bm[i * n + jcol];
            for (int j = 0; j < n; ++j) Bm[i * n + j] -= f * Bm[jcol * n + j];
            for (int k = 0; k < m; ++k) R[i * m + k] -= f * R[jcol * m + k];
        }
    }
    return true;
}

// Block elimination and back substitution. Returns false on a singular or non-finite system.
bool solve(int n, int nj, const std::vector<double>& A, const std::vector<double>& B,
           const std::vector<double>& D, const std::vector<double>& G, std::vector<double>& dc, Pivot pivot,
           bool archival = false) {
    const std::size_t nn = static_cast<std::size_t>(n) * n;
    dc.assign(static_cast<std::size_t>(n) * nj, 0.0);
    for (const auto* v : {&A, &B, &D, &G})
        for (double x : *v)
            if (!std::isfinite(x)) return false;
    const int np1 = n + 1;
    std::vector<double> E(static_cast<std::size_t>(nj) * n * np1), Am(nn), Bm(nn), R(static_cast<std::size_t>(n) * np1);
    auto Ej = [&](int j) { return E.data() + static_cast<std::size_t>(j) * n * np1; };
    auto blk = [&](const std::vector<double>& M, int j) { return M.data() + static_cast<std::size_t>(j) * nn; };
    auto block_solve = [&](int m) {
        return pivot == Pivot::legacy ? solve_legacy(n, m, Bm.data(), R.data(), archival) : solve_partial(n, m, Bm.data(), R.data());
    };
    // node 0: B0 [E0 | e0] = [D0 | G0]
    std::copy(blk(B, 0), blk(B, 0) + nn, Bm.begin());
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n; ++k) R[i * np1 + k] = blk(D, 0)[i * n + k];
        R[i * np1 + n] = G[i];
    }
    if (!block_solve(np1)) return false;
    for (int k = 0; k < n; ++k) {
        Ej(0)[k * np1 + n] = R[k * np1 + n];
        for (int l = 0; l < n; ++l) Ej(0)[k * np1 + l] = -R[k * np1 + l];
    }
    for (int j = 1; j < nj; ++j) {
        std::copy(blk(A, j), blk(A, j) + nn, Am.begin());
        std::copy(blk(B, j), blk(B, j) + nn, Bm.begin());
        for (int i = 0; i < n; ++i)
            for (int k = 0; k < n; ++k) R[i * np1 + k] = blk(D, j)[i * n + k];
        const double* E1 = Ej(j - 1);
        for (int i = 0; i < n; ++i) {
            R[i * np1 + n] = -G[static_cast<std::size_t>(j) * n + i];
            for (int l = 0; l < n; ++l) {
                R[i * np1 + n] += Am[i * n + l] * E1[l * np1 + n];
                for (int k = 0; k < n; ++k) Bm[i * n + k] += Am[i * n + l] * E1[l * np1 + k];
            }
        }
        if (!block_solve(np1)) return false;
        for (int k = 0; k < n; ++k)
            for (int q = 0; q < np1; ++q) Ej(j)[k * np1 + q] = -R[k * np1 + q];
    }
    for (int k = 0; k < n; ++k) dc[static_cast<std::size_t>(nj - 1) * n + k] = Ej(nj - 1)[k * np1 + n];
    for (int j = nj - 2; j >= 0; --j) {
        const double* Ec = Ej(j);
        double* x = dc.data() + static_cast<std::size_t>(j) * n;
        const double* xn = x + n;
        for (int k = 0; k < n; ++k) {
            x[k] = Ec[k * np1 + n];
            for (int l = 0; l < n; ++l) x[k] += Ec[k * np1 + l] * xn[l];
        }
    }
    for (double x : dc)
        if (!std::isfinite(x)) return false;
    return true;
}

}  // namespace band

// ============================== input ==============================
double r32(double x) { return static_cast<double>(static_cast<float>(x)); }

// x**m for an integer m as gfortran evaluates it (libgcc __powidf2: repeated squaring, then a
// reciprocal for negative m), so that faithful mode matches the original bit for bit.
double powi(double x, int m) {
    unsigned n = m < 0 ? static_cast<unsigned>(-m) : static_cast<unsigned>(m);
    double y = (n % 2) ? x : 1.0;
    while (n >>= 1) {
        x = x * x;
        if (n % 2) y = y * x;
    }
    return m < 0 ? 1.0 / y : y;
}

// Redlich-Kister coefficients of the uniform-particle model's OCP fit
const double AK_UNIFORM[11] = {-0.2018059457910574, 0.1123408808723528, -0.0483699097647364, 0.0231624989428732,
                               -0.0377897311905149, -0.3307806975105846, 0.2392976745148739, 0.7787126945566982,
                               -0.2599451275866008, -0.5898456896544948, 0.0520147453263591};
// Redlich-Kister coefficients of the agglomerate model's OCP fit
const double AK_AGG[12] = {-0.255139064974728, 0.0691287746986728, -0.1178158454270744, -0.0444434841626702,
                           0.243569591966704, 0.0775338354167729, -1.0934643144519782, -0.8893166395840808,
                           1.7690915896916977, 1.8213923583001588, -1.2074949744867922, -1.3952076583801158};

// ---- corrected mode: log variables and Scharfetter-Gummel fluxes (docs/model.md section 8) ----
double bern(double x) {  // B(x) = x/(e^x - 1), with its series near 0
    if (std::abs(x) < 1.0e-3) return 1.0 - x / 2.0 + x * x / 12.0;
    if (x > 700.0) return x * std::exp(-x);
    return x / (std::exp(x) - 1.0);
}
double bern_p(double x) {  // B'(x) = B(x) (1 - B(-x)) / x, with its series near 0
    if (std::abs(x) < 1.0e-3) return -0.5 + x / 6.0 - x * x * x / 180.0;
    return bern(x) * (1.0 - bern(-x)) / x;
}
double softplus(double x) { return std::max(x, 0.0) + std::log(1.0 + std::exp(-std::abs(x))); }
double sigm(double x) { return 0.5 * (1.0 + std::tanh(0.5 * x)); }

struct Params {
    std::string particle_model = "uniform";
    // agglomerate model (particle_model = 'agglomerate')
    int nja = 33;
    double time_mod = 20.0, R_agg = 1.0e-4 * 5.0, R_xtal = 200.0e-7, eps_agg = 0.2, D_agg = -1.0, tortuosity_e = -1.0;
    double mass_loading = 0.020, percent_active = 0.95, mol_vol = 0.0476881609, c0_init = 1.0e-3;
    double tau_agg = -1.0, sigma_agg = -1.0, dt_s = 1.0;  // corrected: < 0 means Bruggeman / sigma
    double kappa_bg = 1.0e-8;  // corrected: background (solvent) ionic conductivity [S/cm]
    double L_cath_um = 24.0, L_sep = 25.0e-4;  // cathode thickness in um
    int nj = 43, sep_node = 22;
    double eps = 0.5, eps_AM = 0.4, eps_sep = 0.39, tau_sep = 4.8, bruggeman = -0.5;
    double D = 2.0e-6, t_plus = 0.25, c_bulk = 1.0e-3, z_plus = 1.0, z_minus = -1.0;
    double sigma = -1.0, M = 96.46, rho = 4.6, Q_th = 0.150, R_p = 200.0e-7;
    double k_rxn = -1.0;  // < 0: default for the chosen mode
    double alpha_a = 0.5, alpha_c = 0.5, k_Li = 1.0e-6, c_Li_ref = 1.0e-3;
    double R = 8.314, T = 298.0, F = 96485.0;
    double C_rate = 1.0, phi1_init = 4.1, phi2_init = 0.0, cs_init = 1.0e-5, t_max = 36000.0;
    int n_steps = 36000;
    double V_min = 2.5, V_max = 4.2;       // corrected-mode cutoffs (D-3)
    double fd_step = 1.0e-6;
    double newton_tol = 1.0e-10;           // corrected-mode Newton tolerance (D-7)
    int newton_max_iter = 25;
    std::string mode = "faithful";
    std::string steps;                     // corrected-mode protocol (docs/protocol.md)
    int cycles = 1;
    double write_interval = 18.0;          // [s], corrected mode
    std::string file = "Time_Voltage.txt";
};

std::string lower(std::string s) {
    for (auto& ch : s) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    return s;
}

// Minimal Fortran-namelist reader: &group name = value, ... / ; '!' comments; d exponents.
Params read_input(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("input file not found: " + path);
    std::string text, line;
    while (std::getline(in, line)) {
        std::string out;
        char quote = 0;
        for (char ch : line) {
            if (quote) { if (ch == quote) quote = 0; }
            else if (ch == '\'' || ch == '"') quote = ch;
            else if (ch == '!') break;
            out += ch;
        }
        text += out + "\n";
    }
    auto parse = [&](Params& p) {
    std::map<std::string, double*> reals = {
        {"l_cath_um", &p.L_cath_um}, {"l_sep", &p.L_sep}, {"eps", &p.eps}, {"eps_am", &p.eps_AM},
        {"eps_sep", &p.eps_sep}, {"tau_sep", &p.tau_sep}, {"bruggeman", &p.bruggeman}, {"d", &p.D},
        {"t_plus", &p.t_plus}, {"c_bulk", &p.c_bulk}, {"z_plus", &p.z_plus}, {"z_minus", &p.z_minus},
        {"sigma", &p.sigma}, {"m", &p.M}, {"rho", &p.rho}, {"q_th", &p.Q_th}, {"r_p", &p.R_p},
        {"k_rxn", &p.k_rxn}, {"alpha_a", &p.alpha_a}, {"alpha_c", &p.alpha_c}, {"k_li", &p.k_Li},
        {"c_li_ref", &p.c_Li_ref}, {"r", &p.R}, {"t", &p.T}, {"f", &p.F}, {"c_rate", &p.C_rate},
        {"phi1_init", &p.phi1_init}, {"phi2_init", &p.phi2_init}, {"cs_init", &p.cs_init},
        {"t_max", &p.t_max}, {"fd_step", &p.fd_step}, {"v_min", &p.V_min}, {"v_max", &p.V_max},
        {"newton_tol", &p.newton_tol}, {"write_interval", &p.write_interval}, {"time_mod", &p.time_mod},
        {"r_agg", &p.R_agg}, {"r_xtal", &p.R_xtal}, {"eps_agg", &p.eps_agg}, {"d_agg", &p.D_agg},
        {"tortuosity_e", &p.tortuosity_e}, {"mass_loading", &p.mass_loading}, {"percent_active", &p.percent_active},
        {"mol_vol", &p.mol_vol}, {"c0_init", &p.c0_init}, {"tau_agg", &p.tau_agg}, {"sigma_agg", &p.sigma_agg},
        {"dt_s", &p.dt_s}, {"kappa_bg", &p.kappa_bg}};
    std::map<std::string, int*> ints = {{"nj", &p.nj}, {"sep_node", &p.sep_node}, {"n_steps", &p.n_steps},
                                        {"newton_max_iter", &p.newton_max_iter}, {"cycles", &p.cycles}, {"nja", &p.nja}};
    std::map<std::string, std::string*> strs = {{"mode", &p.mode}, {"file", &p.file}, {"steps", &p.steps},
                                                {"particle_model", &p.particle_model}};
    const std::regex group(R"(&(\w+)([\s\S]*?)/)");
    const std::regex entry(R"((\w+)\s*=\s*('[^']*'|"[^"]*"|[^,\s/]+))");
    for (std::sregex_iterator g(text.begin(), text.end(), group), end; g != end; ++g) {
        const std::string body = (*g)[2];
        for (std::sregex_iterator e(body.begin(), body.end(), entry); e != end; ++e) {
            const std::string name = lower((*e)[1]);
            std::string val = (*e)[2];
            if (strs.count(name)) {
                *strs[name] = (val.front() == '\'' || val.front() == '"') ? val.substr(1, val.size() - 2) : val;
            } else {
                std::replace(val.begin(), val.end(), 'd', 'e');
                std::replace(val.begin(), val.end(), 'D', 'e');
                if (ints.count(name)) *ints[name] = std::stoi(val);
                else if (reals.count(name)) *reals[name] = std::stod(val);
                else throw std::runtime_error("unknown name in input: " + name);
            }
        }
    }
    };
    Params p;
    parse(p);
    if (p.particle_model == "agglomerate") {  // defaults of the original conductivity-study template, then the file
        Params q;
        q.particle_model = "agglomerate";
        q.nj = 75; q.L_cath_um = 100.0; q.eps = 0.4; q.tau_sep = 4.0; q.D = 2.89e-6; q.t_plus = 0.375;
        q.sigma = 0.1; q.M = 96.46; q.rho = 4.7; q.Q_th = 0.155; q.phi1_init = 4.3; q.t_max = 72000.0;
        q.V_min = 3.0; q.V_max = 4.4;
        parse(q);
        p = q;
    }
    return p;
}

// ============================== model ==============================
constexpr int NV = 4, IC = 0, IP1 = 1, IP2 = 2, ICS = 3;
constexpr double THETA_REG = 1.0e-6;  // D-13 regularization threshold

// ============================== protocol ==============================
enum class Kind { cc, cv, rest };
struct Step {
    Kind kind;
    double C = 0.0, V = 0.0, t = -1.0, Vmin = 0.0, Vmax = 0.0, Imin = -1.0;  // t, Imin < 0: unset
};

std::vector<Step> parse_protocol(const std::string& text, double C_rate, double V_min, double V_max, int cycles) {
    std::vector<Step> steps;
    auto blank = [](const std::string& x) { return x.find_first_not_of(" \t") == std::string::npos; };
    if (blank(text)) {
        Step st{Kind::cc};
        st.C = C_rate; st.Vmin = V_min; st.Vmax = V_max;
        return {st};
    }
    std::stringstream parts(text);
    std::string part;
    int n = 0;
    while (std::getline(parts, part, ';')) {
        if (blank(part)) continue;
        ++n;
        auto fail = [&](const std::string& msg) { throw std::runtime_error("protocol step " + std::to_string(n) + ": " + msg); };
        std::istringstream words(part);
        std::string w;
        words >> w;
        Step st{Kind::cc};
        const std::string kind = lower(w);
        if (kind == "cc") st.kind = Kind::cc;
        else if (kind == "cv") st.kind = Kind::cv;
        else if (kind == "rest") st.kind = Kind::rest;
        else fail("unknown step type " + w);
        st.Vmin = V_min; st.Vmax = V_max;
        bool hasC = false, hasV = false;
        while (words >> w) {
            const auto eq = w.find('=');
            if (eq == std::string::npos) fail("expected key=value, got " + w);
            const std::string key = lower(w.substr(0, eq));
            std::string val = w.substr(eq + 1);
            std::replace(val.begin(), val.end(), 'd', 'e');
            std::replace(val.begin(), val.end(), 'D', 'e');
            double v;
            try { v = std::stod(val); } catch (...) { fail("bad number " + val); }
            if (st.kind == Kind::cc) {
                if (key == "c") { st.C = v; hasC = true; }
                else if (key == "t") st.t = v;
                else if (key == "vmin") st.Vmin = v;
                else if (key == "vmax") st.Vmax = v;
                else fail("cc does not take " + key);
            } else if (st.kind == Kind::cv) {
                if (key == "v") { st.V = v; hasV = true; }
                else if (key == "t") st.t = v;
                else if (key == "imin") st.Imin = v;
                else fail("cv does not take " + key);
            } else {
                if (key != "t") fail("rest does not take " + key);
                st.t = v;
            }
            if (key == "t" && v <= 0) fail("t must be positive");
        }
        if (st.kind == Kind::cc && !hasC) fail("cc needs C=");
        if (st.kind == Kind::cv && !hasV) fail("cv needs V=");
        if (st.kind == Kind::cv && st.t < 0 && st.Imin < 0) fail("cv needs t= or Imin= to end");
        if (st.kind == Kind::rest && st.t < 0) fail("rest needs t=");
        steps.push_back(st);
    }
    std::vector<Step> all;
    for (int k = 0; k < std::max(1, cycles); ++k) all.insert(all.end(), steps.begin(), steps.end());
    return all;
}
using Mat = double[NV][NV];

struct Model {
    Params p;
    bool faithful;
    double lit36, x_max, spec_a, tortuosity, i_spec, eps_sep_face, phi1_sign, mass_area, i_1C;
    double vf_AM, L_cath;  // active volume fraction; cathode thickness [cm]
    double dplus = 0.0, dminus = 0.0;  // ion diffusivities (corrected mode)
    double i_app;  // applied current density [A/cm2]; changes step by step in corrected mode
    bool full_current;
    double dcat_s, dan_s, ucat_s, uan_s, dcat_c, dan_c, ucat_c, uan_c;
    int s;  // 0-based interface node
    std::vector<double> dx, aW, aE, bW, bE;

    explicit Model(Params in) : p(std::move(in)) {
        if (p.mode == "faithful") faithful = true;
        else if (p.mode == "corrected") faithful = false;
        else throw std::runtime_error("mode must be 'faithful' or 'corrected'");
        if (faithful) {  // single-precision literals of the original (D-4)
            p.R = r32(p.R); p.c_bulk = r32(p.c_bulk); p.Q_th = r32(p.Q_th); p.M = r32(p.M); p.rho = r32(p.rho);
            p.phi1_init = r32(p.phi1_init); p.eps_sep = r32(p.eps_sep); p.eps_AM = r32(p.eps_AM);
            p.tau_sep = r32(p.tau_sep);
            p.C_rate = r32(p.C_rate);
            // the original's fitted values are not distributed: a faithful run must supply them
            if (p.k_rxn < 0 || p.sigma < 0) throw std::runtime_error("faithful uniform-particle runs need k_rxn and sigma in the input");
            lit36 = r32(3.6);
            x_max = r32(0.55);
            eps_sep_face = p.eps;   // D-2
            phi1_sign = -1.0;       // D-1
            full_current = false;   // D-11
        } else {
            if (p.k_rxn < 0) p.k_rxn = 2.5e-6;  // generic defaults (docs/parameters.md)
            if (p.sigma < 0) p.sigma = 0.1;
            lit36 = 3.6;
            x_max = 0.55;
            eps_sep_face = p.eps_sep;
            phi1_sign = 1.0;
            full_current = true;
        }
        L_cath = p.L_cath_um * 1.0e-4;  // the original writes 24 * 1.0d-4
        vf_AM = p.eps_AM;
        spec_a = 3 * vf_AM / p.R_p;
        tortuosity = std::pow(p.eps, p.bruggeman);
        i_spec = p.Q_th * p.C_rate;
        i_app = i_spec * L_cath * vf_AM * p.rho;
        mass_area = L_cath * vf_AM * p.rho;
        i_1C = p.Q_th * mass_area;
        const double t_an = 1.0 - p.t_plus;
        const double d_cat = p.D * (1.0 + (t_an / p.t_plus)) / (2.0 * t_an / p.t_plus);
        const double d_an = d_cat * t_an / p.t_plus;
        dplus = d_cat; dminus = d_an;
        const double u_cat = d_cat / (p.R * p.T), u_an = d_an / (p.R * p.T);
        dcat_s = d_cat / p.tau_sep; dan_s = d_an / p.tau_sep; ucat_s = u_cat / p.tau_sep; uan_s = u_an / p.tau_sep;
        dcat_c = d_cat / tortuosity; dan_c = d_an / tortuosity; ucat_c = u_cat / tortuosity; uan_c = u_an / tortuosity;

        const int nj = p.nj;
        s = p.sep_node - 1;
        dx.assign(nj, 0.0); aW.assign(nj, 0.0); aE.assign(nj, 0.0); bW.assign(nj, 0.0); bE.assign(nj, 0.0);
        const double h_sep = p.L_sep / static_cast<double>(p.sep_node - 2);
        const double h_cat = L_cath / static_cast<double>(nj - p.sep_node - 1);
        for (int j = 1; j < s; ++j) dx[j] = h_sep;
        for (int j = s + 1; j < nj - 1; ++j) dx[j] = h_cat;
        for (int j = 1; j < nj; ++j) { aW[j] = dx[j - 1] / (dx[j - 1] + dx[j]); bW[j] = 2.0 / (dx[j - 1] + dx[j]); }
        for (int j = 0; j < nj - 1; ++j) { aE[j] = dx[j] / (dx[j + 1] + dx[j]); bE[j] = 2.0 / (dx[j] + dx[j + 1]); }
    }

    // ---- kinetics ----
    double cs_max() const { return x_max * (p.rho / p.M); }
    // Redlich-Kister sum; accumulated in single precision in faithful mode (the original's Vint, D-4)
    double rk_sum(double th) const {
        float vint = 0.0f;
        double vd = 0.0;
        for (int k = 0; k <= 10; ++k) {
            const double ak = faithful ? r32(AK_UNIFORM[k]) : AK_UNIFORM[k];
            const double term = ak * (powi(2 * th - 1, k + 1) - (2 * th * k * (1 - th)) / powi(2 * th - 1, 1 - k));
            if (faithful) vint = static_cast<float>(static_cast<double>(vint) + term);
            else vd = vd + term;
        }
        return faithful ? static_cast<double>(vint) : vd;
    }
    // U = U_ref + (RT/F) ln[(c/c_bulk)(1-theta)/theta] + Redlich-Kister sum (docs/model.md section 2)
    double ocp(double c, double cs) const {
        const double th = (cs / (p.rho / p.M)) / x_max;
        const double u_ref = faithful ? r32(3.8685682447595453) : 3.8685682447595453;
        return u_ref + p.R * p.T / p.F * std::log(c / p.c_bulk * (1.0 - th) / th) + rk_sum(th);
    }
    // dU/dc and dU/dcs (corrected mode)
    void ocp_slopes(double c, double cs, double& du_dc, double& du_dcs) const {
        const double th = (cs / (p.rho / p.M)) / x_max, x = 2 * th - 1, rtf = p.R * p.T / p.F;
        double drk = 0.0;
        for (int k = 0; k <= 10; ++k) {
            double term = 2.0 * (2 * k + 1) * powi(x, k);
            if (k >= 2) term = term - 4.0 * k * (k - 1) * th * (1 - th) * powi(x, k - 2);
            drk = drk + AK_UNIFORM[k] * term;
        }
        du_dc = rtf / c;
        du_dcs = (rtf * (-1.0 / (1.0 - th) - 1.0 / th) + drk) / ((p.rho / p.M) * x_max);
    }
    double rate(double c, double cs, double p1, double p2) const {
        const double eta = p1 - p2 - ocp(c, cs);
        double i0 = p.F * p.k_rxn * std::pow(c, p.alpha_a) * std::pow(cs_max() - cs, p.alpha_a) * std::pow(cs, p.alpha_c);
        if (faithful) i0 = r32(i0);  // D-4
        return i0 * (std::exp(p.alpha_a * p.F * eta / (p.R * p.T)) - std::exp(-(p.alpha_c * p.F * eta / (p.R * p.T))));
    }
    // rate and finite-difference derivatives w.r.t. (c, phi1, phi2, cs) (D-6)
    // rate and exact derivatives; corrected mode (fixes D-6)
    // x^alpha, replaced below delta by a C1 quadratic with g(0) = 0 and a finite slope (D-13)
    static void power_reg(double x, double alpha, double delta, double& g, double& dg) {
        if (x < delta) {
            const double u = x / delta;
            g = std::pow(delta, alpha) * ((2.0 - alpha) * u + (alpha - 1.0) * u * u);
            dg = std::pow(delta, alpha - 1.0) * ((2.0 - alpha) + 2.0 * (alpha - 1.0) * u);
        } else {
            g = std::pow(x, alpha);
            dg = alpha * std::pow(x, alpha - 1.0);
        }
    }
    void rate_derivs_exact(double c, double cs, double p1, double p2, double& i, double di[NV]) const {
        const double rt = p.R * p.T, aa = p.alpha_a * p.F / rt, bb = p.alpha_c * p.F / rt;
        const double eta = p1 - p2 - ocp(c, cs);
        double du_dc, du_dcs;
        ocp_slopes(c, cs, du_dc, du_dcs);
        double gv, dgv, gs, dgs;
        power_reg(cs_max() - cs, p.alpha_a, THETA_REG * cs_max(), gv, dgv);
        power_reg(cs, p.alpha_c, THETA_REG * cs_max(), gs, dgs);
        const double pre = p.F * p.k_rxn * std::pow(c, p.alpha_a);
        const double i0 = pre * gv * gs, di0 = pre * (gs * -dgv + gv * dgs);
        const double ea = std::exp(aa * eta), ec = std::exp(-bb * eta);
        i = i0 * (ea - ec);
        const double di_deta = i0 * (aa * ea + bb * ec);
        di[IC] = p.alpha_a * i / c - di_deta * du_dc;
        di[ICS] = di0 * (ea - ec) - di_deta * du_dcs;
        di[IP1] = di_deta;
        di[IP2] = -di_deta;
    }
    void rate_derivs(double c, double cs, double p1, double p2, double& i, double di[NV]) const {
        if (!faithful) { rate_derivs_exact(c, cs, p1, p2, i, di); return; }
        const double h = p.fd_step;
        i = rate(c, cs, p1, p2);
        di[IC] = c <= h ? (rate(c + h, cs, p1, p2) - i) / h
                        : (rate(c + h, cs, p1, p2) - rate(c - h, cs, p1, p2)) / (2.0 * h);
        di[ICS] = cs <= h ? (rate(c, cs + h, p1, p2) - i) / h
                          : (rate(c, cs + h, p1, p2) - rate(c, cs - h, p1, p2)) / (2.0 * h);
        di[IP1] = (rate(c, cs, p1 + h, p2) - rate(c, cs, p1 - h, p2)) / (2.0 * h);
        di[IP2] = (rate(c, cs, p1, p2 + h) - rate(c, cs, p1, p2 - h)) / (2.0 * h);
    }

    // ---- corrected-mode time step ----
    // storage coefficients: the time-derivative part of each row is T * (c - c_old)
    std::vector<double> time_terms(double dt) const {
        const int nj = p.nj;
        std::vector<double> T(static_cast<std::size_t>(nj) * NV, 0.0);
        T[ICS] = -(vf_AM / dt);
        for (int j = 1; j < s; ++j) { T[j * NV + IC] = -(p.eps_sep / dt * dx[j]); T[j * NV + ICS] = -((1.0 - p.eps_sep) / dt); }
        for (int j = s + 1; j < nj - 1; ++j) T[j * NV + IC] = -((p.eps / dt) * dx[j]);
        for (int j = s; j < nj; ++j) T[j * NV + ICS] = -(vf_AM / dt);
        return T;
    }
    // largest step <= 1 keeping 0 < c and 0 < cs < cs_max (at most 90 % of the way to a bound)
    // and changing no potential by more than 0.1 V (Butler-Volmer exponentials make Newton overshoot)
    double bounded_step(const std::vector<double>& c, const std::vector<double>& d) const {
        const double keep = 0.9, csm = cs_max(), max_dphi = 0.1;  // potential change cap per iteration [V]
        double lam = 1.0, dphi = 0.0;
        for (int j = 0; j < p.nj; ++j)
            dphi = std::max({dphi, std::abs(d[j * NV + IP1]), std::abs(d[j * NV + IP2])});
        if (dphi > max_dphi) lam = max_dphi / dphi;
        for (int j = 0; j < p.nj; ++j) {
            const double dcc = d[j * NV + IC], dcs = d[j * NV + ICS];
            if (dcc < 0) lam = std::min(lam, keep * (c[j * NV + IC] - 0.0) / -dcc);
            if (dcs < 0) lam = std::min(lam, keep * (c[j * NV + ICS] - 0.0) / -dcs);
            if (dcs > 0) lam = std::min(lam, keep * (csm - c[j * NV + ICS]) / dcs);
        }
        return std::max(lam, 0.0);
    }

    // ---- assembly ----
    void face(double e, double dcat, double ucat, double dan, double uan, double cf, double gphi2, Mat& dd, Mat& ff) const {
        const double k = p.z_plus * p.z_plus * ucat + p.z_minus * p.z_minus * uan;
        dd[IC][IC] = -(e * dcat);
        ff[IC][IC] = -(e * p.z_plus * ucat * p.F * gphi2);
        dd[IC][IP2] = -(e * p.z_plus * ucat * p.F * cf);
        dd[IP2][IC] = -(e * p.F * (p.z_plus * dcat + p.z_minus * dan));
        ff[IP2][IC] = -(e * (p.F * p.F) * k * gphi2);
        dd[IP2][IP2] = -(e * (p.F * p.F) * k * cf);
    }
    double current(const Mat& ff, const Mat& dd, const double* cf, const double* gf) const {
        double v = ff[IP2][IP2] * cf[IP2] + dd[IP2][IP2] * gf[IP2];
        if (full_current) v = v + dd[IP2][IC] * gf[IC];  // D-11
        return v;
    }
    void solid_row(double i, const double di[NV], double dt, Mat& rj, double* g) const {
        for (int k = 0; k < NV; ++k) rj[ICS][k] = -(spec_a * di[k] / p.F);
        rj[ICS][ICS] = -(spec_a * di[ICS] / p.F) - 1.0 * vf_AM / dt;
        g[ICS] = +(spec_a * i / p.F);
    }

    void assemble(const std::vector<double>& c, double dt, std::vector<double>& A, std::vector<double>& B,
                  std::vector<double>& Dm, std::vector<double>& G) const {
        const int nj = p.nj;
        const double F = p.F, a = spec_a, eps = p.eps, eps_sep = p.eps_sep, sig = p.sigma;
        std::fill(A.begin(), A.end(), 0.0); std::fill(B.begin(), B.end(), 0.0);
        std::fill(Dm.begin(), Dm.end(), 0.0); std::fill(G.begin(), G.end(), 0.0);
        auto C = [&](int j, int k) { return c[static_cast<std::size_t>(j) * NV + k]; };
        for (int j = 0; j < nj; ++j) {
            Mat dW{}, dE{}, fW{}, fE{}, rj{};
            double g[NV] = {0, 0, 0, 0}, cW[NV] = {0, 0, 0, 0}, cE[NV] = {0, 0, 0, 0}, gW[NV] = {0, 0, 0, 0},
                   gE[NV] = {0, 0, 0, 0}, i, di[NV];
            for (int k = 0; k < NV; ++k) {
                if (j > 0) { cW[k] = aW[j] * C(j, k) + (1.0 - aW[j]) * C(j - 1, k); gW[k] = bW[j] * (C(j, k) - C(j - 1, k)); }
                if (j < nj - 1) { cE[k] = aE[j] * C(j + 1, k) + (1.0 - aE[j]) * C(j, k); gE[k] = bE[j] * (C(j + 1, k) - C(j, k)); }
            }
            rate_derivs(C(j, IC), C(j, ICS), C(j, IP1), C(j, IP2), i, di);
            const double dxj = dx[j];

            if (j == 0) {  // Li-foil face
                dE[IC][IC] = -(eps_sep_face * dcat_s);
                fE[IC][IC] = -(eps_sep_face * p.z_plus * ucat_s * F * gE[IP2]);
                dE[IC][IP2] = -(eps_sep_face * p.z_plus * ucat_s * F * cE[IC]);
                g[IC] = -i_app / F + (dE[IC][IC] * gE[IC] + fE[IC][IC] * cE[IC]);
                rj[ICS][ICS] = 0.0 - 1.0 * vf_AM / dt;
                dE[IP1][IP1] = -(1.0 - eps_sep_face) * sig;
                g[IP1] = phi1_sign * dE[IP1][IP1] * gE[IP1];
                rj[IP2][IP2] = 1.0;
                g[IP2] = 0.0 - C(j, IP2);
            } else if (j < s) {  // separator interior
                face(eps_sep, dcat_s, ucat_s, dan_s, uan_s, cW[IC], gW[IP2], dW, fW);
                face(eps_sep, dcat_s, ucat_s, dan_s, uan_s, cE[IC], gE[IP2], dE, fE);
                rj[IC][IC] = -(eps_sep / dt * dxj);
                g[IC] = 0.0 - (fW[IC][IC] * cW[IC] + dW[IC][IC] * gW[IC]) + (fE[IC][IC] * cE[IC] + dE[IC][IC] * gE[IC]);
                rj[ICS][ICS] = -(1.0 * (1.0 - eps_sep) / dt);
                dW[IP1][IP1] = -(1.0 - eps_sep) * sig;
                dE[IP1][IP1] = -(1.0 - eps_sep) * sig;
                g[IP1] = 0.0 - (fW[IP1][IP1] * cW[IP1] + dW[IP1][IP1] * gW[IP1]) + (fE[IP1][IP1] * cE[IP1] + dE[IP1][IP1] * gE[IP1]);
                g[IP2] = 0.0 - current(fW, dW, cW, gW) + current(fE, dE, cE, gE);
            } else if (j == s) {  // separator/cathode interface
                face(eps_sep_face, dcat_s, ucat_s, dan_s, uan_s, cW[IC], gW[IP2], dW, fW);
                face(eps, dcat_c, ucat_c, dan_c, uan_c, cE[IC], gE[IP2], dE, fE);
                g[IC] = 0.0 - (fW[IC][IC] * cW[IC] + dW[IC][IC] * gW[IC]) + (fE[IC][IC] * cE[IC] + dE[IC][IC] * gE[IC]);
                dW[IP1][IP1] = -(1.0 - eps_sep_face) * sig;
                dE[IP1][IP1] = -(1.0 - eps) * sig;
                g[IP1] = 0.0 - (fW[IP1][IP1] * cW[IP1] + dW[IP1][IP1] * gW[IP1]) + (fE[IP1][IP1] * cE[IP1] + dE[IP1][IP1] * gE[IP1]);
                g[IP2] = 0.0 - current(fW, dW, cW, gW) + current(fE, dE, cE, gE);
                solid_row(i, di, dt, rj, g);
            } else if (j < nj - 1) {  // cathode interior
                face(eps, dcat_c, ucat_c, dan_c, uan_c, cW[IC], gW[IP2], dW, fW);
                face(eps, dcat_c, ucat_c, dan_c, uan_c, cE[IC], gE[IP2], dE, fE);
                for (int k = 0; k < NV; ++k) {
                    rj[IC][k] = (a * di[k] / F) * dxj;
                    rj[IP1][k] = -((a * di[k]) * dxj);
                    rj[IP2][k] = (a * di[k]) * dxj;
                }
                rj[IC][IC] = (a * di[IC] / F) * dxj - (eps / dt) * dxj;
                g[IC] = -((a * i / F) * dxj) - (fW[IC][IC] * cW[IC] + dW[IC][IC] * gW[IC]) + (fE[IC][IC] * cE[IC] + dE[IC][IC] * gE[IC]);
                dW[IP1][IP1] = -(1.0 - eps) * sig;
                dE[IP1][IP1] = -(1.0 - eps) * sig;
                g[IP1] = (a * i) * dxj - (fW[IP1][IP1] * cW[IP1] + dW[IP1][IP1] * gW[IP1]) + (fE[IP1][IP1] * cE[IP1] + dE[IP1][IP1] * gE[IP1]);
                g[IP2] = -((a * i) * dxj) - current(fW, dW, cW, gW) + current(fE, dE, cE, gE);
                solid_row(i, di, dt, rj, g);
            } else {  // current collector
                face(eps, dcat_c, ucat_c, dan_c, uan_c, cW[IC], gW[IP2], dW, fW);
                g[IC] = 0.0 - dW[IC][IC] * gW[IC] - fW[IC][IC] * cW[IC];
                dW[IP1][IP1] = -(1.0 - eps) * sig;
                g[IP1] = i_app - dW[IP1][IP1] * gW[IP1];
                g[IP2] = 0.0 - dW[IP2][IC] * gW[IC] - fW[IP2][IC] * cW[IC];
                solid_row(i, di, dt, rj, g);
            }

            // control-volume blocks: A dc[j-1] + B dc[j] + D dc[j+1] = G
            const std::size_t o = static_cast<std::size_t>(j) * NV * NV;
            for (int r = 0; r < NV; ++r) {
                for (int k = 0; k < NV; ++k) {
                    const std::size_t q = o + r * NV + k;
                    if (j == 0) {
                        B[q] = rj[r][k] - (1.0 - aE[j]) * fE[r][k] + bE[j] * dE[r][k];
                        Dm[q] = -(aE[j] * fE[r][k]) - bE[j] * dE[r][k];
                    } else if (j == nj - 1) {
                        A[q] = (1.0 - aW[j]) * fW[r][k] - bW[j] * dW[r][k];
                        B[q] = rj[r][k] + bW[j] * dW[r][k] + aW[j] * fW[r][k];
                    } else {
                        A[q] = (1.0 - aW[j]) * fW[r][k] - bW[j] * dW[r][k];
                        B[q] = rj[r][k] + bW[j] * dW[r][k] + aW[j] * fW[r][k] - (1.0 - aE[j]) * fE[r][k] + bE[j] * dE[r][k];
                        Dm[q] = -(aE[j] * fE[r][k]) - bE[j] * dE[r][k];
                    }
                }
                G[static_cast<std::size_t>(j) * NV + r] = g[r];
            }
        }
    }

    // ---- corrected mode (log variables; docs/model.md section 8) ----
    // Butler-Volmer: i_n and d i_n / d(u, phi1, phi2, s); agg selects the agglomerate OCP fit and c_s,max
    void log_rate(const double* x, bool agg, double csmax_agg, double& i, double di[NV]) const {
        const double rtf = p.R * p.T / p.F, ff = 1.0 / rtf, th = sigm(x[ICS]), xx = 2 * th - 1;
        const double uref = agg ? 3.8637058886774844 : 3.8685682447595453, csm = agg ? csmax_agg : cs_max();
        const int nk = agg ? 12 : 11;
        double rk = 0.0, drk = 0.0;
        for (int k = 0; k < nk; ++k) {
            const double a = agg ? AK_AGG[k] : AK_UNIFORM[k];
            rk = rk + a * (powi(xx, k + 1) - (2 * th * k * (1 - th)) / powi(xx, 1 - k));
            double term = 2.0 * (2 * k + 1) * powi(xx, k);
            if (k >= 2) term = term - 4.0 * k * (k - 1) * th * (1 - th) * powi(xx, k - 2);
            drk = drk + a * term;
        }
        const double dU_ds = -rtf + drk * th * (1 - th);
        const double ln_i0 = std::log(p.F * p.k_rxn * std::pow(p.c_bulk, p.alpha_a) * std::pow(csm, p.alpha_a + p.alpha_c))
                             + p.alpha_a * x[IC] - p.alpha_a * softplus(x[ICS]) - p.alpha_c * softplus(-x[ICS]);
        const double eta = x[IP1] - x[IP2] - (uref + rtf * (x[IC] - x[ICS]) + rk);
        const double ea = std::exp(ln_i0 + p.alpha_a * ff * eta), ec = std::exp(ln_i0 - p.alpha_c * ff * eta);
        i = ea - ec;
        const double di_deta = p.alpha_a * ff * ea + p.alpha_c * ff * ec;
        di[IC] = i * p.alpha_a - di_deta * rtf;
        di[IP1] = di_deta;
        di[IP2] = -di_deta;
        di[ICS] = i * (-p.alpha_a * th + p.alpha_c * (1 - th)) - di_deta * dU_ds;
    }
    // face fluxes (N+, i1, i2) from xa to xb (Scharfetter-Gummel) and derivatives dFa, dFb [3][NV];
    // g = eps/(tau h), gs = (1-eps) sigma/h
    void sg_face(const double* xa, const double* xb, double g, double gs, double* Fv, double* dFa, double* dFb) const {
        const double F = p.F, ff = F / (p.R * p.T);
        const double ca = p.c_bulk * std::exp(xa[IC]), cb = p.c_bulk * std::exp(xb[IC]), d = ff * (xb[IP2] - xa[IP2]);
        const double Bp = bern(d), Bm = bern(-d), dBp = bern_p(d), dBm = bern_p(-d);
        const double Np = g * dplus * (Bp * ca - Bm * cb), Nm = g * dminus * (Bm * ca - Bp * cb);
        const double dNp_dd = g * dplus * (dBp * ca + dBm * cb), dNm_dd = g * dminus * (-dBm * ca - dBp * cb);
        const double kb = g * p.kappa_bg;
        Fv[0] = Np;
        Fv[1] = -gs * (xb[IP1] - xa[IP1]);
        Fv[2] = F * (Np - Nm) - kb * (xb[IP2] - xa[IP2]);
        for (int z = 0; z < 3 * NV; ++z) { dFa[z] = 0.0; dFb[z] = 0.0; }
        dFa[0 * NV + IC] = g * dplus * Bp * ca;
        dFb[0 * NV + IC] = -g * dplus * Bm * cb;
        dFa[0 * NV + IP2] = -ff * dNp_dd;
        dFb[0 * NV + IP2] = ff * dNp_dd;
        dFa[1 * NV + IP1] = gs;
        dFb[1 * NV + IP1] = -gs;
        dFa[2 * NV + IC] = F * (dFa[0 * NV + IC] - g * dminus * Bm * ca);
        dFb[2 * NV + IC] = F * (dFb[0 * NV + IC] + g * dminus * Bp * cb);
        dFa[2 * NV + IP2] = -F * ff * (dNp_dd - dNm_dd) + kb;
        dFb[2 * NV + IP2] = F * ff * (dNp_dd - dNm_dd) - kb;
    }
    // electrode residual and blocks (G = -R) without agglomerate sources; particles: the uniform model's
    // particles at nodes s..nj-1 (else the S column is held fixed)
    void log_assemble(const std::vector<double>& x, const std::vector<double>& xold, double dt, double Ia, bool particles,
                      std::vector<double>& A, std::vector<double>& B, std::vector<double>& Dm, std::vector<double>& G) const {
        const int nj = p.nj;
        const double F = p.F;
        std::fill(A.begin(), A.end(), 0.0); std::fill(B.begin(), B.end(), 0.0); std::fill(Dm.begin(), Dm.end(), 0.0);
        std::vector<double> R(static_cast<std::size_t>(nj) * NV, 0.0), Fv(3 * static_cast<std::size_t>(nj - 1)),
            dFa(3 * NV * static_cast<std::size_t>(nj - 1)), dFb(dFa.size());
        for (int k = 0; k < nj - 1; ++k) {
            const double h = (dx[k] + dx[k + 1]) / 2.0;
            const bool sep = k < s;
            const double g = (sep ? p.eps_sep / p.tau_sep : p.eps / tortuosity) / h;
            const double gs = (sep ? 1.0 - p.eps_sep : 1.0 - p.eps) * p.sigma / h;
            sg_face(&x[static_cast<std::size_t>(k) * NV], &x[static_cast<std::size_t>(k + 1) * NV], g, gs, &Fv[3 * k],
                    &dFa[3 * NV * k], &dFb[3 * NV * k]);
        }
        auto X = [&](int j, int v) { return x[static_cast<std::size_t>(j) * NV + v]; };
        auto blk = [&](std::vector<double>& M, int j, int r, int v) -> double& { return M[(static_cast<std::size_t>(j) * NV + r) * NV + v]; };
        auto FA = [&](int k, int f, int v) { return dFa[3 * NV * k + f * NV + v]; };
        auto FB = [&](int k, int f, int v) { return dFb[3 * NV * k + f * NV + v]; };
        auto Rr = [&](int j, int r) -> double& { return R[static_cast<std::size_t>(j) * NV + r]; };
        // foil face: N+ = I/F, zero electronic current, phi2 = 0 (the gauge)
        Rr(0, IC) = Fv[0] - Ia / F;
        for (int v = 0; v < NV; ++v) { blk(B, 0, IC, v) = FA(0, 0, v); blk(Dm, 0, IC, v) = FB(0, 0, v); }
        Rr(0, IP1) = X(1, IP1) - X(0, IP1);
        blk(B, 0, IP1, IP1) = -1.0; blk(Dm, 0, IP1, IP1) = 1.0;
        Rr(0, IP2) = X(0, IP2);
        blk(B, 0, IP2, IP2) = 1.0;
        // separator, interface and cathode
        for (int j = 1; j < nj - 1; ++j) {
            for (int f = 0; f < 3; ++f) {
                Rr(j, f) = Fv[3 * j + f] - Fv[3 * (j - 1) + f];
                for (int v = 0; v < NV; ++v) {
                    blk(B, j, f, v) += FA(j, f, v) - FB(j - 1, f, v);
                    blk(Dm, j, f, v) += FB(j, f, v);
                    blk(A, j, f, v) += -FA(j - 1, f, v);
                }
            }
            const double e = j < s ? p.eps_sep : p.eps;
            const double cj = p.c_bulk * std::exp(X(j, IC)), cjo = p.c_bulk * std::exp(xold[static_cast<std::size_t>(j) * NV + IC]);
            Rr(j, IC) += e * dx[j] * (cj - cjo) / dt;
            blk(B, j, IC, IC) += e * dx[j] * cj / dt;
        }
        // collector: no salt flux, no ionic current, electronic current I
        const int n = nj - 1, kl = nj - 2;
        Rr(n, IC) = -Fv[3 * kl];
        Rr(n, IP1) = Fv[3 * kl + 1] - Ia;
        Rr(n, IP2) = Fv[3 * kl + 2];
        for (int v = 0; v < NV; ++v) {
            blk(A, n, IC, v) = -FA(kl, 0, v); blk(B, n, IC, v) = -FB(kl, 0, v);
            blk(A, n, IP1, v) = FA(kl, 1, v); blk(B, n, IP1, v) = FB(kl, 1, v);
            blk(A, n, IP2, v) = FA(kl, 2, v); blk(B, n, IP2, v) = FB(kl, 2, v);
        }
        // the S column: fixed, or the uniform model's particles
        for (int j = 0; j < nj; ++j) {
            Rr(j, ICS) = X(j, ICS) - xold[static_cast<std::size_t>(j) * NV + ICS];
            for (int v = 0; v < NV; ++v) blk(B, j, ICS, v) = 0.0;
            blk(B, j, ICS, ICS) = 1.0;
        }
        if (particles) {
            for (int j = s; j < nj; ++j) {
                double i, di[NV];
                log_rate(&x[static_cast<std::size_t>(j) * NV], false, 0.0, i, di);
                const double a = spec_a;
                Rr(j, IC) -= a * i * dx[j] / F;
                Rr(j, IP1) += a * i * dx[j];
                Rr(j, IP2) -= a * i * dx[j];
                for (int v = 0; v < NV; ++v) {
                    blk(B, j, IC, v) -= a * di[v] * dx[j] / F;
                    blk(B, j, IP1, v) += a * di[v] * dx[j];
                    blk(B, j, IP2, v) -= a * di[v] * dx[j];
                }
                const double th = sigm(X(j, ICS)), tho = sigm(xold[static_cast<std::size_t>(j) * NV + ICS]);
                Rr(j, ICS) = vf_AM * cs_max() * (th - tho) / dt + a * i / F;
                for (int v = 0; v < NV; ++v) blk(B, j, ICS, v) = a * di[v] / F;
                blk(B, j, ICS, ICS) += vf_AM * cs_max() * th * (1 - th) / dt;
            }
        }
        for (std::size_t q = 0; q < G.size(); ++q) G[q] = -R[q];
    }
};

// Newton convergence: the scaled update is below tol, or it has stagnated at the round-off floor
// (within 1e3*tol and down by less than half since the last iteration). Near a full or empty particle
// the Jacobian is ill-conditioned and the linear solve's round-off can leave the update just above tol.
bool converged(double upd, double prev, double tol) {
    return upd <= tol || (upd <= 1.0e3 * tol && upd >= 0.5 * prev);
}

// Newton step length <= 1 limiting |du| <= 1, |dphi| <= 0.1 V and |ds| <= 2 per iteration
double lbounded(const std::vector<double>& d) {
    const double caps[NV] = {1.0, 0.1, 0.1, 2.0};
    double lam = 1.0;
    for (int k = 0; k < NV; ++k) {
        double mx = 0.0;
        for (std::size_t q = k; q < d.size(); q += NV) mx = std::max(mx, std::abs(d[q]));
        if (mx > caps[k]) lam = std::min(lam, caps[k] / mx);
    }
    return lam;
}

// the Newton update in the physical variables: max of e^u |du|, |dphi| and theta(1-theta) |ds|
double phys_update(const std::vector<double>& x, const std::vector<double>& d, bool frozen_s) {
    double u = 0.0;
    for (std::size_t q = 0; q < x.size(); q += NV) {
        u = std::max({u, std::exp(x[q + IC]) * std::abs(d[q + IC]), std::abs(d[q + IP1]), std::abs(d[q + IP2])});
        if (!frozen_s) {
            const double th = sigm(x[q + ICS]);
            u = std::max(u, th * (1 - th) * std::abs(d[q + ICS]));
        }
    }
    return u;
}

// scale every equation by the largest entry of its row in B (the solution is unchanged)
void equilibrate_blocks(int nj, std::vector<double>& A, std::vector<double>& B, std::vector<double>& Dm, std::vector<double>& G) {
    for (int j = 0; j < nj; ++j)
        for (int r = 0; r < NV; ++r) {
            const std::size_t o = (static_cast<std::size_t>(j) * NV + r) * NV;
            double sc = 0.0;
            for (int k = 0; k < NV; ++k) sc = std::max(sc, std::abs(B[o + k]));
            if (sc == 0.0) sc = 1.0;
            for (int k = 0; k < NV; ++k) { A[o + k] /= sc; B[o + k] /= sc; Dm[o + k] /= sc; }
            G[static_cast<std::size_t>(j) * NV + r] /= sc;
        }
}

// ============================== output ==============================
std::string fixed12(double v) {
    char b[64];
    if (std::isnan(v)) std::snprintf(b, sizeof b, "%12s", "NaN");
    else std::snprintf(b, sizeof b, "%12.5f", v);
    return b;
}
std::string sci15(double v) {
    char b[64];
    if (std::isnan(v)) std::snprintf(b, sizeof b, "%15s", "NaN");
    else std::snprintf(b, sizeof b, "%15.5E", v);
    return b;
}
std::string fit(const std::string& s, int w) {  // Fortran A<w>: right-justify, keep leftmost w characters
    const std::string t = s.substr(0, w);
    return std::string(w - t.size(), ' ') + t;
}


// ============================== agglomerate model ==============================
// Porous electrode whose particles are porous spherical agglomerates of uniform crystals
// (docs/model.md section 4). Faithful mode reproduces NMC111_agg.f95: single-precision
// constants, one linearized solve per step, the two scales solved one after the other.


struct Agglomerate {
    Params p;
    int nj, na, s;
    double PIr, pi4, pi43, pi43b, vol_agg, area_agg, v_AM, a_e, a_agg, i_final, i_specific, L_cath;
    double tmx, dt_nom, thr_dep, thr_dep2, thr_x, lit36, dcat, dan, u0, ucat, uan, uagg, D_agg, k_rxn;
    long long nsteps;
    std::vector<double> dx, aW, aE, bW, bE, dxa, xa, aWa, aEa, bWa, bEa, AWs, AEs, dV;
    std::vector<double> ce, dce, ca, dca, ie;          // ce (nj x 4), ca (nj x na x 4)
    std::vector<double> Ae, Be, De, Ge, Aa, Ba, Da, Ga;

    static void faces(const std::vector<double>& d, std::vector<double>& aw, std::vector<double>& ae,
                      std::vector<double>& bw, std::vector<double>& be) {
        const int n = static_cast<int>(d.size());
        aw.assign(n, 0.0); ae.assign(n, 0.0); bw.assign(n, 0.0); be.assign(n, 0.0);
        for (int j = 1; j < n; ++j) { aw[j] = d[j - 1] / (d[j - 1] + d[j]); bw[j] = 2.0 / (d[j - 1] + d[j]); }
        for (int j = 0; j < n - 1; ++j) { ae[j] = d[j] / (d[j + 1] + d[j]); be[j] = 2.0 / (d[j] + d[j + 1]); }
    }

    explicit Agglomerate(Params in) : p(std::move(in)) {
        if (p.mode != "faithful")
            throw std::runtime_error("the agglomerate model is available in faithful mode only (corrected mode: work in progress)");
        // single-precision literals of the original (D-4)
        p.R = r32(p.R); p.c_bulk = r32(p.c_bulk); p.c0_init = r32(p.c0_init); p.eps_sep = r32(p.eps_sep);
        p.M = r32(p.M); p.rho = r32(p.rho); p.mol_vol = r32(p.mol_vol); p.Q_th = r32(p.Q_th); p.eps_agg = r32(p.eps_agg);
        p.mass_loading = r32(p.mass_loading); p.percent_active = r32(p.percent_active); p.phi1_init = r32(p.phi1_init);
        p.eps = r32(p.eps); p.sigma = r32(p.sigma); p.C_rate = r32(p.C_rate); p.time_mod = r32(p.time_mod);
        // the original's fitted values are not distributed: a faithful run must supply them
        if (p.D_agg < 0 || p.k_rxn < 0 || p.tortuosity_e < 0)
            throw std::runtime_error("faithful agglomerate runs need D_agg, k_rxn and tortuosity_e in the input");
        D_agg = p.D_agg;
        k_rxn = p.k_rxn;
        PIr = r32(3.141592654);
        lit36 = r32(3.6);
        L_cath = static_cast<double>(static_cast<float>(p.L_cath_um) / 10000.0f);   // THICKNESS/10000.0
        tmx = static_cast<double>(static_cast<float>(p.t_max));                       // tmax is an implicitly REAL parameter
        nsteps = static_cast<long long>(3.6e3 * p.C_rate * p.time_mod);
        dt_nom = static_cast<double>(static_cast<float>(tmx) / static_cast<float>(nsteps));
        thr_dep = r32(0.0001); thr_dep2 = r32(0.00001); thr_x = r32(0.55);
        v_AM = p.percent_active * p.mass_loading / (p.rho * L_cath);
        a_e = 3.0 * v_AM / p.R_agg;
        a_agg = 3.0 * (1.0 - p.eps_agg) / p.R_xtal;
        i_specific = p.Q_th * p.C_rate;
        i_final = i_specific * L_cath * v_AM * p.rho;
        const double t_an = 1.0 - p.t_plus, diff_e = p.eps * p.D / p.tortuosity_e;
        dcat = diff_e * (1.0 + (t_an / p.t_plus)) / (2.0 * t_an / p.t_plus);
        dan = dcat * t_an / p.t_plus;
        u0 = p.D / (p.R * p.T); ucat = dcat / (p.R * p.T); uan = dan / (p.R * p.T); uagg = D_agg / (p.R * p.T);
        nj = p.nj; na = p.nja; s = p.sep_node - 1;
        dx.assign(nj, 0.0);
        const double h_sep = p.L_sep / static_cast<double>(static_cast<float>(p.sep_node - 2));
        const double h_cat = L_cath / static_cast<double>(static_cast<float>(nj - p.sep_node - 1));
        for (int j = 1; j < s; ++j) dx[j] = h_sep;
        for (int j = s + 1; j < nj - 1; ++j) dx[j] = h_cat;
        faces(dx, aW, aE, bW, bE);
        const double h_c = p.R_agg / static_cast<double>(static_cast<float>(na - 2));
        xa.assign(na, 0.0); dxa.assign(na, 0.0);
        for (int j = 1; j < na - 1; ++j) { xa[j] = h_c * static_cast<double>(static_cast<float>(j)) - h_c / 2.0; dxa[j] = h_c; }
        xa[na - 1] = p.R_agg;
        faces(dxa, aWa, aEa, bWa, bEa);
        const float pf = static_cast<float>(PIr);
        pi4 = static_cast<double>(4.0f * pf);
        pi43 = static_cast<double>((4.0f * pf) / 3.0f);
        pi43b = static_cast<double>((4.0f / 3.0f) * pf);
        AWs.assign(na, 0.0); AEs.assign(na, 0.0); dV.assign(na, 0.0);
        for (int j = 0; j < na; ++j) {
            const double rW = xa[j] - dxa[j] / 2.0, rE = xa[j] + dxa[j] / 2.0;
            AWs[j] = pi4 * std::pow(rW, 2.0);
            AEs[j] = pi4 * std::pow(rE, 2.0);
            dV[j] = pi43 * (std::pow(rE, 3.0) - std::pow(rW, 3.0));
        }
        vol_agg = pi43b * powi(p.R_agg, 3);
        area_agg = pi4 * powi(p.R_agg, 2);
        ce.assign(static_cast<std::size_t>(nj) * NV, 0.0); dce.assign(ce.size(), 0.0); ie.assign(nj, 0.0);
        ca.assign(static_cast<std::size_t>(nj) * na * NV, 0.0);
        for (int j = 0; j < nj; ++j) {
            ce[j * NV + IC] = p.c_bulk; ce[j * NV + IP1] = p.phi1_init; ce[j * NV + IP2] = 0.0; ce[j * NV + ICS] = p.cs_init;
            for (int r = 0; r < na; ++r) {
                double* x = &ca[(static_cast<std::size_t>(j) * na + r) * NV];
                x[IC] = p.c0_init; x[ICS] = p.cs_init; x[IP1] = p.phi1_init; x[IP2] = 0.0;
            }
        }
        Ae.assign(static_cast<std::size_t>(nj) * NV * NV, 0.0); Be = Ae; De = Ae; Ge.assign(static_cast<std::size_t>(nj) * NV, 0.0);
        Aa.assign(static_cast<std::size_t>(na) * NV * NV, 0.0); Ba = Aa; Da = Aa; Ga.assign(static_cast<std::size_t>(na) * NV, 0.0);
    }

    double* C(int j) { return &ce[static_cast<std::size_t>(j) * NV]; }
    double* CA(int l, int r) { return &ca[(static_cast<std::size_t>(l) * na + r) * NV]; }

    // ---- kinetics ----
    double ocp(double cc, double cs) const {
        const double xish_max = p.M * p.Q_th * 3600 / p.F;
        const double th = (cs / p.mol_vol) / xish_max;
        float vint = 0.0f;
        for (int k = 0; k <= 11; ++k)
            vint = static_cast<float>(static_cast<double>(vint) +
                                      r32(AK_AGG[k]) * (powi(2 * th - 1, k + 1) - (2 * th * k * (1 - th)) / powi(2 * th - 1, 1 - k)));
        return r32(3.8637058886774844) + p.R * p.T / p.F * std::log(cc / p.c_bulk * (1.0 - th) / th) + static_cast<double>(vint);
    }
    double ex(double cc, double cs) const {
        const double cimax = p.mol_vol * p.M * p.Q_th * 3600 / p.F;
        return p.F * k_rxn * std::pow(cc, p.alpha_a) * std::pow(cimax - cs, p.alpha_a) * std::pow(cs, p.alpha_c);
    }
    double rate(double cc, double cs, double p1, double p2) const {
        const double eta = p1 - p2 - ocp(cc, cs);
        const double e = r32(ex(cc, cs));   // ex_curr is REAL (D-4)
        return e * (std::exp(p.alpha_a * p.F * eta / (p.R * p.T)) - std::exp(-(p.alpha_c * p.F * eta / (p.R * p.T))));
    }
    void rate_fd(double cc, double cs, double p1, double p2, double& i, double di[NV]) const {
        const double h = p.fd_step;
        i = rate(cc, cs, p1, p2);
        di[IC] = cc <= h ? (rate(cc + h, cs, p1, p2) - i) / h : (rate(cc + h, cs, p1, p2) - rate(cc - h, cs, p1, p2)) / (2.0 * h);
        di[ICS] = cs <= h ? (rate(cc, cs + h, p1, p2) - i) / h : (rate(cc, cs + h, p1, p2) - rate(cc, cs - h, p1, p2)) / (2.0 * h);
        di[IP1] = (rate(cc, cs, p1 + h, p2) - rate(cc, cs, p1 - h, p2)) / (2.0 * h);
        di[IP2] = (rate(cc, cs, p1, p2 + h) - rate(cc, cs, p1, p2 - h)) / (2.0 * h);
    }

    // control-volume coefficients of node j -> blocks (A dc[j-1] + B dc[j] + D dc[j+1] = G)
    static void blocks(int j, int n, const Mat& rj, const Mat& dW, const Mat& dE, const Mat& fW, const Mat& fE,
                       double aw, double ae, double bw, double be, double* A, double* B, double* D) {
        for (int r = 0; r < NV; ++r)
            for (int k = 0; k < NV; ++k) {
                const int q = r * NV + k;
                A[q] = 0.0; D[q] = 0.0;
                if (j == 0) {
                    B[q] = rj[r][k] - (1.0 - ae) * fE[r][k] + be * dE[r][k];
                    D[q] = -(ae * fE[r][k]) - be * dE[r][k];
                } else if (j == n - 1) {
                    A[q] = (1.0 - aw) * fW[r][k] - bw * dW[r][k];
                    B[q] = rj[r][k] + bw * dW[r][k] + aw * fW[r][k];
                } else {
                    A[q] = (1.0 - aw) * fW[r][k] - bw * dW[r][k];
                    B[q] = rj[r][k] + bw * dW[r][k] + aw * fW[r][k] - (1.0 - ae) * fE[r][k] + be * dE[r][k];
                    D[q] = -(ae * fE[r][k]) - be * dE[r][k];
                }
            }
    }

    // electrode scale (the kinetics use the agglomerate-surface crystal concentration)
    void electrode(double h, double I) {
        const double F = p.F, F2 = F * F, eps = p.eps, es = p.eps_sep, sig = p.sigma, a = a_e;
        const double Ds = p.D / p.tau_sep, us = u0 / p.tau_sep, ksep = us + us, ksep_b = u0 + u0, kcat = ucat + uan;
        const int sn = s;
        for (int j = 0; j < nj; ++j) {
            Mat dW{}, dE{}, fW{}, fE{}, rj{};
            double g[NV] = {0, 0, 0, 0}, cW[NV] = {0, 0, 0, 0}, cE[NV] = {0, 0, 0, 0}, gW[NV] = {0, 0, 0, 0}, gE[NV] = {0, 0, 0, 0}, di[NV];
            const double* c0 = C(j);
            for (int k = 0; k < NV; ++k) {
                if (j > 0) { cW[k] = aW[j] * c0[k] + (1.0 - aW[j]) * C(j - 1)[k]; gW[k] = bW[j] * (c0[k] - C(j - 1)[k]); }
                if (j < nj - 1) { cE[k] = aE[j] * C(j + 1)[k] + (1.0 - aE[j]) * c0[k]; gE[k] = bE[j] * (C(j + 1)[k] - c0[k]); }
            }
            double i;
            rate_fd(c0[IC], CA(j, na - 1)[ICS], c0[IP1], c0[IP2], i, di);
            ie[j] = i;
            auto solid = [&]() {
                for (int k = 0; k < NV; ++k) rj[ICS][k] = -(a * di[k] / F);
                rj[ICS][ICS] = -(a * di[ICS] / F) - 1.0 * v_AM / h;
                g[ICS] = +(a * i / F);
            };
            if (j == 0) {
                dE[IC][IC] = -(es * Ds);
                fE[IC][IC] = -(es * 1.0 * us * F * gE[IP2]);
                dE[IC][IP2] = -(es * 1.0 * us * F * cE[IC]);
                g[IC] = -I / F + (dE[IC][IC] * gE[IC] + fE[IC][IC] * cE[IC]);
                rj[ICS][ICS] = 0.0 - 1.0 * (1.0 - es) / h;
                dE[IP1][IP1] = -(1.0 - es) * sig;
                g[IP1] = 0.0 - dE[IP1][IP1] * gE[IP1];  // D-1
                rj[IP2][IP2] = 1.0;
                g[IP2] = 0.0 - c0[IP2];
            } else if (j < sn) {
                dW[IC][IC] = -(es * Ds); dE[IC][IC] = -(es * Ds);
                fW[IC][IC] = -(es * 1.0 * us * F * gW[IP2]); fE[IC][IC] = -(es * 1.0 * us * F * gE[IP2]);
                dW[IC][IP2] = -(es * 1.0 * us * F * cW[IC]); dE[IC][IP2] = -(es * 1.0 * us * F * cE[IC]);
                rj[IC][IC] = -(es / h * dx[j]);
                g[IC] = 0.0 - (fW[IC][IC] * cW[IC] + dW[IC][IC] * gW[IC]) + (fE[IC][IC] * cE[IC] + dE[IC][IC] * gE[IC]);
                rj[ICS][ICS] = -(1.0 * (1.0 - es) / h);
                dW[IP1][IP1] = -(1.0 - es) * sig; dE[IP1][IP1] = -(1.0 - es) * sig;
                g[IP1] = 0.0 - (fW[IP1][IP1] * cW[IP1] + dW[IP1][IP1] * gW[IP1]) + (fE[IP1][IP1] * cE[IP1] + dE[IP1][IP1] * gE[IP1]);
                dW[IP2][IC] = -(es * F * (1.0 * Ds + (-1.0) * Ds)); dE[IP2][IC] = dW[IP2][IC];
                fW[IP2][IC] = -((es / p.tau_sep) * F2 * ksep_b * gW[IP2]); fE[IP2][IC] = -((es / p.tau_sep) * F2 * ksep_b * gE[IP2]);
                dW[IP2][IP2] = -(es * F2 * ksep * cW[IC]); dE[IP2][IP2] = -(es * F2 * ksep * cE[IC]);
                g[IP2] = 0.0 - (fW[IP2][IP2] * cW[IP2] + dW[IP2][IP2] * gW[IP2]) + (fE[IP2][IP2] * cE[IP2] + dE[IP2][IP2] * gE[IP2]);
            } else if (j == sn) {
                dW[IC][IC] = -(es * Ds); dE[IC][IC] = -dcat;
                fW[IC][IC] = -(es * 1.0 * us * F * gW[IP2]); fE[IC][IC] = -(1.0 * ucat * F * gE[IP2]);
                dW[IC][IP2] = -(es * 1.0 * us * F * cW[IC]); dE[IC][IP2] = -(1.0 * ucat * F * cE[IC]);
                g[IC] = 0.0 - (fW[IC][IC] * cW[IC] + dW[IC][IC] * gW[IC]) + (fE[IC][IC] * cE[IC] + dE[IC][IC] * gE[IC]);
                dW[IP1][IP1] = -(1.0 - eps) * sig; dE[IP1][IP1] = -(1.0 - eps) * sig;
                g[IP1] = 0.0 - (fW[IP1][IP1] * cW[IP1] + dW[IP1][IP1] * gW[IP1]) + (fE[IP1][IP1] * cE[IP1] + dE[IP1][IP1] * gE[IP1]);
                dW[IP2][IC] = -(es * F * (1.0 * Ds + (-1.0) * Ds)); dE[IP2][IC] = -(F * (1.0 * dcat + (-1.0) * dan));
                fW[IP2][IC] = -(es * F2 * ksep * gW[IP2]); fE[IP2][IC] = -(F2 * kcat * gE[IP2]);
                dW[IP2][IP2] = -(es * F2 * ksep * cW[IC]); dE[IP2][IP2] = -(F2 * kcat * cE[IC]);
                g[IP2] = 0.0 - (fW[IP2][IP2] * cW[IP2] + dW[IP2][IP2] * gW[IP2]) + (fE[IP2][IP2] * cE[IP2] + dE[IP2][IP2] * gE[IP2]);
                solid();
            } else if (j < nj - 1) {
                dW[IC][IC] = -dcat; dE[IC][IC] = -dcat;
                fW[IC][IC] = -(1.0 * ucat * F * gW[IP2]); fE[IC][IC] = -(1.0 * ucat * F * gE[IP2]);
                dW[IC][IP2] = -(1.0 * ucat * F * cW[IC]); dE[IC][IP2] = -(1.0 * ucat * F * cE[IC]);
                for (int k = 0; k < NV; ++k) {
                    rj[IC][k] = (a * di[k] / F) * dx[j];
                    rj[IP1][k] = -((a * di[k]) * dx[j]);
                    rj[IP2][k] = (a * di[k]) * dx[j];
                }
                rj[IC][IC] = (a * di[IC] / F) * dx[j] - (eps / h) * dx[j];
                g[IC] = -((a * i / F) * dx[j]) - (fW[IC][IC] * cW[IC] + dW[IC][IC] * gW[IC]) + (fE[IC][IC] * cE[IC] + dE[IC][IC] * gE[IC]);
                dW[IP1][IP1] = -(1.0 - eps) * sig; dE[IP1][IP1] = -(1.0 - eps) * sig;
                g[IP1] = (a * i) * dx[j] - (fW[IP1][IP1] * cW[IP1] + dW[IP1][IP1] * gW[IP1]) + (fE[IP1][IP1] * cE[IP1] + dE[IP1][IP1] * gE[IP1]);
                dW[IP2][IC] = -(F * (1.0 * dcat + (-1.0) * dan)); dE[IP2][IC] = dW[IP2][IC];
                fW[IP2][IC] = -(F2 * kcat * gW[IP2]); fE[IP2][IC] = -(F2 * kcat * gE[IP2]);
                dW[IP2][IP2] = -(F2 * kcat * cW[IC]); dE[IP2][IP2] = -(F2 * kcat * cE[IC]);
                g[IP2] = -((a * i) * dx[j]) - (fW[IP2][IP2] * cW[IP2] + dW[IP2][IP2] * gW[IP2]) + (fE[IP2][IP2] * cE[IP2] + dE[IP2][IP2] * gE[IP2]);  // D-11
                solid();
            } else {
                dW[IC][IC] = -dcat;
                fW[IC][IC] = -(1.0 * ucat * F * gW[IP2]);
                dW[IC][IP2] = -(1.0 * ucat * F * cW[IC]);
                g[IC] = 0.0 - dW[IC][IC] * gW[IC] - fW[IC][IC] * cW[IC];
                dW[IP1][IP1] = -(1.0 - eps) * sig;
                g[IP1] = I - dW[IP1][IP1] * gW[IP1];
                dW[IP2][IC] = -(F * (1.0 * dcat + (-1.0) * dan));
                fW[IP2][IC] = -(F2 * kcat * gW[IP2]);
                dW[IP2][IP2] = -(F2 * kcat * cW[IC]);
                g[IP2] = 0.0 - dW[IP2][IC] * gW[IC] - fW[IP2][IC] * cW[IC];
                solid();
            }
            const std::size_t o = static_cast<std::size_t>(j) * NV * NV;
            blocks(j, nj, rj, dW, dE, fW, fE, aW[j], aE[j], bW[j], bE[j], &Ae[o], &Be[o], &De[o]);
            for (int k = 0; k < NV; ++k) Ge[static_cast<std::size_t>(j) * NV + k] = g[k];
        }
    }

    // one agglomerate (at electrode node l): spherical control volumes
    void particle(int l, double h, double dcs_dt) {
        const double F = p.F, F2 = F * F, ea = p.eps_agg, sig = p.sigma, a = a_agg, kagg = uagg + uagg;
        const double c_spec = dcs_dt * F / p.rho;
        const double i_agg = c_spec * vol_agg * (1 - ea) * p.rho / area_agg;
        for (int j = 0; j < na; ++j) {
            Mat dW{}, dE{}, fW{}, fE{}, rj{};
            double g[NV] = {0, 0, 0, 0}, cW[NV] = {0, 0, 0, 0}, cE[NV] = {0, 0, 0, 0}, gW[NV] = {0, 0, 0, 0}, gE[NV] = {0, 0, 0, 0}, di[NV];
            const double* x = CA(l, j);
            for (int k = 0; k < NV; ++k) {
                if (j > 0) { cW[k] = aWa[j] * x[k] + (1.0 - aWa[j]) * CA(l, j - 1)[k]; gW[k] = bWa[j] * (x[k] - CA(l, j - 1)[k]); }
                if (j < na - 1) { cE[k] = aEa[j] * CA(l, j + 1)[k] + (1.0 - aEa[j]) * x[k]; gE[k] = bEa[j] * (CA(l, j + 1)[k] - x[k]); }
            }
            double i;
            rate_fd(x[IC], x[ICS], x[IP1], x[IP2], i, di);
            for (int k = 0; k < NV; ++k) rj[ICS][k] = -(a * di[k] / F);
            g[ICS] = +(a * i / F);
            if (j == 0) {
                dE[IC][IC] = 1.0;
                g[IC] = 0.0 - dE[IC][IC] * gE[IC];  // D-1 (center)
                dE[IP1][IP1] = -(1.0 - ea) * sig;
                g[IP1] = 0.0 - dE[IP1][IP1] * gE[IP1];
                dE[IP2][IP2] = 1.0;
                g[IP2] = 0.0 - dE[IP2][IP2] * gE[IP2];
                rj[ICS][ICS] = -(a * di[ICS] / F) - 1.0 * (1.0 - ea) / h;
            } else if (j < na - 1) {
                dW[IC][IC] = AWs[j] * (-(ea * D_agg)); dE[IC][IC] = AEs[j] * (-(ea * D_agg));
                for (int k = 0; k < NV; ++k) {
                    rj[IC][k] = (a * di[k] / F) * dV[j];
                    rj[IP1][k] = -((a * di[k]) * dV[j]);
                    rj[IP2][k] = (a * di[k]) * dV[j];
                }
                rj[IC][IC] = (a * di[IC] / F) * dV[j] - (ea / h) * dV[j];
                g[IC] = -((a * i / F) * dV[j]) - (fW[IC][IC] * cW[IC] + dW[IC][IC] * gW[IC]) + (fE[IC][IC] * cE[IC] + dE[IC][IC] * gE[IC]);
                rj[ICS][ICS] = -(a * di[ICS] / F) - (1.0 - ea) / h;
                dW[IP1][IP1] = AWs[j] * (-(1.0 - ea) * sig); dE[IP1][IP1] = AEs[j] * (-(1.0 - ea) * sig);
                g[IP1] = (a * i) * dV[j] - (fW[IP1][IP1] * cW[IP1] + dW[IP1][IP1] * gW[IP1]) + (fE[IP1][IP1] * cE[IP1] + dE[IP1][IP1] * gE[IP1]);
                dW[IP2][IC] = AWs[j] * (-(ea * F * (1.0 * D_agg + (-1.0) * D_agg)));
                dE[IP2][IC] = AEs[j] * (-(ea * F * (1.0 * D_agg + (-1.0) * D_agg)));
                fW[IP2][IC] = AWs[j] * (-(ea * F2 * kagg * gW[IP2])); fE[IP2][IC] = AEs[j] * (-(ea * F2 * kagg * gE[IP2]));
                dW[IP2][IP2] = AWs[j] * (-(ea * F2 * kagg * cW[IC])); dE[IP2][IP2] = AEs[j] * (-(ea * F2 * kagg * cE[IC]));
                fW[IP2][IP2] = AWs[j] * 0.0; fE[IP2][IP2] = AEs[j] * 0.0;
                g[IP2] = -((a * i) * dV[j]) - (fW[IP2][IP2] * cW[IP2] + dW[IP2][IP2] * gW[IP2]) + (fE[IP2][IP2] * cE[IP2] + dE[IP2][IP2] * gE[IP2]);
            } else {
                rj[IC][IC] = 1.0;
                g[IC] = p.c_bulk - x[IC];  // D-17
                rj[ICS][ICS] = -(a * di[ICS] / F) - (1.0 - ea) / h;
                dW[IP1][IP1] = -(1.0 - ea) * sig;
                g[IP1] = i_agg - dW[IP1][IP1] * gW[IP1];
                rj[IP2][IP2] = 1.0;
                g[IP2] = C(l)[IP2] - x[IP2];
            }
            const std::size_t o = static_cast<std::size_t>(j) * NV * NV;
            blocks(j, na, rj, dW, dE, fW, fE, aWa[j], aEa[j], bWa[j], bEa[j], &Aa[o], &Ba[o], &Da[o]);
            for (int k = 0; k < NV; ++k) Ga[static_cast<std::size_t>(j) * NV + k] = g[k];
        }
    }

    int run() {
        std::FILE* out = std::fopen(p.file.c_str(), "w");
        if (!out) throw std::runtime_error("cannot write " + p.file);
        double t = 0.0, mAhg = 0.0, dt = dt_nom, ramp = 1.0, current = 0.0;
        const double write_every = dt_nom, to_electrons = 1.0 / p.rho * p.M;
        long last_write = 0;  // uninitialized in the original (D-19); 0 as on the HPC cluster
        const char state = 'D';
        std::string exit_reason = "max_steps";
        long long nsolve = 0;
        auto write_row = [&](bool header) {
            if (header) {
                const char* h1[] = {"State", "Time", "Voltage", "mAhg", "Equivalence", "Solid_Conc", "current_density", "iloc",
                                    "Solution_Pot", "c0", "cs_edge", "cs", "eta_contact", "U", "eta_rxn", "i0", "current", "ramp_current"};
                const char* h2[] = {"CDR", "hours", "Volts", "mAh/g", "LixNMC", "LixNMC", "A/cm2", "A/cm2", "Volts", "mol/cm3",
                                    "mol/cm3", "mol/cm3", "Volts", "Volts", "Volts", "A/cm2", "A/cm2", ","};
                for (auto* h : {h1, h2}) {
                    std::string l = fit(h[0], 5) + " " + fit(h[1], 12) + " " + fit(h[2], 12);
                    for (int k = 3; k < 18; ++k) l += " " + fit(h[k], 15);
                    std::fprintf(out, "%s\n", l.c_str());
                }
            }
            const double p1 = C(nj - 1)[IP1], p2 = C(nj - 1)[IP2], c0 = C(nj - 1)[IC], cs = CA(nj - 1, na - 1)[ICS];
            const double u = ocp(c0, cs), eta_contact = current * 0.0;
            const double vals[16] = {mAhg, mAhg * p.M * lit36 / p.F, cs * to_electrons, i_final, rate(c0, cs, p1, p2), p2, c0, cs,
                                     C(nj - 1)[ICS], eta_contact, u, p1 - p2 - u, ex(c0, cs), current, ramp, 0.0};
            std::string l = std::string(4, ' ') + state + " " + fixed12(t / 3600.0) + " " + fixed12(p1 - eta_contact);
            for (int k = 0; k < 15; ++k) l += " " + sci15(vals[k]);
            std::fprintf(out, "%s\n", l.c_str());
        };
        auto solve_faithful = [&](int n, std::vector<double>& A, std::vector<double>& B, std::vector<double>& D,
                                  std::vector<double>& G, std::vector<double>& x) {
            if (!band::solve(NV, n, A, B, D, G, x, band::Pivot::legacy, true))
                std::fill(x.begin(), x.end(), std::numeric_limits<double>::quiet_NaN());
        };
        for (long long it = 1; it <= nsteps; ++it) {
            if (it == 1) write_row(true);
            else if ((t - last_write) >= write_every) { write_row(false); last_write = static_cast<long>(t - dt); }
            else if (it >= nsteps) write_row(false);
            if (C(nj - 1)[IP1] >= 99.0 && state == 'C') { write_row(false); exit_reason = "end_of_charge"; break; }
            else if (CA(s, na - 1)[ICS] * to_electrons >= thr_x) { write_row(false); exit_reason = "x_limit"; break; }
            else if (std::isnan(dce[IC])) { write_row(false); exit_reason = "nan"; break; }
            else if (t >= 99.0 * 3600.0) { write_row(false); exit_reason = "max_time"; break; }
            if (C(nj - 1)[IC] <= thr_dep) dt = static_cast<double>(static_cast<float>(tmx) / (static_cast<float>(nsteps) * 10.0f));
            else if (C(nj - 1)[IC] <= thr_dep2) dt = static_cast<double>(static_cast<float>(tmx) / (static_cast<float>(nsteps) * 100.0f));
            else dt = dt_nom;
            t = t + dt;
            if (ramp == 1.0) { current = i_final / 50.0; ramp = 2.0; }
            else if (ramp == 2.0 && std::abs(current * 1.5) < std::abs(i_final)) current = current * 1.5;
            else { current = i_final; ramp = 0.0; }
            electrode(dt, current);
            mAhg = mAhg + 1000.0 * i_specific * dt / 3600.0;
            solve_faithful(nj, Ae, Be, De, Ge, dce);
            for (std::size_t k = 0; k < ce.size(); ++k) ce[k] = ce[k] + dce[k];
            for (int l = s; l < nj; ++l) {
                particle(l, dt, dce[static_cast<std::size_t>(l) * NV + ICS] / dt);
                solve_faithful(na, Aa, Ba, Da, Ga, dca);
                double* x = CA(l, 0);
                for (std::size_t k = 0; k < dca.size(); ++k) x[k] = x[k] + dca[k];
            }
            nsolve = it;
        }
        std::fclose(out);
        std::printf("faithful agglomerate run, C-rate %.17g: exit %s after %lld steps; wrote %s\n", p.C_rate, exit_reason.c_str(),
                    nsolve, p.file.c_str());
        return 0;
    }
};

}  // namespace

// ============================== agglomerate model, corrected mode ==============================
// Double porosity with one reaction description (docs/model.md sections 6 and 8): the electrode scale
// (macro-pores, no active material of its own) and porous agglomerates at each cathode volume, which
// draw salt, ionic and electronic current through their surfaces. Log variables (u, phi1, phi2, s) and
// Scharfetter-Gummel fluxes, as in Model::log_assemble. Each Newton iteration condenses the
// agglomerates onto the electrode's diagonal blocks (a fully coupled step). The protocol driver in
// main() is shared with the uniform model.
struct AggCorrected {
    Model& m;
    const Params& p;
    int nj, s, na, nl;
    double v_AM, v_agg, s_agg, a_x, x_max_a, csmax_a, sig_a;
    std::vector<double> xa, dxa, dV, a_area, a_g, a_gs;
    std::vector<double> Aq, Bq, Dq, Gq;   // stacked agglomerate blocks

    explicit AggCorrected(Model& mm) : m(mm), p(mm.p) {
        nj = p.nj; s = m.s; na = p.nja; nl = nj - s - 2;
        const double L_cath = m.L_cath;
        m.mass_area = p.percent_active * p.mass_loading;  // = L_cath*v_AM*rho
        m.i_1C = p.Q_th * m.mass_area;
        v_AM = p.percent_active * p.mass_loading / (p.rho * L_cath);
        v_agg = v_AM / (1.0 - p.eps_agg);
        s_agg = 3.0 * v_agg / p.R_agg;
        a_x = 3.0 * (1.0 - p.eps_agg) / p.R_xtal;
        if (p.eps + v_agg > 1.0 + 1e-12) throw std::runtime_error("porosity + agglomerate volume fraction exceed 1");
        x_max_a = p.M * p.Q_th * 3600.0 / p.F;
        csmax_a = p.mol_vol * x_max_a;
        const double tau_a = p.tau_agg < 0 ? std::pow(p.eps_agg, -0.5) : p.tau_agg;
        sig_a = p.sigma_agg < 0 ? p.sigma : p.sigma_agg;
        // radial grid (zero-volume center and surface nodes)
        xa.assign(na, 0.0); dxa.assign(na, 0.0); dV.assign(na, 0.0);
        a_area.assign(na - 1, 0.0); a_g.assign(na - 1, 0.0); a_gs.assign(na - 1, 0.0);
        const double hr = p.R_agg / static_cast<double>(na - 2), pi = std::acos(-1.0);
        for (int j = 1; j < na - 1; ++j) { xa[j] = hr * static_cast<double>(j) - hr / 2.0; dxa[j] = hr; }
        xa[na - 1] = p.R_agg;
        for (int j = 0; j < na; ++j) {
            double rW = xa[j] - dxa[j] / 2.0, rE = xa[j] + dxa[j] / 2.0;
            if (j == 0) rW = rE = 0.0;
            if (j == na - 1) rW = rE = p.R_agg;
            dV[j] = 4.0 * pi / 3.0 * (rE * rE * rE - rW * rW * rW);
        }
        for (int j = 0; j < na - 1; ++j) {  // face j joins nodes j and j+1
            const double rf = j == 0 ? 0.0 : xa[j] + dxa[j] / 2.0, h = (dxa[j] + dxa[j + 1]) / 2.0;
            a_area[j] = 4.0 * pi * rf * rf;
            a_g[j] = p.eps_agg / tau_a / h;
            a_gs[j] = (1.0 - p.eps_agg) * sig_a / h;
        }
        const std::size_t nq = static_cast<std::size_t>(nl) * na;
        Aq.assign(nq * NV * NV, 0.0); Bq.assign(Aq.size(), 0.0); Dq.assign(Aq.size(), 0.0); Gq.assign(nq * NV, 0.0);
    }

    std::vector<double> initial() const {
        std::vector<double> ca(static_cast<std::size_t>(nl) * na * NV);
        const double th0 = p.cs_init / csmax_a;
        for (std::size_t q = 0; q < ca.size() / NV; ++q) {
            ca[q * NV + IC] = std::log(p.c0_init / p.c_bulk); ca[q * NV + IP1] = p.phi1_init; ca[q * NV + IP2] = 0.0;
            ca[q * NV + ICS] = std::log(th0 / (1.0 - th0));
        }
        return ca;
    }

    // stacked blocks and G = -R of every agglomerate; the surface takes the electrode's u, phi1, phi2
    void assemble(const std::vector<double>& c, const std::vector<double>& ca, const std::vector<double>& ca_old, double h) {
        const double F = p.F, ea = p.eps_agg;
        std::vector<double> Fv(3 * static_cast<std::size_t>(na - 1)), dFa(3 * NV * static_cast<std::size_t>(na - 1)),
            dFb(dFa.size());
        for (int l = 0; l < nl; ++l) {
            const double* xl = &ca[static_cast<std::size_t>(l) * na * NV];
            for (int k = 0; k < na - 1; ++k)
                m.sg_face(xl + k * NV, xl + (k + 1) * NV, a_g[k], a_gs[k], &Fv[3 * k], &dFa[3 * NV * k], &dFb[3 * NV * k]);
            auto FA = [&](int k, int f, int v) { return dFa[3 * NV * k + f * NV + v]; };
            auto FB = [&](int k, int f, int v) { return dFb[3 * NV * k + f * NV + v]; };
            for (int j = 0; j < na; ++j) {
                const std::size_t q = static_cast<std::size_t>(l) * na + j, o = q * NV * NV;
                double Rr[NV] = {0, 0, 0, 0}, i, di[NV];
                for (int z = 0; z < NV * NV; ++z) { Aq[o + z] = 0.0; Bq[o + z] = 0.0; Dq[o + z] = 0.0; }
                auto At = [&](int r, int v) -> double& { return Aq[o + r * NV + v]; };
                auto Bt = [&](int r, int v) -> double& { return Bq[o + r * NV + v]; };
                auto Dt = [&](int r, int v) -> double& { return Dq[o + r * NV + v]; };
                const double* x = xl + j * NV;
                m.log_rate(x, true, csmax_a, i, di);
                if (j == 0) {  // center: zero gradient
                    for (int col = IC; col <= IP2; ++col) { Rr[col] = xl[NV + col] - x[col]; Bt(col, col) = -1.0; Dt(col, col) = 1.0; }
                } else if (j == na - 1) {  // surface: the electrode's u, phi1, phi2
                    for (int col = IC; col <= IP2; ++col) {
                        Rr[col] = x[col] - c[static_cast<std::size_t>(s + 1 + l) * NV + col];
                        Bt(col, col) = 1.0;
                    }
                } else {
                    const double V = dV[j];
                    for (int f = 0; f < 3; ++f) {
                        Rr[f] = a_area[j] * Fv[3 * j + f] - a_area[j - 1] * Fv[3 * (j - 1) + f];
                        for (int v = 0; v < NV; ++v) {
                            Bt(f, v) = a_area[j] * FA(j, f, v) - a_area[j - 1] * FB(j - 1, f, v);
                            Dt(f, v) = a_area[j] * FB(j, f, v);
                            At(f, v) = -a_area[j - 1] * FA(j - 1, f, v);
                        }
                    }
                    const double cj = p.c_bulk * std::exp(x[IC]), cjo = p.c_bulk * std::exp(ca_old[q * NV + IC]);
                    Rr[IC] += ea * V * (cj - cjo) / h - a_x * i * V / F;
                    Bt(IC, IC) += ea * V * cj / h;
                    for (int v = 0; v < NV; ++v) {
                        Bt(IC, v) -= a_x * di[v] * V / F;
                        Bt(IP1, v) += a_x * di[v] * V;
                        Bt(IP2, v) -= a_x * di[v] * V;
                    }
                    Rr[IP1] += a_x * i * V;
                    Rr[IP2] -= a_x * i * V;
                }
                // crystals, every node: (1 - eps_agg) dcs/dt = -a i/F
                const double th = sigm(x[ICS]), tho = sigm(ca_old[q * NV + ICS]);
                Rr[ICS] = (1.0 - ea) * csmax_a * (th - tho) / h + a_x * i / F;
                for (int v = 0; v < NV; ++v) Bt(ICS, v) = a_x * di[v] / F;
                Bt(ICS, ICS) += (1.0 - ea) * csmax_a * th * (1 - th) / h;
                for (int r = 0; r < NV; ++r) Gq[q * NV + r] = -Rr[r];
            }
        }
    }

    // flux into agglomerate l through r = R: q = (N+, i1, i2), and dq/dx at its last two nodes
    void surface(int l, const std::vector<double>& ca, double qv[3], double Jq[3][2 * NV]) const {
        const double* xin = &ca[(static_cast<std::size_t>(l) * na + na - 2) * NV];
        double Fv[3], dFa[3][NV], dFb[3][NV];
        m.sg_face(xin, xin + NV, a_g[na - 2], a_gs[na - 2], Fv, &dFa[0][0], &dFb[0][0]);
        for (int r = 0; r < 3; ++r) {
            qv[r] = -Fv[r];
            for (int v = 0; v < NV; ++v) { Jq[r][v] = -dFa[r][v]; Jq[r][NV + v] = -dFb[r][v]; }
        }
    }

    // one backward-Euler step of length h for both scales (condensed Newton); c, ca updated on success
    bool newton(std::vector<double>& c, std::vector<double>& ca, double h, std::vector<double>& A, std::vector<double>& B,
                std::vector<double>& Dm, std::vector<double>& G) {
        const std::vector<double> c_old = c, ca_old = ca;
        const std::size_t nq = static_cast<std::size_t>(nl) * na;
        std::vector<double> x0, Zc[3], E(nq * NV), dce, dca(nq * NV), rsurf(3 * static_cast<std::size_t>(nl));
        double prev = std::numeric_limits<double>::infinity();
        for (int it = 0; it < p.newton_max_iter; ++it) {
            // agglomerates: the update at fixed surface values and the three surface responses
            assemble(c, ca, ca_old, h);
            for (std::size_t q = 0; q < nq; ++q)
                for (int r = 0; r < NV; ++r) {
                    const std::size_t o = (q * NV + r) * NV;
                    double sc = 0.0;
                    for (int k = 0; k < NV; ++k) sc = std::max(sc, std::abs(Bq[o + k]));
                    if (sc == 0.0) sc = 1.0;
                    for (int k = 0; k < NV; ++k) { Aq[o + k] /= sc; Bq[o + k] /= sc; Dq[o + k] /= sc; }
                    Gq[q * NV + r] /= sc;
                    if (static_cast<int>(q % na) == na - 1 && r <= IP2) rsurf[3 * (q / na) + r] = 1.0 / sc;
                }
            if (!band::solve(NV, static_cast<int>(nq), Aq, Bq, Dq, Gq, x0, band::Pivot::partial)) break;
            bool good = true;
            for (int col = 0; col < 3 && good; ++col) {
                std::fill(E.begin(), E.end(), 0.0);
                for (int l = 0; l < nl; ++l) E[(static_cast<std::size_t>(l) * na + na - 1) * NV + col] = rsurf[3 * l + col];
                good = band::solve(NV, static_cast<int>(nq), Aq, Bq, Dq, E, Zc[col], band::Pivot::partial);
            }
            if (!good) break;
            // electrode with the condensed agglomerate sources: R_e + w q
            m.log_assemble(c, c_old, h, m.i_app, false, A, B, Dm, G);
            for (int l = 0; l < nl; ++l) {
                const int j = s + 1 + l;
                const double w = s_agg * m.dx[j];
                double qv[3], Jq[3][2 * NV];
                surface(l, ca, qv, Jq);
                const std::size_t b0 = (static_cast<std::size_t>(l) * na + na - 2) * NV;
                for (int r = 0; r < 3; ++r) {
                    double jx = 0.0;
                    for (int k = 0; k < 2 * NV; ++k) jx += Jq[r][k] * x0[b0 + k];
                    G[static_cast<std::size_t>(j) * NV + r] += -w * qv[r] - w * jx;
                    for (int col = 0; col < 3; ++col) {
                        double jz = 0.0;
                        for (int k = 0; k < 2 * NV; ++k) jz += Jq[r][k] * Zc[col][b0 + k];
                        B[(static_cast<std::size_t>(j) * NV + r) * NV + col] += w * jz;
                    }
                }
            }
            equilibrate_blocks(nj, A, B, Dm, G);
            if (!band::solve(NV, nj, A, B, Dm, G, dce, band::Pivot::partial)) break;
            for (int l = 0; l < nl; ++l) {
                const std::size_t je = static_cast<std::size_t>(s + 1 + l) * NV;
                for (int k = 0; k < na; ++k) {
                    const std::size_t q = (static_cast<std::size_t>(l) * na + k) * NV;
                    for (int v = 0; v < NV; ++v)
                        dca[q + v] = x0[q + v] + Zc[0][q + v] * dce[je + IC] + Zc[1][q + v] * dce[je + IP1] + Zc[2][q + v] * dce[je + IP2];
                }
            }
            const double lam = std::min(lbounded(dce), lbounded(dca));
            double raw = 0.0;
            for (double v : dce) raw = std::max(raw, std::abs(v));
            for (double v : dca) raw = std::max(raw, std::abs(v));
            for (std::size_t q = 0; q < c.size(); ++q) c[q] = c[q] + lam * dce[q];
            for (std::size_t q = 0; q < ca.size(); ++q) ca[q] = ca[q] + lam * dca[q];
            const double upd = std::max(phys_update(c, dce, true), phys_update(ca, dca, false));
            if (!std::isfinite(raw) || raw > 1.0e3) break;
            if (lam == 1.0 && converged(upd, prev, p.newton_tol)) return true;
            prev = upd;
        }
        c = c_old;
        ca = ca_old;
        return false;
    }
};

int main(int argc, char** argv) {
    try {
        Params input = read_input(argc > 1 ? argv[1] : "nmc.nml");
        if (input.particle_model == "agglomerate" && input.mode == "faithful") return Agglomerate(input).run();
        if (input.particle_model != "uniform" && input.particle_model != "agglomerate")
            throw std::runtime_error("particle_model must be 'uniform' or 'agglomerate'");
        const bool aggc = input.particle_model == "agglomerate";
        if (aggc) {
            input.eps_AM = 0.0;  // the electrode scale has no active material of its own
            if (input.k_rxn < 0) input.k_rxn = 2.5e-6;
        }
        Model m(input);
        std::unique_ptr<AggCorrected> agg;
        if (aggc) agg = std::make_unique<AggCorrected>(m);
        std::vector<double> ca;
        if (aggc) ca = agg->initial();
        const Params& p = m.p;
        const int nj = p.nj;
        std::vector<Step> steps;
        if (!m.faithful) steps = parse_protocol(p.steps, p.C_rate, p.V_min, p.V_max, p.cycles);
        std::FILE* out = std::fopen(p.file.c_str(), "w");
        if (!out) throw std::runtime_error("cannot write " + p.file);

        std::vector<double> c(static_cast<std::size_t>(nj) * NV), dc(c.size(), 0.0), A(c.size() * NV), B(A.size()),
            Dm(A.size()), G(c.size());
        for (int j = 0; j < nj; ++j) {
            c[j * NV + IC] = p.c_bulk; c[j * NV + IP1] = p.phi1_init; c[j * NV + IP2] = p.phi2_init; c[j * NV + ICS] = p.cs_init;
        }
        double t = 0.0, mAhg = 0.0, dt = p.t_max / static_cast<double>(p.n_steps);
        if (aggc) dt = p.dt_s;
        if (!m.faithful) {  // corrected mode: u = ln(c/c_bulk) and the particles' log-odds
            const double th0 = p.cs_init / m.cs_max();
            for (int j = 0; j < nj; ++j) {
                c[j * NV + IC] = 0.0;
                if (!aggc) c[j * NV + ICS] = std::log(th0 / (1.0 - th0));
            }
        }
        std::string exit_reason = "max_steps";
        int nsolve = 0;

        auto header_lines = [&](bool extended) {
            std::vector<std::string> h1 = {"State", "Time", "Voltage", "Equivalence", "Anode_Eta", "anode_exchange_c", "Edge_c0"};
            std::vector<std::string> h2 = {"CDR", "hours", "Volts", "electron_equivs", "mV", "mA/cm2", "mol/cm3"};
            if (extended) {
                for (const char* x : {"Current", "Step", "Li_Nernst"}) h1.push_back(x);
                for (const char* x : {"mA/cm2", "#", "mV"}) h2.push_back(x);
                if (aggc) { h1.push_back("x_front"); h1.push_back("c_collector"); h2.push_back("LixNMC"); h2.push_back("mol/cm3"); }
            }
            for (auto* h : {&h1, &h2}) {
                std::string l = fit((*h)[0], 5) + " " + fit((*h)[1], 12) + " " + fit((*h)[2], 12);
                for (std::size_t k = 3; k < h->size(); ++k) l += " " + fit((*h)[k], 15);
                std::fprintf(out, "%s\n", l.c_str());
            }
        };
        // electrolyte concentration at the foil (corrected mode stores u = ln(c/c_bulk))
        auto c_foil = [&]() { return m.faithful ? c[IC] : p.c_bulk * std::exp(c[IC]); };
        auto i0_li = [&]() { return p.F * p.k_Li * std::pow(c_foil(), 0.5) * std::pow(p.c_Li_ref, 0.5); };
        // Li counter-electrode overpotential; corrected mode: symmetric Butler-Volmer (D-12)
        auto li_eta = [&](char state) {
            const double alpha = 0.5;
            if (!m.faithful) return -(p.R * p.T / (alpha * p.F)) * std::asinh(m.i_app / (2.0 * i0_li()));
            if (state == 'C') return 0.5 * std::log(m.i_app / i0_li()) / (alpha * p.F / (p.R * p.T));
            if (state == 'D') return -(0.5 * std::log(m.i_app / i0_li())) / (alpha * p.F / (p.R * p.T));
            return 0.0;
        };
        // corrected mode: Nernst potential of the lithium foil, (RT/F) ln(c/c_Li_ref)
        auto li_nernst = [&]() { return p.R * p.T / p.F * std::log(c_foil() / p.c_Li_ref); };
        // corrected mode: voltage against the lithium foil (0 V). The solver fixes the gauge with phi2 = 0
        // at the foil face; the equations depend only on potential differences, so the foil-referenced
        // potentials are the solved ones minus U_Li + eta_Li (D-16).
        auto cell_voltage = [&]() { return c[(nj - 1) * NV + IP1] + li_eta('D') - li_nernst(); };

        if (m.faithful) {
            // ---------------- the original program's constant-current discharge ----------------
            const char state = p.C_rate < 0 ? 'C' : 'D';
            const double write_every = p.t_max / p.n_steps / 200;
            long last_write = 0;  // an integer in the original (D-4)
            auto write_row = [&](bool header) {
                if (header) header_lines(false);
                const double eta = li_eta(state);
                std::fprintf(out, "%5c %s %s %s %s %s %s\n", state, fixed12(t / static_cast<double>(3600)).c_str(),
                             fixed12(c[(nj - 1) * NV + IP1] + eta).c_str(), sci15(mAhg * p.M * m.lit36 / p.F).c_str(),
                             sci15(eta * 1.0e3).c_str(), sci15(i0_li() * 1.0e3).c_str(), sci15(c[IC]).c_str());
            };
            for (int it = 1; it <= p.n_steps; ++it) {
                if (it == 1) {
                    write_row(true);
                } else if ((t - last_write) / 3600 >= write_every) {
                    write_row(false);
                    last_write = static_cast<long>(t - dt);
                } else if (it >= p.n_steps) {
                    write_row(false);
                } else if (c[(nj - 1) * NV + IP1] >= 99.0 && state == 'C') {
                    write_row(false); exit_reason = "end_of_charge"; break;
                } else if (std::isnan(dc[IC])) {
                    write_row(false); exit_reason = "nan"; break;
                } else if (t >= 99.0 * 3600.0) {
                    write_row(false); exit_reason = "max_time"; break;
                }
                if (state == 'D') mAhg = mAhg + 1000.0 * m.i_spec * dt / 3600.0;
                else mAhg = mAhg - 1000.0 * m.i_spec * dt / 3600.0;
                m.assemble(c, dt, A, B, Dm, G);
                if (!band::solve(NV, nj, A, B, Dm, G, dc, band::Pivot::legacy))
                    std::fill(dc.begin(), dc.end(), std::numeric_limits<double>::quiet_NaN());
                for (std::size_t k = 0; k < c.size(); ++k) c[k] = c[k] + dc[k];
                nsolve = it;
                dt = p.t_max / static_cast<double>(p.n_steps);
                t = t + dt;
            }
        } else {
            // ---------------- corrected mode: protocol of cc / cv / rest steps ----------------
            auto write_row = [&](bool header, int step) {
                if (header) header_lines(true);
                const char st = m.i_app > 0 ? 'D' : (m.i_app < 0 ? 'C' : 'R');
                const double eta = li_eta(st);
                std::fprintf(out, "%5c %s %s %s %s %s %s %s %15d %s", st, fixed12(t / 3600.0).c_str(),
                             fixed12(cell_voltage()).c_str(), sci15(mAhg * p.M * 3.6 / p.F).c_str(),
                             sci15(eta * 1.0e3).c_str(), sci15(i0_li() * 1.0e3).c_str(), sci15(c_foil()).c_str(),
                             sci15(m.i_app * 1.0e3).c_str(), step, sci15(li_nernst() * 1.0e3).c_str());
                if (aggc)
                    std::fprintf(out, " %s %s",
                                 sci15(agg->csmax_a * sigm(ca[static_cast<std::size_t>(p.nja - 1) * NV + ICS]) / p.mol_vol).c_str(),
                                 sci15(p.c_bulk * std::exp(c[(nj - 1) * NV + IC])).c_str());
                std::fprintf(out, "\n");
            };
            // one backward-Euler step of length h at the current m.i_app, solved with Newton (log variables)
            auto newton_step = [&](double h) {
                if (aggc) return agg->newton(c, ca, h, A, B, Dm, G);
                const std::vector<double> c_old = c;
                double prev = std::numeric_limits<double>::infinity();
                for (int k = 0; k < p.newton_max_iter; ++k) {
                    m.log_assemble(c, c_old, h, m.i_app, true, A, B, Dm, G);
                    equilibrate_blocks(nj, A, B, Dm, G);
                    if (!band::solve(NV, nj, A, B, Dm, G, dc, band::Pivot::partial)) break;
                    const double lam = lbounded(dc);
                    double raw = 0.0;
                    for (double v : dc) raw = std::max(raw, std::abs(v));
                    for (std::size_t q = 0; q < c.size(); ++q) c[q] = c[q] + lam * dc[q];
                    const double upd = phys_update(c, dc, false);
                    if (!std::isfinite(raw) || raw > 1.0e3) break;
                    if (lam == 1.0 && converged(upd, prev, p.newton_tol)) return true;
                    prev = upd;
                }
                c = c_old;
                return false;
            };
            // advance by dt; with check, locate a voltage-cutoff crossing to within 0.1 mV
            auto advance = [&](double dtt, double& t_done, bool& stopped, bool check, double vlo, double vhi) {
                const double min_dt = 1.0e-10, event_dv = 1.0e-4, event_min_dt = 1.0e-12;
                const int max_failures = 200;  // Newton failures allowed within one step
                int failures = 0;
                const std::vector<double> c_begin = c, ca_begin = ca;
                double hh = dtt;
                t_done = 0.0;
                stopped = false;
                while (t_done < dtt) {
                    hh = std::min(hh, dtt - t_done);
                    const std::vector<double> c_save = c, ca_save = ca;
                    if (!newton_step(hh)) {
                        if (hh / 2 < min_dt || ++failures >= max_failures) {
                            c = c_begin;  // give up: report the state at the start of the step
                            ca = ca_begin;
                            return false;
                        }
                        hh = hh / 2;
                        continue;
                    }
                    if (check) {
                        const double vv = cell_voltage(), mg = std::min(vv - vlo, vhi - vv);
                        if (mg < 0.0) {
                            if (mg < -event_dv && hh / 2 >= event_min_dt) { c = c_save; ca = ca_save; hh = hh / 2; continue; }
                            t_done = t_done + hh;
                            stopped = true;
                            return true;
                        }
                    }
                    t_done = t_done + hh;
                    hh = 2.0 * hh;  // grow back after a success (up to dt, by the min above)
                }
                return true;
            };
            // the physical limit the state has reached, reported as the exit reason when a step cannot be
            // solved: electrolyte below 1e-3*c_bulk anywhere, or particles within 1e-3 of full or empty
            auto limit_reason = [&]() -> std::string {
                double cmin = std::numeric_limits<double>::infinity(), thmin = cmin, thmax = -cmin;
                for (int j = 0; j < nj; ++j) cmin = std::min(cmin, p.c_bulk * std::exp(c[j * NV + IC]));
                if (aggc) {
                    for (std::size_t q = 0; q < ca.size() / NV; ++q) {
                        cmin = std::min(cmin, p.c_bulk * std::exp(ca[q * NV + IC]));
                        thmin = std::min(thmin, sigm(ca[q * NV + ICS]));
                        thmax = std::max(thmax, sigm(ca[q * NV + ICS]));
                    }
                } else {
                    for (int j = m.s; j < nj; ++j) {
                        thmin = std::min(thmin, sigm(c[j * NV + ICS]));
                        thmax = std::max(thmax, sigm(c[j * NV + ICS]));
                    }
                }
                if (cmin < 1.0e-3 * p.c_bulk) return "electrolyte_depleted";
                if (thmax > 1.0 - 1.0e-3) return "particles_full";
                if (thmin < 1.0e-3) return "particles_empty";
                return "solver_fail";
            };
            // one constant-voltage time step: find I with V(I) = V_set (see simulate.cv_step in Python)
            auto cv_step = [&](double h, double V_set, double& I) {
                const double tol = 1.0e-9, inf = std::numeric_limits<double>::infinity();
                const std::vector<double> c_start = c, ca_start = ca;
                auto f = [&](double Itry, bool& good) {
                    c = c_start;
                    ca = ca_start;
                    m.i_app = Itry;
                    good = newton_step(h);
                    if (!good) return Itry < 0 ? inf : -inf;
                    return cell_voltage() - V_set;
                };
                bool good;
                double fI = f(I, good);
                if (good && std::abs(fI) <= tol) return true;
                double grow = std::max(std::abs(I), 1.0e-2 * m.i_1C), a = 0, b = 0, fa = 0, fb = 0;
                bool have_a = false, have_b = false;
                for (int it = 0; it < 60; ++it) {
                    if (fI > 0) {
                        a = I; fa = fI; have_a = true;
                        if (have_b) break;
                        I = I < 0 ? 0.0 : I + grow;
                    } else {
                        b = I; fb = fI; have_b = true;
                        if (have_a) break;
                        I = I > 0 ? 0.0 : I - grow;
                    }
                    grow *= 2.0;
                    fI = f(I, good);
                    if (good && std::abs(fI) <= tol) return true;
                }
                if (!(have_a && have_b)) { c = c_start; ca = ca_start; return false; }
                int side = 0;
                for (int it = 0; it < 200; ++it) {
                    if (std::isfinite(fa) && std::isfinite(fb)) {
                        I = (a * fb - b * fa) / (fb - fa);
                        if (!(a < I && I < b)) I = 0.5 * (a + b);
                    } else {
                        I = 0.5 * (a + b);
                    }
                    fI = f(I, good);
                    if (std::abs(fI) <= tol || (b - a) <= 1.0e-14 * m.i_1C) {
                        if (!good) { c = c_start; ca = ca_start; }
                        return good;
                    }
                    if (fI > 0) {
                        a = I; fa = fI;
                        if (side == 1 && std::isfinite(fb)) fb *= 0.5;
                        side = 1;
                    } else {
                        b = I; fb = fI;
                        if (side == -1 && std::isfinite(fa)) fa *= 0.5;
                        side = -1;
                    }
                }
                c = c_start;
                ca = ca_start;
                return false;
            };

            double I = steps[0].kind == Kind::cc ? steps[0].C * m.i_1C : 0.0;
            m.i_app = I;
            write_row(true, 1);
            double last_write = t;
            std::string reason;
            int n_done = 0;
            bool finished = true;
            for (std::size_t k = 0; k < steps.size() && finished; ++k) {
                const Step& st = steps[k];
                const int step_no = static_cast<int>(k) + 1;
                double t_step = 0.0;
                if (st.kind == Kind::cc) I = st.C * m.i_1C;
                else if (st.kind == Kind::rest) I = 0.0;
                while (true) {
                    const double h = st.t < 0 ? dt : std::min(dt, st.t - t_step);
                    double h_done = 0.0;
                    bool stopped = false, ok;
                    std::string why;
                    if (st.kind == Kind::cv) {
                        ok = cv_step(h, st.V, I);
                        h_done = h;
                        stopped = st.Imin >= 0 && std::abs(I) <= st.Imin * m.i_1C;
                        why = "current_limit";
                    } else {
                        m.i_app = I;
                        ok = advance(h, h_done, stopped, st.kind == Kind::cc, st.Vmin, st.Vmax);
                        why = stopped && cell_voltage() <= st.Vmin ? "cutoff_low" : "cutoff_high";
                    }
                    m.i_app = I;
                    if (!ok) { write_row(false, step_no); exit_reason = limit_reason(); finished = false; break; }
                    mAhg = mAhg + 1000.0 * (I / m.mass_area) * h_done / 3600.0;
                    t = t + h_done;
                    t_step = t_step + h_done;
                    ++n_done;
                    bool bad = false;
                    for (double x : c) bad = bad || std::isnan(x);
                    for (double x : ca) bad = bad || std::isnan(x);
                    if (bad) { write_row(false, step_no); exit_reason = "nan"; finished = false; break; }
                    if (stopped || (st.t >= 0 && t_step >= st.t * (1.0 - 1.0e-12))) {
                        write_row(false, step_no);
                        last_write = t;
                        reason = stopped ? why : "duration";
                        break;
                    }
                    if (t - last_write >= p.write_interval) { write_row(false, step_no); last_write = t; }
                    if (t >= 99.0 * 3600.0) { write_row(false, step_no); exit_reason = "max_time"; finished = false; break; }
                }
            }
            if (finished) exit_reason = steps.size() == 1 ? reason : "end_of_protocol";
            nsolve = n_done;
        }
        std::fclose(out);
        std::printf("%s run, C-rate %.17g: exit %s after %d steps; wrote %s\n", p.mode.c_str(), p.C_rate,
                    exit_reason.c_str(), nsolve, p.file.c_str());
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 2;
    }
}
