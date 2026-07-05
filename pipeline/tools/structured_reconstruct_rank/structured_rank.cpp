#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

static constexpr uint64_t P = 2305843009213693951ULL;
using Poly = std::vector<uint64_t>;

static uint64_t addm(uint64_t a, uint64_t b) {
    uint64_t c = a + b;
    return c >= P ? c - P : c;
}

static uint64_t subm(uint64_t a, uint64_t b) {
    return a >= b ? a - b : P - (b - a);
}

static uint64_t mulm(uint64_t a, uint64_t b) {
    return static_cast<uint64_t>((static_cast<__uint128_t>(a) * b) % P);
}

static uint64_t powm(uint64_t a, uint64_t e) {
    uint64_t r = 1;
    while (e) {
        if (e & 1) r = mulm(r, a);
        e >>= 1;
        if (e) a = mulm(a, a);
    }
    return r;
}

static Poly fit(Poly v, int L) {
    v.resize(L, 0);
    return v;
}

static Poly addp(const Poly& a, const Poly& b, int L) {
    Poly o(L);
    for (int i = 0; i < L; ++i) o[i] = addm(a[i], b[i]);
    return o;
}

static Poly subp(const Poly& a, const Poly& b, int L) {
    Poly o(L);
    for (int i = 0; i < L; ++i) o[i] = subm(a[i], b[i]);
    return o;
}

static Poly scalep(const Poly& a, uint64_t c, int L) {
    Poly o(L);
    for (int i = 0; i < L; ++i) o[i] = mulm(a[i], c);
    return o;
}

static Poly shiftp(const Poly& a, int k, int L) {
    Poly o(L);
    if (k < L) {
        for (int i = 0; i + k < L; ++i) o[i + k] = a[i];
    }
    return o;
}

static Poly mulp(const Poly& a, const Poly& b, int L) {
    Poly o(L);
    for (int i = 0; i < L; ++i) {
        if (!a[i]) continue;
        for (int j = 0; i + j < L; ++j) {
            if (b[j]) o[i + j] = addm(o[i + j], mulm(a[i], b[j]));
        }
    }
    return o;
}

static Poly powp(Poly a, int e, int L) {
    Poly o(L);
    o[0] = 1;
    while (e) {
        if (e & 1) o = mulp(o, a, L);
        e >>= 1;
        if (e) a = mulp(a, a, L);
    }
    return o;
}

static Poly invs(const Poly& a, int L) {
    Poly o(L);
    uint64_t inv0 = powm(a[0], P - 2);
    o[0] = inv0;
    for (int k = 1; k < L; ++k) {
        uint64_t s = 0;
        for (int i = 1; i <= k && i < static_cast<int>(a.size()); ++i) {
            s = addm(s, mulm(a[i], o[k - i]));
        }
        o[k] = mulm(subm(0, s), inv0);
    }
    return o;
}

static Poly divs(const Poly& a, const Poly& b, int L) {
    return mulp(a, invs(b, L), L);
}

static std::vector<std::string> split(const std::string& line) {
    std::istringstream is(line);
    std::vector<std::string> out;
    std::string x;
    while (is >> x) out.push_back(x);
    return out;
}

static uint64_t to_mod(const std::string& s) {
    return std::stoull(s) % P;
}

struct Term {
    std::string master;
    uint64_t c = 0;
    std::vector<std::pair<std::string, int>> n;
    std::vector<std::pair<std::string, int>> d;
};

struct Input {
    Poly p0;
    std::vector<std::pair<std::string, Poly>> red;
    std::unordered_map<std::string, Poly> ms;
    std::map<std::string, int> degree;
    std::unordered_map<std::string, Poly> factor;
    std::vector<Term> terms;
};

static Input read_flat(const std::string& path, int L) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    Input x;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        auto t = split(line);
        if (t.empty()) continue;
        if (t[0] == "P") {
            uint64_t p = std::stoull(t[1]);
            if (p != P) throw std::runtime_error("unexpected prime");
        } else if (t[0] == "P0") {
            for (size_t i = 1; i < t.size(); ++i) x.p0.push_back(to_mod(t[i]));
            x.p0 = fit(std::move(x.p0), L);
        } else if (t[0] == "RED") {
            Poly a;
            for (size_t i = 2; i < t.size(); ++i) a.push_back(to_mod(t[i]));
            x.red.push_back({t[1], fit(std::move(a), L)});
        } else if (t[0] == "MS") {
            Poly a;
            for (size_t i = 2; i < t.size(); ++i) a.push_back(to_mod(t[i]));
            x.ms[t[1]] = fit(std::move(a), L);
        } else if (t[0] == "FAC") {
            std::string id = t[1];
            int deg = std::stoi(t[3]);
            Poly a;
            for (size_t i = 4; i < t.size(); ++i) a.push_back(to_mod(t[i]));
            x.degree[id] = deg;
            x.factor[id] = fit(std::move(a), L);
        } else if (t[0] == "TERM") {
            Term term;
            term.master = t[1];
            term.c = to_mod(t[2]);
            int nc = std::stoi(t[3]);
            int dc = std::stoi(t[4]);
            for (int i = 0; i < nc; ++i) {
                if (!std::getline(in, line)) throw std::runtime_error("short N list");
                auto u = split(line);
                if (u[0] != "N") throw std::runtime_error("expected N line");
                term.n.push_back({u[1], std::stoi(u[2])});
            }
            for (int i = 0; i < dc; ++i) {
                if (!std::getline(in, line)) throw std::runtime_error("short D list");
                auto u = split(line);
                if (u[0] != "D") throw std::runtime_error("expected D line");
                term.d.push_back({u[1], std::stoi(u[2])});
            }
            x.terms.push_back(std::move(term));
        }
    }
    if (x.p0.empty()) throw std::runtime_error("missing P0");
    return x;
}

static int rank_mod(std::vector<std::vector<uint64_t>>& a) {
    int rows = static_cast<int>(a.size());
    if (!rows) return 0;
    int cols = static_cast<int>(a[0].size());
    int r = 0;
    for (int c = 0; c < cols && r < rows; ++c) {
        int piv = -1;
        for (int i = r; i < rows; ++i) {
            if (a[i][c]) {
                piv = i;
                break;
            }
        }
        if (piv < 0) continue;
        std::swap(a[r], a[piv]);
        uint64_t inv = powm(a[r][c], P - 2);
        for (int j = c; j < cols; ++j) a[r][j] = mulm(a[r][j], inv);
        for (int i = 0; i < rows; ++i) {
            if (i == r || !a[i][c]) continue;
            uint64_t f = a[i][c];
            for (int j = c; j < cols; ++j) a[i][j] = subm(a[i][j], mulm(f, a[r][j]));
        }
        ++r;
    }
    return r;
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: structured_rank <structured_input.flat> <N>\n";
        return 2;
    }
    const std::string path = argv[1];
    const int N = std::stoi(argv[2]);
    const int L = N + 1;

    auto t0 = std::chrono::steady_clock::now();
    Input in = read_flat(path, L);

    Poly C(L, 0);
    for (const auto& rt : in.red) {
        C = addp(C, mulp(rt.second, in.ms.at(rt.first), L), L);
    }
    C = divs(C, in.p0, L);

    std::set<std::string> used;
    std::map<std::string, int> max_beta;
    for (const auto& term : in.terms) {
        for (const auto& p : term.n) used.insert(p.first);
        for (const auto& p : term.d) {
            used.insert(p.first);
            max_beta[p.first] = std::max(max_beta[p.first], p.second);
        }
    }

    std::unordered_map<std::string, Poly> finv;
    for (const auto& id : used) finv[id] = invs(in.factor.at(id), L);

    Poly D(L, 0);
    D[0] = 1;
    for (const auto& kv : max_beta) {
        D = mulp(D, powp(in.factor.at(kv.first), kv.second, L), L);
    }
    Poly lhs = mulp(C, D, L);

    std::vector<Poly> term_series;
    term_series.reserve(in.terms.size());
    for (const auto& term : in.terms) {
        Poly prod(L, 0);
        prod[0] = term.c;
        prod = mulp(prod, in.ms.at(term.master), L);
        for (const auto& p : term.n) prod = mulp(prod, powp(in.factor.at(p.first), p.second, L), L);

        std::map<std::string, int> beta;
        for (const auto& p : term.d) beta[p.first] += p.second;
        for (const auto& kv : max_beta) {
            int q = kv.second - beta[kv.first];
            if (q) prod = mulp(prod, powp(in.factor.at(kv.first), q, L), L);
        }
        term_series.push_back(std::move(prod));
    }

    Poly eq = lhs;
    for (const auto& s : term_series) eq = subp(eq, s, L);
    int nonzero = 0;
    for (uint64_t v : eq) {
        if (v) ++nonzero;
    }

    std::vector<std::pair<std::string, int>> vars;
    for (const auto& id : used) {
        for (int k = 1; k <= in.degree.at(id) && k < L; ++k) vars.push_back({id, k});
    }
    for (const auto& term : in.terms) vars.push_back({"c:" + term.master, -1});

    std::unordered_map<std::string, Poly> lhs_over;
    for (const auto& kv : max_beta) lhs_over[kv.first] = mulp(lhs, finv.at(kv.first), L);

    std::vector<std::unordered_map<std::string, Poly>> term_over(term_series.size());
    for (size_t j = 0; j < term_series.size(); ++j) {
        for (const auto& id : used) term_over[j][id] = mulp(term_series[j], finv.at(id), L);
    }

    std::vector<std::vector<uint64_t>> mat(L, std::vector<uint64_t>(vars.size()));
    for (size_t v = 0; v < vars.size(); ++v) {
        Poly col(L, 0);
        const auto& name = vars[v].first;
        int k = vars[v].second;
        if (name.rfind("c:", 0) == 0) {
            std::string master = name.substr(2);
            size_t j = 0;
            while (j < in.terms.size() && in.terms[j].master != master) ++j;
            if (j == in.terms.size()) throw std::runtime_error("missing constant term " + master);
            uint64_t neg_inv_c = subm(0, powm(in.terms[j].c, P - 2));
            col = scalep(term_series[j], neg_inv_c, L);
        } else {
            auto mb = max_beta.find(name);
            int maxb_val = mb == max_beta.end() ? 0 : mb->second;
            if (mb != max_beta.end()) {
                col = addp(col, scalep(shiftp(lhs_over.at(name), k, L), maxb_val, L), L);
            }
            for (size_t j = 0; j < in.terms.size(); ++j) {
                int alpha = 0;
                int beta = 0;
                for (const auto& p : in.terms[j].n) {
                    if (p.first == name) alpha += p.second;
                }
                for (const auto& p : in.terms[j].d) {
                    if (p.first == name) beta += p.second;
                }
                int q = maxb_val - beta;
                int coeff = alpha + q;
                if (coeff) {
                    col = subp(col, scalep(shiftp(term_over[j].at(name), k, L), static_cast<uint64_t>(coeff), L), L);
                }
            }
        }
        for (int i = 0; i < L; ++i) mat[i][v] = col[i];
    }

    int r = rank_mod(mat);
    auto t1 = std::chrono::steady_clock::now();
    double sec = std::chrono::duration<double>(t1 - t0).count();

    std::cout << "input " << path << "\n";
    std::cout << "N " << N << " equations " << L << " unknowns " << vars.size()
              << " rank " << r << " full " << (r == static_cast<int>(vars.size()) ? "true" : "false") << "\n";
    std::cout << "known_solution_equation_nonzero " << nonzero << "\n";
    std::cout << "elapsed_sec " << sec << "\n";
    return r == static_cast<int>(vars.size()) && nonzero == 0 ? 0 : 1;
}
