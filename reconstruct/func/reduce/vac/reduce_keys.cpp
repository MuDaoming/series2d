#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>

#include "../../../include/ffType.hpp"
#include "../../../include/parser.hpp"
#include "../../../include/fbi_reducer.hpp"
#include "firefly/FFInt.hpp"

namespace {

struct Config {
    int B = -1;
    int N = -1;
    unsigned long long p = 0;
    unsigned long long d = 0;
    unsigned long long a = 0;
    unsigned long long b = 0;
};

std::string trim(const std::string& s) {
    const char* ws = " \t\n\r";
    const auto b = s.find_first_not_of(ws);
    if (b == std::string::npos) return "";
    const auto e = s.find_last_not_of(ws);
    return s.substr(b, e - b + 1);
}

Config parseConfigFile(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open config file: " + path);
    Config cfg;
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) continue;
        const auto pos = line.find('=');
        if (pos == std::string::npos) continue;
        const std::string key = trim(line.substr(0, pos));
        const std::string val = trim(line.substr(pos + 1));
        if (key == "B") cfg.B = std::stoi(val);
        else if (key == "N") cfg.N = std::stoi(val);
        else if (key == "p") cfg.p = std::stoull(val);
        else if (key == "d") cfg.d = std::stoull(val);
        else if (key == "a") cfg.a = std::stoull(val);
        else if (key == "b") cfg.b = std::stoull(val);
    }
    if (cfg.B <= 0 || cfg.N <= 0 || cfg.p == 0 || cfg.d == 0) {
        throw std::runtime_error("config missing required keys: B,N,p,d");
    }
    return cfg;
}

template<typename T>
Matrix<T> parseSMatrix(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open S file: " + path);
    std::string content((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    Parser<T> parser;
    return parser.parseMatrix(content);
}

template<typename T>
std::vector<std::vector<Polynomial<T>>> toPolynomialTopS(const Matrix<T>& ratS) {
    std::vector<std::vector<Polynomial<T>>> polyS(ratS.size(), std::vector<Polynomial<T>>(ratS[0].size()));
    for (size_t i = 0; i < ratS.size(); ++i) {
        for (size_t j = 0; j < ratS[i].size(); ++j) {
            const auto& r = ratS[i][j];
            if (!r.denominator.isConstant()) {
                throw std::runtime_error("S entry has non-constant denominator");
            }
            T den0 = r.denominator.getCoeff(0, 0);
            if (den0 == T(0)) throw std::runtime_error("S entry has zero denominator");
            Polynomial<T> p = r.numerator;
            p *= (T(1) / den0);
            polyS[i][j] = p;
        }
    }
    return polyS;
}

template<typename T>
Polynomial<T> buildShiftedUPolynomial(const T& a, const T& b) {
    Polynomial<T> U;
    const T one = T(1);
    const T two = T(2);
    const T c00 = a - a * a + b - two * a * b + a * a * b - b * b + two * a * b * b - a * a * b * b;
    const T c10 = one - two * a - two * b + two * a * b + two * b * b - two * a * b * b;
    const T c20 = T(-1) + b - b * b;
    const T c01 = one - two * a + a * a - two * b + T(4) * a * b - two * a * a * b;
    const T c11 = T(-2) + two * a + T(4) * b - T(4) * a * b;
    const T c21 = one - two * b;
    const T c02 = T(-1) + two * a - a * a;
    const T c12 = two - two * a;
    const T c22 = T(-1);
    U.addMonomial(c00, Power(0, 0)); U.addMonomial(c10, Power(1, 0)); U.addMonomial(c20, Power(2, 0));
    U.addMonomial(c01, Power(0, 1)); U.addMonomial(c11, Power(1, 1)); U.addMonomial(c21, Power(2, 1));
    U.addMonomial(c02, Power(0, 2)); U.addMonomial(c12, Power(1, 2)); U.addMonomial(c22, Power(2, 2));
    return U;
}

struct Key { std::vector<int> nu; firefly::FFInt delta; };

std::vector<Key> parseKeys(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open key file: " + path);
    std::vector<Key> out;
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) continue;
        auto semi = line.find(';');
        if (semi == std::string::npos) continue;
        std::string left = trim(line.substr(0, semi));
        std::string right = trim(line.substr(semi + 1));
        if (left.size() < 2 || left.front() != '{' || left.back() != '}') {
            throw std::runtime_error("Invalid key line: " + line);
        }
        left = left.substr(1, left.size() - 2);
        std::vector<int> nu;
        std::stringstream ss(left);
        std::string tok;
        while (std::getline(ss, tok, ',')) {
            tok = trim(tok);
            if (!tok.empty()) nu.push_back(std::stoi(tok));
        }
        out.push_back({nu, firefly::FFInt(std::stoll(right))});
    }
    return out;
}

int toSigned(const firefly::FFInt& x) {
    long long v = static_cast<long long>(x.n);
    long long p = static_cast<long long>(firefly::FFInt::p);
    if (v > p / 2) v -= p;
    return static_cast<int>(v);
}

int deltaPowU(int L, const firefly::FFInt& D_in,
              int nuTotT, const firefly::FFInt& DT,
              int nuTotS, const firefly::FFInt& DS) {
    int offsetT = toSigned(DT - D_in);
    int offsetS = toSigned(DS - D_in);
    bool parityT = (offsetT % 2 == 0);
    bool parityS = (offsetS % 2 == 0);
    int DbarT = offsetT + (parityT ? 0 : 1);
    int DbarS = offsetS + (parityS ? 0 : 1);
    int DbarDiff = DbarT - DbarS;
    return (nuTotT - nuTotS) - (L + 1) * DbarDiff / 2;
}

firefly::FFInt powInt(firefly::FFInt base, int exp) {
    if (exp == 0) return firefly::FFInt(1);
    if (exp < 0) return firefly::FFInt(1) / powInt(base, -exp);
    firefly::FFInt res(1);
    int e = exp;
    while (e > 0) {
        if (e & 1) res = res * base;
        e >>= 1;
        if (e) base = base * base;
    }
    return res;
}

int nuTot(const std::vector<int>& nu) { int s = 0; for (int v : nu) s += v; return s; }

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 5) {
            std::cerr << "Usage: " << argv[0] << " <S_path> <config_path> <keys_path> <output_path>\n";
            return 1;
        }
        const std::string sPath = argv[1];
        const std::string configPath = argv[2];
        const std::string keysPath = argv[3];
        const std::string outPath = argv[4];

        const Config cfg = parseConfigFile(configPath);
        firefly::FFInt::set_new_prime(static_cast<uint64_t>(cfg.p));

        Matrix<firefly::FFInt> topSRat = parseSMatrix<firefly::FFInt>(sPath);
        auto topSPoly = toPolynomialTopS(topSRat);
        std::vector<std::vector<firefly::FFInt>> numTopS(topSPoly.size(), std::vector<firefly::FFInt>(topSPoly[0].size()));
        const firefly::FFInt evalX(0);
        const firefly::FFInt evalY(0);
        for (size_t i = 0; i < topSPoly.size(); ++i) {
            for (size_t j = 0; j < topSPoly[i].size(); ++j) {
                numTopS[i][j] = topSPoly[i][j].evaluate(evalX, evalY);
            }
        }

        FBIReducer<firefly::FFInt> reducer(numTopS, cfg.N, cfg.B, firefly::FFInt(cfg.d));
        auto keys = parseKeys(keysPath);

        auto shiftedU = buildShiftedUPolynomial(firefly::FFInt(cfg.a), firefly::FFInt(cfg.b));
        firefly::FFInt U0 = shiftedU.evaluate(firefly::FFInt(0), firefly::FFInt(0));

        const auto& masterNus = reducer.getMasterNus();
        const auto& masterDeltas = reducer.getMasterDeltas();

        std::ofstream out(outPath);
        if (!out) throw std::runtime_error("Cannot open output file: " + outPath);

        out << "p=" << cfg.p << " a=" << cfg.a << " b=" << cfg.b << " d=" << cfg.d << "\n";
        out << "S_eval_point={0,0}\n";
        out << "U(a,b)=" << U0.n << "\n";
        out << "num_master=" << masterNus.size() << "\n\n";

        for (size_t idx = 0; idx < keys.size(); ++idx) {
            const auto& nu = keys[idx].nu;
            const auto delta = keys[idx].delta;
            const auto coeffRaw = reducer.getReductionCoeff(nu, delta);

            out << "KEY " << idx + 1 << ": nu={";
            for (size_t i = 0; i < nu.size(); ++i) { if (i) out << ","; out << nu[i]; }
            out << "}; delta=" << delta.n << "\n";

            out << "RAW={";
            for (size_t j = 0; j < coeffRaw.size(); ++j) {
                if (j) out << ",";
                out << coeffRaw[j].n;
            }
            out << "}\n";

            const int nuT = nuTot(nu);
            out << "REDEF={";
            for (size_t j = 0; j < coeffRaw.size(); ++j) {
                if (j) out << ",";
                int dp = deltaPowU(2, firefly::FFInt(cfg.d), nuT, delta, nuTot(masterNus[j]), masterDeltas[j]);
                firefly::FFInt c = coeffRaw[j] * powInt(U0, dp);
                out << c.n;
            }
            out << "}\n\n";
        }

        std::cout << "Done: " << outPath << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
