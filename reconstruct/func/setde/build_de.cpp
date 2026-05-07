#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <filesystem>
#include <stdexcept>
#include <cstdlib>

#include "../../include/ffType.hpp"
#include "../../include/parser.hpp"
#include "../../include/de_interpolater.hpp"
#include "firefly/FFInt.hpp"

namespace {

struct Config {
    int B = -1;
    int N = -1;
    unsigned long long p = 0;
    unsigned long long d = 0; // delta
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
    if (!in) {
        throw std::runtime_error("Cannot open config file: " + path);
    }

    Config cfg;
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) continue;

        auto pos = line.find('=');
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
        throw std::runtime_error("config missing required keys: B, N, p, d");
    }
    return cfg;
}

template<typename T>
Matrix<T> parseSMatrix(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Cannot open S file: " + path);
    }
    std::string content((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    Parser<T> parser;
    return parser.parseMatrix(content);
}

template<typename T>
void writeRationalMatrix(const std::string& path, const Matrix<T>& mat) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Cannot open output file: " + path);
    }

    out << "{";
    for (size_t i = 0; i < mat.size(); ++i) {
        if (i) out << ",\n ";
        out << "{";
        for (size_t j = 0; j < mat[i].size(); ++j) {
            if (j) out << ", ";
            out << mat[i][j].toString();
        }
        out << "}";
    }
    out << "}\n";
}

template<typename T>
std::vector<std::vector<Polynomial<T>>> toPolynomialTopS(const Matrix<T>& ratS) {
    std::vector<std::vector<Polynomial<T>>> polyS(
        ratS.size(), std::vector<Polynomial<T>>(ratS[0].size()));

    for (size_t i = 0; i < ratS.size(); ++i) {
        for (size_t j = 0; j < ratS[i].size(); ++j) {
            const auto& r = ratS[i][j];
            if (!r.denominator.isConstant()) {
                throw std::runtime_error("S contains non-constant rational denominator at entry (" +
                                         std::to_string(i) + "," + std::to_string(j) + ")");
            }
            T den0 = r.denominator.getCoeff(0, 0);
            if (den0 == T(0)) {
                throw std::runtime_error("S contains zero denominator at entry (" +
                                         std::to_string(i) + "," + std::to_string(j) + ")");
            }
            Polynomial<T> p = r.numerator;
            p *= (T(1) / den0);
            polyS[i][j] = p;
        }
    }
    return polyS;
}

template<typename T>
Polynomial<T> buildShiftedUPolynomial(const T& a, const T& b) {
    // Urect(X,Y)=X-X^2+Y-2XY+X^2Y-Y^2+2XY^2-X^2Y^2
    // Use shifted form Urect(X+a,Y+b) with explicit coefficients.
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

    U.addMonomial(c00, Power(0, 0));
    U.addMonomial(c10, Power(1, 0));
    U.addMonomial(c20, Power(2, 0));
    U.addMonomial(c01, Power(0, 1));
    U.addMonomial(c11, Power(1, 1));
    U.addMonomial(c21, Power(2, 1));
    U.addMonomial(c02, Power(0, 2));
    U.addMonomial(c12, Power(1, 2));
    U.addMonomial(c22, Power(2, 2));
    return U;
}

} // namespace

int main(int argc, char** argv) {
    try {
        const std::string workDir = (argc >= 2) ? argv[1] : "./db";
        const std::string configPath = workDir + "/config";
        const std::string sPathNormalized = workDir + "/S_normalized";
        const std::string sPathRaw = workDir + "/S";
        const std::string sPath = std::filesystem::exists(sPathNormalized) ? sPathNormalized : sPathRaw;

        Config cfg = parseConfigFile(configPath);
        firefly::FFInt::set_new_prime(static_cast<uint64_t>(cfg.p));

        Matrix<firefly::FFInt> topSRat = parseSMatrix<firefly::FFInt>(sPath);
        auto topSPoly = toPolynomialTopS(topSRat);

        const int expected = cfg.B + cfg.N;
        if (static_cast<int>(topSRat.size()) != expected || static_cast<int>(topSRat[0].size()) != expected) {
            std::ostringstream oss;
            oss << "S dimension mismatch: expected " << expected << "x" << expected
                << ", got " << topSRat.size() << "x" << topSRat[0].size();
            throw std::runtime_error(oss.str());
        }

        Polynomial<firefly::FFInt> UPoly = buildShiftedUPolynomial(
            firefly::FFInt(cfg.a), firefly::FFInt(cfg.b));
        DEInterpolater<firefly::FFInt> interpolater(
            topSPoly, UPoly, cfg.N, cfg.B, firefly::FFInt(cfg.d), static_cast<uint64_t>(cfg.p));

        // Optional one-point numeric debug mode:
        // export SETDE_POINT_X=... SETDE_POINT_Y=... and run build_de.
        if (const char* xEnv = std::getenv("SETDE_POINT_X")) {
            if (const char* yEnv = std::getenv("SETDE_POINT_Y")) {
                firefly::FFInt X0(std::stoull(xEnv));
                firefly::FFInt Y0(std::stoull(yEnv));
                auto dUxPoly = UPoly.derivativeX();
                auto dUyPoly = UPoly.derivativeY();
                auto U0 = UPoly.evaluate(X0, Y0);
                auto dUx0 = dUxPoly.evaluate(X0, Y0);
                auto dUy0 = dUyPoly.evaluate(X0, Y0);
                auto dbgBuilder = interpolater.createDEBuilder(X0, Y0);
                const auto& masterNus = dbgBuilder.getReducer().getMasterNus();
                std::cout << "DEBUG_MASTER_TOTALS=";
                for (size_t i = 0; i < masterNus.size(); ++i) {
                    int s = 0;
                    for (int v : masterNus[i]) s += v;
                    if (i) std::cout << ",";
                    std::cout << s;
                }
                std::cout << "\n";
                for (size_t i = 0; i < masterNus.size(); ++i) {
                    std::cout << "DEBUG_MASTER_NU_" << (i + 1) << "=";
                    for (size_t j = 0; j < masterNus[i].size(); ++j) {
                        if (j) std::cout << ",";
                        std::cout << masterNus[i][j];
                    }
                    std::cout << "\n";
                }
                std::vector<std::vector<firefly::FFInt>> dbgAX, dbgAY;
                dbgBuilder.buildDEMatrices(dbgAX, dbgAY);
                std::cout << "DEBUG_POINT X=" << xEnv << " Y=" << yEnv << "\n";
                std::cout << "DEBUG_U=" << U0.n << " DEBUG_dUdX=" << dUx0.n
                          << " DEBUG_dUdY=" << dUy0.n << "\n";
                if (!dbgAX.empty() && !dbgAX[0].empty() && !dbgAY.empty() && !dbgAY[0].empty()) {
                    std::cout << "DEBUG_NUM_AX11=" << dbgAX[0][0].n << "\n";
                    std::cout << "DEBUG_NUM_AY11=" << dbgAY[0][0].n << "\n";
                }
            }
        }

        Matrix<firefly::FFInt> AX, AY;
        interpolater.buildDEMatrices(AX, AY);

        writeRationalMatrix(workDir + "/AX", AX);
        writeRationalMatrix(workDir + "/AY", AY);

        std::cout << "Built DE matrices successfully.\n";
        std::cout << "Output: " << workDir << "/AX and " << workDir << "/AY\n";
        std::cout << "Firefly probes: " << interpolater.getTotalProbes() << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
