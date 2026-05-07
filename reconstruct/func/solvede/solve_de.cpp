#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <cctype>

#include "../../include/ffType.hpp"
#include "../../include/parser.hpp"
#include "../../include/diffeq.hpp"
#include "firefly/FFInt.hpp"

namespace {

struct Config {
    unsigned long long p = 0;
};

std::string trim(const std::string& s) {
    const char* ws = " \t\n\r";
    const auto b = s.find_first_not_of(ws);
    if (b == std::string::npos) return "";
    const auto e = s.find_last_not_of(ws);
    return s.substr(b, e - b + 1);
}

Config parseConfig(const std::string& path) {
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
        if (key == "p") cfg.p = std::stoull(val);
    }
    if (cfg.p == 0) throw std::runtime_error("config missing key: p");
    return cfg;
}

template<typename T>
Matrix<T> parseMatrixFile(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open matrix file: " + path);
    std::string content((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    Parser<T> parser;
    return parser.parseMatrix(content);
}

std::vector<firefly::FFInt> parseF0File(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open f0 file: " + path);
    std::string s((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    s = trim(s);
    if (s.empty()) throw std::runtime_error("f0 file is empty");
    if (s.front() == '{' && s.back() == '}') s = s.substr(1, s.size() - 2);

    std::vector<firefly::FFInt> out;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        tok = trim(tok);
        if (tok.empty()) continue;
        long long v = std::stoll(tok);
        out.emplace_back(v);
    }
    return out;
}

void writeSolution(const std::string& outPath, const std::vector<Series<firefly::FFInt>>& sol) {
    std::ofstream out(outPath);
    if (!out) throw std::runtime_error("Cannot open output file: " + outPath);
    out << "{";
    for (size_t i = 0; i < sol.size(); ++i) {
        if (i) out << ",\n ";
        out << sol[i].toString();
    }
    out << "}\n";
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 4 || argc > 5) {
            std::cerr << "Usage: " << argv[0] << " <work_dir> <f0_file> <deg> [out_file]\n";
            return 1;
        }

        const std::string workDir = argv[1];
        const std::string f0Path = argv[2];
        const int deg = std::stoi(argv[3]);
        const std::string outPath = (argc == 5) ? argv[4] : (workDir + "/solution");

        const Config cfg = parseConfig(workDir + "/config");
        firefly::FFInt::set_new_prime(static_cast<uint64_t>(cfg.p));

        auto AX = parseMatrixFile<firefly::FFInt>(workDir + "/AX");
        auto AY = parseMatrixFile<firefly::FFInt>(workDir + "/AY");
        auto f0 = parseF0File(f0Path);

        if (AX.size() != AY.size() || AX.empty()) {
            throw std::runtime_error("AX/AY size mismatch or empty");
        }
        if (f0.size() != AX.size()) {
            std::ostringstream oss;
            oss << "f0 size mismatch: expected " << AX.size() << ", got " << f0.size();
            throw std::runtime_error(oss.str());
        }

        DiffSystem<firefly::FFInt> sys(AX, AY);
        std::vector<Series<firefly::FFInt>> sol;
        sys.solve(sol, f0, deg);
        writeSolution(outPath, sol);

        std::cout << "Solved DE system successfully.\n";
        std::cout << "Output: " << outPath << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}

