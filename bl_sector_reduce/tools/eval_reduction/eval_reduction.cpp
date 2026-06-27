#include <cctype>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "ff_type.hpp"
#include "label.hpp"
#include "polynomial_1d.hpp"

namespace {

std::string trim(const std::string& s) {
    size_t l = 0;
    while (l < s.size() && std::isspace(static_cast<unsigned char>(s[l]))) ++l;
    size_t r = s.size();
    while (r > l && std::isspace(static_cast<unsigned char>(s[r - 1]))) --r;
    return s.substr(l, r - l);
}

std::vector<FlintMod> parsePoly(const std::string& raw) {
    const size_t l = raw.find('{');
    const size_t r = raw.rfind('}');
    if (l == std::string::npos || r == std::string::npos || r <= l) {
        throw std::runtime_error("Invalid polynomial: " + raw);
    }
    std::vector<FlintMod> coeffs;
    std::stringstream ss(raw.substr(l + 1, r - l - 1));
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        tok = trim(tok);
        if (!tok.empty()) coeffs.emplace_back(static_cast<unsigned long long>(std::stoull(tok)));
    }
    return coeffs;
}

std::map<std::string, std::string> parseHeader(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) throw std::runtime_error("Cannot open reduction file: " + path);
    std::map<std::string, std::string> kv;
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) continue;
        if (line[0] != '#') break;
        const size_t eq = line.find('=');
        if (eq == std::string::npos) continue;
        kv[trim(line.substr(1, eq - 1))] = trim(line.substr(eq + 1));
    }
    return kv;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <reduction_path> <object_label> <delta_value> <nu_size>\n";
        return 1;
    }
    const std::string reductionPath = argv[1];
    const std::string objectText = argv[2];
    const std::string deltaText = argv[3];
    const int nuSize = std::stoi(argv[4]);

    try {
        const auto header = parseHeader(reductionPath);
        auto pIt = header.find("p");
        if (pIt == header.end()) throw std::runtime_error("Missing # p header");
        FlintMod::set_modulus(static_cast<mp_limb_t>(std::stoull(pIt->second)));
        const FlintMod delta(static_cast<unsigned long long>(std::stoull(deltaText)));
        const std::string object = objectLabelToString(parseObjectLabel(objectText, nuSize));

        std::ifstream in(reductionPath);
        if (!in.is_open()) throw std::runtime_error("Cannot open reduction file: " + reductionPath);

        bool inBlock = false;
        bool active = false;
        bool zero = false;
        Polynomial1D<FlintMod> den;
        std::vector<std::pair<std::string, Polynomial1D<FlintMod>>> terms;
        std::string line;

        auto flush = [&]() {
            if (!active) return false;
            if (zero) {
                std::cout << object << " = 0\n";
                return true;
            }
            const FlintMod denVal = den.eval(delta);
            std::cout << object << " = ";
            bool first = true;
            for (const auto& term : terms) {
                const FlintMod coeff = term.second.eval(delta) / denVal;
                if (!first) std::cout << " + ";
                std::cout << coeff << "*" << term.first;
                first = false;
            }
            std::cout << "\n";
            return true;
        };

        while (std::getline(in, line)) {
            line = trim(line);
            if (line == "[sector_reductions]") {
                inBlock = true;
                continue;
            }
            if (line == "[global_reductions]") break;
            if (!inBlock) continue;
            if (line.empty()) {
                if (flush()) return 0;
                active = false;
                zero = false;
                terms.clear();
                den = Polynomial1D<FlintMod>();
                continue;
            }
            if (line.rfind("object=", 0) == 0) {
                const std::string val = objectLabelToString(parseObjectLabel(line.substr(7), nuSize));
                active = val == object;
                zero = false;
                terms.clear();
                den = Polynomial1D<FlintMod>();
                continue;
            }
            if (!active) continue;
            if (line == "zero") {
                zero = true;
            } else if (line.rfind("den=", 0) == 0) {
                den = Polynomial1D<FlintMod>(parsePoly(line.substr(4)));
            } else if (line.rfind("term ", 0) == 0) {
                const size_t eq = line.find('=');
                if (eq == std::string::npos) throw std::runtime_error("Invalid term line: " + line);
                const std::string label =
                    objectLabelToString(parseObjectLabel(line.substr(5, eq - 5), nuSize));
                terms.emplace_back(label, Polynomial1D<FlintMod>(parsePoly(line.substr(eq + 1))));
            }
        }
        if (flush()) return 0;
        throw std::runtime_error("Object not found in reductions: " + object);
    } catch (const std::exception& e) {
        std::cerr << "[eval_reduction] Error: " << e.what() << "\n";
        return 1;
    }
}
