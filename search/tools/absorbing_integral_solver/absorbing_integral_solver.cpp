#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "ff_type.hpp"
#include "integral_reduction_searcher.hpp"
#include "io.hpp"
#include "relation_formatter.hpp"
#include "relation_types.hpp"

namespace {

struct BasisTerm {
    int integralId = -1;
    int deltaPower = 0;
    FlintMod coeff;
};

struct BasisRelation {
    std::vector<BasisTerm> terms;
};

std::string trim(const std::string& s) {
    size_t l = 0;
    while (l < s.size() && std::isspace(static_cast<unsigned char>(s[l]))) ++l;
    size_t r = s.size();
    while (r > l && std::isspace(static_cast<unsigned char>(s[r - 1]))) --r;
    return s.substr(l, r - l);
}

int inferNuSize(const std::string& gPath) {
    std::ifstream in(gPath);
    if (!in.is_open()) throw std::runtime_error("Cannot open G file: " + gPath);
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) continue;
        const size_t l = line.find('{');
        const size_t r = line.rfind('}');
        if (l == std::string::npos || r == std::string::npos || r <= l) {
            throw std::runtime_error("Invalid G line: " + line);
        }
        std::stringstream ss(line.substr(l + 1, r - l - 1));
        std::string tok;
        int n = 0;
        while (std::getline(ss, tok, ',')) {
            if (!trim(tok).empty()) ++n;
        }
        return n;
    }
    throw std::runtime_error("G file is empty");
}

std::unordered_map<std::string, std::string> parseHeader(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) throw std::runtime_error("Cannot open basis: " + path);
    std::unordered_map<std::string, std::string> kv;
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

std::vector<BasisRelation> parseRelations(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) throw std::runtime_error("Cannot open basis: " + path);
    std::vector<BasisRelation> relations;
    bool inRelations = false;
    BasisRelation* current = nullptr;
    int termsLeft = 0;
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) continue;
        if (line == "[relations]") {
            inRelations = true;
            continue;
        }
        if (!inRelations) continue;
        if (line.front() == '[') break;

        if (line.rfind("relation ", 0) == 0) {
            std::istringstream iss(line);
            std::string word;
            int relId = -1;
            std::string leadWord, maxWord, termsWord;
            int leadI = -1, leadK = -1, maxK = -1;
            iss >> word >> relId >> leadWord >> leadI >> leadK >> maxWord >> maxK >> termsWord >> termsLeft;
            if (!iss || word != "relation" || leadWord != "lead" || maxWord != "max" || termsWord != "terms") {
                throw std::runtime_error("Bad relation header: " + line);
            }
            relations.emplace_back();
            current = &relations.back();
            current->terms.reserve(static_cast<size_t>(termsLeft));
            continue;
        }

        if (current == nullptr || termsLeft <= 0) {
            throw std::runtime_error("Term outside relation: " + line);
        }
        std::istringstream iss(line);
        BasisTerm term;
        unsigned long long coeff = 0;
        iss >> term.integralId >> term.deltaPower >> coeff;
        if (!iss) throw std::runtime_error("Bad relation term: " + line);
        term.coeff = FlintMod(coeff);
        current->terms.push_back(term);
        --termsLeft;
    }
    if (relations.empty()) {
        throw std::runtime_error("No relations found in absorbing basis");
    }
    return relations;
}

FlintMod powDelta(FlintMod x, int k) {
    FlintMod out(1);
    for (int i = 0; i < k; ++i) out = out * x;
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <G_path> <absorbing_relation_path> <delta_value> <output_path>\n";
        return 1;
    }
    try {
        const std::string gPath = argv[1];
        const std::string basisPath = argv[2];

        auto kv = parseHeader(basisPath);
        mp_limb_t p = 2305843009213693951ULL;
        if (auto it = kv.find("p"); it != kv.end()) {
            p = static_cast<mp_limb_t>(std::stoull(it->second));
        }
        FlintMod::set_modulus(p);
        const FlintMod delta(static_cast<unsigned long long>(std::stoull(argv[3])));

        const int nuSize = inferNuSize(gPath);
        const auto G = parseSearchTargetFile(gPath, nuSize);
        const auto basis = parseRelations(basisPath);

        std::vector<FlintMod> deltaPowers(1, FlintMod(1));
        int maxK = 0;
        for (const auto& rel : basis) {
            for (const auto& term : rel.terms) maxK = std::max(maxK, term.deltaPower);
        }
        for (int k = 1; k <= maxK; ++k) deltaPowers.push_back(deltaPowers.back() * delta);

        std::vector<IntegralRelation<FlintMod>> integralRelations;
        integralRelations.reserve(basis.size());
        for (const auto& rel : basis) {
            std::vector<FlintMod> coeffs(G.size(), FlintMod(0));
            for (const auto& term : rel.terms) {
                if (term.integralId < 0 || term.integralId >= static_cast<int>(G.size())) {
                    throw std::runtime_error("Relation integral id out of range");
                }
                coeffs[static_cast<size_t>(term.integralId)] +=
                    term.coeff * deltaPowers[static_cast<size_t>(term.deltaPower)];
            }
            bool nonzero = false;
            for (const auto& c : coeffs) {
                if (c != FlintMod(0)) {
                    nonzero = true;
                    break;
                }
            }
            if (!nonzero) continue;
            IntegralRelation<FlintMod> outRel;
            outRel.integrals = G;
            outRel.coeffs = std::move(coeffs);
            integralRelations.push_back(std::move(outRel));
        }

        std::ofstream out(argv[4]);
        if (!out.is_open()) throw std::runtime_error("Cannot open output: " + std::string(argv[4]));
        out << "# p = " << p << "\n";
        out << "# delta = " << delta << "\n";
        out << "# |G| = " << G.size() << "\n";
        out << "# absorbing relations = " << basis.size() << "\n\n";

        if (integralRelations.empty()) {
            out << "# integral variables = 0\n";
            out << "# integral pivot columns = 0\n";
            out << "# integral free columns = 0\n";
            out << "#MIs\n\n[reductions]\n\n[relations]\n\n[integral_rref]\n";
            return 0;
        }

        IntegralReductionSearcher<FlintMod> searcher(integralRelations);
        const auto result = searcher.search();
        RelationFormatter<FlintMod>::writeIntegralReductionSummary(out, result);
        RelationFormatter<FlintMod>::writeIntegralMasterBasis(out, result);
        out << "\n";
        RelationFormatter<FlintMod>::writeIntegralReductions(out, result);
        out << "\n";
        RelationFormatter<FlintMod>::writeIntegralRelations(out, integralRelations);
        out << "\n";
        RelationFormatter<FlintMod>::writeIntegralRREF(out, result);

        std::cerr << "[absorbing_integral_solver] G size = " << G.size()
                  << ", relations = " << integralRelations.size()
                  << ", free columns = " << result.freeColumns.size() << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
