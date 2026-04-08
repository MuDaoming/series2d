#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "coefficient_relation_expander.hpp"
#include "ff_type.hpp"
#include "fi_reduction_searcher.hpp"
#include "io.hpp"
#include "relation_formatter.hpp"
#include "relation_types.hpp"

namespace {

std::string trim(const std::string& s) {
    size_t l = 0;
    while (l < s.size() && std::isspace(static_cast<unsigned char>(s[l]))) ++l;
    size_t r = s.size();
    while (r > l && std::isspace(static_cast<unsigned char>(s[r - 1]))) --r;
    return s.substr(l, r - l);
}

int inferNuSize(const std::string& gPath) {
    std::ifstream in(gPath);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open G file: " + gPath);
    }

    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) continue;

        const size_t l = line.find('{');
        const size_t r = line.rfind('}');
        if (l == std::string::npos || r == std::string::npos || r <= l) {
            throw std::runtime_error("Invalid nu line in G file: " + line);
        }

        std::string body = line.substr(l + 1, r - l - 1);
        std::stringstream ss(body);
        std::string tok;
        int cnt = 0;
        while (std::getline(ss, tok, ',')) {
            tok = trim(tok);
            if (!tok.empty()) ++cnt;
        }
        if (cnt <= 0) {
            throw std::runtime_error("Empty nu entry in G file");
        }
        return cnt;
    }

    throw std::runtime_error("G file is empty: " + gPath);
}

std::unordered_map<std::string, std::string> parseHeaderKV(const std::string& polyPath) {
    std::ifstream in(polyPath);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open poly_relation file: " + polyPath);
    }

    std::unordered_map<std::string, std::string> kv;
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) continue;
        if (line[0] != '#') {
            break;
        }
        const size_t eq = line.find('=');
        if (eq == std::string::npos) continue;

        std::string key = trim(line.substr(1, eq - 1));
        std::string val = trim(line.substr(eq + 1));
        if (!key.empty() && !val.empty()) {
            kv[key] = val;
        }
    }
    return kv;
}

std::vector<std::vector<FlintMod>> parseRREFRows(const std::string& polyPath) {
    std::ifstream in(polyPath);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open poly_relation file: " + polyPath);
    }

    std::vector<std::vector<FlintMod>> rows;
    bool inRREF = false;
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) continue;

        if (!inRREF) {
            if (line == "[rref]") {
                inRREF = true;
            }
            continue;
        }

        if (line.front() != '{') {
            continue;
        }

        const size_t l = line.find('{');
        const size_t r = line.rfind('}');
        if (l == std::string::npos || r == std::string::npos || r <= l) {
            throw std::runtime_error("Invalid rref row: " + line);
        }

        std::vector<FlintMod> row;
        std::string body = line.substr(l + 1, r - l - 1);
        std::stringstream ss(body);
        std::string tok;
        while (std::getline(ss, tok, ',')) {
            tok = trim(tok);
            if (tok.empty()) {
                row.emplace_back(0ULL);
            } else {
                row.emplace_back(static_cast<unsigned long long>(std::stoull(tok)));
            }
        }
        if (!row.empty()) {
            rows.push_back(std::move(row));
        }
    }

    if (rows.empty()) {
        throw std::runtime_error("No [rref] block rows found in poly_relation");
    }
    return rows;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 3 && argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <G_path> <poly_relation_path> [output_path]\n";
        std::cerr << "Default output_path: fi_solution\n";
        return 1;
    }

    const std::string gPath = argv[1];
    const std::string polyPath = argv[2];
    const std::string outPath = (argc == 4) ? argv[3] : "fi_solution";

    try {
        const int nuSize = inferNuSize(gPath);
        const auto G = parseSearchTargetFile(gPath, nuSize);
        const auto kv = parseHeaderKV(polyPath);

        mp_limb_t p = 2305843009213693951ULL;
        auto pIt = kv.find("p");
        if (pIt != kv.end()) {
            p = static_cast<mp_limb_t>(std::stoull(pIt->second));
        } else {
            std::cerr << "[fi_solver] Warning: '# p = ...' not found. Use default p=" << p << "\n";
        }
        FlintMod::set_modulus(p);

        int m = -1;
        auto mIt = kv.find("m");
        if (mIt != kv.end()) {
            m = std::stoi(mIt->second);
        }

        int numVarsFromHeader = -1;
        auto vIt = kv.find("variables");
        if (vIt != kv.end()) {
            numVarsFromHeader = std::stoi(vIt->second);
        }

        if (m < 0) {
            if (numVarsFromHeader > 0 &&
                static_cast<int>(G.size()) > 0 &&
                numVarsFromHeader % static_cast<int>(G.size()) == 0) {
                m = numVarsFromHeader / static_cast<int>(G.size()) - 1;
            } else {
                throw std::runtime_error("Cannot determine m from poly_relation header");
            }
        }

        std::vector<RelationVariable> variables;
        variables.reserve(G.size() * static_cast<size_t>(m + 1));
        for (const auto& g : G) {
            for (int k = 0; k <= m; ++k) {
                RelationVariable var;
                var.integral = g;
                var.k = k;
                variables.push_back(std::move(var));
            }
        }
        std::sort(variables.begin(), variables.end(), RelationVariableMoreComplexFirst{});

        auto rref = parseRREFRows(polyPath);
        const int nCols = static_cast<int>(rref.front().size());
        for (const auto& row : rref) {
            if (static_cast<int>(row.size()) != nCols) {
                throw std::runtime_error("RREF row size mismatch in poly_relation");
            }
        }
        if (nCols != static_cast<int>(variables.size())) {
            throw std::runtime_error(
                "Variable count mismatch: from G,m is " + std::to_string(variables.size()) +
                ", but RREF has " + std::to_string(nCols) + " columns");
        }

        RelationSearchResult<FlintMod> stage1;
        stage1.variables = std::move(variables);
        stage1.rrefMatrix = std::move(rref);
        stage1.pivotColumns.assign(stage1.rrefMatrix.size(), -1);
        std::vector<bool> isPivotCol(static_cast<size_t>(nCols), false);
        for (size_t i = 0; i < stage1.rrefMatrix.size(); ++i) {
            for (int col = 0; col < nCols; ++col) {
                if (stage1.rrefMatrix[i][static_cast<size_t>(col)] != FlintMod(0)) {
                    stage1.pivotColumns[i] = col;
                    isPivotCol[static_cast<size_t>(col)] = true;
                    break;
                }
            }
        }
        for (int col = 0; col < nCols; ++col) {
            if (!isPivotCol[static_cast<size_t>(col)]) {
                stage1.freeColumns.push_back(col);
            }
        }

        CoefficientRelationExpander<FlintMod> expander;
        const auto assignments = expander.expandAssignments(stage1);
        const auto fiRelations = expander.buildFIRelations(assignments);

        std::ofstream out(outPath);
        if (!out.is_open()) {
            throw std::runtime_error("Cannot open output file: " + outPath);
        }

        out << "# p = " << p << "\n";
        out << "# m = " << m << "\n";
        out << "# |G| = " << G.size() << "\n\n";
        RelationFormatter<FlintMod>::writeFIRelations(out, fiRelations);
        out << "\n";

        if (fiRelations.empty()) {
            out << "# FI variables = 0\n";
            out << "# FI pivot columns = 0\n";
            out << "# FI free columns = 0\n";
            out << "#MIs\n\n";
            out << "[fi_reductions]\n\n";
            out << "[fi_rref]\n";

            std::cerr << "[fi_solver] G size = " << G.size()
                      << ", m = " << m
                      << ", FI relations = 0"
                      << ", FI free columns = 0 (no free columns in stage I)"
                      << "\n";
            return 0;
        }

        FIReductionSearcher<FlintMod> fiSearcher(fiRelations);
        const auto fiResult = fiSearcher.search();

        RelationFormatter<FlintMod>::writeFIReductionSummary(out, fiResult);
        RelationFormatter<FlintMod>::writeFIMasterBasis(out, fiResult);
        out << "\n";
        RelationFormatter<FlintMod>::writeFIReductions(out, fiResult);
        out << "\n";
        RelationFormatter<FlintMod>::writeFIRREF(out, fiResult);

        std::cerr << "[fi_solver] G size = " << G.size()
                  << ", m = " << m
                  << ", FI relations = " << fiRelations.size()
                  << ", FI free columns = " << fiResult.freeColumns.size()
                  << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
