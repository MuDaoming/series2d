#include <cctype>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "absorbing_relation_searcher.hpp"
#include "ff_type.hpp"
#include "io.hpp"
#include "relation_certifier.hpp"
#include "relation_types.hpp"

namespace {

struct Config {
    int nuSize = 0;
    int deg = 0;
    int m = 0;
    int ncheck = 1;
    mp_limb_t p = 0;
};

std::string trim(const std::string& s) {
    size_t l = 0;
    while (l < s.size() && std::isspace(static_cast<unsigned char>(s[l]))) ++l;
    size_t r = s.size();
    while (r > l && std::isspace(static_cast<unsigned char>(s[r - 1]))) --r;
    return s.substr(l, r - l);
}

Config parseConfig(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) throw std::runtime_error("Cannot open config: " + path);
    std::unordered_map<std::string, std::string> kv;
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) continue;
        const size_t eq = line.find('=');
        if (eq == std::string::npos) continue;
        kv[trim(line.substr(0, eq))] = trim(line.substr(eq + 1));
    }
    auto need = [&](const std::string& key) -> std::string {
        auto it = kv.find(key);
        if (it == kv.end()) throw std::runtime_error("Missing key in config: " + key);
        return it->second;
    };
    Config cfg;
    cfg.nuSize = std::stoi(need("N"));
    cfg.deg = std::stoi(need("deg"));
    cfg.m = std::stoi(need("m"));
    if (auto it = kv.find("ncheck"); it != kv.end()) cfg.ncheck = std::stoi(it->second);
    cfg.p = static_cast<mp_limb_t>(std::stoull(need("p")));
    return cfg;
}

std::vector<std::pair<std::string, std::string>> parseSeriesList(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) throw std::runtime_error("Cannot open series list: " + path);
    std::vector<std::pair<std::string, std::string>> out;
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::string series, target;
        iss >> series >> target;
        if (series.empty() || target.empty()) throw std::runtime_error("Bad series list line: " + line);
        out.emplace_back(series, target);
    }
    if (out.empty()) throw std::runtime_error("series list is empty");
    return out;
}

std::vector<SeriesSample<FlintMod>> parseSeriesForG(
    const std::string& seriesPath,
    const std::string& targetPath,
    const std::vector<IntegralLabel>& G,
    int deg,
    int bcIndex) {
    const int nuSize = static_cast<int>(G.front().nu.size());
    const auto fullTarget = parseSearchTargetFile(targetPath, nuSize);
    std::map<IntegralLabel, int, IntegralLabelLess> labelToG;
    for (int i = 0; i < static_cast<int>(G.size()); ++i) labelToG[G[i]] = i;
    std::vector<int> fullToG(fullTarget.size(), -1);
    for (int i = 0; i < static_cast<int>(fullTarget.size()); ++i) {
        auto it = labelToG.find(fullTarget[i]);
        if (it != labelToG.end()) fullToG[i] = it->second;
    }

    std::ifstream in(seriesPath);
    if (!in.is_open()) throw std::runtime_error("Cannot open series: " + seriesPath);
    std::vector<SeriesSample<FlintMod>> result(G.size());
    for (int i = 0; i < static_cast<int>(G.size()); ++i) {
        result[i].label.integral = G[i];
        result[i].label.bcIndex = bcIndex;
    }
    std::string line;
    int lineIdx = 0;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) continue;
        if (lineIdx >= static_cast<int>(fullTarget.size())) {
            throw std::runtime_error("Too many lines in series: " + seriesPath);
        }
        const int gIdx = fullToG[lineIdx];
        if (gIdx >= 0) {
            const size_t l = line.find('{');
            const size_t r = line.rfind('}');
            if (l == std::string::npos || r == std::string::npos || r <= l) {
                throw std::runtime_error("Bad series line in: " + seriesPath);
            }
            std::string body = line.substr(l + 1, r - l - 1);
            std::stringstream ss(body);
            std::string tok;
            auto& coeffs = result[gIdx].coeffs;
            coeffs.reserve(deg + 1);
            while (static_cast<int>(coeffs.size()) <= deg && std::getline(ss, tok, ',')) {
                tok = trim(tok);
                if (!tok.empty()) coeffs.emplace_back(static_cast<unsigned long long>(std::stoull(tok)));
            }
            if (static_cast<int>(coeffs.size()) != deg + 1) {
                throw std::runtime_error("Series has fewer than deg+1 coefficients: " + seriesPath);
            }
        }
        ++lineIdx;
    }
    if (lineIdx != static_cast<int>(fullTarget.size())) {
        throw std::runtime_error("Series line count mismatch: " + seriesPath);
    }
    for (int i = 0; i < static_cast<int>(G.size()); ++i) {
        if (result[i].coeffs.empty()) {
            throw std::runtime_error("G integral not found in target: " + integralLabelToString(G[i]));
        }
    }
    return result;
}

void writeResult(
    const std::string& path,
    const Config& cfg,
    const AbsorbingSearchResult<FlintMod>& result) {
    std::ofstream out(path);
    if (!out.is_open()) throw std::runtime_error("Cannot open output: " + path);
    out << "# format = absorbing_poly_relation_v1\n";
    out << "# p = " << cfg.p << "\n";
    out << "# max_m = " << result.maxDeltaDegreeM << "\n";
    out << "# train_deg = " << result.trainDegree << "\n";
    out << "# check_start = " << result.checkStart << "\n";
    out << "# check_end = " << result.checkEnd << "\n";
    out << "# integrals = " << result.integrals.size() << "\n";
    out << "# relations = " << result.relations.size() << "\n\n";

    out << "[integrals]\n";
    for (int i = 0; i < static_cast<int>(result.integrals.size()); ++i) {
        out << i << " " << integralLabelToString(result.integrals[i]) << "\n";
    }
    out << "\n[relations]\n";
    for (int r = 0; r < static_cast<int>(result.relations.size()); ++r) {
        const auto& rel = result.relations[r];
        out << "relation " << r
            << " lead " << rel.leadIntegralId << " " << rel.leadDeltaPower
            << " max " << rel.maxDeltaPower
            << " terms " << rel.terms.size() << "\n";
        for (const auto& term : rel.terms) {
            out << term.integralId << " " << term.deltaPower << " " << term.coeff << "\n";
        }
    }
    out << "\n[leading_rules]\n";
    for (const auto& rule : result.leadingRules) {
        out << rule.integralId << " "
            << rule.leadDeltaPower << " "
            << rule.relationId << " "
            << rule.relationMaxDeltaPower << "\n";
    }
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <config_path> <G_path> <series_list_path> <output_path>\n";
        return 1;
    }
    try {
        const Config cfg = parseConfig(argv[1]);
        const DegreeWindow window = makeDegreeWindow(cfg.deg, cfg.ncheck);
        FlintMod::set_modulus(cfg.p);

        const auto G = parseSearchTargetFile(argv[2], cfg.nuSize);
        const auto seriesList = parseSeriesList(argv[3]);

        SearchInput<FlintMod> input;
        input.degreeD = cfg.deg;
        input.maxDeltaDegreeM = cfg.m;
        input.numFBIMasters = static_cast<int>(seriesList.size());
        input.targets = G;
        for (int b = 0; b < static_cast<int>(seriesList.size()); ++b) {
            auto samples = parseSeriesForG(seriesList[b].first, seriesList[b].second, G, cfg.deg, b);
            input.samples.insert(input.samples.end(), samples.begin(), samples.end());
        }

        AbsorbingSearchOptions options;
        options.trainDegree = window.trainDeg;
        options.checkStart = window.checkStart;
        options.checkEnd = window.checkEnd;
        AbsorbingRelationSearcher<FlintMod> searcher(input, options);
        const auto result = searcher.search();
        writeResult(argv[4], cfg, result);
        std::cerr << "[absorbing_poly_searcher] G size = " << G.size()
                  << ", relations = " << result.relations.size()
                  << ", max_m = " << cfg.m << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
