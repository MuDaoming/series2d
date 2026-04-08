// poly_relation_searcher.cpp
//
// Stage I of IBP search: find polynomial relations among FI(delta).
//
// Usage:
//   ./stage1_runner <config> <G_file> <series_list> <output>
//
// config format:
//   N   = 6                      # nu vector length
//   deg = 200                    # how many orders to use (must be <= series deg)
//   m   = 4                      # polynomial degree bound for searched relations
//   p   = 2305843009213693951    # prime modulus
//
// G_file:
//   one integral per line, same format as expand target, e.g.
//     {1, 1, 1, 1, 1, 1}
//   G must be a subset of each series' own target.
//
// series_list:
//   one boundary condition per line:
//     path/to/series  path/to/target/of/series
//
// output:
//   summary + [relations] + [rref]

#include <cctype>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "ff_type.hpp"
#include "io.hpp"
#include "relation_formatter.hpp"
#include "relation_searcher.hpp"
#include "relation_types.hpp"

// ---------------------------------------------------------------------------
// Config
// ---------------------------------------------------------------------------

struct Stage1Config {
    int nuSize = 0;
    int deg    = 0;
    int m      = 0;
    mp_limb_t p = 0;
};

static std::string s1Trim(const std::string& s) {
    size_t l = 0;
    while (l < s.size() && std::isspace(static_cast<unsigned char>(s[l]))) ++l;
    size_t r = s.size();
    while (r > l && std::isspace(static_cast<unsigned char>(s[r - 1]))) --r;
    return s.substr(l, r - l);
}

static Stage1Config parseConfig(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open())
        throw std::runtime_error("Cannot open config: " + path);

    std::unordered_map<std::string, std::string> kv;
    std::string line;
    while (std::getline(in, line)) {
        line = s1Trim(line);
        if (line.empty()) continue;
        size_t eq = line.find('=');
        if (eq == std::string::npos) continue;
        kv[s1Trim(line.substr(0, eq))] = s1Trim(line.substr(eq + 1));
    }

    auto need = [&](const std::string& key) -> std::string {
        auto it = kv.find(key);
        if (it == kv.end())
            throw std::runtime_error("Missing key in config: " + key);
        return it->second;
    };

    Stage1Config cfg;
    cfg.nuSize = std::stoi(need("N"));
    cfg.deg    = std::stoi(need("deg"));
    cfg.m      = std::stoi(need("m"));
    cfg.p      = static_cast<mp_limb_t>(std::stoull(need("p")));
    return cfg;
}

// ---------------------------------------------------------------------------
// Series list file: each line is "seriesPath targetPath"
// ---------------------------------------------------------------------------

static std::vector<std::pair<std::string, std::string>>
parseSeriesListFile(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open())
        throw std::runtime_error("Cannot open series list: " + path);

    std::vector<std::pair<std::string, std::string>> entries;
    std::string line;
    while (std::getline(in, line)) {
        line = s1Trim(line);
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::string seriesPath, targetPath;
        iss >> seriesPath >> targetPath;
        if (seriesPath.empty() || targetPath.empty())
            throw std::runtime_error("Bad series_list line: " + line);
        entries.emplace_back(std::move(seriesPath), std::move(targetPath));
    }
    if (entries.empty())
        throw std::runtime_error("series_list file is empty: " + path);
    return entries;
}

// ---------------------------------------------------------------------------
// Read a series file, filtering to only the integrals in G and truncating
// to the first deg+1 coefficients.
// ---------------------------------------------------------------------------

static std::vector<SeriesSample<FlintMod>>
parseSeriesForG(const std::string& seriesPath,
                const std::string& seriesTargetPath,
                const std::vector<IntegralLabel>& G,
                int deg,
                int bcIndex)
{
    // 1. Read the full target for this series.
    const int nuSize = static_cast<int>(G[0].nu.size());
    const auto fullTarget = parseSearchTargetFile(seriesTargetPath, nuSize);

    // 2. Build set: nu -> index in G.
    std::map<std::vector<int>, int> nuToGIdx;
    for (int i = 0; i < static_cast<int>(G.size()); ++i)
        nuToGIdx[G[i].nu] = i;

    // 3. Map each fullTarget position to a G index (-1 if absent).
    std::vector<int> fullToG(fullTarget.size(), -1);
    for (int i = 0; i < static_cast<int>(fullTarget.size()); ++i) {
        auto it = nuToGIdx.find(fullTarget[i].nu);
        if (it != nuToGIdx.end())
            fullToG[i] = it->second;
    }

    // 4. Read series file line by line; extract relevant rows.
    std::ifstream in(seriesPath);
    if (!in.is_open())
        throw std::runtime_error("Cannot open series file: " + seriesPath);

    std::vector<SeriesSample<FlintMod>> result(G.size());
    for (int i = 0; i < static_cast<int>(G.size()); ++i) {
        result[i].label.integral = G[i];
        result[i].label.bcIndex  = bcIndex;
    }

    std::string line;
    int lineIdx = 0;
    while (std::getline(in, line)) {
        line = s1Trim(line);
        if (line.empty()) continue;
        if (lineIdx >= static_cast<int>(fullTarget.size()))
            throw std::runtime_error("Too many lines in series file: " + seriesPath);

        const int gIdx = fullToG[lineIdx];
        if (gIdx >= 0) {
            // Parse the line: {c0,c1,...}
            size_t l = line.find('{');
            size_t r = line.rfind('}');
            if (l == std::string::npos || r == std::string::npos || r <= l)
                throw std::runtime_error("Bad series line in: " + seriesPath);

            std::string body = line.substr(l + 1, r - l - 1);
            std::istringstream iss(body);
            std::string tok;
            int count = 0;
            auto& coeffs = result[gIdx].coeffs;
            coeffs.reserve(deg + 1);
            while (count <= deg && std::getline(iss, tok, ',')) {
                tok = s1Trim(tok);
                if (tok.empty()) continue;
                coeffs.emplace_back(
                    static_cast<unsigned long long>(std::stoull(tok)));
                ++count;
            }
            if (static_cast<int>(coeffs.size()) != deg + 1)
                throw std::runtime_error(
                    "Series in " + seriesPath + " has fewer than deg+1=" +
                    std::to_string(deg + 1) + " coefficients");
        }
        ++lineIdx;
    }

    if (lineIdx != static_cast<int>(fullTarget.size()))
        throw std::runtime_error("Series line count mismatch in: " + seriesPath);

    // Verify all G integrals were matched.
    for (int i = 0; i < static_cast<int>(G.size()); ++i) {
        if (result[i].coeffs.empty())
            throw std::runtime_error(
                "Integral in G not found in series target " + seriesTargetPath +
                ": " + nuToString(G[i].nu));
    }

    return result;
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <config> <G_file> <series_list> <output>\n";
        return 1;
    }

    const std::string configPath     = argv[1];
    const std::string gPath          = argv[2];
    const std::string seriesListPath = argv[3];
    const std::string outputPath     = argv[4];

    try {
        const Stage1Config cfg = parseConfig(configPath);
        FlintMod::set_modulus(cfg.p);

        const auto G          = parseSearchTargetFile(gPath, cfg.nuSize);
        const auto seriesList = parseSeriesListFile(seriesListPath);
        const int numBC       = static_cast<int>(seriesList.size());

        std::cerr << "[stage1] G size = " << G.size()
                  << ", numBC = " << numBC
                  << ", deg = " << cfg.deg
                  << ", m = " << cfg.m << "\n";

        SearchInput<FlintMod> input;
        input.degreeD         = cfg.deg;
        input.maxDeltaDegreeM = cfg.m;
        input.numFBIMasters   = numBC;
        input.targets         = G;

        for (int b = 0; b < numBC; ++b) {
            std::cerr << "[stage1] reading bc" << (b + 1)
                      << ": " << seriesList[b].first << "\n";
            auto samples = parseSeriesForG(
                seriesList[b].first,
                seriesList[b].second,
                G, cfg.deg, b);
            input.samples.insert(input.samples.end(),
                                 samples.begin(), samples.end());
        }

        std::cerr << "[stage1] building and solving A^(delta)...\n";
        RelationSearcher<FlintMod> searcher(input);
        const auto result = searcher.search();

        std::ofstream out(outputPath);
        if (!out.is_open())
            throw std::runtime_error("Cannot open output: " + outputPath);

        out << "# p = " << cfg.p << "\n";
        out << "# m = " << cfg.m << "\n";
        RelationFormatter<FlintMod>::writeSummary(out, result);
        out << "\n";
        RelationFormatter<FlintMod>::writeRelations(out, result);
        out << "\n";
        RelationFormatter<FlintMod>::writeRREF(out, result);

        std::cerr << "[stage1] done. free columns = "
                  << result.freeColumns.size() << "\n";
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
