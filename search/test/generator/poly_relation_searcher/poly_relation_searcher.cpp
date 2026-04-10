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
//   ncheck = 1                   # optional, default 1
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
#include "relation_matrix_builder.hpp"
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
    int ncheck = 1;
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
    {
        auto it = kv.find("ncheck");
        if (it != kv.end()) cfg.ncheck = std::stoi(it->second);
    }
    if (cfg.ncheck < 0)
        throw std::runtime_error("ncheck must be >= 0");
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

static std::vector<SeriesSample<FlintMod>>
truncateSamplesForTrain(const std::vector<SeriesSample<FlintMod>>& samples,
                        int trainDeg) {
    std::vector<SeriesSample<FlintMod>> out = samples;
    for (auto& s : out) {
        if (static_cast<int>(s.coeffs.size()) < trainDeg + 1) {
            throw std::runtime_error("series coeff count smaller than trainDeg+1");
        }
        s.coeffs.resize(trainDeg + 1);
    }
    return out;
}

static bool checkNullspaceShrink(
    const std::vector<std::vector<FlintMod>>& trainRREF,
    const std::vector<int>& pivotColumns,
    const std::vector<std::vector<FlintMod>>& checkRows) {
    const FlintMod zero(0ULL);
    if (trainRREF.empty() || pivotColumns.empty()) {
        for (const auto& row : checkRows) {
            for (const auto& x : row) {
                if (x != zero) return true;
            }
        }
        return false;
    }

    for (const auto& rawRow : checkRows) {
        std::vector<FlintMod> row = rawRow;
        for (int r = 0; r < static_cast<int>(trainRREF.size()); ++r) {
            const int pc = pivotColumns[r];
            if (pc < 0 || pc >= static_cast<int>(row.size())) continue;
            const FlintMod factor = row[pc];
            if (factor == zero) continue;
            for (int c = pc; c < static_cast<int>(row.size()); ++c) {
                row[c] = row[c] - factor * trainRREF[r][c];
            }
        }
        for (const auto& x : row) {
            if (x != zero) {
                return true;
            }
        }
    }
    return false;
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
        if (cfg.ncheck > cfg.deg) {
            throw std::runtime_error("ncheck must satisfy ncheck <= deg");
        }
        const int trainDeg = cfg.deg - cfg.ncheck;
        FlintMod::set_modulus(cfg.p);

        const auto G          = parseSearchTargetFile(gPath, cfg.nuSize);
        const auto seriesList = parseSeriesListFile(seriesListPath);
        const int numBC       = static_cast<int>(seriesList.size());

        std::cerr << "[stage1] G size = " << G.size()
                  << ", numBC = " << numBC
                  << ", deg = " << cfg.deg
                  << ", m = " << cfg.m
                  << ", ncheck = " << cfg.ncheck << "\n";

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

        SearchInput<FlintMod> trainInput = input;
        trainInput.degreeD = trainDeg;
        trainInput.samples = truncateSamplesForTrain(input.samples, trainDeg);

        std::cerr << "[stage1] solving train window: [0.." << trainDeg << "]\n";
        RelationSearcher<FlintMod> searcher(trainInput);
        const auto result = searcher.search();

        bool isNullShrink = false;
        if (cfg.ncheck > 0) {
            std::cerr << "[stage1] checking held-out window: ["
                      << (trainDeg + 1) << ".." << cfg.deg << "]\n";
            RelationMatrixBuilder<FlintMod> fullBuilder(input);
            const auto fullMatrix = fullBuilder.buildMatrix(result.variables);
            std::vector<std::vector<FlintMod>> checkRows;
            checkRows.reserve(static_cast<size_t>(numBC * cfg.ncheck));
            for (int b = 0; b < numBC; ++b) {
                for (int n = trainDeg + 1; n <= cfg.deg; ++n) {
                    const int rowIdx = b * (cfg.deg + 1) + n;
                    checkRows.push_back(fullMatrix[rowIdx]);
                }
            }
            isNullShrink = checkNullspaceShrink(
                result.rrefMatrix, result.pivotColumns, checkRows);
        } else {
            std::cerr << "[stage1] ncheck=0, skip held-out check.\n";
        }

        std::ofstream out(outputPath);
        if (!out.is_open())
            throw std::runtime_error("Cannot open output: " + outputPath);

        out << "# p = " << cfg.p << "\n";
        out << "# m = " << cfg.m << "\n";
        out << "# ncheck = " << cfg.ncheck << "\n";
        out << "# train_deg = " << trainDeg << "\n";
        out << "# nullspace_shrink = " << (isNullShrink ? 1 : 0) << "\n";
        if (isNullShrink) {
            out << "[status]\n";
            out << "no-certified-solution\n";
            std::cerr << "[stage1] shrink detected on held-out rows; no certified solution.\n";
            out << "\n";
        }
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
