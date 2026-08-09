#include <cctype>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>

#include "ff_type.hpp"
#include "formatter.hpp"
#include "io.hpp"
#include "reducer.hpp"
#include "sector_map.hpp"

namespace {

std::string trimLocal(const std::string& s) {
    size_t l = 0;
    while (l < s.size() && std::isspace(static_cast<unsigned char>(s[l]))) ++l;
    size_t r = s.size();
    while (r > l && std::isspace(static_cast<unsigned char>(s[r - 1]))) --r;
    return s.substr(l, r - l);
}

template<typename T>
Polynomial1D<T> parsePolynomialLocal(const std::string& text) {
    const size_t l = text.find('{');
    const size_t r = text.rfind('}');
    if (l == std::string::npos || r == std::string::npos || r <= l) {
        throw std::runtime_error("Invalid polynomial: " + text);
    }
    std::vector<T> coeffs;
    std::stringstream ss(text.substr(l + 1, r - l - 1));
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        tok = trimLocal(tok);
        if (!tok.empty()) coeffs.emplace_back(static_cast<unsigned long long>(std::stoull(tok)));
    }
    if (coeffs.empty()) coeffs.emplace_back(T(0));
    return Polynomial1D<T>(std::move(coeffs));
}

template<typename T>
std::vector<SectorReduction<T>> parseSeedReductions(const std::string& path, int nuSize) {
    std::ifstream in(path);
    if (!in.is_open()) throw std::runtime_error("Cannot open seed reduction file: " + path);

    std::vector<SectorReduction<T>> out;
    SectorReduction<T> current;
    bool inBlock = false;
    bool keep = false;

    auto flush = [&]() {
        if (inBlock && keep && !current.failed && !current.isZero) {
            out.push_back(std::move(current));
        }
        current = SectorReduction<T>{};
        inBlock = false;
        keep = false;
    };

    std::string line;
    while (std::getline(in, line)) {
        line = trimLocal(line);
        if (line.empty()) {
            flush();
            continue;
        }
        if (line[0] == '#' || line == "[sector_reductions]") continue;
        if (line == "[global_reductions]") {
            flush();
            break;
        }
        if (line.rfind("sector=", 0) == 0) {
            flush();
            inBlock = true;
            current.sector = parseSectorId(line.substr(7), nuSize);
            continue;
        }
        if (!inBlock) continue;
        if (line.rfind("object=", 0) == 0) {
            current.object = parseObjectLabel(line.substr(7), nuSize);
        } else if (line == "status=success") {
            keep = true;
        } else if (line == "status=failed") {
            current.failed = true;
            keep = false;
        } else if (line == "zero") {
            current.isZero = true;
        } else if (line.rfind("free_master=", 0) == 0) {
            current.isFreeMaster = line.substr(12) == "1";
        } else if (line.rfind("den=", 0) == 0) {
            current.denominator = parsePolynomialLocal<T>(line.substr(4));
        } else if (line.rfind("term ", 0) == 0) {
            const size_t eq = line.find('=');
            if (eq == std::string::npos) {
                throw std::runtime_error("Invalid seed term line: " + line);
            }
            ReductionTerm<T> term;
            term.master = parseObjectLabel(line.substr(5, eq - 5), nuSize);
            term.numerator = parsePolynomialLocal<T>(line.substr(eq + 1));
            current.terms.push_back(std::move(term));
        }
    }
    flush();
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 4 && argc != 5 && argc != 6) {
        std::cerr << "Usage: " << argv[0]
                  << " <config_path> <sector_series_list_path> <output_path> [object_list_path] [seed_reduction_path]\n";
        return 1;
    }

    const std::string configPath = argv[1];
    const std::string sectorListPath = argv[2];
    const std::string outputPath = argv[3];
    const bool hasObjectList = argc >= 5;
    const bool hasSeedList = argc == 6;
    const std::string objectListPath = (argc >= 5) ? argv[4] : "";
    const std::string seedReductionPath = hasSeedList ? argv[5] : "";

    try {
        BLSectorConfig config = parseBLSectorConfig(configPath);
        FlintMod::set_modulus(config.prime);

        const auto entries = parseSectorSeriesList(sectorListPath, config.nuSize);
        SeriesStore<FlintMod> series;
        MasterData masters;
        loadSectorData(entries, config.degreeD, config.nuSize, series, masters);

        const auto sectors = series.sectors();
        SectorTree tree = config.sectorMapPath.empty()
            ? SectorTree(sectors)
            : SectorTree(sectors, SectorMap::fromFile(config.sectorMapPath, config.nuSize));

        std::vector<ObjectLabel> objects =
            hasObjectList ? parseObjectList(objectListPath, config.nuSize) : series.objects();

        BLSectorReducer<FlintMod> reducer(config, tree, series, masters);
        auto writeSnapshot = [&](const std::vector<SectorReduction<FlintMod>>& reductions) {
            std::ofstream out(outputPath);
            if (!out.is_open()) {
                throw std::runtime_error("Cannot open output file: " + outputPath);
            }
            writeReductions(out, config, tree, masters, reductions);
        };
        const auto seedReductions = hasSeedList
            ? parseSeedReductions<FlintMod>(seedReductionPath, config.nuSize)
            : std::vector<SectorReduction<FlintMod>>{};
        writeSnapshot(seedReductions);
        const auto reductions = reducer.reduceAll(objects, seedReductions, writeSnapshot);
        writeSnapshot(reductions);
    } catch (const std::exception& e) {
        std::cerr << "[bl_sector_reducer] Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
