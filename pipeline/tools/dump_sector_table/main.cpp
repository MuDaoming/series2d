#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../../expand/include/fi_pipeline.hpp"

namespace {

std::string vecToString(const std::vector<int>& values) {
    std::ostringstream out;
    out << "{";
    for (size_t i = 0; i < values.size(); ++i) {
        if (i) out << ",";
        out << values[i];
    }
    out << "}";
    return out.str();
}

std::string boolString(bool value) {
    return value ? "1" : "0";
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <S_path> <config_path> <output_path>\n";
        return 1;
    }

    try {
        namespace fs = std::filesystem;
        fs::path sPath = fs::path(argv[1]);
        fs::path configPath = fs::path(argv[2]);
        fs::path outputPath = fs::path(argv[3]);

        if (!sPath.is_absolute()) sPath = fs::current_path() / sPath;
        if (!configPath.is_absolute()) configPath = fs::current_path() / configPath;
        if (!outputPath.is_absolute()) outputPath = fs::current_path() / outputPath;
        if (outputPath.has_parent_path()) {
            fs::create_directories(outputPath.parent_path());
        }

        InputConfig cfg = parseConfigFile(configPath.string());
        GiNaC::symbol X("X"), Y("Y");
        std::vector<std::vector<GiNaC::ex>> topS =
            parseMatrixFile(sPath.string(), X, Y);
        if (topS.empty() || topS.size() != topS[0].size()) {
            throw std::runtime_error("S matrix must be non-empty and square");
        }
        if (static_cast<int>(topS.size()) != cfg.B + cfg.N) {
            throw std::runtime_error("S dimension mismatch: expected B+N");
        }

        Family<GiNaC::ex, GiNaC::ex, GiNaC::ex> family(topS, cfg.N, cfg.B, X, Y);
        if (!cfg.sectorMapPath.empty()) {
            family.setSectorMap(SectorMap::fromFile(cfg.sectorMapPath, cfg.N));
        }

        std::ofstream out(outputPath.string());
        if (!out.is_open()) {
            throw std::runtime_error("Cannot open output file: " + outputPath.string());
        }

        out << "sector\tcase\trepresentative\tis_representative\thas_cut\tprop_count\tbranch_support\n";
        const auto& branchIndices = family.getBranchIndices();
        const int total = 1 << cfg.N;
        for (int idx = total - 1; idx >= 1; --idx) {
            std::vector<int> sector = family.secvecFromIdx(idx);
            const int caseNum = family.getCase(sector);
            if (caseNum < 0) continue;

            std::vector<int> representative = family.canonicalizeSector(sector);
            const bool isRepresentative = (sector == representative);
            const bool hasCut = isRepresentative && caseNum == 0;

            std::vector<int> branchSupport(cfg.B, 0);
            int propCount = 0;
            for (int i = 0; i < cfg.N; ++i) {
                if (!sector[i]) continue;
                ++propCount;
                branchSupport[branchIndices[i]] = 1;
            }

            out << vecToString(sector) << '\t'
                << caseNum << '\t'
                << vecToString(representative) << '\t'
                << boolString(isRepresentative) << '\t'
                << boolString(hasCut) << '\t'
                << propCount << '\t'
                << vecToString(branchSupport) << '\n';
        }

        std::cout << "Wrote sector table to: " << outputPath.string() << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
