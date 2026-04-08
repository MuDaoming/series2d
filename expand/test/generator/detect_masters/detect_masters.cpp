#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../../include/fi_pipeline.hpp"

int main(int argc, char** argv) {
    if (argc != 3 && argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <S_path> <config_path> [output_path]\n";
        std::cerr << "Default output_path: masters\n";
        return 1;
    }

    try {
        namespace fs = std::filesystem;
        fs::path sPath = fs::path(argv[1]);
        fs::path configPath = fs::path(argv[2]);
        fs::path outputPath = (argc == 4) ? fs::path(argv[3]) : fs::path("masters");

        if (!sPath.is_absolute()) sPath = fs::current_path() / sPath;
        if (!configPath.is_absolute()) configPath = fs::current_path() / configPath;
        if (!outputPath.is_absolute()) outputPath = fs::current_path() / outputPath;
        if (outputPath.has_parent_path()) {
            fs::create_directories(outputPath.parent_path());
        }

        InputConfig cfg = parseConfigFile(configPath.string());
        GiNaC::symbol X("X"), Y("Y");
        std::vector<std::vector<GiNaC::ex>> topS = parseMatrixFile(sPath.string(), X, Y);
        if (topS.empty() || topS.size() != topS[0].size()) {
            throw std::runtime_error("S matrix must be non-empty and square");
        }
        if (static_cast<int>(topS.size()) != cfg.B + cfg.N) {
            throw std::runtime_error("S dimension mismatch: expected B+N");
        }

        Family<GiNaC::ex, GiNaC::ex, GiNaC::ex> family(topS, cfg.N, cfg.B, X, Y);
        const std::vector<int>& masterIdxs = family.getMasterIdxs();

        std::ofstream out(outputPath.string());
        if (!out.is_open()) {
            throw std::runtime_error("Cannot open output file: " + outputPath.string());
        }

        for (int idx : masterIdxs) {
            std::vector<int> nu = family.secvecFromIdx(idx);
            out << "{";
            for (size_t i = 0; i < nu.size(); ++i) {
                out << nu[i];
                if (i + 1 < nu.size()) out << ",";
            }
            out << "}\n";
        }

        std::cout << "Wrote " << masterIdxs.size()
                  << " master nus to: " << outputPath.string() << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
