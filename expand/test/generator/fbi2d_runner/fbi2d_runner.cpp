#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "../../include/converter.hpp"
#include "../../include/integrand_expander.hpp"
#include "../../include/io.hpp"
#include "../../include/series_solver.hpp"

#include <ginac/ginac.h>

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <S_path> <config_path> <target_path> <output_path>\n";
        return 1;
    }

    try {
        const std::string sPath = argv[1];
        const std::string configPath = argv[2];
        const std::string targetPath = argv[3];
        const std::string outputPath = argv[4];

        InputConfig cfg = parseConfigFile(configPath);
        TargetConfig targetCfg = parseTargetFile(targetPath, cfg.N);
        FlintMod::set_modulus(cfg.p);

        GiNaC::symbol X("X"), Y("Y");
        auto topS = parseMatrixFile(sPath, X, Y);
        if (topS.empty() || topS.size() != topS[0].size()) {
            throw std::runtime_error("S matrix must be non-empty and square");
        }
        if (static_cast<int>(topS.size()) != cfg.B + cfg.N) {
            throw std::runtime_error("S dimension mismatch: expected B+N");
        }

        const int numBranch = cfg.B;
        const int numProps = static_cast<int>(topS.size()) - numBranch;
        Family<GiNaC::ex, GiNaC::ex, GiNaC::ex> ginacFamily(topS, numProps, numBranch, X, Y);
        auto family = convertFamily(ginacFamily, X, Y);

        const FlintMod feynmanD(cfg.d);
        const FlintMod shiftA(cfg.a);
        const FlintMod shiftB(cfg.b);
        family.setMasterDelta(feynmanD);

        SeriesSolver<Rational<FlintMod>, Polynomial<FlintMod>, FlintMod> solver(family, cfg.deg);
        solver.setReduceMode(cfg.reduceMode);
        if (cfg.bc.size() != static_cast<size_t>(solver.getNumMaster())) {
            throw std::runtime_error("config.bc size mismatch with number of masters");
        }
        for (int i = 0; i < solver.getNumMaster(); ++i) {
            solver.setMasterBoundary(i, FlintMod(cfg.bc[i]));
        }

        const int numLoops = 2;
        const FlintMod fbiDelta = FlintMod(numLoops) * feynmanD / FlintMod(2);
        IntegrandExpander<Rational<FlintMod>, Polynomial<FlintMod>, FlintMod> expander(
            solver, numLoops, cfg.deg, feynmanD, shiftA, shiftB);
        Redefinition<Polynomial<FlintMod>, FlintMod> redef(numLoops, fbiDelta, expander.getShiftedU());
        solver.setRedefinition(&redef);

        std::ofstream out(outputPath);
        if (!out.is_open()) {
            throw std::runtime_error("Cannot open output file: " + outputPath);
        }

        out << "{\n";
        for (size_t i = 0; i < targetCfg.nus.size(); ++i) {
            const auto& nu = targetCfg.nus[i];
            const auto& series = solver.getFBISeries(nu, fbiDelta, cfg.deg);
            out << "  {" << series.toString() << "}";
            if (i + 1 < targetCfg.nus.size()) out << ",";
            out << "\n";
        }
        out << "}\n";

        std::cout << "FBI 2D series finished.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}

