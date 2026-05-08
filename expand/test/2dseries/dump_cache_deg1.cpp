#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../include/converter.hpp"
#include "../../include/integrand_expander.hpp"
#include "../../include/io.hpp"
#include "../../include/series_solver.hpp"

#include <ginac/ginac.h>

static void writeText(const std::string& path, const std::string& s) {
    std::ofstream out(path);
    if (!out.is_open()) throw std::runtime_error("Cannot open output file: " + path);
    out << s;
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <work_dir>\n";
        return 1;
    }

    try {
        const std::string workDir = argv[1];
        const std::string sPath = workDir + "/S";
        const std::string configPath = workDir + "/config";
        const std::string targetPath = workDir + "/target";

        InputConfig cfg = parseConfigFile(configPath);
        cfg.deg = 1;
        FlintMod::set_modulus(cfg.p);

        GiNaC::symbol X("X"), Y("Y");
        auto topS = parseMatrixFile(sPath, X, Y);
        if (topS.empty() || topS.size() != topS[0].size()) {
            throw std::runtime_error("S matrix must be non-empty and square");
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

        TargetConfig targetCfg = parseTargetFile(targetPath, cfg.N);
        for (const auto& nu : targetCfg.nus) {
            (void)solver.getFBISeries(nu, fbiDelta, cfg.deg);
        }

        std::ostringstream cacheOss;
        {
            auto* old = std::cout.rdbuf(cacheOss.rdbuf());
            solver.printAllCache();
            std::cout.rdbuf(old);
        }
        const std::string cacheText = cacheOss.str();
        writeText(workDir + "/cache_all_deg1.txt", cacheText);

        std::istringstream in(cacheText);
        std::ostringstream keyOut;
        std::string line;
        while (std::getline(in, line)) {
            auto p = line.find("nu = {");
            if (p == std::string::npos) continue;
            auto q = line.find("}, delta = ", p);
            if (q == std::string::npos) continue;
            std::string nu = line.substr(p + 6, q - (p + 6));
            std::string delta = line.substr(q + 11);
            keyOut << "{" << nu << "};" << delta << "\n";
        }
        writeText(workDir + "/cache_keys_deg1.txt", keyOut.str());

        std::cout << "Done. Wrote:\n"
                  << "  " << workDir + "/cache_all_deg1.txt" << "\n"
                  << "  " << workDir + "/cache_keys_deg1.txt" << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
