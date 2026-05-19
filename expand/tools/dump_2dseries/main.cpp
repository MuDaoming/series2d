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

namespace {

struct CliOptions {
    std::string sPath;
    std::string configPath;
    std::string targetPath;
    std::string outputPath;
};

void writeText(const std::string& path, const std::string& s) {
    std::ofstream out(path);
    if (!out.is_open()) throw std::runtime_error("Cannot open output file: " + path);
    out << s;
}

CliOptions parseArgs(int argc, char** argv) {
    if (argc != 5) {
        throw std::runtime_error(
            "Usage: dump_2dseries <S_path> <config_path> <target_path> <output_path>");
    }

    CliOptions opts;
    opts.sPath = argv[1];
    opts.configPath = argv[2];
    opts.targetPath = argv[3];
    opts.outputPath = argv[4];
    return opts;
}

} // namespace

int main(int argc, char** argv) {
    try {
        const CliOptions opts = parseArgs(argc, argv);

        InputConfig cfg = parseConfigFile(opts.configPath);
        TargetConfig targetCfg = parseTargetFile(opts.targetPath, cfg.N);
        FlintMod::set_modulus(cfg.p);

        GiNaC::symbol X("X"), Y("Y");
        auto topS = parseMatrixFile(opts.sPath, X, Y);
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

        if (cfg.print2DMode == "target") {
            std::ofstream out(opts.outputPath);
            if (!out.is_open()) {
                throw std::runtime_error("Cannot open output file: " + opts.outputPath);
            }

            out << "{\n";
            for (size_t i = 0; i < targetCfg.targets.size(); ++i) {
                const auto& nu = targetCfg.targets[i].nu;
                const auto& series = solver.getFBISeries(nu, fbiDelta, cfg.deg);
                out << "  {" << series.toString() << "}";
                if (i + 1 < targetCfg.targets.size()) out << ",";
                out << "\n";
            }
            out << "}\n";

            std::cout << "Wrote " << opts.outputPath << "\n";
            return 0;
        }

        for (const auto& tag : targetCfg.targets) {
            (void)solver.getFBISeries(tag.nu, fbiDelta, cfg.deg);
        }

        std::ostringstream cacheOss;
        {
            auto* old = std::cout.rdbuf(cacheOss.rdbuf());
            solver.printAllCache();
            std::cout.rdbuf(old);
        }
        const std::string cacheText = cacheOss.str();

        writeText(opts.outputPath, cacheText);

        std::string keyOutPath = opts.outputPath + ".keys";

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
        writeText(keyOutPath, keyOut.str());

        std::cout << "Wrote " << opts.outputPath << "\n";
        std::cout << "Wrote " << keyOutPath << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
