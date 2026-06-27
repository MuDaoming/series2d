#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../include/converter.hpp"
#include "../../include/io.hpp"

#include <ginac/ginac.h>

namespace {

struct CliOptions {
    std::string sPath;
    std::string configPath;
    std::string outputPrefix;
};

std::string trim(const std::string& s) {
    const char* ws = " \t\n\r";
    size_t b = s.find_first_not_of(ws);
    if (b == std::string::npos) return "";
    size_t e = s.find_last_not_of(ws);
    return s.substr(b, e - b + 1);
}

CliOptions parseArgs(int argc, char** argv) {
    if (argc != 4) {
        throw std::runtime_error(
            "Usage: dump_cz <S_path> <config_path> <output_prefix>");
    }
    CliOptions opts;
    opts.sPath = argv[1];
    opts.configPath = argv[2];
    opts.outputPrefix = argv[3];
    return opts;
}

std::string sectorSuffix(const std::vector<int>& nu) {
    std::ostringstream oss;
    oss << "sec_";
    for (size_t i = 0; i < nu.size(); ++i) {
        if (i) oss << "_";
        oss << nu[i];
    }
    return oss.str();
}

} // namespace

int main(int argc, char** argv) {
    try {
        CliOptions opts = parseArgs(argc, argv);
        InputConfig cfg = parseConfigFile(opts.configPath);
        FlintMod::set_modulus(cfg.p);

        GiNaC::symbol X("X"), Y("Y");
        auto topS = parseMatrixFile(opts.sPath, X, Y);
        if (topS.empty() || topS.size() != topS[0].size()) {
            throw std::runtime_error("S matrix must be non-empty and square");
        }

        std::vector<int> nu = cfg.sector.empty() ? std::vector<int>(cfg.N, 1) : cfg.sector;
        if (static_cast<int>(nu.size()) != cfg.N) {
            throw std::runtime_error("config.sector size mismatch with N");
        }

        Family<GiNaC::ex, GiNaC::ex, GiNaC::ex> ginacFamily(topS, cfg.N, cfg.B, X, Y);
        if (!cfg.sectorMapPath.empty()) {
            ginacFamily.setSectorMap(SectorMap::fromFile(cfg.sectorMapPath, cfg.N));
        }
        const auto* ginacSector = ginacFamily.getSector(nu);
        if (!ginacSector) {
            throw std::runtime_error("sector was not found before convert");
        }

        auto family = convertFamily(ginacFamily, X, Y);
        const auto* sector = family.getSector(nu);
        if (!sector) {
            throw std::runtime_error("sector was not found after convert");
        }

        const std::string suffix = sectorSuffix(nu);
        const std::string prefix = opts.outputPrefix;
        const std::string textPath = prefix + ".txt";
        const std::string wlPath = prefix + "_data.wl";

        std::ofstream out(textPath);
        if (!out) {
            throw std::runtime_error("Cannot open output file: " + textPath);
        }

        out << "sector = {";
        for (int i = 0; i < cfg.N; ++i) {
            if (i) out << ",";
            out << nu[i];
        }
        out << "}\n";
        out << "\n[before_convert: Sector<GiNaC::ex, GiNaC::ex>]\n";
        out << "case = " << ginacSector->getCase() << "\n";
        out << "dimNull = " << ginacSector->getDimNull() << "\n";
        out << "z0 = " << ginacSector->getZ0() << "\n";
        out << "C = " << ginacSector->getCSum() << "\n";
        out << "denoCandZ = " << ginacSector->getDenoCandZ() << "\n";
        out << "numeC = " << ginacSector->getNumeC() << "\n";
        out << "z = {";
        for (int i = 0; i < cfg.N; ++i) {
            if (i) out << ", ";
            out << ginacSector->getZ(i);
        }
        out << "}\n";
        out << "numeZ = {";
        for (int i = 0; i < cfg.N; ++i) {
            if (i) out << ", ";
            out << ginacSector->getNumeZ(i);
        }
        out << "}\n";

        out << "\n[after_convert: Sector<Rational<FlintMod>, Polynomial<FlintMod>>]\n";
        out << "case = " << sector->getCase() << "\n";
        out << "dimNull = " << sector->getDimNull() << "\n";
        out << "z0 = " << sector->getZ0() << "\n";
        out << "C = " << sector->getCSum().toString() << "\n";
        out << "denoCandZ = " << sector->getDenoCandZ().toString() << "\n";
        out << "numeC = " << sector->getNumeC().toString() << "\n";
        out << "z = {";
        for (int i = 0; i < cfg.N; ++i) {
            if (i) out << ", ";
            out << sector->getZ(i).toString();
        }
        out << "}\n";
        out << "numeZ = {";
        for (int i = 0; i < cfg.N; ++i) {
            if (i) out << ", ";
            out << sector->getNumeZ(i).toString();
        }
        out << "}\n";

        std::ofstream wl(wlPath);
        if (!wl) {
            throw std::runtime_error("Cannot open output file: " + wlPath);
        }
        wl << "cppBeforeC = " << ginacSector->getCSum() << ";\n";
        wl << "cppBeforeZ = {";
        for (int i = 0; i < cfg.N; ++i) {
            if (i) wl << ", ";
            wl << ginacSector->getZ(i);
        }
        wl << "};\n";
        wl << "cppBeforeCandZ = {";
        for (int i = 0; i < cfg.B; ++i) {
            if (i) wl << ", ";
            wl << ginacSector->getC(i);
        }
        for (int i = 0; i < cfg.N; ++i) {
            wl << ", " << ginacSector->getZ(i);
        }
        wl << "};\n";

        wl << "cppAfterC = " << sector->getCSum().toString() << ";\n";
        wl << "cppAfterZ = {";
        for (int i = 0; i < cfg.N; ++i) {
            if (i) wl << ", ";
            wl << sector->getZ(i).toString();
        }
        wl << "};\n";
        wl << "cppAfterCandZ = {";
        for (int i = 0; i < cfg.B; ++i) {
            if (i) wl << ", ";
            wl << sector->getC(i).toString();
        }
        for (int i = 0; i < cfg.N; ++i) {
            wl << ", " << sector->getZ(i).toString();
        }
        wl << "};\n";

        std::cout << "Wrote " << textPath << "\n";
        std::cout << "Wrote " << wlPath << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
