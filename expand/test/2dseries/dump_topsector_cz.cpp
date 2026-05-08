#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../include/converter.hpp"
#include "../../include/io.hpp"

#include <ginac/ginac.h>

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <work_dir>\n";
        return 1;
    }

    try {
        const std::string workDir = argv[1];
        InputConfig cfg = parseConfigFile(workDir + "/config");
        FlintMod::set_modulus(cfg.p);

        GiNaC::symbol X("X"), Y("Y");
        auto topS = parseMatrixFile(workDir + "/S", X, Y);
        if (topS.empty() || topS.size() != topS[0].size()) {
            throw std::runtime_error("S matrix must be non-empty and square");
        }

        Family<GiNaC::ex, GiNaC::ex, GiNaC::ex> ginacFamily(topS, cfg.N, cfg.B, X, Y);
        std::vector<int> topNu(cfg.N, 1);
        const auto* ginacSector = ginacFamily.getSector(topNu);
        if (!ginacSector) {
            throw std::runtime_error("topsector was not found before convert");
        }

        auto family = convertFamily(ginacFamily, X, Y);
        const auto* sector = family.getSector(topNu);
        if (!sector) {
            throw std::runtime_error("topsector was not found after convert");
        }

        std::ofstream out(workDir + "/topsector_C_z.txt");
        if (!out) {
            throw std::runtime_error("Cannot open output file");
        }

        out << "topsector = {";
        for (int i = 0; i < cfg.N; ++i) {
            if (i) out << ",";
            out << 1;
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

        std::ofstream wl(workDir + "/topsector_C_z_data.wl");
        if (!wl) {
            throw std::runtime_error("Cannot open wl output file");
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

        std::cout << "Wrote " << workDir << "/topsector_C_z.txt\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
