#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>

#include "../../../expand/include/ff_type.hpp"
#include "../../../expand/include/io.hpp"
#include "../../include/symmetry_finder.hpp"
#include "../../include/symmetry_formatter.hpp"

namespace {

GiNaC::numeric reconstructRational(mp_limb_t residue, mp_limb_t modulus) {
    using Signed = long long;
    using Wide = __int128_t;

    const Signed m = static_cast<Signed>(modulus);
    Signed a = static_cast<Signed>(residue % modulus);
    if (a > m / 2) a -= m;

    const Signed bound =
        static_cast<Signed>(std::floor(std::sqrt(static_cast<long double>(m) / 2.0L)));

    Signed r0 = m;
    Signed r1 = a;
    Signed t0 = 0;
    Signed t1 = 1;
    while (std::llabs(r1) > bound) {
        const Signed q = r0 / r1;
        const Signed nextR = r0 - q * r1;
        const Signed nextT = t0 - q * t1;
        r0 = r1;
        r1 = nextR;
        t0 = t1;
        t1 = nextT;
    }

    Signed numerator = r1;
    Signed denominator = t1;
    if (denominator < 0) {
        numerator = -numerator;
        denominator = -denominator;
    }
    if (denominator == 0 || std::llabs(numerator) > bound ||
        denominator > bound || std::gcd(std::llabs(numerator), denominator) != 1) {
        throw std::runtime_error("Finite-field value has no small rational reconstruction");
    }

    Wide check = static_cast<Wide>(numerator) -
                 static_cast<Wide>(residue) * denominator;
    check %= static_cast<Wide>(modulus);
    if (check < 0) check += modulus;
    if (check != 0) {
        throw std::runtime_error("Rational reconstruction verification failed");
    }
    return GiNaC::numeric(numerator, denominator);
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
        fs::path sPath = fs::absolute(argv[1]);
        fs::path configPath = fs::absolute(argv[2]);
        fs::path outputPath = fs::absolute(argv[3]);
        if (outputPath.has_parent_path()) {
            fs::create_directories(outputPath.parent_path());
        }

        const InputConfig config = parseConfigFile(configPath.string());
        GiNaC::symbol X("X"), Y("Y");
        const auto topS = parseMatrixFile(sPath.string(), X, Y);
        const GiNaC::numeric shiftA = reconstructRational(config.a, config.p);
        const GiNaC::numeric shiftB = reconstructRational(config.b, config.p);

        SymmetryFinder finder(
            topS, config.N, config.B, X, Y, shiftA, shiftB);
        const auto orbits = finder.findOrbits();

        std::ofstream out(outputPath);
        if (!out.is_open()) {
            throw std::runtime_error("Cannot open output file");
        }
        out << "# shiftA=" << shiftA << "\n";
        out << "# shiftB=" << shiftB << "\n";
        writeSymmetryReport(out, orbits);

        size_t mappings = 0;
        for (const auto& orbit : orbits) {
            mappings += orbit.mappingsToRepresentative.size();
        }
        std::cout << "Found " << orbits.size() << " sector orbits and "
                  << mappings << " nontrivial sector mappings.\n";
        std::cout << "Wrote: " << outputPath << "\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << "\n";
        return 1;
    }
}

