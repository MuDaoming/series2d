#pragma once

#include <sstream>

inline std::string permutationToString(const std::vector<int>& permutation,
                                       bool oneBased = true) {
    std::ostringstream out;
    out << "{";
    bool first = true;
    for (size_t i = 0; i < permutation.size(); ++i) {
        if (permutation[i] < 0) continue;
        if (!first) out << ",";
        const int offset = oneBased ? 1 : 0;
        out << i + offset << "->" << permutation[i] + offset;
        first = false;
    }
    out << "}";
    return out.str();
}

inline void writeOneMapping(std::ostream& out, const SectorMapping& mapping) {
    out << "source=" << sectorToString(mapping.source)
        << " target=" << sectorToString(mapping.target)
        << " sigma=" << permutationToString(mapping.sourceToTarget)
        << " tau=" << permutationToString(mapping.branchPermutation)
        << " Xmap=" << mapping.transformedX
        << " Ymap=" << mapping.transformedY
        << " verified=" << (mapping.verified ? "true" : "false")
        << "\n";
}

inline void writeSymmetryReport(
    std::ostream& out,
    const std::vector<SymmetryOrbit>& orbits) {
    out << "[orbits]\n";
    for (const SymmetryOrbit& orbit : orbits) {
        out << "representative=" << sectorToString(orbit.representative)
            << " members={";
        for (size_t i = 0; i < orbit.members.size(); ++i) {
            if (i != 0) out << ",";
            out << sectorToString(orbit.members[i]);
        }
        out << "}\n";
    }

    out << "\n[sector_mappings]\n";
    for (const SymmetryOrbit& orbit : orbits) {
        for (const SectorMapping& mapping : orbit.mappingsToRepresentative) {
            writeOneMapping(out, mapping);
        }
    }

    out << "\n[automorphisms]\n";
    for (const SymmetryOrbit& orbit : orbits) {
        for (const SectorMapping& mapping : orbit.automorphisms) {
            writeOneMapping(out, mapping);
        }
    }
}

