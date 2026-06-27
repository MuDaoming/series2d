#pragma once

#include <ginac/ginac.h>

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct SectorId {
    std::vector<int> bits;

    bool operator==(const SectorId& other) const { return bits == other.bits; }
    bool operator<(const SectorId& other) const { return bits < other.bits; }
};

inline int sectorIndex(const SectorId& sector) {
    int result = 0;
    for (int bit : sector.bits) {
        if (bit != 0 && bit != 1) {
            throw std::invalid_argument("Sector bits must be 0 or 1");
        }
        result = (result << 1) | bit;
    }
    return result;
}

inline std::string sectorToString(const SectorId& sector) {
    std::ostringstream out;
    out << "{";
    for (size_t i = 0; i < sector.bits.size(); ++i) {
        if (i != 0) out << ",";
        out << sector.bits[i];
    }
    out << "}";
    return out.str();
}

inline std::vector<int> activePropagators(const SectorId& sector) {
    std::vector<int> result;
    for (size_t i = 0; i < sector.bits.size(); ++i) {
        if (sector.bits[i] == 1) result.push_back(static_cast<int>(i));
    }
    return result;
}

struct CanonicalWitness {
    std::vector<int> branchPermutation;
    std::vector<int> orderedProps;
    std::string key;
};

struct SectorCanonicalForm {
    SectorId sector;
    std::string key;
    std::vector<CanonicalWitness> witnesses;
};

struct SectorMapping {
    SectorId source;
    SectorId target;
    std::vector<int> sourceToTarget;
    std::vector<int> branchPermutation;
    GiNaC::ex transformedX;
    GiNaC::ex transformedY;
    bool verified = false;
};

struct SymmetryOrbit {
    SectorId representative;
    std::vector<SectorId> members;
    std::vector<SectorMapping> mappingsToRepresentative;
    std::vector<SectorMapping> automorphisms;
};

