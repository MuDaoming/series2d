#pragma once

#include <ginac/ginac.h>

#include <map>
#include <vector>

#include "sector_canonicalizer.hpp"
#include "symmetry_types.hpp"

class SymmetryFinder {
public:
    SymmetryFinder(const std::vector<std::vector<GiNaC::ex>>& topS,
                   int numProps,
                   int numBranches,
                   const GiNaC::symbol& X,
                   const GiNaC::symbol& Y,
                   const GiNaC::ex& shiftA,
                   const GiNaC::ex& shiftB);

    std::vector<SectorId> enumerateValidSectors() const;
    std::vector<SymmetryOrbit> findOrbits() const;
    bool verify(SectorMapping& mapping) const;

    const std::vector<int>& branchOfProp() const { return branchOfProp_; }

private:
    int numProps_;
    int numBranches_;
    const GiNaC::symbol& X_;
    const GiNaC::symbol& Y_;
    GiNaC::ex shiftA_;
    GiNaC::ex shiftB_;
    std::vector<int> branchOfProp_;
    std::vector<std::vector<GiNaC::ex>> R_;
    SectorCanonicalizer canonicalizer_;
    BranchTransformBuilder transforms_;

    static SectorId sectorFromIndex(int index, int numProps);
    SectorMapping mappingFromWitnesses(
        const SectorCanonicalForm& source,
        const CanonicalWitness& sourceWitness,
        const SectorCanonicalForm& target,
        const CanonicalWitness& targetWitness) const;
};

#include "../src/symmetry_finder.tpp"

