#pragma once

#include <ginac/ginac.h>

#include <string>
#include <vector>

#include "branch_transform.hpp"
#include "symmetry_types.hpp"

class SectorCanonicalizer {
public:
    SectorCanonicalizer(const std::vector<std::vector<GiNaC::ex>>& R,
                        const std::vector<int>& branchOfProp,
                        const GiNaC::symbol& X,
                        const GiNaC::symbol& Y,
                        const GiNaC::ex& shiftA,
                        const GiNaC::ex& shiftB);

    SectorCanonicalForm canonicalize(const SectorId& sector) const;

    static std::string expressionKey(const GiNaC::ex& value);

private:
    struct PartialCandidate {
        std::vector<int> orderedProps;
        std::vector<bool> used;
        std::vector<std::string> prefixTokens;
    };

    const std::vector<std::vector<GiNaC::ex>>& R_;
    const std::vector<int>& branchOfProp_;
    const GiNaC::symbol& X_;
    const GiNaC::symbol& Y_;
    BranchTransformBuilder transforms_;

    std::vector<CanonicalWitness> canonicalizeForTau(
        const SectorId& sector,
        const std::vector<int>& tau) const;

    std::vector<std::vector<GiNaC::ex>> transformedR(
        const std::pair<GiNaC::ex, GiNaC::ex>& transform) const;

    static std::string encodeTokens(const std::vector<std::string>& tokens);
};

#include "../src/sector_canonicalizer.tpp"
