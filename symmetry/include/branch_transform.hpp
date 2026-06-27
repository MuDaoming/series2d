#pragma once

#include <ginac/ginac.h>

#include <utility>
#include <vector>

class BranchTransformBuilder {
public:
    BranchTransformBuilder(const GiNaC::symbol& X,
                           const GiNaC::symbol& Y,
                           const GiNaC::ex& shiftA,
                           const GiNaC::ex& shiftB);

    std::pair<GiNaC::ex, GiNaC::ex> inducedTransform(
        const std::vector<int>& tau) const;

    std::vector<std::vector<int>> allBranchPermutations() const;

private:
    const GiNaC::symbol& X_;
    const GiNaC::symbol& Y_;
    GiNaC::ex shiftA_;
    GiNaC::ex shiftB_;
};

#include "../src/branch_transform.tpp"

