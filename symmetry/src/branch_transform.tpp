#pragma once

#include <algorithm>
#include <stdexcept>

inline BranchTransformBuilder::BranchTransformBuilder(
    const GiNaC::symbol& X,
    const GiNaC::symbol& Y,
    const GiNaC::ex& shiftA,
    const GiNaC::ex& shiftB)
    : X_(X), Y_(Y), shiftA_(shiftA), shiftB_(shiftB) {}

inline std::pair<GiNaC::ex, GiNaC::ex>
BranchTransformBuilder::inducedTransform(const std::vector<int>& tau) const {
    if (tau.size() != 3) {
        throw std::invalid_argument("Branch transform currently requires B=3");
    }
    std::vector<int> sorted = tau;
    std::sort(sorted.begin(), sorted.end());
    if (sorted != std::vector<int>({0, 1, 2})) {
        throw std::invalid_argument("tau is not a permutation of {0,1,2}");
    }

    const GiNaC::ex x1 = X_ + shiftA_;
    const GiNaC::ex oneMinusX1 = 1 - x1;
    const GiNaC::ex x2 = oneMinusX1 * (Y_ + shiftB_);
    const GiNaC::ex x3 = oneMinusX1 * (1 - Y_ - shiftB_);
    const std::vector<GiNaC::ex> branches = {x1, x2, x3};

    const GiNaC::ex newX = GiNaC::normal(branches[tau[0]] - shiftA_);
    const GiNaC::ex newY =
        GiNaC::normal(branches[tau[1]] / (1 - branches[tau[0]]) - shiftB_);
    return {newX, newY};
}

inline std::vector<std::vector<int>>
BranchTransformBuilder::allBranchPermutations() const {
    std::vector<std::vector<int>> result;
    std::vector<int> tau = {0, 1, 2};
    do {
        result.push_back(tau);
    } while (std::next_permutation(tau.begin(), tau.end()));
    return result;
}

