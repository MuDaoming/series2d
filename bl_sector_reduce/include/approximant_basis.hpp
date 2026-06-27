#pragma once

#include <vector>

#include "polynomial_1d.hpp"

template<typename T>
struct ApproximantRequest {
    std::vector<T> target;
    std::vector<std::vector<T>> masters;
    int maxDegree = 0;
    int workOrder = 0;
};

template<typename T>
struct ApproximantResult {
    bool success = false;
    std::vector<Polynomial1D<T>> polynomials;
};

template<typename T>
class ApproximantBasisSolver {
public:
    using PolyVec = std::vector<Polynomial1D<T>>;

    ApproximantResult<T> solve(const ApproximantRequest<T>& request) const;
    std::vector<PolyVec> basis(const ApproximantRequest<T>& request) const;

private:
    T discrepancy(const std::vector<std::vector<T>>& f,
                  const PolyVec& row,
                  int order) const;

    int rowDegree(const PolyVec& row) const;
    bool satisfiesBound(const PolyVec& row, int maxDegree) const;
    bool verify(const std::vector<std::vector<T>>& f,
                const PolyVec& row,
                int workOrder) const;
};

#include "../src/approximant_basis.tpp"
