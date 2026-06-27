#pragma once

#include <stdexcept>

template<typename T>
ApproximantResult<T> ApproximantBasisSolver<T>::solve(
    const ApproximantRequest<T>& request) const {
    const auto resultBasis = basis(request);
    const int r = static_cast<int>(request.masters.size());
    const int n = r + 1;
    std::vector<std::vector<T>> f(n);
    f[0] = request.target;
    for (int j = 0; j < r; ++j) {
        f[j + 1].resize(request.masters[j].size());
        for (size_t k = 0; k < request.masters[j].size(); ++k) {
            f[j + 1][k] = -request.masters[j][k];
        }
    }

    for (const auto& row : resultBasis) {
        if (row[0].isZero()) continue;
        if (!satisfiesBound(row, request.maxDegree)) continue;
        if (!verify(f, row, request.workOrder)) continue;
        ApproximantResult<T> result;
        result.success = true;
        result.polynomials = row;
        return result;
    }

    return {};
}

template<typename T>
std::vector<typename ApproximantBasisSolver<T>::PolyVec>
ApproximantBasisSolver<T>::basis(const ApproximantRequest<T>& request) const {
    if (request.maxDegree < 0) throw std::runtime_error("maxDegree must be nonnegative");
    if (request.workOrder <= 0) throw std::runtime_error("workOrder must be positive");
    const int r = static_cast<int>(request.masters.size());
    const int n = r + 1;

    std::vector<std::vector<T>> f(n);
    f[0] = request.target;
    for (int j = 0; j < r; ++j) {
        f[j + 1].resize(request.masters[j].size());
        for (size_t k = 0; k < request.masters[j].size(); ++k) {
            f[j + 1][k] = -request.masters[j][k];
        }
    }

    std::vector<PolyVec> basis(n, PolyVec(n));
    for (int i = 0; i < n; ++i) {
        basis[i][i].setCoeff(0, T(1));
    }

    for (int k = 0; k < request.workOrder; ++k) {
        std::vector<T> disc(n, T(0));
        int pivot = -1;
        int pivotDeg = 0;
        for (int i = 0; i < n; ++i) {
            disc[i] = discrepancy(f, basis[i], k);
            if (disc[i] == T(0)) continue;
            const int deg = rowDegree(basis[i]);
            if (pivot < 0 || deg < pivotDeg) {
                pivot = i;
                pivotDeg = deg;
            }
        }
        if (pivot < 0) continue;

        const PolyVec pivotRow = basis[pivot];
        const T pivotDisc = disc[pivot];
        for (int i = 0; i < n; ++i) {
            if (i == pivot || disc[i] == T(0)) continue;
            const T factor = disc[i] / pivotDisc;
            for (int j = 0; j < n; ++j) {
                basis[i][j] = basis[i][j] - pivotRow[j] * factor;
            }
        }
        for (int j = 0; j < n; ++j) {
            basis[pivot][j].shiftByX();
        }
    }

    return basis;
}

template<typename T>
T ApproximantBasisSolver<T>::discrepancy(const std::vector<std::vector<T>>& f,
                                         const PolyVec& row,
                                         int order) const {
    T out(0);
    for (int j = 0; j < static_cast<int>(f.size()); ++j) {
        for (int d = 0; d <= order; ++d) {
            const T pc = row[j].coeff(d);
            if (pc == T(0)) continue;
            const int idx = order - d;
            if (idx < static_cast<int>(f[j].size())) {
                out += pc * f[j][idx];
            }
        }
    }
    return out;
}

template<typename T>
int ApproximantBasisSolver<T>::rowDegree(const PolyVec& row) const {
    int d = -1;
    for (const auto& p : row) d = std::max(d, p.degree());
    return d;
}

template<typename T>
bool ApproximantBasisSolver<T>::satisfiesBound(const PolyVec& row, int maxDegree) const {
    for (const auto& p : row) {
        if (p.degree() > maxDegree) return false;
    }
    return true;
}

template<typename T>
bool ApproximantBasisSolver<T>::verify(const std::vector<std::vector<T>>& f,
                                       const PolyVec& row,
                                       int workOrder) const {
    for (int k = 0; k < workOrder; ++k) {
        if (discrepancy(f, row, k) != T(0)) return false;
    }
    return true;
}
