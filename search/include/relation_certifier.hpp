#pragma once

#include <stdexcept>
#include <vector>

#include "relation_matrix_builder.hpp"
#include "relation_types.hpp"

struct DegreeWindow {
    int trainDeg = 0;
    int checkStart = 0;
    int checkEnd = -1;
};

inline DegreeWindow makeDegreeWindow(int effectiveDeg, int ncheck) {
    if (effectiveDeg < 0) {
        throw std::runtime_error("effectiveDeg must be >= 0");
    }
    if (ncheck < 0) {
        throw std::runtime_error("ncheck must be >= 0");
    }
    if (ncheck > effectiveDeg) {
        throw std::runtime_error("ncheck must satisfy ncheck <= effectiveDeg");
    }
    DegreeWindow window;
    window.trainDeg = effectiveDeg - ncheck;
    window.checkStart = window.trainDeg + 1;
    window.checkEnd = effectiveDeg;
    return window;
}

template<typename T>
std::vector<SeriesSample<T>> truncateSamplesForDegree(
    const std::vector<SeriesSample<T>>& samples,
    int degree) {
    std::vector<SeriesSample<T>> out = samples;
    for (auto& sample : out) {
        if (static_cast<int>(sample.coeffs.size()) < degree + 1) {
            throw std::runtime_error("series coeff count smaller than requested degree+1");
        }
        sample.coeffs.resize(degree + 1);
    }
    return out;
}

template<typename T>
bool verifyRelationOnWindow(
    const std::vector<IntegralLabel>& labels,
    const std::vector<std::vector<T>>& polynomials,
    const std::vector<SeriesSample<T>>& samples,
    int numBC,
    int startDeg,
    int endDeg) {
    if (startDeg > endDeg) {
        return true;
    }
    if (labels.size() != polynomials.size()) {
        throw std::runtime_error("labels/polynomials size mismatch");
    }
    const int numLabels = static_cast<int>(labels.size());
    if (numBC <= 0 || static_cast<int>(samples.size()) != numBC * numLabels) {
        throw std::runtime_error("samples size mismatch in relation verification");
    }

    const T zero(0ULL);
    for (int bc = 0; bc < numBC; ++bc) {
        for (int n = startDeg; n <= endDeg; ++n) {
            T sum(0ULL);
            for (int i = 0; i < numLabels; ++i) {
                const auto& poly = polynomials[i];
                const auto& coeffs = samples[bc * numLabels + i].coeffs;
                for (int k = 0; k < static_cast<int>(poly.size()); ++k) {
                    if (n >= k) {
                        if (n - k >= static_cast<int>(coeffs.size())) {
                            throw std::runtime_error("series coeff count too short for check window");
                        }
                        sum += poly[k] * coeffs[n - k];
                    }
                }
            }
            if (sum != zero) {
                return false;
            }
        }
    }
    return true;
}

template<typename T>
std::vector<std::vector<T>> buildCheckRows(
    const SearchInput<T>& fullInput,
    const std::vector<RelationVariable>& variables,
    int checkStart,
    int checkEnd) {
    std::vector<std::vector<T>> rows;
    if (checkStart > checkEnd) {
        return rows;
    }
    RelationMatrixBuilder<T> fullBuilder(fullInput);
    return fullBuilder.buildRowsForDegreeWindow(variables, checkStart, checkEnd);
}

template<typename T>
bool checkNullspaceShrink(
    const std::vector<std::vector<T>>& trainRREF,
    const std::vector<int>& pivotColumns,
    const std::vector<std::vector<T>>& checkRows) {
    const T zero(0ULL);
    if (trainRREF.empty() || pivotColumns.empty()) {
        for (const auto& row : checkRows) {
            for (const auto& x : row) {
                if (x != zero) return true;
            }
        }
        return false;
    }

    for (const auto& rawRow : checkRows) {
        std::vector<T> row = rawRow;
        for (int r = 0; r < static_cast<int>(trainRREF.size()); ++r) {
            const int pc = pivotColumns[r];
            if (pc < 0 || pc >= static_cast<int>(row.size())) continue;
            const T factor = row[pc];
            if (factor == zero) continue;
            for (int c = pc; c < static_cast<int>(row.size()); ++c) {
                row[c] = row[c] - factor * trainRREF[r][c];
            }
        }
        for (const auto& x : row) {
            if (x != zero) {
                return true;
            }
        }
    }
    return false;
}
