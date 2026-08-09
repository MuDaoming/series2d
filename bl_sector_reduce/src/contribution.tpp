#pragma once

#include <stdexcept>
#include <sstream>
#include <vector>

namespace contribution_detail {

template<typename T>
std::vector<T> trim(std::vector<T> coeffs) {
    while (coeffs.size() > 1 && coeffs.back() == T(0)) coeffs.pop_back();
    if (coeffs.empty()) coeffs.push_back(T(0));
    return coeffs;
}

template<typename T>
std::vector<T> multiplyCoeffs(const std::vector<T>& lhs,
                              const std::vector<T>& rhs) {
    if (lhs.empty() || rhs.empty()) return std::vector<T>{T(0)};
    std::vector<T> out(lhs.size() + rhs.size() - 1, T(0));
    for (size_t i = 0; i < lhs.size(); ++i) {
        for (size_t j = 0; j < rhs.size(); ++j) {
            out[i + j] += lhs[i] * rhs[j];
        }
    }
    return trim(std::move(out));
}

template<typename T>
bool equalCoeffs(const std::vector<T>& lhs, const std::vector<T>& rhs) {
    return trim(lhs) == trim(rhs);
}

template<typename T>
bool equalRationalPolynomial(const Polynomial1D<T>& lhsNumerator,
                             const Polynomial1D<T>& lhsDenominator,
                             const Polynomial1D<T>& rhsNumerator,
                             const Polynomial1D<T>& rhsDenominator) {
    const auto lhsCross = multiplyCoeffs(lhsNumerator.coeffs(),
                                         rhsDenominator.coeffs());
    const auto rhsCross = multiplyCoeffs(rhsNumerator.coeffs(),
                                         lhsDenominator.coeffs());
    return equalCoeffs(lhsCross, rhsCross);
}

template<typename T>
int duplicateStatus(const ObjectLabel& master,
                    const Polynomial1D<T>& numerator,
                    const Polynomial1D<T>& denominator,
                    const std::vector<SeenReductionTerm<T>>& seenTerms) {
    for (const auto& seen : seenTerms) {
        if (!equalObjectLabel(master, seen.master)) continue;
        if (equalRationalPolynomial(numerator, denominator,
                                    seen.numerator, seen.denominator)) {
            return 1;
        }
        return -1;
    }
    return 0;
}

}  // namespace contribution_detail

template<typename T>
ContributionBuilder<T>::ContributionBuilder(const SeriesStore<T>& series,
                                             const SectorTree& tree)
    : series_(series), tree_(tree) {}

template<typename T>
std::vector<T> ContributionBuilder<T>::buildContribution(
    const ObjectLabel& object,
    int sectorIndex,
    const std::vector<SectorReduction<T>>& knownReductions) const {
    const SectorId& sector = tree_.sectorAt(sectorIndex);
    std::vector<T> contribution = series_.getSeries(sector, object);
    const int degreeD = static_cast<int>(contribution.size()) - 1;
    const auto ancestors = tree_.ancestorsOf(sectorIndex);
    std::vector<SeenReductionTerm<T>> seenTerms;

    for (int ancestorIdx : ancestors) {
        const SectorId& ancestor = tree_.sectorAt(ancestorIdx);
        for (const auto& red : knownReductions) {
            if (!equalSectorId(red.sector, ancestor)) continue;
            if (!equalObjectLabel(red.object, object)) continue;
            if (red.isZero) continue;
            const auto eval = evaluateReductionAtSector(red, sector, degreeD, seenTerms);
            for (int i = 0; i <= degreeD; ++i) contribution[i] -= eval[i];
        }
    }
    return contribution;
}

template<typename T>
std::vector<T> ContributionBuilder<T>::evaluateReductionAtSector(
    const SectorReduction<T>& reduction,
    const SectorId& sector,
    int degreeD,
    std::vector<SeenReductionTerm<T>>& seenTerms) const {
    std::vector<T> numeratorSeries(degreeD + 1, T(0));
    for (const auto& term : reduction.terms) {
        const int duplicate = contribution_detail::duplicateStatus(
            term.master, term.numerator, reduction.denominator, seenTerms);
        if (duplicate == 1) {
            continue;
        }
        if (duplicate == -1) {
            std::ostringstream oss;
            oss << "Conflicting duplicate master contribution for "
                << objectLabelToString(term.master)
                << " while evaluating sector " << sectorIdToString(sector)
                << ". The same master appears in multiple ancestor reductions "
                << "with different rational coefficients.";
            throw std::runtime_error(oss.str());
        }
        seenTerms.push_back(SeenReductionTerm<T>{
            term.master,
            term.numerator,
            reduction.denominator
        });
        const auto& masterSeries = series_.getSeries(sector, term.master);
        if (static_cast<int>(masterSeries.size()) < degreeD + 1) {
            std::ostringstream oss;
            oss << "Master series is shorter than target contribution: master="
                << objectLabelToString(term.master)
                << " sector=" << sectorIdToString(sector)
                << " have_D=" << static_cast<int>(masterSeries.size()) - 1
                << " need_D=" << degreeD;
            throw std::runtime_error(oss.str());
        }
        const auto prod = multiplyPolySeries(term.numerator, masterSeries, degreeD);
        for (int i = 0; i <= degreeD; ++i) numeratorSeries[i] += prod[i];
    }
    return divideSeriesByPoly(numeratorSeries, reduction.denominator, degreeD);
}

template<typename T>
bool isZeroSeries(const std::vector<T>& series) {
    for (const auto& x : series) {
        if (x != T(0)) return false;
    }
    return true;
}
