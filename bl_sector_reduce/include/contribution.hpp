#pragma once

#include <string>
#include <vector>

#include "label.hpp"
#include "polynomial_1d.hpp"
#include "sector_tree.hpp"
#include "series_store.hpp"

template<typename T>
struct ReductionTerm {
    ObjectLabel master;
    Polynomial1D<T> numerator;
};

template<typename T>
struct SectorReduction {
    SectorId sector;
    ObjectLabel object;
    Polynomial1D<T> denominator;
    std::vector<ReductionTerm<T>> terms;
    bool isFreeMaster = false;
    bool isZero = false;
    bool failed = false;
    std::string failureReason;
    std::vector<T> residualPrefix;
};

template<typename T>
struct SeenReductionTerm {
    ObjectLabel master;
    Polynomial1D<T> numerator;
    Polynomial1D<T> denominator;
};

template<typename T>
class ContributionBuilder {
public:
    ContributionBuilder(const SeriesStore<T>& series,
                        const SectorTree& tree);

    std::vector<T> buildContribution(
        const ObjectLabel& object,
        int sectorIndex,
        const std::vector<SectorReduction<T>>& knownReductions) const;

private:
    const SeriesStore<T>& series_;
    const SectorTree& tree_;

    std::vector<T> evaluateReductionAtSector(
        const SectorReduction<T>& reduction,
        const SectorId& sector,
        std::vector<SeenReductionTerm<T>>& seenTerms) const;
};

template<typename T>
bool isZeroSeries(const std::vector<T>& series);

#include "../src/contribution.tpp"
