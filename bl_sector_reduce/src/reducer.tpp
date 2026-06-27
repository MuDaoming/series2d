#pragma once

#include <sstream>
#include <stdexcept>

template<typename T>
BLSectorReducer<T>::BLSectorReducer(const BLSectorConfig& config,
                                    const SectorTree& tree,
                                    const SeriesStore<T>& series,
                                    const MasterData& masters)
    : config_(config), tree_(tree), series_(series), masters_(masters) {}

template<typename T>
std::vector<SectorReduction<T>> BLSectorReducer<T>::reduceAll(
    const std::vector<ObjectLabel>& objects,
    const std::function<void(const std::vector<SectorReduction<T>>&)>& onSectorDone) const {
    std::vector<SectorReduction<T>> reductions;
    ContributionBuilder<T> contrib(series_, tree_);

    for (const auto& object : objects) {
        for (int sectorIndex : tree_.processingOrder()) {
            const SectorId& sector = tree_.sectorAt(sectorIndex);
            if (!series_.hasSeries(sector, object)) {
                if (onSectorDone) onSectorDone(reductions);
                continue;
            }
            auto current = contrib.buildContribution(object, sectorIndex, reductions);
            if (isZeroSeries(current)) {
                SectorReduction<T> zero;
                zero.sector = sector;
                zero.object = object;
                zero.isZero = true;
                reductions.push_back(std::move(zero));
                if (onSectorDone) onSectorDone(reductions);
                continue;
            }

            if (masters_.isMaster(sector, object)) {
                reductions.push_back(makeFreeMasterReduction(object, sector));
                if (onSectorDone) onSectorDone(reductions);
                continue;
            }

            const auto& masterLabels = masters_.mastersFor(sector);
            const int r = static_cast<int>(masterLabels.size());
            const int mSupported = supportedDegree(r);
            if (mSupported < 0) {
                throw std::runtime_error("Insufficient series degree for object " +
                                         objectLabelToString(object) + " sector " +
                                         sectorIdToString(sector));
            }
            const int mMax = std::min(config_.maxDegree, mSupported);
            std::vector<std::vector<T>> masterSeries;
            masterSeries.reserve(masterLabels.size());
            for (const auto& m : masterLabels) {
                masterSeries.push_back(series_.getSeries(sector, m));
            }

            bool found = false;
            for (int m : degreeSchedule(mMax)) {
                const int workOrder =
                    (r + 1) * (m + 1) + config_.safetyOrder + config_.certOrder;
                ApproximantRequest<T> req;
                req.target = current;
                req.masters = masterSeries;
                req.maxDegree = m;
                req.workOrder = workOrder;
                const auto result = solver_.solve(req);
                if (!result.success) continue;
                reductions.push_back(makeReduction(object, sector, masterLabels, result.polynomials));
                found = true;
                break;
            }
            if (!found) {
                std::ostringstream oss;
                oss << "Unresolved sector contribution: object=" << objectLabelToString(object)
                    << " sector=" << sectorIdToString(sector)
                    << " masters=" << r
                    << " D=" << config_.degreeD
                    << " m_max=" << mMax
                    << " K_safety=" << config_.safetyOrder
                    << " K_cert=" << config_.certOrder;
                SectorReduction<T> failed;
                failed.sector = sector;
                failed.object = object;
                failed.failed = true;
                failed.failureReason = oss.str();
                const size_t prefixSize = std::min<size_t>(10, current.size());
                failed.residualPrefix.assign(current.begin(), current.begin() + prefixSize);
                reductions.push_back(std::move(failed));
                if (onSectorDone) onSectorDone(reductions);
                throw std::runtime_error(oss.str());
            }
            if (onSectorDone) onSectorDone(reductions);
        }
    }
    return reductions;
}

template<typename T>
int BLSectorReducer<T>::supportedDegree(int numMasters) const {
    const int denom = numMasters + 1;
    const int available = config_.degreeD + 1 - config_.safetyOrder - config_.certOrder;
    if (available < denom) return -1;
    return available / denom - 1;
}

template<typename T>
std::vector<int> BLSectorReducer<T>::degreeSchedule(int maxDegree) const {
    std::vector<int> out;
    if (maxDegree < 0) return out;
    out.push_back(0);
    int m = 1;
    while (m <= maxDegree) {
        out.push_back(m);
        if (m > maxDegree / 2) break;
        m *= 2;
    }
    if (out.back() != maxDegree) out.push_back(maxDegree);
    return out;
}

template<typename T>
SectorReduction<T> BLSectorReducer<T>::makeFreeMasterReduction(
    const ObjectLabel& object,
    const SectorId& sector) const {
    SectorReduction<T> red;
    red.sector = sector;
    red.object = object;
    red.denominator = Polynomial1D<T>(std::vector<T>{T(1)});
    red.isFreeMaster = true;
    ReductionTerm<T> term;
    term.master = object;
    term.numerator = Polynomial1D<T>(std::vector<T>{T(1)});
    red.terms.push_back(std::move(term));
    return red;
}

template<typename T>
SectorReduction<T> BLSectorReducer<T>::makeReduction(
    const ObjectLabel& object,
    const SectorId& sector,
    const std::vector<ObjectLabel>& masterLabels,
    const std::vector<Polynomial1D<T>>& polynomials) const {
    if (polynomials.empty() || polynomials[0].isZero()) {
        throw std::runtime_error("Cannot build reduction with zero P0");
    }
    if (polynomials.size() != masterLabels.size() + 1) {
        throw std::runtime_error("Approximant polynomial count mismatch");
    }
    SectorReduction<T> red;
    red.sector = sector;
    red.object = object;
    red.denominator = polynomials[0];
    for (size_t i = 0; i < masterLabels.size(); ++i) {
        if (polynomials[i + 1].isZero()) continue;
        ReductionTerm<T> term;
        term.master = masterLabels[i];
        term.numerator = polynomials[i + 1];
        red.terms.push_back(std::move(term));
    }
    return red;
}
