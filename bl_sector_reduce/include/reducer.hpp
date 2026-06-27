#pragma once

#include <functional>
#include <vector>

#include "approximant_basis.hpp"
#include "bl_config.hpp"
#include "contribution.hpp"
#include "master_data.hpp"
#include "series_store.hpp"

template<typename T>
class BLSectorReducer {
public:
    BLSectorReducer(const BLSectorConfig& config,
                    const SectorTree& tree,
                    const SeriesStore<T>& series,
                    const MasterData& masters);

    std::vector<SectorReduction<T>> reduceAll(
        const std::vector<ObjectLabel>& objects,
        const std::function<void(const std::vector<SectorReduction<T>>&)>& onSectorDone = {}) const;

private:
    BLSectorConfig config_;
    const SectorTree& tree_;
    const SeriesStore<T>& series_;
    const MasterData& masters_;
    ApproximantBasisSolver<T> solver_;

    int supportedDegree(int numMasters) const;
    std::vector<int> degreeSchedule(int maxDegree) const;
    SectorReduction<T> makeFreeMasterReduction(const ObjectLabel& object,
                                                const SectorId& sector) const;
    SectorReduction<T> makeReduction(const ObjectLabel& object,
                                      const SectorId& sector,
                                      const std::vector<ObjectLabel>& masterLabels,
                                      const std::vector<Polynomial1D<T>>& polynomials) const;
};

#include "../src/reducer.tpp"
