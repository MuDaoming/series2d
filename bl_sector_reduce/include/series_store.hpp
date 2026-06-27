#pragma once

#include <map>
#include <vector>

#include "label.hpp"
#include "sector_tree.hpp"

template<typename T>
class SeriesStore {
public:
    void addSeries(const SectorId& sector,
                   const ObjectLabel& object,
                   std::vector<T> coeffs);

    const std::vector<T>& getSeries(const SectorId& sector,
                                    const ObjectLabel& object) const;

    bool hasSeries(const SectorId& sector,
                   const ObjectLabel& object) const;

    int degree() const;
    std::vector<ObjectLabel> objects() const;
    std::vector<SectorId> sectors() const;

private:
    using ObjectMap = std::map<ObjectLabel, std::vector<T>, ObjectLabelLess>;
    std::map<SectorId, ObjectMap, SectorIdLess> data_;
    int degree_ = -1;
};

#include "../src/series_store.tpp"

