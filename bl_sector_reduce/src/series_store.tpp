#pragma once

#include <algorithm>
#include <stdexcept>

template<typename T>
void SeriesStore<T>::addSeries(const SectorId& sector,
                               const ObjectLabel& object,
                               std::vector<T> coeffs) {
    if (coeffs.empty()) throw std::runtime_error("Cannot add empty series");
    data_[sector][object] = std::move(coeffs);
}

template<typename T>
const std::vector<T>& SeriesStore<T>::getSeries(const SectorId& sector,
                                                const ObjectLabel& object) const {
    auto sit = data_.find(sector);
    if (sit == data_.end()) {
        throw std::runtime_error("Missing sector series: " + sectorIdToString(sector));
    }
    auto oit = sit->second.find(object);
    if (oit == sit->second.end()) {
        throw std::runtime_error("Missing object series: " + objectLabelToString(object) +
                                 " at sector " + sectorIdToString(sector));
    }
    return oit->second;
}

template<typename T>
bool SeriesStore<T>::hasSeries(const SectorId& sector,
                               const ObjectLabel& object) const {
    auto sit = data_.find(sector);
    if (sit == data_.end()) return false;
    return sit->second.find(object) != sit->second.end();
}

template<typename T>
int SeriesStore<T>::degree() const {
    int out = -1;
    for (const auto& sec : data_) {
        for (const auto& obj : sec.second) {
            out = std::max(out, static_cast<int>(obj.second.size()) - 1);
        }
    }
    return out;
}

template<typename T>
int SeriesStore<T>::degree(const SectorId& sector,
                           const ObjectLabel& object) const {
    const auto& s = getSeries(sector, object);
    return static_cast<int>(s.size()) - 1;
}

template<typename T>
std::vector<ObjectLabel> SeriesStore<T>::objects() const {
    std::map<ObjectLabel, bool, ObjectLabelLess> seen;
    for (const auto& sec : data_) {
        for (const auto& obj : sec.second) seen[obj.first] = true;
    }
    std::vector<ObjectLabel> out;
    for (const auto& x : seen) out.push_back(x.first);
    return out;
}

template<typename T>
std::vector<SectorId> SeriesStore<T>::sectors() const {
    std::vector<SectorId> out;
    for (const auto& x : data_) out.push_back(x.first);
    return out;
}
