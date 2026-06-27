#pragma once

#include <stdexcept>

template<typename T>
void SeriesStore<T>::addSeries(const SectorId& sector,
                               const ObjectLabel& object,
                               std::vector<T> coeffs) {
    if (coeffs.empty()) throw std::runtime_error("Cannot add empty series");
    if (degree_ < 0) {
        degree_ = static_cast<int>(coeffs.size()) - 1;
    } else if (static_cast<int>(coeffs.size()) != degree_ + 1) {
        throw std::runtime_error("Series degree mismatch for " + objectLabelToString(object));
    }
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
    return degree_;
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

