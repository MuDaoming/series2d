#pragma once

#include <stdexcept>

void MasterData::setMasters(const SectorId& sector, std::vector<ObjectLabel> masters) {
    masters_[sector] = std::move(masters);
}

const std::vector<ObjectLabel>& MasterData::mastersFor(const SectorId& sector) const {
    auto it = masters_.find(sector);
    if (it == masters_.end()) {
        throw std::runtime_error("Missing master list for sector " + sectorIdToString(sector));
    }
    return it->second;
}

bool MasterData::isMaster(const SectorId& sector, const ObjectLabel& object) const {
    const auto& ms = mastersFor(sector);
    for (const auto& m : ms) {
        if (equalObjectLabel(m, object)) return true;
    }
    return false;
}

