#pragma once

#include <map>
#include <vector>

#include "label.hpp"
#include "sector_tree.hpp"

class MasterData {
public:
    void setMasters(const SectorId& sector, std::vector<ObjectLabel> masters);
    const std::vector<ObjectLabel>& mastersFor(const SectorId& sector) const;
    bool isMaster(const SectorId& sector, const ObjectLabel& object) const;

private:
    std::map<SectorId, std::vector<ObjectLabel>, SectorIdLess> masters_;
};

#include "../src/master_data.tpp"

