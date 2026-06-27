#pragma once

#include <iosfwd>
#include <vector>

#include "bl_config.hpp"
#include "contribution.hpp"
#include "master_data.hpp"
#include "sector_tree.hpp"

template<typename T>
void writeReductions(std::ostream& out,
                     const BLSectorConfig& config,
                     const SectorTree& tree,
                     const MasterData& masters,
                     const std::vector<SectorReduction<T>>& reductions);

template<typename T>
std::string polynomialToString(const Polynomial1D<T>& poly);

#include "../src/formatter.tpp"

