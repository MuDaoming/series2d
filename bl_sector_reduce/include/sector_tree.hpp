#pragma once

#include <string>
#include <vector>

struct SectorId {
    std::vector<int> bits;
};

bool equalSectorId(const SectorId& lhs, const SectorId& rhs);
bool sectorContains(const SectorId& lhs, const SectorId& rhs);
int sectorPopcount(const SectorId& sector);
std::string sectorIdToString(const SectorId& sector);
SectorId parseSectorId(const std::string& text, int expectedNuSize);

struct SectorIdLess {
    bool operator()(const SectorId& lhs, const SectorId& rhs) const;
};

class SectorTree {
public:
    explicit SectorTree(const std::vector<SectorId>& sectors);

    int rootIndex() const;
    const std::vector<int>& processingOrder() const;
    int parentOf(int sectorIndex) const;
    std::vector<int> ancestorsOf(int sectorIndex) const;
    int indexOf(const SectorId& sector) const;
    const SectorId& sectorAt(int idx) const;
    int size() const;

private:
    std::vector<SectorId> sectors_;
    std::vector<int> parents_;
    std::vector<int> order_;
    int rootIndex_ = -1;

    void buildOrder();
};

#include "../src/sector_tree.tpp"
