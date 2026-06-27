#pragma once

#include <string>
#include <utility>
#include <vector>

struct SectorMapEntry {
    std::vector<int> source;
    std::vector<int> target;
    std::vector<int> sourceToTarget;
};

class SectorMap {
public:
    SectorMap() = default;
    explicit SectorMap(int numProps) : numProps_(numProps) {}

    static SectorMap fromFile(const std::string& path, int numProps);

    bool empty() const { return entries_.empty(); }
    int numProps() const { return numProps_; }
    const std::vector<SectorMapEntry>& entries() const { return entries_; }

    bool hasSource(const std::vector<int>& sector) const;
    std::vector<int> canonicalizeSector(const std::vector<int>& sector) const;
    std::vector<int> canonicalizeNu(const std::vector<int>& nu) const;

private:
    int numProps_ = 0;
    std::vector<SectorMapEntry> entries_;

    const SectorMapEntry* findEntry(const std::vector<int>& sector) const;
    void addEntry(const std::vector<int>& source,
                  const std::vector<int>& target,
                  const std::vector<std::pair<int, int>>& explicitMap);
};

#include "../src/sector_map.tpp"
