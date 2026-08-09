#pragma once

#include <string>
#include <vector>

#include "bl_config.hpp"
#include "master_data.hpp"
#include "series_store.hpp"

struct SectorSeriesEntry {
    SectorId sector;
    std::string seriesPath;
    std::string targetPath;
    std::string masterPath;
    int degreeD = -1;
};

BLSectorConfig parseBLSectorConfig(const std::string& path);
std::vector<SectorSeriesEntry> parseSectorSeriesList(const std::string& path,
                                                     int expectedNuSize);
std::vector<ObjectLabel> parseObjectList(const std::string& path, int expectedNuSize);
std::vector<ObjectLabel> parseTargetFile(const std::string& path, int expectedNuSize);

template<typename T>
void loadSectorData(const std::vector<SectorSeriesEntry>& entries,
                    int& degreeD,
                    int expectedNuSize,
                    SeriesStore<T>& series,
                    MasterData& masters);

#include "../src/io.tpp"
