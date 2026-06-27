#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "ff_type.hpp"
#include "formatter.hpp"
#include "io.hpp"
#include "reducer.hpp"

int main(int argc, char** argv) {
    if (argc != 4 && argc != 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <config_path> <sector_series_list_path> <output_path> [object_list_path]\n";
        return 1;
    }

    const std::string configPath = argv[1];
    const std::string sectorListPath = argv[2];
    const std::string outputPath = argv[3];
    const bool hasObjectList = argc == 5;
    const std::string objectListPath = hasObjectList ? argv[4] : "";

    try {
        const BLSectorConfig config = parseBLSectorConfig(configPath);
        FlintMod::set_modulus(config.prime);

        const auto entries = parseSectorSeriesList(sectorListPath, config.nuSize);
        SeriesStore<FlintMod> series;
        MasterData masters;
        loadSectorData(entries, config.degreeD, config.nuSize, series, masters);

        const auto sectors = series.sectors();
        SectorTree tree(sectors);

        std::vector<ObjectLabel> objects =
            hasObjectList ? parseObjectList(objectListPath, config.nuSize) : series.objects();

        BLSectorReducer<FlintMod> reducer(config, tree, series, masters);
        auto writeSnapshot = [&](const std::vector<SectorReduction<FlintMod>>& reductions) {
            std::ofstream out(outputPath);
            if (!out.is_open()) {
                throw std::runtime_error("Cannot open output file: " + outputPath);
            }
            writeReductions(out, config, tree, masters, reductions);
        };
        writeSnapshot({});
        const auto reductions = reducer.reduceAll(objects, writeSnapshot);
        writeSnapshot(reductions);
    } catch (const std::exception& e) {
        std::cerr << "[bl_sector_reducer] Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
