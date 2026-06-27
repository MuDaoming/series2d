#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "approximant_basis.hpp"
#include "ff_type.hpp"
#include "formatter.hpp"
#include "io.hpp"

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0]
                  << " <config> <sector_series_list> <sector> <object> <output>\n";
        return 1;
    }

    try {
        const BLSectorConfig config = parseBLSectorConfig(argv[1]);
        FlintMod::set_modulus(config.prime);
        const SectorId wantedSector = parseSectorId(argv[3], config.nuSize);
        const ObjectLabel object = parseObjectLabel(argv[4], config.nuSize);

        auto entries = parseSectorSeriesList(argv[2], config.nuSize);
        entries.erase(
            std::remove_if(entries.begin(), entries.end(), [&](const SectorSeriesEntry& entry) {
                return !equalSectorId(entry.sector, wantedSector);
            }),
            entries.end());
        if (entries.size() != 1) {
            throw std::runtime_error("Requested sector must match exactly one input entry");
        }

        SeriesStore<FlintMod> series;
        MasterData masters;
        loadSectorData(entries, config.degreeD, config.nuSize, series, masters);
        const auto& masterLabels = masters.mastersFor(wantedSector);
        const int r = static_cast<int>(masterLabels.size());
        const int workOrder =
            (r + 1) * (config.maxDegree + 1) + config.safetyOrder + config.certOrder;
        if (workOrder > config.degreeD + 1) {
            throw std::runtime_error("Requested work order exceeds available series length");
        }

        ApproximantRequest<FlintMod> request;
        request.target = series.getSeries(wantedSector, object);
        request.maxDegree = config.maxDegree;
        request.workOrder = workOrder;
        for (const auto& master : masterLabels) {
            request.masters.push_back(series.getSeries(wantedSector, master));
        }

        ApproximantBasisSolver<FlintMod> solver;
        const auto basis = solver.basis(request);
        std::ofstream out(argv[5]);
        if (!out.is_open()) throw std::runtime_error("Cannot open output file");

        out << "# sector = " << sectorIdToString(wantedSector) << "\n";
        out << "# object = " << objectLabelToString(object) << "\n";
        out << "# masters = " << r << "\n";
        out << "# rows = " << basis.size() << "\n";
        out << "# m = " << config.maxDegree << "\n";
        out << "# K = " << workOrder << "\n\n";
        for (int i = 0; i < static_cast<int>(masterLabels.size()); ++i) {
            out << "# M" << i + 1 << " = " << objectLabelToString(masterLabels[i]) << "\n";
        }
        out << "\n";

        for (int i = 0; i < static_cast<int>(basis.size()); ++i) {
            int maxDegree = -1;
            for (const auto& poly : basis[i]) maxDegree = std::max(maxDegree, poly.degree());
            out << "[row " << i << "]\n";
            out << "max_degree=" << maxDegree << "\n";
            out << "P0_nonzero=" << (basis[i][0].isZero() ? 0 : 1) << "\n";
            for (int j = 0; j < static_cast<int>(basis[i].size()); ++j) {
                out << "P" << j << "=" << polynomialToString(basis[i][j]) << "\n";
            }
            out << "\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "[dump_approximant_basis] Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
