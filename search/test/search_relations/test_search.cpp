#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "fi_pipeline.hpp"
#include "search_pipeline.hpp"

int main() {
    try {
        const std::string base = ".";
        const std::string seriesPath = base + "/series";
        const std::string relationPath = base + "/relations";

        runFI1DSeriesPipeline(
            base + "/S",
            base + "/config",
            base + "/target",
            seriesPath);

        runRelationSearchPipeline(
            std::vector<std::string>{seriesPath},
            base + "/config",
            base + "/target",
            base + "/maxsearchdeg",
            relationPath);

        std::cout << "search test pipeline finished.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
