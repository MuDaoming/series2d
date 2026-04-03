#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "fi_pipeline.hpp"
#include "coefficient_relation_expander.hpp"
#include "ff_type.hpp"
#include "fi_reduction_searcher.hpp"
#include "io.hpp"
#include "relation_formatter.hpp"
#include "relation_searcher.hpp"

int main() {
    try {
        const std::string inputBase = "../search_relations";
        const std::string outputBase = ".";

        runFI1DSeriesPipeline(
            inputBase + "/S",
            inputBase + "/config",
            inputBase + "/target",
            outputBase + "/series");

        SearchConfig cfg = parseSearchConfigFile(inputBase + "/config");
        FlintMod::set_modulus(cfg.p);
        const int maxSearchDeg = parseMaxSearchDegreeFile(inputBase + "/maxsearchdeg");
        auto targets = parseSearchTargetFile(inputBase + "/target", cfg.nuSize);

        SearchInput<FlintMod> input;
        input.degreeD = cfg.degreeD;
        input.maxDeltaDegreeM = maxSearchDeg;
        input.numFBIMasters = cfg.numFBIMasters;
        input.targets = targets;
        auto samples = parseSeriesFile<FlintMod>(outputBase + "/series", targets, cfg.degreeD, 0);
        input.samples.insert(input.samples.end(), samples.begin(), samples.end());

        RelationSearcher<FlintMod> relationSearcher(input);
        auto relationResult = relationSearcher.search();

        CoefficientRelationExpander<FlintMod> expander;
        auto assignments = expander.expandAssignments(relationResult);
        auto fiRelations = expander.buildFIRelations(assignments);

        FIReductionSearcher<FlintMod> fiSearcher(fiRelations);
        auto fiReduction = fiSearcher.search();

        std::ofstream out(outputBase + "/fi_reduction_output");
        if (!out.is_open()) {
            throw std::runtime_error("Cannot open output file: " + outputBase + "/fi_reduction_output");
        }
        RelationFormatter<FlintMod>::writeFIReductionSummary(out, fiReduction);
        out << "\n";
        RelationFormatter<FlintMod>::writeFIMasterBasis(out, fiReduction);
        out << "\n";
        RelationFormatter<FlintMod>::writeFIReductions(out, fiReduction);

        std::cout << "FI solve pipeline finished.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
