#include <filesystem>
#include <fstream>

void runRelationSearchPipeline(const std::vector<std::string>& seriesPaths,
                               const std::string& configPath,
                               const std::string& targetPath,
                               const std::string& maxSearchDegPath,
                               const std::string& outputPath) {
    namespace fs = std::filesystem;

    fs::path cfgP = fs::path(configPath);
    fs::path tgtP = fs::path(targetPath);
    fs::path degP = fs::path(maxSearchDegPath);
    fs::path outP = fs::path(outputPath);
    if (!cfgP.is_absolute()) cfgP = fs::current_path() / cfgP;
    if (!tgtP.is_absolute()) tgtP = fs::current_path() / tgtP;
    if (!degP.is_absolute()) degP = fs::current_path() / degP;
    if (!outP.is_absolute()) outP = fs::current_path() / outP;
    if (outP.has_parent_path()) {
        fs::create_directories(outP.parent_path());
    }

    SearchConfig cfg = parseSearchConfigFile(cfgP.string());
    FlintMod::set_modulus(cfg.p);
    const int maxSearchDeg = parseMaxSearchDegreeFile(degP.string());
    auto targets = parseSearchTargetFile(tgtP.string(), cfg.nuSize);

    if (static_cast<int>(seriesPaths.size()) != cfg.numFBIMasters) {
        throw std::runtime_error(
            "number of series files must equal numFBIMasters");
    }

    SearchInput<FlintMod> input;
    input.degreeD = cfg.degreeD;
    input.maxDeltaDegreeM = maxSearchDeg;
    input.numFBIMasters = cfg.numFBIMasters;
    input.targets = targets;

    for (int bcIndex = 0; bcIndex < cfg.numFBIMasters; ++bcIndex) {
        fs::path seriesP = fs::path(seriesPaths[bcIndex]);
        if (!seriesP.is_absolute()) seriesP = fs::current_path() / seriesP;
        auto samples = parseSeriesFile<FlintMod>(seriesP.string(), targets, cfg.degreeD, bcIndex);
        input.samples.insert(input.samples.end(), samples.begin(), samples.end());
    }

    RelationSearcher<FlintMod> searcher(input);
    auto result = searcher.search();
    CoefficientRelationExpander<FlintMod> expander;
    auto assignments = expander.expandAssignments(result);
    auto fiRelations = expander.buildFIRelations(assignments);
    FIReductionSearcher<FlintMod> fiSearcher(fiRelations);
    auto fiReduction = fiSearcher.search();

    std::ofstream out(outP.string());
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open output file: " + outP.string());
    }

    RelationFormatter<FlintMod>::writeSummary(out, result);
    out << "\n";
    RelationFormatter<FlintMod>::writeRelations(out, result);
    out << "\n";
    RelationFormatter<FlintMod>::writeAssignments(out, assignments);
    out << "\n";
    RelationFormatter<FlintMod>::writeFIRelations(out, fiRelations);
    out << "\n";
    RelationFormatter<FlintMod>::writeFIReductionSummary(out, fiReduction);
    out << "\n";
    RelationFormatter<FlintMod>::writeFIReductions(out, fiReduction);
    out << "\n";
    RelationFormatter<FlintMod>::writeFIRREF(out, fiReduction);
    out << "\n";
    RelationFormatter<FlintMod>::writeRREF(out, result);
}
