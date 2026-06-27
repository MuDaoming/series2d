#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "ff_type.hpp"
#include "io.hpp"
#include "relation_certifier.hpp"
#include "relation_matrix_builder.hpp"
#include "relation_searcher.hpp"
#include "relation_types.hpp"

struct DEConfig {
    int nuSize = 0;
    int degreeD = 0;
    int maxM = 0;
    int ncheck = 1;
    mp_limb_t p = 0;
};

static std::string trim(const std::string& s) {
    size_t l = 0;
    while (l < s.size() && std::isspace(static_cast<unsigned char>(s[l]))) ++l;
    size_t r = s.size();
    while (r > l && std::isspace(static_cast<unsigned char>(s[r - 1]))) --r;
    return s.substr(l, r - l);
}

static DEConfig parseConfig(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open config: " + path);
    }

    std::unordered_map<std::string, std::string> kv;
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) continue;
        const size_t eq = line.find('=');
        if (eq == std::string::npos) continue;
        kv[trim(line.substr(0, eq))] = trim(line.substr(eq + 1));
    }

    auto need = [&](const std::string& key) -> std::string {
        auto it = kv.find(key);
        if (it == kv.end()) {
            throw std::runtime_error("Missing key in config: " + key);
        }
        return it->second;
    };

    DEConfig cfg;
    cfg.nuSize = std::stoi(need("N"));
    cfg.degreeD = std::stoi(need("deg"));
    cfg.maxM = std::stoi(need("m"));
    {
        auto it = kv.find("ncheck");
        if (it != kv.end()) cfg.ncheck = std::stoi(it->second);
    }
    if (cfg.ncheck < 0) {
        throw std::runtime_error("ncheck must be >= 0");
    }
    cfg.p = static_cast<mp_limb_t>(std::stoull(need("p")));
    return cfg;
}

static std::filesystem::path resolvePath(
    const std::filesystem::path& baseDir,
    const std::string& raw) {
    std::filesystem::path p(raw);
    if (p.is_absolute()) {
        return p;
    }
    return baseDir / p;
}

static std::vector<std::pair<std::filesystem::path, std::filesystem::path>>
parseSeriesList(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open series list: " + path);
    }

    const auto baseDir = std::filesystem::absolute(path).parent_path();
    std::vector<std::pair<std::filesystem::path, std::filesystem::path>> entries;
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::string seriesPath;
        std::string targetPath;
        iss >> seriesPath >> targetPath;
        if (seriesPath.empty() || targetPath.empty()) {
            throw std::runtime_error("Bad series list line: " + line);
        }
        entries.emplace_back(resolvePath(baseDir, seriesPath),
                             resolvePath(baseDir, targetPath));
    }
    if (entries.empty()) {
        throw std::runtime_error("series list is empty: " + path);
    }
    return entries;
}

static std::vector<IntegralLabel> parseMastersFile(
    const std::string& path,
    int expectedNuSize) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open masters file: " + path);
    }

    std::vector<IntegralLabel> masters;
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        const size_t l = line.find('{');
        const size_t r = line.rfind('}');
        if (l == std::string::npos || r == std::string::npos || r <= l) {
            throw std::runtime_error("Invalid master label: " + line);
        }

        IntegralLabel label;
        const std::string prefix = trim(line.substr(0, l));
        if (prefix.empty()) {
            label.head = IntegralHead::FI;
        } else {
            const size_t lb = prefix.find('[');
            if (lb == std::string::npos) {
                label.head = parseSearchIntegralHead(prefix);
            } else {
                const size_t rb = prefix.rfind(']');
                if (rb == std::string::npos || rb <= lb || rb + 1 != prefix.size()) {
                    throw std::runtime_error("Invalid boundary list in master label: " + line);
                }
                label.head = parseSearchIntegralHead(trim(prefix.substr(0, lb)));
                label.boundaries = parseSearchBoundaryList(prefix.substr(lb + 1, rb - lb - 1));
            }
        }

        std::string body = line.substr(l + 1, r - l - 1);
        std::stringstream ss(body);
        std::string tok;
        while (std::getline(ss, tok, ',')) {
            tok = trim(tok);
            if (!tok.empty()) {
                label.nu.push_back(std::stoi(tok));
            }
        }
        if (static_cast<int>(label.nu.size()) != expectedNuSize) {
            throw std::runtime_error("nu length mismatch in master label: " + line);
        }
        validateSearchIntegralLabel(label, line);
        masters.push_back(std::move(label));
    }

    if (masters.empty()) {
        throw std::runtime_error("No masters in file: " + path);
    }
    return masters;
}

static IntegralLabel shiftedLabel(const IntegralLabel& label, int shift) {
    IntegralLabel out = label;
    for (int& x : out.nu) {
        x += shift;
    }
    return out;
}

static std::vector<FlintMod> parseSeriesLine(const std::string& line) {
    const size_t l = line.find('{');
    const size_t r = line.rfind('}');
    if (l == std::string::npos || r == std::string::npos || r <= l) {
        throw std::runtime_error("Bad series line: " + line);
    }
    std::vector<FlintMod> coeffs;
    std::string body = line.substr(l + 1, r - l - 1);
    std::istringstream iss(body);
    std::string tok;
    while (std::getline(iss, tok, ',')) {
        tok = trim(tok);
        if (!tok.empty()) {
            coeffs.emplace_back(static_cast<unsigned long long>(std::stoull(tok)));
        }
    }
    return coeffs;
}

using SeriesByLabel = std::map<IntegralLabel, std::vector<FlintMod>, IntegralLabelLess>;

static SeriesByLabel readSeriesForLabels(
    const std::filesystem::path& seriesPath,
    const std::filesystem::path& targetPath,
    const std::vector<IntegralLabel>& labels,
    int nuSize,
    int degreeD) {
    const auto target = parseSearchTargetFile(targetPath.string(), nuSize);

    std::map<IntegralLabel, int, IntegralLabelLess> wanted;
    for (int i = 0; i < static_cast<int>(labels.size()); ++i) {
        wanted[labels[i]] = i;
    }

    std::vector<int> targetToWanted(target.size(), -1);
    for (int i = 0; i < static_cast<int>(target.size()); ++i) {
        auto it = wanted.find(target[i]);
        if (it != wanted.end()) {
            targetToWanted[i] = it->second;
        }
    }

    std::ifstream in(seriesPath);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open series file: " + seriesPath.string());
    }

    SeriesByLabel result;
    std::string line;
    int lineIdx = 0;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) continue;
        if (lineIdx >= static_cast<int>(target.size())) {
            throw std::runtime_error("Too many series lines in: " + seriesPath.string());
        }
        const int wantedIdx = targetToWanted[lineIdx];
        if (wantedIdx >= 0) {
            auto coeffs = parseSeriesLine(line);
            if (static_cast<int>(coeffs.size()) < degreeD + 1) {
                throw std::runtime_error("Series shorter than deg+1 in: " + seriesPath.string());
            }
            coeffs.resize(degreeD + 1);
            result[labels[wantedIdx]] = std::move(coeffs);
        }
        ++lineIdx;
    }

    if (lineIdx != static_cast<int>(target.size())) {
        throw std::runtime_error("Series line count mismatch in: " + seriesPath.string());
    }

    for (const auto& label : labels) {
        if (result.find(label) == result.end()) {
            throw std::runtime_error(
                "Missing master series: " + integralLabelToString(label));
        }
    }
    return result;
}

static bool sameVariable(const RelationVariable& var, const IntegralLabel& label) {
    return equalIntegralLabel(var.integral, label);
}

static void validateTargetsInMasters(
    const std::vector<IntegralLabel>& deTargets,
    const std::vector<IntegralLabel>& masters) {
    std::map<IntegralLabel, bool, IntegralLabelLess> masterSet;
    for (const auto& master : masters) {
        masterSet[master] = true;
    }
    for (const auto& target : deTargets) {
        if (masterSet.find(target) == masterSet.end()) {
            throw std::runtime_error(
                "G_path contains a DE target not present in masters_path: " +
                integralLabelToString(target));
        }
    }
}

static int findChosenFreeColumn(
    const RelationSearchResult<FlintMod>& result,
    const IntegralLabel& derivativeLabel) {
    for (int freeCol : result.freeColumns) {
        if (freeCol < 0 || freeCol >= static_cast<int>(result.variables.size())) continue;
        if (sameVariable(result.variables[freeCol], derivativeLabel)) {
            return freeCol;
        }
        for (size_t row = 0; row < result.pivotColumns.size(); ++row) {
            const int pivotCol = result.pivotColumns[row];
            if (pivotCol < 0) continue;
            if (!sameVariable(result.variables[pivotCol], derivativeLabel)) continue;
            if (result.rrefMatrix[row][freeCol] != FlintMod(0ULL)) {
                return freeCol;
            }
        }
    }
    return -1;
}

static std::vector<FlintMod> buildCoefficientVector(
    const RelationSearchResult<FlintMod>& result,
    int chosenFreeColumn) {
    std::vector<FlintMod> coeffs(result.variables.size(), FlintMod(0ULL));
    coeffs[chosenFreeColumn] = FlintMod(1ULL);
    for (size_t row = 0; row < result.pivotColumns.size(); ++row) {
        const int pivotCol = result.pivotColumns[row];
        if (pivotCol < 0) continue;
        coeffs[pivotCol] = FlintMod(0ULL) - result.rrefMatrix[row][chosenFreeColumn];
    }
    return coeffs;
}

static std::vector<FlintMod> polynomialForLabel(
    const std::vector<RelationVariable>& variables,
    const std::vector<FlintMod>& coeffs,
    const IntegralLabel& label,
    int m) {
    std::vector<FlintMod> poly(static_cast<size_t>(m + 1), FlintMod(0ULL));
    for (size_t i = 0; i < variables.size(); ++i) {
        if (sameVariable(variables[i], label)) {
            poly[variables[i].k] = coeffs[i];
        }
    }
    return poly;
}

static bool nonzeroPolynomial(const std::vector<FlintMod>& poly) {
    for (const auto& c : poly) {
        if (c != FlintMod(0ULL)) return true;
    }
    return false;
}

static std::string polynomialToString(const std::vector<FlintMod>& poly) {
    std::ostringstream out;
    bool first = true;
    for (int k = 0; k < static_cast<int>(poly.size()); ++k) {
        if (poly[k] == FlintMod(0ULL)) continue;
        if (!first) {
            out << " + ";
        }
        out << poly[k];
        if (k == 1) {
            out << "*delta";
        } else if (k > 1) {
            out << "*delta^" << k;
        }
        first = false;
    }
    if (first) {
        return "0";
    }
    return out.str();
}

static std::vector<int> buildExponentialMSchedule(int maxM) {
    if (maxM < 0) {
        throw std::runtime_error("m must be >= 0");
    }

    std::vector<int> schedule;
    schedule.push_back(0);
    if (maxM == 0) {
        return schedule;
    }

    for (int m = 1; m < maxM;) {
        schedule.push_back(m);
        if (m > maxM / 2) {
            break;
        }
        m *= 2;
    }
    if (schedule.back() != maxM) {
        schedule.push_back(maxM);
    }
    return schedule;
}

int main(int argc, char** argv) {
    if (argc != 6 && argc != 7) {
        std::cerr << "Usage: " << argv[0]
                  << " <config_path> <G_path> <series_list_path> <masters_path> <output_path> [shift]\n";
        return 1;
    }

    const std::string configPath = argv[1];
    const std::string seriesListPath = argv[3];
    const std::string mastersPath = argv[4];
    const std::filesystem::path outputPath(argv[5]);
    const int shift = (argc == 7) ? std::stoi(argv[6]) : 100;

    try {
        const DEConfig cfg = parseConfig(configPath);
        FlintMod::set_modulus(cfg.p);

        const auto seriesList = parseSeriesList(seriesListPath);
        const int numBC = static_cast<int>(seriesList.size());
        const int derivativeDegree = cfg.degreeD - 1;
        if (derivativeDegree < 0) {
            throw std::runtime_error("config deg must be positive");
        }
        const DegreeWindow window = makeDegreeWindow(derivativeDegree, cfg.ncheck);

        const auto deTargets = parseSearchTargetFile(argv[2], cfg.nuSize);
        const auto masters = parseMastersFile(mastersPath, cfg.nuSize);
        validateTargetsInMasters(deTargets, masters);

        std::vector<SeriesByLabel> masterSeriesByBC;
        masterSeriesByBC.reserve(numBC);
        for (const auto& entry : seriesList) {
            masterSeriesByBC.push_back(readSeriesForLabels(
                entry.first, entry.second, masters, cfg.nuSize, cfg.degreeD));
        }

        if (outputPath.has_parent_path()) {
            std::filesystem::create_directories(outputPath.parent_path());
        }
        std::ofstream out(outputPath);
        if (!out.is_open()) {
            throw std::runtime_error("Cannot open output: " + outputPath.string());
        }

        out << "# p = " << cfg.p << "\n";
        out << "# deg = " << derivativeDegree << "\n";
        out << "# ncheck = " << cfg.ncheck << "\n";
        out << "# train_deg = " << window.trainDeg << "\n";
        out << "# m_max = " << cfg.maxM << "\n";
        out << "# m_schedule = exponential\n";
        out << "# shift = " << shift << "\n";
        out << "# de_targets = " << deTargets.size() << "\n";
        out << "# masters = " << masters.size() << "\n\n";
        out.flush();

        const auto mSchedule = buildExponentialMSchedule(cfg.maxM);

        for (size_t masterIdx = 0; masterIdx < deTargets.size(); ++masterIdx) {
            const auto& master = deTargets[masterIdx];
            const std::string masterName = integralLabelToString(master);
            std::cout << "[search_de] target " << (masterIdx + 1)
                      << "/" << deTargets.size() << ": "
                      << masterName << "\n";

            const IntegralLabel dLabel = shiftedLabel(master, shift);
            std::vector<IntegralLabel> labels;
            labels.reserve(masters.size() + 1);
            labels.push_back(dLabel);
            labels.insert(labels.end(), masters.begin(), masters.end());

            std::vector<SeriesSample<FlintMod>> samples;
            samples.reserve(static_cast<size_t>(numBC) * labels.size());
            for (int bc = 0; bc < numBC; ++bc) {
                const auto& original = masterSeriesByBC[bc].at(master);
                SeriesSample<FlintMod> dSample;
                dSample.label.integral = dLabel;
                dSample.label.bcIndex = bc;
                dSample.coeffs.reserve(derivativeDegree + 1);
                for (int n = 0; n <= derivativeDegree; ++n) {
                    dSample.coeffs.push_back(
                        FlintMod(static_cast<unsigned long long>(n + 1)) *
                        original[n + 1]);
                }
                samples.push_back(std::move(dSample));

                for (const auto& mlabel : masters) {
                    SeriesSample<FlintMod> sample;
                    sample.label.integral = mlabel;
                    sample.label.bcIndex = bc;
                    sample.coeffs = masterSeriesByBC[bc].at(mlabel);
                    sample.coeffs.resize(derivativeDegree + 1);
                    samples.push_back(std::move(sample));
                }
            }

            bool found = false;
            int foundM = -1;
            std::vector<FlintMod> foundCoeffVector;
            RelationSearchResult<FlintMod> foundResult;
            int failedCheckM = -1;

            for (int m : mSchedule) {
                std::cout << "[search_de]   m = " << m << "\n";
                SearchInput<FlintMod> trainInput;
                trainInput.degreeD = window.trainDeg;
                trainInput.maxDeltaDegreeM = m;
                trainInput.numFBIMasters = numBC;
                trainInput.targets = labels;
                trainInput.samples = samples;

                RelationSearcher<FlintMod> trainSearcher(trainInput);
                const auto trainResult = trainSearcher.search();
                const int chosenFree = findChosenFreeColumn(trainResult, dLabel);
                if (chosenFree < 0) {
                    continue;
                }

                SearchInput<FlintMod> fullInput;
                fullInput.degreeD = derivativeDegree;
                fullInput.maxDeltaDegreeM = m;
                fullInput.numFBIMasters = numBC;
                fullInput.targets = labels;
                fullInput.samples = samples;

                const auto checkRows = buildCheckRows(
                    fullInput,
                    trainResult.variables,
                    window.checkStart,
                    window.checkEnd);
                if (checkNullspaceShrink(
                        trainResult.rrefMatrix,
                        trainResult.pivotColumns,
                        checkRows)) {
                    failedCheckM = m;
                    continue;
                }

                found = true;
                foundM = m;
                foundCoeffVector = buildCoefficientVector(trainResult, chosenFree);
                foundResult = trainResult;
                break;
            }

            out << "## d(" << masterName << ")\n";
            out << "# shifted label: " << integralLabelToString(dLabel) << "\n";
            if (!found) {
                if (failedCheckM >= 0) {
                    out << "# train nullspace failed held-out check, last failed m = "
                        << failedCheckM << "\n";
                }
                out << "not found for m <= " << cfg.maxM << "\n\n";
                out.flush();
                std::cout << "[search_de]   not found for m <= "
                          << cfg.maxM << "\n";
                continue;
            }
            std::cout << "[search_de]   found at m = " << foundM << "\n";

            std::vector<std::vector<FlintMod>> polys;
            polys.reserve(labels.size());
            for (const auto& label : labels) {
                polys.push_back(polynomialForLabel(
                    foundResult.variables, foundCoeffVector, label, foundM));
            }

            const bool trainVerified = verifyRelationOnWindow(
                labels, polys, samples, numBC, 0, window.trainDeg);
            const bool checkVerified = verifyRelationOnWindow(
                labels, polys, samples, numBC, window.checkStart, window.checkEnd);
            out << "# m = " << foundM << "\n";
            out << "# train_verified = " << (trainVerified ? 1 : 0) << "\n";
            out << "# check_verified = " << (checkVerified ? 1 : 0) << "\n";
            out << "# certified = " << ((trainVerified && checkVerified) ? 1 : 0) << "\n";

            bool first = true;
            for (size_t i = 0; i < labels.size(); ++i) {
                if (!nonzeroPolynomial(polys[i])) continue;
                if (!first) {
                    out << " +\n";
                }
                out << "(" << polynomialToString(polys[i]) << ")*";
                if (i == 0) {
                    out << "d(" << integralLabelToString(master) << ")";
                } else {
                    out << integralLabelToString(labels[i]);
                }
                first = false;
            }
            if (first) {
                out << "0";
            }
            out << " = 0\n\n";
            out.flush();
        }

        std::cout << "[search_de] wrote " << outputPath << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
