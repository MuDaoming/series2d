#include <algorithm>
#include <map>
#include <stdexcept>

template<typename T>
RelationMatrixBuilder<T>::RelationMatrixBuilder(const SearchInput<T>& input)
    : input_(input) {}

template<typename T>
std::vector<RelationVariable> RelationMatrixBuilder<T>::buildVariables() const {
    std::vector<RelationVariable> vars;
    vars.reserve(input_.targets.size() * static_cast<size_t>(input_.maxDeltaDegreeM + 1));

    for (const auto& target : input_.targets) {
        for (int k = 0; k <= input_.maxDeltaDegreeM; ++k) {
            RelationVariable var;
            var.integral = target;
            var.k = k;
            vars.push_back(std::move(var));
        }
    }

    std::sort(vars.begin(), vars.end(), RelationVariableMoreComplexFirst{});
    return vars;
}

template<typename T>
const SeriesSample<T>& RelationMatrixBuilder<T>::findSample(
    const IntegralLabel& label,
    int bcIndex) const {
    for (const auto& sample : input_.samples) {
        if (sample.label.bcIndex == bcIndex &&
            equalIntegralLabel(sample.label.integral, label)) {
            return sample;
        }
    }
    throw std::runtime_error(
        "Missing series sample for bcIndex=" + std::to_string(bcIndex) +
        ", integral=" + integralLabelToString(label));
}

template<typename T>
std::vector<std::vector<T>> RelationMatrixBuilder<T>::buildMatrix(
    const std::vector<RelationVariable>& variables) const {
    return buildRowsForDegreeWindow(variables, 0, input_.degreeD);
}

template<typename T>
std::vector<std::vector<T>> RelationMatrixBuilder<T>::buildRowsForDegreeWindow(
    const std::vector<RelationVariable>& variables,
    int startDegree,
    int endDegree) const {
    std::vector<std::vector<T>> matrix;
    if (startDegree > endDegree) {
        return matrix;
    }
    if (startDegree < 0 || endDegree > input_.degreeD) {
        throw std::runtime_error("requested degree window is outside input degree range");
    }

    std::map<IntegralLabel, int, IntegralLabelLess> targetIndex;
    for (int i = 0; i < static_cast<int>(input_.targets.size()); ++i) {
        targetIndex[input_.targets[i]] = i;
    }

    std::vector<int> variableTargetIndex(variables.size(), -1);
    for (int col = 0; col < static_cast<int>(variables.size()); ++col) {
        auto it = targetIndex.find(variables[col].integral);
        if (it == targetIndex.end()) {
            throw std::runtime_error(
                "Relation variable integral is not in input targets: " +
                integralLabelToString(variables[col].integral));
        }
        variableTargetIndex[col] = it->second;
    }

    std::vector<std::vector<const SeriesSample<T>*>> samplesByBCAndTarget(
        input_.numFBIMasters,
        std::vector<const SeriesSample<T>*>(input_.targets.size(), nullptr));
    for (int bcIndex = 0; bcIndex < input_.numFBIMasters; ++bcIndex) {
        for (int target = 0; target < static_cast<int>(input_.targets.size()); ++target) {
            samplesByBCAndTarget[bcIndex][target] =
                &findSample(input_.targets[target], bcIndex);
        }
    }

    matrix.reserve(static_cast<size_t>(
        input_.numFBIMasters * (endDegree - startDegree + 1)));

    for (int bcIndex = 0; bcIndex < input_.numFBIMasters; ++bcIndex) {
        for (int n = startDegree; n <= endDegree; ++n) {
            std::vector<T> row(variables.size(), T(0));
            for (size_t col = 0; col < variables.size(); ++col) {
                const auto& var = variables[col];
                if (n < var.k) {
                    continue;
                }
                const auto& sample =
                    *samplesByBCAndTarget[bcIndex][variableTargetIndex[col]];
                if (n - var.k >= static_cast<int>(sample.coeffs.size())) {
                    throw std::runtime_error("series coeff count too short for matrix row");
                }
                row[col] = sample.coeffs[n - var.k];
            }
            matrix.push_back(std::move(row));
        }
    }

    return matrix;
}
