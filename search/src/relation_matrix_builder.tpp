#include <algorithm>

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
    std::vector<std::vector<T>> matrix;
    matrix.reserve(static_cast<size_t>(input_.numFBIMasters * (input_.degreeD + 1)));

    for (int bcIndex = 0; bcIndex < input_.numFBIMasters; ++bcIndex) {
        for (int n = 0; n <= input_.degreeD; ++n) {
            std::vector<T> row(variables.size(), T(0));
            for (size_t col = 0; col < variables.size(); ++col) {
                const auto& var = variables[col];
                if (n < var.k) {
                    continue;
                }
                const auto& sample = findSample(var.integral, bcIndex);
                row[col] = sample.coeffs[n - var.k];
            }
            matrix.push_back(std::move(row));
        }
    }

    return matrix;
}
