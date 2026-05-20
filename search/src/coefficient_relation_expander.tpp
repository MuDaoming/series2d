#include <algorithm>

template<typename T>
std::vector<CoefficientAssignment<T>> CoefficientRelationExpander<T>::expandAssignments(
    const RelationSearchResult<T>& result) const {
    std::vector<CoefficientAssignment<T>> assignments;
    assignments.reserve(result.freeColumns.size());

    for (int freeCol : result.freeColumns) {
        CoefficientAssignment<T> assignment;
        assignment.variables = result.variables;
        assignment.values.assign(result.variables.size(), T(0));
        assignment.chosenFreeColumn = freeCol;

        if (freeCol >= 0 && freeCol < static_cast<int>(assignment.values.size())) {
            assignment.values[freeCol] = T(1);
        }

        for (int row = static_cast<int>(result.pivotColumns.size()) - 1; row >= 0; --row) {
            const int pivotCol = result.pivotColumns[row];
            if (pivotCol < 0) {
                continue;
            }

            T sum = T(0);
            for (size_t col = static_cast<size_t>(pivotCol + 1); col < result.variables.size(); ++col) {
                sum = sum + result.rrefMatrix[row][col] * assignment.values[col];
            }
            assignment.values[pivotCol] = T(0) - sum;
        }

        assignments.push_back(std::move(assignment));
    }

    return assignments;
}

template<typename T>
std::vector<IntegralRelation<T>> CoefficientRelationExpander<T>::buildIntegralRelations(
    const std::vector<CoefficientAssignment<T>>& assignments,
    const T& deltaValue) const {
    std::vector<IntegralRelation<T>> relations;
    relations.reserve(assignments.size());

    if (assignments.empty()) {
        return relations;
    }

    const auto& variables = assignments.front().variables;
    const size_t numVars = variables.size();

    int maxK = 0;
    for (const auto& var : variables) {
        maxK = std::max(maxK, var.k);
    }

    std::vector<T> deltaPowers(static_cast<size_t>(maxK + 1), T(1));
    for (int k = 1; k <= maxK; ++k) {
        deltaPowers[static_cast<size_t>(k)] =
            deltaPowers[static_cast<size_t>(k - 1)] * deltaValue;
    }

    std::vector<T> varWeights;
    varWeights.reserve(numVars);
    for (const auto& var : variables) {
        varWeights.push_back(deltaPowers[static_cast<size_t>(var.k)]);
    }

    std::vector<IntegralLabel> integrals;
    std::vector<size_t> varToIntegralIdx(numVars, 0);
    for (size_t i = 0; i < numVars; ++i) {
        const auto& integral = variables[i].integral;
        auto it = std::find_if(
            integrals.begin(),
            integrals.end(),
            [&](const IntegralLabel& label) {
                return equalIntegralLabel(label, integral);
            });

        if (it == integrals.end()) {
            varToIntegralIdx[i] = integrals.size();
            integrals.push_back(integral);
        } else {
            varToIntegralIdx[i] = static_cast<size_t>(it - integrals.begin());
        }
    }

    for (const auto& assignment : assignments) {
        IntegralRelation<T> relation;
        relation.integrals = integrals;
        relation.coeffs.assign(integrals.size(), T(0));

        for (size_t i = 0; i < numVars; ++i) {
            const T& value = assignment.values[i];
            if (value == T(0)) {
                continue;
            }

            const size_t idx = varToIntegralIdx[i];
            relation.coeffs[idx] = relation.coeffs[idx] + value * varWeights[i];
        }
        relations.push_back(std::move(relation));
    }

    return relations;
}
