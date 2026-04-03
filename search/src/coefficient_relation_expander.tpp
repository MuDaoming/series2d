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
std::vector<FIRelation<T>> CoefficientRelationExpander<T>::buildFIRelations(
    const std::vector<CoefficientAssignment<T>>& assignments) const {
    std::vector<FIRelation<T>> relations;
    relations.reserve(assignments.size());

    for (const auto& assignment : assignments) {
        FIRelation<T> relation;
        for (size_t i = 0; i < assignment.variables.size(); ++i) {
            const auto& integral = assignment.variables[i].integral;
            auto it = std::find_if(
                relation.integrals.begin(),
                relation.integrals.end(),
                [&](const IntegralLabel& label) {
                    return equalNu(label.nu, integral.nu);
                });

            if (it == relation.integrals.end()) {
                relation.integrals.push_back(integral);
                relation.coeffs.push_back(assignment.values[i]);
            } else {
                const size_t idx = static_cast<size_t>(it - relation.integrals.begin());
                relation.coeffs[idx] = relation.coeffs[idx] + assignment.values[i];
            }
        }
        relations.push_back(std::move(relation));
    }

    return relations;
}
