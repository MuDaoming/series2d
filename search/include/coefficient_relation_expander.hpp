#pragma once

#include <vector>

#include "relation_types.hpp"

template<typename T>
class CoefficientRelationExpander {
public:
    std::vector<CoefficientAssignment<T>> expandAssignments(
        const RelationSearchResult<T>& result) const;

    std::vector<FIRelation<T>> buildFIRelations(
        const std::vector<CoefficientAssignment<T>>& assignments,
        const T& deltaValue) const;
};

#include "../src/coefficient_relation_expander.tpp"
