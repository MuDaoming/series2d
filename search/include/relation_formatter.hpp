#pragma once

#include <ostream>

#include "relation_types.hpp"

template<typename T>
class RelationFormatter {
public:
    static void writeSummary(
        std::ostream& out,
        const RelationSearchResult<T>& result);

    static void writeRelations(
        std::ostream& out,
        const RelationSearchResult<T>& result);

    static void writeRREF(
        std::ostream& out,
        const RelationSearchResult<T>& result);

    static void writeAssignments(
        std::ostream& out,
        const std::vector<CoefficientAssignment<T>>& assignments);

    static void writeIntegralRelations(
        std::ostream& out,
        const std::vector<IntegralRelation<T>>& relations);

    static void writeIntegralReductionSummary(
        std::ostream& out,
        const IntegralReductionResult<T>& result);

    static void writeIntegralMasterBasis(
        std::ostream& out,
        const IntegralReductionResult<T>& result);

    static void writeIntegralReductions(
        std::ostream& out,
        const IntegralReductionResult<T>& result);

    static void writeIntegralRREF(
        std::ostream& out,
        const IntegralReductionResult<T>& result);
};

#include "../src/relation_formatter.tpp"
