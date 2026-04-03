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

    static void writeFIRelations(
        std::ostream& out,
        const std::vector<FIRelation<T>>& relations);

    static void writeFIReductionSummary(
        std::ostream& out,
        const FIReductionResult<T>& result);

    static void writeFIMasterBasis(
        std::ostream& out,
        const FIReductionResult<T>& result);

    static void writeFIReductions(
        std::ostream& out,
        const FIReductionResult<T>& result);

    static void writeFIRREF(
        std::ostream& out,
        const FIReductionResult<T>& result);
};

#include "../src/relation_formatter.tpp"
