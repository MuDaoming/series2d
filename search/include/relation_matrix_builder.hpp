#pragma once

#include <vector>

#include "relation_types.hpp"

template<typename T>
class RelationMatrixBuilder {
public:
    explicit RelationMatrixBuilder(const SearchInput<T>& input);

    std::vector<RelationVariable> buildVariables() const;
    std::vector<std::vector<T>> buildMatrix(
        const std::vector<RelationVariable>& variables) const;

private:
    const SearchInput<T>& input_;

    const SeriesSample<T>& findSample(const IntegralLabel& label, int bcIndex) const;
};

#include "../src/relation_matrix_builder.tpp"
