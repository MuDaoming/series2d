#pragma once

#include <vector>

#include "relation_types.hpp"

template<typename T>
class IntegralReductionBuilder {
public:
    explicit IntegralReductionBuilder(const std::vector<IntegralRelation<T>>& relations);

    std::vector<IntegralLabel> buildIntegralVariables() const;
    std::vector<std::vector<T>> buildMatrix(
        const std::vector<IntegralLabel>& integrals) const;

private:
    const std::vector<IntegralRelation<T>>& relations_;
};

#include "../src/integral_reduction_builder.tpp"
