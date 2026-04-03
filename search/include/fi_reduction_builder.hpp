#pragma once

#include <vector>

#include "relation_types.hpp"

template<typename T>
class FIReductionBuilder {
public:
    explicit FIReductionBuilder(const std::vector<FIRelation<T>>& relations);

    std::vector<IntegralLabel> buildIntegralVariables() const;
    std::vector<std::vector<T>> buildMatrix(
        const std::vector<IntegralLabel>& integrals) const;

private:
    const std::vector<FIRelation<T>>& relations_;
};

#include "../src/fi_reduction_builder.tpp"
