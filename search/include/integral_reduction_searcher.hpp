#pragma once

#include <vector>

#include "integral_reduction_builder.hpp"
#include "linear.hpp"
#include "relation_types.hpp"

template<typename T>
class IntegralReductionSearcher {
public:
    explicit IntegralReductionSearcher(const std::vector<IntegralRelation<T>>& relations);

    IntegralReductionResult<T> search() const;

private:
    const std::vector<IntegralRelation<T>>& relations_;
};

#include "../src/integral_reduction_searcher.tpp"
