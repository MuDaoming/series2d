#pragma once

#include <vector>

#include "fi_reduction_builder.hpp"
#include "linear.hpp"
#include "relation_types.hpp"

template<typename T>
class FIReductionSearcher {
public:
    explicit FIReductionSearcher(const std::vector<FIRelation<T>>& relations);

    FIReductionResult<T> search() const;

private:
    const std::vector<FIRelation<T>>& relations_;
};

#include "../src/fi_reduction_searcher.tpp"
