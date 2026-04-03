#pragma once

#include "linear.hpp"
#include "relation_matrix_builder.hpp"
#include "relation_types.hpp"

template<typename T>
class RelationSearcher {
public:
    explicit RelationSearcher(const SearchInput<T>& input);

    RelationSearchResult<T> search() const;

private:
    const SearchInput<T>& input_;
};

#include "../src/relation_searcher.tpp"
