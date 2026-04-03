template<typename T>
RelationSearcher<T>::RelationSearcher(const SearchInput<T>& input)
    : input_(input) {}

template<typename T>
RelationSearchResult<T> RelationSearcher<T>::search() const {
    RelationMatrixBuilder<T> builder(input_);
    RelationSearchResult<T> result;
    result.variables = builder.buildVariables();

    auto matrix = builder.buildMatrix(result.variables);
    LinearSystem<T> system(matrix);
    system.eliminate();

    result.rrefMatrix = system.getRREFMatrix();
    result.pivotColumns = system.getPivotColumns();
    result.freeColumns = system.getFreeVariableColumns();
    return result;
}
