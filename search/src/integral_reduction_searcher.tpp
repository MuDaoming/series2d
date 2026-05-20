template<typename T>
IntegralReductionSearcher<T>::IntegralReductionSearcher(const std::vector<IntegralRelation<T>>& relations)
    : relations_(relations) {}

template<typename T>
IntegralReductionResult<T> IntegralReductionSearcher<T>::search() const {
    IntegralReductionBuilder<T> builder(relations_);
    IntegralReductionResult<T> result;
    result.integrals = builder.buildIntegralVariables();

    auto matrix = builder.buildMatrix(result.integrals);
    LinearSystem<T> system(matrix);
    system.eliminate();

    result.rrefMatrix = system.getRREFMatrix();
    result.pivotColumns = system.getPivotColumns();
    result.freeColumns = system.getFreeVariableColumns();
    return result;
}
