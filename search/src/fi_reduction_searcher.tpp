template<typename T>
FIReductionSearcher<T>::FIReductionSearcher(const std::vector<FIRelation<T>>& relations)
    : relations_(relations) {}

template<typename T>
FIReductionResult<T> FIReductionSearcher<T>::search() const {
    FIReductionBuilder<T> builder(relations_);
    FIReductionResult<T> result;
    result.integrals = builder.buildIntegralVariables();

    auto matrix = builder.buildMatrix(result.integrals);
    LinearSystem<T> system(matrix);
    system.eliminate();

    result.rrefMatrix = system.getRREFMatrix();
    result.pivotColumns = system.getPivotColumns();
    result.freeColumns = system.getFreeVariableColumns();
    return result;
}
