#include <algorithm>

template<typename T>
FIReductionBuilder<T>::FIReductionBuilder(const std::vector<FIRelation<T>>& relations)
    : relations_(relations) {}

template<typename T>
std::vector<IntegralLabel> FIReductionBuilder<T>::buildIntegralVariables() const {
    std::vector<IntegralLabel> integrals;
    for (const auto& relation : relations_) {
        for (const auto& integral : relation.integrals) {
            auto it = std::find_if(
                integrals.begin(),
                integrals.end(),
                [&](const IntegralLabel& label) {
                    return equalNu(label.nu, integral.nu);
                });
            if (it == integrals.end()) {
                integrals.push_back(integral);
            }
        }
    }

    std::sort(integrals.begin(), integrals.end(), IntegralLabelMoreComplexFirst{});
    return integrals;
}

template<typename T>
std::vector<std::vector<T>> FIReductionBuilder<T>::buildMatrix(
    const std::vector<IntegralLabel>& integrals) const {
    std::vector<std::vector<T>> matrix;
    matrix.reserve(relations_.size());

    for (const auto& relation : relations_) {
        std::vector<T> row(integrals.size(), T(0));
        for (size_t col = 0; col < integrals.size(); ++col) {
            auto it = std::find_if(
                relation.integrals.begin(),
                relation.integrals.end(),
                [&](const IntegralLabel& label) {
                    return equalNu(label.nu, integrals[col].nu);
                });
            if (it != relation.integrals.end()) {
                const size_t idx = static_cast<size_t>(it - relation.integrals.begin());
                row[col] = relation.coeffs[idx];
            }
        }
        matrix.push_back(std::move(row));
    }

    return matrix;
}
