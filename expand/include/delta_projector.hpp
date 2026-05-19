#pragma once

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "integral_tag.hpp"
#include "series.hpp"

template<typename ST>
struct DeltaProjectionConfig {
    ST shiftA;
    ST shiftB;
    int degree;
};

template<typename ST>
class DeltaProjector {
private:
    DeltaProjectionConfig<ST> config_;
    std::vector<ST> intX_;
    std::vector<ST> intY_;
    std::vector<ST> xUpper_;
    std::vector<ST> xLower_;
    std::vector<ST> yUpper_;
    std::vector<ST> yLower_;

public:
    explicit DeltaProjector(const DeltaProjectionConfig<ST>& config);

    std::vector<ST> project(const Series<ST>& series, const IntegralTag& tag) const;

private:
    static std::vector<ST> powers(ST base, int degree);
    static void validateTag(const IntegralTag& tag);
    const std::vector<ST>& boundaryWeights(const BoundaryTag& tag) const;
};

#include "../src/delta_projector.tpp"
