#pragma once

#include <algorithm>
#include <stdexcept>
#include <vector>
#include "series.hpp"

template<typename ST>
struct IntegrationConfig {
    ST shiftA_;
    ST shiftB_;
    int degree_;
};

template<typename ST>
class SeriesIntegrator {
private:
    IntegrationConfig<ST> config_;

public:
    explicit SeriesIntegrator(const IntegrationConfig<ST>& config);

    // Convert 2D series to 1D integrated coefficients.
    std::vector<ST> integrate(const Series<ST>& series) const;

private:
    ST monomialWeight(int xPower, int yPower) const;
    static ST powInt(ST base, int exp);
};

#include "../src/series_integrator.tpp"
