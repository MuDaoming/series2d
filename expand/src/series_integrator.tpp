template<typename ST>
SeriesIntegrator<ST>::SeriesIntegrator(const IntegrationConfig<ST>& config)
    : config_(config) {
    if (config_.degree_ < 0) {
        throw std::invalid_argument("Integration degree must be non-negative");
    }
}

template<typename ST>
ST SeriesIntegrator<ST>::powInt(ST base, int exp) {
    ST result = ST(1);
    while (exp > 0) {
        if (exp & 1) {
            result *= base;
        }
        exp >>= 1;
        if (exp > 0) {
            base *= base;
        }
    }
    return result;
}

template<typename ST>
ST SeriesIntegrator<ST>::monomialWeight(int xPower, int yPower) const {
    const ST a = config_.shiftA_;
    const ST b = config_.shiftB_;

    const ST one = ST(1);
    const ST negA = ST(0) - a;
    const ST negB = ST(0) - b;
    const ST oneMinusA = one - a;
    const ST oneMinusB = one - b;

    const int px = xPower + 1;
    const int py = yPower + 1;

    const ST weightX = (powInt(oneMinusA, px) - powInt(negA, px)) / ST(px);
    const ST weightY = (powInt(oneMinusB, py) - powInt(negB, py)) / ST(py);
    return weightX * weightY;
}

template<typename ST>
std::vector<ST> SeriesIntegrator<ST>::integrate(const Series<ST>& series) const {
    std::vector<ST> coeffs(config_.degree_ + 1, ST(0));
    const int maxDeg = std::min(config_.degree_, series.getDeg());

    for (int d = 0; d <= maxDeg; ++d) {
        ST sum = ST(0);
        for (int p = 0; p <= d; ++p) {
            int q = d - p;
            ST c = series.getCoeff(p, q);
            if (c == ST(0)) continue;
            sum += c * monomialWeight(p, q);
        }
        coeffs[d] = sum;
    }

    return coeffs;
}
