template<typename ST>
std::vector<ST> DeltaProjector<ST>::powers(ST base, int degree) {
    std::vector<ST> result(static_cast<size_t>(degree + 1), ST(1));
    for (int i = 1; i <= degree; ++i) {
        result[static_cast<size_t>(i)] = result[static_cast<size_t>(i - 1)] * base;
    }
    return result;
}

template<typename ST>
DeltaProjector<ST>::DeltaProjector(const DeltaProjectionConfig<ST>& config)
    : config_(config) {
    if (config_.degree < 0) {
        throw std::invalid_argument("Projection degree must be non-negative");
    }

    const int degree = config_.degree;
    const ST one(1);
    const ST negA = ST(0) - config_.shiftA;
    const ST negB = ST(0) - config_.shiftB;
    const ST oneMinusA = one - config_.shiftA;
    const ST oneMinusB = one - config_.shiftB;

    xLower_ = powers(negA, degree + 1);
    xUpper_ = powers(oneMinusA, degree + 1);
    yLower_ = powers(negB, degree + 1);
    yUpper_ = powers(oneMinusB, degree + 1);

    intX_.assign(static_cast<size_t>(degree + 1), ST(0));
    intY_.assign(static_cast<size_t>(degree + 1), ST(0));
    for (int k = 0; k <= degree; ++k) {
        intX_[static_cast<size_t>(k)] =
            (xUpper_[static_cast<size_t>(k + 1)] - xLower_[static_cast<size_t>(k + 1)]) / ST(k + 1);
        intY_[static_cast<size_t>(k)] =
            (yUpper_[static_cast<size_t>(k + 1)] - yLower_[static_cast<size_t>(k + 1)]) / ST(k + 1);
    }
}

template<typename ST>
void DeltaProjector<ST>::validateTag(const IntegralTag& tag) {
    if (tag.head == IntegralHead::FI) {
        if (!tag.boundaries.empty()) {
            throw std::invalid_argument("FI target must not have boundary tags");
        }
        return;
    }
    if (tag.head == IntegralHead::BFI) {
        if (tag.boundaries.size() != 1) {
            throw std::invalid_argument("BFI target must have exactly one boundary tag");
        }
        return;
    }
    if (tag.boundaries.size() != 2) {
        throw std::invalid_argument("BBFI target must have exactly two boundary tags");
    }
    bool hasX = false;
    bool hasY = false;
    for (const auto& b : tag.boundaries) {
        hasX = hasX || b.axis == BoundaryAxis::X;
        hasY = hasY || b.axis == BoundaryAxis::Y;
    }
    if (!hasX || !hasY) {
        throw std::invalid_argument("BBFI target must have one X boundary and one Y boundary");
    }
}

template<typename ST>
const std::vector<ST>& DeltaProjector<ST>::boundaryWeights(const BoundaryTag& tag) const {
    if (tag.axis == BoundaryAxis::X) {
        return tag.side == BoundarySide::U ? xUpper_ : xLower_;
    }
    return tag.side == BoundarySide::U ? yUpper_ : yLower_;
}

template<typename ST>
std::vector<ST> DeltaProjector<ST>::project(const Series<ST>& series, const IntegralTag& tag) const {
    validateTag(tag);

    const int outputDegree = config_.degree;
    std::vector<ST> coeffs(static_cast<size_t>(outputDegree + 1), ST(0));
    const int maxDeg = std::min(config_.degree, series.getDeg());

    for (int d = 0; d <= maxDeg; ++d) {
        for (int p = 0; p <= d; ++p) {
            const int q = d - p;
            const ST c = series.getCoeff(p, q);
            if (c == ST(0)) continue;

            ST weight(0);
            int outDeg = 0;
            if (tag.head == IntegralHead::FI) {
                outDeg = p + q + 2;
                if (outDeg > outputDegree) continue;
                weight = intX_[static_cast<size_t>(p)] * intY_[static_cast<size_t>(q)];
            } else if (tag.head == IntegralHead::BFI) {
                outDeg = p + q + 1;
                if (outDeg > outputDegree) continue;
                const auto& b = tag.boundaries[0];
                if (b.axis == BoundaryAxis::X) {
                    weight = boundaryWeights(b)[static_cast<size_t>(p)] * intY_[static_cast<size_t>(q)];
                } else {
                    weight = intX_[static_cast<size_t>(p)] * boundaryWeights(b)[static_cast<size_t>(q)];
                }
            } else {
                outDeg = p + q;
                if (outDeg > outputDegree) continue;
                ST xWeight(1);
                ST yWeight(1);
                for (const auto& b : tag.boundaries) {
                    if (b.axis == BoundaryAxis::X) {
                        xWeight = boundaryWeights(b)[static_cast<size_t>(p)];
                    } else {
                        yWeight = boundaryWeights(b)[static_cast<size_t>(q)];
                    }
                }
                weight = xWeight * yWeight;
            }

            coeffs[static_cast<size_t>(outDeg)] += c * weight;
        }
    }
    return coeffs;
}
