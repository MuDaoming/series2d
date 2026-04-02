template<typename RT, typename PT, typename ST>
IntegrandExpander<RT, PT, ST>::IntegrandExpander(SeriesSolver<RT, PT, ST>& solver,
                                                  int numLoops,
                                                  int targetDeg,
                                                  const ST& feynmanD,
                                                  const ST& shiftA,
                                                  const ST& shiftB)
    : solver_(solver),
      numLoops_(numLoops),
      targetDeg_(targetDeg),
      feynmanD_(feynmanD),
      fbiDelta_(computeFBIDelta()),
      shiftA_(shiftA),
      shiftB_(shiftB),
      gamma_(computeGamma()),
      shiftedU_(buildShiftedU()),
      uPowerSeriesCache_(targetDeg),
      uPowerSeriesCached_(false) {
    if (numLoops_ <= 0) {
        throw std::invalid_argument("numLoops must be positive");
    }
    if (targetDeg_ < 0) {
        throw std::invalid_argument("targetDeg must be non-negative");
    }
}

template<typename RT, typename PT, typename ST>
ST IntegrandExpander<RT, PT, ST>::computeFBIDelta() const {
    return ST(numLoops_) * feynmanD_ / ST(2);
}

template<typename RT, typename PT, typename ST>
ST IntegrandExpander<RT, PT, ST>::computeGamma() const {
    return -(ST(numLoops_ + 1) * feynmanD_) / ST(2);
}

template<typename RT, typename PT, typename ST>
PT IntegrandExpander<RT, PT, ST>::makeConstantPoly(const ST& c) {
    PT p;
    p.addMonomial(c, Power(0, 0));
    return p;
}

template<typename RT, typename PT, typename ST>
PT IntegrandExpander<RT, PT, ST>::makeMonomialPoly(const ST& c, int xPow, int yPow) {
    PT p;
    p.addMonomial(c, Power(xPow, yPow));
    return p;
}

template<typename RT, typename PT, typename ST>
PT IntegrandExpander<RT, PT, ST>::multiplyPoly(const PT& a, const PT& b) {
    PT result;
    for (const auto& [pa, ca] : a) {
        for (const auto& [pb, cb] : b) {
            result.addMonomial(ca * cb, Power(pa.x_power + pb.x_power, pa.y_power + pb.y_power));
        }
    }
    return result;
}

template<typename RT, typename PT, typename ST>
PT IntegrandExpander<RT, PT, ST>::powPoly(PT base, int exp) {
    if (exp < 0) {
        throw std::invalid_argument("Polynomial negative exponent is not supported");
    }
    PT result = makeConstantPoly(ST(1));
    while (exp > 0) {
        if (exp & 1) {
            result = multiplyPoly(result, base);
        }
        exp >>= 1;
        if (exp > 0) {
            base = multiplyPoly(base, base);
        }
    }
    return result;
}

template<typename RT, typename PT, typename ST>
PT IntegrandExpander<RT, PT, ST>::applyShift(const PT& xrYrPoly) const {
    // Shift from (Xr, Yr) to (Xs, Ys):
    // Xr -> Xs + a, Yr -> Ys + b.
    const PT xVar = makeMonomialPoly(ST(1), 1, 0);
    const PT yVar = makeMonomialPoly(ST(1), 0, 1);
    const PT xShift = xVar + makeConstantPoly(shiftA_);
    const PT yShift = yVar + makeConstantPoly(shiftB_);

    PT result;
    for (const auto& [power, coeff] : xrYrPoly) {
        PT term = powPoly(xShift, power.x_power);
        term = multiplyPoly(term, powPoly(yShift, power.y_power));
        term *= coeff;
        result += term;
    }
    return result;
}

template<typename RT, typename PT, typename ST>
PT IntegrandExpander<RT, PT, ST>::buildShiftedU() const {
    // Local variables
    const PT xVar = makeMonomialPoly(ST(1), 1, 0);
    const PT yVar = makeMonomialPoly(ST(1), 0, 1);

    // Shifted local coordinates
    const PT Xs = xVar + makeConstantPoly(shiftA_);
    const PT Y0 = yVar + makeConstantPoly(shiftB_);

    // Mapping from simplex variables to unit square variables then shifted
    const PT one = makeConstantPoly(ST(1));
    const PT oneMinusXs = one + (Xs * ST(-1));
    const PT Ys = multiplyPoly(oneMinusXs, Y0);
    const PT Zs = one + (Xs * ST(-1)) + (Ys * ST(-1));

    // U = X*Y + Y*Z + Z*X
    PT Uxy = multiplyPoly(Xs, Ys);
    PT Uyz = multiplyPoly(Ys, Zs);
    PT Uzx = multiplyPoly(Zs, Xs);
    return Uxy + Uyz + Uzx;
}

template<typename RT, typename PT, typename ST>
PT IntegrandExpander<RT, PT, ST>::buildFIPolynomial(const std::vector<int>& nu) const {
    const auto& family = solver_.getFamily();
    const int numBranch = family.getNumBranch();
    const int numProps = family.getNumProps();
    if (numBranch != 3) {
        throw std::invalid_argument("buildFIPolynomial currently expects exactly 3 branches");
    }
    if (static_cast<int>(nu.size()) != numProps) {
        throw std::invalid_argument("nu size must equal number of propagators");
    }
    const auto& branchIndices = family.getBranchIndices();
    if (static_cast<int>(branchIndices.size()) != numProps) {
        throw std::runtime_error("branchIndices size mismatch");
    }

    std::vector<int> nuBranch(numBranch, 0);
    for (int i = 0; i < numProps; ++i) {
        int b = branchIndices[i];
        if (b < 0 || b >= numBranch) {
            throw std::runtime_error("invalid branch index in family");
        }
        nuBranch[b] += nu[i];
    }

    const int nuX = nuBranch[0];
    const int nuY = nuBranch[1];
    const int nuZ = nuBranch[2];

    const int ex = nuX - 1;
    const int ey = nuY - 1;
    const int ez = nuZ - 1;
    if (ex < 0 || ey < 0 || ez < 0) {
        throw std::invalid_argument("nuX, nuY, nuZ must be >= 1 for polynomial construction");
    }

    // Build FI polynomial first in (Xr, Yr):
    // X0 = Xr, Y0 = (1-Xr)Yr, Z0 = 1-Xr-(1-Xr)Yr.
    const PT xVar = makeMonomialPoly(ST(1), 1, 0);  // Xr
    const PT yVar = makeMonomialPoly(ST(1), 0, 1);  // Yr
    const PT one = makeConstantPoly(ST(1));
    const PT X0 = xVar;
    const PT Y0 = multiplyPoly(one + (xVar * ST(-1)), yVar);
    const PT Z0 = one + (xVar * ST(-1)) + (Y0 * ST(-1));

    // Jacobian from simplex -> square mapping
    PT J = one + (xVar * ST(-1));

    PT W = powPoly(X0, ex);
    W = multiplyPoly(W, powPoly(Y0, ey));
    W = multiplyPoly(W, powPoly(Z0, ez));

    // Integer U power in FI integrand: U^nu with nu = sum_i nu_i.
    int nuSum = 0;
    for (int v : nu) nuSum += v;
    if (nuSum < 0) {
        throw std::invalid_argument("sum(nu) must be non-negative");
    }

    PT baseXr = multiplyPoly(J, W);
    if (nuSum == 0) {
        return applyShift(baseXr);
    }
    // Build U in (Xr,Yr), then multiply U^nu, then shift once.
    const PT Uxy = multiplyPoly(X0, Y0);
    const PT Uyz = multiplyPoly(Y0, Z0);
    const PT Uzx = multiplyPoly(Z0, X0);
    PT uXr = Uxy + Uyz + Uzx;
    PT baseWithUXr = multiplyPoly(baseXr, powPoly(uXr, nuSum));
    return applyShift(baseWithUXr);
}

template<typename RT, typename PT, typename ST>
void IntegrandExpander<RT, PT, ST>::multiplySeries(Series<ST>& result,
                                                   const Series<ST>& a,
                                                   const Series<ST>& b) {
    if (result.getDeg() != a.getDeg() || result.getDeg() != b.getDeg()) {
        throw std::invalid_argument("Series degrees must match in multiplySeries");
    }
    result.setZero();
    const int deg = result.getDeg();
    for (int d = 0; d <= deg; ++d) {
        for (int p = 0; p <= d; ++p) {
            int q = d - p;
            ST sum = ST(0);
            for (int i = 0; i <= p; ++i) {
                for (int j = 0; j <= q; ++j) {
                    sum += a.getCoeff(i, j) * b.getCoeff(p - i, q - j);
                }
            }
            result.setCoeff(p, q, sum);
        }
    }
}

template<typename RT, typename PT, typename ST>
Series<ST> IntegrandExpander<RT, PT, ST>::expandUPower() const {
    // Solve for F = U^gamma:
    // U * dF/dX = gamma * (dU/dX) * F
    // U * dF/dY = gamma * (dU/dY) * F
    // with F(0,0)=1
    Series<ST> F(targetDeg_);
    F.setZero();
    F.setCoeff(0, 0, ST(1));

    const PT dUx = shiftedU_.derivativeX();
    const PT dUy = shiftedU_.derivativeY();
    const ST U00 = shiftedU_.getCoeff(0, 0);
    if (U00 == ST(0)) {
        throw std::runtime_error("U(0,0)=0, cannot use PDE recurrence for U^gamma");
    }

    for (int deg = 1; deg <= targetDeg_; ++deg) {
        // p > 0, use X equation.
        for (int p = 1; p <= deg; ++p) {
            int q = deg - p;

            ST rhs = gamma_ * SeriesSolver<RT, PT, ST>::polySeriesCoeff(dUx, F, p - 1, q);
            ST known = ST(0);

            for (const auto& [power, coeff] : shiftedU_) {
                const int a = power.x_power;
                const int b = power.y_power;
                if (a == 0 && b == 0) continue;
                if (p >= a && q >= b) {
                    int px = p - a;
                    int qy = q - b;
                    if (px > 0) {
                        known += coeff * ST(px) * F.getCoeff(px, qy);
                    }
                }
            }

            F.setCoeff(p, q, (rhs - known) / (U00 * ST(p)));
        }

        // p = 0, use Y equation.
        int q = deg;
        ST rhs = gamma_ * SeriesSolver<RT, PT, ST>::polySeriesCoeff(dUy, F, 0, q - 1);
        ST known = ST(0);

        for (const auto& [power, coeff] : shiftedU_) {
            const int a = power.x_power;
            const int b = power.y_power;
            if (a == 0 && b == 0) continue;
            if (a == 0 && q >= b) {
                int qy = q - b;
                if (qy > 0) {
                    known += coeff * ST(qy) * F.getCoeff(0, qy);
                }
            }
        }

        F.setCoeff(0, q, (rhs - known) / (U00 * ST(q)));
    }

    return F;
}

template<typename RT, typename PT, typename ST>
const Series<ST>& IntegrandExpander<RT, PT, ST>::getUPowerSeries() const {
    if (!uPowerSeriesCached_) {
        uPowerSeriesCache_ = expandUPower();
        uPowerSeriesCached_ = true;
    }
    return uPowerSeriesCache_;
}

template<typename RT, typename PT, typename ST>
void IntegrandExpander<RT, PT, ST>::clearCache() {
    uPowerSeriesCache_.setZero();
    uPowerSeriesCached_ = false;
}

template<typename RT, typename PT, typename ST>
Series<ST> IntegrandExpander<RT, PT, ST>::getFI2DSeries(const std::vector<int>& nu) const {
    PT poly = buildFIPolynomial(nu);
    const Series<ST>& fbi = solver_.getFBISeries(nu, fbiDelta_);

    Series<ST> tmp(targetDeg_);
    Series<ST>::mulPoly(tmp, fbi, poly);

    Series<ST> result(targetDeg_);
    multiplySeries(result, getUPowerSeries(), tmp);
    return result;
}
