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
      shiftedU_(buildShiftedU()) {
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

// ============================================================================
// 多项式辅助函数
// ============================================================================

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

// ============================================================================
// 坐标平移
// ============================================================================

template<typename RT, typename PT, typename ST>
PT IntegrandExpander<RT, PT, ST>::applyShift(const PT& xrYrPoly) const {
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

// ============================================================================
// 构建平移后的 U 多项式
// ============================================================================

template<typename RT, typename PT, typename ST>
PT IntegrandExpander<RT, PT, ST>::buildShiftedU() const {
    const PT xVar = makeMonomialPoly(ST(1), 1, 0);
    const PT yVar = makeMonomialPoly(ST(1), 0, 1);

    const PT Xs = xVar + makeConstantPoly(shiftA_);
    const PT Y0 = yVar + makeConstantPoly(shiftB_);

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

// ============================================================================
// 构建 FI 多项式: J · W（不含 U^{nuTot}）
// ============================================================================

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

    std::vector<int> nuBranch(numBranch, 0);
    for (int i = 0; i < numProps; ++i) {
        nuBranch[branchIndices[i]] += nu[i];
    }

    const int ex = nuBranch[0] - 1;
    const int ey = nuBranch[1] - 1;
    const int ez = nuBranch[2] - 1;
    if (ex < 0 || ey < 0 || ez < 0) {
        throw std::invalid_argument("nuX, nuY, nuZ must be >= 1");
    }

    const PT xVar = makeMonomialPoly(ST(1), 1, 0);
    const PT yVar = makeMonomialPoly(ST(1), 0, 1);
    const PT one = makeConstantPoly(ST(1));
    const PT X0 = xVar;
    const PT Y0 = multiplyPoly(one + (xVar * ST(-1)), yVar);
    const PT Z0 = one + (xVar * ST(-1)) + (Y0 * ST(-1));
    PT J = one + (xVar * ST(-1));  // Jacobian

    PT W = powPoly(X0, ex);
    W = multiplyPoly(W, powPoly(Y0, ey));
    W = multiplyPoly(W, powPoly(Z0, ez));

    // J · W（不含 U）
    return applyShift(multiplyPoly(J, W));
}

// ============================================================================
// FI 被积函数的二维级数
// ============================================================================

template<typename RT, typename PT, typename ST>
Series<ST> IntegrandExpander<RT, PT, ST>::getIntegrand2DSeries(const std::vector<int>& nu) const {
    // FI = J·W · Ĩ（所有 U 已被吸收进 Ĩ）
    const Series<ST>& fbi = solver_.getFBISeries(nu, fbiDelta_, targetDeg_);
    PT poly = buildFIPolynomial(nu);
    Series<ST> result(targetDeg_);
    Series<ST>::mulPoly(result, fbi, poly);
    return result;
}
