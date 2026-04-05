#pragma once

#include <stdexcept>
#include <vector>
#include "series_solver.hpp"

/**
 * IntegrandExpander - 被积函数展开器
 * 
 * 将重定义后的 FBI 级数 Ĩ 组合成 Feynman 积分的被积函数:
 *   FI = J · W · Ĩ
 * 其中 J 是 Jacobian, W = X^{e_X} · Y^{e_Y} · Z^{e_Z}
 */
template<typename RT, typename PT, typename ST>
class IntegrandExpander {
private:
    SeriesSolver<RT, PT, ST>& solver_;

    int numLoops_;
    int targetDeg_;
    ST feynmanD_;
    ST fbiDelta_;
    ST shiftA_;
    ST shiftB_;

    PT shiftedU_;

public:
    IntegrandExpander(SeriesSolver<RT, PT, ST>& solver,
                      int numLoops,
                      int targetDeg,
                      const ST& feynmanD,
                      const ST& shiftA,
                      const ST& shiftB);

    /// 构建 FI 被积函数的二维级数: FI = J·W · Ĩ_nu
    Series<ST> getFI2DSeries(const std::vector<int>& nu) const;

    const PT& getShiftedU() const { return shiftedU_; }
    ST getFeynmanD() const { return feynmanD_; }
    ST getFBIDelta() const { return fbiDelta_; }

private:
    ST computeFBIDelta() const;

    PT buildShiftedU() const;
    PT buildFIPolynomial(const std::vector<int>& nu) const;
    PT applyShift(const PT& xrYrPoly) const;

    static PT makeConstantPoly(const ST& c);
    static PT makeMonomialPoly(const ST& c, int xPow, int yPow);
    static PT multiplyPoly(const PT& a, const PT& b);
    static PT powPoly(PT base, int exp);
};

#include "../src/integrand_expander.tpp"
