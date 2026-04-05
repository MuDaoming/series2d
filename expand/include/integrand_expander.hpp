#pragma once

#include <chrono>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>
#include "series_solver.hpp"

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

    ST gamma_;
    PT shiftedU_;

    mutable Series<ST> uPowerSeriesCache_;
    mutable bool uPowerSeriesCached_;

    // Redefinition 支持
    const Redefinition<PT, ST>* redef_ = nullptr;
    ST gammaRedef_;    // gamma_base = gamma - (L+1)*D_in/2
    mutable Series<ST> uPowerRedefCache_;
    mutable bool uPowerRedefCached_;

public:
    IntegrandExpander(SeriesSolver<RT, PT, ST>& solver,
                      int numLoops,
                      int targetDeg,
                      const ST& feynmanD,
                      const ST& shiftA,
                      const ST& shiftB);

    // Build FI integrand 2D series:
    // P(X,Y) * U(X,Y)^gamma * I_nu^Delta(X,Y),
    // where P already includes the integer power U^nu.
    Series<ST> getFI2DSeries(const std::vector<int>& nu) const;

    const PT& getShiftedU() const { return shiftedU_; }
    ST getGamma() const { return gamma_; }
    ST getFeynmanD() const { return feynmanD_; }
    ST getFBIDelta() const { return fbiDelta_; }
    const Series<ST>& getUPowerSeries() const;
    const Series<ST>& getUPowerRedefSeries() const;

    void setRedefinition(const Redefinition<PT, ST>* redef);
    void clearCache();

private:
    ST computeFBIDelta() const;
    ST computeGamma() const;

    PT buildShiftedU() const;
    PT buildFIPolynomial(const std::vector<int>& nu) const;
    PT applyShift(const PT& xrYrPoly) const;

    Series<ST> expandUPower() const;

    static PT makeConstantPoly(const ST& c);
    static PT makeMonomialPoly(const ST& c, int xPow, int yPow);
    static PT multiplyPoly(const PT& a, const PT& b);
    static PT powPoly(PT base, int exp);
    static void multiplySeries(Series<ST>& result, const Series<ST>& a, const Series<ST>& b);
};

#include "../src/integrand_expander.tpp"
