#pragma once

#include <vector>

template<typename T>
class Polynomial1D {
public:
    Polynomial1D();
    explicit Polynomial1D(int degree);
    explicit Polynomial1D(std::vector<T> coeffs);

    int degree() const;
    bool isZero() const;
    const std::vector<T>& coeffs() const;
    T coeff(int i) const;
    void setCoeff(int i, const T& value);
    void trim();
    void shiftByX();
    T eval(const T& x) const;

private:
    std::vector<T> coeffs_;
};

template<typename T>
Polynomial1D<T> operator+(const Polynomial1D<T>& lhs, const Polynomial1D<T>& rhs);

template<typename T>
Polynomial1D<T> operator-(const Polynomial1D<T>& lhs, const Polynomial1D<T>& rhs);

template<typename T>
Polynomial1D<T> operator*(const Polynomial1D<T>& lhs, const T& rhs);

template<typename T>
Polynomial1D<T> operator*(const T& lhs, const Polynomial1D<T>& rhs);

template<typename T>
std::vector<T> multiplyPolySeries(const Polynomial1D<T>& poly,
                                  const std::vector<T>& series,
                                  int degreeD);

template<typename T>
std::vector<T> divideSeriesByPoly(const std::vector<T>& numerator,
                                  const Polynomial1D<T>& denominator,
                                  int degreeD);

#include "../src/polynomial_1d.tpp"
