#pragma once

#include <stdexcept>

template<typename T>
Polynomial1D<T>::Polynomial1D() = default;

template<typename T>
Polynomial1D<T>::Polynomial1D(int degree) : coeffs_(degree + 1, T(0)) {}

template<typename T>
Polynomial1D<T>::Polynomial1D(std::vector<T> coeffs) : coeffs_(std::move(coeffs)) {
    trim();
}

template<typename T>
int Polynomial1D<T>::degree() const {
    for (int i = static_cast<int>(coeffs_.size()) - 1; i >= 0; --i) {
        if (coeffs_[i] != T(0)) return i;
    }
    return -1;
}

template<typename T>
bool Polynomial1D<T>::isZero() const {
    return degree() < 0;
}

template<typename T>
const std::vector<T>& Polynomial1D<T>::coeffs() const {
    return coeffs_;
}

template<typename T>
T Polynomial1D<T>::coeff(int i) const {
    if (i < 0 || i >= static_cast<int>(coeffs_.size())) return T(0);
    return coeffs_[i];
}

template<typename T>
void Polynomial1D<T>::setCoeff(int i, const T& value) {
    if (i < 0) throw std::runtime_error("Negative polynomial index");
    if (i >= static_cast<int>(coeffs_.size())) coeffs_.resize(i + 1, T(0));
    coeffs_[i] = value;
}

template<typename T>
void Polynomial1D<T>::trim() {
    while (!coeffs_.empty() && coeffs_.back() == T(0)) coeffs_.pop_back();
}

template<typename T>
void Polynomial1D<T>::shiftByX() {
    if (coeffs_.empty()) return;
    coeffs_.insert(coeffs_.begin(), T(0));
}

template<typename T>
T Polynomial1D<T>::eval(const T& x) const {
    T result(0);
    for (int i = static_cast<int>(coeffs_.size()) - 1; i >= 0; --i) {
        result = result * x + coeffs_[i];
    }
    return result;
}

template<typename T>
Polynomial1D<T> operator+(const Polynomial1D<T>& lhs, const Polynomial1D<T>& rhs) {
    const int n = std::max(lhs.degree(), rhs.degree()) + 1;
    std::vector<T> out(n, T(0));
    for (int i = 0; i < n; ++i) out[i] = lhs.coeff(i) + rhs.coeff(i);
    return Polynomial1D<T>(std::move(out));
}

template<typename T>
Polynomial1D<T> operator-(const Polynomial1D<T>& lhs, const Polynomial1D<T>& rhs) {
    const int n = std::max(lhs.degree(), rhs.degree()) + 1;
    std::vector<T> out(n, T(0));
    for (int i = 0; i < n; ++i) out[i] = lhs.coeff(i) - rhs.coeff(i);
    return Polynomial1D<T>(std::move(out));
}

template<typename T>
Polynomial1D<T> operator*(const Polynomial1D<T>& lhs, const T& rhs) {
    std::vector<T> out = lhs.coeffs();
    for (auto& x : out) x *= rhs;
    return Polynomial1D<T>(std::move(out));
}

template<typename T>
Polynomial1D<T> operator*(const T& lhs, const Polynomial1D<T>& rhs) {
    return rhs * lhs;
}

template<typename T>
std::vector<T> multiplyPolySeries(const Polynomial1D<T>& poly,
                                  const std::vector<T>& series,
                                  int degreeD) {
    std::vector<T> out(degreeD + 1, T(0));
    for (int i = 0; i <= degreeD; ++i) {
        const T pc = poly.coeff(i);
        if (pc == T(0)) continue;
        for (int j = 0; j + i <= degreeD && j < static_cast<int>(series.size()); ++j) {
            out[i + j] += pc * series[j];
        }
    }
    return out;
}

template<typename T>
std::vector<T> divideSeriesByPoly(const std::vector<T>& numerator,
                                  const Polynomial1D<T>& denominator,
                                  int degreeD) {
    const T d0 = denominator.coeff(0);
    if (d0 == T(0)) {
        throw std::runtime_error("Series division requires nonzero denominator constant term");
    }
    std::vector<T> out(degreeD + 1, T(0));
    for (int n = 0; n <= degreeD; ++n) {
        T rhs = n < static_cast<int>(numerator.size()) ? numerator[n] : T(0);
        for (int i = 1; i <= n; ++i) {
            rhs -= denominator.coeff(i) * out[n - i];
        }
        out[n] = rhs / d0;
    }
    return out;
}
