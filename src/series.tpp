#include <vector>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <string>


template<typename T>
Series<T>::Series(int deg) : deg_(deg) {
    int totalSize = (deg + 1) * (deg + 2) / 2;
    coefficients_.resize(totalSize, T(0));
}


template<typename T>
void Series<T>::setZero() {
    std::fill(coefficients_.begin(), coefficients_.end(), T(0));
}

template<typename T>
T Series<T>::getCoeff(int i, int j) const {
    if (i < 0 || j < 0 || i + j > deg_) {
        return T(0);  // 超出范围的系数为0
    }
    return coefficients_[getIndex(i, j)];
}

template<typename T>
void Series<T>::setCoeff(int i, int j, const T& coeff) {
    if (i < 0 || j < 0 || i + j > deg_) {
        throw std::out_of_range("Index out of range");
    }
    coefficients_[getIndex(i, j)] = coeff;
}

template<typename T>
int Series<T>::getDeg() const {
    return deg_;
}

template<typename T>
int Series<T>::size() const {
    return coefficients_.size();
}

template<typename T>
const std::vector<T>& Series<T>::getCoefficients() const {
    return coefficients_;
}

template<typename T>
std::string Series<T>::toString() const {
    std::ostringstream oss;
    bool first = true;
    
    // 按照总度数递增的顺序遍历所有项
    for (int totalDeg = 0; totalDeg <= deg_; ++totalDeg) {
        for (int i = totalDeg; i >= 0; --i) {
            int j = totalDeg - i;
            
            T coeff = getCoeff(i, j);
            
            // 跳过零系数
            if (coeff == T(0)) continue;
            
            // 处理符号
            if (!first) {
                if (coeff > T(0)) {
                    oss << " + ";
                } else {
                    oss << " - ";
                    coeff = -coeff;  // 取绝对值
                }
            } else {
                if (coeff < T(0)) {
                    oss << "-";
                    coeff = -coeff;
                }
                first = false;
            }
            
            // 处理系数
            bool needCoeff = true;
            if (coeff == T(1)) {
                // 系数为1时，如果是常数项需要显示，否则可以省略
                if (i == 0 && j == 0) {
                    oss << "1";
                    needCoeff = false;
                } else {
                    needCoeff = false;  // 省略系数1
                }
            } else {
                oss << coeff;
                needCoeff = false;
            }
            
            // 处理X的幂次
            if (i > 0) {
                if (needCoeff) oss << "*";
                if (i == 1) {
                    oss << "X";
                } else {
                    oss << "X^" << i;
                }
            }
            
            // 处理Y的幂次
            if (j > 0) {
                if (needCoeff || i > 0) oss << "*";
                if (j == 1) {
                    oss << "Y";
                } else {
                    oss << "Y^" << j;
                }
            }
            
            // 如果系数为1且是常数项，已经处理过了
            if (i == 0 && j == 0 && needCoeff) {
                oss << coeff;
            }
        }
    }
    
    // 如果所有系数都为零
    if (first) {
        return "0";
    }
    
    return oss.str();
}

template<typename T>
void Series<T>::addSeries(Series<T>& result, const Series<T>& series1, const Series<T>& series2) {
    if (result.getDeg() != series1.getDeg() || result.getDeg() != series2.getDeg()) {
        throw std::invalid_argument("Series degrees must be equal for addition");
    }

    // 对一维数组进行加法即可
    for (size_t idx = 0; idx < result.size(); ++idx) {
        result.coefficients_[idx] = series1.coefficients_[idx] + series2.coefficients_[idx];
    }
}

template<typename T>
void Series<T>::mulPoly(Series<T>& result, const Series<T>& series, const Polynomial<T>& poly) {
    if (result.getDeg() != series.getDeg()) {
        throw std::invalid_argument("Result series degree must match input series degree for multiplication");
    }
    int poly_deg = poly.getDegree();
    int series_deg = series.getDeg();

    for (int deg = series_deg; deg >= 0; --deg) {
        for (int m = deg; m >= 0; --m) {
            int n = deg - m;
            int l = deg * (deg + 1) / 2 + m;

            T sum = T(0);
            
            int i_min = (n - poly_deg > 0) ? (n - poly_deg) : 0;
            int i_max = (series_deg < n) ? series_deg : n;
            int j_min = (m - poly_deg > 0) ? (m - poly_deg) : 0;
            int j_max = (series_deg < m) ? series_deg : m;

            int numOfTerms = (i_max - i_min + 1) * (j_max - j_min + 1);

            
            int numOfMono = poly.getNumOfMonomials();

            if (numOfTerms < numOfMono)
            {   
                // 循环的范围为 i,j
                for (int i = i_min; i <= i_max; ++i) {
                    // j 的范围：Max[0, m - poly.getDegree()], Min[this->getDeg(), m]
                    for (int j = j_min; j <= j_max; ++j) {
                        T series_coeff = series.getCoeff(i, j);
                        T poly_coeff = poly.getCoeff(n - i, m - j);
                        sum += series_coeff * poly_coeff;
                
                    }
                }
            }
            else
            {
                // 循环的范围为多项式的单项式
                for (const auto& [power, coeff] : poly) {
                    int i = power.x_power;
                    int j = power.y_power;
                    if (i <= n && j <= m) {
                        T series_coeff = series.getCoeff(n - i, m - j);
                        sum += series_coeff * coeff;
                    }
                }
            }

            result.coefficients_[l] = sum;
        }
    }
}

template<typename T>
void Series<T>::divPoly(Series<T>& result, const Series<T>& series, const Polynomial<T>& poly) {
    // 验证度数一致性
    if (result.getDeg() != series.getDeg()) {
        throw std::invalid_argument("Result and series must have same degree");
    }
    
    // 检查多项式常数项是否为零
    T b00 = poly.getCoeff(0, 0);
    if (b00 == T(0)) {
        throw std::invalid_argument("Error: denominator has zero constant term!");
    }
    
    // 初始化常数项
    result.setCoeff(0, 0, series.getCoeff(0, 0) / b00);

    int poly_deg = poly.getDegree();
    int series_deg = series.getDeg();

    // 按照一维序列从前往后计算（与.wl算法对应）
    for (int deg = 1; deg <= series_deg; ++deg) {
        for (int m = 0; m <= deg; ++m) {
            int n = deg - m;
            int l = deg * (deg + 1) / 2 + m;
            
            // 计算 result[n, m] = (a1[n,m] - Sum[a2[i,j] * result[n-i, m-j]]) / b00
            // 其中 (i,j) != (0,0)
            T numerator_term = series.getCoeff(n, m);
            T sum = T(0);
            
            // 内部循环保持与 .wl 算法相同的起止点
            int i_min = (n - series_deg > 0) ? (n - series_deg) : 0;
            int i_max = (poly_deg < n) ? poly_deg : n;
            int j_min = (m - series_deg > 0) ? (m - series_deg) : 0;
            int j_max = (poly_deg < m) ? poly_deg : m;

            int numOfTerms = (i_max - i_min + 1) * (j_max - j_min + 1);
            int numOfMono = poly.getNumOfMonomials();

            if (numOfTerms < numOfMono)
            {   
                // 循环的范围为 i,j
                for (int i = i_min; i <= i_max; ++i) {
                    for (int j = j_min; j <= j_max; ++j) {
                        // 跳过 (0,0) 项
                        if (i == 0 && j == 0) continue;
                        
                        T poly_coeff = poly.getCoeff(i, j);
                        T result_coeff = result.getCoeff(n - i, m - j);
                        sum += poly_coeff * result_coeff;
                    }
                }
            }
            else
            {
                // 循环的范围为多项式的单项式
                for (const auto& [power, coeff] : poly) {
                    int i = power.x_power;
                    int j = power.y_power;
                    // 跳过 (0,0) 项
                    if (i == 0 && j == 0) continue;
                    
                    if (i <= n && j <= m) {
                        T result_coeff = result.getCoeff(n - i, m - j);
                        sum += coeff * result_coeff;
                    }
                }
            }

            result.coefficients_[l] = (numerator_term - sum) / b00;
        }
    }
}

template<typename T>
Series<T> Series<T>::operator+(const Series<T>& other) const {
    Series<T> result(deg_);
    Series<T>::addSeries(result, *this, other);
    return result;
}


template<typename T>
Series<T>& Series<T>::operator+=(const Series<T>& other) {
    Series<T>::addSeries(*this, *this, other);
    return *this;
}

template<typename T>
Series<T> Series<T>::operator*(const Polynomial<T>& poly) const {
    Series<T> result(deg_);
    Series<T>::mulPoly(result, *this, poly);
    return result;
}

template<typename T>
Series<T>& Series<T>::operator*=(const Polynomial<T>& poly) {
    Series<T>::mulPoly(*this, *this, poly);
    return *this;
}

template<typename T>
Series<T> Series<T>::operator/(const Polynomial<T>& poly) const {
    Series<T> result(deg_);
    Series<T>::divPoly(result, *this, poly);
    return result;
}

template<typename T>
Series<T>& Series<T>::operator/=(const Polynomial<T>& poly) {
    Series<T>::divPoly(*this, *this, poly);
    return *this;
}

template<typename T>
void Series<T>::divRat(Series<T>& result, const Series<T>& series, const Rational<T>& rational) {
    if (result.getDeg() != series.getDeg()) {
        throw std::invalid_argument("Result and series must have same degree");
    }

    Series<T>::mulPoly(result, series, rational.getNumerator());
    Series<T>::divPoly(result, result, rational.getDenominator());

}

template<typename T>
void Series<T>::mulRat(Series<T>& result, const Series<T>& series, const Rational<T>& rational) {
    if (result.getDeg() != series.getDeg()) {
        throw std::invalid_argument("Result and series must have same degree");
    }

    Series<T>::mulPoly(result, series, rational.numerator);
    Series<T>::divPoly(result, result, rational.denominator);

}





