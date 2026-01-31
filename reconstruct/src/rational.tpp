#include <sstream>

// ===== Power类实现 =====

bool Power::operator<(const Power& other) const {
    // 先比较总度数 i+j
    int this_degree = x_power + y_power;
    int other_degree = other.x_power + other.y_power;
    if (this_degree != other_degree) {
        return this_degree < other_degree;
    }
    // 度数相同时，比较 x_power (i)
    return x_power < other.x_power;
}

std::string Power::toString() const {
    if (isConstant()) {
        return "1";
    }
    
    std::ostringstream oss;
    bool first = true;
    
    if (x_power != 0) {
        if (x_power == 1) {
            oss << "X";
        } else {
            oss << "X^" << x_power;
        }
        first = false;
    }
    
    if (y_power != 0) {
        if (!first) {
            oss << "*";
        }
        if (y_power == 1) {
            oss << "Y";
        } else {
            oss << "Y^" << y_power;
        }
    }
    
    return oss.str();
}

// ===== Polynomial模板类实现 =====

template<typename T>
void Polynomial<T>::addMonomial(const T& coeff, const Power& power) {
    auto it = poly_.find(power);
    if (it != poly_.end()) {
        it->second += coeff;  // 合并同类项
        // 如果系数为0，删除该项
        if (it->second == T{0}) {
            poly_.erase(it);
            // 重新计算参数
            updateDegreeAndCount();
        }
    } else if (coeff != T{0}) {
        // 非零系数才添加
        poly_[power] = coeff;
        // 更新单项式数量
        numOfMono_ = poly_.size();
        // 更新度：需要重新计算因为unordered_map无序
        int current_deg = power.x_power + power.y_power;
        if (deg_ < current_deg) {
            deg_ = current_deg;
        }
    }
}

// 辅助函数：重新计算度和单项式数量
template<typename T>
void Polynomial<T>::updateDegreeAndCount() {
    numOfMono_ = poly_.size();
    if (poly_.empty()) {
        deg_ = -1;
    } else {
        // 遍历所有项找最大度
        deg_ = 0;
        for (const auto& term : poly_) {
            int current_deg = term.first.x_power + term.first.y_power;
            if (current_deg > deg_) {
                deg_ = current_deg;
            }
        }
    }
}


template<typename T>
bool Polynomial<T>::isEmpty() const {
    return poly_.empty();
}

template<typename T>
bool Polynomial<T>::isConstant() const {
    return poly_.size() <= 1 && 
           (poly_.empty() || poly_.begin()->first.isConstant());
}

template<typename T>
size_t Polynomial<T>::size() const {
    return poly_.size();
}

template<typename T>
std::string Polynomial<T>::toString() const {
    if (isEmpty()) {
        return "0";
    }
    
    std::ostringstream oss;
    bool first = true;
    
    for (const auto& [power, coeff] : poly_) {
        if (!first) {
            oss << " + ";
        }
        
        if (power.isConstant()) {
            oss << coeff;
        } else {
            if (coeff == T{1}) {
                oss << power.toString();
            } else {
                oss << coeff << "*" << power.toString();
            }
        }
        
        first = false;
    }
    
    return oss.str();
}

template<typename T>
T Polynomial<T>::evaluate(const T& X, const T& Y) const {
    if (isEmpty()) {
        return T{0};
    }
    
    T result{0};
    for (const auto& [power, coeff] : poly_) {
        // 计算 X^x_power
        T x_term{1};
        for (int i = 0; i < power.x_power; ++i) {
            x_term *= X;
        }
        
        // 计算 Y^y_power
        T y_term{1};
        for (int i = 0; i < power.y_power; ++i) {
            y_term *= Y;
        }
        
        // 累加 coeff * X^x_power * Y^y_power
        result += coeff * x_term * y_term;
    }
    
    return result;
}

template<typename T>
Polynomial<T> Polynomial<T>::derivativeX() const {
    Polynomial<T> result;
    
    // 对每一项求导：d/dX (c * X^i * Y^j) = i * c * X^(i-1) * Y^j
    for (const auto& [power, coeff] : poly_) {
        if (power.x_power > 0) {
            // 新系数 = i * 原系数
            T new_coeff = T(power.x_power) * coeff;
            // 新幂次：X的指数减1
            Power new_power(power.x_power - 1, power.y_power);
            result.addMonomial(new_coeff, new_power);
        }
        // 如果 x_power == 0，该项对X的导数为0，不添加
    }
    
    return result;
}

template<typename T>
Polynomial<T> Polynomial<T>::derivativeY() const {
    Polynomial<T> result;
    
    // 对每一项求导：d/dY (c * X^i * Y^j) = j * c * X^i * Y^(j-1)
    for (const auto& [power, coeff] : poly_) {
        if (power.y_power > 0) {
            // 新系数 = j * 原系数
            T new_coeff = T(power.y_power) * coeff;
            // 新幂次：Y的指数减1
            Power new_power(power.x_power, power.y_power - 1);
            result.addMonomial(new_coeff, new_power);
        }
        // 如果 y_power == 0，该项对Y的导数为0，不添加
    }
    
    return result;
}

template<typename T>
Polynomial<T> Polynomial<T>::operator+(const Polynomial<T>& other) const {
    Polynomial<T> result = *this;
    result += other;
    return result;
}

template<typename T>
Polynomial<T>& Polynomial<T>::operator+=(const Polynomial<T>& other) {
    // 将other的每一项添加到当前多项式
    for (const auto& [power, coeff] : other.poly_) {
        addMonomial(coeff, power);
    }
    return *this;
}

template<typename T>
Polynomial<T> Polynomial<T>::operator*(const T& scalar) const {
    Polynomial<T> result = *this;
    result *= scalar;
    return result;
}

template<typename T>
Polynomial<T>& Polynomial<T>::operator*=(const T& scalar) {
    if (scalar == T(0)) {
        // 乘以0结果为零多项式
        poly_.clear();
        deg_ = -1;
        numOfMono_ = 0;
        return *this;
    }
    
    // 将每一项的系数乘以标量
    for (auto& [power, coeff] : poly_) {
        coeff *= scalar;
    }
    
    return *this;
}

// ===== Rational模板类实现 =====

template<typename T>
Rational<T>::Rational() {
    // 分母默认为1
    denominator.addMonomial(T{1}, Power(0, 0));
}

template<typename T>
Rational<T>::Rational(const Polynomial<T>& num) 
    : numerator(num) {
    // 分母为1
    denominator.addMonomial(T{1}, Power(0, 0));
}

template<typename T>
Rational<T>::Rational(const Polynomial<T>& num, const Polynomial<T>& den) 
    : numerator(num), denominator(den) {
    // 检查分母不为零
    if (denominator.isEmpty()) {
        throw std::invalid_argument("Denominator cannot be zero polynomial");
    }
}

template<typename T>
bool Rational<T>::isEmpty() const {
    return numerator.isEmpty();
}

template<typename T>
bool Rational<T>::isConstant() const {
    return numerator.isConstant() && denominator.isConstant();
}

template<typename T>
std::string Rational<T>::toString() const {
    if (isEmpty()) {
        return "0";
    }
    
    std::string num_str = numerator.toString();
    std::string den_str = denominator.toString();
    
    // 如果分母是1，只返回分子
    if (den_str == "1") {
        return num_str;
    }
    
    // 否则返回分数形式
    return "(" + num_str + ") / (" + den_str + ")";
}

// ===== Rational类运算符 =====

/// 有理函数加法: (a/b) + (c/d) = (a*d + b*c) / (b*d)
template<typename T>
Rational<T> Rational<T>::operator+(const Rational<T>& other) const {
    // 计算新的分子: a*d + b*c
    Polynomial<T> newNumerator;
    
    // a*d
    for (const auto& [power, coeff] : numerator) {
        for (const auto& [power2, coeff2] : other.denominator) {
            Power newPower(power.x_power + power2.x_power, power.y_power + power2.y_power);
            newNumerator.addMonomial(coeff * coeff2, newPower);
        }
    }
    
    // + b*c
    for (const auto& [power, coeff] : denominator) {
        for (const auto& [power2, coeff2] : other.numerator) {
            Power newPower(power.x_power + power2.x_power, power.y_power + power2.y_power);
            newNumerator.addMonomial(coeff * coeff2, newPower);
        }
    }
    
    // 计算新的分母: b*d
    Polynomial<T> newDenominator;
    for (const auto& [power, coeff] : denominator) {
        for (const auto& [power2, coeff2] : other.denominator) {
            Power newPower(power.x_power + power2.x_power, power.y_power + power2.y_power);
            newDenominator.addMonomial(coeff * coeff2, newPower);
        }
    }
    
    return Rational<T>(newNumerator, newDenominator);
}

/// 有理函数乘法: (a/b) * (c/d) = (a*c) / (b*d)
template<typename T>
Rational<T> Rational<T>::operator*(const Rational<T>& other) const {
    // 计算新的分子: a*c
    Polynomial<T> newNumerator;
    for (const auto& [power, coeff] : numerator) {
        for (const auto& [power2, coeff2] : other.numerator) {
            Power newPower(power.x_power + power2.x_power, power.y_power + power2.y_power);
            newNumerator.addMonomial(coeff * coeff2, newPower);
        }
    }
    
    // 计算新的分母: b*d
    Polynomial<T> newDenominator;
    for (const auto& [power, coeff] : denominator) {
        for (const auto& [power2, coeff2] : other.denominator) {
            Power newPower(power.x_power + power2.x_power, power.y_power + power2.y_power);
            newDenominator.addMonomial(coeff * coeff2, newPower);
        }
    }
    
    return Rational<T>(newNumerator, newDenominator);
}

// // ===== 显式实例化 =====
// template class Polynomial<int>;
// template class Polynomial<long long>;
// template class Rational<int>;
// template class Rational<long long>;
