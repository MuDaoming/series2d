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
            if (coeff >= T{0}) {
                oss << " + ";
            } else {
                oss << " - ";
            }
        } else if (coeff < T{0}) {
            oss << "-";
        }
        
        // 计算绝对值系数
        T abs_coeff = (coeff < T{0}) ? -coeff : coeff;
        
        if (power.isConstant()) {
            oss << abs_coeff;
        } else {
            if (abs_coeff == T{1}) {
                oss << power.toString();
            } else {
                oss << abs_coeff << "*" << power.toString();
            }
        }
        
        first = false;
    }
    
    return oss.str();
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

// // ===== 显式实例化 =====
// template class Polynomial<int>;
// template class Polynomial<long long>;
// template class Rational<int>;
// template class Rational<long long>;
