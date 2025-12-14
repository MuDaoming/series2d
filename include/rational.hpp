#pragma once

#include <map>
#include <unordered_map>
#include <vector>
#include <string>
#include <stdexcept>

/// 幂次对，表示 X^x_power * Y^y_power
struct Power {
    int x_power;    ///< X的指数
    int y_power;    ///< Y的指数
    
    /// 默认构造函数（常数项）
    Power(int x = 0, int y = 0) : x_power(x), y_power(y) {}
    
    /// 字典序比较，用于std::map排序 (保留以备后用)
    bool operator<(const Power& other) const;
    
    /// 相等比较，用于std::unordered_map
    bool operator==(const Power& other) const {
        return x_power == other.x_power && y_power == other.y_power;
    }
    
    /// 转换为字符串（调试用）
    std::string toString() const;
    
    /// 检查是否为常数项
    bool isConstant() const {
        return x_power == 0 && y_power == 0;
    }
};

static inline uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

struct PowerHash {
    size_t operator()(const Power& p) const noexcept {
        uint64_t x = ((uint64_t)p.x_power << 32) ^ (uint64_t)p.y_power;
        return splitmix64(x);
    }
};


// /// Power的哈希函数，用于std::unordered_map
// struct PowerHash {
//     std::size_t operator()(const Power& p) const {
//         // 使用简单的哈希组合
//         return std::hash<int>()(p.x_power) ^ (std::hash<int>()(p.y_power) << 1);
//     }
// };

/// 多项式模板类
/// @tparam T 系数类型（可以是整数、有理数、有限域元素等）
template<typename T>
class Polynomial {
private:
    std::unordered_map<Power, T, PowerHash> poly_;  ///< 幂次 -> 系数映射，无序哈希表
    int deg_;                                       ///< 多项式的度，定义为max(i+j)，其中幂次为X^i*Y^j
    size_t numOfMono_;                              ///< 单项式的数量

public:
    /// 默认构造函数（零多项式）
    Polynomial() : deg_(-1), numOfMono_(0) {poly_.max_load_factor(0.5f); }
    
    /// 添加单项式
    /// @param coeff 系数
    /// @param power 幂次对
    void addMonomial(const T& coeff, const Power& power);

    /// 检查是否为零多项式
    bool isEmpty() const;
    
    /// 检查是否为常数
    bool isConstant() const;
    
    /// 获取项的数量
    size_t size() const;
    
    /// 获取多项式的度
    /// @return 多项式的度，零多项式返回-1
    int getDegree() const { return deg_; }
    
    /// 获取单项式数量
    /// @return 单项式的数量
    size_t getNumOfMonomials() const { return numOfMono_; }
    
    /// 转换为字符串（调试用）
    std::string toString() const;

    /// 获取多项式的迭代器（只读访问）
    auto begin() const -> decltype(poly_.begin()) {
        return poly_.begin();
    }
    
    auto end() const -> decltype(poly_.end()) {
        return poly_.end();
    }
    
    /// 获取指定幂次的系数
    T getCoeff(const Power& power) const {
        auto it = poly_.find(power);
        return (it != poly_.end()) ? it->second : T(0);
    }
    
    /// 获取指定幂次的系数 (通过x,y指数)
    T getCoeff(int x, int y) const {
        return getCoeff(Power(x, y));
    }

private:
    /// 辅助函数：重新计算度和单项式数量
    void updateDegreeAndCount();

};

/// 有理函数模板类
/// @tparam T 系数类型
template<typename T>
struct Rational {
    Polynomial<T> numerator;      ///< 分子
    Polynomial<T> denominator;    ///< 分母
    
    /// 默认构造函数（零有理函数）
    Rational();
    
    /// 单参数构造函数（多项式作为分子，分母为1）
    Rational(const Polynomial<T>& num);
    
    /// 完整构造函数
    Rational(const Polynomial<T>& num, const Polynomial<T>& den);

    /// 检查是否为空
    bool isEmpty() const;

    /// 检查是否为常数
    bool isConstant() const;
    
    /// 转换为字符串（调试用）
    std::string toString() const;
};

/// 矩阵类型别名
template<typename T>
using Matrix = std::vector<std::vector<Rational<T>>>;

#include "../src/rational.tpp"
