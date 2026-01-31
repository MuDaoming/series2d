#pragma once

#include <vector>
#include "rational.hpp"
#include "series.hpp"

/// 二维微分方程组求解器类
/// 用于求解形如：
/// d f / d X = AX f
/// d f / d Y = AY f
/// 其中 f 是未知函数向量，AX 和 AY 是有理函数系数矩阵
/// @tparam T 系数类型（可以是整数、有理数等）
template<typename T>
class DiffSystem {
private:
    const std::vector<std::vector<Rational<T>>>& AX_;  ///< X方向系数矩阵的引用
    const std::vector<std::vector<Rational<T>>>& AY_;  ///< Y方向系数矩阵的引用

public:
    /// 构造函数
    /// @param AX X方向系数矩阵
    /// @param AY Y方向系数矩阵
    DiffSystem(const std::vector<std::vector<Rational<T>>>& AX,
               const std::vector<std::vector<Rational<T>>>& AY);

    /// 求解微分方程组
    /// @param result 输出的解向量
    /// @param f0 初值条件向量
    /// @param deg 级数的度数
    void solve(std::vector<Series<T>>& result, const std::vector<T>& f0, int deg);

    /// 获取系统规模（方程个数）
    /// @return 方程个数
    int getSystemSize() const;

private:
    /// 计算耦合项
    /// @param g1 X方向耦合项的结果（输出参数）
    /// @param g2 Y方向耦合项的结果（输出参数）
    /// @param f 当前解向量
    /// @param i 当前方程索引
    void getSub(Series<T>& g1, Series<T>& g2, const std::vector<Series<T>>& f, int i);

    /// 求解标准形式的微分方程（输出参数版本）
    /// @param result 输出的结果级数
    /// @param R1 X方向的有理函数系数
    /// @param g1 X方向的右端项级数（会被原地修改）
    /// @param R2 Y方向的有理函数系数  
    /// @param g2 Y方向的右端项级数（会被原地修改）
    /// @param f0 初值条件
    void solveStandardDE(Series<T>& result, const Rational<T>& R1, Series<T>& g1,
                        const Rational<T>& R2, Series<T>& g2,
                        const T& f0);

    /// 验证系统的一致性
    /// @throws std::invalid_argument 如果系统不一致
    void validateSystem() const;
};

#include "../src/diffeq.tpp"
