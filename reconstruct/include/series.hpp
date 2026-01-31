#pragma once

#include <vector>
#include <string>
#include <stdexcept>
#include "rational.hpp"

/// 两变量幂级数类，表示 sum coefficient_{ij}X^iY^j, i+j<=deg
/// 存储方式：一维数组 + 度数
/// 一维数组的存储顺序：(0,0), (1,0), (0,1), (2,0), (1,1), (0,2), ...., (deg,0),...,(0,deg)
/// @tparam T 系数类型（可以是整数、double、有理数等）
template<typename T>
class Series {
private:
    std::vector<T> coefficients_;  ///< 系数数组，按特定顺序存储
    int deg_;                      ///< 级数的度数


public:
    /// 将二维坐标(i,j)转换为一维索引
    /// @param i X的指数
    /// @param j Y的指数
    /// @return 一维索引
    inline int getIndex(int i, int j) const {
        // 存储顺序：(0,0), (1,0), (0,1), (2,0), (1,1), (0,2), ...., (deg,0),...,(0,deg)
        // 按照总度数 n = i + j 分组，每组内按照 i 递减, j 递增的顺序
        int deg = i + j;
        return (deg * (deg + 1)) / 2 + j;
    }
    /// 构造函数
    /// @param deg 级数的度数
    Series(int deg = 0);

    void setZero();


    /// 获取指定项的系数
    /// @param i X的指数
    /// @param j Y的指数
    /// @return 系数值
    T getCoeff(int i, int j) const;

    /// 设置指定项的系数
    /// @param i X的指数
    /// @param j Y的指数
    /// @param coeff 系数值
    void setCoeff(int i, int j, const T& coeff);

    /// 获取级数的度数
    /// @return 度数
    int getDeg() const;

    /// 获取系数向量的大小
    /// @return 系数向量大小
    int size() const;

    /// 级数加法
    /// @param other 另一个级数
    /// @return 相加后的级数
    Series<T> operator+(const Series<T>& other) const;

    /// 级数加法（原地操作）
    /// @param other 另一个级数
    /// @return 自身的引用
    Series<T>& operator+=(const Series<T>& other);

    /// 获取系数向量（只读）
    /// @return 系数向量的常量引用
    const std::vector<T>& getCoefficients() const;

    /// 级数与多项式乘法
    /// @param poly 多项式
    /// @return 相乘后的级数
    Series<T> operator*(const Polynomial<T>& poly) const;

    /// 级数与多项式乘法（原地操作）
    /// @param poly 多项式
    /// @return 自身的引用
    Series<T>& operator*=(const Polynomial<T>& poly);

    /// 级数与多项式除法
    /// @param poly 多项式
    /// @return 相除后的级数
    Series<T> operator/(const Polynomial<T>& poly) const;

    /// 级数与多项式除法（原地操作）
    /// @param poly 多项式
    /// @return 自身的引用
    Series<T>& operator/=(const Polynomial<T>& poly);

    /// 级数与有理函数乘法
    /// @param rational 有理函数
    /// @return 相乘后的级数
    Series<T> operator*(const Rational<T>& rational) const;

    /// 级数与有理函数乘法（原地操作）
    /// @param rational 有理函数
    /// @return 自身的引用
    Series<T>& operator*=(const Rational<T>& rational);

    /// 级数与有理函数除法
    /// @param rational 有理函数
    /// @return 相除后的级数
    Series<T> operator/(const Rational<T>& rational) const;

    /// 级数与有理函数除法（原地操作）
    /// @param rational 有理函数
    /// @return 自身的引用
    Series<T>& operator/=(const Rational<T>& rational);

    /// 级数与级数乘法
    /// @param other 另一个级数
    /// @return 相乘后的级数
    Series<T> operator*(const Series<T>& other) const;

    /// 级数与级数乘法（原地操作）
    /// @param other 另一个级数
    /// @return 自身的引用
    Series<T>& operator*=(const Series<T>& other);

    /// 打印级数为数学表达式
    /// @return 级数的字符串表示，格式为 sum a_{ij}X^iY^j
    std::string toString() const;

    /// 静态函数：级数与多项式乘法（使用预分配空间）
    /// @param result 结果级数（预分配）
    /// @param series 输入级数
    /// @param poly 多项式
    /// @throws std::invalid_argument 如果result的度数与series不匹配
    static void mulPoly(Series<T>& result, const Series<T>& series, const Polynomial<T>& poly);

    /// 静态函数：级数与多项式除法（使用预分配空间）
    /// @param result 结果级数（预分配）
    /// @param series 输入级数
    /// @param poly 多项式
    /// @throws std::invalid_argument 如果result的度数与series不匹配
    static void divPoly(Series<T>& result, const Series<T>& series, const Polynomial<T>& poly);

    /// 静态函数：级数与有理函数乘法（使用预分配空间）
    /// @param result 结果级数（预分配）
    /// @param series 输入级数
    /// @param rational 有理函数
    /// @throws std::invalid_argument 如果result的度数与series不匹配
    static void mulRat(Series<T>& result, const Series<T>& series, const Rational<T>& rational);

    /// 静态函数：级数与有理函数除法（使用预分配空间）
    /// @param result 结果级数（预分配）
    /// @param series 输入级数
    /// @param rational 有理函数
    /// @throws std::invalid_argument 如果result的度数与series不匹配
    static void divRat(Series<T>& result, const Series<T>& series, const Rational<T>& rational);

    /// 静态函数：级数与级数加法（使用预分配空间）
    /// @param result 结果级数（预分配）
    /// @param series1 第一个级数
    /// @param series2 第二个级数
    /// @throws std::invalid_argument 如果result的度数不匹配
    static void addSeries(Series<T>& result, const Series<T>& series1, const Series<T>& series2);
};

#include "../src/series.tpp"
