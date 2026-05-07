#pragma once

#include "rational.hpp"
#include "de_builder.hpp"

// 使用 #define 访问 Firefly 私有成员以便自定义质数
#define private public
#define protected public
#include "firefly/Reconstructor.hpp"
#undef protected
#undef private

#include <vector>
#include <memory>

// 前向声明
template<typename T>
class DEInterpolater;

namespace firefly {
    /// BlackBox 类：用于插值微分方程矩阵
    /// 在给定(X, Y)时，返回AX和AY矩阵的所有元素（扁平化）
    template<typename T>
    class DEBlackBox : public BlackBoxBase<DEBlackBox<T>> {
    public:
        /// 构造函数
        explicit DEBlackBox(const DEInterpolater<T>* interpolater)
            : interpolater_(interpolater) {}
        
        /// BlackBox 评估函数
        /// values[0] = X, values[1] = Y
        /// 返回AX和AY矩阵的所有元素（扁平化）
        template<typename FFIntTemp>
        std::vector<FFIntTemp> operator()(const std::vector<FFIntTemp>& values);
        
        /// 质数变化时调用
        inline void prime_changed() {
            // 质数已经自动在 FFInt::p 中更新
        }
        
    private:
        const DEInterpolater<T>* interpolater_;  ///< 指向插值器
    };
}

/// DEInterpolater 类：符号构建微分方程，返回有理函数形式的系数矩阵
/// 使用Firefly插值重构得到微分方程系数矩阵关于(X,Y)的有理函数表达式
/// @tparam T 有限域类型（如 FlintMod）
template<typename T>
class DEInterpolater {
public:
    /// 构造函数
    /// @param topS topS 矩阵，每个元素是 Polynomial<T>
    /// @param numProps 传播子数量
    /// @param numBranch Branch 数量
    /// @param delta delta 参数
    /// @param prime 有限域的素数模数
    DEInterpolater(const std::vector<std::vector<Polynomial<T>>>& topS,
                   const Polynomial<T>& UPoly,
                   int numProps,
                   int numBranch,
                   T delta,
                   uint64_t prime);
    
    /// 构建符号微分方程矩阵
    /// @param AX 输出：X方向微分方程矩阵（有理函数形式）
    /// @param AY 输出：Y方向微分方程矩阵（有理函数形式）
    void buildDEMatrices(std::vector<std::vector<Rational<T>>>& AX,
                        std::vector<std::vector<Rational<T>>>& AY);
    
    /// 获取总探测点数（最近一次调用的统计）
    uint32_t getTotalProbes() const { return total_probes_; }
    
    /// 在数值点创建DEBuilder（供BlackBox使用）
    DEBuilder<T> createDEBuilder(const T& X, const T& Y) const;

private:
    std::vector<std::vector<Polynomial<T>>> topS_;      ///< 多项式形式的topS矩阵
    Polynomial<T> UPoly_;                               ///< 符号形式的 U(X,Y)
    std::vector<std::vector<Polynomial<T>>> dRdX_;      ///< 符号形式的dR/dX矩阵
    std::vector<std::vector<Polynomial<T>>> dRdY_;      ///< 符号形式的dR/dY矩阵
    Polynomial<T> dUdXPoly_;                            ///< 符号形式的 dU/dX
    Polynomial<T> dUdYPoly_;                            ///< 符号形式的 dU/dY
    int numProps_;                                       ///< 传播子数量
    int numBranch_;                                      ///< Branch数量
    T delta_;                                            ///< delta参数
    uint64_t prime_;                                     ///< 有限域素数
    int numMasterFBI_;                                   ///< master FBI数量
    uint32_t total_probes_;                             ///< 总探测点数（最近一次）
    
    /// 在数值点求值topS
    std::vector<std::vector<T>> evaluateTopS(const T& X, const T& Y) const;
    
    /// 在数值点求值dRdX
    std::vector<std::vector<T>> evaluateDRdX(const T& X, const T& Y) const;
    
    /// 在数值点求值dRdY
    std::vector<std::vector<T>> evaluateDRdY(const T& X, const T& Y) const;
    T evaluateU(const T& X, const T& Y) const;
    T evaluateDUdX(const T& X, const T& Y) const;
    T evaluateDUdY(const T& X, const T& Y) const;
    
    /// 计算符号导数矩阵
    void computeSymbolicDerivatives();
    
    /// 类型转换：PolynomialFF → Polynomial<T>
    Polynomial<T> convertToPolynomial(const firefly::PolynomialFF& ff_poly);
    
    /// 类型转换：RationalFunctionFF → Rational<T>
    Rational<T> convertToRational(const firefly::RationalFunctionFF& ff_result);
    
    friend class firefly::DEBlackBox<T>;
};

#include "../src/de_interpolater.tpp"
