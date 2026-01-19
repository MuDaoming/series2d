#pragma once

#include "rational.hpp"
#include "fbi_reducer.hpp"

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
class FBIInterpolater;

namespace firefly {
    /// BlackBox 类：用于批量插值多个FBI的约化系数
    /// 在给定(X, Y)时，返回所有FBI的所有约化系数（扁平化）
    template<typename T>
    class FBIBlackBox : public BlackBoxBase<FBIBlackBox<T>> {
    public:
        /// 构造函数
        FBIBlackBox(const FBIInterpolater<T>* interpolater,
                    const std::vector<std::pair<std::vector<int>, T>>& fbi_list)
            : interpolater_(interpolater), fbi_list_(fbi_list) {}
        
        /// BlackBox 评估函数
        /// values[0] = X, values[1] = Y
        /// 返回所有FBI的约化系数向量（扁平化）
        template<typename FFIntTemp>
        std::vector<FFIntTemp> operator()(const std::vector<FFIntTemp>& values);
        
        /// 质数变化时调用
        inline void prime_changed() {
            // 质数已经自动在 FFInt::p 中更新
        }
        
    private:
        const FBIInterpolater<T>* interpolater_;                     ///< 指向插值器
        std::vector<std::pair<std::vector<int>, T>> fbi_list_;       ///< 所有FBI的(nu, delta)对
    };
}

/// FBIInterpolater 类：符号约化，返回有理函数形式的约化系数
/// 使用批量插值重构得到每个master FBI系数关于(X,Y)的有理函数表达式
/// @tparam T 有限域类型（如 FlintMod）
template<typename T>
class FBIInterpolater {
public:
    /// 构造函数
    /// @param topS topS 矩阵，每个元素是 Polynomial<T>
    /// @param numProps 传播子数量
    /// @param numBranch Branch 数量
    /// @param delta delta 参数
    /// @param prime 有限域的素数模数
    FBIInterpolater(const std::vector<std::vector<Polynomial<T>>>& topS,
                    int numProps,
                    int numBranch,
                    T delta,
                    uint64_t prime);
    
    /// 批量接口：一次性获取多个FBI的约化系数（所有函数一起插值）
    /// @param fbi_list 所有FBI的(nu, delta)对列表
    /// @return 约化系数向量的向量，顺序与fbi_list一致
    std::vector<std::vector<Rational<T>>>
    getReductionCoeff(const std::vector<std::pair<std::vector<int>, T>>& fbi_list);
    
    /// 获取总探测点数（最近一次调用的统计）
    uint32_t getTotalProbes() const { return total_probes_; }
    
    /// 打印统计信息
    void printStats() const;
    
    /// 在数值点创建FBIReducer（供BlackBox使用）
    FBIReducer<T> createReducer(const T& X, const T& Y) const;

private:
    std::vector<std::vector<Polynomial<T>>> topS_;  ///< 多项式形式的topS矩阵
    int numProps_;                                   ///< 传播子数量
    int numBranch_;                                  ///< Branch数量
    T delta_;                                        ///< delta参数
    uint64_t prime_;                                 ///< 有限域素数
    size_t numMasterFBIs_;                          ///< master FBI数量
    uint32_t total_probes_;                         ///< 总探测点数（最近一次）
    
    /// 在数值点求值topS
    std::vector<std::vector<T>> evaluateTopS(const T& X, const T& Y) const;
    
    /// 类型转换：PolynomialFF → Polynomial<T>
    Polynomial<T> convertToPolynomial(const firefly::PolynomialFF& ff_poly);
    
    /// 类型转换：RationalFunctionFF → Rational<T>
    Rational<T> convertToRational(const firefly::RationalFunctionFF& ff_result);
    
    friend class firefly::FBIBlackBox<T>;
};

#include "../src/fbi_interpolater.tpp"
