#pragma once

#include "polyFamilyNumReducer.hpp"
#include "rational.hpp"

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
class PolyFamilyReducer;

namespace firefly {
    /// BlackBox 类：用于单个FBI的约化系数插值
    /// 在给定(X, Y)时，返回所有 numMasterFBIs 个约化系数
    template<typename T>
    class PolyFamilyBlackBox : public BlackBoxBase<PolyFamilyBlackBox<T>> {
    public:
        /// 构造函数
        PolyFamilyBlackBox(std::shared_ptr<PolyFamilyNumReducer<T>> reducer,
                           const std::vector<int>& nu,
                           const T& delta)
            : reducer_(reducer), nu_(nu), delta_(delta) {}
        
        /// BlackBox 评估函数
        /// values[0] = X, values[1] = Y
        /// 返回所有 master FBI 的约化系数向量
        template<typename FFIntTemp>
        std::vector<FFIntTemp> operator()(const std::vector<FFIntTemp>& values) {
            // 从 values 中提取 X 和 Y
            T X, Y;
            if constexpr (std::is_same_v<FFIntTemp, FFInt>) {
                X = T(values[0].n);
                Y = T(values[1].n);
            } else {
                X = T(values[0][0].n);
                Y = T(values[1][0].n);
            }
            
            // 获取约化系数向量
            std::vector<T> coeffs = reducer_->getReductionCoeff(nu_, delta_, X, Y);
            
            // 转换为 FFIntTemp 类型
            std::vector<FFIntTemp> result;
            result.reserve(coeffs.size());
            for (const auto& c : coeffs) {
                // FFInt可以直接构造
                result.emplace_back(c.n);
            }
            
            return result;
        }
        
        /// 质数变化时调用
        inline void prime_changed() {
            // 质数已经自动在 FFInt::p 中更新
        }
        
    private:
        std::shared_ptr<PolyFamilyNumReducer<T>> reducer_;  ///< 数值约化器
        std::vector<int> nu_;                               ///< FBI 指数向量
        T delta_;                                            ///< delta 参数
    };
    
    /// BatchBlackBox 类：用于批量插值多个FBI的约化系数
    /// 在给定(X, Y)时，返回所有FBI的所有约化系数（扁平化）
    template<typename T>
    class PolyFamilyBatchBlackBox : public BlackBoxBase<PolyFamilyBatchBlackBox<T>> {
    public:
        /// 构造函数
        PolyFamilyBatchBlackBox(std::shared_ptr<PolyFamilyNumReducer<T>> reducer,
                                const std::vector<std::pair<std::vector<int>, T>>& fbi_list)
            : reducer_(reducer), fbi_list_(fbi_list) {}
        
        /// BlackBox 评估函数
        /// values[0] = X, values[1] = Y
        /// 返回所有FBI的约化系数向量（扁平化）
        template<typename FFIntTemp>
        std::vector<FFIntTemp> operator()(const std::vector<FFIntTemp>& values) {
            // 从 values 中提取 X 和 Y
            T X, Y;
            if constexpr (std::is_same_v<FFIntTemp, FFInt>) {
                X = T(values[0].n);
                Y = T(values[1].n);
            } else {
                X = T(values[0][0].n);
                Y = T(values[1][0].n);
            }
            
            // 存储所有FBI的约化系数
            std::vector<FFIntTemp> result;
            
            // 遍历所有FBI
            for (const auto& [nu, delta] : fbi_list_) {
                // 获取当前FBI的约化系数
                std::vector<T> coeffs = reducer_->getReductionCoeff(nu, delta, X, Y);
                
                // 添加到结果向量
                for (const auto& c : coeffs) {
                    result.emplace_back(c.n);
                }
            }
            
            return result;
        }
        
        /// 质数变化时调用
        inline void prime_changed() {
            // 质数已经自动在 FFInt::p 中更新
        }
        
    private:
        std::shared_ptr<PolyFamilyNumReducer<T>> reducer_;           ///< 数值约化器
        std::vector<std::pair<std::vector<int>, T>> fbi_list_;       ///< 所有FBI的(nu, delta)对
    };
}

/// PolyFamilyReducer 类：符号约化，返回有理函数形式的约化系数
/// 使用批量插值重构得到每个master FBI系数关于(X,Y)的有理函数表达式
/// @tparam T 有限域类型（如 FlintMod）
template<typename T>
class PolyFamilyReducer {
public:
    /// 构造函数
    /// @param polyTopS topS 矩阵，每个元素是 Polynomial<T>
    /// @param numProps 传播子数量
    /// @param numBranch Branch 数量
    /// @param prime 有限域的素数模数
    PolyFamilyReducer(const std::vector<std::vector<Polynomial<T>>>& polyTopS,
                      int numProps,
                      int numBranch,
                      uint64_t prime);
    
    /// 主要接口：获取符号约化系数（有理函数形式）
    /// @param nu FBI 的指数向量
    /// @param delta delta 参数
    /// @return 约化系数向量，每个元素是 Rational<T>（关于X,Y的有理函数）
    std::vector<Rational<T>> getReductionCoeff(const std::vector<int>& nu, const T& delta);
    
    /// 批量接口：一次性获取多个FBI的约化系数（所有函数一起插值）
    /// @param all_nus 所有FBI的指数向量
    /// @param delta delta 参数
    /// @return 扁平化的约化系数向量：[FBI1的系数..., FBI2的系数..., ...]
    std::vector<Rational<T>> getBatchReductionCoeff(
        const std::vector<std::pair<std::vector<int>, T>>& fbi_list);
    
    /// Debug版本：获取符号约化系数，并输出每个系数使用的(X,Y)点数
    /// @param nu FBI 的指数向量
    /// @param delta delta 参数
    /// @return 约化系数向量，每个元素是 Rational<T>（关于X,Y的有理函数）
    std::vector<Rational<T>> getReductionCoeffDebug(const std::vector<int>& nu, const T& delta);
    
    /// 获取总探测点数（最近一次调用的统计）
    uint32_t getTotalProbes() const { return total_probes_; }
    
    /// 获取使用的(X,Y)点数量（numReducer的cache大小）
    size_t getCacheSize() const { return numReducer_->getCacheSize(); }
    
    /// 打印统计信息
    void printStats() const;

private:
    std::shared_ptr<PolyFamilyNumReducer<T>> numReducer_;  ///< 数值约化器（shared_ptr便于传递给BlackBox）
    uint64_t prime_;                                        ///< 有限域素数
    size_t numMasterFBIs_;                                  ///< master FBI 数量
    uint32_t total_probes_;                                 ///< 总探测点数（最近一次）
    
    /// 类型转换：PolynomialFF → Polynomial<T>
    Polynomial<T> convertToPolynomial(const firefly::PolynomialFF& ff_poly);
    
    /// 类型转换：RationalFunctionFF → Rational<T>
    Rational<T> convertToRational(const firefly::RationalFunctionFF& ff_result);
};

#include "../src/polyFamilyReducer.tpp"
