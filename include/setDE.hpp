#pragma once

#include <vector>
#include <map>
#include <memory>
#include "family.hpp"
#include "polyFamilyNumReducer.hpp"
#include "rational.hpp"

// 使用 #define 访问 Firefly 私有成员以便自定义质数
#define private public
#define protected public
#include "firefly/Reconstructor.hpp"
#undef protected
#undef private

/// DENumBuilder：在数值点(X,Y)下构建微分方程
/// 负责计算所有MFBI对X和Y的导数，并通过数值约化得到微分方程系数矩阵
template<typename T>
class DENumBuilder {
public:
    /// 构造函数
    /// @param topS 在数值点(X,Y)下的S矩阵（N×N）
    /// @param dRdX 在数值点下的dR/dX矩阵
    /// @param dRdY 在数值点下的dR/dY矩阵
    /// @param numReducer FBI数值约化器（在数值点(X,Y)下）
    /// @param family FBI族
    /// @param delta 工作维度
    DENumBuilder(const std::vector<std::vector<T>>& topS,
                 const std::vector<std::vector<T>>& dRdX,
                 const std::vector<std::vector<T>>& dRdY,
                 std::shared_ptr<PolyFamilyNumReducer<T>> numReducer,
                 const Family<T>* family,
                 T delta,
                 const T& X,
                 const T& Y);
    
    /// 计算微分方程矩阵（数值形式）
    /// @param AX 输出：X方向微分方程矩阵
    /// @param AY 输出：Y方向微分方程矩阵
    void buildDEMatrices(std::vector<std::vector<T>>& AX,
                        std::vector<std::vector<T>>& AY);
    
    // 为测试暂时设为public
    std::map<std::pair<std::vector<int>, T>, T> 
    computeMasterDerivativeX(int masterIdx);
    
    std::map<std::pair<std::vector<int>, T>, T> 
    computeMasterDerivativeY(int masterIdx);

private:
    const Family<T>* family_;
    std::shared_ptr<PolyFamilyNumReducer<T>> numReducer_;
    T delta_;
    T X_, Y_;  // 保存数值点
    int numMasterFBI_;
    int numProps_;
    int numBranch_;
    
    std::vector<std::vector<T>> topS_;
    std::vector<std::vector<T>> dRdX_;
    std::vector<std::vector<T>> dRdY_;
    
    std::map<std::pair<std::vector<int>, T>, T>
    computeFBIDerivative(const std::vector<int>& nu,
                        const std::vector<std::vector<T>>& dR);
};

/// DEBuilder：符号构建微分方程，使用Firefly插值重构
template<typename T>
class DEBuilder {
public:
    /// 构造函数
    DEBuilder(const std::vector<std::vector<Polynomial<T>>>& polyTopS,
              const Family<T>* family,
              T delta,
              uint64_t prime);
    
    /// 构建微分方程矩阵（符号形式）
    void buildDEMatrices(std::vector<std::vector<Rational<T>>>& AX,
                        std::vector<std::vector<Rational<T>>>& AY);
    
    /// 辅助：在数值点(X,Y)创建DENumBuilder
    std::unique_ptr<DENumBuilder<T>> createNumBuilder(const T& X, const T& Y);

private:
    std::vector<std::vector<Polynomial<T>>> S_;
    std::vector<std::vector<Polynomial<T>>> dRdX_;
    std::vector<std::vector<Polynomial<T>>> dRdY_;
    
    const Family<T>* family_;
    std::shared_ptr<PolyFamilyNumReducer<T>> numReducer_;
    T delta_;
    uint64_t prime_;
    int numMasterFBI_;
    int numProps_;
    int numBranch_;
    
    Rational<T> convertToRational(const firefly::RationalFunctionFF& ff_result);
};

// BlackBox单独定义，简化设计
namespace firefly {
    template<typename T>
    class DEMatrixBlackBox : public BlackBoxBase<DEMatrixBlackBox<T>> {
    public:
        DEMatrixBlackBox(DEBuilder<T>* builder) : builder_(builder) {}
        
        template<typename FFIntTemp>
        std::vector<FFIntTemp> operator()(const std::vector<FFIntTemp>& values) {
            // 提取X和Y
            T X, Y;
            if constexpr (std::is_same_v<FFIntTemp, FFInt>) {
                X = T(values[0].n);
                Y = T(values[1].n);
            } else {
                X = T(values[0][0].n);
                Y = T(values[1][0].n);
            }
            
            // 创建数值builder并计算
            auto numBuilder = builder_->createNumBuilder(X, Y);
            std::vector<std::vector<T>> AX, AY;
            numBuilder->buildDEMatrices(AX, AY);
            
            // 扁平化结果
            std::vector<FFIntTemp> result;
            int M = AX.size();
            result.reserve(2 * M * M);
            
            for (int i = 0; i < M; ++i) {
                for (int j = 0; j < M; ++j) {
                    result.emplace_back(AX[i][j].n);
                }
            }
            for (int i = 0; i < M; ++i) {
                for (int j = 0; j < M; ++j) {
                    result.emplace_back(AY[i][j].n);
                }
            }
            
            return result;
        }
        
        inline void prime_changed() {}
        
    private:
        DEBuilder<T>* builder_;
    };
}

#include "../src/setDE.tpp"
