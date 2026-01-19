#pragma once

#include <vector>
#include <map>
#include <memory>
#include "fbi_reducer.hpp"

/// DEBuilder：在数值点(X,Y)下构建微分方程矩阵
/// 持有FBIReducer对象，负责计算MFBI对X和Y的导数并构建微分方程系数矩阵
template<typename T>
class DEBuilder {
public:
    /// 构造函数
    /// @param topS 在数值点(X,Y)下的S矩阵（N+B × N+B）
    /// @param dRdX 在数值点(X,Y)下的dR/dX矩阵（N×N）
    /// @param dRdY 在数值点(X,Y)下的dR/dY矩阵（N×N）
    /// @param numProps 传播子数量N
    /// @param numBranch 分支数量B
    /// @param delta 工作维度
    DEBuilder(const std::vector<std::vector<T>>& topS,
              const std::vector<std::vector<T>>& dRdX,
              const std::vector<std::vector<T>>& dRdY,
              int numProps,
              int numBranch,
              T delta);
    
    /// 构建微分方程矩阵（数值形式）
    /// @param AX 输出：X方向微分方程矩阵
    /// @param AY 输出：Y方向微分方程矩阵
    void buildDEMatrices(std::vector<std::vector<T>>& AX,
                        std::vector<std::vector<T>>& AY);
    
    /// 访问FBIReducer
    inline const FBIReducer<T>& getReducer() const { return reducer_; }
    
    // 为测试暂时设为public
    std::map<std::pair<std::vector<int>, T>, T> 
    computeMasterDerivativeX(int masterIdx);
    
    std::map<std::pair<std::vector<int>, T>, T> 
    computeMasterDerivativeY(int masterIdx);

private:
    FBIReducer<T> reducer_;                  // 持有约化器
    std::vector<std::vector<T>> dRdX_;      // dR/dX 矩阵（N×N）
    std::vector<std::vector<T>> dRdY_;      // dR/dY 矩阵（N×N）
    
    int numMasterFBI_;
    int numProps_;
    int numBranch_;
    
    /// 辅助方法：计算导数矩阵
    void computeDerivatives(const std::vector<std::vector<T>>& topS);
    
    /// 计算FBI对变量的导数
    std::map<std::pair<std::vector<int>, T>, T>
    computeFBIDerivative(const std::vector<int>& nu,
                        const std::vector<std::vector<T>>& dR);
};

#include "../src/de_builder.tpp"
