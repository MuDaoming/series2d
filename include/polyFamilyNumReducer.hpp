#pragma once

#include <vector>
#include <map>
#include <utility>
#include <memory>
#include "rational.hpp"
#include "family.hpp"
#include "FBIReducer.hpp"

/// PolyFamilyNumReducer 类：topS 的元素是关于 X 和 Y 的多项式
/// 提供数值约化功能：在给定的 (X,Y) 值处对 topS 求值，创建对应的 Family 和 Reducer
/// @tparam T 有限域类型（如 FlintMod）
template<typename T>
class PolyFamilyNumReducer {
public:
    /// 构造函数
    /// @param polyTopS topS 矩阵，每个元素是 Polynomial<T>
    /// @param numProps 传播子数量
    /// @param numBranch Branch 数量
    PolyFamilyNumReducer(const std::vector<std::vector<Polynomial<T>>>& polyTopS, 
                         int numProps, 
                         int numBranch);
    
    /// 析构函数：清理缓存
    ~PolyFamilyNumReducer();
    
    /// 设置当前工作点 (X, Y)
    /// 注意：会重新创建 Family 和 Reducer，清除 FBIReducer 的内部缓存
    /// @param X X 的值
    /// @param Y Y 的值
    void setCurrentPoint(const T& X, const T& Y);
    
    /// 主要接口：在当前点进行数值约化
    /// 注意：必须先调用 setCurrentPoint 设置工作点
    /// @param nu FBI 的指数向量
    /// @param delta delta 参数
    /// @return 约化系数向量（对应所有 master FBI）
    std::vector<T> getReductionCoeff(const std::vector<int>& nu, T delta);
    
    /// 获取 numProps
    inline int getNumProps() const { return numProps_; }
    
    /// 获取 numBranch
    inline int getNumBranch() const { return numBranch_; }
    
    /// 获取 master FBI 的数量
    size_t getNumMasterFBIs() const;
    
    /// 检查是否已设置当前点
    inline bool hasCurrentPoint() const { return hasCurrentPoint_; }
    
    /// 获取当前的 FBIReducer（如果存在）
    FBIReducer<T>* getCurrentReducer() { return curReducer_.get(); }

private:
    std::vector<std::vector<Polynomial<T>>> polyTopS_;  ///< 多项式 topS 矩阵
    int numProps_;                                      ///< 传播子数量
    int numBranch_;                                     ///< Branch 数量
    
    /// 当前工作点状态
    T curX_;                                            ///< 当前 X 值
    T curY_;                                            ///< 当前 Y 值
    bool hasCurrentPoint_;                              ///< 是否已设置当前点
    std::unique_ptr<Family<T>> curFamily_;              ///< 当前点的 Family
    std::unique_ptr<FBIReducer<T>> curReducer_;         ///< 当前点的 Reducer
    
    /// 在 (X, Y) 处对 polyTopS 求值，得到数值 topS
    /// @param X X 的值
    /// @param Y Y 的值
    /// @return 数值 topS 矩阵
    std::vector<std::vector<T>> evaluateTopS(const T& X, const T& Y) const;
};

#include "../src/polyFamilyNumReducer.tpp"
