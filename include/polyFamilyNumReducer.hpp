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
    
    /// 主要接口：在给定 (X,Y) 值处进行数值约化
    /// @param nu FBI 的指数向量
    /// @param delta delta 参数
    /// @param X X 的值
    /// @param Y Y 的值
    /// @return 约化系数向量（对应所有 master FBI）
    std::vector<T> getReductionCoeff(const std::vector<int>& nu, T delta, const T& X, const T& Y);
    
    /// 获取 numProps
    inline int getNumProps() const { return numProps_; }
    
    /// 获取 numBranch
    inline int getNumBranch() const { return numBranch_; }
    
    /// 获取 master FBI 的数量
    size_t getNumMasterFBIs() const;
    
    /// 清除缓存
    void clearCache();
    
    /// 获取缓存大小
    size_t getCacheSize() const { return familyCache_.size(); }
    
    /// 打印缓存信息
    void printCacheInfo() const;

private:
    std::vector<std::vector<Polynomial<T>>> polyTopS_;  ///< 多项式 topS 矩阵
    int numProps_;                                      ///< 传播子数量
    int numBranch_;                                     ///< Branch 数量
    
    /// 缓存：(X, Y) -> (Family, FBIReducer)
    /// 使用 unique_ptr 自动管理内存
    std::map<std::pair<T, T>, std::pair<std::unique_ptr<Family<T>>, std::unique_ptr<FBIReducer<T>>>> familyCache_;
    
    /// 在 (X, Y) 处对 polyTopS 求值，得到数值 topS
    /// @param X X 的值
    /// @param Y Y 的值
    /// @return 数值 topS 矩阵
    std::vector<std::vector<T>> evaluateTopS(const T& X, const T& Y) const;
    
    /// 获取或创建 Family 和 Reducer
    /// @param X X 的值
    /// @param Y Y 的值
    /// @return Family 和 FBIReducer 的指针对
    std::pair<Family<T>*, FBIReducer<T>*> getOrCreateReducer(const T& X, const T& Y);
};

#include "../src/polyFamilyNumReducer.tpp"
