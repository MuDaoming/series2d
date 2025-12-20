#pragma once

#include <vector>
#include <map>
#include <tuple>
#include <iostream>
#include "FBI.hpp"
#include "family.hpp"

// FBI 约化器
template<typename T>
class FBIReducer {
private:
    const Family<T>* family;           // 绑定的 family
    
    // 缓存键：(sectorNu, nu向量, delta) -> 约化系数
    using CacheKey = std::tuple<int, std::vector<int>, T>;
    std::map<CacheKey, std::vector<T>> cache;
    
    // 辅助函数：创建缓存键
    CacheKey makeKey(int sectorNu, const std::vector<int>& nu, T delta) const;
    
    // 辅助函数：比较两个 nu 向量是否相等
    bool nuEqual(const std::vector<int>& nu1, const std::vector<int>& nu2) const;
    
    // 辅助函数：找到 nu 中最大元素的索引
    std::pair<int, int> findMaxIndex(const std::vector<int>& nu) const;
    
    // 辅助函数：检查 nu 是否全部等于 sectorNu
    bool isMasterFBI(const std::vector<int>& nu, int sectorNu, const Family<T>* fam) const;
    
    // 辅助函数：从 nu 向量获取 sectorNu
    int getSectorNu(const std::vector<int>& nu) const;
    
    // Case 0 约化函数
    std::vector<T> reduceCase0(int sectorNu, const std::vector<int>& nu, T delta);
    
    // Case 0：步骤1 - IBP 恒等式（从 delta 到 delta+1）
    std::vector<T> case0IBP(int sectorNu, const std::vector<int>& nu, T delta);
    
    // Case 0：步骤2 - 维度变换（从 delta+1 到 delta）
    std::vector<T> case0DimensionShift(int sectorNu, const std::vector<int>& nu, 
                                        const std::vector<T>& coeffDeltaPlus1, T delta);
    
    // 通用约化函数（根据 case 调用相应的约化方法）
    std::vector<T> reduceInternal(int sectorNu, const std::vector<int>& nu, T delta);

public:
    // 构造函数
    FBIReducer(const Family<T>* fam);
    
    // 主接口：获取约化系数（自动使用缓存）
    const std::vector<T>& getReductionCoeff(const std::vector<int>& nu, T delta);
    
    // 主接口：为 FBI 对象获取约化系数
    const std::vector<T>& getReductionCoeff(const FBI<T>& fbi);
    
    // 构造 MFBIs（特定维度的所有 master FBI）
    std::vector<FBI<T>> buildMFBIs(T delta) const;
    
    // 清空缓存
    void clearCache();
    
    // 获取缓存大小
    size_t getCacheSize() const;
    
    // Getter
    inline const Family<T>* getFamily() const { return family; }
};

#include "../src/FBIReducer.tpp"
