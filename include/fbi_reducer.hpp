#pragma once

#include <vector>
#include <map>
#include <tuple>
#include <iostream>
#include "family.hpp"

// FBI 约化器
template<typename T>
class FBIReducer {
public:
    // 构造函数：持有 Family 对象，传入 delta
    FBIReducer(const std::vector<std::vector<T>>& topS, 
               int numProps, 
               int numBranch,
               T delta)
        : family_(topS, numProps, numBranch, delta), 
          delta_(delta) 
    {
        // 初始化 Master FBI 信息
        numMaster_ = family_.getMasterIdxs().size();
        masterDeltas_.assign(numMaster_, delta);  // 所有 Master FBI 用同一个 delta
        
        // 提取每个 Master FBI 的 nu
        masterNus_.reserve(numMaster_);
        for (int idx : family_.getMasterIdxs()) {
            masterNus_.push_back(family_.secvecFromIdx(idx));
        }
    }

    const std::vector<T>& getReductionCoeff(const std::vector<int>& nu, T delta); 

    inline const Family<T>& getFamily() const { return family_; }
    inline int getNumMaster() const { return numMaster_; }
    inline const std::vector<T>& getMasterDeltas() const { return masterDeltas_; }
    inline const std::vector<std::vector<int>>& getMasterNus() const { return masterNus_; }
    
    void clearCache();
    size_t getCacheSize() const;
    void printCache() const;

private:
    Family<T> family_;                           // 持有 Family 对象
    T delta_;                                     // 工作维度
    
    // Master FBI 信息
    int numMaster_;                              // Master FBI 数量
    std::vector<T> masterDeltas_;               // 每个 Master FBI 的 delta
    std::vector<std::vector<int>> masterNus_;   // 每个 Master FBI 的 nu
    
    // 缓存
    using CacheKey = std::tuple<std::vector<int>, T>;
    mutable std::map<CacheKey, std::vector<T>> cache_;
    
    // main reduction function
    std::vector<T> reduceFBI(const std::vector<int>& nu, T delta);

    // case 0
    std::vector<T> reduceCase0(const std::vector<int>& nu, T delta);
    std::vector<T> case0IBP(const std::vector<int>& nu, T delta);
    std::vector<T> case0DimensionShift(const std::vector<int>& nu, T delta);
    std::vector<T> case0CornerDimensionShift(const std::vector<int>& nu, T delta, T targetDelta);
    std::vector<T> case0CornerDimensionShiftUp(const std::vector<int>& nu, T delta, T targetDelta);
    std::vector<T> case0CornerDimensionShiftDown(const std::vector<int>& nu, T delta, T targetDelta);

    // case 1
    std::vector<T> reduceCase1(const std::vector<int>& nu, T delta);

    // case 2
    std::vector<T> reduceCase2(const std::vector<int>& nu, T delta);

    // case 3
    std::vector<T> reduceCase3(const std::vector<int>& nu, T delta);

    // auxiliary functions
    CacheKey makeKey(const std::vector<int>& nu, T delta) const {return std::make_tuple(nu, delta); };
    
    // 判断nu是否为corner（所有元素都是0或1）
    bool isCorner(const std::vector<int>& nu) const {
        for (int val : nu) {
            if (val != 0 && val != 1) {
                return false;
            }
        }
        return true;
    }
    
    std::pair<int, int> findMaxIndex(const std::vector<int>& nu) const {
        int maxIdxInTopSector = 0;
        int maxVal = nu[0];
        for (size_t i = 1; i < nu.size(); i++) {
            if (nu[i] > maxVal) {
                maxVal = nu[i];
                maxIdxInTopSector = i;
            }
        }

        int maxIdxInCurrentSector = 0;
        for (int i = 0; i < maxIdxInTopSector; i++) {
            if (nu[i] > 0) {
                maxIdxInCurrentSector++;
            }
        }
        return std::pair(maxIdxInTopSector, maxIdxInCurrentSector);
    };
};

#include "../src/fbi_reducer.tpp"
