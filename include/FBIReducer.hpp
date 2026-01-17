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
    FBIReducer(const Family<T>* family) : family_(family) {};

    const std::vector<T>& getReductionCoeff(const std::vector<int>& nu, T delta); 

    inline const Family<T>* getFamily() const { return family_; }
    void clearCache();
    size_t getCacheSize() const;
    void printCache() const;

private:
    const Family<T>* family_;   
    using CacheKey = std::tuple<std::vector<int>, T>;
    std::map<CacheKey, std::vector<T>> cache_;
    
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

#include "../src/FBIReducer.tpp"
