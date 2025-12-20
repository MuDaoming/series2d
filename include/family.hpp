#pragma once

#include <vector>
#include <algorithm>
#include "sector.hpp"

template<typename T>
class Family {
private:
    // 基本参数
    std::vector<std::vector<T>> topS;  // 顶层 Gram 矩阵
    int numProps;                       // 传播子数量 N
    int numBranch;                      // Branch 数量 B
    
    // Sector 相关
    std::vector<Sector<T>> sectors;     // 所有有效的 sector
    std::vector<int> sectorNus;         // 每个 sector 对应的 nu（整数表示）
    std::vector<std::vector<int>> cases; // 长度为4，cases[i]存储属于case i的sector的nu
    
    // 辅助函数：nu 和整数的转换
    std::vector<int> intToNu(int n) const;          // 整数转 nu 向量
    int nuToInt(const std::vector<int>& nu) const;  // nu 向量转整数
    
    // 构造所有有效的 sectors
    void findSectors();

public:
    // 构造函数
    Family(const std::vector<std::vector<T>>& top_s, int n_props, int n_branch);
    
    // Getter 方法
    inline const std::vector<std::vector<T>>& getTopS() const { return topS; }
    inline int getNumProps() const { return numProps; }
    inline int getNumBranch() const { return numBranch; }
    inline const std::vector<Sector<T>>& getSectors() const { return sectors; }
    inline const std::vector<int>& getSectorNus() const { return sectorNus; }
    inline const std::vector<std::vector<int>>& getCases() const { return cases; }
    inline int getNumSectors() const { return sectors.size(); }
    
    // 根据 nu（整数）获取对应的 sector
    const Sector<T>* getSectorByNu(int nu) const;
};

#include "../src/family.tpp"
