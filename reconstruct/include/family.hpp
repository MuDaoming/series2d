#pragma once

#include <vector>
#include <algorithm>
#include "sector.hpp"

template<typename T>
class Family {

public:
    // 构造函数
    Family(const std::vector<std::vector<T>>& topS, int numProps, int numBranch, T targetDelta = T(0))
        : topS_(topS), numProps_(numProps), numBranch_(numBranch), targetDelta_(targetDelta)
    {
        constructBranchIndices();
        findSectors();
    }

    inline int getNumProps() const { return numProps_; }
    inline int getNumBranch() const { return numBranch_; }
    inline int getCase(std::vector<int> nu) const;
    inline int getIndexOfMaster(std::vector<int> nu) const;
    const Sector<T>* getSectorByIdx(int idx) const;
    const Sector<T>* getSectorBySecvec(std::vector<int> secvec) const;
    const Sector<T>* getSector(std::vector<int> nu) const {
        std::vector<int> secvec(nu.size());
        for (size_t i = 0; i < nu.size(); i++) { secvec[i] = (nu[i] > 0) ? 1 : 0; }
        return getSectorBySecvec(secvec);
    }
    
    std::vector<int> secvecFromIdx(int n) const;          // 整数转 secvec 向量
    inline int idxFromSecvec(const std::vector<int>& secvec) const;  // secvec 向量转整数
    inline int nBranch(std::vector<int> const& nu) const;          // 计算 nu 对应的 branch 数量
    inline int nProps(std::vector<int> const& nu) const;           // 计算 nu 中非零元素个数
    inline int nuSum(std::vector<int> const& nu) const;            // 计算 nu 所有元素之和
    inline int getNumMaster() const { return masterIdxs_.size(); } // 获取 master FBI 数量
    inline bool isMaster(const std::vector<int>& nu) const;   // 判断给定 nu 是否对应 master FBI
    inline const std::vector<int>& getMasterIdxs() const { return masterIdxs_; } // 获取所有 master FBI 的索引
    inline int getCaseByIdx(int idx) const { return (idx >= 0 && idx < (int)cases_.size()) ? cases_[idx] : -1; } // 根据 idx 获取对应的 case
    inline const std::vector<std::vector<T>>& getTopS() const { return topS_; } // 获取topS矩阵
    inline const std::vector<T>& getMasterDeltas() const { return masterDeltas_; } // 获取所有主积分的delta值
    inline T getMasterDelta(int masterIdx) const { return masterDeltas_[masterIdx]; } // 获取指定主积分的delta值
    inline void setMasterDelta(T targetDelta) { masterDeltas_.assign(masterIdxs_.size(), targetDelta); } // 设置所有主积分的delta为统一值

private:

    std::vector<std::vector<T>> topS_;    // 顶层 Gram 矩阵
    int numProps_;                        // 传播子数量 N
    int numBranch_;                       // Branch 数量 B
    T targetDelta_;                       // 主积分的目标维度

    std::vector<int> branchIndices_;      // size N, 每个传播子对应的 branch 索引, 值在 [0, B-1] 之间
    std::vector<Sector<T>> sectors_;      // 所有有效的 sector
    std::vector<int> sectorIdxs_;         // 每个 sector 对应的 nu（整数表示）
    std::vector<int> cases_; // size 2^N, 存储每个 sector 的 case 类型, 0,1,2,3 分别对应四种 case, -1 表示无效 sector
    std::vector<int> masterIdxs_;    // 所有 master FBI 对应的 sectorIdx
    std::vector<T> masterDeltas_;    // 所有 master FBI 对应的 delta 值


    void constructBranchIndices();
    bool getSubS(const std::vector<int>& nu, std::vector<std::vector<T>>& subS) const;
    void findSectors();

};

#include "../src/family.tpp"
