/**
 * Family class - 三参数模板版本
 * 管理所有sector，从topS构造各个子sector
 * 
 * 模板参数：
 *   - RT: 有理函数类型（如 GiNaC::ex 或 Rational<FlintMod>）
 *   - PT: 多项式类型（如 GiNaC::ex 或 Polynomial<FlintMod>）
 *   - ST: 标量类型（如 GiNaC::ex 或 FlintMod，用于delta等数值）
 */

#pragma once

#include <vector>
#include <algorithm>
#include <sstream>
#include "sector.hpp"
#include "sector_map.hpp"

template<typename RT, typename PT, typename ST>
class Family {

    // 声明友元用于转换
    template<typename, typename, typename> friend class Family;

public:
    // 构造函数（带X和Y符号变量，自动计算dRdX和dRdY）
    template<typename Symbol>
    Family(const std::vector<std::vector<PT>>& topS, int numProps, int numBranch,
           const Symbol& X, const Symbol& Y)
        : topS_(topS), numProps_(numProps), numBranch_(numBranch),
          sectorMap_(numProps)
    {
        constructBranchIndices();
        computeDerivatives(X, Y);
        findSectors();
        setMasterDelta(ST(0));  // delta 使用 ST 类型
    }

    // 转换构造函数：从另一种类型的Family转换
    // PolyConverter: 将OtherPT转换为PT的函数对象
    // RatConverter: 将OtherRT转换为RT的函数对象
    // ScalarConverter: 将OtherST转换为ST的函数对象
    template<typename OtherRT, typename OtherPT, typename OtherST,
             typename PolyConverter, typename RatConverter, typename ScalarConverter>
    Family(const Family<OtherRT, OtherPT, OtherST>& other,
           PolyConverter polyConv, RatConverter ratConv, ScalarConverter scalarConv)
        : numProps_(other.numProps_), numBranch_(other.numBranch_),
          branchIndices_(other.branchIndices_),
          sectorIdxs_(other.sectorIdxs_),
          cases_(other.cases_),
          masterIdxs_(other.masterIdxs_),
          sectorMap_(other.sectorMap_)
    {
        // 转换topS_矩阵
        topS_.resize(other.topS_.size());
        for (size_t i = 0; i < other.topS_.size(); ++i) {
            topS_[i].resize(other.topS_[i].size());
            for (size_t j = 0; j < other.topS_[i].size(); ++j) {
                topS_[i][j] = polyConv(other.topS_[i][j]);
            }
        }
        
        // 转换dRdX_和dRdY_
        dRdX_.resize(other.dRdX_.size());
        dRdY_.resize(other.dRdY_.size());
        for (size_t i = 0; i < other.dRdX_.size(); ++i) {
            dRdX_[i].resize(other.dRdX_[i].size());
            dRdY_[i].resize(other.dRdY_[i].size());
            for (size_t j = 0; j < other.dRdX_[i].size(); ++j) {
                dRdX_[i][j] = polyConv(other.dRdX_[i][j]);
                dRdY_[i][j] = polyConv(other.dRdY_[i][j]);
            }
        }
        
        // 转换sectors_
        sectors_.reserve(other.sectors_.size());
        for (size_t i = 0; i < other.sectors_.size(); ++i) {
            try {
                sectors_.emplace_back(other.sectors_[i], polyConv, ratConv);
            } catch (const std::exception& e) {
                int secIdx = (i < other.sectorIdxs_.size()) ? other.sectorIdxs_[i] : -1;
                std::ostringstream oss;
                oss << "Family convert failed at sector[" << i << "]";
                if (secIdx >= 0) {
                    oss << ", idx=" << secIdx << ", secvec={";
                    std::vector<int> secvec = other.secvecFromIdx(secIdx);
                    for (size_t k = 0; k < secvec.size(); ++k) {
                        oss << secvec[k];
                        if (k + 1 < secvec.size()) oss << ",";
                    }
                    oss << "}";
                }
                oss << ", reason: " << e.what();
                throw std::runtime_error(oss.str());
            }
        }
        
        // 转换masterDeltas_
        masterDeltas_.resize(other.masterDeltas_.size());
        for (size_t i = 0; i < other.masterDeltas_.size(); ++i) {
            masterDeltas_[i] = scalarConv(other.masterDeltas_[i]);
        }
    }

    inline int getNumProps() const { return numProps_; }
    inline int getNumBranch() const { return numBranch_; }
    inline int getCase(std::vector<int> nu) const;
    inline int getIndexOfMaster(std::vector<int> nu) const;
    const Sector<RT, PT>* getSectorByIdx(int idx) const;
    const Sector<RT, PT>* getSectorBySecvec(std::vector<int> secvec) const;
    const Sector<RT, PT>* getSector(std::vector<int> nu) const {
        std::vector<int> secvec(nu.size());
        for (size_t i = 0; i < nu.size(); i++) { secvec[i] = (nu[i] > 0) ? 1 : 0; }
        return getSectorBySecvec(secvec);
    }
    
    std::vector<int> secvecFromIdx(int n) const;
    inline int idxFromSecvec(const std::vector<int>& secvec) const;
    inline int nBranch(std::vector<int> const& nu) const;
    inline int nProps(std::vector<int> const& nu) const;
    inline int nuSum(std::vector<int> const& nu) const;
    inline int getNumMaster() const { return masterIdxs_.size(); }
    inline bool isMaster(const std::vector<int>& nu) const;
    inline std::vector<int> canonicalizeNu(const std::vector<int>& nu) const {
        return sectorMap_.canonicalizeNu(nu);
    }
    inline std::vector<int> canonicalizeSector(const std::vector<int>& sector) const {
        return sectorMap_.canonicalizeSector(sector);
    }
    inline const SectorMap& getSectorMap() const { return sectorMap_; }
    void setSectorMap(const SectorMap& sectorMap);
    inline const std::vector<int>& getMasterIdxs() const { return masterIdxs_; }
    inline int getCaseByIdx(int idx) const { return (idx >= 0 && idx < (int)cases_.size()) ? cases_[idx] : -1; }
    inline const std::vector<std::vector<PT>>& getTopS() const { return topS_; }
    inline const std::vector<int>& getBranchIndices() const { return branchIndices_; }
    inline const std::vector<ST>& getMasterDeltas() const { return masterDeltas_; }
    inline ST getMasterDelta(int masterIdx) const { return masterDeltas_[masterIdx]; }
    inline void setMasterDelta(ST Delta) { masterDeltas_.assign(masterIdxs_.size(), Delta); }

    // 获取所有有效的sector
    const std::vector<Sector<RT, PT>>& getSectors() const { return sectors_; }
    const std::vector<int>& getSectorIdxs() const { return sectorIdxs_; }

    // 获取dR/dX和dR/dY矩阵（N×N，多项式矩阵）
    const std::vector<std::vector<PT>>& getDRdX() const { return dRdX_; }
    const std::vector<std::vector<PT>>& getDRdY() const { return dRdY_; }

private:

    std::vector<std::vector<PT>> topS_;  // 多项式矩阵
    int numProps_;
    int numBranch_;

    std::vector<int> branchIndices_;
    std::vector<Sector<RT, PT>> sectors_;
    std::vector<int> sectorIdxs_;
    std::vector<int> cases_;
    std::vector<int> masterIdxs_;
    std::vector<ST> masterDeltas_;  // delta 使用 ST 类型
    SectorMap sectorMap_;

    // R矩阵的导数（N×N，多项式矩阵）
    std::vector<std::vector<PT>> dRdX_;  // dR/dX，R是topS的右下角N×N子矩阵
    std::vector<std::vector<PT>> dRdY_;  // dR/dY

    void constructBranchIndices();
    template<typename Symbol>
    void computeDerivatives(const Symbol& X, const Symbol& Y);  // 计算dRdX和dRdY
    bool getSubS(const std::vector<int>& nu, std::vector<std::vector<PT>>& subS) const;
    void findSectors();
    void rebuildMasters();

};

#include "../src/family.tpp"
