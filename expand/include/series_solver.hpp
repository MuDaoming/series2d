/**
 * SeriesSolver - 级数求解器
 * 统一处理微分方程求解和IBP约化，直接计算重定义后FBI (Ĩ = U^{pow_U} · I) 的二维幂级数展开
 * 
 * 核心思想：将所有问题统一为 LRR (Linear Recurrence Relation) 形式:
 *   D(X,Y) · g(X,Y) = sum_i N_i(X,Y) · f_i(X,Y)
 * 
 * 模板参数：
 *   - RT: 有理函数类型（如 GiNaC::ex 或 Rational<FlintMod>）
 *   - PT: 多项式类型（如 GiNaC::ex 或 Polynomial<FlintMod>）
 *   - ST: 标量类型（如 GiNaC::ex 或 FlintMod，用于级数系数和delta）
 */

#pragma once

#include <vector>
#include <map>
#include <tuple>
#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <iostream>
#include <string>
#include "family.hpp"
#include "series.hpp"
#include "rational.hpp"

// ============================================================================
// Redefinition 辅助结构体
// ============================================================================

template<typename PT, typename ST>
struct Redefinition {
    ST D_in;           // 主积分 delta = L * D_Feynman / 2
    int L;             // 圈数
    PT shiftedU;       // 平移后的 U(X,Y) 多项式
    PT dUdX;           // ∂U/∂X
    PT dUdY;           // ∂U/∂Y

    Redefinition(int L_, const ST& D_in_, const PT& U_)
        : D_in(D_in_), L(L_), shiftedU(U_),
          dUdX(U_.derivativeX()), dUdY(U_.derivativeY()) {}

    /// 将 Z_p 中的值转为小有符号整数（|val| << p）
    static int toSignedInt(const ST& val) {
        mp_limb_t v = val.get_value();
        mp_limb_t p = ST::get_modulus();
        if (v > p / 2) return static_cast<int64_t>(v) - static_cast<int64_t>(p);
        return static_cast<int64_t>(v);
    }

    /// D 与 D_in 是否同奇偶
    bool sameParity(const ST& D) const {
        int offset = toSignedInt(D - D_in);
        return (offset % 2 == 0);
    }

    int etaMinus(const ST& D) const { return sameParity(D) ? 0 : 2; }
    int etaPlus(const ST& D) const { return sameParity(D) ? -2 : 0; }

    /// 计算幂差 Δp_U（整数）
    int deltaPowU(int nuTotT, const ST& DT, int nuTotS, const ST& DS) const {
        int offsetT = toSignedInt(DT - D_in);
        int offsetS = toSignedInt(DS - D_in);
        bool parityT = (offsetT % 2 == 0);
        bool parityS = (offsetS % 2 == 0);
        int DbarT = offsetT + (parityT ? 0 : 1);
        int DbarS = offsetS + (parityS ? 0 : 1);
        int DbarDiff = DbarT - DbarS;
        return (nuTotT - nuTotS) - (L + 1) * DbarDiff / 2;
    }

    /// pow_U 在 Z_p 中的标量值（微分方程 dlog 系数）
    ST powUScalar(int nuTot, const ST& D) const {
        int offset = toSignedInt(D - D_in);
        bool parity = (offset % 2 == 0);
        int Dbar_offset = offset + (parity ? 0 : 1);
        return ST(nuTot) - ST(L + 1) * (D_in + ST(Dbar_offset)) / ST(2);
    }
};

// ============================================================================
// SeriesSolver 类
// ============================================================================
template<typename RT, typename PT, typename ST>
class SeriesSolver {
public:
    enum class ReduceMode {
        Normal,
        Cut
    };
    // ========================================================================
    // 构造与基本接口
    // ========================================================================
    
    SeriesSolver(Family<RT, PT, ST>& family, int targetDeg);
    
    void setMasterBoundary(int masterIdx, const ST& value);
    void setAllMasterBoundary(const ST& value);
    void setReduceMode(ReduceMode mode) { reduceMode_ = mode; }
    void setReduceMode(const std::string& modeName);
    void setCutSector(const std::vector<int>& sector);
    ReduceMode getReduceMode() const { return reduceMode_; }
    
    /// 设置重定义（必须在调用 getFBISeries() 之前调用）
    void setRedefinition(const Redefinition<PT, ST>* redef);

    /// 获取FBI的级数展开（递归调用,会触发约化），并确保至少到 needDeg
    const Series<ST>& getFBISeries(const std::vector<int>& nu, const ST& delta, int needDeg);
    
    // ========================================================================
    // 访问器
    // ========================================================================
    
    int getTargetDeg() const { return targetDeg_; }
    int getNumMaster() const { return numMaster_; }
    const Family<RT, PT, ST>& getFamily() const { return family_; }
    const Series<ST>& getMasterSeries(int masterIdx) const;
    size_t getCacheSize() const { return cache_.size(); }
    void clearCache();
    void printCacheInfo() const;
    void printAllCache() const;

    // ========================================================================
    // 辅助函数：多项式与级数的卷积
    // ========================================================================
    
    static ST polySeriesCoeff(const PT& poly, const Series<ST>& series, int p, int q);
    
    static ST sumPolySeriesCoeff(const std::vector<PT>& polys, 
                                  const std::vector<const Series<ST>*>& series,
                                  int p, int q);
    
    // ========================================================================
    // LRR求解器
    // ========================================================================
    
    static void solveLRRAtDeg(Series<ST>& g, const PT& D,
                               const std::vector<PT>& polys,
                               const std::vector<const Series<ST>*>& series,
                               int deg);

    // ========================================================================
    // 多项式辅助函数（静态）
    // ========================================================================
    static PT multiplyPolys(const PT& a, const PT& b);
    static PT powPolyExpand(const PT& base, int exp);

private:
    // ========================================================================
    // 成员变量
    // ========================================================================
    
    Family<RT, PT, ST>& family_;
    int targetDeg_;
    int numMaster_;
    int numProps_;
    int numBranch_;
    
    std::vector<std::vector<int>> masterNus_;
    std::vector<ST> masterDeltas_;
    std::vector<ST> masterBoundary_;
    
    // FBI缓存键
    struct CacheKey {
        std::vector<int> nu;
        ST delta;
        bool operator<(const CacheKey& other) const {
            if (nu != other.nu) return nu < other.nu;
            return delta < other.delta;
        }
    };
    
    mutable std::map<CacheKey, Series<ST>> cache_;
    mutable std::map<CacheKey, int> cacheCurrentDeg_;
    ReduceMode reduceMode_ = ReduceMode::Normal;
    std::vector<int> cutSector_;
    
    // ========================================================================
    // Redefinition 相关
    // ========================================================================
    
    const Redefinition<PT, ST>* redef_ = nullptr;
    std::vector<std::vector<PT>> dRdXModified_;   // dRdX * U^L
    std::vector<std::vector<PT>> dRdYModified_;   // dRdY * U^L
    std::vector<ST> masterPowU_;                   // 每个主积分的 pow_U
    
    // ========================================================================
    // 内部辅助方法
    // ========================================================================
    
    CacheKey makeKey(const std::vector<int>& nu, const ST& delta) const {
        return CacheKey{nu, delta};
    }
    
    bool isCorner(const std::vector<int>& nu) const;
    bool sourceContainsCutSector(const std::vector<int>& sourceNu) const;
    bool shouldDropInCut(const std::vector<int>& sourceNu) const;
    std::pair<int, int> findMaxIndex(const std::vector<int>& nu) const;
    
    static int nuTotSum(const std::vector<int>& nu) {
        int s = 0; for (int v : nu) s += v; return s;
    }

    /// 对约化函数应用 ratio 因子（通分 D、预乘分子 polys）
    void applyRatioFactors(PT& D,
                           std::vector<PT>& polys,
                           const std::vector<int>& deltaPs) const;

    // ========================================================================
    // 主求解逻辑
    // ========================================================================
    
    void solveMasterAtDeg(int masterIdx, int deg);
    void solveMasterCoeffX(int masterIdx, int p, int q);
    void solveMasterCoeffY(int masterIdx, int q);
    
    // ========================================================================
    // 约化逻辑
    // ========================================================================
    
    void reduceFBIAtDeg(Series<ST>& result, const std::vector<int>& nu, 
                        const ST& delta, int deg);
    
    void reduceCase0AtDeg(Series<ST>& result, const std::vector<int>& nu, 
                          const ST& delta, int deg);
    void case0IBPAtDeg(Series<ST>& result, const std::vector<int>& nu, 
                       const ST& delta, int deg);
    void case0DimShiftDownAtDeg(Series<ST>& result, const std::vector<int>& nu, 
                                 const ST& delta, int deg);
    void case0DimShiftUpAtDeg(Series<ST>& result, const std::vector<int>& nu, 
                               const ST& delta, int deg);
    void reduceCase1AtDeg(Series<ST>& result, const std::vector<int>& nu, 
                          const ST& delta, int deg);
    void reduceCase2AtDeg(Series<ST>& result, const std::vector<int>& nu, 
                          const ST& delta, int deg);
    void reduceCase3AtDeg(Series<ST>& result, const std::vector<int>& nu, 
                          const ST& delta, int deg);
};

#include "../src/series_solver.tpp"
