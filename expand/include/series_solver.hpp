/**
 * SeriesSolver - 级数求解器
 * 统一处理微分方程求解和IBP约化，直接计算FBI的二维幂级数展开
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
#include <stdexcept>
#include <iostream>
#include "family.hpp"
#include "series.hpp"
#include "rational.hpp"

// ============================================================================
// SeriesSolver 类
// ============================================================================
template<typename RT, typename PT, typename ST>
class SeriesSolver {
public:
    // ========================================================================
    // 构造与基本接口
    // ========================================================================
    
    /**
     * 构造函数
     * @param family Family对象的引用
     * @param targetDeg 目标级数度数
     */
    SeriesSolver(Family<RT, PT, ST>& family, int targetDeg);
    
    /**
     * 设置主积分边界条件（零阶系数）
     * @param masterIdx 主积分索引
     * @param value 零阶系数值
     */
    void setMasterBoundary(int masterIdx, const ST& value);
    
    /**
     * 设置所有主积分边界条件为相同值
     * @param value 零阶系数值
     */
    void setAllMasterBoundary(const ST& value);
    
    /**
     * 主求解函数：按度数递推求解所有主积分的级数
     */
    void solve();
    
    /**
     * 获取FBI的级数展开（递归调用,会触发约化）
     * @param nu 传播子指数向量
     * @param delta 维度参数
     * @return 级数的常量引用
     */
    const Series<ST>& getFBISeries(const std::vector<int>& nu, const ST& delta);
    
    // ========================================================================
    // 访问器
    // ========================================================================
    
    int getCurrentDeg() const { return currentDeg_; }
    int getTargetDeg() const { return targetDeg_; }
    int getNumMaster() const { return numMaster_; }
    const Series<ST>& getMasterSeries(int masterIdx) const;
    size_t getCacheSize() const { return cache_.size(); }
    void clearCache();
    void printCacheInfo() const;
    void printAllCache() const;  // 打印所有缓存的(nu,delta,series)

    // ========================================================================
    // 辅助函数：多项式与级数的卷积
    // ========================================================================
    
    /**
     * 计算 [P · f]_{p,q}：多项式P与级数f卷积后的(p,q)项系数
     * 
     * [P·f]_{p,q} = Σ_{a,b} P_{ab} · f_{p-a, q-b}
     * 
     * @param poly 多项式 P
     * @param series 级数 f
     * @param p X的指数
     * @param q Y的指数
     * @return 卷积的(p,q)项系数
     */
    static ST polySeriesCoeff(const PT& poly, const Series<ST>& series, int p, int q);
    
    /**
     * 计算 Σ_i [N_i · f_i]_{p,q}：多个多项式与级数卷积的求和
     * 
     * @param polys 多项式向量 [N_0, N_1, ...]
     * @param series 级数指针向量 [f_0, f_1, ...]
     * @param p X的指数
     * @param q Y的指数
     * @return Σ_i [N_i · f_i]_{p,q}
     */
    static ST sumPolySeriesCoeff(const std::vector<PT>& polys, 
                                  const std::vector<const Series<ST>*>& series,
                                  int p, int q);
    
    // ========================================================================
    // LRR求解器
    // ========================================================================
    
    /**
     * 使用LRR递推求解级数g在度数deg的所有系数
     * 
     * LRR形式: D · g = Σ_i N_i · f_i
     * 
     * 递推公式: g_{p,q} = (1/D_00) · ([Σ N_i·f_i]_{pq} - Σ_{(a,b)≠(0,0)} D_{ab}·g_{p-a,q-b})
     * 
     * 前提条件:
     *   - g 的 deg < deg 的系数已知
     *   - f_i 的 deg <= deg 的系数已知
     *   - D_00 ≠ 0
     * 
     * @param g 输出级数（会被修改）
     * @param D 分母多项式
     * @param polys 分子多项式向量 [N_0, N_1, ...]
     * @param series 级数指针向量 [f_0, f_1, ...]
     * @param deg 要求解的度数
     */
    static void solveLRRAtDeg(Series<ST>& g, const PT& D,
                               const std::vector<PT>& polys,
                               const std::vector<const Series<ST>*>& series,
                               int deg);

private:
    // ========================================================================
    // 成员变量
    // ========================================================================
    
    Family<RT, PT, ST>& family_;    // Family对象引用
    int targetDeg_;                  // 目标级数度数
    int currentDeg_;                 // 当前已求解的度数
    int numMaster_;                  // 主积分数量
    int numProps_;                   // 传播子数量
    int numBranch_;                  // 分支数量
    
    // 主积分信息
    std::vector<std::vector<int>> masterNus_;     // 每个主积分的nu
    std::vector<ST> masterDeltas_;                // 每个主积分的delta
    
    // 主积分边界条件
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
    
    // 统一的FBI级数缓存（包括主积分）
    mutable std::map<CacheKey, Series<ST>> cache_;
    
    // 每个缓存级数已经约化到的度数
    mutable std::map<CacheKey, int> cacheCurrentDeg_;
    
    // ========================================================================
    // 内部辅助方法
    // ========================================================================
    
    CacheKey makeKey(const std::vector<int>& nu, const ST& delta) const {
        return CacheKey{nu, delta};
    }
    
    bool isCorner(const std::vector<int>& nu) const;
    std::pair<int, int> findMaxIndex(const std::vector<int>& nu) const;
    
    // ========================================================================
    // 主求解逻辑
    // ========================================================================
    
    /**
     * 求解度数deg的所有主积分系数
     */
    void solveAtDeg(int deg);
    
    /**
     * 使用微分方程计算主积分在(p,q)的系数
     * 
     * 微分方程: ∂f/∂X = Σ_{i,j} (-1/2)·dR_{ij}/dX·factor·FBI(ν+e_i+e_j, δ+1)
     * 递推: f_{p,q} = [rhs]_{p-1,q} / p  (p>0)
     *       f_{0,q} = [rhs_Y]_{0,q-1} / q  (p=0, q>0)
     */
    void solveMasterCoeffX(int masterIdx, int p, int q);
    void solveMasterCoeffY(int masterIdx, int q);
    
    // ========================================================================
    // 约化逻辑
    // ========================================================================
    
    /**
     * 约化并求解FBI级数（整体求解,按度数递推）
     */
    void reduceFBI(Series<ST>& result, const std::vector<int>& nu, const ST& delta);
    
    /**
     * 使用LRR约化FBI在度数deg的系数
     * 
     * 根据Case类型选择不同的约化公式:
     * - Case 0: IBP约化 + 维度迁移
     * - Case 1: (2Δ-ν-B)·I = -Σz_α·I_{ν-e_α}^{Δ-1}
     * - Case 2: C·I = Σz_α·I_{ν-e_α}^Δ
     * - Case 3: z_β·I = -Σz_α·I_{ν+e_β-e_α}^Δ
     */
    void reduceFBIAtDeg(Series<ST>& result, const std::vector<int>& nu, 
                        const ST& delta, int deg);
    
    // Case 0 细分
    void reduceCase0AtDeg(Series<ST>& result, const std::vector<int>& nu, 
                          const ST& delta, int deg);
    void case0IBPAtDeg(Series<ST>& result, const std::vector<int>& nu, 
                       const ST& delta, int deg);
    void case0DimShiftDownAtDeg(Series<ST>& result, const std::vector<int>& nu, 
                                 const ST& delta, int deg);
    void case0DimShiftUpAtDeg(Series<ST>& result, const std::vector<int>& nu, 
                               const ST& delta, int deg);
    
    // Case 1-3
    void reduceCase1AtDeg(Series<ST>& result, const std::vector<int>& nu, 
                          const ST& delta, int deg);
    void reduceCase2AtDeg(Series<ST>& result, const std::vector<int>& nu, 
                          const ST& delta, int deg);
    void reduceCase3AtDeg(Series<ST>& result, const std::vector<int>& nu, 
                          const ST& delta, int deg);
};

#include "../src/series_solver.tpp"
