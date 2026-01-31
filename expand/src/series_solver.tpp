/**
 * SeriesSolver 实现文件
 * 统一处理微分方程求解和IBP约化，直接计算FBI的二维幂级数展开
 */

// ============================================================================
// 构造函数与基本接口
// ============================================================================

template<typename RT, typename PT, typename ST>
SeriesSolver<RT, PT, ST>::SeriesSolver(Family<RT, PT, ST>& family, int targetDeg)
    : family_(family), targetDeg_(targetDeg), currentDeg_(-1)
{
    numProps_ = family_.getNumProps();
    numBranch_ = family_.getNumBranch();
    numMaster_ = family_.getNumMaster();
    
    // 提取每个主积分的 nu 和 delta
    masterNus_.reserve(numMaster_);
    masterDeltas_ = family_.getMasterDeltas();
    
    for (int idx : family_.getMasterIdxs()) {
        masterNus_.push_back(family_.secvecFromIdx(idx));
    }
    
    // 初始化主积分边界条件
    masterBoundary_.resize(numMaster_, ST(0));
}

// ============================================================================
// 主求解函数
// ============================================================================

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::solve() {
    // 初始化主积分零阶系数
    for (int k = 0; k < numMaster_; ++k) {
        auto key = makeKey(masterNus_[k], masterDeltas_[k]);
        cache_[key].setCoeff(0, 0, masterBoundary_[k]);
    }
    currentDeg_ = 0;
    
    // 逐度数递推
    for (int deg = 1; deg <= targetDeg_; ++deg) {
        solveAtDeg(deg);
        currentDeg_ = deg;
    }
}

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::solveAtDeg(int deg) {
    // 对于度数deg，递推所有主积分的系数
    // 使用微分方程: ∂f/∂X 或 ∂f/∂Y
    
    for (int p = 0; p <= deg; ++p) {
        int q = deg - p;
        
        // 对每个主积分
        for (int k = 0; k < numMaster_; ++k) {
            if (p > 0) {
                // 使用X方向微分方程
                solveMasterCoeffX(k, p, q);
            } else if (q > 0) {
                // p = 0, q > 0：使用Y方向微分方程
                solveMasterCoeffY(k, q);
            }
            // p = 0, q = 0 的情况已经在初始化时处理
        }
    }
}

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::solveMasterCoeffX(int masterIdx, int p, int q) {
    // 微分方程: ∂f/∂X = rhs
    // 其中 rhs = Σ_{i,j} (-1/2) * dR_{ij}/dX * factor_{ij} * I_{ν+e_i+e_j}^{Δ+1}
    // 递推: p * f_{p,q} = [rhs]_{p-1,q}
    // 即: f_{p,q} = [rhs]_{p-1,q} / p
    
    const std::vector<int>& nu = masterNus_[masterIdx];
    const ST& delta = masterDeltas_[masterIdx];
    const auto& dRdX = family_.getDRdX();
    
    ST rhsCoeff = ST(0);
    
    // 计算 [rhs]_{p-1,q} = Σ_{i,j} (-1/2) * [dR_{ij}/dX * I_{ν+e_i+e_j}^{Δ+1}]_{p-1,q} * factor_{ij}
    for (int i = 0; i < numProps_; ++i) {
        for (int j = 0; j < numProps_; ++j) {
            // 计算 factor_{ij}
            ST factor_ij;
            if (i != j) {
                factor_ij = ST(nu[i]) * ST(nu[j]);
            } else {
                factor_ij = ST(nu[i]) * ST(nu[i] + 1);
            }
            
            // 跳过 factor = 0 or dRdX = 0 的情况
            if (isZero(factor_ij) || isZero(dRdX[i][j])) continue;
            
            // 构造 ν + e_i + e_j
            std::vector<int> nuShifted = nu;
            nuShifted[i]++;
            nuShifted[j]++;
            
            // 获取 I_{ν+e_i+e_j}^{Δ+1}
            const Series<ST>& seriesShifted = getFBISeries(nuShifted, delta + ST(1));
            
            // 计算 [dR_{ij}/dX * I_{ν+e_i+e_j}^{Δ+1}]_{p-1,q}
            ST convCoeff = polySeriesCoeff(dRdX[i][j], seriesShifted, p - 1, q);
            
            // 累加 (-1/2) * factor_ij * convCoeff
            // 注意: 在有限域中 -1/2 需要特殊处理
            // -1/2 = -(p+1)/2 mod p，但这里简化用 ST(-1)/ST(2)
            rhsCoeff -= factor_ij * convCoeff / ST(2);
        }
    }
    
    // f_{p,q} = rhsCoeff / p
    auto key = makeKey(nu, delta);
    cache_[key].setCoeff(p, q, rhsCoeff / ST(p));
}

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::solveMasterCoeffY(int masterIdx, int q) {
    // 微分方程: ∂f/∂Y = rhs
    // 其中 rhs = Σ_{i,j} (-1/2) * dR_{ij}/dY * factor_{ij} * I_{ν+e_i+e_j}^{Δ+1}
    // 递推: q * f_{0,q} = [rhs]_{0,q-1}
    // 即: f_{0,q} = [rhs]_{0,q-1} / q
    
    const std::vector<int>& nu = masterNus_[masterIdx];
    const ST& delta = masterDeltas_[masterIdx];
    const auto& dRdY = family_.getDRdY();
    
    ST rhsCoeff = ST(0);
    
    // 计算 [rhs]_{0,q-1} = Σ_{i,j} (-1/2) * [dR_{ij}/dY * I_{ν+e_i+e_j}^{Δ+1}]_{0,q-1} * factor_{ij}
    for (int i = 0; i < numProps_; ++i) {
        for (int j = 0; j < numProps_; ++j) {
            // 计算 factor_{ij}
            ST factor_ij;
            if (i != j) {
                factor_ij = ST(nu[i]) * ST(nu[j]);
            } else {
                factor_ij = ST(nu[i]) * ST(nu[i] + 1);
            }
            
            // 跳过 factor = 0 or dRdY = 0 的情况
            if (isZero(factor_ij) || isZero(dRdY[i][j])) continue;

            // 构造 ν + e_i + e_j
            std::vector<int> nuShifted = nu;
            nuShifted[i]++;
            nuShifted[j]++;
            
            // 获取 I_{ν+e_i+e_j}^{Δ+1}
            const Series<ST>& seriesShifted = getFBISeries(nuShifted, delta + ST(1));
            
            // 计算 [dR_{ij}/dY * I_{ν+e_i+e_j}^{Δ+1}]_{0,q-1}
            ST convCoeff = polySeriesCoeff(dRdY[i][j], seriesShifted, 0, q - 1);
            
            // 累加 (-1/2) * factor_ij * convCoeff
            rhsCoeff -= factor_ij * convCoeff / ST(2);
        }
    }
    
    // f_{0,q} = rhsCoeff / q
    auto key = makeKey(nu, delta);
    cache_[key].setCoeff(0, q, rhsCoeff / ST(q));
}

// ============================================================================
// FBI级数获取（递归调用，触发约化）
// ============================================================================

template<typename RT, typename PT, typename ST>
const Series<ST>& SeriesSolver<RT, PT, ST>::getFBISeries(const std::vector<int>& nu, const ST& delta) {
    auto key = makeKey(nu, delta);
    auto it = cache_.find(key);
    
    if (it != cache_.end()) {
        // 缓存存在，检查是否需要继续约化到更高度数
        int cachedDeg = cacheCurrentDeg_[key];
        if (cachedDeg < currentDeg_) {
            // 需要继续约化
            for (int deg = cachedDeg + 1; deg <= currentDeg_ && deg <= targetDeg_; ++deg) {
                reduceFBIAtDeg(it->second, nu, delta, deg);
            }
            cacheCurrentDeg_[key] = std::min(currentDeg_, targetDeg_);
        }
        return it->second;
    }
    
    // 不在缓存中，需要创建并约化
    cache_[key] = Series<ST>(targetDeg_);
    for (int deg = 0; deg <= currentDeg_ && deg <= targetDeg_; ++deg) {
        reduceFBIAtDeg(cache_[key], nu, delta, deg);
    }
    cacheCurrentDeg_[key] = std::min(currentDeg_, targetDeg_);
    
    return cache_[key];
}

// ============================================================================
// 约化主逻辑
// ============================================================================

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::reduceFBI(Series<ST>& result, const std::vector<int>& nu, const ST& delta) {
    // 逐度数约化
    result.setCoeff(0, 0, ST(0));  // 初始化零阶
    
    for (int deg = 0; deg <= currentDeg_ && deg <= targetDeg_; ++deg) {
        reduceFBIAtDeg(result, nu, delta, deg);
    }
}

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::reduceFBIAtDeg(Series<ST>& result, const std::vector<int>& nu, 
                                                const ST& delta, int deg) {
    // 根据Case类型选择约化方法
    int caseType = family_.getCase(nu);
    
    switch (caseType) {
        case 0:
            reduceCase0AtDeg(result, nu, delta, deg);
            break;
        case 1:
            reduceCase1AtDeg(result, nu, delta, deg);
            break;
        case 2:
            reduceCase2AtDeg(result, nu, delta, deg);
            break;
        case 3:
            reduceCase3AtDeg(result, nu, delta, deg);
            break;
        default:
            throw std::runtime_error("Invalid case type in reduceFBIAtDeg");
    }
}

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::setMasterBoundary(int masterIdx, const ST& value) {
    if (masterIdx < 0 || masterIdx >= numMaster_) {
        throw std::runtime_error("Invalid master index in setMasterBoundary");
    }
    masterBoundary_[masterIdx] = value;
    
    // 初始化主积分的cache_和cacheCurrentDeg_
    auto key = makeKey(masterNus_[masterIdx], masterDeltas_[masterIdx]);
    if (cache_.find(key) == cache_.end()) {
        cache_[key] = Series<ST>(targetDeg_);
    }
    // 主积分通过微分方程求解，设为targetDeg_以避免被约化
    cacheCurrentDeg_[key] = targetDeg_;
}

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::setAllMasterBoundary(const ST& value) {
    for (int i = 0; i < numMaster_; ++i) {
        setMasterBoundary(i, value);
    }
}

template<typename RT, typename PT, typename ST>
const Series<ST>& SeriesSolver<RT, PT, ST>::getMasterSeries(int masterIdx) const {
    if (masterIdx < 0 || masterIdx >= numMaster_) {
        throw std::runtime_error("Invalid master index in getMasterSeries");
    }
    auto key = makeKey(masterNus_[masterIdx], masterDeltas_[masterIdx]);
    auto it = cache_.find(key);
    if (it != cache_.end()) {
        return it->second;
    }
    throw std::runtime_error("Master series not found in cache");
}

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::clearCache() {
    cache_.clear();
}

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::printCacheInfo() const {
    std::cout << "\n========== SeriesSolver Cache Info ==========\n";
    std::cout << "Total cached FBI series: " << cache_.size() << "\n";
    std::cout << "Current degree: " << currentDeg_ << " / " << targetDeg_ << "\n";
    std::cout << "============================================\n\n";
}

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::printAllCache() const {
    std::cout << "\n========== All Cached FBI Series ==========\n";
    std::cout << "Total entries: " << cache_.size() << "\n\n";
    
    int count = 0;
    for (const auto& [key, series] : cache_) {
        const auto& nu = key.nu;
        const auto& delta = key.delta;
        
        std::cout << "Entry " << (++count) << ":\n";
        std::cout << "  nu = {";
        for (size_t i = 0; i < nu.size(); ++i) {
            std::cout << nu[i];
            if (i < nu.size() - 1) std::cout << ",";
        }
        std::cout << "}, delta = " << delta << "\n";
        
        int deg = series.getDeg();
        std::cout << "  Series (deg=" << deg << "):\n";
        for (int d = 0; d <= deg; ++d) {
            std::cout << "    deg=" << d << ": ";
            for (int p = 0; p <= d; ++p) {
                int q = d - p;
                ST coeff = series.getCoeff(p, q);
                if (p > 0) std::cout << ", ";
                std::cout << "(" << p << "," << q << ")=" << coeff;
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }
    std::cout << "============================================\n\n";
}

// ============================================================================
// 辅助函数
// ============================================================================

template<typename RT, typename PT, typename ST>
bool SeriesSolver<RT, PT, ST>::isCorner(const std::vector<int>& nu) const {
    for (int val : nu) {
        if (val != 0 && val != 1) {
            return false;
        }
    }
    return true;
}

template<typename RT, typename PT, typename ST>
std::pair<int, int> SeriesSolver<RT, PT, ST>::findMaxIndex(const std::vector<int>& nu) const {
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
    return std::make_pair(maxIdxInTopSector, maxIdxInCurrentSector);
}

// ============================================================================
// 多项式与级数的卷积（静态方法）
// ============================================================================

template<typename RT, typename PT, typename ST>
ST SeriesSolver<RT, PT, ST>::polySeriesCoeff(const PT& poly, const Series<ST>& series, int p, int q) {
    ST result = ST(0);
    
    // 遍历多项式的所有单项式
    for (auto it = poly.begin(); it != poly.end(); ++it) {
        int a = it->first.x_power;
        int b = it->first.y_power;
        const ST& polyCoeff = it->second;
        
        // 卷积公式: [P·f]_{p,q} = Σ_{a,b} P_{ab} · f_{p-a, q-b}
        if (p >= a && q >= b) {
            result += polyCoeff * series.getCoeff(p - a, q - b);
        }
    }
    
    return result;
}

template<typename RT, typename PT, typename ST>
ST SeriesSolver<RT, PT, ST>::sumPolySeriesCoeff(
    const std::vector<PT>& polys,
    const std::vector<const Series<ST>*>& series,
    int p, int q)
{
    ST result = ST(0);
    
    for (size_t i = 0; i < polys.size() && i < series.size(); ++i) {
        result += polySeriesCoeff(polys[i], *series[i], p, q);
    }
    
    return result;
}

// ============================================================================
// LRR求解器（静态方法）
// ============================================================================

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::solveLRRAtDeg(
    Series<ST>& g,
    const PT& D,
    const std::vector<PT>& polys,
    const std::vector<const Series<ST>*>& series,
    int deg)
{
    ST D00 = D.getCoeff(0, 0);
    
    if (isZero(D00)) {
        throw std::runtime_error("D00 is zero in solveLRRAtDeg");
    }
    
    // 遍历度数为 deg 的所有项 (p, q) 满足 p + q = deg
    for (int p = 0; p <= deg; ++p) {
        int q = deg - p;
        
        // 计算右侧：Σ_i [N_i · f_i]_{p,q}
        ST rhs = sumPolySeriesCoeff(polys, series, p, q);
        
        // 减去 Σ_{(a,b)≠(0,0)} D_{ab} · g_{p-a,q-b}
        for (auto it = D.begin(); it != D.end(); ++it) {
            int a = it->first.x_power;
            int b = it->first.y_power;
            const ST& Dab = it->second;
            
            if ((a != 0 || b != 0) && p >= a && q >= b) {
                rhs -= Dab * g.getCoeff(p - a, q - b);
            }
        }
        
        // g_{p,q} = rhs / D00
        g.setCoeff(p, q, rhs / D00);
    }
}

// ============================================================================
// Case 0 约化实现
// ============================================================================

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::reduceCase0AtDeg(Series<ST>& result, const std::vector<int>& nu, 
                                                  const ST& delta, int deg) {
    // 判断是否为corner（所有元素都是0或1）
    bool corner = isCorner(nu);
    
    if (!corner) {
        // 如果不是corner，执行IBP
        case0IBPAtDeg(result, nu, delta, deg);
    } else {
        // 是corner，检查是否为主积分
        int mfbiIndex = family_.getIndexOfMaster(nu);
        
        if (mfbiIndex >= 0) {
            // 是主积分
            ST targetDelta = masterDeltas_[mfbiIndex];
            
            if (delta == targetDelta) {
                // delta相同，直接从cache_复制主积分系数
                auto masterKey = makeKey(masterNus_[mfbiIndex], masterDeltas_[mfbiIndex]);
                const Series<ST>& masterSeries = cache_.at(masterKey);
                for (int p = 0; p <= deg; ++p) {
                    int q = deg - p;
                    result.setCoeff(p, q, masterSeries.getCoeff(p, q));
                }
            } else {
                // 比较向上和向下移动的距离
                ST upDist = targetDelta - delta;
                ST downDist = delta - targetDelta;
                
                if (upDist <= downDist) {
                    // 向上递推路径更短
                    case0DimShiftUpAtDeg(result, nu, delta, deg);
                } else {
                    // 向下递推路径更短
                    case0DimShiftDownAtDeg(result, nu, delta, deg);
                }
            }
        } else {
            throw std::runtime_error("Corner integral in case0 but not a master FBI");
        }
    }
}

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::case0IBPAtDeg(Series<ST>& result, const std::vector<int>& nuPlus, 
                                               const ST& delta, int deg) {
    // IBP约化公式:
    // ν_max · denoInvS · I_{nuPlus}^Δ = Σ_j numeInvS[row][j] · rhs_j
    // 其中 rhs = (-I_ν^{Δ-1}, ..., +I_{ν-e_α}^{Δ-1}, ...)
    // LRR形式: D = denoInvS * ν_max, N_j = ±numeInvS[row][j], f_j = I...

    auto [maxIndex, maxIndexCur] = findMaxIndex(nuPlus);
    std::vector<int> nu = nuPlus;
    nu[maxIndex]--;
    
    int numBranchCur = family_.nBranch(nu);
    const auto* sector = family_.getSector(nu);
    const auto& denoInvS = sector->getDenoInvS();
    const auto& numeInvS = sector->getNumeInvS();
    
    int row = numBranchCur + maxIndexCur;
    
    // 分母多项式 D = denoInvS[row] * ν_max
    PT D = denoInvS[row] * ST(nu[maxIndex]);
    
    // 构造分子多项式和级数指针列表
    std::vector<PT> polys;
    std::vector<const Series<ST>*> seriesPtrs;
    
    // 获取 I_{nu}^{delta-1}
    const Series<ST>& seriesNuDeltaMinus1 = getFBISeries(nu, delta - ST(1));
    
    // 前B项: -numeInvS[row][j] * I_nu^{delta-1}
    for (int j = 0; j < numBranchCur; ++j) {
        polys.push_back(numeInvS[row][j] * ST(-1));
        seriesPtrs.push_back(&seriesNuDeltaMinus1);
    }
    
    // 后N项: +numeInvS[row][numBranch+iCur] * I_{nu-e_i}^{delta-1}
    int iCur = -1;
    for (int i = 0; i < numProps_; ++i) {
        if (nu[i] > 0) {
            iCur++;
            std::vector<int> nuMinusEi = nu;
            nuMinusEi[i]--;
            
            if (family_.nBranch(nuMinusEi) != numBranchCur) continue;
            
            polys.push_back(numeInvS[row][numBranchCur + iCur]);
            seriesPtrs.push_back(&getFBISeries(nuMinusEi, delta - ST(1)));
        }
    }

    solveLRRAtDeg(result, D, polys, seriesPtrs, deg);
}

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::case0DimShiftDownAtDeg(Series<ST>& result, const std::vector<int>& nu, 
                                                        const ST& delta, int deg) {
    // 向下维度迁移公式:
    // (2Δ-ν-B)·z_0 · I_ν^Δ = C·I_ν^{Δ-1} - Σ_α z_α·I_{ν-e_α}^{Δ-1}
    // LRR形式: D = denoC * (2Δ-ν-B) * z_0, N_0 = numeC, N_i = -numeZ_i
    
    int nuSum = family_.nuSum(nu);
    int numBranchCur = family_.nBranch(nu);
    const auto* sector = family_.getSector(nu);
    
    const PT& denoC = sector->getDenoCandZ();
    const PT& numeC = sector->getNumeC();
    int z0 = sector->getZ0();
    
    // 计算因子 (2Δ - ν - B) * z0
    ST factor = ST(z0) * (ST(2) * delta - ST(nuSum) - ST(numBranchCur));
    
    // D = denoC * factor
    PT D = denoC * factor;
    
    // 构造分子多项式和级数指针列表
    std::vector<PT> polys;
    std::vector<const Series<ST>*> seriesPtrs;
    
    // 第一项: numeC * I_ν^{Δ-1}
    const Series<ST>& seriesDeltaMinus1 = getFBISeries(nu, delta - ST(1));
    polys.push_back(numeC);
    seriesPtrs.push_back(&seriesDeltaMinus1);
    
    // 后续项: -numeZ_i * I_{ν-e_i}^{Δ-1}
    int iCur = -1;
    for (int i = 0; i < numProps_; ++i) {
        if (nu[i] > 0) {
            iCur++;
            std::vector<int> nuMinusEi = nu;
            nuMinusEi[i]--;
            
            if (family_.nBranch(nuMinusEi) != numBranchCur) continue;
            
            polys.push_back(sector->getNumeZ(iCur) * ST(-1));
            seriesPtrs.push_back(&getFBISeries(nuMinusEi, delta - ST(1)));
        }
    }
    
    solveLRRAtDeg(result, D, polys, seriesPtrs, deg);
}

// ============================================================================
// Case 1 约化实现
// ============================================================================

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::reduceCase1AtDeg(Series<ST>& result, const std::vector<int>& nu, 
                                                  const ST& delta, int deg) {
    // Case 1公式 (dimNull=0, C=0):
    // (2Δ - ν - B) · I_ν^Δ = -Σ_α z_α · I_{ν-e_α}^{Δ-1}
    // 
    // LRR形式: D · g = Σ_i N_i · f_i
    // 其中: D = denoC * (2Δ - ν - B) （作为常数多项式）
    //       N_i = -numeZ_i
    //       f_i = I_{ν-e_α}^{Δ-1}
    
    int nuSum = family_.nuSum(nu);
    int numBranchCur = family_.nBranch(nu);
    const auto* sector = family_.getSector(nu);
    
    // 计算因子 (2Δ - ν - B)
    ST factor = ST(2) * delta - ST(nuSum) - ST(numBranchCur);
    
    const PT& denoC = sector->getDenoCandZ();
    
    // 构造分母多项式 D = denoC * factor
    PT D = denoC * factor;
    
    // 构造分子多项式和级数指针列表
    std::vector<PT> polys;
    std::vector<const Series<ST>*> seriesPtrs;
    
    int iCur = -1;
    for (int i = 0; i < numProps_; ++i) {
        if (nu[i] > 0) {
            iCur++;
            std::vector<int> nuMinusEi = nu;
            nuMinusEi[i]--;
            
            if (family_.nBranch(nuMinusEi) != numBranchCur) {
                continue;
            }
            
            // N_i = -numeZ_i
            PT Ni = sector->getNumeZ(iCur) * ST(-1);
            polys.push_back(Ni);
            
            // f_i = I_{ν-e_α}^{Δ-1}
            const Series<ST>& seriesNuMinusEi = getFBISeries(nuMinusEi, delta - ST(1));
            seriesPtrs.push_back(&seriesNuMinusEi);
        }
    }
    
    // 使用LRR求解器
    solveLRRAtDeg(result, D, polys, seriesPtrs, deg);
}

// ============================================================================
// Case 2 约化实现
// ============================================================================

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::reduceCase2AtDeg(Series<ST>& result, const std::vector<int>& nu, 
                                                  const ST& delta, int deg) {
    // Case 2公式 (dimNull>0, C≠0):
    // numeC · I_ν^Δ = Σ_α numeZ_α · I_{ν-e_α}^Δ
    // LRR形式: D = numeC, N_i = numeZ_i, f_i = I_{ν-e_α}^Δ
    
    int numBranchCur = family_.nBranch(nu);
    const auto* sector = family_.getSector(nu);
    const PT& numeC = sector->getNumeC();
    
    // 构造分子多项式和级数指针列表
    std::vector<PT> polys;
    std::vector<const Series<ST>*> seriesPtrs;
    
    int iCur = -1;
    for (int i = 0; i < numProps_; ++i) {
        if (nu[i] > 0) {
            iCur++;
            std::vector<int> nuMinusEi = nu;
            nuMinusEi[i]--;
            
            if (family_.nBranch(nuMinusEi) != numBranchCur) continue;
            
            polys.push_back(sector->getNumeZ(iCur));
            seriesPtrs.push_back(&getFBISeries(nuMinusEi, delta));
        }
    }
    
    solveLRRAtDeg(result, numeC, polys, seriesPtrs, deg);
}

// ============================================================================
// Case 3 约化实现
// ============================================================================

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::reduceCase3AtDeg(Series<ST>& result, const std::vector<int>& nu, 
                                                  const ST& delta, int deg) {
    // Case 3公式 (dimNull>0, C=0):
    // z_β · I_ν^Δ = -Σ_{α≠β} z_α · I_{ν+e_β-e_α}^Δ
    // LRR形式: D = numeZ_β, N_i = -numeZ_α, f_i = I_{ν+e_β-e_α}^Δ
    
    int numBranchCur = family_.nBranch(nu);
    const auto* sector = family_.getSector(nu);
    
    // 找到非零的 z_β
    int jCur = -1;
    int j = -1;
    PT numeZ_j;
    
    int iCurTemp = -1;
    for (int i = 0; i < numProps_; ++i) {
        if (nu[i] > 0) {
            iCurTemp++;
            const PT& numeZ_i = sector->getNumeZ(iCurTemp);
            if (!isZero(numeZ_i)) {
                jCur = iCurTemp;
                j = i;
                numeZ_j = numeZ_i;
                break;
            }
        }
    }
    
    if (j < 0) {
        for (int p = 0; p <= deg; ++p) {
            result.setCoeff(p, deg - p, ST(0));
        }
        return;
    }
    
    // 构造分子多项式和级数指针列表
    std::vector<PT> polys;
    std::vector<const Series<ST>*> seriesPtrs;
    
    int iCur = -1;
    for (int i = 0; i < numProps_; ++i) {
        if (nu[i] > 0) {
            iCur++;
            if (iCur == jCur) continue;  // 跳过 β
            
            std::vector<int> nuShifted = nu;
            nuShifted[i]--;
            nuShifted[j]++;
            
            if (family_.nBranch(nuShifted) != numBranchCur) continue;
            
            // N_i = -numeZ_α
            polys.push_back(sector->getNumeZ(iCur) * ST(-1));
            seriesPtrs.push_back(&getFBISeries(nuShifted, delta));
        }
    }
    
    solveLRRAtDeg(result, numeZ_j, polys, seriesPtrs, deg);
}

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::case0DimShiftUpAtDeg(Series<ST>& result, const std::vector<int>& nu, 
                                                      const ST& delta, int deg) {
    // 向上维度迁移公式:
    // C · I_ν^Δ = (2(Δ+1)-ν-B)·z_0·I_ν^{Δ+1} + Σ_α z_α·I_{ν-e_α}^Δ
    // 
    // 由于使用 numeC = C * denoCandZ, numeZ = z * denoCandZ，需要统一分母：
    // numeC · I_ν^Δ = (2(Δ+1)-ν-B)·z_0·denoCandZ·I_ν^{Δ+1} + Σ_α numeZ_α·I_{ν-e_α}^Δ
    // 
    // LRR形式: D = numeC, N_0 = (2(Δ+1)-ν-B)*z_0*denoCandZ (常数乘以多项式), N_i = numeZ_i
    
    int nuSum = family_.nuSum(nu);
    int numBranchCur = family_.nBranch(nu);
    const auto* sector = family_.getSector(nu);
    
    const PT& numeC = sector->getNumeC();
    const PT& denoCandZ = sector->getDenoCandZ();
    int z0 = sector->getZ0();
    
    // 计算因子 (2(Δ+1) - ν - B) * z0
    ST factor = ST(z0) * (ST(2) * (delta + ST(1)) - ST(nuSum) - ST(numBranchCur));
    
    // 构造分子多项式和级数指针列表
    std::vector<PT> polys;
    std::vector<const Series<ST>*> seriesPtrs;
    
    // 第一项: factor * denoCandZ * I_ν^{Δ+1} (factor乘以denoCandZ多项式)
    const Series<ST>& seriesDeltaPlus1 = getFBISeries(nu, delta + ST(1));
    PT factorPoly = denoCandZ * factor;  // 乘以denoCandZ以统一分母
    polys.push_back(factorPoly);
    seriesPtrs.push_back(&seriesDeltaPlus1);
    
    // 后续项: numeZ_i * I_{ν-e_i}^Δ
    int iCur = -1;
    for (int i = 0; i < numProps_; ++i) {
        if (nu[i] > 0) {
            iCur++;
            std::vector<int> nuMinusEi = nu;
            nuMinusEi[i]--;
            
            if (family_.nBranch(nuMinusEi) != numBranchCur) continue;
            
            polys.push_back(sector->getNumeZ(iCur));
            seriesPtrs.push_back(&getFBISeries(nuMinusEi, delta));
        }
    }
    
    solveLRRAtDeg(result, numeC, polys, seriesPtrs, deg);
}

