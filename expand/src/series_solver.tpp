/**
 * SeriesSolver 实现文件
 * 统一处理微分方程求解和IBP约化，直接计算重定义后FBI的二维幂级数展开
 * 
 * 所有 FBI 均为重定义后的 Ĩ = U^{pow_U} · I
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
    
    masterNus_.reserve(numMaster_);
    masterDeltas_ = family_.getMasterDeltas();
    
    for (int idx : family_.getMasterIdxs()) {
        masterNus_.push_back(family_.secvecFromIdx(idx));
    }
    
    masterBoundary_.resize(numMaster_, ST(0));
}

// ============================================================================
// 主求解函数
// ============================================================================

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::solve() {
    if (!redef_) {
        throw std::runtime_error("Redefinition must be set before calling solve()");
    }
    
    // 初始化主积分零阶系数
    for (int k = 0; k < numMaster_; ++k) {
        auto key = makeKey(masterNus_[k], masterDeltas_[k]);
        cache_[key].setCoeff(0, 0, masterBoundary_[k]);
        cacheCurrentDeg_[key] = 0;
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
    for (int p = 0; p <= deg; ++p) {
        int q = deg - p;
        for (int k = 0; k < numMaster_; ++k) {
            if (p > 0) {
                solveMasterCoeffX(k, p, q);
            } else if (q > 0) {
                solveMasterCoeffY(k, q);
            }
        }
    }
}

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::solveMasterAtDeg(int masterIdx, int deg) {
    if (masterIdx < 0 || masterIdx >= numMaster_) {
        throw std::runtime_error("Invalid master index in solveMasterAtDeg");
    }
    if (deg < 0 || deg > targetDeg_) {
        throw std::runtime_error("Invalid degree in solveMasterAtDeg");
    }

    const auto key = makeKey(masterNus_[masterIdx], masterDeltas_[masterIdx]);
    if (cache_.find(key) == cache_.end()) {
        cache_[key] = Series<ST>(targetDeg_);
    }

    int cachedDeg = -1;
    auto itDeg = cacheCurrentDeg_.find(key);
    if (itDeg != cacheCurrentDeg_.end()) {
        cachedDeg = itDeg->second;
    }
    if (cachedDeg >= deg) {
        return;
    }

    if (cachedDeg < 0) {
        cache_[key].setCoeff(0, 0, masterBoundary_[masterIdx]);
        cachedDeg = 0;
        cacheCurrentDeg_[key] = 0;
    }
    if (deg == 0) {
        return;
    }

    int prevCurrentDeg = currentDeg_;
    for (int d = cachedDeg + 1; d <= deg; ++d) {
        currentDeg_ = d;
        for (int p = 0; p <= d; ++p) {
            int q = d - p;
            if (p > 0) {
                solveMasterCoeffX(masterIdx, p, q);
            } else if (q > 0) {
                solveMasterCoeffY(masterIdx, q);
            }
        }
        cacheCurrentDeg_[key] = d;
    }
    currentDeg_ = prevCurrentDeg;
}

// ============================================================================
// 微分方程求解（重定义后）
// ============================================================================

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::solveMasterCoeffX(int masterIdx, int p, int q) {
    // 重定义后的微分方程:
    // U · ∂_X Ĩ = pow_U · (∂_X U) · Ĩ + Σ (-1/2) · (dR/dX · U^L) · factor · Ĩ_source
    // 
    // 递推:
    // Ĩ_{p,q} = (rhsCoeff + dlogCoeff - lhsCorr) / (U_00 · p)
    
    const std::vector<int>& nu = masterNus_[masterIdx];
    const ST& delta = masterDeltas_[masterIdx];
    const auto& dRdX = family_.getDRdX();
    
    // 计算 rhs: Σ (-1/2) · factor · [dRdXModified · Ĩ_source]_{p-1,q}
    ST rhsCoeff = ST(0);
    for (int i = 0; i < numProps_; ++i) {
        if (nu[i] <= 0) continue;
        for (int j = 0; j < numProps_; ++j) {
            if (nu[j] <= 0) continue;
            ST factor_ij;
            if (i != j) {
                factor_ij = ST(nu[i]) * ST(nu[j]);
            } else {
                factor_ij = ST(nu[i]) * ST(nu[i] + 1);
            }
            
            if (isZero(factor_ij) || isZero(dRdX[i][j])) continue;
            
            std::vector<int> nuShifted = nu;
            nuShifted[i]++;
            nuShifted[j]++;
            
            const Series<ST>& seriesShifted = getFBISeries(nuShifted, delta + ST(1), p + q - 1);
            ST convCoeff = polySeriesCoeff(dRdXModified_[i][j], seriesShifted, p - 1, q);
            rhsCoeff -= factor_ij * convCoeff / ST(2);
        }
    }
    
    auto key = makeKey(nu, delta);
    const Series<ST>& masterSeries = cache_.at(key);
    ST powU = masterPowU_[masterIdx];
    
    // dlogCoeff = +pow_U * [dUdX · Ĩ]_{p-1,q}
    ST dlogCoeff = powU * polySeriesCoeff(redef_->dUdX, masterSeries, p - 1, q);
    
    // lhsCorrection = Σ_{(a,b)≠(0,0)} U_{ab} · (p-a) · Ĩ_{p-a, q-b}
    ST lhsCorr = ST(0);
    for (auto it = redef_->shiftedU.begin(); it != redef_->shiftedU.end(); ++it) {
        int a = it->first.x_power;
        int b = it->first.y_power;
        if (a == 0 && b == 0) continue;
        if (p >= a && q >= b) {
            lhsCorr += it->second * ST(p - a) * masterSeries.getCoeff(p - a, q - b);
        }
    }
    
    ST U00 = redef_->shiftedU.getCoeff(0, 0);
    cache_[key].setCoeff(p, q, (rhsCoeff + dlogCoeff - lhsCorr) / (U00 * ST(p)));
}

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::solveMasterCoeffY(int masterIdx, int q) {
    // 重定义后的微分方程 (Y 方向):
    // U · ∂_Y Ĩ = pow_U · (∂_Y U) · Ĩ + Σ (-1/2) · (dR/dY · U^L) · factor · Ĩ_source
    // 
    // 递推 (p=0):
    // Ĩ_{0,q} = (rhsCoeff + dlogCoeff - lhsCorr) / (U_00 · q)
    
    const std::vector<int>& nu = masterNus_[masterIdx];
    const ST& delta = masterDeltas_[masterIdx];
    const auto& dRdY = family_.getDRdY();
    
    // 计算 rhs: Σ (-1/2) · factor · [dRdYModified · Ĩ_source]_{0,q-1}
    ST rhsCoeff = ST(0);
    for (int i = 0; i < numProps_; ++i) {
        if (nu[i] <= 0) continue;
        for (int j = 0; j < numProps_; ++j) {
            if (nu[j] <= 0) continue;
            ST factor_ij;
            if (i != j) {
                factor_ij = ST(nu[i]) * ST(nu[j]);
            } else {
                factor_ij = ST(nu[i]) * ST(nu[i] + 1);
            }
            
            if (isZero(factor_ij) || isZero(dRdY[i][j])) continue;
            
            std::vector<int> nuShifted = nu;
            nuShifted[i]++;
            nuShifted[j]++;
            
            const Series<ST>& seriesShifted = getFBISeries(nuShifted, delta + ST(1), q - 1);
            ST convCoeff = polySeriesCoeff(dRdYModified_[i][j], seriesShifted, 0, q - 1);
            rhsCoeff -= factor_ij * convCoeff / ST(2);
        }
    }
    
    auto key = makeKey(nu, delta);
    const Series<ST>& masterSeries = cache_.at(key);
    ST powU = masterPowU_[masterIdx];
    
    // dlogCoeff = +pow_U * [dUdY · Ĩ]_{0,q-1}
    ST dlogCoeff = powU * polySeriesCoeff(redef_->dUdY, masterSeries, 0, q - 1);
    
    // lhsCorrection = Σ_{b>0} U_{0,b} · (q-b) · Ĩ_{0, q-b}
    ST lhsCorr = ST(0);
    for (auto it = redef_->shiftedU.begin(); it != redef_->shiftedU.end(); ++it) {
        int a = it->first.x_power;
        int b = it->first.y_power;
        if (a != 0 || b == 0) continue;
        if (q >= b) {
            lhsCorr += it->second * ST(q - b) * masterSeries.getCoeff(0, q - b);
        }
    }
    
    ST U00 = redef_->shiftedU.getCoeff(0, 0);
    cache_[key].setCoeff(0, q, (rhsCoeff + dlogCoeff - lhsCorr) / (U00 * ST(q)));
}

// ============================================================================
// FBI级数获取（递归调用，触发约化）
// ============================================================================

template<typename RT, typename PT, typename ST>
const Series<ST>& SeriesSolver<RT, PT, ST>::getFBISeries(
    const std::vector<int>& nu, const ST& delta, int needDeg) {
    const int targetNeedDeg = std::min(
        targetDeg_, (needDeg >= 0 ? needDeg : std::max(currentDeg_, 0)));
    auto key = makeKey(nu, delta);

    if (cache_.find(key) == cache_.end()) {
        cache_[key] = Series<ST>(targetDeg_);
        cacheCurrentDeg_[key] = -1;
    } else if (cacheCurrentDeg_.find(key) == cacheCurrentDeg_.end()) {
        cacheCurrentDeg_[key] = -1;
    }

    int cachedDeg = cacheCurrentDeg_[key];
    if (cachedDeg >= targetNeedDeg) {
        return cache_[key];
    }

    bool isMasterKey = false;
    int masterIdx = -1;
    if (family_.isMaster(nu)) {
        masterIdx = family_.getIndexOfMaster(nu);
        isMasterKey = (masterIdx >= 0 && delta == masterDeltas_[masterIdx]);
    }

    for (int deg = cachedDeg + 1; deg <= targetNeedDeg; ++deg) {
        if (isMasterKey) {
            solveMasterAtDeg(masterIdx, deg);
        } else {
            reduceFBIAtDeg(cache_[key], nu, delta, deg);
            cacheCurrentDeg_[key] = deg;
        }
    }
    return cache_[key];
}

// ============================================================================
// 约化主逻辑
// ============================================================================

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::reduceFBI(Series<ST>& result, const std::vector<int>& nu, const ST& delta) {
    result.setCoeff(0, 0, ST(0));
    for (int deg = 0; deg <= currentDeg_ && deg <= targetDeg_; ++deg) {
        reduceFBIAtDeg(result, nu, delta, deg);
    }
}

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::reduceFBIAtDeg(Series<ST>& result, const std::vector<int>& nu, 
                                                const ST& delta, int deg) {
    int caseType = family_.getCase(nu);
    switch (caseType) {
        case 0: reduceCase0AtDeg(result, nu, delta, deg); break;
        case 1: reduceCase1AtDeg(result, nu, delta, deg); break;
        case 2: reduceCase2AtDeg(result, nu, delta, deg); break;
        case 3: reduceCase3AtDeg(result, nu, delta, deg); break;
        default: throw std::runtime_error("Invalid case type in reduceFBIAtDeg");
    }
}

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::setMasterBoundary(int masterIdx, const ST& value) {
    if (masterIdx < 0 || masterIdx >= numMaster_) {
        throw std::runtime_error("Invalid master index in setMasterBoundary");
    }
    masterBoundary_[masterIdx] = value;
    
    auto key = makeKey(masterNus_[masterIdx], masterDeltas_[masterIdx]);
    if (cache_.find(key) == cache_.end()) {
        cache_[key] = Series<ST>(targetDeg_);
    }
    cache_[key].setCoeff(0, 0, value);
    cacheCurrentDeg_[key] = 0;
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
    cacheCurrentDeg_.clear();
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
        if (val != 0 && val != 1) return false;
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
        if (nu[i] > 0) maxIdxInCurrentSector++;
    }
    return std::make_pair(maxIdxInTopSector, maxIdxInCurrentSector);
}

// ============================================================================
// 多项式与级数的卷积（静态方法）
// ============================================================================

template<typename RT, typename PT, typename ST>
ST SeriesSolver<RT, PT, ST>::polySeriesCoeff(const PT& poly, const Series<ST>& series, int p, int q) {
    ST result = ST(0);
    for (auto it = poly.begin(); it != poly.end(); ++it) {
        int a = it->first.x_power;
        int b = it->first.y_power;
        const ST& polyCoeff = it->second;
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
    
    for (int p = 0; p <= deg; ++p) {
        int q = deg - p;
        ST rhs = sumPolySeriesCoeff(polys, series, p, q);
        
        for (auto it = D.begin(); it != D.end(); ++it) {
            int a = it->first.x_power;
            int b = it->first.y_power;
            const ST& Dab = it->second;
            if ((a != 0 || b != 0) && p >= a && q >= b) {
                rhs -= Dab * g.getCoeff(p - a, q - b);
            }
        }
        
        g.setCoeff(p, q, rhs / D00);
    }
}

// ============================================================================
// Case 0 约化实现
// ============================================================================

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::reduceCase0AtDeg(Series<ST>& result, const std::vector<int>& nu, 
                                                  const ST& delta, int deg) {
    if (!isCorner(nu)) {
        case0IBPAtDeg(result, nu, delta, deg);
        return;
    }
    int mfbiIndex = family_.getIndexOfMaster(nu);
    if (mfbiIndex < 0) {
        throw std::runtime_error("Corner integral in case0 but not a master FBI");
    }

    ST targetDelta = masterDeltas_[mfbiIndex];
    if (delta == targetDelta) {
        auto masterKey = makeKey(masterNus_[mfbiIndex], masterDeltas_[mfbiIndex]);
        if (cache_.find(masterKey) == cache_.end() ||
            cacheCurrentDeg_.find(masterKey) == cacheCurrentDeg_.end() ||
            cacheCurrentDeg_[masterKey] < deg) {
            const Series<ST>& ensured = getFBISeries(masterNus_[mfbiIndex], masterDeltas_[mfbiIndex], deg);
            (void)ensured;
        }
        const Series<ST>& masterSeries = cache_.at(masterKey);
        for (int p = 0; p <= deg; ++p) {
            result.setCoeff(p, deg - p, masterSeries.getCoeff(p, deg - p));
        }
    } else {
        ST upDist = targetDelta - delta;
        ST downDist = delta - targetDelta;
        if (upDist <= downDist) {
            case0DimShiftUpAtDeg(result, nu, delta, deg);
        } else {
            case0DimShiftDownAtDeg(result, nu, delta, deg);
        }
    }
}

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::case0IBPAtDeg(Series<ST>& result, const std::vector<int>& nuPlus, 
                                               const ST& delta, int deg) {
    auto [maxIndex, maxIndexCur] = findMaxIndex(nuPlus);
    std::vector<int> nu = nuPlus;
    nu[maxIndex]--;
    
    int numBranchCur = family_.nBranch(nu);
    const auto* sector = family_.getSector(nu);
    const auto& denoInvS = sector->getDenoInvS();
    const auto& numeInvS = sector->getNumeInvS();
    
    int row = numBranchCur + maxIndexCur;
    PT D = denoInvS[row] * ST(nu[maxIndex]);
    
    std::vector<PT> polys;
    std::vector<const Series<ST>*> seriesPtrs;
    
    const Series<ST>& seriesNuDeltaMinus1 = getFBISeries(nu, delta - ST(1), deg);
    for (int j = 0; j < numBranchCur; ++j) {
        polys.push_back(numeInvS[row][j] * ST(-1));
        seriesPtrs.push_back(&seriesNuDeltaMinus1);
    }
    
    int iCur = -1;
    for (int i = 0; i < numProps_; ++i) {
        if (nu[i] > 0) {
            iCur++;
            std::vector<int> nuMinusEi = nu;
            nuMinusEi[i]--;
            if (family_.nBranch(nuMinusEi) != numBranchCur) continue;
            polys.push_back(numeInvS[row][numBranchCur + iCur]);
            seriesPtrs.push_back(&getFBISeries(nuMinusEi, delta - ST(1), deg));
        }
    }

    // 应用 ratio 因子
    int nuTotT = nuTotSum(nuPlus);
    int nuTotNu = nuTotSum(nu);
    ST deltaMinus1 = delta - ST(1);
    std::vector<int> deltaPs;
    for (int j = 0; j < numBranchCur; ++j)
        deltaPs.push_back(redef_->deltaPowU(nuTotT, delta, nuTotNu, deltaMinus1));
    int iC2 = -1;
    for (int i = 0; i < numProps_; ++i) {
        if (nu[i] > 0) {
            iC2++;
            std::vector<int> nme = nu; nme[i]--;
            if (family_.nBranch(nme) != numBranchCur) continue;
            deltaPs.push_back(redef_->deltaPowU(nuTotT, delta, nuTotNu - 1, deltaMinus1));
        }
    }
    applyRatioFactors(D, polys, deltaPs);

    solveLRRAtDeg(result, D, polys, seriesPtrs, deg);
}

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::case0DimShiftDownAtDeg(Series<ST>& result, const std::vector<int>& nu, 
                                                        const ST& delta, int deg) {
    int nuSum = family_.nuSum(nu);
    int numBranchCur = family_.nBranch(nu);
    const auto* sector = family_.getSector(nu);
    
    const PT& denoC = sector->getDenoCandZ();
    const PT& numeC = sector->getNumeC();
    int z0 = sector->getZ0();
    
    ST factor = ST(z0) * (ST(2) * delta - ST(nuSum) - ST(numBranchCur));
    PT D = denoC * factor;
    
    std::vector<PT> polys;
    std::vector<const Series<ST>*> seriesPtrs;
    
    const Series<ST>& seriesDeltaMinus1 = getFBISeries(nu, delta - ST(1), deg);
    polys.push_back(numeC);
    seriesPtrs.push_back(&seriesDeltaMinus1);
    
    int iCur = -1;
    for (int i = 0; i < numProps_; ++i) {
        if (nu[i] > 0) {
            iCur++;
            std::vector<int> nuMinusEi = nu;
            nuMinusEi[i]--;
            if (family_.nBranch(nuMinusEi) != numBranchCur) continue;
            polys.push_back(sector->getNumeZ(iCur) * ST(-1));
            seriesPtrs.push_back(&getFBISeries(nuMinusEi, delta - ST(1), deg));
        }
    }
    
    // 应用 ratio 因子
    int nuTotT = nuTotSum(nu);
    ST deltaMinus1 = delta - ST(1);
    std::vector<int> deltaPs;
    deltaPs.push_back(redef_->deltaPowU(nuTotT, delta, nuTotT, deltaMinus1));
    int iC2 = -1;
    for (int i = 0; i < numProps_; ++i) {
        if (nu[i] > 0) {
            iC2++;
            std::vector<int> nme = nu; nme[i]--;
            if (family_.nBranch(nme) != numBranchCur) continue;
            deltaPs.push_back(redef_->deltaPowU(nuTotT, delta, nuTotT - 1, deltaMinus1));
        }
    }
    applyRatioFactors(D, polys, deltaPs);
    
    solveLRRAtDeg(result, D, polys, seriesPtrs, deg);
}

// ============================================================================
// Case 1 约化实现
// ============================================================================

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::reduceCase1AtDeg(Series<ST>& result, const std::vector<int>& nu, 
                                                  const ST& delta, int deg) {
    int nuSum = family_.nuSum(nu);
    int numBranchCur = family_.nBranch(nu);
    const auto* sector = family_.getSector(nu);
    
    ST factor = ST(2) * delta - ST(nuSum) - ST(numBranchCur);
    const PT& denoC = sector->getDenoCandZ();
    PT D = denoC * factor;
    
    std::vector<PT> polys;
    std::vector<const Series<ST>*> seriesPtrs;
    
    int iCur = -1;
    for (int i = 0; i < numProps_; ++i) {
        if (nu[i] > 0) {
            iCur++;
            std::vector<int> nuMinusEi = nu;
            nuMinusEi[i]--;
            if (family_.nBranch(nuMinusEi) != numBranchCur) continue;
            polys.push_back(sector->getNumeZ(iCur) * ST(-1));
            seriesPtrs.push_back(&getFBISeries(nuMinusEi, delta - ST(1), deg));
        }
    }
    
    // 应用 ratio 因子
    int nuTotT = nuTotSum(nu);
    ST deltaMinus1 = delta - ST(1);
    std::vector<int> deltaPs;
    for (size_t idx = 0; idx < polys.size(); ++idx)
        deltaPs.push_back(redef_->deltaPowU(nuTotT, delta, nuTotT - 1, deltaMinus1));
    applyRatioFactors(D, polys, deltaPs);
    
    solveLRRAtDeg(result, D, polys, seriesPtrs, deg);
}

// ============================================================================
// Case 2 约化实现
// ============================================================================

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::reduceCase2AtDeg(Series<ST>& result, const std::vector<int>& nu, 
                                                  const ST& delta, int deg) {
    int numBranchCur = family_.nBranch(nu);
    const auto* sector = family_.getSector(nu);
    const PT& numeC = sector->getNumeC();
    
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
            seriesPtrs.push_back(&getFBISeries(nuMinusEi, delta, deg));
        }
    }
    
    // 应用 ratio 因子
    PT D = numeC;
    int nuTotT = nuTotSum(nu);
    std::vector<int> deltaPs;
    for (size_t idx = 0; idx < polys.size(); ++idx)
        deltaPs.push_back(redef_->deltaPowU(nuTotT, delta, nuTotT - 1, delta));
    applyRatioFactors(D, polys, deltaPs);
    
    solveLRRAtDeg(result, D, polys, seriesPtrs, deg);
}

// ============================================================================
// Case 3 约化实现
// ============================================================================

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::reduceCase3AtDeg(Series<ST>& result, const std::vector<int>& nu, 
                                                  const ST& delta, int deg) {
    int numBranchCur = family_.nBranch(nu);
    const auto* sector = family_.getSector(nu);
    
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
    
    std::vector<PT> polys;
    std::vector<const Series<ST>*> seriesPtrs;
    
    int iCur = -1;
    for (int i = 0; i < numProps_; ++i) {
        if (nu[i] > 0) {
            iCur++;
            if (iCur == jCur) continue;
            
            std::vector<int> nuShifted = nu;
            nuShifted[i]--;
            nuShifted[j]++;
            if (family_.nBranch(nuShifted) != numBranchCur) continue;
            
            polys.push_back(sector->getNumeZ(iCur) * ST(-1));
            seriesPtrs.push_back(&getFBISeries(nuShifted, delta, deg));
        }
    }
    
    solveLRRAtDeg(result, numeZ_j, polys, seriesPtrs, deg);
}

// ============================================================================
// 多项式辅助函数 + Redefinition 支持
// ============================================================================

template<typename RT, typename PT, typename ST>
PT SeriesSolver<RT, PT, ST>::multiplyPolys(const PT& a, const PT& b) {
    PT result;
    for (const auto& [pa, ca] : a) {
        for (const auto& [pb, cb] : b) {
            result.addMonomial(ca * cb, Power(pa.x_power + pb.x_power, pa.y_power + pb.y_power));
        }
    }
    return result;
}

template<typename RT, typename PT, typename ST>
PT SeriesSolver<RT, PT, ST>::powPolyExpand(const PT& base, int exp) {
    PT result;
    result.addMonomial(ST(1), Power(0, 0));
    if (exp <= 0) return result;
    PT cur = base;
    int e = exp;
    while (e > 0) {
        if (e & 1) result = multiplyPolys(result, cur);
        e >>= 1;
        if (e > 0) cur = multiplyPolys(cur, cur);
    }
    return result;
}

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::applyRatioFactors(
    PT& D,
    std::vector<PT>& polys,
    const std::vector<int>& deltaPs) const
{
    if (deltaPs.empty()) return;

    // Ĩ = U^{pow_U} · I
    // 重定义后: D · Ĩ_T = Σ N_i · U^{Δp_i} · Ĩ_S_i
    // 若有 Δp_i < 0, 两边乘以 U^m (m = max(0, -min(Δp_i))):
    // D·U^m · Ĩ_T = Σ (N_i · U^{Δp_i + m}) · Ĩ_S_i
    int dpMin = *std::min_element(deltaPs.begin(), deltaPs.end());
    int m = std::max(0, -dpMin);

    if (m > 0) {
        D = multiplyPolys(D, powPolyExpand(redef_->shiftedU, m));
    }

    for (size_t idx = 0; idx < deltaPs.size() && idx < polys.size(); ++idx) {
        int adjPow = deltaPs[idx] + m;
        if (adjPow > 0) {
            polys[idx] = multiplyPolys(polys[idx], powPolyExpand(redef_->shiftedU, adjPow));
        }
    }
}

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::setRedefinition(const Redefinition<PT, ST>* redef) {
    redef_ = redef;
    if (!redef_) return;

    // 预计算 dRdX * U^L 和 dRdY * U^L
    const auto& dRdX = family_.getDRdX();
    const auto& dRdY = family_.getDRdY();
    PT UL = powPolyExpand(redef_->shiftedU, redef_->L);

    dRdXModified_.resize(numProps_);
    dRdYModified_.resize(numProps_);
    for (int i = 0; i < numProps_; ++i) {
        dRdXModified_[i].resize(numProps_);
        dRdYModified_[i].resize(numProps_);
        for (int j = 0; j < numProps_; ++j) {
            if (!isZero(dRdX[i][j]))
                dRdXModified_[i][j] = multiplyPolys(dRdX[i][j], UL);
            if (!isZero(dRdY[i][j]))
                dRdYModified_[i][j] = multiplyPolys(dRdY[i][j], UL);
        }
    }

    // 预计算每个主积分的 pow_U
    masterPowU_.resize(numMaster_);
    for (int k = 0; k < numMaster_; ++k) {
        int nuTot = nuTotSum(masterNus_[k]);
        masterPowU_[k] = redef_->powUScalar(nuTot, masterDeltas_[k]);
    }
}

template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::case0DimShiftUpAtDeg(Series<ST>& result, const std::vector<int>& nu, 
                                                      const ST& delta, int deg) {
    int nuSum = family_.nuSum(nu);
    int numBranchCur = family_.nBranch(nu);
    const auto* sector = family_.getSector(nu);
    
    const PT& numeC = sector->getNumeC();
    const PT& denoCandZ = sector->getDenoCandZ();
    int z0 = sector->getZ0();
    
    ST factor = ST(z0) * (ST(2) * (delta + ST(1)) - ST(nuSum) - ST(numBranchCur));
    
    std::vector<PT> polys;
    std::vector<const Series<ST>*> seriesPtrs;
    
    const Series<ST>& seriesDeltaPlus1 = getFBISeries(nu, delta + ST(1), deg);
    PT factorPoly = denoCandZ * factor;
    polys.push_back(factorPoly);
    seriesPtrs.push_back(&seriesDeltaPlus1);
    
    int iCur = -1;
    for (int i = 0; i < numProps_; ++i) {
        if (nu[i] > 0) {
            iCur++;
            std::vector<int> nuMinusEi = nu;
            nuMinusEi[i]--;
            if (family_.nBranch(nuMinusEi) != numBranchCur) continue;
            polys.push_back(sector->getNumeZ(iCur));
            seriesPtrs.push_back(&getFBISeries(nuMinusEi, delta, deg));
        }
    }
    
    // 应用 ratio 因子
    PT D = numeC;
    int nuTotT = nuTotSum(nu);
    ST deltaPlus1 = delta + ST(1);
    std::vector<int> deltaPs;
    deltaPs.push_back(redef_->deltaPowU(nuTotT, delta, nuTotT, deltaPlus1));
    int iC2 = -1;
    for (int i = 0; i < numProps_; ++i) {
        if (nu[i] > 0) {
            iC2++;
            std::vector<int> nme = nu; nme[i]--;
            if (family_.nBranch(nme) != numBranchCur) continue;
            deltaPs.push_back(redef_->deltaPowU(nuTotT, delta, nuTotT - 1, delta));
        }
    }
    applyRatioFactors(D, polys, deltaPs);
    
    solveLRRAtDeg(result, D, polys, seriesPtrs, deg);
}
