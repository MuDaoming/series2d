template<typename T>
const std::vector<T>& FBIReducer<T>::getReductionCoeff(const std::vector<int>& nu, T delta) {

    // if in cache
    auto key = makeKey(nu, delta);
    auto it = cache_.find(key);
    if (it != cache_.end()) { return it->second; }
    
    // if not
    auto result = reduceFBI(nu, delta);
    cache_[key] = result;
    return cache_[key];
}

template<typename T>
std::vector<T> FBIReducer<T>::reduceFBI(const std::vector<int>& nu, T delta) {
    int caseNum = family_->getCase(nu);
    if (caseNum == 0) { return reduceCase0(nu, delta);} 
    else if (caseNum == 1) { return reduceCase1(nu, delta); }
    else if (caseNum == 2) { return reduceCase2(nu, delta); }
    else if (caseNum == 3) { return reduceCase3(nu, delta); }
    else {
        throw std::runtime_error("Invalid case number in reduceFBI");
    }
}

template<typename T>
std::vector<T> FBIReducer<T>::reduceCase0(const std::vector<int>& nu, T delta) {
    
    // 判断是否为corner（所有元素都是0或1）
    bool corner = isCorner(nu);
    
    if (!corner) {
        // 4.a: 如果不是corner，执行IBP
        return case0IBP(nu, delta);
    }
    
    // 是corner，检查是否为主积分
    int mfbiIndex = family_->getIndexOfMaster(nu);
    
    if (mfbiIndex == -1) {
        // 不是主积分，但是corner，这不应该发生在case0中
        throw std::runtime_error("Corner integral in case0 but not a master FBI");
    }
    
    // 是主积分，获取目标delta
    T targetDelta = family_->getMasterDelta(mfbiIndex);
    
    if (delta == targetDelta) {
        // 4.c: corner且delta相同，返回单位向量
        int numMFBIs = family_->getNumMaster();
        std::vector<T> result(numMFBIs, T(0));
        result[mfbiIndex] = T(1);
        return result;
    }
    
    // 4.b: corner但delta不一致，执行维度递推
    return case0CornerDimensionShift(nu, delta, targetDelta);
}

// Case0 Corner维度递推包装函数
template<typename T>
std::vector<T> FBIReducer<T>::case0CornerDimensionShift(const std::vector<int>& nu, T delta, T targetDelta) {
    // 在有限域下直接比较判断递推方向
    T upDist = targetDelta - delta;
    T downDist = delta - targetDelta;
    
    if (upDist <= downDist) {
        // 向上递推路径更短
        return case0CornerDimensionShiftUp(nu, delta, targetDelta);
    } else {
        // 向下递推路径更短
        return case0CornerDimensionShiftDown(nu, delta, targetDelta);
    }
}

// Corner向上递推：从delta递推到targetDelta (delta < targetDelta)
// I_ν^delta = ((2*(delta+1) - ν - B) z_0 I_ν^(delta+1)  + Σ z_α I_(ν-e_α)^delta) / C
template<typename T>
std::vector<T> FBIReducer<T>::case0CornerDimensionShiftUp(const std::vector<int>& nu, T delta, T targetDelta) {
    int numBranch = family_->getNumBranch();
    int numProps = family_->getNumProps();
    int nuSum = family_->nuSum(nu);
    const auto* sector = family_->getSector(nu);
    T C = sector->getCSum();
    T z0 = sector->getZ0();
    
    // 递归调用 I_nu^{delta+1}（这是向上递推的目标）
    const auto& coeffDeltaPlus1 = getReductionCoeff(nu, delta + T(1));
    
    // 计算 (2*(delta+1) - nu - B) z_0 I(delta+1)
    T factor = z0 * (T(2) * (delta + T(1)) - T(nuSum) - T(numBranch));
    
    std::vector<T> result(coeffDeltaPlus1.size());
    for (size_t i = 0; i < coeffDeltaPlus1.size(); i++) {
        result[i] = factor * coeffDeltaPlus1[i];
    }
    
    // 加上 Σ z_α I_(ν-e_α)^delta 项
    int iCur = -1;
    for (int i = 0; i < numProps; i++) {
        if (nu[i] > 0) {
            iCur++;
            T z_i = sector->getZ(iCur);
            std::vector<int> nuMinusEi = nu;
            nuMinusEi[i]--;
            
            if (family_->nBranch(nuMinusEi) != numBranch) {
                continue;
            }
            
            const auto& coeffNuMinusEi = getReductionCoeff(nuMinusEi, delta);
            for (size_t j = 0; j < coeffNuMinusEi.size() && j < result.size(); j++) {
                result[j] += z_i * coeffNuMinusEi[j];
            }
        }
    }
    
    // 除以C
    for (size_t i = 0; i < result.size(); i++) {
        result[i] = result[i] / C;
    }
    
    return result;
}

// Corner向下递推：从delta递推到targetDelta (delta > targetDelta)
// I_ν^delta = [C I_ν^(delta-1) - Σ z_α I_(ν-e_α)^(delta-1)] / [(2*delta - ν - B) z_0]
template<typename T>
std::vector<T> FBIReducer<T>::case0CornerDimensionShiftDown(const std::vector<int>& nu, T delta, T targetDelta) {
    int numBranch = family_->getNumBranch();
    int numProps = family_->getNumProps();
    int nuSum = family_->nuSum(nu);
    const auto* sector = family_->getSector(nu);
    T C = sector->getCSum();
    T z0 = sector->getZ0();
    
    // 递归调用 I_nu^{delta-1}
    const auto& coeffDeltaMinus1 = getReductionCoeff(nu, delta - T(1));
    
    // 计算 C I_ν^(delta-1)
    std::vector<T> result(coeffDeltaMinus1.size());
    for (size_t i = 0; i < coeffDeltaMinus1.size(); i++) {
        result[i] = C * coeffDeltaMinus1[i];
    }
    
    // 减去 Σ z_α I_(ν-e_α)^(delta-1) 项
    int iCur = -1;
    for (int i = 0; i < numProps; i++) {
        if (nu[i] > 0) {
            iCur++;
            T z_i = sector->getZ(iCur);
            std::vector<int> nuMinusEi = nu;
            nuMinusEi[i]--;
            
            if (family_->nBranch(nuMinusEi) != numBranch) {
                continue;
            }
            
            const auto& coeffNuMinusEi = getReductionCoeff(nuMinusEi, delta - T(1));
            for (size_t j = 0; j < coeffNuMinusEi.size() && j < result.size(); j++) {
                result[j] -= z_i * coeffNuMinusEi[j];
            }
        }
    }
    
    // 除以 (2*delta - ν - B) z_0
    T factor = z0 * (T(2) * delta - T(nuSum) - T(numBranch));
    for (size_t i = 0; i < result.size(); i++) {
        result[i] = result[i] / factor;
    }
    
    return result;
}


template<typename T>
std::vector<T> FBIReducer<T>::case0IBP(const std::vector<int>& nuPlus, T delta) {
    // 构造 nu = nuPlus - e_maxIdx
    auto [maxIndex, maxIndexCur] = findMaxIndex(nuPlus);
    std::vector<int> nu = nuPlus;
    nu[maxIndex]--;

    // 获取 nu 的约化系数
    const auto& coeffNu = getReductionCoeff(nu, delta - T(1));
    
    // 构造右侧向量 rhs
    int numBranch = family_->getNumBranch();
    int numProps = family_->getNumProps();
    int numPropsCur = family_->nProps(nu);
    int nCur = numBranch + numPropsCur;
    int numMaster = family_->getNumMaster();
    std::vector<std::vector<T>> rhs(nCur, std::vector<T>(numMaster));
    
    // rhs 前 B 个元素：-FBI(nu, delta - 1) 的约化系数
    for (int i = 0; i < numBranch; i++) {
        rhs[i] = coeffNu;
        for (auto& x : rhs[i]) x = -x;
    }
    
    // rhs 后 N 个元素：FBI(nu - e_i, delta - 1) 的约化系数
    int iCur = -1;
    for (int i = 0; i < numProps; i++) {
        std::vector<int> nuMinusEi;
        if (nu[i] > 0) {
            iCur++;
            nuMinusEi = nu;
            nuMinusEi[i]--;
            if (family_->nBranch(nuMinusEi) != numBranch) {
                continue;
            }
            const auto& coeffNuMinusEi = getReductionCoeff(nuMinusEi, delta - T(1));
            rhs[numBranch + iCur] = coeffNuMinusEi; 
        }
    }

    // 使用 invS 求解
    auto sector = family_->getSector(nu);
    const auto& invS = sector->getInvS();
    std::vector<T> solution(numMaster, T(0));
    
    for (int i = 0; i < numMaster; i++) {
        for (int j = 0; j < nCur; j++) {
            solution[i] += invS[numBranch+maxIndexCur][j] * rhs[j][i];
        }
        solution[i] /= nu[maxIndex];
    }

    auto key = makeKey(nuPlus, delta);
    cache_[key] = solution;  
    return solution;
}

template<typename T>
std::vector<T> FBIReducer<T>::case0DimensionShift(const std::vector<int>& nu, T delta) {

    int numBranch = family_->getNumBranch();
    int numProps = family_->getNumProps();
    int nuSum = family_->nuSum(nu);
    const auto& coeffDeltaPlus1 = getReductionCoeff(nu, delta + 1);
    const auto* sector = family_->getSector(nu);
    T C = sector->getCSum();
    T z0 = sector->getZ0();
    T constant = z0 * (T(2) * (delta + 1) - T(nuSum) - T(numBranch));


    std::vector<T> result(coeffDeltaPlus1.size());
    for (size_t i = 0; i < coeffDeltaPlus1.size(); i++) {
        result[i] = constant * coeffDeltaPlus1[i];
    }
    

    int iCur = -1;
    for (int i = 0; i < numProps; i++) {
        std::vector<int> nuMinusEi;
        if (nu[i] > 0) {
            iCur ++;
            T z_i = sector->getZ(iCur);
            nuMinusEi = nu;
            nuMinusEi[i]--;
            if (family_->nBranch(nuMinusEi) != numBranch) {
                continue;
            }
            const auto& coeffNuMinusEi = getReductionCoeff(nuMinusEi, delta);
            for (size_t j = 0; j < coeffNuMinusEi.size() && j < result.size(); j++) {
                result[j] += z_i * coeffNuMinusEi[j];
            }
        }
    }
    

    for (size_t i = 0; i < result.size(); i++) {
        result[i] = result[i] / C;
    }
    
    return result;
}


template<typename T>
std::vector<T> FBIReducer<T>::reduceCase1(const std::vector<int>& nu, T delta) {
    // Formula: (2*Delta - nu_sum - B)*I_nu^{Delta} = -sum_{alpha=1}^{N} z_alpha * I_{nu-e_alpha}^{Delta-1}
    int numBranch = family_->getNumBranch();
    int numProps = family_->getNumProps();
    int nuSum = family_->nuSum(nu);
    const auto* sector = family_->getSector(nu);
    
    T factor = T(2) * delta - T(nuSum) - T(numBranch);
    
    std::vector<T> result;
    result.resize(family_->getNumMaster(), T(0));
    
    int iCur = -1;
    for (int i = 0; i < numProps; i++) {
        if (nu[i] > 0) {
            iCur++;
            T z_i = sector->getZ(iCur);
            std::vector<int> nuMinusEi = nu;
            nuMinusEi[i]--;
            
            if (family_->nBranch(nuMinusEi) != numBranch) {
                continue;
            }
            
            const auto& coeffNuMinusEi = getReductionCoeff(nuMinusEi, delta - 1);
            for (size_t j = 0; j < coeffNuMinusEi.size() && j < result.size(); j++) {
                result[j] += (-z_i) * coeffNuMinusEi[j];
            }
        }
    }
    
    for (size_t i = 0; i < result.size(); i++) {
        result[i] = result[i] / factor;
    }
    
    return result;
}

template<typename T>
std::vector<T> FBIReducer<T>::reduceCase2(const std::vector<int>& nu, T delta) {
// I(nu, delta) = 1/C sum_i z_i I(nu - e_i, delta)
    int numBranch = family_->getNumBranch();
    int numProps = family_->getNumProps();
    const auto* sector = family_->getSector(nu);
    T C = sector->getCSum();

    std::vector<T> result;
    result.resize(family_->getNumMaster(), T(0));

    int iCur = -1;
    for (int i = 0; i < numProps; i++) {
        std::vector<int> nuMinusEi;
        if (nu[i] > 0) {
            iCur ++;
            T z_i = sector->getZ(iCur);
            nuMinusEi = nu;
            nuMinusEi[i]--;
            if (family_->nBranch(nuMinusEi) != numBranch) {
                continue;
            }
            const auto& coeffNuMinusEi = getReductionCoeff(nuMinusEi, delta);
            for (size_t j = 0; j < coeffNuMinusEi.size() && j < result.size(); j++) {
                result[j] += z_i * coeffNuMinusEi[j];
            }
        }
    }

    for (size_t i = 0; i < result.size(); i++) {
        result[i] = result[i] / C;
    }

    return result;
}


template<typename T>
std::vector<T> FBIReducer<T>::reduceCase3(const std::vector<int>& nu, T delta) {
// z_j !=0
// I(nu ,delta) = - sum_{i!=j} z_i / z_j I(nu + e_j - e_i, delta)

    int numBranch = family_->getNumBranch();
    int numProps = family_->getNumProps();
    const auto* sector = family_->getSector(nu);

    // 找到非零 z_j
    int iCur = -1;
    int jCur = -1;
    int j;
    T z_j = T(0);
    for (int i = 0; i < numProps; i++) {
        if (nu[i] > 0) {
            iCur++;
            T z_i = sector->getZ(iCur);
            if (z_i != T(0)) {
                jCur = iCur;
                z_j = z_i;
                j = i;
                break;
            }
        }
    }

    std::vector<T> result;
    result.resize(family_->getNumMaster(), T(0));

    iCur = -1;
    for (int i = 0; i < numProps; i++) {
        std::vector<int> nuMinusEiPlusEj;
        if (nu[i] > 0) {
            iCur++;
            if (iCur == jCur) continue; // skip j
            T z_i = sector->getZ(iCur); // iCur + 1 because we skipped j
            nuMinusEiPlusEj = nu;
            nuMinusEiPlusEj[i]--;
            nuMinusEiPlusEj[j]++;
            if (family_->nBranch(nuMinusEiPlusEj) != numBranch) {
                continue;
            }
            const auto& coeffNuMinusEiPlusEj = getReductionCoeff(nuMinusEiPlusEj, delta);
            for (size_t k = 0; k < coeffNuMinusEiPlusEj.size() && k < result.size(); k++) {
                result[k] += (- z_i / z_j) * coeffNuMinusEiPlusEj[k];
        }
    }
    }

    return result;
}

// 清空缓存
template<typename T>
void FBIReducer<T>::clearCache() {
    cache_.clear();
}

// 获取缓存大小
template<typename T>
size_t FBIReducer<T>::getCacheSize() const {
    return cache_.size();
}

// 打印缓存内容
template<typename T>
void FBIReducer<T>::printCache() const {
    std::cout << "\n========== Cache Contents ==========\n";
    std::cout << "Total entries: " << cache_.size() << "\n\n";
    
    int count = 0;
    for (const auto& [key, coeffs] : cache_) {
        const auto& nu = std::get<0>(key);
        const auto& delta = std::get<1>(key);
        
        std::cout << "Entry " << (++count) << ":\n";
        std::cout << "  nu = {";
        for (size_t i = 0; i < nu.size(); i++) {
            std::cout << nu[i];
            if (i < nu.size() - 1) std::cout << ", ";
        }
        std::cout << "}\n";
        std::cout << "  delta = " << delta << "\n";
        std::cout << "  coeffs (size = " << coeffs.size() << "): [";
        for (size_t i = 0; i < coeffs.size(); i++) {
            std::cout << coeffs[i];
            if (i < coeffs.size() - 1) std::cout << ", ";
        }
        std::cout << "]\n\n";
    }
    std::cout << "====================================\n\n";
}
