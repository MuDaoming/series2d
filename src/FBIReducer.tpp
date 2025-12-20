// FBIReducer 类实现文件

// 构造函数
template<typename T>
FBIReducer<T>::FBIReducer(const Family<T>* fam) : family(fam) {}

// 辅助函数：创建缓存键
template<typename T>
typename FBIReducer<T>::CacheKey FBIReducer<T>::makeKey(int sectorNu, const std::vector<int>& nu, T delta) const {
    return std::make_tuple(sectorNu, nu, delta);
}

// 辅助函数：比较两个 nu 向量是否相等
template<typename T>
bool FBIReducer<T>::nuEqual(const std::vector<int>& nu1, const std::vector<int>& nu2) const {
    if (nu1.size() != nu2.size()) return false;
    for (size_t i = 0; i < nu1.size(); i++) {
        if (nu1[i] != nu2[i]) return false;
    }
    return true;
}

// 辅助函数：找到 nu 中最大元素的索引
// 计算两种索引
// 1. 最大元素在topsector中的索引, 也即长度为numProps
// 2. 最大元素在当前sector中的索引, 也即长度为|nuSector|  
template<typename T>
std::pair<int, int> FBIReducer<T>::findMaxIndex(const std::vector<int>& nu) const {
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
}

// 辅助函数：从 nu 向量获取 sectorNu
template<typename T>
int FBIReducer<T>::getSectorNu(const std::vector<int>& nu) const {
    int sectorNu = 0;
    for (size_t i = 0; i < nu.size(); i++) {
        sectorNu = sectorNu * 2 + (nu[i] > 0 ? 1 : 0);
    }
    return sectorNu;
}

// 辅助函数：检查 nu 是否是 Master FBI
template<typename T>
bool FBIReducer<T>::isMasterFBI(const std::vector<int>& nu, int sectorNu, const Family<T>* fam) const {
    // 将 sectorNu 转换为 nu 向量
    int numProps = fam->getNumProps();
    std::vector<int> sectorNuVec(numProps);
    int temp = sectorNu;
    for (int i = numProps - 1; i >= 0; i--) {
        sectorNuVec[i] = temp % 2;
        temp /= 2;
    }
    
    // 比较 nu 和 sectorNuVec
    return nuEqual(nu, sectorNuVec);
}

// 通用约化函数（根据 case 调用相应的约化方法）
template<typename T>
std::vector<T> FBIReducer<T>::reduceInternal(int sectorNu, const std::vector<int>& nu, T delta) {
    // 获取 sector
    const Sector<T>* sector = family->getSectorByNu(sectorNu);
    if (!sector) {
        throw std::runtime_error("Sector not found for sectorNu");
    }
    
    // 根据 sector 的 case 调用相应的约化方法
    int caseNum = sector->getCase();
    
    if (caseNum == 0) {
        return reduceCase0(sectorNu, nu, delta);
    } else {
        // 其他 case 的实现留待后续添加
        throw std::runtime_error("Reduction for case " + std::to_string(caseNum) + " not implemented yet");
    }
}

// Case 0 约化函数
template<typename T>
std::vector<T> FBIReducer<T>::reduceCase0(int sectorNu, const std::vector<int>& nuPlus, T delta) {
    // 获取所有 case 0 的 sectorNu（MFBIs）
    const auto& cases = family->getCases();
    const auto& case0Sectors = cases[0];
    int numMFBIs = case0Sectors.size();
    
    // 步骤b：检查是否是 MFBI
    if (isMasterFBI(nuPlus, sectorNu, family)) {
        // 返回单位向量（该 FBI 本身就是 MFBI）
        // 找到这个 MFBI 在所有 MFBIs 中的索引
        int mfbiIndex = -1;
        for (size_t i = 0; i < case0Sectors.size(); i++) {
            if (case0Sectors[i] == sectorNu) {
                mfbiIndex = i;
                break;
            }
        }
        
        if (mfbiIndex == -1) {
            throw std::runtime_error("MFBI not found in case 0 sectors");
        }
        
        // 创建单位向量
        std::vector<T> result(numMFBIs, T(0));
        result[mfbiIndex] = T(1);
        return result;
    }
    
    // 步骤c & d：先用 IBP 得到 delta+1，再降维到 delta
    auto coeffDeltaPlus1 = case0IBP(sectorNu, nuPlus, delta);
    return case0DimensionShift(sectorNu, nuPlus, coeffDeltaPlus1, delta);
}

// Case 0：步骤1 - IBP 恒等式（从 delta 到 delta+1）
template<typename T>
std::vector<T> FBIReducer<T>::case0IBP(int sectorNu, const std::vector<int>& nuPlus, T delta) {
    // check coeffs of FBI(nuPlus, delta+1) in cache, if exists, return directly, else compute and store in cache
    auto key = makeKey(sectorNu, nuPlus, delta + 1);
    auto it = cache.find(key);
    if (it != cache.end()) {
        return it->second;
    }

    // 找到 nuPlus 中最大的元素索引
    auto [maxIdxInTopSector, maxIdxInCurrentSector] = findMaxIndex(nuPlus);

    // 构造 nu = nuPlus - e_maxIdx
    std::vector<int> nu = nuPlus;
    nu[maxIdxInTopSector]--;
    
    // 递归获取 FBI(nu, delta) 的约化系数
    const auto& coeffNu = getReductionCoeff(nu, delta);
    
    // 获取 sector 信息
    const Sector<T>* sector = family->getSectorByNu(sectorNu);
    if (!sector) {
        throw std::runtime_error("Sector not found");
    }
    
    // 计算 |nu| (nu 的元素和)
    int nuSum = 0;
    for (int elem : nu) {
        nuSum += elem > 0 ? 1 : 0;
    }
    
    // 使用 invS 计算 FBI(nuPlus, delta+1) 的约化系数
    // FBI(nuPlus, delta+1) = FBI(nu + e_maxIdx, delta+1)
    // = (invS * {-FBI(nu,delta)_1, ..., FBI(nu,delta)_B, 
    //            FBI(nu-e_0,delta), ..., FBI(nu-e_{|nu|-1},delta)})_{B+maxIdx}
    
    int numBranch = family->getNumBranch();
    int numProps = family->getNumProps();
    int n = numBranch + nuSum;  // 矩阵规模是 B + |nu|
    int numMaster = family->getCases()[0].size();
    
    // 构造右侧向量
    // std::vector<T> rhs(n, T(0));
    std::vector<std::vector<T>> rhs(n, std::vector<T>(numMaster));
    
    // 前 B 个元素：-FBI(nu,delta) 在各个 branch 的系数（带负号）
    for (int i = 0; i < numBranch && i < coeffNu.size(); i++) {
        rhs[i] = coeffNu;
        for (auto& x : rhs[i]) x = -x;
    }
    
    // 获取 topS 的 branch 信息
    Sector<T> topSector(family->getTopS(), numProps, numBranch);
    auto topBranch = topSector.getBranch();
    
    // 后 N 个元素：FBI(nu - e_i, delta) 的约化系数
    int iCurrentSector = -1;
    for (int i = 0; i < numProps; i++) {
        std::vector<int> nuMinusEi = nu;
        if (nuMinusEi[i] > 0) {
            nuMinusEi[i]--;
            iCurrentSector++;
            
            // 检查 nuMinusEi 是否会导致某个 branch 被完全删除
            std::vector<int> deletedProps;
            for (int j = 0; j < numProps; j++) {
                if (nuMinusEi[j] == 0) {
                    deletedProps.push_back(j);
                }
            }
            
            // 如果会删除 branch，跳过
            if (Sector<T>::isDeleteBranch(topBranch, deletedProps)) {
                continue;
            }
            
            const auto& coeffNuMinusEi = getReductionCoeff(nuMinusEi, delta);
            // 这里需要适当处理系数
            // 暂时简化处理
            if (!coeffNuMinusEi.empty()) {
                rhs[numBranch + iCurrentSector] = coeffNuMinusEi;  // 简化处理
            }
        }
    }
    
    // 使用 rowOperation (invS) 求解
    const auto& invS = sector->getRowOperation();
    std::vector<T> solution(numMaster, T(0));
    
    for (int i = 0; i < numMaster; i++) {
        for (int j = 0; j < n; j++) {
            solution[i] += invS[numBranch+maxIdxInCurrentSector][j] * rhs[j][i];
        }
        solution[i] /= nu[maxIdxInTopSector];
    }
    
    // solution 只是IBP方程的临时解，需要递归获取 FBI(nuPlus, delta+1) 的真正约化系数
    // 但是 delta+1 的递归会导致问题，所以这里应该直接返回 coeffNu
    // TODO: 这里的逻辑需要重新设计


    // put the computed solution into cache
    cache[key] = solution;  

    return solution;
}

// Case 0：步骤2 - 维度变换（从 delta+1 到 delta）
template<typename T>
std::vector<T> FBIReducer<T>::case0DimensionShift(int sectorNu, const std::vector<int>& nuPlus, 
                                                    const std::vector<T>& coeffDeltaPlus1, T delta) {
    // FBI(nuPlus, delta) = 1/C * (z0*(2*delta-|nuPlus|-B) * FBI(nuPlus, delta+1) 
    //                             + sum_i z_i * FBI(nuPlus - e_i, delta))
    
    std::cerr << "DEBUG: case0DimensionShift called with sectorNu=" << sectorNu << std::endl;
    
    const Sector<T>* sector = family->getSectorByNu(sectorNu);
    if (!sector) {
        throw std::runtime_error("Sector not found");
    }
    
    std::cerr << "DEBUG: Getting sector info..." << std::endl;
    
    T C = sector->getCSum();
    T z0 = sector->getZ0();
    int numBranch = family->getNumBranch();
    int numProps = family->getNumProps();
    
    std::cerr << "DEBUG: C=" << C << ", z0=" << z0 << ", numBranch=" << numBranch << ", numProps=" << numProps << std::endl;
    
    // 计算 |nuPlus| (nuPlus 的和)
    int nuSum = 0;
    for (int n : nuPlus) {
        nuSum += n;
    }
    
    std::cerr << "DEBUG: nuSum=" << nuSum << std::endl;
    
    // 计算常数项：z0 * (2*delta - |nuPlus| - B)
    T constant = z0 * (T(2) * (delta + 1) - T(nuSum) - T(numBranch));
    
    std::cerr << "DEBUG: constant calculated" << std::endl;
    
    // 初始化结果向量：constant * coeffDeltaPlus1
    std::vector<T> result(coeffDeltaPlus1.size());
    for (size_t i = 0; i < coeffDeltaPlus1.size(); i++) {
        result[i] = constant * coeffDeltaPlus1[i];
    }
    
    std::cerr << "DEBUG: Starting loop over numProps" << std::endl;
    
    // 获取 topS 的 branch 信息
    Sector<T> topSector(family->getTopS(), numProps, numBranch);
    auto topBranch = topSector.getBranch();
    
    // 添加 sum_i z_i * FBI(nuPlus - e_i, delta) 的贡献
    int iCurrentSector = 0;
    for (int i = 0; i < numProps; i++) {
        std::cerr << "DEBUG: Processing prop " << i << std::endl;
        T z_i = sector->getZ(iCurrentSector);
        iCurrentSector += (nuPlus[i] > 0) ? 1 : 0;
        
        std::vector<int> nuPlusMinusEi = nuPlus;
        if (nuPlusMinusEi[i] > 0) {
            nuPlusMinusEi[i]--;
            
            // 检查 nuPlusMinusEi 是否会导致某个 branch 被完全删除
            std::vector<int> deletedProps;
            for (int j = 0; j < numProps; j++) {
                if (nuPlusMinusEi[j] == 0) {
                    deletedProps.push_back(j);
                }
            }
            
            // 如果会删除 branch，跳过
            if (Sector<T>::isDeleteBranch(topBranch, deletedProps)) {
                std::cerr << "DEBUG: Skipping prop " << i << " (would delete branch)" << std::endl;
                continue;
            }
            
            std::cerr << "DEBUG: Calling getReductionCoeff recursively for prop " << i << std::endl;
            const auto& coeffNuPlusMinusEi = getReductionCoeff(nuPlusMinusEi, delta);
            std::cerr << "DEBUG: Got coefficients for prop " << i << std::endl;
            
            // 将 z_i * coeffNuPlusMinusEi 加到结果中
            for (size_t j = 0; j < coeffNuPlusMinusEi.size() && j < result.size(); j++) {
                result[j] += z_i * coeffNuPlusMinusEi[j];
            }
        }
    }
    
    // 除以 C
    for (size_t i = 0; i < result.size(); i++) {
        result[i] = result[i] / C;
    }
    
    return result;
}

// 主接口：获取约化系数（自动使用缓存）
template<typename T>
const std::vector<T>& FBIReducer<T>::getReductionCoeff(const std::vector<int>& nu, T delta) {
    // 从 nu 获取 sectorNu
    int sectorNu = getSectorNu(nu);
    
    // 检查缓存
    auto key = makeKey(sectorNu, nu, delta);
    auto it = cache.find(key);
    
    if (it != cache.end()) {
        // 缓存命中
        return it->second;
    }
    
    // 缓存未命中，计算约化系数
    auto result = reduceInternal(sectorNu, nu, delta);
    
    // 存入缓存
    cache[key] = result;
    
    return cache[key];
}

// 主接口：为 FBI 对象获取约化系数
template<typename T>
const std::vector<T>& FBIReducer<T>::getReductionCoeff(const FBI<T>& fbi) {
    return getReductionCoeff(fbi.getNu(), fbi.getDelta());
}

// 构造 MFBIs（特定维度的所有 master FBI）
template<typename T>
std::vector<FBI<T>> FBIReducer<T>::buildMFBIs(T delta) const {
    std::vector<FBI<T>> mfbis;
    
    // 获取所有 case 0 的 sectorNu
    const auto& cases = family->getCases();
    const auto& case0Sectors = cases[0];
    
    // 为每个 case 0 sector 创建对应的 master FBI
    for (int sectorNu : case0Sectors) {
        // 将 sectorNu 转换为 nu 向量
        int numProps = family->getNumProps();
        std::vector<int> nu(numProps);
        int temp = sectorNu;
        for (int i = numProps - 1; i >= 0; i--) {
            nu[i] = temp % 2;
            temp /= 2;
        }
        
        // 创建 FBI
        FBI<T> mfbi(family, sectorNu, nu, delta);
        mfbis.push_back(mfbi);
    }
    
    return mfbis;
}

// 清空缓存
template<typename T>
void FBIReducer<T>::clearCache() {
    cache.clear();
}

// 获取缓存大小
template<typename T>
size_t FBIReducer<T>::getCacheSize() const {
    return cache.size();
}
