#include <stdexcept>
#include <iostream>

// DEBuilder 实现

template<typename T>
DEBuilder<T>::DEBuilder(const std::vector<std::vector<T>>& topS,
                        const std::vector<std::vector<T>>& dRdX,
                        const std::vector<std::vector<T>>& dRdY,
                        int numProps,
                        int numBranch,
                        T delta)
    : reducer_(topS, numProps, numBranch, delta),
      dRdX_(dRdX),
      dRdY_(dRdY),
      numProps_(numProps),
      numBranch_(numBranch) {
    
    numMasterFBI_ = reducer_.getNumMaster();
}

template<typename T>
void DEBuilder<T>::buildDEMatrices(std::vector<std::vector<T>>& AX,
                                    std::vector<std::vector<T>>& AY) {
    // 初始化矩阵
    AX.assign(numMasterFBI_, std::vector<T>(numMasterFBI_, T(0)));
    AY.assign(numMasterFBI_, std::vector<T>(numMasterFBI_, T(0)));
    
    // 对每个MFBI计算导数
    for (int i = 0; i < numMasterFBI_; ++i) {
        // 计算第i个MFBI对X的导数
        auto derivX = computeMasterDerivativeX(i);
        
        // 约化导数中的每个FBI到MFBI
        for (const auto& [key, coeff] : derivX) {
            const auto& [nu, delta] = key;
            // 使用reducer约化FBI(nu, delta)
            auto reduction = reducer_.getReductionCoeff(nu, delta);
            // reduction[j]是FBI(nu,delta)在第j个MFBI上的系数
            for (int j = 0; j < numMasterFBI_; ++j) {
                AX[i][j] += coeff * reduction[j];
            }
        }
        
        // 计算第i个MFBI对Y的导数
        auto derivY = computeMasterDerivativeY(i);
        for (const auto& [key, coeff] : derivY) {
            const auto& [nu, delta] = key;
            auto reduction = reducer_.getReductionCoeff(nu, delta);
            for (int j = 0; j < numMasterFBI_; ++j) {
                AY[i][j] += coeff * reduction[j];
            }
        }
    }
}

template<typename T>
std::map<std::pair<std::vector<int>, T>, T> 
DEBuilder<T>::computeMasterDerivativeX(int masterIdx) {
    // 获取第masterIdx个MFBI的nu
    const auto& masterNus = reducer_.getMasterNus();
    auto nu = masterNus[masterIdx];
    
    // 使用dRdX_计算导数
    return computeFBIDerivative(nu, dRdX_);
}

template<typename T>
std::map<std::pair<std::vector<int>, T>, T> 
DEBuilder<T>::computeMasterDerivativeY(int masterIdx) {
    const auto& masterNus = reducer_.getMasterNus();
    auto nu = masterNus[masterIdx];
    
    return computeFBIDerivative(nu, dRdY_);
}

template<typename T>
std::map<std::pair<std::vector<int>, T>, T>
DEBuilder<T>::computeFBIDerivative(const std::vector<int>& nu,
                                   const std::vector<std::vector<T>>& dR) {
    std::map<std::pair<std::vector<int>, T>, T> result;
    
    // 使用公式：d FBI(nu,delta)/dVar = sum_{i,j} -1/2 * dR[i][j] * factor * FBI(nu+e_i+e_j, delta+1)
    // 其中 factor = (i!=j ? nu[i]*nu[j] : nu[i]*(nu[i]+1))
    
    // 获取工作delta
    const auto& masterDeltas = reducer_.getMasterDeltas();
    T delta = masterDeltas[0];  // 所有master FBI使用相同的delta
    
    for (int i = 0; i < numProps_; ++i) {
        for (int j = 0; j < numProps_; ++j) {
            if (dR[i][j] == T(0)) continue;
            
            // 计算factor
            T factor;
            if (i != j) {
                factor = T(nu[i]) * T(nu[j]);
            } else {
                factor = T(nu[i]) * T(nu[i] + 1);
            }
            
            if (factor == T(0)) continue;
            
            // 构造新的nu: nu + e_i + e_j
            std::vector<int> newNu = nu;
            newNu[i]++;
            newNu[j]++;
            
            // 新的delta
            T newDelta = delta + T(1);
            
            // 系数
            T coeff = T(-1) / T(2) * dR[i][j] * factor;
            
            // 添加到结果
            auto key = std::make_pair(newNu, newDelta);
            result[key] += coeff;
        }
    }
    
    return result;
}
