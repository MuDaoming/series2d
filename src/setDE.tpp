#pragma once

#include <stdexcept>
#include <iostream>

// DENumBuilder 实现

template<typename T>
DENumBuilder<T>::DENumBuilder(const std::vector<std::vector<T>>& topS,
                               const std::vector<std::vector<T>>& dRdX,
                               const std::vector<std::vector<T>>& dRdY,
                               std::shared_ptr<PolyFamilyNumReducer<T>> numReducer,
                               const Family<T>* family,
                               T delta,
                               const T& X,
                               const T& Y)
    : family_(family), numReducer_(numReducer), delta_(delta),
      X_(X), Y_(Y),
      topS_(topS), dRdX_(dRdX), dRdY_(dRdY) {
    
    numMasterFBI_ = family_->getMasterIdxs().size();
    numProps_ = family_->getNumProps();
    numBranch_ = family_->getNumBranch();
    
    // 设置当前工作点
    numReducer_->setCurrentPoint(X, Y);
}

template<typename T>
void DENumBuilder<T>::buildDEMatrices(std::vector<std::vector<T>>& AX,
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
            // 使用数值约化器约化FBI(nu, delta)
            // 注意：已在构造函数中通过 setCurrentPoint 设置了 (X, Y)
            auto reduction = numReducer_->getReductionCoeff(nu, delta);
            // reduction[j]是FBI(nu,delta)在第j个MFBI上的系数
            for (int j = 0; j < numMasterFBI_; ++j) {
                AX[i][j] += coeff * reduction[j];
            }
        }
        
        // 计算第i个MFBI对Y的导数
        auto derivY = computeMasterDerivativeY(i);
        for (const auto& [key, coeff] : derivY) {
            const auto& [nu, delta] = key;
            auto reduction = numReducer_->getReductionCoeff(nu, delta);
            for (int j = 0; j < numMasterFBI_; ++j) {
                AY[i][j] += coeff * reduction[j];
            }
        }
    }
}

template<typename T>
std::map<std::pair<std::vector<int>, T>, T> 
DENumBuilder<T>::computeMasterDerivativeX(int masterIdx) {
    // 获取第masterIdx个MFBI的nu
    const auto& masterIdxs = family_->getMasterIdxs();
    int sectorIdx = masterIdxs[masterIdx];
    auto nu = family_->secvecFromIdx(sectorIdx);
    
    // 使用dRdX_计算导数
    return computeFBIDerivative(nu, dRdX_);
}

template<typename T>
std::map<std::pair<std::vector<int>, T>, T> 
DENumBuilder<T>::computeMasterDerivativeY(int masterIdx) {
    const auto& masterIdxs = family_->getMasterIdxs();
    int sectorIdx = masterIdxs[masterIdx];
    auto nu = family_->secvecFromIdx(sectorIdx);
    
    return computeFBIDerivative(nu, dRdY_);
}

template<typename T>
std::map<std::pair<std::vector<int>, T>, T>
DENumBuilder<T>::computeFBIDerivative(const std::vector<int>& nu,
                                      const std::vector<std::vector<T>>& dR) {
    std::map<std::pair<std::vector<int>, T>, T> result;
    
    // 使用公式：d FBI(nu,delta)/dVar = sum_{i,j} -1/2 * dR[i][j] * factor * FBI(nu+e_i+e_j, delta+1)
    // 其中 factor = (i!=j ? nu[i]*nu[j] : nu[i]*(nu[i]+1))
    
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
            T newDelta = delta_ + T(1);
            
            // 系数
            T coeff = T(-1) / T(2) * dR[i][j] * factor;
            
            // 添加到结果
            auto key = std::make_pair(newNu, newDelta);
            result[key] += coeff;
        }
    }
    
    return result;
}

// DEBuilder 实现

template<typename T>
DEBuilder<T>::DEBuilder(const std::vector<std::vector<Polynomial<T>>>& polyTopS,
                        const Family<T>* family,
                        T delta,
                        uint64_t prime)
    : S_(polyTopS), family_(family), delta_(delta), prime_(prime) {
    
    numMasterFBI_ = family_->getMasterIdxs().size();
    numProps_ = family_->getNumProps();
    numBranch_ = family_->getNumBranch();
    
    // 计算dR/dX和dR/dY（符号形式）
    dRdX_.resize(numProps_, std::vector<Polynomial<T>>(numProps_));
    dRdY_.resize(numProps_, std::vector<Polynomial<T>>(numProps_));
    
    for (int i = 0; i < numProps_; ++i) {
        for (int j = 0; j < numProps_; ++j) {
            dRdX_[i][j] = S_[i + numBranch_][j + numBranch_].derivativeX();
            dRdY_[i][j] = S_[i + numBranch_][j + numBranch_].derivativeY();
        }
    }
    
    // 创建数值约化器
    numReducer_ = std::make_shared<PolyFamilyNumReducer<T>>(S_, numProps_, numBranch_);
}

template<typename T>
std::unique_ptr<DENumBuilder<T>> DEBuilder<T>::createNumBuilder(const T& X, const T& Y) {
    // 求值S、dRdX、dRdY
    auto evaluateMatrix = [&](const std::vector<std::vector<Polynomial<T>>>& mat) {
        std::vector<std::vector<T>> result(mat.size(), std::vector<T>(mat[0].size()));
        for (size_t i = 0; i < mat.size(); ++i) {
            for (size_t j = 0; j < mat[i].size(); ++j) {
                result[i][j] = mat[i][j].evaluate(X, Y);
            }
        }
        return result;
    };
    
    auto topS_val = evaluateMatrix(S_);
    auto dRdX_val = evaluateMatrix(dRdX_);
    auto dRdY_val = evaluateMatrix(dRdY_);
    
    // 重要：为每个点创建独立的 numReducer，避免多线程竞争
    // 在多线程环境下，共享 numReducer_ 会导致 setCurrentPoint 的竞争条件
    auto local_numReducer = std::make_shared<PolyFamilyNumReducer<T>>(S_, numProps_, numBranch_);
    
    return std::make_unique<DENumBuilder<T>>(
        topS_val, dRdX_val, dRdY_val,
        local_numReducer, family_, delta_, X, Y
    );
}

template<typename T>
void DEBuilder<T>::buildDEMatrices(std::vector<std::vector<Rational<T>>>& AX,
                                   std::vector<std::vector<Rational<T>>>& AY) {
    using namespace firefly;
    
    // 初始化矩阵
    AX.assign(numMasterFBI_, std::vector<Rational<T>>(numMasterFBI_));
    AY.assign(numMasterFBI_, std::vector<Rational<T>>(numMasterFBI_));
    
    // 设置质数
    FFInt::set_new_prime(prime_);
    
    // 创建BlackBox
    DEMatrixBlackBox<T> blackbox(this);
    
    // 创建Reconstructor - 使用6个线程
    Reconstructor<DEMatrixBlackBox<T>> reconstructor(
        2,  // n_vars (X, Y)
        6,  // n_threads (6线程)
        blackbox  // black box
    );
    
    // 确保使用我们设置的质数
    FFInt::set_new_prime(prime_);
    blackbox.prime_changed_internal();
    reconstructor.prime_it = 0;  // 重置质数迭代器
    
    // 执行重构
    std::cout << "开始Firefly重构微分方程矩阵..." << std::endl;
    reconstructor.reconstruct(1);  // 在一个质数域重构
    auto results = reconstructor.get_result_ff();  // 获取有限域结果
    std::cout << "Firefly重构完成！" << std::endl;
    
    // 解析结果
    int idx = 0;
    
    // 解析AX矩阵
    for (int i = 0; i < numMasterFBI_; ++i) {
        for (int j = 0; j < numMasterFBI_; ++j) {
            AX[i][j] = convertToRational(results[idx++]);
        }
    }
    
    // 解析AY矩阵
    for (int i = 0; i < numMasterFBI_; ++i) {
        for (int j = 0; j < numMasterFBI_; ++j) {
            AY[i][j] = convertToRational(results[idx++]);
        }
    }
}

template<typename T>
Rational<T> DEBuilder<T>::convertToRational(const firefly::RationalFunctionFF& ff_result) {
    // 转换分子
    Polynomial<T> numerator;
    for (const auto& term : ff_result.numerator.coefs) {
        const auto& powers = term.first;   // vector<uint32_t>
        const auto& coeff_ff = term.second; // FFInt
        
        if (powers.size() != 2) {
            throw std::runtime_error("Unexpected monomial dimension in Firefly result");
        }
        
        T coeff(coeff_ff.n);
        Power p{static_cast<int>(powers[0]), static_cast<int>(powers[1])};
        numerator.addMonomial(coeff, p);
    }
    
    // 转换分母
    Polynomial<T> denominator;
    for (const auto& term : ff_result.denominator.coefs) {
        const auto& powers = term.first;   // vector<uint32_t>
        const auto& coeff_ff = term.second; // FFInt
        
        if (powers.size() != 2) {
            throw std::runtime_error("Unexpected monomial dimension in Firefly result");
        }
        
        T coeff(coeff_ff.n);
        Power p{static_cast<int>(powers[0]), static_cast<int>(powers[1])};
        denominator.addMonomial(coeff, p);
    }
    
    // 如果分母为空，设为1
    if (denominator.isEmpty()) {
        denominator.addMonomial(T(1), Power{0, 0});
    }
    
    return Rational<T>{numerator, denominator};
}
