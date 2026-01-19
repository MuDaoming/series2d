#include <iostream>
#include <thread>
#include "../include/fbi_reducer.hpp"
#include "../include/de_builder.hpp"
// ===== DEBlackBox 实现 =====

namespace firefly {
    template<typename T>
    template<typename FFIntTemp>
    std::vector<FFIntTemp> DEBlackBox<T>::operator()(const std::vector<FFIntTemp>& values) {
        // 从 values 中提取 X 和 Y
        T X, Y;
        if constexpr (std::is_same_v<FFIntTemp, FFInt>) {
            X = T(values[0].n);
            Y = T(values[1].n);
        } else {
            X = T(values[0][0].n);
            Y = T(values[1][0].n);
        }
        
        // 每次创建新的DEBuilder（线程安全）
        DEBuilder<T> deBuilder = interpolater_->createDEBuilder(X, Y);
        
        // 计算微分方程矩阵
        std::vector<std::vector<T>> AX, AY;
        deBuilder.buildDEMatrices(AX, AY);
        
        // 扁平化结果：先AX，后AY
        std::vector<FFIntTemp> result;
        int M = AX.size();
        result.reserve(2 * M * M);
        
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < M; ++j) {
                result.emplace_back(AX[i][j].n);
            }
        }
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < M; ++j) {
                result.emplace_back(AY[i][j].n);
            }
        }
        
        return result;
    }
}

// ===== DEInterpolater 模板类实现 =====

template<typename T>
DEInterpolater<T>::DEInterpolater(const std::vector<std::vector<Polynomial<T>>>& topS,
                                  int numProps,
                                  int numBranch,
                                  T delta,
                                  uint64_t prime)
    : topS_(topS),
      numProps_(numProps),
      numBranch_(numBranch),
      delta_(delta),
      prime_(prime),
      total_probes_(0) {
    
    // 计算符号导数矩阵
    computeSymbolicDerivatives();
    
    // 计算master FBI数量（和FBIInterpolater一样，只创建FBIReducer）
    auto numTopS = evaluateTopS(T(1), T(1));
    // debug - print numTopS
    for (const auto& row : numTopS) {
        for (const auto& val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    FBIReducer<T> tempReducer(numTopS, numProps, numBranch, delta);
    numMasterFBI_ = tempReducer.getNumMaster();
}

template<typename T>
void DEInterpolater<T>::computeSymbolicDerivatives() {
    // 初始化导数矩阵
    dRdX_.resize(numProps_, std::vector<Polynomial<T>>(numProps_));
    dRdY_.resize(numProps_, std::vector<Polynomial<T>>(numProps_));
    
    // topS矩阵的结构：前B行/列是branch，后N行/列是传播子
    // R矩阵是topS的右下角N×N子矩阵
    for (int i = 0; i < numProps_; ++i) {
        for (int j = 0; j < numProps_; ++j) {
            dRdX_[i][j] = topS_[i + numBranch_][j + numBranch_].derivativeX();
            dRdY_[i][j] = topS_[i + numBranch_][j + numBranch_].derivativeY();
        }
    }
}

template<typename T>
std::vector<std::vector<T>> DEInterpolater<T>::evaluateTopS(const T& X, const T& Y) const {
    std::vector<std::vector<T>> result(topS_.size(), std::vector<T>(topS_[0].size()));
    for (size_t i = 0; i < topS_.size(); ++i) {
        for (size_t j = 0; j < topS_[i].size(); ++j) {
            result[i][j] = topS_[i][j].evaluate(X, Y);
        }
    }
    return result;
}

template<typename T>
std::vector<std::vector<T>> DEInterpolater<T>::evaluateDRdX(const T& X, const T& Y) const {
    std::vector<std::vector<T>> result(numProps_, std::vector<T>(numProps_));
    for (int i = 0; i < numProps_; ++i) {
        for (int j = 0; j < numProps_; ++j) {
            result[i][j] = dRdX_[i][j].evaluate(X, Y);
        }
    }
    return result;
}

template<typename T>
std::vector<std::vector<T>> DEInterpolater<T>::evaluateDRdY(const T& X, const T& Y) const {
    std::vector<std::vector<T>> result(numProps_, std::vector<T>(numProps_));
    for (int i = 0; i < numProps_; ++i) {
        for (int j = 0; j < numProps_; ++j) {
            result[i][j] = dRdY_[i][j].evaluate(X, Y);
        }
    }
    return result;
}

template<typename T>
DEBuilder<T> DEInterpolater<T>::createDEBuilder(const T& X, const T& Y) const {
    auto numTopS = evaluateTopS(X, Y);
    auto numDRdX = evaluateDRdX(X, Y);
    auto numDRdY = evaluateDRdY(X, Y);
    return DEBuilder<T>(numTopS, numDRdX, numDRdY, numProps_, numBranch_, delta_);
}

template<typename T>
Polynomial<T> DEInterpolater<T>::convertToPolynomial(const firefly::PolynomialFF& ff_poly) {
    Polynomial<T> result;
    
    // 遍历 PolynomialFF 的所有项
    for (const auto& term : ff_poly.coefs) {
        const auto& powers = term.first;   // vector<uint32_t>: [x_power, y_power]
        const auto& coeff_ff = term.second; // FFInt
        
        // 转换系数：FFInt → T
        T coeff(coeff_ff.n);
        
        // 获取幂次（假设是二元，X 和 Y）
        int x_power = (powers.size() > 0) ? static_cast<int>(powers[0]) : 0;
        int y_power = (powers.size() > 1) ? static_cast<int>(powers[1]) : 0;
        
        // 添加单项式
        result.addMonomial(coeff, Power(x_power, y_power));
    }
    
    return result;
}

template<typename T>
Rational<T> DEInterpolater<T>::convertToRational(const firefly::RationalFunctionFF& ff_result) {
    // 转换分子和分母
    Polynomial<T> numerator = convertToPolynomial(ff_result.numerator);
    Polynomial<T> denominator = convertToPolynomial(ff_result.denominator);
    
    // 检查分母是否为空（零多项式）
    if (denominator.isEmpty()) {
        // 如果分母为空，设置为常数多项式 1
        denominator = Polynomial<T>();
        denominator.addMonomial(T(1), Power(0, 0));
    }
    
    return Rational<T>(numerator, denominator);
}

template<typename T>
void DEInterpolater<T>::buildDEMatrices(std::vector<std::vector<Rational<T>>>& AX,
                                        std::vector<std::vector<Rational<T>>>& AY) {
    std::cout << "\n=== Building DE Matrices via Interpolation ===\n";
    std::cout << "  Master FBI count: " << numMasterFBI_ << "\n";
    std::cout << "  Matrix size: " << numMasterFBI_ << " x " << numMasterFBI_ << "\n";
    std::cout << "  Total functions to interpolate: " << 2 * numMasterFBI_ * numMasterFBI_ << "\n\n";
    
    // 创建 DEBlackBox
    firefly::DEBlackBox<T> bb(this);
    
    // 创建 Reconstructor (2个变量: X 和 Y)
    std::cout << "  Creating Reconstructor...\n";
    firefly::Reconstructor<firefly::DEBlackBox<T>> reconst(
        2,  // n_vars (X, Y)
        std::thread::hardware_concurrency(),  // n_threads
        bb  // black box
    );
    
    // 设置自定义质数
    firefly::FFInt::set_new_prime(prime_);
    bb.prime_changed_internal();
    reconst.prime_it = 0;
    
    std::cout << "  Using prime: " << firefly::FFInt::p << "\n";
    std::cout << "  Threads: " << std::thread::hardware_concurrency() << "\n\n";
    
    // 执行重构（只在一个质数域，不做CRT）
    std::cout << "  Starting reconstruction...\n";
    reconst.reconstruct(1);
    firefly::FFInt::set_new_prime(prime_);
    
    std::cout << "\n  Reconstruction completed!\n";
    
    // 获取结果（有限域上的有理函数）
    std::vector<firefly::RationalFunctionFF> ff_results = reconst.get_result_ff();
    
    size_t expected_size = 2 * numMasterFBI_ * numMasterFBI_;
    std::cout << "  Reconstructed functions: " << ff_results.size() 
              << " (expected: " << expected_size << ")\n";
    
    if (ff_results.size() != expected_size) {
        std::cerr << "  Warning: Result size mismatch!\n";
    }
    
    // 解析结果
    AX.assign(numMasterFBI_, std::vector<Rational<T>>(numMasterFBI_));
    AY.assign(numMasterFBI_, std::vector<Rational<T>>(numMasterFBI_));
    
    size_t idx = 0;
    
    // 解析AX矩阵
    for (int i = 0; i < numMasterFBI_; ++i) {
        for (int j = 0; j < numMasterFBI_; ++j) {
            if (idx < ff_results.size()) {
                AX[i][j] = convertToRational(ff_results[idx++]);
            }
        }
    }
    
    // 解析AY矩阵
    for (int i = 0; i < numMasterFBI_; ++i) {
        for (int j = 0; j < numMasterFBI_; ++j) {
            if (idx < ff_results.size()) {
                AY[i][j] = convertToRational(ff_results[idx++]);
            }
        }
    }
    
    std::cout << "\n=== DE Matrices Building Complete ===\n\n";
}
