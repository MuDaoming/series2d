#include <iostream>
#include <thread>
#include "../include/fbi_reducer.hpp"

// ===== FBIBlackBox 实现 =====

namespace firefly {
    template<typename T>
    template<typename FFIntTemp>
    std::vector<FFIntTemp> FBIBlackBox<T>::operator()(const std::vector<FFIntTemp>& values) {
        // 从 values 中提取 X 和 Y
        T X, Y;
        if constexpr (std::is_same_v<FFIntTemp, FFInt>) {
            X = T(values[0].n);
            Y = T(values[1].n);
        } else {
            X = T(values[0][0].n);
            Y = T(values[1][0].n);
        }
        
        // 每次创建新的FBIReducer（线程安全）
        FBIReducer<T> reducer = interpolater_->createReducer(X, Y);
        
        // 存储所有FBI的约化系数
        std::vector<FFIntTemp> result;
        
        // 遍历所有FBI
        for (const auto& [nu, delta] : fbi_list_) {
            // 获取当前FBI的约化系数
            std::vector<T> coeffs = reducer.getReductionCoeff(nu, delta);
            
            // 添加到结果向量
            for (const auto& c : coeffs) {
                result.emplace_back(c.n);
            }
        }
        
        return result;
    }
}

// ===== FBIInterpolater 模板类实现 =====

template<typename T>
FBIInterpolater<T>::FBIInterpolater(const std::vector<std::vector<Polynomial<T>>>& topS,
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
    
    // 计算master FBI数量（需要先创建一个临时的FBIReducer）
    // 在点(1,1)处求值topS
    auto numTopS = evaluateTopS(T(1), T(1));
    // debug
    // print numTopS
    for (const auto& row : numTopS) {
        for (const auto& val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    FBIReducer<T> tempReducer(numTopS, numProps, numBranch, delta);
    numMasterFBIs_ = tempReducer.getNumMaster();
}

template<typename T>
std::vector<std::vector<T>> FBIInterpolater<T>::evaluateTopS(const T& X, const T& Y) const {
    std::vector<std::vector<T>> result(topS_.size(), std::vector<T>(topS_[0].size()));
    for (size_t i = 0; i < topS_.size(); ++i) {
        for (size_t j = 0; j < topS_[i].size(); ++j) {
            result[i][j] = topS_[i][j].evaluate(X, Y);
        }
    }
    return result;
}

template<typename T>
FBIReducer<T> FBIInterpolater<T>::createReducer(const T& X, const T& Y) const {
    auto numTopS = evaluateTopS(X, Y);
    return FBIReducer<T>(numTopS, numProps_, numBranch_, delta_);
}

template<typename T>
Polynomial<T> FBIInterpolater<T>::convertToPolynomial(const firefly::PolynomialFF& ff_poly) {
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
Rational<T> FBIInterpolater<T>::convertToRational(const firefly::RationalFunctionFF& ff_result) {
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

/// 批量接口实现：一次性插值所有FBI的约化系数
template<typename T>
std::vector<std::vector<Rational<T>>>
FBIInterpolater<T>::getReductionCoeff(
    const std::vector<std::pair<std::vector<int>, T>>& fbi_list) {
    
    std::cout << "\n=== Batch Interpolation for " << fbi_list.size() << " FBIs ===\n";
    std::cout << "  Each FBI has " << numMasterFBIs_ << " coefficients\n";
    std::cout << "  Total functions to interpolate: " << fbi_list.size() * numMasterFBIs_ << "\n\n";
    
    // 创建 FBIBlackBox
    firefly::FBIBlackBox<T> bb(this, fbi_list);
    
    // 创建 Reconstructor (2个变量: X 和 Y)
    std::cout << "  Creating Reconstructor...\n";
    firefly::Reconstructor<firefly::FBIBlackBox<T>> reconst(
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
    std::cout << "  Starting batch reconstruction...\n";
    reconst.reconstruct(1);
    firefly::FFInt::set_new_prime(prime_);

    std::cout << "\n  Batch reconstruction completed!\n";
    
    // 获取结果（有限域上的有理函数）
    std::vector<firefly::RationalFunctionFF> ff_results = reconst.get_result_ff();
    
    size_t expected_size = fbi_list.size() * numMasterFBIs_;
    std::cout << "  Reconstructed functions: " << ff_results.size() 
              << " (expected: " << expected_size << ")\n";
    
    if (ff_results.size() != expected_size) {
        std::cerr << "  Warning: Result size mismatch!\n";
    }
    
    // 转换为 Rational<T> 类型并组织成vector，保持fbi_list的顺序
    std::vector<std::vector<Rational<T>>> result;
    result.reserve(fbi_list.size());
    
    size_t idx = 0;
    for (size_t fbi_idx = 0; fbi_idx < fbi_list.size(); ++fbi_idx) {
        std::vector<Rational<T>> coeffs;
        coeffs.reserve(numMasterFBIs_);
        
        for (size_t i = 0; i < numMasterFBIs_; ++i) {
            if (idx < ff_results.size()) {
                coeffs.push_back(convertToRational(ff_results[idx++]));
            }
        }
        
        result.push_back(coeffs);
    }
    
    std::cout << "\n=== Batch Interpolation Complete ===\n\n";
    
    return result;
}

template<typename T>
void FBIInterpolater<T>::printStats() const {
    std::cout << "FBIInterpolater Statistics:\n";
    std::cout << "  Prime: " << prime_ << "\n";
    std::cout << "  Number of master FBIs: " << numMasterFBIs_ << "\n";
    std::cout << "  Total probes (last call): " << total_probes_ << "\n";
}
