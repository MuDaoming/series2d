#include <iostream>

// ===== PolyFamilyReducer 模板类实现 =====
// 使用 Firefly Reconstructor 进行批量插值重构

template<typename T>
PolyFamilyReducer<T>::PolyFamilyReducer(const std::vector<std::vector<Polynomial<T>>>& polyTopS,
                                         int numProps,
                                         int numBranch,
                                         uint64_t prime)
    : prime_(prime), total_probes_(0) {
    
    // 创建内部的数值约化器（使用 shared_ptr 便于传递给 BlackBox）
    numReducer_ = std::make_shared<PolyFamilyNumReducer<T>>(polyTopS, numProps, numBranch);
    
    // 获取 master FBI 数量
    numMasterFBIs_ = numReducer_->getNumMasterFBIs();
}

template<typename T>
Polynomial<T> PolyFamilyReducer<T>::convertToPolynomial(const firefly::PolynomialFF& ff_poly) {
    Polynomial<T> result;
    
    // 遍历 PolynomialFF 的所有项
    // ff_poly.coefs 是 map<vector<uint32_t>, FFInt>
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
Rational<T> PolyFamilyReducer<T>::convertToRational(const firefly::RationalFunctionFF& ff_result) {
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
std::vector<Rational<T>> PolyFamilyReducer<T>::getReductionCoeff(const std::vector<int>& nu, const T& delta) {
    // 设置质数
    firefly::FFInt::set_new_prime(prime_);
    
    // 创建 BlackBox
    firefly::PolyFamilyBlackBox<T> blackbox(numReducer_, nu, delta);
    
    // 创建 Reconstructor (2 个变量: X 和 Y)
    firefly::Reconstructor<firefly::PolyFamilyBlackBox<T>> reconst(
        2,  // n_vars (X, Y)
        1,  // n_threads (单线程)
        blackbox  // black box
    );
    
    // 确保使用我们设置的质数
    firefly::FFInt::set_new_prime(prime_);
    blackbox.prime_changed_internal();
    reconst.prime_it = 0;  // 重置质数迭代器
    
    // 执行重构（只在一个质数域，不做CRT）
    reconst.reconstruct(1);
    
    // 获取结果（有限域上的有理函数）
    std::vector<firefly::RationalFunctionFF> ff_results = reconst.get_result_ff();
    
    // 检查结果数量
    if (ff_results.size() != numMasterFBIs_) {
        std::cerr << "Warning: Expected " << numMasterFBIs_ << " results, got " << ff_results.size() << std::endl;
    }
    
    // 转换为 Rational<T> 类型
    std::vector<Rational<T>> result;
    result.reserve(ff_results.size());
    for (const auto& ff_result : ff_results) {
        result.push_back(convertToRational(ff_result));
    }
    
    // 更新统计信息（总探测点数）
    // 注意：Reconstructor 不直接提供总探测点数，这里设置为 0
    // 可以通过其他方式统计，例如记录 BlackBox 的调用次数
    total_probes_ = 0;  // 可以根据需要实现更精确的统计
    
    return result;
}

template<typename T>
std::vector<Rational<T>> PolyFamilyReducer<T>::getReductionCoeffDebug(const std::vector<int>& nu, const T& delta) {
    std::cout << "开始批量插值重构 (Debug模式)...\n";
    std::cout << "  Master FBI 数量: " << numMasterFBIs_ << "\n";
    
    // 记录插值前的 cache 大小
    size_t cache_before = numReducer_->getCacheSize();
    
    // 设置质数
    firefly::FFInt::set_new_prime(prime_);
    
    // 创建 BlackBox
    firefly::PolyFamilyBlackBox<T> blackbox(numReducer_, nu, delta);
    
    // 创建 Reconstructor (2 个变量: X 和 Y)
    firefly::Reconstructor<firefly::PolyFamilyBlackBox<T>> reconst(
        2,  // n_vars (X, Y)
        1,  // n_threads (单线程)
        blackbox  // black box
    );
    
    // 确保使用我们设置的质数
    firefly::FFInt::set_new_prime(prime_);
    blackbox.prime_changed_internal();
    reconst.prime_it = 0;  // 重置质数迭代器
    
    std::cout << "  执行 Firefly Reconstructor...\n";
    
    // 执行重构（只在一个质数域，不做CRT）
    reconst.reconstruct(1);
    
    // 记录插值后的 cache 大小
    size_t cache_after = numReducer_->getCacheSize();
    size_t total_points_used = cache_after - cache_before;
    
    std::cout << "  插值完成!\n";
    std::cout << "  使用的 (X,Y) 点数: " << total_points_used << "\n";
    std::cout << "  Cache 大小变化: " << cache_before << " -> " << cache_after << "\n";
    
    // 获取结果（有限域上的有理函数）
    std::vector<firefly::RationalFunctionFF> ff_results = reconst.get_result_ff();
    
    std::cout << "  重构的有理函数数量: " << ff_results.size() << "\n";
    
    // 检查结果数量
    if (ff_results.size() != numMasterFBIs_) {
        std::cerr << "  Warning: Expected " << numMasterFBIs_ << " results, got " << ff_results.size() << std::endl;
    }
    
    // 输出每个系数的简要信息
    for (size_t i = 0; i < ff_results.size(); ++i) {
        const auto& ff_result = ff_results[i];
        std::cout << "    Coefficient " << (i+1) << "/" << ff_results.size()
                  << ": " << ff_result.to_string({"X", "Y"}) << "\n";
    }
    
    // 转换为 Rational<T> 类型
    std::vector<Rational<T>> result;
    result.reserve(ff_results.size());
    for (const auto& ff_result : ff_results) {
        result.push_back(convertToRational(ff_result));
    }
    
    // 更新统计信息
    total_probes_ = 0;  // 可以根据需要实现更精确的统计
    
    return result;
}

/// 批量接口实现：一次性插值所有FBI的约化系数
template<typename T>
std::vector<Rational<T>> PolyFamilyReducer<T>::getBatchReductionCoeff(
    const std::vector<std::pair<std::vector<int>, T>>& fbi_list) {
    
    std::cout << "\n=== Batch Interpolation for " << fbi_list.size() << " FBIs ===\n";
    std::cout << "  Each FBI has " << numMasterFBIs_ << " coefficients\n";
    std::cout << "  Total functions to interpolate: " << fbi_list.size() * numMasterFBIs_ << "\n\n";
    
    // 记录插值前的 cache 大小
    size_t cache_before = numReducer_->getCacheSize();
    
    // 创建 BatchBlackBox
    firefly::PolyFamilyBatchBlackBox<T> bb(numReducer_, fbi_list);
    
    // 创建 Reconstructor (2个变量: X 和 Y)
    std::cout << "  Creating Reconstructor...\n";
    firefly::Reconstructor<firefly::PolyFamilyBatchBlackBox<T>> reconst(
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
    
    // 记录插值后的 cache 大小
    size_t cache_after = numReducer_->getCacheSize();
    size_t total_points_used = cache_after - cache_before;
    
    std::cout << "\n  Batch reconstruction completed!\n";
    std::cout << "  Used (X,Y) points: " << total_points_used << "\n";
    std::cout << "  Cache size: " << cache_before << " -> " << cache_after << "\n";
    
    // 获取结果（有限域上的有理函数）
    std::vector<firefly::RationalFunctionFF> ff_results = reconst.get_result_ff();
    
    size_t expected_size = fbi_list.size() * numMasterFBIs_;
    std::cout << "  Reconstructed functions: " << ff_results.size() 
              << " (expected: " << expected_size << ")\n";
    
    if (ff_results.size() != expected_size) {
        std::cerr << "  Warning: Result size mismatch!\n";
    }
    
    // 转换为 Rational<T> 类型
    std::vector<Rational<T>> result;
    result.reserve(ff_results.size());
    for (const auto& ff_result : ff_results) {
        result.push_back(convertToRational(ff_result));
    }
    
    std::cout << "\n=== Batch Interpolation Complete ===\n\n";
    
    return result;
}

template<typename T>
void PolyFamilyReducer<T>::printStats() const {
    std::cout << "PolyFamilyReducer Statistics:\n";
    std::cout << "  Prime: " << prime_ << "\n";
    std::cout << "  Number of master FBIs: " << numMasterFBIs_ << "\n";
    std::cout << "  Cache size (X,Y points): " << numReducer_->getCacheSize() << "\n";
    std::cout << "  Total probes (last call): " << total_probes_ << "\n";
}
