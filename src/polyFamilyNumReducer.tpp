#include <iostream>

// ===== PolyFamilyNumReducer 模板类实现 =====

template<typename T>
PolyFamilyNumReducer<T>::PolyFamilyNumReducer(const std::vector<std::vector<Polynomial<T>>>& polyTopS, 
                                               int numProps, 
                                               int numBranch)
    : polyTopS_(polyTopS), numProps_(numProps), numBranch_(numBranch) {
    // 检查 polyTopS 的维度
    if (polyTopS_.size() != static_cast<size_t>(numProps + numBranch)) {
        throw std::invalid_argument("polyTopS row size must equal numProps + numBranch");
    }
    for (const auto& row : polyTopS_) {
        if (row.size() != static_cast<size_t>(numProps + numBranch)) {
            throw std::invalid_argument("polyTopS must be a square matrix");
        }
    }
}

template<typename T>
PolyFamilyNumReducer<T>::~PolyFamilyNumReducer() {
    // unique_ptr 会自动清理，但我们可以显式清空
    clearCache();
}

template<typename T>
std::vector<std::vector<T>> PolyFamilyNumReducer<T>::evaluateTopS(const T& X, const T& Y) const {
    int size = numProps_ + numBranch_;
    std::vector<std::vector<T>> numericTopS(size, std::vector<T>(size));
    
    // 对每个多项式在 (X, Y) 处求值
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            numericTopS[i][j] = polyTopS_[i][j].evaluate(X, Y);
        }
    }
    
    return numericTopS;
}

template<typename T>
std::pair<Family<T>*, FBIReducer<T>*> PolyFamilyNumReducer<T>::getOrCreateReducer(const T& X, const T& Y) {
    auto key = std::make_pair(X, Y);
    
    // 检查缓存中是否已存在
    auto it = familyCache_.find(key);
    if (it != familyCache_.end()) {
        // 缓存命中
        return std::make_pair(it->second.first.get(), it->second.second.get());
    }
    
    // 缓存未命中，创建新的 Family 和 Reducer
    auto numericTopS = evaluateTopS(X, Y);
    
    // 创建 Family（使用 unique_ptr 管理）
    auto family = std::make_unique<Family<T>>(numericTopS, numProps_, numBranch_);
    
    // 创建 FBIReducer（使用 unique_ptr 管理）
    auto reducer = std::make_unique<FBIReducer<T>>(family.get());
    
    // 获取原始指针用于返回
    Family<T>* family_ptr = family.get();
    FBIReducer<T>* reducer_ptr = reducer.get();
    
    // 存入缓存
    familyCache_[key] = std::make_pair(std::move(family), std::move(reducer));
    
    return std::make_pair(family_ptr, reducer_ptr);
}

template<typename T>
std::vector<T> PolyFamilyNumReducer<T>::getReductionCoeff(const std::vector<int>& nu, T delta, const T& X, const T& Y) {
    // 获取或创建对应 (X, Y) 的 Reducer
    auto [family_ptr, reducer_ptr] = getOrCreateReducer(X, Y);
    
    // 调用 Reducer 的 getReductionCoeff 方法
    return reducer_ptr->getReductionCoeff(nu, delta);
}

template<typename T>
size_t PolyFamilyNumReducer<T>::getNumMasterFBIs() const {
    // 创建一个临时的 Family 来获取 master FBI 数量
    // 使用任意的 (X, Y) 值，例如 (0, 0)
    T zero(0);
    auto numericTopS = evaluateTopS(zero, zero);
    Family<T> temp_family(numericTopS, numProps_, numBranch_);
    return temp_family.getNumMaster();
}

template<typename T>
void PolyFamilyNumReducer<T>::clearCache() {
    familyCache_.clear();
}

template<typename T>
void PolyFamilyNumReducer<T>::printCacheInfo() const {
    std::cout << "PolyFamily Cache Information:\n";
    std::cout << "  Total cached (X, Y) points: " << familyCache_.size() << "\n";
    
    if (!familyCache_.empty()) {
        std::cout << "  Cached points:\n";
        int idx = 0;
        for (const auto& [key, value] : familyCache_) {
            std::cout << "    [" << idx << "] (X=" << key.first 
                      << ", Y=" << key.second << ")\n";
            idx++;
        }
    }
}
