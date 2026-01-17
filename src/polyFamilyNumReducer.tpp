#include <iostream>

// ===== PolyFamilyNumReducer 模板类实现 =====

template<typename T>
PolyFamilyNumReducer<T>::PolyFamilyNumReducer(const std::vector<std::vector<Polynomial<T>>>& polyTopS, 
                                               int numProps, 
                                               int numBranch)
    : polyTopS_(polyTopS), numProps_(numProps), numBranch_(numBranch), 
      curX_(0), curY_(0), hasCurrentPoint_(false) {
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
    // unique_ptr 会自动清理
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
void PolyFamilyNumReducer<T>::setCurrentPoint(const T& X, const T& Y) {
    curX_ = X;
    curY_ = Y;
    hasCurrentPoint_ = true;
    
    // 在新点处求值 topS
    auto numericTopS = evaluateTopS(X, Y);
    
    // 创建新的 Family 和 Reducer（使用 unique_ptr 自动管理）
    curFamily_ = std::make_unique<Family<T>>(numericTopS, numProps_, numBranch_);
    curReducer_ = std::make_unique<FBIReducer<T>>(curFamily_.get());
}

template<typename T>
std::vector<T> PolyFamilyNumReducer<T>::getReductionCoeff(const std::vector<int>& nu, T delta) {
    if (!hasCurrentPoint_) {
        throw std::runtime_error("PolyFamilyNumReducer: Must call setCurrentPoint before getReductionCoeff!");
    }
    
    // 调用当前 Reducer 的 getReductionCoeff 方法
    return curReducer_->getReductionCoeff(nu, delta);
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
