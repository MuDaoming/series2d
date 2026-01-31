// public getter functions
template<typename T>
int Family<T>::getCase(std::vector<int> nu) const
{
    // checke size
    if (nu.size() != static_cast<size_t>(numProps_)) {
        throw std::invalid_argument("nu size does not match numProps_");
    }
    std::vector<int> secvec(numProps_);
    for (size_t i = 0; i < nu.size(); i++) {
        secvec[i] = (nu[i] > 0) ? 1 : 0;
    }   
    int idx = idxFromSecvec(secvec);
    return cases_[idx];
}

template<typename T>
int Family<T>::getIndexOfMaster(std::vector<int> nu) const
{
    if (!isMaster(nu)) {
        return -1;
    }
    int idx = idxFromSecvec(nu);
    auto it = std::find(masterIdxs_.begin(), masterIdxs_.end(), idx);
    if (it != masterIdxs_.end()) {
        return std::distance(masterIdxs_.begin(), it);
    } else {
        throw std::runtime_error("Master nu index not found");
    }
}

template<typename T>
const Sector<T>* Family<T>::getSectorByIdx(int idx) const
{
    for (size_t i = 0; i < sectorIdxs_.size(); i++) {
        if (sectorIdxs_[i] == idx) {
            return &sectors_[i];
        }
    }
    return nullptr;  // 未找到
}

template<typename T>
const Sector<T>* Family<T>::getSectorBySecvec(std::vector<int> secvec) const
{
    return getSectorByIdx(idxFromSecvec(secvec));
}

// public auxiliary functions

template<typename T>
std::vector<int> Family<T>::secvecFromIdx(int n) const
{
    // check if n is in valid range
    if (n < 0 || n >= (1 << numProps_)) {
        throw std::out_of_range("nu index out of range");
    }
    std::vector<int> nu(numProps_);
    for (int i = numProps_ - 1; i >= 0; i--) {
        nu[i] = n % 2;
        n /= 2;
    }
    return nu;
}

template<typename T>
int Family<T>::idxFromSecvec(const std::vector<int>& secvec) const
{
    // check if length of secvec matches numProps_
    if (secvec.size() != static_cast<size_t>(numProps_)) {
        throw std::invalid_argument("secvec size does not match numProps_");
    }

    int idx = 0;
    for (int i = 0; i < numProps_; i++) {
        idx = (idx << 1) | secvec[i];
    }
    return idx;
}

template<typename T>
int Family<T>::nBranch(std::vector<int> const& nu) const
{
    // check if length of nu matches numProps_
    if (nu.size() != static_cast<size_t>(numProps_)) {
        throw std::invalid_argument("nu size does not match numProps_");
    }

    std::vector<bool> branchPresent(numBranch_, false);
    for (int i = 0; i < numProps_; ++i) {
        if (nu[i] > 0) {
            branchPresent[branchIndices_[i]] = true;
        }
    }

    return std::count(branchPresent.begin(), branchPresent.end(), true);
}

template<typename T>
int Family<T>::nProps(std::vector<int> const& nu) const
{
    // check if length of nu matches numProps_
    if (nu.size() != static_cast<size_t>(numProps_)) {
        throw std::invalid_argument("nu size does not match numProps_");
    }
    int nprops = 0;
    for (const auto& val : nu) {
        nprops += val > 0 ? 1 : 0;
    }
    return nprops;
}

template<typename T>
int Family<T>::nuSum(std::vector<int> const& nu) const
{
    // check if length of nu matches numProps_
    if (nu.size() != static_cast<size_t>(numProps_)) {
        throw std::invalid_argument("nu size does not match numProps_");
    }
    int sum = 0;
    for (const auto& val : nu) {
        sum += val;
    }
    return sum;
}

template<typename T>
bool Family<T>::isMaster(const std::vector<int>& nu) const
{
    for (int val : nu) {
        if (val != 0 && val != 1) {
            return false;
        }
    }

    int idx = idxFromSecvec(nu);

    return cases_[idx] == 0;
}

// private member functions

template<typename T>
void Family<T>::constructBranchIndices()
{
    branchIndices_.resize(numProps_);
    for (int i = 0; i < numProps_; i++) {
        for (int j = 0; j < numBranch_; j++) {
            if (topS_[numBranch_ + i][j] == T(1)) {
                branchIndices_[i] = j;
                continue;
            }
        }
    }
}

template<typename T>
bool Family<T>::getSubS(const std::vector<int>& nu, std::vector<std::vector<T>>& subS) const
{
    // check if length of nu matches numProps_
    if (nu.size() != static_cast<size_t>(numProps_)) {
        throw std::invalid_argument("nu size does not match numProps_");
    }
    // check nBranch(nu) == numBranch_
    if (nBranch(nu) != numBranch_) {
        return false;
    }

    int n = numBranch_ + numProps_;
    int sub_n = numBranch_ + nProps(nu);      // nBranch(nu) always equals numBranch_
    subS.resize(sub_n, std::vector<T>(sub_n,0));

    // 构造 subS
    int sub_i = 0;
    int sub_j = numBranch_;
    for (int i = 0; i < numBranch_; i++) {
        sub_i = i;
        sub_j = numBranch_;
        for (int j = numBranch_; j < numBranch_ + numProps_; j++) {
            if (nu[j - numBranch_] == 1) {
                subS[sub_i][sub_j] = topS_[i][j];
                subS[sub_j][sub_i] = topS_[i][j];
                sub_j++;
            }
        }
    }
    sub_i = numBranch_;
    for (int i = numBranch_; i < numBranch_ + numProps_; i++) {
        sub_j = numBranch_;
        if (nu[i - numBranch_] == 1) {
            for (int j = numBranch_; j < numBranch_ + numProps_; j++) {
                if (nu[j - numBranch_] == 1) {
                    if(i==4 && j==3)
                    {
                        int debug_var = 0;
                    }
                    subS[sub_i][sub_j] = topS_[i][j];
                    subS[sub_j][sub_i] = topS_[i][j];
                    sub_j++;
                }
            }
            sub_i++;
        }
    }
    return true;
}



// 构造所有有效的 sectors
template<typename T>
void Family<T>::findSectors() {
    int total = 1 << numProps_;  // 2^numProps

    cases_.resize(total);

    for (int i = total - 1; i >=0; i--) {
        std::vector<int> secvec = secvecFromIdx(i);
        std::vector<std::vector<T>> subS;
        bool isSector = getSubS(secvec, subS);
        if (isSector) {
            Sector<T> tempSector(std::move(subS), nProps(secvec), numBranch_);
            int caseNum = tempSector.getCase();
            cases_[i] = caseNum;
            if (caseNum == 0) masterIdxs_.push_back(i);
            sectors_.push_back(std::move(tempSector));
            sectorIdxs_.push_back(i);
        }
        else {
            cases_[i] = -1;
        }
    }
    
    // 初始化masterDeltas_，使用targetDelta_
    masterDeltas_.assign(masterIdxs_.size(), targetDelta_);
}
