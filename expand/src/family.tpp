/**
 * Family class 实现
 * 管理所有sector，从topS构造各个子sector
 */

// 计算dR/dX和dR/dY矩阵
// R是topS的右下角N×N子矩阵：R[i][j] = topS_[numBranch_+i][numBranch_+j]
template<typename RT, typename PT, typename ST>
template<typename Symbol>
void Family<RT, PT, ST>::computeDerivatives(const Symbol& X, const Symbol& Y) {
    // 初始化dRdX和dRdY为N×N矩阵
    dRdX_.resize(numProps_, std::vector<PT>(numProps_));
    dRdY_.resize(numProps_, std::vector<PT>(numProps_));
    
    // 计算每个R[i][j]对X和Y的偏导数（R是多项式矩阵，导数也是多项式）
    for (int i = 0; i < numProps_; ++i) {
        for (int j = 0; j < numProps_; ++j) {
            // R[i][j] = topS_[numBranch_+i][numBranch_+j]
            const PT& Rij = topS_[numBranch_ + i][numBranch_ + j];
            dRdX_[i][j] = diff(Rij, X);
            dRdY_[i][j] = diff(Rij, Y);
        }
    }
}

// public getter functions
template<typename RT, typename PT, typename ST>
int Family<RT, PT, ST>::getCase(std::vector<int> nu) const
{
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

template<typename RT, typename PT, typename ST>
int Family<RT, PT, ST>::getIndexOfMaster(std::vector<int> nu) const
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

template<typename RT, typename PT, typename ST>
const Sector<RT, PT>* Family<RT, PT, ST>::getSectorByIdx(int idx) const
{
    for (size_t i = 0; i < sectorIdxs_.size(); i++) {
        if (sectorIdxs_[i] == idx) {
            return &sectors_[i];
        }
    }
    return nullptr;
}

template<typename RT, typename PT, typename ST>
const Sector<RT, PT>* Family<RT, PT, ST>::getSectorBySecvec(std::vector<int> secvec) const
{
    return getSectorByIdx(idxFromSecvec(secvec));
}

// public auxiliary functions

template<typename RT, typename PT, typename ST>
std::vector<int> Family<RT, PT, ST>::secvecFromIdx(int n) const
{
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

template<typename RT, typename PT, typename ST>
int Family<RT, PT, ST>::idxFromSecvec(const std::vector<int>& secvec) const
{
    if (secvec.size() != static_cast<size_t>(numProps_)) {
        throw std::invalid_argument("secvec size does not match numProps_");
    }

    int idx = 0;
    for (int i = 0; i < numProps_; i++) {
        idx = (idx << 1) | secvec[i];
    }
    return idx;
}

template<typename RT, typename PT, typename ST>
int Family<RT, PT, ST>::nBranch(std::vector<int> const& nu) const
{
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

template<typename RT, typename PT, typename ST>
int Family<RT, PT, ST>::nProps(std::vector<int> const& nu) const
{
    if (nu.size() != static_cast<size_t>(numProps_)) {
        throw std::invalid_argument("nu size does not match numProps_");
    }
    int nprops = 0;
    for (const auto& val : nu) {
        nprops += val > 0 ? 1 : 0;
    }
    return nprops;
}

template<typename RT, typename PT, typename ST>
int Family<RT, PT, ST>::nuSum(std::vector<int> const& nu) const
{
    if (nu.size() != static_cast<size_t>(numProps_)) {
        throw std::invalid_argument("nu size does not match numProps_");
    }
    int sum = 0;
    for (const auto& val : nu) {
        sum += val;
    }
    return sum;
}

template<typename RT, typename PT, typename ST>
bool Family<RT, PT, ST>::isMaster(const std::vector<int>& nu) const
{
    for (int val : nu) {
        if (val != 0 && val != 1) {
            return false;
        }
    }

    int idx = idxFromSecvec(nu);

    return std::find(masterIdxs_.begin(), masterIdxs_.end(), idx) !=
           masterIdxs_.end();
}

template<typename RT, typename PT, typename ST>
void Family<RT, PT, ST>::setSectorMap(const SectorMap& sectorMap)
{
    if (sectorMap.numProps() != 0 && sectorMap.numProps() != numProps_) {
        throw std::invalid_argument("SectorMap N does not match Family numProps");
    }
    for (const auto& entry : sectorMap.entries()) {
        if (cases_[idxFromSecvec(entry.source)] < 0 ||
            cases_[idxFromSecvec(entry.target)] < 0) {
            throw std::invalid_argument(
                "SectorMap source and target must be valid Family sectors");
        }
    }
    sectorMap_ = sectorMap;
    rebuildMasters();
    masterDeltas_.assign(masterIdxs_.size(), ST(0));
}

// private member functions

template<typename RT, typename PT, typename ST>
void Family<RT, PT, ST>::constructBranchIndices()
{
    branchIndices_.resize(numProps_);
    for (int i = 0; i < numProps_; i++) {
        for (int j = 0; j < numBranch_; j++) {
            // 检查 topS_[numBranch_ + i][j] == PT(1)
            PT diff = topS_[numBranch_ + i][j] - PT(1);
            if (isZero(diff)) {
                branchIndices_[i] = j;
                break;  // 找到了，跳出内层循环
            }
        }
    }
}

template<typename RT, typename PT, typename ST>
bool Family<RT, PT, ST>::getSubS(const std::vector<int>& nu, std::vector<std::vector<PT>>& subS) const
{
    if (nu.size() != static_cast<size_t>(numProps_)) {
        throw std::invalid_argument("nu size does not match numProps_");
    }
    if (nBranch(nu) != numBranch_) {
        return false;
    }

    int n = numBranch_ + numProps_;
    int sub_n = numBranch_ + nProps(nu);
    subS.resize(sub_n, std::vector<PT>(sub_n, PT(0)));

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
template<typename RT, typename PT, typename ST>
void Family<RT, PT, ST>::findSectors() {
    int total = 1 << numProps_;  // 2^numProps

    cases_.resize(total);

    for (int i = total - 1; i >= 0; i--) {
        std::vector<int> secvec = secvecFromIdx(i);
        std::vector<std::vector<PT>> subS;
        bool isSector = getSubS(secvec, subS);
        if (isSector) {
            Sector<RT, PT> tempSector(std::move(subS), nProps(secvec), numBranch_);
            int caseNum = tempSector.getCase();
            cases_[i] = caseNum;
            sectors_.push_back(std::move(tempSector));
            sectorIdxs_.push_back(i);
        }
        else {
            cases_[i] = -1;
        }
    }
    rebuildMasters();
}

template<typename RT, typename PT, typename ST>
void Family<RT, PT, ST>::rebuildMasters() {
    masterIdxs_.clear();
    for (int idx = static_cast<int>(cases_.size()) - 1; idx >= 0; --idx) {
        if (cases_[idx] != 0) continue;
        const std::vector<int> sector = secvecFromIdx(idx);
        if (sectorMap_.hasSource(sector)) continue;
        masterIdxs_.push_back(idx);
    }
}
