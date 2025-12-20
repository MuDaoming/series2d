// Family 类构造函数
template<typename T>
Family<T>::Family(const std::vector<std::vector<T>>& top_s, int n_props, int n_branch)
    : topS(top_s), numProps(n_props), numBranch(n_branch) {
    // 通过 findSectors 构造所有有效的 sectors
    findSectors();
}

// 整数转 nu 向量
template<typename T>
std::vector<int> Family<T>::intToNu(int n) const {
    std::vector<int> nu(numProps);
    for (int i = numProps - 1; i >= 0; i--) {
        nu[i] = n % 2;
        n /= 2;
    }
    return nu;
}

// nu 向量转整数
template<typename T>
int Family<T>::nuToInt(const std::vector<int>& nu) const {
    int result = 0;
    for (int i = 0; i < numProps; i++) {
        result = result * 2 + nu[i];
    }
    return result;
}

// 构造所有有效的 sectors
template<typename T>
void Family<T>::findSectors() {
    int total = 1 << numProps;  // 2^numProps
    
    // 初始化 cases，4个case
    cases.resize(4);
    
    // 首先从 topS 获取原始的 branch 信息（所有传播子都在的情况）
    // 创建一个全1的 nu 向量来获取完整的 branch 信息
    std::vector<int> fullNu(numProps, 1);
    Sector<T> topSector(topS, numProps, numBranch);
    std::vector<std::vector<int>> topBranch = topSector.getBranch();
    
    // 从 1 到 total-1（跳过全0的情况）
    for (int i = 1; i < total; i++) {
        std::vector<int> nu = intToNu(i);
        
        // 构造被删除的传播子列表（nu 中为 0 的）
        std::vector<int> deletedProps;
        for (int j = 0; j < numProps; j++) {
            if (nu[j] == 0) {
                deletedProps.push_back(j);
            }
        }
        
        // 使用 topBranch 检查是否有 branch 被完全删除
        bool branchDeleted = Sector<T>::isDeleteBranch(topBranch, deletedProps);
        
        // 如果没有 branch 被完全删除，则这是一个有效的 sector
        if (!branchDeleted) {
            // 从 topS 和 nu 得到 subS
            std::vector<std::vector<T>> subS = Sector<T>::getSubS(topS, nu, numBranch);
            
            // 创建 sector
            Sector<T> tempSector(subS, std::count(nu.begin(), nu.end(), 1), numBranch);
            
            // 获取该 sector 的 case 并添加到对应的 cases 中
            int caseNum = tempSector.getCase();
            cases[caseNum].push_back(i);
            
            sectors.push_back(std::move(tempSector));
            sectorNus.push_back(i);
        }
    }
}

// 根据 nu（整数）获取对应的 sector
template<typename T>
const Sector<T>* Family<T>::getSectorByNu(int nu) const {
    for (size_t i = 0; i < sectorNus.size(); i++) {
        if (sectorNus[i] == nu) {
            return &sectors[i];
        }
    }
    return nullptr;  // 未找到
}
