#pragma once

#include <vector>
#include "family.hpp"

template<typename T>
class FBI { 
private:
    const Family<T>* family;    // 指向所属的 family（所有FBI可共享）
    int sectorNu;               // 所属 sector 的 nu（整数形式）
    int numProps;               // N（传播子数量）
    std::vector<int> nu;        // FBI 的 nu 向量，长度为 N，元素非负
    T delta;                    // L * Dimension / 2

public:
    // 构造函数
    FBI(const Family<T>* fam, int sec_nu, const std::vector<int>& fbi_nu, T d)
        : family(fam), sectorNu(sec_nu), numProps(fbi_nu.size()), nu(fbi_nu), delta(d) {}
    
    // Getter 方法
    inline const Family<T>* getFamily() const { return family; }
    inline int getSectorNu() const { return sectorNu; }
    inline int getNumProps() const { return numProps; }
    inline const std::vector<int>& getNu() const { return nu; }
    inline const T& getDelta() const { return delta; }
    
    // 获取所属的 Sector
    const Sector<T>* getSector() const {
        return family->getSectorByNu(sectorNu);
    }
    
    // 检查 nu 的有效性（元素非负）
    bool isValidNu() const {
        for (int n : nu) {
            if (n < 0) return false;
        }
        return true;
    }
    
    // 计算 nu 的总和（可能用于排序或其他目的）
    int getNuSum() const {
        int sum = 0;
        for (int n : nu) {
            sum += n;
        }
        return sum;
    }
};

#include "../src/FBI.tpp"
