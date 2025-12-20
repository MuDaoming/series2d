#pragma once

#include <vector>
#include <stdexcept>

template<typename T>
class Sector { 
private:
    // from construction
    std::vector<std::vector<T>> S; // size (N + B)*(N + B)
    int numProps;  // N
    int numBranch; // B
    std::vector<int> nu; // sector nu vector, size N, elements are 0 or 1
    int Idx; // sector nu as integer


    // RREF, stored the following variables in row reduction
    std::vector<std::vector<T>> reducedS;
    std::vector<std::vector<T>> rowOperation; // note that if dimNull == 0, rowOperation is invS
    int dimNull;
    T z0;                    // z0 = 1 if dimNull == 0, else z0 = 0


    // if dimNull == 0; S (C_1,...,C_B,z_1,...,z_N) = (1_1,...,1_B,0,...,0)
    // if dimNull == 1; S (C_1,...,C_B,z_1,...,z_N) = (0,...,0,0,...,0), and there is only one solution up to scaling
    // if dimNull >=2; S (C_1,...,C_B,z_1,...,z_N) = (0,...,0,0,...,0), and there are multiple solutions, choose a solution such that C is as nonzero as possible
    std::vector<T> candz; // size N + B
    T C;                     // C = C_1 + ... + C_B

public: 
    std::vector<int> nuFromIdx(int idx, int numProps) const {
        std::vector<int> nu(numProps);
        for (int i = numProps - 1; i >= 0; i--) {
            nu[i] = idx % 2;
            idx /= 2;
        }
        return nu;
    }
    int idxFromNu(const std::vector<int>& nu) const {
        int idx = 0;
        for (size_t i = 0; i < nu.size(); i++) {
            idx = idx * 2 + (nu[i] > 0 ? 1 : 0);
        }
        return idx;
    }
    // Constructor
    Sector(const std::vector<std::vector<T>>& matrix, int n_props, int n_branch);
    
    // RREF
    void rowReduce();          // call rowReduce at the end of construction

    
    std::vector<std::vector<T>> findNullSpace() const; // finds basis for null space if dimNull >=1
    std::vector<T> solveLinear(const std::vector<T>& b) const; // if dimNull == 0, solves S x = b via x = invS * b

    void solveCandZ(); // solveCandZ by calling solveLinear or findNullSpace depending on dimNull
    
    inline const T& getCSum() const {
        return C;
    }

    inline const T& getZ0() const {
        return z0;
    }

    inline const T& getC(int i) const { 
        // if (i < 0 || i >= numBranch) throw std::out_of_range("C index out of range");
        return candz[i]; 
    }

    inline const T& getZ(int i) const { 
        // if (i < 0 || i >= numProps) throw std::out_of_range("Z index out of range");
        return candz[numBranch + i];
    }

    inline const int getDimNull() const {
        return dimNull;
    }

    // 添加用于测试的getter方法
    inline const std::vector<std::vector<T>>& getReducedS() const {
        return reducedS;
    }

    inline const std::vector<std::vector<T>>& getRowOperation() const {
        return rowOperation;
    }

    inline const std::vector<T>& getCandZ() const {
        return candz;
    }

    inline int getCase() const {
        // dimNull == 0 && C != 0 case 0
        if (dimNull == 0 && C != T(0)) return 0;
        // dimNull == 0 && C == 0 case 1
        else if (dimNull == 0 && C == T(0)) return 1;
        // dimNull !=0 && C !=0 case 2
        else if (dimNull != 0 && C != T(0)) return 2;
        // dimNull !=0 && C ==0 case 3
        else return 3;
    }

    std::vector<std::vector<int>> getBranch() const;

    // 判断删除指定的传播子后，是否有某个branch的所有元素都被删除
    static bool isDeleteBranch(const std::vector<std::vector<int>>& branch, const std::vector<int>& props);
    
    // 从 topS 和 nu 得到 subS（静态方法）
    static std::vector<std::vector<T>> getSubS(const std::vector<std::vector<T>>& topS, 
                                                const std::vector<int>& nu, 
                                                int numBranch);

};

#include "../src/sector.tpp"
