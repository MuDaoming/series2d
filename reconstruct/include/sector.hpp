#pragma once

#include <vector>
#include <stdexcept>

template<typename T>
class Sector { 

public: 
    // Constructor
    Sector(const std::vector<std::vector<T>>& S, int numProps, int numBranch) : S_(S), numProps_(numProps), numBranch_(numBranch) {
    candz_.resize(numBranch_ + numProps_);
    rowReduce();
    solveCandZ();
    }

    inline const T& getCSum() const { return C_; }
    inline const T& getZ0() const { return z0_; }
    inline const T& getC(int i) const { return candz_[i]; }
    inline const T& getZ(int i) const { return candz_[numBranch_ + i]; }
    inline const int getDimNull() const { return dimNull_; }
    inline const std::vector<std::vector<T>>& getInvS() const { if(dimNull_ == 0) {return rowOperation_;} else {throw std::runtime_error("InvS does not exist when dimNull > 0");} }
    inline int getCase() const {
        if (dimNull_ == 0 && C_ != T(0)) return 0;
        else if (dimNull_ == 0 && C_ == T(0)) return 1;
        else if (dimNull_ != 0 && C_ != T(0)) return 2;
        else return 3;
    }


private:
    // from construction
    std::vector<std::vector<T>> S_; // size (N + B)*(N + B)
    int numProps_;  // N
    int numBranch_; // B

    // RREF, stored the following variables in row reduction
    std::vector<std::vector<T>> reducedS_;
    std::vector<std::vector<T>> rowOperation_; // note that if dimNull == 0, rowOperation is invS
    int dimNull_;
    T z0_;

    // solved candz
    std::vector<T> candz_;     // size N + B
    T C_;                     // C = C_1 + ... + C_B

    // RREF
    void rowReduce(); // call rowReduce at the end of construction
    std::vector<std::vector<T>> findNullSpace() const; // finds basis for null space if dimNull >=1
    std::vector<T> solveLinear(const std::vector<T>& b) const; // if dimNull == 0, solves S x = b via x = invS * b
    void solveCandZ(); // solveCandZ by calling solveLinear or findNullSpace depending on dimNull

};



#include "../src/sector.tpp"
