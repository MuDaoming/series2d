/**
 * Sector class - 双参数模板版本
 * 处理子矩阵S，进行RREF行化简，计算C和z系数
 * 
 * 模板参数：
 *   - RT: 有理函数类型（如 GiNaC::ex 或 Rational<FlintMod>）
 *   - PT: 多项式类型（如 GiNaC::ex 或 Polynomial<FlintMod>）
 * 
 * RT和PT都需要支持:
 *   - 从int构造（0和1）
 *   - +, -, *, / 运算符
 */

#pragma once

#include <vector>
#include <sstream>
#include <stdexcept>
#include <type_traits>

// 前向声明GiNaC::ex（避免在未使用GiNaC时引入头文件）
namespace GiNaC { class ex; }

// 前向声明Rational和Polynomial
template<typename T> struct Rational;
template<typename T> struct Polynomial;

// 类型萃取：判断是否为GiNaC::ex
template<typename T>
struct is_ginac_ex : std::is_same<T, GiNaC::ex> {};

// 类型萃取：判断是否为Rational<T>
template<typename T>
struct is_rational_type : std::false_type {};

template<typename T>
struct is_rational_type<Rational<T>> : std::true_type {};

// 类型萃取：判断是否为Polynomial<T>
template<typename T>
struct is_polynomial_type : std::false_type {};

template<typename T>
struct is_polynomial_type<Polynomial<T>> : std::true_type {};

// 零判断函数
template<typename T>
inline bool isZero(const T& val) {
    if constexpr (is_ginac_ex<T>::value) {
        // GiNaC::ex类型：使用normal().is_zero()
        return GiNaC::normal(val).is_zero();
    } else if constexpr (is_rational_type<T>::value) {
        // Rational<T>类型：检查分子的单项式数量是否为0
        return val.numerator.getNumOfMonomials() == 0;
    } else if constexpr (is_polynomial_type<T>::value) {
        // Polynomial<T>类型：检查单项式数量是否为0
        return val.getNumOfMonomials() == 0;
    } else {
        // 数值类型：使用 == T(0)
        return val == T(0);
    }
}

// normalize函数
template<typename T>
inline T normalize(const T& val) {
    if constexpr (is_ginac_ex<T>::value) {
        // GiNaC::ex类型：使用normal()约分
        return GiNaC::normal(val);
    } else {
        // 数值类型：直接返回
        return val;
    }
}

// lcm函数（最小公倍数）
template<typename T>
inline T lcm(const T& a, const T& b) {
    if constexpr (is_ginac_ex<T>::value) {
        // GiNaC::ex类型：使用GiNaC的lcm函数
        return GiNaC::lcm(a, b);
    } else {
        // 数值类型：对于有限域等类型，直接返回乘积（或根据具体类型实现）
        // 这里假设数值类型的lcm就是乘积（适用于互质情况或有限域）
        return a * b;
    }
}

// 计算向量中所有元素的lcm
template<typename T>
inline T lcmOfVector(const std::vector<T>& vec) {
    if (vec.empty()) {
        return T(1);
    }
    T result = vec[0];
    for (size_t i = 1; i < vec.size(); ++i) {
        if (!isZero(vec[i])) {
            result = lcm(result, vec[i]);
        }
    }
    return normalize(result);
}

// 获取有理函数的分母
template<typename T>
inline T getDenominator(const T& val) {
    if constexpr (is_ginac_ex<T>::value) {
        // GiNaC::ex类型：使用denom()获取分母
        return GiNaC::denom(GiNaC::normal(val));
    } else {
        // 数值类型：分母为1
        return T(1);
    }
}

// 获取有理函数的分子
template<typename T>
inline T getNumerator(const T& val) {
    if constexpr (is_ginac_ex<T>::value) {
        // GiNaC::ex类型：使用numer()获取分子
        return GiNaC::numer(GiNaC::normal(val));
    } else {
        // 数值类型：分子就是自身
        return val;
    }
}

// 计算向量中所有元素分母的lcm
template<typename T>
inline T lcmOfDenominators(const std::vector<T>& vec) {
    if (vec.empty()) {
        return T(1);
    }
    T result = getDenominator(vec[0]);
    for (size_t i = 1; i < vec.size(); ++i) {
        T denom = getDenominator(vec[i]);
        if (!isZero(denom)) {
            result = lcm(result, denom);
        }
    }
    return normalize(result);
}

// 通用求导函数（对指定符号变量求偏导）
template<typename T, typename Symbol>
inline T diff(const T& val, const Symbol& var) {
    if constexpr (is_ginac_ex<T>::value) {
        // GiNaC::ex类型：使用GiNaC的diff函数
        return GiNaC::diff(val, var);
    } else {
        // 数值类型：导数为0（常数对变量的导数）
        return T(0);
    }
}

/**
 * Sector类
 * @tparam RT 有理函数类型（用于candz_, C_, invS_等）
 * @tparam PT 多项式类型（用于S_矩阵元素）
 */
template<typename RT, typename PT>
class Sector {

    // 声明友元用于转换
    template<typename, typename> friend class Sector;
    template<typename, typename, typename> friend class Family;

public: 
    // 默认构造函数（用于转换）
    Sector() : numProps_(0), numBranch_(0), dimNull_(0), z0_(0) {}

    // Constructor
    Sector(const std::vector<std::vector<PT>>& S, int numProps, int numBranch) 
        : S_(S), numProps_(numProps), numBranch_(numBranch) 
    {
        candz_.resize(numBranch_ + numProps_);
        rowReduce();
        solveCandZ();
    }

    // 转换构造函数：从另一种类型的Sector转换
    // Converter是一个函数对象，用于类型转换
    template<typename OtherRT, typename OtherPT, typename PolyConverter, typename RatConverter>
    Sector(const Sector<OtherRT, OtherPT>& other, 
           PolyConverter polyConv, RatConverter ratConv)
        : numProps_(other.numProps_), numBranch_(other.numBranch_),
          dimNull_(other.dimNull_), z0_(other.z0_)
    {
        // 转换S_矩阵
        S_.resize(other.S_.size());
        for (size_t i = 0; i < other.S_.size(); ++i) {
            S_[i].resize(other.S_[i].size());
            for (size_t j = 0; j < other.S_[i].size(); ++j) {
                S_[i][j] = polyConv(other.S_[i][j]);
            }
        }
        
        // 转换reducedS_
        reducedS_.resize(other.reducedS_.size());
        for (size_t i = 0; i < other.reducedS_.size(); ++i) {
            reducedS_[i].resize(other.reducedS_[i].size());
            for (size_t j = 0; j < other.reducedS_[i].size(); ++j) {
                try {
                    reducedS_[i][j] = ratConv(other.reducedS_[i][j]);
                } catch (const std::exception& e) {
                    std::ostringstream oss;
                    oss << "Sector convert failed at reducedS[" << i << "][" << j << "], expr="
                        << other.reducedS_[i][j] << ", reason: " << e.what();
                    throw std::runtime_error(oss.str());
                }
            }
        }
        
        // 转换rowOperation_
        rowOperation_.resize(other.rowOperation_.size());
        for (size_t i = 0; i < other.rowOperation_.size(); ++i) {
            rowOperation_[i].resize(other.rowOperation_[i].size());
            for (size_t j = 0; j < other.rowOperation_[i].size(); ++j) {
                try {
                    rowOperation_[i][j] = ratConv(other.rowOperation_[i][j]);
                } catch (const std::exception& e) {
                    std::ostringstream oss;
                    oss << "Sector convert failed at rowOperation[" << i << "][" << j << "], expr="
                        << other.rowOperation_[i][j] << ", reason: " << e.what();
                    throw std::runtime_error(oss.str());
                }
            }
        }
        
        // 转换candz_
        candz_.resize(other.candz_.size());
        for (size_t i = 0; i < other.candz_.size(); ++i) {
            try {
                candz_[i] = ratConv(other.candz_[i]);
            } catch (const std::exception& e) {
                std::ostringstream oss;
                oss << "Sector convert failed at candz[" << i << "], expr="
                    << other.candz_[i] << ", reason: " << e.what();
                throw std::runtime_error(oss.str());
            }
        }
        try {
            C_ = ratConv(other.C_);
        } catch (const std::exception& e) {
            std::ostringstream oss;
            oss << "Sector convert failed at C, expr=" << other.C_ << ", reason: " << e.what();
            throw std::runtime_error(oss.str());
        }
        
        // 转换分子分母形式
        denoCandZ_ = polyConv(other.denoCandZ_);
        numeCandZ_.resize(other.numeCandZ_.size());
        for (size_t i = 0; i < other.numeCandZ_.size(); ++i) {
            numeCandZ_[i] = polyConv(other.numeCandZ_[i]);
        }
        numeC_ = polyConv(other.numeC_);
        
        // 转换denoRowOperation_和numeRowOperation_
        denoRowOperation_.resize(other.denoRowOperation_.size());
        for (size_t i = 0; i < other.denoRowOperation_.size(); ++i) {
            denoRowOperation_[i] = polyConv(other.denoRowOperation_[i]);
        }
        numeRowOperation_.resize(other.numeRowOperation_.size());
        for (size_t i = 0; i < other.numeRowOperation_.size(); ++i) {
            numeRowOperation_[i].resize(other.numeRowOperation_[i].size());
            for (size_t j = 0; j < other.numeRowOperation_[i].size(); ++j) {
                numeRowOperation_[i][j] = polyConv(other.numeRowOperation_[i][j]);
            }
        }
    }

    inline const RT& getCSum() const { return C_; }
    inline int getZ0() const { return z0_; }  // z0 是 0 或 1
    inline const RT& getC(int i) const { return candz_[i]; }
    inline const RT& getZ(int i) const { return candz_[numBranch_ + i]; }
    inline int getDimNull() const { return dimNull_; }
    inline const std::vector<std::vector<RT>>& getInvS() const { 
        if(dimNull_ == 0) {
            return rowOperation_;
        } else {
            throw std::runtime_error("InvS does not exist when dimNull > 0");
        } 
    }
    inline int getCase() const {
        if (dimNull_ == 0 && !isZero(C_)) return 0;
        else if (dimNull_ == 0 && isZero(C_)) return 1;
        else if (dimNull_ != 0 && !isZero(C_)) return 2;
        else return 3;
    }

    // 获取candz的引用（用于调试输出）
    const std::vector<RT>& getCandZ() const { return candz_; }
    int getNumBranch() const { return numBranch_; }
    int getNumProps() const { return numProps_; }

    // 获取分母和分子形式的candz
    inline const PT& getDenoCandZ() const { return denoCandZ_; }
    inline const std::vector<PT>& getNumeCandZ() const { return numeCandZ_; }
    inline const PT& getNumeC() const { return numeC_; }
    inline PT getNumeC(int i) const { return numeCandZ_[i]; }  // C_i * denoCandZ
    inline PT getNumeZ(int i) const { return numeCandZ_[numBranch_ + i]; }  // z_i * denoCandZ

    // 获取分母和分子形式的rowOperation（当dimNull == 0时为invS）
    inline const std::vector<PT>& getDenoRowOperation() const { return denoRowOperation_; }
    inline const std::vector<std::vector<PT>>& getNumeRowOperation() const { return numeRowOperation_; }
    inline const std::vector<std::vector<PT>>& getNumeInvS() const {
        if(dimNull_ == 0) {
            return numeRowOperation_;
        } else {
            throw std::runtime_error("NumeInvS does not exist when dimNull > 0");
        }
    }
    inline const std::vector<PT>& getDenoInvS() const {
        if(dimNull_ == 0) {
            return denoRowOperation_;
        } else {
            throw std::runtime_error("DenoInvS does not exist when dimNull > 0");
        }
    }

private:
    // from construction
    std::vector<std::vector<PT>> S_; // size (N + B)*(N + B), 多项式矩阵
    int numProps_;  // N
    int numBranch_; // B

    // RREF, stored the following variables in row reduction
    std::vector<std::vector<RT>> reducedS_;      // RREF后的矩阵（有理函数）
    std::vector<std::vector<RT>> rowOperation_;  // 行变换矩阵（dimNull=0时为invS）
    int dimNull_;
    int z0_;  // z0 是 0 或 1（整数）

    // solved candz
    std::vector<RT> candz_;     // size N + B，有理函数
    RT C_;                      // C = C_1 + ... + C_B，有理函数

    // 分母和分子形式的candz（避免有理函数运算）
    PT denoCandZ_;                      // candz_所有元素分母的lcm
    std::vector<PT> numeCandZ_;         // numeCandZ_[i] = candz_[i] * denoCandZ_
    PT numeC_;                          // numeC_ = C_ * denoCandZ_

    // 分母和分子形式的rowOperation（每行单独处理）
    std::vector<PT> denoRowOperation_;              // 每行的分母lcm
    std::vector<std::vector<PT>> numeRowOperation_; // numeRowOperation_[i][j] = rowOperation_[i][j] * denoRowOperation_[i]

    // RREF
    void rowReduce();
    std::vector<std::vector<RT>> findNullSpace() const;
    std::vector<RT> solveLinear(const std::vector<RT>& b) const;
    void solveCandZ();

};

#include "../src/sector.tpp"
