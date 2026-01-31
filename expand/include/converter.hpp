/**
 * converter.hpp - GiNaC到FlintMod类型的转换函数
 * 
 * 提供以下转换功能：
 *   - GiNaC::ex (多项式) -> Polynomial<FlintMod>
 *   - GiNaC::ex (有理函数) -> Rational<FlintMod>
 *   - GiNaC::ex (标量) -> FlintMod
 *   - Family<ex, ex, ex> -> Family<Rational<FlintMod>, Polynomial<FlintMod>, FlintMod>
 */

#pragma once

#include <sstream>
#include <ginac/ginac.h>
#include "ffType.hpp"
#include "rational.hpp"
#include "family.hpp"

/**
 * 将GiNaC::ex多项式转换为Polynomial<FlintMod>
 * @param expr GiNaC表达式（必须是关于X,Y的多项式）
 * @param X X符号变量
 * @param Y Y符号变量
 * @return Polynomial<FlintMod>
 */
inline Polynomial<FlintMod> exToPolynomial(const GiNaC::ex& expr, 
                                            const GiNaC::symbol& X, 
                                            const GiNaC::symbol& Y) {
    using namespace GiNaC;
    
    Polynomial<FlintMod> result;
    
    // 展开表达式
    ex expanded = expand(expr);
    
    // 获取X和Y的最大度数
    int maxDegX = expanded.degree(X);
    int maxDegY = expanded.degree(Y);
    
    // 遍历所有可能的(i,j)组合
    for (int i = 0; i <= maxDegX; ++i) {
        for (int j = 0; j <= maxDegY; ++j) {
            // 提取X^i * Y^j 的系数
            ex coeff = expanded.coeff(X, i).coeff(Y, j);
            
            // 如果系数是数值（不包含X或Y）
            if (is_a<numeric>(coeff)) {
                numeric num = ex_to<numeric>(coeff);
                // 转换为整数（如果是有理数，需要在有限域中处理）
                if (num.is_integer()) {
                    long val = num.to_long();
                    if (val != 0) {
                        result.addMonomial(FlintMod(val), Power(i, j));
                    }
                } else if (num.is_rational()) {
                    // 有理数系数：numer/denom，在有限域中计算
                    long numer = num.numer().to_long();
                    long denom = num.denom().to_long();
                    FlintMod val = FlintMod(numer) / FlintMod(denom);
                    result.addMonomial(val, Power(i, j));
                }
            } else if (!coeff.is_zero()) {
                // 系数应该是纯数值，如果不是则抛出异常
                std::ostringstream oss;
                oss << coeff;
                throw std::runtime_error("Coefficient is not numeric: " + oss.str());
            }
        }
    }
    
    return result;
}

/**
 * 将GiNaC::ex有理函数转换为Rational<FlintMod>
 * @param expr GiNaC表达式（有理函数）
 * @param X X符号变量
 * @param Y Y符号变量
 * @return Rational<FlintMod>
 */
inline Rational<FlintMod> exToRational(const GiNaC::ex& expr,
                                        const GiNaC::symbol& X,
                                        const GiNaC::symbol& Y) {
    using namespace GiNaC;
    
    // 规范化表达式
    ex normalized = normal(expr);
    
    // 获取分子和分母
    ex numer_ex = numer(normalized);
    ex denom_ex = denom(normalized);
    
    // 转换分子和分母为Polynomial<FlintMod>
    Polynomial<FlintMod> numer_poly = exToPolynomial(numer_ex, X, Y);
    Polynomial<FlintMod> denom_poly = exToPolynomial(denom_ex, X, Y);
    
    return Rational<FlintMod>(numer_poly, denom_poly);
}

/**
 * 将GiNaC::ex标量转换为FlintMod
 * @param expr GiNaC表达式（必须是数值）
 * @return FlintMod
 */
inline FlintMod exToScalar(const GiNaC::ex& expr) {
    using namespace GiNaC;
    
    if (!is_a<numeric>(expr)) {
        std::ostringstream oss;
        oss << expr;
        throw std::runtime_error("Expression is not numeric: " + oss.str());
    }
    
    numeric num = ex_to<numeric>(expr);
    if (num.is_integer()) {
        return FlintMod(num.to_long());
    } else if (num.is_rational()) {
        return FlintMod(num.numer().to_long()) / FlintMod(num.denom().to_long());
    }
    
    return FlintMod(0);
}

/**
 * 转换器类：封装X和Y符号变量的转换函数
 */
class GiNaCToFlintConverter {
public:
    GiNaCToFlintConverter(const GiNaC::symbol& X, const GiNaC::symbol& Y)
        : X_(X), Y_(Y) {}
    
    // 多项式转换器
    Polynomial<FlintMod> convertPoly(const GiNaC::ex& expr) const {
        return exToPolynomial(expr, X_, Y_);
    }
    
    // 有理函数转换器
    Rational<FlintMod> convertRat(const GiNaC::ex& expr) const {
        return exToRational(expr, X_, Y_);
    }
    
    // 标量转换器
    FlintMod convertScalar(const GiNaC::ex& expr) const {
        return exToScalar(expr);
    }
    
private:
    const GiNaC::symbol& X_;
    const GiNaC::symbol& Y_;
};

/**
 * 从Family<GiNaC::ex, GiNaC::ex, GiNaC::ex>
 * 转换为Family<Rational<FlintMod>, Polynomial<FlintMod>, FlintMod>
 * 
 * 注意：两者都是含有符号变量X,Y的表达式，
 * 只是存储类型不同（GiNaC vs FlintMod多项式）
 * 
 * @param ginacFamily 使用GiNaC表示的Family
 * @param X X符号变量
 * @param Y Y符号变量
 * @return Family<Rational<FlintMod>, Polynomial<FlintMod>, FlintMod>
 */
inline Family<Rational<FlintMod>, Polynomial<FlintMod>, FlintMod> 
convertFamily(const Family<GiNaC::ex, GiNaC::ex, GiNaC::ex>& ginacFamily,
              const GiNaC::symbol& X,
              const GiNaC::symbol& Y) {
    GiNaCToFlintConverter converter(X, Y);
    
    // 使用lambda包装converter方法
    auto polyConv = [&converter](const GiNaC::ex& e) { return converter.convertPoly(e); };
    auto ratConv = [&converter](const GiNaC::ex& e) { return converter.convertRat(e); };
    auto scalarConv = [&converter](const GiNaC::ex& e) { return converter.convertScalar(e); };
    
    // 调用Family的转换构造函数
    return Family<Rational<FlintMod>, Polynomial<FlintMod>, FlintMod>(
        ginacFamily, polyConv, ratConv, scalarConv);
}
