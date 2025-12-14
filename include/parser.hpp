#pragma once

#include "token.hpp"
#include "rational.hpp"
#include "lexer.hpp"
#include <string>
#include <stdexcept>
#include <vector>
#include <type_traits>
#include <string_view>
#include <charconv>

// ========== 独立工具函数 ==========

/// 将字符串转换为数字类型T
/// @tparam T 目标数字类型
/// @param text 要解析的字符串视图
/// @return 转换后的数字
/// @throws std::invalid_argument 如果字符串格式无效
/// @throws std::out_of_range 如果数字超出类型范围
template<typename T>
T stringToNumber(std::string_view text);

/// 解析异常类
class ParseError : public std::runtime_error {
private:
    int line_;
    int column_;
    
public:
    ParseError(const std::string& message, int line, int column)
        : std::runtime_error(message), line_(line), column_(column) {}
    
    int getLine() const { return line_; }
    int getColumn() const { return column_; }
    
    std::string getFormattedMessage() const {
        return "Parse error at line " + std::to_string(line_) + 
               ", column " + std::to_string(column_) + ": " + what();
    }
};

/// 语法分析器模板类
/// @tparam T 系数类型（整数、有限域等）
template<typename T>
class Parser {
private:
    Lexer lexer_;              ///< 词法分析器
    std::vector<Token> tokens_; ///< Token列表
    size_t token_index_;       ///< 当前Token索引
    Token current_token_;      ///< 当前token

    // ========== 语法分析方法 ==========
    
    /// 解析有理函数
    /// Grammar: Rational ::= Polynomial ('/' Polynomial)?
    Rational<T> parseRational();
    
    /// 解析多项式
    /// Grammar: Polynomial ::= ['(' ] Monomial (('+' | '-') Monomial)* [ ')' ]
    Polynomial<T> parsePolynomial();
    
    /// 解析单项式
    /// Grammar: Monomial ::= ['-'] Factor ('*' Factor)*
    /// 返回系数和幂次对
    std::pair<T, Power> parseMonomial();
    
    /// 解析因子
    /// Grammar: Factor ::= Number | Variable ('^' Number)?
    std::pair<T, Power> parseFactor();
    
    /// 解析系数（数字）
    T parseCoefficient();
    
    
    // ========== 辅助函数 ==========
    
    /// 将token转换为系数
    T stringToCoefficient(const Token& token);
    
    /// 从 long long 转换为类型 T
    T convertFromLongLong(long long value);
    
    /// 期望特定类型的token
    void expect(TokenType expected);
    
    /// 消费当前token，如果类型匹配
    bool consume(TokenType expected);
    
    /// 检查当前token类型
    bool check(TokenType type) const;
    
    /// 获取下一个token
    void nextToken();
    
    /// 抛出解析错误
    [[noreturn]] void error(const std::string& message);

public:
    /// 构造函数
    Parser();
    
    /// 解析矩阵的主入口
    /// Grammar: Matrix ::= '{' '{' RationalList '}' (',' '{' RationalList '}')* '}'
    Matrix<T> parseMatrix(const std::string& input);
};

#include "../src/parser.tpp"
