#pragma once

#include "token.hpp"
#include <vector>
#include <string>

/// 词法分析器类
class Lexer {
private:
    const char* input_;            ///< 输入字符串
    size_t pos_;                   ///< 当前解析位置
    size_t length_;                ///< 输入长度
    size_t line_;                  ///< 当前行号
    size_t column_;                ///< 当前列号
    std::vector<Token> tokens_;    ///< 已解析的tokens

public:
    /// 构造函数
    Lexer();

    /// 对输入字符串进行词法分析
    /// @param input 输入字符串
    /// @return 解析得到的Token向量
    std::vector<Token> tokenize(const std::string& input);

private:
    /// 获取下一个token
    Token nextToken();
    
    /// 解析数字
    Token readNumber();
    
    /// 解析变量名
    Token readVariable();
    
    /// 跳过空白字符
    void skipWhitespace();
    
    /// 获取当前字符
    char getCurrentChar() const;
    
    /// 前进一个字符
    void advance();
};

#include "../src/lexer.tpp"
