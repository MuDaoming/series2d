#pragma once

#include <string>
#include <string_view>

/// Token类型枚举
enum class TokenType {
    NUMBER,        ///< 数字：123, -456
    VARIABLE,      ///< 变量：X, Y  
    OPERATOR,      ///< 运算符：+, -, *, ^, /
    LPAREN,        ///< 左括号：(
    RPAREN,        ///< 右括号：)
    LBRACE,        ///< 左大括号：{
    RBRACE,        ///< 右大括号：}
    COMMA,         ///< 逗号：,
    EOF_TOKEN,     ///< 文件结束
    INVALID        ///< 无效token
};

/// 词法单元（零拷贝设计）
struct Token {
    TokenType type;         ///< Token类型
    const char* start;      ///< 指向原始字符串位置
    size_t length;          ///< Token长度
    int line;              ///< 行号（1-based）
    int column;            ///< 列号（1-based）
    
    /// 默认构造函数
    Token() : type(TokenType::EOF_TOKEN), start(nullptr), length(0), line(0), column(0) {}
    
    /// 构造函数
    Token(TokenType t, const char* s, size_t len, int l, int c) 
        : type(t), start(s), length(len), line(l), column(c) {}
    
    /// 获取Token文本（零拷贝）
    std::string_view getText() const {
        return start ? std::string_view(start, length) : std::string_view{};
    }
    
    /// 检查Token是否为空或无效
    bool isEmpty() const {
        return start == nullptr || length == 0 || type == TokenType::INVALID;
    }
    
    /// 检查Token是否为EOF
    bool isEOF() const {
        return type == TokenType::EOF_TOKEN;
    }
    
    /// 检查Token是否为特定类型
    bool is(TokenType t) const {
        return type == t;
    }
    
    /// 转换为字符串（仅调试用）
    std::string toString() const;
    
    /// 获取Token的位置信息
    std::string getPositionString() const {
        return "line " + std::to_string(line) + ", column " + std::to_string(column);
    }
};

/// Token类型转字符串（调试用）
const char* tokenTypeToString(TokenType type);

/// 检查字符是否为运算符
bool isOperatorChar(char c);

/// 字符转换为对应的运算符Token类型
TokenType charToOperatorType(char c);


#include "../src/token.tpp"