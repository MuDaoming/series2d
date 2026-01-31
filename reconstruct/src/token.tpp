#include <sstream>

std::string Token::toString() const {
    std::ostringstream oss;
    oss << "Token(" << tokenTypeToString(type) << ", \"" << getText() << "\", " 
        << line << ":" << column << ")";
    return oss.str();
}

const char* tokenTypeToString(TokenType type) {
    switch (type) {
        case TokenType::NUMBER:    return "NUMBER";
        case TokenType::VARIABLE:  return "VARIABLE";
        case TokenType::OPERATOR:  return "OPERATOR";
        case TokenType::LPAREN:    return "LPAREN";
        case TokenType::RPAREN:    return "RPAREN";
        case TokenType::LBRACE:    return "LBRACE";
        case TokenType::RBRACE:    return "RBRACE";
        case TokenType::COMMA:     return "COMMA";
        case TokenType::EOF_TOKEN: return "EOF";
        case TokenType::INVALID:   return "INVALID";
        default:                   return "UNKNOWN";
    }
}

bool isOperatorChar(char c) {
    return c == '+' || c == '-' || c == '*' || c == '/' || c == '^';
}

TokenType charToOperatorType(char c) {
    switch (c) {
        case '+':
        case '-':
        case '*':
        case '/':
        case '^':
            return TokenType::OPERATOR;
        case '(':
            return TokenType::LPAREN;
        case ')':
            return TokenType::RPAREN;
        case '{':
            return TokenType::LBRACE;
        case '}':
            return TokenType::RBRACE;
        case ',':
            return TokenType::COMMA;
        default:
            return TokenType::INVALID;
    }
}
