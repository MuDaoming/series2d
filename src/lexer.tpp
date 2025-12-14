#include <cctype>
#include <stdexcept>

Lexer::Lexer() : pos_(0), line_(1), column_(1) {}

std::vector<Token> Lexer::tokenize(const std::string& input) {
    input_ = input.c_str();
    pos_ = 0;
    length_ = input.length();
    line_ = 1;
    column_ = 1;
    tokens_.clear();
    
    while (pos_ < length_) {
        skipWhitespace();
        if (pos_ >= length_) break;
        
        Token token = nextToken();
        if (token.type == TokenType::EOF_TOKEN) break;
        tokens_.push_back(token);
    }
    
    tokens_.emplace_back(TokenType::EOF_TOKEN, nullptr, 0, line_, column_);
    return tokens_;
}

Token Lexer::nextToken() {
    skipWhitespace();
    
    if (pos_ >= length_) {
        return Token(TokenType::EOF_TOKEN, nullptr, 0, line_, column_);
    }
    
    char c = getCurrentChar();
    const char* start = input_ + pos_;
    int startLine = line_;
    int startCol = column_;
    
    // 数字
    if (std::isdigit(c)) {
        return readNumber();
    }
    
    // 变量
    if (std::isalpha(c)) {
        return readVariable();
    }
    
    // 单字符token
    advance();
    return Token(charToOperatorType(c), start, 1, startLine, startCol);
}

Token Lexer::readNumber() {
    const char* start = input_ + pos_;
    int startLine = line_;
    int startCol = column_;
    size_t len = 0;
    
    while (pos_ < length_ && std::isdigit(getCurrentChar())) {
        advance();
        len++;
    }
    
    return Token(TokenType::NUMBER, start, len, startLine, startCol);
}

Token Lexer::readVariable() {
    const char* start = input_ + pos_;
    int startLine = line_;
    int startCol = column_;
    size_t len = 0;
    
    while (pos_ < length_ && std::isalnum(getCurrentChar())) {
        advance();
        len++;
    }
    
    return Token(TokenType::VARIABLE, start, len, startLine, startCol);
}

void Lexer::skipWhitespace() {
    while (pos_ < length_) {
        char c = getCurrentChar();
        if (c == ' ' || c == '\t') {
            advance();
        } else if (c == '\n') {
            line_++;
            column_ = 1;
            pos_++;
        } else if (c == '\r') {
            // 处理 \r\n 或单独的 \r
            pos_++;
            if (pos_ < length_ && input_[pos_] == '\n') {
                pos_++;
            }
            line_++;
            column_ = 1;
        } else {
            break;
        }
    }
}

char Lexer::getCurrentChar() const {
    return (pos_ < length_) ? input_[pos_] : '\0';
}

void Lexer::advance() {
    if (pos_ < length_) {
        pos_++;
        column_++;
    }
}
