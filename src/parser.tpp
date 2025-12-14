

template<typename T>
Parser<T>::Parser() : token_index_(0) {}

// ========== 独立工具函数实现 ==========

template<typename T>
T stringToNumber(std::string_view text) {
    if (text.empty()) {
        throw std::invalid_argument("Empty string cannot be converted to number");
    }
    
    // 1. 对于自带整数类型，直接使用 std::from_chars
    if constexpr (std::is_integral_v<T>) {
        T result;
        auto [ptr, ec] = std::from_chars(text.data(), text.data() + text.size(), result);
        
        if (ec == std::errc::invalid_argument) {
            throw std::invalid_argument("Invalid number format: " + std::string(text));
        } else if (ec == std::errc::result_out_of_range) {
            throw std::out_of_range("Number out of range for type: " + std::string(text));
        } else if (ptr != text.data() + text.size()) {
            throw std::invalid_argument("Invalid number format (extra characters): " + std::string(text));
        }
        
        return result;
    }
    // 2. 对于自带浮点数类型，先读取为long long，再转化成浮点类型
    else if constexpr (std::is_floating_point_v<T>) {
        long long intermediate;
        auto [ptr, ec] = std::from_chars(text.data(), text.data() + text.size(), intermediate);
        
        if (ec == std::errc::invalid_argument) {
            throw std::invalid_argument("Invalid number format: " + std::string(text));
        } else if (ec == std::errc::result_out_of_range) {
            throw std::out_of_range("Number out of range for long long: " + std::string(text));
        } else if (ptr != text.data() + text.size()) {
            throw std::invalid_argument("Invalid number format (extra characters): " + std::string(text));
        }
        
        return static_cast<T>(intermediate);
    }
    // 3. 对于其它类型，先读取为unsigned long long，再转换为T()对应类型
    else {
        unsigned long long intermediate;
        auto [ptr, ec] = std::from_chars(text.data(), text.data() + text.size(), intermediate);
        
        if (ec == std::errc::invalid_argument) {
            throw std::invalid_argument("Invalid number format: " + std::string(text));
        } else if (ec == std::errc::result_out_of_range) {
            throw std::out_of_range("Number out of range for unsigned long long: " + std::string(text));
        } else if (ptr != text.data() + text.size()) {
            throw std::invalid_argument("Invalid number format (extra characters): " + std::string(text));
        }
        
        return T(intermediate);
    }
}

template<typename T>
Matrix<T> Parser<T>::parseMatrix(const std::string& input) {
    // 使用词法分析器进行分词
    tokens_ = lexer_.tokenize(input);
    token_index_ = 0;
    
    // if (tokens_.empty()) {
    //     throw ParseError("Empty input", 1, 1);
    // }
    
    current_token_ = tokens_[0];
    
    Matrix<T> matrix;
    
    // 解析最外层大括号
    expect(TokenType::LBRACE);
    
    while (!check(TokenType::RBRACE)) {
        // 解析内层大括号
        expect(TokenType::LBRACE);
        
        std::vector<Rational<T>> row;
        
        // 解析第一个有理函数
        auto element = parseRational();
        row.push_back(element);
        
        // 解析后续有理函数
        while (check(TokenType::COMMA)) {
            nextToken(); // 消费逗号
            auto element = parseRational();
            row.push_back(element);
        }
        
        expect(TokenType::RBRACE);
        matrix.push_back(row);
        
        // 检查是否有更多行
        if (check(TokenType::COMMA)) {
            nextToken(); // 消费逗号
        }
    }
    
    // 期望最外层的 '}'
    expect(TokenType::RBRACE);
    expect(TokenType::EOF_TOKEN);
    
    return matrix;
}


template<typename T>
Rational<T> Parser<T>::parseRational() {
    auto numerator = parsePolynomial();
    
    if (check(TokenType::OPERATOR) && current_token_.getText() == "/") {
        consume(TokenType::OPERATOR); // consume '/'
        auto denominator = parsePolynomial();
        return Rational<T>(numerator, denominator);
    }
    
    // 分母为1的情况
    Polynomial<T> denominator;
    denominator.addMonomial(T{1}, Power(0, 0));
    return Rational<T>(numerator, denominator);
}

template<typename T>
Polynomial<T> Parser<T>::parsePolynomial() {
    Polynomial<T> poly;
    
    // 检查是否有开括号
    bool hasParentheses = false;
    if (check(TokenType::LPAREN)) {
        hasParentheses = true;
        consume(TokenType::LPAREN);
    }
    
    // 解析第一个单项式
    auto firstMono = parseMonomial();
    poly.addMonomial(firstMono.first, firstMono.second);
    
    // 解析后续单项式
    while (check(TokenType::OPERATOR) && 
           (current_token_.getText() == "+" || current_token_.getText() == "-")) {
        
        bool isNegative = (current_token_.getText() == "-");
        consume(TokenType::OPERATOR);
        
        auto mono = parseMonomial();
        T coeff = isNegative ? -mono.first : mono.first;
        poly.addMonomial(coeff, mono.second);
    }
    
    // 如果有开括号，检查闭括号
    if (hasParentheses) {
        expect(TokenType::RPAREN);
    }
    
    return poly;
}

template<typename T>
std::pair<T, Power> Parser<T>::parseMonomial() {
    T coefficient = T{1};
    Power power(0, 0);
    
    // 处理可选的负号
    bool isNegative = false;
    if (check(TokenType::OPERATOR) && current_token_.getText() == "-") {
        isNegative = true;
        consume(TokenType::OPERATOR);
    }
    
    // 解析第一个因子
    auto firstFactor = parseFactor();
    coefficient *= firstFactor.first;
    power.x_power += firstFactor.second.x_power;
    power.y_power += firstFactor.second.y_power;
    
    // 解析后续因子
    while (check(TokenType::OPERATOR) && current_token_.getText() == "*") {
        consume(TokenType::OPERATOR);
        auto factor = parseFactor();
        coefficient *= factor.first;
        power.x_power += factor.second.x_power;
        power.y_power += factor.second.y_power;
    }
    
    if (isNegative) {
        coefficient = -coefficient;
    }
    
    return {coefficient, power};
}

template<typename T>
std::pair<T, Power> Parser<T>::parseFactor() {
    if (check(TokenType::NUMBER)) {
        T coeff = stringToCoefficient(current_token_);
        consume(TokenType::NUMBER);
        return {coeff, Power(0, 0)};
    }
    
    if (check(TokenType::VARIABLE)) {
        std::string_view varName = current_token_.getText();
        consume(TokenType::VARIABLE);
        
        int exponent = 1;
        if (check(TokenType::OPERATOR) && current_token_.getText() == "^") {
            consume(TokenType::OPERATOR);
            if (!check(TokenType::NUMBER)) {
                error("Expected number after '^'");
            }
            // 指数直接使用stringToNumber解析，避免字符串构造
            exponent = stringToNumber<int>(current_token_.getText());
            consume(TokenType::NUMBER);
        }
        
        if (varName == "X") {
            return {T{1}, Power(exponent, 0)};
        } else if (varName == "Y") {
            return {T{1}, Power(0, exponent)};
        } else {
            error("Unknown variable: " + std::string(varName));
        }
    }
    
    error("Expected number or variable");
}

template<typename T>
T Parser<T>::parseCoefficient() {
    expect(TokenType::NUMBER);
    return stringToCoefficient(tokens_[token_index_ - 1]);
}

template<typename T>
T Parser<T>::convertFromLongLong(long long value) {
        T result = T(value);
        return result;
}

template<typename T>
T Parser<T>::stringToCoefficient(const Token& token) {
    std::string_view text = token.getText();
    if (text.empty()) {
        error("Empty number token");
    }
    
    try {
        return stringToNumber<T>(text);
    } catch (const std::invalid_argument& e) {
        error(e.what());
    } catch (const std::out_of_range& e) {
        error(e.what());
    }
}

template<typename T>
void Parser<T>::expect(TokenType expected) {
    if (!check(expected)) {
        error("Expected " + std::string(tokenTypeToString(expected)) + 
              " but got " + std::string(tokenTypeToString(current_token_.type)));
    }
    nextToken();
}

template<typename T>
bool Parser<T>::consume(TokenType expected) {
    if (check(expected)) {
        nextToken();
        return true;
    }
    return false;
}

template<typename T>
bool Parser<T>::check(TokenType type) const {
    return current_token_.type == type;
}

template<typename T>
void Parser<T>::nextToken() {
    if (token_index_ + 1 < tokens_.size()) {
        token_index_++;
        current_token_ = tokens_[token_index_];
    } else {
        current_token_ = Token{TokenType::EOF_TOKEN, nullptr, 0, 
                              current_token_.line, current_token_.column + 1};
    }
}

template<typename T>
void Parser<T>::error(const std::string& message) {
    std::string formatted_message = "Parse error at line " + 
        std::to_string(current_token_.line) + ", column " + 
        std::to_string(current_token_.column) + ": " + message;
    throw ParseError(formatted_message, current_token_.line, current_token_.column);
}

// // 显式模板实例化
// template class Parser<int>;
// template class Parser<long long>;
