// 静态成员定义
nmod_t FlintMod::mod_ctx;
bool FlintMod::ctx_initialized = false;

// 构造函数
FlintMod::FlintMod() : value(0) {}

FlintMod::FlintMod(const FlintMod& other) : value(other.value) {}

// 模板构造函数支持所有整数类型
template<typename IntType>
FlintMod::FlintMod(IntType val, std::enable_if_t<std::is_integral_v<IntType>>*) {
    if (!ctx_initialized) {
        throw std::runtime_error("FlintMod context not initialized. Call set_modulus first.");
    }
    
    if constexpr (std::is_same_v<IntType, mp_limb_t>) {
        // mp_limb_t 类型直接赋值
        value = val;
    } else if constexpr (std::is_unsigned_v<IntType>) {
        // 无符号类型转换为 unsigned long long 然后调用 nmod_set_si
        value = nmod_set_si(static_cast<unsigned long long>(val), mod_ctx);
    } else {
        // 有符号类型转换为 long long 然后调用 nmod_set_si
        value = nmod_set_si(static_cast<long long>(val), mod_ctx);
    }
}

// 设置模数（静态方法）
void FlintMod::set_modulus(mp_limb_t p) {
    nmod_init(&mod_ctx, p);
    ctx_initialized = true;
}

// 基本运算
FlintMod FlintMod::operator+(const FlintMod& other) const {
    FlintMod result;
    result.value = nmod_add(value, other.value, mod_ctx);
    return result;
}

FlintMod FlintMod::operator-(const FlintMod& other) const {
    FlintMod result;
    result.value = nmod_sub(value, other.value, mod_ctx);
    return result;
}

FlintMod FlintMod::operator*(const FlintMod& other) const {
    FlintMod result;
    result.value = nmod_mul(value, other.value, mod_ctx);
    return result;
}

FlintMod FlintMod::operator/(const FlintMod& other) const {
    FlintMod result;
    result.value = nmod_div(value, other.value, mod_ctx);
    return result;
}

FlintMod FlintMod::operator-() const {
    FlintMod result;
    result.value = nmod_neg(value, mod_ctx);
    return result;
}

FlintMod& FlintMod::operator+=(const FlintMod& other) {
    value = nmod_add(value, other.value, mod_ctx);
    return *this;
}

FlintMod& FlintMod::operator-=(const FlintMod& other) {
    value = nmod_sub(value, other.value, mod_ctx);
    return *this;
}

FlintMod& FlintMod::operator*=(const FlintMod& other) {
    value = nmod_mul(value, other.value, mod_ctx);
    return *this;
}

FlintMod& FlintMod::operator/=(const FlintMod& other) {
    value = nmod_div(value, other.value, mod_ctx);
    return *this;
}

FlintMod& FlintMod::operator=(const FlintMod& other) {
    value = other.value;
    return *this;
}

// 比较运算
bool FlintMod::operator==(const FlintMod& other) const {
    return value == other.value;
}

bool FlintMod::operator!=(const FlintMod& other) const {
    return value != other.value;
}

bool FlintMod::operator<(const FlintMod& other) const {
    return value < other.value;
}

bool FlintMod::operator>(const FlintMod& other) const {
    return value > other.value;
}

bool FlintMod::operator<=(const FlintMod& other) const {
    return value <= other.value;
}

bool FlintMod::operator>=(const FlintMod& other) const {
    return value >= other.value;
}

// 获取原始值
mp_limb_t FlintMod::get_value() const { 
    return value; 
}

// 获取模数
mp_limb_t FlintMod::get_modulus() {
    if (!ctx_initialized) {
        throw std::runtime_error("FlintMod context not initialized. Call set_modulus first.");
    }
    return mod_ctx.n;
}

// 输出支持
std::ostream& operator<<(std::ostream& os, const FlintMod& fm) {
    os << fm.value;
    return os;
}

// 转换为字符串
std::string FlintMod::to_string() const {
    return std::to_string(value);
}
