#pragma once

#include <flint/nmod.h>
#include <flint/fmpz.h>
#include <iostream>
#include <string>
#include <stdexcept>
#include <type_traits>

/// Flint有限域包装类
class FlintMod {
private:
    static nmod_t mod_ctx;
    static bool ctx_initialized;
    mp_limb_t value;
    
public:
    // 构造函数
    FlintMod();
    FlintMod(const FlintMod& other);
    explicit FlintMod(const std::string& val);
    
    // 模板构造函数支持所有整数类型
    template<typename IntType>
    FlintMod(IntType val, std::enable_if_t<std::is_integral_v<IntType>>* = nullptr);
    
    // 设置模数（静态方法）
    static void set_modulus(mp_limb_t p);
    
    // 基本运算
    FlintMod operator+(const FlintMod& other) const;
    FlintMod operator-(const FlintMod& other) const;
    FlintMod operator*(const FlintMod& other) const;
    FlintMod operator/(const FlintMod& other) const;
    FlintMod operator-() const;
    
    FlintMod& operator+=(const FlintMod& other);
    FlintMod& operator-=(const FlintMod& other);
    FlintMod& operator*=(const FlintMod& other);
    FlintMod& operator/=(const FlintMod& other);
    FlintMod& operator=(const FlintMod& other);
    
    // 比较运算
    bool operator==(const FlintMod& other) const;
    bool operator!=(const FlintMod& other) const;
    bool operator<(const FlintMod& other) const;
    bool operator>(const FlintMod& other) const;
    bool operator<=(const FlintMod& other) const;
    bool operator>=(const FlintMod& other) const;
    
    // 获取原始值
    mp_limb_t get_value() const;
    
    // 获取模数（静态方法）
    static mp_limb_t get_modulus();
    
    // 输出支持
    friend std::ostream& operator<<(std::ostream& os, const FlintMod& fm);
    
    // 转换为字符串
    std::string to_string() const;
};

#include "../src/ff_type.tpp"
