#include <flint/fmpz.h>
#include <flint/nmod.h>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using ModValue = mp_limb_t;

ModValue fromFmpzString(const std::string& s, ModValue p) {
    fmpz_t z;
    fmpz_init(z);
    if (fmpz_set_str(z, s.c_str(), 10) != 0) {
        fmpz_clear(z);
        throw std::runtime_error("fmpz_set_str failed for: " + s);
    }
    ModValue value = fmpz_fdiv_ui(z, p);
    fmpz_clear(z);
    return value;
}

ModValue fromDecimalLoop(const std::string& s, ModValue p) {
    bool neg = false;
    size_t i = 0;
    if (i < s.size() && (s[i] == '-' || s[i] == '+')) {
        neg = s[i] == '-';
        ++i;
    }

    unsigned __int128 rem = 0;
    for (; i < s.size(); ++i) {
        if (s[i] < '0' || s[i] > '9') {
            throw std::runtime_error("bad digit in: " + s);
        }
        rem = (rem * 10 + static_cast<unsigned>(s[i] - '0')) % p;
    }

    ModValue value = static_cast<ModValue>(rem);
    if (neg && value != 0) value = p - value;
    return value;
}

} // namespace

int main() {
    const ModValue p = 2305843009213693951ULL;
    nmod_t ctx;
    nmod_init(&ctx, p);

    const std::vector<std::string> tests = {
        "0",
        "1",
        "-1",
        "9223372036854775807",
        "9223372036854775808",
        "11601456425899298122",
        "-22037830804544671200",
        "1570192927044706665000",
        "-12647814950176488287696",
        "999999999999999999999999999999999999999999999999999"
    };

    for (const auto& s : tests) {
        ModValue viaFmpz = fromFmpzString(s, p);
        ModValue viaLoop = fromDecimalLoop(s, p);
        std::cout << s << "\n";
        std::cout << "  fmpz_fdiv_ui = " << viaFmpz << "\n";
        std::cout << "  decimal_loop = " << viaLoop << "\n";
        std::cout << "  equal = " << (viaFmpz == viaLoop ? "yes" : "no") << "\n";
    }

    fmpz_t numer;
    fmpz_t denom;
    fmpz_init(numer);
    fmpz_init(denom);
    fmpz_set_str(numer, "1570192927044706665000", 10);
    fmpz_set_str(denom, "-12647814950176488287696", 10);
    ModValue n = fmpz_fdiv_ui(numer, p);
    ModValue d = fmpz_fdiv_ui(denom, p);
    ModValue q = nmod_div(n, d, ctx);
    std::cout << "rational_sample\n";
    std::cout << "  numer_mod = " << n << "\n";
    std::cout << "  denom_mod = " << d << "\n";
    std::cout << "  numer/denom mod p = " << q << "\n";
    fmpz_clear(numer);
    fmpz_clear(denom);

    return 0;
}
