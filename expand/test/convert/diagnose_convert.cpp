#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../include/converter.hpp"
#include "../../include/io.hpp"

#include <ginac/ginac.h>

namespace {

std::string toString(const GiNaC::ex& e) {
    std::ostringstream oss;
    oss << e;
    return oss.str();
}

std::string toString(const GiNaC::numeric& n) {
    std::ostringstream oss;
    oss << n;
    return oss.str();
}

std::string stripSign(std::string s) {
    if (!s.empty() && (s[0] == '-' || s[0] == '+')) {
        s.erase(s.begin());
    }
    return s;
}

bool outsideLongRange(const std::string& signedDecimal) {
    std::string s = stripSign(signedDecimal);
    s.erase(0, s.find_first_not_of('0'));
    if (s.empty()) return false;

    const bool neg = !signedDecimal.empty() && signedDecimal[0] == '-';
    const std::string limit = neg ? "9223372036854775808" : "9223372036854775807";
    if (s.size() != limit.size()) return s.size() > limit.size();
    return s > limit;
}

FlintMod decimalToMod(const std::string& signedDecimal) {
    bool neg = false;
    size_t pos = 0;
    if (pos < signedDecimal.size() && (signedDecimal[pos] == '-' || signedDecimal[pos] == '+')) {
        neg = signedDecimal[pos] == '-';
        ++pos;
    }

    const mp_limb_t p = FlintMod::get_modulus();
    unsigned __int128 rem = 0;
    for (; pos < signedDecimal.size(); ++pos) {
        unsigned char ch = static_cast<unsigned char>(signedDecimal[pos]);
        if (!std::isdigit(ch)) {
            throw std::runtime_error("Non-decimal integer: " + signedDecimal);
        }
        rem = (rem * 10 + (signedDecimal[pos] - '0')) % p;
    }

    mp_limb_t val = static_cast<mp_limb_t>(rem);
    if (neg && val != 0) {
        val = p - val;
    }
    return FlintMod(val);
}

FlintMod numericToModSafe(const GiNaC::numeric& num) {
    if (num.is_integer()) {
        return decimalToMod(toString(num));
    }
    if (num.is_rational()) {
        return decimalToMod(toString(num.numer())) / decimalToMod(toString(num.denom()));
    }
    throw std::runtime_error("Coefficient is not rational: " + toString(num));
}

Polynomial<FlintMod> exToPolynomialSafe(const GiNaC::ex& expr,
                                         const GiNaC::symbol& X,
                                         const GiNaC::symbol& Y) {
    using namespace GiNaC;

    Polynomial<FlintMod> result;
    ex expanded = expand(expr);
    int maxDegX = expanded.degree(X);
    int maxDegY = expanded.degree(Y);

    for (int i = 0; i <= maxDegX; ++i) {
        for (int j = 0; j <= maxDegY; ++j) {
            ex coeff = expanded.coeff(X, i).coeff(Y, j);
            if (coeff.is_zero()) continue;
            if (!is_a<numeric>(coeff)) {
                throw std::runtime_error("Non-numeric coefficient at X^" + std::to_string(i) +
                                         " Y^" + std::to_string(j) + ": " + toString(coeff));
            }
            result.addMonomial(numericToModSafe(ex_to<numeric>(coeff)), Power(i, j));
        }
    }

    return result;
}

Rational<FlintMod> exToRationalSafe(const GiNaC::ex& expr,
                                     const GiNaC::symbol& X,
                                     const GiNaC::symbol& Y) {
    using namespace GiNaC;
    ex normalized = normal(expr);
    return Rational<FlintMod>(
        exToPolynomialSafe(numer(normalized), X, Y),
        exToPolynomialSafe(denom(normalized), X, Y));
}

struct PolyStats {
    int terms = 0;
    int outsideLong = 0;
    std::vector<std::string> examples;
};

void recordNumericIssue(PolyStats& stats,
                        const std::string& label,
                        const GiNaC::numeric& num) {
    auto checkOne = [&](const std::string& part, const GiNaC::numeric& n) {
        std::string s = toString(n);
        if (!outsideLongRange(s)) return;
        ++stats.outsideLong;
        if (stats.examples.size() < 5) {
            std::ostringstream oss;
            oss << label << " " << part << "=" << s;
            try {
                oss << " to_long=" << n.to_long();
            } catch (const std::exception& e) {
                oss << " to_long threw: " << e.what();
            }
            stats.examples.push_back(oss.str());
        }
    };

    if (num.is_integer()) {
        checkOne("coeff", num);
    } else if (num.is_rational()) {
        checkOne("numer", num.numer());
        checkOne("denom", num.denom());
    }
}

PolyStats analyzePoly(const GiNaC::ex& expr,
                      const GiNaC::symbol& X,
                      const GiNaC::symbol& Y) {
    using namespace GiNaC;

    PolyStats stats;
    ex expanded = expand(expr);
    int maxDegX = expanded.degree(X);
    int maxDegY = expanded.degree(Y);

    for (int i = 0; i <= maxDegX; ++i) {
        for (int j = 0; j <= maxDegY; ++j) {
            ex coeff = expanded.coeff(X, i).coeff(Y, j);
            if (coeff.is_zero()) continue;
            ++stats.terms;
            if (!is_a<numeric>(coeff)) {
                if (stats.examples.size() < 5) {
                    stats.examples.push_back("non_numeric X^" + std::to_string(i) +
                                             " Y^" + std::to_string(j) + ": " + toString(coeff));
                }
                continue;
            }
            recordNumericIssue(stats,
                               "X^" + std::to_string(i) + " Y^" + std::to_string(j),
                               ex_to<numeric>(coeff));
        }
    }

    return stats;
}

struct RatStats {
    PolyStats numerator;
    PolyStats denominator;
};

RatStats analyzeRat(const GiNaC::ex& expr,
                    const GiNaC::symbol& X,
                    const GiNaC::symbol& Y) {
    using namespace GiNaC;
    ex normalized = normal(expr);
    return {analyzePoly(numer(normalized), X, Y),
            analyzePoly(denom(normalized), X, Y)};
}

FlintMod evalRat(const Rational<FlintMod>& r, const FlintMod& x, const FlintMod& y) {
    FlintMod den = r.denominator.evaluate(x, y);
    if (den == FlintMod(0)) {
        throw std::runtime_error("zero denominator at test point");
    }
    return r.numerator.evaluate(x, y) / den;
}

bool sameAtSamplePoints(const Rational<FlintMod>& a,
                        const Rational<FlintMod>& b,
                        std::string& detail) {
    const std::vector<std::pair<long, long>> points = {
        {0, 0}, {1, 0}, {0, 1}, {1, 1}, {2, 3}, {5, 7}, {11, 13}
    };

    int checked = 0;
    for (const auto& [px, py] : points) {
        try {
            FlintMod av = evalRat(a, FlintMod(px), FlintMod(py));
            FlintMod bv = evalRat(b, FlintMod(px), FlintMod(py));
            ++checked;
            if (av != bv) {
                std::ostringstream oss;
                oss << "DIFF at (" << px << "," << py << "): current=" << av
                    << " safe=" << bv;
                detail = oss.str();
                return false;
            }
        } catch (const std::exception&) {
            continue;
        }
    }

    detail = "OK on " + std::to_string(checked) + " sample points";
    return true;
}

void reportExpr(std::ostream& out,
                const std::string& name,
                const GiNaC::ex& expr,
                const GiNaC::symbol& X,
                const GiNaC::symbol& Y) {
    out << "\n[" << name << "]\n";
    RatStats stats = analyzeRat(expr, X, Y);
    out << "numerator_terms = " << stats.numerator.terms
        << ", numerator_outside_long = " << stats.numerator.outsideLong << "\n";
    for (const auto& e : stats.numerator.examples) {
        out << "  example: " << e << "\n";
    }
    out << "denominator_terms = " << stats.denominator.terms
        << ", denominator_outside_long = " << stats.denominator.outsideLong << "\n";
    for (const auto& e : stats.denominator.examples) {
        out << "  example: " << e << "\n";
    }

    Rational<FlintMod> current = exToRational(expr, X, Y);
    Rational<FlintMod> safe = exToRationalSafe(expr, X, Y);
    std::string detail;
    bool same = sameAtSamplePoints(current, safe, detail);
    out << "source_exToRational_vs_safe = " << (same ? "OK" : "DIFF")
        << " (" << detail << ")\n";
}

} // namespace

int main(int argc, char** argv) {
    const std::string workDir = argc > 1 ? argv[1] : "expand/test/convert";
    const std::string inputDir = argc > 2 ? argv[2] : "expand/test/2dseries/dp_planar";

    try {
        InputConfig cfg = parseConfigFile(inputDir + "/config");
        FlintMod::set_modulus(cfg.p);

        GiNaC::symbol X("X"), Y("Y");
        auto topS = parseMatrixFile(inputDir + "/S", X, Y);
        Family<GiNaC::ex, GiNaC::ex, GiNaC::ex> ginacFamily(topS, cfg.N, cfg.B, X, Y);

        std::vector<int> topNu(cfg.N, 1);
        const auto* sector = ginacFamily.getSector(topNu);
        if (!sector) {
            throw std::runtime_error("topsector was not found");
        }

        std::ofstream wl(workDir + "/ginac_topsector_C_z.wl");
        if (!wl) {
            throw std::runtime_error("Cannot open ginac_topsector_C_z.wl");
        }
        wl << "ginacC = " << sector->getCSum() << ";\n";
        wl << "ginacCandZ = {";
        for (int i = 0; i < cfg.B; ++i) {
            if (i) wl << ", ";
            wl << sector->getC(i);
        }
        for (int i = 0; i < cfg.N; ++i) {
            wl << ", " << sector->getZ(i);
        }
        wl << "};\n";

        std::ofstream out(workDir + "/diagnose_convert_report.txt");
        if (!out) {
            throw std::runtime_error("Cannot open diagnose_convert_report.txt");
        }

        out << "input_dir = " << inputDir << "\n";
        out << "p = " << cfg.p << "\n";
        out << "N = " << cfg.N << ", B = " << cfg.B << "\n";
        out << "The comparison below calls the existing source exToRational and compares it\n";
        out << "against a local diagnostic converter that reduces decimal GiNaC integers\n";
        out << "mod p without going through long/to_long().\n";

        for (int i = 0; i < cfg.B; ++i) {
            reportExpr(out, "c" + std::to_string(i + 1), sector->getC(i), X, Y);
        }
        reportExpr(out, "C_sum", sector->getCSum(), X, Y);
        for (int i = 0; i < cfg.N; ++i) {
            reportExpr(out, "z" + std::to_string(i + 1), sector->getZ(i), X, Y);
        }

        std::cout << "Wrote " << workDir << "/diagnose_convert_report.txt\n";
        std::cout << "Wrote " << workDir << "/ginac_topsector_C_z.wl\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
