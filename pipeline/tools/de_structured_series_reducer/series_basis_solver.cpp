#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using u128 = __uint128_t;

struct BasisInput {
    u64 p = 0;
    int order = -1;
    std::vector<u64> target;
    std::vector<std::string> labels;
    std::vector<std::vector<u64>> columns;
};

static u64 addm(u64 a, u64 b, u64 p) {
    const u64 c = a + b;
    return c >= p ? c - p : c;
}

static u64 subm(u64 a, u64 b, u64 p) {
    return a >= b ? a - b : p - (b - a);
}

static u64 mulm(u64 a, u64 b, u64 p) {
    return static_cast<u64>((static_cast<u128>(a) * b) % p);
}

static u64 powm(u64 a, u64 n, u64 p) {
    u64 out = 1;
    while (n) {
        if (n & 1) out = mulm(out, a, p);
        n >>= 1;
        if (n) a = mulm(a, a, p);
    }
    return out;
}

static BasisInput read_input(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open input: " + path);
    BasisInput data;
    std::string tag;
    std::size_t count = 0;
    in >> tag >> data.p;
    if (tag != "P" || data.p < 2) throw std::runtime_error("bad P header");
    in >> tag >> data.order;
    if (tag != "ORDER" || data.order < 0) throw std::runtime_error("bad ORDER header");
    in >> tag >> count;
    if (tag != "TARGET" || count != static_cast<std::size_t>(data.order + 1)) {
        throw std::runtime_error("bad TARGET header");
    }
    data.target.resize(count);
    for (u64& value : data.target) in >> value;
    std::size_t column_count = 0;
    in >> tag >> column_count;
    if (tag != "COLUMNS") throw std::runtime_error("bad COLUMNS header");
    data.labels.reserve(column_count);
    data.columns.reserve(column_count);
    for (std::size_t j = 0; j < column_count; ++j) {
        std::string label;
        std::size_t length = 0;
        in >> tag >> label >> length;
        if (tag != "COLUMN" || length != count) {
            throw std::runtime_error("bad COLUMN record at index " + std::to_string(j));
        }
        std::vector<u64> column(length);
        for (u64& value : column) in >> value;
        data.labels.push_back(std::move(label));
        data.columns.push_back(std::move(column));
    }
    if (!in) throw std::runtime_error("truncated input");
    return data;
}

struct SolveResult {
    bool consistent = false;
    bool unique = false;
    int rank = 0;
    int equations_used = 0;
    std::vector<u64> solution;
};

static SolveResult solve_prefix(const BasisInput& data, int equation_count) {
    const int n = static_cast<int>(data.columns.size());
    const int m = std::min(equation_count, data.order + 1);
    std::vector<std::vector<u64>> a(
        m, std::vector<u64>(static_cast<std::size_t>(n) + 1, 0));
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) a[i][j] = data.columns[j][i];
        a[i][n] = data.target[i];
    }

    std::vector<int> pivot_columns;
    int row = 0;
    for (int col = 0; col < n && row < m; ++col) {
        int pivot = row;
        while (pivot < m && a[pivot][col] == 0) ++pivot;
        if (pivot == m) continue;
        std::swap(a[row], a[pivot]);
        const u64 inverse = powm(a[row][col], data.p - 2, data.p);
        for (int j = col; j <= n; ++j) a[row][j] = mulm(a[row][j], inverse, data.p);
        for (int i = 0; i < m; ++i) {
            if (i == row || a[i][col] == 0) continue;
            const u64 factor = a[i][col];
            for (int j = col; j <= n; ++j) {
                a[i][j] = subm(a[i][j], mulm(factor, a[row][j], data.p), data.p);
            }
        }
        pivot_columns.push_back(col);
        ++row;
    }

    bool consistent = true;
    for (int i = row; i < m; ++i) {
        bool zero = true;
        for (int j = 0; j < n; ++j) zero = zero && (a[i][j] == 0);
        if (zero && a[i][n] != 0) {
            consistent = false;
            break;
        }
    }

    SolveResult out;
    out.consistent = consistent;
    out.unique = consistent && row == n;
    out.rank = row;
    out.equations_used = m;
    out.solution.assign(n, 0);
    if (consistent) {
        for (int i = 0; i < row; ++i) out.solution[pivot_columns[i]] = a[i][n];
    }
    return out;
}

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "usage: series_basis_solver INPUT SOLUTION TRAIN_EQUATIONS\n";
        return 2;
    }
    try {
        const BasisInput data = read_input(argv[1]);
        const int requested = std::stoi(argv[3]);
        SolveResult result = solve_prefix(data, requested);
        std::cout << "columns=" << data.columns.size()
                  << " equations_used=" << result.equations_used
                  << " rank=" << result.rank
                  << " consistent=" << result.consistent
                  << " unique=" << result.unique << "\n";
        if (!result.consistent) return 3;
        if (!result.unique) return 4;

        int mismatch_count = 0;
        int first_mismatch = -1;
        for (int i = 0; i <= data.order; ++i) {
            u64 value = 0;
            for (std::size_t j = 0; j < data.columns.size(); ++j) {
                value = addm(value, mulm(result.solution[j], data.columns[j][i], data.p), data.p);
            }
            if (value != data.target[i]) {
                ++mismatch_count;
                if (first_mismatch < 0) first_mismatch = i;
            }
        }
        std::cout << "verification_order=" << data.order
                  << " mismatch_count=" << mismatch_count
                  << " first_mismatch=" << first_mismatch << "\n";

        std::ofstream out(argv[2]);
        if (!out) throw std::runtime_error("cannot open solution output");
        out << "P " << data.p << "\n";
        out << "UNKNOWN_COUNT " << result.solution.size() << "\n";
        out << "TRAIN_EQUATIONS " << result.equations_used << "\n";
        out << "VERIFICATION_ORDER " << data.order << "\n";
        out << "MISMATCH_COUNT " << mismatch_count << "\n";
        for (std::size_t j = 0; j < result.solution.size(); ++j) {
            out << "COEFFICIENT " << data.labels[j] << " " << result.solution[j] << "\n";
        }
        return mismatch_count == 0 ? 0 : 5;
    } catch (const std::exception& error) {
        std::cerr << "series_basis_solver: " << error.what() << "\n";
        return 1;
    }
}
