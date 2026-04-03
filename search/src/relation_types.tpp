int countProps(const std::vector<int>& nu) {
    int props = 0;
    for (int x : nu) {
        if (x != 0) {
            ++props;
        }
    }
    return props;
}

int countDots(const std::vector<int>& nu) {
    int sum = 0;
    for (int x : nu) {
        sum += x;
    }
    return sum - countProps(nu);
}

bool equalNu(const std::vector<int>& lhs, const std::vector<int>& rhs) {
    if (lhs.size() != rhs.size()) {
        return false;
    }
    for (size_t i = 0; i < lhs.size(); ++i) {
        if (lhs[i] != rhs[i]) {
            return false;
        }
    }
    return true;
}

std::string nuToString(const std::vector<int>& nu) {
    std::string s = "{";
    for (size_t i = 0; i < nu.size(); ++i) {
        s += std::to_string(nu[i]);
        if (i + 1 < nu.size()) {
            s += ",";
        }
    }
    s += "}";
    return s;
}

std::string relationVariableToString(const RelationVariable& var) {
    return "c^" + nuToString(var.integral.nu) + "_" + std::to_string(var.k);
}

std::string fiVariableToString(const IntegralLabel& label) {
    return "FI^" + nuToString(label.nu);
}

bool RelationVariableMoreComplexFirst::operator()(
    const RelationVariable& lhs,
    const RelationVariable& rhs) const {
    const int lhsProps = countProps(lhs.integral.nu);
    const int rhsProps = countProps(rhs.integral.nu);
    if (lhsProps != rhsProps) {
        return lhsProps > rhsProps;
    }

    const int lhsDots = countDots(lhs.integral.nu);
    const int rhsDots = countDots(rhs.integral.nu);
    if (lhsDots != rhsDots) {
        return lhsDots > rhsDots;
    }

    if (lhs.k != rhs.k) {
        return lhs.k > rhs.k;
    }

    const size_t n = std::min(lhs.integral.nu.size(), rhs.integral.nu.size());
    for (size_t i = 0; i < n; ++i) {
        if (lhs.integral.nu[i] != rhs.integral.nu[i]) {
            return lhs.integral.nu[i] > rhs.integral.nu[i];
        }
    }
    return lhs.integral.nu.size() > rhs.integral.nu.size();
}

bool IntegralLabelMoreComplexFirst::operator()(
    const IntegralLabel& lhs,
    const IntegralLabel& rhs) const {
    const int lhsProps = countProps(lhs.nu);
    const int rhsProps = countProps(rhs.nu);
    if (lhsProps != rhsProps) {
        return lhsProps > rhsProps;
    }

    const int lhsDots = countDots(lhs.nu);
    const int rhsDots = countDots(rhs.nu);
    if (lhsDots != rhsDots) {
        return lhsDots > rhsDots;
    }

    const size_t n = std::min(lhs.nu.size(), rhs.nu.size());
    for (size_t i = 0; i < n; ++i) {
        if (lhs.nu[i] != rhs.nu[i]) {
            return lhs.nu[i] > rhs.nu[i];
        }
    }
    return lhs.nu.size() > rhs.nu.size();
}
