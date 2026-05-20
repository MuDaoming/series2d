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

bool equalBoundaryTag(const BoundaryTag& lhs, const BoundaryTag& rhs) {
    return lhs.axis == rhs.axis && lhs.side == rhs.side;
}

bool equalIntegralLabel(const IntegralLabel& lhs, const IntegralLabel& rhs) {
    if (lhs.head != rhs.head || !equalNu(lhs.nu, rhs.nu) ||
        lhs.boundaries.size() != rhs.boundaries.size()) {
        return false;
    }
    for (size_t i = 0; i < lhs.boundaries.size(); ++i) {
        if (!equalBoundaryTag(lhs.boundaries[i], rhs.boundaries[i])) {
            return false;
        }
    }
    return true;
}

static int headComplexity(IntegralHead head) {
    switch (head) {
        case IntegralHead::FI: return 2;
        case IntegralHead::BFI: return 1;
        case IntegralHead::BBFI: return 0;
    }
    throw std::runtime_error("Unknown IntegralHead");
}

static int headLexRank(IntegralHead head) {
    switch (head) {
        case IntegralHead::FI: return 0;
        case IntegralHead::BFI: return 1;
        case IntegralHead::BBFI: return 2;
    }
    throw std::runtime_error("Unknown IntegralHead");
}

static int boundaryRank(const BoundaryTag& tag) {
    int axis = tag.axis == BoundaryAxis::X ? 0 : 1;
    int side = tag.side == BoundarySide::U ? 0 : 1;
    return 2 * axis + side;
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

std::string integralHeadToString(IntegralHead head) {
    switch (head) {
        case IntegralHead::FI: return "FI";
        case IntegralHead::BFI: return "BFI";
        case IntegralHead::BBFI: return "BBFI";
    }
    throw std::runtime_error("Unknown IntegralHead");
}

std::string boundaryTagToString(const BoundaryTag& tag) {
    std::string s;
    s += (tag.axis == BoundaryAxis::X) ? "X" : "Y";
    s += (tag.side == BoundarySide::U) ? "U" : "D";
    return s;
}

std::string integralLabelToString(const IntegralLabel& label) {
    std::string s = integralHeadToString(label.head);
    if (!label.boundaries.empty()) {
        s += "[";
        for (size_t i = 0; i < label.boundaries.size(); ++i) {
            s += boundaryTagToString(label.boundaries[i]);
            if (i + 1 < label.boundaries.size()) {
                s += ",";
            }
        }
        s += "]";
    }
    s += nuToString(label.nu);
    return s;
}

std::string relationVariableToString(const RelationVariable& var) {
    return "c^" + integralLabelToString(var.integral) + "_" + std::to_string(var.k);
}

std::string integralVariableToString(const IntegralLabel& label) {
    return integralLabelToString(label);
}

bool IntegralLabelLess::operator()(const IntegralLabel& lhs,
                                   const IntegralLabel& rhs) const {
    if (lhs.head != rhs.head) {
        return headLexRank(lhs.head) < headLexRank(rhs.head);
    }

    if (lhs.boundaries.size() != rhs.boundaries.size()) {
        return lhs.boundaries.size() < rhs.boundaries.size();
    }
    for (size_t i = 0; i < lhs.boundaries.size(); ++i) {
        const int lr = boundaryRank(lhs.boundaries[i]);
        const int rr = boundaryRank(rhs.boundaries[i]);
        if (lr != rr) {
            return lr < rr;
        }
    }

    return std::lexicographical_compare(
        lhs.nu.begin(), lhs.nu.end(), rhs.nu.begin(), rhs.nu.end());
}

bool RelationVariableMoreComplexFirst::operator()(
    const RelationVariable& lhs,
    const RelationVariable& rhs) const {
    IntegralLabelMoreComplexFirst labelCmp;
    if (labelCmp(lhs.integral, rhs.integral)) {
        return true;
    }
    if (labelCmp(rhs.integral, lhs.integral)) {
        return false;
    }
    return lhs.k > rhs.k;
}

bool IntegralLabelMoreComplexFirst::operator()(
    const IntegralLabel& lhs,
    const IntegralLabel& rhs) const {
    const int lhsHead = headComplexity(lhs.head);
    const int rhsHead = headComplexity(rhs.head);
    if (lhsHead != rhsHead) {
        return lhsHead > rhsHead;
    }

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
    if (lhs.nu.size() != rhs.nu.size()) {
        return lhs.nu.size() > rhs.nu.size();
    }

    const size_t nb = std::min(lhs.boundaries.size(), rhs.boundaries.size());
    for (size_t i = 0; i < nb; ++i) {
        const int lr = boundaryRank(lhs.boundaries[i]);
        const int rr = boundaryRank(rhs.boundaries[i]);
        if (lr != rr) {
            return lr > rr;
        }
    }
    return lhs.boundaries.size() > rhs.boundaries.size();
}
