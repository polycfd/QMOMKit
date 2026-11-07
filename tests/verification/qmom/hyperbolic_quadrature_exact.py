"""Generate the references of hyperbolic_quadrature_exact.cpp.

For n = 1 to 10 and three shapes of weighted points, the moments
m_0, ..., m_{2n} rounded to double, and the nodes and weights of the HyQMOM
quadrature of those doubles: the recurrence and the closure in exact rational
arithmetic, then the Gauss quadrature of J_{n+1} with mpmath at 60 digits.
Asserts that every sigma(k, k) exceeds its rounding bound tau_k a
thousandfold. Prints the C++ arrays stored in the test. Requires mpmath.

    python3 hyperbolic_quadrature_exact.py
"""

from fractions import Fraction as F
from math import comb

import mpmath as mp

mp.mp.dps = 60
EPS = F(1, 2**52)


def shapes(n):
    # Sums of more than n + 1 weighted points, with moments rounded to double:
    # a symmetric binomial quadrature, a skewed geometric one with weights
    # 625 (3/4)^j, and two binomial clusters of masses 3/4 and 1/4, shifted
    # off the origin.
    p = 2 * n + 2
    symmetric = ([F(j - n - 1, 2) for j in range(p + 1)], [F(comb(p, j), 2**p) for j in range(p + 1)])
    skewed = ([F(j - n - 1, 2) for j in range(p + 1)], [2500 * F(3, 4)**j / 4 for j in range(p + 1)])
    bimodal = ([1 + c + F(j, 2) - F(n, 4) for c in (-2, 2) for j in range(n + 1)],
               [mass * F(comb(n, j), 2**n) for mass in (F(3, 4), F(1, 4)) for j in range(n + 1)])
    return {"symmetric": symmetric, "skewed": skewed, "bimodal": bimodal}


def reference(m, n):
    # The recurrence of m_0, ..., m_{2n} in exact arithmetic,
    # docs/developer/qmom.md, with its absolute-value bound, then the closure
    # a_n = mean(a_0, ..., a_{n-1}).
    sigma = {(-1, l): F(0) for l in range(2 * n + 2)}
    bound = dict(sigma)
    for l in range(2 * n + 1):
        sigma[(0, l)] = m[l]
        bound[(0, l)] = abs(m[l])
    a, b = [m[1] / m[0]], [F(0)]
    margin = None
    for k in range(1, n + 1):
        for l in range(k, 2 * n - k + 1):
            sigma[(k, l)] = sigma[(k - 1, l + 1)] - a[k - 1] * sigma[(k - 1, l)] - b[k - 1] * sigma[(k - 2, l)]
            bound[(k, l)] = bound[(k - 1, l + 1)] + abs(a[k - 1]) * bound[(k - 1, l)] + abs(b[k - 1]) * bound[(k - 2, l)]
        ratio = sigma[(k, k)] / (100 * EPS * bound[(k, k)])
        margin = ratio if margin is None else min(margin, ratio)
        b.append(sigma[(k, k)] / sigma[(k - 1, k - 1)])
        if k < n:
            a.append(sigma[(k, k + 1)] / sigma[(k, k)] - sigma[(k - 1, k)] / sigma[(k - 1, k - 1)])
    a.append(sum(a) / n)
    # The Gauss quadrature of J_{n+1}, diagonalized at 60 digits.
    q = lambda v: mp.mpf(v.numerator) / v.denominator
    J = mp.matrix(n + 1, n + 1)
    for i in range(n + 1):
        J[i, i] = q(a[i])
    for i in range(1, n + 1):
        J[i - 1, i] = J[i, i - 1] = mp.sqrt(q(b[i]))
    E, Q = mp.eigsy(J)
    order = sorted(range(n + 1), key=lambda i: E[i])
    return [E[i] for i in order], [q(m[0]) * Q[0, i] ** 2 for i in order], margin


def main():
    fmt = lambda values: ", ".join(repr(float(v)) for v in values)
    for name in ["symmetric", "skewed", "bimodal"]:
        print(f"const std::vector<HyperbolicReference> {name}_references{{")
        for n in range(1, 11):
            x, w = shapes(n)[name]
            m = [F(float(sum(wi * xi**k for xi, wi in zip(x, w)))) for k in range(2 * n + 1)]
            nodes, weights, margin = reference(m, n)
            assert margin > 1000, (name, n, float(margin))
            print(f"    {{{n}, {{{fmt(m)}}}, {{{fmt(nodes)}}}, {{{fmt(weights)}}}}},")
        print("};")


if __name__ == "__main__":
    main()
