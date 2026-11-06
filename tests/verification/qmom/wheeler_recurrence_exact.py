"""Generate the references of wheeler_recurrence_exact.cpp.

For six quadratures whose nodes and weights are short binary fractions, the
moments m_0, ..., m_{2n-1}, which are exact in double, and the coefficients
a_k, b_k and sigma(k, k) of their Wheeler recurrence in exact rational
arithmetic. For each quadrature, prints its name and moments, then a_k, b_k
and sigma(k, k) to 17 significant digits, then its largest |node|.

    python3 wheeler_recurrence_exact.py
"""

from fractions import Fraction as F

quadratures = {
    "two_points": ([F(-1, 2), F(5, 4)], [F(1, 4), F(3, 4)]),
    "three_points": ([F(-5, 4), F(1, 2), F(2)], [F(1, 4), F(3, 8), F(3, 8)]),
    "four_points": ([F(-3, 2), F(-1, 4), F(3, 4), F(5, 2)], [F(1, 8), F(3, 8), F(1, 4), F(1, 4)]),
    "five_points": ([F(-2), F(-3, 4), F(1, 4), F(3, 2), F(13, 4)],
                    [F(1, 16), F(1, 4), F(3, 8), F(3, 16), F(1, 8)]),
    "symmetric_three_points": ([F(-1), F(0), F(1)], [F(1, 4), F(1, 2), F(1, 4)]),
    "three_points_at_number_density_1024": ([F(-5, 4), F(1, 2), F(2)], [F(256), F(384), F(384)]),
}


def recurrence(m, n):
    # The Wheeler recurrence of m_0, ..., m_{2n-1}, docs/developer/qmom.md.
    sigma = {(-1, l): F(0) for l in range(2 * n + 1)}
    for l in range(2 * n):
        sigma[(0, l)] = m[l]
    a, b = [m[1] / m[0]], [F(0)]
    for k in range(1, n):
        for l in range(k, 2 * n - k):
            sigma[(k, l)] = sigma[(k - 1, l + 1)] - a[k - 1] * sigma[(k - 1, l)] - b[k - 1] * sigma[(k - 2, l)]
        a.append(sigma[(k, k + 1)] / sigma[(k, k)] - sigma[(k - 1, k)] / sigma[(k - 1, k - 1)])
        b.append(sigma[(k, k)] / sigma[(k - 1, k - 1)])
    return a, b, [sigma[(k, k)] for k in range(n)]


for name, (x, w) in quadratures.items():
    n = len(x)
    m = [sum(wi * xi**k for xi, wi in zip(x, w)) for k in range(2 * n)]
    assert all(F(float(v)) == v for v in m), name
    a, b, s = recurrence(m, n)
    print(name, [repr(float(v)) for v in m])
    for values in (a, b, s):
        print([f"{float(v):.17g}" for v in values])
    print(repr(float(max(abs(v) for v in x))))
