"""Generate the references of hyperbolic_characteristic_speeds_exact.cpp.

For the skewed and bimodal shapes of hyperbolic_quadrature_exact.py, which it
imports from beside this file, n = 1 to 10, the moments m_0, ..., m_{2n}
rounded to double and the characteristic speeds of the moment equations closed
by HyQMOM: the roots of the characteristic polynomial of their Jacobian, from
the closure in exact rational arithmetic and its exact derivatives, with mpmath
at 60 digits. Also the moments of the normal cases, whose speeds the test
evaluates from closed-form roots, and the moments and speeds of the boundary
cases. Prints the C++ arrays stored in the test. Requires mpmath.

    python3 hyperbolic_characteristic_speeds_exact.py
"""

from fractions import Fraction as F
from math import comb

import mpmath as mp

from hyperbolic_quadrature_exact import shapes

mp.mp.dps = 60


class D:
    # A dual number over the rationals: a value and its exact gradient in
    # m_0, ..., m_{2n}.
    def __init__(self, v, g):
        self.v, self.g = v, g

    def __add__(s, o):
        return D(s.v + o.v, [x + y for x, y in zip(s.g, o.g)])

    def __sub__(s, o):
        return D(s.v - o.v, [x - y for x, y in zip(s.g, o.g)])

    def __mul__(s, o):
        return D(s.v * o.v, [s.v * y + o.v * x for x, y in zip(s.g, o.g)])

    def __truediv__(s, o):
        return D(s.v / o.v, [(x * o.v - s.v * y) / (o.v * o.v) for x, y in zip(s.g, o.g)])


def closure(m, n):
    # m_{2n+1} of the HyQMOM closure and its exact gradient: the recurrence
    # of m_0, ..., m_{2n}, docs/developer/qmom.md, the closure
    # a_n = mean(a_0, ..., a_{n-1}), then the reverse step of Fox & Laurent
    # 2022, Algorithm 2.
    size = 2 * n + 1
    zero = D(F(0), [F(0)] * size)
    s = {(-1, l): zero for l in range(size + 1)}
    for l in range(size):
        s[(0, l)] = D(m[l], [F(int(j == l)) for j in range(size)])
    a, b = [s[(0, 1)] / s[(0, 0)]], [zero]
    for k in range(1, n + 1):
        for l in range(k, 2 * n - k + 1):
            s[(k, l)] = s[(k - 1, l + 1)] - a[k - 1] * s[(k - 1, l)] - b[k - 1] * s[(k - 2, l)]
        b.append(s[(k, k)] / s[(k - 1, k - 1)])
        if k < n:
            a.append(s[(k, k + 1)] / s[(k, k)] - s[(k - 1, k)] / s[(k - 1, k - 1)])
    total = zero
    for value in a:
        total = total + value
    a.append(total / D(F(n), [F(0)] * size))
    s[(n, n + 1)] = s[(n, n)] * (a[n] + s[(n - 1, n)] / s[(n - 1, n - 1)])
    for k in range(n - 1, -1, -1):
        s[(k, 2 * n - k + 1)] = s[(k + 1, 2 * n - k)] + a[k] * s[(k, 2 * n - k)] + b[k] * s[(k - 1, 2 * n - k)]
    return s[(0, size)], [v.v for v in a], [v.v for v in b]


def q(v):
    return mp.mpf(v.numerator) / v.denominator


def eigenvalues(a, off, size):
    J = mp.matrix(size, size)
    for i in range(size):
        J[i, i] = q(a[i])
    for i in range(1, size):
        J[i - 1, i] = J[i, i - 1] = mp.sqrt(q(off[i]))
    E, Q = mp.eigsy(J)
    return sorted(E[i] for i in range(size)), E, Q


def speeds(m, n):
    # The eigenvalues of the Jacobian of the closed system, rows e_{k+1} for
    # k < 2n and the gradient of m_{2n+1}: the roots of
    # lambda^{2n+1} - sum_j (d m_{2n+1} / d m_j) lambda^j.
    m_next, a, b = closure(m, n)
    coefficients = [mp.mpf(1)] + [-q(c) for c in reversed(m_next.g)]
    roots = mp.polyroots(coefficients, maxsteps=4000, extraprec=800)
    assert all(abs(mp.im(r)) < mp.mpf(10)**-40 for r in roots)
    roots = sorted(mp.re(r) for r in roots)
    scale = max(abs(r) for r in roots)
    # They are the eigenvalues of J_n and K_{n+1} (eqs. 4.17 and 4.18), which
    # interlace (Theorem 4.8): K_{n+1} at even positions, J_n at odd ones.
    beta = b[:n] + [F(2 * n + 1, n) * b[n]]
    j_roots, _, _ = eigenvalues(a, b, n)
    k_roots, _, _ = eigenvalues(a, beta, n + 1)
    assert max(abs(x - y) for x, y in zip(roots[1::2], j_roots)) < mp.mpf(10)**-40 * scale
    assert max(abs(x - y) for x, y in zip(roots[0::2], k_roots)) < mp.mpf(10)**-40 * scale
    assert all(x < y for x, y in zip(roots, roots[1:]))
    # The closure is the moment m_{2n+1} of the HyQMOM quadrature, the Gauss
    # quadrature of J_{n+1}.
    _, E, Q = eigenvalues(a, b, n + 1)
    quadrature_moment = sum(q(m[0]) * Q[0, i]**2 * E[i]**(2 * n + 1) for i in range(n + 1))
    assert abs(quadrature_moment - q(m_next.v)) < mp.mpf(10)**-40 * max(1, abs(quadrature_moment))
    gap = min(y - x for x, y in zip(roots, roots[1:])) / scale
    return roots, gap


def exact(values):
    # Moments that double holds exactly.
    assert all(F(float(v)) == v for v in values)
    return values


def normal_moments(m0, mu, sigma, n):
    # m_k = m0 E[(mu + sigma Z)^k], with E[Z^j] = (j - 1)!! for even j.
    double_factorial = lambda j: 1 if j <= 0 else j * double_factorial(j - 2)
    return exact([m0 * sum(comb(k, j) * mu**(k - j) * sigma**j * double_factorial(j - 1) for j in range(0, k + 1, 2)) for k in range(2 * n + 1)])


def main():
    fmt = lambda values: ", ".join(repr(float(v)) for v in values)
    for name in ["skewed", "bimodal"]:
        print(f"const std::vector<SpeedReference> {name}_references{{")
        for n in range(1, 11):
            x, w = shapes(n)[name]
            m = [F(float(sum(wi * xi**k for xi, wi in zip(x, w)))) for k in range(2 * n + 1)]
            roots, gap = speeds(m, n)
            print(f"    {{{n}, {{{fmt(m)}}}, {{{fmt(roots)}}}}},  // gap {mp.nstr(gap, 3)}")
        print("};")

    print("const std::vector<NormalReference> normal_references{")
    for m0, mu, sigma in [(F(2), F(3, 4), F(3, 2)), (F(1, 2), F(-5, 4), F(1, 2))]:
        for n in range(1, 4):
            print(f"    {{{n}, {float(mu)!r}, {float(sigma)!r}, {{{fmt(normal_moments(m0, mu, sigma, n))}}}}},")
    print("};")

    print("const std::vector<BoundaryReference> boundary_references{")
    for n, x, w in [(2, [F(-5, 4), F(11, 4)], [F(1, 4), F(3, 4)]), (3, [F(-3, 2), F(1, 4), F(11, 4)], [F(1, 4), F(1, 2), F(1, 4)]),
                    (3, [F(-5, 4), F(11, 4)], [F(1, 4), F(3, 4)]), (3, [F(1, 2)], [F(2)])]:
        m = exact([sum(wi * xi**k for xi, wi in zip(x, w)) for k in range(2 * n + 1)])
        print(f"    {{{n}, {{{fmt(m)}}}, {{{fmt(sorted(x + x + [sum(x) / len(x)]))}}}}},")
    print("};")


if __name__ == "__main__":
    main()
