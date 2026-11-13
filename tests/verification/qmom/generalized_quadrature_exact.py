"""Generate the references of generalized_quadrature_exact.cpp.

GQMOM with the normal continuation (Fox, Laurent & Passalacqua 2023,
Definition 1, eqs. 21 and 22, and eq. 36). For each case, the moments are
rounded to double, and the reference is computed from those doubles: the
recurrence of docs/developer/qmom.md in exact rational arithmetic, the
continuation in exact rational arithmetic, and the Gauss quadrature of the
Jacobi matrix with mpmath at 50 digits. Where a closed form exists, the script
asserts that the reference agrees with it.

The script also asserts that the decisions QMOMKit takes on the moments are not
marginal. Writes generalized_quadrature_exact.inc beside this file; the
tolerances are set in the test. It then prints the numbers that the test
header and README quote: the margin sigma(8, 8) / tau_8 of the shifted
Gaussian at N = 8, which the set leaves out, and the largest change of the
nodes, relative to their largest magnitude, when every moment moves by one ulp
in a random direction, over five draws, for the trimodal case at N = 3, 6
and 9. Runs serially, in
about eight minutes. Requires mpmath.

    python3 generalized_quadrature_exact.py
"""

import random
from fractions import Fraction as F
from math import comb, inf, nextafter
from pathlib import Path

import mpmath as mp

EPS = F(1, 2**52)


def recurrence(m, pairs):
    # a_0, ..., a_{pairs-1} and b_0 = 0, b_1, ..., b_pairs from the moments
    # m_0, ..., m_{2 pairs}, or a_0, ..., a_{pairs-1} and b_1, ..., b_{pairs-1}
    # from m_0, ..., m_{2 pairs - 1}; sigma(k, k) with tau_k, QMOMKit's
    # rounding bound 100 eps times the absolute-value recurrence.
    count = len(m)
    sigma = {(-1, l): F(0) for l in range(count + 1)}
    bound = dict(sigma)
    for l in range(count):
        sigma[(0, l)] = m[l]
        bound[(0, l)] = abs(m[l])
    a, b, diagonal = [m[1] / m[0]], [F(0)], [None]
    for k in range(1, pairs + 1):
        if 2 * k > count - 1:
            break
        for l in range(k, count - k):
            sigma[(k, l)] = sigma[(k - 1, l + 1)] - a[k - 1] * sigma[(k - 1, l)] - b[k - 1] * sigma[(k - 2, l)]
            bound[(k, l)] = bound[(k - 1, l + 1)] + abs(a[k - 1]) * bound[(k - 1, l)] + abs(b[k - 1]) * bound[(k - 2, l)]
        b.append(sigma[(k, k)] / sigma[(k - 1, k - 1)])
        diagonal.append((sigma[(k, k)], 100 * EPS * bound[(k, k)]))
        if k + 1 < count - k:
            a.append(sigma[(k, k + 1)] / sigma[(k, k)] - sigma[(k - 1, k)] / sigma[(k - 1, k - 1)])
    return a, b, diagonal


def gauss(a, b, m0, size):
    # The Gauss quadrature of J_size, diagonal a_0, ..., a_{size-1} and
    # off-diagonal sqrt(b_1), ..., sqrt(b_{size-1}): ascending nodes, weights m0
    # times the squared first eigenvector components.
    q = lambda v: mp.mpf(v.numerator) / v.denominator
    J = mp.matrix(size, size)
    for i in range(size):
        J[i, i] = q(a[i])
    for i in range(1, size):
        J[i - 1, i] = J[i, i - 1] = mp.sqrt(q(b[i]))
    E, Q = mp.eigsy(J)
    order = sorted(range(size), key=lambda i: E[i])
    return [E[i] for i in order], [q(m0) * Q[0, i] ** 2 for i in order]


def generalized(m, n, N, nu, k):
    # Definition 1 continued from the k-node quadrature, k <= N: every a_i,
    # i >= k, the mean of a_0, ..., a_{k-1}, b_k kept, b_i = (i/k)^nu b_k.
    a, b, _ = recurrence(m, n)
    mean = sum(a[:k]) / k
    ac = a[:k] + [mean] * (N - k)
    bc = b[: k + 1] + [F(i, k) ** nu * b[k] for i in range(k + 1, N)]
    return gauss(ac, bc, m[0], N)


def classical(m, k):
    # The classical k-node quadrature of m_0, ..., m_{2k-1}.
    a, b, _ = recurrence(m[: 2 * k], k)
    return gauss(a, b, m[0], k)


def evaluate(case, m):
    mp.mp.dps = 50
    if case["kind"] == "continued":
        return generalized(m, case["n"], case["N"], case["nu"], case["base"])
    return classical(m, case["base"])


def moments_of(points, weights, count):
    return [F(float(sum(w * x**k for x, w in zip(points, weights)))) for k in range(count)]


def exact_moments(m):
    # The moments of the case, asserted exact in double.
    for v in m:
        assert F(float(v)) == v, v
    return list(m)


def gaussian_moments(m0, mu, sd, count):
    # Raw moments of m0 N(mu, sd^2): E[Z^j] = (j - 1)!! for even j.
    def central(j):
        if j % 2:
            return F(0)
        r = F(1)
        for i in range(1, j, 2):
            r *= i
        return r

    return [m0 * sum(comb(k, j) * mu ** (k - j) * sd**j * central(j) for j in range(k + 1)) for k in range(count)]


def trimodal_moments(count):
    # Fig. 5 (section 4.1): weights 1/3, means -3, 0, 4, standard deviations
    # 1/10, 1/20, 1/30, rounded to double.
    exact = [sum(gaussian_moments(F(1, 3), mu, sd, count)[k] for mu, sd in ((-3, F(1, 10)), (0, F(1, 20)), (4, F(1, 30)))) for k in range(count)]
    return [F(float(v)) for v in exact]


def hermite(m0, mu, sd, N):
    # Gauss-Hermite rule of m0 N(mu, sd^2), from the Hermite recurrence of
    # eq. 20: a_i = mu, b_i = i sd^2.
    mp.mp.dps = 50
    return gauss([mu] * N, [F(0)] + [i * sd**2 for i in range(1, N)], m0, N)


def chebyshev(m0, a0, b1, N):
    # Footnote 2: nu = 0 and N = 1 give the Chebyshev polynomials of the second
    # kind on a0 +- 2 sqrt(b1), nodes a0 + 2 sqrt(b1) cos(j pi / (N_G + 1)) and
    # weights m0 2 sin^2(j pi / (N_G + 1)) / (N_G + 1).
    mp.mp.dps = 50
    q = lambda v: mp.mpf(v.numerator) / v.denominator
    t = [mp.pi * j / (N + 1) for j in range(N, 0, -1)]
    return [q(a0) + 2 * mp.sqrt(q(b1)) * mp.cos(x) for x in t], [q(m0) * 2 * mp.sin(x) ** 2 / (N + 1) for x in t]


def cases():
    out = []

    def add(group, name, m, n, N, nu, kind, base, reproduced, check_fixed=True, fixed_refused=False, rmin=1.0e-10, closed=None):
        out.append(dict(group=group, name=name, m=m, n=n, N=N, nu=nu, kind=kind, base=base, reproduced=reproduced, check_fixed=check_fixed,
                        fixed_refused=fixed_refused, rmin=rmin, closed=closed))

    # V1: Gaussian moments and nu = 1 give the Gauss-Hermite rule. Their
    # moments are exact in double and their recurrence coefficients, a_k = mu
    # and b_k = k sd^2, small dyadic numbers, so rounding in the recurrence
    # plays little part: these cases mainly test the continuation and the
    # eigensolve. The shifted Gaussian stops at N = 5: at N = 8, sigma(8, 8)
    # exceeds tau_8 by less than the thousandfold margin of assert_decisions,
    # as main() prints.
    for mu, sd, orders in ((F(0), F(1), (1, 3, 5, 8)), (F(3), F(1, 2), (1, 3, 5))):
        for n in orders:
            m = exact_moments(gaussian_moments(F(2), mu, sd, 2 * n + 1))
            for N in (3 * n, 101):
                add("hermite", f"mu{float(mu):g}_sd{float(sd):g}_n{n}_N{N}", m, n, N, 1, "continued", n, 2 * n,
                    closed=lambda m0=F(2), mu=mu, sd=sd, N=N: hermite(m0, mu, sd, N))
    # V2: N = 1 and nu = 0 give the Chebyshev rule of the second kind.
    m = exact_moments([F(2), F(1), F(1)])
    for N in (3, 10, 101):
        add("chebyshev", f"N{N}", m, 1, N, 0, "continued", 1, 2, closed=lambda N=N: chebyshev(F(2), F(1, 2), F(1, 4), N))
    # V3: the paper's figures. Fig. 1: standardized Gaussian moments, N = 5,
    # N_G = 101, nu = 0 and 2. Fig. 5: the trimodal case, N = 3, 6, 9, N_G = 101,
    # nu = 1, and N = 3 with nu = -1.
    m = exact_moments(gaussian_moments(F(1), F(0), F(1), 11))
    for nu in (0, 2):
        add("figures", f"figure1_nu{nu}", m, 5, 101, nu, "continued", 5, 10)
    for n, nu in ((3, 1), (6, 1), (9, 1), (3, -1)):
        add("figures", f"trimodal_n{n}_nu{nu}", trimodal_moments(2 * n + 1), n, 101, nu, "continued", n, 2 * n)
    # V4: moments on the boundary, of dyadic points. Three points at N = 3 give
    # their classical quadrature in both modes; two points at N = 4 give it in
    # adaptive mode, and fixed mode is not_realizable.
    m = exact_moments(moments_of([F(-1), F(1, 2), F(2)], [F(1, 4), F(1, 2), F(1, 4)], 7))
    add("boundary", "three_points_n3", m, 3, 9, 1, "classical", 3, 6)
    m = exact_moments(moments_of([F(-1, 2), F(3, 2)], [F(3, 4), F(1, 4)], 9))
    add("boundary", "two_points_n4", m, 4, 12, 1, "classical", 2, 8, fixed_refused=True)
    # V5: a point of weight 1e-8 makes the classical 3-node rule fail
    # rmin = 1e-6; the 2-node rule passes and is continued.
    m = moments_of([F(-1), F(1), F(3), F(-2)], [F(1, 2), F(1, 2), F(1, 10**8), F(1, 10**9)], 7)
    add("reduction", "rmin_n3", m, 3, 9, 1, "continued", 2, 4, check_fixed=False, rmin=1.0e-6)
    return out


def assert_decisions(case):
    # The decisions QMOMKit takes on the moments in double, with margins.
    m, n = case["m"], case["n"]
    mp.mp.dps = 50
    if case["group"] in ("hermite", "chebyshev", "figures"):
        _, _, diagonal = recurrence(m, n)
        for k in range(1, n + 1):
            assert diagonal[k][0] > 1000 * diagonal[k][1], (case["name"], k)
        if n > 1:
            x, w = classical(m, n)
            assert min(w) / max(w) > 10 * case["rmin"], case["name"]
            gaps = [x[i] - x[i - 1] for i in range(1, n)]
            spread = (x[-1] - x[0]) if n > 2 else max(abs(x[0]), abs(x[1]))
            assert min(gaps) / spread > 1.0e-9, case["name"]
    if case["group"] == "reduction":
        _, _, diagonal = recurrence(m, n)
        x, w = classical(m, 3)
        assert min(w) / max(w) < case["rmin"] / 10, case["name"]
        x, w = classical(m, 2)
        assert min(w) / max(w) > 10 * case["rmin"], case["name"]
        assert diagonal[2][0] > 1000 * diagonal[2][1], case["name"]


def one_ulp_change(m, n, N, draws=5):
    # The largest change of the nodes of Definition 1 with nu = 1, relative to
    # their largest magnitude, when every moment, a double, moves to its
    # neighbor above or below at random.
    mp.mp.dps = 50
    x, _ = generalized(m, n, N, 1, n)
    scale = max(abs(v) for v in x)
    rng = random.Random(n)
    change = mp.mpf(0)
    for _ in range(draws):
        p = [F(nextafter(float(v), inf if rng.random() < 0.5 else -inf)) for v in m]
        y, _ = generalized(p, n, N, 1, n)
        change = max(change, max(abs(u - v) for u, v in zip(y, x)) / scale)
    return change


CASES = cases()


def main():
    for case in CASES:
        assert_decisions(case)
    fmt = lambda values: ", ".join(repr(float(v)) for v in values)
    lines = ["// Generated by generalized_quadrature_exact.py; do not edit.", ""]
    for group in ("hermite", "chebyshev", "figures", "boundary", "reduction"):
        lines.append(f"const std::vector<GeneralizedReference> {group}_references{{")
        for case in CASES:
            if case["group"] != group:
                continue
            x, w = evaluate(case, case["m"])
            if case["closed"] is not None:
                cx, cw = case["closed"]()
                scale = max(abs(v) for v in cx)
                assert max(abs(u - v) for u, v in zip(x, cx)) < mp.mpf(10) ** -40 * scale, case["name"]
                assert max(abs(u - v) for u, v in zip(w, cw)) < mp.mpf(10) ** -40 * sum(cw), case["name"]
            lines.append(f"    {{\"{case['name']}\", {case['n']}, {case['N']}, {float(case['nu'])!r}, {case['rmin']!r}, {case['base']}, "
                         f"{case['reproduced']}, {str(case['check_fixed']).lower()}, {str(case['fixed_refused']).lower()},")
            lines.append(f"     {{{fmt(case['m'])}}},")
            lines.append(f"     {{{fmt(x)}}},")
            lines.append(f"     {{{fmt(w)}}}}},")
        lines.append("};")
        lines.append("")
    Path(__file__).with_suffix(".inc").write_text("\n".join(lines))

    _, _, diagonal = recurrence(exact_moments(gaussian_moments(F(2), F(3), F(1, 2), 17)), 8)
    print(f"Shifted Gaussian at N = 8, left out: sigma(8, 8) / tau_8 = {float(diagonal[8][0] / diagonal[8][1]):.0f}")
    print("Change of the nodes for a change of one ulp in the moments:")
    for n in (3, 6, 9):
        print(f"  trimodal, N = {n}, N_G = 101: {float(one_ulp_change(trimodal_moments(2 * n + 1), n, 101)):.1e}")


if __name__ == "__main__":
    main()
