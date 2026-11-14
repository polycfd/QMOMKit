"""Generate the references of generalized_quadrature_exact.cpp.

GQMOM with the normal continuation (Fox, Laurent & Passalacqua 2023,
Definition 1, eqs. 21 and 22, and eq. 36) and the gamma continuation
(Definition 2, eqs. 23, 25 and 26). For each case, the moments are rounded to
double, and the reference is computed from those doubles: the recurrence of
docs/developer/qmom.md in exact rational arithmetic, the continuation in exact
rational arithmetic, and the Gauss quadrature of the Jacobi matrix with mpmath
at 50 digits. Where a closed form exists, the script asserts that the reference
agrees with it.

The script also asserts that the decisions QMOMKit takes on the moments are not
marginal. Writes generalized_quadrature_exact.inc beside this file; the
tolerances are set in the test. It then prints the numbers that the test
header and README quote: the margin sigma(8, 8) / tau_8 of the shifted
Gaussian at N = 8, which the set leaves out, and the largest change of the
nodes, relative to their largest magnitude, when every moment moves by one ulp
in a random direction, over five draws, for the trimodal case at N = 3, 6
and 9, with nu = 1, and for the gamma mixture at N = 3, 4 and 5. Runs serially,
in about fourteen minutes. Requires mpmath.

    python3 generalized_quadrature_exact.py
"""

import random
from fractions import Fraction as F
from math import comb, inf, nextafter
from pathlib import Path
from sys import float_info

import mpmath as mp

EPS = F(1, 2**52)
MAX_DOUBLE = float_info.max


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


def gamma_recurrence(m, n, N, lower, k):
    # Definition 2 continued from the k-node quadrature, k <= N, on
    # [lower, inf), as the paper states it: the moments M of the coordinate
    # minus lower, alpha = M_1^2 / (M_2 M_0 - M_1^2) - 1, zeta_1, ..., zeta_{2k}
    # of the recurrence of M (eq. 23), the continued zeta_i, i > 2k (eq. 25),
    # and a_i = zeta_{2i} + zeta_{2i+1}, i >= k, and b_i = zeta_{2i-1} zeta_{2i},
    # i > k (eq. 26). The Jacobi matrix of the coordinate adds lower to every
    # a_i.
    M = [sum(comb(j, i) * (-lower) ** (j - i) * m[i] for i in range(j + 1)) for j in range(len(m))]
    a, b, _ = recurrence(M, n)
    alpha = M[1] ** 2 / (M[2] * M[0] - M[1] ** 2) - 1
    zeta = [None] * (2 * N + 1)
    zeta[1] = a[0]
    for j in range(1, k + 1):
        zeta[2 * j] = b[j] / zeta[2 * j - 1]
        if j < k:
            zeta[2 * j + 1] = a[j] - zeta[2 * j]
    for i in range(k + 1, N + 1):
        zeta[2 * i - 1] = (i + alpha) / (k + alpha) * zeta[2 * k - 1]
        zeta[2 * i] = F(i, k) * zeta[2 * k]
    ac = [v + lower for v in a[:k]] + [zeta[2 * i] + zeta[2 * i + 1] + lower for i in range(k, N)]
    bc = b[: k + 1] + [zeta[2 * i - 1] * zeta[2 * i] for i in range(k + 1, N)]
    return ac, bc


def gamma(m, n, N, lower, k):
    ac, bc = gamma_recurrence(m, n, N, lower, k)
    assert all(v > 0 for v in bc[1:]), "continued b must be positive"
    return gauss(ac, bc, m[0], N)


def classical(m, k):
    # The classical k-node quadrature of m_0, ..., m_{2k-1}.
    a, b, _ = recurrence(m[: 2 * k], k)
    return gauss(a, b, m[0], k)


def lower_of(case):
    return case["bounds"][0] if case["bounds"] is not None else F(0)


def evaluate(case, m):
    mp.mp.dps = 50
    if case["kind"] == "classical":
        return classical(m, case["base"])
    if case["continuation"] == "gamma":
        return gamma(m, case["n"], case["N"], lower_of(case), case["base"])
    return generalized(m, case["n"], case["N"], case["nu"], case["base"])


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


def gamma_moments(m0, alpha, beta, lower, count):
    # Raw moments of m0 times the gamma density of shape alpha + 1 and rate
    # beta, shifted to [lower, inf): with y the unshifted coordinate,
    # E[y^j] = (alpha + 1) ... (alpha + j) / beta^j.
    y = [F(1)]
    for j in range(1, count):
        y.append(y[-1] * (alpha + j) / beta)
    return [m0 * sum(comb(k, j) * lower ** (k - j) * y[j] for j in range(k + 1)) for k in range(count)]


def gamma_mixture_moments(count):
    # 1/3 of the gamma density with alpha = 2, beta = 1 and 2/3 of that with
    # alpha = 9, beta = 1/2, rounded to double.
    exact = [u + v for u, v in zip(gamma_moments(F(1, 3), F(2), F(1), F(0), count), gamma_moments(F(2, 3), F(9), F(1, 2), F(0), count))]
    return [F(float(v)) for v in exact]


def laguerre(m0, alpha, beta, lower, N):
    # Gauss rule of the shifted gamma number density, from the recurrence of
    # the generalized Laguerre polynomials (eq. 24): a_i = (2i + alpha + 1) /
    # beta + lower, b_i = i (i + alpha) / beta^2.
    mp.mp.dps = 50
    return gauss([F(2 * i + 1) / beta + alpha / beta + lower for i in range(N)], [F(0)] + [i * (i + alpha) / beta**2 for i in range(1, N)], m0, N)


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

    def add(group, name, m, n, N, nu, kind, base, reproduced, check_fixed=True, fixed_refused=False, rmin=1.0e-10, closed=None,
            continuation="normal", bounds=None):
        out.append(dict(group=group, name=name, m=m, n=n, N=N, nu=nu, kind=kind, base=base, reproduced=reproduced, check_fixed=check_fixed,
                        fixed_refused=fixed_refused, rmin=rmin, closed=closed, continuation=continuation, bounds=bounds))

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
    # V6: gamma moments with the gamma continuation and the default bounds
    # give the Gauss rule of the generalized Laguerre polynomials. Fig. 2:
    # alpha = 10, beta = 1, N = 2, 4, 6, N_G = 101. Then m0 = 2, alpha = 1/2,
    # beta = 2, N = 1, 3, 5, N_G = 10 and 101. Their moments are exact in
    # double.
    for m0, alpha, beta, orders, counts, prefix in ((F(1), F(10), F(1), (2, 4, 6), (101,), "figure2"),
                                                    (F(2), F(1, 2), F(2), (1, 3, 5), (10, 101), "alpha0.5_beta2")):
        for n in orders:
            m = exact_moments(gamma_moments(m0, alpha, beta, F(0), 2 * n + 1))
            for N in counts:
                add("laguerre", f"{prefix}_n{n}_N{N}", m, n, N, 1, "continued", n, 2 * n, continuation="gamma",
                    closed=lambda m0=m0, alpha=alpha, beta=beta, N=N: laguerre(m0, alpha, beta, F(0), N))
    # V7: the gamma density with alpha = 1 and beta = 1 shifted to [lower, inf),
    # lower = -1 and 4, with the bounds [lower, the largest double].
    for lower in (F(-1), F(4)):
        for n in (1, 3, 5):
            m = exact_moments(gamma_moments(F(1), F(1), F(1), lower, 2 * n + 1))
            add("shifted_laguerre", f"lower{float(lower):g}_n{n}_N101", m, n, 101, 1, "continued", n, 2 * n, continuation="gamma",
                bounds=(lower, F(MAX_DOUBLE)), closed=lambda lower=lower: laguerre(F(1), F(1), F(1), lower, 101))
    # V8: a mixture of two gamma number densities, 1/3 of alpha = 2, beta = 1
    # and 2/3 of alpha = 9, beta = 1/2, moments rounded to double, with the
    # default bounds. Not a gamma density, so alpha from m_0, m_1 and m_2
    # differs from the alpha of the last known zeta pair: an alpha taken from
    # that pair moves the nodes by more than 3e-2 of their largest magnitude in
    # every case.
    for n in (3, 4, 5):
        m = gamma_mixture_moments(2 * n + 1)
        for N in (20, 101):
            add("gamma_mixture", f"mixture_n{n}_N{N}", m, n, N, 1, "continued", n, 2 * n, continuation="gamma")
    # V9: moments on the boundary, of dyadic points in (0, inf), as V4.
    m = exact_moments(moments_of([F(1, 2), F(1), F(3)], [F(1, 4), F(1, 2), F(1, 4)], 7))
    add("gamma_boundary", "three_points_n3", m, 3, 9, 1, "classical", 3, 6, continuation="gamma")
    m = exact_moments(moments_of([F(1, 2), F(5, 2)], [F(3, 4), F(1, 4)], 9))
    add("gamma_boundary", "two_points_n4", m, 4, 12, 1, "classical", 2, 8, fixed_refused=True, continuation="gamma")
    # V10: the moments of points -1, 1 and 2 at N = 3 with the default bounds.
    # Their classical 3-node rule has the node -1, and the continuation from
    # the 2-node rule has negative nodes, so adaptive mode continues the
    # 1-node rule and fixed mode is not_realizable.
    m = exact_moments(moments_of([F(-1), F(1), F(2)], [F(1, 4), F(1, 2), F(1, 4)], 7))
    add("gamma_reduction", "outside_n3", m, 3, 9, 1, "continued", 1, 2, fixed_refused=True, continuation="gamma")
    return out


def assert_decisions(case):
    # The decisions QMOMKit takes on the moments in double, with margins.
    m, n = case["m"], case["n"]
    mp.mp.dps = 50
    if case["group"] in ("hermite", "chebyshev", "figures", "laguerre", "shifted_laguerre", "gamma_mixture"):
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
    if case["group"] == "gamma_reduction":
        # Node count 3: the classical rule, on the boundary, has the node -1.
        _, _, diagonal = recurrence(m, n)
        assert abs(diagonal[3][0]) <= diagonal[3][1], case["name"]
        x, _ = classical(m, 3)
        assert x[0] < -0.5, case["name"]
        # Node count 2: the classical rule is interior and accepted, and its
        # continuation has nodes clearly below zero.
        x, w = classical(m, 2)
        assert min(w) / max(w) > 10 * case["rmin"], case["name"]
        assert diagonal[2][0] > 1000 * diagonal[2][1], case["name"]
        x, _ = gamma(m[:5], 2, case["N"], F(0), 2)
        assert x[0] < -1.0e-6 * max(abs(v) for v in x), case["name"]
        # Node count 1: the continued rule follows in main().
        assert diagonal[1][0] > 1000 * diagonal[1][1], case["name"]


def assert_inside(case, x):
    # The continued nodes lie above the lower bound with a margin far beyond
    # the tolerance 1e-12 max |xi_i| of abscissas_inside_bounds.
    if case["continuation"] == "gamma" and case["kind"] == "continued":
        lower = lower_of(case)
        assert x[0] - mp.mpf(lower.numerator) / lower.denominator > 1.0e-9 * max(abs(v) for v in x), case["name"]
        assert case["bounds"] is None or float(x[-1]) < case["bounds"][1], case["name"]


def one_ulp_change(m, n, continued, draws=5):
    # The largest change of the nodes continued(moments), relative to their
    # largest magnitude, when every moment, a double, moves to its neighbor
    # above or below at random.
    mp.mp.dps = 50
    x, _ = continued(m)
    scale = max(abs(v) for v in x)
    rng = random.Random(n)
    change = mp.mpf(0)
    for _ in range(draws):
        p = [F(nextafter(float(v), inf if rng.random() < 0.5 else -inf)) for v in m]
        y, _ = continued(p)
        change = max(change, max(abs(u - v) for u, v in zip(y, x)) / scale)
    return change


CASES = cases()


def main():
    for case in CASES:
        assert_decisions(case)
    fmt = lambda values: ", ".join(repr(float(v)) for v in values)
    lines = ["// Generated by generalized_quadrature_exact.py; do not edit.", ""]
    for group in ("hermite", "chebyshev", "figures", "boundary", "reduction", "laguerre", "shifted_laguerre", "gamma_mixture", "gamma_boundary",
                  "gamma_reduction"):
        lines.append(f"const std::vector<GeneralizedReference> {group}_references{{")
        for case in CASES:
            if case["group"] != group:
                continue
            x, w = evaluate(case, case["m"])
            assert_inside(case, x)
            if case["closed"] is not None:
                cx, cw = case["closed"]()
                scale = max(abs(v) for v in cx)
                assert max(abs(u - v) for u, v in zip(x, cx)) < mp.mpf(10) ** -40 * scale, case["name"]
                assert max(abs(u - v) for u, v in zip(w, cw)) < mp.mpf(10) ** -40 * sum(cw), case["name"]
            lines.append(f"    {{\"{case['name']}\", {case['n']}, {case['N']}, {float(case['nu'])!r}, {case['rmin']!r}, {case['base']}, "
                         f"{case['reproduced']}, {str(case['check_fixed']).lower()}, {str(case['fixed_refused']).lower()},")
            lines.append(f"     {{{fmt(case['m'])}}},")
            lines.append(f"     {{{fmt(x)}}},")
            if case["continuation"] == "normal":
                lines.append(f"     {{{fmt(w)}}}}},")
            else:
                bounds = "std::nullopt" if case["bounds"] is None else f"qmomkit::AbscissaBounds{{{fmt(case['bounds'])}}}"
                lines.append(f"     {{{fmt(w)}}},")
                lines.append(f"     qmomkit::GeneralizedContinuation::gamma, {bounds}}},")
        lines.append("};")
        lines.append("")
    Path(__file__).with_suffix(".inc").write_text("\n".join(lines))

    _, _, diagonal = recurrence(exact_moments(gaussian_moments(F(2), F(3), F(1, 2), 17)), 8)
    print(f"Shifted Gaussian at N = 8, left out: sigma(8, 8) / tau_8 = {float(diagonal[8][0] / diagonal[8][1]):.0f}")
    print("Change of the nodes for a change of one ulp in the moments:")
    for n in (3, 6, 9):
        change = one_ulp_change(trimodal_moments(2 * n + 1), n, lambda p, n=n: generalized(p, n, 101, 1, n))
        print(f"  trimodal, N = {n}, N_G = 101: {float(change):.1e}")
    for n in (3, 4, 5):
        change = one_ulp_change(gamma_mixture_moments(2 * n + 1), n, lambda p, n=n: gamma(p, n, 101, F(0), n))
        print(f"  gamma mixture, N = {n}, N_G = 101: {float(change):.1e}")


if __name__ == "__main__":
    main()
