"""Generate the references of generalized_quadrature_exact.cpp.

GQMOM with the normal continuation (Fox, Laurent & Passalacqua 2023,
Definition 1, eqs. 21 and 22, and eq. 36), the gamma continuation
(Definition 2, eqs. 23, 25 and 26), the lognormal continuation (Definition 3,
eqs. 23, 26, 29 and 30) and the beta continuation (Definition 4, eqs. 23, 26
and 31 to 35). For each case, the moments are rounded to double, and the
reference is computed from those doubles: the recurrence of
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
    # Definition 1 of Fox, Laurent & Passalacqua (2023), continued from the
    # k-node quadrature, k <= N: every a_i, i >= k, the mean of a_0, ...,
    # a_{k-1}, b_k kept, b_i = (i/k)^nu b_k.
    a, b, _ = recurrence(m, n)
    mean = sum(a[:k]) / k
    ac = a[:k] + [mean] * (N - k)
    bc = b[: k + 1] + [F(i, k) ** nu * b[k] for i in range(k + 1, N)]
    return gauss(ac, bc, m[0], N)


def gamma_recurrence(m, n, N, lower, k):
    # Definition 2 of Fox, Laurent & Passalacqua (2023), continued from the
    # k-node quadrature, k <= N, on [lower, inf), as the paper states it: the
    # moments M of the coordinate minus lower,
    # alpha = M_1^2 / (M_2 M_0 - M_1^2) - 1, zeta_1, ..., zeta_{2k} of the
    # recurrence of M (eq. 23), the continued zeta_i, i > 2k (eq. 25), and
    # a_i = zeta_{2i} + zeta_{2i+1}, i >= k, and b_i = zeta_{2i-1} zeta_{2i},
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


def lognormal_recurrence(m, n, N, lower, k):
    # Definition 3 of Fox, Laurent & Passalacqua (2023), continued from the
    # k-node quadrature, k <= N, on [lower, inf), as the paper states it: the
    # moments M of the coordinate minus lower, eta^2 = M_2 M_0 / M_1^2
    # (eq. 29), zeta_1, ..., zeta_{2k} of the recurrence of M (eq. 23), the
    # continued zeta_i, i > 2k (eq. 30), and a_i = zeta_{2i} + zeta_{2i+1},
    # i >= k, and b_i = zeta_{2i-1} zeta_{2i}, i > k (eq. 26). Eq. 30 has only
    # even powers of eta, so eta^2 keeps the arithmetic exact. The Jacobi matrix
    # of the coordinate adds lower to every a_i.
    M = [sum(comb(j, i) * (-lower) ** (j - i) * m[i] for i in range(j + 1)) for j in range(len(m))]
    a, b, _ = recurrence(M, n)
    eta2 = M[2] * M[0] / M[1] ** 2
    zeta = [None] * (2 * N + 1)
    zeta[1] = a[0]
    for j in range(1, k + 1):
        zeta[2 * j] = b[j] / zeta[2 * j - 1]
        if j < k:
            zeta[2 * j + 1] = a[j] - zeta[2 * j]
    for i in range(k + 1, N + 1):
        zeta[2 * i - 1] = eta2 ** (2 * (i - k)) * zeta[2 * k - 1]
        zeta[2 * i] = eta2 ** (i - k) * (eta2 ** i - 1) / (eta2 ** k - 1) * zeta[2 * k]
    ac = [v + lower for v in a[:k]] + [zeta[2 * i] + zeta[2 * i + 1] + lower for i in range(k, N)]
    bc = b[: k + 1] + [zeta[2 * i - 1] * zeta[2 * i] for i in range(k + 1, N)]
    return ac, bc


def lognormal(m, n, N, lower, k):
    ac, bc = lognormal_recurrence(m, n, N, lower, k)
    assert all(v > 0 for v in bc[1:]), "continued b must be positive"
    return gauss(ac, bc, m[0], N)


def beta_recurrence(m, n, N, lower, upper, k):
    # Definition 4 of Fox, Laurent & Passalacqua (2023), continued from the
    # k-node quadrature, k <= N, on [lower, upper], as the paper states it for
    # (0, 1) after the linear change of variables y = (x - lower) / (upper -
    # lower): the moments M of y, zeta_1, ..., zeta_{2k} of the recurrence of M
    # (eq. 23), the canonical moments p_i with zeta_i = p_i (1 - p_{i-1}) and
    # p_0 = 0 (eq. 31), alpha and beta from p_1 and p_2 (eq. 33), the canonical
    # moments p^J_i of the beta number density (eq. 32), the continued p_i,
    # i > 2k (eqs. 34 and 35), zeta_i = p_i (1 - p_{i-1}) for i > 2k, and
    # a_i = zeta_{2i} + zeta_{2i+1}, i >= k, and b_i = zeta_{2i-1} zeta_{2i},
    # i > k (eq. 26). The Jacobi matrix of x takes lower + (upper - lower) a_i
    # and (upper - lower)^2 b_i. Also returns, for each continued p_i, what
    # decides its case of eqs. 34 and 35, and the margin of that decision. The
    # first case is decided by its first condition, p_l <= p^J_l, with the
    # larger amount by which either condition holds, or by its second
    # condition alone, p^J_l >= p^J_j, with the smaller amount by which the
    # first fails and the second holds. The second case has the smaller amount
    # by which both conditions fail.
    width = upper - lower
    M = [sum(comb(j, i) * (-lower) ** (j - i) * m[i] for i in range(j + 1)) / width**j for j in range(len(m))]
    a, b, _ = recurrence(M, n)
    zeta = [None] * (2 * N + 1)
    zeta[1] = a[0]
    for j in range(1, k + 1):
        zeta[2 * j] = b[j] / zeta[2 * j - 1]
        if j < k:
            zeta[2 * j + 1] = a[j] - zeta[2 * j]
    p = [F(0)] * (2 * N + 1)
    for i in range(1, 2 * k + 1):
        p[i] = zeta[i] / (1 - p[i - 1])
    alpha = (1 - p[1] - 2 * p[2] + p[1] * p[2]) / p[2]
    beta = (p[1] - p[2] - p[1] * p[2]) / p[2]

    def jacobi(j):
        i = (j + 1) // 2
        return (beta + i) / (2 * i + alpha + beta) if j % 2 else F(i) / (2 * i + 1 + alpha + beta)

    choices = []
    for j in range(2 * k + 1, 2 * N + 1):
        last = 2 * k - 1 if j % 2 else 2 * k
        if p[last] <= jacobi(last):
            p[j] = jacobi(j) / jacobi(last) * p[last]
            choices.append((j, "first condition", max(jacobi(last) - p[last], jacobi(last) - jacobi(j))))
        elif jacobi(last) >= jacobi(j):
            p[j] = jacobi(j) / jacobi(last) * p[last]
            choices.append((j, "second condition", min(p[last] - jacobi(last), jacobi(last) - jacobi(j))))
        else:
            p[j] = (p[last] * (1 - jacobi(j)) + jacobi(j) - jacobi(last)) / (1 - jacobi(last))
            choices.append((j, "second case", min(p[last] - jacobi(last), jacobi(j) - jacobi(last))))
        zeta[j] = p[j] * (1 - p[j - 1])
    ay = a[:k] + [zeta[2 * i] + zeta[2 * i + 1] for i in range(k, N)]
    by = b[: k + 1] + [zeta[2 * i - 1] * zeta[2 * i] for i in range(k + 1, N)]
    return [lower + width * v for v in ay], [width**2 * v for v in by], choices


def beta(m, n, N, lower, upper, k):
    ac, bc, _ = beta_recurrence(m, n, N, lower, upper, k)
    assert all(v > 0 for v in bc[1:]), "continued b must be positive"
    return gauss(ac, bc, m[0], N)


def classical(m, k):
    # The classical k-node quadrature of m_0, ..., m_{2k-1}.
    a, b, _ = recurrence(m[: 2 * k], k)
    return gauss(a, b, m[0], k)


def lower_of(case):
    return case["bounds"][0] if case["bounds"] is not None else F(0)


def upper_of(case):
    # The upper bound, [0, 1] for the beta continuation without bounds.
    if case["bounds"] is not None:
        return case["bounds"][1]
    return F(1) if case["continuation"] == "beta" else F(MAX_DOUBLE)


def evaluate(case, m):
    mp.mp.dps = 50
    if case["kind"] == "classical":
        return classical(m, case["base"])
    if case["continuation"] == "gamma":
        return gamma(m, case["n"], case["N"], lower_of(case), case["base"])
    if case["continuation"] == "lognormal":
        return lognormal(m, case["n"], case["N"], lower_of(case), case["base"])
    if case["continuation"] == "beta":
        return beta(m, case["n"], case["N"], lower_of(case), upper_of(case), case["base"])
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
    # Fig. 5 of Fox, Laurent & Passalacqua (2023, section 4.1): weights 1/3,
    # means -3, 0, 4, standard deviations 1/10, 1/20, 1/30, rounded to double.
    exact = [sum(gaussian_moments(F(1, 3), mu, sd, count)[k] for mu, sd in ((-3, F(1, 10)), (0, F(1, 20)), (4, F(1, 30)))) for k in range(count)]
    return [F(float(v)) for v in exact]


def hermite(m0, mu, sd, N):
    # Gauss-Hermite rule of m0 N(mu, sd^2), from the Hermite recurrence
    # (Fox, Laurent & Passalacqua 2023, eq. 20): a_i = mu, b_i = i sd^2.
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
    # the generalized Laguerre polynomials (Fox, Laurent & Passalacqua 2023,
    # eq. 24): a_i = (2i + alpha + 1) / beta + lower, b_i = i (i + alpha) /
    # beta^2.
    mp.mp.dps = 50
    return gauss([F(2 * i + 1) / beta + alpha / beta + lower for i in range(N)], [F(0)] + [i * (i + alpha) / beta**2 for i in range(1, N)], m0, N)


def lognormal_moments(m0, scale, eta, lower, count):
    # Raw moments of m0 times the lognormal density (Fox, Laurent & Passalacqua
    # 2023, eq. 27) with exp(mu) = scale and eta = exp(sigma^2 / 2), shifted to
    # [lower, inf): with y the unshifted coordinate, E[y^j] = scale^j eta^(j^2).
    y = [scale ** j * eta ** (j * j) for j in range(count)]
    return [m0 * sum(comb(k, j) * lower ** (k - j) * y[j] for j in range(k + 1)) for k in range(count)]


def lognormal_mixture_moments(count):
    # 1/2 of the lognormal density with scale 1 and eta = 101/100 and 1/2 of
    # that with scale 5/4 and eta = 51/50, rounded to double.
    exact = [u + v for u, v in zip(lognormal_moments(F(1, 2), F(1), F(101, 100), F(0), count),
                                   lognormal_moments(F(1, 2), F(5, 4), F(51, 50), F(0), count))]
    return [F(float(v)) for v in exact]


def stieltjes_wigert(m0, scale, eta, lower, N):
    # Gauss rule of the shifted lognormal number density, from the zeta of the
    # Stieltjes-Wigert polynomials (Fox, Laurent & Passalacqua 2023, eq. 28):
    # zeta_{2i-1} = scale eta^(4i-3) and
    # zeta_{2i} = scale eta^(2i-1) (eta^(2i) - 1), and eq. 23 with every a_i
    # shifted by lower.
    mp.mp.dps = 50
    zeta = [None]
    for i in range(1, N + 1):
        zeta += [scale * eta ** (4 * i - 3), scale * eta ** (2 * i - 1) * (eta ** (2 * i) - 1)]
    a = [zeta[1] + lower] + [zeta[2 * i] + zeta[2 * i + 1] + lower for i in range(1, N)]
    b = [F(0)] + [zeta[2 * i - 1] * zeta[2 * i] for i in range(1, N)]
    return gauss(a, b, m0, N)


def beta_moments(m0, alpha, beta, lower, upper, count):
    # Raw moments of m0 times the beta number density x^beta (1 - x)^alpha on
    # (0, 1), in the notation of Fox, Laurent & Passalacqua (2023, eq. 32),
    # mapped to [lower, upper]: with y the coordinate on (0, 1),
    # E[y^j] = E[y^(j-1)] (beta + j) / (alpha + beta + 1 + j).
    y = [F(1)]
    for j in range(1, count):
        y.append(y[-1] * (beta + j) / (alpha + beta + 1 + j))
    width = upper - lower
    return [m0 * sum(comb(k, j) * lower ** (k - j) * width**j * y[j] for j in range(k + 1)) for k in range(count)]


def beta_mixture_moments(parts, count):
    # The moments of a mixture of beta number densities on (0, 1), parts of
    # (weight, alpha, beta), rounded to double.
    exact = [sum(beta_moments(w, alpha, beta, F(0), F(1), count)[k] for w, alpha, beta in parts) for k in range(count)]
    return [F(float(v)) for v in exact]


def jacobi(m0, alpha, beta, lower, upper, N):
    # Gauss rule of the beta number density x^beta (1 - x)^alpha mapped to
    # [lower, upper], from the recurrence of the monic Jacobi polynomials
    # P^(alpha, beta) on (-1, 1) (DLMF 18.9.2, Table 18.9.1), with weight
    # (1 - t)^alpha (1 + t)^beta, under t = 2 (x - lower) / (upper - lower) - 1,
    # with b_1 from its own form, which holds for alpha + beta = -1 as well.
    # Independent of eqs. 31 to 35.
    mp.mp.dps = 50
    width = upper - lower
    a, b = [], [F(0)]
    for n in range(N):
        s = 2 * n + alpha + beta
        diagonal = (beta**2 - alpha**2) / (s * (s + 2)) if s != 0 else (beta - alpha) / (alpha + beta + 2)
        a.append(lower + width * (1 + diagonal) / 2)
        if n == 1:
            b.append(width**2 * (1 + alpha) * (1 + beta) / ((2 + alpha + beta) ** 2 * (3 + alpha + beta)))
        elif n > 1:
            b.append(width**2 * n * (n + alpha) * (n + beta) * (n + alpha + beta) / (s**2 * (s + 1) * (s - 1)))
    return gauss(a, b, m0, N)


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
    # V3: the figures of Fox, Laurent & Passalacqua (2023). Fig. 1: standardized
    # Gaussian moments, N = 5, N_G = 101, nu = 0 and 2. Fig. 5: the trimodal
    # case, N = 3, 6, 9, N_G = 101, nu = 1, and N = 3 with nu = -1.
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
    # give the Gauss rule of the generalized Laguerre polynomials. Fig. 2 of
    # Fox, Laurent & Passalacqua (2023): alpha = 10, beta = 1, N = 2, 4, 6,
    # N_G = 101. Then m0 = 2, alpha = 1/2, beta = 2, N = 1, 3, 5, N_G = 10 and
    # 101. Their moments are exact in double.
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
    # differs from the alpha of the last known zeta pair.
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
    # V11: lognormal moments with the lognormal continuation and the default
    # bounds give the Gauss rule of the Stieltjes-Wigert polynomials: eta = 9/8
    # with m0 = 1 and scale 1, and eta = 5/4 with m0 = 2 and scale 2, N = 1, 2,
    # N_G = 10 and 20. Their moments are exact in double.
    for m0, scale, eta, prefix in ((F(1), F(1), F(9, 8), "eta1.125"), (F(2), F(2), F(5, 4), "eta1.25")):
        for n in (1, 2):
            m = exact_moments(lognormal_moments(m0, scale, eta, F(0), 2 * n + 1))
            for N in (10, 20):
                add("stieltjes_wigert", f"{prefix}_n{n}_N{N}", m, n, N, 1, "continued", n, 2 * n, continuation="lognormal",
                    closed=lambda m0=m0, scale=scale, eta=eta, N=N: stieltjes_wigert(m0, scale, eta, F(0), N))
    # V12: the lognormal density with eta = 5/4, m0 = 2 and scale 2 shifted to
    # [lower, inf), lower = -1 and 4, with the bounds [lower, the largest
    # double], N = 1, 2, N_G = 20.
    for lower in (F(-1), F(4)):
        for n in (1, 2):
            m = exact_moments(lognormal_moments(F(2), F(2), F(5, 4), lower, 2 * n + 1))
            add("shifted_stieltjes_wigert", f"lower{float(lower):g}_n{n}_N20", m, n, 20, 1, "continued", n, 2 * n, continuation="lognormal",
                bounds=(lower, F(MAX_DOUBLE)), closed=lambda lower=lower: stieltjes_wigert(F(2), F(2), F(5, 4), lower, 20))
    # V13: a mixture of two lognormal number densities, moments rounded to
    # double, with the default bounds, N = 3, N_G = 20. Not a lognormal
    # density, so eta from m_0, m_1 and m_2 differs from the eta of the last
    # known zeta pair.
    add("lognormal_mixture", "lognormal_mixture_n3_N20", lognormal_mixture_moments(7), 3, 20, 1, "continued", 3, 6, continuation="lognormal")
    # V14: as V10 with the lognormal continuation.
    m = exact_moments(moments_of([F(-1), F(1), F(2)], [F(1, 4), F(1, 2), F(1, 4)], 7))
    add("lognormal_reduction", "outside_n3", m, 3, 9, 1, "continued", 1, 2, fixed_refused=True, continuation="lognormal")
    # V15: beta moments with the beta continuation and the default bounds give
    # the Gauss rule of the Jacobi polynomials: alpha = 1/2 and beta = -1/2
    # with m0 = 1, and alpha = 3/2 and beta = 1/2 with m0 = 2, N = 1, 2, 3,
    # N_G = 10 and 51. Their moments are exact in double.
    for m0, alpha, beta_, prefix in ((F(1), F(1, 2), F(-1, 2), "alpha0.5_beta-0.5"), (F(2), F(3, 2), F(1, 2), "alpha1.5_beta0.5")):
        for n in (1, 2, 3):
            m = exact_moments(beta_moments(m0, alpha, beta_, F(0), F(1), 2 * n + 1))
            for N in (10, 51):
                add("jacobi", f"{prefix}_n{n}_N{N}", m, n, N, 1, "continued", n, 2 * n, continuation="beta",
                    closed=lambda m0=m0, alpha=alpha, beta_=beta_, N=N: jacobi(m0, alpha, beta_, F(0), F(1), N))
    # V16: the beta density with alpha = 3/2, beta = 1/2 and m0 = 2 mapped to
    # [lower, upper] = [-1, 3] and [1/4, 3/4], with these bounds, N = 1, 2,
    # N_G = 20.
    for lower, upper in ((F(-1), F(3)), (F(1, 4), F(3, 4))):
        for n in (1, 2):
            m = exact_moments(beta_moments(F(2), F(3, 2), F(1, 2), lower, upper, 2 * n + 1))
            add("scaled_jacobi", f"lower{float(lower):g}_upper{float(upper):g}_n{n}_N20", m, n, 20, 1, "continued", n, 2 * n,
                continuation="beta", bounds=(lower, upper),
                closed=lambda lower=lower, upper=upper: jacobi(F(2), F(3, 2), F(1, 2), lower, upper, 20))
    # V17: four mixtures of beta number densities, moments rounded to double,
    # with the default bounds, N = 2, N_G = 20. Not beta densities, so alpha
    # and beta from p_1 and p_2 differ from those of the last known pair.
    # Every continued p of the first takes the first case of eqs. 34 and 35 by
    # its first condition, and every continued p of the second the second
    # case. The third takes the first case by its second condition for odd p,
    # and the fourth for even p. main() asserts that the four together decide
    # by each condition and take the second case for odd and even canonical
    # moments.
    for name, parts in (("beta_mixture_first_case_n2_N20", ((F(1, 2), F(2), F(5)), (F(1, 2), F(6), F(1)))),
                        ("beta_mixture_second_case_n2_N20", ((F(1, 4), F(3), F(3)), (F(3, 4), F(30), F(10)))),
                        ("beta_mixture_second_condition_odd_n2_N20", ((F(1, 2), F(1, 4), F(3)), (F(1, 2), F(5), F(5)))),
                        ("beta_mixture_second_condition_even_n2_N20", ((F(3, 4), F(-1, 4), F(2)), (F(1, 4), F(10), F(-1, 2))))):
        add("beta_mixture", name, beta_mixture_moments(parts, 5), 2, 20, 1, "continued", 2, 4, continuation="beta")
    # V18: moments on the boundary, of dyadic points in (0, 1), as V4.
    m = exact_moments(moments_of([F(1, 8), F(1, 2), F(3, 4)], [F(1, 4), F(1, 2), F(1, 4)], 7))
    add("beta_boundary", "three_points_n3", m, 3, 9, 1, "classical", 3, 6, continuation="beta")
    m = exact_moments(moments_of([F(1, 4), F(3, 4)], [F(3, 4), F(1, 4)], 9))
    add("beta_boundary", "two_points_n4", m, 4, 12, 1, "classical", 2, 8, fixed_refused=True, continuation="beta")
    # V19: the moments of three dyadic points, one outside [0, 1], at N = 3
    # with the default bounds. Their classical 3-node rule has that point as a
    # node. The continuation from the 2-node rule has a node above 1 for the
    # points 1/4, 1/2 and 5/4, and a negative b for the points -1/8, 1/4 and
    # 3/4, so adaptive mode continues the 1-node rule, and fixed mode is
    # not_realizable.
    for name, points in (("above_n3", [F(1, 4), F(1, 2), F(5, 4)]), ("negative_b_n3", [F(-1, 8), F(1, 4), F(3, 4)])):
        m = exact_moments(moments_of(points, [F(1, 4), F(1, 2), F(1, 4)], 7))
        add("beta_reduction", name, m, 3, 9, 1, "continued", 1, 2, fixed_refused=True, continuation="beta")
    return out


def assert_decisions(case):
    # The decisions QMOMKit takes on the moments in double, with margins.
    m, n = case["m"], case["n"]
    mp.mp.dps = 50
    if case["group"] in ("hermite", "chebyshev", "figures", "laguerre", "shifted_laguerre", "gamma_mixture", "stieltjes_wigert",
                         "shifted_stieltjes_wigert", "lognormal_mixture", "jacobi", "scaled_jacobi", "beta_mixture"):
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
    if case["group"] in ("gamma_reduction", "lognormal_reduction"):
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
        continued = gamma if case["continuation"] == "gamma" else lognormal
        x, _ = continued(m[:5], 2, case["N"], F(0), 2)
        assert x[0] < -1.0e-6 * max(abs(v) for v in x), case["name"]
        # Node count 1: the continued rule follows in main().
        assert diagonal[1][0] > 1000 * diagonal[1][1], case["name"]
    if case["group"] == "beta_mixture":
        # Every case of eqs. 34 and 35 is taken with a margin far beyond the
        # rounding of the canonical moments in double.
        _, _, choices = beta_recurrence(m, n, case["N"], F(0), F(1), n)
        assert all(margin > 1.0e-6 for _, _, margin in choices), case["name"]
    if case["group"] == "beta_reduction":
        # Node count 3: the classical rule, on the boundary, has the node
        # outside [0, 1].
        _, _, diagonal = recurrence(m, n)
        assert abs(diagonal[3][0]) <= diagonal[3][1], case["name"]
        x, _ = classical(m, 3)
        assert x[0] < -0.1 or x[-1] > 1.1, case["name"]
        # Node count 2: the classical rule is interior and accepted, and its
        # continuation has a node clearly above 1 or a clearly negative b.
        x, w = classical(m, 2)
        assert min(w) / max(w) > 10 * case["rmin"], case["name"]
        assert diagonal[2][0] > 1000 * diagonal[2][1], case["name"]
        ac, bc, _ = beta_recurrence(m[:5], 2, case["N"], F(0), F(1), 2)
        if case["name"] == "negative_b_n3":
            assert min(bc[1:]) < -0.1 * max(bc[1:]), case["name"]
        else:
            assert all(v > 0 for v in bc[1:]), case["name"]
            x, _ = gauss(ac, bc, m[0], case["N"])
            assert x[-1] > 1.1, case["name"]
        # Node count 1: the continued rule follows in main().
        assert diagonal[1][0] > 1000 * diagonal[1][1], case["name"]


def assert_inside(case, x):
    # The continued nodes lie above the lower bound with a margin far beyond
    # the tolerance 1e-12 max |xi_i| of abscissas_inside_bounds, and below the
    # upper bound, for the beta continuation with the same margin.
    if case["continuation"] in ("gamma", "lognormal", "beta") and case["kind"] == "continued":
        q = lambda v: mp.mpf(v.numerator) / v.denominator
        margin = 1.0e-9 * max(abs(v) for v in x)
        assert x[0] - q(lower_of(case)) > margin, case["name"]
        if case["continuation"] == "beta":
            assert q(upper_of(case)) - x[-1] > margin, case["name"]
        else:
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
    # The beta mixtures decide by each condition of the first case of eqs. 34
    # and 35 and take the second case, for both the odd and the even canonical
    # moments.
    taken = set()
    for case in CASES:
        if case["group"] == "beta_mixture":
            _, _, choices = beta_recurrence(case["m"], case["n"], case["N"], F(0), F(1), case["n"])
            taken |= {(j % 2, which) for j, which, _ in choices}
    assert taken == {(parity, decision) for parity in (0, 1) for decision in ("first condition", "second condition", "second case")}, taken
    fmt = lambda values: ", ".join(repr(float(v)) for v in values)
    lines = ["// Generated by generalized_quadrature_exact.py; do not edit.", ""]
    for group in ("hermite", "chebyshev", "figures", "boundary", "reduction", "laguerre", "shifted_laguerre", "gamma_mixture", "gamma_boundary",
                  "gamma_reduction", "stieltjes_wigert", "shifted_stieltjes_wigert", "lognormal_mixture", "lognormal_reduction", "jacobi",
                  "scaled_jacobi", "beta_mixture", "beta_boundary", "beta_reduction"):
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
                lines.append(f"     qmomkit::GeneralizedContinuation::{case['continuation']}, {bounds}}},")
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
