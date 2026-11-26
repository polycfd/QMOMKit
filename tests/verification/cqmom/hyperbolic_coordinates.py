"""Generate the references of hyperbolic_coordinates.cpp.

Six sets of weighted points in two dimensions, chosen so that every raw
moment is exactly a double. For each set, the raw moment tensor over the
moment extents of its hierarchy and, where the case compares nodes, the
terminal nodes of the CQMOM quadrature reconstructed from those moments: the
HyQMOM rule of a coordinate from its recurrence and closure in exact rational
arithmetic, then the Gauss quadrature of J_{n+1} with mpmath at 60 digits.

For the two cases whose conditional moments are polynomials of the hyperbolic
parent coordinate, the script solves the Vandermonde systems of the child
moments at the HyQMOM nodes at 60 digits and asserts that the conditional
moments equal the polynomials of the parent coordinate that the points were
built from. For every HyQMOM rule, it asserts
that each sigma(k, k) exceeds the rounding bound QMOMKit applies to it a
thousandfold, so that the moments lie well inside moment space. Prints the C++
arrays stored in the test. Requires mpmath.

    python3 hyperbolic_coordinates.py
"""

from fractions import Fraction as F

import mpmath as mp

mp.mp.dps = 60
EPS = F(1, 2**52)


def q(value):
    return mp.mpf(value.numerator) / value.denominator


def tensor(points, extents):
    # The raw moments sum_j W_j x_j^p y_j^r, flattened with the last index
    # fastest, as qmomkit::flatten_multi_index. Each must be exactly a double.
    values = []
    for p in range(extents[0]):
        for r in range(extents[1]):
            value = sum(w * x**p * y**r for w, x, y in points)
            assert F(float(value)) == value, (p, r)
            values.append(value)
    return values


def moment(values, extents, p, r):
    return values[p * extents[1] + r]


def hyqmom(m, n, scale):
    # The recurrence of m_0, ..., m_{2n} in exact arithmetic,
    # docs/developer/qmom.md, with its absolute-value bound, then the closure
    # a_n = mean(a_0, ..., a_{n-1}). scale is the coordinate scale s of
    # conditional moments, whose rounding floor is m_0 s^(2k), and zero for
    # moments that are data. Returns the Gauss quadrature of J_{n+1} and the
    # smallest ratio of sigma(k, k) to its rounding bound.
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
        floor = abs(m[0]) * scale ** (2 * k)
        ratio = sigma[(k, k)] / (100 * EPS * max(bound[(k, k)], floor))
        margin = ratio if margin is None else min(margin, ratio)
        b.append(sigma[(k, k)] / sigma[(k - 1, k - 1)])
        if k < n:
            a.append(sigma[(k, k + 1)] / sigma[(k, k)] - sigma[(k - 1, k)] / sigma[(k - 1, k - 1)])
    a.append(sum(a) / n)
    J = mp.matrix(n + 1, n + 1)
    for i in range(n + 1):
        J[i, i] = q(a[i])
    for i in range(1, n + 1):
        J[i - 1, i] = J[i, i - 1] = mp.sqrt(q(b[i]))
    E, Q = mp.eigsy(J)
    order = sorted(range(n + 1), key=lambda i: E[i])
    return [E[i] for i in order], [q(m[0]) * Q[0, i] ** 2 for i in order], margin


def conditional_moments(nodes, weights, rhs):
    # The solution c of sum_a w_a x_a^p c_a = rhs[p], p = 0, ..., N - 1, the
    # square Vandermonde system of the child moments at the parent nodes.
    size = len(nodes)
    V = mp.matrix(size, size)
    for p in range(size):
        for a in range(size):
            V[p, a] = weights[a] * nodes[a] ** p
    return mp.lu_solve(V, mp.matrix([q(value) for value in rhs]))


def scale(values, extents, axis):
    # The coordinate scale sqrt(m_2 / m_0) of the marginal of the axis,
    # src/coordinate_scales.cpp.
    m0 = moment(values, extents, 0, 0)
    m2 = moment(values, extents, 2, 0) if axis == 0 else moment(values, extents, 0, 2)
    return F(float(mp.sqrt(q(m2 / m0))))


def polynomial_mean():
    # x hyperbolic with n = 2 as parent, y classical with one node, and
    # y = f(x) = 1/2 + x/4 - x^2/2 at every point. The child reads only
    # moments of y order 1, so a spread of y about f(x) would not change them.
    f = lambda x: F(1, 2) + x / 4 - x**2 / 2
    xs = [F(-1), F(-1, 2), F(1, 4), F(1), F(2)]
    ws = [F(1, 8), F(1, 4), F(1, 4), F(1, 4), F(1, 8)]
    extents = (5, 2)
    values = tensor([(w, x, f(x)) for w, x in zip(ws, xs)], extents)
    nodes, weights, margin = hyqmom([moment(values, extents, p, 0) for p in range(5)], 2, 0)
    assert margin > 1000, float(margin)
    c = conditional_moments(nodes, weights, [moment(values, extents, p, 1) for p in range(3)])
    for a, node in enumerate(nodes):
        assert abs(c[a] - (q(F(1, 2)) + node / 4 - node**2 / 2)) < mp.mpf(10) ** -50
    terminals = [((a, 0), weights[a], (nodes[a], c[a])) for a in range(3)]
    return values, terminals


def polynomial_moments():
    # x hyperbolic with n = 3 as parent, y classical with two nodes, and the
    # two points (x, f(x) - 1/2) and (x, f(x) + 1/2) of equal weight at every
    # parent point, with f(x) = 1/4 + x/2. The conditional moments of y up to
    # order 3 are polynomials of degree at most 3 in x.
    f = lambda x: F(1, 4) + x / 2
    sd = F(1, 2)
    xs = [F(-3, 2), F(-1), F(-1, 4), F(1, 2), F(1), F(2)]
    ws = [F(1, 16), F(3, 16), F(1, 4), F(1, 4), F(1, 8), F(1, 8)]
    extents = (7, 4)
    points = [(w / 2, x, f(x) + sign * sd) for w, x in zip(ws, xs) for sign in (-1, 1)]
    values = tensor(points, extents)
    nodes, weights, margin = hyqmom([moment(values, extents, p, 0) for p in range(7)], 3, 0)
    assert margin > 1000, float(margin)
    for r in range(1, 4):
        c = conditional_moments(nodes, weights, [moment(values, extents, p, r) for p in range(4)])
        for a, node in enumerate(nodes):
            mean = q(F(1, 4)) + node / 2
            expected = ((mean - q(sd)) ** r + (mean + q(sd)) ** r) / 2
            assert abs(c[a] - expected) < mp.mpf(10) ** -50
    terminals = []
    for a, node in enumerate(nodes):
        mean = q(F(1, 4)) + node / 2
        terminals.append(((a, 0), weights[a] / 2, (node, mean - q(sd))))
        terminals.append(((a, 1), weights[a] / 2, (node, mean + q(sd))))
    return values, terminals


def classical_parent(branches):
    # x classical with n = 2 at exactly two points, so that its nodes are those
    # points, and y hyperbolic with n = 2 as child. Each branch is
    # (x, W, [(local weight, y)]), with local weights summing to one.
    extents = (4, 5)
    points = [(W * v, x, y) for x, W, ys in branches for v, y in ys]
    values = tensor(points, extents)
    s = scale(values, extents, 1)
    terminals = []
    for i, (x, W, ys) in enumerate(branches):
        m = [sum(v * y**r for v, y in ys) for r in range(5)]
        if len(ys) > 2:
            nodes, weights, margin = hyqmom(m, 2, s)
            assert margin > 1000, float(margin)
        else:
            nodes, weights = [q(y) for _, y in ys], [q(v) for v, _ in ys]
        terminals += [((i, a), q(W) * weights[a], (q(x), nodes[a])) for a in range(len(nodes))]
    return values, terminals


def last_coordinate():
    return classical_parent(
        [
            (F(-1, 2), F(3, 8), [(F(1, 4), F(-1)), (F(1, 4), F(0)), (F(1, 4), F(1, 2)), (F(1, 4), F(2))]),
            (F(3, 2), F(5, 8), [(F(1, 8), F(-2)), (F(1, 4), F(-1, 2)), (F(1, 4), F(1)), (F(1, 4), F(3, 2)), (F(1, 8), F(3))]),
        ]
    )


def boundary_branch():
    # The first branch of last_coordinate, and a second branch with exactly
    # two points, whose conditional moments lie on the boundary of moment
    # space for n = 2.
    return classical_parent(
        [
            (F(-1, 2), F(3, 8), [(F(1, 4), F(-1)), (F(1, 4), F(0)), (F(1, 4), F(1, 2)), (F(1, 4), F(2))]),
            (F(3, 2), F(5, 8), [(F(1, 4), F(-1)), (F(3, 4), F(2))]),
        ]
    )


def parent_on_the_boundary():
    # x hyperbolic with n = 2 at exactly two points, so that its moments lie on
    # the boundary of moment space and adaptive mode returns those points as the
    # classical rule of two nodes. y classical with two nodes, at two points for
    # each x. The children are conditioned on the two nodes with the parent
    # orders 0 and 1 of the layout, which also holds the parent order 2. The
    # terminal nodes are the points.
    branches = [
        (F(-1, 2), F(3, 8), [(F(1, 2), F(-1)), (F(1, 2), F(1))]),
        (F(3, 2), F(5, 8), [(F(1, 4), F(0)), (F(3, 4), F(2))]),
    ]
    extents = (5, 4)
    values = tensor([(W * v, x, y) for x, W, ys in branches for v, y in ys], extents)
    terminals = [((i, a), q(W * v), (q(x), q(y))) for i, (x, W, ys) in enumerate(branches) for a, (v, y) in enumerate(ys)]
    return values, terminals


def joint_moments():
    # x and y hyperbolic with n = 2, conditioning order {1, 0}, so y is the
    # parent. Three values of x at each of four values of y, with the mean of
    # x not a polynomial of degree 2 or less in y. Asserts that the y marginal
    # and the conditional moments of x at the three HyQMOM nodes of y lie well
    # inside moment space.
    rows = [
        (F(-1), [(F(1, 16), F(-1)), (F(1, 8), F(0)), (F(1, 16), F(1))]),
        (F(0), [(F(1, 16), F(-1, 2)), (F(1, 8), F(1, 2)), (F(1, 16), F(3, 2))]),
        (F(1), [(F(1, 16), F(0)), (F(1, 8), F(1)), (F(1, 16), F(2))]),
        (F(2), [(F(1, 16), F(-1)), (F(1, 8), F(1, 2)), (F(1, 16), F(2))]),
    ]
    extents = (5, 5)
    values = tensor([(w, x, y) for y, xs in rows for w, x in xs], extents)
    nodes, weights, margin = hyqmom([moment(values, extents, 0, r) for r in range(5)], 2, 0)
    assert margin > 1000, float(margin)
    s = scale(values, extents, 0)
    c = [conditional_moments(nodes, weights, [moment(values, extents, p, r) for r in range(3)]) for p in range(1, 5)]
    for a in range(3):
        m = [F(1)] + [F(mp.nstr(c[p][a], 55)) for p in range(4)]
        assert hyqmom(m, 2, s)[2] > 1000
    return values, None


def main():
    fmt = lambda values: ", ".join(repr(float(v)) for v in values)
    cases = [
        ("polynomial_mean", polynomial_mean),
        ("polynomial_moments", polynomial_moments),
        ("last_coordinate", last_coordinate),
        ("boundary_branch", boundary_branch),
        ("parent_on_the_boundary", parent_on_the_boundary),
        ("joint_moments", joint_moments),
    ]
    for name, case in cases:
        values, terminals = case()
        print(f"const std::vector<double> {name}_moments{{{fmt(values)}}};")
        if terminals is not None:
            print(f"const std::vector<ReferenceNode> {name}_nodes{{")
            for path, weight, abscissas in terminals:
                print(f"    {{{{{path[0]}, {path[1]}}}, {float(weight)!r}, {{{fmt(abscissas)}}}}},")
            print("};")


if __name__ == "__main__":
    main()
