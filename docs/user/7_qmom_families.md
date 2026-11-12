# QMOM families

QMOM, CQMOM, and the methods derived from them share one idea: a quadrature that
reproduces the moments of a number density function. The number of coordinates
selects QMOM or CQMOM ([Concepts](3_concepts.md#qmom-and-cqmom)), whereas the
*family* decides how the moments of a coordinate become quadrature nodes. The
other pages of this guide assume the default, `classical` QMOM. A specific QMOM
family can be selected on the case, for instance

```cpp
qmom_case.qmom_options.family = qmomkit::QmomFamily::hyperbolic;
```

or in a file with `[phase_space] qmom_family = hyperbolic`.

| Family | Moments for node count $N$ | Nodes | Base node count | Bounds | Used with | Also gives |
| --- | --- | --- | --- | --- | --- | --- |
| `classical` | $m_0, \ldots, m_{2N-1}$ | up to $N$ | the number of nodes | yes | QMOM, CQMOM | — |
| `hyperbolic` | $m_0, \ldots, m_{2N}$ | up to $N + 1$ | $N$ for $N + 1$ nodes, otherwise the number of nodes | no | QMOM | characteristic speeds |
| `generalized` | $m_0, \ldots, m_{2N}$ | up to $N_G > N$, which you set | the number $k \le N$ of classical nodes it continues | no | QMOM | — |

## Hyperbolic

The hyperbolic version of QMOM is HyQMOM
[[Fox & Laurent, 2022](https://www.doi.org/10.1137/21M1406143)]. It is meant for
a coordinate that is the velocity of the population elements, such as the
velocity of the droplets of a spray, when your code also transports the moments
in physical space with that velocity. With classical QMOM, the equations of that
transport are only weakly hyperbolic: their characteristic speeds, the speeds at
which changes in the moments travel, coincide in pairs (Fox & Laurent, 2022,
Theorem 4.3). HyQMOM chooses the moment that closes the transport equations so
that the equations are hyperbolic, with real characteristic speeds that are
distinct inside the realizable set. For the same node count, HyQMOM uses one
moment and one node more than classical QMOM.

### Reconstruction

From $m_0, \ldots, m_{2N}$, so $2N + 1$ moments, HyQMOM reconstructs $N + 1$
quadrature nodes that reproduce these moments. The transport equation of
$m_{2N}$ contains $m_{2N+1}$, and the moment $m_{2N+1}$ computed from the
quadrature closes the transport equations. Fox & Laurent (2022) show that this
moment makes the transport equations hyperbolic for $N \le 9$; QMOMKit's tests
check hyperbolicity for sample moments up to $N = 10$.

Inside the realizable set, the reconstruction returns the $N + 1$ nodes, and
the `rmin` and `eabs` tests do not apply to these nodes. On the boundary of the
realizable set, where at most $N$ nodes carry the moments, or outside the set,
an adaptive reconstruction returns the quadrature that adaptive classical QMOM
reconstructs from $m_0, \ldots, m_{2N-1}$, with the `rmin` and `eabs` tests,
and a fixed reconstruction fails.

The base node count of the $N + 1$ nodes is $N$, because HyQMOM builds them
from the recurrence of the classical $N$-node quadrature of
$m_0, \ldots, m_{2N-1}$, extended by $m_{2N}$ and by the closure that gives
$m_{2N+1}$. When the reconstruction returns a classical quadrature, on the
boundary or outside the realizable set, the base node count is the number of
nodes of that quadrature.

A state of node count $N$ has at most $N + 1$ nodes, `maximum_node_count()`.
Having lost nodes, the state regains at most one node per step (see
[Evolving in time](5_evolving_in_time.md#adaptive-node-counts)), through the
classical quadratures of up to $N$ nodes to the $N + 1$ nodes of HyQMOM.

### Characteristic speeds

Your code needs the characteristic speeds to transport the moments, for example
for the fluxes and the time step of a finite-volume scheme, as in Section 6 of
Fox & Laurent (2022). A state of the hyperbolic family gives these speeds:

```cpp
const qmomkit::HyqmomCharacteristicSpeeds result = qmomkit::compute_hyqmom_characteristic_speeds(state);
if (result.ok())
{
    // The speeds are in ascending order: front(), the first element, is the
    // most negative speed, and back(), the last element, the most positive.
    const double leftmost = result.speeds.front();
    const double rightmost = result.speeds.back();
}
```

The speeds are those of the quadrature that QMOMKit reconstructs from the
state's moments with the state's options, in ascending order: $2N + 1$ speeds
for $N + 1$ nodes, and $2k + 1$ speeds for a classical quadrature of $k \le N$
nodes, namely the abscissas of that quadrature, each twice, and the arithmetic
mean of these abscissas. After advection has brought nodes closer together than
the moments resolve, the speeds are those of fewer nodes than the state holds,
and in fixed mode the computation of the speeds fails. For moments your code
changed, create a new state with `state_candidate(...)`
([Coupling QMOMKit](6_coupling_qmomkit.md#moments-changed-by-your-code)). The
Reference lists the [statuses](8_reference.md#statuses) and
[messages](8_reference.md#troubleshooting).

The [`hyqmom_characteristic_speeds`](../../examples/hyqmom_characteristic_speeds/README.md)
example computes the speeds of Fig. 1 of Fox & Laurent (2022), and
[`hyqmom_1d_drag_bgk-collisions`](../../examples/hyqmom_1d_drag_bgk-collisions/README.md)
evolves a velocity distribution with the hyperbolic family.

## Generalized

The generalized version of QMOM is GQMOM
[[Fox, Laurent & Passalacqua, 2023](https://doi.org/10.1016/j.jaerosci.2022.106096)].
GQMOM reconstructs a quadrature
of many more nodes, $N_G$, that still reproduces $m_0, \ldots, m_{2N}$.
By assuming an underlying type of distributions, e.g. a normal distribution, GQMOM
continues the recurrence of the orthogonal polynomials of the moments beyond
the coefficients the moments determine, so that the quadrature represents a
continuous number density rather than $N$ points. More nodes can improve, for instance, the
evaluation of source terms without more moments, but $N$ still controls the
accuracy of the moment method.

$N_G$ has no default and must exceed $N$:

```cpp
qmom_case.qmom_options.family = qmomkit::QmomFamily::generalized;
qmom_case.qmom_options.generalized.node_count = 3 * qmom_case.node_count;
```

or in a `.qmom` configuration file with `qmom_family = generalized` and
`generalized_node_counts`. Fox et al. (2023, Section 4.1) find $N_G = 3N$
adequate in most cases, and more for highly localized source terms.

At present, the only available distribution type for the quadrature continuation is `normal`, for a number
density on the real line, so the family takes no bounds. Its parameter $\nu$,
`generalized.nu` or `generalized_nu`, 1 by default, controls the tails: $\nu = 1$
gives Gaussian tails and, for the moments of a normal distribution, its
Gauss-Hermite quadrature, $\nu = 2$ exponential tails, and $\nu = 0$ a number
density on a bounded interval (Fox et al., 2023, Fig. 1). A negative $\nu$
drives the continued $b_i$ toward zero, and eventually several modes appear
(Fox et al., 2023, Section 3.1.2). A strongly negative $\nu$ can drive them
below what a double resolves, for example $\nu = -20$ with $N = 1$ and
$N_G = 101$: many nodes then coincide and carry zero weight, while the moments
are still reproduced.

### Reconstruction

An adaptive reconstruction first reconstructs a classical quadrature of
$k \le N$ nodes, as adaptive classical QMOM does from $m_0, \ldots, m_{2N-1}$,
with the `rmin` and `eabs` tests. Where $m_0, \ldots, m_{2k}$ lie inside the
realizable set, it continues the quadrature to $N_G$ nodes, which reproduce
$m_0, \ldots, m_{2k}$. On the boundary of the realizable set, the moments are
those of $k$ points, GQMOM is classical QMOM, and the result is the classical
$k$-node quadrature, as it is below the boundary. A fixed reconstruction takes
$k = N$, so it returns $N_G$ nodes, or the classical $N$ nodes on the boundary,
and fails where the moments do not support $N$ nodes or lie below the boundary.
The base node count of a result is $k$. `rmin` and `eabs` do not apply to the
continued nodes, whose weights in the tails are often tiny.

Each reconstruction solves an eigenproblem of $N_G \times N_G$ entries, so its
memory grows as $N_G^2$ and its time faster. An $N_G$ too large for the memory
makes `compute_qmom(...)` throw `std::bad_alloc` or `std::length_error`.

A state has at most $N_G$ nodes, `maximum_node_count()`. Having lost nodes, it
regains at most one base node per step and keeps $N_G$ nodes while its
moments are inside the realizable set.

The [`gqmom_trimodal_distribution`](../../examples/gqmom_trimodal_distribution/README.md)
example reconstructs the trimodal number density of Section 4.1 of Fox et al.
(2023) with GQMOM and with classical QMOM.
