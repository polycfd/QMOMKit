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

| Family | Moments for node count $N$ | Nodes | Bounds | Used with | Also gives |
| --- | --- | --- | --- | --- | --- |
| `classical` | $m_0, \ldots, m_{2N-1}$ | up to $N$ | yes | QMOM, CQMOM | — |
| `hyperbolic` | $m_0, \ldots, m_{2N}$ | up to $N + 1$ | no | QMOM | characteristic speeds |

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
