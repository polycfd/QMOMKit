# QMOM families

QMOM, CQMOM and the methods derived from them share one idea: a quadrature that
reproduces the moments of a number density function. The number of coordinates
selects QMOM or CQMOM ([Concepts](3_concepts.md#qmom-and-cqmom)); the *family*
decides how the moments of a coordinate become nodes. The rest of this guide
describes the default, `classical`. Select another on the case,

```cpp
qmom_case.qmom_options.family = qmomkit::QmomFamily::hyperbolic;
```

or in a file with `[phase_space] qmom_family = hyperbolic`.

| Family | Moments for node count $N$ | Nodes | Bounds | Methods | Also gives |
| --- | --- | --- | --- | --- | --- |
| `classical` | $m_0, \ldots, m_{2N-1}$ | up to $N$ | yes | QMOM, CQMOM | — |
| `hyperbolic` | $m_0, \ldots, m_{2N}$ | up to $N + 1$ | no | HyQMOM | characteristic speeds |

## Hyperbolic

The hyperbolic version of QMOM is HyQMOM
[[Fox & Laurent, 2022](https://www.doi.org/10.1137/21M1406143)]. It is meant for
a coordinate that is the velocity of the elements, such as the velocity of the
droplets of a spray, when your code also transports the moments in physical
space with that velocity. With the classical family, the equations of that
transport are only weakly hyperbolic: their characteristic speeds, the speeds at
which changes in the moments travel, coincide in pairs (Fox & Laurent, 2022,
Theorem 4.3). HyQMOM chooses the moment that closes them so that they are
hyperbolic, with real speeds that are distinct inside the realizable set, and it
uses one moment and one node more for the same node count. Fox & Laurent also
expect its nodes to evaluate phase-space terms, such as drag, at least as well
as classical QMOM with $N + 1$ nodes (section 7). For a coordinate on a
half-line or an interval, such as a size, they see no obvious advantage, and the
family takes no bounds.

### Reconstruction

From $m_0, \ldots, m_{2N}$, so $2N + 1$ moments, HyQMOM reconstructs $N + 1$
nodes that reproduce them; in a file, `[initial] maximum_order` defaults to
$2N$. Their $m_{2N+1}$ is the moment that closes the transport equations. Fox &
Laurent show that it makes them hyperbolic for $N \le 9$; QMOMKit's tests check
it for sample moments up to $N = 10$.

Inside the realizable set, the reconstruction returns the $N + 1$ nodes without
the `rmin` and `eabs` tests. On its boundary, where at most $N$ nodes carry the
moments, or outside it, an adaptive reconstruction returns the classical
quadrature of $m_0, \ldots, m_{2N-1}$, and a fixed one fails.

A state of node count $N$ has at most $N + 1$ nodes, `maximum_node_count()`.
Having lost nodes, it regains at most one per step
([Evolving in time](5_evolving_in_time.md#adaptive-node-counts)), through the
classical quadratures of up to $N$ nodes to the $N + 1$ nodes of HyQMOM.

### Characteristic speeds

A code that transports the moments needs the characteristic speeds, for example
for the fluxes and the time step of a finite-volume scheme, as in section 6 of
Fox & Laurent (2022). A state of the hyperbolic family gives them:

```cpp
const qmomkit::HyqmomCharacteristicSpeeds result = qmomkit::compute_hyqmom_characteristic_speeds(state);
if (result.ok())
{
    const double slowest = result.speeds.front();
    const double fastest = result.speeds.back();
}
```

The speeds are ascending and belong to the quadrature that the state's moments
reconstruct to with its options: $2N + 1$ for $N + 1$ nodes, and for a classical
quadrature of $k \le N$ nodes $2k + 1$, its abscissas, each twice, and their
mean. Where advection has brought nodes closer than the moments resolve, they
are the speeds of fewer nodes than the state holds, and in fixed mode they fail.
For moments your code changed, create a new state with `state_candidate(...)`
([Coupling QMOMKit](6_coupling_qmomkit.md#moments-changed-by-your-code)). The
[Reference](8_reference.md#statuses) lists the statuses and messages.

The [`hyqmom_characteristic_speeds`](../../examples/hyqmom_characteristic_speeds/README.md)
example computes the speeds of Fig. 1 of Fox & Laurent (2022), and
[`hyqmom_1d_drag_bgk-collisions`](../../examples/hyqmom_1d_drag_bgk-collisions/README.md)
evolves a velocity distribution with the family.
