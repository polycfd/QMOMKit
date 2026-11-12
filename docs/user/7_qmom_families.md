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

| Family | Use it for | Moments for node count $N$ | Nodes | Bounds | Used with |
| --- | --- | --- | --- | --- | --- |
| `classical` | the default, any coordinate | $m_0, \ldots, m_{2N-1}$ | up to $N$ | yes | QMOM, CQMOM |
| `hyperbolic` | a velocity coordinate whose moments your code transports in physical space | $m_0, \ldots, m_{2N}$ | up to $N + 1$ | no | QMOM |
| `generalized` | source terms that $N$ nodes evaluate poorly | $m_0, \ldots, m_{2N}$ | up to $N_G > N$ | no | QMOM |

How each family reconstructs its quadrature is described in the developer
docs' [QMOM families](../developer/qmom_families.md).

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

Inside the realizable set, HyQMOM returns $N + 1$ nodes with base node count
$N$. On the boundary of the realizable set or outside it, an adaptive
reconstruction returns the classical quadrature of $m_0, \ldots, m_{2N-1}$, and
a fixed reconstruction fails.

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

The speeds are those of the quadrature reconstructed from the state's moments:
$2N + 1$ speeds for $N + 1$ nodes, and $2k + 1$ for a classical quadrature of
$k \le N$ nodes. For moments your code changed, create a new state with
`state_candidate(...)`
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
From the same moments $m_0, \ldots, m_{2N}$, GQMOM reconstructs many more
nodes, $N_G$, which still reproduce these moments and represent a continuous
number density rather than $N$ points. More nodes can improve the evaluation
of source terms without more moments, but $N$ still controls the accuracy of
the moment method.

$N_G$ has no default and must exceed $N$:

```cpp
qmom_case.qmom_options.family = qmomkit::QmomFamily::generalized;
qmom_case.qmom_options.generalized.node_count = 3 * qmom_case.node_count;
```

or in a file with `[phase_space] qmom_family = generalized` and
`generalized_node_counts`.
Fox et al. (2023, Section 4.1) find $N_G = 3N$ adequate in most cases, and more
for highly localized source terms.

The only continuation is `normal`, for a number density on the real line, so
the family takes no bounds. Its parameter $\nu$, `generalized.nu` or
`generalized_nu`, 1 by default, controls the tails: $\nu = 1$ gives Gaussian
tails, $\nu = 2$ exponential tails, $\nu = 0$ a number density on a bounded
interval (Fox et al., 2023, Fig. 1), and a negative $\nu$ can give several
modes (Section 3.1.2).

An adaptive reconstruction first finds a classical quadrature of $k \le N$
nodes, with the `rmin` and `eabs` tests. Where $m_0, \ldots, m_{2k}$ lie inside
the realizable set, GQMOM continues it to $N_G$ nodes, which reproduce these
moments and have base node count $k$. `rmin` and `eabs` do not apply to the
$N_G$ nodes. Where these moments lie on or below the boundary, the result is
the classical quadrature itself. A fixed reconstruction takes $k = N$, so it
returns $N_G$ nodes or, on the boundary, the classical $N$ nodes, and fails
otherwise. Each reconstruction solves an eigenproblem of
$N_G \times N_G$ entries, so keep $N_G$ modest.
