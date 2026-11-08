# QMOM families

A QMOM family, `QmomOptions::family`
([`include/qmomkit/qmom.hpp`](../../include/qmomkit/qmom.hpp)), decides how the
moments of one coordinate become nodes. `classical`, the default, is the Wheeler
inversion of [QMOM](qmom.md): $m_0,\ldots,m_{2N-1}$ for up to $N$ nodes.
`hyperbolic` is HyQMOM: $m_0,\ldots,m_{2N}$ for up to $N+1$ nodes. Why a user
would choose a family is in the user guide's
[QMOM families](../user/7_qmom_families.md); this page describes where a family
enters the code and the mathematics of each family other than `classical`.

## Where a family enters the code

| Concern | Code | Behavior |
| --- | --- | --- |
| Moments read | `internal::qmom_required_moment_count` ([`src/internal/qmom_family.hpp`](../../src/internal/qmom_family.hpp)) | $2N$ or $2N+1$ for node count $N$, the one definition cases, state candidates, the evolution and the configuration use. |
| Maximum node count | `QmomState::maximum_node_count()` ([`include/qmomkit/qmom.hpp`](../../include/qmomkit/qmom.hpp)) | $N$ or $N+1$ for the case's node count $N$. |
| Reconstruction | `compute_qmom_from_moments` ([`src/public_api.cpp`](../../src/public_api.cpp)) | `classical` takes the paths of [QMOM](qmom.md), `hyperbolic` calls `compute_hyperbolic_quadrature` with `node_count` as $n$; abscissa bounds with `hyperbolic` are `invalid_input` before any reconstruction. |
| Evolution | `qmom_reconstruction_request` ([`src/phase_space_evolution.cpp`](../../src/phase_space_evolution.cpp)) | The count rule of [Phase-space evolution](phase_space_evolution.md#node-counts): the maximum count asks for the family's quadrature, a count below the maximum for a classical one. |
| Configuration | `parse_qmom_family`, `infer_qmom_initial_maximum_order` ([`src/configuration.cpp`](../../src/configuration.cpp)) | `[phase_space] qmom_family`, and the default `[initial] maximum_order` from the moments the family reads ([Configuration and output](configuration_and_output.md#mapping-to-objects)). |
| CQMOM | `compute_cqmom(const CqmomCase&)`, `internal::compute_cqmom_internal` ([`src/public_api.cpp`](../../src/public_api.cpp)) | Every family other than `classical` is `invalid_input` ([CQMOM](cqmom.md)). |
| Characteristic speeds | `compute_hyqmom_characteristic_speeds` ([`src/public_api.cpp`](../../src/public_api.cpp)) | `invalid_input` for every family other than `hyperbolic`, for which they are the [speeds](#characteristic-speeds) of the reconstruction of the stored moments. |

Every `switch` over `QmomFamily` has a case per family and no `default`, so the
compiler's `-Wswitch` names each place a new family must handle.

## Hyperbolic

HyQMOM [[Fox & Laurent, 2022](https://www.doi.org/10.1137/21M1406143)] closes
$m_{2n+1}$ from $m_0,\ldots,m_{2n}$ on the real line so that the moment system
of free transport is globally hyperbolic, as shown for $n\le9$ (Fox & Laurent,
2022, Theorems 4.6 and 4.8). The internal `compute_hyperbolic_quadrature`
([`src/internal/hyperbolic_quadrature.hpp`](../../src/internal/hyperbolic_quadrature.hpp))
returns the $(n+1)$-node quadrature of that closure, which reproduces
$m_0,\ldots,m_{2n}$ and the closed $m_{2n+1}$. It uses $s=0$, as
`compute_quadrature(...)` does. HyQMOM is defined on the real line, so abscissa
bounds give the message `The hyperbolic QMOM family takes no abscissa bounds.`

### Reconstruction

`compute_recurrence` runs the table of
[Wheeler inversion](qmom.md#wheeler-inversion) on $2n+1$ moments: one more
column in every row $k<n$, and row $n$, whose only entry is $\sigma_{n,n}$, with
$b_n=\sigma_{n,n}/\sigma_{n-1,n-1}$ and $\tau_n$. The entries for $k<n$ are
those of the $2n$-moment table. $a_n$ would need $m_{2n+1}$; the closure (Fox &
Laurent, 2022, eq. 4.16, or eq. 4.13 in standardized moments) sets it:

```math
a_n=\frac{1}{n}\sum_{k=0}^{n-1}a_k .
```

The quadrature is the Gauss quadrature of $J_{n+1}$, with diagonal
$a_0,\ldots,a_n$ and off-diagonal $\sqrt{b_1},\ldots,\sqrt{b_n}$, the matrix
associated with the closure (Fox & Laurent, 2022, Section 7). For moments
symmetric about their mean, every $a_k$ is the mean, the closure is exact, and
the quadrature is the classical $(n+1)$-node Gauss quadrature of
$m_0,\ldots,m_{2n}$ and the symmetric $m_{2n+1}$. As $b_n\to0$, $J_{n+1}$ splits
into $J_n$ and the entry $a_n$: the quadrature tends to the classical $n$-node
quadrature and a node at $a_n$ of vanishing weight. On the boundary of the
realizable set the moments determine the distribution, a sum of at most $n$
weighted points (Fox & Laurent, 2022, Section 2.4), and the classical quadrature
is that distribution.

The input checks are those of
[Input checks and the cutoff](qmom.md#input-checks-and-the-cutoff) with $2n+1$
moments, messages naming `node_count`, the message `Moment sequence has fewer
than 2*n+1 entries.`, and no one-node shortcut: $n=1$ gives the two nodes of
$J_2$. The rounding test runs to $k=n$:

| Moments | Adaptive | Fixed |
| --- | --- | --- |
| Interior: $\sigma_{k,k}>\tau_k$ for $k=1,\ldots,n$ | Quadrature of $J_{n+1}$, without `rmin` and `eabs` | Quadrature of $J_{n+1}$ |
| On or below the boundary: first $k$ with $\sigma_{k,k}\le\tau_k$ | The classical selection from $k$, exactly what `compute_quadrature(moments, n, options)` returns | `supports at most k nodes` or the negative $b$ |

Before the quadrature of $J_{n+1}$, both modes check $a_0,\ldots,a_n$ and
$b_1,\ldots,b_n$ as in [Fixed node count](qmom.md#fixed-node-count), and a
failed eigensolve is the result. If $\sigma_{n,n}$ overflows, so does $\tau_n$.
Adaptive mode then takes the classical selection, since $\sigma_{n,n}\le\tau_n$
holds; fixed mode skips the non-finite $\tau_n$ and reports a non-finite
recurrence, even where `compute_quadrature(moments, n, options)`, which does not
read $m_{2n}$, would succeed.

Adaptive mode applies `rmin` and `eabs` only to the classical quadratures.
Applied to the quadrature of $J_{n+1}$, they would reduce interior moments to a
quadrature outside the closure: the classical $n$-node quadrature does not
reproduce $m_{2n}$, and its closure is QMOM's, which is only weakly hyperbolic
(Fox & Laurent, 2022, Theorem 4.3). They are not needed either. The rounding
test already excludes nodes made of rounding, as in fixed mode, and near the
boundary the node of small weight tends to $a_n$, a mean of well-determined
coefficients, with a weight that vanishes continuously. So every interior result
is the HyQMOM quadrature, and every boundary result its limit. A reduced result
is always a classical quadrature of at most $n$ nodes, and only on or below the
boundary.

### Characteristic speeds

The moment system $\partial_t m_k+\partial_x m_{k+1}=0$, $k=0,\ldots,2n$, closed
by HyQMOM has the characteristic polynomial $P_{2n+1}=Q_nR_{n+1}$, with $Q_n$
the $n$-th monic orthogonal polynomial of the recurrence and

```math
R_{n+1}(\lambda)=(\lambda-a_n)\,Q_n(\lambda)-\beta_nQ_{n-1}(\lambda),\qquad \beta_n=\frac{2n+1}{n}\,b_n
```

(Fox & Laurent, 2022, Theorem 4.6, eqs. 4.11, 4.13 and 4.16), shown for $n\le9$
and postulated beyond; the verification checks it for sample moments up to
$n=10$ ([QMOM verification](../../tests/verification/qmom/README.md#hyqmom)). Its
$2n+1$ roots, the characteristic speeds, are the eigenvalues of $J_n$, diagonal
$a_0,\ldots,a_{n-1}$ and off-diagonal $\sqrt{b_1},\ldots,\sqrt{b_{n-1}}$, and of
$K_{n+1}$, which is $J_{n+1}$ with $\sqrt{\beta_n}$ in place of $\sqrt{b_n}$
(Fox & Laurent, 2022, eqs. 4.17 and 4.18); $\beta_n$ enters only the speeds,
never the quadrature. For $\beta_n>0$ the roots of $R_{n+1}$ strictly separate
those of $Q_n$ (Fox & Laurent, 2022, Theorem 4.8), so the speeds are distinct.

The public `compute_hyqmom_characteristic_speeds(state)`
([`src/public_api.cpp`](../../src/public_api.cpp)) returns the speeds of the
quadrature that `compute_qmom(...)` reconstructs from the state's stored moments
with its node count and options:

| Reconstruction of the stored moments | Speeds |
| --- | --- |
| $n+1$ nodes | The eigenvalues of $J_n$ and $K_{n+1}$: $2n+1$ values |
| A classical quadrature of $k\le n$ nodes | Its $k$ abscissas, the roots of $Q_k$, each twice, and its closure $a_k$, the mean of $a_0,\ldots,a_{k-1}$, which is the trace of $J_k$ over $k$ and so the mean of the abscissas: $2k+1$ values |

The second row is the $k$-node system on its own boundary, $\beta_k=0$, where
$R_{k+1}(\lambda)=(\lambda-a_k)\,Q_k(\lambda)$ (Fox & Laurent, 2022, Remark
4.9). At $k=n$ it is the first row with $\beta_n=0$, so the speeds are
continuous as the moments approach the boundary. Fox & Laurent (2022) define no
speeds on the boundary for $k<n$; returning those of the $k$-node system is
QMOMKit's choice.

The speeds count only the nodes that the stored moments resolve.
[Node advection](phase_space_evolution.md#node-advection) moves the nodes
without reconstructing them, so a state can hold more nodes than its moments
resolve, for example after drag has brought its nodes together; its speeds are
then those of the fewer nodes a reconstruction would give. For a state that
holds a classical quadrature, the speeds reproduce its abscissas only to
rounding, since they come from a reconstruction of its stored moments, not from
its nodes.

The internal `compute_hyperbolic_characteristic_speeds(moments, n, options)` and
`compute_hyperbolic_quadrature` share `hyperbolic_recurrence`: the input checks,
the recurrence of $m_0,\ldots,m_{2n}$ with the closure, and the decision of
[Reconstruction](#reconstruction). The speeds therefore cannot decide
differently from the reconstruction. A failed reconstruction is returned with
its status and message and no speeds, so fixed mode is `not_realizable` when the
moments no longer support $n+1$ nodes. For $n+1$ nodes, the speeds take $J_n$
and $K_{n+1}$ from the recurrence in place of the eigensolve of $J_{n+1}$; an
eigensolve that does not converge or a speed beyond the range of a double is
`not_realizable`. Each call records one `wheeler_calls`
([Performance counters](performance_counters.md)).

### Where the code lives

| Function | File | Role |
| --- | --- | --- |
| `compute_hyperbolic_quadrature` | [`src/wheeler.cpp`](../../src/wheeler.cpp) | Internal HyQMOM reconstruction, $s=0$. |
| `compute_hyqmom_characteristic_speeds` | [`src/public_api.cpp`](../../src/public_api.cpp) | Public HyQMOM characteristic speeds of a state: family dispatch, then the internal speeds of the stored moments with the state's node count and options. |
| `compute_hyperbolic_characteristic_speeds` | [`src/wheeler.cpp`](../../src/wheeler.cpp) | Internal speeds of $m_0,\ldots,m_{2n}$, those of the quadrature `hyperbolic_recurrence` decides. |
| `hyperbolic_recurrence`, `hyperbolic_closure` | [`src/wheeler.cpp`](../../src/wheeler.cpp) | The input checks, the recurrence with the closure and the decision that the reconstruction and the speeds share; the closure $a_n$. |

## Adding a family

A family is a value of `QmomFamily` with:

1. its moment count in `internal::qmom_required_moment_count` and its maximum in
   `QmomState::maximum_node_count()`;
2. its reconstruction, dispatched in `compute_qmom_from_moments`, with the
   options and bounds it accepts and what it returns where the moments do not
   support its maximum count. The reconstruction also sets the base node count
   of every successful result, because the evolution takes the node count of a
   state's own reconstructions from it;
3. what its maximum count reconstructs in the evolution's count rule, and how a
   state grows to it
   ([Phase-space evolution](phase_space_evolution.md#node-counts));
4. its `qmom_family` value and message in `parse_qmom_family`, and its default
   `[initial] maximum_order`;
5. its support or refusal in CQMOM;
6. unit tests of each of these, a verification case against an independent
   reference ([Testing](testing.md)), and public-header coverage;
7. a section here, and its row and section in the user guide's
   [QMOM families](../user/7_qmom_families.md).
