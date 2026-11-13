# QMOM families

A QMOM family, `QmomOptions::family`
([`include/qmomkit/qmom.hpp`](../../include/qmomkit/qmom.hpp)), decides how the
moments of one coordinate become a quadrature. `classical`, the default, is
QMOM with the Wheeler inversion of [QMOM](qmom.md): $m_0,\ldots,m_{2N-1}$ for
up to $N$ nodes. `hyperbolic` is HyQMOM: $m_0,\ldots,m_{2N}$ for up to $N+1$
nodes. `generalized` is GQMOM: $m_0,\ldots,m_{2N}$ for up to $N_G>N$ nodes,
with $N_G$ set in `QmomOptions::generalized`. This page
describes where a family enters the code and the mathematics of each family
other than `classical`.

## Where a family enters the code

| Concern | Code | Behavior |
| --- | --- | --- |
| Moments read | `internal::qmom_required_moment_count` ([`src/internal/qmom_family.hpp`](../../src/internal/qmom_family.hpp)) | $2N$ for `classical` and $2N+1$ otherwise, for node count $N$, the one definition cases, state candidates, the evolution and the configuration use. |
| Maximum node count | `QmomState::maximum_node_count()` ([`include/qmomkit/qmom.hpp`](../../include/qmomkit/qmom.hpp)) | $N$, $N+1$ or $N_G$ for the case's node count $N$. |
| Reconstruction | `compute_qmom_from_moments` ([`src/public_api.cpp`](../../src/public_api.cpp)) | `classical` takes the paths of [QMOM](qmom.md), `hyperbolic` calls `compute_hyperbolic_quadrature` and `generalized` calls `compute_generalized_quadrature`, through `compute_bounded_generalized_quadrature` for the gamma continuation, each with `node_count` as $N$. Abscissa bounds with `hyperbolic` or the normal continuation are `invalid_input` before any reconstruction. |
| Evolution | `holds_full_quadrature`, `full_quadrature_request`, `reduced_quadrature_request`, `current_reconstruction_request` ([`src/phase_space_evolution.cpp`](../../src/phase_space_evolution.cpp)) | The count rule of [Phase-space evolution](phase_space_evolution.md#node-counts). The family's full quadrature is requested from the family. A reduced quadrature is requested from the classical family for `classical` and `hyperbolic`, and from the family itself for `generalized`. |
| Configuration | `parse_qmom_family`, `parse_generalized_options`, `infer_qmom_initial_maximum_order` ([`src/configuration.cpp`](../../src/configuration.cpp)) | `[phase_space] qmom_family`, the `generalized_node_counts`, `generalized_continuation` and `generalized_nu` of the generalized family, and the default `[initial] maximum_order` from the moments the family reads ([Configuration and output](configuration_and_output.md#mapping-to-objects)). |
| CQMOM | `compute_cqmom(const CqmomCase&)`, `internal::compute_cqmom_internal` ([`src/public_api.cpp`](../../src/public_api.cpp)) | Every family other than `classical` is `invalid_input` ([CQMOM](cqmom.md)). |
| Characteristic speeds | `compute_hyqmom_characteristic_speeds` ([`src/public_api.cpp`](../../src/public_api.cpp)) | `invalid_input` for every family other than `hyperbolic`, for which they are the [speeds](#characteristic-speeds) of the reconstruction of the stored moments. |

Every `switch` over `QmomFamily` has a case per family and no `default`, so the
compiler's `-Wswitch` names each place a new family must handle.

## Hyperbolic

HyQMOM [[Fox & Laurent, 2022](https://www.doi.org/10.1137/21M1406143)] closes
$m_{2N+1}$ from $m_0,\ldots,m_{2N}$ on the real line so that the moment system
of free transport is globally hyperbolic, as shown for $N\le9$ (Fox & Laurent,
2022, Theorems 4.6 and 4.8). The internal `compute_hyperbolic_quadrature`
([`src/internal/hyperbolic_quadrature.hpp`](../../src/internal/hyperbolic_quadrature.hpp))
returns the $(N+1)$-node quadrature of that closure, which reproduces
$m_0,\ldots,m_{2N}$ and the closed $m_{2N+1}$. It uses $s=0$, as
`compute_quadrature(...)` does. HyQMOM is defined on the real line, so abscissa
bounds give the message `The hyperbolic QMOM family takes no abscissa bounds.`

### Reconstruction

`compute_recurrence` runs the table of
[Wheeler inversion](qmom.md#wheeler-inversion) on $2N+1$ moments: one more
column in every row $k<N$, and row $N$, whose only entry is $\sigma_{N,N}$, with
$b_N=\sigma_{N,N}/\sigma_{N-1,N-1}$ and $\tau_N$. The entries for $k<N$ are
those of the $2N$-moment table. $a_N$ would need $m_{2N+1}$; the closure (Fox &
Laurent, 2022, eq. 4.16, or eq. 4.13 in standardized moments) sets it:

```math
a_N=\frac{1}{N}\sum_{k=0}^{N-1}a_k .
```

The quadrature is the Gauss quadrature of $J_{N+1}$, with diagonal
$a_0,\ldots,a_N$ and off-diagonal $\sqrt{b_1},\ldots,\sqrt{b_N}$, the matrix
associated with the closure (Fox & Laurent, 2022, Section 7). For moments
symmetric about their mean, every $a_k$ is the mean, the closure is exact, and
the quadrature is the classical $(N+1)$-node Gauss quadrature of
$m_0,\ldots,m_{2N}$ and the symmetric $m_{2N+1}$. As $b_N\to0$, $J_{N+1}$ splits
into $J_N$ and the entry $a_N$: the quadrature tends to the classical $N$-node
quadrature and a node at $a_N$ of vanishing weight. On the boundary of the
realizable set the moments determine the distribution, a sum of at most $N$
weighted points (Fox & Laurent, 2022, Section 2.4), and the classical quadrature
is that distribution.

The input checks are those of
[Input checks and the cutoff](qmom.md#input-checks-and-the-cutoff) with $2N+1$
moments, messages naming `node_count`, the message `Moment sequence has fewer
than 2*n+1 entries.`, and no one-node shortcut: $N=1$ gives the two nodes of
$J_2$. The rounding test runs to $k=N$:

| Moments | Adaptive | Fixed |
| --- | --- | --- |
| Interior: $\sigma_{k,k}>\tau_k$ for $k=1,\ldots,N$ | Quadrature of $J_{N+1}$, without `rmin` and `eabs` | Quadrature of $J_{N+1}$ |
| On or below the boundary: first $k$ with $\sigma_{k,k}\le\tau_k$ | The classical selection from $k$, exactly what `compute_quadrature(moments, N, options)` returns | `supports at most k nodes` or the negative $b$ |

Before the quadrature of $J_{N+1}$, both modes check $a_0,\ldots,a_N$ and
$b_1,\ldots,b_N$ as in [Fixed node count](qmom.md#fixed-node-count), and a
failed eigensolve is the result. If $\sigma_{N,N}$ overflows, so does $\tau_N$.
Adaptive mode then takes the classical selection, since $\sigma_{N,N}\le\tau_N$
holds; fixed mode skips the non-finite $\tau_N$ and reports a non-finite
recurrence, even where `compute_quadrature(moments, N, options)`, which does not
read $m_{2N}$, would succeed.

Adaptive mode applies `rmin` and `eabs` only to the classical quadratures.
Applied to the quadrature of $J_{N+1}$, they would reduce interior moments to a
quadrature outside the closure: the classical $N$-node quadrature does not
reproduce $m_{2N}$, and its closure is QMOM's, which is only weakly hyperbolic
(Fox & Laurent, 2022, Theorem 4.3). They are not needed either. The rounding
test already excludes nodes made of rounding, as in fixed mode, and near the
boundary the node of small weight tends to $a_N$, a mean of well-determined
coefficients, with a weight that vanishes continuously. So every interior result
is the HyQMOM quadrature, and every boundary result its limit. A reduced result
is always a classical quadrature of at most $N$ nodes, and only on or below the
boundary.

### Characteristic speeds

The moment system $\partial_t m_k+\partial_x m_{k+1}=0$, $k=0,\ldots,2N$, closed
by HyQMOM has the characteristic polynomial $P_{2N+1}=Q_NR_{N+1}$, with $Q_N$
the $N$-th monic orthogonal polynomial of the recurrence and

```math
R_{N+1}(\lambda)=(\lambda-a_N)\,Q_N(\lambda)-\beta_NQ_{N-1}(\lambda),\qquad \beta_N=\frac{2N+1}{N}\,b_N
```

(Fox & Laurent, 2022, Theorem 4.6, eqs. 4.11, 4.13 and 4.16), shown for $N\le9$
and postulated beyond; the verification checks it for sample moments up to
$N=10$ ([QMOM verification](../../tests/verification/qmom/README.md#hyqmom)). Its
$2N+1$ roots, the characteristic speeds, are the eigenvalues of $J_N$, diagonal
$a_0,\ldots,a_{N-1}$ and off-diagonal $\sqrt{b_1},\ldots,\sqrt{b_{N-1}}$, and of
$K_{N+1}$, which is $J_{N+1}$ with $\sqrt{\beta_N}$ in place of $\sqrt{b_N}$
(Fox & Laurent, 2022, eqs. 4.17 and 4.18); $\beta_N$ enters only the speeds,
never the quadrature. For $\beta_N>0$ the roots of $R_{N+1}$ strictly separate
those of $Q_N$ (Fox & Laurent, 2022, Theorem 4.8), so the speeds are distinct.

The public `compute_hyqmom_characteristic_speeds(state)`
([`src/public_api.cpp`](../../src/public_api.cpp)) returns the speeds of the
quadrature that `compute_qmom(...)` reconstructs from the state's stored moments
with its node count and options:

| Reconstruction of the stored moments | Speeds |
| --- | --- |
| $N+1$ nodes | The eigenvalues of $J_N$ and $K_{N+1}$: $2N+1$ values |
| A classical quadrature of $k\le N$ nodes | Its $k$ abscissas, the roots of $Q_k$, each twice, and its closure $a_k$, the mean of $a_0,\ldots,a_{k-1}$, which is the trace of $J_k$ over $k$ and so the mean of the abscissas: $2k+1$ values |

The second row is the $k$-node system on its own boundary, $\beta_k=0$, where
$R_{k+1}(\lambda)=(\lambda-a_k)\,Q_k(\lambda)$ (Fox & Laurent, 2022, Remark
4.9). At $k=N$ it is the first row with $\beta_N=0$, so the speeds are
continuous as the moments approach the boundary. Fox & Laurent (2022) define no
speeds on the boundary for $k<N$; returning those of the $k$-node system is
QMOMKit's choice.

The speeds count only the nodes that the stored moments resolve.
[Node advection](phase_space_evolution.md#node-advection) moves the nodes
without reconstructing them, so a state can hold more nodes than its moments
resolve, for example after drag has brought its nodes together; its speeds are
then those of the fewer nodes a reconstruction would give. For a state that
holds a classical quadrature, the speeds reproduce its abscissas only to
rounding, since they come from a reconstruction of its stored moments, not from
its nodes.

The internal `compute_hyperbolic_characteristic_speeds(moments, N, options)` and
`compute_hyperbolic_quadrature` share `hyperbolic_recurrence`: the input checks,
the recurrence of $m_0,\ldots,m_{2N}$ with the closure, and the decision of
[Reconstruction](#reconstruction). The speeds therefore cannot decide
differently from the reconstruction. A failed reconstruction is returned with
its status and message and no speeds, so fixed mode is `not_realizable` when the
moments no longer support $N+1$ nodes. For $N+1$ nodes, the speeds take $J_N$
and $K_{N+1}$ from the recurrence in place of the eigensolve of $J_{N+1}$; an
eigensolve that does not converge or a speed beyond the range of a double is
`not_realizable`. Each call records one `wheeler_calls`
([Performance counters](performance_counters.md)).

### Where the code lives

| Function | File | Role |
| --- | --- | --- |
| `compute_hyperbolic_quadrature` | [`src/wheeler.cpp`](../../src/wheeler.cpp) | Internal HyQMOM reconstruction, $s=0$. |
| `compute_hyqmom_characteristic_speeds` | [`src/public_api.cpp`](../../src/public_api.cpp) | Public HyQMOM characteristic speeds of a state: family dispatch, then the internal speeds of the stored moments with the state's node count and options. |
| `compute_hyperbolic_characteristic_speeds` | [`src/wheeler.cpp`](../../src/wheeler.cpp) | Internal speeds of $m_0,\ldots,m_{2N}$, those of the quadrature `hyperbolic_recurrence` decides. |
| `hyperbolic_recurrence`, `hyperbolic_closure` | [`src/wheeler.cpp`](../../src/wheeler.cpp) | The input checks, the recurrence with the closure and the decision that the reconstruction and the speeds share; the closure $a_N$. |

## Generalized

GQMOM [[Fox, Laurent & Passalacqua, 2023](https://doi.org/10.1016/j.jaerosci.2022.106096)]
reconstructs $m_0,\ldots,m_{2N}$ like QMOM and continues the recurrence of their
orthogonal polynomials beyond the coefficients the moments determine, so that
the quadrature has $N_G>N$ nodes and still reproduces $m_0,\ldots,m_{2N}$ (Fox
et al., 2023, eq. 14). The continuation chooses a number density among the many
with these moments, and its higher moments lie in the interior of moment space.
The internal `compute_generalized_quadrature`
([`src/internal/generalized_quadrature.hpp`](../../src/internal/generalized_quadrature.hpp))
returns the quadrature with its base node count and uses $s=0$, as
`compute_quadrature(...)` does. `QmomOptions::generalized` holds $N_G$, which
has no default, the continuation, and the parameter $\nu$ of the normal
continuation. The `normal` continuation is for a number density on the real
line, so abscissa bounds give the message
`The normal continuation of the generalized QMOM family takes no abscissa bounds.`
The `gamma` continuation is for a number density on $[\ell,\infty)$, with
$\ell$ the lower abscissa bound, and its nodes must lie within the bounds
([Gamma continuation](#gamma-continuation)).

### Reconstruction

`compute_recurrence` runs the table of
[Wheeler inversion](qmom.md#wheeler-inversion) on $2N+1$ moments. Its entries
for $k<N$ do not depend on $m_{2N}$, so they are those of the $2N$-moment table,
and row $N$, whose only entry is $\sigma_{N,N}$, adds $b_N$ and $\tau_N$. First,
a classical quadrature of $k\le N$ nodes is found:

| Mode | Classical quadrature |
| --- | --- |
| Adaptive | The classical selection from the first $k\le N$ with $\sigma_{k,k}\le\tau_k$, or from $N$, with `rmin` and `eabs`: exactly what `compute_quadrature(moments, N, options)` returns |
| Fixed | The checks of [Fixed node count](qmom.md#fixed-node-count) for $N$ nodes, then the $N$-node rule |

Then, since $\sigma_{k,k}=m_0b_1\cdots b_k$, the test $\sigma_{k,k}>\tau_k$
decides whether $b_k$ is resolvably positive. If it is, the $k$-node quadrature
is continued to $N_G$ nodes. Fixed mode gives `not_realizable` with the message
`Wheeler recurrence produced non-finite values.` where $\sigma_{k,k}$ or
$\tau_k$ is not finite. In adaptive mode, an overflowed $\sigma_{k,k}=\pm\infty$
has an infinite $\tau_k$, so $\sigma_{k,k}\le\tau_k$ holds and the result is
the classical $k$-node quadrature, whose coefficients are finite. Whether the
overflow gives $\pm\infty$ or NaN can depend on the build, for example on
whether it fuses multiply-adds. A NaN fails $\sigma_{k,k}\le\tau_k$, and
adaptive mode then reports the non-finite $b_k$ of the continuation as
`not_realizable`. If $\sigma_{k,k}$ is zero within $\tau_k$, the moments lie on
the boundary of moment space at $k$ nodes, where GQMOM is QMOM (Fox et al.,
2023, Section 5), and the result is the classical $k$-node quadrature. Below
$-\tau_k$, $b_k$ is negative and the moments lie below the boundary: adaptive
mode returns the classical $k$-node quadrature, as the classical selection
does, and fixed mode is `not_realizable` with the message of a negative $b$.
The base node count is $k$, so fixed mode always returns base node count $N$:
$N_G$ nodes in the interior and the classical $N$ nodes on the boundary.

The normal continuation (Fox et al., 2023, Definition 1, eqs. 21 and 22)
continues the recurrence of $k$ nodes with

```math
a_i=\frac{1}{k}\sum_{j=0}^{k-1}a_j\quad(i\ge k),\qquad b_i=\Big(\frac{i}{k}\Big)^{\nu}b_k\quad(i>k).
```

Eq. 22 writes $a_i$ as the mean of $a_0,\ldots,a_{i-1}$, which equals the mean
of $a_0,\ldots,a_{k-1}$ for every $i\ge k$, so the code computes it once. For
Gaussian moments and $\nu=1$, these are the coefficients of the Hermite
polynomials (eq. 20), and the quadrature is the Gauss-Hermite rule. Larger
$\nu$ gives heavier tails, and $\nu=0$ with $k=1$ gives the Chebyshev
polynomials of the second kind on a finite interval (Fox et al., 2023, Fig. 1
and footnote 2). The quadrature is the Gauss quadrature of $J_{N_G}$, with
diagonal $a_0,\ldots,a_{N_G-1}$ and off-diagonal
$\sqrt{b_1},\ldots,\sqrt{b_{N_G-1}}$, from the dense eigensolver of
[Linear algebra](linear_algebra.md). Its memory grows as $N_G^2$ and its time
faster. Since $N_G$ has no cap, an $N_G$ too large for the memory makes
`compute_qmom(...)` throw `std::bad_alloc` or `std::length_error`. Before the
eigensolve, the continued coefficients are checked as in
[Fixed node count](qmom.md#fixed-node-count), so an overflow, for example of
$(i/k)^\nu$ for a large $\nu$, gives `not_realizable` with the message
`Wheeler recurrence produced non-finite values.`

A strongly negative $\nu$ makes the continued $b_i$ tiny. Where their square
roots fall below the resolution of the diagonal $a_i$, many continued nodes
coincide in double precision and their weights are zero, for example at
$k=1$, $N_G=101$ and $\nu=-20$. The result is still a success and reproduces
$m_0,\ldots,m_{2k}$ to rounding. This is a shortcoming of double precision,
not of the method, and QMOMKit leaves it as it is.

### Gamma continuation

The gamma continuation (Fox et al., 2023, Definition 2, eqs. 23–26) is for a
number density on $[\ell,\infty)$, with $\ell$ the lower abscissa bound.
Subtracting $\ell$ from every $a_i$ shifts the recurrence exactly to
$[0,\infty)$, where it is written with
$a_0=\zeta_1$, $a_i=\zeta_{2i}+\zeta_{2i+1}$ and $b_i=\zeta_{2i-1}\zeta_{2i}$
(eq. 23). The recurrence of $k$ nodes gives $\zeta_1,\ldots,\zeta_{2k}$, and
Definition 2's $\alpha=m_1^2/(m_2m_0-m_1^2)-1$ is $\zeta_1/\zeta_2-1$ in the
shifted moments. The continuation (eq. 25, with $k$ in place of $n$, as for the
normal continuation) is

```math
\zeta_{2i-1}=\frac{i+\alpha}{k+\alpha}\,\zeta_{2k-1},\qquad \zeta_{2i}=\frac{i}{k}\,\zeta_{2k}\qquad(i>k),
```

and eq. 26 gives $a_i=\zeta_{2i}+\zeta_{2i+1}$ for $i\ge k$ and
$b_i=\zeta_{2i-1}\zeta_{2i}$ for $i>k$, shifted back by $\ell$. Eq. 26 is
stated for $i>n$, but $a_n$ needs the continued $\zeta_{2n+1}$ as well, as the
paper's selection of every $\zeta_i$ with $i>2n$ implies. For the moments of a
gamma number density, these are the coefficients of the Laguerre polynomials
(eq. 24), from every $k$.

The abscissa bounds $[\ell,u]$ are the support, $[0,\text{DBL\_MAX}]$ without
bounds. `compute_bounded_generalized_quadrature` selects as bounded classical
QMOM does ([Abscissa bounds](qmom.md#abscissa-bounds)): a reconstruction is
accepted when it succeeds and its nodes pass `abscissas_inside_bounds`, the
check of bounded classical QMOM. Fixed mode tries node count $N$ only. Adaptive
mode retries below the base node count of a rejected reconstruction, down to
one node, so the continuation runs once per base node count. Moments that are
not those of a number density on $[\ell,\infty)$ give a non-positive $\zeta$,
so the continued recurrence fails the check of a negative or non-finite $b$, or
its nodes leave the bounds. Both reject the reconstruction. The continuation
does not read $\nu$.

`rmin` and `eabs` act on the classical quadratures only. The continued nodes
carry the tails of the number density, whose weights can be far below `rmin`
(Fox et al., 2023, Fig. 1), and the test on $b_k$ already excludes nodes made of
rounding.

$N_G$ is carried in the options, independently of $N$ and $k$, so a reduced
quadrature is continued to the same $N_G$ nodes as the full one. A
reconstruction with node count $k$ of the same moments that finds the same $k$
nodes returns the identical result, which is what the evolution relies on
([Phase-space evolution](phase_space_evolution.md#node-counts)).

The input checks, in this order:

| Condition | Status | Message |
| --- | --- | --- |
| abscissa bounds with the normal continuation | `invalid_input` | `The normal continuation of the generalized QMOM family takes no abscissa bounds.` |
| abscissa bounds with a non-finite value or $\ell\ge u$, with the gamma continuation | `invalid_input` | `Bounded-support interval must have finite lower < upper values.` |
| non-finite $\nu$ with the normal continuation | `invalid_input` | `The generalized parameter nu must be finite.` |
| $N_G\le N$, which includes an unset $N_G=0$ | `invalid_input` | `The generalized node count must exceed node_count.` |
| a continuation outside `GeneralizedContinuation` | `invalid_input` | `Unsupported generalized continuation.` |
| the checks of [Input checks and the cutoff](qmom.md#input-checks-and-the-cutoff) with $2N+1$ moments and no one-node shortcut, with the message `Moment sequence has fewer than 2*n+1 entries.` | | |

### Where the code lives

| Function | File | Role |
| --- | --- | --- |
| `compute_generalized_quadrature` | [`src/wheeler.cpp`](../../src/wheeler.cpp) | Internal GQMOM reconstruction, $s=0$, with its base node count. |
| `is_supported_continuation`, `normal_continuation`, `gamma_continuation` | [`src/wheeler.cpp`](../../src/wheeler.cpp) | The check of the continuation; the normal continuation of eqs. 21 and 22; the gamma continuation of eqs. 23–26. |
| `compute_bounded_generalized_quadrature`, `abscissas_inside_bounds` | [`src/wheeler.cpp`](../../src/wheeler.cpp) | The selection within abscissa bounds for the gamma continuation, with the node check of bounded classical QMOM. |
| `input_failure_through_m2n` | [`src/wheeler.cpp`](../../src/wheeler.cpp) | The input checks of $m_0,\ldots,m_{2N}$, shared with HyQMOM. |

## Adding a family

A family is a value of `QmomFamily` with:

1. its moment count in `internal::qmom_required_moment_count` and its maximum in
   `QmomState::maximum_node_count()`;
2. its reconstruction, dispatched in `compute_qmom_from_moments`, with the
   options and bounds it accepts and what it returns where the moments do not
   support its maximum count. The reconstruction also sets the base node count
   of every successful result, because the evolution takes the node count of a
   state's own reconstructions from it;
3. its full quadrature and its reduced quadratures, as listed in
   [Phase-space evolution](phase_space_evolution.md#node-counts).
   `full_quadrature_request` requests the full quadrature from the family's
   reconstruction with node count $N$. `reduced_quadrature_request` requests a
   reduced quadrature with base node count $k$, from the classical family for
   `classical` and `hyperbolic` and from the family itself for `generalized`.
   A reduced state grows back by at most one base node per step up to $N$, and
   where the full quadrature has more than $N$ nodes, one more accepted trial
   restores it;
4. its `qmom_family` value and message in `parse_qmom_family`, the keys of its
   settings, and its default `[initial] maximum_order`;
5. its support or refusal in CQMOM;
6. unit tests of each of these, a verification case against an independent
   reference ([Testing](testing.md)), and public-header coverage;
7. a section here, and its row and section in the user guide's
   [QMOM families](../user/7_qmom_families.md).
