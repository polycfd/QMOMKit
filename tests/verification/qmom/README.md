# One-Dimensional QMOM Analytical Verification

This layer verifies the public `qmomkit::compute_qmom` API from inputs loaded
through the public configuration API. The configuration files live in
`tests/verification/qmom/` and contain only QMOM metadata and raw input moments.
Expected quadrature weights and abscissas are defined independently in
`tests/verification/qmom/expected.hpp`.

## Manufactured Discrete Cases

| Identifier | Configuration | Exact weights | Exact abscissas |
| --- | --- | --- | --- |
| `single_dirac` | `single_dirac.qmom` | `2.5` | `-0.7` |
| `two_node_discrete` | `two_node_discrete.qmom` | `1/3, 2/3` | `1, 3` |
| `three_node_discrete` | `three_node_discrete.qmom` | `0.2, 0.3, 0.5` | `-1.5, 0.4, 2.2` |

These cases verify exact inversion of already-discrete raw moment sequences,
including preservation of non-unit mass.

## Continuous Distribution Cases

Uniform distribution on `[0, 1]`:

```text
f(xi) = 1
m_k = 1 / (k + 1)
```

The expected quadrature is Gauss-Legendre on `[-1, 1]` mapped to `[0, 1]` by
`xi = (x + 1) / 2`; weights are scaled by `1/2`. The verification includes
`uniform_n2.qmom` and `uniform_n3.qmom`.

Exponential distribution on `[0, infinity)`:

```text
f(xi) = exp(-xi)
m_k = k!
```

The expected quadrature is Gauss-Laguerre with `alpha = 0`. The verification
includes `exponential_n2.qmom` and `exponential_n3.qmom`.

Standard normal distribution:

```text
f(xi) = exp(-xi^2 / 2) / sqrt(2*pi)
m_(2k) = (2k - 1)!!
m_(2k+1) = 0
```

The expected quadrature uses the probabilists' Gauss-Hermite convention for the
standard normal weight. Do not compare these nodes directly to physicists'
Gauss-Hermite nodes for `exp(-x^2)`: the standard-normal/probabilists nodes are
larger by a factor of `sqrt(2)` after converting the weight convention. The
verification includes `gaussian_n2.qmom` and `gaussian_n3.qmom`.

Each advection case uses exactly the first `2N` raw moments, `m_0` through
`m_(2N-1)`. For these three advection families, `m_0 = 1` because the chosen
weight functions are normalized probability densities, not because QMOMKit
normalizes raw moments.

## Checks

The reusable executable `qmom_analytical_verification`:

- loads each input moment sequence from `.qmom` through
  `qmomkit::config::load_qmom`;
- calls public `qmomkit::compute_qmom(...)`;
- sorts recovered nodes by abscissa before direct comparison;
- compares weights and nodes against independent classical quadrature values;
- reconstructs every configured raw moment from the recovered quadrature;
- reports maximum absolute and relative reconstruction errors;
- verifies all weights are positive;
- verifies `sum(weights) = m_0` without assuming unit mass.

`qmom_compute_failures.cpp` verifies direct public `compute_qmom(...)`
failure contracts without going through configuration parsing. These cases are
executable documentation for user code that constructs `MomentSet` and
`node_count` directly:

- `m0 < 0` is a non-realizable moment sequence and must return
  `QmomStatus::not_realizable` with an explanatory message;
- a negative-variance moment sequence is mathematically impossible and must
  return `QmomStatus::not_realizable`;
- too few raw moments for the requested node count must return
  `QmomStatus::invalid_input`;
- requesting zero nodes must return `QmomStatus::invalid_input`.

The failure verifier checks `ok()`, status, and message content. Structurally
invalid calls must not expose fallback quadrature data. Non-realizable calls
may carry diagnostic fallback values, but they must never be reported as
successful quadrature.

Distribution-specific checks:

- uniform nodes lie in `[0, 1]`, with nodes and weights symmetric about `0.5`;
- exponential nodes are positive;
- standard-normal nodes are symmetric about zero, paired weights are equal, and
  odd reconstructed moments vanish within tolerance.

## Moment-Space Boundary Cases

`moment_space_boundary.cpp` covers the regime
[qmom_wheeler.md](../../../docs/developer/qmom_wheeler.md) identifies as most
sensitive: moment sequences approaching the boundary of realizable moment
space. The cases above are either comfortably interior, such as the classical
quadrature rules, or exactly degenerate; these cover the approach itself.

Every reference is constructed from known nodes and weights, so the expected
moments are defined independently of the inversion under test. Two families are
used:

| Family | Construction | Boundary approached by |
| --- | --- | --- |
| Node coalescence | three-node rule with the upper two nodes `separation` apart | `separation -> 0` |
| Weight collapse | three-node rule with a small middle weight | middle weight `-> 0` |

A third construction takes the six raw moments of an exact two-node rule. For a
three-node request that sequence lies exactly on the boundary: it is realizable,
but only by two atoms. Perturbing `m_4` alone then moves it strictly inside or
strictly outside while every other moment stays on the boundary.

The registered cases are:

- `coalescing_nodes_preserve_supplied_moments`: over a separation sweep, every
  accepted rule reproduces the moments its node count represents, carries
  non-negative weights summing to `m_0`, and returns sorted abscissas. A
  rejection reports `not_realizable`, because the input is structurally valid
  and realizable by construction.
- `adaptive_reduction_preserves_represented_moments`: the same sweep with
  adaptive inversion is accepted at every point, including where fixed-count
  inversion is entitled to reject.
- `accepted_reduction_respects_weight_ratio`: every accepted multi-node rule
  satisfies `min(w) / max(w) > rmin`, which is the documented adaptive
  acceptance criterion.
- `boundary_crossing_matches_realizability`: inside the boundary,
  `compute_qmom` succeeds and `check_realizability` reports strict
  realizability; outside it, both report non-realizability. Two independent
  implementations of the same mathematical fact must agree on which side a
  sequence lies.
- `fixed_count_does_not_reduce_on_the_boundary`: a fixed three-node request on
  the boundary must report `not_realizable` rather than quietly returning the
  two-node rule, while the same sequence is accepted once reduction is
  permitted.

These cases deliberately assert mathematical invariants and documented
contracts, not the separation or weight at which behavior changes. That
threshold depends on floating-point conditioning and is not a contract; pinning
it would make the suite platform-sensitive without testing anything QMOMKit
promises.

## Tolerance Policy

The current tolerance is `1e-12` for direct weights, abscissas, reconstructed
moments, mass, and distribution-specific symmetry checks. Relative error uses a
unit denominator floor so exact-zero references, such as Gaussian odd moments
and the center node of the three-node Hermite rule, are governed by the absolute
tolerance. Do not loosen this tolerance to mask a failure; first identify
whether the discrepancy comes from the input moments, reference convention, or
Wheeler inversion.

## Running

Build and run all QMOM analytical cases:

```bash
cd tests
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_PREFIX_PATH="../install"
cmake --build build
ctest --test-dir build -R verification_qmom_analytical --output-on-failure
```

Run one registered case:

```bash
ctest --test-dir build -R verification_qmom_analytical_uniform_n2 --output-on-failure
ctest --test-dir build -R verification_qmom_analytical_exponential_n3 --output-on-failure
ctest --test-dir build -R verification_qmom_analytical_gaussian_n3 --output-on-failure
ctest --test-dir build -R verification_qmom_compute_failures --output-on-failure
```

Run the executable manually from the `tests` directory:

```bash
./build/qmom_analytical_verification verification/qmom/uniform_n2.qmom uniform_n2
./build/qmom_analytical_verification verification/qmom/exponential_n2.qmom exponential_n2
./build/qmom_analytical_verification verification/qmom/gaussian_n2.qmom gaussian_n2
```

CTest passes absolute source paths, so all cases run correctly when CTest uses
the build directory as its working directory.
