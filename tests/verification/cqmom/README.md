# Two-Dimensional CQMOM Analytical Verification

These cases verify the public CQMOM path from configuration loading through
`qmomkit::compute_cqmom(...)`. Inputs are loaded from `.qmom` files in this
directory through `qmomkit::config::load_cqmom_dense`; the executables do
not overwrite or regenerate configured raw moments.

## Raw Moment Reconstruction

`moment_reconstruction.cpp` verifies public moment reconstruction APIs on
independently known quadratures:

```text
reconstruct_raw_moments(Quadrature1D, moment_count)
reconstruct_raw_moments(CqmomQuadrature, RawMomentRequest)
reconstruct_cqmom_moments(CqmomQuadrature, CqmomMomentLayout)
```

Selected reconstruction is intentionally not a dense tensor view. A
`RawMomentRequest` is an ordered list, so returned values must preserve
non-monotonic request order and duplicates. The verifier requests:

```text
{2,0}, {0,0}, {1,1}, {0,2}, {2,0}, {3,1}
```

and checks the values against both an independent analytical product
quadrature and the corresponding entries in a dense reconstruction. This makes
the ordering rule explicit: the repeated `{2,0}` entry must appear twice with
the same value, and no sorting or deduplication may occur.

The `reconstruct_cqmom_moments(...)` check uses the public lean
`CqmomMomentLayout` from `load_cqmom_lean(...)` and verifies values align with
`layout.required_indices`. Failure-path checks cover an unsuccessful 1D
quadrature, a selected CQMOM request with the wrong rank, and a CQMOM layout
whose hierarchy does not match the quadrature. These public APIs report invalid
reconstruction input by throwing exceptions with explanatory messages rather
than by returning a status object.

## Public Workflow Consistency

`cqmom_public_workflows.cpp` verifies the recommended user-facing CQMOM loading
paths on the independent-product case:

```text
load_cqmom_dense(...) -> CqmomDenseCase -> compute_cqmom(...)
load_cqmom_lean(...)  -> CqmomLeanCase  -> compute_cqmom(...)
load_cqmom(...)       -> CqmomCase      -> compute_cqmom(...)
```

The default `CqmomCase` workflow is the normal user path because it stores the
configured CQMOM setup and avoids exposing dense or lean moment storage. The
dense workflow is useful when callers intentionally provide a full rectangular
raw-moment tensor. The lean workflow is useful when callers intentionally
manage only the minimal CQMOM moment layout required by recursive
reconstruction.

All three workflows are expected to recover the same analytical terminal
quadrature and reconstruct the same raw moments for the same mathematical
problem. A failure indicates a public representation mismatch, lost
configuration policy, or an unintended difference in diagnostics rather than a
new recursive quadrature formula.

The explicit lean checks also verify the expert representation itself:

- the required multi-index order is the documented fixed 2x2 CQMOM order
  `{0,0}`, `{1,0}`, `{2,0}`, `{3,0}`, `{0,1}`, `{0,2}`, `{0,3}`, `{1,1}`,
  `{1,2}`, `{1,3}`;
- `CqmomMomentData::values` has one value per required index;
- the lean representation stores fewer selected moments than the dense 4x4
  tensor;
- every selected value matches the independent-product analytical moment.

These checks document that the public lean workflow does not require a dense
tensor at compute time: `compute_cqmom(CqmomLeanCase)` receives the selected
required moment representation and must produce the same mathematical
quadrature as the dense and default workflows.

## Public Compute Failure Contracts

Direct public `compute_cqmom(...)` failure/status regression coverage lives in
`tests/component_regression/cqmom_compute_failures.cpp`. It checks invalid
public requests and status mapping across default, dense, and lean public case
objects. This verification directory keeps the analytical and public workflow
success cases.

The non-realizable case uses a finite two-dimensional moment tensor with
`M[0,0] < 0`. Configuration parsing succeeds because the input is structurally
valid, but CQMOM reconstruction cannot succeed because the root one-dimensional
moment sequence has negative total number density. The verifier checks
consistent behavior for:

```text
compute_cqmom(CqmomCase)
compute_cqmom(CqmomDenseCase)
compute_cqmom(CqmomLeanCase)
```

All three workflows must fail with `CqmomStatus::qmom_failed`, preserve an
explanatory message, and expose public diagnostics identifying the root
operation and non-success QMOM status.

The same verifier also covers invalid public compute inputs that should be
rejected before successful quadrature is exposed:

- dense tensor extent rank mismatch;
- dense tensor value-count mismatch;
- non-finite dense raw moments;
- dense extents that omit required CQMOM moments;
- lean moment data with missing required values;
- lean moment data with non-finite values;
- lean moment data whose required-index layout does not match its hierarchy;
- adaptive Wheeler reductions requested without enabling
  `adaptive_branch_counts`.

These failures are expected to return `CqmomStatus::invalid_input` with useful
messages. Diagnostic wording is not required to be identical across public
representations, but the user-facing behavior must be consistent: `ok()` is
false, status is non-success, and no invalid quadrature is reported as
successful.

## Case A: Independent Product

The parent variable is:

```text
x = {1, 3}
w_x = {0.4, 0.6}
```

The conditional variable is independent of `x`:

```text
y = {-1, 2}
w_y = {0.25, 0.75}
```

The raw moments are products:

```text
M[i,j] = m_x[i] * m_y[j]
```

with fixed `node_counts = {2, 2}`, `conditioning_order = {0, 1}`,
`qmom.adaptive = false`, `adaptive_branch_counts = false`, and
`extents = {4, 4}`. Every parent branch must recover the same conditional
`y` rule. Terminal/global weights are products `w_x * w_y`.

## Case B: Deterministic Affine Correlation

The parent distribution is the same as Case A, and:

```text
y = 2*x + 1
M[i,j] = sum_q w_q * x_q^i * (2*x_q + 1)^j
```

Thus:

```text
x = 1 -> y = 3
x = 3 -> y = 7
```

Each branch has exactly one conditional node with local conditional weight
`1`. The configured hierarchy still requests two `y` nodes, but
`adaptive_branch_counts = true`, so CQMOM automatically uses adaptive Wheeler
and reduces each conditional branch from two requested nodes to one actual
node. This reduction is expected success, not failure. The fixture also stores
`qmom.adaptive = true`, but standalone QMOM remains independently configurable.

## Case C: Branch-Dependent Conditional Rules

The parent distribution is:

```text
x = {1, 3}
w_x = {0.35, 0.65}
```

Conditional on `x = 1`:

```text
y = {-0.5}
w_y = {1}
```

Conditional on `x = 3`:

```text
y = {1, 4}
w_y = {0.3, 0.7}
```

The raw moments are generated from the terminal/global quadrature:

```text
M[i,j] = sum_terminal W * x^i * y^j
```

The first branch adaptively reduces to one child, while the second keeps two
children. The result is a ragged terminal tree with paths `{0,0}`, `{1,0}`,
and `{1,1}`.

## Weight Semantics

Local conditional weights are branch-local Wheeler weights. Accumulated or
terminal/global weights multiply all local weights along a complete
conditioning path:

```text
W_global = w_x * w_y|x
```

Fixed results use rectangular phase-space-order storage. Adaptive branch-count
results expose terminal arrays and `ConditioningPath` values in deterministic
traversal order; abscissas are still stored by phase-space dimension.

## Adaptive Controls

`qmom.adaptive` controls standalone one-dimensional QMOM. For CQMOM,
`adaptive_branch_counts` controls whether recursive CQMOM accepts branch-local
node counts smaller than the configured maximum hierarchy. When adaptive branch
counts are enabled, recursive CQMOM automatically uses adaptive Wheeler for each
marginal and conditional inversion. Recursive CQMOM still rejects adaptive
Wheeler reductions in fixed rectangular mode.

## Checks

The reusable executable `cqmom_2d_analytical_verification`:

- loads every input from the corresponding `.qmom`;
- verifies configured node counts, conditioning order, `qmom.adaptive`,
  `adaptive_branch_counts`, and extents;
- calls public `qmomkit::compute_cqmom(...)`;
- compares parent nodes and parent weights;
- compares local conditional nodes and weights;
- compares terminal/global weights and phase-space coordinates;
- checks branch-specific terminal paths and requested/actual diagnostics;
- verifies all terminal weights are positive and total weight equals `M[0,0]`;
- reconstructs the complete configured raw-moment tensor with
  `qmomkit::reconstruct_raw_moments(...)`.

Direct quadrature and path checks are required in addition to moment
reconstruction because the branch structures are analytically known.

The public workflow executable `cqmom_public_workflows_verification` reuses the
independent-product reference. It checks successful computation, terminal node
and weight accuracy, raw-moment reconstruction, branch diagnostics, and
dense/default/lean equivalence. The expert lean workflow is covered as part of
`verification_cqmom_public_workflows`.

## Tolerance Policy

The comparison tolerance is `1e-12` for nodes, weights, and reconstructed raw
moments. Relative moment errors use a unit denominator floor so exact or
near-zero moments are governed by absolute error rather than unstable relative
ratios.

Some high-order moments in these cases are hundreds to thousands in magnitude.
An absolute reconstruction error near `1e-11` can therefore correspond to a
unit-floored relative error near machine precision, about `1e-15`. The
verification reports both values so scale-dependent roundoff is visible without
weakening near-zero checks. Do not loosen tolerances without first identifying
whether the discrepancy comes from the configuration, analytical reference,
conditioning solve, Wheeler inversion, or raw-moment reconstruction.

## Running

Configure and build the standalone test project first:

```bash
cd tests
cmake -S . -B build \
    -DCMAKE_PREFIX_PATH="../install"
cmake --build build
```

Run all three CQMOM analytical cases:

```bash
ctest --test-dir build -R verification_cqmom_ --output-on-failure
```

Run individual cases:

```bash
ctest --test-dir build -R verification_cqmom_public_workflows --output-on-failure
ctest --test-dir build -R verification_cqmom_moment_reconstruction_selected --output-on-failure
ctest --test-dir build -R verification_cqmom_independent_product --output-on-failure
ctest --test-dir build -R verification_cqmom_deterministic_affine --output-on-failure
ctest --test-dir build -R verification_cqmom_branch_dependent --output-on-failure
```

Run manually from the standalone test project after building `tests/build`:

```bash
./build/cqmom_public_workflows_verification verification/cqmom/independent_product.qmom
./build/moment_reconstruction_verification verification/cqmom/independent_product.qmom
./build/cqmom_2d_analytical_verification verification/cqmom/independent_product.qmom independent_product
./build/cqmom_2d_analytical_verification verification/cqmom/deterministic_affine.qmom deterministic_affine
./build/cqmom_2d_analytical_verification verification/cqmom/branch_dependent.qmom branch_dependent
```

# Higher-Dimensional CQMOM Analytical Verification

The high-dimensional executable
`cqmom_high_dimensional_analytical` extends the same configuration-driven
policy to deeper recursive traversals and non-identity conditioning orders.
The configured raw moments are generated analytically from exact terminal
quadratures:

```text
M[k_1,...,k_D] = sum_q W_q product_d xi_{q,d}^{k_d}
```

The C++ verifier keeps independent analytical reference data only: terminal
coordinates, global weights, local branch weights, conditioning paths, branch
counts, and affine coefficients.

## Affine Chain 5D

The 5D affine chain starts from:

```text
xi_1 = {-0.5, 1.5}, weights = {0.4, 0.6}
xi_2 = 2 xi_1 + 1
xi_3 = -0.5 xi_2 + 0.75
xi_4 = 1.25 xi_3 - 0.2
xi_5 = -0.8 xi_4 + 0.5
```

It uses `conditioning_order = {0,1,2,3,4}` and
`adaptive_branch_counts = true`. The first level has two nodes; CQMOM
automatically uses adaptive Wheeler, so every later branch reduces to one
deterministic child with local weight `1`. The terminal weights remain
`{0.4, 0.6}` and the terminal coordinates are stored in phase-space dimension
order. The fixture also records `qmom.adaptive = true`; that option controls
standalone QMOM behavior independently.

This case validates recursion depth and repeated adaptive reduction. Its
configured tensor has `4^5 = 1024` moments and two terminal points.

## Independent Products D=3..6

The independent product sweep uses two-node one-dimensional rules in every
dimension. The first six rules are:

```text
dim 0: nodes {-1.0, 1.5}, weights {0.35, 0.65}
dim 1: nodes { 0.2, 1.2}, weights {0.60, 0.40}
dim 2: nodes {-0.5, 0.75}, weights {0.45, 0.55}
dim 3: nodes { 2.0, 3.5}, weights {0.70, 0.30}
dim 4: nodes {-2.0, 0.5}, weights {0.25, 0.75}
dim 5: nodes { 1.0, 2.25}, weights {0.80, 0.20}
```

For each dimension count, terminal weights are products of the selected
one-dimensional weights. Fixed branch counts are used:
`qmom.adaptive = false` and `adaptive_branch_counts = false`.

The non-identity conditioning orders are:

```text
D=3: {1,0,2}
D=4: {2,0,3,1}
D=5: {4,1,3,0,2}
D=6: {5,2,0,4,1,3}
```

These cases validate conditioning-order handling, phase-space-order terminal
coordinates, identical conditional rules across independent branches, positive
weights, total mass, and full raw-moment reconstruction. The computational
sizes are `4^3 = 64`, `4^4 = 256`, `4^5 = 1024`, and `4^6 = 4096`
configured moments with `2^D` terminal points.

## Mixed Affine And Independent 4D

The mixed case uses:

```text
xi_1 = {-0.5, 1.0}, weights = {0.45, 0.55}
xi_2 = 2 xi_1 + 1
xi_3 = {-1.0, 2.0}, weights = {0.30, 0.70}
xi_4 = -xi_2 + 0.5 xi_3
```

with `conditioning_order = {0,2,1,3}` and
`adaptive_branch_counts = true`. CQMOM automatically uses adaptive Wheeler for
the marginal and conditional inversions. The exact terminal quadrature has four
points: each `xi_1` branch paired with each independent `xi_3` node. The first
two levels have two children, and the affine descendants reduce to one child
per parent path. Local weights are `{w_xi1, w_xi3, 1, 1}` along a complete
path; accumulated terminal weights are `w_xi1 * w_xi3`. The fixture also
records `qmom.adaptive = true`; that option remains a standalone QMOM setting.

This case validates coupled conditional reconstruction with a non-identity
conditioning order and a ragged adaptive tree. It has `4^4 = 256` configured
moments and four terminal points.

## High-Dimensional Checks

The high-dimensional executable:

- loads every input from the corresponding `.qmom`;
- verifies configured node counts, conditioning order, `qmom.adaptive`,
  `adaptive_branch_counts`, and extents;
- checks configured moments against independently generated analytical raw
  moments without modifying the request;
- calls public `qmomkit::compute_cqmom(...)`;
- compares terminal coordinates, global weights, local branch weights, branch
  counts, and paths where adaptive paths are uniquely defined;
- verifies independent-product conditional rules are identical across parent
  branches for each phase-space dimension;
- verifies affine relations directly;
- verifies positive weights and `sum(W) = M[0,...,0]`;
- reconstructs every configured raw moment using
  `qmomkit::reconstruct_raw_moments(...)`.

Adaptive cases compare exact `ConditioningPath` values. Fixed independent
cases compare unordered physical terminal quadrature sets by nearest
coordinate matching, because rectangular storage is deterministic but tiny
roundoff can otherwise perturb lexicographic sorting of equal coordinates.

The high-dimensional tolerance is `1e-11` for nodes, local weights, global
weights, and unit-floored relative raw-moment errors. Absolute errors are used
for exact or near-zero quantities. The tolerance is looser than the 2D cases
because the 6D fixed product case exercises several nested Vandermonde solves
and 4096 moment entries; observed relative moment errors remain near machine
precision.

## High-Dimensional Running

Run all high-dimensional cases:

```bash
ctest --test-dir build -R 'verification_cqmom_(affine_chain_5d|independent_product_[3456]d|mixed_affine_independent_4d)' --output-on-failure
```

Run individual cases:

```bash
ctest --test-dir build -R verification_cqmom_affine_chain_5d --output-on-failure
ctest --test-dir build -R verification_cqmom_independent_product_3d --output-on-failure
ctest --test-dir build -R verification_cqmom_independent_product_4d --output-on-failure
ctest --test-dir build -R verification_cqmom_independent_product_5d --output-on-failure
ctest --test-dir build -R verification_cqmom_independent_product_6d --output-on-failure
ctest --test-dir build -R verification_cqmom_mixed_affine_independent_4d --output-on-failure
```
