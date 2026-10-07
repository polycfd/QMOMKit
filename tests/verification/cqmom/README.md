# Two-Dimensional CQMOM Analytical Verification

These cases verify the public CQMOM path from configuration loading through
`qmomkit::compute_cqmom(...)`. Inputs are loaded from `.qmom` files in this
directory through `qmomkit::config::load_cqmom` and
`qmomkit::config::load_cqmom_prepared`; the executables do not overwrite or
regenerate configured raw moments. The recursive-inversion, bounded and
branch-growth studies described under
[Recursive Inversion Studies](#recursive-inversion-studies) build their moments
in code instead, as do the disparate-units automatic-scaling case and the
one-dimensional all-orders reconstruction check.

## Raw Moment Reconstruction

`moment_reconstruction.cpp` verifies public moment reconstruction APIs on
independently known quadratures:

```text
reconstruct_raw_moments(Quadrature1D, moment_count)
reconstruct_raw_moments(CqmomQuadrature, RawMomentRequest)
reconstruct_cqmom_moments(CqmomQuadrature, CqmomMomentLayout)
```

Selected reconstruction is intentionally not a view into the full tensor. A
`RawMomentRequest` is an ordered list, so returned values must preserve
non-monotonic request order and duplicates. The verifier requests:

```text
{2,0}, {0,0}, {1,1}, {0,2}, {2,0}, {3,1}
```

and checks the values against both an independent analytical product quadrature
and the corresponding entries of a full-tensor reconstruction. This makes the
ordering rule explicit: the repeated `{2,0}` entry must appear twice with the
same value, and no sorting or deduplication may occur.

The `reconstruct_cqmom_moments(...)` check uses the `CqmomMomentLayout` of the
case loaded with `load_cqmom(...)` and verifies values align with
`layout.required_indices`. Failure-path checks cover an unsuccessful 1D
quadrature, a selected CQMOM request with the wrong rank, and a CQMOM layout
whose hierarchy does not match the quadrature. These public APIs report invalid
reconstruction input by throwing exceptions with explanatory messages rather
than by returning a status object.

## Public Workflow Consistency

`cqmom_public_workflows.cpp` verifies the recommended user-facing CQMOM loading
paths on the independent-product case:

```text
load_cqmom(...) -> CqmomCase -> compute_cqmom(...)
load_cqmom_prepared(...)      -> CqmomPreparedCase     -> compute_cqmom(...)
```

The default `CqmomPreparedCase` workflow is the normal user path: it builds the
recursive branch structure of its hierarchy once and shares it with every
reconstruction and state. An open `CqmomCase` suits callers that build or edit
a case in code; its own `compute_cqmom(...)` rebuilds that structure on every
call.

Both workflows are expected to recover the same analytical terminal
quadrature and reconstruct the same raw moments for the same mathematical
problem. A failure indicates a mismatch between the two case types, lost
configuration policy, or an unintended difference in diagnostics rather than a
new recursive quadrature formula.

The `CqmomCase` checks also verify the moment data the loader stores:

- the required multi-index order is the documented fixed 2x2 CQMOM order
  `{0,0}`, `{1,0}`, `{2,0}`, `{3,0}`, `{0,1}`, `{0,2}`, `{0,3}`, `{1,1}`,
  `{1,2}`, `{1,3}`;
- `CqmomMomentData::values` has one value per required index;
- the case stores fewer moments than the full 4x4 tensor;
- every stored value matches the independent-product analytical moment.

These checks document that `compute_cqmom(CqmomCase)` needs only the required
moments, not the full tensor, and must produce the same quadrature as the
default workflow.

## Public Compute Failure Contracts

Direct public `compute_cqmom(...)` failure/status regression coverage lives in
`tests/component_regression/cqmom_compute_failures.cpp`. It checks invalid
public requests and status mapping across `CqmomPreparedCase` and `CqmomCase`.
This verification directory keeps the analytical and public workflow success
cases.

The non-realizable case uses a finite two-dimensional moment tensor with
`M[0,0] < 0`. Configuration parsing succeeds because the input is structurally
valid, but CQMOM reconstruction cannot succeed because the root one-dimensional
moment sequence has negative total number density. The verifier checks
consistent behavior for:

```text
compute_cqmom(CqmomPreparedCase)
compute_cqmom(CqmomCase)
```

Both workflows must fail with `CqmomStatus::qmom_failed`, preserve an
explanatory message, and expose public diagnostics identifying the root
operation and non-success QMOM status.

The same verifier also covers structural mistakes in a full raw moment tensor,
which `extract_cqmom_moments(...)` rejects by throwing before a case exists:

- tensor extent rank mismatch (`std::invalid_argument`);
- tensor value-count mismatch (`std::invalid_argument`);
- extents that omit a required CQMOM moment (`std::out_of_range`);

and invalid public compute inputs that should be rejected before successful
quadrature is exposed:

- `CqmomCase` moment data with missing required values;
- `CqmomCase` moment data with non-finite values;
- `CqmomCase` moment data whose required-index layout does not match its
  hierarchy;
- adaptive Wheeler reductions requested without enabling
  `adaptive_branch_counts`.

The compute failures are expected to return `CqmomStatus::invalid_input` with
useful messages. Diagnostic wording is not required to be identical across the
two case types, but the user-facing behavior must be consistent: `ok()` is
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
results expose terminal arrays and a `ConditioningPathTable` in deterministic
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
- reconstructs the configured required CQMOM moments with
  `qmomkit::reconstruct_cqmom_moments(...)` and compares them using the
  unit-floored relative error.

Direct quadrature and path checks are required in addition to moment
reconstruction because the branch structures are analytically known.

The public workflow executable `cqmom_public_workflows_verification` reuses the
independent-product reference. It checks successful computation, terminal node
and weight accuracy, raw-moment reconstruction, branch diagnostics, and
agreement between the `CqmomCase` and `CqmomPreparedCase` workflows.

## Tolerance Policy

The comparison tolerance is `1e-12` for nodes, weights, and reconstructed raw
moments. Relative moment errors use a unit denominator floor so exact or
near-zero moments are governed by absolute error rather than unstable relative
ratios.

Some high-order moments in these cases are hundreds to thousands in magnitude.
An absolute reconstruction error near `1e-11` can therefore correspond to a
unit-floored relative error near machine precision, about `1e-15`. On failure
the verification reports both values so scale-dependent roundoff is visible without
weakening near-zero checks. Do not loosen tolerances without first identifying
whether the discrepancy comes from the configuration, analytical reference,
conditioning solve, Wheeler inversion, or raw-moment reconstruction.

## Running

Configure and build the test suite from the repository root, in Debug as
described in [Testing](../../../docs/developer/testing.md#debug-build-requirement):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DQMOMKIT_BUILD_TESTS=ON
cmake --build build --parallel
```

Run the three two-dimensional analytical cases:

```bash
ctest --test-dir build -R '^verification_cqmom_(independent_product|deterministic_affine|branch_dependent)$' --output-on-failure
```

Run individual cases:

```bash
ctest --test-dir build -R '^verification_cqmom_public_workflows$' --output-on-failure
ctest --test-dir build -R '^verification_cqmom_moment_reconstruction_cqmom_selected$' --output-on-failure
ctest --test-dir build -R '^verification_cqmom_independent_product$' --output-on-failure
ctest --test-dir build -R '^verification_cqmom_deterministic_affine$' --output-on-failure
ctest --test-dir build -R '^verification_cqmom_branch_dependent$' --output-on-failure
```

Run the executables directly from the repository root after building:

```bash
./build/tests/verification/cqmom/cqmom_public_workflows_verification tests/verification/cqmom/independent_product.qmom
./build/tests/verification/cqmom/moment_reconstruction_verification tests/verification/cqmom/independent_product.qmom
./build/tests/verification/cqmom/cqmom_2d_analytical_verification tests/verification/cqmom/independent_product.qmom independent_product
./build/tests/verification/cqmom/cqmom_2d_analytical_verification tests/verification/cqmom/deterministic_affine.qmom deterministic_affine
./build/tests/verification/cqmom/cqmom_2d_analytical_verification tests/verification/cqmom/branch_dependent.qmom branch_dependent
```

## Recursive Inversion Studies

Four more executables in this directory verify the recursive inversion on
quadratures written down in advance. They build their moments in code with a
local helper that sums weight times position to the required power, so the
reference does not depend on the library's reconstruction.

- `recursive_quadrature.cpp` (`verification_recursive_quadrature_manufactured_*`):
  the recursion recovers the prescribed quadrature in one to three dimensions,
  for identity and non-identity conditioning orders, non-uniform node counts,
  and adaptive branches that return fewer nodes than requested. The
  one-dimensional case reduces to plain Wheeler. Rebuilt moments are compared
  with per-case tolerances from `1e-10` to `1e-7`, and branch node counts and
  diagnostic record counts must match exactly.
- `recursive_cqmom.cpp` (`verification_cqmom_recursive_*`): six fixed cases,
  from one to four dimensions with identity and non-identity conditioning
  orders and non-uniform node counts, must recover the prescribed weights and
  positions; the five dimension sweeps `dimension_sweep_D2` to `D6` must
  reproduce the moments they were given, without a node comparison; and one
  adaptive 4D case must return the prescribed ragged structure, with the right
  node count on each branch, terminal paths, and local and global weights.
  Values use per-case tolerances, `1e-12` for the sweeps.
- `quadrature_bounds.cpp` (`verification_cqmom_quadrature_bounds_*`): a branch
  whose positions lie inside a coordinate's bounds is accepted unchanged, one
  with a position outside is rejected rather than clipped, and a bound applies
  to its phase-space coordinate whatever its conditioning level. Accepted values
  agree to `1e-12`.
- `adaptive_branch_growth.cpp` (`verification_cqmom_adaptive_branch_growth_*`):
  during a moment-source step, a branch that starts with one conditional node
  grows to two when the moments move to a state that needs two. Terminal paths
  and per-branch node counts must match the target exactly, the diagnostics
  must record the growth, and values agree to `5e-11`.

Run them, from the same build as above, with:

```bash
ctest --test-dir build -R '^verification_(recursive_quadrature_manufactured_|cqmom_(recursive_|quadrature_bounds_|adaptive_branch_growth_))' --output-on-failure
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
- checks the configured required moments against independently generated
  analytical raw moments without modifying the request;
- calls public `qmomkit::compute_cqmom(...)`;
- compares terminal coordinates, global weights, local branch weights, branch
  counts, and paths where adaptive paths are uniquely defined;
- verifies independent-product conditional rules are identical across parent
  branches for each phase-space dimension;
- verifies affine relations directly;
- verifies positive weights and `sum(W) = M[0,...,0]`;
- reconstructs the full tensor of the configured extents with
  `qmomkit::reconstruct_raw_moments(...)` and compares it with the analytical
  tensor.

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

## Automatic Scaling CQMOM Verification

`automatic_scaling.cpp` verifies automatic reconstruction scaling through real
configuration loading for a prescribed 3D ragged terminal rule.  The hierarchy
requests `{2,2,2}` nodes in phase-space order but recurses in the non-identity
order `{2,0,1}`.  Its five paths are `{0,0,0}`, `{0,0,1}`, `{1,0,0}`,
`{1,1,0}`, and `{1,1,1}`: the level-one branches have active counts one and
two, and the final branches have counts two, one, and two.

The fixtures `automatic_scaling_reference_3d.qmom` and
`automatic_scaling_rescaled_3d.qmom` prescribe the same terminal rule in
physical coordinate units `(1,1,1)` and `(1e6,1e-6,1e3)`. Every complete
`4 x 4 x 4` tensor is independently summed from the terminal rule. The
physical rescaled tensor has `M'_alpha = (s0^a0 s1^a1 s2^a2) M_alpha`; the
automatic-scaling normalization checked by this test is therefore

```text
Mhat_a0,a1,a2 = M'_a0,a1,a2 / (s0^a0 s1^a1 s2^a2)
```

when read back in canonical units, irrespective of the recursive order. The
verifier loads the configured case, confirms the configured
automatic-scaling/adaptive policy, reconstructs it, and checks paths, active
branch counts, global weights, and physical coordinates against the
terminal-rule oracle, and the 22 required CQMOM moments against the configured
moments, which are summed independently from the terminal rule. This is
scaling-enabled interaction and covariance evidence: it covers adaptive
topology, non-identity conditioning, configuration propagation, and
physical-unit restoration. It does not claim scaling is necessary for the
`(1e6,1e-6,1e3)` representation.

The tolerance is `2e-10` for unit-floored relative node and moment errors and
absolute weights. It leaves cross-platform margin for nested adaptive
inversions; observed errors are near machine precision.

`verification_cqmom_automatic_scaling_disparate_units` uses the same
independently prescribed five-terminal rule through the direct public CQMOM
route, but represents its coordinates in units `(1e8,1e-8,1e4)`. The
complete `4 x 4 x 4` raw tensor is explicitly summed from those physical
nodes and positive weights, so the terminal rule remains the oracle. The
physical factors create mixed moments from approximately `1e-24` to `1e36`;
automatic scaling computes each coordinate scale from its marginal second
moment and returns the normalized problem to O(1) coordinates.

With scaling on, the verifier requires the five terminals, their paths,
physical nodes, weights, and the required CQMOM moments to agree with the oracle
to `1e-10` (unit-floored relative errors for nodes and moments, absolute error
for weights). The per-branch active counts are not checked in this case. With
scaling off, the same physical moment tensor must recover the same rule to the
same tolerances: adaptive node reduction compares each recurrence entry with a
roundoff bound in the same units, so it does not depend on the units of the
coordinates.

## High-Dimensional Running

Run all high-dimensional cases from the same build:

```bash
ctest --test-dir build -R '^verification_cqmom_(affine_chain_5d|independent_product_[3456]d|mixed_affine_independent_4d)$' --output-on-failure
```

Run individual cases:

```bash
ctest --test-dir build -R '^verification_cqmom_affine_chain_5d$' --output-on-failure
ctest --test-dir build -R '^verification_cqmom_independent_product_3d$' --output-on-failure
ctest --test-dir build -R '^verification_cqmom_independent_product_4d$' --output-on-failure
ctest --test-dir build -R '^verification_cqmom_independent_product_5d$' --output-on-failure
ctest --test-dir build -R '^verification_cqmom_independent_product_6d$' --output-on-failure
ctest --test-dir build -R '^verification_cqmom_mixed_affine_independent_4d$' --output-on-failure
ctest --test-dir build -R '^verification_cqmom_automatic_scaling_' --output-on-failure
ctest --test-dir build -R '^verification_cqmom_automatic_scaling_disparate_units$' --output-on-failure
```
