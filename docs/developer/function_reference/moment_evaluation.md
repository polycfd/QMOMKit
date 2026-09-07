# Moment evaluation and instantaneous rates

The numerical definitions are in
[moment_reconstruction.cpp](../../../src/moment_reconstruction.cpp),
[advection_moment_rates.cpp](../../../src/advection_moment_rates.cpp), and
[moment_source_rates.cpp](../../../src/moment_source_rates.cpp). Public declarations
are in the corresponding [installed headers](../../../include/qmomkit).
See [numerical conventions](../numerical_conventions.md) and
[phase-space design](../phase_space_evolution_design.md). Public functions are in
`qmomkit`; the per-value arithmetic helpers are local, and trusted reconstruction
and advection-rate helpers are in `qmomkit::internal`.

## Raw moments — Level A

`reconstruct_raw_moments` computes
`M_alpha = sum_j W_j product_d x[d,j]^alpha[d]` from terminal quadrature values.
The QMOM overload uses orders zero through a requested maximum. CQMOM overloads
produce a dense tensor for requested extents or a selected vector of raw indices.
Selected requests retain order and duplicates; the algorithm does not silently
normalize by total weight.

Local `integer_power` uses exponentiation by squaring. Its loop consumes exponent
bits, multiplying the result for set bits and squaring the base only when another
bit remains. This avoids unnecessary final squaring, but does not provide scaled
or arbitrary-precision powers. `reconstruct_cqmom_moment_value` loops over nodes,
then coordinates, accumulates each weighted monomial, and checks finite products
and running sums. Summation is ordinary double arithmetic, without compensation.
An overflowing intermediate can fail even when a differently ordered or exact
sum would be finite.

Input validation rejects malformed/empty quadratures, nonfinite storage, and
rank mismatches. The QMOM result overload also requires success status. A public
CQMOM quadrature object does not carry its producing result's status, so its caller
must first check `CqmomResult`. Invalid input throws `invalid_argument`;
nonfinite reconstructed arithmetic throws `overflow_error`.

`reconstruct_cqmom_moments_impl` validates the required layout/hierarchy and
quadrature structure, then evaluates each required index. Its public wrapper
`reconstruct_cqmom_moments` checks canonical indices; internal
`reconstruct_trusted_cqmom_moments` reuses that previously established invariant.
Trusted mode still checks numerical storage. This function is central to accepted
state canonicalization after node motion or source inversion.

Evidence: [moment reconstruction tests](../../../tests/unit/moment_reconstruction.cpp)
use specified nodes/weights, selected duplicates, non-identity conditioning,
ragged output, invalid requests, and overflow cases. These establish exact finite
examples and failure branches, not cancellation or high-order accuracy limits.

## Generalized moments — Level A

`evaluate_generalized_moment` uses the same terminal sum with real exponents.
Local `generalized_power` rejects a negative abscissa with a noninteger exponent;
`is_integer_exponent` means exact `trunc(exponent) == exponent`, with no rounding
tolerance. Otherwise it calls `pow` and rejects nonfinite output. Thus zero raised
to a negative exponent fails, and a zero-weight node does not bypass an invalid
power. Every exponent must be finite; standalone generalized evaluation permits
negative exponents when the node arithmetic is defined and finite.

The CQMOM coordinate loop validates exponent rank, evaluates each factor, then
checks product and accumulation. No support restriction is inferred from coordinate
names. A finite terminal approximation neither establishes integrability of a
continuous singular integrand nor proves an underlying distribution is supported
on the positive half-line. Source-context generalized access imposes the additional
nonnegative-exponent restriction described in [sources](moment_sources.md).

Evidence: fractional/integer negative-node, zero-node, malformed exponent, and
overflow assertions in [moment reconstruction tests](../../../tests/unit/moment_reconstruction.cpp).
No operating envelope for exponents approaching an integer from either side,
near-zero nodes with negative powers, or signed cancellation was identified.

## Chain-rule advection rates — Level A

Local `reconstruct_advection_moment_rate_value` and internal
`reconstruct_advection_moment_rates` evaluate

`dM_alpha/dt = sum_j W_j sum_d alpha[d] v[d,j] x[d,j]^(alpha[d]-1) product_(k!=d) x[k,j]^alpha[k]`.

The loops respectively visit requested moments, terminal nodes, derivative
coordinates, and monomial coordinates. Terms with zero `alpha[d]` are skipped
before decrementing the exponent. The implementation uses powers rather than
`M/x`, so zero abscissas cause no spurious division by zero. Zeroth-moment
advection rate is zero because its derivative-coordinate sum is empty. Weights
are held fixed; this generic identity does not imply conservation of nonzero-order
moments under arbitrary caller velocity.

`validate_advection_moment_rate_storage` checks the quadrature, velocity count,
coordinate rank, and finite values. Every term/product/sum is checked; overflow
throws. Complexity grows with requested moment count, terminal nodes, and rank;
repeated powers are recomputed rather than globally cached.

Public `evaluate_qmom_advection_moment_rates` and
`evaluate_cqmom_advection_moment_rates` reconstruct a quadrature, evaluate terminal
velocities, and then apply that identity. Internal terminal-velocity evaluation
constructs a coordinate vector for each node and invokes the caller at the supplied
time. Wrong size/nonfinite velocity or a caught callback exception becomes a
velocity failure; reconstruction failure is separately reported. Final arithmetic
overflow can propagate as an exception rather than a returned velocity status.
These functions produce an instantaneous rate; they do not advance the state.

Evidence: [moment reconstruction tests](../../../tests/unit/moment_reconstruction.cpp)
compare advection derivatives with independently perturbed coordinates;
[advection rate tests](../../../tests/unit/advection_moment_rates.cpp) exercise
reconstruction, callback failures, and dense/lean interfaces. Fixed finite-difference
examples do not measure loss of accuracy at very high order or cancellation.

## Source-rate wrappers — Level B

`evaluate_qmom_moment_source_rates` and `evaluate_cqmom_moment_source_rates` create
borrowed source contexts and validate the returned rate shape/finiteness. They
report invalid request, callback failure, or success without integration.
Where a boundary callback overload is provided, source evaluation precedes
boundary evaluation; a failed source prevents calling the boundary. The two rate
vectors are added and the combined result checked. The boundary contribution is
an additive moment rate, not a spatial boundary condition.

Standalone contexts do not receive reconstruction metadata for untracked moments.
QMOM's instantaneous source wrapper checks nonempty state, finite time, and callback
presence but does not itself reject every nonfinite input moment; the CQMOM route
checks finite values. A callback that ignores a bad QMOM input can therefore return
a finite rate. Integration has stronger initial-state validation. This is an API
boundary distinction, not a claim that source evaluation repairs invalid moments.

Evidence: [source-rate tests](../../../tests/unit/moment_source_rates.cpp) exercise
shape, callback, boundary-addition, and error propagation; no targeted ignored-NaN
QMOM input assertion was identified.

## Local validators and adapters — Level C

| Family | Purpose / semantics | Evidence |
|---|---|---|
| `validate_reconstructed_value`, `validate_advection_moment_rate_value`, `validate_generalized_value`, `validate_generalized_exponent` | Distinguish malformed exponent input from nonfinite arithmetic, using operation-specific exception text. | Domain/overflow unit cases |
| Quadrature validation message builders | Supply context-specific error wording to shared structural checks; do not change the checks. | Invalid-input tests indirectly |
| Rate invalid/failure result factories and `validate_*` helpers | Preserve requested output shape/status and validate context/rate structure before success. | Rate API failure tests |
| Dense/lean and scalar/multidimensional forwarding overloads | Adapt the caller representation to the same mathematical sum/context operation. | Dense/lean parity and public-header coverage |

Exact overload names and definitions are listed in the [inventory](coverage_inventory.md).
