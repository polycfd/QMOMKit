# Moment evaluation

The numerical definitions are in
[moment_reconstruction.cpp](../../../src/moment_reconstruction.cpp). Public
declarations are in the corresponding [installed
headers](../../../include/qmomkit). See [numerical
conventions](../numerical_conventions.md) and [phase-space
design](../phase_space_evolution_design.md). Public functions are in `qmomkit`;
the per-value arithmetic helpers are local, and trusted reconstruction helpers
are in `qmomkit::internal`.

## Raw moments — Level A

`reconstruct_raw_moments` computes `M_alpha = sum_j W_j product_d
x[d,j]^alpha[d]` from terminal quadrature values. The QMOM overload uses orders
zero through a requested maximum. CQMOM overloads produce a dense tensor for
requested extents or a selected vector of raw indices. Selected requests retain
order and duplicates; the algorithm does not silently normalize by total
weight.

Local `integer_power` uses exponentiation by squaring. Its loop consumes
exponent bits, multiplying the result for set bits and squaring the base only
when another bit remains. This avoids unnecessary final squaring, but does not
provide scaled or arbitrary-precision powers. `reconstruct_cqmom_moment_value`
loops over nodes, then coordinates, accumulates each weighted monomial, and
checks finite products and running sums. Summation is ordinary double
arithmetic, without compensation. An overflowing intermediate can fail even
when a differently ordered or exact sum would be finite.

Input validation rejects malformed/empty quadratures, nonfinite storage, and
rank mismatches. The QMOM result overload also requires success status. A
public CQMOM quadrature object does not carry its producing result's status, so
its caller must first check `CqmomResult`. Invalid input throws
`invalid_argument`; nonfinite reconstructed arithmetic throws `overflow_error`.

`reconstruct_cqmom_moments_impl` validates the required layout/hierarchy and
quadrature structure, then evaluates each required index. Its public wrapper
`reconstruct_cqmom_moments` checks canonical indices; internal
`reconstruct_trusted_cqmom_moments` reuses that previously established
invariant. Trusted mode still checks numerical storage. This function is
central to accepted state canonicalization after node motion or source
inversion.

Evidence: [moment reconstruction
tests](../../../tests/unit/moment_reconstruction.cpp) use specified
nodes/weights, selected duplicates, non-identity conditioning, ragged output,
invalid requests, and overflow cases. These establish exact finite examples and
failure branches, not cancellation or high-order accuracy limits.

## Generalized moments — Level A

`evaluate_generalized_moment` uses the same terminal sum with real exponents.
Local `generalized_power` rejects a negative abscissa with a noninteger
exponent; `is_integer_exponent` means exact `trunc(exponent) == exponent`, with
no rounding tolerance. Otherwise it calls `pow` and rejects nonfinite output.
Thus zero raised to a negative exponent fails, and a zero-weight node does not
bypass an invalid power. Every exponent must be finite; standalone generalized
evaluation permits negative exponents when the node arithmetic is defined and
finite.

The CQMOM coordinate loop validates exponent rank, evaluates each factor, then
checks product and accumulation. No support restriction is inferred from
coordinate names. A finite terminal approximation neither establishes
integrability of a continuous singular integrand nor proves an underlying
distribution is supported on the positive half-line. Source-context generalized
access imposes the additional nonnegative-exponent restriction described in
[sources](moment_sources.md).

Evidence: fractional/integer negative-node, zero-node, malformed exponent, and
overflow assertions in [moment reconstruction tests](../../../tests/unit/moment_reconstruction.cpp).
No operating envelope for exponents approaching an integer from either side,
near-zero nodes with negative powers, or signed cancellation was identified.
