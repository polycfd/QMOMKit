#pragma once

namespace qmomkit
{

/**
 * User-facing Vandermonde solver selection for CQMOM conditional moment solves.
 *
 * `automatic` preserves QMOMKit's default behavior: try the
 * Vandermonde-specific Björck-Pereyra algorithm with residual validation and
 * use the dense solver as fallback when validation or the fast solve fails.
 * `pivoted_qr` is an explicit policy and is not selected by `automatic`.
 */
enum class VandermondeSolver
{
    automatic,
    dense,
    pivoted_qr,
    bjorck_pereyra,
    bjorck_pereyra_with_dense_fallback
};

} // namespace qmomkit
