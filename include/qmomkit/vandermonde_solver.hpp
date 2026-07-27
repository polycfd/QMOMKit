#pragma once

namespace qmomkit
{

/**
 * User-facing Vandermonde solver selection for CQMOM conditional moment solves.
 *
 * `automatic` preserves QMOMKit's default behavior: try the
 * Vandermonde-specific Björck-Pereyra algorithm with residual validation and
 * use the dense solver as fallback when validation or the fast solve fails.
 */
enum class VandermondeSolver
{
    automatic,
    dense,
    bjorck_pereyra,
    bjorck_pereyra_with_dense_fallback
};

} // namespace qmomkit
