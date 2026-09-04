#pragma once

namespace qmomkit
{

/**
 * User-facing Vandermonde solver selection for CQMOM conditional moment solves.
 *
 * `automatic` is the default policy and uses the Björck-Pereyra solver.
 * `bjorck_pereyra_with_dense_fallback` retries with the dense solver when the
 * Björck-Pereyra solve fails or is rejected.
 * `pivoted_qr` and `svd` select their respective direct solvers.
 */
enum class VandermondeSolver
{
    automatic,
    dense,
    pivoted_qr,
    svd,
    bjorck_pereyra,
    bjorck_pereyra_with_dense_fallback
};

} // namespace qmomkit
