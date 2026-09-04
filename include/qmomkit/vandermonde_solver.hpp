#pragma once

namespace qmomkit
{

/**
 * User-facing Vandermonde solver selection for CQMOM conditional moment solves.
 *
 * `automatic` is the default policy and uses the Björck-Pereyra solver.
 * `pivoted_qr` and `svd` select their respective direct solvers.
 */
enum class VandermondeSolver
{
    automatic,
    dense,
    pivoted_qr,
    svd,
    bjorck_pereyra
};

} // namespace qmomkit
