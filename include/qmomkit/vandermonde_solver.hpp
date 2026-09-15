#pragma once

namespace qmomkit
{

/**
 * User-facing Vandermonde solver selection for CQMOM conditional moment solves.
 *
 * `automatic` is the default policy and uses the Björck-Pereyra solver.
 * `dense`,`pivoted_qr` and `svd` select direct solvers.
 * `moore_penrose` selects a truncated SVD pseudoinverse solver.
 */
enum class VandermondeSolver
{
    automatic,
    dense,
    pivoted_qr,
    svd,
    moore_penrose,
    bjorck_pereyra
};

} // namespace qmomkit
