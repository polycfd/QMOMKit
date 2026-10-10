#pragma once

namespace qmomkit
{

/**
 * User-facing Vandermonde solver selection for CQMOM conditional moment solves.
 *
 * `automatic`, the default, uses the Björck-Pereyra solver.
 * `bjorck_pereyra` selects that solver explicitly.
 * `dense`, `pivoted_qr` and `svd` select direct solvers.
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
