#pragma once

#include "linear_algebra.hpp"

#include <vector>

namespace qmomkit::internal
{

/**
 * Solve the unweighted transposed Vandermonde system V^T Y = B using the
 * Bjorck-Pereyra algorithm, where V_{i,p} = nodes[i]^p and B is the rhs
 * argument.
 *
 * Rows of rhs are powers 0..N-1, rows of the returned solution are nodes, and
 * each rhs column is solved independently. This kernel does not know about
 * quadrature weights, conditional moments, or solver policy.
 */
[[nodiscard]] LinearSolveResult solve_bjorck_pereyra(const std::vector<double>& nodes, DenseMatrix rhs, double singular_tolerance);

} // namespace qmomkit::internal
