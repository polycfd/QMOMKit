#pragma once

#include "internal/linear_algebra.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

// The weighted conditional Vandermonde matrix A(k, j) = w_j x_j^k, which the
// library solves without forming it. Tests build it independently as the
// reference the solvers are compared against.
namespace qmomkit::test
{

[[nodiscard]] inline internal::DenseMatrix make_conditional_vandermonde_matrix(const std::vector<double>& nodes, const std::vector<double>& weights)
{
    const std::size_t n = nodes.size();
    internal::DenseMatrix matrix(n, n);
    for (std::size_t power = 0; power < n; ++power)
    {
        for (std::size_t node = 0; node < n; ++node)
        {
            matrix(power, node) = std::pow(nodes[node], static_cast<int>(power)) * weights[node];
        }
    }
    return matrix;
}

} // namespace qmomkit::test
