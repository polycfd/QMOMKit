#include "bjorck_pereyra.hpp"

#include "numeric_utils.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>

namespace qmomkit::internal
{
namespace
{

[[nodiscard]] bool matrix_all_finite(const DenseMatrix& matrix)
{
    for (std::size_t row = 0; row < matrix.rows(); ++row)
    {
        for (std::size_t col = 0; col < matrix.cols(); ++col)
        {
            if (!std::isfinite(matrix(row, col)))
            {
                return false;
            }
        }
    }
    return true;
}

[[nodiscard]] LinearSolveResult failure_result(const std::size_t rows, const std::size_t cols, const std::string& message)
{
    return LinearSolveResult{DenseMatrix(rows, cols), false, message};
}

[[nodiscard]] double node_scale(const double lhs, const double rhs)
{
    return std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

[[nodiscard]] bool denominator_is_singular(const double denominator, const double lhs, const double rhs, const double singular_tolerance)
{
    return std::abs(denominator) <= singular_tolerance * node_scale(lhs, rhs);
}

} // namespace

LinearSolveResult solve_bjorck_pereyra(const std::vector<double>& nodes, DenseMatrix rhs, const double singular_tolerance)
{
    const std::size_t n = nodes.size();
    const std::size_t rhs_cols = rhs.cols();
    if (n == 0)
    {
        return failure_result(n, rhs_cols, "Transposed Vandermonde system must not be empty.");
    }
    if (rhs.rows() != n)
    {
        return failure_result(n, rhs_cols, "Transposed Vandermonde dimensions are inconsistent.");
    }
    if (!std::isfinite(singular_tolerance) || singular_tolerance <= 0.0)
    {
        return failure_result(n, rhs_cols, "Transposed Vandermonde singular tolerance is invalid.");
    }
    if (!qmomkit::internal::all_finite(nodes) || !matrix_all_finite(rhs))
    {
        return failure_result(n, rhs_cols, "Transposed Vandermonde inputs must be finite.");
    }

    // Apply the transposed elementary operations of the Bjorck-Pereyra
    // Vandermonde solve. The result is Y in V^T Y = rhs.
    for (std::size_t k = 0; k + 1 < n; ++k)
    {
        for (std::size_t reverse_i = n - 1; reverse_i > k; --reverse_i)
        {
            const std::size_t i = reverse_i - 1;
            for (std::size_t col = 0; col < rhs_cols; ++col)
            {
                rhs(i + 1, col) -= nodes[k] * rhs(i, col);
            }
        }
    }

    for (std::size_t reverse_k = n - 1; reverse_k > 0; --reverse_k)
    {
        const std::size_t k = reverse_k - 1;
        for (std::size_t i = k + 1; i < n; ++i)
        {
            const std::size_t denominator_node = i - k - 1;
            const double denominator = nodes[i] - nodes[denominator_node];
            if (denominator_is_singular(denominator, nodes[i], nodes[denominator_node], singular_tolerance))
            {
                return failure_result(n, rhs_cols, "Transposed Vandermonde system is singular.");
            }
            for (std::size_t col = 0; col < rhs_cols; ++col)
            {
                const double value = rhs(i, col) / denominator;
                rhs(i - 1, col) -= value;
                rhs(i, col) = value;
            }
        }
    }

    if (!matrix_all_finite(rhs))
    {
        return failure_result(n, rhs_cols, "Transposed Vandermonde solve produced non-finite values.");
    }

    return LinearSolveResult{std::move(rhs), true, ""};
}

} // namespace qmomkit::internal
