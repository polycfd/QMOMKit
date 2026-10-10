#include "bjorck_pereyra.hpp"

#include "numeric_utils.hpp"
#include "performance.hpp"

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

// The system counts as singular when its condition number reaches
// 1 / singular_tolerance; the dense solver's pivot test uses the same
// tolerance. With the nodes divided by the power of two above the largest |x|,
// every entry of V^T lies in [-1, 1], so ||V^T||_inf <= n, and
// Gautschi's bound ||V^-T||_inf <= max_i prod_{j != i} (1 + |x_j|) / |x_i - x_j|
// (Numer. Math. 4, 117-123, 1962, Theorem 1; an equality when the nodes share
// a sign) gives the condition number up to the factor n. The estimate does not
// depend on the unit of the coordinate, covers several nodes that are close
// together as well as a single close pair, and is infinite for coincident
// nodes.
[[nodiscard]] bool is_singular(const std::vector<double>& nodes, const double singular_tolerance)
{
    double largest_node = 0.0;
    for (const double node : nodes)
    {
        largest_node = std::max(largest_node, std::abs(node));
    }
    int exponent = 0;
    if (largest_node > 0.0)
    {
        static_cast<void>(std::frexp(largest_node, &exponent));
    }
    // (1 + |x_j| / s) / (|x_i - x_j| / s) = (s + |x_j|) / |x_i - x_j|.
    const double scale = std::ldexp(1.0, exponent);

    const std::size_t n = nodes.size();
    double inverse_bound = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        double product = 1.0;
        for (std::size_t j = 0; j < n; ++j)
        {
            if (j != i)
            {
                product *= (scale + std::abs(nodes[j])) / std::abs(nodes[i] - nodes[j]);
            }
        }
        inverse_bound = std::max(inverse_bound, product);
    }
    return !(static_cast<double>(n) * inverse_bound * singular_tolerance < 1.0);
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
    if (is_singular(nodes, singular_tolerance))
    {
        return failure_result(n, rhs_cols, "Transposed Vandermonde system is singular.");
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
                performance::record_bjorck_pereyra_transform_update();
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
            for (std::size_t col = 0; col < rhs_cols; ++col)
            {
                const double value = rhs(i, col) / denominator;
                rhs(i - 1, col) -= value;
                rhs(i, col) = value;
                performance::record_bjorck_pereyra_transform_update();
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
