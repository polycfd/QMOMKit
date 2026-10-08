#include "linear_algebra.hpp"

#include "performance.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
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

} // namespace

double matrix_infinity_norm(const DenseMatrix& matrix)
{
    double norm = 0.0;
    for (std::size_t row = 0; row < matrix.rows(); ++row)
    {
        double row_sum = 0.0;
        for (std::size_t col = 0; col < matrix.cols(); ++col)
        {
            row_sum += std::abs(matrix(row, col));
        }
        norm = std::max(norm, row_sum);
    }
    return norm;
}

LinearSolveResult solve_dense(DenseMatrix matrix, DenseMatrix rhs, const double singular_tolerance, const std::string& singular_message)
{
    performance::record_dense_linear_solve();
    const std::size_t n = matrix.rows();
    const std::size_t rhs_cols = rhs.cols();
    // Pivots are compared with the matrix itself, so scaling the whole matrix,
    // for example by the number density it carries, does not change the result.
    const double scale = matrix_infinity_norm(matrix);

    for (std::size_t pivot_col = 0; pivot_col < n; ++pivot_col)
    {
        std::size_t pivot_row = pivot_col;
        double max_abs = std::abs(matrix(pivot_col, pivot_col));
        for (std::size_t row = pivot_col + 1; row < n; ++row)
        {
            const double candidate = std::abs(matrix(row, pivot_col));
            if (candidate > max_abs)
            {
                max_abs = candidate;
                pivot_row = row;
            }
        }

        if (max_abs <= singular_tolerance * scale)
        {
            return LinearSolveResult{DenseMatrix(n, rhs_cols), false, singular_message};
        }

        if (pivot_row != pivot_col)
        {
            for (std::size_t col = 0; col < n; ++col)
            {
                std::swap(matrix(pivot_col, col), matrix(pivot_row, col));
            }
            for (std::size_t col = 0; col < rhs_cols; ++col)
            {
                std::swap(rhs(pivot_col, col), rhs(pivot_row, col));
            }
        }

        for (std::size_t row = pivot_col + 1; row < n; ++row)
        {
            const double factor = matrix(row, pivot_col) / matrix(pivot_col, pivot_col);
            matrix(row, pivot_col) = 0.0;
            for (std::size_t col = pivot_col + 1; col < n; ++col)
            {
                matrix(row, col) -= factor * matrix(pivot_col, col);
            }
            for (std::size_t col = 0; col < rhs_cols; ++col)
            {
                rhs(row, col) -= factor * rhs(pivot_col, col);
            }
        }
    }

    DenseMatrix solution(n, rhs_cols);
    for (std::size_t reverse_row = n; reverse_row > 0; --reverse_row)
    {
        const std::size_t row = reverse_row - 1;
        for (std::size_t rhs_col = 0; rhs_col < rhs_cols; ++rhs_col)
        {
            double value = rhs(row, rhs_col);
            for (std::size_t col = row + 1; col < n; ++col)
            {
                value -= matrix(row, col) * solution(col, rhs_col);
            }
            solution(row, rhs_col) = value / matrix(row, row);
        }
    }

    if (!matrix_all_finite(solution))
    {
        return LinearSolveResult{DenseMatrix(n, rhs_cols), false, "Dense linear solve produced non-finite solution values."};
    }

    return LinearSolveResult{std::move(solution), true, ""};
}

DenseMatrix make_identity(const std::size_t n)
{
    DenseMatrix identity(n, n);
    for (std::size_t i = 0; i < n; ++i)
    {
        identity(i, i) = 1.0;
    }
    return identity;
}

double residual_norm(const DenseMatrix& matrix, const DenseMatrix& solution, const DenseMatrix& rhs)
{
    performance::record_residual_check();
    double sum = 0.0;
    for (std::size_t row = 0; row < matrix.rows(); ++row)
    {
        for (std::size_t col = 0; col < rhs.cols(); ++col)
        {
            double reconstructed = 0.0;
            for (std::size_t k = 0; k < matrix.cols(); ++k)
            {
                reconstructed += matrix(row, k) * solution(k, col);
            }
            const double difference = reconstructed - rhs(row, col);
            sum += difference * difference;
        }
    }
    return std::sqrt(sum);
}

double condition_estimate(const DenseMatrix& matrix, const double singular_tolerance)
{
    performance::record_condition_estimate();
    LinearSolveResult inverse = solve_dense(matrix, make_identity(matrix.rows()), singular_tolerance, "");
    if (!inverse.success)
    {
        return std::numeric_limits<double>::infinity();
    }
    return matrix_infinity_norm(matrix) * matrix_infinity_norm(inverse.solution);
}

namespace
{

// The exponent shift that places the largest entry in [1, 2); zero for a zero
// matrix or one with a non-finite entry, which is then left unscaled.
[[nodiscard]] int power_of_two_shift(const DenseMatrix& matrix)
{
    double largest = 0.0;
    for (std::size_t row = 0; row < matrix.rows(); ++row)
    {
        for (std::size_t col = 0; col < matrix.cols(); ++col)
        {
            const double value = std::abs(matrix(row, col));
            if (!std::isfinite(value))
            {
                return 0;
            }
            largest = std::max(largest, value);
        }
    }
    if (largest == 0.0)
    {
        return 0;
    }
    int exponent = 0;
    static_cast<void>(std::frexp(largest, &exponent));
    return 1 - exponent;
}

void rotate(DenseMatrix& matrix, DenseMatrix& eigenvectors, const std::size_t p, const std::size_t q)
{
    if (matrix(p, q) == 0.0)
    {
        return;
    }

    const double tau = (matrix(q, q) - matrix(p, p)) / (2.0 * matrix(p, q));
    const double sign = tau >= 0.0 ? 1.0 : -1.0;
    const double t = sign / (std::abs(tau) + std::sqrt(1.0 + tau * tau));
    const double c = 1.0 / std::sqrt(1.0 + t * t);
    const double s = t * c;
    const double app = matrix(p, p);
    const double aqq = matrix(q, q);
    const double apq = matrix(p, q);

    matrix(p, p) = app - t * apq;
    matrix(q, q) = aqq + t * apq;
    matrix(p, q) = 0.0;
    matrix(q, p) = 0.0;

    for (std::size_t k = 0; k < matrix.rows(); ++k)
    {
        if (k != p && k != q)
        {
            const double akp = matrix(k, p);
            const double akq = matrix(k, q);
            matrix(k, p) = c * akp - s * akq;
            matrix(p, k) = matrix(k, p);
            matrix(k, q) = s * akp + c * akq;
            matrix(q, k) = matrix(k, q);
        }
    }

    for (std::size_t k = 0; k < eigenvectors.rows(); ++k)
    {
        const double vkp = eigenvectors(k, p);
        const double vkq = eigenvectors(k, q);
        eigenvectors(k, p) = c * vkp - s * vkq;
        eigenvectors(k, q) = s * vkp + c * vkq;
    }
}

} // namespace

EigenDecomposition symmetric_jacobi_eigen(DenseMatrix matrix, const JacobiOptions& options)
{
    const std::size_t n = matrix.rows();
    DenseMatrix eigenvectors(n, n);
    for (std::size_t i = 0; i < n; ++i)
    {
        eigenvectors(i, i) = 1.0;
    }

    // Scaling by a power of two is exact, except for entries so far below the
    // largest that they leave the normal range, which lie far below the
    // tolerance. Every rotation parameter is a ratio of entries, so the
    // rotations are those of the unscaled matrix while the tolerance becomes
    // relative to its largest entry.
    const int shift = power_of_two_shift(matrix);
    for (std::size_t row = 0; row < n; ++row)
    {
        for (std::size_t col = 0; col < n; ++col)
        {
            matrix(row, col) = std::ldexp(matrix(row, col), shift);
        }
    }

    const std::size_t max_rotations = options.max_rotations > 0 ? options.max_rotations : 10 * n * n;
    bool converged = false;
    for (std::size_t rotation = 0;; ++rotation)
    {
        // Convergence is tested by the largest off-diagonal magnitude found
        // below. A separate off-diagonal norm check would be redundant: the
        // norm bounds that maximum from above, so it can only stop the
        // iteration where the maximum already does.
        std::size_t p = 0;
        std::size_t q = 0;
        double max_value = 0.0;
        for (std::size_t row = 0; row < n; ++row)
        {
            for (std::size_t col = row + 1; col < n; ++col)
            {
                const double value = std::abs(matrix(row, col));
                if (value > max_value)
                {
                    max_value = value;
                    p = row;
                    q = col;
                }
            }
        }

        if (max_value <= options.tolerance)
        {
            converged = true;
            break;
        }
        if (rotation == max_rotations)
        {
            break;
        }
        rotate(matrix, eigenvectors, p, q);
    }

    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&matrix](const std::size_t lhs, const std::size_t rhs) { return matrix(lhs, lhs) < matrix(rhs, rhs); });

    EigenDecomposition result{std::vector<double>(n, 0.0), DenseMatrix(n, n)};
    for (std::size_t sorted_col = 0; sorted_col < n; ++sorted_col)
    {
        const std::size_t original_col = order[sorted_col];
        result.eigenvalues[sorted_col] = std::ldexp(matrix(original_col, original_col), -shift);
        for (std::size_t row = 0; row < n; ++row)
        {
            result.eigenvectors(row, sorted_col) = eigenvectors(row, original_col);
        }
    }

    result.converged = converged;
    return result;
}

} // namespace qmomkit::internal
