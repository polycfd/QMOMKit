#include "jacobi_svd.hpp"

#include "linear_algebra.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <utility>

namespace qmomkit::internal
{
namespace
{

[[nodiscard]] bool finite_matrix(const DenseMatrix& matrix)
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

[[nodiscard]] bool finite_values(const std::vector<double>& values)
{
    return std::all_of(values.begin(), values.end(), [](const double value) { return std::isfinite(value); });
}

[[nodiscard]] double max_abs_entry(const DenseMatrix& matrix)
{
    double max_abs = 0.0;
    for (std::size_t row = 0; row < matrix.rows(); ++row)
    {
        for (std::size_t col = 0; col < matrix.cols(); ++col)
        {
            max_abs = std::max(max_abs, std::abs(matrix(row, col)));
        }
    }
    return max_abs;
}

void scale_matrix(DenseMatrix& matrix, const double scale)
{
    for (std::size_t row = 0; row < matrix.rows(); ++row)
    {
        for (std::size_t col = 0; col < matrix.cols(); ++col)
        {
            matrix(row, col) *= scale;
        }
    }
}

[[nodiscard]] double stable_column_norm(const DenseMatrix& matrix, const std::size_t col)
{
    double norm = 0.0;
    for (std::size_t row = 0; row < matrix.rows(); ++row)
    {
        norm = std::hypot(norm, matrix(row, col));
    }
    return norm;
}

[[nodiscard]] double stable_column_dot(const DenseMatrix& matrix, const std::size_t lhs, const std::size_t rhs)
{
    long double dot = 0.0L;
    for (std::size_t row = 0; row < matrix.rows(); ++row)
    {
        dot += static_cast<long double>(matrix(row, lhs)) * static_cast<long double>(matrix(row, rhs));
    }
    return static_cast<double>(dot);
}

void rotate_columns(DenseMatrix& matrix, const std::size_t p, const std::size_t q, const double c, const double s)
{
    for (std::size_t row = 0; row < matrix.rows(); ++row)
    {
        const double bp = matrix(row, p);
        const double bq = matrix(row, q);
        matrix(row, p) = c * bp - s * bq;
        matrix(row, q) = s * bp + c * bq;
    }
}

void copy_column(const DenseMatrix& source, const std::size_t source_col, DenseMatrix& destination, const std::size_t destination_col)
{
    for (std::size_t row = 0; row < source.rows(); ++row)
    {
        destination(row, destination_col) = source(row, source_col);
    }
}

[[nodiscard]] double resolved_jacobi_tolerance(const std::size_t rows, const std::size_t cols, const double requested)
{
    if (requested != 0.0)
    {
        return requested;
    }
    const double dimension = static_cast<double>(std::max(rows, cols));
    return std::max(1.0e-14, 10.0 * std::numeric_limits<double>::epsilon() * dimension);
}

[[nodiscard]] std::size_t resolved_max_sweeps(const std::size_t cols, const std::size_t requested)
{
    if (requested != 0)
    {
        return requested;
    }
    return std::max<std::size_t>(32, 8 * cols);
}

[[nodiscard]] int resolved_scale_exponent(const double max_abs)
{
    int exponent = 0;
    (void)std::frexp(max_abs, &exponent);
    const int requested = -exponent;
    const int max_finite_exponent = std::numeric_limits<double>::max_exponent - 1;
    const int min_finite_exponent = std::numeric_limits<double>::min_exponent - std::numeric_limits<double>::digits;
    return std::max(min_finite_exponent, std::min(max_finite_exponent, requested));
}

[[nodiscard]] bool valid_options(const SvdOptions& options)
{
    return std::isfinite(options.relative_rank_tolerance) && options.relative_rank_tolerance > 0.0 && std::isfinite(options.jacobi_tolerance) &&
           options.jacobi_tolerance >= 0.0;
}

[[nodiscard]] SvdFactorization factorization_failure_result(DenseMatrix matrix, const SvdOptions& options, std::string message)
{
    SvdFactorization result;
    result.rows = matrix.rows();
    result.cols = matrix.cols();
    result.left_singular_vectors = DenseMatrix(matrix.rows(), matrix.cols());
    result.singular_values.assign(matrix.cols(), 0.0);
    result.right_singular_vectors = make_identity(matrix.cols());
    result.relative_rank_tolerance = options.relative_rank_tolerance;
    result.jacobi_tolerance = options.jacobi_tolerance;
    result.status = SvdStatus::invalid_input;
    result.message = std::move(message);
    return result;
}

[[nodiscard]] SvdSolveResult solve_failure_result(const std::size_t rows, const std::size_t cols, const std::size_t numerical_rank,
                                                  const SvdStatus status, std::string message)
{
    return SvdSolveResult{DenseMatrix(rows, cols), numerical_rank, status, std::move(message)};
}

[[nodiscard]] std::size_t numerical_rank_from_singular_values(const std::vector<double>& singular_values, const double relative_rank_tolerance)
{
    if (singular_values.empty() || singular_values.front() == 0.0)
    {
        return 0;
    }

    const double threshold = relative_rank_tolerance * singular_values.front();
    std::size_t rank = 0;
    while (rank < singular_values.size() && singular_values[rank] > threshold)
    {
        ++rank;
    }
    return rank;
}

[[nodiscard]] bool singular_values_are_nonnegative_descending(const std::vector<double>& singular_values)
{
    for (std::size_t i = 0; i < singular_values.size(); ++i)
    {
        if (singular_values[i] < 0.0)
        {
            return false;
        }
        if (i > 0 && singular_values[i - 1] < singular_values[i])
        {
            return false;
        }
    }
    return true;
}

[[nodiscard]] bool valid_svd_status(const SvdStatus status)
{
    switch (status)
    {
    case SvdStatus::success:
    case SvdStatus::invalid_input:
    case SvdStatus::rank_deficient:
    case SvdStatus::no_convergence:
        return true;
    }
    return false;
}

[[nodiscard]] bool rank_metadata_is_consistent(const SvdFactorization& factorization)
{
    const std::size_t computed_rank = numerical_rank_from_singular_values(factorization.singular_values, factorization.relative_rank_tolerance);
    if (factorization.numerical_rank != computed_rank)
    {
        return false;
    }

    switch (factorization.status)
    {
    case SvdStatus::success:
        return computed_rank == factorization.cols;
    case SvdStatus::rank_deficient:
        return computed_rank < factorization.cols;
    case SvdStatus::no_convergence:
    case SvdStatus::invalid_input:
        return true;
    }
    return false;
}

[[nodiscard]] std::size_t diagnostic_rank_from_metadata(const SvdFactorization& factorization)
{
    if (factorization.singular_values.size() == factorization.cols && finite_values(factorization.singular_values) &&
        singular_values_are_nonnegative_descending(factorization.singular_values) && std::isfinite(factorization.relative_rank_tolerance) &&
        factorization.relative_rank_tolerance > 0.0)
    {
        return numerical_rank_from_singular_values(factorization.singular_values, factorization.relative_rank_tolerance);
    }
    return 0;
}

[[nodiscard]] bool full_rank_singular_values_are_safe(const SvdFactorization& factorization)
{
    if (factorization.cols == 0 || factorization.singular_values.empty())
    {
        return false;
    }
    const double threshold = factorization.relative_rank_tolerance * factorization.singular_values.front();
    return factorization.singular_values[factorization.cols - 1] > threshold;
}

[[nodiscard]] bool factorization_metadata_is_consistent(const SvdFactorization& factorization)
{
    return factorization.rows > 0 && factorization.cols > 0 && factorization.rows >= factorization.cols &&
           factorization.left_singular_vectors.rows() == factorization.rows && factorization.left_singular_vectors.cols() == factorization.cols &&
           factorization.right_singular_vectors.rows() == factorization.cols && factorization.right_singular_vectors.cols() == factorization.cols &&
           factorization.singular_values.size() == factorization.cols && finite_matrix(factorization.left_singular_vectors) &&
           finite_matrix(factorization.right_singular_vectors) && finite_values(factorization.singular_values) &&
           singular_values_are_nonnegative_descending(factorization.singular_values) && valid_svd_status(factorization.status) &&
           std::isfinite(factorization.relative_rank_tolerance) && factorization.relative_rank_tolerance > 0.0 &&
           std::isfinite(factorization.jacobi_tolerance) && factorization.jacobi_tolerance > 0.0 &&
           factorization.numerical_rank <= factorization.cols && rank_metadata_is_consistent(factorization) &&
           (factorization.status != SvdStatus::success || full_rank_singular_values_are_safe(factorization));
}

[[nodiscard]] bool pair_is_converged(const double alpha, const double beta, const double gamma, const double tolerance)
{
    if (alpha == 0.0 || beta == 0.0 || gamma == 0.0)
    {
        return true;
    }
    return std::abs(gamma) <= tolerance * std::sqrt(alpha * beta);
}

[[nodiscard]] double jacobi_tangent(const double alpha, const double beta, const double gamma)
{
    const double zeta = (beta - alpha) / (2.0 * gamma);
    if (std::isfinite(zeta))
    {
        const double sign = zeta >= 0.0 ? 1.0 : -1.0;
        return sign / (std::abs(zeta) + std::hypot(1.0, zeta));
    }

    const double denominator = beta - alpha;
    if (denominator != 0.0 && std::isfinite(denominator))
    {
        return gamma / denominator;
    }
    return std::copysign(1.0, gamma);
}

[[nodiscard]] bool rotate_pair(DenseMatrix& working, DenseMatrix& right_singular_vectors, const std::size_t p, const std::size_t q,
                               const double tolerance)
{
    const double alpha = stable_column_dot(working, p, p);
    const double beta = stable_column_dot(working, q, q);
    const double gamma = stable_column_dot(working, p, q);
    if (pair_is_converged(alpha, beta, gamma, tolerance))
    {
        return false;
    }

    const double t = jacobi_tangent(alpha, beta, gamma);
    const double c = 1.0 / std::hypot(1.0, t);
    const double s = t * c;
    rotate_columns(working, p, q, c, s);
    rotate_columns(right_singular_vectors, p, q, c, s);
    return true;
}

void sort_singular_triplets(SvdFactorization& factorization)
{
    std::vector<std::size_t> order(factorization.cols);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&factorization](const std::size_t lhs, const std::size_t rhs)
                     { return factorization.singular_values[lhs] > factorization.singular_values[rhs]; });

    std::vector<double> sorted_singular_values(factorization.cols, 0.0);
    DenseMatrix sorted_left(factorization.rows, factorization.cols);
    DenseMatrix sorted_right(factorization.cols, factorization.cols);
    for (std::size_t sorted_col = 0; sorted_col < factorization.cols; ++sorted_col)
    {
        const std::size_t original_col = order[sorted_col];
        sorted_singular_values[sorted_col] = factorization.singular_values[original_col];
        copy_column(factorization.left_singular_vectors, original_col, sorted_left, sorted_col);
        copy_column(factorization.right_singular_vectors, original_col, sorted_right, sorted_col);
    }

    factorization.singular_values = std::move(sorted_singular_values);
    factorization.left_singular_vectors = std::move(sorted_left);
    factorization.right_singular_vectors = std::move(sorted_right);
}

[[nodiscard]] SvdFactorization completed_factorization(DenseMatrix working, DenseMatrix right_singular_vectors, const SvdOptions& options,
                                                       const double jacobi_tolerance, const std::size_t sweeps, const int scale_exponent,
                                                       const SvdStatus convergence_status)
{
    const std::size_t rows = working.rows();
    const std::size_t cols = working.cols();

    SvdFactorization result;
    result.rows = rows;
    result.cols = cols;
    result.left_singular_vectors = DenseMatrix(rows, cols);
    result.singular_values.assign(cols, 0.0);
    result.right_singular_vectors = std::move(right_singular_vectors);
    result.sweeps = sweeps;
    result.relative_rank_tolerance = options.relative_rank_tolerance;
    result.jacobi_tolerance = jacobi_tolerance;

    std::vector<double> scaled_singular_values(cols, 0.0);
    for (std::size_t col = 0; col < cols; ++col)
    {
        scaled_singular_values[col] = stable_column_norm(working, col);
        result.singular_values[col] = std::ldexp(scaled_singular_values[col], -scale_exponent);
        if (scaled_singular_values[col] > 0.0)
        {
            for (std::size_t row = 0; row < rows; ++row)
            {
                result.left_singular_vectors(row, col) = working(row, col) / scaled_singular_values[col];
            }
        }
    }

    sort_singular_triplets(result);
    result.numerical_rank = numerical_rank_from_singular_values(result.singular_values, options.relative_rank_tolerance);

    if (convergence_status == SvdStatus::no_convergence)
    {
        result.status = SvdStatus::no_convergence;
        result.message = "SVD Jacobi iteration did not converge.";
    }
    else if (result.numerical_rank < cols)
    {
        result.status = SvdStatus::rank_deficient;
        result.message = "SVD factorization is rank deficient.";
    }

    if (!finite_matrix(result.left_singular_vectors) || !finite_matrix(result.right_singular_vectors) || !finite_values(result.singular_values))
    {
        result.status = SvdStatus::invalid_input;
        result.message = "SVD factorization produced non-finite values.";
    }

    return result;
}

} // namespace

SvdFactorization factorize_svd(DenseMatrix matrix, const SvdOptions& options)
{
    const std::size_t rows = matrix.rows();
    const std::size_t cols = matrix.cols();
    if (rows == 0 || cols == 0)
    {
        return factorization_failure_result(std::move(matrix), options, "SVD factorization matrix must not be empty.");
    }
    if (rows < cols)
    {
        return factorization_failure_result(std::move(matrix), options, "SVD factorization requires rows >= columns.");
    }
    if (!valid_options(options))
    {
        return factorization_failure_result(std::move(matrix), options, "SVD options are invalid.");
    }
    if (!finite_matrix(matrix))
    {
        return factorization_failure_result(std::move(matrix), options, "SVD factorization matrix entries must be finite.");
    }

    const double jacobi_tolerance = resolved_jacobi_tolerance(rows, cols, options.jacobi_tolerance);
    const std::size_t max_sweeps = resolved_max_sweeps(cols, options.max_sweeps);
    const double max_abs = max_abs_entry(matrix);
    DenseMatrix right_singular_vectors = make_identity(cols);

    if (max_abs == 0.0)
    {
        return completed_factorization(std::move(matrix), std::move(right_singular_vectors), options, jacobi_tolerance, 0, 0, SvdStatus::success);
    }

    const int scale_exponent = resolved_scale_exponent(max_abs);
    const double scale = std::ldexp(1.0, scale_exponent);
    scale_matrix(matrix, scale);

    std::size_t completed_sweeps = 0;
    bool converged = false;
    for (std::size_t sweep = 0; sweep < max_sweeps; ++sweep)
    {
        bool rotated = false;
        for (std::size_t p = 0; p < cols; ++p)
        {
            for (std::size_t q = p + 1; q < cols; ++q)
            {
                rotated = rotate_pair(matrix, right_singular_vectors, p, q, jacobi_tolerance) || rotated;
            }
        }
        ++completed_sweeps;

        if (!finite_matrix(matrix) || !finite_matrix(right_singular_vectors))
        {
            return factorization_failure_result(std::move(matrix), options, "SVD Jacobi iteration produced non-finite values.");
        }
        if (!rotated)
        {
            converged = true;
            break;
        }
    }

    return completed_factorization(std::move(matrix), std::move(right_singular_vectors), options, jacobi_tolerance, completed_sweeps, scale_exponent,
                                   converged ? SvdStatus::success : SvdStatus::no_convergence);
}

SvdSolveResult solve_svd(const SvdFactorization& factorization, DenseMatrix rhs)
{
    const std::size_t rhs_cols = rhs.cols();
    if (!factorization_metadata_is_consistent(factorization))
    {
        return solve_failure_result(factorization.cols, rhs_cols, diagnostic_rank_from_metadata(factorization), SvdStatus::invalid_input,
                                    "SVD factorization metadata is inconsistent.");
    }
    if (rhs.rows() != factorization.rows)
    {
        return solve_failure_result(factorization.cols, rhs_cols, factorization.numerical_rank, SvdStatus::invalid_input,
                                    "SVD solve right-hand side row count is inconsistent.");
    }
    if (!finite_matrix(rhs))
    {
        return solve_failure_result(factorization.cols, rhs_cols, factorization.numerical_rank, SvdStatus::invalid_input,
                                    "SVD solve right-hand side entries must be finite.");
    }
    if (factorization.status == SvdStatus::rank_deficient || factorization.numerical_rank < factorization.cols)
    {
        return solve_failure_result(factorization.cols, rhs_cols, factorization.numerical_rank, SvdStatus::rank_deficient,
                                    "SVD solve requires full column rank.");
    }
    if (factorization.status == SvdStatus::no_convergence)
    {
        return solve_failure_result(factorization.cols, rhs_cols, factorization.numerical_rank, SvdStatus::no_convergence,
                                    "SVD solve requires a converged factorization.");
    }
    if (factorization.status != SvdStatus::success)
    {
        return solve_failure_result(factorization.cols, rhs_cols, factorization.numerical_rank, factorization.status, factorization.message);
    }

    DenseMatrix projected(factorization.cols, rhs_cols);
    for (std::size_t row = 0; row < factorization.cols; ++row)
    {
        for (std::size_t rhs_col = 0; rhs_col < rhs_cols; ++rhs_col)
        {
            double value = 0.0;
            for (std::size_t k = 0; k < factorization.rows; ++k)
            {
                value += factorization.left_singular_vectors(k, row) * rhs(k, rhs_col);
            }
            projected(row, rhs_col) = value / factorization.singular_values[row];
        }
    }

    DenseMatrix solution(factorization.cols, rhs_cols);
    for (std::size_t row = 0; row < factorization.cols; ++row)
    {
        for (std::size_t rhs_col = 0; rhs_col < rhs_cols; ++rhs_col)
        {
            double value = 0.0;
            for (std::size_t k = 0; k < factorization.cols; ++k)
            {
                value += factorization.right_singular_vectors(row, k) * projected(k, rhs_col);
            }
            solution(row, rhs_col) = value;
        }
    }

    if (!finite_matrix(solution))
    {
        return solve_failure_result(factorization.cols, rhs_cols, factorization.numerical_rank, SvdStatus::invalid_input,
                                    "SVD solve produced non-finite solution values.");
    }

    return SvdSolveResult{std::move(solution), factorization.numerical_rank, SvdStatus::success, ""};
}

SvdSolveResult solve_svd(DenseMatrix matrix, DenseMatrix rhs, const SvdOptions& options)
{
    SvdFactorization factorization = factorize_svd(std::move(matrix), options);
    if (factorization.status == SvdStatus::invalid_input)
    {
        return solve_failure_result(factorization.cols, rhs.cols(), factorization.numerical_rank, factorization.status, factorization.message);
    }
    return solve_svd(factorization, std::move(rhs));
}

SvdSolveResult solve_svd_pseudoinverse(const SvdFactorization& factorization, DenseMatrix rhs, const double rcond)
{
    const std::size_t rhs_cols = rhs.cols();
    if (!factorization_metadata_is_consistent(factorization))
    {
        return solve_failure_result(factorization.cols, rhs_cols, diagnostic_rank_from_metadata(factorization), SvdStatus::invalid_input,
                                    "SVD factorization metadata is inconsistent.");
    }
    if (rhs.rows() != factorization.rows)
    {
        return solve_failure_result(factorization.cols, rhs_cols, factorization.numerical_rank, SvdStatus::invalid_input,
                                    "SVD solve right-hand side row count is inconsistent.");
    }
    if (!finite_matrix(rhs))
    {
        return solve_failure_result(factorization.cols, rhs_cols, factorization.numerical_rank, SvdStatus::invalid_input,
                                    "SVD solve right-hand side entries must be finite.");
    }
    if (!std::isfinite(rcond) || rcond < 0.0 || rcond >= 1.0)
    {
        return solve_failure_result(factorization.cols, rhs_cols, factorization.numerical_rank, SvdStatus::invalid_input,
                                    "SVD pseudoinverse rcond must be finite and satisfy 0 <= rcond < 1.");
    }
    if (factorization.status == SvdStatus::no_convergence)
    {
        return solve_failure_result(factorization.cols, rhs_cols, factorization.numerical_rank, SvdStatus::no_convergence,
                                    "SVD pseudoinverse requires a converged factorization.");
    }
    if (factorization.status == SvdStatus::invalid_input)
    {
        return solve_failure_result(factorization.cols, rhs_cols, factorization.numerical_rank, factorization.status, factorization.message);
    }

    const double cutoff = rcond * factorization.singular_values.front();
    DenseMatrix projected(factorization.cols, rhs_cols);
    std::size_t retained_rank = 0;
    for (std::size_t row = 0; row < factorization.cols; ++row)
    {
        const double singular_value = factorization.singular_values[row];
        if (singular_value <= cutoff)
        {
            continue;
        }
        ++retained_rank;
        for (std::size_t rhs_col = 0; rhs_col < rhs_cols; ++rhs_col)
        {
            double value = 0.0;
            for (std::size_t k = 0; k < factorization.rows; ++k)
            {
                value += factorization.left_singular_vectors(k, row) * rhs(k, rhs_col);
            }
            projected(row, rhs_col) = value / singular_value;
        }
    }

    DenseMatrix solution(factorization.cols, rhs_cols);
    for (std::size_t row = 0; row < factorization.cols; ++row)
    {
        for (std::size_t rhs_col = 0; rhs_col < rhs_cols; ++rhs_col)
        {
            double value = 0.0;
            for (std::size_t k = 0; k < factorization.cols; ++k)
            {
                value += factorization.right_singular_vectors(row, k) * projected(k, rhs_col);
            }
            solution(row, rhs_col) = value;
        }
    }

    if (!finite_matrix(solution))
    {
        return solve_failure_result(factorization.cols, rhs_cols, retained_rank, SvdStatus::invalid_input,
                                    "SVD pseudoinverse solve produced non-finite solution values.");
    }
    return SvdSolveResult{std::move(solution), retained_rank, SvdStatus::success, ""};
}

} // namespace qmomkit::internal
