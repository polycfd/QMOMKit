#include "householder_qr.hpp"

#include <algorithm>
#include <cmath>
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

[[nodiscard]] double stable_column_norm(const DenseMatrix& matrix, const std::size_t first_row, const std::size_t col)
{
    double norm = 0.0;
    for (std::size_t row = first_row; row < matrix.rows(); ++row)
    {
        norm = std::hypot(norm, matrix(row, col));
    }
    return norm;
}

void swap_columns(DenseMatrix& matrix, const std::size_t lhs, const std::size_t rhs)
{
    if (lhs == rhs)
    {
        return;
    }
    for (std::size_t row = 0; row < matrix.rows(); ++row)
    {
        std::swap(matrix(row, lhs), matrix(row, rhs));
    }
}

[[nodiscard]] PivotedQrFactorization qr_failure_result(DenseMatrix matrix, const double relative_rank_tolerance, std::string message)
{
    PivotedQrFactorization result;
    result.rows = matrix.rows();
    result.cols = matrix.cols();
    result.factors = std::move(matrix);
    result.relative_rank_tolerance = relative_rank_tolerance;
    result.status = QrStatus::invalid_input;
    result.message = std::move(message);
    return result;
}

[[nodiscard]] QrSolveResult solve_failure_result(const std::size_t rows, const std::size_t cols, const std::size_t numerical_rank,
                                                 const QrStatus status, std::string message)
{
    return QrSolveResult{DenseMatrix(rows, cols), numerical_rank, status, std::move(message)};
}

[[nodiscard]] std::size_t numerical_rank_from_factor(const DenseMatrix& factors, const std::size_t cols, const double relative_rank_tolerance)
{
    if (cols == 0 || factors.rows() == 0 || std::abs(factors(0, 0)) == 0.0)
    {
        return 0;
    }

    const double threshold = relative_rank_tolerance * std::abs(factors(0, 0));
    std::size_t rank = 0;
    while (rank < cols && std::abs(factors(rank, rank)) > threshold)
    {
        ++rank;
    }
    return rank;
}

[[nodiscard]] bool valid_permutation(const std::vector<std::size_t>& permutation, const std::size_t cols)
{
    if (permutation.size() != cols)
    {
        return false;
    }
    std::vector<bool> seen(cols, false);
    for (const std::size_t value : permutation)
    {
        if (value >= cols || seen[value])
        {
            return false;
        }
        seen[value] = true;
    }
    return true;
}

[[nodiscard]] bool valid_qr_status(const QrStatus status)
{
    switch (status)
    {
    case QrStatus::success:
    case QrStatus::invalid_input:
    case QrStatus::rank_deficient:
        return true;
    }
    return false;
}

[[nodiscard]] bool finite_values(const std::vector<double>& values)
{
    return std::all_of(values.begin(), values.end(), [](const double value) { return std::isfinite(value); });
}

[[nodiscard]] bool rank_metadata_is_consistent(const PivotedQrFactorization& factorization)
{
    const std::size_t computed_rank = numerical_rank_from_factor(factorization.factors, factorization.cols, factorization.relative_rank_tolerance);
    if (factorization.numerical_rank != computed_rank)
    {
        return false;
    }

    switch (factorization.status)
    {
    case QrStatus::success:
        return computed_rank == factorization.cols;
    case QrStatus::rank_deficient:
        return computed_rank < factorization.cols;
    case QrStatus::invalid_input:
        return true;
    }
    return false;
}

[[nodiscard]] std::size_t diagnostic_rank_from_metadata(const PivotedQrFactorization& factorization)
{
    if (factorization.rows == factorization.factors.rows() && factorization.cols == factorization.factors.cols() &&
        factorization.rows >= factorization.cols && factorization.cols > 0 && std::isfinite(factorization.relative_rank_tolerance) &&
        factorization.relative_rank_tolerance > 0.0 && finite_matrix(factorization.factors))
    {
        return numerical_rank_from_factor(factorization.factors, factorization.cols, factorization.relative_rank_tolerance);
    }
    return 0;
}

[[nodiscard]] bool factorization_metadata_is_consistent(const PivotedQrFactorization& factorization)
{
    return factorization.rows > 0 && factorization.cols > 0 && factorization.rows == factorization.factors.rows() &&
           factorization.cols == factorization.factors.cols() && factorization.rows >= factorization.cols &&
           factorization.householder_coefficients.size() == factorization.cols && finite_values(factorization.householder_coefficients) &&
           valid_permutation(factorization.column_permutation, factorization.cols) && valid_qr_status(factorization.status) &&
           std::isfinite(factorization.relative_rank_tolerance) && factorization.relative_rank_tolerance > 0.0 &&
           factorization.numerical_rank <= factorization.cols && finite_matrix(factorization.factors) && rank_metadata_is_consistent(factorization);
}

void apply_householder_to_rhs(const PivotedQrFactorization& factorization, const std::size_t reflector, DenseMatrix& rhs)
{
    const double coefficient = factorization.householder_coefficients[reflector];
    if (coefficient == 0.0)
    {
        return;
    }

    for (std::size_t col = 0; col < rhs.cols(); ++col)
    {
        double projection = rhs(reflector, col);
        for (std::size_t row = reflector + 1; row < factorization.rows; ++row)
        {
            projection += factorization.factors(row, reflector) * rhs(row, col);
        }
        projection *= coefficient;

        rhs(reflector, col) -= projection;
        for (std::size_t row = reflector + 1; row < factorization.rows; ++row)
        {
            rhs(row, col) -= factorization.factors(row, reflector) * projection;
        }
    }
}

} // namespace

PivotedQrFactorization factorize_qr(DenseMatrix matrix, const double relative_rank_tolerance)
{
    const std::size_t rows = matrix.rows();
    const std::size_t cols = matrix.cols();
    if (rows == 0 || cols == 0)
    {
        return qr_failure_result(std::move(matrix), relative_rank_tolerance, "QR factorization matrix must not be empty.");
    }
    if (rows < cols)
    {
        return qr_failure_result(std::move(matrix), relative_rank_tolerance, "QR factorization requires rows >= columns.");
    }
    if (!std::isfinite(relative_rank_tolerance) || relative_rank_tolerance <= 0.0)
    {
        return qr_failure_result(std::move(matrix), relative_rank_tolerance, "QR relative rank tolerance is invalid.");
    }
    if (!finite_matrix(matrix))
    {
        return qr_failure_result(std::move(matrix), relative_rank_tolerance, "QR factorization matrix entries must be finite.");
    }

    PivotedQrFactorization result;
    result.rows = rows;
    result.cols = cols;
    result.factors = std::move(matrix);
    result.householder_coefficients.assign(cols, 0.0);
    result.column_permutation.resize(cols);
    std::iota(result.column_permutation.begin(), result.column_permutation.end(), 0);
    result.relative_rank_tolerance = relative_rank_tolerance;

    for (std::size_t step = 0; step < cols; ++step)
    {
        std::size_t pivot_col = step;
        double pivot_norm = stable_column_norm(result.factors, step, step);
        for (std::size_t col = step + 1; col < cols; ++col)
        {
            const double candidate_norm = stable_column_norm(result.factors, step, col);
            if (candidate_norm > pivot_norm)
            {
                pivot_norm = candidate_norm;
                pivot_col = col;
            }
        }

        if (pivot_col != step)
        {
            swap_columns(result.factors, step, pivot_col);
            std::swap(result.column_permutation[step], result.column_permutation[pivot_col]);
        }

        const double norm = stable_column_norm(result.factors, step, step);
        if (norm == 0.0)
        {
            result.householder_coefficients[step] = 0.0;
            continue;
        }

        const double x0 = result.factors(step, step);
        const double alpha = -std::copysign(norm, x0);
        const double coefficient = (alpha - x0) / alpha;
        const double scale = 1.0 / (x0 - alpha);
        for (std::size_t row = step + 1; row < rows; ++row)
        {
            result.factors(row, step) *= scale;
        }
        result.factors(step, step) = alpha;
        result.householder_coefficients[step] = coefficient;

        for (std::size_t col = step + 1; col < cols; ++col)
        {
            double projection = result.factors(step, col);
            for (std::size_t row = step + 1; row < rows; ++row)
            {
                projection += result.factors(row, step) * result.factors(row, col);
            }
            projection *= coefficient;

            result.factors(step, col) -= projection;
            for (std::size_t row = step + 1; row < rows; ++row)
            {
                result.factors(row, col) -= result.factors(row, step) * projection;
            }
        }
    }

    result.numerical_rank = numerical_rank_from_factor(result.factors, cols, relative_rank_tolerance);
    if (result.numerical_rank < cols)
    {
        result.status = QrStatus::rank_deficient;
        result.message = "QR factorization is rank deficient.";
    }
    return result;
}

QrSolveResult solve_qr(const PivotedQrFactorization& factorization, DenseMatrix rhs)
{
    const std::size_t rhs_cols = rhs.cols();
    if (!factorization_metadata_is_consistent(factorization))
    {
        return solve_failure_result(factorization.cols, rhs_cols, diagnostic_rank_from_metadata(factorization), QrStatus::invalid_input,
                                    "QR factorization metadata is inconsistent.");
    }
    if (rhs.rows() != factorization.rows)
    {
        return solve_failure_result(factorization.cols, rhs_cols, factorization.numerical_rank, QrStatus::invalid_input,
                                    "QR solve right-hand side row count is inconsistent.");
    }
    if (!finite_matrix(rhs))
    {
        return solve_failure_result(factorization.cols, rhs_cols, factorization.numerical_rank, QrStatus::invalid_input,
                                    "QR solve right-hand side entries must be finite.");
    }
    if (factorization.status == QrStatus::rank_deficient || factorization.numerical_rank < factorization.cols)
    {
        return solve_failure_result(factorization.cols, rhs_cols, factorization.numerical_rank, QrStatus::rank_deficient,
                                    "QR solve requires full column rank.");
    }
    if (factorization.status != QrStatus::success)
    {
        return solve_failure_result(factorization.cols, rhs_cols, factorization.numerical_rank, factorization.status, factorization.message);
    }

    for (std::size_t reflector = 0; reflector < factorization.cols; ++reflector)
    {
        apply_householder_to_rhs(factorization, reflector, rhs);
    }

    DenseMatrix pivoted_solution(factorization.cols, rhs_cols);
    for (std::size_t reverse_row = factorization.cols; reverse_row > 0; --reverse_row)
    {
        const std::size_t row = reverse_row - 1;
        for (std::size_t rhs_col = 0; rhs_col < rhs_cols; ++rhs_col)
        {
            double value = rhs(row, rhs_col);
            for (std::size_t col = row + 1; col < factorization.cols; ++col)
            {
                value -= factorization.factors(row, col) * pivoted_solution(col, rhs_col);
            }
            pivoted_solution(row, rhs_col) = value / factorization.factors(row, row);
        }
    }

    DenseMatrix solution(factorization.cols, rhs_cols);
    for (std::size_t pivoted_col = 0; pivoted_col < factorization.cols; ++pivoted_col)
    {
        const std::size_t original_col = factorization.column_permutation[pivoted_col];
        for (std::size_t rhs_col = 0; rhs_col < rhs_cols; ++rhs_col)
        {
            solution(original_col, rhs_col) = pivoted_solution(pivoted_col, rhs_col);
        }
    }

    return QrSolveResult{std::move(solution), factorization.numerical_rank, QrStatus::success, ""};
}

QrSolveResult solve_qr(DenseMatrix matrix, DenseMatrix rhs, const double relative_rank_tolerance)
{
    PivotedQrFactorization factorization = factorize_qr(std::move(matrix), relative_rank_tolerance);
    if (factorization.status == QrStatus::invalid_input)
    {
        return solve_failure_result(factorization.cols, rhs.cols(), factorization.numerical_rank, factorization.status, factorization.message);
    }
    return solve_qr(factorization, std::move(rhs));
}

} // namespace qmomkit::internal
