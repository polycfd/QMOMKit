#include "internal/householder_qr.hpp"

#include "internal/linear_algebra.hpp"
#include "internal/vandermonde.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace
{

constexpr double default_tolerance = 1.0e-12;

qmomkit::internal::DenseMatrix make_matrix(const std::vector<std::vector<double>>& rows)
{
    const std::size_t row_count = rows.size();
    const std::size_t col_count = rows.empty() ? 0 : rows.front().size();
    qmomkit::internal::DenseMatrix matrix(row_count, col_count);
    for (std::size_t row = 0; row < row_count; ++row)
    {
        assert(rows[row].size() == col_count);
        for (std::size_t col = 0; col < col_count; ++col)
        {
            matrix(row, col) = rows[row][col];
        }
    }
    return matrix;
}

qmomkit::internal::DenseMatrix multiply(const qmomkit::internal::DenseMatrix& lhs, const qmomkit::internal::DenseMatrix& rhs)
{
    assert(lhs.cols() == rhs.rows());
    qmomkit::internal::DenseMatrix product(lhs.rows(), rhs.cols());
    for (std::size_t row = 0; row < lhs.rows(); ++row)
    {
        for (std::size_t col = 0; col < rhs.cols(); ++col)
        {
            double value = 0.0;
            for (std::size_t k = 0; k < lhs.cols(); ++k)
            {
                value += lhs(row, k) * rhs(k, col);
            }
            product(row, col) = value;
        }
    }
    return product;
}

qmomkit::internal::DenseMatrix transpose(const qmomkit::internal::DenseMatrix& matrix)
{
    qmomkit::internal::DenseMatrix transposed(matrix.cols(), matrix.rows());
    for (std::size_t row = 0; row < matrix.rows(); ++row)
    {
        for (std::size_t col = 0; col < matrix.cols(); ++col)
        {
            transposed(col, row) = matrix(row, col);
        }
    }
    return transposed;
}

qmomkit::internal::DenseMatrix residual_matrix(const qmomkit::internal::DenseMatrix& matrix, const qmomkit::internal::DenseMatrix& solution,
                                               const qmomkit::internal::DenseMatrix& rhs)
{
    qmomkit::internal::DenseMatrix residual = multiply(matrix, solution);
    for (std::size_t row = 0; row < residual.rows(); ++row)
    {
        for (std::size_t col = 0; col < residual.cols(); ++col)
        {
            residual(row, col) -= rhs(row, col);
        }
    }
    return residual;
}

double frobenius_norm(const qmomkit::internal::DenseMatrix& matrix)
{
    double norm = 0.0;
    for (std::size_t row = 0; row < matrix.rows(); ++row)
    {
        for (std::size_t col = 0; col < matrix.cols(); ++col)
        {
            norm = std::hypot(norm, matrix(row, col));
        }
    }
    return norm;
}

void assert_near(const double actual, const double expected, const double tolerance)
{
    assert(std::abs(actual - expected) <= tolerance);
}

void assert_matrix_near(const qmomkit::internal::DenseMatrix& actual, const qmomkit::internal::DenseMatrix& expected, const double tolerance)
{
    assert(actual.rows() == expected.rows());
    assert(actual.cols() == expected.cols());
    for (std::size_t row = 0; row < actual.rows(); ++row)
    {
        for (std::size_t col = 0; col < actual.cols(); ++col)
        {
            assert_near(actual(row, col), expected(row, col), tolerance);
        }
    }
}

void assert_scaled_residual(const qmomkit::internal::DenseMatrix& matrix, const qmomkit::internal::DenseMatrix& solution,
                            const qmomkit::internal::DenseMatrix& rhs, const double tolerance)
{
    const double residual = qmomkit::internal::residual_norm(matrix, solution, rhs);
    const double scale = std::max(1.0, qmomkit::internal::matrix_infinity_norm(matrix) * qmomkit::internal::matrix_infinity_norm(solution) +
                                           qmomkit::internal::matrix_infinity_norm(rhs));
    assert(residual <= tolerance * scale);
}

qmomkit::internal::DenseMatrix make_rhs(const qmomkit::internal::DenseMatrix& matrix, const qmomkit::internal::DenseMatrix& solution)
{
    return multiply(matrix, solution);
}

void apply_reflector_left(const qmomkit::internal::PivotedQrFactorization& factorization, const std::size_t reflector,
                          qmomkit::internal::DenseMatrix& matrix)
{
    const double coefficient = factorization.householder_coefficients[reflector];
    if (coefficient == 0.0)
    {
        return;
    }

    for (std::size_t col = 0; col < matrix.cols(); ++col)
    {
        double projection = matrix(reflector, col);
        for (std::size_t row = reflector + 1; row < factorization.rows; ++row)
        {
            projection += factorization.factors(row, reflector) * matrix(row, col);
        }
        projection *= coefficient;

        matrix(reflector, col) -= projection;
        for (std::size_t row = reflector + 1; row < factorization.rows; ++row)
        {
            matrix(row, col) -= factorization.factors(row, reflector) * projection;
        }
    }
}

qmomkit::internal::DenseMatrix explicit_q(const qmomkit::internal::PivotedQrFactorization& factorization)
{
    qmomkit::internal::DenseMatrix q = qmomkit::internal::make_identity(factorization.rows);
    for (std::size_t reverse_reflector = factorization.cols; reverse_reflector > 0; --reverse_reflector)
    {
        apply_reflector_left(factorization, reverse_reflector - 1, q);
    }
    return q;
}

qmomkit::internal::DenseMatrix explicit_r(const qmomkit::internal::PivotedQrFactorization& factorization)
{
    qmomkit::internal::DenseMatrix r(factorization.rows, factorization.cols);
    for (std::size_t row = 0; row < factorization.cols; ++row)
    {
        for (std::size_t col = row; col < factorization.cols; ++col)
        {
            r(row, col) = factorization.factors(row, col);
        }
    }
    return r;
}

qmomkit::internal::DenseMatrix pivoted_columns(const qmomkit::internal::DenseMatrix& matrix, const std::vector<std::size_t>& column_permutation)
{
    qmomkit::internal::DenseMatrix pivoted(matrix.rows(), matrix.cols());
    for (std::size_t pivoted_col = 0; pivoted_col < matrix.cols(); ++pivoted_col)
    {
        const std::size_t original_col = column_permutation[pivoted_col];
        for (std::size_t row = 0; row < matrix.rows(); ++row)
        {
            pivoted(row, pivoted_col) = matrix(row, original_col);
        }
    }
    return pivoted;
}

qmomkit::internal::DenseMatrix make_patterned_solution(const std::size_t rows, const std::size_t cols)
{
    qmomkit::internal::DenseMatrix solution(rows, cols);
    for (std::size_t row = 0; row < rows; ++row)
    {
        for (std::size_t col = 0; col < cols; ++col)
        {
            const double sign = (row + col) % 2 == 0 ? 1.0 : -1.0;
            solution(row, col) = sign * (0.25 + 0.125 * static_cast<double>(row + 1) + 0.0625 * static_cast<double>(col + 1));
        }
    }
    return solution;
}

void test_square_full_rank_solve()
{
    const qmomkit::internal::DenseMatrix matrix = make_matrix({{3.0, -1.0, 2.0}, {2.0, 4.0, -1.0}, {1.0, 0.5, 5.0}});
    const qmomkit::internal::DenseMatrix expected = make_matrix({{1.0, -2.0}, {0.5, 3.0}, {-1.5, 0.25}});
    const qmomkit::internal::DenseMatrix rhs = make_rhs(matrix, expected);

    const qmomkit::internal::QrSolveResult result = qmomkit::internal::solve_qr(matrix, rhs);

    assert(result.status == qmomkit::internal::QrStatus::success);
    assert(result.numerical_rank == matrix.cols());
    assert_matrix_near(result.solution, expected, 1.0e-12);
    assert_scaled_residual(matrix, result.solution, rhs, 1.0e-13);
}

void test_tall_exact_least_squares_system()
{
    const qmomkit::internal::DenseMatrix matrix = make_matrix({{1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}, {2.0, -1.0}});
    const qmomkit::internal::DenseMatrix expected = make_matrix({{1.25}, {-0.75}});
    const qmomkit::internal::DenseMatrix rhs = make_rhs(matrix, expected);

    const qmomkit::internal::QrSolveResult result = qmomkit::internal::solve_qr(matrix, rhs);

    assert(result.status == qmomkit::internal::QrStatus::success);
    assert(result.numerical_rank == matrix.cols());
    assert_matrix_near(result.solution, expected, 1.0e-12);
    assert_scaled_residual(matrix, result.solution, rhs, 1.0e-13);
}

void test_tall_inconsistent_least_squares_system()
{
    const qmomkit::internal::DenseMatrix matrix = make_matrix({{1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}});
    const qmomkit::internal::DenseMatrix rhs = make_matrix({{1.0}, {2.0}, {2.0}});
    const qmomkit::internal::DenseMatrix expected = make_matrix({{2.0 / 3.0}, {5.0 / 3.0}});

    const qmomkit::internal::QrSolveResult result = qmomkit::internal::solve_qr(matrix, rhs);
    const qmomkit::internal::DenseMatrix residual = residual_matrix(matrix, result.solution, rhs);
    const qmomkit::internal::DenseMatrix normal_residual = multiply(transpose(matrix), residual);

    assert(result.status == qmomkit::internal::QrStatus::success);
    assert_matrix_near(result.solution, expected, 1.0e-12);
    assert_near(frobenius_norm(residual), std::sqrt(1.0 / 3.0), 1.0e-12);
    assert(frobenius_norm(normal_residual) <= 1.0e-13);
}

void test_nontrivial_column_pivot_restores_solution_order()
{
    const qmomkit::internal::DenseMatrix matrix = make_matrix({{1.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 2.0}, {0.5, 0.0, 0.0}});
    const qmomkit::internal::DenseMatrix expected = make_matrix({{4.0}, {-2.0}, {3.0}});
    const qmomkit::internal::DenseMatrix rhs = make_rhs(matrix, expected);

    const qmomkit::internal::PivotedQrFactorization factorization = qmomkit::internal::factorize_qr(matrix);
    const qmomkit::internal::QrSolveResult result = qmomkit::internal::solve_qr(factorization, rhs);

    assert(factorization.status == qmomkit::internal::QrStatus::success);
    assert(factorization.column_permutation[0] == 1);
    assert(result.status == qmomkit::internal::QrStatus::success);
    assert_matrix_near(result.solution, expected, 1.0e-12);
}

void test_reusable_factorization()
{
    const qmomkit::internal::DenseMatrix matrix = make_matrix({{2.0, -1.0}, {1.0, 3.0}, {4.0, 0.5}});
    const qmomkit::internal::DenseMatrix expected_first = make_matrix({{1.0}, {2.0}});
    const qmomkit::internal::DenseMatrix expected_second = make_matrix({{-3.0, 0.25}, {0.5, -1.5}});
    const qmomkit::internal::DenseMatrix rhs_first = make_rhs(matrix, expected_first);
    const qmomkit::internal::DenseMatrix rhs_second = make_rhs(matrix, expected_second);

    const qmomkit::internal::PivotedQrFactorization factorization = qmomkit::internal::factorize_qr(matrix);
    const qmomkit::internal::DenseMatrix factors_before = factorization.factors;
    const std::vector<double> coefficients_before = factorization.householder_coefficients;
    const std::vector<std::size_t> permutation_before = factorization.column_permutation;

    const qmomkit::internal::QrSolveResult first = qmomkit::internal::solve_qr(factorization, rhs_first);
    const qmomkit::internal::QrSolveResult second = qmomkit::internal::solve_qr(factorization, rhs_second);

    assert(first.status == qmomkit::internal::QrStatus::success);
    assert(second.status == qmomkit::internal::QrStatus::success);
    assert_matrix_near(first.solution, expected_first, 1.0e-12);
    assert_matrix_near(second.solution, expected_second, 1.0e-12);
    assert_matrix_near(factorization.factors, factors_before, 0.0);
    assert(factorization.householder_coefficients == coefficients_before);
    assert(factorization.column_permutation == permutation_before);
}

void test_rank_deficient_matrix_refuses_solve()
{
    const qmomkit::internal::DenseMatrix matrix = make_matrix({{1.0, 2.0, 1.0}, {3.0, 4.0, 3.0}, {5.0, 6.0, 5.0}, {7.0, 8.0, 7.0}});
    const qmomkit::internal::DenseMatrix rhs = make_matrix({{1.0}, {2.0}, {3.0}, {4.0}});

    const qmomkit::internal::PivotedQrFactorization factorization = qmomkit::internal::factorize_qr(matrix);
    const qmomkit::internal::QrSolveResult result = qmomkit::internal::solve_qr(factorization, rhs);

    assert(factorization.status == qmomkit::internal::QrStatus::rank_deficient);
    assert(factorization.numerical_rank == 2);
    assert(factorization.rows == matrix.rows());
    assert(factorization.cols == matrix.cols());
    assert(factorization.factors.rows() == matrix.rows());
    assert(factorization.factors.cols() == matrix.cols());
    assert(factorization.householder_coefficients.size() == matrix.cols());
    assert(factorization.column_permutation.size() == matrix.cols());
    assert(result.status == qmomkit::internal::QrStatus::rank_deficient);
    assert(result.numerical_rank == 2);
    assert(result.solution.rows() == matrix.cols());
    assert(result.solution.cols() == rhs.cols());
}

void test_zero_matrix_reports_rank_zero()
{
    const qmomkit::internal::DenseMatrix matrix(3, 2);

    const qmomkit::internal::PivotedQrFactorization factorization = qmomkit::internal::factorize_qr(matrix);

    assert(factorization.status == qmomkit::internal::QrStatus::rank_deficient);
    assert(factorization.numerical_rank == 0);
}

void test_nearly_dependent_columns_follow_relative_tolerance()
{
    const qmomkit::internal::DenseMatrix full_rank = make_matrix({{1.0, 1.0}, {0.0, 1.0e-6}, {0.0, 0.0}});
    const qmomkit::internal::DenseMatrix rank_deficient = make_matrix({{1.0, 1.0}, {0.0, 1.0e-10}, {0.0, 0.0}});

    const qmomkit::internal::PivotedQrFactorization full = qmomkit::internal::factorize_qr(full_rank, 1.0e-8);
    const qmomkit::internal::PivotedQrFactorization deficient = qmomkit::internal::factorize_qr(rank_deficient, 1.0e-8);

    assert(full.status == qmomkit::internal::QrStatus::success);
    assert(full.numerical_rank == 2);
    assert(deficient.status == qmomkit::internal::QrStatus::rank_deficient);
    assert(deficient.numerical_rank == 1);
}

void test_rank_decision_is_scaling_invariant()
{
    const qmomkit::internal::DenseMatrix matrix = make_matrix({{1.0, 1.0}, {0.0, 1.0e-7}, {0.0, 0.0}});
    qmomkit::internal::DenseMatrix small = matrix;
    qmomkit::internal::DenseMatrix large = matrix;
    for (std::size_t row = 0; row < matrix.rows(); ++row)
    {
        for (std::size_t col = 0; col < matrix.cols(); ++col)
        {
            small(row, col) *= 1.0e-80;
            large(row, col) *= 1.0e80;
        }
    }

    const qmomkit::internal::PivotedQrFactorization unscaled = qmomkit::internal::factorize_qr(matrix, 1.0e-8);
    const qmomkit::internal::PivotedQrFactorization small_scaled = qmomkit::internal::factorize_qr(small, 1.0e-8);
    const qmomkit::internal::PivotedQrFactorization large_scaled = qmomkit::internal::factorize_qr(large, 1.0e-8);

    assert(unscaled.numerical_rank == small_scaled.numerical_rank);
    assert(unscaled.numerical_rank == large_scaled.numerical_rank);
    assert(unscaled.status == small_scaled.status);
    assert(unscaled.status == large_scaled.status);
}

void test_invalid_input_reports_status()
{
    qmomkit::internal::DenseMatrix wide(2, 3);
    qmomkit::internal::PivotedQrFactorization factorization = qmomkit::internal::factorize_qr(wide);
    assert(factorization.status == qmomkit::internal::QrStatus::invalid_input);

    factorization = qmomkit::internal::factorize_qr(make_matrix({{1.0}}), 0.0);
    assert(factorization.status == qmomkit::internal::QrStatus::invalid_input);

    factorization = qmomkit::internal::factorize_qr(make_matrix({{std::numeric_limits<double>::infinity()}}));
    assert(factorization.status == qmomkit::internal::QrStatus::invalid_input);

    const qmomkit::internal::DenseMatrix matrix = make_matrix({{1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}});
    factorization = qmomkit::internal::factorize_qr(matrix);

    qmomkit::internal::QrSolveResult solve = qmomkit::internal::solve_qr(factorization, qmomkit::internal::DenseMatrix(2, 1));
    assert(solve.status == qmomkit::internal::QrStatus::invalid_input);

    solve = qmomkit::internal::solve_qr(factorization, make_matrix({{1.0}, {std::numeric_limits<double>::quiet_NaN()}, {2.0}}));
    assert(solve.status == qmomkit::internal::QrStatus::invalid_input);

    qmomkit::internal::PivotedQrFactorization inconsistent = factorization;
    inconsistent.column_permutation = {0, 0};
    solve = qmomkit::internal::solve_qr(inconsistent, make_matrix({{1.0}, {2.0}, {3.0}}));
    assert(solve.status == qmomkit::internal::QrStatus::invalid_input);
}

void test_malformed_full_rank_factorization_rejects_zero_diagonal()
{
    const qmomkit::internal::DenseMatrix matrix = make_matrix({{2.0, -1.0}, {1.0, 3.0}, {4.0, 0.5}});
    const qmomkit::internal::DenseMatrix expected = make_matrix({{1.0}, {-2.0}});
    const qmomkit::internal::DenseMatrix rhs = make_rhs(matrix, expected);
    qmomkit::internal::PivotedQrFactorization factorization = qmomkit::internal::factorize_qr(matrix);
    qmomkit::internal::PivotedQrFactorization malformed = factorization;

    malformed.status = qmomkit::internal::QrStatus::success;
    malformed.numerical_rank = malformed.cols;
    malformed.factors(malformed.cols - 1, malformed.cols - 1) = 0.0;
    const qmomkit::internal::DenseMatrix malformed_factors_before = malformed.factors;
    const std::vector<double> malformed_coefficients_before = malformed.householder_coefficients;
    const std::vector<std::size_t> malformed_permutation_before = malformed.column_permutation;

    const qmomkit::internal::QrSolveResult result = qmomkit::internal::solve_qr(malformed, rhs);

    assert(result.status == qmomkit::internal::QrStatus::invalid_input);
    assert(result.numerical_rank < malformed.cols);
    assert(result.solution.rows() == malformed.cols);
    assert(result.solution.cols() == rhs.cols());
    for (std::size_t row = 0; row < result.solution.rows(); ++row)
    {
        for (std::size_t col = 0; col < result.solution.cols(); ++col)
        {
            assert(std::isfinite(result.solution(row, col)));
            assert(result.solution(row, col) == 0.0);
        }
    }
    assert_matrix_near(malformed.factors, malformed_factors_before, 0.0);
    assert(malformed.householder_coefficients == malformed_coefficients_before);
    assert(malformed.column_permutation == malformed_permutation_before);
}

void test_factorization_identity()
{
    const qmomkit::internal::DenseMatrix matrix = make_matrix({{1.0, 4.0, -1.0}, {2.0, 0.0, 3.0}, {-1.0, 1.0, 2.0}, {0.5, -2.0, 1.0}});

    const qmomkit::internal::PivotedQrFactorization factorization = qmomkit::internal::factorize_qr(matrix);
    const qmomkit::internal::DenseMatrix q = explicit_q(factorization);
    const qmomkit::internal::DenseMatrix r = explicit_r(factorization);
    const qmomkit::internal::DenseMatrix ap = pivoted_columns(matrix, factorization.column_permutation);
    const qmomkit::internal::DenseMatrix qr = multiply(q, r);
    const qmomkit::internal::DenseMatrix qtq = multiply(transpose(q), q);
    const qmomkit::internal::DenseMatrix identity = qmomkit::internal::make_identity(q.rows());

    assert(factorization.status == qmomkit::internal::QrStatus::success);
    assert_matrix_near(qr, ap, 1.0e-12);
    assert_matrix_near(qtq, identity, 1.0e-12);
}

void assert_vandermonde_case_matches_references(const std::vector<double>& abscissas, const std::vector<double>& weights,
                                                const qmomkit::internal::DenseMatrix& expected_solution, const double tolerance)
{
    const qmomkit::internal::DenseMatrix system = qmomkit::internal::make_conditional_vandermonde_matrix(abscissas, weights);
    const qmomkit::internal::DenseMatrix rhs = make_rhs(system, expected_solution);
    const qmomkit::internal::LinearSolveResult dense = qmomkit::internal::solve_dense(system, rhs, 1.0e-14, "singular");
    const qmomkit::internal::VandermondeSolveOptions options{qmomkit::VandermondeSolver::bjorck_pereyra, 1.0e-14, false, 1.0e-10};
    const qmomkit::internal::VandermondeSolveResult bjorck_pereyra =
        qmomkit::internal::solve_conditional_vandermonde_system(abscissas, weights, rhs, options, "singular");
    const qmomkit::internal::QrSolveResult qr = qmomkit::internal::solve_qr(system, rhs);

    assert(dense.success);
    assert(bjorck_pereyra.success);
    assert(qr.status == qmomkit::internal::QrStatus::success);
    assert_matrix_near(qr.solution, expected_solution, tolerance);
    assert_matrix_near(qr.solution, dense.solution, tolerance);
    assert_matrix_near(qr.solution, bjorck_pereyra.solution, tolerance);
    assert_scaled_residual(system, qr.solution, rhs, tolerance);
}

void test_vandermonde_shaped_systems()
{
    assert_vandermonde_case_matches_references({0.25, 1.0, 2.0}, {0.5, 1.0, 1.5}, make_patterned_solution(3, 2), 1.0e-12);
    assert_vandermonde_case_matches_references({-2.0, 0.25, 1.5, 4.0}, {0.5, 1.25, -0.75, 2.0}, make_patterned_solution(4, 3), 1.0e-10);
    assert_vandermonde_case_matches_references({1.0, 1.0 + 1.0e-5, 1.0 + 3.0e-5}, {0.75, -1.5, 2.25}, make_patterned_solution(3, 4), 5.0e-5);
}

} // namespace

int main()
{
    test_square_full_rank_solve();
    test_tall_exact_least_squares_system();
    test_tall_inconsistent_least_squares_system();
    test_nontrivial_column_pivot_restores_solution_order();
    test_reusable_factorization();
    test_rank_deficient_matrix_refuses_solve();
    test_zero_matrix_reports_rank_zero();
    test_nearly_dependent_columns_follow_relative_tolerance();
    test_rank_decision_is_scaling_invariant();
    test_invalid_input_reports_status();
    test_malformed_full_rank_factorization_rejects_zero_diagonal();
    test_factorization_identity();
    test_vandermonde_shaped_systems();

    return 0;
}
