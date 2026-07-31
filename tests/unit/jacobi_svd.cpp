#include "internal/jacobi_svd.hpp"

#include "internal/householder_qr.hpp"
#include "internal/linear_algebra.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace
{

constexpr double default_tolerance = 1.0e-11;

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

double matrix_difference_norm(const qmomkit::internal::DenseMatrix& lhs, const qmomkit::internal::DenseMatrix& rhs)
{
    assert(lhs.rows() == rhs.rows());
    assert(lhs.cols() == rhs.cols());
    double norm = 0.0;
    for (std::size_t row = 0; row < lhs.rows(); ++row)
    {
        for (std::size_t col = 0; col < lhs.cols(); ++col)
        {
            norm = std::hypot(norm, lhs(row, col) - rhs(row, col));
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
    const double scale = std::max({1.0, frobenius_norm(actual), frobenius_norm(expected)});
    assert(matrix_difference_norm(actual, expected) <= tolerance * scale);
}

qmomkit::internal::DenseMatrix make_rhs(const qmomkit::internal::DenseMatrix& matrix, const qmomkit::internal::DenseMatrix& solution)
{
    return multiply(matrix, solution);
}

qmomkit::internal::DenseMatrix reconstruct(const qmomkit::internal::SvdFactorization& factorization)
{
    qmomkit::internal::DenseMatrix u_sigma(factorization.rows, factorization.cols);
    for (std::size_t col = 0; col < factorization.cols; ++col)
    {
        for (std::size_t row = 0; row < factorization.rows; ++row)
        {
            u_sigma(row, col) = factorization.left_singular_vectors(row, col) * factorization.singular_values[col];
        }
    }
    return multiply(u_sigma, transpose(factorization.right_singular_vectors));
}

void assert_reconstructs(const qmomkit::internal::DenseMatrix& matrix, const qmomkit::internal::SvdFactorization& factorization,
                         const double tolerance)
{
    assert(factorization.left_singular_vectors.rows() == matrix.rows());
    assert(factorization.left_singular_vectors.cols() == matrix.cols());
    assert(factorization.right_singular_vectors.rows() == matrix.cols());
    assert(factorization.right_singular_vectors.cols() == matrix.cols());
    assert(factorization.singular_values.size() == matrix.cols());
    assert_matrix_near(reconstruct(factorization), matrix, tolerance);
}

void assert_right_orthogonal(const qmomkit::internal::SvdFactorization& factorization, const double tolerance)
{
    assert_matrix_near(multiply(transpose(factorization.right_singular_vectors), factorization.right_singular_vectors),
                       qmomkit::internal::make_identity(factorization.cols), tolerance);
}

void assert_active_left_orthogonal(const qmomkit::internal::SvdFactorization& factorization, const double tolerance)
{
    qmomkit::internal::DenseMatrix active(factorization.rows, factorization.numerical_rank);
    for (std::size_t col = 0; col < factorization.numerical_rank; ++col)
    {
        for (std::size_t row = 0; row < factorization.rows; ++row)
        {
            active(row, col) = factorization.left_singular_vectors(row, col);
        }
    }
    assert_matrix_near(multiply(transpose(active), active), qmomkit::internal::make_identity(factorization.numerical_rank), tolerance);
}

void assert_finite_matrix(const qmomkit::internal::DenseMatrix& matrix)
{
    for (std::size_t row = 0; row < matrix.rows(); ++row)
    {
        for (std::size_t col = 0; col < matrix.cols(); ++col)
        {
            assert(std::isfinite(matrix(row, col)));
        }
    }
}

void assert_finite_factorization(const qmomkit::internal::SvdFactorization& factorization)
{
    assert_finite_matrix(factorization.left_singular_vectors);
    assert_finite_matrix(factorization.right_singular_vectors);
    for (const double value : factorization.singular_values)
    {
        assert(std::isfinite(value));
    }
}

double condition_number(const qmomkit::internal::SvdFactorization& factorization)
{
    if (factorization.numerical_rank < factorization.cols || factorization.singular_values.empty() || factorization.singular_values.back() == 0.0)
    {
        return std::numeric_limits<double>::infinity();
    }
    return factorization.singular_values.front() / factorization.singular_values.back();
}

void test_factorization_identity_and_orthogonality()
{
    const qmomkit::internal::DenseMatrix matrix = make_matrix({{1.0, 4.0, -1.0}, {2.0, 0.0, 3.0}, {-1.0, 1.0, 2.0}, {0.5, -2.0, 1.0}});

    const qmomkit::internal::SvdFactorization factorization = qmomkit::internal::factorize_svd(matrix);

    assert(factorization.status == qmomkit::internal::SvdStatus::success);
    assert(factorization.numerical_rank == matrix.cols());
    assert_reconstructs(matrix, factorization, default_tolerance);
    assert_right_orthogonal(factorization, default_tolerance);
    assert_active_left_orthogonal(factorization, default_tolerance);
}

void test_known_singular_values()
{
    const qmomkit::internal::SvdFactorization diagonal =
        qmomkit::internal::factorize_svd(make_matrix({{5.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, 0.5}}));
    assert(diagonal.status == qmomkit::internal::SvdStatus::success);
    assert_near(diagonal.singular_values[0], 5.0, 1.0e-13);
    assert_near(diagonal.singular_values[1], 2.0, 1.0e-13);
    assert_near(diagonal.singular_values[2], 0.5, 1.0e-13);

    const qmomkit::internal::SvdFactorization scaled_orthogonal =
        qmomkit::internal::factorize_svd(make_matrix({{3.0, 0.0}, {0.0, 0.0}, {0.0, 4.0}, {0.0, 0.0}}));
    assert(scaled_orthogonal.status == qmomkit::internal::SvdStatus::success);
    assert_near(scaled_orthogonal.singular_values[0], 4.0, 1.0e-13);
    assert_near(scaled_orthogonal.singular_values[1], 3.0, 1.0e-13);

    const qmomkit::internal::SvdFactorization two_by_two = qmomkit::internal::factorize_svd(make_matrix({{1.0, 1.0}, {0.0, 1.0}}));
    const double root_five = std::sqrt(5.0);
    assert(two_by_two.status == qmomkit::internal::SvdStatus::success);
    assert_near(two_by_two.singular_values[0], std::sqrt((3.0 + root_five) / 2.0), 1.0e-13);
    assert_near(two_by_two.singular_values[1], std::sqrt((3.0 - root_five) / 2.0), 1.0e-13);

    const qmomkit::internal::DenseMatrix rank_one_matrix = make_matrix({{1.0, 2.0}, {2.0, 4.0}, {3.0, 6.0}});
    const qmomkit::internal::SvdFactorization rank_one = qmomkit::internal::factorize_svd(rank_one_matrix);
    assert(rank_one.status == qmomkit::internal::SvdStatus::rank_deficient);
    assert(rank_one.numerical_rank == 1);
    assert_near(rank_one.singular_values[0], std::sqrt(70.0), 1.0e-12);
    assert_near(rank_one.singular_values[1], 0.0, 1.0e-12);
    assert_reconstructs(rank_one_matrix, rank_one, default_tolerance);
}

void test_ordering_and_repeated_singular_values()
{
    const qmomkit::internal::DenseMatrix ordered_matrix = make_matrix({{0.0, 3.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 2.0}});
    const qmomkit::internal::SvdFactorization ordered = qmomkit::internal::factorize_svd(ordered_matrix);
    assert(ordered.status == qmomkit::internal::SvdStatus::success);
    assert(ordered.singular_values[0] >= ordered.singular_values[1]);
    assert(ordered.singular_values[1] >= ordered.singular_values[2]);
    assert_near(ordered.singular_values[0], 3.0, 1.0e-13);
    assert_near(ordered.singular_values[1], 2.0, 1.0e-13);
    assert_near(ordered.singular_values[2], 1.0, 1.0e-13);
    assert_reconstructs(ordered_matrix, ordered, default_tolerance);

    const qmomkit::internal::DenseMatrix repeated_matrix = make_matrix({{1.0, 0.0}, {0.0, 1.0}, {0.0, 0.0}});
    const qmomkit::internal::SvdFactorization repeated = qmomkit::internal::factorize_svd(repeated_matrix);
    assert(repeated.status == qmomkit::internal::SvdStatus::success);
    assert_near(repeated.singular_values[0], 1.0, 1.0e-13);
    assert_near(repeated.singular_values[1], 1.0, 1.0e-13);
    assert_reconstructs(repeated_matrix, repeated, default_tolerance);
    assert_right_orthogonal(repeated, default_tolerance);
    assert_active_left_orthogonal(repeated, default_tolerance);
}

void test_square_solve_matches_dense()
{
    const qmomkit::internal::DenseMatrix matrix = make_matrix({{3.0, -1.0, 2.0}, {2.0, 4.0, -1.0}, {1.0, 0.5, 5.0}});
    const qmomkit::internal::DenseMatrix expected = make_matrix({{1.0, -2.0}, {0.5, 3.0}, {-1.5, 0.25}});
    const qmomkit::internal::DenseMatrix rhs = make_rhs(matrix, expected);

    const qmomkit::internal::SvdSolveResult svd = qmomkit::internal::solve_svd(matrix, rhs);
    const qmomkit::internal::LinearSolveResult dense = qmomkit::internal::solve_dense(matrix, rhs, 1.0e-14, "singular");

    assert(svd.status == qmomkit::internal::SvdStatus::success);
    assert(dense.success);
    assert_matrix_near(svd.solution, expected, default_tolerance);
    assert_matrix_near(svd.solution, dense.solution, default_tolerance);
}

void test_tall_exact_and_inconsistent_least_squares()
{
    const qmomkit::internal::DenseMatrix exact_matrix = make_matrix({{1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}, {2.0, -1.0}});
    const qmomkit::internal::DenseMatrix exact_expected = make_matrix({{1.25}, {-0.75}});
    const qmomkit::internal::DenseMatrix exact_rhs = make_rhs(exact_matrix, exact_expected);
    const qmomkit::internal::SvdSolveResult exact = qmomkit::internal::solve_svd(exact_matrix, exact_rhs);
    assert(exact.status == qmomkit::internal::SvdStatus::success);
    assert_matrix_near(exact.solution, exact_expected, default_tolerance);

    const qmomkit::internal::DenseMatrix inconsistent_matrix = make_matrix({{1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}});
    const qmomkit::internal::DenseMatrix inconsistent_rhs = make_matrix({{1.0}, {2.0}, {2.0}});
    const qmomkit::internal::DenseMatrix inconsistent_expected = make_matrix({{2.0 / 3.0}, {5.0 / 3.0}});
    const qmomkit::internal::SvdSolveResult inconsistent = qmomkit::internal::solve_svd(inconsistent_matrix, inconsistent_rhs);
    const qmomkit::internal::DenseMatrix residual = multiply(inconsistent_matrix, inconsistent.solution);
    qmomkit::internal::DenseMatrix adjusted_residual = residual;
    for (std::size_t row = 0; row < adjusted_residual.rows(); ++row)
    {
        adjusted_residual(row, 0) -= inconsistent_rhs(row, 0);
    }
    const qmomkit::internal::DenseMatrix normal_residual = multiply(transpose(inconsistent_matrix), adjusted_residual);

    assert(inconsistent.status == qmomkit::internal::SvdStatus::success);
    assert_matrix_near(inconsistent.solution, inconsistent_expected, default_tolerance);
    assert(frobenius_norm(normal_residual) <= 1.0e-12);
}

void test_multiple_rhs_and_qr_comparison()
{
    const qmomkit::internal::DenseMatrix square = make_matrix({{2.0, -1.0}, {1.0, 3.0}});
    const qmomkit::internal::DenseMatrix square_expected = make_matrix({{1.0, 0.25}, {-2.0, 3.0}});
    const qmomkit::internal::DenseMatrix square_rhs = make_rhs(square, square_expected);
    const qmomkit::internal::SvdSolveResult square_svd = qmomkit::internal::solve_svd(square, square_rhs);
    const qmomkit::internal::QrSolveResult square_qr = qmomkit::internal::solve_qr(square, square_rhs);
    assert(square_svd.status == qmomkit::internal::SvdStatus::success);
    assert(square_qr.status == qmomkit::internal::QrStatus::success);
    assert_matrix_near(square_svd.solution, square_expected, default_tolerance);
    assert_matrix_near(square_svd.solution, square_qr.solution, default_tolerance);

    const qmomkit::internal::DenseMatrix tall = make_matrix({{2.0, -1.0}, {1.0, 3.0}, {4.0, 0.5}});
    const qmomkit::internal::DenseMatrix tall_expected = make_matrix({{-3.0, 0.25}, {0.5, -1.5}});
    const qmomkit::internal::DenseMatrix tall_rhs = make_rhs(tall, tall_expected);
    const qmomkit::internal::SvdSolveResult tall_svd = qmomkit::internal::solve_svd(tall, tall_rhs);
    const qmomkit::internal::QrSolveResult tall_qr = qmomkit::internal::solve_qr(tall, tall_rhs);
    assert(tall_svd.status == qmomkit::internal::SvdStatus::success);
    assert(tall_qr.status == qmomkit::internal::QrStatus::success);
    assert_matrix_near(tall_svd.solution, tall_expected, default_tolerance);
    assert_matrix_near(tall_svd.solution, tall_qr.solution, default_tolerance);
}

void test_rank_deficiency_cases_refuse_solve()
{
    const std::vector<qmomkit::internal::DenseMatrix> matrices{
        make_matrix({{1.0, 2.0}, {3.0, 6.0}, {5.0, 10.0}}),
        make_matrix({{1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}}),
        make_matrix({{1.0, 2.0}, {2.0, 4.0}, {3.0, 6.0}}),
        qmomkit::internal::DenseMatrix(3, 2),
    };

    for (const qmomkit::internal::DenseMatrix& matrix : matrices)
    {
        const qmomkit::internal::SvdFactorization factorization = qmomkit::internal::factorize_svd(matrix);
        const qmomkit::internal::SvdSolveResult solve = qmomkit::internal::solve_svd(factorization, qmomkit::internal::DenseMatrix(matrix.rows(), 1));
        assert(factorization.status == qmomkit::internal::SvdStatus::rank_deficient);
        assert(factorization.numerical_rank < matrix.cols());
        assert(factorization.left_singular_vectors.rows() == matrix.rows());
        assert(factorization.right_singular_vectors.cols() == matrix.cols());
        assert(factorization.singular_values.size() == matrix.cols());
        assert(solve.status == qmomkit::internal::SvdStatus::rank_deficient);
        assert(solve.solution.rows() == matrix.cols());
        assert(solve.solution.cols() == 1);
    }
}

void test_near_rank_deficiency_and_condition_numbers()
{
    const qmomkit::internal::SvdOptions options{1.0e-8, 0.0, 0};
    const qmomkit::internal::SvdFactorization full = qmomkit::internal::factorize_svd(make_matrix({{1.0, 0.0}, {0.0, 1.0e-6}}), options);
    const qmomkit::internal::SvdFactorization deficient = qmomkit::internal::factorize_svd(make_matrix({{1.0, 0.0}, {0.0, 1.0e-10}}), options);
    assert(full.status == qmomkit::internal::SvdStatus::success);
    assert(full.numerical_rank == 2);
    assert(deficient.status == qmomkit::internal::SvdStatus::rank_deficient);
    assert(deficient.numerical_rank == 1);

    const qmomkit::internal::SvdFactorization conditioned = qmomkit::internal::factorize_svd(make_matrix({{6.0, 0.0}, {0.0, 2.0}}));
    const qmomkit::internal::SvdFactorization scaled_conditioned = qmomkit::internal::factorize_svd(make_matrix({{600.0, 0.0}, {0.0, 200.0}}));
    assert_near(condition_number(conditioned), 3.0, 1.0e-13);
    assert_near(condition_number(scaled_conditioned), 3.0, 1.0e-13);
    assert(condition_number(deficient) == std::numeric_limits<double>::infinity());
}

void test_uniform_scaling_large_and_small()
{
    const qmomkit::internal::DenseMatrix matrix = make_matrix({{3.0, 1.0}, {0.0, 2.0}, {-1.0, 4.0}});
    qmomkit::internal::DenseMatrix large = matrix;
    qmomkit::internal::DenseMatrix small = matrix;
    for (std::size_t row = 0; row < matrix.rows(); ++row)
    {
        for (std::size_t col = 0; col < matrix.cols(); ++col)
        {
            large(row, col) *= 1.0e100;
            small(row, col) *= 1.0e-100;
        }
    }

    const qmomkit::internal::SvdFactorization reference = qmomkit::internal::factorize_svd(matrix);
    const qmomkit::internal::SvdFactorization large_factorization = qmomkit::internal::factorize_svd(large);
    const qmomkit::internal::SvdFactorization small_factorization = qmomkit::internal::factorize_svd(small);

    assert(reference.status == qmomkit::internal::SvdStatus::success);
    assert(large_factorization.status == qmomkit::internal::SvdStatus::success);
    assert(small_factorization.status == qmomkit::internal::SvdStatus::success);
    assert_finite_factorization(large_factorization);
    assert_finite_factorization(small_factorization);
    assert(reference.numerical_rank == large_factorization.numerical_rank);
    assert(reference.numerical_rank == small_factorization.numerical_rank);
    for (std::size_t i = 0; i < reference.singular_values.size(); ++i)
    {
        assert_near(large_factorization.singular_values[i] / 1.0e100, reference.singular_values[i], 5.0e-12);
        assert_near(small_factorization.singular_values[i] / 1.0e-100, reference.singular_values[i], 5.0e-12);
    }
    assert_reconstructs(large, large_factorization, 1.0e-11);
    assert_reconstructs(small, small_factorization, 1.0e-11);
}

void test_nonconvergence_is_reported_and_refused()
{
    const qmomkit::internal::DenseMatrix matrix = make_matrix({{1.0, 1.0}, {0.25, 2.0}, {-1.0, 0.5}});
    const qmomkit::internal::SvdOptions options{1.0e-14, 1.0e-16, 1};

    const qmomkit::internal::SvdFactorization factorization = qmomkit::internal::factorize_svd(matrix, options);
    const qmomkit::internal::SvdSolveResult solve = qmomkit::internal::solve_svd(factorization, qmomkit::internal::DenseMatrix(matrix.rows(), 1));

    assert(factorization.status == qmomkit::internal::SvdStatus::no_convergence);
    assert(factorization.sweeps == 1);
    assert(factorization.left_singular_vectors.rows() == matrix.rows());
    assert(factorization.singular_values.size() == matrix.cols());
    assert(solve.status == qmomkit::internal::SvdStatus::no_convergence);
}

void test_invalid_factorization_input()
{
    qmomkit::internal::SvdFactorization factorization = qmomkit::internal::factorize_svd(qmomkit::internal::DenseMatrix(0, 1));
    assert(factorization.status == qmomkit::internal::SvdStatus::invalid_input);

    factorization = qmomkit::internal::factorize_svd(qmomkit::internal::DenseMatrix(2, 3));
    assert(factorization.status == qmomkit::internal::SvdStatus::invalid_input);

    factorization = qmomkit::internal::factorize_svd(make_matrix({{std::numeric_limits<double>::infinity()}}));
    assert(factorization.status == qmomkit::internal::SvdStatus::invalid_input);

    factorization = qmomkit::internal::factorize_svd(make_matrix({{1.0}}), {0.0, 0.0, 0});
    assert(factorization.status == qmomkit::internal::SvdStatus::invalid_input);

    factorization = qmomkit::internal::factorize_svd(make_matrix({{1.0}}), {1.0e-14, -1.0, 0});
    assert(factorization.status == qmomkit::internal::SvdStatus::invalid_input);
}

void test_invalid_solve_input_and_malformed_metadata()
{
    const qmomkit::internal::DenseMatrix matrix = make_matrix({{2.0, -1.0}, {1.0, 3.0}, {4.0, 0.5}});
    const qmomkit::internal::SvdFactorization factorization = qmomkit::internal::factorize_svd(matrix);

    qmomkit::internal::SvdSolveResult solve = qmomkit::internal::solve_svd(factorization, qmomkit::internal::DenseMatrix(2, 1));
    assert(solve.status == qmomkit::internal::SvdStatus::invalid_input);

    solve = qmomkit::internal::solve_svd(factorization, make_matrix({{1.0}, {std::numeric_limits<double>::quiet_NaN()}, {2.0}}));
    assert(solve.status == qmomkit::internal::SvdStatus::invalid_input);

    qmomkit::internal::SvdFactorization malformed = factorization;
    malformed.left_singular_vectors = qmomkit::internal::DenseMatrix(2, 2);
    solve = qmomkit::internal::solve_svd(malformed, qmomkit::internal::DenseMatrix(3, 1));
    assert(solve.status == qmomkit::internal::SvdStatus::invalid_input);

    malformed = factorization;
    malformed.right_singular_vectors = qmomkit::internal::DenseMatrix(3, 3);
    solve = qmomkit::internal::solve_svd(malformed, qmomkit::internal::DenseMatrix(3, 1));
    assert(solve.status == qmomkit::internal::SvdStatus::invalid_input);

    malformed = factorization;
    malformed.singular_values.pop_back();
    solve = qmomkit::internal::solve_svd(malformed, qmomkit::internal::DenseMatrix(3, 1));
    assert(solve.status == qmomkit::internal::SvdStatus::invalid_input);

    malformed = factorization;
    malformed.singular_values[0] = -1.0;
    solve = qmomkit::internal::solve_svd(malformed, qmomkit::internal::DenseMatrix(3, 1));
    assert(solve.status == qmomkit::internal::SvdStatus::invalid_input);

    malformed = factorization;
    std::swap(malformed.singular_values[0], malformed.singular_values[1]);
    solve = qmomkit::internal::solve_svd(malformed, qmomkit::internal::DenseMatrix(3, 1));
    assert(solve.status == qmomkit::internal::SvdStatus::invalid_input);

    malformed = factorization;
    malformed.numerical_rank = 0;
    solve = qmomkit::internal::solve_svd(malformed, qmomkit::internal::DenseMatrix(3, 1));
    assert(solve.status == qmomkit::internal::SvdStatus::invalid_input);

    malformed = factorization;
    malformed.singular_values[1] = 0.0;
    malformed.numerical_rank = malformed.cols;
    malformed.status = qmomkit::internal::SvdStatus::success;
    solve = qmomkit::internal::solve_svd(malformed, qmomkit::internal::DenseMatrix(3, 1));
    assert(solve.status == qmomkit::internal::SvdStatus::invalid_input);
}

void test_zero_column_rhs_and_reuse_immutability()
{
    const qmomkit::internal::DenseMatrix matrix = make_matrix({{2.0, -1.0}, {1.0, 3.0}, {4.0, 0.5}});
    const qmomkit::internal::SvdFactorization factorization = qmomkit::internal::factorize_svd(matrix);
    const qmomkit::internal::DenseMatrix left_before = factorization.left_singular_vectors;
    const qmomkit::internal::DenseMatrix right_before = factorization.right_singular_vectors;
    const std::vector<double> singular_values_before = factorization.singular_values;

    const qmomkit::internal::SvdSolveResult zero_rhs = qmomkit::internal::solve_svd(factorization, qmomkit::internal::DenseMatrix(matrix.rows(), 0));
    assert(zero_rhs.status == qmomkit::internal::SvdStatus::success);
    assert(zero_rhs.solution.rows() == matrix.cols());
    assert(zero_rhs.solution.cols() == 0);

    const qmomkit::internal::DenseMatrix expected_first = make_matrix({{1.0}, {2.0}});
    const qmomkit::internal::DenseMatrix expected_second = make_matrix({{-3.0, 0.25}, {0.5, -1.5}});
    const qmomkit::internal::SvdSolveResult first = qmomkit::internal::solve_svd(factorization, make_rhs(matrix, expected_first));
    const qmomkit::internal::SvdSolveResult second = qmomkit::internal::solve_svd(factorization, make_rhs(matrix, expected_second));

    assert(first.status == qmomkit::internal::SvdStatus::success);
    assert(second.status == qmomkit::internal::SvdStatus::success);
    assert_matrix_near(first.solution, expected_first, default_tolerance);
    assert_matrix_near(second.solution, expected_second, default_tolerance);
    assert_matrix_near(factorization.left_singular_vectors, left_before, 0.0);
    assert_matrix_near(factorization.right_singular_vectors, right_before, 0.0);
    assert(factorization.singular_values == singular_values_before);
}

} // namespace

int main()
{
    test_factorization_identity_and_orthogonality();
    test_known_singular_values();
    test_ordering_and_repeated_singular_values();
    test_square_solve_matches_dense();
    test_tall_exact_and_inconsistent_least_squares();
    test_multiple_rhs_and_qr_comparison();
    test_rank_deficiency_cases_refuse_solve();
    test_near_rank_deficiency_and_condition_numbers();
    test_uniform_scaling_large_and_small();
    test_nonconvergence_is_reported_and_refused();
    test_invalid_factorization_input();
    test_invalid_solve_input_and_malformed_metadata();
    test_zero_column_rhs_and_reuse_immutability();

    return 0;
}
