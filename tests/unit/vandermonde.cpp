#include "internal/vandermonde.hpp"

#include "internal/linear_algebra.hpp"
#include "internal/performance.hpp"
#include "test_case_dispatch.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace
{

constexpr double singular_tolerance = 1.0e-14;
constexpr const char* singular_message = "Weighted Vandermonde system is singular.";

void assert_near(const double actual, const double expected, const double tolerance)
{
    assert(std::abs(actual - expected) <= tolerance);
}

void assert_close(const double actual, const double expected, const double absolute_tolerance, const double relative_tolerance)
{
    assert(std::abs(actual - expected) <= absolute_tolerance + relative_tolerance * std::abs(expected));
}

qmomkit::internal::DenseMatrix make_solution(const std::vector<std::vector<double>>& rows)
{
    const std::size_t row_count = rows.size();
    const std::size_t col_count = rows.empty() ? 0 : rows.front().size();
    qmomkit::internal::DenseMatrix solution(row_count, col_count);
    for (std::size_t row = 0; row < row_count; ++row)
    {
        assert(rows[row].size() == col_count);
        for (std::size_t col = 0; col < col_count; ++col)
        {
            solution(row, col) = rows[row][col];
        }
    }
    return solution;
}

qmomkit::internal::DenseMatrix make_rhs(const std::vector<double>& abscissas, const std::vector<double>& weights,
                                        const qmomkit::internal::DenseMatrix& solution)
{
    const std::size_t n = abscissas.size();
    qmomkit::internal::DenseMatrix rhs(n, solution.cols());
    for (std::size_t order = 0; order < n; ++order)
    {
        for (std::size_t col = 0; col < solution.cols(); ++col)
        {
            double value = 0.0;
            double power = 1.0;
            for (std::size_t node = 0; node < n; ++node)
            {
                if (order == 0)
                {
                    power = 1.0;
                }
                else
                {
                    power = std::pow(abscissas[node], static_cast<int>(order));
                }
                value += weights[node] * power * solution(node, col);
            }
            rhs(order, col) = value;
        }
    }
    return rhs;
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

void assert_matrix_close(const qmomkit::internal::DenseMatrix& actual, const qmomkit::internal::DenseMatrix& expected,
                         const double absolute_tolerance, const double relative_tolerance)
{
    assert(actual.rows() == expected.rows());
    assert(actual.cols() == expected.cols());
    for (std::size_t row = 0; row < actual.rows(); ++row)
    {
        for (std::size_t col = 0; col < actual.cols(); ++col)
        {
            assert_close(actual(row, col), expected(row, col), absolute_tolerance, relative_tolerance);
        }
    }
}

void assert_matches_dense_reference(const std::vector<double>& abscissas, const std::vector<double>& weights,
                                    const qmomkit::internal::DenseMatrix& expected_solution, const double solution_tolerance,
                                    const double residual_tolerance)
{
    const qmomkit::internal::DenseMatrix rhs = make_rhs(abscissas, weights, expected_solution);
    const qmomkit::internal::DenseMatrix system = qmomkit::internal::make_conditional_vandermonde_matrix(abscissas, weights);

    const qmomkit::internal::LinearSolveResult dense = qmomkit::internal::solve_dense(system, rhs, singular_tolerance, singular_message);
    assert(dense.success);

    const qmomkit::internal::VandermondeSolveOptions options{qmomkit::VandermondeSolver::bjorck_pereyra, singular_tolerance, false, 1.0e-10};
    const qmomkit::internal::VandermondeSolveResult structured =
        qmomkit::internal::solve_conditional_vandermonde_system(abscissas, weights, rhs, options, singular_message);
    assert(structured.success);

    assert_matrix_near(structured.solution, dense.solution, solution_tolerance);
    assert_matrix_near(structured.solution, expected_solution, solution_tolerance);
    assert(qmomkit::internal::residual_norm(system, structured.solution, rhs) <= residual_tolerance);
}

void assert_wrapper_matches_dense_reference(const std::vector<double>& abscissas, const std::vector<double>& weights,
                                            const qmomkit::internal::DenseMatrix& expected_solution, const qmomkit::VandermondeSolver solver,
                                            const qmomkit::VandermondeSolver expected_selected_solver,
                                            const qmomkit::internal::VandermondeAlgorithm expected_algorithm_used, const double solution_tolerance)
{
    const qmomkit::internal::DenseMatrix rhs = make_rhs(abscissas, weights, expected_solution);
    const qmomkit::internal::DenseMatrix system = qmomkit::internal::make_conditional_vandermonde_matrix(abscissas, weights);
    const qmomkit::internal::LinearSolveResult dense = qmomkit::internal::solve_dense(system, rhs, singular_tolerance, singular_message);
    assert(dense.success);

    const qmomkit::internal::VandermondeSolveOptions options{solver, singular_tolerance, true, 1.0e-10};
    const qmomkit::internal::VandermondeSolveResult result =
        qmomkit::internal::solve_conditional_vandermonde_system(abscissas, weights, rhs, options, singular_message);
    assert(result.success);
    assert(result.requested_solver == solver);
    assert(result.selected_solver == expected_selected_solver);
    assert(result.algorithm_used == expected_algorithm_used);
    assert_matrix_near(result.solution, dense.solution, solution_tolerance);
}

void assert_pivoted_qr_success(const std::vector<double>& abscissas, const std::vector<double>& weights,
                               const qmomkit::internal::DenseMatrix& expected_solution, const bool validate_residual, const double residual_tolerance,
                               const std::size_t expected_residual_checks, const double solution_tolerance)
{
    const qmomkit::internal::DenseMatrix rhs = make_rhs(abscissas, weights, expected_solution);
    const qmomkit::internal::DenseMatrix system = qmomkit::internal::make_conditional_vandermonde_matrix(abscissas, weights);
    const qmomkit::internal::LinearSolveResult dense = qmomkit::internal::solve_dense(system, rhs, singular_tolerance, singular_message);
    assert(dense.success);

    qmomkit::internal::performance::CounterSession session;
    qmomkit::internal::VandermondeSolveResult result{qmomkit::internal::DenseMatrix(0, 0),
                                                     false,
                                                     "",
                                                     qmomkit::VandermondeSolver::pivoted_qr,
                                                     qmomkit::VandermondeSolver::pivoted_qr,
                                                     qmomkit::internal::VandermondeAlgorithm::pivoted_qr};
    {
        qmomkit::internal::performance::CounterScope scope(session);
        const qmomkit::internal::VandermondeSolveOptions options{qmomkit::VandermondeSolver::pivoted_qr, singular_tolerance, validate_residual,
                                                                 residual_tolerance};
        result = qmomkit::internal::solve_conditional_vandermonde_system(abscissas, weights, rhs, options, singular_message);
    }

    assert(result.success);
    assert(result.requested_solver == qmomkit::VandermondeSolver::pivoted_qr);
    assert(result.selected_solver == qmomkit::VandermondeSolver::pivoted_qr);
    assert(result.algorithm_used == qmomkit::internal::VandermondeAlgorithm::pivoted_qr);
    assert_matrix_near(result.solution, dense.solution, solution_tolerance);

    const qmomkit::internal::performance::PerformanceCounters counters = session.snapshot();
#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)
    assert(counters.linear_algebra.vandermonde_solves == 1);
    assert(counters.linear_algebra.pivoted_qr_solves == 1);
    assert(counters.linear_algebra.svd_solves == 0);
    assert(counters.linear_algebra.bjorck_pereyra_solves == 0);
    assert(counters.linear_algebra.dense_linear_solves == 0);
    assert(counters.linear_algebra.residual_checks == expected_residual_checks);
#else
    (void)expected_residual_checks;
    assert(counters.linear_algebra.vandermonde_solves == 0);
    assert(counters.linear_algebra.pivoted_qr_solves == 0);
    assert(counters.linear_algebra.svd_solves == 0);
    assert(counters.linear_algebra.bjorck_pereyra_solves == 0);
    assert(counters.linear_algebra.dense_linear_solves == 0);
    assert(counters.linear_algebra.residual_checks == 0);
#endif
}

void assert_svd_success(const std::vector<double>& abscissas, const std::vector<double>& weights,
                        const qmomkit::internal::DenseMatrix& expected_solution, const bool validate_residual, const double residual_tolerance,
                        const std::size_t expected_residual_checks, const double absolute_tolerance, const double relative_tolerance)
{
    const qmomkit::internal::DenseMatrix rhs = make_rhs(abscissas, weights, expected_solution);
    const qmomkit::internal::DenseMatrix system = qmomkit::internal::make_conditional_vandermonde_matrix(abscissas, weights);
    const qmomkit::internal::LinearSolveResult dense = qmomkit::internal::solve_dense(system, rhs, singular_tolerance, singular_message);
    assert(dense.success);

    qmomkit::internal::performance::CounterSession session;
    qmomkit::internal::VandermondeSolveResult result{
        qmomkit::internal::DenseMatrix(0, 0),        false, "", qmomkit::VandermondeSolver::svd, qmomkit::VandermondeSolver::svd,
        qmomkit::internal::VandermondeAlgorithm::svd};
    {
        qmomkit::internal::performance::CounterScope scope(session);
        const qmomkit::internal::VandermondeSolveOptions options{qmomkit::VandermondeSolver::svd, singular_tolerance, validate_residual,
                                                                 residual_tolerance};
        result = qmomkit::internal::solve_conditional_vandermonde_system(abscissas, weights, rhs, options, singular_message);
    }

    assert(result.success);
    assert(result.requested_solver == qmomkit::VandermondeSolver::svd);
    assert(result.selected_solver == qmomkit::VandermondeSolver::svd);
    assert(result.algorithm_used == qmomkit::internal::VandermondeAlgorithm::svd);
    assert_matrix_close(result.solution, dense.solution, absolute_tolerance, relative_tolerance);
    assert_matrix_close(result.solution, expected_solution, absolute_tolerance, relative_tolerance);

    const qmomkit::internal::performance::PerformanceCounters counters = session.snapshot();
#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)
    assert(counters.linear_algebra.vandermonde_solves == 1);
    assert(counters.linear_algebra.svd_solves == 1);
    assert(counters.linear_algebra.pivoted_qr_solves == 0);
    assert(counters.linear_algebra.bjorck_pereyra_solves == 0);
    assert(counters.linear_algebra.dense_linear_solves == 0);
    assert(counters.linear_algebra.residual_checks == expected_residual_checks);
#else
    (void)expected_residual_checks;
    assert(counters.linear_algebra.vandermonde_solves == 0);
    assert(counters.linear_algebra.svd_solves == 0);
    assert(counters.linear_algebra.pivoted_qr_solves == 0);
    assert(counters.linear_algebra.bjorck_pereyra_solves == 0);
    assert(counters.linear_algebra.dense_linear_solves == 0);
    assert(counters.linear_algebra.residual_checks == 0);
#endif
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

std::vector<double> make_centered_nodes(const std::size_t n)
{
    std::vector<double> nodes(n, 0.0);
    if (n == 1)
    {
        nodes.front() = 0.25;
        return nodes;
    }
    for (std::size_t i = 0; i < n; ++i)
    {
        nodes[i] = -0.8 + 1.6 * static_cast<double>(i) / static_cast<double>(n - 1);
    }
    return nodes;
}

std::vector<double> make_mild_weights(const std::size_t n)
{
    std::vector<double> weights(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        weights[i] = 0.75 + 0.125 * static_cast<double>((i % 5) + 1);
    }
    return weights;
}

void test_one_node_success()
{
    assert_matches_dense_reference({2.5}, {-4.0}, make_solution({{3.0, -7.0, 0.25}}), 1.0e-14, 1.0e-14);
}

void test_well_conditioned_uniform_weights()
{
    assert_matches_dense_reference({0.0, 1.0, 3.0}, {1.0, 1.0, 1.0}, make_solution({{1.0, 2.0}, {-3.0, 4.5}, {5.0, -6.0}}), 1.0e-12, 1.0e-12);
}

void test_mixed_sign_abscissas()
{
    assert_matches_dense_reference({-2.0, 0.25, 1.5, 4.0}, {0.5, 1.25, -0.75, 2.0},
                                   make_solution({{1.0, -2.0, 3.0}, {0.5, 4.0, -1.0}, {-3.0, 2.5, 0.0}, {2.0, -1.0, 6.0}}), 1.0e-11, 1.0e-10);
}

void test_clustered_nodes_match_dense_reference()
{
    assert_matches_dense_reference({1.0, 1.0 + 1.0e-5, 1.0 + 3.0e-5}, {0.75, -1.5, 2.25}, make_solution({{1.0, 0.25}, {-0.5, 2.0}, {3.0, -4.0}}),
                                   2.0e-5, 1.0e-11);
}

void test_nearly_repeated_nodes_match_dense_reference()
{
    assert_matches_dense_reference({-1.0, -1.0 + 1.0e-10}, {0.5, 0.75}, make_solution({{2.0, -1.0}, {4.0, 3.0}}), 1.0e-5, 1.0e-14);
}

void test_larger_node_count_smoke_cases()
{
    assert_matches_dense_reference(make_centered_nodes(8), make_mild_weights(8), make_patterned_solution(8, 4), 1.0e-9, 1.0e-12);
    assert_matches_dense_reference(make_centered_nodes(16), make_mild_weights(16), make_patterned_solution(16, 3), 1.0e-5, 1.0e-10);
}

void test_zero_rhs_columns_success()
{
    const std::vector<double> abscissas{-1.0, 0.5, 2.0};
    const std::vector<double> weights{0.25, 0.5, 0.75};
    const qmomkit::internal::DenseMatrix rhs(abscissas.size(), 0);
    const qmomkit::internal::VandermondeSolveOptions options{qmomkit::VandermondeSolver::bjorck_pereyra, singular_tolerance, false, 1.0e-10};
    const qmomkit::internal::VandermondeSolveResult structured =
        qmomkit::internal::solve_conditional_vandermonde_system(abscissas, weights, rhs, options, singular_message);
    assert(structured.success);
    assert(structured.solution.rows() == abscissas.size());
    assert(structured.solution.cols() == 0);
}

void test_large_dynamic_range_rhs()
{
    // The weighted system has an infinity-norm condition number of about
    // 2e8. Keep its widely scaled right-hand sides, but use small manufactured
    // coefficients that remain resolvable at that conditioning level.
    const std::vector<double> abscissas{-0.75, -0.1, 0.4, 0.9};
    const std::vector<double> weights{1.0e-3, -2.0, 5.0e2, -1.0e5};
    const qmomkit::internal::DenseMatrix expected =
        make_solution({{1.0e-4, -1.0e12, 3.0}, {-2.0e-4, 4.0e9, -5.0}, {3.0e-4, -6.0e6, 7.0}, {-4.0e-4, 8.0e3, -9.0}});
    const qmomkit::internal::DenseMatrix rhs = make_rhs(abscissas, weights, expected);
    const qmomkit::internal::DenseMatrix system = qmomkit::internal::make_conditional_vandermonde_matrix(abscissas, weights);
    const qmomkit::internal::VandermondeSolveOptions options{qmomkit::VandermondeSolver::bjorck_pereyra, singular_tolerance, false, 1.0e-10};
    const qmomkit::internal::VandermondeSolveResult result =
        qmomkit::internal::solve_conditional_vandermonde_system(abscissas, weights, rhs, options, singular_message);

    assert(result.success);
    assert(result.solution.rows() == expected.rows());
    assert(result.solution.cols() == expected.cols());

    constexpr double relative_solution_tolerance = 1.0e-6;
    constexpr double near_zero_absolute_tolerance = 1.0e-12;
    double expected_solution_infinity_norm = 0.0;
    for (std::size_t row = 0; row < expected.rows(); ++row)
    {
        for (std::size_t column = 0; column < expected.cols(); ++column)
        {
            const double expected_value = expected(row, column);
            const double solution_tolerance = near_zero_absolute_tolerance + relative_solution_tolerance * std::abs(expected_value);
            assert(std::abs(result.solution(row, column) - expected_value) <= solution_tolerance);
            expected_solution_infinity_norm = std::max(expected_solution_infinity_norm, std::abs(expected_value));
        }
    }

    double system_infinity_norm = 0.0;
    double rhs_infinity_norm = 0.0;
    for (std::size_t row = 0; row < system.rows(); ++row)
    {
        double row_sum = 0.0;
        for (std::size_t column = 0; column < system.cols(); ++column)
        {
            row_sum += std::abs(system(row, column));
        }
        system_infinity_norm = std::max(system_infinity_norm, row_sum);
    }
    for (std::size_t row = 0; row < rhs.rows(); ++row)
    {
        for (std::size_t column = 0; column < rhs.cols(); ++column)
        {
            rhs_infinity_norm = std::max(rhs_infinity_norm, std::abs(rhs(row, column)));
        }
    }

    const double residual_scale = system_infinity_norm * expected_solution_infinity_norm + rhs_infinity_norm;
    assert(qmomkit::internal::residual_norm(system, result.solution, rhs) <= 1.0e-10 * residual_scale);
}

void test_wrapper_forced_modes_match_dense_reference()
{
    const std::vector<double> abscissas{-0.75, -0.1, 0.4, 0.9};
    const std::vector<double> weights{0.25, 0.5, -0.75, 1.25};
    const qmomkit::internal::DenseMatrix expected_solution = make_patterned_solution(4, 8);

    assert_wrapper_matches_dense_reference(abscissas, weights, expected_solution, qmomkit::VandermondeSolver::dense,
                                           qmomkit::VandermondeSolver::dense, qmomkit::internal::VandermondeAlgorithm::dense, 1.0e-12);
    assert_wrapper_matches_dense_reference(abscissas, weights, expected_solution, qmomkit::VandermondeSolver::bjorck_pereyra,
                                           qmomkit::VandermondeSolver::bjorck_pereyra, qmomkit::internal::VandermondeAlgorithm::bjorck_pereyra,
                                           1.0e-11);
    assert_wrapper_matches_dense_reference(abscissas, weights, expected_solution, qmomkit::VandermondeSolver::pivoted_qr,
                                           qmomkit::VandermondeSolver::pivoted_qr, qmomkit::internal::VandermondeAlgorithm::pivoted_qr, 1.0e-12);
    assert_wrapper_matches_dense_reference(abscissas, weights, expected_solution, qmomkit::VandermondeSolver::svd, qmomkit::VandermondeSolver::svd,
                                           qmomkit::internal::VandermondeAlgorithm::svd, 1.0e-11);
    assert_wrapper_matches_dense_reference(abscissas, weights, expected_solution, qmomkit::VandermondeSolver::moore_penrose,
                                           qmomkit::VandermondeSolver::moore_penrose, qmomkit::internal::VandermondeAlgorithm::moore_penrose,
                                           1.0e-11);
}

void test_wrapper_pivoted_qr_solves_weighted_system_without_residual_by_default()
{
    assert_pivoted_qr_success({-0.75, -0.1, 0.4, 0.9}, {0.25, 0.5, -0.75, 1.25}, make_patterned_solution(4, 8), false, -1.0, 0, 1.0e-12);
}

void test_wrapper_pivoted_qr_honors_optional_residual_validation()
{
    assert_pivoted_qr_success({1.0, 1.0 + 1.0e-5, 1.0 + 3.0e-5}, {0.75, -1.5, 2.25}, make_solution({{1.0, 0.25}, {-0.5, 2.0}, {3.0, -4.0}}), true,
                              1.0e-10, 1, 2.0e-5);
}

void test_wrapper_pivoted_qr_zero_rhs_columns_success()
{
    const std::vector<double> abscissas{-1.0, 0.5, 2.0};
    const std::vector<double> weights{0.25, 0.5, 0.75};
    const qmomkit::internal::DenseMatrix rhs(abscissas.size(), 0);
    const qmomkit::internal::VandermondeSolveOptions options{qmomkit::VandermondeSolver::pivoted_qr, singular_tolerance, false, 1.0e-10};
    const qmomkit::internal::VandermondeSolveResult result =
        qmomkit::internal::solve_conditional_vandermonde_system(abscissas, weights, rhs, options, singular_message);
    assert(result.success);
    assert(result.selected_solver == qmomkit::VandermondeSolver::pivoted_qr);
    assert(result.algorithm_used == qmomkit::internal::VandermondeAlgorithm::pivoted_qr);
    assert(result.solution.rows() == abscissas.size());
    assert(result.solution.cols() == 0);
}

void test_wrapper_svd_solves_weighted_system_without_residual_by_default()
{
    assert_svd_success({-0.75, -0.1, 0.4, 0.9}, {0.25, 0.5, -0.75, 1.25}, make_patterned_solution(4, 8), false, -1.0, 0, 1.0e-11, 1.0e-11);
}

void test_wrapper_svd_honors_optional_residual_validation()
{
    // The clustered abscissas below give a weighted Vandermonde system with a
    // condition number near 6.7e10, so the solution can only be pinned down to
    // roughly eps * cond = 1.5e-5. Compare against that floor with a safety
    // margin: a tighter bound would make the test depend on which rounding a
    // given toolchain happens to produce.
    constexpr double ill_conditioned_solution_tolerance = 1.0e-3;
    assert_svd_success({1.0, 1.0 + 1.0e-5, 1.0 + 3.0e-5}, {0.75, -1.5, 2.25}, make_solution({{1.0, 0.25}, {-0.5, 2.0}, {3.0, -4.0}}), true, 1.0e-10,
                       1, ill_conditioned_solution_tolerance, 1.0e-10);
}

void test_wrapper_svd_zero_rhs_columns_success()
{
    const std::vector<double> abscissas{-1.0, 0.5, 2.0};
    const std::vector<double> weights{0.25, 0.5, 0.75};
    const qmomkit::internal::DenseMatrix rhs(abscissas.size(), 0);
    const qmomkit::internal::VandermondeSolveOptions options{qmomkit::VandermondeSolver::svd, singular_tolerance, false, 1.0e-10};
    const qmomkit::internal::VandermondeSolveResult result =
        qmomkit::internal::solve_conditional_vandermonde_system(abscissas, weights, rhs, options, singular_message);
    assert(result.success);
    assert(result.requested_solver == qmomkit::VandermondeSolver::svd);
    assert(result.selected_solver == qmomkit::VandermondeSolver::svd);
    assert(result.algorithm_used == qmomkit::internal::VandermondeAlgorithm::svd);
    assert(result.solution.rows() == abscissas.size());
    assert(result.solution.cols() == 0);
}

void test_wrapper_moore_penrose_solves_rank_deficient_system_with_minimum_norm_solution()
{
    const std::vector<double> abscissas{1.0, 1.0};
    const std::vector<double> weights{0.5, 0.5};
    const qmomkit::internal::DenseMatrix rhs = make_solution({{1.0}, {1.0}});
    qmomkit::internal::performance::CounterSession session;
    qmomkit::internal::VandermondeSolveResult result{qmomkit::internal::DenseMatrix(0, 0),
                                                     false,
                                                     "",
                                                     qmomkit::VandermondeSolver::moore_penrose,
                                                     qmomkit::VandermondeSolver::moore_penrose,
                                                     qmomkit::internal::VandermondeAlgorithm::moore_penrose};
    {
        qmomkit::internal::performance::CounterScope scope(session);
        const qmomkit::internal::VandermondeSolveOptions options{qmomkit::VandermondeSolver::moore_penrose, singular_tolerance, true, 1.0e-10,
                                                                 1.0e-14};
        result = qmomkit::internal::solve_conditional_vandermonde_system(abscissas, weights, rhs, options, singular_message);
    }

    assert(result.success);
    assert(result.requested_solver == qmomkit::VandermondeSolver::moore_penrose);
    assert(result.selected_solver == qmomkit::VandermondeSolver::moore_penrose);
    assert(result.algorithm_used == qmomkit::internal::VandermondeAlgorithm::moore_penrose);
    assert_matrix_near(result.solution, make_solution({{1.0}, {1.0}}), 1.0e-12);

    const qmomkit::internal::performance::PerformanceCounters counters = session.snapshot();
#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)
    assert(counters.linear_algebra.vandermonde_solves == 1);
    assert(counters.linear_algebra.svd_solves == 1);
    assert(counters.linear_algebra.pivoted_qr_solves == 0);
    assert(counters.linear_algebra.bjorck_pereyra_solves == 0);
    assert(counters.linear_algebra.dense_linear_solves == 0);
    assert(counters.linear_algebra.residual_checks == 1);
#else
    assert(counters.linear_algebra.vandermonde_solves == 0);
    assert(counters.linear_algebra.svd_solves == 0);
    assert(counters.linear_algebra.pivoted_qr_solves == 0);
    assert(counters.linear_algebra.bjorck_pereyra_solves == 0);
    assert(counters.linear_algebra.dense_linear_solves == 0);
    assert(counters.linear_algebra.residual_checks == 0);
#endif
}

void test_wrapper_moore_penrose_rcond_truncates_full_rank_system()
{
    const std::vector<double> abscissas{0.0, 1.0};
    const std::vector<double> weights{1.0, 1.0};
    const qmomkit::internal::DenseMatrix rhs = make_solution({{1.0}, {0.0}});
    const qmomkit::internal::VandermondeSolveOptions low_rcond{qmomkit::VandermondeSolver::moore_penrose, singular_tolerance, false, 1.0e-10, 0.1};
    const qmomkit::internal::VandermondeSolveOptions high_rcond{qmomkit::VandermondeSolver::moore_penrose, singular_tolerance, false, 1.0e-10, 0.5};
    const qmomkit::internal::VandermondeSolveResult low =
        qmomkit::internal::solve_conditional_vandermonde_system(abscissas, weights, rhs, low_rcond, singular_message);
    const qmomkit::internal::VandermondeSolveResult high =
        qmomkit::internal::solve_conditional_vandermonde_system(abscissas, weights, rhs, high_rcond, singular_message);

    assert(low.success);
    assert(high.success);
    assert_matrix_near(low.solution, make_solution({{1.0}, {0.0}}), 1.0e-12);
    assert_matrix_near(high.solution, make_solution({{0.276393202250021}, {0.447213595499958}}), 1.0e-12);
}

void test_wrapper_automatic_selects_plain_bp_for_small_wide_system()
{
    const std::vector<double> abscissas = make_centered_nodes(8);
    const std::vector<double> weights = make_mild_weights(8);
    assert_wrapper_matches_dense_reference(abscissas, weights, make_patterned_solution(8, 256), qmomkit::VandermondeSolver::automatic,
                                           qmomkit::VandermondeSolver::bjorck_pereyra, qmomkit::internal::VandermondeAlgorithm::bjorck_pereyra,
                                           1.0e-9);
}

void test_wrapper_automatic_selects_plain_bp_for_large_system()
{
    const std::vector<double> abscissas = make_centered_nodes(16);
    const std::vector<double> weights = make_mild_weights(16);
    assert_wrapper_matches_dense_reference(abscissas, weights, make_patterned_solution(16, 8), qmomkit::VandermondeSolver::automatic,
                                           qmomkit::VandermondeSolver::bjorck_pereyra, qmomkit::internal::VandermondeAlgorithm::bjorck_pereyra,
                                           1.0e-5);
}

void test_wrapper_automatic_matches_plain_bp_without_forced_residual_validation()
{
    const std::vector<double> abscissas{-0.75, -0.1, 0.4, 0.9};
    const std::vector<double> weights{0.25, 0.5, -0.75, 1.25};
    const qmomkit::internal::DenseMatrix expected_solution = make_patterned_solution(4, 8);
    const qmomkit::internal::DenseMatrix rhs = make_rhs(abscissas, weights, expected_solution);
    const qmomkit::internal::VandermondeSolveOptions options{qmomkit::VandermondeSolver::automatic, singular_tolerance, false, 1.0e-10};

    qmomkit::internal::performance::CounterSession session;
    qmomkit::internal::VandermondeSolveResult result{
        qmomkit::internal::DenseMatrix(0, 0),          false, "", qmomkit::VandermondeSolver::automatic, qmomkit::VandermondeSolver::automatic,
        qmomkit::internal::VandermondeAlgorithm::dense};
    {
        qmomkit::internal::performance::CounterScope scope(session);
        result = qmomkit::internal::solve_conditional_vandermonde_system(abscissas, weights, rhs, options, singular_message);
    }

    assert(result.success);
    assert(result.requested_solver == qmomkit::VandermondeSolver::automatic);
    assert(result.selected_solver == qmomkit::VandermondeSolver::bjorck_pereyra);
    assert(result.algorithm_used == qmomkit::internal::VandermondeAlgorithm::bjorck_pereyra);
    assert_matrix_near(result.solution, expected_solution, 1.0e-11);

    const qmomkit::internal::performance::PerformanceCounters counters = session.snapshot();
#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)
    assert(counters.linear_algebra.vandermonde_solves == 1);
    assert(counters.linear_algebra.bjorck_pereyra_solves == 1);
    assert(counters.linear_algebra.dense_linear_solves == 0);
    assert(counters.linear_algebra.residual_checks == 0);
#else
    assert(counters.linear_algebra.vandermonde_solves == 0);
    assert(counters.linear_algebra.bjorck_pereyra_solves == 0);
    assert(counters.linear_algebra.dense_linear_solves == 0);
    assert(counters.linear_algebra.residual_checks == 0);
#endif
}

void assert_singular_failure(const std::vector<double>& abscissas, const std::vector<double>& weights, const qmomkit::internal::DenseMatrix& rhs,
                             const std::string& expected_message)
{
    const qmomkit::internal::VandermondeSolveOptions options{qmomkit::VandermondeSolver::bjorck_pereyra, singular_tolerance, false, 1.0e-10};
    const qmomkit::internal::VandermondeSolveResult structured =
        qmomkit::internal::solve_conditional_vandermonde_system(abscissas, weights, rhs, options, singular_message);
    assert(!structured.success);
    assert(structured.message == expected_message);
}

struct InvalidPolicyExpectation
{
    qmomkit::VandermondeSolver requested_solver;
    qmomkit::VandermondeSolver selected_solver;
    qmomkit::internal::VandermondeAlgorithm algorithm_used;
    bool success;
    const char* message;
    std::size_t vandermonde_solves;
    std::size_t bjorck_pereyra_solves;
    std::size_t pivoted_qr_solves;
    std::size_t svd_solves;
    std::size_t dense_linear_solves;
    std::size_t residual_checks;
};

void assert_counter_expectations(const qmomkit::internal::performance::PerformanceCounters& counters, const InvalidPolicyExpectation& expected)
{
#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)
    assert(counters.linear_algebra.vandermonde_solves == expected.vandermonde_solves);
    assert(counters.linear_algebra.bjorck_pereyra_solves == expected.bjorck_pereyra_solves);
    assert(counters.linear_algebra.pivoted_qr_solves == expected.pivoted_qr_solves);
    assert(counters.linear_algebra.svd_solves == expected.svd_solves);
    assert(counters.linear_algebra.dense_linear_solves == expected.dense_linear_solves);
    assert(counters.linear_algebra.residual_checks == expected.residual_checks);
#else
    (void)expected;
    assert(counters.linear_algebra.vandermonde_solves == 0);
    assert(counters.linear_algebra.bjorck_pereyra_solves == 0);
    assert(counters.linear_algebra.pivoted_qr_solves == 0);
    assert(counters.linear_algebra.svd_solves == 0);
    assert(counters.linear_algebra.dense_linear_solves == 0);
    assert(counters.linear_algebra.residual_checks == 0);
#endif
}

void assert_invalid_policy_behavior(const InvalidPolicyExpectation& expected)
{
    qmomkit::internal::performance::CounterSession session;
    qmomkit::internal::VandermondeSolveResult result{
        qmomkit::internal::DenseMatrix(0, 0), false, "", expected.requested_solver, expected.selected_solver, expected.algorithm_used};
    {
        qmomkit::internal::performance::CounterScope scope(session);
        const std::vector<double> abscissas{1.0, std::numeric_limits<double>::infinity()};
        const std::vector<double> weights{1.0, 1.0};
        const qmomkit::internal::DenseMatrix rhs(2, 1);
        const qmomkit::internal::VandermondeSolveOptions options{expected.requested_solver, singular_tolerance, false, 1.0e-10};
        result = qmomkit::internal::solve_conditional_vandermonde_system(abscissas, weights, rhs, options, singular_message);
    }

    assert(result.requested_solver == expected.requested_solver);
    assert(result.selected_solver == expected.selected_solver);
    assert(result.algorithm_used == expected.algorithm_used);
    assert(result.success == expected.success);
    assert(result.message == expected.message);
    assert_counter_expectations(session.snapshot(), expected);
}

void test_repeated_nodes_fail()
{
    assert_singular_failure({1.0, 1.0}, {0.5, 0.5}, make_solution({{1.0, 2.0}, {3.0, 4.0}}), singular_message);
}

void test_nonuniform_weights()
{
    assert_matches_dense_reference({-3.0, -0.5, 2.0}, {1.0e-3, -2.0, 7.5}, make_solution({{2.0}, {-4.0}, {8.0}}), 1.0e-10, 1.0e-12);
}

void test_small_weight_failure()
{
    assert_singular_failure({-1.0, 2.0}, {1.0, 1.0e-16}, make_solution({{1.0}, {2.0}}), singular_message);
}

void test_zero_weight_failure()
{
    assert_singular_failure({-1.0, 2.0}, {1.0, 0.0}, make_solution({{1.0}, {2.0}}), singular_message);
}

void assert_pivoted_qr_failure(const std::vector<double>& abscissas, const std::vector<double>& weights, const qmomkit::internal::DenseMatrix& rhs,
                               const bool validate_residual, const double residual_tolerance)
{
    qmomkit::internal::performance::CounterSession session;
    qmomkit::internal::VandermondeSolveResult result{qmomkit::internal::DenseMatrix(0, 0),
                                                     true,
                                                     "",
                                                     qmomkit::VandermondeSolver::pivoted_qr,
                                                     qmomkit::VandermondeSolver::pivoted_qr,
                                                     qmomkit::internal::VandermondeAlgorithm::pivoted_qr};
    {
        qmomkit::internal::performance::CounterScope scope(session);
        const qmomkit::internal::VandermondeSolveOptions options{qmomkit::VandermondeSolver::pivoted_qr, singular_tolerance, validate_residual,
                                                                 residual_tolerance};
        result = qmomkit::internal::solve_conditional_vandermonde_system(abscissas, weights, rhs, options, singular_message);
    }

    assert(!result.success);
    assert(result.message == singular_message);
    assert(result.requested_solver == qmomkit::VandermondeSolver::pivoted_qr);
    assert(result.selected_solver == qmomkit::VandermondeSolver::pivoted_qr);
    assert(result.algorithm_used == qmomkit::internal::VandermondeAlgorithm::pivoted_qr);

    const qmomkit::internal::performance::PerformanceCounters counters = session.snapshot();
#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)
    assert(counters.linear_algebra.vandermonde_solves == 1);
    assert(counters.linear_algebra.pivoted_qr_solves == 1);
    assert(counters.linear_algebra.svd_solves == 0);
    assert(counters.linear_algebra.bjorck_pereyra_solves == 0);
    assert(counters.linear_algebra.dense_linear_solves == 0);
    assert(counters.linear_algebra.residual_checks == (validate_residual ? 1 : 0));
#else
    assert(counters.linear_algebra.vandermonde_solves == 0);
    assert(counters.linear_algebra.pivoted_qr_solves == 0);
    assert(counters.linear_algebra.svd_solves == 0);
    assert(counters.linear_algebra.bjorck_pereyra_solves == 0);
    assert(counters.linear_algebra.dense_linear_solves == 0);
    assert(counters.linear_algebra.residual_checks == 0);
#endif
}

void assert_svd_failure(const std::vector<double>& abscissas, const std::vector<double>& weights, const qmomkit::internal::DenseMatrix& rhs,
                        const bool validate_residual, const double residual_tolerance, const std::string& expected_message,
                        const std::size_t expected_residual_checks)
{
    qmomkit::internal::performance::CounterSession session;
    qmomkit::internal::VandermondeSolveResult result{
        qmomkit::internal::DenseMatrix(0, 0),        true, "", qmomkit::VandermondeSolver::svd, qmomkit::VandermondeSolver::svd,
        qmomkit::internal::VandermondeAlgorithm::svd};
    {
        qmomkit::internal::performance::CounterScope scope(session);
        const qmomkit::internal::VandermondeSolveOptions options{qmomkit::VandermondeSolver::svd, singular_tolerance, validate_residual,
                                                                 residual_tolerance};
        result = qmomkit::internal::solve_conditional_vandermonde_system(abscissas, weights, rhs, options, singular_message);
    }

    assert(!result.success);
    assert(result.message == expected_message);
    assert(result.requested_solver == qmomkit::VandermondeSolver::svd);
    assert(result.selected_solver == qmomkit::VandermondeSolver::svd);
    assert(result.algorithm_used == qmomkit::internal::VandermondeAlgorithm::svd);

    const qmomkit::internal::performance::PerformanceCounters counters = session.snapshot();
#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)
    assert(counters.linear_algebra.vandermonde_solves == 1);
    assert(counters.linear_algebra.svd_solves == 1);
    assert(counters.linear_algebra.pivoted_qr_solves == 0);
    assert(counters.linear_algebra.bjorck_pereyra_solves == 0);
    assert(counters.linear_algebra.dense_linear_solves == 0);
    assert(counters.linear_algebra.residual_checks == expected_residual_checks);
#else
    (void)expected_residual_checks;
    assert(counters.linear_algebra.vandermonde_solves == 0);
    assert(counters.linear_algebra.svd_solves == 0);
    assert(counters.linear_algebra.pivoted_qr_solves == 0);
    assert(counters.linear_algebra.bjorck_pereyra_solves == 0);
    assert(counters.linear_algebra.dense_linear_solves == 0);
    assert(counters.linear_algebra.residual_checks == 0);
#endif
}

void test_wrapper_pivoted_qr_reports_rank_deficiency()
{
    assert_pivoted_qr_failure({1.0, 1.0}, {0.5, 0.5}, make_solution({{1.0, 2.0}, {3.0, 4.0}}), false, 1.0e-10);
    assert_pivoted_qr_failure({-1.0, 2.0}, {1.0, 1.0e-16}, make_solution({{1.0}, {2.0}}), false, 1.0e-10);
    assert_pivoted_qr_failure({-1.0, 2.0}, {1.0, 0.0}, make_solution({{1.0}, {2.0}}), false, 1.0e-10);
}

void test_wrapper_svd_reports_rank_deficiency()
{
    assert_svd_failure({1.0, 1.0}, {0.5, 0.5}, make_solution({{1.0, 2.0}, {3.0, 4.0}}), false, 1.0e-10, singular_message, 0);
    assert_svd_failure({-1.0, 2.0}, {1.0, 1.0e-16}, make_solution({{1.0}, {2.0}}), false, 1.0e-10, singular_message, 0);
    assert_svd_failure({-1.0, 2.0}, {1.0, 0.0}, make_solution({{1.0}, {2.0}}), false, 1.0e-10, singular_message, 0);
}

void test_wrapper_pivoted_qr_residual_rejection_fails()
{
    const std::vector<double> abscissas{-0.75, -0.1, 0.4, 0.9};
    const std::vector<double> weights{0.25, 0.5, -0.75, 1.25};
    const qmomkit::internal::DenseMatrix rhs = make_rhs(abscissas, weights, make_patterned_solution(4, 8));
    assert_pivoted_qr_failure(abscissas, weights, rhs, true, -1.0);
}

void test_wrapper_svd_residual_rejection_fails()
{
    const std::vector<double> abscissas{-0.75, -0.1, 0.4, 0.9};
    const std::vector<double> weights{0.25, 0.5, -0.75, 1.25};
    const qmomkit::internal::DenseMatrix rhs = make_rhs(abscissas, weights, make_patterned_solution(4, 8));
    assert_svd_failure(abscissas, weights, rhs, true, -1.0, singular_message, 1);
}

void test_invalid_inputs_fail()
{
    assert_singular_failure({}, {}, qmomkit::internal::DenseMatrix(0, 1), "Weighted Vandermonde system must not be empty.");
    assert_singular_failure({1.0, 2.0}, {1.0}, qmomkit::internal::DenseMatrix(2, 1), "Weighted Vandermonde dimensions are inconsistent.");
    assert_singular_failure({1.0, std::numeric_limits<double>::infinity()}, {1.0, 1.0}, qmomkit::internal::DenseMatrix(2, 1),
                            "Weighted Vandermonde inputs must be finite.");
}

void test_invalid_inputs_preserve_policy_diagnostics_and_counters()
{
    assert_invalid_policy_behavior({qmomkit::VandermondeSolver::dense, qmomkit::VandermondeSolver::dense,
                                    qmomkit::internal::VandermondeAlgorithm::dense, false, singular_message, 1, 0, 0, 0, 1, 0});
    assert_invalid_policy_behavior({qmomkit::VandermondeSolver::pivoted_qr, qmomkit::VandermondeSolver::pivoted_qr,
                                    qmomkit::internal::VandermondeAlgorithm::pivoted_qr, false, "QR factorization matrix entries must be finite.", 1,
                                    0, 1, 0, 0, 0});
    assert_invalid_policy_behavior({qmomkit::VandermondeSolver::svd, qmomkit::VandermondeSolver::svd, qmomkit::internal::VandermondeAlgorithm::svd,
                                    false, "SVD factorization matrix entries must be finite.", 1, 0, 0, 1, 0, 0});
    assert_invalid_policy_behavior({qmomkit::VandermondeSolver::bjorck_pereyra, qmomkit::VandermondeSolver::bjorck_pereyra,
                                    qmomkit::internal::VandermondeAlgorithm::bjorck_pereyra, false, "Weighted Vandermonde inputs must be finite.", 1,
                                    1, 0, 0, 0, 0});
    assert_invalid_policy_behavior({qmomkit::VandermondeSolver::automatic, qmomkit::VandermondeSolver::bjorck_pereyra,
                                    qmomkit::internal::VandermondeAlgorithm::bjorck_pereyra, false, "Weighted Vandermonde inputs must be finite.", 1,
                                    1, 0, 0, 0, 0});
}

void test_dense_success_requires_finite_solution()
{
    // This scalar weighted system passes dense elimination's absolute pivot
    // threshold, but its exact solution 5e313 is not representable as double.
    const qmomkit::internal::VandermondeSolveOptions options{qmomkit::VandermondeSolver::dense, singular_tolerance, false, 1.0e-10};
    const qmomkit::internal::VandermondeSolveResult result =
        qmomkit::internal::solve_conditional_vandermonde_system({0.0}, {2.0e-14}, make_solution({{1.0e300}}), options, singular_message);

    assert(!result.success || std::isfinite(result.solution(0, 0)));
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"one_node_success", test_one_node_success},
            {"well_conditioned_uniform_weights", test_well_conditioned_uniform_weights},
            {"mixed_sign_abscissas", test_mixed_sign_abscissas},
            {"clustered_nodes_match_dense_reference", test_clustered_nodes_match_dense_reference},
            {"nearly_repeated_nodes_match_dense_reference", test_nearly_repeated_nodes_match_dense_reference},
            {"larger_node_count_smoke_cases", test_larger_node_count_smoke_cases},
            {"zero_rhs_columns_success", test_zero_rhs_columns_success},
            {"large_dynamic_range_rhs", test_large_dynamic_range_rhs},
            {"wrapper_forced_modes_match_dense_reference", test_wrapper_forced_modes_match_dense_reference},
            {"wrapper_pivoted_qr_solves_weighted_system_without_residual_by_default",
             test_wrapper_pivoted_qr_solves_weighted_system_without_residual_by_default},
            {"wrapper_pivoted_qr_honors_optional_residual_validation", test_wrapper_pivoted_qr_honors_optional_residual_validation},
            {"wrapper_pivoted_qr_zero_rhs_columns_success", test_wrapper_pivoted_qr_zero_rhs_columns_success},
            {"wrapper_svd_solves_weighted_system_without_residual_by_default", test_wrapper_svd_solves_weighted_system_without_residual_by_default},
            {"wrapper_svd_honors_optional_residual_validation", test_wrapper_svd_honors_optional_residual_validation},
            {"wrapper_svd_zero_rhs_columns_success", test_wrapper_svd_zero_rhs_columns_success},
            {"wrapper_moore_penrose_solves_rank_deficient_system_with_minimum_norm_solution",
             test_wrapper_moore_penrose_solves_rank_deficient_system_with_minimum_norm_solution},
            {"wrapper_moore_penrose_rcond_truncates_full_rank_system", test_wrapper_moore_penrose_rcond_truncates_full_rank_system},
            {"wrapper_automatic_selects_plain_bp_for_small_wide_system", test_wrapper_automatic_selects_plain_bp_for_small_wide_system},
            {"wrapper_automatic_selects_plain_bp_for_large_system", test_wrapper_automatic_selects_plain_bp_for_large_system},
            {"wrapper_automatic_matches_plain_bp_without_forced_residual_validation",
             test_wrapper_automatic_matches_plain_bp_without_forced_residual_validation},
            {"repeated_nodes_fail", test_repeated_nodes_fail},
            {"nonuniform_weights", test_nonuniform_weights},
            {"small_weight_failure", test_small_weight_failure},
            {"zero_weight_failure", test_zero_weight_failure},
            {"wrapper_pivoted_qr_reports_rank_deficiency", test_wrapper_pivoted_qr_reports_rank_deficiency},
            {"wrapper_svd_reports_rank_deficiency", test_wrapper_svd_reports_rank_deficiency},
            {"wrapper_pivoted_qr_residual_rejection_fails", test_wrapper_pivoted_qr_residual_rejection_fails},
            {"wrapper_svd_residual_rejection_fails", test_wrapper_svd_residual_rejection_fails},
            {"invalid_inputs_fail", test_invalid_inputs_fail},
            {"invalid_inputs_preserve_policy_diagnostics_and_counters", test_invalid_inputs_preserve_policy_diagnostics_and_counters},
            {"dense_success_requires_finite_solution", test_dense_success_requires_finite_solution},
        });
}
