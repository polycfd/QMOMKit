#include "internal/vandermonde.hpp"

#include "internal/linear_algebra.hpp"

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
    const qmomkit::internal::DenseMatrix system = qmomkit::internal::make_weighted_vandermonde(abscissas, weights);

    const qmomkit::internal::LinearSolveResult dense = qmomkit::internal::solve_dense(system, rhs, singular_tolerance, singular_message);
    assert(dense.success);

    const qmomkit::internal::LinearSolveResult structured =
        qmomkit::internal::solve_weighted_transposed_vandermonde_bjorck_pereyra(abscissas, weights, rhs, singular_tolerance, singular_message);
    assert(structured.success);

    assert_matrix_near(structured.solution, dense.solution, solution_tolerance);
    assert_matrix_near(structured.solution, expected_solution, solution_tolerance);
    assert(qmomkit::internal::residual_norm(system, structured.solution, rhs) <= residual_tolerance);
}

void assert_matches_dense_reference_relative(const std::vector<double>& abscissas, const std::vector<double>& weights,
                                             const qmomkit::internal::DenseMatrix& expected_solution, const double absolute_tolerance,
                                             const double relative_tolerance, const double residual_tolerance)
{
    const qmomkit::internal::DenseMatrix rhs = make_rhs(abscissas, weights, expected_solution);
    const qmomkit::internal::DenseMatrix system = qmomkit::internal::make_weighted_vandermonde(abscissas, weights);

    const qmomkit::internal::LinearSolveResult dense = qmomkit::internal::solve_dense(system, rhs, singular_tolerance, singular_message);
    assert(dense.success);

    const qmomkit::internal::LinearSolveResult structured =
        qmomkit::internal::solve_weighted_transposed_vandermonde_bjorck_pereyra(abscissas, weights, rhs, singular_tolerance, singular_message);
    assert(structured.success);

    assert_matrix_close(structured.solution, dense.solution, absolute_tolerance, relative_tolerance);
    assert_matrix_close(structured.solution, expected_solution, absolute_tolerance, relative_tolerance);
    assert(qmomkit::internal::residual_norm(system, structured.solution, rhs) <= residual_tolerance);
}

void assert_wrapper_matches_dense_reference(const std::vector<double>& abscissas, const std::vector<double>& weights,
                                            const qmomkit::internal::DenseMatrix& expected_solution, const qmomkit::VandermondeSolver solver,
                                            const qmomkit::VandermondeSolver expected_selected_solver,
                                            const qmomkit::internal::VandermondeAlgorithm expected_algorithm_used, const double solution_tolerance)
{
    const qmomkit::internal::DenseMatrix rhs = make_rhs(abscissas, weights, expected_solution);
    const qmomkit::internal::DenseMatrix system = qmomkit::internal::make_weighted_vandermonde(abscissas, weights);
    const qmomkit::internal::LinearSolveResult dense = qmomkit::internal::solve_dense(system, rhs, singular_tolerance, singular_message);
    assert(dense.success);

    const qmomkit::internal::VandermondeSolveOptions options{solver, singular_tolerance, true, 1.0e-10};
    const qmomkit::internal::VandermondeSolveResult result =
        qmomkit::internal::solve_weighted_vandermonde(abscissas, weights, rhs, options, singular_message);
    assert(result.success);
    assert(result.requested_solver == solver);
    assert(result.selected_solver == expected_selected_solver);
    assert(result.algorithm_used == expected_algorithm_used);
    assert(!result.dense_fallback_occurred);
    assert_matrix_near(result.solution, dense.solution, solution_tolerance);
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
    const qmomkit::internal::LinearSolveResult structured =
        qmomkit::internal::solve_weighted_transposed_vandermonde_bjorck_pereyra(abscissas, weights, rhs, singular_tolerance, singular_message);
    assert(structured.success);
    assert(structured.solution.rows() == abscissas.size());
    assert(structured.solution.cols() == 0);
}

void test_large_dynamic_range_rhs()
{
    assert_matches_dense_reference_relative(
        {-0.75, -0.1, 0.4, 0.9}, {1.0e-3, -2.0, 5.0e2, -1.0e5},
        make_solution({{1.0e-12, -1.0e12, 3.0}, {-2.0e-9, 4.0e9, -5.0}, {3.0e-6, -6.0e6, 7.0}, {-4.0e-3, 8.0e3, -9.0}}), 1.0e-2, 1.0e-10, 1.0e-1);
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
    assert_wrapper_matches_dense_reference(abscissas, weights, expected_solution, qmomkit::VandermondeSolver::bjorck_pereyra_with_dense_fallback,
                                           qmomkit::VandermondeSolver::bjorck_pereyra_with_dense_fallback,
                                           qmomkit::internal::VandermondeAlgorithm::bjorck_pereyra, 1.0e-11);
}

void test_wrapper_automatic_selects_bp_fallback_for_small_wide_system()
{
    const std::vector<double> abscissas = make_centered_nodes(8);
    const std::vector<double> weights = make_mild_weights(8);
    assert_wrapper_matches_dense_reference(abscissas, weights, make_patterned_solution(8, 256), qmomkit::VandermondeSolver::automatic,
                                           qmomkit::VandermondeSolver::bjorck_pereyra_with_dense_fallback,
                                           qmomkit::internal::VandermondeAlgorithm::bjorck_pereyra, 1.0e-9);
}

void test_wrapper_automatic_selects_bp_fallback_for_large_system()
{
    const std::vector<double> abscissas = make_centered_nodes(16);
    const std::vector<double> weights = make_mild_weights(16);
    assert_wrapper_matches_dense_reference(abscissas, weights, make_patterned_solution(16, 8), qmomkit::VandermondeSolver::automatic,
                                           qmomkit::VandermondeSolver::bjorck_pereyra_with_dense_fallback,
                                           qmomkit::internal::VandermondeAlgorithm::bjorck_pereyra, 1.0e-5);
}

void assert_singular_failure(const std::vector<double>& abscissas, const std::vector<double>& weights, const qmomkit::internal::DenseMatrix& rhs,
                             const std::string& expected_message)
{
    const qmomkit::internal::LinearSolveResult structured =
        qmomkit::internal::solve_weighted_transposed_vandermonde_bjorck_pereyra(abscissas, weights, rhs, singular_tolerance, singular_message);
    assert(!structured.success);
    assert(structured.message == expected_message);
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

void test_wrapper_bp_fallback_uses_dense_on_repeated_nodes()
{
    const std::vector<double> abscissas{1.0, 1.0};
    const std::vector<double> weights{0.5, 0.5};
    const qmomkit::internal::DenseMatrix rhs = make_rhs(abscissas, weights, make_solution({{1.0}, {2.0}}));
    const qmomkit::internal::VandermondeSolveOptions options{qmomkit::VandermondeSolver::bjorck_pereyra_with_dense_fallback, singular_tolerance, true,
                                                             1.0e-10};
    const qmomkit::internal::VandermondeSolveResult result =
        qmomkit::internal::solve_weighted_vandermonde(abscissas, weights, rhs, options, singular_message);
    assert(!result.success);
    assert(result.requested_solver == qmomkit::VandermondeSolver::bjorck_pereyra_with_dense_fallback);
    assert(result.selected_solver == qmomkit::VandermondeSolver::bjorck_pereyra_with_dense_fallback);
    assert(result.algorithm_used == qmomkit::internal::VandermondeAlgorithm::dense);
    assert(result.dense_fallback_occurred);
    assert(result.message == singular_message);
}

void test_invalid_inputs_fail()
{
    assert_singular_failure({}, {}, qmomkit::internal::DenseMatrix(0, 1), "Weighted Vandermonde system must not be empty.");
    assert_singular_failure({1.0, 2.0}, {1.0}, qmomkit::internal::DenseMatrix(2, 1), "Weighted Vandermonde dimensions are inconsistent.");
    assert_singular_failure({1.0, std::numeric_limits<double>::infinity()}, {1.0, 1.0}, qmomkit::internal::DenseMatrix(2, 1),
                            "Weighted Vandermonde inputs must be finite.");
}

} // namespace

int main()
{
    test_one_node_success();
    test_well_conditioned_uniform_weights();
    test_mixed_sign_abscissas();
    test_clustered_nodes_match_dense_reference();
    test_nearly_repeated_nodes_match_dense_reference();
    test_larger_node_count_smoke_cases();
    test_zero_rhs_columns_success();
    test_large_dynamic_range_rhs();
    test_wrapper_forced_modes_match_dense_reference();
    test_wrapper_automatic_selects_bp_fallback_for_small_wide_system();
    test_wrapper_automatic_selects_bp_fallback_for_large_system();
    test_repeated_nodes_fail();
    test_nonuniform_weights();
    test_small_weight_failure();
    test_zero_weight_failure();
    test_wrapper_bp_fallback_uses_dense_on_repeated_nodes();
    test_invalid_inputs_fail();

    return 0;
}
