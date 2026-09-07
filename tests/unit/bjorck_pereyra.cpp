#include "internal/bjorck_pereyra.hpp"

#include "internal/performance.hpp"
#include "test_case_dispatch.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace
{

constexpr double singular_tolerance = 1.0e-14;

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

qmomkit::internal::DenseMatrix make_transposed_vandermonde_rhs(const std::vector<double>& nodes, const qmomkit::internal::DenseMatrix& solution)
{
    const std::size_t n = nodes.size();
    qmomkit::internal::DenseMatrix rhs(n, solution.cols());
    for (std::size_t power = 0; power < n; ++power)
    {
        for (std::size_t col = 0; col < solution.cols(); ++col)
        {
            double value = 0.0;
            for (std::size_t node = 0; node < n; ++node)
            {
                value += std::pow(nodes[node], static_cast<int>(power)) * solution(node, col);
            }
            rhs(power, col) = value;
        }
    }
    return rhs;
}

void assert_near(const double actual, const double expected, const double tolerance)
{
    assert(std::abs(actual - expected) <= tolerance);
}

void assert_close(const double actual, const double expected, const double absolute_tolerance, const double relative_tolerance)
{
    assert(std::abs(actual - expected) <= absolute_tolerance + relative_tolerance * std::abs(expected));
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

void assert_kernel_recovers_solution(const std::vector<double>& nodes, const qmomkit::internal::DenseMatrix& expected_solution,
                                     const double tolerance)
{
    const qmomkit::internal::DenseMatrix rhs = make_transposed_vandermonde_rhs(nodes, expected_solution);
    const qmomkit::internal::LinearSolveResult result = qmomkit::internal::solve_bjorck_pereyra(nodes, rhs, singular_tolerance);
    assert(result.success);
    assert_matrix_near(result.solution, expected_solution, tolerance);
}

void test_exact_transposed_vandermonde_solve()
{
    assert_kernel_recovers_solution({0.25, 1.5, 3.0}, make_solution({{1.0}, {-2.0}, {4.0}}), 1.0e-13);
}

void test_multiple_rhs_columns()
{
    assert_kernel_recovers_solution({0.25, 1.0, 2.0}, make_solution({{1.0, -2.0, 3.0}, {0.5, 4.0, -1.0}, {-3.0, 2.5, 0.0}}), 1.0e-12);
}

void test_zero_rhs_columns_success()
{
    const std::vector<double> nodes{-1.0, 0.5, 2.0};
    const qmomkit::internal::DenseMatrix rhs(nodes.size(), 0);
    const qmomkit::internal::LinearSolveResult result = qmomkit::internal::solve_bjorck_pereyra(nodes, rhs, singular_tolerance);
    assert(result.success);
    assert(result.solution.rows() == nodes.size());
    assert(result.solution.cols() == 0);
}

void test_mixed_sign_nodes()
{
    assert_kernel_recovers_solution({-2.0, 0.0, 1.5, 4.0}, make_solution({{1.0, -2.0}, {0.5, 4.0}, {-3.0, 2.5}, {2.0, -1.0}}), 1.0e-11);
}

void test_clustered_distinct_nodes()
{
    const std::vector<double> nodes{1.0, 1.0 + 1.0e-5, 1.0 + 3.0e-5};
    const qmomkit::internal::DenseMatrix expected = make_solution({{1.0, 0.25}, {-0.5, 2.0}, {3.0, -4.0}});
    const qmomkit::internal::DenseMatrix rhs = make_transposed_vandermonde_rhs(nodes, expected);
    const qmomkit::internal::LinearSolveResult result = qmomkit::internal::solve_bjorck_pereyra(nodes, rhs, singular_tolerance);
    assert(result.success);
    assert_matrix_close(result.solution, expected, 2.0e-5, 1.0e-12);
}

void test_repeated_nodes_fail()
{
    const qmomkit::internal::DenseMatrix rhs = make_solution({{1.0, 2.0}, {3.0, 4.0}});
    const qmomkit::internal::LinearSolveResult result = qmomkit::internal::solve_bjorck_pereyra({1.0, 1.0}, rhs, singular_tolerance);
    assert(!result.success);
    assert(result.message == "Transposed Vandermonde system is singular.");
}

void test_invalid_inputs_fail()
{
    qmomkit::internal::LinearSolveResult result =
        qmomkit::internal::solve_bjorck_pereyra({}, qmomkit::internal::DenseMatrix(0, 1), singular_tolerance);
    assert(!result.success);
    assert(result.message == "Transposed Vandermonde system must not be empty.");

    result = qmomkit::internal::solve_bjorck_pereyra({1.0, 2.0}, qmomkit::internal::DenseMatrix(3, 1), singular_tolerance);
    assert(!result.success);
    assert(result.message == "Transposed Vandermonde dimensions are inconsistent.");

    result = qmomkit::internal::solve_bjorck_pereyra({1.0, 2.0}, qmomkit::internal::DenseMatrix(2, 1), 0.0);
    assert(!result.success);
    assert(result.message == "Transposed Vandermonde singular tolerance is invalid.");

    result = qmomkit::internal::solve_bjorck_pereyra({1.0, std::numeric_limits<double>::infinity()}, qmomkit::internal::DenseMatrix(2, 1),
                                                     singular_tolerance);
    assert(!result.success);
    assert(result.message == "Transposed Vandermonde inputs must be finite.");

    qmomkit::internal::DenseMatrix rhs(2, 1);
    rhs(0, 0) = 1.0;
    rhs(1, 0) = std::numeric_limits<double>::quiet_NaN();
    result = qmomkit::internal::solve_bjorck_pereyra({1.0, 2.0}, rhs, singular_tolerance);
    assert(!result.success);
    assert(result.message == "Transposed Vandermonde inputs must be finite.");
}

qmomkit::internal::DenseMatrix make_benign_solution(const std::size_t n, const std::size_t rhs_cols)
{
    qmomkit::internal::DenseMatrix solution(n, rhs_cols);
    for (std::size_t row = 0; row < n; ++row)
    {
        for (std::size_t col = 0; col < rhs_cols; ++col)
        {
            solution(row, col) = static_cast<double>((row + 1) * (col + 1)) / static_cast<double>(n + rhs_cols);
        }
    }
    return solution;
}

std::vector<double> make_well_separated_nodes(const std::size_t n)
{
    std::vector<double> nodes(n);
    for (std::size_t node = 0; node < n; ++node)
    {
        nodes[node] = static_cast<double>(node) / static_cast<double>(n);
    }
    return nodes;
}

void assert_bjorck_pereyra_transform_updates(const std::size_t n, const std::size_t rhs_cols)
{
    const std::vector<double> nodes = make_well_separated_nodes(n);
    const qmomkit::internal::DenseMatrix rhs = make_transposed_vandermonde_rhs(nodes, make_benign_solution(n, rhs_cols));

    qmomkit::internal::performance::CounterSession session;
    qmomkit::internal::LinearSolveResult result{qmomkit::internal::DenseMatrix(0, 0), false, ""};
    {
        qmomkit::internal::performance::CounterScope scope(session);
        result = qmomkit::internal::solve_bjorck_pereyra(nodes, rhs, singular_tolerance);
    }

    assert(result.success);
#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)
    assert(session.snapshot().linear_algebra.bjorck_pereyra_transform_updates == n * (n - 1) * rhs_cols);
#endif
}

void test_transform_update_count_scaling()
{
    constexpr std::size_t fixed_rhs_cols = 3;
    for (const std::size_t n : {std::size_t{4}, std::size_t{8}, std::size_t{16}, std::size_t{32}})
    {
        assert_bjorck_pereyra_transform_updates(n, fixed_rhs_cols);
    }

    assert_bjorck_pereyra_transform_updates(8, 1);
    assert_bjorck_pereyra_transform_updates(8, 5);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(argc, argv,
                                             {
                                                 {"exact_transposed_vandermonde_solve", test_exact_transposed_vandermonde_solve},
                                                 {"multiple_rhs_columns", test_multiple_rhs_columns},
                                                 {"zero_rhs_columns_success", test_zero_rhs_columns_success},
                                                 {"mixed_sign_nodes", test_mixed_sign_nodes},
                                                 {"clustered_distinct_nodes", test_clustered_distinct_nodes},
                                                 {"repeated_nodes_fail", test_repeated_nodes_fail},
                                                 {"invalid_inputs_fail", test_invalid_inputs_fail},
                                                 {"transform_update_count_scaling", test_transform_update_count_scaling},
                                             });
}
