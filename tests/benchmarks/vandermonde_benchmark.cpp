#include "../conditional_vandermonde_matrix.hpp"
#include "internal/bjorck_pereyra.hpp"
#include "internal/linear_algebra.hpp"
#include "internal/vandermonde.hpp"
#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <benchmark/benchmark.h>
#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace
{

constexpr double singular_tolerance = 1.0e-300;
constexpr const char* singular_message = "Weighted Vandermonde system is singular.";

enum class NumericalCase
{
    separated = 0,
    mixed_sign = 1,
    clustered = 2,
    tiny_weights = 3
};

struct VandermondeFixture
{
    std::vector<double> abscissas;
    std::vector<double> weights;
    qmomkit::internal::DenseMatrix rhs;
    qmomkit::internal::DenseMatrix system;
};

std::vector<double> make_abscissas(const std::size_t n, const NumericalCase numerical_case)
{
    std::vector<double> abscissas(n, 0.0);
    if (n == 1)
    {
        abscissas.front() = numerical_case == NumericalCase::mixed_sign ? -0.25 : 0.25;
        return abscissas;
    }

    for (std::size_t i = 0; i < n; ++i)
    {
        const double fraction = static_cast<double>(i) / static_cast<double>(n - 1);
        switch (numerical_case)
        {
        case NumericalCase::separated:
        case NumericalCase::tiny_weights:
            abscissas[i] = 0.05 + 0.90 * fraction;
            break;
        case NumericalCase::mixed_sign:
            abscissas[i] = -0.90 + 1.80 * fraction;
            break;
        case NumericalCase::clustered:
            abscissas[i] = 1.0 + 1.0e-4 * static_cast<double>(i);
            break;
        }
    }
    return abscissas;
}

std::vector<double> make_weights(const std::size_t n, const NumericalCase numerical_case)
{
    std::vector<double> weights(n, 1.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        switch (numerical_case)
        {
        case NumericalCase::tiny_weights:
            weights[i] = i == n / 2 ? 1.0e-10 : 1.0 + 0.01 * static_cast<double>(i);
            break;
        case NumericalCase::mixed_sign:
        default:
            weights[i] = 0.75 + 0.05 * static_cast<double>((i % 7) + 1);
            break;
        }
    }
    return weights;
}

qmomkit::internal::DenseMatrix make_solution(const std::size_t rows, const std::size_t cols)
{
    qmomkit::internal::DenseMatrix solution(rows, cols);
    for (std::size_t row = 0; row < rows; ++row)
    {
        for (std::size_t col = 0; col < cols; ++col)
        {
            const double sign = (row + col) % 2 == 0 ? 1.0 : -1.0;
            solution(row, col) = sign * (0.125 + 0.01 * static_cast<double>(row + 1) + 0.0025 * static_cast<double>(col + 1));
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
            for (std::size_t node = 0; node < n; ++node)
            {
                value += weights[node] * std::pow(abscissas[node], static_cast<int>(order)) * solution(node, col);
            }
            rhs(order, col) = value;
        }
    }
    return rhs;
}

VandermondeFixture make_fixture(const std::size_t n, const std::size_t rhs_cols, const NumericalCase numerical_case)
{
    std::vector<double> abscissas = make_abscissas(n, numerical_case);
    std::vector<double> weights = make_weights(n, numerical_case);
    qmomkit::internal::DenseMatrix rhs = make_rhs(abscissas, weights, make_solution(n, rhs_cols));
    qmomkit::internal::DenseMatrix system = qmomkit::test::make_conditional_vandermonde_matrix(abscissas, weights);
    return VandermondeFixture{std::move(abscissas), std::move(weights), std::move(rhs), std::move(system)};
}

qmomkit::internal::LinearSolveResult solve_weighted_bjorck_pereyra_for_benchmark(const std::vector<double>& abscissas,
                                                                                 const std::vector<double>& weights,
                                                                                 qmomkit::internal::DenseMatrix rhs)
{
    double weight_scale = 1.0;
    for (const double weight : weights)
    {
        weight_scale = std::max(weight_scale, std::abs(weight));
    }
    for (const double weight : weights)
    {
        if (std::abs(weight) <= singular_tolerance * weight_scale)
        {
            return qmomkit::internal::LinearSolveResult{qmomkit::internal::DenseMatrix(abscissas.size(), rhs.cols()), false, singular_message};
        }
    }

    qmomkit::internal::LinearSolveResult result = qmomkit::internal::solve_bjorck_pereyra(abscissas, std::move(rhs), singular_tolerance);
    if (!result.success)
    {
        return result;
    }
    for (std::size_t row = 0; row < weights.size(); ++row)
    {
        for (std::size_t col = 0; col < result.solution.cols(); ++col)
        {
            result.solution(row, col) /= weights[row];
        }
    }
    return result;
}

bool skip_if_invalid(benchmark::State& state, const VandermondeFixture& fixture)
{
    const qmomkit::internal::LinearSolveResult dense =
        qmomkit::internal::solve_dense(fixture.system, fixture.rhs, singular_tolerance, singular_message);
    const qmomkit::internal::LinearSolveResult structured =
        solve_weighted_bjorck_pereyra_for_benchmark(fixture.abscissas, fixture.weights, fixture.rhs);
    if (!dense.success || !structured.success)
    {
        state.SkipWithError("Vandermonde benchmark fixture is singular under the configured tolerance.");
        return true;
    }
    return false;
}

void run_dense_solve(benchmark::State& state, const std::size_t n, const std::size_t rhs_cols, const NumericalCase numerical_case)
{
    const VandermondeFixture fixture = make_fixture(n, rhs_cols, numerical_case);
    if (skip_if_invalid(state, fixture))
    {
        return;
    }

    for (auto _ : state)
    {
        qmomkit::internal::LinearSolveResult result =
            qmomkit::internal::solve_dense(fixture.system, fixture.rhs, singular_tolerance, singular_message);
        benchmark::DoNotOptimize(result.solution.rows());
        benchmark::DoNotOptimize(result.solution.cols());
        benchmark::DoNotOptimize(result.success);
    }
}

void set_scaling_counters(benchmark::State& state, const std::size_t n, const std::size_t rhs_cols)
{
    state.counters["N"] = static_cast<double>(n);
    state.counters["RHS"] = static_cast<double>(rhs_cols);
    state.counters["N3"] = static_cast<double>(n) * static_cast<double>(n) * static_cast<double>(n);
    state.counters["N2_RHS"] = static_cast<double>(n) * static_cast<double>(n) * static_cast<double>(rhs_cols);
}

void set_solver_selection_counters(benchmark::State& state, const qmomkit::internal::VandermondeSolveResult& result)
{
    state.counters["selected_dense"] = result.selected_solver == qmomkit::VandermondeSolver::dense ? 1.0 : 0.0;
    state.counters["selected_qr"] = result.selected_solver == qmomkit::VandermondeSolver::pivoted_qr ? 1.0 : 0.0;
    state.counters["selected_svd"] = result.selected_solver == qmomkit::VandermondeSolver::svd ? 1.0 : 0.0;
    state.counters["selected_bp"] = result.selected_solver == qmomkit::VandermondeSolver::bjorck_pereyra ? 1.0 : 0.0;
    state.counters["used_dense"] = result.algorithm_used == qmomkit::internal::VandermondeAlgorithm::dense ? 1.0 : 0.0;
    state.counters["used_qr"] = result.algorithm_used == qmomkit::internal::VandermondeAlgorithm::pivoted_qr ? 1.0 : 0.0;
    state.counters["used_svd"] = result.algorithm_used == qmomkit::internal::VandermondeAlgorithm::svd ? 1.0 : 0.0;
    state.counters["used_bp"] = result.algorithm_used == qmomkit::internal::VandermondeAlgorithm::bjorck_pereyra ? 1.0 : 0.0;
}

void run_bjorck_pereyra_solve(benchmark::State& state, const std::size_t n, const std::size_t rhs_cols, const NumericalCase numerical_case)
{
    const VandermondeFixture fixture = make_fixture(n, rhs_cols, numerical_case);
    if (skip_if_invalid(state, fixture))
    {
        return;
    }

    for (auto _ : state)
    {
        qmomkit::internal::LinearSolveResult result = solve_weighted_bjorck_pereyra_for_benchmark(fixture.abscissas, fixture.weights, fixture.rhs);
        benchmark::DoNotOptimize(result.solution.rows());
        benchmark::DoNotOptimize(result.solution.cols());
        benchmark::DoNotOptimize(result.success);
    }
}

void run_dense_scaling_solve(benchmark::State& state, const std::size_t n, const std::size_t rhs_cols, const NumericalCase numerical_case)
{
    set_scaling_counters(state, n, rhs_cols);
    run_dense_solve(state, n, rhs_cols, numerical_case);
}

void run_bjorck_pereyra_scaling_solve(benchmark::State& state, const std::size_t n, const std::size_t rhs_cols, const NumericalCase numerical_case)
{
    set_scaling_counters(state, n, rhs_cols);
    run_bjorck_pereyra_solve(state, n, rhs_cols, numerical_case);
}

void run_weighted_wrapper_solve(benchmark::State& state, const std::size_t n, const std::size_t rhs_cols, const NumericalCase numerical_case,
                                const qmomkit::VandermondeSolver solver)
{
    set_scaling_counters(state, n, rhs_cols);
    const VandermondeFixture fixture = make_fixture(n, rhs_cols, numerical_case);
    if (skip_if_invalid(state, fixture))
    {
        return;
    }

    const qmomkit::internal::VandermondeSolveOptions options{solver, singular_tolerance};
    const qmomkit::internal::VandermondeSolveResult setup_check =
        qmomkit::internal::solve_conditional_vandermonde_system(fixture.abscissas, fixture.weights, fixture.rhs, options, singular_message);
    if (!setup_check.success)
    {
        state.SkipWithError("Vandermonde solver benchmark fixture failed setup.");
        return;
    }
    set_solver_selection_counters(state, setup_check);

    for (auto _ : state)
    {
        const qmomkit::internal::VandermondeSolveResult result =
            qmomkit::internal::solve_conditional_vandermonde_system(fixture.abscissas, fixture.weights, fixture.rhs, options, singular_message);
        benchmark::DoNotOptimize(result.solution.rows());
        benchmark::DoNotOptimize(result.solution.cols());
    }
}

void run_direct_solve(benchmark::State& state, const std::size_t n, const std::size_t rhs_cols, const NumericalCase numerical_case,
                      const qmomkit::VandermondeSolver solver)
{
    switch (solver)
    {
    case qmomkit::VandermondeSolver::dense:
        run_dense_scaling_solve(state, n, rhs_cols, numerical_case);
        return;
    case qmomkit::VandermondeSolver::bjorck_pereyra:
        run_bjorck_pereyra_scaling_solve(state, n, rhs_cols, numerical_case);
        return;
    case qmomkit::VandermondeSolver::pivoted_qr:
    case qmomkit::VandermondeSolver::svd:
    case qmomkit::VandermondeSolver::moore_penrose:
    case qmomkit::VandermondeSolver::automatic:
        run_weighted_wrapper_solve(state, n, rhs_cols, numerical_case, solver);
        return;
    }
}

void bm_vandermonde_scaling_matrix_size(benchmark::State& state, const qmomkit::VandermondeSolver solver)
{
    run_direct_solve(state, static_cast<std::size_t>(state.range(0)), 16, NumericalCase::separated, solver);
}

void bm_vandermonde_scaling_rhs_width(benchmark::State& state, const qmomkit::VandermondeSolver solver)
{
    run_direct_solve(state, 32, static_cast<std::size_t>(state.range(0)), NumericalCase::separated, solver);
}

void bm_vandermonde_robustness(benchmark::State& state, const NumericalCase numerical_case, const qmomkit::VandermondeSolver solver)
{
    run_direct_solve(state, 16, 16, numerical_case, solver);
}

void bm_vandermonde_solver_selection(benchmark::State& state, const qmomkit::VandermondeSolver solver)
{
    run_weighted_wrapper_solve(state, static_cast<std::size_t>(state.range(0)), static_cast<std::size_t>(state.range(1)), NumericalCase::separated,
                               solver);
}

BENCHMARK_CAPTURE(bm_vandermonde_scaling_matrix_size, dense, qmomkit::VandermondeSolver::dense)
    ->ArgName("N")
    ->Args({4})
    ->Args({8})
    ->Args({16})
    ->Args({32})
    ->Args({64})
    ->Args({128});
BENCHMARK_CAPTURE(bm_vandermonde_scaling_matrix_size, bjorck_pereyra, qmomkit::VandermondeSolver::bjorck_pereyra)
    ->ArgName("N")
    ->Args({4})
    ->Args({8})
    ->Args({16})
    ->Args({32})
    ->Args({64})
    ->Args({128});

BENCHMARK_CAPTURE(bm_vandermonde_scaling_rhs_width, dense, qmomkit::VandermondeSolver::dense)->ArgName("RHS")->Args({1})->Args({16})->Args({256});
BENCHMARK_CAPTURE(bm_vandermonde_scaling_rhs_width, bjorck_pereyra, qmomkit::VandermondeSolver::bjorck_pereyra)
    ->ArgName("RHS")
    ->Args({1})
    ->Args({16})
    ->Args({256});

BENCHMARK_CAPTURE(bm_vandermonde_robustness, clustered_nodes_dense, NumericalCase::clustered, qmomkit::VandermondeSolver::dense);
BENCHMARK_CAPTURE(bm_vandermonde_robustness, clustered_nodes_bjorck_pereyra, NumericalCase::clustered, qmomkit::VandermondeSolver::bjorck_pereyra);
BENCHMARK_CAPTURE(bm_vandermonde_robustness, mixed_sign_nodes_dense, NumericalCase::mixed_sign, qmomkit::VandermondeSolver::dense);
BENCHMARK_CAPTURE(bm_vandermonde_robustness, mixed_sign_nodes_bjorck_pereyra, NumericalCase::mixed_sign, qmomkit::VandermondeSolver::bjorck_pereyra);
BENCHMARK_CAPTURE(bm_vandermonde_robustness, tiny_weights_dense, NumericalCase::tiny_weights, qmomkit::VandermondeSolver::dense);
BENCHMARK_CAPTURE(bm_vandermonde_robustness, tiny_weights_bjorck_pereyra, NumericalCase::tiny_weights, qmomkit::VandermondeSolver::bjorck_pereyra);

BENCHMARK_CAPTURE(bm_vandermonde_solver_selection, dense, qmomkit::VandermondeSolver::dense)
    ->ArgNames({"N", "RHS"})
    ->Args({8, 256})
    ->Args({16, 8})
    ->Args({32, 64})
    ->Args({64, 256});
BENCHMARK_CAPTURE(bm_vandermonde_solver_selection, bjorck_pereyra, qmomkit::VandermondeSolver::bjorck_pereyra)
    ->ArgNames({"N", "RHS"})
    ->Args({8, 256})
    ->Args({16, 8})
    ->Args({32, 64})
    ->Args({64, 256});
BENCHMARK_CAPTURE(bm_vandermonde_solver_selection, automatic, qmomkit::VandermondeSolver::automatic)
    ->ArgNames({"N", "RHS"})
    ->Args({8, 256})
    ->Args({16, 8})
    ->Args({32, 64})
    ->Args({64, 256});

} // namespace
