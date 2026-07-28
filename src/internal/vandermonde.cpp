#include "vandermonde.hpp"

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

[[nodiscard]] bool result_is_finite(const LinearSolveResult& result)
{
    return !result.success || matrix_all_finite(result.solution);
}

[[nodiscard]] LinearSolveResult failure_result(const std::size_t rows, const std::size_t cols, const std::string& message)
{
    return LinearSolveResult{DenseMatrix(rows, cols), false, message};
}

[[nodiscard]] VandermondeSolver select_automatic_solver()
{
    // Automatic favors the structured Vandermonde solve, while retaining the
    // dense path as a residual-checked fallback for difficult node sets.
    return VandermondeSolver::bjorck_pereyra_with_dense_fallback;
}

[[nodiscard]] VandermondeSolver resolve_solver(const VandermondeSolver requested_solver, const std::size_t n, const std::size_t rhs_cols)
{
    switch (requested_solver)
    {
    case VandermondeSolver::automatic:
        (void)n;
        (void)rhs_cols;
        return select_automatic_solver();
    case VandermondeSolver::dense:
    case VandermondeSolver::bjorck_pereyra:
    case VandermondeSolver::bjorck_pereyra_with_dense_fallback:
        return requested_solver;
    }
    return VandermondeSolver::dense;
}

[[nodiscard]] VandermondeSolveResult vandermonde_result(LinearSolveResult result, const VandermondeSolver requested_solver,
                                                        const VandermondeSolver selected_solver, const VandermondeAlgorithm algorithm_used,
                                                        const bool dense_fallback_occurred)
{
    return VandermondeSolveResult{std::move(result.solution), result.success,  std::move(result.message),
                                  requested_solver,           selected_solver, algorithm_used,
                                  dense_fallback_occurred};
}

[[nodiscard]] double node_scale(const double lhs, const double rhs)
{
    return std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

[[nodiscard]] bool denominator_is_singular(const double denominator, const double lhs, const double rhs, const double singular_tolerance)
{
    return std::abs(denominator) <= singular_tolerance * node_scale(lhs, rhs);
}

[[nodiscard]] bool residual_is_accepted(const DenseMatrix& system, const DenseMatrix& solution, const DenseMatrix& rhs, const double tolerance)
{
    const double residual = residual_norm(system, solution, rhs);
    if (!std::isfinite(residual))
    {
        return false;
    }
    const double scale = std::max({1.0, matrix_infinity_norm(system) * matrix_infinity_norm(solution), matrix_infinity_norm(rhs)});
    return residual <= tolerance * scale;
}

[[nodiscard]] bool residual_validation_required(const VandermondeSolver requested_solver, const VandermondeSolver selected_solver)
{
    return requested_solver == VandermondeSolver::automatic || selected_solver == VandermondeSolver::bjorck_pereyra_with_dense_fallback;
}

} // namespace

LinearSolveResult solve_weighted_transposed_vandermonde_bjorck_pereyra(const std::vector<double>& abscissas, const std::vector<double>& weights,
                                                                       DenseMatrix rhs, const double singular_tolerance,
                                                                       const std::string& singular_message)
{
    const std::size_t n = abscissas.size();
    const std::size_t rhs_cols = rhs.cols();
    if (n == 0)
    {
        return failure_result(n, rhs_cols, "Weighted Vandermonde system must not be empty.");
    }
    if (weights.size() != n || rhs.rows() != n)
    {
        return failure_result(n, rhs_cols, "Weighted Vandermonde dimensions are inconsistent.");
    }
    if (!std::isfinite(singular_tolerance) || singular_tolerance <= 0.0)
    {
        return failure_result(n, rhs_cols, "Weighted Vandermonde singular tolerance is invalid.");
    }
    if (!qmomkit::internal::all_finite(abscissas) || !qmomkit::internal::all_finite(weights) || !matrix_all_finite(rhs))
    {
        return failure_result(n, rhs_cols, "Weighted Vandermonde inputs must be finite.");
    }

    // In CQMOM this system maps mixed parent/target moments to conditional
    // target moments C at each parent node:
    //     (V^T diag(w)) C = RHS.
    double weight_scale = 1.0;
    for (const double weight : weights)
    {
        weight_scale = std::max(weight_scale, std::abs(weight));
    }
    for (const double weight : weights)
    {
        if (std::abs(weight) <= singular_tolerance * weight_scale)
        {
            return failure_result(n, rhs_cols, singular_message);
        }
    }

    // Apply the transposed elementary operations of the Björck-Pereyra
    // Vandermonde solve. The result is Y in V^T Y = rhs.
    for (std::size_t k = 0; k + 1 < n; ++k)
    {
        for (std::size_t reverse_i = n - 1; reverse_i > k; --reverse_i)
        {
            const std::size_t i = reverse_i - 1;
            for (std::size_t col = 0; col < rhs_cols; ++col)
            {
                rhs(i + 1, col) -= abscissas[k] * rhs(i, col);
            }
        }
    }

    for (std::size_t reverse_k = n - 1; reverse_k > 0; --reverse_k)
    {
        const std::size_t k = reverse_k - 1;
        for (std::size_t i = k + 1; i < n; ++i)
        {
            const std::size_t denominator_node = i - k - 1;
            const double denominator = abscissas[i] - abscissas[denominator_node];
            if (denominator_is_singular(denominator, abscissas[i], abscissas[denominator_node], singular_tolerance))
            {
                return failure_result(n, rhs_cols, singular_message);
            }
            for (std::size_t col = 0; col < rhs_cols; ++col)
            {
                const double value = rhs(i, col) / denominator;
                rhs(i - 1, col) -= value;
                rhs(i, col) = value;
            }
        }
    }

    for (std::size_t row = 0; row < n; ++row)
    {
        // The Björck-Pereyra stages solve V^T Y = RHS; dividing by weights
        // converts Y into the conditional-moment unknowns C.
        for (std::size_t col = 0; col < rhs_cols; ++col)
        {
            rhs(row, col) /= weights[row];
        }
    }

    if (!matrix_all_finite(rhs))
    {
        return failure_result(n, rhs_cols, "Weighted Vandermonde solve produced non-finite values.");
    }

    return LinearSolveResult{std::move(rhs), true, ""};
}

VandermondeSolveResult solve_weighted_vandermonde(const std::vector<double>& abscissas, const std::vector<double>& weights, DenseMatrix rhs,
                                                  const VandermondeSolveOptions& options, const std::string& singular_message)
{
    const VandermondeSolver requested_solver = options.solver;
    const VandermondeSolver selected_solver = resolve_solver(requested_solver, abscissas.size(), rhs.cols());

    // Solver choice is part of diagnostics: forced dense uses Gaussian
    // elimination, forced Björck-Pereyra reports that result directly, and
    // automatic/fallback validates the structured solve before using dense.
    if (selected_solver == VandermondeSolver::dense)
    {
        performance::record_vandermonde_solve(performance::VandermondeSolveKind::dense);
        DenseMatrix system = make_weighted_vandermonde(abscissas, weights);
        return vandermonde_result(solve_dense(std::move(system), std::move(rhs), options.singular_tolerance, singular_message), requested_solver,
                                  selected_solver, VandermondeAlgorithm::dense, false);
    }

    performance::record_vandermonde_solve(performance::VandermondeSolveKind::bjorck_pereyra);
    LinearSolveResult bjorck_pereyra =
        solve_weighted_transposed_vandermonde_bjorck_pereyra(abscissas, weights, rhs, options.singular_tolerance, singular_message);
    if (bjorck_pereyra.success && result_is_finite(bjorck_pereyra))
    {
        const bool validate_residual = options.validate_residual || residual_validation_required(requested_solver, selected_solver);
        if (!validate_residual)
        {
            return vandermonde_result(bjorck_pereyra, requested_solver, selected_solver, VandermondeAlgorithm::bjorck_pereyra, false);
        }
        const DenseMatrix system = make_weighted_vandermonde(abscissas, weights);
        if (residual_is_accepted(system, bjorck_pereyra.solution, rhs, options.residual_tolerance))
        {
            return vandermonde_result(bjorck_pereyra, requested_solver, selected_solver, VandermondeAlgorithm::bjorck_pereyra, false);
        }
    }

    if (selected_solver == VandermondeSolver::bjorck_pereyra)
    {
        return vandermonde_result(bjorck_pereyra, requested_solver, selected_solver, VandermondeAlgorithm::bjorck_pereyra, false);
    }

    DenseMatrix system = make_weighted_vandermonde(abscissas, weights);
    performance::record_vandermonde_fallback();
    performance::record_vandermonde_solve(performance::VandermondeSolveKind::dense);
    return vandermonde_result(solve_dense(std::move(system), std::move(rhs), options.singular_tolerance, singular_message), requested_solver,
                              selected_solver, VandermondeAlgorithm::dense, true);
}

} // namespace qmomkit::internal
