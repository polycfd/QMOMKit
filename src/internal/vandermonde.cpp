#include "vandermonde.hpp"

#include "bjorck_pereyra.hpp"
#include "householder_qr.hpp"
#include "jacobi_svd.hpp"
#include "numeric_utils.hpp"
#include "performance.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
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

[[nodiscard]] LinearSolveResult failure_result(const std::size_t rows, const std::size_t cols, const std::string& message)
{
    return LinearSolveResult{DenseMatrix(rows, cols), false, message};
}

[[nodiscard]] VandermondeSolver select_automatic_solver()
{
    return VandermondeSolver::bjorck_pereyra;
}

[[nodiscard]] VandermondeSolver resolve_solver(const VandermondeSolver requested_solver)
{
    switch (requested_solver)
    {
    case VandermondeSolver::automatic:
        return select_automatic_solver();
    case VandermondeSolver::dense:
    case VandermondeSolver::pivoted_qr:
    case VandermondeSolver::svd:
    case VandermondeSolver::moore_penrose:
    case VandermondeSolver::bjorck_pereyra:
        return requested_solver;
    }
    return VandermondeSolver::dense;
}

[[nodiscard]] VandermondeSolveResult vandermonde_result(LinearSolveResult result, const VandermondeSolver requested_solver,
                                                        const VandermondeSolver selected_solver)
{
    return VandermondeSolveResult{std::move(result.solution), result.success, std::move(result.message), requested_solver, selected_solver};
}

[[nodiscard]] LinearSolveResult qr_linear_solve_result(QrSolveResult result, const std::size_t rows, const std::size_t cols,
                                                       const std::string& singular_message)
{
    if (result.status == QrStatus::success)
    {
        return LinearSolveResult{std::move(result.solution), true, ""};
    }

    if (result.status == QrStatus::rank_deficient)
    {
        return failure_result(rows, cols, singular_message);
    }

    const std::string message = result.message.empty() ? singular_message : result.message;
    return failure_result(rows, cols, message);
}

[[nodiscard]] LinearSolveResult svd_linear_solve_result(SvdSolveResult result, const std::size_t rows, const std::size_t cols,
                                                        const std::string& singular_message)
{
    if (result.status == SvdStatus::success)
    {
        return LinearSolveResult{std::move(result.solution), true, ""};
    }

    if (result.status == SvdStatus::rank_deficient)
    {
        return failure_result(rows, cols, singular_message);
    }

    const std::string message = result.message.empty() ? singular_message : result.message;
    return failure_result(rows, cols, message);
}

// The checks every solver path applies before solving. A weight counts as zero
// relative to the largest weight, so the test does not depend on the number
// density the weights carry.
[[nodiscard]] std::optional<LinearSolveResult> weighted_vandermonde_input_failure(const std::vector<double>& nodes,
                                                                                  const std::vector<double>& weights, const DenseMatrix& rhs,
                                                                                  const double singular_tolerance,
                                                                                  const std::string& singular_message)
{
    const std::size_t n = nodes.size();
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
    if (!qmomkit::internal::all_finite(nodes) || !qmomkit::internal::all_finite(weights) || !matrix_all_finite(rhs))
    {
        return failure_result(n, rhs_cols, "Weighted Vandermonde inputs must be finite.");
    }
    double weight_scale = 0.0;
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
    return std::nullopt;
}

// A(k, j) = w_j x_j^k factors as S V D with S = diag(s^k), D = diag(w) and V the
// unweighted Vandermonde matrix of x / s, where s is the power of two above the
// largest |x|. The direct solvers solve V Y = S^-1 R and then take C = D^-1 Y,
// separating the weights as the Bjorck-Pereyra path does. Their singularity,
// rank and truncation decisions then depend on neither the number density nor
// the unit of the nodes: the first row of V is all ones and every entry lies in
// [-1, 1]. Both scalings are by powers of two and therefore exact.
struct ScaledVandermondeSystem
{
    DenseMatrix system;
    DenseMatrix rhs;
};

[[nodiscard]] ScaledVandermondeSystem make_scaled_vandermonde_system(const std::vector<double>& nodes, DenseMatrix rhs)
{
    const std::size_t n = nodes.size();
    double largest_node = 0.0;
    for (const double node : nodes)
    {
        largest_node = std::max(largest_node, std::abs(node));
    }
    int exponent = 0;
    if (largest_node > 0.0)
    {
        static_cast<void>(std::frexp(largest_node, &exponent));
    }

    DenseMatrix system(n, n);
    for (std::size_t node = 0; node < n; ++node)
    {
        const double scaled_node = std::ldexp(nodes[node], -exponent);
        double power = 1.0;
        for (std::size_t row = 0; row < n; ++row)
        {
            system(row, node) = power;
            power *= scaled_node;
        }
    }
    for (std::size_t row = 0; row < n; ++row)
    {
        for (std::size_t col = 0; col < rhs.cols(); ++col)
        {
            rhs(row, col) = std::ldexp(rhs(row, col), -exponent * static_cast<int>(row));
        }
    }
    return ScaledVandermondeSystem{std::move(system), std::move(rhs)};
}

// C = D^-1 Y. A solution that overflows on the division is a failure.
[[nodiscard]] LinearSolveResult divide_by_weights(LinearSolveResult result, const std::vector<double>& weights)
{
    if (!result.success)
    {
        return result;
    }
    for (std::size_t row = 0; row < result.solution.rows(); ++row)
    {
        for (std::size_t col = 0; col < result.solution.cols(); ++col)
        {
            result.solution(row, col) /= weights[row];
        }
    }
    if (!matrix_all_finite(result.solution))
    {
        return failure_result(result.solution.rows(), result.solution.cols(), "Weighted Vandermonde solve produced non-finite values.");
    }
    return result;
}

} // namespace

ConditionalVandermondeDiagnostics conditional_vandermonde_diagnostics(const std::vector<double>& nodes, const std::vector<double>& weights,
                                                                      DenseMatrix rhs, const DenseMatrix& solution, const double singular_tolerance)
{
    const ScaledVandermondeSystem scaled = make_scaled_vandermonde_system(nodes, std::move(rhs));
    DenseMatrix scaled_solution = solution;
    for (std::size_t row = 0; row < scaled_solution.rows(); ++row)
    {
        for (std::size_t col = 0; col < scaled_solution.cols(); ++col)
        {
            scaled_solution(row, col) *= weights[row];
        }
    }
    const double residual = residual_norm(scaled.system, scaled_solution, scaled.rhs);
    const double scale = std::max(matrix_infinity_norm(scaled.system) * matrix_infinity_norm(scaled_solution), matrix_infinity_norm(scaled.rhs));
    return ConditionalVandermondeDiagnostics{scale > 0.0 ? residual / scale : residual, condition_estimate(scaled.system, singular_tolerance)};
}

VandermondeSolveResult solve_conditional_vandermonde_system(const std::vector<double>& nodes, const std::vector<double>& weights, DenseMatrix rhs,
                                                            const VandermondeSolveOptions& options, const std::string& singular_message)
{
    const std::size_t n = nodes.size();
    const std::size_t rhs_cols = rhs.cols();
    const VandermondeSolver requested_solver = options.solver;
    const VandermondeSolver selected_solver = resolve_solver(requested_solver);

    // The direct solvers share the input checks and the scaled system; only the
    // kernel differs.
    const auto solve_direct = [&](const auto& solve_scaled_system)
    {
        if (std::optional<LinearSolveResult> failure =
                weighted_vandermonde_input_failure(nodes, weights, rhs, options.singular_tolerance, singular_message))
        {
            return vandermonde_result(std::move(*failure), requested_solver, selected_solver);
        }
        const ScaledVandermondeSystem scaled = make_scaled_vandermonde_system(nodes, std::move(rhs));
        return vandermonde_result(divide_by_weights(solve_scaled_system(scaled), weights), requested_solver, selected_solver);
    };

    if (selected_solver == VandermondeSolver::dense)
    {
        performance::record_vandermonde_solve(performance::VandermondeSolveKind::dense);
        return solve_direct([&](const ScaledVandermondeSystem& scaled)
                            { return solve_dense(scaled.system, scaled.rhs, options.singular_tolerance, singular_message); });
    }

    if (selected_solver == VandermondeSolver::pivoted_qr)
    {
        performance::record_vandermonde_solve(performance::VandermondeSolveKind::pivoted_qr);
        return solve_direct([&](const ScaledVandermondeSystem& scaled)
                            { return qr_linear_solve_result(solve_qr(scaled.system, scaled.rhs), n, rhs_cols, singular_message); });
    }

    if (selected_solver == VandermondeSolver::svd)
    {
        performance::record_vandermonde_solve(performance::VandermondeSolveKind::svd);
        return solve_direct([&](const ScaledVandermondeSystem& scaled)
                            { return svd_linear_solve_result(solve_svd(scaled.system, scaled.rhs), n, rhs_cols, singular_message); });
    }

    if (selected_solver == VandermondeSolver::moore_penrose)
    {
        performance::record_vandermonde_solve(performance::VandermondeSolveKind::svd);
        // The truncated minimum-norm solution is minimal in the scaled unknowns Y,
        // and rcond is relative to the singular values of the scaled system.
        return solve_direct(
            [&](const ScaledVandermondeSystem& scaled)
            {
                return svd_linear_solve_result(solve_svd_pseudoinverse(factorize_svd(scaled.system), scaled.rhs, options.moore_penrose_rcond), n,
                                               rhs_cols, singular_message);
            });
    }

    performance::record_vandermonde_solve(performance::VandermondeSolveKind::bjorck_pereyra);
    LinearSolveResult bjorck_pereyra{DenseMatrix(n, rhs_cols), true, ""};
    if (std::optional<LinearSolveResult> failure =
            weighted_vandermonde_input_failure(nodes, weights, rhs, options.singular_tolerance, singular_message))
    {
        bjorck_pereyra = std::move(*failure);
    }
    if (bjorck_pereyra.success)
    {
        // `rhs` is dead after this call: the remaining work reads only the
        // returned solution, and a Bjorck-Pereyra failure becomes a failure
        // result rather than a retry with the original right-hand side. Moving
        // hands the kernel the caller's matrix instead of copying it, which is
        // one n-by-rhs_cols allocation and copy per conditional solve.
        LinearSolveResult kernel_result = solve_bjorck_pereyra(nodes, std::move(rhs), options.singular_tolerance);
        bjorck_pereyra = kernel_result.success ? std::move(kernel_result) : failure_result(n, rhs_cols, singular_message);
    }

    if (bjorck_pereyra.success)
    {
        for (std::size_t row = 0; row < n; ++row)
        {
            // The BP kernel solves V^T Y = RHS; dividing by weights converts Y
            // into the conditional-moment unknowns C for A*C = RHS.
            for (std::size_t col = 0; col < rhs_cols; ++col)
            {
                bjorck_pereyra.solution(row, col) /= weights[row];
            }
        }

        if (!matrix_all_finite(bjorck_pereyra.solution))
        {
            bjorck_pereyra = failure_result(n, rhs_cols, "Weighted Vandermonde solve produced non-finite values.");
        }
    }

    return vandermonde_result(std::move(bjorck_pereyra), requested_solver, selected_solver);
}

} // namespace qmomkit::internal
