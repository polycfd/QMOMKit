#include "vandermonde.hpp"

#include "bjorck_pereyra.hpp"
#include "householder_qr.hpp"
#include "jacobi_svd.hpp"
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
                                                        const VandermondeSolver selected_solver, const VandermondeAlgorithm algorithm_used)
{
    return VandermondeSolveResult{std::move(result.solution), result.success,  std::move(result.message),
                                  requested_solver,           selected_solver, algorithm_used};
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

} // namespace

DenseMatrix make_conditional_vandermonde_matrix(const std::vector<double>& nodes, const std::vector<double>& weights)
{
    const std::size_t n = nodes.size();
    DenseMatrix matrix(n, n);
    for (std::size_t power = 0; power < n; ++power)
    {
        for (std::size_t node = 0; node < n; ++node)
        {
            matrix(power, node) = std::pow(nodes[node], static_cast<int>(power)) * weights[node];
        }
    }
    return matrix;
}

VandermondeSolveResult solve_conditional_vandermonde_system(const std::vector<double>& nodes, const std::vector<double>& weights, DenseMatrix rhs,
                                                            const VandermondeSolveOptions& options, const std::string& singular_message)
{
    const std::size_t n = nodes.size();
    const std::size_t rhs_cols = rhs.cols();
    const VandermondeSolver requested_solver = options.solver;
    const VandermondeSolver selected_solver = resolve_solver(requested_solver);

    if (selected_solver == VandermondeSolver::dense)
    {
        performance::record_vandermonde_solve(performance::VandermondeSolveKind::dense);
        if (weights.size() != n || rhs.rows() != n)
        {
            return vandermonde_result(failure_result(n, rhs_cols, "Weighted Vandermonde dimensions are inconsistent."), requested_solver,
                                      selected_solver, VandermondeAlgorithm::dense);
        }
        DenseMatrix system = make_conditional_vandermonde_matrix(nodes, weights);
        return vandermonde_result(solve_dense(std::move(system), std::move(rhs), options.singular_tolerance, singular_message), requested_solver,
                                  selected_solver, VandermondeAlgorithm::dense);
    }

    if (selected_solver == VandermondeSolver::pivoted_qr)
    {
        performance::record_vandermonde_solve(performance::VandermondeSolveKind::pivoted_qr);
        if (weights.size() != n || rhs.rows() != n)
        {
            return vandermonde_result(failure_result(n, rhs_cols, "Weighted Vandermonde dimensions are inconsistent."), requested_solver,
                                      selected_solver, VandermondeAlgorithm::pivoted_qr);
        }

        const DenseMatrix system = make_conditional_vandermonde_matrix(nodes, weights);
        LinearSolveResult pivoted_qr = qr_linear_solve_result(solve_qr(system, rhs), n, rhs_cols, singular_message);
        if (!pivoted_qr.success)
        {
            return vandermonde_result(std::move(pivoted_qr), requested_solver, selected_solver, VandermondeAlgorithm::pivoted_qr);
        }

        if (options.validate_residual && !residual_is_accepted(system, pivoted_qr.solution, rhs, options.residual_tolerance))
        {
            return vandermonde_result(failure_result(n, rhs_cols, singular_message), requested_solver, selected_solver,
                                      VandermondeAlgorithm::pivoted_qr);
        }

        return vandermonde_result(std::move(pivoted_qr), requested_solver, selected_solver, VandermondeAlgorithm::pivoted_qr);
    }

    if (selected_solver == VandermondeSolver::svd)
    {
        performance::record_vandermonde_solve(performance::VandermondeSolveKind::svd);
        if (weights.size() != n || rhs.rows() != n)
        {
            return vandermonde_result(failure_result(n, rhs_cols, "Weighted Vandermonde dimensions are inconsistent."), requested_solver,
                                      selected_solver, VandermondeAlgorithm::svd);
        }

        const DenseMatrix system = make_conditional_vandermonde_matrix(nodes, weights);
        LinearSolveResult svd = svd_linear_solve_result(solve_svd(system, rhs), n, rhs_cols, singular_message);
        if (!svd.success)
        {
            return vandermonde_result(std::move(svd), requested_solver, selected_solver, VandermondeAlgorithm::svd);
        }

        if (options.validate_residual && !residual_is_accepted(system, svd.solution, rhs, options.residual_tolerance))
        {
            return vandermonde_result(failure_result(n, rhs_cols, singular_message), requested_solver, selected_solver, VandermondeAlgorithm::svd);
        }

        return vandermonde_result(std::move(svd), requested_solver, selected_solver, VandermondeAlgorithm::svd);
    }

    if (selected_solver == VandermondeSolver::moore_penrose)
    {
        performance::record_vandermonde_solve(performance::VandermondeSolveKind::svd);
        if (weights.size() != n || rhs.rows() != n)
        {
            return vandermonde_result(failure_result(n, rhs_cols, "Weighted Vandermonde dimensions are inconsistent."), requested_solver,
                                      selected_solver, VandermondeAlgorithm::moore_penrose);
        }

        const DenseMatrix system = make_conditional_vandermonde_matrix(nodes, weights);
        LinearSolveResult moore_penrose =
            svd_linear_solve_result(solve_svd_pseudoinverse(factorize_svd(system), rhs, options.moore_penrose_rcond), n, rhs_cols, singular_message);
        if (!moore_penrose.success)
        {
            return vandermonde_result(std::move(moore_penrose), requested_solver, selected_solver, VandermondeAlgorithm::moore_penrose);
        }

        if (options.validate_residual && !residual_is_accepted(system, moore_penrose.solution, rhs, options.residual_tolerance))
        {
            return vandermonde_result(failure_result(n, rhs_cols, singular_message), requested_solver, selected_solver,
                                      VandermondeAlgorithm::moore_penrose);
        }

        return vandermonde_result(std::move(moore_penrose), requested_solver, selected_solver, VandermondeAlgorithm::moore_penrose);
    }

    performance::record_vandermonde_solve(performance::VandermondeSolveKind::bjorck_pereyra);
    LinearSolveResult bjorck_pereyra{DenseMatrix(n, rhs_cols), true, ""};
    if (n == 0)
    {
        bjorck_pereyra = failure_result(n, rhs_cols, "Weighted Vandermonde system must not be empty.");
    }
    else if (weights.size() != n || rhs.rows() != n)
    {
        bjorck_pereyra = failure_result(n, rhs_cols, "Weighted Vandermonde dimensions are inconsistent.");
    }
    else if (!std::isfinite(options.singular_tolerance) || options.singular_tolerance <= 0.0)
    {
        bjorck_pereyra = failure_result(n, rhs_cols, "Weighted Vandermonde singular tolerance is invalid.");
    }
    else if (!qmomkit::internal::all_finite(nodes) || !qmomkit::internal::all_finite(weights) || !matrix_all_finite(rhs))
    {
        bjorck_pereyra = failure_result(n, rhs_cols, "Weighted Vandermonde inputs must be finite.");
    }

    if (bjorck_pereyra.success)
    {
        // A weight counts as zero relative to the largest weight, so the test
        // does not depend on the number density the weights carry.
        double weight_scale = 0.0;
        for (const double weight : weights)
        {
            weight_scale = std::max(weight_scale, std::abs(weight));
        }
        for (const double weight : weights)
        {
            if (std::abs(weight) <= options.singular_tolerance * weight_scale)
            {
                bjorck_pereyra = failure_result(n, rhs_cols, singular_message);
                break;
            }
        }
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

    return vandermonde_result(std::move(bjorck_pereyra), requested_solver, selected_solver, VandermondeAlgorithm::bjorck_pereyra);
}

} // namespace qmomkit::internal
