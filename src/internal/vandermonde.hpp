#pragma once

#include "linear_algebra.hpp"
#include "qmomkit/vandermonde_solver.hpp"

#include <string>
#include <vector>

namespace qmomkit::internal
{

enum class VandermondeAlgorithm
{
    dense,
    pivoted_qr,
    svd,
    moore_penrose,
    bjorck_pereyra
};

struct VandermondeSolveOptions
{
    VandermondeSolver solver = VandermondeSolver::automatic;
    double singular_tolerance = 1.0e-14;
    double moore_penrose_rcond = 1.0e-14;
};

struct VandermondeSolveResult
{
    DenseMatrix solution;
    bool success = true;
    std::string message;
    VandermondeSolver requested_solver = VandermondeSolver::automatic;
    VandermondeSolver selected_solver = VandermondeSolver::dense;
    VandermondeAlgorithm algorithm_used = VandermondeAlgorithm::dense;
};

/**
 * Construct the conditional Vandermonde matrix A = V^T diag(w), with
 * A(power, node) = weights[node] * nodes[node]^power.
 */
[[nodiscard]] DenseMatrix make_conditional_vandermonde_matrix(const std::vector<double>& nodes, const std::vector<double>& weights);

/**
 * Solve the complete conditional system A*C = rhs, where
 * A(power, node) = weights[node] * nodes[node]^power.
 */
[[nodiscard]] VandermondeSolveResult solve_conditional_vandermonde_system(const std::vector<double>& nodes, const std::vector<double>& weights,
                                                                          DenseMatrix rhs, const VandermondeSolveOptions& options,
                                                                          const std::string& singular_message);

struct ConditionalVandermondeDiagnostics
{
    double relative_residual = 0.0;
    double condition_estimate = 0.0;
};

/**
 * Diagnostics of a successful solve of A*C = rhs, evaluated on the system the
 * direct solvers decide on, V*Y = S^-1*rhs with A = S*V*D and Y = D*C: the
 * infinity-norm condition estimate of V, and the Frobenius norm of
 * V*Y - S^-1*rhs relative to max(||V|| ||Y||, ||S^-1*rhs||) in the infinity
 * norm. Neither depends on the number density, nor on a change of the unit of
 * the nodes by a power of two. Another change of unit changes the condition
 * estimate by less than a factor 2^(n-1), because the nodes are normalized
 * only up to a factor of two.
 */
[[nodiscard]] ConditionalVandermondeDiagnostics conditional_vandermonde_diagnostics(const std::vector<double>& nodes,
                                                                                    const std::vector<double>& weights, DenseMatrix rhs,
                                                                                    const DenseMatrix& solution, double singular_tolerance);

} // namespace qmomkit::internal
