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
    bjorck_pereyra
};

struct VandermondeSolveOptions
{
    VandermondeSolver solver = VandermondeSolver::automatic;
    double singular_tolerance = 1.0e-14;
    bool validate_residual = false;
    double residual_tolerance = 1.0e-10;
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

} // namespace qmomkit::internal
