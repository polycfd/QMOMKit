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
    bool dense_fallback_occurred = false;
};

[[nodiscard]] LinearSolveResult solve_weighted_transposed_vandermonde_bjorck_pereyra(const std::vector<double>& abscissas,
                                                                                     const std::vector<double>& weights, DenseMatrix rhs,
                                                                                     double singular_tolerance, const std::string& singular_message);

[[nodiscard]] VandermondeSolveResult solve_weighted_vandermonde(const std::vector<double>& abscissas, const std::vector<double>& weights,
                                                                DenseMatrix rhs, const VandermondeSolveOptions& options,
                                                                const std::string& singular_message);

} // namespace qmomkit::internal
