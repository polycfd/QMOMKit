#pragma once

#include "qmomkit/qmomkit.hpp"

#include <array>
#include <cstddef>
#include <optional>
#include <string>

namespace cqmom_2d_borderline_conditioning
{

struct SolverChoice
{
    qmomkit::VandermondeSolver solver = qmomkit::VandermondeSolver::dense;
    const char* label = "dense";
};

struct MatrixDiagnostics
{
    double sigma_max = 0.0;
    double sigma_min = 0.0;
    double ratio = 0.0;
    std::size_t rank = 0;
    std::optional<double> condition_number;
};

struct RunRecord
{
    std::string mode;
    double delta = 0.0;
    std::string conditioning_order;
    std::string requested_solver;
    std::string selected_solver;
    bool success = false;
    std::string status;
    std::string active_counts;
    std::string diagnostic_level;
    std::string parent_path;
    std::optional<double> condition_estimate_inf;
    MatrixDiagnostics matrix;
    std::optional<double> residual_norm;
    std::optional<double> moment_error;
    std::string failed_dimension;
    std::string message;
};

inline const std::array<SolverChoice, 4> solver_choices{{
    {qmomkit::VandermondeSolver::dense, "dense"},
    {qmomkit::VandermondeSolver::bjorck_pereyra, "bjorck_pereyra"},
    {qmomkit::VandermondeSolver::pivoted_qr, "pivoted_qr"},
    {qmomkit::VandermondeSolver::svd, "svd"},
}};

inline const std::array<std::size_t, 4> report_solver_indices{{0, 1, 2, 3}};

inline const std::array<double, 8> delta_values{{1.0e-1, 1.0e-6, 1.0e-10, 1.0e-12, 1.0e-13, 2.5e-14, 2.0e-14, 0.0}};

} // namespace cqmom_2d_borderline_conditioning
