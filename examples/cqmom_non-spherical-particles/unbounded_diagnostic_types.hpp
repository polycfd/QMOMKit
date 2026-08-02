#pragma once

#include "qmomkit/qmomkit.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace cqmom_non_spherical_particles
{

struct Matrix
{
    std::size_t rows = 0;
    std::size_t cols = 0;
    std::vector<double> values;

    Matrix() = default;
    Matrix(std::size_t row_count, std::size_t col_count);

    double& operator()(std::size_t row, std::size_t col);
    double operator()(std::size_t row, std::size_t col) const;
};

struct EigenResult
{
    std::vector<double> eigenvalues;
    Matrix eigenvectors;
};

struct MatrixDiagnostic
{
    std::string family;
    std::size_t k = 0;
    std::size_t l = 0;
    Matrix matrix;
    double determinant = 0.0;
    std::vector<double> eigenvalues;
    double minimum_eigenvalue = 0.0;
    double infinity_norm = 0.0;
    double effective_tolerance = 0.0;
    double signed_margin = 0.0;
    std::string classification;
};

struct ConditionalSolveDiagnostic
{
    std::vector<double> parent_nodes;
    std::vector<double> parent_weights;
    std::vector<double> parent_global_weights;
    Matrix vandermonde;
    Matrix rhs;
    Matrix solved;
    double condition_estimate = 0.0;
    double residual_norm = 0.0;
    double maximum_input_scale = 0.0;
    double maximum_output_scale = 0.0;
};

struct StageMoments
{
    std::string label;
    double time = 0.0;
    std::vector<double> values;
    std::vector<double> branch_moments;
};

struct SourceReplay
{
    std::optional<qmomkit::CqmomMomentLayout> layout;
    std::vector<StageMoments> stages;
    double max_difference_from_public_source_update = 0.0;
};

struct OffendingNode
{
    std::size_t terminal_node = 0;
    qmomkit::ConditioningPath path;
    std::vector<double> state;
    double terminal_weight = 0.0;
    double parent_global_weight = 0.0;
    double chi_conditional_weight = 0.0;
    std::vector<double> chi_nodes;
    std::vector<double> chi_conditional_weights;
    std::vector<double> chi_global_weights;
};

struct BoundedAcceptedBranch
{
    qmomkit::ConditioningPath parent_path;
    std::vector<double> chi_nodes;
    std::vector<double> conditional_weights;
    std::vector<double> global_weights;
    std::vector<double> canonical_moments;
    std::vector<double> retained_prefix_absolute_residuals;
    std::vector<double> retained_prefix_relative_residuals;
    double parent_global_weight = 0.0;
};

struct BoundaryAdvance
{
    bool ok = false;
    std::string message;
    std::optional<qmomkit::CqmomState> pre_source_state;
    std::optional<qmomkit::CqmomQuadrature> pre_source_quadrature;
    std::optional<qmomkit::CqmomState> boundary_state;
    std::optional<qmomkit::CqmomQuadrature> boundary_quadrature;
    std::optional<qmomkit::CqmomState> final_state;
    std::optional<qmomkit::CqmomQuadrature> final_quadrature;
    qmomkit::TransportDiagnostics source_diagnostics;
    qmomkit::TransportStage failed_stage = qmomkit::TransportStage::none;
    double source_start_time = 0.0;
    double boundary_time = 0.0;
    std::string stage_label;
};

struct DiagnosticReport
{
    std::size_t macro_step = 0;
    double time = 0.0;
    std::string operator_stage_label;
    OffendingNode offending;
    qmomkit::CqmomDiagnostics unbounded_branch;
    std::optional<qmomkit::CqmomDiagnostics> bounded_branch;
    std::optional<qmomkit::CqmomNodeCountTransition> bounded_transition;
    std::optional<BoundedAcceptedBranch> bounded_accepted_branch;
    ConditionalSolveDiagnostic conditional_solve;
    SourceReplay source_replay;
    std::vector<double> branch_moments;
    std::vector<double> reconstructed_moments;
    qmomkit::Quadrature1D wheeler;
    double configured_rmin = 0.0;
    std::vector<double> recurrence_a;
    std::vector<double> recurrence_b;
    Matrix sigma;
    std::vector<MatrixDiagnostic> matrices;
};

} // namespace cqmom_non_spherical_particles
