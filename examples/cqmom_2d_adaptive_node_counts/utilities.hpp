#pragma once

#include "qmomkit/conditioning_path.hpp"
#include "qmomkit/cqmom.hpp"
#include "qmomkit/moment_reconstruction.hpp"
#include "qmomkit/phase_space_evolution.hpp"

#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cqmom_2d_adaptive_node_counts
{

struct MomentDiagnostics
{
    double mass = 0.0;
    double size_mean = 0.0;
    double size_variance = 0.0;
    double velocity_mean = 0.0;
    double velocity_variance = 0.0;
};

struct VelocityBranchId
{
    std::size_t conditioning_level = 0;
    std::size_t phase_space_dimension = 0;
    qmomkit::ConditioningPath parent_path;
};

struct VelocityBranchCount
{
    VelocityBranchId id;
    std::size_t configured_maximum = 0;
    std::size_t active_count = 0;
};

struct HistorySample
{
    double time = 0.0;
    bool source_active = false;
    MomentDiagnostics moments;
    std::vector<VelocityBranchCount> velocity_branch_counts;
    qmomkit::CqmomQuadrature quadrature;
};

struct ActiveCountTransition
{
    double time = 0.0;
    VelocityBranchId branch;
    std::size_t configured_maximum = 0;
    std::optional<std::size_t> before;
    std::optional<std::size_t> after;
    qmomkit::CqmomNodeCountTransitionKind kind = qmomkit::CqmomNodeCountTransitionKind::grew;
    std::vector<double> velocity_moments;
    long double rank_indicator = 0.0L;
};

void print_configuration_errors(const std::string& error);

void ensure_output_directory(const std::filesystem::path& path);

[[nodiscard]] HistorySample make_history_sample(double time, bool source_active, const qmomkit::CqmomQuadrature& quadrature,
                                                const std::vector<qmomkit::CqmomActiveBranchCount>& active_branch_counts,
                                                std::size_t velocity_dimension);

void append_active_count_transitions(double time, const qmomkit::CqmomQuadrature& previous_quadrature,
                                     const qmomkit::CqmomQuadrature& current_quadrature,
                                     const std::vector<qmomkit::CqmomNodeCountTransition>& native_transitions,
                                     std::vector<ActiveCountTransition>& transitions, std::size_t velocity_dimension);

[[nodiscard]] bool validate_history(const std::vector<HistorySample>& history, double initial_velocity_mean, double relaxation_rate,
                                    double equilibrium_velocity, double pulse_start, double pulse_end, std::string& message);

[[nodiscard]] bool validate_required_velocity_moments(const qmomkit::CqmomMomentLayout& layout, std::size_t velocity_dimension, std::string& message);

[[nodiscard]] bool validate_active_count_lifecycle(const std::vector<HistorySample>& history, const std::vector<ActiveCountTransition>& transitions,
                                                   double pulse_start, double pulse_end, std::string& message);

void print_result_summary(const qmomkit::SelectedRawMomentResult& initial_moments, const qmomkit::SelectedRawMomentResult& final_moments,
                          const std::vector<HistorySample>& output_history, const std::vector<ActiveCountTransition>& active_count_transitions,
                          const qmomkit::CqmomQuadrature& final_quadrature, const std::filesystem::path& output_directory,
                          const std::vector<std::size_t>& conditioning_order, const std::vector<std::size_t>& node_counts, double relaxation_rate,
                          double equilibrium_velocity, double pulse_start, double pulse_end, double source_rate, double velocity_displacement,
                          double time_step, double final_time, std::size_t required_velocity_order, const std::string& source_validation_message,
                          const std::string& moment_layout_message, bool wrote_moments, bool wrote_quadrature);

} // namespace cqmom_2d_adaptive_node_counts
