#include "utilities.hpp"

#include "qmomkit/qmom.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <system_error>

namespace cqmom_2d_adaptive_node_counts
{
namespace
{

constexpr std::size_t max_console_moments = 12;
constexpr std::size_t source_rate_checked_order = 5;

void print_section_heading(const char* title)
{
    std::cout << '\n' << title << '\n';
    for (const char* character = title; *character != '\0'; ++character)
    {
        std::cout << '-';
    }
    std::cout << '\n';
}

std::string format_multi_index(const qmomkit::MultiIndex& order)
{
    std::ostringstream stream;
    stream << '[';
    for (std::size_t dimension = 0; dimension < order.size(); ++dimension)
    {
        if (dimension > 0)
        {
            stream << ',';
        }
        stream << order[dimension];
    }
    stream << ']';
    return stream.str();
}

std::string format_path(const qmomkit::ConditioningPath& path)
{
    if (path.empty())
    {
        return "root";
    }
    std::ostringstream stream;
    stream << '[';
    for (std::size_t level = 0; level < path.size(); ++level)
    {
        if (level > 0)
        {
            stream << ',';
        }
        stream << path.node_index_at_level(level);
    }
    stream << ']';
    return stream.str();
}

std::string format_terminal_path(const qmomkit::ConditioningPathView path)
{
    std::ostringstream stream;
    stream << '[';
    for (std::size_t level = 0; level < path.size(); ++level)
    {
        if (level > 0)
        {
            stream << ',';
        }
        stream << path.node_index_at_level(level);
    }
    stream << ']';
    return stream.str();
}

std::string format_terminal_paths(const qmomkit::ConditioningPathTable& paths)
{
    std::ostringstream stream;
    for (std::size_t path = 0; path < paths.size(); ++path)
    {
        if (path > 0)
        {
            stream << ' ';
        }
        stream << format_terminal_path(paths[path]);
    }
    return stream.str();
}

std::string format_index_vector(const std::vector<std::size_t>& values)
{
    std::ostringstream stream;
    stream << '[';
    for (std::size_t value = 0; value < values.size(); ++value)
    {
        if (value > 0)
        {
            stream << ',';
        }
        stream << values[value];
    }
    stream << ']';
    return stream.str();
}

void print_selected_moments(const char* label, const qmomkit::SelectedRawMomentResult& moments)
{
    print_section_heading(label);
    if (!moments.ok())
    {
        std::cout << "  unavailable: " << moments.message << '\n';
        return;
    }

    const std::size_t printed_count = std::min(max_console_moments, moments.values.size());
    std::size_t label_width = 0;
    for (std::size_t index = 0; index < printed_count; ++index)
    {
        label_width = std::max(label_width, std::string("M" + format_multi_index(moments.indices[index])).size());
    }
    for (std::size_t index = 0; index < printed_count; ++index)
    {
        const std::string moment_label = "M" + format_multi_index(moments.indices[index]);
        std::cout << std::left << std::setw(static_cast<int>(label_width + 2)) << moment_label << std::right << std::scientific
                  << std::setprecision(6) << moments.values[index] << '\n';
    }
    if (moments.values.size() > printed_count)
    {
        std::cout << "... " << (moments.values.size() - printed_count) << " more moments omitted from the console summary\n";
    }
}

void print_quadrature(const qmomkit::CqmomQuadrature& quadrature)
{
    print_section_heading("Final quadrature");
    std::cout << " node  " << std::left << std::setw(10) << "path" << std::right << std::setw(16) << "weight" << std::setw(16) << "d" << std::setw(16)
              << "u" << '\n';
    for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
    {
        std::cout << std::setw(5) << node << "  " << std::left << std::setw(10) << format_terminal_path(quadrature.terminal_paths()[node])
                  << std::right << std::scientific << std::setprecision(6) << std::setw(16) << quadrature.weights()[node] << std::setw(16)
                  << quadrature.abscissas()[0][node] << std::setw(16) << quadrature.abscissas()[1][node] << '\n';
    }
}

double selected_moment(const qmomkit::SelectedRawMomentResult& moments, const qmomkit::MultiIndex& index)
{
    for (std::size_t value = 0; value < moments.values.size(); ++value)
    {
        if (moments.indices[value] == index)
        {
            return moments.values[value];
        }
    }
    throw std::runtime_error("Requested diagnostic moment is unavailable.");
}

MomentDiagnostics compute_moment_diagnostics(const qmomkit::CqmomQuadrature& quadrature)
{
    const qmomkit::RawMomentRequest request{{{0, 0}, {1, 0}, {2, 0}, {0, 1}, {0, 2}}};
    const qmomkit::SelectedRawMomentResult moments = qmomkit::reconstruct_raw_moments(quadrature, request);
    if (!moments.ok())
    {
        throw std::runtime_error("Could not reconstruct diagnostic moments: " + moments.message);
    }

    MomentDiagnostics diagnostics;
    diagnostics.mass = selected_moment(moments, {0, 0});
    diagnostics.size_mean = selected_moment(moments, {1, 0}) / diagnostics.mass;
    diagnostics.velocity_mean = selected_moment(moments, {0, 1}) / diagnostics.mass;
    diagnostics.size_variance = selected_moment(moments, {2, 0}) / diagnostics.mass - diagnostics.size_mean * diagnostics.size_mean;
    diagnostics.velocity_variance = selected_moment(moments, {0, 2}) / diagnostics.mass - diagnostics.velocity_mean * diagnostics.velocity_mean;
    return diagnostics;
}

bool branch_less(const VelocityBranchId& lhs, const VelocityBranchId& rhs)
{
    if (lhs.conditioning_level != rhs.conditioning_level)
    {
        return lhs.conditioning_level < rhs.conditioning_level;
    }
    if (lhs.phase_space_dimension != rhs.phase_space_dimension)
    {
        return lhs.phase_space_dimension < rhs.phase_space_dimension;
    }
    return lhs.parent_path < rhs.parent_path;
}

std::vector<VelocityBranchCount> velocity_branch_counts(const std::vector<qmomkit::CqmomActiveBranchCount>& active_branch_counts,
                                                        const std::size_t velocity_dimension)
{
    std::vector<VelocityBranchCount> counts;
    for (const qmomkit::CqmomActiveBranchCount& active : active_branch_counts)
    {
        if (active.phase_space_dimension == velocity_dimension)
        {
            counts.push_back(VelocityBranchCount{VelocityBranchId{active.conditioning_level, active.phase_space_dimension, active.parent_path},
                                                 active.configured_node_count, active.active_node_count});
        }
    }
    std::sort(counts.begin(), counts.end(),
              [](const VelocityBranchCount& lhs, const VelocityBranchCount& rhs) { return branch_less(lhs.id, rhs.id); });
    return counts;
}

std::vector<double> branch_velocity_moments(const qmomkit::CqmomQuadrature& quadrature, const VelocityBranchId& branch,
                                            const std::size_t maximum_order, const std::size_t velocity_dimension)
{
    std::vector<double> moments(maximum_order + 1, 0.0);
    for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
    {
        const qmomkit::ConditioningPathView path = quadrature.terminal_paths()[node];
        if (path.size() <= branch.conditioning_level)
        {
            continue;
        }

        bool same_parent = true;
        for (std::size_t level = 0; level < branch.conditioning_level; ++level)
        {
            if (path.node_index_at_level(level) != branch.parent_path.node_index_at_level(level))
            {
                same_parent = false;
                break;
            }
        }
        if (!same_parent)
        {
            continue;
        }

        double power = 1.0;
        const double velocity = quadrature.abscissas()[velocity_dimension][node];
        for (std::size_t order = 0; order <= maximum_order; ++order)
        {
            moments[order] += quadrature.weights()[node] * power;
            power *= velocity;
        }
    }
    return moments;
}

long double centered_rank_indicator(const std::vector<double>& moments)
{
    if (moments.size() < 5 || moments[0] == 0.0)
    {
        return 0.0L;
    }

    const long double m0 = static_cast<long double>(moments[0]);
    const long double mean = static_cast<long double>(moments[1]) / m0;
    long double centered[5] = {};
    for (std::size_t order = 0; order <= 4; ++order)
    {
        for (std::size_t k = 0; k <= order; ++k)
        {
            const long double sign = ((order - k) % 2 == 0) ? 1.0L : -1.0L;
            long double binomial = 1.0L;
            for (std::size_t i = 1; i <= k; ++i)
            {
                binomial *= static_cast<long double>(order + 1 - i);
                binomial /= static_cast<long double>(i);
            }
            centered[order] += sign * binomial * std::pow(mean, static_cast<int>(order - k)) * static_cast<long double>(moments[k]) / m0;
        }
    }
    return centered[2] * centered[4] - centered[3] * centered[3] - centered[2] * centered[2] * centered[2];
}

const char* qmom_status_name(const qmomkit::QmomStatus status)
{
    switch (status)
    {
    case qmomkit::QmomStatus::success:
        return "success";
    case qmomkit::QmomStatus::invalid_input:
        return "invalid_input";
    case qmomkit::QmomStatus::not_realizable:
        return "not_realizable";
    }
    return "unknown";
}

std::string format_branch_counts(const std::vector<VelocityBranchCount>& counts)
{
    std::ostringstream stream;
    for (std::size_t branch = 0; branch < counts.size(); ++branch)
    {
        if (branch > 0)
        {
            stream << ", ";
        }
        stream << format_path(counts[branch].id.parent_path) << ": " << counts[branch].active_count << '/' << counts[branch].configured_maximum;
    }
    return stream.str();
}

void print_variance_history(const std::vector<HistorySample>& history)
{
    print_section_heading("Variance and active-count history");
    std::cout << "Active counts are shown as parent-path: active/configured.\n";
    std::cout << std::right << std::setw(12) << "time" << "  " << std::left << std::setw(3) << "src" << "  " << std::right << std::setw(14)
              << "mean(u)" << std::setw(16) << "var(u)" << "  " << std::left << std::setw(28) << "active u counts" << "terminal paths\n";
    for (const HistorySample& sample : history)
    {
        std::cout << std::right << std::scientific << std::setprecision(3) << std::setw(12) << sample.time << "  " << std::left << std::setw(3)
                  << (sample.source_active ? "on" : "off") << "  " << std::right << std::scientific << std::setprecision(6) << std::setw(14)
                  << sample.moments.velocity_mean << std::setw(16) << sample.moments.velocity_variance << "  " << std::left << std::setw(28)
                  << format_branch_counts(sample.velocity_branch_counts) << format_terminal_paths(sample.quadrature.terminal_paths()) << '\n';
    }
}

void print_active_count_transitions(const std::vector<ActiveCountTransition>& transitions)
{
    print_section_heading("Active velocity-count transitions");
    if (transitions.empty())
    {
        std::cout << "none\n";
        return;
    }

    std::cout << std::right << std::setw(12) << "time" << "  " << std::left << std::setw(10) << "branch" << std::setw(10) << "change" << std::setw(10)
              << "kind" << std::right << std::setw(8) << "maximum" << '\n';
    for (const ActiveCountTransition& transition : transitions)
    {
        std::ostringstream change;
        if (transition.before.has_value())
        {
            change << *transition.before;
        }
        else
        {
            change << "created";
        }
        change << " -> ";
        if (transition.after.has_value())
        {
            change << *transition.after;
        }
        else
        {
            change << "removed";
        }
        const char* kind = transition.kind == qmomkit::CqmomNodeCountTransitionKind::reduced   ? "reduced"
                           : transition.kind == qmomkit::CqmomNodeCountTransitionKind::grew    ? "grew"
                           : transition.kind == qmomkit::CqmomNodeCountTransitionKind::created ? "created"
                                                                                               : "removed";
        std::cout << std::right << std::scientific << std::setprecision(3) << std::setw(12) << transition.time << "  " << std::left << std::setw(10)
                  << format_path(transition.branch.parent_path) << std::setw(10) << change.str() << std::setw(10) << kind << std::right
                  << std::setw(8) << transition.configured_maximum << '\n';
    }
}

void print_three_node_growth_diagnostics(const std::vector<ActiveCountTransition>& transitions)
{
    std::vector<ActiveCountTransition> growth_transitions;
    for (const ActiveCountTransition& transition : transitions)
    {
        if (transition.before == std::optional<std::size_t>{2} && transition.after == std::optional<std::size_t>{3})
        {
            growth_transitions.push_back(transition);
        }
    }
    if (growth_transitions.empty())
    {
        return;
    }

    print_section_heading("Three-node growth diagnostic");
    std::cout << "R2 > 0 indicates moment information beyond a two-node representation.\n";
    std::cout << std::right << std::setw(12) << "time" << "  " << std::left << std::setw(10) << "branch" << std::right << std::setw(14) << "R2"
              << "  " << std::left << std::setw(18) << "Wheeler ceiling 2" << "Wheeler ceiling 3\n";
    for (const ActiveCountTransition& transition : growth_transitions)
    {
        qmomkit::QmomOptions options;
        options.adaptive = true;
        const qmomkit::QmomResult current_ceiling =
            qmomkit::compute_qmom(qmomkit::MomentSet{transition.velocity_moments}, *transition.before, options);
        const qmomkit::QmomResult growth_ceiling = qmomkit::compute_qmom(qmomkit::MomentSet{transition.velocity_moments}, *transition.after, options);
        std::ostringstream current_status;
        current_status << qmom_status_name(current_ceiling.status) << ' ' << current_ceiling.abscissas.size();
        std::ostringstream growth_status;
        growth_status << qmom_status_name(growth_ceiling.status) << ' ' << growth_ceiling.abscissas.size();
        std::cout << std::right << std::scientific << std::setprecision(3) << std::setw(12) << transition.time << "  " << std::left << std::setw(10)
                  << format_path(transition.branch.parent_path) << std::right << std::scientific << std::setprecision(6) << std::setw(14)
                  << static_cast<double>(transition.rank_indicator) << "  " << std::left << std::setw(18) << current_status.str()
                  << growth_status.str() << '\n';

        std::cout << "  moments";
        for (std::size_t order = 0; order < transition.velocity_moments.size(); ++order)
        {
            std::cout << "  m" << order << '=' << std::scientific << std::setprecision(6) << transition.velocity_moments[order];
        }
        std::cout << '\n';
    }
}

} // namespace

void print_configuration_errors(const std::string& error)
{
    std::cerr << error << '\n';
}

void ensure_output_directory(const std::filesystem::path& path)
{
    std::error_code error;
    std::filesystem::create_directories(path, error);
    if (error)
    {
        throw std::runtime_error("Could not create output directory " + path.string() + ": " + error.message());
    }
}

HistorySample make_history_sample(const double time, const bool source_active, const qmomkit::CqmomQuadrature& quadrature,
                                  const std::vector<qmomkit::CqmomActiveBranchCount>& active_branch_counts, const std::size_t velocity_dimension)
{
    return HistorySample{time, source_active, compute_moment_diagnostics(quadrature),
                         velocity_branch_counts(active_branch_counts, velocity_dimension), quadrature};
}

void append_active_count_transitions(const double time, const qmomkit::CqmomQuadrature& previous_quadrature,
                                     const qmomkit::CqmomQuadrature& current_quadrature,
                                     const std::vector<qmomkit::CqmomNodeCountTransition>& native_transitions,
                                     std::vector<ActiveCountTransition>& transitions, const std::size_t velocity_dimension)
{
    for (const qmomkit::CqmomNodeCountTransition& native : native_transitions)
    {
        if (native.phase_space_dimension != velocity_dimension)
        {
            continue;
        }

        const VelocityBranchId branch{native.conditioning_level, native.phase_space_dimension, native.parent_path};
        const std::size_t diagnostic_order = std::max(std::size_t{4}, 2 * native.configured_maximum_node_count - 1);
        const qmomkit::CqmomQuadrature& diagnostic_quadrature = native.active_node_count_after.has_value() ? current_quadrature : previous_quadrature;
        const std::vector<double> moments = branch_velocity_moments(diagnostic_quadrature, branch, diagnostic_order, velocity_dimension);
        transitions.push_back(ActiveCountTransition{time, branch, native.configured_maximum_node_count, native.active_node_count_before,
                                                    native.active_node_count_after, native.kind, moments, centered_rank_indicator(moments)});
    }
}

bool validate_history(const std::vector<HistorySample>& history, const double initial_velocity_mean, const double relaxation_rate,
                      const double equilibrium_velocity, const double pulse_start, const double pulse_end, std::string& message)
{
    if (history.empty())
    {
        message = "No history samples were recorded.";
        return false;
    }

    const double initial_mass = history.front().moments.mass;
    double maximum_mass_change = 0.0;
    double maximum_mean_error = 0.0;
    double variance_before_pulse = std::numeric_limits<double>::quiet_NaN();
    double maximum_pulse_variance = -std::numeric_limits<double>::infinity();

    for (const HistorySample& sample : history)
    {
        maximum_mass_change = std::max(maximum_mass_change, std::abs(sample.moments.mass - initial_mass));
        const double expected_mean = equilibrium_velocity + (initial_velocity_mean - equilibrium_velocity) * std::exp(-relaxation_rate * sample.time);
        maximum_mean_error = std::max(maximum_mean_error, std::abs(sample.moments.velocity_mean - expected_mean));

        if (sample.time < pulse_start)
        {
            variance_before_pulse = sample.moments.velocity_variance;
        }
        if (sample.time >= pulse_start && sample.time <= pulse_end + 1.0e-12)
        {
            maximum_pulse_variance = std::max(maximum_pulse_variance, sample.moments.velocity_variance);
        }
    }

    if (maximum_mass_change > 1.0e-8)
    {
        message = "Total weight changed by more than 1e-8.";
        return false;
    }
    if (maximum_mean_error > 5.0e-5)
    {
        message = "Velocity mean deviated from the advection relaxation reference by more than 5e-5.";
        return false;
    }
    if (!std::isfinite(variance_before_pulse) || maximum_pulse_variance <= variance_before_pulse + 1.0e-4)
    {
        message = "Velocity variance did not increase clearly during the source pulse.";
        return false;
    }

    message = "mass conserved, velocity mean follows advection relaxation, and pulse increases velocity variance";
    return true;
}

bool validate_required_velocity_moments(const qmomkit::CqmomMomentLayout& layout, const std::size_t velocity_dimension, std::string& message)
{
    const std::size_t required_velocity_order = 2 * layout.hierarchy.node_count(velocity_dimension) - 1;
    for (std::size_t order = 0; order <= required_velocity_order; ++order)
    {
        qmomkit::MultiIndex index(layout.hierarchy.number_of_dimensions(), 0);
        index[velocity_dimension] = order;
        if (std::find(layout.required_indices.begin(), layout.required_indices.end(), index) == layout.required_indices.end())
        {
            std::ostringstream stream;
            stream << "The CQMOM source layout does not expose required velocity moment M";
            stream << '[';
            for (std::size_t dimension = 0; dimension < index.size(); ++dimension)
            {
                if (dimension > 0)
                {
                    stream << ',';
                }
                stream << index[dimension];
            }
            stream << "].";
            message = stream.str();
            return false;
        }
    }

    std::ostringstream stream;
    stream << "highest order: " << required_velocity_order;
    message = stream.str();
    return true;
}

bool validate_active_count_lifecycle(const std::vector<HistorySample>& history, const std::vector<ActiveCountTransition>& transitions,
                                     const double pulse_start, const double pulse_end, std::string& message)
{
    if (history.empty())
    {
        message = "No history samples were recorded.";
        return false;
    }

    bool saw_pre_pulse_reduction = false;
    bool saw_pulse_growth = false;
    bool saw_post_pulse_reduction = false;
    bool three_node_growth_required = false;
    bool saw_two_to_three_growth = false;
    constexpr double time_tolerance = 1.0e-12;

    for (const HistorySample& sample : history)
    {
        for (const VelocityBranchCount& branch : sample.velocity_branch_counts)
        {
            if (branch.configured_maximum >= 3)
            {
                three_node_growth_required = true;
            }
        }
    }

    for (const ActiveCountTransition& transition : transitions)
    {
        if (!transition.before.has_value() || !transition.after.has_value())
        {
            continue;
        }

        if (*transition.after > *transition.before)
        {
            if (transition.time + time_tolerance < pulse_start)
            {
                message = "A velocity branch grew before the variance-injection pulse.";
                return false;
            }
            if (transition.time > pulse_end + time_tolerance)
            {
                message = "A velocity branch grew after the variance-injection pulse.";
                return false;
            }
            if (*transition.after > std::min(*transition.before + 1, transition.configured_maximum))
            {
                message = "A velocity branch grew by more than one active node.";
                return false;
            }
            saw_pulse_growth = true;
            if (*transition.before == 2 && *transition.after == 3)
            {
                saw_two_to_three_growth = true;
            }
        }
        else if (*transition.after < *transition.before)
        {
            if (transition.time + time_tolerance < pulse_start)
            {
                saw_pre_pulse_reduction = true;
            }
            else if (transition.time > pulse_end + time_tolerance)
            {
                saw_post_pulse_reduction = true;
            }
        }
    }

    if (!saw_pre_pulse_reduction)
    {
        message = "No velocity branch adaptively reduced before the pulse.";
        return false;
    }
    if (!saw_pulse_growth)
    {
        message = "No collapsed velocity branch regrew during the pulse.";
        return false;
    }
    if (three_node_growth_required && !saw_two_to_three_growth)
    {
        message = "No velocity branch demonstrated controlled 2 -> 3 growth.";
        return false;
    }
    if (!saw_post_pulse_reduction)
    {
        message = "No velocity branch reduced after the pulse.";
        return false;
    }

    message = three_node_growth_required ? "velocity branches reduced, regrew by one node per step through three nodes, then reduced again"
                                         : "velocity branches reduced, regrew by one node per step, then reduced again";
    return true;
}

void print_result_summary(const qmomkit::SelectedRawMomentResult& initial_moments, const qmomkit::SelectedRawMomentResult& final_moments,
                          const std::vector<HistorySample>& output_history, const std::vector<ActiveCountTransition>& active_count_transitions,
                          const qmomkit::CqmomQuadrature& final_quadrature, const std::filesystem::path& output_directory,
                          const std::vector<std::size_t>& conditioning_order, const std::vector<std::size_t>& node_counts,
                          const double relaxation_rate, const double equilibrium_velocity, const double pulse_start, const double pulse_end,
                          const double source_rate, const double velocity_displacement, const double time_step, const double final_time,
                          const std::size_t required_velocity_order, const std::string& source_validation_message,
                          const std::string& moment_layout_message, const bool wrote_moments, const bool wrote_quadrature)
{
    std::cout << "2D CQMOM adaptive node counts\n";
    std::cout << "-----------------------------\n";
    std::cout << std::left << std::setw(22) << "Phase space" << "(d, u)\n";
    std::cout << std::setw(22) << "Conditioning order" << format_index_vector(conditioning_order) << '\n';
    std::cout << std::setw(22) << "Node counts" << format_index_vector(node_counts) << '\n';
    std::cout << std::setw(22) << "Relaxation rate" << std::scientific << std::setprecision(6) << relaxation_rate << '\n';
    std::cout << std::left << std::setw(22) << "Equilibrium velocity" << std::scientific << std::setprecision(6) << equilibrium_velocity << '\n';
    std::cout << std::left << std::setw(22) << "Pulse interval" << '[' << std::scientific << std::setprecision(6) << pulse_start << ", " << pulse_end
              << ")\n";
    std::cout << std::left << std::setw(22) << "Moment-source rate" << std::scientific << std::setprecision(6) << source_rate << '\n';
    std::cout << std::left << std::setw(22) << "Velocity displacement" << std::scientific << std::setprecision(6) << velocity_displacement << '\n';
    std::cout << std::left << std::setw(22) << "Time step" << std::scientific << std::setprecision(6) << time_step << '\n';
    std::cout << std::left << std::setw(22) << "Final time" << std::scientific << std::setprecision(6) << final_time << '\n';

    print_selected_moments("Initial moments", initial_moments);

    print_selected_moments("Final moments", final_moments);

    print_variance_history(output_history);

    print_active_count_transitions(active_count_transitions);

    print_three_node_growth_diagnostics(active_count_transitions);

    print_quadrature(final_quadrature);

    print_section_heading("Validation");
    std::cout << "[pass] Mass conservation\n";
    std::cout << "[pass] Velocity-mean relaxation\n";
    std::cout << "[pass] Variance injection\n";
    std::cout << "[pass] Adaptive collapse and regrowth\n";
    std::cout << "[pass] Source moments through order " << source_rate_checked_order << "   " << source_validation_message << '\n';
    std::cout << "[pass] Required velocity moment layout  " << moment_layout_message << '\n';
    if (required_velocity_order > source_rate_checked_order)
    {
        std::cout << "[note] Source self-check covers order " << source_rate_checked_order << "; required order is " << required_velocity_order
                  << '\n';
    }

    if (wrote_moments || wrote_quadrature)
    {
        print_section_heading("Output files");
        std::cout << std::left << std::setw(10) << "directory" << output_directory.string() << '\n';
        std::cout << std::setw(10) << "files";
        bool first_file = true;
        if (wrote_moments)
        {
            std::cout << "initial_moments.txt, moment_history.txt";
            first_file = false;
        }
        if (wrote_quadrature)
        {
            if (!first_file)
            {
                std::cout << ", ";
            }
            std::cout << "initial_quadrature.txt, quadrature_history.txt";
        }
        std::cout << '\n';
    }
    else
    {
        print_section_heading("Output files");
        std::cout << "disabled by configuration\n";
    }
}

} // namespace cqmom_2d_adaptive_node_counts
