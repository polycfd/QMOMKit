#include "utilities.hpp"

#include "qmomkit/qmom.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <system_error>

namespace cqmom_2d_variance_reinjection
{
namespace
{

constexpr std::size_t max_console_moments = 12;

void print_multi_index(const qmomkit::MultiIndex& order)
{
    std::cout << '[';
    for (std::size_t dimension = 0; dimension < order.size(); ++dimension)
    {
        if (dimension > 0)
        {
            std::cout << ',';
        }
        std::cout << order[dimension];
    }
    std::cout << ']';
}

void print_path(const qmomkit::ConditioningPath& path)
{
    std::cout << '[';
    for (std::size_t level = 0; level < path.size(); ++level)
    {
        if (level > 0)
        {
            std::cout << ',';
        }
        std::cout << path.node_index_at_level(level);
    }
    std::cout << ']';
}

void print_paths(const std::vector<qmomkit::ConditioningPath>& paths)
{
    for (std::size_t path = 0; path < paths.size(); ++path)
    {
        if (path > 0)
        {
            std::cout << ',';
        }
        print_path(paths[path]);
    }
}

void print_conditioning_order(const std::vector<std::size_t>& conditioning_order)
{
    std::cout << '[';
    for (std::size_t dimension = 0; dimension < conditioning_order.size(); ++dimension)
    {
        if (dimension > 0)
        {
            std::cout << ',';
        }
        std::cout << conditioning_order[dimension];
    }
    std::cout << ']';
}

void print_selected_moments(const char* label, const qmomkit::SelectedRawMomentResult& moments)
{
    std::cout << label << '\n';
    if (!moments.ok())
    {
        std::cout << "  unavailable: " << moments.message << '\n';
        return;
    }

    const std::size_t printed_count = std::min(max_console_moments, moments.values.size());
    for (std::size_t index = 0; index < printed_count; ++index)
    {
        std::cout << "M";
        print_multi_index(moments.indices[index]);
        std::cout << " = " << moments.values[index] << '\n';
    }
    if (moments.values.size() > printed_count)
    {
        std::cout << "... " << (moments.values.size() - printed_count) << " more moments omitted from the console summary\n";
    }
}

void print_quadrature(const qmomkit::CqmomQuadrature& quadrature)
{
    std::cout << "node  path       weight          d               u\n";
    for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
    {
        std::cout << node << "     ";
        print_path(quadrature.terminal_paths()[node]);
        std::cout << "   " << quadrature.weights()[node] << "   " << quadrature.abscissas()[0][node] << "   " << quadrature.abscissas()[1][node]
                  << '\n';
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

bool same_branch(const VelocityBranchId& lhs, const VelocityBranchId& rhs)
{
    return lhs.conditioning_level == rhs.conditioning_level && lhs.phase_space_dimension == rhs.phase_space_dimension &&
           lhs.parent_path == rhs.parent_path;
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

std::optional<VelocityBranchCount> find_branch_count(const std::vector<VelocityBranchCount>& counts, const VelocityBranchId& branch)
{
    const auto found =
        std::find_if(counts.begin(), counts.end(), [&branch](const VelocityBranchCount& count) { return same_branch(count.id, branch); });
    if (found == counts.end())
    {
        return std::nullopt;
    }
    return *found;
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
        const qmomkit::ConditioningPath& path = quadrature.terminal_paths()[node];
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

void print_branch_id(const VelocityBranchId& branch)
{
    std::cout << "level " << branch.conditioning_level << " dim " << branch.phase_space_dimension << " parent ";
    print_path(branch.parent_path);
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
    case qmomkit::QmomStatus::no_adaptive_quadrature:
        return "no_adaptive_quadrature";
    }
    return "unknown";
}

void print_branch_counts(const std::vector<VelocityBranchCount>& counts)
{
    std::cout << '[';
    for (std::size_t branch = 0; branch < counts.size(); ++branch)
    {
        if (branch > 0)
        {
            std::cout << ',';
        }
        std::cout << "level " << counts[branch].id.conditioning_level << " dim " << counts[branch].id.phase_space_dimension << " parent ";
        print_path(counts[branch].id.parent_path);
        std::cout << ':' << counts[branch].active_count << '/' << counts[branch].configured_maximum;
    }
    std::cout << ']';
}

void print_variance_history(const std::vector<HistorySample>& history)
{
    std::cout << "\nVariance and active-count history:\n";
    std::cout << "time  source  velocity_mean  velocity_variance  active_u_branch_counts  active_paths\n";
    for (const HistorySample& sample : history)
    {
        std::cout << sample.time << "  " << (sample.source_active ? "on " : "off") << "  " << sample.moments.velocity_mean << "  "
                  << sample.moments.velocity_variance << "  ";
        print_branch_counts(sample.velocity_branch_counts);
        std::cout << "  ";
        print_paths(sample.quadrature.terminal_paths());
        std::cout << '\n';
    }
}

void print_active_count_transitions(const std::vector<ActiveCountTransition>& transitions)
{
    std::cout << "\nObserved active velocity-count transitions:\n";
    if (transitions.empty())
    {
        std::cout << "none\n";
        return;
    }

    for (const ActiveCountTransition& transition : transitions)
    {
        std::cout << transition.time << "  ";
        print_branch_id(transition.branch);
        std::cout << "  " << transition.before << " -> " << transition.after << " / " << transition.configured_maximum;
        if (!transition.velocity_moments.empty())
        {
            std::cout << "  R2 = " << static_cast<double>(transition.rank_indicator);
        }
        std::cout << '\n';
        if (transition.before == 2 && transition.after == 3)
        {
            std::cout << "    velocity moments:";
            for (std::size_t order = 0; order < transition.velocity_moments.size(); ++order)
            {
                std::cout << " m" << order << '=' << transition.velocity_moments[order];
            }
            std::cout << '\n';
            qmomkit::QmomOptions options;
            options.adaptive = true;
            const qmomkit::QmomResult current_ceiling =
                qmomkit::compute_qmom(qmomkit::MomentSet{transition.velocity_moments}, transition.before, options);
            const qmomkit::QmomResult growth_ceiling =
                qmomkit::compute_qmom(qmomkit::MomentSet{transition.velocity_moments}, transition.after, options);
            std::cout << "    Wheeler diagnostic: ceiling " << transition.before << " returned " << current_ceiling.abscissas.size() << " status "
                      << qmom_status_name(current_ceiling.status) << "; ceiling " << transition.after << " returned "
                      << growth_ceiling.abscissas.size() << " status " << qmom_status_name(growth_ceiling.status) << '\n';
        }
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

void append_active_count_transitions(const HistorySample& previous, const HistorySample& current, std::vector<ActiveCountTransition>& transitions,
                                     const std::size_t velocity_dimension)
{
    for (const VelocityBranchCount& after : current.velocity_branch_counts)
    {
        const std::optional<VelocityBranchCount> before = find_branch_count(previous.velocity_branch_counts, after.id);
        if (!before.has_value() || before->active_count == after.active_count)
        {
            continue;
        }

        const std::size_t diagnostic_order = std::max(std::size_t{4}, 2 * after.configured_maximum - 1);
        const std::vector<double> moments = branch_velocity_moments(current.quadrature, after.id, diagnostic_order, velocity_dimension);
        transitions.push_back(ActiveCountTransition{current.time, after.id, after.configured_maximum, before->active_count, after.active_count,
                                                    moments, centered_rank_indicator(moments)});
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
        message = "Velocity mean deviated from the continuous relaxation reference by more than 5e-5.";
        return false;
    }
    if (!std::isfinite(variance_before_pulse) || maximum_pulse_variance <= variance_before_pulse + 1.0e-4)
    {
        message = "Velocity variance did not increase clearly during the source pulse.";
        return false;
    }

    message = "mass conserved, velocity mean follows continuous relaxation, and pulse increases velocity variance";
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
    stream << "source layout exposes velocity moments through order " << required_velocity_order;
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
        if (transition.after > transition.before)
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
            if (transition.after > std::min(transition.before + 1, transition.configured_maximum))
            {
                message = "A velocity branch grew by more than one active node.";
                return false;
            }
            saw_pulse_growth = true;
            if (transition.before == 2 && transition.after == 3)
            {
                saw_two_to_three_growth = true;
            }
        }
        else if (transition.after < transition.before)
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
                          const std::vector<std::size_t>& conditioning_order, const double relaxation_rate, const double equilibrium_velocity,
                          const double pulse_start, const double pulse_end, const double source_rate, const double velocity_displacement,
                          const double time_step, const double final_time, const bool wrote_moments, const bool wrote_quadrature)
{
    std::cout << std::setprecision(12);
    std::cout << "2D CQMOM variance reinjection example\n";
    std::cout << "phase-space ordering: (d,u)\n";
    std::cout << "conditioning_order = ";
    print_conditioning_order(conditioning_order);
    std::cout << '\n';
    std::cout << "relaxation_rate = " << relaxation_rate << '\n';
    std::cout << "equilibrium_velocity = " << equilibrium_velocity << '\n';
    std::cout << "pulse_start = " << pulse_start << '\n';
    std::cout << "pulse_end = " << pulse_end << '\n';
    std::cout << "source_rate = " << source_rate << '\n';
    std::cout << "velocity_displacement = " << velocity_displacement << '\n';
    std::cout << "time_step = " << time_step << '\n';
    std::cout << "final_time = " << final_time << '\n';

    std::cout << '\n';
    print_selected_moments("Initial moments:", initial_moments);

    std::cout << '\n';
    print_selected_moments("Final moments:", final_moments);

    print_variance_history(output_history);

    print_active_count_transitions(active_count_transitions);

    std::cout << "\nFinal quadrature:\n";
    print_quadrature(final_quadrature);

    if (wrote_moments || wrote_quadrature)
    {
        std::cout << "\nOutputs written to:\n";
        if (wrote_moments)
        {
            std::cout << (output_directory / "initial_moments.txt").string() << '\n';
            std::cout << (output_directory / "moment_history.txt").string() << '\n';
        }
        if (wrote_quadrature)
        {
            std::cout << (output_directory / "initial_quadrature.txt").string() << '\n';
            std::cout << (output_directory / "quadrature_history.txt").string() << '\n';
        }
    }
    else
    {
        std::cout << "\nOutput disabled by configuration.\n";
    }
}

} // namespace cqmom_2d_variance_reinjection
