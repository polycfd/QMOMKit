#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

constexpr double initial_time = 0.0;
constexpr double final_time = 8.0;
constexpr double time_step = 0.002;
constexpr double output_interval = 0.01;
constexpr std::size_t velocity_dimension = 1;
constexpr double relaxation_rate = 5.0;
constexpr double equilibrium_velocity = 0.0;
constexpr double pulse_start = 4.0;
constexpr double pulse_end = 4.05;
constexpr double source_rate = 20.0;
constexpr double velocity_displacement = 0.10;
constexpr double mass_tolerance = 1.0e-8;
constexpr double mean_tolerance = 5.0e-5;

struct Sample
{
    double time = 0.0;
    double mass = 0.0;
    double mean_velocity = 0.0;
    double velocity_variance = 0.0;
    std::array<std::size_t, 2> branch_counts{};
};

bool source_pulse_active(const double time)
{
    return time + 1.0e-12 >= pulse_start && time < pulse_end - 1.0e-12;
}

double binomial_coefficient(const std::size_t n, const std::size_t k)
{
    const std::size_t reduced = std::min(k, n - k);
    double coefficient = 1.0;
    for (std::size_t factor = 1; factor <= reduced; ++factor)
    {
        coefficient *= static_cast<double>(n + 1 - factor);
        coefficient /= static_cast<double>(factor);
    }
    return coefficient;
}

qmomkit::PhaseSpaceVelocity velocity()
{
    return [](double, const std::vector<double>& coordinates)
    {
        if (coordinates.size() != 2)
        {
            throw std::invalid_argument("The adaptive paper study requires phase-space coordinates ordered as (d,u).");
        }
        return std::vector<double>{0.0, -relaxation_rate * (coordinates[velocity_dimension] - equilibrium_velocity)};
    };
}

qmomkit::CqmomMomentSource variance_injection_source()
{
    return [](const double time, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> derivative(context.size(), 0.0);
        if (!source_pulse_active(time))
        {
            return qmomkit::CqmomMomentRate{std::move(derivative)};
        }
        for (std::size_t entry = 0; entry < context.size(); ++entry)
        {
            const qmomkit::MultiIndex& index = context.index(entry);
            if (index.size() != 2)
            {
                throw std::invalid_argument("The adaptive paper source requires 2D raw moments ordered as (d,u).");
            }
            double redistribution_change = 0.0;
            for (std::size_t displacement_order = 2; displacement_order <= index[velocity_dimension]; displacement_order += 2)
            {
                qmomkit::MultiIndex lower = index;
                lower[velocity_dimension] -= displacement_order;
                const std::optional<double> lower_moment = context.value(lower);
                if (!lower_moment.has_value())
                {
                    throw std::invalid_argument("The variance-injection source requires a lower-order velocity moment.");
                }
                redistribution_change += binomial_coefficient(index[velocity_dimension], displacement_order) *
                                         std::pow(velocity_displacement, static_cast<double>(displacement_order)) * *lower_moment;
            }
            derivative[entry] = source_rate * redistribution_change;
        }
        return qmomkit::CqmomMomentRate{std::move(derivative)};
    };
}

bool finite_quadrature(const qmomkit::CqmomQuadrature& quadrature)
{
    const auto finite_vector = [](const std::vector<double>& values)
    { return std::all_of(values.begin(), values.end(), [](const double value) { return std::isfinite(value); }); };
    if (!finite_vector(quadrature.weights()))
    {
        return false;
    }
    return std::all_of(quadrature.abscissas().begin(), quadrature.abscissas().end(), finite_vector);
}

std::array<std::size_t, 2> velocity_branch_counts(const std::vector<qmomkit::CqmomActiveBranchCount>& active_branches)
{
    std::array<std::size_t, 2> counts{};
    std::array<bool, 2> found{};
    for (const qmomkit::CqmomActiveBranchCount& branch : active_branches)
    {
        if (branch.conditioning_level != 1 || branch.phase_space_dimension != velocity_dimension || branch.parent_path.size() != 1)
        {
            continue;
        }
        const std::size_t branch_index = branch.parent_path.node_index_at_level(0);
        if (branch_index >= counts.size())
        {
            throw std::runtime_error("The adaptive paper study expected exactly two first-level size branches.");
        }
        counts[branch_index] = branch.active_node_count;
        found[branch_index] = true;
    }
    if (!found[0] || !found[1])
    {
        throw std::runtime_error("The adaptive paper study could not find both velocity-conditioned size branches.");
    }
    return counts;
}

Sample make_sample(const double time, const qmomkit::CqmomQuadrature& quadrature, const std::vector<qmomkit::CqmomActiveBranchCount>& active_branches)
{
    const qmomkit::RawMomentRequest request{{{0, 0}, {0, 1}, {0, 2}}};
    const qmomkit::SelectedRawMomentResult moments = qmomkit::reconstruct_raw_moments(quadrature, request);
    if (!moments.ok() || moments.values.size() != 3 ||
        !std::all_of(moments.values.begin(), moments.values.end(), [](const double value) { return std::isfinite(value); }))
    {
        throw std::runtime_error("The adaptive paper study could not reconstruct finite diagnostic moments.");
    }
    if (moments.values[0] <= 0.0)
    {
        throw std::runtime_error("The adaptive paper study produced a non-positive total mass.");
    }
    const double mean = moments.values[1] / moments.values[0];
    const double variance = moments.values[2] / moments.values[0] - mean * mean;
    if (!std::isfinite(mean) || !std::isfinite(variance))
    {
        throw std::runtime_error("The adaptive paper study produced non-finite velocity statistics.");
    }
    return Sample{time, moments.values[0], mean, variance, velocity_branch_counts(active_branches)};
}

bool transitions_match_counts(const std::array<std::size_t, 2>& before, const std::array<std::size_t, 2>& after,
                              const std::vector<qmomkit::CqmomNodeCountTransition>& transitions)
{
    for (const qmomkit::CqmomNodeCountTransition& transition : transitions)
    {
        if (transition.conditioning_level != 1 || transition.phase_space_dimension != velocity_dimension || transition.parent_path.size() != 1)
        {
            continue;
        }
        const std::size_t branch = transition.parent_path.node_index_at_level(0);
        if (branch >= before.size() || !transition.active_node_count_before.has_value() || !transition.active_node_count_after.has_value() ||
            *transition.active_node_count_before != before[branch] || *transition.active_node_count_after != after[branch])
        {
            return false;
        }
    }
    return true;
}

void write_csv(const std::vector<Sample>& samples)
{
    std::ofstream output("adaptive_branch_quadrature.csv");
    if (!output)
    {
        throw std::runtime_error("Could not open adaptive_branch_quadrature.csv for writing.");
    }
    output << "time,total_mass,mean_velocity,velocity_variance,branch_0_active_nodes,branch_1_active_nodes\n";
    output << std::scientific << std::setprecision(17);
    for (const Sample& sample : samples)
    {
        output << sample.time << ',' << sample.mass << ',' << sample.mean_velocity << ',' << sample.velocity_variance << ','
               << sample.branch_counts[0] << ',' << sample.branch_counts[1] << '\n';
    }
}

void print_summary(const std::array<std::size_t, 2>& initial_counts, const std::array<std::size_t, 2>& minimum_counts,
                   const std::array<std::size_t, 2>& maximum_after_reinjection, const std::array<double, 2>& collapse_times,
                   const double maximum_mass_error, const double maximum_mean_error)
{
    std::cout << "Adaptive branch-specific quadrature\n";
    std::cout << "Stable branch ordering: branch 0 is first-level size path [0], branch 1 is path [1].\n";
    std::cout << "initial active counts: " << initial_counts[0] << ", " << initial_counts[1] << '\n';
    std::cout << "minimum active counts: " << minimum_counts[0] << ", " << minimum_counts[1] << '\n';
    std::cout << "maximum counts after reinjection: " << maximum_after_reinjection[0] << ", " << maximum_after_reinjection[1] << '\n';
    std::cout << "collapse time range: [" << collapse_times[0] << ", " << collapse_times[1] << "]\n";
    std::cout << "reinjection time range: [" << pulse_start << ", " << pulse_end << ")\n";
    std::cout << std::scientific << std::setprecision(6) << "maximum mass error: " << maximum_mass_error << '\n';
    std::cout << "maximum mean-velocity error: " << maximum_mean_error << '\n';
}

} // namespace

int main()
{
    try
    {
        const auto configured = qmomkit::config::load_cqmom("configuration.qmom");
        if (!configured.ok())
        {
            std::cerr << configured.error() << '\n';
            return 1;
        }
        const qmomkit::CqmomCase case_ = configured.value();
        if (case_.hierarchy().node_counts() != std::vector<std::size_t>({2, 3}) ||
            case_.hierarchy().conditioning_order() != std::vector<std::size_t>({0, 1}) || !case_.adaptive_branch_counts())
        {
            throw std::runtime_error(
                "The adaptive paper configuration must define adaptive CQMOM with node counts {2,3} and conditioning order {0,1}.");
        }

        const qmomkit::CqmomResult initial_reconstruction = qmomkit::compute_cqmom(case_);
        if (!initial_reconstruction.ok() || !finite_quadrature(initial_reconstruction.quadrature))
        {
            std::cerr << "Initial CQMOM reconstruction failed: " << initial_reconstruction.message << '\n';
            return 1;
        }

        qmomkit::CqmomState state = case_.initial_state();
        qmomkit::CqmomQuadrature quadrature = initial_reconstruction.quadrature;
        std::vector<Sample> output_samples;
        Sample previous = make_sample(initial_time, quadrature, state.active_branch_counts());
        output_samples.push_back(previous);
        const double initial_mass = previous.mass;
        const std::array<std::size_t, 2> initial_counts = previous.branch_counts;
        std::array<std::size_t, 2> minimum_counts = initial_counts;
        std::array<std::size_t, 2> maximum_after_reinjection{};
        std::array<double, 2> collapse_times{std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
        double maximum_mass_error = 0.0;
        double maximum_mean_error = 0.0;
        bool transition_consistent = true;

        qmomkit::CqmomPhaseSpaceEvolutionOptions options;
        options.node_options.scheme = qmomkit::NodeAdvectionScheme::rk23;
        options.node_options.adaptive_options.relative_tolerance = 1.0e-12;
        options.node_options.adaptive_options.absolute_tolerance = 1.0e-14;
        options.node_options.adaptive_options.initial_step = 0.01;
        options.node_options.adaptive_options.max_step = 0.01;
        options.source_options.ssp_rk3_options.estimate_error = true;
        options.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;

        double time = initial_time;
        double next_output_time = output_interval;
        while (time < final_time - 1.0e-12)
        {
            const double step = std::min(time_step, final_time - time);
            const qmomkit::CqmomPhaseSpaceEvolutionResult advanced =
                qmomkit::advance_cqmom_phase_space(state, time, step, velocity(), variance_injection_source(), options);
            if (!advanced.ok() || !advanced.quadrature.has_value() || !finite_quadrature(*advanced.quadrature) ||
                !std::all_of(advanced.state.moment_values().begin(), advanced.state.moment_values().end(),
                             [](const double value) { return std::isfinite(value); }))
            {
                std::cerr << "CQMOM phase-space evolution failed: " << advanced.message << '\n';
                return 1;
            }

            state = advanced.state;
            quadrature = *advanced.quadrature;
            time = advanced.final_time;
            const Sample current = make_sample(time, quadrature, state.active_branch_counts());
            transition_consistent = transition_consistent && transitions_match_counts(previous.branch_counts, current.branch_counts,
                                                                                      advanced.diagnostics.cqmom_node_count_transitions);
            for (std::size_t branch = 0; branch < current.branch_counts.size(); ++branch)
            {
                minimum_counts[branch] = std::min(minimum_counts[branch], current.branch_counts[branch]);
                if (current.branch_counts[branch] < 3 && !std::isfinite(collapse_times[branch]))
                {
                    collapse_times[branch] = time;
                }
                if (time >= pulse_start - 1.0e-12)
                {
                    maximum_after_reinjection[branch] = std::max(maximum_after_reinjection[branch], current.branch_counts[branch]);
                }
            }
            maximum_mass_error = std::max(maximum_mass_error, std::abs(current.mass - initial_mass));
            const double expected_mean =
                equilibrium_velocity + (output_samples.front().mean_velocity - equilibrium_velocity) * std::exp(-relaxation_rate * current.time);
            maximum_mean_error = std::max(maximum_mean_error, std::abs(current.mean_velocity - expected_mean));

            if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
            {
                output_samples.push_back(current);
                while (next_output_time <= time + 1.0e-12)
                {
                    next_output_time += output_interval;
                }
            }
            previous = current;
        }

        write_csv(output_samples);
        print_summary(initial_counts, minimum_counts, maximum_after_reinjection, collapse_times, maximum_mass_error, maximum_mean_error);

        const bool initial_topology = initial_counts[0] == 3 && initial_counts[1] == 3;
        const bool collapsed = minimum_counts[0] < 3 && minimum_counts[1] < 3 && std::isfinite(collapse_times[0]) && std::isfinite(collapse_times[1]);
        const bool regrew = maximum_after_reinjection[0] == 3 && maximum_after_reinjection[1] == 3;
        if (!initial_topology || !collapsed || !regrew || !transition_consistent || maximum_mass_error > mass_tolerance ||
            maximum_mean_error > mean_tolerance)
        {
            std::cerr << "The adaptive branch-specific lifecycle verification failed.\n";
            return 1;
        }
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
