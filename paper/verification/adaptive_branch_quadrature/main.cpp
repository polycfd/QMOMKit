#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace
{

constexpr double initial_time = 0.0;
constexpr double final_time = 8.0;
constexpr double time_step = 0.002;
constexpr std::size_t velocity_dimension = 1;
constexpr double relaxation_rate = 5.0;
constexpr double equilibrium_velocity = 0.5;
constexpr double pulse_start = 4.0;
constexpr double pulse_end = 4.05;
constexpr double source_rate = 20.0;
constexpr double velocity_displacement = 0.10;
constexpr double number_density_tolerance = 1.0e-8;
constexpr double mean_tolerance = 5.0e-5;

struct Sample
{
    double time = 0.0;
    double number_density = 0.0;
    double mean_velocity = 0.0;
    double velocity_variance = 0.0;
    std::vector<std::size_t> branch_counts;
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
    // Velocities relax toward the equilibrium velocity; sizes do not move.
    return qmomkit::make_phase_space_velocity(
        [](double, const std::vector<double>& coordinates)
        { return std::vector<double>{0.0, -relaxation_rate * (coordinates[velocity_dimension] - equilibrium_velocity)}; });
}

// During the pulse, each velocity distribution is replaced at rate source_rate
// by two copies of itself shifted by -velocity_displacement and
// +velocity_displacement. This injects velocity variance without changing the
// mean, so a collapsed velocity branch can regrow.
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
            const qmomkit::MultiIndexView index = context.index(entry);
            double redistribution_change = 0.0;
            for (std::size_t displacement_order = 2; displacement_order <= index[velocity_dimension]; displacement_order += 2)
            {
                qmomkit::MultiIndex lower = index.to_multi_index();
                lower[velocity_dimension] -= displacement_order;
                redistribution_change += binomial_coefficient(index[velocity_dimension], displacement_order) *
                                         std::pow(velocity_displacement, static_cast<double>(displacement_order)) * context.moment(lower);
            }
            derivative[entry] = source_rate * redistribution_change;
        }
        return qmomkit::CqmomMomentRate{std::move(derivative)};
    };
}

// Active node count of the velocity branch conditioned on each size node.
std::vector<std::size_t> velocity_branch_counts(const std::vector<qmomkit::CqmomActiveBranchCount>& active_branches,
                                                const std::size_t configured_size_node_count)
{
    std::vector<std::size_t> counts(configured_size_node_count, 0);
    for (const qmomkit::CqmomActiveBranchCount& branch : active_branches)
    {
        if (branch.phase_space_dimension == velocity_dimension)
        {
            counts[branch.parent_path.node_index_at_level(0)] = branch.active_node_count;
        }
    }
    return counts;
}

// The number density and velocity statistics reported in the paper figure and
// used for the analytical checks below, from the quadrature of the state.
Sample make_sample(const double time, const qmomkit::CqmomState& state, const std::size_t configured_size_node_count)
{
    const qmomkit::RawMomentRequest request{{{0, 0}, {0, 1}, {0, 2}}};
    const qmomkit::SelectedRawMomentResult moments = qmomkit::reconstruct_raw_moments(state.quadrature(), request);
    if (!moments.ok())
    {
        throw std::runtime_error("The velocity moment reconstruction failed: " + moments.message);
    }
    const double mean = moments.values[1] / moments.values[0];
    const double variance = moments.values[2] / moments.values[0] - mean * mean;
    return Sample{time, moments.values[0], mean, variance, velocity_branch_counts(state.active_branch_counts(), configured_size_node_count)};
}

void write_csv(const std::vector<Sample>& samples)
{
    std::ofstream output("adaptive_branch_quadrature.csv");
    if (!output)
    {
        throw std::runtime_error("Could not open adaptive_branch_quadrature.csv for writing.");
    }
    output << "time,total_mass,mean_velocity,velocity_variance,branch_index,active_nodes\n";
    output << std::scientific << std::setprecision(17);
    for (const Sample& sample : samples)
    {
        for (std::size_t branch = 0; branch < sample.branch_counts.size(); ++branch)
        {
            output << sample.time << ',' << sample.number_density << ',' << sample.mean_velocity << ',' << sample.velocity_variance << ',' << branch
                   << ',' << sample.branch_counts[branch] << '\n';
        }
    }
}

} // namespace

int main()
{
    try
    {
        const std::filesystem::path configuration_file = "configuration.qmom";
        auto configured_case = qmomkit::config::load_cqmom_prepared(configuration_file);
        if (!configured_case)
        {
            std::cerr << configured_case.error() << '\n';
            return 1;
        }
        auto configured_options = qmomkit::config::load_cqmom_phase_space_evolution_options(configuration_file);
        if (!configured_options)
        {
            std::cerr << configured_options.error() << '\n';
            return 1;
        }
        auto configured_output = qmomkit::config::load_output_options(configuration_file);
        if (!configured_output)
        {
            std::cerr << configured_output.error() << '\n';
            return 1;
        }
        const qmomkit::CqmomPreparedCase cqmom_case = std::move(configured_case).value();
        const std::size_t configured_size_node_count = cqmom_case.hierarchy().node_count(0);
        const qmomkit::CqmomPhaseSpaceEvolutionOptions options = std::move(configured_options).value();
        const qmomkit::OutputOptions output = std::move(configured_output).value();

        qmomkit::CqmomState state = cqmom_case.initial_state();
        const Sample initial_sample = make_sample(initial_time, state, configured_size_node_count);
        std::vector<Sample> output_samples{initial_sample};
        const std::vector<std::size_t> initial_counts = initial_sample.branch_counts;
        std::vector<std::size_t> minimum_counts = initial_counts;
        std::vector<std::size_t> maximum_after_injection(configured_size_node_count, 0);
        double maximum_number_density_error = 0.0;
        double maximum_mean_error = 0.0;

        const qmomkit::PhaseSpaceVelocity advection_velocity = velocity();
        const qmomkit::CqmomMomentSource injection_source = variance_injection_source();
        const double output_interval = output.interval.value_or(final_time - initial_time);
        double time = initial_time;
        double next_output_time = initial_time + output_interval;
        while (time < final_time - 1.0e-12)
        {
            const double step = std::min(time_step, final_time - time);
            qmomkit::CqmomPhaseSpaceEvolutionResult result =
                qmomkit::advance_cqmom_phase_space(state, time, step, advection_velocity, injection_source, options);
            if (!result.ok())
            {
                std::cerr << "CQMOM phase-space evolution failed: " << result.message << '\n';
                return 1;
            }
            state = std::move(result.state);
            time = result.final_time;

            const Sample current = make_sample(time, state, configured_size_node_count);
            for (std::size_t branch = 0; branch < current.branch_counts.size(); ++branch)
            {
                minimum_counts[branch] = std::min(minimum_counts[branch], current.branch_counts[branch]);
                if (time >= pulse_start - 1.0e-12)
                {
                    maximum_after_injection[branch] = std::max(maximum_after_injection[branch], current.branch_counts[branch]);
                }
            }
            maximum_number_density_error = std::max(maximum_number_density_error, std::abs(current.number_density - initial_sample.number_density));
            const double exact_mean =
                equilibrium_velocity + (initial_sample.mean_velocity - equilibrium_velocity) * std::exp(-relaxation_rate * current.time);
            maximum_mean_error = std::max(maximum_mean_error, std::abs(current.mean_velocity - exact_mean));

            if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
            {
                output_samples.push_back(current);
                while (next_output_time <= time + 1.0e-12)
                {
                    next_output_time += output_interval;
                }
            }
        }

        write_csv(output_samples);

        // Every velocity branch must lose nodes as relaxation collapses it, and
        // regrow after the variance injection; number density and mean velocity
        // must follow their exact values.
        bool collapsed = true;
        bool regrew = true;
        for (std::size_t branch = 0; branch < configured_size_node_count; ++branch)
        {
            collapsed = collapsed && minimum_counts[branch] < initial_counts[branch];
            regrew = regrew && maximum_after_injection[branch] > minimum_counts[branch];
        }
        if (!collapsed || !regrew || maximum_number_density_error > number_density_tolerance || maximum_mean_error > mean_tolerance)
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
