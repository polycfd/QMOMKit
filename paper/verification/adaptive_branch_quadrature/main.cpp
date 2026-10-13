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
// The source and the reduction of the velocity node count both preserve number
// density and mean velocity exactly, so their errors come only from the node
// advection, whose relative tolerance is 1e-12, and from round-off.
constexpr double number_density_tolerance = 1.0e-11;
constexpr double mean_tolerance = 1.0e-11;

struct Sample
{
    double time = 0.0;
    double number_density = 0.0;
    double mean_velocity = 0.0;
    std::size_t velocity_node_count = 0;
};

// The source sees the pulse only strictly inside (pulse_start, pulse_end), so the
// step that ends at pulse_start, whose last stage is evaluated there, injects
// nothing.
bool source_pulse_active(const double time)
{
    return time > pulse_start + 1.0e-12 && time < pulse_end - 1.0e-12;
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
// by the average of two copies of itself shifted by -velocity_displacement and
// +velocity_displacement. This injects velocity variance without changing the
// number density or the mean, so a collapsed velocity branch can regrow.
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

// The number density and mean velocity reported in the paper figure and used
// for the analytical checks below, from the quadrature of the state. The size
// coordinate has a single node, so the terminal nodes are the active nodes of
// the one velocity branch.
Sample make_sample(const double time, const qmomkit::CqmomState& state)
{
    const qmomkit::RawMomentRequest request{{{0, 0}, {0, 1}}};
    const qmomkit::SelectedRawMomentResult moments = qmomkit::reconstruct_raw_moments(state.quadrature(), request);
    if (!moments.ok())
    {
        throw std::runtime_error("The velocity moment reconstruction failed: " + moments.message);
    }
    return Sample{time, moments.values[0], moments.values[1] / moments.values[0], state.quadrature().total_nodes()};
}

void write_csv(const std::vector<Sample>& samples)
{
    std::ofstream output("adaptive_branch_quadrature.csv");
    if (!output)
    {
        throw std::runtime_error("Could not open adaptive_branch_quadrature.csv for writing.");
    }
    output << "time,number_density,mean_velocity,active_velocity_nodes\n";
    output << std::scientific << std::setprecision(17);
    for (const Sample& sample : samples)
    {
        output << sample.time << ',' << sample.number_density << ',' << sample.mean_velocity << ',' << sample.velocity_node_count << '\n';
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
        const qmomkit::CqmomPhaseSpaceEvolutionOptions options = std::move(configured_options).value();
        const qmomkit::OutputOptions output = std::move(configured_output).value();

        qmomkit::CqmomState state = cqmom_case.initial_state();
        const Sample initial_sample = make_sample(initial_time, state);
        std::vector<Sample> output_samples{initial_sample};
        std::size_t minimum_node_count = initial_sample.velocity_node_count;
        std::size_t maximum_node_count_after_injection = 0;
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

            const Sample current = make_sample(time, state);
            minimum_node_count = std::min(minimum_node_count, current.velocity_node_count);
            if (time >= pulse_start - 1.0e-12)
            {
                maximum_node_count_after_injection = std::max(maximum_node_count_after_injection, current.velocity_node_count);
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

        // The velocity branch must lose nodes as relaxation collapses it, and
        // regrow after the variance injection; number density and mean velocity
        // must follow their exact values.
        const bool collapsed = minimum_node_count < initial_sample.velocity_node_count;
        const bool regrew = maximum_node_count_after_injection > minimum_node_count;
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
