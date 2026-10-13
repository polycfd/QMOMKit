#include "breakage_model.hpp"
#include "fluid_field.hpp"
#include "monte_carlo/evolution.hpp"
#include "monte_carlo/moments.hpp"
#include "monte_carlo/output.hpp"
#include "particle_model.hpp"
#include "particles_mc.hpp"
#include "qmomkit/configuration.hpp"
#include "qmomkit/output.hpp"
#include "sampling.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <path-to-configuration.qmom>\n";
        return 1;
    }
    const std::filesystem::path configuration_file = argv[1];
    const double initial_time = 0.0;
    const double final_time = 1.0;
    const double time_step = 0.005;

    // The hierarchy sizes the moment history; the node-advection options set the
    // accuracy of the particle trajectories.
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
    const qmomkit::AdaptiveRkOptions advection_options = configured_options.value().node_options.adaptive_options;
    const qmomkit::OutputOptions output = std::move(configured_output).value();

    // The same sample the QMOMKit run starts from, here kept as particles of
    // equal statistical weight.
    const std::optional<cqmom_5d_non_spherical_particles::SampledInitialEnsemble> sampled =
        cqmom_5d_non_spherical_particles::sample_initial_ensemble(configuration_file, configured_case.value().hierarchy());
    if (!sampled)
    {
        return 1;
    }
    qmomkit::monte_carlo::Ensemble ensemble = sampled->sampling.ensemble;
    const std::size_t initial_particle_count = ensemble.size();
    const double particle_weight = sampled->sampling.number_density / static_cast<double>(initial_particle_count);
    const cqmom_5d_non_spherical_particles::BinaryFragmentationParameters breakage_parameters =
        cqmom_5d_non_spherical_particles::binary_fragmentation_parameters();
    std::mt19937_64 breakage_generator(cqmom_5d_non_spherical_particles::breakage_seed);

    // carrier_samples is the host-side carrier flow, owned here by the driver and
    // handed to the callback below by reference. The time loop publishes a new
    // host level once per outer step at the same instants as the CQMOM driver, so
    // both runs reconstruct the same carrier flow.
    cqmom_5d_non_spherical_particles::CarrierFlowSamples carrier_samples;
    carrier_samples.current = cqmom_5d_non_spherical_particles::carrier_flow_velocity(initial_time);
    carrier_samples.sample_time = initial_time;
    const qmomkit::PhaseSpaceVelocity velocity = cqmom_5d_non_spherical_particles::particle_phase_space_velocity(
        cqmom_5d_non_spherical_particles::particle_model_parameters(), carrier_samples);

    const qmomkit::RawMomentRequest printed_moments = cqmom_5d_non_spherical_particles::printed_moments();
    const auto estimate_printed_moments = [&]()
    {
        const qmomkit::monte_carlo::SelectedEnsembleMomentEstimate estimate = qmomkit::monte_carlo::compute_requested_ensemble_moments(
            ensemble, particle_weight * static_cast<double>(ensemble.size()), printed_moments);
        if (!estimate.ok())
        {
            throw std::runtime_error(estimate.message);
        }
        return estimate.values;
    };

    try
    {
        // The estimates and the writer throw when they fail.
        const std::vector<double> initial_moments = estimate_printed_moments();
        std::optional<qmomkit::MomentWriter> moment_writer;
        if (output.write_moments)
        {
            std::filesystem::create_directories(output.directory);
            moment_writer.emplace(
                qmomkit::monte_carlo::make_moment_history_writer(output.directory / "moment_history_mc.txt", sampled->moment_extents, output));
        }
        const auto write_output = [&](const double time)
        {
            if (moment_writer)
            {
                qmomkit::monte_carlo::write_moment_history_row(
                    *moment_writer, time,
                    cqmom_5d_non_spherical_particles::compute_ensemble_moments(ensemble, particle_weight, sampled->moment_extents), ensemble,
                    particle_weight * static_cast<double>(ensemble.size()), output);
            }
        };
        write_output(initial_time);

        const double output_interval = output.interval.value_or(final_time - initial_time);
        double next_output_time = initial_time + output_interval;
        double time = initial_time;
        while (time < final_time - 1.0e-12)
        {
            const double step = std::min(time_step, final_time - time);
            const double half_step = 0.5 * step;
            if (time > initial_time)
            {
                carrier_samples.previous = carrier_samples.current;
                carrier_samples.current = cqmom_5d_non_spherical_particles::carrier_flow_velocity(time);
                carrier_samples.sample_interval = time - carrier_samples.sample_time;
                carrier_samples.sample_time = time;
            }

            // Strang splitting, as in the QMOMKit run: an advection half-step of the
            // particles, stochastic binary breakage over the whole step, and a
            // second advection half-step.
            qmomkit::monte_carlo::EnsembleEvolutionResult first_advection =
                qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, time, time + half_step, advection_options);
            if (!first_advection.ok())
            {
                std::cerr << "First advection particle half-step failed: " << first_advection.message << '\n';
                return 1;
            }
            ensemble = std::move(first_advection.ensemble);

            cqmom_5d_non_spherical_particles::apply_binary_breakage(ensemble, step, breakage_generator, breakage_parameters);

            qmomkit::monte_carlo::EnsembleEvolutionResult second_advection =
                qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, time + half_step, time + step, advection_options);
            if (!second_advection.ok())
            {
                std::cerr << "Second advection particle half-step failed: " << second_advection.message << '\n';
                return 1;
            }
            ensemble = std::move(second_advection.ensemble);
            time += step;

            if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
            {
                write_output(time);
                next_output_time += output_interval;
            }
        }

        cqmom_5d_non_spherical_particles::print_moments(std::cout, "Monte Carlo run", initial_moments, estimate_printed_moments());
        std::cout << "Particles: " << initial_particle_count << " -> " << ensemble.size() << '\n';
        if (moment_writer)
        {
            std::cout << "Output written to " << output.directory.string() << "/\n";
        }
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
