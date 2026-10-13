#include "monte_carlo/evolution.hpp"
#include "monte_carlo/moments.hpp"
#include "monte_carlo/output.hpp"
#include "particles_mc.hpp"
#include "physics.hpp"
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
#include <stdexcept>
#include <utility>
#include <vector>

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";
    const double initial_time = 0.0;
    const double final_time = 13.9;
    const double time_step = 0.01;

    // The hierarchy decides which moments the sample estimates.
    auto configured_case = qmomkit::config::load_cqmom_prepared(configuration_file);
    if (!configured_case)
    {
        std::cerr << configured_case.error() << '\n';
        return 1;
    }
    auto configured_output = qmomkit::config::load_output_options(configuration_file);
    if (!configured_output)
    {
        std::cerr << configured_output.error() << '\n';
        return 1;
    }
    const qmomkit::OutputOptions output = std::move(configured_output).value();

    // The same sample the QMOMKit run starts from, here kept as bubbles of equal
    // statistical weight.
    const std::optional<qmomkit::monte_carlo::TensorSampledInitialization> sampled =
        cqmom_2d_bubbles_cavitation::sample_initial_distribution(configuration_file, configured_case.value().hierarchy());
    if (!sampled)
    {
        return 1;
    }
    qmomkit::monte_carlo::Ensemble ensemble = sampled->sampling.ensemble;
    const double particle_weight = sampled->sampling.number_density / static_cast<double>(ensemble.size());
    const std::vector<std::size_t> moment_extents = sampled->moments.moments.extents;

    const qmomkit::PhaseSpaceVelocity velocity = cqmom_2d_bubbles_cavitation::advection_velocity();
    const qmomkit::AdaptiveRkOptions advection_options;
    const qmomkit::RawMomentRequest printed_moments = cqmom_2d_bubbles_cavitation::printed_moments();
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
                qmomkit::monte_carlo::make_moment_history_writer(output.directory / "moment_history_mc.txt", moment_extents, output));
        }
        const auto write_output = [&](const double time)
        {
            if (moment_writer)
            {
                qmomkit::monte_carlo::write_moment_history_row(
                    *moment_writer, time, cqmom_2d_bubbles_cavitation::compute_particle_moments(ensemble, particle_weight, moment_extents), ensemble,
                    particle_weight * static_cast<double>(ensemble.size()), output);
            }
        };
        write_output(initial_time);

        const double output_interval = output.interval.value_or(final_time - initial_time);
        double next_output_time = initial_time + output_interval;
        double time = initial_time;
        while (time < final_time - 1.0e-12)
        {
            // Each bubble follows its own radius and wall-velocity trajectory.
            const double step = std::min(time_step, final_time - time);
            qmomkit::monte_carlo::EnsembleEvolutionResult advection =
                qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, time, time + step, advection_options);
            if (!advection.ok())
            {
                std::cerr << "Bubble advection failed: " << advection.message << '\n';
                return 1;
            }
            ensemble = std::move(advection.ensemble);
            time += step;

            if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
            {
                write_output(time);
                next_output_time += output_interval;
            }
        }

        cqmom_2d_bubbles_cavitation::print_moments(std::cout, "Monte Carlo run", initial_moments, estimate_printed_moments());
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
