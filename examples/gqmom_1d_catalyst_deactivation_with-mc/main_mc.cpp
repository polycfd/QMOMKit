#include "monte_carlo/evolution.hpp"
#include "monte_carlo/moments.hpp"
#include "monte_carlo/output.hpp"
#include "physics.hpp"
#include "qmomkit/configuration.hpp"
#include "qmomkit/output.hpp"
#include "sampling.hpp"

#include <cmath>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <optional>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";
    const int step_count = 100;
    const double time_step = 0.1;

    auto configured_output = qmomkit::config::load_output_options(configuration_file);
    if (!configured_output)
    {
        std::cerr << configured_output.error() << '\n';
        return 1;
    }
    const qmomkit::OutputOptions output = std::move(configured_output).value();

    // The same sample the QMOMKit runs start from, here kept as particles of
    // equal statistical weight.
    std::optional<gqmom_1d_catalyst_deactivation_with_mc::SamplingSetup> sampled =
        gqmom_1d_catalyst_deactivation_with_mc::sample_initial_distribution(configuration_file);
    if (!sampled)
    {
        return 1;
    }
    qmomkit::monte_carlo::Ensemble particles = std::move(sampled->sampling.ensemble);
    const std::size_t initial_particle_count = particles.size();
    const double particle_weight = sampled->sampling.number_density / static_cast<double>(initial_particle_count);
    const std::size_t maximum_order = sampled->moments.moments.values().size() - 1;
    const std::vector<double> initial_moments = sampled->moments.moments.values();
    const qmomkit::PhaseSpaceVelocity velocity = gqmom_1d_catalyst_deactivation_with_mc::coking_velocity();
    std::mt19937_64 removal_generator(gqmom_1d_catalyst_deactivation_with_mc::removal_seed);

    try
    {
        // The writer throws when a file cannot be written, and the coking steps
        // and current_moments with the message of a failure the companion
        // reports.
        std::optional<qmomkit::MomentWriter> moment_writer;
        if (output.write_moments)
        {
            std::filesystem::create_directories(output.directory);
            moment_writer.emplace(
                qmomkit::monte_carlo::make_moment_history_writer(output.directory / "moment_history_mc.txt", initial_moments.size(), output));
        }
        const auto current_moments = [&]()
        {
            qmomkit::monte_carlo::OneDimensionalMomentEstimate estimate =
                qmomkit::monte_carlo::compute_ensemble_moments(particles, particle_weight * static_cast<double>(particles.size()), maximum_order);
            if (!estimate.ok())
            {
                throw std::runtime_error(estimate.message);
            }
            return estimate.moments;
        };
        const auto write_output = [&](const double time)
        {
            if (moment_writer)
            {
                qmomkit::monte_carlo::write_moment_history_row(*moment_writer, time, current_moments(), particles,
                                                               particle_weight * static_cast<double>(particles.size()), output);
            }
        };
        write_output(0.0);

        // Write every output interval of [output], and at the end without one.
        const int steps_per_output = output.interval ? static_cast<int>(std::lround(*output.interval / time_step)) : step_count;
        for (int step = 1; step <= step_count; ++step)
        {
            // The Strang-split step of the QMOMKit runs, applied to the particles:
            // coking over half the step, stochastic removal over the step, and
            // coking over the other half.
            const double time = (step - 1) * time_step;
            qmomkit::monte_carlo::EnsembleEvolutionResult first_half =
                qmomkit::monte_carlo::evolve_ensemble(particles, velocity, time, time + (0.5 * time_step));
            if (!first_half.ok())
            {
                throw std::runtime_error(first_half.message);
            }
            particles = std::move(first_half.ensemble);
            gqmom_1d_catalyst_deactivation_with_mc::apply_removal(particles, time_step, removal_generator);
            qmomkit::monte_carlo::EnsembleEvolutionResult second_half =
                qmomkit::monte_carlo::evolve_ensemble(particles, velocity, time + (0.5 * time_step), time + time_step);
            if (!second_half.ok())
            {
                throw std::runtime_error(second_half.message);
            }
            particles = std::move(second_half.ensemble);

            if (step % steps_per_output == 0 || step == step_count)
            {
                write_output(step * time_step);
            }
        }

        // The initial and final raw moments side by side, as main.cpp prints them.
        const std::vector<double> final_moments = current_moments().values();
        std::cout << "Monte Carlo run, raw moments of the coked-fraction distribution:\n";
        std::cout << "      " << std::setw(14) << "initial" << std::setw(14) << "final" << '\n';
        for (std::size_t order = 0; order < initial_moments.size(); ++order)
        {
            std::cout << "  M" << order << "  " << std::setw(14) << initial_moments[order] << std::setw(14) << final_moments[order] << '\n';
        }
        std::cout << "Particles: " << initial_particle_count << " -> " << particles.size() << '\n';
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
