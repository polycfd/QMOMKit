#include "monte_carlo/evolution.hpp"
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
#include <iomanip>
#include <iostream>
#include <optional>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";
    const double initial_time = 0.0;
    const double final_time = 1.0;
    const double time_step = 0.01;

    const std::optional<std::string> configuration_text = qmom_1d_breakage_advection_mc_initialization::read_configuration_text(configuration_file);
    if (!configuration_text.has_value())
    {
        return 1;
    }

    auto initial_condition = qmomkit::config::parse_qmom_initial_condition(*configuration_text, configuration_file.string());
    if (!initial_condition)
    {
        qmom_1d_breakage_advection_mc_initialization::print_configuration_errors(initial_condition.error());
        return 1;
    }

    auto output_configured = qmomkit::config::load_output_options(configuration_file);
    if (!output_configured)
    {
        qmom_1d_breakage_advection_mc_initialization::print_configuration_errors(output_configured.error());
        return 1;
    }
    const qmomkit::OutputOptions output_options = std::move(output_configured).value();

    const std::optional<qmom_1d_breakage_advection_mc_initialization::SamplingSetup> sampled =
        qmom_1d_breakage_advection_mc_initialization::sample_initial_distribution(std::move(initial_condition).value());
    if (!sampled.has_value())
    {
        return 1;
    }
    if (!sampled->sampling.ok())
    {
        std::cerr << sampled->sampling.message << '\n';
        return 1;
    }
    if (!sampled->moments.ok())
    {
        std::cerr << sampled->moments.message << '\n';
        return 1;
    }

    const std::size_t maximum_order = sampled->moments.moments.values().size() - 1;

    std::vector<qmom_1d_breakage_advection_mc_initialization::Particle> particles;
    try
    {
        particles = qmom_1d_breakage_advection_mc_initialization::make_particles(sampled->sampling.ensemble);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    const double particle_weight = sampled->sampling.number_density / static_cast<double>(sampled->sampling.ensemble.size());
    std::size_t next_particle_id = particles.size();
    std::mt19937_64 breakage_generator(qmom_1d_breakage_advection_mc_initialization::breakage_seed);

    const qmomkit::PhaseSpaceVelocity velocity = qmom_1d_breakage_advection_mc_initialization::advection_velocity();
    const qmomkit::AdaptiveRkOptions advection_options;

    std::vector<double> final_moments =
        qmom_1d_breakage_advection_mc_initialization::compute_particle_moments(particles, particle_weight, maximum_order);
    const std::vector<double> initial_moments = final_moments;
    std::vector<std::tuple<double, std::size_t, double>> particle_history;
    if (output_options.write_moments)
    {
        for (const qmom_1d_breakage_advection_mc_initialization::Particle& particle : particles)
        {
            particle_history.emplace_back(initial_time, particle.id, particle.diameter);
        }
    }

    std::optional<qmomkit::MomentWriter> moment_history_writer;
    try
    {
        if (output_options.write_moments)
        {
            std::filesystem::create_directories(output_options.directory);
            moment_history_writer.emplace(qmomkit::MomentWriter::for_qmom(output_options.directory / "moment_history_mc.txt", final_moments.size()));
            moment_history_writer->write(initial_time, qmomkit::MomentSet(final_moments));
        }
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    const double output_interval = output_options.interval.value_or(final_time - initial_time);
    double time = initial_time;
    double next_output_time = initial_time + output_interval;
    while (time < final_time - 1.0e-12)
    {
        const double current_step = std::min(time_step, final_time - time);
        const double half_step = 0.5 * current_step;

        // Strang stage 1: advection half-step evolves particle diameters.
        qmomkit::monte_carlo::EnsembleEvolutionResult first_advection = qmomkit::monte_carlo::evolve_ensemble(
            qmom_1d_breakage_advection_mc_initialization::to_ensemble(particles), velocity, time, time + half_step, advection_options);
        if (!first_advection.ok())
        {
            std::cerr << "First advection particle half-step failed: " << first_advection.message << '\n';
            return 1;
        }

        try
        {
            qmom_1d_breakage_advection_mc_initialization::update_particle_diameters(particles, first_advection.ensemble);
        }
        catch (const std::exception& error)
        {
            std::cerr << error.what() << '\n';
            return 1;
        }

        // Strang stage 2: stochastic binary fragmentation acts on particles.
        qmom_1d_breakage_advection_mc_initialization::apply_binary_breakage(particles, current_step, breakage_generator, next_particle_id);

        // Strang stage 3: advection half-step evolves the post-breakage particles.
        qmomkit::monte_carlo::EnsembleEvolutionResult second_advection = qmomkit::monte_carlo::evolve_ensemble(
            qmom_1d_breakage_advection_mc_initialization::to_ensemble(particles), velocity, time + half_step, time + current_step, advection_options);
        if (!second_advection.ok())
        {
            std::cerr << "Second advection particle half-step failed: " << second_advection.message << '\n';
            return 1;
        }

        try
        {
            qmom_1d_breakage_advection_mc_initialization::update_particle_diameters(particles, second_advection.ensemble);
        }
        catch (const std::exception& error)
        {
            std::cerr << error.what() << '\n';
            return 1;
        }

        time += current_step;

        if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
        {
            final_moments = qmom_1d_breakage_advection_mc_initialization::compute_particle_moments(particles, particle_weight, maximum_order);
            if (moment_history_writer.has_value())
            {
                try
                {
                    moment_history_writer->write(time, qmomkit::MomentSet(final_moments));
                    for (const qmom_1d_breakage_advection_mc_initialization::Particle& particle : particles)
                    {
                        particle_history.emplace_back(time, particle.id, particle.diameter);
                    }
                }
                catch (const std::exception& error)
                {
                    std::cerr << error.what() << '\n';
                    return 1;
                }
            }
            while (next_output_time <= time + 1.0e-12)
            {
                next_output_time += output_interval;
            }
        }
    }

    const std::filesystem::path particle_history_file = output_options.directory / "particle_history.txt";

    try
    {
        if (output_options.write_moments)
        {
            qmom_1d_breakage_advection_mc_initialization::write_particle_history(particle_history_file, particle_history);
        }
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    std::cout << std::setprecision(12);
    qmom_1d_breakage_advection_mc_initialization::print_common_run_summary(
        std::cout, "Monte Carlo particles", configuration_file, output_options.directory,
        qmom_1d_breakage_advection_mc_initialization::initialization_sample_count, qmom_1d_breakage_advection_mc_initialization::initialization_seed,
        initial_time, final_time, time_step, output_interval);
    std::cout << "Representation: sampled particles with constant statistical weight\n";
    std::cout << "Breakage seed: " << qmom_1d_breakage_advection_mc_initialization::breakage_seed << '\n';
    std::cout << "Fragmentation rate: " << qmom_1d_breakage_advection_mc_initialization::fragmentation_rate << '\n';

    qmom_1d_breakage_advection_mc_initialization::print_moments(std::cout, "Initial moments", initial_moments);
    qmom_1d_breakage_advection_mc_initialization::print_moments(std::cout, "Final moments", final_moments);

    std::cout << "\nFinal ensemble at t = " << time << ":\n";
    std::cout << "initial_particle_count = " << qmom_1d_breakage_advection_mc_initialization::initialization_sample_count << '\n';
    std::cout << "final_particle_count = " << particles.size() << '\n';
    if (output_options.write_moments)
    {
        std::cout << "\nCompleted Monte Carlo particle workflow. Outputs written to " << output_options.directory << '\n';
    }
    else
    {
        std::cout << "\nCompleted Monte Carlo particle workflow. Output disabled by configuration.\n";
    }

    return 0;
}
