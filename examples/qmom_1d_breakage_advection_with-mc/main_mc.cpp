#include "monte_carlo/evolution.hpp"
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
#include <random>
#include <tuple>
#include <utility>
#include <vector>

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";
    const double initial_time = 0.0;
    const double final_time = 1.0;
    const double time_step = 0.01;

    auto configured_output = qmomkit::config::load_output_options(configuration_file);
    if (!configured_output)
    {
        std::cerr << configured_output.error() << '\n';
        return 1;
    }
    const qmomkit::OutputOptions output = std::move(configured_output).value();

    // The same sample the QMOMKit run starts from, here kept as particles of
    // equal statistical weight.
    const std::optional<qmom_1d_breakage_advection_with_mc::SamplingSetup> sampled =
        qmom_1d_breakage_advection_with_mc::sample_initial_distribution(configuration_file);
    if (!sampled)
    {
        return 1;
    }
    std::vector<qmom_1d_breakage_advection_with_mc::Particle> particles =
        qmom_1d_breakage_advection_with_mc::make_particles(sampled->sampling.ensemble);
    const std::size_t initial_particle_count = particles.size();
    const double particle_weight = sampled->sampling.number_density / static_cast<double>(initial_particle_count);
    const std::size_t maximum_order = sampled->moments.moments.values().size() - 1;
    std::size_t next_particle_id = particles.size();
    std::mt19937_64 breakage_generator(qmom_1d_breakage_advection_with_mc::breakage_seed);

    const qmomkit::PhaseSpaceVelocity velocity = qmom_1d_breakage_advection_with_mc::advection_velocity();
    const qmomkit::AdaptiveRkOptions advection_options;
    const std::vector<double> initial_moments =
        qmom_1d_breakage_advection_with_mc::compute_particle_moments(particles, particle_weight, maximum_order);

    try
    {
        // The writers throw when a file cannot be written. The particle history is
        // kept with the moment history and written at the end.
        std::optional<qmomkit::MomentWriter> moment_writer;
        std::vector<std::tuple<double, std::size_t, double>> particle_history;
        if (output.write_moments)
        {
            std::filesystem::create_directories(output.directory);
            moment_writer.emplace(
                qmomkit::monte_carlo::make_moment_history_writer(output.directory / "moment_history_mc.txt", initial_moments.size(), output));
        }
        const auto write_output = [&](const double time)
        {
            if (moment_writer)
            {
                const std::vector<double> moments =
                    qmom_1d_breakage_advection_with_mc::compute_particle_moments(particles, particle_weight, maximum_order);
                qmomkit::monte_carlo::write_moment_history_row(*moment_writer, time, qmomkit::MomentSet(moments),
                                                               qmom_1d_breakage_advection_with_mc::to_ensemble(particles),
                                                               particle_weight * static_cast<double>(particles.size()), output);
                for (const qmom_1d_breakage_advection_with_mc::Particle& particle : particles)
                {
                    particle_history.emplace_back(time, particle.id, particle.volume);
                }
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

            // Strang splitting, as in the QMOMKit run: an advection half-step of the
            // particle volumes, stochastic binary fragmentation over the whole
            // step, and a second advection half-step.
            const qmomkit::monte_carlo::EnsembleEvolutionResult first_advection = qmomkit::monte_carlo::evolve_ensemble(
                qmom_1d_breakage_advection_with_mc::to_ensemble(particles), velocity, time, time + half_step, advection_options);
            if (!first_advection.ok())
            {
                std::cerr << "First advection particle half-step failed: " << first_advection.message << '\n';
                return 1;
            }
            qmom_1d_breakage_advection_with_mc::update_particle_volumes(particles, first_advection.ensemble);

            qmom_1d_breakage_advection_with_mc::apply_binary_breakage(particles, step, breakage_generator, next_particle_id);

            const qmomkit::monte_carlo::EnsembleEvolutionResult second_advection = qmomkit::monte_carlo::evolve_ensemble(
                qmom_1d_breakage_advection_with_mc::to_ensemble(particles), velocity, time + half_step, time + step, advection_options);
            if (!second_advection.ok())
            {
                std::cerr << "Second advection particle half-step failed: " << second_advection.message << '\n';
                return 1;
            }
            qmom_1d_breakage_advection_with_mc::update_particle_volumes(particles, second_advection.ensemble);

            time += step;
            if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
            {
                write_output(time);
                next_output_time += output_interval;
            }
        }
        if (output.write_moments)
        {
            qmom_1d_breakage_advection_with_mc::write_particle_history(output.directory / "particle_history.txt", particle_history);
        }

        const std::vector<double> final_moments =
            qmom_1d_breakage_advection_with_mc::compute_particle_moments(particles, particle_weight, maximum_order);
        qmom_1d_breakage_advection_with_mc::print_moments(std::cout, "Monte Carlo run", initial_moments, final_moments);
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
