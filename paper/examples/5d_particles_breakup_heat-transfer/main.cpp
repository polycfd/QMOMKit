#include "monte_carlo/evolution.hpp"
#include "monte_carlo/output.hpp"
#include "physics.hpp"
#include "qmomkit/qmomkit.hpp"
#include "sampling.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <random>
#include <utility>
#include <vector>

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";
    const double initial_time = 0.0;
    const double final_time = 0.1;
    const double time_step = 0.001;
    constexpr unsigned long long breakup_seed = 321;

    // The case, the evolution options and the output options come from the same
    // file. A faulty file is reported as a message naming the key and line.
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
    const qmomkit::CqmomPhaseSpaceEvolutionOptions options = std::move(configured_options).value();
    const qmomkit::OutputOptions output = std::move(configured_output).value();

    // The CQMOM run starts from the moments of a finite sample of the configured
    // distribution, and the Monte Carlo reference advances that same sample.
    const std::optional<qmomkit::monte_carlo::TensorSampledInitialization> sampled =
        paper_5d_particles_breakup_heat_transfer::sample_initial_distribution(configuration_file, configured_case.value().hierarchy());
    if (!sampled)
    {
        return 1;
    }
    auto sampled_case = qmomkit::config::with_initial_moments(configured_case.value(), sampled->moments.moments);
    if (!sampled_case)
    {
        std::cerr << sampled_case.error() << '\n';
        return 1;
    }
    const qmomkit::CqmomPreparedCase cqmom_case = std::move(sampled_case).value();
    qmomkit::monte_carlo::Ensemble ensemble = sampled->sampling.ensemble;
    const double particle_weight = sampled->sampling.number_density / static_cast<double>(ensemble.size());
    const std::vector<std::size_t> moment_extents = sampled->moments.moments.extents;
    std::mt19937_64 generator(breakup_seed);

    const qmomkit::PhaseSpaceVelocity velocity = paper_5d_particles_breakup_heat_transfer::particle_velocity();
    const qmomkit::CqmomMomentSource source = paper_5d_particles_breakup_heat_transfer::breakup_source();
    const qmomkit::AdaptiveRkOptions& mc_advection_options = options.node_options.adaptive_options;

    try
    {
        // initial_state() throws when the sampled moments admit no quadrature, and
        // the writers and the Monte Carlo estimates throw when they fail.
        qmomkit::CqmomState state = cqmom_case.initial_state();

        std::optional<qmomkit::MomentWriter> moment_writer;
        std::optional<qmomkit::MomentWriter> mc_moment_writer;
        std::optional<qmomkit::QuadratureWriter> quadrature_writer;
        if (output.write_moments || output.write_quadrature)
        {
            std::filesystem::create_directories(output.directory);
        }
        if (output.write_moments)
        {
            // Both histories hold the moments the state tracks, followed by the
            // generalized moments configured in [output].
            moment_writer.emplace(qmomkit::MomentWriter::for_cqmom(output.directory / "moment_history.txt", state, output));
            mc_moment_writer.emplace(qmomkit::monte_carlo::make_moment_history_writer(output.directory / "moment_history_mc.txt", moment_extents,
                                                                                      state.moment_indices(), output));
        }
        if (output.write_quadrature)
        {
            quadrature_writer.emplace(
                qmomkit::QuadratureWriter::for_cqmom(output.directory / "quadrature_history.txt", state.hierarchy().number_of_dimensions()));
        }
        const auto write_output = [&](const double time)
        {
            if (moment_writer)
            {
                moment_writer->write(time, state, state.quadrature());
                qmomkit::monte_carlo::write_moment_history_row(
                    *mc_moment_writer, time, paper_5d_particles_breakup_heat_transfer::ensemble_moments(ensemble, particle_weight, moment_extents),
                    ensemble, particle_weight * static_cast<double>(ensemble.size()), output);
            }
            if (quadrature_writer)
            {
                quadrature_writer->write(time, state.quadrature());
            }
        };
        write_output(initial_time);

        const double output_interval = output.interval.value_or(final_time - initial_time);
        double next_output_time = initial_time + output_interval;
        double time = initial_time;
        while (time < final_time - 1.0e-12)
        {
            // Strang splitting for both runs: an advection half-step, breakup over
            // the whole step, and a second advection half-step. The particles break
            // up stochastically.
            const double step = std::min(time_step, final_time - time);
            qmomkit::CqmomPhaseSpaceEvolutionResult result = qmomkit::advance_cqmom_phase_space(state, time, step, velocity, source, options);
            if (!result.ok())
            {
                std::cerr << "CQMOM phase-space evolution failed: " << result.message << '\n';
                return 1;
            }
            state = std::move(result.state);

            const double half_step = 0.5 * step;
            qmomkit::monte_carlo::EnsembleEvolutionResult first_advection =
                qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, time, time + half_step, mc_advection_options);
            if (!first_advection.ok())
            {
                std::cerr << "First advection particle half-step failed: " << first_advection.message << '\n';
                return 1;
            }
            ensemble = std::move(first_advection.ensemble);

            paper_5d_particles_breakup_heat_transfer::apply_breakup(ensemble, step, generator);

            qmomkit::monte_carlo::EnsembleEvolutionResult second_advection =
                qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, time + half_step, time + step, mc_advection_options);
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

        if (moment_writer || quadrature_writer)
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
