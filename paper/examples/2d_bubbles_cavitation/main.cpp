#include "monte_carlo/evolution.hpp"
#include "particles_mc.hpp"
#include "physics.hpp"
#include "qmomkit/qmomkit.hpp"
#include "sampling.hpp"

#include <algorithm>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <utility>
#include <vector>

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";
    const double initial_time = 0.0;
    const double final_time = 13.9;
    const double time_step = 0.139;

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
        paper_2d_bubbles_cavitation::sample_initial_distribution(configuration_file, configured_case.value().hierarchy());
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

    // The bubbles only move in phase space, so each step is the advection stage
    // alone, for the quadrature nodes and for the sampled bubbles alike.
    const qmomkit::PhaseSpaceVelocity velocity = paper_2d_bubbles_cavitation::advection_velocity();
    const qmomkit::AdaptiveRkOptions& mc_advection_options = options.node_options.adaptive_options;

    try
    {
        // initial_state() throws when the sampled moments admit no quadrature, and
        // the writers throw when a file cannot be written.
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
            // The QMOMKit history holds the moments the state tracks; the Monte
            // Carlo history holds the full moment tensor of the sample.
            moment_writer.emplace(qmomkit::MomentWriter::for_cqmom(output.directory / "moment_history.txt", state, output));
            mc_moment_writer.emplace(qmomkit::MomentWriter::for_cqmom(output.directory / "moment_history_mc.txt", moment_extents));
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
                mc_moment_writer->write(time, paper_2d_bubbles_cavitation::compute_particle_moments(ensemble, particle_weight, moment_extents));
            }
            if (quadrature_writer)
            {
                quadrature_writer->write(time, state.quadrature());
            }
        };
        write_output(initial_time);

        double time = initial_time;
        while (time < final_time - 1.0e-12)
        {
            const double step = std::min(time_step, final_time - time);
            qmomkit::CqmomPhaseSpaceEvolutionResult result = qmomkit::advance_cqmom_advection(state, time, step, velocity, options);
            if (!result.ok())
            {
                std::cerr << "CQMOM advection failed: " << result.message << '\n';
                return 1;
            }
            state = std::move(result.state);

            qmomkit::monte_carlo::EnsembleEvolutionResult mc_advection =
                qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, time, time + step, mc_advection_options);
            if (!mc_advection.ok())
            {
                std::cerr << "Bubble advection failed: " << mc_advection.message << '\n';
                return 1;
            }
            ensemble = std::move(mc_advection.ensemble);

            time += step;
            write_output(time);
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
