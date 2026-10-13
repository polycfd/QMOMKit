#include "physics.hpp"
#include "qmomkit/qmomkit.hpp"
#include "sampling.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <utility>

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";
    const double initial_time = 0.0;
    const double final_time = 0.1;
    const double time_step = 0.001;

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

    // Both the deterministic run and its Monte Carlo reference start from one
    // sampled ensemble, so any later difference between them comes from the
    // evolution rather than from a different starting population. Only the
    // initial moments are replaced; the configured hierarchy, supports and
    // reconstruction settings are kept.
    const std::optional<qmomkit::monte_carlo::TensorSampledInitialization> sampled =
        cqmom_4d_particles_breakage_recoil::sample_initial_population(configuration_file, configured_case.value().hierarchy());
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

    // particle_velocity() relaxes each particle's velocity (u, v, w) towards the
    // gas by drag, and breakup_source() replaces breaking particles by two
    // daughters with diameter-scaled recoil.
    const qmomkit::PhaseSpaceVelocity velocity = cqmom_4d_particles_breakage_recoil::particle_velocity();
    const qmomkit::CqmomMomentSource source = cqmom_4d_particles_breakage_recoil::breakup_source();
    const qmomkit::RawMomentRequest printed_moments = cqmom_4d_particles_breakage_recoil::printed_moments();

    try
    {
        // initial_state() throws when the sampled moments admit no quadrature, and
        // the writers throw when a file cannot be written.
        qmomkit::CqmomState state = cqmom_case.initial_state();
        const qmomkit::SelectedRawMomentResult initial_moments = qmomkit::reconstruct_raw_moments(state.quadrature(), printed_moments);
        if (!initial_moments.ok())
        {
            std::cerr << "Initial moment reconstruction failed: " << initial_moments.message << '\n';
            return 1;
        }

        std::optional<qmomkit::MomentWriter> moment_writer;
        std::optional<qmomkit::QuadratureWriter> quadrature_writer;
        if (output.write_moments || output.write_quadrature)
        {
            std::filesystem::create_directories(output.directory);
        }
        if (output.write_moments)
        {
            moment_writer.emplace(qmomkit::MomentWriter::for_cqmom(output.directory / "moment_history.txt", state, output));
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
        while (time < final_time - 1.0e-14)
        {
            // One step of the configured splitting scheme in the four-dimensional
            // phase space (d, u, v, w).
            qmomkit::CqmomPhaseSpaceEvolutionResult result =
                qmomkit::advance_cqmom_phase_space(state, time, std::min(time_step, final_time - time), velocity, source, options);
            if (!result.ok())
            {
                std::cerr << "CQMOM phase-space evolution failed: " << result.message << '\n';
                return 1;
            }
            state = std::move(result.state);
            time = result.final_time;

            if (time + 1.0e-14 >= next_output_time || time >= final_time - 1.0e-14)
            {
                write_output(time);
                next_output_time += output_interval;
            }
        }

        const qmomkit::SelectedRawMomentResult final_moments = qmomkit::reconstruct_raw_moments(state.quadrature(), printed_moments);
        if (!final_moments.ok())
        {
            std::cerr << "Final moment reconstruction failed: " << final_moments.message << '\n';
            return 1;
        }
        cqmom_4d_particles_breakage_recoil::print_moments(std::cout, "QMOMKit run", initial_moments.values, final_moments.values);
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
