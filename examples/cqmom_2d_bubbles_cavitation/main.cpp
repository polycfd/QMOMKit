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
    const double final_time = 13.9;
    const double time_step = 0.01;

    // The case, the advection options and the output options come from the same
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

    // Starting from the moments of the sample rather than from the configured
    // analytical moments is what makes the two runs comparable: both begin from
    // the same finite population. The case keeps its hierarchy and settings.
    const std::optional<qmomkit::monte_carlo::TensorSampledInitialization> sampled =
        cqmom_2d_bubbles_cavitation::sample_initial_distribution(configuration_file, configured_case.value().hierarchy());
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

    const qmomkit::PhaseSpaceVelocity velocity = cqmom_2d_bubbles_cavitation::advection_velocity();
    const qmomkit::RawMomentRequest printed_moments = cqmom_2d_bubbles_cavitation::printed_moments();

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
        while (time < final_time - 1.0e-12)
        {
            // The deterministic half of the comparison. The bubbles only move in
            // phase space and there is no moment source, so the step is the
            // advection stage alone rather than a split step. The companion
            // executable moves the sampled bubbles directly.
            qmomkit::CqmomPhaseSpaceEvolutionResult result =
                qmomkit::advance_cqmom_advection(state, time, std::min(time_step, final_time - time), velocity, options);
            if (!result.ok())
            {
                std::cerr << "CQMOM advection failed: " << result.message << '\n';
                return 1;
            }
            state = std::move(result.state);
            time = result.final_time;

            if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
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
        cqmom_2d_bubbles_cavitation::print_moments(std::cout, "QMOMKit run", initial_moments.values, final_moments.values);
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
