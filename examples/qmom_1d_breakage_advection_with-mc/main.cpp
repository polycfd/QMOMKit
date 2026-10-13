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
#include <vector>

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";
    const double initial_time = 0.0;
    const double final_time = 1.0;
    const double time_step = 0.01;

    // The case, the evolution options and the output options come from the same
    // file. A faulty file is reported as a message naming the key and line.
    auto configured_case = qmomkit::config::load_qmom(configuration_file);
    if (!configured_case)
    {
        std::cerr << configured_case.error() << '\n';
        return 1;
    }
    auto configured_options = qmomkit::config::load_qmom_phase_space_evolution_options(configuration_file);
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
    const qmomkit::QmomCase qmom_case = std::move(configured_case).value();
    const qmomkit::QmomPhaseSpaceEvolutionOptions options = std::move(configured_options).value();
    const qmomkit::OutputOptions output = std::move(configured_output).value();

    // Starting from the moments of the sample rather than from the configured
    // analytical moments is what makes the two runs comparable: both begin from
    // the same finite population. The case keeps the configured node count and
    // inversion options.
    const std::optional<qmom_1d_breakage_advection_with_mc::SamplingSetup> sampled =
        qmom_1d_breakage_advection_with_mc::sample_initial_distribution(configuration_file);
    if (!sampled)
    {
        return 1;
    }
    const qmomkit::QmomCase sampled_case{sampled->moments.moments, qmom_case.node_count, qmom_case.qmom_options};
    const std::vector<double> initial_moments = sampled->moments.moments.values();

    const qmomkit::PhaseSpaceVelocity velocity = qmom_1d_breakage_advection_with_mc::advection_velocity();
    const qmomkit::QmomMomentSource source = qmom_1d_breakage_advection_with_mc::binary_fragmentation_source();

    try
    {
        // initial_state() throws when the sampled moments admit no quadrature, and
        // the writers throw when a file cannot be written.
        qmomkit::QmomState state = sampled_case.initial_state();

        std::optional<qmomkit::MomentWriter> moment_writer;
        std::optional<qmomkit::QuadratureWriter> quadrature_writer;
        if (output.write_moments || output.write_quadrature)
        {
            std::filesystem::create_directories(output.directory);
        }
        if (output.write_moments)
        {
            moment_writer.emplace(qmomkit::MomentWriter::for_qmom(output.directory / "moment_history.txt", state, output));
        }
        if (output.write_quadrature)
        {
            quadrature_writer.emplace(qmomkit::QuadratureWriter::for_qmom(output.directory / "quadrature_history.txt", state.maximum_node_count()));
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
            // The deterministic half of the comparison. The Monte Carlo companion
            // advances the same physics over the same interval.
            qmomkit::QmomPhaseSpaceEvolutionResult result =
                qmomkit::advance_qmom_phase_space(state, time, std::min(time_step, final_time - time), velocity, source, options);
            if (!result.ok())
            {
                std::cerr << "QMOM phase-space evolution failed: " << result.message << '\n';
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

        qmom_1d_breakage_advection_with_mc::print_moments(std::cout, "QMOMKit run", initial_moments, state.moment_values());
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
