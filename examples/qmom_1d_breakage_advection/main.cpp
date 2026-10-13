#include "physics.hpp"
#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <utility>

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

    // The physics is supplied as two callbacks defined in physics.cpp: a velocity
    // that moves the quadrature nodes, and a fragmentation source that returns
    // moment rates. QMOMKit composes them; it never interprets the coordinate.
    const qmomkit::PhaseSpaceVelocity velocity = qmom_1d_breakage_advection::advection_velocity();
    const qmomkit::QmomMomentSource source = qmom_1d_breakage_advection::binary_fragmentation_source();

    try
    {
        // initial_state() throws when the configured moments admit no quadrature,
        // and the writers throw when a file cannot be written.
        qmomkit::QmomState state = qmom_case.initial_state();
        const double initial_number_density = state.moment_values()[0];
        const double initial_mean_diameter = state.moment_values()[1] / initial_number_density;

        std::optional<qmomkit::MomentWriter> moment_writer;
        std::optional<qmomkit::QuadratureWriter> quadrature_writer;
        if (output.write_moments || output.write_quadrature)
        {
            std::filesystem::create_directories(output.directory);
        }
        if (output.write_moments)
        {
            moment_writer.emplace(qmomkit::MomentWriter::for_qmom(output.directory / "moments.txt", state, output));
        }
        if (output.write_quadrature)
        {
            quadrature_writer.emplace(qmomkit::QuadratureWriter::for_qmom(output.directory / "quadrature.txt", state.maximum_node_count()));
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
            // One step of the configured splitting scheme: advection of the nodes
            // and the fragmentation source on the moments.
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

        const double final_number_density = state.moment_values()[0];
        std::cout << "Number density M0:    " << initial_number_density << " -> " << final_number_density << '\n';
        std::cout << "Mean diameter M1/M0:  " << initial_mean_diameter << " -> " << state.moment_values()[1] / final_number_density << '\n';
        std::cout << "Final quadrature at t = " << time << ":\n";
        const qmomkit::QmomResult& quadrature = state.quadrature();
        for (std::size_t node = 0; node < quadrature.weights.size(); ++node)
        {
            std::cout << "  weight " << quadrature.weights[node] << " at d = " << quadrature.abscissas[node] << '\n';
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
