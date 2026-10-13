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
    const double final_time = 0.5;
    const double time_step = 0.1;

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

    const qmomkit::PhaseSpaceVelocity velocity = qmom_1d_explicit_strang_splitting::advection_velocity();
    const qmomkit::QmomMomentSource source = qmom_1d_explicit_strang_splitting::variance_injection_source();

    try
    {
        // initial_state() throws when the configured moments admit no quadrature,
        // and the writers throw when a file cannot be written.
        qmomkit::QmomState state = qmom_case.initial_state();

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

        std::cout << "initial active nodes: " << state.active_node_count() << '\n';

        double time = initial_time;
        while (time < final_time - 1.0e-12)
        {
            const double step = std::min(time_step, final_time - time);
            const double half_step = 0.5 * step;
            // advance_qmom_phase_space(...) would compose these three stages for
            // you. Driving them by hand is worth it only when the application needs
            // to do its own work between them - inspect an intermediate state,
            // couple to another solver, or apply a different acceptance rule. The
            // composition below is exactly what the `strang` scheme does internally.
            const qmomkit::QmomPhaseSpaceEvolutionResult first_advection = qmomkit::advance_qmom_advection(state, time, half_step, velocity, options);
            if (!first_advection.ok())
            {
                std::cerr << "QMOM first advection stage failed: " << first_advection.message << '\n';
                return 1;
            }

            const qmomkit::QmomPhaseSpaceEvolutionResult moment_source =
                qmomkit::advance_qmom_moment_source(first_advection.state, time, step, source, options);
            if (!moment_source.ok())
            {
                std::cerr << "QMOM moment-source stage failed: " << moment_source.message << '\n';
                return 1;
            }

            // The work between stages that motivates driving them by hand: here,
            // reporting when the reconstruction after the source changes the node count.
            const std::size_t nodes_before_source = first_advection.state.active_node_count();
            const std::size_t nodes_after_source = moment_source.state.active_node_count();
            if (nodes_before_source != nodes_after_source)
            {
                std::cout << "t = " << time + half_step << ": moment-source reconstruction changed active nodes " << nodes_before_source << " -> "
                          << nodes_after_source << '\n';
            }

            qmomkit::QmomPhaseSpaceEvolutionResult second_advection =
                qmomkit::advance_qmom_advection(moment_source.state, time + half_step, half_step, velocity, options);
            if (!second_advection.ok())
            {
                std::cerr << "QMOM second advection stage failed: " << second_advection.message << '\n';
                return 1;
            }

            state = std::move(second_advection.state);
            time += step;
            write_output(time);
        }

        std::cout << "final active nodes: " << state.active_node_count() << '\n';
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
