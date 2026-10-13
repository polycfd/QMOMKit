#include "physics.hpp"
#include "qmomkit/qmomkit.hpp"
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
    const double final_time = 8.0;
    const double time_step = 0.002;

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
    const qmomkit::CqmomPreparedCase cqmom_case = std::move(configured_case).value();
    const qmomkit::CqmomPhaseSpaceEvolutionOptions options = std::move(configured_options).value();
    const qmomkit::OutputOptions output = std::move(configured_output).value();

    // The physics lives in physics.cpp: an advection velocity that relaxes the
    // velocity distribution toward a single value, and a moment source that
    // injects velocity variance in a short pulse. Relaxation lets a velocity
    // branch lose a node; the pulse lets it grow one back.
    const cqmom_2d_adaptive_node_counts::VarianceInjectionParameters parameters;
    const qmomkit::PhaseSpaceVelocity velocity = cqmom_2d_adaptive_node_counts::advection_velocity(parameters);
    const qmomkit::CqmomMomentSource source = cqmom_2d_adaptive_node_counts::variance_injection_source(parameters);

    try
    {
        // The configuration enables adaptive node counts: each branch uses as many
        // nodes as its moments support. initial_state() throws when the moments
        // admit no quadrature, and the writers throw when a file cannot be
        // written.
        qmomkit::CqmomState state = cqmom_case.initial_state();

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

        std::cout << "Node-count changes of the velocity branches:\n";
        const double output_interval = output.interval.value_or(final_time - initial_time);
        double next_output_time = initial_time + output_interval;
        double time = initial_time;
        while (time < final_time - 1.0e-12)
        {
            qmomkit::CqmomPhaseSpaceEvolutionResult result =
                qmomkit::advance_cqmom_phase_space(state, time, std::min(time_step, final_time - time), velocity, source, options);
            if (!result.ok())
            {
                std::cerr << "CQMOM phase-space evolution failed: " << result.message << '\n';
                return 1;
            }
            state = std::move(result.state);
            time = result.final_time;

            // The diagnostics of each step report where a branch changed its node
            // count, which is what this example exists to show.
            for (const qmomkit::CqmomNodeCountTransition& transition : result.diagnostics.cqmom_node_count_transitions)
            {
                if (transition.phase_space_dimension == cqmom_2d_adaptive_node_counts::velocity_dimension)
                {
                    cqmom_2d_adaptive_node_counts::print_transition(std::cout, time, transition);
                }
            }

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
