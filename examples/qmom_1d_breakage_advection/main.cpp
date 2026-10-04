#include "physics.hpp"
#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <iomanip>
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

    // Three loads from one file. load_qmom(...) validates every recognised key
    // but returns only the case, so the output and evolution settings have to be
    // requested separately.
    auto configured = qmomkit::config::load_qmom(configuration_file);
    if (!configured)
    {
        std::cerr << configured.error() << '\n';
        return 1;
    }

    auto output_configured = qmomkit::config::load_output_options(configuration_file);
    if (!output_configured)
    {
        std::cerr << output_configured.error() << '\n';
        return 1;
    }
    auto phase_space_evolution_configured = qmomkit::config::load_qmom_phase_space_evolution_options(configuration_file);
    if (!phase_space_evolution_configured)
    {
        std::cerr << phase_space_evolution_configured.error() << '\n';
        return 1;
    }

    const qmomkit::QmomCase qmom_case = std::move(configured).value();
    const qmomkit::OutputOptions output_options = std::move(output_configured).value();
    const std::vector<double> initial_moments = qmom_case.moments.values();
    // initial_state() inverts the configured moments, so it throws when they do
    // not admit a quadrature at the requested node count. Constructing the state
    // inside a try block keeps that failure separate from a configuration error.
    std::optional<qmomkit::QmomState> initial_state_storage;
    try
    {
        initial_state_storage.emplace(qmom_case.initial_state());
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    qmomkit::QmomState state = *initial_state_storage;
    const qmomkit::QmomResult initial_quadrature = state.quadrature();

    // The physics is supplied as two callbacks defined in physics.cpp: a velocity
    // that moves the quadrature abscissas, and a source that changes the moments
    // directly. QMOMKit composes them; it never interprets the coordinates.
    const qmomkit::PhaseSpaceVelocity velocity = qmom_1d_breakage_advection::advection_velocity();
    const qmomkit::QmomMomentSource source = qmom_1d_breakage_advection::binary_fragmentation_source();
    qmomkit::QmomPhaseSpaceEvolutionOptions phase_space_evolution_options = std::move(phase_space_evolution_configured).value();

    std::optional<qmomkit::MomentWriter> moment_writer;
    std::optional<qmomkit::QuadratureWriter> quadrature_writer;
    try
    {
        if (output_options.write_moments || output_options.write_quadrature)
        {
            std::filesystem::create_directories(output_options.directory);
        }
        if (output_options.write_moments)
        {
            moment_writer.emplace(qmomkit::MomentWriter::for_qmom(output_options.directory / "moments.txt", state, output_options));
            moment_writer->write(initial_time, state, state.quadrature());
        }
        if (output_options.write_quadrature)
        {
            quadrature_writer.emplace(
                qmomkit::QuadratureWriter::for_qmom(output_options.directory / "quadrature.txt", state.quadrature().weights.size()));
            quadrature_writer->write(initial_time, state.quadrature());
        }
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    double time = initial_time;
    const double output_interval = output_options.interval.value_or(final_time - initial_time);
    double next_output_time = initial_time + output_interval;
    while (time < final_time - 1.0e-12)
    {
        const double current_step = std::min(time_step, final_time - time);
        // One complete step: advection and moment source composed by the
        // configured splitting scheme. The returned state replaces the current
        // one only after ok() confirms the step was accepted.
        const qmomkit::QmomPhaseSpaceEvolutionResult phase_space_evolved =
            qmomkit::advance_qmom_phase_space(state, time, current_step, velocity, source, phase_space_evolution_options);
        if (!phase_space_evolved.ok())
        {
            std::cerr << "QMOM phase-space evolution failed: " << phase_space_evolved.message << '\n';
            return 1;
        }

        state = phase_space_evolved.state;
        time = phase_space_evolved.final_time;

        if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
        {
            try
            {
                if (moment_writer.has_value())
                {
                    moment_writer->write(time, state, state.quadrature());
                }
                if (quadrature_writer.has_value())
                {
                    quadrature_writer->write(time, state.quadrature());
                }
            }
            catch (const std::exception& error)
            {
                std::cerr << error.what() << '\n';
                return 1;
            }
            while (next_output_time <= time + 1.0e-12)
            {
                next_output_time += output_interval;
            }
        }
    }

    std::cout << std::setprecision(12);
    std::cout << "1D QMOM phase-space advection + moment-source evolution example\n";

    std::cout << '\n';
    std::cout << "Initial moments:\n";
    for (std::size_t order = 0; order < initial_moments.size(); ++order)
    {
        std::cout << "M" << order << " = " << initial_moments[order] << '\n';
    }

    std::cout << '\n';
    std::cout << "Final moments:\n";
    for (std::size_t order = 0; order < state.moment_values().size(); ++order)
    {
        std::cout << "M" << order << " = " << state.moment_values()[order] << '\n';
    }

    std::cout << '\n';
    std::cout << "Initial quadrature:\n";
    std::cout << "node  weight        abscissa\n";
    for (std::size_t node = 0; node < initial_quadrature.weights.size(); ++node)
    {
        std::cout << node << "     " << initial_quadrature.weights[node] << "   " << initial_quadrature.abscissas[node] << '\n';
    }

    std::cout << '\n';
    std::cout << "Final quadrature:\n";
    std::cout << "node  weight        abscissa\n";
    for (std::size_t node = 0; node < state.quadrature().weights.size(); ++node)
    {
        std::cout << node << "     " << state.quadrature().weights[node] << "   " << state.quadrature().abscissas[node] << '\n';
    }
    std::cout << "final_time = " << time << '\n';

    if (output_options.write_moments || output_options.write_quadrature)
    {
        std::cout << "\nOutputs written to:\n";
        if (output_options.write_moments)
        {
            std::cout << (output_options.directory / "moments.txt").string() << '\n';
        }
        if (output_options.write_quadrature)
        {
            std::cout << (output_options.directory / "quadrature.txt").string() << '\n';
        }
    }
    else
    {
        std::cout << "\nOutput disabled by configuration.\n";
    }

    return 0;
}
