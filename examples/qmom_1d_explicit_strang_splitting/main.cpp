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
    const double final_time = 0.5;
    const double time_step = 0.1;
    const std::size_t maximum_node_count = 2;

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
    qmomkit::QmomPhaseSpaceEvolutionOptions phase_space_evolution_options = std::move(phase_space_evolution_configured).value();
    if (qmom_case.node_count != maximum_node_count)
    {
        std::cerr << "The explicit Strang example expects configuration node_counts = 2.\n";
        return 1;
    }

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
    const qmomkit::PhaseSpaceVelocity velocity = qmom_1d_explicit_strang_splitting::advection_velocity();
    const qmomkit::MomentSource source = qmom_1d_explicit_strang_splitting::variance_injection_source();

    std::optional<qmomkit::QuadratureWriter> quadrature_history_writer;
    try
    {
        if (output_options.write_quadrature)
        {
            std::filesystem::create_directories(output_options.directory);
            quadrature_history_writer.emplace(
                qmomkit::QuadratureWriter::for_qmom(output_options.directory / "quadrature_history.txt", maximum_node_count));
            quadrature_history_writer->write(initial_time, state.quadrature());
        }
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    std::cout << std::setprecision(12);
    std::cout << "initial active nodes: " << state.active_node_count() << '\n';

    double time = initial_time;
    while (time < final_time - 1.0e-12)
    {
        const double step = std::min(time_step, final_time - time);
        const double half_step = 0.5 * step;
        // advance_qmom_phase_space(...) would compose these three stages for you.
        // Driving them by hand is worth it only when the application needs to do
        // its own work between them - inspect an intermediate state, couple to
        // another solver, or apply a different acceptance rule. The composition
        // below is exactly what the `strang` scheme does internally.
        qmomkit::QmomState trial_state = state;

        // First advection half-step.
        const qmomkit::QmomPhaseSpaceEvolutionResult first_advection =
            qmomkit::advance_qmom_advection(trial_state, time, half_step, velocity, phase_space_evolution_options);
        if (!first_advection.ok())
        {
            std::cerr << "QMOM first advection stage failed: " << first_advection.message << '\n';
            return 1;
        }
        trial_state = first_advection.state;

        // Full moment-source stage.
        const std::size_t before_moment_source_nodes = trial_state.active_node_count();
        const qmomkit::QmomPhaseSpaceEvolutionResult moment_source =
            qmomkit::advance_qmom_moment_source(trial_state, time + half_step, step, source, phase_space_evolution_options);
        if (!moment_source.ok())
        {
            std::cerr << "QMOM moment-source stage failed: " << moment_source.message << '\n';
            return 1;
        }
        trial_state = moment_source.state;
        const std::size_t after_moment_source_nodes = trial_state.active_node_count();

        // Second advection half-step.
        const qmomkit::QmomPhaseSpaceEvolutionResult second_advection =
            qmomkit::advance_qmom_advection(trial_state, time + half_step, half_step, velocity, phase_space_evolution_options);
        if (!second_advection.ok())
        {
            std::cerr << "QMOM second advection stage failed: " << second_advection.message << '\n';
            return 1;
        }
        trial_state = second_advection.state;

        state = trial_state;
        time += step;

        if (before_moment_source_nodes != after_moment_source_nodes)
        {
            std::cout << "t = " << time - half_step << ": moment-source reconstruction changed active nodes " << before_moment_source_nodes << " -> "
                      << after_moment_source_nodes << '\n';
        }

        try
        {
            if (quadrature_history_writer.has_value())
            {
                quadrature_history_writer->write(time, state.quadrature());
            }
        }
        catch (const std::exception& error)
        {
            std::cerr << error.what() << '\n';
            return 1;
        }
    }

    std::cout << "final active nodes: " << state.active_node_count() << '\n';

    return 0;
}
