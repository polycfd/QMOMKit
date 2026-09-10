#include "physics.hpp"
#include "qmomkit/qmomkit.hpp"
#include "sampling.hpp"
#include "utilities.hpp"

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

    const std::optional<std::string> configuration_text = qmom_1d_breakage_advection_with_mc::read_configuration_text(configuration_file);
    if (!configuration_text.has_value())
    {
        return 1;
    }

    auto configured = qmomkit::config::load_qmom(configuration_file);
    if (!configured)
    {
        qmom_1d_breakage_advection_with_mc::print_configuration_errors(configured.error());
        return 1;
    }

    auto output_configured = qmomkit::config::load_output_options(configuration_file);
    if (!output_configured)
    {
        qmom_1d_breakage_advection_with_mc::print_configuration_errors(output_configured.error());
        return 1;
    }
    auto phase_space_evolution_configured = qmomkit::config::load_qmom_phase_space_evolution_options(configuration_file);
    if (!phase_space_evolution_configured)
    {
        qmom_1d_breakage_advection_with_mc::print_configuration_errors(phase_space_evolution_configured.error());
        return 1;
    }

    const qmomkit::QmomCase qmom_case = std::move(configured).value();
    const qmomkit::OutputOptions output_options = std::move(output_configured).value();
    auto initial_condition = qmomkit::config::parse_qmom_initial_condition(*configuration_text, configuration_file.string());
    if (!initial_condition)
    {
        qmom_1d_breakage_advection_with_mc::print_configuration_errors(initial_condition.error());
        return 1;
    }

    const std::optional<qmom_1d_breakage_advection_with_mc::SamplingSetup> sampled =
        qmom_1d_breakage_advection_with_mc::sample_initial_distribution(std::move(initial_condition).value());
    if (!sampled.has_value())
    {
        return 1;
    }
    if (!sampled->sampling.ok())
    {
        std::cerr << sampled->sampling.message << '\n';
        return 1;
    }

    const qmomkit::monte_carlo::OneDimensionalMomentEstimate& sampled_moments = sampled->moments;
    if (!sampled_moments.ok())
    {
        std::cerr << sampled_moments.message << '\n';
        return 1;
    }

    const std::vector<double> initial_moments = sampled_moments.moments.values();
    std::vector<double> moments = initial_moments;
    const qmomkit::QmomCase sampled_case{sampled_moments.moments, qmom_case.node_count, qmom_case.options};
    std::optional<qmomkit::QmomState> initial_state_storage;
    try
    {
        initial_state_storage.emplace(sampled_case.initial_state());
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    qmomkit::QmomState state = *initial_state_storage;
    moments = state.moments();

    const qmomkit::QmomResult initial_quadrature = state.quadrature();
    qmomkit::QmomResult quadrature = initial_quadrature;
    const qmomkit::PhaseSpaceVelocity velocity = qmom_1d_breakage_advection_with_mc::advection_velocity();
    const qmomkit::MomentSource source = qmom_1d_breakage_advection_with_mc::binary_fragmentation_source();
    qmomkit::QmomPhaseSpaceEvolutionOptions phase_space_evolution_options = std::move(phase_space_evolution_configured).value();
    phase_space_evolution_options.qmom_options = sampled_case.options;
    std::optional<qmomkit::MomentWriter> moment_history_writer;
    std::optional<qmomkit::QuadratureWriter> quadrature_history_writer;
    try
    {
        if (output_options.write_moments || output_options.write_quadrature)
        {
            std::filesystem::create_directories(output_options.directory);
        }
        if (output_options.write_moments)
        {
            qmomkit::MomentWriter initial_moment_writer =
                qmomkit::MomentWriter::for_qmom(output_options.directory / "initial_moments.txt", state, output_options);
            initial_moment_writer.write(initial_time, state, quadrature);

            moment_history_writer.emplace(qmomkit::MomentWriter::for_qmom(output_options.directory / "moment_history.txt", state, output_options));
            moment_history_writer->write(initial_time, state, quadrature);
        }
        if (output_options.write_quadrature)
        {
            qmomkit::QuadratureWriter initial_quadrature_writer =
                qmomkit::QuadratureWriter::for_qmom(output_options.directory / "initial_quadrature.txt", quadrature.weights.size());
            initial_quadrature_writer.write(initial_time, quadrature);

            quadrature_history_writer.emplace(
                qmomkit::QuadratureWriter::for_qmom(output_options.directory / "quadrature_history.txt", quadrature.weights.size()));
            quadrature_history_writer->write(initial_time, quadrature);
        }
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    const double output_interval = output_options.interval.value_or(final_time - initial_time);
    double time = initial_time;
    double next_output_time = initial_time + output_interval;
    while (time < final_time - 1.0e-12)
    {
        const double current_step = std::min(time_step, final_time - time);
        const qmomkit::QmomPhaseSpaceEvolutionResult phase_space_evolved =
            qmomkit::advance_qmom_phase_space(state, sampled_case.node_count, time, current_step, velocity, source, phase_space_evolution_options);
        if (!phase_space_evolved.ok())
        {
            std::cerr << "QMOM phase-space evolution failed: " << phase_space_evolved.message << '\n';
            return 1;
        }

        state = phase_space_evolved.state;
        moments = state.moments();
        quadrature = state.quadrature();
        time = phase_space_evolved.final_time;

        if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
        {
            try
            {
                if (moment_history_writer.has_value())
                {
                    moment_history_writer->write(time, state, quadrature);
                }
                if (quadrature_history_writer.has_value())
                {
                    quadrature_history_writer->write(time, quadrature);
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

    qmom_1d_breakage_advection_with_mc::print_common_run_summary(
        std::cout, "QMOM quadrature", configuration_file, output_options.directory, qmom_1d_breakage_advection_with_mc::initialization_sample_count,
        qmom_1d_breakage_advection_with_mc::initialization_seed, initial_time, final_time, time_step, output_interval);
    std::cout << "Representation: sampled moments reconstructed as QMOM nodes and weights\n";

    qmom_1d_breakage_advection_with_mc::print_moments(std::cout, "Initial moments", initial_moments);

    std::cout << "\nInitial quadrature:\n";
    std::cout << "node  weight        abscissa\n";
    for (std::size_t node = 0; node < initial_quadrature.weights.size(); ++node)
    {
        std::cout << node << "     " << initial_quadrature.weights[node] << "   " << initial_quadrature.abscissas[node] << '\n';
    }

    qmom_1d_breakage_advection_with_mc::print_moments(std::cout, "Final moments", moments);

    std::cout << "\nFinal quadrature at t = " << time << ":\n";
    std::cout << "node  initial_abscissa  final_abscissa\n";
    for (std::size_t node = 0; node < quadrature.abscissas.size(); ++node)
    {
        std::cout << node << "     " << initial_quadrature.abscissas[node] << "         " << quadrature.abscissas[node] << '\n';
    }
    if (output_options.write_moments || output_options.write_quadrature)
    {
        std::cout << "\nCompleted QMOM workflow. Outputs written to " << output_options.directory << '\n';
    }
    else
    {
        std::cout << "\nCompleted QMOM workflow. Output disabled by configuration.\n";
    }

    return 0;
}
