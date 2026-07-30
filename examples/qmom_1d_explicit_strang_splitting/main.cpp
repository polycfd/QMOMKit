#include "qmomkit/qmomkit.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>
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

    auto transport_configured = qmomkit::config::load_qmom_transport_options(configuration_file);
    if (!transport_configured)
    {
        std::cerr << transport_configured.error() << '\n';
        return 1;
    }

    const qmomkit::QmomCase qmom_case = std::move(configured).value();
    const qmomkit::OutputOptions output_options = std::move(output_configured).value();
    qmomkit::QmomTransportOptions transport_options = std::move(transport_configured).value();
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
    const qmomkit::PhaseSpaceVelocity velocity = [](double, const std::vector<double>& coordinates)
    {
        constexpr double relaxation_rate = 500.0;
        if (coordinates.size() != 1)
        {
            throw std::invalid_argument("The explicit Strang example expects one phase-space coordinate.");
        }
        return std::vector<double>{-relaxation_rate * coordinates[0]};
    };
    const qmomkit::MomentSource source = [](const double time, const std::vector<double>& current_moments)
    {
        constexpr double source_start = 0.16;
        constexpr double source_end = 0.25;
        constexpr double variance_injection_rate = 12.0;

        std::vector<double> derivative(current_moments.size(), 0.0);
        if (time >= source_start && time < source_end && derivative.size() > 2)
        {
            derivative[2] = variance_injection_rate;
        }
        return derivative;
    };

    std::optional<qmomkit::QuadratureWriter> quadrature_history_writer;
    try
    {
        if (output_options.write_quadrature)
        {
            qmom_1d_explicit_strang_splitting::ensure_output_directory(output_options.directory);
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
        qmomkit::QmomState trial_state = state;

        // First continuous half-step.
        const qmomkit::QmomTransportResult first_continuous =
            qmomkit::advance_qmom_continuous(trial_state, time, half_step, velocity, transport_options);
        if (!first_continuous.ok())
        {
            std::cerr << "QMOM first continuous stage failed: " << first_continuous.message << '\n';
            return 1;
        }
        trial_state = first_continuous.state;

        // Full discontinuous source step.
        const std::size_t before_discontinuous_nodes = trial_state.active_node_count();
        const qmomkit::QmomTransportResult discontinuous =
            qmomkit::advance_qmom_discontinuous(trial_state, maximum_node_count, time + half_step, step, source, transport_options);
        if (!discontinuous.ok())
        {
            std::cerr << "QMOM discontinuous stage failed: " << discontinuous.message << '\n';
            return 1;
        }
        trial_state = discontinuous.state;
        const std::size_t after_discontinuous_nodes = trial_state.active_node_count();

        // Second continuous half-step.
        const qmomkit::QmomTransportResult second_continuous =
            qmomkit::advance_qmom_continuous(trial_state, time + half_step, half_step, velocity, transport_options);
        if (!second_continuous.ok())
        {
            std::cerr << "QMOM second continuous stage failed: " << second_continuous.message << '\n';
            return 1;
        }
        trial_state = second_continuous.state;

        state = trial_state;
        time += step;

        if (before_discontinuous_nodes != after_discontinuous_nodes)
        {
            std::cout << "t = " << time - half_step << ": discontinuous reconstruction changed active nodes " << before_discontinuous_nodes << " -> "
                      << after_discontinuous_nodes << '\n';
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
