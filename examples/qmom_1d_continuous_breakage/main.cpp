#include "physics.hpp"
#include "qmomkit/qmomkit.hpp"
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

    auto configured = qmomkit::config::load_qmom(configuration_file);
    if (!configured)
    {
        qmom_1d_continuous_breakage::print_configuration_errors(configured.error());
        return 1;
    }

    auto output_configured = qmomkit::config::load_output_options(configuration_file);
    if (!output_configured)
    {
        qmom_1d_continuous_breakage::print_configuration_errors(output_configured.error());
        return 1;
    }
    auto transport_configured = qmomkit::config::load_qmom_transport_options(configuration_file);
    if (!transport_configured)
    {
        qmom_1d_continuous_breakage::print_configuration_errors(transport_configured.error());
        return 1;
    }

    const qmomkit::QmomCase qmom_case = std::move(configured).value();
    const qmomkit::OutputOptions output_options = std::move(output_configured).value();
    const std::vector<double> initial_moments = qmom_case.moments.values();
    std::vector<double> moments = initial_moments;
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
    moments = state.moments();
    const qmomkit::QmomResult initial_quadrature = state.quadrature();

    qmomkit::Quadrature1D quadrature{initial_quadrature.abscissas, initial_quadrature.weights};
    const qmomkit::PhaseSpaceVelocity velocity = qmom_1d_continuous_breakage::continuous_velocity();
    const qmomkit::MomentSource source = qmom_1d_continuous_breakage::binary_fragmentation_source();
    qmomkit::QmomTransportOptions transport_options = std::move(transport_configured).value();

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
            moment_writer.emplace(qmomkit::MomentWriter::for_qmom(output_options.directory / "moments.txt", moments.size()));
            moment_writer->write(initial_time, qmomkit::MomentSet(moments));
        }
        if (output_options.write_quadrature)
        {
            quadrature_writer.emplace(qmomkit::QuadratureWriter::for_qmom(output_options.directory / "quadrature.txt", quadrature.weights.size()));
            quadrature_writer->write(initial_time, quadrature);
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
        const qmomkit::QmomTransportResult transported =
            qmomkit::advance_transport_qmom(state, qmom_case.node_count, time, current_step, velocity, source, transport_options);
        if (!transported.ok())
        {
            std::cerr << "QMOM transport failed: " << transported.message << '\n';
            return 1;
        }

        state = transported.state;
        moments = state.moments();
        quadrature = qmomkit::Quadrature1D{state.quadrature().abscissas, state.quadrature().weights, qmomkit::QuadratureStatus::success, ""};
        time = transported.final_time;

        if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
        {
            try
            {
                if (moment_writer.has_value())
                {
                    moment_writer->write(time, qmomkit::MomentSet(moments));
                }
                if (quadrature_writer.has_value())
                {
                    quadrature_writer->write(time, quadrature);
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
    std::cout << "1D QMOM continuous + discontinuous evolution example\n";

    std::cout << '\n';
    std::cout << "Initial moments:\n";
    for (std::size_t order = 0; order < initial_moments.size(); ++order)
    {
        std::cout << "M" << order << " = " << initial_moments[order] << '\n';
    }

    std::cout << '\n';
    std::cout << "Final moments:\n";
    for (std::size_t order = 0; order < moments.size(); ++order)
    {
        std::cout << "M" << order << " = " << moments[order] << '\n';
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
    for (std::size_t node = 0; node < quadrature.weights.size(); ++node)
    {
        std::cout << node << "     " << quadrature.weights[node] << "   " << quadrature.abscissas[node] << '\n';
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
