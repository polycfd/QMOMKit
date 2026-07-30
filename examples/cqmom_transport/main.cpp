#include "physics.hpp"
#include "qmomkit/qmomkit.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <exception>
#include <filesystem>
#include <iostream>
#include <utility>
#include <vector>

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";
    const double initial_time = 0.0;
    const double final_time = 1.0;
    const double time_step = 0.01;

    try
    {
        auto configured = qmomkit::config::load_cqmom(configuration_file);
        if (!configured)
        {
            cqmom_transport::print_configuration_errors(configured.error());
            return 1;
        }

        auto output_configured = qmomkit::config::load_output_options(configuration_file);
        if (!output_configured)
        {
            cqmom_transport::print_configuration_errors(output_configured.error());
            return 1;
        }
        auto transport_configured = qmomkit::config::load_cqmom_transport_options(configuration_file);
        if (!transport_configured)
        {
            cqmom_transport::print_configuration_errors(transport_configured.error());
            return 1;
        }

        // CqmomCase is the immutable configuration and problem definition.
        const qmomkit::CqmomCase cqmom_case = std::move(configured).value();
        const qmomkit::OutputOptions output_options = std::move(output_configured).value();

        // compute_cqmom(...) computes the initial quadrature for the case.
        const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(cqmom_case);
        if (!cqmom.ok())
        {
            std::cerr << "CQMOM inversion failed: " << cqmom.message << '\n';
            return 1;
        }

        // CqmomState is the runtime state that transport evolves.
        qmomkit::CqmomState state = cqmom_case.initial_state();
        qmomkit::CqmomQuadrature quadrature = cqmom.quadrature;
        const qmomkit::CqmomMomentLayout moment_output_layout = qmomkit::create_cqmom_moment_layout(cqmom_case.hierarchy());
        const qmomkit::RawMomentRequest moment_output_request{moment_output_layout.required_indices};

        // User physics enters through callbacks.
        const qmomkit::PhaseSpaceVelocity velocity = cqmom_transport::continuous_velocity();
        const qmomkit::CqmomMomentSource source = cqmom_transport::binary_fragmentation_source();
        const qmomkit::CqmomMomentBoundaryContribution boundary = cqmom_transport::injection_boundary_contribution();

        // User-facing transport controls live in CqmomTransportOptions.
        qmomkit::CqmomTransportOptions transport_options = std::move(transport_configured).value();
        std::cout << qmomkit::format_cqmom_settings(cqmom_case.moment_storage(), cqmom_case.adaptive_branch_counts(), cqmom_case.vandermonde_solver(),
                                                    cqmom_case.diagnostic_level());
        std::cout << "Transport scheme: ";
        switch (transport_options.scheme)
        {
        case qmomkit::TransportScheme::strang:
            std::cout << "strang";
            break;
        case qmomkit::TransportScheme::cd:
            std::cout << "cd";
            break;
        case qmomkit::TransportScheme::dc:
            std::cout << "dc";
            break;
        }
        std::cout << '\n';
        std::cout << qmomkit::format_node_advection_options(transport_options.node_options);

        std::vector<std::pair<double, qmomkit::CqmomQuadrature>> quadrature_history{{initial_time, quadrature}};

        const double output_interval = output_options.interval.value_or(final_time - initial_time);
        double time = initial_time;
        double next_output_time = initial_time + output_interval;
        while (time < final_time - 1.0e-12)
        {
            const double current_step = std::min(time_step, final_time - time);
            const qmomkit::CqmomTransportResult transported =
                qmomkit::advance_transport_cqmom(state, time, current_step, velocity, source, boundary, transport_options);
            if (!transported.ok())
            {
                std::cerr << "CQMOM transport failed: " << transported.message << '\n';
                return 1;
            }
            if (!transported.quadrature.has_value())
            {
                std::cerr << "CQMOM transport did not return a quadrature.\n";
                return 1;
            }

            state = transported.state;
            quadrature = *transported.quadrature;
            time = transported.final_time;

            if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
            {
                quadrature_history.emplace_back(time, quadrature);
                while (next_output_time <= time + 1.0e-12)
                {
                    next_output_time += output_interval;
                }
            }
        }

        if (output_options.write_quadrature)
        {
            cqmom_transport::ensure_output_directory(output_options.directory);
        }
        if (output_options.write_moments)
        {
            cqmom_transport::ensure_output_directory(output_options.directory);
            qmomkit::MomentWriter::for_cqmom(output_options.directory / "initial_moments.txt", moment_output_request)
                .write(initial_time, cqmom.quadrature);

            qmomkit::MomentWriter moment_history_writer =
                qmomkit::MomentWriter::for_cqmom(output_options.directory / "moment_history.txt", moment_output_request);
            for (const auto& [history_time, history_quadrature] : quadrature_history)
            {
                moment_history_writer.write(history_time, history_quadrature);
            }
        }
        if (output_options.write_quadrature)
        {
            qmomkit::QuadratureWriter::for_cqmom(output_options.directory / "initial_quadrature.txt", cqmom.quadrature.number_of_dimensions())
                .write(initial_time, cqmom.quadrature);

            qmomkit::QuadratureWriter quadrature_history_writer = qmomkit::QuadratureWriter::for_cqmom(
                output_options.directory / "quadrature_history.txt", quadrature_history.front().second.number_of_dimensions());
            for (const auto& [history_time, history_quadrature] : quadrature_history)
            {
                quadrature_history_writer.write(history_time, history_quadrature);
            }
        }
        const qmomkit::SelectedRawMomentResult initial_moments = qmomkit::reconstruct_raw_moments(cqmom.quadrature, moment_output_request);
        if (!initial_moments.ok())
        {
            std::cerr << "Initial CQMOM moment reconstruction failed: " << initial_moments.message << '\n';
            return 1;
        }
        const qmomkit::SelectedRawMomentResult final_moments = qmomkit::reconstruct_raw_moments(quadrature, moment_output_request);
        if (!final_moments.ok())
        {
            std::cerr << "Final CQMOM moment reconstruction failed: " << final_moments.message << '\n';
            return 1;
        }
        cqmom_transport::print_result_summary(initial_moments, final_moments, quadrature, output_options.directory, output_options.write_moments,
                                              output_options.write_quadrature);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
