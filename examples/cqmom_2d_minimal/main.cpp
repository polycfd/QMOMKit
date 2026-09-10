#include "physics.hpp"
#include "qmomkit/qmomkit.hpp"
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

    try
    {
        auto configured = qmomkit::config::load_cqmom(configuration_file);
        if (!configured)
        {
            cqmom_2d_minimal::print_configuration_errors(configured.error());
            return 1;
        }

        auto output_configured = qmomkit::config::load_output_options(configuration_file);
        if (!output_configured)
        {
            cqmom_2d_minimal::print_configuration_errors(output_configured.error());
            return 1;
        }
        auto phase_space_evolution_configured = qmomkit::config::load_cqmom_phase_space_evolution_options(configuration_file);
        if (!phase_space_evolution_configured)
        {
            cqmom_2d_minimal::print_configuration_errors(phase_space_evolution_configured.error());
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

        // CqmomState is the runtime state advanced by phase-space evolution.
        qmomkit::CqmomState state = cqmom_case.initial_state();
        qmomkit::CqmomQuadrature quadrature = cqmom.quadrature;
        const qmomkit::CqmomMomentLayout moment_output_layout = qmomkit::create_cqmom_moment_layout(cqmom_case.hierarchy());
        const qmomkit::RawMomentRequest moment_output_request{moment_output_layout.required_indices};

        // User physics enters through callbacks.
        const qmomkit::PhaseSpaceVelocity velocity = cqmom_2d_minimal::advection_velocity();
        const qmomkit::CqmomMomentSource source = cqmom_2d_minimal::binary_fragmentation_source();
        const qmomkit::CqmomMomentBoundaryContribution boundary = cqmom_2d_minimal::injection_boundary_contribution();

        // User-facing phase-space evolution controls live in CqmomPhaseSpaceEvolutionOptions.
        qmomkit::CqmomPhaseSpaceEvolutionOptions phase_space_evolution_options = std::move(phase_space_evolution_configured).value();
        std::cout << qmomkit::format_cqmom_settings(cqmom_case.moment_storage(), cqmom_case.adaptive_branch_counts(), cqmom_case.vandermonde_solver(),
                                                    cqmom_case.diagnostic_level());
        std::cout << "Phase-space evolution scheme: ";
        switch (phase_space_evolution_options.scheme)
        {
        case qmomkit::PhaseSpaceEvolutionScheme::strang:
            std::cout << "strang";
            break;
        case qmomkit::PhaseSpaceEvolutionScheme::advsrc:
            std::cout << "advsrc";
            break;
        case qmomkit::PhaseSpaceEvolutionScheme::srcadv:
            std::cout << "srcadv";
            break;
        }
        std::cout << '\n';
        std::cout << qmomkit::format_node_advection_options(phase_space_evolution_options.node_options);

        if (output_options.write_moments || output_options.write_quadrature)
        {
            cqmom_2d_minimal::ensure_output_directory(output_options.directory);
        }
        std::optional<qmomkit::MomentWriter> moment_history_writer;
        if (output_options.write_moments)
        {
            qmomkit::MomentWriter::for_cqmom(output_options.directory / "initial_moments.txt", state, output_options)
                .write(initial_time, state, quadrature);
            moment_history_writer.emplace(qmomkit::MomentWriter::for_cqmom(output_options.directory / "moment_history.txt", state, output_options));
            moment_history_writer->write(initial_time, state, quadrature);
        }
        std::vector<std::pair<double, qmomkit::CqmomQuadrature>> quadrature_history{{initial_time, quadrature}};

        const double output_interval = output_options.interval.value_or(final_time - initial_time);
        double time = initial_time;
        double next_output_time = initial_time + output_interval;
        while (time < final_time - 1.0e-12)
        {
            const double current_step = std::min(time_step, final_time - time);
            const qmomkit::CqmomPhaseSpaceEvolutionResult phase_space_evolved =
                qmomkit::advance_cqmom_phase_space(state, time, current_step, velocity, source, boundary, phase_space_evolution_options);
            if (!phase_space_evolved.ok())
            {
                std::cerr << "CQMOM phase-space evolution failed: " << phase_space_evolved.message << '\n';
                return 1;
            }
            if (!phase_space_evolved.quadrature.has_value())
            {
                std::cerr << "CQMOM phase-space evolution did not return a quadrature.\n";
                return 1;
            }

            state = phase_space_evolved.state;
            quadrature = *phase_space_evolved.quadrature;
            time = phase_space_evolved.final_time;

            if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
            {
                quadrature_history.emplace_back(time, quadrature);
                if (moment_history_writer.has_value())
                {
                    moment_history_writer->write(time, state, quadrature);
                }
                while (next_output_time <= time + 1.0e-12)
                {
                    next_output_time += output_interval;
                }
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
        cqmom_2d_minimal::print_result_summary(initial_moments, final_moments, quadrature, output_options.directory, output_options.write_moments,
                                               output_options.write_quadrature);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
