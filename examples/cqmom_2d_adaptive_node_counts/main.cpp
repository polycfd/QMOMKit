#include "physics.hpp"
#include "qmomkit/qmomkit.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <string>
#include <utility>
#include <vector>

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";
    const double initial_time = 0.0;
    const double final_time = 8.0;
    const double time_step = 0.002;

    try
    {
        auto configured = qmomkit::config::load_cqmom(configuration_file);
        if (!configured)
        {
            cqmom_2d_adaptive_node_counts::print_configuration_errors(configured.error());
            return 1;
        }
        auto output_configured = qmomkit::config::load_output_options(configuration_file);
        if (!output_configured)
        {
            cqmom_2d_adaptive_node_counts::print_configuration_errors(output_configured.error());
            return 1;
        }
        auto phase_space_evolution_configured = qmomkit::config::load_cqmom_phase_space_evolution_options(configuration_file);
        if (!phase_space_evolution_configured)
        {
            cqmom_2d_adaptive_node_counts::print_configuration_errors(phase_space_evolution_configured.error());
            return 1;
        }

        const qmomkit::CqmomCase cqmom_case = std::move(configured).value();
        const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(cqmom_case);
        if (!cqmom.ok())
        {
            std::cerr << "CQMOM inversion failed: " << cqmom.message << '\n';
            return 1;
        }

        qmomkit::CqmomState state = cqmom_case.initial_state();
        qmomkit::CqmomQuadrature quadrature = cqmom.quadrature;
        const qmomkit::CqmomMomentLayout moment_output_layout = qmomkit::create_cqmom_moment_layout(cqmom_case.hierarchy());
        const qmomkit::RawMomentRequest moment_output_request{moment_output_layout.required_indices};
        const cqmom_2d_adaptive_node_counts::VarianceInjectionParameters source_parameters =
            cqmom_2d_adaptive_node_counts::default_variance_injection_parameters();
        const std::size_t velocity_dimension = cqmom_2d_adaptive_node_counts::velocity_phase_space_dimension();
        std::string source_validation_message;
        if (!cqmom_2d_adaptive_node_counts::validate_variance_injection_source_rates(source_parameters, source_validation_message))
        {
            std::cerr << "Variance-injection source validation failed: " << source_validation_message << '\n';
            return 1;
        }
        std::string moment_layout_message;
        if (!cqmom_2d_adaptive_node_counts::validate_required_velocity_moments(moment_output_layout, velocity_dimension, moment_layout_message))
        {
            std::cerr << "Adaptive node-count moment-layout validation failed: " << moment_layout_message << '\n';
            return 1;
        }
        const qmomkit::PhaseSpaceVelocity velocity = cqmom_2d_adaptive_node_counts::advection_velocity(source_parameters);
        const qmomkit::CqmomMomentSource source = cqmom_2d_adaptive_node_counts::variance_injection_source(source_parameters);
        qmomkit::CqmomPhaseSpaceEvolutionOptions phase_space_evolution_options = std::move(phase_space_evolution_configured).value();
        const qmomkit::OutputOptions output_options = std::move(output_configured).value();

        if (output_options.write_moments || output_options.write_quadrature)
        {
            cqmom_2d_adaptive_node_counts::ensure_output_directory(output_options.directory);
        }
        std::optional<qmomkit::MomentWriter> moment_history_writer;
        if (output_options.write_moments)
        {
            qmomkit::MomentWriter::for_cqmom(output_options.directory / "initial_moments.txt", state, output_options)
                .write(initial_time, state, quadrature);
            moment_history_writer.emplace(qmomkit::MomentWriter::for_cqmom(output_options.directory / "moment_history.txt", state, output_options));
            moment_history_writer->write(initial_time, state, quadrature);
        }
        std::vector<cqmom_2d_adaptive_node_counts::HistorySample> step_history{cqmom_2d_adaptive_node_counts::make_history_sample(
            initial_time, cqmom_2d_adaptive_node_counts::source_pulse_active(initial_time, source_parameters), quadrature,
            state.active_branch_counts(), velocity_dimension)};
        std::vector<cqmom_2d_adaptive_node_counts::HistorySample> output_history{step_history.front()};
        std::vector<cqmom_2d_adaptive_node_counts::ActiveCountTransition> active_count_transitions;

        const double output_interval = output_options.interval.value_or(final_time - initial_time);
        double time = initial_time;
        double next_output_time = initial_time + output_interval;
        while (time < final_time - 1.0e-12)
        {
            const double current_step = std::min(time_step, final_time - time);
            const qmomkit::CqmomPhaseSpaceEvolutionResult phase_space_evolved =
                qmomkit::advance_cqmom_phase_space(state, time, current_step, velocity, source, phase_space_evolution_options);
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

            step_history.push_back(
                cqmom_2d_adaptive_node_counts::make_history_sample(time, cqmom_2d_adaptive_node_counts::source_pulse_active(time, source_parameters),
                                                                   quadrature, state.active_branch_counts(), velocity_dimension));
            cqmom_2d_adaptive_node_counts::append_active_count_transitions(
                time, step_history[step_history.size() - 2].quadrature, step_history.back().quadrature,
                phase_space_evolved.diagnostics.cqmom_node_count_transitions, active_count_transitions, velocity_dimension);

            if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
            {
                output_history.push_back(step_history.back());
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
                output_options.directory / "quadrature_history.txt", output_history.front().quadrature.number_of_dimensions());
            for (const cqmom_2d_adaptive_node_counts::HistorySample& sample : output_history)
            {
                quadrature_history_writer.write(sample.time, sample.quadrature);
            }
        }

        const qmomkit::SelectedRawMomentResult initial_moments = qmomkit::reconstruct_raw_moments(cqmom.quadrature, moment_output_request);
        const qmomkit::SelectedRawMomentResult final_moments = qmomkit::reconstruct_raw_moments(quadrature, moment_output_request);
        std::string validation_message;
        if (!cqmom_2d_adaptive_node_counts::validate_history(step_history, step_history.front().moments.velocity_mean,
                                                             source_parameters.relaxation_rate, source_parameters.equilibrium_velocity,
                                                             source_parameters.pulse_start, source_parameters.pulse_end, validation_message))
        {
            std::cerr << "Adaptive node-count validation failed: " << validation_message << '\n';
            return 1;
        }
        std::string lifecycle_message;
        if (!cqmom_2d_adaptive_node_counts::validate_active_count_lifecycle(step_history, active_count_transitions, source_parameters.pulse_start,
                                                                            source_parameters.pulse_end, lifecycle_message))
        {
            std::cerr << "Adaptive node-count lifecycle validation failed: " << lifecycle_message << '\n';
            return 1;
        }
        cqmom_2d_adaptive_node_counts::print_result_summary(
            initial_moments, final_moments, output_history, active_count_transitions, quadrature, output_options.directory,
            cqmom_case.hierarchy().conditioning_order(), cqmom_case.hierarchy().node_counts(), source_parameters.relaxation_rate,
            source_parameters.equilibrium_velocity, source_parameters.pulse_start, source_parameters.pulse_end, source_parameters.source_rate,
            source_parameters.velocity_displacement, time_step, final_time, 2 * cqmom_case.hierarchy().node_count(velocity_dimension) - 1,
            source_validation_message, moment_layout_message, output_options.write_moments, output_options.write_quadrature);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
