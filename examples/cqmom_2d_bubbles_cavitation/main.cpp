#include "physics.hpp"
#include "qmomkit/qmomkit.hpp"
#include "sampling.hpp"
#include "utilities.hpp"

#include <algorithm>
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
    const double final_time = 13.9;
    const double time_step = 0.01;

    try
    {
        // A .qmom file can hold the case, the evolution controls and the output
        // settings together, but they are loaded separately: the case loaders
        // validate every recognised key yet return only the case, so anything
        // the driver needs must be asked for explicitly.
        auto output_configured = qmomkit::config::load_output_options(configuration_file);
        if (!output_configured)
        {
            cqmom_2d_bubbles_cavitation::print_configuration_errors(output_configured.error());
            return 1;
        }
        auto phase_space_evolution_configured = qmomkit::config::load_cqmom_phase_space_evolution_options(configuration_file);
        if (!phase_space_evolution_configured)
        {
            cqmom_2d_bubbles_cavitation::print_configuration_errors(phase_space_evolution_configured.error());
            return 1;
        }
        const qmomkit::OutputOptions output_options = std::move(output_configured).value();

        // The initial condition is parsed separately from the case because this
        // example does not use it directly: it samples from it instead, below.
        const std::optional<std::string> configuration_text = cqmom_2d_bubbles_cavitation::read_configuration_text(configuration_file);
        if (!configuration_text.has_value())
            return 1;

        auto configured_cqmom = qmomkit::config::parse_cqmom_prepared(*configuration_text, configuration_file.string());
        if (!configured_cqmom)
        {
            cqmom_2d_bubbles_cavitation::print_configuration_errors(configured_cqmom.error());
            return 1;
        }

        auto initial_condition = qmomkit::config::parse_cqmom_initial_condition(*configuration_text, configuration_file.string());
        if (!initial_condition)
        {
            cqmom_2d_bubbles_cavitation::print_configuration_errors(initial_condition.error());
            return 1;
        }

        // Draw a sample ensemble from the configured distribution and compute its
        // raw moments. Starting the deterministic run from sampled moments is what
        // makes it directly comparable with the Monte Carlo companion.
        const std::optional<qmomkit::monte_carlo::TensorSampledInitialization> sampled =
            cqmom_2d_bubbles_cavitation::sample_initial_distribution(std::move(initial_condition).value(), configured_cqmom.value().hierarchy());
        if (!sampled.has_value())
            return 1;
        if (!sampled->ok())
        {
            std::cerr << sampled->sampling.message << sampled->moments.message << '\n';
            return 1;
        }

        // Swap the sampled moments into the loaded case. This keeps the configured
        // hierarchy, supports and reconstruction settings and replaces only the
        // initial values, which avoids rewriting the configuration file.
        auto sampled_case_result = qmomkit::config::with_initial_moments(configured_cqmom.value(), std::move(sampled->moments.moments));
        if (!sampled_case_result)
        {
            cqmom_2d_bubbles_cavitation::print_configuration_errors(sampled_case_result.error());
            return 1;
        }
        const qmomkit::CqmomPreparedCase cqmom_case = std::move(sampled_case_result).value();

        // The case is the fixed setup; the state is what advances in time.
        qmomkit::CqmomState state = cqmom_case.initial_state();
        const qmomkit::CqmomQuadrature initial_quadrature = state.quadrature();
        const qmomkit::CqmomMomentLayout moment_output_layout = qmomkit::create_cqmom_moment_layout(cqmom_case.hierarchy());
        const qmomkit::RawMomentRequest moment_output_request{moment_output_layout.required_indices.to_multi_indices()};
        const qmomkit::PhaseSpaceVelocity velocity = cqmom_2d_bubbles_cavitation::advection_velocity();

        qmomkit::CqmomPhaseSpaceEvolutionOptions phase_space_evolution_options = std::move(phase_space_evolution_configured).value();
        std::cout << qmomkit::format_node_advection_options(phase_space_evolution_options.node_options);

        if (output_options.write_moments || output_options.write_quadrature)
        {
            cqmom_2d_bubbles_cavitation::ensure_output_directory(output_options.directory);
        }
        std::optional<qmomkit::MomentWriter> moment_history_writer;
        if (output_options.write_moments)
        {
            qmomkit::MomentWriter::for_cqmom(output_options.directory / "initial_moments.txt", state, output_options)
                .write(initial_time, state, state.quadrature());
            moment_history_writer.emplace(qmomkit::MomentWriter::for_cqmom(output_options.directory / "moment_history.txt", state, output_options));
            moment_history_writer->write(initial_time, state, state.quadrature());
        }
        std::optional<qmomkit::QuadratureWriter> quadrature_history_writer;
        if (output_options.write_quadrature)
        {
            qmomkit::QuadratureWriter::for_cqmom(output_options.directory / "initial_quadrature.txt", initial_quadrature.number_of_dimensions())
                .write(initial_time, initial_quadrature);
            quadrature_history_writer.emplace(
                qmomkit::QuadratureWriter::for_cqmom(output_options.directory / "quadrature_history.txt", initial_quadrature.number_of_dimensions()));
            quadrature_history_writer->write(initial_time, initial_quadrature);
        }

        // This example models advection only, so the time loop calls
        // advance_cqmom_advection(...) rather than the full
        // advance_cqmom_phase_space(...) step. There is no moment source to
        // compose with, so there is nothing for a splitting scheme to do.
        const double output_interval = output_options.interval.value_or(final_time - initial_time);
        double time = initial_time;
        double next_output_time = initial_time + output_interval;
        while (time < final_time - 1.0e-12)
        {
            const double current_step = std::min(time_step, final_time - time);
            const qmomkit::CqmomPhaseSpaceEvolutionResult phase_space_evolved =
                qmomkit::advance_cqmom_advection(state, time, current_step, velocity, phase_space_evolution_options);
            if (!phase_space_evolved.ok())
            {
                std::cerr << "CQMOM phase-space evolution failed: " << phase_space_evolved.message << '\n';
                return 1;
            }
            state = phase_space_evolved.state;
            time = phase_space_evolved.final_time;

            if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
            {
                if (quadrature_history_writer.has_value())
                {
                    quadrature_history_writer->write(time, state.quadrature());
                }
                if (moment_history_writer.has_value())
                {
                    moment_history_writer->write(time, state, state.quadrature());
                }
                while (next_output_time <= time + 1.0e-12)
                {
                    next_output_time += output_interval;
                }
            }
        }

        // Reconstruct the same moment set from the initial and final quadratures
        // so the summary compares like with like.
        const qmomkit::SelectedRawMomentResult initial_moments = qmomkit::reconstruct_raw_moments(initial_quadrature, moment_output_request);
        if (!initial_moments.ok())
        {
            std::cerr << "Initial CQMOM moment reconstruction failed: " << initial_moments.message << '\n';
            return 1;
        }
        const qmomkit::SelectedRawMomentResult final_moments = qmomkit::reconstruct_raw_moments(state.quadrature(), moment_output_request);
        if (!final_moments.ok())
        {
            std::cerr << "Final CQMOM moment reconstruction failed: " << final_moments.message << '\n';
            return 1;
        }
        cqmom_2d_bubbles_cavitation::print_result_summary(
            initial_moments, final_moments, state.quadrature(), output_options.directory, cqmom_2d_bubbles_cavitation::initialization_sample_count,
            cqmom_2d_bubbles_cavitation::initialization_seed, time_step, final_time, output_options.write_moments, output_options.write_quadrature);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
