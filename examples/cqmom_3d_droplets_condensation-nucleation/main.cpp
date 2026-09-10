#include "droplet_model.hpp"
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
    const double final_time = 1.0;
    const double time_step = 0.01;

    try
    {
        auto output_configured = qmomkit::config::load_output_options(configuration_file);
        if (!output_configured)
        {
            cqmom_3d_droplets_condensation_nucleation::print_configuration_errors(output_configured.error());
            return 1;
        }
        auto phase_space_evolution_configured = qmomkit::config::load_cqmom_phase_space_evolution_options(configuration_file);
        if (!phase_space_evolution_configured)
        {
            cqmom_3d_droplets_condensation_nucleation::print_configuration_errors(phase_space_evolution_configured.error());
            return 1;
        }
        const qmomkit::OutputOptions output_options = std::move(output_configured).value();

        const std::optional<std::string> configuration_text = cqmom_3d_droplets_condensation_nucleation::read_configuration_text(configuration_file);
        if (!configuration_text.has_value())
            return 1;

        auto configured_cqmom = qmomkit::config::parse_cqmom(*configuration_text, configuration_file.string());
        if (!configured_cqmom)
        {
            cqmom_3d_droplets_condensation_nucleation::print_configuration_errors(configured_cqmom.error());
            return 1;
        }

        auto initial_condition = qmomkit::config::parse_cqmom_initial_condition(*configuration_text, configuration_file.string());
        if (!initial_condition)
        {
            cqmom_3d_droplets_condensation_nucleation::print_configuration_errors(initial_condition.error());
            return 1;
        }

        const std::optional<qmomkit::monte_carlo::TensorSampledInitialization> sampled =
            cqmom_3d_droplets_condensation_nucleation::sample_initial_distribution(std::move(initial_condition).value(),
                                                                                   configured_cqmom.value().hierarchy());
        if (!sampled.has_value())
            return 1;
        if (!sampled->ok())
        {
            std::cerr << sampled->sampling.message << sampled->moments.message << '\n';
            return 1;
        }

        auto sampled_case_result = qmomkit::config::with_initial_moments(configured_cqmom.value(), std::move(sampled->moments.moments));
        if (!sampled_case_result)
        {
            cqmom_3d_droplets_condensation_nucleation::print_configuration_errors(sampled_case_result.error());
            return 1;
        }
        const qmomkit::CqmomCase cqmom_case = std::move(sampled_case_result).value();

        const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(cqmom_case);
        if (!cqmom.ok())
        {
            std::cerr << "CQMOM inversion failed: " << cqmom.message << '\n';
            return 1;
        }

        const qmomkit::RawMomentRequest initial_number_request{std::vector<qmomkit::MultiIndex>{{0, 0, 0}}};
        const qmomkit::SelectedRawMomentResult initial_number_moment = qmomkit::reconstruct_raw_moments(cqmom.quadrature, initial_number_request);
        if (!initial_number_moment.ok() || initial_number_moment.values.size() != 1)
        {
            std::cerr << "Could not obtain initial M[0,0,0] for the nucleation model.\n";
            return 1;
        }
        const cqmom_3d_droplets_condensation_nucleation::NucleationModel nucleation_model =
            cqmom_3d_droplets_condensation_nucleation::make_nucleation_model(initial_number_moment.values.front());

        qmomkit::CqmomState state = cqmom_case.initial_state();
        qmomkit::CqmomQuadrature quadrature = cqmom.quadrature;
        const qmomkit::CqmomMomentLayout moment_output_layout = qmomkit::create_cqmom_moment_layout(cqmom_case.hierarchy());
        const qmomkit::RawMomentRequest moment_output_request{moment_output_layout.required_indices};
        const qmomkit::PhaseSpaceVelocity velocity = cqmom_3d_droplets_condensation_nucleation::droplet_phase_space_velocity();

        qmomkit::CqmomPhaseSpaceEvolutionOptions phase_space_evolution_options = std::move(phase_space_evolution_configured).value();
        if (phase_space_evolution_options.scheme != qmomkit::PhaseSpaceEvolutionScheme::strang)
        {
            std::cerr << "The condensation-nucleation example requires configured Strang splitting.\n";
            return 1;
        }
        cqmom_3d_droplets_condensation_nucleation::validate_nucleation_interval(initial_time, final_time);
        const qmomkit::CqmomMomentSource zero_moment_source = [](double, const qmomkit::CqmomSourceContext& context)
        { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };
        const auto nucleation_boundary = cqmom_3d_droplets_condensation_nucleation::nucleation_boundary_contribution(nucleation_model);
        cqmom_3d_droplets_condensation_nucleation::print_cqmom_run_configuration(cqmom_case, phase_space_evolution_options);

        if (output_options.write_moments || output_options.write_quadrature)
        {
            cqmom_3d_droplets_condensation_nucleation::ensure_output_directory(output_options.directory);
        }
        std::optional<qmomkit::MomentWriter> moment_history_writer;
        if (output_options.write_moments)
        {
            qmomkit::MomentWriter::for_cqmom(output_options.directory / "initial_moments.txt", state).write(initial_time, state);
            moment_history_writer.emplace(qmomkit::MomentWriter::for_cqmom(output_options.directory / "moment_history.txt", state));
            moment_history_writer->write(initial_time, state);
        }
        std::vector<std::pair<double, qmomkit::CqmomQuadrature>> quadrature_history{{initial_time, quadrature}};

        const double output_interval = output_options.interval.value_or(final_time - initial_time);
        double time = initial_time;
        double next_output_time = initial_time + output_interval;
        while (time < final_time - 1.0e-12)
        {
            const double current_step =
                cqmom_3d_droplets_condensation_nucleation::limit_step_to_nucleation_interval(time, std::min(time_step, final_time - time));
            const bool nucleation_active = cqmom_3d_droplets_condensation_nucleation::nucleation_is_active(time);
            const qmomkit::CqmomPhaseSpaceEvolutionResult phase_space_evolved =
                nucleation_active ? qmomkit::advance_cqmom_phase_space(state, time, current_step, velocity, zero_moment_source, nucleation_boundary,
                                                                       phase_space_evolution_options)
                                  : qmomkit::advance_cqmom_advection(state, time, current_step, velocity, phase_space_evolution_options);
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
                    moment_history_writer->write(time, state);
                }
                while (next_output_time <= time + 1.0e-12)
                {
                    next_output_time += output_interval;
                }
            }
        }

        if (output_options.write_quadrature)
        {
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
        cqmom_3d_droplets_condensation_nucleation::print_result_summary(initial_moments, final_moments, quadrature, output_options.directory,
                                                                        cqmom_3d_droplets_condensation_nucleation::initialization_sample_count,
                                                                        cqmom_3d_droplets_condensation_nucleation::initialization_seed, time_step,
                                                                        final_time, output_options.write_moments, output_options.write_quadrature);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
