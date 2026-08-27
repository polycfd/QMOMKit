#include "droplet_model.hpp"
#include "monte_carlo/evolution.hpp"
#include "monte_carlo/moments.hpp"
#include "qmomkit/configuration.hpp"
#include "qmomkit/cqmom.hpp"
#include "qmomkit/output.hpp"
#include "qmomkit/qmomkit.hpp"
#include "sampling.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <utility>

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";
    const double initial_time = 0.0;
    const double final_time = 1.0;
    const double time_step = 0.01;

    try
    {
        const std::optional<std::string> configuration_text = cqmom_3d_condensing_nucleating_droplets::read_configuration_text(configuration_file);
        if (!configuration_text.has_value())
        {
            return 1;
        }

        auto configured_cqmom = qmomkit::config::parse_cqmom(*configuration_text, configuration_file.string());
        if (!configured_cqmom)
        {
            cqmom_3d_condensing_nucleating_droplets::print_configuration_errors(configured_cqmom.error());
            return 1;
        }

        auto initial_condition = qmomkit::config::parse_cqmom_initial_condition(*configuration_text, configuration_file.string());
        if (!initial_condition)
        {
            cqmom_3d_condensing_nucleating_droplets::print_configuration_errors(initial_condition.error());
            return 1;
        }

        auto output_configured = qmomkit::config::load_output_options(configuration_file);
        if (!output_configured)
        {
            cqmom_3d_condensing_nucleating_droplets::print_configuration_errors(output_configured.error());
            return 1;
        }
        auto phase_space_evolution_configured = qmomkit::config::load_cqmom_phase_space_evolution_options(configuration_file);
        if (!phase_space_evolution_configured)
        {
            cqmom_3d_condensing_nucleating_droplets::print_configuration_errors(phase_space_evolution_configured.error());
            return 1;
        }
        const qmomkit::CqmomCase cqmom_case = std::move(configured_cqmom).value();
        const qmomkit::OutputOptions output_options = std::move(output_configured).value();
        const qmomkit::CqmomPhaseSpaceEvolutionOptions phase_space_evolution_options = std::move(phase_space_evolution_configured).value();
        const qmomkit::CqmomMomentLayout moment_layout = qmomkit::create_cqmom_moment_layout(cqmom_case.hierarchy());
        const qmomkit::RawMomentRequest moment_request{moment_layout.required_indices};

        const std::optional<qmomkit::monte_carlo::TensorSampledInitialization> sampled =
            cqmom_3d_condensing_nucleating_droplets::sample_initial_distribution(std::move(initial_condition).value(), cqmom_case.hierarchy());
        if (!sampled.has_value())
        {
            return 1;
        }
        if (!sampled->ok())
        {
            std::cerr << sampled->sampling.message << sampled->moments.message << '\n';
            return 1;
        }

        qmomkit::monte_carlo::Ensemble ensemble = sampled->sampling.ensemble;
        double number_density = sampled->sampling.number_density;
        const cqmom_3d_condensing_nucleating_droplets::NucleationModel nucleation_model =
            cqmom_3d_condensing_nucleating_droplets::make_nucleation_model(number_density);
        const auto velocity = cqmom_3d_condensing_nucleating_droplets::droplet_phase_space_velocity();
        const qmomkit::AdaptiveRkOptions advection_options = phase_space_evolution_options.node_options.adaptive_options;
        if (phase_space_evolution_options.scheme != qmomkit::PhaseSpaceEvolutionScheme::strang)
        {
            std::cerr << "The condensation-nucleation example requires configured Strang splitting.\n";
            return 1;
        }
        cqmom_3d_condensing_nucleating_droplets::validate_nucleation_interval(initial_time, final_time);

        cqmom_3d_condensing_nucleating_droplets::print_monte_carlo_run_configuration(cqmom_case, phase_space_evolution_options.node_options);

        const auto adapt_monte_carlo_result = [](const qmomkit::monte_carlo::SelectedEnsembleMomentEstimate& estimate)
        {
            if (!estimate.ok())
            {
                throw std::runtime_error("Monte Carlo moment extraction failed: " + estimate.message);
            }
            return qmomkit::SelectedRawMomentResult{estimate.indices,
                                                    estimate.values,
                                                    qmomkit::RawMomentStatus::success,
                                                    "",
                                                    qmomkit::invalid_component_index,
                                                    qmomkit::invalid_dimension_index,
                                                    {}};
        };

        qmomkit::SelectedRawMomentResult final_moments =
            adapt_monte_carlo_result(qmomkit::monte_carlo::compute_requested_ensemble_moments(ensemble, number_density, moment_request));
        const qmomkit::SelectedRawMomentResult initial_moments = final_moments;

        std::optional<qmomkit::MomentWriter> moment_history_writer;
        if (output_options.write_moments)
        {
            cqmom_3d_condensing_nucleating_droplets::ensure_output_directory(output_options.directory);
            moment_history_writer.emplace(qmomkit::MomentWriter::for_cqmom(output_options.directory / "moment_history_mc.txt", moment_layout));
            moment_history_writer->write(initial_time, qmomkit::create_cqmom_moment_data(moment_layout, final_moments));
        }

        const double output_interval = output_options.interval.value_or(final_time - initial_time);
        double time = initial_time;
        double next_output_time = initial_time + output_interval;
        std::size_t birth_event_index = 0;
        while (time < final_time - 1.0e-12)
        {
            const double current_step =
                cqmom_3d_condensing_nucleating_droplets::limit_step_to_nucleation_interval(time, std::min(time_step, final_time - time));
            const bool nucleation_active = cqmom_3d_condensing_nucleating_droplets::nucleation_is_active(time);
            const double source_time = time + 0.5 * current_step;
            if (nucleation_active)
            {
                qmomkit::monte_carlo::EnsembleEvolutionResult first_advection =
                    qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, time, source_time, advection_options);
                if (!first_advection.ok())
                {
                    std::cerr << "First advection half-step failed: " << first_advection.message << '\n';
                    return 1;
                }
                const qmomkit::monte_carlo::PopulationBirthResult birth = cqmom_3d_condensing_nucleating_droplets::inject_nucleated_droplets(
                    nucleation_model, first_advection.ensemble, number_density, time, current_step,
                    cqmom_3d_condensing_nucleating_droplets::initialization_seed + static_cast<unsigned long long>(birth_event_index++));
                if (!birth.ok())
                {
                    std::cerr << "Monte Carlo population injection failed: " << birth.message << '\n';
                    return 1;
                }
                number_density = birth.number_density;
                qmomkit::monte_carlo::EnsembleEvolutionResult second_advection =
                    qmomkit::monte_carlo::evolve_ensemble(birth.ensemble, velocity, source_time, time + current_step, advection_options);
                if (!second_advection.ok())
                {
                    std::cerr << "Second advection half-step failed: " << second_advection.message << '\n';
                    return 1;
                }
                ensemble = std::move(second_advection.ensemble);
            }
            else
            {
                qmomkit::monte_carlo::EnsembleEvolutionResult advection =
                    qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, time, time + current_step, advection_options);
                if (!advection.ok())
                {
                    std::cerr << "Advection particle step failed: " << advection.message << '\n';
                    return 1;
                }
                ensemble = std::move(advection.ensemble);
            }
            time += current_step;

            if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
            {
                final_moments =
                    adapt_monte_carlo_result(qmomkit::monte_carlo::compute_requested_ensemble_moments(ensemble, number_density, moment_request));
                if (moment_history_writer.has_value())
                {
                    moment_history_writer->write(time, qmomkit::create_cqmom_moment_data(moment_layout, final_moments));
                }
                while (next_output_time <= time + 1.0e-12)
                {
                    next_output_time += output_interval;
                }
            }
        }

        if (!moment_history_writer.has_value())
        {
            final_moments =
                adapt_monte_carlo_result(qmomkit::monte_carlo::compute_requested_ensemble_moments(ensemble, number_density, moment_request));
        }
        cqmom_3d_condensing_nucleating_droplets::print_monte_carlo_result_summary(
            initial_moments, final_moments, output_options.directory, ensemble.size(),
            cqmom_3d_condensing_nucleating_droplets::initialization_sample_count, cqmom_3d_condensing_nucleating_droplets::initialization_seed,
            time_step, final_time, output_options.write_moments);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
