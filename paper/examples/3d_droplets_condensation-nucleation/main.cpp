#include "droplet_model.hpp"
#include "monte_carlo/evolution.hpp"
#include "monte_carlo/moments.hpp"
#include "qmomkit/qmomkit.hpp"
#include "sampling.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <filesystem>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";
    const double initial_time = 0.0;
    constexpr double final_time = 1.0;
    constexpr double time_step = 0.01;

    try
    {
        auto output_configured = qmomkit::config::load_output_options(configuration_file);
        if (!output_configured)
        {
            paper_3d_droplets_condensation_nucleation::print_configuration_errors(output_configured.error());
            return 1;
        }
        auto phase_space_evolution_configured = qmomkit::config::load_cqmom_phase_space_evolution_options(configuration_file);
        if (!phase_space_evolution_configured)
        {
            paper_3d_droplets_condensation_nucleation::print_configuration_errors(phase_space_evolution_configured.error());
            return 1;
        }
        const qmomkit::OutputOptions output_options = std::move(output_configured).value();

        const std::optional<std::string> configuration_text = paper_3d_droplets_condensation_nucleation::read_configuration_text(configuration_file);
        if (!configuration_text.has_value())
        {
            return 1;
        }

        auto configured_cqmom = qmomkit::config::parse_cqmom_prepared(*configuration_text, configuration_file.string());
        if (!configured_cqmom)
        {
            paper_3d_droplets_condensation_nucleation::print_configuration_errors(configured_cqmom.error());
            return 1;
        }

        auto initial_condition = qmomkit::config::parse_cqmom_initial_condition(*configuration_text, configuration_file.string());
        if (!initial_condition)
        {
            paper_3d_droplets_condensation_nucleation::print_configuration_errors(initial_condition.error());
            return 1;
        }

        std::optional<qmomkit::monte_carlo::TensorSampledInitialization> sampled =
            paper_3d_droplets_condensation_nucleation::sample_initial_distribution(std::move(initial_condition).value(),
                                                                                   configured_cqmom.value().hierarchy());
        if (!sampled.has_value())
        {
            return 1;
        }
        if (!sampled->ok())
        {
            std::cerr << sampled->sampling.message << sampled->moments.message << '\n';
            return 1;
        }

        // The paper comparison retains the canonical common initialization: these
        // CQMOM moments and the independently evolved MC ensemble come from one sample.
        qmomkit::monte_carlo::Ensemble ensemble = sampled->sampling.ensemble;
        double number_density = sampled->sampling.number_density;
        auto sampled_case_result = qmomkit::config::with_initial_moments(configured_cqmom.value(), std::move(sampled->moments.moments));
        if (!sampled_case_result)
        {
            paper_3d_droplets_condensation_nucleation::print_configuration_errors(sampled_case_result.error());
            return 1;
        }
        const qmomkit::CqmomPreparedCase cqmom_prepared_case = std::move(sampled_case_result).value();

        const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(cqmom_prepared_case);
        if (!cqmom.ok())
        {
            std::cerr << "CQMOM inversion failed: " << cqmom.message << '\n';
            return 1;
        }

        const qmomkit::RawMomentRequest initial_number_request{std::vector<qmomkit::MultiIndex>{{0, 0, 0}}};
        const qmomkit::SelectedRawMomentResult initial_number_moment = qmomkit::reconstruct_raw_moments(cqmom.quadrature, initial_number_request);
        if (!initial_number_moment.ok() || initial_number_moment.values.size() != 1 || !std::isfinite(initial_number_moment.values.front()) ||
            initial_number_moment.values.front() <= 0.0)
        {
            std::cerr << "Could not obtain a positive initial M_000 for the nucleation flux.\n";
            return 1;
        }
        const paper_3d_droplets_condensation_nucleation::NucleationModel nucleation_model =
            paper_3d_droplets_condensation_nucleation::make_nucleation_model(initial_number_moment.values.front());

        qmomkit::CqmomState state = cqmom_prepared_case.initial_state();
        qmomkit::CqmomQuadrature quadrature = cqmom.quadrature;
        const qmomkit::CqmomMomentLayout moment_output_layout = qmomkit::create_cqmom_moment_layout(cqmom_prepared_case.hierarchy());
        const qmomkit::RawMomentRequest moment_output_request{moment_output_layout.required_indices.to_multi_indices()};
        const qmomkit::PhaseSpaceVelocity velocity = paper_3d_droplets_condensation_nucleation::droplet_phase_space_velocity();
        qmomkit::CqmomPhaseSpaceEvolutionOptions phase_space_evolution_options = std::move(phase_space_evolution_configured).value();
        if (phase_space_evolution_options.node_options.scheme != qmomkit::NodeAdvectionScheme::rk23)
        {
            std::cerr << "The paper Monte Carlo reference requires RK23 node advection.\n";
            return 1;
        }
        if (phase_space_evolution_options.scheme != qmomkit::PhaseSpaceEvolutionScheme::strang)
        {
            std::cerr << "The condensing-droplet nucleation example requires the configured Strang splitting scheme.\n";
            return 1;
        }
        const qmomkit::AdaptiveRkOptions mc_advection_options = phase_space_evolution_options.node_options.adaptive_options;

        paper_3d_droplets_condensation_nucleation::validate_nucleation_interval(initial_time, final_time);

        const qmomkit::CqmomMomentSource zero_moment_source = [](double, const qmomkit::CqmomSourceContext& context)
        { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };
        const qmomkit::CqmomMomentBoundaryContribution nucleation_boundary =
            paper_3d_droplets_condensation_nucleation::nucleation_boundary_contribution(nucleation_model);

        paper_3d_droplets_condensation_nucleation::print_cqmom_run_configuration(cqmom_prepared_case, phase_space_evolution_options);

        paper_3d_droplets_condensation_nucleation::ensure_output_directory(output_options.directory);

        const auto estimate_mc_moments = [&](const qmomkit::RawMomentRequest& request)
        {
            const qmomkit::monte_carlo::SelectedEnsembleMomentEstimate estimate =
                qmomkit::monte_carlo::compute_requested_ensemble_moments(ensemble, number_density, request);
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

        std::vector<std::pair<double, qmomkit::CqmomQuadrature>> quadrature_history{{initial_time, quadrature}};
        std::optional<qmomkit::MomentWriter> mc_moment_history_writer;
        qmomkit::SelectedRawMomentResult mc_final_moments = estimate_mc_moments(moment_output_request);
        const qmomkit::SelectedRawMomentResult mc_initial_moments = mc_final_moments;
        if (output_options.write_moments)
        {
            mc_moment_history_writer.emplace(
                qmomkit::MomentWriter::for_cqmom(output_options.directory / "moment_history_mc.txt", moment_output_layout));
            mc_moment_history_writer->write(initial_time, qmomkit::create_cqmom_moment_data(moment_output_layout, mc_final_moments));
        }

        // Each physical birth event needs a distinct deterministic sample.
        std::size_t birth_event_index = 0;
        double time = initial_time;
        while (time < final_time)
        {
            const double current_step =
                paper_3d_droplets_condensation_nucleation::limit_step_to_nucleation_interval(time, std::min(time_step, final_time - time));
            const bool nucleation_active = paper_3d_droplets_condensation_nucleation::nucleation_is_active(time);
            const qmomkit::CqmomPhaseSpaceEvolutionResult phase_space_evolved =
                nucleation_active ? qmomkit::advance_cqmom_phase_space(state, time, current_step, velocity, zero_moment_source, nucleation_boundary,
                                                                       phase_space_evolution_options)
                                  : qmomkit::advance_cqmom_advection(state, time, current_step, velocity, phase_space_evolution_options);
            if (!phase_space_evolved.ok())
            {
                std::cerr << "CQMOM phase-space evolution failed: " << phase_space_evolved.message << '\n';
                return 1;
            }
            state = phase_space_evolved.state;
            quadrature = phase_space_evolved.state.quadrature();
            if (std::abs(phase_space_evolved.final_time - (time + current_step)) > 1.0e-12 * std::max(1.0, std::abs(time + current_step)))
            {
                std::cerr << "CQMOM phase-space evolution did not reach the prescribed step endpoint.\n";
                return 1;
            }

            if (nucleation_active)
            {
                const double source_time = time + 0.5 * current_step;
                qmomkit::monte_carlo::EnsembleEvolutionResult mc_first_half_advection =
                    qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, time, source_time, mc_advection_options);
                if (!mc_first_half_advection.ok())
                {
                    std::cerr << "First-half Monte Carlo advection failed: " << mc_first_half_advection.message << '\n';
                    return 1;
                }
                ensemble = std::move(mc_first_half_advection.ensemble);
                const unsigned long long birth_seed =
                    paper_3d_droplets_condensation_nucleation::initialization_seed + static_cast<unsigned long long>(birth_event_index++);
                const qmomkit::monte_carlo::PopulationBirthResult birth = paper_3d_droplets_condensation_nucleation::inject_nucleated_droplets(
                    nucleation_model, ensemble, number_density, time, current_step, birth_seed);
                if (!birth.ok())
                {
                    std::cerr << "Monte Carlo population injection failed: " << birth.message << '\n';
                    return 1;
                }
                ensemble = birth.ensemble;
                number_density = birth.number_density;
                qmomkit::monte_carlo::EnsembleEvolutionResult mc_second_half_advection =
                    qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, source_time, time + current_step, mc_advection_options);
                if (!mc_second_half_advection.ok())
                {
                    std::cerr << "Second-half Monte Carlo advection failed: " << mc_second_half_advection.message << '\n';
                    return 1;
                }
                ensemble = std::move(mc_second_half_advection.ensemble);
            }
            else
            {
                qmomkit::monte_carlo::EnsembleEvolutionResult mc_advection =
                    qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, time, time + current_step, mc_advection_options);
                if (!mc_advection.ok())
                {
                    std::cerr << "Monte Carlo advection failed: " << mc_advection.message << '\n';
                    return 1;
                }
                ensemble = std::move(mc_advection.ensemble);
            }

            time += current_step;
            quadrature_history.emplace_back(time, quadrature);
            mc_final_moments = estimate_mc_moments(moment_output_request);
            if (mc_moment_history_writer.has_value())
            {
                mc_moment_history_writer->write(time, qmomkit::create_cqmom_moment_data(moment_output_layout, mc_final_moments));
            }
        }

        if (output_options.write_moments)
        {
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
            qmomkit::QuadratureWriter quadrature_history_writer = qmomkit::QuadratureWriter::for_cqmom(
                output_options.directory / "quadrature_history.txt", quadrature_history.front().second.number_of_dimensions());
            for (const auto& [history_time, history_quadrature] : quadrature_history)
            {
                quadrature_history_writer.write(history_time, history_quadrature);
            }
        }

        const qmomkit::SelectedRawMomentResult initial_moments = qmomkit::reconstruct_raw_moments(cqmom.quadrature, moment_output_request);
        const qmomkit::SelectedRawMomentResult final_moments = qmomkit::reconstruct_raw_moments(quadrature, moment_output_request);
        if (!initial_moments.ok() || !final_moments.ok())
        {
            std::cerr << "CQMOM moment reconstruction failed while producing the summary.\n";
            return 1;
        }
        paper_3d_droplets_condensation_nucleation::print_result_summary(initial_moments, final_moments, quadrature, output_options.directory,
                                                                        paper_3d_droplets_condensation_nucleation::initialization_sample_count,
                                                                        paper_3d_droplets_condensation_nucleation::initialization_seed, time_step,
                                                                        final_time, output_options.write_moments, output_options.write_quadrature);
        paper_3d_droplets_condensation_nucleation::print_monte_carlo_result_summary(
            mc_initial_moments, mc_final_moments, output_options.directory, ensemble.size(),
            paper_3d_droplets_condensation_nucleation::initialization_sample_count, paper_3d_droplets_condensation_nucleation::initialization_seed,
            time_step, final_time, output_options.write_moments);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
