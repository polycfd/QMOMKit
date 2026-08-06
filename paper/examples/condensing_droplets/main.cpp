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
#include <fstream>
#include <iomanip>
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
    const double final_time = 1.0;
    const double time_step = 0.01;

    try
    {
        auto output_configured = qmomkit::config::load_output_options(configuration_file);
        if (!output_configured)
        {
            cqmom_condensing_droplets::print_configuration_errors(output_configured.error());
            return 1;
        }
        auto phase_space_evolution_configured = qmomkit::config::load_cqmom_phase_space_evolution_options(configuration_file);
        if (!phase_space_evolution_configured)
        {
            cqmom_condensing_droplets::print_configuration_errors(phase_space_evolution_configured.error());
            return 1;
        }
        const qmomkit::OutputOptions output_options = std::move(output_configured).value();

        const std::optional<std::string> configuration_text = cqmom_condensing_droplets::read_configuration_text(configuration_file);
        if (!configuration_text.has_value())
        {
            return 1;
        }

        auto configured_cqmom = qmomkit::config::parse_cqmom(*configuration_text, configuration_file.string());
        if (!configured_cqmom)
        {
            cqmom_condensing_droplets::print_configuration_errors(configured_cqmom.error());
            return 1;
        }

        auto initial_condition = qmomkit::config::parse_cqmom_initial_condition(*configuration_text, configuration_file.string());
        if (!initial_condition)
        {
            cqmom_condensing_droplets::print_configuration_errors(initial_condition.error());
            return 1;
        }

        std::optional<qmomkit::monte_carlo::TensorSampledInitialization> sampled =
            cqmom_condensing_droplets::sample_initial_distribution(std::move(initial_condition).value(), configured_cqmom.value().hierarchy());
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
        const double number_density = sampled->sampling.number_density;
        auto sampled_case_result = qmomkit::config::with_initial_moments(configured_cqmom.value(), std::move(sampled->moments.moments));
        if (!sampled_case_result)
        {
            cqmom_condensing_droplets::print_configuration_errors(sampled_case_result.error());
            return 1;
        }
        const qmomkit::CqmomCase cqmom_case = std::move(sampled_case_result).value();

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
        const qmomkit::PhaseSpaceVelocity velocity = cqmom_condensing_droplets::droplet_phase_space_velocity();
        qmomkit::CqmomPhaseSpaceEvolutionOptions phase_space_evolution_options = std::move(phase_space_evolution_configured).value();
        if (phase_space_evolution_options.node_options.scheme != qmomkit::NodeAdvectionScheme::rk23)
        {
            std::cerr << "The paper Monte Carlo reference requires RK23 node advection.\n";
            return 1;
        }
        const qmomkit::AdaptiveRkOptions mc_advection_options = phase_space_evolution_options.node_options.adaptive_options;

        cqmom_condensing_droplets::print_cqmom_run_configuration(cqmom_case, phase_space_evolution_options);

        // The configured phase-space order is (r, w, T). The requested
        // Sauter-mean quantity is M_300/M_200 without a factor of two.
        const std::array<qmomkit::MultiIndex, 5> comparison_indices{{{0, 0, 0}, {3, 0, 0}, {2, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
        const std::array<const char*, 4> comparison_labels{{"M_000", "M_300/M_200", "M_010/M_000", "M_001/M_000"}};
        const qmomkit::RawMomentRequest comparison_request{std::vector<qmomkit::MultiIndex>(comparison_indices.begin(), comparison_indices.end())};

        cqmom_condensing_droplets::ensure_output_directory(output_options.directory);
        std::ofstream comparison_file(output_options.directory / "moment_comparison.csv");
        if (!comparison_file)
        {
            throw std::runtime_error("Could not open paper moment-comparison output.");
        }
        comparison_file << "time,quantity,qmomkit,monte_carlo\n";
        comparison_file << std::setprecision(17);

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
        const auto checked_ratio = [](const double numerator, const double denominator, const char* quantity, const double time)
        {
            if (!std::isfinite(denominator) || std::abs(denominator) <= std::numeric_limits<double>::min())
            {
                throw std::runtime_error("Invalid denominator for " + std::string(quantity) + " at time " + std::to_string(time) + ": " +
                                         std::to_string(denominator));
            }
            const double value = numerator / denominator;
            if (!std::isfinite(value))
            {
                throw std::runtime_error("Non-finite value for " + std::string(quantity) + " at time " + std::to_string(time));
            }
            return value;
        };
        const auto calculate_comparison_values = [&](const std::vector<double>& moments, const double time)
        {
            return std::array<double, 4>{moments[0], checked_ratio(moments[1], moments[2], "M_300/M_200", time),
                                         checked_ratio(moments[3], moments[0], "M_010/M_000", time),
                                         checked_ratio(moments[4], moments[0], "M_001/M_000", time)};
        };
        const auto write_comparison = [&](const double time, const qmomkit::CqmomQuadrature& current_quadrature)
        {
            const qmomkit::SelectedRawMomentResult qmomkit_moments = qmomkit::reconstruct_raw_moments(current_quadrature, comparison_request);
            if (!qmomkit_moments.ok())
            {
                throw std::runtime_error("CQMOM paper moment reconstruction failed: " + qmomkit_moments.message);
            }
            const qmomkit::SelectedRawMomentResult mc_moments = estimate_mc_moments(comparison_request);
            const std::array<double, 4> qmomkit_values = calculate_comparison_values(qmomkit_moments.values, time);
            const std::array<double, 4> mc_values = calculate_comparison_values(mc_moments.values, time);
            for (std::size_t index = 0; index < comparison_labels.size(); ++index)
            {
                comparison_file << time << ',' << comparison_labels[index] << ',' << qmomkit_values[index] << ',' << mc_values[index] << '\n';
            }
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
        write_comparison(initial_time, quadrature);

        const double output_interval = output_options.interval.value_or(final_time - initial_time);
        double time = initial_time;
        double next_output_time = initial_time + output_interval;
        while (time < final_time - 1.0e-12)
        {
            const double current_step = std::min(time_step, final_time - time);
            const double step_start_time = time;
            const qmomkit::CqmomPhaseSpaceEvolutionResult phase_space_evolved =
                qmomkit::advance_cqmom_advection(state, step_start_time, current_step, velocity, phase_space_evolution_options);
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

            qmomkit::monte_carlo::EnsembleEvolutionResult mc_advection =
                qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, step_start_time, step_start_time + current_step, mc_advection_options);
            if (!mc_advection.ok())
            {
                std::cerr << "Advection particle step failed: " << mc_advection.message << '\n';
                return 1;
            }
            ensemble = std::move(mc_advection.ensemble);

            if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
            {
                quadrature_history.emplace_back(time, quadrature);
                mc_final_moments = estimate_mc_moments(moment_output_request);
                if (mc_moment_history_writer.has_value())
                {
                    mc_moment_history_writer->write(time, qmomkit::create_cqmom_moment_data(moment_output_layout, mc_final_moments));
                }
                write_comparison(time, quadrature);
                while (next_output_time <= time + 1.0e-12)
                {
                    next_output_time += output_interval;
                }
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
        cqmom_condensing_droplets::print_result_summary(
            initial_moments, final_moments, quadrature, output_options.directory, cqmom_condensing_droplets::initialization_sample_count,
            cqmom_condensing_droplets::initialization_seed, time_step, final_time, output_options.write_moments, output_options.write_quadrature);
        cqmom_condensing_droplets::print_monte_carlo_result_summary(
            mc_initial_moments, mc_final_moments, output_options.directory, ensemble.size(), cqmom_condensing_droplets::initialization_sample_count,
            cqmom_condensing_droplets::initialization_seed, time_step, final_time, output_options.write_moments);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
