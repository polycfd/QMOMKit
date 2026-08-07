#include "droplet_model.hpp"
#include "monte_carlo/birth.hpp"
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
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{

struct NucleationParameters
{
    double radius = 50.0e-5;
    double vertical_velocity = 0.0;
    double temperature = 298.15;
    double start_time = 0.20;
    double end_time = 0.40;
    double relative_population_increase = 1.0;
};

constexpr NucleationParameters nucleation_parameters{};

[[nodiscard]] double integer_power(const double base, const std::size_t exponent) noexcept
{
    double result = 1.0;
    for (std::size_t power = 0; power < exponent; ++power)
    {
        result *= base;
    }
    return result;
}

} // namespace

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";
    const double initial_time = 0.0;
    constexpr double final_time = 1.0;
    constexpr std::size_t outer_step_count = 100;

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
        const double outer_time_step = final_time / static_cast<double>(outer_step_count);

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
        double number_density = sampled->sampling.number_density;
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

        const qmomkit::RawMomentRequest initial_number_request{std::vector<qmomkit::MultiIndex>{{0, 0, 0}}};
        const qmomkit::SelectedRawMomentResult initial_number_moment = qmomkit::reconstruct_raw_moments(cqmom.quadrature, initial_number_request);
        if (!initial_number_moment.ok() || initial_number_moment.values.size() != 1 || !std::isfinite(initial_number_moment.values.front()) ||
            initial_number_moment.values.front() <= 0.0)
        {
            std::cerr << "Could not obtain a positive initial M_000 for the nucleation flux.\n";
            return 1;
        }
        const double initial_number_density = initial_number_moment.values.front();
        const double nucleation_duration = nucleation_parameters.end_time - nucleation_parameters.start_time;
        const double nucleation_flux = nucleation_parameters.relative_population_increase * initial_number_density / nucleation_duration;
        if (!std::isfinite(nucleation_flux) || nucleation_flux <= 0.0)
        {
            std::cerr << "Invalid nucleation flux.\n";
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
        if (phase_space_evolution_options.scheme != qmomkit::PhaseSpaceEvolutionScheme::strang)
        {
            std::cerr << "The condensing-droplet nucleation example requires the configured Strang splitting scheme.\n";
            return 1;
        }
        const qmomkit::AdaptiveRkOptions mc_advection_options = phase_space_evolution_options.node_options.adaptive_options;

        const auto aligned_outer_step = [&](const double time)
        {
            const double unrounded_step = (time - initial_time) / outer_time_step;
            const auto outer_step = static_cast<std::size_t>(std::llround(unrounded_step));
            if (std::abs(unrounded_step - static_cast<double>(outer_step)) > 1.0e-12 * std::max(1.0, std::abs(unrounded_step)))
            {
                throw std::runtime_error("The nucleation interval does not align with the application outer time grid.");
            }
            return outer_step;
        };
        const std::size_t nucleation_start_step = aligned_outer_step(nucleation_parameters.start_time);
        const std::size_t nucleation_end_step = aligned_outer_step(nucleation_parameters.end_time);
        if (nucleation_end_step <= nucleation_start_step || nucleation_end_step > outer_step_count)
        {
            std::cerr << "Invalid nucleation interval on the application outer time grid.\n";
            return 1;
        }

        const qmomkit::CqmomMomentSource zero_moment_source = [](double, const qmomkit::CqmomSourceContext& context)
        { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };
        const qmomkit::monte_carlo::BirthCoordinateSampler sample_nucleated_droplet = [](std::mt19937_64&)
        {
            return qmomkit::monte_carlo::Sample{nucleation_parameters.radius, nucleation_parameters.vertical_velocity,
                                                nucleation_parameters.temperature};
        };

        cqmom_condensing_droplets::print_cqmom_run_configuration(cqmom_case, phase_space_evolution_options);

        // The configured phase-space order is (r, w, T). The requested
        // Sauter-mean quantity is M_300/M_200 without a factor of two.
        const std::array<qmomkit::MultiIndex, 7> comparison_indices{{{0, 0, 0}, {3, 0, 0}, {2, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 0}, {1, 0, 1}}};
        const std::array<const char*, 6> comparison_labels{{"M_000", "M_300/M_200", "M_010/M_000", "M_001/M_000", "sigma_r^2", "cov_r_T"}};
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
            const double mean_radius = checked_ratio(moments[5], moments[0], "M_100/M_000", time);
            const double mean_temperature = checked_ratio(moments[4], moments[0], "M_001/M_000", time);
            const double radius_variance = checked_ratio(moments[2], moments[0], "M_200/M_000", time) - mean_radius * mean_radius;
            const double radius_temperature_covariance = checked_ratio(moments[6], moments[0], "M_101/M_000", time) - mean_radius * mean_temperature;
            if (!std::isfinite(radius_variance) || !std::isfinite(radius_temperature_covariance))
            {
                throw std::runtime_error("Non-finite derived comparison quantity at time " + std::to_string(time));
            }
            return std::array<double, 6>{moments[0],
                                         checked_ratio(moments[1], moments[2], "M_300/M_200", time),
                                         checked_ratio(moments[3], moments[0], "M_010/M_000", time),
                                         mean_temperature,
                                         radius_variance,
                                         radius_temperature_covariance};
        };
        const auto write_comparison = [&](const double time, const qmomkit::CqmomQuadrature& current_quadrature)
        {
            const qmomkit::SelectedRawMomentResult qmomkit_moments = qmomkit::reconstruct_raw_moments(current_quadrature, comparison_request);
            if (!qmomkit_moments.ok())
            {
                throw std::runtime_error("CQMOM paper moment reconstruction failed: " + qmomkit_moments.message);
            }
            const qmomkit::SelectedRawMomentResult mc_moments = estimate_mc_moments(comparison_request);
            const std::array<double, 6> qmomkit_values = calculate_comparison_values(qmomkit_moments.values, time);
            const std::array<double, 6> mc_values = calculate_comparison_values(mc_moments.values, time);
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

        for (std::size_t outer_step = 0; outer_step < outer_step_count; ++outer_step)
        {
            const double step_start_time = initial_time + static_cast<double>(outer_step) * outer_time_step;
            const double output_time =
                outer_step + 1 == outer_step_count ? final_time : initial_time + static_cast<double>(outer_step + 1) * outer_time_step;
            const double current_step = output_time - step_start_time;
            const bool nucleation_active = outer_step >= nucleation_start_step && outer_step < nucleation_end_step;
            const qmomkit::CqmomMomentBoundaryContribution boundary_contribution =
                [nucleation_active, nucleation_flux](double, const qmomkit::CqmomSourceContext& context)
            {
                std::vector<double> rates(context.size(), 0.0);
                if (!nucleation_active)
                {
                    return qmomkit::CqmomMomentRate{std::move(rates)};
                }
                for (std::size_t index = 0; index < context.size(); ++index)
                {
                    const qmomkit::MultiIndex& order = context.index(index);
                    if (order.size() != 3)
                    {
                        throw std::runtime_error("The condensing-droplet nucleation contribution requires the phase-space order (r,w,T).");
                    }
                    rates[index] = nucleation_flux * integer_power(nucleation_parameters.radius, order[0]) *
                                   integer_power(nucleation_parameters.vertical_velocity, order[1]) *
                                   integer_power(nucleation_parameters.temperature, order[2]);
                }
                return qmomkit::CqmomMomentRate{std::move(rates)};
            };
            const qmomkit::CqmomPhaseSpaceEvolutionResult phase_space_evolved =
                nucleation_active ? qmomkit::advance_cqmom_phase_space(state, step_start_time, current_step, velocity, zero_moment_source,
                                                                       boundary_contribution, phase_space_evolution_options)
                                  : qmomkit::advance_cqmom_advection(state, step_start_time, current_step, velocity, phase_space_evolution_options);
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
            if (std::abs(phase_space_evolved.final_time - output_time) > 1.0e-12 * std::max(1.0, std::abs(output_time)))
            {
                std::cerr << "CQMOM phase-space evolution did not reach the prescribed outer-step endpoint.\n";
                return 1;
            }

            if (nucleation_active)
            {
                const double source_time = step_start_time + 0.5 * current_step;
                qmomkit::monte_carlo::EnsembleEvolutionResult mc_first_half_advection =
                    qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, step_start_time, source_time, mc_advection_options);
                if (!mc_first_half_advection.ok())
                {
                    std::cerr << "First-half Monte Carlo advection failed: " << mc_first_half_advection.message << '\n';
                    return 1;
                }
                ensemble = std::move(mc_first_half_advection.ensemble);
                const std::size_t ensemble_size_before_birth = ensemble.size();
                const unsigned long long birth_seed = cqmom_condensing_droplets::initialization_seed + static_cast<unsigned long long>(outer_step);
                const qmomkit::monte_carlo::PopulationBirthResult birth = qmomkit::monte_carlo::inject_population(
                    ensemble, number_density, nucleation_flux * current_step, sample_nucleated_droplet, birth_seed);
                if (!birth.ok())
                {
                    std::cerr << "Monte Carlo population injection failed: " << birth.message << '\n';
                    return 1;
                }
                for (std::size_t sample = ensemble_size_before_birth; sample < birth.ensemble.size(); ++sample)
                {
                    if (birth.ensemble[sample] != qmomkit::monte_carlo::Sample{nucleation_parameters.radius, nucleation_parameters.vertical_velocity,
                                                                               nucleation_parameters.temperature})
                    {
                        std::cerr << "Monte Carlo population injection returned an unexpected birth coordinate.\n";
                        return 1;
                    }
                }
                ensemble = birth.ensemble;
                number_density = birth.number_density;
                qmomkit::monte_carlo::EnsembleEvolutionResult mc_second_half_advection =
                    qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, source_time, output_time, mc_advection_options);
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
                    qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, step_start_time, output_time, mc_advection_options);
                if (!mc_advection.ok())
                {
                    std::cerr << "Monte Carlo advection failed: " << mc_advection.message << '\n';
                    return 1;
                }
                ensemble = std::move(mc_advection.ensemble);
            }

            quadrature_history.emplace_back(output_time, quadrature);
            mc_final_moments = estimate_mc_moments(moment_output_request);
            if (mc_moment_history_writer.has_value())
            {
                mc_moment_history_writer->write(output_time, qmomkit::create_cqmom_moment_data(moment_output_layout, mc_final_moments));
            }
            write_comparison(output_time, quadrature);
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
        cqmom_condensing_droplets::print_result_summary(initial_moments, final_moments, quadrature, output_options.directory,
                                                        cqmom_condensing_droplets::initialization_sample_count,
                                                        cqmom_condensing_droplets::initialization_seed, outer_time_step, final_time,
                                                        output_options.write_moments, output_options.write_quadrature);
        cqmom_condensing_droplets::print_monte_carlo_result_summary(
            mc_initial_moments, mc_final_moments, output_options.directory, ensemble.size(), cqmom_condensing_droplets::initialization_sample_count,
            cqmom_condensing_droplets::initialization_seed, outer_time_step, final_time, output_options.write_moments);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
