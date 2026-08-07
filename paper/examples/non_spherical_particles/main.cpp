#include "breakage_model.hpp"
#include "monte_carlo/evolution.hpp"
#include "monte_carlo/moments.hpp"
#include "particle_model.hpp"
#include "particles_mc.hpp"
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

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";
    const double initial_time = 0.0;
    constexpr double final_time = 1.0;
    constexpr std::size_t outer_step_count = 100;

    try
    {
        const cqmom_non_spherical_particles::ParticleModelParameters particle_parameters = cqmom_non_spherical_particles::particle_model_parameters();
        cqmom_non_spherical_particles::validate_particle_model_parameters(particle_parameters);
        const cqmom_non_spherical_particles::BinaryFragmentationParameters breakage_parameters =
            cqmom_non_spherical_particles::binary_fragmentation_parameters();
        cqmom_non_spherical_particles::validate_binary_fragmentation_parameters(breakage_parameters);
        const qmomkit::PhaseSpaceVelocity velocity = cqmom_non_spherical_particles::particle_phase_space_velocity(particle_parameters);

        auto configured = qmomkit::config::load_cqmom(configuration_file);
        if (!configured)
        {
            cqmom_non_spherical_particles::print_configuration_errors(configured.error());
            return 1;
        }
        auto output_configured = qmomkit::config::load_output_options(configuration_file);
        if (!output_configured)
        {
            cqmom_non_spherical_particles::print_configuration_errors(output_configured.error());
            return 1;
        }
        auto phase_space_evolution_configured = qmomkit::config::load_cqmom_phase_space_evolution_options(configuration_file);
        if (!phase_space_evolution_configured)
        {
            cqmom_non_spherical_particles::print_configuration_errors(phase_space_evolution_configured.error());
            return 1;
        }

        const std::optional<std::string> configuration_text = cqmom_non_spherical_particles::read_configuration_text(configuration_file);
        if (!configuration_text.has_value())
        {
            return 1;
        }
        auto initial_condition = qmomkit::config::parse_cqmom_initial_condition(*configuration_text, configuration_file.string());
        if (!initial_condition)
        {
            cqmom_non_spherical_particles::print_configuration_errors(initial_condition.error());
            return 1;
        }

        std::optional<qmomkit::monte_carlo::TensorSampledInitialization> sampled =
            cqmom_non_spherical_particles::sample_initial_distribution(std::move(initial_condition).value(), configured.value().hierarchy());
        if (!sampled.has_value())
        {
            return 1;
        }
        if (!sampled->ok())
        {
            std::cerr << sampled->sampling.message << sampled->moments.message << '\n';
            return 1;
        }

        // The inherited CQMOM and MC paths intentionally start from this one
        // fixed-seed ensemble; the paper CSV evaluates raw moments from it directly.
        qmomkit::monte_carlo::Ensemble ensemble = sampled->sampling.ensemble;
        const double particle_weight = sampled->sampling.number_density / static_cast<double>(ensemble.size());
        const std::vector<std::size_t> moment_extents = sampled->moments.moments.extents;
        auto sampled_case_result = qmomkit::config::with_initial_moments(configured.value(), std::move(sampled->moments.moments));
        if (!sampled_case_result)
        {
            cqmom_non_spherical_particles::print_configuration_errors(sampled_case_result.error());
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
        const qmomkit::RawMomentRequest file_moment_request{moment_output_layout.required_indices};
        const qmomkit::RawMomentRequest terminal_moment_request{
            cqmom_non_spherical_particles::select_terminal_moment_indices(moment_output_layout.required_indices)};
        const qmomkit::CqmomMomentSource source = cqmom_non_spherical_particles::binary_breakage_source(cqmom_case, breakage_parameters);
        const qmomkit::OutputOptions output_options = std::move(output_configured).value();
        const double outer_time_step = final_time / static_cast<double>(outer_step_count);
        const qmomkit::CqmomPhaseSpaceEvolutionOptions phase_space_evolution_options = std::move(phase_space_evolution_configured).value();
        if (phase_space_evolution_options.node_options.scheme != qmomkit::NodeAdvectionScheme::rk23)
        {
            std::cerr << "The paper Monte Carlo reference requires RK23 node advection.\n";
            return 1;
        }
        const qmomkit::AdaptiveRkOptions mc_advection_options = phase_space_evolution_options.node_options.adaptive_options;
        std::cout << qmomkit::format_cqmom_settings(cqmom_case.moment_storage(), cqmom_case.adaptive_branch_counts(), cqmom_case.vandermonde_solver(),
                                                    cqmom_case.diagnostic_level());
        std::cout << qmomkit::format_node_advection_options(phase_space_evolution_options.node_options);

        // The configured phase-space order is (d, chi, u, v, w). Theta_x,
        // Theta_y, and Theta_z are the componentwise raw-moment quantities
        // M_00200/M_00000 - (M_00100/M_00000)^2 and its y/z analogues.
        const std::array<qmomkit::MultiIndex, 13> comparison_indices{{
            {0, 0, 0, 0, 0},
            {3, 0, 0, 0, 0},
            {2, 0, 0, 0, 0},
            {0, 1, 0, 0, 0},
            {3, 0, 1, 0, 0},
            {3, 0, 0, 1, 0},
            {3, 0, 0, 0, 1},
            {0, 0, 2, 0, 0},
            {0, 0, 1, 0, 0},
            {0, 0, 0, 2, 0},
            {0, 0, 0, 1, 0},
            {0, 0, 0, 0, 2},
            {0, 0, 0, 0, 1},
        }};
        const std::array<const char*, 9> comparison_labels{{"M_00000", "M_30000/M_20000", "M_01000/M_00000", "M_30100/M_00000", "M_30010/M_00000",
                                                            "M_30001/M_00000", "Theta_x", "Theta_y", "Theta_z"}};
        const qmomkit::RawMomentRequest comparison_request{std::vector<qmomkit::MultiIndex>(comparison_indices.begin(), comparison_indices.end())};

        cqmom_non_spherical_particles::ensure_output_directory(output_options.directory);
        std::ofstream comparison_file(output_options.directory / "moment_comparison.csv");
        if (!comparison_file)
        {
            throw std::runtime_error("Could not open paper moment-comparison output.");
        }
        comparison_file << "time,quantity,qmomkit,monte_carlo\n";
        comparison_file << std::setprecision(17);

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
            const double m00000 = moments[0];
            const double mean_chi = checked_ratio(moments[3], m00000, "M_01000/M_00000", time);
            const double mean_u = checked_ratio(moments[8], m00000, "M_00100/M_00000", time);
            const double mean_v = checked_ratio(moments[10], m00000, "M_00010/M_00000", time);
            const double mean_w = checked_ratio(moments[12], m00000, "M_00001/M_00000", time);
            return std::array<double, 9>{m00000,
                                         checked_ratio(moments[1], moments[2], "M_30000/M_20000", time),
                                         mean_chi,
                                         checked_ratio(moments[4], m00000, "M_30100/M_00000", time),
                                         checked_ratio(moments[5], m00000, "M_30010/M_00000", time),
                                         checked_ratio(moments[6], m00000, "M_30001/M_00000", time),
                                         checked_ratio(moments[7], m00000, "Theta_x", time) - mean_u * mean_u,
                                         checked_ratio(moments[9], m00000, "Theta_y", time) - mean_v * mean_v,
                                         checked_ratio(moments[11], m00000, "Theta_z", time) - mean_w * mean_w};
        };
        const auto estimate_mc_moments = [&]()
        {
            const double number_density = particle_weight * static_cast<double>(ensemble.size());
            const qmomkit::monte_carlo::SelectedEnsembleMomentEstimate estimate =
                qmomkit::monte_carlo::compute_requested_ensemble_moments(ensemble, number_density, comparison_request);
            if (!estimate.ok())
            {
                throw std::runtime_error("Monte Carlo paper moment extraction failed: " + estimate.message);
            }
            return estimate.values;
        };
        const auto write_comparison = [&](const double time, const qmomkit::CqmomQuadrature& current_quadrature)
        {
            const qmomkit::SelectedRawMomentResult qmomkit_moments = qmomkit::reconstruct_raw_moments(current_quadrature, comparison_request);
            if (!qmomkit_moments.ok())
            {
                throw std::runtime_error("CQMOM paper moment reconstruction failed: " + qmomkit_moments.message);
            }
            const std::array<double, 9> qmomkit_values = calculate_comparison_values(qmomkit_moments.values, time);
            const std::array<double, 9> mc_values = calculate_comparison_values(estimate_mc_moments(), time);
            for (std::size_t index = 0; index < comparison_labels.size(); ++index)
            {
                comparison_file << time << ',' << comparison_labels[index] << ',' << qmomkit_values[index] << ',' << mc_values[index] << '\n';
            }
        };

        std::vector<std::pair<double, qmomkit::CqmomQuadrature>> quadrature_history{{initial_time, quadrature}};
        qmomkit::RawMomentTensor mc_final_moments =
            cqmom_non_spherical_particles::compute_ensemble_moments(ensemble, particle_weight, moment_extents);
        const qmomkit::RawMomentTensor mc_initial_moments = mc_final_moments;
        std::optional<qmomkit::MomentWriter> mc_moment_history_writer;
        if (output_options.write_moments)
        {
            mc_moment_history_writer.emplace(qmomkit::MomentWriter::for_cqmom(output_options.directory / "moment_history_mc.txt", moment_extents));
            mc_moment_history_writer->write(initial_time, mc_final_moments);
        }
        write_comparison(initial_time, quadrature);

        std::mt19937_64 breakage_generator(cqmom_non_spherical_particles::breakage_seed);
        for (std::size_t outer_step = 0; outer_step < outer_step_count; ++outer_step)
        {
            const double step_start_time = initial_time + static_cast<double>(outer_step) * outer_time_step;
            const double output_time = outer_step + 1 == outer_step_count
                                           ? final_time
                                           : initial_time + static_cast<double>(outer_step + 1) * outer_time_step;
            const double current_step = output_time - step_start_time;
            const qmomkit::CqmomPhaseSpaceEvolutionResult phase_space_evolved =
                qmomkit::advance_cqmom_phase_space(state, step_start_time, current_step, velocity, source, phase_space_evolution_options);
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

            const double half_step = 0.5 * current_step;
            qmomkit::monte_carlo::EnsembleEvolutionResult first_advection =
                qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, step_start_time, step_start_time + half_step, mc_advection_options);
            if (!first_advection.ok())
            {
                std::cerr << "First advection particle half-step failed: " << first_advection.message << '\n';
                return 1;
            }
            ensemble = std::move(first_advection.ensemble);
            cqmom_non_spherical_particles::apply_binary_breakage(ensemble, current_step, breakage_generator, breakage_parameters);
            qmomkit::monte_carlo::EnsembleEvolutionResult second_advection = qmomkit::monte_carlo::evolve_ensemble(
                ensemble, velocity, step_start_time + half_step, output_time, mc_advection_options);
            if (!second_advection.ok())
            {
                std::cerr << "Second advection particle half-step failed: " << second_advection.message << '\n';
                return 1;
            }
            ensemble = std::move(second_advection.ensemble);

            quadrature_history.emplace_back(output_time, quadrature);
            mc_final_moments = cqmom_non_spherical_particles::compute_ensemble_moments(ensemble, particle_weight, moment_extents);
            if (mc_moment_history_writer.has_value())
            {
                mc_moment_history_writer->write(output_time, mc_final_moments);
            }
            write_comparison(output_time, quadrature);
        }

        if (output_options.write_moments)
        {
            qmomkit::MomentWriter::for_cqmom(output_options.directory / "initial_moments.txt", file_moment_request)
                .write(initial_time, cqmom.quadrature);
            qmomkit::MomentWriter moment_history_writer =
                qmomkit::MomentWriter::for_cqmom(output_options.directory / "moment_history.txt", file_moment_request);
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

        const qmomkit::SelectedRawMomentResult initial_moments = qmomkit::reconstruct_raw_moments(cqmom.quadrature, terminal_moment_request);
        const qmomkit::SelectedRawMomentResult final_moments = qmomkit::reconstruct_raw_moments(quadrature, terminal_moment_request);
        if (!initial_moments.ok() || !final_moments.ok())
        {
            std::cerr << "CQMOM moment reconstruction failed while producing the summary.\n";
            return 1;
        }
        cqmom_non_spherical_particles::print_result_summary(
            initial_moments, final_moments, quadrature, output_options.directory, cqmom_non_spherical_particles::initialization_sample_count,
            cqmom_non_spherical_particles::initialization_seed, outer_time_step, final_time, output_options.write_moments,
            output_options.write_quadrature);
        cqmom_non_spherical_particles::print_monte_carlo_result_summary(
            mc_initial_moments, mc_final_moments, terminal_moment_request.indices, output_options.directory, ensemble.size(),
            cqmom_non_spherical_particles::initialization_sample_count, cqmom_non_spherical_particles::initialization_seed, outer_time_step,
            final_time, output_options.write_moments);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
