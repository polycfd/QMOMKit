#include "monte_carlo/evolution.hpp"
#include "monte_carlo/moments.hpp"
#include "particles_mc.hpp"
#include "physics.hpp"
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
    constexpr double final_time = 13.9;
    constexpr double time_step = 0.139;

    try
    {
        auto output_configured = qmomkit::config::load_output_options(configuration_file);
        if (!output_configured)
        {
            paper_2d_bubbles_cavitation::print_configuration_errors(output_configured.error());
            return 1;
        }
        auto phase_space_evolution_configured = qmomkit::config::load_cqmom_phase_space_evolution_options(configuration_file);
        if (!phase_space_evolution_configured)
        {
            paper_2d_bubbles_cavitation::print_configuration_errors(phase_space_evolution_configured.error());
            return 1;
        }
        const qmomkit::OutputOptions output_options = std::move(output_configured).value();
        const std::optional<std::string> configuration_text = paper_2d_bubbles_cavitation::read_configuration_text(configuration_file);
        if (!configuration_text.has_value())
        {
            return 1;
        }

        auto configured_cqmom = qmomkit::config::parse_cqmom(*configuration_text, configuration_file.string());
        if (!configured_cqmom)
        {
            paper_2d_bubbles_cavitation::print_configuration_errors(configured_cqmom.error());
            return 1;
        }

        auto initial_condition = qmomkit::config::parse_cqmom_initial_condition(*configuration_text, configuration_file.string());
        if (!initial_condition)
        {
            paper_2d_bubbles_cavitation::print_configuration_errors(initial_condition.error());
            return 1;
        }

        std::optional<qmomkit::monte_carlo::TensorSampledInitialization> sampled =
            paper_2d_bubbles_cavitation::sample_initial_distribution(std::move(initial_condition).value(), configured_cqmom.value().hierarchy());
        if (!sampled.has_value())
        {
            return 1;
        }
        if (!sampled->ok())
        {
            std::cerr << sampled->sampling.message << sampled->moments.message << '\n';
            return 1;
        }

        // The paper output compares the CQMOM state with the same finite MC
        // realization that supplied its initial raw moments.
        qmomkit::monte_carlo::Ensemble ensemble = sampled->sampling.ensemble;
        const double particle_weight = sampled->sampling.number_density / static_cast<double>(ensemble.size());
        const std::vector<std::size_t> moment_extents = sampled->moments.moments.extents;
        auto sampled_case_result = qmomkit::config::with_initial_moments(configured_cqmom.value(), std::move(sampled->moments.moments));
        if (!sampled_case_result)
        {
            paper_2d_bubbles_cavitation::print_configuration_errors(sampled_case_result.error());
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
        const qmomkit::PhaseSpaceVelocity velocity = paper_2d_bubbles_cavitation::advection_velocity();
        qmomkit::CqmomPhaseSpaceEvolutionOptions phase_space_evolution_options = std::move(phase_space_evolution_configured).value();
        if (phase_space_evolution_options.node_options.scheme != qmomkit::NodeAdvectionScheme::rk23)
        {
            std::cerr << "The paper Monte Carlo reference requires RK23 node advection.\n";
            return 1;
        }
        const qmomkit::AdaptiveRkOptions mc_advection_options = phase_space_evolution_options.node_options.adaptive_options;
        std::cout << qmomkit::format_node_advection_options(phase_space_evolution_options.node_options);

        // The configured phase-space order is (R, U): population-normalized
        // first, mixed, and second raw moments are compared directly.
        const std::array<qmomkit::MultiIndex, 6> comparison_indices{{{0, 0}, {1, 0}, {0, 1}, {1, 1}, {2, 0}, {0, 2}}};
        const std::array<const char*, 6> comparison_labels{{"M_00", "M_10/M_00", "M_01/M_00", "M_11/M_00", "M_20/M_00", "M_02/M_00"}};
        const qmomkit::RawMomentRequest comparison_request{std::vector<qmomkit::MultiIndex>(comparison_indices.begin(), comparison_indices.end())};

        paper_2d_bubbles_cavitation::ensure_output_directory(output_options.directory);
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
            return std::array<double, 6>{moments[0],
                                         checked_ratio(moments[1], moments[0], "M_10/M_00", time),
                                         checked_ratio(moments[2], moments[0], "M_01/M_00", time),
                                         checked_ratio(moments[3], moments[0], "M_11/M_00", time),
                                         checked_ratio(moments[4], moments[0], "M_20/M_00", time),
                                         checked_ratio(moments[5], moments[0], "M_02/M_00", time)};
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
            const std::array<double, 6> qmomkit_values = calculate_comparison_values(qmomkit_moments.values, time);
            const std::array<double, 6> mc_values = calculate_comparison_values(estimate_mc_moments(), time);
            for (std::size_t index = 0; index < comparison_labels.size(); ++index)
            {
                comparison_file << time << ',' << comparison_labels[index] << ',' << qmomkit_values[index] << ',' << mc_values[index] << '\n';
            }
        };

        std::vector<std::pair<double, qmomkit::CqmomQuadrature>> quadrature_history{{initial_time, quadrature}};
        qmomkit::RawMomentTensor mc_final_moments = paper_2d_bubbles_cavitation::compute_particle_moments(ensemble, particle_weight, moment_extents);
        const qmomkit::RawMomentTensor mc_initial_moments = mc_final_moments;
        std::optional<qmomkit::MomentWriter> mc_moment_history_writer;
        if (output_options.write_moments)
        {
            mc_moment_history_writer.emplace(qmomkit::MomentWriter::for_cqmom(output_options.directory / "moment_history_mc.txt", moment_extents));
            mc_moment_history_writer->write(initial_time, mc_final_moments);
        }
        write_comparison(initial_time, quadrature);

        double time = initial_time;
        while (time < final_time)
        {
            const double current_step = std::min(time_step, final_time - time);
            const qmomkit::CqmomPhaseSpaceEvolutionResult phase_space_evolved =
                qmomkit::advance_cqmom_advection(state, time, current_step, velocity, phase_space_evolution_options);
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
            if (std::abs(phase_space_evolved.final_time - (time + current_step)) > 1.0e-12 * std::max(1.0, std::abs(time + current_step)))
            {
                std::cerr << "CQMOM phase-space evolution did not reach the prescribed step endpoint.\n";
                return 1;
            }

            qmomkit::monte_carlo::EnsembleEvolutionResult mc_advection =
                qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, time, time + current_step, mc_advection_options);
            if (!mc_advection.ok())
            {
                std::cerr << "Advection particle step failed: " << mc_advection.message << '\n';
                return 1;
            }
            ensemble = std::move(mc_advection.ensemble);

            time += current_step;
            quadrature_history.emplace_back(time, quadrature);
            mc_final_moments = paper_2d_bubbles_cavitation::compute_particle_moments(ensemble, particle_weight, moment_extents);
            if (mc_moment_history_writer.has_value())
            {
                mc_moment_history_writer->write(time, mc_final_moments);
            }
            write_comparison(time, quadrature);
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
        const qmomkit::SelectedRawMomentResult final_moments = qmomkit::reconstruct_raw_moments(quadrature, moment_output_request);
        if (!initial_moments.ok() || !final_moments.ok())
        {
            std::cerr << "CQMOM moment reconstruction failed while producing the summary.\n";
            return 1;
        }
        paper_2d_bubbles_cavitation::print_result_summary(
            initial_moments, final_moments, quadrature, output_options.directory, paper_2d_bubbles_cavitation::initialization_sample_count,
            paper_2d_bubbles_cavitation::initialization_seed, time_step, final_time, output_options.write_moments, output_options.write_quadrature);
        paper_2d_bubbles_cavitation::print_monte_carlo_result_summary(
            mc_initial_moments, mc_final_moments, moment_output_layout.required_indices, output_options.directory, ensemble.size(),
            paper_2d_bubbles_cavitation::initialization_sample_count, paper_2d_bubbles_cavitation::initialization_seed, time_step, final_time,
            output_options.write_moments);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
