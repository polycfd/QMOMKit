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

        auto configured_cqmom = qmomkit::config::parse_cqmom_prepared(*configuration_text, configuration_file.string());
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
        const qmomkit::CqmomPreparedCase cqmom_case = std::move(sampled_case_result).value();

        const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(cqmom_case);
        if (!cqmom.ok())
        {
            std::cerr << "CQMOM inversion failed: " << cqmom.message << '\n';
            return 1;
        }

        qmomkit::CqmomState state = cqmom_case.initial_state();
        qmomkit::CqmomQuadrature quadrature = cqmom.quadrature;
        const qmomkit::CqmomMomentLayout moment_output_layout = qmomkit::create_cqmom_moment_layout(cqmom_case.hierarchy());
        const qmomkit::RawMomentRequest moment_output_request{moment_output_layout.required_indices.to_multi_indices()};
        const qmomkit::PhaseSpaceVelocity velocity = paper_2d_bubbles_cavitation::advection_velocity();
        qmomkit::CqmomPhaseSpaceEvolutionOptions phase_space_evolution_options = std::move(phase_space_evolution_configured).value();
        if (phase_space_evolution_options.node_options.scheme != qmomkit::NodeAdvectionScheme::rk23)
        {
            std::cerr << "The paper Monte Carlo reference requires RK23 node advection.\n";
            return 1;
        }
        const qmomkit::AdaptiveRkOptions mc_advection_options = phase_space_evolution_options.node_options.adaptive_options;
        std::cout << qmomkit::format_node_advection_options(phase_space_evolution_options.node_options);

        paper_2d_bubbles_cavitation::ensure_output_directory(output_options.directory);

        std::vector<std::pair<double, qmomkit::CqmomQuadrature>> quadrature_history{{initial_time, quadrature}};
        qmomkit::RawMomentTensor mc_final_moments = paper_2d_bubbles_cavitation::compute_particle_moments(ensemble, particle_weight, moment_extents);
        const qmomkit::RawMomentTensor mc_initial_moments = mc_final_moments;
        std::optional<qmomkit::MomentWriter> mc_moment_history_writer;
        if (output_options.write_moments)
        {
            mc_moment_history_writer.emplace(qmomkit::MomentWriter::for_cqmom(output_options.directory / "moment_history_mc.txt", moment_extents));
            mc_moment_history_writer->write(initial_time, mc_final_moments);
        }

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
            mc_initial_moments, mc_final_moments, moment_output_layout.required_indices.to_multi_indices(), output_options.directory, ensemble.size(),
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
