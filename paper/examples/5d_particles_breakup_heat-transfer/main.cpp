#include "monte_carlo/evolution.hpp"
#include "physics.hpp"
#include "qmomkit/qmomkit.hpp"
#include "sampling.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

int main()
{
    using namespace paper_5d_particles_breakup_heat_transfer;
    constexpr unsigned long long breakup_seed = 321;
    constexpr double initial_time = 0.0;
    constexpr double final_time = 0.1;
    constexpr double time_step = 0.001;

    try
    {
        const std::filesystem::path configuration_file = "configuration.qmom";
        const auto configured = qmomkit::config::load_cqmom_prepared(configuration_file);
        const auto evolution_configured = qmomkit::config::load_cqmom_phase_space_evolution_options(configuration_file);
        const auto output_configured = qmomkit::config::load_output_options(configuration_file);
        const auto configuration_text = read_configuration_text(configuration_file);
        if (!configured || !evolution_configured || !output_configured || !configuration_text)
        {
            std::cerr << "Could not load configuration.qmom.\n";
            return 1;
        }

        const auto initial_condition = qmomkit::config::parse_cqmom_initial_condition(*configuration_text, configuration_file.string());
        if (!initial_condition)
        {
            std::cerr << initial_condition.error() << '\n';
            return 1;
        }
        const auto sampled = sample_initial_population(initial_condition.value(), configured.value().hierarchy());
        if (!sampled || !sampled->ok())
        {
            std::cerr << "Sampled initialization failed.\n";
            return 1;
        }

        qmomkit::monte_carlo::Ensemble ensemble = sampled->sampling.ensemble;
        const double particle_weight = sampled->sampling.number_density / static_cast<double>(ensemble.size());
        qmomkit::RawMomentTensor mc_final_moments = ensemble_moments(ensemble, particle_weight, sampled->moments.moments.extents);
        const qmomkit::RawMomentTensor mc_initial_moments = mc_final_moments;
        if (maximum_tensor_difference(sampled->moments.moments, mc_initial_moments) != 0.0)
        {
            throw std::runtime_error("The sampled CQMOM and Monte Carlo initial moments differ.");
        }

        const auto sampled_case_result = qmomkit::config::with_initial_moments(configured.value(), sampled->moments.moments);
        if (!sampled_case_result)
        {
            std::cerr << sampled_case_result.error() << '\n';
            return 1;
        }
        const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(sampled_case_result.value());
        if (!initial.ok())
        {
            std::cerr << "Initial CQMOM inversion failed: " << initial.message << '\n';
            return 1;
        }

        qmomkit::CqmomState state = sampled_case_result.value().initial_state();
        qmomkit::CqmomQuadrature quadrature = initial.quadrature;
        const qmomkit::CqmomMomentLayout moment_output_layout = qmomkit::create_cqmom_moment_layout(sampled_case_result.value().hierarchy());
        const qmomkit::RawMomentRequest moment_output_request{moment_output_layout.required_indices.to_multi_indices()};
        const qmomkit::CqmomMomentSource source = breakup_source();
        const qmomkit::CqmomPhaseSpaceEvolutionOptions evolution_options = evolution_configured.value();
        if (evolution_options.node_options.scheme != qmomkit::NodeAdvectionScheme::rk23)
        {
            std::cerr << "The paper Monte Carlo reference requires RK23 node advection.\n";
            return 1;
        }
        const qmomkit::AdaptiveRkOptions mc_advection_options = evolution_options.node_options.adaptive_options;
        const qmomkit::OutputOptions output_options = output_configured.value();

        ensure_output_directory(output_options.directory);
        std::ofstream mc_generalized_moment_history_file(output_options.directory / "generalized_moment_history_mc.txt");
        if (!mc_generalized_moment_history_file)
        {
            throw std::runtime_error("Could not open Monte Carlo generalized-moment history output.");
        }
        mc_generalized_moment_history_file << "# time M[4/3,0,0,0,0] M[4/3,1,0,0,0]\n";
        mc_generalized_moment_history_file << std::setprecision(17);
        // The QMOMKit generalized moments are written through the native moment
        // writer; the ensemble has no such writer, so its two fractional-order
        // columns are recorded here.
        const auto write_mc_generalized_moment_history = [&](const double time)
        {
            const std::vector<double> mc_values = paper_moments(ensemble, particle_weight);
            if (mc_values.size() != paper_moment_labels().size())
            {
                throw std::runtime_error("Paper moment output columns are inconsistent.");
            }
            mc_generalized_moment_history_file << time << ' ' << mc_values[6] << ' ' << mc_values[8] << '\n';
            if (!mc_generalized_moment_history_file)
            {
                throw std::runtime_error("Failed while writing paper moment output.");
            }
        };

        std::vector<std::pair<double, qmomkit::CqmomQuadrature>> quadrature_history{{initial_time, quadrature}};
        std::optional<qmomkit::MomentWriter> mc_moment_history_writer;
        if (output_options.write_moments)
        {
            mc_moment_history_writer.emplace(
                qmomkit::MomentWriter::for_cqmom(output_options.directory / "moment_history_mc.txt", mc_final_moments.extents));
            mc_moment_history_writer->write(initial_time, mc_final_moments);
        }
        write_mc_generalized_moment_history(initial_time);

        std::mt19937_64 generator(breakup_seed);
        const double output_interval = output_options.interval.value_or(final_time - initial_time);
        double time = initial_time;
        double next_output = initial_time + output_interval;
        while (time < final_time - 1.0e-14)
        {
            const double current_step = std::min(time_step, final_time - time);
            const qmomkit::CqmomPhaseSpaceEvolutionResult phase_space_evolved =
                qmomkit::advance_cqmom_phase_space(state, time, current_step, particle_velocity(), source, evolution_options);
            if (!phase_space_evolved.ok())
            {
                std::cerr << "CQMOM phase-space evolution failed: " << phase_space_evolved.message << '\n';
                return 1;
            }
            if (std::abs(phase_space_evolved.final_time - (time + current_step)) > 1.0e-12 * std::max(1.0, std::abs(time + current_step)))
            {
                std::cerr << "CQMOM phase-space evolution did not reach the prescribed step endpoint.\n";
                return 1;
            }
            state = phase_space_evolved.state;
            quadrature = phase_space_evolved.state.quadrature();

            qmomkit::monte_carlo::EnsembleEvolutionResult first_advection =
                qmomkit::monte_carlo::evolve_ensemble(ensemble, particle_velocity(), time, time + 0.5 * current_step, mc_advection_options);
            if (!first_advection.ok())
            {
                std::cerr << "First MC advection half-step failed: " << first_advection.message << '\n';
                return 1;
            }
            ensemble = std::move(first_advection.ensemble);
            apply_breakup(ensemble, current_step, generator);
            qmomkit::monte_carlo::EnsembleEvolutionResult second_advection = qmomkit::monte_carlo::evolve_ensemble(
                ensemble, particle_velocity(), time + 0.5 * current_step, time + current_step, mc_advection_options);
            if (!second_advection.ok())
            {
                std::cerr << "Second MC advection half-step failed: " << second_advection.message << '\n';
                return 1;
            }
            ensemble = std::move(second_advection.ensemble);

            time += current_step;
            if (time + 1.0e-14 >= next_output || time >= final_time - 1.0e-14)
            {
                quadrature_history.emplace_back(time, quadrature);
                mc_final_moments = ensemble_moments(ensemble, particle_weight, mc_final_moments.extents);
                if (mc_moment_history_writer.has_value())
                {
                    mc_moment_history_writer->write(time, mc_final_moments);
                }
                write_mc_generalized_moment_history(time);
                while (next_output <= time + 1.0e-14)
                {
                    next_output += output_interval;
                }
            }
        }

        if (output_options.write_moments)
        {
            qmomkit::MomentWriter::for_cqmom(output_options.directory / "initial_moments.txt", moment_output_request)
                .write(initial_time, initial.quadrature);
            qmomkit::MomentWriter moment_history_writer =
                qmomkit::MomentWriter::for_cqmom(output_options.directory / "moment_history.txt", moment_output_request);
            for (const auto& [history_time, history_quadrature] : quadrature_history)
            {
                moment_history_writer.write(history_time, history_quadrature);
            }
        }
        if (output_options.write_quadrature)
        {
            qmomkit::QuadratureWriter::for_cqmom(output_options.directory / "initial_quadrature.txt", initial.quadrature.number_of_dimensions())
                .write(initial_time, initial.quadrature);
            qmomkit::QuadratureWriter quadrature_history_writer =
                qmomkit::QuadratureWriter::for_cqmom(output_options.directory / "quadrature_history.txt", initial.quadrature.number_of_dimensions());
            for (const auto& [history_time, history_quadrature] : quadrature_history)
            {
                quadrature_history_writer.write(history_time, history_quadrature);
            }
        }

        std::cout << "5D spherical breakup paper comparison completed.\n";
        std::cout << "Outputs written to " << output_options.directory << ".\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
