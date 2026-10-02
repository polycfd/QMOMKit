#include "monte_carlo/evolution.hpp"
#include "monte_carlo/output.hpp"
#include "particles_mc.hpp"
#include "physics.hpp"
#include "qmomkit/configuration.hpp"
#include "qmomkit/output.hpp"
#include "sampling.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
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
        const std::optional<std::string> configuration_text = cqmom_2d_bubbles_cavitation::read_configuration_text(configuration_file);
        if (!configuration_text.has_value())
        {
            return 1;
        }

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

        auto output_configured = qmomkit::config::load_output_options(configuration_file);
        if (!output_configured)
        {
            cqmom_2d_bubbles_cavitation::print_configuration_errors(output_configured.error());
            return 1;
        }
        const qmomkit::OutputOptions output_options = std::move(output_configured).value();

        const std::optional<qmomkit::monte_carlo::TensorSampledInitialization> sampled =
            cqmom_2d_bubbles_cavitation::sample_initial_distribution(std::move(initial_condition).value(), configured_cqmom.value().hierarchy());
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
        const double particle_weight = sampled->sampling.number_density / static_cast<double>(sampled->sampling.ensemble.size());
        const std::vector<std::size_t> moment_extents = sampled->moments.moments.extents;
        const qmomkit::CqmomMomentLayout moment_output_layout = qmomkit::create_cqmom_moment_layout(configured_cqmom.value().hierarchy());
        const qmomkit::PhaseSpaceVelocity velocity = cqmom_2d_bubbles_cavitation::advection_velocity();
        const qmomkit::AdaptiveRkOptions advection_options;

        qmomkit::RawMomentTensor final_moments = cqmom_2d_bubbles_cavitation::compute_particle_moments(ensemble, particle_weight, moment_extents);
        const qmomkit::RawMomentTensor initial_moments = final_moments;

        std::optional<qmomkit::MomentWriter> moment_history_writer;
        if (output_options.write_moments)
        {
            cqmom_2d_bubbles_cavitation::ensure_output_directory(output_options.directory);
            moment_history_writer.emplace(
                qmomkit::monte_carlo::make_moment_history_writer(output_options.directory / "moment_history_mc.txt", moment_extents, output_options));
            qmomkit::monte_carlo::write_moment_history_row(*moment_history_writer, initial_time, final_moments, ensemble,
                                                           particle_weight * static_cast<double>(ensemble.size()), output_options);
        }

        const double output_interval = output_options.interval.value_or(final_time - initial_time);
        double time = initial_time;
        double next_output_time = initial_time + output_interval;
        while (time < final_time - 1.0e-12)
        {
            const double current_step = std::min(time_step, final_time - time);
            qmomkit::monte_carlo::EnsembleEvolutionResult advection =
                qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, time, time + current_step, advection_options);
            if (!advection.ok())
            {
                std::cerr << "Advection particle step failed: " << advection.message << '\n';
                return 1;
            }
            ensemble = std::move(advection.ensemble);
            time += current_step;

            if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
            {
                final_moments = cqmom_2d_bubbles_cavitation::compute_particle_moments(ensemble, particle_weight, moment_extents);
                if (moment_history_writer.has_value())
                {
                    qmomkit::monte_carlo::write_moment_history_row(*moment_history_writer, time, final_moments, ensemble,
                                                                   particle_weight * static_cast<double>(ensemble.size()), output_options);
                }
                while (next_output_time <= time + 1.0e-12)
                {
                    next_output_time += output_interval;
                }
            }
        }

        if (!moment_history_writer.has_value())
        {
            final_moments = cqmom_2d_bubbles_cavitation::compute_particle_moments(ensemble, particle_weight, moment_extents);
        }
        cqmom_2d_bubbles_cavitation::print_monte_carlo_result_summary(
            initial_moments, final_moments, moment_output_layout.required_indices.to_multi_indices(), output_options.directory, ensemble.size(),
            cqmom_2d_bubbles_cavitation::initialization_sample_count, cqmom_2d_bubbles_cavitation::initialization_seed, time_step, final_time,
            output_options.write_moments);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
