#include "droplet_model.hpp"
#include "monte_carlo/evolution.hpp"
#include "qmomkit/configuration.hpp"
#include "qmomkit/cqmom.hpp"
#include "qmomkit/output.hpp"
#include "qmomkit/qmomkit.hpp"
#include "sampling.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <utility>

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";
    const double initial_time = 0.0;
    const double final_time = 1.0;
    const double time_step = 0.01;

    try
    {
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

        auto output_configured = qmomkit::config::load_output_options(configuration_file);
        if (!output_configured)
        {
            cqmom_condensing_droplets::print_configuration_errors(output_configured.error());
            return 1;
        }
        auto transport_configured = qmomkit::config::load_cqmom_transport_options(configuration_file);
        if (!transport_configured)
        {
            cqmom_condensing_droplets::print_configuration_errors(transport_configured.error());
            return 1;
        }
        const qmomkit::CqmomCase cqmom_case = std::move(configured_cqmom).value();
        const qmomkit::OutputOptions output_options = std::move(output_configured).value();
        const qmomkit::CqmomTransportOptions transport_options = std::move(transport_configured).value();
        const qmomkit::CqmomMomentLayout moment_layout = qmomkit::create_cqmom_moment_layout(cqmom_case.hierarchy());
        const qmomkit::RawMomentRequest moment_request{moment_layout.required_indices};

        const std::optional<qmomkit::monte_carlo::TensorSampledInitialization> sampled =
            cqmom_condensing_droplets::sample_initial_distribution(std::move(initial_condition).value(), cqmom_case.hierarchy());
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
        const auto velocity = cqmom_condensing_droplets::droplet_phase_space_velocity();
        const qmomkit::AdaptiveRk23Options continuous_options = transport_options.node_options.rk23_options;

        cqmom_condensing_droplets::print_monte_carlo_run_configuration(cqmom_case, transport_options.node_options);

        qmomkit::SelectedRawMomentResult final_moments =
            cqmom_condensing_droplets::compute_ensemble_moments(ensemble, particle_weight, moment_request);
        const qmomkit::SelectedRawMomentResult initial_moments = final_moments;

        std::optional<qmomkit::MomentWriter> moment_history_writer;
        if (output_options.write_moments)
        {
            cqmom_condensing_droplets::ensure_output_directory(output_options.directory);
            moment_history_writer.emplace(qmomkit::MomentWriter::for_cqmom(output_options.directory / "moment_history_mc.txt", moment_layout));
            moment_history_writer->write(initial_time, qmomkit::create_cqmom_moment_data(moment_layout, final_moments));
        }

        const double output_interval = output_options.interval.value_or(final_time - initial_time);
        double time = initial_time;
        double next_output_time = initial_time + output_interval;
        while (time < final_time - 1.0e-12)
        {
            const double current_step = std::min(time_step, final_time - time);
            qmomkit::monte_carlo::EnsembleEvolutionResult continuous =
                qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, time, time + current_step, continuous_options);
            if (!continuous.ok())
            {
                std::cerr << "Continuous particle step failed: " << continuous.message << '\n';
                return 1;
            }
            ensemble = std::move(continuous.ensemble);
            time = continuous.final_time;

            if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
            {
                final_moments = cqmom_condensing_droplets::compute_ensemble_moments(ensemble, particle_weight, moment_request);
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
            final_moments = cqmom_condensing_droplets::compute_ensemble_moments(ensemble, particle_weight, moment_request);
        }
        cqmom_condensing_droplets::print_monte_carlo_result_summary(
            initial_moments, final_moments, output_options.directory, ensemble.size(), cqmom_condensing_droplets::initialization_sample_count,
            cqmom_condensing_droplets::initialization_seed, time_step, final_time, output_options.write_moments);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
