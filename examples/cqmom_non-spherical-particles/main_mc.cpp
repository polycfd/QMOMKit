#include "breakage_model.hpp"
#include "monte_carlo/evolution.hpp"
#include "particle_model.hpp"
#include "particles_mc.hpp"
#include "qmomkit/configuration.hpp"
#include "qmomkit/cqmom.hpp"
#include "qmomkit/output.hpp"
#include "sampling.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <random>
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
        const cqmom_non_spherical_particles::ParticleModelParameters particle_parameters = cqmom_non_spherical_particles::particle_model_parameters();
        cqmom_non_spherical_particles::validate_particle_model_parameters(particle_parameters);
        const cqmom_non_spherical_particles::BinaryFragmentationParameters breakage_parameters =
            cqmom_non_spherical_particles::binary_fragmentation_parameters();
        cqmom_non_spherical_particles::validate_binary_fragmentation_parameters(breakage_parameters);
        const qmomkit::PhaseSpaceVelocity velocity = cqmom_non_spherical_particles::particle_phase_space_velocity(particle_parameters);

        const std::optional<std::string> configuration_text = cqmom_non_spherical_particles::read_configuration_text(configuration_file);
        if (!configuration_text.has_value())
        {
            return 1;
        }

        auto configured = qmomkit::config::load_cqmom(configuration_file);
        if (!configured)
        {
            cqmom_non_spherical_particles::print_configuration_errors(configured.error());
            return 1;
        }
        const qmomkit::CqmomCase configured_case = std::move(configured).value();
        const qmomkit::CqmomMomentLayout terminal_moment_layout = qmomkit::create_cqmom_moment_layout(configured_case.hierarchy());
        // Reuse the CQMOM hierarchy only to select and order console moments. Dense
        // Monte Carlo file output retains its native tensor-storage ordering.
        const std::vector<qmomkit::MultiIndex> terminal_moment_indices =
            cqmom_non_spherical_particles::select_terminal_moment_indices(terminal_moment_layout.required_indices);

        auto initial_condition = qmomkit::config::parse_cqmom_initial_condition(*configuration_text, configuration_file.string());
        if (!initial_condition)
        {
            cqmom_non_spherical_particles::print_configuration_errors(initial_condition.error());
            return 1;
        }

        auto output_configured = qmomkit::config::load_output_options(configuration_file);
        if (!output_configured)
        {
            cqmom_non_spherical_particles::print_configuration_errors(output_configured.error());
            return 1;
        }
        const qmomkit::OutputOptions output_options = std::move(output_configured).value();

        const std::optional<qmomkit::monte_carlo::TensorSampledInitialization> sampled =
            cqmom_non_spherical_particles::sample_initial_distribution(std::move(initial_condition).value(), configured_case.hierarchy());
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
        // Each sampled particle carries n0 / N so the ensemble represents
        // unnormalized raw moments of the particle number-density distribution.
        const double particle_weight = sampled->sampling.number_density / static_cast<double>(sampled->sampling.ensemble.size());
        const std::vector<std::size_t> moment_extents = sampled->moments.moments.extents;
        const qmomkit::AdaptiveRkOptions continuous_options;
        std::mt19937_64 breakage_generator(cqmom_non_spherical_particles::breakage_seed);

        qmomkit::RawMomentTensor final_moments = cqmom_non_spherical_particles::compute_ensemble_moments(ensemble, particle_weight, moment_extents);
        const qmomkit::RawMomentTensor initial_moments = final_moments;

        std::optional<qmomkit::MomentWriter> moment_history_writer;
        if (output_options.write_moments)
        {
            cqmom_non_spherical_particles::ensure_output_directory(output_options.directory);
            moment_history_writer.emplace(qmomkit::MomentWriter::for_cqmom(output_options.directory / "moment_history_mc.txt", moment_extents));
            moment_history_writer->write(initial_time, final_moments);
        }

        const double output_interval = output_options.interval.value_or(final_time - initial_time);
        double time = initial_time;
        double next_output_time = initial_time + output_interval;
        while (time < final_time - 1.0e-12)
        {
            // The MC macro step mirrors the operator roles of the CQMOM transport:
            // continuous stages advance particle trajectories, while discontinuous
            // stages apply stochastic fragmentation and may increase the ensemble size.
            const double current_step = std::min(time_step, final_time - time);
            const double half_step = 0.5 * current_step;

            qmomkit::monte_carlo::EnsembleEvolutionResult first_continuous =
                qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, time, time + half_step, continuous_options);
            if (!first_continuous.ok())
            {
                std::cerr << "First continuous particle half-step failed: " << first_continuous.message << '\n';
                return 1;
            }
            ensemble = std::move(first_continuous.ensemble);

            cqmom_non_spherical_particles::apply_binary_breakage(ensemble, current_step, breakage_generator, breakage_parameters);

            qmomkit::monte_carlo::EnsembleEvolutionResult second_continuous =
                qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, time + half_step, time + current_step, continuous_options);
            if (!second_continuous.ok())
            {
                std::cerr << "Second continuous particle half-step failed: " << second_continuous.message << '\n';
                return 1;
            }
            ensemble = std::move(second_continuous.ensemble);
            time += current_step;

            if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
            {
                final_moments = cqmom_non_spherical_particles::compute_ensemble_moments(ensemble, particle_weight, moment_extents);
                if (moment_history_writer.has_value())
                {
                    moment_history_writer->write(time, final_moments);
                }
                while (next_output_time <= time + 1.0e-12)
                {
                    next_output_time += output_interval;
                }
            }
        }

        if (!moment_history_writer.has_value())
        {
            final_moments = cqmom_non_spherical_particles::compute_ensemble_moments(ensemble, particle_weight, moment_extents);
        }
        cqmom_non_spherical_particles::print_monte_carlo_result_summary(
            initial_moments, final_moments, terminal_moment_indices, output_options.directory, ensemble.size(),
            cqmom_non_spherical_particles::initialization_sample_count, cqmom_non_spherical_particles::initialization_seed, time_step, final_time,
            output_options.write_moments);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
