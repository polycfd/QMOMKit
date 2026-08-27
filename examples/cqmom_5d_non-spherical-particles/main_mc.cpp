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
#include <string>
#include <utility>
#include <vector>

namespace
{
constexpr double initial_time = 0.0;
constexpr double final_time = 1.0;
constexpr double time_step = 0.005;
} // namespace

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <path-to-configuration.qmom>\n";
        return 1;
    }

    const std::filesystem::path configuration_file = argv[1];

    try
    {
        using namespace cqmom_5d_non_spherical_particles;

        auto configured = qmomkit::config::load_cqmom(configuration_file);
        if (!configured)
        {
            print_configuration_errors(configured.error());
            return 1;
        }
        const qmomkit::CqmomCase configured_case = std::move(configured).value();

        auto output_configured = qmomkit::config::load_output_options(configuration_file);
        if (!output_configured)
        {
            print_configuration_errors(output_configured.error());
            return 1;
        }
        const qmomkit::OutputOptions output_options = std::move(output_configured).value();

        auto evolution_configured = qmomkit::config::load_cqmom_phase_space_evolution_options(configuration_file);
        if (!evolution_configured)
        {
            print_configuration_errors(evolution_configured.error());
            return 1;
        }
        const qmomkit::AdaptiveRkOptions advection_options = std::move(evolution_configured).value().node_options.adaptive_options;

        const std::optional<std::string> configuration_text = read_configuration_text(configuration_file);
        if (!configuration_text.has_value())
        {
            return 1;
        }
        auto initial_condition = qmomkit::config::parse_cqmom_initial_condition(*configuration_text, configuration_file.string());
        if (!initial_condition)
        {
            print_configuration_errors(initial_condition.error());
            return 1;
        }
        const auto sampled = sample_initial_distribution(std::move(initial_condition).value(), configured_case.hierarchy());
        if (!sampled.has_value())
        {
            return 1;
        }
        if (!sampled->ok())
        {
            std::cerr << sampled->sampling.message << sampled->moments.message << '\n';
            return 1;
        }

        const qmomkit::CqmomMomentLayout moment_layout = qmomkit::create_cqmom_moment_layout(configured_case.hierarchy());
        const std::vector<qmomkit::MultiIndex> terminal_moment_indices = select_terminal_moment_indices(moment_layout.required_indices);
        qmomkit::monte_carlo::Ensemble ensemble = sampled->sampling.ensemble;
        const double particle_weight = sampled->sampling.number_density / static_cast<double>(ensemble.size());
        const std::vector<std::size_t> moment_extents = sampled->moments.moments.extents;

        const ParticleModelParameters particle_parameters = particle_model_parameters();
        validate_particle_model_parameters(particle_parameters);
        const BinaryFragmentationParameters breakage_parameters = binary_fragmentation_parameters();
        validate_binary_fragmentation_parameters(breakage_parameters);
        const qmomkit::PhaseSpaceVelocity velocity = particle_phase_space_velocity(particle_parameters);
        std::mt19937_64 breakage_generator(breakage_seed);

        qmomkit::RawMomentTensor final_moments = compute_ensemble_moments(ensemble, particle_weight, moment_extents);
        const qmomkit::RawMomentTensor initial_moments = final_moments;

        std::optional<qmomkit::MomentWriter> moment_history_writer;
        if (output_options.write_moments)
        {
            ensure_output_directory(output_options.directory);
            moment_history_writer.emplace(qmomkit::MomentWriter::for_cqmom(output_options.directory / "moment_history_mc.txt", moment_extents));
            moment_history_writer->write(initial_time, final_moments);
        }

        const double output_interval = output_options.interval.value_or(final_time - initial_time);
        double time = initial_time;
        double next_output_time = initial_time + output_interval;
        while (time < final_time - 1.0e-12)
        {
            const double step = std::min(time_step, final_time - time);
            const double half_step = 0.5 * step;

            qmomkit::monte_carlo::EnsembleEvolutionResult first_advection =
                qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, time, time + half_step, advection_options);
            if (!first_advection.ok())
            {
                std::cerr << "First advection particle half-step failed: " << first_advection.message << '\n';
                return 1;
            }
            ensemble = std::move(first_advection.ensemble);

            apply_binary_breakage(ensemble, step, breakage_generator, breakage_parameters);

            qmomkit::monte_carlo::EnsembleEvolutionResult second_advection =
                qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, time + half_step, time + step, advection_options);
            if (!second_advection.ok())
            {
                std::cerr << "Second advection particle half-step failed: " << second_advection.message << '\n';
                return 1;
            }
            ensemble = std::move(second_advection.ensemble);
            time += step;

            if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
            {
                final_moments = compute_ensemble_moments(ensemble, particle_weight, moment_extents);
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
            final_moments = compute_ensemble_moments(ensemble, particle_weight, moment_extents);
        }
        print_monte_carlo_result_summary(initial_moments, final_moments, terminal_moment_indices, output_options.directory, ensemble.size(),
                                         initialization_sample_count, initialization_seed, time_step, final_time, output_options.write_moments);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
