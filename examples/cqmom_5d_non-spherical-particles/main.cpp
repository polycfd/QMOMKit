#include "breakage_model.hpp"
#include "particle_model.hpp"
#include "qmomkit/qmomkit.hpp"
#include "sampling.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <stdexcept>
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

        // The phase space here is (d, chi, u, v, w) with three nodes per
        // dimension, so the quadrature carries 3^5 = 243 terminal nodes. The
        // callbacks below see all of them at once, packed node-major.
        auto configured = qmomkit::config::load_cqmom_prepared(configuration_file);
        if (!configured)
        {
            print_configuration_errors(configured.error());
            return 1;
        }
        auto output_configured = qmomkit::config::load_output_options(configuration_file);
        if (!output_configured)
        {
            print_configuration_errors(output_configured.error());
            return 1;
        }
        auto evolution_configured = qmomkit::config::load_cqmom_phase_space_evolution_options(configuration_file);
        if (!evolution_configured)
        {
            print_configuration_errors(evolution_configured.error());
            return 1;
        }

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
        // CQMOM consumes only the moments its hierarchy requires, so the sampler
        // is asked for exactly those indices rather than for a dense tensor
        // whose remaining entries the case never reads. The lean replacement
        // then seeds the case without that tensor ever being built. For the
        // shipped configuration this is 606 moments instead of 7776; see the
        // README for the figures that belong to that hierarchy.
        const qmomkit::CqmomMomentLayout moment_layout = qmomkit::create_cqmom_moment_layout(configured.value().hierarchy());
        const qmomkit::RawMomentRequest required_moments{moment_layout.required_indices.to_multi_indices()};
        const auto sampled = sample_requested_initial_moments(std::move(initial_condition).value(), required_moments);
        if (!sampled.has_value())
        {
            return 1;
        }
        if (!sampled->ok())
        {
            std::cerr << sampled->sampling.message << sampled->moments.message << '\n';
            return 1;
        }

        auto sampled_case_result = qmomkit::config::with_initial_moments(
            configured.value(),
            qmomkit::create_cqmom_moment_data(moment_layout, qmomkit::make_selected_raw_moment_result(required_moments, sampled->moments.values)));
        if (!sampled_case_result)
        {
            print_configuration_errors(sampled_case_result.error());
            return 1;
        }
        qmomkit::CqmomPreparedCase cqmom_prepared_case = std::move(sampled_case_result).value();
        const ParticleModelParameters particle_parameters = particle_model_parameters();
        validate_particle_model_parameters(particle_parameters);
        const BinaryFragmentationParameters breakage_parameters = binary_fragmentation_parameters();
        validate_binary_fragmentation_parameters(breakage_parameters);
        // The model parameters are validated before the callbacks are built, so a
        // bad coefficient is reported here rather than as a velocity failure deep
        // inside the first step.
        const qmomkit::PhaseSpaceVelocity velocity = particle_phase_space_velocity(particle_parameters);
        const qmomkit::CqmomMomentSource source = binary_breakage_source(breakage_parameters);

        const qmomkit::OutputOptions output_options = std::move(output_configured).value();
        const qmomkit::CqmomPhaseSpaceEvolutionOptions evolution_options = std::move(evolution_configured).value();
        // Replacement preserves the hierarchy, so moment_layout built above is
        // this case's layout too and is reused rather than rebuilt here. The
        // hierarchy requires more moments than this example reports, so the
        // summary asks for a chosen subset rather than the whole layout.
        const qmomkit::RawMomentRequest terminal_moment_request{select_terminal_moment_indices(moment_layout.required_indices)};

        std::cout << qmomkit::format_cqmom_settings(cqmom_prepared_case.adaptive_branch_counts(), cqmom_prepared_case.vandermonde_solver(),
                                                    cqmom_prepared_case.diagnostic_level());
        std::cout << "Phase-space evolution scheme: ";
        switch (evolution_options.scheme)
        {
        case qmomkit::PhaseSpaceEvolutionScheme::strang:
            std::cout << "strang";
            break;
        case qmomkit::PhaseSpaceEvolutionScheme::advsrc:
            std::cout << "advsrc";
            break;
        case qmomkit::PhaseSpaceEvolutionScheme::srcadv:
            std::cout << "srcadv";
            break;
        }
        std::cout << '\n' << qmomkit::format_node_advection_options(evolution_options.node_options);

        qmomkit::CqmomState state = cqmom_prepared_case.initial_state();
        const qmomkit::CqmomQuadrature initial_quadrature = state.quadrature();
        if (output_options.write_quadrature || output_options.write_moments)
        {
            ensure_output_directory(output_options.directory);
        }
        std::optional<qmomkit::MomentWriter> moment_history_writer;
        if (output_options.write_moments)
        {
            qmomkit::MomentWriter::for_cqmom(output_options.directory / "initial_moments.txt", state, output_options)
                .write(initial_time, state, state.quadrature());
            moment_history_writer.emplace(qmomkit::MomentWriter::for_cqmom(output_options.directory / "moment_history.txt", state, output_options));
            moment_history_writer->write(initial_time, state, state.quadrature());
        }
        std::optional<qmomkit::QuadratureWriter> quadrature_history_writer;
        if (output_options.write_quadrature)
        {
            qmomkit::QuadratureWriter::for_cqmom(output_options.directory / "initial_quadrature.txt", initial_quadrature.number_of_dimensions())
                .write(initial_time, initial_quadrature);
            quadrature_history_writer.emplace(
                qmomkit::QuadratureWriter::for_cqmom(output_options.directory / "quadrature_history.txt", initial_quadrature.number_of_dimensions()));
            quadrature_history_writer->write(initial_time, initial_quadrature);
        }
        const double output_interval = output_options.interval.value_or(final_time - initial_time);
        double next_output_time = initial_time + output_interval;
        double time = initial_time;

        while (time < final_time - 1.0e-12)
        {
            const double step = std::min(time_step, final_time - time);
            const qmomkit::CqmomPhaseSpaceEvolutionResult evolved =
                // One complete step. With 243 terminal nodes the packed velocity
                // form matters: the callback is evaluated once per stage for all
                // nodes together rather than once per node.
                qmomkit::advance_cqmom_phase_space(state, time, step, velocity, source, evolution_options);
            if (!evolved.ok())
            {
                std::cerr << "CQMOM phase-space evolution failed: " << evolved.message << '\n';
                return 1;
            }
            state = evolved.state;
            time = evolved.final_time;
            if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
            {
                if (quadrature_history_writer.has_value())
                {
                    quadrature_history_writer->write(time, state.quadrature());
                }
                if (moment_history_writer.has_value())
                {
                    moment_history_writer->write(time, state, state.quadrature());
                }
                while (next_output_time <= time + 1.0e-12)
                {
                    next_output_time += output_interval;
                }
            }
        }

        const qmomkit::SelectedRawMomentResult initial_moments = qmomkit::reconstruct_raw_moments(initial_quadrature, terminal_moment_request);
        if (!initial_moments.ok())
        {
            throw std::runtime_error("Initial CQMOM moment reconstruction failed: " + initial_moments.message);
        }
        const qmomkit::SelectedRawMomentResult final_moments = qmomkit::reconstruct_raw_moments(state.quadrature(), terminal_moment_request);
        if (!final_moments.ok())
        {
            throw std::runtime_error("Final CQMOM moment reconstruction failed: " + final_moments.message);
        }
        print_result_summary(initial_moments, final_moments, state.quadrature(), output_options.directory, initialization_sample_count,
                             initialization_seed, time_step, final_time, output_options.write_moments, output_options.write_quadrature);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
