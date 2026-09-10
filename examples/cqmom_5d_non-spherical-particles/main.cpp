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

        auto configured = qmomkit::config::load_cqmom(configuration_file);
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
        const auto sampled = sample_initial_distribution(std::move(initial_condition).value(), configured.value().hierarchy());
        if (!sampled.has_value())
        {
            return 1;
        }
        if (!sampled->ok())
        {
            std::cerr << sampled->sampling.message << sampled->moments.message << '\n';
            return 1;
        }

        auto sampled_case = qmomkit::config::with_initial_moments(configured.value(), std::move(sampled->moments.moments));
        if (!sampled_case)
        {
            print_configuration_errors(sampled_case.error());
            return 1;
        }
        qmomkit::CqmomCase cqmom_case = std::move(sampled_case).value();
        const qmomkit::CqmomResult initial_cqmom = qmomkit::compute_cqmom(cqmom_case);
        if (!initial_cqmom.ok())
        {
            std::cerr << "CQMOM inversion failed: " << initial_cqmom.message << '\n';
            return 1;
        }

        const ParticleModelParameters particle_parameters = particle_model_parameters();
        validate_particle_model_parameters(particle_parameters);
        const BinaryFragmentationParameters breakage_parameters = binary_fragmentation_parameters();
        validate_binary_fragmentation_parameters(breakage_parameters);
        const qmomkit::PhaseSpaceVelocity velocity = particle_phase_space_velocity(particle_parameters);
        const qmomkit::CqmomMomentSource source = binary_breakage_source(breakage_parameters);

        const qmomkit::OutputOptions output_options = std::move(output_configured).value();
        const qmomkit::CqmomPhaseSpaceEvolutionOptions evolution_options = std::move(evolution_configured).value();
        const qmomkit::CqmomMomentLayout moment_layout = qmomkit::create_cqmom_moment_layout(cqmom_case.hierarchy());
        const qmomkit::RawMomentRequest terminal_moment_request{select_terminal_moment_indices(moment_layout.required_indices)};

        std::cout << qmomkit::format_cqmom_settings(cqmom_case.moment_storage(), cqmom_case.adaptive_branch_counts(), cqmom_case.vandermonde_solver(),
                                                    cqmom_case.diagnostic_level());
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

        qmomkit::CqmomState state = cqmom_case.initial_state();
        qmomkit::CqmomQuadrature quadrature = initial_cqmom.quadrature;
        if (output_options.write_quadrature || output_options.write_moments)
        {
            ensure_output_directory(output_options.directory);
        }
        std::optional<qmomkit::MomentWriter> moment_history_writer;
        if (output_options.write_moments)
        {
            qmomkit::MomentWriter::for_cqmom(output_options.directory / "initial_moments.txt", state).write(initial_time, state);
            moment_history_writer.emplace(qmomkit::MomentWriter::for_cqmom(output_options.directory / "moment_history.txt", state));
            moment_history_writer->write(initial_time, state);
        }
        std::vector<std::pair<double, qmomkit::CqmomQuadrature>> quadrature_history{{initial_time, quadrature}};
        const double output_interval = output_options.interval.value_or(final_time - initial_time);
        double next_output_time = initial_time + output_interval;
        double time = initial_time;

        while (time < final_time - 1.0e-12)
        {
            const double step = std::min(time_step, final_time - time);
            const qmomkit::CqmomPhaseSpaceEvolutionResult evolved =
                qmomkit::advance_cqmom_phase_space(state, time, step, velocity, source, evolution_options);
            if (!evolved.ok())
            {
                std::cerr << "CQMOM phase-space evolution failed: " << evolved.message << '\n';
                return 1;
            }
            if (!evolved.quadrature.has_value())
            {
                std::cerr << "CQMOM phase-space evolution did not return a quadrature.\n";
                return 1;
            }

            state = evolved.state;
            quadrature = *evolved.quadrature;
            time = evolved.final_time;
            if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
            {
                quadrature_history.emplace_back(time, quadrature);
                if (moment_history_writer.has_value())
                {
                    moment_history_writer->write(time, state);
                }
                while (next_output_time <= time + 1.0e-12)
                {
                    next_output_time += output_interval;
                }
            }
        }

        if (output_options.write_quadrature)
        {
            qmomkit::QuadratureWriter::for_cqmom(output_options.directory / "initial_quadrature.txt", initial_cqmom.quadrature.number_of_dimensions())
                .write(initial_time, initial_cqmom.quadrature);
            qmomkit::QuadratureWriter writer = qmomkit::QuadratureWriter::for_cqmom(output_options.directory / "quadrature_history.txt",
                                                                                    quadrature_history.front().second.number_of_dimensions());
            for (const auto& [history_time, history_quadrature] : quadrature_history)
            {
                writer.write(history_time, history_quadrature);
            }
        }

        const qmomkit::SelectedRawMomentResult initial_moments = qmomkit::reconstruct_raw_moments(initial_cqmom.quadrature, terminal_moment_request);
        if (!initial_moments.ok())
        {
            throw std::runtime_error("Initial CQMOM moment reconstruction failed: " + initial_moments.message);
        }
        const qmomkit::SelectedRawMomentResult final_moments = qmomkit::reconstruct_raw_moments(quadrature, terminal_moment_request);
        if (!final_moments.ok())
        {
            throw std::runtime_error("Final CQMOM moment reconstruction failed: " + final_moments.message);
        }
        print_result_summary(initial_moments, final_moments, quadrature, output_options.directory, initialization_sample_count, initialization_seed,
                             time_step, final_time, output_options.write_moments, output_options.write_quadrature);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
