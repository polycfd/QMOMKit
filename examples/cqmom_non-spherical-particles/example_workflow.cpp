#include "example_workflow.hpp"

#include "sampling.hpp"
#include "utilities.hpp"

#include <cstddef>
#include <iomanip>
#include <iostream>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <utility>

namespace cqmom_non_spherical_particles
{
namespace
{

void print_path(const qmomkit::ConditioningPath& path)
{
    std::cout << '[';
    for (std::size_t level = 0; level < path.size(); ++level)
    {
        if (level > 0)
        {
            std::cout << ',';
        }
        std::cout << path.node_index_at_level(level);
    }
    std::cout << ']';
}

std::string path_string(const qmomkit::ConditioningPath& path)
{
    std::ostringstream text;
    text << '[';
    for (std::size_t level = 0; level < path.size(); ++level)
    {
        if (level > 0)
        {
            text << ',';
        }
        text << path.node_index_at_level(level);
    }
    text << ']';
    return text.str();
}

std::string support_event_key(const qmomkit::CqmomDiagnostics& branch)
{
    return branch.bounded_candidate_attempt_type + ":" + std::to_string(branch.phase_space_dimension) + ":" + path_string(branch.parent_path) + ":" +
           std::to_string(branch.bounded_candidate_node_count) + ":" + std::to_string(branch.bounded_accepted_fallback_node_count);
}

void print_event_header(const char* label, const double time)
{
    std::cout << std::setprecision(12) << label << '\n';
    std::cout << "  time: " << time << '\n';
}

} // namespace

std::optional<PreparedCqmomExample> prepare_cqmom_example(const std::filesystem::path& configuration_file)
{
    ParticleModelParameters particle_parameters = particle_model_parameters();
    validate_particle_model_parameters(particle_parameters);
    BinaryFragmentationParameters breakage_parameters = binary_fragmentation_parameters();
    validate_binary_fragmentation_parameters(breakage_parameters);
    qmomkit::PhaseSpaceVelocity velocity = particle_phase_space_velocity(particle_parameters);

    auto configured = qmomkit::config::load_cqmom(configuration_file);
    if (!configured)
    {
        print_configuration_errors(configured.error());
        return std::nullopt;
    }
    auto output_configured = qmomkit::config::load_output_options(configuration_file);
    if (!output_configured)
    {
        print_configuration_errors(output_configured.error());
        return std::nullopt;
    }
    auto phase_space_evolution_configured = qmomkit::config::load_cqmom_phase_space_evolution_options(configuration_file);
    if (!phase_space_evolution_configured)
    {
        print_configuration_errors(phase_space_evolution_configured.error());
        return std::nullopt;
    }

    const std::optional<std::string> configuration_text = read_configuration_text(configuration_file);
    if (!configuration_text.has_value())
    {
        return std::nullopt;
    }
    auto initial_condition = qmomkit::config::parse_cqmom_initial_condition(*configuration_text, configuration_file.string());
    if (!initial_condition)
    {
        print_configuration_errors(initial_condition.error());
        return std::nullopt;
    }
    // The CQMOM workflows consume deterministic finite-sample moments rather than
    // analytical initialization moments so the bounded run, unbounded diagnostic,
    // and MC comparison start from one reproducible sampled realization.
    const std::optional<qmomkit::monte_carlo::TensorSampledInitialization> sampled =
        sample_initial_distribution(std::move(initial_condition).value(), configured.value().hierarchy());
    if (!sampled.has_value())
    {
        return std::nullopt;
    }
    if (!sampled->ok())
    {
        std::cerr << sampled->sampling.message << sampled->moments.message << '\n';
        return std::nullopt;
    }

    auto sampled_case_result = qmomkit::config::with_initial_moments(configured.value(), std::move(sampled->moments.moments));
    if (!sampled_case_result)
    {
        print_configuration_errors(sampled_case_result.error());
        return std::nullopt;
    }
    qmomkit::CqmomCase cqmom_case = std::move(sampled_case_result).value();
    const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(cqmom_case);
    if (!cqmom.ok())
    {
        std::cerr << "CQMOM inversion failed: " << cqmom.message << '\n';
        return std::nullopt;
    }

    qmomkit::CqmomMomentLayout moment_output_layout = qmomkit::create_cqmom_moment_layout(cqmom_case.hierarchy());
    qmomkit::RawMomentRequest file_moment_request{moment_output_layout.required_indices};
    const std::vector<qmomkit::MultiIndex> terminal_moment_indices = select_terminal_moment_indices(moment_output_layout.required_indices);
    qmomkit::RawMomentRequest terminal_moment_request{terminal_moment_indices};
    qmomkit::CqmomMomentSource source = binary_breakage_source(cqmom_case, breakage_parameters);
    qmomkit::CqmomState state = cqmom_case.initial_state();

    return PreparedCqmomExample{configuration_file,
                                particle_parameters,
                                breakage_parameters,
                                std::move(cqmom_case),
                                std::move(output_configured).value(),
                                std::move(phase_space_evolution_configured).value(),
                                std::move(velocity),
                                std::move(source),
                                std::move(state),
                                cqmom.quadrature,
                                std::move(moment_output_layout),
                                std::move(file_moment_request),
                                std::move(terminal_moment_request)};
}

void print_phase_space_evolution_settings(const PreparedCqmomExample& example)
{
    std::cout << qmomkit::format_cqmom_settings(example.cqmom_case.moment_storage(), example.cqmom_case.adaptive_branch_counts(),
                                                example.cqmom_case.vandermonde_solver(), example.cqmom_case.diagnostic_level());
    std::cout << "Phase-space evolution scheme: ";
    switch (example.phase_space_evolution_options.scheme)
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
    std::cout << '\n';
    std::cout << qmomkit::format_node_advection_options(example.phase_space_evolution_options.node_options);
}

void print_bounded_support_event_summaries(const qmomkit::PhaseSpaceEvolutionDiagnostics& diagnostics, const double time)
{
    static std::set<std::string> printed_events;
    for (const qmomkit::CqmomDiagnostics& branch : diagnostics.post_source_cqmom_diagnostics)
    {
        if (!branch.bounded_support_checked || !branch.bounded_support_interval_present || branch.bounded_accepted_fallback_node_count == 0 ||
            branch.bounded_accepted_fallback_node_count >= branch.bounded_candidate_node_count)
        {
            continue;
        }
        if (!printed_events.insert(support_event_key(branch)).second)
        {
            continue;
        }

        const std::size_t previous_count = branch.bounded_previous_active_node_count.value_or(branch.bounded_candidate_node_count);
        if (branch.bounded_candidate_attempt_type == "growth")
        {
            print_event_header("Support-infeasible growth rejected:", time);
            std::cout << "  dimension: " << branch.phase_space_dimension << '\n';
            std::cout << "  branch: ";
            print_path(branch.parent_path);
            std::cout << '\n';
            std::cout << "  preserved count: " << previous_count << '\n';
            std::cout << "  rejected growth count: " << branch.bounded_candidate_node_count << '\n';
        }
        else
        {
            print_event_header("Support-driven active-count reduction:", time);
            std::cout << "  dimension: " << branch.phase_space_dimension << '\n';
            std::cout << "  branch: ";
            print_path(branch.parent_path);
            std::cout << '\n';
            std::cout << "  previous count: " << previous_count << '\n';
            std::cout << "  rejected count: " << branch.bounded_candidate_node_count << '\n';
            std::cout << "  accepted count: " << branch.bounded_accepted_fallback_node_count << '\n';
        }
        std::cout << "  interval: [" << branch.bounded_support_interval.lower << ", " << branch.bounded_support_interval.upper << "]\n";
    }
}

void write_cqmom_outputs_and_summary(const PreparedCqmomExample& example,
                                     const std::vector<std::pair<double, qmomkit::CqmomQuadrature>>& quadrature_history,
                                     const qmomkit::CqmomQuadrature& final_quadrature, const double final_time)
{
    if (example.output_options.write_quadrature || example.output_options.write_moments)
    {
        ensure_output_directory(example.output_options.directory);
    }
    if (example.output_options.write_moments)
    {
        qmomkit::MomentWriter::for_cqmom(example.output_options.directory / "initial_moments.txt", example.file_moment_request)
            .write(cqmom_initial_time, example.initial_quadrature);

        qmomkit::MomentWriter moment_history_writer =
            qmomkit::MomentWriter::for_cqmom(example.output_options.directory / "moment_history.txt", example.file_moment_request);
        for (const auto& [history_time, history_quadrature] : quadrature_history)
        {
            moment_history_writer.write(history_time, history_quadrature);
        }
    }
    if (example.output_options.write_quadrature)
    {
        qmomkit::QuadratureWriter::for_cqmom(example.output_options.directory / "initial_quadrature.txt",
                                             example.initial_quadrature.number_of_dimensions())
            .write(cqmom_initial_time, example.initial_quadrature);

        qmomkit::QuadratureWriter quadrature_history_writer = qmomkit::QuadratureWriter::for_cqmom(
            example.output_options.directory / "quadrature_history.txt", quadrature_history.front().second.number_of_dimensions());
        for (const auto& [history_time, history_quadrature] : quadrature_history)
        {
            quadrature_history_writer.write(history_time, history_quadrature);
        }
    }

    const qmomkit::SelectedRawMomentResult initial_moments =
        qmomkit::reconstruct_raw_moments(example.initial_quadrature, example.terminal_moment_request);
    if (!initial_moments.ok())
    {
        throw std::runtime_error("Initial CQMOM moment reconstruction failed: " + initial_moments.message);
    }
    const qmomkit::SelectedRawMomentResult final_moments = qmomkit::reconstruct_raw_moments(final_quadrature, example.terminal_moment_request);
    if (!final_moments.ok())
    {
        throw std::runtime_error("Final CQMOM moment reconstruction failed: " + final_moments.message);
    }
    print_result_summary(initial_moments, final_moments, final_quadrature, example.output_options.directory, initialization_sample_count,
                         initialization_seed, cqmom_time_step, final_time, example.output_options.write_moments,
                         example.output_options.write_quadrature);
}

} // namespace cqmom_non_spherical_particles
