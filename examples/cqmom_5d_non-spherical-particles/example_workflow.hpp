#pragma once

#include "breakage_model.hpp"
#include "particle_model.hpp"
#include "qmomkit/qmomkit.hpp"

#include <filesystem>
#include <optional>
#include <utility>
#include <vector>

namespace cqmom_5d_non_spherical_particles
{

constexpr double cqmom_initial_time = 0.0;
constexpr double cqmom_final_time = 1.0;
constexpr double cqmom_time_step = 0.005;

/**
 * Shared setup for the bounded and intentionally unbounded CQMOM examples.
 *
 * Both workflows load their own `.qmom` file, then replace the distributional
 * `[initial]` section with raw moments from the same fixed-seed finite ensemble.
 * This keeps the comparison focused on the chi support declaration rather than
 * on sampling variation.
 */
struct PreparedCqmomExample
{
    std::filesystem::path configuration_file;
    ParticleModelParameters particle_parameters;
    BinaryFragmentationParameters breakage_parameters;
    qmomkit::CqmomCase cqmom_case;
    qmomkit::OutputOptions output_options;
    qmomkit::CqmomPhaseSpaceEvolutionOptions phase_space_evolution_options;
    qmomkit::PhaseSpaceVelocity velocity;
    qmomkit::CqmomMomentSource source;
    qmomkit::CqmomState state;
    qmomkit::CqmomQuadrature initial_quadrature;
    qmomkit::CqmomMomentLayout moment_output_layout;
    qmomkit::RawMomentRequest file_moment_request;
    qmomkit::RawMomentRequest terminal_moment_request;
};

/**
 * Load a CQMOM example configuration, sample the configured initialization with
 * the example seed, and build the initial state, quadrature, phase-space evolution callback,
 * and breakage source used by the CQMOM drivers.
 *
 * Returns `std::nullopt` after printing user-facing configuration or
 * initialization errors.
 */
[[nodiscard]] std::optional<PreparedCqmomExample> prepare_cqmom_example(const std::filesystem::path& configuration_file);

/** Print the selected CQMOM reconstruction and phase-space evolution settings. */
void print_phase_space_evolution_settings(const PreparedCqmomExample& example);

/**
 * Print compact, deduplicated support decisions for the normal bounded run.
 *
 * The ordinary example keeps the terminal focused on topology changes only; the
 * full determinant/eigenvalue/Vandermonde explanation belongs to the unbounded
 * diagnostic report.
 */
void print_bounded_support_event_summaries(const qmomkit::PhaseSpaceEvolutionDiagnostics& diagnostics, double time);

/** Write the configured CQMOM output files and the final console summary. */
void write_cqmom_outputs_and_summary(const PreparedCqmomExample& example,
                                     const std::vector<std::pair<double, qmomkit::CqmomQuadrature>>& quadrature_history,
                                     const qmomkit::CqmomQuadrature& final_quadrature, double final_time);

} // namespace cqmom_5d_non_spherical_particles
