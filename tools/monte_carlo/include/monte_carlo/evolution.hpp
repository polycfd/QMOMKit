#pragma once

#include "monte_carlo/ensemble.hpp"
#include "qmomkit/node_advection.hpp"
#include "qmomkit/runge_kutta.hpp"

#include <cstddef>
#include <limits>
#include <string>

namespace qmomkit::monte_carlo
{

inline constexpr std::size_t invalid_sample_index = std::numeric_limits<std::size_t>::max();

/**
 * Status for independently evolving every sample in a finite ensemble.
 */
enum class EnsembleEvolutionStatus
{
    success,
    invalid_input,
    velocity_failed,
    rk23_failed
};

/**
 * Result of advection sample-coordinate evolution.
 *
 * On success, `ensemble` contains one evolved sample for each input sample. On
 * failure, no partially evolved ensemble is returned; diagnostics identify the
 * failing sample and RK23 status when available.
 */
struct EnsembleEvolutionResult
{
    Ensemble ensemble;
    EnsembleEvolutionStatus status = EnsembleEvolutionStatus::success;
    std::size_t failed_sample_index = invalid_sample_index;
    qmomkit::AdaptiveRkStatus rk23_status = qmomkit::AdaptiveRkStatus::success;
    std::string message;
    // Totals accumulated over completed samples and, for RK23 failures, the
    // failed sample's reported counters.
    std::size_t total_accepted_steps = 0;
    std::size_t total_rejected_steps = 0;
    std::size_t total_function_evaluations = 0;
    // Maximum rejected-step count reported by any individual sample.
    std::size_t maximum_rejected_steps = 0;
    // Requested final time on success; failure time for the failing sample.
    double final_time = 0.0;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == EnsembleEvolutionStatus::success;
    }
};

/**
 * Evolve each sample independently through the supplied phase-space velocity.
 *
 * Each sample is presented to the callback as a single packed node, so
 * `coordinates` holds that sample's phase-space coordinates and `velocities`
 * must be filled with one derivative per coordinate. Leaving an entry unwritten
 * is a contract violation that debug builds report as a velocity failure, and
 * that release builds cannot detect. Samples are treated as generic coordinate
 * realizations; no moment source, boundary contribution, or operator-split
 * phase-space evolution is applied.
 */
[[nodiscard]] EnsembleEvolutionResult evolve_ensemble(const Ensemble& ensemble, const qmomkit::PhaseSpaceVelocity& velocity, double t0, double t1,
                                                      const qmomkit::AdaptiveRkOptions& options = {});

} // namespace qmomkit::monte_carlo
