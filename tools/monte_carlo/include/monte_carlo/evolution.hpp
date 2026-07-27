#pragma once

#include "monte_carlo/ensemble.hpp"
#include "monte_carlo/explicit_distribution.hpp"
#include "monte_carlo/status.hpp"
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
 * Result of continuous sample-coordinate evolution.
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
    qmomkit::AdaptiveRk23Status rk23_status = qmomkit::AdaptiveRk23Status::success;
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
 * Result of evolving the nodes of an explicit weighted quadrature.
 *
 * The weights are unchanged by continuous phase-space motion. The evolved
 * distribution therefore stores the advanced node coordinates with the original
 * absolute weights and number-density convention.
 */
struct ExplicitQuadratureEvolutionResult
{
    ExplicitDiscreteDistribution initial_distribution;
    ExplicitDiscreteDistribution evolved_distribution;
    ComparisonStatus status = ComparisonStatus::success;
    std::size_t failed_node_index = static_cast<std::size_t>(-1);
    qmomkit::AdaptiveRk23Status rk23_status = qmomkit::AdaptiveRk23Status::success;
    std::string message;
    std::size_t total_accepted_steps = 0;
    std::size_t total_rejected_steps = 0;
    std::size_t total_function_evaluations = 0;
    std::size_t maximum_rejected_steps = 0;
    double final_time = 0.0;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == ComparisonStatus::success;
    }
};

/**
 * Evolve each sample independently through the supplied phase-space velocity.
 *
 * The callback must return one derivative per phase-space coordinate. Samples
 * are treated as generic coordinate realizations; no moment source, boundary
 * contribution, or transport splitting is applied.
 */
[[nodiscard]] EnsembleEvolutionResult evolve_ensemble(const Ensemble& ensemble, const qmomkit::PhaseSpaceVelocity& velocity, double t0, double t1,
                                                      const qmomkit::AdaptiveRk23Options& options = {});

/**
 * Evolve one-dimensional quadrature nodes with QMOMKit node advection.
 *
 * This helper is used to compare continuous node motion with continuous sample
 * motion under the same velocity callback.
 */
[[nodiscard]] qmomkit::QmomNodeAdvectionResult evolve_quadrature_nodes(const qmomkit::Quadrature1D& quadrature,
                                                                       const qmomkit::PhaseSpaceVelocity& velocity, double t0, double t1,
                                                                       const qmomkit::AdaptiveRk23Options& options = {});

/**
 * Evolve the nodes of an explicit weighted quadrature distribution.
 *
 * The returned distribution is a deterministic reference for moment comparison:
 * raw moments are reconstructed from the evolved nodes and unchanged absolute
 * weights.
 */
[[nodiscard]] ExplicitQuadratureEvolutionResult evolve_explicit_quadrature(const ExplicitDiscreteDistribution& distribution,
                                                                           const qmomkit::PhaseSpaceVelocity& velocity, double t0, double t1,
                                                                           const qmomkit::AdaptiveRk23Options& options);

} // namespace qmomkit::monte_carlo
