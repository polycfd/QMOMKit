#pragma once

#include "qmomkit/cqmom.hpp"
#include "qmomkit/quadrature.hpp"
#include "qmomkit/runge_kutta.hpp"

#include <cstddef>
#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace qmomkit
{

/**
 * Phase-space velocity callback, evaluated once for all quadrature nodes.
 *
 * `coordinates` holds `nodes * dimensions` values packed node-major, so node
 * `i` occupies `[i * dimensions, (i + 1) * dimensions)` and within a node the
 * values are ordered by phase-space dimension. `velocities` is supplied with
 * the same size and must be filled with dxi/dt in the same layout.
 *
 * The callback must not resize `velocities`; the library sizes it and reports
 * a size mismatch as a velocity failure. Every entry must be written: a
 * partially filled buffer is a contract violation that debug builds report as
 * a velocity failure, and that release builds cannot detect. No physical
 * meaning is assigned to any coordinate by the library.
 *
 * Evaluating all nodes together keeps the cost of one right-hand-side
 * evaluation independent of the terminal node count, which is the product of
 * the per-dimension node counts.
 */
using PhaseSpaceVelocity = std::function<void(double, const std::vector<double>&, std::vector<double>&, std::size_t)>;

/**
 * Convenience callback shape evaluating one quadrature node at a time.
 *
 * The coordinate vector is ordered by phase-space dimension. The returned
 * vector must have the same size and represents dxi/dt at the supplied time.
 *
 * Convert with `make_phase_space_velocity` where a `PhaseSpaceVelocity` is
 * required. The conversion calls this callback once per node and therefore
 * costs one call and one returned vector per node per right-hand-side
 * evaluation; prefer writing a `PhaseSpaceVelocity` directly when advancing
 * many nodes.
 */
using PerNodePhaseSpaceVelocity = std::function<std::vector<double>(double, const std::vector<double>&)>;

/**
 * Build a `PhaseSpaceVelocity` from a per-node callback.
 *
 * The result evaluates `velocity` once per node in ascending node order and
 * preserves its per-node validation, so a wrong-sized or non-finite per-node
 * return is reported exactly as the per-node contract describes. An empty
 * callback converts to an empty callback.
 */
[[nodiscard]] PhaseSpaceVelocity make_phase_space_velocity(PerNodePhaseSpaceVelocity velocity);

enum class NodeAdvectionStatus
{
    success,
    invalid_input,
    velocity_failed,
    integration_failed
};

enum class NodeAdvectionScheme
{
    rk23,
    rk45,
    ssp_rk3,
};

struct NodeAdvectionOptions
{
    NodeAdvectionScheme scheme = NodeAdvectionScheme::rk23;
    AdaptiveRkOptions adaptive_options{};
};

struct NodeAdvectionDiagnostics
{
    NodeAdvectionScheme scheme = NodeAdvectionScheme::rk23;
    std::size_t step_count = 0;
    std::size_t rejected_step_count = 0;
    std::size_t rhs_evaluations = 0;
    double final_time = 0.0;
    double last_step = 0.0;
    std::optional<double> estimated_error;
    // The status of the RK23 or RK45 integration. Empty for SSP-RK3, and empty
    // when the integration did not complete because the velocity callback
    // threw; `rhs_evaluations` then counts the evaluations made before that.
    std::optional<AdaptiveRkStatus> adaptive_status;
    std::string message;
};

struct QmomNodeAdvectionResult
{
    Quadrature1D quadrature;
    NodeAdvectionStatus status = NodeAdvectionStatus::success;
    NodeAdvectionDiagnostics diagnostics;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == NodeAdvectionStatus::success;
    }
};

struct CqmomNodeAdvectionResult
{
    CqmomQuadrature quadrature;
    NodeAdvectionStatus status = NodeAdvectionStatus::success;
    NodeAdvectionDiagnostics diagnostics;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == NodeAdvectionStatus::success;
    }
};

/**
 * Advanced composable operation: advance one-dimensional QMOM abscissas
 * through a phase-space velocity field. Weights are unchanged.
 */
[[nodiscard]] QmomNodeAdvectionResult advance_nodes(const Quadrature1D& quadrature, double t, double time_step, const PhaseSpaceVelocity& velocity,
                                                    const NodeAdvectionOptions& options = {});

/**
 * Advanced composable operation: advance one-dimensional public QMOM output
 * through a phase-space velocity field. Weights are unchanged.
 */
[[nodiscard]] QmomNodeAdvectionResult advance_nodes(const QmomResult& quadrature, double t, double time_step, const PhaseSpaceVelocity& velocity,
                                                    const NodeAdvectionOptions& options = {});

/**
 * Advanced composable operation: advance multidimensional CQMOM terminal
 * abscissas through a phase-space velocity field. Weights and
 * terminal paths are unchanged. Fixed rectangular and adaptive ragged
 * quadratures are both represented through their terminal-node arrays.
 */
[[nodiscard]] CqmomNodeAdvectionResult advance_nodes(const CqmomQuadrature& quadrature, double t, double time_step,
                                                     const PhaseSpaceVelocity& velocity, const NodeAdvectionOptions& options = {});

/**
 * Format effective node-advection settings for startup logs or diagnostics.
 */
[[nodiscard]] std::string format_node_advection_options(const NodeAdvectionOptions& options);

} // namespace qmomkit
