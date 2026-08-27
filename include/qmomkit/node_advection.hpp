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
 * Phase-space velocity callback for one quadrature node.
 *
 * The coordinate vector is ordered by phase-space dimension. The returned
 * vector must have the same size and represents dxi/dt at the supplied time.
 * No physical meaning is assigned to any coordinate by the library.
 */
using PhaseSpaceVelocity = std::function<std::vector<double>(double, const std::vector<double>&)>;

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
    AdaptiveRkOptions adaptive_options;
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
