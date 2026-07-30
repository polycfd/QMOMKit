#pragma once

#include "qmomkit/cqmom.hpp"
#include "qmomkit/moment_reconstruction.hpp"
#include "qmomkit/moment_source.hpp"
#include "qmomkit/node_advection.hpp"
#include "qmomkit/qmom.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace qmomkit
{

/**
 * Operator-splitting scheme for one transport macro step.
 *
 * `strang`: Strang splitting, continuous(dt/2) -> discontinuous(dt) -> continuous(dt/2)
 * `cd`: continuous-then-discontinuous Lie splitting, continuous(dt) -> discontinuous(dt)
 * `dc`: discontinuous-then-continuous Lie splitting, discontinuous(dt) -> continuous(dt)
 *
 * Continuous transport advects quadrature nodes through phase-space velocity.
 * Discontinuous transport advances raw moments through MomentSource plus any
 * additive MomentBoundaryContribution supplied by the caller.
 */
enum class TransportScheme
{
    strang,
    cd,
    dc
};

enum class TransportStage
{
    none,
    initial_inversion,
    first_continuous_half_step,
    first_moment_reconstruction,
    discontinuous_source_step,
    post_source_inversion,
    second_continuous_half_step,
    final_moment_reconstruction
};

enum class TransportStatus
{
    success,
    invalid_input,
    inversion_failed,
    continuous_failed,
    moment_reconstruction_failed,
    source_failed
};

enum class CqmomNodeCountTransitionKind
{
    reduced,
    grew,
    created,
    removed
};

struct MomentReconstructionDiagnostics
{
    bool success = false;
    std::size_t moment_count = 0;
    std::string message;
};

struct CqmomNodeCountTransition
{
    TransportStage stage = TransportStage::none;
    std::size_t conditioning_level = 0;
    ConditioningPath parent_path;
    std::size_t phase_space_dimension = 0;
    std::size_t configured_maximum_node_count = 0;
    std::optional<std::size_t> active_node_count_before;
    std::optional<std::size_t> active_node_count_after;
    CqmomNodeCountTransitionKind kind = CqmomNodeCountTransitionKind::grew;
    bool rank_reduced = false;
    bool rank_grew = false;
};

struct TransportDiagnostics
{
    TransportScheme scheme = TransportScheme::strang;
    TransportStage failed_stage = TransportStage::none;
    CqmomMomentStorage cqmom_moment_storage = CqmomMomentStorage::lean;
    CqmomDiagnosticLevel cqmom_diagnostic_level = CqmomDiagnosticLevel::detailed;
    VandermondeSolver cqmom_vandermonde_solver = VandermondeSolver::automatic;
    QmomStatus post_source_qmom_status = QmomStatus::success;
    CqmomStatus initial_cqmom_status = CqmomStatus::success;
    CqmomStatus post_source_cqmom_status = CqmomStatus::success;
    std::vector<CqmomDiagnostics> initial_cqmom_diagnostics;
    std::vector<CqmomDiagnostics> post_source_cqmom_diagnostics;
    NodeAdvectionStatus first_continuous_status = NodeAdvectionStatus::success;
    NodeAdvectionStatus second_continuous_status = NodeAdvectionStatus::success;
    NodeAdvectionDiagnostics first_continuous_diagnostics;
    NodeAdvectionDiagnostics second_continuous_diagnostics;
    MomentReconstructionDiagnostics first_reconstruction_diagnostics;
    MomentReconstructionDiagnostics final_reconstruction_diagnostics;
    MomentSourceUpdateStatus source_status = MomentSourceUpdateStatus::success;
    MomentSourceUpdateDiagnostics source_diagnostics;
    std::vector<CqmomNodeCountTransition> cqmom_node_count_transitions;
    std::size_t total_node_advection_steps = 0;
    std::size_t total_node_advection_rejected_steps = 0;
    std::size_t total_node_advection_rhs_evaluations = 0;
    std::size_t source_retries = 0;
    std::string message;
};

struct QmomTransportOptions
{
    QmomOptions qmom_options;
    NodeAdvectionOptions node_options;
    MomentSourceUpdateOptions source_options;
    TransportScheme scheme = TransportScheme::strang;
};

struct CqmomTransportOptions
{
    NodeAdvectionOptions node_options;
    MomentSourceUpdateOptions source_options;
    TransportScheme scheme = TransportScheme::strang;
};

struct QmomTransportResult
{
    QmomState state;
    TransportStatus status = TransportStatus::success;
    TransportDiagnostics diagnostics;
    double final_time = 0.0;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == TransportStatus::success;
    }
};

struct CqmomTransportResult
{
    CqmomState state;
    std::optional<CqmomQuadrature> quadrature;
    TransportStatus status = TransportStatus::success;
    TransportDiagnostics diagnostics;
    double final_time = 0.0;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == TransportStatus::success;
    }
};

/**
 * Advance an accepted one-dimensional QMOM state by one transport macro step.
 *
 * Continuous transport advects quadrature nodes through the phase-space
 * velocity field, then reconstructs raw moments. Discontinuous transport
 * advances raw moments through MomentSource plus any additive boundary
 * contribution supplied by a boundary-aware overload.
 *
 * TransportScheme::strang applies second-order Strang splitting:
 * continuous(dt/2) -> discontinuous(dt) -> continuous(dt/2).
 * TransportScheme::cd applies first-order continuous-then-discontinuous Lie
 * splitting: continuous(dt) -> discontinuous(dt). TransportScheme::dc applies
 * first-order discontinuous-then-continuous Lie splitting:
 * discontinuous(dt) -> continuous(dt). In first-order schemes, diagnostic
 * "half-step" stage names refer to the
 * first/only continuous substep for backward compatibility.
 *
 * `maximum_node_count` is the configured fixed-count size or adaptive ceiling.
 * In adaptive mode the accepted state's quadrature size is the active count.
 */
[[nodiscard]] QmomTransportResult advance_transport_qmom(const QmomState& state, std::size_t maximum_node_count, double t, double time_step,
                                                         const PhaseSpaceVelocity& velocity, const MomentSource& source,
                                                         const QmomTransportOptions& options = {});

/**
 * Advance an accepted one-dimensional QMOM state using a generic moment source
 * plus an additive phase-space boundary exchange
 * contribution. The source and boundary rates are combined in the
 * discontinuous moment update; boundary exchange does not add a separate
 * transport stage or alter the selected operator-splitting scheme.
 */
[[nodiscard]] QmomTransportResult advance_transport_qmom(const QmomState& state, std::size_t maximum_node_count, double t, double time_step,
                                                         const PhaseSpaceVelocity& velocity, const MomentSource& source,
                                                         const MomentBoundaryContribution& boundary, const QmomTransportOptions& options);

/**
 * Advance a CQMOM transport state by one transport macro step.
 *
 * Reconstruction policies are taken from the state. The options configure only
 * transport-time controls such as node advection, source integration, and
 * the operator-splitting scheme. Continuous transport advects CQMOM quadrature
 * nodes and reconstructs state moments; discontinuous transport advances
 * moments through CqmomMomentSource plus any boundary contribution supplied by
 * a boundary-aware overload.
 *
 * TransportScheme::strang applies second-order Strang splitting:
 * continuous(dt/2) -> discontinuous(dt) -> continuous(dt/2).
 * TransportScheme::cd applies first-order continuous-then-discontinuous Lie
 * splitting: continuous(dt) -> discontinuous(dt). TransportScheme::dc applies
 * first-order discontinuous-then-continuous Lie splitting:
 * discontinuous(dt) -> continuous(dt).
 *
 * In adaptive branch-count mode, `CqmomState` carries current branch-specific
 * active counts separately from the hierarchy's configured maximum counts.
 * Continuous stages preserve the active topology, while accepted source stages
 * may grow active branches one node at a time toward the configured maximum.
 */
[[nodiscard]] CqmomTransportResult advance_transport_cqmom(const CqmomState& state, double t, double time_step, const PhaseSpaceVelocity& velocity,
                                                           const CqmomMomentSource& source, const CqmomTransportOptions& options = {});

/**
 * Advance a user-facing CQMOM runtime state using a CQMOM-native source
 * callback plus an additive phase-space boundary exchange contribution. The
 * callbacks remain runtime physics inputs and are not stored in the state. The
 * source and boundary rates are combined in the discontinuous moment update.
 */
[[nodiscard]] CqmomTransportResult advance_transport_cqmom(const CqmomState& state, double t, double time_step, const PhaseSpaceVelocity& velocity,
                                                           const CqmomMomentSource& source, const CqmomMomentBoundaryContribution& boundary,
                                                           const CqmomTransportOptions& options);

} // namespace qmomkit
