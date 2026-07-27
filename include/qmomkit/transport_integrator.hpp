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
 * Operator ordering for one transport macro step.
 *
 * `strang`: continuous(dt/2) -> discontinuous(dt) -> continuous(dt/2)
 * `cd`: continuous(dt) -> discontinuous(dt)
 * `dc`: discontinuous(dt) -> continuous(dt)
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

struct MomentReconstructionDiagnostics
{
    bool success = false;
    std::size_t moment_count = 0;
    std::string message;
};

struct TransportDiagnostics
{
    TransportScheme scheme = TransportScheme::strang;
    TransportStage failed_stage = TransportStage::none;
    QmomStatus initial_qmom_status = QmomStatus::success;
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

struct CqmomDenseTransportOptions
{
    ConditionalMomentOptions conditional_options;
    QmomOptions qmom_options;
    NodeAdvectionOptions node_options;
    MomentSourceUpdateOptions source_options;
    bool adaptive_branch_counts = false;
    TransportScheme scheme = TransportScheme::strang;
};

/**
 * Lean CQMOM transport state.
 *
 * The state stores only the raw moments required by the fixed recursive CQMOM
 * layout. Quadrature is reconstructed as a transport-stage artifact and is not
 * cached in the state.
 */
struct CqmomLeanTransportState
{
    CqmomMomentData moments;
};

struct CqmomLeanTransportOptions
{
    ConditionalMomentOptions conditional_options;
    QmomOptions qmom_options;
    NodeAdvectionOptions node_options;
    MomentSourceUpdateOptions source_options;
    TransportScheme scheme = TransportScheme::strang;
};

struct CqmomStateTransportOptions
{
    NodeAdvectionOptions node_options;
    MomentSourceUpdateOptions source_options;
    TransportScheme scheme = TransportScheme::strang;
};

struct QmomTransportResult
{
    std::vector<double> moments;
    QmomResult quadrature;
    TransportStatus status = TransportStatus::success;
    TransportDiagnostics diagnostics;
    double final_time = 0.0;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == TransportStatus::success;
    }
};

struct CqmomDenseTransportResult
{
    RawMomentTensor moments;
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

struct CqmomLeanTransportResult
{
    CqmomLeanTransportState state;
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

struct CqmomStateTransportResult
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
 * Advance one-dimensional raw moments by one hybrid transport macro step.
 *
 * Continuous transport advects quadrature nodes through the phase-space
 * velocity field, then reconstructs raw moments. Discontinuous transport
 * advances raw moments through MomentSource plus any additive boundary
 * contribution supplied by a boundary-aware overload.
 *
 * TransportScheme::strang applies
 * continuous(dt/2) -> discontinuous(dt) -> continuous(dt/2).
 * TransportScheme::cd applies continuous(dt) -> discontinuous(dt).
 * TransportScheme::dc applies discontinuous(dt) -> continuous(dt). In
 * first-order schemes, diagnostic "half-step" stage names refer to the
 * first/only continuous substep for backward compatibility.
 */
[[nodiscard]] QmomTransportResult advance_transport_qmom(const std::vector<double>& moments, std::size_t node_count, double t, double time_step,
                                                         const PhaseSpaceVelocity& velocity, const MomentSource& source,
                                                         const QmomTransportOptions& options = {});

/**
 * Advance one-dimensional raw moments by one hybrid transport macro step using
 * a generic moment source plus an additive phase-space boundary exchange
 * contribution. The source and boundary rates are combined in the
 * discontinuous moment update; boundary exchange does not add a separate
 * transport stage or alter the selected TransportScheme ordering.
 */
[[nodiscard]] QmomTransportResult advance_transport_qmom(const std::vector<double>& moments, std::size_t node_count, double t, double time_step,
                                                         const PhaseSpaceVelocity& velocity, const MomentSource& source,
                                                         const MomentBoundaryContribution& boundary, const QmomTransportOptions& options);

/**
 * Advance dense multidimensional raw moments by one hybrid transport macro
 * step while preserving phase-space-dimension ordering.
 *
 * Continuous transport advects CQMOM quadrature nodes through phase-space
 * velocity, then reconstructs dense raw moments. Discontinuous transport
 * advances the dense raw-moment vector through MomentSource plus any additive
 * MomentBoundaryContribution supplied by a boundary-aware overload.
 *
 * TransportScheme::strang applies
 * continuous(dt/2) -> discontinuous(dt) -> continuous(dt/2).
 * TransportScheme::cd applies continuous(dt) -> discontinuous(dt).
 * TransportScheme::dc applies discontinuous(dt) -> continuous(dt). In
 * first-order schemes, diagnostic "half-step" stage names refer to the
 * first/only continuous substep for backward compatibility.
 */
[[nodiscard]] CqmomDenseTransportResult advance_transport_cqmom(const RawMomentTensor& moments, const PhaseSpaceHierarchy& hierarchy, double t,
                                                                double time_step, const PhaseSpaceVelocity& velocity, const MomentSource& source,
                                                                const CqmomDenseTransportOptions& options = {});

/**
 * Advance dense multidimensional raw moments by one hybrid transport macro step
 * using a generic moment source plus an additive phase-space boundary exchange
 * contribution. The source and boundary rates are combined in the
 * discontinuous moment update; boundary exchange does not add a separate
 * transport stage or alter the selected TransportScheme ordering. The dense
 * workflow intentionally owns full raw-tensor storage; representation-neutral
 * CQMOM workflows should prefer the CqmomState overload.
 */
[[nodiscard]] CqmomDenseTransportResult advance_transport_cqmom(const RawMomentTensor& moments, const PhaseSpaceHierarchy& hierarchy, double t,
                                                                double time_step, const PhaseSpaceVelocity& velocity, const MomentSource& source,
                                                                const MomentBoundaryContribution& boundary,
                                                                const CqmomDenseTransportOptions& options);

/**
 * Advance lean fixed-layout CQMOM moments by one hybrid transport macro step
 * using a CQMOM-native source callback.
 *
 * Continuous transport advects CQMOM quadrature nodes through phase-space
 * velocity, then reconstructs required fixed-layout moments. Discontinuous
 * transport advances those moments through CqmomMomentSource plus any additive
 * CqmomMomentBoundaryContribution supplied by a boundary-aware overload.
 *
 * TransportScheme::strang applies
 * continuous(dt/2) -> discontinuous(dt) -> continuous(dt/2).
 * TransportScheme::cd applies continuous(dt) -> discontinuous(dt).
 * TransportScheme::dc applies discontinuous(dt) -> continuous(dt). Adaptive
 * branch evolution is intentionally unsupported for the explicit lean
 * transport API.
 */
[[nodiscard]] CqmomLeanTransportResult advance_transport_cqmom(const CqmomLeanTransportState& state, double t, double time_step,
                                                               const PhaseSpaceVelocity& velocity, const CqmomMomentSource& source,
                                                               const CqmomLeanTransportOptions& options = {});

/**
 * Advance lean fixed-layout CQMOM moments using a CQMOM-native source callback
 * plus an additive phase-space boundary exchange contribution. Both callbacks
 * receive the representation-neutral `CqmomSourceContext`, and their rates are
 * combined in the discontinuous moment update.
 */
[[nodiscard]] CqmomLeanTransportResult advance_transport_cqmom(const CqmomLeanTransportState& state, double t, double time_step,
                                                               const PhaseSpaceVelocity& velocity, const CqmomMomentSource& source,
                                                               const CqmomMomentBoundaryContribution& boundary,
                                                               const CqmomLeanTransportOptions& options);

/**
 * Advance a user-facing CQMOM runtime state by one Hybrid Transport macro step.
 *
 * Reconstruction policies are taken from the state. The options configure only
 * transport-time controls such as node advection, source integration, and
 * operator-splitting scheme. Continuous transport advects CQMOM quadrature
 * nodes and reconstructs state moments; discontinuous transport advances
 * moments through CqmomMomentSource plus any boundary contribution supplied by
 * a boundary-aware overload.
 *
 * TransportScheme::strang applies
 * continuous(dt/2) -> discontinuous(dt) -> continuous(dt/2).
 * TransportScheme::cd applies continuous(dt) -> discontinuous(dt).
 * TransportScheme::dc applies discontinuous(dt) -> continuous(dt).
 */
[[nodiscard]] CqmomStateTransportResult advance_transport_cqmom(const CqmomState& state, double t, double time_step,
                                                                const PhaseSpaceVelocity& velocity, const CqmomMomentSource& source,
                                                                const CqmomStateTransportOptions& options = {});

/**
 * Advance a user-facing CQMOM runtime state using a CQMOM-native source
 * callback plus an additive phase-space boundary exchange contribution. The
 * callbacks remain runtime physics inputs and are not stored in the state. The
 * source and boundary rates are combined in the discontinuous moment update.
 */
[[nodiscard]] CqmomStateTransportResult advance_transport_cqmom(const CqmomState& state, double t, double time_step,
                                                                const PhaseSpaceVelocity& velocity, const CqmomMomentSource& source,
                                                                const CqmomMomentBoundaryContribution& boundary,
                                                                const CqmomStateTransportOptions& options);

/**
 * Advance lean fixed-layout CQMOM moments with a zero source callback.
 *
 * This is the same hybrid transport composition as the CQMOM-native source
 * overload, with the discontinuous moment update supplied by a zero
 * CqmomMomentSource.
 *
 * TransportScheme::strang applies
 * continuous(dt/2) -> discontinuous(dt) -> continuous(dt/2).
 * TransportScheme::cd applies continuous(dt) -> discontinuous(dt).
 * TransportScheme::dc applies discontinuous(dt) -> continuous(dt).
 *
 * Adaptive branch evolution is intentionally unsupported for the explicit lean
 * transport API.
 */
[[nodiscard]] CqmomLeanTransportResult advance_transport_cqmom(const CqmomLeanTransportState& state, double t, double time_step,
                                                               const PhaseSpaceVelocity& velocity, const CqmomLeanTransportOptions& options = {});

} // namespace qmomkit
