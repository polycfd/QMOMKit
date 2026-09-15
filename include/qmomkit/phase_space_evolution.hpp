#pragma once

#include "qmomkit/cqmom.hpp"
#include "qmomkit/moment_reconstruction.hpp"
#include "qmomkit/moment_source.hpp"
#include "qmomkit/node_advection.hpp"
#include "qmomkit/qmom.hpp"

#include <cstddef>
#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace qmomkit
{

/**
 * Operator-splitting scheme for one phase-space evolution macro step.
 *
 * `strang`: Strang splitting, advection(dt/2) -> moment-source(dt) -> advection(dt/2)
 * `advsrc`: advection-then-moment-source Lie splitting, advection(dt) -> moment-source(dt)
 * `srcadv`: moment-source-then-advection Lie splitting, moment-source(dt) -> advection(dt)
 *
 * Phase-space advection advects quadrature nodes through phase-space velocity.
 * Moment-source evolution advances raw moments through a `QmomMomentSource` plus
 * any additive `QmomMomentBoundaryContribution` supplied by the caller.
 */
enum class PhaseSpaceEvolutionScheme
{
    strang,
    advsrc,
    srcadv
};

enum class PhaseSpaceEvolutionStage
{
    none,
    initial_inversion,
    first_advection_half_step,
    first_moment_reconstruction,
    moment_source_step,
    post_source_inversion,
    second_advection_half_step,
    final_moment_reconstruction
};

enum class PhaseSpaceEvolutionStatus
{
    success,
    invalid_input,
    inversion_failed,
    advection_failed,
    moment_reconstruction_failed,
    moment_source_failed
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
    PhaseSpaceEvolutionStage stage = PhaseSpaceEvolutionStage::none;
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

struct PhaseSpaceEvolutionDiagnostics
{
    PhaseSpaceEvolutionScheme scheme = PhaseSpaceEvolutionScheme::strang;
    PhaseSpaceEvolutionStage failed_stage = PhaseSpaceEvolutionStage::none;
    CqmomMomentStorage cqmom_moment_storage = CqmomMomentStorage::lean;
    CqmomDiagnosticLevel cqmom_diagnostic_level = CqmomDiagnosticLevel::none;
    VandermondeSolver cqmom_vandermonde_solver = VandermondeSolver::automatic;
    QmomStatus post_source_qmom_status = QmomStatus::success;
    CqmomStatus initial_cqmom_status = CqmomStatus::success;
    CqmomStatus post_source_cqmom_status = CqmomStatus::success;
    std::vector<CqmomDiagnostics> initial_cqmom_diagnostics;
    std::vector<CqmomDiagnostics> post_source_cqmom_diagnostics;
    NodeAdvectionStatus first_advection_status = NodeAdvectionStatus::success;
    NodeAdvectionStatus second_advection_status = NodeAdvectionStatus::success;
    NodeAdvectionDiagnostics first_advection_diagnostics;
    NodeAdvectionDiagnostics second_advection_diagnostics;
    MomentReconstructionDiagnostics first_reconstruction_diagnostics;
    MomentReconstructionDiagnostics final_reconstruction_diagnostics;
    MomentSourceUpdateStatus moment_source_status = MomentSourceUpdateStatus::success;
    MomentSourceUpdateDiagnostics moment_source_diagnostics;
    std::vector<CqmomNodeCountTransition> cqmom_node_count_transitions;
    std::size_t total_node_advection_steps = 0;
    std::size_t total_node_advection_rejected_steps = 0;
    std::size_t total_node_advection_rhs_evaluations = 0;
    std::size_t moment_source_retries = 0;
    std::string message;
};

/** Per-call controls for one-dimensional QMOM phase-space evolution.
 *
 * Reconstruction options and the configured maximum node count are carried by
 * QmomState.
 */
struct QmomPhaseSpaceEvolutionOptions
{
    NodeAdvectionOptions node_options;
    MomentSourceUpdateOptions source_options;
    PhaseSpaceEvolutionScheme scheme = PhaseSpaceEvolutionScheme::strang;
    RealizabilityPolicy realizability_policy;
};

struct CqmomPhaseSpaceEvolutionOptions
{
    NodeAdvectionOptions node_options;
    MomentSourceUpdateOptions source_options;
    PhaseSpaceEvolutionScheme scheme = PhaseSpaceEvolutionScheme::strang;
};

struct QmomPhaseSpaceEvolutionResult
{
    QmomState state;
    PhaseSpaceEvolutionStatus status = PhaseSpaceEvolutionStatus::success;
    PhaseSpaceEvolutionDiagnostics diagnostics;
    double final_time = 0.0;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == PhaseSpaceEvolutionStatus::success;
    }
};

struct CqmomPhaseSpaceEvolutionResult
{
    CqmomState state;
    std::optional<CqmomQuadrature> quadrature;
    PhaseSpaceEvolutionStatus status = PhaseSpaceEvolutionStatus::success;
    PhaseSpaceEvolutionDiagnostics diagnostics;
    double final_time = 0.0;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == PhaseSpaceEvolutionStatus::success;
    }
};

/**
 * Batch-index-aware phase-space velocity callback.
 *
 * `batch_index` is the zero-based position of the state in the collection
 * supplied to a batch phase-space evolution operation. Coordinates and
 * velocities follow the packed `PhaseSpaceVelocity` layout.
 */
using BatchPhaseSpaceVelocity = std::function<void(std::size_t batch_index, double time, const std::vector<double>& coordinates,
                                                   std::vector<double>& velocities, std::size_t dimensions)>;

/**
 * Batch-index-aware one-dimensional moment-source callback.
 *
 * `batch_index` is the zero-based position of the state in the collection
 * supplied to a batch phase-space evolution operation.
 */
using BatchQmomMomentSource = std::function<std::vector<double>(std::size_t batch_index, double time, const QmomSourceContext& context)>;

/**
 * Batch-index-aware one-dimensional boundary-contribution callback.
 *
 * `batch_index` is the zero-based position of the state in the collection
 * supplied to a batch phase-space evolution operation.
 */
using BatchQmomMomentBoundaryContribution = std::function<QmomMomentRate(std::size_t batch_index, double time, const std::vector<double>& moments)>;

/**
 * Batch-index-aware CQMOM moment-source callback.
 *
 * `batch_index` is the zero-based position of the state in the collection
 * supplied to a batch phase-space evolution operation. The CQMOM source context remains
 * non-owning and valid only during the callback invocation.
 */
using BatchCqmomMomentSource = std::function<CqmomMomentRate(std::size_t batch_index, double time, const CqmomSourceContext& context)>;

/**
 * Batch-index-aware CQMOM boundary-contribution callback.
 *
 * `batch_index` is the zero-based position of the state in the collection
 * supplied to a batch phase-space evolution operation. The CQMOM source context remains
 * non-owning and valid only during the callback invocation.
 */
using BatchCqmomMomentBoundaryContribution = std::function<CqmomMomentRate(std::size_t batch_index, double time, const CqmomSourceContext& context)>;

/**
 * Advance an accepted one-dimensional QMOM state through one phase-space
 * advection stage.
 *
 * The stage advects the accepted quadrature, preserves the active node count,
 * does not perform QMOM inversion, and returns a canonical state whose full
 * moment vector is reconstructed from the advanced quadrature. On failure, the
 * returned state is the unchanged input state.
 */
[[nodiscard]] QmomPhaseSpaceEvolutionResult advance_qmom_advection(const QmomState& state, double t, double time_step,
                                                                   const PhaseSpaceVelocity& velocity,
                                                                   const QmomPhaseSpaceEvolutionOptions& options = {});

/**
 * Advance an accepted one-dimensional QMOM state through one moment-source
 * evolution stage.
 *
 * The stage advances the full stored moment vector, applies fixed or adaptive
 * post-source QMOM reconstruction, may reduce the active node count, and may
 * grow it by at most one node toward the state's configured maximum. On success it returns
 * a canonical accepted state; on failure, it returns the unchanged input state.
 */
[[nodiscard]] QmomPhaseSpaceEvolutionResult advance_qmom_moment_source(const QmomState& state, double t, double time_step,
                                                                       const QmomMomentSource& source,
                                                                       const QmomPhaseSpaceEvolutionOptions& options = {});

/**
 * Advance an accepted one-dimensional QMOM state through one moment-source
 * evolution stage using a generic moment source plus an additive phase-space
 * boundary exchange contribution.
 */
[[nodiscard]] QmomPhaseSpaceEvolutionResult advance_qmom_moment_source(const QmomState& state, double t, double time_step,
                                                                       const QmomMomentSource& source, const QmomMomentBoundaryContribution& boundary,
                                                                       const QmomPhaseSpaceEvolutionOptions& options);

/**
 * Advance an accepted one-dimensional QMOM state by one phase-space evolution
 * macro step.
 *
 * Phase-space advection advects quadrature nodes through the phase-space
 * velocity field, then reconstructs raw moments. Moment-source evolution
 * advances raw moments through a `QmomMomentSource` plus any additive boundary
 * contribution supplied by a boundary-aware overload.
 *
 * `PhaseSpaceEvolutionScheme::strang` applies second-order Strang splitting:
 * advection(dt/2) -> moment-source(dt) -> advection(dt/2).
 * `PhaseSpaceEvolutionScheme::advsrc` applies first-order advection-then-moment-source Lie
 * splitting: advection(dt) -> moment-source(dt). `PhaseSpaceEvolutionScheme::srcadv` applies
 * first-order moment-source-then-advection Lie splitting:
 * moment-source(dt) -> advection(dt). In first-order schemes, the single
 * advection substep's diagnostics land in the `first_*` fields for `advsrc`
 * and in the `second_*`/`final_*` fields for `srcadv`.
 *
 * The state's configured maximum node count and reconstruction options govern
 * reconstruction. In adaptive mode the accepted state's quadrature size is the
 * active count.
 */
[[nodiscard]] QmomPhaseSpaceEvolutionResult advance_qmom_phase_space(const QmomState& state, double t, double time_step,
                                                                     const PhaseSpaceVelocity& velocity, const QmomMomentSource& source,
                                                                     const QmomPhaseSpaceEvolutionOptions& options = {});

/**
 * Advance an accepted one-dimensional QMOM state using a generic moment source
 * plus an additive phase-space boundary exchange
 * contribution. The source and boundary rates are combined in the
 * moment-source update; boundary exchange does not add a separate
 * phase-space evolution stage or alter the selected operator-splitting scheme.
 */
[[nodiscard]] QmomPhaseSpaceEvolutionResult advance_qmom_phase_space(const QmomState& state, double t, double time_step,
                                                                     const PhaseSpaceVelocity& velocity, const QmomMomentSource& source,
                                                                     const QmomMomentBoundaryContribution& boundary,
                                                                     const QmomPhaseSpaceEvolutionOptions& options);

/**
 * Advance an accepted CQMOM state through one phase-space advection stage.
 *
 * The stage advects accepted terminal quadrature nodes, preserves terminal
 * topology and active branch counts, and reconstructs the state moments from
 * the advanced quadrature without a post-stage CQMOM inversion. On failure, the
 * returned state is the unchanged input state.
 */
[[nodiscard]] CqmomPhaseSpaceEvolutionResult advance_cqmom_advection(const CqmomState& state, double t, double time_step,
                                                                     const PhaseSpaceVelocity& velocity,
                                                                     const CqmomPhaseSpaceEvolutionOptions& options = {});

/**
 * Advance an accepted CQMOM state through one moment-source evolution stage.
 *
 * The stage advances the canonical CQMOM phase-space evolution moment layout, applies the
 * configured recursive CQMOM reconstruction policy, may reduce active branch
 * counts, and may grow retained adaptive branches by at most one node toward
 * the configured hierarchy maximum. On failure, the returned state is the
 * unchanged input state.
 */
[[nodiscard]] CqmomPhaseSpaceEvolutionResult advance_cqmom_moment_source(const CqmomState& state, double t, double time_step,
                                                                         const CqmomMomentSource& source,
                                                                         const CqmomPhaseSpaceEvolutionOptions& options = {});

/**
 * Advance an accepted CQMOM state through one moment-source evolution stage using
 * a CQMOM-native source plus an additive phase-space boundary exchange
 * contribution.
 */
[[nodiscard]] CqmomPhaseSpaceEvolutionResult advance_cqmom_moment_source(const CqmomState& state, double t, double time_step,
                                                                         const CqmomMomentSource& source,
                                                                         const CqmomMomentBoundaryContribution& boundary,
                                                                         const CqmomPhaseSpaceEvolutionOptions& options);

/**
 * Advance a CQMOM phase-space evolution state by one macro step.
 *
 * Reconstruction policies are taken from the state. The options configure only
 * phase-space evolution controls such as node advection, moment-source integration, and
 * the operator-splitting scheme. Phase-space advection advects CQMOM quadrature
 * nodes and reconstructs state moments; moment-source evolution advances
 * moments through a `CqmomMomentSource` plus any boundary contribution supplied by
 * a boundary-aware overload.
 *
 * `PhaseSpaceEvolutionScheme::strang` applies second-order Strang splitting:
 * advection(dt/2) -> moment-source(dt) -> advection(dt/2).
 * `PhaseSpaceEvolutionScheme::advsrc` applies first-order advection-then-moment-source Lie
 * splitting: advection(dt) -> moment-source(dt). `PhaseSpaceEvolutionScheme::srcadv` applies
 * first-order moment-source-then-advection Lie splitting:
 * moment-source(dt) -> advection(dt).
 *
 * In adaptive branch-count mode, `CqmomState` carries current branch-specific
 * active counts separately from the hierarchy's configured maximum counts.
 * Phase-space advection stages preserve the active topology, while accepted moment-source stages
 * may grow active branches one node at a time toward the configured maximum.
 */
[[nodiscard]] CqmomPhaseSpaceEvolutionResult advance_cqmom_phase_space(const CqmomState& state, double t, double time_step,
                                                                       const PhaseSpaceVelocity& velocity, const CqmomMomentSource& source,
                                                                       const CqmomPhaseSpaceEvolutionOptions& options = {});

/**
 * Advance a user-facing CQMOM runtime state using a CQMOM-native source
 * callback plus an additive phase-space boundary exchange contribution. The
 * callbacks remain runtime physics inputs and are not stored in the state. The
 * source and boundary rates are combined in the moment-source update.
 */
[[nodiscard]] CqmomPhaseSpaceEvolutionResult advance_cqmom_phase_space(const CqmomState& state, double t, double time_step,
                                                                       const PhaseSpaceVelocity& velocity, const CqmomMomentSource& source,
                                                                       const CqmomMomentBoundaryContribution& boundary,
                                                                       const CqmomPhaseSpaceEvolutionOptions& options);

/**
 * Serial batch convenience wrappers for accepted QMOM and CQMOM phase-space
 * evolution
 * states.
 *
 * Every wrapper reserves one result per input state and delegates in
 * increasing input order to the corresponding scalar public operation. Input
 * states are not modified, adaptive behavior remains independent per entry,
 * and an empty input returns an empty result collection without validation or
 * callback invocation. Time, timestep, and evolution options are shared by the
 * batch; each QMOM state supplies its own reconstruction metadata.
 *
 * Ordinary batch functions reuse one scalar callback for every entry.
 * Batch-index-aware functions have a `_batch_with_index` suffix and bind the
 * zero-based input position into each callback invocation; that position is not
 * a spatial or host-application identifier. Distinct names select the callback
 * family without callable-conversion ambiguity. QMOMKit does not provide
 * parallel execution or spatial organization through this API, and hosts may
 * instead call scalar operations inside their own parallel loops.
 *
 * An unsuccessful scalar result is appended unchanged. Exceptions that escape
 * a scalar call escape unchanged and stop processing; no partial result
 * collection is returned.
 */

// Ordinary-callback advection-stage batch wrappers.
[[nodiscard]] std::vector<QmomPhaseSpaceEvolutionResult> advance_qmom_advection_batch(const std::vector<QmomState>& states, double time,
                                                                                      double time_step, const PhaseSpaceVelocity& velocity,
                                                                                      const QmomPhaseSpaceEvolutionOptions& options = {});

[[nodiscard]] std::vector<CqmomPhaseSpaceEvolutionResult> advance_cqmom_advection_batch(const std::vector<CqmomState>& states, double time,
                                                                                        double time_step, const PhaseSpaceVelocity& velocity,
                                                                                        const CqmomPhaseSpaceEvolutionOptions& options = {});

// Ordinary-callback moment-source-evolution batch wrappers.
[[nodiscard]] std::vector<QmomPhaseSpaceEvolutionResult> advance_qmom_moment_source_batch(const std::vector<QmomState>& states, double time,
                                                                                          double time_step, const QmomMomentSource& source,
                                                                                          const QmomPhaseSpaceEvolutionOptions& options = {});

[[nodiscard]] std::vector<CqmomPhaseSpaceEvolutionResult> advance_cqmom_moment_source_batch(const std::vector<CqmomState>& states, double time,
                                                                                            double time_step, const CqmomMomentSource& source,
                                                                                            const CqmomPhaseSpaceEvolutionOptions& options = {});

[[nodiscard]] std::vector<QmomPhaseSpaceEvolutionResult> advance_qmom_moment_source_batch(const std::vector<QmomState>& states, double time,
                                                                                          double time_step, const QmomMomentSource& source,
                                                                                          const QmomMomentBoundaryContribution& boundary,
                                                                                          const QmomPhaseSpaceEvolutionOptions& options);

[[nodiscard]] std::vector<CqmomPhaseSpaceEvolutionResult> advance_cqmom_moment_source_batch(const std::vector<CqmomState>& states, double time,
                                                                                            double time_step, const CqmomMomentSource& source,
                                                                                            const CqmomMomentBoundaryContribution& boundary,
                                                                                            const CqmomPhaseSpaceEvolutionOptions& options);

// Ordinary-callback complete phase-space evolution batch wrappers.
[[nodiscard]] std::vector<QmomPhaseSpaceEvolutionResult> advance_qmom_phase_space_batch(const std::vector<QmomState>& states, double time,
                                                                                        double time_step, const PhaseSpaceVelocity& velocity,
                                                                                        const QmomMomentSource& source,
                                                                                        const QmomPhaseSpaceEvolutionOptions& options = {});

[[nodiscard]] std::vector<CqmomPhaseSpaceEvolutionResult> advance_cqmom_phase_space_batch(const std::vector<CqmomState>& states, double time,
                                                                                          double time_step, const PhaseSpaceVelocity& velocity,
                                                                                          const CqmomMomentSource& source,
                                                                                          const CqmomPhaseSpaceEvolutionOptions& options = {});

[[nodiscard]] std::vector<QmomPhaseSpaceEvolutionResult> advance_qmom_phase_space_batch(const std::vector<QmomState>& states, double time,
                                                                                        double time_step, const PhaseSpaceVelocity& velocity,
                                                                                        const QmomMomentSource& source,
                                                                                        const QmomMomentBoundaryContribution& boundary,
                                                                                        const QmomPhaseSpaceEvolutionOptions& options);

[[nodiscard]] std::vector<CqmomPhaseSpaceEvolutionResult> advance_cqmom_phase_space_batch(const std::vector<CqmomState>& states, double time,
                                                                                          double time_step, const PhaseSpaceVelocity& velocity,
                                                                                          const CqmomMomentSource& source,
                                                                                          const CqmomMomentBoundaryContribution& boundary,
                                                                                          const CqmomPhaseSpaceEvolutionOptions& options);

// Batch-index-aware advection-stage batch wrappers.
[[nodiscard]] std::vector<QmomPhaseSpaceEvolutionResult> advance_qmom_advection_batch_with_index(const std::vector<QmomState>& states, double time,
                                                                                                 double time_step,
                                                                                                 const BatchPhaseSpaceVelocity& velocity,
                                                                                                 const QmomPhaseSpaceEvolutionOptions& options = {});

[[nodiscard]] std::vector<CqmomPhaseSpaceEvolutionResult>
advance_cqmom_advection_batch_with_index(const std::vector<CqmomState>& states, double time, double time_step,
                                         const BatchPhaseSpaceVelocity& velocity, const CqmomPhaseSpaceEvolutionOptions& options = {});

// Batch-index-aware moment-source-evolution batch wrappers.
[[nodiscard]] std::vector<QmomPhaseSpaceEvolutionResult>
advance_qmom_moment_source_batch_with_index(const std::vector<QmomState>& states, double time, double time_step, const BatchQmomMomentSource& source,
                                            const QmomPhaseSpaceEvolutionOptions& options = {});

[[nodiscard]] std::vector<CqmomPhaseSpaceEvolutionResult>
advance_cqmom_moment_source_batch_with_index(const std::vector<CqmomState>& states, double time, double time_step,
                                             const BatchCqmomMomentSource& source, const CqmomPhaseSpaceEvolutionOptions& options = {});

[[nodiscard]] std::vector<QmomPhaseSpaceEvolutionResult>
advance_qmom_moment_source_batch_with_index(const std::vector<QmomState>& states, double time, double time_step, const BatchQmomMomentSource& source,
                                            const BatchQmomMomentBoundaryContribution& boundary, const QmomPhaseSpaceEvolutionOptions& options);

[[nodiscard]] std::vector<CqmomPhaseSpaceEvolutionResult>
advance_cqmom_moment_source_batch_with_index(const std::vector<CqmomState>& states, double time, double time_step,
                                             const BatchCqmomMomentSource& source, const BatchCqmomMomentBoundaryContribution& boundary,
                                             const CqmomPhaseSpaceEvolutionOptions& options);

// Batch-index-aware complete phase-space evolution batch wrappers.
[[nodiscard]] std::vector<QmomPhaseSpaceEvolutionResult>
advance_qmom_phase_space_batch_with_index(const std::vector<QmomState>& states, double time, double time_step,
                                          const BatchPhaseSpaceVelocity& velocity, const BatchQmomMomentSource& source,
                                          const QmomPhaseSpaceEvolutionOptions& options = {});

[[nodiscard]] std::vector<CqmomPhaseSpaceEvolutionResult>
advance_cqmom_phase_space_batch_with_index(const std::vector<CqmomState>& states, double time, double time_step,
                                           const BatchPhaseSpaceVelocity& velocity, const BatchCqmomMomentSource& source,
                                           const CqmomPhaseSpaceEvolutionOptions& options = {});

[[nodiscard]] std::vector<QmomPhaseSpaceEvolutionResult>
advance_qmom_phase_space_batch_with_index(const std::vector<QmomState>& states, double time, double time_step,
                                          const BatchPhaseSpaceVelocity& velocity, const BatchQmomMomentSource& source,
                                          const BatchQmomMomentBoundaryContribution& boundary, const QmomPhaseSpaceEvolutionOptions& options);

[[nodiscard]] std::vector<CqmomPhaseSpaceEvolutionResult>
advance_cqmom_phase_space_batch_with_index(const std::vector<CqmomState>& states, double time, double time_step,
                                           const BatchPhaseSpaceVelocity& velocity, const BatchCqmomMomentSource& source,
                                           const BatchCqmomMomentBoundaryContribution& boundary, const CqmomPhaseSpaceEvolutionOptions& options);

} // namespace qmomkit
