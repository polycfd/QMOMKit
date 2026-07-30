#include "qmomkit/transport_integrator.hpp"

#include "internal/cqmom_active_counts.hpp"
#include "internal/cqmom_case.hpp"
#include "internal/cqmom_validation.hpp"
#include "internal/moment_reconstruction.hpp"
#include "internal/numeric_utils.hpp"
#include "internal/performance.hpp"
#include "internal/qmom_state_access.hpp"
#include "qmomkit/multi_index.hpp"
#include "qmomkit/quadrature.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace qmomkit
{

// Operator splitting composes each macro step from two numerical operators:
// C: continuous phase-space transport, which advects quadrature nodes and
//    reconstructs raw moments from the advected quadrature.
// D: discontinuous moment evolution, which combines MomentSource with any
//    MomentBoundaryContribution, advances moments with SSP-RK3, then
//    reconstructs quadrature to check realizability before more node motion.
//
// Supported splitting schemes:
//   strang: Strang splitting, C(dt/2) -> D(dt) -> C(dt/2)
//   cd:     continuous-then-discontinuous Lie splitting, C(dt) -> D(dt)
//   dc:     discontinuous-then-continuous Lie splitting, D(dt) -> C(dt)

namespace
{

struct DenseCqmomTransportResult
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

struct LeanCqmomTransportState
{
    CqmomMomentData moments;
};

struct LeanCqmomTransportResult
{
    LeanCqmomTransportState state;
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

void validate_time_input(const double t, const double time_step)
{
    if (!std::isfinite(t) || !std::isfinite(time_step) || !std::isfinite(t + time_step))
    {
        throw std::invalid_argument("Transport time input contains non-finite values.");
    }
    if (time_step <= 0.0)
    {
        throw std::invalid_argument("Transport timestep must be positive.");
    }
}

void validate_transport_scheme(const TransportScheme scheme)
{
    switch (scheme)
    {
    case TransportScheme::strang:
    case TransportScheme::cd:
    case TransportScheme::dc:
        return;
    }
    throw std::invalid_argument("Unsupported transport scheme.");
}

template <typename ApplyContinuous, typename ApplySource>
[[nodiscard]] bool execute_transport_sequence(const TransportScheme scheme, const double t, const double time_step,
                                              ApplyContinuous&& apply_continuous, ApplySource&& apply_source)
{
    const double half_step = 0.5 * time_step;
    switch (scheme)
    {
    // strang: C(dt/2) -> D(dt) -> C(dt/2)
    case TransportScheme::strang:
        if (!apply_continuous(t, half_step, false))
        {
            return false;
        }
        if (!apply_source(t + half_step))
        {
            return false;
        }
        return apply_continuous(t + half_step, half_step, true);
    // cd: C(dt) -> D(dt)
    case TransportScheme::cd:
        if (!apply_continuous(t, time_step, false))
        {
            return false;
        }
        return apply_source(t);
    // dc: D(dt) -> C(dt)
    case TransportScheme::dc:
        if (!apply_source(t))
        {
            return false;
        }
        return apply_continuous(t, time_step, true);
    }
    throw std::invalid_argument("Unsupported transport scheme.");
}

void validate_qmom_input(const QmomState& state, const std::size_t maximum_node_count, const QmomOptions& qmom_options, const double t,
                         const double time_step, const PhaseSpaceVelocity& velocity, const MomentSource& source, const TransportScheme scheme)
{
    validate_time_input(t, time_step);
    validate_transport_scheme(scheme);
    if (maximum_node_count == 0)
    {
        throw std::invalid_argument("QMOM transport requires at least one configured node.");
    }
    if (state.moments().size() < 2 * maximum_node_count)
    {
        throw std::invalid_argument("QMOM transport state has too few moments for the configured maximum node count.");
    }
    if (!internal::all_finite(state.moments()))
    {
        throw std::invalid_argument("QMOM transport state moments must be finite.");
    }
    if (state.active_node_count() > maximum_node_count)
    {
        throw std::invalid_argument("QMOM transport state active node count exceeds the configured maximum.");
    }
    if (!qmom_options.adaptive && state.active_node_count() != maximum_node_count)
    {
        throw std::invalid_argument("Fixed-count QMOM transport state must contain the configured node count.");
    }
    if (!velocity)
    {
        throw std::invalid_argument("Transport phase-space velocity callback must be provided.");
    }
    if (!source)
    {
        throw std::invalid_argument("Transport moment-source callback must be provided.");
    }
}

void validate_boundary(const MomentBoundaryContribution* const boundary)
{
    if (boundary != nullptr && !*boundary)
    {
        throw std::invalid_argument("Moment-boundary contribution callback must be provided.");
    }
}

void validate_dense_cqmom_input(const RawMomentTensor& moments, const PhaseSpaceHierarchy& hierarchy, const double t, const double time_step,
                                const PhaseSpaceVelocity& velocity, const CqmomMomentSource& source, const TransportScheme scheme)
{
    validate_time_input(t, time_step);
    validate_transport_scheme(scheme);
    internal::validate_dense_cqmom_moment_tensor(moments, hierarchy, "CQMOM transport moment tensor rank must match the phase-space hierarchy.",
                                                 "CQMOM transport moment tensor size does not match its extents.",
                                                 "CQMOM transport moments must be finite.");
    if (!velocity)
    {
        throw std::invalid_argument("Transport phase-space velocity callback must be provided.");
    }
    if (!source)
    {
        throw std::invalid_argument("CQMOM source callback must be provided.");
    }
}

void validate_cqmom_boundary(const CqmomMomentBoundaryContribution* const boundary)
{
    if (boundary != nullptr && !*boundary)
    {
        throw std::invalid_argument("CQMOM boundary contribution callback must be provided.");
    }
}

void validate_trusted_lean_cqmom_input(const LeanCqmomTransportState& state, const double t, const double time_step,
                                       const PhaseSpaceVelocity& velocity, const TransportScheme scheme)
{
    validate_time_input(t, time_step);
    validate_transport_scheme(scheme);
    // CqmomState owns already validated layout metadata. The trusted path still
    // checks the evolving values that transport mutates between reconstructions.
    internal::validate_cqmom_moment_values_match_layout(
        state.moments, {"", "Lean CQMOM transport moments must match the fixed layout size.", "Lean CQMOM transport moments must be finite."});
    if (!velocity)
    {
        throw std::invalid_argument("Transport phase-space velocity callback must be provided.");
    }
}

void validate_trusted_lean_cqmom_input(const LeanCqmomTransportState& state, const double t, const double time_step,
                                       const PhaseSpaceVelocity& velocity, const CqmomMomentSource& source, const TransportScheme scheme)
{
    validate_trusted_lean_cqmom_input(state, t, time_step, velocity, scheme);
    if (!source)
    {
        throw std::invalid_argument("Lean CQMOM transport moment-source callback must be provided.");
    }
}

[[nodiscard]] MomentReconstructionDiagnostics successful_reconstruction(const std::size_t moment_count)
{
    return MomentReconstructionDiagnostics{true, moment_count, ""};
}

void record_continuous_diagnostics(TransportDiagnostics& diagnostics, const bool final_continuous, const NodeAdvectionStatus status,
                                   const NodeAdvectionDiagnostics& continuous_diagnostics)
{
    if (final_continuous)
    {
        diagnostics.second_continuous_status = status;
        diagnostics.second_continuous_diagnostics = continuous_diagnostics;
    }
    else
    {
        diagnostics.first_continuous_status = status;
        diagnostics.first_continuous_diagnostics = continuous_diagnostics;
    }
}

template <typename SourceUpdateResult>
void record_source_update_diagnostics(TransportDiagnostics& diagnostics, const SourceUpdateResult& source_update)
{
    diagnostics.source_status = source_update.status;
    diagnostics.source_diagnostics = source_update.diagnostics;
    diagnostics.source_retries = source_update.diagnostics.retries;
}

[[nodiscard]] QmomResult qmom_from_quadrature(const Quadrature1D& quadrature)
{
    QmomStatus status = QmomStatus::success;
    switch (quadrature.status)
    {
    case QuadratureStatus::success:
        status = QmomStatus::success;
        break;
    case QuadratureStatus::invalid_input:
        status = QmomStatus::invalid_input;
        break;
    case QuadratureStatus::not_realizable:
        status = QmomStatus::not_realizable;
        break;
    case QuadratureStatus::no_adaptive_quadrature:
        status = QmomStatus::no_adaptive_quadrature;
        break;
    }
    return QmomResult{quadrature.abscissas, quadrature.weights, status, quadrature.message};
}

void add_continuous_totals(TransportDiagnostics& diagnostics)
{
    diagnostics.total_node_advection_steps =
        diagnostics.first_continuous_diagnostics.step_count + diagnostics.second_continuous_diagnostics.step_count;
    diagnostics.total_node_advection_rejected_steps =
        diagnostics.first_continuous_diagnostics.rejected_step_count + diagnostics.second_continuous_diagnostics.rejected_step_count;
    diagnostics.total_node_advection_rhs_evaluations =
        diagnostics.first_continuous_diagnostics.rhs_evaluations + diagnostics.second_continuous_diagnostics.rhs_evaluations;
}

void finalize_diagnostics(TransportDiagnostics& diagnostics)
{
    add_continuous_totals(diagnostics);
    diagnostics.source_retries = diagnostics.source_diagnostics.retries;
}

void record_failure_diagnostics(TransportDiagnostics& diagnostics, const TransportStage stage, const std::string& message)
{
    diagnostics.failed_stage = stage;
    diagnostics.message = message;
    finalize_diagnostics(diagnostics);
}

[[nodiscard]] QmomTransportResult qmom_failure(QmomState state, const TransportStatus status, const TransportStage stage,
                                               TransportDiagnostics diagnostics, const double final_time, std::string message)
{
    record_failure_diagnostics(diagnostics, stage, message);
    return QmomTransportResult{std::move(state), status, std::move(diagnostics), final_time, std::move(message)};
}

[[nodiscard]] QmomResult select_post_source_qmom(const std::vector<double>& post_source_moments, const std::size_t active_node_count,
                                                 const std::size_t maximum_node_count, const QmomOptions& options)
{
    if (!options.adaptive)
    {
        QmomResult fixed = compute_qmom(MomentSet(post_source_moments), maximum_node_count, options);
        if (fixed.ok() && fixed.weights.size() != maximum_node_count)
        {
            fixed.status = QmomStatus::not_realizable;
            fixed.message = "Fixed-count QMOM reconstruction did not return the configured node count.";
        }
        return fixed;
    }

    QmomResult current = compute_qmom(MomentSet(post_source_moments), active_node_count, options);
    if (!current.ok())
    {
        return current;
    }
    if (current.weights.size() < active_node_count || active_node_count >= maximum_node_count)
    {
        return current;
    }

    QmomResult growth = compute_qmom(MomentSet(post_source_moments), active_node_count + 1, options);
    if (growth.ok() && growth.weights.size() == active_node_count + 1)
    {
        return growth;
    }
    return current;
}

[[nodiscard]] DenseCqmomTransportResult cqmom_failure(RawMomentTensor moments, std::optional<CqmomQuadrature> quadrature,
                                                      const TransportStatus status, const TransportStage stage, TransportDiagnostics diagnostics,
                                                      const double final_time, std::string message)
{
    record_failure_diagnostics(diagnostics, stage, message);
    return DenseCqmomTransportResult{std::move(moments), std::move(quadrature), status, std::move(diagnostics), final_time, std::move(message)};
}

[[nodiscard]] LeanCqmomTransportResult lean_cqmom_failure(LeanCqmomTransportState state, std::optional<CqmomQuadrature> quadrature,
                                                          const TransportStatus status, const TransportStage stage, TransportDiagnostics diagnostics,
                                                          const double final_time, std::string message)
{
    record_failure_diagnostics(diagnostics, stage, message);
    return LeanCqmomTransportResult{std::move(state), std::move(quadrature), status, std::move(diagnostics), final_time, std::move(message)};
}

[[nodiscard]] std::string neutral_cqmom_message(std::string message)
{
    constexpr std::string_view lean_prefix = "Lean CQMOM ";
    constexpr std::string_view neutral_prefix = "CQMOM ";
    std::size_t position = message.find(lean_prefix);
    while (position != std::string::npos)
    {
        message.replace(position, lean_prefix.size(), neutral_prefix);
        position = message.find(lean_prefix, position + neutral_prefix.size());
    }
    return message;
}

[[nodiscard]] TransportDiagnostics neutral_cqmom_diagnostics(TransportDiagnostics diagnostics)
{
    diagnostics.message = neutral_cqmom_message(std::move(diagnostics.message));
    diagnostics.source_diagnostics.message = neutral_cqmom_message(std::move(diagnostics.source_diagnostics.message));
    diagnostics.first_reconstruction_diagnostics.message = neutral_cqmom_message(std::move(diagnostics.first_reconstruction_diagnostics.message));
    diagnostics.final_reconstruction_diagnostics.message = neutral_cqmom_message(std::move(diagnostics.final_reconstruction_diagnostics.message));
    return diagnostics;
}

struct ActiveBranchKey
{
    std::size_t conditioning_level = 0;
    ConditioningPath parent_path;
    std::size_t phase_space_dimension = 0;
};

[[nodiscard]] bool operator<(const ActiveBranchKey& lhs, const ActiveBranchKey& rhs) noexcept
{
    if (lhs.conditioning_level != rhs.conditioning_level)
    {
        return lhs.conditioning_level < rhs.conditioning_level;
    }
    if (lhs.parent_path != rhs.parent_path)
    {
        return lhs.parent_path < rhs.parent_path;
    }
    return lhs.phase_space_dimension < rhs.phase_space_dimension;
}

[[nodiscard]] ActiveBranchKey active_branch_key(const CqmomActiveBranchCount& count)
{
    return {count.conditioning_level, count.parent_path, count.phase_space_dimension};
}

[[nodiscard]] bool path_has_prefix(const ConditioningPath& path, const ConditioningPath& prefix)
{
    if (prefix.size() > path.size())
    {
        return false;
    }
    for (std::size_t level = 0; level < prefix.size(); ++level)
    {
        if (path.node_index_at_level(level) != prefix.node_index_at_level(level))
        {
            return false;
        }
    }
    return true;
}

[[nodiscard]] bool has_changed_strict_ancestor(const CqmomActiveBranchCount& count,
                                               const std::map<ActiveBranchKey, CqmomNodeCountTransitionKind>& changed_keys)
{
    for (const auto& [changed_key, kind] : changed_keys)
    {
        (void)kind;
        if (changed_key.conditioning_level >= count.conditioning_level)
        {
            continue;
        }
        if (path_has_prefix(count.parent_path, changed_key.parent_path))
        {
            return true;
        }
    }
    return false;
}

void record_cqmom_node_count_transitions(TransportDiagnostics& diagnostics, const TransportStage stage,
                                         const std::vector<CqmomActiveBranchCount>& before_counts,
                                         const std::vector<CqmomActiveBranchCount>& after_counts)
{
    std::map<ActiveBranchKey, CqmomActiveBranchCount> before_by_key;
    std::map<ActiveBranchKey, CqmomActiveBranchCount> after_by_key;
    for (const CqmomActiveBranchCount& before : before_counts)
    {
        before_by_key.emplace(active_branch_key(before), before);
    }
    for (const CqmomActiveBranchCount& after : after_counts)
    {
        after_by_key.emplace(active_branch_key(after), after);
    }

    std::map<ActiveBranchKey, CqmomNodeCountTransitionKind> changed_keys;
    for (const auto& [key, before] : before_by_key)
    {
        const auto after = after_by_key.find(key);
        if (after == after_by_key.end())
        {
            changed_keys.emplace(key, CqmomNodeCountTransitionKind::removed);
        }
        else if (before.active_node_count != after->second.active_node_count)
        {
            changed_keys.emplace(key, before.active_node_count > after->second.active_node_count ? CqmomNodeCountTransitionKind::reduced
                                                                                                 : CqmomNodeCountTransitionKind::grew);
        }
    }
    for (const auto& [key, after] : after_by_key)
    {
        if (before_by_key.find(key) == before_by_key.end())
        {
            changed_keys.emplace(key, CqmomNodeCountTransitionKind::created);
        }
    }

    const auto append_transition = [&](const CqmomActiveBranchCount& reference, const std::optional<std::size_t> before,
                                       const std::optional<std::size_t> after, const CqmomNodeCountTransitionKind kind)
    {
        CqmomNodeCountTransition transition;
        transition.stage = stage;
        transition.conditioning_level = reference.conditioning_level;
        transition.parent_path = reference.parent_path;
        transition.phase_space_dimension = reference.phase_space_dimension;
        transition.configured_maximum_node_count = reference.configured_node_count;
        transition.active_node_count_before = before;
        transition.active_node_count_after = after;
        transition.kind = kind;
        transition.rank_reduced = kind == CqmomNodeCountTransitionKind::reduced || kind == CqmomNodeCountTransitionKind::removed;
        transition.rank_grew = kind == CqmomNodeCountTransitionKind::grew || kind == CqmomNodeCountTransitionKind::created;
        diagnostics.cqmom_node_count_transitions.push_back(std::move(transition));
    };

    for (const auto& [key, before] : before_by_key)
    {
        const auto after = after_by_key.find(key);
        const bool ancestor_changed = has_changed_strict_ancestor(before, changed_keys);
        if (after == after_by_key.end() || ancestor_changed)
        {
            append_transition(before, before.active_node_count, std::nullopt, CqmomNodeCountTransitionKind::removed);
        }
        else if (before.active_node_count > after->second.active_node_count)
        {
            append_transition(after->second, before.active_node_count, after->second.active_node_count, CqmomNodeCountTransitionKind::reduced);
        }
        else if (before.active_node_count < after->second.active_node_count)
        {
            append_transition(after->second, before.active_node_count, after->second.active_node_count, CqmomNodeCountTransitionKind::grew);
        }
    }
    for (const auto& [key, after] : after_by_key)
    {
        const auto before = before_by_key.find(key);
        if (before == before_by_key.end() || has_changed_strict_ancestor(after, changed_keys))
        {
            append_transition(after, std::nullopt, after.active_node_count, CqmomNodeCountTransitionKind::created);
        }
    }
}

} // namespace

// Public QMOM overloads dispatch here; this helper applies the selected
// C/D numerical transport sequence to an accepted one-dimensional state.
static QmomTransportResult advance_transport_qmom_impl(const QmomState& state, const std::size_t maximum_node_count, const double t,
                                                       const double time_step, const PhaseSpaceVelocity& velocity, const MomentSource& source,
                                                       const MomentBoundaryContribution* const boundary, const QmomTransportOptions& options)
{
    internal::performance::record_transport_macro_step();
    TransportDiagnostics diagnostics;
    diagnostics.scheme = options.scheme;
    try
    {
        validate_qmom_input(state, maximum_node_count, options.qmom_options, t, time_step, velocity, source, options.scheme);
        validate_boundary(boundary);
    }
    catch (const std::exception& error)
    {
        return qmom_failure(state, TransportStatus::invalid_input, TransportStage::none, diagnostics, t, error.what());
    }

    QmomState current_state = state;
    std::optional<QmomTransportResult> failure;

    const auto apply_continuous = [&](const double start_time, const double step, const bool final_continuous)
    {
        internal::performance::record_continuous_operator_call();
        const QmomNodeAdvectionResult continuous = advance_nodes(current_state.quadrature(), start_time, step, velocity, options.node_options);
        record_continuous_diagnostics(diagnostics, final_continuous, continuous.status, continuous.diagnostics);
        if (!continuous.ok())
        {
            failure = qmom_failure(state, TransportStatus::continuous_failed,
                                   final_continuous ? TransportStage::second_continuous_half_step : TransportStage::first_continuous_half_step,
                                   diagnostics, start_time, continuous.message);
            return false;
        }

        try
        {
            internal::performance::record_moment_reconstruction_call();
            current_state = internal::QmomStateAccess::from_quadrature(qmom_from_quadrature(continuous.quadrature), state.moments().size());
            if (final_continuous)
            {
                diagnostics.final_reconstruction_diagnostics = successful_reconstruction(current_state.moments().size());
            }
            else
            {
                diagnostics.first_reconstruction_diagnostics = successful_reconstruction(current_state.moments().size());
            }
        }
        catch (const std::exception& error)
        {
            if (final_continuous)
            {
                diagnostics.final_reconstruction_diagnostics = {false, state.moments().size(), error.what()};
            }
            else
            {
                diagnostics.first_reconstruction_diagnostics = {false, state.moments().size(), error.what()};
            }
            failure = qmom_failure(state, TransportStatus::moment_reconstruction_failed,
                                   final_continuous ? TransportStage::final_moment_reconstruction : TransportStage::first_moment_reconstruction,
                                   diagnostics, start_time + step, error.what());
            return false;
        }
        return true;
    };

    const auto apply_source = [&](const double start_time)
    {
        internal::performance::record_discontinuous_operator_call();
        // D-step: MomentSource and optional boundary contribution share the
        // existing SSP-RK3 source update path.
        const MomentSourceUpdateResult source_update =
            boundary == nullptr
                ? advance_moments_with_source(current_state.moments(), start_time, time_step, source, options.source_options)
                : advance_moments_with_source(current_state.moments(), start_time, time_step, source, *boundary, options.source_options);
        record_source_update_diagnostics(diagnostics, source_update);
        if (!source_update.ok())
        {
            failure = qmom_failure(state, TransportStatus::source_failed, TransportStage::discontinuous_source_step, diagnostics, start_time,
                                   source_update.message);
            return false;
        }

        internal::performance::record_post_source_inversion();
        // Reconstruct quadrature after the moment update so realizability is
        // checked before any later continuous stage can use these moments.
        QmomResult selected_quadrature =
            select_post_source_qmom(source_update.moments, current_state.active_node_count(), maximum_node_count, options.qmom_options);
        diagnostics.post_source_qmom_status = selected_quadrature.status;
        if (!selected_quadrature.ok())
        {
            failure = qmom_failure(state, TransportStatus::inversion_failed, TransportStage::post_source_inversion, diagnostics, start_time,
                                   selected_quadrature.message);
            return false;
        }
        try
        {
            internal::performance::record_moment_reconstruction_call();
            current_state = internal::QmomStateAccess::from_quadrature(std::move(selected_quadrature), state.moments().size());
        }
        catch (const std::exception& error)
        {
            std::string message = std::string("QMOM post-source moment reconstruction failed: ") + error.what();
            diagnostics.final_reconstruction_diagnostics = {false, state.moments().size(), message};
            failure = qmom_failure(state, TransportStatus::moment_reconstruction_failed, TransportStage::post_source_inversion, diagnostics,
                                   start_time, std::move(message));
            return false;
        }
        return true;
    };

    if (!execute_transport_sequence(options.scheme, t, time_step, apply_continuous, apply_source))
    {
        return std::move(*failure);
    }

    finalize_diagnostics(diagnostics);
    return QmomTransportResult{std::move(current_state), TransportStatus::success, std::move(diagnostics), t + time_step, ""};
}

QmomTransportResult advance_transport_qmom(const QmomState& state, const std::size_t maximum_node_count, const double t, const double time_step,
                                           const PhaseSpaceVelocity& velocity, const MomentSource& source, const QmomTransportOptions& options)
{
    return advance_transport_qmom_impl(state, maximum_node_count, t, time_step, velocity, source, nullptr, options);
}

QmomTransportResult advance_transport_qmom(const QmomState& state, const std::size_t maximum_node_count, const double t, const double time_step,
                                           const PhaseSpaceVelocity& velocity, const MomentSource& source, const MomentBoundaryContribution& boundary,
                                           const QmomTransportOptions& options)
{
    return advance_transport_qmom_impl(state, maximum_node_count, t, time_step, velocity, source, &boundary, options);
}

// Dense-backed CqmomState transport applies the selected C/D sequence while the
// state owns full raw-tensor storage internally.
static DenseCqmomTransportResult advance_transport_cqmom_dense_impl(const CqmomDenseCase& dense_case, const double t, const double time_step,
                                                                    const PhaseSpaceVelocity& velocity, const CqmomMomentSource& source,
                                                                    const CqmomMomentBoundaryContribution* const boundary,
                                                                    const CqmomTransportOptions& options,
                                                                    std::vector<CqmomActiveBranchCount> initial_active_branch_counts)
{
    internal::performance::record_transport_macro_step();
    TransportDiagnostics diagnostics;
    diagnostics.scheme = options.scheme;
    diagnostics.cqmom_moment_storage = CqmomMomentStorage::dense;
    diagnostics.cqmom_diagnostic_level = dense_case.diagnostic_level;
    diagnostics.cqmom_vandermonde_solver = dense_case.vandermonde_solver;
    try
    {
        validate_dense_cqmom_input(dense_case.moments, dense_case.hierarchy, t, time_step, velocity, source, options.scheme);
        validate_cqmom_boundary(boundary);
    }
    catch (const std::exception& error)
    {
        return cqmom_failure(dense_case.moments, std::nullopt, TransportStatus::invalid_input, TransportStage::none, diagnostics, t, error.what());
    }

    std::vector<CqmomActiveBranchCount> current_active_branch_counts = std::move(initial_active_branch_counts);
    const std::vector<CqmomActiveBranchCount> before_initial_active_branch_counts = current_active_branch_counts;
    CqmomResult initial_quadrature = internal::compute_trusted_cqmom(dense_case, current_active_branch_counts);
    diagnostics.initial_cqmom_status = initial_quadrature.status;
    diagnostics.initial_cqmom_diagnostics = initial_quadrature.diagnostics;
    if (!initial_quadrature.ok())
    {
        return cqmom_failure(dense_case.moments, initial_quadrature.quadrature, TransportStatus::inversion_failed, TransportStage::initial_inversion,
                             diagnostics, t, initial_quadrature.message);
    }

    RawMomentTensor current_moments = dense_case.moments;
    CqmomQuadrature current_quadrature = initial_quadrature.quadrature;
    current_active_branch_counts = internal::active_branch_counts_from_quadrature(initial_quadrature.quadrature);
    record_cqmom_node_count_transitions(diagnostics, TransportStage::initial_inversion, before_initial_active_branch_counts,
                                        current_active_branch_counts);
    const CqmomMomentLayout source_layout = create_cqmom_moment_layout(dense_case.hierarchy);
    std::optional<DenseCqmomTransportResult> failure;

    const auto apply_continuous = [&](const double start_time, const double step, const bool final_continuous)
    {
        internal::performance::record_continuous_operator_call();
        const CqmomNodeAdvectionResult continuous = advance_nodes(current_quadrature, start_time, step, velocity, options.node_options);
        record_continuous_diagnostics(diagnostics, final_continuous, continuous.status, continuous.diagnostics);
        if (!continuous.ok())
        {
            failure = cqmom_failure(std::move(current_moments), current_quadrature, TransportStatus::continuous_failed,
                                    final_continuous ? TransportStage::second_continuous_half_step : TransportStage::first_continuous_half_step,
                                    diagnostics, start_time, continuous.message);
            return false;
        }

        try
        {
            internal::performance::record_moment_reconstruction_call();
            current_moments = reconstruct_raw_moments(continuous.quadrature, dense_case.moments.extents);
            if (final_continuous)
            {
                diagnostics.final_reconstruction_diagnostics = successful_reconstruction(current_moments.values.size());
            }
            else
            {
                diagnostics.first_reconstruction_diagnostics = successful_reconstruction(current_moments.values.size());
            }
            current_quadrature = continuous.quadrature;
            current_active_branch_counts = internal::active_branch_counts_from_quadrature(current_quadrature);
        }
        catch (const std::exception& error)
        {
            if (final_continuous)
            {
                diagnostics.final_reconstruction_diagnostics = {false, dense_case.moments.values.size(), error.what()};
            }
            else
            {
                diagnostics.first_reconstruction_diagnostics = {false, dense_case.moments.values.size(), error.what()};
            }
            failure = cqmom_failure(std::move(current_moments), continuous.quadrature, TransportStatus::moment_reconstruction_failed,
                                    final_continuous ? TransportStage::final_moment_reconstruction : TransportStage::first_moment_reconstruction,
                                    diagnostics, start_time + step, error.what());
            return false;
        }
        return true;
    };

    const auto apply_source = [&](const double start_time)
    {
        internal::performance::record_discontinuous_operator_call();
        // D-step: CQMOM-native source and boundary callbacks update only the
        // canonical transport layout. Dense-only tensor entries are closure
        // values reconstructed from the accepted post-source quadrature below.
        CqmomMomentData source_moments = extract_cqmom_moments(current_moments, source_layout);
        const CqmomMomentSourceUpdateResult source_update =
            boundary == nullptr ? advance_cqmom_moments_with_source(source_moments, start_time, time_step, source, options.source_options)
                                : advance_cqmom_moments_with_source(source_moments, start_time, time_step, source, *boundary, options.source_options);
        record_source_update_diagnostics(diagnostics, source_update);
        if (!source_update.ok())
        {
            failure = cqmom_failure(std::move(current_moments), current_quadrature, TransportStatus::source_failed,
                                    TransportStage::discontinuous_source_step, diagnostics, start_time, source_update.message);
            return false;
        }

        const std::vector<CqmomActiveBranchCount> before_source_active_branch_counts = current_active_branch_counts;
        const bool allow_active_branch_count_growth = dense_case.adaptive_branch_counts;
        // Reconstruct quadrature after the moment update so CQMOM
        // realizability is checked before any later continuous stage.
        internal::performance::record_post_source_inversion();
        CqmomLeanCase post_source_case{source_update.moments,       dense_case.conditional_options,
                                       dense_case.qmom_options,     dense_case.adaptive_branch_counts,
                                       dense_case.diagnostic_level, dense_case.vandermonde_solver};
        CqmomResult post_source_quadrature =
            internal::compute_trusted_cqmom(post_source_case, current_active_branch_counts, allow_active_branch_count_growth);
        diagnostics.post_source_cqmom_status = post_source_quadrature.status;
        diagnostics.post_source_cqmom_diagnostics = post_source_quadrature.diagnostics;
        if (!post_source_quadrature.ok())
        {
            failure = cqmom_failure(std::move(current_moments), post_source_quadrature.quadrature, TransportStatus::inversion_failed,
                                    TransportStage::post_source_inversion, diagnostics, start_time,
                                    neutral_cqmom_message(std::move(post_source_quadrature.message)));
            return false;
        }
        current_quadrature = post_source_quadrature.quadrature;
        try
        {
            internal::performance::record_moment_reconstruction_call();
            current_moments = reconstruct_raw_moments(current_quadrature, dense_case.moments.extents);
        }
        catch (const std::exception& error)
        {
            std::string message = std::string("Dense CQMOM post-source moment reconstruction failed: ") + error.what();
            diagnostics.final_reconstruction_diagnostics = {false, dense_case.moments.values.size(), message};
            failure = cqmom_failure(std::move(current_moments), current_quadrature, TransportStatus::moment_reconstruction_failed,
                                    TransportStage::post_source_inversion, diagnostics, start_time, std::move(message));
            return false;
        }
        current_active_branch_counts = internal::active_branch_counts_from_quadrature(current_quadrature);
        record_cqmom_node_count_transitions(diagnostics, TransportStage::post_source_inversion, before_source_active_branch_counts,
                                            current_active_branch_counts);
        return true;
    };

    if (!execute_transport_sequence(options.scheme, t, time_step, apply_continuous, apply_source))
    {
        return std::move(*failure);
    }

    finalize_diagnostics(diagnostics);
    return DenseCqmomTransportResult{
        std::move(current_moments), std::move(current_quadrature), TransportStatus::success, std::move(diagnostics), t + time_step, ""};
}

namespace
{

// Lean-backed CqmomState transport applies the selected C/D sequence to the
// fixed CQMOM moment layout and rebuilds quadrature at stage boundaries.
LeanCqmomTransportResult advance_transport_cqmom_lean_impl(const CqmomLeanCase& lean_case, const double t, const double time_step,
                                                           const PhaseSpaceVelocity& velocity, const CqmomMomentSource& source,
                                                           const CqmomMomentBoundaryContribution* const boundary,
                                                           const CqmomTransportOptions& options,
                                                           std::vector<CqmomActiveBranchCount> initial_active_branch_counts = {})
{
    internal::performance::record_transport_macro_step();
    TransportDiagnostics diagnostics;
    diagnostics.scheme = options.scheme;
    diagnostics.cqmom_moment_storage = CqmomMomentStorage::lean;
    diagnostics.cqmom_diagnostic_level = lean_case.diagnostic_level;
    diagnostics.cqmom_vandermonde_solver = lean_case.vandermonde_solver;
    const LeanCqmomTransportState state{lean_case.moments};
    try
    {
        validate_trusted_lean_cqmom_input(state, t, time_step, velocity, source, options.scheme);
        validate_cqmom_boundary(boundary);
    }
    catch (const std::exception& error)
    {
        return lean_cqmom_failure(state, std::nullopt, TransportStatus::invalid_input, TransportStage::none, diagnostics, t, error.what());
    }

    const auto compute_stage_cqmom = [&](const CqmomMomentData& moments, const std::vector<CqmomActiveBranchCount>& active_branch_counts,
                                         const bool allow_active_branch_count_growth = false)
    {
        const CqmomLeanCase stage_case{moments,
                                       lean_case.conditional_options,
                                       lean_case.qmom_options,
                                       lean_case.adaptive_branch_counts,
                                       lean_case.diagnostic_level,
                                       lean_case.vandermonde_solver};
        return internal::compute_trusted_cqmom(stage_case, active_branch_counts, allow_active_branch_count_growth);
    };

    std::vector<CqmomActiveBranchCount> current_active_branch_counts = std::move(initial_active_branch_counts);
    const std::vector<CqmomActiveBranchCount> before_initial_active_branch_counts = current_active_branch_counts;
    CqmomResult initial_quadrature = compute_stage_cqmom(state.moments, current_active_branch_counts);
    diagnostics.initial_cqmom_status = initial_quadrature.status;
    diagnostics.initial_cqmom_diagnostics = initial_quadrature.diagnostics;
    if (!initial_quadrature.ok())
    {
        return lean_cqmom_failure(state, initial_quadrature.quadrature, TransportStatus::inversion_failed, TransportStage::initial_inversion,
                                  diagnostics, t, initial_quadrature.message);
    }
    current_active_branch_counts = internal::active_branch_counts_from_quadrature(initial_quadrature.quadrature);
    record_cqmom_node_count_transitions(diagnostics, TransportStage::initial_inversion, before_initial_active_branch_counts,
                                        current_active_branch_counts);

    LeanCqmomTransportState current_state = state;
    CqmomQuadrature current_quadrature = initial_quadrature.quadrature;
    std::optional<LeanCqmomTransportResult> failure;

    const auto apply_continuous = [&](const double start_time, const double step, const bool final_continuous)
    {
        internal::performance::record_continuous_operator_call();
        const CqmomNodeAdvectionResult continuous = advance_nodes(current_quadrature, start_time, step, velocity, options.node_options);
        record_continuous_diagnostics(diagnostics, final_continuous, continuous.status, continuous.diagnostics);
        if (!continuous.ok())
        {
            failure = lean_cqmom_failure(std::move(current_state), current_quadrature, TransportStatus::continuous_failed,
                                         final_continuous ? TransportStage::second_continuous_half_step : TransportStage::first_continuous_half_step,
                                         diagnostics, start_time, continuous.message);
            return false;
        }

        try
        {
            internal::performance::record_moment_reconstruction_call();
            // The trusted lean-backed path preserves the existing fixed layout
            // and reconstructs only updated moment values after node advection.
            current_state.moments = internal::reconstruct_trusted_cqmom_moments(continuous.quadrature, current_state.moments.layout);
            if (final_continuous)
            {
                diagnostics.final_reconstruction_diagnostics = successful_reconstruction(current_state.moments.values.size());
            }
            else
            {
                diagnostics.first_reconstruction_diagnostics = successful_reconstruction(current_state.moments.values.size());
            }
            current_quadrature = continuous.quadrature;
            current_active_branch_counts = internal::active_branch_counts_from_quadrature(current_quadrature);
        }
        catch (const std::exception& error)
        {
            if (final_continuous)
            {
                diagnostics.final_reconstruction_diagnostics = {false, current_state.moments.values.size(), error.what()};
            }
            else
            {
                diagnostics.first_reconstruction_diagnostics = {false, current_state.moments.values.size(), error.what()};
            }
            failure = lean_cqmom_failure(std::move(current_state), continuous.quadrature, TransportStatus::moment_reconstruction_failed,
                                         final_continuous ? TransportStage::final_moment_reconstruction : TransportStage::first_moment_reconstruction,
                                         diagnostics, start_time + step, error.what());
            return false;
        }
        return true;
    };

    const auto apply_source = [&](const double start_time)
    {
        internal::performance::record_discontinuous_operator_call();
        // D-step: CqmomMomentSource and optional boundary contribution share
        // the existing SSP-RK3 source update path.
        const CqmomMomentSourceUpdateResult source_update =
            boundary == nullptr
                ? advance_cqmom_moments_with_source(current_state.moments, start_time, time_step, source, options.source_options)
                : advance_cqmom_moments_with_source(current_state.moments, start_time, time_step, source, *boundary, options.source_options);
        record_source_update_diagnostics(diagnostics, source_update);
        if (!source_update.ok())
        {
            failure = lean_cqmom_failure(std::move(current_state), current_quadrature, TransportStatus::source_failed,
                                         TransportStage::discontinuous_source_step, diagnostics, start_time, source_update.message);
            return false;
        }

        current_state.moments = source_update.moments;
        const std::vector<CqmomActiveBranchCount> before_source_active_branch_counts = current_active_branch_counts;
        const bool allow_active_branch_count_growth = lean_case.adaptive_branch_counts;
        // Reconstruct quadrature after the moment update so CQMOM
        // realizability is checked before any later continuous stage.
        internal::performance::record_post_source_inversion();
        CqmomResult post_source_quadrature =
            compute_stage_cqmom(current_state.moments, current_active_branch_counts, allow_active_branch_count_growth);
        diagnostics.post_source_cqmom_status = post_source_quadrature.status;
        diagnostics.post_source_cqmom_diagnostics = post_source_quadrature.diagnostics;
        if (!post_source_quadrature.ok())
        {
            failure = lean_cqmom_failure(std::move(current_state), post_source_quadrature.quadrature, TransportStatus::inversion_failed,
                                         TransportStage::post_source_inversion, diagnostics, start_time, post_source_quadrature.message);
            return false;
        }
        current_quadrature = post_source_quadrature.quadrature;
        current_active_branch_counts = internal::active_branch_counts_from_quadrature(current_quadrature);
        record_cqmom_node_count_transitions(diagnostics, TransportStage::post_source_inversion, before_source_active_branch_counts,
                                            current_active_branch_counts);
        return true;
    };

    if (!execute_transport_sequence(options.scheme, t, time_step, apply_continuous, apply_source))
    {
        return std::move(*failure);
    }

    finalize_diagnostics(diagnostics);
    return LeanCqmomTransportResult{
        std::move(current_state), std::move(current_quadrature), TransportStatus::success, std::move(diagnostics), t + time_step, ""};
}

} // namespace

static CqmomTransportResult advance_transport_cqmom_state_impl(const CqmomState& state, double t, double time_step,
                                                               const PhaseSpaceVelocity& velocity, const CqmomMomentSource& source,
                                                               const CqmomMomentBoundaryContribution* boundary, const CqmomTransportOptions& options);

CqmomTransportResult advance_transport_cqmom(const CqmomState& state, const double t, const double time_step, const PhaseSpaceVelocity& velocity,
                                             const CqmomMomentSource& source, const CqmomTransportOptions& options)
{
    return advance_transport_cqmom_state_impl(state, t, time_step, velocity, source, nullptr, options);
}

CqmomTransportResult advance_transport_cqmom(const CqmomState& state, const double t, const double time_step, const PhaseSpaceVelocity& velocity,
                                             const CqmomMomentSource& source, const CqmomMomentBoundaryContribution& boundary,
                                             const CqmomTransportOptions& options)
{
    return advance_transport_cqmom_state_impl(state, t, time_step, velocity, source, &boundary, options);
}

static CqmomTransportResult advance_transport_cqmom_state_impl(const CqmomState& state, const double t, const double time_step,
                                                               const PhaseSpaceVelocity& velocity, const CqmomMomentSource& source,
                                                               const CqmomMomentBoundaryContribution* const boundary,
                                                               const CqmomTransportOptions& options)
{
    const std::vector<CqmomActiveBranchCount>& active_branch_counts = internal::cqmom_state_active_branch_counts(state);

    if (const CqmomDenseCase* dense_case = internal::cqmom_state_dense_case(state))
    {
        DenseCqmomTransportResult dense_result =
            advance_transport_cqmom_dense_impl(*dense_case, t, time_step, velocity, source, boundary, options, active_branch_counts);
        if (!dense_result.ok())
        {
            return CqmomTransportResult{state,
                                        std::move(dense_result.quadrature),
                                        dense_result.status,
                                        std::move(dense_result.diagnostics),
                                        dense_result.final_time,
                                        std::move(dense_result.message)};
        }

        CqmomDenseCase updated_case{dense_case->hierarchy,         std::move(dense_result.moments),    dense_case->conditional_options,
                                    dense_case->qmom_options,      dense_case->adaptive_branch_counts, dense_case->diagnostic_level,
                                    dense_case->vandermonde_solver};
        std::vector<CqmomActiveBranchCount> updated_active_branch_counts;
        if (dense_result.quadrature.has_value())
        {
            updated_active_branch_counts = internal::active_branch_counts_from_quadrature(*dense_result.quadrature);
        }
        return CqmomTransportResult{internal::create_cqmom_state(std::move(updated_case), std::move(updated_active_branch_counts)),
                                    std::move(dense_result.quadrature),
                                    dense_result.status,
                                    std::move(dense_result.diagnostics),
                                    dense_result.final_time,
                                    std::move(dense_result.message)};
    }

    const CqmomLeanCase* lean_case = internal::cqmom_state_lean_case(state);
    if (lean_case == nullptr)
    {
        TransportDiagnostics diagnostics;
        diagnostics.scheme = options.scheme;
        return CqmomTransportResult{state, std::nullopt, TransportStatus::invalid_input, std::move(diagnostics), t, "Unknown CQMOM state storage."};
    }

    LeanCqmomTransportResult lean_result =
        advance_transport_cqmom_lean_impl(*lean_case, t, time_step, velocity, source, boundary, options, active_branch_counts);

    if (!lean_result.ok())
    {
        TransportDiagnostics diagnostics = neutral_cqmom_diagnostics(std::move(lean_result.diagnostics));
        std::string message = neutral_cqmom_message(std::move(lean_result.message));
        return CqmomTransportResult{
            state, std::move(lean_result.quadrature), lean_result.status, std::move(diagnostics), lean_result.final_time, std::move(message)};
    }

    CqmomLeanCase updated_case{lean_result.state.moments,         lean_case->conditional_options, lean_case->qmom_options,
                               lean_case->adaptive_branch_counts, lean_case->diagnostic_level,    lean_case->vandermonde_solver};
    std::vector<CqmomActiveBranchCount> updated_active_branch_counts;
    if (lean_result.quadrature.has_value())
    {
        updated_active_branch_counts = internal::active_branch_counts_from_quadrature(*lean_result.quadrature);
    }
    TransportDiagnostics diagnostics = neutral_cqmom_diagnostics(std::move(lean_result.diagnostics));
    std::string message = neutral_cqmom_message(std::move(lean_result.message));
    return CqmomTransportResult{internal::create_cqmom_state(std::move(updated_case), std::move(updated_active_branch_counts)),
                                std::move(lean_result.quadrature),
                                lean_result.status,
                                std::move(diagnostics),
                                lean_result.final_time,
                                std::move(message)};
}

} // namespace qmomkit
