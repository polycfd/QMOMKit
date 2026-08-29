#include "qmomkit/phase_space_evolution.hpp"

#include "internal/bounded_quadrature.hpp"
#include "internal/cqmom_active_counts.hpp"
#include "internal/cqmom_case.hpp"
#include "internal/cqmom_validation.hpp"
#include "internal/moment_reconstruction.hpp"
#include "internal/moment_source.hpp"
#include "internal/numeric_utils.hpp"
#include "internal/performance.hpp"
#include "internal/qmom_state_access.hpp"
#include "internal/recursive_quadrature.hpp"
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
#include <type_traits>
#include <utility>
#include <vector>

namespace qmomkit
{

// Operator splitting composes each macro step from two numerical operations:
// phase-space advection, which advects quadrature nodes and reconstructs raw
// moments from the advected quadrature; and moment-source evolution, which
// combines the caller-supplied moment source with any boundary contribution,
// advances moments with SSP-RK3, then reconstructs quadrature to check
// realizability before more node motion.
//
// Supported splitting schemes:
//   strang: Strang splitting, advection(dt/2) -> moment-source(dt) -> advection(dt/2)
//   advsrc: advection-then-moment-source Lie splitting, advection(dt) -> moment-source(dt)
//   srcadv: moment-source-then-advection Lie splitting, moment-source(dt) -> advection(dt)

namespace
{

template <typename State, typename Apply>
[[nodiscard]] auto apply_batch(const std::vector<State>& states, Apply&& apply)
    -> std::vector<std::invoke_result_t<Apply&, std::size_t, const State&>>
{
    using Result = std::invoke_result_t<Apply&, std::size_t, const State&>;
    std::vector<Result> results;
    results.reserve(states.size());
    for (std::size_t batch_index = 0; batch_index < states.size(); ++batch_index)
    {
        results.push_back(apply(batch_index, states[batch_index]));
    }
    return results;
}

void validate_time_input(const double t, const double time_step)
{
    if (!std::isfinite(t) || !std::isfinite(time_step) || !std::isfinite(t + time_step))
    {
        throw std::invalid_argument("Phase-space evolution time input contains non-finite values.");
    }
    if (time_step <= 0.0)
    {
        throw std::invalid_argument("Phase-space evolution timestep must be positive.");
    }
}

void validate_phase_space_evolution_scheme(const PhaseSpaceEvolutionScheme scheme)
{
    switch (scheme)
    {
    case PhaseSpaceEvolutionScheme::strang:
    case PhaseSpaceEvolutionScheme::advsrc:
    case PhaseSpaceEvolutionScheme::srcadv:
        return;
    }
    throw std::invalid_argument("Unsupported phase-space evolution scheme.");
}

template <typename ApplyAdvection, typename ApplySource>
[[nodiscard]] bool execute_phase_space_evolution_sequence(const PhaseSpaceEvolutionScheme scheme, const double t, const double time_step,
                                                          ApplyAdvection&& apply_advection, ApplySource&& apply_source)
{
    const double half_step = 0.5 * time_step;
    switch (scheme)
    {
    // strang: C(dt/2) -> D(dt) -> C(dt/2)
    case PhaseSpaceEvolutionScheme::strang:
        if (!apply_advection(t, half_step, false))
        {
            return false;
        }
        if (!apply_source(t + half_step))
        {
            return false;
        }
        return apply_advection(t + half_step, half_step, true);
    // advsrc: C(dt) -> D(dt)
    case PhaseSpaceEvolutionScheme::advsrc:
        if (!apply_advection(t, time_step, false))
        {
            return false;
        }
        return apply_source(t);
    // srcadv: D(dt) -> C(dt)
    case PhaseSpaceEvolutionScheme::srcadv:
        if (!apply_source(t))
        {
            return false;
        }
        return apply_advection(t, time_step, true);
    }
    throw std::invalid_argument("Unsupported phase-space evolution scheme.");
}

void validate_qmom_input(const QmomState& state, const std::size_t maximum_node_count, const QmomOptions& qmom_options, const double t,
                         const double time_step, const PhaseSpaceVelocity& velocity, const MomentSource& source,
                         const PhaseSpaceEvolutionScheme scheme)
{
    validate_time_input(t, time_step);
    validate_phase_space_evolution_scheme(scheme);
    if (maximum_node_count == 0)
    {
        throw std::invalid_argument("QMOM phase-space evolution requires at least one configured node.");
    }
    if (state.moments().size() < 2 * maximum_node_count)
    {
        throw std::invalid_argument("QMOM phase-space evolution state has too few moments for the configured maximum node count.");
    }
    if (!internal::all_finite(state.moments()))
    {
        throw std::invalid_argument("QMOM phase-space evolution state moments must be finite.");
    }
    if (state.active_node_count() > maximum_node_count)
    {
        throw std::invalid_argument("QMOM phase-space evolution state active node count exceeds the configured maximum.");
    }
    if (!qmom_options.adaptive && state.active_node_count() != maximum_node_count)
    {
        throw std::invalid_argument("Fixed-count QMOM phase-space evolution state must contain the configured node count.");
    }
    if (!velocity)
    {
        throw std::invalid_argument("A phase-space velocity callback must be provided for phase-space evolution.");
    }
    if (!source)
    {
        throw std::invalid_argument("A moment-source callback must be provided for phase-space evolution.");
    }
}

void validate_boundary(const MomentBoundaryContribution* const boundary)
{
    if (boundary != nullptr && !*boundary)
    {
        throw std::invalid_argument("Moment-boundary contribution callback must be provided.");
    }
}

void validate_cqmom_boundary(const CqmomMomentBoundaryContribution* const boundary)
{
    if (boundary != nullptr && !*boundary)
    {
        throw std::invalid_argument("CQMOM boundary contribution callback must be provided.");
    }
}

[[nodiscard]] MomentReconstructionDiagnostics successful_reconstruction(const std::size_t moment_count)
{
    return MomentReconstructionDiagnostics{true, moment_count, ""};
}

void record_advection_diagnostics(PhaseSpaceEvolutionDiagnostics& diagnostics, const bool final_advection, const NodeAdvectionStatus status,
                                  const NodeAdvectionDiagnostics& advection_diagnostics)
{
    if (final_advection)
    {
        diagnostics.second_advection_status = status;
        diagnostics.second_advection_diagnostics = advection_diagnostics;
    }
    else
    {
        diagnostics.first_advection_status = status;
        diagnostics.first_advection_diagnostics = advection_diagnostics;
    }
}

template <typename SourceUpdateResult>
void record_source_update_diagnostics(PhaseSpaceEvolutionDiagnostics& diagnostics, const SourceUpdateResult& source_update)
{
    diagnostics.moment_source_status = source_update.status;
    diagnostics.moment_source_diagnostics = source_update.diagnostics;
    diagnostics.moment_source_retries = source_update.diagnostics.retries;
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
    }
    return QmomResult{quadrature.abscissas, quadrature.weights, status, quadrature.message};
}

void add_advection_totals(PhaseSpaceEvolutionDiagnostics& diagnostics)
{
    diagnostics.total_node_advection_steps = diagnostics.first_advection_diagnostics.step_count + diagnostics.second_advection_diagnostics.step_count;
    diagnostics.total_node_advection_rejected_steps =
        diagnostics.first_advection_diagnostics.rejected_step_count + diagnostics.second_advection_diagnostics.rejected_step_count;
    diagnostics.total_node_advection_rhs_evaluations =
        diagnostics.first_advection_diagnostics.rhs_evaluations + diagnostics.second_advection_diagnostics.rhs_evaluations;
}

void finalize_diagnostics(PhaseSpaceEvolutionDiagnostics& diagnostics)
{
    add_advection_totals(diagnostics);
    diagnostics.moment_source_retries = diagnostics.moment_source_diagnostics.retries;
}

void record_failure_diagnostics(PhaseSpaceEvolutionDiagnostics& diagnostics, const PhaseSpaceEvolutionStage stage, const std::string& message)
{
    diagnostics.failed_stage = stage;
    diagnostics.message = message;
    finalize_diagnostics(diagnostics);
}

[[nodiscard]] QmomPhaseSpaceEvolutionResult qmom_failure(QmomState state, const PhaseSpaceEvolutionStatus status,
                                                         const PhaseSpaceEvolutionStage stage, PhaseSpaceEvolutionDiagnostics diagnostics,
                                                         const double final_time, std::string message)
{
    record_failure_diagnostics(diagnostics, stage, message);
    return QmomPhaseSpaceEvolutionResult{std::move(state), status, std::move(diagnostics), final_time, std::move(message)};
}

[[nodiscard]] QmomResult select_post_source_qmom(const std::vector<double>& post_source_moments, const std::size_t active_node_count,
                                                 const std::size_t maximum_node_count, const QmomOptions& options,
                                                 const std::optional<AbscissaBounds>& abscissa_bounds)
{
    if (abscissa_bounds.has_value())
    {
        const internal::BoundedQuadratureSelectionRequest request{options.adaptive ? active_node_count : maximum_node_count, maximum_node_count,
                                                                  options.adaptive ? std::optional<std::size_t>{active_node_count} : std::nullopt,
                                                                  options.adaptive, options.adaptive};
        const internal::BoundedQuadratureSelection selection = internal::select_bounded_quadrature(
            post_source_moments, request, *abscissa_bounds, WheelerOptions{options.adaptive, options.rmin, options.eabs, options.cutoff});
        return qmom_from_quadrature(selection.quadrature);
    }

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

void validate_qmom_advection_input(const QmomState& state, const double t, const double time_step, const PhaseSpaceVelocity& velocity)
{
    validate_time_input(t, time_step);
    if (!internal::all_finite(state.moments()))
    {
        throw std::invalid_argument("QMOM phase-space evolution state moments must be finite.");
    }
    if (!velocity)
    {
        throw std::invalid_argument("A phase-space velocity callback must be provided for phase-space evolution.");
    }
}

void validate_qmom_moment_source_input(const QmomState& state, const std::size_t maximum_node_count, const QmomOptions& qmom_options, const double t,
                                       const double time_step, const MomentSource& source, const MomentBoundaryContribution* const boundary)
{
    validate_time_input(t, time_step);
    if (maximum_node_count == 0)
    {
        throw std::invalid_argument("QMOM phase-space evolution requires at least one configured node.");
    }
    if (state.moments().size() < 2 * maximum_node_count)
    {
        throw std::invalid_argument("QMOM phase-space evolution state has too few moments for the configured maximum node count.");
    }
    if (!internal::all_finite(state.moments()))
    {
        throw std::invalid_argument("QMOM phase-space evolution state moments must be finite.");
    }
    if (state.active_node_count() > maximum_node_count)
    {
        throw std::invalid_argument("QMOM phase-space evolution state active node count exceeds the configured maximum.");
    }
    if (!qmom_options.adaptive && state.active_node_count() != maximum_node_count)
    {
        throw std::invalid_argument("Fixed-count QMOM phase-space evolution state must contain the configured node count.");
    }
    if (!source)
    {
        throw std::invalid_argument("A moment-source callback must be provided for phase-space evolution.");
    }
    validate_boundary(boundary);
}

[[nodiscard]] bool advance_qmom_advection_stage(const QmomState& state, const QmomState& rollback_state, const double t, const double time_step,
                                                const PhaseSpaceVelocity& velocity, const QmomPhaseSpaceEvolutionOptions& options,
                                                const bool final_advection, PhaseSpaceEvolutionDiagnostics& diagnostics, QmomState& accepted_state,
                                                std::optional<QmomPhaseSpaceEvolutionResult>& failure)
{
    internal::performance::record_advection_operator_call();
    const QmomNodeAdvectionResult advection = advance_nodes(state.quadrature(), t, time_step, velocity, options.node_options);
    record_advection_diagnostics(diagnostics, final_advection, advection.status, advection.diagnostics);
    if (!advection.ok())
    {
        failure =
            qmom_failure(rollback_state, PhaseSpaceEvolutionStatus::advection_failed,
                         final_advection ? PhaseSpaceEvolutionStage::second_advection_half_step : PhaseSpaceEvolutionStage::first_advection_half_step,
                         diagnostics, t, advection.message);
        return false;
    }

    try
    {
        internal::performance::record_moment_reconstruction_call();
        accepted_state =
            internal::QmomStateAccess::from_quadrature(qmom_from_quadrature(advection.quadrature), state.moments().size(), state.abscissa_bounds());
        if (final_advection)
        {
            diagnostics.final_reconstruction_diagnostics = successful_reconstruction(accepted_state.moments().size());
        }
        else
        {
            diagnostics.first_reconstruction_diagnostics = successful_reconstruction(accepted_state.moments().size());
        }
    }
    catch (const std::exception& error)
    {
        if (final_advection)
        {
            diagnostics.final_reconstruction_diagnostics = {false, state.moments().size(), error.what()};
        }
        else
        {
            diagnostics.first_reconstruction_diagnostics = {false, state.moments().size(), error.what()};
        }
        failure = qmom_failure(rollback_state, PhaseSpaceEvolutionStatus::moment_reconstruction_failed,
                               final_advection ? PhaseSpaceEvolutionStage::final_moment_reconstruction
                                               : PhaseSpaceEvolutionStage::first_moment_reconstruction,
                               diagnostics, t + time_step, error.what());
        return false;
    }
    return true;
}

[[nodiscard]] bool advance_qmom_moment_source_stage(const QmomState& state, const QmomState& rollback_state, const std::size_t maximum_node_count,
                                                    const double t, const double time_step, const MomentSource& source,
                                                    const MomentBoundaryContribution* const boundary, const QmomPhaseSpaceEvolutionOptions& options,
                                                    PhaseSpaceEvolutionDiagnostics& diagnostics, QmomState& accepted_state,
                                                    std::optional<QmomPhaseSpaceEvolutionResult>& failure)
{
    internal::performance::record_moment_source_operator_call();
    MomentSourceUpdateOptions source_options = options.source_options;
    const MomentStateConverter to_moment_set = [](const std::vector<double>& values) { return MomentSet(values); };
    if (options.realizability_policy.check_mode == RealizabilityCheckMode::FinalOnly)
    {
        const StateAcceptancePredicate policy_accept = make_realizability_acceptance_predicate(options.realizability_policy, to_moment_set);
        const StateAcceptancePredicate caller_accept = source_options.ssp_rk3_options.accept;
        source_options.ssp_rk3_options.accept = [caller_accept, policy_accept](const std::vector<double>& candidate)
        { return (!caller_accept || caller_accept(candidate)) && policy_accept(candidate); };
    }
    else if (options.realizability_policy.check_mode == RealizabilityCheckMode::EveryStage)
    {
        const StageAcceptancePredicate policy_accept = make_realizability_stage_acceptance_predicate(options.realizability_policy, to_moment_set);
        const StageAcceptancePredicate caller_accept = source_options.ssp_rk3_options.stage_accept;
        source_options.ssp_rk3_options.stage_accept = [caller_accept, policy_accept](const SspRk3Stage stage, const std::vector<double>& candidate)
        { return (!caller_accept || caller_accept(stage, candidate)) && policy_accept(stage, candidate); };
    }
    const MomentSourceUpdateResult source_update =
        boundary == nullptr
            ? internal::advance_trusted_qmom_moments_with_source(state.moments(), t, time_step, source, source_options, state.active_node_count(),
                                                                 options.qmom_options, state.abscissa_bounds())
            : internal::advance_trusted_qmom_moments_with_source(state.moments(), t, time_step, source, *boundary, source_options,
                                                                 state.active_node_count(), options.qmom_options, state.abscissa_bounds());
    record_source_update_diagnostics(diagnostics, source_update);
    if (!source_update.ok())
    {
        failure = qmom_failure(rollback_state, PhaseSpaceEvolutionStatus::moment_source_failed, PhaseSpaceEvolutionStage::moment_source_step,
                               diagnostics, t, source_update.message);
        return false;
    }

    internal::performance::record_post_source_inversion();
    QmomResult selected_quadrature =
        select_post_source_qmom(source_update.moments, state.active_node_count(), maximum_node_count, options.qmom_options, state.abscissa_bounds());
    diagnostics.post_source_qmom_status = selected_quadrature.status;
    if (!selected_quadrature.ok())
    {
        failure = qmom_failure(rollback_state, PhaseSpaceEvolutionStatus::inversion_failed, PhaseSpaceEvolutionStage::post_source_inversion,
                               diagnostics, t, selected_quadrature.message);
        return false;
    }
    try
    {
        internal::performance::record_moment_reconstruction_call();
        accepted_state = internal::QmomStateAccess::from_quadrature(std::move(selected_quadrature), state.moments().size(), state.abscissa_bounds());
    }
    catch (const std::exception& error)
    {
        std::string message = std::string("QMOM post-source moment reconstruction failed: ") + error.what();
        diagnostics.final_reconstruction_diagnostics = {false, state.moments().size(), message};
        failure = qmom_failure(rollback_state, PhaseSpaceEvolutionStatus::moment_reconstruction_failed,
                               PhaseSpaceEvolutionStage::post_source_inversion, diagnostics, t, std::move(message));
        return false;
    }
    return true;
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

[[nodiscard]] PhaseSpaceEvolutionDiagnostics neutral_cqmom_diagnostics(PhaseSpaceEvolutionDiagnostics diagnostics)
{
    diagnostics.message = neutral_cqmom_message(std::move(diagnostics.message));
    diagnostics.moment_source_diagnostics.message = neutral_cqmom_message(std::move(diagnostics.moment_source_diagnostics.message));
    diagnostics.first_reconstruction_diagnostics.message = neutral_cqmom_message(std::move(diagnostics.first_reconstruction_diagnostics.message));
    diagnostics.final_reconstruction_diagnostics.message = neutral_cqmom_message(std::move(diagnostics.final_reconstruction_diagnostics.message));
    return diagnostics;
}

[[nodiscard]] CqmomPhaseSpaceEvolutionResult neutral_cqmom_result(CqmomPhaseSpaceEvolutionResult result)
{
    result.message = neutral_cqmom_message(std::move(result.message));
    result.diagnostics = neutral_cqmom_diagnostics(std::move(result.diagnostics));
    return result;
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

void record_cqmom_node_count_transitions(PhaseSpaceEvolutionDiagnostics& diagnostics, const PhaseSpaceEvolutionStage stage,
                                         const std::vector<CqmomActiveBranchCount>& before_counts,
                                         const std::vector<CqmomActiveBranchCount>& after_counts)
{
    if (!before_counts.empty() || !after_counts.empty())
    {
        internal::performance::record_active_count_transition_build();
    }

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
// operator-split phase-space evolution sequence to an accepted one-dimensional state.
static QmomPhaseSpaceEvolutionResult advance_qmom_phase_space_impl(const QmomState& state, const std::size_t maximum_node_count, const double t,
                                                                   const double time_step, const PhaseSpaceVelocity& velocity,
                                                                   const MomentSource& source, const MomentBoundaryContribution* const boundary,
                                                                   const QmomPhaseSpaceEvolutionOptions& options)
{
    internal::performance::record_phase_space_evolution_macro_step();
    PhaseSpaceEvolutionDiagnostics diagnostics;
    diagnostics.scheme = options.scheme;
    try
    {
        validate_qmom_input(state, maximum_node_count, options.qmom_options, t, time_step, velocity, source, options.scheme);
        validate_boundary(boundary);
    }
    catch (const std::exception& error)
    {
        return qmom_failure(state, PhaseSpaceEvolutionStatus::invalid_input, PhaseSpaceEvolutionStage::none, diagnostics, t, error.what());
    }

    QmomState current_state = state;
    std::optional<QmomPhaseSpaceEvolutionResult> failure;

    const auto apply_advection = [&](const double start_time, const double step, const bool final_advection)
    {
        return advance_qmom_advection_stage(current_state, state, start_time, step, velocity, options, final_advection, diagnostics, current_state,
                                            failure);
    };

    const auto apply_source = [&](const double start_time)
    {
        return advance_qmom_moment_source_stage(current_state, state, maximum_node_count, start_time, time_step, source, boundary, options,
                                                diagnostics, current_state, failure);
    };

    if (!execute_phase_space_evolution_sequence(options.scheme, t, time_step, apply_advection, apply_source))
    {
        return std::move(*failure);
    }

    finalize_diagnostics(diagnostics);
    return QmomPhaseSpaceEvolutionResult{std::move(current_state), PhaseSpaceEvolutionStatus::success, std::move(diagnostics), t + time_step, ""};
}

QmomPhaseSpaceEvolutionResult advance_qmom_advection(const QmomState& state, const double t, const double time_step,
                                                     const PhaseSpaceVelocity& velocity, const QmomPhaseSpaceEvolutionOptions& options)
{
    PhaseSpaceEvolutionDiagnostics diagnostics;
    diagnostics.scheme = options.scheme;
    try
    {
        validate_qmom_advection_input(state, t, time_step, velocity);
    }
    catch (const std::exception& error)
    {
        return qmom_failure(state, PhaseSpaceEvolutionStatus::invalid_input, PhaseSpaceEvolutionStage::none, diagnostics, t, error.what());
    }

    QmomState accepted_state = state;
    std::optional<QmomPhaseSpaceEvolutionResult> failure;
    if (!advance_qmom_advection_stage(state, state, t, time_step, velocity, options, false, diagnostics, accepted_state, failure))
    {
        return std::move(*failure);
    }

    finalize_diagnostics(diagnostics);
    return QmomPhaseSpaceEvolutionResult{std::move(accepted_state), PhaseSpaceEvolutionStatus::success, std::move(diagnostics), t + time_step, ""};
}

static QmomPhaseSpaceEvolutionResult advance_qmom_moment_source_impl(const QmomState& state, const std::size_t maximum_node_count, const double t,
                                                                     const double time_step, const MomentSource& source,
                                                                     const MomentBoundaryContribution* const boundary,
                                                                     const QmomPhaseSpaceEvolutionOptions& options)
{
    PhaseSpaceEvolutionDiagnostics diagnostics;
    diagnostics.scheme = options.scheme;
    try
    {
        validate_qmom_moment_source_input(state, maximum_node_count, options.qmom_options, t, time_step, source, boundary);
    }
    catch (const std::exception& error)
    {
        return qmom_failure(state, PhaseSpaceEvolutionStatus::invalid_input, PhaseSpaceEvolutionStage::none, diagnostics, t, error.what());
    }

    QmomState accepted_state = state;
    std::optional<QmomPhaseSpaceEvolutionResult> failure;
    if (!advance_qmom_moment_source_stage(state, state, maximum_node_count, t, time_step, source, boundary, options, diagnostics, accepted_state,
                                          failure))
    {
        return std::move(*failure);
    }

    finalize_diagnostics(diagnostics);
    return QmomPhaseSpaceEvolutionResult{std::move(accepted_state), PhaseSpaceEvolutionStatus::success, std::move(diagnostics), t + time_step, ""};
}

QmomPhaseSpaceEvolutionResult advance_qmom_moment_source(const QmomState& state, const std::size_t maximum_node_count, const double t,
                                                         const double time_step, const MomentSource& source,
                                                         const QmomPhaseSpaceEvolutionOptions& options)
{
    return advance_qmom_moment_source_impl(state, maximum_node_count, t, time_step, source, nullptr, options);
}

QmomPhaseSpaceEvolutionResult advance_qmom_moment_source(const QmomState& state, const std::size_t maximum_node_count, const double t,
                                                         const double time_step, const MomentSource& source,
                                                         const MomentBoundaryContribution& boundary, const QmomPhaseSpaceEvolutionOptions& options)
{
    return advance_qmom_moment_source_impl(state, maximum_node_count, t, time_step, source, &boundary, options);
}

QmomPhaseSpaceEvolutionResult advance_qmom_phase_space(const QmomState& state, const std::size_t maximum_node_count, const double t,
                                                       const double time_step, const PhaseSpaceVelocity& velocity, const MomentSource& source,
                                                       const QmomPhaseSpaceEvolutionOptions& options)
{
    return advance_qmom_phase_space_impl(state, maximum_node_count, t, time_step, velocity, source, nullptr, options);
}

QmomPhaseSpaceEvolutionResult advance_qmom_phase_space(const QmomState& state, const std::size_t maximum_node_count, const double t,
                                                       const double time_step, const PhaseSpaceVelocity& velocity, const MomentSource& source,
                                                       const MomentBoundaryContribution& boundary, const QmomPhaseSpaceEvolutionOptions& options)
{
    return advance_qmom_phase_space_impl(state, maximum_node_count, t, time_step, velocity, source, &boundary, options);
}

namespace
{

void initialize_cqmom_diagnostics(PhaseSpaceEvolutionDiagnostics& diagnostics, const CqmomDenseCase& dense_case)
{
    diagnostics.cqmom_moment_storage = CqmomMomentStorage::dense;
    diagnostics.cqmom_diagnostic_level = dense_case.diagnostic_level;
    diagnostics.cqmom_vandermonde_solver = dense_case.vandermonde_solver;
}

void initialize_cqmom_diagnostics(PhaseSpaceEvolutionDiagnostics& diagnostics, const CqmomLeanCase& lean_case)
{
    diagnostics.cqmom_moment_storage = CqmomMomentStorage::lean;
    diagnostics.cqmom_diagnostic_level = lean_case.diagnostic_level;
    diagnostics.cqmom_vandermonde_solver = lean_case.vandermonde_solver;
}

void initialize_cqmom_diagnostics(PhaseSpaceEvolutionDiagnostics& diagnostics, const CqmomState& state)
{
    if (const CqmomDenseCase* dense_case = internal::CqmomConstructionAccess::state_dense_case(state))
    {
        initialize_cqmom_diagnostics(diagnostics, *dense_case);
        return;
    }
    if (const CqmomLeanCase* lean_case = internal::CqmomConstructionAccess::state_lean_case(state))
    {
        initialize_cqmom_diagnostics(diagnostics, *lean_case);
    }
}

void validate_dense_cqmom_moments_for_phase_space_evolution(const CqmomDenseCase& dense_case)
{
    internal::validate_dense_cqmom_moment_tensor(
        dense_case.moments, dense_case.hierarchy, "CQMOM phase-space evolution moment tensor rank must match the phase-space hierarchy.",
        "CQMOM phase-space evolution moment tensor size does not match its extents.", "CQMOM phase-space evolution moments must be finite.");
}

void validate_lean_cqmom_moments_for_phase_space_evolution(const CqmomLeanCase& lean_case)
{
    internal::validate_cqmom_moment_values_match_layout(lean_case.moments,
                                                        {"", "Lean CQMOM phase-space evolution moments must match the fixed layout size.",
                                                         "Lean CQMOM phase-space evolution moments must be finite."});
}

void validate_cqmom_state_for_phase_space_evolution(const CqmomState& state)
{
    if (const CqmomDenseCase* dense_case = internal::CqmomConstructionAccess::state_dense_case(state))
    {
        validate_dense_cqmom_moments_for_phase_space_evolution(*dense_case);
        return;
    }
    if (const CqmomLeanCase* lean_case = internal::CqmomConstructionAccess::state_lean_case(state))
    {
        validate_lean_cqmom_moments_for_phase_space_evolution(*lean_case);
        return;
    }
    throw std::invalid_argument("Unknown CQMOM state storage.");
}

void validate_cqmom_advection_input(const CqmomState& state, const double t, const double time_step, const PhaseSpaceVelocity& velocity)
{
    validate_time_input(t, time_step);
    validate_cqmom_state_for_phase_space_evolution(state);
    if (!velocity)
    {
        throw std::invalid_argument("A phase-space velocity callback must be provided for phase-space evolution.");
    }
}

void validate_cqmom_moment_source_input(const CqmomState& state, const double t, const double time_step, const CqmomMomentSource& source,
                                        const CqmomMomentBoundaryContribution* const boundary)
{
    validate_time_input(t, time_step);
    validate_cqmom_state_for_phase_space_evolution(state);
    if (!source)
    {
        throw std::invalid_argument("CQMOM source callback must be provided.");
    }
    validate_cqmom_boundary(boundary);
}

void validate_cqmom_phase_space_evolution_input(const CqmomState& state, const double t, const double time_step, const PhaseSpaceVelocity& velocity,
                                                const CqmomMomentSource& source, const CqmomMomentBoundaryContribution* const boundary,
                                                const PhaseSpaceEvolutionScheme scheme)
{
    validate_time_input(t, time_step);
    validate_phase_space_evolution_scheme(scheme);
    validate_cqmom_state_for_phase_space_evolution(state);
    if (!velocity)
    {
        throw std::invalid_argument("A phase-space velocity callback must be provided for phase-space evolution.");
    }
    if (!source)
    {
        throw std::invalid_argument("CQMOM source callback must be provided.");
    }
    validate_cqmom_boundary(boundary);
}

struct CqmomStageOutcome
{
    CqmomState state;
    std::optional<CqmomQuadrature> diagnostic_quadrature;
    PhaseSpaceEvolutionStatus status = PhaseSpaceEvolutionStatus::success;
    PhaseSpaceEvolutionDiagnostics diagnostics;
    double final_time = 0.0;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == PhaseSpaceEvolutionStatus::success;
    }
};

struct CqmomPhaseSpaceEvolutionWorkingState
{
    std::optional<CqmomQuadrature> reconstructed_quadrature;
    const CqmomQuadrature* accepted_quadrature = nullptr;
    std::vector<CqmomActiveBranchCount> reconstructed_active_branch_counts;
    const std::vector<CqmomActiveBranchCount>* active_branch_counts = nullptr;
};

[[nodiscard]] std::vector<CqmomActiveBranchCount> take_cqmom_advection_active_counts(CqmomPhaseSpaceEvolutionWorkingState& working_state)
{
    if (working_state.reconstructed_quadrature.has_value())
    {
        return std::move(working_state.reconstructed_active_branch_counts);
    }
    return *working_state.active_branch_counts;
}

[[nodiscard]] CqmomStageOutcome cqmom_state_failure(const CqmomState& rollback_state, std::optional<CqmomQuadrature> quadrature,
                                                    const PhaseSpaceEvolutionStatus status, const PhaseSpaceEvolutionStage stage,
                                                    PhaseSpaceEvolutionDiagnostics diagnostics, const double final_time, std::string message)
{
    record_failure_diagnostics(diagnostics, stage, message);
    return CqmomStageOutcome{rollback_state, std::move(quadrature), status, std::move(diagnostics), final_time, std::move(message)};
}

[[nodiscard]] CqmomStageOutcome cqmom_success(CqmomState state, PhaseSpaceEvolutionDiagnostics diagnostics, const double final_time)
{
    return CqmomStageOutcome{std::move(state), std::nullopt, PhaseSpaceEvolutionStatus::success, std::move(diagnostics), final_time, ""};
}

[[nodiscard]] CqmomPhaseSpaceEvolutionResult materialize_cqmom_phase_space_evolution_result(CqmomStageOutcome outcome)
{
    if (!outcome.ok())
    {
        return neutral_cqmom_result(CqmomPhaseSpaceEvolutionResult{std::move(outcome.state), std::move(outcome.diagnostic_quadrature), outcome.status,
                                                                   std::move(outcome.diagnostics), outcome.final_time, std::move(outcome.message)});
    }

    const std::optional<CqmomQuadrature>& accepted_quadrature = internal::CqmomConstructionAccess::state_accepted_quadrature(outcome.state);
    if (!accepted_quadrature.has_value())
    {
        throw std::logic_error("Successful CQMOM phase-space evolution state has no accepted quadrature cache.");
    }

    CqmomQuadrature public_quadrature = *accepted_quadrature;
    finalize_diagnostics(outcome.diagnostics);
    return neutral_cqmom_result(CqmomPhaseSpaceEvolutionResult{std::move(outcome.state), std::move(public_quadrature),
                                                               PhaseSpaceEvolutionStatus::success, std::move(outcome.diagnostics), outcome.final_time,
                                                               ""});
}

[[nodiscard]] CqmomState dense_state_from_accepted_quadrature(const CqmomDenseCase& base_case, RawMomentTensor moments,
                                                              std::vector<CqmomActiveBranchCount> active_counts, CqmomQuadrature&& quadrature,
                                                              const internal::PreparedFixedCqmomMetadataPtr& prepared_fixed_metadata)
{
    CqmomDenseCase updated_case{base_case.hierarchy,
                                std::move(moments),
                                base_case.conditional_options,
                                base_case.qmom_options,
                                base_case.adaptive_branch_counts,
                                base_case.diagnostic_level,
                                base_case.vandermonde_solver,
                                base_case.coordinate_support};
    return internal::CqmomConstructionAccess::create_state(std::move(updated_case), std::move(active_counts), std::move(quadrature),
                                                           prepared_fixed_metadata);
}

[[nodiscard]] CqmomState lean_state_from_accepted_quadrature(const CqmomLeanCase& base_case, CqmomMomentData moments,
                                                             std::vector<CqmomActiveBranchCount> active_counts, CqmomQuadrature&& quadrature,
                                                             const internal::PreparedFixedCqmomMetadataPtr& prepared_fixed_metadata)
{
    CqmomLeanCase updated_case{
        std::move(moments),         base_case.conditional_options, base_case.qmom_options,      base_case.adaptive_branch_counts,
        base_case.diagnostic_level, base_case.vandermonde_solver,  base_case.coordinate_support};
    return internal::CqmomConstructionAccess::create_state(std::move(updated_case), std::move(active_counts), std::move(quadrature),
                                                           prepared_fixed_metadata);
}

[[nodiscard]] bool prepare_cqmom_phase_space_evolution_working_state(const CqmomState& state, const CqmomState& rollback_state, const double t,
                                                                     PhaseSpaceEvolutionDiagnostics& diagnostics,
                                                                     CqmomPhaseSpaceEvolutionWorkingState& working_state,
                                                                     std::optional<CqmomStageOutcome>& failure)
{
    const std::optional<CqmomQuadrature>& accepted_quadrature = internal::CqmomConstructionAccess::state_accepted_quadrature(state);
    if (accepted_quadrature.has_value())
    {
        working_state.accepted_quadrature = &*accepted_quadrature;
        working_state.active_branch_counts = &internal::CqmomConstructionAccess::state_active_branch_counts(state);
        return true;
    }

    const std::vector<CqmomActiveBranchCount>& before_initial_active_branch_counts =
        internal::CqmomConstructionAccess::state_active_branch_counts(state);
    working_state.reconstructed_active_branch_counts = before_initial_active_branch_counts;
    const CqmomDenseCase* dense_case = internal::CqmomConstructionAccess::state_dense_case(state);
    const internal::PreparedFixedCqmomMetadata* prepared_fixed_metadata =
        internal::CqmomConstructionAccess::state_prepared_fixed_metadata(state).get();
    CqmomResult initial_quadrature = dense_case != nullptr
                                         ? internal::compute_cqmom_internal(*dense_case, working_state.reconstructed_active_branch_counts, false,
                                                                            prepared_fixed_metadata, !dense_case->adaptive_branch_counts)
                                         : internal::compute_cqmom_internal(*internal::CqmomConstructionAccess::state_lean_case(state),
                                                                            working_state.reconstructed_active_branch_counts, false,
                                                                            prepared_fixed_metadata, !state.adaptive_branch_counts());
    diagnostics.initial_cqmom_status = initial_quadrature.status;
    diagnostics.initial_cqmom_diagnostics = initial_quadrature.diagnostics;
    if (!initial_quadrature.ok())
    {
        failure = cqmom_state_failure(rollback_state, std::move(initial_quadrature.quadrature), PhaseSpaceEvolutionStatus::inversion_failed,
                                      PhaseSpaceEvolutionStage::initial_inversion, diagnostics, t,
                                      neutral_cqmom_message(std::move(initial_quadrature.message)));
        return false;
    }

    working_state.reconstructed_quadrature = std::move(initial_quadrature.quadrature);
    working_state.accepted_quadrature = &*working_state.reconstructed_quadrature;
    working_state.reconstructed_active_branch_counts = internal::active_branch_counts_from_quadrature(*working_state.accepted_quadrature);
    working_state.active_branch_counts = &working_state.reconstructed_active_branch_counts;
    record_cqmom_node_count_transitions(diagnostics, PhaseSpaceEvolutionStage::initial_inversion, before_initial_active_branch_counts,
                                        *working_state.active_branch_counts);
    return true;
}

[[nodiscard]] CqmomStageOutcome advance_cqmom_advection_stage(const CqmomState& state, const CqmomState& rollback_state, const double t,
                                                              const double time_step, const PhaseSpaceVelocity& velocity,
                                                              const CqmomPhaseSpaceEvolutionOptions& options, const bool final_advection,
                                                              PhaseSpaceEvolutionDiagnostics& diagnostics)
{
    CqmomPhaseSpaceEvolutionWorkingState working_state;
    std::optional<CqmomStageOutcome> failure;
    if (!prepare_cqmom_phase_space_evolution_working_state(state, rollback_state, t, diagnostics, working_state, failure))
    {
        return std::move(*failure);
    }

    internal::performance::record_advection_operator_call();
    CqmomNodeAdvectionResult advection = advance_nodes(*working_state.accepted_quadrature, t, time_step, velocity, options.node_options);
    record_advection_diagnostics(diagnostics, final_advection, advection.status, advection.diagnostics);
    if (!advection.ok())
    {
        return cqmom_state_failure(rollback_state, *working_state.accepted_quadrature, PhaseSpaceEvolutionStatus::advection_failed,
                                   final_advection ? PhaseSpaceEvolutionStage::second_advection_half_step
                                                   : PhaseSpaceEvolutionStage::first_advection_half_step,
                                   diagnostics, t, advection.message);
    }

    if (const CqmomDenseCase* dense_case = internal::CqmomConstructionAccess::state_dense_case(state))
    {
        const internal::PreparedFixedCqmomMetadataPtr& prepared_fixed_metadata =
            internal::CqmomConstructionAccess::state_prepared_fixed_metadata(state);
        try
        {
            internal::performance::record_moment_reconstruction_call();
            RawMomentTensor moments = reconstruct_raw_moments(advection.quadrature, dense_case->moments.extents);
            if (final_advection)
            {
                diagnostics.final_reconstruction_diagnostics = successful_reconstruction(moments.values.size());
            }
            else
            {
                diagnostics.first_reconstruction_diagnostics = successful_reconstruction(moments.values.size());
            }
            std::vector<CqmomActiveBranchCount> active_counts = take_cqmom_advection_active_counts(working_state);
            return cqmom_success(dense_state_from_accepted_quadrature(*dense_case, std::move(moments), std::move(active_counts),
                                                                      std::move(advection.quadrature), prepared_fixed_metadata),
                                 diagnostics, t + time_step);
        }
        catch (const std::exception& error)
        {
            if (final_advection)
            {
                diagnostics.final_reconstruction_diagnostics = {false, dense_case->moments.values.size(), error.what()};
            }
            else
            {
                diagnostics.first_reconstruction_diagnostics = {false, dense_case->moments.values.size(), error.what()};
            }
            return cqmom_state_failure(rollback_state, advection.quadrature, PhaseSpaceEvolutionStatus::moment_reconstruction_failed,
                                       final_advection ? PhaseSpaceEvolutionStage::final_moment_reconstruction
                                                       : PhaseSpaceEvolutionStage::first_moment_reconstruction,
                                       diagnostics, t + time_step, error.what());
        }
    }

    const CqmomLeanCase* lean_case = internal::CqmomConstructionAccess::state_lean_case(state);
    const internal::PreparedFixedCqmomMetadataPtr& prepared_fixed_metadata = internal::CqmomConstructionAccess::state_prepared_fixed_metadata(state);
    try
    {
        internal::performance::record_moment_reconstruction_call();
        CqmomMomentData moments = internal::reconstruct_trusted_cqmom_moments(advection.quadrature, lean_case->moments.layout);
        if (final_advection)
        {
            diagnostics.final_reconstruction_diagnostics = successful_reconstruction(moments.values.size());
        }
        else
        {
            diagnostics.first_reconstruction_diagnostics = successful_reconstruction(moments.values.size());
        }
        std::vector<CqmomActiveBranchCount> active_counts = take_cqmom_advection_active_counts(working_state);
        return cqmom_success(lean_state_from_accepted_quadrature(*lean_case, std::move(moments), std::move(active_counts),
                                                                 std::move(advection.quadrature), prepared_fixed_metadata),
                             diagnostics, t + time_step);
    }
    catch (const std::exception& error)
    {
        if (final_advection)
        {
            diagnostics.final_reconstruction_diagnostics = {false, lean_case->moments.values.size(), error.what()};
        }
        else
        {
            diagnostics.first_reconstruction_diagnostics = {false, lean_case->moments.values.size(), error.what()};
        }
        return cqmom_state_failure(rollback_state, advection.quadrature, PhaseSpaceEvolutionStatus::moment_reconstruction_failed,
                                   final_advection ? PhaseSpaceEvolutionStage::final_moment_reconstruction
                                                   : PhaseSpaceEvolutionStage::first_moment_reconstruction,
                                   diagnostics, t + time_step, error.what());
    }
}

[[nodiscard]] CqmomStageOutcome advance_cqmom_moment_source_stage(const CqmomState& state, const CqmomState& rollback_state, const double t,
                                                                  const double time_step, const CqmomMomentSource& source,
                                                                  const CqmomMomentBoundaryContribution* const boundary,
                                                                  const CqmomPhaseSpaceEvolutionOptions& options,
                                                                  PhaseSpaceEvolutionDiagnostics& diagnostics)
{
    CqmomPhaseSpaceEvolutionWorkingState working_state;
    std::optional<CqmomStageOutcome> failure;
    if (!prepare_cqmom_phase_space_evolution_working_state(state, rollback_state, t, diagnostics, working_state, failure))
    {
        return std::move(*failure);
    }

    internal::performance::record_moment_source_operator_call();
    if (const CqmomDenseCase* dense_case = internal::CqmomConstructionAccess::state_dense_case(state))
    {
        const internal::PreparedFixedCqmomMetadataPtr& prepared_fixed_metadata =
            internal::CqmomConstructionAccess::state_prepared_fixed_metadata(state);
        std::optional<CqmomMomentLayout> generated_source_layout;
        const CqmomMomentLayout* source_layout = prepared_fixed_metadata == nullptr ? nullptr : &prepared_fixed_metadata->moment_layout;
        if (source_layout == nullptr)
        {
            generated_source_layout.emplace(create_cqmom_moment_layout(dense_case->hierarchy));
            source_layout = &*generated_source_layout;
        }
        CqmomMomentData source_moments = extract_cqmom_moments(dense_case->moments, *source_layout);
        const CqmomLeanCase source_reconstruction_case{source_moments,
                                                       dense_case->conditional_options,
                                                       dense_case->qmom_options,
                                                       dense_case->adaptive_branch_counts,
                                                       dense_case->diagnostic_level,
                                                       dense_case->vandermonde_solver,
                                                       dense_case->coordinate_support};
        CqmomMomentSourceUpdateResult source_update =
            boundary == nullptr
                ? internal::advance_trusted_cqmom_moments_with_source(source_moments, t, time_step, source, options.source_options,
                                                                      source_reconstruction_case, *working_state.active_branch_counts,
                                                                      prepared_fixed_metadata.get())
                : internal::advance_trusted_cqmom_moments_with_source(source_moments, t, time_step, source, *boundary, options.source_options,
                                                                      source_reconstruction_case, *working_state.active_branch_counts,
                                                                      prepared_fixed_metadata.get());
        record_source_update_diagnostics(diagnostics, source_update);
        if (!source_update.ok())
        {
            return cqmom_state_failure(rollback_state, *working_state.accepted_quadrature, PhaseSpaceEvolutionStatus::moment_source_failed,
                                       PhaseSpaceEvolutionStage::moment_source_step, diagnostics, t, source_update.message);
        }

        const std::vector<CqmomActiveBranchCount>& before_source_active_branch_counts = *working_state.active_branch_counts;
        internal::performance::record_post_source_inversion();
        CqmomLeanCase post_source_case{std::move(source_update.moments),   dense_case->conditional_options, dense_case->qmom_options,
                                       dense_case->adaptive_branch_counts, dense_case->diagnostic_level,    dense_case->vandermonde_solver,
                                       dense_case->coordinate_support};
        CqmomResult post_source_quadrature =
            internal::compute_cqmom_internal(post_source_case, *working_state.active_branch_counts, dense_case->adaptive_branch_counts,
                                             prepared_fixed_metadata.get(), !dense_case->adaptive_branch_counts);
        diagnostics.post_source_cqmom_status = post_source_quadrature.status;
        diagnostics.post_source_cqmom_diagnostics = post_source_quadrature.diagnostics;
        if (!post_source_quadrature.ok())
        {
            return cqmom_state_failure(rollback_state, post_source_quadrature.quadrature, PhaseSpaceEvolutionStatus::inversion_failed,
                                       PhaseSpaceEvolutionStage::post_source_inversion, diagnostics, t,
                                       neutral_cqmom_message(std::move(post_source_quadrature.message)));
        }

        try
        {
            internal::performance::record_moment_reconstruction_call();
            RawMomentTensor moments = reconstruct_raw_moments(post_source_quadrature.quadrature, dense_case->moments.extents);
            std::vector<CqmomActiveBranchCount> after_source_active_branch_counts =
                internal::active_branch_counts_from_quadrature(post_source_quadrature.quadrature);
            record_cqmom_node_count_transitions(diagnostics, PhaseSpaceEvolutionStage::post_source_inversion, before_source_active_branch_counts,
                                                after_source_active_branch_counts);
            return cqmom_success(dense_state_from_accepted_quadrature(*dense_case, std::move(moments), std::move(after_source_active_branch_counts),
                                                                      std::move(post_source_quadrature.quadrature), prepared_fixed_metadata),
                                 diagnostics, t + time_step);
        }
        catch (const std::exception& error)
        {
            std::string message = std::string("Dense CQMOM post-source moment reconstruction failed: ") + error.what();
            diagnostics.final_reconstruction_diagnostics = {false, dense_case->moments.values.size(), message};
            return cqmom_state_failure(rollback_state, post_source_quadrature.quadrature, PhaseSpaceEvolutionStatus::moment_reconstruction_failed,
                                       PhaseSpaceEvolutionStage::post_source_inversion, diagnostics, t, std::move(message));
        }
    }

    const CqmomLeanCase* lean_case = internal::CqmomConstructionAccess::state_lean_case(state);
    const internal::PreparedFixedCqmomMetadataPtr& prepared_fixed_metadata = internal::CqmomConstructionAccess::state_prepared_fixed_metadata(state);
    CqmomMomentSourceUpdateResult source_update =
        boundary == nullptr
            ? internal::advance_trusted_cqmom_moments_with_source(lean_case->moments, t, time_step, source, options.source_options, *lean_case,
                                                                  *working_state.active_branch_counts, prepared_fixed_metadata.get())
            : internal::advance_trusted_cqmom_moments_with_source(lean_case->moments, t, time_step, source, *boundary, options.source_options,
                                                                  *lean_case, *working_state.active_branch_counts, prepared_fixed_metadata.get());
    record_source_update_diagnostics(diagnostics, source_update);
    if (!source_update.ok())
    {
        return cqmom_state_failure(rollback_state, *working_state.accepted_quadrature, PhaseSpaceEvolutionStatus::moment_source_failed,
                                   PhaseSpaceEvolutionStage::moment_source_step, diagnostics, t, source_update.message);
    }

    const std::vector<CqmomActiveBranchCount>& before_source_active_branch_counts = *working_state.active_branch_counts;
    internal::performance::record_post_source_inversion();
    CqmomLeanCase post_source_case{std::move(source_update.moments),  lean_case->conditional_options, lean_case->qmom_options,
                                   lean_case->adaptive_branch_counts, lean_case->diagnostic_level,    lean_case->vandermonde_solver,
                                   lean_case->coordinate_support};
    CqmomResult post_source_quadrature =
        internal::compute_cqmom_internal(post_source_case, *working_state.active_branch_counts, lean_case->adaptive_branch_counts,
                                         prepared_fixed_metadata.get(), !lean_case->adaptive_branch_counts);
    diagnostics.post_source_cqmom_status = post_source_quadrature.status;
    diagnostics.post_source_cqmom_diagnostics = post_source_quadrature.diagnostics;
    if (!post_source_quadrature.ok())
    {
        return cqmom_state_failure(rollback_state, post_source_quadrature.quadrature, PhaseSpaceEvolutionStatus::inversion_failed,
                                   PhaseSpaceEvolutionStage::post_source_inversion, diagnostics, t,
                                   neutral_cqmom_message(std::move(post_source_quadrature.message)));
    }

    try
    {
        internal::performance::record_moment_reconstruction_call();
        CqmomMomentData moments = internal::reconstruct_trusted_cqmom_moments(post_source_quadrature.quadrature, lean_case->moments.layout);
        std::vector<CqmomActiveBranchCount> after_source_active_branch_counts =
            internal::active_branch_counts_from_quadrature(post_source_quadrature.quadrature);
        record_cqmom_node_count_transitions(diagnostics, PhaseSpaceEvolutionStage::post_source_inversion, before_source_active_branch_counts,
                                            after_source_active_branch_counts);
        return cqmom_success(lean_state_from_accepted_quadrature(*lean_case, std::move(moments), std::move(after_source_active_branch_counts),
                                                                 std::move(post_source_quadrature.quadrature), prepared_fixed_metadata),
                             diagnostics, t + time_step);
    }
    catch (const std::exception& error)
    {
        std::string message = std::string("Lean CQMOM post-source moment reconstruction failed: ") + error.what();
        diagnostics.final_reconstruction_diagnostics = {false, lean_case->moments.values.size(), message};
        return cqmom_state_failure(rollback_state, post_source_quadrature.quadrature, PhaseSpaceEvolutionStatus::moment_reconstruction_failed,
                                   PhaseSpaceEvolutionStage::post_source_inversion, diagnostics, t, std::move(message));
    }
}

CqmomStageOutcome advance_cqmom_advection_impl(const CqmomState& state, const CqmomState& rollback_state, const double t, const double time_step,
                                               const PhaseSpaceVelocity& velocity, const CqmomPhaseSpaceEvolutionOptions& options,
                                               const bool final_advection, PhaseSpaceEvolutionDiagnostics& diagnostics)
{
    return advance_cqmom_advection_stage(state, rollback_state, t, time_step, velocity, options, final_advection, diagnostics);
}

CqmomStageOutcome advance_cqmom_moment_source_impl(const CqmomState& state, const CqmomState& rollback_state, const double t, const double time_step,
                                                   const CqmomMomentSource& source, const CqmomMomentBoundaryContribution* const boundary,
                                                   const CqmomPhaseSpaceEvolutionOptions& options, PhaseSpaceEvolutionDiagnostics& diagnostics)
{
    return advance_cqmom_moment_source_stage(state, rollback_state, t, time_step, source, boundary, options, diagnostics);
}

CqmomPhaseSpaceEvolutionResult advance_cqmom_phase_space_state_impl(const CqmomState& state, const double t, const double time_step,
                                                                    const PhaseSpaceVelocity& velocity, const CqmomMomentSource& source,
                                                                    const CqmomMomentBoundaryContribution* const boundary,
                                                                    const CqmomPhaseSpaceEvolutionOptions& options)
{
    internal::performance::record_phase_space_evolution_macro_step();
    PhaseSpaceEvolutionDiagnostics diagnostics;
    diagnostics.scheme = options.scheme;
    initialize_cqmom_diagnostics(diagnostics, state);
    try
    {
        validate_cqmom_phase_space_evolution_input(state, t, time_step, velocity, source, boundary, options.scheme);
    }
    catch (const std::exception& error)
    {
        return materialize_cqmom_phase_space_evolution_result(cqmom_state_failure(state, std::nullopt, PhaseSpaceEvolutionStatus::invalid_input,
                                                                                  PhaseSpaceEvolutionStage::none, diagnostics, t, error.what()));
    }

    const CqmomState* current_state = &state;
    std::optional<CqmomState> owned_current_state;
    std::optional<CqmomStageOutcome> failure;
    const auto apply_advection = [&](const double start_time, const double step, const bool final_advection)
    {
        CqmomStageOutcome result =
            advance_cqmom_advection_impl(*current_state, state, start_time, step, velocity, options, final_advection, diagnostics);
        if (!result.ok())
        {
            failure = std::move(result);
            return false;
        }
        owned_current_state = std::move(result.state);
        current_state = &*owned_current_state;
        return true;
    };

    const auto apply_source = [&](const double start_time)
    {
        CqmomStageOutcome result =
            advance_cqmom_moment_source_impl(*current_state, state, start_time, time_step, source, boundary, options, diagnostics);
        if (!result.ok())
        {
            failure = std::move(result);
            return false;
        }
        owned_current_state = std::move(result.state);
        current_state = &*owned_current_state;
        return true;
    };

    if (!execute_phase_space_evolution_sequence(options.scheme, t, time_step, apply_advection, apply_source))
    {
        return materialize_cqmom_phase_space_evolution_result(std::move(*failure));
    }

    return materialize_cqmom_phase_space_evolution_result(cqmom_success(std::move(*owned_current_state), std::move(diagnostics), t + time_step));
}

} // namespace

CqmomPhaseSpaceEvolutionResult advance_cqmom_advection(const CqmomState& state, const double t, const double time_step,
                                                       const PhaseSpaceVelocity& velocity, const CqmomPhaseSpaceEvolutionOptions& options)
{
    PhaseSpaceEvolutionDiagnostics diagnostics;
    diagnostics.scheme = options.scheme;
    initialize_cqmom_diagnostics(diagnostics, state);
    try
    {
        validate_cqmom_advection_input(state, t, time_step, velocity);
    }
    catch (const std::exception& error)
    {
        return materialize_cqmom_phase_space_evolution_result(cqmom_state_failure(state, std::nullopt, PhaseSpaceEvolutionStatus::invalid_input,
                                                                                  PhaseSpaceEvolutionStage::none, diagnostics, t, error.what()));
    }
    return materialize_cqmom_phase_space_evolution_result(
        advance_cqmom_advection_impl(state, state, t, time_step, velocity, options, false, diagnostics));
}

namespace
{

CqmomPhaseSpaceEvolutionResult advance_cqmom_moment_source_public_impl(const CqmomState& state, const double t, const double time_step,
                                                                       const CqmomMomentSource& source,
                                                                       const CqmomMomentBoundaryContribution* const boundary,
                                                                       const CqmomPhaseSpaceEvolutionOptions& options)
{
    PhaseSpaceEvolutionDiagnostics diagnostics;
    diagnostics.scheme = options.scheme;
    initialize_cqmom_diagnostics(diagnostics, state);
    try
    {
        validate_cqmom_moment_source_input(state, t, time_step, source, boundary);
    }
    catch (const std::exception& error)
    {
        return materialize_cqmom_phase_space_evolution_result(cqmom_state_failure(state, std::nullopt, PhaseSpaceEvolutionStatus::invalid_input,
                                                                                  PhaseSpaceEvolutionStage::none, diagnostics, t, error.what()));
    }
    return materialize_cqmom_phase_space_evolution_result(
        advance_cqmom_moment_source_impl(state, state, t, time_step, source, boundary, options, diagnostics));
}

} // namespace

CqmomPhaseSpaceEvolutionResult advance_cqmom_moment_source(const CqmomState& state, const double t, const double time_step,
                                                           const CqmomMomentSource& source, const CqmomPhaseSpaceEvolutionOptions& options)
{
    return advance_cqmom_moment_source_public_impl(state, t, time_step, source, nullptr, options);
}

CqmomPhaseSpaceEvolutionResult advance_cqmom_moment_source(const CqmomState& state, const double t, const double time_step,
                                                           const CqmomMomentSource& source, const CqmomMomentBoundaryContribution& boundary,
                                                           const CqmomPhaseSpaceEvolutionOptions& options)
{
    return advance_cqmom_moment_source_public_impl(state, t, time_step, source, &boundary, options);
}

CqmomPhaseSpaceEvolutionResult advance_cqmom_phase_space(const CqmomState& state, const double t, const double time_step,
                                                         const PhaseSpaceVelocity& velocity, const CqmomMomentSource& source,
                                                         const CqmomPhaseSpaceEvolutionOptions& options)
{
    return advance_cqmom_phase_space_state_impl(state, t, time_step, velocity, source, nullptr, options);
}

CqmomPhaseSpaceEvolutionResult advance_cqmom_phase_space(const CqmomState& state, const double t, const double time_step,
                                                         const PhaseSpaceVelocity& velocity, const CqmomMomentSource& source,
                                                         const CqmomMomentBoundaryContribution& boundary,
                                                         const CqmomPhaseSpaceEvolutionOptions& options)
{
    return advance_cqmom_phase_space_state_impl(state, t, time_step, velocity, source, &boundary, options);
}

std::vector<QmomPhaseSpaceEvolutionResult> advance_qmom_advection_batch(const std::vector<QmomState>& states, const double time,
                                                                        const double time_step, const PhaseSpaceVelocity& velocity,
                                                                        const QmomPhaseSpaceEvolutionOptions& options)
{
    return apply_batch(states,
                       [&](const std::size_t, const QmomState& state) { return advance_qmom_advection(state, time, time_step, velocity, options); });
}

std::vector<CqmomPhaseSpaceEvolutionResult> advance_cqmom_advection_batch(const std::vector<CqmomState>& states, const double time,
                                                                          const double time_step, const PhaseSpaceVelocity& velocity,
                                                                          const CqmomPhaseSpaceEvolutionOptions& options)
{
    return apply_batch(states, [&](const std::size_t, const CqmomState& state)
                       { return advance_cqmom_advection(state, time, time_step, velocity, options); });
}

std::vector<QmomPhaseSpaceEvolutionResult> advance_qmom_moment_source_batch(const std::vector<QmomState>& states,
                                                                            const std::size_t maximum_node_count, const double time,
                                                                            const double time_step, const MomentSource& source,
                                                                            const QmomPhaseSpaceEvolutionOptions& options)
{
    return apply_batch(states, [&](const std::size_t, const QmomState& state)
                       { return advance_qmom_moment_source(state, maximum_node_count, time, time_step, source, options); });
}

std::vector<CqmomPhaseSpaceEvolutionResult> advance_cqmom_moment_source_batch(const std::vector<CqmomState>& states, const double time,
                                                                              const double time_step, const CqmomMomentSource& source,
                                                                              const CqmomPhaseSpaceEvolutionOptions& options)
{
    return apply_batch(states, [&](const std::size_t, const CqmomState& state)
                       { return advance_cqmom_moment_source(state, time, time_step, source, options); });
}

std::vector<QmomPhaseSpaceEvolutionResult> advance_qmom_moment_source_batch(const std::vector<QmomState>& states,
                                                                            const std::size_t maximum_node_count, const double time,
                                                                            const double time_step, const MomentSource& source,
                                                                            const MomentBoundaryContribution& boundary,
                                                                            const QmomPhaseSpaceEvolutionOptions& options)
{
    return apply_batch(states, [&](const std::size_t, const QmomState& state)
                       { return advance_qmom_moment_source(state, maximum_node_count, time, time_step, source, boundary, options); });
}

std::vector<CqmomPhaseSpaceEvolutionResult> advance_cqmom_moment_source_batch(const std::vector<CqmomState>& states, const double time,
                                                                              const double time_step, const CqmomMomentSource& source,
                                                                              const CqmomMomentBoundaryContribution& boundary,
                                                                              const CqmomPhaseSpaceEvolutionOptions& options)
{
    return apply_batch(states, [&](const std::size_t, const CqmomState& state)
                       { return advance_cqmom_moment_source(state, time, time_step, source, boundary, options); });
}

std::vector<QmomPhaseSpaceEvolutionResult> advance_qmom_phase_space_batch(const std::vector<QmomState>& states, const std::size_t maximum_node_count,
                                                                          const double time, const double time_step,
                                                                          const PhaseSpaceVelocity& velocity, const MomentSource& source,
                                                                          const QmomPhaseSpaceEvolutionOptions& options)
{
    return apply_batch(states, [&](const std::size_t, const QmomState& state)
                       { return advance_qmom_phase_space(state, maximum_node_count, time, time_step, velocity, source, options); });
}

std::vector<CqmomPhaseSpaceEvolutionResult> advance_cqmom_phase_space_batch(const std::vector<CqmomState>& states, const double time,
                                                                            const double time_step, const PhaseSpaceVelocity& velocity,
                                                                            const CqmomMomentSource& source,
                                                                            const CqmomPhaseSpaceEvolutionOptions& options)
{
    return apply_batch(states, [&](const std::size_t, const CqmomState& state)
                       { return advance_cqmom_phase_space(state, time, time_step, velocity, source, options); });
}

std::vector<QmomPhaseSpaceEvolutionResult> advance_qmom_phase_space_batch(const std::vector<QmomState>& states, const std::size_t maximum_node_count,
                                                                          const double time, const double time_step,
                                                                          const PhaseSpaceVelocity& velocity, const MomentSource& source,
                                                                          const MomentBoundaryContribution& boundary,
                                                                          const QmomPhaseSpaceEvolutionOptions& options)
{
    return apply_batch(states, [&](const std::size_t, const QmomState& state)
                       { return advance_qmom_phase_space(state, maximum_node_count, time, time_step, velocity, source, boundary, options); });
}

std::vector<CqmomPhaseSpaceEvolutionResult> advance_cqmom_phase_space_batch(const std::vector<CqmomState>& states, const double time,
                                                                            const double time_step, const PhaseSpaceVelocity& velocity,
                                                                            const CqmomMomentSource& source,
                                                                            const CqmomMomentBoundaryContribution& boundary,
                                                                            const CqmomPhaseSpaceEvolutionOptions& options)
{
    return apply_batch(states, [&](const std::size_t, const CqmomState& state)
                       { return advance_cqmom_phase_space(state, time, time_step, velocity, source, boundary, options); });
}

std::vector<QmomPhaseSpaceEvolutionResult> advance_qmom_advection_batch_with_index(const std::vector<QmomState>& states, const double time,
                                                                                   const double time_step, const BatchPhaseSpaceVelocity& velocity,
                                                                                   const QmomPhaseSpaceEvolutionOptions& options)
{
    return apply_batch(
        states,
        [&](const std::size_t batch_index, const QmomState& state)
        {
            const PhaseSpaceVelocity local_velocity =
                velocity ? PhaseSpaceVelocity{[&velocity, batch_index](const double callback_time, const std::vector<double>& phase_space_location)
                                              { return velocity(batch_index, callback_time, phase_space_location); }}
                         : PhaseSpaceVelocity{};
            return advance_qmom_advection(state, time, time_step, local_velocity, options);
        });
}

std::vector<CqmomPhaseSpaceEvolutionResult> advance_cqmom_advection_batch_with_index(const std::vector<CqmomState>& states, const double time,
                                                                                     const double time_step, const BatchPhaseSpaceVelocity& velocity,
                                                                                     const CqmomPhaseSpaceEvolutionOptions& options)
{
    return apply_batch(
        states,
        [&](const std::size_t batch_index, const CqmomState& state)
        {
            const PhaseSpaceVelocity local_velocity =
                velocity ? PhaseSpaceVelocity{[&velocity, batch_index](const double callback_time, const std::vector<double>& phase_space_location)
                                              { return velocity(batch_index, callback_time, phase_space_location); }}
                         : PhaseSpaceVelocity{};
            return advance_cqmom_advection(state, time, time_step, local_velocity, options);
        });
}

std::vector<QmomPhaseSpaceEvolutionResult> advance_qmom_moment_source_batch_with_index(const std::vector<QmomState>& states,
                                                                                       const std::size_t maximum_node_count, const double time,
                                                                                       const double time_step, const BatchMomentSource& source,
                                                                                       const QmomPhaseSpaceEvolutionOptions& options)
{
    return apply_batch(states,
                       [&](const std::size_t batch_index, const QmomState& state)
                       {
                           const MomentSource local_source =
                               source ? MomentSource{[&source, batch_index](const double callback_time, const QmomSourceContext& context)
                                                     { return source(batch_index, callback_time, context); }}
                                      : MomentSource{};
                           return advance_qmom_moment_source(state, maximum_node_count, time, time_step, local_source, options);
                       });
}

std::vector<CqmomPhaseSpaceEvolutionResult> advance_cqmom_moment_source_batch_with_index(const std::vector<CqmomState>& states, const double time,
                                                                                         const double time_step, const BatchCqmomMomentSource& source,
                                                                                         const CqmomPhaseSpaceEvolutionOptions& options)
{
    return apply_batch(states,
                       [&](const std::size_t batch_index, const CqmomState& state)
                       {
                           const CqmomMomentSource local_source =
                               source ? CqmomMomentSource{[&source, batch_index](const double callback_time, const CqmomSourceContext& context)
                                                          { return source(batch_index, callback_time, context); }}
                                      : CqmomMomentSource{};
                           return advance_cqmom_moment_source(state, time, time_step, local_source, options);
                       });
}

std::vector<QmomPhaseSpaceEvolutionResult> advance_qmom_moment_source_batch_with_index(const std::vector<QmomState>& states,
                                                                                       const std::size_t maximum_node_count, const double time,
                                                                                       const double time_step, const BatchMomentSource& source,
                                                                                       const BatchMomentBoundaryContribution& boundary,
                                                                                       const QmomPhaseSpaceEvolutionOptions& options)
{
    return apply_batch(
        states,
        [&](const std::size_t batch_index, const QmomState& state)
        {
            const MomentSource local_source = source
                                                  ? MomentSource{[&source, batch_index](const double callback_time, const QmomSourceContext& context)
                                                                 { return source(batch_index, callback_time, context); }}
                                                  : MomentSource{};
            const MomentBoundaryContribution local_boundary =
                boundary ? MomentBoundaryContribution{[&boundary, batch_index](const double callback_time, const std::vector<double>& moments)
                                                      { return boundary(batch_index, callback_time, moments); }}
                         : MomentBoundaryContribution{};
            return advance_qmom_moment_source(state, maximum_node_count, time, time_step, local_source, local_boundary, options);
        });
}

std::vector<CqmomPhaseSpaceEvolutionResult> advance_cqmom_moment_source_batch_with_index(const std::vector<CqmomState>& states, const double time,
                                                                                         const double time_step, const BatchCqmomMomentSource& source,
                                                                                         const BatchCqmomMomentBoundaryContribution& boundary,
                                                                                         const CqmomPhaseSpaceEvolutionOptions& options)
{
    return apply_batch(
        states,
        [&](const std::size_t batch_index, const CqmomState& state)
        {
            const CqmomMomentSource local_source =
                source ? CqmomMomentSource{[&source, batch_index](const double callback_time, const CqmomSourceContext& context)
                                           { return source(batch_index, callback_time, context); }}
                       : CqmomMomentSource{};
            const CqmomMomentBoundaryContribution local_boundary =
                boundary ? CqmomMomentBoundaryContribution{[&boundary, batch_index](const double callback_time, const CqmomSourceContext& context)
                                                           { return boundary(batch_index, callback_time, context); }}
                         : CqmomMomentBoundaryContribution{};
            return advance_cqmom_moment_source(state, time, time_step, local_source, local_boundary, options);
        });
}

std::vector<QmomPhaseSpaceEvolutionResult> advance_qmom_phase_space_batch_with_index(const std::vector<QmomState>& states,
                                                                                     const std::size_t maximum_node_count, const double time,
                                                                                     const double time_step, const BatchPhaseSpaceVelocity& velocity,
                                                                                     const BatchMomentSource& source,
                                                                                     const QmomPhaseSpaceEvolutionOptions& options)
{
    return apply_batch(
        states,
        [&](const std::size_t batch_index, const QmomState& state)
        {
            const PhaseSpaceVelocity local_velocity =
                velocity ? PhaseSpaceVelocity{[&velocity, batch_index](const double callback_time, const std::vector<double>& phase_space_location)
                                              { return velocity(batch_index, callback_time, phase_space_location); }}
                         : PhaseSpaceVelocity{};
            const MomentSource local_source = source
                                                  ? MomentSource{[&source, batch_index](const double callback_time, const QmomSourceContext& context)
                                                                 { return source(batch_index, callback_time, context); }}
                                                  : MomentSource{};
            return advance_qmom_phase_space(state, maximum_node_count, time, time_step, local_velocity, local_source, options);
        });
}

std::vector<CqmomPhaseSpaceEvolutionResult> advance_cqmom_phase_space_batch_with_index(const std::vector<CqmomState>& states, const double time,
                                                                                       const double time_step,
                                                                                       const BatchPhaseSpaceVelocity& velocity,
                                                                                       const BatchCqmomMomentSource& source,
                                                                                       const CqmomPhaseSpaceEvolutionOptions& options)
{
    return apply_batch(
        states,
        [&](const std::size_t batch_index, const CqmomState& state)
        {
            const PhaseSpaceVelocity local_velocity =
                velocity ? PhaseSpaceVelocity{[&velocity, batch_index](const double callback_time, const std::vector<double>& phase_space_location)
                                              { return velocity(batch_index, callback_time, phase_space_location); }}
                         : PhaseSpaceVelocity{};
            const CqmomMomentSource local_source =
                source ? CqmomMomentSource{[&source, batch_index](const double callback_time, const CqmomSourceContext& context)
                                           { return source(batch_index, callback_time, context); }}
                       : CqmomMomentSource{};
            return advance_cqmom_phase_space(state, time, time_step, local_velocity, local_source, options);
        });
}

std::vector<QmomPhaseSpaceEvolutionResult>
advance_qmom_phase_space_batch_with_index(const std::vector<QmomState>& states, const std::size_t maximum_node_count, const double time,
                                          const double time_step, const BatchPhaseSpaceVelocity& velocity, const BatchMomentSource& source,
                                          const BatchMomentBoundaryContribution& boundary, const QmomPhaseSpaceEvolutionOptions& options)
{
    return apply_batch(
        states,
        [&](const std::size_t batch_index, const QmomState& state)
        {
            const PhaseSpaceVelocity local_velocity =
                velocity ? PhaseSpaceVelocity{[&velocity, batch_index](const double callback_time, const std::vector<double>& phase_space_location)
                                              { return velocity(batch_index, callback_time, phase_space_location); }}
                         : PhaseSpaceVelocity{};
            const MomentSource local_source = source
                                                  ? MomentSource{[&source, batch_index](const double callback_time, const QmomSourceContext& context)
                                                                 { return source(batch_index, callback_time, context); }}
                                                  : MomentSource{};
            const MomentBoundaryContribution local_boundary =
                boundary ? MomentBoundaryContribution{[&boundary, batch_index](const double callback_time, const std::vector<double>& moments)
                                                      { return boundary(batch_index, callback_time, moments); }}
                         : MomentBoundaryContribution{};
            return advance_qmom_phase_space(state, maximum_node_count, time, time_step, local_velocity, local_source, local_boundary, options);
        });
}

std::vector<CqmomPhaseSpaceEvolutionResult>
advance_cqmom_phase_space_batch_with_index(const std::vector<CqmomState>& states, const double time, const double time_step,
                                           const BatchPhaseSpaceVelocity& velocity, const BatchCqmomMomentSource& source,
                                           const BatchCqmomMomentBoundaryContribution& boundary, const CqmomPhaseSpaceEvolutionOptions& options)
{
    return apply_batch(
        states,
        [&](const std::size_t batch_index, const CqmomState& state)
        {
            const PhaseSpaceVelocity local_velocity =
                velocity ? PhaseSpaceVelocity{[&velocity, batch_index](const double callback_time, const std::vector<double>& phase_space_location)
                                              { return velocity(batch_index, callback_time, phase_space_location); }}
                         : PhaseSpaceVelocity{};
            const CqmomMomentSource local_source =
                source ? CqmomMomentSource{[&source, batch_index](const double callback_time, const CqmomSourceContext& context)
                                           { return source(batch_index, callback_time, context); }}
                       : CqmomMomentSource{};
            const CqmomMomentBoundaryContribution local_boundary =
                boundary ? CqmomMomentBoundaryContribution{[&boundary, batch_index](const double callback_time, const CqmomSourceContext& context)
                                                           { return boundary(batch_index, callback_time, context); }}
                         : CqmomMomentBoundaryContribution{};
            return advance_cqmom_phase_space(state, time, time_step, local_velocity, local_source, local_boundary, options);
        });
}

} // namespace qmomkit
