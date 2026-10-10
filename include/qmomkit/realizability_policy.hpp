#pragma once

#include "qmomkit/realizability.hpp"
#include "qmomkit/runge_kutta.hpp"

#include <functional>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace qmomkit
{

/**
 * Mode for optional runtime realizability checks during SSP-RK3 integration.
 *
 * `None` installs no acceptance predicate. `FinalOnly` checks only complete
 * candidate states, preserving the original SSP-RK3 acceptance behavior.
 * `EveryStage` checks stage-1, stage-2, and final candidate states.
 */
enum class RealizabilityCheckMode
{
    None,
    EveryStage,
    FinalOnly
};

/**
 * Opt-in runtime policy for using one-dimensional realizability as a candidate
 * state acceptance predicate.
 *
 * The policy is separate from the generic Runge-Kutta solvers. Callers own any
 * conversion from their numerical state representation to `MomentSet`.
 *
 * A policy accepts strictly realizable moments, boundary moments when
 * `accept_boundary` is set, and moments the check cannot classify
 * (`numerically_indeterminate`) when `accept_indeterminate` is set. The latter
 * arises, for example, for distributions that are narrow compared with their
 * distance from the origin, at high orders, where the moments carry too few
 * digits to decide, or with a coarse tolerance. Such moments may lie just
 * outside the realizable set within the precision of the data; accepting them
 * leaves them to the reconstruction, which reduces the node count in adaptive
 * mode and may fail in fixed mode.
 */
struct RealizabilityPolicy
{
    RealizabilityCheckMode check_mode = RealizabilityCheckMode::None;
    MomentSupport support = MomentSupport::real_line;
    std::optional<BoundedInterval> interval;
    RealizabilityOptions options;
    bool accept_boundary = true;
    bool accept_indeterminate = false;
};

using MomentStateConverter = std::function<MomentSet(const std::vector<double>&)>;

/**
 * Whether `policy` accepts moments that the check classified as `status`. The
 * check mode decides when a check runs, not what it accepts.
 */
[[nodiscard]] inline bool accepts_realizability_status(const RealizabilityPolicy& policy, const RealizabilityStatus status) noexcept
{
    switch (status)
    {
    case RealizabilityStatus::strictly_realizable:
        return true;
    case RealizabilityStatus::boundary:
        return policy.accept_boundary;
    case RealizabilityStatus::numerically_indeterminate:
        return policy.accept_indeterminate;
    case RealizabilityStatus::not_realizable:
    case RealizabilityStatus::invalid_input:
        return false;
    }
    return false;
}

/**
 * Evaluate a policy against a candidate generic state.
 *
 * Disabled policies accept without invoking the converter. Enabled policies
 * convert exactly once, run `check_realizability(...)`, and accept the result
 * as `accepts_realizability_status(...)` does. Converter exceptions are not
 * caught; they follow the acceptance-predicate exception behavior of the
 * calling integrator or wrapper.
 */
[[nodiscard]] inline bool accepts_realizability_policy(const RealizabilityPolicy& policy, const std::vector<double>& state,
                                                       const MomentStateConverter& converter)
{
    if (policy.check_mode == RealizabilityCheckMode::None)
    {
        return true;
    }
    if (!converter)
    {
        throw std::invalid_argument("Realizability policy requires a state-to-moment converter.");
    }

    const MomentSet moments = converter(state);
    const RealizabilityResult result = (policy.support == MomentSupport::bounded_interval && policy.interval.has_value())
                                           ? check_realizability(moments, *policy.interval, policy.options)
                                           : check_realizability(moments, policy.support, policy.options);
    return accepts_realizability_status(policy, result.status);
}

/**
 * Adapt a realizability policy to the existing SSP-RK3 acceptance callback.
 *
 * Only FinalOnly policies produce a predicate for the final-candidate
 * acceptance callback. None and EveryStage policies return an empty predicate; use
 * make_realizability_stage_acceptance_predicate for stage-aware integration.
 * Converter exceptions propagate through the returned predicate.
 */
[[nodiscard]] inline StateAcceptancePredicate make_realizability_acceptance_predicate(const RealizabilityPolicy& policy,
                                                                                      MomentStateConverter converter)
{
    if (policy.check_mode != RealizabilityCheckMode::FinalOnly)
    {
        return {};
    }

    return [policy, converter = std::move(converter)](const std::vector<double>& state)
    { return accepts_realizability_policy(policy, state, converter); };
}

/**
 * Adapt a realizability policy to the SSP-RK3 stage-aware acceptance callback.
 *
 * None returns an empty predicate. FinalOnly checks only SspRk3Stage::Final.
 * EveryStage checks every state supplied by the SSP-RK3 integrator.
 */
[[nodiscard]] inline StageAcceptancePredicate make_realizability_stage_acceptance_predicate(const RealizabilityPolicy& policy,
                                                                                            MomentStateConverter converter)
{
    if (policy.check_mode == RealizabilityCheckMode::None)
    {
        return {};
    }

    return [policy, converter = std::move(converter)](const SspRk3Stage stage, const std::vector<double>& state)
    {
        if (policy.check_mode == RealizabilityCheckMode::FinalOnly && stage != SspRk3Stage::Final)
        {
            return true;
        }
        return accepts_realizability_policy(policy, state, converter);
    };
}

} // namespace qmomkit
