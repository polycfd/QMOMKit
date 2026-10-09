#pragma once

#include "qmomkit/runge_kutta.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace qmomkit::internal
{

constexpr std::size_t maximum_embedded_rk_stages = 7;

struct EmbeddedRkTableau
{
    const char* name = "";
    std::size_t stage_count = 0;
    std::size_t higher_order = 0;
    std::size_t embedded_order = 0;
    std::array<double, maximum_embedded_rk_stages> c{};
    std::array<std::array<double, maximum_embedded_rk_stages>, maximum_embedded_rk_stages> a{};
    std::array<double, maximum_embedded_rk_stages> b{};
    std::array<double, maximum_embedded_rk_stages> b_embedded{};
    bool fsal = false;
};

/**
 * Buffers and cached slope carried across the attempts of one integration.
 *
 * The controller owns a single workspace for the whole integration and hands it
 * to every step attempt, so the stage vectors are allocated on the first attempt
 * and reused afterwards instead of being allocated per attempt. A workspace is
 * an ordinary local of the controller rather than shared state, so concurrent
 * integrations each work on their own.
 *
 * After an attempt the workspace also holds that attempt's outcome:
 * `higher_state` and `embedded_state` are the two candidate solutions and
 * `error` their componentwise difference.
 */
struct EmbeddedRkWorkspace
{
    /**
     * Stage derivatives, of which the first `stage_count` entries are in use.
     * Slot zero is the slope at the starting state; for a first-same-as-last
     * tableau the last slot in use is the slope at the accepted state.
     */
    std::array<std::vector<double>, maximum_embedded_rk_stages> derivatives;
    /** State that each stage after the first evaluates the right-hand side at. */
    std::vector<double> stage_state;
    std::vector<double> higher_state;
    std::vector<double> embedded_state;
    std::vector<double> error;
    /**
     * Whether `derivatives[0]` holds the slope at the starting state and time of
     * the next attempt. Keeping the cache flag beside the buffer it describes is
     * what lets the controller carry a slope from one attempt to the next
     * without moving it out of the workspace and reallocating it.
     */
    bool has_first_derivative = false;
};

[[nodiscard]] const EmbeddedRkTableau& embedded_rk_tableau(EmbeddedRkMethod method);

/**
 * Take one attempt of an embedded Runge-Kutta step, leaving the outcome in
 * `workspace`.
 *
 * A workspace whose `has_first_derivative` is set supplies the slope at `state`
 * and `t` instead of it being evaluated; validity of that cache is the caller's
 * invariant, exactly as it was when the slope was passed in explicitly. On
 * return the flag is set, because `derivatives[0]` is then the slope at `state`
 * and `t` whichever way it was obtained, so a retry of the same step with a
 * smaller `step` reuses it. A controller that accepts the attempt must instead
 * either install the outgoing slope of a first-same-as-last tableau or clear the
 * flag, since the cache does not belong to the state it just moved to.
 */
void take_embedded_rk_step(const EmbeddedRkTableau& tableau, const OdeRhs& rhs, double t, const std::vector<double>& state, double step,
                           std::size_t& function_evaluations, EmbeddedRkWorkspace& workspace);

[[nodiscard]] double embedded_rk_controller_exponent(const EmbeddedRkTableau& tableau);

// Why `options` cannot drive an embedded Runge-Kutta integration, or nullptr
// when they can. Phase-space evolution uses it to check node-advection options
// before any stage runs.
[[nodiscard]] const char* adaptive_rk_options_error(const AdaptiveRkOptions& options);

} // namespace qmomkit::internal
