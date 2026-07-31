#pragma once

#include <cstddef>
#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace qmomkit
{

using OdeRhs = std::function<std::vector<double>(double, const std::vector<double>&)>;
using StateAcceptancePredicate = std::function<bool(const std::vector<double>&)>;

enum class SspRk3Stage
{
    Stage1,
    Stage2,
    Final
};

using StageAcceptancePredicate = std::function<bool(SspRk3Stage, const std::vector<double>&)>;

enum class SspRk3Status
{
    success,
    rejected_by_predicate,
    minimum_timestep_reached
};

/**
 * Controls for optional SSP-RK3 candidate acceptance.
 *
 * If accept is supplied, each candidate step is tested after the full SSP-RK3
 * update. Rejected candidates are discarded, the timestep is halved, and the
 * step is retried from the original state. If stage_accept is supplied, it is
 * called for stage-1, stage-2, and final candidate states during each attempt.
 * The predicates are entirely generic: they must not encode CQMOM-specific
 * behavior inside the solver.
 */
struct SspRk3Options
{
    bool estimate_error = false;
    StateAcceptancePredicate accept;
    StageAcceptancePredicate stage_accept;
    std::size_t max_retries = 0;
    double min_time_step = 0.0;
};

/**
 * Result of one explicit SSP-RK3 step.
 *
 * The optional error estimate is ||y_new - y_test|| / ||y_new||, where y_test
 * is the embedded second-order state used only for timestep estimation.
 */
struct SspRk3Result
{
    std::vector<double> state;
    std::optional<double> error_estimate;
    SspRk3Status status = SspRk3Status::success;
    std::size_t retries = 0;
    double time_step_used = 0.0;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == SspRk3Status::success;
    }
};

/**
 * Advance a vector-valued state by one deterministic SSP-RK3 step.
 *
 * The RHS callback must return a vector with the same size as the input state.
 * This routine is generic and does not assume any CQMOM, particle, or physics
 * semantics.
 */
[[nodiscard]] SspRk3Result ssp_rk3_step(const std::vector<double>& state, double t, double time_step, const OdeRhs& rhs, bool estimate_error = false);

/**
 * Advance a vector-valued state by one SSP-RK3 step with optional generic
 * acceptance/rejection.
 */
[[nodiscard]] SspRk3Result ssp_rk3_step(const std::vector<double>& state, double t, double time_step, const OdeRhs& rhs,
                                        const SspRk3Options& options);

/**
 * Adapt a timestep using the SSP-RK3 error controller.
 */
[[nodiscard]] double adapt_ssp_rk3_timestep(double current_dt, double timestep_error, double error_tolerance, double min_dt, double max_dt);

enum class AdaptiveRk23Status
{
    success,
    invalid_input,
    step_size_underflow,
    max_steps_exceeded
};

/**
 * Controls for the generic adaptive RK23 integrator.
 *
 * The coefficients and error control are the Bogacki-Shampine 3(2) pair used by
 * SciPy's solve_ivp(method="RK23"). Scalar tolerances are applied componentwise
 * with scale_i = atol + max(|y_i|, |y_new_i|) * rtol.
 */
struct AdaptiveRk23Options
{
    double relative_tolerance = 1.0e-6;
    double absolute_tolerance = 1.0e-9;
    std::optional<double> initial_step;
    double max_step = 0.0;
    double min_step = 0.0;
    std::size_t max_steps = 1000000000;
};

struct AdaptiveRk23Result
{
    std::vector<double> state;
    double t = 0.0;
    AdaptiveRk23Status status = AdaptiveRk23Status::success;
    std::string message;
    std::size_t accepted_steps = 0;
    std::size_t rejected_steps = 0;
    std::size_t function_evaluations = 0;
    double last_step = 0.0;
};

/**
 * Integrate y' = rhs(t, y) from t0 to t_bound with adaptive RK23.
 */
[[nodiscard]] AdaptiveRk23Result integrate_rk23(const std::vector<double>& initial_state, double t0, double t_bound, const OdeRhs& rhs,
                                                const AdaptiveRk23Options& options = {});

} // namespace qmomkit
