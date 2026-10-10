#include "qmomkit/runge_kutta.hpp"

#include "internal/embedded_runge_kutta.hpp"
#include "internal/numeric_utils.hpp"
#include "internal/ode_rhs_validation.hpp"
#include "internal/performance.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace qmomkit
{
namespace
{

void require_same_size(const std::vector<double>& lhs, const std::vector<double>& rhs, const char* message)
{
    if (lhs.size() != rhs.size())
    {
        throw std::invalid_argument(message);
    }
}

[[nodiscard]] double euclidean_norm(const std::vector<double>& values)
{
    double sum = 0.0;
    for (const double value : values)
    {
        sum += value * value;
    }
    return std::sqrt(sum);
}

[[nodiscard]] std::vector<double> add_scaled(const std::vector<double>& state, const double scale, const std::vector<double>& derivative)
{
    require_same_size(state, derivative, "Vector sizes do not match.");
    std::vector<double> result(state.size(), 0.0);
    for (std::size_t i = 0; i < state.size(); ++i)
    {
        result[i] = state[i] + scale * derivative[i];
    }
    return result;
}

[[nodiscard]] double sign(const double value)
{
    return value >= 0.0 ? 1.0 : -1.0;
}

[[nodiscard]] double validate_or_default_max_step(const AdaptiveRkOptions& options, const double interval)
{
    if (options.max_step > 0.0)
    {
        return options.max_step;
    }
    return std::abs(interval);
}

// The automatic minimum is ten spacings of the double at t. An explicit
// min_step is never smaller: a shorter step can round to no step at all, or
// back to the step just rejected, and the controller would never finish.
[[nodiscard]] double minimum_step(const double t, const double direction, const AdaptiveRkOptions& options)
{
    const double toward_infinity = direction > 0.0 ? std::numeric_limits<double>::infinity() : -std::numeric_limits<double>::infinity();
    return std::max(options.min_step, 10.0 * std::abs(std::nextafter(t, toward_infinity) - t));
}

struct SspRk3Attempt
{
    SspRk3Result result;
    bool accepted = true;
};

[[nodiscard]] double error_norm(const std::vector<double>& error, const std::vector<double>& y, const std::vector<double>& y_new,
                                const AdaptiveRkOptions& options)
{
    if (error.empty())
    {
        return 0.0;
    }

    // Root-mean-square of the tolerance-scaled error, accumulated in index
    // order so the summation matches a scaled buffer reduced term by term.
    double sum_of_squares = 0.0;
    for (std::size_t i = 0; i < error.size(); ++i)
    {
        const double scale = options.absolute_tolerance + std::max(std::abs(y[i]), std::abs(y_new[i])) * options.relative_tolerance;
        const double scaled = error[i] / scale;
        sum_of_squares += scaled * scaled;
    }
    return std::sqrt(sum_of_squares) / std::sqrt(static_cast<double>(error.size()));
}

[[nodiscard]] AdaptiveRkResult invalid_embedded_rk_result(std::vector<double> state, const double t, std::string message)
{
    return AdaptiveRkResult{std::move(state), t, AdaptiveRkStatus::invalid_input, std::move(message), 0, 0, 0, 0.0};
}

[[nodiscard]] bool stage_accepted(const StageAcceptancePredicate& accept, const SspRk3Stage stage, const std::vector<double>& state)
{
    return !accept || accept(stage, state);
}

[[nodiscard]] SspRk3Attempt try_ssp_rk3_step(const std::vector<double>& state, const double t, const double time_step, const OdeRhs& rhs,
                                             const SspRk3Options& options)
{
    // Perform exactly one fixed-timestep SSP-RK3 attempt. The caller owns any
    // retry policy and discards rejected attempts by restarting from state.
    // One derivative buffer serves all three stages: each stage's derivative is
    // consumed into the next stage state before the following stage overwrites
    // it, so the attempt evaluates three stages against one allocation.
    std::vector<double> derivative;

    internal::performance::record_source_rk_stage();
    internal::rhs_checked(rhs, t, state, derivative);
    const std::vector<double> stage_1 = add_scaled(state, time_step, derivative);
    if (!stage_accepted(options.stage_accept, SspRk3Stage::Stage1, stage_1))
    {
        return SspRk3Attempt{SspRk3Result{state, std::nullopt, SspRk3Status::success, 0, time_step}, false};
    }

    internal::performance::record_source_rk_stage();
    internal::rhs_checked(rhs, t + time_step, stage_1, derivative);
    std::vector<double> test_state(state.size(), 0.0);
    std::vector<double> stage_2(state.size(), 0.0);
    for (std::size_t i = 0; i < state.size(); ++i)
    {
        test_state[i] = 0.5 * (state[i] + (stage_1[i] + time_step * derivative[i]));
        stage_2[i] = 0.75 * state[i] + 0.25 * (stage_1[i] + time_step * derivative[i]);
    }
    if (!stage_accepted(options.stage_accept, SspRk3Stage::Stage2, stage_2))
    {
        return SspRk3Attempt{SspRk3Result{state, std::nullopt, SspRk3Status::success, 0, time_step}, false};
    }

    internal::performance::record_source_rk_stage();
    internal::rhs_checked(rhs, t + 0.5 * time_step, stage_2, derivative);
    std::vector<double> new_state(state.size(), 0.0);
    for (std::size_t i = 0; i < state.size(); ++i)
    {
        new_state[i] = (state[i] + 2.0 * (stage_2[i] + time_step * derivative[i])) / 3.0;
    }
    if (!internal::all_finite(new_state))
    {
        throw std::invalid_argument("SSP-RK3 produced a non-finite candidate state.");
    }

    std::optional<double> error_estimate;
    if (options.estimate_error)
    {
        std::vector<double> difference(state.size(), 0.0);
        for (std::size_t i = 0; i < state.size(); ++i)
        {
            difference[i] = new_state[i] - test_state[i];
        }
        const double denominator = euclidean_norm(new_state);
        error_estimate = denominator == 0.0 ? euclidean_norm(difference) : euclidean_norm(difference) / denominator;
    }

    const bool accepted = stage_accepted(options.stage_accept, SspRk3Stage::Final, new_state) && (!options.accept || options.accept(new_state));
    return SspRk3Attempt{SspRk3Result{std::move(new_state), error_estimate, SspRk3Status::success, 0, time_step}, accepted};
}

} // namespace

SspRk3Result ssp_rk3_step(const std::vector<double>& state, const double t, const double time_step, const OdeRhs& rhs, const SspRk3Options& options)
{
    const auto invalid_input = [&state](std::string message)
    { return SspRk3Result{state, std::nullopt, SspRk3Status::invalid_input, 0, 0.0, std::move(message)}; };
    if (state.empty())
    {
        return invalid_input("SSP-RK3 state must not be empty.");
    }
    if (!internal::all_finite(state) || !std::isfinite(t) || !std::isfinite(time_step))
    {
        return invalid_input("SSP-RK3 input contains non-finite values.");
    }
    // A step already shorter than min_time_step is attempted, just not halved,
    // so that a caller retrying with a shorter step is not refused for it.
    if (time_step <= 0.0 || !std::isfinite(options.min_time_step) || options.min_time_step < 0.0)
    {
        return invalid_input("SSP-RK3 timestep controls must be positive.");
    }

    double candidate_dt = time_step;
    std::size_t retries = 0;
    while (true)
    {
        SspRk3Attempt attempt = try_ssp_rk3_step(state, t, candidate_dt, rhs, options);
        SspRk3Result candidate = std::move(attempt.result);
        candidate.retries = retries;

        if (attempt.accepted)
        {
            return candidate;
        }

        if (retries >= options.max_retries)
        {
            return SspRk3Result{state, std::nullopt, SspRk3Status::rejected_by_predicate, retries, candidate_dt};
        }

        const double next_dt = 0.5 * candidate_dt;
        if (!std::isfinite(next_dt) || next_dt <= 0.0 || next_dt >= candidate_dt || (options.min_time_step > 0.0 && next_dt < options.min_time_step))
        {
            return SspRk3Result{state, std::nullopt, SspRk3Status::minimum_timestep_reached, retries, candidate_dt};
        }

        candidate_dt = next_dt;
        internal::performance::record_source_retry();
        ++retries;
    }
}

const char* internal::adaptive_rk_options_error(const AdaptiveRkOptions& options)
{
    if (!std::isfinite(options.relative_tolerance) || !std::isfinite(options.absolute_tolerance) || options.relative_tolerance <= 0.0 ||
        options.absolute_tolerance <= 0.0)
    {
        return "Embedded Runge-Kutta tolerances must be finite and positive.";
    }
    if (options.max_steps == 0)
    {
        return "Embedded Runge-Kutta max_steps must be positive.";
    }
    if ((!std::isfinite(options.initial_step) || options.initial_step <= 0.0) || !std::isfinite(options.max_step) || options.max_step < 0.0 ||
        !std::isfinite(options.min_step) || options.min_step < 0.0)
    {
        return "Invalid embedded Runge-Kutta step-size option.";
    }
    // Only two explicit limits can contradict each other. A default max_step
    // is the interval itself, and an interval shorter than min_step is taken
    // as one step, so that a caller retrying with a shorter interval is not
    // refused for it.
    if (options.max_step > 0.0 && options.min_step > options.max_step)
    {
        return "Embedded Runge-Kutta min_step must not exceed max_step.";
    }
    return nullptr;
}

AdaptiveRkResult integrate_embedded_rk(const std::vector<double>& initial_state, const double t0, const double t_bound, const OdeRhs& rhs,
                                       const EmbeddedRkMethod method, const AdaptiveRkOptions& options)
{
    const internal::EmbeddedRkTableau* tableau = nullptr;
    try
    {
        tableau = &internal::embedded_rk_tableau(method);
    }
    catch (const std::exception& error)
    {
        return invalid_embedded_rk_result(initial_state, t0, error.what());
    }

    if (initial_state.empty())
    {
        return invalid_embedded_rk_result(initial_state, t0, "Embedded Runge-Kutta state must not be empty.");
    }
    if (!internal::all_finite(initial_state) || !std::isfinite(t0) || !std::isfinite(t_bound))
    {
        return invalid_embedded_rk_result(initial_state, t0, "Invalid embedded Runge-Kutta input.");
    }
    // Validated before the empty-interval return below, so invalid options are
    // reported whether or not the interval is empty.
    if (const char* error = internal::adaptive_rk_options_error(options); error != nullptr)
    {
        return invalid_embedded_rk_result(initial_state, t0, error);
    }
    if (t0 == t_bound)
    {
        return AdaptiveRkResult{initial_state, t0, AdaptiveRkStatus::success, "", 0, 0, 0, 0.0};
    }

    const double direction = sign(t_bound - t0);
    const double max_step = validate_or_default_max_step(options, t_bound - t0);
    if (max_step <= 0.0 || !std::isfinite(max_step))
    {
        return invalid_embedded_rk_result(initial_state, t0, "Invalid embedded Runge-Kutta maximum step.");
    }

    std::vector<double> y = initial_state;
    double t = t0;
    std::size_t function_evaluations = 0;
    double h_abs = std::min(std::abs(options.initial_step), max_step);
    // One workspace serves every attempt of this integration, so the stage
    // vectors are allocated on the first attempt rather than on each one. It is
    // a local and not shared state, so two threads integrating at the same time
    // never meet in it.
    internal::EmbeddedRkWorkspace workspace;

    AdaptiveRkResult result;
    result.state = y;
    result.final_time = t;
    result.status = AdaptiveRkStatus::success;
    result.function_evaluations = function_evaluations;

    constexpr double safety = 0.9;
    constexpr double min_factor = 0.2;
    constexpr double max_factor = 10.0;
    const double error_exponent = internal::embedded_rk_controller_exponent(*tableau);

    while (direction * (t_bound - t) > 0.0)
    {
        if (result.accepted_steps >= options.max_steps)
        {
            result.state = y;
            result.final_time = t;
            result.status = AdaptiveRkStatus::max_steps_exceeded;
            result.message = std::string(tableau->name) + " exceeded max_steps.";
            result.function_evaluations = function_evaluations;
            return result;
        }

        const double min_step = minimum_step(t, direction, options);
        h_abs = std::min(std::max(h_abs, min_step), max_step);

        bool step_accepted = false;
        bool step_rejected = false;
        while (!step_accepted)
        {
            const double remaining_interval = std::abs(t_bound - t);
            const double effective_lower_bound = std::min(min_step, remaining_interval);
            if (h_abs < effective_lower_bound)
            {
                h_abs = effective_lower_bound;
            }
            const double requested_h_abs = h_abs;

            double h = direction * h_abs;
            double t_new = t + h;
            if (direction * (t_new - t_bound) > 0.0)
            {
                t_new = t_bound;
            }
            h = t_new - t;
            h_abs = std::abs(h);

            internal::take_embedded_rk_step(*tableau, rhs, t, y, h, function_evaluations, workspace);
            const bool step_is_finite = internal::all_finite(workspace.higher_state) && internal::all_finite(workspace.embedded_state) &&
                                        internal::all_finite(workspace.error);
            const double norm =
                step_is_finite ? error_norm(workspace.error, y, workspace.higher_state, options) : std::numeric_limits<double>::infinity();
            const bool error_is_acceptable = std::isfinite(norm) && norm < 1.0;
            const bool step_is_at_effective_lower_bound = requested_h_abs <= effective_lower_bound || h_abs <= effective_lower_bound;
            const bool lower_bound_acceptance = step_is_finite && std::isfinite(norm) && step_is_at_effective_lower_bound;

            if (error_is_acceptable || lower_bound_acceptance)
            {
                double factor = max_factor;
                if (norm > 0.0)
                {
                    factor = std::min(max_factor, safety * std::pow(norm, error_exponent));
                }
                if (step_rejected)
                {
                    factor = std::min(1.0, factor);
                }

                result.last_step = h;
                ++result.accepted_steps;
                internal::performance::record_embedded_rk_step(true);
                t = t_new;
                // Swapping rather than assigning hands the retired state vector
                // back to the workspace, which overwrites it on the next attempt.
                y.swap(workspace.higher_state);
                // The cached slope belongs to the state the step started from, so
                // accepting invalidates it. A first-same-as-last tableau evaluates
                // its last stage at the accepted state, which is the slope the next
                // step starts from, so that stage becomes the new cache.
                if (tableau->fsal)
                {
                    workspace.derivatives[0].swap(workspace.derivatives[tableau->stage_count - 1]);
                }
                else
                {
                    workspace.has_first_derivative = false;
                }
                h_abs *= factor;
                step_accepted = true;
            }
            else
            {
                double factor = min_factor;
                if (std::isfinite(norm) && norm > 0.0)
                {
                    factor = std::max(min_factor, safety * std::pow(norm, error_exponent));
                }
                h_abs *= factor;
                step_rejected = true;
                // A rejected attempt leaves t and y where they were, so the slope
                // the attempt evaluated there is still the one the retry needs and
                // the workspace keeps it.
                ++result.rejected_steps;
                internal::performance::record_embedded_rk_step(false);
                if (step_is_at_effective_lower_bound)
                {
                    // Reaching this branch at the lower bound means the attempt
                    // was not merely inaccurate: a finite candidate with a finite
                    // norm is accepted there regardless of the error. What is
                    // left is a non-finite state or error estimate, so the
                    // message says that rather than blaming the step size the
                    // status name refers to.
                    result.state = y;
                    result.final_time = t;
                    result.status = AdaptiveRkStatus::step_size_underflow;
                    result.message = std::string(tableau->name) + " produced a non-finite state or error estimate at the minimum step.";
                    result.function_evaluations = function_evaluations;
                    return result;
                }
            }
        }
    }

    result.state = std::move(y);
    result.final_time = t;
    result.status = AdaptiveRkStatus::success;
    result.function_evaluations = function_evaluations;
    return result;
}

} // namespace qmomkit
