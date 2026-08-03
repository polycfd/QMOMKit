#include "qmomkit/runge_kutta.hpp"

#include "internal/embedded_runge_kutta.hpp"
#include "internal/numeric_utils.hpp"
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

[[nodiscard]] double rms_norm(const std::vector<double>& values)
{
    if (values.empty())
    {
        return 0.0;
    }
    return euclidean_norm(values) / std::sqrt(static_cast<double>(values.size()));
}

[[nodiscard]] std::vector<double> rhs_checked(const OdeRhs& rhs, const double t, const std::vector<double>& state)
{
    std::vector<double> derivative = rhs(t, state);
    require_same_size(derivative, state, "RHS returned a vector with the wrong size.");
    if (!internal::all_finite(derivative))
    {
        throw std::invalid_argument("RHS returned a non-finite value.");
    }
    return derivative;
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

[[nodiscard]] double minimum_step(const double t, const double direction, const AdaptiveRkOptions& options)
{
    if (options.min_step > 0.0)
    {
        return options.min_step;
    }
    const double toward_infinity = direction > 0.0 ? std::numeric_limits<double>::infinity() : -std::numeric_limits<double>::infinity();
    return 10.0 * std::abs(std::nextafter(t, toward_infinity) - t);
}

struct SspRk3Attempt
{
    SspRk3Result result;
    bool accepted = true;
};

[[nodiscard]] double error_norm(const std::vector<double>& error, const std::vector<double>& y, const std::vector<double>& y_new,
                                const AdaptiveRkOptions& options)
{
    std::vector<double> scaled(error.size(), 0.0);
    for (std::size_t i = 0; i < error.size(); ++i)
    {
        const double scale = options.absolute_tolerance + std::max(std::abs(y[i]), std::abs(y_new[i])) * options.relative_tolerance;
        scaled[i] = error[i] / scale;
    }
    return rms_norm(scaled);
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
    internal::performance::record_source_rk_stage();
    const std::vector<double> k1 = rhs_checked(rhs, t, state);
    const std::vector<double> stage_1 = add_scaled(state, time_step, k1);
    if (!stage_accepted(options.stage_accept, SspRk3Stage::Stage1, stage_1))
    {
        return SspRk3Attempt{SspRk3Result{state, std::nullopt, SspRk3Status::success, 0, time_step}, false};
    }

    internal::performance::record_source_rk_stage();
    const std::vector<double> k2 = rhs_checked(rhs, t + time_step, stage_1);
    std::vector<double> test_state(state.size(), 0.0);
    std::vector<double> stage_2(state.size(), 0.0);
    for (std::size_t i = 0; i < state.size(); ++i)
    {
        test_state[i] = 0.5 * (state[i] + (stage_1[i] + time_step * k2[i]));
        stage_2[i] = 0.75 * state[i] + 0.25 * (stage_1[i] + time_step * k2[i]);
    }
    if (!stage_accepted(options.stage_accept, SspRk3Stage::Stage2, stage_2))
    {
        return SspRk3Attempt{SspRk3Result{state, std::nullopt, SspRk3Status::success, 0, time_step}, false};
    }

    internal::performance::record_source_rk_stage();
    const std::vector<double> k3 = rhs_checked(rhs, t + 0.5 * time_step, stage_2);
    std::vector<double> new_state(state.size(), 0.0);
    for (std::size_t i = 0; i < state.size(); ++i)
    {
        new_state[i] = (state[i] + 2.0 * (stage_2[i] + time_step * k3[i])) / 3.0;
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

SspRk3Result ssp_rk3_step(const std::vector<double>& state, const double t, const double time_step, const OdeRhs& rhs, const bool estimate_error)
{
    SspRk3Options options;
    options.estimate_error = estimate_error;
    return ssp_rk3_step(state, t, time_step, rhs, options);
}

SspRk3Result ssp_rk3_step(const std::vector<double>& state, const double t, const double time_step, const OdeRhs& rhs, const SspRk3Options& options)
{
    if (state.empty())
    {
        throw std::invalid_argument("SSP-RK3 state must not be empty.");
    }
    if (!internal::all_finite(state) || !std::isfinite(t) || !std::isfinite(time_step))
    {
        throw std::invalid_argument("SSP-RK3 input contains non-finite values.");
    }
    if (time_step <= 0.0 || options.min_time_step < 0.0 || (options.min_time_step > 0.0 && time_step < options.min_time_step))
    {
        throw std::invalid_argument("SSP-RK3 timestep controls must be positive.");
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
        if (options.min_time_step > 0.0 && next_dt < options.min_time_step)
        {
            return SspRk3Result{state, std::nullopt, SspRk3Status::minimum_timestep_reached, retries, candidate_dt};
        }

        candidate_dt = next_dt;
        internal::performance::record_source_retry();
        ++retries;
    }
}

double adapt_ssp_rk3_timestep(const double current_dt, const double timestep_error, const double error_tolerance, const double min_dt,
                              const double max_dt)
{
    if (current_dt <= 0.0 || timestep_error <= 0.0 || error_tolerance <= 0.0 || min_dt <= 0.0 || max_dt < min_dt)
    {
        throw std::invalid_argument("Invalid SSP-RK3 timestep adaptation input.");
    }

    const double error_fraction = std::sqrt(0.5 * error_tolerance / timestep_error);
    const double time_step_factor = std::min(std::max(error_fraction, 0.3), 2.0);
    const double new_dt = 0.9 * time_step_factor * current_dt;
    return std::min(std::max(new_dt, min_dt), max_dt);
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
    if (!internal::all_finite(initial_state) || !std::isfinite(t0) || !std::isfinite(t_bound) || !std::isfinite(options.relative_tolerance) ||
        !std::isfinite(options.absolute_tolerance) || options.relative_tolerance <= 0.0 || options.absolute_tolerance <= 0.0 ||
        options.max_steps == 0)
    {
        return invalid_embedded_rk_result(initial_state, t0, "Invalid embedded Runge-Kutta input.");
    }
    if (t0 == t_bound)
    {
        return AdaptiveRkResult{initial_state, t0, AdaptiveRkStatus::success, "", 0, 0, 0, 0.0};
    }
    if ((!std::isfinite(options.initial_step) || options.initial_step <= 0.0) || !std::isfinite(options.max_step) || options.max_step < 0.0 ||
        !std::isfinite(options.min_step) || options.min_step < 0.0)
    {
        return invalid_embedded_rk_result(initial_state, t0, "Invalid embedded Runge-Kutta step-size option.");
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
    std::optional<std::vector<double>> first_derivative;

    AdaptiveRkResult result;
    result.state = y;
    result.t = t;
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
            result.t = t;
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

            internal::EmbeddedRkStep step = internal::take_embedded_rk_step(*tableau, rhs, t, y, h, first_derivative, function_evaluations);
            const bool step_is_finite =
                internal::all_finite(step.higher_state) && internal::all_finite(step.embedded_state) && internal::all_finite(step.error);
            const double norm = step_is_finite ? error_norm(step.error, y, step.higher_state, options) : std::numeric_limits<double>::infinity();
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
                y = std::move(step.higher_state);
                first_derivative = tableau->fsal ? std::optional<std::vector<double>>{std::move(step.final_derivative)} : std::nullopt;
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
                first_derivative = std::move(step.first_derivative);
                ++result.rejected_steps;
                internal::performance::record_embedded_rk_step(false);
                if (step_is_at_effective_lower_bound)
                {
                    result.state = y;
                    result.t = t;
                    result.status = AdaptiveRkStatus::step_size_underflow;
                    result.message = std::string(tableau->name) + " required a step smaller than min_step.";
                    result.function_evaluations = function_evaluations;
                    return result;
                }
            }
        }
    }

    result.state = std::move(y);
    result.t = t;
    result.status = AdaptiveRkStatus::success;
    result.function_evaluations = function_evaluations;
    return result;
}

AdaptiveRkResult integrate_rk23(const std::vector<double>& initial_state, const double t0, const double t_bound, const OdeRhs& rhs,
                                const AdaptiveRkOptions& options)
{
    return integrate_embedded_rk(initial_state, t0, t_bound, rhs, EmbeddedRkMethod::rk23, options);
}

} // namespace qmomkit
