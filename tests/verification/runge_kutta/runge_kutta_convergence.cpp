// Convergence order of the time integrators.
//
// What is verified
//   Two different things, deliberately kept apart.
//
//   Asymptotic convergence order, with the step size forced. The controller is
//   pinned to a single step (initial = min = max) and its tolerances set to
//   1.0e-15 so it cannot adapt, giving a genuine halving refinement. SSP-RK3,
//   RK23 and RK45 are each checked this way, along with the embedded error
//   estimate of the two pairs.
//
//   Controller behaviour, with the step size free. The requested tolerance is
//   tightened tenfold per run and the controller chooses its own steps. What is
//   checked here is that the achieved error stays within the tolerance asked
//   for, that it falls at every tightening, that the work does not decrease,
//   and that the error scales roughly in proportion to the tolerance.
//
// How the reference is constructed
//   The test problem is y' = -y from t = 0, whose exact solution exp(-t) is
//   known in closed form. The reference is that value, evaluated directly, so
//   it does not depend on any integrator.
//
// Error metric
//   Absolute error of the final value against exp(-t). Slopes are taken on
//   log-log axes: against the step size for the fixed-step studies, where the
//   slope is a convergence order, and against the requested tolerance for the
//   adaptive studies, where it is not.
//
// Pass criteria
//   Every slope is taken from the finest pair of runs alone - the estimate
//   between the second smallest and the smallest control value. Convergence
//   order is an asymptotic statement, so a coarser pair, or an average over the
//   sweep, would report a pre-asymptotic value and let it gate the result.
//
//   Fixed step, dt halved exactly:
//     SSP-RK3          observed order > 2.85, five percent below nominal 3
//     RK23             solution and embedded order both 3.0 +/- 0.15
//     RK45             solution and embedded order both 5.0 +/- 0.15
//   Free step, tolerance tightened tenfold:
//     RK23 and RK45    achieved error within the requested tolerance, and the
//                      error-versus-tolerance slope between 0.85 and 1.25
//
//   RK23 and RK45 are each checked on two numbers that come out equal - 3 for
//   RK23, 5 for RK45 - but that measure different things:
//
//     final error order          how fast the answer at the end of the whole
//       run approaches the exact value as the step is halved. This is the order
//       of the method: it is off by roughly h^(p+1) on each step and takes 1/h
//       steps, leaving roughly h^p at the end.
//
//     single step estimate order how fast the scheme's own error estimate for
//       one step shrinks as the step is halved. That estimate is the gap
//       between the two answers the scheme produces for the same step, and its
//       size is set by the less accurate of the two, which is off by roughly
//       h^p on a single step.
//
//   So both land on p, one by accumulating over a whole run and one from a
//   single step. Neither is the accuracy of the less accurate answer over a
//   whole run, which would be p-1; the scheme never forms that, and this study
//   does not measure it.
//
//   The adaptive studies deliberately assert no order. Measuring error against
//   the average accepted step there does not give the order of the method, it
//   gives the controller's step exponent, because the controller sets the step
//   from the tolerance as tol^(1/(q+1)) with q the embedded order. That
//   coincides with the nominal order only when the accepted step count is large
//   enough for its integer rounding to wash out. At these tolerances RK45 takes
//   5 to 14 steps and the measurement reads about 7 rather than 5; the fixed
//   step study above reads 5.03 on the same method.
//
// Output
//   Nothing on success. When a check fails, the reason and the full refinement
//   table go to stderr, because a bare pass/fail cannot distinguish an order
//   that grazed its threshold from one that collapsed.
//
// Reproduce
//   ctest -R verification_.*(ssp_rk3|rk23|rk45).*convergence
//
// Build modes
//   Both. Checks go through QMOMKIT_VERIFY_*, not assert().

#include "../verification_helpers.hpp"
#include "internal/embedded_runge_kutta.hpp"
#include "qmomkit/runge_kutta.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

namespace
{

struct FixedStepRun
{
    std::size_t steps = 0;
    double dt = 0.0;
    double error = 0.0;
    double observed_order = std::numeric_limits<double>::quiet_NaN();
};

struct AdaptiveRun
{
    double tolerance = 0.0;
    std::size_t accepted_steps = 0;
    std::size_t rejected_steps = 0;
    std::size_t function_evaluations = 0;
    double average_step = 0.0;
    double error = 0.0;
    // Slope of achieved error against requested tolerance. Not a convergence
    // order: the controller chooses the steps, so there is no refined h here.
    double tolerance_slope = std::numeric_limits<double>::quiet_NaN();
};

struct EmbeddedFixedRun
{
    std::size_t steps = 0;
    double dt = 0.0;
    // Error left at the end of the whole run.
    double final_error = 0.0;
    // The error estimate the scheme forms for a single step, which is the gap
    // between its two answers for that step.
    double single_step_estimate = 0.0;
    double final_error_order = std::numeric_limits<double>::quiet_NaN();
    double single_step_estimate_order = std::numeric_limits<double>::quiet_NaN();
};

qmomkit::OdeRhs exponential_decay_rhs()
{
    return [](double, const std::vector<double>& state, std::vector<double>& derivative) { derivative[0] = -state[0]; };
}

// Slope of one quantity against another on log-log axes. It is a convergence
// order only when the control variable is a refined step size; against a
// requested tolerance it measures how faithfully the controller delivers what
// it was asked for.
double log_slope(const double coarse_value, const double fine_value, const double coarse_control, const double fine_control)
{
    QMOMKIT_VERIFY_CHECK(coarse_value > 0.0);
    QMOMKIT_VERIFY_CHECK(fine_value > 0.0);
    QMOMKIT_VERIFY_CHECK(coarse_control > fine_control);
    return std::log(coarse_value / fine_value) / std::log(coarse_control / fine_control);
}

// Refinement tables are built only when a check fails. A bare "expression was
// false" cannot distinguish an order that grazed its threshold from one that
// collapsed, and that difference is the whole diagnosis; on a passing run the
// same table tells a reader nothing and CTest hides it anyway.
[[noreturn]] void fail_with_table(const std::string& reason, const std::string& header, const std::string& rows)
{
    std::cerr << reason << '\n' << header << '\n' << rows;
    std::exit(EXIT_FAILURE);
}

std::string fixed_step_rows(const std::vector<FixedStepRun>& runs)
{
    std::ostringstream out;
    out << std::setprecision(17);
    for (const FixedStepRun& run : runs)
    {
        out << run.steps << ' ' << run.dt << ' ' << run.error << ' ' << run.observed_order << '\n';
    }
    return out.str();
}

std::string adaptive_rows(const std::vector<AdaptiveRun>& runs)
{
    std::ostringstream out;
    out << std::setprecision(17);
    for (const AdaptiveRun& run : runs)
    {
        out << run.tolerance << ' ' << run.accepted_steps << ' ' << run.rejected_steps << ' ' << run.function_evaluations << ' ' << run.average_step
            << ' ' << run.error << ' ' << run.tolerance_slope << '\n';
    }
    return out.str();
}

std::string embedded_rows(const std::vector<EmbeddedFixedRun>& runs)
{
    std::ostringstream out;
    out << std::setprecision(17);
    for (const EmbeddedFixedRun& run : runs)
    {
        out << run.steps << ' ' << run.dt << ' ' << run.final_error << ' ' << run.single_step_estimate << ' ' << run.final_error_order << ' '
            << run.single_step_estimate_order << '\n';
    }
    return out.str();
}

std::vector<FixedStepRun> verify_ssp_rk3_convergence()
{
    // The fixed-step study integrates y' = -y from t=0 to T=1, whose exact
    // solution is y(T)=exp(-T). Repeated public SSP-RK3 steps should therefore
    // show third-order global convergence as dt is halved.
    constexpr double final_time = 1.0;
    const double exact = std::exp(-final_time);
    const std::vector<std::size_t> step_counts{5, 10, 20, 40, 80};
    std::vector<FixedStepRun> runs;
    qmomkit::SspRk3Options options;
    options.estimate_error = true;

    for (const std::size_t steps : step_counts)
    {
        const double dt = final_time / static_cast<double>(steps);
        double t = 0.0;
        std::vector<double> state{1.0};
        for (std::size_t step = 0; step < steps; ++step)
        {
            const qmomkit::SspRk3Result result = qmomkit::ssp_rk3_step(state, t, dt, exponential_decay_rhs(), options);
            QMOMKIT_VERIFY_CHECK(result.ok());
            QMOMKIT_VERIFY_CHECK(result.error_estimate.has_value());
            state = result.state;
            t += dt;
        }
        QMOMKIT_VERIFY_NEAR(t, final_time, 1.0e-14);
        runs.push_back(FixedStepRun{steps, dt, std::abs(state[0] - exact), std::numeric_limits<double>::quiet_NaN()});
    }

    for (std::size_t run = 1; run < runs.size(); ++run)
    {
        QMOMKIT_VERIFY_CHECK(runs[run].error < runs[run - 1].error);
        runs[run].observed_order = log_slope(runs[run - 1].error, runs[run].error, runs[run - 1].dt, runs[run].dt);
    }

    // Only the finest pair: the order is asymptotic, so the estimate from the
    // second smallest to the smallest step is the one that carries meaning.
    if (!(runs.back().observed_order > 2.85))
    {
        fail_with_table("SSP-RK3 fixed-step order " + std::to_string(runs.back().observed_order) + " did not exceed 2.85.",
                        "steps dt error observed_order", fixed_step_rows(runs));
    }
    return runs;
}

std::vector<AdaptiveRun> verify_rk23_convergence()
{
    // The adaptive study uses the same independent exact solution, but asks
    // RK23 to choose its own step sizes under successively tighter tolerances.
    // The average accepted step is not prescribed; it is used only to measure
    // the observed global order of the public adaptive integration result.
    constexpr double final_time = 2.0;
    const double exact = std::exp(-final_time);
    const std::vector<double> tolerances{1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7};
    std::vector<AdaptiveRun> runs;

    for (const double tolerance : tolerances)
    {
        qmomkit::AdaptiveRkOptions options;
        options.relative_tolerance = tolerance;
        options.absolute_tolerance = tolerance * 1.0e-3;
        options.initial_step = 0.1;
        options.max_step = final_time;

        const qmomkit::AdaptiveRkResult result =
            qmomkit::integrate_embedded_rk({1.0}, 0.0, final_time, exponential_decay_rhs(), qmomkit::EmbeddedRkMethod::rk23, options);
        QMOMKIT_VERIFY_CHECK(result.status == qmomkit::AdaptiveRkStatus::success);
        QMOMKIT_VERIFY_NEAR(result.final_time, final_time, 1.0e-14);
        QMOMKIT_VERIFY_CHECK(result.accepted_steps > 0);
        QMOMKIT_VERIFY_CHECK(result.function_evaluations >= 1 + 3 * result.accepted_steps);

        runs.push_back(AdaptiveRun{tolerance, result.accepted_steps, result.rejected_steps, result.function_evaluations,
                                   final_time / static_cast<double>(result.accepted_steps), std::abs(result.state[0] - exact),
                                   std::numeric_limits<double>::quiet_NaN()});
    }

    for (std::size_t run = 1; run < runs.size(); ++run)
    {
        QMOMKIT_VERIFY_CHECK(runs[run].error < runs[run - 1].error);
        QMOMKIT_VERIFY_CHECK(runs[run].accepted_steps >= runs[run - 1].accepted_steps);
        runs[run].tolerance_slope = log_slope(runs[run - 1].error, runs[run].error, runs[run - 1].tolerance, runs[run].tolerance);
        // The controller must deliver at least the accuracy it was asked for.
        QMOMKIT_VERIFY_CHECK(runs[run].error <= runs[run].tolerance);
    }

    // Tightening the requested tolerance tenfold should shrink the achieved
    // error about tenfold. The order of the method is established by the
    // fixed-step study instead, where a step size actually exists to refine.
    if (!(runs.back().tolerance_slope > 0.85 && runs.back().tolerance_slope < 1.25))
    {
        fail_with_table("RK23 error did not track the requested tolerance: slope " + std::to_string(runs.back().tolerance_slope) +
                            " is outside [0.85, 1.25].",
                        "tolerance accepted_steps rejected_steps function_evaluations average_step error tolerance_slope", adaptive_rows(runs));
    }
    return runs;
}

std::vector<AdaptiveRun> verify_rk45_adaptive_convergence()
{
    // Dormand-Prince RK45 advances the accepted fifth-order solution while
    // selecting its own steps from the embedded fourth-order estimate.  The
    // exact exponential solution supplies an independent error reference.
    constexpr double final_time = 2.0;
    const double exact = std::exp(-final_time);
    const std::vector<double> tolerances{1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7};
    std::vector<AdaptiveRun> runs;

    for (const double tolerance : tolerances)
    {
        qmomkit::AdaptiveRkOptions options;
        options.relative_tolerance = tolerance;
        options.absolute_tolerance = tolerance * 1.0e-3;
        options.initial_step = 0.1;
        options.max_step = final_time;

        const qmomkit::AdaptiveRkResult result =
            qmomkit::integrate_embedded_rk({1.0}, 0.0, final_time, exponential_decay_rhs(), qmomkit::EmbeddedRkMethod::rk45, options);
        QMOMKIT_VERIFY_CHECK(result.status == qmomkit::AdaptiveRkStatus::success);
        QMOMKIT_VERIFY_NEAR(result.final_time, final_time, 1.0e-14);
        QMOMKIT_VERIFY_CHECK(result.accepted_steps > 0);
        QMOMKIT_VERIFY_CHECK(result.function_evaluations >= 1 + 6 * result.accepted_steps);

        runs.push_back(AdaptiveRun{tolerance, result.accepted_steps, result.rejected_steps, result.function_evaluations,
                                   final_time / static_cast<double>(result.accepted_steps), std::abs(result.state[0] - exact),
                                   std::numeric_limits<double>::quiet_NaN()});
    }

    for (std::size_t run = 1; run < runs.size(); ++run)
    {
        QMOMKIT_VERIFY_CHECK(runs[run].error < runs[run - 1].error);
        QMOMKIT_VERIFY_CHECK(runs[run].accepted_steps >= runs[run - 1].accepted_steps);
        runs[run].tolerance_slope = log_slope(runs[run - 1].error, runs[run].error, runs[run - 1].tolerance, runs[run].tolerance);
        // The controller must deliver at least the accuracy it was asked for.
        QMOMKIT_VERIFY_CHECK(runs[run].error <= runs[run].tolerance);
    }

    // Same contract as RK23: error proportional to the requested tolerance. The
    // fifth-order claim belongs to the fixed-step study.
    if (!(runs.back().tolerance_slope > 0.85 && runs.back().tolerance_slope < 1.25))
    {
        fail_with_table("RK45 error did not track the requested tolerance: slope " + std::to_string(runs.back().tolerance_slope) +
                            " is outside [0.85, 1.25].",
                        "tolerance accepted_steps rejected_steps function_evaluations average_step error tolerance_slope", adaptive_rows(runs));
    }
    return runs;
}

std::vector<EmbeddedFixedRun> verify_embedded_fixed_step_convergence(const qmomkit::EmbeddedRkMethod method, const double expected_final_error_order,
                                                                     const double expected_single_step_estimate_order)
{
    constexpr double final_time = 1.0;
    const double exact = std::exp(-final_time);
    const std::vector<std::size_t> step_counts{8, 16, 32, 64};
    const qmomkit::internal::EmbeddedRkTableau& tableau = qmomkit::internal::embedded_rk_tableau(method);
    std::vector<EmbeddedFixedRun> runs;

    for (const std::size_t steps : step_counts)
    {
        const double dt = final_time / static_cast<double>(steps);
        qmomkit::AdaptiveRkOptions options;
        options.initial_step = dt;
        options.max_step = dt;
        options.min_step = dt;
        options.relative_tolerance = 1.0e-15;
        options.absolute_tolerance = 1.0e-15;
        const qmomkit::AdaptiveRkResult result = qmomkit::integrate_embedded_rk({1.0}, 0.0, final_time, exponential_decay_rhs(), method, options);
        QMOMKIT_VERIFY_CHECK(result.status == qmomkit::AdaptiveRkStatus::success);
        QMOMKIT_VERIFY_CHECK(result.accepted_steps == steps);

        std::size_t function_evaluations = 0;
        qmomkit::internal::EmbeddedRkWorkspace one_step;
        qmomkit::internal::take_embedded_rk_step(tableau, exponential_decay_rhs(), 0.0, {1.0}, dt, function_evaluations, one_step);
        QMOMKIT_VERIFY_CHECK(function_evaluations == tableau.stage_count);
        runs.push_back(EmbeddedFixedRun{steps, dt, std::abs(result.state[0] - exact), std::abs(one_step.error[0])});
    }

    for (std::size_t run = 1; run < runs.size(); ++run)
    {
        runs[run].final_error_order = log_slope(runs[run - 1].final_error, runs[run].final_error, runs[run - 1].dt, runs[run].dt);
        runs[run].single_step_estimate_order =
            log_slope(runs[run - 1].single_step_estimate, runs[run].single_step_estimate, runs[run - 1].dt, runs[run].dt);
    }

    const bool final_ok =
        runs.back().final_error_order > expected_final_error_order - 0.15 && runs.back().final_error_order < expected_final_error_order + 0.15;
    const bool estimate_ok = runs.back().single_step_estimate_order > expected_single_step_estimate_order - 0.15 &&
                             runs.back().single_step_estimate_order < expected_single_step_estimate_order + 0.15;
    if (!final_ok || !estimate_ok)
    {
        fail_with_table("Fixed-step order outside its band: final error order " + std::to_string(runs.back().final_error_order) + " expected " +
                            std::to_string(expected_final_error_order) + ", single step estimate order " +
                            std::to_string(runs.back().single_step_estimate_order) + " expected " +
                            std::to_string(expected_single_step_estimate_order) + ", both +/- 0.15.",
                        "steps dt final_error single_step_estimate final_error_order single_step_estimate_order", embedded_rows(runs));
    }
    return runs;
}

} // namespace

int main(const int argc, char** argv)
{
    if (argc > 2)
    {
        std::cerr << "Usage: test_runge_kutta_convergence [ssp_rk3|rk23_adaptive|rk23_fixed_step|rk45_adaptive|rk45_fixed_step]\n";
        return 2;
    }

    const std::string selector = argc == 2 ? argv[1] : "";
    if (selector.empty() || selector == "ssp_rk3")
    {
        verify_ssp_rk3_convergence();
        if (!selector.empty())
        {
            return 0;
        }
    }
    if (selector.empty() || selector == "rk23_adaptive")
    {
        verify_rk23_convergence();
        if (!selector.empty())
        {
            return 0;
        }
    }
    if (selector.empty() || selector == "rk23_fixed_step")
    {
        verify_embedded_fixed_step_convergence(qmomkit::EmbeddedRkMethod::rk23, 3.0, 3.0);
        if (!selector.empty())
        {
            return 0;
        }
    }
    if (selector.empty() || selector == "rk45_adaptive")
    {
        verify_rk45_adaptive_convergence();
        if (!selector.empty())
        {
            return 0;
        }
    }
    if (selector.empty() || selector == "rk45_fixed_step")
    {
        verify_embedded_fixed_step_convergence(qmomkit::EmbeddedRkMethod::rk45, 5.0, 5.0);
        return 0;
    }

    std::cerr << "Unknown Runge-Kutta convergence selector: " << selector << '\n';
    return 2;
}
