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
    double observed_order = std::numeric_limits<double>::quiet_NaN();
};

struct EmbeddedFixedRun
{
    std::size_t steps = 0;
    double dt = 0.0;
    double solution_error = 0.0;
    double embedded_difference = 0.0;
    double solution_order = std::numeric_limits<double>::quiet_NaN();
    double embedded_order = std::numeric_limits<double>::quiet_NaN();
};

qmomkit::OdeRhs exponential_decay_rhs()
{
    return [](double, const std::vector<double>& state) { return std::vector<double>{-state[0]}; };
}

double observed_order(const double coarse_error, const double fine_error, const double coarse_step, const double fine_step)
{
    QMOMKIT_VERIFY_CHECK(coarse_error > 0.0);
    QMOMKIT_VERIFY_CHECK(fine_error > 0.0);
    QMOMKIT_VERIFY_CHECK(coarse_step > fine_step);
    return std::log(coarse_error / fine_error) / std::log(coarse_step / fine_step);
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

    for (const std::size_t steps : step_counts)
    {
        const double dt = final_time / static_cast<double>(steps);
        double t = 0.0;
        std::vector<double> state{1.0};
        for (std::size_t step = 0; step < steps; ++step)
        {
            const qmomkit::SspRk3Result result = qmomkit::ssp_rk3_step(state, t, dt, exponential_decay_rhs(), true);
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
        runs[run].observed_order = observed_order(runs[run - 1].error, runs[run].error, runs[run - 1].dt, runs[run].dt);
    }

    // The coarsest pair may still contain pre-asymptotic behavior, but the
    // last refinements should be clearly third order for this smooth scalar
    // problem.
    QMOMKIT_VERIFY_CHECK(runs[runs.size() - 2].observed_order > 2.85);
    QMOMKIT_VERIFY_CHECK(runs[runs.size() - 1].observed_order > 2.85);
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

        const qmomkit::AdaptiveRkResult result = qmomkit::integrate_rk23({1.0}, 0.0, final_time, exponential_decay_rhs(), options);
        QMOMKIT_VERIFY_CHECK(result.status == qmomkit::AdaptiveRkStatus::success);
        QMOMKIT_VERIFY_NEAR(result.t, final_time, 1.0e-14);
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
        runs[run].observed_order = observed_order(runs[run - 1].error, runs[run].error, runs[run - 1].average_step, runs[run].average_step);
    }

    // Bogacki-Shampine RK23 advances the accepted solution with the third-order
    // formula while controlling an embedded second-order error estimate. The
    // finest tolerance refinements should therefore show a clear high-order
    // trend when measured against the accepted average step size.
    QMOMKIT_VERIFY_CHECK(runs[runs.size() - 2].observed_order > 2.2);
    QMOMKIT_VERIFY_CHECK(runs[runs.size() - 1].observed_order > 2.2);
    return runs;
}

std::vector<EmbeddedFixedRun> verify_embedded_fixed_step_convergence(const qmomkit::EmbeddedRkMethod method, const double expected_solution_order,
                                                                     const double expected_embedded_order)
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
        const qmomkit::internal::EmbeddedRkStep one_step =
            qmomkit::internal::take_embedded_rk_step(tableau, exponential_decay_rhs(), 0.0, {1.0}, dt, std::nullopt, function_evaluations);
        QMOMKIT_VERIFY_CHECK(function_evaluations == tableau.stage_count);
        runs.push_back(EmbeddedFixedRun{steps, dt, std::abs(result.state[0] - exact), std::abs(one_step.error[0])});
    }

    for (std::size_t run = 1; run < runs.size(); ++run)
    {
        runs[run].solution_order = observed_order(runs[run - 1].solution_error, runs[run].solution_error, runs[run - 1].dt, runs[run].dt);
        runs[run].embedded_order = observed_order(runs[run - 1].embedded_difference, runs[run].embedded_difference, runs[run - 1].dt, runs[run].dt);
    }

    QMOMKIT_VERIFY_CHECK(runs.back().solution_order > expected_solution_order - 0.15);
    QMOMKIT_VERIFY_CHECK(runs.back().solution_order < expected_solution_order + 0.15);
    QMOMKIT_VERIFY_CHECK(runs.back().embedded_order > expected_embedded_order - 0.15);
    QMOMKIT_VERIFY_CHECK(runs.back().embedded_order < expected_embedded_order + 0.15);
    return runs;
}

void print_ssp_table(const std::vector<FixedStepRun>& runs)
{
    std::cout << "SSP-RK3 convergence for y'=-y, y(0)=1, T=1\n";
    std::cout << "steps dt error observed_order\n";
    for (const FixedStepRun& run : runs)
    {
        std::cout << run.steps << ' ' << run.dt << ' ' << run.error << ' ' << run.observed_order << '\n';
    }
}

void print_rk23_table(const std::vector<AdaptiveRun>& runs)
{
    std::cout << "Adaptive RK23 convergence for y'=-y, y(0)=1, T=2\n";
    std::cout << "tolerance accepted_steps rejected_steps function_evaluations average_step error observed_order\n";
    for (const AdaptiveRun& run : runs)
    {
        std::cout << run.tolerance << ' ' << run.accepted_steps << ' ' << run.rejected_steps << ' ' << run.function_evaluations << ' '
                  << run.average_step << ' ' << run.error << ' ' << run.observed_order << '\n';
    }
}

void print_embedded_fixed_step_table(const char* const method, const std::vector<EmbeddedFixedRun>& runs)
{
    std::cout << method << " fixed-step convergence for y'=-y, y(0)=1, T=1\n";
    std::cout << "steps dt solution_error embedded_difference solution_order embedded_order\n";
    for (const EmbeddedFixedRun& run : runs)
    {
        std::cout << run.steps << ' ' << run.dt << ' ' << run.solution_error << ' ' << run.embedded_difference << ' ' << run.solution_order << ' '
                  << run.embedded_order << '\n';
    }
}

} // namespace

int main()
{
    std::cout << std::setprecision(17);
    const std::vector<FixedStepRun> ssp_runs = verify_ssp_rk3_convergence();
    const std::vector<AdaptiveRun> rk23_runs = verify_rk23_convergence();
    const std::vector<EmbeddedFixedRun> rk23_fixed_runs = verify_embedded_fixed_step_convergence(qmomkit::EmbeddedRkMethod::rk23, 3.0, 3.0);
    const std::vector<EmbeddedFixedRun> rk45_fixed_runs = verify_embedded_fixed_step_convergence(qmomkit::EmbeddedRkMethod::rk45, 5.0, 5.0);
    print_ssp_table(ssp_runs);
    print_rk23_table(rk23_runs);
    print_embedded_fixed_step_table("RK23", rk23_fixed_runs);
    print_embedded_fixed_step_table("RK45", rk45_fixed_runs);
    return 0;
}
