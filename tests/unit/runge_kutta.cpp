#include "qmomkit/runge_kutta.hpp"

#include "internal/embedded_runge_kutta.hpp"
#include "test_case_dispatch.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

void assert_near(const double actual, const double expected, const double tolerance)
{
    assert(std::abs(actual - expected) <= tolerance);
}

void assert_vector_near(const std::vector<double>& actual, const std::vector<double>& expected, const double tolerance)
{
    assert(actual.size() == expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        assert_near(actual[i], expected[i], tolerance);
    }
}

qmomkit::OdeRhs exponential_decay_rhs();

double reference_scalar_embedded_step(const qmomkit::EmbeddedRkMethod method, const double value, const double lambda, const double step,
                                      const bool higher)
{
    if (method == qmomkit::EmbeddedRkMethod::rk23)
    {
        const double k1 = lambda * value;
        const double k2 = lambda * (value + step * (1.0 / 2.0) * k1);
        const double k3 = lambda * (value + step * (3.0 / 4.0) * k2);
        const double high = value + step * ((2.0 / 9.0) * k1 + (1.0 / 3.0) * k2 + (4.0 / 9.0) * k3);
        const double k4 = lambda * high;
        return higher ? high : value + step * ((7.0 / 24.0) * k1 + (1.0 / 4.0) * k2 + (1.0 / 3.0) * k3 + (1.0 / 8.0) * k4);
    }

    const double k1 = lambda * value;
    const double k2 = lambda * (value + step * ((1.0 / 5.0) * k1));
    const double k3 = lambda * (value + step * ((3.0 / 40.0) * k1 + (9.0 / 40.0) * k2));
    const double k4 = lambda * (value + step * ((44.0 / 45.0) * k1 - (56.0 / 15.0) * k2 + (32.0 / 9.0) * k3));
    const double k5 = lambda * (value + step * ((19372.0 / 6561.0) * k1 - (25360.0 / 2187.0) * k2 + (64448.0 / 6561.0) * k3 - (212.0 / 729.0) * k4));
    const double k6 =
        lambda *
        (value + step * ((9017.0 / 3168.0) * k1 - (355.0 / 33.0) * k2 + (46732.0 / 5247.0) * k3 + (49.0 / 176.0) * k4 - (5103.0 / 18656.0) * k5));
    const double high =
        value + step * ((35.0 / 384.0) * k1 + (500.0 / 1113.0) * k3 + (125.0 / 192.0) * k4 - (2187.0 / 6784.0) * k5 + (11.0 / 84.0) * k6);
    const double k7 = lambda * high;
    return higher ? high
                  : value + step * ((5179.0 / 57600.0) * k1 + (7571.0 / 16695.0) * k3 + (393.0 / 640.0) * k4 - (92097.0 / 339200.0) * k5 +
                                    (187.0 / 2100.0) * k6 + (1.0 / 40.0) * k7);
}

void assert_tableau_invariants(const qmomkit::internal::EmbeddedRkTableau& tableau, const std::size_t expected_stages,
                               const std::size_t expected_higher_order, const std::size_t expected_embedded_order)
{
    assert(tableau.stage_count == expected_stages);
    assert(tableau.higher_order == expected_higher_order);
    assert(tableau.embedded_order == expected_embedded_order);
    assert(tableau.fsal);

    double higher_sum = 0.0;
    double embedded_sum = 0.0;
    for (std::size_t row = 0; row < tableau.stage_count; ++row)
    {
        double row_sum = 0.0;
        for (std::size_t column = 0; column < tableau.stage_count; ++column)
        {
            if (column >= row)
            {
                assert(tableau.a[row][column] == 0.0);
            }
            row_sum += tableau.a[row][column];
        }
        assert_near(row_sum, tableau.c[row], 2.0e-15);
        higher_sum += tableau.b[row];
        embedded_sum += tableau.b_embedded[row];
        assert(tableau.a[tableau.stage_count - 1][row] == tableau.b[row]);
    }
    assert_near(higher_sum, 1.0, 2.0e-15);
    assert_near(embedded_sum, 1.0, 2.0e-15);
}

void test_embedded_rk_tableaus_match_primary_papers()
{
    const qmomkit::internal::EmbeddedRkTableau& rk23 = qmomkit::internal::embedded_rk_tableau(qmomkit::EmbeddedRkMethod::rk23);
    assert_tableau_invariants(rk23, 4, 3, 2);
    assert((rk23.c == std::array<double, 7>{0.0, 1.0 / 2.0, 3.0 / 4.0, 1.0, 0.0, 0.0, 0.0}));
    assert((rk23.b == std::array<double, 7>{2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0, 0.0, 0.0, 0.0, 0.0}));
    assert((rk23.b_embedded == std::array<double, 7>{7.0 / 24.0, 1.0 / 4.0, 1.0 / 3.0, 1.0 / 8.0, 0.0, 0.0, 0.0}));
    assert(rk23.a[1][0] == 1.0 / 2.0);
    assert(rk23.a[2][1] == 3.0 / 4.0);
    assert(rk23.a[3][0] == 2.0 / 9.0);
    assert(rk23.a[3][1] == 1.0 / 3.0);
    assert(rk23.a[3][2] == 4.0 / 9.0);

    const qmomkit::internal::EmbeddedRkTableau& rk45 = qmomkit::internal::embedded_rk_tableau(qmomkit::EmbeddedRkMethod::rk45);
    assert_tableau_invariants(rk45, 7, 5, 4);
    assert((rk45.c == std::array<double, 7>{0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0}));
    assert((rk45.b == std::array<double, 7>{35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0}));
    assert((rk45.b_embedded ==
            std::array<double, 7>{5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0}));
    assert(rk45.a[1][0] == 1.0 / 5.0);
    assert(rk45.a[2][0] == 3.0 / 40.0 && rk45.a[2][1] == 9.0 / 40.0);
    assert(rk45.a[3][0] == 44.0 / 45.0 && rk45.a[3][1] == -56.0 / 15.0 && rk45.a[3][2] == 32.0 / 9.0);
    assert(rk45.a[4][0] == 19372.0 / 6561.0 && rk45.a[4][1] == -25360.0 / 2187.0 && rk45.a[4][2] == 64448.0 / 6561.0 &&
           rk45.a[4][3] == -212.0 / 729.0);
    assert(rk45.a[5][0] == 9017.0 / 3168.0 && rk45.a[5][1] == -355.0 / 33.0 && rk45.a[5][2] == 46732.0 / 5247.0 && rk45.a[5][3] == 49.0 / 176.0 &&
           rk45.a[5][4] == -5103.0 / 18656.0);
    assert(rk45.a[6][0] == 35.0 / 384.0 && rk45.a[6][1] == 0.0 && rk45.a[6][2] == 500.0 / 1113.0 && rk45.a[6][3] == 125.0 / 192.0 &&
           rk45.a[6][4] == -2187.0 / 6784.0 && rk45.a[6][5] == 11.0 / 84.0);

    assert_near(qmomkit::internal::embedded_rk_controller_exponent(rk23), -1.0 / 3.0, 0.0);
    assert_near(qmomkit::internal::embedded_rk_controller_exponent(rk45), -1.0 / 5.0, 0.0);
}

void test_embedded_rk_one_steps_match_independent_scalar_formulae()
{
    for (const qmomkit::EmbeddedRkMethod method : {qmomkit::EmbeddedRkMethod::rk23, qmomkit::EmbeddedRkMethod::rk45})
    {
        const qmomkit::internal::EmbeddedRkTableau& tableau = qmomkit::internal::embedded_rk_tableau(method);
        std::size_t evaluations = 0;
        qmomkit::internal::EmbeddedRkWorkspace forward;
        qmomkit::internal::take_embedded_rk_step(
            tableau,
            [](double, const std::vector<double>& state, std::vector<double>& derivative)
            {
                derivative[0] = state[0];
                derivative[1] = -2.0 * state[1];
            },
            0.25, {1.5, -0.5}, 0.2, evaluations, forward);
        const std::vector<double> expected_high{reference_scalar_embedded_step(method, 1.5, 1.0, 0.2, true),
                                                reference_scalar_embedded_step(method, -0.5, -2.0, 0.2, true)};
        const std::vector<double> expected_low{reference_scalar_embedded_step(method, 1.5, 1.0, 0.2, false),
                                               reference_scalar_embedded_step(method, -0.5, -2.0, 0.2, false)};
        assert_vector_near(forward.higher_state, expected_high, 2.0e-15);
        assert_vector_near(forward.embedded_state, expected_low, 2.0e-15);
        assert_vector_near(forward.error, {expected_high[0] - expected_low[0], expected_high[1] - expected_low[1]}, 2.0e-15);
        assert(evaluations == tableau.stage_count);

        evaluations = 0;
        qmomkit::internal::EmbeddedRkWorkspace backward;
        qmomkit::internal::take_embedded_rk_step(
            tableau, [](double, const std::vector<double>& state, std::vector<double>& derivative) { derivative[0] = state[0]; }, 0.25, {1.5}, -0.2,
            evaluations, backward);
        assert_near(backward.higher_state[0], reference_scalar_embedded_step(method, 1.5, 1.0, -0.2, true), 2.0e-15);
    }
}

void test_embedded_rk_one_steps_use_stage_abscissae_for_time_dependent_rhs()
{
    constexpr double initial_time = 0.3;
    constexpr double initial_value = 1.0;
    constexpr double step = 0.2;
    const double expected = initial_value + initial_time * step + 0.5 * step * step;

    for (const qmomkit::EmbeddedRkMethod method : {qmomkit::EmbeddedRkMethod::rk23, qmomkit::EmbeddedRkMethod::rk45})
    {
        const qmomkit::internal::EmbeddedRkTableau& tableau = qmomkit::internal::embedded_rk_tableau(method);
        std::size_t evaluations = 0;
        qmomkit::internal::EmbeddedRkWorkspace result;
        qmomkit::internal::take_embedded_rk_step(
            tableau, [](const double time, const std::vector<double>&, std::vector<double>& derivative) { derivative[0] = time; }, initial_time,
            {initial_value}, step, evaluations, result);

        // y' = t has y(t + h) = y(t) + t h + h^2 / 2.  Both members of
        // these pairs must use their declared stage abscissae to recover it.
        assert_near(result.higher_state[0], expected, 2.0e-15);
        assert_near(result.embedded_state[0], expected, 2.0e-15);
        assert_near(result.error[0], 0.0, 2.0e-15);
        assert(evaluations == tableau.stage_count);
    }
}

void test_embedded_rk_workspace_reuse_is_transparent()
{
    for (const qmomkit::EmbeddedRkMethod method : {qmomkit::EmbeddedRkMethod::rk23, qmomkit::EmbeddedRkMethod::rk45})
    {
        const qmomkit::internal::EmbeddedRkTableau& tableau = qmomkit::internal::embedded_rk_tableau(method);
        const qmomkit::OdeRhs rhs = [](double, const std::vector<double>& state, std::vector<double>& derivative)
        {
            derivative[0] = state[0];
            derivative[1] = -2.0 * state[1];
        };
        const std::vector<double> start{1.5, -0.5};

        std::size_t fresh_evaluations = 0;
        qmomkit::internal::EmbeddedRkWorkspace fresh;
        qmomkit::internal::take_embedded_rk_step(tableau, rhs, 0.25, start, 0.2, fresh_evaluations, fresh);
        assert(fresh_evaluations == tableau.stage_count);
        assert(fresh.has_first_derivative);

        // A workspace that already served an unrelated step, and whose cache the
        // caller therefore invalidated, must take the step exactly as a fresh one
        // does: reused storage may not carry values between attempts.
        std::size_t reused_evaluations = 0;
        qmomkit::internal::EmbeddedRkWorkspace reused;
        qmomkit::internal::take_embedded_rk_step(tableau, rhs, 0.1, {2.0, 1.0}, 0.3, reused_evaluations, reused);
        reused.has_first_derivative = false;
        reused_evaluations = 0;
        qmomkit::internal::take_embedded_rk_step(tableau, rhs, 0.25, start, 0.2, reused_evaluations, reused);
        assert(reused_evaluations == tableau.stage_count);
        assert_vector_near(reused.higher_state, fresh.higher_state, 0.0);
        assert_vector_near(reused.embedded_state, fresh.embedded_state, 0.0);
        assert_vector_near(reused.error, fresh.error, 0.0);

        // Retaking the step the workspace just took is the rejected-attempt retry:
        // the cached slope at the unchanged starting state replaces exactly one
        // evaluation and reproduces the same step bit for bit.
        std::size_t retry_evaluations = 0;
        qmomkit::internal::take_embedded_rk_step(tableau, rhs, 0.25, start, 0.2, retry_evaluations, reused);
        assert(retry_evaluations == tableau.stage_count - 1);
        assert(reused.has_first_derivative);
        assert_vector_near(reused.higher_state, fresh.higher_state, 0.0);
        assert_vector_near(reused.embedded_state, fresh.embedded_state, 0.0);
        assert_vector_near(reused.error, fresh.error, 0.0);
    }
}

void test_adaptive_initial_step_default_and_override()
{
    qmomkit::AdaptiveRkOptions options;
    assert(options.initial_step == qmomkit::default_adaptive_initial_step);
    assert(options.initial_step == 1.0e-6);

    options.initial_step = 2.5e-5;
    assert(options.initial_step == 2.5e-5);
}

void test_embedded_rk_default_initial_step_is_used_for_both_methods()
{
    for (const qmomkit::EmbeddedRkMethod method : {qmomkit::EmbeddedRkMethod::rk23, qmomkit::EmbeddedRkMethod::rk45})
    {
        std::vector<double> rhs_times;
        const double initial_time = 2.0;
        const qmomkit::AdaptiveRkResult result =
            qmomkit::integrate_embedded_rk({0.0}, initial_time, initial_time + 2.0 * qmomkit::default_adaptive_initial_step,
                                           [&rhs_times](const double time, const std::vector<double>&, std::vector<double>& derivative)
                                           {
                                               rhs_times.push_back(time);
                                               derivative[0] = 1.0;
                                           },
                                           method, {});

        assert(result.status == qmomkit::AdaptiveRkStatus::success);
        assert_near(result.final_time, initial_time + 2.0 * qmomkit::default_adaptive_initial_step, 1.0e-15);
        assert_near(result.state[0], 2.0 * qmomkit::default_adaptive_initial_step, 1.0e-15);
        assert(rhs_times.size() >= 2);
        assert_near(rhs_times[0], initial_time, 0.0);
        assert_near(rhs_times[1] - initial_time, qmomkit::internal::embedded_rk_tableau(method).c[1] * qmomkit::default_adaptive_initial_step,
                    1.0e-15);
    }
}

void test_embedded_rk_fsal_reuses_accepted_endpoint_derivative()
{
    for (const qmomkit::EmbeddedRkMethod method : {qmomkit::EmbeddedRkMethod::rk23, qmomkit::EmbeddedRkMethod::rk45})
    {
        qmomkit::AdaptiveRkOptions options;
        options.initial_step = 0.5;
        options.max_step = 0.5;
        const qmomkit::AdaptiveRkResult result = qmomkit::integrate_embedded_rk(
            {1.0}, 0.0, 1.0, [](double, const std::vector<double>& state, std::vector<double>& derivative) { derivative[0] = 0.0 * state[0]; },
            method, options);
        assert(result.status == qmomkit::AdaptiveRkStatus::success);
        assert(result.accepted_steps == 2);
        assert(result.function_evaluations == 2 * qmomkit::internal::embedded_rk_tableau(method).stage_count - 1);
    }
}

void test_embedded_rk_adaptive_controls_for_both_methods()
{
    for (const qmomkit::EmbeddedRkMethod method : {qmomkit::EmbeddedRkMethod::rk23, qmomkit::EmbeddedRkMethod::rk45})
    {
        const qmomkit::AdaptiveRkResult empty_interval = qmomkit::integrate_embedded_rk({1.0}, 0.25, 0.25, exponential_decay_rhs(), method, {});
        assert(empty_interval.status == qmomkit::AdaptiveRkStatus::success);
        assert(empty_interval.function_evaluations == 0);

        qmomkit::AdaptiveRkOptions accepted;
        accepted.initial_step = 0.2;
        accepted.max_step = 0.2;
        const qmomkit::AdaptiveRkResult first_trial = qmomkit::integrate_embedded_rk(
            {2.0}, 0.0, 0.2, [](double, const std::vector<double>&, std::vector<double>& derivative) { derivative[0] = 3.0; }, method, accepted);
        assert(first_trial.status == qmomkit::AdaptiveRkStatus::success);
        assert(first_trial.accepted_steps == 1);
        assert(first_trial.rejected_steps == 0);
        assert_near(first_trial.final_time, 0.2, 1.0e-15);
        assert_near(first_trial.state[0], 2.6, 1.0e-15);

        qmomkit::AdaptiveRkOptions rejected;
        rejected.initial_step = 1.0;
        rejected.max_step = 1.0;
        rejected.relative_tolerance = 1.0e-9;
        rejected.absolute_tolerance = 1.0e-12;
        const qmomkit::AdaptiveRkResult recovered = qmomkit::integrate_embedded_rk(
            {1.0}, 0.0, 1.0, [](double, const std::vector<double>& state, std::vector<double>& derivative) { derivative[0] = -10.0 * state[0]; },
            method, rejected);
        assert(recovered.status == qmomkit::AdaptiveRkStatus::success);
        assert(recovered.rejected_steps > 1);
        const std::size_t stage_count = qmomkit::internal::embedded_rk_tableau(method).stage_count;
        assert(recovered.function_evaluations == stage_count + (recovered.accepted_steps + recovered.rejected_steps - 1) * (stage_count - 1));
        assert_near(recovered.final_time, 1.0, 1.0e-14);

        qmomkit::AdaptiveRkOptions clamped;
        clamped.initial_step = 1.0;
        clamped.max_step = 0.25;
        const qmomkit::AdaptiveRkResult maximum_step = qmomkit::integrate_embedded_rk(
            {1.0}, 0.0, 1.0, [](double, const std::vector<double>& state, std::vector<double>& derivative) { derivative[0] = 0.0 * state[0]; },
            method, clamped);
        assert(maximum_step.status == qmomkit::AdaptiveRkStatus::success);
        assert(maximum_step.accepted_steps == 4);
        assert_near(maximum_step.final_time, 1.0, 1.0e-15);

        qmomkit::AdaptiveRkOptions nonfinite_error;
        nonfinite_error.initial_step = 2.0;
        nonfinite_error.max_step = 2.0;
        nonfinite_error.min_step = 0.5;
        const qmomkit::AdaptiveRkResult underflow = qmomkit::integrate_embedded_rk(
            {0.0}, 0.0, 2.0, [](double, const std::vector<double>&, std::vector<double>& derivative)
            { derivative[0] = std::numeric_limits<double>::max(); }, method, nonfinite_error);
        assert(underflow.status == qmomkit::AdaptiveRkStatus::step_size_underflow);
        assert(underflow.rejected_steps > 0);

        bool nonfinite_derivative_rejected = false;
        try
        {
            (void)qmomkit::integrate_embedded_rk(
                {1.0}, 0.0, 1.0, [](double, const std::vector<double>&, std::vector<double>& derivative)
                { derivative[0] = std::numeric_limits<double>::infinity(); }, method, accepted);
        }
        catch (const std::invalid_argument&)
        {
            nonfinite_derivative_rejected = true;
        }
        assert(nonfinite_derivative_rejected);
    }
}

void test_embedded_rk_public_backward_integration()
{
    for (const qmomkit::EmbeddedRkMethod method : {qmomkit::EmbeddedRkMethod::rk23, qmomkit::EmbeddedRkMethod::rk45})
    {
        qmomkit::AdaptiveRkOptions options;
        options.initial_step = 0.3;
        options.max_step = 0.3;
        std::vector<double> rhs_times;
        const qmomkit::AdaptiveRkResult result = qmomkit::integrate_embedded_rk(
            {2.0}, 2.0, 1.0,
            [&rhs_times](const double time, const std::vector<double>&, std::vector<double>& derivative)
            {
                rhs_times.push_back(time);
                derivative[0] = 1.0;
            },
            method, options);

        assert(result.status == qmomkit::AdaptiveRkStatus::success);
        assert_near(result.final_time, 1.0, 0.0);
        assert_near(result.state[0], 1.0, 1.0e-15);
        assert(rhs_times.size() >= 2);
        assert(rhs_times[1] < rhs_times[0]);
    }
}

qmomkit::OdeRhs exponential_decay_rhs()
{
    return [](double, const std::vector<double>& y, std::vector<double>& derivative) { derivative[0] = -y[0]; };
}

void test_ssp_rk3_accept_all_predicate()
{
    qmomkit::SspRk3Options options;
    options.accept = [](const std::vector<double>&) { return true; };
    options.max_retries = 4;

    const qmomkit::SspRk3Result result = qmomkit::ssp_rk3_step({2.0}, 0.0, 0.1, exponential_decay_rhs(), options);

    assert(result.ok());
    assert(result.status == qmomkit::SspRk3Status::success);
    assert(result.retries == 0);
    assert_near(result.time_step_used, 0.1, 1.0e-15);
    assert_vector_near(result.state, {1.8096666666666668}, 1.0e-15);
}

void test_ssp_rk3_reject_negative_values_then_success()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&, std::vector<double>& derivative) { derivative[0] = -2.0; };
    qmomkit::SspRk3Options options;
    options.accept = [](const std::vector<double>& state) { return state[0] >= 0.0; };
    options.max_retries = 2;

    const qmomkit::SspRk3Result result = qmomkit::ssp_rk3_step({1.0}, 0.0, 1.0, rhs, options);

    assert(result.ok());
    assert(result.retries == 1);
    assert_near(result.time_step_used, 0.5, 1.0e-15);
    assert_vector_near(result.state, {0.0}, 1.0e-15);
}

void test_ssp_rk3_multiple_retries()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&, std::vector<double>& derivative) { derivative[0] = -8.0; };
    qmomkit::SspRk3Options options;
    options.accept = [](const std::vector<double>& state) { return state[0] >= 0.0; };
    options.max_retries = 4;

    const qmomkit::SspRk3Result result = qmomkit::ssp_rk3_step({1.0}, 0.0, 1.0, rhs, options);

    assert(result.ok());
    assert(result.retries == 3);
    assert_near(result.time_step_used, 0.125, 1.0e-15);
    assert_vector_near(result.state, {0.0}, 1.0e-15);
}

void test_ssp_rk3_failure_after_max_retries()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&, std::vector<double>& derivative) { derivative[0] = -8.0; };
    qmomkit::SspRk3Options options;
    options.accept = [](const std::vector<double>& state) { return state[0] >= 0.0; };
    options.max_retries = 2;

    const qmomkit::SspRk3Result result = qmomkit::ssp_rk3_step({1.0}, 0.0, 1.0, rhs, options);

    assert(!result.ok());
    assert(result.status == qmomkit::SspRk3Status::rejected_by_predicate);
    assert(result.retries == 2);
    assert_near(result.time_step_used, 0.25, 1.0e-15);
    assert_vector_near(result.state, {1.0}, 1.0e-15);
}

void test_ssp_rk3_failure_after_reaching_minimum_timestep()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&, std::vector<double>& derivative) { derivative[0] = -8.0; };
    qmomkit::SspRk3Options options;
    options.accept = [](const std::vector<double>& state) { return state[0] >= 0.0; };
    options.max_retries = 10;
    options.min_time_step = 0.2;

    const qmomkit::SspRk3Result result = qmomkit::ssp_rk3_step({1.0}, 0.0, 1.0, rhs, options);

    assert(!result.ok());
    assert(result.status == qmomkit::SspRk3Status::minimum_timestep_reached);
    assert(result.retries == 2);
    assert_near(result.time_step_used, 0.25, 1.0e-15);
    assert_vector_near(result.state, {1.0}, 1.0e-15);
}

void test_ssp_rk3_rejects_initial_step_below_minimum()
{
    qmomkit::SspRk3Options options;
    options.min_time_step = 0.2;

    bool threw = false;
    try
    {
        (void)qmomkit::ssp_rk3_step({1.0}, 0.0, 0.1, exponential_decay_rhs(), options);
    }
    catch (const std::invalid_argument&)
    {
        threw = true;
    }

    assert(threw);
}

void test_ssp_rk3_rejects_retry_that_underflows_to_zero()
{
    std::size_t rhs_calls = 0;
    qmomkit::SspRk3Options options;
    options.accept = [](const std::vector<double>&) { return false; };
    options.max_retries = 1;

    const qmomkit::SspRk3Result result = qmomkit::ssp_rk3_step(
        {0.0}, 0.0, std::numeric_limits<double>::denorm_min(),
        [&rhs_calls](double, const std::vector<double>&, std::vector<double>& derivative)
        {
            ++rhs_calls;
            derivative[0] = 0.0;
        },
        options);

    assert(!result.ok());
    assert(result.status == qmomkit::SspRk3Status::minimum_timestep_reached);
    assert(result.retries == 0);
    assert(result.time_step_used == std::numeric_limits<double>::denorm_min());
    assert(rhs_calls == 3);
}

void test_ssp_rk3_without_stage_callback_accepts_unchecked_candidate()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&, std::vector<double>& derivative)
    {
        derivative[0] = 0.0;
        derivative[1] = 0.0;
        derivative[2] = -2.0;
    };

    const qmomkit::SspRk3Result result = qmomkit::ssp_rk3_step({1.0, 0.0, 1.0}, 0.0, 1.0, rhs, qmomkit::SspRk3Options{});

    assert(result.ok());
    assert(result.retries == 0);
    assert_vector_near(result.state, {1.0, 0.0, -1.0}, 1.0e-15);
}

void test_ssp_rk3_final_only_callback_invoked_once()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&, std::vector<double>& derivative)
    {
        derivative[0] = 0.0;
        derivative[1] = 0.0;
        derivative[2] = -2.0;
    };
    int calls = 0;
    qmomkit::SspRk3Options options;
    options.accept = [&](const std::vector<double>& state)
    {
        ++calls;
        assert_vector_near(state, {1.0, 0.0, -1.0}, 1.0e-15);
        return true;
    };

    const qmomkit::SspRk3Result result = qmomkit::ssp_rk3_step({1.0, 0.0, 1.0}, 0.0, 1.0, rhs, options);

    assert(result.ok());
    assert(result.retries == 0);
    assert(calls == 1);
}

void test_ssp_rk3_stage_callback_invoked_for_each_stage()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&, std::vector<double>& derivative) { derivative[0] = 1.0; };
    std::vector<qmomkit::SspRk3Stage> stages;
    std::vector<double> states;
    qmomkit::SspRk3Options options;
    options.stage_accept = [&](const qmomkit::SspRk3Stage stage, const std::vector<double>& state)
    {
        stages.push_back(stage);
        states.push_back(state[0]);
        return true;
    };

    const qmomkit::SspRk3Result result = qmomkit::ssp_rk3_step({0.0}, 0.0, 1.0, rhs, options);

    assert(result.ok());
    assert(result.retries == 0);
    assert(stages.size() == 3);
    assert(stages[0] == qmomkit::SspRk3Stage::Stage1);
    assert(stages[1] == qmomkit::SspRk3Stage::Stage2);
    assert(stages[2] == qmomkit::SspRk3Stage::Final);
    assert_vector_near(states, {1.0, 0.5, 1.0}, 1.0e-15);
}

void test_ssp_rk3_failed_stage_check_retries_with_smaller_timestep()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&, std::vector<double>& derivative) { derivative[0] = -2.0; };
    int stage_1_calls = 0;
    qmomkit::SspRk3Options options;
    options.stage_accept = [&](const qmomkit::SspRk3Stage stage, const std::vector<double>& state)
    {
        if (stage == qmomkit::SspRk3Stage::Stage1)
        {
            ++stage_1_calls;
        }
        return state[0] >= 0.0;
    };
    options.max_retries = 2;

    const qmomkit::SspRk3Result result = qmomkit::ssp_rk3_step({1.0}, 0.0, 1.0, rhs, options);

    assert(result.ok());
    assert(result.retries == 1);
    assert(stage_1_calls == 2);
    assert_near(result.time_step_used, 0.5, 1.0e-15);
    assert_vector_near(result.state, {0.0}, 1.0e-15);
}

void test_ssp_rk3_stage2_rejection_skips_k3_and_retries_from_original_state()
{
    std::vector<double> rhs_times;
    std::vector<double> rhs_states;
    const qmomkit::OdeRhs rhs = [&](const double time, const std::vector<double>& state, std::vector<double>& derivative)
    {
        rhs_times.push_back(time);
        rhs_states.push_back(state[0]);
        derivative[0] = -0.8;
    };

    std::vector<qmomkit::SspRk3Stage> checked_stages;
    std::vector<double> checked_states;
    qmomkit::SspRk3Options options;
    options.stage_accept = [&](const qmomkit::SspRk3Stage stage, const std::vector<double>& state)
    {
        checked_stages.push_back(stage);
        checked_states.push_back(state[0]);
        if (stage == qmomkit::SspRk3Stage::Stage2)
        {
            return state[0] >= 0.75;
        }
        return true;
    };
    options.max_retries = 2;

    const qmomkit::SspRk3Result result = qmomkit::ssp_rk3_step({1.0}, 0.0, 1.0, rhs, options);

    assert(result.ok());
    assert(result.retries == 1);
    assert_near(result.time_step_used, 0.5, 1.0e-15);
    assert_vector_near(result.state, {0.6}, 1.0e-15);

    assert(rhs_times.size() == 5);
    assert_vector_near(rhs_times, {0.0, 1.0, 0.0, 0.5, 0.25}, 1.0e-15);
    assert_vector_near(rhs_states, {1.0, 0.2, 1.0, 0.6, 0.8}, 1.0e-15);

    assert(checked_stages.size() == 5);
    assert(checked_stages[0] == qmomkit::SspRk3Stage::Stage1);
    assert(checked_stages[1] == qmomkit::SspRk3Stage::Stage2);
    assert(checked_stages[2] == qmomkit::SspRk3Stage::Stage1);
    assert(checked_stages[3] == qmomkit::SspRk3Stage::Stage2);
    assert(checked_stages[4] == qmomkit::SspRk3Stage::Final);
    assert_vector_near(checked_states, {0.2, 0.6, 0.6, 0.8, 0.6}, 1.0e-15);
}

void test_ssp_rk3_final_stage_rejection_short_circuits_final_accept()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&, std::vector<double>& derivative) { derivative[0] = 1.0; };
    int final_accept_calls = 0;
    qmomkit::SspRk3Options options;
    options.stage_accept = [](const qmomkit::SspRk3Stage stage, const std::vector<double>&) { return stage != qmomkit::SspRk3Stage::Final; };
    options.accept = [&](const std::vector<double>&)
    {
        ++final_accept_calls;
        return true;
    };

    const qmomkit::SspRk3Result result = qmomkit::ssp_rk3_step({0.0}, 0.0, 1.0, rhs, options);

    assert(!result.ok());
    assert(result.status == qmomkit::SspRk3Status::rejected_by_predicate);
    assert(result.retries == 0);
    assert_near(result.time_step_used, 1.0, 1.0e-15);
    assert_vector_near(result.state, {0.0}, 1.0e-15);
    assert(final_accept_calls == 0);
}

void test_ssp_rk3_stage_rejection_exhaustion_preserves_original_state()
{
    const std::vector<double> initial{2.0};
    qmomkit::SspRk3Options options;
    options.stage_accept = [](qmomkit::SspRk3Stage, const std::vector<double>&) { return false; };
    options.max_retries = 1;

    const qmomkit::SspRk3Result result = qmomkit::ssp_rk3_step(
        initial, 0.0, 0.8, [](double, const std::vector<double>&, std::vector<double>& derivative) { derivative[0] = 3.0; }, options);

    assert(!result.ok());
    assert(result.status == qmomkit::SspRk3Status::rejected_by_predicate);
    assert(result.retries == 1);
    assert_near(result.time_step_used, 0.4, 1.0e-15);
    assert_vector_near(result.state, initial, 1.0e-15);
}

void test_ssp_rk3_stage_rejection_minimum_timestep_preserves_original_state()
{
    const std::vector<double> initial{2.0};
    qmomkit::SspRk3Options options;
    options.stage_accept = [](qmomkit::SspRk3Stage, const std::vector<double>&) { return false; };
    options.max_retries = 4;
    options.min_time_step = 0.5;

    const qmomkit::SspRk3Result result = qmomkit::ssp_rk3_step(
        initial, 0.0, 0.8, [](double, const std::vector<double>&, std::vector<double>& derivative) { derivative[0] = 3.0; }, options);

    assert(!result.ok());
    assert(result.status == qmomkit::SspRk3Status::minimum_timestep_reached);
    assert(result.retries == 0);
    assert_near(result.time_step_used, 0.8, 1.0e-15);
    assert_vector_near(result.state, initial, 1.0e-15);
}

void test_ssp_rk3_rhs_stage_audit_preserves_times_and_retry_states()
{
    std::vector<double> rhs_times;
    std::vector<double> rhs_states;
    const qmomkit::OdeRhs rhs = [&](const double time, const std::vector<double>& state, std::vector<double>& derivative)
    {
        rhs_times.push_back(time);
        rhs_states.push_back(state[0]);
        derivative[0] = 2.0;
    };

    int final_stage_calls = 0;
    qmomkit::SspRk3Options options;
    options.stage_accept = [&](const qmomkit::SspRk3Stage stage, const std::vector<double>&)
    {
        if (stage == qmomkit::SspRk3Stage::Final)
        {
            ++final_stage_calls;
            return final_stage_calls > 1;
        }
        return true;
    };
    options.max_retries = 1;

    const qmomkit::SspRk3Result result = qmomkit::ssp_rk3_step({1.0}, 0.25, 0.4, rhs, options);

    assert(result.ok());
    assert(result.retries == 1);
    assert_near(result.time_step_used, 0.2, 1.0e-15);
    assert_vector_near(result.state, {1.4}, 1.0e-15);

    assert(rhs_times.size() == 6);
    assert_vector_near(rhs_times, {0.25, 0.65, 0.45, 0.25, 0.45, 0.35}, 1.0e-15);
    assert_vector_near(rhs_states, {1.0, 1.8, 1.4, 1.0, 1.4, 1.2}, 1.0e-15);
}

void test_rk23_rejects_overly_large_initial_step()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>& y, std::vector<double>& derivative) { derivative[0] = -10.0 * y[0]; };
    qmomkit::AdaptiveRkOptions options;
    options.initial_step = 1.0;
    options.max_step = 1.0;
    options.relative_tolerance = 1.0e-9;
    options.absolute_tolerance = 1.0e-12;

    const qmomkit::AdaptiveRkResult result = qmomkit::integrate_embedded_rk({1.0}, 0.0, 1.0, rhs, qmomkit::EmbeddedRkMethod::rk23, options);

    assert(result.status == qmomkit::AdaptiveRkStatus::success);
    assert(result.rejected_steps > 0);
    assert_near(result.state[0], std::exp(-10.0), 5.0e-9);
}

void test_rk23_zero_rhs_has_no_timestep_pathology()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&, std::vector<double>& derivative)
    { std::fill(derivative.begin(), derivative.end(), 0.0); };
    qmomkit::AdaptiveRkOptions options;
    options.initial_step = 0.1;
    options.relative_tolerance = 1.0e-8;
    options.absolute_tolerance = 1.0e-12;

    const qmomkit::AdaptiveRkResult result =
        qmomkit::integrate_embedded_rk({0.0, 2.0, -3.0}, 0.0, 1.0, rhs, qmomkit::EmbeddedRkMethod::rk23, options);

    assert(result.status == qmomkit::AdaptiveRkStatus::success);
    assert_near(result.final_time, 1.0, 1.0e-14);
    assert_vector_near(result.state, {0.0, 2.0, -3.0}, 1.0e-15);
    assert(result.rejected_steps == 0);
}

void test_rk23_invalid_options_fail_before_loop()
{
    qmomkit::AdaptiveRkOptions options;
    options.absolute_tolerance = 0.0;

    const qmomkit::AdaptiveRkResult zero_atol =
        qmomkit::integrate_embedded_rk({0.0}, 0.0, 1.0, exponential_decay_rhs(), qmomkit::EmbeddedRkMethod::rk23, options);
    assert(zero_atol.status == qmomkit::AdaptiveRkStatus::invalid_input);

    for (const double initial_step :
         {0.0, -1.0, std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
    {
        options.absolute_tolerance = 1.0e-9;
        options.initial_step = initial_step;
        for (const qmomkit::EmbeddedRkMethod method : {qmomkit::EmbeddedRkMethod::rk23, qmomkit::EmbeddedRkMethod::rk45})
        {
            const qmomkit::AdaptiveRkResult result = qmomkit::integrate_embedded_rk({1.0}, 0.0, 1.0, exponential_decay_rhs(), method, options);
            assert(result.status == qmomkit::AdaptiveRkStatus::invalid_input);
        }
    }
}

// An empty interval needs no step, but invalid step options are still invalid:
// they are reported whether or not the interval is empty.
void test_invalid_step_options_are_rejected_for_an_empty_interval()
{
    for (const qmomkit::EmbeddedRkMethod method : {qmomkit::EmbeddedRkMethod::rk23, qmomkit::EmbeddedRkMethod::rk45})
    {
        qmomkit::AdaptiveRkOptions negative_initial_step;
        negative_initial_step.initial_step = -1.0;
        const qmomkit::AdaptiveRkResult initial_step_result =
            qmomkit::integrate_embedded_rk({1.0}, 0.5, 0.5, exponential_decay_rhs(), method, negative_initial_step);
        assert(initial_step_result.status == qmomkit::AdaptiveRkStatus::invalid_input);
        assert(initial_step_result.final_time == 0.5);

        qmomkit::AdaptiveRkOptions non_finite_min_step;
        non_finite_min_step.min_step = std::numeric_limits<double>::quiet_NaN();
        assert(qmomkit::integrate_embedded_rk({1.0}, 0.5, 0.5, exponential_decay_rhs(), method, non_finite_min_step).status ==
               qmomkit::AdaptiveRkStatus::invalid_input);

        qmomkit::AdaptiveRkOptions negative_max_step;
        negative_max_step.max_step = -1.0;
        assert(qmomkit::integrate_embedded_rk({1.0}, 0.5, 0.5, exponential_decay_rhs(), method, negative_max_step).status ==
               qmomkit::AdaptiveRkStatus::invalid_input);

        qmomkit::AdaptiveRkOptions min_step_above_max_step;
        min_step_above_max_step.min_step = 0.5;
        min_step_above_max_step.max_step = 0.1;
        assert(qmomkit::integrate_embedded_rk({1.0}, 0.5, 0.5, exponential_decay_rhs(), method, min_step_above_max_step).status ==
               qmomkit::AdaptiveRkStatus::invalid_input);

        // Valid options over an empty interval return the input state unchanged.
        const qmomkit::AdaptiveRkResult empty = qmomkit::integrate_embedded_rk({1.0}, 0.5, 0.5, exponential_decay_rhs(), method, {});
        assert(empty.status == qmomkit::AdaptiveRkStatus::success);
        assert(empty.state == std::vector<double>{1.0});
        assert(empty.final_time == 0.5);
        assert(empty.function_evaluations == 0);
    }
}

void test_rk23_rejects_min_step_exceeding_max_step()
{
    for (const qmomkit::EmbeddedRkMethod method : {qmomkit::EmbeddedRkMethod::rk23, qmomkit::EmbeddedRkMethod::rk45})
    {
        qmomkit::AdaptiveRkOptions invalid_options;
        invalid_options.min_step = 5.0;
        invalid_options.max_step = 2.0;
        const qmomkit::AdaptiveRkResult invalid_result =
            qmomkit::integrate_embedded_rk({1.0}, 0.0, 100.0, exponential_decay_rhs(), method, invalid_options);
        assert(invalid_result.status == qmomkit::AdaptiveRkStatus::invalid_input);
        assert(!invalid_result.message.empty());

        // min_step == max_step (a fixed effective step) remains valid.
        qmomkit::AdaptiveRkOptions boundary_options;
        boundary_options.min_step = 2.0;
        boundary_options.max_step = 2.0;
        boundary_options.initial_step = 2.0;
        const qmomkit::AdaptiveRkResult boundary_result =
            qmomkit::integrate_embedded_rk({1.0}, 0.0, 4.0, exponential_decay_rhs(), method, boundary_options);
        assert(boundary_result.status == qmomkit::AdaptiveRkStatus::success);

        // An auto-derived max_step (options.max_step == 0.0, defaulting to the
        // integration interval length) must also be checked against an explicit
        // min_step, not just an explicitly configured max_step.
        qmomkit::AdaptiveRkOptions auto_max_step_options;
        auto_max_step_options.min_step = 5.0;
        const qmomkit::AdaptiveRkResult auto_max_step_result =
            qmomkit::integrate_embedded_rk({1.0}, 0.0, 2.0, exponential_decay_rhs(), method, auto_max_step_options);
        assert(auto_max_step_result.status == qmomkit::AdaptiveRkStatus::invalid_input);
    }
}

void test_rk23_rejects_nonfinite_trial_state()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&, std::vector<double>& derivative)
    { derivative[0] = std::numeric_limits<double>::max(); };
    qmomkit::AdaptiveRkOptions options;
    options.initial_step = 2.0;
    options.max_step = 2.0;
    options.min_step = 0.5;

    const qmomkit::AdaptiveRkResult result = qmomkit::integrate_embedded_rk({0.0}, 0.0, 2.0, rhs, qmomkit::EmbeddedRkMethod::rk23, options);

    assert(result.status == qmomkit::AdaptiveRkStatus::step_size_underflow);
    assert(result.rejected_steps > 0);
    assert(std::isfinite(result.state[0]));
}

void test_rhs_that_resizes_the_derivative_buffer_is_rejected()
{
    // The integrator sizes the buffer and the callback fills it. A callback that
    // assigns a differently sized vector instead is the fill-style spelling of
    // the old "returned the wrong size", and must fail the same way.
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&, std::vector<double>& derivative) { derivative.assign(3, 1.0); };

    // The message matters as much as the throw. Downstream code would also fail
    // on the mis-sized vector, but it would blame a vector-size mismatch or the
    // integrator's own cached derivative rather than the callback that resized
    // the buffer, so both drivers are pinned to the accurate message here.
    const std::string expected = "RHS returned a vector with the wrong size.";

    std::string ssp_message;
    try
    {
        (void)qmomkit::ssp_rk3_step({1.0}, 0.0, 0.1, rhs, {});
    }
    catch (const std::invalid_argument& error)
    {
        ssp_message = error.what();
    }
    assert(ssp_message == expected);

    for (const qmomkit::EmbeddedRkMethod method : {qmomkit::EmbeddedRkMethod::rk23, qmomkit::EmbeddedRkMethod::rk45})
    {
        std::string embedded_message;
        try
        {
            (void)qmomkit::integrate_embedded_rk({1.0}, 0.0, 0.1, rhs, method, {});
        }
        catch (const std::invalid_argument& error)
        {
            embedded_message = error.what();
        }
        assert(embedded_message == expected);
    }
}

void test_rhs_that_leaves_a_component_unwritten_is_rejected_in_debug_builds()
{
    // A callback that writes only part of the buffer would otherwise integrate
    // whatever the previous stage left there, and neither the size nor the
    // finiteness of the buffer can reveal it. Debug builds mark the buffer so
    // that an unwritten entry is reported as the missing write it is; release
    // builds cannot detect it and are not asked to.
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&, std::vector<double>& derivative) { derivative[0] = 1.0; };

    bool threw = false;
    try
    {
        (void)qmomkit::integrate_embedded_rk({1.0, 2.0}, 0.0, 0.1, rhs, qmomkit::EmbeddedRkMethod::rk23, {});
    }
    catch (const std::invalid_argument&)
    {
        threw = true;
    }
#ifndef NDEBUG
    assert(threw);
#else
    // A release build cannot detect the missing write and is not asked to. What
    // it does with the unwritten component is deliberately not pinned here: that
    // is the consequence of a contract violation, not a promised behaviour.
    assert(!threw);
#endif
}

void test_rk23_initial_request_below_minimum_uses_minimum_step()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&, std::vector<double>& derivative) { derivative[0] = 1.0; };
    qmomkit::AdaptiveRkOptions options;
    options.initial_step = 1.0e-3;
    options.max_step = 0.05;
    options.min_step = 0.05;

    const qmomkit::AdaptiveRkResult result = qmomkit::integrate_embedded_rk({1.0}, 0.0, 0.05, rhs, qmomkit::EmbeddedRkMethod::rk23, options);

    assert(result.status == qmomkit::AdaptiveRkStatus::success);
    assert(result.accepted_steps == 1);
    assert(result.rejected_steps == 0);
    assert_near(result.final_time, 0.05, 1.0e-15);
    assert_near(result.last_step, 0.05, 1.0e-15);
    assert_near(result.state[0], 1.05, 1.0e-15);
}

// At t = 1e6 the spacing of the double is about 1e-10, so an explicit
// min_step of 1e-12 or 1e-20 is a fraction of it. Such a step can round to no
// step at all, or back to the one just rejected, so the minimum step must not
// fall below the automatic one there.
void test_embedded_rk_min_step_below_time_resolution_terminates()
{
    constexpr double t0 = 1.0e6;
    for (const qmomkit::EmbeddedRkMethod method : {qmomkit::EmbeddedRkMethod::rk23, qmomkit::EmbeddedRkMethod::rk45})
    {
        // dx/dt = -0.5/x drives x from 0.5 to zero at t0 + 0.25, where no step
        // meets the tolerance; the integration must still return.
        std::size_t evaluations = 0;
        const qmomkit::OdeRhs singular = [&evaluations](double, const std::vector<double>& y, std::vector<double>& derivative)
        {
            if (++evaluations > 100000)
            {
                throw std::runtime_error("The integration did not terminate.");
            }
            derivative[0] = -0.5 / y[0];
        };
        qmomkit::AdaptiveRkOptions options;
        options.min_step = 1.0e-12;
        options.max_steps = 1000;
        const qmomkit::AdaptiveRkResult stalled = qmomkit::integrate_embedded_rk({0.5}, t0, t0 + 1.0, singular, method, options);
        assert(stalled.status != qmomkit::AdaptiveRkStatus::success);

        // A smooth problem started at a step far below the spacing must make
        // progress rather than accept steps of zero length.
        options.min_step = 1.0e-20;
        options.initial_step = 1.0e-20;
        const qmomkit::AdaptiveRkResult smooth = qmomkit::integrate_embedded_rk({1.0}, t0, t0 + 1.0, exponential_decay_rhs(), method, options);
        assert(smooth.status == qmomkit::AdaptiveRkStatus::success);
        assert(smooth.final_time == t0 + 1.0);
    }
}

void test_rk23_accepts_finite_update_at_minimum_step()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>& y, std::vector<double>& derivative) { derivative[0] = 1000.0 * y[0]; };
    qmomkit::AdaptiveRkOptions options;
    options.initial_step = 1.0e-3;
    options.max_step = 1.0e-3;
    options.min_step = 1.0e-3;
    options.relative_tolerance = 1.0e-15;
    options.absolute_tolerance = 1.0e-15;

    const qmomkit::AdaptiveRkResult result = qmomkit::integrate_embedded_rk({1.0}, 0.0, 1.0e-3, rhs, qmomkit::EmbeddedRkMethod::rk23, options);

    assert(result.status == qmomkit::AdaptiveRkStatus::success);
    assert(result.accepted_steps == 1);
    assert(result.rejected_steps == 0);
    assert_near(result.final_time, 1.0e-3, 1.0e-18);
    assert_near(result.last_step, 1.0e-3, 1.0e-18);
    assert_near(result.state[0], 8.0 / 3.0, 1.0e-14);
}

void test_rk23_rejected_request_clamps_to_minimum_step()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>& y, std::vector<double>& derivative)
    {
        derivative[0] = y[0];
        derivative[1] = y[1];
        derivative[2] = y[2];
    };
    qmomkit::AdaptiveRkOptions options;
    options.initial_step = 0.04;
    options.max_step = 0.04;
    options.min_step = 0.02;
    options.relative_tolerance = 1.0e-12;
    options.absolute_tolerance = 1.0e-14;

    const qmomkit::AdaptiveRkResult result =
        qmomkit::integrate_embedded_rk({0.5, 1.0, 2.0}, 0.0, 0.04, rhs, qmomkit::EmbeddedRkMethod::rk23, options);

    assert(result.status == qmomkit::AdaptiveRkStatus::success);
    assert(result.accepted_steps == 2);
    assert(result.rejected_steps > 0);
    assert_near(result.final_time, 0.04, 1.0e-15);
    assert_near(result.last_step, 0.02, 1.0e-15);
}

void test_rk23_accepts_final_remainder_below_minimum_step()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>& y, std::vector<double>& derivative) { derivative[0] = 1000.0 * y[0]; };
    qmomkit::AdaptiveRkOptions options;
    options.initial_step = 1.0e-3;
    options.max_step = 1.0e-3;
    options.min_step = 1.0e-3;
    options.relative_tolerance = 1.0e-15;
    options.absolute_tolerance = 1.0e-15;

    const qmomkit::AdaptiveRkResult result = qmomkit::integrate_embedded_rk({1.0}, 0.0, 1.5e-3, rhs, qmomkit::EmbeddedRkMethod::rk23, options);

    assert(result.status == qmomkit::AdaptiveRkStatus::success);
    assert(result.accepted_steps == 2);
    assert(result.rejected_steps == 0);
    assert_near(result.final_time, 1.5e-3, 1.0e-18);
    assert_near(result.last_step, 5.0e-4, 1.0e-18);
    assert_near(result.state[0], (8.0 / 3.0) * (79.0 / 48.0), 1.0e-14);
}

void test_rk23_accepts_complete_interval_below_minimum_step()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>& y, std::vector<double>& derivative) { derivative[0] = 1000.0 * y[0]; };
    qmomkit::AdaptiveRkOptions options;
    options.initial_step = 1.0e-3;
    options.max_step = 1.0e-3;
    options.min_step = 1.0e-3;
    options.relative_tolerance = 1.0e-15;
    options.absolute_tolerance = 1.0e-15;

    const qmomkit::AdaptiveRkResult result = qmomkit::integrate_embedded_rk({1.0}, 0.0, 5.0e-4, rhs, qmomkit::EmbeddedRkMethod::rk23, options);

    assert(result.status == qmomkit::AdaptiveRkStatus::success);
    assert(result.accepted_steps == 1);
    assert(result.rejected_steps == 0);
    assert_near(result.final_time, 5.0e-4, 1.0e-18);
    assert_near(result.last_step, 5.0e-4, 1.0e-18);
    assert_near(result.state[0], 79.0 / 48.0, 1.0e-14);
}

void test_rk23_reports_maximum_step_count_failure()
{
    qmomkit::AdaptiveRkOptions options;
    options.initial_step = 0.1;
    options.max_step = 0.1;
    options.max_steps = 1;

    const qmomkit::AdaptiveRkResult result =
        qmomkit::integrate_embedded_rk({1.0}, 0.0, 1.0, exponential_decay_rhs(), qmomkit::EmbeddedRkMethod::rk23, options);

    assert(result.status == qmomkit::AdaptiveRkStatus::max_steps_exceeded);
    assert(!result.message.empty());
    assert(result.accepted_steps == 1);
    assert(result.final_time > 0.0);
    assert(result.final_time < 1.0);
    assert(std::isfinite(result.state[0]));
}

void test_rk23_max_step_changes_observable_diagnostics()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&, std::vector<double>& derivative)
    { std::fill(derivative.begin(), derivative.end(), 0.0); };

    qmomkit::AdaptiveRkOptions large_step;
    large_step.initial_step = 0.5;
    large_step.max_step = 0.5;

    qmomkit::AdaptiveRkOptions small_step;
    small_step.initial_step = 0.1;
    small_step.max_step = 0.1;

    const qmomkit::AdaptiveRkResult large = qmomkit::integrate_embedded_rk({1.0}, 0.0, 1.0, rhs, qmomkit::EmbeddedRkMethod::rk23, large_step);
    const qmomkit::AdaptiveRkResult small = qmomkit::integrate_embedded_rk({1.0}, 0.0, 1.0, rhs, qmomkit::EmbeddedRkMethod::rk23, small_step);

    assert(large.status == qmomkit::AdaptiveRkStatus::success);
    assert(small.status == qmomkit::AdaptiveRkStatus::success);
    assert(large.accepted_steps == 2);
    assert(small.accepted_steps >= 10);
    assert(large.function_evaluations == 7);
    assert(small.function_evaluations > large.function_evaluations);
    assert_vector_near(large.state, small.state, 1.0e-15);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"embedded_rk_tableaus_match_primary_papers", test_embedded_rk_tableaus_match_primary_papers},
            {"embedded_rk_one_steps_match_independent_scalar_formulae", test_embedded_rk_one_steps_match_independent_scalar_formulae},
            {"embedded_rk_one_steps_use_stage_abscissae_for_time_dependent_rhs",
             test_embedded_rk_one_steps_use_stage_abscissae_for_time_dependent_rhs},
            {"embedded_rk_workspace_reuse_is_transparent", test_embedded_rk_workspace_reuse_is_transparent},
            {"adaptive_initial_step_default_and_override", test_adaptive_initial_step_default_and_override},
            {"embedded_rk_default_initial_step_is_used_for_both_methods", test_embedded_rk_default_initial_step_is_used_for_both_methods},
            {"embedded_rk_fsal_reuses_accepted_endpoint_derivative", test_embedded_rk_fsal_reuses_accepted_endpoint_derivative},
            {"embedded_rk_adaptive_controls_for_both_methods", test_embedded_rk_adaptive_controls_for_both_methods},
            {"embedded_rk_public_backward_integration", test_embedded_rk_public_backward_integration},
            {"ssp_rk3_accept_all_predicate", test_ssp_rk3_accept_all_predicate},
            {"ssp_rk3_reject_negative_values_then_success", test_ssp_rk3_reject_negative_values_then_success},
            {"ssp_rk3_multiple_retries", test_ssp_rk3_multiple_retries},
            {"ssp_rk3_failure_after_max_retries", test_ssp_rk3_failure_after_max_retries},
            {"ssp_rk3_failure_after_reaching_minimum_timestep", test_ssp_rk3_failure_after_reaching_minimum_timestep},
            {"ssp_rk3_rejects_initial_step_below_minimum", test_ssp_rk3_rejects_initial_step_below_minimum},
            {"ssp_rk3_rejects_retry_that_underflows_to_zero", test_ssp_rk3_rejects_retry_that_underflows_to_zero},
            {"ssp_rk3_without_stage_callback_accepts_unchecked_candidate", test_ssp_rk3_without_stage_callback_accepts_unchecked_candidate},
            {"ssp_rk3_final_only_callback_invoked_once", test_ssp_rk3_final_only_callback_invoked_once},
            {"ssp_rk3_stage_callback_invoked_for_each_stage", test_ssp_rk3_stage_callback_invoked_for_each_stage},
            {"ssp_rk3_failed_stage_check_retries_with_smaller_timestep", test_ssp_rk3_failed_stage_check_retries_with_smaller_timestep},
            {"ssp_rk3_stage2_rejection_skips_k3_and_retries_from_original_state",
             test_ssp_rk3_stage2_rejection_skips_k3_and_retries_from_original_state},
            {"ssp_rk3_final_stage_rejection_short_circuits_final_accept", test_ssp_rk3_final_stage_rejection_short_circuits_final_accept},
            {"ssp_rk3_stage_rejection_exhaustion_preserves_original_state", test_ssp_rk3_stage_rejection_exhaustion_preserves_original_state},
            {"ssp_rk3_stage_rejection_minimum_timestep_preserves_original_state",
             test_ssp_rk3_stage_rejection_minimum_timestep_preserves_original_state},
            {"ssp_rk3_rhs_stage_audit_preserves_times_and_retry_states", test_ssp_rk3_rhs_stage_audit_preserves_times_and_retry_states},
            {"rk23_rejects_overly_large_initial_step", test_rk23_rejects_overly_large_initial_step},
            {"rk23_zero_rhs_has_no_timestep_pathology", test_rk23_zero_rhs_has_no_timestep_pathology},
            {"rk23_invalid_options_fail_before_loop", test_rk23_invalid_options_fail_before_loop},
            {"invalid_step_options_are_rejected_for_an_empty_interval", test_invalid_step_options_are_rejected_for_an_empty_interval},
            {"rk23_rejects_min_step_exceeding_max_step", test_rk23_rejects_min_step_exceeding_max_step},
            {"rk23_rejects_nonfinite_trial_state", test_rk23_rejects_nonfinite_trial_state},
            {"rhs_that_resizes_the_derivative_buffer_is_rejected", test_rhs_that_resizes_the_derivative_buffer_is_rejected},
            {"rhs_that_leaves_a_component_unwritten_is_rejected_in_debug_builds",
             test_rhs_that_leaves_a_component_unwritten_is_rejected_in_debug_builds},
            {"rk23_initial_request_below_minimum_uses_minimum_step", test_rk23_initial_request_below_minimum_uses_minimum_step},
            {"embedded_rk_min_step_below_time_resolution_terminates", test_embedded_rk_min_step_below_time_resolution_terminates},
            {"rk23_accepts_finite_update_at_minimum_step", test_rk23_accepts_finite_update_at_minimum_step},
            {"rk23_rejected_request_clamps_to_minimum_step", test_rk23_rejected_request_clamps_to_minimum_step},
            {"rk23_accepts_final_remainder_below_minimum_step", test_rk23_accepts_final_remainder_below_minimum_step},
            {"rk23_accepts_complete_interval_below_minimum_step", test_rk23_accepts_complete_interval_below_minimum_step},
            {"rk23_reports_maximum_step_count_failure", test_rk23_reports_maximum_step_count_failure},
            {"rk23_max_step_changes_observable_diagnostics", test_rk23_max_step_changes_observable_diagnostics},
        });
}
