#include "qmomkit/runge_kutta.hpp"

#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
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

qmomkit::OdeRhs exponential_decay_rhs()
{
    return [](double, const std::vector<double>& y) { return std::vector<double>{-y[0]}; };
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
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&) { return std::vector<double>{-2.0}; };
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
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&) { return std::vector<double>{-8.0}; };
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
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&) { return std::vector<double>{-8.0}; };
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
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&) { return std::vector<double>{-8.0}; };
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

void test_ssp_rk3_without_stage_callback_accepts_unchecked_candidate()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&) { return std::vector<double>{0.0, 0.0, -2.0}; };

    const qmomkit::SspRk3Result result = qmomkit::ssp_rk3_step({1.0, 0.0, 1.0}, 0.0, 1.0, rhs);

    assert(result.ok());
    assert(result.retries == 0);
    assert_vector_near(result.state, {1.0, 0.0, -1.0}, 1.0e-15);
}

void test_ssp_rk3_final_only_callback_invoked_once()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&) { return std::vector<double>{0.0, 0.0, -2.0}; };
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
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&) { return std::vector<double>{1.0}; };
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
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&) { return std::vector<double>{-2.0}; };
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

void test_rk23_rejects_overly_large_initial_step()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>& y) { return std::vector<double>{-10.0 * y[0]}; };
    qmomkit::AdaptiveRk23Options options;
    options.initial_step = 1.0;
    options.max_step = 1.0;
    options.relative_tolerance = 1.0e-9;
    options.absolute_tolerance = 1.0e-12;

    const qmomkit::AdaptiveRk23Result result = qmomkit::integrate_rk23({1.0}, 0.0, 1.0, rhs, options);

    assert(result.status == qmomkit::AdaptiveRk23Status::success);
    assert(result.rejected_steps > 0);
    assert_near(result.state[0], std::exp(-10.0), 5.0e-9);
}

void test_rk23_zero_rhs_has_no_timestep_pathology()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>& y) { return std::vector<double>(y.size(), 0.0); };
    qmomkit::AdaptiveRk23Options options;
    options.initial_step = 0.1;
    options.relative_tolerance = 1.0e-8;
    options.absolute_tolerance = 1.0e-12;

    const qmomkit::AdaptiveRk23Result result = qmomkit::integrate_rk23({0.0, 2.0, -3.0}, 0.0, 1.0, rhs, options);

    assert(result.status == qmomkit::AdaptiveRk23Status::success);
    assert_near(result.t, 1.0, 1.0e-14);
    assert_vector_near(result.state, {0.0, 2.0, -3.0}, 1.0e-15);
    assert(result.rejected_steps == 0);
}

void test_rk23_invalid_options_fail_before_loop()
{
    qmomkit::AdaptiveRk23Options options;
    options.absolute_tolerance = 0.0;

    const qmomkit::AdaptiveRk23Result zero_atol = qmomkit::integrate_rk23({0.0}, 0.0, 1.0, exponential_decay_rhs(), options);
    assert(zero_atol.status == qmomkit::AdaptiveRk23Status::invalid_input);

    options.absolute_tolerance = 1.0e-9;
    options.initial_step = -1.0;
    const qmomkit::AdaptiveRk23Result negative_initial_step = qmomkit::integrate_rk23({1.0}, 0.0, 1.0, exponential_decay_rhs(), options);
    assert(negative_initial_step.status == qmomkit::AdaptiveRk23Status::invalid_input);
}

void test_rk23_rejects_nonfinite_trial_state()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&) { return std::vector<double>{std::numeric_limits<double>::max()}; };
    qmomkit::AdaptiveRk23Options options;
    options.initial_step = 2.0;
    options.max_step = 2.0;
    options.min_step = 0.5;

    const qmomkit::AdaptiveRk23Result result = qmomkit::integrate_rk23({0.0}, 0.0, 2.0, rhs, options);

    assert(result.status == qmomkit::AdaptiveRk23Status::step_size_underflow);
    assert(result.rejected_steps > 0);
    assert(std::isfinite(result.state[0]));
}

void test_rk23_reports_maximum_step_count_failure()
{
    qmomkit::AdaptiveRk23Options options;
    options.initial_step = 0.1;
    options.max_step = 0.1;
    options.max_steps = 1;

    const qmomkit::AdaptiveRk23Result result = qmomkit::integrate_rk23({1.0}, 0.0, 1.0, exponential_decay_rhs(), options);

    assert(result.status == qmomkit::AdaptiveRk23Status::max_steps_exceeded);
    assert(!result.message.empty());
    assert(result.accepted_steps == 1);
    assert(result.t > 0.0);
    assert(result.t < 1.0);
    assert(std::isfinite(result.state[0]));
}

} // namespace

int main()
{
    test_ssp_rk3_accept_all_predicate();
    test_ssp_rk3_reject_negative_values_then_success();
    test_ssp_rk3_multiple_retries();
    test_ssp_rk3_failure_after_max_retries();
    test_ssp_rk3_failure_after_reaching_minimum_timestep();
    test_ssp_rk3_rejects_initial_step_below_minimum();
    test_ssp_rk3_without_stage_callback_accepts_unchecked_candidate();
    test_ssp_rk3_final_only_callback_invoked_once();
    test_ssp_rk3_stage_callback_invoked_for_each_stage();
    test_ssp_rk3_failed_stage_check_retries_with_smaller_timestep();
    test_rk23_rejects_overly_large_initial_step();
    test_rk23_zero_rhs_has_no_timestep_pathology();
    test_rk23_invalid_options_fail_before_loop();
    test_rk23_rejects_nonfinite_trial_state();
    test_rk23_reports_maximum_step_count_failure();

    return 0;
}
