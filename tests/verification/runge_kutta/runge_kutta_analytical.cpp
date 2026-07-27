#include "../verification_helpers.hpp"
#include "qmomkit/runge_kutta.hpp"

#include <cmath>
#include <vector>

namespace
{

void assert_near(const double actual, const double expected, const double tolerance)
{
    QMOMKIT_VERIFY_NEAR(actual, expected, tolerance);
}

void assert_vector_near(const std::vector<double>& actual, const std::vector<double>& expected, const double tolerance)
{
    QMOMKIT_VERIFY_CHECK(actual.size() == expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        assert_near(actual[i], expected[i], tolerance);
    }
}

qmomkit::OdeRhs exponential_decay_rhs()
{
    return [](double, const std::vector<double>& state) { return std::vector<double>{-state[0]}; };
}

void test_ssp_rk3_constant_rhs()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&) { return std::vector<double>{1.0, -2.0, 0.5}; };
    const qmomkit::SspRk3Result result = qmomkit::ssp_rk3_step({0.0, 3.0, -1.0}, 1.5, 0.25, rhs, true);

    assert_vector_near(result.state, {0.25, 2.5, -0.875}, 1.0e-15);
    QMOMKIT_VERIFY_CHECK(result.error_estimate.has_value());
    assert_near(*result.error_estimate, 0.0, 1.0e-15);
}

void test_rk23_exponential_decay()
{
    qmomkit::AdaptiveRk23Options options;
    options.relative_tolerance = 1.0e-7;
    options.absolute_tolerance = 1.0e-10;

    const qmomkit::AdaptiveRk23Result result = qmomkit::integrate_rk23({1.0}, 0.0, 2.0, exponential_decay_rhs(), options);

    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::AdaptiveRk23Status::success);
    assert_near(result.t, 2.0, 1.0e-14);
    assert_near(result.state[0], std::exp(-2.0), 5.0e-7);
    QMOMKIT_VERIFY_CHECK(result.accepted_steps > 1);
    QMOMKIT_VERIFY_CHECK(result.function_evaluations > result.accepted_steps);
}

void test_rk23_constant_rhs()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&) { return std::vector<double>{2.0, -3.0}; };
    qmomkit::AdaptiveRk23Options options;
    options.initial_step = 0.75;
    options.max_step = 0.75;

    const qmomkit::AdaptiveRk23Result result = qmomkit::integrate_rk23({1.0, -2.0}, 0.0, 1.5, rhs, options);

    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::AdaptiveRk23Status::success);
    assert_vector_near(result.state, {4.0, -6.5}, 1.0e-13);
}

void test_rk23_vector_coupled_linear_system()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>& state) { return std::vector<double>{state[1], -state[0]}; };
    qmomkit::AdaptiveRk23Options options;
    options.relative_tolerance = 1.0e-8;
    options.absolute_tolerance = 1.0e-10;

    const double final_time = 1.0;
    const qmomkit::AdaptiveRk23Result result = qmomkit::integrate_rk23({1.0, 0.0}, 0.0, final_time, rhs, options);

    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::AdaptiveRk23Status::success);
    assert_near(result.state[0], std::cos(final_time), 2.0e-7);
    assert_near(result.state[1], -std::sin(final_time), 2.0e-7);
}

} // namespace

int main()
{
    test_ssp_rk3_constant_rhs();
    test_rk23_exponential_decay();
    test_rk23_constant_rhs();
    test_rk23_vector_coupled_linear_system();
    return 0;
}
