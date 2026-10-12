// Time integrators against problems with known closed-form solutions.
//
// What is verified
//   That each integrator reproduces the exact trajectory of test problems
//   whose solution can be written down, within the tolerance the integrator
//   was asked to achieve, and that the adaptive controller reports step and
//   evaluation counts consistent with the work it did.
//
// How the reference is constructed
//   The exact solution is evaluated directly from its closed form at the final
//   time, independently of any integrator. Requested controller tolerances are
//   stated per case (relative 1.0e-7 to 1.0e-10, absolute 1.0e-10 to 1.0e-12)
//   and the comparison tolerance is chosen to sit above the error the
//   controller was asked to deliver. A constant right-hand side is integrated
//   exactly, so those cases are compared at roundoff.
//
// Error metric and pass criteria
//   Maximum absolute deviation of each state component from the exact solution
//   at the final time, against the per-case tolerance. A result that merely
//   converged is not enough: it must land within the accuracy requested.
//
// Reproduce
//   ctest -R verification_.*runge_kutta
//
// Build modes
//   Both. Checks go through QMOMKIT_VERIFY_*, not assert().

#include "../verification_helpers.hpp"
#include "qmomkit/runge_kutta.hpp"
#include "test_case_dispatch.hpp"

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
    return [](double, const std::vector<double>& state, std::vector<double>& derivative) { derivative[0] = -state[0]; };
}

void test_ssp_rk3_constant_rhs()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&, std::vector<double>& derivative)
    {
        derivative[0] = 1.0;
        derivative[1] = -2.0;
        derivative[2] = 0.5;
    };
    qmomkit::SspRk3Options options;
    options.estimate_error = true;
    const qmomkit::SspRk3Result result = qmomkit::ssp_rk3_step({0.0, 3.0, -1.0}, 1.5, 0.25, rhs, options);

    assert_vector_near(result.state, {0.25, 2.5, -0.875}, 1.0e-15);
    QMOMKIT_VERIFY_CHECK(result.error_estimate.has_value());
    assert_near(*result.error_estimate, 0.0, 1.0e-15);
}

void test_rk23_exponential_decay()
{
    qmomkit::AdaptiveRkOptions options;
    options.relative_tolerance = 1.0e-7;
    options.absolute_tolerance = 1.0e-10;
    options.initial_step = 0.1;

    const qmomkit::AdaptiveRkResult result =
        qmomkit::integrate_embedded_rk({1.0}, 0.0, 2.0, exponential_decay_rhs(), qmomkit::EmbeddedRkMethod::rk23, options);

    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::AdaptiveRkStatus::success);
    assert_near(result.final_time, 2.0, 1.0e-14);
    assert_near(result.state[0], std::exp(-2.0), 5.0e-7);
    QMOMKIT_VERIFY_CHECK(result.accepted_steps > 1);
    QMOMKIT_VERIFY_CHECK(result.function_evaluations > result.accepted_steps);
}

void test_rk23_constant_rhs()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&, std::vector<double>& derivative)
    {
        derivative[0] = 2.0;
        derivative[1] = -3.0;
    };
    qmomkit::AdaptiveRkOptions options;
    options.initial_step = 0.75;
    options.max_step = 0.75;

    const qmomkit::AdaptiveRkResult result = qmomkit::integrate_embedded_rk({1.0, -2.0}, 0.0, 1.5, rhs, qmomkit::EmbeddedRkMethod::rk23, options);

    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::AdaptiveRkStatus::success);
    assert_vector_near(result.state, {4.0, -6.5}, 1.0e-13);
}

void test_rk23_vector_coupled_linear_system()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>& state, std::vector<double>& derivative)
    {
        derivative[0] = state[1];
        derivative[1] = -state[0];
    };
    qmomkit::AdaptiveRkOptions options;
    options.relative_tolerance = 1.0e-8;
    options.absolute_tolerance = 1.0e-10;
    options.initial_step = 0.1;

    const double final_time = 1.0;
    const qmomkit::AdaptiveRkResult result =
        qmomkit::integrate_embedded_rk({1.0, 0.0}, 0.0, final_time, rhs, qmomkit::EmbeddedRkMethod::rk23, options);

    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::AdaptiveRkStatus::success);
    assert_near(result.state[0], std::cos(final_time), 2.0e-7);
    assert_near(result.state[1], -std::sin(final_time), 2.0e-7);
}

void test_rk45_vector_coupled_linear_system()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>& state, std::vector<double>& derivative)
    {
        derivative[0] = state[1];
        derivative[1] = -state[0];
    };
    qmomkit::AdaptiveRkOptions options;
    options.relative_tolerance = 1.0e-10;
    options.absolute_tolerance = 1.0e-12;
    options.initial_step = 0.25;

    constexpr double final_time = 1.0;
    const qmomkit::AdaptiveRkResult result =
        qmomkit::integrate_embedded_rk({1.0, 0.0}, 0.0, final_time, rhs, qmomkit::EmbeddedRkMethod::rk45, options);

    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::AdaptiveRkStatus::success);
    assert_near(result.final_time, final_time, 1.0e-14);
    assert_near(result.state[0], std::cos(final_time), 2.0e-10);
    assert_near(result.state[1], -std::sin(final_time), 2.0e-10);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(argc, argv,
                                             {
                                                 {"ssp_rk3_constant_rhs", test_ssp_rk3_constant_rhs},
                                                 {"rk23_exponential_decay", test_rk23_exponential_decay},
                                                 {"rk23_constant_rhs", test_rk23_constant_rhs},
                                                 {"rk23_vector_coupled_linear_system", test_rk23_vector_coupled_linear_system},
                                                 {"rk45_vector_coupled_linear_system", test_rk45_vector_coupled_linear_system},
                                             });
}
