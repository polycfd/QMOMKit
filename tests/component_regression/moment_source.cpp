#include "../verification/verification_helpers.hpp"
#include "qmomkit/qmomkit.hpp"

#include <cmath>
#include <cstddef>
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

double rk3_linear_factor(const double lambda, const double time_step)
{
    const double z = lambda * time_step;
    return 1.0 + z + 0.5 * z * z + (z * z * z) / 6.0;
}

void test_constant_source_term()
{
    const std::vector<double> initial{0.0, 3.0, -1.0};
    const std::vector<double> source_values{1.0, -2.0, 0.5};
    constexpr double time_step = 0.25;

    const qmomkit::MomentSourceUpdateResult result =
        qmomkit::advance_moments_with_source(initial, 1.5, time_step,
                                             [source_values](double, const qmomkit::QmomSourceContext& context)
                                             {
                                                 QMOMKIT_VERIFY_CHECK(context.size() == source_values.size());
                                                 return source_values;
                                             });

    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::MomentSourceUpdateStatus::success);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.ssp_rk3_status == qmomkit::SspRk3Status::success);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.accepted);
    QMOMKIT_VERIFY_CHECK(!result.diagnostics.rejected);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.retries == 0);
    assert_near(result.diagnostics.time_step_used, time_step, 1.0e-15);
    assert_vector_near(result.moments, {0.25, 2.5, -0.875}, 1.0e-15);
}

void test_linear_source_term()
{
    const std::vector<double> initial{2.0, -1.0};
    constexpr double lambda = -0.75;
    constexpr double time_step = 0.4;

    const qmomkit::MomentSourceUpdateResult result = qmomkit::advance_moments_with_source(initial, 0.0, time_step,
                                                                                          [](double, const qmomkit::QmomSourceContext& context)
                                                                                          {
                                                                                              std::vector<double> derivative(context.size(), 0.0);
                                                                                              for (std::size_t i = 0; i < context.size(); ++i)
                                                                                              {
                                                                                                  derivative[i] = lambda * context.value(i);
                                                                                              }
                                                                                              return derivative;
                                                                                          });

    const double factor = rk3_linear_factor(lambda, time_step);
    QMOMKIT_VERIFY_CHECK(result.ok());
    assert_vector_near(result.moments, {initial[0] * factor, initial[1] * factor}, 1.0e-15);
}

void test_vector_valued_coupled_source_term()
{
    const std::vector<double> initial{1.0, 0.0};
    constexpr double time_step = 0.2;

    const qmomkit::MomentSourceUpdateResult result =
        qmomkit::advance_moments_with_source(initial, 0.0, time_step, [](double, const qmomkit::QmomSourceContext& context)
                                             { return std::vector<double>{context.value(1), -context.value(0)}; });

    QMOMKIT_VERIFY_CHECK(result.ok());
    assert_vector_near(result.moments, {1.0 - 0.5 * time_step * time_step, -time_step + (time_step * time_step * time_step) / 6.0}, 1.0e-15);
}

} // namespace

int main()
{
    test_constant_source_term();
    test_linear_source_term();
    test_vector_valued_coupled_source_term();
    return 0;
}
