#include "internal/moment_source_rates.hpp"

#include "internal/performance.hpp"
#include "qmomkit/qmomkit.hpp"
#include "test_case_dispatch.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>
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

bool same_layout(const qmomkit::CqmomMomentLayout& lhs, const qmomkit::CqmomMomentLayout& rhs)
{
    return lhs.hierarchy.node_counts() == rhs.hierarchy.node_counts() && lhs.hierarchy.conditioning_order() == rhs.hierarchy.conditioning_order() &&
           lhs.required_indices == rhs.required_indices;
}

qmomkit::CqmomMomentData make_cqmom_moment_data(const qmomkit::PhaseSpaceHierarchy& hierarchy)
{
    qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    std::vector<double> values(layout.required_indices.size(), 0.0);
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        const qmomkit::MultiIndex& index = layout.required_indices[i];
        values[i] = 1.0 + static_cast<double>(index[0]) + 10.0 * static_cast<double>(index.size() > 1 ? index[1] : 0);
    }
    return qmomkit::CqmomMomentData{layout, values};
}

void test_qmom_success_ordering_and_single_invocation()
{
    const std::vector<double> moments{1.0, 2.0, 3.0};
    std::size_t calls = 0;
    const qmomkit::MomentSource source = [&calls](const double time, const qmomkit::QmomSourceContext& state)
    {
        ++calls;
        assert_near(time, 0.25, 0.0);
        assert_vector_near(state.values(), {1.0, 2.0, 3.0}, 0.0);
        return std::vector<double>{state.value(2), state.value(0), state.value(1)};
    };

    const qmomkit::QmomMomentSourceRateResult result = qmomkit::internal::evaluate_qmom_moment_source_rates(0.25, moments, source);

    assert(result.ok());
    assert(result.status == qmomkit::MomentSourceRateStatus::success);
    assert_vector_near(result.rates, {3.0, 1.0, 2.0}, 0.0);
    assert(calls == 1);
}

void test_qmom_failure_paths()
{
    const std::vector<double> moments{1.0, 2.0};

    const qmomkit::QmomMomentSourceRateResult missing = qmomkit::internal::evaluate_qmom_moment_source_rates(0.0, moments, {});
    assert(!missing.ok());
    assert(missing.status == qmomkit::MomentSourceRateStatus::invalid_input);

    const qmomkit::QmomMomentSourceRateResult wrong_size = qmomkit::internal::evaluate_qmom_moment_source_rates(
        0.0, moments, [](double, const qmomkit::QmomSourceContext& state) { return std::vector<double>(state.size() + 1, 0.0); });
    assert(!wrong_size.ok());
    assert(wrong_size.status == qmomkit::MomentSourceRateStatus::invalid_rate);
    assert(wrong_size.message == "Moment-source callback rate values must match the moment state size.");

    const qmomkit::QmomMomentSourceRateResult non_finite =
        qmomkit::internal::evaluate_qmom_moment_source_rates(0.0, moments,
                                                             [](double, const qmomkit::QmomSourceContext& state)
                                                             {
                                                                 std::vector<double> rate(state.size(), 0.0);
                                                                 rate.back() = std::numeric_limits<double>::quiet_NaN();
                                                                 return rate;
                                                             });
    assert(!non_finite.ok());
    assert(non_finite.status == qmomkit::MomentSourceRateStatus::invalid_rate);
    assert(non_finite.message == "Moment-source callback rate values must be finite.");

    const qmomkit::QmomMomentSourceRateResult thrown = qmomkit::internal::evaluate_qmom_moment_source_rates(
        0.0, moments, [](double, const qmomkit::QmomSourceContext&) -> std::vector<double> { throw std::runtime_error("source model failed"); });
    assert(!thrown.ok());
    assert(thrown.status == qmomkit::MomentSourceRateStatus::source_failed);
    assert(thrown.message.find("source model failed") != std::string::npos);
}

void test_qmom_non_finite_input_is_rejected_before_callback()
{
    for (const double non_finite :
         {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()})
    {
        std::size_t calls = 0;
        const qmomkit::QmomMomentSourceRateResult result =
            qmomkit::evaluate_qmom_moment_source_rates(0.0, {1.0, non_finite, 3.0},
                                                       [&calls](double, const qmomkit::QmomSourceContext& context)
                                                       {
                                                           ++calls;
                                                           return std::vector<double>(context.size(), 2.0);
                                                       });

        assert(!result.ok());
        assert(result.status == qmomkit::MomentSourceRateStatus::invalid_input);
        assert(result.rates.empty());
        assert(calls == 0);
    }
}

void test_cqmom_identity_context_ordering_and_layout_preservation()
{
    const qmomkit::CqmomMomentData moments = make_cqmom_moment_data(qmomkit::PhaseSpaceHierarchy({2, 2}));
    std::size_t calls = 0;

    const qmomkit::CqmomMomentSource source = [&moments, &calls](const double time, const qmomkit::CqmomSourceContext& context)
    {
        ++calls;
        assert_near(time, 0.5, 0.0);
        assert(context.size() == moments.values.size());
        assert(context.hierarchy().node_counts() == moments.layout.hierarchy.node_counts());
        assert(context.hierarchy().conditioning_order() == moments.layout.hierarchy.conditioning_order());

        std::vector<double> rates(context.size(), 0.0);
        for (std::size_t i = 0; i < context.size(); ++i)
        {
            assert(context.index(i) == moments.layout.required_indices[i]);
            assert_near(context.value(i), moments.values[i], 0.0);
            rates[i] = static_cast<double>(i) + 0.25;
        }
        return qmomkit::CqmomMomentRate{rates};
    };

    const qmomkit::CqmomMomentSourceRateResult result = qmomkit::internal::evaluate_cqmom_moment_source_rates(0.5, moments, source);

    assert(result.ok());
    assert(result.rates.has_value());
    assert(same_layout(result.rates->layout, moments.layout));
    for (std::size_t i = 0; i < result.rates->values.size(); ++i)
    {
        assert_near(result.rates->values[i], static_cast<double>(i) + 0.25, 0.0);
    }
    assert(calls == 1);
}

void test_cqmom_non_identity_ordering()
{
    const qmomkit::CqmomMomentData moments = make_cqmom_moment_data(qmomkit::PhaseSpaceHierarchy({2, 2}, {1, 0}));
    const qmomkit::CqmomMomentSource source = [&moments](double, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> rates(context.size(), 0.0);
        for (std::size_t i = 0; i < context.size(); ++i)
        {
            assert(context.index(i) == moments.layout.required_indices[i]);
            rates[i] = 100.0 + static_cast<double>(context.index(i)[0]) + 10.0 * static_cast<double>(context.index(i)[1]);
        }
        return qmomkit::CqmomMomentRate{rates};
    };

    const qmomkit::CqmomMomentSourceRateResult result = qmomkit::internal::evaluate_cqmom_moment_source_rates(0.0, moments, source);

    assert(result.ok());
    assert(result.rates.has_value());
    assert(same_layout(result.rates->layout, moments.layout));
    for (std::size_t i = 0; i < moments.layout.required_indices.size(); ++i)
    {
        const qmomkit::MultiIndex& index = moments.layout.required_indices[i];
        assert_near(result.rates->values[i], 100.0 + static_cast<double>(index[0]) + 10.0 * static_cast<double>(index[1]), 0.0);
    }
}

void test_cqmom_failure_paths()
{
    const qmomkit::CqmomMomentData moments = make_cqmom_moment_data(qmomkit::PhaseSpaceHierarchy({2, 2}));

    const qmomkit::CqmomMomentSourceRateResult missing = qmomkit::internal::evaluate_cqmom_moment_source_rates(0.0, moments, {});
    assert(!missing.ok());
    assert(missing.status == qmomkit::MomentSourceRateStatus::invalid_input);
    assert(!missing.rates.has_value());

    qmomkit::CqmomMomentData malformed = moments;
    malformed.values.pop_back();
    const qmomkit::CqmomMomentSourceRateResult malformed_result =
        qmomkit::internal::evaluate_cqmom_moment_source_rates(0.0, malformed, [](double, const qmomkit::CqmomSourceContext& context)
                                                              { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; });
    assert(!malformed_result.ok());
    assert(malformed_result.status == qmomkit::MomentSourceRateStatus::invalid_input);
    assert(malformed_result.message == "CQMOM moment-source stage values must match the fixed layout size.");

    const qmomkit::CqmomMomentSourceRateResult wrong_size =
        qmomkit::internal::evaluate_cqmom_moment_source_rates(0.0, moments, [](double, const qmomkit::CqmomSourceContext& context)
                                                              { return qmomkit::CqmomMomentRate{std::vector<double>(context.size() + 1, 0.0)}; });
    assert(!wrong_size.ok());
    assert(wrong_size.status == qmomkit::MomentSourceRateStatus::invalid_rate);
    assert(wrong_size.message == "CQMOM moment-source rate values must match the source context size.");

    const qmomkit::CqmomMomentSourceRateResult non_finite =
        qmomkit::internal::evaluate_cqmom_moment_source_rates(0.0, moments,
                                                              [](double, const qmomkit::CqmomSourceContext& context)
                                                              {
                                                                  std::vector<double> rate(context.size(), 0.0);
                                                                  rate.front() = std::numeric_limits<double>::infinity();
                                                                  return qmomkit::CqmomMomentRate{rate};
                                                              });
    assert(!non_finite.ok());
    assert(non_finite.status == qmomkit::MomentSourceRateStatus::invalid_rate);
    assert(non_finite.message == "CQMOM moment-source rate values must be finite.");

    const qmomkit::CqmomMomentSourceRateResult thrown =
        qmomkit::internal::evaluate_cqmom_moment_source_rates(0.0, moments, [](double, const qmomkit::CqmomSourceContext&) -> qmomkit::CqmomMomentRate
                                                              { throw std::runtime_error("CQMOM source model failed"); });
    assert(!thrown.ok());
    assert(thrown.status == qmomkit::MomentSourceRateStatus::source_failed);
    assert(thrown.message.find("CQMOM source model failed") != std::string::npos);
}

void test_direct_evaluators_do_not_integrate_retry_or_reconstruct()
{
    const std::vector<double> qmom_moments{1.0, 2.0};
    const qmomkit::CqmomMomentData cqmom_moments = make_cqmom_moment_data(qmomkit::PhaseSpaceHierarchy({2, 2}));
    bool acceptability_callback_invoked = false;
    std::size_t qmom_calls = 0;
    std::size_t cqmom_calls = 0;

    qmomkit::internal::performance::CounterSession session;
    {
        qmomkit::internal::performance::CounterScope scope(session);
        const qmomkit::QmomMomentSourceRateResult qmom_result =
            qmomkit::internal::evaluate_qmom_moment_source_rates(0.0, qmom_moments,
                                                                 [&qmom_calls](double, const qmomkit::QmomSourceContext& moments)
                                                                 {
                                                                     ++qmom_calls;
                                                                     return std::vector<double>(moments.size(), 0.0);
                                                                 });
        const qmomkit::CqmomMomentSourceRateResult cqmom_result =
            qmomkit::internal::evaluate_cqmom_moment_source_rates(0.0, cqmom_moments,
                                                                  [&cqmom_calls](double, const qmomkit::CqmomSourceContext& context)
                                                                  {
                                                                      ++cqmom_calls;
                                                                      return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)};
                                                                  });
        assert(qmom_result.ok());
        assert(cqmom_result.ok());
    }

    const qmomkit::internal::performance::PerformanceCounters counters = session.snapshot();
    assert(counters.time_integration.source_rk_stages == 0);
    assert(counters.time_integration.source_rk_retries == 0);
    assert(counters.reconstruction.qmom_reconstructions == 0);
    assert(counters.reconstruction.cqmom_reconstructions == 0);
    assert(counters.reconstruction.wheeler_calls == 0);
    assert(counters.reconstruction.moment_layout_builds == 0);
    assert(qmom_calls == 1);
    assert(cqmom_calls == 1);
    assert(!acceptability_callback_invoked);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"qmom_success_ordering_and_single_invocation", test_qmom_success_ordering_and_single_invocation},
            {"qmom_failure_paths", test_qmom_failure_paths},
            {"qmom_non_finite_input_is_rejected_before_callback", test_qmom_non_finite_input_is_rejected_before_callback},
            {"cqmom_identity_context_ordering_and_layout_preservation", test_cqmom_identity_context_ordering_and_layout_preservation},
            {"cqmom_non_identity_ordering", test_cqmom_non_identity_ordering},
            {"cqmom_failure_paths", test_cqmom_failure_paths},
            {"direct_evaluators_do_not_integrate_retry_or_reconstruct", test_direct_evaluators_do_not_integrate_retry_or_reconstruct},
        });
}
