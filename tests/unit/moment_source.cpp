#include "qmomkit/qmomkit.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
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

double rk3_linear_factor(const double lambda, const double time_step)
{
    const double z = lambda * time_step;
    return 1.0 + z + 0.5 * z * z + (z * z * z) / 6.0;
}

bool same_layout(const qmomkit::CqmomMomentLayout& lhs, const qmomkit::CqmomMomentLayout& rhs)
{
    return lhs.hierarchy.node_counts() == rhs.hierarchy.node_counts() && lhs.hierarchy.conditioning_order() == rhs.hierarchy.conditioning_order() &&
           lhs.required_indices == rhs.required_indices;
}

qmomkit::CqmomMomentData create_moment_data()
{
    qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(qmomkit::PhaseSpaceHierarchy{{2, 1}, {1, 0}});
    std::vector<double> values(layout.required_indices.size(), 0.0);
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        values[i] = 1.0 + 0.25 * static_cast<double>(i);
    }
    return qmomkit::CqmomMomentData{layout, values};
}

void test_accept_all_predicate()
{
    qmomkit::MomentSourceUpdateOptions options;
    options.ssp_rk3_options.accept = [](const std::vector<double>&) { return true; };
    options.ssp_rk3_options.max_retries = 4;

    const qmomkit::MomentSourceUpdateResult result =
        qmomkit::advance_moments_with_source({1.0}, 0.0, 0.2, [](double, const std::vector<double>&) { return std::vector<double>{3.0}; }, options);

    assert(result.ok());
    assert(result.diagnostics.accepted);
    assert(!result.diagnostics.rejected);
    assert(result.diagnostics.retries == 0);
    assert_near(result.moments[0], 1.6, 1.0e-15);
}

void test_reject_negative_moments_predicate_accepts_nonnegative_update()
{
    qmomkit::MomentSourceUpdateOptions options;
    options.ssp_rk3_options.accept = [](const std::vector<double>& state)
    {
        for (const double value : state)
        {
            if (value < 0.0)
            {
                return false;
            }
        }
        return true;
    };
    options.ssp_rk3_options.max_retries = 2;

    const qmomkit::MomentSourceUpdateResult result = qmomkit::advance_moments_with_source(
        {0.3, 1.0}, 0.0, 0.2, [](double, const std::vector<double>&) { return std::vector<double>{-0.5, -1.0}; }, options);

    assert(result.ok());
    assert(result.diagnostics.accepted);
    assert(result.diagnostics.retries == 0);
    assert_vector_near(result.moments, {0.2, 0.8}, 1.0e-15);
}

void test_retry_behavior_first_step_rejected_then_halved_step_succeeds()
{
    qmomkit::MomentSourceUpdateOptions options;
    options.ssp_rk3_options.accept = [](const std::vector<double>& state) { return state[0] >= 0.0; };
    options.ssp_rk3_options.max_retries = 2;

    const qmomkit::MomentSourceUpdateResult result =
        qmomkit::advance_moments_with_source({1.0}, 0.0, 1.0, [](double, const std::vector<double>&) { return std::vector<double>{-2.0}; }, options);

    assert(result.ok());
    assert(result.diagnostics.accepted);
    assert(result.diagnostics.rejected);
    assert(result.diagnostics.retries == 1);
    assert_near(result.diagnostics.time_step_used, 0.5, 1.0e-15);
    assert_near(result.moments[0], 0.0, 1.0e-15);
}

void test_failure_behavior_predicate_rejects_all_attempts_and_preserves_original_state()
{
    const std::vector<double> initial{1.0, -2.0};

    qmomkit::MomentSourceUpdateOptions options;
    options.ssp_rk3_options.accept = [](const std::vector<double>&) { return false; };
    options.ssp_rk3_options.max_retries = 2;

    const qmomkit::MomentSourceUpdateResult result = qmomkit::advance_moments_with_source(
        initial, 0.0, 1.0, [](double, const std::vector<double>& moments) { return std::vector<double>(moments.size(), 0.0); }, options);

    assert(!result.ok());
    assert(result.status == qmomkit::MomentSourceUpdateStatus::rejected_by_predicate);
    assert(result.diagnostics.ssp_rk3_status == qmomkit::SspRk3Status::rejected_by_predicate);
    assert(!result.diagnostics.accepted);
    assert(result.diagnostics.rejected);
    assert(result.diagnostics.retries == 2);
    assert_near(result.diagnostics.time_step_used, 0.25, 1.0e-15);
    assert(result.moments == initial);
}

void test_failure_behavior_minimum_timestep_preserves_original_state()
{
    const std::vector<double> initial{1.0};

    qmomkit::MomentSourceUpdateOptions options;
    options.ssp_rk3_options.accept = [](const std::vector<double>& state) { return state[0] >= 0.0; };
    options.ssp_rk3_options.max_retries = 10;
    options.ssp_rk3_options.min_time_step = 0.2;

    const qmomkit::MomentSourceUpdateResult result = qmomkit::advance_moments_with_source(
        initial, 0.0, 1.0, [](double, const std::vector<double>&) { return std::vector<double>{-8.0}; }, options);

    assert(!result.ok());
    assert(result.status == qmomkit::MomentSourceUpdateStatus::minimum_timestep_reached);
    assert(result.diagnostics.ssp_rk3_status == qmomkit::SspRk3Status::minimum_timestep_reached);
    assert(!result.diagnostics.accepted);
    assert(result.diagnostics.rejected);
    assert(result.diagnostics.retries == 2);
    assert_near(result.diagnostics.time_step_used, 0.25, 1.0e-15);
    assert(result.moments == initial);
}

void test_cqmom_moment_data_zero_source_preserves_state()
{
    const qmomkit::CqmomMomentData initial = create_moment_data();

    const qmomkit::CqmomMomentSourceUpdateResult result =
        qmomkit::advance_cqmom_moments_with_source(initial, 0.1, 0.2, [](double, const qmomkit::CqmomSourceContext& moments)
                                                   { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; });

    assert(result.ok());
    assert(same_layout(result.moments.layout, initial.layout));
    assert_vector_near(result.moments.values, initial.values, 1.0e-15);
}

void test_cqmom_source_context_exposes_moment_metadata()
{
    const qmomkit::CqmomMomentData initial = create_moment_data();
    bool inspected_context = false;

    const qmomkit::CqmomMomentSourceUpdateResult result = qmomkit::advance_cqmom_moments_with_source(
        initial, 0.1, 0.2,
        [&initial, &inspected_context](const double time, const qmomkit::CqmomSourceContext& moments)
        {
            assert_near(time, moments.time(), 1.0e-15);
            if (!inspected_context)
            {
                assert_near(moments.time(), 0.1, 1.0e-15);
                assert(moments.size() == initial.values.size());
                assert(moments.moment_count() == initial.values.size());
                assert(moments.dimension_count() == initial.layout.hierarchy.number_of_dimensions());
                assert(moments.hierarchy().node_counts() == initial.layout.hierarchy.node_counts());
                assert(moments.hierarchy().conditioning_order() == initial.layout.hierarchy.conditioning_order());
                assert(moments.indices() == initial.layout.required_indices);
                assert(moments.values() == initial.values);
                for (std::size_t i = 0; i < moments.size(); ++i)
                {
                    assert(moments.index(i) == initial.layout.required_indices[i]);
                    assert_near(moments.value(i), initial.values[i], 1.0e-15);
                    const std::optional<std::size_t> position = moments.position(initial.layout.required_indices[i]);
                    assert(position.has_value());
                    assert(*position == i);
                    const std::optional<double> value = moments.value(initial.layout.required_indices[i]);
                    assert(value.has_value());
                    assert_near(*value, initial.values[i], 1.0e-15);
                }
                assert(!moments.position(qmomkit::MultiIndex{99, 99}).has_value());
                assert(!moments.value(qmomkit::MultiIndex{99, 99}).has_value());
                inspected_context = true;
            }
            return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)};
        });

    assert(result.ok());
    assert(inspected_context);
}

void test_cqmom_source_context_validates_indices_and_values()
{
    const qmomkit::CqmomMomentData initial = create_moment_data();

    bool rejected_mismatched_values = false;
    try
    {
        static_cast<void>(qmomkit::CqmomSourceContext{0.0, initial.layout.hierarchy, initial.layout.required_indices, {1.0}});
    }
    catch (const std::invalid_argument& error)
    {
        rejected_mismatched_values = true;
        assert(std::string(error.what()).find("indices and values") != std::string::npos);
    }
    assert(rejected_mismatched_values);

    std::vector<qmomkit::MultiIndex> bad_indices = initial.layout.required_indices;
    bad_indices.front().push_back(0);
    bool rejected_bad_rank = false;
    try
    {
        static_cast<void>(qmomkit::CqmomSourceContext{0.0, initial.layout.hierarchy, bad_indices, initial.values});
    }
    catch (const std::invalid_argument& error)
    {
        rejected_bad_rank = true;
        assert(std::string(error.what()).find("dimension count") != std::string::npos);
    }
    assert(rejected_bad_rank);
}

void test_cqmom_moment_data_linear_source_matches_ssp_rk3_polynomial()
{
    const qmomkit::CqmomMomentData initial = create_moment_data();
    constexpr double lambda = -0.7;
    constexpr double time_step = 0.2;

    const qmomkit::CqmomMomentSourceUpdateResult result =
        qmomkit::advance_cqmom_moments_with_source(initial, 0.0, time_step,
                                                   [lambda](double, const qmomkit::CqmomSourceContext& moments)
                                                   {
                                                       std::vector<double> rate(moments.size(), 0.0);
                                                       for (std::size_t i = 0; i < moments.size(); ++i)
                                                       {
                                                           rate[i] = lambda * moments.value(i);
                                                       }
                                                       return qmomkit::CqmomMomentRate{rate};
                                                   });

    std::vector<double> expected = initial.values;
    const double factor = rk3_linear_factor(lambda, time_step);
    for (double& value : expected)
    {
        value *= factor;
    }

    assert(result.ok());
    assert_vector_near(result.moments.values, expected, 1.0e-14);
}

void test_cqmom_moment_data_rate_size_mismatch_fails()
{
    const qmomkit::CqmomMomentData initial = create_moment_data();

    const qmomkit::CqmomMomentSourceUpdateResult result =
        qmomkit::advance_cqmom_moments_with_source(initial, 0.0, 0.1, [](double, const qmomkit::CqmomSourceContext& moments)
                                                   { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size() + 1, 0.0)}; });

    assert(!result.ok());
    assert(result.status == qmomkit::MomentSourceUpdateStatus::source_failed);
    assert_vector_near(result.moments.values, initial.values, 1.0e-15);
}

void test_cqmom_moment_data_non_finite_source_fails()
{
    const qmomkit::CqmomMomentData initial = create_moment_data();

    const qmomkit::CqmomMomentSourceUpdateResult result =
        qmomkit::advance_cqmom_moments_with_source(initial, 0.0, 0.1,
                                                   [](double, const qmomkit::CqmomSourceContext& moments)
                                                   {
                                                       std::vector<double> rate(moments.size(), 0.0);
                                                       rate.front() = std::numeric_limits<double>::infinity();
                                                       return qmomkit::CqmomMomentRate{rate};
                                                   });

    assert(!result.ok());
    assert(result.status == qmomkit::MomentSourceUpdateStatus::source_failed);
    assert_vector_near(result.moments.values, initial.values, 1.0e-15);
}

void test_cqmom_moment_data_rk_failure_propagates_and_preserves_state()
{
    const qmomkit::CqmomMomentData initial = create_moment_data();

    qmomkit::MomentSourceUpdateOptions options;
    options.ssp_rk3_options.accept = [](const std::vector<double>&) { return false; };
    options.ssp_rk3_options.max_retries = 1;

    const qmomkit::CqmomMomentSourceUpdateResult result = qmomkit::advance_cqmom_moments_with_source(
        initial, 0.0, 0.1, [](double, const qmomkit::CqmomSourceContext& moments)
        { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; }, options);

    assert(!result.ok());
    assert(result.status == qmomkit::MomentSourceUpdateStatus::rejected_by_predicate);
    assert(result.diagnostics.ssp_rk3_status == qmomkit::SspRk3Status::rejected_by_predicate);
    assert(result.diagnostics.rejected);
    assert_vector_near(result.moments.values, initial.values, 1.0e-15);
}

void test_boundary_contribution_adds_to_moment_source()
{
    const std::vector<double> initial{1.0, 2.0, -1.0};
    constexpr double time_step = 0.25;

    const qmomkit::MomentSource source = [](double, const std::vector<double>& moments) { return std::vector<double>{1.0, 0.0, -0.5 * moments[2]}; };
    const qmomkit::MomentBoundaryContribution boundary = [](double, const std::vector<double>& moments)
    { return qmomkit::MomentBoundaryRate{std::vector<double>(moments.size(), 0.5)}; };

    const qmomkit::MomentSourceUpdateResult result = qmomkit::advance_moments_with_source(initial, 0.0, time_step, source, boundary, {});

    assert(result.ok());
    assert_vector_near(result.moments, {1.375, 2.125, -0.76497395833333337}, 1.0e-15);
}

void test_boundary_contribution_wrong_size_fails_with_boundary_message()
{
    const std::vector<double> initial{1.0, 2.0};
    const qmomkit::MomentSource source = [](double, const std::vector<double>& moments) { return std::vector<double>(moments.size(), 0.0); };
    const qmomkit::MomentBoundaryContribution boundary = [](double, const std::vector<double>& moments)
    { return qmomkit::MomentBoundaryRate{std::vector<double>(moments.size() + 1, 0.0)}; };

    const qmomkit::MomentSourceUpdateResult result = qmomkit::advance_moments_with_source(initial, 0.0, 0.1, source, boundary, {});

    assert(!result.ok());
    assert(result.status == qmomkit::MomentSourceUpdateStatus::source_failed);
    assert(result.message.find("boundary contribution") != std::string::npos);
    assert_vector_near(result.moments, initial, 1.0e-15);
}

void test_boundary_contribution_non_finite_fails_with_boundary_message()
{
    const std::vector<double> initial{1.0, 2.0};
    const qmomkit::MomentSource source = [](double, const std::vector<double>& moments) { return std::vector<double>(moments.size(), 0.0); };
    const qmomkit::MomentBoundaryContribution boundary = [](double, const std::vector<double>& moments)
    {
        std::vector<double> rate(moments.size(), 0.0);
        rate.front() = std::numeric_limits<double>::quiet_NaN();
        return qmomkit::MomentBoundaryRate{rate};
    };

    const qmomkit::MomentSourceUpdateResult result = qmomkit::advance_moments_with_source(initial, 0.0, 0.1, source, boundary, {});

    assert(!result.ok());
    assert(result.status == qmomkit::MomentSourceUpdateStatus::source_failed);
    assert(result.message.find("boundary contribution") != std::string::npos);
    assert_vector_near(result.moments, initial, 1.0e-15);
}

void test_boundary_contribution_callback_exception_fails_with_boundary_message()
{
    const std::vector<double> initial{1.0};
    const qmomkit::MomentSource source = [](double, const std::vector<double>& moments) { return std::vector<double>(moments.size(), 0.0); };
    const qmomkit::MomentBoundaryContribution boundary = [](double, const std::vector<double>&) -> qmomkit::MomentBoundaryRate
    { throw std::runtime_error("user model rejected the stage"); };

    const qmomkit::MomentSourceUpdateResult result = qmomkit::advance_moments_with_source(initial, 0.0, 0.1, source, boundary, {});

    assert(!result.ok());
    assert(result.status == qmomkit::MomentSourceUpdateStatus::source_failed);
    assert(result.message.find("boundary contribution") != std::string::npos);
    assert(result.message.find("user model rejected the stage") != std::string::npos);
    assert(result.diagnostics.message == result.message);
    assert_vector_near(result.moments, initial, 1.0e-15);
}

void test_source_callback_exception_keeps_source_message()
{
    const std::vector<double> initial{1.0};
    const qmomkit::MomentSource source = [](double, const std::vector<double>&) -> std::vector<double>
    { throw std::runtime_error("source model rejected the stage"); };
    const qmomkit::MomentBoundaryContribution boundary = [](double, const std::vector<double>& moments)
    { return qmomkit::MomentBoundaryRate{std::vector<double>(moments.size(), 0.0)}; };

    const qmomkit::MomentSourceUpdateResult result = qmomkit::advance_moments_with_source(initial, 0.0, 0.1, source, boundary, {});

    assert(!result.ok());
    assert(result.status == qmomkit::MomentSourceUpdateStatus::source_failed);
    assert(result.message.find("source model rejected the stage") != std::string::npos);
    assert(result.message.find("boundary contribution") == std::string::npos);
    assert_vector_near(result.moments, initial, 1.0e-15);
}

void test_cqmom_boundary_contribution_uses_source_context_ordering()
{
    const qmomkit::CqmomMomentData initial = create_moment_data();
    constexpr double time_step = 0.2;
    bool inspected_context = false;

    const qmomkit::CqmomMomentSource source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };
    const qmomkit::CqmomMomentBoundaryContribution boundary = [&initial, &inspected_context](double, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> rate(context.size(), 0.0);
        for (std::size_t i = 0; i < context.size(); ++i)
        {
            assert(context.index(i) == initial.layout.required_indices[i]);
            rate[i] = 0.1 + static_cast<double>(context.index(i)[0]) + 2.0 * static_cast<double>(context.index(i)[1]);
        }
        inspected_context = true;
        return qmomkit::CqmomMomentRate{rate};
    };

    const qmomkit::CqmomMomentSourceUpdateResult result = qmomkit::advance_cqmom_moments_with_source(initial, 0.0, time_step, source, boundary, {});

    std::vector<double> expected = initial.values;
    for (std::size_t i = 0; i < expected.size(); ++i)
    {
        expected[i] += time_step * (0.1 + static_cast<double>(initial.layout.required_indices[i][0]) +
                                    2.0 * static_cast<double>(initial.layout.required_indices[i][1]));
    }
    assert(result.ok());
    assert(inspected_context);
    assert_vector_near(result.moments.values, expected, 1.0e-15);
}

} // namespace

int main()
{
    test_accept_all_predicate();
    test_reject_negative_moments_predicate_accepts_nonnegative_update();
    test_retry_behavior_first_step_rejected_then_halved_step_succeeds();
    test_failure_behavior_predicate_rejects_all_attempts_and_preserves_original_state();
    test_failure_behavior_minimum_timestep_preserves_original_state();
    test_cqmom_moment_data_zero_source_preserves_state();
    test_cqmom_source_context_exposes_moment_metadata();
    test_cqmom_source_context_validates_indices_and_values();
    test_cqmom_moment_data_linear_source_matches_ssp_rk3_polynomial();
    test_cqmom_moment_data_rate_size_mismatch_fails();
    test_cqmom_moment_data_non_finite_source_fails();
    test_cqmom_moment_data_rk_failure_propagates_and_preserves_state();
    test_boundary_contribution_adds_to_moment_source();
    test_boundary_contribution_wrong_size_fails_with_boundary_message();
    test_boundary_contribution_non_finite_fails_with_boundary_message();
    test_boundary_contribution_callback_exception_fails_with_boundary_message();
    test_source_callback_exception_keeps_source_message();
    test_cqmom_boundary_contribution_uses_source_context_ordering();
    return 0;
}
