#include "qmomkit/qmomkit.hpp"
#include "test_case_dispatch.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace
{

static_assert(std::is_constructible_v<qmomkit::CqmomSourceContext, double, const qmomkit::PhaseSpaceHierarchy&, const qmomkit::MultiIndexTable&,
                                      const std::vector<double>&>);
static_assert(!std::is_constructible_v<qmomkit::CqmomSourceContext, const qmomkit::PhaseSpaceHierarchy&, const std::vector<qmomkit::MultiIndex>&,
                                       const std::vector<double>&>);

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

void set_coordinate_marginal(qmomkit::CqmomMomentData& moments, const std::size_t coordinate, const std::vector<double>& marginal)
{
    for (std::size_t position = 0; position < moments.layout.required_indices.size(); ++position)
    {
        const qmomkit::MultiIndexView index = moments.layout.required_indices[position];
        bool pure_coordinate = true;
        for (std::size_t dimension = 0; dimension < index.size(); ++dimension)
        {
            if (dimension != coordinate && index[dimension] != 0)
            {
                pure_coordinate = false;
                break;
            }
        }
        if (pure_coordinate)
        {
            assert(index[coordinate] < marginal.size());
            moments.values[position] = marginal[index[coordinate]];
        }
    }
}

qmomkit::CqmomMomentData create_coordinate_policy_moment_data(const std::size_t dimensions)
{
    std::vector<std::size_t> node_counts(dimensions, 2);
    std::vector<std::size_t> conditioning_order(dimensions, 0);
    for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
    {
        conditioning_order[dimension] = dimension;
    }
    qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout({node_counts, conditioning_order});
    return qmomkit::CqmomMomentData{layout, std::vector<double>(layout.required_indices.size(), 1.0)};
}

qmomkit::CqmomMomentSource zero_cqmom_source()
{
    return [](double, const qmomkit::CqmomSourceContext& moments) { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };
}

void test_accept_all_predicate()
{
    qmomkit::MomentSourceUpdateOptions options;
    options.ssp_rk3_options.accept = [](const std::vector<double>&) { return true; };
    options.ssp_rk3_options.max_retries = 4;

    const qmomkit::QmomMomentSourceUpdateResult result = qmomkit::advance_qmom_moments_with_source(
        {1.0}, 0.0, 0.2, [](double, const qmomkit::QmomSourceContext&) { return qmomkit::QmomMomentRate{std::vector<double>{3.0}}; }, options);

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

    const qmomkit::QmomMomentSourceUpdateResult result = qmomkit::advance_qmom_moments_with_source(
        {0.3, 1.0}, 0.0, 0.2, [](double, const qmomkit::QmomSourceContext&) { return qmomkit::QmomMomentRate{std::vector<double>{-0.5, -1.0}}; },
        options);

    assert(result.ok());
    assert(result.diagnostics.accepted);
    assert(result.diagnostics.retries == 0);
    assert_vector_near(result.moments, {0.2, 0.8}, 1.0e-15);
}

void test_retry_behavior_completes_requested_interval_from_unchanged_state()
{
    qmomkit::MomentSourceUpdateOptions options;
    std::size_t acceptance_calls = 0;
    options.ssp_rk3_options.accept = [&acceptance_calls](const std::vector<double>&)
    {
        ++acceptance_calls;
        return acceptance_calls != 1;
    };
    options.ssp_rk3_options.max_retries = 2;
    options.ssp_rk3_options.min_time_step = 0.25;

    std::vector<double> initial_stage_values;
    const qmomkit::QmomMomentSourceUpdateResult result = qmomkit::advance_qmom_moments_with_source(
        {0.0}, 0.0, 1.0,
        [&initial_stage_values](const double stage_time, const qmomkit::QmomSourceContext& state)
        {
            if (stage_time == 0.0)
            {
                initial_stage_values.push_back(state.value(0));
            }
            return qmomkit::QmomMomentRate{std::vector<double>{1.0}};
        },
        options);

    assert(result.ok());
    assert(result.diagnostics.accepted);
    assert(result.diagnostics.rejected);
    assert(result.diagnostics.retries == 1);
    assert_near(result.diagnostics.time_step_used, 0.5, 1.0e-15);
    assert_vector_near(initial_stage_values, {0.0, 0.0}, 1.0e-15);
    assert_near(result.moments[0], 1.0, 1.0e-15);
}

void test_repeated_reduced_steps_complete_requested_interval()
{
    qmomkit::MomentSourceUpdateOptions options;
    options.ssp_rk3_options.min_time_step = 0.1;
    options.ssp_rk3_options.max_retries = 1;
    std::size_t acceptance_calls = 0;
    options.ssp_rk3_options.accept = [&acceptance_calls](const std::vector<double>&)
    {
        ++acceptance_calls;
        return acceptance_calls == 2 || acceptance_calls == 4 || acceptance_calls == 6 || acceptance_calls == 7;
    };

    const qmomkit::QmomMomentSourceUpdateResult result = qmomkit::advance_qmom_moments_with_source(
        {0.0}, 0.0, 1.0, [](double, const qmomkit::QmomSourceContext&) { return qmomkit::QmomMomentRate{std::vector<double>{1.0}}; }, options);

    assert(result.ok());
    assert(result.diagnostics.rejected);
    assert(result.diagnostics.retries == 3);
    assert_near(result.diagnostics.time_step_used, 0.125, 1.0e-15);
    assert_near(result.moments[0], 1.0, 1.0e-15);
}

void test_failure_after_accepted_retry_preserves_original_state()
{
    const std::vector<double> initial{1.0, -2.0};

    qmomkit::MomentSourceUpdateOptions options;
    std::size_t acceptance_calls = 0;
    options.ssp_rk3_options.accept = [&acceptance_calls](const std::vector<double>&)
    {
        ++acceptance_calls;
        return acceptance_calls == 2;
    };
    options.ssp_rk3_options.max_retries = 1;
    options.ssp_rk3_options.min_time_step = 0.3;

    const qmomkit::QmomMomentSourceUpdateResult result = qmomkit::advance_qmom_moments_with_source(
        initial, 0.0, 1.0,
        [](double, const qmomkit::QmomSourceContext& moments) { return qmomkit::QmomMomentRate{std::vector<double>(moments.size(), 1.0)}; }, options);

    assert(!result.ok());
    assert(result.status == qmomkit::MomentSourceUpdateStatus::minimum_timestep_reached);
    assert(result.diagnostics.ssp_rk3_status == qmomkit::SspRk3Status::minimum_timestep_reached);
    assert(!result.diagnostics.accepted);
    assert(result.diagnostics.rejected);
    assert(result.diagnostics.retries == 1);
    assert_near(result.diagnostics.time_step_used, 0.5, 1.0e-15);
    assert(result.moments == initial);
}

void test_failure_behavior_minimum_timestep_preserves_original_state()
{
    const std::vector<double> initial{1.0};

    qmomkit::MomentSourceUpdateOptions options;
    options.ssp_rk3_options.accept = [](const std::vector<double>& state) { return state[0] >= 0.0; };
    options.ssp_rk3_options.max_retries = 10;
    options.ssp_rk3_options.min_time_step = 0.2;

    const qmomkit::QmomMomentSourceUpdateResult result = qmomkit::advance_qmom_moments_with_source(
        initial, 0.0, 1.0, [](double, const qmomkit::QmomSourceContext&) { return qmomkit::QmomMomentRate{std::vector<double>{-8.0}}; }, options);

    assert(!result.ok());
    assert(result.status == qmomkit::MomentSourceUpdateStatus::minimum_timestep_reached);
    assert(result.diagnostics.ssp_rk3_status == qmomkit::SspRk3Status::minimum_timestep_reached);
    assert(!result.diagnostics.accepted);
    assert(result.diagnostics.rejected);
    assert(result.diagnostics.retries == 2);
    assert_near(result.diagnostics.time_step_used, 0.25, 1.0e-15);
    assert(result.moments == initial);
}

void test_reduced_accepted_step_must_reduce_remaining_interval()
{
    qmomkit::MomentSourceUpdateOptions options;
    std::size_t acceptance_calls = 0;
    std::size_t source_calls = 0;
    options.ssp_rk3_options.accept = [&acceptance_calls](const std::vector<double>&)
    {
        ++acceptance_calls;
        return acceptance_calls % 55 == 0;
    };
    options.ssp_rk3_options.max_retries = 54;

    const std::vector<double> initial{0.0};
    const qmomkit::QmomMomentSourceUpdateResult result = qmomkit::advance_qmom_moments_with_source(
        initial, 0.0, 1.0,
        [&source_calls](double, const qmomkit::QmomSourceContext& context)
        {
            ++source_calls;
            return qmomkit::QmomMomentRate{std::vector<double>(context.size(), 0.0)};
        },
        options);

    assert(!result.ok());
    assert(result.status == qmomkit::MomentSourceUpdateStatus::source_failed);
    assert(result.moments == initial);
    assert(result.diagnostics.ssp_rk3_status == qmomkit::SspRk3Status::success);
    assert(result.diagnostics.accepted);
    assert(result.diagnostics.rejected);
    assert(result.diagnostics.retries == 54);
    assert(result.diagnostics.time_step_used == std::ldexp(1.0, -54));
    assert(result.message.find("representable interval progress") != std::string::npos);
    assert(source_calls == 165);
}

void test_nan_minimum_timestep_is_invalid()
{
    qmomkit::MomentSourceUpdateOptions options;
    options.ssp_rk3_options.min_time_step = std::numeric_limits<double>::quiet_NaN();
    std::size_t source_calls = 0;

    const qmomkit::QmomMomentSourceUpdateResult result = qmomkit::advance_qmom_moments_with_source(
        {1.0}, 0.0, 0.1,
        [&source_calls](double, const qmomkit::QmomSourceContext& context)
        {
            ++source_calls;
            return qmomkit::QmomMomentRate{std::vector<double>(context.size(), 0.0)};
        },
        options);

    assert(!result.ok());
    assert(result.status == qmomkit::MomentSourceUpdateStatus::invalid_input);
    assert(result.moments == std::vector<double>{1.0});
    assert(source_calls == 0);
}

void test_qmom_rejects_coordinate_realizability_policies()
{
    // coordinate_realizability_policies is a CQMOM-only option (see
    // MomentSourceUpdateOptions); the QMOM plain-vector path must reject it
    // rather than silently ignore it.
    qmomkit::RealizabilityPolicy policy;
    policy.check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;
    qmomkit::MomentSourceUpdateOptions options;
    options.coordinate_realizability_policies = {policy};
    std::size_t source_calls = 0;

    const qmomkit::QmomMomentSourceUpdateResult result = qmomkit::advance_qmom_moments_with_source(
        {1.0}, 0.0, 0.1,
        [&source_calls](double, const qmomkit::QmomSourceContext& context)
        {
            ++source_calls;
            return qmomkit::QmomMomentRate{std::vector<double>(context.size(), 0.0)};
        },
        options);

    assert(!result.ok());
    assert(result.status == qmomkit::MomentSourceUpdateStatus::invalid_input);
    assert(result.message.find("coordinate_realizability_policies") != std::string::npos);
    assert(result.moments == std::vector<double>{1.0});
    assert(source_calls == 0);

    // A vector that is merely sized but holds no configured policy (every entry
    // std::nullopt) is inert, matching the CQMOM path's own treatment of
    // unconfigured entries, and must not be rejected.
    qmomkit::MomentSourceUpdateOptions inert_options;
    inert_options.coordinate_realizability_policies.resize(3);
    const qmomkit::QmomMomentSourceUpdateResult inert_result = qmomkit::advance_qmom_moments_with_source(
        {1.0}, 0.0, 0.1, [](double, const qmomkit::QmomSourceContext& context)
        { return qmomkit::QmomMomentRate{std::vector<double>(context.size(), 0.0)}; }, inert_options);
    assert(inert_result.ok());
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
    constexpr double initial_time = 1.2345;
    bool inspected_context = false;

    const qmomkit::CqmomMomentSourceUpdateResult result = qmomkit::advance_cqmom_moments_with_source(
        initial, initial_time, 0.2,
        [&initial, &inspected_context](const double time, const qmomkit::CqmomSourceContext& moments)
        {
            assert_near(time, moments.time(), 1.0e-15);
            if (!inspected_context)
            {
                assert_near(time, initial_time, 1.0e-15);
                assert_near(moments.time(), initial_time, 1.0e-15);
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

    // A ragged index set can no longer be built at all: the flat table enforces
    // one rank for every row, so the mismatch is caught before any context sees it.
    std::vector<qmomkit::MultiIndex> ragged_indices = initial.layout.required_indices.to_multi_indices();
    ragged_indices.front().push_back(0);
    bool rejected_ragged = false;
    try
    {
        static_cast<void>(qmomkit::MultiIndexTable(ragged_indices));
    }
    catch (const std::invalid_argument& error)
    {
        rejected_ragged = true;
        assert(std::string(error.what()).find("same rank") != std::string::npos);
    }
    assert(rejected_ragged);

    // A uniformly wrong rank is well-formed as a table and must be rejected by
    // the context against the hierarchy dimension count.
    std::vector<qmomkit::MultiIndex> bad_rank_indices = initial.layout.required_indices.to_multi_indices();
    for (qmomkit::MultiIndex& index : bad_rank_indices)
    {
        index.push_back(0);
    }
    const qmomkit::MultiIndexTable bad_rank_table(bad_rank_indices);
    bool rejected_bad_rank = false;
    try
    {
        static_cast<void>(qmomkit::CqmomSourceContext{0.0, initial.layout.hierarchy, bad_rank_table, initial.values});
    }
    catch (const std::invalid_argument& error)
    {
        rejected_bad_rank = true;
        assert(std::string(error.what()).find("dimension count") != std::string::npos);
    }
    assert(rejected_bad_rank);
}

void test_cqmom_source_context_terminal_nodes_require_reconstruction_metadata()
{
    const qmomkit::CqmomMomentData initial = create_moment_data();
    const qmomkit::CqmomSourceContext context{0.0, initial.layout.hierarchy, initial.layout.required_indices, initial.values};
    bool rejected = false;
    try
    {
        context.for_each_terminal_node([](const double, const std::vector<double>&) {});
    }
    catch (const std::out_of_range& error)
    {
        rejected = true;
        assert(std::string(error.what()).find("source-stage quadrature") != std::string::npos);
    }
    assert(rejected);
}

void test_cqmom_source_context_generalized_moments_require_reconstruction_metadata()
{
    const qmomkit::CqmomMomentData initial = create_moment_data();
    const qmomkit::CqmomSourceContext context{0.0, initial.layout.hierarchy, initial.layout.required_indices, initial.values};
    bool rejected = false;
    try
    {
        static_cast<void>(context.generalized_moment({1.0, 0.0}));
    }
    catch (const std::out_of_range& error)
    {
        rejected = true;
        assert(std::string(error.what()).find("source-stage quadrature") != std::string::npos);
    }
    assert(rejected);
}

void test_cqmom_moment_data_linear_source_matches_ssp_rk3_polynomial()
{
    const qmomkit::CqmomMomentData initial = create_moment_data();
    constexpr double lambda = -0.7;
    constexpr double time_step = 0.2;

    const qmomkit::CqmomMomentSourceUpdateResult result =
        qmomkit::advance_cqmom_moments_with_source(initial, 0.0, time_step,
                                                   [](double, const qmomkit::CqmomSourceContext& moments)
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

void test_cqmom_coordinate_final_only_policy_rejects_final_candidate()
{
    qmomkit::CqmomMomentData initial = create_coordinate_policy_moment_data(2);
    set_coordinate_marginal(initial, 0, {1.0, 0.0, 1.0, 0.0});

    qmomkit::RealizabilityPolicy policy;
    policy.check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;
    policy.support = qmomkit::MomentSupport::real_line;

    qmomkit::MomentSourceUpdateOptions options;
    options.coordinate_realizability_policies = {policy};
    std::size_t source_calls = 0;
    std::size_t caller_accept_calls = 0;
    options.ssp_rk3_options.accept = [&caller_accept_calls](const std::vector<double>&)
    {
        ++caller_accept_calls;
        return true;
    };

    const qmomkit::CqmomMomentSourceUpdateResult result = qmomkit::advance_cqmom_moments_with_source(
        initial, 0.0, 0.2,
        [&source_calls](double, const qmomkit::CqmomSourceContext& moments)
        {
            ++source_calls;
            std::vector<double> rate(moments.size(), 0.0);
            const std::optional<std::size_t> position = moments.position({2, 0});
            assert(position.has_value());
            rate[*position] = -10.0;
            return qmomkit::CqmomMomentRate{rate};
        },
        options);

    assert(!result.ok());
    assert(result.status == qmomkit::MomentSourceUpdateStatus::rejected_by_predicate);
    assert(source_calls == 3);
    assert(caller_accept_calls == 1);
    assert_vector_near(result.moments.values, initial.values, 1.0e-15);
}

void test_cqmom_coordinate_every_stage_policy_rejects_stage_candidate()
{
    qmomkit::CqmomMomentData initial = create_coordinate_policy_moment_data(2);
    set_coordinate_marginal(initial, 0, {1.0, 0.0, 1.0, 0.0});

    qmomkit::RealizabilityPolicy policy;
    policy.check_mode = qmomkit::RealizabilityCheckMode::EveryStage;
    policy.support = qmomkit::MomentSupport::real_line;

    qmomkit::MomentSourceUpdateOptions options;
    options.coordinate_realizability_policies = {policy};
    std::size_t source_calls = 0;
    const qmomkit::CqmomMomentSourceUpdateResult result = qmomkit::advance_cqmom_moments_with_source(
        initial, 0.0, 0.2,
        [&source_calls](double, const qmomkit::CqmomSourceContext& moments)
        {
            ++source_calls;
            std::vector<double> rate(moments.size(), 0.0);
            const std::optional<std::size_t> position = moments.position({2, 0});
            assert(position.has_value());
            rate[*position] = -10.0;
            return qmomkit::CqmomMomentRate{rate};
        },
        options);

    assert(!result.ok());
    assert(result.status == qmomkit::MomentSourceUpdateStatus::rejected_by_predicate);
    assert(source_calls == 1);
    assert_vector_near(result.moments.values, initial.values, 1.0e-15);
}

void test_cqmom_coordinate_policies_check_configured_marginals_only()
{
    qmomkit::CqmomMomentData initial = create_coordinate_policy_moment_data(4);
    set_coordinate_marginal(initial, 0, {1.0, 0.0, 1.0, 0.0});
    set_coordinate_marginal(initial, 1, {1.0, 1.0, 1.0, 1.0});
    set_coordinate_marginal(initial, 2, {1.0, 0.5, 0.25, 0.125});
    set_coordinate_marginal(initial, 3, {1.0, 0.0, -1.0, 0.0});

    qmomkit::RealizabilityPolicy real_line;
    real_line.check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;
    real_line.support = qmomkit::MomentSupport::real_line;
    qmomkit::RealizabilityPolicy nonnegative = real_line;
    nonnegative.support = qmomkit::MomentSupport::nonnegative;
    qmomkit::RealizabilityPolicy bounded = real_line;
    bounded.support = qmomkit::MomentSupport::bounded_interval;
    bounded.interval = qmomkit::BoundedInterval{0.0, 1.0};

    qmomkit::MomentSourceUpdateOptions options;
    options.coordinate_realizability_policies = {real_line, nonnegative, bounded};
    const qmomkit::CqmomMomentSourceUpdateResult result = qmomkit::advance_cqmom_moments_with_source(initial, 0.0, 0.1, zero_cqmom_source(), options);

    assert(result.ok());
    assert_vector_near(result.moments.values, initial.values, 1.0e-15);
}

void test_cqmom_coordinate_policies_preserve_caller_predicates()
{
    qmomkit::CqmomMomentData initial = create_coordinate_policy_moment_data(2);
    set_coordinate_marginal(initial, 0, {1.0, 0.0, 1.0, 0.0});
    set_coordinate_marginal(initial, 1, {1.0, 1.0, 1.0, 1.0});

    qmomkit::RealizabilityPolicy final_only;
    final_only.check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;
    final_only.support = qmomkit::MomentSupport::real_line;
    qmomkit::MomentSourceUpdateOptions final_options;
    final_options.coordinate_realizability_policies = {final_only};
    std::size_t accept_calls = 0;
    final_options.ssp_rk3_options.accept = [&accept_calls](const std::vector<double>&)
    {
        ++accept_calls;
        return false;
    };
    const qmomkit::CqmomMomentSourceUpdateResult final_result =
        qmomkit::advance_cqmom_moments_with_source(initial, 0.0, 0.1, zero_cqmom_source(), final_options);
    assert(!final_result.ok());
    assert(final_result.status == qmomkit::MomentSourceUpdateStatus::rejected_by_predicate);
    assert(accept_calls == 1);

    qmomkit::RealizabilityPolicy every_stage;
    every_stage.check_mode = qmomkit::RealizabilityCheckMode::EveryStage;
    every_stage.support = qmomkit::MomentSupport::nonnegative;
    qmomkit::MomentSourceUpdateOptions stage_options;
    stage_options.coordinate_realizability_policies = {std::nullopt, every_stage};
    std::size_t stage_accept_calls = 0;
    stage_options.ssp_rk3_options.stage_accept = [&stage_accept_calls](qmomkit::SspRk3Stage, const std::vector<double>&)
    {
        ++stage_accept_calls;
        return false;
    };
    const qmomkit::CqmomMomentSourceUpdateResult stage_result =
        qmomkit::advance_cqmom_moments_with_source(initial, 0.0, 0.1, zero_cqmom_source(), stage_options);
    assert(!stage_result.ok());
    assert(stage_result.status == qmomkit::MomentSourceUpdateStatus::rejected_by_predicate);
    assert(stage_accept_calls == 1);
}

void test_cqmom_coordinate_policy_rejects_out_of_range_coordinate()
{
    const qmomkit::CqmomMomentData initial = create_coordinate_policy_moment_data(2);
    qmomkit::RealizabilityPolicy policy;
    policy.check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;

    qmomkit::MomentSourceUpdateOptions options;
    options.coordinate_realizability_policies.resize(3);
    options.coordinate_realizability_policies[2] = policy;
    const qmomkit::CqmomMomentSourceUpdateResult result = qmomkit::advance_cqmom_moments_with_source(initial, 0.0, 0.1, zero_cqmom_source(), options);

    assert(!result.ok());
    assert(result.status == qmomkit::MomentSourceUpdateStatus::invalid_input);
    assert(result.message.find("coordinate realizability policy index") != std::string::npos);
}

void test_cqmom_coordinate_policy_retry_completes_source_interval()
{
    qmomkit::CqmomMomentData initial = create_coordinate_policy_moment_data(2);
    set_coordinate_marginal(initial, 0, {1.0, 0.0, 1.0, 0.0});
    set_coordinate_marginal(initial, 1, {1.0, 0.0, -1.0, 0.0});

    qmomkit::RealizabilityPolicy policy;
    policy.check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;
    policy.support = qmomkit::MomentSupport::real_line;

    qmomkit::MomentSourceUpdateOptions options;
    options.coordinate_realizability_policies = {policy};
    options.ssp_rk3_options.max_retries = 1;
    options.ssp_rk3_options.min_time_step = 0.25;

    const qmomkit::CqmomMomentSourceUpdateResult result = qmomkit::advance_cqmom_moments_with_source(
        initial, 0.0, 1.0,
        [](double, const qmomkit::CqmomSourceContext& moments)
        {
            std::vector<double> rate(moments.size(), 0.0);
            const std::optional<std::size_t> position = moments.position({2, 0});
            assert(position.has_value());
            rate[*position] = -3.0 * moments.value(*position);
            return qmomkit::CqmomMomentRate{rate};
        },
        options);

    assert(result.ok());
    assert(result.diagnostics.retries == 1);
    assert_near(result.diagnostics.time_step_used, 0.5, 1.0e-15);
    assert_near(result.moments.moment({2, 0}), 1.0 / 256.0, 1.0e-15);
}

void test_boundary_contribution_adds_to_moment_source()
{
    const std::vector<double> initial{1.0, 2.0, -1.0};
    constexpr double time_step = 0.25;

    const qmomkit::QmomMomentSource source = [](double, const qmomkit::QmomSourceContext& moments)
    { return qmomkit::QmomMomentRate{std::vector<double>{1.0, 0.0, -0.5 * moments.value(2)}}; };
    const qmomkit::QmomMomentBoundaryContribution boundary = [](double, const qmomkit::QmomSourceContext& moments)
    { return qmomkit::QmomMomentRate{std::vector<double>(moments.size(), 0.5)}; };

    const qmomkit::QmomMomentSourceUpdateResult result = qmomkit::advance_qmom_moments_with_source(initial, 0.0, time_step, source, boundary, {});

    assert(result.ok());
    assert_vector_near(result.moments, {1.375, 2.125, -0.76497395833333337}, 1.0e-15);
}

void test_boundary_contribution_wrong_size_fails_with_boundary_message()
{
    const std::vector<double> initial{1.0, 2.0};
    const qmomkit::QmomMomentSource source = [](double, const qmomkit::QmomSourceContext& moments)
    { return qmomkit::QmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };
    const qmomkit::QmomMomentBoundaryContribution boundary = [](double, const qmomkit::QmomSourceContext& moments)
    { return qmomkit::QmomMomentRate{std::vector<double>(moments.size() + 1, 0.0)}; };

    const qmomkit::QmomMomentSourceUpdateResult result = qmomkit::advance_qmom_moments_with_source(initial, 0.0, 0.1, source, boundary, {});

    assert(!result.ok());
    assert(result.status == qmomkit::MomentSourceUpdateStatus::source_failed);
    assert(result.message.find("boundary contribution") != std::string::npos);
    assert_vector_near(result.moments, initial, 1.0e-15);
}

void test_boundary_contribution_non_finite_fails_with_boundary_message()
{
    const std::vector<double> initial{1.0, 2.0};
    const qmomkit::QmomMomentSource source = [](double, const qmomkit::QmomSourceContext& moments)
    { return qmomkit::QmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };
    const qmomkit::QmomMomentBoundaryContribution boundary = [](double, const qmomkit::QmomSourceContext& moments)
    {
        std::vector<double> rate(moments.size(), 0.0);
        rate.front() = std::numeric_limits<double>::quiet_NaN();
        return qmomkit::QmomMomentRate{rate};
    };

    const qmomkit::QmomMomentSourceUpdateResult result = qmomkit::advance_qmom_moments_with_source(initial, 0.0, 0.1, source, boundary, {});

    assert(!result.ok());
    assert(result.status == qmomkit::MomentSourceUpdateStatus::source_failed);
    assert(result.message.find("boundary contribution") != std::string::npos);
    assert_vector_near(result.moments, initial, 1.0e-15);
}

void test_boundary_contribution_callback_exception_fails_with_boundary_message()
{
    const std::vector<double> initial{1.0};
    const qmomkit::QmomMomentSource source = [](double, const qmomkit::QmomSourceContext& moments)
    { return qmomkit::QmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };
    const qmomkit::QmomMomentBoundaryContribution boundary = [](double, const qmomkit::QmomSourceContext&) -> qmomkit::QmomMomentRate
    { throw std::runtime_error("user model rejected the stage"); };

    const qmomkit::QmomMomentSourceUpdateResult result = qmomkit::advance_qmom_moments_with_source(initial, 0.0, 0.1, source, boundary, {});

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
    const qmomkit::QmomMomentSource source = [](double, const qmomkit::QmomSourceContext&) -> qmomkit::QmomMomentRate
    { throw std::runtime_error("source model rejected the stage"); };
    const qmomkit::QmomMomentBoundaryContribution boundary = [](double, const qmomkit::QmomSourceContext& moments)
    { return qmomkit::QmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };

    const qmomkit::QmomMomentSourceUpdateResult result = qmomkit::advance_qmom_moments_with_source(initial, 0.0, 0.1, source, boundary, {});

    assert(!result.ok());
    assert(result.status == qmomkit::MomentSourceUpdateStatus::source_failed);
    assert(result.message.find("source model rejected the stage") != std::string::npos);
    assert(result.message.find("boundary contribution") == std::string::npos);
    assert_vector_near(result.moments, initial, 1.0e-15);
}

void test_cqmom_boundary_contribution_uses_source_context_ordering()
{
    const qmomkit::CqmomMomentData initial = create_moment_data();
    constexpr double initial_time = 1.2345;
    constexpr double time_step = 0.2;
    bool inspected_context = false;
    std::vector<double> source_times;
    std::vector<double> boundary_times;

    const qmomkit::CqmomMomentSource source = [&source_times](const double time, const qmomkit::CqmomSourceContext& moments)
    {
        assert_near(time, moments.time(), 1.0e-15);
        source_times.push_back(time);
        return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)};
    };
    const qmomkit::CqmomMomentBoundaryContribution boundary =
        [&initial, &inspected_context, &boundary_times](const double time, const qmomkit::CqmomSourceContext& context)
    {
        assert_near(time, context.time(), 1.0e-15);
        boundary_times.push_back(time);
        std::vector<double> rate(context.size(), 0.0);
        for (std::size_t i = 0; i < context.size(); ++i)
        {
            assert(context.index(i) == initial.layout.required_indices[i]);
            rate[i] = 0.1 + static_cast<double>(context.index(i)[0]) + 2.0 * static_cast<double>(context.index(i)[1]);
        }
        inspected_context = true;
        return qmomkit::CqmomMomentRate{rate};
    };

    const qmomkit::CqmomMomentSourceUpdateResult result =
        qmomkit::advance_cqmom_moments_with_source(initial, initial_time, time_step, source, boundary, {});

    std::vector<double> expected = initial.values;
    for (std::size_t i = 0; i < expected.size(); ++i)
    {
        expected[i] += time_step * (0.1 + static_cast<double>(initial.layout.required_indices[i][0]) +
                                    2.0 * static_cast<double>(initial.layout.required_indices[i][1]));
    }
    assert(result.ok());
    assert(inspected_context);
    assert(!source_times.empty());
    assert(source_times == boundary_times);
    assert_near(source_times.front(), initial_time, 1.0e-15);
    assert_vector_near(result.moments.values, expected, 1.0e-15);
}

void test_cqmom_source_context_index_lookup_matches_a_scan()
{
    // position() is backed by an ordered view of the entries rather than a
    // scan. It must still agree with a scan entry for entry, report a missing
    // index as absent, and -- because a caller-supplied table may repeat an
    // index -- report the first occurrence of a repeated one.
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    qmomkit::MultiIndexTable indices;
    for (const qmomkit::MultiIndex& index : std::vector<qmomkit::MultiIndex>{{0, 0}, {1, 0}, {0, 1}, {2, 1}, {1, 0}, {3, 2}, {0, 1}})
    {
        indices.push_back(index);
    }
    const std::vector<double> values{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
    const qmomkit::CqmomSourceContext context(0.0, hierarchy, indices, values);

    for (std::size_t entry = 0; entry < indices.size(); ++entry)
    {
        const qmomkit::MultiIndexView index = indices[entry];
        const auto scanned = std::find(indices.begin(), indices.end(), index);
        const std::size_t expected = static_cast<std::size_t>(std::distance(indices.begin(), scanned));

        const std::optional<std::size_t> position = context.position(index);
        assert(position.has_value());
        assert(*position == expected);

        const std::optional<double> value = context.value(index);
        assert(value.has_value());
        assert_near(*value, values[expected], 1.0e-15);
    }

    // The repeated entries resolve to their first occurrence, not their last.
    assert(*context.position(qmomkit::MultiIndex{1, 0}) == 1);
    assert(*context.position(qmomkit::MultiIndex{0, 1}) == 2);

    // The table above is too short to discriminate: a sort of a few entries is
    // stable in practice, so first-occurrence order survives by luck rather than
    // by the ordering being built to preserve it. Repeat the check on a table
    // long enough to be sorted by an unstable algorithm, with every index
    // duplicated and the copies far apart.
    qmomkit::MultiIndexTable repeated;
    std::vector<qmomkit::MultiIndex> distinct;
    for (std::size_t first = 0; first < 10; ++first)
    {
        for (std::size_t second = 0; second < 10; ++second)
        {
            distinct.push_back(qmomkit::MultiIndex{first, second});
        }
    }
    for (auto entry = distinct.rbegin(); entry != distinct.rend(); ++entry)
    {
        repeated.push_back(*entry);
    }
    for (const qmomkit::MultiIndex& entry : distinct)
    {
        repeated.push_back(entry);
    }
    const std::vector<double> repeated_values(repeated.size(), 1.0);
    const qmomkit::CqmomSourceContext repeated_context(0.0, hierarchy, repeated, repeated_values);

    for (const qmomkit::MultiIndex& entry : distinct)
    {
        const auto scanned = std::find(repeated.begin(), repeated.end(), qmomkit::MultiIndexView(entry));
        const std::size_t expected = static_cast<std::size_t>(std::distance(repeated.begin(), scanned));
        const std::optional<std::size_t> position = repeated_context.position(entry);
        assert(position.has_value());
        assert(*position == expected);
    }

    // A query whose rank does not match the table is absent, not a match.
    assert(!context.position(qmomkit::MultiIndex{0, 0, 0}).has_value());

    // Absent indices, including ones that sort before and after every entry.
    assert(!context.position(qmomkit::MultiIndex{9, 9}).has_value());
    assert(!context.position(qmomkit::MultiIndex{2, 0}).has_value());
    assert(!context.value(qmomkit::MultiIndex{9, 9}).has_value());
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"cqmom_source_context_index_lookup_matches_a_scan", test_cqmom_source_context_index_lookup_matches_a_scan},
            {"accept_all_predicate", test_accept_all_predicate},
            {"reject_negative_moments_predicate_accepts_nonnegative_update", test_reject_negative_moments_predicate_accepts_nonnegative_update},
            {"retry_behavior_completes_requested_interval_from_unchanged_state",
             test_retry_behavior_completes_requested_interval_from_unchanged_state},
            {"repeated_reduced_steps_complete_requested_interval", test_repeated_reduced_steps_complete_requested_interval},
            {"failure_after_accepted_retry_preserves_original_state", test_failure_after_accepted_retry_preserves_original_state},
            {"failure_behavior_minimum_timestep_preserves_original_state", test_failure_behavior_minimum_timestep_preserves_original_state},
            {"reduced_accepted_step_must_reduce_remaining_interval", test_reduced_accepted_step_must_reduce_remaining_interval},
            {"nan_minimum_timestep_is_invalid", test_nan_minimum_timestep_is_invalid},
            {"qmom_rejects_coordinate_realizability_policies", test_qmom_rejects_coordinate_realizability_policies},
            {"cqmom_moment_data_zero_source_preserves_state", test_cqmom_moment_data_zero_source_preserves_state},
            {"cqmom_source_context_exposes_moment_metadata", test_cqmom_source_context_exposes_moment_metadata},
            {"cqmom_source_context_validates_indices_and_values", test_cqmom_source_context_validates_indices_and_values},
            {"cqmom_source_context_terminal_nodes_require_reconstruction_metadata",
             test_cqmom_source_context_terminal_nodes_require_reconstruction_metadata},
            {"cqmom_source_context_generalized_moments_require_reconstruction_metadata",
             test_cqmom_source_context_generalized_moments_require_reconstruction_metadata},
            {"cqmom_moment_data_linear_source_matches_ssp_rk3_polynomial", test_cqmom_moment_data_linear_source_matches_ssp_rk3_polynomial},
            {"cqmom_moment_data_rate_size_mismatch_fails", test_cqmom_moment_data_rate_size_mismatch_fails},
            {"cqmom_moment_data_non_finite_source_fails", test_cqmom_moment_data_non_finite_source_fails},
            {"cqmom_moment_data_rk_failure_propagates_and_preserves_state", test_cqmom_moment_data_rk_failure_propagates_and_preserves_state},
            {"cqmom_coordinate_final_only_policy_rejects_final_candidate", test_cqmom_coordinate_final_only_policy_rejects_final_candidate},
            {"cqmom_coordinate_every_stage_policy_rejects_stage_candidate", test_cqmom_coordinate_every_stage_policy_rejects_stage_candidate},
            {"cqmom_coordinate_policies_check_configured_marginals_only", test_cqmom_coordinate_policies_check_configured_marginals_only},
            {"cqmom_coordinate_policies_preserve_caller_predicates", test_cqmom_coordinate_policies_preserve_caller_predicates},
            {"cqmom_coordinate_policy_rejects_out_of_range_coordinate", test_cqmom_coordinate_policy_rejects_out_of_range_coordinate},
            {"cqmom_coordinate_policy_retry_completes_source_interval", test_cqmom_coordinate_policy_retry_completes_source_interval},
            {"boundary_contribution_adds_to_moment_source", test_boundary_contribution_adds_to_moment_source},
            {"boundary_contribution_wrong_size_fails_with_boundary_message", test_boundary_contribution_wrong_size_fails_with_boundary_message},
            {"boundary_contribution_non_finite_fails_with_boundary_message", test_boundary_contribution_non_finite_fails_with_boundary_message},
            {"boundary_contribution_callback_exception_fails_with_boundary_message",
             test_boundary_contribution_callback_exception_fails_with_boundary_message},
            {"source_callback_exception_keeps_source_message", test_source_callback_exception_keeps_source_message},
            {"cqmom_boundary_contribution_uses_source_context_ordering", test_cqmom_boundary_contribution_uses_source_context_ordering},
        });
}
