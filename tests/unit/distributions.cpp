#include "qmomkit/distributions.hpp"

#include "qmomkit/cqmom.hpp"
#include "test_case_dispatch.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace
{

void assert_near(const double actual, const double expected, const double tolerance)
{
    assert(std::abs(actual - expected) <= tolerance);
}

void assert_vector_equal(const std::vector<std::size_t>& actual, const std::vector<std::size_t>& expected)
{
    assert(actual == expected);
}

void test_selected_raw_moment_result_preserves_request_order()
{
    const qmomkit::RawMomentRequest request{{{3}, {0}, {5}}};
    const qmomkit::SelectedRawMomentResult result = qmomkit::make_selected_raw_moment_result(request, {30.0, 0.0, 50.0});

    assert(result.ok());
    assert(result.indices == request.indices);
    assert(result.values.size() == request.indices.size());
    assert_vector_equal(result.indices[0], {3});
    assert_vector_equal(result.indices[1], {0});
    assert_vector_equal(result.indices[2], {5});
    assert_near(result.values[0], 30.0, 0.0);
    assert_near(result.values[1], 0.0, 0.0);
    assert_near(result.values[2], 50.0, 0.0);
}

void test_selected_raw_moment_result_accepts_noncontiguous_multidimensional_indices()
{
    const qmomkit::RawMomentRequest request{{{1, 0}, {0, 4}, {3, 2}}};
    const qmomkit::SelectedRawMomentResult result = qmomkit::make_selected_raw_moment_result(request, {10.0, 4.0, 32.0});

    assert(result.ok());
    assert(result.indices == request.indices);
    assert(result.values.size() == 3);
    assert_vector_equal(result.indices[0], {1, 0});
    assert_vector_equal(result.indices[1], {0, 4});
    assert_vector_equal(result.indices[2], {3, 2});
    assert_near(result.values[0], 10.0, 0.0);
    assert_near(result.values[1], 4.0, 0.0);
    assert_near(result.values[2], 32.0, 0.0);
}

void test_selected_raw_moment_result_preserves_duplicate_indices()
{
    const qmomkit::RawMomentRequest request{{{2}, {2}, {0}}};
    const qmomkit::SelectedRawMomentResult result = qmomkit::make_selected_raw_moment_result(request, {4.0, 4.0, 1.0});

    assert(result.ok());
    assert(result.indices == request.indices);
    assert(result.values.size() == 3);
    assert_vector_equal(result.indices[0], {2});
    assert_vector_equal(result.indices[1], {2});
    assert_near(result.values[0], result.values[1], 0.0);
}

void test_selected_raw_moment_empty_request_is_valid()
{
    const qmomkit::SelectedRawMomentResult result = qmomkit::make_selected_raw_moment_result(qmomkit::RawMomentRequest{}, {});

    assert(result.ok());
    assert(result.indices.empty());
    assert(result.values.empty());
}

void test_selected_raw_moment_invalid_requests()
{
    const qmomkit::SelectedRawMomentResult mismatched_count = qmomkit::make_selected_raw_moment_result({{{0}, {1}}}, {1.0});
    assert(!mismatched_count.ok());
    assert(mismatched_count.status == qmomkit::RawMomentStatus::invalid_input);

    const qmomkit::SelectedRawMomentResult zero_rank = qmomkit::make_selected_raw_moment_result({{{}}}, {1.0});
    assert(!zero_rank.ok());
    assert(zero_rank.status == qmomkit::RawMomentStatus::invalid_input);

    const qmomkit::SelectedRawMomentResult rank_mismatch = qmomkit::make_selected_raw_moment_result({{{0}, {1, 2}}}, {1.0, 2.0});
    assert(!rank_mismatch.ok());
    assert(rank_mismatch.status == qmomkit::RawMomentStatus::invalid_input);
    assert_vector_equal(rank_mismatch.failed_order, {1, 2});

    const qmomkit::SelectedRawMomentResult order_too_large =
        qmomkit::make_selected_raw_moment_result({{{1, std::numeric_limits<std::size_t>::max()}}}, {1.0});
    assert(!order_too_large.ok());
    assert(order_too_large.status == qmomkit::RawMomentStatus::invalid_input);
    assert(order_too_large.failed_dimension == 1);
    assert_vector_equal(order_too_large.failed_order, {1, std::numeric_limits<std::size_t>::max()});

    const qmomkit::SelectedRawMomentResult nonfinite_value =
        qmomkit::make_selected_raw_moment_result({{{3}}}, {std::numeric_limits<double>::infinity()});
    assert(!nonfinite_value.ok());
    assert(nonfinite_value.status == qmomkit::RawMomentStatus::overflow);
    assert_vector_equal(nonfinite_value.failed_order, {3});
}

void assert_selected_matches_all_orders(const qmomkit::DistributionSpecification& distribution, const qmomkit::RawMomentRequest& request,
                                        const double total_weight, const std::size_t maximum_order, const double tolerance)
{
    const qmomkit::RawMomentResult all_orders = qmomkit::compute_raw_moments(distribution, maximum_order, total_weight);
    const qmomkit::SelectedRawMomentResult selected = qmomkit::compute_raw_moments(distribution, request, total_weight);

    assert(all_orders.ok());
    assert(selected.ok());
    assert(selected.indices == request.indices);
    assert(selected.values.size() == request.indices.size());
    for (std::size_t index = 0; index < request.indices.size(); ++index)
    {
        assert(request.indices[index].size() == 1);
        // Relative, because high-order moments reach 1e13 and more, where an
        // absolute bound would demand bit equality between the two paths.
        const double expected = all_orders.moments[request.indices[index][0]];
        assert(std::abs(selected.values[index] - expected) <= tolerance * std::abs(expected));
    }
}

void test_selected_distribution_contiguous_orders_match_all_orders()
{
    const qmomkit::RawMomentRequest request{{{0}, {1}, {2}, {3}}};
    assert_selected_matches_all_orders(qmomkit::NormalDistribution{0.5, 1.25}, request, 2.0, 3, 1.0e-13);
    assert_selected_matches_all_orders(qmomkit::LogNormalDistribution{0.1, 0.2}, request, 2.0, 3, 1.0e-13);
    assert_selected_matches_all_orders(qmomkit::BetaDistribution{2.0, 5.0}, request, 2.0, 3, 1.0e-13);
    assert_selected_matches_all_orders(qmomkit::GammaDistribution{1.5, 0.75}, request, 2.0, 3, 1.0e-13);
    assert_selected_matches_all_orders(qmomkit::ChiDistribution{4.0}, request, 2.0, 3, 1.0e-13);
    assert_selected_matches_all_orders(qmomkit::ChiSquaredDistribution{4.0}, request, 2.0, 3, 1.0e-13);
    assert_selected_matches_all_orders(qmomkit::WeibullDistribution{1.5, 0.75}, request, 2.0, 3, 1.0e-13);
}

void test_selected_distribution_noncontiguous_orders_match_all_orders()
{
    const qmomkit::RawMomentRequest request{{{0}, {2}, {5}}};
    assert_selected_matches_all_orders(qmomkit::NormalDistribution{0.25, 0.75}, request, 1.5, 5, 1.0e-12);
    assert_selected_matches_all_orders(qmomkit::LogNormalDistribution{-0.1, 0.15}, request, 1.5, 5, 1.0e-12);
    assert_selected_matches_all_orders(qmomkit::BetaDistribution{3.0, 4.0}, request, 1.5, 5, 1.0e-12);
    assert_selected_matches_all_orders(qmomkit::GammaDistribution{2.0, 0.5}, request, 1.5, 5, 1.0e-12);
    assert_selected_matches_all_orders(qmomkit::ChiDistribution{5.0}, request, 1.5, 5, 1.0e-12);
    assert_selected_matches_all_orders(qmomkit::ChiSquaredDistribution{5.0}, request, 1.5, 5, 1.0e-12);
    assert_selected_matches_all_orders(qmomkit::WeibullDistribution{2.0, 0.5}, request, 1.5, 5, 1.0e-12);
}

void test_selected_distribution_reordered_and_duplicate_orders_match_all_orders()
{
    const qmomkit::RawMomentRequest request{{{5}, {0}, {2}, {2}, {1}}};
    assert_selected_matches_all_orders(qmomkit::NormalDistribution{0.25, 0.75}, request, 1.25, 5, 1.0e-12);
    assert_selected_matches_all_orders(qmomkit::LogNormalDistribution{-0.1, 0.15}, request, 1.25, 5, 1.0e-12);
    assert_selected_matches_all_orders(qmomkit::BetaDistribution{3.0, 4.0}, request, 1.25, 5, 1.0e-12);
    assert_selected_matches_all_orders(qmomkit::GammaDistribution{2.0, 0.5}, request, 1.25, 5, 1.0e-12);
    assert_selected_matches_all_orders(qmomkit::ChiDistribution{5.0}, request, 1.25, 5, 1.0e-12);
    assert_selected_matches_all_orders(qmomkit::ChiSquaredDistribution{5.0}, request, 1.25, 5, 1.0e-12);
    assert_selected_matches_all_orders(qmomkit::WeibullDistribution{2.0, 0.5}, request, 1.25, 5, 1.0e-12);
}

void test_selected_distribution_high_order_gaps_match_all_orders()
{
    const qmomkit::RawMomentRequest request{{{0}, {1}, {12}}};
    assert_selected_matches_all_orders(qmomkit::NormalDistribution{0.1, 0.5}, request, 1.0, 12, 1.0e-12);
    assert_selected_matches_all_orders(qmomkit::LogNormalDistribution{-0.05, 0.1}, request, 1.0, 12, 1.0e-12);
    assert_selected_matches_all_orders(qmomkit::BetaDistribution{2.0, 3.0}, request, 1.0, 12, 1.0e-12);
    assert_selected_matches_all_orders(qmomkit::GammaDistribution{1.5, 0.25}, request, 1.0, 12, 1.0e-12);
    assert_selected_matches_all_orders(qmomkit::ChiDistribution{3.0}, request, 1.0, 12, 1.0e-12);
    assert_selected_matches_all_orders(qmomkit::ChiSquaredDistribution{3.0}, request, 1.0, 12, 1.0e-12);
    assert_selected_matches_all_orders(qmomkit::WeibullDistribution{2.0, 0.25}, request, 1.0, 12, 1.0e-12);
}

void test_selected_distribution_zero_total_weight()
{
    const qmomkit::RawMomentRequest request{{{3}, {0}, {5}}};
    const qmomkit::SelectedRawMomentResult normal = qmomkit::compute_raw_moments(qmomkit::NormalDistribution{1.0e308, 1.0}, request, 0.0);
    assert(normal.ok());
    assert(normal.indices == request.indices);
    for (double value : normal.values)
    {
        assert_near(value, 0.0, 0.0);
    }

    const qmomkit::SelectedRawMomentResult gamma = qmomkit::compute_raw_moments(qmomkit::GammaDistribution{1.0, 0.0}, request, 0.0);
    assert(!gamma.ok());
    assert(gamma.status == qmomkit::RawMomentStatus::invalid_input);
}

void test_selected_distribution_invalid_request_rank()
{
    const qmomkit::SelectedRawMomentResult zero_rank =
        qmomkit::compute_raw_moments(qmomkit::NormalDistribution{0.0, 1.0}, qmomkit::RawMomentRequest{{{}}});
    assert(!zero_rank.ok());
    assert(zero_rank.status == qmomkit::RawMomentStatus::invalid_input);

    const qmomkit::SelectedRawMomentResult multi_rank =
        qmomkit::compute_raw_moments(qmomkit::NormalDistribution{0.0, 1.0}, qmomkit::RawMomentRequest{{{0, 1}}});
    assert(!multi_rank.ok());
    assert(multi_rank.status == qmomkit::RawMomentStatus::invalid_input);
    assert_vector_equal(multi_rank.failed_order, {0, 1});
}

void test_selected_distribution_overflow_and_underflow_failures()
{
    const qmomkit::SelectedRawMomentResult overflow =
        qmomkit::compute_raw_moments(qmomkit::LogNormalDistribution{100.0, 10.0}, qmomkit::RawMomentRequest{{{20}}});
    assert(!overflow.ok());
    assert(overflow.status == qmomkit::RawMomentStatus::overflow);
    assert(overflow.failed_dimension == 0);
    assert_vector_equal(overflow.failed_order, {20});

    const qmomkit::SelectedRawMomentResult underflow =
        qmomkit::compute_raw_moments(qmomkit::LogNormalDistribution{-1000.0, 0.1}, qmomkit::RawMomentRequest{{{1}}});
    assert(!underflow.ok());
    assert(underflow.status == qmomkit::RawMomentStatus::underflow);
    assert(underflow.failed_dimension == 0);
    assert_vector_equal(underflow.failed_order, {1});
}

void test_variant_dispatch_and_nonunit_weight()
{
    const qmomkit::DistributionSpecification distribution = qmomkit::NormalDistribution{2.0, 3.0};
    const qmomkit::RawMomentResult result = qmomkit::compute_raw_moments(distribution, 2, 4.0);

    assert(result.ok());
    assert(result.moments.size() == 3);
    assert_near(result.moments[0], 4.0, 0.0);
    assert_near(result.moments[1], 8.0, 1.0e-14);
    assert_near(result.moments[2], 52.0, 1.0e-14);
}

void test_maximum_order_zero()
{
    const qmomkit::RawMomentResult result = qmomkit::compute_raw_moments(qmomkit::GammaDistribution{2.0, 3.0}, 0, 7.5);

    assert(result.ok());
    assert(result.moments.size() == 1);
    assert_near(result.moments[0], 7.5, 1.0e-14);
}

void test_zero_total_weight_returns_zero_sequence()
{
    const qmomkit::RawMomentResult log_normal = qmomkit::compute_raw_moments(qmomkit::LogNormalDistribution{-1000.0, 0.5}, 4, 0.0);

    assert(log_normal.ok());
    assert(log_normal.moments.size() == 5);
    for (double moment : log_normal.moments.values())
    {
        assert_near(moment, 0.0, 0.0);
    }

    const qmomkit::RawMomentResult normal = qmomkit::compute_raw_moments(qmomkit::NormalDistribution{1.0e308, 1.0}, 2, 0.0);
    assert(normal.ok());
    assert(normal.moments.size() == 3);
    for (double moment : normal.moments.values())
    {
        assert_near(moment, 0.0, 0.0);
    }
}

void test_invalid_total_weight()
{
    const qmomkit::RawMomentResult negative = qmomkit::compute_raw_moments(qmomkit::BetaDistribution{1.0, 1.0}, 2, -1.0);
    assert(!negative.ok());
    assert(negative.status == qmomkit::RawMomentStatus::invalid_input);
    assert(!negative.message.empty());

    const qmomkit::RawMomentResult nonfinite =
        qmomkit::compute_raw_moments(qmomkit::BetaDistribution{1.0, 1.0}, 2, std::numeric_limits<double>::infinity());
    assert(!nonfinite.ok());
    assert(nonfinite.status == qmomkit::RawMomentStatus::invalid_input);

    const qmomkit::RawMomentResult not_a_number =
        qmomkit::compute_raw_moments(qmomkit::BetaDistribution{1.0, 1.0}, 2, std::numeric_limits<double>::quiet_NaN());
    assert(!not_a_number.ok());
    assert(not_a_number.status == qmomkit::RawMomentStatus::invalid_input);
}

void test_invalid_distribution_parameters()
{
    assert(qmomkit::compute_raw_moments(qmomkit::NormalDistribution{0.0, 0.0}, 2).status == qmomkit::RawMomentStatus::invalid_input);
    assert(qmomkit::compute_raw_moments(qmomkit::NormalDistribution{std::numeric_limits<double>::infinity(), 1.0}, 2).status ==
           qmomkit::RawMomentStatus::invalid_input);
    assert(qmomkit::compute_raw_moments(qmomkit::LogNormalDistribution{0.0, -1.0}, 2).status == qmomkit::RawMomentStatus::invalid_input);
    assert(qmomkit::compute_raw_moments(qmomkit::LogNormalDistribution{std::numeric_limits<double>::quiet_NaN(), 1.0}, 2).status ==
           qmomkit::RawMomentStatus::invalid_input);
    assert(qmomkit::compute_raw_moments(qmomkit::BetaDistribution{0.0, 1.0}, 2).status == qmomkit::RawMomentStatus::invalid_input);
    assert(qmomkit::compute_raw_moments(qmomkit::BetaDistribution{1.0, std::numeric_limits<double>::infinity()}, 2).status ==
           qmomkit::RawMomentStatus::invalid_input);
    assert(qmomkit::compute_raw_moments(qmomkit::GammaDistribution{1.0, 0.0}, 2).status == qmomkit::RawMomentStatus::invalid_input);
    assert(qmomkit::compute_raw_moments(qmomkit::GammaDistribution{std::numeric_limits<double>::quiet_NaN(), 1.0}, 2).status ==
           qmomkit::RawMomentStatus::invalid_input);
    assert(qmomkit::compute_raw_moments(qmomkit::ChiDistribution{0.0}, 2).status == qmomkit::RawMomentStatus::invalid_input);
    assert(qmomkit::compute_raw_moments(qmomkit::ChiDistribution{std::numeric_limits<double>::infinity()}, 2).status ==
           qmomkit::RawMomentStatus::invalid_input);
    assert(qmomkit::compute_raw_moments(qmomkit::ChiSquaredDistribution{-1.0}, 2).status == qmomkit::RawMomentStatus::invalid_input);
    assert(qmomkit::compute_raw_moments(qmomkit::ChiSquaredDistribution{std::numeric_limits<double>::quiet_NaN()}, 2).status ==
           qmomkit::RawMomentStatus::invalid_input);
    assert(qmomkit::compute_raw_moments(qmomkit::WeibullDistribution{1.0, std::numeric_limits<double>::quiet_NaN()}, 2).status ==
           qmomkit::RawMomentStatus::invalid_input);
    assert(qmomkit::compute_raw_moments(qmomkit::WeibullDistribution{std::numeric_limits<double>::infinity(), 1.0}, 2).status ==
           qmomkit::RawMomentStatus::invalid_input);
}

// Distributions paired with the first order whose ln E[X^k] exceeds
// ln(DBL_MAX) ~ 709.8. The neighbouring orders lie at least ln 2 from that
// bound, so rounding cannot move it:
//   gamma(1, DBL_MAX / 2)  E[X^k] = k! (DBL_MAX / 2)^k                        k = 2
//   chi(3)                 E[X^k] = 2^(k/2) Gamma((3 + k) / 2) / Gamma(3/2)   k = 300
//   chi-squared(2)         E[X^k] = 2^k k!                                    k = 151
//   Weibull(0.1, 2)        E[X^k] = 2^k Gamma(1 + 10 k)                       k = 17
std::vector<std::pair<qmomkit::DistributionSpecification, std::size_t>> overflowing_distributions()
{
    return {
        {qmomkit::GammaDistribution{1.0, std::numeric_limits<double>::max() / 2.0}, 2},
        {qmomkit::ChiDistribution{3.0}, 300},
        {qmomkit::ChiSquaredDistribution{2.0}, 151},
        {qmomkit::WeibullDistribution{0.1, 2.0}, 17},
    };
}

void test_overflow_status()
{
    const qmomkit::RawMomentResult log_normal = qmomkit::compute_raw_moments(qmomkit::LogNormalDistribution{100.0, 10.0}, 20);

    assert(!log_normal.ok());
    assert(log_normal.status == qmomkit::RawMomentStatus::overflow);
    // ln E[X^k] = 100 k + 50 k^2 first exceeds ln(DBL_MAX) ~ 709.8 at k = 3.
    assert(log_normal.failed_order == 3);
    assert(!log_normal.message.empty());

    for (const auto& [distribution, failed_order] : overflowing_distributions())
    {
        const qmomkit::RawMomentResult result = qmomkit::compute_raw_moments(distribution, 2 * failed_order);
        assert(result.status == qmomkit::RawMomentStatus::overflow);
        assert(result.failed_order == failed_order);
        assert(qmomkit::compute_raw_moments(distribution, failed_order - 1).ok());
    }
}

// The selected-order path reports the same first overflowing order, and the
// order just below it is still representable.
void test_selected_overflow_reports_the_failed_order()
{
    for (const auto& [distribution, failed_order] : overflowing_distributions())
    {
        const qmomkit::SelectedRawMomentResult overflowing =
            qmomkit::compute_raw_moments(distribution, qmomkit::RawMomentRequest{{{0}, {failed_order}}});
        assert(overflowing.status == qmomkit::RawMomentStatus::overflow);
        assert_vector_equal(overflowing.failed_order, {failed_order});

        const qmomkit::SelectedRawMomentResult representable =
            qmomkit::compute_raw_moments(distribution, qmomkit::RawMomentRequest{{{0}, {failed_order - 1}}});
        assert(representable.ok());
    }
}

// A moment overflows exactly when it exceeds the largest double and underflows
// exactly when it rounds to zero; a subnormal moment is returned as it is. For
// gamma(1, theta) the first moment is total_weight * theta, which the product
// forms with a single rounding.
void test_product_moments_use_the_whole_double_range()
{
    const double half_maximum = std::numeric_limits<double>::max() / 2.0;
    const qmomkit::RawMomentResult near_maximum = qmomkit::compute_raw_moments(qmomkit::GammaDistribution{1.0, half_maximum}, 1, 1.9);
    assert(near_maximum.ok());
    assert(near_maximum.moments[1] == 1.9 * half_maximum);
    const qmomkit::RawMomentResult above_maximum = qmomkit::compute_raw_moments(qmomkit::GammaDistribution{1.0, half_maximum}, 1, 2.1);
    assert(above_maximum.status == qmomkit::RawMomentStatus::overflow);
    assert(above_maximum.failed_order == 1);

    const qmomkit::RawMomentResult subnormal = qmomkit::compute_raw_moments(qmomkit::GammaDistribution{1.0, 1.0e-300}, 1, 1.0e-10);
    assert(subnormal.ok());
    assert(subnormal.moments[1] > 0.0 && subnormal.moments[1] < std::numeric_limits<double>::min());
    const qmomkit::RawMomentResult below_minimum = qmomkit::compute_raw_moments(qmomkit::GammaDistribution{1.0, 1.0e-300}, 1, 1.0e-30);
    assert(below_minimum.status == qmomkit::RawMomentStatus::underflow);
    assert(below_minimum.failed_order == 1);
}

void test_underflow_status()
{
    const qmomkit::RawMomentResult log_normal = qmomkit::compute_raw_moments(qmomkit::LogNormalDistribution{-1000.0, 0.1}, 1);

    assert(!log_normal.ok());
    assert(log_normal.status == qmomkit::RawMomentStatus::underflow);
    assert(log_normal.failed_order == 1);
    assert(!log_normal.message.empty());

    const double tiny_scale = std::numeric_limits<double>::denorm_min();
    assert(qmomkit::compute_raw_moments(qmomkit::GammaDistribution{1.0, tiny_scale}, 2).status == qmomkit::RawMomentStatus::underflow);
    assert(qmomkit::compute_raw_moments(qmomkit::WeibullDistribution{1.0, tiny_scale}, 2).status == qmomkit::RawMomentStatus::underflow);
}

void test_normal_zero_mean_large_variance_does_not_spuriously_report_non_finite()
{
    // A zero-mean normal distribution has exactly zero odd raw moments. With a large
    // enough (but still finite) standard deviation, evaluating the recurrence in the
    // wrong multiplication order previously produced `0 * inf = NaN` at order 3 instead
    // of the true value 0.0 (see normal_expectations in src/distributions.cpp).
    const double large_standard_deviation = 1.2065474487291362e154;
    const double variance = large_standard_deviation * large_standard_deviation;
    assert(std::isfinite(variance));
    assert(!std::isfinite(2.0 * variance)); // confirms this case exercises the previously-spurious overflow.

    const qmomkit::NormalDistribution distribution{0.0, large_standard_deviation};

    const qmomkit::RawMomentResult all_orders = qmomkit::compute_raw_moments(distribution, 3);
    assert(all_orders.ok());
    assert_near(all_orders.moments[0], 1.0, 0.0);
    assert_near(all_orders.moments[1], 0.0, 0.0);
    assert_near(all_orders.moments[2], variance, 0.0);
    assert_near(all_orders.moments[3], 0.0, 0.0);

    const qmomkit::SelectedRawMomentResult selected = qmomkit::compute_raw_moments(distribution, qmomkit::RawMomentRequest{{{3}}});
    assert(selected.ok());
    assert_near(selected.values[0], 0.0, 0.0);

    // A genuine overflow (the true value truly exceeds double precision) at a higher
    // order must still be reported as such, not masked by the fix above.
    const qmomkit::RawMomentResult genuine_overflow = qmomkit::compute_raw_moments(distribution, 4);
    assert(!genuine_overflow.ok());
    assert(genuine_overflow.status == qmomkit::RawMomentStatus::overflow);
}

// The variance enters the normal moments from order 2 and the log-normal
// moments from order 1, so a variance that overflows must fail only requests
// that reach those orders.
void test_overflowing_variance_fails_only_the_orders_that_use_it()
{
    const qmomkit::NormalDistribution wide_normal{1.0, 1.0e200};
    const qmomkit::RawMomentResult normal_first_order = qmomkit::compute_raw_moments(wide_normal, 1);
    assert(normal_first_order.ok());
    assert_near(normal_first_order.moments[0], 1.0, 0.0);
    assert_near(normal_first_order.moments[1], 1.0, 0.0);
    assert(qmomkit::compute_raw_moments(wide_normal, 0).ok());
    const qmomkit::SelectedRawMomentResult normal_selected = qmomkit::compute_raw_moments(wide_normal, qmomkit::RawMomentRequest{{{1}, {0}}});
    assert(normal_selected.ok());
    assert_near(normal_selected.values[0], 1.0, 0.0);
    assert_near(normal_selected.values[1], 1.0, 0.0);
    const qmomkit::RawMomentResult normal_second_order = qmomkit::compute_raw_moments(wide_normal, 2);
    assert(normal_second_order.status == qmomkit::RawMomentStatus::overflow);
    assert(normal_second_order.failed_order == 2);

    const qmomkit::LogNormalDistribution wide_log_normal{0.0, 1.0e200};
    const qmomkit::RawMomentResult log_normal_zeroth_order = qmomkit::compute_raw_moments(wide_log_normal, 0);
    assert(log_normal_zeroth_order.ok());
    assert_near(log_normal_zeroth_order.moments[0], 1.0, 0.0);
    // Zero total weight returns zeros without evaluating any moment.
    const qmomkit::RawMomentResult log_normal_zero_weight = qmomkit::compute_raw_moments(wide_log_normal, 3, 0.0);
    assert(log_normal_zero_weight.ok());
    for (const double moment : log_normal_zero_weight.moments.values())
    {
        assert(moment == 0.0);
    }
    const qmomkit::SelectedRawMomentResult log_normal_selected_zero_weight =
        qmomkit::compute_raw_moments(wide_log_normal, qmomkit::RawMomentRequest{{{0}, {3}}}, 0.0);
    assert(log_normal_selected_zero_weight.ok());
    assert_near(log_normal_selected_zero_weight.values[0], 0.0, 0.0);
    assert_near(log_normal_selected_zero_weight.values[1], 0.0, 0.0);
    const qmomkit::RawMomentResult log_normal_first_order = qmomkit::compute_raw_moments(wide_log_normal, 2);
    assert(log_normal_first_order.status == qmomkit::RawMomentStatus::overflow);
    assert(log_normal_first_order.failed_order == 1);
    const qmomkit::SelectedRawMomentResult log_normal_selected = qmomkit::compute_raw_moments(wide_log_normal, qmomkit::RawMomentRequest{{{0}, {2}}});
    assert(log_normal_selected.status == qmomkit::RawMomentStatus::overflow);
    assert_vector_equal(log_normal_selected.failed_order, {2});
}

// Even-order normal moments are strictly positive, and odd-order moments are
// nonzero whenever the mean is, so such a moment that rounds to zero is an
// underflow, not an exact zero.
void test_normal_moments_that_underflow_are_reported()
{
    // The variance underflows, so with a zero mean every even moment rounds to zero.
    const qmomkit::NormalDistribution narrow{0.0, 1.0e-200};
    const qmomkit::RawMomentResult narrow_moments = qmomkit::compute_raw_moments(narrow, 4);
    assert(narrow_moments.status == qmomkit::RawMomentStatus::underflow);
    assert(narrow_moments.failed_order == 2);
    const qmomkit::SelectedRawMomentResult narrow_selected = qmomkit::compute_raw_moments(narrow, qmomkit::RawMomentRequest{{{2}}});
    assert(narrow_selected.status == qmomkit::RawMomentStatus::underflow);
    assert_vector_equal(narrow_selected.failed_order, {2});
    // The odd moments of a zero-mean normal are exactly zero, not an underflow.
    const qmomkit::RawMomentResult narrow_first_order = qmomkit::compute_raw_moments(narrow, 1);
    assert(narrow_first_order.ok());
    assert(narrow_first_order.moments[1] == 0.0);
    // With a nonzero mean the same variance is negligible next to the mean, so
    // every moment is representable and exactly the power of the mean.
    const qmomkit::RawMomentResult offset_narrow = qmomkit::compute_raw_moments(qmomkit::NormalDistribution{1.0, 1.0e-200}, 4);
    assert(offset_narrow.ok());
    for (const double moment : offset_narrow.moments.values())
    {
        assert(moment == 1.0);
    }

    // An even-order moment underflows in the recurrence.
    const qmomkit::RawMomentResult even_order = qmomkit::compute_raw_moments(qmomkit::NormalDistribution{0.0, 1.0e-100}, 4);
    assert(even_order.status == qmomkit::RawMomentStatus::underflow);
    assert(even_order.failed_order == 4);

    // An odd-order moment of a nonzero-mean normal underflows in the recurrence.
    const qmomkit::RawMomentResult odd_order = qmomkit::compute_raw_moments(qmomkit::NormalDistribution{1.0e-160, 1.0e-160}, 3);
    assert(odd_order.status == qmomkit::RawMomentStatus::underflow);
    assert(odd_order.failed_order == 3);

    // The weighting underflows.
    const qmomkit::RawMomentResult weighted = qmomkit::compute_raw_moments(qmomkit::NormalDistribution{0.0, 1.0e-100}, 2, 1.0e-300);
    assert(weighted.status == qmomkit::RawMomentStatus::underflow);
    assert(weighted.failed_order == 2);
    const qmomkit::SelectedRawMomentResult weighted_selected =
        qmomkit::compute_raw_moments(qmomkit::NormalDistribution{0.0, 1.0e-100}, qmomkit::RawMomentRequest{{{2}}}, 1.0e-300);
    assert(weighted_selected.status == qmomkit::RawMomentStatus::underflow);
    assert_vector_equal(weighted_selected.failed_order, {2});
    const qmomkit::RawMomentResult weighted_mean = qmomkit::compute_raw_moments(qmomkit::NormalDistribution{1.0e-200, 1.0}, 1, 1.0e-200);
    assert(weighted_mean.status == qmomkit::RawMomentStatus::underflow);
    assert(weighted_mean.failed_order == 1);
}

void test_empty_mixture_is_invalid()
{
    const qmomkit::RawMomentResult result = qmomkit::compute_raw_moments(qmomkit::DistributionMixture{}, 2);

    assert(!result.ok());
    assert(result.status == qmomkit::RawMomentStatus::invalid_input);
    assert(!result.message.empty());
}

void test_mixture_maximum_order_zero()
{
    const qmomkit::DistributionMixture mixture{{{qmomkit::GammaDistribution{2.0, 3.0}, 2.0}, {qmomkit::BetaDistribution{1.0, 1.0}, 0.5}}};

    const qmomkit::RawMomentResult result = qmomkit::compute_raw_moments(mixture, 0);

    assert(result.ok());
    assert(result.moments.size() == 1);
    assert_near(result.moments[0], 2.5, 0.0);
}

void test_one_component_mixture_matches_distribution_result()
{
    const qmomkit::DistributionMixture mixture{{{qmomkit::GammaDistribution{2.0, 3.0}, 4.0}}};
    const qmomkit::RawMomentResult mixture_result = qmomkit::compute_raw_moments(mixture, 3);
    const qmomkit::RawMomentResult distribution_result = qmomkit::compute_raw_moments(qmomkit::GammaDistribution{2.0, 3.0}, 3, 4.0);

    assert(mixture_result.ok());
    assert(distribution_result.ok());
    assert(mixture_result.moments.size() == distribution_result.moments.size());
    for (std::size_t order = 0; order < mixture_result.moments.size(); ++order)
    {
        assert_near(mixture_result.moments[order], distribution_result.moments[order], 1.0e-13);
    }
}

void test_heterogeneous_mixture_dispatch()
{
    const qmomkit::DistributionMixture mixture{{{qmomkit::GammaDistribution{2.0, 3.0}, 2.0}, {qmomkit::BetaDistribution{2.0, 3.0}, 5.0}}};
    const qmomkit::RawMomentResult result = qmomkit::compute_raw_moments(mixture, 1);

    assert(result.ok());
    assert_near(result.moments[0], 7.0, 1.0e-14);
    assert_near(result.moments[1], 2.0 * 6.0 + 5.0 * 0.4, 1.0e-13);
}

void test_zero_density_component_is_skipped()
{
    const qmomkit::DistributionMixture mixture{{{qmomkit::NormalDistribution{0.0, 0.0}, 0.0}, {qmomkit::GammaDistribution{1.0, 2.0}, 3.0}}};
    const qmomkit::RawMomentResult result = qmomkit::compute_raw_moments(mixture, 2);

    assert(result.ok());
    assert_near(result.moments[0], 3.0, 0.0);
    assert_near(result.moments[1], 6.0, 1.0e-14);
    assert_near(result.moments[2], 24.0, 1.0e-13);
}

void test_zero_density_extreme_component_is_skipped()
{
    const qmomkit::DistributionMixture mixture{
        {{qmomkit::NormalDistribution{std::numeric_limits<double>::max(), 1.0}, 0.0}, {qmomkit::BetaDistribution{1.0, 1.0}, 2.0}}};
    const qmomkit::RawMomentResult result = qmomkit::compute_raw_moments(mixture, 2);

    assert(result.ok());
    assert_near(result.moments[0], 2.0, 0.0);
    assert_near(result.moments[1], 1.0, 1.0e-14);
    assert_near(result.moments[2], 2.0 / 3.0, 1.0e-14);
}

void test_invalid_mixture_densities()
{
    const qmomkit::DistributionMixture negative_density{{{qmomkit::NormalDistribution{0.0, 1.0}, -1.0}}};
    const qmomkit::RawMomentResult negative_result = qmomkit::compute_raw_moments(negative_density, 2);
    assert(!negative_result.ok());
    assert(negative_result.status == qmomkit::RawMomentStatus::invalid_input);
    assert(negative_result.failed_component == 0);

    const qmomkit::DistributionMixture nan_density{{{qmomkit::NormalDistribution{0.0, 1.0}, std::numeric_limits<double>::quiet_NaN()}}};
    assert(qmomkit::compute_raw_moments(nan_density, 2).status == qmomkit::RawMomentStatus::invalid_input);

    const qmomkit::DistributionMixture infinite_density{{{qmomkit::NormalDistribution{0.0, 1.0}, std::numeric_limits<double>::infinity()}}};
    assert(qmomkit::compute_raw_moments(infinite_density, 2).status == qmomkit::RawMomentStatus::invalid_input);
}

void test_all_zero_density_mixture_is_invalid()
{
    const qmomkit::DistributionMixture mixture{{{qmomkit::NormalDistribution{0.0, 0.0}, 0.0}, {qmomkit::LogNormalDistribution{0.0, 0.0}, 0.0}}};
    const qmomkit::RawMomentResult result = qmomkit::compute_raw_moments(mixture, 2);

    assert(!result.ok());
    assert(result.status == qmomkit::RawMomentStatus::invalid_input);
}

void test_positive_component_failure_reports_index()
{
    const qmomkit::DistributionMixture mixture{{{qmomkit::NormalDistribution{0.0, 1.0}, 1.0}, {qmomkit::NormalDistribution{0.0, 0.0}, 2.0}}};
    const qmomkit::RawMomentResult result = qmomkit::compute_raw_moments(mixture, 2);

    assert(!result.ok());
    assert(result.status == qmomkit::RawMomentStatus::invalid_input);
    assert(result.failed_component == 1);
    assert(!result.message.empty());
}

void test_surrounded_positive_component_failure_reports_index()
{
    const qmomkit::DistributionMixture mixture{
        {{qmomkit::GammaDistribution{1.0, 1.0}, 2.0}, {qmomkit::LogNormalDistribution{0.0, 0.0}, 3.0}, {qmomkit::BetaDistribution{1.0, 1.0}, 4.0}}};
    const qmomkit::RawMomentResult result = qmomkit::compute_raw_moments(mixture, 2);

    assert(!result.ok());
    assert(result.status == qmomkit::RawMomentStatus::invalid_input);
    assert(result.failed_component == 1);
    assert(!result.message.empty());
}

void test_mixture_accumulation_overflow_reports_index()
{
    const double large_mean = 0.75 * std::numeric_limits<double>::max();
    const qmomkit::DistributionMixture mixture{
        {{qmomkit::NormalDistribution{large_mean, 1.0}, 1.0}, {qmomkit::NormalDistribution{large_mean, 1.0}, 1.0}}};
    const qmomkit::RawMomentResult result = qmomkit::compute_raw_moments(mixture, 1);

    assert(!result.ok());
    assert(result.status == qmomkit::RawMomentStatus::overflow);
    assert(result.failed_order == 1);
    assert(result.failed_component == 1);
}

void assert_selected_mixture_matches_all_orders(const qmomkit::DistributionMixture& mixture, const qmomkit::RawMomentRequest& request,
                                                const std::size_t maximum_order, const double tolerance)
{
    const qmomkit::RawMomentResult all_orders = qmomkit::compute_raw_moments(mixture, maximum_order);
    const qmomkit::SelectedRawMomentResult selected = qmomkit::compute_raw_moments(mixture, request);

    assert(all_orders.ok());
    assert(selected.ok());
    assert(selected.indices == request.indices);
    assert(selected.values.size() == request.indices.size());
    for (std::size_t index = 0; index < request.indices.size(); ++index)
    {
        assert(request.indices[index].size() == 1);
        assert_near(selected.values[index], all_orders.moments[request.indices[index][0]], tolerance);
    }
}

void test_selected_one_component_mixture_matches_all_orders()
{
    const qmomkit::DistributionMixture mixture{{{qmomkit::GammaDistribution{2.0, 3.0}, 4.0}}};
    assert_selected_mixture_matches_all_orders(mixture, {{{0}, {2}, {4}}}, 4, 1.0e-12);
}

void test_selected_heterogeneous_mixture_matches_all_orders()
{
    const qmomkit::DistributionMixture mixture{
        {{qmomkit::GammaDistribution{2.0, 3.0}, 2.0}, {qmomkit::BetaDistribution{2.0, 3.0}, 5.0}, {qmomkit::NormalDistribution{1.0, 0.5}, 1.5}}};
    assert_selected_mixture_matches_all_orders(mixture, {{{0}, {1}, {3}}}, 3, 1.0e-12);
}

void test_selected_mixture_preserves_reordered_and_duplicate_requests()
{
    const qmomkit::DistributionMixture mixture{{{qmomkit::GammaDistribution{1.5, 0.5}, 2.0}, {qmomkit::ChiSquaredDistribution{4.0}, 0.75}}};
    assert_selected_mixture_matches_all_orders(mixture, {{{4}, {0}, {2}, {2}, {1}}}, 4, 1.0e-11);
}

void test_selected_mixture_zero_density_component_is_skipped()
{
    const qmomkit::DistributionMixture mixture{{{qmomkit::NormalDistribution{0.0, 0.0}, 0.0}, {qmomkit::GammaDistribution{1.0, 2.0}, 3.0}}};
    const qmomkit::RawMomentRequest request{{{2}, {0}, {1}}};
    const qmomkit::SelectedRawMomentResult result = qmomkit::compute_raw_moments(mixture, request);

    assert(result.ok());
    assert(result.indices == request.indices);
    assert_near(result.values[0], 24.0, 1.0e-13);
    assert_near(result.values[1], 3.0, 0.0);
    assert_near(result.values[2], 6.0, 1.0e-14);
}

void test_selected_mixture_empty_request_is_valid_for_valid_mixture()
{
    const qmomkit::DistributionMixture mixture{{{qmomkit::GammaDistribution{1.0, 2.0}, 3.0}}};
    const qmomkit::SelectedRawMomentResult result = qmomkit::compute_raw_moments(mixture, qmomkit::RawMomentRequest{});

    assert(result.ok());
    assert(result.indices.empty());
    assert(result.values.empty());
}

void test_selected_mixture_invalid_component_failure_reports_index()
{
    const qmomkit::DistributionMixture mixture{{{qmomkit::GammaDistribution{1.0, 1.0}, 2.0}, {qmomkit::LogNormalDistribution{0.0, 0.0}, 3.0}}};
    const qmomkit::SelectedRawMomentResult result = qmomkit::compute_raw_moments(mixture, {{{0}, {1}}});

    assert(!result.ok());
    assert(result.status == qmomkit::RawMomentStatus::invalid_input);
    assert(result.failed_component == 1);
    assert(result.failed_dimension == 0);
    assert_vector_equal(result.failed_order, {0});
    assert(!result.message.empty());
}

void test_selected_mixture_overflow_propagates_component_and_order()
{
    const qmomkit::DistributionMixture mixture{{{qmomkit::GammaDistribution{1.0, 1.0}, 2.0}, {qmomkit::LogNormalDistribution{100.0, 10.0}, 3.0}}};
    const qmomkit::SelectedRawMomentResult result = qmomkit::compute_raw_moments(mixture, {{{20}}});

    assert(!result.ok());
    assert(result.status == qmomkit::RawMomentStatus::overflow);
    assert(result.failed_component == 1);
    assert(result.failed_dimension == 0);
    assert_vector_equal(result.failed_order, {20});
}

void test_selected_mixture_underflow_propagates_component_and_order()
{
    const qmomkit::DistributionMixture mixture{{{qmomkit::LogNormalDistribution{-1000.0, 0.1}, 3.0}}};
    const qmomkit::SelectedRawMomentResult result = qmomkit::compute_raw_moments(mixture, {{{1}}});

    assert(!result.ok());
    assert(result.status == qmomkit::RawMomentStatus::underflow);
    assert(result.failed_component == 0);
    assert(result.failed_dimension == 0);
    assert_vector_equal(result.failed_order, {1});
}

void test_selected_mixture_accumulation_overflow_reports_index()
{
    const double large_mean = 0.75 * std::numeric_limits<double>::max();
    const qmomkit::DistributionMixture mixture{
        {{qmomkit::NormalDistribution{large_mean, 1.0}, 1.0}, {qmomkit::NormalDistribution{large_mean, 1.0}, 1.0}}};
    const qmomkit::SelectedRawMomentResult result = qmomkit::compute_raw_moments(mixture, {{{1}}});

    assert(!result.ok());
    assert(result.status == qmomkit::RawMomentStatus::overflow);
    assert(result.failed_component == 1);
    assert(result.failed_dimension == 0);
    assert_vector_equal(result.failed_order, {1});
}

void test_independent_product_zero_dimension_is_invalid()
{
    const qmomkit::TensorMomentResult result = qmomkit::compute_raw_moments(qmomkit::IndependentProductDistribution{}, std::vector<std::size_t>{});

    assert(!result.ok());
    assert(result.status == qmomkit::RawMomentStatus::invalid_input);
    assert(!result.message.empty());
}

void test_independent_product_rank_mismatch_is_invalid()
{
    const qmomkit::IndependentProductDistribution product{{qmomkit::NormalDistribution{0.0, 1.0}}, 1.0};
    const qmomkit::TensorMomentResult result = qmomkit::compute_raw_moments(product, {1, 1});

    assert(!result.ok());
    assert(result.status == qmomkit::RawMomentStatus::invalid_input);
}

void test_independent_product_maximum_order_overflow_guard()
{
    const qmomkit::IndependentProductDistribution product{{qmomkit::NormalDistribution{0.0, 1.0}}, 1.0};
    const qmomkit::TensorMomentResult result = qmomkit::compute_raw_moments(product, {std::numeric_limits<std::size_t>::max()});

    assert(!result.ok());
    assert(result.status == qmomkit::RawMomentStatus::invalid_input);
    assert(result.failed_dimension == 0);
    assert_vector_equal(result.failed_order, {std::numeric_limits<std::size_t>::max()});

    // Each maximum order is valid on its own, but the full tensor of raw
    // moments they span has more entries than std::size_t can count.
    const qmomkit::IndependentProductDistribution product_2d{{qmomkit::NormalDistribution{0.0, 1.0}, qmomkit::NormalDistribution{0.0, 1.0}}, 1.0};
    const qmomkit::TensorMomentResult tensor_overflow = qmomkit::compute_raw_moments(product_2d, {std::numeric_limits<std::size_t>::max() - 1, 1});
    assert(!tensor_overflow.ok());
    assert(tensor_overflow.status == qmomkit::RawMomentStatus::invalid_input);
    assert(tensor_overflow.message == "The extents of the full tensor of raw moments are invalid: Index extents overflow size_t.");
}

void test_independent_product_invalid_global_number_density()
{
    const qmomkit::IndependentProductDistribution product{{qmomkit::NormalDistribution{0.0, 1.0}}, -1.0};
    assert(qmomkit::compute_raw_moments(product, {1}).status == qmomkit::RawMomentStatus::invalid_input);

    const qmomkit::IndependentProductDistribution nan_product{{qmomkit::NormalDistribution{0.0, 1.0}}, std::numeric_limits<double>::quiet_NaN()};
    assert(qmomkit::compute_raw_moments(nan_product, {1}).status == qmomkit::RawMomentStatus::invalid_input);

    const qmomkit::IndependentProductDistribution infinite_product{{qmomkit::NormalDistribution{0.0, 1.0}}, std::numeric_limits<double>::infinity()};
    assert(qmomkit::compute_raw_moments(infinite_product, {1}).status == qmomkit::RawMomentStatus::invalid_input);
}

void test_independent_product_invalid_marginal_reports_dimension()
{
    const qmomkit::IndependentProductDistribution product{{qmomkit::NormalDistribution{0.0, 1.0}, qmomkit::GammaDistribution{0.0, 1.0}}, 2.0};
    const qmomkit::TensorMomentResult result = qmomkit::compute_raw_moments(product, {1, 1});

    assert(!result.ok());
    assert(result.status == qmomkit::RawMomentStatus::invalid_input);
    assert(result.failed_dimension == 1);
    assert_vector_equal(result.failed_order, {0, 0});
    assert(!result.message.empty());
}

void test_independent_product_marginal_failure_reports_order()
{
    const qmomkit::IndependentProductDistribution product{{qmomkit::NormalDistribution{0.0, 1.0}, qmomkit::LogNormalDistribution{100.0, 10.0}}, 1.0};
    const qmomkit::TensorMomentResult result = qmomkit::compute_raw_moments(product, {1, 20});

    assert(!result.ok());
    assert(result.status == qmomkit::RawMomentStatus::overflow);
    assert(result.failed_dimension == 1);
    assert(result.failed_order.size() == 2);
    assert(result.failed_order[0] == 0);
    assert(result.failed_order[1] > 0);
}

void test_independent_product_maximum_orders_zero()
{
    const qmomkit::IndependentProductDistribution product{{qmomkit::NormalDistribution{10.0, 2.0}, qmomkit::GammaDistribution{3.0, 4.0}}, 7.5};
    const qmomkit::TensorMomentResult result = qmomkit::compute_raw_moments(product, {0, 0});

    assert(result.ok());
    assert_vector_equal(result.moments.extents, {1, 1});
    assert(result.moments.values.size() == 1);
    assert_near(result.moments.values[0], 7.5, 0.0);
}

void test_independent_product_zero_global_density_skips_marginals()
{
    const qmomkit::IndependentProductDistribution product{{qmomkit::NormalDistribution{0.0, 0.0}, qmomkit::GammaDistribution{0.0, 1.0}}, 0.0};
    const qmomkit::TensorMomentResult result = qmomkit::compute_raw_moments(product, {2, 3});

    assert(result.ok());
    assert_vector_equal(result.moments.extents, {3, 4});
    assert(result.moments.values.size() == 12);
    for (double value : result.moments.values)
    {
        assert_near(value, 0.0, 0.0);
    }
}

void test_independent_product_zero_global_density_still_validates_extents()
{
    const std::size_t large_order = std::numeric_limits<std::size_t>::max() / 2;
    const qmomkit::IndependentProductDistribution product{{qmomkit::NormalDistribution{0.0, 0.0}, qmomkit::GammaDistribution{0.0, 1.0}}, 0.0};
    const qmomkit::TensorMomentResult result = qmomkit::compute_raw_moments(product, {large_order, 1});

    assert(!result.ok());
    assert(result.status == qmomkit::RawMomentStatus::invalid_input);
    assert(!result.message.empty());
}

void test_independent_product_tensor_size_overflow()
{
    const std::size_t large_order = std::numeric_limits<std::size_t>::max() / 2;
    const qmomkit::IndependentProductDistribution product{{qmomkit::NormalDistribution{0.0, 1.0}, qmomkit::NormalDistribution{0.0, 1.0}}, 1.0};
    const qmomkit::TensorMomentResult result = qmomkit::compute_raw_moments(product, {large_order, 1});

    assert(!result.ok());
    assert(result.status == qmomkit::RawMomentStatus::invalid_input);
    assert(!result.message.empty());
}

void test_independent_product_row_major_flattening_spots()
{
    const qmomkit::IndependentProductDistribution product{{qmomkit::NormalDistribution{2.0, 1.0}, qmomkit::GammaDistribution{1.0, 1.0}}, 3.0};
    const qmomkit::TensorMomentResult result = qmomkit::compute_raw_moments(product, {1, 2});

    assert(result.ok());
    assert_vector_equal(result.moments.extents, {2, 3});
    assert(result.moments.values.size() == 6);
    assert_near(result.moments.values[qmomkit::flatten_multi_index({0, 0}, result.moments.extents)], 3.0, 0.0);
    assert_near(result.moments.values[qmomkit::flatten_multi_index({0, 1}, result.moments.extents)], 3.0, 1.0e-14);
    assert_near(result.moments.values[qmomkit::flatten_multi_index({1, 0}, result.moments.extents)], 6.0, 1.0e-14);
    assert_near(result.moments.values[qmomkit::flatten_multi_index({1, 2}, result.moments.extents)], 12.0, 1.0e-14);
}

void test_independent_product_mixture_marginal_is_normalized()
{
    const qmomkit::NormalDistribution normal{-0.4, 0.7};
    const qmomkit::DistributionMixture mixture{
        {{qmomkit::LogNormalDistribution{-0.3, 0.25}, 1.5}, {qmomkit::LogNormalDistribution{0.4, 0.45}, 3.25}}};
    const double number_density = 2.75;
    const qmomkit::IndependentProductDistribution product{{qmomkit::DistributionSpecification{normal}, mixture}, number_density};

    const qmomkit::RawMomentResult normal_moments = qmomkit::compute_raw_moments(normal, 1, 1.0);
    const qmomkit::RawMomentResult mixture_moments = qmomkit::compute_raw_moments(mixture, 2);
    const qmomkit::TensorMomentResult result = qmomkit::compute_raw_moments(product, {1, 2});

    assert(normal_moments.ok());
    assert(mixture_moments.ok());
    assert(result.ok());
    assert_near(result.moments.values[qmomkit::flatten_multi_index({0, 0}, result.moments.extents)], number_density, 0.0);
    assert_near(result.moments.values[qmomkit::flatten_multi_index({1, 0}, result.moments.extents)],
                number_density * normal_moments.moments[1] / normal_moments.moments[0], 1.0e-14);
    assert_near(result.moments.values[qmomkit::flatten_multi_index({0, 1}, result.moments.extents)],
                number_density * mixture_moments.moments[1] / mixture_moments.moments[0], 1.0e-14);
    assert_near(result.moments.values[qmomkit::flatten_multi_index({1, 1}, result.moments.extents)],
                number_density * normal_moments.moments[1] / normal_moments.moments[0] * mixture_moments.moments[1] / mixture_moments.moments[0],
                1.0e-14);
    assert_near(result.moments.values[qmomkit::flatten_multi_index({0, 2}, result.moments.extents)],
                number_density * mixture_moments.moments[2] / mixture_moments.moments[0], 1.0e-14);
}

void assert_selected_product_matches_moment_tensor(const qmomkit::IndependentProductDistribution& product, const qmomkit::RawMomentRequest& request,
                                                   const std::vector<std::size_t>& maximum_orders, const double tolerance)
{
    const qmomkit::TensorMomentResult tensor_result = qmomkit::compute_raw_moments(product, maximum_orders);
    const qmomkit::SelectedRawMomentResult selected = qmomkit::compute_raw_moments(product, request);

    assert(tensor_result.ok());
    assert(selected.ok());
    assert(selected.indices == request.indices);
    assert(selected.values.size() == request.indices.size());
    for (std::size_t index = 0; index < request.indices.size(); ++index)
    {
        assert_near(selected.values[index],
                    tensor_result.moments.values[qmomkit::flatten_multi_index(request.indices[index], tensor_result.moments.extents)], tolerance);
    }
}

void test_selected_independent_product_one_dimensional_matches_moment_tensor()
{
    const qmomkit::IndependentProductDistribution product{{qmomkit::GammaDistribution{1.5, 2.0}}, 2.5};
    assert_selected_product_matches_moment_tensor(product, {{{0}, {3}, {1}}}, {3}, 1.0e-12);
}

void test_selected_independent_product_two_dimensional_nonrectangular_matches_moment_tensor()
{
    const qmomkit::IndependentProductDistribution product{{qmomkit::NormalDistribution{2.0, 1.0}, qmomkit::GammaDistribution{1.0, 1.0}}, 3.0};
    assert_selected_product_matches_moment_tensor(product, {{{0, 0}, {2, 1}, {1, 3}}}, {2, 3}, 1.0e-12);
}

void test_selected_independent_product_mixture_marginal_matches_moment_tensor()
{
    const qmomkit::DistributionMixture mixture{
        {{qmomkit::LogNormalDistribution{-0.3, 0.25}, 1.5}, {qmomkit::LogNormalDistribution{0.4, 0.45}, 3.25}}};
    const qmomkit::IndependentProductDistribution product{{qmomkit::DistributionSpecification{qmomkit::NormalDistribution{-0.4, 0.7}}, mixture},
                                                          2.75};
    assert_selected_product_matches_moment_tensor(product, {{{1, 0}, {0, 1}, {1, 1}, {0, 2}}}, {1, 2}, 1.0e-12);
}

void test_selected_independent_product_three_dimensional_matches_moment_tensor()
{
    const qmomkit::IndependentProductDistribution product{
        {qmomkit::NormalDistribution{1.0, 0.5}, qmomkit::BetaDistribution{2.0, 3.0}, qmomkit::ChiSquaredDistribution{4.0}}, 1.25};
    assert_selected_product_matches_moment_tensor(product, {{{0, 0, 0}, {2, 1, 3}, {1, 2, 0}}}, {2, 2, 3}, 1.0e-10);
}

void test_selected_independent_product_preserves_reordered_and_duplicate_indices()
{
    const qmomkit::IndependentProductDistribution product{{qmomkit::NormalDistribution{2.0, 1.0}, qmomkit::GammaDistribution{1.0, 1.0}}, 3.0};
    assert_selected_product_matches_moment_tensor(product, {{{1, 2}, {0, 0}, {1, 2}, {0, 1}}}, {1, 2}, 1.0e-12);
}

void test_selected_independent_product_zero_global_density_skips_marginals()
{
    const qmomkit::IndependentProductDistribution product{{qmomkit::NormalDistribution{0.0, 0.0}, qmomkit::GammaDistribution{0.0, 1.0}}, 0.0};
    const qmomkit::RawMomentRequest request{{{2, 3}, {0, 0}}};
    const qmomkit::SelectedRawMomentResult result = qmomkit::compute_raw_moments(product, request);

    assert(result.ok());
    assert(result.indices == request.indices);
    assert(result.values.size() == request.indices.size());
    for (double value : result.values)
    {
        assert_near(value, 0.0, 0.0);
    }
}

void test_selected_independent_product_empty_request_is_valid()
{
    const qmomkit::IndependentProductDistribution product{{qmomkit::NormalDistribution{0.0, 1.0}, qmomkit::GammaDistribution{1.0, 2.0}}, 3.0};
    const qmomkit::SelectedRawMomentResult result = qmomkit::compute_raw_moments(product, qmomkit::RawMomentRequest{});

    assert(result.ok());
    assert(result.indices.empty());
    assert(result.values.empty());
}

void test_selected_independent_product_invalid_rank_reports_order()
{
    const qmomkit::IndependentProductDistribution product{{qmomkit::NormalDistribution{0.0, 1.0}, qmomkit::GammaDistribution{1.0, 2.0}}, 3.0};
    const qmomkit::SelectedRawMomentResult result = qmomkit::compute_raw_moments(product, qmomkit::RawMomentRequest{{{0}}});

    assert(!result.ok());
    assert(result.status == qmomkit::RawMomentStatus::invalid_input);
    assert(result.failed_dimension == qmomkit::invalid_dimension_index);
    assert_vector_equal(result.failed_order, {0});
}

void test_selected_independent_product_oversized_order_reports_dimension()
{
    const qmomkit::IndependentProductDistribution product{{qmomkit::NormalDistribution{0.0, 1.0}, qmomkit::GammaDistribution{1.0, 2.0}}, 3.0};
    const qmomkit::SelectedRawMomentResult result = qmomkit::compute_raw_moments(product, {{{1, std::numeric_limits<std::size_t>::max()}}});

    assert(!result.ok());
    assert(result.status == qmomkit::RawMomentStatus::invalid_input);
    assert(result.failed_dimension == 1);
    assert_vector_equal(result.failed_order, {1, std::numeric_limits<std::size_t>::max()});
}

void test_selected_independent_product_marginal_failure_reports_dimension()
{
    const qmomkit::IndependentProductDistribution product{{qmomkit::NormalDistribution{0.0, 1.0}, qmomkit::GammaDistribution{0.0, 1.0}}, 2.0};
    const qmomkit::SelectedRawMomentResult result = qmomkit::compute_raw_moments(product, {{{1, 2}}});

    assert(!result.ok());
    assert(result.status == qmomkit::RawMomentStatus::invalid_input);
    assert(result.failed_dimension == 1);
    assert_vector_equal(result.failed_order, {1, 2});
    assert(!result.message.empty());
}

void test_selected_independent_product_marginal_overflow_reports_full_order()
{
    const qmomkit::IndependentProductDistribution product{{qmomkit::NormalDistribution{0.0, 1.0}, qmomkit::LogNormalDistribution{100.0, 10.0}}, 1.0};
    const qmomkit::SelectedRawMomentResult result = qmomkit::compute_raw_moments(product, {{{0, 20}}});

    assert(!result.ok());
    assert(result.status == qmomkit::RawMomentStatus::overflow);
    assert(result.failed_dimension == 1);
    assert_vector_equal(result.failed_order, {0, 20});
}

void test_selected_independent_product_product_overflow_reports_full_order()
{
    const qmomkit::IndependentProductDistribution product{
        {qmomkit::NormalDistribution{0.5 * std::numeric_limits<double>::max(), 1.0}, qmomkit::NormalDistribution{3.0, 1.0}}, 1.0};
    const qmomkit::SelectedRawMomentResult result = qmomkit::compute_raw_moments(product, {{{1, 1}}});

    assert(!result.ok());
    assert(result.status == qmomkit::RawMomentStatus::overflow);
    assert(result.failed_dimension == qmomkit::invalid_dimension_index);
    assert_vector_equal(result.failed_order, {1, 1});
}

void test_selected_independent_product_product_underflow_reports_full_order()
{
    const qmomkit::IndependentProductDistribution product{{qmomkit::NormalDistribution{std::numeric_limits<double>::denorm_min(), 1.0},
                                                           qmomkit::NormalDistribution{std::numeric_limits<double>::denorm_min(), 1.0}},
                                                          1.0};
    const qmomkit::SelectedRawMomentResult result = qmomkit::compute_raw_moments(product, {{{1, 1}}});

    assert(!result.ok());
    assert(result.status == qmomkit::RawMomentStatus::underflow);
    assert(result.failed_dimension == qmomkit::invalid_dimension_index);
    assert_vector_equal(result.failed_order, {1, 1});
}

void test_selected_independent_product_matches_cqmom_layout_subselection()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::IndependentProductDistribution product{{qmomkit::GammaDistribution{1.5, 0.5}, qmomkit::BetaDistribution{2.0, 3.0}}, 2.0};

    const qmomkit::TensorMomentResult tensor_result = qmomkit::compute_raw_moments(product, {3, 3});
    const qmomkit::CqmomMomentData required_moment_data = qmomkit::extract_cqmom_moments(tensor_result.moments, layout);
    const qmomkit::SelectedRawMomentResult selected =
        qmomkit::compute_raw_moments(product, qmomkit::RawMomentRequest{layout.required_indices.to_multi_indices()});

    assert(tensor_result.ok());
    assert(selected.ok());
    assert(selected.indices == layout.required_indices.to_multi_indices());
    assert(selected.values.size() == required_moment_data.values.size());
    for (std::size_t index = 0; index < selected.values.size(); ++index)
    {
        assert_near(selected.values[index], required_moment_data.values[index], 1.0e-12);
    }
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"selected_raw_moment_result_preserves_request_order", test_selected_raw_moment_result_preserves_request_order},
            {"selected_raw_moment_result_accepts_noncontiguous_multidimensional_indices",
             test_selected_raw_moment_result_accepts_noncontiguous_multidimensional_indices},
            {"selected_raw_moment_result_preserves_duplicate_indices", test_selected_raw_moment_result_preserves_duplicate_indices},
            {"selected_raw_moment_empty_request_is_valid", test_selected_raw_moment_empty_request_is_valid},
            {"selected_raw_moment_invalid_requests", test_selected_raw_moment_invalid_requests},
            {"selected_distribution_contiguous_orders_match_all_orders", test_selected_distribution_contiguous_orders_match_all_orders},
            {"selected_distribution_noncontiguous_orders_match_all_orders", test_selected_distribution_noncontiguous_orders_match_all_orders},
            {"selected_distribution_reordered_and_duplicate_orders_match_all_orders",
             test_selected_distribution_reordered_and_duplicate_orders_match_all_orders},
            {"selected_distribution_high_order_gaps_match_all_orders", test_selected_distribution_high_order_gaps_match_all_orders},
            {"selected_distribution_zero_total_weight", test_selected_distribution_zero_total_weight},
            {"selected_distribution_invalid_request_rank", test_selected_distribution_invalid_request_rank},
            {"selected_distribution_overflow_and_underflow_failures", test_selected_distribution_overflow_and_underflow_failures},
            {"variant_dispatch_and_nonunit_weight", test_variant_dispatch_and_nonunit_weight},
            {"maximum_order_zero", test_maximum_order_zero},
            {"zero_total_weight_returns_zero_sequence", test_zero_total_weight_returns_zero_sequence},
            {"invalid_total_weight", test_invalid_total_weight},
            {"invalid_distribution_parameters", test_invalid_distribution_parameters},
            {"overflow_status", test_overflow_status},
            {"selected_overflow_reports_the_failed_order", test_selected_overflow_reports_the_failed_order},
            {"product_moments_use_the_whole_double_range", test_product_moments_use_the_whole_double_range},
            {"overflowing_variance_fails_only_the_orders_that_use_it", test_overflowing_variance_fails_only_the_orders_that_use_it},
            {"normal_moments_that_underflow_are_reported", test_normal_moments_that_underflow_are_reported},
            {"underflow_status", test_underflow_status},
            {"normal_zero_mean_large_variance_does_not_spuriously_report_non_finite",
             test_normal_zero_mean_large_variance_does_not_spuriously_report_non_finite},
            {"empty_mixture_is_invalid", test_empty_mixture_is_invalid},
            {"mixture_maximum_order_zero", test_mixture_maximum_order_zero},
            {"one_component_mixture_matches_distribution_result", test_one_component_mixture_matches_distribution_result},
            {"heterogeneous_mixture_dispatch", test_heterogeneous_mixture_dispatch},
            {"zero_density_component_is_skipped", test_zero_density_component_is_skipped},
            {"zero_density_extreme_component_is_skipped", test_zero_density_extreme_component_is_skipped},
            {"invalid_mixture_densities", test_invalid_mixture_densities},
            {"all_zero_density_mixture_is_invalid", test_all_zero_density_mixture_is_invalid},
            {"positive_component_failure_reports_index", test_positive_component_failure_reports_index},
            {"surrounded_positive_component_failure_reports_index", test_surrounded_positive_component_failure_reports_index},
            {"mixture_accumulation_overflow_reports_index", test_mixture_accumulation_overflow_reports_index},
            {"selected_one_component_mixture_matches_all_orders", test_selected_one_component_mixture_matches_all_orders},
            {"selected_heterogeneous_mixture_matches_all_orders", test_selected_heterogeneous_mixture_matches_all_orders},
            {"selected_mixture_preserves_reordered_and_duplicate_requests", test_selected_mixture_preserves_reordered_and_duplicate_requests},
            {"selected_mixture_zero_density_component_is_skipped", test_selected_mixture_zero_density_component_is_skipped},
            {"selected_mixture_empty_request_is_valid_for_valid_mixture", test_selected_mixture_empty_request_is_valid_for_valid_mixture},
            {"selected_mixture_invalid_component_failure_reports_index", test_selected_mixture_invalid_component_failure_reports_index},
            {"selected_mixture_overflow_propagates_component_and_order", test_selected_mixture_overflow_propagates_component_and_order},
            {"selected_mixture_underflow_propagates_component_and_order", test_selected_mixture_underflow_propagates_component_and_order},
            {"selected_mixture_accumulation_overflow_reports_index", test_selected_mixture_accumulation_overflow_reports_index},
            {"independent_product_zero_dimension_is_invalid", test_independent_product_zero_dimension_is_invalid},
            {"independent_product_rank_mismatch_is_invalid", test_independent_product_rank_mismatch_is_invalid},
            {"independent_product_maximum_order_overflow_guard", test_independent_product_maximum_order_overflow_guard},
            {"independent_product_invalid_global_number_density", test_independent_product_invalid_global_number_density},
            {"independent_product_invalid_marginal_reports_dimension", test_independent_product_invalid_marginal_reports_dimension},
            {"independent_product_marginal_failure_reports_order", test_independent_product_marginal_failure_reports_order},
            {"independent_product_maximum_orders_zero", test_independent_product_maximum_orders_zero},
            {"independent_product_zero_global_density_skips_marginals", test_independent_product_zero_global_density_skips_marginals},
            {"independent_product_zero_global_density_still_validates_extents", test_independent_product_zero_global_density_still_validates_extents},
            {"independent_product_tensor_size_overflow", test_independent_product_tensor_size_overflow},
            {"independent_product_row_major_flattening_spots", test_independent_product_row_major_flattening_spots},
            {"independent_product_mixture_marginal_is_normalized", test_independent_product_mixture_marginal_is_normalized},
            {"selected_independent_product_one_dimensional_matches_moment_tensor",
             test_selected_independent_product_one_dimensional_matches_moment_tensor},
            {"selected_independent_product_two_dimensional_nonrectangular_matches_moment_tensor",
             test_selected_independent_product_two_dimensional_nonrectangular_matches_moment_tensor},
            {"selected_independent_product_mixture_marginal_matches_moment_tensor",
             test_selected_independent_product_mixture_marginal_matches_moment_tensor},
            {"selected_independent_product_three_dimensional_matches_moment_tensor",
             test_selected_independent_product_three_dimensional_matches_moment_tensor},
            {"selected_independent_product_preserves_reordered_and_duplicate_indices",
             test_selected_independent_product_preserves_reordered_and_duplicate_indices},
            {"selected_independent_product_zero_global_density_skips_marginals",
             test_selected_independent_product_zero_global_density_skips_marginals},
            {"selected_independent_product_empty_request_is_valid", test_selected_independent_product_empty_request_is_valid},
            {"selected_independent_product_invalid_rank_reports_order", test_selected_independent_product_invalid_rank_reports_order},
            {"selected_independent_product_oversized_order_reports_dimension", test_selected_independent_product_oversized_order_reports_dimension},
            {"selected_independent_product_marginal_failure_reports_dimension", test_selected_independent_product_marginal_failure_reports_dimension},
            {"selected_independent_product_marginal_overflow_reports_full_order",
             test_selected_independent_product_marginal_overflow_reports_full_order},
            {"selected_independent_product_product_overflow_reports_full_order",
             test_selected_independent_product_product_overflow_reports_full_order},
            {"selected_independent_product_product_underflow_reports_full_order",
             test_selected_independent_product_product_underflow_reports_full_order},
            {"selected_independent_product_matches_cqmom_layout_subselection", test_selected_independent_product_matches_cqmom_layout_subselection},
        });
}
