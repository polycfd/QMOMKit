#include "qmomkit/distributions.hpp"

#include "qmomkit/cqmom.hpp"
#include "test_case_dispatch.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
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

void assert_selected_matches_dense(const qmomkit::DistributionSpecification& distribution, const qmomkit::RawMomentRequest& request,
                                   const double total_weight, const std::size_t maximum_order, const double tolerance)
{
    const qmomkit::RawMomentResult dense = qmomkit::compute_raw_moments(distribution, maximum_order, total_weight);
    const qmomkit::SelectedRawMomentResult selected = qmomkit::compute_raw_moments(distribution, request, total_weight);

    assert(dense.ok());
    assert(selected.ok());
    assert(selected.indices == request.indices);
    assert(selected.values.size() == request.indices.size());
    for (std::size_t index = 0; index < request.indices.size(); ++index)
    {
        assert(request.indices[index].size() == 1);
        assert_near(selected.values[index], dense.moments[request.indices[index][0]], tolerance);
    }
}

void test_selected_distribution_contiguous_orders_match_dense()
{
    const qmomkit::RawMomentRequest request{{{0}, {1}, {2}, {3}}};
    assert_selected_matches_dense(qmomkit::NormalDistribution{0.5, 1.25}, request, 2.0, 3, 1.0e-13);
    assert_selected_matches_dense(qmomkit::LogNormalDistribution{0.1, 0.2}, request, 2.0, 3, 1.0e-13);
    assert_selected_matches_dense(qmomkit::BetaDistribution{2.0, 5.0}, request, 2.0, 3, 1.0e-13);
    assert_selected_matches_dense(qmomkit::GammaDistribution{1.5, 0.75}, request, 2.0, 3, 1.0e-13);
    assert_selected_matches_dense(qmomkit::ChiDistribution{4.0}, request, 2.0, 3, 1.0e-13);
    assert_selected_matches_dense(qmomkit::ChiSquaredDistribution{4.0}, request, 2.0, 3, 1.0e-13);
    assert_selected_matches_dense(qmomkit::WeibullDistribution{1.5, 0.75}, request, 2.0, 3, 1.0e-13);
}

void test_selected_distribution_noncontiguous_orders_match_dense()
{
    const qmomkit::RawMomentRequest request{{{0}, {2}, {5}}};
    assert_selected_matches_dense(qmomkit::NormalDistribution{0.25, 0.75}, request, 1.5, 5, 1.0e-12);
    assert_selected_matches_dense(qmomkit::LogNormalDistribution{-0.1, 0.15}, request, 1.5, 5, 1.0e-12);
    assert_selected_matches_dense(qmomkit::BetaDistribution{3.0, 4.0}, request, 1.5, 5, 1.0e-12);
    assert_selected_matches_dense(qmomkit::GammaDistribution{2.0, 0.5}, request, 1.5, 5, 1.0e-12);
    assert_selected_matches_dense(qmomkit::ChiDistribution{5.0}, request, 1.5, 5, 1.0e-12);
    assert_selected_matches_dense(qmomkit::ChiSquaredDistribution{5.0}, request, 1.5, 5, 1.0e-12);
    assert_selected_matches_dense(qmomkit::WeibullDistribution{2.0, 0.5}, request, 1.5, 5, 1.0e-12);
}

void test_selected_distribution_reordered_and_duplicate_orders_match_dense()
{
    const qmomkit::RawMomentRequest request{{{5}, {0}, {2}, {2}, {1}}};
    assert_selected_matches_dense(qmomkit::NormalDistribution{0.25, 0.75}, request, 1.25, 5, 1.0e-12);
    assert_selected_matches_dense(qmomkit::LogNormalDistribution{-0.1, 0.15}, request, 1.25, 5, 1.0e-12);
    assert_selected_matches_dense(qmomkit::BetaDistribution{3.0, 4.0}, request, 1.25, 5, 1.0e-12);
    assert_selected_matches_dense(qmomkit::GammaDistribution{2.0, 0.5}, request, 1.25, 5, 1.0e-12);
    assert_selected_matches_dense(qmomkit::ChiDistribution{5.0}, request, 1.25, 5, 1.0e-12);
    assert_selected_matches_dense(qmomkit::ChiSquaredDistribution{5.0}, request, 1.25, 5, 1.0e-12);
    assert_selected_matches_dense(qmomkit::WeibullDistribution{2.0, 0.5}, request, 1.25, 5, 1.0e-12);
}

void test_selected_distribution_high_order_gaps_match_dense()
{
    const qmomkit::RawMomentRequest request{{{0}, {1}, {12}}};
    assert_selected_matches_dense(qmomkit::NormalDistribution{0.1, 0.5}, request, 1.0, 12, 1.0e-9);
    assert_selected_matches_dense(qmomkit::LogNormalDistribution{-0.05, 0.1}, request, 1.0, 12, 1.0e-9);
    assert_selected_matches_dense(qmomkit::BetaDistribution{2.0, 3.0}, request, 1.0, 12, 1.0e-12);
    assert_selected_matches_dense(qmomkit::GammaDistribution{1.5, 0.25}, request, 1.0, 12, 1.0e-9);
    assert_selected_matches_dense(qmomkit::ChiDistribution{3.0}, request, 1.0, 12, 1.0e-8);
    assert_selected_matches_dense(qmomkit::ChiSquaredDistribution{3.0}, request, 1.0, 12, 1.0e-6);
    assert_selected_matches_dense(qmomkit::WeibullDistribution{2.0, 0.25}, request, 1.0, 12, 1.0e-9);
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

void test_overflow_status()
{
    const qmomkit::RawMomentResult log_normal = qmomkit::compute_raw_moments(qmomkit::LogNormalDistribution{100.0, 10.0}, 20);

    assert(!log_normal.ok());
    assert(log_normal.status == qmomkit::RawMomentStatus::overflow);
    assert(log_normal.failed_order > 0);
    assert(!log_normal.message.empty());

    const double large_scale = std::numeric_limits<double>::max() / 2.0;
    assert(qmomkit::compute_raw_moments(qmomkit::GammaDistribution{1.0, large_scale}, 2).status == qmomkit::RawMomentStatus::overflow);
    assert(qmomkit::compute_raw_moments(qmomkit::ChiDistribution{3.0}, 2048).status == qmomkit::RawMomentStatus::overflow);
    assert(qmomkit::compute_raw_moments(qmomkit::ChiSquaredDistribution{2.0}, 1024).status == qmomkit::RawMomentStatus::overflow);
    assert(qmomkit::compute_raw_moments(qmomkit::WeibullDistribution{0.1, 2.0}, 100).status == qmomkit::RawMomentStatus::overflow);
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

void assert_selected_mixture_matches_dense(const qmomkit::DistributionMixture& mixture, const qmomkit::RawMomentRequest& request,
                                           const std::size_t maximum_order, const double tolerance)
{
    const qmomkit::RawMomentResult dense = qmomkit::compute_raw_moments(mixture, maximum_order);
    const qmomkit::SelectedRawMomentResult selected = qmomkit::compute_raw_moments(mixture, request);

    assert(dense.ok());
    assert(selected.ok());
    assert(selected.indices == request.indices);
    assert(selected.values.size() == request.indices.size());
    for (std::size_t index = 0; index < request.indices.size(); ++index)
    {
        assert(request.indices[index].size() == 1);
        assert_near(selected.values[index], dense.moments[request.indices[index][0]], tolerance);
    }
}

void test_selected_one_component_mixture_matches_dense_subselection()
{
    const qmomkit::DistributionMixture mixture{{{qmomkit::GammaDistribution{2.0, 3.0}, 4.0}}};
    assert_selected_mixture_matches_dense(mixture, {{{0}, {2}, {4}}}, 4, 1.0e-12);
}

void test_selected_heterogeneous_mixture_matches_dense_subselection()
{
    const qmomkit::DistributionMixture mixture{
        {{qmomkit::GammaDistribution{2.0, 3.0}, 2.0}, {qmomkit::BetaDistribution{2.0, 3.0}, 5.0}, {qmomkit::NormalDistribution{1.0, 0.5}, 1.5}}};
    assert_selected_mixture_matches_dense(mixture, {{{0}, {1}, {3}}}, 3, 1.0e-12);
}

void test_selected_mixture_preserves_reordered_and_duplicate_requests()
{
    const qmomkit::DistributionMixture mixture{{{qmomkit::GammaDistribution{1.5, 0.5}, 2.0}, {qmomkit::ChiSquaredDistribution{4.0}, 0.75}}};
    assert_selected_mixture_matches_dense(mixture, {{{4}, {0}, {2}, {2}, {1}}}, 4, 1.0e-11);
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

void assert_selected_product_matches_dense(const qmomkit::IndependentProductDistribution& product, const qmomkit::RawMomentRequest& request,
                                           const std::vector<std::size_t>& maximum_orders, const double tolerance)
{
    const qmomkit::TensorMomentResult dense = qmomkit::compute_raw_moments(product, maximum_orders);
    const qmomkit::SelectedRawMomentResult selected = qmomkit::compute_raw_moments(product, request);

    assert(dense.ok());
    assert(selected.ok());
    assert(selected.indices == request.indices);
    assert(selected.values.size() == request.indices.size());
    for (std::size_t index = 0; index < request.indices.size(); ++index)
    {
        assert_near(selected.values[index], dense.moments.values[qmomkit::flatten_multi_index(request.indices[index], dense.moments.extents)],
                    tolerance);
    }
}

void test_selected_independent_product_one_dimensional_matches_dense()
{
    const qmomkit::IndependentProductDistribution product{{qmomkit::GammaDistribution{1.5, 2.0}}, 2.5};
    assert_selected_product_matches_dense(product, {{{0}, {3}, {1}}}, {3}, 1.0e-12);
}

void test_selected_independent_product_two_dimensional_nonrectangular_matches_dense()
{
    const qmomkit::IndependentProductDistribution product{{qmomkit::NormalDistribution{2.0, 1.0}, qmomkit::GammaDistribution{1.0, 1.0}}, 3.0};
    assert_selected_product_matches_dense(product, {{{0, 0}, {2, 1}, {1, 3}}}, {2, 3}, 1.0e-12);
}

void test_selected_independent_product_mixture_marginal_matches_dense()
{
    const qmomkit::DistributionMixture mixture{
        {{qmomkit::LogNormalDistribution{-0.3, 0.25}, 1.5}, {qmomkit::LogNormalDistribution{0.4, 0.45}, 3.25}}};
    const qmomkit::IndependentProductDistribution product{{qmomkit::DistributionSpecification{qmomkit::NormalDistribution{-0.4, 0.7}}, mixture},
                                                          2.75};
    assert_selected_product_matches_dense(product, {{{1, 0}, {0, 1}, {1, 1}, {0, 2}}}, {1, 2}, 1.0e-12);
}

void test_selected_independent_product_three_dimensional_matches_dense()
{
    const qmomkit::IndependentProductDistribution product{
        {qmomkit::NormalDistribution{1.0, 0.5}, qmomkit::BetaDistribution{2.0, 3.0}, qmomkit::ChiSquaredDistribution{4.0}}, 1.25};
    assert_selected_product_matches_dense(product, {{{0, 0, 0}, {2, 1, 3}, {1, 2, 0}}}, {2, 2, 3}, 1.0e-10);
}

void test_selected_independent_product_preserves_reordered_and_duplicate_indices()
{
    const qmomkit::IndependentProductDistribution product{{qmomkit::NormalDistribution{2.0, 1.0}, qmomkit::GammaDistribution{1.0, 1.0}}, 3.0};
    assert_selected_product_matches_dense(product, {{{1, 2}, {0, 0}, {1, 2}, {0, 1}}}, {1, 2}, 1.0e-12);
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

    const qmomkit::TensorMomentResult dense = qmomkit::compute_raw_moments(product, {3, 3});
    const qmomkit::CqmomMomentData required_moment_data = qmomkit::extract_cqmom_moments(dense.moments, layout);
    const qmomkit::SelectedRawMomentResult selected = qmomkit::compute_raw_moments(product, qmomkit::RawMomentRequest{layout.required_indices});

    assert(dense.ok());
    assert(selected.ok());
    assert(selected.indices == layout.required_indices);
    assert(selected.values.size() == required_moment_data.values.size());
    for (std::size_t index = 0; index < selected.values.size(); ++index)
    {
        assert_near(selected.values[index], required_moment_data.values[index], 1.0e-12);
    }
}

} // namespace

void run_all_distribution_tests()
{
    test_selected_raw_moment_result_preserves_request_order();
    test_selected_raw_moment_result_accepts_noncontiguous_multidimensional_indices();
    test_selected_raw_moment_result_preserves_duplicate_indices();
    test_selected_raw_moment_empty_request_is_valid();
    test_selected_raw_moment_invalid_requests();
    test_selected_distribution_contiguous_orders_match_dense();
    test_selected_distribution_noncontiguous_orders_match_dense();
    test_selected_distribution_reordered_and_duplicate_orders_match_dense();
    test_selected_distribution_high_order_gaps_match_dense();
    test_selected_distribution_zero_total_weight();
    test_selected_distribution_invalid_request_rank();
    test_selected_distribution_overflow_and_underflow_failures();
    test_variant_dispatch_and_nonunit_weight();
    test_maximum_order_zero();
    test_zero_total_weight_returns_zero_sequence();
    test_invalid_total_weight();
    test_invalid_distribution_parameters();
    test_overflow_status();
    test_underflow_status();
    test_empty_mixture_is_invalid();
    test_mixture_maximum_order_zero();
    test_one_component_mixture_matches_distribution_result();
    test_heterogeneous_mixture_dispatch();
    test_zero_density_component_is_skipped();
    test_zero_density_extreme_component_is_skipped();
    test_invalid_mixture_densities();
    test_all_zero_density_mixture_is_invalid();
    test_positive_component_failure_reports_index();
    test_surrounded_positive_component_failure_reports_index();
    test_mixture_accumulation_overflow_reports_index();
    test_selected_one_component_mixture_matches_dense_subselection();
    test_selected_heterogeneous_mixture_matches_dense_subselection();
    test_selected_mixture_preserves_reordered_and_duplicate_requests();
    test_selected_mixture_zero_density_component_is_skipped();
    test_selected_mixture_empty_request_is_valid_for_valid_mixture();
    test_selected_mixture_invalid_component_failure_reports_index();
    test_selected_mixture_overflow_propagates_component_and_order();
    test_selected_mixture_underflow_propagates_component_and_order();
    test_selected_mixture_accumulation_overflow_reports_index();
    test_independent_product_zero_dimension_is_invalid();
    test_independent_product_rank_mismatch_is_invalid();
    test_independent_product_maximum_order_overflow_guard();
    test_independent_product_invalid_global_number_density();
    test_independent_product_invalid_marginal_reports_dimension();
    test_independent_product_marginal_failure_reports_order();
    test_independent_product_maximum_orders_zero();
    test_independent_product_zero_global_density_skips_marginals();
    test_independent_product_zero_global_density_still_validates_extents();
    test_independent_product_tensor_size_overflow();
    test_independent_product_row_major_flattening_spots();
    test_independent_product_mixture_marginal_is_normalized();
    test_selected_independent_product_one_dimensional_matches_dense();
    test_selected_independent_product_two_dimensional_nonrectangular_matches_dense();
    test_selected_independent_product_mixture_marginal_matches_dense();
    test_selected_independent_product_three_dimensional_matches_dense();
    test_selected_independent_product_preserves_reordered_and_duplicate_indices();
    test_selected_independent_product_zero_global_density_skips_marginals();
    test_selected_independent_product_empty_request_is_valid();
    test_selected_independent_product_invalid_rank_reports_order();
    test_selected_independent_product_oversized_order_reports_dimension();
    test_selected_independent_product_marginal_failure_reports_dimension();
    test_selected_independent_product_marginal_overflow_reports_full_order();
    test_selected_independent_product_product_overflow_reports_full_order();
    test_selected_independent_product_product_underflow_reports_full_order();
    test_selected_independent_product_matches_cqmom_layout_subselection();
}

void run_selected_result_tests()
{
    test_selected_raw_moment_result_preserves_request_order();
    test_selected_raw_moment_result_accepts_noncontiguous_multidimensional_indices();
    test_selected_raw_moment_result_preserves_duplicate_indices();
    test_selected_raw_moment_empty_request_is_valid();
    test_selected_raw_moment_invalid_requests();
    test_selected_distribution_contiguous_orders_match_dense();
    test_selected_distribution_noncontiguous_orders_match_dense();
    test_selected_distribution_reordered_and_duplicate_orders_match_dense();
    test_selected_distribution_high_order_gaps_match_dense();
    test_selected_distribution_zero_total_weight();
    test_selected_distribution_invalid_request_rank();
    test_selected_distribution_overflow_and_underflow_failures();
    test_variant_dispatch_and_nonunit_weight();
    test_maximum_order_zero();
    test_zero_total_weight_returns_zero_sequence();
    test_invalid_total_weight();
    test_invalid_distribution_parameters();
    test_overflow_status();
    test_underflow_status();
}

void run_mixture_tests()
{
    test_empty_mixture_is_invalid();
    test_mixture_maximum_order_zero();
    test_one_component_mixture_matches_distribution_result();
    test_heterogeneous_mixture_dispatch();
    test_zero_density_component_is_skipped();
    test_zero_density_extreme_component_is_skipped();
    test_invalid_mixture_densities();
    test_all_zero_density_mixture_is_invalid();
    test_positive_component_failure_reports_index();
    test_surrounded_positive_component_failure_reports_index();
    test_mixture_accumulation_overflow_reports_index();
    test_selected_one_component_mixture_matches_dense_subselection();
    test_selected_heterogeneous_mixture_matches_dense_subselection();
    test_selected_mixture_preserves_reordered_and_duplicate_requests();
    test_selected_mixture_zero_density_component_is_skipped();
    test_selected_mixture_empty_request_is_valid_for_valid_mixture();
    test_selected_mixture_invalid_component_failure_reports_index();
    test_selected_mixture_overflow_propagates_component_and_order();
    test_selected_mixture_underflow_propagates_component_and_order();
    test_selected_mixture_accumulation_overflow_reports_index();
}

void run_independent_product_tests()
{
    test_independent_product_zero_dimension_is_invalid();
    test_independent_product_rank_mismatch_is_invalid();
    test_independent_product_maximum_order_overflow_guard();
    test_independent_product_invalid_global_number_density();
    test_independent_product_invalid_marginal_reports_dimension();
    test_independent_product_marginal_failure_reports_order();
    test_independent_product_maximum_orders_zero();
    test_independent_product_zero_global_density_skips_marginals();
    test_independent_product_zero_global_density_still_validates_extents();
    test_independent_product_tensor_size_overflow();
    test_independent_product_row_major_flattening_spots();
    test_independent_product_mixture_marginal_is_normalized();
    test_selected_independent_product_one_dimensional_matches_dense();
    test_selected_independent_product_two_dimensional_nonrectangular_matches_dense();
    test_selected_independent_product_mixture_marginal_matches_dense();
    test_selected_independent_product_three_dimensional_matches_dense();
    test_selected_independent_product_preserves_reordered_and_duplicate_indices();
    test_selected_independent_product_zero_global_density_skips_marginals();
    test_selected_independent_product_empty_request_is_valid();
    test_selected_independent_product_invalid_rank_reports_order();
    test_selected_independent_product_oversized_order_reports_dimension();
    test_selected_independent_product_marginal_failure_reports_dimension();
    test_selected_independent_product_marginal_overflow_reports_full_order();
    test_selected_independent_product_product_overflow_reports_full_order();
    test_selected_independent_product_product_underflow_reports_full_order();
    test_selected_independent_product_matches_cqmom_layout_subselection();
}

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
            {"selected_distribution_contiguous_orders_match_dense", test_selected_distribution_contiguous_orders_match_dense},
            {"selected_distribution_noncontiguous_orders_match_dense", test_selected_distribution_noncontiguous_orders_match_dense},
            {"selected_distribution_reordered_and_duplicate_orders_match_dense",
             test_selected_distribution_reordered_and_duplicate_orders_match_dense},
            {"selected_distribution_high_order_gaps_match_dense", test_selected_distribution_high_order_gaps_match_dense},
            {"selected_distribution_zero_total_weight", test_selected_distribution_zero_total_weight},
            {"selected_distribution_invalid_request_rank", test_selected_distribution_invalid_request_rank},
            {"selected_distribution_overflow_and_underflow_failures", test_selected_distribution_overflow_and_underflow_failures},
            {"variant_dispatch_and_nonunit_weight", test_variant_dispatch_and_nonunit_weight},
            {"maximum_order_zero", test_maximum_order_zero},
            {"zero_total_weight_returns_zero_sequence", test_zero_total_weight_returns_zero_sequence},
            {"invalid_total_weight", test_invalid_total_weight},
            {"invalid_distribution_parameters", test_invalid_distribution_parameters},
            {"overflow_status", test_overflow_status},
            {"underflow_status", test_underflow_status},
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
            {"selected_one_component_mixture_matches_dense_subselection", test_selected_one_component_mixture_matches_dense_subselection},
            {"selected_heterogeneous_mixture_matches_dense_subselection", test_selected_heterogeneous_mixture_matches_dense_subselection},
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
            {"selected_independent_product_one_dimensional_matches_dense", test_selected_independent_product_one_dimensional_matches_dense},
            {"selected_independent_product_two_dimensional_nonrectangular_matches_dense",
             test_selected_independent_product_two_dimensional_nonrectangular_matches_dense},
            {"selected_independent_product_mixture_marginal_matches_dense", test_selected_independent_product_mixture_marginal_matches_dense},
            {"selected_independent_product_three_dimensional_matches_dense", test_selected_independent_product_three_dimensional_matches_dense},
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
