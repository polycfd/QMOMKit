#include "qmomkit/qmomkit.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace
{

void assert_near(const double actual, const double expected, const double tolerance)
{
    assert(std::abs(actual - expected) <= tolerance);
}

void assert_near_relative(const double actual, const double expected, const double tolerance)
{
    const double scale = std::max(1.0, std::abs(expected));
    assert(std::abs(actual - expected) <= tolerance * scale);
}

void assert_vector_near(const std::vector<double>& actual, const std::vector<double>& expected, const double tolerance)
{
    assert(actual.size() == expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        assert_near(actual[i], expected[i], tolerance);
    }
}

void assert_vector_near_relative(const std::vector<double>& actual, const std::vector<double>& expected, const double tolerance)
{
    assert(actual.size() == expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        assert_near_relative(actual[i], expected[i], tolerance);
    }
}

void assert_sorted_ascending(const std::vector<double>& values)
{
    for (std::size_t i = 1; i < values.size(); ++i)
    {
        assert(values[i - 1] <= values[i]);
    }
}

std::vector<double> reconstruct_moments(const qmomkit::Quadrature1D& quadrature, const std::size_t n_moments)
{
    std::vector<double> moments(n_moments, 0.0);
    for (std::size_t moment_order = 0; moment_order < n_moments; ++moment_order)
    {
        for (std::size_t node = 0; node < quadrature.weights.size(); ++node)
        {
            moments[moment_order] += quadrature.weights[node] * std::pow(quadrature.abscissas[node], moment_order);
        }
    }
    return moments;
}

void test_adaptive_near_duplicate_nodes_collapse()
{
    const qmomkit::MomentSet moments{1.0, 1.0000000000005, 1.000000000001, 1.0000000000015001};
    qmomkit::WheelerOptions options;
    options.adaptive = true;

    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments, 2, options);

    assert(quadrature.ok());
    assert_vector_near(quadrature.abscissas, {1.0000000000005}, 1.0e-14);
    assert_vector_near(quadrature.weights, {1.0}, 1.0e-14);
}

void test_adaptive_tiny_weight_node_collapse()
{
    const qmomkit::MomentSet moments{1.000000000001, 1.999999999999, 4.000000000001, 7.999999999999};
    qmomkit::WheelerOptions options;
    options.adaptive = true;

    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments, 2, options);

    assert(quadrature.ok());
    assert_vector_near(quadrature.abscissas, {1.9999999999969997}, 1.0e-12);
    assert_vector_near(quadrature.weights, {1.000000000001}, 1.0e-12);
}

void test_scaled_weights_transform_quadrature()
{
    const qmomkit::MomentSet moments{37.5, 41.374999999999993, 199.31249999999997, 411.24874999999992, 1400.1131249999999, 3497.4708874999988};

    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments, 3);

    assert(quadrature.ok());
    assert_sorted_ascending(quadrature.abscissas);
    assert_vector_near(quadrature.abscissas, {-1.7000000000000013, 0.19999999999999335, 2.799999999999998}, 1.0e-12);
    assert_vector_near(quadrature.weights, {11.249999999999989, 5.0000000000000142, 21.250000000000011}, 1.0e-11);
    assert_vector_near_relative(reconstruct_moments(quadrature, moments.size()), moments.values(), 1.0e-13);
}

void test_scaled_abscissas_transform_quadrature()
{
    const qmomkit::MomentSet moments{3.0, 8.2750000000000004, 99.65625, 514.06093750000002, 4375.353515625, 27323.991308593748};

    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments, 3);

    assert(quadrature.ok());
    assert_sorted_ascending(quadrature.abscissas);
    assert_vector_near(quadrature.abscissas, {-4.2500000000000036, 0.49999999999999301, 7.0000000000000009}, 1.0e-11);
    assert_vector_near(quadrature.weights, {0.89999999999999947, 0.4000000000000008, 1.6999999999999997}, 1.0e-11);
    assert_vector_near_relative(reconstruct_moments(quadrature, moments.size()), moments.values(), 1.0e-13);
}

void test_adaptive_exact_degenerate_sequence_reduces_before_finite_check()
{
    const qmomkit::MomentSet moments{1.0, 1.0, 1.0, 1.0};
    qmomkit::WheelerOptions options;
    options.adaptive = true;

    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments, 2, options);

    assert(quadrature.ok());
    assert_vector_near(quadrature.abscissas, {1.0}, 1.0e-14);
    assert_vector_near(quadrature.weights, {1.0}, 1.0e-14);
}

void test_adaptive_exact_large_abscissa_degenerate_sequence_reduces()
{
    const qmomkit::MomentSet moments{1.0, 7.0, 49.0, 343.0};
    qmomkit::WheelerOptions options;
    options.adaptive = true;

    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments, 2, options);

    assert(quadrature.ok());
    assert_vector_near(quadrature.abscissas, {7.0}, 1.0e-14);
    assert_vector_near(quadrature.weights, {1.0}, 1.0e-14);
}

void test_adaptive_roundoff_positive_degenerate_sequence_reduces()
{
    const qmomkit::MomentSet moments{1.0000000000000007, -5.4825828376551026e-16, 2.193033135062041e-15, 0.0};
    qmomkit::WheelerOptions options;
    options.adaptive = true;

    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments, 2, options);

    assert(quadrature.ok());
    assert_vector_near(quadrature.abscissas, {-5.4825828376550987e-16}, 1.0e-14);
    assert_vector_near(quadrature.weights, {1.0000000000000007}, 1.0e-14);
}

void test_adaptive_threshold_rejection_falls_back_to_one_node()
{
    const qmomkit::MomentSet moments{1.0, 1.25, 3.25, 5.75};
    qmomkit::WheelerOptions options;
    options.adaptive = true;
    options.rmin = 1.0;

    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments, 2, options);

    assert(quadrature.ok());
    assert(quadrature.status == qmomkit::QuadratureStatus::success);
    assert_vector_near(quadrature.abscissas, {1.25}, 1.0e-14);
    assert_vector_near(quadrature.weights, {1.0}, 1.0e-14);
}

void test_adaptive_requested_two_node_negative_variance_recovers_one_node()
{
    const qmomkit::MomentSet moments{1.0, 0.0, -1.0, 0.0};
    qmomkit::WheelerOptions options;
    options.adaptive = true;

    // For a two-node request, the failing sigma diagonal is the variance
    // numerator m0*m2 - m1*m1. The adaptive retry uses the lower-order prefix
    // and is not expected to reproduce the discarded inconsistent m2.
    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments, 2, options);

    assert(quadrature.ok());
    assert(quadrature.status == qmomkit::QuadratureStatus::success);
    assert(quadrature.abscissas.size() == 1);
    assert(quadrature.weights.size() == 1);
    assert(std::isfinite(quadrature.abscissas.front()));
    assert(std::isfinite(quadrature.weights.front()));
    assert_vector_near(reconstruct_moments(quadrature, 2), {moments[0], moments[1]}, 1.0e-14);
}

void test_adaptive_requested_three_node_high_rank_sigma_recovers_two_nodes()
{
    const qmomkit::MomentSet moments{1.0, 0.8, 2.8, 4.4, 0.0, 0.0};
    qmomkit::WheelerOptions options;
    options.adaptive = true;

    // The first four moments support a two-node quadrature, but the higher-rank
    // sigma condition needed by the requested three-node inversion fails.
    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments, 3, options);

    assert(quadrature.ok());
    assert(quadrature.status == qmomkit::QuadratureStatus::success);
    assert(quadrature.abscissas.size() == 2);
    assert(quadrature.weights.size() == 2);
    assert_vector_near_relative(reconstruct_moments(quadrature, 4), {1.0, 0.8, 2.8, 4.4}, 1.0e-13);
}

void test_nonadaptive_singular_sequence_reports_failure()
{
    const qmomkit::MomentSet moments{1.0, 1.0, 1.0, 1.0};

    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments, 2);

    assert(!quadrature.ok());
    assert(quadrature.status == qmomkit::QuadratureStatus::not_realizable);
    assert(!quadrature.message.empty());
    assert_vector_near(quadrature.abscissas, {1.0}, 1.0e-14);
    assert_vector_near(quadrature.weights, {1.0}, 1.0e-14);
}

void test_fixed_count_negative_variance_remains_strict_failure()
{
    const qmomkit::MomentSet moments{1.0, 0.0, -1.0, 0.0};

    // Without adaptive mode, a negative two-node variance numerator remains a
    // not-realizable failure even though diagnostic one-node data is available.
    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments, 2);

    assert(!quadrature.ok());
    assert(quadrature.status == qmomkit::QuadratureStatus::not_realizable);
    assert(!quadrature.message.empty());
    assert_vector_near(quadrature.abscissas, {0.0}, 1.0e-14);
    assert_vector_near(quadrature.weights, {1.0}, 1.0e-14);
}

void test_adaptive_scaled_requested_two_node_negative_variance_recovers_one_node()
{
    const qmomkit::MomentSet moments{1.0, 100.0, 9999.0, 999900.0};
    qmomkit::WheelerOptions options;
    options.adaptive = true;

    // The variance numerator is materially negative at large abscissa scale.
    // The one-node adaptive retry still reproduces only m0 and m1.
    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments, 2, options);

    assert(quadrature.ok());
    assert(quadrature.status == qmomkit::QuadratureStatus::success);
    assert_vector_near(quadrature.abscissas, {100.0}, 1.0e-14);
    assert_vector_near(quadrature.weights, {1.0}, 1.0e-14);
    assert_vector_near(reconstruct_moments(quadrature, 2), {moments[0], moments[1]}, 1.0e-14);
}

void test_adaptive_negative_mass_remains_failure()
{
    const qmomkit::MomentSet moments{-1.0, 0.0, 1.0, 0.0};
    qmomkit::WheelerOptions options;
    options.adaptive = true;

    // Adaptive rank reduction must not reinterpret a negative mass moment as a
    // valid lower-rank quadrature.
    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments, 2, options);

    assert(!quadrature.ok());
    assert(quadrature.status == qmomkit::QuadratureStatus::not_realizable);
    assert(!quadrature.message.empty());
}

void test_invalid_input_requires_two_moments_per_node()
{
    const qmomkit::MomentSet moments{1.0, 0.0, 1.0};

    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments, 2);

    assert(!quadrature.ok());
    assert(quadrature.status == qmomkit::QuadratureStatus::invalid_input);
    assert(!quadrature.message.empty());
    assert(quadrature.abscissas.empty());
    assert(quadrature.weights.empty());
}

void test_zero_requested_nodes_reports_invalid_input()
{
    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(qmomkit::MomentSet{1.0, 0.0}, 0);

    assert(!quadrature.ok());
    assert(quadrature.status == qmomkit::QuadratureStatus::invalid_input);
    assert(!quadrature.message.empty());
    assert(quadrature.abscissas.empty());
    assert(quadrature.weights.empty());
}

void test_nonfinite_moment_reports_invalid_input()
{
    const qmomkit::MomentSet moments{1.0, std::numeric_limits<double>::quiet_NaN()};
    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments, 1);

    assert(!quadrature.ok());
    assert(quadrature.status == qmomkit::QuadratureStatus::invalid_input);
    assert(!quadrature.message.empty());
    assert(quadrature.abscissas.empty());
    assert(quadrature.weights.empty());
}

void test_invalid_wheeler_options_fail_before_inversion()
{
    qmomkit::WheelerOptions options;
    options.cutoff = 0.0;
    qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(qmomkit::MomentSet{1.0, 0.0}, 1, options);

    assert(!quadrature.ok());
    assert(quadrature.status == qmomkit::QuadratureStatus::invalid_input);
    assert(!quadrature.message.empty());
    assert(quadrature.abscissas.empty());
    assert(quadrature.weights.empty());

    options.cutoff = 1.0e-30;
    options.rmin = -1.0;
    quadrature = qmomkit::compute_quadrature(qmomkit::MomentSet{1.0, 0.0}, 1, options);

    assert(!quadrature.ok());
    assert(quadrature.status == qmomkit::QuadratureStatus::invalid_input);
    assert(!quadrature.message.empty());

    options.rmin = 1.0e-10;
    options.eabs = std::numeric_limits<double>::quiet_NaN();
    quadrature = qmomkit::compute_quadrature(qmomkit::MomentSet{1.0, 0.0}, 1, options);

    assert(!quadrature.ok());
    assert(quadrature.status == qmomkit::QuadratureStatus::invalid_input);
    assert(!quadrature.message.empty());
}

} // namespace

int main()
{
    test_adaptive_near_duplicate_nodes_collapse();
    test_adaptive_tiny_weight_node_collapse();
    test_scaled_weights_transform_quadrature();
    test_scaled_abscissas_transform_quadrature();
    test_adaptive_exact_degenerate_sequence_reduces_before_finite_check();
    test_adaptive_exact_large_abscissa_degenerate_sequence_reduces();
    test_adaptive_roundoff_positive_degenerate_sequence_reduces();
    test_adaptive_threshold_rejection_falls_back_to_one_node();
    test_adaptive_requested_two_node_negative_variance_recovers_one_node();
    test_adaptive_requested_three_node_high_rank_sigma_recovers_two_nodes();
    test_nonadaptive_singular_sequence_reports_failure();
    test_fixed_count_negative_variance_remains_strict_failure();
    test_adaptive_scaled_requested_two_node_negative_variance_recovers_one_node();
    test_adaptive_negative_mass_remains_failure();
    test_invalid_input_requires_two_moments_per_node();
    test_zero_requested_nodes_reports_invalid_input();
    test_nonfinite_moment_reports_invalid_input();
    test_invalid_wheeler_options_fail_before_inversion();

    return 0;
}
