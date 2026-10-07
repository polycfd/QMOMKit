#include "internal/bounded_quadrature.hpp"
#include "qmomkit/qmomkit.hpp"
#include "test_case_dispatch.hpp"

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

qmomkit::MomentSet moments_from_nodes(const std::vector<double>& abscissas, const std::vector<double>& weights)
{
    return qmomkit::MomentSet(reconstruct_moments(qmomkit::Quadrature1D{abscissas, weights}, 2 * weights.size()));
}

void test_adaptive_two_node_span_separation_uses_larger_abscissa_magnitude()
{
    // Per docs/developer/qmom_wheeler.md, the two-node separation ratio is the
    // node distance divided by the larger abscissa magnitude: unlike the n>2
    // pairwise min/max-distance ratio, this is not translation invariant. The
    // same absolute distance (0.01) is accepted near the origin, where it is
    // large relative to the abscissa magnitudes, but rejected once both
    // abscissas are shifted far from zero, where it becomes comparatively tiny.
    qmomkit::WheelerOptions options;
    options.adaptive = true;
    options.eabs = 0.005;

    const qmomkit::Quadrature1D near_origin = qmomkit::compute_quadrature(moments_from_nodes({0.0, 0.01}, {0.5, 0.5}), 2, options);
    const qmomkit::Quadrature1D translated = qmomkit::compute_quadrature(moments_from_nodes({10.0, 10.01}, {0.5, 0.5}), 2, options);

    assert(near_origin.ok());
    assert(translated.ok());
    assert(near_origin.weights.size() == 2);
    assert(translated.weights.size() == 1);
    assert(std::isfinite(translated.abscissas.front()));
    assert(std::isfinite(translated.weights.front()));
    assert_near(translated.abscissas.front(), 10.005, 1.0e-10);
    assert_near(translated.weights.front(), 1.0, 1.0e-10);
}

void test_adaptive_two_node_span_separation_is_scale_invariant()
{
    qmomkit::WheelerOptions options;
    options.adaptive = true;
    options.eabs = 0.5;

    const qmomkit::Quadrature1D reference = qmomkit::compute_quadrature(moments_from_nodes({0.0, 0.01}, {0.5, 0.5}), 2, options);
    const qmomkit::Quadrature1D scaled = qmomkit::compute_quadrature(moments_from_nodes({0.0, 1.0}, {0.5, 0.5}), 2, options);

    assert(reference.ok());
    assert(scaled.ok());
    assert(reference.weights.size() == 2);
    assert(scaled.weights.size() == 2);
}

void test_adaptive_two_node_relative_separation_accepts_small_symmetric_pair()
{
    qmomkit::WheelerOptions options;
    options.adaptive = true;
    options.eabs = 0.1;

    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments_from_nodes({-1.0e-6, 1.0e-6}, {0.5, 0.5}), 2, options);

    assert(quadrature.ok());
    assert(quadrature.weights.size() == 2);
    assert_vector_near(quadrature.abscissas, {-1.0e-6, 1.0e-6}, 1.0e-14);
}

void test_adaptive_two_node_zero_span_reduces_without_nonfinite_values()
{
    qmomkit::WheelerOptions options;
    options.adaptive = true;

    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments_from_nodes({0.0, 0.0}, {0.5, 0.5}), 2, options);

    assert(quadrature.ok());
    assert(quadrature.weights.size() == 1);
    assert(quadrature.abscissas.size() == 1);
    assert(std::isfinite(quadrature.weights.front()));
    assert(std::isfinite(quadrature.abscissas.front()));
}

void test_adaptive_density_above_cutoff_is_not_reduced_by_rmin()
{
    qmomkit::WheelerOptions options;
    options.adaptive = true;

    // These scales straddle the default rmin value while remaining above the
    // absolute cutoff. The wide support keeps the unchanged sigma criterion
    // well away from its finite-precision reduction boundary.
    for (const double density_scale : {5.0e-11, 2.0e-10})
    {
        const qmomkit::Quadrature1D quadrature =
            qmomkit::compute_quadrature(moments_from_nodes({100.0, 200.0}, {0.5 * density_scale, 0.5 * density_scale}), 2, options);

        assert(quadrature.ok());
        assert(quadrature.weights.size() == 2);
        assert_vector_near_relative(reconstruct_moments(quadrature, 4),
                                    moments_from_nodes({100.0, 200.0}, {0.5 * density_scale, 0.5 * density_scale}).values(), 1.0e-13);
    }
}

void test_cutoff_remains_the_absolute_population_gate()
{
    qmomkit::WheelerOptions options;
    options.adaptive = true;

    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments_from_nodes({100.0, 200.0}, {0.25e-30, 0.25e-30}), 2, options);

    assert(!quadrature.ok());
    assert(quadrature.status == qmomkit::QuadratureStatus::not_realizable);
    assert(quadrature.message == "Moment sequence has m0 < cutoff.");
}

void test_adaptive_three_node_relative_separation_uses_quadrature_span()
{
    qmomkit::WheelerOptions options;
    options.adaptive = true;
    // The global pairwise span gives 0.1 / 1.1 < eabs, while normalization
    // by the largest abscissa magnitude or largest adjacent spacing gives 0.1.
    options.eabs = 0.095;

    const qmomkit::Quadrature1D quadrature =
        qmomkit::compute_quadrature(moments_from_nodes({-1.0, -0.9, 0.1}, {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}), 3, options);

    assert(quadrature.ok());
    assert(quadrature.weights.size() == 2);
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
    // Conditional moments of a single point at zero, on a coordinate of unit
    // scale, with the rounding of the solve that produced them: the variance
    // 2.2e-15 is rounding relative to that scale. As data, the same moments
    // describe a genuine spread of about 5e-8 and keep two nodes.
    const qmomkit::MomentSet moments{1.0000000000000007, -5.4825828376551026e-16, 2.193033135062041e-15, 0.0};
    qmomkit::WheelerOptions options;
    options.adaptive = true;

    assert(qmomkit::compute_quadrature(moments, 2, options).weights.size() == 2);
    const qmomkit::Quadrature1D quadrature = qmomkit::internal::compute_conditional_quadrature(moments, 2, options, 1.0);

    assert(quadrature.ok());
    assert_vector_near(quadrature.abscissas, {-5.4825828376550987e-16}, 1.0e-14);
    assert_vector_near(quadrature.weights, {1.0000000000000007}, 1.0e-14);
}

void test_adaptive_node_count_does_not_depend_on_units_or_number_density()
{
    // sigma(k, k) and the roundoff bound it is compared with both scale as m0
    // times the coordinate to the power 2(k - 1), so neither the unit of the
    // coordinate nor the number density may change the adaptive node count.
    // The conditioning of the cluster limits the nodes to about 1e-5, in fixed
    // mode as well, far below their separation of 0.03.
    const std::vector<double> reference_points{0.95, 0.98, 1.01, 1.04};
    const std::vector<double> reference_weights{0.2, 0.3, 0.3, 0.2};
    qmomkit::WheelerOptions options;
    options.adaptive = true;
    for (const double unit : {1.0e-9, 1.0e-6, 1.0, 1.0e4, 1.0e9})
    {
        for (const double number_density : {1.0e-14, 1.0, 1.0e10})
        {
            std::vector<double> points = reference_points;
            std::vector<double> weights = reference_weights;
            for (std::size_t node = 0; node < points.size(); ++node)
            {
                points[node] *= unit;
                weights[node] *= number_density;
            }
            const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments_from_nodes(points, weights), 4, options);
            assert(quadrature.ok());
            assert(quadrature.weights.size() == 4);
            for (std::size_t node = 0; node < points.size(); ++node)
            {
                assert_near(quadrature.abscissas[node] / unit, reference_points[node], 1.0e-4);
            }
        }
    }
}

void test_adaptive_degenerate_sequences_reduce_at_any_unit_and_location()
{
    // Two points requested as four nodes, and one point requested as three,
    // must reduce to the points they are, in any unit and wherever they sit:
    // the rounding of their moments stays within the bound, which grows with
    // the cancellation far from the origin.
    qmomkit::WheelerOptions options;
    options.adaptive = true;
    for (const double unit : {1.0e-6, 1.0, 1.0e4})
    {
        for (const double location : {0.0, 10.0, 100.0})
        {
            const std::vector<double> two_points{unit * (location - 0.5), unit * (location + 0.5)};
            const qmomkit::Quadrature1D two =
                qmomkit::compute_quadrature(qmomkit::MomentSet(reconstruct_moments(qmomkit::Quadrature1D{two_points, {0.4, 0.6}}, 8)), 4, options);
            assert(two.ok());
            assert(two.weights.size() == 2);

            const std::vector<double> one_point{unit * (location + 0.3)};
            const qmomkit::Quadrature1D one =
                qmomkit::compute_quadrature(qmomkit::MomentSet(reconstruct_moments(qmomkit::Quadrature1D{one_point, {1.0}}, 6)), 3, options);
            assert(one.ok());
            assert(one.weights.size() == 1);
        }
    }
}

void test_adaptive_relative_weight_threshold_rejection_falls_back_to_one_node()
{
    qmomkit::WheelerOptions options;
    options.adaptive = true;
    options.rmin = 0.25;

    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments_from_nodes({-1.0, 2.0}, {0.2, 0.8}), 2, options);

    assert(quadrature.ok());
    assert(quadrature.status == qmomkit::QuadratureStatus::success);
    assert_vector_near(quadrature.abscissas, {1.4}, 1.0e-14);
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
    qmomkit::WheelerOptions options;
    options.adaptive = false;

    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments, 2, options);

    assert(!quadrature.ok());
    assert(quadrature.status == qmomkit::QuadratureStatus::not_realizable);
    assert(!quadrature.message.empty());
    assert_vector_near(quadrature.abscissas, {1.0}, 1.0e-14);
    assert_vector_near(quadrature.weights, {1.0}, 1.0e-14);
}

void test_fixed_count_negative_variance_remains_strict_failure()
{
    const qmomkit::MomentSet moments{1.0, 0.0, -1.0, 0.0};

    qmomkit::WheelerOptions options;
    options.adaptive = false;

    // Without adaptive mode, a negative two-node variance numerator remains a
    // not-realizable failure even though diagnostic one-node data is available.
    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments, 2, options);

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

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"adaptive_two_node_span_separation_uses_larger_abscissa_magnitude",
             test_adaptive_two_node_span_separation_uses_larger_abscissa_magnitude},
            {"adaptive_two_node_span_separation_is_scale_invariant", test_adaptive_two_node_span_separation_is_scale_invariant},
            {"adaptive_two_node_relative_separation_accepts_small_symmetric_pair",
             test_adaptive_two_node_relative_separation_accepts_small_symmetric_pair},
            {"adaptive_two_node_zero_span_reduces_without_nonfinite_values", test_adaptive_two_node_zero_span_reduces_without_nonfinite_values},
            {"adaptive_density_above_cutoff_is_not_reduced_by_rmin", test_adaptive_density_above_cutoff_is_not_reduced_by_rmin},
            {"cutoff_remains_the_absolute_population_gate", test_cutoff_remains_the_absolute_population_gate},
            {"adaptive_three_node_relative_separation_uses_quadrature_span", test_adaptive_three_node_relative_separation_uses_quadrature_span},
            {"adaptive_near_duplicate_nodes_collapse", test_adaptive_near_duplicate_nodes_collapse},
            {"adaptive_tiny_weight_node_collapse", test_adaptive_tiny_weight_node_collapse},
            {"scaled_weights_transform_quadrature", test_scaled_weights_transform_quadrature},
            {"scaled_abscissas_transform_quadrature", test_scaled_abscissas_transform_quadrature},
            {"adaptive_exact_degenerate_sequence_reduces_before_finite_check", test_adaptive_exact_degenerate_sequence_reduces_before_finite_check},
            {"adaptive_exact_large_abscissa_degenerate_sequence_reduces", test_adaptive_exact_large_abscissa_degenerate_sequence_reduces},
            {"adaptive_roundoff_positive_degenerate_sequence_reduces", test_adaptive_roundoff_positive_degenerate_sequence_reduces},
            {"adaptive_node_count_does_not_depend_on_units_or_number_density", test_adaptive_node_count_does_not_depend_on_units_or_number_density},
            {"adaptive_degenerate_sequences_reduce_at_any_unit_and_location", test_adaptive_degenerate_sequences_reduce_at_any_unit_and_location},
            {"adaptive_relative_weight_threshold_rejection_falls_back_to_one_node",
             test_adaptive_relative_weight_threshold_rejection_falls_back_to_one_node},
            {"adaptive_requested_two_node_negative_variance_recovers_one_node", test_adaptive_requested_two_node_negative_variance_recovers_one_node},
            {"adaptive_requested_three_node_high_rank_sigma_recovers_two_nodes",
             test_adaptive_requested_three_node_high_rank_sigma_recovers_two_nodes},
            {"nonadaptive_singular_sequence_reports_failure", test_nonadaptive_singular_sequence_reports_failure},
            {"fixed_count_negative_variance_remains_strict_failure", test_fixed_count_negative_variance_remains_strict_failure},
            {"adaptive_scaled_requested_two_node_negative_variance_recovers_one_node",
             test_adaptive_scaled_requested_two_node_negative_variance_recovers_one_node},
            {"adaptive_negative_mass_remains_failure", test_adaptive_negative_mass_remains_failure},
            {"invalid_input_requires_two_moments_per_node", test_invalid_input_requires_two_moments_per_node},
            {"zero_requested_nodes_reports_invalid_input", test_zero_requested_nodes_reports_invalid_input},
            {"nonfinite_moment_reports_invalid_input", test_nonfinite_moment_reports_invalid_input},
            {"invalid_wheeler_options_fail_before_inversion", test_invalid_wheeler_options_fail_before_inversion},
        });
}
