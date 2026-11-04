#include "internal/bounded_quadrature.hpp"
#include "internal/hyperbolic_quadrature.hpp"
#include "qmomkit/qmomkit.hpp"
#include "test_case_dispatch.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <tuple>
#include <utility>
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
    // Per docs/developer/qmom.md, the two-node separation ratio is the
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
    // times the coordinate to the power 2k, so neither the unit of the
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
    // numerator m0*m2 - m1*m1. The adaptive reduction keeps the one-node rule
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
    // The one-node rule of the adaptive reduction reproduces only m0 and m1.
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

void test_fixed_count_refuses_sequences_supporting_fewer_nodes()
{
    // Moments of a rule with fewer weighted points than requested lie on the
    // boundary, where a recurrence entry is zero up to roundoff. A fixed count
    // must refuse them whatever sign that roundoff takes, which changes with
    // the unit and the number density, and name the count the moments support.
    qmomkit::WheelerOptions fixed;
    fixed.adaptive = false;
    for (const double unit : {1.0e-6, 0.37, 1.0, 3.0, 1.0e6})
    {
        for (const double number_density : {1.0e-12, 1.0, 1.0e12})
        {
            const qmomkit::Quadrature1D dirac{{0.4 * unit}, {number_density}};
            const qmomkit::Quadrature1D two_points{{-0.75 * unit, 1.25 * unit}, {0.4 * number_density, 0.6 * number_density}};
            for (const auto& [rule, requested, message] :
                 {std::tuple{dirac, std::size_t{2}, "supports at most 1 node:"}, std::tuple{two_points, std::size_t{3}, "supports at most 2 nodes:"},
                  std::tuple{two_points, std::size_t{4}, "supports at most 2 nodes:"}})
            {
                const qmomkit::Quadrature1D result =
                    qmomkit::compute_quadrature(qmomkit::MomentSet(reconstruct_moments(rule, 2 * requested)), requested, fixed);
                assert(result.status == qmomkit::QuadratureStatus::not_realizable);
                assert(result.message.find(message) != std::string::npos);
            }

            // A rule with as many points as requested is still recovered.
            const std::vector<double> abscissas{-0.75 * unit, 0.25 * unit, 1.25 * unit};
            const std::vector<double> weights{0.3 * number_density, 0.3 * number_density, 0.4 * number_density};
            const qmomkit::Quadrature1D recovered = qmomkit::compute_quadrature(moments_from_nodes(abscissas, weights), 3, fixed);
            assert(recovered.ok());
            assert_vector_near_relative(recovered.abscissas, abscissas, 1.0e-12);
            assert_vector_near_relative(recovered.weights, weights, 1.0e-12);
        }
    }
}

// The rounding bound tau_k of docs/developer/qmom.md at coordinate scale zero,
// from the full sigma table of the moments and its absolute-value counterpart.
// Row j of both is stored at index j + 1, so that row -1 is the zero row 0.
double documented_rounding_bound(const std::vector<double>& moments, const std::size_t k)
{
    const std::size_t count = moments.size();
    std::vector<std::vector<double>> sigma(k + 2, std::vector<double>(count, 0.0));
    std::vector<std::vector<double>> bound = sigma;
    std::vector<double> a(k + 1, 0.0);
    std::vector<double> b(k + 1, 0.0);
    for (std::size_t l = 0; l < count; ++l)
    {
        sigma[1][l] = moments[l];
        bound[1][l] = std::abs(moments[l]);
    }
    a[0] = moments[1] / moments[0];
    for (std::size_t j = 1; j <= k; ++j)
    {
        for (std::size_t l = j; l + j < count; ++l)
        {
            sigma[j + 1][l] = sigma[j][l + 1] - a[j - 1] * sigma[j][l] - b[j - 1] * sigma[j - 1][l];
            bound[j + 1][l] = bound[j][l + 1] + std::abs(a[j - 1]) * bound[j][l] + std::abs(b[j - 1]) * bound[j - 1][l];
        }
        a[j] = sigma[j + 1][j + 1] / sigma[j + 1][j] - sigma[j][j] / sigma[j][j - 1];
        b[j] = sigma[j + 1][j] / sigma[j][j - 1];
    }
    return 100.0 * std::numeric_limits<double>::epsilon() * bound[k + 1][k];
}

void test_fixed_count_reports_the_documented_rounding_bound()
{
    // Three weighted points requested as four nodes are refused at sigma(3, 3),
    // whose bound tau_3 takes every term of the bound recurrence, b_2 included,
    // from rows 1 and 2. The refusal prints tau_3, which must be the documented
    // bound, near the origin and at 10, where it carries the cancellation.
    qmomkit::WheelerOptions fixed;
    fixed.adaptive = false;
    for (const double location : {0.0, 10.0})
    {
        const qmomkit::Quadrature1D points{{location - 0.75, location + 0.25, location + 1.25}, {0.3, 0.3, 0.4}};
        const std::vector<double> moments = reconstruct_moments(points, 8);
        const qmomkit::Quadrature1D result = qmomkit::compute_quadrature(qmomkit::MomentSet(moments), 4, fixed);
        assert(result.status == qmomkit::QuadratureStatus::not_realizable);
        assert(result.message.find("supports at most 3 nodes: sigma(3, 3) = ") != std::string::npos);

        const std::string marker = "within the tolerance ";
        const std::size_t position = result.message.find(marker);
        assert(position != std::string::npos);
        const double printed = std::stod(result.message.substr(position + marker.size()));
        const double expected = documented_rounding_bound(moments, 3);
        assert(std::abs(printed - expected) <= 1.0e-12 * expected);
    }
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

void assert_failed_recurrence(const qmomkit::WheelerRecurrence& recurrence, const qmomkit::QuadratureStatus status, const std::string& message)
{
    assert(!recurrence.ok());
    assert(recurrence.status == status);
    assert(recurrence.message == message);
    assert(recurrence.a.empty());
    assert(recurrence.b.empty());
    assert(recurrence.sigma.empty());
    assert(recurrence.tolerance.empty());
}

void test_recurrence_reports_input_failures()
{
    using qmomkit::QuadratureStatus;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const qmomkit::MomentSet two_points = moments_from_nodes({-1.0, 2.0}, {0.4, 0.6});

    assert_failed_recurrence(qmomkit::compute_wheeler_recurrence(two_points, 0), QuadratureStatus::invalid_input, "node_count must be positive.");
    assert_failed_recurrence(qmomkit::compute_wheeler_recurrence(two_points, std::numeric_limits<std::size_t>::max()),
                             QuadratureStatus::invalid_input, "node_count is too large to derive the required moment count.");
    assert_failed_recurrence(qmomkit::compute_wheeler_recurrence(qmomkit::MomentSet{1.0, 0.8, 2.8}, 2), QuadratureStatus::invalid_input,
                             "Moment sequence has fewer than 2*n entries.");
    assert_failed_recurrence(qmomkit::compute_wheeler_recurrence(qmomkit::MomentSet{1.0, 0.8, nan, 4.4}, 2), QuadratureStatus::invalid_input,
                             "Moment sequence contains non-finite values.");
    assert_failed_recurrence(qmomkit::compute_wheeler_recurrence(qmomkit::MomentSet{-1.0, 0.8, 2.8, 4.4}, 2), QuadratureStatus::not_realizable,
                             "Moment sequence is not realizable: m0 < 0.");
    for (const double zero : {0.0, -0.0})
    {
        assert_failed_recurrence(qmomkit::compute_wheeler_recurrence(qmomkit::MomentSet{zero, 0.8, 2.8, 4.4}, 2), QuadratureStatus::not_realizable,
                                 "Moment sequence has m0 = 0.");
    }

    // Moments beyond m_{2n-1} are not read.
    assert(qmomkit::compute_wheeler_recurrence(qmomkit::MomentSet{1.0, 0.8, 2.8, 4.4, nan}, 2).ok());
}

void test_recurrence_entry_zero_and_lengths()
{
    const std::vector<double> points{-1.25, 0.5, 1.0, 2.75};
    const std::vector<double> weights{0.25, 0.375, 0.125, 0.25};
    const qmomkit::MomentSet moments = moments_from_nodes(points, weights);
    for (std::size_t n = 1; n <= points.size(); ++n)
    {
        const qmomkit::WheelerRecurrence recurrence = qmomkit::compute_wheeler_recurrence(moments, n);
        assert(recurrence.ok());
        assert(recurrence.message.empty());
        assert(recurrence.a.size() == n);
        assert(recurrence.b.size() == n);
        assert(recurrence.sigma.size() == n);
        assert(recurrence.tolerance.size() == n);
        assert(recurrence.a[0] == moments[1] / moments[0]);
        assert(recurrence.b[0] == 0.0);
        assert(recurrence.sigma[0] == moments[0]);
        assert(recurrence.tolerance[0] == 0.0);
    }
}

void test_recurrence_sigma_is_m0_times_the_product_of_b()
{
    // sigma(k, k) = m0 b_1 ... b_k, up to the rounding of the product.
    const qmomkit::MomentSet moments = moments_from_nodes({-1.25, 0.5, 1.0, 2.75}, {0.25, 0.375, 0.125, 0.25});
    const qmomkit::WheelerRecurrence recurrence = qmomkit::compute_wheeler_recurrence(moments, 4);
    assert(recurrence.ok());
    double product = moments[0];
    for (std::size_t k = 1; k < 4; ++k)
    {
        product *= recurrence.b[k];
        assert(recurrence.sigma[k] > 0.0);
        assert(std::abs(recurrence.sigma[k] - product) <= 1.0e-14 * std::abs(product));
    }
}

void test_recurrence_tolerance_is_the_documented_bound()
{
    // tau_k = 100 epsilon sigma-bar(k, k) at coordinate scale zero, as the
    // documented formula evaluates it independently. At scale 1e-3 the bound
    // lies far below m0, so a floor m0 s^(2k) with s of order one would show.
    const std::vector<double> points{-1.25, 0.5, 1.0, 2.75};
    const std::vector<double> weights{0.25, 0.375, 0.125, 0.25};
    for (const double scale : {1.0e-3, 1.0, 1.0e3})
    {
        std::vector<double> scaled = points;
        for (double& point : scaled)
        {
            point *= scale;
        }
        const qmomkit::MomentSet moments = moments_from_nodes(scaled, weights);
        const qmomkit::WheelerRecurrence recurrence = qmomkit::compute_wheeler_recurrence(moments, points.size());
        assert(recurrence.ok());
        for (std::size_t k = 1; k < points.size(); ++k)
        {
            const double expected = documented_rounding_bound(moments.values(), k);
            assert(std::abs(recurrence.tolerance[k] - expected) <= 1.0e-12 * expected);
        }
    }
}

// The usable count of adaptive inversion: the first k >= 1 with
// sigma(k, k) <= tau_k, or n if there is none. For distinct points of positive
// weight and rmin = eabs = 0, it is the count returned.
std::size_t count_from_recurrence(const qmomkit::WheelerRecurrence& recurrence)
{
    for (std::size_t k = 1; k < recurrence.sigma.size(); ++k)
    {
        if (recurrence.sigma[k] <= recurrence.tolerance[k])
        {
            return k;
        }
    }
    return recurrence.sigma.size();
}

// The k of a fixed-count refusal: the first k >= 1 with a finite tau_k and
// |sigma(k, k)| <= tau_k, or zero if there is none.
std::size_t refusal_from_recurrence(const qmomkit::WheelerRecurrence& recurrence)
{
    for (std::size_t k = 1; k < recurrence.sigma.size(); ++k)
    {
        if (std::isfinite(recurrence.tolerance[k]) && std::abs(recurrence.sigma[k]) <= recurrence.tolerance[k])
        {
            return k;
        }
    }
    return 0;
}

void test_recurrence_decides_the_node_count()
{
    // Two nodes merge, alone (two requested, sigma(1, 1) crosses tau_1) or
    // beside a third (three requested, sigma(2, 2) crosses tau_2), with the
    // separation shrinking through the threshold. The adaptive count and the
    // fixed-count refusal must follow the recurrence on both sides of it, for
    // each configuration.
    qmomkit::WheelerOptions adaptive;
    adaptive.rmin = 0.0;
    adaptive.eabs = 0.0;
    qmomkit::WheelerOptions fixed;
    fixed.adaptive = false;

    for (const auto& [base_points, weights] : {std::pair<std::vector<double>, std::vector<double>>{{1.0, 1.0}, {0.4, 0.6}},
                                               std::pair<std::vector<double>, std::vector<double>>{{-1.0, 1.0, 1.0}, {0.3, 0.3, 0.4}}})
    {
        bool merged = false;
        bool separated = false;
        for (int exponent = 1; exponent <= 12; ++exponent)
        {
            std::vector<double> points = base_points;
            points.back() += std::pow(10.0, -exponent);
            const std::size_t n = points.size();
            const qmomkit::MomentSet moments = moments_from_nodes(points, weights);
            const qmomkit::WheelerRecurrence recurrence = qmomkit::compute_wheeler_recurrence(moments, n);
            assert(recurrence.ok());

            const std::size_t count = count_from_recurrence(recurrence);
            if (count < n)
            {
                merged = true;
            }
            else
            {
                separated = true;
            }
            const qmomkit::Quadrature1D adaptive_rule = qmomkit::compute_quadrature(moments, n, adaptive);
            assert(adaptive_rule.ok());
            assert(adaptive_rule.weights.size() == count);

            const std::size_t refusal = refusal_from_recurrence(recurrence);
            const qmomkit::Quadrature1D fixed_rule = qmomkit::compute_quadrature(moments, n, fixed);
            if (refusal == 0)
            {
                assert(fixed_rule.ok());
                assert(fixed_rule.weights.size() == n);
            }
            else
            {
                const std::string expected = "Moment sequence supports at most " + std::to_string(refusal) + (refusal == 1 ? " node:" : " nodes:");
                assert(fixed_rule.status == qmomkit::QuadratureStatus::not_realizable);
                assert(fixed_rule.message.rfind(expected, 0) == 0);
            }
        }
        assert(merged);
        assert(separated);
    }
}

// The first count moments of the weighted points, exact in double for short
// binary fractions.
qmomkit::MomentSet moments_of(const std::vector<double>& abscissas, const std::vector<double>& weights, const std::size_t count)
{
    return qmomkit::MomentSet(reconstruct_moments(qmomkit::Quadrature1D{abscissas, weights}, count));
}

void assert_identical(const qmomkit::Quadrature1D& actual, const qmomkit::Quadrature1D& expected)
{
    assert(actual.abscissas == expected.abscissas);
    assert(actual.weights == expected.weights);
    assert(actual.status == expected.status);
    assert(actual.message == expected.message);
}

void test_hyperbolic_reports_input_failures()
{
    using qmomkit::QuadratureStatus;
    using qmomkit::internal::compute_hyperbolic_quadrature;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const qmomkit::MomentSet moments{1.0, 0.25, 1.0, 0.5, 2.0};
    const qmomkit::WheelerOptions options;

    qmomkit::WheelerOptions invalid;
    invalid.rmin = -1.0;
    const qmomkit::Quadrature1D invalid_options = compute_hyperbolic_quadrature(moments, 2, invalid);
    assert(invalid_options.status == QuadratureStatus::invalid_input);
    assert(invalid_options.message == "Invalid Wheeler numerical option.");
    assert(invalid_options.abscissas.empty() && invalid_options.weights.empty());

    const std::vector<std::pair<qmomkit::Quadrature1D, std::string>> invalid_inputs{
        {compute_hyperbolic_quadrature(moments, 0, options), "node_count must be positive."},
        {compute_hyperbolic_quadrature(moments, std::numeric_limits<std::size_t>::max(), options),
         "node_count is too large to derive the required moment count."},
        {compute_hyperbolic_quadrature(qmomkit::MomentSet{1.0, 0.25, 1.0, 0.5}, 2, options), "Moment sequence has fewer than 2*n+1 entries."},
        {compute_hyperbolic_quadrature(qmomkit::MomentSet{1.0, 0.25, 1.0, 0.5, nan}, 2, options), "Moment sequence contains non-finite values."},
    };
    for (const auto& [result, message] : invalid_inputs)
    {
        assert(result.status == QuadratureStatus::invalid_input);
        assert(result.message == message);
        assert(result.abscissas.empty() && result.weights.empty());
    }

    // m0 < 0 and m0 < cutoff fail as in compute_quadrature, with its diagnostic
    // vectors.
    for (const qmomkit::MomentSet& failing : {qmomkit::MomentSet{-1.0, 0.25, 1.0, 0.5, 2.0}, qmomkit::MomentSet{1.0e-31, 0.0, 1.0e-31, 0.0, 1.0e-31}})
    {
        const qmomkit::Quadrature1D result = compute_hyperbolic_quadrature(failing, 2, options);
        assert(result.status == QuadratureStatus::not_realizable);
        assert_identical(result, qmomkit::compute_quadrature(failing, 2, options));
    }

    // m_{2n+1} is not read.
    assert(compute_hyperbolic_quadrature(qmomkit::MomentSet{1.0, 0.25, 1.0, 0.5, 2.0, nan}, 2, options).ok());
}

void test_hyperbolic_on_the_boundary_returns_the_classical_rule()
{
    // k <= n weighted points put m_0, ..., m_{2n} on the boundary at k. Adaptive
    // mode returns exactly the classical selection for n nodes, which still
    // applies rmin: from three points on, the weight ratio 1/3 is below
    // rmin = 0.5, and the selection reduces. Fixed mode refuses at k.
    const std::vector<double> points{-1.25, 0.5, 1.0, 2.75};
    const std::vector<double> weights{0.25, 0.375, 0.125, 0.25};
    qmomkit::WheelerOptions strict;
    strict.rmin = 0.5;
    qmomkit::WheelerOptions fixed;
    fixed.adaptive = false;
    for (std::size_t n = 2; n <= 4; ++n)
    {
        for (std::size_t k = 1; k <= n; ++k)
        {
            const std::vector<double> x(points.begin(), points.begin() + static_cast<std::ptrdiff_t>(k));
            const std::vector<double> w(weights.begin(), weights.begin() + static_cast<std::ptrdiff_t>(k));
            const qmomkit::MomentSet moments = moments_of(x, w, 2 * n + 1);
            const qmomkit::Quadrature1D result = qmomkit::internal::compute_hyperbolic_quadrature(moments, n, {});
            assert(result.ok());
            assert(result.weights.size() == k);
            assert_identical(result, qmomkit::compute_quadrature(moments, n, {}));

            const qmomkit::Quadrature1D reduced = qmomkit::internal::compute_hyperbolic_quadrature(moments, n, strict);
            assert(reduced.ok());
            assert_identical(reduced, qmomkit::compute_quadrature(moments, n, strict));
            if (k >= 3)
            {
                assert(reduced.weights.size() < k);
            }

            const qmomkit::Quadrature1D refused = qmomkit::internal::compute_hyperbolic_quadrature(moments, n, fixed);
            const std::string expected = "Moment sequence supports at most " + std::to_string(k) + (k == 1 ? " node:" : " nodes:");
            assert(refused.status == qmomkit::QuadratureStatus::not_realizable);
            assert(refused.message.rfind(expected, 0) == 0);
        }
    }
}

void test_hyperbolic_below_the_boundary_returns_the_classical_rule()
{
    // sigma(n, n) depends on m_{2n} with coefficient one, so lowering m_{2n} by
    // 2 sigma(n, n) makes it -sigma(n, n), far below -tau_n. Adaptive mode
    // returns the classical rule of m_0, ..., m_{2n-1}; fixed mode reports the
    // negative b_n.
    const std::vector<double> points{-1.25, 0.5, 1.0, 2.75};
    const std::vector<double> weights{0.25, 0.375, 0.125, 0.25};
    qmomkit::WheelerOptions fixed;
    fixed.adaptive = false;
    for (std::size_t n = 1; n <= 3; ++n)
    {
        std::vector<double> values = moments_of(points, weights, 2 * n + 2).values();
        values[2 * n] -= 2.0 * qmomkit::compute_wheeler_recurrence(qmomkit::MomentSet(values), n + 1).sigma[n];
        const qmomkit::MomentSet lowered(values);
        const qmomkit::WheelerRecurrence recurrence = qmomkit::compute_wheeler_recurrence(lowered, n + 1);
        assert(recurrence.sigma[n] < -recurrence.tolerance[n]);

        const qmomkit::Quadrature1D result = qmomkit::internal::compute_hyperbolic_quadrature(lowered, n, {});
        assert(result.ok());
        assert_identical(result, qmomkit::compute_quadrature(lowered, n, {}));

        const qmomkit::Quadrature1D refused = qmomkit::internal::compute_hyperbolic_quadrature(lowered, n, fixed);
        assert(refused.status == qmomkit::QuadratureStatus::not_realizable);
        assert(refused.message.rfind("Moment sequence has negative recurrence b = ", 0) == 0);
    }
}

void test_hyperbolic_overflow_of_the_extra_row()
{
    // With finite moments {1, 1e200, 1e308}, sigma(1, 1) = m2 - a0 m1 overflows to
    // -inf, and so does tau_1 to inf. Adaptive mode then finds sigma(1, 1) <=
    // tau_1 and returns the classical rule; fixed mode skips the non-finite
    // tau_1 and reports the non-finite recurrence, where the classical fixed
    // rule of n = 1 succeeds without reading m2.
    const qmomkit::MomentSet moments{1.0, 1.0e200, 1.0e308};
    const qmomkit::WheelerRecurrence recurrence = qmomkit::compute_wheeler_recurrence(qmomkit::MomentSet{1.0, 1.0e200, 1.0e308, 0.0}, 2);
    assert(std::isinf(recurrence.sigma[1]) && std::isinf(recurrence.tolerance[1]));

    const qmomkit::Quadrature1D result = qmomkit::internal::compute_hyperbolic_quadrature(moments, 1, {});
    assert(result.ok());
    assert_identical(result, qmomkit::compute_quadrature(moments, 1, {}));

    qmomkit::WheelerOptions fixed;
    fixed.adaptive = false;
    assert(qmomkit::compute_quadrature(moments, 1, fixed).ok());
    const qmomkit::Quadrature1D refused = qmomkit::internal::compute_hyperbolic_quadrature(moments, 1, fixed);
    assert(refused.status == qmomkit::QuadratureStatus::not_realizable);
    assert(refused.message == "Wheeler recurrence produced non-finite values.");
}

void test_hyperbolic_keeps_a_small_weight_in_the_interior()
{
    // Points -1 and 1 with weight 1/2 and 2 with weight 1e-12 are interior at
    // n = 2, sigma(2, 2) about 200 tau_2. The J_3 rule has a weight ratio of
    // about 2e-11, below the default rmin, and adaptive mode returns it as
    // fixed mode does.
    const qmomkit::MomentSet moments = moments_of({-1.0, 1.0, 2.0}, {0.5, 0.5, 1.0e-12}, 6);
    const qmomkit::WheelerRecurrence recurrence = qmomkit::compute_wheeler_recurrence(moments, 3);
    assert(recurrence.sigma[1] > recurrence.tolerance[1]);
    assert(recurrence.sigma[2] > 100.0 * recurrence.tolerance[2]);

    const qmomkit::WheelerOptions options;
    const qmomkit::Quadrature1D result = qmomkit::internal::compute_hyperbolic_quadrature(moments, 2, options);
    assert(result.ok());
    assert(result.weights.size() == 3);
    const auto [smallest, largest] = std::minmax_element(result.weights.begin(), result.weights.end());
    assert(*smallest / *largest < options.rmin);

    qmomkit::WheelerOptions fixed;
    fixed.adaptive = false;
    assert_identical(result, qmomkit::internal::compute_hyperbolic_quadrature(moments, 2, fixed));
}

// Whether a rule passes the weight-ratio and separation tests of adaptive
// inversion (docs/developer/qmom.md, "Adaptive reduction").
bool passes_acceptance(const qmomkit::Quadrature1D& rule, const qmomkit::WheelerOptions& options)
{
    const auto [smallest_weight, largest_weight] = std::minmax_element(rule.weights.begin(), rule.weights.end());
    const std::vector<double>& x = rule.abscissas;
    double separation = 0.0;
    if (x.size() == 2)
    {
        separation = std::abs(x[1] - x[0]) / std::max(std::abs(x[0]), std::abs(x[1]));
    }
    else
    {
        double smallest_gap = x[1] - x[0];
        for (std::size_t i = 2; i < x.size(); ++i)
        {
            smallest_gap = std::min(smallest_gap, x[i] - x[i - 1]);
        }
        separation = smallest_gap / (x.back() - x.front());
    }
    return *smallest_weight / *largest_weight > options.rmin && separation > options.eabs;
}

void test_hyperbolic_does_not_apply_rmin_or_eabs()
{
    // With rmin = 0.99 and eabs = 0.9, the classical rules of these moments
    // reduce, and the J_{n+1} rule of the same interior moments fails the
    // tests too, but is returned whole, as fixed mode returns it. The points
    // lie away from the origin, so that the two equal weights of J_2 still
    // fail the two-node separation test.
    const std::vector<double> points{6.75, 8.5, 9.0, 10.75, 11.0};
    const std::vector<double> weights{0.25, 0.375, 0.125, 0.25, 0.0625};
    qmomkit::WheelerOptions strict;
    strict.rmin = 0.99;
    strict.eabs = 0.9;
    qmomkit::WheelerOptions fixed;
    fixed.adaptive = false;
    for (std::size_t n = 1; n <= 3; ++n)
    {
        const qmomkit::MomentSet moments = moments_of(points, weights, 2 * n + 2);
        assert(qmomkit::compute_quadrature(moments, n + 1, strict).weights.size() < n + 1);

        const qmomkit::Quadrature1D result = qmomkit::internal::compute_hyperbolic_quadrature(moments, n, strict);
        assert(result.ok());
        assert(result.weights.size() == n + 1);
        assert(!passes_acceptance(result, strict));
        assert_identical(result, qmomkit::internal::compute_hyperbolic_quadrature(moments, n, fixed));
    }
}

void test_hyperbolic_n_equal_one_gives_two_nodes_or_refuses_zero_variance()
{
    // Two nodes at mu -+ s with weights m0 / 2 have a_1 = a_0, so J_2 returns
    // them; zero variance is the boundary at one node.
    const qmomkit::MomentSet moments = moments_of({0.25, 1.25}, {2.5, 2.5}, 3);
    qmomkit::WheelerOptions fixed;
    fixed.adaptive = false;
    for (const qmomkit::WheelerOptions& options : {qmomkit::WheelerOptions{}, fixed})
    {
        const qmomkit::Quadrature1D result = qmomkit::internal::compute_hyperbolic_quadrature(moments, 1, options);
        assert(result.ok());
        assert_vector_near(result.abscissas, {0.25, 1.25}, 1.0e-15);
        assert_vector_near(result.weights, {2.5, 2.5}, 1.0e-14);
    }

    const qmomkit::MomentSet zero_variance{2.0, 1.5, 1.125};
    const qmomkit::Quadrature1D one_node = qmomkit::internal::compute_hyperbolic_quadrature(zero_variance, 1, {});
    assert(one_node.ok());
    assert_identical(one_node, qmomkit::compute_quadrature(zero_variance, 1, {}));
    const qmomkit::Quadrature1D refused = qmomkit::internal::compute_hyperbolic_quadrature(zero_variance, 1, fixed);
    assert(refused.status == qmomkit::QuadratureStatus::not_realizable);
    assert(refused.message.rfind("Moment sequence supports at most 1 node:", 0) == 0);
}

// The standardized moments S_3, S_4 and S_5 of a rule.
std::vector<double> standardized_moments(const qmomkit::Quadrature1D& rule)
{
    double mass = 0.0;
    double mean = 0.0;
    for (std::size_t i = 0; i < rule.weights.size(); ++i)
    {
        mass += rule.weights[i];
        mean += rule.weights[i] * rule.abscissas[i];
    }
    mean /= mass;
    std::vector<double> central(6, 0.0);
    for (std::size_t order = 2; order < central.size(); ++order)
    {
        for (std::size_t i = 0; i < rule.weights.size(); ++i)
        {
            central[order] += rule.weights[i] * std::pow(rule.abscissas[i] - mean, static_cast<double>(order)) / mass;
        }
    }
    return {central[3] / std::pow(central[2], 1.5), central[4] / (central[2] * central[2]), central[5] / std::pow(central[2], 2.5)};
}

void test_hyperbolic_rule_closes_the_next_moment_as_fox_and_laurent_state()
{
    // The closure gives S_3 = 0 for n = 1 (Fox & Laurent 2022, Section 5.1) and
    // S_5 = S_3 (5 S_4 - 3 S_3^2 - 1) / 2 for n = 2 (Section 5.2). The rules
    // reproduce m_0, ..., m_{2n}.
    const std::vector<double> points{-1.25, 0.5, 1.0, 2.75};
    const std::vector<double> weights{0.25, 0.375, 0.125, 0.25};
    for (std::size_t n = 1; n <= 2; ++n)
    {
        const qmomkit::MomentSet moments = moments_of(points, weights, 2 * n + 1);
        const qmomkit::Quadrature1D rule = qmomkit::internal::compute_hyperbolic_quadrature(moments, n, {});
        assert(rule.ok());
        assert(rule.weights.size() == n + 1);
        assert_vector_near_relative(reconstruct_moments(rule, 2 * n + 1), moments.values(), 1.0e-13);

        const std::vector<double> standardized = standardized_moments(rule);
        if (n == 1)
        {
            assert_near(standardized[0], 0.0, 1.0e-13);
        }
        else
        {
            const double s3 = standardized[0];
            const double s4 = standardized[1];
            assert_near_relative(standardized[2], 0.5 * s3 * (5.0 * s4 - 3.0 * s3 * s3 - 1.0), 1.0e-12);
        }
    }
}

void test_hyperbolic_speeds_are_covariant_for_a_narrow_distribution_far_from_the_origin()
{
    // The speeds of the points shift + scale y are shift + scale times those of
    // the points y. Six points of spread 2^-10 around 4 make every Lanczos step
    // cancel most of diag(x) q_k; removing the projections twice keeps the
    // speeds accurate to the rounding of the shift, 1e-12 of the spread, where
    // removing them once loses the spread entirely.
    const std::vector<double> points{-1.5, -0.5, 0.0, 0.5, 1.0, 2.0};
    const std::vector<double> weights{0.125, 0.25, 0.125, 0.25, 0.125, 0.125};
    const double shift = 4.0;
    const double scale = std::ldexp(1.0, -10);
    std::vector<double> narrow;
    for (const double point : points)
    {
        narrow.push_back(shift + (scale * point));
    }
    const qmomkit::internal::HyperbolicCharacteristicSpeeds reference = qmomkit::internal::compute_hyperbolic_characteristic_speeds(points, weights);
    const qmomkit::internal::HyperbolicCharacteristicSpeeds speeds = qmomkit::internal::compute_hyperbolic_characteristic_speeds(narrow, weights);
    assert(reference.status == qmomkit::QuadratureStatus::success && speeds.status == qmomkit::QuadratureStatus::success);
    assert(speeds.speeds.size() == 11 && reference.speeds.size() == 11);
    for (std::size_t i = 0; i < speeds.speeds.size(); ++i)
    {
        assert_near(speeds.speeds[i], shift + (scale * reference.speeds[i]), 1.0e-10 * scale);
    }
}

void test_hyperbolic_speeds_report_a_recurrence_that_vanishes_or_overflows()
{
    // Three points at zero give b_1 = 0 exactly at the first Lanczos step,
    // before the recurrence of n = 2 is complete; points near 1e200 overflow
    // b_1.
    const qmomkit::internal::HyperbolicCharacteristicSpeeds vanishing =
        qmomkit::internal::compute_hyperbolic_characteristic_speeds({0.0, 0.0, 0.0}, {0.25, 0.5, 0.25});
    assert(vanishing.status == qmomkit::QuadratureStatus::not_realizable);
    assert(vanishing.message == "Characteristic speed recurrence of the nodes vanished.");
    assert(vanishing.speeds.empty());

    const qmomkit::internal::HyperbolicCharacteristicSpeeds overflowing =
        qmomkit::internal::compute_hyperbolic_characteristic_speeds({1.0e200, 2.0e200}, {0.5, 0.5});
    assert(overflowing.status == qmomkit::QuadratureStatus::not_realizable);
    assert(overflowing.message == "Characteristic speed recurrence produced non-finite values.");
    assert(overflowing.speeds.empty());
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
            {"fixed_count_refuses_sequences_supporting_fewer_nodes", test_fixed_count_refuses_sequences_supporting_fewer_nodes},
            {"fixed_count_reports_the_documented_rounding_bound", test_fixed_count_reports_the_documented_rounding_bound},
            {"invalid_wheeler_options_fail_before_inversion", test_invalid_wheeler_options_fail_before_inversion},
            {"recurrence_reports_input_failures", test_recurrence_reports_input_failures},
            {"recurrence_entry_zero_and_lengths", test_recurrence_entry_zero_and_lengths},
            {"recurrence_sigma_is_m0_times_the_product_of_b", test_recurrence_sigma_is_m0_times_the_product_of_b},
            {"recurrence_tolerance_is_the_documented_bound", test_recurrence_tolerance_is_the_documented_bound},
            {"recurrence_decides_the_node_count", test_recurrence_decides_the_node_count},
            {"hyperbolic_reports_input_failures", test_hyperbolic_reports_input_failures},
            {"hyperbolic_on_the_boundary_returns_the_classical_rule", test_hyperbolic_on_the_boundary_returns_the_classical_rule},
            {"hyperbolic_below_the_boundary_returns_the_classical_rule", test_hyperbolic_below_the_boundary_returns_the_classical_rule},
            {"hyperbolic_overflow_of_the_extra_row", test_hyperbolic_overflow_of_the_extra_row},
            {"hyperbolic_keeps_a_small_weight_in_the_interior", test_hyperbolic_keeps_a_small_weight_in_the_interior},
            {"hyperbolic_does_not_apply_rmin_or_eabs", test_hyperbolic_does_not_apply_rmin_or_eabs},
            {"hyperbolic_n_equal_one_gives_two_nodes_or_refuses_zero_variance", test_hyperbolic_n_equal_one_gives_two_nodes_or_refuses_zero_variance},
            {"hyperbolic_rule_closes_the_next_moment_as_fox_and_laurent_state", test_hyperbolic_rule_closes_the_next_moment_as_fox_and_laurent_state},
            {"hyperbolic_speeds_are_covariant_for_a_narrow_distribution_far_from_the_origin",
             test_hyperbolic_speeds_are_covariant_for_a_narrow_distribution_far_from_the_origin},
            {"hyperbolic_speeds_report_a_recurrence_that_vanishes_or_overflows",
             test_hyperbolic_speeds_report_a_recurrence_that_vanishes_or_overflows},
        });
}
