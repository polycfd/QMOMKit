// Wheeler's inversion against sets of weighted points chosen in advance.
//
// What is verified
//   That compute_quadrature returns exactly the positions and weights the
//   moments were built from, for one to six and for eight nodes; that
//   multiplying every weight by a common factor leaves the positions unchanged
//   and scales the weights by that factor; that the same points multiplied by
//   powers of ten from 1e-20 to 1e20, or given in meters, are recovered as
//   accurately as at unit scale; that, for one to three nodes, reconstructing
//   the moments from the result returns the input; and that presenting the
//   same points in a different input order changes neither the moments nor the
//   sorted output.
//
// How the reference is constructed
//   Chosen, not computed. Each case starts from a list of positions and
//   weights, and the moments handed to the inversion are built from them by
//   construct_moments, which sums weight times position to the required power.
//   A set of n points is exactly determined by its first 2n moments, so the
//   inversion has one correct answer: the input points.
//
// Error metric and pass criteria
//   Maximum absolute deviation of each position and weight from the chosen
//   value, and of each rebuilt moment from the input moment. Tolerances are
//   per case: 1.0e-14 for the one-node and input-order cases and 1.0e-12 for
//   the two- and three-node references, whose values are order one; rebuilt
//   moments use 1.0e-13 for one node, 1.0e-12 for two and 1.0e-11 for three,
//   whose highest moment is about 67. Measured deviations are 1e-15 to 1e-13.
//
//   From four nodes on, positions are compared relative to the largest
//   |position| and weights relative to the total weight. The map from 2n
//   moments to n points is increasingly ill-conditioned as n grows, so the
//   rounding of the constructed input moments limits the attainable accuracy:
//   measured deviations are about 1e-15 for four nodes, 2e-14 for six and
//   1.6e-12 for eight. The tolerances, 1e-12 up to six nodes and 1e-10 for
//   eight, leave a margin of about fifty or more over those deviations while
//   still rejecting any misplaced point or weight.
//
//   The scaled points use the same relative metric. Three points at every
//   fourth power of ten from 1e-20 to 1e20 are recovered to within 5.2e-14,
//   against a tolerance of 1e-12; four droplet diameters between 20 and 200
//   micrometers, given in meters, to within 2.9e-13, against 2e-12, since the
//   spread of a factor ten makes that rule less well conditioned.
//
// Reproduce
//   ctest -R verification_wheeler_analytical
//
// Build modes
//   Both. QMOMKIT_VERIFY_CHECK and QMOMKIT_VERIFY_NEAR report and exit rather
//   than relying on assert(), so NDEBUG does not weaken this study.

#include "../verification_helpers.hpp"
#include "qmomkit/qmomkit.hpp"
#include "test_case_dispatch.hpp"

#include <algorithm>
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

void assert_sorted_ascending(const std::vector<double>& values)
{
    for (std::size_t i = 1; i < values.size(); ++i)
    {
        QMOMKIT_VERIFY_CHECK(values[i - 1] <= values[i]);
    }
}

std::vector<double> construct_moments(const std::vector<double>& abscissas, const std::vector<double>& weights, const std::size_t moment_count)
{
    QMOMKIT_VERIFY_CHECK(abscissas.size() == weights.size());
    std::vector<double> moments(moment_count, 0.0);
    for (std::size_t order = 0; order < moment_count; ++order)
    {
        for (std::size_t node = 0; node < weights.size(); ++node)
        {
            moments[order] += weights[node] * std::pow(abscissas[node], order);
        }
    }
    return moments;
}

std::vector<double> reconstruct_moments(const qmomkit::Quadrature1D& quadrature, const std::size_t moment_count)
{
    std::vector<double> moments(moment_count, 0.0);
    for (std::size_t order = 0; order < moment_count; ++order)
    {
        for (std::size_t node = 0; node < quadrature.weights.size(); ++node)
        {
            moments[order] += quadrature.weights[node] * std::pow(quadrature.abscissas[node], order);
        }
    }
    return moments;
}

void test_one_node_reference()
{
    const qmomkit::MomentSet moments{3.0, 7.5, 18.75, 46.875};
    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments, 1);

    QMOMKIT_VERIFY_CHECK(quadrature.ok());
    assert_vector_near(quadrature.abscissas, {2.5}, 1.0e-14);
    assert_vector_near(quadrature.weights, {3.0}, 1.0e-14);
}

void test_one_node_constructed_moment_reconstruction()
{
    const std::vector<double> expected_abscissas{-4.25};
    const std::vector<double> expected_weights{7.0};
    const qmomkit::MomentSet moments(construct_moments(expected_abscissas, expected_weights, 2));
    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments, 1);

    QMOMKIT_VERIFY_CHECK(quadrature.ok());
    assert_vector_near(quadrature.abscissas, expected_abscissas, 1.0e-14);
    assert_vector_near(quadrature.weights, expected_weights, 1.0e-14);
    assert_vector_near(reconstruct_moments(quadrature, moments.size()), moments.values(), 1.0e-13);
}

void test_two_node_reference()
{
    const qmomkit::MomentSet moments{1.0, 1.25, 3.25, 5.75};
    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments, 2);

    QMOMKIT_VERIFY_CHECK(quadrature.ok());
    assert_sorted_ascending(quadrature.abscissas);
    assert_vector_near(quadrature.abscissas, {-1.0, 2.0}, 1.0e-12);
    assert_vector_near(quadrature.weights, {0.25, 0.75}, 1.0e-12);
    assert_vector_near(reconstruct_moments(quadrature, moments.size()), moments.values(), 1.0e-12);
}

void test_three_node_reference()
{
    const qmomkit::MomentSet moments{1.0, 0.7499999999999999, 3.625, 6.5625, 27.53125, 66.51562499999999};
    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments, 3);

    QMOMKIT_VERIFY_CHECK(quadrature.ok());
    assert_sorted_ascending(quadrature.abscissas);
    assert_vector_near(quadrature.abscissas, {-2.0, 0.5, 3.0}, 1.0e-12);
    assert_vector_near(quadrature.weights, {0.2, 0.5, 0.3}, 1.0e-12);
    assert_vector_near(reconstruct_moments(quadrature, moments.size()), moments.values(), 1.0e-11);
}

void test_input_node_order_does_not_change_moments_or_sorted_output()
{
    const std::vector<double> ordered_moments = construct_moments({-2.0, 4.0}, {0.7, 1.3}, 4);
    const std::vector<double> shuffled_moments = construct_moments({4.0, -2.0}, {1.3, 0.7}, 4);
    assert_vector_near(ordered_moments, shuffled_moments, 1.0e-14);

    const qmomkit::Quadrature1D ordered = qmomkit::compute_quadrature(qmomkit::MomentSet(ordered_moments), 2);
    const qmomkit::Quadrature1D shuffled = qmomkit::compute_quadrature(qmomkit::MomentSet(shuffled_moments), 2);

    QMOMKIT_VERIFY_CHECK(ordered.ok());
    QMOMKIT_VERIFY_CHECK(shuffled.ok());
    assert_sorted_ascending(ordered.abscissas);
    assert_sorted_ascending(shuffled.abscissas);
    assert_vector_near(ordered.abscissas, {-2.0, 4.0}, 1.0e-14);
    assert_vector_near(ordered.weights, {0.7, 1.3}, 1.0e-14);
    assert_vector_near(shuffled.abscissas, ordered.abscissas, 1.0e-14);
    assert_vector_near(shuffled.weights, ordered.weights, 1.0e-14);
}

// Compares a recovered quadrature with the chosen points: positions relative to
// the largest |position|, weights relative to the total weight.
void assert_recovers_points(const std::vector<double>& abscissas, const std::vector<double>& weights, const double tolerance)
{
    const std::size_t n = abscissas.size();
    const qmomkit::MomentSet moments(construct_moments(abscissas, weights, 2 * n));
    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments, n);

    QMOMKIT_VERIFY_CHECK(quadrature.ok());
    QMOMKIT_VERIFY_CHECK(quadrature.abscissas.size() == n);
    assert_sorted_ascending(quadrature.abscissas);
    double position_scale = 0.0;
    double total_weight = 0.0;
    for (std::size_t node = 0; node < n; ++node)
    {
        position_scale = std::max(position_scale, std::abs(abscissas[node]));
        total_weight += weights[node];
    }
    for (std::size_t node = 0; node < n; ++node)
    {
        assert_near(quadrature.abscissas[node], abscissas[node], tolerance * position_scale);
        assert_near(quadrature.weights[node], weights[node], tolerance * total_weight);
    }
}

void test_four_node_reference()
{
    assert_recovers_points({-1.5, -0.5, 0.75, 2.0}, {0.1, 0.4, 0.3, 0.2}, 1.0e-12);
}

void test_five_node_reference()
{
    assert_recovers_points({-2.0, -1.0, 0.0, 1.5, 3.0}, {0.05, 0.2, 0.35, 0.25, 0.15}, 1.0e-12);
}

void test_six_node_reference()
{
    assert_recovers_points({-2.5, -1.5, -0.5, 0.5, 1.5, 2.5}, {0.08, 0.17, 0.25, 0.22, 0.18, 0.1}, 1.0e-12);
    assert_recovers_points({-0.9, -0.6, -0.2, 0.1, 0.5, 0.8}, {0.08, 0.17, 0.25, 0.22, 0.18, 0.1}, 1.0e-12);
}

// The same eight points with total weight 1 and with total weight 2500: the
// positions do not depend on the number density, and the weights scale with
// it.
void test_eight_node_reference_at_two_number_densities()
{
    const std::vector<double> abscissas{-0.95, -0.7, -0.4, -0.1, 0.15, 0.45, 0.7, 0.9};
    const std::vector<double> weights{0.06, 0.1, 0.14, 0.2, 0.18, 0.14, 0.11, 0.07};
    assert_recovers_points(abscissas, weights, 1.0e-10);

    std::vector<double> scaled_weights = weights;
    for (double& weight : scaled_weights)
    {
        weight *= 2500.0;
    }
    assert_recovers_points(abscissas, scaled_weights, 1.0e-10);
}

// The inversion must not depend on the unit of the coordinate: the same points
// multiplied by any power of ten are recovered with the same relative accuracy.
void test_points_at_any_scale_are_recovered()
{
    for (int exponent = -20; exponent <= 20; exponent += 4)
    {
        const double scale = std::pow(10.0, exponent);
        assert_recovers_points({scale, 2.0 * scale, 4.0 * scale}, {0.2, 0.5, 0.3}, 1.0e-12);
    }
}

// Droplet diameters between 20 and 200 micrometers, in meters, as an
// application provides them.
void test_droplet_diameters_in_metres_are_recovered()
{
    assert_recovers_points({20.0e-6, 50.0e-6, 100.0e-6, 200.0e-6}, {0.1, 0.4, 0.35, 0.15}, 2.0e-12);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"one_node_reference", test_one_node_reference},
            {"one_node_constructed_moment_reconstruction", test_one_node_constructed_moment_reconstruction},
            {"two_node_reference", test_two_node_reference},
            {"three_node_reference", test_three_node_reference},
            {"four_node_reference", test_four_node_reference},
            {"five_node_reference", test_five_node_reference},
            {"six_node_reference", test_six_node_reference},
            {"eight_node_reference_at_two_number_densities", test_eight_node_reference_at_two_number_densities},
            {"points_at_any_scale_are_recovered", test_points_at_any_scale_are_recovered},
            {"droplet_diameters_in_metres_are_recovered", test_droplet_diameters_in_metres_are_recovered},
            {"input_node_order_does_not_change_moments_or_sorted_output", test_input_node_order_does_not_change_moments_or_sorted_output},
        });
}
