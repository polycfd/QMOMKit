#include "../verification_helpers.hpp"
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

} // namespace

int main()
{
    test_one_node_reference();
    test_one_node_constructed_moment_reconstruction();
    test_two_node_reference();
    test_three_node_reference();
    test_input_node_order_does_not_change_moments_or_sorted_output();
    return 0;
}
