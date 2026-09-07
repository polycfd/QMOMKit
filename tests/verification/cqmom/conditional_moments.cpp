#include "../verification_helpers.hpp"
#include "qmomkit/qmomkit.hpp"
#include "test_case_dispatch.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{

void assert_near(const double actual, const double expected, const double tolerance)
{
    QMOMKIT_VERIFY_NEAR(actual, expected, tolerance);
}

std::vector<double> make_mixed_moments(const std::vector<std::size_t>& extents, const std::size_t parent_dimension,
                                       const std::size_t target_dimension, const std::vector<double>& parent_abscissas,
                                       const std::vector<double>& parent_weights, const std::vector<std::vector<double>>& conditional_moments)
{
    QMOMKIT_VERIFY_CHECK(parent_abscissas.size() == parent_weights.size());
    QMOMKIT_VERIFY_CHECK(parent_abscissas.size() == conditional_moments.size());

    std::vector<double> moments(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat_index = 0; flat_index < moments.size(); ++flat_index)
    {
        const qmomkit::MultiIndex moment_index = qmomkit::unflatten_multi_index(flat_index, extents);
        const std::size_t parent_order = moment_index[parent_dimension];
        const std::size_t target_order = moment_index[target_dimension];

        double value = 0.0;
        for (std::size_t node = 0; node < parent_abscissas.size(); ++node)
        {
            value +=
                parent_weights[node] * std::pow(parent_abscissas[node], static_cast<int>(parent_order)) * conditional_moments[node][target_order];
        }
        moments[flat_index] = value;
    }
    return moments;
}

void assert_conditional_moments(const qmomkit::ConditionalMomentResult& result, const std::vector<std::vector<double>>& expected,
                                const double tolerance)
{
    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.parent_node_count == expected.size());
    QMOMKIT_VERIFY_CHECK(result.target_moment_count == expected.front().size());
    for (std::size_t node = 0; node < expected.size(); ++node)
    {
        for (std::size_t order = 0; order < expected[node].size(); ++order)
        {
            assert_near(result.moment(node, order), expected[node][order], tolerance);
        }
    }
}

void test_two_parent_nodes_conditional_orders_zero_to_three()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<double> parent_abscissas{-1.0, 2.0};
    const std::vector<double> parent_weights{0.25, 0.75};
    const std::vector<std::vector<double>> expected{{1.0, 3.0, 9.0, 27.0}, {1.0, -2.0, 4.0, -8.0}};
    const std::vector<std::size_t> extents{2, 4};

    const qmomkit::ConditionalMomentRequest request{
        hierarchy, 0, parent_abscissas, parent_weights, {extents, make_mixed_moments(extents, 0, 1, parent_abscissas, parent_weights, expected)}, {}};

    const qmomkit::ConditionalMomentResult result = qmomkit::reconstruct_conditional_moments(request);

    QMOMKIT_VERIFY_CHECK(result.parent_dimension == 0);
    QMOMKIT_VERIFY_CHECK(result.target_dimension == 1);
    assert_conditional_moments(result, expected, 1.0e-13);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.residual_norm.has_value());
    QMOMKIT_VERIFY_CHECK(*result.diagnostics.residual_norm <= 1.0e-14);
    QMOMKIT_VERIFY_CHECK(!result.diagnostics.conditioning_warning);
}

void test_analytically_known_2d_quadrature_reference_case()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<double> parent_abscissas{0.5, 1.5};
    const std::vector<double> parent_weights{0.4, 0.6};
    const std::vector<std::vector<double>> expected{{1.0, 1.25, 1.5625, 1.953125}, {1.0, -0.75, 0.5625, -0.421875}};
    const std::vector<std::size_t> extents{2, 4};

    const qmomkit::ConditionalMomentRequest request{
        hierarchy, 0, parent_abscissas, parent_weights, {extents, make_mixed_moments(extents, 0, 1, parent_abscissas, parent_weights, expected)}, {}};

    const qmomkit::ConditionalMomentResult result = qmomkit::reconstruct_conditional_moments(request);

    assert_conditional_moments(result, expected, 1.0e-13);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.residual_norm.has_value());
    QMOMKIT_VERIFY_CHECK(*result.diagnostics.residual_norm <= 1.0e-14);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"two_parent_nodes_conditional_orders_zero_to_three", test_two_parent_nodes_conditional_orders_zero_to_three},
            {"analytically_known_2d_quadrature_reference_case", test_analytically_known_2d_quadrature_reference_case},
        });
}
