#include "qmomkit/qmomkit.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{

void assert_near(const double actual, const double expected, const double tolerance)
{
    assert(std::abs(actual - expected) <= tolerance);
}

std::vector<double> make_mixed_moments(const std::vector<std::size_t>& extents, const std::size_t parent_dimension,
                                       const std::size_t target_dimension, const std::vector<double>& parent_abscissas,
                                       const std::vector<double>& parent_weights, const std::vector<std::vector<double>>& conditional_moments)
{
    assert(parent_abscissas.size() == parent_weights.size());
    assert(parent_abscissas.size() == conditional_moments.size());

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
    assert(result.ok());
    assert(result.parent_node_count == expected.size());
    assert(result.target_moment_count == expected.front().size());
    for (std::size_t node = 0; node < expected.size(); ++node)
    {
        for (std::size_t order = 0; order < expected[node].size(); ++order)
        {
            assert_near(result.moment(node, order), expected[node][order], tolerance);
        }
    }
}

void test_non_identity_conditioning_order_uses_phase_space_moment_axes()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 1, 2}, {2, 0, 1});
    const std::vector<double> parent_abscissas{-2.0, 0.25};
    const std::vector<double> parent_weights{0.3, 0.7};
    const std::vector<std::vector<double>> expected{{1.0, 4.0, 16.0, 64.0}, {1.0, -1.5, 2.25, -3.375}};
    const std::vector<std::size_t> extents{4, 1, 2};

    const qmomkit::ConditionalMomentRequest request{
        hierarchy, 0, parent_abscissas, parent_weights, {extents, {1.0, -0.425, 0.15, -2.6625, 6.375, -9.20625, 16.8375, -38.990625}}, {}};

    const qmomkit::ConditionalMomentResult result = qmomkit::reconstruct_conditional_moments(request);

    assert(result.parent_dimension == 2);
    assert(result.target_dimension == 0);
    assert_conditional_moments(result, expected, 1.0e-12);
    assert(result.diagnostics.residual_norm.has_value());
    assert(*result.diagnostics.residual_norm <= 1.0e-13);
}

void test_unused_phase_space_axes_are_read_at_zero_order()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 1, 2}, {2, 0, 1});
    const std::vector<double> parent_abscissas{-2.0, 0.25};
    const std::vector<double> parent_weights{0.3, 0.7};
    const std::vector<std::vector<double>> expected{{1.0, 4.0, 16.0, 64.0}, {1.0, -1.5, 2.25, -3.375}};
    const std::vector<std::size_t> extents{4, 3, 2};
    std::vector<double> raw_moments(qmomkit::product_extents(extents), 1000.0);

    for (std::size_t parent_order = 0; parent_order < parent_abscissas.size(); ++parent_order)
    {
        for (std::size_t target_order = 0; target_order < expected.front().size(); ++target_order)
        {
            double value = 0.0;
            for (std::size_t node = 0; node < parent_abscissas.size(); ++node)
            {
                value += parent_weights[node] * std::pow(parent_abscissas[node], static_cast<int>(parent_order)) * expected[node][target_order];
            }
            raw_moments[qmomkit::flatten_multi_index({target_order, 0, parent_order}, extents)] = value;
        }
    }

    const qmomkit::ConditionalMomentRequest request{hierarchy, 0, parent_abscissas, parent_weights, {extents, raw_moments}, {}};

    const qmomkit::ConditionalMomentResult result = qmomkit::reconstruct_conditional_moments(request);

    assert(result.parent_dimension == 2);
    assert(result.target_dimension == 0);
    assert_conditional_moments(result, expected, 1.0e-12);
}

void test_singular_parent_abscissas_fail()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{2, 4};
    const std::vector<double> parent_abscissas{1.0, 1.0};
    const std::vector<double> parent_weights{0.5, 0.5};
    const std::vector<std::vector<double>> expected{{1.0, 2.0, 4.0, 8.0}, {1.0, -1.0, 1.0, -1.0}};

    const qmomkit::ConditionalMomentRequest request{
        hierarchy, 0, parent_abscissas, parent_weights, {extents, make_mixed_moments(extents, 0, 1, parent_abscissas, parent_weights, expected)}, {}};

    const qmomkit::ConditionalMomentResult result = qmomkit::reconstruct_conditional_moments(request);

    assert(!result.ok());
    assert(result.status == qmomkit::ConditionalMomentStatus::singular_system);
}

void test_nearly_singular_parent_abscissas_warn()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{2, 4};
    const std::vector<double> parent_abscissas{1.0, 1.0 + 1.0e-10};
    const std::vector<double> parent_weights{0.5, 0.5};
    const std::vector<std::vector<double>> expected{{1.0, 0.5, 0.25, 0.125}, {1.0, -0.5, 0.25, -0.125}};
    qmomkit::ConditionalMomentOptions options;
    options.conditioning_warning_threshold = 1.0e8;

    const qmomkit::ConditionalMomentRequest request{
        hierarchy, 0, parent_abscissas, parent_weights, {extents, make_mixed_moments(extents, 0, 1, parent_abscissas, parent_weights, expected)},
        options};

    const qmomkit::ConditionalMomentResult result = qmomkit::reconstruct_conditional_moments(request);

    assert_conditional_moments(result, expected, 1.0e-6);
    assert(result.diagnostics.conditioning_warning);
    assert(result.diagnostics.condition_estimate.has_value());
    assert(*result.diagnostics.condition_estimate > options.conditioning_warning_threshold);
}

void test_invalid_inputs_report_failure()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    qmomkit::ConditionalMomentRequest request{hierarchy, 1, {0.0, 1.0}, {0.5, 0.5}, {{2, 4}, std::vector<double>(8, 0.0)}, {}};
    qmomkit::ConditionalMomentResult result = qmomkit::reconstruct_conditional_moments(request);
    assert(!result.ok());
    assert(result.status == qmomkit::ConditionalMomentStatus::invalid_input);

    request.conditioning_level = 0;
    request.parent_weights = {0.5};
    result = qmomkit::reconstruct_conditional_moments(request);
    assert(!result.ok());
    assert(result.status == qmomkit::ConditionalMomentStatus::invalid_input);

    request.parent_weights = {0.5, 0.5};
    request.raw_moments.extents = {1, 4};
    request.raw_moments.values.assign(4, 0.0);
    result = qmomkit::reconstruct_conditional_moments(request);
    assert(!result.ok());
    assert(result.status == qmomkit::ConditionalMomentStatus::invalid_input);

    request.conditioning_level = std::numeric_limits<std::size_t>::max();
    request.parent_weights = {0.5, 0.5};
    request.raw_moments.extents = {2, 4};
    request.raw_moments.values.assign(8, 0.0);
    result = qmomkit::reconstruct_conditional_moments(request);
    assert(!result.ok());
    assert(result.status == qmomkit::ConditionalMomentStatus::invalid_input);
}

} // namespace

int main()
{
    test_non_identity_conditioning_order_uses_phase_space_moment_axes();
    test_unused_phase_space_axes_are_read_at_zero_order();
    test_singular_parent_abscissas_fail();
    test_nearly_singular_parent_abscissas_warn();
    test_invalid_inputs_report_failure();

    return 0;
}
