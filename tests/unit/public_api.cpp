#include "internal/performance.hpp"
#include "internal/reconstruction_scaling.hpp"
#include "qmomkit/qmomkit.hpp"
#include "test_case_dispatch.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace
{

void assert_near(const double actual, const double expected, const double tolerance)
{
    assert(std::abs(actual - expected) <= tolerance);
}

void assert_relative_near(const double actual, const double expected, const double tolerance)
{
    assert(std::abs(actual - expected) <= tolerance * std::max(1.0, std::abs(expected)));
}

bool contains_index(const std::vector<qmomkit::MultiIndex>& indices, const qmomkit::MultiIndex& index)
{
    return std::find(indices.begin(), indices.end(), index) != indices.end();
}

void append_unique(std::vector<qmomkit::MultiIndex>& indices, qmomkit::MultiIndex index)
{
    if (!contains_index(indices, index))
    {
        indices.push_back(std::move(index));
    }
}

void assert_indices_equal(const std::vector<qmomkit::MultiIndex>& actual, const std::vector<qmomkit::MultiIndex>& expected)
{
    assert(actual == expected);
}

void assert_unique_indices(const std::vector<qmomkit::MultiIndex>& indices)
{
    for (std::size_t i = 0; i < indices.size(); ++i)
    {
        for (std::size_t j = i + 1; j < indices.size(); ++j)
        {
            assert(indices[i] != indices[j]);
        }
    }
}

template <typename Function> void assert_throws_out_of_range(Function function)
{
    bool threw = false;
    try
    {
        function();
    }
    catch (const std::out_of_range&)
    {
        threw = true;
    }
    assert(threw);
}

template <typename Function> void assert_throws_invalid_argument(Function function)
{
    bool threw = false;
    try
    {
        function();
    }
    catch (const std::invalid_argument&)
    {
        threw = true;
    }
    assert(threw);
}

std::vector<double> raw_moments_1d(const std::vector<double>& abscissas, const std::vector<double>& weights, const std::size_t moment_count)
{
    std::vector<double> moments(moment_count, 0.0);
    for (std::size_t order = 0; order < moment_count; ++order)
    {
        for (std::size_t node = 0; node < abscissas.size(); ++node)
        {
            moments[order] += weights[node] * std::pow(abscissas[node], static_cast<int>(order));
        }
    }
    return moments;
}

void test_moment_set_public_container_behavior()
{
    const qmomkit::MomentSet empty;
    assert(empty.empty());
    assert(empty.size() == 0);
    assert(empty.values().empty());

    const qmomkit::MomentSet moments{1.0, -2.5, 4.25};
    assert(!moments.empty());
    assert(moments.size() == 3);
    assert(moments.values() == std::vector<double>({1.0, -2.5, 4.25}));
    assert_near(moments[0], 1.0, 0.0);
    assert_near(moments[2], 4.25, 0.0);
    assert_throws_out_of_range([&moments]() { (void)moments[3]; });
}

std::vector<double> indexed_values(const std::vector<std::size_t>& extents)
{
    std::vector<double> values(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat = 0; flat < values.size(); ++flat)
    {
        values[flat] = 10.0 + static_cast<double>(flat);
    }
    return values;
}

std::vector<double> indexed_count_values(const std::size_t count)
{
    std::vector<double> values(count, 0.0);
    for (std::size_t index = 0; index < values.size(); ++index)
    {
        values[index] = 10.0 + static_cast<double>(index);
    }
    return values;
}

std::vector<double> raw_moments_from_grid(const qmomkit::PhaseSpaceHierarchy& hierarchy, const std::vector<std::size_t>& extents,
                                          const std::vector<double>& weights, const std::vector<std::vector<double>>& abscissas)
{
    std::vector<double> moments(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat_moment = 0; flat_moment < moments.size(); ++flat_moment)
    {
        const qmomkit::MultiIndex moment_index = qmomkit::unflatten_multi_index(flat_moment, extents);
        double value = 0.0;
        for (const qmomkit::MultiIndex& node_index : qmomkit::enumerate_quadrature_node_multi_indices(hierarchy))
        {
            const std::size_t flat_node = hierarchy.flatten_node_index(node_index);
            double contribution = weights[flat_node];
            for (std::size_t dimension = 0; dimension < hierarchy.number_of_dimensions(); ++dimension)
            {
                contribution *= std::pow(abscissas[dimension][flat_node], static_cast<int>(moment_index[dimension]));
            }
            value += contribution;
        }
        moments[flat_moment] = value;
    }
    return moments;
}

std::vector<double> raw_moments_from_terminals(const std::vector<std::size_t>& extents, const std::vector<double>& weights,
                                               const std::vector<std::vector<double>>& abscissas)
{
    std::vector<double> moments(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat_moment = 0; flat_moment < moments.size(); ++flat_moment)
    {
        const qmomkit::MultiIndex moment_index = qmomkit::unflatten_multi_index(flat_moment, extents);
        for (std::size_t node = 0; node < weights.size(); ++node)
        {
            double contribution = weights[node];
            for (std::size_t dimension = 0; dimension < abscissas.size(); ++dimension)
            {
                contribution *= std::pow(abscissas[dimension][node], static_cast<int>(moment_index[dimension]));
            }
            moments[flat_moment] += contribution;
        }
    }
    return moments;
}

void assert_same_quadrature(const qmomkit::CqmomQuadrature& lhs, const qmomkit::CqmomQuadrature& rhs, const double tolerance)
{
    assert(lhs.number_of_dimensions() == rhs.number_of_dimensions());
    assert(lhs.total_nodes() == rhs.total_nodes());
    assert(lhs.adaptive_branch_counts_enabled() == rhs.adaptive_branch_counts_enabled());
    assert(lhs.has_rectangular_storage() == rhs.has_rectangular_storage());
    assert(lhs.terminal_paths() == rhs.terminal_paths());
    for (std::size_t node = 0; node < lhs.weights().size(); ++node)
    {
        assert_near(lhs.weights()[node], rhs.weights()[node], tolerance);
    }
    for (std::size_t dimension = 0; dimension < lhs.number_of_dimensions(); ++dimension)
    {
        for (std::size_t node = 0; node < lhs.total_nodes(); ++node)
        {
            assert_near(lhs.abscissas()[dimension][node], rhs.abscissas()[dimension][node], tolerance);
        }
    }
}

void assert_same_cqmom_result(const qmomkit::CqmomResult& lhs, const qmomkit::CqmomResult& rhs, const double tolerance)
{
    assert(lhs.status == rhs.status);
    assert(lhs.ok() == rhs.ok());
    assert(lhs.message == rhs.message);
    assert(lhs.failed_phase_space_dimension == rhs.failed_phase_space_dimension);
    assert(lhs.diagnostics.size() == rhs.diagnostics.size());
    for (std::size_t diagnostic = 0; diagnostic < lhs.diagnostics.size(); ++diagnostic)
    {
        assert(lhs.diagnostics[diagnostic].conditioning_level == rhs.diagnostics[diagnostic].conditioning_level);
        assert(lhs.diagnostics[diagnostic].parent_path == rhs.diagnostics[diagnostic].parent_path);
        assert(lhs.diagnostics[diagnostic].phase_space_dimension == rhs.diagnostics[diagnostic].phase_space_dimension);
        assert(lhs.diagnostics[diagnostic].requested_node_count == rhs.diagnostics[diagnostic].requested_node_count);
        assert(lhs.diagnostics[diagnostic].actual_node_count == rhs.diagnostics[diagnostic].actual_node_count);
        assert(lhs.diagnostics[diagnostic].adaptive_branch_counts_enabled == rhs.diagnostics[diagnostic].adaptive_branch_counts_enabled);
        assert(lhs.diagnostics[diagnostic].adaptive_node_count_reduced == rhs.diagnostics[diagnostic].adaptive_node_count_reduced);
        assert(lhs.diagnostics[diagnostic].qmom_status == rhs.diagnostics[diagnostic].qmom_status);
        assert(lhs.diagnostics[diagnostic].conditional_moment_status == rhs.diagnostics[diagnostic].conditional_moment_status);
        assert(lhs.diagnostics[diagnostic].conditional_moment_diagnostics.residual_norm.has_value() ==
               rhs.diagnostics[diagnostic].conditional_moment_diagnostics.residual_norm.has_value());
        assert(lhs.diagnostics[diagnostic].conditional_moment_diagnostics.condition_estimate.has_value() ==
               rhs.diagnostics[diagnostic].conditional_moment_diagnostics.condition_estimate.has_value());
        assert(lhs.diagnostics[diagnostic].conditional_moment_diagnostics.conditioning_warning ==
               rhs.diagnostics[diagnostic].conditional_moment_diagnostics.conditioning_warning);
    }
    assert_same_quadrature(lhs.quadrature, rhs.quadrature, tolerance);
}

qmomkit::CqmomMomentData lean_moments_from_dense(const qmomkit::PhaseSpaceHierarchy& hierarchy, const qmomkit::RawMomentTensor& dense)
{
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    return qmomkit::extract_cqmom_moments(dense, layout);
}

void assert_public_lean_cqmom_matches_dense(const qmomkit::PhaseSpaceHierarchy& hierarchy)
{
    const std::vector<std::size_t> extents = hierarchy.moment_extents();
    std::vector<double> weights(hierarchy.total_nodes(), 0.0);
    std::vector<std::vector<double>> abscissas(hierarchy.number_of_dimensions(), std::vector<double>(hierarchy.total_nodes(), 0.0));

    for (const qmomkit::MultiIndex& node_index : qmomkit::enumerate_quadrature_node_multi_indices(hierarchy))
    {
        const std::size_t flat_node = hierarchy.flatten_node_index(node_index);
        double weight = 1.0;
        for (std::size_t dimension = 0; dimension < hierarchy.number_of_dimensions(); ++dimension)
        {
            weight *= 0.25 + 0.06 * static_cast<double>(dimension + node_index[dimension] + 1);
            abscissas[dimension][flat_node] = -1.0 + 0.35 * static_cast<double>(dimension) + 0.8 * static_cast<double>(node_index[dimension]);
        }
        weights[flat_node] = weight;
    }

    const qmomkit::RawMomentTensor dense_moments{extents, raw_moments_from_grid(hierarchy, extents, weights, abscissas)};
    const qmomkit::CqmomResult dense = qmomkit::compute_cqmom({hierarchy, dense_moments, {}, {}});
    const qmomkit::CqmomMomentData lean_moments = lean_moments_from_dense(hierarchy, dense_moments);
    const qmomkit::CqmomResult lean = qmomkit::compute_cqmom({lean_moments, {}, {}});

    assert(dense.ok());
    assert(lean.ok());
    assert_same_cqmom_result(lean, dense, 1.0e-10);
}

void test_cqmom_moment_layout_two_dimensional_reduced_set()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);

    assert(layout.hierarchy.number_of_dimensions() == 2);
    assert(layout.hierarchy.node_counts() == std::vector<std::size_t>({2, 2}));
    assert(layout.hierarchy.conditioning_order() == std::vector<std::size_t>({0, 1}));
    assert_indices_equal(layout.required_indices,
                         std::vector<qmomkit::MultiIndex>({{0, 0}, {1, 0}, {2, 0}, {3, 0}, {0, 1}, {0, 2}, {0, 3}, {1, 1}, {1, 2}, {1, 3}}));
    assert_unique_indices(layout.required_indices);
}

void test_cqmom_moment_layout_three_dimensional_recursive_indices()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);

    std::vector<qmomkit::MultiIndex> expected;
    for (std::size_t i = 0; i < 4; ++i)
    {
        append_unique(expected, {i, 0, 0});
    }
    for (std::size_t i = 0; i < 2; ++i)
    {
        for (std::size_t j = 0; j < 4; ++j)
        {
            append_unique(expected, {i, j, 0});
        }
    }
    for (std::size_t i = 0; i < 2; ++i)
    {
        for (std::size_t j = 0; j < 2; ++j)
        {
            for (std::size_t k = 0; k < 4; ++k)
            {
                append_unique(expected, {i, j, k});
            }
        }
    }

    assert(expected.size() == 22);
    assert_indices_equal(layout.required_indices, expected);
    assert_unique_indices(layout.required_indices);
}

void test_cqmom_moment_layout_non_identity_conditioning_order()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2}, {2, 0, 1});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);

    std::vector<qmomkit::MultiIndex> expected;
    for (std::size_t order_dim2 = 0; order_dim2 < 4; ++order_dim2)
    {
        append_unique(expected, {0, 0, order_dim2});
    }
    for (std::size_t order_dim2 = 0; order_dim2 < 2; ++order_dim2)
    {
        for (std::size_t order_dim0 = 0; order_dim0 < 4; ++order_dim0)
        {
            append_unique(expected, {order_dim0, 0, order_dim2});
        }
    }
    for (std::size_t order_dim2 = 0; order_dim2 < 2; ++order_dim2)
    {
        for (std::size_t order_dim0 = 0; order_dim0 < 2; ++order_dim0)
        {
            for (std::size_t order_dim1 = 0; order_dim1 < 4; ++order_dim1)
            {
                append_unique(expected, {order_dim0, order_dim1, order_dim2});
            }
        }
    }

    assert(layout.hierarchy.conditioning_order() == std::vector<std::size_t>({2, 0, 1}));
    assert(expected.size() == 22);
    assert_indices_equal(layout.required_indices, expected);
    assert_unique_indices(layout.required_indices);
}

void test_cqmom_moment_layout_unequal_node_counts()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 3, 4});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);

    assert(layout.hierarchy.node_counts() == std::vector<std::size_t>({2, 3, 4}));
    assert(layout.required_indices.size() == 56);
    assert(contains_index(layout.required_indices, {3, 0, 0}));
    assert(!contains_index(layout.required_indices, {4, 0, 0}));
    assert(contains_index(layout.required_indices, {1, 5, 0}));
    assert(!contains_index(layout.required_indices, {2, 5, 0}));
    assert(contains_index(layout.required_indices, {1, 2, 7}));
    assert(!contains_index(layout.required_indices, {1, 3, 7}));
    assert_unique_indices(layout.required_indices);
}

void test_select_cqmom_moments_extracts_required_dense_values()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const std::vector<std::size_t> extents{4, 4};
    const qmomkit::RawMomentTensor dense{extents, indexed_values(extents)};

    const qmomkit::CqmomMomentData selected = qmomkit::extract_cqmom_moments(dense, layout);

    assert(selected.layout.required_indices == layout.required_indices);
    assert(selected.values.size() == layout.required_indices.size());
    for (std::size_t offset = 0; offset < layout.required_indices.size(); ++offset)
    {
        const qmomkit::MultiIndex& index = layout.required_indices[offset];
        const double expected = dense.values[qmomkit::flatten_multi_index(index, dense.extents)];
        assert(selected.values[offset] == expected);
        assert(selected.contains_required_index(index));
        assert(selected.moment(index) == expected);
    }
    assert(!selected.contains_required_index({2, 2}));
    assert_throws_out_of_range([&selected]() { (void)selected.moment({2, 2}); });
}

void test_cqmom_moment_data_lookup_reflects_public_value_updates()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const std::vector<std::size_t> extents{4, 4};
    qmomkit::CqmomMomentData selected = qmomkit::extract_cqmom_moments({extents, indexed_values(extents)}, layout);

    assert(selected.moment({1, 0}) == 14.0);
    selected.values[1] = 1234.0;
    assert(selected.moment({1, 0}) == 1234.0);
    assert(!selected.contains_required_index({0}));
}

void test_select_cqmom_moments_respects_non_identity_conditioning_order()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2}, {2, 0, 1});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const std::vector<std::size_t> extents{4, 4, 4};
    const qmomkit::RawMomentTensor dense{extents, indexed_values(extents)};

    const qmomkit::CqmomMomentData selected = qmomkit::extract_cqmom_moments(dense, layout);

    assert(selected.values.size() == 22);
    for (std::size_t offset = 0; offset < layout.required_indices.size(); ++offset)
    {
        const qmomkit::MultiIndex& index = layout.required_indices[offset];
        assert(selected.values[offset] == dense.values[qmomkit::flatten_multi_index(index, dense.extents)]);
    }
    assert(selected.moment({0, 0, 3}) == dense.values[qmomkit::flatten_multi_index({0, 0, 3}, dense.extents)]);
    assert(selected.moment({3, 0, 1}) == dense.values[qmomkit::flatten_multi_index({3, 0, 1}, dense.extents)]);
    assert(selected.moment({1, 3, 1}) == dense.values[qmomkit::flatten_multi_index({1, 3, 1}, dense.extents)]);
}

void test_select_cqmom_moments_handles_unequal_node_counts()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 3, 4});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const std::vector<std::size_t> extents{4, 6, 8};
    const qmomkit::RawMomentTensor dense{extents, indexed_values(extents)};

    const qmomkit::CqmomMomentData selected = qmomkit::extract_cqmom_moments(dense, layout);

    assert(selected.values.size() == 56);
    assert(selected.moment({3, 0, 0}) == dense.values[qmomkit::flatten_multi_index({3, 0, 0}, dense.extents)]);
    assert(selected.moment({1, 5, 0}) == dense.values[qmomkit::flatten_multi_index({1, 5, 0}, dense.extents)]);
    assert(selected.moment({1, 2, 7}) == dense.values[qmomkit::flatten_multi_index({1, 2, 7}, dense.extents)]);
}

void test_select_cqmom_moments_reports_missing_dense_moment()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const std::vector<std::size_t> extents{3, 4};
    const qmomkit::RawMomentTensor dense{extents, indexed_values(extents)};

    assert_throws_out_of_range([&dense, &layout]() { (void)qmomkit::extract_cqmom_moments(dense, layout); });
}

void test_create_cqmom_moment_data_converts_valid_selected_result()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::RawMomentRequest request{layout.required_indices};
    const qmomkit::SelectedRawMomentResult selected =
        qmomkit::make_selected_raw_moment_result(request, indexed_count_values(layout.required_indices.size()));

    const qmomkit::CqmomMomentData data = qmomkit::create_cqmom_moment_data(layout, selected);

    assert(data.layout.required_indices == layout.required_indices);
    assert(data.values == selected.values);
    assert(data.moment(layout.required_indices.front()) == selected.values.front());
}

void test_create_cqmom_moment_data_rejects_failed_selected_result()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    qmomkit::SelectedRawMomentResult selected =
        qmomkit::make_selected_raw_moment_result({layout.required_indices}, indexed_count_values(layout.required_indices.size()));
    selected.status = qmomkit::RawMomentStatus::invalid_input;
    selected.message = "failed upstream";

    assert_throws_invalid_argument([&layout, &selected]() { (void)qmomkit::create_cqmom_moment_data(layout, selected); });
}

void test_create_cqmom_moment_data_rejects_index_order_mismatch()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    qmomkit::RawMomentRequest request{layout.required_indices};
    std::swap(request.indices.front(), request.indices.back());
    const qmomkit::SelectedRawMomentResult selected =
        qmomkit::make_selected_raw_moment_result(request, indexed_count_values(layout.required_indices.size()));

    assert_throws_invalid_argument([&layout, &selected]() { (void)qmomkit::create_cqmom_moment_data(layout, selected); });
}

void test_create_cqmom_moment_data_rejects_missing_values()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    qmomkit::SelectedRawMomentResult selected =
        qmomkit::make_selected_raw_moment_result({layout.required_indices}, indexed_count_values(layout.required_indices.size()));
    selected.values.pop_back();

    assert_throws_invalid_argument([&layout, &selected]() { (void)qmomkit::create_cqmom_moment_data(layout, selected); });
}

void test_create_cqmom_moment_data_rejects_extra_values()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    qmomkit::SelectedRawMomentResult selected =
        qmomkit::make_selected_raw_moment_result({layout.required_indices}, indexed_count_values(layout.required_indices.size()));
    selected.values.push_back(123.0);

    assert_throws_invalid_argument([&layout, &selected]() { (void)qmomkit::create_cqmom_moment_data(layout, selected); });
}

void test_create_cqmom_moment_data_rejects_non_finite_values()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    qmomkit::SelectedRawMomentResult selected =
        qmomkit::make_selected_raw_moment_result({layout.required_indices}, indexed_count_values(layout.required_indices.size()));
    selected.values.front() = std::numeric_limits<double>::infinity();

    assert_throws_invalid_argument([&layout, &selected]() { (void)qmomkit::create_cqmom_moment_data(layout, selected); });
}

void test_create_cqmom_moment_data_from_independent_product_selected_moments()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::IndependentProductDistribution product{{qmomkit::GammaDistribution{1.5, 0.5}, qmomkit::BetaDistribution{2.0, 3.0}}, 2.0};

    const qmomkit::SelectedRawMomentResult selected = qmomkit::compute_raw_moments(product, qmomkit::RawMomentRequest{layout.required_indices});
    const qmomkit::CqmomMomentData data = qmomkit::create_cqmom_moment_data(layout, selected);
    const qmomkit::TensorMomentResult dense = qmomkit::compute_raw_moments(product, {3, 3});
    const qmomkit::CqmomMomentData required_moment_data = qmomkit::extract_cqmom_moments(dense.moments, layout);

    assert(selected.ok());
    assert(dense.ok());
    assert(data.layout.required_indices == layout.required_indices);
    assert(data.values.size() == required_moment_data.values.size());
    for (std::size_t index = 0; index < data.values.size(); ++index)
    {
        assert_near(data.values[index], required_moment_data.values[index], 1.0e-12);
    }
}

void test_lean_cqmom_public_api_matches_dense_two_dimensional()
{
    assert_public_lean_cqmom_matches_dense(qmomkit::PhaseSpaceHierarchy({2, 2}));
}

void test_lean_cqmom_public_api_matches_dense_three_dimensional()
{
    assert_public_lean_cqmom_matches_dense(qmomkit::PhaseSpaceHierarchy({2, 2, 2}));
}

void test_lean_cqmom_public_api_matches_dense_non_identity()
{
    assert_public_lean_cqmom_matches_dense(qmomkit::PhaseSpaceHierarchy({2, 2, 2}, {2, 0, 1}));
}

void test_lean_cqmom_public_api_matches_dense_unequal_node_counts()
{
    assert_public_lean_cqmom_matches_dense(qmomkit::PhaseSpaceHierarchy({2, 3, 4}));
}

void test_lean_cqmom_public_api_matches_dense_adaptive()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2});
    const std::vector<std::size_t> extents{4, 4, 4};
    const std::vector<double> weights{0.05, 0.05, 0.30, 0.60};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, -1.0, 1.0}, {-2.0, -2.0, 2.0, 3.0}, {-1.0, 1.0, 0.0, 2.0}};
    const qmomkit::RawMomentTensor dense_moments{extents, raw_moments_from_terminals(extents, weights, abscissas)};
    const qmomkit::QmomOptions qmom_options{true, 1.0e-10, 1.0e-10, 1.0e-12};

    qmomkit::CqmomDenseCase dense_case{hierarchy, dense_moments, {}, qmom_options, true};
    dense_case.diagnostic_level = qmomkit::CqmomDiagnosticLevel::basic;
    qmomkit::CqmomLeanCase lean_case{lean_moments_from_dense(hierarchy, dense_moments), {}, qmom_options, true};
    lean_case.diagnostic_level = qmomkit::CqmomDiagnosticLevel::basic;

    const qmomkit::CqmomResult dense = qmomkit::compute_cqmom(dense_case);
    const qmomkit::CqmomResult lean = qmomkit::compute_cqmom(lean_case);

    assert(dense.ok());
    assert(lean.ok());
    assert(dense.quadrature.adaptive_branch_counts_enabled());
    assert(!dense.quadrature.has_rectangular_storage());
    assert_same_cqmom_result(lean, dense, 1.0e-10);
}

void test_lean_cqmom_public_api_rejects_missing_values()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    qmomkit::CqmomMomentData moment_data = lean_moments_from_dense(hierarchy, {extents, indexed_values(extents)});
    moment_data.values.pop_back();

    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(qmomkit::CqmomLeanCase{moment_data});

    assert(!result.ok());
    assert(result.status == qmomkit::CqmomStatus::invalid_input);
    assert(result.message == "CQMOM moment data values must match the layout size.");
}

void test_lean_cqmom_public_api_rejects_non_finite_values()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    qmomkit::CqmomMomentData moment_data = lean_moments_from_dense(hierarchy, {extents, indexed_values(extents)});
    moment_data.values.front() = std::numeric_limits<double>::infinity();

    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(qmomkit::CqmomLeanCase{moment_data});

    assert(!result.ok());
    assert(result.status == qmomkit::CqmomStatus::invalid_input);
    assert(result.message == "CQMOM moment data values must be finite.");
}

void test_lean_cqmom_public_api_rejects_invalid_options()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const qmomkit::CqmomMomentData moment_data = lean_moments_from_dense(hierarchy, {extents, indexed_values(extents)});
    qmomkit::CqmomLeanCase lean_case{moment_data};
    lean_case.conditional_options.singular_tolerance = 0.0;

    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(lean_case);

    assert(!result.ok());
    assert(result.status == qmomkit::CqmomStatus::invalid_input);
    assert(!result.message.empty());
}

void test_lean_cqmom_public_api_rejects_incompatible_layout()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2});
    const qmomkit::PhaseSpaceHierarchy incompatible_hierarchy({2, 2, 2}, {2, 0, 1});
    const qmomkit::CqmomMomentData moment_data = lean_moments_from_dense(incompatible_hierarchy, {{4, 4, 4}, indexed_values({4, 4, 4})});
    qmomkit::CqmomLeanCase lean_case{moment_data};
    lean_case.moments.layout.hierarchy = hierarchy;

    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(lean_case);

    assert(!result.ok());
    assert(result.status == qmomkit::CqmomStatus::invalid_input);
    assert(result.message == "CQMOM moment data layout must match the recursive hierarchy.");
}

void test_cqmom_public_api_reports_missing_required_dense_moment()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{3, 4};
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom({hierarchy, {extents, indexed_values(extents)}, {}, {}});

    assert(!result.ok());
    assert(result.status == qmomkit::CqmomStatus::invalid_input);
    assert(!result.message.empty());
}

void test_cqmom_public_api_preserves_dense_invalid_input_messages()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});

    const qmomkit::CqmomResult rank =
        qmomkit::compute_cqmom(qmomkit::CqmomDenseCase{hierarchy, qmomkit::RawMomentTensor{{4}, std::vector<double>(4, 1.0)}, {}, {}});
    assert(!rank.ok());
    assert(rank.status == qmomkit::CqmomStatus::invalid_input);
    assert(rank.message == "Raw moment extents must match the hierarchy dimension.");

    const qmomkit::CqmomResult count =
        qmomkit::compute_cqmom(qmomkit::CqmomDenseCase{hierarchy, qmomkit::RawMomentTensor{{4, 4}, std::vector<double>(15, 1.0)}, {}, {}});
    assert(!count.ok());
    assert(count.status == qmomkit::CqmomStatus::invalid_input);
    assert(count.message == "Raw moment data size does not match its extents.");

    qmomkit::RawMomentTensor nonfinite{{4, 4}, std::vector<double>(16, 1.0)};
    nonfinite.values.front() = std::numeric_limits<double>::infinity();
    const qmomkit::CqmomResult finite = qmomkit::compute_cqmom(qmomkit::CqmomDenseCase{hierarchy, nonfinite, {}, {}});
    assert(!finite.ok());
    assert(finite.status == qmomkit::CqmomStatus::invalid_input);
    assert(finite.message == "Raw moments must be finite.");

    const qmomkit::CqmomResult extents =
        qmomkit::compute_cqmom(qmomkit::CqmomDenseCase{hierarchy, qmomkit::RawMomentTensor{{3, 4}, indexed_values({3, 4})}, {}, {}});
    assert(!extents.ok());
    assert(extents.status == qmomkit::CqmomStatus::invalid_input);
    assert(extents.message == "Raw moment extents must contain at least 2*N orders per dimension.");
}

template <typename Case> void assert_cqmom_rejects_unrepresentable_required_moment_extent(const Case& case_)
{
    bool leaked_length_error = false;
    try
    {
        const qmomkit::CqmomResult result = qmomkit::compute_cqmom(case_);
        assert(!result.ok());
        assert(result.status == qmomkit::CqmomStatus::invalid_input);
        assert(result.message == "CQMOM node count is too large to derive the required moment extent.");
    }
    catch (const std::length_error&)
    {
        leaked_length_error = true;
    }
    assert(!leaked_length_error);
}

void test_cqmom_public_api_rejects_unrepresentable_required_moment_extent_before_storage()
{
    const std::size_t node_count = std::numeric_limits<std::size_t>::max() / 2 + 1;
    const qmomkit::PhaseSpaceHierarchy hierarchy({node_count});

    const qmomkit::CqmomDenseCase dense_case{hierarchy, {{1}, {1.0}}, {}, {}};
    assert_cqmom_rejects_unrepresentable_required_moment_extent(dense_case);

    const qmomkit::CqmomMomentData lean_moments{{hierarchy, {}}, {}};
    const qmomkit::CqmomLeanCase lean_case{lean_moments};
    assert_cqmom_rejects_unrepresentable_required_moment_extent(lean_case);
}

void test_qmom_public_api_matches_wheeler()
{
    // Public-wrapper consistency: the unscaled public result must preserve the underlying Wheeler result; analytical Wheeler correctness is verified
    // separately.
    const std::vector<double> moments = raw_moments_1d({-0.5, 1.75}, {0.25, 0.75}, 4);

    const qmomkit::QmomResult qmom = qmomkit::compute_qmom(qmomkit::MomentSet(moments), 2);
    const qmomkit::Quadrature1D wheeler = qmomkit::compute_quadrature(qmomkit::MomentSet(moments), 2);

    assert(qmom.ok());
    assert(wheeler.ok());
    assert(qmom.abscissas.size() == 2);
    for (std::size_t node = 0; node < qmom.abscissas.size(); ++node)
    {
        assert_near(qmom.abscissas[node], wheeler.abscissas[node], 1.0e-13);
        assert_near(qmom.weights[node], wheeler.weights[node], 1.0e-13);
    }
}

void test_qmom_rejects_node_count_with_unrepresentable_moment_requirement()
{
    const std::size_t impossible_node_count = std::numeric_limits<std::size_t>::max() / 2 + 1;
    const qmomkit::QmomResult result = qmomkit::compute_qmom(qmomkit::MomentSet{1.0, 0.0}, impossible_node_count);

    assert(!result.ok());
    assert(result.status == qmomkit::QmomStatus::invalid_input);
    assert(result.message.find("too large") != std::string::npos);
}

void test_qmom_automatic_scaling_restores_physical_abscissas_and_bounds()
{
    const std::vector<double> abscissas{2.0e50, 5.0e50};
    const std::vector<double> weights{0.25, 0.75};
    const qmomkit::MomentSet moments(raw_moments_1d(abscissas, weights, 4));
    qmomkit::QmomOptions options;
    options.automatic_scaling = true;

    const qmomkit::QmomResult result = qmomkit::compute_qmom(moments, 2, options, qmomkit::AbscissaBounds{1.0e50, 6.0e50});

    assert(result.ok());
    assert(result.abscissas.size() == abscissas.size());
    for (std::size_t node = 0; node < abscissas.size(); ++node)
    {
        assert_relative_near(result.abscissas[node], abscissas[node], 1.0e-12);
        assert_near(result.weights[node], weights[node], 1.0e-12);
    }
}

void test_qmom_automatic_scaling_is_covariant_and_uses_existing_scale_moments()
{
    const std::vector<double> reference_abscissas{-0.5, 1.75};
    const std::vector<double> weights{0.25, 0.75};
    qmomkit::QmomOptions options;
    options.automatic_scaling = true;
    const qmomkit::QmomResult reference = qmomkit::compute_qmom(qmomkit::MomentSet(raw_moments_1d(reference_abscissas, weights, 4)), 2, options);
    assert(reference.ok());

    for (const double factor : std::vector<double>{1.0e-9, 1.0e-6, 1.0, 1.0e3, 1.0e6})
    {
        std::vector<double> scaled_abscissas = reference_abscissas;
        for (double& abscissa : scaled_abscissas)
        {
            abscissa *= factor;
        }
        const qmomkit::QmomResult result = qmomkit::compute_qmom(qmomkit::MomentSet(raw_moments_1d(scaled_abscissas, weights, 4)), 2, options);
        assert(result.ok());
        assert(result.weights.size() == reference.weights.size());
        for (std::size_t node = 0; node < result.weights.size(); ++node)
        {
            assert_near(result.weights[node], reference.weights[node], 1.0e-12);
            assert_near(result.abscissas[node] / factor, reference.abscissas[node], 1.0e-12);
        }
    }

    const qmomkit::QmomResult mean_fallback = qmomkit::compute_qmom(qmomkit::MomentSet({2.0, 6.0}), 1, options);
    assert(mean_fallback.ok());
    assert_near(mean_fallback.weights.front(), 2.0, 1.0e-13);
    assert_near(mean_fallback.abscissas.front(), 3.0, 1.0e-13);

    const qmomkit::QmomResult identity_fallback = qmomkit::compute_qmom(qmomkit::MomentSet({2.0, 0.0}), 1, options);
    assert(identity_fallback.ok());
    assert_near(identity_fallback.weights.front(), 2.0, 1.0e-13);
    assert_near(identity_fallback.abscissas.front(), 0.0, 1.0e-13);

    qmomkit::QmomOptions adaptive_options = options;
    adaptive_options.adaptive = true;
    const qmomkit::QmomResult adaptive_reference = qmomkit::compute_qmom(qmomkit::MomentSet(raw_moments_1d({0.25}, {1.0}, 4)), 2, adaptive_options);
    assert(adaptive_reference.ok());
    for (const double factor : std::vector<double>{1.0e-9, 1.0e6})
    {
        const qmomkit::QmomResult adaptive_result =
            qmomkit::compute_qmom(qmomkit::MomentSet(raw_moments_1d({0.25 * factor}, {1.0}, 4)), 2, adaptive_options);
        assert(adaptive_result.ok());
        assert(adaptive_result.weights.size() == adaptive_reference.weights.size());
        assert_near(adaptive_result.abscissas.front() / factor, adaptive_reference.abscissas.front(), 1.0e-12);
    }
}

void test_qmom_automatic_scaling_rejects_nonzero_underflow_before_reconstruction()
{
    const double denormal = std::numeric_limits<double>::denorm_min();
    const std::vector<double> moments{denormal, 1.0, 1.0, denormal};
    const auto scaling = qmomkit::internal::make_qmom_scaling_context(moments, std::nullopt);
    assert(!scaling.has_value());

    qmomkit::QmomOptions options;
    options.automatic_scaling = true;
    const qmomkit::QmomResult scaled = qmomkit::compute_qmom(qmomkit::MomentSet(moments), 2, options);
    options.automatic_scaling = false;
    const qmomkit::QmomResult unscaled = qmomkit::compute_qmom(qmomkit::MomentSet(moments), 2, options);
    assert(scaled.status == unscaled.status);
    assert(scaled.message == unscaled.message);
    assert(scaled.abscissas == unscaled.abscissas);
    assert(scaled.weights == unscaled.weights);
}

void test_cqmom_automatic_scaling_falls_back_atomically_when_preparation_is_unsafe()
{
    // Finite normal second moments produce 1e-100 coordinate scales. The
    // finite normal mixed M[1,1] then overflows when divided by both scales,
    // so preparation must be rejected before recursive CQMOM sees transformed
    // data. This is not a successful physical quadrature oracle; the contract
    // is that the scaled request returns the untouched unscaled failure result,
    // including physical diagnostics and empty quadrature.
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {1, 0});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    qmomkit::CqmomMomentData moments{layout, std::vector<double>(layout.required_indices.size(), 1.0)};
    const auto set_required = [&moments](const qmomkit::MultiIndex& index, const double value)
    {
        const auto found = std::find(moments.layout.required_indices.begin(), moments.layout.required_indices.end(), index);
        assert(found != moments.layout.required_indices.end());
        moments.values[static_cast<std::size_t>(std::distance(moments.layout.required_indices.begin(), found))] = value;
    };
    set_required({0, 0}, 1.0);
    set_required({2, 0}, 1.0e-200);
    set_required({0, 2}, 1.0e-200);
    set_required({1, 1}, 1.0e200);

    const auto scaling = qmomkit::internal::make_cqmom_scaling_context(moments, {});
    assert(!scaling.has_value());

    qmomkit::QmomOptions options;
    options.automatic_scaling = true;
    const qmomkit::CqmomResult scaled = qmomkit::compute_cqmom(qmomkit::CqmomLeanCase{moments, {}, options});
    options.automatic_scaling = false;
    const qmomkit::CqmomResult unscaled = qmomkit::compute_cqmom(qmomkit::CqmomLeanCase{moments, {}, options});

    assert(scaled.status == unscaled.status);
    assert(scaled.message == unscaled.message);
    assert(scaled.failed_phase_space_dimension == unscaled.failed_phase_space_dimension);
    assert(scaled.quadrature.weights() == unscaled.quadrature.weights());
    assert(scaled.quadrature.abscissas() == unscaled.quadrature.abscissas());
    assert(scaled.diagnostics.size() == unscaled.diagnostics.size());
    for (std::size_t diagnostic = 0; diagnostic < scaled.diagnostics.size(); ++diagnostic)
    {
        const qmomkit::CqmomDiagnostics& scaled_diagnostic = scaled.diagnostics[diagnostic];
        const qmomkit::CqmomDiagnostics& unscaled_diagnostic = unscaled.diagnostics[diagnostic];
        assert(scaled_diagnostic.conditioning_level == unscaled_diagnostic.conditioning_level);
        assert(scaled_diagnostic.parent_path == unscaled_diagnostic.parent_path);
        assert(scaled_diagnostic.phase_space_dimension == unscaled_diagnostic.phase_space_dimension);
        assert(scaled_diagnostic.requested_node_count == unscaled_diagnostic.requested_node_count);
        assert(scaled_diagnostic.actual_node_count == unscaled_diagnostic.actual_node_count);
        assert(scaled_diagnostic.qmom_status == unscaled_diagnostic.qmom_status);
        assert(scaled_diagnostic.conditional_moment_status == unscaled_diagnostic.conditional_moment_status);
        assert(scaled_diagnostic.message == unscaled_diagnostic.message);
        assert(scaled_diagnostic.bounded_returned_nodes == unscaled_diagnostic.bounded_returned_nodes);
        assert(scaled_diagnostic.bounded_returned_weights == unscaled_diagnostic.bounded_returned_weights);
    }
}

void test_one_dimensional_cqmom_agrees_with_qmom()
{
    const std::vector<double> moments = raw_moments_1d({-1.0, 2.0}, {0.4, 0.6}, 4);
    const qmomkit::QmomResult qmom = qmomkit::compute_qmom(qmomkit::MomentSet(moments), 2);

    const qmomkit::PhaseSpaceHierarchy hierarchy({2});
    const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom({hierarchy, {{4}, moments}, {}, {}});

    assert(qmom.ok());
    assert(cqmom.ok());
    assert(cqmom.quadrature.number_of_dimensions() == 1);
    assert(cqmom.quadrature.total_nodes() == 2);

    for (std::size_t node = 0; node < 2; ++node)
    {
        assert_near(cqmom.quadrature.weight({node}), qmom.weights[node], 1.0e-13);
        assert_near(cqmom.quadrature.abscissa(0, {node}), qmom.abscissas[node], 1.0e-13);
    }
}

void test_cqmom_public_api_two_dimensional_manufactured_case()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {1, 0});
    const std::vector<std::size_t> extents{4, 4};

    const std::vector<double> x{-0.75, 1.25};
    const std::vector<double> y{0.2, 1.4};
    const std::vector<double> x_weights{0.35, 0.65};
    const std::vector<double> y_weights{0.8, 0.2};

    std::vector<double> weights(hierarchy.total_nodes(), 0.0);
    std::vector<std::vector<double>> abscissas(2, std::vector<double>(hierarchy.total_nodes(), 0.0));
    for (std::size_t i = 0; i < 2; ++i)
    {
        for (std::size_t j = 0; j < 2; ++j)
        {
            const std::size_t flat = hierarchy.flatten_node_index({i, j});
            weights[flat] = x_weights[i] * y_weights[j];
            abscissas[0][flat] = x[i];
            abscissas[1][flat] = y[j];
        }
    }

    const std::vector<double> moments = raw_moments_from_grid(hierarchy, extents, weights, abscissas);
    const qmomkit::CqmomDenseCase case_{hierarchy, {extents, moments}, {}, {}, false, qmomkit::CqmomDiagnosticLevel::detailed};
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(case_);

    assert(result.ok());
    assert(result.status == qmomkit::CqmomStatus::success);
    assert(!result.diagnostics.empty());
    for (const qmomkit::CqmomDiagnostics& diagnostics : result.diagnostics)
    {
        assert(diagnostics.qmom_status == qmomkit::QmomStatus::success);
        assert(diagnostics.conditional_moment_status == qmomkit::ConditionalMomentStatus::success);
        if (diagnostics.conditioning_level + 1 < hierarchy.number_of_dimensions())
        {
            assert(diagnostics.conditional_moment_diagnostics.residual_norm.has_value());
            assert(diagnostics.conditional_moment_diagnostics.condition_estimate.has_value());
        }
    }

    for (const qmomkit::MultiIndex& node_index : qmomkit::enumerate_quadrature_node_multi_indices(hierarchy))
    {
        const std::size_t flat = hierarchy.flatten_node_index(node_index);
        assert_near(result.quadrature.weight(node_index), weights[flat], 1.0e-12);
        assert_near(result.quadrature.abscissa(0, node_index), abscissas[0][flat], 1.0e-12);
        assert_near(result.quadrature.abscissa(1, node_index), abscissas[1][flat], 1.0e-12);
    }

    const std::vector<double> first_node = result.quadrature.node({0, 0});
    const std::vector<double> last_node = result.quadrature.node({1, 1});
    assert(first_node.size() == 2);
    assert(last_node.size() == 2);
    assert_near(first_node[0], x[0], 1.0e-12);
    assert_near(first_node[1], y[0], 1.0e-12);
    assert_near(last_node[0], x[1], 1.0e-12);
    assert_near(last_node[1], y[1], 1.0e-12);
}

void test_cqmom_automatic_scaling_restores_each_physical_coordinate()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {1, 0});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> x{2.0e50, 5.0e50};
    const std::vector<double> y{-4.0e-50, 3.0e-50};
    const std::vector<double> x_weights{0.25, 0.75};
    const std::vector<double> y_weights{0.6, 0.4};
    std::vector<double> weights(hierarchy.total_nodes(), 0.0);
    std::vector<std::vector<double>> abscissas(2, std::vector<double>(hierarchy.total_nodes(), 0.0));
    for (std::size_t i = 0; i < 2; ++i)
    {
        for (std::size_t j = 0; j < 2; ++j)
        {
            const std::size_t flat = hierarchy.flatten_node_index({i, j});
            weights[flat] = x_weights[i] * y_weights[j];
            abscissas[0][flat] = x[i];
            abscissas[1][flat] = y[j];
        }
    }

    qmomkit::QmomOptions options;
    options.automatic_scaling = true;
    qmomkit::CqmomDenseCase case_{hierarchy, {extents, raw_moments_from_grid(hierarchy, extents, weights, abscissas)},
                                  {},        options,
                                  false,     qmomkit::CqmomDiagnosticLevel::detailed};
    case_.coordinate_support = {qmomkit::AbscissaBounds{1.0e50, 6.0e50}, qmomkit::AbscissaBounds{-5.0e-50, 4.0e-50}};
    qmomkit::CqmomLeanCase lean_case{qmomkit::extract_cqmom_moments(case_.moments, qmomkit::create_cqmom_moment_layout(hierarchy)),
                                     {},
                                     options,
                                     false,
                                     qmomkit::CqmomDiagnosticLevel::detailed};
    lean_case.coordinate_support = case_.coordinate_support;

    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(case_);
    const qmomkit::CqmomResult lean_result = qmomkit::compute_cqmom(lean_case);

    const std::vector<double> reference_x{2.0, 5.0};
    const std::vector<double> reference_y{-4.0, 3.0};
    std::vector<std::vector<double>> reference_abscissas(2, std::vector<double>(hierarchy.total_nodes(), 0.0));
    for (std::size_t i = 0; i < 2; ++i)
    {
        for (std::size_t j = 0; j < 2; ++j)
        {
            const std::size_t flat = hierarchy.flatten_node_index({i, j});
            reference_abscissas[0][flat] = reference_x[i];
            reference_abscissas[1][flat] = reference_y[j];
        }
    }
    const qmomkit::CqmomDenseCase reference_case{hierarchy, {extents, raw_moments_from_grid(hierarchy, extents, weights, reference_abscissas)},
                                                 {},        options,
                                                 false,     qmomkit::CqmomDiagnosticLevel::detailed};
    const qmomkit::CqmomResult reference = qmomkit::compute_cqmom(reference_case);

    assert(result.ok());
    assert(lean_result.ok());
    assert(reference.ok());
    for (const qmomkit::MultiIndex& node_index : qmomkit::enumerate_quadrature_node_multi_indices(hierarchy))
    {
        const std::size_t flat = hierarchy.flatten_node_index(node_index);
        assert_near(result.quadrature.weight(node_index), weights[flat], 1.0e-12);
        assert_relative_near(result.quadrature.abscissa(0, node_index), abscissas[0][flat], 1.0e-11);
        assert_relative_near(result.quadrature.abscissa(1, node_index), abscissas[1][flat], 1.0e-11);
        assert_near(lean_result.quadrature.weight(node_index), result.quadrature.weight(node_index), 1.0e-12);
        assert_relative_near(lean_result.quadrature.abscissa(0, node_index), result.quadrature.abscissa(0, node_index), 1.0e-12);
        assert_relative_near(lean_result.quadrature.abscissa(1, node_index), result.quadrature.abscissa(1, node_index), 1.0e-12);
        assert_near(result.quadrature.weight(node_index), reference.quadrature.weight(node_index), 1.0e-12);
        assert_relative_near(result.quadrature.abscissa(0, node_index) / 1.0e50, reference.quadrature.abscissa(0, node_index), 1.0e-11);
        assert_relative_near(result.quadrature.abscissa(1, node_index) / 1.0e-50, reference.quadrature.abscissa(1, node_index), 1.0e-11);
    }
    assert(!result.diagnostics.empty());
    for (const qmomkit::CqmomDiagnostics& diagnostics : result.diagnostics)
    {
        if (!diagnostics.bounded_support_interval_present)
        {
            continue;
        }
        const qmomkit::AbscissaBounds expected = case_.coordinate_support[diagnostics.phase_space_dimension].value();
        assert_relative_near(diagnostics.bounded_support_interval.lower, expected.lower, 1.0e-12);
        assert_relative_near(diagnostics.bounded_support_interval.upper, expected.upper, 1.0e-12);
    }
}

void test_cqmom_public_api_reports_conditional_moment_failure()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.5, 0.5};
    const std::vector<std::vector<double>> abscissas{{-1.0, 1.0}, {2.0, -1.0}};
    qmomkit::ConditionalMomentOptions conditional_options;
    conditional_options.singular_tolerance = 1.0;
    const qmomkit::CqmomDenseCase case_{hierarchy, {extents, raw_moments_from_terminals(extents, weights, abscissas)}, conditional_options, {}};

    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(case_);

    assert(!result.ok());
    assert(result.status == qmomkit::CqmomStatus::conditional_moment_failed);
    assert(result.failed_phase_space_dimension == 0);
    // The public contract makes the quadrature diagnostic-only whenever the
    // overall CQMOM result fails; callers must gate it on result.ok().
    assert(!result.message.empty());
    assert(result.message.find("Vandermonde") != std::string::npos);
    assert(!result.diagnostics.empty());
    const auto failed_diagnostic =
        std::find_if(result.diagnostics.begin(), result.diagnostics.end(), [](const qmomkit::CqmomDiagnostics& diagnostics)
                     { return diagnostics.conditional_moment_status == qmomkit::ConditionalMomentStatus::singular_system; });
    assert(failed_diagnostic != result.diagnostics.end());
    assert(failed_diagnostic->conditioning_level == 0);
    assert(failed_diagnostic->phase_space_dimension == 0);
    assert(failed_diagnostic->parent_path.empty());
    assert(!failed_diagnostic->message.empty());
}

void test_cqmom_public_api_exposes_adaptive_branch_counts()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.1, 0.3, 0.6};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 1.0}, {-2.0, 2.0, 3.0}};
    const std::vector<double> moments = raw_moments_from_terminals(extents, weights, abscissas);

    const qmomkit::CqmomDenseCase case_{
        hierarchy, {extents, moments}, {}, {true, 1.0e-10, 1.0e-10, 1.0e-12}, true, qmomkit::CqmomDiagnosticLevel::basic};
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(case_);

    assert(result.ok());
    assert(result.quadrature.adaptive_branch_counts_enabled());
    assert(!result.quadrature.has_rectangular_storage());
    assert(result.quadrature.total_nodes() == 3);
    assert(result.quadrature.terminal_paths().size() == 3);

    bool saw_reduced_branch = false;
    for (const qmomkit::CqmomDiagnostics& diagnostics : result.diagnostics)
    {
        assert(diagnostics.adaptive_branch_counts_enabled);
        if (diagnostics.conditioning_level == 1 && diagnostics.parent_path == qmomkit::ConditioningPath({1}))
        {
            saw_reduced_branch = true;
            assert(diagnostics.requested_node_count == 2);
            assert(diagnostics.actual_node_count == 1);
            assert(diagnostics.adaptive_node_count_reduced);
        }
    }
    assert(saw_reduced_branch);
}

void test_cqmom_adaptive_root_branch_uses_relative_rmin_only()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const double density_scale = 0.05;
    const std::vector<double> weights(4, 0.25 * density_scale);
    const std::vector<std::vector<double>> abscissas{{10.0, 10.0, 20.0, 20.0}, {-1.0, 1.0, -1.0, 1.0}};
    const qmomkit::QmomOptions qmom_options{true, 0.1, 1.0e-10, 1.0e-30};
    const qmomkit::CqmomDenseCase case_{
        hierarchy, {extents, raw_moments_from_terminals(extents, weights, abscissas)}, {}, qmom_options, true, qmomkit::CqmomDiagnosticLevel::basic};

    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(case_);

    assert(result.ok());
    assert(result.quadrature.total_nodes() == 4);
    assert_near(std::accumulate(result.quadrature.weights().begin(), result.quadrature.weights().end(), 0.0), density_scale, 1.0e-14);

    bool saw_unreduced_root = false;
    for (const qmomkit::CqmomDiagnostics& diagnostics : result.diagnostics)
    {
        if (diagnostics.conditioning_level == 0)
        {
            saw_unreduced_root = true;
            assert(diagnostics.actual_node_count == 2);
            assert(!diagnostics.adaptive_node_count_reduced);
        }
    }
    assert(saw_unreduced_root);
}

void test_cqmom_public_api_rejects_adaptive_wheeler_without_adaptive_branch_counts()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.1, 0.2, 0.3, 0.4};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 1.0, 1.0}, {-2.0, 2.0, -2.0, 2.0}};
    const std::vector<double> moments = raw_moments_from_terminals(extents, weights, abscissas);
    const qmomkit::QmomOptions adaptive_qmom{true, 1.0e-10, 1.0e-10, 1.0e-12};

    const qmomkit::CqmomResult fixed = qmomkit::compute_cqmom({hierarchy, {extents, moments}, {}, {}, false});
    assert(fixed.ok());
    assert(!fixed.quadrature.adaptive_branch_counts_enabled());
    assert(fixed.quadrature.has_rectangular_storage());

    const qmomkit::CqmomResult rejected = qmomkit::compute_cqmom({hierarchy, {extents, moments}, {}, adaptive_qmom, false});
    assert(!rejected.ok());
    assert(rejected.status == qmomkit::CqmomStatus::invalid_input);
    assert(rejected.failed_phase_space_dimension == 0);
    assert(rejected.diagnostics.empty());
    assert(rejected.message.find("Adaptive branch-specific node counts") != std::string::npos);
    assert(rejected.message.find("adaptive_branch_counts") != std::string::npos);

    const qmomkit::CqmomResult adaptive = qmomkit::compute_cqmom({hierarchy, {extents, moments}, {}, adaptive_qmom, true});
    assert(adaptive.ok());
    assert(adaptive.quadrature.adaptive_branch_counts_enabled());
    assert(!adaptive.quadrature.has_rectangular_storage());
    assert(adaptive.quadrature.total_nodes() == weights.size());
}

void test_cqmom_diagnostic_levels_control_expensive_diagnostics()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.2, 0.3, 0.1, 0.4};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 1.5, 1.5}, {0.25, 2.0, 0.25, 2.0}};
    const std::vector<double> moments = raw_moments_from_grid(hierarchy, extents, weights, abscissas);

    qmomkit::CqmomDenseCase default_case{hierarchy, {extents, moments}, {}, {}};
    assert(default_case.diagnostic_level == qmomkit::CqmomDiagnosticLevel::none);

    qmomkit::internal::performance::CounterSession default_session;
    const qmomkit::CqmomResult default_result = [&]()
    {
        qmomkit::internal::performance::CounterScope scope(default_session);
        return qmomkit::compute_cqmom(default_case);
    }();
    assert(default_result.ok());
    assert(default_result.diagnostics.empty());
    assert(default_session.snapshot().linear_algebra.condition_estimates == 0);

    qmomkit::CqmomDenseCase detailed_case = default_case;
    detailed_case.diagnostic_level = qmomkit::CqmomDiagnosticLevel::detailed;
    const qmomkit::CqmomResult detailed_result = qmomkit::compute_cqmom(detailed_case);
    assert(detailed_result.ok());
    assert_same_quadrature(default_result.quadrature, detailed_result.quadrature, 1.0e-13);
    assert(!detailed_result.diagnostics.empty());
    assert(detailed_result.diagnostics.front().conditional_moment_diagnostics.residual_norm.has_value());
    assert(detailed_result.diagnostics.front().conditional_moment_diagnostics.condition_estimate.has_value());

    qmomkit::CqmomDenseCase basic_case = default_case;
    basic_case.diagnostic_level = qmomkit::CqmomDiagnosticLevel::basic;
    const qmomkit::CqmomResult basic_result = qmomkit::compute_cqmom(basic_case);
    assert(basic_result.ok());
    assert_same_quadrature(default_result.quadrature, basic_result.quadrature, 1.0e-13);
    assert(basic_result.diagnostics.size() == detailed_result.diagnostics.size());
    for (const qmomkit::CqmomDiagnostics& diagnostics : basic_result.diagnostics)
    {
        assert(diagnostics.qmom_status == qmomkit::QmomStatus::success);
        assert(diagnostics.conditional_moment_status == qmomkit::ConditionalMomentStatus::success);
        assert(!diagnostics.conditional_moment_diagnostics.residual_norm.has_value());
        assert(!diagnostics.conditional_moment_diagnostics.condition_estimate.has_value());
        assert(!diagnostics.conditional_moment_diagnostics.conditioning_warning);
    }

    qmomkit::CqmomDenseCase none_case = default_case;
    none_case.diagnostic_level = qmomkit::CqmomDiagnosticLevel::none;
    const qmomkit::CqmomResult none_result = qmomkit::compute_cqmom(none_case);
    assert(none_result.ok());
    assert_same_quadrature(default_result.quadrature, none_result.quadrature, 1.0e-13);
    assert(none_result.diagnostics.empty());
}

void test_cqmom_none_diagnostics_preserves_failure_context()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    qmomkit::CqmomDenseCase dense_case{hierarchy, {{4, 4}, std::vector<double>(16, 0.0)}, {}, {}};
    dense_case.diagnostic_level = qmomkit::CqmomDiagnosticLevel::none;

    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(dense_case);

    assert(!result.ok());
    assert(result.status == qmomkit::CqmomStatus::qmom_failed);
    assert(!result.message.empty());
    assert(!result.diagnostics.empty());
    assert(result.diagnostics.front().qmom_status != qmomkit::QmomStatus::success);
    assert(!result.diagnostics.front().message.empty());
}

void test_cqmom_diagnostics_default_bounded_support_placeholder()
{
    const qmomkit::CqmomDiagnostics diagnostics;
    assert(!diagnostics.bounded_support_checked);
    assert(!diagnostics.bounded_support_interval_present);
    assert_near(diagnostics.bounded_support_interval.lower, 0.0, 0.0);
    assert_near(diagnostics.bounded_support_interval.upper, 0.0, 0.0);
}

void test_cqmom_storage_settings_format()
{
    assert(qmomkit::to_string(qmomkit::CqmomMomentStorage::lean) == "lean");
    assert(qmomkit::to_string(qmomkit::CqmomMomentStorage::dense) == "dense");

    const std::string formatted = qmomkit::format_cqmom_settings(qmomkit::CqmomMomentStorage::dense, true, qmomkit::VandermondeSolver::dense,
                                                                 qmomkit::CqmomDiagnosticLevel::detailed);
    assert(formatted.find("moment storage = dense") != std::string::npos);
    assert(formatted.find("adaptive branch counts = enabled") != std::string::npos);
    assert(formatted.find("Vandermonde solver = dense") != std::string::npos);
    assert(formatted.find("diagnostic level = detailed") != std::string::npos);

    const std::string pivoted_qr = qmomkit::format_cqmom_settings(qmomkit::CqmomMomentStorage::lean, false, qmomkit::VandermondeSolver::pivoted_qr,
                                                                  qmomkit::CqmomDiagnosticLevel::basic);
    assert(pivoted_qr.find("Vandermonde solver = pivoted_qr") != std::string::npos);

    const std::string svd = qmomkit::format_cqmom_settings(qmomkit::CqmomMomentStorage::lean, false, qmomkit::VandermondeSolver::svd,
                                                           qmomkit::CqmomDiagnosticLevel::basic);
    assert(svd.find("Vandermonde solver = svd") != std::string::npos);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"moment_set_public_container_behavior", test_moment_set_public_container_behavior},
            {"cqmom_moment_layout_two_dimensional_reduced_set", test_cqmom_moment_layout_two_dimensional_reduced_set},
            {"cqmom_moment_layout_three_dimensional_recursive_indices", test_cqmom_moment_layout_three_dimensional_recursive_indices},
            {"cqmom_moment_layout_non_identity_conditioning_order", test_cqmom_moment_layout_non_identity_conditioning_order},
            {"cqmom_moment_layout_unequal_node_counts", test_cqmom_moment_layout_unequal_node_counts},
            {"select_cqmom_moments_extracts_required_dense_values", test_select_cqmom_moments_extracts_required_dense_values},
            {"cqmom_moment_data_lookup_reflects_public_value_updates", test_cqmom_moment_data_lookup_reflects_public_value_updates},
            {"select_cqmom_moments_respects_non_identity_conditioning_order", test_select_cqmom_moments_respects_non_identity_conditioning_order},
            {"select_cqmom_moments_handles_unequal_node_counts", test_select_cqmom_moments_handles_unequal_node_counts},
            {"select_cqmom_moments_reports_missing_dense_moment", test_select_cqmom_moments_reports_missing_dense_moment},
            {"create_cqmom_moment_data_converts_valid_selected_result", test_create_cqmom_moment_data_converts_valid_selected_result},
            {"create_cqmom_moment_data_rejects_failed_selected_result", test_create_cqmom_moment_data_rejects_failed_selected_result},
            {"create_cqmom_moment_data_rejects_index_order_mismatch", test_create_cqmom_moment_data_rejects_index_order_mismatch},
            {"create_cqmom_moment_data_rejects_missing_values", test_create_cqmom_moment_data_rejects_missing_values},
            {"create_cqmom_moment_data_rejects_extra_values", test_create_cqmom_moment_data_rejects_extra_values},
            {"create_cqmom_moment_data_rejects_non_finite_values", test_create_cqmom_moment_data_rejects_non_finite_values},
            {"create_cqmom_moment_data_from_independent_product_selected_moments",
             test_create_cqmom_moment_data_from_independent_product_selected_moments},
            {"lean_cqmom_public_api_matches_dense_two_dimensional", test_lean_cqmom_public_api_matches_dense_two_dimensional},
            {"lean_cqmom_public_api_matches_dense_three_dimensional", test_lean_cqmom_public_api_matches_dense_three_dimensional},
            {"lean_cqmom_public_api_matches_dense_non_identity", test_lean_cqmom_public_api_matches_dense_non_identity},
            {"lean_cqmom_public_api_matches_dense_unequal_node_counts", test_lean_cqmom_public_api_matches_dense_unequal_node_counts},
            {"lean_cqmom_public_api_matches_dense_adaptive", test_lean_cqmom_public_api_matches_dense_adaptive},
            {"lean_cqmom_public_api_rejects_missing_values", test_lean_cqmom_public_api_rejects_missing_values},
            {"lean_cqmom_public_api_rejects_non_finite_values", test_lean_cqmom_public_api_rejects_non_finite_values},
            {"lean_cqmom_public_api_rejects_invalid_options", test_lean_cqmom_public_api_rejects_invalid_options},
            {"lean_cqmom_public_api_rejects_incompatible_layout", test_lean_cqmom_public_api_rejects_incompatible_layout},
            {"cqmom_public_api_reports_missing_required_dense_moment", test_cqmom_public_api_reports_missing_required_dense_moment},
            {"cqmom_public_api_preserves_dense_invalid_input_messages", test_cqmom_public_api_preserves_dense_invalid_input_messages},
            {"cqmom_public_api_rejects_unrepresentable_required_moment_extent_before_storage",
             test_cqmom_public_api_rejects_unrepresentable_required_moment_extent_before_storage},
            {"qmom_public_api_matches_wheeler", test_qmom_public_api_matches_wheeler},
            {"qmom_rejects_node_count_with_unrepresentable_moment_requirement", test_qmom_rejects_node_count_with_unrepresentable_moment_requirement},
            {"qmom_automatic_scaling_restores_physical_abscissas_and_bounds", test_qmom_automatic_scaling_restores_physical_abscissas_and_bounds},
            {"qmom_automatic_scaling_is_covariant_and_uses_existing_scale_moments",
             test_qmom_automatic_scaling_is_covariant_and_uses_existing_scale_moments},
            {"qmom_automatic_scaling_rejects_nonzero_underflow_before_reconstruction",
             test_qmom_automatic_scaling_rejects_nonzero_underflow_before_reconstruction},
            {"cqmom_automatic_scaling_falls_back_atomically_when_preparation_is_unsafe",
             test_cqmom_automatic_scaling_falls_back_atomically_when_preparation_is_unsafe},
            {"one_dimensional_cqmom_agrees_with_qmom", test_one_dimensional_cqmom_agrees_with_qmom},
            {"cqmom_public_api_two_dimensional_manufactured_case", test_cqmom_public_api_two_dimensional_manufactured_case},
            {"cqmom_automatic_scaling_restores_each_physical_coordinate", test_cqmom_automatic_scaling_restores_each_physical_coordinate},
            {"cqmom_public_api_reports_conditional_moment_failure", test_cqmom_public_api_reports_conditional_moment_failure},
            {"cqmom_public_api_exposes_adaptive_branch_counts", test_cqmom_public_api_exposes_adaptive_branch_counts},
            {"cqmom_adaptive_root_branch_uses_relative_rmin_only", test_cqmom_adaptive_root_branch_uses_relative_rmin_only},
            {"cqmom_public_api_rejects_adaptive_wheeler_without_adaptive_branch_counts",
             test_cqmom_public_api_rejects_adaptive_wheeler_without_adaptive_branch_counts},
            {"cqmom_diagnostic_levels_control_expensive_diagnostics", test_cqmom_diagnostic_levels_control_expensive_diagnostics},
            {"cqmom_none_diagnostics_preserves_failure_context", test_cqmom_none_diagnostics_preserves_failure_context},
            {"cqmom_diagnostics_default_bounded_support_placeholder", test_cqmom_diagnostics_default_bounded_support_placeholder},
            {"cqmom_storage_settings_format", test_cqmom_storage_settings_format},
        });
}
