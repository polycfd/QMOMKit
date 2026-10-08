#include "internal/recursive_quadrature.hpp"

#include "internal/coordinate_scales.hpp"
#include "qmomkit/qmomkit.hpp"
#include "test_case_dispatch.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{

// These tests construct a full tensor of raw moments; reconstruction takes only
// the moments the hierarchy requires. The request is fixed unless the Wheeler
// options ask for adaptive reconstruction.
qmomkit::internal::RecursiveQuadratureMomentDataRequest
moment_data_request(const qmomkit::PhaseSpaceHierarchy& hierarchy, const qmomkit::RawMomentTensor& moment_tensor,
                    const qmomkit::ConditionalMomentOptions& conditional_options = {},
                    const qmomkit::WheelerOptions& wheeler_options = qmomkit::WheelerOptions{false},
                    const qmomkit::CqmomDiagnosticLevel diagnostic_level = qmomkit::CqmomDiagnosticLevel::none,
                    const qmomkit::VandermondeSolver solver = qmomkit::VandermondeSolver::automatic,
                    const qmomkit::CqmomAbscissaBounds& coordinate_support = {})
{
    return {hierarchy,
            qmomkit::extract_cqmom_moments(moment_tensor, qmomkit::create_cqmom_moment_layout(hierarchy)),
            conditional_options,
            wheeler_options,
            diagnostic_level,
            solver,
            {},
            false,
            coordinate_support};
}

void assert_vector_equal(const std::vector<std::size_t>& actual, const std::vector<std::size_t>& expected)
{
    assert(actual == expected);
}

void assert_vector_near(const std::vector<double>& actual, const std::vector<double>& expected)
{
    assert(actual.size() == expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        assert(actual[i] == expected[i]);
    }
}

void assert_near(const double actual, const double expected, const double tolerance)
{
    assert(std::abs(actual - expected) <= tolerance);
}

std::vector<double> raw_moments_from_quadrature(const qmomkit::internal::RecursiveQuadrature& quadrature, const std::vector<std::size_t>& extents)
{
    std::vector<double> moments(qmomkit::product_extents(extents), 0.0);
    const qmomkit::PhaseSpaceHierarchy& hierarchy = quadrature.hierarchy();
    for (std::size_t flat_moment = 0; flat_moment < moments.size(); ++flat_moment)
    {
        const qmomkit::MultiIndex moment_index = qmomkit::unflatten_multi_index(flat_moment, extents);
        double value = 0.0;
        for (const qmomkit::MultiIndex& node_index : qmomkit::enumerate_quadrature_node_multi_indices(hierarchy))
        {
            double contribution = quadrature.weight(node_index);
            for (std::size_t dimension = 0; dimension < hierarchy.number_of_dimensions(); ++dimension)
            {
                contribution *= std::pow(quadrature.abscissa(dimension, node_index), static_cast<int>(moment_index[dimension]));
            }
            value += contribution;
        }
        moments[flat_moment] = value;
    }
    return moments;
}

std::vector<double> raw_moments_from_terminal_quadrature(const std::vector<std::vector<double>>& abscissas, const std::vector<double>& weights,
                                                         const std::vector<std::size_t>& extents)
{
    std::vector<double> moments(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat_moment = 0; flat_moment < moments.size(); ++flat_moment)
    {
        const qmomkit::MultiIndex moment_index = qmomkit::unflatten_multi_index(flat_moment, extents);
        double value = 0.0;
        for (std::size_t node = 0; node < weights.size(); ++node)
        {
            double contribution = weights[node];
            for (std::size_t dimension = 0; dimension < abscissas.size(); ++dimension)
            {
                contribution *= std::pow(abscissas[dimension][node], static_cast<int>(moment_index[dimension]));
            }
            value += contribution;
        }
        moments[flat_moment] = value;
    }
    return moments;
}

std::vector<double> indexed_values(const std::vector<std::size_t>& extents)
{
    std::vector<double> values(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat = 0; flat < values.size(); ++flat)
    {
        values[flat] = 100.0 + static_cast<double>(flat);
    }
    return values;
}

struct BranchMomentSliceFixture
{
    std::shared_ptr<const qmomkit::internal::BranchMomentLayout> layout;
    qmomkit::internal::BranchMomentSlice slice;
};

BranchMomentSliceFixture make_branch_slice_from_tensor(const qmomkit::PhaseSpaceHierarchy& hierarchy, const qmomkit::RawMomentTensor& moment_tensor)
{
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const std::vector<std::size_t> moment_extents = hierarchy.moment_extents();
    std::vector<std::size_t> ordered_extents;
    ordered_extents.reserve(hierarchy.number_of_dimensions());
    for (const std::size_t dimension : hierarchy.conditioning_order())
    {
        ordered_extents.push_back(moment_extents[dimension]);
    }

    std::vector<qmomkit::MultiIndex> ordered_indices;
    std::vector<double> values;
    ordered_indices.reserve(layout.required_indices.size());
    values.reserve(layout.required_indices.size());
    for (const qmomkit::MultiIndexView phase_space_index : layout.required_indices)
    {
        qmomkit::MultiIndex ordered_index(hierarchy.number_of_dimensions(), 0);
        for (std::size_t axis = 0; axis < hierarchy.number_of_dimensions(); ++axis)
        {
            ordered_index[axis] = phase_space_index[hierarchy.conditioning_order()[axis]];
        }
        ordered_indices.push_back(std::move(ordered_index));
        values.push_back(moment_tensor.values[qmomkit::flatten_multi_index(phase_space_index, moment_tensor.extents)]);
    }

    std::shared_ptr<const qmomkit::internal::BranchMomentLayout> branch_layout =
        qmomkit::internal::create_branch_moment_layout(std::move(ordered_extents), std::move(ordered_indices));
    return {branch_layout, {branch_layout.get(), std::move(values)}};
}

void assert_matrix_value(const qmomkit::internal::DenseMatrix& matrix, const std::size_t row, const std::size_t col, const double expected)
{
    assert(matrix(row, col) == expected);
}

void assert_raw_moments_reconstructed(const qmomkit::internal::RecursiveQuadratureResult& result, const std::vector<double>& expected_moments,
                                      const std::vector<std::size_t>& extents, const double tolerance)
{
    assert(result.ok());
    const std::vector<double> reconstructed = raw_moments_from_quadrature(result.quadrature, extents);
    assert(reconstructed.size() == expected_moments.size());
    for (std::size_t i = 0; i < expected_moments.size(); ++i)
    {
        assert_near(reconstructed[i], expected_moments[i], tolerance);
    }
}

void assert_quadrature_matches(const qmomkit::internal::RecursiveQuadrature& actual, const qmomkit::internal::RecursiveQuadrature& expected,
                               const double tolerance)
{
    assert(actual.hierarchy().node_counts() == expected.hierarchy().node_counts());
    assert(actual.hierarchy().conditioning_order() == expected.hierarchy().conditioning_order());
    for (const qmomkit::MultiIndex& node_index : qmomkit::enumerate_quadrature_node_multi_indices(expected.hierarchy()))
    {
        assert_near(actual.weight(node_index), expected.weight(node_index), tolerance);
        for (std::size_t dimension = 0; dimension < expected.number_of_dimensions(); ++dimension)
        {
            assert_near(actual.abscissa(dimension, node_index), expected.abscissa(dimension, node_index), tolerance);
        }
    }
}

void assert_terminal_quadrature_matches(const qmomkit::internal::RecursiveQuadrature& actual, const qmomkit::internal::RecursiveQuadrature& expected,
                                        const double tolerance)
{
    assert(actual.hierarchy().node_counts() == expected.hierarchy().node_counts());
    assert(actual.hierarchy().conditioning_order() == expected.hierarchy().conditioning_order());
    assert(actual.weights().size() == expected.weights().size());
    assert(actual.abscissas().size() == expected.abscissas().size());
    assert(actual.terminal_paths() == expected.terminal_paths());
    for (std::size_t node = 0; node < expected.weights().size(); ++node)
    {
        assert_near(actual.weights()[node], expected.weights()[node], tolerance);
    }
    for (std::size_t dimension = 0; dimension < expected.abscissas().size(); ++dimension)
    {
        assert(actual.abscissas()[dimension].size() == expected.abscissas()[dimension].size());
        for (std::size_t node = 0; node < expected.abscissas()[dimension].size(); ++node)
        {
            assert_near(actual.abscissas()[dimension][node], expected.abscissas()[dimension][node], tolerance);
        }
    }
}

void assert_recursive_results_equivalent(const qmomkit::internal::RecursiveQuadratureResult& actual,
                                         const qmomkit::internal::RecursiveQuadratureResult& expected, const double tolerance)
{
    assert(actual.status == expected.status);
    assert(actual.ok() == expected.ok());
    assert(actual.failed_path == expected.failed_path);
    assert(actual.failed_conditioning_level == expected.failed_conditioning_level);
    assert(actual.diagnostics.size() == expected.diagnostics.size());
    for (std::size_t i = 0; i < expected.diagnostics.size(); ++i)
    {
        assert(actual.diagnostics[i].conditioning_level == expected.diagnostics[i].conditioning_level);
        assert(actual.diagnostics[i].path == expected.diagnostics[i].path);
        assert(actual.diagnostics[i].phase_space_dimension == expected.diagnostics[i].phase_space_dimension);
        assert(actual.diagnostics[i].requested_node_count == expected.diagnostics[i].requested_node_count);
        assert(actual.diagnostics[i].actual_node_count == expected.diagnostics[i].actual_node_count);
        assert(actual.diagnostics[i].wheeler_status == expected.diagnostics[i].wheeler_status);
        assert(actual.diagnostics[i].conditional_moment_status == expected.diagnostics[i].conditional_moment_status);
    }

    if (actual.ok())
    {
        if (expected.quadrature.has_rectangular_storage())
        {
            assert_quadrature_matches(actual.quadrature, expected.quadrature, tolerance);
        }
        else
        {
            assert_terminal_quadrature_matches(actual.quadrature, expected.quadrature, tolerance);
        }
    }
}

void assert_recursive_numerics_equivalent(const qmomkit::internal::RecursiveQuadratureResult& actual,
                                          const qmomkit::internal::RecursiveQuadratureResult& expected, const double tolerance)
{
    assert(actual.status == expected.status);
    assert(actual.ok() == expected.ok());
    assert(actual.failed_path == expected.failed_path);
    assert(actual.failed_conditioning_level == expected.failed_conditioning_level);

    if (actual.ok())
    {
        if (expected.quadrature.has_rectangular_storage())
        {
            assert_quadrature_matches(actual.quadrature, expected.quadrature, tolerance);
        }
        else
        {
            assert_terminal_quadrature_matches(actual.quadrature, expected.quadrature, tolerance);
        }
    }
}

std::vector<qmomkit::CqmomDiagnosticLevel> diagnostic_levels()
{
    return {qmomkit::CqmomDiagnosticLevel::none, qmomkit::CqmomDiagnosticLevel::basic, qmomkit::CqmomDiagnosticLevel::detailed};
}

qmomkit::CqmomAbscissaBounds bounded_support_for_dimension(const std::size_t dimensions, const std::size_t dimension,
                                                           const qmomkit::AbscissaBounds interval)
{
    qmomkit::CqmomAbscissaBounds support(dimensions);
    support.at(dimension) = interval;
    return support;
}

void assert_cqmom_results_numerically_equivalent(const qmomkit::CqmomResult& actual, const qmomkit::CqmomResult& expected, const double tolerance)
{
    assert(actual.status == expected.status);
    assert(actual.ok() == expected.ok());
    assert(actual.failed_phase_space_dimension == expected.failed_phase_space_dimension);
    if (!actual.ok())
    {
        return;
    }

    assert(actual.quadrature.hierarchy().node_counts() == expected.quadrature.hierarchy().node_counts());
    assert(actual.quadrature.hierarchy().conditioning_order() == expected.quadrature.hierarchy().conditioning_order());
    assert(actual.quadrature.adaptive_branch_counts_enabled() == expected.quadrature.adaptive_branch_counts_enabled());
    assert(actual.quadrature.has_rectangular_storage() == expected.quadrature.has_rectangular_storage());
    assert(actual.quadrature.terminal_paths() == expected.quadrature.terminal_paths());
    assert(actual.quadrature.weights().size() == expected.quadrature.weights().size());
    assert(actual.quadrature.abscissas().size() == expected.quadrature.abscissas().size());
    for (std::size_t node = 0; node < expected.quadrature.weights().size(); ++node)
    {
        assert_near(actual.quadrature.weights()[node], expected.quadrature.weights()[node], tolerance);
    }
    for (std::size_t dimension = 0; dimension < expected.quadrature.abscissas().size(); ++dimension)
    {
        assert(actual.quadrature.abscissas()[dimension].size() == expected.quadrature.abscissas()[dimension].size());
        for (std::size_t node = 0; node < expected.quadrature.abscissas()[dimension].size(); ++node)
        {
            assert_near(actual.quadrature.abscissas()[dimension][node], expected.quadrature.abscissas()[dimension][node], tolerance);
        }
    }
}

void assert_terminal_coordinate_inside_interval(const qmomkit::internal::RecursiveQuadratureResult& result, const std::size_t dimension,
                                                const qmomkit::AbscissaBounds interval, const double tolerance)
{
    assert(result.ok());
    assert(dimension < result.quadrature.abscissas().size());
    for (const double node : result.quadrature.abscissas()[dimension])
    {
        assert(node >= interval.lower - tolerance);
        assert(node <= interval.upper + tolerance);
    }
}

void assert_terminal_weights_positive(const qmomkit::internal::RecursiveQuadratureResult& result)
{
    assert(result.ok());
    for (const double weight : result.quadrature.weights())
    {
        assert(weight > 0.0);
    }
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

template <typename Function> void assert_throws_logic_error(Function function)
{
    bool threw = false;
    try
    {
        function();
    }
    catch (const std::logic_error&)
    {
        threw = true;
    }
    assert(threw);
}

void assert_ragged_tree_consistent(const qmomkit::internal::RecursiveQuadrature& quadrature)
{
    assert(quadrature.adaptive_branch_counts_enabled());
    assert(!quadrature.has_rectangular_storage());
    assert(!quadrature.tree_nodes().empty());
    assert(quadrature.terminal_node_ids().size() == quadrature.total_nodes());
    assert(quadrature.terminal_paths().size() == quadrature.total_nodes());

    const qmomkit::internal::RecursiveQuadrature::TreeNode& root = quadrature.tree_nodes().front();
    assert(root.parent_id == std::numeric_limits<std::size_t>::max());
    assert(root.path.empty());
    assert(!root.terminal);

    for (std::size_t node_id = 0; node_id < quadrature.tree_nodes().size(); ++node_id)
    {
        const qmomkit::internal::RecursiveQuadrature::TreeNode& node = quadrature.tree_nodes()[node_id];
        assert(node.id == node_id);
        assert(node.abscissas_by_dimension.size() == quadrature.number_of_dimensions());
        if (node_id == 0)
        {
            continue;
        }

        assert(node.parent_id < quadrature.tree_nodes().size());
        const qmomkit::internal::RecursiveQuadrature::TreeNode& parent = quadrature.tree_nodes()[node.parent_id];
        assert(node.path.size() == parent.path.size() + 1);
        for (std::size_t level = 0; level < parent.path.size(); ++level)
        {
            assert(node.path.node_index_at_level(level) == parent.path.node_index_at_level(level));
        }
        assert(node.path.node_index_at_level(node.path.size() - 1) == node.local_node_index);
        assert_near(node.accumulated_weight, parent.accumulated_weight * node.local_weight, 1.0e-12);

        bool listed_by_parent = false;
        for (const std::size_t child_id : parent.children)
        {
            listed_by_parent = listed_by_parent || child_id == node_id;
        }
        assert(listed_by_parent);
    }

    for (std::size_t terminal = 0; terminal < quadrature.total_nodes(); ++terminal)
    {
        const std::size_t node_id = quadrature.terminal_node_ids()[terminal];
        assert(node_id < quadrature.tree_nodes().size());
        const qmomkit::internal::RecursiveQuadrature::TreeNode& node = quadrature.tree_nodes()[node_id];
        assert(node.terminal);
        assert(node.path == quadrature.terminal_paths()[terminal]);
        assert_near(node.accumulated_weight, quadrature.weights()[terminal], 1.0e-12);
        for (std::size_t dimension = 0; dimension < quadrature.number_of_dimensions(); ++dimension)
        {
            assert_near(node.abscissas_by_dimension[dimension], quadrature.abscissas()[dimension][terminal], 1.0e-12);
        }
    }
}

void test_flat_multi_index_round_trip()
{
    const std::vector<std::size_t> extents{2, 3, 4};
    const std::size_t total = qmomkit::product_extents(extents);

    assert(total == 24);
    for (std::size_t flat_index = 0; flat_index < total; ++flat_index)
    {
        const qmomkit::MultiIndex multi_index = qmomkit::unflatten_multi_index(flat_index, extents);
        assert(qmomkit::flatten_multi_index(multi_index, extents) == flat_index);
    }

    assert_vector_equal(qmomkit::unflatten_multi_index(0, extents), {0, 0, 0});
    assert_vector_equal(qmomkit::unflatten_multi_index(1, extents), {0, 0, 1});
    assert_vector_equal(qmomkit::unflatten_multi_index(4, extents), {0, 1, 0});
    assert_vector_equal(qmomkit::unflatten_multi_index(23, extents), {1, 2, 3});
    assert(qmomkit::flatten_multi_index({1, 2, 3}, extents) == 23);
}

void test_enumerate_multi_indices_row_major()
{
    const std::vector<qmomkit::MultiIndex> indices = qmomkit::enumerate_multi_indices({2, 3});

    assert(indices.size() == 6);
    assert_vector_equal(indices[0], {0, 0});
    assert_vector_equal(indices[1], {0, 1});
    assert_vector_equal(indices[2], {0, 2});
    assert_vector_equal(indices[3], {1, 0});
    assert_vector_equal(indices[5], {1, 2});
}

void test_branch_local_layout_identity_wheeler_and_rhs_offsets()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moment_tensor{{4, 4}, indexed_values({4, 4})};
    const BranchMomentSliceFixture fixture = make_branch_slice_from_tensor(hierarchy, moment_tensor);
    const qmomkit::internal::BranchMomentSlice& slice = fixture.slice;

    assert(slice.layout);
    assert_vector_equal(slice.layout->extents, {4, 4});
    assert(slice.values.size() == slice.layout->indices.size());
    assert_vector_near(qmomkit::internal::extract_wheeler_moments(slice), {100.0, 104.0, 108.0, 112.0});

    const qmomkit::internal::DenseMatrix rhs = qmomkit::internal::make_branch_rhs_matrix(slice, 2);
    assert(rhs.rows() == 2);
    assert(rhs.cols() == 4);
    assert_matrix_value(rhs, 0, 0, 100.0);
    assert_matrix_value(rhs, 0, 1, 101.0);
    assert_matrix_value(rhs, 0, 2, 102.0);
    assert_matrix_value(rhs, 0, 3, 103.0);
    assert_matrix_value(rhs, 1, 0, 104.0);
    assert_matrix_value(rhs, 1, 1, 105.0);
    assert_matrix_value(rhs, 1, 2, 106.0);
    assert_matrix_value(rhs, 1, 3, 107.0);
}

void test_branch_local_layout_non_identity_wheeler_and_rhs_offsets()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {1, 0});
    const qmomkit::RawMomentTensor moment_tensor{{4, 4}, indexed_values({4, 4})};
    const BranchMomentSliceFixture fixture = make_branch_slice_from_tensor(hierarchy, moment_tensor);
    const qmomkit::internal::BranchMomentSlice& slice = fixture.slice;

    assert(slice.layout);
    assert_vector_equal(slice.layout->extents, {4, 4});
    assert_vector_near(qmomkit::internal::extract_wheeler_moments(slice), {100.0, 101.0, 102.0, 103.0});

    const qmomkit::internal::DenseMatrix rhs = qmomkit::internal::make_branch_rhs_matrix(slice, 2);
    assert(rhs.rows() == 2);
    assert(rhs.cols() == 4);
    assert_matrix_value(rhs, 0, 0, 100.0);
    assert_matrix_value(rhs, 0, 1, 104.0);
    assert_matrix_value(rhs, 0, 2, 108.0);
    assert_matrix_value(rhs, 0, 3, 112.0);
    assert_matrix_value(rhs, 1, 0, 101.0);
    assert_matrix_value(rhs, 1, 1, 105.0);
    assert_matrix_value(rhs, 1, 2, 109.0);
    assert_matrix_value(rhs, 1, 3, 113.0);
}

void test_branch_local_layout_unequal_node_counts_and_adaptive_child_case()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 3});
    const qmomkit::RawMomentTensor moment_tensor{{4, 6}, indexed_values({4, 6})};
    const BranchMomentSliceFixture fixture = make_branch_slice_from_tensor(hierarchy, moment_tensor);
    const qmomkit::internal::BranchMomentSlice& slice = fixture.slice;

    assert(slice.layout);
    assert_vector_equal(slice.layout->extents, {4, 6});
    assert_vector_near(qmomkit::internal::extract_wheeler_moments(slice), {100.0, 106.0, 112.0, 118.0});

    const qmomkit::internal::DenseMatrix fixed_rhs = qmomkit::internal::make_branch_rhs_matrix(slice, 2);
    assert(fixed_rhs.rows() == 2);
    assert(fixed_rhs.cols() == 6);
    assert_matrix_value(fixed_rhs, 0, 5, 105.0);
    assert_matrix_value(fixed_rhs, 1, 5, 111.0);

    const qmomkit::internal::DenseMatrix adaptive_rhs = qmomkit::internal::make_branch_rhs_matrix(slice, 1);
    assert(adaptive_rhs.rows() == 1);
    assert(adaptive_rhs.cols() == 6);
    assert_matrix_value(adaptive_rhs, 0, 0, 100.0);
    assert_matrix_value(adaptive_rhs, 0, 5, 105.0);
    const qmomkit::internal::BranchMomentChildLayout& child_case = qmomkit::internal::child_layout_for_parent_node_count(*slice.layout, 1);
    assert(child_case.child_layout);
    assert_vector_equal(child_case.child_layout->extents, {6});
    assert(child_case.child_layout->indices.size() == 6);
}

void test_moment_data_recursive_rejects_missing_values()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    qmomkit::CqmomMomentData moment_data = qmomkit::extract_cqmom_moments({{4, 4}, indexed_values({4, 4})}, layout);
    moment_data.values.pop_back();

    const qmomkit::internal::RecursiveQuadratureResult result =
        qmomkit::internal::compute_recursive_quadrature_from_moment_data({hierarchy, moment_data, {}, {}});

    assert(!result.ok());
    assert(result.status == qmomkit::internal::RecursiveQuadratureStatus::invalid_input);
    assert(result.message == "CQMOM moment data values must match the layout size.");
}

void test_moment_data_recursive_rejects_non_finite_values()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    qmomkit::CqmomMomentData moment_data = qmomkit::extract_cqmom_moments({{4, 4}, indexed_values({4, 4})}, layout);
    moment_data.values[0] = std::numeric_limits<double>::infinity();

    const qmomkit::internal::RecursiveQuadratureResult result =
        qmomkit::internal::compute_recursive_quadrature_from_moment_data({hierarchy, moment_data, {}, {}});

    assert(!result.ok());
    assert(result.status == qmomkit::internal::RecursiveQuadratureStatus::invalid_input);
    assert(result.message == "CQMOM moment data values must be finite.");
}

void test_moment_data_recursive_rejects_incompatible_layout()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2});
    const qmomkit::PhaseSpaceHierarchy incompatible_hierarchy({2, 2, 2}, {2, 0, 1});
    const qmomkit::CqmomMomentLayout incompatible_layout = qmomkit::create_cqmom_moment_layout(incompatible_hierarchy);
    const qmomkit::CqmomMomentData moment_data = qmomkit::extract_cqmom_moments({{4, 4, 4}, indexed_values({4, 4, 4})}, incompatible_layout);

    const qmomkit::internal::RecursiveQuadratureResult result =
        qmomkit::internal::compute_recursive_quadrature_from_moment_data({hierarchy, moment_data, {}, {}});

    assert(!result.ok());
    assert(result.status == qmomkit::internal::RecursiveQuadratureStatus::invalid_input);
    assert(result.message == "CQMOM moment data layout must match the recursive hierarchy.");
}

void test_recursive_bounded_support_accepts_feasible_branch()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.2, 0.3, 0.1, 0.4};
    const std::vector<std::vector<double>> abscissas{{0.2, 0.2, 0.8, 0.8}, {-1.0, 1.0, -1.0, 1.0}};
    const std::vector<double> moments = raw_moments_from_terminal_quadrature(abscissas, weights, extents);
    const qmomkit::CqmomAbscissaBounds support = bounded_support_for_dimension(2, 0, {0.0, 1.0});

    const qmomkit::internal::RecursiveQuadratureResult result = qmomkit::internal::compute_recursive_quadrature_from_moment_data(
        moment_data_request(hierarchy, {extents, moments}, {}, qmomkit::WheelerOptions{false}, qmomkit::CqmomDiagnosticLevel::detailed,
                            qmomkit::VandermondeSolver::dense, support));

    assert(result.ok());
    assert(result.quadrature.has_rectangular_storage());
    for (const qmomkit::MultiIndex& node_index : qmomkit::enumerate_quadrature_node_multi_indices(hierarchy))
    {
        const double x = result.quadrature.abscissa(0, node_index);
        assert(x >= -1.0e-12);
        assert(x <= 1.0 + 1.0e-12);
    }
    assert(!result.diagnostics.empty());
    assert(result.diagnostics.front().bounded_support_checked);
    assert(result.diagnostics.front().bounded_candidate_node_count == 2);
}

void test_recursive_shifted_bounded_support_accepts_feasible_fixed_branch()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 1});
    const std::vector<std::size_t> extents{4, 2};
    const qmomkit::AbscissaBounds interval{2.0, 5.0};
    const std::vector<double> weights{0.4, 0.6};
    const std::vector<std::vector<double>> abscissas{{2.25, 4.5}, {0.0, 0.0}};
    const std::vector<double> moments = raw_moments_from_terminal_quadrature(abscissas, weights, extents);
    const qmomkit::CqmomAbscissaBounds support = bounded_support_for_dimension(2, 0, interval);

    const qmomkit::internal::RecursiveQuadratureResult result = qmomkit::internal::compute_recursive_quadrature_from_moment_data(
        moment_data_request(hierarchy, {extents, moments}, {}, qmomkit::WheelerOptions{false}, qmomkit::CqmomDiagnosticLevel::detailed,
                            qmomkit::VandermondeSolver::dense, support));

    assert(result.ok());
    assert(result.quadrature.has_rectangular_storage());
    assert_terminal_coordinate_inside_interval(result, 0, interval, 1.0e-12);
    assert_terminal_weights_positive(result);
    assert_raw_moments_reconstructed(result, moments, extents, 1.0e-10);

    const std::vector<double> reconstructed_prefix =
        raw_moments_from_terminal_quadrature({result.quadrature.abscissas()[0]}, result.quadrature.weights(), {4});
    const std::vector<double> expected_prefix = raw_moments_from_terminal_quadrature({abscissas[0]}, weights, {4});
    for (std::size_t order = 0; order < expected_prefix.size(); ++order)
    {
        assert_near(reconstructed_prefix[order], expected_prefix[order], 1.0e-10);
    }
}

void test_recursive_shifted_bounded_support_rejects_fixed_and_reduces_adaptive()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 1});
    const std::vector<std::size_t> extents{4, 2};
    const qmomkit::AbscissaBounds interval{2.0, 5.0};
    const std::vector<double> weights{0.5, 0.5};
    const std::vector<std::vector<double>> abscissas{{3.0, 7.0}, {0.0, 0.0}};
    const std::vector<double> moments = raw_moments_from_terminal_quadrature(abscissas, weights, extents);
    const qmomkit::CqmomAbscissaBounds support = bounded_support_for_dimension(2, 0, interval);

    const qmomkit::internal::RecursiveQuadratureResult fixed = qmomkit::internal::compute_recursive_quadrature_from_moment_data(
        moment_data_request(hierarchy, {extents, moments}, {}, qmomkit::WheelerOptions{false}, qmomkit::CqmomDiagnosticLevel::detailed,
                            qmomkit::VandermondeSolver::dense, support));
    assert(!fixed.ok());
    assert(fixed.failed_conditioning_level == 0);
    assert(!fixed.diagnostics.empty());
    assert(fixed.diagnostics.front().bounded_support_checked);
    assert(fixed.diagnostics.front().bounded_returned_nodes.size() == 2);
    assert(std::any_of(fixed.diagnostics.front().bounded_returned_nodes.begin(), fixed.diagnostics.front().bounded_returned_nodes.end(),
                       [interval](const double node) { return node > interval.upper; }));
    assert(fixed.message.find("returned nodes outside") != std::string::npos);

    const qmomkit::WheelerOptions adaptive_wheeler{true, 1.0e-10, 1.0e-10, 1.0e-12};
    const qmomkit::internal::RecursiveQuadratureResult adaptive = qmomkit::internal::compute_recursive_quadrature_from_moment_data(
        moment_data_request(hierarchy, {extents, moments}, {}, adaptive_wheeler, qmomkit::CqmomDiagnosticLevel::detailed,
                            qmomkit::VandermondeSolver::dense, support));
    assert(adaptive.ok());
    assert(adaptive.quadrature.total_nodes() == 1);
    assert_terminal_coordinate_inside_interval(adaptive, 0, interval, 1.0e-12);
    assert(!adaptive.diagnostics.empty());
    assert(adaptive.diagnostics.front().bounded_candidate_node_count == 2);
    assert(adaptive.diagnostics.front().bounded_accepted_fallback_node_count == 1);
    assert(adaptive.diagnostics.front().adaptive_node_count_reduced);
}

void test_recursive_bounded_support_rejects_fixed_infeasible_count()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.2, 0.3, 0.1, 0.4};
    const std::vector<std::vector<double>> abscissas{{2.0, 2.0, 3.0, 3.0}, {-1.0, 1.0, -1.0, 1.0}};
    const std::vector<double> moments = raw_moments_from_terminal_quadrature(abscissas, weights, extents);
    const qmomkit::CqmomAbscissaBounds support = bounded_support_for_dimension(2, 0, {0.0, 1.0});

    const qmomkit::internal::RecursiveQuadratureResult result = qmomkit::internal::compute_recursive_quadrature_from_moment_data(
        moment_data_request(hierarchy, {extents, moments}, {}, qmomkit::WheelerOptions{false}, qmomkit::CqmomDiagnosticLevel::detailed,
                            qmomkit::VandermondeSolver::dense, support));

    assert(!result.ok());
    assert(result.status == qmomkit::internal::RecursiveQuadratureStatus::wheeler_failed);
    assert(result.failed_conditioning_level == 0);
    assert(result.diagnostics.size() == 1);
    assert(result.diagnostics.front().bounded_support_checked);
    assert(result.diagnostics.front().bounded_candidate_attempt_type == "fixed");
    assert(result.diagnostics.front().bounded_candidate_node_count == 2);
    assert(result.diagnostics.front().bounded_returned_nodes.size() == 2);
    assert(std::all_of(result.diagnostics.front().bounded_returned_nodes.begin(), result.diagnostics.front().bounded_returned_nodes.end(),
                       [](const double node) { return node > 1.0; }));
    assert(result.message.find("returned nodes outside") != std::string::npos);
}

void test_recursive_bounded_support_adaptive_falls_back_to_supported_count()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({3, 1});
    const std::vector<std::size_t> extents{6, 2};
    const std::vector<double> valid_prefix{1.0, 0.5, 0.3125, 0.21875};
    std::vector<double> moments(qmomkit::product_extents(extents), 0.0);
    for (std::size_t order = 0; order < valid_prefix.size(); ++order)
    {
        moments[qmomkit::flatten_multi_index({order, 0}, extents)] = valid_prefix[order];
    }
    moments[qmomkit::flatten_multi_index({4, 0}, extents)] = 10.0;
    moments[qmomkit::flatten_multi_index({5, 0}, extents)] = 10.0;
    const qmomkit::WheelerOptions wheeler_options{true, 1.0e-10, 1.0e-10, 1.0e-12};
    const qmomkit::CqmomAbscissaBounds support = bounded_support_for_dimension(2, 0, {0.0, 1.0});

    const qmomkit::internal::RecursiveQuadratureResult result = qmomkit::internal::compute_recursive_quadrature_from_moment_data(moment_data_request(
        hierarchy, {extents, moments}, {}, wheeler_options, qmomkit::CqmomDiagnosticLevel::detailed, qmomkit::VandermondeSolver::dense, support));

    assert(result.ok());
    assert(result.quadrature.total_nodes() == 2);
    assert(!result.diagnostics.empty());
    const qmomkit::internal::RecursiveQuadratureBranchDiagnostics& diagnostic = result.diagnostics.front();
    assert(diagnostic.bounded_support_checked);
    assert(diagnostic.bounded_candidate_node_count == 3);
    assert(diagnostic.bounded_accepted_fallback_node_count == 2);
    assert(diagnostic.adaptive_node_count_reduced);
}

void test_recursive_bounded_support_adaptive_two_node_uses_span_separation()
{
    // Per docs/developer/qmom_wheeler.md, the two-node separation ratio is the
    // node distance divided by the larger abscissa magnitude: 0.01 / 1.01 ~
    // 0.0099. That is below eabs = 0.1, so this candidate is correctly reduced
    // to one node instead of being kept at two.
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 1});
    const std::vector<std::size_t> extents{4, 2};
    const std::vector<double> moments = raw_moments_from_terminal_quadrature({{1.0, 1.01}, {0.0, 0.0}}, {0.5, 0.5}, extents);
    const qmomkit::WheelerOptions wheeler_options{true, 1.0e-10, 0.1, 1.0e-12};
    const qmomkit::CqmomAbscissaBounds support = bounded_support_for_dimension(2, 0, {0.0, 2.0});

    const qmomkit::internal::RecursiveQuadratureResult result = qmomkit::internal::compute_recursive_quadrature_from_moment_data(moment_data_request(
        hierarchy, {extents, moments}, {}, wheeler_options, qmomkit::CqmomDiagnosticLevel::detailed, qmomkit::VandermondeSolver::dense, support));

    assert(result.ok());
    assert(result.quadrature.total_nodes() == 1);
    assert(!result.diagnostics.empty());
    assert(result.diagnostics.front().bounded_candidate_node_count == 2);
    assert(result.diagnostics.front().bounded_accepted_fallback_node_count == 1);
    assert(result.diagnostics.front().adaptive_node_count_reduced);
}

void test_recursive_bounded_support_adaptive_applies_relative_weight_check()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 1});
    const std::vector<std::size_t> extents{4, 2};
    const std::vector<double> moments = raw_moments_from_terminal_quadrature({{0.2, 0.8}, {0.0, 0.0}}, {0.99, 0.01}, extents);
    const qmomkit::WheelerOptions wheeler_options{true, 0.1, 1.0e-10, 1.0e-12};
    const qmomkit::CqmomAbscissaBounds support = bounded_support_for_dimension(2, 0, {0.0, 1.0});

    const qmomkit::internal::RecursiveQuadratureResult result = qmomkit::internal::compute_recursive_quadrature_from_moment_data(moment_data_request(
        hierarchy, {extents, moments}, {}, wheeler_options, qmomkit::CqmomDiagnosticLevel::detailed, qmomkit::VandermondeSolver::dense, support));

    assert(result.ok());
    assert(result.quadrature.total_nodes() == 1);
    assert(!result.diagnostics.empty());
    assert(result.diagnostics.front().bounded_candidate_node_count == 2);
    assert(result.diagnostics.front().bounded_accepted_fallback_node_count == 1);
    assert(result.diagnostics.front().adaptive_node_count_reduced);
}

void test_recursive_bounded_support_rejects_growth_and_preserves_current_count()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({3, 1});
    const std::vector<std::size_t> extents{6, 2};
    const std::vector<double> valid_prefix{1.0, 0.5, 0.3125, 0.21875};
    std::vector<double> moments(qmomkit::product_extents(extents), 0.0);
    for (std::size_t order = 0; order < valid_prefix.size(); ++order)
    {
        moments[qmomkit::flatten_multi_index({order, 0}, extents)] = valid_prefix[order];
    }
    moments[qmomkit::flatten_multi_index({4, 0}, extents)] = 10.0;
    moments[qmomkit::flatten_multi_index({5, 0}, extents)] = 10.0;
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData moment_data = qmomkit::extract_cqmom_moments({extents, moments}, layout);
    const std::vector<qmomkit::CqmomActiveBranchCount> active_counts{
        {0, {}, 0, 3, 2}, {1, qmomkit::ConditioningPath({0}), 1, 1, 1}, {1, qmomkit::ConditioningPath({1}), 1, 1, 1}};
    const qmomkit::WheelerOptions wheeler_options{true, 1.0e-10, 1.0e-10, 1.0e-12};
    const qmomkit::CqmomAbscissaBounds support = bounded_support_for_dimension(2, 0, {0.0, 1.0});

    const qmomkit::internal::RecursiveQuadratureResult result =
        qmomkit::internal::compute_recursive_quadrature_from_moment_data({hierarchy,
                                                                          moment_data,
                                                                          {},
                                                                          wheeler_options,
                                                                          qmomkit::CqmomDiagnosticLevel::detailed,
                                                                          qmomkit::VandermondeSolver::dense,
                                                                          active_counts,
                                                                          true,
                                                                          support});

    assert(result.ok());
    assert(result.quadrature.total_nodes() == 2);
    const qmomkit::internal::RecursiveQuadratureBranchDiagnostics& diagnostic = result.diagnostics.front();
    assert(diagnostic.bounded_support_checked);
    assert(diagnostic.bounded_candidate_attempt_type == "growth");
    assert(diagnostic.bounded_candidate_node_count == 3);
    assert(diagnostic.bounded_accepted_fallback_node_count == 2);
    assert(!diagnostic.adaptive_node_count_reduced);
}

void test_recursive_bounded_support_uses_phase_space_coordinate_index()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {1, 0});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.2, 0.3, 0.1, 0.4};
    const std::vector<std::vector<double>> abscissas{{-1.0, 1.0, -1.0, 1.0}, {0.2, 0.2, 0.8, 0.8}};
    const std::vector<double> moments = raw_moments_from_terminal_quadrature(abscissas, weights, extents);
    const qmomkit::CqmomAbscissaBounds support = bounded_support_for_dimension(2, 1, {0.0, 1.0});

    const qmomkit::internal::RecursiveQuadratureResult result = qmomkit::internal::compute_recursive_quadrature_from_moment_data(
        moment_data_request(hierarchy, {extents, moments}, {}, qmomkit::WheelerOptions{false}, qmomkit::CqmomDiagnosticLevel::detailed,
                            qmomkit::VandermondeSolver::dense, support));

    assert(result.ok());
    assert(!result.diagnostics.empty());
    assert(result.diagnostics.front().conditioning_level == 0);
    assert(result.diagnostics.front().phase_space_dimension == 1);
    assert(result.diagnostics.front().bounded_support_checked);
}

void test_fixed_public_terminal_paths_align_with_physical_storage()
{
    // The manufactured terminals are deliberately distinct.  The vectors are
    // in physical row-major order: (x, y) = 00, 01, 10, 11.
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {1, 0});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.2, 0.3, 0.1, 0.4};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}};
    const std::vector<double> moments = raw_moments_from_terminal_quadrature(abscissas, weights, extents);

    const std::vector<qmomkit::ConditioningPath> identity_paths{qmomkit::ConditioningPath({0, 0}), qmomkit::ConditioningPath({0, 1}),
                                                                qmomkit::ConditioningPath({1, 0}), qmomkit::ConditioningPath({1, 1})};
    const qmomkit::PhaseSpaceHierarchy identity_hierarchy({2, 2});
    const qmomkit::CqmomResult identity = qmomkit::compute_cqmom(
        qmomkit::CqmomCase{qmomkit::extract_cqmom_moments({extents, moments}, qmomkit::create_cqmom_moment_layout(identity_hierarchy)),
                           {},
                           qmomkit::QmomOptions{false}});
    assert(identity.ok());
    assert(identity.quadrature.terminal_paths().to_conditioning_paths() == identity_paths);

    const qmomkit::CqmomCase cqmom_case{
        qmomkit::extract_cqmom_moments({extents, moments}, qmomkit::create_cqmom_moment_layout(hierarchy)), {}, qmomkit::QmomOptions{false}};
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(cqmom_case);
    assert(result.ok());
    assert(result.quadrature.has_rectangular_storage());

    // Row-major storage requires paths in conditioning-level order (y, x):
    // physical 00, 01, 10, 11 map to paths 00, 10, 01, 11.
    const std::vector<qmomkit::ConditioningPath> expected_paths{qmomkit::ConditioningPath({0, 0}), qmomkit::ConditioningPath({1, 0}),
                                                                qmomkit::ConditioningPath({0, 1}), qmomkit::ConditioningPath({1, 1})};

    // The quadrature values themselves retain their independently supplied
    // physical-node identities and reconstruct every input raw moment.
    for (std::size_t node = 0; node < weights.size(); ++node)
    {
        assert_near(result.quadrature.weights()[node], weights[node], 1.0e-12);
        assert_near(result.quadrature.abscissas()[0][node], abscissas[0][node], 1.0e-12);
        assert_near(result.quadrature.abscissas()[1][node], abscissas[1][node], 1.0e-12);
    }
    const qmomkit::RawMomentTensor reconstructed = qmomkit::reconstruct_raw_moments(result.quadrature, extents);
    for (std::size_t moment = 0; moment < moments.size(); ++moment)
    {
        assert_near(reconstructed.values[moment], moments[moment], 1.0e-12);
    }

    assert(result.quadrature.terminal_paths().to_conditioning_paths() == expected_paths);
}

void test_recursive_bounded_support_applies_to_conditional_phase_space_dimension()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {0, 1});
    const std::vector<std::size_t> extents{4, 4};
    const qmomkit::AbscissaBounds interval{2.0, 5.0};
    const std::vector<double> weights{0.2, 0.3, 0.1, 0.4};
    const std::vector<std::vector<double>> abscissas{{-3.0, -3.0, 8.0, 8.0}, {2.25, 4.5, 2.5, 4.75}};
    const std::vector<double> moments = raw_moments_from_terminal_quadrature(abscissas, weights, extents);
    const qmomkit::CqmomAbscissaBounds support = bounded_support_for_dimension(2, 1, interval);

    const qmomkit::internal::RecursiveQuadratureResult result = qmomkit::internal::compute_recursive_quadrature_from_moment_data(
        moment_data_request(hierarchy, {extents, moments}, {}, qmomkit::WheelerOptions{false}, qmomkit::CqmomDiagnosticLevel::detailed,
                            qmomkit::VandermondeSolver::dense, support));

    assert(result.ok());
    assert_terminal_coordinate_inside_interval(result, 1, interval, 1.0e-12);
    assert(result.quadrature.abscissas()[0].front() < interval.lower);
    assert(result.quadrature.abscissas()[0].back() > interval.upper);
    assert_raw_moments_reconstructed(result, moments, extents, 1.0e-9);

    bool checked_conditional_dimension = false;
    for (const qmomkit::internal::RecursiveQuadratureBranchDiagnostics& diagnostic : result.diagnostics)
    {
        if (diagnostic.conditioning_level == 1 && diagnostic.phase_space_dimension == 1)
        {
            assert(diagnostic.bounded_support_checked);
            checked_conditional_dimension = true;
        }
        if (diagnostic.phase_space_dimension == 0)
        {
            assert(!diagnostic.bounded_support_checked);
        }
    }
    assert(checked_conditional_dimension);
}

void test_recursive_bounded_support_returned_node_tolerance_and_no_clipping()
{
    const qmomkit::AbscissaBounds interval{2.0, 5.0};
    // The acceptance margin is relative to the larger bound magnitude.
    const double tolerance = 1.0e-12 * std::max(std::abs(interval.lower), std::abs(interval.upper));
    const double weight = 0.01;
    const qmomkit::PhaseSpaceHierarchy hierarchy({1});
    const qmomkit::CqmomAbscissaBounds support = bounded_support_for_dimension(1, 0, interval);

    const auto run = [&](const double node, const double node_weight)
    {
        const std::vector<std::size_t> extents{2};
        const std::vector<double> moments = raw_moments_from_terminal_quadrature({{node}}, {node_weight}, extents);
        return qmomkit::internal::compute_recursive_quadrature_from_moment_data(
            moment_data_request(hierarchy, {extents, moments}, {}, qmomkit::WheelerOptions{false}, qmomkit::CqmomDiagnosticLevel::detailed,
                                qmomkit::VandermondeSolver::dense, support));
    };

    for (const double node : {interval.lower, interval.upper, interval.lower - 0.5 * tolerance, interval.upper + 0.5 * tolerance})
    {
        const qmomkit::internal::RecursiveQuadratureResult result = run(node, weight);
        assert(result.ok());
        assert(result.quadrature.weights().size() == 1);
        assert(result.quadrature.abscissas()[0].size() == 1);
        assert(std::isfinite(result.quadrature.weights()[0]));
        assert(std::isfinite(result.quadrature.abscissas()[0][0]));
        assert_near(result.quadrature.weights()[0], weight, 0.0);
        assert_near(result.quadrature.abscissas()[0][0], node, 0.01 * tolerance);
        if (node < interval.lower)
        {
            assert(result.quadrature.abscissas()[0][0] < interval.lower);
            assert(interval.lower - result.quadrature.abscissas()[0][0] <= tolerance);
        }
        if (node > interval.upper)
        {
            assert(result.quadrature.abscissas()[0][0] > interval.upper);
            assert(result.quadrature.abscissas()[0][0] - interval.upper <= tolerance);
        }
    }

    const qmomkit::internal::RecursiveQuadratureResult lower_rejected = run(interval.lower - 2.0 * tolerance, weight);
    const qmomkit::internal::RecursiveQuadratureResult upper_rejected = run(interval.upper + 2.0 * tolerance, weight);
    assert(!lower_rejected.ok());
    assert(!upper_rejected.ok());
    assert(!lower_rejected.diagnostics.empty());
    assert(!upper_rejected.diagnostics.empty());
    assert(lower_rejected.diagnostics.front().bounded_support_checked);
    assert(upper_rejected.diagnostics.front().bounded_support_checked);
    assert(lower_rejected.diagnostics.front().bounded_returned_nodes.size() == 1);
    assert(upper_rejected.diagnostics.front().bounded_returned_nodes.size() == 1);
    assert(lower_rejected.message.find("returned nodes outside") != std::string::npos);
    assert(upper_rejected.message.find("returned nodes outside") != std::string::npos);
}

void test_one_dimensional_hierarchy_reduces_to_stage_1_ordering()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({3});

    assert(hierarchy.number_of_dimensions() == 1);
    assert(hierarchy.total_nodes() == 3);
    assert_vector_equal(hierarchy.node_counts(), {3});
    assert_vector_equal(hierarchy.conditioning_order(), {0});
    assert_vector_equal(hierarchy.moment_extents(), {6});

    const std::vector<qmomkit::MultiIndex> moment_indices = qmomkit::enumerate_multi_indices(hierarchy.moment_extents());
    const std::vector<qmomkit::MultiIndex> node_indices = qmomkit::enumerate_quadrature_node_multi_indices(hierarchy);

    assert(moment_indices.size() == 6);
    assert(node_indices.size() == 3);
    assert_vector_equal(moment_indices[0], {0});
    assert_vector_equal(moment_indices[5], {5});
    assert_vector_equal(node_indices[0], {0});
    assert_vector_equal(node_indices[2], {2});

    qmomkit::internal::RecursiveQuadrature quadrature(hierarchy);
    quadrature.set_weight({0}, 0.25);
    quadrature.set_abscissa(0, {0}, -1.0);
    quadrature.set_weight({1}, 0.75);
    quadrature.set_abscissa(0, {1}, 2.0);
    quadrature.set_weight({2}, 0.5);
    quadrature.set_abscissa(0, {2}, 3.5);

    assert(quadrature.weight({0}) == 0.25);
    assert(quadrature.abscissa(0, {1}) == 2.0);
    assert_vector_near(quadrature.node({2}), {3.5});
}

void test_two_dimensional_hierarchy_indexing()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 3});

    assert(hierarchy.total_nodes() == 6);
    assert_vector_equal(hierarchy.moment_extents(), {4, 6});
    assert(hierarchy.flatten_node_index({0, 0}) == 0);
    assert(hierarchy.flatten_node_index({0, 1}) == 1);
    assert(hierarchy.flatten_node_index({0, 2}) == 2);
    assert(hierarchy.flatten_node_index({1, 0}) == 3);
    assert(hierarchy.flatten_node_index({1, 2}) == 5);
    assert_vector_equal(hierarchy.unflatten_node_index(4), {1, 1});

    const std::vector<qmomkit::MultiIndex> node_indices = qmomkit::enumerate_quadrature_node_multi_indices(hierarchy);
    assert(node_indices.size() == 6);
    assert_vector_equal(node_indices[0], {0, 0});
    assert_vector_equal(node_indices[5], {1, 2});

    const std::vector<qmomkit::MultiIndex> moment_indices = qmomkit::enumerate_multi_indices(hierarchy.moment_extents());
    assert(moment_indices.size() == 24);
    assert_vector_equal(moment_indices[0], {0, 0});
    assert_vector_equal(moment_indices[23], {3, 5});
}

void test_three_dimensional_hierarchy_indexing()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2});

    assert(hierarchy.total_nodes() == 8);
    assert(hierarchy.flatten_node_index({0, 0, 0}) == 0);
    assert(hierarchy.flatten_node_index({0, 0, 1}) == 1);
    assert(hierarchy.flatten_node_index({0, 1, 0}) == 2);
    assert(hierarchy.flatten_node_index({1, 0, 0}) == 4);
    assert(hierarchy.flatten_node_index({1, 1, 1}) == 7);
    assert_vector_equal(hierarchy.unflatten_node_index(6), {1, 1, 0});

    const std::vector<qmomkit::MultiIndex> moment_indices = qmomkit::enumerate_multi_indices(hierarchy.moment_extents());
    assert(moment_indices.size() == 64);
    assert_vector_equal(moment_indices.front(), {0, 0, 0});
    assert_vector_equal(moment_indices.back(), {3, 3, 3});
}

void test_non_uniform_node_counts()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 3, 2});

    assert(hierarchy.number_of_dimensions() == 3);
    assert(hierarchy.total_nodes() == 12);
    assert_vector_equal(hierarchy.moment_extents(), {4, 6, 4});
    assert(hierarchy.flatten_node_index({0, 0, 0}) == 0);
    assert(hierarchy.flatten_node_index({0, 1, 0}) == 2);
    assert(hierarchy.flatten_node_index({1, 0, 0}) == 6);
    assert(hierarchy.flatten_node_index({1, 2, 1}) == 11);

    qmomkit::internal::RecursiveQuadrature quadrature(hierarchy);
    for (const qmomkit::MultiIndex& node_index : qmomkit::enumerate_quadrature_node_multi_indices(hierarchy))
    {
        const double flat = static_cast<double>(hierarchy.flatten_node_index(node_index));
        quadrature.set_weight(node_index, flat + 1.0);
        quadrature.set_abscissa(0, node_index, 10.0 + flat);
        quadrature.set_abscissa(1, node_index, 20.0 + flat);
        quadrature.set_abscissa(2, node_index, 30.0 + flat);
    }

    assert(quadrature.weight({1, 2, 1}) == 12.0);
    assert_vector_near(quadrature.node({1, 2, 1}), {21.0, 31.0, 41.0});
}

void test_identity_conditioning_order_matches_phase_space_order()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 3, 2}, {0, 1, 2});

    assert_vector_equal(hierarchy.conditioning_order(), {0, 1, 2});
    assert(hierarchy.flatten_node_index({0, 0, 0}) == 0);
    assert(hierarchy.flatten_node_index({0, 0, 1}) == 1);
    assert(hierarchy.flatten_node_index({0, 1, 0}) == 2);
    assert(hierarchy.flatten_node_index({1, 0, 0}) == 6);
    assert_vector_equal(hierarchy.unflatten_node_index(11), {1, 2, 1});

    const std::vector<qmomkit::MultiIndex> node_indices = qmomkit::enumerate_quadrature_node_multi_indices(hierarchy);
    assert_vector_equal(node_indices[0], {0, 0, 0});
    assert_vector_equal(node_indices[1], {0, 0, 1});
    assert_vector_equal(node_indices[2], {0, 1, 0});
    assert_vector_equal(node_indices[11], {1, 2, 1});
}

void test_non_identity_conditioning_order_controls_recursive_traversal()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 3, 2}, {2, 0, 1});

    assert_vector_equal(hierarchy.conditioning_order(), {2, 0, 1});
    // Flat node indexing remains in phase-space dimension order.
    assert(hierarchy.flatten_node_index({0, 0, 0}) == 0);
    assert(hierarchy.flatten_node_index({0, 0, 1}) == 1);
    assert(hierarchy.flatten_node_index({0, 1, 0}) == 2);
    assert(hierarchy.flatten_node_index({1, 0, 0}) == 6);
    assert_vector_equal(hierarchy.unflatten_node_index(2), {0, 1, 0});

    // Recursive traversal follows conditioning levels: phase-space dimension 2,
    // then 0, then 1.
    const std::vector<qmomkit::MultiIndex> node_indices = qmomkit::enumerate_quadrature_node_multi_indices(hierarchy);
    assert_vector_equal(node_indices[0], {0, 0, 0});
    assert_vector_equal(node_indices[1], {0, 1, 0});
    assert_vector_equal(node_indices[2], {0, 2, 0});
    assert_vector_equal(node_indices[3], {1, 0, 0});
    assert_vector_equal(node_indices[6], {0, 0, 1});
    assert_vector_equal(node_indices[11], {1, 2, 1});
}

void test_reversed_conditioning_order_controls_recursive_traversal()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 3, 2}, {2, 1, 0});

    assert_vector_equal(hierarchy.conditioning_order(), {2, 1, 0});
    // Flat node indexing remains in phase-space dimension order.
    assert(hierarchy.flatten_node_index({0, 0, 0}) == 0);
    assert(hierarchy.flatten_node_index({0, 0, 1}) == 1);
    assert(hierarchy.flatten_node_index({0, 1, 0}) == 2);
    assert(hierarchy.flatten_node_index({1, 0, 0}) == 6);
    assert(hierarchy.flatten_node_index({1, 2, 1}) == 11);
    assert_vector_equal(hierarchy.unflatten_node_index(1), {0, 0, 1});
    assert_vector_equal(hierarchy.unflatten_node_index(6), {1, 0, 0});

    const std::vector<qmomkit::MultiIndex> node_indices = qmomkit::enumerate_quadrature_node_multi_indices(hierarchy);
    assert_vector_equal(node_indices[0], {0, 0, 0});
    assert_vector_equal(node_indices[1], {1, 0, 0});
    assert_vector_equal(node_indices[2], {0, 1, 0});
    assert_vector_equal(node_indices[6], {0, 0, 1});
    assert_vector_equal(node_indices[11], {1, 2, 1});
}

// Terminal paths all reach the deepest conditioning level, so they are stored
// contiguously: one allocation per table instead of one per path, and one
// allocation when a quadrature carries them through node advection.
void test_conditioning_path_table_round_trips_and_views_rows()
{
    const std::vector<qmomkit::ConditioningPath> paths{qmomkit::ConditioningPath({0, 1}), qmomkit::ConditioningPath({1, 0}),
                                                       qmomkit::ConditioningPath({1, 1})};
    const qmomkit::ConditioningPathTable table(paths);

    assert(table.size() == 3);
    assert(table.levels() == 2);
    assert(!table.empty());
    assert(table.node_indices() == std::vector<std::size_t>({0, 1, 1, 0, 1, 1}));

    for (std::size_t path = 0; path < paths.size(); ++path)
    {
        const qmomkit::ConditioningPathView view = table[path];
        assert(view.size() == 2);
        assert(view == qmomkit::ConditioningPathView(paths[path]));
        assert(view.to_conditioning_path() == paths[path]);
        for (std::size_t level = 0; level < view.size(); ++level)
        {
            assert(view.node_index_at_level(level) == paths[path].node_index_at_level(level));
        }
    }
    assert(table.to_conditioning_paths() == paths);

    // Rows are adjacent in memory: this is what removes the per-path allocation.
    assert(table[1].data() == table[0].data() + table.levels());
}

void test_conditioning_path_table_rejects_uneven_depth_and_checks_bounds()
{
    const qmomkit::ConditioningPathTable flat(2, std::vector<std::size_t>{1, 2, 3, 4});
    assert(flat.size() == 2);
    assert(flat[1] == qmomkit::ConditioningPathView(qmomkit::ConditioningPath({3, 4})));

    assert_throws_invalid_argument([] { (void)qmomkit::ConditioningPathTable(2, std::vector<std::size_t>{1, 2, 3}); });
    assert_throws_invalid_argument([] { (void)qmomkit::ConditioningPathTable(0, std::vector<std::size_t>{1}); });
    assert_throws_invalid_argument([] { (void)qmomkit::ConditioningPathTable({qmomkit::ConditioningPath({0, 1}), qmomkit::ConditioningPath({0})}); });

    assert_throws_out_of_range([&flat] { (void)flat.at(2); });
    assert_throws_out_of_range([&flat] { (void)flat[0].node_index_at_level(2); });
}

void test_conditioning_path_table_grows_iterates_and_compares()
{
    qmomkit::ConditioningPathTable table;
    assert(table.empty());
    assert(table.levels() == 0);
    assert(table.begin() == table.end());

    table.reserve(2);
    table.push_back(qmomkit::ConditioningPath({2, 3}));
    assert(table.levels() == 2);
    const qmomkit::ConditioningPathTable source({qmomkit::ConditioningPath({4, 5})});
    table.push_back(source[0]);
    assert(table.size() == 2);
    assert_throws_invalid_argument([&table] { table.push_back(qmomkit::ConditioningPath({1, 2, 3})); });

    std::size_t visited = 0;
    for (const qmomkit::ConditioningPathView row : table)
    {
        assert(row.size() == 2);
        ++visited;
    }
    assert(visited == 2);

    assert(table == qmomkit::ConditioningPathTable({qmomkit::ConditioningPath({2, 3}), qmomkit::ConditioningPath({4, 5})}));
    assert(table != qmomkit::ConditioningPathTable({qmomkit::ConditioningPath({2, 3})}));

    // clear() keeps the depth so the table can be refilled with the same shape.
    table.clear();
    assert(table.empty());
    assert(table.levels() == 2);
}

// A view maps to phase-space order exactly as the owning path does; the owning
// path delegates to it so the mapping exists once.
void test_conditioning_path_view_maps_to_phase_space_node_index()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 3, 4}, {2, 0, 1});
    const qmomkit::ConditioningPath path({3, 1, 2});
    const qmomkit::ConditioningPathTable table({path});

    assert(table[0].phase_space_node_index(hierarchy) == path.phase_space_node_index(hierarchy));
    assert(qmomkit::recursive_path_to_node_index(hierarchy, table[0]) == qmomkit::recursive_path_to_node_index(hierarchy, path));
}

void test_conditioning_path_maps_to_phase_space_node_index()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 3, 2}, {2, 0, 1});
    const qmomkit::ConditioningPath partial_path({1, 0});
    const qmomkit::ConditioningPath full_path({1, 0, 2});

    assert(partial_path.size() == 2);
    assert(partial_path.node_index_at_level(0) == 1);
    assert(partial_path.node_index_at_level(1) == 0);
    assert_vector_equal(partial_path.phase_space_node_index(hierarchy), {0, 0, 1});
    assert_vector_equal(full_path.phase_space_node_index(hierarchy), {0, 2, 1});
    assert_vector_equal(qmomkit::recursive_path_to_node_index(hierarchy, full_path), {0, 2, 1});
    assert(hierarchy.flatten_node_index(qmomkit::recursive_path_to_node_index(hierarchy, full_path)) == hierarchy.flatten_node_index({0, 2, 1}));
    assert(partial_path.with_child(2) == full_path);
    assert(qmomkit::ConditioningPath({0, 2}) < qmomkit::ConditioningPath({1, 0}));
    assert(!partial_path.to_string().empty());
    assert(!partial_path.describe(hierarchy).empty());
}

void test_invalid_hierarchy_inputs()
{
    assert_throws_invalid_argument([] { (void)qmomkit::PhaseSpaceHierarchy({}); });
    assert_throws_invalid_argument([] { (void)qmomkit::PhaseSpaceHierarchy({2, 0}); });
    assert_throws_invalid_argument([] { (void)qmomkit::PhaseSpaceHierarchy({2, 3}, {0}); });
    assert_throws_invalid_argument([] { (void)qmomkit::PhaseSpaceHierarchy({2, 3, 2}, {0, 0, 1}); });
    assert_throws_invalid_argument([] { (void)qmomkit::PhaseSpaceHierarchy({2, 3, 2}, {0, 1, 1}); });
    assert_throws_invalid_argument([] { (void)qmomkit::PhaseSpaceHierarchy({2, 3, 2}, {0, 1, 3}); });
}

void test_invalid_index_inputs()
{
    assert_throws_invalid_argument([] { (void)qmomkit::flatten_multi_index({0}, {2, 3}); });
    assert_throws_out_of_range([] { (void)qmomkit::flatten_multi_index({2, 0}, {2, 3}); });
    assert_throws_out_of_range([] { (void)qmomkit::unflatten_multi_index(6, {2, 3}); });
    assert_throws_invalid_argument([] { (void)qmomkit::product_extents({std::numeric_limits<std::size_t>::max(), std::size_t{2}}); });
    assert_throws_invalid_argument([] { (void)qmomkit::flatten_multi_index({0, 0}, {std::numeric_limits<std::size_t>::max(), std::size_t{2}}); });
    assert_throws_invalid_argument([] { (void)qmomkit::PhaseSpaceHierarchy({std::numeric_limits<std::size_t>::max()}).moment_extents(); });

    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 3});
    assert_throws_invalid_argument([&hierarchy] { (void)hierarchy.flatten_node_index({0}); });
    assert_throws_out_of_range([&hierarchy] { (void)hierarchy.flatten_node_index({0, 3}); });
    assert_throws_out_of_range([&hierarchy] { (void)hierarchy.unflatten_node_index(6); });

    qmomkit::internal::RecursiveQuadrature quadrature(hierarchy);
    assert_throws_out_of_range([&quadrature] { (void)quadrature.abscissa(2, {0, 0}); });
    assert_throws_out_of_range([&quadrature] { quadrature.set_abscissa(2, {0, 0}, 0.0); });
    assert_throws_invalid_argument([&hierarchy] { (void)qmomkit::internal::RecursiveQuadrature(hierarchy, {1.0}, {{0.0}}); });
    assert_throws_invalid_argument([&hierarchy] { (void)qmomkit::recursive_path_to_node_index(hierarchy, qmomkit::ConditioningPath({1})); });
    assert_throws_out_of_range([&hierarchy] { (void)qmomkit::ConditioningPath({3, 0}).phase_space_node_index(hierarchy); });
}

void test_recursive_failure_diagnostics_identify_failing_branch()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<double> parent_abscissas{-1.0, 1.0};
    const std::vector<double> parent_weights{0.5, 0.5};
    const std::vector<std::size_t> extents{4, 4};
    std::vector<double> moments(qmomkit::product_extents(extents), 0.0);
    const std::vector<std::vector<double>> conditional_moments{{1.0, 0.0, 1.0, 0.0}, {1.0, 0.0, -1.0, 0.0}};
    for (std::size_t flat_index = 0; flat_index < moments.size(); ++flat_index)
    {
        const qmomkit::MultiIndex moment_index = qmomkit::unflatten_multi_index(flat_index, extents);
        for (std::size_t parent_node = 0; parent_node < parent_abscissas.size(); ++parent_node)
        {
            moments[flat_index] += parent_weights[parent_node] * std::pow(parent_abscissas[parent_node], static_cast<int>(moment_index[0])) *
                                   conditional_moments[parent_node][moment_index[1]];
        }
    }

    const qmomkit::internal::RecursiveQuadratureResult result = qmomkit::internal::compute_recursive_quadrature_from_moment_data(
        moment_data_request(hierarchy, {extents, moments}, {}, qmomkit::WheelerOptions{false}));

    assert(!result.ok());
    assert(result.status == qmomkit::internal::RecursiveQuadratureStatus::wheeler_failed);
    assert(result.failed_conditioning_level == 1);
    assert(result.failed_path == qmomkit::ConditioningPath({1}));
    assert(!result.diagnostics.empty());
    assert(result.diagnostics.back().path == qmomkit::ConditioningPath({1}));
    assert(result.diagnostics.back().conditioning_level == 1);
    assert(result.diagnostics.back().wheeler_status == qmomkit::QuadratureStatus::not_realizable);
}

void test_recursive_dense_and_bp_equivalence_fixed()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2}, {2, 0, 1});
    const std::vector<std::size_t> extents{4, 4, 4};
    const std::vector<double> weights{0.05, 0.15, 0.20, 0.10, 0.12, 0.18, 0.08, 0.12};
    const std::vector<std::vector<double>> abscissas{{-1.5, 0.25, -0.75, 1.75, -1.25, 0.75, -0.10, 2.00},
                                                     {-2.0, -1.0, 0.5, 1.5, -1.5, -0.5, 1.0, 2.0},
                                                     {-0.80, -0.80, -0.20, -0.20, 0.60, 0.60, 1.20, 1.20}};
    const std::vector<double> moments = raw_moments_from_terminal_quadrature(abscissas, weights, extents);

    const qmomkit::internal::RecursiveQuadratureMomentDataRequest dense_request =
        moment_data_request(hierarchy, {extents, moments}, {}, qmomkit::WheelerOptions{false}, qmomkit::CqmomDiagnosticLevel::detailed,
                            qmomkit::VandermondeSolver::dense);
    const qmomkit::internal::RecursiveQuadratureMomentDataRequest bp_request =
        moment_data_request(hierarchy, {extents, moments}, {}, qmomkit::WheelerOptions{false}, qmomkit::CqmomDiagnosticLevel::detailed,
                            qmomkit::VandermondeSolver::bjorck_pereyra);

    const qmomkit::internal::RecursiveQuadratureResult dense = qmomkit::internal::compute_recursive_quadrature_from_moment_data(dense_request);
    const qmomkit::internal::RecursiveQuadratureResult bp = qmomkit::internal::compute_recursive_quadrature_from_moment_data(bp_request);

    assert_recursive_results_equivalent(bp, dense, 1.0e-9);
    assert(bp.ok());
    assert(!bp.diagnostics.empty());
    for (const qmomkit::internal::RecursiveQuadratureBranchDiagnostics& diagnostics : bp.diagnostics)
    {
        if (diagnostics.conditioning_level + 1 < hierarchy.number_of_dimensions())
        {
            assert(diagnostics.conditional_vandermonde_requested_solver == qmomkit::VandermondeSolver::bjorck_pereyra);
            assert(diagnostics.conditional_vandermonde_selected_solver == qmomkit::VandermondeSolver::bjorck_pereyra);
            assert(diagnostics.conditional_vandermonde_algorithm_used == qmomkit::internal::VandermondeAlgorithm::bjorck_pereyra);
        }
    }
}

void test_recursive_dense_and_pivoted_qr_equivalence_fixed()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2}, {2, 0, 1});
    const std::vector<std::size_t> extents{4, 4, 4};
    const std::vector<double> weights{0.05, 0.15, 0.20, 0.10, 0.12, 0.18, 0.08, 0.12};
    const std::vector<std::vector<double>> abscissas{{-1.5, 0.25, -0.75, 1.75, -1.25, 0.75, -0.10, 2.00},
                                                     {-2.0, -1.0, 0.5, 1.5, -1.5, -0.5, 1.0, 2.0},
                                                     {-0.80, -0.80, -0.20, -0.20, 0.60, 0.60, 1.20, 1.20}};
    const std::vector<double> moments = raw_moments_from_terminal_quadrature(abscissas, weights, extents);

    const qmomkit::internal::RecursiveQuadratureMomentDataRequest dense_request =
        moment_data_request(hierarchy, {extents, moments}, {}, qmomkit::WheelerOptions{false}, qmomkit::CqmomDiagnosticLevel::detailed,
                            qmomkit::VandermondeSolver::dense);
    const qmomkit::internal::RecursiveQuadratureMomentDataRequest qr_request =
        moment_data_request(hierarchy, {extents, moments}, {}, qmomkit::WheelerOptions{false}, qmomkit::CqmomDiagnosticLevel::detailed,
                            qmomkit::VandermondeSolver::pivoted_qr);

    const qmomkit::internal::RecursiveQuadratureResult dense = qmomkit::internal::compute_recursive_quadrature_from_moment_data(dense_request);
    const qmomkit::internal::RecursiveQuadratureResult qr = qmomkit::internal::compute_recursive_quadrature_from_moment_data(qr_request);

    assert_recursive_results_equivalent(qr, dense, 1.0e-9);
    assert(qr.ok());
    assert(!qr.diagnostics.empty());
    for (const qmomkit::internal::RecursiveQuadratureBranchDiagnostics& diagnostics : qr.diagnostics)
    {
        if (diagnostics.conditioning_level + 1 < hierarchy.number_of_dimensions())
        {
            assert(diagnostics.conditional_vandermonde_requested_solver == qmomkit::VandermondeSolver::pivoted_qr);
            assert(diagnostics.conditional_vandermonde_selected_solver == qmomkit::VandermondeSolver::pivoted_qr);
            assert(diagnostics.conditional_vandermonde_algorithm_used == qmomkit::internal::VandermondeAlgorithm::pivoted_qr);
        }
    }
}

void test_recursive_dense_and_svd_equivalence_fixed()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2}, {2, 0, 1});
    const std::vector<std::size_t> extents{4, 4, 4};
    const std::vector<double> weights{0.05, 0.15, 0.20, 0.10, 0.12, 0.18, 0.08, 0.12};
    const std::vector<std::vector<double>> abscissas{{-1.5, 0.25, -0.75, 1.75, -1.25, 0.75, -0.10, 2.00},
                                                     {-2.0, -1.0, 0.5, 1.5, -1.5, -0.5, 1.0, 2.0},
                                                     {-0.80, -0.80, -0.20, -0.20, 0.60, 0.60, 1.20, 1.20}};
    const std::vector<double> moments = raw_moments_from_terminal_quadrature(abscissas, weights, extents);

    const qmomkit::internal::RecursiveQuadratureMomentDataRequest dense_request =
        moment_data_request(hierarchy, {extents, moments}, {}, qmomkit::WheelerOptions{false}, qmomkit::CqmomDiagnosticLevel::detailed,
                            qmomkit::VandermondeSolver::dense);
    const qmomkit::internal::RecursiveQuadratureMomentDataRequest svd_request = moment_data_request(
        hierarchy, {extents, moments}, {}, qmomkit::WheelerOptions{false}, qmomkit::CqmomDiagnosticLevel::detailed, qmomkit::VandermondeSolver::svd);

    const qmomkit::internal::RecursiveQuadratureResult dense = qmomkit::internal::compute_recursive_quadrature_from_moment_data(dense_request);
    const qmomkit::internal::RecursiveQuadratureResult svd = qmomkit::internal::compute_recursive_quadrature_from_moment_data(svd_request);

    assert_recursive_results_equivalent(svd, dense, 1.0e-9);
    assert(svd.ok());
    assert(!svd.diagnostics.empty());
    for (const qmomkit::internal::RecursiveQuadratureBranchDiagnostics& diagnostics : svd.diagnostics)
    {
        if (diagnostics.conditioning_level + 1 < hierarchy.number_of_dimensions())
        {
            assert(diagnostics.conditional_vandermonde_requested_solver == qmomkit::VandermondeSolver::svd);
            assert(diagnostics.conditional_vandermonde_selected_solver == qmomkit::VandermondeSolver::svd);
            assert(diagnostics.conditional_vandermonde_algorithm_used == qmomkit::internal::VandermondeAlgorithm::svd);
        }
    }
}

void test_recursive_dense_and_moore_penrose_equivalence_fixed()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2}, {2, 0, 1});
    const std::vector<std::size_t> extents{4, 4, 4};
    const std::vector<double> weights{0.05, 0.15, 0.20, 0.10, 0.12, 0.18, 0.08, 0.12};
    const std::vector<std::vector<double>> abscissas{{-1.5, 0.25, -0.75, 1.75, -1.25, 0.75, -0.10, 2.00},
                                                     {-2.0, -1.0, 0.5, 1.5, -1.5, -0.5, 1.0, 2.0},
                                                     {-0.80, -0.80, -0.20, -0.20, 0.60, 0.60, 1.20, 1.20}};
    const std::vector<double> moments = raw_moments_from_terminal_quadrature(abscissas, weights, extents);

    const qmomkit::internal::RecursiveQuadratureMomentDataRequest dense_request =
        moment_data_request(hierarchy, {extents, moments}, {}, qmomkit::WheelerOptions{false}, qmomkit::CqmomDiagnosticLevel::detailed,
                            qmomkit::VandermondeSolver::dense);
    qmomkit::ConditionalMomentOptions conditional_options;
    conditional_options.moore_penrose_rcond = 1.0e-14;
    const qmomkit::internal::RecursiveQuadratureMomentDataRequest moore_penrose_request =
        moment_data_request(hierarchy, {extents, moments}, conditional_options, qmomkit::WheelerOptions{false},
                            qmomkit::CqmomDiagnosticLevel::detailed, qmomkit::VandermondeSolver::moore_penrose);

    const qmomkit::internal::RecursiveQuadratureResult dense = qmomkit::internal::compute_recursive_quadrature_from_moment_data(dense_request);
    const qmomkit::internal::RecursiveQuadratureResult moore_penrose =
        qmomkit::internal::compute_recursive_quadrature_from_moment_data(moore_penrose_request);

    assert_recursive_results_equivalent(moore_penrose, dense, 1.0e-9);
    assert(moore_penrose.ok());
    for (const qmomkit::internal::RecursiveQuadratureBranchDiagnostics& diagnostics : moore_penrose.diagnostics)
    {
        if (diagnostics.conditioning_level + 1 < hierarchy.number_of_dimensions())
        {
            assert(diagnostics.conditional_vandermonde_requested_solver == qmomkit::VandermondeSolver::moore_penrose);
            assert(diagnostics.conditional_vandermonde_selected_solver == qmomkit::VandermondeSolver::moore_penrose);
            assert(diagnostics.conditional_vandermonde_algorithm_used == qmomkit::internal::VandermondeAlgorithm::moore_penrose);
        }
    }
}

void test_recursive_moore_penrose_rcond_changes_conditional_reconstruction()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 1});
    const std::vector<std::size_t> extents{4, 2};
    const std::vector<std::vector<double>> abscissas{{0.0, 1.0}, {0.0, 2.0}};
    const std::vector<double> weights{1.0, 1.0};
    const std::vector<double> moments = raw_moments_from_terminal_quadrature(abscissas, weights, extents);

    qmomkit::ConditionalMomentOptions low_options;
    low_options.moore_penrose_rcond = 0.1;
    qmomkit::ConditionalMomentOptions high_options = low_options;
    high_options.moore_penrose_rcond = 0.5;
    const qmomkit::internal::RecursiveQuadratureMomentDataRequest low_request =
        moment_data_request(hierarchy, {extents, moments}, low_options, qmomkit::WheelerOptions{false}, qmomkit::CqmomDiagnosticLevel::detailed,
                            qmomkit::VandermondeSolver::moore_penrose);
    const qmomkit::internal::RecursiveQuadratureMomentDataRequest high_request =
        moment_data_request(hierarchy, {extents, moments}, high_options, qmomkit::WheelerOptions{false}, qmomkit::CqmomDiagnosticLevel::detailed,
                            qmomkit::VandermondeSolver::moore_penrose);

    const qmomkit::internal::RecursiveQuadratureResult low = qmomkit::internal::compute_recursive_quadrature_from_moment_data(low_request);
    const qmomkit::internal::RecursiveQuadratureResult high = qmomkit::internal::compute_recursive_quadrature_from_moment_data(high_request);

    // The root conditional system is solved with the parent nodes divided by
    // s = 2, so the truncated system is V = [[1, 1], [0, 1/2]] with
    // right-hand side [[2, 2], [1/2, 1]]. Its singular-value ratio is
    // (9 - sqrt(65)) / 4, so rcond = 0.1 retains both directions and
    // rcond = 0.5 retains only the larger direction. That direction is (1, t)
    // with t = (1 + sqrt(65)) / 8 and squared singular value 1 + t; projecting
    // V^T times the right-hand side onto it gives the retained solution.
    const double t = (1.0 + std::sqrt(65.0)) / 8.0;
    const double retained_weight = (2.0 + 2.25 * t) / ((1.0 + t) * (1.0 + t * t));
    const double retained_abscissa = (2.0 + 2.5 * t) / (2.0 + 2.25 * t);

    assert(low.ok());
    assert(high.ok());
    assert_near(low.quadrature.weight({0, 0}), 1.0, 1.0e-12);
    assert_near(low.quadrature.weight({1, 0}), 1.0, 1.0e-12);
    assert_near(high.quadrature.weight({0, 0}), retained_weight, 1.0e-12);
    assert_near(high.quadrature.weight({1, 0}), t * retained_weight, 1.0e-12);
    assert_near(high.quadrature.abscissa(1, {0, 0}), retained_abscissa, 1.0e-12);
    assert_near(high.quadrature.abscissa(1, {1, 0}), retained_abscissa, 1.0e-12);

    assert(!low.diagnostics.empty());
    assert(!high.diagnostics.empty());
    assert(low.diagnostics.front().conditional_vandermonde_requested_solver == qmomkit::VandermondeSolver::moore_penrose);
    assert(high.diagnostics.front().conditional_vandermonde_requested_solver == qmomkit::VandermondeSolver::moore_penrose);
}

void test_coordinate_scales_fall_back_from_rms_to_mean_to_one()
{
    // With two nodes per coordinate the second marginal moments are required,
    // so each scale is the RMS sqrt(m2 / m0) of its coordinate.
    const qmomkit::PhaseSpaceHierarchy two_nodes({2, 2});
    const std::vector<std::size_t> two_node_extents{4, 4};
    const std::vector<double> rms_scales = qmomkit::internal::cqmom_coordinate_scales(qmomkit::extract_cqmom_moments(
        {two_node_extents, raw_moments_from_terminal_quadrature({{1.0, 3.0}, {-2.0, 2.0}}, {0.5, 0.5}, two_node_extents)},
        qmomkit::create_cqmom_moment_layout(two_nodes)));
    assert(rms_scales.size() == 2);
    assert_near(rms_scales[0], std::sqrt(5.0), 1.0e-15);
    assert_near(rms_scales[1], 2.0, 1.0e-15);

    // One node per coordinate requires no second moment, so the scale falls
    // back to |m1| / m0, and to one where the mean is zero.
    const qmomkit::PhaseSpaceHierarchy one_node({1, 1});
    const std::vector<std::size_t> one_node_extents{2, 2};
    const std::vector<double> fallback_scales = qmomkit::internal::cqmom_coordinate_scales(
        qmomkit::extract_cqmom_moments({one_node_extents, raw_moments_from_terminal_quadrature({{-3.0}, {0.0}}, {2.0}, one_node_extents)},
                                       qmomkit::create_cqmom_moment_layout(one_node)));
    assert(fallback_scales.size() == 2);
    assert_near(fallback_scales[0], 3.0, 1.0e-15);
    assert(fallback_scales[1] == 1.0);
}

void test_recursive_bp_success_independent_of_diagnostic_level()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2}, {2, 0, 1});
    const std::vector<std::size_t> extents{4, 4, 4};
    const std::vector<double> weights{0.05, 0.15, 0.20, 0.10, 0.12, 0.18, 0.08, 0.12};
    const std::vector<std::vector<double>> abscissas{{-1.5, 0.25, -0.75, 1.75, -1.25, 0.75, -0.10, 2.00},
                                                     {-2.0, -1.0, 0.5, 1.5, -1.5, -0.5, 1.0, 2.0},
                                                     {-0.80, -0.80, -0.20, -0.20, 0.60, 0.60, 1.20, 1.20}};
    const std::vector<double> moments = raw_moments_from_terminal_quadrature(abscissas, weights, extents);

    std::vector<qmomkit::internal::RecursiveQuadratureResult> results;
    for (const qmomkit::CqmomDiagnosticLevel level : diagnostic_levels())
    {
        results.push_back(qmomkit::internal::compute_recursive_quadrature_from_moment_data(moment_data_request(
            hierarchy, {extents, moments}, {}, qmomkit::WheelerOptions{false}, level, qmomkit::VandermondeSolver::bjorck_pereyra)));
    }

    assert(results.front().ok());
    for (const qmomkit::internal::RecursiveQuadratureResult& result : results)
    {
        assert_recursive_numerics_equivalent(result, results.front(), 1.0e-9);
    }
    assert(results[0].diagnostics.empty());
    for (std::size_t i = 1; i < results.size(); ++i)
    {
        assert(!results[i].diagnostics.empty());
        for (const qmomkit::internal::RecursiveQuadratureBranchDiagnostics& diagnostics : results[i].diagnostics)
        {
            if (diagnostics.conditioning_level + 1 < hierarchy.number_of_dimensions())
            {
                assert(diagnostics.conditional_vandermonde_algorithm_used == qmomkit::internal::VandermondeAlgorithm::bjorck_pereyra);
            }
        }
    }
}

void test_recursive_automatic_uses_plain_bp_for_small_branches()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {1, 0});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.08, 0.12, 0.56, 0.24};
    const std::vector<std::vector<double>> abscissas{{-2.0, 1.0, -0.5, 2.25}, {-1.5, -1.5, 1.0, 1.0}};
    const std::vector<double> moments = raw_moments_from_terminal_quadrature(abscissas, weights, extents);

    const qmomkit::internal::RecursiveQuadratureMomentDataRequest dense_request =
        moment_data_request(hierarchy, {extents, moments}, {}, qmomkit::WheelerOptions{false}, qmomkit::CqmomDiagnosticLevel::detailed,
                            qmomkit::VandermondeSolver::dense);
    const qmomkit::internal::RecursiveQuadratureMomentDataRequest automatic_request =
        moment_data_request(hierarchy, {extents, moments}, {}, qmomkit::WheelerOptions{false}, qmomkit::CqmomDiagnosticLevel::detailed,
                            qmomkit::VandermondeSolver::automatic);

    const qmomkit::internal::RecursiveQuadratureResult dense = qmomkit::internal::compute_recursive_quadrature_from_moment_data(dense_request);
    const qmomkit::internal::RecursiveQuadratureResult automatic =
        qmomkit::internal::compute_recursive_quadrature_from_moment_data(automatic_request);

    assert_recursive_results_equivalent(automatic, dense, 1.0e-12);
    assert(automatic.ok());
    assert(!automatic.diagnostics.empty());
    for (const qmomkit::internal::RecursiveQuadratureBranchDiagnostics& diagnostics : automatic.diagnostics)
    {
        if (diagnostics.conditioning_level + 1 < hierarchy.number_of_dimensions())
        {
            assert(diagnostics.conditional_vandermonde_requested_solver == qmomkit::VandermondeSolver::automatic);
            assert(diagnostics.conditional_vandermonde_selected_solver == qmomkit::VandermondeSolver::bjorck_pereyra);
            assert(diagnostics.conditional_vandermonde_algorithm_used == qmomkit::internal::VandermondeAlgorithm::bjorck_pereyra);
        }
    }
}

void test_recursive_automatic_independent_of_diagnostic_level()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {1, 0});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.08, 0.12, 0.56, 0.24};
    const std::vector<std::vector<double>> abscissas{{-2.0, 1.0, -0.5, 2.25}, {-1.5, -1.5, 1.0, 1.0}};
    const std::vector<double> moments = raw_moments_from_terminal_quadrature(abscissas, weights, extents);

    std::vector<qmomkit::internal::RecursiveQuadratureResult> results;
    for (const qmomkit::CqmomDiagnosticLevel level : diagnostic_levels())
    {
        results.push_back(qmomkit::internal::compute_recursive_quadrature_from_moment_data(
            moment_data_request(hierarchy, {extents, moments}, {}, qmomkit::WheelerOptions{false}, level, qmomkit::VandermondeSolver::automatic)));
    }

    assert(results.front().ok());
    for (const qmomkit::internal::RecursiveQuadratureResult& result : results)
    {
        assert_recursive_numerics_equivalent(result, results.front(), 1.0e-12);
    }
}

void test_recursive_dense_and_bp_equivalence_adaptive()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2});
    const std::vector<std::size_t> extents{4, 4, 4};
    const std::vector<double> weights{0.05, 0.05, 0.30, 0.60};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, -1.0, 1.0}, {-2.0, -2.0, 2.0, 3.0}, {-1.0, 1.0, 0.0, 2.0}};
    const std::vector<double> moments = raw_moments_from_terminal_quadrature(abscissas, weights, extents);
    const qmomkit::WheelerOptions wheeler_options{true, 1.0e-10, 1.0e-10, 1.0e-12};

    const qmomkit::internal::RecursiveQuadratureMomentDataRequest dense_request = moment_data_request(
        hierarchy, {extents, moments}, {}, wheeler_options, qmomkit::CqmomDiagnosticLevel::basic, qmomkit::VandermondeSolver::dense);
    const qmomkit::internal::RecursiveQuadratureMomentDataRequest bp_request = moment_data_request(
        hierarchy, {extents, moments}, {}, wheeler_options, qmomkit::CqmomDiagnosticLevel::basic, qmomkit::VandermondeSolver::bjorck_pereyra);

    const qmomkit::internal::RecursiveQuadratureResult dense = qmomkit::internal::compute_recursive_quadrature_from_moment_data(dense_request);
    const qmomkit::internal::RecursiveQuadratureResult bp = qmomkit::internal::compute_recursive_quadrature_from_moment_data(bp_request);

    assert_recursive_results_equivalent(bp, dense, 1.0e-9);
    assert(bp.ok());
    assert_ragged_tree_consistent(bp.quadrature);
}

void test_recursive_adaptive_reduces_roundoff_deterministic_mixed_affine_branch()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2, 2}, {0, 2, 1, 3});
    const std::vector<std::size_t> extents{4, 4, 4, 4};
    const std::vector<double> weights{0.45 * 0.3, 0.45 * 0.7, 0.55 * 0.3, 0.55 * 0.7};
    const std::vector<std::vector<double>> abscissas{{-0.5, -0.5, 1.0, 1.0}, {0.0, 0.0, 3.0, 3.0}, {-1.0, 2.0, -1.0, 2.0}, {-0.5, 1.0, -3.5, -2.0}};
    const std::vector<double> moments = raw_moments_from_terminal_quadrature(abscissas, weights, extents);
    const qmomkit::WheelerOptions wheeler_options{true, 1.0e-10, 1.0e-10, 1.0e-30};

    const qmomkit::internal::RecursiveQuadratureMomentDataRequest automatic_request = moment_data_request(
        hierarchy, {extents, moments}, {}, wheeler_options, qmomkit::CqmomDiagnosticLevel::detailed, qmomkit::VandermondeSolver::automatic);
    const qmomkit::internal::RecursiveQuadratureMomentDataRequest dense_request = moment_data_request(
        hierarchy, {extents, moments}, {}, wheeler_options, qmomkit::CqmomDiagnosticLevel::detailed, qmomkit::VandermondeSolver::dense);

    const qmomkit::internal::RecursiveQuadratureResult automatic =
        qmomkit::internal::compute_recursive_quadrature_from_moment_data(automatic_request);
    const qmomkit::internal::RecursiveQuadratureResult dense = qmomkit::internal::compute_recursive_quadrature_from_moment_data(dense_request);

    assert_recursive_results_equivalent(automatic, dense, 1.0e-9);
    assert(automatic.ok());
    assert_ragged_tree_consistent(automatic.quadrature);

    bool found_branch = false;
    for (const qmomkit::internal::RecursiveQuadratureBranchDiagnostics& diagnostics : automatic.diagnostics)
    {
        if (diagnostics.conditioning_level == 2 && diagnostics.path == qmomkit::ConditioningPath({0, 0}))
        {
            found_branch = true;
            assert(diagnostics.phase_space_dimension == 1);
            assert(diagnostics.requested_node_count == 2);
            assert(diagnostics.actual_node_count == 1);
            assert(diagnostics.adaptive_node_count_reduced);
            assert(diagnostics.wheeler_status == qmomkit::QuadratureStatus::success);
        }
    }
    assert(found_branch);
}

void test_recursive_adaptive_reduces_roundoff_branch_with_bounded_support()
{
    // The same rule with non-binding bounds on every coordinate, so that the
    // bounded selection inverts the conditional moments. The branch whose
    // conditional distribution is a single point must still reduce despite
    // the rounding the conditional solve leaves in its moments.
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2, 2}, {0, 2, 1, 3});
    const std::vector<std::size_t> extents{4, 4, 4, 4};
    const std::vector<double> weights{0.45 * 0.3, 0.45 * 0.7, 0.55 * 0.3, 0.55 * 0.7};
    const std::vector<std::vector<double>> abscissas{{-0.5, -0.5, 1.0, 1.0}, {0.0, 0.0, 3.0, 3.0}, {-1.0, 2.0, -1.0, 2.0}, {-0.5, 1.0, -3.5, -2.0}};
    const std::vector<double> moments = raw_moments_from_terminal_quadrature(abscissas, weights, extents);
    const qmomkit::WheelerOptions wheeler_options{true, 1.0e-10, 1.0e-10, 1.0e-30};
    const qmomkit::CqmomAbscissaBounds support(4, qmomkit::AbscissaBounds{-10.0, 10.0});

    const qmomkit::internal::RecursiveQuadratureResult result = qmomkit::internal::compute_recursive_quadrature_from_moment_data(moment_data_request(
        hierarchy, {extents, moments}, {}, wheeler_options, qmomkit::CqmomDiagnosticLevel::detailed, qmomkit::VandermondeSolver::automatic, support));
    assert(result.ok());

    bool found_branch = false;
    for (const qmomkit::internal::RecursiveQuadratureBranchDiagnostics& diagnostics : result.diagnostics)
    {
        if (diagnostics.conditioning_level == 2 && diagnostics.path == qmomkit::ConditioningPath({0, 0}))
        {
            found_branch = true;
            assert(diagnostics.bounded_support_checked);
            assert(diagnostics.actual_node_count == 1);
        }
    }
    assert(found_branch);
}

// An adaptive request for a 3D rule of six weighted points whose second branch
// holds a single y point, so that branch reduces to one node.
qmomkit::internal::RecursiveQuadratureMomentDataRequest single_point_branch_3d_request()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2});
    const std::vector<std::size_t> extents{4, 4, 4};
    const std::vector<double> weights{0.4 * 0.5 * 0.3, 0.4 * 0.5 * 0.7, 0.4 * 0.5 * 0.3, 0.4 * 0.5 * 0.7, 0.6 * 0.3, 0.6 * 0.7};
    const std::vector<std::vector<double>> abscissas{
        {-1.0, -1.0, -1.0, -1.0, 1.0, 1.0}, {-1.0, -1.0, 1.0, 1.0, 0.0, 0.0}, {-0.5, 0.5, -0.5, 0.5, -0.5, 0.5}};
    return moment_data_request(hierarchy, {extents, raw_moments_from_terminal_quadrature(abscissas, weights, extents)}, {},
                               qmomkit::WheelerOptions{true, 1.0e-10, 1.0e-10, 1.0e-12});
}

void test_recursive_adaptive_records_active_counts_in_parent_path_order()
{
    // The traversal visits branches level by level, but the recorded counts
    // list each branch before its descendants, the order in which the terminal
    // paths visit them.
    qmomkit::internal::RecursiveQuadratureMomentDataRequest request = single_point_branch_3d_request();
    const qmomkit::internal::RecursiveQuadratureResult unrecorded = qmomkit::internal::compute_recursive_quadrature_from_moment_data(request);
    assert(unrecorded.ok());
    assert(unrecorded.active_branch_counts.empty());

    request.record_active_branch_counts = true;
    const qmomkit::internal::RecursiveQuadratureResult result = qmomkit::internal::compute_recursive_quadrature_from_moment_data(request);
    assert(result.ok());
    assert(result.quadrature.total_nodes() == 6);
    const std::vector<qmomkit::CqmomActiveBranchCount> expected{
        {0, qmomkit::ConditioningPath{}, 0, 2, 2},       {1, qmomkit::ConditioningPath({0}), 1, 2, 2},
        {2, qmomkit::ConditioningPath({0, 0}), 2, 2, 2}, {2, qmomkit::ConditioningPath({0, 1}), 2, 2, 2},
        {1, qmomkit::ConditioningPath({1}), 1, 2, 1},    {2, qmomkit::ConditioningPath({1, 0}), 2, 2, 2}};
    assert(result.active_branch_counts.size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i)
    {
        assert(result.active_branch_counts[i].conditioning_level == expected[i].conditioning_level);
        assert(result.active_branch_counts[i].parent_path == expected[i].parent_path);
        assert(result.active_branch_counts[i].phase_space_dimension == expected[i].phase_space_dimension);
        assert(result.active_branch_counts[i].configured_node_count == expected[i].configured_node_count);
        assert(result.active_branch_counts[i].active_node_count == expected[i].active_node_count);
    }

    // A fixed reconstruction has no active counts to record.
    const std::vector<double> product_weights{0.06, 0.14, 0.06, 0.14, 0.09, 0.21, 0.09, 0.21};
    const std::vector<std::vector<double>> product_abscissas{
        {-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0}, {-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0}, {-0.5, 0.5, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5}};
    const std::vector<std::size_t> extents{4, 4, 4};
    qmomkit::internal::RecursiveQuadratureMomentDataRequest fixed_request =
        moment_data_request(request.hierarchy, {extents, raw_moments_from_terminal_quadrature(product_abscissas, product_weights, extents)});
    fixed_request.record_active_branch_counts = true;
    const qmomkit::internal::RecursiveQuadratureResult fixed = qmomkit::internal::compute_recursive_quadrature_from_moment_data(fixed_request);
    assert(fixed.ok());
    assert(fixed.active_branch_counts.empty());
}

void test_recursive_adaptive_validates_supplied_active_counts()
{
    // The rule constrained by the active counts its own reconstruction
    // records. Counts in any order describe the same topology; counts that do
    // not describe one record per reachable branch are rejected.
    qmomkit::internal::RecursiveQuadratureMomentDataRequest request = single_point_branch_3d_request();
    request.record_active_branch_counts = true;
    const qmomkit::internal::RecursiveQuadratureResult recorded = qmomkit::internal::compute_recursive_quadrature_from_moment_data(request);
    assert(recorded.ok());
    const std::vector<qmomkit::CqmomActiveBranchCount> counts = recorded.active_branch_counts;
    assert(counts.size() == 6);

    const auto reconstruct_with = [&request](std::vector<qmomkit::CqmomActiveBranchCount> supplied)
    {
        qmomkit::internal::RecursiveQuadratureMomentDataRequest constrained = request;
        constrained.active_branch_counts = std::move(supplied);
        return qmomkit::internal::compute_recursive_quadrature_from_moment_data(constrained);
    };
    const auto assert_rejected = [&reconstruct_with](std::vector<qmomkit::CqmomActiveBranchCount> supplied, const std::string& message)
    {
        const qmomkit::internal::RecursiveQuadratureResult result = reconstruct_with(std::move(supplied));
        assert(result.status == qmomkit::internal::RecursiveQuadratureStatus::invalid_input);
        assert(result.message.find(message) != std::string::npos);
    };

    const qmomkit::internal::RecursiveQuadratureResult reversed = reconstruct_with({counts.rbegin(), counts.rend()});
    assert(reversed.ok());
    assert(reversed.quadrature.total_nodes() == recorded.quadrature.total_nodes());
    assert(reversed.active_branch_counts.size() == counts.size());
    for (std::size_t i = 0; i < counts.size(); ++i)
    {
        assert(reversed.active_branch_counts[i].parent_path == counts[i].parent_path);
        assert(reversed.active_branch_counts[i].active_node_count == counts[i].active_node_count);
    }

    // counts[4] is the single-node branch [1], so [1, 1] is not a branch.
    std::vector<qmomkit::CqmomActiveBranchCount> stale = counts;
    stale.push_back({2, qmomkit::ConditioningPath({1, 1}), 2, 2, 2});
    assert_rejected(stale, "unreachable or stale branches");

    std::vector<qmomkit::CqmomActiveBranchCount> missing = counts;
    missing.erase(missing.begin() + 3);
    assert_rejected(missing, "missing a reachable branch");

    // Without its parent [1], the record for [1, 0] is unreachable as well; the
    // missing branch is reported.
    std::vector<qmomkit::CqmomActiveBranchCount> orphan = counts;
    orphan.erase(orphan.begin() + 4);
    assert_rejected(orphan, "missing a reachable branch");

    std::vector<qmomkit::CqmomActiveBranchCount> duplicate = counts;
    duplicate.push_back(counts[1]);
    assert_rejected(duplicate, "duplicate entries for a branch");

    std::vector<qmomkit::CqmomActiveBranchCount> over_capacity = counts;
    over_capacity[1].active_node_count = 3;
    assert_rejected(over_capacity, "outside the configured node-count capacity");
}

void test_recursive_reconstruction_does_not_depend_on_number_density()
{
    // The same four weighted points at unit number density and at 1e-16: the
    // nodes must agree and the weights must scale with the number density, in
    // fixed and in adaptive mode. The branch-local weighted Vandermonde solve
    // sees parent weights of 5e-17 at the lower density.
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.0, 0.0, 3.0}};
    const std::vector<double> weights{0.2, 0.3, 0.1, 0.4};
    constexpr double number_density = 1.0e-16;
    std::vector<double> dilute_weights = weights;
    for (double& weight : dilute_weights)
    {
        weight *= number_density;
    }

    for (const bool adaptive : {false, true})
    {
        const qmomkit::WheelerOptions wheeler_options{adaptive, 1.0e-10, 1.0e-10, 1.0e-30};
        const qmomkit::internal::RecursiveQuadratureResult reference = qmomkit::internal::compute_recursive_quadrature_from_moment_data(
            moment_data_request(hierarchy, {extents, raw_moments_from_terminal_quadrature(abscissas, weights, extents)}, {}, wheeler_options));
        const qmomkit::internal::RecursiveQuadratureResult dilute = qmomkit::internal::compute_recursive_quadrature_from_moment_data(
            moment_data_request(hierarchy, {extents, raw_moments_from_terminal_quadrature(abscissas, dilute_weights, extents)}, {}, wheeler_options));
        assert(reference.ok());
        assert(dilute.ok());
        assert(dilute.quadrature.weights().size() == weights.size());
        assert(dilute.quadrature.terminal_paths() == reference.quadrature.terminal_paths());
        for (std::size_t node = 0; node < weights.size(); ++node)
        {
            assert_near(dilute.quadrature.weights()[node] / number_density, reference.quadrature.weights()[node], 1.0e-12);
            for (std::size_t dimension = 0; dimension < abscissas.size(); ++dimension)
            {
                assert_near(dilute.quadrature.abscissas()[dimension][node], reference.quadrature.abscissas()[dimension][node], 1.0e-12);
            }
        }
    }
}

void test_recursive_bp_failure_independent_of_diagnostic_level()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<double> parent_abscissas{-1.0, 1.0};
    const std::vector<double> parent_weights{0.5, 0.5};
    const std::vector<std::size_t> extents{4, 4};
    std::vector<double> moments(qmomkit::product_extents(extents), 0.0);
    const std::vector<std::vector<double>> conditional_moments{{1.0, 0.0, 1.0, 0.0}, {1.0, 0.0, -1.0, 0.0}};
    for (std::size_t flat_index = 0; flat_index < moments.size(); ++flat_index)
    {
        const qmomkit::MultiIndex moment_index = qmomkit::unflatten_multi_index(flat_index, extents);
        for (std::size_t parent_node = 0; parent_node < parent_abscissas.size(); ++parent_node)
        {
            moments[flat_index] += parent_weights[parent_node] * std::pow(parent_abscissas[parent_node], static_cast<int>(moment_index[0])) *
                                   conditional_moments[parent_node][moment_index[1]];
        }
    }

    std::vector<qmomkit::internal::RecursiveQuadratureResult> results;
    for (const qmomkit::CqmomDiagnosticLevel level : diagnostic_levels())
    {
        results.push_back(qmomkit::internal::compute_recursive_quadrature_from_moment_data(moment_data_request(
            hierarchy, {extents, moments}, {}, qmomkit::WheelerOptions{false}, level, qmomkit::VandermondeSolver::bjorck_pereyra)));
    }

    assert(!results.front().ok());
    for (const qmomkit::internal::RecursiveQuadratureResult& result : results)
    {
        assert_recursive_numerics_equivalent(result, results.front(), 0.0);
        assert(result.status == qmomkit::internal::RecursiveQuadratureStatus::wheeler_failed);
        assert(result.failed_conditioning_level == 1);
    }
}

// Diagnostics are observational: raising the diagnostic level may compute extra
// residuals and condition estimates, but the accepted quadrature must be
// identical, not merely close. This is checked on a three-dimensional
// non-identity conditioning order with the Bjorck-Pereyra solver, where the
// extra diagnostic work touches every conditional branch.
void test_public_cqmom_bp_success_independent_of_diagnostic_level()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2}, {2, 0, 1});
    const std::vector<std::size_t> extents{4, 4, 4};
    const std::vector<double> weights{0.05, 0.15, 0.20, 0.10, 0.12, 0.18, 0.08, 0.12};
    const std::vector<std::vector<double>> abscissas{{-1.5, 0.25, -0.75, 1.75, -1.25, 0.75, -0.10, 2.00},
                                                     {-2.0, -1.0, 0.5, 1.5, -1.5, -0.5, 1.0, 2.0},
                                                     {-0.80, -0.80, -0.20, -0.20, 0.60, 0.60, 1.20, 1.20}};
    const qmomkit::RawMomentTensor moments{extents, raw_moments_from_terminal_quadrature(abscissas, weights, extents)};
    const qmomkit::CqmomMomentData required_moments = qmomkit::extract_cqmom_moments(moments, qmomkit::create_cqmom_moment_layout(hierarchy));

    std::vector<qmomkit::CqmomResult> results;
    for (const qmomkit::CqmomDiagnosticLevel level : diagnostic_levels())
    {
        results.push_back(qmomkit::compute_cqmom(
            qmomkit::CqmomCase{required_moments, {}, qmomkit::QmomOptions{false}, level, qmomkit::VandermondeSolver::bjorck_pereyra}));
    }

    assert(results.front().ok());
    for (const qmomkit::CqmomResult& result : results)
    {
        assert_cqmom_results_numerically_equivalent(result, results.front(), 0.0);
    }
}

// A flat multi-index table stores equal-rank indices contiguously, so copying
// one costs a single allocation instead of one per row. These tests pin the
// contract that moment layouts and terminal-path tables will rely on.
void test_multi_index_table_round_trips_and_views_rows()
{
    const std::vector<qmomkit::MultiIndex> rows{{0, 1, 2}, {3, 4, 5}, {6, 7, 8}};
    const qmomkit::MultiIndexTable table(rows);

    assert(table.size() == 3);
    assert(table.rank() == 3);
    assert(!table.empty());
    assert(table.values() == std::vector<std::size_t>({0, 1, 2, 3, 4, 5, 6, 7, 8}));

    for (std::size_t row = 0; row < rows.size(); ++row)
    {
        const qmomkit::MultiIndexView view = table[row];
        assert(view.size() == 3);
        assert(view == rows[row]);
        assert(rows[row] == view);
        assert(view.to_multi_index() == rows[row]);
        for (std::size_t axis = 0; axis < view.size(); ++axis)
        {
            assert(view[axis] == rows[row][axis]);
        }
    }
    assert(table.to_multi_indices() == rows);
}

// Contiguity is the point of the type: rows are adjacent in memory, which is
// what removes the per-row allocation and keeps traversal sequential.
void test_multi_index_table_storage_is_contiguous_and_row_major()
{
    const qmomkit::MultiIndexTable table(std::vector<qmomkit::MultiIndex>{{1, 2}, {3, 4}});
    assert(table[1].data() == table[0].data() + table.rank());
    assert(table.values().data() == table[0].data());
}

void test_multi_index_table_accepts_flat_storage_and_rejects_ragged_input()
{
    const qmomkit::MultiIndexTable flat(2, std::vector<std::size_t>{1, 2, 3, 4});
    assert(flat.size() == 2);
    assert(flat[1] == qmomkit::MultiIndex({3, 4}));

    assert_throws_invalid_argument([] { (void)qmomkit::MultiIndexTable(2, std::vector<std::size_t>{1, 2, 3}); });
    assert_throws_invalid_argument([] { (void)qmomkit::MultiIndexTable(0, std::vector<std::size_t>{1}); });
    assert_throws_invalid_argument([] { (void)qmomkit::MultiIndexTable(std::vector<qmomkit::MultiIndex>{{1, 2}, {3}}); });

    const qmomkit::MultiIndexTable table(std::vector<qmomkit::MultiIndex>{{1, 2}});
    bool threw = false;
    try
    {
        (void)table.at(1);
    }
    catch (const std::out_of_range&)
    {
        threw = true;
    }
    assert(threw);
}

void test_multi_index_table_grows_and_compares()
{
    qmomkit::MultiIndexTable table;
    assert(table.empty());
    assert(table.size() == 0);
    assert(table.rank() == 0);

    table.reserve(2);
    table.push_back(qmomkit::MultiIndex{5, 6});
    assert(table.rank() == 2);
    const qmomkit::MultiIndexTable source(std::vector<qmomkit::MultiIndex>{{7, 8}});
    table.push_back(source[0]);
    assert(table.size() == 2);
    assert(table[1] == qmomkit::MultiIndex({7, 8}));
    assert_throws_invalid_argument([&table] { table.push_back(qmomkit::MultiIndex{1, 2, 3}); });

    assert(table == qmomkit::MultiIndexTable(std::vector<qmomkit::MultiIndex>{{5, 6}, {7, 8}}));
    assert(table != qmomkit::MultiIndexTable(std::vector<qmomkit::MultiIndex>{{5, 6}}));
    // Equal flat storage with a different rank is a different table.
    assert(qmomkit::MultiIndexTable(1, std::vector<std::size_t>{1, 2}) != qmomkit::MultiIndexTable(2, std::vector<std::size_t>{1, 2}));

    table.clear();
    assert(table.empty());
}

// Range-based for and std::find are how moment layouts are traversed, so the
// row iterator must support both.
void test_multi_index_table_iterates_and_finds_rows()
{
    const std::vector<qmomkit::MultiIndex> rows{{1, 2}, {3, 4}, {5, 6}};
    const qmomkit::MultiIndexTable table(rows);

    std::size_t visited = 0;
    for (const qmomkit::MultiIndexView row : table)
    {
        assert(row == rows[visited]);
        ++visited;
    }
    assert(visited == rows.size());

    const auto found = std::find(table.begin(), table.end(), qmomkit::MultiIndex{3, 4});
    assert(found != table.end());
    assert(static_cast<std::size_t>(std::distance(table.begin(), found)) == 1);
    assert(std::find(table.begin(), table.end(), qmomkit::MultiIndex{9, 9}) == table.end());

    const qmomkit::MultiIndexTable empty_table;
    assert(empty_table.begin() == empty_table.end());
}

void test_multi_index_table_matches_enumerate_multi_indices()
{
    const std::vector<std::size_t> extents{2, 3};
    const std::vector<qmomkit::MultiIndex> enumerated = qmomkit::enumerate_multi_indices(extents);
    const qmomkit::MultiIndexTable table(enumerated);
    assert(table.size() == qmomkit::product_extents(extents));
    for (std::size_t row = 0; row < table.size(); ++row)
    {
        assert(table[row] == enumerated[row]);
        assert(qmomkit::flatten_multi_index(table[row].to_multi_index(), extents) == row);
    }
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"flat_multi_index_round_trip", test_flat_multi_index_round_trip},
            {"enumerate_multi_indices_row_major", test_enumerate_multi_indices_row_major},
            {"branch_local_layout_identity_wheeler_and_rhs_offsets", test_branch_local_layout_identity_wheeler_and_rhs_offsets},
            {"branch_local_layout_non_identity_wheeler_and_rhs_offsets", test_branch_local_layout_non_identity_wheeler_and_rhs_offsets},
            {"branch_local_layout_unequal_node_counts_and_adaptive_child_case", test_branch_local_layout_unequal_node_counts_and_adaptive_child_case},
            {"moment_data_recursive_rejects_missing_values", test_moment_data_recursive_rejects_missing_values},
            {"moment_data_recursive_rejects_non_finite_values", test_moment_data_recursive_rejects_non_finite_values},
            {"moment_data_recursive_rejects_incompatible_layout", test_moment_data_recursive_rejects_incompatible_layout},
            {"recursive_bounded_support_accepts_feasible_branch", test_recursive_bounded_support_accepts_feasible_branch},
            {"recursive_shifted_bounded_support_accepts_feasible_fixed_branch", test_recursive_shifted_bounded_support_accepts_feasible_fixed_branch},
            {"recursive_shifted_bounded_support_rejects_fixed_and_reduces_adaptive",
             test_recursive_shifted_bounded_support_rejects_fixed_and_reduces_adaptive},
            {"recursive_bounded_support_rejects_fixed_infeasible_count", test_recursive_bounded_support_rejects_fixed_infeasible_count},
            {"recursive_bounded_support_adaptive_falls_back_to_supported_count",
             test_recursive_bounded_support_adaptive_falls_back_to_supported_count},
            {"recursive_bounded_support_adaptive_two_node_uses_span_separation",
             test_recursive_bounded_support_adaptive_two_node_uses_span_separation},
            {"recursive_bounded_support_adaptive_applies_relative_weight_check",
             test_recursive_bounded_support_adaptive_applies_relative_weight_check},
            {"recursive_bounded_support_rejects_growth_and_preserves_current_count",
             test_recursive_bounded_support_rejects_growth_and_preserves_current_count},
            {"recursive_bounded_support_uses_phase_space_coordinate_index", test_recursive_bounded_support_uses_phase_space_coordinate_index},
            {"fixed_public_terminal_paths_align_with_physical_storage", test_fixed_public_terminal_paths_align_with_physical_storage},
            {"recursive_bounded_support_applies_to_conditional_phase_space_dimension",
             test_recursive_bounded_support_applies_to_conditional_phase_space_dimension},
            {"recursive_bounded_support_returned_node_tolerance_and_no_clipping",
             test_recursive_bounded_support_returned_node_tolerance_and_no_clipping},
            {"one_dimensional_hierarchy_reduces_to_stage_1_ordering", test_one_dimensional_hierarchy_reduces_to_stage_1_ordering},
            {"two_dimensional_hierarchy_indexing", test_two_dimensional_hierarchy_indexing},
            {"three_dimensional_hierarchy_indexing", test_three_dimensional_hierarchy_indexing},
            {"non_uniform_node_counts", test_non_uniform_node_counts},
            {"identity_conditioning_order_matches_phase_space_order", test_identity_conditioning_order_matches_phase_space_order},
            {"non_identity_conditioning_order_controls_recursive_traversal", test_non_identity_conditioning_order_controls_recursive_traversal},
            {"reversed_conditioning_order_controls_recursive_traversal", test_reversed_conditioning_order_controls_recursive_traversal},
            {"conditioning_path_table_round_trips_and_views_rows", test_conditioning_path_table_round_trips_and_views_rows},
            {"conditioning_path_table_rejects_uneven_depth_and_checks_bounds", test_conditioning_path_table_rejects_uneven_depth_and_checks_bounds},
            {"conditioning_path_table_grows_iterates_and_compares", test_conditioning_path_table_grows_iterates_and_compares},
            {"conditioning_path_view_maps_to_phase_space_node_index", test_conditioning_path_view_maps_to_phase_space_node_index},
            {"conditioning_path_maps_to_phase_space_node_index", test_conditioning_path_maps_to_phase_space_node_index},
            {"invalid_hierarchy_inputs", test_invalid_hierarchy_inputs},
            {"invalid_index_inputs", test_invalid_index_inputs},
            {"recursive_failure_diagnostics_identify_failing_branch", test_recursive_failure_diagnostics_identify_failing_branch},
            {"recursive_dense_and_bp_equivalence_fixed", test_recursive_dense_and_bp_equivalence_fixed},
            {"recursive_dense_and_pivoted_qr_equivalence_fixed", test_recursive_dense_and_pivoted_qr_equivalence_fixed},
            {"recursive_dense_and_svd_equivalence_fixed", test_recursive_dense_and_svd_equivalence_fixed},
            {"recursive_dense_and_moore_penrose_equivalence_fixed", test_recursive_dense_and_moore_penrose_equivalence_fixed},
            {"recursive_moore_penrose_rcond_changes_conditional_reconstruction",
             test_recursive_moore_penrose_rcond_changes_conditional_reconstruction},
            {"coordinate_scales_fall_back_from_rms_to_mean_to_one", test_coordinate_scales_fall_back_from_rms_to_mean_to_one},
            {"recursive_bp_success_independent_of_diagnostic_level", test_recursive_bp_success_independent_of_diagnostic_level},
            {"recursive_automatic_uses_plain_bp_for_small_branches", test_recursive_automatic_uses_plain_bp_for_small_branches},
            {"recursive_automatic_independent_of_diagnostic_level", test_recursive_automatic_independent_of_diagnostic_level},
            {"recursive_dense_and_bp_equivalence_adaptive", test_recursive_dense_and_bp_equivalence_adaptive},
            {"recursive_adaptive_reduces_roundoff_deterministic_mixed_affine_branch",
             test_recursive_adaptive_reduces_roundoff_deterministic_mixed_affine_branch},
            {"recursive_adaptive_reduces_roundoff_branch_with_bounded_support", test_recursive_adaptive_reduces_roundoff_branch_with_bounded_support},
            {"recursive_adaptive_records_active_counts_in_parent_path_order", test_recursive_adaptive_records_active_counts_in_parent_path_order},
            {"recursive_adaptive_validates_supplied_active_counts", test_recursive_adaptive_validates_supplied_active_counts},
            {"recursive_reconstruction_does_not_depend_on_number_density", test_recursive_reconstruction_does_not_depend_on_number_density},
            {"recursive_bp_failure_independent_of_diagnostic_level", test_recursive_bp_failure_independent_of_diagnostic_level},
            {"public_cqmom_bp_success_independent_of_diagnostic_level", test_public_cqmom_bp_success_independent_of_diagnostic_level},
            {"multi_index_table_round_trips_and_views_rows", test_multi_index_table_round_trips_and_views_rows},
            {"multi_index_table_storage_is_contiguous_and_row_major", test_multi_index_table_storage_is_contiguous_and_row_major},
            {"multi_index_table_accepts_flat_storage_and_rejects_ragged_input", test_multi_index_table_accepts_flat_storage_and_rejects_ragged_input},
            {"multi_index_table_grows_and_compares", test_multi_index_table_grows_and_compares},
            {"multi_index_table_iterates_and_finds_rows", test_multi_index_table_iterates_and_finds_rows},
            {"multi_index_table_matches_enumerate_multi_indices", test_multi_index_table_matches_enumerate_multi_indices},
        });
}
