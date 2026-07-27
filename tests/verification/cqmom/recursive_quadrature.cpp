#include "internal/recursive_quadrature.hpp"

#include "../verification_helpers.hpp"
#include "qmomkit/qmomkit.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{

void assert_vector_equal(const std::vector<std::size_t>& actual, const std::vector<std::size_t>& expected)
{
    QMOMKIT_VERIFY_CHECK(actual == expected);
}

void assert_vector_near(const std::vector<double>& actual, const std::vector<double>& expected)
{
    QMOMKIT_VERIFY_CHECK(actual.size() == expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        QMOMKIT_VERIFY_CHECK(actual[i] == expected[i]);
    }
}

void assert_near(const double actual, const double expected, const double tolerance)
{
    QMOMKIT_VERIFY_NEAR(actual, expected, tolerance);
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

void assert_raw_moments_reconstructed(const qmomkit::internal::RecursiveQuadratureResult& result, const std::vector<double>& expected_moments,
                                      const std::vector<std::size_t>& extents, const double tolerance)
{
    QMOMKIT_VERIFY_CHECK(result.ok());
    const std::vector<double> reconstructed = raw_moments_from_quadrature(result.quadrature, extents);
    QMOMKIT_VERIFY_CHECK(reconstructed.size() == expected_moments.size());
    for (std::size_t i = 0; i < expected_moments.size(); ++i)
    {
        assert_near(reconstructed[i], expected_moments[i], tolerance);
    }
}

void assert_terminal_moments_reconstructed(const qmomkit::internal::RecursiveQuadratureResult& result, const std::vector<double>& expected_moments,
                                           const std::vector<std::size_t>& extents, const double tolerance)
{
    QMOMKIT_VERIFY_CHECK(result.ok());
    const std::vector<double> reconstructed =
        raw_moments_from_terminal_quadrature(result.quadrature.abscissas(), result.quadrature.weights(), extents);
    QMOMKIT_VERIFY_CHECK(reconstructed.size() == expected_moments.size());
    for (std::size_t i = 0; i < expected_moments.size(); ++i)
    {
        assert_near(reconstructed[i], expected_moments[i], tolerance);
    }
}

void assert_quadrature_matches(const qmomkit::internal::RecursiveQuadrature& actual, const qmomkit::internal::RecursiveQuadrature& expected,
                               const double tolerance)
{
    QMOMKIT_VERIFY_CHECK(actual.hierarchy().node_counts() == expected.hierarchy().node_counts());
    QMOMKIT_VERIFY_CHECK(actual.hierarchy().conditioning_order() == expected.hierarchy().conditioning_order());
    for (const qmomkit::MultiIndex& node_index : qmomkit::enumerate_quadrature_node_multi_indices(expected.hierarchy()))
    {
        assert_near(actual.weight(node_index), expected.weight(node_index), tolerance);
        for (std::size_t dimension = 0; dimension < expected.number_of_dimensions(); ++dimension)
        {
            assert_near(actual.abscissa(dimension, node_index), expected.abscissa(dimension, node_index), tolerance);
        }
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
    QMOMKIT_VERIFY_CHECK(threw);
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
    QMOMKIT_VERIFY_CHECK(threw);
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
    QMOMKIT_VERIFY_CHECK(threw);
}

void assert_ragged_tree_consistent(const qmomkit::internal::RecursiveQuadrature& quadrature)
{
    QMOMKIT_VERIFY_CHECK(quadrature.adaptive_branch_counts_enabled());
    QMOMKIT_VERIFY_CHECK(!quadrature.has_rectangular_storage());
    QMOMKIT_VERIFY_CHECK(!quadrature.tree_nodes().empty());
    QMOMKIT_VERIFY_CHECK(quadrature.terminal_node_ids().size() == quadrature.total_nodes());
    QMOMKIT_VERIFY_CHECK(quadrature.terminal_paths().size() == quadrature.total_nodes());

    const qmomkit::internal::RecursiveQuadrature::TreeNode& root = quadrature.tree_nodes().front();
    QMOMKIT_VERIFY_CHECK(root.parent_id == std::numeric_limits<std::size_t>::max());
    QMOMKIT_VERIFY_CHECK(root.path.empty());
    QMOMKIT_VERIFY_CHECK(!root.terminal);

    for (std::size_t node_id = 0; node_id < quadrature.tree_nodes().size(); ++node_id)
    {
        const qmomkit::internal::RecursiveQuadrature::TreeNode& node = quadrature.tree_nodes()[node_id];
        QMOMKIT_VERIFY_CHECK(node.id == node_id);
        QMOMKIT_VERIFY_CHECK(node.abscissas_by_dimension.size() == quadrature.number_of_dimensions());
        if (node_id == 0)
        {
            continue;
        }

        QMOMKIT_VERIFY_CHECK(node.parent_id < quadrature.tree_nodes().size());
        const qmomkit::internal::RecursiveQuadrature::TreeNode& parent = quadrature.tree_nodes()[node.parent_id];
        QMOMKIT_VERIFY_CHECK(node.path.size() == parent.path.size() + 1);
        for (std::size_t level = 0; level < parent.path.size(); ++level)
        {
            QMOMKIT_VERIFY_CHECK(node.path.node_index_at_level(level) == parent.path.node_index_at_level(level));
        }
        QMOMKIT_VERIFY_CHECK(node.path.node_index_at_level(node.path.size() - 1) == node.local_node_index);
        assert_near(node.accumulated_weight, parent.accumulated_weight * node.local_weight, 1.0e-12);

        bool listed_by_parent = false;
        for (const std::size_t child_id : parent.children)
        {
            listed_by_parent = listed_by_parent || child_id == node_id;
        }
        QMOMKIT_VERIFY_CHECK(listed_by_parent);
    }

    for (std::size_t terminal = 0; terminal < quadrature.total_nodes(); ++terminal)
    {
        const std::size_t node_id = quadrature.terminal_node_ids()[terminal];
        QMOMKIT_VERIFY_CHECK(node_id < quadrature.tree_nodes().size());
        const qmomkit::internal::RecursiveQuadrature::TreeNode& node = quadrature.tree_nodes()[node_id];
        QMOMKIT_VERIFY_CHECK(node.terminal);
        QMOMKIT_VERIFY_CHECK(node.path == quadrature.terminal_paths()[terminal]);
        assert_near(node.accumulated_weight, quadrature.weights()[terminal], 1.0e-12);
        for (std::size_t dimension = 0; dimension < quadrature.number_of_dimensions(); ++dimension)
        {
            assert_near(node.abscissas_by_dimension[dimension], quadrature.abscissas()[dimension][terminal], 1.0e-12);
        }
    }
}

void test_recursive_1d_reduces_to_wheeler()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2});
    const std::vector<double> moments{1.0, 0.7, 1.9, 3.7};

    const qmomkit::Quadrature1D wheeler = qmomkit::compute_quadrature(qmomkit::MomentSet(moments), 2);
    const qmomkit::internal::RecursiveQuadratureResult result = qmomkit::internal::compute_recursive_quadrature({hierarchy, {{4}, moments}, {}, {}});

    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.diagnostics.size() == 1);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.front().conditioning_level == 0);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.front().path.empty());
    for (std::size_t node = 0; node < 2; ++node)
    {
        assert_near(result.quadrature.weight({node}), wheeler.weights[node], 1.0e-12);
        assert_near(result.quadrature.abscissa(0, {node}), wheeler.abscissas[node], 1.0e-12);
    }
}

void test_recursive_2d_identity_conditioning_order()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    qmomkit::internal::RecursiveQuadrature expected(hierarchy);

    const std::vector<double> x{-1.0, 2.0};
    const std::vector<double> wx{0.25, 0.75};
    const std::vector<std::vector<double>> y{{-3.0, 1.0}, {-2.0, 4.0}};
    const std::vector<std::vector<double>> wy{{0.2, 0.8}, {0.6, 0.4}};
    for (std::size_t i = 0; i < 2; ++i)
    {
        for (std::size_t j = 0; j < 2; ++j)
        {
            expected.set_weight({i, j}, wx[i] * wy[i][j]);
            expected.set_abscissa(0, {i, j}, x[i]);
            expected.set_abscissa(1, {i, j}, y[i][j]);
        }
    }

    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> moments = raw_moments_from_quadrature(expected, extents);
    const qmomkit::internal::RecursiveQuadratureResult result =
        qmomkit::internal::compute_recursive_quadrature({hierarchy, {extents, moments}, {}, {}});

    assert_raw_moments_reconstructed(result, moments, extents, 1.0e-10);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.size() == 3);
    QMOMKIT_VERIFY_CHECK(result.quadrature.number_of_dimensions() == 2);
}

void test_recursive_2d_reversed_conditioning_order()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {1, 0});
    qmomkit::internal::RecursiveQuadrature expected(hierarchy);

    const std::vector<double> y{-0.5, 1.25};
    const std::vector<double> wy{0.2, 0.8};
    const std::vector<std::vector<double>> x{{-4.0, -1.0}, {0.25, 2.0}};
    const std::vector<std::vector<double>> wx{{0.35, 0.65}, {0.7, 0.3}};
    for (std::size_t y_node = 0; y_node < 2; ++y_node)
    {
        for (std::size_t x_node = 0; x_node < 2; ++x_node)
        {
            const qmomkit::MultiIndex node_index{x_node, y_node};
            expected.set_weight(node_index, wy[y_node] * wx[y_node][x_node]);
            expected.set_abscissa(0, node_index, x[y_node][x_node]);
            expected.set_abscissa(1, node_index, y[y_node]);
        }
    }

    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> moments = raw_moments_from_quadrature(expected, extents);
    const qmomkit::internal::RecursiveQuadratureResult result =
        qmomkit::internal::compute_recursive_quadrature({hierarchy, {extents, moments}, {}, {}});

    assert_raw_moments_reconstructed(result, moments, extents, 1.0e-10);
    QMOMKIT_VERIFY_CHECK(result.diagnostics[0].phase_space_dimension == 1);
    QMOMKIT_VERIFY_CHECK(result.diagnostics[1].phase_space_dimension == 0);
}

void test_recursive_3d_non_identity_conditioning_order()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 3, 2}, {2, 0, 1});
    qmomkit::internal::RecursiveQuadrature expected(hierarchy);

    const std::vector<double> z{-1.5, 0.75};
    const std::vector<double> wz{0.45, 0.55};
    const std::vector<std::vector<double>> x{{-3.0, 0.5}, {-2.0, 4.0}};
    const std::vector<std::vector<double>> wx{{0.25, 0.75}, {0.4, 0.6}};
    for (std::size_t z_node = 0; z_node < 2; ++z_node)
    {
        for (std::size_t x_node = 0; x_node < 2; ++x_node)
        {
            const std::vector<double> y{-2.0 + static_cast<double>(z_node), 0.25 + 0.5 * static_cast<double>(x_node),
                                        2.5 + static_cast<double>(z_node) + static_cast<double>(x_node)};
            const std::vector<double> wy{0.2 + 0.05 * static_cast<double>(x_node), 0.3, 0.5 - 0.05 * static_cast<double>(x_node)};
            for (std::size_t y_node = 0; y_node < 3; ++y_node)
            {
                const qmomkit::MultiIndex node_index{x_node, y_node, z_node};
                expected.set_weight(node_index, wz[z_node] * wx[z_node][x_node] * wy[y_node]);
                expected.set_abscissa(0, node_index, x[z_node][x_node]);
                expected.set_abscissa(1, node_index, y[y_node]);
                expected.set_abscissa(2, node_index, z[z_node]);
            }
        }
    }

    const std::vector<std::size_t> extents{4, 6, 4};
    const std::vector<double> moments = raw_moments_from_quadrature(expected, extents);
    const qmomkit::internal::RecursiveQuadratureResult result =
        qmomkit::internal::compute_recursive_quadrature({hierarchy, {extents, moments}, {}, {}});

    assert_quadrature_matches(result.quadrature, expected, 1.0e-8);
    assert_raw_moments_reconstructed(result, moments, extents, 1.0e-7);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.front().phase_space_dimension == 2);
}

void test_recursive_non_uniform_fixed_node_counts()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 3});
    qmomkit::internal::RecursiveQuadrature expected(hierarchy);

    const std::vector<double> x{-2.0, 0.5};
    const std::vector<double> wx{0.3, 0.7};
    const std::vector<std::vector<double>> y{{-2.0, 0.0, 3.0}, {-1.0, 1.0, 4.0}};
    const std::vector<std::vector<double>> wy{{0.2, 0.5, 0.3}, {0.25, 0.25, 0.5}};
    for (std::size_t i = 0; i < 2; ++i)
    {
        for (std::size_t j = 0; j < 3; ++j)
        {
            expected.set_weight({i, j}, wx[i] * wy[i][j]);
            expected.set_abscissa(0, {i, j}, x[i]);
            expected.set_abscissa(1, {i, j}, y[i][j]);
        }
    }

    const std::vector<std::size_t> extents{4, 6};
    const std::vector<double> moments = raw_moments_from_quadrature(expected, extents);
    const qmomkit::internal::RecursiveQuadratureResult result =
        qmomkit::internal::compute_recursive_quadrature({hierarchy, {extents, moments}, {}, {}});

    assert_raw_moments_reconstructed(result, moments, extents, 1.0e-8);
    QMOMKIT_VERIFY_CHECK(result.quadrature.total_nodes() == 6);
}

void test_recursive_adaptive_2d_branch_specific_reduction()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.1, 0.3, 0.6};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 1.0}, {-2.0, 2.0, 3.0}};
    const std::vector<double> moments = raw_moments_from_terminal_quadrature(abscissas, weights, extents);

    const qmomkit::internal::RecursiveQuadratureResult result =
        qmomkit::internal::compute_recursive_quadrature({hierarchy, {extents, moments}, {}, {true, 1.0e-10, 1.0e-10, 1.0e-12}, true});

    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.quadrature.adaptive_branch_counts_enabled());
    QMOMKIT_VERIFY_CHECK(!result.quadrature.has_rectangular_storage());
    QMOMKIT_VERIFY_CHECK(result.quadrature.total_nodes() == 3);
    QMOMKIT_VERIFY_CHECK(result.quadrature.terminal_node_ids().size() == 3);
    QMOMKIT_VERIFY_CHECK(result.quadrature.terminal_paths().size() == 3);
    QMOMKIT_VERIFY_CHECK(result.quadrature.terminal_paths()[0] == qmomkit::ConditioningPath({0, 0}));
    QMOMKIT_VERIFY_CHECK(result.quadrature.terminal_paths()[1] == qmomkit::ConditioningPath({0, 1}));
    QMOMKIT_VERIFY_CHECK(result.quadrature.terminal_paths()[2] == qmomkit::ConditioningPath({1, 0}));
    assert_ragged_tree_consistent(result.quadrature);
    assert_throws_logic_error([&result] { (void)result.quadrature.weight({0, 0}); });
    assert_throws_logic_error([&result] { (void)result.quadrature.abscissa(0, {0, 0}); });
    assert_terminal_moments_reconstructed(result, moments, extents, 1.0e-9);

    bool saw_reduced_branch = false;
    for (const qmomkit::internal::RecursiveQuadratureBranchDiagnostics& diagnostics : result.diagnostics)
    {
        QMOMKIT_VERIFY_CHECK(diagnostics.requested_node_count == 2);
        if (diagnostics.path == qmomkit::ConditioningPath({1}) && diagnostics.conditioning_level == 1)
        {
            saw_reduced_branch = true;
            QMOMKIT_VERIFY_CHECK(diagnostics.actual_node_count == 1);
            QMOMKIT_VERIFY_CHECK(diagnostics.adaptive_node_count_reduced);
        }
    }
    QMOMKIT_VERIFY_CHECK(saw_reduced_branch);
}

void test_recursive_adaptive_3d_branch_specific_reduction()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2});
    const std::vector<std::size_t> extents{4, 4, 4};
    const std::vector<double> weights{0.05, 0.05, 0.3, 0.6};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, -1.0, 1.0}, {-2.0, -2.0, 2.0, 3.0}, {-1.0, 1.0, 0.0, 2.0}};
    const std::vector<double> moments = raw_moments_from_terminal_quadrature(abscissas, weights, extents);

    const qmomkit::internal::RecursiveQuadratureResult result =
        qmomkit::internal::compute_recursive_quadrature({hierarchy, {extents, moments}, {}, {true, 1.0e-10, 1.0e-10, 1.0e-12}, true});

    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.quadrature.total_nodes() == 4);
    QMOMKIT_VERIFY_CHECK(result.quadrature.tree_nodes().size() > result.quadrature.total_nodes());
    assert_ragged_tree_consistent(result.quadrature);
    assert_terminal_moments_reconstructed(result, moments, extents, 1.0e-8);

    bool saw_y_reduction = false;
    bool saw_z_reduction = false;
    for (const qmomkit::internal::RecursiveQuadratureBranchDiagnostics& diagnostics : result.diagnostics)
    {
        if (diagnostics.conditioning_level == 1 && diagnostics.path == qmomkit::ConditioningPath({1}))
        {
            saw_y_reduction = true;
            QMOMKIT_VERIFY_CHECK(diagnostics.actual_node_count == 1);
        }
        if (diagnostics.conditioning_level == 2 && diagnostics.actual_node_count == 1)
        {
            saw_z_reduction = true;
            QMOMKIT_VERIFY_CHECK(diagnostics.adaptive_node_count_reduced);
        }
    }
    QMOMKIT_VERIFY_CHECK(saw_y_reduction);
    QMOMKIT_VERIFY_CHECK(saw_z_reduction);
}

void test_recursive_adaptive_non_identity_conditioning_order()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {1, 0});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.25, 0.15, 0.6};
    const std::vector<std::vector<double>> abscissas{{-2.0, 2.0, 3.0}, {-1.0, -1.0, 1.0}};
    const std::vector<double> moments = raw_moments_from_terminal_quadrature(abscissas, weights, extents);

    const qmomkit::internal::RecursiveQuadratureResult result =
        qmomkit::internal::compute_recursive_quadrature({hierarchy, {extents, moments}, {}, {true, 1.0e-10, 1.0e-10, 1.0e-12}, true});

    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.quadrature.total_nodes() == 3);
    assert_ragged_tree_consistent(result.quadrature);
    assert_terminal_moments_reconstructed(result, moments, extents, 1.0e-9);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.front().phase_space_dimension == 1);

    bool saw_reduced_branch = false;
    for (const qmomkit::internal::RecursiveQuadratureBranchDiagnostics& diagnostics : result.diagnostics)
    {
        if (diagnostics.conditioning_level == 1 && diagnostics.path == qmomkit::ConditioningPath({1}))
        {
            saw_reduced_branch = true;
            QMOMKIT_VERIFY_CHECK(diagnostics.phase_space_dimension == 0);
            QMOMKIT_VERIFY_CHECK(diagnostics.requested_node_count == 2);
            QMOMKIT_VERIFY_CHECK(diagnostics.actual_node_count == 1);
        }
    }
    QMOMKIT_VERIFY_CHECK(saw_reduced_branch);
}

} // namespace

int main()
{
    test_recursive_1d_reduces_to_wheeler();
    test_recursive_2d_identity_conditioning_order();
    test_recursive_2d_reversed_conditioning_order();
    test_recursive_3d_non_identity_conditioning_order();
    test_recursive_non_uniform_fixed_node_counts();
    test_recursive_adaptive_2d_branch_specific_reduction();
    test_recursive_adaptive_3d_branch_specific_reduction();
    test_recursive_adaptive_non_identity_conditioning_order();

    return 0;
}
