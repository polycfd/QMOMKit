// The recursive CQMOM inversion, against quadratures written down in advance.
//
// What is verified
//   That the recursion recovers the quadrature the moments came from, for one
//   to three dimensions, for identity and non-identity conditioning orders, for
//   non-uniform node counts per axis, and for adaptive cases where individual
//   branches must return fewer nodes than requested. Also that the one
//   dimensional case reduces to plain Wheeler, and that a branch which collapses
//   does so for either conditioning order.
//
// How the reference is constructed
//   Chosen, not computed. Each case builds the answer first - positions and
//   weights placed directly into a RecursiveQuadrature - and then evaluates its
//   raw moments with a local helper that sums weight times position to the
//   required power. That helper is written in this file and does not call the
//   library's reconstruction, so the reference does not depend on the code path
//   under test. Feeding those moments back in must return the quadrature the
//   case started from.
//
// Error metric and pass criteria
//   Maximum absolute deviation of the rebuilt moments from the prescribed ones,
//   with the branch node counts and diagnostic record counts required to match
//   exactly. Tolerances are per case: 1.0e-10 for the well-separated fixed
//   cases, loosening to 1.0e-8 and 1.0e-7 where the case deliberately places
//   positions close together or drives a branch to collapse, since the solve
//   behind those is correspondingly worse conditioned.
//
// Reproduce
//   ctest -R verification_recursive_quadrature
//
// Build modes
//   Both. Failures are reported through an exit code, not assert().

#include "internal/recursive_quadrature.hpp"

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

void assert_terminal_paths_equal(const qmomkit::internal::RecursiveQuadrature& quadrature, const std::vector<qmomkit::ConditioningPath>& expected)
{
    QMOMKIT_VERIFY_CHECK(quadrature.terminal_paths().size() == expected.size());
    for (std::size_t path = 0; path < expected.size(); ++path)
    {
        QMOMKIT_VERIFY_CHECK(quadrature.terminal_paths()[path] == expected[path]);
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

std::vector<double> product_moments_from_1d_sequences(const std::vector<double>& x_moments, const std::vector<double>& y_moments,
                                                      const std::vector<std::size_t>& extents)
{
    std::vector<double> moments(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat_moment = 0; flat_moment < moments.size(); ++flat_moment)
    {
        const qmomkit::MultiIndex moment_index = qmomkit::unflatten_multi_index(flat_moment, extents);
        moments[flat_moment] = x_moments.at(moment_index[0]) * y_moments.at(moment_index[1]);
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
    // Structural dimensional-reduction consistency: 1D recursive CQMOM must match Wheeler; independent Wheeler correctness is covered by separate
    // analytical tests.
    const qmomkit::PhaseSpaceHierarchy hierarchy({2});
    const std::vector<double> moments{1.0, 0.7, 1.9, 3.7};

    const qmomkit::Quadrature1D wheeler = qmomkit::compute_quadrature(qmomkit::MomentSet(moments), 2);
    const qmomkit::internal::RecursiveQuadratureResult result = qmomkit::internal::compute_recursive_quadrature_from_moment_data(
        moment_data_request(hierarchy, {{4}, moments}, {}, qmomkit::WheelerOptions{false}, qmomkit::CqmomDiagnosticLevel::basic));

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
    const qmomkit::internal::RecursiveQuadratureResult result = qmomkit::internal::compute_recursive_quadrature_from_moment_data(
        moment_data_request(hierarchy, {extents, moments}, {}, qmomkit::WheelerOptions{false}, qmomkit::CqmomDiagnosticLevel::basic));

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
    const qmomkit::internal::RecursiveQuadratureResult result = qmomkit::internal::compute_recursive_quadrature_from_moment_data(
        moment_data_request(hierarchy, {extents, moments}, {}, qmomkit::WheelerOptions{false}, qmomkit::CqmomDiagnosticLevel::basic));

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
    const qmomkit::internal::RecursiveQuadratureResult result = qmomkit::internal::compute_recursive_quadrature_from_moment_data(
        moment_data_request(hierarchy, {extents, moments}, {}, qmomkit::WheelerOptions{false}, qmomkit::CqmomDiagnosticLevel::basic));

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
    const qmomkit::internal::RecursiveQuadratureResult result = qmomkit::internal::compute_recursive_quadrature_from_moment_data(
        moment_data_request(hierarchy, {extents, moments}, {}, qmomkit::WheelerOptions{false}, qmomkit::CqmomDiagnosticLevel::basic));

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

    const qmomkit::internal::RecursiveQuadratureResult result = qmomkit::internal::compute_recursive_quadrature_from_moment_data(
        moment_data_request(hierarchy, {extents, moments}, {}, {true, 1.0e-10, 1.0e-10, 1.0e-12}, qmomkit::CqmomDiagnosticLevel::basic));

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

    const qmomkit::internal::RecursiveQuadratureResult result = qmomkit::internal::compute_recursive_quadrature_from_moment_data(
        moment_data_request(hierarchy, {extents, moments}, {}, {true, 1.0e-10, 1.0e-10, 1.0e-12}, qmomkit::CqmomDiagnosticLevel::basic));

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

    const qmomkit::internal::RecursiveQuadratureResult result = qmomkit::internal::compute_recursive_quadrature_from_moment_data(
        moment_data_request(hierarchy, {extents, moments}, {}, {true, 1.0e-10, 1.0e-10, 1.0e-12}, qmomkit::CqmomDiagnosticLevel::basic));

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

void test_recursive_adaptive_velocity_collapse_reduces_for_both_conditioning_orders()
{
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> size_moments{1.0, 0.0, 1.0, 0.0};
    const std::vector<double> collapsing_velocity_moments{1.0, 0.0, -1.0, 0.0};
    const std::vector<double> moments = product_moments_from_1d_sequences(size_moments, collapsing_velocity_moments, extents);
    const qmomkit::WheelerOptions wheeler_options{true, 1.0e-10, 1.0e-10, 1.0e-12};

    // With order [0,1], the same physical velocity dimension is conditional.
    // Each size root keeps one velocity child after the negative-variance retry.
    const qmomkit::internal::RecursiveQuadratureResult conditional =
        qmomkit::internal::compute_recursive_quadrature_from_moment_data(moment_data_request(
            qmomkit::PhaseSpaceHierarchy({2, 2}, {0, 1}), {extents, moments}, {}, wheeler_options, qmomkit::CqmomDiagnosticLevel::basic));
    QMOMKIT_VERIFY_CHECK(conditional.ok());
    QMOMKIT_VERIFY_CHECK(conditional.quadrature.total_nodes() == 2);
    assert_ragged_tree_consistent(conditional.quadrature);
    assert_terminal_paths_equal(conditional.quadrature, {qmomkit::ConditioningPath({0, 0}), qmomkit::ConditioningPath({1, 0})});

    std::size_t conditional_reduced_branches = 0;
    for (const qmomkit::internal::RecursiveQuadratureBranchDiagnostics& diagnostics : conditional.diagnostics)
    {
        if (diagnostics.conditioning_level == 1 && diagnostics.phase_space_dimension == 1)
        {
            QMOMKIT_VERIFY_CHECK(diagnostics.requested_node_count == 2);
            QMOMKIT_VERIFY_CHECK(diagnostics.actual_node_count == 1);
            QMOMKIT_VERIFY_CHECK(diagnostics.adaptive_node_count_reduced);
            ++conditional_reduced_branches;
        }
    }
    QMOMKIT_VERIFY_CHECK(conditional_reduced_branches == 2);

    // With order [1,0], the velocity dimension is marginal. The root velocity
    // quadrature reduces to one node, and that root retains two size children.
    const qmomkit::internal::RecursiveQuadratureResult marginal =
        qmomkit::internal::compute_recursive_quadrature_from_moment_data(moment_data_request(
            qmomkit::PhaseSpaceHierarchy({2, 2}, {1, 0}), {extents, moments}, {}, wheeler_options, qmomkit::CqmomDiagnosticLevel::basic));
    QMOMKIT_VERIFY_CHECK(marginal.ok());
    QMOMKIT_VERIFY_CHECK(marginal.quadrature.total_nodes() == 2);
    assert_ragged_tree_consistent(marginal.quadrature);
    assert_terminal_paths_equal(marginal.quadrature, {qmomkit::ConditioningPath({0, 0}), qmomkit::ConditioningPath({0, 1})});
    QMOMKIT_VERIFY_CHECK(!marginal.diagnostics.empty());
    QMOMKIT_VERIFY_CHECK(marginal.diagnostics.front().conditioning_level == 0);
    QMOMKIT_VERIFY_CHECK(marginal.diagnostics.front().path.empty());
    QMOMKIT_VERIFY_CHECK(marginal.diagnostics.front().phase_space_dimension == 1);
    QMOMKIT_VERIFY_CHECK(marginal.diagnostics.front().requested_node_count == 2);
    QMOMKIT_VERIFY_CHECK(marginal.diagnostics.front().actual_node_count == 1);
    QMOMKIT_VERIFY_CHECK(marginal.diagnostics.front().adaptive_node_count_reduced);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"recursive_1d_reduces_to_wheeler", test_recursive_1d_reduces_to_wheeler},
            {"recursive_2d_identity_conditioning_order", test_recursive_2d_identity_conditioning_order},
            {"recursive_2d_reversed_conditioning_order", test_recursive_2d_reversed_conditioning_order},
            {"recursive_3d_non_identity_conditioning_order", test_recursive_3d_non_identity_conditioning_order},
            {"recursive_non_uniform_fixed_node_counts", test_recursive_non_uniform_fixed_node_counts},
            {"recursive_adaptive_2d_branch_specific_reduction", test_recursive_adaptive_2d_branch_specific_reduction},
            {"recursive_adaptive_3d_branch_specific_reduction", test_recursive_adaptive_3d_branch_specific_reduction},
            {"recursive_adaptive_non_identity_conditioning_order", test_recursive_adaptive_non_identity_conditioning_order},
            {"recursive_adaptive_velocity_collapse_reduces_for_both_conditioning_orders",
             test_recursive_adaptive_velocity_collapse_reduces_for_both_conditioning_orders},
        });
}
