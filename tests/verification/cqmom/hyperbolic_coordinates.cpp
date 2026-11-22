// CQMOM with hyperbolic coordinates, against references computed from the
// mathematics in exact rational arithmetic and with mpmath at 60 digits.
//
// What is verified
//   That a hyperbolic coordinate of a CQMOM hierarchy is reconstructed by
//   HyQMOM on every branch and that its children are conditioned on its n + 1
//   nodes with the parent orders 0, ..., n, in two dimensions:
//   - hyperbolic_parent_with_a_polynomial_conditional_mean: a hyperbolic
//     parent with n = 2 and a classical child with one node, where the child
//     coordinate is the quadratic f(x) = 1/2 + x/4 - x^2/2 of the parent at
//     each of five weighted points. The child moments M_{p,1}, p = 0, 1, 2,
//     involve parent moments up to m_4, which the HyQMOM rule reproduces, so
//     the terminal nodes are (x_a, f(x_a)) with the HyQMOM nodes x_a and
//     weights w_a of the parent marginal;
//   - hyperbolic_parent_with_polynomial_conditional_moments: a hyperbolic
//     parent with n = 3 and a classical child with two nodes, where each of
//     six parent points carries the child values f(x) - 1/2 and f(x) + 1/2
//     with equal weights, f(x) = 1/4 + x/2. The conditional moments up to
//     order 3 are polynomials of degree at most 3 in x, so the terminal nodes
//     are (x_a, f(x_a) - 1/2) and (x_a, f(x_a) + 1/2) with weights w_a / 2;
//   - hyperbolic_last_coordinate: a classical parent with exactly two points,
//     so that its nodes are those points, and a hyperbolic child with n = 2
//     whose conditional distributions at the two points are four and five
//     weighted points. Each branch returns the HyQMOM rule of its conditional
//     moments;
//   - joint_moments_are_reproduced: two hyperbolic coordinates with n = 2 and
//     conditioning order {1, 0}, for twelve correlated weighted points. The
//     terminal quadrature reproduces every moment of the layout, including
//     those with parent order 2 of the hyperbolic parent;
//   - boundary_branch_adaptive_returns_the_classical_rule and
//     boundary_branch_fixed_fails: as hyperbolic_last_coordinate, but with
//     two points in the second branch, whose conditional moments support
//     only n = 2 nodes. Adaptive mode returns those two points there, and
//     fixed mode fails as qmom_failed in the child coordinate;
//   - hyperbolic_parent_on_the_boundary_conditions_on_its_classical_nodes: a
//     hyperbolic parent with n = 2 at exactly two points, whose moments lie on
//     the boundary, and a classical child with two nodes at two points for each
//     parent point. Adaptive mode returns the two parent points as the
//     classical rule of two nodes, and the children, conditioned on those
//     nodes with the parent orders 0 and 1 of a layout that also holds the
//     parent order 2, return their points.
//   Every case with nodes runs in fixed and in adaptive mode, except the cases
//   on the boundary, which run in adaptive mode or in the mode they name.
//
// How the reference is constructed
//   The script hyperbolic_coordinates.py, beside this file, generated the
//   moments and references stored here. Its points have short binary
//   coordinates and weights, and it asserts that every raw moment is exactly
//   a double, so the moments stored here are the exact moments of the points.
//   A HyQMOM rule is computed from the recurrence of its moments and the
//   closure a_n = (a_0 + ... + a_{n-1}) / n in exact rational arithmetic,
//   then from the eigenvalues and first eigenvector components of J_{n+1}
//   with mpmath at 60 digits, rounded to double. For the two cases whose
//   conditional moments are polynomials of the hyperbolic parent coordinate,
//   the script also solves the Vandermonde systems of the child moments at the
//   HyQMOM nodes at 60 digits and asserts that the conditional moments equal
//   the polynomials of the parent that the points were built from, to 1e-50. It asserts that every sigma(k, k) of every
//   HyQMOM rule, the joint-moment branches included, exceeds the rounding
//   bound QMOMKit applies to it, with the coordinate scale for conditional
//   moments, a thousandfold. The smallest margin is about 1e12.
//
// Error metric and pass criteria
//   Positions relative to the largest |position| of their coordinate in the
//   reference, weights relative to the total weight, and moments relative to
//   max(|moment|, 1), with the tolerance 1e-12 of the two-dimensional cases
//   of README.md. The conditional moments pass through Vandermonde solves of
//   at most four nodes, whose rounding stays far below that tolerance.
//
// Reproduce
//   ctest -R "^verification_cqmom_hyperbolic_coordinates_"
//   The references: python3 hyperbolic_coordinates.py, with mpmath, prints
//   the arrays stored below.
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

constexpr double tolerance = 1.0e-12;

// A terminal node: its path in conditioning order, its weight and its
// abscissas in phase-space order.
struct ReferenceNode
{
    std::vector<std::size_t> path;
    double weight;
    std::vector<double> abscissas;
};

const std::vector<double> polynomial_mean_moments{1.0,        0.1015625,        0.3125,       -0.154296875,      0.953125, -0.44482421875,
                                                  1.09765625, -0.9119873046875, 2.3916015625, -1.964324951171875};
const std::vector<ReferenceNode> polynomial_mean_nodes{
    {{0, 0}, 0.3074990456417663, {-0.8111801444298873, -0.031801649466118194}},
    {{1, 0}, 0.48903186964566037, {0.39060605115736324, 0.5213649691889665}},
    {{2, 0}, 0.20346908471257336, {1.8229713535464966, -0.7058694395389489}},
};
const std::vector<double> polynomial_moments_moments{1.0,
                                                     0.328125,
                                                     0.609375,
                                                     0.56396484375,
                                                     0.15625,
                                                     0.5546875,
                                                     0.4951171875,
                                                     0.9871826171875,
                                                     1.03125,
                                                     0.634765625,
                                                     1.172119140625,
                                                     1.492706298828125,
                                                     0.75390625,
                                                     1.51123046875,
                                                     1.76458740234375,
                                                     3.1724777221679688,
                                                     2.6455078125,
                                                     2.3966064453125,
                                                     3.9514312744140625,
                                                     5.715913772583008,
                                                     3.470458984375,
                                                     5.381805419921875,
                                                     7.059505462646484,
                                                     11.964820384979248,
                                                     9.02838134765625,
                                                     9.692878723144531,
                                                     15.018082618713379,
                                                     22.855962872505188};
const std::vector<ReferenceNode> polynomial_moments_nodes{
    {{0, 0}, 0.09389008092361735, {-1.2850943811923872, -0.8925471905961936}},
    {{0, 1}, 0.09389008092361735, {-1.2850943811923872, 0.10745280940380643}},
    {{1, 0}, 0.19336382241992936, {-0.24991091315038974, -0.37495545657519486}},
    {{1, 1}, 0.19336382241992936, {-0.24991091315038974, 0.6250445434248051}},
    {{2, 0}, 0.15231292853325068, {0.822905151916855, 0.16145257595842746}},
    {{2, 1}, 0.15231292853325068, {0.822905151916855, 1.1614525759584275}},
    {{3, 0}, 0.060433168123202624, {2.014907625462493, 0.7574538127312465}},
    {{3, 1}, 0.060433168123202624, {2.014907625462493, 1.7574538127312465}},
};
const std::vector<double> last_coordinate_moments{
    1.0, 0.53125,   2.0546875,   2.81640625,   10.134765625,   0.75,   0.515625,   2.09765625,   2.888671875,   12.0029296875,
    1.5, 0.9140625, 3.638671875, 5.0009765625, 19.60400390625, 2.0625, 1.30078125, 5.2119140625, 7.16748046875, 28.606201171875};
const std::vector<ReferenceNode> last_coordinate_nodes{
    {{0, 0}, 0.126046937593809, {-0.5, -0.88004229201117}},   {{0, 1}, 0.1644418743386979, {-0.5, 0.4720934807708678}},
    {{0, 2}, 0.08451118806749311, {-0.5, 2.057948811240302}}, {{1, 0}, 0.11759280027111996, {1.5, -1.80964466581535}},
    {{1, 1}, 0.3482433816092918, {1.5, 0.5456845955970924}},  {{1, 2}, 0.15916381811958827, {1.5, 2.597293403551591}},
};
const std::vector<double> boundary_branch_moments{
    1.0, 0.921875,   2.5234375,   4.26171875,   9.255859375,    0.75,   1.1015625,   2.80078125,   5.056640625,    10.6845703125,
    1.5, 1.79296875, 4.693359375, 8.2529296875, 17.62646484375, 2.0625, 2.619140625, 6.7939453125, 12.04541015625, 25.639892578125};
const std::vector<ReferenceNode> boundary_branch_nodes{
    {{0, 0}, 0.126046937593809, {-0.5, -0.88004229201117}},
    {{0, 1}, 0.1644418743386979, {-0.5, 0.4720934807708678}},
    {{0, 2}, 0.08451118806749311, {-0.5, 2.057948811240302}},
    {{1, 0}, 0.15625, {1.5, -1.0}},
    {{1, 1}, 0.46875, {1.5, 2.0}},
};
const std::vector<double> parent_on_the_boundary_moments{1.0,     0.9375,   2.25,     3.75,       0.75,     1.40625,  2.625,
                                                         5.625,   1.5,      2.109375, 4.3125,     8.4375,   2.0625,   3.1640625,
                                                         6.28125, 12.65625, 3.1875,   4.74609375, 9.515625, 18.984375};
const std::vector<ReferenceNode> parent_on_the_boundary_nodes{
    {{0, 0}, 0.1875, {-0.5, -1.0}},
    {{0, 1}, 0.1875, {-0.5, 1.0}},
    {{1, 0}, 0.15625, {1.5, 0.0}},
    {{1, 1}, 0.46875, {1.5, 2.0}},
};
const std::vector<double> joint_moments_moments{1.0,  0.5,     1.5,       2.0,      4.5,     0.5,    0.5,      0.75,    1.25,
                                                2.25, 1.03125, 0.9375,    1.875,    3.0,     6.0,    1.296875, 1.53125, 2.4375,
                                                4.25, 7.875,   2.6484375, 3.140625, 5.53125, 9.5625, 18.375};

qmomkit::PhaseSpaceHierarchy hierarchy_of(std::vector<std::size_t> node_counts, std::vector<std::size_t> conditioning_order,
                                          const std::vector<qmomkit::QmomFamily>& families)
{
    qmomkit::HierarchyFamilies hierarchy_families;
    for (const qmomkit::QmomFamily family : families)
    {
        hierarchy_families.coordinates.push_back({family});
    }
    return {std::move(node_counts), std::move(conditioning_order), hierarchy_families};
}

// The case of the required moments of the hierarchy, selected from the full
// tensor over its moment extents.
qmomkit::CqmomPreparedCase case_of(const qmomkit::PhaseSpaceHierarchy& hierarchy, const std::vector<double>& moments, const bool adaptive)
{
    qmomkit::CqmomCase cqmom_case{
        qmomkit::extract_cqmom_moments({hierarchy.moment_extents(), moments}, qmomkit::create_cqmom_moment_layout(hierarchy))};
    cqmom_case.qmom_options.adaptive = adaptive;
    return qmomkit::CqmomPreparedCase{cqmom_case};
}

void verify_terminal_nodes(const qmomkit::CqmomResult& result, const std::vector<ReferenceNode>& reference)
{
    QMOMKIT_VERIFY_CHECK(result.ok());
    const qmomkit::CqmomQuadrature& quadrature = result.quadrature;
    const qmomkit::PhaseSpaceHierarchy& hierarchy = quadrature.hierarchy();
    QMOMKIT_VERIFY_CHECK(quadrature.total_nodes() == reference.size());

    double total_weight = 0.0;
    std::vector<double> position_scales(hierarchy.number_of_dimensions(), 0.0);
    for (const ReferenceNode& node : reference)
    {
        total_weight += node.weight;
        for (std::size_t dimension = 0; dimension < position_scales.size(); ++dimension)
        {
            position_scales[dimension] = std::max(position_scales[dimension], std::abs(node.abscissas[dimension]));
        }
    }

    // Fixed results store the nodes by phase-space node index, adaptive
    // results by terminal path.
    const std::vector<qmomkit::ConditioningPath> paths = quadrature.terminal_paths().to_conditioning_paths();
    for (const ReferenceNode& node : reference)
    {
        const qmomkit::ConditioningPath path(node.path);
        std::size_t position = 0;
        if (quadrature.has_rectangular_storage())
        {
            position = hierarchy.flatten_node_index(path.phase_space_node_index(hierarchy));
        }
        else
        {
            const auto found = std::find(paths.begin(), paths.end(), path);
            QMOMKIT_VERIFY_CHECK(found != paths.end());
            position = static_cast<std::size_t>(found - paths.begin());
        }
        QMOMKIT_VERIFY_NEAR(quadrature.weights()[position], node.weight, tolerance * total_weight);
        for (std::size_t dimension = 0; dimension < position_scales.size(); ++dimension)
        {
            QMOMKIT_VERIFY_NEAR(quadrature.abscissas()[dimension][position], node.abscissas[dimension], tolerance * position_scales[dimension]);
        }
    }
}

void verify_both_modes(const qmomkit::PhaseSpaceHierarchy& hierarchy, const std::vector<double>& moments, const std::vector<ReferenceNode>& reference)
{
    for (const bool adaptive : {false, true})
    {
        verify_terminal_nodes(qmomkit::compute_cqmom(case_of(hierarchy, moments, adaptive)), reference);
    }
}

void verify_hyperbolic_parent_with_a_polynomial_conditional_mean()
{
    verify_both_modes(hierarchy_of({2, 1}, {0, 1}, {qmomkit::QmomFamily::hyperbolic, qmomkit::QmomFamily::classical}), polynomial_mean_moments,
                      polynomial_mean_nodes);
}

void verify_hyperbolic_parent_with_polynomial_conditional_moments()
{
    verify_both_modes(hierarchy_of({3, 2}, {0, 1}, {qmomkit::QmomFamily::hyperbolic, qmomkit::QmomFamily::classical}), polynomial_moments_moments,
                      polynomial_moments_nodes);
}

void verify_hyperbolic_last_coordinate()
{
    verify_both_modes(hierarchy_of({2, 2}, {0, 1}, {qmomkit::QmomFamily::classical, qmomkit::QmomFamily::hyperbolic}), last_coordinate_moments,
                      last_coordinate_nodes);
}

void verify_hyperbolic_parent_on_the_boundary_conditions_on_its_classical_nodes()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy = hierarchy_of({2, 2}, {0, 1}, {qmomkit::QmomFamily::hyperbolic, qmomkit::QmomFamily::classical});
    verify_terminal_nodes(qmomkit::compute_cqmom(case_of(hierarchy, parent_on_the_boundary_moments, true)), parent_on_the_boundary_nodes);
}

void verify_joint_moments_are_reproduced()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy = hierarchy_of({2, 2}, {1, 0}, {qmomkit::QmomFamily::hyperbolic, qmomkit::QmomFamily::hyperbolic});
    for (const bool adaptive : {false, true})
    {
        const qmomkit::CqmomPreparedCase prepared_case = case_of(hierarchy, joint_moments_moments, adaptive);
        const qmomkit::CqmomResult result = qmomkit::compute_cqmom(prepared_case);
        QMOMKIT_VERIFY_CHECK(result.ok());
        QMOMKIT_VERIFY_CHECK(result.quadrature.total_nodes() == 9);
        const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
        const qmomkit::CqmomMomentData expected = qmomkit::extract_cqmom_moments({hierarchy.moment_extents(), joint_moments_moments}, layout);
        const qmomkit::CqmomMomentData reconstructed = qmomkit::reconstruct_cqmom_moments(result.quadrature, layout);
        QMOMKIT_VERIFY_CHECK(reconstructed.values.size() == expected.values.size());
        for (std::size_t position = 0; position < expected.values.size(); ++position)
        {
            QMOMKIT_VERIFY_NEAR(reconstructed.values[position], expected.values[position],
                                tolerance * std::max(std::abs(expected.values[position]), 1.0));
        }
    }
}

void verify_boundary_branch_adaptive_returns_the_classical_rule()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy = hierarchy_of({2, 2}, {0, 1}, {qmomkit::QmomFamily::classical, qmomkit::QmomFamily::hyperbolic});
    verify_terminal_nodes(qmomkit::compute_cqmom(case_of(hierarchy, boundary_branch_moments, true)), boundary_branch_nodes);
}

void verify_boundary_branch_fixed_fails()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy = hierarchy_of({2, 2}, {0, 1}, {qmomkit::QmomFamily::classical, qmomkit::QmomFamily::hyperbolic});
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(case_of(hierarchy, boundary_branch_moments, false));
    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::CqmomStatus::qmom_failed);
    QMOMKIT_VERIFY_CHECK(result.failed_phase_space_dimension == 1);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"hyperbolic_parent_with_a_polynomial_conditional_mean", verify_hyperbolic_parent_with_a_polynomial_conditional_mean},
            {"hyperbolic_parent_with_polynomial_conditional_moments", verify_hyperbolic_parent_with_polynomial_conditional_moments},
            {"hyperbolic_last_coordinate", verify_hyperbolic_last_coordinate},
            {"joint_moments_are_reproduced", verify_joint_moments_are_reproduced},
            {"boundary_branch_adaptive_returns_the_classical_rule", verify_boundary_branch_adaptive_returns_the_classical_rule},
            {"boundary_branch_fixed_fails", verify_boundary_branch_fixed_fails},
            {"hyperbolic_parent_on_the_boundary_conditions_on_its_classical_nodes",
             verify_hyperbolic_parent_on_the_boundary_conditions_on_its_classical_nodes},
        });
}
