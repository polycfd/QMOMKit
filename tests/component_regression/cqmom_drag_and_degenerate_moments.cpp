/*
 * WHAT IS TESTED
 *   That three-dimensional CQMOM reconstruction handles two moment sets that
 *   are harder to invert than a rule whose nodes are far apart: one whose
 *   velocities have been changed by drag, and one with two velocity nodes
 *   almost on top of each other.
 *
 * HOW IT IS TESTED
 *   Both cases start from the same explicit list of eight quadrature nodes and
 *   compute their raw moments from those nodes, so each moment set is
 *   realizable by construction and exact to round-off.
 *
 *   The first case applies drag to every node for one time step: each node's
 *   velocity moves toward the surrounding fluid velocity by the exponential
 *   factor exp(-rate * time step). Both the drag rate and the fluid velocity
 *   depend on the node's other two coordinates, so every node moves by a
 *   different amount, and the velocities end up closer together than the ones
 *   the case was built from.
 *
 *   The second case instead puts two velocity nodes of one branch 1.5e-6
 *   apart. Moments of two nodes that close are almost the moments of a single
 *   node, so that branch's moment sequence carries very little information
 *   about where the two nodes actually are.
 *
 * PASS CRITERIA
 *   Both cases reconstruct successfully with eight terminal nodes and reproduce
 *   every required moment to 1e-12, absolute or relative.
 *
 *   The drag case additionally returns the nodes it was built from, to 1e-10.
 *   The case with two nodes 1.5e-6 apart deliberately does not check that:
 *   their positions are not determined by the moments, so returning the
 *   original pair is not required and any pair reproducing the moments is a
 *   correct answer. Reproducing the moments is the whole of the contract there.
 *
 *   One branch of the drag case is also inverted on its own with
 *   one-dimensional QMOM, which shows that the moment sequence handed to that
 *   branch is usable by itself, so a failure inside the recursive
 *   reconstruction could not be blamed on the branch data.
 *
 * WHEN IT IS TESTED
 *   Debug and Release. The checks are explicit and do not rely on assert().
 */

#include "../verification/verification_helpers.hpp"
#include "qmomkit/qmomkit.hpp"
#include "test_case_dispatch.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace
{

constexpr std::size_t coordinate_count = 3;

struct TerminalNode
{
    double weight = 0.0;
    std::array<double, coordinate_count> coordinate{};
    qmomkit::ConditioningPath path;
};

[[nodiscard]] double integer_power(const double base, const std::size_t exponent)
{
    double value = 1.0;
    for (std::size_t step = 0; step < exponent; ++step)
    {
        value *= base;
    }
    return value;
}

// Eight nodes on a 2x2x2 conditioning tree, with branch-dependent conditional
// positions so the case is not an independent product.
[[nodiscard]] std::vector<TerminalNode> generating_nodes()
{
    const std::array<double, 2> d_nodes{0.82, 1.66};
    const std::array<double, 2> d_weights{0.43, 0.57};
    const std::array<std::array<double, 2>, 2> chi_nodes{{{{0.30, 0.94}}, {{0.48, 1.32}}}};
    const std::array<std::array<double, 2>, 2> chi_weights{{{{0.37, 0.63}}, {{0.59, 0.41}}}};
    const std::array<std::array<std::array<double, 2>, 2>, 2> u_nodes{{{{{{-0.35, 1.15}}, {{0.15, 1.80}}}}, {{{{-0.85, 0.90}}, {{0.42, 2.35}}}}}};
    const std::array<std::array<std::array<double, 2>, 2>, 2> u_weights{{{{{{0.62, 0.38}}, {{0.44, 0.56}}}}, {{{{0.71, 0.29}}, {{0.52, 0.48}}}}}};

    std::vector<TerminalNode> nodes;
    for (std::size_t i = 0; i < 2; ++i)
    {
        for (std::size_t j = 0; j < 2; ++j)
        {
            for (std::size_t k = 0; k < 2; ++k)
            {
                nodes.push_back(TerminalNode{d_weights[i] * chi_weights[i][j] * u_weights[i][j][k],
                                             {d_nodes[i], chi_nodes[i][j], u_nodes[i][j][k]},
                                             qmomkit::ConditioningPath{{i, j, k}}});
            }
        }
    }
    return nodes;
}

// Both the drag rate and the fluid velocity depend on a node's first two
// coordinates, so every node's velocity changes by a different amount and the
// velocities do not simply shift together.
[[nodiscard]] double drag_rate(const double d, const double chi)
{
    return (1.0 + 0.35 * chi + 0.08 * chi * chi) / (d * d);
}

[[nodiscard]] double fluid_velocity(const double d, const double chi)
{
    return 0.35 + 0.25 * chi - 0.10 * d;
}

[[nodiscard]] std::vector<TerminalNode> drag_transformed_nodes()
{
    constexpr double time_step = 0.7;
    std::vector<TerminalNode> nodes = generating_nodes();
    for (TerminalNode& node : nodes)
    {
        const double rate = drag_rate(node.coordinate[0], node.coordinate[1]);
        const double target = fluid_velocity(node.coordinate[0], node.coordinate[1]);
        const double initial_velocity = node.coordinate[2];
        node.coordinate[2] = target + (initial_velocity - target) * std::exp(-rate * time_step);
        // Drag moves a velocity toward the fluid velocity and never past it, so
        // the result lies between where the node started and the fluid
        // velocity. This guards the arithmetic in this file, not QMOMKit.
        QMOMKIT_VERIFY_CHECK(node.coordinate[2] >= std::min(initial_velocity, target));
        QMOMKIT_VERIFY_CHECK(node.coordinate[2] <= std::max(initial_velocity, target));
    }
    return nodes;
}

[[nodiscard]] std::vector<TerminalNode> near_degenerate_nodes()
{
    std::vector<TerminalNode> nodes = generating_nodes();
    for (TerminalNode& node : nodes)
    {
        if (node.path.node_indices_by_level() == qmomkit::MultiIndex{1, 1, 0})
        {
            node.coordinate[2] = 0.7500000;
        }
        if (node.path.node_indices_by_level() == qmomkit::MultiIndex{1, 1, 1})
        {
            node.coordinate[2] = 0.7500015;
        }
    }
    return nodes;
}

[[nodiscard]] qmomkit::RawMomentTensor raw_moments_from_nodes(const std::vector<TerminalNode>& nodes, const std::vector<std::size_t>& extents)
{
    std::vector<double> values(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat = 0; flat < values.size(); ++flat)
    {
        const qmomkit::MultiIndex order = qmomkit::unflatten_multi_index(flat, extents);
        double moment = 0.0;
        for (const TerminalNode& node : nodes)
        {
            double contribution = node.weight;
            for (std::size_t dimension = 0; dimension < order.size(); ++dimension)
            {
                contribution *= integer_power(node.coordinate[dimension], order[dimension]);
            }
            moment += contribution;
        }
        values[flat] = moment;
    }
    return {extents, values};
}

// One options object for both the recursive reconstruction and the direct
// branch inversion below, so the direct check runs under exactly the numerical
// settings the recursive path uses.
[[nodiscard]] qmomkit::QmomOptions reconstruction_options()
{
    qmomkit::QmomOptions options;
    options.adaptive = false;
    options.rmin = 1.0e-12;
    options.eabs = 1.0e-12;
    return options;
}

[[nodiscard]] qmomkit::CqmomMomentData required_moments_of(const std::vector<TerminalNode>& nodes, const qmomkit::PhaseSpaceHierarchy& hierarchy)
{
    return qmomkit::extract_cqmom_moments(raw_moments_from_nodes(nodes, hierarchy.moment_extents()), qmomkit::create_cqmom_moment_layout(hierarchy));
}

[[nodiscard]] qmomkit::CqmomResult reconstruct(const qmomkit::CqmomMomentData& moments)
{
    return qmomkit::compute_cqmom(qmomkit::CqmomCase{moments, {}, reconstruction_options(), false});
}

// Reconstructing the moments from the accepted quadrature must return the
// moments the case was built from, whether or not the node positions
// themselves are determined.
void check_reproduces_required_moments(const qmomkit::CqmomMomentData& expected, const qmomkit::CqmomResult& result)
{
    const qmomkit::CqmomMomentData produced = qmomkit::reconstruct_cqmom_moments(result.quadrature, expected.layout);
    QMOMKIT_VERIFY_CHECK(produced.values.size() == expected.values.size());
    for (std::size_t entry = 0; entry < expected.values.size(); ++entry)
    {
        QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(produced.values[entry], expected.values[entry], 1.0e-12, 1.0e-12);
    }
}

// Each generating node must be matched by exactly one accepted node, so a
// reconstruction that collapsed two nodes onto one cannot pass. Matching each
// accepted node to its nearest unclaimed generating node is unambiguous here
// because the generating nodes are far apart relative to the tolerance; this
// helper is therefore not used for the near-degenerate case.
void check_recovers_generating_nodes(const std::vector<TerminalNode>& nodes, const qmomkit::CqmomResult& result, const double tolerance)
{
    QMOMKIT_VERIFY_CHECK(result.quadrature.weights().size() == nodes.size());
    std::vector<bool> claimed(nodes.size(), false);
    for (std::size_t terminal = 0; terminal < result.quadrature.weights().size(); ++terminal)
    {
        std::size_t closest = nodes.size();
        double closest_distance = std::numeric_limits<double>::infinity();
        for (std::size_t candidate = 0; candidate < nodes.size(); ++candidate)
        {
            if (claimed[candidate])
            {
                continue;
            }
            double distance = 0.0;
            for (std::size_t dimension = 0; dimension < coordinate_count; ++dimension)
            {
                distance = std::max(distance, std::abs(result.quadrature.abscissas()[dimension][terminal] - nodes[candidate].coordinate[dimension]));
            }
            if (distance < closest_distance)
            {
                closest_distance = distance;
                closest = candidate;
            }
        }
        QMOMKIT_VERIFY_CHECK(closest < nodes.size());
        claimed[closest] = true;
        QMOMKIT_VERIFY_NEAR(closest_distance, 0.0, tolerance);
        QMOMKIT_VERIFY_NEAR(result.quadrature.weights()[terminal], nodes[closest].weight, tolerance);
    }
}

void test_drag_transformed_moment_set_reconstructs()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2});
    const std::vector<TerminalNode> nodes = drag_transformed_nodes();
    const qmomkit::CqmomMomentData moments = required_moments_of(nodes, hierarchy);
    const qmomkit::CqmomResult result = reconstruct(moments);

    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.quadrature.total_nodes() == nodes.size());
    check_reproduces_required_moments(moments, result);
    check_recovers_generating_nodes(nodes, result, 1.0e-10);
}

// The conditional moment sequence a branch hands to one-dimensional QMOM must
// be realizable on its own, so a failure inside recursive reconstruction can be
// attributed to the recursion rather than to the branch data.
void test_drag_transformed_branch_sequence_inverts_directly()
{
    const std::vector<TerminalNode> nodes = drag_transformed_nodes();
    const qmomkit::ConditioningPath branch{{0, 1}};

    double branch_mass = 0.0;
    std::vector<TerminalNode> branch_nodes;
    for (const TerminalNode& node : nodes)
    {
        const qmomkit::MultiIndex indices = node.path.node_indices_by_level();
        const qmomkit::MultiIndex branch_indices = branch.node_indices_by_level();
        if (indices.size() >= branch_indices.size() && std::equal(branch_indices.begin(), branch_indices.end(), indices.begin()))
        {
            branch_nodes.push_back(node);
            branch_mass += node.weight;
        }
    }
    QMOMKIT_VERIFY_CHECK(!branch_nodes.empty());
    QMOMKIT_VERIFY_CHECK(branch_mass > 0.0);

    constexpr std::size_t moment_count = 4;
    std::vector<double> conditional_moments(moment_count, 0.0);
    for (const TerminalNode& node : branch_nodes)
    {
        const double conditional_weight = node.weight / branch_mass;
        for (std::size_t order = 0; order < moment_count; ++order)
        {
            conditional_moments[order] += conditional_weight * integer_power(node.coordinate[2], order);
        }
    }

    const qmomkit::QmomOptions qmom_options = reconstruction_options();
    const qmomkit::WheelerOptions wheeler_options{qmom_options.adaptive, qmom_options.rmin, qmom_options.eabs, qmom_options.cutoff};
    const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(qmomkit::MomentSet(conditional_moments), 2, wheeler_options);
    QMOMKIT_VERIFY_CHECK(quadrature.ok());
    QMOMKIT_VERIFY_CHECK(quadrature.weights.size() == 2);
}

void test_near_degenerate_branch_moment_set_reconstructs()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2});
    const std::vector<TerminalNode> nodes = near_degenerate_nodes();
    const qmomkit::CqmomMomentData moments = required_moments_of(nodes, hierarchy);
    const qmomkit::CqmomResult result = reconstruct(moments);

    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.quadrature.total_nodes() == nodes.size());
    check_reproduces_required_moments(moments, result);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"drag_transformed_moment_set_reconstructs", test_drag_transformed_moment_set_reconstructs},
            {"drag_transformed_branch_sequence_inverts_directly", test_drag_transformed_branch_sequence_inverts_directly},
            {"near_degenerate_branch_moment_set_reconstructs", test_near_degenerate_branch_moment_set_reconstructs},
        });
}
