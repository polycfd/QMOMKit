// QMOM inversion when the caller restricts the positions to an interval.
//
// What is verified
//   That a fixed-node inversion is accepted when its positions already lie
//   inside the requested interval; that it is rejected, rather than quietly
//   moved onto the boundary, when a position falls outside; that the
//   adaptive path returns the largest node count whose positions all fit; and
//   that whether a position on or just beyond a bound is accepted does not
//   depend on the unit of the coordinate.
//
// How the reference is constructed
//   Chosen, not computed. Each case starts from positions and weights picked
//   to sit either inside or outside the interval under test, and the moments
//   are built from them. The expected outcome is therefore known before the
//   inversion runs: which node count survives, and where the surviving nodes
//   must be.
//
// Error metric and pass criteria
//   Status must match the expected outcome exactly. Where a quadrature is
//   expected, positions and weights must agree with the chosen values to
//   1.0e-12 absolute. The rejection case additionally checks that the returned
//   position is unchanged, which is what distinguishes rejection from
//   clipping. The unit case repeats the same rules at coordinate scales 1e-6,
//   1 and 1e6 and requires the same decision at every scale: a single position
//   half the acceptance margin beyond a bound is accepted and one twice the
//   margin beyond is rejected, where the margin is 1e-12 times the larger bound
//   magnitude; rules with a position exactly on a bound are accepted.
//
// Reproduce
//   ctest -R verification_qmom_.*bounded
//
// Build modes
//   Both. Checks go through QMOMKIT_VERIFY_*, not assert().

#include "../verification_helpers.hpp"
#include "qmomkit/qmom.hpp"
#include "test_case_dispatch.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace
{

constexpr double tolerance = 1.0e-12;

std::vector<double> raw_moments(const std::vector<double>& weights, const std::vector<double>& nodes)
{
    std::vector<double> moments(2 * nodes.size(), 0.0);
    for (std::size_t order = 0; order < moments.size(); ++order)
    {
        for (std::size_t node = 0; node < nodes.size(); ++node)
        {
            moments[order] += weights[node] * std::pow(nodes[node], static_cast<int>(order));
        }
    }
    return moments;
}

void verify_fixed_bounded_qmom_accepts_nodes_inside_interval()
{
    const std::vector<double> weights{0.4, 0.6};
    const std::vector<double> nodes{2.25, 4.5};
    const qmomkit::MomentSet moments(raw_moments(weights, nodes));
    const qmomkit::AbscissaBounds bounds{2.0, 5.0};

    const qmomkit::QmomResult bounded = qmomkit::compute_qmom(moments, 2, qmomkit::QmomOptions{false}, bounds);
    const qmomkit::QmomResult unbounded = qmomkit::compute_qmom(moments, 2, qmomkit::QmomOptions{false});

    QMOMKIT_VERIFY_CHECK(bounded.ok());
    QMOMKIT_VERIFY_CHECK(unbounded.ok());
    QMOMKIT_VERIFY_CHECK(bounded.abscissas == unbounded.abscissas);
    QMOMKIT_VERIFY_CHECK(bounded.weights == unbounded.weights);
}

void verify_fixed_bounded_qmom_rejects_nodes_outside_interval_without_clipping()
{
    const std::vector<double> weights{0.5, 0.5};
    const std::vector<double> nodes{0.25, 1.25};
    const qmomkit::QmomResult result =
        qmomkit::compute_qmom(qmomkit::MomentSet(raw_moments(weights, nodes)), 2, qmomkit::QmomOptions{false}, qmomkit::AbscissaBounds{0.0, 1.0});

    QMOMKIT_VERIFY_CHECK(!result.ok());
    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::QmomStatus::not_realizable);
    QMOMKIT_VERIFY_CHECK(result.abscissas.size() == nodes.size());
    QMOMKIT_VERIFY_NEAR(result.abscissas[1], nodes[1], tolerance);
    QMOMKIT_VERIFY_CHECK(result.abscissas[1] > 1.0);
}

void verify_adaptive_bounded_qmom_uses_highest_admissible_count()
{
    const std::vector<double> weights{0.5, 0.5};
    const std::vector<double> nodes{0.25, 1.25};
    qmomkit::QmomOptions options;
    options.adaptive = true;

    const qmomkit::QmomResult result =
        qmomkit::compute_qmom(qmomkit::MomentSet(raw_moments(weights, nodes)), 2, options, qmomkit::AbscissaBounds{0.0, 1.0});

    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.weights.size() == 1);
    QMOMKIT_VERIFY_CHECK(result.abscissas.size() == 1);
    QMOMKIT_VERIFY_NEAR(result.weights[0], 1.0, tolerance);
    QMOMKIT_VERIFY_NEAR(result.abscissas[0], 0.75, tolerance);
}

void verify_bounded_qmom_acceptance_does_not_depend_on_the_unit()
{
    // Rules with a position exactly on a bound. Rounding moves that position by
    // an amount proportional to the coordinate, so a margin fixed in absolute
    // terms would reject some of them at a large scale.
    struct Rule
    {
        std::vector<double> weights;
        std::vector<double> nodes;
    };
    const std::vector<Rule> rules{{{0.3, 0.7}, {0.0, 0.6}},
                                  {{0.6, 0.4}, {0.35, 1.0}},
                                  {{0.2, 0.5, 0.3}, {0.0, 0.4, 0.8}},
                                  {{0.25, 0.35, 0.4}, {0.15, 0.55, 1.0}},
                                  {{0.1, 0.3, 0.4, 0.2}, {0.0, 0.3, 0.6, 0.85}},
                                  {{0.2, 0.3, 0.3, 0.2}, {0.1, 0.4, 0.7, 1.0}}};

    for (const double scale : {1.0e-6, 1.0, 1.0e6})
    {
        const qmomkit::AbscissaBounds bounds{0.0, scale};
        const double margin = 1.0e-12 * std::max(std::abs(bounds.lower), std::abs(bounds.upper));
        const auto single_position = [&bounds](const double node)
        { return qmomkit::compute_qmom(qmomkit::MomentSet(raw_moments({0.3}, {node})), 1, qmomkit::QmomOptions{false}, bounds); };

        QMOMKIT_VERIFY_CHECK(single_position(bounds.lower - 0.5 * margin).ok());
        QMOMKIT_VERIFY_CHECK(single_position(bounds.upper + 0.5 * margin).ok());
        for (const double outside : {bounds.lower - 2.0 * margin, bounds.upper + 2.0 * margin})
        {
            const qmomkit::QmomResult rejected = single_position(outside);
            QMOMKIT_VERIFY_CHECK(!rejected.ok());
            QMOMKIT_VERIFY_CHECK(rejected.status == qmomkit::QmomStatus::not_realizable);
        }

        for (const Rule& rule : rules)
        {
            std::vector<double> nodes = rule.nodes;
            for (double& node : nodes)
            {
                node *= scale;
            }
            const qmomkit::QmomResult result =
                qmomkit::compute_qmom(qmomkit::MomentSet(raw_moments(rule.weights, nodes)), nodes.size(), qmomkit::QmomOptions{false}, bounds);
            QMOMKIT_VERIFY_CHECK(result.ok());
            QMOMKIT_VERIFY_CHECK(result.abscissas.size() == nodes.size());
        }
    }
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"fixed_bounded_qmom_accepts_nodes_inside_interval", verify_fixed_bounded_qmom_accepts_nodes_inside_interval},
            {"fixed_bounded_qmom_rejects_nodes_outside_interval_without_clipping",
             verify_fixed_bounded_qmom_rejects_nodes_outside_interval_without_clipping},
            {"adaptive_bounded_qmom_uses_highest_admissible_count", verify_adaptive_bounded_qmom_uses_highest_admissible_count},
            {"bounded_qmom_acceptance_does_not_depend_on_the_unit", verify_bounded_qmom_acceptance_does_not_depend_on_the_unit},
        });
}
