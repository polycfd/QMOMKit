// Two properties of the Wheeler inversion that must hold for every input,
// checked on seeded random rules of weighted points, many of them on the
// boundary of moment space.
//
// What is verified
//   1. No decision depends on the unit of the coordinate or on the number
//      density; only the absolute cutoff on m_0 does. Scaling the coordinate
//      by 2^p and the number density by 2^q scales every moment exactly, so
//      the result must scale exactly as well: the same status and node count,
//      the abscissas times 2^p and the weights times 2^q, bit for bit. This
//      holds for the classical and the hyperbolic family: the recurrence
//      coefficients, their mean, which is the HyQMOM closure, and the rounding
//      bounds tau_k scale exactly.
//   2. The inversion never returns more nodes than the weighted points the
//      moments are built from: a fixed request for more nodes is refused, and
//      adaptive inversion returns at most as many nodes as there are points.
//      This holds for both families; the hyperbolic one requests n + 1 nodes
//      for node count n.
//
// How the reference is constructed
//   No reference rule is needed. Case 1 compares the inversion with itself,
//   case 2 compares node counts with the number of points the moments are
//   built from. Rules of one to eight points are drawn with std::mt19937_64,
//   whose output the standard fixes, mapped to [0, 1) by its upper 53 bits;
//   std::pow may round the spreads and weights differently on another
//   platform, which neither property depends on. The points are clustered,
//   with spreads from 1 down to 1.0e-3 around centers in [-4, 4], and their
//   weights span six decades. Up to two nodes more than the points are
//   requested, so many sequences lie on the boundary of moment space, and both
//   fixed-count refusals and adaptive reductions occur; the cases check that
//   they do, so the verification cannot pass vacuously. The hyperbolic family
//   is given n from 1 to the point count less one and m_0, ..., m_{2n}, so
//   that the moments lie in the interior unless the points cluster below the
//   rounding of the recurrence; its case checks that rules of n + 1 nodes and
//   classical fallbacks both occur. For case 2, the hyperbolic family is given
//   rules of one to eight points and n from 1 to the point count plus one, so
//   that it requests up to two nodes more than the points, as the classical
//   family does.
//
// Error metric and pass criteria
//   Exact equality, no tolerance. Scalings that would take a moment out of the
//   normal double range are skipped, because they are not exact.
//
// Reproduce
//   ctest -R verification_qmom_scaling_and_support
//
// Build modes
//   Both. Checks go through QMOMKIT_VERIFY_*, not assert().

#include "../verification_helpers.hpp"
#include "qmomkit/qmom.hpp"
#include "test_case_dispatch.hpp"

#include <cmath>
#include <cstddef>
#include <random>
#include <vector>

namespace
{

// The two classical cases draw the same rules from the same seed; the two
// hyperbolic cases draw their own from it.
constexpr int rule_count = 2000;
constexpr unsigned long long seed = 2026;

struct Rule
{
    std::size_t point_count = 0;
    std::size_t requested_node_count = 0;
    std::vector<double> moments;
};

double uniform(std::mt19937_64& generator)
{
    return static_cast<double>(generator() >> 11U) * 0x1.0p-53;
}

// The first moment_count moments of point_count clustered weighted points.
std::vector<double> draw_moments(std::mt19937_64& generator, const std::size_t point_count, const std::size_t moment_count)
{
    const double center = 8.0 * (uniform(generator) - 0.5);
    const double spread = std::pow(10.0, -3.0 * uniform(generator));
    std::vector<double> abscissas(point_count);
    std::vector<double> weights(point_count);
    for (std::size_t point = 0; point < point_count; ++point)
    {
        abscissas[point] = center + spread * (uniform(generator) - 0.5);
        weights[point] = std::pow(10.0, -6.0 * uniform(generator));
    }
    std::vector<double> moments(moment_count, 0.0);
    for (std::size_t point = 0; point < point_count; ++point)
    {
        double power = 1.0;
        for (double& moment : moments)
        {
            moment += weights[point] * power;
            power *= abscissas[point];
        }
    }
    return moments;
}

// A clustered rule of one to eight weighted points and the 2N moments of a
// request for N nodes, N between the point count and two more.
Rule draw_rule(std::mt19937_64& generator)
{
    Rule rule;
    rule.point_count = 1 + static_cast<std::size_t>(generator() % 8U);
    rule.requested_node_count = rule.point_count + static_cast<std::size_t>(generator() % 3U);
    rule.moments = draw_moments(generator, rule.point_count, 2 * rule.requested_node_count);
    return rule;
}

// A clustered rule of two to eight weighted points and the 2n + 1 moments of a
// hyperbolic request with node count n, n from 1 to the point count less one.
Rule draw_hyperbolic_rule(std::mt19937_64& generator)
{
    Rule rule;
    rule.point_count = 2 + static_cast<std::size_t>(generator() % 7U);
    rule.requested_node_count = 1 + static_cast<std::size_t>(generator() % (rule.point_count - 1));
    rule.moments = draw_moments(generator, rule.point_count, 2 * rule.requested_node_count + 1);
    return rule;
}

// A clustered rule of one to eight weighted points and the 2n + 1 moments of a
// hyperbolic request with node count n, n from 1 to the point count plus one,
// so that the n + 1 requested nodes exceed the points by up to two.
Rule draw_hyperbolic_boundary_rule(std::mt19937_64& generator)
{
    Rule rule;
    rule.point_count = 1 + static_cast<std::size_t>(generator() % 8U);
    rule.requested_node_count = 1 + static_cast<std::size_t>(generator() % (rule.point_count + 1));
    rule.moments = draw_moments(generator, rule.point_count, 2 * rule.requested_node_count + 1);
    return rule;
}

// Inverts the rule at every representable power-of-two scaling and requires
// the result of base, scaled exactly; returns the number of scaled inversions.
std::size_t verify_scaled_inversions(const Rule& rule, const qmomkit::QmomOptions& options, const qmomkit::QmomResult& base)
{
    std::size_t scaled_inversions = 0;
    for (const int coordinate_exponent : {-40, -7, 13, 60})
    {
        for (const int density_exponent : {-50, 0, 30})
        {
            std::vector<double> scaled(rule.moments.size());
            bool exact = true;
            for (std::size_t order = 0; order < scaled.size(); ++order)
            {
                scaled[order] = std::ldexp(rule.moments[order], density_exponent + coordinate_exponent * static_cast<int>(order));
                exact = exact && (rule.moments[order] == 0.0 || std::fpclassify(scaled[order]) == FP_NORMAL);
            }
            if (!exact)
            {
                continue;
            }
            const qmomkit::QmomResult result = qmomkit::compute_qmom(qmomkit::MomentSet{scaled}, rule.requested_node_count, options);
            ++scaled_inversions;
            QMOMKIT_VERIFY_CHECK(result.status == base.status);
            QMOMKIT_VERIFY_CHECK(result.weights.size() == base.weights.size());
            for (std::size_t node = 0; node < base.weights.size(); ++node)
            {
                QMOMKIT_VERIFY_CHECK(result.abscissas[node] == std::ldexp(base.abscissas[node], coordinate_exponent));
                QMOMKIT_VERIFY_CHECK(result.weights[node] == std::ldexp(base.weights[node], density_exponent));
            }
        }
    }
    return scaled_inversions;
}

void verify_power_of_two_scaling_is_exact()
{
    std::mt19937_64 generator(seed);
    std::size_t scaled_inversions = 0;
    std::size_t refusals = 0;
    std::size_t reductions = 0;
    for (int trial = 0; trial < rule_count; ++trial)
    {
        const Rule rule = draw_rule(generator);
        for (const bool adaptive : {false, true})
        {
            qmomkit::QmomOptions options;
            options.adaptive = adaptive;
            const qmomkit::QmomResult base = qmomkit::compute_qmom(qmomkit::MomentSet{rule.moments}, rule.requested_node_count, options);
            refusals += !adaptive && !base.ok();
            reductions += base.ok() && base.weights.size() < rule.requested_node_count;
            scaled_inversions += verify_scaled_inversions(rule, options, base);
        }
    }
    QMOMKIT_VERIFY_CHECK(scaled_inversions > 0);
    QMOMKIT_VERIFY_CHECK(refusals > 0);
    QMOMKIT_VERIFY_CHECK(reductions > 0);
}

void verify_hyperbolic_power_of_two_scaling_is_exact()
{
    std::mt19937_64 generator(seed);
    std::size_t scaled_inversions = 0;
    std::size_t hyperbolic_rules = 0;
    std::size_t classical_fallbacks = 0;
    for (int trial = 0; trial < rule_count; ++trial)
    {
        const Rule rule = draw_hyperbolic_rule(generator);
        for (const bool adaptive : {false, true})
        {
            qmomkit::QmomOptions options;
            options.adaptive = adaptive;
            options.family = qmomkit::QmomFamily::hyperbolic;
            const qmomkit::QmomResult base = qmomkit::compute_qmom(qmomkit::MomentSet{rule.moments}, rule.requested_node_count, options);
            hyperbolic_rules += base.ok() && base.weights.size() == rule.requested_node_count + 1;
            classical_fallbacks += base.ok() && base.weights.size() <= rule.requested_node_count;
            scaled_inversions += verify_scaled_inversions(rule, options, base);
        }
    }
    QMOMKIT_VERIFY_CHECK(scaled_inversions > 0);
    QMOMKIT_VERIFY_CHECK(hyperbolic_rules > 0);
    QMOMKIT_VERIFY_CHECK(classical_fallbacks > 0);
}

void verify_no_more_nodes_than_weighted_points()
{
    std::mt19937_64 generator(seed);
    std::size_t boundary_requests = 0;
    for (int trial = 0; trial < rule_count; ++trial)
    {
        const Rule rule = draw_rule(generator);
        boundary_requests += rule.requested_node_count > rule.point_count;

        qmomkit::QmomOptions fixed_options;
        fixed_options.adaptive = false;
        const qmomkit::QmomResult fixed = qmomkit::compute_qmom(qmomkit::MomentSet{rule.moments}, rule.requested_node_count, fixed_options);
        if (rule.requested_node_count > rule.point_count)
        {
            QMOMKIT_VERIFY_CHECK(fixed.status == qmomkit::QmomStatus::not_realizable);
        }

        const qmomkit::QmomResult adaptive =
            qmomkit::compute_qmom(qmomkit::MomentSet{rule.moments}, rule.requested_node_count, qmomkit::QmomOptions{});
        QMOMKIT_VERIFY_CHECK(adaptive.ok());
        QMOMKIT_VERIFY_CHECK(adaptive.weights.size() <= rule.point_count);
    }
    QMOMKIT_VERIFY_CHECK(boundary_requests > 0);
}

void verify_hyperbolic_no_more_nodes_than_weighted_points()
{
    std::mt19937_64 generator(seed);
    std::size_t boundary_requests = 0;
    for (int trial = 0; trial < rule_count; ++trial)
    {
        const Rule rule = draw_hyperbolic_boundary_rule(generator);
        const bool beyond_points = rule.requested_node_count + 1 > rule.point_count;
        boundary_requests += beyond_points;

        qmomkit::QmomOptions fixed_options;
        fixed_options.adaptive = false;
        fixed_options.family = qmomkit::QmomFamily::hyperbolic;
        const qmomkit::QmomResult fixed = qmomkit::compute_qmom(qmomkit::MomentSet{rule.moments}, rule.requested_node_count, fixed_options);
        if (beyond_points)
        {
            QMOMKIT_VERIFY_CHECK(fixed.status == qmomkit::QmomStatus::not_realizable);
        }

        qmomkit::QmomOptions adaptive_options;
        adaptive_options.family = qmomkit::QmomFamily::hyperbolic;
        const qmomkit::QmomResult adaptive = qmomkit::compute_qmom(qmomkit::MomentSet{rule.moments}, rule.requested_node_count, adaptive_options);
        QMOMKIT_VERIFY_CHECK(adaptive.ok());
        QMOMKIT_VERIFY_CHECK(adaptive.weights.size() <= rule.point_count);
    }
    QMOMKIT_VERIFY_CHECK(boundary_requests > 0);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"power_of_two_scaling_is_exact", verify_power_of_two_scaling_is_exact},
            {"hyperbolic_power_of_two_scaling_is_exact", verify_hyperbolic_power_of_two_scaling_is_exact},
            {"no_more_nodes_than_weighted_points", verify_no_more_nodes_than_weighted_points},
            {"hyperbolic_no_more_nodes_than_weighted_points", verify_hyperbolic_no_more_nodes_than_weighted_points},
        });
}
