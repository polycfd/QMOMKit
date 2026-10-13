// Two properties of the Wheeler inversion that must hold for every input,
// checked on seeded random rules of weighted points, many of them on the
// boundary of moment space.
//
// What is verified
//   1. No decision depends on the unit of the coordinate or on the number
//      density; only the absolute cutoff on m_0 does. Scaling the coordinate
//      by 2^p and the number density by 2^q scales every moment exactly, so
//      the result must scale exactly as well: the same status and node count,
//      the abscissas times 2^p and the weights times 2^q, bit for bit.
//   2. The inversion never returns more nodes than the weighted points the
//      moments are built from: a fixed request for more nodes is refused, and
//      adaptive inversion returns at most as many nodes as there are points.
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
//   they do, so the verification cannot pass vacuously.
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

// Both cases draw the same rules from the same seed.
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

// A clustered rule of one to eight weighted points and the 2N moments of a
// request for N nodes, N between the point count and two more.
Rule draw_rule(std::mt19937_64& generator)
{
    Rule rule;
    rule.point_count = 1 + static_cast<std::size_t>(generator() % 8U);
    rule.requested_node_count = rule.point_count + static_cast<std::size_t>(generator() % 3U);
    const double center = 8.0 * (uniform(generator) - 0.5);
    const double spread = std::pow(10.0, -3.0 * uniform(generator));
    std::vector<double> abscissas(rule.point_count);
    std::vector<double> weights(rule.point_count);
    for (std::size_t point = 0; point < rule.point_count; ++point)
    {
        abscissas[point] = center + spread * (uniform(generator) - 0.5);
        weights[point] = std::pow(10.0, -6.0 * uniform(generator));
    }
    rule.moments.assign(2 * rule.requested_node_count, 0.0);
    for (std::size_t point = 0; point < rule.point_count; ++point)
    {
        double power = 1.0;
        for (double& moment : rule.moments)
        {
            moment += weights[point] * power;
            power *= abscissas[point];
        }
    }
    return rule;
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
        }
    }
    QMOMKIT_VERIFY_CHECK(scaled_inversions > 0);
    QMOMKIT_VERIFY_CHECK(refusals > 0);
    QMOMKIT_VERIFY_CHECK(reductions > 0);
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

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(argc, argv,
                                             {
                                                 {"power_of_two_scaling_is_exact", verify_power_of_two_scaling_is_exact},
                                                 {"no_more_nodes_than_weighted_points", verify_no_more_nodes_than_weighted_points},
                                             });
}
