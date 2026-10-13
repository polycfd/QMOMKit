// CQMOM reconstruction under an exact change of the unit of every coordinate
// and of the number density, on seeded random rules of weighted points.
//
// What is verified
//   No decision of the recursive reconstruction depends on the unit of a
//   coordinate or on the number density; only the absolute cutoff on m_0
//   does. Scaling coordinate d by 2^p_d and the number density by 2^q scales
//   the moment M_k by 2^(q + sum_d p_d k_d) exactly, so the result must scale
//   exactly as well: the same status and terminal node count, the abscissas
//   of coordinate d times 2^p_d and the weights times 2^q, bit for bit. This
//   covers the Wheeler inversions, the conditional Vandermonde solves with
//   every solver, the coordinate scales of adaptive inversion, and failures
//   as well as accepted rules.
//
// How the reference is constructed
//   No reference rule is needed: the reconstruction is compared with itself.
//   Rules in two and three dimensions with one to three nodes per coordinate
//   and a random conditioning order are drawn with std::mt19937_64, whose
//   output the standard fixes, mapped to [0, 1) by its upper 53 bits. Each
//   parent point has one up to the requested number of children in the next
//   coordinate of the conditioning order, clustered with spreads from 1 down to
//   1.0e-2 around centers in [-2, 2], with weights spanning three decades per
//   level, so failed fixed requests, adaptive reductions and non-identity
//   conditioning orders all occur; the case checks that they do, so it cannot
//   pass vacuously. Four random
//   scalings with p_d in [-20, 20] and q in [-30, 30] are applied to each.
//
// Error metric and pass criteria
//   Exact equality, no tolerance. Scalings that would take a moment out of the
//   normal double range are skipped, because they are not exact.
//
// Reproduce
//   ctest -R verification_cqmom_power_of_two_scaling_results_scale_exactly
//
// Build modes
//   Both. Checks go through QMOMKIT_VERIFY_*, not assert().

#include "../verification_helpers.hpp"
#include "qmomkit/cqmom.hpp"
#include "test_case_dispatch.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <random>
#include <utility>
#include <vector>

namespace
{

constexpr int rule_count = 300;

struct WeightedPoint
{
    std::vector<double> coordinates;
    double weight = 1.0;
};

double uniform(std::mt19937_64& generator)
{
    return static_cast<double>(generator() >> 11U) * 0x1.0p-53;
}

// A conditioning order drawn with a Fisher-Yates shuffle of its own, because
// std::shuffle may draw differently on another platform.
std::vector<std::size_t> draw_conditioning_order(std::mt19937_64& generator, const std::size_t dimension_count)
{
    std::vector<std::size_t> order(dimension_count);
    for (std::size_t dimension = 0; dimension < dimension_count; ++dimension)
    {
        order[dimension] = dimension;
    }
    for (std::size_t last = dimension_count - 1; last > 0; --last)
    {
        std::swap(order[last], order[static_cast<std::size_t>(generator() % (last + 1))]);
    }
    return order;
}

// A rule with the conditional structure CQMOM represents: every parent point
// has between one and the requested number of children in the next coordinate
// of the conditioning order.
std::vector<WeightedPoint> draw_rule(std::mt19937_64& generator, const std::vector<std::size_t>& node_counts,
                                     const std::vector<std::size_t>& conditioning_order)
{
    std::vector<WeightedPoint> points{WeightedPoint{std::vector<double>(node_counts.size(), 0.0), 1.0}};
    for (const std::size_t dimension : conditioning_order)
    {
        std::vector<WeightedPoint> children;
        for (const WeightedPoint& parent : points)
        {
            const std::size_t child_count = 1 + static_cast<std::size_t>(generator() % node_counts[dimension]);
            const double center = 4.0 * (uniform(generator) - 0.5);
            const double spread = std::pow(10.0, -2.0 * uniform(generator));
            for (std::size_t child = 0; child < child_count; ++child)
            {
                WeightedPoint point = parent;
                point.coordinates[dimension] = center + spread * (uniform(generator) - 0.5);
                point.weight *= std::pow(10.0, -3.0 * uniform(generator));
                children.push_back(point);
            }
        }
        points = std::move(children);
    }
    return points;
}

std::vector<double> required_moments(const std::vector<WeightedPoint>& points, const std::vector<qmomkit::MultiIndex>& indices)
{
    std::vector<double> moments(indices.size(), 0.0);
    for (std::size_t entry = 0; entry < indices.size(); ++entry)
    {
        for (const WeightedPoint& point : points)
        {
            double term = point.weight;
            for (std::size_t dimension = 0; dimension < point.coordinates.size(); ++dimension)
            {
                for (std::size_t power = 0; power < indices[entry][dimension]; ++power)
                {
                    term *= point.coordinates[dimension];
                }
            }
            moments[entry] += term;
        }
    }
    return moments;
}

void verify_results_scale_exactly()
{
    std::mt19937_64 generator(7);
    std::size_t scaled_reconstructions = 0;
    std::size_t failures = 0;
    std::size_t reductions = 0;
    std::size_t reordered_rules = 0;
    for (int trial = 0; trial < rule_count; ++trial)
    {
        std::vector<std::size_t> node_counts(2 + static_cast<std::size_t>(generator() % 2U));
        std::size_t full_node_count = 1;
        for (std::size_t& node_count : node_counts)
        {
            node_count = 1 + static_cast<std::size_t>(generator() % 3U);
            full_node_count *= node_count;
        }
        const std::vector<std::size_t> conditioning_order = draw_conditioning_order(generator, node_counts.size());
        reordered_rules += !std::is_sorted(conditioning_order.begin(), conditioning_order.end());
        const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(qmomkit::PhaseSpaceHierarchy(node_counts, conditioning_order));
        const std::vector<qmomkit::MultiIndex> indices = layout.required_indices.to_multi_indices();
        const std::vector<double> moments = required_moments(draw_rule(generator, node_counts, conditioning_order), indices);

        for (const bool adaptive : {false, true})
        {
            for (const qmomkit::VandermondeSolver solver :
                 {qmomkit::VandermondeSolver::automatic, qmomkit::VandermondeSolver::dense, qmomkit::VandermondeSolver::pivoted_qr,
                  qmomkit::VandermondeSolver::svd, qmomkit::VandermondeSolver::moore_penrose, qmomkit::VandermondeSolver::bjorck_pereyra})
            {
                qmomkit::CqmomCase cqmom_case{qmomkit::CqmomMomentData{layout, moments}};
                cqmom_case.qmom_options.adaptive = adaptive;
                cqmom_case.vandermonde_solver = solver;
                const qmomkit::CqmomResult base = qmomkit::compute_cqmom(cqmom_case);
                failures += !base.ok();
                reductions += base.ok() && base.quadrature.total_nodes() < full_node_count;

                for (int scaling = 0; scaling < 4; ++scaling)
                {
                    std::vector<int> coordinate_exponents(node_counts.size());
                    for (int& exponent : coordinate_exponents)
                    {
                        exponent = static_cast<int>(generator() % 41U) - 20;
                    }
                    const int density_exponent = static_cast<int>(generator() % 61U) - 30;

                    std::vector<double> scaled(moments.size());
                    bool exact = true;
                    for (std::size_t entry = 0; entry < moments.size(); ++entry)
                    {
                        int exponent = density_exponent;
                        for (std::size_t dimension = 0; dimension < node_counts.size(); ++dimension)
                        {
                            exponent += coordinate_exponents[dimension] * static_cast<int>(indices[entry][dimension]);
                        }
                        scaled[entry] = std::ldexp(moments[entry], exponent);
                        exact = exact && (moments[entry] == 0.0 || std::fpclassify(scaled[entry]) == FP_NORMAL);
                    }
                    if (!exact)
                    {
                        continue;
                    }

                    qmomkit::CqmomCase scaled_case = cqmom_case;
                    scaled_case.moments = qmomkit::CqmomMomentData{layout, scaled};
                    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(scaled_case);
                    ++scaled_reconstructions;
                    QMOMKIT_VERIFY_CHECK(result.status == base.status);
                    QMOMKIT_VERIFY_CHECK(result.quadrature.total_nodes() == base.quadrature.total_nodes());
                    for (std::size_t node = 0; node < base.quadrature.total_nodes(); ++node)
                    {
                        QMOMKIT_VERIFY_CHECK(result.quadrature.weights()[node] == std::ldexp(base.quadrature.weights()[node], density_exponent));
                        for (std::size_t dimension = 0; dimension < node_counts.size(); ++dimension)
                        {
                            QMOMKIT_VERIFY_CHECK(result.quadrature.abscissas()[dimension][node] ==
                                                 std::ldexp(base.quadrature.abscissas()[dimension][node], coordinate_exponents[dimension]));
                        }
                    }
                }
            }
        }
    }
    QMOMKIT_VERIFY_CHECK(scaled_reconstructions > 0);
    QMOMKIT_VERIFY_CHECK(failures > 0);
    QMOMKIT_VERIFY_CHECK(reductions > 0);
    QMOMKIT_VERIFY_CHECK(reordered_rules > 0);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(argc, argv,
                                             {
                                                 {"results_scale_exactly", verify_results_scale_exactly},
                                             });
}
