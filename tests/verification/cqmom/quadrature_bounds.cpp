// CQMOM inversion when the caller restricts a coordinate to an interval.
//
// What is verified
//   That a branch whose positions already lie inside the requested interval is
//   accepted unchanged; that one with a position outside is rejected rather
//   than moved onto the boundary; and that a bound placed on one phase-space
//   coordinate applies to that coordinate rather than to whichever conditioning
//   level happens to hold it.
//
// How the reference is constructed
//   Chosen, not computed. Each case starts from a terminal rule written in
//   phase-space coordinate order and builds the moments from it, so the
//   expected positions and weights are known before the inversion runs. The
//   coordinate-specific case matters because phase-space order and conditioning
//   order are distinct: a bound on y must follow y even when y is conditioned
//   first.
//
// Error metric and pass criteria
//   Status must match the expected outcome; accepted positions and weights must
//   agree with the chosen rule to 1.0e-12 absolute. The rejection case checks
//   the returned position is unchanged, which is what separates rejection from
//   clipping.
//
// Reproduce
//   ctest -R "^verification_cqmom_quadrature_bounds_"
//
// Build modes
//   Both. Checks go through QMOMKIT_VERIFY_*, not assert().

#include "../verification_helpers.hpp"
#include "qmomkit/cqmom.hpp"
#include "qmomkit/moment_reconstruction.hpp"
#include "qmomkit/multi_index.hpp"
#include "test_case_dispatch.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

namespace
{

constexpr double tolerance = 1.0e-12;

std::vector<double> raw_moments(const std::vector<double>& weights, const std::vector<double>& nodes)
{
    std::vector<double> moments(4, 0.0);
    for (std::size_t order = 0; order < moments.size(); ++order)
    {
        for (std::size_t node = 0; node < nodes.size(); ++node)
        {
            moments[order] += weights[node] * std::pow(nodes[node], static_cast<int>(order));
        }
    }
    return moments;
}

qmomkit::RawMomentTensor raw_moments_from_terminals(const std::vector<std::size_t>& extents, const std::vector<double>& weights,
                                                    const std::vector<std::vector<double>>& abscissas)
{
    std::vector<double> moments(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat = 0; flat < moments.size(); ++flat)
    {
        const qmomkit::MultiIndex index = qmomkit::unflatten_multi_index(flat, extents);
        for (std::size_t node = 0; node < weights.size(); ++node)
        {
            double contribution = weights[node];
            for (std::size_t dimension = 0; dimension < abscissas.size(); ++dimension)
            {
                contribution *= std::pow(abscissas[dimension][node], static_cast<int>(index[dimension]));
            }
            moments[flat] += contribution;
        }
    }
    return {extents, moments};
}

qmomkit::CqmomCase bounded_cqmom_case(const std::vector<double>& moments, const qmomkit::AbscissaBounds interval)
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2});
    qmomkit::CqmomCase cqmom_case{
        qmomkit::extract_cqmom_moments(qmomkit::RawMomentTensor{{4}, moments}, qmomkit::create_cqmom_moment_layout(hierarchy))};
    cqmom_case.qmom_options.adaptive = false;
    cqmom_case.adaptive_branch_counts = false;
    cqmom_case.abscissa_bounds = {interval};
    cqmom_case.diagnostic_level = qmomkit::CqmomDiagnosticLevel::basic;
    return cqmom_case;
}

void verify_accepted_nodes_inside_interval()
{
    const std::vector<double> weights{0.4, 0.6};
    const std::vector<double> nodes{2.25, 4.5};
    const qmomkit::AbscissaBounds interval{2.0, 5.0};
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(bounded_cqmom_case(raw_moments(weights, nodes), interval));

    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.quadrature.weights().size() == weights.size());
    QMOMKIT_VERIFY_CHECK(result.quadrature.abscissas().size() == 1);
    QMOMKIT_VERIFY_CHECK(result.quadrature.abscissas()[0].size() == nodes.size());
    for (std::size_t node = 0; node < nodes.size(); ++node)
    {
        QMOMKIT_VERIFY_NEAR(result.quadrature.weights()[node], weights[node], tolerance);
        QMOMKIT_VERIFY_NEAR(result.quadrature.abscissas()[0][node], nodes[node], tolerance);
        QMOMKIT_VERIFY_CHECK(result.quadrature.abscissas()[0][node] >= interval.lower);
        QMOMKIT_VERIFY_CHECK(result.quadrature.abscissas()[0][node] <= interval.upper);
    }
}

void verify_rejected_nodes_outside_interval()
{
    const std::vector<double> weights{0.5, 0.5};
    const std::vector<double> nodes{0.25, 1.25};
    const qmomkit::AbscissaBounds interval{0.0, 1.0};
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(bounded_cqmom_case(raw_moments(weights, nodes), interval));

    QMOMKIT_VERIFY_CHECK(!result.ok());
    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::CqmomStatus::qmom_failed);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.size() == 1);
    const qmomkit::CqmomDiagnostics& diagnostic = result.diagnostics.front();
    QMOMKIT_VERIFY_CHECK(diagnostic.bounded_support_checked);
    QMOMKIT_VERIFY_CHECK(diagnostic.bounded_returned_nodes.size() == nodes.size());
    QMOMKIT_VERIFY_CHECK(diagnostic.bounded_returned_weights.size() == weights.size());
    for (std::size_t node = 0; node < nodes.size(); ++node)
    {
        QMOMKIT_VERIFY_NEAR(diagnostic.bounded_returned_nodes[node], nodes[node], tolerance);
        QMOMKIT_VERIFY_NEAR(diagnostic.bounded_returned_weights[node], weights[node], tolerance);
    }
    QMOMKIT_VERIFY_CHECK(diagnostic.bounded_returned_nodes[1] > interval.upper);
}

void verify_coordinate_specific_bounds_in_recursive_cqmom()
{
    // This terminal rule is prescribed in phase-space coordinate order (x, y),
    // while CQMOM conditions in the non-identity order (y, x).  The disjoint
    // supports make a conditioning-level/coordinate mix-up fail reconstruction.
    const std::vector<double> weights{0.12, 0.18, 0.28, 0.42};
    const std::vector<std::vector<double>> abscissas{{1.25, 3.5, 1.25, 3.5}, {-1.5, -1.5, -0.25, -0.25}};
    const qmomkit::AbscissaBounds x_bounds{1.0, 4.0};
    const qmomkit::AbscissaBounds y_bounds{-2.0, 0.0};
    const qmomkit::RawMomentTensor expected = raw_moments_from_terminals({4, 4}, weights, abscissas);

    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {1, 0});
    qmomkit::CqmomCase cqmom_case{qmomkit::extract_cqmom_moments(expected, qmomkit::create_cqmom_moment_layout(hierarchy))};
    cqmom_case.abscissa_bounds = {x_bounds, y_bounds};
    cqmom_case.diagnostic_level = qmomkit::CqmomDiagnosticLevel::basic;
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(cqmom_case);

    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.quadrature.total_nodes() == weights.size());
    QMOMKIT_VERIFY_CHECK(result.diagnostics.size() == 3);

    const qmomkit::RawMomentTensor reconstructed = qmomkit::reconstruct_raw_moments(result.quadrature, expected.extents);
    QMOMKIT_VERIFY_CHECK(reconstructed.values.size() == expected.values.size());
    for (std::size_t moment = 0; moment < expected.values.size(); ++moment)
    {
        QMOMKIT_VERIFY_NEAR(reconstructed.values[moment], expected.values[moment], tolerance);
    }

    bool checked_y_root = false;
    std::size_t checked_x_children = 0;
    for (const qmomkit::CqmomDiagnostics& diagnostic : result.diagnostics)
    {
        QMOMKIT_VERIFY_CHECK(diagnostic.bounded_support_checked);
        QMOMKIT_VERIFY_CHECK(diagnostic.bounded_support_interval_present);
        if (diagnostic.phase_space_dimension == 1)
        {
            QMOMKIT_VERIFY_NEAR(diagnostic.bounded_support_interval.lower, y_bounds.lower, tolerance);
            QMOMKIT_VERIFY_NEAR(diagnostic.bounded_support_interval.upper, y_bounds.upper, tolerance);
            checked_y_root = true;
        }
        else
        {
            QMOMKIT_VERIFY_CHECK(diagnostic.phase_space_dimension == 0);
            QMOMKIT_VERIFY_NEAR(diagnostic.bounded_support_interval.lower, x_bounds.lower, tolerance);
            QMOMKIT_VERIFY_NEAR(diagnostic.bounded_support_interval.upper, x_bounds.upper, tolerance);
            ++checked_x_children;
        }
    }
    QMOMKIT_VERIFY_CHECK(checked_y_root);
    QMOMKIT_VERIFY_CHECK(checked_x_children == 2);

    for (std::size_t node = 0; node < result.quadrature.total_nodes(); ++node)
    {
        const double x = result.quadrature.abscissas()[0][node];
        const double y = result.quadrature.abscissas()[1][node];
        QMOMKIT_VERIFY_CHECK(x >= x_bounds.lower && x <= x_bounds.upper);
        QMOMKIT_VERIFY_CHECK(y >= y_bounds.lower && y <= y_bounds.upper);
    }
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"accepted_nodes_inside_interval", verify_accepted_nodes_inside_interval},
            {"rejected_nodes_outside_interval", verify_rejected_nodes_outside_interval},
            {"coordinate_specific_bounds_in_recursive_cqmom", verify_coordinate_specific_bounds_in_recursive_cqmom},
        });
}
