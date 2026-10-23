// One-dimensional inversion as a moment sequence approaches the edge of the
// set of moments a real distribution can produce.
//
// What is verified
//   That an inversion accepted near that edge still satisfies the properties
//   every quadrature must have - non-negative weights that sum to the zeroth
//   moment, positions in increasing order - and still reproduces the moments
//   its node count represents, whatever node count the adaptive path settles
//   on; that compute_qmom and check_realizability agree on which side of the
//   edge a sequence lies; that a fixed node count is refused, not reduced, on
//   the edge; and that a well-conditioned anchor recovers the chosen points.
//   The neighboring files cover only sequences that are
//   comfortably interior (classical rules) or exactly degenerate; these cover
//   the approach itself, where the inversion is most sensitive to rounding.
//
// How the reference is constructed
//   Chosen, not computed. Each case starts from positions and weights and the
//   moments are built from them: two positions placed a controlled distance
//   apart, a middle weight made small, or the six moments of a two-point rule
//   with m_4 moved up or down. Shrinking the distance or the weight walks the
//   sequence towards the edge. The sweeps compare moments, not positions and
//   weights: close to the edge those are ill-conditioned functions of the
//   moments (see the pass criteria), so they are compared with the chosen
//   points only at the anchor, two nodes 1e-2 apart.
//
// Error metric
//   Maximum deviation of the rebuilt moments from the supplied moments,
//   absolute or relative to |m_k| where that exceeds one. At the anchor, also
//   the absolute deviation of each position and weight from the chosen point;
//   the positions are of order one and the weights sum to m_0 = 1.
//
// Pass criteria
//   Rebuilt moments and the sum of the weights within 1.0e-12, plus the
//   quadrature properties above on every accepted result. The separation at
//   which the adaptive path starts returning fewer nodes is deliberately not
//   pinned: that threshold depends on floating-point conditioning rather than
//   on any documented contract, so asserting it would turn a numerical
//   accident into a requirement.
//
//   At the anchor, positions within 100 eps / d^2 and weights within
//   100 eps / d^3, with d = 1e-2 the separation of the close pair. Rounding the
//   moments fixes the centroid of the close pair to about eps, but splits the
//   pair into its two nodes only through the pair's third central moment, which
//   is of order d^3. The split, and with it each weight, is therefore uncertain
//   by about eps / d^3, and each position, which depends on the split through a
//   factor d, by about eps / d^2. Measured from d = 1e-1 to 1e-4, the factors
//   in front are 0.5 to 1.1 for positions and below 0.5 for weights, so the
//   tolerances leave a margin of about a hundred.
//
// Reproduce
//   ctest -R verification_qmom_.*(coalescing|adaptive_reduction|boundary)
//
// Build modes
//   Both. Checks go through QMOMKIT_VERIFY_*, not assert().

#include "../verification_helpers.hpp"
#include "qmomkit/moment_reconstruction.hpp"
#include "qmomkit/qmom.hpp"
#include "qmomkit/realizability.hpp"
#include "test_case_dispatch.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace
{

constexpr double tolerance = 1.0e-12;
// Fixed-count requests, which deliver the requested node count or fail.
const qmomkit::QmomOptions fixed_options{false};

std::vector<double> raw_moments(const std::vector<double>& weights, const std::vector<double>& abscissas, const std::size_t count)
{
    std::vector<double> moments(count, 0.0);
    for (std::size_t order = 0; order < count; ++order)
    {
        for (std::size_t node = 0; node < weights.size(); ++node)
        {
            moments[order] += weights[node] * std::pow(abscissas[node], static_cast<int>(order));
        }
    }
    return moments;
}

// Invariants that must hold for every accepted rule, whatever node count the
// inversion settled on. An accepted rule of n nodes represents the first 2n
// supplied moments; it makes no claim about the moments above that, which is
// why the reproduction check is bounded by the returned count.
void check_accepted_rule_invariants(const qmomkit::QmomResult& result, const std::vector<double>& moments, const std::size_t requested_nodes)
{
    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(!result.weights.empty());
    QMOMKIT_VERIFY_CHECK(result.weights.size() == result.abscissas.size());
    QMOMKIT_VERIFY_CHECK(result.weights.size() <= requested_nodes);
    QMOMKIT_VERIFY_CHECK(std::is_sorted(result.abscissas.begin(), result.abscissas.end()));

    double mass = 0.0;
    for (const double weight : result.weights)
    {
        QMOMKIT_VERIFY_CHECK(std::isfinite(weight));
        QMOMKIT_VERIFY_CHECK(weight >= 0.0);
        mass += weight;
    }
    QMOMKIT_VERIFY_NEAR(mass, moments[0], tolerance);

    for (const double abscissa : result.abscissas)
    {
        QMOMKIT_VERIFY_CHECK(std::isfinite(abscissa));
    }

    const std::size_t represented = 2 * result.weights.size();
    const std::vector<double> reconstructed =
        qmomkit::reconstruct_raw_moments(qmomkit::Quadrature1D{result.abscissas, result.weights}, std::min(represented, moments.size()));
    for (std::size_t order = 0; order < reconstructed.size(); ++order)
    {
        QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(reconstructed[order], moments[order], tolerance, tolerance);
    }
}

// Six moments of an exact two-node rule. For a three-node request this sequence
// sits exactly on the boundary: it is realizable, but only by two weighted
// points.
std::vector<double> two_node_boundary_moments()
{
    return raw_moments({0.4, 0.6}, {-0.75, 1.25}, 6);
}

// A three-node rule whose upper two nodes are `separation` apart. The rule stays
// exactly realizable for every positive separation and degenerates continuously
// to a two-node rule as the separation vanishes.
const std::vector<double> coalescing_weights{0.5, 0.3, 0.2};

std::vector<double> coalescing_positions(const double separation)
{
    return {-1.0, 0.5, 0.5 + separation};
}

std::vector<double> coalescing_moments(const double separation)
{
    return raw_moments(coalescing_weights, coalescing_positions(separation), 6);
}

void verify_coalescing_nodes_preserve_supplied_moments()
{
    // Anchor: a separation this large is comfortably interior, so the requested
    // count must be delivered in full. The sweep below makes no such claim.
    const double anchor_separation = 1.0e-2;
    const std::vector<double> interior = coalescing_moments(anchor_separation);
    const qmomkit::QmomResult anchored = qmomkit::compute_qmom(qmomkit::MomentSet(interior), 3, fixed_options);
    QMOMKIT_VERIFY_CHECK(anchored.weights.size() == 3);
    check_accepted_rule_invariants(anchored, interior, 3);

    // The anchor recovers the chosen points, within the rounding the pair's
    // conditioning allows (see the pass criteria in the file header).
    const double epsilon = std::numeric_limits<double>::epsilon();
    const double position_tolerance = 100.0 * epsilon / (anchor_separation * anchor_separation);
    const double weight_tolerance = position_tolerance / anchor_separation;
    const std::vector<double> chosen_positions = coalescing_positions(anchor_separation);
    for (std::size_t node = 0; node < 3; ++node)
    {
        QMOMKIT_VERIFY_NEAR(anchored.abscissas[node], chosen_positions[node], position_tolerance);
        QMOMKIT_VERIFY_NEAR(anchored.weights[node], coalescing_weights[node], weight_tolerance);
    }

    for (const double separation : {1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8})
    {
        const std::vector<double> moments = coalescing_moments(separation);
        const qmomkit::QmomResult result = qmomkit::compute_qmom(qmomkit::MomentSet(moments), 3, fixed_options);
        if (result.ok())
        {
            check_accepted_rule_invariants(result, moments, 3);
        }
        else
        {
            // The input is structurally valid and its moments are realizable by
            // construction, so a rejection here reports conditioning, not a
            // malformed request.
            QMOMKIT_VERIFY_CHECK(result.status == qmomkit::QmomStatus::not_realizable);
        }
    }
}

void verify_adaptive_reduction_preserves_represented_moments()
{
    qmomkit::QmomOptions options;
    options.adaptive = true;

    for (const double separation : {1.0e-1, 1.0e-3, 1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8, 1.0e-10})
    {
        const std::vector<double> moments = coalescing_moments(separation);
        const qmomkit::QmomResult result = qmomkit::compute_qmom(qmomkit::MomentSet(moments), 3, options);

        // Adaptive inversion accepts a reduced rule instead of failing, so every
        // point of this sweep must be accepted even where the fixed-count sweep
        // above is entitled to reject.
        check_accepted_rule_invariants(result, moments, 3);
    }
}

void verify_accepted_reduction_respects_weight_ratio()
{
    qmomkit::QmomOptions options;
    options.adaptive = true;

    // `rmin` is documented as the minimum accepted ratio between the smallest
    // and largest weight. Checking it on the returned rule pins the contract
    // without hard-coding the middle weight at which reduction happens.
    for (const double middle_weight : {1.0e-2, 1.0e-4, 1.0e-6, 1.0e-8, 1.0e-10, 1.0e-11, 1.0e-12, 1.0e-14})
    {
        const std::vector<double> moments = raw_moments({0.5, middle_weight, 0.5}, {-1.0, 0.2, 1.4}, 6);
        const qmomkit::QmomResult result = qmomkit::compute_qmom(qmomkit::MomentSet(moments), 3, options);
        check_accepted_rule_invariants(result, moments, 3);

        if (result.weights.size() > 1)
        {
            const auto extremes = std::minmax_element(result.weights.begin(), result.weights.end());
            QMOMKIT_VERIFY_CHECK(*extremes.first / *extremes.second > options.rmin);
        }
    }
}

void verify_boundary_crossing_matches_realizability()
{
    const std::vector<double> boundary = two_node_boundary_moments();

    // Perturbing the fourth moment alone moves the sequence across the boundary
    // while leaving every other moment on it. The offsets stay well outside the
    // realizability checker's tolerance band, where its classification is
    // definite rather than indeterminate.
    for (const double offset : {1.0e-1, 1.0e-3, 1.0e-6})
    {
        std::vector<double> inside = boundary;
        inside[4] += offset;
        const qmomkit::QmomResult accepted = qmomkit::compute_qmom(qmomkit::MomentSet(inside), 3, fixed_options);
        check_accepted_rule_invariants(accepted, inside, 3);
        QMOMKIT_VERIFY_CHECK(accepted.weights.size() == 3);

        const qmomkit::RealizabilityResult classification =
            qmomkit::check_realizability(qmomkit::MomentSet(inside), qmomkit::MomentSupport::real_line);
        QMOMKIT_VERIFY_CHECK(classification.status == qmomkit::RealizabilityStatus::strictly_realizable);

        std::vector<double> outside = boundary;
        outside[4] -= offset;
        const qmomkit::QmomResult rejected = qmomkit::compute_qmom(qmomkit::MomentSet(outside), 3, fixed_options);
        QMOMKIT_VERIFY_CHECK(!rejected.ok());
        QMOMKIT_VERIFY_CHECK(rejected.status == qmomkit::QmomStatus::not_realizable);

        const qmomkit::RealizabilityResult rejected_classification =
            qmomkit::check_realizability(qmomkit::MomentSet(outside), qmomkit::MomentSupport::real_line);
        QMOMKIT_VERIFY_CHECK(rejected_classification.status == qmomkit::RealizabilityStatus::not_realizable);
        QMOMKIT_VERIFY_CHECK(!rejected_classification.realizable());
    }
}

void verify_fixed_count_does_not_reduce_on_the_boundary()
{
    const std::vector<double> boundary = two_node_boundary_moments();

    // The sequence is realizable, but only by two weighted points. A fixed
    // three-node request cannot be satisfied and must say so rather than
    // quietly returning the two-node rule or a third node made of roundoff:
    // fixed-count reconstruction does not silently reduce. On the boundary the
    // third recurrence entry is zero up to roundoff, so the refusal must not
    // depend on the sign that roundoff takes in a particular unit or number
    // density.
    for (const double unit : {1.0e-3, 0.5, 1.0, 3.0, 1.0e3})
    {
        for (const double number_density : {1.0e-10, 1.0, 1.0e10})
        {
            const std::vector<double> scaled = raw_moments({0.4 * number_density, 0.6 * number_density}, {-0.75 * unit, 1.25 * unit}, 6);
            const qmomkit::QmomResult fixed = qmomkit::compute_qmom(qmomkit::MomentSet(scaled), 3, fixed_options);
            QMOMKIT_VERIFY_CHECK(!fixed.ok());
            QMOMKIT_VERIFY_CHECK(fixed.status == qmomkit::QmomStatus::not_realizable);
            QMOMKIT_VERIFY_CHECK(fixed.message.find("supports at most 2 nodes") != std::string::npos);
        }
    }

    // The same sequence is accepted once reduction is permitted, and the reduced
    // rule reproduces the four moments two nodes represent.
    qmomkit::QmomOptions options;
    options.adaptive = true;
    const qmomkit::QmomResult adaptive = qmomkit::compute_qmom(qmomkit::MomentSet(boundary), 3, options);
    check_accepted_rule_invariants(adaptive, boundary, 3);
    QMOMKIT_VERIFY_CHECK(adaptive.weights.size() < 3);

    // Two weighted points cannot also carry the fifth and sixth moments of a
    // sequence that sits on the boundary, so the accepted rule is not a
    // three-node substitute.
    const qmomkit::RealizabilityResult classification = qmomkit::check_realizability(qmomkit::MomentSet(boundary), qmomkit::MomentSupport::real_line);
    QMOMKIT_VERIFY_CHECK(classification.on_boundary());
    QMOMKIT_VERIFY_CHECK(classification.realizable());
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"coalescing_nodes_preserve_supplied_moments", verify_coalescing_nodes_preserve_supplied_moments},
            {"adaptive_reduction_preserves_represented_moments", verify_adaptive_reduction_preserves_represented_moments},
            {"accepted_reduction_respects_weight_ratio", verify_accepted_reduction_respects_weight_ratio},
            {"boundary_crossing_matches_realizability", verify_boundary_crossing_matches_realizability},
            {"fixed_count_does_not_reduce_on_the_boundary", verify_fixed_count_does_not_reduce_on_the_boundary},
        });
}
