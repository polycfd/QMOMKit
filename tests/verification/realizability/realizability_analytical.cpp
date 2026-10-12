// Deciding whether a moment sequence could have come from a real distribution.
//
// What is verified
//   That sequences known to be realizable are accepted and sequences known not
//   to be are rejected, and that in both cases the diagnostics identify which
//   order and which matrix decided it, with a pivot whose sign and magnitude
//   are consistent with that verdict.
//
// How the reference is constructed
//   By construction on both sides. Realizable cases are built from actual
//   weighted points, so a distribution producing those moments demonstrably
//   exists. Unrealizable cases are built by taking such a sequence and
//   violating a specific requirement - a negative zeroth moment, a negative
//   variance, a Hankel determinant driven through zero - so the expected
//   verdict and the order at which it should be reached are both known before
//   the check runs.
//
//   The remaining cases rest on two standard facts. A measure of n points
//   inside the support is strictly realizable up to order 2n - 1 and lies on
//   the boundary from order 2n on, where the (n + 1) x (n + 1) Hankel matrix
//   becomes singular. A point at an end of the support contributes nothing to
//   the localizing matrix for that end, so the sequence is on the boundary once
//   that matrix is larger than the number of the other points; a single point
//   at an end is on the boundary from order one. And on an interval [a, b]
//   every distribution satisfies the
//   Bhatia-Davis bound variance <= (b - mean)(mean - a), so a larger variance
//   is not realizable.
//
//   Realizability is invariant under an affine change of coordinate, so a
//   sequence built from a location c and a spread s has the classification of
//   the one at c = 0 and s = 1. Double-precision input resolves it only up to
//   a shift c / s that falls with the order; each case requires its verdict up
//   to a shift at least three times below the first one, on a grid of 0 to
//   1e8, at which any spread from 1e-9 to 1e9 gave a different answer.
//
//   Weighted points on a lattice of 1/4 below 16 in magnitude, with weights
//   that are multiples of 1/8, have raw moments up to order 8 that are exact
//   in double precision, so such sequences are realizable as given, and on the
//   boundary from order 2k on for k points.
//
// Error metric and pass criteria
//   Verdict and reported order must match exactly. The reported pivot must
//   exceed the tolerance actually used for a strictly realizable sequence, and
//   fall below its negative for a rejected one, which is what distinguishes a
//   decision from a coin flip near the boundary. The classification runs at
//   its default relative tolerance of 1.0e-12, except in the near-boundary
//   cases, which set 1.0e-8 so that moving m_4 by 1e-10 falls inside the
//   tolerance band and moving it by 1e-6 falls outside. Beyond its resolved
//   shift, a case may only be boundary or numerically indeterminate instead of
//   its own verdict. An exact point sequence may never be rejected, nor called
//   strictly realizable on the boundary.
//
// Reproduce
//   ctest -R verification_realizability
//
// Build modes
//   Both. Failures are reported through an exit code, not assert().

#include "../verification_helpers.hpp"
#include "test_case_dispatch.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <qmomkit/distributions.hpp>
#include <qmomkit/realizability.hpp>
#include <string>
#include <utility>
#include <vector>

namespace
{

using qmomkit::BetaDistribution;
using qmomkit::BoundedInterval;
using qmomkit::GammaDistribution;
using qmomkit::LogNormalDistribution;
using qmomkit::MomentSet;
using qmomkit::MomentSupport;
using qmomkit::RealizabilityOptions;
using qmomkit::RealizabilityStatus;

bool contains_text(const std::string& text, const std::string& fragment)
{
    return text.find(fragment) != std::string::npos;
}

[[nodiscard]] MomentSet discrete_moments(const std::vector<double>& weights, const std::vector<double>& nodes, const std::size_t maximum_order)
{
    std::vector<double> moments(maximum_order + 1, 0.0);
    for (std::size_t order = 0; order <= maximum_order; ++order)
    {
        for (std::size_t node = 0; node < nodes.size(); ++node)
        {
            moments[order] += weights[node] * std::pow(nodes[node], static_cast<double>(order));
        }
    }
    return MomentSet(std::move(moments));
}

[[nodiscard]] MomentSet standard_normal_moments(const std::size_t maximum_order)
{
    std::vector<double> moments(maximum_order + 1, 0.0);
    moments[0] = 1.0;
    if (maximum_order >= 2)
    {
        for (std::size_t order = 2; order <= maximum_order; ++order)
        {
            moments[order] = static_cast<double>(order - 1) * moments[order - 2];
        }
    }
    return MomentSet(std::move(moments));
}

[[nodiscard]] MomentSet scaled_moments(const MomentSet& moments, const double scale)
{
    std::vector<double> scaled = moments.values();
    for (double& value : scaled)
    {
        value *= scale;
    }
    return MomentSet(std::move(scaled));
}

[[nodiscard]] MomentSet uniform_01_moments(const std::size_t maximum_order)
{
    std::vector<double> moments(maximum_order + 1, 0.0);
    for (std::size_t order = 0; order <= maximum_order; ++order)
    {
        moments[order] = 1.0 / static_cast<double>(order + 1);
    }
    return MomentSet(std::move(moments));
}

// Moments of the uniform distribution on [lower, upper], with unit number density.
[[nodiscard]] MomentSet uniform_moments(const double lower, const double upper, const std::size_t maximum_order)
{
    std::vector<double> moments(maximum_order + 1, 0.0);
    for (std::size_t order = 0; order <= maximum_order; ++order)
    {
        const double k = static_cast<double>(order);
        moments[order] = (std::pow(upper, k + 1.0) - std::pow(lower, k + 1.0)) / ((k + 1.0) * (upper - lower));
    }
    return MomentSet(std::move(moments));
}

void expect_status(const MomentSet& moments, const MomentSupport support, const RealizabilityStatus status, const std::string& name)
{
    const auto result = qmomkit::check_realizability(moments, support);
    if (result.status != status)
    {
        std::cerr << name << " expected status " << static_cast<int>(status) << " but got " << static_cast<int>(result.status) << ": "
                  << result.message << '\n';
    }
    QMOMKIT_VERIFY_CHECK(result.status == status);
}

void expect_status(const MomentSet& moments, const BoundedInterval interval, const RealizabilityStatus status, const std::string& name)
{
    const auto result = qmomkit::check_realizability(moments, interval);
    if (result.status != status)
    {
        std::cerr << name << " expected status " << static_cast<int>(status) << " but got " << static_cast<int>(result.status) << ": "
                  << result.message << '\n';
    }
    QMOMKIT_VERIFY_CHECK(result.status == status);
}

void verify_success_diagnostics(const qmomkit::RealizabilityResult& result, const RealizabilityStatus expected_status,
                                const std::size_t expected_checked_order, const std::string& support_fragment, const std::string& name)
{
    // A successful public realizability result should still tell users what
    // was checked. The "failed_matrix" field is also used as the limiting
    // public matrix diagnostic on successful classifications.
    QMOMKIT_VERIFY_CHECK(result.status == expected_status);
    QMOMKIT_VERIFY_CHECK(result.realizable());
    QMOMKIT_VERIFY_CHECK(result.on_boundary() == (expected_status == RealizabilityStatus::boundary));
    QMOMKIT_VERIFY_CHECK(result.checked_order == expected_checked_order);
    QMOMKIT_VERIFY_CHECK(!result.failed_matrix.empty());
    QMOMKIT_VERIFY_CHECK(result.failed_matrix_order <= expected_checked_order);
    QMOMKIT_VERIFY_CHECK(std::isfinite(result.minimum_pivot));
    QMOMKIT_VERIFY_CHECK(result.tolerance_used > 0.0);
    QMOMKIT_VERIFY_CHECK(!result.message.empty());
    QMOMKIT_VERIFY_CHECK(contains_text(result.message, support_fragment));
    if (expected_status == RealizabilityStatus::strictly_realizable && result.minimum_pivot <= result.tolerance_used)
    {
        std::cerr << name << " did not report a strictly positive limiting pivot. pivot=" << result.minimum_pivot
                  << " tolerance=" << result.tolerance_used << '\n';
        QMOMKIT_VERIFY_CHECK(result.minimum_pivot > result.tolerance_used);
    }
}

void verify_failure_diagnostics(const qmomkit::RealizabilityResult& result, const std::size_t expected_checked_order,
                                const std::string& expected_matrix, const std::size_t expected_matrix_order,
                                const std::string& expected_message_fragment)
{
    // Non-realizability diagnostics should identify the public mathematical
    // condition that failed, not merely report a final status code.
    QMOMKIT_VERIFY_CHECK(result.status == RealizabilityStatus::not_realizable);
    QMOMKIT_VERIFY_CHECK(!result.realizable());
    QMOMKIT_VERIFY_CHECK(!result.on_boundary());
    QMOMKIT_VERIFY_CHECK(result.checked_order == expected_checked_order);
    QMOMKIT_VERIFY_CHECK(result.failed_matrix == expected_matrix);
    QMOMKIT_VERIFY_CHECK(result.failed_matrix_order == expected_matrix_order);
    QMOMKIT_VERIFY_CHECK(result.minimum_pivot < -result.tolerance_used || expected_matrix == "m0");
    QMOMKIT_VERIFY_CHECK(result.tolerance_used > 0.0);
    QMOMKIT_VERIFY_CHECK(contains_text(result.message, expected_message_fragment));
}

void verify_diagnostic_success_cases()
{
    // Hamburger: the standard normal moments are strictly realizable on the
    // full real line. The moment matrix diagnostics should report a successful
    // positive limiting check through order five.
    verify_success_diagnostics(qmomkit::check_realizability(standard_normal_moments(5), MomentSupport::real_line),
                               RealizabilityStatus::strictly_realizable, 5, "real_line", "standard normal Hamburger diagnostics");

    // Stieltjes: gamma moments come from a positive density on [0, infinity).
    const auto gamma = qmomkit::compute_raw_moments(GammaDistribution{2.0, 1.0}, 5);
    QMOMKIT_VERIFY_CHECK(gamma.ok());
    verify_success_diagnostics(qmomkit::check_realizability(gamma.moments, MomentSupport::nonnegative), RealizabilityStatus::strictly_realizable, 5,
                               "nonnegative", "gamma Stieltjes diagnostics");

    // Hausdorff: beta moments come from a positive density on [0,1]. The
    // bounded-interval overload should preserve interval-aware diagnostics.
    const auto beta = qmomkit::compute_raw_moments(BetaDistribution{2.0, 3.0}, 5);
    QMOMKIT_VERIFY_CHECK(beta.ok());
    verify_success_diagnostics(qmomkit::check_realizability(beta.moments, BoundedInterval{0.0, 1.0}), RealizabilityStatus::strictly_realizable, 5,
                               "bounded_interval", "beta Hausdorff diagnostics");
}

void verify_diagnostic_failure_cases()
{
    // A negative number density is support-independent: no nonnegative
    // measure can have m0 < 0, so the public failure should name m0 directly.
    const auto negative_mass = qmomkit::check_realizability(MomentSet{-1.0, 0.0, 1.0}, MomentSupport::real_line);
    verify_failure_diagnostics(negative_mass, 2, "m0", 0, "m0 is negative");
    QMOMKIT_VERIFY_NEAR(negative_mass.minimum_pivot, -1.0, 1.0e-15);

    // Hamburger positivity: m2 < m1^2/m0 gives a negative variance. The
    // failing public condition is the first nontrivial Hankel matrix H_1.
    verify_failure_diagnostics(qmomkit::check_realizability(MomentSet{1.0, 0.0, -1.0}, MomentSupport::real_line), 2, "H_1", 1, "negative eigenvalue");

    // Stieltjes localizing positivity: a negative first moment cannot come
    // from support [0, infinity), even though the ordinary H_1 matrix is
    // positive semidefinite. The shifted matrix S_0 should identify the
    // support-specific failure.
    verify_failure_diagnostics(qmomkit::check_realizability(MomentSet{1.0, -1.0, 2.0}, MomentSupport::nonnegative), 2, "S_0", 0,
                               "negative eigenvalue");

    // Hausdorff upper-bound support: m1 > m0 on [0,1] violates the upper
    // localizing matrix for an order-one moment sequence.
    verify_failure_diagnostics(qmomkit::check_realizability(MomentSet{1.0, 1.2}, BoundedInterval{0.0, 1.0}), 1, "L_upper_0", 0,
                               "negative eigenvalue");
}

void verify_boundary_diagnostics()
{
    // The zero measure is an analytical boundary case. The public diagnostic
    // identifies the zero-measure boundary directly without reporting a matrix
    // failure.
    const auto zero = qmomkit::check_realizability(MomentSet{0.0, 0.0, 0.0}, MomentSupport::real_line);
    QMOMKIT_VERIFY_CHECK(zero.status == RealizabilityStatus::boundary);
    QMOMKIT_VERIFY_CHECK(zero.realizable());
    QMOMKIT_VERIFY_CHECK(zero.on_boundary());
    QMOMKIT_VERIFY_CHECK(zero.checked_order == 2);
    QMOMKIT_VERIFY_CHECK(zero.failed_matrix == "zero_measure");
    QMOMKIT_VERIFY_CHECK(zero.failed_matrix_order == 0);
    QMOMKIT_VERIFY_NEAR(zero.minimum_pivot, 0.0, 0.0);
    QMOMKIT_VERIFY_CHECK(zero.tolerance_used > 0.0);
    QMOMKIT_VERIFY_CHECK(contains_text(zero.message, "zero measure"));

    // Tolerance-controlled near-boundary cases. Two points at -1 and 1 put H_2
    // exactly on the boundary at order 4; lowering m_4 by 1e-10 keeps its zero
    // eigenvalue within a loose tolerance, so its sign is indeterminate. A
    // variance of -1e-6 or 1e-6 is decided at its own scale. The tolerance used
    // refers to the standardized matrices, against which the limiting
    // eigenvalue is compared.
    RealizabilityOptions loose;
    loose.tolerance = 1.0e-8;

    const auto tiny_negative = qmomkit::check_realizability(MomentSet{1.0, 0.0, 1.0, 0.0, 1.0 - 1.0e-10}, MomentSupport::real_line, loose);
    QMOMKIT_VERIFY_CHECK(tiny_negative.status == RealizabilityStatus::numerically_indeterminate);
    QMOMKIT_VERIFY_CHECK(tiny_negative.checked_order == 4);
    QMOMKIT_VERIFY_CHECK(tiny_negative.failed_matrix == "H_2");
    QMOMKIT_VERIFY_CHECK(tiny_negative.failed_matrix_order == 2);
    QMOMKIT_VERIFY_CHECK(tiny_negative.minimum_pivot < 0.0);
    QMOMKIT_VERIFY_CHECK(std::abs(tiny_negative.minimum_pivot) < tiny_negative.tolerance_used);
    QMOMKIT_VERIFY_CHECK(contains_text(tiny_negative.message, "indeterminate"));

    const auto clear_negative = qmomkit::check_realizability(MomentSet{1.0, 0.0, -1.0e-6}, MomentSupport::real_line, loose);
    verify_failure_diagnostics(clear_negative, 2, "H_1", 1, "negative eigenvalue");

    const auto clear_positive = qmomkit::check_realizability(MomentSet{1.0, 0.0, 1.0e-6}, MomentSupport::real_line, loose);
    QMOMKIT_VERIFY_CHECK(clear_positive.status == RealizabilityStatus::strictly_realizable);
    QMOMKIT_VERIFY_CHECK(clear_positive.checked_order == 2);
    QMOMKIT_VERIFY_CHECK(clear_positive.failed_matrix == "H_1");
    QMOMKIT_VERIFY_CHECK(clear_positive.minimum_pivot > clear_positive.tolerance_used);
}

void verify_hamburger()
{
    // Hamburger checks cover moments on the full real line: strictly
    // realizable normal moments, boundary weighted-point measures, and basic
    // positivity/range violations.
    expect_status(standard_normal_moments(5), MomentSupport::real_line, RealizabilityStatus::strictly_realizable, "normal order 5");
    expect_status(discrete_moments({2.0}, {-1.5}, 5), MomentSupport::real_line, RealizabilityStatus::boundary, "single Dirac order 5");
    expect_status(discrete_moments({0.25, 0.75}, {-1.0, 2.0}, 5), MomentSupport::real_line, RealizabilityStatus::boundary, "two-point order 5");
    expect_status(MomentSet{1.0, 0.0, -1.0}, MomentSupport::real_line, RealizabilityStatus::not_realizable, "variance violation");
    expect_status(MomentSet{1.0, 0.0, 0.0, 1.0}, MomentSupport::real_line, RealizabilityStatus::not_realizable, "odd range violation");
}

void verify_stieltjes()
{
    // Stieltjes checks add nonnegative-support constraints. Gamma and
    // lognormal moments are analytical positive-density references, while the
    // failure cases isolate ordinary and shifted-matrix violations.
    const auto gamma = qmomkit::compute_raw_moments(GammaDistribution{2.0, 1.0}, 5);
    QMOMKIT_VERIFY_CHECK(gamma.ok());
    expect_status(gamma.moments, MomentSupport::nonnegative, RealizabilityStatus::strictly_realizable, "gamma order 5");

    const auto lognormal = qmomkit::compute_raw_moments(LogNormalDistribution{0.0, 0.25}, 4);
    QMOMKIT_VERIFY_CHECK(lognormal.ok());
    expect_status(lognormal.moments, MomentSupport::nonnegative, RealizabilityStatus::strictly_realizable, "lognormal order 4");

    expect_status(discrete_moments({1.0}, {0.0}, 5), MomentSupport::nonnegative, RealizabilityStatus::boundary, "Dirac zero");
    expect_status(discrete_moments({1.0}, {2.0}, 5), MomentSupport::nonnegative, RealizabilityStatus::boundary, "positive atom");
    expect_status(MomentSet{1.0, -1.0, 2.0}, MomentSupport::nonnegative, RealizabilityStatus::not_realizable, "negative first moment");
    expect_status(MomentSet{1.0, 1.0, 1.0, -1.0}, MomentSupport::nonnegative, RealizabilityStatus::not_realizable, "shifted matrix failure");
    expect_status(MomentSet{1.0, 0.0, 1.0}, MomentSupport::nonnegative, RealizabilityStatus::not_realizable, "singular shifted range failure");
    expect_status(MomentSet{1.0, 0.0, 0.0, 1.0}, MomentSupport::nonnegative, RealizabilityStatus::not_realizable, "odd ordinary range failure");
}

void verify_hausdorff()
{
    // Hausdorff checks use the bounded interval [0,1]. Beta and uniform
    // moments are analytical references, and the boundary/failure cases protect
    // lower- and upper-localizing behavior.
    const BoundedInterval unit{0.0, 1.0};
    const auto beta = qmomkit::compute_raw_moments(BetaDistribution{2.0, 3.0}, 5);
    QMOMKIT_VERIFY_CHECK(beta.ok());
    expect_status(beta.moments, unit, RealizabilityStatus::strictly_realizable, "beta order 5");
    expect_status(uniform_01_moments(4), unit, RealizabilityStatus::strictly_realizable, "uniform order 4");

    expect_status(discrete_moments({1.0}, {0.0}, 5), unit, RealizabilityStatus::boundary, "lower bound atom");
    expect_status(discrete_moments({1.0}, {1.0}, 5), unit, RealizabilityStatus::boundary, "upper bound atom");
    expect_status(discrete_moments({0.5, 0.5}, {0.0, 0.5}, 5), unit, RealizabilityStatus::boundary, "two-point with boundary atom");
    expect_status(MomentSet{1.0, 1.2}, unit, RealizabilityStatus::not_realizable, "mean outside interval");
    expect_status(MomentSet{1.0, 0.0, 1.0}, unit, RealizabilityStatus::not_realizable, "interval localizing failure");
    expect_status(MomentSet{1.0, 0.25, 0.125, 0.0625}, unit, RealizabilityStatus::boundary, "odd order boundary");
}

void verify_hausdorff_on_other_intervals()
{
    // The [0,1] cases above, on an interval containing zero and on one away
    // from it. The variance cases sit 0.25 above and below the Bhatia-Davis
    // bound for a mean at the center.
    for (const BoundedInterval interval : {BoundedInterval{-1.0, 2.0}, BoundedInterval{2.0, 5.0}})
    {
        const double a = interval.lower;
        const double b = interval.upper;
        const double centre = 0.5 * (a + b);
        const double largest_variance = (b - centre) * (centre - a);

        expect_status(uniform_moments(a, b, 4), interval, RealizabilityStatus::strictly_realizable, "uniform order 4");
        expect_status(uniform_moments(a, b, 5), interval, RealizabilityStatus::strictly_realizable, "uniform order 5");
        expect_status(discrete_moments({1.5}, {centre}, 1), interval, RealizabilityStatus::strictly_realizable, "interior atom order 1");
        expect_status(discrete_moments({1.5}, {centre}, 2), interval, RealizabilityStatus::boundary, "interior atom order 2");
        expect_status(discrete_moments({0.3, 0.7}, {a + 0.5, b - 1.0}, 3), interval, RealizabilityStatus::strictly_realizable,
                      "two interior atoms order 3");
        expect_status(discrete_moments({0.3, 0.7}, {a + 0.5, b - 1.0}, 4), interval, RealizabilityStatus::boundary, "two interior atoms order 4");
        expect_status(discrete_moments({1.0}, {a}, 3), interval, RealizabilityStatus::boundary, "lower end atom");
        expect_status(discrete_moments({1.0}, {b}, 3), interval, RealizabilityStatus::boundary, "upper end atom");
        expect_status(MomentSet{1.0, b + 0.2}, interval, RealizabilityStatus::not_realizable, "mean above the interval");
        expect_status(MomentSet{1.0, a - 0.2}, interval, RealizabilityStatus::not_realizable, "mean below the interval");
        expect_status(MomentSet{1.0, centre, centre * centre + largest_variance - 0.25}, interval, RealizabilityStatus::strictly_realizable,
                      "variance below the Bhatia-Davis bound");
        expect_status(MomentSet{1.0, centre, centre * centre + largest_variance + 0.25}, interval, RealizabilityStatus::not_realizable,
                      "variance above the Bhatia-Davis bound");
    }
}

void verify_even_order_boundaries()
{
    // n interior points: strictly realizable at order 2n - 1, boundary at 2n.
    expect_status(discrete_moments({0.4, 0.6}, {-0.75, 1.25}, 3), MomentSupport::real_line, RealizabilityStatus::strictly_realizable,
                  "two points order 3");
    expect_status(discrete_moments({0.4, 0.6}, {-0.75, 1.25}, 4), MomentSupport::real_line, RealizabilityStatus::boundary, "two points order 4");
    expect_status(discrete_moments({0.25, 0.5, 0.25}, {-1.0, 0.5, 2.0}, 5), MomentSupport::real_line, RealizabilityStatus::strictly_realizable,
                  "three points order 5");
    expect_status(discrete_moments({0.25, 0.5, 0.25}, {-1.0, 0.5, 2.0}, 6), MomentSupport::real_line, RealizabilityStatus::boundary,
                  "three points order 6");
    expect_status(discrete_moments({0.4, 0.6}, {0.5, 2.0}, 3), MomentSupport::nonnegative, RealizabilityStatus::strictly_realizable,
                  "two positive points order 3");
    expect_status(discrete_moments({0.4, 0.6}, {0.5, 2.0}, 4), MomentSupport::nonnegative, RealizabilityStatus::boundary,
                  "two positive points order 4");
    expect_status(discrete_moments({0.4, 0.6}, {0.25, 0.75}, 4), BoundedInterval{0.0, 1.0}, RealizabilityStatus::boundary,
                  "two interior points order 4");

    // The localizing matrix x H does not see the point at zero, so it has rank
    // one here and is singular from order 3, where it is 2 x 2.
    expect_status(discrete_moments({0.4, 0.6}, {0.0, 2.0}, 3), MomentSupport::nonnegative, RealizabilityStatus::boundary, "zero point order 3");
    expect_status(discrete_moments({0.4, 0.6}, {0.0, 2.0}, 4), MomentSupport::nonnegative, RealizabilityStatus::boundary, "zero point order 4");
}

// Raw moments of the normal distribution with the given mean and standard
// deviation, from m_k = mean m_(k-1) + (k - 1) deviation^2 m_(k-2).
[[nodiscard]] MomentSet normal_moments(const double mean, const double deviation, const std::size_t maximum_order)
{
    std::vector<double> moments{1.0, mean};
    for (std::size_t order = 2; order <= maximum_order; ++order)
    {
        moments.push_back(mean * moments[order - 1] + static_cast<double>(order - 1) * deviation * deviation * moments[order - 2]);
    }
    return MomentSet(std::move(moments));
}

// Adds offset * spread^k to m_k, a fixed change of the standardized moment.
[[nodiscard]] MomentSet with_offset(const MomentSet& moments, const std::size_t order, const double offset, const double spread)
{
    std::vector<double> values = moments.values();
    values[order] += offset * std::pow(spread, static_cast<double>(order));
    return MomentSet(std::move(values));
}

// A sequence built from a location c and a spread s, whose classification
// does not depend on either. `resolved_shift` is the largest c / s at which
// double-precision input still resolves it; `positive_shift` marks sequences
// whose points must stay to the right of the origin.
struct InvarianceCase
{
    std::string name;
    std::function<qmomkit::RealizabilityResult(double, double)> check;
    RealizabilityStatus expected = RealizabilityStatus::strictly_realizable;
    double resolved_shift = 0.0;
    bool positive_shift = false;
};

[[nodiscard]] std::vector<InvarianceCase> invariance_cases()
{
    const auto real_line = [](const MomentSet& moments) { return qmomkit::check_realizability(moments, MomentSupport::real_line); };
    const auto nonnegative = [](const MomentSet& moments) { return qmomkit::check_realizability(moments, MomentSupport::nonnegative); };
    const auto three_points = [](const double c, const double s, const std::size_t order)
    { return discrete_moments({0.25, 0.5, 0.25}, {c - s, c, c + s}, order); };
    const auto two_points = [](const double c, const double s, const std::size_t order)
    { return discrete_moments({0.5, 0.5}, {c - s, c + s}, order); };
    const auto four_points = [](const double c, const double s, const std::size_t order)
    { return discrete_moments({0.1, 0.4, 0.4, 0.1}, {c - 1.5 * s, c - 0.5 * s, c + 0.5 * s, c + 1.5 * s}, order); };
    constexpr auto strict = RealizabilityStatus::strictly_realizable;
    constexpr auto boundary = RealizabilityStatus::boundary;
    constexpr auto not_realizable = RealizabilityStatus::not_realizable;

    return {
        {"real line, three points, order 5", [=](double c, double s) { return real_line(three_points(c, s, 5)); }, strict, 1.0e3},
        {"real line, three points, order 6", [=](double c, double s) { return real_line(three_points(c, s, 6)); }, boundary, 1.0e2},
        {"real line, normal, order 6", [=](double c, double s) { return real_line(normal_moments(c, s, 6)); }, strict, 1.0e2},
        {"real line, normal, order 8", [=](double c, double s) { return real_line(normal_moments(c, s, 8)); }, strict, 1.0e1},
        {"real line, four points, order 8", [=](double c, double s) { return real_line(four_points(c, s, 8)); }, boundary, 1.0e1},
        {"real line, two points, order 4", [=](double c, double s) { return real_line(two_points(c, s, 4)); }, boundary, 1.0e3},
        // The next column lies in the range of the singular H_2.
        {"real line, two points, order 5", [=](double c, double s) { return real_line(two_points(c, s, 5)); }, boundary, 1.0e3},
        {"real line, variance s^2", [=](double c, double s) { return real_line(MomentSet{1.0, c, c * c + s * s}); }, strict, 1.0e3},
        {"real line, variance -s^2", [=](double c, double s) { return real_line(MomentSet{1.0, c, c * c - s * s}); }, not_realizable, 1.0e3},
        // H_2 has a negative eigenvalue.
        {"real line, two points, m4 lowered", [=](double c, double s) { return real_line(with_offset(two_points(c, s, 4), 4, -0.5, s)); },
         not_realizable, 1.0e3},
        {"real line, three points, m6 lowered", [=](double c, double s) { return real_line(with_offset(three_points(c, s, 6), 6, -0.1, s)); },
         not_realizable, 1.0e1},
        // The next column leaves the range of the singular H_2.
        {"real line, two points, m5 raised", [=](double c, double s) { return real_line(with_offset(two_points(c, s, 5), 5, 0.5, s)); },
         not_realizable, 1.0e2},
        // H_2 is positive semidefinite but not a flat extension of H_1.
        {"real line, one point, m4 raised", [=](double c, double s) { return real_line(with_offset(discrete_moments({1.0}, {c}, 4), 4, 0.5, s)); },
         not_realizable, 1.0e1},
        {"real line, two points, m6 raised", [=](double c, double s) { return real_line(with_offset(two_points(c, s, 6), 6, 0.5, s)); },
         not_realizable, 2.0},
        {"nonnegative, two points, order 3", [=](double c, double s) { return nonnegative(two_points(c, s, 3)); }, strict, 1.0e3, true},
        {"nonnegative, three points, order 5", [=](double c, double s) { return nonnegative(three_points(c, s, 5)); }, strict, 1.0e3, true},
        {"nonnegative, three points, order 6", [=](double c, double s) { return nonnegative(three_points(c, s, 6)); }, boundary, 1.0e2, true},
        {"interval, three points, order 5", [=](double c, double s)
         { return qmomkit::check_realizability(three_points(c, s, 5), BoundedInterval{c - 2.0 * s, c + 2.0 * s}); }, strict, 1.0e2},
        {"interval, three points, order 6", [=](double c, double s)
         { return qmomkit::check_realizability(three_points(c, s, 6), BoundedInterval{c - 2.0 * s, c + 2.0 * s}); }, boundary, 1.0e2},
        {"interval, points at both ends, order 3",
         [=](double c, double s) { return qmomkit::check_realizability(two_points(c, s, 3), BoundedInterval{c - s, c + s}); }, boundary, 1.0e3},
        {"interval, points outside, order 2", [=](double c, double s)
         { return qmomkit::check_realizability(two_points(c, s, 2), BoundedInterval{c - 0.5 * s, c + 0.5 * s}); }, not_realizable, 1.0e3},
    };
}

constexpr double invariance_spreads[] = {1.0e-9, 1.0e-3, 1.0, 1.0e3, 1.0e9};

void verify_shift_and_scale_preserve_classification()
{
    // Realizability does not depend on where the distribution sits or on the
    // unit of its coordinate, and the relative tolerance applies to the
    // standardized moments, so every classification must hold for any spread
    // and for a location of up to the resolved shift in units of the spread.
    for (const InvarianceCase& test_case : invariance_cases())
    {
        for (const double shift : {0.0, 2.0, 10.0, 1.0e2, 1.0e3})
        {
            if (shift > test_case.resolved_shift || (test_case.positive_shift && shift == 0.0))
            {
                continue;
            }
            for (const double spread : invariance_spreads)
            {
                const auto result = test_case.check(shift * spread, spread);
                if (result.status != test_case.expected)
                {
                    std::cerr << test_case.name << " at shift " << shift << " and spread " << spread << '\n';
                }
                QMOMKIT_VERIFY_CHECK(result.status == test_case.expected);
            }
        }
    }
}

void verify_resolution_limit_never_reverses_classification()
{
    // Beyond the resolved shift, the rounding of the raw input moments hides
    // the standardized moments, and the check may answer boundary or
    // numerically indeterminate instead, but never the opposite verdict: a
    // realizable sequence is never rejected, and a non-realizable one or a
    // boundary one is never strictly realizable. At 10^8 spreads the variance
    // is below the rounding of m_2, so the input moments are those of a single
    // point.
    for (const InvarianceCase& test_case : invariance_cases())
    {
        for (const double shift : {0.0, 2.0, 10.0, 30.0, 1.0e2, 3.0e2, 1.0e3, 3.0e3, 1.0e4, 1.0e5, 1.0e6, 1.0e8})
        {
            if (test_case.positive_shift && shift == 0.0)
            {
                continue;
            }
            for (const double spread : invariance_spreads)
            {
                const RealizabilityStatus status = test_case.check(shift * spread, spread).status;
                const bool undecided = status == RealizabilityStatus::boundary || status == RealizabilityStatus::numerically_indeterminate;
                if (status != test_case.expected && !undecided)
                {
                    std::cerr << test_case.name << " at shift " << shift << " and spread " << spread << '\n';
                }
                QMOMKIT_VERIFY_CHECK(status == test_case.expected || undecided);
            }
        }
    }
}

void verify_exact_point_sequences_are_never_rejected()
{
    // Up to four points at multiples of 1/4 below 16 in magnitude, with weights
    // that are multiples of 1/8, have raw moments up to order 8 that are exact
    // in double precision: each term is an integer below 2^51 times a power of
    // two, and so is the sum of up to four terms. Every such sequence is
    // realizable, so no verdict may reject it. From order 2k on, k points are
    // on the boundary and may not be called strictly realizable either. The
    // points cluster within 3/4 of a center anywhere in that range, so many
    // sequences lie tens of spreads from the origin.
    std::uint64_t state = 0x2545F4914F6CDD1DULL;
    const auto next = [&state](const std::uint64_t range)
    {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return (state >> 33U) % range;
    };
    std::size_t boundary_sequences = 0;
    for (std::size_t sample = 0; sample < 50000; ++sample)
    {
        const std::size_t point_count = 1 + next(4);
        std::vector<double> points;
        std::vector<double> weights;
        const double centre = static_cast<double>(next(121)) - 60.0;
        while (points.size() < point_count)
        {
            const double point = (centre + static_cast<double>(next(7)) - 3.0) / 4.0;
            if (std::find(points.begin(), points.end(), point) == points.end())
            {
                points.push_back(point);
                weights.push_back(static_cast<double>(1 + next(8)) / 8.0);
            }
        }
        const std::size_t maximum_order = 2 + next(7);
        std::vector<double> moments(maximum_order + 1, 0.0);
        for (std::size_t point = 0; point < points.size(); ++point)
        {
            double term = weights[point];
            for (std::size_t order = 0; order <= maximum_order; ++order)
            {
                moments[order] += term;
                term *= points[point];
            }
        }

        const RealizabilityStatus status = qmomkit::check_realizability(MomentSet(moments), MomentSupport::real_line).status;
        const bool on_boundary = maximum_order >= 2 * point_count;
        boundary_sequences += on_boundary ? 1 : 0;
        if (status == RealizabilityStatus::not_realizable || (on_boundary && status == RealizabilityStatus::strictly_realizable))
        {
            std::cerr << point_count << " points, order " << maximum_order << ", sample " << sample << '\n';
        }
        QMOMKIT_VERIFY_CHECK(status != RealizabilityStatus::not_realizable);
        QMOMKIT_VERIFY_CHECK(!on_boundary || status != RealizabilityStatus::strictly_realizable);
    }
    QMOMKIT_VERIFY_CHECK(boundary_sequences > 20000);
}

void verify_near_boundary_and_scaling()
{
    // Near-boundary cases make tolerance semantics explicit. The relative
    // tolerance applies to the standardized moment matrices. Two points at -1
    // and 1 put H_2 exactly on the boundary at order 4; moving m_4 by 1e-10
    // keeps its zero eigenvalue within the loose tolerance, and by 1e-6 moves
    // it decisively to either side.
    RealizabilityOptions loose;
    loose.tolerance = 1.0e-8;
    const auto two_points = [](const double fourth_moment) { return MomentSet{1.0, 0.0, 1.0, 0.0, fourth_moment}; };
    QMOMKIT_VERIFY_CHECK(qmomkit::check_realizability(two_points(1.0), MomentSupport::real_line, loose).status == RealizabilityStatus::boundary);
    QMOMKIT_VERIFY_CHECK(qmomkit::check_realizability(two_points(1.0 - 1.0e-10), MomentSupport::real_line, loose).status ==
                         RealizabilityStatus::numerically_indeterminate);
    QMOMKIT_VERIFY_CHECK(qmomkit::check_realizability(two_points(1.0 + 1.0e-10), MomentSupport::real_line, loose).status ==
                         RealizabilityStatus::numerically_indeterminate);
    QMOMKIT_VERIFY_CHECK(qmomkit::check_realizability(two_points(1.0 - 1.0e-6), MomentSupport::real_line, loose).status ==
                         RealizabilityStatus::not_realizable);
    QMOMKIT_VERIFY_CHECK(qmomkit::check_realizability(two_points(1.0 + 1.0e-6), MomentSupport::real_line, loose).status ==
                         RealizabilityStatus::strictly_realizable);

    // A variance of 1e-10 or 1e-6 is decided at its own scale: negative is not
    // realizable, positive is strictly realizable.
    for (const double variance : {1.0e-10, 1.0e-6})
    {
        QMOMKIT_VERIFY_CHECK(qmomkit::check_realizability(MomentSet{1.0, 0.0, -variance}, MomentSupport::real_line, loose).status ==
                             RealizabilityStatus::not_realizable);
        QMOMKIT_VERIFY_CHECK(qmomkit::check_realizability(MomentSet{1.0, 0.0, variance}, MomentSupport::real_line, loose).status ==
                             RealizabilityStatus::strictly_realizable);
    }

    // Positive number-density scaling preserves a strict classification.
    const MomentSet strict = standard_normal_moments(4);
    for (const double scale : {1.0, 1.0e-14, 1.0e-50, 1.0e-100, 1.0e-200, 1.0e50})
    {
        const auto result = qmomkit::check_realizability(scaled_moments(strict, scale), MomentSupport::real_line);
        QMOMKIT_VERIFY_CHECK(result.status == RealizabilityStatus::strictly_realizable);
    }

    // Non-diagonal matrices exercise Jacobi rotations. These boundary and
    // non-realizable cases must retain their classifications at low density.
    const MomentSet boundary{1.0, 1.0, 1.0};
    const MomentSet non_realizable{1.0, 2.0, 1.0};
    for (const double scale : {1.0, 1.0e-50})
    {
        QMOMKIT_VERIFY_CHECK(qmomkit::check_realizability(scaled_moments(boundary, scale), MomentSupport::real_line).status ==
                             RealizabilityStatus::boundary);
        QMOMKIT_VERIFY_CHECK(qmomkit::check_realizability(scaled_moments(non_realizable, scale), MomentSupport::real_line).status ==
                             RealizabilityStatus::not_realizable);
    }

    // The H_1 row sum overflows at this scale although its entries and
    // eigenvalues remain finite. Relative classification must remain unchanged.
    const MomentSet norm_overflow_non_realizable{1.0, 0.9, 0.1};
    const double norm_overflow_scale = 0.6 * std::numeric_limits<double>::max();
    QMOMKIT_VERIFY_CHECK(qmomkit::check_realizability(norm_overflow_non_realizable, MomentSupport::real_line).status ==
                         RealizabilityStatus::not_realizable);
    QMOMKIT_VERIFY_CHECK(
        qmomkit::check_realizability(scaled_moments(norm_overflow_non_realizable, norm_overflow_scale), MomentSupport::real_line).status ==
        RealizabilityStatus::not_realizable);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"diagnostic_success_cases", verify_diagnostic_success_cases},
            {"diagnostic_failure_cases", verify_diagnostic_failure_cases},
            {"boundary_diagnostics", verify_boundary_diagnostics},
            {"hamburger", verify_hamburger},
            {"stieltjes", verify_stieltjes},
            {"hausdorff", verify_hausdorff},
            {"hausdorff_on_other_intervals", verify_hausdorff_on_other_intervals},
            {"even_order_boundaries", verify_even_order_boundaries},
            {"shift_and_scale_preserve_classification", verify_shift_and_scale_preserve_classification},
            {"resolution_limit_never_reverses_classification", verify_resolution_limit_never_reverses_classification},
            {"exact_point_sequences_are_never_rejected", verify_exact_point_sequences_are_never_rejected},
            {"near_boundary_and_scaling", verify_near_boundary_and_scaling},
        });
}
