#include "../verification_helpers.hpp"

#include <cmath>
#include <iostream>
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

[[nodiscard]] MomentSet uniform_01_moments(const std::size_t maximum_order)
{
    std::vector<double> moments(maximum_order + 1, 0.0);
    for (std::size_t order = 0; order <= maximum_order; ++order)
    {
        moments[order] = 1.0 / static_cast<double>(order + 1);
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
    // Negative mass is support-independent: no nonnegative measure can have
    // m0 < 0, so the public failure should name m0 directly.
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

    // Tolerance-controlled near-boundary cases use m0=1, m1=0, and a tiny
    // positive/negative variance. With a loose absolute tolerance, both tiny
    // signs are classified as numerically indeterminate; larger perturbations
    // move deterministically to strict or non-realizable classifications.
    RealizabilityOptions loose;
    loose.tolerance = 1.0e-8;
    loose.relative_tolerance = false;

    const auto tiny_negative = qmomkit::check_realizability(MomentSet{1.0, 0.0, -1.0e-10}, MomentSupport::real_line, loose);
    QMOMKIT_VERIFY_CHECK(tiny_negative.status == RealizabilityStatus::numerically_indeterminate);
    QMOMKIT_VERIFY_CHECK(tiny_negative.checked_order == 2);
    QMOMKIT_VERIFY_CHECK(tiny_negative.failed_matrix == "H_1");
    QMOMKIT_VERIFY_CHECK(tiny_negative.failed_matrix_order == 1);
    QMOMKIT_VERIFY_NEAR(tiny_negative.tolerance_used, loose.tolerance, 0.0);
    QMOMKIT_VERIFY_CHECK(contains_text(tiny_negative.message, "indeterminate"));

    const auto clear_negative = qmomkit::check_realizability(MomentSet{1.0, 0.0, -1.0e-6}, MomentSupport::real_line, loose);
    verify_failure_diagnostics(clear_negative, 2, "H_1", 1, "negative eigenvalue");
    QMOMKIT_VERIFY_NEAR(clear_negative.tolerance_used, loose.tolerance, 0.0);

    const auto clear_positive = qmomkit::check_realizability(MomentSet{1.0, 0.0, 1.0e-6}, MomentSupport::real_line, loose);
    QMOMKIT_VERIFY_CHECK(clear_positive.status == RealizabilityStatus::strictly_realizable);
    QMOMKIT_VERIFY_CHECK(clear_positive.checked_order == 2);
    QMOMKIT_VERIFY_CHECK(clear_positive.failed_matrix == "H_1");
    QMOMKIT_VERIFY_CHECK(clear_positive.minimum_pivot > clear_positive.tolerance_used);
    QMOMKIT_VERIFY_NEAR(clear_positive.tolerance_used, loose.tolerance, 0.0);
}

void verify_hamburger()
{
    // Hamburger checks cover moments on the full real line: strictly
    // realizable advection moments, boundary atomic measures, and basic
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

void verify_near_boundary_and_scaling()
{
    // Near-boundary cases make tolerance semantics explicit, and the scaling
    // loop verifies that relative tolerance handling does not change the
    // realizability classification of a strictly realizable sequence.
    RealizabilityOptions loose;
    loose.tolerance = 1.0e-8;
    QMOMKIT_VERIFY_CHECK(qmomkit::check_realizability(MomentSet{1.0, 0.0, -1.0e-10}, MomentSupport::real_line, loose).status ==
                         RealizabilityStatus::numerically_indeterminate);
    QMOMKIT_VERIFY_CHECK(qmomkit::check_realizability(MomentSet{1.0, 0.0, 1.0e-10}, MomentSupport::real_line, loose).status ==
                         RealizabilityStatus::numerically_indeterminate);
    QMOMKIT_VERIFY_CHECK(qmomkit::check_realizability(MomentSet{1.0, 0.0, -1.0e-6}, MomentSupport::real_line, loose).status ==
                         RealizabilityStatus::not_realizable);
    QMOMKIT_VERIFY_CHECK(qmomkit::check_realizability(MomentSet{1.0, 0.0, 1.0e-6}, MomentSupport::real_line, loose).status ==
                         RealizabilityStatus::strictly_realizable);

    const auto base = standard_normal_moments(4).values();
    for (const double scale : {1.0e-6, 1.0, 1.0e6})
    {
        std::vector<double> scaled = base;
        for (double& value : scaled)
        {
            value *= scale;
        }
        const auto result = qmomkit::check_realizability(MomentSet(std::move(scaled)), MomentSupport::real_line);
        QMOMKIT_VERIFY_CHECK(result.status == RealizabilityStatus::strictly_realizable);
    }
}

} // namespace

int main()
{
    verify_hamburger();
    verify_stieltjes();
    verify_hausdorff();
    verify_near_boundary_and_scaling();
    verify_diagnostic_success_cases();
    verify_diagnostic_failure_cases();
    verify_boundary_diagnostics();
    std::cout << "QMOMKit one-dimensional realizability analytical verification passed\n";
    return 0;
}
