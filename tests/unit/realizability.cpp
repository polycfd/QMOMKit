#include "test_case_dispatch.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <qmomkit/realizability.hpp>
#include <qmomkit/realizability_policy.hpp>
#include <qmomkit/runge_kutta.hpp>
#include <utility>
#include <vector>

namespace
{

using qmomkit::BoundedInterval;
using qmomkit::MomentSet;
using qmomkit::MomentSupport;
using qmomkit::RealizabilityOptions;
using qmomkit::RealizabilityStatus;

void test_invalid_inputs()
{
    assert(qmomkit::check_realizability(MomentSet{}, MomentSupport::real_line).status == RealizabilityStatus::invalid_input);
    assert(qmomkit::check_realizability(MomentSet{1.0, std::numeric_limits<double>::quiet_NaN()}, MomentSupport::real_line).status ==
           RealizabilityStatus::invalid_input);
    assert(qmomkit::check_realizability(MomentSet{1.0, std::numeric_limits<double>::infinity()}, MomentSupport::real_line).status ==
           RealizabilityStatus::invalid_input);

    RealizabilityOptions bad_options;
    bad_options.tolerance = 0.0;
    assert(qmomkit::check_realizability(MomentSet{1.0}, MomentSupport::real_line, bad_options).status == RealizabilityStatus::invalid_input);

    assert(qmomkit::check_realizability(MomentSet{1.0, 0.5}, BoundedInterval{1.0, 1.0}).status == RealizabilityStatus::invalid_input);
    assert(qmomkit::check_realizability(MomentSet{1.0, 0.5}, BoundedInterval{2.0, 1.0}).status == RealizabilityStatus::invalid_input);
    assert(qmomkit::check_realizability(MomentSet{0.0, 0.0}, BoundedInterval{1.0, 1.0}).status == RealizabilityStatus::invalid_input);
    assert(qmomkit::check_realizability(MomentSet{1.0, 0.5}, MomentSupport::bounded_interval).status == RealizabilityStatus::invalid_input);
}

void test_m0_policy()
{
    assert(qmomkit::check_realizability(MomentSet{-1.0, 0.0}, MomentSupport::real_line).status == RealizabilityStatus::not_realizable);

    const auto zero = qmomkit::check_realizability(MomentSet{0.0, 0.0, 0.0}, MomentSupport::nonnegative);
    assert(zero.status == RealizabilityStatus::boundary);
    assert(zero.realizable());
    assert(zero.on_boundary());

    const auto inconsistent_zero = qmomkit::check_realizability(MomentSet{0.0, 1.0}, MomentSupport::real_line);
    assert(inconsistent_zero.status == RealizabilityStatus::not_realizable);
    assert(!inconsistent_zero.realizable());
}

void test_matrix_diagnostics()
{
    const auto hamburger = qmomkit::check_realizability(MomentSet{1.0, 0.0, -1.0}, MomentSupport::real_line);
    assert(hamburger.status == RealizabilityStatus::not_realizable);
    assert(hamburger.checked_order == 2);
    assert(hamburger.failed_matrix == "H_1");
    assert(hamburger.failed_matrix_order == 1);
    assert(hamburger.minimum_pivot < 0.0);

    const auto stieltjes = qmomkit::check_realizability(MomentSet{1.0, -1.0}, MomentSupport::nonnegative);
    assert(stieltjes.status == RealizabilityStatus::not_realizable);
    assert(stieltjes.checked_order == 1);
    assert(stieltjes.failed_matrix == "S_0");

    const auto hausdorff = qmomkit::check_realizability(MomentSet{1.0, 2.0}, BoundedInterval{0.0, 1.0});
    assert(hausdorff.status == RealizabilityStatus::not_realizable);
    assert(hausdorff.checked_order == 1);
    assert(hausdorff.failed_matrix == "L_upper_0");
}

void test_status_helpers()
{
    const auto strict = qmomkit::check_realizability(MomentSet{1.0, 0.0, 1.0}, MomentSupport::real_line);
    assert(strict.status == RealizabilityStatus::strictly_realizable);
    assert(strict.realizable());
    assert(!strict.on_boundary());

    const auto boundary = qmomkit::check_realizability(MomentSet{1.0, 2.0, 4.0}, MomentSupport::real_line);
    assert(boundary.status == RealizabilityStatus::boundary);
    assert(boundary.realizable());
    assert(boundary.on_boundary());

    // Two points at -1 and 1 put H_2 exactly on the boundary at order 4; moving
    // m_4 by 1e-10 moves its zero eigenvalue by less than the loose tolerance,
    // so neither direction can be decided.
    RealizabilityOptions loose;
    loose.tolerance = 1.0e-8;
    const auto indeterminate = qmomkit::check_realizability(MomentSet{1.0, 0.0, 1.0, 0.0, 1.0 - 1.0e-10}, MomentSupport::real_line, loose);
    assert(indeterminate.status == RealizabilityStatus::numerically_indeterminate);
    assert(!indeterminate.realizable());
    assert(!indeterminate.on_boundary());

    const auto small_positive = qmomkit::check_realizability(MomentSet{1.0, 0.0, 1.0, 0.0, 1.0 + 1.0e-10}, MomentSupport::real_line, loose);
    assert(small_positive.status == RealizabilityStatus::numerically_indeterminate);
    assert(!small_positive.on_boundary());

    // The relative tolerance applies to the standardized moments, so a tiny
    // variance is decided at its own scale: negative is not realizable, and
    // positive is a narrow but strictly realizable distribution.
    assert(qmomkit::check_realizability(MomentSet{1.0, 0.0, -1.0e-10}, MomentSupport::real_line, loose).status ==
           RealizabilityStatus::not_realizable);
    assert(qmomkit::check_realizability(MomentSet{1.0, 0.0, 1.0e-10}, MomentSupport::real_line, loose).status ==
           RealizabilityStatus::strictly_realizable);
}

void test_even_order_hamburger_singular_kernel_consistency()
{
    // m2 = integral(x^2 dn) = 0 forces any nonnegative representing
    // measure on the real line to be concentrated at zero, hence m4 = 0.
    // This exact sequence is therefore not Hamburger-realizable despite the
    // PSD Hankel matrix diag(1, 0, 1).
    const MomentSet inconsistent{1.0, 0.0, 0.0, 0.0, 1.0};
    const auto result = qmomkit::check_realizability(inconsistent, MomentSupport::real_line);

    assert(result.status == RealizabilityStatus::not_realizable);
    assert(!result.realizable());
    qmomkit::RealizabilityPolicy policy;
    policy.check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;
    policy.accept_boundary = true;
    assert(
        !qmomkit::accepts_realizability_policy(policy, inconsistent.values(), [](const std::vector<double>& values) { return MomentSet(values); }));
}

void test_every_eigenvalue_is_classified()
{
    // Two points at -1e-5 and 1e-5 give H_2 with eigenvalues of about 1, 1e-10
    // and 0. Under an absolute tolerance of 1e-8 the eigenvalue 1e-10 is
    // neither zero nor positive, so the sequence is indeterminate although the
    // smallest eigenvalue alone would make it boundary.
    const MomentSet moments{1.0, 0.0, 1.0e-10, 0.0, 1.0e-20, 0.0};
    RealizabilityOptions absolute;
    absolute.tolerance = 1.0e-8;
    absolute.relative_tolerance = false;
    const auto result = qmomkit::check_realizability(moments, MomentSupport::real_line, absolute);
    assert(result.status == RealizabilityStatus::numerically_indeterminate);
    assert(result.failed_matrix == "H_2");

    // Standardized, the same points sit at -1 and 1, and every eigenvalue is
    // decided.
    assert(qmomkit::check_realizability(moments, MomentSupport::real_line).status == RealizabilityStatus::boundary);
}

// Raw moments of the normal distribution with mean c and unit variance, from
// m_k = c m_(k-1) + (k - 1) m_(k-2); exact for an integer mean while they
// stay below 2^53.
[[nodiscard]] MomentSet normal_moments(const double mean, const std::size_t maximum_order)
{
    std::vector<double> moments{1.0, mean};
    for (std::size_t order = 2; order <= maximum_order; ++order)
    {
        moments.push_back(mean * moments[order - 1] + static_cast<double>(order - 1) * moments[order - 2]);
    }
    return MomentSet(std::move(moments));
}

void test_flat_extension_within_data_uncertainty_is_indeterminate()
{
    // At a mean of 250 standard deviations the rounding of the raw input
    // moments makes the smallest eigenvalues of H_2 and H_3 zero within the
    // data uncertainty, although both matrices are positive definite as given.
    // The rank mismatch therefore does not decide.
    const auto result = qmomkit::check_realizability(normal_moments(250.0, 6), MomentSupport::real_line);
    assert(result.status == RealizabilityStatus::numerically_indeterminate);
    assert(result.failed_matrix == "H_3_flat_extension");
    assert(result.minimum_pivot > 0.0);
}

void test_lost_information_is_indeterminate()
{
    // At a mean of 10^6 standard deviations the order-6 central moments carry
    // no reliable digit.
    const auto result = qmomkit::check_realizability(normal_moments(1.0e6, 6), MomentSupport::real_line);
    assert(result.status == RealizabilityStatus::numerically_indeterminate);
    assert(result.failed_matrix == "H_3");
}

void test_ssp_rk3_predicate_can_be_explicitly_wrapped()
{
    const qmomkit::OdeRhs rhs = [](double, const std::vector<double>&, std::vector<double>& derivative)
    { std::fill(derivative.begin(), derivative.end(), 0.0); };
    qmomkit::SspRk3Options options;
    options.accept = [](const std::vector<double>& state)
    {
        const auto result = qmomkit::check_realizability(MomentSet(state), MomentSupport::nonnegative);
        return result.status == RealizabilityStatus::strictly_realizable || result.status == RealizabilityStatus::boundary;
    };

    const auto result = qmomkit::ssp_rk3_step(std::vector<double>{1.0, 1.0, 2.0}, 0.0, 0.1, rhs, options);
    assert(result.status == qmomkit::SspRk3Status::success);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"invalid_inputs", test_invalid_inputs},
            {"m0_policy", test_m0_policy},
            {"matrix_diagnostics", test_matrix_diagnostics},
            {"status_helpers", test_status_helpers},
            {"even_order_hamburger_singular_kernel_consistency", test_even_order_hamburger_singular_kernel_consistency},
            {"every_eigenvalue_is_classified", test_every_eigenvalue_is_classified},
            {"flat_extension_within_data_uncertainty_is_indeterminate", test_flat_extension_within_data_uncertainty_is_indeterminate},
            {"lost_information_is_indeterminate", test_lost_information_is_indeterminate},
            {"ssp_rk3_predicate_can_be_explicitly_wrapped", test_ssp_rk3_predicate_can_be_explicitly_wrapped},
        });
}
