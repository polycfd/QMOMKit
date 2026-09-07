#include <cassert>
#include <cmath>
#include <limits>
#include <qmomkit/realizability.hpp>
#include <qmomkit/realizability_policy.hpp>
#include <qmomkit/runge_kutta.hpp>
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

    RealizabilityOptions loose;
    loose.tolerance = 1.0e-8;
    const auto indeterminate = qmomkit::check_realizability(MomentSet{1.0, 0.0, -1.0e-10}, MomentSupport::real_line, loose);
    assert(indeterminate.status == RealizabilityStatus::numerically_indeterminate);
    assert(!indeterminate.realizable());
    assert(!indeterminate.on_boundary());

    const auto small_positive = qmomkit::check_realizability(MomentSet{1.0, 0.0, 1.0e-10}, MomentSupport::real_line, loose);
    assert(small_positive.status == RealizabilityStatus::numerically_indeterminate);
    assert(!small_positive.on_boundary());
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

void test_ssp_rk3_predicate_can_be_explicitly_wrapped()
{
    const auto rhs = [](double, const std::vector<double>&) { return std::vector<double>{0.0, 0.0, 0.0}; };
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

int main()
{
    test_invalid_inputs();
    test_m0_policy();
    test_matrix_diagnostics();
    test_status_helpers();
    test_even_order_hamburger_singular_kernel_consistency();
    test_ssp_rk3_predicate_can_be_explicitly_wrapped();
    return 0;
}
