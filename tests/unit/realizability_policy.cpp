#include "test_case_dispatch.hpp"

#include <cassert>
#include <qmomkit/realizability_policy.hpp>
#include <stdexcept>
#include <vector>

namespace
{

qmomkit::MomentSet identity_moments(const std::vector<double>& state)
{
    return qmomkit::MomentSet(state);
}

void test_default_policy_is_disabled()
{
    const qmomkit::RealizabilityPolicy policy;
    assert(policy.check_mode == qmomkit::RealizabilityCheckMode::None);

    const qmomkit::StateAcceptancePredicate predicate = qmomkit::make_realizability_acceptance_predicate(policy, identity_moments);
    assert(!predicate);
    const qmomkit::StageAcceptancePredicate stage_predicate = qmomkit::make_realizability_stage_acceptance_predicate(policy, identity_moments);
    assert(!stage_predicate);
}

void test_disabled_policy_does_not_invoke_converter()
{
    const qmomkit::RealizabilityPolicy policy;
    bool converter_called = false;
    const bool accepted = qmomkit::accepts_realizability_policy(policy, {1.0, 0.0, -1.0},
                                                                [&](const std::vector<double>& state)
                                                                {
                                                                    converter_called = true;
                                                                    return qmomkit::MomentSet(state);
                                                                });

    assert(accepted);
    assert(!converter_called);
}

void test_final_only_policy_invokes_converter_once()
{
    qmomkit::RealizabilityPolicy policy;
    policy.check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;

    int converter_calls = 0;
    const qmomkit::StateAcceptancePredicate predicate = qmomkit::make_realizability_acceptance_predicate(policy,
                                                                                                         [&](const std::vector<double>& state)
                                                                                                         {
                                                                                                             ++converter_calls;
                                                                                                             return qmomkit::MomentSet(state);
                                                                                                         });

    assert(predicate);
    assert(predicate({1.0, 0.0, 1.0}));
    assert(converter_calls == 1);
}

void test_every_stage_policy_uses_stage_predicate()
{
    qmomkit::RealizabilityPolicy policy;
    policy.check_mode = qmomkit::RealizabilityCheckMode::EveryStage;

    int converter_calls = 0;
    const qmomkit::StageAcceptancePredicate predicate = qmomkit::make_realizability_stage_acceptance_predicate(policy,
                                                                                                               [&](const std::vector<double>& state)
                                                                                                               {
                                                                                                                   ++converter_calls;
                                                                                                                   return qmomkit::MomentSet(state);
                                                                                                               });

    assert(predicate);
    assert(predicate(qmomkit::SspRk3Stage::Stage1, {1.0, 0.0, 1.0}));
    assert(predicate(qmomkit::SspRk3Stage::Stage2, {1.0, 0.0, 1.0}));
    assert(predicate(qmomkit::SspRk3Stage::Final, {1.0, 0.0, 1.0}));
    assert(converter_calls == 3);
}

void test_support_and_interval_are_forwarded()
{
    qmomkit::RealizabilityPolicy real_line_policy;
    real_line_policy.check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;
    real_line_policy.support = qmomkit::MomentSupport::real_line;

    assert(qmomkit::accepts_realizability_policy(real_line_policy, {1.0, 2.0}, identity_moments));

    qmomkit::RealizabilityPolicy interval_policy;
    interval_policy.check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;
    interval_policy.support = qmomkit::MomentSupport::bounded_interval;
    interval_policy.interval = qmomkit::BoundedInterval{0.0, 1.0};

    assert(!qmomkit::accepts_realizability_policy(interval_policy, {1.0, 2.0}, identity_moments));
}

void test_explicit_bounded_interval_policy_accepts_realizable_moments()
{
    qmomkit::RealizabilityPolicy policy;
    policy.check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;
    policy.support = qmomkit::MomentSupport::bounded_interval;
    policy.interval = qmomkit::BoundedInterval{0.0, 1.0};

    assert(qmomkit::accepts_realizability_policy(policy, {1.0, 0.5, 1.0 / 3.0}, identity_moments));
}

void test_strict_boundary_and_failure_statuses()
{
    qmomkit::RealizabilityPolicy policy;
    policy.check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;
    policy.support = qmomkit::MomentSupport::real_line;

    assert(qmomkit::accepts_realizability_policy(policy, {1.0, 0.0, 1.0}, identity_moments));
    assert(qmomkit::accepts_realizability_policy(policy, {1.0, 2.0, 4.0}, identity_moments));

    policy.accept_boundary = false;
    assert(qmomkit::accepts_realizability_policy(policy, {1.0e-14, 0.0, 1.0e-14}, identity_moments));
    assert(!qmomkit::accepts_realizability_policy(policy, {1.0, 2.0, 4.0}, identity_moments));
    assert(!qmomkit::accepts_realizability_policy(policy, {1.0, 0.0, -1.0}, identity_moments));
}

void test_numerically_indeterminate_is_rejected()
{
    qmomkit::RealizabilityPolicy policy;
    policy.check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;
    policy.support = qmomkit::MomentSupport::real_line;
    policy.options.tolerance = 1.0e-8;

    assert(!qmomkit::accepts_realizability_policy(policy, {1.0, 0.0, 1.0e-10}, identity_moments));
}

void test_converter_failures_propagate()
{
    qmomkit::RealizabilityPolicy policy;
    policy.check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;

    bool threw = false;
    try
    {
        (void)qmomkit::accepts_realizability_policy(policy, {1.0}, [](const std::vector<double>&) -> qmomkit::MomentSet
                                                    { throw std::runtime_error("converter failed"); });
    }
    catch (const std::runtime_error&)
    {
        threw = true;
    }
    assert(threw);
}

void test_ssp_rk3_accepts_realizable_candidate()
{
    qmomkit::RealizabilityPolicy policy;
    policy.check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;
    policy.support = qmomkit::MomentSupport::nonnegative;

    qmomkit::SspRk3Options options;
    options.accept = qmomkit::make_realizability_acceptance_predicate(policy, identity_moments);

    const auto rhs = [](double, const std::vector<double>&) { return std::vector<double>{0.0, 0.0, 0.0}; };
    const qmomkit::SspRk3Result result = qmomkit::ssp_rk3_step({1.0, 1.0, 2.0}, 0.0, 0.1, rhs, options);

    assert(result.ok());
    assert(result.retries == 0);
}

void test_ssp_rk3_rejects_and_retries_with_realizability_predicate()
{
    qmomkit::RealizabilityPolicy policy;
    policy.check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;
    policy.support = qmomkit::MomentSupport::real_line;

    int converter_calls = 0;
    qmomkit::SspRk3Options options;
    options.accept = qmomkit::make_realizability_acceptance_predicate(policy,
                                                                      [&](const std::vector<double>& state)
                                                                      {
                                                                          ++converter_calls;
                                                                          return qmomkit::MomentSet(state);
                                                                      });
    options.max_retries = 2;

    const auto rhs = [](double, const std::vector<double>&) { return std::vector<double>{0.0, 0.0, -2.0}; };
    const qmomkit::SspRk3Result result = qmomkit::ssp_rk3_step({1.0, 0.0, 1.0}, 0.0, 1.0, rhs, options);

    assert(result.ok());
    assert(result.retries == 1);
    assert(converter_calls == 2);
    assert(result.state[2] == 0.0);
}

void test_ssp_rk3_disabled_policy_installs_no_predicate()
{
    qmomkit::RealizabilityPolicy policy;
    int converter_calls = 0;

    qmomkit::SspRk3Options options;
    options.accept = qmomkit::make_realizability_acceptance_predicate(policy,
                                                                      [&](const std::vector<double>& state)
                                                                      {
                                                                          ++converter_calls;
                                                                          return qmomkit::MomentSet(state);
                                                                      });

    const auto rhs = [](double, const std::vector<double>&) { return std::vector<double>{0.0, 0.0, -2.0}; };
    const qmomkit::SspRk3Result result = qmomkit::ssp_rk3_step({1.0, 0.0, 1.0}, 0.0, 1.0, rhs, options);

    assert(result.ok());
    assert(result.retries == 0);
    assert(converter_calls == 0);
    assert(result.state[2] == -1.0);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"default_policy_is_disabled", test_default_policy_is_disabled},
            {"disabled_policy_does_not_invoke_converter", test_disabled_policy_does_not_invoke_converter},
            {"final_only_policy_invokes_converter_once", test_final_only_policy_invokes_converter_once},
            {"every_stage_policy_uses_stage_predicate", test_every_stage_policy_uses_stage_predicate},
            {"support_and_interval_are_forwarded", test_support_and_interval_are_forwarded},
            {"explicit_bounded_interval_policy_accepts_realizable_moments", test_explicit_bounded_interval_policy_accepts_realizable_moments},
            {"strict_boundary_and_failure_statuses", test_strict_boundary_and_failure_statuses},
            {"numerically_indeterminate_is_rejected", test_numerically_indeterminate_is_rejected},
            {"converter_failures_propagate", test_converter_failures_propagate},
            {"ssp_rk3_accepts_realizable_candidate", test_ssp_rk3_accepts_realizable_candidate},
            {"ssp_rk3_rejects_and_retries_with_realizability_predicate", test_ssp_rk3_rejects_and_retries_with_realizability_predicate},
            {"ssp_rk3_disabled_policy_installs_no_predicate", test_ssp_rk3_disabled_policy_installs_no_predicate},
        });
}
