// Realizability policies and retries in the moment-source step.
//
// What is verified
//   That a QMOM realizability policy in `FinalOnly` mode checks the result of
//   the SSP-RK3 source step but not its stages, and in `EveryStage` mode also
//   its first stage; that a rejected step is halved, counted as a retry and
//   the rest of the interval completed when retries are allowed; that it fails
//   with the documented status when they are not, or when halving would go
//   below the minimum time step; and that a failed step returns the state and
//   time it started from, not the state after its first advection half-step.
//
// How the reference is constructed
//   By hand. The moments (1, 0, m2, 0) of two points at -sqrt(m2) and sqrt(m2)
//   relax under dm2/dt = -r m2 with m0, m1 and m3 fixed, over one Strang step
//   of 1.5. The accepted cases use zero velocity and start at t = 0; the
//   failing ones use dxi/dt = 0.1 xi and start at t = 1, so that the first
//   advection half-step changes the moments and a defaulted final time differs
//   from the start. One SSP-RK3 step of size dt multiplies m2 by P(-r dt),
//   where P(z) = 1 + z + z^2/2 + z^3/6, and its first stage gives
//   m2 (1 - r dt). The real-line check rejects a negative m2:
//   - FinalOnly with r = 2 and one retry allowed: the full step gives
//     P(-3) = -2, which FinalOnly rejects; each half step of 0.75 has the
//     negative first stage 1 - 1.5, which FinalOnly ignores, and the result
//     P(-1.5) = 0.0625, so m2 = P(-1.5)^2, with one retry;
//   - EveryStage with r = 1: the full step has the negative first stage
//     1 - 1.5 but the realizable result P(-1.5); EveryStage rejects the first
//     stage and, with one retry allowed, halves the step to 0.75 and completes
//     the interval with a second step of 0.75, both accepted, so
//     m2 = P(-0.75)^2, with one retry;
//   - EveryStage without retries fails with rejected_by_predicate, since the
//     first stage is negative for any m2 > 0;
//   - EveryStage with a minimum time step of 1.0 fails with
//     minimum_timestep_reached, since the halved step 0.75 is below it.
//
// Error metric and pass criteria
//   Statuses, failed stage and retry counts must match exactly. The accepted
//   m2 and the nodes +-sqrt(m2) with weights 0.5 must match the reference to
//   1.0e-12, relative to max(|reference|, 1); observed errors are of order
//   1e-16. A failed step must return the moments of its input state exactly
//   and its start time.
//
// Reproduce
//   ctest -R verification_phase_space_evolution_source_retries_
//
// Build modes
//   Both. Checks go through QMOMKIT_VERIFY_*, not assert().

#include "../verification_helpers.hpp"
#include "qmomkit/qmomkit.hpp"
#include "test_case_dispatch.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

namespace
{

constexpr double tolerance = 1.0e-12;
constexpr double time_step = 1.5;

double ssp_rk3_factor(const double z)
{
    return 1.0 + z + (z * z / 2.0) + (z * z * z / 6.0);
}

qmomkit::QmomState initial_state()
{
    qmomkit::QmomCase qmom_case;
    qmom_case.moments = qmomkit::MomentSet({1.0, 0.0, 1.0, 0.0});
    qmom_case.node_count = 2;
    qmom_case.qmom_options.adaptive = false;
    return qmom_case.initial_state();
}

struct StepSettings
{
    qmomkit::RealizabilityCheckMode check_mode = qmomkit::RealizabilityCheckMode::EveryStage;
    double decay_rate = 1.0; // dm2/dt = -decay_rate m2
    std::size_t max_retries = 0;
    double min_time_step = 0.0;
    double velocity_rate = 0.0; // dxi/dt = velocity_rate xi
    double start_time = 0.0;
};

qmomkit::QmomPhaseSpaceEvolutionResult advance(const qmomkit::QmomState& state, const StepSettings& settings)
{
    const double velocity_rate = settings.velocity_rate;
    const qmomkit::PhaseSpaceVelocity velocity = qmomkit::make_phase_space_velocity([velocity_rate](double, const std::vector<double>& xi)
                                                                                    { return std::vector<double>{velocity_rate * xi[0]}; });
    const double decay_rate = settings.decay_rate;
    const qmomkit::QmomMomentSource relaxation = [decay_rate](double, const qmomkit::QmomSourceContext& context)
    {
        qmomkit::QmomMomentRate rate{std::vector<double>(context.size(), 0.0)};
        rate.values[2] = -decay_rate * context.value(2);
        return rate;
    };

    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    options.realizability_policy.check_mode = settings.check_mode;
    options.realizability_policy.support = qmomkit::MomentSupport::real_line;
    options.source_options.ssp_rk3_options.max_retries = settings.max_retries;
    options.source_options.ssp_rk3_options.min_time_step = settings.min_time_step;
    return qmomkit::advance_qmom_phase_space(state, settings.start_time, time_step, velocity, relaxation, options);
}

void verify_accepted_step(const qmomkit::QmomPhaseSpaceEvolutionResult& result, const double expected_m2, const std::size_t expected_retries)
{
    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.diagnostics.moment_source_retries == expected_retries);
    QMOMKIT_VERIFY_NEAR(result.final_time, time_step, 1.0e-15);

    const std::vector<double>& moments = result.state.moment_values();
    QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(moments[0], 1.0, 0.0, tolerance);
    QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(moments[1], 0.0, 0.0, tolerance);
    QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(moments[2], expected_m2, 0.0, tolerance);
    QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(moments[3], 0.0, 0.0, tolerance);

    const qmomkit::QmomResult& quadrature = result.state.quadrature();
    QMOMKIT_VERIFY_CHECK(quadrature.abscissas.size() == 2);
    const double spread = std::sqrt(expected_m2);
    QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(quadrature.abscissas[0], -spread, 0.0, tolerance);
    QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(quadrature.abscissas[1], spread, 0.0, tolerance);
    QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(quadrature.weights[0], 0.5, 0.0, tolerance);
    QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(quadrature.weights[1], 0.5, 0.0, tolerance);
}

void verify_failed_step(const qmomkit::QmomState& input, const StepSettings& settings, const qmomkit::MomentSourceUpdateStatus expected_source_status)
{
    const qmomkit::QmomPhaseSpaceEvolutionResult result = advance(input, settings);
    QMOMKIT_VERIFY_CHECK(!result.ok());
    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::PhaseSpaceEvolutionStatus::moment_source_failed);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::moment_source_step);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.moment_source_status == expected_source_status);
    QMOMKIT_VERIFY_CHECK(result.state.moment_values() == input.moment_values());
    QMOMKIT_VERIFY_CHECK(result.final_time == settings.start_time);
}

void verify_final_only_policy_checks_only_the_step_result()
{
    // With dm2/dt = -2 m2 the full step gives P(-3) = -2 and is rejected; each
    // half step has the negative first stage 1 - 1.5, which FinalOnly ignores,
    // and the result P(-1.5).
    const double half_step_factor = ssp_rk3_factor(-time_step);
    verify_accepted_step(advance(initial_state(), {qmomkit::RealizabilityCheckMode::FinalOnly, 2.0, 1}), half_step_factor * half_step_factor, 1);
}

void verify_every_stage_policy_halves_a_rejected_step()
{
    const double half_step_factor = ssp_rk3_factor(-0.5 * time_step);
    verify_accepted_step(advance(initial_state(), {qmomkit::RealizabilityCheckMode::EveryStage, 1.0, 1}), half_step_factor * half_step_factor, 1);
}

void verify_rejected_step_without_retries_fails()
{
    verify_failed_step(initial_state(), {qmomkit::RealizabilityCheckMode::EveryStage, 1.0, 0, 0.0, 0.1, 1.0},
                       qmomkit::MomentSourceUpdateStatus::rejected_by_predicate);
}

void verify_minimum_time_step_stops_halving()
{
    verify_failed_step(initial_state(), {qmomkit::RealizabilityCheckMode::EveryStage, 1.0, 4, 1.0, 0.1, 1.0},
                       qmomkit::MomentSourceUpdateStatus::minimum_timestep_reached);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"final_only_policy_checks_only_the_step_result", verify_final_only_policy_checks_only_the_step_result},
            {"every_stage_policy_halves_a_rejected_step", verify_every_stage_policy_halves_a_rejected_step},
            {"rejected_step_without_retries_fails", verify_rejected_step_without_retries_fails},
            {"minimum_time_step_stops_halving", verify_minimum_time_step_stops_halving},
        });
}
