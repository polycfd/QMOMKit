// Realizability policies and retries in the moment-source step.
//
// What is verified
//   That a QMOM realizability policy in `FinalOnly` mode checks only the result
//   of the SSP-RK3 source step and in `EveryStage` mode also its first stage;
//   that a rejected step is halved, counted as a retry and the rest of the
//   interval completed when retries are allowed; that it fails with the
//   documented status when they are not, or when halving would go below the
//   minimum time step; and that a failed step returns the state it started
//   from, not the state after its first advection half-step.
//
// How the reference is constructed
//   By hand. The moments (1, 0, m2, 0) of two points at -sqrt(m2) and sqrt(m2)
//   relax under dm2/dt = -m2 with m0, m1 and m3 fixed. The accepted cases use
//   zero velocity; the failing ones use dxi/dt = 0.1 xi, so that the first
//   advection half-step changes the moments before the source step fails.
//   One SSP-RK3 step of size dt multiplies m2 by P(-dt), where
//   P(z) = 1 + z + z^2/2 + z^3/6, and its first stage gives m2 (1 - dt). With
//   dt = 1.5 the first stage is negative for any m2 > 0, which the real-line
//   check rejects, while the step result P(-1.5) = 0.0625 is realizable:
//   - FinalOnly checks only the result and accepts the full step, so
//     m2 = P(-1.5), with no retry;
//   - EveryStage rejects the first stage and, with one retry allowed, halves
//     the step to 0.75 and completes the interval with a second step of 0.75,
//     both accepted, so m2 = P(-0.75)^2, with one retry;
//   - EveryStage without retries fails with rejected_by_predicate;
//   - EveryStage with a minimum time step of 1.0 fails with
//     minimum_timestep_reached, since the halved step 0.75 is below it.
//
// Error metric and pass criteria
//   Statuses, failed stage and retry counts must match exactly. The accepted
//   m2 and the nodes +-sqrt(m2) with weights 0.5 must match the reference to
//   1.0e-12, relative to max(|reference|, 1); observed errors are of order
//   1e-16. A failed step must return the moments of its input state exactly
//   and the final time 0.
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

qmomkit::QmomPhaseSpaceEvolutionResult advance(const qmomkit::QmomState& state, const qmomkit::RealizabilityCheckMode check_mode,
                                               const std::size_t max_retries, const double min_time_step, const double velocity_rate)
{
    const qmomkit::PhaseSpaceVelocity velocity = qmomkit::make_phase_space_velocity([velocity_rate](double, const std::vector<double>& xi)
                                                                                    { return std::vector<double>{velocity_rate * xi[0]}; });
    const qmomkit::QmomMomentSource relaxation = [](double, const qmomkit::QmomSourceContext& context)
    {
        qmomkit::QmomMomentRate rate{std::vector<double>(context.size(), 0.0)};
        rate.values[2] = -context.value(2);
        return rate;
    };

    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    options.realizability_policy.check_mode = check_mode;
    options.realizability_policy.support = qmomkit::MomentSupport::real_line;
    options.source_options.ssp_rk3_options.max_retries = max_retries;
    options.source_options.ssp_rk3_options.min_time_step = min_time_step;
    return qmomkit::advance_qmom_phase_space(state, 0.0, time_step, velocity, relaxation, options);
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

void verify_failed_step(const qmomkit::QmomState& input, const qmomkit::QmomPhaseSpaceEvolutionResult& result,
                        const qmomkit::MomentSourceUpdateStatus expected_source_status)
{
    QMOMKIT_VERIFY_CHECK(!result.ok());
    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::PhaseSpaceEvolutionStatus::moment_source_failed);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::moment_source_step);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.moment_source_status == expected_source_status);
    QMOMKIT_VERIFY_CHECK(result.state.moment_values() == input.moment_values());
    QMOMKIT_VERIFY_CHECK(result.final_time == 0.0);
}

void verify_final_only_policy_accepts_the_full_step()
{
    verify_accepted_step(advance(initial_state(), qmomkit::RealizabilityCheckMode::FinalOnly, 0, 0.0, 0.0), ssp_rk3_factor(-time_step), 0);
}

void verify_every_stage_policy_halves_a_rejected_step()
{
    const double half_step_factor = ssp_rk3_factor(-0.5 * time_step);
    verify_accepted_step(advance(initial_state(), qmomkit::RealizabilityCheckMode::EveryStage, 1, 0.0, 0.0), half_step_factor * half_step_factor, 1);
}

void verify_rejected_step_without_retries_fails()
{
    const qmomkit::QmomState state = initial_state();
    verify_failed_step(state, advance(state, qmomkit::RealizabilityCheckMode::EveryStage, 0, 0.0, 0.1),
                       qmomkit::MomentSourceUpdateStatus::rejected_by_predicate);
}

void verify_minimum_time_step_stops_halving()
{
    const qmomkit::QmomState state = initial_state();
    verify_failed_step(state, advance(state, qmomkit::RealizabilityCheckMode::EveryStage, 4, 1.0, 0.1),
                       qmomkit::MomentSourceUpdateStatus::minimum_timestep_reached);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(argc, argv,
                                             {
                                                 {"final_only_policy_accepts_the_full_step", verify_final_only_policy_accepts_the_full_step},
                                                 {"every_stage_policy_halves_a_rejected_step", verify_every_stage_policy_halves_a_rejected_step},
                                                 {"rejected_step_without_retries_fails", verify_rejected_step_without_retries_fails},
                                                 {"minimum_time_step_stops_halving", verify_minimum_time_step_stops_halving},
                                             });
}
