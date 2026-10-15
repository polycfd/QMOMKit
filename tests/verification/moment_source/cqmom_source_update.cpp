// The CQMOM moment-source update called directly, advance_cqmom_moments_with_source.
//
// What is verified
//   That one call advances the required moments of a CQMOM layout by exactly
//   the SSP-RK3 step of the given rate, entry by entry in layout order, and
//   keeps the layout; that a boundary contribution is added to the source rate
//   and both are evaluated at the SSP-RK3 stage times; and that a realizability
//   policy for one coordinate checks that coordinate's marginal moments, in
//   FinalOnly mode the step result and in EveryStage mode also its stages,
//   halves a rejected step when retries are allowed, fails without them and
//   returns the input moments, also after an accepted part of the interval,
//   and leaves the other coordinates unchecked.
//
// How the reference is constructed
//   By hand, from the weighted points the moments are built from, all with
//   node counts {2, 2}:
//   - relaxation: the points (-0.5, -1.0), (-0.5, 0.8), (1.2, -1.0),
//     (1.2, 0.8) with weights 0.20, 0.30, 0.15, 0.35 relax toward the same
//     points with weights 0.25, 0.20, 0.20, 0.35,
//     dM_a/dt = -r (M_a - M*_a). For this linear rate one SSP-RK3 step of size
//     dt is exactly M* + P(-r dt) (M - M*), with
//     P(z) = 1 + z + z^2/2 + z^3/6, so n steps give P(-r dt)^n. The rate reads
//     each entry's multi-index from the source context, so values and indices
//     must stay aligned for the reference to hold.
//   - source and boundary: a constant source s_a = 0.5 c^a and a boundary
//     birth g_a(t) = (0.8 - 0.6 t + 1.5 t^2) b^a at c = (0.7, -0.2) and
//     b = (0.3, -0.4), where x^a is the monomial of the multi-index a. SSP-RK3
//     weighs its stage times t, t + dt and t + dt/2 like Simpson's rule, which
//     integrates a quadratic exactly, so one step adds
//     s_a dt + b^a (0.8 dt - 0.3 (t1^2 - t0^2) + 0.5 (t1^3 - t0^3)).
//   - coordinate policies: independent coordinates, the first at -0.5 and 1.2
//     with weights 0.4 and 0.6, the second at -1 and 1 with weight 0.5 each,
//     so M_(k,j) = A_k B_j with B = (1, 0, 1, 0). The rate -r M_a on every
//     entry with second exponent 2 multiplies B_2, the second moment of the
//     second coordinate's marginal, by exactly P(-r dt) per step, and its
//     first stage by 1 - r dt. The real-line check rejects a negative B_2:
//     - FinalOnly, r = 2, one retry: the full step of 1.5 gives P(-3) = -2 and
//       is rejected; each half step has the negative first stage 1 - 1.5,
//       which FinalOnly ignores, and the result P(-1.5), so B_2 = P(-1.5)^2;
//     - EveryStage, r = 1, one retry: the full step has the negative first
//       stage 1 - 1.5 and is rejected; each half step of 0.75 is accepted, so
//       B_2 = P(-0.75)^2;
//     - EveryStage without retries fails with rejected_by_predicate;
//     - EveryStage with r = 1 before t = 0.75 and r = 2 from then on, one
//       retry and a minimum time step of 0.5: the full step has the negative
//       first stage 1 - 1.5; the half step to 0.75 has the stages 0.25,
//       0.71875 and 0.453125 and is accepted; the rest of the interval has the
//       negative first stage 1 - 2 * 0.75, and halving it would go below the
//       minimum time step, so the update fails with minimum_timestep_reached
//       and returns its input, not the moments at 0.75;
//     - the FinalOnly policy on the first coordinate instead accepts the full
//       step, B_2 = P(-3) = -2, since the first coordinate's marginal does not
//       change.
//
// Error metric and pass criteria
//   Every required moment within 1.0e-13 of the reference, relative to
//   max(|reference|, 1); the references are exact up to rounding, and the
//   observed errors are of order 1e-16. Statuses and retry counts must match
//   exactly, and a failed update must return its input moments exactly.
//
// Reproduce
//   ctest -R verification_moment_source_cqmom_update_
//
// Build modes
//   Both. Checks go through QMOMKIT_VERIFY_*, not assert().

#include "../verification_helpers.hpp"
#include "qmomkit/qmomkit.hpp"
#include "test_case_dispatch.hpp"

#include <cstddef>
#include <optional>
#include <vector>

namespace
{

constexpr double tolerance = 1.0e-13;

struct WeightedPoint
{
    double first = 0.0;
    double second = 0.0;
    double weight = 0.0;
};

double ssp_rk3_factor(const double z)
{
    return 1.0 + z + (z * z / 2.0) + (z * z * z / 6.0);
}

double integer_power(const double base, const std::size_t exponent)
{
    double value = 1.0;
    for (std::size_t power = 0; power < exponent; ++power)
    {
        value *= base;
    }
    return value;
}

double monomial(const double first, const double second, const qmomkit::MultiIndexView index)
{
    return integer_power(first, index[0]) * integer_power(second, index[1]);
}

double moment_of(const std::vector<WeightedPoint>& points, const qmomkit::MultiIndexView index)
{
    double moment = 0.0;
    for (const WeightedPoint& point : points)
    {
        moment += point.weight * monomial(point.first, point.second, index);
    }
    return moment;
}

const qmomkit::CqmomMomentLayout& layout()
{
    static const qmomkit::CqmomMomentLayout instance = qmomkit::create_cqmom_moment_layout(qmomkit::PhaseSpaceHierarchy({2, 2}));
    return instance;
}

qmomkit::CqmomMomentData moments_of(const std::vector<WeightedPoint>& points)
{
    std::vector<double> values;
    for (const qmomkit::MultiIndexView index : layout().required_indices)
    {
        values.push_back(moment_of(points, index));
    }
    return qmomkit::CqmomMomentData{layout(), values};
}

void verify_moments(const qmomkit::CqmomMomentData& actual, const std::vector<double>& expected)
{
    QMOMKIT_VERIFY_CHECK(actual.layout.required_indices == layout().required_indices);
    QMOMKIT_VERIFY_CHECK(actual.values.size() == expected.size());
    for (std::size_t entry = 0; entry < expected.size(); ++entry)
    {
        QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(actual.values[entry], expected[entry], 0.0, tolerance);
    }
}

void verify_accepted(const qmomkit::CqmomMomentSourceUpdateResult& result, const std::size_t expected_retries)
{
    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.diagnostics.accepted);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.retries == expected_retries);
}

const std::vector<WeightedPoint> initial_points{{-0.5, -1.0, 0.20}, {-0.5, 0.8, 0.30}, {1.2, -1.0, 0.15}, {1.2, 0.8, 0.35}};
const std::vector<WeightedPoint> target_points{{-0.5, -1.0, 0.25}, {-0.5, 0.8, 0.20}, {1.2, -1.0, 0.20}, {1.2, 0.8, 0.35}};

void verify_relaxation_matches_the_ssp_rk3_step()
{
    constexpr double rate = 0.45;
    constexpr double time_step = 0.06;
    const qmomkit::CqmomMomentSource relaxation = [](double, const qmomkit::CqmomSourceContext& context)
    {
        qmomkit::CqmomMomentRate result{std::vector<double>(context.size(), 0.0)};
        for (std::size_t entry = 0; entry < context.size(); ++entry)
        {
            result.values[entry] = -rate * (context.values()[entry] - moment_of(target_points, context.index(entry)));
        }
        return result;
    };

    const qmomkit::CqmomMomentData initial = moments_of(initial_points);
    const qmomkit::CqmomMomentData target = moments_of(target_points);
    const double factor = ssp_rk3_factor(-rate * time_step);

    // One step from t = 0.1.
    const qmomkit::CqmomMomentSourceUpdateResult one_step = qmomkit::advance_cqmom_moments_with_source(initial, 0.1, time_step, relaxation);
    verify_accepted(one_step, 0);
    QMOMKIT_VERIFY_CHECK(one_step.diagnostics.time_step_used == time_step);
    std::vector<double> expected(initial.values.size());
    for (std::size_t entry = 0; entry < expected.size(); ++entry)
    {
        expected[entry] = target.values[entry] + factor * (initial.values[entry] - target.values[entry]);
    }
    verify_moments(one_step.moments, expected);

    // Five consecutive steps.
    qmomkit::CqmomMomentData moments = initial;
    double time = 0.1;
    for (int step = 0; step < 5; ++step)
    {
        qmomkit::CqmomMomentSourceUpdateResult result = qmomkit::advance_cqmom_moments_with_source(moments, time, time_step, relaxation);
        verify_accepted(result, 0);
        moments = std::move(result.moments);
        time += time_step;
    }
    const double five_step_factor = factor * factor * factor * factor * factor;
    for (std::size_t entry = 0; entry < expected.size(); ++entry)
    {
        expected[entry] = target.values[entry] + five_step_factor * (initial.values[entry] - target.values[entry]);
    }
    verify_moments(moments, expected);
}

void verify_source_and_boundary_are_integrated_exactly()
{
    constexpr double source_first = 0.7;
    constexpr double source_second = -0.2;
    constexpr double birth_first = 0.3;
    constexpr double birth_second = -0.4;
    const qmomkit::CqmomMomentSource constant_source = [](double, const qmomkit::CqmomSourceContext& context)
    {
        qmomkit::CqmomMomentRate result{std::vector<double>(context.size(), 0.0)};
        for (std::size_t entry = 0; entry < context.size(); ++entry)
        {
            result.values[entry] = 0.5 * monomial(source_first, source_second, context.index(entry));
        }
        return result;
    };
    const qmomkit::CqmomMomentBoundaryContribution birth = [](const double time, const qmomkit::CqmomSourceContext& context)
    {
        const double flux = 0.8 - 0.6 * time + 1.5 * time * time;
        qmomkit::CqmomMomentRate result{std::vector<double>(context.size(), 0.0)};
        for (std::size_t entry = 0; entry < context.size(); ++entry)
        {
            result.values[entry] = flux * monomial(birth_first, birth_second, context.index(entry));
        }
        return result;
    };

    constexpr double start = 0.2;
    constexpr double time_step = 0.3;
    constexpr double end = start + time_step;
    const qmomkit::CqmomMomentData initial = moments_of(initial_points);
    const qmomkit::CqmomMomentSourceUpdateResult result =
        qmomkit::advance_cqmom_moments_with_source(initial, start, time_step, constant_source, birth, qmomkit::MomentSourceUpdateOptions{});
    verify_accepted(result, 0);

    const double integrated_flux = 0.8 * time_step - 0.3 * (end * end - start * start) + 0.5 * (end * end * end - start * start * start);
    std::vector<double> expected(initial.values.size());
    for (std::size_t entry = 0; entry < expected.size(); ++entry)
    {
        const qmomkit::MultiIndexView index = layout().required_indices[entry];
        expected[entry] = initial.values[entry] + 0.5 * monomial(source_first, source_second, index) * time_step +
                          integrated_flux * monomial(birth_first, birth_second, index);
    }
    verify_moments(result.moments, expected);
}

// Independent coordinates: the first at -0.5 and 1.2, the second at -1 and 1.
const std::vector<WeightedPoint> independent_points{{-0.5, -1.0, 0.2}, {-0.5, 1.0, 0.2}, {1.2, -1.0, 0.3}, {1.2, 1.0, 0.3}};
constexpr double policy_time_step = 1.5;

struct PolicySettings
{
    std::size_t policy_coordinate = 1;
    qmomkit::RealizabilityCheckMode check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;
    double rate = 1.0;
    std::size_t max_retries = 0;
    // From this time on, the rate is `late_rate`.
    double switch_time = 1.0e300;
    double late_rate = 0.0;
    double min_time_step = 0.0;
};

qmomkit::CqmomMomentSourceUpdateResult advance_with_policy(const PolicySettings& settings)
{
    // dM_a/dt = -r(t) M_a on every entry whose second exponent is 2.
    const qmomkit::CqmomMomentSource decay = [settings](const double time, const qmomkit::CqmomSourceContext& context)
    {
        const double rate = time < settings.switch_time ? settings.rate : settings.late_rate;
        qmomkit::CqmomMomentRate result{std::vector<double>(context.size(), 0.0)};
        for (std::size_t entry = 0; entry < context.size(); ++entry)
        {
            if (context.index(entry)[1] == 2)
            {
                result.values[entry] = -rate * context.values()[entry];
            }
        }
        return result;
    };

    qmomkit::RealizabilityPolicy policy;
    policy.check_mode = settings.check_mode;
    policy.support = qmomkit::MomentSupport::real_line;
    qmomkit::MomentSourceUpdateOptions options;
    options.coordinate_realizability_policies.resize(2);
    options.coordinate_realizability_policies[settings.policy_coordinate] = policy;
    options.ssp_rk3_options.max_retries = settings.max_retries;
    options.ssp_rk3_options.min_time_step = settings.min_time_step;
    return qmomkit::advance_cqmom_moments_with_source(moments_of(independent_points), 0.0, policy_time_step, decay, options);
}

// The moments of the independent points with B_2 replaced by `second_moment`.
std::vector<double> independent_moments_with(const double second_moment)
{
    const qmomkit::CqmomMomentData initial = moments_of(independent_points);
    std::vector<double> expected = initial.values;
    for (std::size_t entry = 0; entry < expected.size(); ++entry)
    {
        if (layout().required_indices[entry][1] == 2)
        {
            expected[entry] *= second_moment;
        }
    }
    return expected;
}

void verify_final_only_coordinate_policy_halves_a_rejected_step()
{
    const qmomkit::CqmomMomentSourceUpdateResult result = advance_with_policy({1, qmomkit::RealizabilityCheckMode::FinalOnly, 2.0, 1});
    verify_accepted(result, 1);
    const double half_step_factor = ssp_rk3_factor(-2.0 * 0.5 * policy_time_step);
    verify_moments(result.moments, independent_moments_with(half_step_factor * half_step_factor));
}

void verify_every_stage_coordinate_policy_halves_a_rejected_step()
{
    const qmomkit::CqmomMomentSourceUpdateResult result = advance_with_policy({1, qmomkit::RealizabilityCheckMode::EveryStage, 1.0, 1});
    verify_accepted(result, 1);
    const double half_step_factor = ssp_rk3_factor(-0.5 * policy_time_step);
    verify_moments(result.moments, independent_moments_with(half_step_factor * half_step_factor));
}

void verify_rejected_step_without_retries_returns_the_input()
{
    const qmomkit::CqmomMomentSourceUpdateResult result = advance_with_policy({1, qmomkit::RealizabilityCheckMode::EveryStage, 1.0, 0});
    QMOMKIT_VERIFY_CHECK(!result.ok());
    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::MomentSourceUpdateStatus::rejected_by_predicate);
    QMOMKIT_VERIFY_CHECK(result.moments.values == moments_of(independent_points).values);
}

void verify_failed_update_discards_an_accepted_part_of_the_interval()
{
    const qmomkit::CqmomMomentSourceUpdateResult result =
        advance_with_policy({1, qmomkit::RealizabilityCheckMode::EveryStage, 1.0, 1, 0.75, 2.0, 0.5});
    QMOMKIT_VERIFY_CHECK(!result.ok());
    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::MomentSourceUpdateStatus::minimum_timestep_reached);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.retries == 1);
    QMOMKIT_VERIFY_CHECK(result.moments.values == moments_of(independent_points).values);
}

void verify_policy_checks_only_its_own_coordinate()
{
    const qmomkit::CqmomMomentSourceUpdateResult result = advance_with_policy({0, qmomkit::RealizabilityCheckMode::FinalOnly, 2.0, 1});
    verify_accepted(result, 0);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.time_step_used == policy_time_step);
    verify_moments(result.moments, independent_moments_with(ssp_rk3_factor(-2.0 * policy_time_step)));
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"relaxation_matches_the_ssp_rk3_step", verify_relaxation_matches_the_ssp_rk3_step},
            {"source_and_boundary_are_integrated_exactly", verify_source_and_boundary_are_integrated_exactly},
            {"final_only_coordinate_policy_halves_a_rejected_step", verify_final_only_coordinate_policy_halves_a_rejected_step},
            {"every_stage_coordinate_policy_halves_a_rejected_step", verify_every_stage_coordinate_policy_halves_a_rejected_step},
            {"rejected_step_without_retries_returns_the_input", verify_rejected_step_without_retries_returns_the_input},
            {"failed_update_discards_an_accepted_part_of_the_interval", verify_failed_update_discards_an_accepted_part_of_the_interval},
            {"policy_checks_only_its_own_coordinate", verify_policy_checks_only_its_own_coordinate},
        });
}
