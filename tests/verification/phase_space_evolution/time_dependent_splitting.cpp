// Splitting order when both operators depend on time.
//
// What is verified
//   That `strang` reaches second order and `advsrc` and `srcadv` first order
//   when the velocity and the source depend explicitly on time and do not
//   commute. Each stage must then evaluate its callback on the right time
//   interval: a source stage started at t + dt/2 instead of t, for example,
//   reduces Strang splitting to first order. The other noncommuting cases use
//   callbacks that do not depend on time, so they cannot detect such an error.
//
// How the reference is constructed
//   By hand. Two weighted points (0.4 at 0.5, 0.6 at 2) move under
//   dxi/dt = xi / (1 + t) + (1 + t) cos t, whose solution is
//   xi(t) = (1 + t) (xi(0) + sin t). The first term is the node velocity; the
//   second is a drift written as the moment source dm_k/dt = k c(t) m_{k-1}
//   with c(t) = (1 + t) cos t. The two operators do not commute, the weights
//   are constant, and the moments m_0 to m_3 at t = 1 follow from the exact
//   positions. Node advection uses RK45 at tolerances 1e-13 relative and 1e-15
//   absolute, so the splitting error dominates.
//
// Error metric and pass criteria
//   The error of a run is the largest error of m_0 to m_3 at t = 1, relative
//   to max(|reference|, 1). Each scheme runs with 10, 20, 40 and 80 steps; the
//   error must fall at every refinement, and the observed order
//   log2(e_dt / e_dt/2) on the two finest pairs must be at least 1.75 and 1.85
//   for `strang` and at least 0.75 and 0.85 for `advsrc` and `srcadv`, the
//   thresholds of the other noncommuting QMOM case. The observed orders are
//   2.07 and 2.03 for `strang` and within 0.03 of 1 for the others.
//
// Reproduce
//   ctest -R verification_phase_space_evolution_time_dependent_splitting_
//
// Build modes
//   Both. Checks go through QMOMKIT_VERIFY_*, not assert().

#include "../verification_helpers.hpp"
#include "qmomkit/qmomkit.hpp"
#include "test_case_dispatch.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace
{

constexpr double final_time = 1.0;
constexpr std::array<std::size_t, 4> step_counts{10, 20, 40, 80};
constexpr std::array<double, 2> weights{0.4, 0.6};
constexpr std::array<double, 2> initial_nodes{0.5, 2.0};

std::vector<double> exact_moments(const double time)
{
    std::vector<double> moments(4, 0.0);
    for (std::size_t order = 0; order < moments.size(); ++order)
    {
        for (std::size_t node = 0; node < weights.size(); ++node)
        {
            moments[order] += weights[node] * std::pow((1.0 + time) * (initial_nodes[node] + std::sin(time)), static_cast<int>(order));
        }
    }
    return moments;
}

double run_error(const qmomkit::PhaseSpaceEvolutionScheme scheme, const std::size_t steps)
{
    qmomkit::QmomCase qmom_case;
    qmom_case.moments = qmomkit::MomentSet(exact_moments(0.0));
    qmom_case.node_count = 2;
    qmom_case.qmom_options.adaptive = false;
    qmomkit::QmomState state = qmom_case.initial_state();

    const qmomkit::PhaseSpaceVelocity velocity = qmomkit::make_phase_space_velocity([](const double time, const std::vector<double>& xi)
                                                                                    { return std::vector<double>{xi[0] / (1.0 + time)}; });
    const qmomkit::QmomMomentSource drift = [](const double time, const qmomkit::QmomSourceContext& context)
    {
        qmomkit::QmomMomentRate rate{std::vector<double>(context.size(), 0.0)};
        for (std::size_t order = 1; order < context.size(); ++order)
        {
            rate.values[order] = static_cast<double>(order) * (1.0 + time) * std::cos(time) * context.value(order - 1);
        }
        return rate;
    };

    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    options.scheme = scheme;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::rk45;
    options.node_options.adaptive_options.relative_tolerance = 1.0e-13;
    options.node_options.adaptive_options.absolute_tolerance = 1.0e-15;
    options.node_options.adaptive_options.initial_step = 1.0e-4;

    const double time_step = final_time / static_cast<double>(steps);
    for (std::size_t step = 0; step < steps; ++step)
    {
        const qmomkit::QmomPhaseSpaceEvolutionResult result =
            qmomkit::advance_qmom_phase_space(state, static_cast<double>(step) * time_step, time_step, velocity, drift, options);
        QMOMKIT_VERIFY_CHECK(result.ok());
        state = result.state;
    }

    const std::vector<double> reference = exact_moments(final_time);
    double error = 0.0;
    for (std::size_t order = 0; order < reference.size(); ++order)
    {
        error = std::max(error, std::abs(state.moment_values()[order] - reference[order]) / std::max(1.0, std::abs(reference[order])));
    }
    return error;
}

void verify_scheme_order(const qmomkit::PhaseSpaceEvolutionScheme scheme, const double second_finest_minimum, const double finest_minimum)
{
    std::array<double, step_counts.size()> errors{};
    for (std::size_t run = 0; run < step_counts.size(); ++run)
    {
        errors[run] = run_error(scheme, step_counts[run]);
        QMOMKIT_VERIFY_CHECK(std::isfinite(errors[run]));
        if (run > 0)
        {
            QMOMKIT_VERIFY_CHECK(errors[run] < errors[run - 1]);
        }
    }
    const std::size_t finest = step_counts.size() - 1;
    QMOMKIT_VERIFY_CHECK(std::log2(errors[finest - 2] / errors[finest - 1]) >= second_finest_minimum);
    QMOMKIT_VERIFY_CHECK(std::log2(errors[finest - 1] / errors[finest]) >= finest_minimum);
}

void verify_strang_is_second_order()
{
    verify_scheme_order(qmomkit::PhaseSpaceEvolutionScheme::strang, 1.75, 1.85);
}

void verify_lie_schemes_are_first_order()
{
    verify_scheme_order(qmomkit::PhaseSpaceEvolutionScheme::advsrc, 0.75, 0.85);
    verify_scheme_order(qmomkit::PhaseSpaceEvolutionScheme::srcadv, 0.75, 0.85);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(argc, argv,
                                             {
                                                 {"strang_is_second_order", verify_strang_is_second_order},
                                                 {"lie_schemes_are_first_order", verify_lie_schemes_are_first_order},
                                             });
}
