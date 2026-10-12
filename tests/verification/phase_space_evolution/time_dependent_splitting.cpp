// Splitting with a velocity and a source that both depend on time.
//
// What is verified
//   That `strang`, `advsrc` and `srcadv` reproduce their exact split solutions
//   when the velocity and the source depend explicitly on time and do not
//   commute, which requires every stage to evaluate its callback on the right
//   time interval, and that `strang` then reaches second order. A stage on the
//   wrong interval moves a scheme away from its split solution; a Lie scheme
//   stays first order, so only this comparison detects it there, while a
//   source stage started at t + dt/2 instead of t also reduces Strang splitting
//   to first order. The other noncommuting cases use callbacks that do not
//   depend on time, so they cannot detect such an error.
//
// How the reference is constructed
//   By hand. Two weighted points (0.4 at 0.5, 0.6 at 2) move under
//   dxi/dt = xi / (1 + t) + (1 + t) cos t, whose solution is
//   xi(t) = (1 + t) (xi(0) + sin t). The first term is the node velocity; the
//   second is a drift written as the moment source dm_k/dt = k c(t) m_{k-1}
//   with c(t) = (1 + t) cos t. The two operators do not commute, and the
//   weights are constant. Each operator alone also has a closed-form flow over
//   [s, s + h]: advection maps xi to xi (1 + s + h) / (1 + s), and the drift
//   adds F(s + h) - F(s) with F(t) = (1 + t) sin t + cos t. Composing these
//   flows in the order and on the intervals of a scheme gives its exact split
//   solution. The moments m_0 to m_3 at t = 1 follow from the positions.
//
// Error metric and pass criteria
//   Errors are the largest error of m_0 to m_3 at t = 1, relative to
//   max(|reference|, 1). Each scheme runs with 10, 20, 40 and 80 steps, with
//   RK45 node advection at tolerances 1e-13 relative and 1e-15 absolute.
//   - Against the exact split solution, the deviation is the error of the
//     third-order SSP-RK3 source step, so it must fall at every refinement,
//     converge at an observed order log2(e_dt / e_dt/2) of at least 2.8 on the
//     two finest pairs and stay at most 1e-6 at 80 steps. The observed orders
//     are 2.97 to 3.01, and the deviation at 80 steps is 1e-7; a stage on the
//     wrong interval adds a deviation of order dt.
//   - For `strang`, against the exact solution, the error must fall at every
//     refinement, and the observed order must be at least 1.75 and 1.85 on the
//     two finest pairs, the thresholds of the other noncommuting QMOM case.
//     The observed orders are 2.07 and 2.03. This also checks the Strang split
//     model against the true solution; the Lie split solutions are first order
//     by construction.
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
constexpr double split_deviation_bound = 1.0e-6;

std::vector<double> moments_of(const std::array<double, 2>& nodes)
{
    std::vector<double> moments(4, 0.0);
    for (std::size_t order = 0; order < moments.size(); ++order)
    {
        for (std::size_t node = 0; node < weights.size(); ++node)
        {
            moments[order] += weights[node] * std::pow(nodes[node], static_cast<int>(order));
        }
    }
    return moments;
}

std::vector<double> exact_moments(const double time)
{
    return moments_of({(1.0 + time) * (initial_nodes[0] + std::sin(time)), (1.0 + time) * (initial_nodes[1] + std::sin(time))});
}

// Exact flows of the two operators over [start, start + step].
double advection_flow(const double xi, const double start, const double step)
{
    return xi * (1.0 + start + step) / (1.0 + start);
}

double drift_flow(const double xi, const double start, const double step)
{
    const auto antiderivative = [](const double time) { return ((1.0 + time) * std::sin(time)) + std::cos(time); };
    return xi + antiderivative(start + step) - antiderivative(start);
}

// One step of the scheme with the exact flows, on the intervals the scheme
// assigns to its stages.
double split_step(const qmomkit::PhaseSpaceEvolutionScheme scheme, const double xi, const double start, const double step)
{
    switch (scheme)
    {
    case qmomkit::PhaseSpaceEvolutionScheme::strang:
        return advection_flow(drift_flow(advection_flow(xi, start, 0.5 * step), start, step), start + (0.5 * step), 0.5 * step);
    case qmomkit::PhaseSpaceEvolutionScheme::advsrc:
        return drift_flow(advection_flow(xi, start, step), start, step);
    case qmomkit::PhaseSpaceEvolutionScheme::srcadv:
        return advection_flow(drift_flow(xi, start, step), start, step);
    }
    return xi;
}

double relative_error(const std::vector<double>& moments, const std::vector<double>& reference)
{
    double error = 0.0;
    for (std::size_t order = 0; order < reference.size(); ++order)
    {
        error = std::max(error, std::abs(moments[order] - reference[order]) / std::max(1.0, std::abs(reference[order])));
    }
    return error;
}

// Errors of one run against the exact solution and against the exact split
// solution of its scheme.
struct RunErrors
{
    double exact = 0.0;
    double split = 0.0;
};

RunErrors run(const qmomkit::PhaseSpaceEvolutionScheme scheme, const std::size_t steps)
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

    std::array<double, 2> split_nodes = initial_nodes;
    const double time_step = final_time / static_cast<double>(steps);
    for (std::size_t step = 0; step < steps; ++step)
    {
        const double start = static_cast<double>(step) * time_step;
        const qmomkit::QmomPhaseSpaceEvolutionResult result = qmomkit::advance_qmom_phase_space(state, start, time_step, velocity, drift, options);
        QMOMKIT_VERIFY_CHECK(result.ok());
        state = result.state;
        for (double& node : split_nodes)
        {
            node = split_step(scheme, node, start, time_step);
        }
    }
    return {relative_error(state.moment_values(), exact_moments(final_time)), relative_error(state.moment_values(), moments_of(split_nodes))};
}

using Errors = std::array<double, step_counts.size()>;

// Runs every step count; returns the errors against the exact solution and
// against the exact split solution.
std::array<Errors, 2> run_all(const qmomkit::PhaseSpaceEvolutionScheme scheme)
{
    std::array<Errors, 2> errors{};
    for (std::size_t index = 0; index < step_counts.size(); ++index)
    {
        const RunErrors run_errors = run(scheme, step_counts.at(index));
        QMOMKIT_VERIFY_CHECK(std::isfinite(run_errors.exact) && std::isfinite(run_errors.split));
        errors.at(0).at(index) = run_errors.exact;
        errors.at(1).at(index) = run_errors.split;
    }
    return errors;
}

// Checks that the errors fall at every refinement and reach the given observed
// orders on the two finest pairs.
void verify_order(const Errors& errors, const double second_finest, const double finest)
{
    for (std::size_t index = 1; index < errors.size(); ++index)
    {
        QMOMKIT_VERIFY_CHECK(errors.at(index) < errors.at(index - 1));
    }
    const std::size_t last = errors.size() - 1;
    QMOMKIT_VERIFY_CHECK(std::log2(errors.at(last - 2) / errors.at(last - 1)) >= second_finest);
    QMOMKIT_VERIFY_CHECK(std::log2(errors.at(last - 1) / errors.at(last)) >= finest);
}

void verify_schemes_match_their_exact_split_solutions()
{
    for (const qmomkit::PhaseSpaceEvolutionScheme scheme :
         {qmomkit::PhaseSpaceEvolutionScheme::strang, qmomkit::PhaseSpaceEvolutionScheme::advsrc, qmomkit::PhaseSpaceEvolutionScheme::srcadv})
    {
        const Errors split_deviations = run_all(scheme).at(1);
        verify_order(split_deviations, 2.8, 2.8);
        QMOMKIT_VERIFY_CHECK(split_deviations.back() <= split_deviation_bound);
    }
}

void verify_strang_is_second_order()
{
    verify_order(run_all(qmomkit::PhaseSpaceEvolutionScheme::strang).at(0), 1.75, 1.85);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(argc, argv,
                                             {
                                                 {"schemes_match_their_exact_split_solutions", verify_schemes_match_their_exact_split_solutions},
                                                 {"strang_is_second_order", verify_strang_is_second_order},
                                             });
}
