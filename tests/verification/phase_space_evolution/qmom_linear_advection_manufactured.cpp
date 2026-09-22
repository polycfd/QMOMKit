// QMOM phase-space evolution against a problem whose exact answer is known in
// closed form, and a measurement of the splitting order.
//
// What is verified
//   direct_advection - that advecting the quadrature and then rebuilding the
//     moments reproduces the exact solution, and that advection alone moves
//     the positions while leaving the weights untouched.
//   phase_space_source - that combining advection with a source term converges
//     to the exact solution as the time step shrinks, at the order each
//     splitting scheme is supposed to achieve.
//
// How the reference is constructed
//   The velocity is chosen as u(x) = x, which is the case that can be solved
//   on paper. Each quadrature position then satisfies dx/dt = x, so it grows
//   as x(t) = x(0)*exp(t), and advection never changes a weight. The raw
//   moment of order k is a weighted sum of the k-th powers of the positions,
//   so it follows
//
//       M_k(t) = M_k(0) * exp(k*t).
//
//   Adding a source that injects 0.1 into every moment per unit time gives
//
//       dM_k/dt = k*M_k + 0.1,
//
//   whose solution is M_0(t) = M_0(0) + 0.1*t and, for k >= 1,
//       M_k(t) = M_k(0)*exp(k*t) + 0.1*(exp(k*t) - 1)/k.
//
//   That closed form is the reference. The starting moments come from three
//   chosen positions and weights, so the initial state is exact as well and
//   nothing here is taken from a previous QMOMKit run.
//
// Error metric
//   Maximum moment error relative to max(|expected|, 1). The observed order of
//   convergence is the base-2 logarithm of the ratio of errors between
//   successive halvings of the time step.
//
// Pass criteria
//   direct_advection - positions within 2.0e-12 absolute of x(0)*exp(0.1),
//     weights bit-for-bit unchanged, rebuilt moments within 2.0e-11 relative.
//   phase_space_source - the error must fall at every refinement, and over the
//     last two refinements the observed order must exceed 1.7 for Strang and
//     lie between 0.7 and 1.3 for the two first-order orderings. Those bounds
//     are the textbook orders of the schemes: second order for Strang, first
//     order for a single advection-then-source or source-then-advection pass.
//     They are stated as ranges because an observed order measured from two
//     finite steps is itself only approximate.
//
// Reproduce
//   ctest -R verification_phase_space_evolution_qmom_linear_advection
//
// Build modes
//   Both. The checks use QMOMKIT_VERIFY_CHECK, which reports and exits rather
//   than relying on assert(), so NDEBUG does not weaken this study.

#include "../verification_helpers.hpp"
#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

namespace
{

constexpr std::size_t node_count = 3;
constexpr double final_time = 0.1;

double power(const double value, const std::size_t exponent)
{
    double result = 1.0;
    for (std::size_t i = 0; i < exponent; ++i)
    {
        result *= value;
    }
    return result;
}

std::vector<double> initial_moments(const std::size_t count)
{
    constexpr std::array<double, node_count> nodes{0.5, 1.5, 3.0};
    constexpr std::array<double, node_count> weights{0.2, 0.3, 0.5};
    std::vector<double> moments(count, 0.0);
    for (std::size_t order = 0; order < count; ++order)
    {
        for (std::size_t node = 0; node < node_count; ++node)
        {
            moments[order] += weights[node] * power(nodes[node], order);
        }
    }
    return moments;
}

qmomkit::QmomOptions qmom_options()
{
    qmomkit::QmomOptions options;
    options.adaptive = false;
    options.automatic_scaling = false;
    return options;
}

qmomkit::QmomState initial_state(const std::size_t count)
{
    return qmomkit::QmomCase{qmomkit::MomentSet(initial_moments(count)), node_count, qmom_options()}.initial_state();
}

qmomkit::PhaseSpaceVelocity velocity()
{
    return [](double, const std::vector<double>& coordinates, std::vector<double>& velocities, std::size_t)
    { std::copy(coordinates.begin(), coordinates.end(), velocities.begin()); };
}

double error(const std::vector<double>& actual, const std::vector<double>& expected)
{
    QMOMKIT_VERIFY_CHECK(actual.size() == expected.size());
    double result = 0.0;
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        result = std::max(result, std::abs(actual[i] - expected[i]) / std::max(1.0, std::abs(expected[i])));
    }
    return result;
}

qmomkit::QmomPhaseSpaceEvolutionOptions options(const qmomkit::PhaseSpaceEvolutionScheme scheme)
{
    qmomkit::QmomPhaseSpaceEvolutionOptions result;
    result.node_options.scheme = qmomkit::NodeAdvectionScheme::rk45;
    result.node_options.adaptive_options.relative_tolerance = 1.0e-12;
    result.node_options.adaptive_options.absolute_tolerance = 1.0e-14;
    result.node_options.adaptive_options.initial_step = 0.0025;
    result.node_options.adaptive_options.max_step = 0.0025;
    result.scheme = scheme;
    return result;
}

void verify_direct_phase_space_advection()
{
    qmomkit::QmomState state = initial_state(11);
    const qmomkit::QmomResult initial_quadrature = state.quadrature();
    qmomkit::QmomPhaseSpaceEvolutionOptions advection_options = options(qmomkit::PhaseSpaceEvolutionScheme::strang);
    advection_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    constexpr double time_step = 0.00025;
    for (std::size_t step = 0; step < 400; ++step)
    {
        const qmomkit::QmomPhaseSpaceEvolutionResult result =
            qmomkit::advance_qmom_advection(state, step * time_step, time_step, velocity(), advection_options);
        QMOMKIT_VERIFY_CHECK(result.ok());
        state = result.state;
    }
    for (std::size_t node = 0; node < node_count; ++node)
    {
        QMOMKIT_VERIFY_CHECK(std::abs(state.quadrature().abscissas[node] - initial_quadrature.abscissas[node] * std::exp(final_time)) < 2.0e-12);
        QMOMKIT_VERIFY_CHECK(state.quadrature().weights[node] == initial_quadrature.weights[node]);
    }
    std::vector<double> expected = initial_moments(11);
    for (std::size_t order = 0; order < expected.size(); ++order)
    {
        expected[order] *= std::exp(static_cast<double>(order) * final_time);
    }
    const qmomkit::Quadrature1D final_quadrature{state.quadrature().abscissas, state.quadrature().weights, qmomkit::QuadratureStatus::success, ""};
    QMOMKIT_VERIFY_CHECK(error(qmomkit::reconstruct_raw_moments(final_quadrature, expected.size()), expected) < 2.0e-11);
}

qmomkit::QmomMomentSource injection_source()
{
    return [](double, const qmomkit::QmomSourceContext& context) { return qmomkit::QmomMomentRate{std::vector<double>(context.size(), 0.1)}; };
}

std::vector<double> exact_injection_solution()
{
    std::vector<double> result = initial_moments(6);
    result[0] += 0.1 * final_time;
    for (std::size_t order = 1; order < result.size(); ++order)
    {
        const double k = static_cast<double>(order);
        result[order] = result[order] * std::exp(k * final_time) + 0.1 * (std::exp(k * final_time) - 1.0) / k;
    }
    return result;
}

std::vector<double> advance(const double time_step, const qmomkit::PhaseSpaceEvolutionScheme scheme)
{
    qmomkit::QmomState state = initial_state(6);
    const std::size_t steps = static_cast<std::size_t>(std::lround(final_time / time_step));
    for (std::size_t step = 0; step < steps; ++step)
    {
        const qmomkit::QmomPhaseSpaceEvolutionResult result =
            qmomkit::advance_qmom_phase_space(state, step * time_step, time_step, velocity(), injection_source(), options(scheme));
        QMOMKIT_VERIFY_CHECK(result.ok());
        QMOMKIT_VERIFY_CHECK(result.diagnostics.moment_source_diagnostics.accepted);
        state = result.state;
    }
    return state.moments();
}

void verify_phase_space_source_evolution()
{
    const std::vector<double> time_steps{0.02, 0.01, 0.005, 0.0025};
    const std::array<qmomkit::PhaseSpaceEvolutionScheme, 3> schemes{
        qmomkit::PhaseSpaceEvolutionScheme::strang, qmomkit::PhaseSpaceEvolutionScheme::advsrc, qmomkit::PhaseSpaceEvolutionScheme::srcadv};
    std::array<std::vector<double>, 3> errors;
    for (std::size_t scheme = 0; scheme < schemes.size(); ++scheme)
    {
        for (const double time_step : time_steps)
        {
            errors[scheme].push_back(error(advance(time_step, schemes[scheme]), exact_injection_solution()));
        }
        for (std::size_t run = 1; run < errors[scheme].size(); ++run)
        {
            QMOMKIT_VERIFY_CHECK(errors[scheme][run] < errors[scheme][run - 1]);
        }
    }
    for (std::size_t run = errors[0].size() - 2; run < errors[0].size(); ++run)
    {
        const double strang_order = std::log(errors[0][run - 1] / errors[0][run]) / std::log(2.0);
        const double advsrc_order = std::log(errors[1][run - 1] / errors[1][run]) / std::log(2.0);
        const double srcadv_order = std::log(errors[2][run - 1] / errors[2][run]) / std::log(2.0);
        QMOMKIT_VERIFY_CHECK(strang_order > 1.7);
        QMOMKIT_VERIFY_CHECK(advsrc_order > 0.7 && advsrc_order < 1.3);
        QMOMKIT_VERIFY_CHECK(srcadv_order > 0.7 && srcadv_order < 1.3);
    }
}

} // namespace

int main(const int argc, char** argv)
{
    if (argc != 2)
    {
        return 2;
    }
    if (std::string(argv[1]) == "direct_advection")
    {
        verify_direct_phase_space_advection();
    }
    else if (std::string(argv[1]) == "phase_space_source")
    {
        verify_phase_space_source_evolution();
    }
    else
    {
        return 2;
    }
    return 0;
}
