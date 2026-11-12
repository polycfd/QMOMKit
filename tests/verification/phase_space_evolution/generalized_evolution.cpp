// GQMOM states under the source step and under advection: closed-form moment
// evolutions of states that hold N_G nodes.
//
// What is verified
//   That a generalized QMOM state with node count N, N_G = 3N nodes and its
//   2N + 1 stored moments evolves as the count rule of
//   docs/developer/phase_space_evolution.md requires:
//   - decay on the nodes (N = 2, 3): the source S_k = -c sum_i w_i xi_i^k,
//     evaluated on the N_G nodes of each stage quadrature, keeps N_G nodes and the
//     base node count N, and the stored moments decay as the closed form;
//   - constant velocity (N = 2, 3): advection by a constant velocity shifts
//     the N_G nodes and the distribution, and keeps both node counts;
//   - growth from one point (adaptive, N = 3): relaxation toward the standard
//     normal distribution continues the state to N_G nodes after the first
//     step, and its base node count after step j is exactly min(1 + j, N);
//   - fixed boundary start (N = 2, 3): a fixed state of N points holds their
//     classical N nodes with base node count N, and after one relaxation step
//     its N_G nodes.
//
// How the reference is constructed
//   Closed form. Both sources are linear and autonomous in the moments, so one
//   SSP-RK3 step of size dt multiplies the deviation from equilibrium by
//   R(z) = 1 + z + z^2/2 + z^3/6, with z = -c dt for the decay, whose
//   equilibrium is zero, and z = -dt/theta for the relaxation. The decay is
//   evaluated on the nodes, and every stage quadrature reproduces
//   m_0, ..., m_{2N}, so its rate is -c m_k to rounding. A constant velocity v
//   moves every node by v t, so m_k(t) = sum_j C(k, j) (v t)^(k-j) m_j(0).
//   Every quadrature a generalized state holds reproduces m_0 and m_1: the
//   classical k-node quadrature reproduces m_0, ..., m_{2k-1} and a continued
//   one m_0, ..., m_{2k}, with k >= 1. Each moment relaxes on its own, so m_0
//   and m_1 follow the closed form while the base node count grows. The normal
//   start has m_0 = 2, mean 0.75 and variance 0.5.
//
// Error metric and pass criteria
//   The node counts must match exactly. Values must match the reference
//   relative to max(|reference|, 1): positions, m_0 and m_1 to 1.0e-12, and
//   the stored moments, up to order 6, to 1.0e-11. The only error between the
//   closed form and the evolution is the rounding of at most ten steps, each
//   with its stage reconstructions of the moments of a normal distribution
//   with N <= 3. The observed deviations are 1.7e-16 to 1.0e-15 for
//   positions, m_0 and m_1, and 1.3e-15 to 3.7e-15 for the stored moments, so
//   the tolerances leave a factor of at least 270. A wrong count rule or a
//   source evaluated on wrong nodes changes the counts, or the moments by far
//   more.
//
// Reproduce
//   ctest -R verification_phase_space_evolution_generalized_evolution_
//
// Build modes
//   Both. Checks go through QMOMKIT_VERIFY_*, not assert().

#include "../verification_helpers.hpp"
#include "qmomkit/qmomkit.hpp"
#include "test_case_dispatch.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace
{

constexpr double tolerance = 1.0e-12;
constexpr double moment_tolerance = 1.0e-11;
constexpr double time_step = 0.5;
constexpr double decay_rate = 0.4;
constexpr double theta = 1.0;
constexpr double velocity = 0.5;

// R(z) of one SSP-RK3 step for a linear rate with factor z = rate * dt.
double amplification(const double z)
{
    return 1.0 + z + (z * z / 2.0) + (z * z * z / 6.0);
}

std::vector<double> normal_moments(const double mass, const double mean, const double variance, const std::size_t count)
{
    std::vector<double> moments(count, 0.0);
    moments[0] = mass;
    moments[1] = mass * mean;
    for (std::size_t order = 2; order < count; ++order)
    {
        moments[order] = (mean * moments[order - 1]) + (static_cast<double>(order - 1) * variance * moments[order - 2]);
    }
    return moments;
}

std::vector<double> point_moments(const std::vector<double>& abscissas, const std::vector<double>& weights, const std::size_t count)
{
    std::vector<double> moments(count, 0.0);
    for (std::size_t order = 0; order < count; ++order)
    {
        for (std::size_t point = 0; point < abscissas.size(); ++point)
        {
            moments[order] += weights[point] * std::pow(abscissas[point], static_cast<int>(order));
        }
    }
    return moments;
}

qmomkit::QmomOptions generalized_options(const std::size_t n, const bool adaptive)
{
    qmomkit::QmomOptions options;
    options.adaptive = adaptive;
    options.family = qmomkit::QmomFamily::generalized;
    options.generalized.node_count = 3 * n;
    return options;
}

qmomkit::QmomMomentSource relaxation_towards(const std::vector<double>& equilibrium)
{
    return [equilibrium](double, const qmomkit::QmomSourceContext& context)
    {
        std::vector<double> rate(context.size());
        for (std::size_t order = 0; order < rate.size(); ++order)
        {
            rate[order] = (equilibrium[order] - context.value(order)) / theta;
        }
        return qmomkit::QmomMomentRate{std::move(rate)};
    };
}

void verify_near(const double actual, const double expected, const double relative_tolerance = tolerance)
{
    QMOMKIT_VERIFY_NEAR(actual, expected, relative_tolerance * std::max(std::abs(expected), 1.0));
}

void verify_decay_on_the_nodes()
{
    constexpr std::size_t steps = 10;
    const qmomkit::QmomMomentSource decay = [](double, const qmomkit::QmomSourceContext& context)
    {
        std::vector<double> rate(context.size(), 0.0);
        context.for_each_node(
            [&rate](const double weight, const double abscissa)
            {
                double power = 1.0;
                for (double& value : rate)
                {
                    value -= decay_rate * weight * power;
                    power *= abscissa;
                }
            });
        return qmomkit::QmomMomentRate{std::move(rate)};
    };
    for (const std::size_t n : {std::size_t{2}, std::size_t{3}})
    {
        const std::vector<double> initial = normal_moments(2.0, 0.75, 0.5, (2 * n) + 1);
        qmomkit::QmomState state = qmomkit::QmomCase{qmomkit::MomentSet(initial), n, generalized_options(n, true)}.initial_state();
        QMOMKIT_VERIFY_CHECK(state.active_node_count() == 3 * n && state.base_node_count() == n);
        for (std::size_t step = 1; step <= steps; ++step)
        {
            const qmomkit::QmomPhaseSpaceEvolutionResult result =
                qmomkit::advance_qmom_moment_source(state, static_cast<double>(step - 1) * time_step, time_step, decay);
            QMOMKIT_VERIFY_CHECK(result.ok());
            QMOMKIT_VERIFY_CHECK(result.state.active_node_count() == 3 * n && result.state.base_node_count() == n);
            const double factor = std::pow(amplification(-decay_rate * time_step), static_cast<double>(step));
            for (std::size_t order = 0; order < initial.size(); ++order)
            {
                verify_near(result.state.moment_values()[order], factor * initial[order], moment_tolerance);
            }
            state = result.state;
        }
    }
}

void verify_constant_velocity_shifts_the_distribution()
{
    const qmomkit::PhaseSpaceVelocity constant = [](double, const std::vector<double>&, std::vector<double>& rates, const std::size_t)
    { std::fill(rates.begin(), rates.end(), velocity); };
    for (const std::size_t n : {std::size_t{2}, std::size_t{3}})
    {
        const std::vector<double> initial = normal_moments(2.0, 0.75, 0.5, (2 * n) + 1);
        const qmomkit::QmomState state = qmomkit::QmomCase{qmomkit::MomentSet(initial), n, generalized_options(n, true)}.initial_state();
        const qmomkit::QmomPhaseSpaceEvolutionResult result = qmomkit::advance_qmom_advection(state, 0.0, time_step, constant);
        QMOMKIT_VERIFY_CHECK(result.ok());
        QMOMKIT_VERIFY_CHECK(result.state.active_node_count() == 3 * n && result.state.base_node_count() == n);
        const double shift = velocity * time_step;
        for (std::size_t node = 0; node < 3 * n; ++node)
        {
            verify_near(result.state.quadrature().abscissas[node], state.quadrature().abscissas[node] + shift);
        }
        const std::vector<double> expected = normal_moments(2.0, 0.75 + shift, 0.5, initial.size());
        for (std::size_t order = 0; order < initial.size(); ++order)
        {
            verify_near(result.state.moment_values()[order], expected[order], moment_tolerance);
        }
    }
}

void verify_growth_from_one_point()
{
    constexpr std::size_t n = 3;
    const std::vector<double> equilibrium = normal_moments(1.0, 0.0, 1.0, (2 * n) + 1);
    const std::vector<double> initial = point_moments({0.5}, {2.0}, (2 * n) + 1);
    const qmomkit::QmomMomentSource source = relaxation_towards(equilibrium);
    qmomkit::QmomState state = qmomkit::QmomCase{qmomkit::MomentSet(initial), n, generalized_options(n, true)}.initial_state();
    QMOMKIT_VERIFY_CHECK(state.active_node_count() == 1 && state.base_node_count() == 1);
    const double z = -time_step / theta;
    for (std::size_t step = 1; step <= n + 2; ++step)
    {
        const qmomkit::QmomPhaseSpaceEvolutionResult result =
            qmomkit::advance_qmom_moment_source(state, static_cast<double>(step - 1) * time_step, time_step, source);
        QMOMKIT_VERIFY_CHECK(result.ok());
        QMOMKIT_VERIFY_CHECK(result.state.active_node_count() == 3 * n);
        QMOMKIT_VERIFY_CHECK(result.state.base_node_count() == std::min(1 + step, n));
        const double factor = std::pow(amplification(z), static_cast<double>(step));
        for (const std::size_t order : {std::size_t{0}, std::size_t{1}})
        {
            verify_near(result.state.moment_values()[order], equilibrium[order] + (factor * (initial[order] - equilibrium[order])));
        }
        state = result.state;
    }
}

void verify_fixed_boundary_start()
{
    for (const std::size_t n : {std::size_t{2}, std::size_t{3}})
    {
        const std::vector<double> equilibrium = normal_moments(1.0, 0.0, 1.0, (2 * n) + 1);
        std::vector<double> abscissas(n);
        for (std::size_t point = 0; point < n; ++point)
        {
            abscissas[point] = 0.5 + (0.75 * static_cast<double>(point));
        }
        const std::vector<double> initial = point_moments(abscissas, std::vector<double>(n, 2.0 / static_cast<double>(n)), (2 * n) + 1);
        const qmomkit::QmomState state = qmomkit::QmomCase{qmomkit::MomentSet(initial), n, generalized_options(n, false)}.initial_state();
        QMOMKIT_VERIFY_CHECK(state.active_node_count() == n && state.base_node_count() == n);
        const qmomkit::QmomPhaseSpaceEvolutionResult result =
            qmomkit::advance_qmom_moment_source(state, 0.0, time_step, relaxation_towards(equilibrium));
        QMOMKIT_VERIFY_CHECK(result.ok());
        QMOMKIT_VERIFY_CHECK(result.state.active_node_count() == 3 * n && result.state.base_node_count() == n);
        const double factor = amplification(-time_step / theta);
        for (std::size_t order = 0; order < initial.size(); ++order)
        {
            verify_near(result.state.moment_values()[order], equilibrium[order] + (factor * (initial[order] - equilibrium[order])), moment_tolerance);
        }
    }
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(argc, argv,
                                             {
                                                 {"decay_on_the_nodes", verify_decay_on_the_nodes},
                                                 {"constant_velocity_shifts_the_distribution", verify_constant_velocity_shifts_the_distribution},
                                                 {"growth_from_one_point", verify_growth_from_one_point},
                                                 {"fixed_boundary_start", verify_fixed_boundary_start},
                                             });
}
