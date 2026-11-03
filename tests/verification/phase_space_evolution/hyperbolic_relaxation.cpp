// HyQMOM states under the source step alone: linear relaxation of the stored
// moments towards those of the standard normal distribution.
//
// What is verified
//   That a hyperbolic QMOM state with node count n and its 2n + 1 stored
//   moments evolves under the source step with zero velocity as the count
//   rule of docs/developer/phase_space_evolution.md requires:
//   - interior start (n = 1 to 4): the n + 1 nodes of HyQMOM after every
//     step, stored moments equal to the reference, and, once the deviation
//     from equilibrium is below rounding, the (n + 1)-point Gauss-Hermite rule
//     of the standard normal distribution, for which the closure is exact;
//   - growth from the boundary (adaptive, n = 2 to 4, one point or n points,
//     whose state has the classical one or n nodes): the active count after
//     step j is exactly min(k + j, n + 1), with k the number of points, and
//     m_0 and m_1 equal the reference.
//
// How the reference is constructed
//   Closed form. The rate S_k = (M^eq_k - M_k) / theta is linear and
//   autonomous, so one SSP-RK3 step of size dt multiplies M - M^eq by
//   R(z) = 1 + z + z^2/2 + z^3/6, z = -dt/theta, and after j steps
//   M = M^eq + R(z)^j (M(0) - M^eq), with no time-integration error. With
//   0 < dt/theta <= 1, each stage is a convex combination of the moments and
//   M^eq, so realizable moments stay realizable and become interior after one
//   step. Every rule reproduces m_0 and m_1, and each moment relaxes on its own,
//   so m_0 and m_1 follow the closed form even when a classical rule replaces
//   the higher stored moments by its own; the HyQMOM rule reproduces all of
//   m_0, ..., m_{2n}. M^eq_k is the k-th moment of the standard normal
//   distribution, and the normal start M_0 = 2, mean 0.75 and variance 0.5
//   follows from M_k = mean M_{k-1} + (k - 1) variance M_{k-2}. The
//   Gauss-Hermite nodes are the roots of He_{n+1} in closed form, and the
//   weights (n + 1)! / ((n + 1) He_n(x_i))^2.
//
// Error metric and pass criteria
//   The active node count must match exactly. Values must match the reference
//   relative to max(|reference|, 1): positions, weights, m_0 and m_1 to
//   1.0e-12, as for the exact references of this directory, with observed
//   errors below 1.5e-15; the stored moments of the interior case to 1.0e-11.
//   Those are the moments of the accepted rule, up to order 8 for n = 4, and
//   carry the conditioning of the map from nodes to moments: observed 2.4e-15
//   up to n = 2, 1.3e-14 for n = 3 and 1.4e-13 for n = 4, a margin of 70.
//
// Reproduce
//   ctest -R verification_phase_space_evolution_hyperbolic_relaxation_
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
constexpr double theta = 1.0;
constexpr double time_step = 0.5;

// R(z) of one SSP-RK3 step for the rate (M^eq - M) / theta.
double amplification()
{
    const double z = -time_step / theta;
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

std::vector<double> reference_moments(const std::vector<double>& initial, const std::vector<double>& equilibrium, const std::size_t steps)
{
    const double factor = std::pow(amplification(), static_cast<double>(steps));
    std::vector<double> moments(initial.size());
    for (std::size_t order = 0; order < moments.size(); ++order)
    {
        moments[order] = equilibrium[order] + (factor * (initial[order] - equilibrium[order]));
    }
    return moments;
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

qmomkit::QmomOptions hyperbolic_options()
{
    qmomkit::QmomOptions options;
    options.family = qmomkit::QmomFamily::hyperbolic;
    return options;
}

void verify_near(const double actual, const double expected, const double relative_tolerance = tolerance)
{
    QMOMKIT_VERIFY_NEAR(actual, expected, relative_tolerance * std::max(std::abs(expected), 1.0));
}

// He_k(x) by its three-term recurrence.
double hermite(const std::size_t degree, const double x)
{
    double previous = 1.0;
    double current = x;
    if (degree == 0)
    {
        return previous;
    }
    for (std::size_t k = 1; k < degree; ++k)
    {
        const double next = (x * current) - (static_cast<double>(k) * previous);
        previous = current;
        current = next;
    }
    return current;
}

// The roots of He_count in ascending order, count = 2 to 5.
std::vector<double> gauss_hermite_nodes(const std::size_t count)
{
    switch (count)
    {
    case 2:
        return {-1.0, 1.0};
    case 3:
        return {-std::sqrt(3.0), 0.0, std::sqrt(3.0)};
    case 4:
        return {-std::sqrt(3.0 + std::sqrt(6.0)), -std::sqrt(3.0 - std::sqrt(6.0)), std::sqrt(3.0 - std::sqrt(6.0)), std::sqrt(3.0 + std::sqrt(6.0))};
    default:
        return {-std::sqrt(5.0 + std::sqrt(10.0)), -std::sqrt(5.0 - std::sqrt(10.0)), 0.0, std::sqrt(5.0 - std::sqrt(10.0)),
                std::sqrt(5.0 + std::sqrt(10.0))};
    }
}

void verify_interior_moments_keep_n_plus_one_nodes()
{
    // 80 steps leave R(z)^80 = 3e-18 of the initial deviation.
    constexpr std::size_t steps = 80;
    for (std::size_t n = 1; n <= 4; ++n)
    {
        const std::size_t count = (2 * n) + 1;
        const std::vector<double> equilibrium = normal_moments(1.0, 0.0, 1.0, count);
        const std::vector<double> initial = normal_moments(2.0, 0.75, 0.5, count);
        const qmomkit::QmomMomentSource source = relaxation_towards(equilibrium);
        qmomkit::QmomState state = qmomkit::QmomCase{qmomkit::MomentSet(initial), n, hyperbolic_options()}.initial_state();
        QMOMKIT_VERIFY_CHECK(state.active_node_count() == n + 1);
        for (std::size_t step = 1; step <= steps; ++step)
        {
            const qmomkit::QmomPhaseSpaceEvolutionResult result =
                qmomkit::advance_qmom_moment_source(state, static_cast<double>(step - 1) * time_step, time_step, source);
            QMOMKIT_VERIFY_CHECK(result.ok());
            QMOMKIT_VERIFY_CHECK(result.state.active_node_count() == n + 1);
            QMOMKIT_VERIFY_CHECK(result.state.maximum_node_count() == n + 1);
            const std::vector<double> expected = reference_moments(initial, equilibrium, step);
            for (std::size_t order = 0; order < count; ++order)
            {
                verify_near(result.state.moment_values()[order], expected[order], moment_tolerance);
            }
            state = result.state;
        }

        const std::size_t nodes = n + 1;
        const std::vector<double> abscissas = gauss_hermite_nodes(nodes);
        double factorial = 1.0;
        for (std::size_t k = 2; k <= nodes; ++k)
        {
            factorial *= static_cast<double>(k);
        }
        for (std::size_t node = 0; node < nodes; ++node)
        {
            const double he = static_cast<double>(nodes) * hermite(n, abscissas[node]);
            verify_near(state.quadrature().abscissas[node], abscissas[node]);
            verify_near(state.quadrature().weights[node], factorial / (he * he));
        }
    }
}

void verify_boundary_moments_grow_to_n_plus_one_nodes()
{
    for (std::size_t n = 2; n <= 4; ++n)
    {
        const std::size_t count = (2 * n) + 1;
        const std::vector<double> equilibrium = normal_moments(1.0, 0.0, 1.0, count);
        const qmomkit::QmomMomentSource source = relaxation_towards(equilibrium);
        for (const std::size_t points : {std::size_t{1}, n})
        {
            std::vector<double> abscissas(points);
            for (std::size_t point = 0; point < points; ++point)
            {
                abscissas[point] = 0.5 + (0.75 * static_cast<double>(point));
            }
            const std::vector<double> initial = point_moments(abscissas, std::vector<double>(points, 2.0 / static_cast<double>(points)), count);
            qmomkit::QmomState state = qmomkit::QmomCase{qmomkit::MomentSet(initial), n, hyperbolic_options()}.initial_state();
            QMOMKIT_VERIFY_CHECK(state.active_node_count() == points);
            for (std::size_t step = 1; step <= n + 3 - points; ++step)
            {
                const qmomkit::QmomPhaseSpaceEvolutionResult result =
                    qmomkit::advance_qmom_moment_source(state, static_cast<double>(step - 1) * time_step, time_step, source);
                QMOMKIT_VERIFY_CHECK(result.ok());
                QMOMKIT_VERIFY_CHECK(result.state.active_node_count() == std::min(points + step, n + 1));
                QMOMKIT_VERIFY_CHECK(result.state.maximum_node_count() == n + 1);
                const std::vector<double> expected = reference_moments(initial, equilibrium, step);
                verify_near(result.state.moment_values()[0], expected[0]);
                verify_near(result.state.moment_values()[1], expected[1]);
                state = result.state;
            }
        }
    }
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(argc, argv,
                                             {
                                                 {"interior_moments_keep_n_plus_one_nodes", verify_interior_moments_keep_n_plus_one_nodes},
                                                 {"boundary_moments_grow_to_n_plus_one_nodes", verify_boundary_moments_grow_to_n_plus_one_nodes},
                                             });
}
