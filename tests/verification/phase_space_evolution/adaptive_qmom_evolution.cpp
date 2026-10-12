// Adaptive QMOM evolution: the node count changes only in the reconstruction
// after the source step, growing by one node or dropping a vanishing one.
//
// What is verified
//   That an adaptive QMOM state grows from one to two nodes as soon as its
//   moments support two, and drops to one node when the weight of one of two
//   nodes vanishes; that the accepted quadrature after every step is the rule
//   the moments describe; and that the stored moments are those of that rule.
//
// How the reference is constructed
//   Chosen, not computed. The weights move linearly in time between two rules
//   on the fixed nodes 0.5 and 1.5, so every intermediate state is the
//   two-point rule with weights (1 - t) w_a + t w_b. Its moments change at the
//   constant rate M(w_b) - M(w_a), which SSP-RK3 integrates exactly, and the
//   velocity is zero, so the reference after each step is that rule:
//   - growth: from 1 at 1.5 to 0.5 at 0.5 and 0.5 at 1.5; the weight at 0.5 is
//     0.5 t, so two nodes from the first step on;
//   - reduction: the reverse; the weight at 0.5 is 0.5 (1 - t), so two nodes
//     until it vanishes at t = 1, then one node at 1.5 with weight 1.
//
// Error metric and pass criteria
//   The active node count after each step must equal the expected count
//   exactly. Positions, weights and moments must match the reference to
//   1.0e-12, relative to max(|reference|, 1). Observed errors are of order
//   1e-15, the rounding of the source step and the reconstruction.
//
// Reproduce
//   ctest -R verification_phase_space_evolution_adaptive_qmom_
//
// Build modes
//   Both. Checks go through QMOMKIT_VERIFY_*, not assert().

#include "../verification_helpers.hpp"
#include "qmomkit/qmomkit.hpp"
#include "test_case_dispatch.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace
{

constexpr double tolerance = 1.0e-12;
constexpr double time_step = 0.25;
constexpr std::array<double, 2> nodes{0.5, 1.5};

std::vector<double> raw_moments(const std::vector<double>& weights)
{
    std::vector<double> moments(4, 0.0);
    for (std::size_t order = 0; order < moments.size(); ++order)
    {
        for (std::size_t node = 0; node < nodes.size(); ++node)
        {
            moments[order] += weights[node] * std::pow(nodes[node], static_cast<int>(order));
        }
    }
    return moments;
}

std::vector<double> weights_at(const std::vector<double>& initial_weights, const std::vector<double>& final_weights, const double time)
{
    return {((1.0 - time) * initial_weights[0]) + (time * final_weights[0]), ((1.0 - time) * initial_weights[1]) + (time * final_weights[1])};
}

// Runs four steps of 0.25 to t = 1 and compares the state after each step with
// the two-point rule on {0.5, 1.5}, without the node whose weight is zero.
void verify_linear_weight_path(const std::vector<double>& initial_weights, const std::vector<double>& final_weights,
                               const std::size_t initial_node_count, const std::vector<std::size_t>& expected_node_counts)
{
    const std::vector<double> initial_moments = raw_moments(initial_weights);
    const std::vector<double> final_moments = raw_moments(final_weights);
    std::vector<double> rate(initial_moments.size());
    for (std::size_t order = 0; order < rate.size(); ++order)
    {
        rate[order] = final_moments[order] - initial_moments[order];
    }

    qmomkit::QmomCase qmom_case;
    qmom_case.moments = qmomkit::MomentSet(initial_moments);
    qmom_case.node_count = 2;
    qmom_case.qmom_options.adaptive = true;
    qmomkit::QmomState state = qmom_case.initial_state();
    QMOMKIT_VERIFY_CHECK(state.active_node_count() == initial_node_count);

    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const qmomkit::QmomMomentSource constant_rate = [&rate](double, const qmomkit::QmomSourceContext& context)
    { return qmomkit::QmomMomentRate{std::vector<double>(rate.begin(), rate.begin() + static_cast<std::ptrdiff_t>(context.size()))}; };

    double time = 0.0;
    for (const std::size_t expected_node_count : expected_node_counts)
    {
        const qmomkit::QmomPhaseSpaceEvolutionResult result = qmomkit::advance_qmom_phase_space(state, time, time_step, zero_velocity, constant_rate);
        QMOMKIT_VERIFY_CHECK(result.ok());
        state = result.state;
        time = result.final_time;
        QMOMKIT_VERIFY_CHECK(state.active_node_count() == expected_node_count);

        const std::vector<double> weights = weights_at(initial_weights, final_weights, time);
        std::vector<double> expected_nodes;
        std::vector<double> expected_weights;
        for (std::size_t node = 0; node < nodes.size(); ++node)
        {
            if (weights[node] != 0.0)
            {
                expected_nodes.push_back(nodes[node]);
                expected_weights.push_back(weights[node]);
            }
        }
        const qmomkit::QmomResult& quadrature = state.quadrature();
        QMOMKIT_VERIFY_CHECK(quadrature.abscissas.size() == expected_nodes.size());
        for (std::size_t node = 0; node < expected_nodes.size(); ++node)
        {
            QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(quadrature.abscissas[node], expected_nodes[node], 0.0, tolerance);
            QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(quadrature.weights[node], expected_weights[node], 0.0, tolerance);
        }
        const std::vector<double> expected_moments = raw_moments(weights);
        for (std::size_t order = 0; order < expected_moments.size(); ++order)
        {
            QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(state.moment_values()[order], expected_moments[order], 0.0, tolerance);
        }
    }
    QMOMKIT_VERIFY_NEAR(time, 1.0, 1.0e-15);
}

void verify_grows_by_one_node()
{
    // From one point at 1.5 to equal weights at 0.5 and 1.5: the reconstruction
    // tries one more node after the first step and keeps it.
    verify_linear_weight_path({0.0, 1.0}, {0.5, 0.5}, 1, {2, 2, 2, 2});
}

void verify_drops_a_vanishing_node()
{
    // From equal weights at 0.5 and 1.5 to one point at 1.5: two nodes while
    // the weight at 0.5 is positive, one node once it has vanished.
    verify_linear_weight_path({0.5, 0.5}, {0.0, 1.0}, 2, {2, 2, 2, 1});
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(argc, argv,
                                             {
                                                 {"grows_by_one_node", verify_grows_by_one_node},
                                                 {"drops_a_vanishing_node", verify_drops_a_vanishing_node},
                                             });
}
