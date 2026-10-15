// Abscissa bounds during phase-space evolution.
//
// What is verified
//   That abscissa bounds apply wherever evolution reconstructs a quadrature
//   from moments, and nowhere else: advection moves a node out of the bounds
//   without clipping or rejecting it, in QMOM and CQMOM; a fixed-count
//   post-source reconstruction gives the exact rule while the nodes stay
//   inside the bounds and fails once one leaves them, returning the input
//   state and time; an adaptive one falls back to the one-node rule inside the
//   bounds; the source-stage reconstructions behind the source context apply
//   the bounds in the same way, in QMOM and CQMOM; and CQMOM bounds on one
//   coordinate leave the other coordinate unbounded.
//
// How the reference is constructed
//   By hand, from weighted points with bounds [0, 1] on the first coordinate.
//   - Advection: a constant velocity v moves every node exactly to x + v dt,
//     which the embedded RK23 integrator reproduces; the QMOM and CQMOM cases
//     below move nodes at 0.2 and 0.6 to 0.7 and 1.1.
//   - Source steps: the points move with constant velocities, and the source
//     returns the time derivative of their moments as an explicit function of
//     time, dm_k/dt = sum_j w_j k v_j (x_j + v_j (t - t0))^(k - 1). SSP-RK3
//     weighs its stage times t, t + dt and t + dt/2 like Simpson's rule, which
//     integrates a cubic exactly, so a step gives the moments of the moved
//     points whenever the rate is at most cubic in time: in QMOM, with m0 to
//     m3, it is at most quadratic; in CQMOM with node counts {2, 2}, where
//     every required index has a total order of at most 4, at most cubic.
//     The moved points are then the exact rule, and the one-node rule is
//     (m1/m0, m0).
//   - QMOM: nodes 0.2 and 0.6 with weights 0.4 and 0.6 move at 0.1 and 0.2
//     from t = 0.5: to 0.3 and 0.8 after a step of 1, inside the bounds, and
//     to 0.45 and 1.1 after a step of 2.5, outside them, where the one-node
//     rule is 0.84 with weight 1.
//   - Source stage: the QMOM nodes, advected by 0.5 to 0.7 and 1.1. A source
//     that visits the stage nodes and returns a zero rate then sees the
//     one-node rule 0.94 in every stage of an adaptive state, and fails the
//     stage of a fixed one; so does a CQMOM source that visits the terminal
//     nodes after the first coordinate is advected out of the bounds. Each
//     source step starts at the time the advection reached.
//   - CQMOM: the points (0.2, 0.05), (0.2, 0.95), (0.6, 0.1), (0.6, 0.8)
//     with weights 0.15, 0.2, 0.3, 0.35 and fixed node counts, all inside
//     [0, 1] in both coordinates. The first coordinate moves at 0.25 for 2, to
//     0.7 and 1.1, by advection or by the source; or both coordinates move at
//     (0.1, 0.15) for 2, the first to 0.4 and 0.8, inside its bounds, the
//     second to 0.35 and 1.25 in the first branch and 0.4 and 1.1 in the
//     second, so each branch has one node outside [0, 1].
//
// Error metric and pass criteria
//   Nodes, weights and moments within 1.0e-12 of the reference, relative to
//   max(|reference|, 1); the references are exact up to rounding, and the
//   observed errors are at most 4e-15. Statuses, stages and node counts must
//   match exactly, and a failed step must return its input state exactly and
//   its start time.
//
// Reproduce
//   ctest -R verification_phase_space_evolution_abscissa_bounds_
//
// Build modes
//   Both. Checks go through QMOMKIT_VERIFY_*, not assert().

#include "../verification_helpers.hpp"
#include "qmomkit/qmomkit.hpp"
#include "test_case_dispatch.hpp"

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

namespace
{

constexpr double tolerance = 1.0e-12;
constexpr qmomkit::AbscissaBounds unit_bounds{0.0, 1.0};

double integer_power(const double base, const std::size_t exponent)
{
    double value = 1.0;
    for (std::size_t power = 0; power < exponent; ++power)
    {
        value *= base;
    }
    return value;
}

// QMOM: nodes moving with constant velocities from the start time.
struct MovingNode
{
    double weight = 0.0;
    double position = 0.0;
    double velocity = 0.0;
};

constexpr double qmom_start = 0.5;
const std::vector<MovingNode> qmom_nodes{{0.4, 0.2, 0.1}, {0.6, 0.6, 0.2}};

std::vector<double> qmom_moments(const std::vector<MovingNode>& nodes, const double elapsed)
{
    std::vector<double> moments(4, 0.0);
    for (std::size_t order = 0; order < moments.size(); ++order)
    {
        for (const MovingNode& node : nodes)
        {
            moments[order] += node.weight * integer_power(node.position + node.velocity * elapsed, order);
        }
    }
    return moments;
}

qmomkit::QmomState qmom_state(const std::vector<MovingNode>& nodes, const bool adaptive)
{
    return qmomkit::QmomCase{qmomkit::MomentSet(qmom_moments(nodes, 0.0)), 2, qmomkit::QmomOptions{adaptive}, unit_bounds}.initial_state();
}

qmomkit::QmomPhaseSpaceEvolutionResult move_qmom_nodes(const qmomkit::QmomState& state, const double time_step)
{
    const qmomkit::QmomMomentSource translation = [](const double time, const qmomkit::QmomSourceContext& context)
    {
        qmomkit::QmomMomentRate rate{std::vector<double>(context.size(), 0.0)};
        for (std::size_t order = 1; order < rate.values.size(); ++order)
        {
            for (const MovingNode& node : qmom_nodes)
            {
                rate.values[order] += node.weight * static_cast<double>(order) * node.velocity *
                                      integer_power(node.position + node.velocity * (time - qmom_start), order - 1);
            }
        }
        return rate;
    };
    return qmomkit::advance_qmom_moment_source(state, qmom_start, time_step, translation);
}

void verify_qmom_quadrature(const qmomkit::QmomResult& quadrature, const std::vector<double>& weights, const std::vector<double>& abscissas)
{
    QMOMKIT_VERIFY_CHECK(quadrature.weights.size() == weights.size());
    QMOMKIT_VERIFY_CHECK(quadrature.abscissas.size() == abscissas.size());
    for (std::size_t node = 0; node < weights.size(); ++node)
    {
        QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(quadrature.weights[node], weights[node], 0.0, tolerance);
        QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(quadrature.abscissas[node], abscissas[node], 0.0, tolerance);
    }
}

void verify_qmom_moments(const qmomkit::QmomState& state, const std::vector<double>& expected)
{
    QMOMKIT_VERIFY_CHECK(state.moment_values().size() == expected.size());
    for (std::size_t order = 0; order < expected.size(); ++order)
    {
        QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(state.moment_values()[order], expected[order], 0.0, tolerance);
    }
}

void verify_qmom_failed_step(const qmomkit::QmomPhaseSpaceEvolutionResult& result, const qmomkit::QmomState& input, const double start_time)
{
    QMOMKIT_VERIFY_CHECK(result.state.moment_values() == input.moment_values());
    QMOMKIT_VERIFY_CHECK(result.state.quadrature().abscissas == input.quadrature().abscissas);
    QMOMKIT_VERIFY_CHECK(result.state.quadrature().weights == input.quadrature().weights);
    QMOMKIT_VERIFY_CHECK(result.final_time == start_time);
}

qmomkit::PhaseSpaceVelocity constant_velocity(const std::vector<double>& velocity)
{
    return qmomkit::make_phase_space_velocity([velocity](double, const std::vector<double>&) { return velocity; });
}

// The QMOM nodes advected by 0.5 from t = 0 to 1, to 0.7 and 1.1, outside the
// bounds.
qmomkit::QmomPhaseSpaceEvolutionResult advect_out_of_the_bounds(const bool adaptive)
{
    return qmomkit::advance_qmom_advection(qmom_state(qmom_nodes, adaptive), 0.0, 1.0, constant_velocity({0.5}));
}

void verify_qmom_advection_leaves_the_bounds_without_clipping()
{
    const qmomkit::QmomPhaseSpaceEvolutionResult result = advect_out_of_the_bounds(false);
    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.state.abscissa_bounds().has_value());
    verify_qmom_quadrature(result.state.quadrature(), {0.4, 0.6}, {0.7, 1.1});
    verify_qmom_moments(result.state, qmom_moments({{0.4, 0.7, 0.0}, {0.6, 1.1, 0.0}}, 0.0));
}

void verify_qmom_fixed_source_step_inside_the_bounds_is_exact()
{
    const qmomkit::QmomPhaseSpaceEvolutionResult result = move_qmom_nodes(qmom_state(qmom_nodes, false), 1.0);
    QMOMKIT_VERIFY_CHECK(result.ok());
    verify_qmom_quadrature(result.state.quadrature(), {0.4, 0.6}, {0.3, 0.8});
    verify_qmom_moments(result.state, qmom_moments(qmom_nodes, 1.0));
}

void verify_qmom_fixed_source_step_leaving_the_bounds_fails()
{
    const qmomkit::QmomState initial = qmom_state(qmom_nodes, false);
    const qmomkit::QmomPhaseSpaceEvolutionResult result = move_qmom_nodes(initial, 2.5);
    QMOMKIT_VERIFY_CHECK(!result.ok());
    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::PhaseSpaceEvolutionStatus::inversion_failed);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::post_source_inversion);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::success);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.post_source_qmom_status == qmomkit::QmomStatus::not_realizable);
    verify_qmom_failed_step(result, initial, qmom_start);
}

void verify_qmom_adaptive_source_step_leaving_the_bounds_keeps_one_node()
{
    const qmomkit::QmomPhaseSpaceEvolutionResult result = move_qmom_nodes(qmom_state(qmom_nodes, true), 2.5);
    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.state.active_node_count() == 1);
    verify_qmom_quadrature(result.state.quadrature(), {1.0}, {0.84});
    verify_qmom_moments(result.state, qmom_moments({{1.0, 0.84, 0.0}}, 0.0));
}

struct StageNode
{
    double weight = 0.0;
    double abscissa = 0.0;
};

// A zero rate that records the nodes of every source-stage quadrature.
qmomkit::QmomMomentSource node_recorder(std::vector<std::vector<StageNode>>& stages)
{
    return [&stages](double, const qmomkit::QmomSourceContext& context)
    {
        std::vector<StageNode> nodes;
        context.for_each_node([&nodes](const double weight, const double abscissa) { nodes.push_back({weight, abscissa}); });
        stages.push_back(nodes);
        return qmomkit::QmomMomentRate{std::vector<double>(context.size(), 0.0)};
    };
}

void verify_qmom_fixed_source_stage_rejects_nodes_outside_the_bounds()
{
    const qmomkit::QmomPhaseSpaceEvolutionResult advected = advect_out_of_the_bounds(false);
    QMOMKIT_VERIFY_CHECK(advected.ok());
    std::vector<std::vector<StageNode>> stages;
    const qmomkit::QmomPhaseSpaceEvolutionResult result =
        qmomkit::advance_qmom_moment_source(advected.state, advected.final_time, 0.1, node_recorder(stages));
    QMOMKIT_VERIFY_CHECK(!result.ok());
    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::PhaseSpaceEvolutionStatus::moment_source_failed);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::moment_source_step);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::source_failed);
    QMOMKIT_VERIFY_CHECK(stages.empty());
    verify_qmom_failed_step(result, advected.state, advected.final_time);
}

void verify_qmom_adaptive_source_stage_sees_one_node_inside_the_bounds()
{
    const qmomkit::QmomPhaseSpaceEvolutionResult advected = advect_out_of_the_bounds(true);
    QMOMKIT_VERIFY_CHECK(advected.ok());
    std::vector<std::vector<StageNode>> stages;
    const qmomkit::QmomPhaseSpaceEvolutionResult result =
        qmomkit::advance_qmom_moment_source(advected.state, advected.final_time, 0.1, node_recorder(stages));
    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(stages.size() == 3);
    for (const std::vector<StageNode>& nodes : stages)
    {
        QMOMKIT_VERIFY_CHECK(nodes.size() == 1);
        QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(nodes.front().weight, 1.0, 0.0, tolerance);
        QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(nodes.front().abscissa, 0.94, 0.0, tolerance);
    }
    QMOMKIT_VERIFY_CHECK(result.state.active_node_count() == 1);
    verify_qmom_quadrature(result.state.quadrature(), {1.0}, {0.94});
}

// CQMOM: points in the order of the reconstruction, by the first coordinate,
// then by the second within each first-coordinate node.
struct WeightedPoint
{
    double weight = 0.0;
    std::array<double, 2> position{};
};

constexpr double cqmom_start = 0.25;
constexpr double cqmom_time_step = 2.0;
const std::vector<WeightedPoint> cqmom_points{{0.15, {0.2, 0.05}}, {0.2, {0.2, 0.95}}, {0.3, {0.6, 0.1}}, {0.35, {0.6, 0.8}}};

std::vector<WeightedPoint> moved_points(const std::array<double, 2>& velocity, const double elapsed)
{
    std::vector<WeightedPoint> points = cqmom_points;
    for (WeightedPoint& point : points)
    {
        point.position[0] += velocity[0] * elapsed;
        point.position[1] += velocity[1] * elapsed;
    }
    return points;
}

const qmomkit::CqmomMomentLayout& cqmom_layout()
{
    static const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(qmomkit::PhaseSpaceHierarchy({2, 2}));
    return layout;
}

std::vector<double> cqmom_moments(const std::vector<WeightedPoint>& points)
{
    std::vector<double> values;
    for (const qmomkit::MultiIndexView index : cqmom_layout().required_indices)
    {
        double value = 0.0;
        for (const WeightedPoint& point : points)
        {
            value += point.weight * integer_power(point.position[0], index[0]) * integer_power(point.position[1], index[1]);
        }
        values.push_back(value);
    }
    return values;
}

qmomkit::CqmomState cqmom_state()
{
    qmomkit::CqmomCase cqmom_case{qmomkit::CqmomMomentData{cqmom_layout(), cqmom_moments(cqmom_points)}};
    cqmom_case.qmom_options.adaptive = false;
    cqmom_case.abscissa_bounds = qmomkit::CqmomAbscissaBounds{unit_bounds, std::nullopt};
    return qmomkit::CqmomPreparedCase{cqmom_case}.initial_state();
}

qmomkit::CqmomPhaseSpaceEvolutionResult move_cqmom_points(const qmomkit::CqmomState& state, const std::array<double, 2>& velocity)
{
    // d/dt of x^a y^b along the moving points is a v0 x^(a-1) y^b + b v1 x^a y^(b-1).
    const qmomkit::CqmomMomentSource translation = [velocity](const double time, const qmomkit::CqmomSourceContext& context)
    {
        const std::vector<WeightedPoint> points = moved_points(velocity, time - cqmom_start);
        qmomkit::CqmomMomentRate rate{std::vector<double>(context.size(), 0.0)};
        for (std::size_t entry = 0; entry < context.size(); ++entry)
        {
            const qmomkit::MultiIndexView index = context.index(entry);
            for (const WeightedPoint& point : points)
            {
                const double first = point.position[0];
                const double second = point.position[1];
                if (index[0] > 0)
                {
                    rate.values[entry] += point.weight * static_cast<double>(index[0]) * velocity[0] * integer_power(first, index[0] - 1) *
                                          integer_power(second, index[1]);
                }
                if (index[1] > 0)
                {
                    rate.values[entry] += point.weight * static_cast<double>(index[1]) * velocity[1] * integer_power(first, index[0]) *
                                          integer_power(second, index[1] - 1);
                }
            }
        }
        return rate;
    };
    return qmomkit::advance_cqmom_moment_source(state, cqmom_start, cqmom_time_step, translation);
}

void verify_cqmom_state(const qmomkit::CqmomState& state, const std::vector<WeightedPoint>& points)
{
    const qmomkit::CqmomQuadrature& quadrature = state.quadrature();
    QMOMKIT_VERIFY_CHECK(quadrature.total_nodes() == points.size());
    for (std::size_t node = 0; node < points.size(); ++node)
    {
        QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(quadrature.weights()[node], points[node].weight, 0.0, tolerance);
        for (std::size_t coordinate = 0; coordinate < 2; ++coordinate)
        {
            QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(quadrature.abscissas()[coordinate][node], points[node].position.at(coordinate), 0.0, tolerance);
        }
    }
    const std::vector<double> expected = cqmom_moments(points);
    QMOMKIT_VERIFY_CHECK(state.moment_values().size() == expected.size());
    for (std::size_t entry = 0; entry < expected.size(); ++entry)
    {
        QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(state.moment_values()[entry], expected[entry], 0.0, tolerance);
    }
}

void verify_cqmom_failed_step(const qmomkit::CqmomPhaseSpaceEvolutionResult& result, const qmomkit::CqmomState& input, const double start_time)
{
    QMOMKIT_VERIFY_CHECK(result.state.moment_values() == input.moment_values());
    QMOMKIT_VERIFY_CHECK(result.state.quadrature().abscissas() == input.quadrature().abscissas());
    QMOMKIT_VERIFY_CHECK(result.state.quadrature().weights() == input.quadrature().weights());
    QMOMKIT_VERIFY_CHECK(result.final_time == start_time);
}

void verify_cqmom_advection_leaves_the_bounds_without_clipping()
{
    const qmomkit::CqmomPhaseSpaceEvolutionResult result =
        qmomkit::advance_cqmom_advection(cqmom_state(), cqmom_start, cqmom_time_step, constant_velocity({0.25, 0.0}));
    QMOMKIT_VERIFY_CHECK(result.ok());
    verify_cqmom_state(result.state, moved_points({0.25, 0.0}, cqmom_time_step));
}

void verify_cqmom_fixed_source_stage_rejects_nodes_outside_the_bounds()
{
    const qmomkit::CqmomPhaseSpaceEvolutionResult advected =
        qmomkit::advance_cqmom_advection(cqmom_state(), cqmom_start, cqmom_time_step, constant_velocity({0.25, 0.0}));
    QMOMKIT_VERIFY_CHECK(advected.ok());
    std::size_t visits = 0;
    const qmomkit::CqmomMomentSource node_visitor = [&visits](double, const qmomkit::CqmomSourceContext& context)
    {
        context.for_each_terminal_node([&visits](double, const std::vector<double>&) { ++visits; });
        return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)};
    };
    const qmomkit::CqmomPhaseSpaceEvolutionResult result =
        qmomkit::advance_cqmom_moment_source(advected.state, advected.final_time, 0.1, node_visitor);
    QMOMKIT_VERIFY_CHECK(!result.ok());
    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::PhaseSpaceEvolutionStatus::moment_source_failed);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::moment_source_step);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::source_failed);
    QMOMKIT_VERIFY_CHECK(visits == 0);
    verify_cqmom_failed_step(result, advected.state, advected.final_time);
}

void verify_cqmom_bounded_coordinate_rejects_a_node_outside()
{
    const qmomkit::CqmomState initial = cqmom_state();
    const qmomkit::CqmomPhaseSpaceEvolutionResult result = move_cqmom_points(initial, {0.25, 0.0});
    QMOMKIT_VERIFY_CHECK(!result.ok());
    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::PhaseSpaceEvolutionStatus::inversion_failed);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::post_source_inversion);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::success);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.post_source_cqmom_status == qmomkit::CqmomStatus::qmom_failed);
    verify_cqmom_failed_step(result, initial, cqmom_start);
}

void verify_cqmom_unbounded_coordinate_accepts_nodes_outside()
{
    const qmomkit::CqmomPhaseSpaceEvolutionResult result = move_cqmom_points(cqmom_state(), {0.1, 0.15});
    QMOMKIT_VERIFY_CHECK(result.ok());
    verify_cqmom_state(result.state, moved_points({0.1, 0.15}, cqmom_time_step));
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"qmom_advection_leaves_the_bounds_without_clipping", verify_qmom_advection_leaves_the_bounds_without_clipping},
            {"qmom_fixed_source_step_inside_the_bounds_is_exact", verify_qmom_fixed_source_step_inside_the_bounds_is_exact},
            {"qmom_fixed_source_step_leaving_the_bounds_fails", verify_qmom_fixed_source_step_leaving_the_bounds_fails},
            {"qmom_adaptive_source_step_leaving_the_bounds_keeps_one_node", verify_qmom_adaptive_source_step_leaving_the_bounds_keeps_one_node},
            {"qmom_fixed_source_stage_rejects_nodes_outside_the_bounds", verify_qmom_fixed_source_stage_rejects_nodes_outside_the_bounds},
            {"qmom_adaptive_source_stage_sees_one_node_inside_the_bounds", verify_qmom_adaptive_source_stage_sees_one_node_inside_the_bounds},
            {"cqmom_advection_leaves_the_bounds_without_clipping", verify_cqmom_advection_leaves_the_bounds_without_clipping},
            {"cqmom_fixed_source_stage_rejects_nodes_outside_the_bounds", verify_cqmom_fixed_source_stage_rejects_nodes_outside_the_bounds},
            {"cqmom_bounded_coordinate_rejects_a_node_outside", verify_cqmom_bounded_coordinate_rejects_a_node_outside},
            {"cqmom_unbounded_coordinate_accepts_nodes_outside", verify_cqmom_unbounded_coordinate_accepts_nodes_outside},
        });
}
