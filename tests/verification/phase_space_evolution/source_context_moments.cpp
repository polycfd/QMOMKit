// Moments that a moment-source callback requests beyond the tracked ones.
//
// What is verified
//   That the source context of a phase-space evolution step evaluates the
//   moments a callback requests but the state does not track, from the
//   quadrature of the current SSP-RK3 stage: integer orders above the tracked
//   ones through `QmomSourceContext::moment(...)`, multi-indices outside the
//   moment layout through `CqmomSourceContext::moment(...)`, and real orders
//   through `generalized_moment(...)` in both. Moments evaluated from another
//   quadrature, such as the one at the start of the step, make the step first
//   order.
//
// How the reference is constructed
//   By hand. The velocity is zero, and the source is the moment form of a node
//   motion that needs such moments, so the weights stay constant and the
//   nodes follow closed-form trajectories:
//   - QMOM, 0.4 at 0.5 and 0.6 at 1 to t = 0.5: dxi/dt = xi^2, written as
//     dm_k/dt = k m_{k+1}, so xi(t) = xi(0) / (1 - xi(0) t). With two nodes
//     the state tracks m_0 to m_3, and the rate of m_3 needs m_4.
//   - QMOM, 0.4 at 0.5 and 0.6 at 2 to t = 1: dxi/dt = sqrt(xi), written as
//     dm_k/dt = k m_{k - 1/2}, so xi(t) = (sqrt(xi(0)) + t / 2)^2.
//   - CQMOM, two nodes per coordinate, 0.15 at (0.5, 1), 0.2 at (0.5, 3),
//     0.3 at (1.5, 0.5) and 0.35 at (1.5, 2) to t = 0.5: xi_0 is constant and
//     dxi_1/dt = xi_0^2, written as dM_{i,j}/dt = j M_{i+2,j-1}, or
//     dxi_1/dt = sqrt(xi_0), written as dM_{i,j}/dt = j M_{i+1/2,j-1}; xi_1
//     moves by xi_0^2 t or sqrt(xi_0) t. The layout holds M_{i,0} for i <= 3
//     and M_{i,j} for i <= 1, j <= 3, so M_{2,1}, M_{2,2}, M_{3,1} and
//     M_{3,2} are untracked.
//   The moments at the final time follow from the positions.
//
// Error metric and pass criteria
//   Errors are the largest error of the tracked moments at the final time,
//   relative to max(|reference|, 1).
//   - QMOM: two nodes reproduce the four tracked stage moments, and the motion
//     keeps a measure of two points, so the only error is that of the
//     third-order SSP-RK3 source step. With 20, 40, 80 and 160 steps the
//     error must fall at every refinement, converge at an observed order
//     log2(e_dt / e_dt/2) of at least 2.8 on the two finest pairs, and stay at
//     most 1e-5 at 160 steps. The observed orders on those pairs are 2.94 to
//     3.00, and the errors at 160 steps are 1.7e-6 and 2.1e-8. The bound
//     leaves a factor of 6 over the larger error and lies three orders of
//     magnitude below the errors of 7e-3 and 6e-2 that the quadrature at the
//     start of the step gives.
//   - CQMOM: exact up to rounding. The xi_0 nodes do not move, and the stage
//     quadrature reproduces the conditional moments of orders 0 to 3 at each
//     of them, while the rates need orders up to 2 only. The moment equations
//     are then linear, and each rate lowers j by one, so their matrix A has
//     A^4 = 0 and the SSP-RK3 amplification I + hA + (hA)^2/2 + (hA)^3/6
//     equals the exact flow exp(hA). Four steps of 0.125 must match the
//     reference to 1.0e-12; observed errors are below 1e-14.
//
// Reproduce
//   ctest -R verification_phase_space_evolution_source_context_moments_
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
#include <functional>
#include <vector>

namespace
{

constexpr std::array<std::size_t, 4> step_counts{20, 40, 80, 160};
constexpr double minimum_order = 2.8;
constexpr double finest_error_bound = 1.0e-5;
constexpr double cqmom_tolerance = 1.0e-12;

constexpr std::array<double, 2> qmom_weights{0.4, 0.6};

// Position at the given time of a node that started at the given position.
using Trajectory = std::function<double(double position, double time)>;

std::vector<double> qmom_moments(const std::vector<double>& nodes)
{
    std::vector<double> moments(4, 0.0);
    for (std::size_t order = 0; order < moments.size(); ++order)
    {
        for (std::size_t node = 0; node < nodes.size(); ++node)
        {
            moments[order] += qmom_weights[node] * std::pow(nodes[node], static_cast<int>(order));
        }
    }
    return moments;
}

double relative_error(const std::vector<double>& moments, const std::vector<double>& reference)
{
    double error = 0.0;
    for (std::size_t entry = 0; entry < reference.size(); ++entry)
    {
        error = std::max(error, std::abs(moments[entry] - reference[entry]) / std::max(1.0, std::abs(reference[entry])));
    }
    return error;
}

qmomkit::PhaseSpaceVelocity zero_velocity()
{
    return qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
}

// Runs every step count to the final time and checks that the error falls at
// third order towards the reference.
void verify_qmom_convergence(const std::vector<double>& initial_nodes, const Trajectory& trajectory, const double final_time,
                             const qmomkit::QmomMomentSource& source)
{
    qmomkit::QmomCase qmom_case;
    qmom_case.moments = qmomkit::MomentSet(qmom_moments(initial_nodes));
    qmom_case.node_count = 2;
    qmom_case.qmom_options.adaptive = false;
    const std::vector<double> reference = qmom_moments({trajectory(initial_nodes[0], final_time), trajectory(initial_nodes[1], final_time)});

    const qmomkit::PhaseSpaceVelocity velocity = zero_velocity();
    std::array<double, step_counts.size()> errors{};
    for (std::size_t index = 0; index < step_counts.size(); ++index)
    {
        qmomkit::QmomState state = qmom_case.initial_state();
        const double time_step = final_time / static_cast<double>(step_counts.at(index));
        for (std::size_t step = 0; step < step_counts.at(index); ++step)
        {
            const qmomkit::QmomPhaseSpaceEvolutionResult result =
                qmomkit::advance_qmom_phase_space(state, static_cast<double>(step) * time_step, time_step, velocity, source);
            QMOMKIT_VERIFY_CHECK(result.ok());
            state = result.state;
        }
        errors.at(index) = relative_error(state.moment_values(), reference);
    }

    for (std::size_t index = 1; index < errors.size(); ++index)
    {
        QMOMKIT_VERIFY_CHECK(errors.at(index) < errors.at(index - 1));
    }
    const std::size_t last = errors.size() - 1;
    QMOMKIT_VERIFY_CHECK(std::log2(errors.at(last - 2) / errors.at(last - 1)) >= minimum_order);
    QMOMKIT_VERIFY_CHECK(std::log2(errors.at(last - 1) / errors.at(last)) >= minimum_order);
    QMOMKIT_VERIFY_CHECK(errors.at(last) <= finest_error_bound);
}

void verify_qmom_untracked_moment()
{
    const qmomkit::QmomMomentSource source = [](double, const qmomkit::QmomSourceContext& context)
    {
        qmomkit::QmomMomentRate rate{std::vector<double>(context.size(), 0.0)};
        for (std::size_t order = 1; order < context.size(); ++order)
        {
            rate.values[order] = static_cast<double>(order) * context.moment(order + 1);
        }
        return rate;
    };
    verify_qmom_convergence({0.5, 1.0}, [](const double position, const double time) { return position / (1.0 - (position * time)); }, 0.5, source);
}

void verify_qmom_generalized_moment()
{
    const qmomkit::QmomMomentSource source = [](double, const qmomkit::QmomSourceContext& context)
    {
        qmomkit::QmomMomentRate rate{std::vector<double>(context.size(), 0.0)};
        for (std::size_t order = 1; order < context.size(); ++order)
        {
            rate.values[order] = static_cast<double>(order) * context.generalized_moment(static_cast<double>(order) - 0.5);
        }
        return rate;
    };
    verify_qmom_convergence(
        {0.5, 2.0}, [](const double position, const double time) { return std::pow(std::sqrt(position) + (0.5 * time), 2); }, 1.0, source);
}

struct WeightedPoint
{
    double weight = 0.0;
    std::array<double, 2> position{};
};

constexpr std::array<WeightedPoint, 4> cqmom_points{{{0.15, {0.5, 1.0}}, {0.2, {0.5, 3.0}}, {0.3, {1.5, 0.5}}, {0.35, {1.5, 2.0}}}};

// Moments of the layout once xi_1 has moved by drift(xi_0) * time.
std::vector<double> cqmom_moments(const qmomkit::CqmomMomentLayout& layout, const std::function<double(double)>& drift, const double time)
{
    std::vector<double> values;
    for (const qmomkit::MultiIndexView index : layout.required_indices)
    {
        double value = 0.0;
        for (const WeightedPoint& point : cqmom_points)
        {
            const double xi_1 = point.position[1] + (drift(point.position[0]) * time);
            value += point.weight * std::pow(point.position[0], static_cast<int>(index[0])) * std::pow(xi_1, static_cast<int>(index[1]));
        }
        values.push_back(value);
    }
    return values;
}

// Runs four steps of 0.125 to t = 0.5 and compares the moments with the
// reference.
void verify_cqmom_exact(const std::function<double(double)>& drift, const qmomkit::CqmomMomentSource& source)
{
    constexpr double final_time = 0.5;
    constexpr std::size_t steps = 4;
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(qmomkit::PhaseSpaceHierarchy({2, 2}));
    qmomkit::CqmomCase cqmom_case{qmomkit::CqmomMomentData{layout, cqmom_moments(layout, drift, 0.0)}};
    cqmom_case.qmom_options.adaptive = false;
    qmomkit::CqmomState state = qmomkit::CqmomPreparedCase(cqmom_case).initial_state();
    const qmomkit::PhaseSpaceVelocity velocity = zero_velocity();

    const double time_step = final_time / static_cast<double>(steps);
    for (std::size_t step = 0; step < steps; ++step)
    {
        const qmomkit::CqmomPhaseSpaceEvolutionResult result =
            qmomkit::advance_cqmom_phase_space(state, static_cast<double>(step) * time_step, time_step, velocity, source);
        QMOMKIT_VERIFY_CHECK(result.ok());
        state = result.state;
    }
    QMOMKIT_VERIFY_CHECK(relative_error(state.moment_values(), cqmom_moments(layout, drift, final_time)) <= cqmom_tolerance);
}

void verify_cqmom_untracked_moment()
{
    const qmomkit::CqmomMomentSource source = [](double, const qmomkit::CqmomSourceContext& context)
    {
        qmomkit::CqmomMomentRate rate{std::vector<double>(context.size(), 0.0)};
        for (std::size_t entry = 0; entry < context.size(); ++entry)
        {
            const qmomkit::MultiIndexView index = context.index(entry);
            if (index[1] > 0)
            {
                const qmomkit::MultiIndex shifted{index[0] + 2, index[1] - 1};
                rate.values[entry] = static_cast<double>(index[1]) * context.moment(shifted);
            }
        }
        return rate;
    };
    verify_cqmom_exact([](const double xi_0) { return xi_0 * xi_0; }, source);
}

void verify_cqmom_generalized_moment()
{
    const qmomkit::CqmomMomentSource source = [](double, const qmomkit::CqmomSourceContext& context)
    {
        qmomkit::CqmomMomentRate rate{std::vector<double>(context.size(), 0.0)};
        for (std::size_t entry = 0; entry < context.size(); ++entry)
        {
            const qmomkit::MultiIndexView index = context.index(entry);
            if (index[1] > 0)
            {
                rate.values[entry] = static_cast<double>(index[1]) *
                                     context.generalized_moment({static_cast<double>(index[0]) + 0.5, static_cast<double>(index[1] - 1)});
            }
        }
        return rate;
    };
    verify_cqmom_exact([](const double xi_0) { return std::sqrt(xi_0); }, source);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(argc, argv,
                                             {
                                                 {"qmom_untracked_moment", verify_qmom_untracked_moment},
                                                 {"qmom_generalized_moment", verify_qmom_generalized_moment},
                                                 {"cqmom_untracked_moment", verify_cqmom_untracked_moment},
                                                 {"cqmom_generalized_moment", verify_cqmom_generalized_moment},
                                             });
}
