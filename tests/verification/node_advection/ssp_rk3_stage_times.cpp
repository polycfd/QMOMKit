// SSP-RK3 node advection under velocities that depend on time only, against
// the exact displacement of the nodes.
//
// What is verified
//   That SSP-RK3 node advection evaluates the velocity at its stage times t,
//   t + dt and t + dt/2, in absolute time, for QMOM and CQMOM. The first two
//   stages enter with the same weight, so these cases cannot tell their times
//   apart. unit_node_advection_non_finite_velocity_diagnostic_reports_callback_context
//   checks that the first is at t. The other node-advection cases use
//   velocities that do not depend on time, and time_dependent_splitting.cpp
//   advects its nodes with RK45, so none of them would detect a wrong SSP-RK3
//   stage time.
//
// How the reference is constructed
//   By hand. The velocity of each coordinate is a cubic polynomial in t and
//   does not depend on the position, so every node moves by its integral,
//   P(t_end) - P(t_start), with P the antiderivative of the cubic. For such a
//   velocity, an SSP-RK3 step from t to t + dt reduces to Simpson's rule,
//   (dt/6) [v(t) + 4 v(t + dt/2) + v(t + dt)], which is exact for cubic
//   polynomials. The steps therefore reproduce the exact displacement apart
//   from rounding. The run starts at t = 0.5, so a stage evaluated at the time
//   elapsed within the step instead of the absolute time also fails. The
//   reference positions are the initial positions plus the displacement, so
//   the CQMOM reconstruction that provides the initial quadrature does not
//   enter the comparison.
//
// Error metric and pass criteria
//   Absolute deviation of each position from the reference after four steps of
//   0.25, at most 1e-12. Positions and displacements are of order one, so the
//   bound is relative as well. The largest deviation observed in the Debug,
//   Release, GCC 12 and GCC 16 builds is 6.7e-16. A stage at t instead of
//   t + dt, at t + dt instead of t + dt/2, or at the time since the start of
//   the step moves the nodes by at least 9.4e-3.
//
// Reproduce
//   ctest -R verification_node_advection_ssp_rk3_stage_times_
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

constexpr double start_time = 0.5;
constexpr double time_step = 0.25;
constexpr std::size_t steps = 4;
constexpr double tolerance = 1.0e-12;

// The coefficients c_0 to c_3 of the velocity c_0 + c_1 t + c_2 t^2 + c_3 t^3.
using Cubic = std::array<double, 4>;

constexpr Cubic velocity_0{0.6, 0.8, -1.1, 0.5};
constexpr Cubic velocity_1{-0.3, 0.9, -0.5, 0.8};

double evaluate(const Cubic& c, const double t)
{
    return c[0] + t * (c[1] + t * (c[2] + t * c[3]));
}

double displacement(const Cubic& c, const double from, const double to)
{
    const auto antiderivative = [&c](const double t) { return t * (c[0] + t * (c[1] / 2.0 + t * (c[2] / 3.0 + t * c[3] / 4.0))); };
    return antiderivative(to) - antiderivative(from);
}

qmomkit::NodeAdvectionOptions ssp_rk3_options()
{
    qmomkit::NodeAdvectionOptions options;
    options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    return options;
}

void test_qmom_nodes_move_by_the_exact_displacement()
{
    const std::vector<double> initial{-0.8, 0.4, 1.7};
    qmomkit::Quadrature1D quadrature{initial, {0.2, 0.35, 0.45}, qmomkit::QuadratureStatus::success, {}};
    const qmomkit::PhaseSpaceVelocity velocity = qmomkit::make_phase_space_velocity([](const double time, const std::vector<double>&)
                                                                                    { return std::vector<double>{evaluate(velocity_0, time)}; });

    double time = start_time;
    for (std::size_t step = 0; step < steps; ++step)
    {
        const qmomkit::QmomNodeAdvectionResult result = qmomkit::advance_nodes(quadrature, time, time_step, velocity, ssp_rk3_options());
        QMOMKIT_VERIFY_CHECK(result.ok());
        quadrature = result.quadrature;
        time += time_step;
    }

    const double moved = displacement(velocity_0, start_time, time);
    QMOMKIT_VERIFY_CHECK(quadrature.abscissas.size() == initial.size());
    for (std::size_t node = 0; node < initial.size(); ++node)
    {
        QMOMKIT_VERIFY_NEAR(quadrature.abscissas[node], initial[node] + moved, tolerance);
    }
}

void test_cqmom_nodes_move_by_the_exact_displacement()
{
    // Four weighted points, two positions in xi_1 at each of two in xi_0.
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.2, 0.3, 0.1, 0.4};
    const std::vector<std::vector<double>> points{{-0.5, -0.5, 1.2, 1.2}, {-1.0, 0.6, 0.2, 1.4}};
    std::vector<double> moments(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat = 0; flat < moments.size(); ++flat)
    {
        const qmomkit::MultiIndex index = qmomkit::unflatten_multi_index(flat, extents);
        for (std::size_t point = 0; point < weights.size(); ++point)
        {
            moments[flat] +=
                weights[point] * std::pow(points[0][point], static_cast<int>(index[0])) * std::pow(points[1][point], static_cast<int>(index[1]));
        }
    }
    const qmomkit::CqmomCase cqmom_case{qmomkit::extract_cqmom_moments(qmomkit::RawMomentTensor{extents, moments},
                                                                       qmomkit::create_cqmom_moment_layout(qmomkit::PhaseSpaceHierarchy({2, 2})))};
    const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(cqmom_case);
    QMOMKIT_VERIFY_CHECK(initial.ok());
    QMOMKIT_VERIFY_CHECK(initial.quadrature.total_nodes() == weights.size());

    const qmomkit::PhaseSpaceVelocity velocity = qmomkit::make_phase_space_velocity(
        [](const double time, const std::vector<double>&) { return std::vector<double>{evaluate(velocity_0, time), evaluate(velocity_1, time)}; });
    qmomkit::CqmomQuadrature quadrature = initial.quadrature;
    double time = start_time;
    for (std::size_t step = 0; step < steps; ++step)
    {
        const qmomkit::CqmomNodeAdvectionResult result = qmomkit::advance_nodes(quadrature, time, time_step, velocity, ssp_rk3_options());
        QMOMKIT_VERIFY_CHECK(result.ok());
        quadrature = result.quadrature;
        time += time_step;
    }

    const std::array<double, 2> moved{displacement(velocity_0, start_time, time), displacement(velocity_1, start_time, time)};
    QMOMKIT_VERIFY_CHECK(quadrature.total_nodes() == weights.size());
    for (std::size_t dimension = 0; dimension < moved.size(); ++dimension)
    {
        for (std::size_t node = 0; node < weights.size(); ++node)
        {
            QMOMKIT_VERIFY_NEAR(quadrature.abscissas()[dimension][node], initial.quadrature.abscissas()[dimension][node] + moved[dimension],
                                tolerance);
        }
    }
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(argc, argv,
                                             {
                                                 {"qmom_nodes_move_by_the_exact_displacement", test_qmom_nodes_move_by_the_exact_displacement},
                                                 {"cqmom_nodes_move_by_the_exact_displacement", test_cqmom_nodes_move_by_the_exact_displacement},
                                             });
}
